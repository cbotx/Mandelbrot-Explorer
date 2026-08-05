// Headless verification + benchmark harness for the Mandelbrot renderer.
//
// For a given view it computes the escape-time buffer two ways:
//   (1) the perturbation / series-approximation / rebasing engine (Compute), and
//   (2) a brute-force full-precision reference (ComputeDirect / accuratePointCompute).
// It then reports how well they agree, plus timings and a checksum for
// regression tracking. No libpng / OpenGL / GLUT dependencies.
//
// Usage: verify [case]
//   case = shallow | deep | ticktock | flake | exterior1000 | parity1000 |
//          deep876 | subpixel | minibrot875 | extref875 | slowpoint | point31 |
//          gui875 | julia | julia-ede | julia-dendrite | julia-critical |
//          expression | expression-oracle | expression-suite |
//          expression-residual | multibrot | backend | gpu | all
//   The 1e1000 cases are excluded from "all" because their 3400-bit GMP oracles
//   are intentionally much more expensive than the regular regression set.

#include <gmp.h>
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <string>
#include <chrono>
#include <algorithm>
#include <future>
#include <thread>
#include <vector>

#include "compute_backend.h"
#include "mandel_perturbation.h"
#include "formula_expression.h"
#include "formula_expression_jit.h"
#include "formula_expression_oracle.h"
#include "test_cases.h"

using Clock = std::chrono::high_resolution_clock;
static double since(Clock::time_point t) {
    return std::chrono::duration_cast<std::chrono::milliseconds>(Clock::now() - t).count() / 1000.0;
}
struct TestCase {
    const char* name;
    std::string cx;
    std::string cy;
    std::string scale;   // decimal magnitude, e.g. "1" or "1e30" expanded to digits
    int mxit;
    int step;            // brute-force oracle samples a (step x step) pixel grid
    // Per-case pass tolerances and optional size override (0 => use main()'s W/H).
    // Fully-exterior stress frames (a sub-pixel nucleus casts no real interior)
    // set maxMisPct = 0 so ANY oracle classification disagreement fails.
    double maxDiff = 0.5;
    double maxMisPct = 1.0;
    double maxMeanDiff = 1e300;
    double maxP99Diff = 1e300;
    int w = 0, h = 0;
    const char* goldenPath = nullptr;
    const char* centerShiftRe = nullptr;
};

// Set-interior pixels carry the exact sentinel -2 (EMPTYPIXEL is -10); exterior
// smooth counts stay > -1.5. Testing "< 0" would misread a slightly-negative
// exterior smooth as interior, so classify against the sentinel band instead.
static inline bool isInterior(float v) { return v < -1.5f; }

static uint32_t checksum(const float* v, int n) {
    // Order-independent-ish FNV over the raw float bits.
    uint32_t h = 2166136261u;
    for (int i = 0; i < n; ++i) {
        uint32_t b;
        std::memcpy(&b, &v[i], 4);
        h = (h ^ b) * 16777619u;
    }
    return h;
}

static bool loadGolden(const char* path, float* out, size_t count) {
    FILE* f = fopen(path, "rb");
    std::string alternate;
    if (!f) {
        alternate = std::string("../") + path;   // verify.exe is normally run from build/
        f = fopen(alternate.c_str(), "rb");
    }

    if (!f) {
        fprintf(stderr, "Cannot open golden file: %s\n", path);
        return false;
    }
    size_t got = fread(out, sizeof(float), count, f);
    int extra = fgetc(f);
    fclose(f);
    if (got != count || extra != EOF) {
        fprintf(stderr, "Golden file has wrong size: %s (expected %zu floats)\n", path, count);
        return false;
    }
    return true;
}

static bool saveGolden(const char* path, const float* values, size_t count) {
    FILE* file = fopen(path, "wb");
    std::string alternate;
    if (!file) {
        alternate = std::string("../") + path;
        file = fopen(alternate.c_str(), "wb");
    }
    if (!file) {
        fprintf(stderr, "Cannot create golden file: %s\n", path);
        return false;
    }
    size_t written = fwrite(values, sizeof(float), count, file);
    bool ok = written == count && fclose(file) == 0;
    if (!ok) fprintf(stderr, "Failed writing golden file: %s\n", path);
    return ok;
}

static int runCase(const TestCase& tc, int W, int H) {
    if (tc.w > 0) W = tc.w;
    if (tc.h > 0) H = tc.h;
    // Precision follows the visible scale; 64 guard bits retain center digits
    // beyond the last visible decimal without making shallow tests use the full
    // textual coordinate length.
    int precision = static_cast<int>(tc.scale.size() * log(10) / log(2)) + 64;
    mpf_set_default_prec(precision);

    mpf_t cre, cim, scale;
    mpf_init_set_str(cre, tc.cx.c_str(), 10);
    mpf_init_set_str(cim, tc.cy.c_str(), 10);
    mpf_init_set_str(scale, tc.scale.c_str(), 10);
    if (tc.centerShiftRe) {
        mpf_t shift;
        mpf_init_set_str(shift, tc.centerShiftRe, 10);
        mpf_add(cre, cre, shift);
        mpf_clear(shift);
    }

    float* itp = new float[W * H];
    float* itd = new float[W * H];
    for (int i = 0; i < W * H; ++i) itp[i] = itd[i] = EMPTYPIXEL;

    Mandel mandel(W, H, tc.mxit, 1, itp);
    mandel.setPrecision(precision);

    // Optional EDE mode (exterior distance estimation) to benchmark the
    // accuracy-safe render path used by the GUI. Applied to both the
    // perturbation engine and the brute-force oracle so values stay comparable.
    int c_method = (getenv("MANDEL_EDE") != nullptr) ? ColoringMethod::EXTERIOR_DIST_EST : 0;

    printf("=== %s  (%dx%d, scale=1e%zu, mxit=%d, prec=%d bits%s)\n",
           tc.name, W, H, tc.scale.size() - 1, tc.mxit, precision,
           c_method ? ", EDE" : "");

    auto t0 = Clock::now();
    mandel.Compute(cre, cim, scale, tc.mxit, c_method);
    double t_pert = since(t0);

    t0 = Clock::now();
    const bool useGolden = tc.goldenPath && getenv("MANDEL_FORCE_ORACLE") == nullptr;
    if (useGolden) {
        if (!loadGolden(tc.goldenPath, itd, (size_t)W * H)) {
            delete[] itp; delete[] itd;
            mpf_clear(cre); mpf_clear(cim); mpf_clear(scale);
            return 2;
        }
    } else if (getenv("MANDEL_NOORACLE") == nullptr) {
        mandel.ComputeDirect(tc.mxit, itd, tc.step, c_method);
    }
    double t_direct = since(t0);

    // Compare only pixels the (possibly sparse) oracle actually computed.
    long n_int_both = 0, n_ext_both = 0, n_class_mismatch = 0, n_sampled = 0;
    long n_engine_empty = 0;
    double max_diff = 0, sum_diff = 0;
    std::vector<double> diffs;
    diffs.reserve((size_t)W * H);
    int worst_i = -1, worst_j = -1;
    int first_class_i = -1, first_class_j = -1;
    for (int i = 0; i < H; ++i) {
        for (int j = 0; j < W; ++j) {
            float d = itd[i * W + j];
            if (d == EMPTYPIXEL) continue;      // not sampled by the oracle
            float p = itp[i * W + j];
            ++n_sampled;
            if (p == EMPTYPIXEL) {
                if (first_class_i < 0) { first_class_i = i; first_class_j = j; }
                ++n_engine_empty;
                ++n_class_mismatch;
                continue;
            }
            bool p_int = isInterior(p);
            bool d_int = isInterior(d);
            if (p_int != d_int) {
                if (first_class_i < 0) { first_class_i = i; first_class_j = j; }
                ++n_class_mismatch;
                continue;
            }
            if (p_int) { ++n_int_both; continue; }
            ++n_ext_both;
            double diff = std::fabs((double)p - (double)d);
            sum_diff += diff;
            diffs.push_back(diff);
            if (diff > max_diff) { max_diff = diff; worst_i = i; worst_j = j; }
        }
    }

    double mean_diff = n_ext_both ? sum_diff / n_ext_both : 0.0;
    double p99_diff = 0.0;
    if (!diffs.empty()) {
        size_t q = (size_t)std::ceil(0.99 * diffs.size()) - 1;
        std::nth_element(diffs.begin(), diffs.begin() + q, diffs.end());
        p99_diff = diffs[q];
    }
    // Interior cost estimate over the FULL image: interior pixels (iter < 0)
    // each ran ~mxit iterations; exterior escape-iteration ~ the smooth value.
    long n_interior = 0; double ext_iters = 0;
    for (int q = 0; q < W * H; ++q) {
        if (itp[q] == EMPTYPIXEL) continue;
        if (isInterior(itp[q])) ++n_interior;
        else ext_iters += itp[q];
    }
    double int_iters = (double)n_interior * tc.mxit;
    double frac_int = (int_iters + ext_iters) > 0 ? int_iters / (int_iters + ext_iters) : 0.0;
    printf("  perturbation : %8.3f s   (full %dx%d image)\n", t_pert, W, H);
    printf("  interior px  : %ld / %d (%.1f%%)   est. iters interior/exterior = %.3g / %.3g  (interior = %.1f%% of all iters)\n",
           n_interior, W * H, 100.0 * n_interior / (W * H), int_iters, ext_iters, 100.0 * frac_int);
    if (useGolden)
        printf("  golden load  : %8.3f s   (%ld pixels from %s)\n", t_direct, n_sampled, tc.goldenPath);
    else
        printf("  direct (ref) : %8.3f s   (%ld pixels sampled, step=%d)\n",
               t_direct, n_sampled, tc.step);
    printf("  interior both: %ld   exterior both: %ld   class mismatch: %ld (%.3f%% of sampled)\n",
           n_int_both, n_ext_both, n_class_mismatch,
           n_sampled ? 100.0 * n_class_mismatch / n_sampled : 0.0);
    if (n_engine_empty)
        printf("  ERROR: engine left %ld sampled pixels EMPTY\n", n_engine_empty);
    if (first_class_i >= 0) {
        int q = first_class_i * W + first_class_j;
        printf("  first class mismatch @ %d,%d: pert=%.4f ref=%.4f\n",
               first_class_i, first_class_j, itp[q], itd[q]);
    }
    printf("  escape-time diff vs reference:  max %.6g  mean %.6g  p99 %.6g",
           max_diff, mean_diff, p99_diff);
    if (worst_i >= 0)
        printf("   (worst @ %d,%d: pert=%.4f ref=%.4f)",
               worst_i, worst_j, itp[worst_i * W + worst_j], itd[worst_i * W + worst_j]);
    printf("\n  checksum(pert)=0x%08x\n", checksum(itp, W * H));

    // Boundary pixels can legitimately disagree, but the bulk must match.
    // Fully-exterior stress frames use maxMisPct = 0: a sub-pixel nucleus casts
    // no true interior, so ANY false-interior pixel (e.g. the periodic-reference
    // "black star") or garbage smooth value (an escaping reference never
    // re-referenced) is a real glitch, not a boundary rounding.
    double mismatch_pct = n_sampled ? 100.0 * n_class_mismatch / n_sampled : 0.0;
    bool ok = n_engine_empty == 0 &&
              (max_diff <= tc.maxDiff) && (mean_diff <= tc.maxMeanDiff) &&
              (p99_diff <= tc.maxP99Diff) && (mismatch_pct <= tc.maxMisPct);
    printf("  => %s\n\n", ok ? "PASS" : "CHECK (see mismatches above)");

    delete[] itp;
    delete[] itd;
    mpf_clear(cre);
    mpf_clear(cim);
    mpf_clear(scale);
    return ok ? 0 : 1;
}

static int runAdaptiveGuiCase() {
    constexpr int W = 4, H = 4, SUB = 5, MXIT = 300000;
    std::string scaleText = "1";
    scaleText.append(875, '0');
    int precision = static_cast<int>(scaleText.size() * log(10) / log(2)) + 64;
    mpf_set_default_prec(precision);

    mpf_t cre, cim, scale;
    mpf_init_set_str(cre, testcases::deep1_x, 10);
    mpf_init_set_str(cim, testcases::deep1_y, 10);
    mpf_init_set_str(scale, scaleText.c_str(), 10);

    const size_t count = (size_t)W * H * SUB * SUB;
    std::vector<float> engine(count, EMPTYPIXEL), golden(count, EMPTYPIXEL);
    Mandel mandel(W, H, MXIT, SUB, engine.data());
    mandel.setPrecision(precision);

    printf("=== gui875 (production GUI adaptive SS5, resolved minibrot875)  "
           "(%dx%d, sub=%d, %zu sample slots)\n", W, H, SUB, count);
    auto t0 = Clock::now();
    mandel.Compute(cre, cim, scale, MXIT, ColoringMethod::SUPER_SAMPLING);
    double elapsed = since(t0);
    if (!loadGolden("tests/golden/minibrot875_gui_ss5.f32", golden.data(), count)) {
        mpf_clear(cre); mpf_clear(cim); mpf_clear(scale);
        return 2;
    }

    long computed = 0, base = 0, sub = 0, classMismatch = 0;
    int firstEmptyX = -1, firstEmptyY = -1;
    double maxDiff = 0.0, sumDiff = 0.0;
    std::vector<double> diffs;
    for (int y = 0; y < H * SUB; ++y) {
        for (int x = 0; x < W * SUB; ++x) {
            size_t q = (size_t)y * W * SUB + x;
            if (engine[q] == EMPTYPIXEL) {
                if (firstEmptyX < 0) { firstEmptyX = x; firstEmptyY = y; }
                continue;   // adaptive SS intentionally leaves unflagged subpixels empty
            }
            ++computed;
            bool isBase = (x % SUB == SUB / 2) && (y % SUB == SUB / 2);
            if (isBase) ++base; else ++sub;
            bool ei = isInterior(engine[q]), gi = isInterior(golden[q]);
            if (ei != gi) { ++classMismatch; continue; }
            if (!ei) {
                double d = std::fabs((double)engine[q] - golden[q]);
                maxDiff = std::max(maxDiff, d);
                sumDiff += d;
                diffs.push_back(d);
            }
        }
    }
    double meanDiff = diffs.empty() ? 0.0 : sumDiff / diffs.size();
    double p99Diff = 0.0;
    if (!diffs.empty()) {
        size_t q = (size_t)std::ceil(0.99 * diffs.size()) - 1;
        std::nth_element(diffs.begin(), diffs.begin() + q, diffs.end());
        p99Diff = diffs[q];
    }
    bool ok = base == W * H && sub > 0 && classMismatch == 0 && maxDiff <= 1.0;
    printf("  GUI engine   : %8.3f s   base=%ld/%d adaptive-subpixels=%ld\n",
           elapsed, base, W * H, sub);
    printf("  computed     : %ld/%zu   class mismatch: %ld\n", computed, count, classMismatch);
    if (firstEmptyX >= 0)
        printf("  first empty  : sample (%d,%d), base=%d\n",
               firstEmptyX, firstEmptyY,
               (firstEmptyX % SUB == SUB / 2) && (firstEmptyY % SUB == SUB / 2));
    printf("  escape-time diff vs golden: max %.6g  mean %.6g  p99 %.6g\n",
           maxDiff, meanDiff, p99Diff);
    printf("  => %s\n\n", ok ? "PASS" : "CHECK (GUI adaptive path mismatch)");

    mpf_clear(cre); mpf_clear(cim); mpf_clear(scale);
    return ok ? 0 : 1;
}

static int runJuliaCase(bool ede, bool dendrite = false, bool critical = false) {
    const int W = critical ? 32 : 96, H = critical ? 24 : 64;
    const int mxit = dendrite ? 1000 : 500;
    mpf_set_default_prec(256);
    mpf_t centerRe, centerIm, scale, fixedRe, fixedIm;
    mpf_init_set_ui(centerRe, 0); mpf_init_set_ui(centerIm, 0);
    mpf_init_set_ui(scale, critical ? 1000000 : 1);
    mpf_init_set_str(fixedRe, dendrite ? "-0.8" : "0", 10);
    mpf_init_set_str(fixedIm, dendrite ? "0.156" : "0", 10);
    std::vector<float> engine((size_t)W * H, EMPTYPIXEL);
    std::vector<float> oracle((size_t)W * H, EMPTYPIXEL);
    Mandel julia(W, H, mxit, 1, engine.data());
    const int method = ede ? ColoringMethod::EXTERIOR_DIST_EST : 0;

    auto t0 = Clock::now();
    // Build the exact viewport grid, then run the oracle before the SIMD engine so
    // memory-safety regressions in either path cannot silently bless each other.
    julia.ComputeJulia(centerRe, centerIm, scale, fixedRe, fixedIm, 1, method);
    julia.ComputeJuliaDirect(fixedRe, fixedIm, mxit, oracle.data(), 1, method);
    double oracleTime = since(t0);
    t0 = Clock::now();
    julia.ComputeJulia(centerRe, centerIm, scale, fixedRe, fixedIm, mxit, method);
    double engineTime = since(t0);

    int mismatch = 0, exterior = 0;
    int firstMismatch = -1;
    double maxDiff = 0.0, sumDiff = 0.0;
    std::vector<double> diffs;
    for (int i = 0; i < W * H; ++i) {
        bool a = isInterior(engine[i]), b = isInterior(oracle[i]);
        if (a != b) {
            if (firstMismatch < 0) firstMismatch = i;
            ++mismatch;
            continue;
        }
        if (!a) {
            ++exterior;
            double diff = std::fabs((double)engine[i] - oracle[i]);
            maxDiff = std::max(maxDiff, diff);
            sumDiff += diff;
            diffs.push_back(diff);
        }
    }
    double meanDiff = exterior ? sumDiff / exterior : 0.0;
    double p99Diff = 0.0;
    if (!diffs.empty()) {
        size_t q = (size_t)std::ceil(0.99 * diffs.size()) - 1;
        std::nth_element(diffs.begin(), diffs.begin() + q, diffs.end());
        p99Diff = diffs[q];
    }
    bool ok = mismatch == 0 &&
              (dendrite ? (maxDiff <= 100.0 && meanDiff <= 0.1 && p99Diff <= 0.5)
                        : (maxDiff <= (ede ? 1e-3 : 0.01)));
    if (critical) ok = ok && exterior == W * H;
    printf("=== quadratic Julia z0-plane (c=%s%s%s)  (%dx%d, mxit=%d)\n",
           dendrite ? "-0.8+0.156i" : "0", critical ? ", critical zoom" : "",
           ede ? ", EDE" : "", W, H, mxit);
    printf("  engine/direct: %.3f / %.3f s\n", engineTime, oracleTime);
    printf("  class mismatch: %d   max diff: %.6g   mean diff: %.6g   p99: %.6g\n",
           mismatch, maxDiff, meanDiff, p99Diff);
    if (critical) printf("  independently expected exterior: %d/%d\n", exterior, W * H);
    if (firstMismatch >= 0)
        printf("  first mismatch: (%d,%d) engine=%.6g oracle=%.6g\n",
               firstMismatch % W, firstMismatch / W,
               engine[firstMismatch], oracle[firstMismatch]);
    printf("  checksum=0x%08x\n", checksum(engine.data(), W * H));
    printf("  => %s\n\n", ok ? "PASS" : "CHECK (Julia direct/oracle mismatch)");

    mpf_clears(centerRe, centerIm, scale, fixedRe, fixedIm, (mpf_ptr)0);
    return ok ? 0 : 1;
}

static int runExpressionCoreCase() {
    using formula::Complex;
    using formula::ExpressionContext;
    using formula::ExpressionDerivativeSeed;
    using formula::ExpressionError;
    using formula::ExpressionProgram;

    int failures = 0;
    auto close = [](Complex a, Complex b, double eps = 1e-12) {
        return std::abs(a - b) <= eps * std::max(1.0, std::abs(b));
    };
    auto compile = [&](ExpressionProgram& program, const char* source) {
        ExpressionError error;
        if (!program.compile(source, &error)) {
            printf("  compile failed @ %zu: %s  [%s]\n",
                   error.position, error.message.c_str(), source);
            ++failures;
            return false;
        }
        return true;
    };

    ExpressionContext context;
    context.z = { 1.0, 2.0 };
    context.c = { 0.5, -0.25 };
    context.z0 = { -0.75, 0.1 };
    context.parameters[0] = { 0.2, -0.3 };
    context.iteration = 7;

    ExpressionProgram quadratic;
    if (compile(quadratic, "z*z + c") &&
        !close(quadratic.evaluate(context), Complex{ -2.5, 3.75 }))
        ++failures;
    if (quadratic.fastPath() != ExpressionProgram::FastPath::IntegerPowerPlusC ||
        quadratic.fastIntegerPower() != 2)
        ++failures;
    ExpressionProgram quadraticSquare, quadraticPower, cubicProduct;
    if (!compile(quadraticSquare, "sqr(z)+c") ||
        !compile(quadraticPower, "z^2+c") ||
        !compile(cubicProduct, "z*z*z+c") ||
        quadraticSquare.fastPath() != ExpressionProgram::FastPath::IntegerPowerPlusC ||
        quadraticSquare.fastIntegerPower() != 2 ||
        quadraticPower.fastPath() != ExpressionProgram::FastPath::None ||
        cubicProduct.fastPath() != ExpressionProgram::FastPath::IntegerPowerPlusC ||
        cubicProduct.fastIntegerPower() != 3)
        ++failures;

    ExpressionProgram precedence;
    if (compile(precedence, "-2^2 + 2^3^2") &&
        !close(precedence.evaluate(context), Complex{ 508.0, 0.0 }, 1e-10))
        ++failures;

    ExpressionProgram functions;
    const char* functionSource =
        "sin(z) + conj(p0) + complex(abs(re(c)), abs(im(c))) + n/10";
    if (compile(functions, functionSource)) {
        Complex expected = std::sin(context.z) + std::conj(context.parameters[0]) +
                           Complex{ std::abs(context.c.real()), std::abs(context.c.imag()) } +
                           Complex{ 0.7, 0.0 };
        if (!close(functions.evaluate(context), expected)) ++failures;
        if (functions.fastPath() != ExpressionProgram::FastPath::None) ++failures;
        if (functions.avx2Compatible()) ++failures;
    }

    ExpressionProgram burningShip;
    if (compile(burningShip,
                "sqr(complex(abs(re(z)), abs(im(z)))) + c")) {
        Complex folded{ std::abs(context.z.real()), std::abs(context.z.imag()) };
        if (!close(burningShip.evaluate(context), folded * folded + context.c))
            ++failures;
        if (burningShip.derivativeCompatible()) ++failures;
        ExpressionDerivativeSeed seed;
        Complex value, derivative;
        if (burningShip.evaluateWithDerivative(context, seed, value, derivative))
            ++failures;
    }

    ExpressionDerivativeSeed quadraticSeed;
    quadraticSeed.z = { 0.3, -0.2 };
    quadraticSeed.c = { -0.1, 0.4 };
    Complex dualValue, dualDerivative;
    if (!quadratic.derivativeCompatible() ||
        !quadratic.evaluateWithDerivative(context, quadraticSeed,
                                          dualValue, dualDerivative) ||
        dualValue != context.z * context.z + context.c ||
        dualDerivative != quadraticSeed.z * context.z +
                          context.z * quadraticSeed.z + quadraticSeed.c)
        ++failures;

    ExpressionProgram analytic;
    if (!compile(analytic, "sin(z*z+c)+exp(p0*z0)+log(z+2)") ||
        !analytic.derivativeCompatible()) {
        ++failures;
    } else {
        double maxDerivativeError = 0.0;
        for (int sample = 0; sample < 1000; ++sample) {
            ExpressionContext base = context;
            base.z = { 0.2 + sample * 1e-5, -0.15 + sample * 2e-6 };
            base.c = { 0.4 - sample * 3e-6, 0.12 };
            base.z0 = { -0.3, 0.25 + sample * 1e-6 };
            base.parameters[0] = { 0.6, -0.2 };
            ExpressionDerivativeSeed direction;
            direction.z = { 0.3, -0.1 };
            direction.c = { -0.2, 0.25 };
            direction.z0 = { 0.1, 0.2 };
            direction.parameters[0] = { -0.15, 0.05 };
            Complex value, derivative;
            if (!analytic.evaluateWithDerivative(base, direction, value, derivative)) {
                ++failures;
                break;
            }
            const double h = 1e-6;
            auto shifted = [&](double sign) {
                ExpressionContext shifted = base;
                shifted.z += sign * h * direction.z;
                shifted.c += sign * h * direction.c;
                shifted.z0 += sign * h * direction.z0;
                shifted.parameters[0] += sign * h * direction.parameters[0];
                return analytic.evaluate(shifted);
            };
            Complex finiteDifference = (shifted(1.0) - shifted(-1.0)) / (2.0 * h);
            double relative = std::abs(derivative - finiteDifference) /
                              std::max(1.0, std::abs(derivative));
            maxDerivativeError = std::max(maxDerivativeError, relative);
        }
        if (!(maxDerivativeError < 2e-9)) {
            printf("  derivative finite-difference max error=%.6g\n",
                   maxDerivativeError);
            ++failures;
        }
    }
    ExpressionProgram derivativeDomain;
    if (!compile(derivativeDomain, "log(z)") ||
        !derivativeDomain.derivativeCompatible()) {
        ++failures;
    } else {
        ExpressionContext zero;
        ExpressionDerivativeSeed seed; seed.z = 1.0;
        Complex value, derivative;
        if (derivativeDomain.evaluateWithDerivative(zero, seed, value, derivative))
            ++failures;
    }
    auto derivativeCase = [&](const char* source, Complex z, Complex seedValue,
                              Complex expectedValue, Complex expectedDerivative,
                              bool shouldSucceed) {
        ExpressionProgram program;
        if (!compile(program, source)) return;
        ExpressionContext dc; dc.z = z;
        ExpressionDerivativeSeed seed; seed.z = seedValue;
        Complex value, derivative;
        bool succeeded = program.evaluateWithDerivative(dc, seed, value, derivative);
        if (succeeded != shouldSucceed) { ++failures; return; }
        if (succeeded &&
            (!close(value, expectedValue, 1e-12) ||
             !close(derivative, expectedDerivative, 1e-12)))
            ++failures;
    };
    derivativeCase("z^2+c", {}, 1.0, {}, {}, true);
    {
        ExpressionProgram powerAtZero;
        if (!compile(powerAtZero, "z^2+c")) {
            ++failures;
        } else {
            ExpressionContext dc; dc.z = {}; dc.c = { 0.3, -0.2 };
            ExpressionDerivativeSeed seed; seed.z = 1.0; seed.c = 1.0;
            Complex value, derivative;
            if (!powerAtZero.evaluateWithDerivative(dc, seed, value, derivative) ||
                value != dc.c || derivative != Complex{ 1.0, 0.0 })
                ++failures;
        }
    }
    derivativeCase("log(z)", { -1.0, 0.0 }, 1.0, {}, {}, false);
    derivativeCase("sqrt(z)", { -1.0, 0.0 }, 1.0, {}, {}, false);
    derivativeCase("z^0.5", { -1.0, 0.0 }, 1.0, {}, {}, false);
    derivativeCase("z^(-1)", { -1.0, 0.0 }, 1.0,
                   { -1.0, 0.0 }, { -1.0, 0.0 }, true);
    derivativeCase("z/z", { 1e-200, 0.0 }, 1.0,
                   { 1.0, 0.0 }, {}, true);
    derivativeCase("tan(z)", { 0.0, 1000.0 }, 1.0,
                   { 0.0, 1.0 }, {}, true);
    derivativeCase("tanh(z)", { 1000.0, 0.0 }, 1.0,
                   { 1.0, 0.0 }, {}, true);
    {
        Complex z{ 0.0, 20.0 };
        Complex expected = 1.0 / (std::cos(z) * std::cos(z));
        derivativeCase("tan(z)", z, 1.0, std::tan(z), expected, true);
        z = { 20.0, 0.0 };
        expected = 1.0 / (std::cosh(z) * std::cosh(z));
        derivativeCase("tanh(z)", z, 1.0, std::tanh(z), expected, true);
        double q = std::exp(-600.0);
        double saturatedDerivative = 4.0 * q / ((1.0 + q) * (1.0 + q));
        derivativeCase("tan(z)", { 0.0, 300.0 }, 1.0,
                       std::tan(Complex{ 0.0, 300.0 }),
                       { saturatedDerivative, 0.0 }, true);
        derivativeCase("tanh(z)", { 300.0, 0.0 }, 1.0,
                       std::tanh(Complex{ 300.0, 0.0 }),
                       { saturatedDerivative, 0.0 }, true);
    }

    ExpressionProgram invalid;
    ExpressionError invalidError;
    if (invalid.compile("z + * c", &invalidError) || invalidError.message.empty())
        ++failures;
    std::string tooLong(ExpressionProgram::MAX_SOURCE + 1, '1');
    if (invalid.compile(tooLong, &invalidError)) ++failures;
    std::string tooMany = "z";
    for (size_t i = 0; i < ExpressionProgram::MAX_INSTRUCTIONS; ++i)
        tooMany += "+z";
    if (invalid.compile(tooMany, &invalidError)) ++failures;
    std::string deepPower = "z";
    for (size_t i = 0; i <= ExpressionProgram::MAX_PARSE_DEPTH; ++i)
        deepPower += "^z";
    if (invalid.compile(deepPower, &invalidError)) ++failures;
    std::string unaryChain(ExpressionProgram::MAX_SOURCE - 1, '-');
    unaryChain += "z";
    if (!invalid.compile(unaryChain, &invalidError)) ++failures;

    ExpressionProgram invalidPolar;
    if (compile(invalidPolar, "polar(-1, 0)")) {
        Complex result = invalidPolar.evaluate(context);
        if (!std::isnan(result.real()) || !std::isnan(result.imag())) ++failures;
    }

    // The immutable bytecode must be safe to evaluate concurrently.
    std::atomic<int> parallelFailures{0};
#pragma omp parallel for schedule(static)
    for (int i = 0; i < 10000; ++i) {
        ExpressionContext local = context;
        local.z = { i * 1e-4, -i * 2e-4 };
        Complex interpreted = quadratic.evaluate(local);
        std::array<Complex, ExpressionProgram::MAX_STACK> scratch;
        Complex reused = quadratic.evaluate(local, scratch.data(), quadratic.stackDepth());
        Complex expected = local.z * local.z + local.c;
        if (!close(interpreted, expected) || !close(reused, expected))
            ++parallelFailures;
    }
    failures += parallelFailures.load();

    ExpressionProgram vectorProgram;
    formula::ExpressionJit4 vectorJit;
    double avxMs = 0.0, jitMs = 0.0, jitRawMs = 0.0;
    int jitMismatch = 0;
    if (!compile(vectorProgram, "z*z+c+p0*z+complex(re(z0),im(c))") ||
        !vectorProgram.avx2Compatible() || !vectorJit.compile(vectorProgram)) {
        ++failures;
    } else {
        MEMORY_BASIC_INFORMATION memory{};
        SIZE_T queried = VirtualQuery(vectorJit.codeAddress(), &memory,
                                      sizeof(memory));
        DWORD executableWritable =
            PAGE_EXECUTE_READWRITE | PAGE_EXECUTE_WRITECOPY;
        if (!vectorJit.usesDualMapping() || queried != sizeof(memory) ||
            !(memory.Protect & (PAGE_EXECUTE | PAGE_EXECUTE_READ |
                                PAGE_EXECUTE_READWRITE | PAGE_EXECUTE_WRITECOPY)) ||
            (memory.Protect & executableWritable))
            ++failures;
        for (int group = 0; group < 2500; ++group) {
            ExpressionContext lanes[4];
            Complex vectorResults[4];
            Complex jitResults[4];
            for (int lane = 0; lane < 4; ++lane) {
                int index = group * 4 + lane;
                lanes[lane] = context;
                lanes[lane].z = { index * 1e-4, -index * 3e-5 };
                lanes[lane].c = { 0.2 + lane * 0.01, -0.3 + group * 1e-6 };
                lanes[lane].z0 = { -0.5 + lane * 0.02, 0.1 };
            }
            if (!vectorProgram.evaluate4(lanes, vectorResults)) {
                ++failures;
                break;
            }
            if (!vectorJit.evaluate(lanes, jitResults)) {
                ++failures;
                break;
            }
            for (int lane = 0; lane < 4; ++lane) {
                Complex scalar = vectorProgram.evaluate(lanes[lane]);
                if (scalar.real() != vectorResults[lane].real() ||
                    scalar.imag() != vectorResults[lane].imag())
                    ++failures;
                if (scalar.real() != jitResults[lane].real() ||
                    scalar.imag() != jitResults[lane].imag()) {
                    if (jitMismatch++ == 0)
                        printf("  first JIT mismatch scalar=(%.17g,%.17g) jit=(%.17g,%.17g)\n",
                               scalar.real(), scalar.imag(),
                               jitResults[lane].real(), jitResults[lane].imag());
                }
            }
        }
        auto sameDoubleBits = [](double a, double b) {
            uint64_t aa, bb;
            std::memcpy(&aa, &a, sizeof(a));
            std::memcpy(&bb, &b, sizeof(b));
            return aa == bb;
        };
        for (const char* source : { "-z", "conj(z)" }) {
            ExpressionProgram signProgram;
            if (!compile(signProgram, source) || !signProgram.avx2Compatible()) {
                ++failures;
                continue;
            }
            ExpressionContext lanes[4]{};
            lanes[0].z = { 0.0, 0.0 };
            lanes[1].z = { -0.0, 0.0 };
            lanes[2].z = { 0.0, -0.0 };
            lanes[3].z = { -0.0, -0.0 };
            Complex vectorResults[4];
            Complex jitResults[4];
            formula::ExpressionJit4 signJit;
            if (!signProgram.evaluate4(lanes, vectorResults)) {
                ++failures;
                continue;
            }
            if (!signJit.compile(signProgram) || !signJit.evaluate(lanes, jitResults)) {
                ++failures;
                continue;
            }
            for (int lane = 0; lane < 4; ++lane) {
                Complex scalar = signProgram.evaluate(lanes[lane]);
                if (!sameDoubleBits(scalar.real(), vectorResults[lane].real()) ||
                    !sameDoubleBits(scalar.imag(), vectorResults[lane].imag()))
                    ++failures;
                if (!sameDoubleBits(scalar.real(), jitResults[lane].real()) ||
                    !sameDoubleBits(scalar.imag(), jitResults[lane].imag()))
                    ++jitMismatch;
            }
            failures += jitMismatch;
        }
        formula::ExpressionJit4 unsupportedJit;
        if (unsupportedJit.compile(functions)) ++failures;
        {
            formula::ExpressionJit4 resetJit;
            if (!resetJit.compile(vectorProgram) || !resetJit.valid()) {
                ++failures;
            } else {
                resetJit.reset();
                if (resetJit.valid() || resetJit.codeAddress() != nullptr) ++failures;
            }
        }

        if (vectorJit.valid()) {
            ExpressionContext lanes[4] = { context, context, context, context };
            Complex outputs[4];
            const int iterations = 300000;
            volatile double sink = 0.0;
            auto begin = Clock::now();
            for (int i = 0; i < iterations; ++i) {
                lanes[0].iteration = i;
                vectorProgram.evaluate4(lanes, outputs);
                sink += outputs[0].real();
            }
            avxMs = std::chrono::duration<double, std::milli>(
                Clock::now() - begin).count();
            begin = Clock::now();
            for (int i = 0; i < iterations; ++i) {
                lanes[0].iteration = i;
                vectorJit.evaluate(lanes, outputs);
                sink += outputs[0].real();
            }
            jitMs = std::chrono::duration<double, std::milli>(
                Clock::now() - begin).count();
            formula::ExpressionJitInput4 input;
            formula::ExpressionJitOutput4 output;
            input.setContexts(lanes);
            begin = Clock::now();
            for (int i = 0; i < iterations; ++i) {
                for (int lane = 0; lane < 4; ++lane)
                    input.vectors[formula::ExpressionJitInput4::N_RE][lane] = (double)i;
                vectorJit.evaluate(input, output);
                sink += output.re[0];
            }
            jitRawMs = std::chrono::duration<double, std::milli>(
                Clock::now() - begin).count();
            if (!std::isfinite(sink)) ++failures;
        }
    }

    // Canonical recurrence through the interpreter must classify the same pixels
    // as a hand-written quadratic Julia loop.
    int classificationMismatch = 0;
    for (int y = 0; y < 40; ++y) {
        for (int x = 0; x < 60; ++x) {
            Complex z0{ -2.0 + 4.0 * x / 59.0, -1.3 + 2.6 * y / 39.0 };
            Complex fixedC{ -0.8, 0.156 };
            ExpressionContext ec; ec.z = ec.z0 = z0; ec.c = fixedC;
            Complex handwritten = z0;
            bool exprEscape = false, handwrittenEscape = false;
            for (int n = 0; n < 1000; ++n) {
                ec.iteration = n;
                ec.z = quadratic.evaluate(ec);
                handwritten = handwritten * handwritten + fixedC;
                exprEscape = std::norm(ec.z) > 16.0;
                handwrittenEscape = std::norm(handwritten) > 16.0;
                if (exprEscape || handwrittenEscape) break;
            }
            if (exprEscape != handwrittenEscape) ++classificationMismatch;
        }
    }
    failures += classificationMismatch;

    // Exercise the actual generic render backend for both supported pixel bindings.
    int renderMismatch = 0;
    const int RW = 60, RH = 40, RMIT = 500;
    std::vector<float> rendered((size_t)RW * RH, EMPTYPIXEL);
    Mandel expressionRenderer(RW, RH, RMIT, 1, rendered.data());
    mpf_t centerRe, centerIm, renderScale;
    mpf_init_set_d(centerRe, -0.5); mpf_init_set_ui(centerIm, 0);
    mpf_init_set_ui(renderScale, 1);
    ExpressionContext fixed;
    fixed.z0 = { 0.0, 0.0 };
    if (!expressionRenderer.ComputeExpression(
            centerRe, centerIm, renderScale, quadratic, fixed,
            FormulaParameter::C, RMIT, 4.0)) {
        ++failures;
    } else {
        const double halfW = 2.0, halfH = halfW * RH / RW;
        const double dx = 2.0 * halfW / (RW - 1);
        const double dy = 2.0 * halfH / (RH - 1);
        for (int y = 0; y < RH; ++y) {
            for (int x = 0; x < RW; ++x) {
                Complex c{ -0.5 - halfW + dx * x, -halfH + dy * y };
                Complex z{};
                float expected = -2.0f;
                for (int n = 0; n < RMIT; ++n) {
                    z = z * z + c;
                    if (std::hypot(z.real(), z.imag()) > 4.0) {
                        expected = (float)(n + 1); break;
                    }
                }
                if (rendered[(size_t)y * RW + x] != expected) ++renderMismatch;
            }
        }
    }
    std::vector<float> specialized = rendered;
    ExpressionProgram genericQuadratic;
    if (!compile(genericQuadratic, "z*z+c+0") ||
        !expressionRenderer.ComputeExpression(
            centerRe, centerIm, renderScale, genericQuadratic, fixed,
            FormulaParameter::C, RMIT, 4.0)) {
        ++failures;
    } else {
        for (size_t i = 0; i < rendered.size(); ++i)
            if (rendered[i] != specialized[i]) ++renderMismatch;
    }
    ExpressionProgram genericCubic;
    if (!expressionRenderer.ComputeExpression(
            centerRe, centerIm, renderScale, cubicProduct, fixed,
            FormulaParameter::C, RMIT, 4.0)) {
        ++failures;
    } else {
        specialized = rendered;
        if (!compile(genericCubic, "z*z*z+c+0") ||
            !expressionRenderer.ComputeExpression(
                centerRe, centerIm, renderScale, genericCubic, fixed,
                FormulaParameter::C, RMIT, 4.0)) {
            ++failures;
        } else {
            for (size_t i = 0; i < rendered.size(); ++i)
                if (rendered[i] != specialized[i]) ++renderMismatch;
        }
    }

    fixed.c = { -0.8, 0.156 };
    mpf_set_ui(centerRe, 0);
    if (!expressionRenderer.ComputeExpression(
            centerRe, centerIm, renderScale, quadratic, fixed,
            FormulaParameter::InitialZ, RMIT, 4.0)) {
        ++failures;
    } else {
        const double halfW = 2.0, halfH = halfW * RH / RW;
        const double dx = 2.0 * halfW / (RW - 1);
        const double dy = 2.0 * halfH / (RH - 1);
        for (int y = 0; y < RH; ++y) {
            for (int x = 0; x < RW; ++x) {
                Complex z{ -halfW + dx * x, -halfH + dy * y };
                float expected = std::hypot(z.real(), z.imag()) > 4.0 ? 0.0f : -2.0f;
                if (expected < 0.0f) {
                    for (int n = 0; n < RMIT; ++n) {
                        z = z * z + fixed.c;
                        if (std::hypot(z.real(), z.imag()) > 4.0) {
                            expected = (float)(n + 1); break;
                        }
                    }
                }
                if (rendered[(size_t)y * RW + x] != expected) ++renderMismatch;
            }
        }
    }
    // Degree-aware Smooth/EDE must agree with the generic automatic derivative.
    fixed = {};
    mpf_set_d(centerRe, -0.5); mpf_set_ui(centerIm, 0); mpf_set_ui(renderScale, 1);
    auto checkExpressionColoring = [&](formula::ExpressionColoring coloring) {
        int mismatches = 0;
        if (!expressionRenderer.ComputeExpression(
                centerRe, centerIm, renderScale, cubicProduct, fixed,
                FormulaParameter::C, RMIT, 4.0, coloring))
            return RW * RH;
        const double halfW = 2.0, halfH = halfW * RH / RW;
        const double dx = 2.0 * halfW / (RW - 1);
        const double dy = 2.0 * halfH / (RH - 1);
        for (int y = 0; y < RH; ++y) {
            for (int x = 0; x < RW; ++x) {
                ExpressionContext dc;
                dc.z = dc.z0 = {};
                dc.c = { -0.5 - halfW + dx * x, -halfH + dy * y };
                Complex derivative{};
                float expected = -2.0f;
                for (int n = 0; n < RMIT; ++n) {
                    ExpressionDerivativeSeed seed;
                    seed.z = derivative;
                    seed.c = 1.0;
                    Complex next, nextDerivative;
                    if (!cubicProduct.evaluateWithDerivative(
                            dc, seed, next, nextDerivative)) {
                        expected = 0.0f;
                        break;
                    }
                    dc.z = next;
                    derivative = nextDerivative;
                    double magnitude = std::hypot(next.real(), next.imag());
                    if (!std::isfinite(magnitude) || magnitude > 4.0) {
                        if (coloring == formula::ExpressionColoring::Smooth &&
                            std::isfinite(magnitude)) {
                            expected = (float)(n + 1 -
                                std::log(std::log(magnitude)) / std::log(3.0));
                        } else if (coloring == formula::ExpressionColoring::Distance &&
                                   std::isfinite(magnitude)) {
                            double denominator = std::abs(derivative) * dx;
                            expected = denominator > 0.0 && std::isfinite(denominator)
                                ? (float)(magnitude * std::log(magnitude) / denominator)
                                : 0.0f;
                        } else {
                            expected = (float)(n + 1);
                        }
                        break;
                    }
                }
                float actual = rendered[(size_t)y * RW + x];
                bool sameClass = isInterior(actual) == isInterior(expected);
                double tolerance = coloring == formula::ExpressionColoring::Distance
                    ? 2e-3 * std::max(1.0, std::fabs((double)expected))
                    : 2e-4;
                if (!sameClass || (!isInterior(expected) &&
                                   std::fabs((double)actual - expected) > tolerance))
                    ++mismatches;
            }
        }
        return mismatches;
    };
    renderMismatch += checkExpressionColoring(formula::ExpressionColoring::Smooth);
    renderMismatch += checkExpressionColoring(formula::ExpressionColoring::Distance);
    fixed.c = { -0.8, 0.156 };
    mpf_set_ui(centerRe, 0); mpf_set_ui(centerIm, 0); mpf_set_d(renderScale, 0.4);
    const double z0HalfWidth = 5.0;
    const double z0HalfHeight = z0HalfWidth * RH / RW;
    const double firstMagnitude = std::hypot(z0HalfWidth, z0HalfHeight);
    if (!expressionRenderer.ComputeExpression(
            centerRe, centerIm, renderScale, cubicProduct, fixed,
            FormulaParameter::InitialZ, RMIT, 4.0,
            formula::ExpressionColoring::Smooth)) {
        ++failures;
    } else {
        float expected = (float)(-std::log(std::log(firstMagnitude)) /
                                 std::log(3.0));
        if (rendered[0] != expected) ++renderMismatch;
    }
    if (!expressionRenderer.ComputeExpression(
            centerRe, centerIm, renderScale, cubicProduct, fixed,
            FormulaParameter::InitialZ, RMIT, 4.0,
            formula::ExpressionColoring::Distance)) {
        ++failures;
    } else {
        double pixelDx = 2.0 * z0HalfWidth / (RW - 1);
        float expected = (float)(firstMagnitude * std::log(firstMagnitude) /
                                 pixelDx);
        if (rendered[0] != expected) ++renderMismatch;
    }
    // Coloring requests with bailout<1 must safely downgrade to raw counts.
    fixed = {};
    mpf_set_d(centerRe, -0.5); mpf_set_ui(centerIm, 0); mpf_set_ui(renderScale, 1);
    if (!expressionRenderer.ComputeExpression(
            centerRe, centerIm, renderScale, cubicProduct, fixed,
            FormulaParameter::C, 1, 0.5,
            formula::ExpressionColoring::Distance)) {
        ++failures;
    } else {
        for (float value : rendered)
            if (value >= 0.0f && value != 0.0f && value != 1.0f)
                ++renderMismatch;
    }

    // Extreme finite bailout values must use overflow-safe magnitude tests.
    ExpressionProgram identity;
    if (compile(identity, "z")) {
        fixed = {};
        fixed.z0 = { 1e308, 1e308 };
        mpf_set_ui(centerRe, 0); mpf_set_ui(centerIm, 0);
        if (!expressionRenderer.ComputeExpression(
                centerRe, centerIm, renderScale, identity, fixed,
                FormulaParameter::C, 1, 1e308) ||
            rendered[0] != 0.0f)
            ++failures;
    }
    expressionRenderer.SetHalt(true);
    if (expressionRenderer.ComputeExpression(
            centerRe, centerIm, renderScale, quadratic, fixed,
            FormulaParameter::C, RMIT, 4.0))
        ++failures;
    expressionRenderer.SetHalt(false);
    mpf_clears(centerRe, centerIm, renderScale, (mpf_ptr)0);
    failures += renderMismatch;

    printf("=== arbitrary expression core\n");
    printf("  quadratic instructions=%zu stack=%zu\n",
           quadratic.instructionCount(), quadratic.stackDepth());
    printf("  parallel failures=%d   classification mismatch=%d   render mismatch=%d\n",
           parallelFailures.load(), classificationMismatch, renderMismatch);
    printf("  evaluator AVX2/JIT-wrapper/JIT-raw: %.3f / %.3f / %.3f ms  raw speedup %.2fx\n",
           avxMs, jitMs, jitRawMs, jitRawMs > 0.0 ? avxMs / jitRawMs : 0.0);
    printf("  JIT mismatches=%d\n", jitMismatch);
    if (vectorJit.valid()) {
        MEMORY_BASIC_INFORMATION memory{};
        VirtualQuery(vectorJit.codeAddress(), &memory, sizeof(memory));
        printf("  JIT W^X dual-mapping=%d execute-protect=0x%lx\n",
               vectorJit.usesDualMapping() ? 1 : 0,
               (unsigned long)memory.Protect);
    }
    printf("  => %s\n\n", failures == 0 ? "PASS" : "CHECK (expression core failure)");
    return failures == 0 ? 0 : 1;
}

static int runExpressionOracleCase() {
    using formula::Complex;
    using formula::ExpressionContext;
    using formula::ExpressionError;
    using formula::ExpressionOracle;
    using formula::ExpressionOracleContext;
    using formula::ExpressionProgram;
    using formula::MpfrComplex;

    ExpressionContext direct;
    direct.z = { 0.7, 0.2 };
    direct.c = { 1.2, 0.3 };
    direct.z0 = { 0.9, -0.4 };
    direct.parameters[0] = { 0.6, 0.1 };
    direct.parameters[1] = { 0.8, 0.25 };
    direct.iteration = 7;

    ExpressionOracleContext hp256(256), hp512(512);
    auto configure = [](ExpressionOracleContext& context) {
        context.z.set("0.7", "0.2");
        context.c.set("1.2", "0.3");
        context.z0.set("0.9", "-0.4");
        context.parameters[0].set("0.6", "0.1");
        context.parameters[1].set("0.8", "0.25");
        context.iteration = 7;
    };
    configure(hp256);
    configure(hp512);

    const char* expressions[] = {
        "z*z+c+p0*n/10",
        "sin(z)+cos(c)+sinh(p0)-tanh(z0)",
        "exp(z/3)+log(c)+log10(p0)+sqrt(z0)",
        "pow(z,c)+conj(p0)+complex(real(z),imag(c))+polar(abs(p1),arg(p1))",
        "abs(z)+norm(c)+arg(p0)"
    };

    int failures = 0;
    double maxDoubleError = 0.0;
    double maxPrecisionDelta = 0.0;
    mpfr_t delta;
    mpfr_init2(delta, 512);
    for (const char* source : expressions) {
        ExpressionProgram program;
        ExpressionError compileError;
        if (!program.compile(source, &compileError)) {
            printf("  oracle compile failed @ %zu: %s [%s]\n",
                   compileError.position, compileError.message.c_str(), source);
            ++failures;
            continue;
        }
        MpfrComplex result256(256), result512(512);
        std::string error256, error512;
        bool ok256 = ExpressionOracle::evaluate(program, hp256, result256, &error256);
        bool ok512 = ExpressionOracle::evaluate(program, hp512, result512, &error512);
        if (!ok256 || !ok512) {
            printf("  oracle domain failure: %s / %s [%s]\n",
                   error256.c_str(), error512.c_str(), source);
            ++failures;
            continue;
        }

        mpfr_sub(delta, result512.re, result256.re, MPFR_RNDN);
        mpfr_abs(delta, delta, MPFR_RNDN);
        double reDelta = mpfr_get_d(delta, MPFR_RNDU);
        mpfr_sub(delta, result512.im, result256.im, MPFR_RNDN);
        mpfr_abs(delta, delta, MPFR_RNDN);
        double imDelta = mpfr_get_d(delta, MPFR_RNDU);
        maxPrecisionDelta = std::max({ maxPrecisionDelta, reDelta, imDelta });
        if (reDelta > 1e-60 || imDelta > 1e-60) ++failures;

        Complex expected = result512.toDouble();
        Complex actual = program.evaluate(direct);
        double relative = std::abs(actual - expected) / std::max(1.0, std::abs(expected));
        maxDoubleError = std::max(maxDoubleError, relative);
        if (!(relative <= 2e-13)) ++failures;
    }
    struct DomainCase { const char* source; };
    const DomainCase invalid[] = { { "log(0)" }, { "1/0" }, { "polar(-1,0)" } };
    int domainAccepted = 0;
    for (const DomainCase& item : invalid) {
        ExpressionProgram program;
        ExpressionError compileError;
        if (!program.compile(item.source, &compileError)) {
            ++failures;
            continue;
        }
        MpfrComplex result(256);
        std::string error;
        if (ExpressionOracle::evaluate(program, hp256, result, &error) ||
            !mpfr_nan_p(result.re) || !mpfr_nan_p(result.im) || error.empty())
            ++domainAccepted;
    }
    failures += domainAccepted;

    ExpressionProgram zeroPower;
    ExpressionError zeroPowerError;
    if (!zeroPower.compile("0^2 + 0^0", &zeroPowerError)) {
        ++failures;
    } else {
        MpfrComplex result(256);
        std::string error;
        if (!ExpressionOracle::evaluate(zeroPower, hp256, result, &error) ||
            std::abs(result.toDouble() - Complex{ 1.0, 0.0 }) > 1e-15)
            ++failures;
    }

    auto evaluateSpecial = [&](const char* source, ExpressionOracleContext& context,
                               Complex expected, double tolerance) {
        ExpressionProgram program;
        ExpressionError compileError;
        if (!program.compile(source, &compileError)) { ++failures; return; }
        MpfrComplex result(context.z.precision());
        std::string error;
        if (!ExpressionOracle::evaluate(program, context, result, &error)) {
            printf("  special domain failure: %s [%s]\n", error.c_str(), source);
            ++failures; return;
        }
        Complex actual = result.toDouble();
        if (std::abs(actual - expected) > tolerance * std::max(1.0, std::abs(expected))) {
            printf("  special mismatch [%s]: actual=(%.17g,%.17g) expected=(%.17g,%.17g)\n",
                   source, actual.real(), actual.imag(), expected.real(), expected.imag());
            ++failures;
        }
    };

    ExpressionOracleContext sqrtContext(256);
    sqrtContext.z.set("-1", "1e-40");
    evaluateSpecial("sqrt(z)", sqrtContext, { 5e-41, 1.0 }, 1e-14);

    ExpressionOracleContext tinyDivision(256);
    mpfr_set_ui_2exp(tinyDivision.z.re, 1, -600000000, MPFR_RNDN);
    mpfr_set_zero(tinyDivision.z.im, 0);
    mpfr_set_ui_2exp(tinyDivision.c.re, 1, -600000000, MPFR_RNDN);
    mpfr_set_zero(tinyDivision.c.im, 0);
    evaluateSpecial("z/c", tinyDivision, { 1.0, 0.0 }, 1e-15);

    ExpressionOracleContext hugeDivision(256);
    mpfr_set_ui(hugeDivision.z.re, 1, MPFR_RNDN);
    mpfr_set_exp(hugeDivision.z.re, mpfr_get_emax() - 1);
    mpfr_set(hugeDivision.z.im, hugeDivision.z.re, MPFR_RNDN);
    hugeDivision.c.set(hugeDivision.z);
    evaluateSpecial("z/c", hugeDivision, { 1.0, 0.0 }, 1e-15);

    ExpressionOracleContext tinyQuotient(256);
    mpfr_set_ui(tinyQuotient.z.re, 1, MPFR_RNDN);
    mpfr_set_exp(tinyQuotient.z.re, mpfr_get_emin());
    mpfr_set_zero(tinyQuotient.z.im, 0);
    tinyQuotient.c.set(0.5, 0.0);
    ExpressionProgram quotientProgram;
    ExpressionError quotientError;
    if (!quotientProgram.compile("z/c", &quotientError)) {
        ++failures;
    } else {
        MpfrComplex actual(256);
        std::string error;
        if (!ExpressionOracle::evaluate(quotientProgram, tinyQuotient, actual, &error)) {
            ++failures;
        } else {
            mpfr_div(delta, tinyQuotient.z.re, tinyQuotient.c.re, MPFR_RNDN);
            if (!mpfr_equal_p(actual.re, delta) || !mpfr_zero_p(actual.im))
                ++failures;
        }

        ExpressionOracleContext mixedQuotient(256);
        mixedQuotient.z.set(1.0, 0.0625);
        mpfr_set_ui(mixedQuotient.c.re, 1, MPFR_RNDN);
        mpfr_set_exp(mixedQuotient.c.re, mpfr_get_emax() - 1);
        mpfr_set_zero(mixedQuotient.c.im, 0);
        if (!quotientProgram.valid()) {
            ++failures;
        } else {
            MpfrComplex actual(256);
            std::string error;
            if (!ExpressionOracle::evaluate(quotientProgram, mixedQuotient, actual, &error) ||
                !mpfr_number_p(actual.re) || mpfr_zero_p(actual.re) ||
                !mpfr_zero_p(actual.im))
                ++failures;
        }
    }

    ExpressionOracleContext hugeSqrt(256);
    mpfr_set_ui(hugeSqrt.z.re, 1, MPFR_RNDN);
    mpfr_set_exp(hugeSqrt.z.re, mpfr_get_emax() - 1);
    mpfr_set_zero(hugeSqrt.z.im, 0);
    ExpressionProgram hugeSqrtProgram;
    ExpressionError hugeSqrtError;
    if (!hugeSqrtProgram.compile("sqrt(z)", &hugeSqrtError)) {
        ++failures;
    } else {
        MpfrComplex actual(256);
        std::string error;
        if (!ExpressionOracle::evaluate(hugeSqrtProgram, hugeSqrt, actual, &error)) {
            printf("  huge sqrt domain failure: %s\n", error.c_str());
            ++failures;
        } else {
            mpfr_sqrt(delta, hugeSqrt.z.re, MPFR_RNDN);
            if (!mpfr_equal_p(actual.re, delta) || !mpfr_zero_p(actual.im)) {
                mpfr_t difference, ulp;
                mpfr_init2(difference, 256); mpfr_init2(ulp, 256);
                mpfr_sub(difference, actual.re, delta, MPFR_RNDN);
                mpfr_abs(difference, difference, MPFR_RNDN);
                mpfr_set_ui(ulp, 1, MPFR_RNDN);
                mpfr_exp_t ulpExponent = mpfr_get_exp(actual.re) -
                    (mpfr_exp_t)mpfr_get_prec(actual.re) + 1;
                mpfr_set_exp(ulp, ulpExponent);
                if (mpfr_cmp(difference, ulp) > 0 || !mpfr_zero_p(actual.im)) {
                    printf("  huge sqrt exceeds 1 ulp; diff-exp=%ld expected-exp=%ld\n",
                           mpfr_zero_p(difference) ? 0L : (long)mpfr_get_exp(difference),
                           (long)mpfr_get_exp(delta));
                    ++failures;
                }

                ExpressionOracleContext skewSqrt(256);
                mpfr_set_ui(skewSqrt.z.re, 1, MPFR_RNDN);
                mpfr_set_exp(skewSqrt.z.re, mpfr_get_emax() - 1);
                mpfr_set_d(skewSqrt.z.im, 0.25, MPFR_RNDN);
                if (!hugeSqrtProgram.valid()) {
                    ++failures;
                } else {
                    MpfrComplex actual(256);
                    std::string error;
                    if (!ExpressionOracle::evaluate(hugeSqrtProgram, skewSqrt, actual, &error) ||
                        mpfr_zero_p(actual.im) || !mpfr_number_p(actual.re) ||
                        !mpfr_number_p(actual.im))
                        ++failures;
                }
                mpfr_clear(difference); mpfr_clear(ulp);
            }
        }
    }

    ExpressionOracleContext largeTan(256);
    largeTan.z.set(1.0, 1e9);
    evaluateSpecial("tan(z)", largeTan, { 0.0, 1.0 }, 1e-14);
    largeTan.z.set(1e9, 1.0);
    evaluateSpecial("tanh(z)", largeTan, { 1.0, 0.0 }, 1e-14);
    ExpressionOracleContext tinyTan(2048);
    mpfr_set_zero(tinyTan.z.re, 0);
    mpfr_set_ui_2exp(tinyTan.z.im, 1, -1000, MPFR_RNDN);
    ExpressionProgram tinyIdentityProgram, tinyTanProgram, tinyTanhProgram;
    ExpressionError tinyError;
    if (!tinyIdentityProgram.compile("z", &tinyError) ||
        !tinyTanProgram.compile("tan(z)", &tinyError) ||
        !tinyTanhProgram.compile("tanh(z)", &tinyError)) {
        ++failures;
    } else {
        MpfrComplex identityResult(2048);
        std::string identityError;
        if (!ExpressionOracle::evaluate(tinyIdentityProgram, tinyTan, identityResult,
                                        &identityError))
            ++failures;
        MpfrComplex tanResult(2048), tanhResult(2048);
        std::string error;
        if (!ExpressionOracle::evaluate(tinyTanProgram, tinyTan, tanResult, &error)) {
            printf("  tiny tan domain failure: %s\n", error.c_str());
            ++failures;
        }
        mpfr_set_ui_2exp(delta, 1, -1000, MPFR_RNDN);
        mpfr_sub(delta, tanResult.im, delta, MPFR_RNDN);
        mpfr_abs(delta, delta, MPFR_RNDN);
        if (mpfr_cmp_ui_2exp(delta, 1, -1900) > 0) {
            printf("  tiny tan error exponent=%ld\n", mpfr_get_exp(delta));
            ++failures;
        }
        mpfr_set_ui_2exp(tinyTan.z.re, 1, -1000, MPFR_RNDN);
        mpfr_set_zero(tinyTan.z.im, 0);
        if (!ExpressionOracle::evaluate(tinyTanhProgram, tinyTan, tanhResult, &error)) {
            printf("  tiny tanh domain failure: %s\n", error.c_str());
            ++failures;
        }
        mpfr_set_ui_2exp(delta, 1, -1000, MPFR_RNDN);
        mpfr_sub(delta, tanhResult.re, delta, MPFR_RNDN);
        mpfr_abs(delta, delta, MPFR_RNDN);
        if (mpfr_cmp_ui_2exp(delta, 1, -1900) > 0) {
            printf("  tiny tanh error exponent=%ld\n", mpfr_get_exp(delta));
            ++failures;
        }
    }
    mpfr_clear(delta);

    struct OrbitStats {
        int classMismatch = 0;
        int maxIterationDiff = 0;
        int oracleDomainErrors = 0;
    };
    auto compareOrbitGrid = [&](const char* source, Complex parameter0,
                                int width, int height, int maxIterations) {
        OrbitStats stats;
        ExpressionProgram program;
        ExpressionError compileError;
        if (!program.compile(source, &compileError)) {
            ++stats.oracleDomainErrors;
            return stats;
        }
        mpfr_t magnitude;
        mpfr_init2(magnitude, 256);
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                Complex pixel{ -2.0 + 4.0 * x / (width - 1),
                               -1.3 + 2.6 * y / (height - 1) };
                ExpressionContext dc;
                dc.z = dc.z0 = {};
                dc.c = pixel;
                dc.parameters[0] = parameter0;
                bool doubleEscaped = false;
                int doubleIteration = maxIterations;
                for (int n = 0; n < maxIterations; ++n) {
                    dc.iteration = n;
                    dc.z = program.evaluate(dc);
                    if (!std::isfinite(dc.z.real()) || !std::isfinite(dc.z.imag()) ||
                        std::hypot(dc.z.real(), dc.z.imag()) > 4.0) {
                        doubleEscaped = true;
                        doubleIteration = n + 1;
                        break;
                    }
                }

                ExpressionOracleContext hc(256);
                hc.z.set(0.0); hc.z0.set(0.0);
                hc.c.set(pixel.real(), pixel.imag());
                hc.parameters[0].set(parameter0.real(), parameter0.imag());
                MpfrComplex next(256);
                bool oracleEscaped = false;
                int oracleIteration = maxIterations;
                for (int n = 0; n < maxIterations; ++n) {
                    hc.iteration = n;
                    std::string error;
                    if (!ExpressionOracle::evaluate(program, hc, next, &error)) {
                        ++stats.oracleDomainErrors;
                        oracleEscaped = true;
                        oracleIteration = n + 1;
                        break;
                    }
                    hc.z.set(next);
                    mpfr_hypot(magnitude, hc.z.re, hc.z.im, MPFR_RNDN);
                    if (mpfr_cmp_d(magnitude, 4.0) > 0) {
                        oracleEscaped = true;
                        oracleIteration = n + 1;
                        break;
                    }
                }
                if (doubleEscaped != oracleEscaped) ++stats.classMismatch;
                if (doubleEscaped && oracleEscaped)
                    stats.maxIterationDiff = std::max(
                        stats.maxIterationDiff,
                        std::abs(doubleIteration - oracleIteration));
            }
        }
        mpfr_clear(magnitude);
        return stats;
    };

    OrbitStats quadraticOrbit = compareOrbitGrid("z*z+c", {}, 30, 20, 300);
    OrbitStats parameterOrbit =
        compareOrbitGrid("z*z+c+p0*z*z*z", { 0.03, -0.02 }, 24, 16, 200);
    failures += quadraticOrbit.classMismatch + quadraticOrbit.oracleDomainErrors;
    failures += parameterOrbit.classMismatch + parameterOrbit.oracleDomainErrors;
    if (quadraticOrbit.maxIterationDiff > 0 || parameterOrbit.maxIterationDiff > 1)
        ++failures;

    printf("=== MPFR expression oracle\n");
    printf("  expressions=%zu max double relative error=%.6g\n",
           sizeof(expressions) / sizeof(expressions[0]), maxDoubleError);
    printf("  256-vs-512 max delta=%.6g domain accepted=%d\n",
           maxPrecisionDelta, domainAccepted);
    printf("  orbit quadratic mismatch=%d iter-diff=%d; parameter mismatch=%d iter-diff=%d\n",
           quadraticOrbit.classMismatch, quadraticOrbit.maxIterationDiff,
           parameterOrbit.classMismatch, parameterOrbit.maxIterationDiff);
    printf("  => %s\n\n", failures == 0 ? "PASS" : "CHECK (MPFR oracle failure)");
    return failures == 0 ? 0 : 1;
}

struct FormulaRegressionCase {
    const char* name;
    const char* source;
    FormulaParameter pixel;
    double centerRe;
    double centerIm;
    double scale;
    formula::Complex fixedZ0;
    formula::Complex fixedC;
    std::array<formula::Complex, 8> parameters{};
    double bailout = 4.0;
    int mxit = 500;
    int width = 64;
    int height = 44;
    const char* goldenPath = nullptr;
    int maxClassMismatch = 0;
    double maxDiff = 0.0;
    double maxMeanDiff = 0.0;
    double maxP99Diff = 0.0;
};

static bool renderExpressionOracle(const FormulaRegressionCase& test,
                                   const formula::ExpressionProgram& program,
                                   std::vector<float>& output) {
    using formula::ExpressionOracle;
    using formula::ExpressionOracleContext;
    using formula::MpfrComplex;
    output.assign((size_t)test.width * test.height, -2.0f);
    mpfr_prec_t precision = 512;
    if (const char* value = getenv("MANDEL_FORMULA_ORACLE_BITS"))
        precision = std::max<mpfr_prec_t>(128, (mpfr_prec_t)atol(value));
    const double halfWidth = 2.0 / test.scale;
    const double halfHeight = halfWidth * test.height / test.width;
    const double dx = 2.0 * halfWidth / (test.width - 1);
    const double dy = 2.0 * halfHeight / (test.height - 1);
    mpfr_t magnitude;
    mpfr_init2(magnitude, precision);
    for (int y = 0; y < test.height; ++y) {
        for (int x = 0; x < test.width; ++x) {
            formula::Complex pixel{
                test.centerRe - halfWidth + dx * x,
                test.centerIm - halfHeight + dy * y
            };
            ExpressionOracleContext context(precision);
            context.z0.set(test.fixedZ0.real(), test.fixedZ0.imag());
            context.c.set(test.fixedC.real(), test.fixedC.imag());
            for (int p = 0; p < 8; ++p)
                context.parameters[p].set(test.parameters[p].real(),
                                          test.parameters[p].imag());
            if (test.pixel == FormulaParameter::C)
                context.c.set(pixel.real(), pixel.imag());
            else
                context.z0.set(pixel.real(), pixel.imag());
            context.z.set(context.z0);

            mpfr_hypot(magnitude, context.z.re, context.z.im, MPFR_RNDN);
            if (mpfr_cmp_d(magnitude, test.bailout) > 0) {
                output[(size_t)y * test.width + x] = 0.0f;
                continue;
            }
            MpfrComplex next(precision);
            for (int n = 0; n < test.mxit; ++n) {
                context.iteration = n;
                std::string error;
                if (!ExpressionOracle::evaluate(program, context, next, &error)) {
                    fprintf(stderr,
                            "Formula oracle domain error in %s @ pixel %d,%d iter %d: %s\n",
                            test.name, x, y, n, error.c_str());
                    mpfr_clear(magnitude);
                    return false;
                }
                context.z.set(next);
                mpfr_hypot(magnitude, context.z.re, context.z.im, MPFR_RNDN);
                if (mpfr_cmp_d(magnitude, test.bailout) > 0) {
                    output[(size_t)y * test.width + x] = (float)(n + 1);
                    break;
                }
            }
        }
    }
    mpfr_clear(magnitude);
    return true;
}

static int runFormulaRegressionCase(const FormulaRegressionCase& test,
                                    bool updateGoldens) {
    formula::ExpressionProgram program;
    formula::ExpressionError compileError;
    if (!program.compile(test.source, &compileError)) {
        printf("=== formula %s\n  compile error @ %zu: %s\n  => CHECK\n\n",
               test.name, compileError.position, compileError.message.c_str());
        return 1;
    }

    formula::ExpressionContext fixed;
    fixed.z0 = test.fixedZ0;
    fixed.c = test.fixedC;
    fixed.parameters = test.parameters;
    std::vector<float> engine((size_t)test.width * test.height, EMPTYPIXEL);
    std::vector<float> golden(engine.size(), EMPTYPIXEL);
    Mandel renderer(test.width, test.height, test.mxit, 1, engine.data());
    mpf_t centerRe, centerIm, scale;
    mpf_init_set_d(centerRe, test.centerRe);
    mpf_init_set_d(centerIm, test.centerIm);
    mpf_init_set_d(scale, test.scale);
    auto begin = Clock::now();
    bool rendered = renderer.ComputeExpression(
        centerRe, centerIm, scale, program, fixed, test.pixel,
        test.mxit, test.bailout);
    double engineTime = since(begin);
    mpf_clears(centerRe, centerIm, scale, (mpf_ptr)0);
    if (!rendered) return 1;

    double oracleTime = 0.0;
    if (updateGoldens) {
        begin = Clock::now();
        if (!renderExpressionOracle(test, program, golden)) return 2;
        oracleTime = since(begin);
        if (!saveGolden(test.goldenPath, golden.data(), golden.size())) return 2;
    } else if (!loadGolden(test.goldenPath, golden.data(), golden.size())) {
        return 2;
    }

    int classMismatch = 0, exteriorBoth = 0, emptyPixels = 0;
    double maxDiff = 0.0, sumDiff = 0.0;
    std::vector<double> differences;
    for (size_t i = 0; i < engine.size(); ++i) {
        if (engine[i] == EMPTYPIXEL || golden[i] == EMPTYPIXEL) {
            ++emptyPixels;
            continue;
        }
        bool engineInterior = isInterior(engine[i]);
        bool goldenInterior = isInterior(golden[i]);
        if (engineInterior != goldenInterior) {
            ++classMismatch;
            continue;
        }
        if (!engineInterior) {
            ++exteriorBoth;
            double difference = std::fabs((double)engine[i] - golden[i]);
            maxDiff = std::max(maxDiff, difference);
            sumDiff += difference;
            differences.push_back(difference);
        }
    }
    double meanDiff = exteriorBoth ? sumDiff / exteriorBoth : 0.0;
    double p99Diff = 0.0;
    if (!differences.empty()) {
        size_t index = (size_t)std::ceil(0.99 * differences.size()) - 1;
        std::nth_element(differences.begin(), differences.begin() + index,
                         differences.end());
        p99Diff = differences[index];
    }
    bool ok = emptyPixels == 0 &&
              classMismatch <= test.maxClassMismatch &&
              maxDiff <= test.maxDiff &&
              meanDiff <= test.maxMeanDiff && p99Diff <= test.maxP99Diff;
    printf("=== formula %s  (%dx%d, mxit=%d, %s-plane)\n",
           test.name, test.width, test.height, test.mxit,
           test.pixel == FormulaParameter::C ? "c" : "z0");
    printf("  expression: %s\n", test.source);
    printf("  engine %.3fs%s", engineTime,
           updateGoldens ? "  oracle " : "  golden load");
    if (updateGoldens) printf("%.3fs", oracleTime);
    printf("\n  empty=%d class mismatch=%d  max=%.6g mean=%.6g p99=%.6g"
           "  checksum=0x%08x\n",
           emptyPixels, classMismatch, maxDiff, meanDiff, p99Diff,
           checksum(engine.data(), (int)engine.size()));
    printf("  => %s\n\n", ok ? "PASS" : "CHECK (formula golden mismatch)");
    return ok ? 0 : 1;
}

static int runFormulaRegressionSuite() {
    std::vector<FormulaRegressionCase> cases;
    FormulaRegressionCase quadratic{
        "quadratic-c", "z*z+c", FormulaParameter::C,
        -0.5, 0.0, 1.0, {}, {}, {}, 4.0, 500, 64, 44,
        "tests/golden/formula_quadratic_c.f32"
    };
    cases.push_back(quadratic);
    FormulaRegressionCase julia = quadratic;
    julia.name = "quadratic-z0";
    julia.pixel = FormulaParameter::InitialZ;
    julia.centerRe = 0.0;
    julia.fixedC = { -0.8, 0.156 };
    julia.goldenPath = "tests/golden/formula_quadratic_z0.f32";
    cases.push_back(julia);
    FormulaRegressionCase cubic = quadratic;
    cubic.name = "cubic-c";
    cubic.source = "z*z*z+c";
    cubic.centerRe = 0.0;
    cubic.goldenPath = "tests/golden/formula_cubic_c.f32";
    cases.push_back(cubic);
    FormulaRegressionCase burning = quadratic;
    burning.name = "burning-ship";
    burning.source = "sqr(complex(abs(re(z)),abs(im(z))))+c";
    burning.centerIm = -0.5;
    burning.mxit = 300;
    burning.goldenPath = "tests/golden/formula_burning_ship.f32";
    burning.maxClassMismatch = 11;
    burning.maxDiff = 100.0;
    burning.maxMeanDiff = 0.25;
    burning.maxP99Diff = 0.0;
    cases.push_back(burning);
    FormulaRegressionCase sine = quadratic;
    sine.name = "sine-c";
    sine.source = "sin(z)+c";
    sine.centerRe = 0.0;
    sine.bailout = 8.0;
    sine.mxit = 200;
    sine.goldenPath = "tests/golden/formula_sine_c.f32";
    cases.push_back(sine);
    FormulaRegressionCase parameter = quadratic;
    parameter.name = "parameter-polynomial";
    parameter.source = "z*z+c+p0*z";
    parameter.parameters[0] = { 0.15, -0.05 };
    parameter.goldenPath = "tests/golden/formula_parameter_poly.f32";
    cases.push_back(parameter);
    FormulaRegressionCase branch = quadratic;
    branch.name = "branch-power";
    branch.source = "exp(p0*log(z))+c";
    branch.fixedZ0 = { 0.3, 0.2 };
    branch.parameters[0] = { 2.5, 0.0 };
    branch.mxit = 150;
    branch.goldenPath = "tests/golden/formula_branch_power.f32";
    cases.push_back(branch);
    FormulaRegressionCase iteration = quadratic;
    iteration.name = "iteration-dependent";
    iteration.source = "z*z+c+0.0001*n";
    iteration.mxit = 200;
    iteration.goldenPath = "tests/golden/formula_iteration.f32";
    cases.push_back(iteration);

    bool update = getenv("MANDEL_UPDATE_FORMULA_GOLDENS") != nullptr;
    int result = 0;
    for (const FormulaRegressionCase& test : cases)
        result |= runFormulaRegressionCase(test, update);
    return result;
}

static std::string integerPowerSource(int power) {
    std::string source = "z";
    for (int i = 1; i < power; ++i) source += "*z";
    source += "+c";
    return source;
}

static bool renderIntegerPowerFrame(
        int power, FormulaParameter pixel,
        formula::ExpressionColoring coloring,
        int width, int height, int mxit,
        bool simd, std::vector<float>& output,
        double* elapsed = nullptr) {
    formula::ExpressionProgram program;
    formula::ExpressionError error;
    if (!program.compile(integerPowerSource(power), &error) ||
        program.fastIntegerPower() != power)
        return false;

    formula::ExpressionContext fixed;
    if (pixel == FormulaParameter::InitialZ)
        fixed.c = { -0.35, 0.62 };
    output.assign((size_t)width * height, EMPTYPIXEL);
    Mandel renderer(width, height, mxit, 1, output.data());
    mpf_t centerRe, centerIm, scale;
    mpf_init_set_d(centerRe, power == 2 && pixel == FormulaParameter::C
        ? -0.5 : 0.0);
    mpf_init_set_ui(centerIm, 0);
    mpf_init_set_ui(scale, 1);
    _putenv_s("MANDEL_EXPR_POWER_SIMD", simd ? "1" : "0");
    auto begin = Clock::now();
    bool okay = renderer.ComputeExpression(
        centerRe, centerIm, scale, program, fixed, pixel,
        mxit, 4.0, coloring);
    if (elapsed) *elapsed = since(begin);
    mpf_clears(centerRe, centerIm, scale, (mpf_ptr)0);
    return okay;
}

static int runMultibrotCase() {
    constexpr int W = 144, H = 96, MXIT = 1200;
    int failures = 0;
    int totalBitMismatches = 0;
    int totalClassMismatches = 0;
    double maximumDifference = 0.0;
    const formula::ExpressionColoring colorings[] = {
        formula::ExpressionColoring::Raw,
        formula::ExpressionColoring::Smooth,
        formula::ExpressionColoring::Distance
    };
    const char* coloringNames[] = { "raw", "smooth", "distance" };

    printf("=== integer Multibrot AVX2 correctness\n");
    for (int power = 2; power <= 8; ++power) {
        for (FormulaParameter pixel :
             { FormulaParameter::C, FormulaParameter::InitialZ }) {
            for (int coloringIndex = 0; coloringIndex < 3; ++coloringIndex) {
                std::vector<float> scalar;
                std::vector<float> simd;
                if (!renderIntegerPowerFrame(
                        power, pixel, colorings[coloringIndex],
                        W, H, MXIT, false, scalar) ||
                    !renderIntegerPowerFrame(
                        power, pixel, colorings[coloringIndex],
                        W, H, MXIT, true, simd)) {
                    ++failures;
                    continue;
                }
                int bitMismatches = 0;
                int classMismatches = 0;
                double maxDifference = 0.0;
                for (size_t i = 0; i < scalar.size(); ++i) {
                    if (std::memcmp(&scalar[i], &simd[i], sizeof(float)) != 0)
                        ++bitMismatches;
                    if (isInterior(scalar[i]) != isInterior(simd[i]))
                        ++classMismatches;
                    else if (!isInterior(scalar[i]))
                        maxDifference = std::max(
                            maxDifference,
                            std::fabs((double)scalar[i] - simd[i]));
                }
                totalBitMismatches += bitMismatches;
                totalClassMismatches += classMismatches;
                maximumDifference =
                    std::max(maximumDifference, maxDifference);
                if (bitMismatches || classMismatches) ++failures;
                printf("  z^%d+c %-2s-plane %-8s bits=%d class=%d max=%.6g\n",
                       power,
                       pixel == FormulaParameter::C ? "c" : "z0",
                       coloringNames[coloringIndex],
                       bitMismatches, classMismatches, maxDifference);
            }
        }
    }

    constexpr int OW = 48, OH = 32, OMIT = 300;
    printf("  MPFR raw-orbit gates:\n");
    for (int power = 3; power <= 8; ++power) {
        std::string source = integerPowerSource(power);
        formula::ExpressionProgram program;
        formula::ExpressionError error;
        if (!program.compile(source, &error)) {
            ++failures;
            continue;
        }
        for (FormulaParameter pixel :
             { FormulaParameter::C, FormulaParameter::InitialZ }) {
            FormulaRegressionCase test{};
            test.name = "multibrot-oracle";
            test.source = source.c_str();
            test.pixel = pixel;
            test.centerRe = 0.0;
            test.centerIm = 0.0;
            test.scale = 1.0;
            if (pixel == FormulaParameter::InitialZ)
                test.fixedC = { -0.35, 0.62 };
            test.bailout = 4.0;
            test.mxit = OMIT;
            test.width = OW;
            test.height = OH;
            std::vector<float> oracle;
            std::vector<float> simd;
            if (!renderExpressionOracle(test, program, oracle) ||
                !renderIntegerPowerFrame(
                    power, pixel, formula::ExpressionColoring::Raw,
                    OW, OH, OMIT, true, simd)) {
                ++failures;
                continue;
            }
            int classMismatches = 0;
            int iterationMismatches = 0;
            for (size_t i = 0; i < oracle.size(); ++i) {
                bool oracleInterior = isInterior(oracle[i]);
                bool simdInterior = isInterior(simd[i]);
                if (oracleInterior != simdInterior)
                    ++classMismatches;
                else if (!oracleInterior && oracle[i] != simd[i])
                    ++iterationMismatches;
            }
            if (classMismatches || iterationMismatches) ++failures;
            printf("    z^%d+c %-2s-plane class=%d iteration=%d\n",
                   power,
                   pixel == FormulaParameter::C ? "c" : "z0",
                   classMismatches, iterationMismatches);
        }
    }

    auto extremeParity = [&](const char* name, const char* centerText,
                             const char* scaleText, double bailout,
                             formula::ExpressionColoring coloring,
                             int mxit) {
        formula::ExpressionProgram program;
        formula::ExpressionError error;
        if (!program.compile("z*z+c", &error)) {
            ++failures;
            return;
        }
        constexpr int EW = 4, EH = 3;
        std::vector<float> scalar((size_t)EW * EH, EMPTYPIXEL);
        std::vector<float> simd((size_t)EW * EH, EMPTYPIXEL);
        formula::ExpressionContext fixed;
        mpf_t centerRe, centerIm, scale;
        mpf_init_set_str(centerRe, centerText, 10);
        mpf_init_set_ui(centerIm, 0);
        mpf_init_set_str(scale, scaleText, 10);
        _putenv_s("MANDEL_EXPR_POWER_SIMD", "0");
        Mandel scalarRenderer(EW, EH, mxit, 1, scalar.data());
        bool scalarOkay = scalarRenderer.ComputeExpression(
            centerRe, centerIm, scale, program, fixed,
            FormulaParameter::C, mxit, bailout, coloring);
        _putenv_s("MANDEL_EXPR_POWER_SIMD", "1");
        Mandel simdRenderer(EW, EH, mxit, 1, simd.data());
        bool simdOkay = simdRenderer.ComputeExpression(
            centerRe, centerIm, scale, program, fixed,
            FormulaParameter::C, mxit, bailout, coloring);
        mpf_clears(centerRe, centerIm, scale, (mpf_ptr)0);
        int mismatches = scalar.size() == simd.size()
            ? (int)!std::equal(
                scalar.begin(), scalar.end(), simd.begin(),
                [](float a, float b) {
                    return std::memcmp(&a, &b, sizeof(float)) == 0;
                })
            : 1;
        if (!scalarOkay || !simdOkay || mismatches) ++failures;
        printf("  extreme %-18s parity=%s\n", name,
               scalarOkay && simdOkay && !mismatches ? "exact" : "FAIL");
    };
    extremeParity(
        "huge bailout", "1.5e200", "1", 1.0e200,
        formula::ExpressionColoring::Raw, 1);
    extremeParity(
        "tiny bailout", "0", "1e200", 1.0e-200,
        formula::ExpressionColoring::Raw, 8);
    extremeParity(
        "subnormal bailout", "1.0000000000000002e-160",
        "1e170", 1.0e-160,
        formula::ExpressionColoring::Raw, 1);
    extremeParity(
        "subnormal dx EDE", "2", "1e309", 4.0,
        formula::ExpressionColoring::Distance, 8);

    constexpr int BW = 960, BH = 540, BMIT = 3000, SAMPLES = 5;
    std::vector<double> scalarTimes;
    std::vector<double> simdTimes;
    std::vector<float> benchmarkOutput;
    for (int sample = 0; sample < SAMPLES; ++sample) {
        double scalar = 0.0, simd = 0.0;
        if (!renderIntegerPowerFrame(
                3, FormulaParameter::C,
                formula::ExpressionColoring::Smooth,
                BW, BH, BMIT, false, benchmarkOutput, &scalar) ||
            !renderIntegerPowerFrame(
                3, FormulaParameter::C,
                formula::ExpressionColoring::Smooth,
                BW, BH, BMIT, true, benchmarkOutput, &simd)) {
            ++failures;
            break;
        }
        scalarTimes.push_back(scalar);
        simdTimes.push_back(simd);
    }
    _putenv_s("MANDEL_EXPR_POWER_SIMD", "");
    std::sort(scalarTimes.begin(), scalarTimes.end());
    std::sort(simdTimes.begin(), simdTimes.end());
    double scalarMedian = scalarTimes.empty()
        ? 0.0 : scalarTimes[scalarTimes.size() / 2];
    double simdMedian = simdTimes.empty()
        ? 0.0 : simdTimes[simdTimes.size() / 2];
    double speedup = simdMedian > 0.0 ? scalarMedian / simdMedian : 0.0;
    if (speedup < 1.05) ++failures;

    formula::ExpressionProgram cubicProgram;
    formula::ExpressionError cubicError;
    cubicProgram.compile("z*z*z+c", &cubicError);
    constexpr int RW = 320, RH = 216, RMIT = 1500, RSAMPLES = 3;
    std::vector<double> genericResidualTimes;
    std::vector<double> cubicResidualTimes;
    std::vector<float> genericResidual;
    std::vector<float> cubicResidual;
    auto renderResidual = [&](bool specialized,
                              formula::ExpressionColoring coloring,
                              std::vector<float>& output,
                              double& elapsed) {
        output.assign((size_t)RW * RH, EMPTYPIXEL);
        Mandel renderer(RW, RH, RMIT, 1, output.data());
        formula::ExpressionContext fixed;
        mpf_t centerRe, centerIm, scale;
        mpf_init_set_d(centerRe, 0.3849001794597505);
        mpf_init_set_ui(centerIm, 0);
        mpf_init_set_d(scale, 1e5);
        _putenv_s(
            "MANDEL_EXPR_RESIDUAL_POWER", specialized ? "1" : "0");
        bool used = false;
        auto begin = Clock::now();
        bool okay = renderer.ComputeExpressionResidual(
            centerRe, centerIm, scale, cubicProgram, fixed,
            FormulaParameter::C, RMIT, 4.0, &used, coloring);
        elapsed = since(begin);
        mpf_clears(centerRe, centerIm, scale, (mpf_ptr)0);
        return okay && used;
    };
    const formula::ExpressionColoring residualColorings[] = {
        formula::ExpressionColoring::Raw,
        formula::ExpressionColoring::Smooth,
        formula::ExpressionColoring::Distance
    };
    const char* residualColoringNames[] = {
        "raw", "smooth", "distance"
    };
    for (int coloringIndex = 0; coloringIndex < 3; ++coloringIndex) {
        std::vector<float> direct((size_t)RW * RH, EMPTYPIXEL);
        std::vector<float> residual;
        Mandel directRenderer(RW, RH, RMIT, 1, direct.data());
        formula::ExpressionContext fixed;
        mpf_t centerRe, centerIm, scale;
        mpf_init_set_d(centerRe, 0.3849001794597505);
        mpf_init_set_ui(centerIm, 0);
        mpf_init_set_d(scale, 1e5);
        _putenv_s("MANDEL_EXPR_POWER_SIMD", "0");
        bool directOkay = directRenderer.ComputeExpression(
            centerRe, centerIm, scale, cubicProgram, fixed,
            FormulaParameter::C, RMIT, 4.0,
            residualColorings[coloringIndex]);
        mpf_clears(centerRe, centerIm, scale, (mpf_ptr)0);
        double residualTime = 0.0;
        bool residualOkay = renderResidual(
            true, residualColorings[coloringIndex],
            residual, residualTime);
        int classMismatches = 0;
        int floorMismatches = 0;
        double maxDifference = 0.0;
        for (size_t i = 0; i < direct.size(); ++i) {
            bool directInterior = isInterior(direct[i]);
            bool residualInterior = isInterior(residual[i]);
            if (directInterior != residualInterior) {
                ++classMismatches;
            } else if (!directInterior) {
                if ((int)direct[i] != (int)residual[i])
                    ++floorMismatches;
                maxDifference = std::max(
                    maxDifference,
                    std::fabs((double)direct[i] - residual[i]));
            }
        }
        if (!directOkay || !residualOkay ||
            classMismatches || floorMismatches ||
            (coloringIndex == 0 && maxDifference != 0.0) ||
            (coloringIndex != 0 && maxDifference > 1e-3))
            ++failures;
        printf("  cubic residual %-8s class=%d floor=%d max=%.6g\n",
               residualColoringNames[coloringIndex],
               classMismatches, floorMismatches, maxDifference);
    }
    for (int sample = 0; sample < RSAMPLES; ++sample) {
        double genericTime = 0.0, specializedTime = 0.0;
        if (!renderResidual(
                false, formula::ExpressionColoring::Raw,
                genericResidual, genericTime) ||
            !renderResidual(
                true, formula::ExpressionColoring::Raw,
                cubicResidual, specializedTime)) {
            ++failures;
            break;
        }
        genericResidualTimes.push_back(genericTime);
        cubicResidualTimes.push_back(specializedTime);
    }
    _putenv_s("MANDEL_EXPR_RESIDUAL_POWER", "");
    int residualBitMismatches = 0;
    if (genericResidual.size() == cubicResidual.size()) {
        for (size_t i = 0; i < genericResidual.size(); ++i) {
            if (std::memcmp(
                    &genericResidual[i], &cubicResidual[i],
                    sizeof(float)) != 0)
                ++residualBitMismatches;
        }
    } else {
        residualBitMismatches = 1;
    }
    std::sort(genericResidualTimes.begin(), genericResidualTimes.end());
    std::sort(cubicResidualTimes.begin(), cubicResidualTimes.end());
    double genericResidualMedian = genericResidualTimes.empty()
        ? 0.0 : genericResidualTimes[genericResidualTimes.size() / 2];
    double cubicResidualMedian = cubicResidualTimes.empty()
        ? 0.0 : cubicResidualTimes[cubicResidualTimes.size() / 2];
    double residualSpeedup = cubicResidualMedian > 0.0
        ? genericResidualMedian / cubicResidualMedian : 0.0;
    if (residualBitMismatches || residualSpeedup < 1.10) ++failures;

    printf("  aggregate bits=%d class=%d max=%.6g\n",
           totalBitMismatches, totalClassMismatches, maximumDifference);
    printf("  cubic %dx%d mxit=%d scalar/AVX2 %.3f/%.3f s speedup %.2fx\n",
           BW, BH, BMIT, scalarMedian, simdMedian, speedup);
    printf("  cubic residual generic/specialized %.3f/%.3f s speedup %.2fx"
           " bits=%d\n",
           genericResidualMedian, cubicResidualMedian,
           residualSpeedup, residualBitMismatches);
    printf("  => %s\n\n",
           failures == 0 ? "PASS"
                         : "CHECK (Multibrot correctness/performance)");
    return failures == 0 ? 0 : 1;
}

static int runExpressionResidualSuite() {
    struct ResidualCase {
        const char* name;
        const char* source;
        FormulaParameter pixel;
        formula::Complex center;
        double scale;
        formula::Complex fixedZ0;
        formula::Complex fixedC;
        formula::Complex p0;
        int mxit;
        double bailout;
        bool expectResidual;
        formula::ExpressionColoring coloring =
            formula::ExpressionColoring::Raw;
    };
    const ResidualCase cases[] = {
        { "quadratic-c", "z*z+c", FormulaParameter::C,
          { -0.75, 0.0 }, 1e4, {}, {}, {}, 1000, 4.0, true },
        { "cubic-c", "z*z*z+c", FormulaParameter::C,
          { 0.3849001794597505, 0.0 }, 1e5,
          {}, {}, {}, 1500, 4.0, true },
        { "sine-c", "sin(z)+c", FormulaParameter::C,
          {}, 4.0, {}, {}, {}, 200, 8.0, true },
        { "burning-c", "sqr(complex(abs(re(z)),abs(im(z))))+c",
          FormulaParameter::C, {}, 4.0, {}, {}, {}, 200, 4.0, true },
        { "parameter-c", "z*z+c+p0*z", FormulaParameter::C,
          {}, 4.0, {}, {}, { 0.15, -0.05 }, 300, 4.0, true },
        { "quadratic-z0", "z*z+c", FormulaParameter::InitialZ,
          {}, 1.0, {}, {}, {}, 500, 4.0, true },
        { "escaping-reference-fallback", "z*z+c", FormulaParameter::C,
          { 0.5, 0.5 }, 10.0, {}, {}, {}, 300, 4.0, false }
        ,
        { "cubic-escaping-smooth", "z*z*z+c", FormulaParameter::C,
          { 1.0, 1.0 }, 10.0, {}, {}, {}, 300, 4.0, false,
          formula::ExpressionColoring::Smooth },
        { "cubic-small-bailout-distance", "z*z*z+c",
          FormulaParameter::C, {}, 4.0, {}, {}, {}, 100, 0.5, true,
          formula::ExpressionColoring::Distance }
    };

    int failures = 0;
    for (const ResidualCase& test : cases) {
        formula::ExpressionProgram program;
        formula::ExpressionError error;
        if (!program.compile(test.source, &error)) {
            ++failures;
            continue;
        }
        formula::ExpressionContext fixed;
        fixed.z0 = test.fixedZ0;
        fixed.c = test.fixedC;
        fixed.parameters[0] = test.p0;
        constexpr int W = 64, H = 44;
        std::vector<float> output((size_t)W * H, EMPTYPIXEL);
        std::vector<float> direct;
        Mandel renderer(W, H, test.mxit, 1, output.data());
        mpf_t centerRe, centerIm, scale;
        mpf_init_set_d(centerRe, test.center.real());
        mpf_init_set_d(centerIm, test.center.imag());
        mpf_init_set_d(scale, test.scale);
        auto begin = Clock::now();
        bool ok = renderer.ComputeExpression(
            centerRe, centerIm, scale, program, fixed, test.pixel,
            test.mxit, test.bailout, test.coloring);
        double directTime = since(begin);
        direct = output;
        bool usedResidual = false;
        begin = Clock::now();
        ok = ok && renderer.ComputeExpressionResidual(
            centerRe, centerIm, scale, program, fixed, test.pixel,
            test.mxit, test.bailout, &usedResidual, test.coloring);
        double residualTime = since(begin);
        mpf_clears(centerRe, centerIm, scale, (mpf_ptr)0);

        int classMismatch = 0, iterationMismatch = 0, empty = 0;
        int maxDifference = 0;
        for (size_t i = 0; i < output.size(); ++i) {
            if (direct[i] == EMPTYPIXEL || output[i] == EMPTYPIXEL) {
                ++empty;
                continue;
            }
            if (isInterior(direct[i]) != isInterior(output[i])) {
                ++classMismatch;
            } else if (!isInterior(direct[i])) {
                int difference = (int)std::fabs(direct[i] - output[i]);
                if (difference) ++iterationMismatch;
                maxDifference = std::max(maxDifference, difference);
            }
        }
        bool passed = ok && empty == 0 &&
                      usedResidual == test.expectResidual &&
                      classMismatch == 0 && iterationMismatch == 0 &&
                      maxDifference == 0;
        if (!passed) ++failures;
        printf("=== residual %s\n", test.name);
        printf("  direct/residual %.3f/%.3f s used=%d expected=%d\n",
               directTime, residualTime, usedResidual ? 1 : 0,
               test.expectResidual ? 1 : 0);
        printf("  empty=%d class mismatch=%d iteration mismatch=%d max=%d\n",
               empty, classMismatch, iterationMismatch, maxDifference);
        printf("  => %s\n\n", passed ? "PASS" : "CHECK (residual mismatch)");
    }
    return failures == 0 ? 0 : 1;
}

static int runBackendCase() {
    constexpr int W = 96, H = 64, MXIT = 500;
    mpf_set_default_prec(256);
    mpf_t centerRe, centerIm, scale, fixedCRe, fixedCIm;
    mpf_init_set_d(centerRe, -0.5); mpf_init_set_ui(centerIm, 0);
    mpf_init_set_ui(scale, 1);
    mpf_init_set_d(fixedCRe, -0.8); mpf_init_set_d(fixedCIm, 0.156);
    int failures = 0;

    auto same = [](const std::vector<float>& a, const std::vector<float>& b) {
        return a.size() == b.size() &&
               std::memcmp(a.data(), b.data(), a.size() * sizeof(float)) == 0;
    };
    std::unique_ptr<IComputeBackend> backend = createComputeBackend("cpu");
    if (!backend || backend->info().name != "CPU" ||
        backend->info().hardwareAccelerated || backend->info().fallback)
        ++failures;

    std::vector<float> direct((size_t)W * H), dispatched((size_t)W * H);
    Mandel directMandel(W, H, MXIT, 1, direct.data());
    Mandel backendMandel(W, H, MXIT, 1, dispatched.data());
    directMandel.Compute(centerRe, centerIm, scale, MXIT, 0);
    ComputeRequest request;
    request.mode = ComputeMode::Mandelbrot;
    request.cpuEngine = &backendMandel;
    request.centerRe = centerRe; request.centerIm = centerIm; request.scale = scale;
    request.width = W; request.height = H; request.sub = 1;
    request.maxIterations = MXIT; request.iterations = dispatched.data();
    backend->resetCancellation();
    if (!backend->compute(request) || !same(direct, dispatched)) ++failures;

    directMandel.ComputeJulia(centerRe, centerIm, scale, fixedCRe, fixedCIm,
                              MXIT, 0);
    request.mode = ComputeMode::Julia;
    request.fixedCRe = fixedCRe; request.fixedCIm = fixedCIm;
    backend->resetCancellation();
    if (!backend->compute(request) || !same(direct, dispatched)) ++failures;

    formula::ExpressionProgram expression;
    formula::ExpressionError compileError;
    formula::ExpressionContext fixed;
    if (!expression.compile("z*z+c+0", &compileError)) {
        ++failures;
    } else {
        directMandel.ComputeExpression(centerRe, centerIm, scale, expression,
                                       fixed, FormulaParameter::C, MXIT, 4.0);
        request.mode = ComputeMode::Expression;
        request.expression = &expression;
        request.expressionFixed = &fixed;
        request.expressionPixel = FormulaParameter::C;
        request.expressionBailout = 4.0;
        request.expressionColoring = formula::ExpressionColoring::Raw;
        request.fixedCRe = request.fixedCIm = nullptr;
        backend->resetCancellation();
        if (!backend->compute(request) || !same(direct, dispatched)) ++failures;
    }

    ComputeRequest invalid;
    if (backend->compute(invalid)) ++failures;
    ComputeRequest invalidMode = request;
    invalidMode.mode = (ComputeMode)99;
    if (backend->compute(invalidMode)) ++failures;
    ComputeRequest invalidSize = request;
    invalidSize.width = 1;
    if (backend->compute(invalidSize)) ++failures;
    mpf_t zeroScale;
    mpf_init_set_ui(zeroScale, 0);
    ComputeRequest invalidScale = request;
    invalidScale.scale = zeroScale;
    if (backend->compute(invalidScale)) ++failures;
    mpf_clear(zeroScale);
    backend->cancel();
    if (backend->compute(request)) ++failures;
    backend->resetCancellation();

    std::unique_ptr<IComputeBackend> fallback =
        createComputeBackend("not-a-backend");
    if (!fallback || !fallback->info().fallback ||
        fallback->info().hardwareAccelerated ||
        fallback->info().detail.find("unknown backend") == std::string::npos)
        ++failures;

    // WARP runs the exact D3D11 compute shader even on CI/VM hosts without a
    // physical GPU. GPU smooth values use 2xFP32 orbit arithmetic and FP32
    // transcendentals, so require exact classification and escape-iteration
    // floors rather than byte identity with the CPU's FP64 smooth values.
    std::unique_ptr<IComputeBackend> warp = createComputeBackend("warp");
    int gpuEmpty = 0, gpuClassMismatch = 0, gpuFloorMismatch = 0;
    double gpuMaxDifference = 0.0;
    int gpuWorst = -1;
    float gpuWorstCpu = 0.0f, gpuWorstGpu = 0.0f;
    int gpuStressClass = 0, gpuStressFloor = 0, gpuLayoutMismatch = 0;
    double gpuStressMaxDifference = 0.0;
    if (!warp || warp->info().fallback ||
        warp->info().name != "D3D11 WARP" ||
        warp->info().hardwareAccelerated) {
        ++failures;
    } else {
        std::vector<float> gpu((size_t)W * H, EMPTYPIXEL);
        Mandel gpuMandel(W, H, MXIT, 1, gpu.data());
        directMandel.Compute(centerRe, centerIm, scale, MXIT, 0);
        request.mode = ComputeMode::Mandelbrot;
        request.cpuEngine = &gpuMandel;
        request.iterations = gpu.data();
        request.expression = nullptr;
        request.expressionFixed = nullptr;
        warp->resetCancellation();
        bool gpuOk = warp->compute(request);
        if (!gpuOk || !warp->lastComputeUsedGpuPath()) ++failures;
        for (size_t i = 0; i < gpu.size(); ++i) {
            if (gpu[i] == EMPTYPIXEL) {
                ++gpuEmpty;
                continue;
            }
            bool cpuInterior = isInterior(direct[i]);
            bool gpuInterior = isInterior(gpu[i]);
            if (cpuInterior != gpuInterior) {
                ++gpuClassMismatch;
            } else if (!cpuInterior) {
                if ((int)direct[i] != (int)gpu[i]) ++gpuFloorMismatch;
                double difference =
                    std::fabs((double)direct[i] - (double)gpu[i]);
                if (difference > gpuMaxDifference) {
                    gpuMaxDifference = difference;
                    gpuWorst = static_cast<int>(i);
                    gpuWorstCpu = direct[i];
                    gpuWorstGpu = gpu[i];
                }
            }
        }
        if (gpuEmpty != 0 || gpuClassMismatch != 0 ||
            gpuFloorMismatch != 0 || gpuMaxDifference > 0.01)
            ++failures;

        // Precision-limit shallow view (the GPU gate is scale <= 1e6).
        mpf_set_str(centerRe, "-0.743643887037151", 10);
        mpf_set_str(centerIm, "0.13182590420533", 10);
        mpf_set_ui(scale, 1000000);
        directMandel.Compute(centerRe, centerIm, scale, MXIT, 0);
        std::fill(gpu.begin(), gpu.end(), EMPTYPIXEL);
        request.centerRe = centerRe;
        request.centerIm = centerIm;
        request.scale = scale;
        warp->resetCancellation();
        if (!warp->compute(request)) {
            ++failures;
        } else if (!warp->lastComputeUsedGpuPath()) {
            ++failures;
        } else {
            for (size_t i = 0; i < gpu.size(); ++i) {
                bool cpuInterior = isInterior(direct[i]);
                bool gpuInterior = isInterior(gpu[i]);
                if (cpuInterior != gpuInterior) {
                    ++gpuStressClass;
                } else if (!cpuInterior) {
                    if ((int)direct[i] != (int)gpu[i]) ++gpuStressFloor;
                    gpuStressMaxDifference = std::max(
                        gpuStressMaxDifference,
                        std::fabs((double)direct[i] - (double)gpu[i]));
                }
            }
            if (gpuStressClass != 0 || gpuStressFloor != 0 ||
                gpuStressMaxDifference > 0.01)
                ++failures;
        }

        // The GUI normally keeps an odd adaptive-AA backing grid even when AA is
        // off. The GPU writes only each pixel's center sample and must preserve
        // every unresolved subpixel as EMPTYPIXEL.
        constexpr int LW = 32, LH = 24, LSUB = 5;
        const size_t layoutCount = (size_t)LW * LH * LSUB * LSUB;
        std::vector<float> layoutCpu(layoutCount, EMPTYPIXEL);
        std::vector<float> layoutGpu(layoutCount, EMPTYPIXEL);
        Mandel layoutCpuMandel(LW, LH, MXIT, LSUB, layoutCpu.data());
        Mandel layoutGpuMandel(LW, LH, MXIT, LSUB, layoutGpu.data());
        mpf_set_d(centerRe, -0.5);
        mpf_set_ui(centerIm, 0);
        mpf_set_ui(scale, 1);
        layoutCpuMandel.Compute(centerRe, centerIm, scale, MXIT, 0);
        ComputeRequest layoutRequest;
        layoutRequest.mode = ComputeMode::Mandelbrot;
        layoutRequest.cpuEngine = &layoutGpuMandel;
        layoutRequest.centerRe = centerRe;
        layoutRequest.centerIm = centerIm;
        layoutRequest.scale = scale;
        layoutRequest.width = LW;
        layoutRequest.height = LH;
        layoutRequest.sub = LSUB;
        layoutRequest.maxIterations = MXIT;
        layoutRequest.iterations = layoutGpu.data();
        warp->resetCancellation();
        if (!warp->compute(layoutRequest)) {
            ++failures;
        } else if (!warp->lastComputeUsedGpuPath()) {
            ++failures;
        } else {
            for (size_t i = 0; i < layoutCount; ++i)
                if (layoutCpu[i] != layoutGpu[i] &&
                    !(layoutCpu[i] >= 0.0f && layoutGpu[i] >= 0.0f &&
                      (int)layoutCpu[i] == (int)layoutGpu[i] &&
                      std::fabs(layoutCpu[i] - layoutGpu[i]) <= 0.01f))
                    ++gpuLayoutMismatch;
            if (gpuLayoutMismatch != 0) ++failures;
        }

        // Unsupported modes must remain byte-identical through the owned CPU
        // fallback instead of returning a partially rendered GPU frame.
        directMandel.ComputeJulia(
            centerRe, centerIm, scale, fixedCRe, fixedCIm, MXIT, 0);
        request.mode = ComputeMode::Julia;
        request.fixedCRe = fixedCRe;
        request.fixedCIm = fixedCIm;
        warp->resetCancellation();
        if (!warp->compute(request) || !same(direct, gpu) ||
            warp->lastComputeUsedGpuPath())
            ++failures;

        // Analytic cardioid/bulb membership is valid only for the conventional
        // escape radius >= 2. A smaller custom bailout must use CPU semantics.
        _putenv_s("MANDEL_BAILOUT", "1.1");
        std::vector<float> bailoutCpu(64, EMPTYPIXEL);
        std::vector<float> bailoutGpu(64, EMPTYPIXEL);
        Mandel bailoutCpuMandel(8, 8, MXIT, 1, bailoutCpu.data());
        Mandel bailoutGpuMandel(8, 8, MXIT, 1, bailoutGpu.data());
        mpf_set_d(centerRe, -1.24);
        mpf_set_ui(centerIm, 0);
        mpf_set_ui(scale, 100);
        bailoutCpuMandel.Compute(centerRe, centerIm, scale, MXIT, 0);
        ComputeRequest bailoutRequest;
        bailoutRequest.mode = ComputeMode::Mandelbrot;
        bailoutRequest.cpuEngine = &bailoutGpuMandel;
        bailoutRequest.centerRe = centerRe;
        bailoutRequest.centerIm = centerIm;
        bailoutRequest.scale = scale;
        bailoutRequest.width = 8;
        bailoutRequest.height = 8;
        bailoutRequest.sub = 1;
        bailoutRequest.maxIterations = MXIT;
        bailoutRequest.iterations = bailoutGpu.data();
        warp->resetCancellation();
        if (!warp->compute(bailoutRequest) ||
            warp->lastComputeUsedGpuPath() ||
            !same(bailoutCpu, bailoutGpu))
            ++failures;
        _putenv_s("MANDEL_BAILOUT", "");

        // Batched AVX2 refinement must preserve the scalar renderer's float
        // bailout-square rounding for custom radii.
        _putenv_s("MANDEL_BAILOUT", "2.1");
        float batchStorage[4] = {};
        Mandel batchMandel(2, 2, 64, 1, batchStorage);
        const double batchRe[] = {
            0.2522435803501477, -0.75, 0.3, -1.2
        };
        const double batchIm[] = {0.0, 0.1, 0.5, 0.0};
        float scalarPoints[4], batchPoints[4];
        for (int i = 0; i < 4; ++i)
            scalarPoints[i] =
                batchMandel.ComputeShallowPoint(batchRe[i], batchIm[i], 64);
        batchMandel.ComputeShallowPoints(
            batchRe, batchIm, 4, batchPoints, 64);
        if (std::memcmp(scalarPoints, batchPoints,
                        sizeof(scalarPoints)) != 0)
            ++failures;
        _putenv_s("MANDEL_BAILOUT", "");

        // The split c0/dx values can each fit while a far endpoint does not.
        // Reject before dispatch rather than feeding an overflowing dsMul.
        std::vector<float> rangeOutput((size_t)W * H, EMPTYPIXEL);
        Mandel rangeMandel(W, H, MXIT, 1, rangeOutput.data());
        mpf_set_d(centerRe, 4.0e9);
        mpf_set_d(centerIm, 2.666666666666667e9);
        mpf_set_d(scale, 5.0e-10);
        ComputeRequest rangeRequest;
        rangeRequest.mode = ComputeMode::Mandelbrot;
        rangeRequest.cpuEngine = &rangeMandel;
        rangeRequest.centerRe = centerRe;
        rangeRequest.centerIm = centerIm;
        rangeRequest.scale = scale;
        rangeRequest.width = W;
        rangeRequest.height = H;
        rangeRequest.sub = 1;
        rangeRequest.maxIterations = MXIT;
        rangeRequest.iterations = rangeOutput.data();
        warp->resetCancellation();
        if (!warp->compute(rangeRequest) ||
            warp->lastComputeUsedGpuPath())
            ++failures;

        mpf_set_d(centerRe, -0.5);
        mpf_set_ui(centerIm, 0);
        mpf_set_ui(scale, 1);
    }

    // Running cancellation: a bounded identity expression would otherwise execute
    // ten million iterations for every pixel.
    constexpr int CW = 8, CH = 8, CMXIT = 10000000;
    std::vector<float> cancelOutput((size_t)CW * CH, EMPTYPIXEL);
    Mandel cancelMandel(CW, CH, CMXIT, 1, cancelOutput.data());
    formula::ExpressionProgram identity;
    identity.compile("z", &compileError);
    formula::ExpressionContext cancelFixed;
    ComputeRequest cancelRequest;
    cancelRequest.mode = ComputeMode::Expression;
    cancelRequest.cpuEngine = &cancelMandel;
    cancelRequest.centerRe = centerRe; cancelRequest.centerIm = centerIm;
    cancelRequest.scale = scale;
    cancelRequest.width = CW; cancelRequest.height = CH; cancelRequest.sub = 1;
    cancelRequest.maxIterations = CMXIT;
    cancelRequest.iterations = cancelOutput.data();
    cancelRequest.expression = &identity;
    cancelRequest.expressionFixed = &cancelFixed;
    cancelRequest.expressionPixel = FormulaParameter::C;
    backend->resetCancellation();
    auto start = Clock::now();
    std::future<bool> running = std::async(std::launch::async, [&] {
        return backend->compute(cancelRequest);
    });
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
    backend->cancel();
    bool cancelResult = running.get();
    double cancelSeconds = since(start);
    backend->resetCancellation();
    int empty = (int)std::count(cancelOutput.begin(), cancelOutput.end(), EMPTYPIXEL);
    if (cancelResult || cancelSeconds > 2.0 || empty == 0) ++failures;

    printf("=== compute backend interface\n");
    printf("  backend=%s detail=%s fallback-test=%s\n",
           backend->info().name.c_str(), backend->info().detail.c_str(),
           fallback->info().detail.c_str());
    printf("  Mandel/Julia/Expression byte parity; cancel %.3fs empty=%d/%d\n",
           cancelSeconds, empty, CW * CH);
    printf("  D3D11 WARP: empty=%d class=%d floor=%d max smooth diff=%.6g\n",
           gpuEmpty, gpuClassMismatch, gpuFloorMismatch, gpuMaxDifference);
    if (gpuWorst >= 0)
        printf("    worst=(%d,%d) cpu=%.9g gpu=%.9g\n",
               gpuWorst % W, gpuWorst / W, gpuWorstCpu, gpuWorstGpu);
    printf("    scale1e6 class=%d floor=%d max diff=%.6g; sub5 layout=%d\n",
           gpuStressClass, gpuStressFloor, gpuStressMaxDifference,
           gpuLayoutMismatch);
    printf("  => %s\n\n", failures == 0 ? "PASS" : "CHECK (backend failure)");

    mpf_clears(centerRe, centerIm, scale, fixedCRe, fixedCIm, (mpf_ptr)0);
    return failures == 0 ? 0 : 1;
}

static int runGpuBenchmarkCase(int width, int height) {
    std::unique_ptr<IComputeBackend> gpu = createComputeBackend("gpu");
    printf("=== D3D11 hardware GPU benchmark\n");
    if (!gpu || gpu->info().fallback || !gpu->info().hardwareAccelerated) {
        printf("  SKIP: %s\n\n",
               gpu ? gpu->info().detail.c_str() : "backend creation failed");
        return 0;
    }

    const int mxit = [] {
        const char* value = getenv("MANDEL_GPU_MXIT");
        return value ? std::clamp(atoi(value), 2, 100000) : 5000;
    }();
    mpf_set_default_prec(128);
    mpf_t centerRe, centerIm, scale;
    mpf_init_set_d(centerRe, -0.5);
    mpf_init_set_ui(centerIm, 0);
    mpf_init_set_ui(scale, 1);

    const size_t count = (size_t)width * height;
    std::vector<float> cpuOutput(count, EMPTYPIXEL);
    std::vector<float> gpuOutput(count, EMPTYPIXEL);
    Mandel cpuMandel(width, height, mxit, 1, cpuOutput.data());
    Mandel gpuMandel(width, height, mxit, 1, gpuOutput.data());
    ComputeRequest request;
    request.mode = ComputeMode::Mandelbrot;
    request.cpuEngine = &gpuMandel;
    request.centerRe = centerRe;
    request.centerIm = centerIm;
    request.scale = scale;
    request.width = width;
    request.height = height;
    request.sub = 1;
    request.maxIterations = mxit;
    request.iterations = gpuOutput.data();

    // Warm resource allocation, shader execution, and the graphics driver's
    // pipeline cache before measuring a steady-state interactive frame.
    gpu->resetCancellation();
    bool warmupOk = gpu->compute(request);
    bool warmupUsedGpu = gpu->lastComputeUsedGpuPath();
    std::fill(gpuOutput.begin(), gpuOutput.end(), EMPTYPIXEL);
    _putenv_s("MANDEL_GPU_PROFILE", "1");
    gpu->resetCancellation();
    auto start = Clock::now();
    bool gpuOk = gpu->compute(request);
    bool measuredUsedGpu = gpu->lastComputeUsedGpuPath();
    double gpuSeconds = since(start);
    _putenv_s("MANDEL_GPU_PROFILE", "");

    _putenv_s("MANDEL_COARSE", "0");
    start = Clock::now();
    cpuMandel.Compute(centerRe, centerIm, scale, mxit, 0);
    double cpuSeconds = since(start);

    int empty = 0, classMismatch = 0, floorMismatch = 0;
    double maxDifference = 0.0;
    for (size_t i = 0; i < count; ++i) {
        if (gpuOutput[i] == EMPTYPIXEL) {
            ++empty;
            continue;
        }
        bool cpuInterior = isInterior(cpuOutput[i]);
        bool gpuInterior = isInterior(gpuOutput[i]);
        if (cpuInterior != gpuInterior) {
            ++classMismatch;
        } else if (!cpuInterior) {
            if ((int)cpuOutput[i] != (int)gpuOutput[i]) ++floorMismatch;
            maxDifference = std::max(
                maxDifference,
                std::fabs((double)cpuOutput[i] - (double)gpuOutput[i]));
        }
    }
    double speedup = gpuSeconds > 0.0 ? cpuSeconds / gpuSeconds : 0.0;
    bool accurate = warmupOk && warmupUsedGpu &&
                    gpuOk && measuredUsedGpu && empty == 0 &&
                    classMismatch == 0 && floorMismatch == 0 &&
                    maxDifference <= 0.01;
    bool fastEnough = speedup >= 5.0;
    printf("  backend=%s  %s\n",
           gpu->info().name.c_str(), gpu->info().detail.c_str());
    printf("  frame=%dx%d mxit=%d CPU/GPU=%.3f/%.3f s speedup=%.2fx\n",
           width, height, mxit, cpuSeconds, gpuSeconds, speedup);
    printf("  empty=%d class=%d floor=%d max smooth diff=%.6g\n",
           empty, classMismatch, floorMismatch, maxDifference);
    printf("  => %s\n\n",
           accurate && fastEnough
               ? "PASS"
               : (!accurate ? "CHECK (accuracy failure)"
                            : "CHECK (below 5x acceptance threshold)"));

    mpf_clears(centerRe, centerIm, scale, (mpf_ptr)0);
    return accurate && fastEnough ? 0 : 1;
}

static std::string pow10(int n) {
    std::string s = "1";
    s.append(n, '0');
    return s;
}

int main(int argc, char** argv) {
    std::string which = (argc > 1) ? argv[1] : "all";
    const int W = (argc > 2) ? atoi(argv[2]) : 120;
    const int H = (argc > 3) ? atoi(argv[3]) : 80;

    // step=1: brute-force the whole image. step>1: sparse oracle (feasible when
    // the location needs a high mxit and deep precision).
    TestCase shallow{ "shallow (whole set)", "-0.5", "0", pow10(0), 5000, 1 };
    std::string deep51_scale = "3831277"; deep51_scale.append(45, '0');   // 3.831277e51
    TestCase deep{ "deep51 (3.8e51 stress, maxit 2M)", testcases::deep51_x, testcases::deep51_y, deep51_scale, 2000000, 32 };
    TestCase ticktock{ "ticktock (1e141, glitch stress)", testcases::ticktock_x, testcases::ticktock_y, pow10(141), 200000, 12 };
    TestCase flake{ "flake (1e157, glitch stress)", testcases::flake_x, testcases::flake_y, pow10(157), 200000, 12 };
    // Exact Misiurewicz parameter c=i. Its critical orbit enters a repelling
    // period-2 cycle; 1e1000 pixel deltas amplify and escape after ~2650
    // iterations, making this a stable, genuinely deep exterior stress case.
    TestCase exterior1000{
        "exterior (Misiurewicz c=i, 1e1000)", "0", "1", pow10(1000), 10000, 1
    };
    // A maxit-sensitive period-1329 component boundary near c=i. This exposed a
    // deep-reference bug where even and odd image sizes followed different
    // reference pixels and flipped almost the entire image classification.
    const std::string parityRe =
        "0.87658692204855777255720127240577294219664089579606052252153341247110316040359567346941066530198264333715350076914731265081285605596159889659482671637665433036047479449476861255305601870749958937987160494800193599269036070384980604718254475114301562016620167318503180246283533354808944647785156562902419144968711947754196044352006391776697598867374433836375101918562145980477898443489741920475118804535699922081074341333337897641607372542163474622122834005970824831688364959496394361968461917806403290834529241567373997027330589016e-500";
    std::string parityIm = "0.";
    parityIm.append(498, '9');
    parityIm +=
        "9737732145267474019035193255007042699493697319030874172320234538271496366171723116782567955351384344445642267031664918871592758856681040534525429832129077075898037974437238709557577566760153377050600768882679447728184707424726714911339143303870647306255931485936680095871239708214775864294277313659941154622160496034912392379929833856520170763139643443634966071512688131637522087272155766755998538950242005247687225131897953568014863667043990248406875087728885765117677234374808520874070003498657884677143809794193877711208425737725547735873602";
    TestCase parity1000{
        "reference parity (period-1329 boundary, 1e1000)",
        parityRe, parityIm, pow10(1000), 10000, 16
    };
    // A 1e876 needle point (floatexp path, high-half multiply range). Deep exterior
    // stress with a different reference orbit than the 1e1000 cases.
    TestCase deep876{ "deep876 (1e876 needle, floatexp)", testcases::deep1_x, testcases::deep1_y, pow10(876), 300000, 24 };

    // Sub-pixel-nucleus stress (regression for the exterior-reference glitch):
    // the deep1 period-14164 minibrot is ~1e-876 wide, so at 1e737 it is ~139
    // orders of magnitude sub-pixel. The whole frame is therefore genuinely
    // EXTERIOR (mpmath confirms the centre escapes at n=41351 at this precision).
    // A full-image oracle catches both failure modes this exposed: a periodic
    // reference that renders the sub-pixel nucleus as a false-interior "black
    // star", and an escaping reference (periodic off / EDE / SAC) that is never
    // re-referenced and returns garbage smooth values. maxMisPct = 0 => any
    // false-interior pixel fails; step = 1 => every pixel is checked.
    TestCase subpixel{ "subpixel (deep1 nucleus 1e-876 @ 1e737, all-exterior)",
                       testcases::deep1_x, testcases::deep1_y, pow10(737), 60000, 1 };
    subpixel.maxMisPct = 0.0; subpixel.maxDiff = 0.5; subpixel.w = 48; subpixel.h = 36;

    // Acceptance gate for the real, resolved period-14164 minibrot at its natural
    // depth. The saved full-image oracle was generated at 1100 decimal digits and
    // its 81-interior / 351-exterior classification was stable at 1000 and 1100
    // digits. A render-precision oracle falsely reports the whole frame interior,
    // so this case MUST use the saved high-precision golden rather than ComputeDirect.
    TestCase minibrot875{ "minibrot875 (resolved period-14164 body + exterior)",
                         testcases::deep1_x, testcases::deep1_y, pow10(875), 300000, 1 };
    minibrot875.maxMisPct = 0.0; minibrot875.maxDiff = 1.0;
    minibrot875.w = 24; minibrot875.h = 18;
    minibrot875.goldenPath = "tests/golden/minibrot875.f32";

    // Shift the 1e875 frame right by 0.55 screen widths. The period-14164 nucleus
    // sits just outside the left edge, so findNucleus returns 0 and the exact view
    // centre is an EXTERIOR full reference (escapes near iteration 161971), while
    // the minibrot body remains visible along the left edge (27/432 interior).
    // This guards the full-reference/deep-zero path independently of the interior
    // periodic-nucleus reference used by minibrot875.
    TestCase extref875{ "extref875 (exterior full reference, visible minibrot)",
                        testcases::deep1_x, testcases::deep1_y, pow10(875), 300000, 1 };
    extref875.maxMisPct = 0.0; extref875.maxDiff = 1.0;
    extref875.w = 24; extref875.h = 18;
    extref875.goldenPath = "tests/golden/minibrot875_extref.f32";
    extref875.centerShiftRe = "2.2e-875";

    // Moderate-zoom high-iteration stress: BLA is ineffective and ~39% of pixels
    // are interior, so SIMD perturbation/interior detection dominates. The current
    // double perturbation has known smooth-value drift on a small chaotic tail; lock
    // classification exactly and bound max/mean/p99 so performance work cannot worsen it.
    TestCase slowpoint{ "slowpoint (7.826796e8, interior-detector stress)",
                        "-1.1758621450236620370", "-0.2447677973532398022",
                        "782679600", 500000, 1 };
    slowpoint.maxMisPct = 0.0; slowpoint.maxDiff = 90000.0;
    slowpoint.maxMeanDiff = 300.0; slowpoint.maxP99Diff = 3000.0;
    slowpoint.w = 48; slowpoint.h = 32;
    slowpoint.goldenPath = "tests/golden/slowpoint_7e8.f32";

    std::string point31Scale = "7354177"; point31Scale.append(25, '0');
    TestCase point31{ "point31 (7.354177e31, SA+BLA+SIMD stress)",
                      "-0.749139567333446841955467474699747367338762518832278501811",
                      "0.040823298514634751035521346975478853963578400940553676068",
                      point31Scale, 500000, 1 };
    point31.maxMisPct = 0.0; point31.maxDiff = 12000.0;
    point31.maxMeanDiff = 200.0; point31.maxP99Diff = 4000.0;
    point31.w = 48; point31.h = 32;
    point31.goldenPath = "tests/golden/point_7e31.f32";

    // Optional reference-footprint sweep: override mxit for all cases.
    if (const char* e = getenv("MANDEL_MXIT")) {
        int m = atoi(e);
        shallow.mxit = deep.mxit = ticktock.mxit = flake.mxit = exterior1000.mxit = parity1000.mxit = m;
    }
    if (const char* e = getenv("MANDEL_ORACLE_STEP")) {
        int step = std::max(1, atoi(e));
        shallow.step = deep.step = ticktock.step = flake.step = exterior1000.step = parity1000.step = step;
        deep876.step = subpixel.step = minibrot875.step = extref875.step =
            slowpoint.step = point31.step = step;
    }
    // Optional scale override (decimal exponent) to find interior-containing frames.
    if (const char* e = getenv("MANDEL_SCALE_EXP")) {
        std::string sc = pow10(atoi(e));
        deep.scale = ticktock.scale = flake.scale = exterior1000.scale = parity1000.scale = sc;
    }

    int rc = 0;
    if (which == "all" || which == "shallow")  rc |= runCase(shallow, W, H);
    if (which == "all" || which == "deep")     rc |= runCase(deep, W, H);
    if (which == "all" || which == "ticktock") rc |= runCase(ticktock, W, H);
    if (which == "all" || which == "flake")    rc |= runCase(flake, W, H);
    if (which == "exterior1000")               rc |= runCase(exterior1000, W, H);
    if (which == "parity1000")                 rc |= runCase(parity1000, W, H);
    if (which == "deep876")                    rc |= runCase(deep876, W, H);
    if (which == "subpixel")                   rc |= runCase(subpixel, W, H);
    if (which == "minibrot875")                rc |= runCase(minibrot875, W, H);
    if (which == "extref875")                  rc |= runCase(extref875, W, H);
    if (which == "slowpoint")                  rc |= runCase(slowpoint, W, H);
    if (which == "point31")                    rc |= runCase(point31, W, H);
    if (which == "gui875")                     rc |= runAdaptiveGuiCase();
    if (which == "julia")                      rc |= runJuliaCase(false);
    if (which == "julia-ede")                  rc |= runJuliaCase(true);
    if (which == "julia-dendrite")             rc |= runJuliaCase(false, true);
    if (which == "julia-critical")             rc |= runJuliaCase(false, true, true);
    if (which == "expression")                 rc |= runExpressionCoreCase();
    if (which == "expression-oracle")          rc |= runExpressionOracleCase();
    if (which == "expression-suite")           rc |= runFormulaRegressionSuite();
    if (which == "expression-residual")        rc |= runExpressionResidualSuite();
    if (which == "multibrot")                  rc |= runMultibrotCase();
    if (which == "backend")                    rc |= runBackendCase();
    if (which == "gpu")
        rc |= runGpuBenchmarkCase(
            argc > 2 ? W : 1920, argc > 3 ? H : 1080);
    return rc;
}
