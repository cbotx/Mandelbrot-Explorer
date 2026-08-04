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
//          expression | expression-oracle | all
//   The 1e1000 cases are excluded from "all" because their 3400-bit GMP oracles
//   are intentionally much more expensive than the regular regression set.

#include <gmp.h>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <string>
#include <chrono>
#include <algorithm>
#include <vector>

#include "mandel_perturbation.h"
#include "formula_expression.h"
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
    }

    ExpressionProgram burningShip;
    if (compile(burningShip,
                "sqr(complex(abs(re(z)), abs(im(z)))) + c")) {
        Complex folded{ std::abs(context.z.real()), std::abs(context.z.imag()) };
        if (!close(burningShip.evaluate(context), folded * folded + context.c))
            ++failures;
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
    return rc;
}
