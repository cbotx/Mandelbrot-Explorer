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
//          expression | expression-coloring | expression-orbit |
//          expression-centered | expression-reference | expression-scaled |
//          expression-taylor | expression-deep-render |
//          expression-oracle |
//          oracle | custom-deep | expression-suite | suite |
//          expression-residual | formula-bench | multibrot | backend | gpu | all
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
#include <atomic>
#include <future>
#include <iterator>
#include <limits>
#include <memory>
#include <new>
#include <set>
#include <thread>
#include <tuple>
#include <vector>

#include "compute_backend.h"
#include "mandel_perturbation.h"
#include "formula_expression.h"
#include "formula_expression_centered.h"
#include "formula_expression_orbit.h"
#include "formula_reference_orbit.h"
#include "formula_scaled_residual.h"
#include "formula_taylor_jet.h"
#include "formula_deep_renderer.h"
#include "formula_expression_jit.h"
#include "formula_expression_oracle.h"
#include "custom_deep_zoom.h"
#include "orbit_coloring.h"
#include "orbit_overlay.h"
#include "orbit_thumbnail.h"
#include "mandel_navigator.h"
#include "test_cases.h"

using Clock = std::chrono::high_resolution_clock;

#ifdef MANDEL_VERIFY_NO_JIT
static constexpr bool VERIFY_JIT = false;
#else
static constexpr bool VERIFY_JIT = true;
#endif
static double since(Clock::time_point t) {
    return std::chrono::duration_cast<std::chrono::milliseconds>(Clock::now() - t).count() / 1000.0;
}
struct TestCase {
    const char* name;
    std::string cx;
    std::string cy;
    std::string scale; // decimal magnitude, e.g. "1" or "1e30" expanded to digits
    int mxit;
    int step; // brute-force oracle samples a (step x step) pixel grid
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
static inline bool isInterior(float v) {
    return v < -1.5f;
}

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
        alternate = std::string("../") + path; // verify.exe is normally run from build/
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

    printf("=== %s  (%dx%d, scale=1e%zu, mxit=%d, prec=%d bits%s)\n", tc.name, W, H, tc.scale.size() - 1, tc.mxit, precision, c_method ? ", EDE" : "");

    auto t0 = Clock::now();
    mandel.Compute(cre, cim, scale, tc.mxit, c_method);
    double t_pert = since(t0);

    t0 = Clock::now();
    const bool useGolden = tc.goldenPath && getenv("MANDEL_FORCE_ORACLE") == nullptr;
    if (useGolden) {
        if (!loadGolden(tc.goldenPath, itd, (size_t)W * H)) {
            delete[] itp;
            delete[] itd;
            mpf_clear(cre);
            mpf_clear(cim);
            mpf_clear(scale);
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
            if (d == EMPTYPIXEL) continue; // not sampled by the oracle
            float p = itp[i * W + j];
            ++n_sampled;
            if (p == EMPTYPIXEL) {
                if (first_class_i < 0) {
                    first_class_i = i;
                    first_class_j = j;
                }
                ++n_engine_empty;
                ++n_class_mismatch;
                continue;
            }
            bool p_int = isInterior(p);
            bool d_int = isInterior(d);
            if (p_int != d_int) {
                if (first_class_i < 0) {
                    first_class_i = i;
                    first_class_j = j;
                }
                ++n_class_mismatch;
                continue;
            }
            if (p_int) {
                ++n_int_both;
                continue;
            }
            ++n_ext_both;
            double diff = std::fabs((double)p - (double)d);
            sum_diff += diff;
            diffs.push_back(diff);
            if (diff > max_diff) {
                max_diff = diff;
                worst_i = i;
                worst_j = j;
            }
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
    long n_interior = 0;
    double ext_iters = 0;
    for (int q = 0; q < W * H; ++q) {
        if (itp[q] == EMPTYPIXEL) continue;
        if (isInterior(itp[q]))
            ++n_interior;
        else
            ext_iters += itp[q];
    }
    double int_iters = (double)n_interior * tc.mxit;
    double frac_int = (int_iters + ext_iters) > 0 ? int_iters / (int_iters + ext_iters) : 0.0;
    printf("  perturbation : %8.3f s   (full %dx%d image)\n", t_pert, W, H);
    printf("  interior px  : %ld / %d (%.1f%%)   est. iters interior/exterior = %.3g / %.3g  (interior = %.1f%% of all iters)\n", n_interior, W * H, 100.0 * n_interior / (W * H), int_iters, ext_iters, 100.0 * frac_int);
    if (useGolden)
        printf("  golden load  : %8.3f s   (%ld pixels from %s)\n", t_direct, n_sampled, tc.goldenPath);
    else
        printf("  direct (ref) : %8.3f s   (%ld pixels sampled, step=%d)\n", t_direct, n_sampled, tc.step);
    printf("  interior both: %ld   exterior both: %ld   class mismatch: %ld (%.3f%% of sampled)\n", n_int_both, n_ext_both, n_class_mismatch, n_sampled ? 100.0 * n_class_mismatch / n_sampled : 0.0);
    if (n_engine_empty) printf("  ERROR: engine left %ld sampled pixels EMPTY\n", n_engine_empty);
    if (first_class_i >= 0) {
        int q = first_class_i * W + first_class_j;
        printf("  first class mismatch @ %d,%d: pert=%.4f ref=%.4f\n", first_class_i, first_class_j, itp[q], itd[q]);
    }
    printf("  escape-time diff vs reference:  max %.6g  mean %.6g  p99 %.6g", max_diff, mean_diff, p99_diff);
    if (worst_i >= 0) printf("   (worst @ %d,%d: pert=%.4f ref=%.4f)", worst_i, worst_j, itp[worst_i * W + worst_j], itd[worst_i * W + worst_j]);
    printf("\n  checksum(pert)=0x%08x\n", checksum(itp, W * H));

    // Boundary pixels can legitimately disagree, but the bulk must match.
    // Fully-exterior stress frames use maxMisPct = 0: a sub-pixel nucleus casts
    // no true interior, so ANY false-interior pixel (e.g. the periodic-reference
    // "black star") or garbage smooth value (an escaping reference never
    // re-referenced) is a real glitch, not a boundary rounding.
    double mismatch_pct = n_sampled ? 100.0 * n_class_mismatch / n_sampled : 0.0;
    bool ok = n_engine_empty == 0 && (max_diff <= tc.maxDiff) && (mean_diff <= tc.maxMeanDiff) && (p99_diff <= tc.maxP99Diff) && (mismatch_pct <= tc.maxMisPct);
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
           "(%dx%d, sub=%d, %zu sample slots)\n",
           W, H, SUB, count);
    auto t0 = Clock::now();
    mandel.Compute(cre, cim, scale, MXIT, ColoringMethod::SUPER_SAMPLING);
    double elapsed = since(t0);
    if (!loadGolden("tests/golden/minibrot875_gui_ss5.f32", golden.data(), count)) {
        mpf_clear(cre);
        mpf_clear(cim);
        mpf_clear(scale);
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
                if (firstEmptyX < 0) {
                    firstEmptyX = x;
                    firstEmptyY = y;
                }
                continue; // adaptive SS intentionally leaves unflagged subpixels empty
            }
            ++computed;
            bool isBase = (x % SUB == SUB / 2) && (y % SUB == SUB / 2);
            if (isBase)
                ++base;
            else
                ++sub;
            bool ei = isInterior(engine[q]), gi = isInterior(golden[q]);
            if (ei != gi) {
                ++classMismatch;
                continue;
            }
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
    printf("  GUI engine   : %8.3f s   base=%ld/%d adaptive-subpixels=%ld\n", elapsed, base, W * H, sub);
    printf("  computed     : %ld/%zu   class mismatch: %ld\n", computed, count, classMismatch);
    if (firstEmptyX >= 0) printf("  first empty  : sample (%d,%d), base=%d\n", firstEmptyX, firstEmptyY, (firstEmptyX % SUB == SUB / 2) && (firstEmptyY % SUB == SUB / 2));
    printf("  escape-time diff vs golden: max %.6g  mean %.6g  p99 %.6g\n", maxDiff, meanDiff, p99Diff);
    printf("  => %s\n\n", ok ? "PASS" : "CHECK (GUI adaptive path mismatch)");

    mpf_clear(cre);
    mpf_clear(cim);
    mpf_clear(scale);
    return ok ? 0 : 1;
}

static int runJuliaCase(bool ede, bool dendrite = false, bool critical = false) {
    const int W = critical ? 32 : 96, H = critical ? 24 : 64;
    const int mxit = dendrite ? 1000 : 500;
    mpf_set_default_prec(256);
    mpf_t centerRe, centerIm, scale, fixedRe, fixedIm;
    mpf_init_set_ui(centerRe, 0);
    mpf_init_set_ui(centerIm, 0);
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
    bool ok = mismatch == 0 && (dendrite ? (maxDiff <= 100.0 && meanDiff <= 0.1 && p99Diff <= 0.5) : (maxDiff <= (ede ? 1e-3 : 0.01)));
    if (critical) ok = ok && exterior == W * H;
    printf("=== quadratic Julia z0-plane (c=%s%s%s)  (%dx%d, mxit=%d)\n", dendrite ? "-0.8+0.156i" : "0", critical ? ", critical zoom" : "", ede ? ", EDE" : "", W, H, mxit);
    printf("  engine/direct: %.3f / %.3f s\n", engineTime, oracleTime);
    printf("  class mismatch: %d   max diff: %.6g   mean diff: %.6g   p99: %.6g\n", mismatch, maxDiff, meanDiff, p99Diff);
    if (critical) printf("  independently expected exterior: %d/%d\n", exterior, W * H);
    if (firstMismatch >= 0) printf("  first mismatch: (%d,%d) engine=%.6g oracle=%.6g\n", firstMismatch % W, firstMismatch / W, engine[firstMismatch], oracle[firstMismatch]);
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
    using formula::ExpressionOrbitPlan;
    using formula::ExpressionProgram;

    int failures = 0;
    auto close = [](Complex a, Complex b, double eps = 1e-12) { return std::abs(a - b) <= eps * std::max(1.0, std::abs(b)); };
    auto sameDoubleBits = [](double a, double b) {
        uint64_t aa, bb;
        std::memcpy(&aa, &a, sizeof(a));
        std::memcpy(&bb, &b, sizeof(b));
        return aa == bb;
    };
    auto sameComplexBits = [&](Complex a, Complex b) { return sameDoubleBits(a.real(), b.real()) && sameDoubleBits(a.imag(), b.imag()); };
    auto sameComplexResult = [&](Complex a, Complex b) {
        const auto sameComponent = [&](double left, double right) { return sameDoubleBits(left, right) || (std::isnan(left) && std::isnan(right)); };
        return sameComponent(a.real(), b.real()) && sameComponent(a.imag(), b.imag());
    };
    auto compile = [&](ExpressionProgram& program, const char* source) {
        ExpressionError error;
        if (!program.compile(source, &error)) {
            printf("  compile failed @ %zu: %s  [%s]\n", error.position, error.message.c_str(), source);
            ++failures;
            return false;
        }
        return true;
    };

    ExpressionContext context;
    context.z = {1.0, 2.0};
    context.c = {0.5, -0.25};
    context.z0 = {-0.75, 0.1};
    context.parameters[0] = {0.2, -0.3};
    context.iteration = 7;

    ExpressionProgram quadratic;
    if (compile(quadratic, "z*z + c") && !close(quadratic.evaluate(context), Complex{-2.5, 3.75})) ++failures;
    if (quadratic.fastPath() != ExpressionProgram::FastPath::IntegerPowerPlusC || quadratic.fastIntegerPower() != 2) ++failures;
    ExpressionProgram quadraticSquare, quadraticPower, cubicProduct;
    if (!compile(quadraticSquare, "sqr(z)+c") || !compile(quadraticPower, "z^2+c") || !compile(cubicProduct, "z*z*z+c") || quadraticSquare.fastPath() != ExpressionProgram::FastPath::IntegerPowerPlusC || quadraticSquare.fastIntegerPower() != 2 || quadraticPower.fastPath() != ExpressionProgram::FastPath::None || cubicProduct.fastPath() != ExpressionProgram::FastPath::IntegerPowerPlusC || cubicProduct.fastIntegerPower() != 3) ++failures;

    ExpressionProgram precedence;
    if (compile(precedence, "-2^2 + 2^3^2") && !close(precedence.evaluate(context), Complex{508.0, 0.0}, 1e-10)) ++failures;

    ExpressionProgram functions;
    const char* functionSource = "sin(z) + conj(p0) + complex(abs(re(c)), abs(im(c))) + n/10";
    if (compile(functions, functionSource)) {
        Complex expected = std::sin(context.z) + std::conj(context.parameters[0]) + Complex{std::abs(context.c.real()), std::abs(context.c.imag())} + Complex{0.7, 0.0};
        if (!close(functions.evaluate(context), expected)) ++failures;
        if (functions.fastPath() != ExpressionProgram::FastPath::None) ++failures;
        if (functions.avx2Compatible()) ++failures;
        if (!functions.batchCompatible()) ++failures;
    }

    ExpressionProgram burningShip;
    if (compile(burningShip, "sqr(complex(abs(re(z)), abs(im(z)))) + c")) {
        Complex folded{std::abs(context.z.real()), std::abs(context.z.imag())};
        if (!close(burningShip.evaluate(context), folded * folded + context.c)) ++failures;
        if (burningShip.derivativeCompatible()) ++failures;
        ExpressionDerivativeSeed seed;
        Complex value, derivative;
        if (burningShip.evaluateWithDerivative(context, seed, value, derivative)) ++failures;
    }

    ExpressionDerivativeSeed quadraticSeed;
    quadraticSeed.z = {0.3, -0.2};
    quadraticSeed.c = {-0.1, 0.4};
    Complex dualValue, dualDerivative;
    if (!quadratic.derivativeCompatible() || !quadratic.evaluateWithDerivative(context, quadraticSeed, dualValue, dualDerivative) || dualValue != context.z * context.z + context.c || dualDerivative != quadraticSeed.z * context.z + context.z * quadraticSeed.z + quadraticSeed.c) ++failures;

    ExpressionProgram analytic;
    if (!compile(analytic, "sin(z*z+c)+exp(p0*z0)+log(z+2)") || !analytic.derivativeCompatible()) {
        ++failures;
    } else {
        double maxDerivativeError = 0.0;
        for (int sample = 0; sample < 1000; ++sample) {
            ExpressionContext base = context;
            base.z = {0.2 + sample * 1e-5, -0.15 + sample * 2e-6};
            base.c = {0.4 - sample * 3e-6, 0.12};
            base.z0 = {-0.3, 0.25 + sample * 1e-6};
            base.parameters[0] = {0.6, -0.2};
            ExpressionDerivativeSeed direction;
            direction.z = {0.3, -0.1};
            direction.c = {-0.2, 0.25};
            direction.z0 = {0.1, 0.2};
            direction.parameters[0] = {-0.15, 0.05};
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
            double relative = std::abs(derivative - finiteDifference) / std::max(1.0, std::abs(derivative));
            maxDerivativeError = std::max(maxDerivativeError, relative);
        }
        if (!(maxDerivativeError < 2e-9)) {
            printf("  derivative finite-difference max error=%.6g\n", maxDerivativeError);
            ++failures;
        }
    }
    ExpressionProgram derivativeDomain;
    if (!compile(derivativeDomain, "log(z)") || !derivativeDomain.derivativeCompatible()) {
        ++failures;
    } else {
        ExpressionContext zero;
        ExpressionDerivativeSeed seed;
        seed.z = 1.0;
        Complex value, derivative;
        if (derivativeDomain.evaluateWithDerivative(zero, seed, value, derivative)) ++failures;
    }
    auto derivativeCase = [&](const char* source, Complex z, Complex seedValue, Complex expectedValue, Complex expectedDerivative, bool shouldSucceed) {
        ExpressionProgram program;
        if (!compile(program, source)) return;
        ExpressionContext dc;
        dc.z = z;
        ExpressionDerivativeSeed seed;
        seed.z = seedValue;
        Complex value, derivative;
        bool succeeded = program.evaluateWithDerivative(dc, seed, value, derivative);
        if (succeeded != shouldSucceed) {
            ++failures;
            return;
        }
        if (succeeded && (!close(value, expectedValue, 1e-12) || !close(derivative, expectedDerivative, 1e-12))) ++failures;
    };
    derivativeCase("z^2+c", {}, 1.0, {}, {}, true);
    {
        ExpressionProgram powerAtZero;
        if (!compile(powerAtZero, "z^2+c")) {
            ++failures;
        } else {
            ExpressionContext dc;
            dc.z = {};
            dc.c = {0.3, -0.2};
            ExpressionDerivativeSeed seed;
            seed.z = 1.0;
            seed.c = 1.0;
            Complex value, derivative;
            if (!powerAtZero.evaluateWithDerivative(dc, seed, value, derivative) || value != dc.c || derivative != Complex{1.0, 0.0}) ++failures;
        }
    }
    derivativeCase("log(z)", {-1.0, 0.0}, 1.0, {}, {}, false);
    derivativeCase("sqrt(z)", {-1.0, 0.0}, 1.0, {}, {}, false);
    derivativeCase("z^0.5", {-1.0, 0.0}, 1.0, {}, {}, false);
    derivativeCase("z^(-1)", {-1.0, 0.0}, 1.0, {-1.0, 0.0}, {-1.0, 0.0}, true);
    derivativeCase("z/z", {1e-200, 0.0}, 1.0, {1.0, 0.0}, {}, true);
    derivativeCase("tan(z)", {0.0, 1000.0}, 1.0, {0.0, 1.0}, {}, true);
    derivativeCase("tanh(z)", {1000.0, 0.0}, 1.0, {1.0, 0.0}, {}, true);
    {
        Complex z{0.0, 20.0};
        Complex expected = 1.0 / (std::cos(z) * std::cos(z));
        derivativeCase("tan(z)", z, 1.0, std::tan(z), expected, true);
        z = {20.0, 0.0};
        expected = 1.0 / (std::cosh(z) * std::cosh(z));
        derivativeCase("tanh(z)", z, 1.0, std::tanh(z), expected, true);
        double q = std::exp(-600.0);
        double saturatedDerivative = 4.0 * q / ((1.0 + q) * (1.0 + q));
        derivativeCase("tan(z)", {0.0, 300.0}, 1.0, std::tan(Complex{0.0, 300.0}), {saturatedDerivative, 0.0}, true);
        derivativeCase("tanh(z)", {300.0, 0.0}, 1.0, std::tanh(Complex{300.0, 0.0}), {saturatedDerivative, 0.0}, true);
    }

    ExpressionProgram invalid;
    ExpressionError invalidError;
    if (invalid.compile("z + * c", &invalidError) || invalidError.message.empty()) ++failures;
    std::string tooLong(ExpressionProgram::MAX_SOURCE + 1, '1');
    if (invalid.compile(tooLong, &invalidError)) ++failures;
    std::string tooMany = "z";
    for (size_t i = 0; i < ExpressionProgram::MAX_INSTRUCTIONS; ++i) tooMany += "+z";
    if (invalid.compile(tooMany, &invalidError)) ++failures;
    std::string deepPower = "z";
    for (size_t i = 0; i <= ExpressionProgram::MAX_PARSE_DEPTH; ++i) deepPower += "^z";
    if (invalid.compile(deepPower, &invalidError)) ++failures;
    std::string unaryChain(ExpressionProgram::MAX_SOURCE - 1, '-');
    unaryChain += "z";
    if (!invalid.compile(unaryChain, &invalidError)) ++failures;

    ExpressionProgram invalidPolar;
    if (compile(invalidPolar, "polar(-1, 0)")) {
        Complex result = invalidPolar.evaluate(context);
        if (!std::isnan(result.real()) || !std::isnan(result.imag())) ++failures;
    }

    int specializationFailures = 0;
    int specializationParityMismatches = 0;
    int specializationFoldCases = 0;
    ExpressionContext specializationFixed;
    specializationFixed.z0 = {-0.45, 0.125};
    specializationFixed.c = {-0.8, 0.156};
    for (int p = 0; p < 8; ++p) specializationFixed.parameters[p] = {0.25 + 0.125 * p, -0.5 + 0.0625 * p};

    for (int p = 0; p < 8; ++p) {
        std::string source = "p" + std::to_string(p);
        ExpressionProgram parameter;
        ExpressionProgram cPlane;
        ExpressionProgram z0Plane;
        if (!compile(parameter, source.c_str()) || !parameter.specialize(specializationFixed, FormulaParameter::C, cPlane) || !parameter.specialize(specializationFixed, FormulaParameter::InitialZ, z0Plane)) {
            ++specializationFailures;
            continue;
        }
        ExpressionContext changed = specializationFixed;
        changed.parameters[p] = {99.0 + p, -77.0 - p};
        if (!sameComplexBits(cPlane.evaluate(changed), specializationFixed.parameters[p]) || !sameComplexBits(z0Plane.evaluate(changed), specializationFixed.parameters[p]) || cPlane.source() != source || z0Plane.source() != source) ++specializationFailures;
    }

    auto checkDynamicVariable = [&](const char* source, FormulaParameter pixel, Complex expected, ExpressionContext runtime) {
        ExpressionProgram original;
        ExpressionProgram specializedProgram;
        if (!compile(original, source) || !original.specialize(specializationFixed, pixel, specializedProgram)) {
            ++specializationFailures;
            return;
        }
        if (!sameComplexBits(specializedProgram.evaluate(runtime), expected) || specializedProgram.source() != source) ++specializationFailures;
    };
    ExpressionContext changed = specializationFixed;
    changed.c = {4.0, -3.0};
    changed.z0 = {2.0, 5.0};
    changed.z = {-6.0, 7.0};
    changed.iteration = -11;
    checkDynamicVariable("c", FormulaParameter::C, changed.c, changed);
    checkDynamicVariable("c", FormulaParameter::InitialZ, specializationFixed.c, changed);
    checkDynamicVariable("z0", FormulaParameter::C, specializationFixed.z0, changed);
    checkDynamicVariable("z0", FormulaParameter::InitialZ, changed.z0, changed);
    checkDynamicVariable("z", FormulaParameter::C, changed.z, changed);
    checkDynamicVariable("z", FormulaParameter::InitialZ, changed.z, changed);
    checkDynamicVariable("n", FormulaParameter::C, {(double)changed.iteration, 0.0}, changed);
    checkDynamicVariable("n", FormulaParameter::InitialZ, {(double)changed.iteration, 0.0}, changed);

    const char* foldSources[] = {"-p0", "p0+p1", "p0-p1", "p0*p1", "p0/p1", "pow(p0,p1)", "sqr(p0)", "sin(p0)", "cos(p0)", "tan(p0)", "sinh(p0)", "cosh(p0)", "tanh(p0)", "exp(p0)", "log(p0)", "log10(p0)", "sqrt(p0)", "abs(p0)", "norm(p0)", "arg(p0)", "conj(p0)", "real(p0)", "imag(p0)", "complex(real(p0),imag(p1))", "polar(abs(p0),arg(p1))", "polar(-abs(p0),real(p1))", "sin(0.5)+exp(1)"};
    for (const char* source : foldSources) {
        ExpressionProgram original;
        ExpressionProgram specializedProgram;
        if (!compile(original, source) || !original.specialize(specializationFixed, FormulaParameter::C, specializedProgram)) {
            ++specializationFailures;
            continue;
        }
        ExpressionContext runtime = specializationFixed;
        runtime.parameters.fill({91.0, -37.0});
        Complex expected = original.evaluate(specializationFixed);
        Complex actual = specializedProgram.evaluate(runtime);
        if (specializedProgram.instructionCount() != 1 || specializedProgram.stackDepth() != 1 || specializedProgram.source() != source || !sameComplexBits(actual, expected)) ++specializationFailures;
        ++specializationFoldCases;
    }

    ExpressionProgram analyzedOriginal;
    ExpressionProgram analyzedRuntime;
    if (!compile(analyzedOriginal, "z*(sin(p0)+exp(p1))+abs(p2)") || !analyzedOriginal.specialize(specializationFixed, FormulaParameter::C, analyzedRuntime) || analyzedRuntime.stackDepth() >= analyzedOriginal.stackDepth() || analyzedRuntime.instructionCount() >= analyzedOriginal.instructionCount() || analyzedOriginal.avx2Compatible() || analyzedOriginal.derivativeCompatible() || !analyzedRuntime.avx2Compatible() || !analyzedRuntime.batchCompatible() || !analyzedRuntime.derivativeCompatible()) ++specializationFailures;

    ExpressionProgram quadraticCPlane;
    ExpressionProgram quadraticZ0Plane;
    ExpressionProgram literalRecurrence;
    ExpressionProgram literalRuntime;
    if (!quadratic.specialize(specializationFixed, FormulaParameter::C, quadraticCPlane) || !quadratic.specialize(specializationFixed, FormulaParameter::InitialZ, quadraticZ0Plane) || quadraticCPlane.fastIntegerPower() != 2 || quadraticZ0Plane.fastIntegerPower() != 2 || !compile(literalRecurrence, "z*z+1") || !literalRecurrence.specialize(specializationFixed, FormulaParameter::InitialZ, literalRuntime) || literalRuntime.fastPath() != ExpressionProgram::FastPath::None) ++specializationFailures;

    ExpressionProgram transactionTarget;
    ExpressionProgram uncompiled;
    ExpressionError specializationError;
    if (!compile(transactionTarget, "z+c")) {
        ++specializationFailures;
    } else {
        const std::string oldSource = transactionTarget.source();
        const size_t oldInstructions = transactionTarget.instructionCount();
        Complex oldValue = transactionTarget.evaluate(context);
        if (uncompiled.specialize(specializationFixed, FormulaParameter::C, transactionTarget, &specializationError) || specializationError.message.empty() || transactionTarget.source() != oldSource || transactionTarget.instructionCount() != oldInstructions || !sameComplexBits(transactionTarget.evaluate(context), oldValue)) ++specializationFailures;
        specializationError = {};
        if (quadratic.specialize(specializationFixed, FormulaParameter::Power, transactionTarget, &specializationError) || specializationError.message.empty() || transactionTarget.source() != oldSource || transactionTarget.instructionCount() != oldInstructions) ++specializationFailures;
    }

    const char* paritySources[] = {"z*z+c+sin(p0)+exp(p1)", "z*z+c+complex(real(p0),imag(p1))*z", "sin(p0)+cos(p1)+tan(p2)+z", "sinh(p0)+cosh(p1)-tanh(p2)+z", "exp(p0)+log(p1)+log10(p2)+sqrt(p3)+z", "abs(p0)+norm(p1)+arg(p2)+conj(p3)+z", "complex(real(p0),imag(p1))+z", "polar(abs(p0),arg(p1))+z", "p0/p1+pow(p2,p3)+z", "sqr(p0-p1)+z", "c+z0+z+n+p0", "pow(p0,p1)+c-z0+z", "log(p0)+sqrt(p1)+z", "1/p0+z", "polar(p0,p1)+z", "-p0*z+c", "conj(complex(real(p0),imag(p1)))*z+c", "p0*p1+p2/p3+z+n"};
    const double infValue = std::numeric_limits<double>::infinity();
    const double nanValue = std::numeric_limits<double>::quiet_NaN();
    const Complex specializationEdges[] = {{-1.0, 0.0}, {-1.0, -0.0}, {0.0, 0.0}, {-0.0, -0.0}, {1.0, 0.0}, {0.0, 1.0}, {infValue, 0.0}, {-infValue, -0.0}, {nanValue, 0.0}, {0.0, nanValue}, {1e308, -1e308}, {1e-308, -1e-308}};
    uint64_t specializationRandom = 0x243f6a8885a308d3ull;
    auto specializationRandomComponent = [&]() {
        specializationRandom ^= specializationRandom >> 12;
        specializationRandom ^= specializationRandom << 25;
        specializationRandom ^= specializationRandom >> 27;
        uint64_t bits = specializationRandom * 2685821657736338717ull;
        int32_t centered = (int32_t)(bits >> 32);
        return 4.0 * (double)centered / 2147483648.0;
    };
    for (const char* source : paritySources) {
        ExpressionProgram original;
        if (!compile(original, source)) {
            ++specializationFailures;
            continue;
        }
        for (FormulaParameter pixel : {FormulaParameter::C, FormulaParameter::InitialZ}) {
            for (int fixedSample = 0; fixedSample < 20; ++fixedSample) {
                ExpressionContext fixedSampleContext;
                if (fixedSample < 12) {
                    fixedSampleContext.c = specializationEdges[(fixedSample + 1) % 12];
                    fixedSampleContext.z0 = specializationEdges[(fixedSample + 3) % 12];
                    for (int p = 0; p < 8; ++p) fixedSampleContext.parameters[p] = specializationEdges[(fixedSample + p) % 12];
                } else {
                    fixedSampleContext.c = {specializationRandomComponent(), specializationRandomComponent()};
                    fixedSampleContext.z0 = {specializationRandomComponent(), specializationRandomComponent()};
                    for (Complex& parameter : fixedSampleContext.parameters) { parameter = {specializationRandomComponent(), specializationRandomComponent()}; }
                }
                ExpressionProgram specializedProgram;
                if (!original.specialize(fixedSampleContext, pixel, specializedProgram)) {
                    ++specializationFailures;
                    continue;
                }
                for (int dynamicSample = 0; dynamicSample < 8; ++dynamicSample) {
                    ExpressionContext reference = fixedSampleContext;
                    Complex pixelValue = dynamicSample < 4 ? specializationEdges[(fixedSample + dynamicSample * 2) % 12] : Complex{specializationRandomComponent(), specializationRandomComponent()};
                    if (pixel == FormulaParameter::C)
                        reference.c = pixelValue;
                    else
                        reference.z0 = pixelValue;
                    reference.z = dynamicSample < 4 ? specializationEdges[(fixedSample + dynamicSample * 3 + 5) % 12] : Complex{specializationRandomComponent(), specializationRandomComponent()};
                    reference.iteration = dynamicSample == 0 ? -1 : dynamicSample * 17 + fixedSample;
                    ExpressionContext runtime = reference;
                    for (int p = 0; p < 8; ++p) runtime.parameters[p] = specializationEdges[(fixedSample + dynamicSample + p + 7) % 12];
                    if (pixel == FormulaParameter::C)
                        runtime.z0 = specializationEdges[(fixedSample + dynamicSample + 9) % 12];
                    else
                        runtime.c = specializationEdges[(fixedSample + dynamicSample + 10) % 12];
                    Complex expected = original.evaluate(reference);
                    Complex actual = specializedProgram.evaluate(runtime);
                    if (!sameComplexBits(expected, actual)) {
                        if (specializationParityMismatches++ == 0) printf("  first specialization parity mismatch [%s] plane=%s fixed=%d dynamic=%d\n", source, pixel == FormulaParameter::C ? "c" : "z0", fixedSample, dynamicSample);
                    }
                }
            }
        }
    }

    int orbitPlanFailures = 0;
    int orbitPlanScalarMismatches = 0;
    int orbitPlanBatchMismatches = 0;
    int orbitPlanJitMismatches = 0;
    {
        ExpressionProgram original;
        ExpressionProgram runtime;
        ExpressionOrbitPlan plan;
        const char* source = "sin(c)+z";
        if (!compile(original, source) || !original.specialize(specializationFixed, FormulaParameter::C, runtime) || !plan.build(runtime) || !plan.profitable() || plan.source() != source || runtime.source() != source || original.source() != source || plan.dependencyMask() != (ExpressionOrbitPlan::DependencyZ | ExpressionOrbitPlan::DependencyC) || plan.invariantCount() != 1 || plan.invariantDependencyMask(0) != ExpressionOrbitPlan::DependencyC || plan.invariantInstructionCount(0) != 2 || plan.invariantOperationCount() != 1 || plan.bodyInstructionCount() != 3 || plan.bodyOperationCount() != 1 || !plan.matches(runtime)) { ++orbitPlanFailures; }

        ExpressionProgram duplicateOriginal;
        ExpressionProgram duplicateRuntime;
        ExpressionOrbitPlan duplicatePlan;
        const char* duplicateSource = "sin(c)+sin(c)+z";
        if (!compile(duplicateOriginal, duplicateSource) || !duplicateOriginal.specialize(specializationFixed, FormulaParameter::C, duplicateRuntime) || !duplicatePlan.build(duplicateRuntime) || duplicatePlan.invariantCount() != 2 || duplicatePlan.invariantOperationCount() != 2 || duplicatePlan.bodyInstructionCount() != 3 || duplicatePlan.bodyOperationCount() != 1 || duplicatePlan.source() != duplicateSource) { ++orbitPlanFailures; }

        ExpressionProgram dynamicOriginal;
        ExpressionProgram dynamicRuntime;
        ExpressionOrbitPlan dynamicPlan;
        const char* dynamicSource = "sin(z)+n+c";
        if (!compile(dynamicOriginal, dynamicSource) || !dynamicOriginal.specialize(specializationFixed, FormulaParameter::C, dynamicRuntime) || !dynamicPlan.build(dynamicRuntime) || dynamicPlan.profitable() || dynamicPlan.invariantCount() != 0 || dynamicPlan.bodyInstructionCount() != dynamicRuntime.instructionCount() || dynamicPlan.dependencyMask() != (ExpressionOrbitPlan::DependencyZ | ExpressionOrbitPlan::DependencyIteration | ExpressionOrbitPlan::DependencyC)) { ++orbitPlanFailures; }

        auto nanWithPayload = [](uint64_t payload) {
            uint64_t bits = 0x7ff8000000000000ull | (payload & 0x0007ffffffffffffull);
            double value;
            std::memcpy(&value, &bits, sizeof(value));
            return value;
        };
        ExpressionContext keyFixed = specializationFixed;
        keyFixed.parameters[0] = {nanWithPayload(0x1234), 0.0};
        keyFixed.parameters[1] = {nanWithPayload(0x5678), 0.0};
        auto checkKeyPlan = [&](const char* keySource, size_t expectedInvariants) {
            ExpressionProgram keyOriginal;
            ExpressionProgram keyRuntime;
            ExpressionOrbitPlan keyPlan;
            if (!compile(keyOriginal, keySource) || !keyOriginal.specialize(keyFixed, FormulaParameter::C, keyRuntime) || !keyPlan.build(keyRuntime) || keyPlan.invariantCount() != expectedInvariants) ++orbitPlanFailures;
        };
        checkKeyPlan("(c+p0)*z+(c+p0)", 1);
        checkKeyPlan("(c+p0)*z+(c+p1)", 2);
        checkKeyPlan("(c+p0)*z+(p0+c)", 2);
        keyFixed.parameters[0] = {0.0, -0.0};
        keyFixed.parameters[1] = {-0.0, 0.0};
        checkKeyPlan("(c+p0)*z+(c+p1)", 2);

        ExpressionOrbitPlan transactionPlan;
        ExpressionProgram transactionRuntime;
        ExpressionProgram uncompiledProgram;
        ExpressionError planError;
        if (!original.specialize(specializationFixed, FormulaParameter::C, transactionRuntime) || !transactionPlan.build(transactionRuntime)) {
            ++orbitPlanFailures;
        } else {
            const std::string oldSource = transactionPlan.source();
            const size_t oldInvariants = transactionPlan.invariantCount();
            if (transactionPlan.build(uncompiledProgram, &planError) || planError.message.empty() || transactionPlan.source() != oldSource || transactionPlan.invariantCount() != oldInvariants || !transactionPlan.matches(transactionRuntime)) { ++orbitPlanFailures; }
        }
    }

    const double payloadNaN = [] {
        uint64_t bits = 0x7ff8000000001234ull;
        double value;
        std::memcpy(&value, &bits, sizeof(value));
        return value;
    }();
    const Complex orbitEdges[] = {{0.0, 0.0}, {-0.0, 0.0}, {0.0, -0.0}, {-1.0, -0.0}, {1.0, 0.0}, {-2.0, 0.0}, {0.0, 1.0}, {1e-300, -1e-300}, {1e300, -1e300}, {std::numeric_limits<double>::infinity(), -0.0}, {-std::numeric_limits<double>::infinity(), 0.0}, {payloadNaN, -0.0}};
    struct OrbitParityCase {
        const char* source;
        FormulaParameter pixel;
        bool expectAvx;
        bool expectJit;
    };
    const OrbitParityCase orbitParityCases[] = {{"sqr(c)+z", FormulaParameter::C, true, true}, {"sin(c)+z", FormulaParameter::C, false, true}, {"sin(c)+sin(c)+z", FormulaParameter::C, false, true}, {"log(c)+sqrt(c)+pow(c,complex(0.5,0))+z", FormulaParameter::C, false, true}, {"sin(c)+sin(z)", FormulaParameter::C, false, false}, {"sin(z0)+z", FormulaParameter::InitialZ, false, true}, {"log(z0)+sqrt(z0)+z", FormulaParameter::InitialZ, false, true}};
    for (const OrbitParityCase& test : orbitParityCases) {
        ExpressionProgram original;
        ExpressionProgram runtime;
        ExpressionOrbitPlan plan;
        if (!compile(original, test.source) || !original.specialize(specializationFixed, test.pixel, runtime) || !plan.build(runtime) || !plan.profitable()) {
            ++orbitPlanFailures;
            continue;
        }
        if (runtime.avx2Compatible() != test.expectAvx) ++orbitPlanFailures;
        for (size_t sample = 0; sample < sizeof(orbitEdges) / sizeof(orbitEdges[0]); ++sample) {
            ExpressionContext lane = specializationFixed;
            if (test.pixel == FormulaParameter::C)
                lane.c = orbitEdges[sample];
            else
                lane.z0 = orbitEdges[sample];
            lane.z = orbitEdges[(sample * 5 + 3) % (sizeof(orbitEdges) / sizeof(orbitEdges[0]))];
            lane.iteration = (int)sample * 13 - 7;
            ExpressionOrbitPlan::Prepared prepared;
            if (!plan.prepare(lane, prepared)) {
                ++orbitPlanFailures;
                continue;
            }
            Complex expected = runtime.evaluate(lane);
            Complex actual = plan.evaluate(lane, prepared);
            if (!sameComplexBits(expected, actual)) ++orbitPlanScalarMismatches;
        }

        for (size_t first = 0; first < sizeof(orbitEdges) / sizeof(orbitEdges[0]); first += 4) {
            ExpressionContext lanes[4];
            ExpressionOrbitPlan::Prepared prepared[4];
            for (int lane = 0; lane < 4; ++lane) {
                size_t index = first + (size_t)lane;
                lanes[lane] = specializationFixed;
                if (test.pixel == FormulaParameter::C)
                    lanes[lane].c = orbitEdges[index];
                else
                    lanes[lane].z0 = orbitEdges[index];
                lanes[lane].z = orbitEdges[(index * 7 + 1) % (sizeof(orbitEdges) / sizeof(orbitEdges[0]))];
                lanes[lane].iteration = (int)index * 17 - 11;
            }
            Complex expected[4], actual[4];
            bool originalOkay = runtime.avx2Compatible() ? runtime.evaluate4(lanes, expected) : runtime.evaluate4Hybrid(lanes, expected);
            bool planOkay = plan.prepare4(lanes, 0x0f, prepared) && (plan.avx2Compatible() ? plan.evaluate4(lanes, prepared, actual) : plan.evaluate4Hybrid(lanes, prepared, actual));
            if (!originalOkay || !planOkay) {
                ++orbitPlanFailures;
                continue;
            }
            for (int lane = 0; lane < 4; ++lane) {
                if (!sameComplexResult(expected[lane], actual[lane])) {
                    if (orbitPlanBatchMismatches == 0) printf("  first orbit-plan batch mismatch source=%s sample=%zu lane=%d expected=(%.17g,%.17g) actual=(%.17g,%.17g)\n", test.source, first, lane, expected[lane].real(), expected[lane].imag(), actual[lane].real(), actual[lane].imag());
                    ++orbitPlanBatchMismatches;
                }
            }

            formula::ExpressionJit4 planJit;
            bool compiled = VERIFY_JIT && planJit.compile(plan);
            if (compiled != (VERIFY_JIT && test.expectJit)) {
                ++orbitPlanFailures;
            } else if (compiled) {
                formula::ExpressionJitInput4 input;
                formula::ExpressionJitInvariantInput4 invariantInput;
                formula::ExpressionJitOutput4 output;
                input.setContexts(lanes);
                for (int lane = 0; lane < 4; ++lane) invariantInput.setPreparedLane(lane, plan, prepared[lane]);
                planJit.evaluate(input, &invariantInput, output);
                for (int lane = 0; lane < 4; ++lane) {
                    Complex jitValue{output.re[lane], output.im[lane]};
                    if (!sameComplexResult(expected[lane], jitValue)) {
                        if (orbitPlanJitMismatches == 0) printf("  first orbit-plan JIT mismatch source=%s sample=%zu lane=%d expected=(%.17g,%.17g) actual=(%.17g,%.17g)\n", test.source, first, lane, expected[lane].real(), expected[lane].imag(), jitValue.real(), jitValue.imag());
                        ++orbitPlanJitMismatches;
                    }
                }
                Complex rejected[4];
                if (!planJit.supports(plan) || planJit.evaluate(lanes, rejected)) ++orbitPlanFailures;
            }
        }
    }

    auto checkOrbitPlanFrame = [&](const char* source, FormulaParameter pixel, bool expectJit) {
        ExpressionProgram original;
        ExpressionProgram runtime;
        ExpressionOrbitPlan plan;
        ExpressionContext fixed = specializationFixed;
        fixed.z0 = {0.1, -0.05};
        fixed.c = {-0.7, 0.2};
        if (!compile(original, source) || !original.specialize(fixed, pixel, runtime) || !plan.build(runtime) || !plan.profitable()) {
            ++orbitPlanFailures;
            return;
        }
        formula::ExpressionJit4 planJit;
        bool jitAvailable = VERIFY_JIT && planJit.compile(plan);
        if (jitAvailable != expectJit) ++orbitPlanFailures;
        constexpr int PW = 37, PH = 23, PMXIT = 90;
        std::vector<float> baseline((size_t)PW * PH);
        std::vector<float> scalarPlan((size_t)PW * PH);
        std::vector<float> batchPlan((size_t)PW * PH);
        std::vector<float> jitPlan((size_t)PW * PH);
        mpf_t re, im, scale;
        mpf_init_set_ui(re, 0);
        mpf_init_set_ui(im, 0);
        mpf_init_set_ui(scale, 1);
        auto render = [&](std::vector<float>& output, bool vector, const formula::ExpressionJit4* activeJit, const ExpressionOrbitPlan* activePlan) {
            std::fill(output.begin(), output.end(), EMPTYPIXEL);
            Mandel renderer(PW, PH, PMXIT, 1, output.data());
            _putenv_s("MANDEL_EXPR_VECTOR", vector ? "1" : "0");
            return renderer.ComputeExpression(re, im, scale, runtime, fixed, pixel, PMXIT, 8.0, formula::ExpressionColoring::Raw, activeJit, activePlan);
        };
        bool okay = render(baseline, false, nullptr, nullptr) && render(scalarPlan, false, nullptr, &plan) && render(batchPlan, true, nullptr, &plan);
        if (jitAvailable) okay = okay && render(jitPlan, true, &planJit, &plan);
        _putenv_s("MANDEL_EXPR_VECTOR", "");
        mpf_clears(re, im, scale, (mpf_ptr)0);
        if (!okay || std::memcmp(baseline.data(), scalarPlan.data(), baseline.size() * sizeof(float)) != 0 || std::memcmp(baseline.data(), batchPlan.data(), baseline.size() * sizeof(float)) != 0 || (jitAvailable && std::memcmp(baseline.data(), jitPlan.data(), baseline.size() * sizeof(float)) != 0)) ++orbitPlanFailures;
    };
    checkOrbitPlanFrame("sin(c)+z", FormulaParameter::C, VERIFY_JIT);
    checkOrbitPlanFrame("sin(z0)+z", FormulaParameter::InitialZ, VERIFY_JIT);
    checkOrbitPlanFrame("sin(c)+sin(z)", FormulaParameter::C, false);
    failures += orbitPlanFailures + orbitPlanScalarMismatches + orbitPlanBatchMismatches + orbitPlanJitMismatches;

    struct HybridOpcodeCase {
        const char* name;
        const char* source;
    };
    const HybridOpcodeCase hybridOpcodeCases[] = {{"divide", "z/p0"}, {"power", "pow(z,p0)"}, {"sin", "sin(z)"}, {"cos", "cos(z)"}, {"tan", "tan(z)"}, {"sinh", "sinh(z)"}, {"cosh", "cosh(z)"}, {"tanh", "tanh(z)"}, {"exp", "exp(z)"}, {"log", "log(z)"}, {"log10", "log10(z)"}, {"sqrt", "sqrt(z)"}, {"abs", "abs(z)"}, {"norm", "norm(z)"}, {"arg", "arg(z)"}, {"polar", "polar(real(z),real(p0))"}};
    uint64_t randomState = 0x9e3779b97f4a7c15ull;
    auto randomComponent = [&]() {
        randomState ^= randomState >> 12;
        randomState ^= randomState << 25;
        randomState ^= randomState >> 27;
        uint64_t bits = randomState * 2685821657736338717ull;
        int32_t centered = (int32_t)(bits >> 32);
        return 3.0 * (double)centered / 2147483648.0;
    };
    int hybridMismatch = 0;
    int hybridProgramsTested = 0;
    auto checkHybridLanes = [&](const char* name, const ExpressionProgram& program, ExpressionContext* lanes) {
        Complex batch[4];
        if (!program.evaluate4Hybrid(lanes, batch)) {
            ++hybridMismatch;
            return;
        }
        for (int lane = 0; lane < 4; ++lane) {
            Complex scalar = program.evaluate(lanes[lane]);
            if (!sameComplexBits(scalar, batch[lane])) {
                if (hybridMismatch++ == 0) printf("  first hybrid mismatch %s lane=%d scalar=(%.17g,%.17g) batch=(%.17g,%.17g)\n", name, lane, scalar.real(), scalar.imag(), batch[lane].real(), batch[lane].imag());
            }
        }
    };
    const double inf = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const Complex edgeValues[] = {{-1.0, 0.0}, {-1.0, -0.0}, {0.0, 0.0}, {-0.0, -0.0}, {1.0, 0.0}, {0.0, 1.0}, {inf, 0.0}, {-inf, -0.0}, {nan, 0.0}, {0.0, nan}, {1e308, -1e308}, {1e-308, -1e-308}};
    const Complex edgeParameters[] = {{0.5, 0.0}, {0.5, -0.0}, {0.0, 0.0}, {-0.0, 0.0}, {-1.0, 0.0}, {0.0, -1.0}, {2.0, 0.0}, {inf, 0.0}, {nan, 0.0}, {1.0, nan}, {1e308, 1e308}, {1e-308, 0.0}};
    for (const HybridOpcodeCase& test : hybridOpcodeCases) {
        ExpressionProgram program;
        if (!compile(program, test.source)) continue;
        ++hybridProgramsTested;
        if (!program.batchCompatible() || program.avx2Compatible()) ++hybridMismatch;
        for (int group = 0; group < 256; ++group) {
            ExpressionContext lanes[4];
            for (int lane = 0; lane < 4; ++lane) {
                lanes[lane] = context;
                lanes[lane].z = {randomComponent(), randomComponent()};
                lanes[lane].c = {randomComponent(), randomComponent()};
                lanes[lane].z0 = {randomComponent(), randomComponent()};
                lanes[lane].parameters[0] = {randomComponent(), randomComponent()};
                lanes[lane].iteration = group * 4 + lane;
            }
            checkHybridLanes(test.name, program, lanes);
        }
        for (size_t first = 0; first < sizeof(edgeValues) / sizeof(edgeValues[0]); first += 4) {
            ExpressionContext lanes[4];
            for (int lane = 0; lane < 4; ++lane) {
                size_t index = first + (size_t)lane;
                lanes[lane] = context;
                lanes[lane].z = edgeValues[index];
                lanes[lane].parameters[0] = edgeParameters[index];
            }
            checkHybridLanes(test.name, program, lanes);
        }
    }
    int vectorTranscendentalMismatch = 0;
    double vectorTranscendentalMaxRelativeError = 0.0;
    for (const char* source : {"sin(z)", "cos(z)", "tan(z)"}) {
        ExpressionProgram program;
        if (!compile(program, source)) {
            ++vectorTranscendentalMismatch;
            continue;
        }
        for (int group = 0; group < 256; ++group) {
            ExpressionContext lanes[4];
            for (int lane = 0; lane < 4; ++lane) {
                lanes[lane] = context;
                lanes[lane].z = {randomComponent(), randomComponent()};
            }
            Complex batch[4];
            if (!program.evaluate4Hybrid(lanes, batch, 1)) {
                ++vectorTranscendentalMismatch;
                continue;
            }
            for (int lane = 0; lane < 4; ++lane) {
                const Complex scalar = program.evaluate(lanes[lane]);
                const double scale = std::max(1.0, std::abs(scalar));
                const double relativeError = std::abs(batch[lane] - scalar) / scale;
                vectorTranscendentalMaxRelativeError = std::max(vectorTranscendentalMaxRelativeError, relativeError);
                if (!std::isfinite(relativeError) || relativeError > 5e-14) ++vectorTranscendentalMismatch;
            }
        }
        ExpressionContext signedZeroLanes[4] = {context, context, context, context};
        signedZeroLanes[0].z = {0.0, 0.0};
        signedZeroLanes[1].z = {-0.0, 0.0};
        signedZeroLanes[2].z = {0.0, -0.0};
        signedZeroLanes[3].z = {-0.0, -0.0};
        Complex signedZeroBatch[4];
        if (!program.evaluate4Hybrid(signedZeroLanes, signedZeroBatch, 1)) {
            ++vectorTranscendentalMismatch;
        } else {
            for (int lane = 0; lane < 4; ++lane)
                if (!sameComplexBits(program.evaluate(signedZeroLanes[lane]), signedZeroBatch[lane])) ++vectorTranscendentalMismatch;
        }
        if (std::strcmp(source, "tan(z)") == 0) {
            constexpr double Pi = 3.141592653589793238462643383279502884;
            ExpressionContext poleLanes[4] = {context, context, context, context};
            poleLanes[0].z = {0.5 * Pi, 8e-7};
            poleLanes[1].z = {-0.5 * Pi, 8e-7};
            poleLanes[2].z = {0.5 * Pi, -8e-7};
            poleLanes[3].z = {1.5 * Pi, 1e-6};
            Complex poleBatch[4];
            if (!program.evaluate4Hybrid(poleLanes, poleBatch, 1)) {
                ++vectorTranscendentalMismatch;
            } else {
                for (int lane = 0; lane < 4; ++lane)
                    if (!sameComplexBits(program.evaluate(poleLanes[lane]), poleBatch[lane])) ++vectorTranscendentalMismatch;
            }
        }
    }
    failures += vectorTranscendentalMismatch;
    {
        ExpressionProgram polarDomain;
        if (compile(polarDomain, "polar(z,p0)")) {
            const ExpressionContext polarCases[8] = {[&] {
                                                         ExpressionContext c;
                                                         c.z = {1.0, 1.0};
                                                         return c;
                                                     }(),
                                                     [&] {
                                                         ExpressionContext c;
                                                         c.z = {-1.0, 0.0};
                                                         return c;
                                                     }(),
                                                     [&] {
                                                         ExpressionContext c;
                                                         c.z = {inf, 0.0};
                                                         return c;
                                                     }(),
                                                     [&] {
                                                         ExpressionContext c;
                                                         c.z = {nan, 0.0};
                                                         return c;
                                                     }(),
                                                     [&] {
                                                         ExpressionContext c;
                                                         c.z = {1.0, 0.0};
                                                         c.parameters[0] = {0.0, 1.0};
                                                         return c;
                                                     }(),
                                                     [&] {
                                                         ExpressionContext c;
                                                         c.z = {1.0, 0.0};
                                                         c.parameters[0] = {inf, 0.0};
                                                         return c;
                                                     }(),
                                                     [&] {
                                                         ExpressionContext c;
                                                         c.z = {-0.0, 0.0};
                                                         c.parameters[0] = {-0.0, 0.0};
                                                         return c;
                                                     }(),
                                                     [&] {
                                                         ExpressionContext c;
                                                         c.z = {0.0, -0.0};
                                                         c.parameters[0] = {3.141592653589793, -0.0};
                                                         return c;
                                                     }()};
            for (int first = 0; first < 8; first += 4) {
                ExpressionContext lanes[4];
                for (int lane = 0; lane < 4; ++lane) lanes[lane] = polarCases[first + lane];
                checkHybridLanes("polar-domain", polarDomain, lanes);
            }
        }
    }
    {
        ExpressionProgram uncompiled;
        Complex outputs[4];
        ExpressionContext lanes[4]{};
        if (uncompiled.batchCompatible() || uncompiled.evaluate4Hybrid(lanes, outputs) || functions.evaluate4Hybrid(nullptr, outputs) || functions.evaluate4Hybrid(lanes, nullptr)) ++hybridMismatch;
    }
    failures += hybridMismatch;

    // The immutable bytecode must be safe to evaluate concurrently.
    std::atomic<int> parallelFailures{0};
#pragma omp parallel for schedule(static)
    for (int i = 0; i < 10000; ++i) {
        ExpressionContext local = context;
        local.z = {i * 1e-4, -i * 2e-4};
        Complex interpreted = quadratic.evaluate(local);
        std::array<Complex, ExpressionProgram::MAX_STACK> scratch;
        Complex reused = quadratic.evaluate(local, scratch.data(), quadratic.stackDepth());
        Complex expected = local.z * local.z + local.c;
        if (!close(interpreted, expected) || !close(reused, expected)) ++parallelFailures;
    }
    failures += parallelFailures.load();

    ExpressionProgram vectorProgram;
    formula::ExpressionJit4 vectorJit;
    double avxMs = 0.0, jitMs = 0.0, jitRawMs = 0.0;
    int jitMismatch = 0;
    if (!compile(vectorProgram, "z*z+c+p0*z+complex(re(z0),im(c))") || !vectorProgram.avx2Compatible()) {
        ++failures;
    } else if (VERIFY_JIT && !vectorJit.compile(vectorProgram)) {
        ++failures;
    } else if (VERIFY_JIT) {
        MEMORY_BASIC_INFORMATION memory{};
        SIZE_T queried = VirtualQuery(vectorJit.codeAddress(), &memory, sizeof(memory));
        DWORD executableWritable = PAGE_EXECUTE_READWRITE | PAGE_EXECUTE_WRITECOPY;
        if (!vectorJit.usesDualMapping() || queried != sizeof(memory) || !(memory.Protect & (PAGE_EXECUTE | PAGE_EXECUTE_READ | PAGE_EXECUTE_READWRITE | PAGE_EXECUTE_WRITECOPY)) || (memory.Protect & executableWritable)) ++failures;
        for (int group = 0; group < 2500; ++group) {
            ExpressionContext lanes[4];
            Complex vectorResults[4];
            Complex jitResults[4];
            for (int lane = 0; lane < 4; ++lane) {
                int index = group * 4 + lane;
                lanes[lane] = context;
                lanes[lane].z = {index * 1e-4, -index * 3e-5};
                lanes[lane].c = {0.2 + lane * 0.01, -0.3 + group * 1e-6};
                lanes[lane].z0 = {-0.5 + lane * 0.02, 0.1};
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
                if (scalar.real() != vectorResults[lane].real() || scalar.imag() != vectorResults[lane].imag()) ++failures;
                if (scalar.real() != jitResults[lane].real() || scalar.imag() != jitResults[lane].imag()) {
                    if (jitMismatch++ == 0) printf("  first JIT mismatch scalar=(%.17g,%.17g) jit=(%.17g,%.17g)\n", scalar.real(), scalar.imag(), jitResults[lane].real(), jitResults[lane].imag());
                }
            }
        }
        for (const char* source : {"-z", "conj(z)"}) {
            ExpressionProgram signProgram;
            if (!compile(signProgram, source) || !signProgram.avx2Compatible()) {
                ++failures;
                continue;
            }
            ExpressionContext lanes[4]{};
            lanes[0].z = {0.0, 0.0};
            lanes[1].z = {-0.0, 0.0};
            lanes[2].z = {0.0, -0.0};
            lanes[3].z = {-0.0, -0.0};
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
                if (!sameDoubleBits(scalar.real(), vectorResults[lane].real()) || !sameDoubleBits(scalar.imag(), vectorResults[lane].imag())) ++failures;
                if (!sameDoubleBits(scalar.real(), jitResults[lane].real()) || !sameDoubleBits(scalar.imag(), jitResults[lane].imag())) ++jitMismatch;
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
            ExpressionContext lanes[4] = {context, context, context, context};
            Complex outputs[4];
            const int iterations = 300000;
            volatile double sink = 0.0;
            auto begin = Clock::now();
            for (int i = 0; i < iterations; ++i) {
                lanes[0].iteration = i;
                vectorProgram.evaluate4(lanes, outputs);
                sink += outputs[0].real();
            }
            avxMs = std::chrono::duration<double, std::milli>(Clock::now() - begin).count();
            begin = Clock::now();
            for (int i = 0; i < iterations; ++i) {
                lanes[0].iteration = i;
                vectorJit.evaluate(lanes, outputs);
                sink += outputs[0].real();
            }
            jitMs = std::chrono::duration<double, std::milli>(Clock::now() - begin).count();
            formula::ExpressionJitInput4 input;
            formula::ExpressionJitOutput4 output;
            input.setContexts(lanes);
            begin = Clock::now();
            for (int i = 0; i < iterations; ++i) {
                for (int lane = 0; lane < 4; ++lane) input.vectors[formula::ExpressionJitInput4::N_RE][lane] = (double)i;
                vectorJit.evaluate(input, output);
                sink += output.re[0];
            }
            jitRawMs = std::chrono::duration<double, std::milli>(Clock::now() - begin).count();
            if (!std::isfinite(sink)) ++failures;
        }
    }

    // Canonical recurrence through the interpreter must classify the same pixels
    // as a hand-written quadratic Julia loop.
    int classificationMismatch = 0;
    for (int y = 0; y < 40; ++y) {
        for (int x = 0; x < 60; ++x) {
            Complex z0{-2.0 + 4.0 * x / 59.0, -1.3 + 2.6 * y / 39.0};
            Complex fixedC{-0.8, 0.156};
            ExpressionContext ec;
            ec.z = ec.z0 = z0;
            ec.c = fixedC;
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
    mpf_init_set_d(centerRe, -0.5);
    mpf_init_set_ui(centerIm, 0);
    mpf_init_set_ui(renderScale, 1);
    ExpressionContext fixed;
    fixed.z0 = {0.0, 0.0};
    if (!expressionRenderer.ComputeExpression(centerRe, centerIm, renderScale, quadratic, fixed, FormulaParameter::C, RMIT, 4.0)) {
        ++failures;
    } else {
        const double halfW = 2.0, halfH = halfW * RH / RW;
        const double dx = 2.0 * halfW / (RW - 1);
        const double dy = 2.0 * halfH / (RH - 1);
        for (int y = 0; y < RH; ++y) {
            for (int x = 0; x < RW; ++x) {
                Complex c{-0.5 - halfW + dx * x, -halfH + dy * y};
                Complex z{};
                float expected = -2.0f;
                for (int n = 0; n < RMIT; ++n) {
                    z = z * z + c;
                    if (std::hypot(z.real(), z.imag()) > 4.0) {
                        expected = (float)(n + 1);
                        break;
                    }
                }
                if (rendered[(size_t)y * RW + x] != expected) ++renderMismatch;
            }
        }
    }
    std::vector<float> specialized = rendered;
    ExpressionProgram genericQuadratic;
    if (!compile(genericQuadratic, "z*z+c+0") || !expressionRenderer.ComputeExpression(centerRe, centerIm, renderScale, genericQuadratic, fixed, FormulaParameter::C, RMIT, 4.0)) {
        ++failures;
    } else {
        for (size_t i = 0; i < rendered.size(); ++i)
            if (rendered[i] != specialized[i]) ++renderMismatch;
    }
    ExpressionProgram genericCubic;
    if (!expressionRenderer.ComputeExpression(centerRe, centerIm, renderScale, cubicProduct, fixed, FormulaParameter::C, RMIT, 4.0)) {
        ++failures;
    } else {
        specialized = rendered;
        if (!compile(genericCubic, "z*z*z+c+0") || !expressionRenderer.ComputeExpression(centerRe, centerIm, renderScale, genericCubic, fixed, FormulaParameter::C, RMIT, 4.0)) {
            ++failures;
        } else {
            for (size_t i = 0; i < rendered.size(); ++i)
                if (rendered[i] != specialized[i]) ++renderMismatch;
        }
    }

    int specializationFrameMismatch = 0;
    auto checkSpecializedFrame = [&](const char* source, FormulaParameter pixel, const ExpressionContext& frameFixed, double centerReal) {
        ExpressionProgram original;
        ExpressionProgram runtime;
        if (!compile(original, source) || !original.specialize(frameFixed, pixel, runtime)) {
            ++specializationFailures;
            return;
        }
        formula::ExpressionJit4 runtimeJit;
        const formula::ExpressionJit4* activeJit = runtime.fastPath() == ExpressionProgram::FastPath::None && runtimeJit.compile(runtime) ? &runtimeJit : nullptr;
        mpf_set_d(centerRe, centerReal);
        mpf_set_ui(centerIm, 0);
        mpf_set_ui(renderScale, 1);
        if (!expressionRenderer.ComputeExpression(centerRe, centerIm, renderScale, original, frameFixed, pixel, 240, 8.0)) {
            ++specializationFailures;
            return;
        }
        std::vector<float> originalFrame = rendered;
        if (!expressionRenderer.ComputeExpression(centerRe, centerIm, renderScale, runtime, frameFixed, pixel, 240, 8.0, formula::ExpressionColoring::Raw, activeJit)) {
            ++specializationFailures;
            return;
        }
        for (size_t i = 0; i < rendered.size(); ++i) {
            if (std::memcmp(&originalFrame[i], &rendered[i], sizeof(float)) != 0) ++specializationFrameMismatch;
        }
    };
    ExpressionContext invariantFixed;
    invariantFixed.z0 = {0.0, 0.0};
    invariantFixed.c = {-0.8, 0.156};
    invariantFixed.parameters[0] = {0.15, -0.2};
    invariantFixed.parameters[1] = {-0.35, 0.1};
    for (FormulaParameter pixel : {FormulaParameter::C, FormulaParameter::InitialZ}) {
        checkSpecializedFrame("z*z+c+sin(p0)+exp(p1)", pixel, invariantFixed, pixel == FormulaParameter::C ? -0.5 : 0.0);
        checkSpecializedFrame("z*z+c+complex(real(p0),imag(p1))*z", pixel, invariantFixed, pixel == FormulaParameter::C ? -0.5 : 0.0);
    }
    checkSpecializedFrame("z*z+c", FormulaParameter::InitialZ, invariantFixed, 0.0);
    renderMismatch += specializationFrameMismatch;
    failures += specializationFailures + specializationParityMismatches;

    fixed.c = {-0.8, 0.156};
    mpf_set_ui(centerRe, 0);
    if (!expressionRenderer.ComputeExpression(centerRe, centerIm, renderScale, quadratic, fixed, FormulaParameter::InitialZ, RMIT, 4.0)) {
        ++failures;
    } else {
        const double halfW = 2.0, halfH = halfW * RH / RW;
        const double dx = 2.0 * halfW / (RW - 1);
        const double dy = 2.0 * halfH / (RH - 1);
        for (int y = 0; y < RH; ++y) {
            for (int x = 0; x < RW; ++x) {
                Complex z{-halfW + dx * x, -halfH + dy * y};
                float expected = std::hypot(z.real(), z.imag()) > 4.0 ? 0.0f : -2.0f;
                if (expected < 0.0f) {
                    for (int n = 0; n < RMIT; ++n) {
                        z = z * z + fixed.c;
                        if (std::hypot(z.real(), z.imag()) > 4.0) {
                            expected = (float)(n + 1);
                            break;
                        }
                    }
                }
                if (rendered[(size_t)y * RW + x] != expected) ++renderMismatch;
            }
        }
    }
    // Degree-aware Smooth/EDE must agree with the generic automatic derivative.
    fixed = {};
    mpf_set_d(centerRe, -0.5);
    mpf_set_ui(centerIm, 0);
    mpf_set_ui(renderScale, 1);
    auto checkExpressionColoring = [&](formula::ExpressionColoring coloring) {
        int mismatches = 0;
        if (!expressionRenderer.ComputeExpression(centerRe, centerIm, renderScale, cubicProduct, fixed, FormulaParameter::C, RMIT, 4.0, coloring)) return RW * RH;
        const double halfW = 2.0, halfH = halfW * RH / RW;
        const double dx = 2.0 * halfW / (RW - 1);
        const double dy = 2.0 * halfH / (RH - 1);
        for (int y = 0; y < RH; ++y) {
            for (int x = 0; x < RW; ++x) {
                ExpressionContext dc;
                dc.z = dc.z0 = {};
                dc.c = {-0.5 - halfW + dx * x, -halfH + dy * y};
                Complex derivative{};
                float expected = -2.0f;
                for (int n = 0; n < RMIT; ++n) {
                    ExpressionDerivativeSeed seed;
                    seed.z = derivative;
                    seed.c = 1.0;
                    Complex next, nextDerivative;
                    if (!cubicProduct.evaluateWithDerivative(dc, seed, next, nextDerivative)) {
                        expected = 0.0f;
                        break;
                    }
                    dc.z = next;
                    derivative = nextDerivative;
                    double magnitude = std::hypot(next.real(), next.imag());
                    if (!std::isfinite(magnitude) || magnitude > 4.0) {
                        if (coloring == formula::ExpressionColoring::Smooth && std::isfinite(magnitude)) {
                            expected = (float)(n + 1 - std::log(std::log(magnitude)) / std::log(3.0));
                        } else if (coloring == formula::ExpressionColoring::Distance && std::isfinite(magnitude)) {
                            double denominator = std::abs(derivative) * dx;
                            expected = denominator > 0.0 && std::isfinite(denominator) ? (float)(magnitude * std::log(magnitude) / denominator) : 0.0f;
                        } else {
                            expected = (float)(n + 1);
                        }
                        break;
                    }
                }
                float actual = rendered[(size_t)y * RW + x];
                bool sameClass = isInterior(actual) == isInterior(expected);
                double tolerance = coloring == formula::ExpressionColoring::Distance ? 2e-3 * std::max(1.0, std::fabs((double)expected)) : 2e-4;
                if (!sameClass || (!isInterior(expected) && std::fabs((double)actual - expected) > tolerance)) ++mismatches;
            }
        }
        return mismatches;
    };
    renderMismatch += checkExpressionColoring(formula::ExpressionColoring::Smooth);
    renderMismatch += checkExpressionColoring(formula::ExpressionColoring::Distance);
    fixed.c = {-0.8, 0.156};
    mpf_set_ui(centerRe, 0);
    mpf_set_ui(centerIm, 0);
    mpf_set_d(renderScale, 0.4);
    const double z0HalfWidth = 5.0;
    const double z0HalfHeight = z0HalfWidth * RH / RW;
    const double firstMagnitude = std::hypot(z0HalfWidth, z0HalfHeight);
    if (!expressionRenderer.ComputeExpression(centerRe, centerIm, renderScale, cubicProduct, fixed, FormulaParameter::InitialZ, RMIT, 4.0, formula::ExpressionColoring::Smooth)) {
        ++failures;
    } else {
        float expected = (float)(-std::log(std::log(firstMagnitude)) / std::log(3.0));
        if (rendered[0] != expected) ++renderMismatch;
    }
    if (!expressionRenderer.ComputeExpression(centerRe, centerIm, renderScale, cubicProduct, fixed, FormulaParameter::InitialZ, RMIT, 4.0, formula::ExpressionColoring::Distance)) {
        ++failures;
    } else {
        double pixelDx = 2.0 * z0HalfWidth / (RW - 1);
        float expected = (float)(firstMagnitude * std::log(firstMagnitude) / pixelDx);
        if (rendered[0] != expected) ++renderMismatch;
    }
    // Coloring requests with bailout<1 must safely downgrade to raw counts.
    fixed = {};
    mpf_set_d(centerRe, -0.5);
    mpf_set_ui(centerIm, 0);
    mpf_set_ui(renderScale, 1);
    if (!expressionRenderer.ComputeExpression(centerRe, centerIm, renderScale, cubicProduct, fixed, FormulaParameter::C, 1, 0.5, formula::ExpressionColoring::Distance)) {
        ++failures;
    } else {
        for (float value : rendered)
            if (value >= 0.0f && value != 0.0f && value != 1.0f) ++renderMismatch;
    }

    // Extreme finite bailout values must use overflow-safe magnitude tests.
    ExpressionProgram identity;
    if (compile(identity, "z")) {
        fixed = {};
        fixed.z0 = {1e308, 1e308};
        mpf_set_ui(centerRe, 0);
        mpf_set_ui(centerIm, 0);
        if (!expressionRenderer.ComputeExpression(centerRe, centerIm, renderScale, identity, fixed, FormulaParameter::C, 1, 1e308) || rendered[0] != 0.0f) ++failures;
    }
    expressionRenderer.SetHalt(true);
    if (expressionRenderer.ComputeExpression(centerRe, centerIm, renderScale, quadratic, fixed, FormulaParameter::C, RMIT, 4.0)) ++failures;
    expressionRenderer.SetHalt(false);
    mpf_clears(centerRe, centerIm, renderScale, (mpf_ptr)0);
    failures += renderMismatch;

    printf("=== arbitrary expression core\n");
    printf("  quadratic instructions=%zu stack=%zu\n", quadratic.instructionCount(), quadratic.stackDepth());
    printf("  parallel failures=%d   classification mismatch=%d   render mismatch=%d\n", parallelFailures.load(), classificationMismatch, renderMismatch);
    printf("  hybrid opcode programs=%d lane mismatches=%d vector-transcendental mismatches=%d max-relative=%.3g\n", hybridProgramsTested, hybridMismatch, vectorTranscendentalMismatch, vectorTranscendentalMaxRelativeError);
    printf("  specialization folds=%d failures=%d scalar-bit mismatches=%d frame-bit mismatches=%d\n", specializationFoldCases, specializationFailures, specializationParityMismatches, specializationFrameMismatch);
    printf("  orbit plan failures=%d scalar/AVX-Hybrid/JIT mismatches=%d/%d/%d\n", orbitPlanFailures, orbitPlanScalarMismatches, orbitPlanBatchMismatches, orbitPlanJitMismatches);
    printf("  evaluator AVX2/JIT-wrapper/JIT-raw: %.3f / %.3f / %.3f ms  raw speedup %.2fx\n", avxMs, jitMs, jitRawMs, jitRawMs > 0.0 ? avxMs / jitRawMs : 0.0);
    printf("  JIT mismatches=%d\n", jitMismatch);
    if (vectorJit.valid()) {
        MEMORY_BASIC_INFORMATION memory{};
        VirtualQuery(vectorJit.codeAddress(), &memory, sizeof(memory));
        printf("  JIT W^X dual-mapping=%d execute-protect=0x%lx\n", vectorJit.usesDualMapping() ? 1 : 0, (unsigned long)memory.Protect);
    }
    printf("  => %s\n\n", failures == 0 ? "PASS" : "CHECK (expression core failure)");
    return failures == 0 ? 0 : 1;
}

static int runExpressionCenteredCase() {
    using formula::Complex;
    using formula::ExpressionCenteredEvaluator;
    using formula::ExpressionCenteredResult;
    using formula::ExpressionCenteredStatus;
    using formula::ExpressionContext;
    using formula::ExpressionDeltaContext;
    using formula::ExpressionError;
    using formula::ExpressionOracle;
    using formula::ExpressionOracleContext;
    using formula::ExpressionProgram;
    using formula::MpfrComplex;

    int failures = 0;
    int opcodeChecks = 0;
    int oracleChecks = 0;
    int tinyChecks = 0;
    int recurrenceChecks = 0;
    int fallbackChecks = 0;
    int signedZeroChecks = 0;
    constexpr double centeredPi = 3.141592653589793238462643383279502884;

    auto sameDoubleBits = [](double left, double right) {
        uint64_t leftBits = 0, rightBits = 0;
        std::memcpy(&leftBits, &left, sizeof(leftBits));
        std::memcpy(&rightBits, &right, sizeof(rightBits));
        return leftBits == rightBits;
    };
    auto sameComplexBits = [&](Complex left, Complex right) { return sameDoubleBits(left.real(), right.real()) && sameDoubleBits(left.imag(), right.imag()); };
    auto magnitude = [](Complex value) { return std::hypot(value.real(), value.imag()); };
    auto isZero = [](Complex value) { return value.real() == 0.0 && value.imag() == 0.0; };
    auto residualClose = [&](Complex actual, Complex expected, double relative = 5e-8) {
        double scale = magnitude(expected);
        double error = magnitude(actual - expected);
        return error <= relative * scale + 16.0 * std::numeric_limits<double>::denorm_min();
    };
    auto compile = [&](ExpressionProgram& program, const char* source) {
        ExpressionError error;
        if (!program.compile(source, &error)) {
            printf("  centered compile failed @ %zu: %s [%s]\n", error.position, error.message.c_str(), source);
            ++failures;
            return false;
        }
        return true;
    };
    auto perturbedContext = [](ExpressionContext base, const ExpressionDeltaContext& delta) {
        base.z += delta.z;
        base.c += delta.c;
        base.z0 += delta.z0;
        for (size_t i = 0; i < base.parameters.size(); ++i) base.parameters[i] += delta.parameters[i];
        return base;
    };
    auto setOracleContext = [](ExpressionOracleContext& output, const ExpressionContext& input) {
        output.z.set(input.z.real(), input.z.imag());
        output.c.set(input.c.real(), input.c.imag());
        output.z0.set(input.z0.real(), input.z0.imag());
        for (size_t i = 0; i < input.parameters.size(); ++i) output.parameters[i].set(input.parameters[i].real(), input.parameters[i].imag());
        output.iteration = input.iteration;
    };
    auto addOracleDelta = [](ExpressionOracleContext& output, const ExpressionDeltaContext& delta) {
        auto add = [](MpfrComplex& value, Complex residual) {
            mpfr_add_d(value.re, value.re, residual.real(), MPFR_RNDN);
            mpfr_add_d(value.im, value.im, residual.imag(), MPFR_RNDN);
        };
        add(output.z, delta.z);
        add(output.c, delta.c);
        add(output.z0, delta.z0);
        for (size_t i = 0; i < delta.parameters.size(); ++i) add(output.parameters[i], delta.parameters[i]);
    };
    auto oracleResidual = [&](const ExpressionProgram& program, const ExpressionContext& base, const ExpressionDeltaContext& delta, Complex& residual, mpfr_prec_t precision = 1536) {
        ExpressionOracleContext reference(precision);
        ExpressionOracleContext perturbed(precision);
        setOracleContext(reference, base);
        setOracleContext(perturbed, base);
        addOracleDelta(perturbed, delta);
        MpfrComplex referenceOutput(precision);
        MpfrComplex perturbedOutput(precision);
        std::string referenceError, perturbedError;
        if (!ExpressionOracle::evaluate(program, reference, referenceOutput, &referenceError) || !ExpressionOracle::evaluate(program, perturbed, perturbedOutput, &perturbedError)) {
            printf("  centered oracle failed: %s / %s [%s]\n", referenceError.c_str(), perturbedError.c_str(), program.source().c_str());
            return false;
        }
        MpfrComplex difference(precision);
        mpfr_sub(difference.re, perturbedOutput.re, referenceOutput.re, MPFR_RNDN);
        mpfr_sub(difference.im, perturbedOutput.im, referenceOutput.im, MPFR_RNDN);
        residual = difference.toDouble();
        return true;
    };
    auto expectStatus = [&](const char* source, ExpressionContext context, ExpressionDeltaContext delta, ExpressionCenteredStatus expected) {
        ExpressionProgram program;
        if (!compile(program, source)) return;
        ExpressionCenteredResult result = ExpressionCenteredEvaluator::evaluate(program, context, delta);
        ++fallbackChecks;
        if (result.status != expected) {
            printf("  centered status [%s]: got %s expected %s\n", source, formula::expressionCenteredStatusName(result.status), formula::expressionCenteredStatusName(expected));
            ++failures;
        }
        if (!sameComplexBits(result.base, program.evaluate(context))) {
            printf("  centered fallback base mismatch [%s]\n", source);
            ++failures;
        }
    };
    auto expectComposedBranch = [&](const char* source, ExpressionContext context, ExpressionDeltaContext delta) {
        ++signedZeroChecks;
        expectStatus(source, context, delta, ExpressionCenteredStatus::BranchUncertain);
    };
    auto expectSuccessParity = [&](const char* source, ExpressionContext context, ExpressionDeltaContext delta) {
        ExpressionProgram program;
        if (!compile(program, source)) return;
        ExpressionCenteredResult result = ExpressionCenteredEvaluator::evaluate(program, context, delta);
        ++signedZeroChecks;
        if (!result.success()) {
            printf("  centered benign signed-zero status [%s]: %s\n", source, formula::expressionCenteredStatusName(result.status));
            ++failures;
        }
        if (!sameComplexBits(result.base, program.evaluate(context))) {
            printf("  centered benign signed-zero base mismatch [%s]\n", source);
            ++failures;
        }
    };

    ExpressionContext nominal;
    nominal.z = {0.7, 0.2};
    nominal.c = {1.2, 0.3};
    nominal.z0 = {0.9, -0.4};
    nominal.parameters[0] = {0.6, 0.1};
    nominal.parameters[1] = {1.3, 0.0};
    nominal.parameters[2] = {0.4, 0.0};
    for (int i = 3; i < 8; ++i) nominal.parameters[i] = {0.2 + 0.07 * i, -0.15 + 0.03 * i};
    nominal.iteration = 7;

    const char* opcodeSources[] = {"2.0", "z+c+z0+p0+p1+p2+p3+p4+p5+p6+p7+n", "-z", "z+c", "z-c", "z*c", "z/c", "pow(z,p0)", "sqr(z)", "sin(z)", "cos(z)", "tan(z)", "sinh(z)", "cosh(z)", "tanh(z)", "exp(z)", "log(z)", "log10(z)", "sqrt(z)", "abs(z)", "norm(z)", "arg(z)", "conj(z)", "real(z)", "imag(z)", "complex(real(z),imag(c))", "polar(p1,p2)"};
    uint64_t randomState = 0x9e3779b97f4a7c15ULL;
    auto randomSigned = [&]() {
        randomState = randomState * 6364136223846793005ULL + 1;
        double unit = (double)((randomState >> 11) & ((1ULL << 53) - 1)) / (double)(1ULL << 53);
        return 2.0 * unit - 1.0;
    };
    for (int sample = 0; sample < 24; ++sample) {
        ExpressionContext context = nominal;
        context.z += Complex{0.08 * randomSigned(), 0.08 * randomSigned()};
        context.c += Complex{0.08 * randomSigned(), 0.08 * randomSigned()};
        context.z0 += Complex{0.08 * randomSigned(), 0.08 * randomSigned()};
        ExpressionDeltaContext delta;
        delta.z = {1e-7 * randomSigned(), 1e-7 * randomSigned()};
        delta.c = {1e-7 * randomSigned(), 1e-7 * randomSigned()};
        delta.z0 = {1e-7 * randomSigned(), 1e-7 * randomSigned()};
        for (int i = 0; i < 8; ++i) { delta.parameters[i] = {1e-7 * randomSigned(), 1e-7 * randomSigned()}; }
        delta.parameters[1] = {1e-7 * randomSigned(), 0.0};
        delta.parameters[2] = {1e-7 * randomSigned(), 0.0};

        for (const char* source : opcodeSources) {
            ExpressionProgram program;
            if (!compile(program, source)) continue;
            ExpressionCenteredResult result = ExpressionCenteredEvaluator::evaluate(program, context, delta);
            ++opcodeChecks;
            if (!result.success()) {
                printf("  centered random status %s [%s]\n", formula::expressionCenteredStatusName(result.status), source);
                ++failures;
                continue;
            }
            Complex directBase = program.evaluate(context);
            if (!sameComplexBits(result.base, directBase)) {
                printf("  centered base bit mismatch [%s]\n", source);
                ++failures;
            }
            Complex directPerturbed = program.evaluate(perturbedContext(context, delta));
            double reconstructionError = magnitude(result.base + result.delta - directPerturbed);
            if (reconstructionError > 2e-12 * std::max(1.0, magnitude(directPerturbed))) {
                printf("  centered reconstruction %.3g [%s]\n", reconstructionError, source);
                ++failures;
            }
            Complex expectedResidual;
            if (!oracleResidual(program, context, delta, expectedResidual, 512)) {
                ++failures;
            } else {
                ++oracleChecks;
                if (!residualClose(result.delta, expectedResidual, 2e-7)) {
                    printf("  centered oracle residual actual=(%.17g,%.17g)"
                           " expected=(%.17g,%.17g) [%s]\n",
                           result.delta.real(), result.delta.imag(), expectedResidual.real(), expectedResidual.imag(), source);
                    ++failures;
                }
            }
        }
    }

    struct TinyCase {
        const char* source;
        Complex base;
    };
    const TinyCase tinyCases[] = {{"exp(z)", {1.0, 0.2}}, {"sin(z)", {0.7, 0.2}}, {"log(z)", {1.2, 0.3}}, {"sqrt(z)", {1.2, 0.3}}, {"sin(z)+c", {0.4, -0.15}}};
    const int tinyExponents[] = {2, 8, 16, 40, 100, 200, 300};
    int directZeroDemonstrations = 0;
    for (const TinyCase& test : tinyCases) {
        ExpressionProgram program;
        if (!compile(program, test.source)) continue;
        for (int exponent : tinyExponents) {
            double amount = std::pow(10.0, -exponent);
            ExpressionContext context = nominal;
            context.z = test.base;
            context.c = {-0.2, 0.05};
            ExpressionDeltaContext delta;
            delta.z = {amount, 0.0};
            ExpressionCenteredResult result = ExpressionCenteredEvaluator::evaluate(program, context, delta);
            Complex expectedResidual;
            ++tinyChecks;
            if (!result.success() || !oracleResidual(program, context, delta, expectedResidual, 2048) || !residualClose(result.delta, expectedResidual, 5e-12)) {
                printf("  centered tiny mismatch 1e-%d [%s] status=%s\n", exponent, test.source, formula::expressionCenteredStatusName(result.status));
                ++failures;
            }
            Complex directDifference = program.evaluate(perturbedContext(context, delta)) - program.evaluate(context);
            if (exponent == 300 && isZero(directDifference) && !isZero(result.delta)) ++directZeroDemonstrations;
        }
    }
    if (directZeroDemonstrations != (int)std::size(tinyCases)) {
        printf("  centered tiny direct-zero demonstrations=%d/%zu\n", directZeroDemonstrations, std::size(tinyCases));
        ++failures;
    }

    struct RecurrenceCase {
        const char* source;
        Complex z;
        Complex c;
        Complex p0;
        ExpressionDeltaContext delta;
        int iterations;
    };
    ExpressionDeltaContext sineDelta;
    sineDelta.z = {1e-200, -2e-200};
    sineDelta.c = {3e-210, 1e-210};
    ExpressionDeltaContext expDelta;
    expDelta.z = {-2e-180, 1e-180};
    expDelta.c = {1e-200, -2e-200};
    ExpressionDeltaContext rationalDelta;
    rationalDelta.z = {2e-220, 1e-220};
    rationalDelta.c = {-1e-210, 1e-210};
    ExpressionDeltaContext powerDelta;
    powerDelta.z = {1e-160, -2e-160};
    powerDelta.c = {1e-190, 1e-190};
    powerDelta.parameters[0] = {2e-180, -1e-180};
    const RecurrenceCase recurrences[] = {{"sin(z)+c", {0.2, 0.1}, {-0.05, 0.02}, {}, sineDelta, 12}, {"exp(z)+c", {-1.0, 0.1}, {-1.0, -0.02}, {}, expDelta, 12}, {"z/(2+z)+c", {0.3, 0.1}, {0.02, -0.01}, {}, rationalDelta, 12}, {"pow(z,p0)+c", {0.8, 0.2}, {0.01, -0.02}, {1.3, 0.1}, powerDelta, 7}};
    for (const RecurrenceCase& test : recurrences) {
        ExpressionProgram program;
        if (!compile(program, test.source)) continue;
        ExpressionContext context = nominal;
        context.z = context.z0 = test.z;
        context.c = test.c;
        context.parameters[0] = test.p0;
        ExpressionDeltaContext delta = test.delta;
        for (int iteration = 0; iteration < test.iterations; ++iteration) {
            context.iteration = iteration;
            Complex expectedResidual;
            ExpressionCenteredResult result = ExpressionCenteredEvaluator::evaluate(program, context, delta);
            ++recurrenceChecks;
            if (!result.success() || !oracleResidual(program, context, delta, expectedResidual, 2048) || !residualClose(result.delta, expectedResidual, 2e-10)) {
                printf("  centered recurrence n=%d [%s] status=%s\n", iteration, test.source, formula::expressionCenteredStatusName(result.status));
                ++failures;
                break;
            }
            context.z = result.base;
            delta.z = result.delta;
        }
    }

    {
        ExpressionContext context;
        context.z = {-0.0, 0.0};
        context.c = {0.0, -0.0};
        context.z0 = {-0.0, -0.0};
        ExpressionDeltaContext delta;
        delta.z = {-0.0, 0.0};
        delta.c = {0.0, -0.0};
        delta.z0 = {-0.0, -0.0};
        ExpressionProgram linear, exponential;
        if (compile(linear, "complex(real(z),imag(c))+conj(z0)") && compile(exponential, "exp(complex(real(z),imag(c)))")) {
            ExpressionCenteredResult linearResult = ExpressionCenteredEvaluator::evaluate(linear, context, delta);
            ExpressionCenteredResult exponentialResult = ExpressionCenteredEvaluator::evaluate(exponential, context, delta);
            if (!linearResult.success() || !exponentialResult.success() || !sameComplexBits(linearResult.base, linear.evaluate(context)) || !sameComplexBits(exponentialResult.base, exponential.evaluate(context))) ++failures;
        }
    }

    {
        ExpressionContext context;
        ExpressionDeltaContext delta;
        context.z = {1e308, 0.0};
        delta.z = {1.0, 0.0};
        ExpressionProgram absolute;
        if (compile(absolute, "abs(z)")) {
            ExpressionCenteredResult result = ExpressionCenteredEvaluator::evaluate(absolute, context, delta);
            if (!result.success() || !(result.delta.real() > 0.5 && result.delta.real() < 1.5)) ++failures;
        }

        context.z = {0.25, 300.0};
        delta.z = {1e-12, 0.0};
        ExpressionProgram tangent;
        if (compile(tangent, "tan(z)")) {
            ExpressionCenteredResult result = ExpressionCenteredEvaluator::evaluate(tangent, context, delta);
            if (!result.success() || isZero(result.delta)) ++failures;
        }
        context.z = {300.0, 0.25};
        ExpressionProgram hyperbolicTangent;
        if (compile(hyperbolicTangent, "tanh(z)")) {
            ExpressionCenteredResult result = ExpressionCenteredEvaluator::evaluate(hyperbolicTangent, context, delta);
            if (!result.success() || isZero(result.delta)) ++failures;
        }
        context.z = {-700.0, 0.0};
        delta.z = {1e-10, 0.0};
        ExpressionProgram exponential;
        if (compile(exponential, "exp(z)")) {
            ExpressionCenteredResult result = ExpressionCenteredEvaluator::evaluate(exponential, context, delta);
            if (!result.success() || isZero(result.delta)) ++failures;
        }
        context.z = {1e-300, 1e-300};
        context.c = {2e-300, -2e-300};
        delta.z = {1e-310, -1e-310};
        delta.c = {2e-310, 1e-310};
        ExpressionProgram division;
        if (compile(division, "z/c")) {
            ExpressionCenteredResult result = ExpressionCenteredEvaluator::evaluate(division, context, delta);
            if (!result.success()) ++failures;
        }
    }

    {
        ExpressionContext context;
        ExpressionDeltaContext delta;
        context.z = {-1.0, 1.0};
        delta.z = {0.0, -2.0};
        expectStatus("log(z)", context, delta, ExpressionCenteredStatus::BranchUncertain);
        expectStatus("sqrt(z)", context, delta, ExpressionCenteredStatus::BranchUncertain);
        expectStatus("arg(z)", context, delta, ExpressionCenteredStatus::BranchUncertain);
        context.parameters[0] = {0.7, 0.2};
        expectStatus("pow(z,p0)", context, delta, ExpressionCenteredStatus::BranchUncertain);

        context.z = {1.0, 1.0};
        context.c = {-0.9999999999999999, 1.0};
        delta = {};
        delta.z = {-8.881784197001252e-16, 5.967448757360216e-16};
        delta.c = {1.4988010832439613e-15, -9.71445146547012e-17};
        expectStatus("log(z*c)", context, delta, ExpressionCenteredStatus::BranchUncertain);

        context.z = {-1.0, -0.0};
        delta = {};
        expectStatus("log(z)", context, delta, ExpressionCenteredStatus::BranchUncertain);
        expectStatus("sqrt(z)", context, delta, ExpressionCenteredStatus::BranchUncertain);
        expectStatus("arg(z)", context, delta, ExpressionCenteredStatus::BranchUncertain);
        expectStatus("pow(z,p0)", context, delta, ExpressionCenteredStatus::BranchUncertain);

        for (const char* source : {"log(conj(z))", "sqrt(conj(z))", "arg(conj(z))", "pow(conj(z),p0)", "log(complex(real(z),imag(z)))", "sqrt(complex(real(z),imag(z)))", "arg(complex(real(z),imag(z)))", "pow(complex(real(z),imag(z)),p0)"}) { expectComposedBranch(source, context, delta); }
        context.z = {1.0, -0.0};
        expectComposedBranch("log(-(complex(real(z),imag(z))+complex(0,-0)))", context, delta);

        context.z = {-1.0, -0.0};
        for (const char* source : {"exp(conj(z))", "log(complex(1,imag(z)))", "sqrt(complex(1,imag(z)))", "pow(complex(1,imag(z)),p0)"}) { expectSuccessParity(source, context, delta); }

        context.z = {1.0, 0.0};
        delta.z = {-2.0, 0.0};
        expectStatus("log(z)", context, delta, ExpressionCenteredStatus::Singular);
        context.c = {1.0, 0.0};
        delta = {};
        delta.c = {-2.0, 0.0};
        expectStatus("z/c", context, delta, ExpressionCenteredStatus::Singular);

        context.z = {centeredPi * 0.5 - 0.1, 0.0};
        delta = {};
        delta.z = {0.2, 0.0};
        expectStatus("tan(z)", context, delta, ExpressionCenteredStatus::Singular);
        context.z = {0.0, centeredPi * 0.5 - 0.1};
        delta.z = {0.0, 0.2};
        expectStatus("tanh(z)", context, delta, ExpressionCenteredStatus::Singular);

        context = nominal;
        context.parameters[1] = {0.1, 0.0};
        context.parameters[2] = {0.3, 0.0};
        delta = {};
        delta.parameters[1] = {-0.2, 0.0};
        expectStatus("polar(p1,p2)", context, delta, ExpressionCenteredStatus::Undefined);
        delta = {};
        delta.parameters[1] = {0.0, 1e-10};
        expectStatus("polar(p1,p2)", context, delta, ExpressionCenteredStatus::Undefined);

        context = {};
        delta = {};
        expectStatus("log(z)", context, delta, ExpressionCenteredStatus::Singular);
        expectStatus("sqrt(z)", context, delta, ExpressionCenteredStatus::Singular);
        expectStatus("pow(z,2)", context, delta, ExpressionCenteredStatus::Singular);
        context.z = {710.0, 0.0};
        expectStatus("exp(z)", context, delta, ExpressionCenteredStatus::NonFinite);
        context.z = {std::numeric_limits<double>::quiet_NaN(), 0.0};
        expectStatus("z+1", context, delta, ExpressionCenteredStatus::NonFinite);

        context.z = {-1.0, 0.0};
        ExpressionProgram exactCut;
        if (compile(exactCut, "log(z)")) {
            ExpressionCenteredResult result = ExpressionCenteredEvaluator::evaluate(exactCut, context, {});
            if (!result.success() || !isZero(result.delta) || !sameComplexBits(result.base, exactCut.evaluate(context))) ++failures;
        }
    }

    int batchChecks = 0;
    int batchMismatches = 0;
    const char* batchSources[] = {"z+c", "z-c", "z*c", "z/c", "pow(z,p0)", "sqr(z)", "sin(z)", "cos(z)", "tan(z)", "sinh(z)", "cosh(z)", "tanh(z)", "exp(z)", "log(z)", "log10(z)", "sqrt(z)", "abs(z)", "norm(z)", "arg(z)", "conj(z)", "real(z)", "imag(z)", "complex(real(z),imag(c))", "polar(real(p1),real(p2))", "sin(z*z+c)+exp(p0*z)/(2+z)"};
    std::array<ExpressionDeltaContext, 4> batchDeltas;
    batchDeltas[0].z = {1e-12, -2e-12};
    batchDeltas[0].c = {-1e-13, 2e-13};
    batchDeltas[0].parameters[0] = {1e-13, -1e-13};
    batchDeltas[1] = batchDeltas[0];
    batchDeltas[1].z *= -0.5;
    batchDeltas[1].c *= 2.0;
    batchDeltas[2] = batchDeltas[0];
    batchDeltas[2].z *= Complex{0.25, -0.75};
    batchDeltas[2].parameters[0] *= -1.0;
    batchDeltas[3] = {};
    for (const char* source : batchSources) {
        ExpressionProgram program;
        if (!compile(program, source)) continue;
        ExpressionCenteredResult batchResults[4];
        if (!ExpressionCenteredEvaluator::evaluate4(program, nominal, batchDeltas.data(), batchResults)) {
            ++batchMismatches;
            continue;
        }
        for (int lane = 0; lane < 4; ++lane) {
            const ExpressionCenteredResult scalar = ExpressionCenteredEvaluator::evaluate(program, nominal, batchDeltas[static_cast<size_t>(lane)]);
            ++batchChecks;
            if (scalar.status != batchResults[lane].status || !sameComplexBits(scalar.base, batchResults[lane].base) || !sameComplexBits(scalar.delta, batchResults[lane].delta)) ++batchMismatches;
        }
    }
    {
        ExpressionProgram canonical;
        ExpressionProgram runtime;
        formula::ExpressionOrbitPlan plan;
        ExpressionError error;
        if (!canonical.compile("sin(c)+z", &error) || !canonical.specialize(nominal, FormulaParameter::C, runtime, &error) || !plan.build(runtime, &error) || !plan.profitable()) {
            ++batchMismatches;
        } else {
            ExpressionCenteredResult batchResults[4];
            if (!ExpressionCenteredEvaluator::evaluate4(plan.bodyProgram(), nominal, batchDeltas.data(), batchResults)) {
                ++batchMismatches;
            } else {
                for (int lane = 0; lane < 4; ++lane) {
                    ++batchChecks;
                    const ExpressionCenteredResult scalar = ExpressionCenteredEvaluator::evaluate(plan.bodyProgram(), nominal, batchDeltas[static_cast<size_t>(lane)]);
                    if (scalar.status != ExpressionCenteredStatus::Unsupported || batchResults[lane].status != scalar.status) ++batchMismatches;
                }
            }
        }
    }
    failures += batchMismatches;

    ExpressionProgram benchmarkProgram;
    double scalarMs = 0.0, centeredMs = 0.0, scalar4Ms = 0.0, batch4Ms = 0.0;
    if (compile(benchmarkProgram, "sin(z*z+c)+exp(p0*z)/(2+z)")) {
        ExpressionContext context = nominal;
        ExpressionDeltaContext delta;
        delta.z = {1e-12, -2e-12};
        delta.c = {-1e-13, 2e-13};
        delta.parameters[0] = {1e-13, -1e-13};
        constexpr int repetitions = 100000;
        Complex scalarSink{}, centeredSink{};
        auto begin = Clock::now();
        for (int i = 0; i < repetitions; ++i) {
            context.iteration = i;
            scalarSink += benchmarkProgram.evaluate(context);
        }
        scalarMs = since(begin) * 1000.0;
        begin = Clock::now();
        for (int i = 0; i < repetitions; ++i) {
            context.iteration = i;
            ExpressionCenteredResult result = ExpressionCenteredEvaluator::evaluate(benchmarkProgram, context, delta);
            centeredSink += result.base + result.delta;
        }
        centeredMs = since(begin) * 1000.0;
        std::array<ExpressionDeltaContext, 4> deltas = batchDeltas;
        begin = Clock::now();
        for (int i = 0; i < repetitions; ++i) {
            context.iteration = i;
            for (int lane = 0; lane < 4; ++lane) {
                ExpressionCenteredResult result = ExpressionCenteredEvaluator::evaluate(benchmarkProgram, context, deltas[static_cast<size_t>(lane)]);
                centeredSink += result.base + result.delta;
            }
        }
        scalar4Ms = since(begin) * 1000.0;
        begin = Clock::now();
        for (int i = 0; i < repetitions; ++i) {
            context.iteration = i;
            ExpressionCenteredResult batch[4];
            if (!ExpressionCenteredEvaluator::evaluate4(benchmarkProgram, context, deltas.data(), batch)) {
                ++failures;
                break;
            }
            for (const ExpressionCenteredResult& result : batch) centeredSink += result.base + result.delta;
        }
        batch4Ms = since(begin) * 1000.0;
        if (!std::isfinite(scalarSink.real()) || !std::isfinite(centeredSink.real())) ++failures;
    }

    printf("=== expression centered residual bytecode v1\n");
    printf("  opcode/base/reconstruction=%d oracle=%d tiny=%d"
           " recurrence=%d fallback=%d signed-zero=%d\n",
           opcodeChecks, oracleChecks, tinyChecks, recurrenceChecks, fallbackChecks, signedZeroChecks);
    printf("  direct subtraction zero while centered nonzero=%d/%zu\n", directZeroDemonstrations, std::size(tinyCases));
    printf("  centered batch parity=%d mismatch=%d\n", batchChecks, batchMismatches);
    printf("  scalar/centered %.3f/%.3f ms overhead %.2fx"
           " (informational)\n",
           scalarMs, centeredMs, scalarMs > 0.0 ? centeredMs / scalarMs : 0.0);
    printf("  centered scalar4/batch4 %.3f/%.3f ms speedup %.2fx\n", scalar4Ms, batch4Ms, batch4Ms > 0.0 ? scalar4Ms / batch4Ms : 0.0);
    printf("  => %s\n\n", failures == 0 ? "PASS" : "CHECK (centered expression failure)");
    return failures == 0 ? 0 : 1;
}

static bool renderExpressionColorReference(const formula::ExpressionProgram& program, const formula::ExpressionContext& fixed, FormulaParameter pixel, double centerRe, double centerIm, double scale, int width, int height, int mxit, double bailout, formula::ExpressionColoring coloring, std::vector<float>& output) {
    if (!program.valid() || width < 2 || height < 2 || !(scale > 0.0) || !(bailout > 0.0)) return false;
    const double halfWidth = 2.0 / scale;
    const double halfHeight = halfWidth * height / width;
    const double startRe = centerRe - halfWidth;
    const double startIm = centerIm - halfHeight;
    const double dx = 2.0 * halfWidth / (width - 1);
    const double dy = 2.0 * halfHeight / (height - 1);
    const int power = program.fastPath() == formula::ExpressionProgram::FastPath::IntegerPowerPlusC ? program.fastIntegerPower() : 0;
    output.assign((size_t)width * height, EMPTYPIXEL);
    std::array<formula::Complex, formula::ExpressionProgram::MAX_STACK> stack;
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            formula::ExpressionContext context = fixed;
            formula::Complex coordinate{startRe + dx * x, startIm + dy * y};
            if (pixel == FormulaParameter::C)
                context.c = coordinate;
            else
                context.z0 = coordinate;
            context.z = context.z0;
            double initialRe = context.z.real();
            double initialIm = context.z.imag();
            if (!std::isfinite(initialRe) || !std::isfinite(initialIm) || std::hypot(initialRe, initialIm) > bailout) {
                output[(size_t)y * width + x] = 0.0f;
                continue;
            }

            orbitcolor::FormulaColorAccum accumulator;
            if (coloring == formula::ExpressionColoring::Feather || coloring == formula::ExpressionColoring::OrbitTrap) { accumulator.init(coloring == formula::ExpressionColoring::OrbitTrap ? orbitcolor::FormulaColorMode::OrbitTrap : orbitcolor::FormulaColorMode::Feather); }
            float result = coloring == formula::ExpressionColoring::OrbitTrap ? accumulator.interior() : -2.0f;
            for (int iteration = 0; iteration < mxit; ++iteration) {
                context.iteration = iteration;
                context.z = program.evaluate(context, stack.data(), program.stackDepth());
                double re = context.z.real();
                double im = context.z.imag();
                double magnitude = std::hypot(re, im);
                if (coloring == formula::ExpressionColoring::Feather || coloring == formula::ExpressionColoring::OrbitTrap) accumulator.push(re, im);
                if (!std::isfinite(re) || !std::isfinite(im) || magnitude > bailout) {
                    if (coloring == formula::ExpressionColoring::Feather || coloring == formula::ExpressionColoring::OrbitTrap) {
                        result = accumulator.escaped(iteration + 1, magnitude, power, bailout);
                    } else {
                        result = (float)(iteration + 1);
                    }
                    break;
                }
                if (iteration + 1 == mxit && coloring == formula::ExpressionColoring::OrbitTrap) result = accumulator.interior();
            }
            output[(size_t)y * width + x] = result;
        }
    }
    return true;
}

static int runExpressionColoringCase() {
    struct Case {
        const char* name;
        const char* source;
        FormulaParameter pixel;
        formula::ExpressionContext fixed;
        double centerRe = 0.0;
        double centerIm = 0.0;
        double scale = 1.0;
        double bailout = 4.0;
        int mxit = 80;
    };

    std::vector<Case> cases;
    Case arithmetic{"arithmetic-jit", "z*z+c+p0*z", FormulaParameter::C};
    arithmetic.fixed.parameters[0] = {0.1, -0.05};
    arithmetic.centerRe = -0.4;
    cases.push_back(arithmetic);
    Case arithmeticZ0 = arithmetic;
    arithmeticZ0.name = "arithmetic-jit-z0";
    arithmeticZ0.pixel = FormulaParameter::InitialZ;
    arithmeticZ0.fixed.c = {-0.35, 0.62};
    arithmeticZ0.centerRe = 0.0;
    cases.push_back(arithmeticZ0);
    Case transcendental{"transcendental-hybrid", "sin(z)+c", FormulaParameter::C};
    transcendental.fixed.z0 = {0.1, -0.05};
    transcendental.bailout = 8.0;
    cases.push_back(transcendental);
    Case invariant{"orbit-plan-jit", "0.5*z+0.1*sin(c)", FormulaParameter::C};
    invariant.mxit = 120;
    cases.push_back(invariant);
    Case powerC{"integer-power-c", "z*z*z+c", FormulaParameter::C};
    cases.push_back(powerC);
    Case powerZ0{"integer-power-z0", "z*z*z+c", FormulaParameter::InitialZ};
    powerZ0.fixed.c = {-0.35, 0.62};
    powerZ0.scale = 0.8;
    cases.push_back(powerZ0);
    Case interior{"interior", "0.5*z", FormulaParameter::C};
    interior.fixed.z0 = {};
    interior.scale = 4.0;
    interior.mxit = 20;
    cases.push_back(interior);
    Case initialEscape{"initial-escape", "z*z+c", FormulaParameter::InitialZ};
    initialEscape.fixed.c = {-0.8, 0.156};
    initialEscape.scale = 0.4;
    cases.push_back(initialEscape);
    Case domain{"domain-nonfinite", "log(z)+c", FormulaParameter::C};
    domain.fixed.z0 = {};
    domain.bailout = 8.0;
    cases.push_back(domain);

    constexpr int width = 7;
    constexpr int height = 5;
    int failures = 0;
    int comparisons = 0;
    int mismatches = 0;
    int sentinelCollisions = 0;
    printf("=== universal expression coloring\n");

    for (const Case& test : cases) {
        formula::ExpressionProgram source;
        formula::ExpressionProgram runtime;
        formula::ExpressionOrbitPlan plan;
        formula::ExpressionError error;
        if (!source.compile(test.source, &error) || !source.specialize(test.fixed, test.pixel, runtime, &error) || !plan.build(runtime, &error)) {
            ++failures;
            continue;
        }
        formula::ExpressionJit4 jit;
        bool jitAvailable = runtime.fastPath() == formula::ExpressionProgram::FastPath::None && (plan.profitable() ? jit.compile(plan) : jit.compile(runtime));
        const formula::ExpressionOrbitPlan* activePlan = plan.profitable() ? &plan : nullptr;

        for (formula::ExpressionColoring coloring : {formula::ExpressionColoring::Raw, formula::ExpressionColoring::Feather, formula::ExpressionColoring::OrbitTrap}) {
            std::vector<float> reference;
            if (!renderExpressionColorReference(runtime, test.fixed, test.pixel, test.centerRe, test.centerIm, test.scale, width, height, test.mxit, test.bailout, coloring, reference)) {
                ++failures;
                continue;
            }
            struct Path {
                bool vector;
                bool refill;
                bool powerSimd;
                bool useJit;
            };
            const Path paths[] = {{false, true, false, false}, {true, true, true, false}, {true, false, true, false}, {true, true, true, true}, {true, false, true, true}};
            for (const Path& path : paths) {
                if (path.useJit && !jitAvailable) continue;
                std::vector<float> actual((size_t)width * height, EMPTYPIXEL);
                Mandel renderer(width, height, test.mxit, 1, actual.data());
                mpf_t centerRe, centerIm, scale;
                mpf_init_set_d(centerRe, test.centerRe);
                mpf_init_set_d(centerIm, test.centerIm);
                mpf_init_set_d(scale, test.scale);
                _putenv_s("MANDEL_EXPR_VECTOR", path.vector ? "1" : "0");
                _putenv_s("MANDEL_EXPR_HYBRID_REFILL", path.refill ? "1" : "0");
                _putenv_s("MANDEL_EXPR_POWER_SIMD", path.powerSimd ? "1" : "0");
                bool okay = renderer.ComputeExpression(centerRe, centerIm, scale, runtime, test.fixed, test.pixel, test.mxit, test.bailout, coloring, path.useJit ? &jit : nullptr, activePlan);
                mpf_clears(centerRe, centerIm, scale, (mpf_ptr)0);
                ++comparisons;
                if (!okay) {
                    ++failures;
                    continue;
                }
                for (size_t index = 0; index < actual.size(); ++index) {
                    if (std::memcmp(&actual[index], &reference[index], sizeof(float)) != 0) ++mismatches;
                    if (coloring == formula::ExpressionColoring::OrbitTrap) {
                        if (!(actual[index] >= 0.0f) || actual[index] == EMPTYPIXEL || actual[index] == -2.0f) ++sentinelCollisions;
                    } else if (coloring == formula::ExpressionColoring::Feather) {
                        if (actual[index] != -2.0f && !(actual[index] >= 0.0f)) ++sentinelCollisions;
                    }
                }
            }
        }
    }

    if (expressionColoringFromMethod(0, false) != formula::ExpressionColoring::Raw || expressionColoringFromMethod(0, true) != formula::ExpressionColoring::Smooth || expressionColoringFromMethod(ColoringMethod::EXTERIOR_DIST_EST, true) != formula::ExpressionColoring::Distance || expressionColoringFromMethod(ColoringMethod::STRIPE_AVERAGE, false) != formula::ExpressionColoring::Feather || expressionColoringFromMethod(ColoringMethod::ORBIT_TRAP, false) != formula::ExpressionColoring::OrbitTrap) ++failures;
    for (int index = 0; index < 7; ++index) {
        bool genericExpected = index == 0 || index == 2 || index == 5;
        bool powerExpected = genericExpected || index == 1;
        if (expressionColoringIndexSupported(index, false) != genericExpected || expressionColoringIndexSupported(index, true) != powerExpected) ++failures;
    }

    auto median = [](std::vector<double> values) {
        std::sort(values.begin(), values.end());
        return values[values.size() / 2];
    };
    formula::ExpressionProgram profileSource;
    formula::ExpressionProgram profileRuntime;
    formula::ExpressionOrbitPlan profilePlan;
    formula::ExpressionContext profileFixed;
    formula::ExpressionError profileError;
    double medians[3] = {};
    const formula::ExpressionColoring profileColorings[] = {formula::ExpressionColoring::Raw, formula::ExpressionColoring::Feather, formula::ExpressionColoring::OrbitTrap};
    if (!profileSource.compile("0.5*z+0.1*sin(c)", &profileError) || !profileSource.specialize(profileFixed, FormulaParameter::C, profileRuntime, &profileError) || !profilePlan.build(profileRuntime, &profileError)) {
        ++failures;
    } else {
        formula::ExpressionJit4 profileJit;
        bool profileJitAvailable = profilePlan.profitable() ? profileJit.compile(profilePlan) : profileJit.compile(profileRuntime);
        for (int colorIndex = 0; colorIndex < 3; ++colorIndex) {
            std::vector<double> times;
            std::vector<float> frame(64 * 48, EMPTYPIXEL);
            for (int sample = 0; sample < 5; ++sample) {
                Mandel renderer(64, 48, 300, 1, frame.data());
                mpf_t re, im, scale;
                mpf_init_set_ui(re, 0);
                mpf_init_set_ui(im, 0);
                mpf_init_set_ui(scale, 1);
                auto begin = Clock::now();
                bool okay = renderer.ComputeExpression(re, im, scale, profileRuntime, profileFixed, FormulaParameter::C, 300, 4.0, profileColorings[colorIndex], profileJitAvailable ? &profileJit : nullptr, profilePlan.profitable() ? &profilePlan : nullptr);
                times.push_back(since(begin));
                mpf_clears(re, im, scale, (mpf_ptr)0);
                if (!okay) ++failures;
            }
            medians[colorIndex] = median(times);
        }
    }

    _putenv_s("MANDEL_EXPR_VECTOR", "");
    _putenv_s("MANDEL_EXPR_HYBRID_REFILL", "");
    _putenv_s("MANDEL_EXPR_POWER_SIMD", "");
    failures += mismatches != 0;
    failures += sentinelCollisions != 0;
    printf("  scalar/AVX2/Hybrid/JIT comparisons=%d mismatches=%d sentinels=%d\n", comparisons, mismatches, sentinelCollisions);
    printf("  repeated medians raw/feather/trap %.4f/%.4f/%.4f s"
           " (colored costs reported, not gated)\n",
           medians[0], medians[1], medians[2]);
    printf("  => %s\n\n", failures == 0 ? "PASS" : "CHECK (expression coloring failure)");
    return failures == 0 ? 0 : 1;
}

static int runExpressionOrbitCase() {
    using formula::Complex;
    using formula::ExpressionContext;
    using formula::ExpressionError;
    using formula::ExpressionOrbitEvaluation;
    using formula::ExpressionOrbitSnapshot;
    using formula::ExpressionProgram;

    int failures = 0;
    auto close = [](Complex a, Complex b, double eps = 1e-12) { return std::abs(a - b) <= eps * std::max(1.0, std::abs(b)); };
    auto makeSnapshot = [&](const char* source, FormulaParameter pixelParameter, const ExpressionContext& fixed, double bailout, ExpressionOrbitSnapshot& snapshot) {
        ExpressionProgram sourceProgram, runtime;
        ExpressionError error;
        if (!sourceProgram.compile(source, &error) || !sourceProgram.specialize(fixed, pixelParameter, runtime, &error)) {
            printf("  orbit compile failed @ %zu: %s  [%s]\n", error.position, error.message.c_str(), source);
            ++failures;
            return false;
        }
        snapshot = {};
        snapshot.program = std::move(runtime);
        snapshot.fixed = fixed;
        snapshot.pixelParameter = pixelParameter;
        snapshot.bailout = bailout;
        return true;
    };
    auto evaluate = [&](const ExpressionOrbitSnapshot& snapshot, Complex pixel, int maxIterations, ExpressionOrbitEvaluation& evaluation) {
        if (!formula::evaluateExpressionOrbit(snapshot, pixel, maxIterations, evaluation)) {
            ++failures;
            return false;
        }
        return true;
    };

    ExpressionContext fixed;
    ExpressionOrbitSnapshot quadratic;
    ExpressionOrbitEvaluation evaluation;
    if (makeSnapshot("z*z+c", FormulaParameter::C, fixed, 4.0, quadratic) && evaluate(quadratic, {1.0, 0.0}, 10, evaluation)) {
        const Complex expected[] = {0.0, 1.0, 2.0, 5.0};
        if (!evaluation.escaped || evaluation.iterations != 3 || evaluation.points.size() != std::size(expected)) {
            ++failures;
        } else {
            for (size_t i = 0; i < std::size(expected); ++i)
                if (!close(evaluation.points[i], expected[i])) ++failures;
        }
    }

    fixed = {};
    fixed.z0 = {99.0, -50.0};
    fixed.c = {-0.25, 0.0};
    ExpressionOrbitSnapshot z0Plane;
    if (makeSnapshot("z*z+c", FormulaParameter::InitialZ, fixed, 10.0, z0Plane) && evaluate(z0Plane, {0.5, 0.0}, 3, evaluation)) {
        const Complex expected[] = {0.5, 0.0, -0.25, -0.1875};
        if (evaluation.escaped || evaluation.iterations != 3 || evaluation.points.size() != std::size(expected)) {
            ++failures;
        } else {
            for (size_t i = 0; i < std::size(expected); ++i)
                if (!close(evaluation.points[i], expected[i])) ++failures;
        }
    }

    fixed = {};
    ExpressionOrbitSnapshot iterationDependent;
    if (makeSnapshot("z+c+n", FormulaParameter::C, fixed, 100.0, iterationDependent) && evaluate(iterationDependent, {1.0, 0.0}, 3, evaluation)) {
        const Complex expected[] = {0.0, 1.0, 3.0, 6.0};
        if (evaluation.points.size() != std::size(expected)) {
            ++failures;
        } else {
            for (size_t i = 0; i < std::size(expected); ++i)
                if (!close(evaluation.points[i], expected[i])) ++failures;
        }
    }

    fixed = {};
    fixed.z0 = {0.2, -0.1};
    ExpressionOrbitSnapshot transcendental;
    Complex transcendentalPixel{0.1, 0.2};
    if (makeSnapshot("sin(z)+c", FormulaParameter::C, fixed, 100.0, transcendental) && evaluate(transcendental, transcendentalPixel, 4, evaluation)) {
        Complex expected = fixed.z0;
        if (evaluation.points.size() != 5 || !close(evaluation.points[0], expected)) {
            ++failures;
        } else {
            for (int n = 0; n < 4; ++n) {
                expected = std::sin(expected) + transcendentalPixel;
                if (!close(evaluation.points[(size_t)n + 1], expected)) ++failures;
            }
        }
    }

    fixed = {};
    ExpressionOrbitSnapshot initialEscape;
    if (makeSnapshot("z+c", FormulaParameter::InitialZ, fixed, 4.0, initialEscape) && evaluate(initialEscape, {5.0, 0.0}, 10, evaluation) && (!evaluation.escaped || evaluation.iterations != 0 || evaluation.points.size() != 1 || evaluation.points[0] != Complex{5.0, 0.0})) ++failures;

    fixed = {};
    ExpressionOrbitSnapshot domain;
    if (makeSnapshot("log(z)", FormulaParameter::InitialZ, fixed, 4.0, domain) && evaluate(domain, {}, 10, evaluation)) {
        if (!evaluation.escaped || evaluation.iterations != 1 || evaluation.points.size() != 2 || (std::isfinite(evaluation.points[1].real()) && std::isfinite(evaluation.points[1].imag()))) ++failures;
    }

    fixed = {};
    ExpressionOrbitSnapshot cancellable;
    if (makeSnapshot("z+1", FormulaParameter::C, fixed, 1e100, cancellable)) {
        int checks = 0;
        if (!formula::evaluateExpressionOrbit(cancellable, {}, 100, evaluation, [&checks] { return ++checks > 3; }) || !evaluation.cancelled || evaluation.escaped || evaluation.iterations != 3 || evaluation.points.size() != 4) ++failures;
    }

    auto requestSnapshot = std::make_shared<const ExpressionOrbitSnapshot>(quadratic);
    quadratic = iterationDependent;
    if (evaluate(*requestSnapshot, {1.0, 0.0}, 10, evaluation) && (evaluation.points.size() != 4 || !close(evaluation.points.back(), {5.0, 0.0}))) ++failures;

    auto waitForOrbit = [](OrbitWorker& worker, OrbitResult& result) {
        auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(2);
        while (std::chrono::steady_clock::now() < deadline) {
            if (worker.takeLatest(result)) return true;
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
        }
        return false;
    };
    mpf_t centerRe, centerIm, scale;
    mpf_init_set_ui(centerRe, 0);
    mpf_init_set_ui(centerIm, 0);
    mpf_init_set_ui(scale, 1);
    FormulaContext customFormula = expressionFormula();
    customFormula.slice.pixel = FormulaParameter::C;
    {
        OrbitWorker worker;
        OrbitResult workerResult;
        auto lifetimeSnapshot = std::make_shared<const ExpressionOrbitSnapshot>(cancellable);
        worker.request(centerRe, centerIm, scale, 1, 1, 3, 3, 3, customFormula, lifetimeSnapshot);
        lifetimeSnapshot.reset();
        cancellable = iterationDependent;
        if (!waitForOrbit(worker, workerResult) || workerResult.iterations != 3 || workerResult.points.size() != 4 || !close({workerResult.points.back().re, workerResult.points.back().im}, {3.0, 0.0})) { ++failures; }

        ExpressionOrbitSnapshot plusOne, plusTwo, identity;
        fixed = {};
        makeSnapshot("z+1", FormulaParameter::C, fixed, 1e100, plusOne);
        makeSnapshot("z+2", FormulaParameter::C, fixed, 1e100, plusTwo);
        makeSnapshot("z", FormulaParameter::C, fixed, 1e100, identity);
        worker.request(centerRe, centerIm, scale, 1, 1, 3, 3, 2048, customFormula, std::make_shared<const ExpressionOrbitSnapshot>(plusOne));
        worker.request(centerRe, centerIm, scale, 1, 1, 3, 3, 2, customFormula, std::make_shared<const ExpressionOrbitSnapshot>(plusTwo));
        if (!waitForOrbit(worker, workerResult) || workerResult.iterations != 2 || workerResult.points.size() != 3 || !close({workerResult.points.back().re, workerResult.points.back().im}, {4.0, 0.0})) { ++failures; }

        worker.request(centerRe, centerIm, scale, 1, 1, 3, 3, 100000, customFormula, std::make_shared<const ExpressionOrbitSnapshot>(identity));
        if (!waitForOrbit(worker, workerResult) || workerResult.iterations != 2048 || workerResult.points.size() != 2049) ++failures;
    }
    mpf_clears(centerRe, centerIm, scale, (mpf_ptr)0);

    {
        OrbitThumbnailPixel pixel;
        if (!classifyOrbitThumbnailPixel(nullptr, {2.0, 0.0}, 2, pixel) || !pixel.escaped || pixel.iterations != 2) ++failures;
        if (!classifyOrbitThumbnailPixel(nullptr, {2.0, 0.0}, 1, pixel) || pixel.escaped || pixel.iterations != 1) ++failures;

        ExpressionOrbitSnapshot finalEscape;
        fixed = {};
        if (makeSnapshot("z+1", FormulaParameter::C, fixed, 1.5, finalEscape)) {
            if (!classifyOrbitThumbnailPixel(&finalEscape, {}, 2, pixel) || !pixel.escaped || pixel.iterations != 2) ++failures;
            if (!classifyOrbitThumbnailPixel(&finalEscape, {}, 1, pixel) || pixel.escaped || pixel.iterations != 1) ++failures;
        }
        if (!classifyOrbitThumbnailPixel(&z0Plane, {0.5, 0.0}, 3, pixel) || pixel.escaped || pixel.iterations != 3) ++failures;
        if (!classifyOrbitThumbnailPixel(&initialEscape, {5.0, 0.0}, 10, pixel) || !pixel.escaped || pixel.iterations != 0) ++failures;
        if (!classifyOrbitThumbnailPixel(&domain, {}, 10, pixel) || !pixel.escaped || pixel.iterations != 1) ++failures;
        ExpressionOrbitSnapshot parameterAndIteration;
        fixed = {};
        fixed.parameters[0] = {2.0, 0.0};
        if (makeSnapshot("z+p0+n", FormulaParameter::C, fixed, 4.5, parameterAndIteration) && (!classifyOrbitThumbnailPixel(&parameterAndIteration, {}, 3, pixel) || !pixel.escaped || pixel.iterations != 2)) ++failures;

        auto waitForThumbnail = [](OrbitThumbnailWorker& worker, OrbitThumbnailResult& result) {
            auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(5);
            while (std::chrono::steady_clock::now() < deadline) {
                if (worker.takeLatest(result)) return true;
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
            }
            return false;
        };
        OrbitThumbnailWorker thumbnailWorker;
        ExpressionOrbitSnapshot slow, latest;
        fixed = {};
        makeSnapshot("sin(z)+c", FormulaParameter::C, fixed, 100.0, slow);
        fixed.c = {-0.25, 0.0};
        makeSnapshot("z*z+c", FormulaParameter::InitialZ, fixed, 8.0, latest);
        uint64_t staleGeneration = thumbnailWorker.request(96, 64, -2.5, 1.0, -1.2, 1.2, 160, std::make_shared<const ExpressionOrbitSnapshot>(slow));
        uint64_t latestGeneration = thumbnailWorker.request(32, 24, -2.5, 1.0, -1.2, 1.2, 40, std::make_shared<const ExpressionOrbitSnapshot>(latest));
        OrbitThumbnailResult thumbnail;
        if (latestGeneration <= staleGeneration || !waitForThumbnail(thumbnailWorker, thumbnail) || thumbnail.generation != latestGeneration || !thumbnail.expression || thumbnail.pixelParameter != FormulaParameter::InitialZ || thumbnail.pixels.size() != (size_t)32 * 24 * 3) ++failures;

        uint64_t mandelGeneration = thumbnailWorker.request(16, 12, -2.5, 1.0, -1.2, 1.2, 20, nullptr);
        if (!waitForThumbnail(thumbnailWorker, thumbnail) || thumbnail.generation != mandelGeneration || thumbnail.expression || thumbnail.pixels.size() != (size_t)16 * 12 * 3) ++failures;
    }

    printf("=== expression orbit overlay\n");
    printf("  c-plane/z0/n/transcendental/escape/domain/cancel/lifetime/generation/cap/thumbnail\n");
    printf("  => %s\n\n", failures == 0 ? "PASS" : "CHECK (expression orbit failure)");
    return failures == 0 ? 0 : 1;
}

static int runExpressionReferenceCase() {
    using formula::Complex;
    using formula::ExpressionContext;
    using formula::ExpressionError;
    using formula::ExpressionOracle;
    using formula::ExpressionOracleContext;
    using formula::ExpressionOracleOperation;
    using formula::ExpressionOracleTrace;
    using formula::ExpressionOracleTraceNode;
    using formula::ExpressionOrbitPlan;
    using formula::ExpressionProgram;
    using formula::ExpressionReferenceBuildRequest;
    using formula::ExpressionReferenceBuildStatus;
    using formula::ExpressionReferenceOrbitResult;
    using formula::ExpressionReferencePrecisionPolicy;
    using formula::ExpressionReferenceTapeNode;
    using formula::MpfrComplex;
    using formula::ScaledComplexShadow;
    using formula::ScaledRealShadow;

    int failures = 0;
    int orbitSamples = 0;
    int companionChecks = 0;
    int branchChecks = 0;
    int defectChecks = 0;

    auto compile = [&](ExpressionProgram& program, const std::string& source) {
        ExpressionError error;
        if (!program.compile(source, &error)) {
            printf("  reference compile failed @ %zu: %s [%s]\n", error.position, error.message.c_str(), source.c_str());
            ++failures;
            return false;
        }
        return true;
    };
    auto specialize = [&](const ExpressionProgram& canonical, const ExpressionContext& fixed, FormulaParameter pixel, ExpressionProgram& runtime) {
        ExpressionError error;
        if (!canonical.specialize(fixed, pixel, runtime, &error)) {
            printf("  reference specialize failed: %s [%s]\n", error.message.c_str(), canonical.source().c_str());
            ++failures;
            return false;
        }
        return true;
    };
    auto buildDecimal = [&](const ExpressionProgram& canonical, const ExpressionProgram& runtime, FormulaParameter pixel, const ExpressionContext& fixed, const std::string& real, const std::string& imaginary, int iterations, mpfr_prec_t bits, ExpressionReferenceOrbitResult& result, double bailout = 1e100, mpfr_prec_t guard = 0) {
        ExpressionReferenceBuildRequest request;
        request.canonicalProgram = &canonical;
        request.runtimeProgram = &runtime;
        request.pixelParameter = pixel;
        request.fixed = fixed;
        request.center.realDecimal = real;
        request.center.imaginaryDecimal = imaginary;
        request.bailout = bailout;
        request.maxIterations = iterations;
        request.precision.requestedBits = bits;
        request.precision.minimumBits = 53;
        request.precision.guardBits = guard;
        request.precision.maximumBits = 8192;
        if (!formula::buildExpressionReferenceOrbit(request, result)) {
            printf("  reference build failed: %s [%s]\n", result.error.c_str(), canonical.source().c_str());
            ++failures;
            return false;
        }
        return true;
    };
    auto sameRealShadow = [](const ScaledRealShadow& left, const ScaledRealShadow& right) {
        uint64_t leftBits = 0, rightBits = 0;
        std::memcpy(&leftBits, &left.mantissa, sizeof(leftBits));
        std::memcpy(&rightBits, &right.mantissa, sizeof(rightBits));
        return leftBits == rightBits && left.exponent == right.exponent;
    };
    auto sameShadow = [&](const ScaledComplexShadow& left, const ScaledComplexShadow& right) { return sameRealShadow(left.re, right.re) && sameRealShadow(left.im, right.im); };
    auto closeMpfr = [](mpfr_srcptr actual, mpfr_srcptr expected, int bits) {
        if (mpfr_nan_p(actual) || mpfr_nan_p(expected)) return mpfr_nan_p(actual) && mpfr_nan_p(expected);
        if (mpfr_inf_p(actual) || mpfr_inf_p(expected)) return mpfr_inf_p(actual) && mpfr_inf_p(expected) && mpfr_signbit(actual) == mpfr_signbit(expected);
        if (mpfr_equal_p(actual, expected)) return true;
        if (mpfr_zero_p(expected)) return mpfr_zero_p(actual);
        mpfr_prec_t precision = std::max(mpfr_get_prec(actual), mpfr_get_prec(expected));
        mpfr_t difference, tolerance;
        mpfr_init2(difference, precision);
        mpfr_init2(tolerance, precision);
        mpfr_sub(difference, actual, expected, MPFR_RNDN);
        mpfr_abs(difference, difference, MPFR_RNDN);
        mpfr_abs(tolerance, expected, MPFR_RNDN);
        mpfr_mul_2si(tolerance, tolerance, -bits, MPFR_RNDU);
        bool close = mpfr_cmp(difference, tolerance) <= 0;
        mpfr_clear(difference);
        mpfr_clear(tolerance);
        return close;
    };
    auto closeMpfrComplex = [&](const MpfrComplex& actual, const MpfrComplex& expected, int bits) { return closeMpfr(actual.re, expected.re, bits) && closeMpfr(actual.im, expected.im, bits); };
    auto configureOracle = [](ExpressionOracleContext& context, const ExpressionContext& fixed, FormulaParameter pixel, const std::string& real, const std::string& imaginary) {
        context.c.set(fixed.c.real(), fixed.c.imag());
        context.z0.set(fixed.z0.real(), fixed.z0.imag());
        for (size_t i = 0; i < fixed.parameters.size(); ++i) { context.parameters[i].set(fixed.parameters[i].real(), fixed.parameters[i].imag()); }
        MpfrComplex center(context.z.precision());
        if (!center.set(real, imaginary)) return false;
        if (pixel == FormulaParameter::C)
            context.c.set(center);
        else
            context.z0.set(center);
        context.z.set(context.z0);
        return true;
    };

    {
        MpfrComplex parsed(256);
        if (!parsed.set("0.1", "1e-500") || !mpfr_number_p(parsed.re) || !mpfr_number_p(parsed.im) || mpfr_zero_p(parsed.re) || mpfr_zero_p(parsed.im) || parsed.set("1e-400000000", "0") || parsed.set("1e400000000", "0")) {
            printf("  finite/inexact decimal parse range handling failed\n");
            ++failures;
        }
    }

    std::string deepCoordinate = "1.";
    deepCoordinate.append(499, '0');
    deepCoordinate.push_back('1');
    ExpressionProgram coordinateCanonical, coordinateRuntimeC;
    ExpressionProgram coordinateRuntimeZ0;
    ExpressionContext coordinateFixed;
    if (compile(coordinateCanonical, "z+c") && specialize(coordinateCanonical, coordinateFixed, FormulaParameter::C, coordinateRuntimeC) && specialize(coordinateCanonical, coordinateFixed, FormulaParameter::InitialZ, coordinateRuntimeZ0)) {
        ExpressionReferenceOrbitResult deepC, baseC, deepZ0, baseZ0;
        bool coordinateOk = buildDecimal(coordinateCanonical, coordinateRuntimeC, FormulaParameter::C, coordinateFixed, deepCoordinate, "1e-500", 1, 1800, deepC) && buildDecimal(coordinateCanonical, coordinateRuntimeC, FormulaParameter::C, coordinateFixed, "1", "0", 1, 1800, baseC) && buildDecimal(coordinateCanonical, coordinateRuntimeZ0, FormulaParameter::InitialZ, coordinateFixed, deepCoordinate, "1e-500", 1, 1800, deepZ0, 4.0) && buildDecimal(coordinateCanonical, coordinateRuntimeZ0, FormulaParameter::InitialZ, coordinateFixed, "1", "0", 1, 1800, baseZ0, 4.0);
        if (coordinateOk) {
            MpfrComplex reconstructedDeep(2048);
            MpfrComplex reconstructedBase(2048);
            if (!formula::reconstructMpfrFromShadows(reconstructedDeep, deepC.pixel, deepC.pixelDefect) || !formula::reconstructMpfrFromShadows(reconstructedBase, baseC.pixel, baseC.pixelDefect) || mpfr_cmp(reconstructedDeep.re, reconstructedBase.re) == 0 || mpfr_zero_p(reconstructedDeep.im) || !deepC.pixelDefect.re.isFinite() || deepC.pixel.im.exponent > -1000 || deepZ0.initialZDefect.re.isZero() || deepZ0.initialZ.im.exponent > -1000 || mpfr_get_d(reconstructedDeep.re, MPFR_RNDN) != mpfr_get_d(reconstructedBase.re, MPFR_RNDN) || mpfr_get_d(reconstructedDeep.im, MPFR_RNDN) != 0.0 || mpfr_get_d(reconstructedDeep.im, MPFR_RNDN) != mpfr_get_d(reconstructedBase.im, MPFR_RNDN)) {
                printf("  e500 coordinate shadow preservation failed\n");
                ++failures;
            }
        }
    }

    if (coordinateRuntimeC.valid()) {
        auto checkMpfInput = [&](const char* label, mpf_srcptr real, mpf_srcptr imaginary) {
            auto significantBits = [](mpf_srcptr value) {
                const mp_size_t limbCount = mpf_size(value);
                if (limbCount == 0) return mpfr_prec_t{0};
                const size_t highBit = mpn_sizeinbase(value->_mp_d, limbCount, 2);
                const mp_bitcnt_t lowBit = mpn_scan1(value->_mp_d, 0);
                return static_cast<mpfr_prec_t>(std::max<size_t>(1, highBit - lowBit));
            };
            const mpfr_prec_t required = static_cast<mpfr_prec_t>(std::max(significantBits(real), significantBits(imaginary)));
            ExpressionReferenceBuildRequest request;
            request.canonicalProgram = &coordinateCanonical;
            request.runtimeProgram = &coordinateRuntimeC;
            request.pixelParameter = FormulaParameter::C;
            request.center.realMpf = real;
            request.center.imaginaryMpf = imaginary;
            request.maxIterations = 1;
            request.bailout = 1e100;
            request.precision.minimumBits = 53;
            request.precision.guardBits = 0;
            request.precision.maximumBits = 4096;
            ExpressionReferenceOrbitResult result;
            MpfrComplex reconstructed(required);
            if (!formula::buildExpressionReferenceOrbit(request, result) || result.precision != required || !formula::reconstructMpfrFromShadows(reconstructed, result.pixel, result.pixelDefect) || mpfr_cmp_f(reconstructed.re, real) != 0 || mpfr_cmp_f(reconstructed.im, imaginary) != 0) {
                printf("  exact-significand mpf transfer failed [%s]\n", label);
                ++failures;
            }
            request.precision.maximumBits = required - 1;
            ExpressionReferenceOrbitResult rejected;
            if (formula::buildExpressionReferenceOrbit(request, rejected) || rejected.status != ExpressionReferenceBuildStatus::PrecisionOutOfRange) {
                printf("  significand precision policy failed [%s]\n", label);
                ++failures;
            }
        };

        mpf_t carry, zero;
        mpz_t integer;
        mpf_init2(carry, 64);
        mpf_init2(zero, 64);
        mpz_init(integer);
        mpz_set_ui(integer, 1);
        mpz_mul_2exp(integer, integer, 65);
        mpz_sub_ui(integer, integer, 1);
        mpf_set_z(carry, integer);
        if (mpf_size(carry) != 2 || mpf_get_prec(carry) >= mpf_size(carry) * GMP_NUMB_BITS) ++failures;
        checkMpfInput("dense carry limb", carry, zero);

        mpf_t rawPrecision;
        mpf_init2(rawPrecision, 256);
        mpz_set_ui(integer, 1);
        mpz_mul_2exp(integer, integer, 192);
        mpz_add_ui(integer, integer, 1);
        mpf_set_z(rawPrecision, integer);
        const mp_bitcnt_t originalPrecision = mpf_get_prec(rawPrecision);
        mpf_set_prec_raw(rawPrecision, 64);
        if (mpf_size(rawPrecision) != 4 || mpf_get_prec(rawPrecision) != 64) ++failures;
        checkMpfInput("raw precision live significand", rawPrecision, zero);
        mpf_set_prec_raw(rawPrecision, originalPrecision);

        mpf_clear(rawPrecision);
        mpz_clear(integer);
        mpf_clears(carry, zero, (mpf_ptr)0);
    }

    struct OrbitFormula {
        const char* source;
        FormulaParameter pixel;
        const char* real;
        const char* imaginary;
    };
    const OrbitFormula orbitFormulas[] = {{"z*z+c", FormulaParameter::C, "-0.2", "0.3"}, {"sin(z)+c", FormulaParameter::C, "0.1", "-0.15"}, {"exp(z/4)+c/8", FormulaParameter::C, "-0.3", "0.1"}, {"(z*z+c)/(z+2)", FormulaParameter::C, "0.2", "0.1"}, {"pow(z+1.5,p0)+c", FormulaParameter::C, "-0.15", "0.12"}, {"z+c+n/10", FormulaParameter::C, "0.05", "-0.03"}, {"sqr(complex(abs(real(z)),abs(imag(z))))+c", FormulaParameter::C, "-0.25", "0.2"}};
    for (const OrbitFormula& item : orbitFormulas) {
        ExpressionContext fixed;
        fixed.z0 = {0.125, -0.0625};
        fixed.parameters[0] = {0.5, 0.2};
        ExpressionProgram canonical, runtime;
        if (!compile(canonical, item.source) || !specialize(canonical, fixed, item.pixel, runtime)) continue;
        ExpressionReferenceOrbitResult reference;
        if (!buildDecimal(canonical, runtime, item.pixel, fixed, item.real, item.imaginary, 7, 256, reference)) continue;
        if (!reference.valid || reference.defectPending || reference.programSemanticHash != runtime.semanticHash() || reference.sampleCount != reference.samples.size()) {
            ++failures;
            continue;
        }

        ExpressionOracleContext same(reference.precision);
        ExpressionOracleContext higher(512);
        if (!configureOracle(same, fixed, item.pixel, item.real, item.imaginary) || !configureOracle(higher, fixed, item.pixel, item.real, item.imaginary)) {
            ++failures;
            continue;
        }
        MpfrComplex sameNext(reference.precision);
        MpfrComplex highNext(512);
        MpfrComplex reconstructed(reference.precision);
        for (size_t i = 0; i < reference.samples.size(); ++i) {
            same.iteration = static_cast<int>(i);
            higher.iteration = static_cast<int>(i);
            std::string sameError, highError;
            bool sameOk = ExpressionOracle::evaluate(runtime, same, sameNext, &sameError);
            bool highOk = ExpressionOracle::evaluate(runtime, higher, highNext, &highError);
            const auto& sample = reference.samples[i];
            const auto& root = reference.tape[static_cast<size_t>(sample.tapeOffset) + sample.rootNode];
            ++orbitSamples;
            if (!sameOk || !highOk || !formula::reconstructMpfrFromShadows(reconstructed, sample.next, sample.rootDefect) || !closeMpfrComplex(reconstructed, sameNext, 94) || !closeMpfrComplex(reconstructed, highNext, 90) || !sameShadow(root.output, sample.next) || !sameShadow(root.outputDefect, sample.rootDefect) || sample.tapeCount != runtime.instructionCount()) {
                printf("  root/tape mismatch at %zu [%s]\n", i, item.source);
                ++failures;
                break;
            }
            same.z.set(sameNext);
            higher.z.set(highNext);
        }
    }

    ExpressionProgram companionProgram;
    if (compile(companionProgram, "sin(z)+cos(z)+sinh(z)+cosh(z)+tan(z)+tanh(z)")) {
        auto companionSource = [](ExpressionOracleOperation operation) {
            switch (operation) {
            case ExpressionOracleOperation::Sin:
            case ExpressionOracleOperation::Tan: return "cos(z)";
            case ExpressionOracleOperation::Cos: return "sin(z)";
            case ExpressionOracleOperation::Sinh:
            case ExpressionOracleOperation::Tanh: return "cosh(z)";
            case ExpressionOracleOperation::Cosh: return "sinh(z)";
            default: return static_cast<const char*>(nullptr);
            }
        };
        ExpressionOracleContext context(512);
        context.z.set("0.3", "0.2");
        MpfrComplex output(512);
        ExpressionOracleTrace trace;
        std::string error;
        if (!ExpressionOracle::evaluateTrace(companionProgram, context, output, trace, &error) || trace.nodes.empty() || !mpfr_equal_p(trace.nodes.back().output.re, output.re) || !mpfr_equal_p(trace.nodes.back().output.im, output.im)) {
            printf("  trace root mismatch: %s\n", error.c_str());
            ++failures;
        } else {
            for (const auto& node : trace.nodes) {
                const char* source = companionSource(node.operation);
                if (!source) continue;
                ExpressionProgram independent;
                if (!compile(independent, source)) continue;
                MpfrComplex expected(512);
                std::string independentError;
                ++companionChecks;
                if (!ExpressionOracle::evaluate(independent, context, expected, &independentError) || !(node.flags & formula::OracleTraceHasCompanion) || !closeMpfrComplex(node.auxiliary, expected, 500)) {
                    printf("  nonlinear companion mismatch [%s]\n", source);
                    ++failures;
                }
            }
        }

        ExpressionContext fixed;
        ExpressionProgram runtime;
        ExpressionReferenceOrbitResult tapeReference;
        if (specialize(companionProgram, fixed, FormulaParameter::InitialZ, runtime) && buildDecimal(companionProgram, runtime, FormulaParameter::InitialZ, fixed, "0.3", "0.2", 1, 512, tapeReference)) {
            for (const auto& node : tapeReference.tape) {
                const char* source = companionSource(node.operation);
                if (!source) continue;
                ExpressionProgram independent;
                if (!compile(independent, source)) continue;
                MpfrComplex expected(512);
                MpfrComplex reconstructed(512);
                std::string independentError;
                ++companionChecks;
                if (!ExpressionOracle::evaluate(independent, context, expected, &independentError) || !formula::reconstructMpfrFromShadows(reconstructed, node.auxiliary, node.auxiliaryDefect) || !closeMpfrComplex(reconstructed, expected, 94)) {
                    printf("  compact companion mismatch [%s]\n", source);
                    ++failures;
                }
            }
        }
    }

    struct AsymptoticCompanionCase {
        const char* source;
        const char* real;
        const char* imaginary;
        const char* logScale;
        const char* phase;
        ExpressionOracleOperation operation;
    };
    const AsymptoticCompanionCase asymptoticCases[] = {{"tan(z)", "1", "1e9", "-1e9", "2", ExpressionOracleOperation::Tan}, {"tanh(z)", "1e9", "0.25", "-1e9", "0.5", ExpressionOracleOperation::Tanh}};
    for (const auto& item : asymptoticCases) {
        ExpressionProgram canonical, runtime;
        ExpressionContext fixed;
        if (!compile(canonical, item.source) || !specialize(canonical, fixed, FormulaParameter::InitialZ, runtime)) continue;

        ExpressionOracleContext context(512);
        context.z.set(item.real, item.imaginary);
        MpfrComplex output(512);
        MpfrComplex expectedAuxiliary(512);
        expectedAuxiliary.set(item.logScale, item.phase);
        ExpressionOracleTrace trace;
        std::string error;
        const ExpressionOracleTraceNode* tracedNode = nullptr;
        if (ExpressionOracle::evaluateTrace(runtime, context, output, trace, &error)) {
            for (const auto& node : trace.nodes)
                if (node.operation == item.operation) tracedNode = &node;
        }
        ++companionChecks;
        if (!tracedNode || !mpfr_number_p(output.re) || !mpfr_number_p(output.im) || !(tracedNode->flags & formula::OracleTraceHasCompanion) || !(tracedNode->flags & formula::OracleTraceHasDenominator) || !(tracedNode->flags & formula::OracleTraceHasAsymptoticLogPhase) || tracedNode->clearance != formula::ExpressionOraclePointClearance::ClearAtPoint || !mpfr_number_p(tracedNode->auxiliary.re) || !mpfr_number_p(tracedNode->auxiliary.im) || !closeMpfrComplex(tracedNode->auxiliary, expectedAuxiliary, 480)) {
            printf("  asymptotic companion trace failed [%s]\n", item.source);
            ++failures;
            continue;
        }

        ExpressionReferenceOrbitResult reference;
        if (!buildDecimal(canonical, runtime, FormulaParameter::InitialZ, fixed, item.real, item.imaginary, 1, 512, reference, 1e100)) continue;
        const ExpressionReferenceTapeNode* compactNode = nullptr;
        for (const auto& node : reference.tape)
            if (node.operation == item.operation) compactNode = &node;
        MpfrComplex reconstructedOutput(512);
        MpfrComplex reconstructedAuxiliary(512);
        ++companionChecks;
        if (!compactNode || !(compactNode->flags & formula::OracleTraceHasAsymptoticLogPhase) || !formula::reconstructMpfrFromShadows(reconstructedOutput, compactNode->output, compactNode->outputDefect) || !formula::reconstructMpfrFromShadows(reconstructedAuxiliary, compactNode->auxiliary, compactNode->auxiliaryDefect) || !mpfr_number_p(reconstructedOutput.re) || !mpfr_number_p(reconstructedOutput.im) || !closeMpfrComplex(reconstructedAuxiliary, expectedAuxiliary, 94)) {
            printf("  asymptotic compact companion failed [%s]\n", item.source);
            ++failures;
        }
    }

    auto checkNonFiniteTangentTrace = [&](const char* source, ExpressionOracleOperation operation, bool folded) {
        ExpressionProgram canonical, runtime;
        ExpressionContext fixed;
        if (folded) fixed.parameters[0] = {0.0, 0.0};
        if (!compile(canonical, source)) return;
        const ExpressionProgram* program = &canonical;
        if (folded) {
            if (!specialize(canonical, fixed, FormulaParameter::InitialZ, runtime)) return;
            program = &runtime;
        }

        ExpressionOracleContext context(256);
        context.z.set("0.25", "0.125");
        if (!folded) {
            if (operation == ExpressionOracleOperation::Tan)
                mpfr_set_inf(context.z.im, 1);
            else
                mpfr_set_inf(context.z.re, 1);
        }
        MpfrComplex output(256);
        ExpressionOracleTrace trace;
        std::string error;
        const bool defined = ExpressionOracle::evaluateTrace(*program, context, output, trace, &error);
        const ExpressionOracleTraceNode* tangentNode = nullptr;
        bool sawFoldedNonFinite = !folded;
        for (const auto& node : trace.nodes) {
            if (node.operation == operation) tangentNode = &node;
            if (folded && node.operation == ExpressionOracleOperation::Constant && (node.flags & formula::OracleTraceUndefined) && (!mpfr_number_p(node.output.re) || !mpfr_number_p(node.output.im))) sawFoldedNonFinite = true;
        }
        if (defined || !tangentNode || !sawFoldedNonFinite || !(tangentNode->flags & formula::OracleTraceUndefined) || (tangentNode->flags & formula::OracleTraceHasAsymptoticLogPhase) || tangentNode->clearance != formula::ExpressionOraclePointClearance::NonFiniteAtPoint) {
            printf("  nonfinite tangent trace classification failed [%s]\n", source);
            ++failures;
        }
    };
    checkNonFiniteTangentTrace("tan(z)", ExpressionOracleOperation::Tan, false);
    checkNonFiniteTangentTrace("tanh(z)", ExpressionOracleOperation::Tanh, false);
    checkNonFiniteTangentTrace("tan(z+complex(0,log(p0)))", ExpressionOracleOperation::Tan, true);
    checkNonFiniteTangentTrace("tanh(z+log(p0))", ExpressionOracleOperation::Tanh, true);

    ExpressionProgram defectCanonical, defectRuntime;
    ExpressionContext defectFixed;
    const char* defectReal = "0.123456789012345678901234567890123456789";
    if (compile(defectCanonical, "z+c") && specialize(defectCanonical, defectFixed, FormulaParameter::C, defectRuntime)) {
        for (mpfr_prec_t bits : {256L, 1024L, 2048L}) {
            ExpressionReferenceOrbitResult reference;
            if (!buildDecimal(defectCanonical, defectRuntime, FormulaParameter::C, defectFixed, defectReal, "1e-500", 1, bits, reference)) continue;
            ExpressionOracleContext expected(reference.precision);
            configureOracle(expected, defectFixed, FormulaParameter::C, defectReal, "1e-500");
            MpfrComplex next(reference.precision);
            std::string error;
            expected.iteration = 0;
            if (!ExpressionOracle::evaluate(defectRuntime, expected, next, &error) || reference.samples.size() != 1) {
                ++failures;
                continue;
            }
            MpfrComplex reconstructed(reference.precision);
            const auto& sample = reference.samples[0];
            ++defectChecks;
            if (reference.precision != bits || sample.rootDefect.re.isZero() || sample.next.im.isZero() || sample.next.im.exponent > -1000 || !formula::reconstructMpfrFromShadows(reconstructed, sample.next, sample.rootDefect) || !closeMpfrComplex(reconstructed, next, 94)) {
                printf("  defect reconstruction failed at %ld bits\n", (long)bits);
                ++failures;
            }
        }
    }

    ExpressionProgram foldedDomainCanonical;
    if (compile(foldedDomainCanonical, "log(p0)+z")) {
        for (const Complex parameter : {Complex{0.0, 0.0}, Complex{2.0, 0.0}}) {
            ExpressionContext fixed;
            fixed.parameters[0] = parameter;
            ExpressionProgram runtime;
            if (!specialize(foldedDomainCanonical, fixed, FormulaParameter::InitialZ, runtime)) continue;
            ExpressionOracleContext context(256);
            context.z.set("0.5", "0");
            MpfrComplex output(256);
            ExpressionOracleTrace trace;
            std::string error;
            const bool expectedDefined = parameter.real() != 0.0;
            const bool defined = ExpressionOracle::evaluate(runtime, context, output, &error);
            std::string traceError;
            const bool traced = ExpressionOracle::evaluateTrace(runtime, context, output, trace, &traceError);
            bool constantUndefined = false;
            bool finiteConstant = false;
            for (const auto& node : trace.nodes) {
                if (node.operation != ExpressionOracleOperation::Constant) continue;
                constantUndefined = (node.flags & formula::OracleTraceUndefined) != 0;
                finiteConstant = mpfr_number_p(node.output.re) && mpfr_number_p(node.output.im);
            }
            if (runtime.instructionCount() != 3 || defined != expectedDefined || traced != expectedDefined || constantUndefined == expectedDefined || finiteConstant != expectedDefined) {
                printf("  folded constant domain handling failed\n");
                ++failures;
            }

            ExpressionReferenceOrbitResult reference;
            if (buildDecimal(foldedDomainCanonical, runtime, FormulaParameter::InitialZ, fixed, "0.5", "0", 1, 256, reference, 4.0)) {
                if (expectedDefined) {
                    if (reference.undefined || reference.escaped) ++failures;
                } else if (!reference.undefined || reference.escaped || reference.undefinedIteration != 1) {
                    printf("  folded domain reference escaped\n");
                    ++failures;
                }
            }
        }
    }

    ExpressionProgram branchCanonical, branchRuntime;
    ExpressionContext branchFixed;
    branchFixed.parameters[0] = {0.5, 0.0};
    if (compile(branchCanonical, "log(z)+sqrt(z)+arg(z)+pow(z,p0)+"
                                 "1/(z-z)+tan(z)+tanh(z)") &&
        specialize(branchCanonical, branchFixed, FormulaParameter::InitialZ, branchRuntime)) {
        ExpressionReferenceOrbitResult branch;
        if (buildDecimal(branchCanonical, branchRuntime, FormulaParameter::InitialZ, branchFixed, "-1", "-0", 1, 512, branch)) {
            bool sawLog = false, sawSqrt = false;
            bool sawArg = false, sawPower = false;
            bool sawDivide = false, sawTan = false;
            bool sawTanh = false;
            for (const auto& node : branch.tape) {
                auto branchNode = [&] {
                    ++branchChecks;
                    return (node.flags & formula::OracleTraceBranchSensitive) && (node.cut == formula::ExpressionOracleCutLocation::NegativeRealLowerLip || node.cut == formula::ExpressionOracleCutLocation::NegativeRealUpperLip) && node.certification == formula::ExpressionOracleCertification::PointOnlyNotCertified;
                };
                switch (node.operation) {
                case ExpressionOracleOperation::Log: sawLog = branchNode(); break;
                case ExpressionOracleOperation::Sqrt: sawSqrt = branchNode(); break;
                case ExpressionOracleOperation::Arg: sawArg = branchNode(); break;
                case ExpressionOracleOperation::Power: sawPower = branchNode(); break;
                case ExpressionOracleOperation::Divide:
                    ++branchChecks;
                    sawDivide = node.clearance == formula::ExpressionOraclePointClearance::ZeroAtPoint && (node.flags & formula::OracleTraceSingularPoint);
                    break;
                case ExpressionOracleOperation::Tan:
                    ++branchChecks;
                    sawTan = node.clearance == formula::ExpressionOraclePointClearance::ClearAtPoint && (node.flags & formula::OracleTraceHasDenominator);
                    break;
                case ExpressionOracleOperation::Tanh:
                    ++branchChecks;
                    sawTanh = node.clearance == formula::ExpressionOraclePointClearance::ClearAtPoint && (node.flags & formula::OracleTraceHasDenominator);
                    break;
                default: break;
                }
            }
            if (!branch.undefined || !(sawLog && sawSqrt && sawArg && sawPower && sawDivide && sawTan && sawTanh)) {
                printf("  branch/domain point metadata failed\n");
                ++failures;
            }
        }
    }

    ExpressionProgram simpleCanonical, simpleRuntimeC;
    ExpressionProgram simpleRuntimeZ0;
    ExpressionContext simpleFixed;
    if (compile(simpleCanonical, "z+1") && specialize(simpleCanonical, simpleFixed, FormulaParameter::C, simpleRuntimeC) && specialize(simpleCanonical, simpleFixed, FormulaParameter::InitialZ, simpleRuntimeZ0)) {
        ExpressionReferenceOrbitResult initialEscape;
        if (buildDecimal(simpleCanonical, simpleRuntimeZ0, FormulaParameter::InitialZ, simpleFixed, "5", "0", 10, 256, initialEscape, 4.0) && (!initialEscape.escaped || initialEscape.escapeIteration != 0 || initialEscape.sampleCount != 0)) ++failures;

        ExpressionReferenceBuildRequest cancellation;
        cancellation.canonicalProgram = &simpleCanonical;
        cancellation.runtimeProgram = &simpleRuntimeC;
        cancellation.pixelParameter = FormulaParameter::C;
        cancellation.center.realDecimal = "0";
        cancellation.maxIterations = 100;
        cancellation.bailout = 1e100;
        cancellation.precision.requestedBits = 256;
        cancellation.precision.minimumBits = 53;
        cancellation.precision.guardBits = 0;
        int cancellationChecks = 0;
        cancellation.shouldCancel = [&] { return ++cancellationChecks > 4; };
        ExpressionReferenceOrbitResult cancelled;
        if (!formula::buildExpressionReferenceOrbit(cancellation, cancelled) || !cancelled.cancelled || cancelled.sampleCount != 3 || cancelled.escaped || cancelled.undefined) ++failures;

        cancellation.maxIterations = 1000000;
        cancellation.memoryLimitBytes = 1;
        cancellation.shouldCancel = [] { return true; };
        ExpressionReferenceOrbitResult immediate;
        if (!formula::buildExpressionReferenceOrbit(cancellation, immediate) || !immediate.cancelled || immediate.sampleCount != 0 || immediate.samples.capacity() != 0 || immediate.tape.capacity() != 0) ++failures;

        cancellation.center.realDecimal = "malformed";
        cancellation.memoryLimitBytes = 0;
        cancellation.shouldCancel = {};
        ExpressionReferenceOrbitResult malformed;
        if (formula::buildExpressionReferenceOrbit(cancellation, malformed) || malformed.status != ExpressionReferenceBuildStatus::InputParseError || malformed.samples.capacity() != 0 || malformed.tape.capacity() != 0) ++failures;

        cancellation.center.realDecimal = "0";
        cancellation.maxIterations = 0;
        cancellation.shouldCancel = {};
        ExpressionReferenceOrbitResult zeroIterations;
        if (!formula::buildExpressionReferenceOrbit(cancellation, zeroIterations) || zeroIterations.sampleCount != 0 || zeroIterations.escaped || zeroIterations.cancelled) ++failures;
    }

    ExpressionProgram mismatchCanonical, mismatchRuntime;
    ExpressionContext mismatchFixed;
    if (compile(mismatchCanonical, "z+c") && compile(mismatchRuntime, "z-c")) {
        ExpressionProgram specializedMismatch;
        specialize(mismatchRuntime, mismatchFixed, FormulaParameter::C, specializedMismatch);
        ExpressionReferenceBuildRequest mismatch;
        mismatch.canonicalProgram = &mismatchCanonical;
        mismatch.runtimeProgram = &specializedMismatch;
        mismatch.center.realDecimal = "0";
        mismatch.maxIterations = 1;
        ExpressionReferenceOrbitResult mismatchResult;
        if (formula::buildExpressionReferenceOrbit(mismatch, mismatchResult) || mismatchResult.status != ExpressionReferenceBuildStatus::ProgramMismatch) ++failures;
    }
    ExpressionReferenceBuildRequest invalid;
    invalid.center.realDecimal = "0";
    ExpressionReferenceOrbitResult invalidResult;
    if (formula::buildExpressionReferenceOrbit(invalid, invalidResult) || invalidResult.status != ExpressionReferenceBuildStatus::InvalidRequest) ++failures;
    if (simpleRuntimeC.valid()) {
        mpf_t mixedReal, mixedImaginary;
        mpf_init_set_ui(mixedReal, 0);
        mpf_init_set_ui(mixedImaginary, 0);
        ExpressionReferenceBuildRequest mixed;
        mixed.runtimeProgram = &simpleRuntimeC;
        mixed.center.realMpf = mixedReal;
        mixed.center.imaginaryMpf = mixedImaginary;
        mixed.center.imaginaryDecimal = "0";
        mixed.maxIterations = 0;
        ExpressionReferenceOrbitResult mixedResult;
        if (formula::buildExpressionReferenceOrbit(mixed, mixedResult) || mixedResult.status != ExpressionReferenceBuildStatus::InvalidRequest) {
            printf("  mixed MPF/imaginary-decimal input accepted\n");
            ++failures;
        }
        mixed.center.realMpf = nullptr;
        mixed.center.imaginaryMpf = nullptr;
        if (formula::buildExpressionReferenceOrbit(mixed, mixedResult) || mixedResult.status != ExpressionReferenceBuildStatus::InvalidRequest) {
            printf("  imaginary-only decimal input accepted\n");
            ++failures;
        }
        mpf_clears(mixedReal, mixedImaginary, (mpf_ptr)0);

        invalid.runtimeProgram = &simpleRuntimeC;
        invalid.pixelParameter = FormulaParameter::Power;
        if (formula::buildExpressionReferenceOrbit(invalid, invalidResult)) ++failures;
        invalid.pixelParameter = FormulaParameter::C;
        invalid.precision.maximumBits = 32;
        if (formula::buildExpressionReferenceOrbit(invalid, invalidResult) || invalidResult.status != ExpressionReferenceBuildStatus::PrecisionOutOfRange) ++failures;
        invalid.precision = {};
        invalid.memoryLimitBytes = 1;
        invalid.maxIterations = 10;
        invalid.precision.requestedBits = ExpressionReferencePrecisionPolicy::ApplicationMaximumBits + 1;
        if (formula::buildExpressionReferenceOrbit(invalid, invalidResult) || invalidResult.status != ExpressionReferenceBuildStatus::ResourceLimit || invalidResult.samples.capacity() != 0 || invalidResult.tape.capacity() != 0) ++failures;

        invalid.precision.requestedBits = 32768;
        invalid.memoryLimitBytes = 4096;
        if (formula::buildExpressionReferenceOrbit(invalid, invalidResult) || invalidResult.status != ExpressionReferenceBuildStatus::ResourceLimit || invalidResult.samples.capacity() != 0 || invalidResult.tape.capacity() != 0) ++failures;
    }

    ExpressionProgram plannedProgram;
    ExpressionOrbitPlan plan;
    if (compile(plannedProgram, "z*(sin(c)+exp(z0))+abs(p0)") && plan.build(plannedProgram) && plan.profitable()) {
        ExpressionReferenceBuildRequest unsupported;
        unsupported.runtimeProgram = &plan.bodyProgram();
        unsupported.center.realDecimal = "0";
        unsupported.maxIterations = 1;
        ExpressionReferenceOrbitResult unsupportedResult;
        if (formula::buildExpressionReferenceOrbit(unsupported, unsupportedResult) || unsupportedResult.status != ExpressionReferenceBuildStatus::UnsupportedProgram) ++failures;
    } else {
        ++failures;
    }

    if (coordinateCanonical.valid() && coordinateRuntimeC.valid()) {
        ExpressionReferenceOrbitResult low, high;
        if (buildDecimal(coordinateCanonical, coordinateRuntimeC, FormulaParameter::C, coordinateFixed, deepCoordinate, "1e-500", 3, 1800, low) && buildDecimal(coordinateCanonical, coordinateRuntimeC, FormulaParameter::C, coordinateFixed, deepCoordinate, "1e-500", 3, 2400, high)) {
            if (low.precision != 1800 || high.precision != 2400 || low.samples.size() != high.samples.size()) {
                ++failures;
            } else {
                for (size_t i = 0; i < low.samples.size(); ++i) {
                    MpfrComplex lowValue(2400), highValue(2400);
                    if (!formula::reconstructMpfrFromShadows(lowValue, low.samples[i].next, low.samples[i].rootDefect) || !formula::reconstructMpfrFromShadows(highValue, high.samples[i].next, high.samples[i].rootDefect) || !closeMpfrComplex(lowValue, highValue, 90)) {
                        printf("  e500 precision convergence failed\n");
                        ++failures;
                        break;
                    }
                }
            }
        }
    }
    if (coordinateCanonical.valid() && coordinateRuntimeC.valid()) {
        ExpressionReferenceOrbitResult positiveE500;
        if (!buildDecimal(coordinateCanonical, coordinateRuntimeC, FormulaParameter::C, coordinateFixed, "1e500", "0", 1, 512, positiveE500) || !positiveE500.pixel.re.isFinite() || positiveE500.pixel.re.exponent < 1600) {
            printf("  positive e500 reference input failed\n");
            ++failures;
        }

        for (const char* outOfRange : {"1e-400000000", "1e400000000"}) {
            ExpressionReferenceBuildRequest request;
            request.canonicalProgram = &coordinateCanonical;
            request.runtimeProgram = &coordinateRuntimeC;
            request.center.realDecimal = outOfRange;
            request.maxIterations = 1000000;
            request.precision.requestedBits = 256;
            request.precision.minimumBits = 53;
            request.precision.guardBits = 0;
            request.memoryLimitBytes = 0;
            ExpressionReferenceOrbitResult rejected;
            if (formula::buildExpressionReferenceOrbit(request, rejected) || rejected.status != ExpressionReferenceBuildStatus::InputParseError || !rejected.samples.empty() || rejected.samples.capacity() != 0 || !rejected.tape.empty() || rejected.tape.capacity() != 0) {
                printf("  decimal range rejection failed [%s]\n", outOfRange);
                ++failures;
            }
        }
    }

    {
        const mpfr_exp_t runtimeExponents[] = {mpfr_get_emin(), mpfr_get_emax()};
        for (mpfr_exp_t exponent : runtimeExponents) {
            MpfrComplex boundary(128);
            MpfrComplex reconstructed(128);
            boundary.set(0.75, 0.0);
            if (mpfr_set_exp(boundary.re, exponent) != 0) {
                ++failures;
                continue;
            }
            ScaledRealShadow shadow;
            if (!formula::makeScaledRealShadow(boundary.re, shadow) || shadow.exponent != static_cast<int64_t>(exponent) || !formula::setMpfrFromScaledShadow(reconstructed.re, shadow) || !mpfr_equal_p(reconstructed.re, boundary.re)) {
                printf("  runtime exponent boundary roundtrip failed\n");
                ++failures;
            }
        }

        MpfrComplex roundedMaximum(128);
        roundedMaximum.set(1.0, 0.0);
        mpfr_nextbelow(roundedMaximum.re);
        if (mpfr_set_exp(roundedMaximum.re, mpfr_get_emax()) != 0) {
            ++failures;
        } else {
            ScaledRealShadow shadow;
            MpfrComplex reconstructed(128);
            if (!formula::makeScaledRealShadow(roundedMaximum.re, shadow) || shadow.exponent != static_cast<int64_t>(mpfr_get_emax()) || !formula::setMpfrFromScaledShadow(reconstructed.re, shadow) || mpfr_cmp(reconstructed.re, roundedMaximum.re) > 0) {
                printf("  inward emax shadow extraction failed\n");
                ++failures;
            }
        }

        MpfrComplex rejected(128);
        ScaledRealShadow outside;
        outside.mantissa = 0.75;
        if (mpfr_get_emin() > std::numeric_limits<int64_t>::min()) {
            outside.exponent = static_cast<int64_t>(mpfr_get_emin()) - 1;
            if (formula::setMpfrFromScaledShadow(rejected.re, outside)) ++failures;
        }
        if (mpfr_get_emax() < std::numeric_limits<int64_t>::max()) {
            outside.exponent = static_cast<int64_t>(mpfr_get_emax()) + 1;
            if (formula::setMpfrFromScaledShadow(rejected.re, outside)) ++failures;
        }

        if (coordinateCanonical.valid() && coordinateRuntimeC.valid()) {
            ExpressionProgram auxiliaryCanonical;
            ExpressionProgram auxiliaryRuntime;
            if (compile(auxiliaryCanonical, "z/c") && specialize(auxiliaryCanonical, coordinateFixed, FormulaParameter::C, auxiliaryRuntime)) {
                mpf_t real, imaginary;
                mpf_init2(real, 128);
                mpf_init2(imaginary, 128);
                mpf_set_ui(imaginary, 0);

                mpfr_get_f(real, roundedMaximum.re, MPFR_RNDN);
                ExpressionReferenceBuildRequest request;
                request.canonicalProgram = &auxiliaryCanonical;
                request.runtimeProgram = &auxiliaryRuntime;
                request.pixelParameter = FormulaParameter::C;
                request.center.realMpf = real;
                request.center.imaginaryMpf = imaginary;
                request.maxIterations = 1;
                request.bailout = std::numeric_limits<double>::max();
                request.precision.minimumBits = 53;
                request.precision.guardBits = 0;
                request.precision.maximumBits = 4096;
                ExpressionReferenceOrbitResult maximum;
                const ExpressionReferenceTapeNode* denominator = nullptr;
                if (formula::buildExpressionReferenceOrbit(request, maximum)) {
                    for (const auto& node : maximum.tape)
                        if (node.operation == ExpressionOracleOperation::Divide) denominator = &node;
                }
                MpfrComplex reconstructed(128);
                if (!maximum.valid || !denominator || !formula::reconstructMpfrFromShadows(reconstructed, denominator->auxiliary, denominator->auxiliaryDefect) || !closeMpfr(reconstructed.re, roundedMaximum.re, 94)) {
                    printf("  near-emax input/tape auxiliary compaction failed\n");
                    ++failures;
                }

                MpfrComplex roundedMinimum(128);
                roundedMinimum.set(1.0, 0.0);
                if (mpfr_set_exp(roundedMinimum.re, mpfr_get_emin()) != 0) {
                    ++failures;
                } else {
                    mpfr_nextabove(roundedMinimum.re);
                    mpfr_get_f(real, roundedMinimum.re, MPFR_RNDN);
                    request.canonicalProgram = &coordinateCanonical;
                    request.runtimeProgram = &coordinateRuntimeC;
                    request.maxIterations = 0;
                    ExpressionReferenceOrbitResult minimum;
                    if (formula::buildExpressionReferenceOrbit(request, minimum) || minimum.valid || minimum.status != ExpressionReferenceBuildStatus::CompactionOutOfRange || minimum.cDefect.re.isNan() || minimum.cDefect.im.isNan()) {
                        printf("  near-emin compaction failure was not explicit\n");
                        ++failures;
                    }
                }
                mpf_clears(real, imaginary, (mpf_ptr)0);
            }
        }
    }

    auto benchmarkReference = [&](const char* source, FormulaParameter pixel, const ExpressionContext& fixed, const char* real, const char* imaginary, double& seconds, double& bytesPerSample) {
        ExpressionProgram canonical, runtime;
        if (!compile(canonical, source) || !specialize(canonical, fixed, pixel, runtime)) return false;
        ExpressionReferenceBuildRequest request;
        request.canonicalProgram = &canonical;
        request.runtimeProgram = &runtime;
        request.pixelParameter = pixel;
        request.fixed = fixed;
        request.center.realDecimal = real;
        request.center.imaginaryDecimal = imaginary;
        request.maxIterations = 10000;
        request.bailout = 1e100;
        request.precision.requestedBits = 256;
        request.precision.minimumBits = 53;
        request.precision.guardBits = 0;
        auto start = Clock::now();
        ExpressionReferenceOrbitResult result;
        bool okay = formula::buildExpressionReferenceOrbit(request, result);
        seconds = since(start);
        bytesPerSample = result.sampleCount ? (double)result.memoryBytes / result.sampleCount : 0.0;
        return okay && result.sampleCount == 10000 && !result.escaped && !result.undefined;
    };
    ExpressionContext arithmeticFixed;
    ExpressionContext transcendentalFixed;
    double arithmeticSeconds = 0.0;
    double arithmeticBytes = 0.0;
    double transcendentalSeconds = 0.0;
    double transcendentalBytes = 0.0;
    if (!benchmarkReference("z*0.9999", FormulaParameter::C, arithmeticFixed, "0", "0", arithmeticSeconds, arithmeticBytes) || !benchmarkReference("sin(z)", FormulaParameter::InitialZ, transcendentalFixed, "0.1", "0.05", transcendentalSeconds, transcendentalBytes)) ++failures;

    printf("=== MPFR expression reference orbit/tape\n");
    printf("  root samples=%d companions=%d branches=%d defects=%d\n", orbitSamples, companionChecks, branchChecks, defectChecks);
    printf("  10k arithmetic: %.3f s %.1f bytes/sample; "
           "transcendental: %.3f s %.1f bytes/sample\n",
           arithmeticSeconds, arithmeticBytes, transcendentalSeconds, transcendentalBytes);
    printf("  branch clearance is point-only/not-certified; "
           "GUI dispatch remains disabled\n");
    printf("  => %s\n\n", failures == 0 ? "PASS" : "CHECK (expression reference failure)");
    return failures == 0 ? 0 : 1;
}

static int runExpressionScaledCase() {
    using formula::Complex;
    using formula::ExpressionContext;
    using formula::ExpressionError;
    using formula::ExpressionOracle;
    using formula::ExpressionOracleContext;
    using formula::ExpressionOracleOperation;
    using formula::ExpressionProgram;
    using formula::ExpressionReferenceBuildRequest;
    using formula::ExpressionReferenceOrbitResult;
    using formula::ExpressionScaledResidualEvaluator;
    using formula::ExpressionScaledResidualInput;
    using formula::ExpressionScaledResidualResult;
    using formula::ExpressionScaledResidualStatus;
    using formula::MpfrComplex;
    using formula::ScaledArithmeticStatus;
    using formula::ScaledComplexValue;
    using formula::ScaledRealBall;
    using formula::ScaledRealValue;

    int failures = 0;
    size_t framePixels = 0;
    size_t fallbackPixels = 0;
    size_t seriesOperations = 0;
    double maximumResidualRelativeError = 0.0;
    std::array<bool, static_cast<size_t>(ExpressionOracleOperation::Polar) + 1> opcodeCoverage{};

    auto elapsed = [](Clock::time_point start) { return std::chrono::duration<double>(Clock::now() - start).count(); };
    auto compile = [&](ExpressionProgram& program, const char* source) {
        ExpressionError error;
        if (!program.compile(source, &error)) {
            printf("  scaled compile failed @ %zu: %s [%s]\n", error.position, error.message.c_str(), source);
            ++failures;
            return false;
        }
        return true;
    };
    auto markCoverage = [&](const ExpressionReferenceOrbitResult& reference) {
        if (reference.samples.empty()) return;
        const auto& sample = reference.samples.front();
        for (size_t node = 0; node < sample.tapeCount; ++node) {
            const auto operation = reference.tape[static_cast<size_t>(sample.tapeOffset) + node].operation;
            opcodeCoverage[static_cast<size_t>(operation)] = true;
        }
    };
    auto sameScaled = [](const ScaledRealValue& left, const ScaledRealValue& right) {
        uint64_t leftBits = 0, rightBits = 0;
        std::memcpy(&leftBits, &left.mantissa, sizeof(leftBits));
        std::memcpy(&rightBits, &right.mantissa, sizeof(rightBits));
        return leftBits == rightBits && left.exponent == right.exponent;
    };
    auto sameScaledComplex = [&](const ScaledComplexValue& left, const ScaledComplexValue& right) { return sameScaled(left.re, right.re) && sameScaled(left.im, right.im); };
    {
        struct DiffAbsCase {
            const char* name;
            ScaledRealBall reference;
            ScaledRealBall residual;
            ExpressionScaledResidualStatus expected;
            bool expectPositiveZero = false;
        };
        auto value = [](double mantissa, int64_t exponent) { return ScaledRealValue{mantissa, exponent}; };
        const DiffAbsCase cases[] = {{"same-positive", {value(0.75, 3), {}}, {value(0.5, 0), {}}, ExpressionScaledResidualStatus::Success},
                                     {"same-negative", {value(-0.75, 3), {}}, {value(0.5, 0), {}}, ExpressionScaledResidualStatus::Success},
                                     {"same-positive-radii", {value(0.75, 3), value(0.5, -2)}, {value(0.5, 0), value(0.5, -3)}, ExpressionScaledResidualStatus::Success},
                                     {"cross-positive-negative", {value(0.5, 2), {}}, {value(-0.75, 3), {}}, ExpressionScaledResidualStatus::Success},
                                     {"cross-negative-positive", {value(-0.5, 2), {}}, {value(0.75, 3), {}}, ExpressionScaledResidualStatus::Success},
                                     {"cross-positive-negative-radii", {value(0.5, 2), value(0.5, -3)}, {value(-0.75, 3), value(0.5, -2)}, ExpressionScaledResidualStatus::Success},
                                     {"subnormal-delta", {value(0.75, 0), {}}, {value(0.5, -1073), {}}, ExpressionScaledResidualStatus::Success},
                                     {"e500-delta", {value(0.75, 0), {}}, {value(-0.5, -1660), {}}, ExpressionScaledResidualStatus::Success},
                                     {"negative-positive-zero-delta", {value(-0.75, 0), {}}, {value(0.0, 0), {}}, ExpressionScaledResidualStatus::Success, true},
                                     {"negative-negative-zero-delta", {value(-0.75, 0), {}}, {value(-0.0, 0), {}}, ExpressionScaledResidualStatus::Success, true},
                                     {"exact-zero-input", {value(-0.0, 0), {}}, {value(0.0, 0), {}}, ExpressionScaledResidualStatus::Success, true},
                                     {"reference-positive-zero", {value(0.0, 0), {}}, {value(0.5, 0), {}}, ExpressionScaledResidualStatus::BranchUncertain},
                                     {"reference-negative-zero", {value(-0.0, 0), {}}, {value(-0.5, 0), {}}, ExpressionScaledResidualStatus::BranchUncertain},
                                     {"pixel-positive-zero", {value(0.5, 1), {}}, {value(-0.5, 1), {}}, ExpressionScaledResidualStatus::BranchUncertain},
                                     {"touching-zero", {value(0.5, 1), value(0.5, 1)}, {value(0.5, 0), {}}, ExpressionScaledResidualStatus::BranchUncertain}};
        mpfr_t reference, residual, endpoint, exact;
        mpfr_t referenceRadius, residualRadius;
        mpfr_t actualReference, actualResidual;
        mpfr_t midpoint, radius, difference;
        mpfr_inits2(4096, reference, residual, endpoint, exact, referenceRadius, residualRadius, actualReference, actualResidual, midpoint, radius, difference, (mpfr_ptr)0);
        for (const DiffAbsCase& test : cases) {
            ScaledRealBall output;
            const ExpressionScaledResidualStatus status = formula::certifiedScaledDiffAbsReal(test.reference, test.residual, output);
            bool okay = status == test.expected;
            if (okay && status == ExpressionScaledResidualStatus::Success) {
                okay = formula::setMpfrFromScaledValue(reference, test.reference.value) && formula::setMpfrFromScaledValue(residual, test.residual.value) && formula::setMpfrFromScaledValue(referenceRadius, test.reference.radius) && formula::setMpfrFromScaledValue(residualRadius, test.residual.radius) && formula::setMpfrFromScaledValue(midpoint, output.value) && formula::setMpfrFromScaledValue(radius, output.radius);
                for (int referenceSide = -1; okay && referenceSide <= 1; ++referenceSide) {
                    mpfr_set(actualReference, reference, MPFR_RNDN);
                    if (referenceSide < 0)
                        mpfr_sub(actualReference, actualReference, referenceRadius, MPFR_RNDN);
                    else if (referenceSide > 0)
                        mpfr_add(actualReference, actualReference, referenceRadius, MPFR_RNDN);
                    for (int residualSide = -1; okay && residualSide <= 1; ++residualSide) {
                        mpfr_set(actualResidual, residual, MPFR_RNDN);
                        if (residualSide < 0)
                            mpfr_sub(actualResidual, actualResidual, residualRadius, MPFR_RNDN);
                        else if (residualSide > 0)
                            mpfr_add(actualResidual, actualResidual, residualRadius, MPFR_RNDN);
                        mpfr_add(endpoint, actualReference, actualResidual, MPFR_RNDN);
                        mpfr_abs(endpoint, endpoint, MPFR_RNDN);
                        mpfr_abs(exact, actualReference, MPFR_RNDN);
                        mpfr_sub(exact, endpoint, exact, MPFR_RNDN);
                        mpfr_sub(difference, exact, midpoint, MPFR_RNDU);
                        mpfr_abs(difference, difference, MPFR_RNDU);
                        okay = mpfr_cmp(difference, radius) <= 0;
                    }
                }
                if (test.expectPositiveZero) okay = okay && output.value.isZero() && !std::signbit(output.value.mantissa);
            }
            if (!okay) {
                printf("  certified diffabs failed [%s] status=%s\n", test.name, formula::expressionScaledResidualStatusName(status));
                ++failures;
            }
        }
        mpfr_clears(reference, residual, endpoint, exact, referenceRadius, residualRadius, actualReference, actualResidual, midpoint, radius, difference, (mpfr_ptr)0);
    }
    auto configureOracle = [](ExpressionOracleContext& context, const ExpressionContext& fixed, FormulaParameter pixel, const std::string& centerReal, const std::string& centerImaginary) {
        context.c.set(fixed.c.real(), fixed.c.imag());
        context.z0.set(fixed.z0.real(), fixed.z0.imag());
        for (size_t parameter = 0; parameter < fixed.parameters.size(); ++parameter) { context.parameters[parameter].set(fixed.parameters[parameter].real(), fixed.parameters[parameter].imag()); }
        MpfrComplex center(context.z.precision());
        if (!center.set(centerReal, centerImaginary)) return false;
        if (pixel == FormulaParameter::C)
            context.c.set(center);
        else
            context.z0.set(center);
        context.z.set(context.z0);
        return true;
    };
    auto outsideScaled = [](const ScaledComplexValue& value, double bailout, bool& outside) {
        ScaledRealValue norm, threshold;
        ScaledArithmeticStatus status = formula::scaledNormSquared(value, norm);
        if (status != ScaledArithmeticStatus::Success) return status;
        status = formula::makeScaledRealValue(bailout * bailout, threshold);
        if (status != ScaledArithmeticStatus::Success) return status;
        outside = formula::compareScaledNonnegative(norm, threshold) > 0;
        return ScaledArithmeticStatus::Success;
    };

    // Arithmetic value contract, including signed zero, MPFR round-trip,
    // unrepresentable binary64 values, and explicit exponent overflow.
    {
        ScaledRealValue negativeZero, positiveZero;
        ScaledRealValue tiny, product, overflow;
        ScaledRealValue boundary, boundaryResult;
        MpfrComplex mpfrValue(2048), reconstructed(2048);
        mpfrValue.set("1e-500", "-1e500");
        ScaledComplexValue compact;
        double converted = 0.0;
        if (formula::makeScaledRealValue(-0.0, negativeZero) != ScaledArithmeticStatus::Success || formula::scaledNegate(negativeZero, positiveZero) != ScaledArithmeticStatus::Success || !std::signbit(negativeZero.mantissa) || std::signbit(positiveZero.mantissa) || formula::makeScaledComplexValue(mpfrValue, compact) != ScaledArithmeticStatus::Success || !formula::setMpfrFromScaledValue(reconstructed, compact) || mpfr_cmp(reconstructed.re, mpfrValue.re) == 0 || formula::scaledValueToDouble(compact.im, converted) || formula::scaledValueToDouble(compact.re, converted)) {
            printf("  scaled value conversion contract failed\n");
            ++failures;
        }
        tiny.mantissa = 0.5;
        tiny.exponent = std::numeric_limits<int64_t>::max();
        if (formula::makeScaledRealValue(2.0, product) != ScaledArithmeticStatus::Success || formula::scaledMultiply(tiny, product, overflow) != ScaledArithmeticStatus::ExponentRange) {
            printf("  scaled exponent overflow was not explicit\n");
            ++failures;
        }
        boundary.mantissa = 0.5;
        boundary.exponent = std::numeric_limits<int64_t>::max();
        if (formula::makeScaledRealValue(1.0, product) != ScaledArithmeticStatus::Success || formula::scaledMultiply(boundary, product, boundaryResult) != ScaledArithmeticStatus::Success || boundaryResult.mantissa != 0.5 || boundaryResult.exponent != boundary.exponent) {
            printf("  scaled maximum exponent renormalization failed\n");
            ++failures;
        }
        boundary.exponent = std::numeric_limits<int64_t>::min();
        if (formula::scaledDivideByDouble(boundary, 1.0, boundaryResult) != ScaledArithmeticStatus::Success || boundaryResult.mantissa != 0.5 || boundaryResult.exponent != boundary.exponent) {
            printf("  scaled minimum exponent renormalization failed\n");
            ++failures;
        }
        if (formula::makeScaledRealValue(2.0, product) != ScaledArithmeticStatus::Success || formula::scaledMultiply(boundary, product, boundaryResult) != ScaledArithmeticStatus::Success || boundaryResult.exponent != std::numeric_limits<int64_t>::min() + 1 || formula::scaledMultiply(product, boundary, boundaryResult) != ScaledArithmeticStatus::Success || boundaryResult.exponent != std::numeric_limits<int64_t>::min() + 1) {
            printf("  scaled minimum exponent product ordering failed\n");
            ++failures;
        }
        boundary.exponent = std::numeric_limits<int64_t>::max();
        if (formula::scaledDivideByDouble(boundary, 1.0, boundaryResult) != ScaledArithmeticStatus::Success || boundaryResult.mantissa != 0.5 || boundaryResult.exponent != boundary.exponent) {
            printf("  scaled maximum exponent division failed\n");
            ++failures;
        }
        ScaledRealValue subnormal;
        subnormal.mantissa = 0.75;
        subnormal.exponent = -1074;
        if (!formula::scaledValueToDouble(subnormal, converted) || converted != std::numeric_limits<double>::denorm_min()) {
            printf("  scaled subnormal conversion failed\n");
            ++failures;
        }
    }

    // Moderate local deltas force adaptive nonlinear terms rather than the
    // e500 linear limit. Compare the returned residual with an independent
    // higher-precision MPFR subtraction.
    struct LocalFormula {
        const char* source;
        bool series;
    };
    const LocalFormula localFormulas[] = {{"sin(z)", true},
                                          {"cos(z)", true},
                                          {"sinh(z)", true},
                                          {"cosh(z)", true},
                                          {"exp(z)", true},
                                          {"complex(real(conj(z)),imag(z))"
                                           "+complex(norm(z),0)-(-z)",
                                           false}};
    for (const LocalFormula& formulaCase : localFormulas) {
        const char* source = formulaCase.source;
        ExpressionProgram canonical, runtime;
        ExpressionContext fixed;
        if (!compile(canonical, source)) continue;
        fixed.c = {0.0, 0.0};
        ExpressionError error;
        if (!canonical.specialize(fixed, FormulaParameter::InitialZ, runtime, &error)) {
            ++failures;
            continue;
        }
        ExpressionReferenceBuildRequest request;
        request.canonicalProgram = &canonical;
        request.runtimeProgram = &runtime;
        request.pixelParameter = FormulaParameter::InitialZ;
        request.fixed = fixed;
        request.center.realDecimal = "0.2";
        request.center.imaginaryDecimal = "0.1";
        request.bailout = 100.0;
        request.maxIterations = 1;
        request.precision.requestedBits = 384;
        request.precision.minimumBits = 53;
        request.precision.guardBits = 0;
        ExpressionReferenceOrbitResult reference;
        ExpressionScaledResidualEvaluator evaluator;
        MpfrComplex delta(512);
        delta.set("0.03125", "-0.015625");
        ScaledComplexValue scaledDelta;
        if (!formula::buildExpressionReferenceOrbit(request, reference) || !evaluator.prepare(runtime, reference) || formula::makeScaledComplexValue(delta, scaledDelta) != ScaledArithmeticStatus::Success) {
            printf("  nonlinear scaled setup failed [%s]\n", source);
            ++failures;
            continue;
        }
        markCoverage(reference);
        ExpressionScaledResidualInput input;
        input.z = scaledDelta;
        input.z0 = scaledDelta;
        input.iteration = 0;
        ExpressionScaledResidualResult evaluated = evaluator.evaluate(0, input);

        ExpressionOracleContext base(512), actual(512);
        configureOracle(base, fixed, FormulaParameter::InitialZ, "0.2", "0.1");
        configureOracle(actual, fixed, FormulaParameter::InitialZ, "0.2", "0.1");
        mpfr_add(actual.z.re, actual.z.re, delta.re, MPFR_RNDN);
        mpfr_add(actual.z.im, actual.z.im, delta.im, MPFR_RNDN);
        actual.z0.set(actual.z);
        MpfrComplex baseOutput(512), actualOutput(512);
        std::string baseError, actualError;
        MpfrComplex expected(512), reconstructed(512);
        mpfr_t magnitude, difference, ratio;
        mpfr_inits2(512, magnitude, difference, ratio, (mpfr_ptr)0);
        bool okay = ExpressionOracle::evaluate(runtime, base, baseOutput, &baseError) && ExpressionOracle::evaluate(runtime, actual, actualOutput, &actualError);
        if (okay) {
            mpfr_sub(expected.re, actualOutput.re, baseOutput.re, MPFR_RNDN);
            mpfr_sub(expected.im, actualOutput.im, baseOutput.im, MPFR_RNDN);
            okay = evaluated.status == ExpressionScaledResidualStatus::Success && evaluated.uncertified == formulaCase.series && (formulaCase.series ? !evaluated.remainderEstimate.isZero() : evaluated.remainderEstimate.isZero()) && formula::setMpfrFromScaledValue(reconstructed, evaluated.residual);
        }
        if (okay) {
            mpfr_sub(reconstructed.re, reconstructed.re, expected.re, MPFR_RNDN);
            mpfr_sub(reconstructed.im, reconstructed.im, expected.im, MPFR_RNDN);
            mpfr_hypot(difference, reconstructed.re, reconstructed.im, MPFR_RNDN);
            mpfr_hypot(magnitude, expected.re, expected.im, MPFR_RNDN);
            if (mpfr_zero_p(magnitude))
                okay = mpfr_zero_p(difference);
            else {
                mpfr_div(ratio, difference, magnitude, MPFR_RNDU);
                okay = mpfr_cmp_d(ratio, 3e-14) <= 0;
            }
        }
        mpfr_clears(magnitude, difference, ratio, (mpfr_ptr)0);
        if (!okay) {
            printf("  nonlinear residual mismatch [%s] status=%s\n", source, formula::expressionScaledResidualStatusName(evaluated.status));
            ++failures;
        }
    }

    // Parameter arguments must survive the trace so p0 and p1 leaves cannot
    // alias. This also covers cancellation and exact primary/defect rebasing.
    {
        ExpressionProgram program;
        ExpressionContext fixed;
        fixed.parameters[0] = {0.25, -0.125};
        fixed.parameters[1] = {-0.5, 0.375};
        if (compile(program, "(z+p0)-p0+p1-p1")) {
            ExpressionReferenceBuildRequest request;
            request.runtimeProgram = &program;
            request.pixelParameter = FormulaParameter::C;
            request.fixed = fixed;
            request.center.realDecimal = "0.123456789012345678901234567890123456789";
            request.center.imaginaryDecimal = "1e-500";
            request.bailout = 100.0;
            request.maxIterations = 2;
            request.precision.requestedBits = 1800;
            request.precision.minimumBits = 53;
            request.precision.guardBits = 0;
            ExpressionReferenceOrbitResult reference;
            ExpressionScaledResidualEvaluator evaluator;
            MpfrComplex delta(1800);
            delta.set("1e-500", "-2e-500");
            ScaledComplexValue scaledDelta;
            bool sawP0 = false, sawP1 = false;
            if (formula::buildExpressionReferenceOrbit(request, reference)) {
                markCoverage(reference);
                for (const auto& node : reference.tape) {
                    if (node.operation != ExpressionOracleOperation::Parameter) continue;
                    sawP0 = sawP0 || node.argument == 0;
                    sawP1 = sawP1 || node.argument == 1;
                }
            }
            if (!reference.valid || !evaluator.prepare(program, reference) || formula::makeScaledComplexValue(delta, scaledDelta) != ScaledArithmeticStatus::Success || !(sawP0 && sawP1)) {
                printf("  parameter trace argument setup failed\n");
                ++failures;
            } else {
                ExpressionScaledResidualInput input;
                input.z = scaledDelta;
                input.iteration = 0;
                ExpressionScaledResidualResult evaluated = evaluator.evaluate(0, input);
                formula::ExpressionScaledPrimaryRelativeState primaryRelative;
                ScaledComplexValue nextExact;
                bool rebaseOkay = evaluated.status == ExpressionScaledResidualStatus::Success && sameScaledComplex(evaluated.residual, scaledDelta) && formula::makeExpressionResidualNextPrimaryState(reference.samples[0], evaluated.residual, primaryRelative) == ScaledArithmeticStatus::Success && formula::resetExpressionResidualToExactSample(reference.samples[1], primaryRelative, nextExact) == ScaledArithmeticStatus::Success && sameScaledComplex(nextExact, scaledDelta);
                if (!rebaseOkay) {
                    printf("  exact/primary delta convention failed\n");
                    ++failures;
                }

                ExpressionReferenceOrbitResult malformed = reference;
                for (auto& node : malformed.tape) {
                    if (node.operation == ExpressionOracleOperation::Parameter) {
                        node.argument = 7;
                        break;
                    }
                }
                ExpressionScaledResidualEvaluator rejected;
                if (rejected.prepare(program, malformed)) {
                    printf("  malformed parameter argument accepted\n");
                    ++failures;
                }
                malformed = reference;
                malformed.tape[malformed.samples[0].tapeOffset + malformed.samples[0].rootNode].leftNode = UINT16_MAX;
                if (rejected.prepare(program, malformed)) {
                    printf("  malformed child layout accepted\n");
                    ++failures;
                }
                malformed = reference;
                ++malformed.programSemanticHash;
                if (rejected.prepare(program, malformed)) {
                    printf("  semantic mismatch accepted\n");
                    ++failures;
                }
                if (evaluator.evaluate(reference.samples.size(), input).status != ExpressionScaledResidualStatus::InvalidTape || ([&] {
                                                                                                                                     input.iteration = 1;
                                                                                                                                     return evaluator.evaluate(0, input).status;
                                                                                                                                 })() != ExpressionScaledResidualStatus::InvalidInput) {
                    printf("  sample range/iteration status failed\n");
                    ++failures;
                }
            }
        }
    }

    {
        ExpressionProgram canonical, runtime;
        ExpressionContext fixed;
        fixed.z0 = {0.125, -0.0625};
        if (compile(canonical, "z+c")) {
            ExpressionError error;
            canonical.specialize(fixed, FormulaParameter::C, runtime, &error);
            ExpressionReferenceBuildRequest request;
            request.canonicalProgram = &canonical;
            request.runtimeProgram = &runtime;
            request.pixelParameter = FormulaParameter::C;
            request.fixed = fixed;
            request.center.realDecimal = "0.123456789012345678901234567890123456789";
            request.center.imaginaryDecimal = "1e-500";
            request.bailout = 100.0;
            request.maxIterations = 2;
            request.precision.requestedBits = 1800;
            request.precision.minimumBits = 53;
            request.precision.guardBits = 0;
            ExpressionReferenceOrbitResult reference;
            ExpressionScaledResidualEvaluator evaluator;
            MpfrComplex delta(1800);
            delta.set("-3e-500", "2e-500");
            ScaledComplexValue scaledDelta;
            if (!formula::buildExpressionReferenceOrbit(request, reference) || !evaluator.prepare(runtime, reference) || formula::makeScaledComplexValue(delta, scaledDelta) != ScaledArithmeticStatus::Success || (reference.samples[0].rootDefect.re.isZero() && reference.samples[0].rootDefect.im.isZero())) {
                printf("  reference-defect scaled setup failed\n");
                ++failures;
            } else {
                ExpressionScaledResidualInput input;
                input.c = scaledDelta;
                input.iteration = 0;
                const ExpressionScaledResidualResult evaluated = evaluator.evaluate(0, input);
                formula::ExpressionScaledPrimaryRelativeState primaryRelative;
                ScaledComplexValue nextExact;
                if (evaluated.status != ExpressionScaledResidualStatus::Success || !sameScaledComplex(evaluated.residual, scaledDelta) || formula::makeExpressionResidualNextPrimaryState(reference.samples[0], evaluated.residual, primaryRelative) != ScaledArithmeticStatus::Success || formula::resetExpressionResidualToExactSample(reference.samples[1], primaryRelative, nextExact) != ScaledArithmeticStatus::Success || !sameScaledComplex(nextExact, scaledDelta)) {
                    printf("  reference-defect propagation failed\n");
                    ++failures;
                }
            }
        }
    }

    auto checkExplicitStatus = [&](const char* source, FormulaParameter pixel, const char* centerReal, const char* centerImaginary, const ExpressionContext& fixed, ExpressionScaledResidualStatus expected, const ScaledComplexValue& delta, bool corruptCompanion = false, bool removeCompanion = false) {
        ExpressionProgram canonical, runtime;
        if (!compile(canonical, source)) return;
        ExpressionError error;
        if (!canonical.specialize(fixed, pixel, runtime, &error)) {
            ++failures;
            return;
        }
        ExpressionReferenceBuildRequest request;
        request.canonicalProgram = &canonical;
        request.runtimeProgram = &runtime;
        request.pixelParameter = pixel;
        request.fixed = fixed;
        request.center.realDecimal = centerReal;
        request.center.imaginaryDecimal = centerImaginary;
        request.bailout = 1e100;
        request.maxIterations = 1;
        request.precision.requestedBits = 512;
        request.precision.minimumBits = 53;
        request.precision.guardBits = 0;
        ExpressionReferenceOrbitResult reference;
        if (!formula::buildExpressionReferenceOrbit(request, reference)) {
            ++failures;
            return;
        }
        if (corruptCompanion || removeCompanion) {
            for (auto& node : reference.tape) {
                if (node.operation != ExpressionOracleOperation::Sin) continue;
                if (removeCompanion) node.flags &= static_cast<uint16_t>(~formula::OracleTraceHasCompanion);
                if (corruptCompanion) node.auxiliary.re.mantissa = std::numeric_limits<double>::quiet_NaN();
                break;
            }
        }
        ExpressionScaledResidualEvaluator evaluator;
        const bool prepared = evaluator.prepare(runtime, reference);
        if (corruptCompanion || removeCompanion) {
            if (prepared || evaluator.preparationStatus() != ExpressionScaledResidualStatus::InvalidTape) {
                printf("  malformed companion accepted [%s]\n", source);
                ++failures;
            }
            return;
        }
        ExpressionScaledResidualInput input;
        input.iteration = 0;
        if (pixel == FormulaParameter::InitialZ) {
            input.z = delta;
            input.z0 = delta;
        } else {
            input.c = delta;
        }
        ExpressionScaledResidualResult evaluated = prepared ? evaluator.evaluate(0, input) : ExpressionScaledResidualResult{};
        if (!prepared || evaluated.status != expected) {
            printf("  explicit status [%s]: got %s expected %s\n", source, formula::expressionScaledResidualStatusName(evaluated.status), formula::expressionScaledResidualStatusName(expected));
            ++failures;
        }
    };
    {
        ScaledComplexValue tinyDelta;
        formula::makeScaledComplexValue(Complex{1e-20, -2e-20}, tinyDelta);
        ExpressionContext fixed;
        fixed.c = {1.0, 0.0};
        checkExplicitStatus("z/c", FormulaParameter::InitialZ, "0.25", "0", fixed, ExpressionScaledResidualStatus::Unsupported, tinyDelta);
        checkExplicitStatus("log(z)", FormulaParameter::InitialZ, "0.25", "0", fixed, ExpressionScaledResidualStatus::BranchUncertain, tinyDelta);
        fixed.c = {0.0, 0.0};
        checkExplicitStatus("z/(c-c)", FormulaParameter::InitialZ, "0.25", "0", fixed, ExpressionScaledResidualStatus::Singular, tinyDelta);
        checkExplicitStatus("sin(z)", FormulaParameter::InitialZ, "0.25", "0", fixed, ExpressionScaledResidualStatus::Success, tinyDelta, true, false);
        checkExplicitStatus("sin(z)", FormulaParameter::InitialZ, "0.25", "0", fixed, ExpressionScaledResidualStatus::Success, tinyDelta, false, true);

        ScaledComplexValue outsideSeries;
        formula::makeScaledComplexValue(Complex{0.125, 0.0}, outsideSeries);
        checkExplicitStatus("sin(z)", FormulaParameter::InitialZ, "0.25", "0", fixed, ExpressionScaledResidualStatus::Unsupported, outsideSeries);
    }

    // A nonfinite stored output is a runtime status, while a finite but
    // impossible arithmetic exponent reports exponent-range.
    {
        ExpressionProgram canonical, runtime;
        ExpressionContext fixed;
        if (compile(canonical, "z*z")) {
            ExpressionError error;
            canonical.specialize(fixed, FormulaParameter::InitialZ, runtime, &error);
            ExpressionReferenceBuildRequest request;
            request.canonicalProgram = &canonical;
            request.runtimeProgram = &runtime;
            request.pixelParameter = FormulaParameter::InitialZ;
            request.center.realDecimal = "0.5";
            request.center.imaginaryDecimal = "0";
            request.bailout = 1e100;
            request.maxIterations = 1;
            request.precision.requestedBits = 256;
            request.precision.minimumBits = 53;
            request.precision.guardBits = 0;
            ExpressionReferenceOrbitResult reference;
            if (formula::buildExpressionReferenceOrbit(request, reference)) {
                ExpressionReferenceOrbitResult nonfinite = reference;
                nonfinite.tape[0].output.re.mantissa = std::numeric_limits<double>::infinity();
                ExpressionScaledResidualEvaluator evaluator;
                ExpressionScaledResidualInput input;
                input.iteration = 0;
                if (!evaluator.prepare(runtime, nonfinite) || evaluator.evaluate(0, input).status != ExpressionScaledResidualStatus::Nonfinite) {
                    printf("  nonfinite tape status failed\n");
                    ++failures;
                }
                if (!evaluator.prepare(runtime, reference)) {
                    ++failures;
                } else {
                    input.z.re.mantissa = 0.5;
                    input.z.re.exponent = std::numeric_limits<int64_t>::max();
                    input.z0 = input.z;
                    if (evaluator.evaluate(0, input).status != ExpressionScaledResidualStatus::ExponentRange) {
                        printf("  evaluator exponent-range status failed\n");
                        ++failures;
                    }
                }
            } else {
                ++failures;
            }
        }
    }

    struct FrameSpec {
        const char* name;
        const char* source;
        FormulaParameter pixel;
        const char* centerReal;
        const char* centerImaginary;
        ExpressionContext fixed;
        double bailout;
        int maxIterations;
        int width;
        int height;
        bool preserveParameters = false;
        bool benchmark = false;
    };
    std::vector<FrameSpec> frames;
    FrameSpec arithmetic{"arithmetic-c", "z*z+c", FormulaParameter::C, "-2", "0", {}, 4.0, 850, 64, 44, false, true};
    frames.push_back(arithmetic);
    FrameSpec sine{"sine-c", "sin(z)+c", FormulaParameter::C, "0", "0", {}, 4.0, 40, 13, 9};
    frames.push_back(sine);
    FrameSpec exponential{"exponential-c", "exp(0.1*z)+c", FormulaParameter::C, "-1", "0", {}, 4.0, 50, 13, 9};
    frames.push_back(exponential);
    FrameSpec parameter{"parameter-invariant-c", "z*z+c+sin(p0)+exp(p1)-1+0*n", FormulaParameter::C, "-2", "0", {}, 4.0, 850, 9, 7, true};
    parameter.fixed.parameters[0] = {};
    parameter.fixed.parameters[1] = {};
    frames.push_back(parameter);
    FrameSpec initialZ{"arithmetic-z0", "2*z", FormulaParameter::InitialZ, "0", "0", {}, 4.0, 1670, 9, 7};
    initialZ.fixed.c = {};
    frames.push_back(initialZ);

    double benchmarkReferenceSeconds = 0.0;
    double benchmarkScaledSeconds = 0.0;
    double benchmarkMpfrSeconds = 0.0;
    double benchmarkMinimumSpeedup = 0.0;
    size_t benchmarkMemory = 0;
    size_t benchmarkBytesPerSample = 0;

    for (const FrameSpec& frame : frames) {
        ExpressionProgram canonical, runtime;
        if (!compile(canonical, frame.source)) continue;
        if (frame.preserveParameters) {
            runtime = canonical;
        } else {
            ExpressionError error;
            if (!canonical.specialize(frame.fixed, frame.pixel, runtime, &error)) {
                printf("  scaled specialize failed: %s [%s]\n", error.message.c_str(), frame.source);
                ++failures;
                continue;
            }
        }

        ExpressionReferenceBuildRequest request;
        request.canonicalProgram = frame.preserveParameters ? nullptr : &canonical;
        request.runtimeProgram = &runtime;
        request.pixelParameter = frame.pixel;
        request.fixed = frame.fixed;
        request.center.realDecimal = frame.centerReal;
        request.center.imaginaryDecimal = frame.centerImaginary;
        request.bailout = frame.bailout;
        request.maxIterations = frame.maxIterations;
        request.precision.viewBits = 1750;
        request.precision.minimumBits = 53;
        request.precision.guardBits = 64;
        request.precision.maximumBits = 4096;
        const Clock::time_point referenceStart = Clock::now();
        ExpressionReferenceOrbitResult reference;
        const bool built = formula::buildExpressionReferenceOrbit(request, reference);
        const double referenceSeconds = elapsed(referenceStart);
        ExpressionScaledResidualEvaluator evaluator;
        if (!built || reference.samples.size() != static_cast<size_t>(frame.maxIterations) || reference.escaped || reference.undefined || !evaluator.prepare(runtime, reference)) {
            printf("  frame reference failed [%s]: %s / %s\n", frame.name, reference.error.c_str(), evaluator.error().c_str());
            ++failures;
            continue;
        }
        markCoverage(reference);

        const size_t pixelCount = static_cast<size_t>(frame.width) * frame.height;
        std::vector<ScaledComplexValue> offsets(pixelCount);
        std::set<std::tuple<int64_t, uint64_t, int64_t, uint64_t>> scaledCoordinates;
        std::set<std::pair<double, double>> doubleCoordinates;
        MpfrComplex center(reference.precision);
        MpfrComplex offset(reference.precision);
        MpfrComplex coordinate(reference.precision);
        mpfr_t zoom;
        mpfr_init2(zoom, reference.precision);
        bool coordinateOkay = center.set(frame.centerReal, frame.centerImaginary) && mpfr_set_str(zoom, "1e500", 10, MPFR_RNDN) == 0;
        for (int y = 0; coordinateOkay && y < frame.height; ++y) {
            for (int x = 0; x < frame.width; ++x) {
                const size_t index = static_cast<size_t>(y) * frame.width + x;
                mpfr_set_si(offset.re, x - frame.width / 2, MPFR_RNDN);
                mpfr_div(offset.re, offset.re, zoom, MPFR_RNDN);
                mpfr_set_si(offset.im, y - frame.height / 2, MPFR_RNDN);
                mpfr_div(offset.im, offset.im, zoom, MPFR_RNDN);
                if (formula::makeScaledComplexValue(offset, offsets[index]) != ScaledArithmeticStatus::Success) {
                    coordinateOkay = false;
                    break;
                }
                mpfr_add(coordinate.re, center.re, offset.re, MPFR_RNDN);
                mpfr_add(coordinate.im, center.im, offset.im, MPFR_RNDN);
                const double directReal = mpfr_get_d(coordinate.re, MPFR_RNDN);
                const double directImaginary = mpfr_get_d(coordinate.im, MPFR_RNDN);
                uint64_t scaledRealBits = 0;
                uint64_t scaledImaginaryBits = 0;
                std::memcpy(&scaledRealBits, &offsets[index].re.mantissa, sizeof(scaledRealBits));
                std::memcpy(&scaledImaginaryBits, &offsets[index].im.mantissa, sizeof(scaledImaginaryBits));
                doubleCoordinates.emplace(directReal, directImaginary);
                scaledCoordinates.emplace(offsets[index].re.exponent, scaledRealBits, offsets[index].im.exponent, scaledImaginaryBits);
            }
        }
        mpfr_clear(zoom);
        if (!coordinateOkay || doubleCoordinates.size() != 1 || scaledCoordinates.size() != pixelCount) {
            printf("  e500 coordinate distinction failed [%s] double/scaled=%zu/%zu\n", frame.name, doubleCoordinates.size(), scaledCoordinates.size());
            ++failures;
            continue;
        }

        struct PixelState {
            int escapeIteration = 0;
            int stateIteration = 0;
            ExpressionScaledResidualStatus status = ExpressionScaledResidualStatus::Success;
            ScaledComplexValue residual;
        };
        auto renderScaled = [&](std::vector<PixelState>& states, size_t& localSeriesOperations) {
            states.assign(pixelCount, {});
            localSeriesOperations = 0;
            for (size_t index = 0; index < pixelCount; ++index) {
                PixelState& pixel = states[index];
                pixel.escapeIteration = frame.maxIterations;
                ScaledComplexValue stateDelta;
                ExpressionScaledResidualInput input;
                if (frame.pixel == FormulaParameter::C) {
                    input.c = offsets[index];
                } else {
                    input.z0 = offsets[index];
                    stateDelta = offsets[index];
                }
                ScaledComplexValue initialBase, initialValue;
                ScaledArithmeticStatus arithmeticStatus = formula::makeScaledComplexValue(reference.initialZ, reference.initialZDefect, initialBase);
                if (arithmeticStatus == ScaledArithmeticStatus::Success) arithmeticStatus = formula::scaledAdd(initialBase, stateDelta, initialValue);
                bool escaped = false;
                if (arithmeticStatus == ScaledArithmeticStatus::Success) arithmeticStatus = outsideScaled(initialValue, frame.bailout, escaped);
                if (arithmeticStatus != ScaledArithmeticStatus::Success) {
                    pixel.status = arithmeticStatus == ScaledArithmeticStatus::ExponentRange ? ExpressionScaledResidualStatus::ExponentRange : ExpressionScaledResidualStatus::Nonfinite;
                    continue;
                }
                if (escaped) {
                    pixel.escapeIteration = 0;
                    pixel.stateIteration = 0;
                    pixel.residual = stateDelta;
                    continue;
                }
                for (int iteration = 0; iteration < frame.maxIterations; ++iteration) {
                    input.z = stateDelta;
                    input.iteration = iteration;
                    ExpressionScaledResidualResult evaluated = evaluator.evaluate(static_cast<size_t>(iteration), input);
                    localSeriesOperations += evaluated.seriesOperationCount;
                    if (evaluated.status != ExpressionScaledResidualStatus::Success) {
                        pixel.status = evaluated.status;
                        break;
                    }
                    ScaledComplexValue outputBase, actualOutput;
                    arithmeticStatus = formula::makeScaledComplexValue(reference.samples[iteration].next, reference.samples[iteration].rootDefect, outputBase);
                    if (arithmeticStatus == ScaledArithmeticStatus::Success) arithmeticStatus = formula::scaledAdd(outputBase, evaluated.residual, actualOutput);
                    if (arithmeticStatus == ScaledArithmeticStatus::Success) arithmeticStatus = outsideScaled(actualOutput, frame.bailout, escaped);
                    if (arithmeticStatus != ScaledArithmeticStatus::Success) {
                        pixel.status = arithmeticStatus == ScaledArithmeticStatus::ExponentRange ? ExpressionScaledResidualStatus::ExponentRange : ExpressionScaledResidualStatus::Nonfinite;
                        break;
                    }
                    pixel.residual = evaluated.residual;
                    pixel.stateIteration = iteration + 1;
                    if (escaped) {
                        pixel.escapeIteration = iteration + 1;
                        break;
                    }
                    if (iteration + 1 < frame.maxIterations) {
                        formula::ExpressionScaledPrimaryRelativeState primaryRelative;
                        arithmeticStatus = formula::makeExpressionResidualNextPrimaryState(reference.samples[iteration], evaluated.residual, primaryRelative);
                        if (arithmeticStatus == ScaledArithmeticStatus::Success) arithmeticStatus = formula::resetExpressionResidualToExactSample(reference.samples[iteration + 1], primaryRelative, stateDelta);
                        if (arithmeticStatus != ScaledArithmeticStatus::Success) {
                            pixel.status = ExpressionScaledResidualStatus::ExponentRange;
                            break;
                        }
                    }
                }
            }
        };

        std::vector<PixelState> scaledStates;
        size_t localSeriesOperations = 0;
        const Clock::time_point scaledStart = Clock::now();
        renderScaled(scaledStates, localSeriesOperations);
        const double scaledSeconds = elapsed(scaledStart);
        size_t secondSeriesOperations = 0;
        std::vector<PixelState> repeatedStates;
        double repeatedScaledSeconds = scaledSeconds;
        if (frame.benchmark) {
            const Clock::time_point repeatStart = Clock::now();
            renderScaled(repeatedStates, secondSeriesOperations);
            repeatedScaledSeconds = elapsed(repeatStart);
        }
        seriesOperations += localSeriesOperations;
        uint32_t scaledChecksum = 0;
        uint32_t repeatedChecksum = 0;
        bool sawInterior = false;
        bool sawExterior = false;
        for (size_t index = 0; index < pixelCount; ++index) {
            scaledChecksum = scaledChecksum * 16777619u ^ static_cast<uint32_t>(scaledStates[index].escapeIteration);
            if (frame.benchmark) repeatedChecksum = repeatedChecksum * 16777619u ^ static_cast<uint32_t>(repeatedStates[index].escapeIteration);
            if (scaledStates[index].status != ExpressionScaledResidualStatus::Success) ++fallbackPixels;
            sawInterior = sawInterior || scaledStates[index].escapeIteration == frame.maxIterations;
            sawExterior = sawExterior || scaledStates[index].escapeIteration < frame.maxIterations;
        }
        framePixels += pixelCount;
        if (frame.benchmark && (scaledChecksum != repeatedChecksum || localSeriesOperations != secondSeriesOperations)) {
            printf("  scaled benchmark repeat mismatch\n");
            ++failures;
        }
        if ((frame.benchmark || frame.pixel == FormulaParameter::InitialZ) && !(sawInterior && sawExterior)) {
            printf("  frame lacks interior/exterior coverage [%s]\n", frame.name);
            ++failures;
        }

        const mpfr_prec_t oraclePrecision = reference.precision + 256;
        ExpressionOracleContext highReferenceContext(oraclePrecision);
        std::vector<MpfrComplex> highReference;
        highReference.reserve(static_cast<size_t>(frame.maxIterations) + 1);
        if (!configureOracle(highReferenceContext, frame.fixed, frame.pixel, frame.centerReal, frame.centerImaginary)) {
            ++failures;
            continue;
        }
        highReference.emplace_back(oraclePrecision);
        highReference.back().set(highReferenceContext.z);
        MpfrComplex highNext(oraclePrecision);
        bool highReferenceOkay = true;
        for (int iteration = 0; iteration < frame.maxIterations; ++iteration) {
            highReferenceContext.iteration = iteration;
            std::string oracleError;
            if (!ExpressionOracle::evaluate(runtime, highReferenceContext, highNext, &oracleError)) {
                highReferenceOkay = false;
                break;
            }
            highReferenceContext.z.set(highNext);
            highReference.emplace_back(oraclePrecision);
            highReference.back().set(highNext);
        }
        if (!highReferenceOkay) {
            printf("  high reference failed [%s]\n", frame.name);
            ++failures;
            continue;
        }

        auto renderOracle = [&](bool compareResiduals, std::vector<int>& escapes, uint32_t& checksum) {
            escapes.assign(pixelCount, frame.maxIterations);
            checksum = 0;
            mpfr_t highZoom, magnitude, residualMagnitude;
            mpfr_t differenceMagnitude, ratio;
            mpfr_inits2(oraclePrecision, highZoom, magnitude, residualMagnitude, differenceMagnitude, ratio, (mpfr_ptr)0);
            mpfr_set_str(highZoom, "1e500", 10, MPFR_RNDN);
            bool okay = true;
            for (int y = 0; okay && y < frame.height; ++y) {
                for (int x = 0; x < frame.width; ++x) {
                    const size_t index = static_cast<size_t>(y) * frame.width + x;
                    ExpressionOracleContext context(oraclePrecision);
                    if (!configureOracle(context, frame.fixed, frame.pixel, frame.centerReal, frame.centerImaginary)) {
                        okay = false;
                        break;
                    }
                    MpfrComplex exactOffset(oraclePrecision);
                    mpfr_set_si(exactOffset.re, x - frame.width / 2, MPFR_RNDN);
                    mpfr_div(exactOffset.re, exactOffset.re, highZoom, MPFR_RNDN);
                    mpfr_set_si(exactOffset.im, y - frame.height / 2, MPFR_RNDN);
                    mpfr_div(exactOffset.im, exactOffset.im, highZoom, MPFR_RNDN);
                    if (frame.pixel == FormulaParameter::C) {
                        mpfr_add(context.c.re, context.c.re, exactOffset.re, MPFR_RNDN);
                        mpfr_add(context.c.im, context.c.im, exactOffset.im, MPFR_RNDN);
                    } else {
                        mpfr_add(context.z0.re, context.z0.re, exactOffset.re, MPFR_RNDN);
                        mpfr_add(context.z0.im, context.z0.im, exactOffset.im, MPFR_RNDN);
                        context.z.set(context.z0);
                    }
                    int escapeIteration = frame.maxIterations;
                    mpfr_hypot(magnitude, context.z.re, context.z.im, MPFR_RNDN);
                    if (mpfr_cmp_d(magnitude, frame.bailout) > 0) {
                        escapeIteration = 0;
                    } else {
                        MpfrComplex next(oraclePrecision);
                        for (int iteration = 0; iteration < frame.maxIterations; ++iteration) {
                            context.iteration = iteration;
                            std::string oracleError;
                            if (!ExpressionOracle::evaluate(runtime, context, next, &oracleError)) {
                                okay = false;
                                break;
                            }
                            context.z.set(next);
                            mpfr_hypot(magnitude, context.z.re, context.z.im, MPFR_RNDN);
                            if (mpfr_cmp_d(magnitude, frame.bailout) > 0) {
                                escapeIteration = iteration + 1;
                                break;
                            }
                        }
                    }
                    if (!okay) break;
                    escapes[index] = escapeIteration;
                    checksum = checksum * 16777619u ^ static_cast<uint32_t>(escapeIteration);
                    if (escapeIteration != scaledStates[index].escapeIteration || scaledStates[index].status != ExpressionScaledResidualStatus::Success) {
                        okay = false;
                        break;
                    }
                    if (!compareResiduals) continue;

                    const int stateIteration = scaledStates[index].stateIteration;
                    MpfrComplex expectedResidual(oraclePrecision);
                    mpfr_sub(expectedResidual.re, context.z.re, highReference[stateIteration].re, MPFR_RNDN);
                    mpfr_sub(expectedResidual.im, context.z.im, highReference[stateIteration].im, MPFR_RNDN);
                    MpfrComplex actualResidual(oraclePrecision);
                    if (!formula::setMpfrFromScaledValue(actualResidual, scaledStates[index].residual)) {
                        okay = false;
                        break;
                    }
                    mpfr_sub(actualResidual.re, actualResidual.re, expectedResidual.re, MPFR_RNDN);
                    mpfr_sub(actualResidual.im, actualResidual.im, expectedResidual.im, MPFR_RNDN);
                    mpfr_hypot(differenceMagnitude, actualResidual.re, actualResidual.im, MPFR_RNDN);
                    mpfr_hypot(residualMagnitude, expectedResidual.re, expectedResidual.im, MPFR_RNDN);
                    double relativeError = 0.0;
                    if (mpfr_zero_p(residualMagnitude)) {
                        if (!mpfr_zero_p(differenceMagnitude)) {
                            okay = false;
                            break;
                        }
                    } else {
                        mpfr_div(ratio, differenceMagnitude, residualMagnitude, MPFR_RNDU);
                        relativeError = mpfr_get_d(ratio, MPFR_RNDU);
                        maximumResidualRelativeError = std::max(maximumResidualRelativeError, relativeError);
                        if (!(relativeError <= 1e-8)) {
                            okay = false;
                            break;
                        }
                    }
                }
            }
            mpfr_clears(highZoom, magnitude, residualMagnitude, differenceMagnitude, ratio, (mpfr_ptr)0);
            return okay;
        };

        std::vector<int> oracleEscapes;
        uint32_t oracleChecksum = 0;
        const Clock::time_point oracleStart = Clock::now();
        const bool oracleOkay = renderOracle(true, oracleEscapes, oracleChecksum);
        const double oracleSeconds = elapsed(oracleStart);
        if (!oracleOkay || oracleChecksum != scaledChecksum) {
            printf("  e500 frame mismatch [%s] checksums=%08x/%08x\n", frame.name, scaledChecksum, oracleChecksum);
            ++failures;
        }

        double repeatedOracleSeconds = oracleSeconds;
        if (frame.benchmark) {
            std::vector<int> repeatedOracle;
            uint32_t repeatedOracleChecksum = 0;
            const Clock::time_point repeatStart = Clock::now();
            const bool repeatOkay = renderOracle(false, repeatedOracle, repeatedOracleChecksum);
            repeatedOracleSeconds = elapsed(repeatStart);
            if (!repeatOkay || repeatedOracleChecksum != oracleChecksum) {
                printf("  MPFR benchmark repeat mismatch\n");
                ++failures;
            }
            benchmarkReferenceSeconds = referenceSeconds;
            benchmarkScaledSeconds = std::max(scaledSeconds, repeatedScaledSeconds);
            benchmarkMpfrSeconds = std::min(oracleSeconds, repeatedOracleSeconds);
            benchmarkMinimumSpeedup = benchmarkScaledSeconds > 0.0 ? benchmarkMpfrSeconds / benchmarkScaledSeconds : 0.0;
            benchmarkMemory = reference.memoryBytes;
            benchmarkBytesPerSample = reference.sampleCount ? reference.memoryBytes / reference.sampleCount : 0;
            if (!(benchmarkMinimumSpeedup > 1.0) || reference.memoryBytes > size_t{8} * 1024 * 1024) {
                printf("  scaled benchmark gate failed speedup=%.2fx memory=%zu\n", benchmarkMinimumSpeedup, reference.memoryBytes);
                ++failures;
            }
        }
        printf("  %-22s ref/scaled/MPFR %.3f/%.3f/%.3f s pixels=%zu\n", frame.name, referenceSeconds, scaledSeconds, oracleSeconds, pixelCount);
    }

    // Initial escape produces no sample, and is deliberately not dispatched
    // into the evaluator.
    {
        ExpressionProgram canonical, runtime;
        ExpressionContext fixed;
        fixed.c = {};
        if (compile(canonical, "z+1")) {
            ExpressionError error;
            canonical.specialize(fixed, FormulaParameter::InitialZ, runtime, &error);
            ExpressionReferenceBuildRequest request;
            request.canonicalProgram = &canonical;
            request.runtimeProgram = &runtime;
            request.pixelParameter = FormulaParameter::InitialZ;
            request.center.realDecimal = "5";
            request.center.imaginaryDecimal = "0";
            request.bailout = 4.0;
            request.maxIterations = 10;
            request.precision.requestedBits = 256;
            request.precision.minimumBits = 53;
            request.precision.guardBits = 0;
            ExpressionReferenceOrbitResult reference;
            ExpressionScaledResidualEvaluator evaluator;
            if (!formula::buildExpressionReferenceOrbit(request, reference) || !reference.escaped || reference.escapeIteration != 0 || !reference.samples.empty() || !evaluator.prepare(runtime, reference) || evaluator.evaluate(0, ExpressionScaledResidualInput{}).status != ExpressionScaledResidualStatus::InvalidTape) {
                printf("  initial escape/sample-range handling failed\n");
                ++failures;
            }
        }
    }

    const size_t coveredOpcodes = static_cast<size_t>(std::count(opcodeCoverage.begin(), opcodeCoverage.end(), true));
    const double fallbackRate = framePixels ? 100.0 * static_cast<double>(fallbackPixels) / static_cast<double>(framePixels) : 100.0;
    printf("=== expression scaled residual e500 prototype\n");
    printf("  direct binary64 coordinates collapse; scaled coordinates distinct in every frame\n");
    printf("  opcode coverage=%zu fallback=%zu/%zu (%.3f%%) series-ops=%zu\n", coveredOpcodes, fallbackPixels, framePixels, fallbackRate, seriesOperations);
    printf("  max reconstructed residual relative error=%.3g\n", maximumResidualRelativeError);
    printf("  64x44 arithmetic reference/scaled/MPFR %.3f/%.3f/%.3f s min speedup %.2fx\n", benchmarkReferenceSeconds, benchmarkScaledSeconds, benchmarkMpfrSeconds, benchmarkMinimumSpeedup);
    printf("  reference memory=%zu bytes (%zu/sample); nonlinear result is explicitly uncertified\n", benchmarkMemory, benchmarkBytesPerSample);
    printf("  scalar-real diffabs is certified; principal branch/pole per-step and GUI dispatch remain disabled\n");
    printf("  => %s\n\n", failures == 0 ? "PASS" : "CHECK (scaled expression failure)");
    return failures == 0 ? 0 : 1;
}

static int runExpressionTaylorCase() {
    using formula::Complex;
    using formula::ExpressionContext;
    using formula::ExpressionDeepRenderRequest;
    using formula::ExpressionDeepRenderResult;
    using formula::ExpressionError;
    using formula::ExpressionOracle;
    using formula::ExpressionOracleContext;
    using formula::ExpressionProgram;
    using formula::ExpressionReferenceBuildRequest;
    using formula::ExpressionReferenceOrbitResult;
    using formula::ExpressionTaylorJetBuilder;
    using formula::ExpressionTaylorJetEvaluation;
    using formula::ExpressionTaylorJetEvaluator;
    using formula::ExpressionTaylorJetRequest;
    using formula::ExpressionTaylorJetResult;
    using formula::ExpressionTaylorJetStatus;
    using formula::MpfrComplex;
    using formula::ScaledArithmeticStatus;
    using formula::ScaledComplexBall;
    using formula::ScaledComplexValue;
    using formula::ScaledRealValue;

    struct ProgramPair {
        ExpressionProgram canonical;
        ExpressionProgram runtime;
    };
    auto compilePair = [](const char* source, FormulaParameter pixel, const ExpressionContext& fixed, ProgramPair& pair) {
        ExpressionError error;
        return pair.canonical.compile(source, &error) && pair.canonical.specialize(fixed, pixel, pair.runtime, &error);
    };
    auto makeScale = [](int64_t binaryExponent) {
        ScaledComplexValue scale;
        scale.re.mantissa = 0.5;
        scale.re.exponent = binaryExponent + 1;
        return scale;
    };
    auto buildReference = [](const ProgramPair& pair, const ExpressionContext& fixed, FormulaParameter pixel, const char* centerReal, const char* centerImaginary, int iterations, mpfr_prec_t viewBits, ExpressionReferenceOrbitResult& reference) {
        ExpressionReferenceBuildRequest request;
        request.canonicalProgram = &pair.canonical;
        request.runtimeProgram = &pair.runtime;
        request.pixelParameter = pixel;
        request.center.realDecimal = centerReal;
        request.center.imaginaryDecimal = centerImaginary;
        request.fixed = fixed;
        request.bailout = 4.0;
        request.maxIterations = iterations;
        request.precision.viewBits = viewBits;
        request.precision.minimumBits = 128;
        request.precision.guardBits = 64;
        request.precision.maximumBits = 8192;
        request.certificationPrecision = viewBits + 256;
        request.memoryLimitBytes = size_t{1} << 30;
        return formula::buildExpressionReferenceOrbit(request, reference);
    };
    auto contains = [](const MpfrComplex& exact, const ScaledComplexBall& ball) {
        MpfrComplex midpoint(exact.precision());
        mpfr_t radius, difference, otherDifference;
        mpfr_inits2(exact.precision(), radius, difference, otherDifference, (mpfr_ptr)0);
        bool okay = formula::setMpfrFromScaledValue(midpoint, ball.value) && formula::setMpfrFromScaledValue(radius, ball.radius);
        if (okay) {
            mpfr_sub(difference, exact.re, midpoint.re, MPFR_RNDD);
            mpfr_sub(otherDifference, exact.re, midpoint.re, MPFR_RNDU);
            mpfr_abs(difference, difference, MPFR_RNDU);
            mpfr_abs(otherDifference, otherDifference, MPFR_RNDU);
            okay = mpfr_cmp(difference, radius) <= 0 && mpfr_cmp(otherDifference, radius) <= 0;
        }
        if (okay) {
            mpfr_sub(difference, exact.im, midpoint.im, MPFR_RNDD);
            mpfr_sub(otherDifference, exact.im, midpoint.im, MPFR_RNDU);
            mpfr_abs(difference, difference, MPFR_RNDU);
            mpfr_abs(otherDifference, otherDifference, MPFR_RNDU);
            okay = mpfr_cmp(difference, radius) <= 0 && mpfr_cmp(otherDifference, radius) <= 0;
        }
        mpfr_clears(radius, difference, otherDifference, (mpfr_ptr)0);
        return okay;
    };
    auto evaluateOracle = [](const ExpressionProgram& program, const ExpressionContext& fixed, FormulaParameter pixel, const char* centerReal, const char* centerImaginary, const ScaledComplexValue& scale, const ScaledComplexBall& parameterOffset, Complex q, int iterations, mpfr_prec_t precision, MpfrComplex& state) {
        ExpressionOracleContext context(precision);
        context.c.set(fixed.c.real(), fixed.c.imag());
        context.z0.set(fixed.z0.real(), fixed.z0.imag());
        for (size_t parameter = 0; parameter < fixed.parameters.size(); ++parameter) context.parameters[parameter].set(fixed.parameters[parameter].real(), fixed.parameters[parameter].imag());
        MpfrComplex center(precision);
        MpfrComplex d(precision);
        MpfrComplex p(precision);
        MpfrComplex coordinate(precision);
        MpfrComplex next(precision);
        if (!center.set(centerReal, centerImaginary) || !formula::setMpfrFromScaledValue(d, scale) || !formula::setMpfrFromScaledValue(p, parameterOffset.value)) return false;
        mpfr_mul_d(coordinate.re, d.re, q.real(), MPFR_RNDN);
        mpfr_add(coordinate.re, coordinate.re, center.re, MPFR_RNDN);
        mpfr_add(coordinate.re, coordinate.re, p.re, MPFR_RNDN);
        mpfr_mul_d(coordinate.im, d.re, q.imag(), MPFR_RNDN);
        mpfr_add(coordinate.im, coordinate.im, center.im, MPFR_RNDN);
        mpfr_add(coordinate.im, coordinate.im, p.im, MPFR_RNDN);
        if (pixel == FormulaParameter::C)
            context.c.set(coordinate);
        else
            context.z0.set(coordinate);
        context.z.set(context.z0);
        for (int iteration = 0; iteration < iterations; ++iteration) {
            context.iteration = iteration;
            std::string error;
            if (!ExpressionOracle::evaluate(program, context, next, &error)) return false;
            context.z.set(next);
        }
        state.set(context.z);
        return mpfr_number_p(state.re) && mpfr_number_p(state.im);
    };
    auto reconstructLanding = [](const ExpressionReferenceOrbitResult& reference, const ExpressionTaylorJetResult& jet, const ExpressionTaylorJetEvaluation& evaluated, ScaledComplexBall& landing) {
        if (jet.landingSample >= reference.samples.size()) return false;
        const auto& sample = reference.samples[jet.landingSample];
        ScaledComplexBall base;
        if (formula::makeScaledComplexValue(jet.landingUsesSampleOutput ? sample.next : sample.z, jet.landingUsesSampleOutput ? sample.rootDefect : sample.zDefect, base.value) != ScaledArithmeticStatus::Success || formula::certifiedScaledAdd(base, evaluated.residual, landing) != ScaledArithmeticStatus::Success) return false;
        return formula::certifyScaledMpfrExponentRange(landing) == ScaledArithmeticStatus::Success;
    };

    int failures = 0;
    int containmentChecks = 0;
    std::vector<Complex> qValues = {{0.0, -0.0}, {-0.0, 0.0}, {0.25, -0.5}, {-0.6, 0.2}, {0.45, 0.45}, {-0.137, -0.419}, {0.382, -0.071}, {-0.493, 0.311}, {0.064, 0.527}, {0.571, -0.233}, {-0.284, -0.506}, {0.198, 0.613}, {-0.623, 0.109}, {0.999, 0.0}, {-0.999, 0.0}, {0.0, 0.999}, {0.0, -0.999}, {0.499, 0.499}, {-0.499, 0.499}, {0.499, -0.499}, {-0.499, -0.499}};
    auto sameScaled = [](const ScaledRealValue& left, const ScaledRealValue& right) { return left.mantissa == right.mantissa && left.exponent == right.exponent && (left.mantissa != 0.0 || std::signbit(left.mantissa) == std::signbit(right.mantissa)); };
    auto sameEvaluation = [&](const ExpressionTaylorJetEvaluation& left, const ExpressionTaylorJetEvaluation& right) { return left.status == right.status && left.valid == right.valid && left.operationCount == right.operationCount && sameScaled(left.residual.value.re, right.residual.value.re) && sameScaled(left.residual.value.im, right.residual.value.im) && sameScaled(left.residual.radius, right.residual.radius); };
    auto verifyBatch = [&](const char* name, int landing, const ExpressionTaylorJetResult& jet, const std::array<ScaledComplexBall, 4>& q) {
        std::array<ExpressionTaylorJetEvaluation, 4> scalar{};
        std::array<ExpressionTaylorJetEvaluation, 4> batch{};
        bool okay = true;
        for (size_t lane = 0; lane < q.size(); ++lane) okay = ExpressionTaylorJetEvaluator::evaluate(jet, q[lane], scalar[lane]) && okay;
        okay = ExpressionTaylorJetEvaluator::evaluateBatch(jet, q.data(), q.size(), batch.data()) && okay;
        for (size_t lane = 0; lane < q.size(); ++lane) okay = sameEvaluation(scalar[lane], batch[lane]) && okay;
        if (!okay) {
            printf("  Taylor batch/scalar mismatch [%s/%d]\n", name, landing);
            ++failures;
        }
    };
    uint64_t randomState = 0x9e3779b97f4a7c15ULL;
    while (qValues.size() < 32) {
        randomState = randomState * 6364136223846793005ULL + 1;
        const double x = (static_cast<double>((randomState >> 11) & 0x1fffff) / 1048576.0 - 1.0) * 0.7;
        randomState = randomState * 6364136223846793005ULL + 1;
        const double y = (static_cast<double>((randomState >> 11) & 0x1fffff) / 1048576.0 - 1.0) * 0.7;
        if (std::abs(x) + std::abs(y) <= 0.9) qValues.emplace_back(x, y);
    }
    struct FormulaCase {
        const char* name;
        const char* source;
        FormulaParameter pixel;
        ExpressionContext fixed;
        const char* centerReal;
        const char* centerImaginary;
        int64_t scaleExponent;
        mpfr_prec_t viewBits;
        bool firstCoefficientIsScale;
    };
    ExpressionContext quadratic;
    ExpressionContext parameterized;
    parameterized.parameters[0] = {0.125, -0.0625};
    ExpressionContext z0Plane;
    z0Plane.c = {-0.35, 0.2};
    ExpressionContext fixedDivision;
    fixedDivision.z0 = {2.0, 0.0};
    const FormulaCase cases[] = {{"quadratic-c", "z*z+c", FormulaParameter::C, quadratic, "-0.25", "0.1", -1800, 1800, true},
                                 {"parameter-poly", "z*z+c+p0*z", FormulaParameter::C, parameterized, "-0.2", "0.15", -1800, 1800, true},
                                 {"generic-cubic", "z*z*z+c", FormulaParameter::C, quadratic, "-0.1", "0.2", -1800, 1800, true},
                                 {"quadratic-z0", "z*z+c", FormulaParameter::InitialZ, z0Plane, "0.1", "-0.2", -1800, 1800, false},
                                 {"quadratic-c-e1000", "z*z+c", FormulaParameter::C, quadratic, "-0.25", "0.1", -3600, 3600, true},
                                 {"conjugate-c-e500", "conj(z)+c", FormulaParameter::C, quadratic, "0", "0", -1800, 1800, false},
                                 {"real-c-e500", "real(z)+c", FormulaParameter::C, quadratic, "0", "0", -1800, 1800, false},
                                 {"imaginary-c-e500", "imag(z)*complex(0,1)+c", FormulaParameter::C, quadratic, "0", "0", -1800, 1800, false},
                                 {"norm-c-e500", "norm(z)+c", FormulaParameter::C, quadratic, "0", "0", -1800, 1800, false},
                                 {"make-complex-c-e500", "complex(real(z),imag(c))+c", FormulaParameter::C, quadratic, "0", "0", -1800, 1800, false},
                                 {"mixed-real-polynomial-e500", "(z*z+c)*conj(z+c)+c", FormulaParameter::C, quadratic, "0", "0", -1800, 1800, false},
                                 {"conjugate-z0-e500", "conj(z)+c", FormulaParameter::InitialZ, z0Plane, "0", "0", -1800, 1800, false},
                                 {"norm-c-e1000", "norm(z)+c", FormulaParameter::C, quadratic, "0", "0", -3600, 3600, false},
                                 {"abs-real-positive", "abs(real(c))", FormulaParameter::C, quadratic, "2", "0", -10, 512, false},
                                 {"abs-real-negative", "abs(real(c))", FormulaParameter::C, quadratic, "-2", "0", -10, 512, false},
                                 {"exp-moderate", "exp(z)-1+c", FormulaParameter::C, quadratic, "0", "0", -3, 512, true},
                                 {"sin-moderate", "sin(z)+c", FormulaParameter::C, quadratic, "0", "0", -3, 512, true},
                                 {"cos-moderate", "cos(z)-1+c", FormulaParameter::C, quadratic, "0", "0", -3, 512, true},
                                 {"sinh-moderate", "sinh(z)+c", FormulaParameter::C, quadratic, "0", "0", -3, 512, true},
                                 {"cosh-moderate", "cosh(z)-1+c", FormulaParameter::C, quadratic, "0", "0", -3, 512, true},
                                 {"standalone-exp", "exp(c)-1", FormulaParameter::C, quadratic, "0", "0", -3, 512, true},
                                 {"standalone-sin", "sin(c)", FormulaParameter::C, quadratic, "0", "0", -3, 512, true},
                                 {"standalone-cos", "cos(c)-1", FormulaParameter::C, quadratic, "0", "0", -3, 512, false},
                                 {"standalone-sinh", "sinh(c)", FormulaParameter::C, quadratic, "0", "0", -3, 512, true},
                                 {"standalone-cosh", "cosh(c)-1", FormulaParameter::C, quadratic, "0", "0", -3, 512, false},
                                 {"composed-sin-e500", "sin(z*z+c)+c", FormulaParameter::C, quadratic, "0", "0", -1800, 1800, false},
                                 {"mixed-hyperbolic-e500", "sinh(0.25*z)+cosh(0.25*z)-1+c", FormulaParameter::C, quadratic, "0", "0", -1800, 1800, true},
                                 {"sin-z0-e500", "sin(z)+c", FormulaParameter::InitialZ, z0Plane, "0", "0", -1800, 1800, false},
                                 {"divide-polynomial", "z/(c+2)+c", FormulaParameter::C, quadratic, "0", "0", -5, 512, false},
                                 {"reciprocal-polynomial", "1/(z+2)+c", FormulaParameter::C, quadratic, "0", "0", -6, 512, false},
                                 {"tangent", "tan(z)+c", FormulaParameter::C, quadratic, "0", "0", -8, 512, false},
                                 {"hyperbolic-tangent", "tanh(z)+c", FormulaParameter::C, quadratic, "0", "0", -5, 512, false},
                                 {"nested-entire-division", "(sin(z)+c)/(cos(z)+2)", FormulaParameter::C, quadratic, "0", "0", -7, 512, false},
                                 {"reciprocal-z0", "1/(z+2)+c", FormulaParameter::InitialZ, z0Plane, "0", "0", -7, 512, false},
                                 {"standalone-log", "log(c+2)", FormulaParameter::C, quadratic, "0", "0", -7, 512, false},
                                 {"standalone-log10", "log10(c+2)", FormulaParameter::C, quadratic, "0", "0", -7, 512, false},
                                 {"standalone-sqrt", "sqrt(c+2)", FormulaParameter::C, quadratic, "0", "0", -7, 512, false},
                                 {"variable-power", "pow(c+2,c+0.25)", FormulaParameter::C, quadratic, "0", "0", -9, 512, false},
                                 {"large-exponent-power", "pow(c+1.25,20+c)*0.001", FormulaParameter::C, quadratic, "0", "0", -10, 512, false},
                                 {"small-exponent-power", "pow(c+2,0.001+c)-1", FormulaParameter::C, quadratic, "0", "0", -10, 512, false},
                                 {"nested-branch-entire-division", "exp(log(c+2))-1+sin(sqrt(c+3))/(c+4)", FormulaParameter::C, quadratic, "0", "0", -9, 512, false},
                                 {"branch-z0", "log(z+2)+c", FormulaParameter::InitialZ, z0Plane, "0", "0", -9, 512, false},
                                 {"fixed-z0-division-e500", "c/z0", FormulaParameter::C, fixedDivision, "0", "0", -1800, 1800, false},
                                 {"arg-c-e500", "arg(c)", FormulaParameter::C, quadratic, "2", "0.5", -1800, 1800, false},
                                 {"arg-z-plus-c-e500", "arg(z+2)+c", FormulaParameter::C, quadratic, "0", "0.25", -1800, 1800, false},
                                 {"polar-real-imag-e500", "polar(real(c)+2,imag(z))", FormulaParameter::C, quadratic, "-1", "0", -1800, 1800, false},
                                 {"nested-polar-arg-e500", "polar(abs(real(c))+1,arg(c+2))", FormulaParameter::C, quadratic, "0.5", "0.25", -1800, 1800, false}};

    for (const FormulaCase& test : cases) {
        ProgramPair pair;
        ExpressionReferenceOrbitResult reference;
        if (!compilePair(test.source, test.pixel, test.fixed, pair) || !buildReference(pair, test.fixed, test.pixel, test.centerReal, test.centerImaginary, 28, test.viewBits, reference)) {
            printf("  Taylor setup failed [%s]\n", test.name);
            ++failures;
            continue;
        }
        const ScaledComplexValue scale = makeScale(test.scaleExponent);
        ExpressionTaylorJetRequest request;
        request.program = &pair.runtime;
        request.reference = &reference;
        request.pixelParameter = test.pixel;
        request.parameterScale = scale;
        request.minimumOrder = 8;
        request.preferredOrder = 12;
        request.maximumOrder = 20;
        request.minimumLanding = 1;
        request.maximumCandidateIteration = 12;
        request.bailout = 4.0;
        request.accuracyBudget = 0x1p-40;
        ExpressionTaylorJetResult jet;
        const bool entireCandidate = pair.runtime.scaledResidualCapability() == formula::ExpressionScaledResidualCapability::CertifiedEntireCandidate;
        const bool meromorphicCandidate = pair.runtime.scaledResidualCapability() == formula::ExpressionScaledResidualCapability::CertifiedMeromorphicCandidate;
        const bool branchCandidate = pair.runtime.scaledResidualCapability() == formula::ExpressionScaledResidualCapability::CertifiedBranchCandidate;
        const bool realCandidate = pair.runtime.scaledResidualCapability() == formula::ExpressionScaledResidualCapability::CertifiedRealCandidate;
        const bool piecewiseCandidate = pair.runtime.scaledResidualCapability() == formula::ExpressionScaledResidualCapability::CertifiedPiecewiseCandidate;
        const bool realBivariateCandidate = realCandidate || piecewiseCandidate;
        const bool built = ExpressionTaylorJetBuilder::build(request, jet);
        size_t expectedCoefficientCount = static_cast<size_t>(jet.order) + 1;
        if (realBivariateCandidate && !formula::expressionTaylorBivariateMonomialCount(jet.order, expectedCoefficientCount)) expectedCoefficientCount = 0;
        if (!built || jet.landingIteration < 1 || !jet.certified || jet.coefficients.size() != expectedCoefficientCount || jet.monomialCount != expectedCoefficientCount || (realBivariateCandidate && (jet.layout != formula::ExpressionTaylorJetLayout::RealBivariate)) || (entireCandidate && (jet.functionSeriesCount == 0 || jet.functionSeriesOperationCount == 0 || jet.maximumFunctionSeriesOrder < 1 || formula::compareScaledNonnegative(jet.maximumFunctionSeriesTail, ScaledRealValue{0.5, -39}) > 0)) || (meromorphicCandidate && (jet.reciprocalCount == 0 || jet.reciprocalOperationCount == 0 || jet.maximumReciprocalOrder < 8 || jet.minimumDenominatorClearance.isZero())) || (branchCandidate && (jet.branchCompositionCount == 0 || jet.branchCompositionOperationCount == 0 || jet.maximumBranchSeriesOrder < 1 || jet.minimumBranchCutClearance.isZero() || jet.minimumBranchZeroClearance.isZero())) ||
            (piecewiseCandidate && (jet.absBranchCount == 0 || jet.minimumFoldClearance.isZero() || jet.absPositiveCellCount + jet.absNegativeCellCount == 0))) {
            printf("  Taylor build failed [%s]: %s/%s landing=%d\n", test.name, formula::expressionTaylorJetStatusName(jet.status), jet.failureReason.c_str(), jet.landingIteration);
            ++failures;
            continue;
        }
        const int checkedLanding = std::min(8, jet.landingIteration);
        for (int landing = 1; landing <= checkedLanding; ++landing) {
            request.maximumCandidateIteration = landing;
            ExpressionTaylorJetResult prefix;
            if (!ExpressionTaylorJetBuilder::build(request, prefix) || prefix.landingIteration != landing) {
                printf("  Taylor prefix build failed [%s/%d]\n", test.name, landing);
                ++failures;
                break;
            }
            if (landing == 1 && test.pixel == FormulaParameter::C && test.firstCoefficientIsScale && (prefix.coefficients[1].re.mantissa != scale.re.mantissa || prefix.coefficients[1].re.exponent != scale.re.exponent || !prefix.coefficients[1].im.isZero())) {
                printf("  first Taylor coefficient mismatch [%s]\n", test.name);
                ++failures;
            }
            for (Complex qValue : qValues) {
                ScaledComplexBall q;
                if (formula::makeScaledComplexValue(qValue, q.value) != ScaledArithmeticStatus::Success || !formula::expressionTaylorQInsideUnitDisk(q)) {
                    ++failures;
                    continue;
                }
                ExpressionTaylorJetEvaluation evaluated;
                ScaledComplexBall landingBall;
                MpfrComplex exact(reference.certificationPrecision);
                const bool okay = ExpressionTaylorJetEvaluator::evaluate(prefix, q, evaluated) && reconstructLanding(reference, prefix, evaluated, landingBall) && evaluateOracle(pair.runtime, test.fixed, test.pixel, test.centerReal, test.centerImaginary, scale, ScaledComplexBall{}, qValue, landing, reference.certificationPrecision, exact) && contains(exact, landingBall);
                ++containmentChecks;
                if (!okay) {
                    printf("  Taylor MPFR containment failed [%s/%d q=(%.2f,%.2f)]\n", test.name, landing, qValue.real(), qValue.imag());
                    ++failures;
                }
            }
            std::array<ScaledComplexBall, 4> batchQ{};
            bool batchQOkay = true;
            for (size_t lane = 0; lane < batchQ.size(); ++lane) batchQOkay = formula::makeScaledComplexValue(qValues[lane], batchQ[lane].value) == ScaledArithmeticStatus::Success && batchQOkay;
            if (batchQOkay)
                verifyBatch(test.name, landing, prefix, batchQ);
            else
                ++failures;

            batchQ = {};
            batchQ[0].value.re = {0.5, -1800};
            batchQ[0].value.im.mantissa = -0.0;
            batchQ[1].value.re.mantissa = -0.0;
            batchQ[1].value.im = {-0.5, -3600};
            batchQ[2].value.re = {-0.5, -1024};
            batchQ[2].value.im = {0.5, -2048};
            batchQ[3].value.re.mantissa = 0.0;
            batchQ[3].value.im.mantissa = -0.0;
            verifyBatch(test.name, landing, prefix, batchQ);
        }
    }

    // Tile-local P+D*q jets remain relative to the global reference.
    for (FormulaParameter pixel : {FormulaParameter::C, FormulaParameter::InitialZ}) {
        ExpressionContext fixed;
        if (pixel == FormulaParameter::InitialZ) fixed.c = {-0.2, 0.1};
        ProgramPair pair;
        ExpressionReferenceOrbitResult reference;
        const ScaledComplexValue scale = makeScale(-8);
        ScaledComplexBall parameterOffset;
        bool okay = compilePair("z*z+c", pixel, fixed, pair) && buildReference(pair, fixed, pixel, "0", "0", 12, 512, reference) && formula::makeScaledComplexValue(Complex{0.125, -0.0625}, parameterOffset.value) == ScaledArithmeticStatus::Success;
        parameterOffset.radius = {0.5, -99};
        ExpressionTaylorJetRequest request;
        request.program = &pair.runtime;
        request.reference = &reference;
        request.pixelParameter = pixel;
        request.parameterOffset = parameterOffset;
        request.parameterScale = scale;
        request.minimumLanding = 1;
        request.maximumCandidateIteration = 4;
        ExpressionTaylorJetResult jet;
        okay = okay && ExpressionTaylorJetBuilder::build(request, jet) && jet.parameterOffset.value.re.mantissa == parameterOffset.value.re.mantissa && jet.parameterOffset.value.re.exponent == parameterOffset.value.re.exponent && jet.parameterOffset.radius.mantissa == parameterOffset.radius.mantissa && jet.parameterOffset.radius.exponent == parameterOffset.radius.exponent;
        for (Complex qValue : {Complex{0.0, 0.0}, Complex{-0.75, 0.0}, Complex{0.0, -0.75}, Complex{0.375, 0.375}}) {
            ScaledComplexBall exactQ;
            ScaledComplexBall scaleBall;
            ScaledComplexBall displacement;
            ScaledComplexBall pixelOffset;
            ScaledComplexBall localQ;
            ScaledComplexBall reconstructed;
            scaleBall.value = scale;
            okay = okay && formula::makeScaledComplexValue(qValue, exactQ.value) == ScaledArithmeticStatus::Success && formula::certifiedScaledMultiply(scaleBall, exactQ, displacement) == ScaledArithmeticStatus::Success && formula::certifiedScaledAdd(parameterOffset, displacement, pixelOffset) == ScaledArithmeticStatus::Success && formula::makeExpressionTaylorLocalQ(pixelOffset, parameterOffset, scale, localQ) && formula::expressionTaylorQInsideUnitDisk(localQ) && formula::certifiedScaledMultiply(scaleBall, localQ, displacement) == ScaledArithmeticStatus::Success && formula::certifiedScaledAdd(parameterOffset, displacement, reconstructed) == ScaledArithmeticStatus::Success;
            ExpressionTaylorJetEvaluation evaluated;
            ScaledComplexBall landingBall;
            MpfrComplex exact(reference.certificationPrecision);
            okay = okay && ExpressionTaylorJetEvaluator::evaluate(jet, localQ, evaluated) && reconstructLanding(reference, jet, evaluated, landingBall) && evaluateOracle(pair.runtime, fixed, pixel, "0", "0", scale, parameterOffset, qValue, jet.landingIteration, reference.certificationPrecision, exact) && contains(exact, landingBall);
            ++containmentChecks;
        }
        if (!okay) {
            printf("  offset Taylor containment failed [%s]\n", pixel == FormulaParameter::C ? "c" : "z0");
            ++failures;
        }
    }

    // Odd dimensions and all frame corners must remain in the normalized disk.
    {
        ScaledRealValue real, imaginary, error;
        formula::makeScaledRealValue(3.0, real);
        formula::makeScaledRealValue(2.0, imaginary);
        formula::makeScaledRealValue(0x1p-50, error);
        ScaledComplexValue scale;
        bool okay = formula::makeExpressionTaylorFrameScale(real, imaginary, error, scale);
        for (double y : {-2.0, 0.0, 2.0})
            for (double x : {-3.0, 0.0, 3.0}) {
                ScaledComplexBall offset;
                okay = okay && formula::makeScaledComplexValue(Complex{x, y}, offset.value) == ScaledArithmeticStatus::Success;
                offset.radius = error;
                ScaledComplexBall q;
                okay = okay && formula::makeExpressionTaylorNormalizedQ(offset, scale, q) && formula::expressionTaylorQInsideUnitDisk(q);
            }
        if (!okay) {
            printf("  Taylor odd-frame normalization failed\n");
            ++failures;
        }
        scale = {};
        scale.re.mantissa = 0.5;
        scale.re.exponent = std::numeric_limits<int64_t>::min();
        ScaledComplexBall boundaryOffset, boundaryQ;
        boundaryOffset.value.re.mantissa = 0.5;
        boundaryOffset.value.re.exponent = std::numeric_limits<int64_t>::min();
        if (formula::makeExpressionTaylorNormalizedQ(boundaryOffset, scale, boundaryQ)) {
            printf("  Taylor exponent-boundary normalization was not rejected\n");
            ++failures;
        }
    }

    // Unsupported non-holomorphic bytecode, malformed tape, policy limits,
    // cancellation, and exponent guards must reject without weakening bounds.
    {
        struct CapabilityCase {
            const char* source;
            formula::ExpressionScaledResidualCapability expected;
        };
        const CapabilityCase capabilityCases[] = {{"z/(c+2)+c", formula::ExpressionScaledResidualCapability::CertifiedMeromorphicCandidate},
                                                  {"tan(z)+c", formula::ExpressionScaledResidualCapability::CertifiedMeromorphicCandidate},
                                                  {"tanh(z)+c", formula::ExpressionScaledResidualCapability::CertifiedMeromorphicCandidate},
                                                  {"log(z)+c", formula::ExpressionScaledResidualCapability::CertifiedBranchCandidate},
                                                  {"log10(z)+c", formula::ExpressionScaledResidualCapability::CertifiedBranchCandidate},
                                                  {"sqrt(z)+c", formula::ExpressionScaledResidualCapability::CertifiedBranchCandidate},
                                                  {"pow(z,2.5)+c", formula::ExpressionScaledResidualCapability::CertifiedBranchCandidate},
                                                  {"conj(z)+c", formula::ExpressionScaledResidualCapability::CertifiedRealCandidate},
                                                  {"real(z)+c", formula::ExpressionScaledResidualCapability::CertifiedRealCandidate},
                                                  {"imag(z)+c", formula::ExpressionScaledResidualCapability::CertifiedRealCandidate},
                                                  {"norm(z)+c", formula::ExpressionScaledResidualCapability::CertifiedRealCandidate},
                                                  {"complex(real(z),imag(c))+c", formula::ExpressionScaledResidualCapability::CertifiedRealCandidate},
                                                  {"abs(real(z))+c", formula::ExpressionScaledResidualCapability::CertifiedPiecewiseCandidate},
                                                  {"abs(complex(real(z),0))+c", formula::ExpressionScaledResidualCapability::CertifiedPiecewiseCandidate},
                                                  {"sqr(complex(abs(real(z)),abs(imag(z))))+c", formula::ExpressionScaledResidualCapability::CertifiedPiecewiseCandidate},
                                                  {"sin(conj(z))+c", formula::ExpressionScaledResidualCapability::Unsupported},
                                                  {"conj(z)/(c+2)+c", formula::ExpressionScaledResidualCapability::Unsupported},
                                                  {"arg(z)+c", formula::ExpressionScaledResidualCapability::CertifiedRealCandidate},
                                                  {"polar(real(c)+2,imag(z))", formula::ExpressionScaledResidualCapability::CertifiedRealCandidate},
                                                  {"polar(abs(real(c))+1,arg(c+2))", formula::ExpressionScaledResidualCapability::CertifiedPiecewiseCandidate},
                                                  {"log(abs(z)+2)+c", formula::ExpressionScaledResidualCapability::BranchSensitive},
                                                  {"log(norm(z)+2)+c", formula::ExpressionScaledResidualCapability::BranchSensitive},
                                                  {"abs(z)+c", formula::ExpressionScaledResidualCapability::Unsupported},
                                                  {"polar(1,z)+c", formula::ExpressionScaledResidualCapability::Unsupported}};
        for (const CapabilityCase& test : capabilityCases) {
            ProgramPair pair;
            ExpressionContext fixed;
            if (!compilePair(test.source, FormulaParameter::C, fixed, pair) || pair.runtime.scaledResidualCapability() != test.expected) {
                printf("  Taylor capability tier mismatch [%s]\n", test.source);
                ++failures;
            }
        }
    }
    {
        ProgramPair pair;
        ExpressionContext fixed;
        ExpressionReferenceOrbitResult reference;
        ExpressionTaylorJetRequest request;
        ExpressionTaylorJetResult jet;
        if (!compilePair("abs(real(c))", FormulaParameter::C, fixed, pair) || !buildReference(pair, fixed, FormulaParameter::C, "0", "0", 8, 512, reference)) {
            printf("  absolute-value fold Taylor setup failed\n");
            ++failures;
        } else {
            request.program = &pair.runtime;
            request.reference = &reference;
            request.pixelParameter = FormulaParameter::C;
            request.parameterScale = makeScale(-4);
            request.minimumOrder = 8;
            request.preferredOrder = 8;
            request.maximumOrder = 8;
            request.minimumLanding = 1;
            request.maximumCandidateIteration = 4;
            request.bailout = 4.0;
            request.accuracyBudget = 0x1p-40;
            if (ExpressionTaylorJetBuilder::build(request, jet) || jet.status != ExpressionTaylorJetStatus::BranchRejected || !jet.foldRejected || jet.foldRejectionIteration != 0 || jet.foldRejectionReason.empty()) {
                printf("  absolute-value fold Taylor rejection failed status=%s iteration=%d reason=%s\n", formula::expressionTaylorJetStatusName(jet.status), jet.foldRejectionIteration, jet.foldRejectionReason.c_str());
                ++failures;
            }
        }
    }
    {
        size_t count = 0;
        size_t qIndex = 0;
        size_t conjugateIndex = 0;
        int qDegree = -1;
        int conjugateDegree = -1;
        const bool okay = formula::expressionTaylorBivariateMonomialCount(8, count) && count == 45 && formula::expressionTaylorBivariateIndex(8, 1, 0, qIndex) && formula::expressionTaylorBivariateIndex(8, 0, 1, conjugateIndex) && qIndex == 1 && conjugateIndex == 2 && formula::expressionTaylorBivariateExponents(8, conjugateIndex, qDegree, conjugateDegree) && qDegree == 0 && conjugateDegree == 1;
        if (!okay) {
            printf("  real-bivariate Taylor index mapping failed\n");
            ++failures;
        }
    }
    {
        auto buildIdentityJet = [&](const char* source, ExpressionTaylorJetResult& jet) {
            ProgramPair pair;
            ExpressionReferenceOrbitResult reference;
            ExpressionContext fixed;
            if (!compilePair(source, FormulaParameter::C, fixed, pair) || !buildReference(pair, fixed, FormulaParameter::C, "0", "0", 2, 512, reference)) return false;
            ExpressionTaylorJetRequest request;
            request.program = &pair.runtime;
            request.reference = &reference;
            request.pixelParameter = FormulaParameter::C;
            request.parameterScale = makeScale(-100);
            request.minimumLanding = 1;
            request.maximumCandidateIteration = 1;
            return ExpressionTaylorJetBuilder::build(request, jet) && jet.layout == formula::ExpressionTaylorJetLayout::RealBivariate;
        };
        auto sameScaled = [](const ScaledRealValue& left, const ScaledRealValue& right) { return left.mantissa == right.mantissa && left.exponent == right.exponent; };
        auto containsZero = [](const ScaledComplexBall& ball) {
            ScaledRealValue real = ball.value.re;
            ScaledRealValue imaginary = ball.value.im;
            real.mantissa = std::abs(real.mantissa);
            imaginary.mantissa = std::abs(imaginary.mantissa);
            return formula::compareScaledNonnegative(real, ball.radius) <= 0 && formula::compareScaledNonnegative(imaginary, ball.radius) <= 0;
        };
        auto symmetryHolds = [&](const ExpressionTaylorJetResult& jet) {
            for (size_t index = 0; index < jet.coefficients.size(); ++index) {
                int qDegree = 0;
                int conjugateDegree = 0;
                size_t swappedIndex = 0;
                if (!formula::expressionTaylorBivariateExponents(jet.order, index, qDegree, conjugateDegree) || !formula::expressionTaylorBivariateIndex(jet.order, conjugateDegree, qDegree, swappedIndex)) return false;
                ScaledComplexBall left{jet.coefficients[index], jet.coefficientRadii[index]};
                ScaledComplexBall right{jet.coefficients[swappedIndex], jet.coefficientRadii[swappedIndex]};
                if (formula::scaledNegate(right.value.im, right.value.im) != ScaledArithmeticStatus::Success) return false;
                ScaledComplexBall difference;
                if (formula::certifiedScaledSubtract(left, right, difference) != ScaledArithmeticStatus::Success || !containsZero(difference)) return false;
            }
            return true;
        };

        size_t qIndex = 0;
        size_t conjugateIndex = 0;
        size_t normIndex = 0;
        bool okay = formula::expressionTaylorBivariateIndex(8, 1, 0, qIndex) && formula::expressionTaylorBivariateIndex(8, 0, 1, conjugateIndex) && formula::expressionTaylorBivariateIndex(8, 1, 1, normIndex);
        const ScaledComplexValue scale = makeScale(-100);
        ScaledRealValue halfScale;
        ScaledRealValue normScale;
        okay = okay && formula::scaledDivideByDouble(scale.re, 2.0, halfScale) == ScaledArithmeticStatus::Success && formula::scaledMultiply(scale.re, scale.re, normScale) == ScaledArithmeticStatus::Success;

        ExpressionTaylorJetResult conjugateJet;
        okay = okay && buildIdentityJet("conj(c)", conjugateJet) && conjugateJet.coefficients[qIndex].isZero() && sameScaled(conjugateJet.coefficients[conjugateIndex].re, scale.re) && conjugateJet.coefficients[conjugateIndex].im.isZero() && std::signbit(conjugateJet.coefficients[conjugateIndex].im.mantissa);
        ScaledComplexBall signedAxisQ;
        ExpressionTaylorJetEvaluation conjugateEvaluation;
        okay = okay && formula::makeScaledComplexValue(Complex{0.5, 0.0}, signedAxisQ.value) == ScaledArithmeticStatus::Success && ExpressionTaylorJetEvaluator::evaluate(conjugateJet, signedAxisQ, conjugateEvaluation) && conjugateEvaluation.residual.value.im.isZero() && std::signbit(conjugateEvaluation.residual.value.im.mantissa);

        ExpressionTaylorJetResult realJet;
        okay = okay && buildIdentityJet("real(c)", realJet) && sameScaled(realJet.coefficients[qIndex].re, halfScale) && sameScaled(realJet.coefficients[conjugateIndex].re, halfScale) && symmetryHolds(realJet);

        ExpressionTaylorJetResult imaginaryJet;
        okay = okay && buildIdentityJet("imag(c)", imaginaryJet) && sameScaled(imaginaryJet.coefficients[qIndex].im, ScaledRealValue{-halfScale.mantissa, halfScale.exponent}) && sameScaled(imaginaryJet.coefficients[conjugateIndex].im, halfScale) && symmetryHolds(imaginaryJet);

        ExpressionTaylorJetResult normJet;
        okay = okay && buildIdentityJet("norm(c)", normJet) && sameScaled(normJet.coefficients[normIndex].re, normScale) && normJet.coefficients[normIndex].im.isZero() && normJet.bivariateConvolutionOperationCount > 0 && symmetryHolds(normJet);

        ExpressionTaylorJetResult makeComplexJet;
        okay = okay && buildIdentityJet("complex(real(c),imag(c))", makeComplexJet) && sameScaled(makeComplexJet.coefficients[qIndex].re, scale.re) && makeComplexJet.coefficients[qIndex].im.isZero() && containsZero(ScaledComplexBall{makeComplexJet.coefficients[conjugateIndex], makeComplexJet.coefficientRadii[conjugateIndex]});
        if (!okay) {
            printf("  real-bivariate coefficient identities failed\n");
            ++failures;
        }
    }

    // Principal-branch clearance is certified over the entire q-frame. Point
    // cut metadata, including signed-zero lip labels, is never sufficient.
    {
        auto checkBranchFrame = [&](const char* name, const char* source, const char* centerReal, const char* centerImaginary, int64_t scaleExponent, bool expectAccepted) {
            ProgramPair pair;
            ExpressionReferenceOrbitResult reference;
            ExpressionContext fixed;
            bool okay = compilePair(source, FormulaParameter::C, fixed, pair) && buildReference(pair, fixed, FormulaParameter::C, centerReal, centerImaginary, 4, 512, reference);
            ExpressionTaylorJetRequest request;
            request.program = &pair.runtime;
            request.reference = &reference;
            request.pixelParameter = FormulaParameter::C;
            request.parameterScale = makeScale(scaleExponent);
            request.minimumLanding = 1;
            request.maximumCandidateIteration = 1;
            request.accuracyBudget = 0x1p-40;
            ExpressionTaylorJetResult jet;
            const bool built = okay && ExpressionTaylorJetBuilder::build(request, jet);
            if (expectAccepted) {
                okay = built && jet.branchCompositionCount > 0 && jet.branchCompositionOperationCount > 0 && jet.maximumBranchSeriesOrder > 0 && !jet.minimumBranchCutClearance.isZero() && !jet.minimumBranchZeroClearance.isZero() && !jet.branchRejected;
            } else {
                okay = okay && !built && jet.status == ExpressionTaylorJetStatus::BranchRejected && jet.branchRejected && (jet.minimumBranchCutClearance.isZero() || jet.minimumBranchZeroClearance.isZero()) && jet.failureReason.find("branch input") != std::string::npos;
            }
            if (!okay) {
                printf("  branch-clearance proof failed [%s] built=%d status=%s reason=%s cut=(%.6g,e%lld) zero=(%.6g,e%lld)\n", name, built ? 1 : 0, formula::expressionTaylorJetStatusName(jet.status), jet.failureReason.c_str(), jet.minimumBranchCutClearance.mantissa, (long long)jet.minimumBranchCutClearance.exponent, jet.minimumBranchZeroClearance.mantissa, (long long)jet.minimumBranchZeroClearance.exponent);
                ++failures;
            }
        };
        checkBranchFrame("log-near-cut-safe", "log(c)*0.001", "-2", "0.04", -5, true);
        checkBranchFrame("log-near-cut-overlap", "log(c)*0.001", "-2", "0.02", -5, false);
        checkBranchFrame("log10-negative-neighborhood", "log10(c)*0.001", "-2", "0.04", -5, true);
        checkBranchFrame("sqrt-near-cut-safe", "sqrt(c)*0.001", "-2", "-0.04", -5, true);
        checkBranchFrame("power-near-cut-overlap", "pow(c,0.25)*0.001", "-2", "-0.02", -5, false);
        checkBranchFrame("exact-upper-lip", "log(c)*0.001", "-2", "0", -20, false);
        checkBranchFrame("exact-lower-lip", "log(c)*0.001", "-2", "-0", -20, false);
        checkBranchFrame("origin", "sqrt(c)*0.001", "0", "0", -20, false);
        checkBranchFrame("zero-overlap", "log(c)*0.001", "0.001", "0", -8, false);

        for (const char* source : {"log(c)", "log10(c)", "sqrt(c)", "pow(c,0.25)"}) {
            ExpressionProgram program;
            ExpressionError error;
            ExpressionOracleContext upper(512);
            ExpressionOracleContext lower(512);
            MpfrComplex upperValue(512);
            MpfrComplex lowerValue(512);
            upper.c.set(-2.0, 0.0);
            lower.c.set(-2.0, -0.0);
            std::string upperError, lowerError;
            const bool okay = program.compile(source, &error) && ExpressionOracle::evaluate(program, upper, upperValue, &upperError) && ExpressionOracle::evaluate(program, lower, lowerValue, &lowerError) && mpfr_sgn(upperValue.im) > 0 && mpfr_sgn(lowerValue.im) < 0 && !mpfr_signbit(upperValue.im) && mpfr_signbit(lowerValue.im);
            if (!okay) {
                printf("  principal signed-zero lip semantics failed [%s]\n", source);
                ++failures;
            }
        }
        {
            ExpressionProgram program;
            ExpressionError error;
            ExpressionOracleContext upper(512);
            ExpressionOracleContext lower(512);
            MpfrComplex upperValue(512);
            MpfrComplex lowerValue(512);
            upper.c.set(-2.0, 0.0);
            lower.c.set(-2.0, -0.0);
            std::string upperError, lowerError;
            const bool okay = program.compile("arg(c)", &error) && ExpressionOracle::evaluate(program, upper, upperValue, &upperError) && ExpressionOracle::evaluate(program, lower, lowerValue, &lowerError) && mpfr_sgn(upperValue.re) > 0 && mpfr_sgn(lowerValue.re) < 0 && mpfr_zero_p(upperValue.im) && mpfr_zero_p(lowerValue.im);
            if (!okay) {
                printf("  arg signed-zero lip semantics failed\n");
                ++failures;
            }
        }

        auto checkArgFrame = [&](const char* name, const char* centerReal, const char* centerImaginary, int64_t scaleExponent, bool expectAccepted) {
            ProgramPair pair;
            ExpressionReferenceOrbitResult reference;
            ExpressionContext fixed;
            bool okay = compilePair("0.001*arg(c)", FormulaParameter::C, fixed, pair) && buildReference(pair, fixed, FormulaParameter::C, centerReal, centerImaginary, 4, 512, reference);
            ExpressionTaylorJetRequest request;
            request.program = &pair.runtime;
            request.reference = &reference;
            request.pixelParameter = FormulaParameter::C;
            request.parameterScale = makeScale(scaleExponent);
            request.minimumLanding = 1;
            request.maximumCandidateIteration = 1;
            ExpressionTaylorJetResult jet;
            const bool built = okay && ExpressionTaylorJetBuilder::build(request, jet);
            if (expectAccepted) {
                okay = built && jet.layout == formula::ExpressionTaylorJetLayout::RealBivariate && jet.argCompositionCount > 0 && jet.branchCompositionCount > 0 && jet.branchCompositionOperationCount > 0 && !jet.minimumBranchCutClearance.isZero() && !jet.minimumBranchZeroClearance.isZero() && !jet.branchRejected;
            } else {
                okay = okay && !built && jet.status == ExpressionTaylorJetStatus::BranchRejected && jet.branchRejected && !jet.argRejectionReason.empty();
            }
            if (!okay) {
                printf("  arg-clearance proof failed [%s] built=%d status=%s reason=%s cut=(%.6g,e%lld) zero=(%.6g,e%lld)\n", name, built ? 1 : 0, formula::expressionTaylorJetStatusName(jet.status), jet.failureReason.c_str(), jet.minimumBranchCutClearance.mantissa, (long long)jet.minimumBranchCutClearance.exponent, jet.minimumBranchZeroClearance.mantissa, (long long)jet.minimumBranchZeroClearance.exponent);
                ++failures;
            }
        };
        checkArgFrame("arg-near-cut-safe", "-2", "0.04", -5, true);
        checkArgFrame("arg-near-cut-overlap", "-2", "0.02", -5, false);
        checkArgFrame("arg-upper-lip", "-2", "0", -20, false);
        checkArgFrame("arg-lower-lip", "-2", "-0", -20, false);
        checkArgFrame("arg-origin", "0", "0", -20, false);

        ProgramPair largePower;
        ExpressionReferenceOrbitResult largeReference;
        ExpressionContext fixed;
        bool largeOkay = compilePair("pow(c+1.001,1000+c)-1", FormulaParameter::C, fixed, largePower) && buildReference(largePower, fixed, FormulaParameter::C, "0", "0", 4, 512, largeReference);
        ExpressionTaylorJetRequest largeRequest;
        largeRequest.program = &largePower.runtime;
        largeRequest.reference = &largeReference;
        largeRequest.pixelParameter = FormulaParameter::C;
        largeRequest.parameterScale = makeScale(-24);
        largeRequest.minimumLanding = 1;
        largeRequest.maximumCandidateIteration = 1;
        ExpressionTaylorJetResult largeJet;
        largeOkay = largeOkay && !ExpressionTaylorJetBuilder::build(largeRequest, largeJet) && largeJet.status == ExpressionTaylorJetStatus::AccuracyBudget;
        if (!largeOkay) {
            printf("  large-exponent power fallback failed status=%s reason=%s\n", formula::expressionTaylorJetStatusName(largeJet.status), largeJet.failureReason.c_str());
            ++failures;
        }
    }

    // Polar requires real-valued operands and a radius that is either exactly
    // zero or strictly positive over the complete q-frame.
    {
        auto checkPolarFrame = [&](const char* name, const char* source, const char* centerReal, const char* centerImaginary, int64_t scaleExponent, bool expectAccepted, bool expectExactZero) {
            ProgramPair pair;
            ExpressionReferenceOrbitResult reference;
            ExpressionContext fixed;
            bool okay = compilePair(source, FormulaParameter::C, fixed, pair) && buildReference(pair, fixed, FormulaParameter::C, centerReal, centerImaginary, 4, std::max<mpfr_prec_t>(512, -scaleExponent + 128), reference);
            ExpressionTaylorJetRequest request;
            request.program = &pair.runtime;
            request.reference = &reference;
            request.pixelParameter = FormulaParameter::C;
            request.parameterScale = makeScale(scaleExponent);
            request.minimumLanding = 1;
            request.maximumCandidateIteration = 1;
            ExpressionTaylorJetResult jet;
            const bool built = okay && ExpressionTaylorJetBuilder::build(request, jet);
            if (expectAccepted) {
                okay = built && jet.layout == formula::ExpressionTaylorJetLayout::RealBivariate && jet.polarCompositionCount > 0 && !jet.polarRejected && (expectExactZero ? jet.minimumPolarRadiusClearance.isZero() : !jet.minimumPolarRadiusClearance.isZero() && jet.functionSeriesCount >= 2);
            } else {
                okay = okay && !built && jet.status == ExpressionTaylorJetStatus::BranchRejected && jet.polarRejected && !jet.polarRejectionReason.empty();
            }
            if (!okay) {
                printf("  polar domain proof failed [%s] built=%d status=%s reason=%s radius=(%.6g,e%lld) count=%llu\n", name, built ? 1 : 0, formula::expressionTaylorJetStatusName(jet.status), jet.failureReason.c_str(), jet.minimumPolarRadiusClearance.mantissa, (long long)jet.minimumPolarRadiusClearance.exponent, (unsigned long long)jet.polarCompositionCount);
                ++failures;
            }
        };
        checkPolarFrame("positive-radius", "polar(real(c)+2,imag(c))", "0", "0", -8, true, false);
        checkPolarFrame("exact-zero-radius", "polar(0,real(c))", "1e500", "-0", -1800, true, true);
        checkPolarFrame("near-zero-safe", "polar(real(c)+0.01,imag(c))", "0", "0", -12, true, false);
        checkPolarFrame("radius-crosses-negative", "polar(real(c)+0.01,imag(c))", "0", "0", -5, false, false);
        {
            ProgramPair pair;
            ExpressionReferenceOrbitResult reference;
            ExpressionContext fixed;
            bool okay = compilePair("polar(1,real(c))", FormulaParameter::C, fixed, pair) && buildReference(pair, fixed, FormulaParameter::C, "1e500", "0", 4, 1928, reference);
            ExpressionTaylorJetRequest request;
            request.program = &pair.runtime;
            request.reference = &reference;
            request.pixelParameter = FormulaParameter::C;
            request.parameterScale = makeScale(-1800);
            request.minimumLanding = 1;
            request.maximumCandidateIteration = 1;
            ExpressionTaylorJetResult jet;
            okay = okay && !ExpressionTaylorJetBuilder::build(request, jet) && (jet.status == ExpressionTaylorJetStatus::AccuracyBudget || jet.status == ExpressionTaylorJetStatus::ExponentRange);
            if (!okay) {
                printf("  polar huge-angle fallback failed status=%s reason=%s\n", formula::expressionTaylorJetStatusName(jet.status), jet.failureReason.c_str());
                ++failures;
            }
        }

        ExpressionProgram nonrealRadius;
        ExpressionProgram nonrealAngle;
        ExpressionError error;
        const bool nonrealOkay = nonrealRadius.compile("polar(c,0)", &error) && nonrealAngle.compile("polar(1,c)", &error) && nonrealRadius.scaledResidualCapability() == formula::ExpressionScaledResidualCapability::Unsupported && nonrealAngle.scaledResidualCapability() == formula::ExpressionScaledResidualCapability::Unsupported;
        if (!nonrealOkay) {
            printf("  polar nonreal operand capability failed\n");
            ++failures;
        }

        ExpressionProgram signedPolar;
        ExpressionOracleContext negativeZeroContext(512);
        ExpressionOracleContext positiveZeroContext(512);
        MpfrComplex negativeZeroOutput(512);
        MpfrComplex positiveZeroOutput(512);
        std::string signedError;
        negativeZeroContext.parameters[0].set(-0.0, 0.0);
        negativeZeroContext.parameters[1].set(-0.0, 0.0);
        positiveZeroContext.parameters[0].set(0.0, 0.0);
        positiveZeroContext.parameters[1].set(-0.0, 0.0);
        const bool signedPolarOkay = signedPolar.compile("polar(p0,p1)", &error) && ExpressionOracle::evaluate(signedPolar, negativeZeroContext, negativeZeroOutput, &signedError) && ExpressionOracle::evaluate(signedPolar, positiveZeroContext, positiveZeroOutput, &signedError) && mpfr_zero_p(negativeZeroOutput.re) && mpfr_signbit(negativeZeroOutput.re) && mpfr_zero_p(negativeZeroOutput.im) && !mpfr_signbit(negativeZeroOutput.im) && mpfr_zero_p(positiveZeroOutput.re) && !mpfr_signbit(positiveZeroOutput.re) && mpfr_zero_p(positiveZeroOutput.im) && mpfr_signbit(positiveZeroOutput.im);
        if (!signedPolarOkay) {
            printf("  polar signed-zero semantics failed\n");
            ++failures;
        }

        ProgramPair resourcePair;
        ExpressionReferenceOrbitResult resourceReference;
        ExpressionContext fixed;
        bool resourceOkay = compilePair("polar(real(c)+1,arg(c+2))", FormulaParameter::C, fixed, resourcePair) && buildReference(resourcePair, fixed, FormulaParameter::C, "0", "0", 8, 512, resourceReference);
        ExpressionTaylorJetRequest resourceRequest;
        resourceRequest.program = &resourcePair.runtime;
        resourceRequest.reference = &resourceReference;
        resourceRequest.pixelParameter = FormulaParameter::C;
        resourceRequest.parameterScale = makeScale(-20);
        resourceRequest.minimumLanding = 1;
        resourceRequest.maximumCandidateIteration = 4;
        ExpressionTaylorJetResult baseline;
        resourceOkay = resourceOkay && ExpressionTaylorJetBuilder::build(resourceRequest, baseline) && baseline.memoryBytes > 1 && baseline.argCompositionCount > 0 && baseline.polarCompositionCount > 0;
        if (resourceOkay) {
            resourceRequest.memoryLimitBytes = baseline.memoryBytes - 1;
            ExpressionTaylorJetResult constrained;
            resourceOkay = !ExpressionTaylorJetBuilder::build(resourceRequest, constrained) && constrained.status == ExpressionTaylorJetStatus::ResourceLimit && constrained.coefficients.empty() && constrained.coefficientRadii.empty();
        }
        resourceRequest.memoryLimitBytes = size_t{1} << 30;
        resourceRequest.shouldCancel = [] { return true; };
        ExpressionTaylorJetResult cancelled;
        resourceOkay = resourceOkay && !ExpressionTaylorJetBuilder::build(resourceRequest, cancelled) && cancelled.status == ExpressionTaylorJetStatus::Cancelled && cancelled.coefficients.empty() && cancelled.coefficientRadii.empty();
        if (!resourceOkay) {
            printf("  Arg/Polar Taylor resource policy failed\n");
            ++failures;
        }
    }

    // Branch workspace accounting includes the normalization, logarithm
    // series, power-product, and nested entire-composition buffers.
    {
        ProgramPair pair;
        ExpressionReferenceOrbitResult reference;
        ExpressionContext fixed;
        bool okay = compilePair("pow(c+2,c+0.25)", FormulaParameter::C, fixed, pair) && buildReference(pair, fixed, FormulaParameter::C, "0", "0", 8, 512, reference);
        ExpressionTaylorJetRequest request;
        request.program = &pair.runtime;
        request.reference = &reference;
        request.pixelParameter = FormulaParameter::C;
        request.parameterScale = makeScale(-10);
        request.minimumLanding = 1;
        request.maximumCandidateIteration = 4;
        ExpressionTaylorJetResult baseline;
        okay = okay && ExpressionTaylorJetBuilder::build(request, baseline) && baseline.memoryBytes > 1 && baseline.branchCompositionCount > 0;
        if (okay) {
            request.memoryLimitBytes = baseline.memoryBytes - 1;
            ExpressionTaylorJetResult constrained;
            okay = !ExpressionTaylorJetBuilder::build(request, constrained) && constrained.status == ExpressionTaylorJetStatus::ResourceLimit && constrained.coefficients.empty() && constrained.coefficientRadii.empty();
        }
        if (!okay) {
            printf("  branch Taylor workspace accounting failed\n");
            ++failures;
        }
    }
    // Real-bivariate workspace accounting includes triangular node storage,
    // retained best candidates, and the quadratic convolution workset.
    {
        ProgramPair pair;
        ExpressionReferenceOrbitResult reference;
        ExpressionContext fixed;
        bool okay = compilePair("(z*z+c)*conj(z+c)+norm(z)", FormulaParameter::C, fixed, pair) && buildReference(pair, fixed, FormulaParameter::C, "0", "0", 24, 512, reference);
        ExpressionTaylorJetRequest request;
        request.program = &pair.runtime;
        request.reference = &reference;
        request.pixelParameter = FormulaParameter::C;
        request.parameterScale = makeScale(-600);
        request.minimumLanding = 1;
        request.maximumCandidateIteration = 12;
        ExpressionTaylorJetResult baseline;
        okay = okay && ExpressionTaylorJetBuilder::build(request, baseline) && baseline.layout == formula::ExpressionTaylorJetLayout::RealBivariate && baseline.memoryBytes > 1 && baseline.bivariateConvolutionOperationCount > 0;
        if (okay) {
            request.memoryLimitBytes = baseline.memoryBytes - 1;
            ExpressionTaylorJetResult constrained;
            okay = !ExpressionTaylorJetBuilder::build(request, constrained) && constrained.status == ExpressionTaylorJetStatus::ResourceLimit && constrained.coefficients.empty() && constrained.coefficientRadii.empty();
        }
        request.memoryLimitBytes = size_t{1} << 30;
        request.maximumBivariateOrder = 7;
        ExpressionTaylorJetResult invalidLow;
        okay = okay && !ExpressionTaylorJetBuilder::build(request, invalidLow) && invalidLow.status == ExpressionTaylorJetStatus::InvalidRequest;
        request.maximumBivariateOrder = formula::ExpressionTaylorMaximumBivariateOrder + 1;
        ExpressionTaylorJetResult invalidHigh;
        okay = okay && !ExpressionTaylorJetBuilder::build(request, invalidHigh) && invalidHigh.status == ExpressionTaylorJetStatus::InvalidRequest;
        request.maximumBivariateOrder = formula::ExpressionTaylorMaximumBivariateOrder;
        request.parameterScale.re.exponent = static_cast<int64_t>(mpfr_get_emax());
        ExpressionTaylorJetResult exponentBoundary;
        okay = okay && !ExpressionTaylorJetBuilder::build(request, exponentBoundary) && exponentBoundary.status == ExpressionTaylorJetStatus::ExponentRange;
        if (!okay) {
            printf("  real-bivariate Taylor workspace policy failed\n");
            ++failures;
        }
    }
    {
        ProgramPair pair;
        ExpressionReferenceOrbitResult reference;
        ExpressionContext fixed;
        bool okay = compilePair("z*z+c", FormulaParameter::C, fixed, pair) && buildReference(pair, fixed, FormulaParameter::C, "0", "0", 20, 512, reference);
        ExpressionTaylorJetRequest request;
        request.program = &pair.runtime;
        request.reference = &reference;
        request.pixelParameter = FormulaParameter::C;
        request.parameterScale = makeScale(-600);
        request.minimumLanding = 1;
        ExpressionTaylorJetResult jet;
        ExpressionReferenceOrbitResult malformed = reference;
        malformed.samples[0].rootNode = UINT16_MAX;
        request.reference = &malformed;
        okay = okay && !ExpressionTaylorJetBuilder::build(request, jet) && jet.status == ExpressionTaylorJetStatus::InvalidTape;
        request.reference = &reference;
        request.minimumOrder = 7;
        okay = okay && !ExpressionTaylorJetBuilder::build(request, jet) && jet.status == ExpressionTaylorJetStatus::InvalidRequest;
        request.minimumOrder = 8;
        request.memoryLimitBytes = 1;
        okay = okay && !ExpressionTaylorJetBuilder::build(request, jet) && jet.status == ExpressionTaylorJetStatus::ResourceLimit;
        request.memoryLimitBytes = size_t{1} << 30;
        request.shouldCancel = [] { return true; };
        okay = okay && !ExpressionTaylorJetBuilder::build(request, jet) && jet.status == ExpressionTaylorJetStatus::Cancelled;
        request.shouldCancel = {};
        request.parameterScale.re.exponent = static_cast<int64_t>(mpfr_get_emax());
        okay = okay && !ExpressionTaylorJetBuilder::build(request, jet) && jet.status == ExpressionTaylorJetStatus::ExponentRange;
        if (!okay) {
            printf("  Taylor rejection policy failed\n");
            ++failures;
        }
    }

    // Meromorphic workspace accounting includes reciprocal/convolution
    // buffers and the retained best rejected-order candidate.
    {
        ProgramPair pair;
        ExpressionReferenceOrbitResult reference;
        ExpressionContext fixed;
        bool okay = compilePair("1/(z+2)+c", FormulaParameter::C, fixed, pair) && buildReference(pair, fixed, FormulaParameter::C, "0", "0", 24, 512, reference);
        ExpressionTaylorJetRequest request;
        request.program = &pair.runtime;
        request.reference = &reference;
        request.pixelParameter = FormulaParameter::C;
        request.parameterScale = makeScale(-8);
        request.minimumLanding = 1;
        request.maximumCandidateIteration = 12;
        ExpressionTaylorJetResult baseline;
        okay = okay && ExpressionTaylorJetBuilder::build(request, baseline) && baseline.memoryBytes > 1 && baseline.reciprocalCount > 0;
        if (okay) {
            request.memoryLimitBytes = baseline.memoryBytes - 1;
            ExpressionTaylorJetResult constrained;
            okay = !ExpressionTaylorJetBuilder::build(request, constrained) && constrained.status == ExpressionTaylorJetStatus::ResourceLimit && constrained.coefficients.empty() && constrained.coefficientRadii.empty();
        }
        if (!okay) {
            printf("  meromorphic Taylor workspace accounting failed\n");
            ++failures;
        }
    }

    // Finite-series tails, companion metadata, and large derivative/reference
    // factors must reject rather than relying on coefficient decay.
    {
        auto expectTailRejection = [&](const char* name, const char* source, int maximumCompositionOrder) {
            ProgramPair pair;
            ExpressionReferenceOrbitResult reference;
            ExpressionContext fixed;
            bool okay = compilePair(source, FormulaParameter::C, fixed, pair) && buildReference(pair, fixed, FormulaParameter::C, "0", "0", 4, 512, reference);
            ExpressionTaylorJetRequest request;
            request.program = &pair.runtime;
            request.reference = &reference;
            request.pixelParameter = FormulaParameter::C;
            request.parameterScale = makeScale(-3);
            request.minimumLanding = 1;
            request.maximumCandidateIteration = 1;
            request.maximumCompositionOrder = maximumCompositionOrder;
            ExpressionTaylorJetResult jet;
            okay = okay && !ExpressionTaylorJetBuilder::build(request, jet) && jet.status == ExpressionTaylorJetStatus::AccuracyBudget && jet.failureReason.find("tail") != std::string::npos;
            if (!okay) {
                printf("  Taylor certified-tail rejection failed [%s] status=%s reason=%s\n", name, formula::expressionTaylorJetStatusName(jet.status), jet.failureReason.c_str());
                ++failures;
            }
        };
        expectTailRejection("high-exp-derivative", "exp(100*c)-1+c", 24);
        expectTailRejection("large-exp-reference", "exp(c+100)*0+c", 24);
        expectTailRejection("large-sin-companion", "sin(c+complex(0,100))*0+c", 24);
        expectTailRejection("bounded-composition-order", "exp(c)-1+c", 1);
    }
    {
        ProgramPair pair;
        ExpressionReferenceOrbitResult reference;
        ExpressionContext fixed;
        bool okay = compilePair("sin(z)+c", FormulaParameter::C, fixed, pair) && buildReference(pair, fixed, FormulaParameter::C, "0", "0", 4, 512, reference);
        if (okay) {
            for (auto& node : reference.tape) {
                if (node.operation == formula::ExpressionOracleOperation::Sin) {
                    node.flags &= ~formula::OracleTraceHasCompanion;
                    break;
                }
            }
            ExpressionTaylorJetRequest request;
            request.program = &pair.runtime;
            request.reference = &reference;
            request.pixelParameter = FormulaParameter::C;
            request.parameterScale = makeScale(-16);
            request.minimumLanding = 1;
            request.maximumCandidateIteration = 1;
            ExpressionTaylorJetResult jet;
            okay = !ExpressionTaylorJetBuilder::build(request, jet) && jet.status == ExpressionTaylorJetStatus::InvalidTape;
        }
        if (!okay) {
            printf("  Taylor companion guard failed\n");
            ++failures;
        }
    }

    // The certified prefix must stop before a known escape.
    {
        ProgramPair pair;
        ExpressionReferenceOrbitResult reference;
        ExpressionContext fixed;
        bool okay = compilePair("z*z+c", FormulaParameter::C, fixed, pair) && buildReference(pair, fixed, FormulaParameter::C, "2", "0", 10, 512, reference);
        ExpressionTaylorJetRequest request;
        request.program = &pair.runtime;
        request.reference = &reference;
        request.pixelParameter = FormulaParameter::C;
        request.parameterScale = makeScale(-600);
        request.minimumLanding = 1;
        request.maximumCandidateIteration = 8;
        ExpressionTaylorJetResult jet;
        okay = okay && ExpressionTaylorJetBuilder::build(request, jet) && jet.landingIteration == 1;
        if (!okay) {
            printf("  Taylor first-escape stop failed landing=%d\n", jet.landingIteration);
            ++failures;
        }
    }
    {
        ProgramPair pair;
        ExpressionReferenceOrbitResult reference;
        ExpressionContext fixed;
        bool okay = compilePair("exp(c)+c", FormulaParameter::C, fixed, pair) && buildReference(pair, fixed, FormulaParameter::C, "1.5", "0", 4, 512, reference);
        ExpressionTaylorJetRequest request;
        request.program = &pair.runtime;
        request.reference = &reference;
        request.pixelParameter = FormulaParameter::C;
        request.parameterScale = makeScale(-600);
        request.minimumLanding = 1;
        request.maximumCandidateIteration = 2;
        ExpressionTaylorJetResult jet;
        okay = okay && !ExpressionTaylorJetBuilder::build(request, jet) && jet.landingIteration == 0 && (jet.status == ExpressionTaylorJetStatus::BailoutUncertain || jet.status == ExpressionTaylorJetStatus::NoCoverage);
        if (!okay) {
            printf("  entire Taylor first-escape stop failed landing=%d status=%s\n", jet.landingIteration, formula::expressionTaylorJetStatusName(jet.status));
            ++failures;
        }
    }
    {
        ProgramPair pair;
        ExpressionReferenceOrbitResult reference;
        ExpressionContext fixed;
        bool okay = compilePair("conj(z)+3", FormulaParameter::C, fixed, pair) && buildReference(pair, fixed, FormulaParameter::C, "0", "0", 4, 512, reference);
        ExpressionTaylorJetRequest request;
        request.program = &pair.runtime;
        request.reference = &reference;
        request.pixelParameter = FormulaParameter::C;
        request.parameterScale = makeScale(-600);
        request.minimumLanding = 1;
        request.maximumCandidateIteration = 3;
        ExpressionTaylorJetResult jet;
        okay = okay && ExpressionTaylorJetBuilder::build(request, jet) && jet.landingIteration == 1 && jet.layout == formula::ExpressionTaylorJetLayout::RealBivariate;
        if (!okay) {
            printf("  real-bivariate Taylor first-escape stop failed landing=%d status=%s\n", jet.landingIteration, formula::expressionTaylorJetStatusName(jet.status));
            ++failures;
        }
    }
    {
        ProgramPair pair;
        ExpressionReferenceOrbitResult reference;
        ExpressionContext fixed;
        bool okay = compilePair("10*log(c+2)", FormulaParameter::C, fixed, pair) && buildReference(pair, fixed, FormulaParameter::C, "0", "0", 4, 512, reference);
        ExpressionTaylorJetRequest request;
        request.program = &pair.runtime;
        request.reference = &reference;
        request.pixelParameter = FormulaParameter::C;
        request.parameterScale = makeScale(-20);
        request.minimumLanding = 1;
        request.maximumCandidateIteration = 2;
        ExpressionTaylorJetResult jet;
        okay = okay && !ExpressionTaylorJetBuilder::build(request, jet) && jet.landingIteration == 0 && (jet.status == ExpressionTaylorJetStatus::BailoutUncertain || jet.status == ExpressionTaylorJetStatus::NoCoverage);
        if (!okay) {
            printf("  branch Taylor first-escape stop failed landing=%d status=%s\n", jet.landingIteration, formula::expressionTaylorJetStatusName(jet.status));
            ++failures;
        }
    }

    // Renderer parity, thread identity, and a practical repeated speed gate.
    ProgramPair renderPair;
    ExpressionContext renderFixed;
    if (!compilePair("sin(z)+c", FormulaParameter::C, renderFixed, renderPair)) {
        ++failures;
    } else {
        auto makeRenderRequest = [&](int width, int height, int iterations, std::vector<float>& output) {
            output.assign(static_cast<size_t>(width) * height, formula::ExpressionDeepEmptyPixel);
            ExpressionDeepRenderRequest request;
            request.canonicalProgram = &renderPair.canonical;
            request.runtimeProgram = &renderPair.runtime;
            request.center.realDecimal = "0";
            request.center.imaginaryDecimal = "0";
            request.scale.decimal = "1e500";
            request.fixed = renderFixed;
            request.pixelParameter = FormulaParameter::C;
            request.width = width;
            request.height = height;
            request.maxIterations = iterations;
            request.bailout = 4.0;
            request.output = output.data();
            request.outputCount = output.size();
            request.precision.minimumBits = 128;
            request.precision.guardBits = 64;
            request.precision.maximumBits = 4096;
            request.threading.tileWidth = 8;
            request.threading.tileHeight = 5;
            return request;
        };
        std::vector<float> enabled, scalarBatch, disabled, single;
        ExpressionDeepRenderRequest enabledRequest = makeRenderRequest(41, 25, 240, enabled);
        ExpressionDeepRenderResult enabledResult;
        bool okay = formula::renderExpressionDeepFrame(enabledRequest, enabledResult);
        ExpressionDeepRenderRequest scalarBatchRequest = makeRenderRequest(41, 25, 240, scalarBatch);
        scalarBatchRequest.taylor.enableBatchEvaluation = false;
        ExpressionDeepRenderResult scalarBatchResult;
        okay = okay && formula::renderExpressionDeepFrame(scalarBatchRequest, scalarBatchResult);
        ExpressionDeepRenderRequest disabledRequest = makeRenderRequest(41, 25, 240, disabled);
        disabledRequest.taylor.enableTaylor = false;
        ExpressionDeepRenderResult disabledResult;
        okay = okay && formula::renderExpressionDeepFrame(disabledRequest, disabledResult);
        ExpressionDeepRenderRequest singleRequest = makeRenderRequest(41, 25, 240, single);
        singleRequest.threading.threads = 1;
        ExpressionDeepRenderResult singleResult;
        okay = okay && formula::renderExpressionDeepFrame(singleRequest, singleResult);
        if (!okay || enabled != scalarBatch || enabled != disabled || enabled != single || !enabledResult.taylorAccepted || !scalarBatchResult.taylorAccepted || enabledResult.taylorAcceptedPixelCount == 0 || enabledResult.taylorCoveredIterations < enabledRequest.taylor.minimumLanding || enabledResult.taylorAttemptedJetCount != 1 || enabledResult.taylorAcceptedJetCount != 1 || enabledResult.taylorRejectedJetCount != 0 || enabledResult.taylorTileSplitCount != 0 || enabledResult.taylorMaximumTileDepth != 0 || enabledResult.taylorAcceptedPixelCoverage != enabled.size()) {
            printf("  Taylor renderer parity failed accepted=%d landing=%d pixels=%llu/%llu\n", enabledResult.taylorAccepted ? 1 : 0, enabledResult.taylorCoveredIterations, (unsigned long long)enabledResult.taylorAcceptedPixelCount, (unsigned long long)enabledResult.taylorFallbackPixelCount);
            ++failures;
        }

        double jetSeconds[2]{};
        double baselineSeconds[2]{};
        ExpressionDeepRenderResult benchmarkTelemetry;
        bool speedOkay = true;
        for (int repeat = 0; repeat < 2; ++repeat) {
            std::vector<float> jetOutput;
            ExpressionDeepRenderRequest jetRequest = makeRenderRequest(160, 100, 360, jetOutput);
            ExpressionDeepRenderResult jetResult;
            const Clock::time_point jetStart = Clock::now();
            speedOkay = speedOkay && formula::renderExpressionDeepFrame(jetRequest, jetResult);
            jetSeconds[repeat] = std::chrono::duration<double>(Clock::now() - jetStart).count();
            if (repeat == 0) benchmarkTelemetry = jetResult;

            std::vector<float> baselineOutput;
            ExpressionDeepRenderRequest baselineRequest = makeRenderRequest(160, 100, 360, baselineOutput);
            baselineRequest.taylor.enableTaylor = false;
            ExpressionDeepRenderResult baselineResult;
            const Clock::time_point baselineStart = Clock::now();
            speedOkay = speedOkay && formula::renderExpressionDeepFrame(baselineRequest, baselineResult);
            baselineSeconds[repeat] = std::chrono::duration<double>(Clock::now() - baselineStart).count();
            speedOkay = speedOkay && jetOutput == baselineOutput && jetResult.taylorAccepted && jetResult.taylorAcceptedPixelCount > 0 && jetSeconds[repeat] < baselineSeconds[repeat];
        }
        printf("  Taylor 160x100 total jet/no-jet %.3f/%.3f and %.3f/%.3f s\n", jetSeconds[0], baselineSeconds[0], jetSeconds[1], baselineSeconds[1]);

        std::vector<float> jetSmall, mpfrOutput;
        ExpressionDeepRenderRequest jetSmallRequest = makeRenderRequest(64, 40, 360, jetSmall);
        ExpressionDeepRenderResult jetSmallResult;
        const Clock::time_point jetSmallStart = Clock::now();
        bool mpfrOkay = formula::renderExpressionDeepFrame(jetSmallRequest, jetSmallResult);
        const double jetSmallSeconds = std::chrono::duration<double>(Clock::now() - jetSmallStart).count();
        mpfrOutput.assign(jetSmall.size(), formula::ExpressionDeepEmptyPixel);
        ExpressionDeepRenderResult mpfrResult;
        double mpfrSeconds = 0.0;
        if (mpfrOkay) {
            ExpressionDeepRenderRequest mpfrRequest = jetSmallRequest;
            mpfrRequest.forceMpfrFallbackForVerification = true;
            mpfrRequest.output = mpfrOutput.data();
            mpfrRequest.outputCount = mpfrOutput.size();
            const Clock::time_point mpfrStart = Clock::now();
            mpfrOkay = formula::renderExpressionDeepFrame(mpfrRequest, mpfrResult);
            mpfrSeconds = std::chrono::duration<double>(Clock::now() - mpfrStart).count();
        } else {
            mpfrOkay = false;
        }
        printf("  Taylor breakdown reference/build/eval/residual %.3f/%.3f/%.3f/%.3f s; 64x40 jet/MPFR %.3f/%.3f s\n", benchmarkTelemetry.referenceSeconds, benchmarkTelemetry.taylorBuildSeconds, benchmarkTelemetry.taylorEvaluationSeconds, benchmarkTelemetry.taylorResidualSeconds, jetSmallSeconds, mpfrSeconds);
        speedOkay = speedOkay && mpfrOkay && jetSmall == mpfrOutput && jetSmallResult.taylorAccepted && mpfrResult.fastPixelCount == 0 && mpfrResult.fallbackPixelCount == mpfrOutput.size() && jetSmallSeconds < mpfrSeconds;
        if (!speedOkay) {
            printf("  Taylor repeated speed gate failed\n");
            ++failures;
        }
    }

    printf("=== certified expression Taylor jet\n");
    printf("  MPFR containment checks=%d\n", containmentChecks);
    printf("  => %s\n", failures == 0 ? "PASS" : "FAIL");
    return failures == 0 ? 0 : 1;
}

static int runExpressionDeepRenderCase() {
    using formula::ExpressionContext;
    using formula::ExpressionDeepFallbackReason;
    using formula::ExpressionDeepRenderPhase;
    using formula::ExpressionDeepRenderRequest;
    using formula::ExpressionDeepRenderResult;
    using formula::ExpressionDeepRenderStatus;
    using formula::ExpressionDeepVerificationFault;
    using formula::ExpressionError;
    using formula::ExpressionOracle;
    using formula::ExpressionOracleContext;
    using formula::ExpressionProgram;
    using formula::MpfrComplex;

    struct ProgramPair {
        ExpressionProgram canonical;
        ExpressionProgram runtime;
    };
    auto compilePair = [](const char* source, FormulaParameter pixel, const ExpressionContext& fixed, ProgramPair& pair) {
        ExpressionError error;
        return pair.canonical.compile(source, &error) && pair.canonical.specialize(fixed, pixel, pair.runtime, &error);
    };
    auto makeRequest = [](const ProgramPair& pair, const ExpressionContext& fixed, FormulaParameter pixel, const char* centerReal, const char* centerImaginary, int width, int height, int iterations, std::vector<float>& output) {
        output.assign(static_cast<size_t>(width) * height, formula::ExpressionDeepEmptyPixel);
        ExpressionDeepRenderRequest request;
        request.canonicalProgram = &pair.canonical;
        request.runtimeProgram = &pair.runtime;
        request.center.realDecimal = centerReal;
        request.center.imaginaryDecimal = centerImaginary;
        request.scale.decimal = "1e500";
        request.fixed = fixed;
        request.pixelParameter = pixel;
        request.width = width;
        request.height = height;
        request.maxIterations = iterations;
        request.bailout = 4.0;
        request.output = output.data();
        request.outputCount = output.size();
        request.precision.minimumBits = 128;
        request.precision.guardBits = 64;
        request.precision.maximumBits = 4096;
        request.memory.fallbackGuardBits = 128;
        request.threading.tileWidth = 4;
        request.threading.tileHeight = 3;
        return request;
    };
    auto renderOracle = [](const ExpressionProgram& program, const ExpressionContext& fixed, FormulaParameter pixelParameter, const char* centerReal, const char* centerImaginary, const char* scaleText, int width, int height, int maxIterations, double bailout, mpfr_prec_t precision, std::vector<float>& output) {
        output.assign(static_cast<size_t>(width) * height, formula::ExpressionDeepEmptyPixel);
        std::atomic_bool okay{true};
#pragma omp parallel for schedule(dynamic, 1)
        for (int y = 0; y < height; ++y) {
            if (!okay.load(std::memory_order_acquire)) continue;
            ExpressionOracleContext context(precision);
            MpfrComplex center(precision);
            MpfrComplex coordinate(precision);
            MpfrComplex next(precision);
            mpfr_t scale, dxHalf, dyHalf, temporary, magnitude;
            mpfr_inits2(precision, scale, dxHalf, dyHalf, temporary, magnitude, (mpfr_ptr)0);
            bool rowOkay = center.set(centerReal, centerImaginary) && mpfr_set_str(scale, scaleText, 10, MPFR_RNDN) == 0;
            if (rowOkay) {
                mpfr_mul_ui(temporary, scale, static_cast<unsigned long>(width - 1), MPFR_RNDN);
                mpfr_ui_div(dxHalf, 2, temporary, MPFR_RNDN);
                mpfr_mul_ui(temporary, scale, static_cast<unsigned long>(width), MPFR_RNDN);
                mpfr_mul_ui(temporary, temporary, static_cast<unsigned long>(height - 1), MPFR_RNDN);
                mpfr_ui_div(dyHalf, static_cast<unsigned long>(height), temporary, MPFR_RNDN);
                mpfr_mul_ui(dyHalf, dyHalf, 2, MPFR_RNDN);
            }
            for (int x = 0; rowOkay && x < width; ++x) {
                context.z.set(fixed.z.real(), fixed.z.imag());
                context.c.set(fixed.c.real(), fixed.c.imag());
                context.z0.set(fixed.z0.real(), fixed.z0.imag());
                for (size_t p = 0; p < fixed.parameters.size(); ++p) { context.parameters[p].set(fixed.parameters[p].real(), fixed.parameters[p].imag()); }
                mpfr_mul_si(coordinate.re, dxHalf, static_cast<long>(2LL * x - (width - 1LL)), MPFR_RNDN);
                mpfr_add(coordinate.re, coordinate.re, center.re, MPFR_RNDN);
                mpfr_mul_si(coordinate.im, dyHalf, static_cast<long>(2LL * y - (height - 1LL)), MPFR_RNDN);
                mpfr_add(coordinate.im, coordinate.im, center.im, MPFR_RNDN);
                if (pixelParameter == FormulaParameter::C)
                    context.c.set(coordinate);
                else
                    context.z0.set(coordinate);
                context.z.set(context.z0);

                float value = formula::ExpressionDeepInteriorPixel;
                mpfr_hypot(magnitude, context.z.re, context.z.im, MPFR_RNDN);
                bool escaped = mpfr_number_p(magnitude) && mpfr_cmp_d(magnitude, bailout) > 0;
                bool defined = mpfr_number_p(magnitude);
                if (escaped) value = 0.0f;
                for (int iteration = 0; defined && !escaped && iteration < maxIterations; ++iteration) {
                    context.iteration = iteration;
                    std::string error;
                    defined = ExpressionOracle::evaluate(program, context, next, &error) && mpfr_number_p(next.re) && mpfr_number_p(next.im);
                    if (!defined) break;
                    context.z.set(next);
                    mpfr_hypot(magnitude, context.z.re, context.z.im, MPFR_RNDN);
                    defined = mpfr_number_p(magnitude);
                    if (defined && mpfr_cmp_d(magnitude, bailout) > 0) {
                        value = static_cast<float>(iteration + 1);
                        escaped = true;
                    }
                }
                if (!defined) {
                    rowOkay = false;
                    break;
                }
                output[static_cast<size_t>(y) * width + x] = value;
            }
            mpfr_clears(scale, dxHalf, dyHalf, temporary, magnitude, (mpfr_ptr)0);
            if (!rowOkay) okay.store(false, std::memory_order_release);
        }
        return okay.load(std::memory_order_acquire);
    };
    auto reasonCount = [](const ExpressionDeepRenderResult& result, ExpressionDeepFallbackReason reason) { return result.fallbackReasonCounts[static_cast<size_t>(reason)]; };

    int failures = 0;
    int exactFrames = 0;
    bool sawInitialEscape = false;
    bool sawInterior = false;
    bool sawExterior = false;
    auto verifyFrame = [&](const char* name, const char* source, FormulaParameter pixel, const ExpressionContext& fixed, const char* centerReal, const char* centerImaginary, int width, int height, int iterations, uint64_t expectedFallback, ExpressionDeepFallbackReason expectedReason, bool expectTaylorAccepted) {
        ProgramPair pair;
        if (!compilePair(source, pixel, fixed, pair)) {
            printf("  compile failed [%s]\n", name);
            ++failures;
            return;
        }
        std::vector<float> actual;
        ExpressionDeepRenderRequest request = makeRequest(pair, fixed, pixel, centerReal, centerImaginary, width, height, iterations, actual);
        if (expectTaylorAccepted) request.taylor.requirePredictedBenefit = false;
        ExpressionDeepRenderResult result;
        if (!formula::renderExpressionDeepFrame(request, result)) {
            printf("  render failed [%s]: %s/%s\n", name, formula::expressionDeepRenderStatusName(result.status), result.error.c_str());
            ++failures;
            return;
        }
        std::vector<float> expected;
        const bool oracleOkay = renderOracle(pair.runtime, fixed, pixel, centerReal, centerImaginary, "1e500", width, height, iterations, request.bailout, result.fallbackPrecision, expected);
        const uint64_t pixelCount = static_cast<uint64_t>(width) * height;
        const bool expectsReference = pair.runtime.scaledResidualCapability() == formula::ExpressionScaledResidualCapability::ExactCenteredArithmetic || pair.runtime.scaledResidualCapability() == formula::ExpressionScaledResidualCapability::CertifiedEntireCandidate || pair.runtime.scaledResidualCapability() == formula::ExpressionScaledResidualCapability::CertifiedMeromorphicCandidate || pair.runtime.scaledResidualCapability() == formula::ExpressionScaledResidualCapability::CertifiedBranchCandidate || pair.runtime.scaledResidualCapability() == formula::ExpressionScaledResidualCapability::CertifiedRealCandidate || pair.runtime.scaledResidualCapability() == formula::ExpressionScaledResidualCapability::CertifiedPiecewiseCandidate;
        const bool meromorphicCandidate = pair.runtime.scaledResidualCapability() == formula::ExpressionScaledResidualCapability::CertifiedMeromorphicCandidate;
        const bool branchCandidate = pair.runtime.scaledResidualCapability() == formula::ExpressionScaledResidualCapability::CertifiedBranchCandidate;
        const bool realCandidate = pair.runtime.scaledResidualCapability() == formula::ExpressionScaledResidualCapability::CertifiedRealCandidate;
        const bool piecewiseCandidate = pair.runtime.scaledResidualCapability() == formula::ExpressionScaledResidualCapability::CertifiedPiecewiseCandidate;
        const bool hasArg = std::string(source).find("arg(") != std::string::npos;
        const bool hasPolar = std::string(source).find("polar(") != std::string::npos;
        const bool allMpfrTaylorCandidate = (pair.runtime.scaledResidualCapability() == formula::ExpressionScaledResidualCapability::CertifiedEntireCandidate || meromorphicCandidate || branchCandidate || realCandidate || piecewiseCandidate) && expectedFallback == pixelCount;
        const bool uncertainReason = expectedReason == ExpressionDeepFallbackReason::UncertifiedSeries || expectedReason == ExpressionDeepFallbackReason::BranchSensitive || expectedReason == ExpressionDeepFallbackReason::BailoutUncertain || expectedReason == ExpressionDeepFallbackReason::CertificationFailure;
        if (!oracleOkay || actual != expected || result.fallbackPixelCount != expectedFallback || result.fastPixelCount + result.fallbackPixelCount != pixelCount || result.uncertainPixelCount != (uncertainReason ? expectedFallback : 0) || result.undefinedPixelCount != 0 || (expectedFallback != 0 && reasonCount(result, expectedReason) != expectedFallback) || result.maxFallbackReasonCount != expectedFallback || (expectedFallback != 0 && (result.fallbackTileCount == 0 || result.maxTileFallbackRate != 1.0)) || result.selectedPrecision < 1700 || (expectsReference ? (allMpfrTaylorCandidate ? result.referenceBytes != 0 : result.referenceBytes == 0 || result.certificationPrecision <= result.selectedPrecision) : result.referenceBytes != 0) ||
            (expectTaylorAccepted && (!result.taylorAccepted || result.taylorAcceptedPixelCount == 0 || result.taylorCoveredIterations != iterations || (hasArg && (result.taylorArgCompositionCount == 0 || result.taylorBranchCompositionCount == 0 || result.taylorMinimumBranchCutClearance.isZero() || result.taylorMinimumBranchZeroClearance.isZero())) || (hasPolar && (result.taylorPolarCompositionCount == 0 || result.taylorFunctionSeriesCount == 0 || result.taylorMinimumPolarRadiusClearance.isZero())) ||
                                      (meromorphicCandidate                  ? result.taylorReciprocalCount == 0 || result.taylorReciprocalOperationCount == 0 || result.taylorMaximumReciprocalOrder < 8 || result.taylorMinimumDenominatorClearance.isZero()
                                       : branchCandidate                     ? result.taylorBranchCompositionCount == 0 || result.taylorBranchCompositionOperationCount == 0 || result.taylorMaximumBranchSeriesOrder < 1 || result.taylorMinimumBranchCutClearance.isZero() || result.taylorMinimumBranchZeroClearance.isZero()
                                       : realCandidate || piecewiseCandidate ? result.taylorLayout != formula::ExpressionTaylorJetLayout::RealBivariate || result.taylorMonomialCount == 0 || (piecewiseCandidate && (result.taylorAbsBranchCount == 0 || result.taylorMinimumFoldClearance.isZero()))
                                                                             : result.taylorFunctionSeriesCount == 0 || result.taylorFunctionSeriesOperationCount == 0 || result.taylorMaximumFunctionSeriesOrder < 1))) ||
            result.rendererBytes == 0) {
            printf("  exact frame mismatch [%s] oracle/equal=%d/%d fast/fallback=%llu/%llu uncertain/undefined=%llu/%llu maxreason=%llu tiles=%llu rate=%.3f precision=%lld/%lld reasons cert/bailout/exhausted=%llu/%llu/%llu Taylor=%s/%s accepted=%d px=%llu landing=%d layout=%s monomials=%zu ref/renderer=%zu/%zu\n", name, oracleOkay ? 1 : 0, actual == expected ? 1 : 0, (unsigned long long)result.fastPixelCount, (unsigned long long)result.fallbackPixelCount, (unsigned long long)result.uncertainPixelCount, (unsigned long long)result.undefinedPixelCount, (unsigned long long)result.maxFallbackReasonCount, (unsigned long long)result.fallbackTileCount, result.maxTileFallbackRate, (long long)result.selectedPrecision, (long long)result.certificationPrecision, (unsigned long long)reasonCount(result, ExpressionDeepFallbackReason::CertificationFailure), (unsigned long long)reasonCount(result, ExpressionDeepFallbackReason::BailoutUncertain),
                   (unsigned long long)reasonCount(result, ExpressionDeepFallbackReason::ReferenceExhausted), formula::expressionTaylorJetStatusName(result.taylorStatus), result.taylorFailureReason.c_str(), result.taylorAccepted ? 1 : 0, (unsigned long long)result.taylorAcceptedPixelCount, result.taylorCoveredIterations, formula::expressionTaylorJetLayoutName(result.taylorLayout), result.taylorMonomialCount, result.referenceBytes, result.rendererBytes);
            ++failures;
            return;
        }
        for (float value : actual) {
            sawInitialEscape = sawInitialEscape || value == 0.0f;
            sawInterior = sawInterior || value == formula::ExpressionDeepInteriorPixel;
            sawExterior = sawExterior || value > 0.0f;
        }
        ++exactFrames;
    };

    ExpressionContext mandelbrotFixed;
    verifyFrame("arithmetic-c-tail", "z*z+c", FormulaParameter::C, mandelbrotFixed, "-2", "0", 15, 9, 850, 135, ExpressionDeepFallbackReason::BailoutUncertain, false);

    ExpressionContext z0Fixed;
    z0Fixed.c = {};
    verifyFrame("arithmetic-z0-tail", "2*z", FormulaParameter::InitialZ, z0Fixed, "0", "0", 13, 7, 1670, 0, ExpressionDeepFallbackReason::InvalidTape, false);
    verifyFrame("initial-escape", "z+1", FormulaParameter::InitialZ, z0Fixed, "5", "0", 9, 5, 20, 0, ExpressionDeepFallbackReason::InvalidTape, false);
    verifyFrame("finite-interior", "z", FormulaParameter::C, mandelbrotFixed, "3", "0", 7, 5, 12, 0, ExpressionDeepFallbackReason::InvalidTape, false);
    verifyFrame("bailout-gate-fallback", "z", FormulaParameter::InitialZ, z0Fixed, "4", "0", 7, 5, 2, 35, ExpressionDeepFallbackReason::BailoutUncertain, false);

    ExpressionContext transcendentalFixed;
    verifyFrame("sine-taylor", "sin(z)+c", FormulaParameter::C, transcendentalFixed, "0", "0", 7, 5, 20, 0, ExpressionDeepFallbackReason::InvalidTape, true);
    verifyFrame("exp-taylor", "exp(0.1*z)+c", FormulaParameter::C, transcendentalFixed, "-1", "0", 7, 5, 20, 0, ExpressionDeepFallbackReason::InvalidTape, true);
    verifyFrame("composed-sine-taylor", "sin(z*z+c)+c", FormulaParameter::C, transcendentalFixed, "0", "0", 7, 5, 20, 0, ExpressionDeepFallbackReason::InvalidTape, true);
    verifyFrame("mixed-hyperbolic-taylor", "sinh(0.25*z)+cosh(0.25*z)-1+c", FormulaParameter::C, transcendentalFixed, "0", "0", 21, 13, 20, 0, ExpressionDeepFallbackReason::InvalidTape, true);
    verifyFrame("sine-z0-taylor", "sin(z)+c", FormulaParameter::InitialZ, z0Fixed, "0", "0", 7, 5, 20, 0, ExpressionDeepFallbackReason::InvalidTape, true);
    verifyFrame("conjugate-taylor", "conj(z)+c", FormulaParameter::C, transcendentalFixed, "0", "0", 21, 13, 20, 0, ExpressionDeepFallbackReason::InvalidTape, true);
    verifyFrame("real-taylor", "real(z)+c", FormulaParameter::C, transcendentalFixed, "0", "0", 21, 13, 20, 0, ExpressionDeepFallbackReason::InvalidTape, true);
    verifyFrame("imaginary-taylor", "imag(z)*complex(0,1)+c", FormulaParameter::C, transcendentalFixed, "0", "0", 21, 13, 20, 0, ExpressionDeepFallbackReason::InvalidTape, true);
    verifyFrame("norm-taylor", "norm(z)+c", FormulaParameter::C, transcendentalFixed, "0", "0", 21, 13, 20, 0, ExpressionDeepFallbackReason::InvalidTape, true);
    verifyFrame("make-complex-taylor", "complex(real(z),imag(c))+c", FormulaParameter::C, transcendentalFixed, "0", "0", 21, 13, 20, 0, ExpressionDeepFallbackReason::InvalidTape, true);
    verifyFrame("mixed-real-polynomial-taylor", "(z*z+c)*conj(z+c)+c", FormulaParameter::C, transcendentalFixed, "0", "0", 21, 13, 20, 0, ExpressionDeepFallbackReason::InvalidTape, true);
    verifyFrame("conjugate-z0-taylor", "conj(z)+c", FormulaParameter::InitialZ, z0Fixed, "0", "0", 21, 13, 20, 0, ExpressionDeepFallbackReason::InvalidTape, true);
    verifyFrame("arg-taylor", "arg(c+2)", FormulaParameter::C, transcendentalFixed, "0", "0.25", 21, 13, 20, 0, ExpressionDeepFallbackReason::InvalidTape, true);
    verifyFrame("arg-plus-c-taylor", "arg(z+2)+c", FormulaParameter::C, transcendentalFixed, "0", "0.25", 21, 13, 20, 0, ExpressionDeepFallbackReason::InvalidTape, true);
    verifyFrame("polar-taylor", "polar(real(c)+2,imag(z))", FormulaParameter::C, transcendentalFixed, "-1", "0", 21, 13, 20, 0, ExpressionDeepFallbackReason::InvalidTape, true);
    verifyFrame("nested-polar-arg-taylor", "polar(abs(real(c))+1,arg(c+2))", FormulaParameter::C, transcendentalFixed, "0.5", "0.25", 21, 13, 20, 0, ExpressionDeepFallbackReason::InvalidTape, true);
    verifyFrame("arg-first-escape", "arg(c+2)+5", FormulaParameter::C, transcendentalFixed, "0", "0", 7, 5, 20, 35, ExpressionDeepFallbackReason::CertificationFailure, false);
    verifyFrame("polar-first-escape", "polar(real(c)+5,imag(c))", FormulaParameter::C, transcendentalFixed, "0", "0", 7, 5, 20, 35, ExpressionDeepFallbackReason::CertificationFailure, false);
    verifyFrame("unsupported-real-entire", "sin(conj(z))+c", FormulaParameter::C, transcendentalFixed, "0", "0", 7, 5, 20, 35, ExpressionDeepFallbackReason::UnsupportedOperation, false);
    verifyFrame("unsupported-real-branch", "log(conj(z)+2)+c", FormulaParameter::C, transcendentalFixed, "0", "0", 7, 5, 20, 35, ExpressionDeepFallbackReason::BranchSensitive, false);
    verifyFrame("rational-cost-fallback", "z/(c+2)+c", FormulaParameter::C, transcendentalFixed, "0", "0", 7, 5, 20, 35, ExpressionDeepFallbackReason::CertificationFailure, false);
    verifyFrame("divide-taylor", "z/(c+2)+c", FormulaParameter::C, transcendentalFixed, "0", "0", 21, 13, 20, 0, ExpressionDeepFallbackReason::InvalidTape, true);
    verifyFrame("reciprocal-taylor", "1/(z+2)+c", FormulaParameter::C, transcendentalFixed, "0", "0", 21, 13, 20, 0, ExpressionDeepFallbackReason::InvalidTape, true);
    verifyFrame("tangent-taylor", "tan(z)+c", FormulaParameter::C, transcendentalFixed, "0", "0", 21, 13, 20, 0, ExpressionDeepFallbackReason::InvalidTape, true);
    verifyFrame("tanh-taylor", "tanh(z)+c", FormulaParameter::C, transcendentalFixed, "0", "0", 21, 13, 20, 0, ExpressionDeepFallbackReason::InvalidTape, true);
    verifyFrame("nested-meromorphic-taylor", "(sin(z)+c)/(cos(z)+2)", FormulaParameter::C, transcendentalFixed, "0", "0", 21, 13, 20, 0, ExpressionDeepFallbackReason::InvalidTape, true);
    verifyFrame("reciprocal-z0-taylor", "1/(z+2)+c", FormulaParameter::InitialZ, z0Fixed, "0", "0", 21, 13, 20, 0, ExpressionDeepFallbackReason::InvalidTape, true);
    ExpressionContext branchFixed;
    branchFixed.z0 = {1.0, 0.0};
    verifyFrame("log-taylor", "log(z+2)+c", FormulaParameter::C, branchFixed, "0", "0", 7, 5, 12, 0, ExpressionDeepFallbackReason::InvalidTape, true);
    verifyFrame("log10-taylor", "log10(z+2)+c", FormulaParameter::C, branchFixed, "0", "0", 7, 5, 12, 0, ExpressionDeepFallbackReason::InvalidTape, true);
    verifyFrame("sqrt-taylor", "0.25*sqrt(z+2)+c", FormulaParameter::C, branchFixed, "0", "0", 7, 5, 12, 0, ExpressionDeepFallbackReason::InvalidTape, true);
    verifyFrame("power-taylor", "0.1*pow(z+2,c+0.25)+c", FormulaParameter::C, branchFixed, "0", "0", 7, 5, 12, 0, ExpressionDeepFallbackReason::InvalidTape, true);
    verifyFrame("nested-branch-taylor", "0.1*exp(log(z+2))+0.1*sin(sqrt(z+3))/(z+4)+c", FormulaParameter::C, branchFixed, "0", "0", 7, 5, 12, 0, ExpressionDeepFallbackReason::InvalidTape, true);

    {
        auto runPiecewiseParity = [&](const char* name, const char* source, FormulaParameter pixel, const ExpressionContext& fixed, const char* centerReal, const char* centerImaginary, const char* scale, bool requireTaylor, bool requireMixedFold) {
            ProgramPair pair;
            if (!compilePair(source, pixel, fixed, pair) || pair.runtime.scaledResidualCapability() != formula::ExpressionScaledResidualCapability::CertifiedPiecewiseCandidate) {
                printf("  piecewise setup failed [%s]\n", name);
                ++failures;
                return;
            }

            std::vector<float> accelerated;
            ExpressionDeepRenderRequest request = makeRequest(pair, fixed, pixel, centerReal, centerImaginary, 9, 5, 8, accelerated);
            request.scale.decimal = scale;
            request.taylor.minimumLanding = 1;
            request.taylor.requirePredictedBenefit = false;
            request.threading.threads = 1;
            ExpressionDeepRenderResult acceleratedResult;
            bool okay = formula::renderExpressionDeepFrame(request, acceleratedResult);

            std::vector<float> noTaylor(accelerated.size(), formula::ExpressionDeepEmptyPixel);
            ExpressionDeepRenderRequest noTaylorRequest = request;
            noTaylorRequest.taylor.enableTaylor = false;
            noTaylorRequest.output = noTaylor.data();
            noTaylorRequest.outputCount = noTaylor.size();
            ExpressionDeepRenderResult noTaylorResult;
            okay = okay && formula::renderExpressionDeepFrame(noTaylorRequest, noTaylorResult);

            std::vector<float> allMpfr(accelerated.size(), formula::ExpressionDeepEmptyPixel);
            ExpressionDeepRenderRequest allMpfrRequest = noTaylorRequest;
            allMpfrRequest.forceMpfrFallbackForVerification = true;
            allMpfrRequest.output = allMpfr.data();
            allMpfrRequest.outputCount = allMpfr.size();
            ExpressionDeepRenderResult allMpfrResult;
            okay = okay && formula::renderExpressionDeepFrame(allMpfrRequest, allMpfrResult);

            std::vector<float> threaded(accelerated.size(), formula::ExpressionDeepEmptyPixel);
            ExpressionDeepRenderRequest threadedRequest = request;
            threadedRequest.threading.threads = 4;
            threadedRequest.output = threaded.data();
            threadedRequest.outputCount = threaded.size();
            ExpressionDeepRenderResult threadedResult;
            okay = okay && formula::renderExpressionDeepFrame(threadedRequest, threadedResult);

            const uint64_t pixels = static_cast<uint64_t>(accelerated.size());
            okay = okay && accelerated == noTaylor && accelerated == allMpfr && accelerated == threaded && noTaylorResult.fastPixelCount > 0 && allMpfrResult.fastPixelCount == 0 && allMpfrResult.fallbackPixelCount == pixels && acceleratedResult.fastPixelCount == threadedResult.fastPixelCount && acceleratedResult.fallbackPixelCount == threadedResult.fallbackPixelCount;
            if (requireTaylor) { okay = okay && acceleratedResult.taylorAccepted && acceleratedResult.taylorAcceptedPixelCount > 0 && acceleratedResult.taylorAbsBranchCount > 0 && !acceleratedResult.taylorMinimumFoldClearance.isZero(); }
            if (requireMixedFold) { okay = okay && acceleratedResult.taylorAttempted && !acceleratedResult.taylorAccepted && acceleratedResult.taylorFoldRejected && acceleratedResult.taylorFoldRejectionIteration == 0 && acceleratedResult.fastPixelCount > 0 && acceleratedResult.fallbackPixelCount > 0 && acceleratedResult.fallbackPixelCount < pixels && reasonCount(acceleratedResult, ExpressionDeepFallbackReason::BranchSensitive) > 0; }
            if (!okay) {
                printf("  piecewise parity failed [%s] fast/fallback=%llu/%llu noTaylor=%llu/%llu Taylor=%d/%s fold=%d@%d abs=%llu +/-=%llu/%llu\n", name, (unsigned long long)acceleratedResult.fastPixelCount, (unsigned long long)acceleratedResult.fallbackPixelCount, (unsigned long long)noTaylorResult.fastPixelCount, (unsigned long long)noTaylorResult.fallbackPixelCount, acceleratedResult.taylorAccepted ? 1 : 0, formula::expressionTaylorJetStatusName(acceleratedResult.taylorStatus), acceleratedResult.taylorFoldRejected ? 1 : 0, acceleratedResult.taylorFoldRejectionIteration, (unsigned long long)acceleratedResult.taylorAbsBranchCount, (unsigned long long)acceleratedResult.taylorAbsPositiveCellCount, (unsigned long long)acceleratedResult.taylorAbsNegativeCellCount);
                ++failures;
            } else {
                ++exactFrames;
            }
        };

        ExpressionContext piecewiseC;
        piecewiseC.z0 = {0.75, 0.5};
        ExpressionContext piecewiseZ0;
        piecewiseZ0.c = {-0.2, 0.1};
        const struct {
            const char* name;
            const char* source;
        } formulas[] = {{"burning", "sqr(complex(abs(real(z)),abs(imag(z))))+c"}, {"celtic", "complex(abs(real(z*z)),imag(z*z))+c"}, {"buffalo", "complex(abs(real(z*z)),-abs(imag(z*z)))+c"}};
        for (const auto& formula : formulas) {
            std::string name = std::string(formula.name) + "-c-moderate";
            runPiecewiseParity(name.c_str(), formula.source, FormulaParameter::C, piecewiseC, "-0.2", "0.1", "64", false, false);
            name = std::string(formula.name) + "-c-e500";
            runPiecewiseParity(name.c_str(), formula.source, FormulaParameter::C, piecewiseC, "-0.2", "0.1", "1e500", false, false);
            name = std::string(formula.name) + "-z0-moderate";
            runPiecewiseParity(name.c_str(), formula.source, FormulaParameter::InitialZ, piecewiseZ0, "0.75", "0.5", "64", false, false);
            name = std::string(formula.name) + "-z0-e500";
            runPiecewiseParity(name.c_str(), formula.source, FormulaParameter::InitialZ, piecewiseZ0, "0.75", "0.5", "1e500", false, false);
        }
        runPiecewiseParity("abs-real-safe-taylor", "abs(real(c))", FormulaParameter::C, transcendentalFixed, "2", "0", "1e500", true, false);
        runPiecewiseParity("burning-exact-fold", "sqr(complex(abs(real(z)),abs(imag(z))))+c", FormulaParameter::InitialZ, z0Fixed, "0.125", "1", "4", false, true);

        ProgramPair benchmarkPair;
        if (!compilePair("abs(real(c))", FormulaParameter::C, transcendentalFixed, benchmarkPair)) {
            printf("  piecewise benchmark setup failed\n");
            ++failures;
        } else {
            std::array<double, 2> acceleratedSeconds{};
            std::array<double, 2> mpfrSeconds{};
            ExpressionDeepRenderResult telemetry;
            ExpressionDeepRenderResult perStepTelemetry;
            bool benchmarkOkay = true;
            for (size_t repeat = 0; repeat < acceleratedSeconds.size(); ++repeat) {
                std::vector<float> accelerated;
                ExpressionDeepRenderRequest acceleratedRequest = makeRequest(benchmarkPair, transcendentalFixed, FormulaParameter::C, "2", "0", 64, 40, 360, accelerated);
                acceleratedRequest.scale.decimal = "1e500";
                acceleratedRequest.taylor.minimumLanding = 1;
                const Clock::time_point acceleratedStart = Clock::now();
                ExpressionDeepRenderResult acceleratedResult;
                benchmarkOkay = benchmarkOkay && formula::renderExpressionDeepFrame(acceleratedRequest, acceleratedResult);
                acceleratedSeconds[repeat] = std::chrono::duration<double>(Clock::now() - acceleratedStart).count();

                std::vector<float> perStep(accelerated.size(), formula::ExpressionDeepEmptyPixel);
                ExpressionDeepRenderRequest perStepRequest = acceleratedRequest;
                perStepRequest.taylor.enableTaylor = false;
                perStepRequest.output = perStep.data();
                perStepRequest.outputCount = perStep.size();
                ExpressionDeepRenderResult perStepResult;
                benchmarkOkay = benchmarkOkay && formula::renderExpressionDeepFrame(perStepRequest, perStepResult);

                std::vector<float> allMpfr(accelerated.size(), formula::ExpressionDeepEmptyPixel);
                ExpressionDeepRenderRequest mpfrRequest = perStepRequest;
                mpfrRequest.forceMpfrFallbackForVerification = true;
                mpfrRequest.output = allMpfr.data();
                mpfrRequest.outputCount = allMpfr.size();
                const Clock::time_point mpfrStart = Clock::now();
                ExpressionDeepRenderResult mpfrResult;
                benchmarkOkay = benchmarkOkay && formula::renderExpressionDeepFrame(mpfrRequest, mpfrResult);
                mpfrSeconds[repeat] = std::chrono::duration<double>(Clock::now() - mpfrStart).count();
                benchmarkOkay = benchmarkOkay && accelerated == perStep && accelerated == allMpfr;
                if (repeat == 0) {
                    telemetry = acceleratedResult;
                    perStepTelemetry = perStepResult;
                }
            }
            benchmarkOkay = benchmarkOkay && telemetry.preflightAttempted && !telemetry.preflightRejectedFast && telemetry.preflightFallbackCount == 0 && telemetry.preflightFastCount == telemetry.preflightSampleCount;
            const double minimumAccelerated = *std::min_element(acceleratedSeconds.begin(), acceleratedSeconds.end());
            const double maximumAccelerated = *std::max_element(acceleratedSeconds.begin(), acceleratedSeconds.end());
            const double minimumMpfr = *std::min_element(mpfrSeconds.begin(), mpfrSeconds.end());
            const double maximumMpfr = *std::max_element(mpfrSeconds.begin(), mpfrSeconds.end());
            const bool stableTiming = minimumAccelerated > 0.0 && minimumMpfr > 0.0 && maximumAccelerated / minimumAccelerated <= 2.0 && maximumMpfr / minimumMpfr <= 2.0;
            const double minimumSpeedup = maximumAccelerated > 0.0 ? minimumMpfr / maximumAccelerated : 0.0;
            if (!benchmarkOkay || !stableTiming || !(minimumSpeedup > 1.0)) {
                printf("  repeated piecewise e500 benchmark failed stable=%d speedup=%.2fx\n", stableTiming ? 1 : 0, minimumSpeedup);
                ++failures;
            }
            const uint64_t perStepPixels = perStepTelemetry.fastPixelCount + perStepTelemetry.fallbackPixelCount;
            const double foldRate = perStepPixels ? static_cast<double>(perStepTelemetry.fallbackReasonCounts[static_cast<size_t>(ExpressionDeepFallbackReason::BranchSensitive)]) / static_cast<double>(perStepPixels) : 0.0;
            printf("  piecewise e500 repeated Taylor/MPFR %.3f/%.3f and %.3f/%.3f s speedup %.2fx stable=%d coverage=%d accepted=%llu per-step fast/fallback=%llu/%llu fold=%.3f%%\n", acceleratedSeconds[0], mpfrSeconds[0], acceleratedSeconds[1], mpfrSeconds[1], minimumSpeedup, stableTiming ? 1 : 0, telemetry.taylorCoveredIterations, (unsigned long long)telemetry.taylorAcceptedPixelCount, (unsigned long long)perStepTelemetry.fastPixelCount, (unsigned long long)perStepTelemetry.fallbackPixelCount, 100.0 * foldRate);
        }
    }

    // Production regression for the reported componentwise-abs deep view.
    // Preflight may only reject certified work; the final image remains the
    // independently evaluated generic MPFR bytecode result.
    {
        ProgramPair pair;
        const char* source = "sqr(complex(abs(real(z)),abs(imag(z))))+c";
        const char* centerReal = "-1.013951002213813310632862698887121834129";
        const char* centerImaginary = "-0.7988691125646760914741501921763298573252";
        const char* scale = "3.245427860252859436180864346639097915255e12";
        if (!compilePair(source, FormulaParameter::C, transcendentalFixed, pair)) {
            printf("  reported diffabs setup failed\n");
            ++failures;
        } else {
            auto configure = [&](std::vector<float>& output) {
                ExpressionDeepRenderRequest request = makeRequest(pair, transcendentalFixed, FormulaParameter::C, centerReal, centerImaginary, 33, 21, 2000, output);
                request.scale.decimal = scale;
                request.bailout = 100.0;
                return request;
            };
            std::vector<float> production;
            ExpressionDeepRenderRequest productionRequest = configure(production);
            productionRequest.threading.threads = 1;
            ExpressionDeepRenderResult productionResult;
            bool okay = formula::renderExpressionDeepFrame(productionRequest, productionResult);

            std::vector<float> threaded;
            ExpressionDeepRenderRequest threadedRequest = configure(threaded);
            threadedRequest.threading.threads = 4;
            ExpressionDeepRenderResult threadedResult;
            okay = okay && formula::renderExpressionDeepFrame(threadedRequest, threadedResult);

            std::vector<float> genericMpfr;
            ExpressionDeepRenderRequest genericRequest = configure(genericMpfr);
            genericRequest.preflight.enable = false;
            genericRequest.taylor.enableTaylor = false;
            genericRequest.forceMpfrFallbackForVerification = true;
            genericRequest.disableSpecializedPiecewiseMpfrForVerification = true;
            ExpressionDeepRenderResult genericResult;
            okay = okay && formula::renderExpressionDeepFrame(genericRequest, genericResult);

            ProgramPair genericPair;
            std::vector<float> genericProduction;
            ExpressionDeepRenderResult genericProductionResult;
            bool genericProductionOkay = compilePair("c+sqr(complex(abs(real(z)),abs(imag(z))))", FormulaParameter::C, transcendentalFixed, genericPair);
            if (genericProductionOkay) {
                ExpressionDeepRenderRequest request = makeRequest(genericPair, transcendentalFixed, FormulaParameter::C, centerReal, centerImaginary, 33, 21, 2000, genericProduction);
                request.scale.decimal = scale;
                request.bailout = 100.0;
                genericProductionOkay = genericPair.runtime.piecewiseQuadraticKind() == formula::ExpressionPiecewiseQuadraticKind::None && formula::renderExpressionDeepFrame(request, genericProductionResult) && genericProduction == genericMpfr && genericProductionResult.preflightRejectedFast && !genericProductionResult.usedSpecializedPiecewiseMpfr;
            }

            const uint64_t pixels = static_cast<uint64_t>(production.size());
            okay = okay && production == threaded && production == genericMpfr && productionResult.preflightAttempted && productionResult.preflightRejectedFast && productionResult.preflightFallbackCount == productionResult.preflightSampleCount && productionResult.preflightSampleCount == 16 && productionResult.fastPixelCount == 0 && productionResult.fallbackPixelCount == pixels && !productionResult.taylorAttempted && productionResult.usedPiecewiseBigFixed && !productionResult.usedSpecializedPiecewiseMpfr && productionResult.preflightFallbackReasonCounts[static_cast<size_t>(ExpressionDeepFallbackReason::BranchSensitive)] == 16 && productionResult.preflightFirstUncertainHistogram[6] == 16 && productionResult.preflightMinimumFirstUncertainIteration >= productionRequest.preflight.earlyRejectMinimumFirstUncertainIteration && threadedResult.preflightRejectedFast && threadedResult.usedPiecewiseBigFixed && threadedResult.preflightFallbackCount == productionResult.preflightFallbackCount &&
                   threadedResult.preflightFirstUncertainHistogram == productionResult.preflightFirstUncertainHistogram && genericResult.fastPixelCount == 0 && genericResult.fallbackPixelCount == pixels && !genericResult.usedSpecializedPiecewiseMpfr && genericProductionOkay;
            if (!okay) {
                printf("  reported diffabs regression failed production=%llu/%llu preflight=%llu/%llu Taylor=%d specialized=%d bigfixed=%d min-uncertain=%u generic=%llu/%llu\n", (unsigned long long)productionResult.fastPixelCount, (unsigned long long)productionResult.fallbackPixelCount, (unsigned long long)productionResult.preflightFallbackCount, (unsigned long long)productionResult.preflightSampleCount, productionResult.taylorAttempted ? 1 : 0, productionResult.usedSpecializedPiecewiseMpfr ? 1 : 0, productionResult.usedPiecewiseBigFixed ? 1 : 0, productionResult.preflightMinimumFirstUncertainIteration, (unsigned long long)genericResult.fastPixelCount, (unsigned long long)genericResult.fallbackPixelCount);
                ++failures;
            } else {
                ++exactFrames;
            }

            std::vector<float> cancelledOutput;
            ExpressionDeepRenderRequest cancelRequest = configure(cancelledOutput);
            std::atomic<bool> inPreflight{false};
            std::atomic<int> preflightPolls{0};
            cancelRequest.progress = [&](ExpressionDeepRenderPhase phase, uint64_t, uint64_t) { inPreflight.store(phase == ExpressionDeepRenderPhase::Preflight, std::memory_order_release); };
            cancelRequest.shouldCancel = [&] { return inPreflight.load(std::memory_order_acquire) && preflightPolls.fetch_add(1, std::memory_order_relaxed) > 4; };
            ExpressionDeepRenderResult cancelResult;
            const bool cancelOkay = formula::renderExpressionDeepFrame(cancelRequest, cancelResult);
            if (cancelOkay || !cancelResult.cancelled || std::count(cancelledOutput.begin(), cancelledOutput.end(), formula::ExpressionDeepEmptyPixel) != static_cast<ptrdiff_t>(cancelledOutput.size())) {
                printf("  reported diffabs preflight cancellation failed\n");
                ++failures;
            }

            const struct {
                const char* centerRe;
                const char* centerIm;
                const char* zoom;
                int iterations;
                double escape;
            } nearby[] = {{centerReal, centerImaginary, "9.99e11", 500, 100.0}, {centerReal, centerImaginary, "1.001e12", 500, 100.0}, {"-1.0139510022135", "-0.7988691125650", scale, 900, 4.0}, {"-0.9057323612582908890917025099695845505893", "-0.8657625017293160215313785051941801184211", "3.270573007328246235696425358963614125045e12", 2000, 100.0}};
            for (const auto& view : nearby) {
                std::vector<float> actual, expected;
                ExpressionDeepRenderRequest actualRequest = makeRequest(pair, transcendentalFixed, FormulaParameter::C, view.centerRe, view.centerIm, 9, 7, view.iterations, actual);
                actualRequest.scale.decimal = view.zoom;
                actualRequest.bailout = view.escape;
                ExpressionDeepRenderResult actualResult;
                ExpressionDeepRenderRequest expectedRequest = actualRequest;
                expected.assign(actual.size(), formula::ExpressionDeepEmptyPixel);
                expectedRequest.output = expected.data();
                expectedRequest.outputCount = expected.size();
                expectedRequest.preflight.enable = false;
                expectedRequest.taylor.enableTaylor = false;
                expectedRequest.forceMpfrFallbackForVerification = true;
                expectedRequest.disableSpecializedPiecewiseMpfrForVerification = true;
                ExpressionDeepRenderResult expectedResult;
                if (!formula::renderExpressionDeepFrame(actualRequest, actualResult) || !formula::renderExpressionDeepFrame(expectedRequest, expectedResult) || actual != expected) {
                    printf("  reported diffabs nearby parity failed zoom=%s\n", view.zoom);
                    ++failures;
                }
            }

            for (double escape : {2.1, 4194304.0}) {
                std::vector<float> actual, expected;
                ExpressionDeepRenderRequest actualRequest = makeRequest(pair, transcendentalFixed, FormulaParameter::C, centerReal, centerImaginary, 9, 7, 64, actual);
                actualRequest.scale.decimal = scale;
                actualRequest.bailout = escape;
                actualRequest.forceMpfrFallbackForVerification = true;
                ExpressionDeepRenderResult actualResult;
                ExpressionDeepRenderRequest expectedRequest = actualRequest;
                expected.assign(actual.size(), formula::ExpressionDeepEmptyPixel);
                expectedRequest.output = expected.data();
                expectedRequest.outputCount = expected.size();
                expectedRequest.preflight.enable = false;
                expectedRequest.taylor.enableTaylor = false;
                expectedRequest.forceMpfrFallbackForVerification = true;
                expectedRequest.disableSpecializedPiecewiseMpfrForVerification = true;
                ExpressionDeepRenderResult expectedResult;
                if (!formula::renderExpressionDeepFrame(actualRequest, actualResult) || !formula::renderExpressionDeepFrame(expectedRequest, expectedResult) || actual != expected || actualResult.usedPiecewiseBigFixed || !actualResult.usedSpecializedPiecewiseMpfr) {
                    printf("  BigFixed piecewise range gate failed bailout=%.17g status=%s/%s bigfixed=%d specialized=%d equal=%d\n", escape, formula::expressionDeepRenderStatusName(actualResult.status), formula::expressionDeepRenderStatusName(expectedResult.status), actualResult.usedPiecewiseBigFixed ? 1 : 0, actualResult.usedSpecializedPiecewiseMpfr ? 1 : 0, actual == expected ? 1 : 0);
                    ++failures;
                }
            }

            {
                std::vector<float> actual, expected;
                ExpressionDeepRenderRequest actualRequest = makeRequest(pair, transcendentalFixed, FormulaParameter::C, centerReal, centerImaginary, 9, 7, 5000, actual);
                actualRequest.scale.decimal = scale;
                actualRequest.bailout = 100.0;
                actualRequest.precision.minimumBits = 53;
                actualRequest.precision.guardBits = 0;
                actualRequest.memory.fallbackGuardBits = 0;
                actualRequest.forceMpfrFallbackForVerification = true;
                ExpressionDeepRenderResult actualResult;
                ExpressionDeepRenderRequest expectedRequest = actualRequest;
                expected.assign(actual.size(), formula::ExpressionDeepEmptyPixel);
                expectedRequest.output = expected.data();
                expectedRequest.outputCount = expected.size();
                expectedRequest.disableSpecializedPiecewiseMpfrForVerification = true;
                ExpressionDeepRenderResult expectedResult;
                if (!formula::renderExpressionDeepFrame(actualRequest, actualResult) || !formula::renderExpressionDeepFrame(expectedRequest, expectedResult) || actual != expected || actualResult.usedPiecewiseBigFixed || !actualResult.usedSpecializedPiecewiseMpfr) {
                    printf("  BigFixed piecewise low-precision gate failed\n");
                    ++failures;
                }
            }

            {
                ExpressionContext tinyInitial = transcendentalFixed;
                tinyInitial.z0 = {1e-150, 1e-150};
                std::vector<float> actual, expected;
                ExpressionDeepRenderRequest actualRequest = makeRequest(pair, tinyInitial, FormulaParameter::C, "-2", "0", 5, 3, 500, actual);
                actualRequest.scale.decimal = "1e12";
                actualRequest.bailout = 2.0;
                actualRequest.forceMpfrFallbackForVerification = true;
                ExpressionDeepRenderResult actualResult;
                ExpressionDeepRenderRequest expectedRequest = actualRequest;
                expected.assign(actual.size(), formula::ExpressionDeepEmptyPixel);
                expectedRequest.output = expected.data();
                expectedRequest.outputCount = expected.size();
                expectedRequest.disableSpecializedPiecewiseMpfrForVerification = true;
                ExpressionDeepRenderResult expectedResult;
                if (!formula::renderExpressionDeepFrame(actualRequest, actualResult) || !formula::renderExpressionDeepFrame(expectedRequest, expectedResult) || actual != expected || actualResult.usedPiecewiseBigFixed || !actualResult.usedSpecializedPiecewiseMpfr) {
                    printf("  BigFixed piecewise underflow fallback failed\n");
                    ++failures;
                }
            }

            {
                ExpressionContext boundaryInitial = transcendentalFixed;
                boundaryInitial.z0 = {2.0, std::ldexp(1.0, -160)};
                std::vector<float> actual, expected;
                ExpressionDeepRenderRequest actualRequest = makeRequest(pair, boundaryInitial, FormulaParameter::C, "0", "0", 3, 3, 4, actual);
                actualRequest.scale.decimal = "1e12";
                actualRequest.bailout = 2.0;
                actualRequest.forceMpfrFallbackForVerification = true;
                ExpressionDeepRenderResult actualResult;
                ExpressionDeepRenderRequest expectedRequest = actualRequest;
                expected.assign(actual.size(), formula::ExpressionDeepEmptyPixel);
                expectedRequest.output = expected.data();
                expectedRequest.outputCount = expected.size();
                expectedRequest.disableSpecializedPiecewiseMpfrForVerification = true;
                ExpressionDeepRenderResult expectedResult;
                if (!formula::renderExpressionDeepFrame(actualRequest, actualResult) || !formula::renderExpressionDeepFrame(expectedRequest, expectedResult) || actual != expected || actualResult.specializedPiecewiseMpfrPixelCount == 0) {
                    printf("  BigFixed piecewise bailout-boundary fallback failed\n");
                    ++failures;
                }
            }
        }
    }

    // A fold or principal-cut stripe must only exclude ambiguous leaf tiles.
    {
        auto runAdaptiveStripe = [&](const char* name, const char* source, FormulaParameter pixel, const char* centerReal, const char* centerImaginary, const char* scale, bool expectFold, bool expectCut, bool expectPole) {
            ProgramPair pair;
            if (!compilePair(source, pixel, transcendentalFixed, pair)) {
                printf("  adaptive stripe setup failed [%s]\n", name);
                ++failures;
                return;
            }
            auto configure = [&](std::vector<float>& output) {
                ExpressionDeepRenderRequest request = makeRequest(pair, transcendentalFixed, pixel, centerReal, centerImaginary, 65, 33, 24, output);
                request.scale.decimal = scale;
                request.taylor.minimumLanding = 1;
                request.taylor.requirePredictedBenefit = false;
                request.taylor.maximumDepth = expectPole ? 9 : 11;
                request.taylor.minimumTileWidth = 1;
                request.taylor.minimumTileHeight = 1;
                request.taylor.maximumJetCount = expectPole ? 512 : 1024;
                request.taylor.maximumRejectedBeforeFirstAcceptance = request.taylor.maximumJetCount;
                return request;
            };

            std::vector<float> adaptive;
            ExpressionDeepRenderRequest adaptiveRequest = configure(adaptive);
            ExpressionDeepRenderResult adaptiveResult;
            bool okay = formula::renderExpressionDeepFrame(adaptiveRequest, adaptiveResult);

            std::vector<float> wholeFrame;
            ExpressionDeepRenderRequest wholeRequest = configure(wholeFrame);
            wholeRequest.taylor.enableTileTaylor = false;
            ExpressionDeepRenderResult wholeResult;
            okay = okay && formula::renderExpressionDeepFrame(wholeRequest, wholeResult);

            std::vector<float> noTaylor;
            ExpressionDeepRenderRequest noTaylorRequest = configure(noTaylor);
            noTaylorRequest.taylor.enableTaylor = false;
            ExpressionDeepRenderResult noTaylorResult;
            okay = okay && formula::renderExpressionDeepFrame(noTaylorRequest, noTaylorResult);

            std::vector<float> allMpfr;
            ExpressionDeepRenderRequest mpfrRequest = configure(allMpfr);
            mpfrRequest.forceMpfrFallbackForVerification = true;
            ExpressionDeepRenderResult mpfrResult;
            okay = okay && formula::renderExpressionDeepFrame(mpfrRequest, mpfrResult);

            std::vector<float> threaded;
            ExpressionDeepRenderRequest threadedRequest = configure(threaded);
            threadedRequest.threading.threads = 4;
            ExpressionDeepRenderResult threadedResult;
            okay = okay && formula::renderExpressionDeepFrame(threadedRequest, threadedResult);

            const uint64_t pixels = static_cast<uint64_t>(adaptive.size());
            okay = okay && adaptive == wholeFrame && adaptive == noTaylor && adaptive == allMpfr && adaptive == threaded && !wholeResult.taylorAccepted && adaptiveResult.taylorAccepted && adaptiveResult.taylorAttemptedJetCount > 1 && adaptiveResult.taylorAcceptedJetCount > 0 && adaptiveResult.taylorRejectedJetCount > 0 && adaptiveResult.taylorTileSplitCount > 0 && adaptiveResult.taylorMaximumTileDepth > 0 && adaptiveResult.taylorAcceptedPixelCoverage > pixels / 2 && adaptiveResult.taylorAcceptedPixelCoverage < pixels && adaptiveResult.taylorWeightedLanding >= 1.0 && adaptiveResult.taylorFallbackPixelCount > 0 && adaptiveResult.taylorRetainedBytes > 0 && adaptiveResult.taylorTileMapHash != 0 && adaptiveResult.taylorTileMapHash == threadedResult.taylorTileMapHash && adaptiveResult.taylorAcceptedPixelCoverage == threadedResult.taylorAcceptedPixelCoverage && adaptiveResult.taylorAttemptedJetCount == threadedResult.taylorAttemptedJetCount &&
                   (!expectFold || adaptiveResult.taylorFoldRejectedJetCount > 0 && adaptiveResult.fallbackPixelCount == 32) && (!expectCut || adaptiveResult.taylorCutRejectedJetCount > 0 && adaptiveResult.fallbackPixelCount == 65) && (!expectPole || adaptiveResult.taylorPoleRejectedJetCount > 0);

            std::vector<float> countLimited;
            ExpressionDeepRenderRequest countRequest = configure(countLimited);
            countRequest.taylor.maximumJetCount = 2;
            countRequest.taylor.maximumRejectedBeforeFirstAcceptance = 2;
            ExpressionDeepRenderResult countResult;
            okay = okay && formula::renderExpressionDeepFrame(countRequest, countResult) && countLimited == adaptive && countResult.taylorAttemptedJetCount <= 2;

            std::vector<float> memoryLimited;
            ExpressionDeepRenderRequest memoryRequest = configure(memoryLimited);
            memoryRequest.taylor.maximumJetMemoryBytes = 1;
            ExpressionDeepRenderResult memoryResult;
            okay = okay && formula::renderExpressionDeepFrame(memoryRequest, memoryResult) && memoryLimited == adaptive && !memoryResult.taylorAccepted;

            if (!okay) {
                printf("  adaptive stripe failed [%s] status=%s/%s threaded=%s/%s whole=%d adaptive=%d attempts/accepted/rejected=%llu/%llu/%llu splits=%llu depth=%d coverage=%llu/%llu landing=%.2f runtime=%llu fallback=%llu fold/cut/pole=%llu/%llu/%llu hash=%llu/%llu limits=%llu/%d\n", name, formula::expressionDeepRenderStatusName(adaptiveResult.status), adaptiveResult.error.c_str(), formula::expressionDeepRenderStatusName(threadedResult.status), threadedResult.error.c_str(), wholeResult.taylorAccepted ? 1 : 0, adaptiveResult.taylorAccepted ? 1 : 0, (unsigned long long)adaptiveResult.taylorAttemptedJetCount, (unsigned long long)adaptiveResult.taylorAcceptedJetCount, (unsigned long long)adaptiveResult.taylorRejectedJetCount, (unsigned long long)adaptiveResult.taylorTileSplitCount, adaptiveResult.taylorMaximumTileDepth, (unsigned long long)adaptiveResult.taylorAcceptedPixelCoverage, (unsigned long long)pixels, adaptiveResult.taylorWeightedLanding,
                       (unsigned long long)adaptiveResult.taylorAcceptedPixelCount, (unsigned long long)adaptiveResult.fallbackPixelCount, (unsigned long long)adaptiveResult.taylorFoldRejectedJetCount, (unsigned long long)adaptiveResult.taylorCutRejectedJetCount, (unsigned long long)adaptiveResult.taylorPoleRejectedJetCount, (unsigned long long)adaptiveResult.taylorTileMapHash, (unsigned long long)threadedResult.taylorTileMapHash, (unsigned long long)countResult.taylorAttemptedJetCount, memoryResult.taylorAccepted ? 1 : 0);
                ++failures;
            } else {
                ++exactFrames;
            }
        };

        runAdaptiveStripe("abs-z0-cross-fold", "abs(real(z0))", FormulaParameter::InitialZ, "0", "0", "8", true, false, false);
        runAdaptiveStripe("log-cross-cut", "0.001*log(c)", FormulaParameter::C, "-2", "0", "10", false, true, false);
        runAdaptiveStripe("divide-near-pole", "0.001/(c+2)", FormulaParameter::C, "-1.901", "0", "10", false, false, true);

        ProgramPair benchmarkPair;
        if (!compilePair("0.001*log(c)", FormulaParameter::C, transcendentalFixed, benchmarkPair)) {
            ++failures;
        } else {
            std::array<double, 2> adaptiveSeconds{};
            std::array<double, 2> wholeFrameSeconds{};
            ExpressionDeepRenderResult benchmarkTelemetry;
            bool benchmarkOkay = true;
            for (size_t repeat = 0; repeat < adaptiveSeconds.size(); ++repeat) {
                std::vector<float> adaptive;
                ExpressionDeepRenderRequest adaptiveRequest = makeRequest(benchmarkPair, transcendentalFixed, FormulaParameter::C, "-2", "0", 64, 40, 120, adaptive);
                adaptiveRequest.taylor.maximumDepth = 5;
                adaptiveRequest.taylor.maximumJetCount = 64;
                ExpressionDeepRenderResult adaptiveResult;
                const Clock::time_point adaptiveStart = Clock::now();
                const bool adaptiveOkay = formula::renderExpressionDeepFrame(adaptiveRequest, adaptiveResult);
                adaptiveSeconds[repeat] = std::chrono::duration<double>(Clock::now() - adaptiveStart).count();

                std::vector<float> wholeFrame;
                ExpressionDeepRenderRequest wholeRequest = makeRequest(benchmarkPair, transcendentalFixed, FormulaParameter::C, "-2", "0", 64, 40, 120, wholeFrame);
                wholeRequest.taylor.enableTileTaylor = false;
                ExpressionDeepRenderResult wholeResult;
                const Clock::time_point wholeStart = Clock::now();
                const bool wholeOkay = formula::renderExpressionDeepFrame(wholeRequest, wholeResult);
                wholeFrameSeconds[repeat] = std::chrono::duration<double>(Clock::now() - wholeStart).count();

                const uint64_t pixels = static_cast<uint64_t>(adaptive.size());
                const bool accepted = adaptiveResult.taylorAccepted;
                benchmarkOkay = benchmarkOkay && adaptiveOkay && wholeOkay && adaptive == wholeFrame && !wholeResult.taylorAccepted && wholeResult.fallbackPixelCount == pixels && (accepted ? adaptiveResult.taylorAcceptedJetCount > 1 && adaptiveResult.taylorAcceptedPixelCoverage > 0 && adaptiveResult.fallbackPixelCount < pixels && adaptiveSeconds[repeat] < wholeFrameSeconds[repeat] : adaptiveResult.taylorAcceptedJetCount == 0 && adaptiveResult.taylorAcceptedPixelCoverage == 0 && adaptiveSeconds[repeat] <= wholeFrameSeconds[repeat] * 1.10);
                if (repeat == 0) benchmarkTelemetry = adaptiveResult;
            }
            printf("  adaptive cross-cut e500 Taylor/whole-frame %.3f/%.3f and %.3f/%.3f s coverage=%llu jets=%llu/%llu fallback=%llu\n", adaptiveSeconds[0], wholeFrameSeconds[0], adaptiveSeconds[1], wholeFrameSeconds[1], (unsigned long long)benchmarkTelemetry.taylorAcceptedPixelCoverage, (unsigned long long)benchmarkTelemetry.taylorAcceptedJetCount, (unsigned long long)benchmarkTelemetry.taylorAttemptedJetCount, (unsigned long long)benchmarkTelemetry.fallbackPixelCount);
            if (!benchmarkOkay) {
                printf("  adaptive cross-cut e500 speed gate failed\n");
                ++failures;
            }
        }
    }

    // Pole clearance is a whole-frame proof, not point-only tape metadata.
    // A frame whose denominator neighborhood crosses zero must release the
    // Taylor/reference resources and render entirely through MPFR; a nearby
    // separated frame must be accepted and remain byte-identical.
    {
        ProgramPair pair;
        const char* source = "0.001/(c+2)";
        if (!compilePair(source, FormulaParameter::C, transcendentalFixed, pair)) {
            ++failures;
        } else {
            auto runPoleFrame = [&](const char* scale, bool expectAccepted) {
                std::vector<float> actual, expected;
                ExpressionDeepRenderRequest request = makeRequest(pair, transcendentalFixed, FormulaParameter::C, "-1.9", "0", 64, 40, 20, actual);
                request.scale.decimal = scale;
                request.taylor.enableTileTaylor = false;
                ExpressionDeepRenderResult result;
                bool okay = formula::renderExpressionDeepFrame(request, result) && renderOracle(pair.runtime, transcendentalFixed, FormulaParameter::C, "-1.9", "0", scale, 64, 40, 20, 4.0, result.fallbackPrecision, expected) && actual == expected && result.taylorAttempted;
                if (expectAccepted) {
                    okay = okay && result.taylorAccepted && result.taylorAcceptedPixelCount > 0 && result.fallbackPixelCount == 0 && result.taylorReciprocalCount > 0 && !result.taylorMinimumDenominatorClearance.isZero();
                } else {
                    const uint64_t pixels = static_cast<uint64_t>(actual.size());
                    okay = okay && !result.taylorAccepted && result.taylorPoleRejected && result.taylorStatus == formula::ExpressionTaylorJetStatus::PoleRejected && result.fastPixelCount == 0 && result.fallbackPixelCount == pixels && result.referenceBytes == 0 && reasonCount(result, ExpressionDeepFallbackReason::CertificationFailure) == pixels;
                }
                if (!okay) {
                    printf("  pole-clearance frame failed scale=%s accepted=%d status=%s reason=%s clearance=(%.17g,e%lld)\n", scale, result.taylorAccepted ? 1 : 0, formula::expressionTaylorJetStatusName(result.taylorStatus), result.taylorFailureReason.c_str(), result.taylorMinimumDenominatorClearance.mantissa, (long long)result.taylorMinimumDenominatorClearance.exponent);
                    ++failures;
                }
            };
            runPoleFrame("10", false);
            runPoleFrame("320", true);
        }
    }

    // Principal-cut clearance is also a whole-frame proof. Rejected frames
    // release every fast resource and use exact MPFR for every pixel.
    {
        ProgramPair pair;
        const char* source = "0.001*log(c)";
        if (!compilePair(source, FormulaParameter::C, transcendentalFixed, pair)) {
            ++failures;
        } else {
            auto runBranchFrame = [&](const char* name, const char* centerImaginary, const char* scale, bool expectAccepted) {
                std::vector<float> actual, expected;
                ExpressionDeepRenderRequest request = makeRequest(pair, transcendentalFixed, FormulaParameter::C, "-2", centerImaginary, 64, 40, 20, actual);
                request.scale.decimal = scale;
                request.taylor.enableTileTaylor = false;
                ExpressionDeepRenderResult result;
                bool okay = formula::renderExpressionDeepFrame(request, result) && renderOracle(pair.runtime, transcendentalFixed, FormulaParameter::C, "-2", centerImaginary, scale, 64, 40, 20, 4.0, result.fallbackPrecision, expected) && actual == expected && result.taylorAttempted;
                if (expectAccepted) {
                    okay = okay && result.taylorAccepted && result.taylorAcceptedPixelCount > 0 && result.fallbackPixelCount == 0 && result.taylorBranchCompositionCount > 0 && !result.taylorMinimumBranchCutClearance.isZero() && !result.taylorMinimumBranchZeroClearance.isZero();
                } else {
                    const uint64_t pixels = static_cast<uint64_t>(actual.size());
                    okay = okay && !result.taylorAccepted && result.taylorBranchRejected && result.taylorStatus == formula::ExpressionTaylorJetStatus::BranchRejected && result.fastPixelCount == 0 && result.fallbackPixelCount == pixels && result.referenceBytes == 0 && reasonCount(result, ExpressionDeepFallbackReason::CertificationFailure) == pixels;
                }
                if (!okay) {
                    printf("  branch-clearance frame failed [%s] accepted=%d status=%s reason=%s cut=(%.17g,e%lld) zero=(%.17g,e%lld)\n", name, result.taylorAccepted ? 1 : 0, formula::expressionTaylorJetStatusName(result.taylorStatus), result.taylorFailureReason.c_str(), result.taylorMinimumBranchCutClearance.mantissa, (long long)result.taylorMinimumBranchCutClearance.exponent, result.taylorMinimumBranchZeroClearance.mantissa, (long long)result.taylorMinimumBranchZeroClearance.exponent);
                    ++failures;
                }
            };
            runBranchFrame("near-cut-overlap", "0.1", "10", false);
            runBranchFrame("near-cut-safe", "0.1", "320", true);
            runBranchFrame("upper-lip", "0", "1e500", false);
            runBranchFrame("lower-lip", "-0", "1e500", false);
        }
    }

    // Arg and Polar use the same all-or-MPFR frame policy: a failed whole-
    // frame proof releases the reference and every fast-path allocation.
    {
        ProgramPair argPair;
        const char* source = "0.001*arg(c)";
        if (!compilePair(source, FormulaParameter::C, transcendentalFixed, argPair)) {
            ++failures;
        } else {
            auto runArgFrame = [&](const char* name, const char* centerImaginary, const char* scale, bool expectAccepted) {
                std::vector<float> actual, expected;
                ExpressionDeepRenderRequest request = makeRequest(argPair, transcendentalFixed, FormulaParameter::C, "-2", centerImaginary, 32, 20, 20, actual);
                request.scale.decimal = scale;
                request.taylor.requirePredictedBenefit = false;
                request.taylor.enableTileTaylor = false;
                ExpressionDeepRenderResult result;
                bool okay = formula::renderExpressionDeepFrame(request, result) && renderOracle(argPair.runtime, transcendentalFixed, FormulaParameter::C, "-2", centerImaginary, scale, 32, 20, 20, 4.0, result.fallbackPrecision, expected) && actual == expected && result.taylorAttempted;
                if (expectAccepted) {
                    okay = okay && result.taylorAccepted && result.fallbackPixelCount == 0 && result.taylorArgCompositionCount > 0 && result.taylorBranchCompositionCount > 0 && !result.taylorMinimumBranchCutClearance.isZero() && !result.taylorMinimumBranchZeroClearance.isZero();
                } else {
                    const uint64_t pixels = static_cast<uint64_t>(actual.size());
                    okay = okay && !result.taylorAccepted && result.taylorBranchRejected && !result.taylorArgRejectionReason.empty() && result.fastPixelCount == 0 && result.fallbackPixelCount == pixels && result.referenceBytes == 0 && reasonCount(result, ExpressionDeepFallbackReason::CertificationFailure) == pixels;
                }
                if (!okay) {
                    printf("  arg frame policy failed [%s] accepted=%d status=%s reason=%s cut=(%.6g,e%lld) zero=(%.6g,e%lld)\n", name, result.taylorAccepted ? 1 : 0, formula::expressionTaylorJetStatusName(result.taylorStatus), result.taylorFailureReason.c_str(), result.taylorMinimumBranchCutClearance.mantissa, (long long)result.taylorMinimumBranchCutClearance.exponent, result.taylorMinimumBranchZeroClearance.mantissa, (long long)result.taylorMinimumBranchZeroClearance.exponent);
                    ++failures;
                }
            };
            runArgFrame("near-cut-overlap", "0.02", "10", false);
            runArgFrame("near-cut-safe", "0.1", "320", true);
            runArgFrame("upper-lip", "0", "1e500", false);
            runArgFrame("lower-lip", "-0", "1e500", false);
        }

        auto runPolarFrame = [&](const char* name, const char* source, const char* scale, bool expectAccepted, bool exactZeroRadius) {
            ProgramPair pair;
            bool setup = compilePair(source, FormulaParameter::C, transcendentalFixed, pair);
            std::vector<float> actual, expected;
            ExpressionDeepRenderRequest request = makeRequest(pair, transcendentalFixed, FormulaParameter::C, "0", "0", 32, 20, 20, actual);
            request.scale.decimal = scale;
            request.taylor.requirePredictedBenefit = false;
            request.taylor.enableTileTaylor = false;
            ExpressionDeepRenderResult result;
            bool okay = setup && formula::renderExpressionDeepFrame(request, result) && renderOracle(pair.runtime, transcendentalFixed, FormulaParameter::C, "0", "0", scale, 32, 20, 20, 4.0, result.fallbackPrecision, expected) && actual == expected && result.taylorAttempted;
            if (expectAccepted) {
                okay = okay && result.taylorAccepted && result.fallbackPixelCount == 0 && result.taylorPolarCompositionCount > 0 && !result.taylorPolarRejected && (exactZeroRadius ? result.taylorMinimumPolarRadiusClearance.isZero() : !result.taylorMinimumPolarRadiusClearance.isZero() && result.taylorFunctionSeriesCount >= 2);
            } else {
                const uint64_t pixels = static_cast<uint64_t>(actual.size());
                okay = okay && !result.taylorAccepted && result.taylorPolarRejected && !result.taylorPolarRejectionReason.empty() && result.fastPixelCount == 0 && result.fallbackPixelCount == pixels && result.referenceBytes == 0 && reasonCount(result, ExpressionDeepFallbackReason::CertificationFailure) == pixels;
            }
            if (!okay) {
                printf("  polar frame policy failed [%s] accepted=%d status=%s reason=%s radius=(%.6g,e%lld) count=%llu rejected=%d reject-reason=%s fast/fallback=%llu/%llu ref=%zu cert-reason=%llu oracle/equal=%d/%d\n", name, result.taylorAccepted ? 1 : 0, formula::expressionTaylorJetStatusName(result.taylorStatus), result.taylorFailureReason.c_str(), result.taylorMinimumPolarRadiusClearance.mantissa, (long long)result.taylorMinimumPolarRadiusClearance.exponent, (unsigned long long)result.taylorPolarCompositionCount, result.taylorPolarRejected ? 1 : 0, result.taylorPolarRejectionReason.c_str(), (unsigned long long)result.fastPixelCount, (unsigned long long)result.fallbackPixelCount, result.referenceBytes, (unsigned long long)reasonCount(result, ExpressionDeepFallbackReason::CertificationFailure), expected.empty() ? 0 : 1, actual == expected ? 1 : 0);
                ++failures;
            }
        };
        runPolarFrame("radius-crosses-negative", "polar(real(c)+0.01,imag(c))", "250", false, false);
        runPolarFrame("near-zero-safe", "polar(real(c)+0.01,imag(c))", "1000", true, false);
        runPolarFrame("exact-zero", "polar(0,real(c))", "1e500", true, true);
    }

    // Independently rebuild the higher-precision point orbit and verify that
    // every retained compact arithmetic value lies inside its stored radius.
    {
        ProgramPair pair;
        if (!compilePair("z*z+c", FormulaParameter::C, mandelbrotFixed, pair)) {
            ++failures;
        } else {
            formula::ExpressionReferenceBuildRequest request;
            request.canonicalProgram = &pair.canonical;
            request.runtimeProgram = &pair.runtime;
            request.pixelParameter = FormulaParameter::C;
            request.center.realDecimal = "-1";
            request.center.imaginaryDecimal = "0";
            request.fixed = mandelbrotFixed;
            request.bailout = 4.0;
            request.maxIterations = 40;
            request.precision.viewBits = 1700;
            request.precision.minimumBits = 128;
            request.precision.guardBits = 64;
            request.precision.maximumBits = 4096;
            request.certificationPrecision = 1892;
            formula::ExpressionReferenceOrbitResult reference;
            bool okay = formula::buildExpressionReferenceOrbit(request, reference);
            auto contains = [&](const formula::MpfrComplex& exact, const formula::ScaledComplexShadow& primary, const formula::ScaledComplexShadow& defect, const formula::ScaledRealValue& radius) {
                formula::ScaledComplexValue midpoint;
                formula::MpfrComplex reconstructed(request.certificationPrecision);
                mpfr_t bound, difference, otherDifference;
                mpfr_inits2(request.certificationPrecision, bound, difference, otherDifference, (mpfr_ptr)0);
                bool contained = formula::makeScaledComplexValue(primary, defect, midpoint) == formula::ScaledArithmeticStatus::Success && formula::setMpfrFromScaledValue(reconstructed, midpoint) && formula::setMpfrFromScaledValue(bound, radius);
                if (contained) {
                    mpfr_sub(difference, exact.re, reconstructed.re, MPFR_RNDD);
                    mpfr_sub(otherDifference, exact.re, reconstructed.re, MPFR_RNDU);
                    mpfr_abs(difference, difference, MPFR_RNDU);
                    mpfr_abs(otherDifference, otherDifference, MPFR_RNDU);
                    contained = mpfr_cmp(difference, bound) <= 0 && mpfr_cmp(otherDifference, bound) <= 0;
                }
                if (contained) {
                    mpfr_sub(difference, exact.im, reconstructed.im, MPFR_RNDD);
                    mpfr_sub(otherDifference, exact.im, reconstructed.im, MPFR_RNDU);
                    mpfr_abs(difference, difference, MPFR_RNDU);
                    mpfr_abs(otherDifference, otherDifference, MPFR_RNDU);
                    contained = mpfr_cmp(difference, bound) <= 0 && mpfr_cmp(otherDifference, bound) <= 0;
                }
                mpfr_clears(bound, difference, otherDifference, (mpfr_ptr)0);
                return contained;
            };
            if (okay) {
                ExpressionOracleContext context(request.certificationPrecision);
                context.c.set("-1", "0");
                context.z0.set(0.0, 0.0);
                context.z.set(context.z0);
                MpfrComplex next(request.certificationPrecision);
                for (size_t sampleIndex = 0; okay && sampleIndex < reference.samples.size(); ++sampleIndex) {
                    const auto& sample = reference.samples[sampleIndex];
                    okay = contains(context.z, sample.z, sample.zDefect, sample.zError);
                    context.iteration = static_cast<int>(sampleIndex);
                    formula::ExpressionOracleTrace trace;
                    std::string error;
                    okay = okay && ExpressionOracle::evaluateTrace(pair.runtime, context, next, trace, &error) && contains(next, sample.next, sample.rootDefect, sample.nextError) && contains(next, sample.next, sample.rootDefect, sample.rootError) && trace.nodes.size() == sample.tapeCount;
                    for (size_t nodeIndex = 0; okay && nodeIndex < trace.nodes.size(); ++nodeIndex) {
                        const auto& exact = trace.nodes[nodeIndex];
                        const auto& compact = reference.tape[static_cast<size_t>(sample.tapeOffset) + nodeIndex];
                        okay = contains(exact.output, compact.output, compact.outputDefect, compact.outputError);
                        if (okay && (compact.flags & (formula::OracleTraceHasCompanion | formula::OracleTraceHasDenominator | formula::OracleTraceHasLogarithmBase))) okay = contains(exact.auxiliary, compact.auxiliary, compact.auxiliaryDefect, compact.auxiliaryError);
                    }
                    context.z.set(next);
                }
            }
            if (!okay) {
                printf("  compact certification containment failed\n");
                ++failures;
            }
        }
    }

    // Deterministic output and byte identity across thread policies.
    {
        ProgramPair pair;
        std::vector<float> single, parallel, repeated;
        if (!compilePair("z*z+c", FormulaParameter::C, mandelbrotFixed, pair)) {
            ++failures;
        } else {
            ExpressionDeepRenderRequest one = makeRequest(pair, mandelbrotFixed, FormulaParameter::C, "-2", "0", 17, 11, 850, single);
            one.threading.threads = 1;
            ExpressionDeepRenderResult oneResult;
            bool okay = formula::renderExpressionDeepFrame(one, oneResult);
            ExpressionDeepRenderRequest many = makeRequest(pair, mandelbrotFixed, FormulaParameter::C, "-2", "0", 17, 11, 850, parallel);
            many.threading.threads = 4;
            ExpressionDeepRenderResult manyResult;
            okay = okay && formula::renderExpressionDeepFrame(many, manyResult);
            ExpressionDeepRenderRequest again = makeRequest(pair, mandelbrotFixed, FormulaParameter::C, "-2", "0", 17, 11, 850, repeated);
            again.threading.threads = 4;
            ExpressionDeepRenderResult againResult;
            okay = okay && formula::renderExpressionDeepFrame(again, againResult);
            if (!okay || single != parallel || parallel != repeated || oneResult.fallbackPixelCount != manyResult.fallbackPixelCount || manyResult.fallbackPixelCount != againResult.fallbackPixelCount) {
                printf("  deterministic thread identity failed\n");
                ++failures;
            }
        }
    }

    // Moderate-scale chaotic polynomial views exercise rigorous ambiguity
    // fallback while preserving every first-escape result.
    {
        struct ViewCase {
            const char* name;
            const char* source;
            const char* centerReal;
            const char* centerImaginary;
            const char* scale;
            int iterations;
        };
        const ViewCase views[] = {{"moderate-mandelbrot", "z*z+c", "-0.743643887037151", "0.13182590420533", "64", 300}, {"moderate-cubic", "z*z*z-0.2*z+c", "-0.1", "0.65", "48", 220}};
        for (const ViewCase& view : views) {
            ProgramPair pair;
            std::vector<float> actual, expected;
            if (!compilePair(view.source, FormulaParameter::C, mandelbrotFixed, pair)) {
                ++failures;
                continue;
            }
            ExpressionDeepRenderRequest request = makeRequest(pair, mandelbrotFixed, FormulaParameter::C, view.centerReal, view.centerImaginary, 19, 13, view.iterations, actual);
            request.scale.decimal = view.scale;
            ExpressionDeepRenderResult result;
            const bool okay = formula::renderExpressionDeepFrame(request, result) && renderOracle(pair.runtime, mandelbrotFixed, FormulaParameter::C, view.centerReal, view.centerImaginary, view.scale, 19, 13, view.iterations, 4.0, result.fallbackPrecision, expected);
            if (!okay || actual != expected || result.fastPixelCount + result.fallbackPixelCount != actual.size()) {
                printf("  adversarial view failed [%s]\n", view.name);
                ++failures;
            }
        }
    }

    // Artificially enlarged intervals must only increase fallback, never
    // commit a classification that disagrees with the MPFR oracle.
    {
        ProgramPair pair;
        std::vector<float> baseline, widened, expected;
        if (!compilePair("z*z+c", FormulaParameter::C, mandelbrotFixed, pair)) {
            ++failures;
        } else {
            ExpressionDeepRenderRequest baseRequest = makeRequest(pair, mandelbrotFixed, FormulaParameter::C, "-1", "0", 21, 13, 180, baseline);
            ExpressionDeepRenderResult baseResult;
            bool okay = formula::renderExpressionDeepFrame(baseRequest, baseResult);
            ExpressionDeepRenderRequest wideRequest = makeRequest(pair, mandelbrotFixed, FormulaParameter::C, "-1", "0", 21, 13, 180, widened);
            wideRequest.verificationErrorInflationBits = 2048;
            ExpressionDeepRenderResult wideResult;
            okay = okay && formula::renderExpressionDeepFrame(wideRequest, wideResult) && renderOracle(pair.runtime, mandelbrotFixed, FormulaParameter::C, "-1", "0", "1e500", 21, 13, 180, 4.0, wideResult.fallbackPrecision, expected);
            if (!okay || baseline != expected || widened != expected || wideResult.fallbackPixelCount == 0 || wideResult.fallbackPixelCount < baseResult.fallbackPixelCount) {
                printf("  widened certification fallback failed\n");
                ++failures;
            }
        }
    }

    // Every certified jet layout must agree with Taylor-off, all-MPFR, and
    // single-thread rendering.
    {
        struct EntireFrame {
            const char* name;
            const char* source;
            FormulaParameter pixel;
            const ExpressionContext* fixed;
            const char* centerReal;
        };
        const EntireFrame frames[] = {{"sin", "sin(z)+c", FormulaParameter::C, &transcendentalFixed, "0"},
                                      {"cos", "cos(z)-1+c", FormulaParameter::C, &transcendentalFixed, "0"},
                                      {"exp", "exp(0.1*z)+c", FormulaParameter::C, &transcendentalFixed, "-1"},
                                      {"sinh", "sinh(z)+c", FormulaParameter::C, &transcendentalFixed, "0"},
                                      {"cosh", "cosh(z)-1+c", FormulaParameter::C, &transcendentalFixed, "0"},
                                      {"sin-z0", "sin(z)+c", FormulaParameter::InitialZ, &z0Fixed, "0"},
                                      {"conjugate", "conj(z)+c", FormulaParameter::C, &transcendentalFixed, "0"},
                                      {"real", "real(z)+c", FormulaParameter::C, &transcendentalFixed, "0"},
                                      {"imaginary", "imag(z)*complex(0,1)+c", FormulaParameter::C, &transcendentalFixed, "0"},
                                      {"norm", "norm(z)+c", FormulaParameter::C, &transcendentalFixed, "0"},
                                      {"make-complex", "complex(real(z),imag(c))+c", FormulaParameter::C, &transcendentalFixed, "0"},
                                      {"mixed-real-polynomial", "(z*z+c)*conj(z+c)+c", FormulaParameter::C, &transcendentalFixed, "0"},
                                      {"conjugate-z0", "conj(z)+c", FormulaParameter::InitialZ, &z0Fixed, "0"},
                                      {"divide", "z/(c+2)+c", FormulaParameter::C, &transcendentalFixed, "0"},
                                      {"reciprocal", "1/(z+2)+c", FormulaParameter::C, &transcendentalFixed, "0"},
                                      {"tan", "tan(z)+c", FormulaParameter::C, &transcendentalFixed, "0"},
                                      {"tanh", "tanh(z)+c", FormulaParameter::C, &transcendentalFixed, "0"},
                                      {"nested", "(sin(z)+c)/(cos(z)+2)", FormulaParameter::C, &transcendentalFixed, "0"},
                                      {"reciprocal-z0", "1/(z+2)+c", FormulaParameter::InitialZ, &z0Fixed, "0"},
                                      {"log", "log(z+2)+c", FormulaParameter::C, &branchFixed, "0"},
                                      {"log10", "log10(z+2)+c", FormulaParameter::C, &branchFixed, "0"},
                                      {"sqrt", "0.25*sqrt(z+2)+c", FormulaParameter::C, &branchFixed, "0"},
                                      {"power", "0.1*pow(z+2,c+0.25)+c", FormulaParameter::C, &branchFixed, "0"},
                                      {"nested-branch", "0.1*exp(log(z+2))+0.1*sin(sqrt(z+3))/(z+4)+c", FormulaParameter::C, &branchFixed, "0"}};
        for (const EntireFrame& frame : frames) {
            ProgramPair pair;
            std::vector<float> enabled, disabled, mpfr, single;
            if (!compilePair(frame.source, frame.pixel, *frame.fixed, pair)) {
                ++failures;
                continue;
            }
            ExpressionDeepRenderRequest enabledRequest = makeRequest(pair, *frame.fixed, frame.pixel, frame.centerReal, "0", 11, 7, 30, enabled);
            enabledRequest.taylor.requirePredictedBenefit = false;
            ExpressionDeepRenderResult enabledResult;
            bool okay = formula::renderExpressionDeepFrame(enabledRequest, enabledResult);

            ExpressionDeepRenderRequest disabledRequest = makeRequest(pair, *frame.fixed, frame.pixel, frame.centerReal, "0", 11, 7, 30, disabled);
            disabledRequest.taylor.enableTaylor = false;
            ExpressionDeepRenderResult disabledResult;
            okay = okay && formula::renderExpressionDeepFrame(disabledRequest, disabledResult);

            ExpressionDeepRenderRequest mpfrRequest = makeRequest(pair, *frame.fixed, frame.pixel, frame.centerReal, "0", 11, 7, 30, mpfr);
            mpfrRequest.forceMpfrFallbackForVerification = true;
            ExpressionDeepRenderResult mpfrResult;
            okay = okay && formula::renderExpressionDeepFrame(mpfrRequest, mpfrResult);

            ExpressionDeepRenderRequest singleRequest = makeRequest(pair, *frame.fixed, frame.pixel, frame.centerReal, "0", 11, 7, 30, single);
            singleRequest.threading.threads = 1;
            singleRequest.taylor.requirePredictedBenefit = false;
            ExpressionDeepRenderResult singleResult;
            okay = okay && formula::renderExpressionDeepFrame(singleRequest, singleResult);
            const uint64_t pixels = static_cast<uint64_t>(enabled.size());
            const bool meromorphic = pair.runtime.scaledResidualCapability() == formula::ExpressionScaledResidualCapability::CertifiedMeromorphicCandidate;
            const bool branch = pair.runtime.scaledResidualCapability() == formula::ExpressionScaledResidualCapability::CertifiedBranchCandidate;
            const bool real = pair.runtime.scaledResidualCapability() == formula::ExpressionScaledResidualCapability::CertifiedRealCandidate;
            if (!okay || enabled != disabled || enabled != mpfr || enabled != single || !enabledResult.taylorAccepted || enabledResult.taylorAcceptedPixelCount == 0 || enabledResult.taylorCoveredIterations != 30 || (meromorphic ? enabledResult.taylorReciprocalCount == 0 || enabledResult.taylorReciprocalOperationCount == 0 || enabledResult.taylorMinimumDenominatorClearance.isZero() : branch ? enabledResult.taylorBranchCompositionCount == 0 || enabledResult.taylorBranchCompositionOperationCount == 0 || enabledResult.taylorMinimumBranchCutClearance.isZero() || enabledResult.taylorMinimumBranchZeroClearance.isZero()
                                                                                                                                                                                                                                                                                                                                                                                              : real     ? enabledResult.taylorLayout != formula::ExpressionTaylorJetLayout::RealBivariate || enabledResult.taylorMonomialCount == 0
                                                                                                                                                                                                                                                                                                                                                                                                         : enabledResult.taylorFunctionSeriesCount == 0) ||
                disabledResult.fastPixelCount != 0 || disabledResult.fallbackPixelCount != pixels || mpfrResult.fastPixelCount != 0 || mpfrResult.fallbackPixelCount != pixels) {
                printf("  entire Taylor parity failed [%s] accepted=%d landing=%d\n", frame.name, enabledResult.taylorAccepted ? 1 : 0, enabledResult.taylorCoveredIterations);
                ++failures;
            }
        }

        ProgramPair benchmarkPair;
        if (!compilePair("sin(z)+c", FormulaParameter::C, transcendentalFixed, benchmarkPair)) {
            ++failures;
        } else {
            bool speedOkay = true;
            double jetSeconds[2]{};
            double mpfrSeconds[2]{};
            for (int repeat = 0; repeat < 2; ++repeat) {
                std::vector<float> jetOutput, mpfrOutput;
                ExpressionDeepRenderRequest jetRequest = makeRequest(benchmarkPair, transcendentalFixed, FormulaParameter::C, "0", "0", 80, 50, 120, jetOutput);
                ExpressionDeepRenderResult jetResult;
                const Clock::time_point jetStart = Clock::now();
                const bool jetOkay = formula::renderExpressionDeepFrame(jetRequest, jetResult);
                jetSeconds[repeat] = std::chrono::duration<double>(Clock::now() - jetStart).count();

                ExpressionDeepRenderRequest mpfrRequest = makeRequest(benchmarkPair, transcendentalFixed, FormulaParameter::C, "0", "0", 80, 50, 120, mpfrOutput);
                mpfrRequest.forceMpfrFallbackForVerification = true;
                ExpressionDeepRenderResult mpfrResult;
                const Clock::time_point mpfrStart = Clock::now();
                const bool mpfrOkay = formula::renderExpressionDeepFrame(mpfrRequest, mpfrResult);
                mpfrSeconds[repeat] = std::chrono::duration<double>(Clock::now() - mpfrStart).count();
                const bool repeatOkay = jetOkay && mpfrOkay && jetOutput == mpfrOutput && jetResult.taylorAccepted && jetResult.taylorAcceptedPixelCount > 0 && jetSeconds[repeat] < mpfrSeconds[repeat];
                speedOkay = speedOkay && repeatOkay;
            }
            printf("  entire Taylor e500 total jet/MPFR %.3f/%.3f and %.3f/%.3f s\n", jetSeconds[0], mpfrSeconds[0], jetSeconds[1], mpfrSeconds[1]);
            if (!speedOkay) {
                printf("  entire Taylor repeated speed gate failed\n");
                ++failures;
            }

            std::vector<float> tinyOutput;
            ExpressionDeepRenderRequest tinyRequest = makeRequest(benchmarkPair, transcendentalFixed, FormulaParameter::C, "0", "0", 3, 3, 9, tinyOutput);
            ExpressionDeepRenderResult tinyResult;
            if (!formula::renderExpressionDeepFrame(tinyRequest, tinyResult) || !tinyResult.taylorAttempted || tinyResult.taylorAccepted || tinyResult.fastPixelCount != 0 || tinyResult.fallbackPixelCount != tinyOutput.size()) {
                printf("  entire Taylor cost auto-disable failed\n");
                ++failures;
            }
        }

        ProgramPair realBenchmark;
        if (!compilePair("(z*z+c)*conj(z+c)+c", FormulaParameter::C, transcendentalFixed, realBenchmark)) {
            ++failures;
        } else {
            bool benchmarkOkay = true;
            bool accepted = false;
            bool acceptanceInitialized = false;
            double jetSeconds[2]{};
            double mpfrSeconds[2]{};
            double buildSeconds = 0.0;
            double evaluationSeconds = 0.0;
            uint64_t fallback = 0;
            size_t monomials = 0;
            uint64_t convolutionOperations = 0;
            for (int repeat = 0; repeat < 2; ++repeat) {
                std::vector<float> jetOutput, mpfrOutput;
                ExpressionDeepRenderRequest jetRequest = makeRequest(realBenchmark, transcendentalFixed, FormulaParameter::C, "0", "0", 80, 50, 120, jetOutput);
                ExpressionDeepRenderResult jetResult;
                const Clock::time_point jetStart = Clock::now();
                const bool jetOkay = formula::renderExpressionDeepFrame(jetRequest, jetResult);
                jetSeconds[repeat] = std::chrono::duration<double>(Clock::now() - jetStart).count();

                ExpressionDeepRenderRequest mpfrRequest = makeRequest(realBenchmark, transcendentalFixed, FormulaParameter::C, "0", "0", 80, 50, 120, mpfrOutput);
                mpfrRequest.forceMpfrFallbackForVerification = true;
                ExpressionDeepRenderResult mpfrResult;
                const Clock::time_point mpfrStart = Clock::now();
                const bool mpfrOkay = formula::renderExpressionDeepFrame(mpfrRequest, mpfrResult);
                mpfrSeconds[repeat] = std::chrono::duration<double>(Clock::now() - mpfrStart).count();

                if (!acceptanceInitialized) {
                    accepted = jetResult.taylorAccepted;
                    acceptanceInitialized = true;
                } else {
                    benchmarkOkay = benchmarkOkay && accepted == jetResult.taylorAccepted;
                }
                buildSeconds = jetResult.taylorBuildSeconds;
                evaluationSeconds = jetResult.taylorEvaluationSeconds;
                fallback = jetResult.fallbackPixelCount;
                monomials = jetResult.taylorMonomialCount;
                convolutionOperations = jetResult.taylorBivariateConvolutionOperationCount;
                const bool repeatOkay = jetOkay && mpfrOkay && jetOutput == mpfrOutput && (jetResult.taylorAccepted ? jetResult.taylorLayout == formula::ExpressionTaylorJetLayout::RealBivariate && jetResult.taylorAcceptedPixelCount > 0 && jetResult.taylorMonomialCount > 0 && jetResult.taylorBivariateConvolutionOperationCount > 0 && jetSeconds[repeat] <= mpfrSeconds[repeat] * 1.15 : jetResult.fastPixelCount == 0 && jetResult.fallbackPixelCount == jetOutput.size());
                benchmarkOkay = benchmarkOkay && repeatOkay;
            }
            printf("  real-bivariate Taylor e500 jet/MPFR %.3f/%.3f and %.3f/%.3f s accepted=%d build/eval=%.3f/%.3f fallback=%llu monomials=%zu convolution=%llu\n", jetSeconds[0], mpfrSeconds[0], jetSeconds[1], mpfrSeconds[1], accepted ? 1 : 0, buildSeconds, evaluationSeconds, (unsigned long long)fallback, monomials, (unsigned long long)convolutionOperations);
            if (!benchmarkOkay) {
                printf("  real-bivariate Taylor repeated speed/auto-disable gate failed\n");
                ++failures;
            }

            std::vector<float> tinyOutput;
            ExpressionDeepRenderRequest tinyRequest = makeRequest(realBenchmark, transcendentalFixed, FormulaParameter::C, "0", "0", 3, 3, 9, tinyOutput);
            ExpressionDeepRenderResult tinyResult;
            if (!formula::renderExpressionDeepFrame(tinyRequest, tinyResult) || !tinyResult.taylorAttempted || tinyResult.taylorAccepted || tinyResult.fastPixelCount != 0 || tinyResult.fallbackPixelCount != tinyOutput.size()) {
                printf("  real-bivariate Taylor cost auto-disable failed\n");
                ++failures;
            }

            std::vector<float> baseline, constrained, mpfr;
            ExpressionDeepRenderRequest baselineRequest = makeRequest(realBenchmark, transcendentalFixed, FormulaParameter::C, "0", "0", 80, 50, 300, baseline);
            baselineRequest.threading.threads = 1;
            ExpressionDeepRenderResult baselineResult;
            bool memoryOkay = formula::renderExpressionDeepFrame(baselineRequest, baselineResult) && baselineResult.taylorAccepted && baselineResult.taylorMemoryBytes > 0;
            ExpressionDeepRenderRequest mpfrRequest = makeRequest(realBenchmark, transcendentalFixed, FormulaParameter::C, "0", "0", 80, 50, 300, mpfr);
            mpfrRequest.forceMpfrFallbackForVerification = true;
            mpfrRequest.threading.threads = 1;
            ExpressionDeepRenderResult mpfrResult;
            memoryOkay = memoryOkay && formula::renderExpressionDeepFrame(mpfrRequest, mpfrResult);
            size_t constrainedLimit = 0;
            if (memoryOkay) {
                const size_t fullPeak = baselineResult.referenceBytes + baselineResult.rendererBytes;
                const size_t reduction = std::max<size_t>(1, baselineResult.taylorMemoryBytes / 2);
                constrainedLimit = fullPeak > reduction ? fullPeak - reduction : 1;
                constrainedLimit = std::max(constrainedLimit, mpfrResult.rendererBytes + 4096);
                memoryOkay = constrainedLimit < fullPeak;
            }
            ExpressionDeepRenderRequest constrainedRequest = makeRequest(realBenchmark, transcendentalFixed, FormulaParameter::C, "0", "0", 80, 50, 300, constrained);
            constrainedRequest.memory.memoryLimitBytes = constrainedLimit;
            constrainedRequest.threading.threads = 1;
            ExpressionDeepRenderResult constrainedResult;
            memoryOkay = memoryOkay && formula::renderExpressionDeepFrame(constrainedRequest, constrainedResult);
            const uint64_t pixels = static_cast<uint64_t>(constrained.size());
            if (!memoryOkay || baseline != mpfr || constrained != mpfr || constrainedResult.taylorAccepted || constrainedResult.fastPixelCount != 0 || constrainedResult.fallbackPixelCount != pixels || constrainedResult.referenceBytes != 0 || constrainedResult.rendererBytes > constrainedLimit) {
                printf("  real-bivariate Taylor memory fallback failed\n");
                ++failures;
            }
        }

        ProgramPair branchBenchmark;
        if (!compilePair("log(z+2)+c", FormulaParameter::C, branchFixed, branchBenchmark)) {
            ++failures;
        } else {
            bool speedOkay = true;
            double jetSeconds[2]{};
            double mpfrSeconds[2]{};
            uint64_t fallback = 0;
            formula::ScaledRealValue cutClearance;
            formula::ScaledRealValue zeroClearance;
            for (int repeat = 0; repeat < 2; ++repeat) {
                std::vector<float> jetOutput, mpfrOutput;
                ExpressionDeepRenderRequest jetRequest = makeRequest(branchBenchmark, branchFixed, FormulaParameter::C, "0", "0", 80, 50, 120, jetOutput);
                ExpressionDeepRenderResult jetResult;
                const Clock::time_point jetStart = Clock::now();
                speedOkay = speedOkay && formula::renderExpressionDeepFrame(jetRequest, jetResult);
                jetSeconds[repeat] = std::chrono::duration<double>(Clock::now() - jetStart).count();

                ExpressionDeepRenderRequest mpfrRequest = makeRequest(branchBenchmark, branchFixed, FormulaParameter::C, "0", "0", 80, 50, 120, mpfrOutput);
                mpfrRequest.forceMpfrFallbackForVerification = true;
                ExpressionDeepRenderResult mpfrResult;
                const Clock::time_point mpfrStart = Clock::now();
                speedOkay = speedOkay && formula::renderExpressionDeepFrame(mpfrRequest, mpfrResult);
                mpfrSeconds[repeat] = std::chrono::duration<double>(Clock::now() - mpfrStart).count();
                fallback = jetResult.fallbackPixelCount;
                cutClearance = jetResult.taylorMinimumBranchCutClearance;
                zeroClearance = jetResult.taylorMinimumBranchZeroClearance;
                speedOkay = speedOkay && jetOutput == mpfrOutput && jetResult.taylorAccepted && jetResult.taylorBranchCompositionCount > 0 && jetResult.taylorAcceptedPixelCount > 0 && jetSeconds[repeat] < mpfrSeconds[repeat];
            }
            printf("  branch Taylor e500 jet/MPFR %.3f/%.3f and %.3f/%.3f s fallback=%llu cut=(%.6g,e%lld) zero=(%.6g,e%lld)\n", jetSeconds[0], mpfrSeconds[0], jetSeconds[1], mpfrSeconds[1], (unsigned long long)fallback, cutClearance.mantissa, (long long)cutClearance.exponent, zeroClearance.mantissa, (long long)zeroClearance.exponent);
            if (!speedOkay) {
                printf("  branch Taylor repeated speed gate failed\n");
                ++failures;
            }
        }

        ProgramPair argPolarBenchmark;
        if (!compilePair("polar(real(c)+1,arg(c+2))", FormulaParameter::C, transcendentalFixed, argPolarBenchmark)) {
            ++failures;
        } else {
            bool speedOkay = true;
            bool accepted = false;
            bool acceptanceInitialized = false;
            double jetSeconds[2]{};
            double mpfrSeconds[2]{};
            uint64_t fallback = 0;
            formula::ScaledRealValue cutClearance;
            formula::ScaledRealValue zeroClearance;
            formula::ScaledRealValue radiusClearance;
            std::vector<float> acceptedOutput;
            for (int repeat = 0; repeat < 2; ++repeat) {
                std::vector<float> jetOutput, mpfrOutput;
                ExpressionDeepRenderRequest jetRequest = makeRequest(argPolarBenchmark, transcendentalFixed, FormulaParameter::C, "0", "0", 80, 50, 120, jetOutput);
                ExpressionDeepRenderResult jetResult;
                const Clock::time_point jetStart = Clock::now();
                const bool jetOkay = formula::renderExpressionDeepFrame(jetRequest, jetResult);
                jetSeconds[repeat] = std::chrono::duration<double>(Clock::now() - jetStart).count();

                ExpressionDeepRenderRequest mpfrRequest = makeRequest(argPolarBenchmark, transcendentalFixed, FormulaParameter::C, "0", "0", 80, 50, 120, mpfrOutput);
                mpfrRequest.forceMpfrFallbackForVerification = true;
                ExpressionDeepRenderResult mpfrResult;
                const Clock::time_point mpfrStart = Clock::now();
                const bool mpfrOkay = formula::renderExpressionDeepFrame(mpfrRequest, mpfrResult);
                mpfrSeconds[repeat] = std::chrono::duration<double>(Clock::now() - mpfrStart).count();

                if (!acceptanceInitialized) {
                    accepted = jetResult.taylorAccepted;
                    acceptanceInitialized = true;
                } else {
                    speedOkay = speedOkay && accepted == jetResult.taylorAccepted;
                }
                fallback = jetResult.fallbackPixelCount;
                cutClearance = jetResult.taylorMinimumBranchCutClearance;
                zeroClearance = jetResult.taylorMinimumBranchZeroClearance;
                radiusClearance = jetResult.taylorMinimumPolarRadiusClearance;
                const bool repeatOkay = jetOkay && mpfrOkay && jetOutput == mpfrOutput && jetResult.taylorAccepted && jetResult.taylorArgCompositionCount > 0 && jetResult.taylorPolarCompositionCount > 0 && jetResult.taylorAcceptedPixelCount > 0 && jetResult.fallbackPixelCount == 0 && !cutClearance.isZero() && !zeroClearance.isZero() && !radiusClearance.isZero() && jetSeconds[repeat] <= mpfrSeconds[repeat] * 1.15;
                speedOkay = speedOkay && repeatOkay;
                if (repeat == 0) acceptedOutput = jetOutput;
            }
            std::vector<float> disabledOutput;
            ExpressionDeepRenderRequest disabledRequest = makeRequest(argPolarBenchmark, transcendentalFixed, FormulaParameter::C, "0", "0", 80, 50, 120, disabledOutput);
            disabledRequest.taylor.enableTaylor = false;
            ExpressionDeepRenderResult disabledResult;
            const bool disabledOkay = formula::renderExpressionDeepFrame(disabledRequest, disabledResult);
            speedOkay = speedOkay && disabledOkay && disabledOutput == acceptedOutput && !disabledResult.taylorAttempted && disabledResult.fastPixelCount == 0 && disabledResult.fallbackPixelCount == disabledOutput.size();
            printf("  Arg/Polar Taylor e500 jet/MPFR %.3f/%.3f and %.3f/%.3f s accepted=%d fallback=%llu cut=(%.6g,e%lld) zero=(%.6g,e%lld) radius=(%.6g,e%lld)\n", jetSeconds[0], mpfrSeconds[0], jetSeconds[1], mpfrSeconds[1], accepted ? 1 : 0, (unsigned long long)fallback, cutClearance.mantissa, (long long)cutClearance.exponent, zeroClearance.mantissa, (long long)zeroClearance.exponent, radiusClearance.mantissa, (long long)radiusClearance.exponent);
            if (!speedOkay) {
                printf("  Arg/Polar Taylor repeated speed gate failed\n");
                ++failures;
            }
        }

        ProgramPair reciprocalBenchmark;
        if (!compilePair("1/(z+2)+c", FormulaParameter::C, transcendentalFixed, reciprocalBenchmark)) {
            ++failures;
        } else {
            bool speedOkay = true;
            double jetSeconds[2]{};
            double mpfrSeconds[2]{};
            uint64_t fallback = 0;
            formula::ScaledRealValue clearance;
            formula::ScaledRealValue reciprocalTail;
            for (int repeat = 0; repeat < 2; ++repeat) {
                std::vector<float> jetOutput, mpfrOutput;
                ExpressionDeepRenderRequest jetRequest = makeRequest(reciprocalBenchmark, transcendentalFixed, FormulaParameter::C, "0", "0", 80, 50, 120, jetOutput);
                ExpressionDeepRenderResult jetResult;
                const Clock::time_point jetStart = Clock::now();
                speedOkay = speedOkay && formula::renderExpressionDeepFrame(jetRequest, jetResult);
                jetSeconds[repeat] = std::chrono::duration<double>(Clock::now() - jetStart).count();

                ExpressionDeepRenderRequest mpfrRequest = makeRequest(reciprocalBenchmark, transcendentalFixed, FormulaParameter::C, "0", "0", 80, 50, 120, mpfrOutput);
                mpfrRequest.forceMpfrFallbackForVerification = true;
                ExpressionDeepRenderResult mpfrResult;
                const Clock::time_point mpfrStart = Clock::now();
                speedOkay = speedOkay && formula::renderExpressionDeepFrame(mpfrRequest, mpfrResult);
                mpfrSeconds[repeat] = std::chrono::duration<double>(Clock::now() - mpfrStart).count();
                fallback = jetResult.fallbackPixelCount;
                clearance = jetResult.taylorMinimumDenominatorClearance;
                reciprocalTail = jetResult.taylorMaximumReciprocalTail;
                speedOkay = speedOkay && jetOutput == mpfrOutput && jetResult.taylorAccepted && jetResult.taylorReciprocalCount > 0 && jetResult.taylorAcceptedPixelCount > 0 && jetSeconds[repeat] < mpfrSeconds[repeat];
            }
            printf("  meromorphic Taylor e500 jet/MPFR %.3f/%.3f and %.3f/%.3f s fallback=%llu clearance=(%.6g,e%lld) max-tail=(%.6g,e%lld)\n", jetSeconds[0], mpfrSeconds[0], jetSeconds[1], mpfrSeconds[1], (unsigned long long)fallback, clearance.mantissa, (long long)clearance.exponent, reciprocalTail.mantissa, (long long)reciprocalTail.exponent);
            if (!speedOkay) {
                printf("  meromorphic Taylor repeated speed gate failed\n");
                ++failures;
            }

            std::vector<float> baseline, constrained, mpfr;
            ExpressionDeepRenderRequest baselineRequest = makeRequest(reciprocalBenchmark, transcendentalFixed, FormulaParameter::C, "0", "0", 80, 50, 300, baseline);
            baselineRequest.threading.threads = 1;
            ExpressionDeepRenderResult baselineResult;
            bool memoryOkay = formula::renderExpressionDeepFrame(baselineRequest, baselineResult) && baselineResult.taylorAccepted && baselineResult.taylorMemoryBytes > 0;
            ExpressionDeepRenderRequest mpfrRequest = makeRequest(reciprocalBenchmark, transcendentalFixed, FormulaParameter::C, "0", "0", 80, 50, 300, mpfr);
            mpfrRequest.forceMpfrFallbackForVerification = true;
            mpfrRequest.threading.threads = 1;
            ExpressionDeepRenderResult mpfrResult;
            memoryOkay = memoryOkay && formula::renderExpressionDeepFrame(mpfrRequest, mpfrResult);
            size_t constrainedLimit = 0;
            if (memoryOkay) {
                const size_t fullPeak = baselineResult.referenceBytes + baselineResult.rendererBytes;
                const size_t reduction = std::max<size_t>(1, baselineResult.taylorMemoryBytes / 2);
                constrainedLimit = fullPeak > reduction ? fullPeak - reduction : 1;
                constrainedLimit = std::max(constrainedLimit, mpfrResult.rendererBytes + 4096);
                memoryOkay = constrainedLimit < fullPeak;
            }
            ExpressionDeepRenderRequest constrainedRequest = makeRequest(reciprocalBenchmark, transcendentalFixed, FormulaParameter::C, "0", "0", 80, 50, 300, constrained);
            constrainedRequest.memory.memoryLimitBytes = constrainedLimit;
            constrainedRequest.threading.threads = 1;
            ExpressionDeepRenderResult constrainedResult;
            memoryOkay = memoryOkay && formula::renderExpressionDeepFrame(constrainedRequest, constrainedResult);
            const uint64_t pixels = static_cast<uint64_t>(constrained.size());
            if (!memoryOkay || baseline != mpfr || constrained != mpfr || constrainedResult.taylorAccepted || constrainedResult.fastPixelCount != 0 || constrainedResult.fallbackPixelCount != pixels || constrainedResult.referenceBytes != 0 || constrainedResult.rendererBytes > constrainedLimit) {
                printf("  meromorphic Taylor memory fallback failed okay=%d status=%s error=%s limit=%zu peak=%zu Taylor=%zu attempted/accepted=%d/%d fast/fallback=%llu/%llu ref=%zu renderer=%zu parity=%d/%d\n", memoryOkay ? 1 : 0, formula::expressionDeepRenderStatusName(constrainedResult.status), constrainedResult.error.c_str(), constrainedLimit, baselineResult.referenceBytes + baselineResult.rendererBytes, baselineResult.taylorMemoryBytes, constrainedResult.taylorAttempted ? 1 : 0, constrainedResult.taylorAccepted ? 1 : 0, (unsigned long long)constrainedResult.fastPixelCount, (unsigned long long)constrainedResult.fallbackPixelCount, constrainedResult.referenceBytes, constrainedResult.rendererBytes, baseline == mpfr ? 1 : 0, constrained == mpfr ? 1 : 0);
                ++failures;
            }
        }
    }

    // The opt-in per-step series path remains informational and is never the
    // default.
    {
        ProgramPair pair;
        std::vector<float> output, expected;
        if (compilePair("sin(z)+c", FormulaParameter::C, transcendentalFixed, pair)) {
            ExpressionDeepRenderRequest request = makeRequest(pair, transcendentalFixed, FormulaParameter::C, "0", "0", 9, 7, 20, output);
            request.allowUncertifiedForBenchmark = true;
            request.taylor.enableTaylor = false;
            ExpressionDeepRenderResult result;
            if (!formula::renderExpressionDeepFrame(request, result)) {
                ++failures;
            } else {
                renderOracle(pair.runtime, transcendentalFixed, FormulaParameter::C, "0", "0", "1e500", 9, 7, 20, 4.0, result.fallbackPrecision, expected);
                size_t mismatches = 0;
                for (size_t i = 0; i < output.size(); ++i) mismatches += output[i] != expected[i];
                printf("  uncertified opt-in informational mismatches=%zu fallback=%llu\n", mismatches, (unsigned long long)result.fallbackPixelCount);
            }
        } else {
            ++failures;
        }
    }

    // A defined MPFR infinity is an escape, not an undefined domain result.
    {
        ProgramPair pair;
        std::vector<float> output;
        if (!compilePair("exp(1e100)", FormulaParameter::C, mandelbrotFixed, pair)) {
            ++failures;
        } else {
            ExpressionDeepRenderRequest request = makeRequest(pair, mandelbrotFixed, FormulaParameter::C, "0", "0", 5, 3, 2, output);
            ExpressionDeepRenderResult result;
            if (!formula::renderExpressionDeepFrame(request, result) || result.undefinedPixelCount != 0 || std::count(output.begin(), output.end(), 1.0f) != static_cast<ptrdiff_t>(output.size())) {
                printf("  MPFR infinity escape policy failed status=%s undefined=%llu ones=%td first=%.9g\n", formula::expressionDeepRenderStatusName(result.status), (unsigned long long)result.undefinedPixelCount, std::count(output.begin(), output.end(), 1.0f), output.empty() ? 0.0 : output.front());
                ++failures;
            }
        }
    }

    // Undefined orbits are not reported as successful pixels.
    {
        ProgramPair pair;
        std::vector<float> output;
        if (!compilePair("1/(z-z)", FormulaParameter::InitialZ, z0Fixed, pair)) {
            ++failures;
        } else {
            ExpressionDeepRenderRequest request = makeRequest(pair, z0Fixed, FormulaParameter::InitialZ, "1", "0", 7, 5, 8, output);
            ExpressionDeepRenderResult result;
            if (formula::renderExpressionDeepFrame(request, result) || result.status != ExpressionDeepRenderStatus::UndefinedPixel || result.undefinedPixelCount != output.size() || std::count(output.begin(), output.end(), formula::ExpressionDeepEmptyPixel) != static_cast<ptrdiff_t>(output.size())) {
                printf("  undefined pixel policy failed\n");
                ++failures;
            }
        }
    }

    // Scaled exponents extend beyond MPFR's runtime range. Near either bound,
    // the certified path must defer before a finite reference neighborhood can
    // overflow, underflow, or turn a later multiplication by zero undefined.
    {
        mpf_t centerReal, centerImaginary, scaleValue;
        mpf_init2(centerReal, 768);
        mpf_init2(centerImaginary, 768);
        mpf_init2(scaleValue, 768);
        mpf_set_ui(centerImaginary, 0);
        mpfr_t value, scale, temporary;
        mpfr_inits2(768, value, scale, temporary, (mpfr_ptr)0);

        auto runRangeCase = [&](const char* name, const char* source, ExpressionDeepFallbackReason expectedReason) {
            ProgramPair pair;
            if (!compilePair(source, FormulaParameter::C, mandelbrotFixed, pair)) {
                printf("  exponent-range compile failed [%s]\n", name);
                ++failures;
                return;
            }
            std::vector<float> actual, expected;
            ExpressionDeepRenderRequest request = makeRequest(pair, mandelbrotFixed, FormulaParameter::C, "0", "0", 3, 3, 2, actual);
            request.center.realDecimal.clear();
            request.center.imaginaryDecimal.clear();
            request.center.realMpf = centerReal;
            request.center.imaginaryMpf = centerImaginary;
            request.scale.decimal.clear();
            request.scale.mpf = scaleValue;
            ExpressionDeepRenderResult actualResult;
            const bool actualOkay = formula::renderExpressionDeepFrame(request, actualResult);

            ExpressionDeepRenderRequest oracleRequest = request;
            oracleRequest.output = nullptr;
            oracleRequest.outputCount = 0;
            if (actualResult.fallbackPrecision > request.precision.guardBits) {
                oracleRequest.precision.requestedBits = actualResult.fallbackPrecision - request.precision.guardBits;
                oracleRequest.precision.maximumBits = actualResult.fallbackPrecision;
            }
            oracleRequest.memory.fallbackGuardBits = 0;
            expected.assign(actual.size(), formula::ExpressionDeepEmptyPixel);
            oracleRequest.output = expected.data();
            oracleRequest.outputCount = expected.size();
            ExpressionDeepRenderResult expectedResult;
            const bool expectedOkay = formula::renderExpressionDeepFrame(oracleRequest, expectedResult);
            const uint64_t pixelCount = static_cast<uint64_t>(actual.size());
            if (actualResult.fallbackPrecision == 0 || actualOkay != expectedOkay || actualResult.status != expectedResult.status || actual != expected || actualResult.undefinedPixelCount != expectedResult.undefinedPixelCount || actualResult.fastPixelCount != 0 || actualResult.fallbackPixelCount != pixelCount || reasonCount(actualResult, expectedReason) != pixelCount) {
                printf("  exponent-range fallback mismatch [%s] status=%s/%s fast/fallback=%llu/%llu range=%llu error=%s\n", name, formula::expressionDeepRenderStatusName(actualResult.status), formula::expressionDeepRenderStatusName(expectedResult.status), (unsigned long long)actualResult.fastPixelCount, (unsigned long long)actualResult.fallbackPixelCount, (unsigned long long)reasonCount(actualResult, ExpressionDeepFallbackReason::ExponentRange), actualResult.error.c_str());
                ++failures;
            }
        };

        const mpfr_exp_t emax = mpfr_get_emax();
        const mpfr_exp_t emin = mpfr_get_emin();

        mpfr_set_ui_2exp(value, 1, emax - 2, MPFR_RNDN);
        mpfr_ui_div(scale, 1, value, MPFR_RNDN);
        mpfr_get_f(centerReal, value, MPFR_RNDN);
        mpfr_get_f(scaleValue, scale, MPFR_RNDN);
        runRangeCase("two-c-times-zero", "(2*c)*0", ExpressionDeepFallbackReason::ExponentRange);
        runRangeCase("huge-divide-unit", "c/c", ExpressionDeepFallbackReason::CertificationFailure);
        runRangeCase("huge-reciprocal", "1/c", ExpressionDeepFallbackReason::CertificationFailure);

        const mpfr_exp_t hugePower = (emax - 2) / 2;
        mpfr_set_ui_2exp(value, 1, hugePower, MPFR_RNDN);
        mpfr_ui_div(scale, 1, value, MPFR_RNDN);
        mpfr_get_f(centerReal, value, MPFR_RNDN);
        mpfr_get_f(scaleValue, scale, MPFR_RNDN);
        runRangeCase("huge-multiply", "c*c", ExpressionDeepFallbackReason::ExponentRange);

        const mpfr_exp_t tinyPower = emin / 2;
        mpfr_set_ui_2exp(value, 1, tinyPower, MPFR_RNDN);
        mpfr_get_f(centerReal, value, MPFR_RNDN);
        mpf_set_ui(scaleValue, 1);
        runRangeCase("tiny-underflow", "(c*c)*0.125", ExpressionDeepFallbackReason::ExponentRange);
        runRangeCase("tiny-divide-unit", "c/c", ExpressionDeepFallbackReason::CertificationFailure);
        runRangeCase("tiny-reciprocal", "1/c", ExpressionDeepFallbackReason::CertificationFailure);

        mpfr_const_log2(value, MPFR_RNDN);
        mpfr_mul_si(temporary, value, static_cast<long>(emax), MPFR_RNDN);
        mpfr_sub_ui(temporary, temporary, 1, MPFR_RNDN);
        mpfr_get_f(centerReal, temporary, MPFR_RNDN);
        mpf_set_ui(scaleValue, 1);
        runRangeCase("exp-overflow-domain", "exp(c)*0", ExpressionDeepFallbackReason::CertificationFailure);

        mpfr_clears(value, scale, temporary, (mpfr_ptr)0);
        mpf_clear(centerReal);
        mpf_clear(centerImaginary);
        mpf_clear(scaleValue);
    }

    // Cancellation leaves all unfinished pixels at EMPTY in both phases.
    auto cancellationCase = [&](const char* source, bool fallbackPhase) {
        const bool branchFast = std::string(source) == "log(z+2)+c";
        const bool realFast = std::string(source) == "conj(z)+c";
        ProgramPair pair;
        const FormulaParameter pixel = FormulaParameter::C;
        if (!compilePair(source, pixel, mandelbrotFixed, pair)) {
            ++failures;
            return;
        }
        std::vector<float> output;
        ExpressionDeepRenderRequest request = makeRequest(pair, mandelbrotFixed, pixel, fallbackPhase ? "0" : (branchFast || realFast) ? "0"
                                                                                                                                       : "-2",
                                                          "0", fallbackPhase ? 96 : (branchFast || realFast) ? 400
                                                                                                             : 48,
                                                          fallbackPhase ? 64 : (branchFast || realFast) ? 240
                                                                                                        : 28,
                                                          fallbackPhase ? 30 : 850, output);
        std::atomic<int> phase{static_cast<int>(ExpressionDeepRenderPhase::Reference)};
        std::atomic<int> polls{0};
        request.progress = [&](ExpressionDeepRenderPhase value, uint64_t, uint64_t) { phase.store(static_cast<int>(value), std::memory_order_release); };
        request.shouldCancel = [&] {
            const ExpressionDeepRenderPhase current = static_cast<ExpressionDeepRenderPhase>(phase.load(std::memory_order_acquire));
            const bool target = fallbackPhase ? current == ExpressionDeepRenderPhase::Fallback : current == ExpressionDeepRenderPhase::Fast;
            return target && polls.fetch_add(1, std::memory_order_relaxed) > (fallbackPhase ? 20 : (branchFast || realFast) ? 0
                                                                                                                            : 300);
        };
        ExpressionDeepRenderResult result;
        const bool okay = formula::renderExpressionDeepFrame(request, result);
        const size_t empty = static_cast<size_t>(std::count(output.begin(), output.end(), formula::ExpressionDeepEmptyPixel));
        if (okay || !result.cancelled || result.status != ExpressionDeepRenderStatus::Cancelled || empty == 0) {
            printf("  %s cancellation failed [%s] empty=%zu\n", fallbackPhase ? "fallback" : "fast", source, empty);
            ++failures;
        }
    };
    cancellationCase("z*z+c", false);
    cancellationCase("sin(z)+c", false);
    cancellationCase("log(z+2)+c", false);
    cancellationCase("conj(z)+c", false);
    cancellationCase("0.001/(c-c)", true);

    // Worker and per-iteration allocation failures must not let any thread
    // bypass an OpenMP worksharing barrier.
    auto workerFaultCase = [&](const char* source, ExpressionDeepVerificationFault fault, const char* name) {
        ProgramPair pair;
        if (!compilePair(source, FormulaParameter::C, mandelbrotFixed, pair)) {
            ++failures;
            return;
        }
        std::vector<float> output;
        const bool largeFastFrame = std::string(source) == "1/(z+2)+c" || std::string(source) == "log(z+2)+c" || std::string(source) == "norm(z)+c";
        ExpressionDeepRenderRequest request = makeRequest(pair, mandelbrotFixed, FormulaParameter::C, largeFastFrame ? "0" : "-1", "0", largeFastFrame ? 80 : 16, largeFastFrame ? 50 : 10, largeFastFrame ? 120 : 40, output);
        request.threading.threads = 4;
        request.verificationFault = fault;
        ExpressionDeepRenderResult result;
        const Clock::time_point start = Clock::now();
        const bool okay = formula::renderExpressionDeepFrame(request, result);
        const double elapsed = std::chrono::duration<double>(Clock::now() - start).count();
        if (okay || result.status != ExpressionDeepRenderStatus::ResourceLimit || elapsed > 10.0) {
            printf("  %s fault handling failed status=%s elapsed=%.3f\n", name, formula::expressionDeepRenderStatusName(result.status), elapsed);
            ++failures;
        }
    };
    workerFaultCase("z*z+c", ExpressionDeepVerificationFault::FastWorkerAllocation, "fast worker allocation");
    workerFaultCase("z*z+c", ExpressionDeepVerificationFault::FastIterationAllocation, "fast iteration allocation");
    workerFaultCase("sin(z)+c", ExpressionDeepVerificationFault::FastWorkerAllocation, "entire fast worker allocation");
    workerFaultCase("sin(z)+c", ExpressionDeepVerificationFault::FastIterationAllocation, "entire fast iteration allocation");
    workerFaultCase("1/(z+2)+c", ExpressionDeepVerificationFault::FastWorkerAllocation, "meromorphic fast worker allocation");
    workerFaultCase("1/(z+2)+c", ExpressionDeepVerificationFault::FastIterationAllocation, "meromorphic fast iteration allocation");
    workerFaultCase("log(z+2)+c", ExpressionDeepVerificationFault::FastWorkerAllocation, "branch fast worker allocation");
    workerFaultCase("log(z+2)+c", ExpressionDeepVerificationFault::FastIterationAllocation, "branch fast iteration allocation");
    workerFaultCase("norm(z)+c", ExpressionDeepVerificationFault::FastWorkerAllocation, "real-bivariate fast worker allocation");
    workerFaultCase("norm(z)+c", ExpressionDeepVerificationFault::FastIterationAllocation, "real-bivariate fast iteration allocation");
    workerFaultCase("0.001/(c-c)", ExpressionDeepVerificationFault::FallbackWorkerAllocation, "fallback worker allocation");
    workerFaultCase("0.001/(c-c)", ExpressionDeepVerificationFault::FallbackIterationAllocation, "fallback iteration allocation");

    // Request, resource, and semantic validation.
    {
        ProgramPair pair;
        std::vector<float> output;
        if (!compilePair("z*z+c", FormulaParameter::C, mandelbrotFixed, pair)) {
            ++failures;
        } else {
            ExpressionDeepRenderRequest request = makeRequest(pair, mandelbrotFixed, FormulaParameter::C, "-1", "0", 7, 5, 20, output);
            ExpressionDeepRenderResult result;
            request.outputCount = output.size() - 1;
            if (formula::renderExpressionDeepFrame(request, result) || result.status != ExpressionDeepRenderStatus::InvalidRequest) ++failures;
            request.outputCount = output.size();
            request.scale.decimal = "-1";
            if (formula::renderExpressionDeepFrame(request, result) || result.status != ExpressionDeepRenderStatus::InvalidRequest) ++failures;
            request.scale.decimal = "1e500";
            request.memory.memoryLimitBytes = 1;
            if (formula::renderExpressionDeepFrame(request, result) || result.status != ExpressionDeepRenderStatus::ResourceLimit) ++failures;
            request.memory.memoryLimitBytes = size_t{1} << 30;
            request.precision.maximumBits = 128;
            if (formula::renderExpressionDeepFrame(request, result) || result.status != ExpressionDeepRenderStatus::PrecisionOutOfRange) ++failures;

            ProgramPair mismatch;
            ExpressionError error;
            mismatch.canonical.compile("z-c", &error);
            mismatch.canonical.specialize(mandelbrotFixed, FormulaParameter::C, mismatch.runtime, &error);
            request.precision.maximumBits = 4096;
            request.canonicalProgram = &pair.canonical;
            request.runtimeProgram = &mismatch.runtime;
            if (formula::renderExpressionDeepFrame(request, result) || result.status != ExpressionDeepRenderStatus::ProgramMismatch) ++failures;

            const char* allFallbackMismatchSources[] = {"sin(z)+c", "log(z+2)+c", "z/(c+2)+c"};
            for (const char* source : allFallbackMismatchSources) {
                ProgramPair allFallback;
                std::vector<float> mismatchOutput;
                if (!compilePair(source, FormulaParameter::C, mandelbrotFixed, allFallback)) {
                    ++failures;
                    continue;
                }
                ExpressionDeepRenderRequest mismatchRequest = makeRequest(allFallback, mandelbrotFixed, FormulaParameter::C, "0", "0", 7, 5, 8, mismatchOutput);
                mismatchRequest.runtimeProgram = &mismatch.runtime;
                ExpressionDeepRenderResult mismatchResult;
                if (formula::renderExpressionDeepFrame(mismatchRequest, mismatchResult) || mismatchResult.status != ExpressionDeepRenderStatus::ProgramMismatch || mismatchResult.fastPixelCount != 0 || mismatchResult.fallbackPixelCount != 0) {
                    printf("  all-fallback mismatch was accepted [%s]\n", source);
                    ++failures;
                }
            }

            request.canonicalProgram = &pair.canonical;
            request.runtimeProgram = &pair.runtime;
            request.precision.maximumBits = 4096;
            request.memory.fallbackGuardBits = 4096;
            request.scale.decimal = "1e500";
            if (!formula::renderExpressionDeepFrame(request, result) || result.fallbackPixelCount != 0) ++failures;
        }
    }

    // A production-sized-enough timing sample: reference is reported
    // separately, and the parallel scaled pass must beat all-pixel MPFR.
    double benchmarkReference = 0.0;
    double benchmarkScaled = 0.0;
    double benchmarkMpfr = 0.0;
    double benchmarkSpeedup = 0.0;
    double benchmarkFastPixelsPerSecond = 0.0;
    double benchmarkFallbackRate = 0.0;
    size_t benchmarkMemory = 0;
    uint64_t benchmarkFallback = 0;
    {
        ProgramPair pair;
        if (!compilePair("z*z+c", FormulaParameter::C, mandelbrotFixed, pair)) {
            ++failures;
        } else {
            constexpr int BW = 64;
            constexpr int BH = 40;
            constexpr int BI = 350;
            std::vector<float> first, second;
            ExpressionDeepRenderRequest firstRequest = makeRequest(pair, mandelbrotFixed, FormulaParameter::C, "-1", "0", BW, BH, BI, first);
            firstRequest.threading.tileWidth = 16;
            firstRequest.threading.tileHeight = 8;
            ExpressionDeepRenderResult firstResult;
            bool okay = formula::renderExpressionDeepFrame(firstRequest, firstResult);
            ExpressionDeepRenderRequest secondRequest = makeRequest(pair, mandelbrotFixed, FormulaParameter::C, "-1", "0", BW, BH, BI, second);
            secondRequest.threading.tileWidth = 16;
            secondRequest.threading.tileHeight = 8;
            ExpressionDeepRenderResult secondResult;
            okay = okay && formula::renderExpressionDeepFrame(secondRequest, secondResult);

            std::vector<float> oracleA, oracleB;
            const Clock::time_point oracleStartA = Clock::now();
            okay = okay && renderOracle(pair.runtime, mandelbrotFixed, FormulaParameter::C, "-1", "0", "1e500", BW, BH, BI, 4.0, firstResult.fallbackPrecision, oracleA);
            const double oracleSecondsA = std::chrono::duration<double>(Clock::now() - oracleStartA).count();
            const Clock::time_point oracleStartB = Clock::now();
            okay = okay && renderOracle(pair.runtime, mandelbrotFixed, FormulaParameter::C, "-1", "0", "1e500", BW, BH, BI, 4.0, firstResult.fallbackPrecision, oracleB);
            const double oracleSecondsB = std::chrono::duration<double>(Clock::now() - oracleStartB).count();
            benchmarkReference = std::max(firstResult.referenceSeconds, secondResult.referenceSeconds);
            benchmarkScaled = std::max(firstResult.fastSeconds, secondResult.fastSeconds);
            benchmarkMpfr = std::min(oracleSecondsA, oracleSecondsB);
            benchmarkSpeedup = benchmarkScaled > 0.0 ? benchmarkMpfr / benchmarkScaled : 0.0;
            benchmarkFastPixelsPerSecond = benchmarkScaled > 0.0 ? static_cast<double>(BW * BH) / benchmarkScaled : 0.0;
            benchmarkMemory = firstResult.referenceBytes + firstResult.rendererBytes;
            benchmarkFallback = firstResult.fallbackPixelCount;
            benchmarkFallbackRate = static_cast<double>(benchmarkFallback) / static_cast<double>(BW * BH);
            if (!okay || first != second || first != oracleA || oracleA != oracleB || firstResult.fallbackPixelCount != 0 || secondResult.fallbackPixelCount != 0 || !(benchmarkSpeedup > 1.0)) {
                printf("  deep renderer benchmark gate failed speedup=%.2fx fallback=%llu\n", benchmarkSpeedup, (unsigned long long)benchmarkFallback);
                ++failures;
            }
        }
    }

    if (!(sawInitialEscape && sawInterior && sawExterior)) {
        printf("  initial/interior/exterior coverage failed\n");
        ++failures;
    }

    printf("=== expression deep frame renderer e500\n");
    printf("  exact frames=%d certified entire/meromorphic/branch/real/piecewise jets plus MPFR fallback\n", exactFrames);
    printf("  64x40 reference/scaled/all-MPFR %.3f/%.3f/%.3f s speedup %.2fx\n", benchmarkReference, benchmarkScaled, benchmarkMpfr, benchmarkSpeedup);
    printf("  certified fast rate=%.0f px/s fallback=%llu (%.3f%%) memory=%zu bytes\n", benchmarkFastPixelsPerSecond, (unsigned long long)benchmarkFallback, benchmarkFallbackRate * 100.0, benchmarkMemory);
    printf("  undefined pixels remain EMPTY\n");
    printf("  => %s\n\n", failures == 0 ? "PASS" : "CHECK (deep renderer failure)");
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
    direct.z = {0.7, 0.2};
    direct.c = {1.2, 0.3};
    direct.z0 = {0.9, -0.4};
    direct.parameters[0] = {0.6, 0.1};
    direct.parameters[1] = {0.8, 0.25};
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

    const char* expressions[] = {"z*z+c+p0*n/10", "sin(z)+cos(c)+sinh(p0)-tanh(z0)", "exp(z/3)+log(c)+log10(p0)+sqrt(z0)", "pow(z,c)+conj(p0)+complex(real(z),imag(c))+polar(abs(p1),arg(p1))", "abs(z)+norm(c)+arg(p0)"};

    int failures = 0;
    double maxDoubleError = 0.0;
    double maxPrecisionDelta = 0.0;
    mpfr_t delta;
    mpfr_init2(delta, 512);
    for (const char* source : expressions) {
        ExpressionProgram program;
        ExpressionError compileError;
        if (!program.compile(source, &compileError)) {
            printf("  oracle compile failed @ %zu: %s [%s]\n", compileError.position, compileError.message.c_str(), source);
            ++failures;
            continue;
        }
        MpfrComplex result256(256), result512(512);
        std::string error256, error512;
        bool ok256 = ExpressionOracle::evaluate(program, hp256, result256, &error256);
        bool ok512 = ExpressionOracle::evaluate(program, hp512, result512, &error512);
        if (!ok256 || !ok512) {
            printf("  oracle domain failure: %s / %s [%s]\n", error256.c_str(), error512.c_str(), source);
            ++failures;
            continue;
        }

        mpfr_sub(delta, result512.re, result256.re, MPFR_RNDN);
        mpfr_abs(delta, delta, MPFR_RNDN);
        double reDelta = mpfr_get_d(delta, MPFR_RNDU);
        mpfr_sub(delta, result512.im, result256.im, MPFR_RNDN);
        mpfr_abs(delta, delta, MPFR_RNDN);
        double imDelta = mpfr_get_d(delta, MPFR_RNDU);
        maxPrecisionDelta = std::max({maxPrecisionDelta, reDelta, imDelta});
        if (reDelta > 1e-60 || imDelta > 1e-60) ++failures;

        Complex expected = result512.toDouble();
        Complex actual = program.evaluate(direct);
        double relative = std::abs(actual - expected) / std::max(1.0, std::abs(expected));
        maxDoubleError = std::max(maxDoubleError, relative);
        if (!(relative <= 2e-13)) ++failures;
    }
    struct DomainCase {
        const char* source;
    };
    const DomainCase invalid[] = {{"log(0)"}, {"1/0"}, {"polar(-1,0)"}};
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
        if (ExpressionOracle::evaluate(program, hp256, result, &error) || !mpfr_nan_p(result.re) || !mpfr_nan_p(result.im) || error.empty()) ++domainAccepted;
    }
    failures += domainAccepted;

    {
        MpfrComplex componentSquares(256);
        mpfr_sqr(componentSquares.re, hp256.z.re, MPFR_RNDN);
        mpfr_sqr(componentSquares.im, hp256.z.im, MPFR_RNDN);
        ExpressionProgram reusable;
        ExpressionError reusableError;
        MpfrComplex result(256);
        if (!reusable.compile("sqr(complex(abs(real(z)),abs(imag(z))))", &reusableError) || !ExpressionOracle::evaluateOrbitStep(reusable, hp256, result, &componentSquares) || !reusable.compile("real(z)", &reusableError) || !ExpressionOracle::evaluateOrbitStep(reusable, hp256, result, &componentSquares) || mpfr_cmp(result.re, hp256.z.re) != 0 || !mpfr_zero_p(result.im) || !reusable.compile("z", &reusableError) || !ExpressionOracle::evaluateOrbitStep(reusable, hp256, result, &componentSquares) || mpfr_cmp(result.re, hp256.z.re) != 0 || mpfr_cmp(result.im, hp256.z.im) != 0) {
            printf("  orbit workspace revision/projection regression failed\n");
            ++failures;
        }
        alignas(ExpressionProgram) unsigned char storage[sizeof(ExpressionProgram)];
        ExpressionProgram* first = ::new (storage) ExpressionProgram();
        const bool firstOkay = first->compile("sqr(complex(abs(real(z)),abs(imag(z))))", &reusableError) && ExpressionOracle::evaluateOrbitStep(*first, hp256, result, &componentSquares);
        first->~ExpressionProgram();
        ExpressionProgram* second = ::new (storage) ExpressionProgram();
        const bool secondOkay = second->compile("real(z)", &reusableError) && ExpressionOracle::evaluateOrbitStep(*second, hp256, result, &componentSquares) && mpfr_cmp(result.re, hp256.z.re) == 0 && mpfr_zero_p(result.im);
        second->~ExpressionProgram();
        if (!firstOkay || !secondOkay) {
            printf("  orbit workspace lifetime identity regression failed\n");
            ++failures;
        }
        ExpressionProgram moveSource;
        const bool moveSetup = moveSource.compile("z", &reusableError);
        ExpressionProgram moveTarget(std::move(moveSource));
        if (!moveSetup || moveSource.valid() || !moveTarget.valid() || !ExpressionOracle::evaluateOrbitStep(moveTarget, hp256, result, &componentSquares) || mpfr_cmp(result.re, hp256.z.re) != 0 || mpfr_cmp(result.im, hp256.z.im) != 0) {
            printf("  moved program state regression failed\n");
            ++failures;
        }
        ExpressionOracle::releaseThreadWorkspace();
    }

    ExpressionProgram zeroPower;
    ExpressionError zeroPowerError;
    if (!zeroPower.compile("0^2 + 0^0", &zeroPowerError)) {
        ++failures;
    } else {
        MpfrComplex result(256);
        std::string error;
        if (!ExpressionOracle::evaluate(zeroPower, hp256, result, &error) || std::abs(result.toDouble() - Complex{1.0, 0.0}) > 1e-15) ++failures;
    }

    auto evaluateSpecial = [&](const char* source, ExpressionOracleContext& context, Complex expected, double tolerance) {
        ExpressionProgram program;
        ExpressionError compileError;
        if (!program.compile(source, &compileError)) {
            ++failures;
            return;
        }
        MpfrComplex result(context.z.precision());
        std::string error;
        if (!ExpressionOracle::evaluate(program, context, result, &error)) {
            printf("  special domain failure: %s [%s]\n", error.c_str(), source);
            ++failures;
            return;
        }
        Complex actual = result.toDouble();
        if (std::abs(actual - expected) > tolerance * std::max(1.0, std::abs(expected))) {
            printf("  special mismatch [%s]: actual=(%.17g,%.17g) expected=(%.17g,%.17g)\n", source, actual.real(), actual.imag(), expected.real(), expected.imag());
            ++failures;
        }
    };

    ExpressionOracleContext sqrtContext(256);
    sqrtContext.z.set("-1", "1e-40");
    evaluateSpecial("sqrt(z)", sqrtContext, {5e-41, 1.0}, 1e-14);

    ExpressionOracleContext tinyDivision(256);
    mpfr_set_ui_2exp(tinyDivision.z.re, 1, -600000000, MPFR_RNDN);
    mpfr_set_zero(tinyDivision.z.im, 0);
    mpfr_set_ui_2exp(tinyDivision.c.re, 1, -600000000, MPFR_RNDN);
    mpfr_set_zero(tinyDivision.c.im, 0);
    evaluateSpecial("z/c", tinyDivision, {1.0, 0.0}, 1e-15);

    ExpressionOracleContext hugeDivision(256);
    mpfr_set_ui(hugeDivision.z.re, 1, MPFR_RNDN);
    mpfr_set_exp(hugeDivision.z.re, mpfr_get_emax() - 1);
    mpfr_set(hugeDivision.z.im, hugeDivision.z.re, MPFR_RNDN);
    hugeDivision.c.set(hugeDivision.z);
    evaluateSpecial("z/c", hugeDivision, {1.0, 0.0}, 1e-15);

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
            if (!mpfr_equal_p(actual.re, delta) || !mpfr_zero_p(actual.im)) ++failures;
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
            if (!ExpressionOracle::evaluate(quotientProgram, mixedQuotient, actual, &error) || !mpfr_number_p(actual.re) || mpfr_zero_p(actual.re) || !mpfr_zero_p(actual.im)) ++failures;
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
                mpfr_init2(difference, 256);
                mpfr_init2(ulp, 256);
                mpfr_sub(difference, actual.re, delta, MPFR_RNDN);
                mpfr_abs(difference, difference, MPFR_RNDN);
                mpfr_set_ui(ulp, 1, MPFR_RNDN);
                mpfr_exp_t ulpExponent = mpfr_get_exp(actual.re) - (mpfr_exp_t)mpfr_get_prec(actual.re) + 1;
                mpfr_set_exp(ulp, ulpExponent);
                if (mpfr_cmp(difference, ulp) > 0 || !mpfr_zero_p(actual.im)) {
                    printf("  huge sqrt exceeds 1 ulp; diff-exp=%ld expected-exp=%ld\n", mpfr_zero_p(difference) ? 0L : (long)mpfr_get_exp(difference), (long)mpfr_get_exp(delta));
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
                    if (!ExpressionOracle::evaluate(hugeSqrtProgram, skewSqrt, actual, &error) || mpfr_zero_p(actual.im) || !mpfr_number_p(actual.re) || !mpfr_number_p(actual.im)) ++failures;
                }
                mpfr_clear(difference);
                mpfr_clear(ulp);
            }
        }
    }

    ExpressionOracleContext largeTan(256);
    largeTan.z.set(1.0, 1e9);
    evaluateSpecial("tan(z)", largeTan, {0.0, 1.0}, 1e-14);
    largeTan.z.set(1e9, 1.0);
    evaluateSpecial("tanh(z)", largeTan, {1.0, 0.0}, 1e-14);
    ExpressionOracleContext tinyTan(2048);
    mpfr_set_zero(tinyTan.z.re, 0);
    mpfr_set_ui_2exp(tinyTan.z.im, 1, -1000, MPFR_RNDN);
    ExpressionProgram tinyIdentityProgram, tinyTanProgram, tinyTanhProgram;
    ExpressionError tinyError;
    if (!tinyIdentityProgram.compile("z", &tinyError) || !tinyTanProgram.compile("tan(z)", &tinyError) || !tinyTanhProgram.compile("tanh(z)", &tinyError)) {
        ++failures;
    } else {
        MpfrComplex identityResult(2048);
        std::string identityError;
        if (!ExpressionOracle::evaluate(tinyIdentityProgram, tinyTan, identityResult, &identityError)) ++failures;
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
    auto compareOrbitGrid = [&](const char* source, Complex parameter0, int width, int height, int maxIterations) {
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
                Complex pixel{-2.0 + 4.0 * x / (width - 1), -1.3 + 2.6 * y / (height - 1)};
                ExpressionContext dc;
                dc.z = dc.z0 = {};
                dc.c = pixel;
                dc.parameters[0] = parameter0;
                bool doubleEscaped = false;
                int doubleIteration = maxIterations;
                for (int n = 0; n < maxIterations; ++n) {
                    dc.iteration = n;
                    dc.z = program.evaluate(dc);
                    if (!std::isfinite(dc.z.real()) || !std::isfinite(dc.z.imag()) || std::hypot(dc.z.real(), dc.z.imag()) > 4.0) {
                        doubleEscaped = true;
                        doubleIteration = n + 1;
                        break;
                    }
                }

                ExpressionOracleContext hc(256);
                hc.z.set(0.0);
                hc.z0.set(0.0);
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
                if (doubleEscaped && oracleEscaped) stats.maxIterationDiff = std::max(stats.maxIterationDiff, std::abs(doubleIteration - oracleIteration));
            }
        }
        mpfr_clear(magnitude);
        return stats;
    };

    OrbitStats quadraticOrbit = compareOrbitGrid("z*z+c", {}, 30, 20, 300);
    OrbitStats parameterOrbit = compareOrbitGrid("z*z+c+p0*z*z*z", {0.03, -0.02}, 24, 16, 200);
    failures += quadraticOrbit.classMismatch + quadraticOrbit.oracleDomainErrors;
    failures += parameterOrbit.classMismatch + parameterOrbit.oracleDomainErrors;
    if (quadraticOrbit.maxIterationDiff > 0 || parameterOrbit.maxIterationDiff > 1) ++failures;

    printf("=== MPFR expression oracle\n");
    printf("  expressions=%zu max double relative error=%.6g\n", sizeof(expressions) / sizeof(expressions[0]), maxDoubleError);
    printf("  256-vs-512 max delta=%.6g domain accepted=%d\n", maxPrecisionDelta, domainAccepted);
    printf("  orbit quadratic mismatch=%d iter-diff=%d; parameter mismatch=%d iter-diff=%d\n", quadraticOrbit.classMismatch, quadraticOrbit.maxIterationDiff, parameterOrbit.classMismatch, parameterOrbit.maxIterationDiff);
    printf("  => %s\n\n", failures == 0 ? "PASS" : "CHECK (MPFR oracle failure)");
    return failures == 0 ? 0 : 1;
}

static int runGenericFormulaProfile() {
    struct ProfileCase {
        const char* name;
        const char* source;
        formula::Complex fixedZ0;
        formula::Complex p0;
        double bailout;
        int mxit;
        formula::Complex p1{};
        bool orbitStage = false;
        bool expectOrbitInvariant = false;
    };
    const ProfileCase cases[] = {{"arithmetic", "z*z+c+p0*z", {}, {0.1, -0.05}, 4.0, 1200}, {"iteration", "z*z+c+0.0001*n", {}, {}, 4.0, 800}, {"components", "sqr(complex(real(z),imag(z)))+c+p0", {}, {0.02, -0.01}, 4.0, 800}, {"sine", "sin(z)+c+complex(1e-100,-1e-101)", {}, {}, 8.0, 400}, {"nested-sine", "sin(sin(z))+c+complex(1e-100,-1e-101)", {}, {}, 100.0, 400}, {"nested-cosine", "cos(cos(z))+c+complex(1e-100,-1e-101)", {}, {}, 100.0, 400}, {"nested-tangent", "tan(tan(z))+c+complex(1e-100,-1e-101)", {}, {}, 100.0, 400}, {"branch-power", "exp(p0*log(z))+c", {0.3, 0.2}, {2.5, 0.0}, 8.0, 300}, {"invariant-fn", "z*z+c+sin(p0)+exp(p1)", {}, {0.15, -0.2}, 8.0, 800, {-0.35, 0.1}}, {"invariant-components", "z*z+c+complex(real(p0),imag(p1))*z", {}, {0.15, -0.2}, 4.0, 1000, {-0.35, 0.1}}, {"orbit-sin-c", "0.5*z+0.1*sin(c)", {}, {}, 4.0, 800, {}, true, true}, {"orbit-sin-c-cse", "0.5*z+0.1*(sin(c)+sin(c))", {}, {}, 4.0, 800, {}, true, true}, {"orbit-control-z-n", "0.5*z+0.1*sin(z)+0.000001*n+c", {}, {}, 4.0, 500, {}, true, false}};
    constexpr int W = 322, H = 216, SAMPLES = 7, PAIR_REPEATS = 5;
    int failures = 0;
    int exactHybridCases = 0;
    int hybridNoRegressionCases = 0;
    int refillNoRegressionCases = 0;
    int unaffectedOrbitExact = 0;
    int unaffectedOrbitQuality = 0;
    int unaffectedOrbitNoRegression = 0;
    int orbitInvariantExact = 0;
    int orbitInvariantImproved = 0;
    int orbitInvariantNoRegression = 0;
    int orbitControlExact = 0;
    int orbitControlNoRegression = 0;
    printf("=== generic formula full-frame profile (%dx%d)\n", W, H);
    const char* vectorTranscendentalBenchmark = getenv("MANDEL_FORMULA_BENCH_VECTOR_TRANSCENDENTALS");
    const char* vectorTranscendentalSetting = vectorTranscendentalBenchmark ? vectorTranscendentalBenchmark : "1";
    const bool vectorTranscendentalsEnabled = std::atoi(vectorTranscendentalSetting) != 0;

    for (const ProfileCase& test : cases) {
        formula::ExpressionProgram program;
        formula::ExpressionError error;
        if (!program.compile(test.source, &error)) {
            ++failures;
            continue;
        }
        formula::ExpressionContext fixed;
        fixed.z0 = test.fixedZ0;
        fixed.parameters[0] = test.p0;
        fixed.parameters[1] = test.p1;
        formula::ExpressionProgram specializedProgram;
        if (!program.specialize(fixed, FormulaParameter::C, specializedProgram, &error)) {
            ++failures;
            continue;
        }
        formula::ExpressionJit4 jit;
        bool jitAvailable = jit.compile(program);
        formula::ExpressionJit4 specializedJit;
        bool specializedJitAvailable = specializedProgram.fastPath() == formula::ExpressionProgram::FastPath::None && specializedJit.compile(specializedProgram);
        formula::ExpressionOrbitPlan orbitPlan;
        if (!orbitPlan.build(specializedProgram, &error)) {
            ++failures;
            continue;
        }
        formula::ExpressionJit4 orbitJit;
        bool orbitJitAvailable = orbitPlan.profitable() && orbitJit.compile(orbitPlan);

        auto render = [&](const formula::ExpressionProgram& activeProgram, bool vector, const formula::ExpressionJit4* activeJit, std::vector<float>& output, double& elapsed, bool refill = true, const formula::ExpressionOrbitPlan* activePlan = nullptr) {
            output.assign((size_t)W * H, EMPTYPIXEL);
            Mandel renderer(W, H, test.mxit, 1, output.data());
            mpf_t centerRe, centerIm, scale;
            mpf_init_set_ui(centerRe, 0);
            mpf_init_set_ui(centerIm, 0);
            mpf_init_set_ui(scale, 1);
            _putenv_s("MANDEL_EXPR_VECTOR", vector ? "1" : "0");
            _putenv_s("MANDEL_EXPR_HYBRID_REFILL", refill ? "1" : "0");
            _putenv_s("MANDEL_EXPR_VECTOR_TRANSCENDENTALS", vectorTranscendentalSetting);
            auto begin = Clock::now();
            bool okay = renderer.ComputeExpression(centerRe, centerIm, scale, activeProgram, fixed, FormulaParameter::C, test.mxit, test.bailout, formula::ExpressionColoring::Raw, activeJit, activePlan);
            elapsed = std::chrono::duration<double>(Clock::now() - begin).count();
            mpf_clears(centerRe, centerIm, scale, (mpf_ptr)0);
            return okay;
        };

        std::vector<double> scalarTimes;
        std::vector<double> batchTimes;
        std::vector<double> fixedBatchTimes;
        std::vector<double> jitTimes;
        std::vector<double> defaultTimes;
        std::vector<double> specializedTimes;
        std::vector<double> specializationRatios;
        std::vector<float> scalar;
        std::vector<float> batch;
        std::vector<float> fixedBatch;
        std::vector<float> native;
        std::vector<float> defaultOutput;
        std::vector<float> specializedOutput;
        for (int sample = 0; sample < SAMPLES; ++sample) {
            double scalarTime = 0.0;
            double batchTime = 0.0;
            double jitTime = 0.0;
            double specializedTime = 0.0;
            if (!render(program, false, nullptr, scalar, scalarTime)) {
                ++failures;
                break;
            }
            scalarTimes.push_back(scalarTime);
            if (program.batchCompatible()) {
                if (!program.avx2Compatible()) {
                    double fixedBatchTime = 0.0;
                    bool rendered = (sample & 1) ? render(program, true, nullptr, fixedBatch, fixedBatchTime, false) && render(program, true, nullptr, batch, batchTime, true) : render(program, true, nullptr, batch, batchTime, true) && render(program, true, nullptr, fixedBatch, fixedBatchTime, false);
                    if (!rendered) {
                        ++failures;
                        break;
                    }
                    fixedBatchTimes.push_back(fixedBatchTime);
                } else if (!render(program, true, nullptr, batch, batchTime)) {
                    ++failures;
                    break;
                }
                batchTimes.push_back(batchTime);
            }
            if (jitAvailable) {
                if (!render(program, true, &jit, native, jitTime)) {
                    ++failures;
                    break;
                }
                jitTimes.push_back(jitTime);
            }
            double defaultTime = 0.0;
            bool paired = true;
            for (int repeat = 0; repeat < PAIR_REPEATS && paired; ++repeat) {
                double defaultOnce = 0.0;
                double specializedOnce = 0.0;
                auto renderDefault = [&]() { return render(specializedProgram, true, specializedJitAvailable ? &specializedJit : nullptr, defaultOutput, defaultOnce); };
                auto renderSpecialized = [&]() { return render(specializedProgram, true, orbitPlan.profitable() && orbitJitAvailable ? &orbitJit : (specializedJitAvailable ? &specializedJit : nullptr), specializedOutput, specializedOnce, true, orbitPlan.profitable() ? &orbitPlan : nullptr); };
                if (!orbitPlan.profitable()) {
                    paired = renderDefault();
                    specializedOnce = defaultOnce;
                    specializedOutput = defaultOutput;
                } else {
                    paired = ((sample + repeat) & 1) ? renderSpecialized() && renderDefault() : renderDefault() && renderSpecialized();
                }
                defaultTime += defaultOnce;
                specializedTime += specializedOnce;
            }
            if (!paired) {
                ++failures;
                break;
            }
            defaultTime /= PAIR_REPEATS;
            specializedTime /= PAIR_REPEATS;
            defaultTimes.push_back(defaultTime);
            specializedTimes.push_back(specializedTime);
            specializationRatios.push_back(defaultTime > 0.0 ? specializedTime / defaultTime : 0.0);
        }
        _putenv_s("MANDEL_EXPR_VECTOR", "");
        _putenv_s("MANDEL_EXPR_HYBRID_REFILL", "");
        _putenv_s("MANDEL_EXPR_VECTOR_TRANSCENDENTALS", "");
        auto median = [](std::vector<double>& values) {
            std::sort(values.begin(), values.end());
            return values.empty() ? 0.0 : values[values.size() / 2];
        };
        double scalarTime = median(scalarTimes);
        double batchTime = median(batchTimes);
        double fixedBatchTime = median(fixedBatchTimes);
        double jitTime = median(jitTimes);
        double defaultTime = median(defaultTimes);
        double specializedTime = median(specializedTimes);
        double specializationRatio = median(specializationRatios);
        int batchMismatches = 0, jitMismatches = 0;
        int refillMismatches = 0, specializedMismatches = 0;
        bool batchExact = !batch.empty() && batch.size() == scalar.size() && std::memcmp(scalar.data(), batch.data(), scalar.size() * sizeof(float)) == 0;
        if (!batchExact) {
            size_t compared = std::min(scalar.size(), batch.size());
            batchMismatches = (int)(std::max(scalar.size(), batch.size()) - compared);
            for (size_t i = 0; i < compared; ++i)
                if (std::memcmp(&scalar[i], &batch[i], sizeof(float)) != 0) ++batchMismatches;
        }
        bool jitExact = native.empty() || (native.size() == scalar.size() && std::memcmp(scalar.data(), native.data(), scalar.size() * sizeof(float)) == 0);
        if (!jitExact) {
            size_t compared = std::min(scalar.size(), native.size());
            jitMismatches = (int)(std::max(scalar.size(), native.size()) - compared);
            for (size_t i = 0; i < compared; ++i)
                if (std::memcmp(&scalar[i], &native[i], sizeof(float)) != 0) ++jitMismatches;
        }
        bool refillExact = fixedBatch.empty() || (fixedBatch.size() == batch.size() && std::memcmp(fixedBatch.data(), batch.data(), batch.size() * sizeof(float)) == 0);
        if (!refillExact) {
            size_t compared = std::min(fixedBatch.size(), batch.size());
            refillMismatches = (int)(std::max(fixedBatch.size(), batch.size()) - compared);
            for (size_t i = 0; i < compared; ++i)
                if (std::memcmp(&fixedBatch[i], &batch[i], sizeof(float)) != 0) ++refillMismatches;
        }
        bool specializedExact = specializedOutput.size() == scalar.size() && std::memcmp(scalar.data(), specializedOutput.data(), scalar.size() * sizeof(float)) == 0;
        if (!specializedExact) {
            size_t compared = std::min(specializedOutput.size(), scalar.size());
            specializedMismatches = (int)(std::max(specializedOutput.size(), scalar.size()) - compared);
            for (size_t i = 0; i < compared; ++i) {
                if (std::memcmp(&scalar[i], &specializedOutput[i], sizeof(float)) != 0) ++specializedMismatches;
            }
        }
        struct Accuracy {
            int classMismatch = 0;
            double maxDifference = 0.0;
            double meanDifference = 0.0;
            double p99Difference = 0.0;
        };
        auto measureAccuracy = [&](const std::vector<float>& candidate) {
            Accuracy accuracy;
            double sumDifference = 0.0;
            std::vector<double> differences;
            const size_t count = std::min(candidate.size(), scalar.size());
            for (size_t index = 0; index < count; ++index) {
                const bool scalarInterior = isInterior(scalar[index]);
                const bool candidateInterior = isInterior(candidate[index]);
                if (scalarInterior != candidateInterior) {
                    ++accuracy.classMismatch;
                } else if (!scalarInterior) {
                    const double difference = std::fabs(static_cast<double>(candidate[index]) - scalar[index]);
                    accuracy.maxDifference = std::max(accuracy.maxDifference, difference);
                    sumDifference += difference;
                    differences.push_back(difference);
                }
            }
            accuracy.classMismatch += static_cast<int>(std::max(candidate.size(), scalar.size()) - count);
            accuracy.meanDifference = differences.empty() ? 0.0 : sumDifference / differences.size();
            if (!differences.empty()) {
                const size_t percentile = static_cast<size_t>(std::ceil(0.99 * differences.size())) - 1;
                std::nth_element(differences.begin(), differences.begin() + percentile, differences.end());
                accuracy.p99Difference = differences[percentile];
            }
            return accuracy;
        };
        const Accuracy batchAccuracy = measureAccuracy(batch);
        const Accuracy specializedAccuracy = measureAccuracy(specializedOutput);
        const bool trigonometricCase = std::strstr(test.source, "sin(") || std::strstr(test.source, "cos(") || std::strstr(test.source, "tan(");
        // This profile compares two chaotic binary64 trajectories, so it is a
        // coarse drift guard. The persisted 512-bit MPFR suite below supplies
        // the strict adoption gate for these same generic/nested opcodes.
        const auto qualityAccepted = [&](const Accuracy& accuracy) { return vectorTranscendentalsEnabled && trigonometricCase && accuracy.classMismatch <= 20 && accuracy.maxDifference <= 150.0 && accuracy.meanDifference <= 0.1 && accuracy.p99Difference == 0.0; };
        const bool batchAccepted = batchExact || qualityAccepted(batchAccuracy);
        const bool specializedAccepted = specializedExact || qualityAccepted(specializedAccuracy);
        if (!batchAccepted || !jitExact || !refillExact || !specializedAccepted) ++failures;
        bool noRegression = specializationRatio > 0.0 && specializationRatio <= 1.02;
        bool improved = specializationRatio > 0.0 && specializationRatio < 1.0;
        if (test.orbitStage) {
            bool classificationOkay = orbitPlan.profitable() == test.expectOrbitInvariant;
            if (std::strcmp(test.name, "orbit-sin-c-cse") == 0 && (orbitPlan.invariantCount() < 2 || orbitPlan.invariantOperationCount() >= 4)) classificationOkay = false;
            if (!classificationOkay) ++failures;
            if (test.expectOrbitInvariant) {
                if (specializedExact) ++orbitInvariantExact;
                if (specializedExact && improved) ++orbitInvariantImproved;
                if (specializedExact && noRegression) ++orbitInvariantNoRegression;
            } else {
                if (specializedExact) ++orbitControlExact;
                if (specializedExact && noRegression) ++orbitControlNoRegression;
            }
        } else {
            if (orbitPlan.profitable()) ++failures;
            if (specializedExact) ++unaffectedOrbitExact;
            if (specializedAccepted) ++unaffectedOrbitQuality;
            if (specializedAccepted && noRegression) ++unaffectedOrbitNoRegression;
        }
        bool acceptanceCase = std::strcmp(test.name, "sine") == 0 || std::strcmp(test.name, "branch-power") == 0;
        if (acceptanceCase && program.batchCompatible() && batchExact && batchTime > 0.0) {
            ++exactHybridCases;
            if (batchTime <= scalarTime * 1.02) ++hybridNoRegressionCases;
            if (fixedBatchTime > 0.0 && batchTime <= fixedBatchTime * 1.02) ++refillNoRegressionCases;
        }
        printf("  %-20s ops=%zu->%zu(-%zu) scalar=%.3fs", test.name, program.instructionCount(), specializedProgram.instructionCount(), program.instructionCount() - specializedProgram.instructionCount(), scalarTime);
        if (batchTime > 0.0)
            printf(" %s=%.3fs(%.2fx)", program.avx2Compatible() ? "AVX2" : "Hybrid", batchTime, scalarTime / batchTime);
        else
            printf(" batch=n/a");
        if (fixedBatchTime > 0.0) printf(" refill=%.2fx", fixedBatchTime / batchTime);
        if (jitTime > 0.0)
            printf(" JIT=%.3fs(%.2fx)", jitTime, scalarTime / jitTime);
        else
            printf(" JIT=n/a");
        printf(" default=%.4fs orbit-%s=%.4fs(%.2fx)", defaultTime, (orbitJitAvailable || specializedJitAvailable) ? "JIT" : (specializedProgram.avx2Compatible() ? "AVX2" : "Hybrid"), specializedTime, specializationRatio > 0.0 ? 1.0 / specializationRatio : 0.0);
        printf(" memcmp=%s/%s/%s/%s mismatch=%d/%d/%d/%d quality=%d/%.0f/%.3g/%.0f\n", batchExact ? "exact" : (batchAccepted ? "quality" : "FAIL"), jitExact ? "exact" : "FAIL", refillExact ? "exact" : "FAIL", specializedExact ? "exact" : (specializedAccepted ? "quality" : "FAIL"), batchMismatches, jitMismatches, refillMismatches, specializedMismatches, specializedAccuracy.classMismatch, specializedAccuracy.maxDifference, specializedAccuracy.meanDifference, specializedAccuracy.p99Difference);
        if (test.orbitStage) printf("    orbit invariants=%zu invariant-ops=%zu body-ops=%zu profitable=%d\n", orbitPlan.invariantCount(), orbitPlan.invariantOperationCount(), orbitPlan.bodyOperationCount(), orbitPlan.profitable() ? 1 : 0);
    }
    bool hybridDefaultGate = exactHybridCases == 2 && hybridNoRegressionCases == 2;
    if (!hybridDefaultGate) ++failures;
    printf("  hybrid default gate: exact=%d/2 non-regression=%d/2\n", exactHybridCases, hybridNoRegressionCases);
    if (refillNoRegressionCases != 2) ++failures;
    printf("  hybrid refill non-regression=%d/2\n", refillNoRegressionCases);
    if (unaffectedOrbitExact < 8 || unaffectedOrbitQuality != 10 || unaffectedOrbitNoRegression != 10) ++failures;
    printf("  unaffected orbit-plan gate: exact=%d/10 quality=%d/10 non-regression=%d/10\n", unaffectedOrbitExact, unaffectedOrbitQuality, unaffectedOrbitNoRegression);
    if (orbitInvariantExact != 2 || (orbitInvariantImproved < 1 && orbitInvariantNoRegression != 2) || orbitControlExact != 1 || orbitControlNoRegression != 1) ++failures;
    printf("  orbit-plan gate: invariant exact=%d/2 improved=%d/2 non-regression=%d/2; control exact=%d/1 non-regression=%d/1\n", orbitInvariantExact, orbitInvariantImproved, orbitInvariantNoRegression, orbitControlExact, orbitControlNoRegression);
    {
        if (VERIFY_JIT) {
            formula::ExpressionProgram identity;
            formula::ExpressionError error;
            formula::ExpressionJit4 jit;
            formula::ExpressionContext contexts[4]{};
            float results[4]{};
            std::atomic_bool halt{true};
            bool cancelled = identity.compile("z", &error) && jit.compile(identity) && !jit.evaluateOrbit(contexts, 3, 10000000, 4.0, results, &halt);
            if (!cancelled) ++failures;
            printf("  orbit cancellation=%s\n", cancelled ? "PASS" : "FAIL");
        } else {
            printf("  orbit cancellation=SKIP (non-JIT verifier)\n");
        }
    }
    printf("  => %s\n\n", failures == 0 ? "PASS" : "CHECK (generic profile mismatch)");
    return failures == 0 ? 0 : 1;
}

static int runGenericPeriodicityCase() {
    constexpr int W = 322;
    constexpr int H = 216;
    constexpr int MXIT = 2000;
    formula::ExpressionProgram canonical;
    formula::ExpressionProgram runtime;
    formula::ExpressionContext fixed;
    formula::ExpressionError error;
    const char* source = "sin(z)+c+complex(1e-100,-1e-101)";
    if (!canonical.compile(source, &error) || !canonical.specialize(fixed, FormulaParameter::C, runtime, &error) || runtime.iterationDependent()) return 1;

    mpf_t centerReal, centerImaginary, scale;
    mpf_init_set_d(centerReal, -0.5);
    mpf_init_set_ui(centerImaginary, 0);
    mpf_init_set_ui(scale, 1);
    auto render = [&](bool periodic, std::vector<float>& output, double& seconds, uint64_t& periodicPixels, uint64_t& iterations) {
        output.assign(static_cast<size_t>(W) * H, EMPTYPIXEL);
        Mandel renderer(W, H, MXIT, 1, output.data());
        _putenv_s("MANDEL_EXPR_PERIODIC", periodic ? "1" : "0");
        const Clock::time_point start = Clock::now();
        const bool okay = renderer.ComputeExpression(centerReal, centerImaginary, scale, runtime, fixed, FormulaParameter::C, MXIT, 100.0, formula::ExpressionColoring::Raw);
        seconds = since(start);
        periodicPixels = renderer.expressionPeriodicPixelCount();
        iterations = renderer.expressionIterationCount();
        return okay;
    };

    std::vector<float> baseline;
    std::vector<float> periodic;
    double baselineSeconds = 0.0;
    double periodicSeconds = 0.0;
    uint64_t baselinePeriodicPixels = 0;
    uint64_t periodicPixels = 0;
    uint64_t baselineIterations = 0;
    uint64_t periodicIterations = 0;
    bool okay = render(false, baseline, baselineSeconds, baselinePeriodicPixels, baselineIterations) && render(true, periodic, periodicSeconds, periodicPixels, periodicIterations);

    formula::ExpressionProgram iterationCanonical;
    formula::ExpressionProgram iterationRuntime;
    const bool iterationSetup = iterationCanonical.compile("sin(z)+c+0*n", &error) && iterationCanonical.specialize(fixed, FormulaParameter::C, iterationRuntime, &error) && iterationRuntime.iterationDependent();
    std::vector<float> iterationOutput(static_cast<size_t>(W) * H, EMPTYPIXEL);
    Mandel iterationRenderer(W, H, 256, 1, iterationOutput.data());
    _putenv_s("MANDEL_EXPR_PERIODIC", "1");
    const bool iterationOkay = iterationSetup && iterationRenderer.ComputeExpression(centerReal, centerImaginary, scale, iterationRuntime, fixed, FormulaParameter::C, 256, 100.0, formula::ExpressionColoring::Raw) && iterationRenderer.expressionPeriodicPixelCount() == 0;

    formula::ExpressionProgram deepCanonical;
    formula::ExpressionProgram deepRuntime;
    std::vector<float> deepPeriodic(35, formula::ExpressionDeepEmptyPixel);
    std::vector<float> deepBaseline(35, formula::ExpressionDeepEmptyPixel);
    const bool deepSetup = deepCanonical.compile("z+0*c", &error) && deepCanonical.specialize(fixed, FormulaParameter::C, deepRuntime, &error) && !deepRuntime.iterationDependent();
    auto renderDeep = [&](bool enable, std::vector<float>& output, formula::ExpressionDeepRenderResult& result) {
        formula::ExpressionDeepRenderRequest request;
        request.canonicalProgram = &deepCanonical;
        request.runtimeProgram = &deepRuntime;
        request.center.realDecimal = "0";
        request.center.imaginaryDecimal = "0";
        request.scale.decimal = "1e500";
        request.fixed = fixed;
        request.pixelParameter = FormulaParameter::C;
        request.width = 7;
        request.height = 5;
        request.maxIterations = 4096;
        request.bailout = 100.0;
        request.output = output.data();
        request.outputCount = output.size();
        request.forceMpfrFallbackForVerification = true;
        request.disableSpecializedPiecewiseMpfrForVerification = true;
        _putenv_s("MANDEL_EXPR_PERIODIC", enable ? "1" : "0");
        return formula::renderExpressionDeepFrame(request, result);
    };
    formula::ExpressionDeepRenderResult deepPeriodicResult;
    formula::ExpressionDeepRenderResult deepBaselineResult;
    const bool deepOkay = deepSetup && renderDeep(true, deepPeriodic, deepPeriodicResult) && renderDeep(false, deepBaseline, deepBaselineResult) && deepPeriodic == deepBaseline && deepPeriodicResult.genericMpfrPeriodicPixelCount == deepPeriodic.size() && deepBaselineResult.genericMpfrPeriodicPixelCount == 0 && deepPeriodicResult.totalIterations < deepBaselineResult.totalIterations;

    formula::ExpressionProgram signedZeroCanonical;
    formula::ExpressionProgram signedZeroRuntime;
    std::vector<float> signedZeroPeriodic(9, formula::ExpressionDeepEmptyPixel);
    std::vector<float> signedZeroBaseline(9, formula::ExpressionDeepEmptyPixel);
    const bool signedZeroSetup = signedZeroCanonical.compile("-arg(z)", &error) && signedZeroCanonical.specialize(fixed, FormulaParameter::C, signedZeroRuntime, &error);
    auto renderSignedZero = [&](bool enable, std::vector<float>& output, formula::ExpressionDeepRenderResult& result) {
        formula::ExpressionDeepRenderRequest request;
        request.canonicalProgram = &signedZeroCanonical;
        request.runtimeProgram = &signedZeroRuntime;
        request.center.realDecimal = "0";
        request.center.imaginaryDecimal = "0";
        request.scale.decimal = "1e500";
        request.fixed = fixed;
        request.pixelParameter = FormulaParameter::C;
        request.width = 3;
        request.height = 3;
        request.maxIterations = 4096;
        request.bailout = 2.0;
        request.output = output.data();
        request.outputCount = output.size();
        request.forceMpfrFallbackForVerification = true;
        request.disableSpecializedPiecewiseMpfrForVerification = true;
        _putenv_s("MANDEL_EXPR_PERIODIC", enable ? "1" : "0");
        return formula::renderExpressionDeepFrame(request, result);
    };
    formula::ExpressionDeepRenderResult signedZeroPeriodicResult;
    formula::ExpressionDeepRenderResult signedZeroBaselineResult;
    const bool signedZeroOkay = signedZeroSetup && renderSignedZero(true, signedZeroPeriodic, signedZeroPeriodicResult) && renderSignedZero(false, signedZeroBaseline, signedZeroBaselineResult) && signedZeroPeriodic == signedZeroBaseline && std::all_of(signedZeroPeriodic.begin(), signedZeroPeriodic.end(), [](float value) { return value == 2.0f; }) && signedZeroPeriodicResult.genericMpfrPeriodicPixelCount == 0;
    _putenv_s("MANDEL_EXPR_PERIODIC", "");
    mpf_clears(centerReal, centerImaginary, scale, (mpf_ptr)0);

    const double speedup = periodicSeconds > 0.0 ? baselineSeconds / periodicSeconds : 0.0;
    okay = okay && iterationOkay && deepOkay && signedZeroOkay && baseline == periodic && baselinePeriodicPixels == 0 && periodicPixels > 0 && periodicIterations < baselineIterations;
    printf("=== generic exact periodicity %dx%d mxit=%d\n", W, H, MXIT);
    printf("  off/on %.3f/%.3f s speedup %.2fx periodic=%llu iterations=%llu/%llu\n", baselineSeconds, periodicSeconds, speedup, (unsigned long long)periodicPixels, (unsigned long long)baselineIterations, (unsigned long long)periodicIterations);
    printf("  output=%s iteration-dependent=%s\n", baseline == periodic ? "exact" : "MISMATCH", iterationOkay ? "disabled" : "FAIL");
    printf("  deep MPFR output=%s periodic=%llu iterations=%llu/%llu\n", deepPeriodic == deepBaseline ? "exact" : "MISMATCH", (unsigned long long)deepPeriodicResult.genericMpfrPeriodicPixelCount, (unsigned long long)deepBaselineResult.totalIterations, (unsigned long long)deepPeriodicResult.totalIterations);
    printf("  deep signed-zero=%s periodic=%llu\n", signedZeroOkay ? "exact" : "FAIL", (unsigned long long)signedZeroPeriodicResult.genericMpfrPeriodicPixelCount);
    printf("  => %s\n\n", okay ? "PASS" : "CHECK (generic periodicity failure)");
    return okay ? 0 : 1;
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

static bool renderExpressionOracle(const FormulaRegressionCase& test, const formula::ExpressionProgram& program, std::vector<float>& output) {
    using formula::ExpressionOracle;
    using formula::ExpressionOracleContext;
    using formula::MpfrComplex;
    output.assign((size_t)test.width * test.height, -2.0f);
    mpfr_prec_t precision = 512;
    if (const char* value = getenv("MANDEL_FORMULA_ORACLE_BITS")) precision = std::max<mpfr_prec_t>(128, (mpfr_prec_t)atol(value));
    const double halfWidth = 2.0 / test.scale;
    const double halfHeight = halfWidth * test.height / test.width;
    const double dx = 2.0 * halfWidth / (test.width - 1);
    const double dy = 2.0 * halfHeight / (test.height - 1);
    mpfr_t magnitude;
    mpfr_init2(magnitude, precision);
    for (int y = 0; y < test.height; ++y) {
        for (int x = 0; x < test.width; ++x) {
            formula::Complex pixel{test.centerRe - halfWidth + dx * x, test.centerIm - halfHeight + dy * y};
            ExpressionOracleContext context(precision);
            context.z0.set(test.fixedZ0.real(), test.fixedZ0.imag());
            context.c.set(test.fixedC.real(), test.fixedC.imag());
            for (int p = 0; p < 8; ++p) context.parameters[p].set(test.parameters[p].real(), test.parameters[p].imag());
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
                    fprintf(stderr, "Formula oracle domain error in %s @ pixel %d,%d iter %d: %s\n", test.name, x, y, n, error.c_str());
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

static int runFormulaRegressionCase(const FormulaRegressionCase& test, bool updateGoldens) {
    formula::ExpressionProgram program;
    formula::ExpressionProgram runtimeProgram;
    formula::ExpressionOrbitPlan orbitPlan;
    formula::ExpressionError compileError;
    if (!program.compile(test.source, &compileError)) {
        printf("=== formula %s\n  compile error @ %zu: %s\n  => CHECK\n\n", test.name, compileError.position, compileError.message.c_str());
        return 1;
    }

    formula::ExpressionContext fixed;
    fixed.z0 = test.fixedZ0;
    fixed.c = test.fixedC;
    fixed.parameters = test.parameters;
    if (!program.specialize(fixed, test.pixel, runtimeProgram, &compileError) || !orbitPlan.build(runtimeProgram, &compileError)) return 1;
    formula::ExpressionJit4 jit;
    const formula::ExpressionJit4* activeJit = nullptr;
    if (runtimeProgram.fastPath() == formula::ExpressionProgram::FastPath::None) {
        bool compiled = orbitPlan.profitable() ? jit.compile(orbitPlan) : jit.compile(runtimeProgram);
        if (compiled) activeJit = &jit;
    }
    std::vector<float> engine((size_t)test.width * test.height, EMPTYPIXEL);
    std::vector<float> golden(engine.size(), EMPTYPIXEL);
    Mandel renderer(test.width, test.height, test.mxit, 1, engine.data());
    mpf_t centerRe, centerIm, scale;
    mpf_init_set_d(centerRe, test.centerRe);
    mpf_init_set_d(centerIm, test.centerIm);
    mpf_init_set_d(scale, test.scale);
    auto begin = Clock::now();
    bool rendered = renderer.ComputeExpression(centerRe, centerIm, scale, runtimeProgram, fixed, test.pixel, test.mxit, test.bailout, formula::ExpressionColoring::Raw, activeJit, orbitPlan.profitable() ? &orbitPlan : nullptr);
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
    int diagnosticCount = 0;
    const bool diagnostics = getenv("MANDEL_FORMULA_DIAGNOSTICS") != nullptr;
    for (size_t i = 0; i < engine.size(); ++i) {
        if (engine[i] == EMPTYPIXEL || golden[i] == EMPTYPIXEL) {
            ++emptyPixels;
            continue;
        }
        bool engineInterior = isInterior(engine[i]);
        bool goldenInterior = isInterior(golden[i]);
        if (engineInterior != goldenInterior) {
            ++classMismatch;
            if (diagnostics && diagnosticCount++ < 16) printf("  mismatch pixel=%zu (%zu,%zu) engine=%.9g GT=%.9g class\n", i, i % static_cast<size_t>(test.width), i / static_cast<size_t>(test.width), engine[i], golden[i]);
            continue;
        }
        if (!engineInterior) {
            ++exteriorBoth;
            double difference = std::fabs((double)engine[i] - golden[i]);
            maxDiff = std::max(maxDiff, difference);
            sumDiff += difference;
            differences.push_back(difference);
            if (diagnostics && difference != 0.0 && diagnosticCount++ < 16) printf("  mismatch pixel=%zu (%zu,%zu) engine=%.9g GT=%.9g diff=%.9g\n", i, i % static_cast<size_t>(test.width), i / static_cast<size_t>(test.width), engine[i], golden[i], difference);
        }
    }
    double meanDiff = exteriorBoth ? sumDiff / exteriorBoth : 0.0;
    double p99Diff = 0.0;
    if (!differences.empty()) {
        size_t index = (size_t)std::ceil(0.99 * differences.size()) - 1;
        std::nth_element(differences.begin(), differences.begin() + index, differences.end());
        p99Diff = differences[index];
    }
    bool ok = emptyPixels == 0 && classMismatch <= test.maxClassMismatch && maxDiff <= test.maxDiff && meanDiff <= test.maxMeanDiff && p99Diff <= test.maxP99Diff;
    printf("=== formula %s  (%dx%d, mxit=%d, %s-plane)\n", test.name, test.width, test.height, test.mxit, test.pixel == FormulaParameter::C ? "c" : "z0");
    printf("  expression: %s\n", test.source);
    printf("  engine %.3fs%s", engineTime, updateGoldens ? "  oracle " : "  golden load");
    if (updateGoldens) printf("%.3fs", oracleTime);
    printf("\n  empty=%d class mismatch=%d  max=%.6g mean=%.6g p99=%.6g"
           "  checksum=0x%08x\n",
           emptyPixels, classMismatch, maxDiff, meanDiff, p99Diff, checksum(engine.data(), (int)engine.size()));
    printf("  => %s\n\n", ok ? "PASS" : "CHECK (formula golden mismatch)");
    return ok ? 0 : 1;
}

static int runFormulaRegressionSuite() {
    std::vector<FormulaRegressionCase> cases;
    FormulaRegressionCase quadratic{"quadratic-c", "z*z+c", FormulaParameter::C, -0.5, 0.0, 1.0, {}, {}, {}, 4.0, 500, 64, 44, "tests/golden/formula_quadratic_c.f32"};
    cases.push_back(quadratic);
    FormulaRegressionCase julia = quadratic;
    julia.name = "quadratic-z0";
    julia.pixel = FormulaParameter::InitialZ;
    julia.centerRe = 0.0;
    julia.fixedC = {-0.8, 0.156};
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
    FormulaRegressionCase genericSine = sine;
    genericSine.name = "generic-sine-constant";
    genericSine.source = "sin(z)+c+complex(1e-100,-1e-101)";
    genericSine.centerRe = -0.5;
    genericSine.bailout = 100.0;
    genericSine.mxit = 2000;
    genericSine.goldenPath = "tests/golden/formula_generic_sine_constant.f32";
    genericSine.maxClassMismatch = 1;
    genericSine.maxDiff = 100.0;
    genericSine.maxMeanDiff = 0.3;
    cases.push_back(genericSine);
    FormulaRegressionCase nestedSine = genericSine;
    nestedSine.name = "generic-nested-sine";
    nestedSine.source = "sin(sin(z))+c+complex(1e-100,-1e-101)";
    nestedSine.mxit = 500;
    nestedSine.goldenPath = "tests/golden/formula_generic_nested_sine.f32";
    nestedSine.maxClassMismatch = 0;
    nestedSine.maxDiff = 50.0;
    nestedSine.maxMeanDiff = 0.1;
    cases.push_back(nestedSine);
    FormulaRegressionCase nestedCosine = genericSine;
    nestedCosine.name = "generic-nested-cosine";
    nestedCosine.source = "cos(cos(z))+c+complex(1e-100,-1e-101)";
    nestedCosine.mxit = 500;
    nestedCosine.goldenPath = "tests/golden/formula_generic_nested_cosine.f32";
    nestedCosine.maxClassMismatch = 0;
    nestedCosine.maxDiff = 50.0;
    nestedCosine.maxMeanDiff = 0.1;
    cases.push_back(nestedCosine);
    FormulaRegressionCase nestedTangent = genericSine;
    nestedTangent.name = "generic-nested-tangent";
    nestedTangent.source = "tan(tan(z))+c+complex(1e-100,-1e-101)";
    nestedTangent.mxit = 500;
    nestedTangent.goldenPath = "tests/golden/formula_generic_nested_tangent.f32";
    nestedTangent.maxClassMismatch = 0;
    nestedTangent.maxDiff = 0.0;
    nestedTangent.maxMeanDiff = 0.0;
    cases.push_back(nestedTangent);
    FormulaRegressionCase parameter = quadratic;
    parameter.name = "parameter-polynomial";
    parameter.source = "z*z+c+p0*z";
    parameter.parameters[0] = {0.15, -0.05};
    parameter.goldenPath = "tests/golden/formula_parameter_poly.f32";
    cases.push_back(parameter);
    FormulaRegressionCase branch = quadratic;
    branch.name = "branch-power";
    branch.source = "exp(p0*log(z))+c";
    branch.fixedZ0 = {0.3, 0.2};
    branch.parameters[0] = {2.5, 0.0};
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
    const char* selectedCase = getenv("MANDEL_FORMULA_CASE");
    int result = 0;
    for (const FormulaRegressionCase& test : cases) {
        if (!selectedCase || std::strcmp(test.name, selectedCase) == 0) result |= runFormulaRegressionCase(test, update);
    }
    return result;
}

static std::string integerPowerSource(int power) {
    std::string source = "z";
    for (int i = 1; i < power; ++i) source += "*z";
    source += "+c";
    return source;
}

static bool renderIntegerPowerFrame(int power, FormulaParameter pixel, formula::ExpressionColoring coloring, int width, int height, int mxit, bool simd, std::vector<float>& output, double* elapsed = nullptr) {
    formula::ExpressionProgram program;
    formula::ExpressionError error;
    if (!program.compile(integerPowerSource(power), &error) || program.fastIntegerPower() != power) return false;

    formula::ExpressionContext fixed;
    if (pixel == FormulaParameter::InitialZ) fixed.c = {-0.35, 0.62};
    output.assign((size_t)width * height, EMPTYPIXEL);
    Mandel renderer(width, height, mxit, 1, output.data());
    mpf_t centerRe, centerIm, scale;
    mpf_init_set_d(centerRe, power == 2 && pixel == FormulaParameter::C ? -0.5 : 0.0);
    mpf_init_set_ui(centerIm, 0);
    mpf_init_set_ui(scale, 1);
    _putenv_s("MANDEL_EXPR_POWER_SIMD", simd ? "1" : "0");
    auto begin = Clock::now();
    bool okay = renderer.ComputeExpression(centerRe, centerIm, scale, program, fixed, pixel, mxit, 4.0, coloring);
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
    const formula::ExpressionColoring colorings[] = {formula::ExpressionColoring::Raw, formula::ExpressionColoring::Smooth, formula::ExpressionColoring::Distance, formula::ExpressionColoring::Feather, formula::ExpressionColoring::OrbitTrap};
    const char* coloringNames[] = {"raw", "smooth", "distance", "feather", "trap"};

    printf("=== integer Multibrot AVX2 correctness\n");
    for (int power = 2; power <= 8; ++power) {
        for (FormulaParameter pixel : {FormulaParameter::C, FormulaParameter::InitialZ}) {
            for (int coloringIndex = 0; coloringIndex < (int)std::size(colorings); ++coloringIndex) {
                std::vector<float> scalar;
                std::vector<float> simd;
                if (!renderIntegerPowerFrame(power, pixel, colorings[coloringIndex], W, H, MXIT, false, scalar) || !renderIntegerPowerFrame(power, pixel, colorings[coloringIndex], W, H, MXIT, true, simd)) {
                    ++failures;
                    continue;
                }
                int bitMismatches = 0;
                int classMismatches = 0;
                double maxDifference = 0.0;
                for (size_t i = 0; i < scalar.size(); ++i) {
                    if (std::memcmp(&scalar[i], &simd[i], sizeof(float)) != 0) ++bitMismatches;
                    if (isInterior(scalar[i]) != isInterior(simd[i]))
                        ++classMismatches;
                    else if (!isInterior(scalar[i]))
                        maxDifference = std::max(maxDifference, std::fabs((double)scalar[i] - simd[i]));
                }
                totalBitMismatches += bitMismatches;
                totalClassMismatches += classMismatches;
                maximumDifference = std::max(maximumDifference, maxDifference);
                if (bitMismatches || classMismatches) ++failures;
                printf("  z^%d+c %-2s-plane %-8s bits=%d class=%d max=%.6g\n", power, pixel == FormulaParameter::C ? "c" : "z0", coloringNames[coloringIndex], bitMismatches, classMismatches, maxDifference);
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
        for (FormulaParameter pixel : {FormulaParameter::C, FormulaParameter::InitialZ}) {
            FormulaRegressionCase test{};
            test.name = "multibrot-oracle";
            test.source = source.c_str();
            test.pixel = pixel;
            test.centerRe = 0.0;
            test.centerIm = 0.0;
            test.scale = 1.0;
            if (pixel == FormulaParameter::InitialZ) test.fixedC = {-0.35, 0.62};
            test.bailout = 4.0;
            test.mxit = OMIT;
            test.width = OW;
            test.height = OH;
            std::vector<float> oracle;
            std::vector<float> simd;
            if (!renderExpressionOracle(test, program, oracle) || !renderIntegerPowerFrame(power, pixel, formula::ExpressionColoring::Raw, OW, OH, OMIT, true, simd)) {
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
            printf("    z^%d+c %-2s-plane class=%d iteration=%d\n", power, pixel == FormulaParameter::C ? "c" : "z0", classMismatches, iterationMismatches);
        }
    }
    {
        formula::ExpressionProgram program;
        formula::ExpressionError error;
        program.compile("z*z*z+c", &error);
        FormulaRegressionCase deep{};
        deep.name = "cubic-deep-oracle";
        deep.source = "z*z*z+c";
        deep.pixel = FormulaParameter::C;
        deep.centerRe = 0.3849001794597505;
        deep.centerIm = 0.0;
        deep.scale = 1e10;
        deep.bailout = 4.0;
        deep.mxit = 3000;
        deep.width = 24;
        deep.height = 16;
        std::vector<float> oracle;
        std::vector<float> seriesOutput((size_t)deep.width * deep.height);
        formula::ExpressionContext fixed;
        mpf_t centerRe, centerIm, scale;
        mpf_init_set_d(centerRe, deep.centerRe);
        mpf_init_set_d(centerIm, deep.centerIm);
        mpf_init_set_d(scale, deep.scale);
        Mandel renderer(deep.width, deep.height, deep.mxit, 1, seriesOutput.data());
        bool used = false;
        int series = 0;
        bool okay = renderExpressionOracle(deep, program, oracle) && renderer.ComputeExpressionResidual(centerRe, centerIm, scale, program, fixed, deep.pixel, deep.mxit, deep.bailout, &used, formula::ExpressionColoring::Raw, &series);
        mpf_clears(centerRe, centerIm, scale, (mpf_ptr)0);
        int classMismatches = 0;
        int iterationMismatches = 0;
        if (okay && oracle.size() == seriesOutput.size()) {
            for (size_t i = 0; i < oracle.size(); ++i) {
                bool a = isInterior(oracle[i]);
                bool b = isInterior(seriesOutput[i]);
                if (a != b)
                    ++classMismatches;
                else if (!a && oracle[i] != seriesOutput[i])
                    ++iterationMismatches;
            }
        } else {
            classMismatches = 1;
        }
        if (!used || series <= 0 || classMismatches || iterationMismatches) ++failures;
        printf("    cubic deep 1e10 class=%d iteration=%d SA=%d\n", classMismatches, iterationMismatches, series);
    }

    auto extremeParity = [&](const char* name, const char* centerText, const char* scaleText, double bailout, formula::ExpressionColoring coloring, int mxit) {
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
        bool scalarOkay = scalarRenderer.ComputeExpression(centerRe, centerIm, scale, program, fixed, FormulaParameter::C, mxit, bailout, coloring);
        _putenv_s("MANDEL_EXPR_POWER_SIMD", "1");
        Mandel simdRenderer(EW, EH, mxit, 1, simd.data());
        bool simdOkay = simdRenderer.ComputeExpression(centerRe, centerIm, scale, program, fixed, FormulaParameter::C, mxit, bailout, coloring);
        mpf_clears(centerRe, centerIm, scale, (mpf_ptr)0);
        int mismatches = scalar.size() == simd.size() ? (int)!std::equal(scalar.begin(), scalar.end(), simd.begin(), [](float a, float b) { return std::memcmp(&a, &b, sizeof(float)) == 0; }) : 1;
        if (!scalarOkay || !simdOkay || mismatches) ++failures;
        printf("  extreme %-18s parity=%s\n", name, scalarOkay && simdOkay && !mismatches ? "exact" : "FAIL");
    };
    extremeParity("huge bailout", "1.5e200", "1", 1.0e200, formula::ExpressionColoring::Raw, 1);
    extremeParity("tiny bailout", "0", "1e200", 1.0e-200, formula::ExpressionColoring::Raw, 8);
    extremeParity("subnormal bailout", "1.0000000000000002e-160", "1e170", 1.0e-160, formula::ExpressionColoring::Raw, 1);
    extremeParity("subnormal dx EDE", "2", "1e309", 4.0, formula::ExpressionColoring::Distance, 8);

    constexpr int BW = 960, BH = 540, BMIT = 3000, SAMPLES = 5;
    std::vector<double> scalarTimes;
    std::vector<double> simdTimes;
    std::vector<float> benchmarkOutput;
    for (int sample = 0; sample < SAMPLES; ++sample) {
        double scalar = 0.0, simd = 0.0;
        if (!renderIntegerPowerFrame(3, FormulaParameter::C, formula::ExpressionColoring::Smooth, BW, BH, BMIT, false, benchmarkOutput, &scalar) || !renderIntegerPowerFrame(3, FormulaParameter::C, formula::ExpressionColoring::Smooth, BW, BH, BMIT, true, benchmarkOutput, &simd)) {
            ++failures;
            break;
        }
        scalarTimes.push_back(scalar);
        simdTimes.push_back(simd);
    }
    _putenv_s("MANDEL_EXPR_POWER_SIMD", "");
    std::sort(scalarTimes.begin(), scalarTimes.end());
    std::sort(simdTimes.begin(), simdTimes.end());
    double scalarMedian = scalarTimes.empty() ? 0.0 : scalarTimes[scalarTimes.size() / 2];
    double simdMedian = simdTimes.empty() ? 0.0 : simdTimes[simdTimes.size() / 2];
    double speedup = simdMedian > 0.0 ? scalarMedian / simdMedian : 0.0;
    if (speedup < 1.05) ++failures;

    formula::ExpressionProgram cubicProgram;
    formula::ExpressionError cubicError;
    cubicProgram.compile("z*z*z+c", &cubicError);
    constexpr int RW = 320, RH = 216, RMIT = 3000, RSAMPLES = 3;
    std::vector<double> genericResidualTimes;
    std::vector<double> cubicResidualTimes;
    std::vector<float> genericResidual;
    std::vector<float> cubicResidual;
    auto renderResidual = [&](bool specialized, formula::ExpressionColoring coloring, std::vector<float>& output, double& elapsed, bool seriesEnabled = true, int* seriesIterations = nullptr, double viewScale = 1e10) {
        output.assign((size_t)RW * RH, EMPTYPIXEL);
        Mandel renderer(RW, RH, RMIT, 1, output.data());
        formula::ExpressionContext fixed;
        mpf_t centerRe, centerIm, scale;
        mpf_init_set_d(centerRe, 0.3849001794597505);
        mpf_init_set_ui(centerIm, 0);
        mpf_init_set_d(scale, viewScale);
        _putenv_s("MANDEL_EXPR_RESIDUAL_POWER", specialized ? "1" : "0");
        _putenv_s("MANDEL_EXPR_CUBIC_SA", seriesEnabled ? "1" : "0");
        bool used = false;
        int series = 0;
        auto begin = Clock::now();
        bool okay = renderer.ComputeExpressionResidual(centerRe, centerIm, scale, cubicProgram, fixed, FormulaParameter::C, RMIT, 4.0, &used, coloring, &series);
        elapsed = since(begin);
        mpf_clears(centerRe, centerIm, scale, (mpf_ptr)0);
        if (seriesIterations) *seriesIterations = series;
        return okay && used;
    };
    const formula::ExpressionColoring residualColorings[] = {formula::ExpressionColoring::Raw, formula::ExpressionColoring::Smooth, formula::ExpressionColoring::Distance};
    const char* residualColoringNames[] = {"raw", "smooth", "distance"};
    for (int coloringIndex = 0; coloringIndex < 3; ++coloringIndex) {
        std::vector<float> direct((size_t)RW * RH, EMPTYPIXEL);
        std::vector<float> residual;
        Mandel directRenderer(RW, RH, RMIT, 1, direct.data());
        formula::ExpressionContext fixed;
        mpf_t centerRe, centerIm, scale;
        mpf_init_set_d(centerRe, 0.3849001794597505);
        mpf_init_set_ui(centerIm, 0);
        mpf_init_set_d(scale, 1e10);
        _putenv_s("MANDEL_EXPR_POWER_SIMD", "0");
        bool directOkay = directRenderer.ComputeExpression(centerRe, centerIm, scale, cubicProgram, fixed, FormulaParameter::C, RMIT, 4.0, residualColorings[coloringIndex]);
        mpf_clears(centerRe, centerIm, scale, (mpf_ptr)0);
        double residualTime = 0.0;
        bool residualOkay = renderResidual(true, residualColorings[coloringIndex], residual, residualTime);
        int classMismatches = 0;
        int floorMismatches = 0;
        double maxDifference = 0.0;
        for (size_t i = 0; i < direct.size(); ++i) {
            bool directInterior = isInterior(direct[i]);
            bool residualInterior = isInterior(residual[i]);
            if (directInterior != residualInterior) {
                ++classMismatches;
            } else if (!directInterior) {
                if ((int)direct[i] != (int)residual[i]) ++floorMismatches;
                maxDifference = std::max(maxDifference, std::fabs((double)direct[i] - residual[i]));
            }
        }
        if (!directOkay || !residualOkay || classMismatches || floorMismatches || (coloringIndex == 0 && maxDifference != 0.0) || (coloringIndex != 0 && maxDifference > 1e-3)) ++failures;
        printf("  cubic residual %-8s class=%d floor=%d max=%.6g\n", residualColoringNames[coloringIndex], classMismatches, floorMismatches, maxDifference);
    }
    for (int sample = 0; sample < RSAMPLES; ++sample) {
        double genericTime = 0.0, specializedTime = 0.0;
        if (!renderResidual(false, formula::ExpressionColoring::Raw, genericResidual, genericTime) || !renderResidual(true, formula::ExpressionColoring::Raw, cubicResidual, specializedTime)) {
            ++failures;
            break;
        }
        genericResidualTimes.push_back(genericTime);
        cubicResidualTimes.push_back(specializedTime);
    }
    std::vector<float> cubicNoSeries;
    std::vector<float> cubicWithSeries;
    double cubicNoSeriesTime = 0.0;
    double cubicWithSeriesTime = 0.0;
    int seriesIterations = 0;
    bool noSeriesOkay = renderResidual(true, formula::ExpressionColoring::Raw, cubicNoSeries, cubicNoSeriesTime, false);
    bool withSeriesOkay = renderResidual(true, formula::ExpressionColoring::Raw, cubicWithSeries, cubicWithSeriesTime, true, &seriesIterations);
    int seriesClassMismatches = 0;
    int seriesFloorMismatches = 0;
    if (cubicNoSeries.size() == cubicWithSeries.size()) {
        for (size_t i = 0; i < cubicNoSeries.size(); ++i) {
            bool a = isInterior(cubicNoSeries[i]);
            bool b = isInterior(cubicWithSeries[i]);
            if (a != b)
                ++seriesClassMismatches;
            else if (!a && (int)cubicNoSeries[i] != (int)cubicWithSeries[i])
                ++seriesFloorMismatches;
        }
    } else {
        seriesClassMismatches = 1;
    }
    double seriesSpeedup = cubicWithSeriesTime > 0.0 ? cubicNoSeriesTime / cubicWithSeriesTime : 0.0;
    if (!noSeriesOkay || !withSeriesOkay || seriesIterations <= 0 || seriesClassMismatches || seriesFloorMismatches || seriesSpeedup < 0.95) ++failures;
    std::vector<float> shallowNoSeries;
    std::vector<float> shallowWithSeries;
    double shallowNoSeriesTime = 0.0;
    double shallowWithSeriesTime = 0.0;
    int shallowSeriesIterations = 0;
    bool shallowNoSeriesOkay = renderResidual(true, formula::ExpressionColoring::Raw, shallowNoSeries, shallowNoSeriesTime, false, nullptr, 1.0);
    bool shallowWithSeriesOkay = renderResidual(true, formula::ExpressionColoring::Raw, shallowWithSeries, shallowWithSeriesTime, true, &shallowSeriesIterations, 1.0);
    int shallowMismatches = 0;
    if (shallowNoSeries.size() == shallowWithSeries.size()) {
        for (size_t i = 0; i < shallowNoSeries.size(); ++i) {
            if (std::memcmp(&shallowNoSeries[i], &shallowWithSeries[i], sizeof(float)) != 0) ++shallowMismatches;
        }
    } else {
        shallowMismatches = 1;
    }
    if (!shallowNoSeriesOkay || !shallowWithSeriesOkay || shallowSeriesIterations != 0 || shallowMismatches) ++failures;

    auto renderResidualZ0 = [&](formula::ExpressionColoring coloring, bool seriesEnabled, std::vector<float>& output, int* seriesIterations) {
        constexpr int ZW = 96, ZH = 64, ZMIT = 300;
        output.assign((size_t)ZW * ZH, EMPTYPIXEL);
        Mandel renderer(ZW, ZH, ZMIT, 1, output.data());
        formula::ExpressionContext fixed;
        fixed.c = {0.0, 0.0};
        mpf_t centerRe, centerIm, scale;
        mpf_init_set_ui(centerRe, 1);
        mpf_init_set_ui(centerIm, 0);
        mpf_init_set_d(scale, 1e10);
        _putenv_s("MANDEL_EXPR_CUBIC_SA", seriesEnabled ? "1" : "0");
        bool used = false;
        int series = 0;
        bool okay = renderer.ComputeExpressionResidual(centerRe, centerIm, scale, cubicProgram, fixed, FormulaParameter::InitialZ, ZMIT, 4.0, &used, coloring, &series);
        mpf_clears(centerRe, centerIm, scale, (mpf_ptr)0);
        if (seriesIterations) *seriesIterations = series;
        return okay && used;
    };
    for (formula::ExpressionColoring coloring : {formula::ExpressionColoring::Raw, formula::ExpressionColoring::Distance}) {
        std::vector<float> z0Off;
        std::vector<float> z0On;
        int z0Series = 0;
        bool offOkay = renderResidualZ0(coloring, false, z0Off, nullptr);
        bool onOkay = renderResidualZ0(coloring, true, z0On, &z0Series);
        int classMismatches = 0;
        int floorMismatches = 0;
        double maxDifference = 0.0;
        for (size_t i = 0; i < z0Off.size(); ++i) {
            bool a = isInterior(z0Off[i]);
            bool b = isInterior(z0On[i]);
            if (a != b)
                ++classMismatches;
            else if (!a) {
                if ((int)z0Off[i] != (int)z0On[i]) ++floorMismatches;
                maxDifference = std::max(maxDifference, std::fabs((double)z0Off[i] - z0On[i]));
            }
        }
        double tolerance = coloring == formula::ExpressionColoring::Raw ? 0.0 : 1e-3;
        if (!offOkay || !onOkay || z0Series < 8 || classMismatches || floorMismatches || maxDifference > tolerance) ++failures;
        printf("  cubic SA z0 %-8s iterations=%d class=%d floor=%d"
               " max=%.6g\n",
               coloring == formula::ExpressionColoring::Raw ? "raw" : "distance", z0Series, classMismatches, floorMismatches, maxDifference);
    }
    {
        constexpr int IW = 48, IH = 32, IMIT = 50;
        auto renderInitialEscape = [&](bool seriesEnabled, std::vector<float>& output, int& seriesIterations) {
            output.assign((size_t)IW * IH, EMPTYPIXEL);
            Mandel renderer(IW, IH, IMIT, 1, output.data());
            formula::ExpressionContext fixed;
            fixed.c = {0.099, 0.0};
            mpf_t centerRe, centerIm, scale;
            mpf_init_set_d(centerRe, 0.1);
            mpf_init_set_ui(centerIm, 0);
            mpf_init_set_d(scale, 1e6);
            _putenv_s("MANDEL_EXPR_CUBIC_SA", seriesEnabled ? "1" : "0");
            bool used = false;
            bool okay = renderer.ComputeExpressionResidual(centerRe, centerIm, scale, cubicProgram, fixed, FormulaParameter::InitialZ, IMIT, 0.100001, &used, formula::ExpressionColoring::Raw, &seriesIterations);
            mpf_clears(centerRe, centerIm, scale, (mpf_ptr)0);
            return okay && used;
        };
        std::vector<float> off;
        std::vector<float> on;
        int offSeries = 0, onSeries = 0;
        bool offOkay = renderInitialEscape(false, off, offSeries);
        bool onOkay = renderInitialEscape(true, on, onSeries);
        int mismatches = 0;
        for (size_t i = 0; i < off.size(); ++i) {
            if (std::memcmp(&off[i], &on[i], sizeof(float)) != 0) ++mismatches;
        }
        if (!offOkay || !onOkay || onSeries < 8 || mismatches) ++failures;
        printf("  cubic SA z0 initial-escape iterations=%d bits=%d\n", onSeries, mismatches);
    }
    {
        FormulaRegressionCase z0Deep{};
        z0Deep.name = "cubic-z0-deep-oracle";
        z0Deep.source = "z*z*z+c";
        z0Deep.pixel = FormulaParameter::InitialZ;
        z0Deep.centerRe = 1.0;
        z0Deep.scale = 1e10;
        z0Deep.fixedC = {};
        z0Deep.bailout = 4.0;
        z0Deep.mxit = 300;
        z0Deep.width = 24;
        z0Deep.height = 16;
        std::vector<float> oracle;
        std::vector<float> residual;
        int z0Series = 0;
        bool okay = renderExpressionOracle(z0Deep, cubicProgram, oracle);
        int classMismatches = 0;
        int iterationMismatches = 0;
        if (okay) {
            residual.resize((size_t)z0Deep.width * z0Deep.height);
            Mandel renderer(z0Deep.width, z0Deep.height, z0Deep.mxit, 1, residual.data());
            formula::ExpressionContext fixed;
            mpf_t centerRe, centerIm, scale;
            mpf_init_set_ui(centerRe, 1);
            mpf_init_set_ui(centerIm, 0);
            mpf_init_set_d(scale, 1e10);
            bool used = false;
            okay = renderer.ComputeExpressionResidual(centerRe, centerIm, scale, cubicProgram, fixed, FormulaParameter::InitialZ, z0Deep.mxit, 4.0, &used, formula::ExpressionColoring::Raw, &z0Series) && used;
            mpf_clears(centerRe, centerIm, scale, (mpf_ptr)0);
        }
        if (okay && oracle.size() == residual.size()) {
            for (size_t i = 0; i < oracle.size(); ++i) {
                bool a = isInterior(oracle[i]);
                bool b = isInterior(residual[i]);
                if (a != b)
                    ++classMismatches;
                else if (!a && oracle[i] != residual[i])
                    ++iterationMismatches;
            }
        } else {
            classMismatches = 1;
        }
        if (classMismatches || iterationMismatches) ++failures;
        printf("    cubic z0 deep class=%d iteration=%d SA=%d\n", classMismatches, iterationMismatches, z0Series);
    }
    _putenv_s("MANDEL_EXPR_RESIDUAL_POWER", "");
    _putenv_s("MANDEL_EXPR_CUBIC_SA", "");
    int residualBitMismatches = 0;
    if (genericResidual.size() == cubicResidual.size()) {
        for (size_t i = 0; i < genericResidual.size(); ++i) {
            if (std::memcmp(&genericResidual[i], &cubicResidual[i], sizeof(float)) != 0) ++residualBitMismatches;
        }
    } else {
        residualBitMismatches = 1;
    }
    std::sort(genericResidualTimes.begin(), genericResidualTimes.end());
    std::sort(cubicResidualTimes.begin(), cubicResidualTimes.end());
    double genericResidualMedian = genericResidualTimes.empty() ? 0.0 : genericResidualTimes[genericResidualTimes.size() / 2];
    double cubicResidualMedian = cubicResidualTimes.empty() ? 0.0 : cubicResidualTimes[cubicResidualTimes.size() / 2];
    double residualSpeedup = cubicResidualMedian > 0.0 ? genericResidualMedian / cubicResidualMedian : 0.0;
    if (residualBitMismatches || residualSpeedup < 1.10) ++failures;

    auto renderDirectCusp = [&](double viewScale, std::vector<float>& output, double& elapsed) {
        output.assign((size_t)RW * RH, EMPTYPIXEL);
        Mandel renderer(RW, RH, RMIT, 1, output.data());
        formula::ExpressionContext fixed;
        mpf_t centerRe, centerIm, scale;
        mpf_init_set_d(centerRe, 0.3849001794597505);
        mpf_init_set_ui(centerIm, 0);
        mpf_init_set_d(scale, viewScale);
        _putenv_s("MANDEL_EXPR_POWER_SIMD", "1");
        auto begin = Clock::now();
        bool okay = renderer.ComputeExpression(centerRe, centerIm, scale, cubicProgram, fixed, FormulaParameter::C, RMIT, 4.0, formula::ExpressionColoring::Smooth);
        elapsed = since(begin);
        mpf_clears(centerRe, centerIm, scale, (mpf_ptr)0);
        return okay;
    };
    const double crossoverScales[] = {1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12};
    double firstProfitableScale = 0.0;
    printf("  cubic production crossover (Smooth):\n");
    for (double viewScale : crossoverScales) {
        std::vector<float> direct;
        std::vector<float> residual;
        double directTime = 0.0, residualTime = 0.0;
        int coverage = 0;
        bool directOkay = renderDirectCusp(viewScale, direct, directTime);
        bool residualOkay = renderResidual(true, formula::ExpressionColoring::Smooth, residual, residualTime, true, &coverage, viewScale);
        int classMismatches = 0;
        int floorMismatches = 0;
        double maxDifference = 0.0;
        for (size_t i = 0; i < direct.size(); ++i) {
            bool a = isInterior(direct[i]);
            bool b = isInterior(residual[i]);
            if (a != b)
                ++classMismatches;
            else if (!a) {
                if ((int)direct[i] != (int)residual[i]) ++floorMismatches;
                maxDifference = std::max(maxDifference, std::fabs((double)direct[i] - residual[i]));
            }
        }
        double crossoverSpeedup = residualTime > 0.0 ? directTime / residualTime : 0.0;
        if (!directOkay || !residualOkay || classMismatches || floorMismatches || maxDifference > 1e-3) ++failures;
        if (firstProfitableScale == 0.0 && crossoverSpeedup >= 1.10) firstProfitableScale = viewScale;
        printf("    scale=%.0e SA=%d direct/residual %.3f/%.3f"
               " speedup=%.2fx class=%d floor=%d max=%.3g\n",
               viewScale, coverage, directTime, residualTime, crossoverSpeedup, classMismatches, floorMismatches, maxDifference);
    }
    if (firstProfitableScale == 0.0) ++failures;
    _putenv_s("MANDEL_EXPR_POWER_SIMD", "");

    int backendDeepMismatches = 0;
    int backendShallowMismatches = 0;
    int backendDisabledMismatches = 0;
    int backendZ0Mismatches = 0;
    int backendPartialMismatches = 0;
    int backendColoredMismatches = 0;
    {
        constexpr int PW = 160, PH = 108, PMIT = 3000;
        formula::ExpressionContext fixed;
        std::vector<float> expected((size_t)PW * PH, EMPTYPIXEL);
        std::vector<float> dispatched((size_t)PW * PH, EMPTYPIXEL);
        Mandel expectedRenderer(PW, PH, PMIT, 1, expected.data());
        Mandel dispatchedRenderer(PW, PH, PMIT, 1, dispatched.data());
        std::unique_ptr<IComputeBackend> backend = createComputeBackend("cpu");
        mpf_t centerRe, centerIm, scale;
        mpf_init_set_d(centerRe, 0.3849001794597505);
        mpf_init_set_ui(centerIm, 0);
        mpf_init_set_d(scale, 1e8);
        bool used = false;
        bool expectedOkay = expectedRenderer.ComputeExpressionResidual(centerRe, centerIm, scale, cubicProgram, fixed, FormulaParameter::C, PMIT, 4.0, &used, formula::ExpressionColoring::Smooth);
        ComputeRequest request;
        request.mode = ComputeMode::Expression;
        request.cpuEngine = &dispatchedRenderer;
        request.centerRe = centerRe;
        request.centerIm = centerIm;
        request.scale = scale;
        request.width = PW;
        request.height = PH;
        request.sub = 1;
        request.maxIterations = PMIT;
        request.iterations = dispatched.data();
        request.expression = &cubicProgram;
        request.expressionFixed = &fixed;
        request.expressionPixel = FormulaParameter::C;
        request.expressionBailout = 4.0;
        request.expressionColoring = formula::ExpressionColoring::Smooth;
        _putenv_s("MANDEL_CUBIC_RESIDUAL", "1");
        _putenv_s("MANDEL_CUBIC_RESIDUAL_SCALE", "1e8");
        backend->resetCancellation();
        bool backendOkay = backend->compute(request);
        for (size_t i = 0; i < expected.size(); ++i) {
            if (std::memcmp(&expected[i], &dispatched[i], sizeof(float)) != 0) ++backendDeepMismatches;
        }
        if (!expectedOkay || !used || !backendOkay || backendDeepMismatches) ++failures;

        mpf_set_d(scale, 1e7);
        expectedOkay = expectedRenderer.ComputeExpression(centerRe, centerIm, scale, cubicProgram, fixed, FormulaParameter::C, PMIT, 4.0, formula::ExpressionColoring::Smooth);
        std::fill(dispatched.begin(), dispatched.end(), EMPTYPIXEL);
        backend->resetCancellation();
        backendOkay = backend->compute(request);
        for (size_t i = 0; i < expected.size(); ++i) {
            if (std::memcmp(&expected[i], &dispatched[i], sizeof(float)) != 0) ++backendShallowMismatches;
        }
        if (!expectedOkay || !backendOkay || backendShallowMismatches) ++failures;

        mpf_set_d(scale, 1e5);
        _putenv_s("MANDEL_CUBIC_RESIDUAL_SCALE", "1e5");
        expectedOkay = expectedRenderer.ComputeExpression(centerRe, centerIm, scale, cubicProgram, fixed, FormulaParameter::C, PMIT, 4.0, formula::ExpressionColoring::Smooth);
        std::fill(dispatched.begin(), dispatched.end(), EMPTYPIXEL);
        backend->resetCancellation();
        backendOkay = backend->compute(request);
        for (size_t i = 0; i < expected.size(); ++i) {
            if (std::memcmp(&expected[i], &dispatched[i], sizeof(float)) != 0) ++backendPartialMismatches;
        }
        if (!expectedOkay || !backendOkay || backendPartialMismatches) ++failures;
        _putenv_s("MANDEL_CUBIC_RESIDUAL_SCALE", "1e8");

        mpf_set_d(scale, 1e8);
        expectedOkay = expectedRenderer.ComputeExpression(centerRe, centerIm, scale, cubicProgram, fixed, FormulaParameter::C, PMIT, 4.0, formula::ExpressionColoring::Smooth);
        std::fill(dispatched.begin(), dispatched.end(), EMPTYPIXEL);
        _putenv_s("MANDEL_EXPR_RESIDUAL_POWER", "0");
        backend->resetCancellation();
        backendOkay = backend->compute(request);
        for (size_t i = 0; i < expected.size(); ++i) {
            if (std::memcmp(&expected[i], &dispatched[i], sizeof(float)) != 0) ++backendDisabledMismatches;
        }
        if (!expectedOkay || !backendOkay || backendDisabledMismatches) ++failures;
        _putenv_s("MANDEL_EXPR_RESIDUAL_POWER", "");

        for (formula::ExpressionColoring coloring : {formula::ExpressionColoring::Feather, formula::ExpressionColoring::OrbitTrap}) {
            request.expressionColoring = coloring;
            expectedOkay = expectedRenderer.ComputeExpression(centerRe, centerIm, scale, cubicProgram, fixed, FormulaParameter::C, PMIT, 4.0, coloring);
            std::fill(dispatched.begin(), dispatched.end(), EMPTYPIXEL);
            backend->resetCancellation();
            backendOkay = backend->compute(request);
            for (size_t i = 0; i < expected.size(); ++i) {
                if (std::memcmp(&expected[i], &dispatched[i], sizeof(float)) != 0) ++backendColoredMismatches;
            }
            if (!expectedOkay || !backendOkay || backendColoredMismatches) ++failures;
        }
        request.expressionColoring = formula::ExpressionColoring::Smooth;

        fixed.c = {0.0, 0.0};
        request.expressionPixel = FormulaParameter::InitialZ;
        mpf_set_ui(centerRe, 1);
        mpf_set_d(scale, 1e10);
        expectedOkay = expectedRenderer.ComputeExpression(centerRe, centerIm, scale, cubicProgram, fixed, FormulaParameter::InitialZ, PMIT, 4.0, formula::ExpressionColoring::Smooth);
        std::fill(dispatched.begin(), dispatched.end(), EMPTYPIXEL);
        backend->resetCancellation();
        backendOkay = backend->compute(request);
        for (size_t i = 0; i < expected.size(); ++i) {
            if (std::memcmp(&expected[i], &dispatched[i], sizeof(float)) != 0) ++backendZ0Mismatches;
        }
        if (!expectedOkay || !backendOkay || backendZ0Mismatches) ++failures;
        _putenv_s("MANDEL_CUBIC_RESIDUAL", "");
        _putenv_s("MANDEL_CUBIC_RESIDUAL_SCALE", "");
        mpf_clears(centerRe, centerIm, scale, (mpf_ptr)0);
    }

    printf("  aggregate bits=%d class=%d max=%.6g\n", totalBitMismatches, totalClassMismatches, maximumDifference);
    printf("  cubic %dx%d mxit=%d scalar/AVX2 %.3f/%.3f s speedup %.2fx\n", BW, BH, BMIT, scalarMedian, simdMedian, speedup);
    printf("  cubic residual generic/specialized %.3f/%.3f s speedup %.2fx"
           " bits=%d\n",
           genericResidualMedian, cubicResidualMedian, residualSpeedup, residualBitMismatches);
    printf("  cubic SA iterations=%d off/on %.3f/%.3f s speedup %.2fx"
           " class=%d floor=%d\n",
           seriesIterations, cubicNoSeriesTime, cubicWithSeriesTime, seriesSpeedup, seriesClassMismatches, seriesFloorMismatches);
    printf("  cubic SA shallow iteration=%d bits=%d\n", shallowSeriesIterations, shallowMismatches);
    printf("  cubic production threshold candidate=%.0e\n", firstProfitableScale);
    printf("  cubic backend deep/shallow/partial/disabled/z0/colored"
           " bits=%d/%d/%d/%d/%d/%d\n",
           backendDeepMismatches, backendShallowMismatches, backendPartialMismatches, backendDisabledMismatches, backendZ0Mismatches, backendColoredMismatches);
    printf("  => %s\n\n", failures == 0 ? "PASS" : "CHECK (Multibrot correctness/performance)");
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
        formula::ExpressionColoring coloring = formula::ExpressionColoring::Raw;
    };
    const ResidualCase cases[] = {
        {"quadratic-c", "z*z+c", FormulaParameter::C, {-0.75, 0.0}, 1e4, {}, {}, {}, 1000, 4.0, true}, {"cubic-c", "z*z*z+c", FormulaParameter::C, {0.3849001794597505, 0.0}, 1e5, {}, {}, {}, 1500, 4.0, true}, {"sine-c", "sin(z)+c", FormulaParameter::C, {}, 4.0, {}, {}, {}, 200, 8.0, true}, {"burning-c", "sqr(complex(abs(re(z)),abs(im(z))))+c", FormulaParameter::C, {}, 4.0, {}, {}, {}, 200, 4.0, true}, {"parameter-c", "z*z+c+p0*z", FormulaParameter::C, {}, 4.0, {}, {}, {0.15, -0.05}, 300, 4.0, true}, {"quadratic-z0", "z*z+c", FormulaParameter::InitialZ, {}, 1.0, {}, {}, {}, 500, 4.0, true}, {"escaping-reference-fallback", "z*z+c", FormulaParameter::C, {0.5, 0.5}, 10.0, {}, {}, {}, 300, 4.0, false}, {"cubic-escaping-smooth", "z*z*z+c", FormulaParameter::C, {1.0, 1.0}, 10.0, {}, {}, {}, 300, 4.0, false, formula::ExpressionColoring::Smooth}, {"cubic-small-bailout-distance", "z*z*z+c", FormulaParameter::C, {}, 4.0, {}, {}, {}, 100, 0.5, true, formula::ExpressionColoring::Distance}};

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
        bool ok = renderer.ComputeExpression(centerRe, centerIm, scale, program, fixed, test.pixel, test.mxit, test.bailout, test.coloring);
        double directTime = since(begin);
        direct = output;
        bool usedResidual = false;
        begin = Clock::now();
        ok = ok && renderer.ComputeExpressionResidual(centerRe, centerIm, scale, program, fixed, test.pixel, test.mxit, test.bailout, &usedResidual, test.coloring);
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
        bool passed = ok && empty == 0 && usedResidual == test.expectResidual && classMismatch == 0 && iterationMismatch == 0 && maxDifference == 0;
        if (!passed) ++failures;
        printf("=== residual %s\n", test.name);
        printf("  direct/residual %.3f/%.3f s used=%d expected=%d\n", directTime, residualTime, usedResidual ? 1 : 0, test.expectResidual ? 1 : 0);
        printf("  empty=%d class mismatch=%d iteration mismatch=%d max=%d\n", empty, classMismatch, iterationMismatch, maxDifference);
        printf("  => %s\n\n", passed ? "PASS" : "CHECK (residual mismatch)");
    }
    return failures == 0 ? 0 : 1;
}

static int runCustomDeepZoomCase() {
    using formula::Complex;
    using formula::CustomDeepZoomOutputAdapter;
    using formula::ExpressionColoring;
    using formula::ExpressionContext;
    using formula::ExpressionError;
    using formula::ExpressionProgram;

    struct PreparedFormula {
        ExpressionProgram source;
        ExpressionProgram runtime;
        ExpressionContext fixed;
    };
    auto prepare = [](const char* source, Complex z0 = {}) {
        PreparedFormula prepared;
        prepared.fixed.z0 = z0;
        ExpressionError error;
        if (!prepared.source.compile(source, &error) || !prepared.source.specialize(prepared.fixed, FormulaParameter::C, prepared.runtime, &error)) { fprintf(stderr, "custom-deep compile failed: %s\n", error.message.c_str()); }
        return prepared;
    };

    int failures = 0;
    float radiusBuffer[4]{};
    Mandel radiusProbe(2, 2, 16, 1, radiusBuffer);
    const double productionRadius = radiusProbe.escapeRadius();
    const double bailout = 4.0;
    PreparedFormula quadratic = prepare("z*z+c");
    mpf_t capabilityCenterRe, capabilityCenterIm;
    mpf_init_set_ui(capabilityCenterRe, 0);
    mpf_init_set_ui(capabilityCenterIm, 0);

    auto plan = [&](const PreparedFormula& prepared, FormulaParameter pixel, double bailout, ExpressionColoring coloring, mpf_srcptr scale = nullptr, mpf_srcptr centerRe = nullptr, mpf_srcptr centerIm = nullptr) {
        if (!centerRe) centerRe = capabilityCenterRe;
        if (!centerIm) centerIm = capabilityCenterIm;
        return formula::makeCustomDeepZoomPlan(prepared.source, prepared.runtime, prepared.fixed, pixel, bailout, coloring, scale, centerRe, centerIm, 48, 32);
    };
    auto expect = [&](bool condition, const char* name) {
        if (!condition) {
            ++failures;
            printf("  FAIL: %s\n", name);
        }
    };

    mpf_t below, at, above;
    mpf_init_set_d(below, 9.99999999999e11);
    mpf_init_set_d(at, formula::CUSTOM_DIRECT_ZOOM_LIMIT);
    mpf_init_set_d(above, 1.000000000001e12);
    auto acceptedBelow = plan(quadratic, FormulaParameter::C, bailout, ExpressionColoring::Smooth, below);
    auto acceptedAt = plan(quadratic, FormulaParameter::C, bailout, ExpressionColoring::Smooth, at);
    auto acceptedAbove = plan(quadratic, FormulaParameter::C, bailout, ExpressionColoring::Smooth, above);
    expect(acceptedBelow.canZoomBeyondDirectLimit() && !acceptedBelow.usesQuadraticPerturbation(), "compatible formula stays direct below crossover");
    expect(acceptedAt.canZoomBeyondDirectLimit() && !acceptedAt.usesQuadraticPerturbation(), "compatible formula stays direct through crossover");
    expect(acceptedAbove.usesQuadraticPerturbation(), "default bailout-4 formula dispatches above crossover");
    PreparedFormula canonicalSquare = prepare("sqr(z)+c");
    expect(plan(canonicalSquare, FormulaParameter::C, bailout, ExpressionColoring::Smooth, above).usesQuadraticPerturbation(), "canonical square opcode accepted");
    expect(plan(canonicalSquare, FormulaParameter::C, bailout, ExpressionColoring::Feather, above).usesQuadraticPerturbation() && plan(canonicalSquare, FormulaParameter::C, bailout, ExpressionColoring::OrbitTrap, above).usesQuadraticPerturbation(), "canonical square colored adapters accepted");
    expect(plan(quadratic, FormulaParameter::C, bailout, ExpressionColoring::Distance, above).usesQuadraticPerturbation(), "distance dispatch accepted");
    expect(!plan(quadratic, FormulaParameter::C, bailout, ExpressionColoring::Raw, above).canZoomBeyondDirectLimit(), "quadratic specialized raw adapter remains unavailable");
    expect(formula::makeExpressionProductionPlan(quadratic.source, quadratic.runtime, quadratic.fixed, FormulaParameter::C, bailout, ExpressionColoring::Raw, above, capabilityCenterRe, capabilityCenterIm, 48, 32).usesGenericCertifiedDeep(), "raw quadratic receives generic certified deep fallback");
    auto featherPlan = plan(quadratic, FormulaParameter::C, bailout, ExpressionColoring::Feather, above);
    auto trapPlan = plan(quadratic, FormulaParameter::C, bailout, ExpressionColoring::OrbitTrap, above);
    expect(featherPlan.usesQuadraticPerturbation() && featherPlan.outputAdapter == CustomDeepZoomOutputAdapter::FeatherExpression, "feather dispatch accepted with z1 adapter");
    expect(trapPlan.usesQuadraticPerturbation() && trapPlan.outputAdapter == CustomDeepZoomOutputAdapter::OrbitTrapExpression, "orbit trap dispatch accepted with z1 adapter");
    expect(!plan(quadratic, FormulaParameter::InitialZ, bailout, ExpressionColoring::Smooth, above).canZoomBeyondDirectLimit(), "z0-plane rejected");
    expect(plan(quadratic, FormulaParameter::C, 2.0, ExpressionColoring::Smooth, above).usesQuadraticPerturbation(), "minimum safe bailout accepted independently of production radius");
    expect(!plan(quadratic, FormulaParameter::C, std::nextafter(2.0, 0.0), ExpressionColoring::Smooth, above).canZoomBeyondDirectLimit(), "tiny bailout rejected");
    expect(plan(quadratic, FormulaParameter::C, formula::customDeepMaxEscapeRadius(), ExpressionColoring::Smooth, above).usesQuadraticPerturbation(), "maximum overflow-safe bailout accepted");
    expect(!plan(quadratic, FormulaParameter::C, std::nextafter(formula::customDeepMaxEscapeRadius(), std::numeric_limits<double>::infinity()), ExpressionColoring::Smooth, above).canZoomBeyondDirectLimit(), "overflowing bailout rejected");
    expect(!plan(quadratic, FormulaParameter::C, std::numeric_limits<double>::infinity(), ExpressionColoring::Smooth, above).canZoomBeyondDirectLimit(), "infinite bailout rejected");
    expect(!plan(quadratic, FormulaParameter::C, std::numeric_limits<double>::quiet_NaN(), ExpressionColoring::Smooth, above).canZoomBeyondDirectLimit(), "NaN bailout rejected");
    expect(!plan(quadratic, FormulaParameter::C, std::numeric_limits<double>::max(), ExpressionColoring::Smooth, above).canZoomBeyondDirectLimit(), "radius with non-finite square rejected");

    PreparedFormula nonzeroZ0 = prepare("z*z+c", {1e-300, 0.0});
    PreparedFormula negativeZeroRe = prepare("z*z+c", {-0.0, 0.0});
    PreparedFormula negativeZeroIm = prepare("z*z+c", {0.0, -0.0});
    expect(!plan(nonzeroZ0, FormulaParameter::C, bailout, ExpressionColoring::Smooth, above).canZoomBeyondDirectLimit(), "nonzero z0 rejected");
    expect(!plan(negativeZeroRe, FormulaParameter::C, bailout, ExpressionColoring::Smooth, above).canZoomBeyondDirectLimit(), "negative-zero real z0 rejected");
    expect(!plan(negativeZeroIm, FormulaParameter::C, bailout, ExpressionColoring::Smooth, above).canZoomBeyondDirectLimit(), "negative-zero imaginary z0 rejected");

    PreparedFormula reordered = prepare("c+z*z");
    PreparedFormula parameterized = prepare("z*z+c+p0");
    parameterized.fixed.parameters[0] = {};
    ExpressionError parameterError;
    parameterized.source.specialize(parameterized.fixed, FormulaParameter::C, parameterized.runtime, &parameterError);
    expect(!plan(reordered, FormulaParameter::C, bailout, ExpressionColoring::Smooth, above).canZoomBeyondDirectLimit(), "algebraically reordered recurrence rejected");
    expect(!plan(parameterized, FormulaParameter::C, bailout, ExpressionColoring::Smooth, above).canZoomBeyondDirectLimit(), "folded custom parameter recurrence rejected");
    PreparedFormula unusedParameters = quadratic;
    unusedParameters.fixed.parameters[0] = {9.0, -4.0};
    expect(plan(unusedParameters, FormulaParameter::C, bailout, ExpressionColoring::Smooth, above).canZoomBeyondDirectLimit(), "unused parameters do not affect semantics");
    mpf_t outsideCenter;
    mpf_init_set_d(outsideCenter, bailout * 2.0);
    expect(!plan(quadratic, FormulaParameter::C, bailout, ExpressionColoring::Smooth, above, outsideCenter, capabilityCenterIm).canZoomBeyondDirectLimit(), "z1 outside Custom bailout rejected");
    mpf_clear(outsideCenter);
    mpf_clears(below, at, above, (mpf_ptr)0);

    {
        MandelNavigator nav(16, 12, 1, 1000, 1.0, 1000.0);
        std::array<Complex, 8> parameters{};
        ExpressionError error;
        expect(nav.SetLocation("-0.75", "0", "1000000000000000"), "Mandel location setup");
        expect(nav.SetExpressionFormula("z*z+c", FormulaParameter::C, {}, {}, parameters, bailout, &error), "compatible formula apply");
        mpf_t re, im, scale;
        mpf_inits(re, im, scale, (mpf_ptr)0);
        nav.GetView(re, im, scale);
        expect(mpf_cmp_d(scale, 1e15) == 0, "compatible apply preserves deep view");
        expect(nav.GetCustomDeepZoomPlan().usesQuadraticPerturbation(), "navigator reports deep dispatch");
        expect(nav.GetExpressionAccelerationText().find("deep quadratic perturbation") != std::string::npos, "deep status text");
        expect(nav.SetExpressionFormula("z*z+c+0", FormulaParameter::C, {}, {}, parameters, bailout, &error), "generic raw formula apply succeeds");
        nav.GetView(re, im, scale);
        expect(mpf_cmp_d(scale, 1e15) == 0, "generic raw apply preserves deep view");
        expect(!nav.GetCustomDeepZoomPlan().usesQuadraticPerturbation() && nav.GetExpressionAccelerationText().find("generic deep") != std::string::npos, "navigator reports generic deep dispatch");
        nav.ZoomIn(8, 6);
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        nav.Update();
        nav.UpdateCoords();
        nav.GetView(re, im, scale);
        expect(mpf_cmp_d(scale, 1e15) > 0, "generic raw interactive zoom crosses cap");
        nav.JumpReset();

        nav.SetCMethod(ColoringMethod::STRIPE_AVERAGE);
        nav.GetView(re, im, scale);
        expect(mpf_cmp_d(scale, formula::CUSTOM_DIRECT_ZOOM_LIMIT) == 0, "generic Feather switch clamps transactionally");
        nav.ZoomIn(8, 6);
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        nav.Update();
        nav.UpdateCoords();
        nav.GetView(re, im, scale);
        expect(mpf_cmp_d(scale, formula::CUSTOM_DIRECT_ZOOM_LIMIT) == 0, "generic Feather interactive zoom remains capped");
        nav.JumpReset();
        nav.SetCMethod(0);
        nav.GetView(re, im, scale);
        expect(mpf_cmp_d(scale, formula::CUSTOM_DIRECT_ZOOM_LIMIT) == 0, "returning to raw does not restore lost depth");
        nav.ZoomIn(8, 6);
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        nav.Update();
        nav.UpdateCoords();
        nav.GetView(re, im, scale);
        expect(mpf_cmp_d(scale, formula::CUSTOM_DIRECT_ZOOM_LIMIT) > 0, "returning to raw permits further zoom");
        nav.JumpReset();

        expect(nav.SetExpressionFormula("z*z+c", FormulaParameter::C, {}, {}, parameters, bailout, &error), "compatible formula reapply");
        nav.ZoomIn(8, 6);
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        nav.Update();
        nav.UpdateCoords();
        nav.GetView(re, im, scale);
        expect(mpf_cmp_d(scale, formula::CUSTOM_DIRECT_ZOOM_LIMIT) > 0, "compatible interactive zoom crosses cap");
        nav.JumpReset();
        expect(nav.SetLocation("8", "0", "1000000000000"), "outside-bailout cap-boundary setup");
        nav.ZoomIn(0, 0);
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        nav.Update();
        nav.UpdateCoords();
        nav.GetView(re, im, scale);
        expect(mpf_cmp_d(scale, formula::CUSTOM_DIRECT_ZOOM_LIMIT) == 0 && mpf_cmp_d(re, 8.0) == 0, "outside-bailout interactive crossing stays capped");
        nav.JumpReset();
        expect(nav.SetLocation("-0.75", "0", "1000000000000000000000000000000"), "compatible paste above cap");
        nav.GetView(re, im, scale);
        mpf_t expectedScale;
        mpf_init2(expectedScale, nav.GetViewPrecision());
        mpf_set_str(expectedScale, "1000000000000000000000000000000", 10);
        expect(mpf_cmp(scale, expectedScale) == 0 && nav.GetViewPrecision() > 120, "compatible SetLocation preserves pasted scale");
        nav.SetCMethod(ColoringMethod::STRIPE_AVERAGE);
        nav.GetView(re, im, scale);
        expect(mpf_cmp(scale, expectedScale) == 0 && nav.GetCustomDeepZoomPlan().usesQuadraticPerturbation(), "Feather switch preserves deep view");
        expect(nav.GetExpressionAccelerationText().find("deep quadratic perturbation") != std::string::npos, "Feather deep status text");
        nav.SetCMethod(ColoringMethod::ORBIT_TRAP);
        nav.GetView(re, im, scale);
        expect(mpf_cmp(scale, expectedScale) == 0 && nav.GetCustomDeepZoomPlan().usesQuadraticPerturbation(), "OrbitTrap switch preserves deep view");

        expect(nav.SetLocation("-0.75", "0", "1000000000000001"), "compatible coloring paste accepted");
        nav.GetView(re, im, scale);
        expect(mpf_cmp_d(scale, 1000000000000001.0) == 0, "OrbitTrap SetLocation preserves deep scale");
        mpf_clear(expectedScale);
        mpf_clears(re, im, scale, (mpf_ptr)0);
    }
    {
        MandelNavigator nav(16, 12, 1, 200, 1.0, 1.0);
        std::array<Complex, 8> parameters{};
        ExpressionError error;
        const std::string deepScale = "1e500";
        const std::string deepCenter = "-0.75" + std::string(497, '0') + "1";
        expect(nav.SetLocation(deepCenter, "0", deepScale), "generic e500 Mandel setup");
        expect(nav.SetExpressionFormula("z*z+c+0", FormulaParameter::C, {}, {}, parameters, bailout, &error), "generic c-plane raw apply");
        mpf_t re, im, scale;
        mpf_inits(re, im, scale, (mpf_ptr)0);
        nav.GetView(re, im, scale);
        expect(mpf_cmp_d(scale, 1e100) > 0 && nav.GetExpressionAccelerationText().find("generic deep") != std::string::npos, "generic c-plane e500 accepted");
        const std::string exactText = nav.GetLocationText(true);
        auto lineValue = [&](const char* key) {
            const size_t begin = exactText.find(key);
            if (begin == std::string::npos) return std::string();
            const size_t valueBegin = begin + strlen(key);
            const size_t end = exactText.find_first_of("\r\n", valueBegin);
            return exactText.substr(valueBegin, end == std::string::npos ? std::string::npos : end - valueBegin);
        };
        const std::string copiedX = lineValue("x: ");
        const std::string copiedY = lineValue("y: ");
        const std::string copiedScale = lineValue("zoom: ");
        MandelNavigator roundTrip(16, 12, 1, 200, 1.0, 1.0);
        expect(roundTrip.SetLocation(copiedX, copiedY, copiedScale, nav.GetViewPrecision()), "exact e500 location text reparses");
        mpf_t copiedRe, copiedIm, copiedZoom;
        mpf_inits(copiedRe, copiedIm, copiedZoom, (mpf_ptr)0);
        roundTrip.GetView(copiedRe, copiedIm, copiedZoom);
        const std::string recopied = roundTrip.GetLocationText(true);
        expect(copiedX.rfind(deepCenter, 0) == 0 && recopied.find("x: " + deepCenter) != std::string::npos && copiedScale.find("e+500") != std::string::npos && recopied.find("e+500") != std::string::npos, "exact e500 coordinate digits roundtrip");

        expect(nav.SetExpressionFormula("z*z+c+0", FormulaParameter::InitialZ, {}, {-0.2, 0.1}, parameters, bailout, &error), "generic z0-plane raw apply");
        nav.GetView(re, im, scale);
        expect(mpf_cmp_d(scale, 1e100) > 0, "generic z0-plane preserves e500 view");
        nav.SetCMethod(ColoringMethod::ORBIT_TRAP);
        nav.GetView(re, im, scale);
        expect(mpf_cmp_d(scale, formula::CUSTOM_DIRECT_ZOOM_LIMIT) == 0, "generic OrbitTrap clamps at direct limit");
        nav.SetCMethod(0);
        expect(nav.SetLocation("8", "0", deepScale), "generic raw accepts arbitrary e500 center");
        nav.GetView(re, im, scale);
        expect(mpf_cmp_d(scale, 1e100) > 0 && mpf_cmp_ui(re, 8) == 0, "generic raw tentative center is not quadratic-gated");
        nav.ZoomIn(0, 0);
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        nav.Update();
        nav.UpdateCoords();
        nav.GetView(re, im, scale);
        expect(mpf_cmp_d(scale, 1e100) > 0, "generic raw wheel tentative view remains e500-capable");
        nav.JumpReset();
        expect(nav.SetExpressionFormula("z*z*z+c", FormulaParameter::C, {}, {}, parameters, bailout, &error), "unsupported smooth generic formula applies");
        nav.GetView(re, im, scale);
        expect(mpf_cmp_d(scale, formula::CUSTOM_DIRECT_ZOOM_LIMIT) == 0, "unsupported deep formula switch clamps transactionally");
        const std::string unsupportedScale = "1" + std::string(20000, '0');
        expect(nav.SetLocation("0", "0", unsupportedScale), "over-limit generic location parses");
        nav.GetView(re, im, scale);
        expect(mpf_cmp_d(scale, formula::CUSTOM_DIRECT_ZOOM_LIMIT) == 0, "generic precision cap prevents an unrenderable view");
        mpf_clears(copiedRe, copiedIm, copiedZoom, re, im, scale, (mpf_ptr)0);
    }
    {
        MandelNavigator julia(16, 12, 1, 1000, 1.0, 1.0);
        julia.SetJuliaMode(true);
        expect(julia.SetLocation("0", "0", "1000000000000001"), "Julia SetLocation accepted");
        mpf_t re, im, scale;
        mpf_inits(re, im, scale, (mpf_ptr)0);
        julia.GetView(re, im, scale);
        expect(mpf_cmp_d(scale, formula::CUSTOM_DIRECT_ZOOM_LIMIT) == 0, "Julia remains capped");
        expect(julia.SetJuliaC("0.3", "-0.2"), "legacy Julia constants without a precision hint remain valid");
        const std::string juliaConstant = "0." + std::string(80, '0') + "1";
        expect(julia.SetJuliaC(juliaConstant, "-0", 256), "Julia precision hint accepted");
        mpf_t cRe, cIm;
        mpf_inits(cRe, cIm, (mpf_ptr)0);
        julia.GetJuliaC(cRe, cIm);
        expect(mpf_get_prec(cRe) == 256 && mpf_get_prec(cIm) == 256, "Julia precision hint applied exactly");
        const std::string juliaCopy = julia.GetLocationText(true);
        const size_t hintAt = juliaCopy.find("precision: ");
        expect(hintAt != std::string::npos && juliaCopy.find("precision: 256", hintAt) != std::string::npos, "Julia copy preserves bounded precision hint");
        mpf_clears(cRe, cIm, (mpf_ptr)0);
        mpf_clears(re, im, scale, (mpf_ptr)0);
    }
    {
        MandelNavigator nav(9, 7, 1, 1000000, 1.0, 1.0);
        std::array<Complex, 8> parameters{};
        ExpressionError error;
        expect(nav.SetLocation("0", "0", "1e500") && nav.SetExpressionFormula("z*z+c+0", FormulaParameter::C, {}, {}, parameters, bailout, &error), "generic lifetime setup");
        nav.StartCompute();
        std::this_thread::sleep_for(std::chrono::milliseconds(2));
        const auto swapStart = std::chrono::steady_clock::now();
        const bool swapped = nav.SetExpressionFormula("sin(z)+c", FormulaParameter::C, {}, {}, parameters, bailout, &error);
        const double swapSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - swapStart).count();
        expect(swapped && swapSeconds < 3.0 && !nav.IsComputing(), "formula reapply interrupts before pointer swap");
        nav.SetMxit(32);
        nav.StartCompute();
        const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(5);
        while (nav.IsComputing() && std::chrono::steady_clock::now() < deadline) std::this_thread::sleep_for(std::chrono::milliseconds(2));
        const GenericDeepInfo info = nav.GetLastGenericDeepInfo();
        expect(!nav.IsComputing() && nav.LastComputeUsedGenericDeepPath() && info.success, "reapplied formula owns stable async request state");
    }
    {
        MandelNavigator fromJulia(16, 12, 1, 1000, 1.0, 1.0);
        std::array<Complex, 8> parameters{};
        ExpressionError error;
        expect(fromJulia.SetLocation("-0.75", "0", "1000000000000000"), "Julia transition Mandel setup");
        fromJulia.SetJuliaMode(true);
        expect(fromJulia.SetExpressionFormula("z*z+c", FormulaParameter::C, {}, {}, parameters, bailout, &error), "compatible formula apply from Julia");
        mpf_t re, im, scale;
        mpf_inits(re, im, scale, (mpf_ptr)0);
        fromJulia.GetView(re, im, scale);
        expect(mpf_cmp_d(scale, 1e15) == 0 && fromJulia.GetCustomDeepZoomPlan().usesQuadraticPerturbation(), "Julia-to-Custom preserves compatible saved deep view");
        mpf_clears(re, im, scale, (mpf_ptr)0);
    }

    struct OverlapResult {
        int classMismatch = 0;
        int floorMismatch = 0;
        int byteMismatch = 0;
        double maxDelta = 0.0;
        double maxRelative = 0.0;
    };
    auto adapterFor = [](ExpressionColoring coloring) { return coloring == ExpressionColoring::Smooth ? CustomDeepZoomOutputAdapter::SmoothExpression : coloring == ExpressionColoring::Distance ? CustomDeepZoomOutputAdapter::DistanceExpression
                                                                                                                                                    : coloring == ExpressionColoring::Feather    ? CustomDeepZoomOutputAdapter::FeatherExpression
                                                                                                                                                    : coloring == ExpressionColoring::OrbitTrap  ? CustomDeepZoomOutputAdapter::OrbitTrapExpression
                                                                                                                                                                                                 : CustomDeepZoomOutputAdapter::None; };
    auto methodFor = [](ExpressionColoring coloring) { return coloring == ExpressionColoring::Distance ? ColoringMethod::EXTERIOR_DIST_EST : coloring == ExpressionColoring::Feather ? ColoringMethod::STRIPE_AVERAGE
                                                                                                                                         : coloring == ExpressionColoring::OrbitTrap ? ColoringMethod::ORBIT_TRAP
                                                                                                                                                                                     : 0; };
    auto compareOverlap = [&](ExpressionColoring coloring, const char* scaleText = "1000000000000") {
        constexpr int W = 48, H = 32, MXIT = 4000;
        std::vector<float> direct((size_t)W * H, EMPTYPIXEL);
        std::vector<float> deep((size_t)W * H, EMPTYPIXEL);
        Mandel directEngine(W, H, MXIT, 1, direct.data());
        Mandel deepEngine(W, H, MXIT, 1, deep.data());
        mpf_t re, im, scale;
        mpf_init2(re, 256);
        mpf_init2(im, 256);
        mpf_init2(scale, 256);
        mpf_set_str(re, "-0.75", 10);
        mpf_set_str(im, "0.1", 10);
        mpf_set_str(scale, scaleText, 10);
        int method = methodFor(coloring);
        bool okay = directEngine.ComputeExpression(re, im, scale, quadratic.runtime, quadratic.fixed, FormulaParameter::C, MXIT, bailout, coloring);
        CustomDeepZoomOutputAdapter adapter = adapterFor(coloring);
        {
            Mandel::ScopedCustomCompute customCompute(deepEngine, bailout, adapter);
            okay = okay && customCompute.active();
            if (customCompute.active()) deepEngine.Compute(re, im, scale, MXIT, method);
        }
        mpf_clears(re, im, scale, (mpf_ptr)0);

        OverlapResult result;
        for (size_t i = 0; i < direct.size(); ++i) {
            if (std::memcmp(&direct[i], &deep[i], sizeof(float)) != 0) ++result.byteMismatch;
            if (coloring != ExpressionColoring::OrbitTrap && isInterior(direct[i]) != isInterior(deep[i])) {
                ++result.classMismatch;
                continue;
            }
            if (coloring == ExpressionColoring::Smooth && !isInterior(direct[i]) && std::floor(direct[i]) != std::floor(deep[i])) { ++result.floorMismatch; }
            double delta = std::fabs((double)direct[i] - deep[i]);
            result.maxDelta = std::max(result.maxDelta, delta);
            double denominator = std::max(1e-30, std::fabs((double)direct[i]));
            result.maxRelative = std::max(result.maxRelative, delta / denominator);
        }
        if (!okay) ++failures;
        return result;
    };

    OverlapResult smooth = compareOverlap(ExpressionColoring::Smooth);
    OverlapResult distance = compareOverlap(ExpressionColoring::Distance);
    OverlapResult feather = compareOverlap(ExpressionColoring::Feather);
    OverlapResult trap = compareOverlap(ExpressionColoring::OrbitTrap);
    OverlapResult featherBelow = compareOverlap(ExpressionColoring::Feather, "999999999999");
    OverlapResult featherAbove = compareOverlap(ExpressionColoring::Feather, "1000000000001");
    OverlapResult featherShallow = compareOverlap(ExpressionColoring::Feather, "100");
    OverlapResult trapBelow = compareOverlap(ExpressionColoring::OrbitTrap, "999999999999");
    OverlapResult trapAbove = compareOverlap(ExpressionColoring::OrbitTrap, "1000000000001");
    OverlapResult trapShallow = compareOverlap(ExpressionColoring::OrbitTrap, "100");
    expect(smooth.classMismatch == 0 && smooth.floorMismatch == 0 && smooth.maxDelta <= 0.05, "smooth overlap class/floor/value");
    expect(distance.classMismatch == 0 && distance.maxDelta <= 0.01, "distance overlap class/value");
    expect(feather.byteMismatch == 0 && featherBelow.byteMismatch == 0 && featherAbove.byteMismatch == 0 && featherShallow.byteMismatch == 0, "feather shallow/crossover is byte-identical");
    expect(trap.byteMismatch == 0 && trapBelow.byteMismatch == 0 && trapAbove.byteMismatch == 0 && trapShallow.byteMismatch == 0, "orbit-trap shallow/crossover is byte-identical");

    {
        constexpr int W = 4, H = 4, MXIT = 16;
        const double upperBailout = formula::customDeepMaxEscapeRadius();
        mpf_t re, im, scale;
        mpf_init2(re, 256);
        mpf_init2(im, 256);
        mpf_init2(scale, 256);
        mpf_set_str(re, "0.65", 10);
        mpf_set_ui(im, 0);
        mpf_set_str(scale, "1000000000001", 10);
        auto compareUpperBailout = [&](ExpressionColoring coloring) {
            int method = coloring == ExpressionColoring::Distance ? ColoringMethod::EXTERIOR_DIST_EST : 0;
            CustomDeepZoomOutputAdapter adapter = coloring == ExpressionColoring::Distance ? CustomDeepZoomOutputAdapter::DistanceExpression : CustomDeepZoomOutputAdapter::SmoothExpression;
            std::vector<float> direct((size_t)W * H, EMPTYPIXEL);
            std::vector<float> highPrecision((size_t)W * H, EMPTYPIXEL);
            Mandel directEngine(W, H, MXIT, 1, direct.data());
            Mandel highPrecisionEngine(W, H, MXIT, 1, highPrecision.data());
            bool okay = directEngine.ComputeExpression(re, im, scale, quadratic.runtime, quadratic.fixed, FormulaParameter::C, MXIT, upperBailout, coloring);
            {
                Mandel::ScopedCustomCompute customCompute(highPrecisionEngine, upperBailout, adapter);
                okay = okay && customCompute.active();
                if (customCompute.active()) { highPrecisionEngine.Compute(re, im, scale, MXIT, method); }
            }
            double maxRelative = 0.0;
            for (size_t i = 0; i < direct.size(); ++i) {
                double expected = direct[i];
                double actual = highPrecision[i];
                if (!(std::isfinite(expected) && std::isfinite(actual) && expected > 0.0 && actual > 0.0)) { return std::numeric_limits<double>::infinity(); }
                maxRelative = std::max(maxRelative, std::fabs(actual - expected) / std::max(1.0, std::fabs(expected)));
            }
            return okay ? maxRelative : std::numeric_limits<double>::infinity();
        };
        double upperSmooth = compareUpperBailout(ExpressionColoring::Smooth);
        double upperDistance = compareUpperBailout(ExpressionColoring::Distance);
        expect(upperSmooth <= 1e-5, "upper bailout Smooth matches direct reference");
        expect(upperDistance <= 1e-5, "upper bailout Distance matches high-precision fallback");
        printf("  upper bailout smooth/distance relative=%.9g/%.9g\n", upperSmooth, upperDistance);
        mpf_clears(re, im, scale, (mpf_ptr)0);
    }

    struct ColoredDeepResult {
        bool okay = false;
        bool usedDeep = false;
        bool restored = false;
        int classMismatch = 0;
        int byteMismatch = 0;
        double maxDelta = 0.0;
        float centerValue = EMPTYPIXEL;
        float centerOracle = EMPTYPIXEL;
    };
    auto compareColoredDeepOracle = [&](ExpressionColoring coloring, const char* centerRe, const char* centerIm, const char* scaleText, int width, int height, int mxit, int step) {
        mp_bitcnt_t precision = std::max<mp_bitcnt_t>(256, strlen(scaleText) * 4 + 96);
        mpf_t re, im, scale;
        mpf_init2(re, precision);
        mpf_init2(im, precision);
        mpf_init2(scale, precision);
        mpf_set_str(re, centerRe, 10);
        mpf_set_str(im, centerIm, 10);
        mpf_set_str(scale, scaleText, 10);
        std::vector<float> deep((size_t)width * height, EMPTYPIXEL);
        std::vector<float> oracle((size_t)width * height, EMPTYPIXEL);
        Mandel engine(width, height, mxit, 1, deep.data());
        engine.setPrecision((int)precision);
        std::unique_ptr<IComputeBackend> backend = createComputeBackend("cpu");
        ComputeRequest request;
        request.mode = ComputeMode::Expression;
        request.cpuEngine = &engine;
        request.centerRe = re;
        request.centerIm = im;
        request.scale = scale;
        request.width = width;
        request.height = height;
        request.sub = 1;
        request.maxIterations = mxit;
        request.coloringMethod = methodFor(coloring);
        request.iterations = deep.data();
        request.expressionSource = &quadratic.source;
        request.expression = &quadratic.runtime;
        request.expressionFixed = &quadratic.fixed;
        request.expressionPixel = FormulaParameter::C;
        request.expressionBailout = bailout;
        request.expressionColoring = coloring;
        backend->resetCancellation();
        ColoredDeepResult result;
        result.okay = backend->compute(request);
        result.usedDeep = backend->lastComputeUsedCustomDeepPath();
        {
            Mandel::ScopedCustomCompute scope(engine, bailout, adapterFor(coloring));
            result.okay = result.okay && scope.active();
            if (scope.active()) engine.ComputeDirect(mxit, oracle.data(), step, methodFor(coloring));
        }
        result.restored = engine.escapeRadius() == productionRadius;
        for (int y = 0; y < height; y += step) {
            for (int x = 0; x < width; x += step) {
                size_t index = (size_t)y * width + x;
                if (coloring != ExpressionColoring::OrbitTrap && isInterior(deep[index]) != isInterior(oracle[index])) {
                    ++result.classMismatch;
                    continue;
                }
                if (std::memcmp(&deep[index], &oracle[index], sizeof(float)) != 0) ++result.byteMismatch;
                result.maxDelta = std::max(result.maxDelta, std::fabs((double)deep[index] - oracle[index]));
            }
        }
        const size_t center = (size_t)(height / 2) * width + width / 2;
        result.centerValue = deep[center];
        result.centerOracle = oracle[center];
        mpf_clears(re, im, scale, (mpf_ptr)0);
        return result;
    };

    ColoredDeepResult featherExterior = compareColoredDeepOracle(ExpressionColoring::Feather, "-0.75", "0.1", "10000000000000", 17, 13, 4000, 1);
    ColoredDeepResult trapExterior = compareColoredDeepOracle(ExpressionColoring::OrbitTrap, "-0.75", "0.1", "10000000000000", 17, 13, 4000, 1);
    ColoredDeepResult trapInterior = compareColoredDeepOracle(ExpressionColoring::OrbitTrap, "0", "0", "10000000000000", 9, 7, 96, 1);
    expect(featherExterior.okay && featherExterior.usedDeep && featherExterior.restored && featherExterior.classMismatch == 0 && featherExterior.byteMismatch == 0, "deep Feather exterior matches GMP oracle");
    expect(trapExterior.okay && trapExterior.usedDeep && trapExterior.restored && trapExterior.byteMismatch == 0, "deep OrbitTrap exterior matches GMP oracle");
    expect(trapInterior.okay && trapInterior.usedDeep && trapInterior.restored && trapInterior.byteMismatch == 0 && trapInterior.centerValue == trapInterior.centerOracle, "deep OrbitTrap interior and reference pixel match GMP");
    if (!(trapInterior.okay && trapInterior.usedDeep && trapInterior.restored && trapInterior.byteMismatch == 0 && trapInterior.centerValue == trapInterior.centerOracle)) {
        printf("  trap interior mismatch bytes/max/center="
               "%d/%.9g/%.9g/%.9g\n",
               trapInterior.byteMismatch, trapInterior.maxDelta, trapInterior.centerValue, trapInterior.centerOracle);
    }

    struct HistoryRegressionResult {
        bool okay = false;
        bool usedDeep = false;
        bool directExact = false;
        bool gmpExact = false;
        bool centerExact = false;
        bool seriesWasBuilt = false;
        bool blaWasRequested = false;
        bool skipsDisabled = false;
        long long gmpFallbackPixels = -1;
        int interiorPixels = 0;
        int gmpByteMismatch = 0;
        double gmpMaxDelta = 0.0;
    };
    auto runHistoryRegression = [&](ExpressionColoring coloring) {
        constexpr int W = 5, H = 5, MXIT = 512;
        mpf_t re, im, scale;
        mpf_init2(re, 256);
        mpf_init2(im, 256);
        mpf_init2(scale, 256);
        mpf_set_str(re, "-0.5", 10);
        mpf_set_str(im, "0.5", 10);
        mpf_set_str(scale, "1000000000001", 10);
        std::vector<float> direct((size_t)W * H, EMPTYPIXEL);
        std::vector<float> deep((size_t)W * H, EMPTYPIXEL);
        std::vector<float> gmp((size_t)W * H, EMPTYPIXEL);
        std::vector<float> classification((size_t)W * H, EMPTYPIXEL);
        Mandel directEngine(W, H, MXIT, 1, direct.data());
        Mandel deepEngine(W, H, MXIT, 1, deep.data());
        directEngine.setPrecision(256);
        deepEngine.setPrecision(256);
        HistoryRegressionResult result;
        result.okay = directEngine.ComputeExpression(re, im, scale, quadratic.runtime, quadratic.fixed, FormulaParameter::C, MXIT, bailout, coloring);

        std::unique_ptr<IComputeBackend> backend = createComputeBackend("cpu");
        ComputeRequest request;
        request.mode = ComputeMode::Expression;
        request.cpuEngine = &deepEngine;
        request.centerRe = re;
        request.centerIm = im;
        request.scale = scale;
        request.width = W;
        request.height = H;
        request.sub = 1;
        request.maxIterations = MXIT;
        request.coloringMethod = methodFor(coloring);
        request.iterations = deep.data();
        request.expressionSource = &quadratic.source;
        request.expression = &quadratic.runtime;
        request.expressionFixed = &quadratic.fixed;
        request.expressionPixel = FormulaParameter::C;
        request.expressionBailout = bailout;
        request.expressionColoring = coloring;
        backend->resetCancellation();
        result.okay = result.okay && backend->compute(request);
        result.usedDeep = backend->lastComputeUsedCustomDeepPath();
        result.seriesWasBuilt = deepEngine.lastCustomHistorySeriesWasBuilt();
        result.blaWasRequested = deepEngine.lastCustomHistoryBlaWasRequested();
        result.skipsDisabled = deepEngine.lastCustomHistorySkipsDisabled();
        result.gmpFallbackPixels = deepEngine.lastDeepGmpFallbackPixels();

        {
            Mandel::ScopedCustomCompute scope(deepEngine, bailout, adapterFor(coloring));
            result.okay = result.okay && scope.active();
            if (scope.active()) { deepEngine.ComputeDirect(MXIT, gmp.data(), 1, methodFor(coloring)); }
        }
        {
            Mandel::ScopedCustomCompute scope(deepEngine, bailout, CustomDeepZoomOutputAdapter::None);
            result.okay = result.okay && scope.active();
            if (scope.active()) deepEngine.ComputeDirect(MXIT, classification.data(), 1, 0);
        }
        result.directExact = std::memcmp(deep.data(), direct.data(), deep.size() * sizeof(float)) == 0;
        result.gmpExact = std::memcmp(deep.data(), gmp.data(), deep.size() * sizeof(float)) == 0;
        for (size_t index = 0; index < deep.size(); ++index) {
            if (std::memcmp(&deep[index], &gmp[index], sizeof(float)) != 0) ++result.gmpByteMismatch;
            result.gmpMaxDelta = std::max(result.gmpMaxDelta, std::fabs((double)deep[index] - gmp[index]));
        }
        const size_t center = (size_t)(H / 2) * W + W / 2;
        result.centerExact = std::memcmp(&deep[center], &direct[center], sizeof(float)) == 0 && std::memcmp(&deep[center], &gmp[center], sizeof(float)) == 0;
        result.interiorPixels = (int)std::count(classification.begin(), classification.end(), -2.0f);
        mpf_clears(re, im, scale, (mpf_ptr)0);
        return result;
    };

    _putenv_s("MANDEL_BLA", "1");
    _putenv_s("MANDEL_FE", "0");
    HistoryRegressionResult featherHistory = runHistoryRegression(ExpressionColoring::Feather);
    HistoryRegressionResult trapHistory = runHistoryRegression(ExpressionColoring::OrbitTrap);
    _putenv_s("MANDEL_FE", "1");
    HistoryRegressionResult trapFloatExpBounded = runHistoryRegression(ExpressionColoring::OrbitTrap);
    _putenv_s("MANDEL_BLA", "");
    _putenv_s("MANDEL_FE", "");

    auto historyOkay = [](const HistoryRegressionResult& result) { return result.okay && result.usedDeep && result.gmpMaxDelta <= 1e-4 && result.centerExact && result.seriesWasBuilt && result.blaWasRequested && result.skipsDisabled && result.interiorPixels > 0; };
    expect(historyOkay(featherHistory) && historyOkay(trapHistory), "Custom Feather/Trap observe every orbit value with SA/BLA disabled");
    expect(historyOkay(trapFloatExpBounded) && trapFloatExpBounded.gmpFallbackPixels == 0, "bounded Custom OrbitTrap avoids cutoff-driven full-frame GMP");
    if (!(historyOkay(featherHistory) && historyOkay(trapHistory) && historyOkay(trapFloatExpBounded) && trapFloatExpBounded.gmpFallbackPixels == 0)) {
        printf("  history feather direct/gmp/SA/BLA/off/int="
               "%d/%d/%d/%d/%d/%d"
               " trap=%d/%d/%d/%d/%d/%d mismatch/max=%d/%.9g"
               " FEtrap=%d/%d/%d/%d/%d/%d mismatch/max=%d/%.9g"
               " fallback=%lld\n",
               featherHistory.directExact, featherHistory.gmpExact, featherHistory.seriesWasBuilt, featherHistory.blaWasRequested, featherHistory.skipsDisabled, featherHistory.interiorPixels, trapHistory.directExact, trapHistory.gmpExact, trapHistory.seriesWasBuilt, trapHistory.blaWasRequested, trapHistory.skipsDisabled, trapHistory.interiorPixels, trapHistory.gmpByteMismatch, trapHistory.gmpMaxDelta, trapFloatExpBounded.directExact, trapFloatExpBounded.gmpExact, trapFloatExpBounded.seriesWasBuilt, trapFloatExpBounded.blaWasRequested, trapFloatExpBounded.skipsDisabled, trapFloatExpBounded.interiorPixels, trapFloatExpBounded.gmpByteMismatch, trapFloatExpBounded.gmpMaxDelta, trapFloatExpBounded.gmpFallbackPixels);
    }

    const std::string floatExpScale = "1" + std::string(1000, '0');
    ColoredDeepResult featherFloatExp = compareColoredDeepOracle(ExpressionColoring::Feather, "-0.5", "0.5", floatExpScale.c_str(), 5, 5, 10000, 4);
    ColoredDeepResult trapFloatExp = compareColoredDeepOracle(ExpressionColoring::OrbitTrap, "-0.5", "0.5", floatExpScale.c_str(), 5, 5, 10000, 4);
    expect(featherFloatExp.okay && featherFloatExp.usedDeep && featherFloatExp.restored && featherFloatExp.classMismatch == 0 && featherFloatExp.maxDelta <= 1e-4, "floatexp Feather matches sampled GMP oracle");
    expect(trapFloatExp.okay && trapFloatExp.usedDeep && trapFloatExp.restored && trapFloatExp.maxDelta <= 1e-4, "floatexp OrbitTrap matches sampled GMP oracle");
    if (!(featherFloatExp.okay && featherFloatExp.usedDeep && featherFloatExp.restored && featherFloatExp.classMismatch == 0 && featherFloatExp.maxDelta <= 1e-4) || !(trapFloatExp.okay && trapFloatExp.usedDeep && trapFloatExp.restored && trapFloatExp.maxDelta <= 1e-4)) {
        printf("  floatexp feather class/bytes/max=%d/%d/%.9g"
               " trap bytes/max=%d/%.9g\n",
               featherFloatExp.classMismatch, featherFloatExp.byteMismatch, featherFloatExp.maxDelta, trapFloatExp.byteMismatch, trapFloatExp.maxDelta);
    }

    auto renderColoredDeep = [&](ExpressionColoring coloring, const char* centerRe, const char* centerIm, const char* scaleText, int width, int height, int mxit, std::vector<float>& output) {
        mp_bitcnt_t precision = std::max<mp_bitcnt_t>(256, strlen(scaleText) * 4 + 96);
        mpf_t re, im, scale;
        mpf_init2(re, precision);
        mpf_init2(im, precision);
        mpf_init2(scale, precision);
        mpf_set_str(re, centerRe, 10);
        mpf_set_str(im, centerIm, 10);
        mpf_set_str(scale, scaleText, 10);
        output.assign((size_t)width * height, EMPTYPIXEL);
        Mandel engine(width, height, mxit, 1, output.data());
        engine.setPrecision((int)precision);
        std::unique_ptr<IComputeBackend> backend = createComputeBackend("cpu");
        ComputeRequest request;
        request.mode = ComputeMode::Expression;
        request.cpuEngine = &engine;
        request.centerRe = re;
        request.centerIm = im;
        request.scale = scale;
        request.width = width;
        request.height = height;
        request.sub = 1;
        request.maxIterations = mxit;
        request.coloringMethod = methodFor(coloring);
        request.iterations = output.data();
        request.expressionSource = &quadratic.source;
        request.expression = &quadratic.runtime;
        request.expressionFixed = &quadratic.fixed;
        request.expressionPixel = FormulaParameter::C;
        request.expressionBailout = bailout;
        request.expressionColoring = coloring;
        backend->resetCancellation();
        bool okay = backend->compute(request) && backend->lastComputeUsedCustomDeepPath() && engine.escapeRadius() == productionRadius;
        mpf_clears(re, im, scale, (mpf_ptr)0);
        return okay;
    };
    std::vector<float> featherBlaOn, featherBlaOff;
    _putenv_s("MANDEL_BLA", "1");
    bool featherBlaOnOkay = renderColoredDeep(ExpressionColoring::Feather, "-0.749139567333446841955467474699747367338762518832278501811", "0.040823298514634751035521346975478853963578400940553676068", "73541770000000000000000000000000", 8, 6, 50000, featherBlaOn);
    _putenv_s("MANDEL_BLA", "0");
    bool featherBlaOffOkay = renderColoredDeep(ExpressionColoring::Feather, "-0.749139567333446841955467474699747367338762518832278501811", "0.040823298514634751035521346975478853963578400940553676068", "73541770000000000000000000000000", 8, 6, 50000, featherBlaOff);
    _putenv_s("MANDEL_BLA", "");
    expect(featherBlaOnOkay && featherBlaOffOkay && featherBlaOn.size() == featherBlaOff.size() && std::memcmp(featherBlaOn.data(), featherBlaOff.data(), featherBlaOn.size() * sizeof(float)) == 0, "deep Feather BLA preserves exact accumulator semantics");

    {
        constexpr int W = 9, H = 7, MXIT = 1000;
        mpf_t re, im, scale;
        mpf_init2(re, 256);
        mpf_init2(im, 256);
        mpf_init2(scale, 256);
        mpf_set_d(re, -0.75);
        mpf_set_d(im, 0.1);
        mpf_set_str(scale, "10000000000000", 10);
        std::vector<float> reused((size_t)W * H, EMPTYPIXEL);
        std::vector<float> fresh((size_t)W * H, EMPTYPIXEL);
        Mandel reusedEngine(W, H, MXIT, 1, reused.data());
        Mandel freshEngine(W, H, MXIT, 1, fresh.data());
        std::unique_ptr<IComputeBackend> backend = createComputeBackend("cpu");
        ComputeRequest request;
        request.mode = ComputeMode::Expression;
        request.cpuEngine = &reusedEngine;
        request.centerRe = re;
        request.centerIm = im;
        request.scale = scale;
        request.width = W;
        request.height = H;
        request.sub = 1;
        request.maxIterations = MXIT;
        request.coloringMethod = ColoringMethod::STRIPE_AVERAGE;
        request.iterations = reused.data();
        request.expressionSource = &quadratic.source;
        request.expression = &quadratic.runtime;
        request.expressionFixed = &quadratic.fixed;
        request.expressionPixel = FormulaParameter::C;
        request.expressionBailout = bailout;
        request.expressionColoring = ExpressionColoring::Feather;
        backend->resetCancellation();
        bool success = backend->compute(request) && backend->lastComputeUsedCustomDeepPath();
        std::vector<float> expectedDeep = reused;
        std::fill(reused.begin(), reused.end(), EMPTYPIXEL);
        request.coloringMethod = ColoringMethod::STRIPE_AVERAGE | ColoringMethod::SUPER_SAMPLING;
        backend->resetCancellation();
        bool orthogonalFlagOkay = backend->compute(request) && backend->lastComputeUsedCustomDeepPath() && std::memcmp(expectedDeep.data(), reused.data(), reused.size() * sizeof(float)) == 0;

        request.coloringMethod = ColoringMethod::ORBIT_TRAP;
        backend->resetCancellation();
        bool rejectedMismatch = !backend->compute(request) && !backend->lastComputeUsedCustomDeepPath();

        bool ordinaryDeepSame = true;
        request.mode = ComputeMode::Mandelbrot;
        request.maxIterations = MXIT;
        for (int method : {ColoringMethod::STRIPE_AVERAGE, ColoringMethod::ORBIT_TRAP}) {
            std::fill(reused.begin(), reused.end(), EMPTYPIXEL);
            std::fill(fresh.begin(), fresh.end(), EMPTYPIXEL);
            freshEngine.Compute(re, im, scale, MXIT, method);
            request.coloringMethod = method;
            backend->resetCancellation();
            ordinaryDeepSame = ordinaryDeepSame && backend->compute(request) && !backend->lastComputeUsedCustomDeepPath() && std::memcmp(fresh.data(), reused.data(), fresh.size() * sizeof(float)) == 0;
        }

        std::fill(reused.begin(), reused.end(), EMPTYPIXEL);
        mpf_set_d(re, -0.5);
        mpf_set_ui(im, 0);
        mpf_set_ui(scale, 1);
        freshEngine.Compute(re, im, scale, 300, ColoringMethod::STRIPE_AVERAGE);
        request.maxIterations = 300;
        request.coloringMethod = ColoringMethod::STRIPE_AVERAGE;
        backend->resetCancellation();
        bool ordinaryOkay = backend->compute(request);
        bool ordinarySame = std::memcmp(fresh.data(), reused.data(), fresh.size() * sizeof(float)) == 0;
        expect(success && orthogonalFlagOkay && rejectedMismatch && ordinaryDeepSame && ordinaryOkay && ordinarySame && reusedEngine.escapeRadius() == productionRadius, "backend validates adapter and restores success/failure state");
        mpf_clears(re, im, scale, (mpf_ptr)0);
    }

    auto adaptExpected = [](std::vector<float>& values, CustomDeepZoomOutputAdapter adapter) {
        const float offset = static_cast<float>(-std::log(std::log(2.0)) / std::log(2.0));
        for (float& value : values) {
            if (value < 0.0f) continue;
            if (adapter == CustomDeepZoomOutputAdapter::SmoothExpression) {
                value += offset;
            } else if (adapter == CustomDeepZoomOutputAdapter::DistanceExpression) {
                value *= 0.5f;
            }
        }
    };
    auto deepParity = [&](const char* name, const char* centerRe, const char* centerIm, const char* scaleText, int mxit) {
        constexpr int W = 20, H = 14;
        mp_bitcnt_t precision = std::max<mp_bitcnt_t>(256, strlen(scaleText) * 4 + 64);
        mpf_t re, im, scale;
        mpf_init2(re, precision);
        mpf_init2(im, precision);
        mpf_init2(scale, precision);
        mpf_set_str(re, centerRe, 10);
        mpf_set_str(im, centerIm, 10);
        mpf_set_str(scale, scaleText, 10);
        std::vector<float> ordinary((size_t)W * H, EMPTYPIXEL);
        std::vector<float> custom((size_t)W * H, EMPTYPIXEL);
        Mandel ordinaryEngine(W, H, mxit, 1, ordinary.data());
        Mandel customEngine(W, H, mxit, 1, custom.data());
        ordinaryEngine.setPrecision((int)precision);
        customEngine.setPrecision((int)precision);
        {
            Mandel::ScopedCustomCompute ordinaryCustom(ordinaryEngine, bailout, CustomDeepZoomOutputAdapter::SmoothExpression);
            if (ordinaryCustom.active()) ordinaryEngine.Compute(re, im, scale, mxit, 0);
        }

        std::unique_ptr<IComputeBackend> backend = createComputeBackend("cpu");
        ComputeRequest request;
        request.mode = ComputeMode::Expression;
        request.cpuEngine = &customEngine;
        request.centerRe = re;
        request.centerIm = im;
        request.scale = scale;
        request.width = W;
        request.height = H;
        request.sub = 1;
        request.maxIterations = mxit;
        request.coloringMethod = 0;
        request.iterations = custom.data();
        request.expressionSource = &quadratic.source;
        request.expression = &quadratic.runtime;
        request.expressionFixed = &quadratic.fixed;
        request.expressionPixel = FormulaParameter::C;
        request.expressionBailout = bailout;
        request.expressionColoring = ExpressionColoring::Smooth;
        backend->resetCancellation();
        bool okay = backend->compute(request);
        bool usedDeep = backend->lastComputeUsedCustomDeepPath();
        bool sameDeep = ordinary.size() == custom.size() && std::memcmp(ordinary.data(), custom.data(), ordinary.size() * sizeof(float)) == 0;
        bool restored = customEngine.escapeRadius() == productionRadius;

        std::fill(ordinary.begin(), ordinary.end(), EMPTYPIXEL);
        std::fill(custom.begin(), custom.end(), EMPTYPIXEL);
        mpf_set_d(re, -0.5);
        mpf_set_ui(im, 0);
        mpf_set_ui(scale, 1);
        const int ordinaryMxit = std::min(mxit, 500);
        ordinaryEngine.Compute(re, im, scale, ordinaryMxit, 0);
        request.mode = ComputeMode::Mandelbrot;
        request.maxIterations = ordinaryMxit;
        backend->resetCancellation();
        bool ordinaryOkay = backend->compute(request);
        bool sameOrdinary = ordinary.size() == custom.size() && std::memcmp(ordinary.data(), custom.data(), ordinary.size() * sizeof(float)) == 0;
        expect(okay && usedDeep && sameDeep && restored && ordinaryOkay && sameOrdinary && !backend->lastComputeUsedCustomDeepPath(), name);
        mpf_clears(re, im, scale, (mpf_ptr)0);
    };
    deepParity("1e13 custom/ordinary byte parity", "-0.743643887037151", "0.13182590420533", "10000000000000", 5000);
    deepParity("1e31 custom/ordinary byte parity", "-0.749139567333446841955467474699747367338762518832278501811", "0.040823298514634751035521346975478853963578400940553676068", "73541770000000000000000000000000", 50000);

    {
        constexpr int W = 12, H = 8, MXIT = 2000000, STEP = 4;
        const char* scaleText = "3831277000000000000000000000000000000000000000000000";
        mp_bitcnt_t precision = 384;
        mpf_t re, im, scale;
        mpf_init2(re, precision);
        mpf_init2(im, precision);
        mpf_init2(scale, precision);
        mpf_set_str(re, testcases::deep51_x, 10);
        mpf_set_str(im, testcases::deep51_y, 10);
        mpf_set_str(scale, scaleText, 10);
        std::vector<float> deep((size_t)W * H, EMPTYPIXEL);
        std::vector<float> oracle((size_t)W * H, EMPTYPIXEL);
        Mandel engine(W, H, MXIT, 1, deep.data());
        engine.setPrecision((int)precision);
        std::unique_ptr<IComputeBackend> backend = createComputeBackend("cpu");
        ComputeRequest request;
        request.mode = ComputeMode::Expression;
        request.cpuEngine = &engine;
        request.centerRe = re;
        request.centerIm = im;
        request.scale = scale;
        request.width = W;
        request.height = H;
        request.sub = 1;
        request.maxIterations = MXIT;
        request.coloringMethod = 0;
        request.iterations = deep.data();
        request.expressionSource = &quadratic.source;
        request.expression = &quadratic.runtime;
        request.expressionFixed = &quadratic.fixed;
        request.expressionPixel = FormulaParameter::C;
        request.expressionBailout = bailout;
        request.expressionColoring = ExpressionColoring::Smooth;
        backend->resetCancellation();
        bool deepOkay = backend->compute(request);
        bool oracleScopeActive = false;
        {
            Mandel::ScopedCustomCompute oracleScope(engine, bailout, CustomDeepZoomOutputAdapter::None);
            oracleScopeActive = oracleScope.active();
            if (oracleScopeActive) engine.ComputeDirect(MXIT, oracle.data(), STEP, 0);
        }
        adaptExpected(oracle, CustomDeepZoomOutputAdapter::SmoothExpression);
        int classMismatch = 0, floorMismatch = 0;
        double maxDelta = 0.0;
        for (int y = 0; y < H; y += STEP) {
            for (int x = 0; x < W; x += STEP) {
                size_t index = (size_t)y * W + x;
                if (isInterior(deep[index]) != isInterior(oracle[index])) {
                    ++classMismatch;
                    continue;
                }
                if (!isInterior(deep[index])) {
                    if (std::floor(deep[index]) != std::floor(oracle[index])) ++floorMismatch;
                    maxDelta = std::max(maxDelta, std::fabs((double)deep[index] - oracle[index]));
                }
            }
        }
        bool oraclePassed = deepOkay && oracleScopeActive && backend->lastComputeUsedCustomDeepPath() && engine.escapeRadius() == productionRadius && classMismatch == 0 && floorMismatch == 0 && maxDelta <= 0.05;
        expect(oraclePassed, "1e51 deep frame matches sampled GMP bailout-4 oracle");
        if (!oraclePassed) { printf("  1e51 oracle class/floor/max=%d/%d/%.9g\n", classMismatch, floorMismatch, maxDelta); }
        mpf_clears(re, im, scale, (mpf_ptr)0);
    }

    {
        constexpr int W = 5, H = 5, MXIT = 50000000;
        mpf_t re, im, scale;
        mpf_init2(re, 256);
        mpf_init2(im, 256);
        mpf_init2(scale, 256);
        mpf_set_str(scale, "1000000000001", 10);
        mpf_set_ui(re, 1);
        mpf_div(re, re, scale);
        mpf_set_ui(im, 1);
        std::vector<float> output((size_t)W * H, EMPTYPIXEL);
        Mandel engine(W, H, MXIT, 1, output.data());
        engine.setPrecision(256);
        std::unique_ptr<IComputeBackend> backend = createComputeBackend("cpu");
        ComputeRequest request;
        request.mode = ComputeMode::Expression;
        request.cpuEngine = &engine;
        request.centerRe = re;
        request.centerIm = im;
        request.scale = scale;
        request.width = W;
        request.height = H;
        request.sub = 1;
        request.maxIterations = MXIT;
        request.coloringMethod = ColoringMethod::ORBIT_TRAP;
        request.iterations = output.data();
        request.expressionSource = &quadratic.source;
        request.expression = &quadratic.runtime;
        request.expressionFixed = &quadratic.fixed;
        request.expressionPixel = FormulaParameter::C;
        request.expressionBailout = bailout;
        request.expressionColoring = ExpressionColoring::OrbitTrap;
        backend->resetCancellation();
        auto future = std::async(std::launch::async, [&] { return backend->compute(request); });
        std::this_thread::sleep_for(std::chrono::milliseconds(50));
        auto cancelStart = std::chrono::steady_clock::now();
        backend->cancel();
        bool completed = future.wait_for(std::chrono::seconds(3)) == std::future_status::ready;
        double cancelSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - cancelStart).count();
        bool result = completed ? future.get() : true;
        int empty = (int)std::count(output.begin(), output.end(), EMPTYPIXEL);
        bool usedCustom = backend->lastComputeUsedCustomDeepPath();
        bool radiusRestored = engine.escapeRadius() == productionRadius;
        expect(completed && cancelSeconds < 3.0 && !result && empty > 0 && usedCustom && radiusRestored, "OrbitTrap GMP fallback cancellation restores state");
        printf("  GMP cancellation completed/result/empty/time/restored="
               "%d/%d/%d/%.3f/%d\n",
               completed ? 1 : 0, result ? 1 : 0, empty, cancelSeconds, radiusRestored ? 1 : 0);
        backend->resetCancellation();

        std::vector<float> ordinary((size_t)W * H, EMPTYPIXEL);
        Mandel ordinaryEngine(W, H, 200, 1, ordinary.data());
        std::fill(output.begin(), output.end(), EMPTYPIXEL);
        mpf_set_d(re, -0.5);
        mpf_set_ui(im, 0);
        mpf_set_ui(scale, 1);
        ordinaryEngine.Compute(re, im, scale, 200, ColoringMethod::ORBIT_TRAP);
        request.mode = ComputeMode::Mandelbrot;
        request.maxIterations = 200;
        request.coloringMethod = ColoringMethod::ORBIT_TRAP;
        bool ordinaryOkay = backend->compute(request);
        bool ordinarySame = std::memcmp(ordinary.data(), output.data(), ordinary.size() * sizeof(float)) == 0;
        expect(ordinaryOkay && ordinarySame && engine.escapeRadius() == productionRadius && !backend->lastComputeUsedCustomDeepPath(), "cancelled Custom trap leaves ordinary trap byte-identical");
        ordinaryEngine.setPrecision(256);
        mpf_set_ui(re, 0);
        mpf_set_ui(im, 0);
        mpf_set_str(scale, "10000000000000", 10);
        ordinaryEngine.Compute(re, im, scale, 200, ColoringMethod::ORBIT_TRAP);
        float ordinaryTrapCenter = ordinary[(H / 2) * W + W / 2];
        expect(ordinaryTrapCenter >= 0.0f, "ordinary deep trap reference keeps its interior color");
        if (ordinaryTrapCenter < 0.0f) printf("  ordinary deep trap center=%.9g\n", ordinaryTrapCenter);
        mpf_clears(re, im, scale, (mpf_ptr)0);
    }

    {
        mpf_t re, im, scale;
        mpf_init2(re, 256);
        mpf_init2(im, 256);
        mpf_init2(scale, 256);
        mpf_set_ui(re, 1);
        mpf_set_ui(im, 0);
        mpf_set_str(scale, "10000000000000", 10);
        formula::ExpressionOrbitSnapshot snapshot;
        snapshot.program = quadratic.runtime;
        snapshot.fixed = quadratic.fixed;
        snapshot.pixelParameter = FormulaParameter::C;
        snapshot.bailout = bailout;
        FormulaContext customFormula = expressionFormula();
        customFormula.slice.pixel = FormulaParameter::C;
        OrbitWorker worker;
        for (ExpressionColoring coloring : {ExpressionColoring::Smooth, ExpressionColoring::Feather, ExpressionColoring::OrbitTrap}) {
            OrbitResult orbit;
            worker.request(re, im, scale, 3, 3, 7, 7, 16, customFormula, std::make_shared<const formula::ExpressionOrbitSnapshot>(snapshot), plan(quadratic, FormulaParameter::C, bailout, coloring, scale, re, im));
            auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(3);
            while (!worker.takeLatest(orbit) && std::chrono::steady_clock::now() < deadline) { std::this_thread::sleep_for(std::chrono::milliseconds(1)); }
            expect(orbit.generation != 0 && orbit.usedGmpQuadratic && orbit.pixelParameter == FormulaParameter::C && orbit.points.size() == (size_t)orbit.iterations + 1 && orbit.escaped && orbit.iterations == 3 && orbit.points.back().re == 5.0f && orbit.points.back().im == 0.0f, "Custom deep GMP orbit uses shared coloring plan");
        }
        mpf_clears(re, im, scale, (mpf_ptr)0);
    }

    printf("=== Custom deep-zoom stage 2\n");
    printf("  Custom bailout=%.9g production bailout=%.9g"
           " direct limit=%.0e\n",
           bailout, productionRadius, formula::CUSTOM_DIRECT_ZOOM_LIMIT);
    printf("  overlap smooth class/floor/max=%d/%d/%.6g"
           " distance class/max/rel=%d/%.6g/%.6g\n",
           smooth.classMismatch, smooth.floorMismatch, smooth.maxDelta, distance.classMismatch, distance.maxDelta, distance.maxRelative);
    printf("  feather/trap crossover max delta=%.6g/%.6g\n", feather.maxDelta, trap.maxDelta);
    printf("  => %s\n\n", failures == 0 ? "PASS" : "CHECK (Custom deep-zoom failure)");
    mpf_clears(capabilityCenterRe, capabilityCenterIm, (mpf_ptr)0);
    return failures == 0 ? 0 : 1;
}

static int runBackendCase() {
    constexpr int W = 96, H = 64, MXIT = 500;
    mpf_set_default_prec(256);
    mpf_t centerRe, centerIm, scale, fixedCRe, fixedCIm;
    mpf_init_set_d(centerRe, -0.5);
    mpf_init_set_ui(centerIm, 0);
    mpf_init_set_ui(scale, 1);
    mpf_init_set_d(fixedCRe, -0.8);
    mpf_init_set_d(fixedCIm, 0.156);
    int failures = 0;

    auto same = [](const std::vector<float>& a, const std::vector<float>& b) { return a.size() == b.size() && std::memcmp(a.data(), b.data(), a.size() * sizeof(float)) == 0; };
    std::unique_ptr<IComputeBackend> backend = createComputeBackend("cpu");
    if (!backend || backend->info().name != "CPU" || backend->info().hardwareAccelerated || backend->info().fallback) ++failures;

    std::vector<float> direct((size_t)W * H), dispatched((size_t)W * H);
    Mandel directMandel(W, H, MXIT, 1, direct.data());
    Mandel backendMandel(W, H, MXIT, 1, dispatched.data());
    directMandel.Compute(centerRe, centerIm, scale, MXIT, 0);
    ComputeRequest request;
    request.mode = ComputeMode::Mandelbrot;
    request.cpuEngine = &backendMandel;
    request.centerRe = centerRe;
    request.centerIm = centerIm;
    request.scale = scale;
    request.width = W;
    request.height = H;
    request.sub = 1;
    request.maxIterations = MXIT;
    request.iterations = dispatched.data();
    backend->resetCancellation();
    if (!backend->compute(request) || !same(direct, dispatched)) ++failures;

    directMandel.ComputeJulia(centerRe, centerIm, scale, fixedCRe, fixedCIm, MXIT, 0);
    request.mode = ComputeMode::Julia;
    request.fixedCRe = fixedCRe;
    request.fixedCIm = fixedCIm;
    backend->resetCancellation();
    if (!backend->compute(request) || !same(direct, dispatched)) ++failures;

    formula::ExpressionProgram expression;
    formula::ExpressionProgram runtimeExpression;
    formula::ExpressionOrbitPlan runtimeExpressionPlan;
    formula::ExpressionJit4 runtimeExpressionJit;
    formula::ExpressionError compileError;
    formula::ExpressionContext fixed;
    fixed.parameters[0] = {0.15, -0.2};
    fixed.parameters[1] = {-0.35, 0.1};
    if (!expression.compile("0.5*z+0.1*sin(c)", &compileError) || !expression.specialize(fixed, FormulaParameter::C, runtimeExpression, &compileError) || !runtimeExpressionPlan.build(runtimeExpression, &compileError)) {
        ++failures;
    } else {
        request.mode = ComputeMode::Expression;
        request.expressionSource = &expression;
        request.expression = &runtimeExpression;
        request.expressionFixed = &fixed;
        request.expressionPlan = runtimeExpressionPlan.profitable() ? &runtimeExpressionPlan : nullptr;
        request.expressionJit = (runtimeExpressionPlan.profitable() ? runtimeExpressionJit.compile(runtimeExpressionPlan) : runtimeExpressionJit.compile(runtimeExpression)) ? &runtimeExpressionJit : nullptr;
        request.expressionPixel = FormulaParameter::C;
        request.expressionBailout = 4.0;
        request.fixedCRe = request.fixedCIm = nullptr;
        for (formula::ExpressionColoring coloring : {formula::ExpressionColoring::Raw, formula::ExpressionColoring::Feather, formula::ExpressionColoring::OrbitTrap}) {
            if (!directMandel.ComputeExpression(centerRe, centerIm, scale, runtimeExpression, fixed, FormulaParameter::C, MXIT, 4.0, coloring, request.expressionJit, request.expressionPlan)) ++failures;
            request.expressionColoring = coloring;
            backend->resetCancellation();
            if (!backend->compute(request) || !same(direct, dispatched) || backend->lastComputeUsedCustomDeepPath()) ++failures;
        }
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

    std::unique_ptr<IComputeBackend> fallback = createComputeBackend("not-a-backend");
    if (!fallback || !fallback->info().fallback || fallback->info().hardwareAccelerated || fallback->info().detail.find("unknown backend") == std::string::npos) ++failures;

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
    if (!warp || warp->info().fallback || warp->info().name != "D3D11 WARP" || warp->info().hardwareAccelerated) {
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
                double difference = std::fabs((double)direct[i] - (double)gpu[i]);
                if (difference > gpuMaxDifference) {
                    gpuMaxDifference = difference;
                    gpuWorst = static_cast<int>(i);
                    gpuWorstCpu = direct[i];
                    gpuWorstGpu = gpu[i];
                }
            }
        }
        if (gpuEmpty != 0 || gpuClassMismatch != 0 || gpuFloorMismatch != 0 || gpuMaxDifference > 0.01) ++failures;

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
                    gpuStressMaxDifference = std::max(gpuStressMaxDifference, std::fabs((double)direct[i] - (double)gpu[i]));
                }
            }
            if (gpuStressClass != 0 || gpuStressFloor != 0 || gpuStressMaxDifference > 0.01) ++failures;
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
                if (layoutCpu[i] != layoutGpu[i] && !(layoutCpu[i] >= 0.0f && layoutGpu[i] >= 0.0f && (int)layoutCpu[i] == (int)layoutGpu[i] && std::fabs(layoutCpu[i] - layoutGpu[i]) <= 0.01f)) ++gpuLayoutMismatch;
            if (gpuLayoutMismatch != 0) ++failures;
        }

        // Unsupported modes must remain byte-identical through the owned CPU
        // fallback instead of returning a partially rendered GPU frame.
        directMandel.ComputeJulia(centerRe, centerIm, scale, fixedCRe, fixedCIm, MXIT, 0);
        request.mode = ComputeMode::Julia;
        request.fixedCRe = fixedCRe;
        request.fixedCIm = fixedCIm;
        warp->resetCancellation();
        if (!warp->compute(request) || !same(direct, gpu) || warp->lastComputeUsedGpuPath()) ++failures;

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
        if (!warp->compute(bailoutRequest) || warp->lastComputeUsedGpuPath() || !same(bailoutCpu, bailoutGpu)) ++failures;
        _putenv_s("MANDEL_BAILOUT", "");

        // Batched AVX2 refinement must preserve the scalar renderer's float
        // bailout-square rounding for custom radii.
        _putenv_s("MANDEL_BAILOUT", "2.1");
        float batchStorage[4] = {};
        Mandel batchMandel(2, 2, 64, 1, batchStorage);
        const double batchRe[] = {0.2522435803501477, -0.75, 0.3, -1.2};
        const double batchIm[] = {0.0, 0.1, 0.5, 0.0};
        float scalarPoints[4], batchPoints[4];
        for (int i = 0; i < 4; ++i) scalarPoints[i] = batchMandel.ComputeShallowPoint(batchRe[i], batchIm[i], 64);
        batchMandel.ComputeShallowPoints(batchRe, batchIm, 4, batchPoints, 64);
        if (std::memcmp(scalarPoints, batchPoints, sizeof(scalarPoints)) != 0) ++failures;
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
        if (!warp->compute(rangeRequest) || warp->lastComputeUsedGpuPath()) ++failures;

        mpf_set_d(centerRe, -0.5);
        mpf_set_ui(centerIm, 0);
        mpf_set_ui(scale, 1);
    }

    // Running cancellation through the hoisted lane-refill path.
    constexpr int CW = 8, CH = 8, CMXIT = 10000000;
    std::vector<float> cancelOutput((size_t)CW * CH, EMPTYPIXEL);
    Mandel cancelMandel(CW, CH, CMXIT, 1, cancelOutput.data());
    formula::ExpressionProgram cancelExpression;
    formula::ExpressionProgram cancelRuntime;
    formula::ExpressionOrbitPlan cancelPlan;
    formula::ExpressionJit4 cancelJit;
    cancelExpression.compile("0.5*z+0.1*sin(c)", &compileError);
    formula::ExpressionContext cancelFixed;
    bool cancelSetup = cancelExpression.specialize(cancelFixed, FormulaParameter::C, cancelRuntime, &compileError) && cancelPlan.build(cancelRuntime, &compileError) && cancelPlan.profitable() && (!VERIFY_JIT || cancelJit.compile(cancelPlan));
    ComputeRequest cancelRequest;
    cancelRequest.mode = ComputeMode::Expression;
    cancelRequest.cpuEngine = &cancelMandel;
    cancelRequest.centerRe = centerRe;
    cancelRequest.centerIm = centerIm;
    cancelRequest.scale = scale;
    cancelRequest.width = CW;
    cancelRequest.height = CH;
    cancelRequest.sub = 1;
    cancelRequest.maxIterations = CMXIT;
    cancelRequest.iterations = cancelOutput.data();
    cancelRequest.expression = &cancelRuntime;
    cancelRequest.expressionFixed = &cancelFixed;
    cancelRequest.expressionPlan = &cancelPlan;
    cancelRequest.expressionJit = VERIFY_JIT ? &cancelJit : nullptr;
    cancelRequest.expressionPixel = FormulaParameter::C;
    cancelRequest.expressionColoring = formula::ExpressionColoring::OrbitTrap;
    backend->resetCancellation();
    auto start = Clock::now();
    std::future<bool> running = std::async(std::launch::async, [&] { return backend->compute(cancelRequest); });
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
    backend->cancel();
    bool cancelResult = running.get();
    double cancelSeconds = since(start);
    backend->resetCancellation();
    int empty = (int)std::count(cancelOutput.begin(), cancelOutput.end(), EMPTYPIXEL);
    if (!cancelSetup || cancelResult || cancelSeconds > 2.0 || empty == 0) ++failures;

    printf("=== compute backend interface\n");
    printf("  backend=%s detail=%s fallback-test=%s\n", backend->info().name.c_str(), backend->info().detail.c_str(), fallback->info().detail.c_str());
    printf("  Mandel/Julia/Expression byte parity; cancel %.3fs empty=%d/%d\n", cancelSeconds, empty, CW * CH);
    printf("  D3D11 WARP: empty=%d class=%d floor=%d max smooth diff=%.6g\n", gpuEmpty, gpuClassMismatch, gpuFloorMismatch, gpuMaxDifference);
    if (gpuWorst >= 0) printf("    worst=(%d,%d) cpu=%.9g gpu=%.9g\n", gpuWorst % W, gpuWorst / W, gpuWorstCpu, gpuWorstGpu);
    printf("    scale1e6 class=%d floor=%d max diff=%.6g; sub5 layout=%d\n", gpuStressClass, gpuStressFloor, gpuStressMaxDifference, gpuLayoutMismatch);
    printf("  => %s\n\n", failures == 0 ? "PASS" : "CHECK (backend failure)");

    mpf_clears(centerRe, centerIm, scale, fixedCRe, fixedCIm, (mpf_ptr)0);
    return failures == 0 ? 0 : 1;
}

static int runGenericDeepBackendCase() {
    using formula::ExpressionContext;
    using formula::ExpressionDeepRenderRequest;
    using formula::ExpressionDeepRenderResult;
    using formula::ExpressionError;
    using formula::ExpressionProgram;

    struct Case {
        const char* name;
        const char* source;
        FormulaParameter pixel;
        formula::Complex fixedC;
        double centerRe;
        double centerIm;
    };
    const Case cases[] = {{"arithmetic-c", "z*z+c+0", FormulaParameter::C, {}, -0.75, 0.0},
                          {"arithmetic-z0", "z*z+c+0", FormulaParameter::InitialZ, {-0.1, 0.2}, 0.0, 0.0},
                          {"sin-exp-c", "0.1*sin(z)+0.01*exp(z)+c", FormulaParameter::C, {}, 0.0, 0.0},
                          {"division-tan-c", "0.1*z/(c+3)+0.01*tan(z)+c", FormulaParameter::C, {}, 0.0, 0.0},
                          {"branch-c",
                           "0.01*log(z+3)+0.01*sqrt(z+4)+"
                           "0.01*pow(z+2,complex(0.5,0))+c",
                           FormulaParameter::C,
                           {},
                           0.0,
                           0.0},
                          {"branch-cut-fallback", "0.001*log(c)", FormulaParameter::C, {}, -1.0, 0.0},
                          {"real-conj-norm-c",
                           "0.1*conj(z)+0.01*norm(z)+"
                           "complex(0.1*real(c),0.1*imag(c))",
                           FormulaParameter::C,
                           {},
                           0.0,
                           0.0},
                          {"burning-ship-c", "0.2*sqr(complex(abs(real(z)),abs(imag(z))))+c", FormulaParameter::C, {}, -0.4, -0.2},
                          {"arg-polar-c",
                           "0.1*arg(z+2)+"
                           "0.1*polar(abs(real(c))+1,arg(c+2))",
                           FormulaParameter::C,
                           {},
                           0.0,
                           0.0},
                          {"sin-z0", "sin(z)+c+0", FormulaParameter::InitialZ, {-0.2, 0.1}, 0.0, 0.0}};

    constexpr int W = 7;
    constexpr int H = 5;
    constexpr int MXIT = 24;
    int failures = 0;
    mpf_t centerRe, centerIm, scale;
    mpf_init2(centerRe, 2048);
    mpf_init2(centerIm, 2048);
    mpf_init2(scale, 2048);
    mpf_set_str(scale, "1e500", 10);
    GenericDeepInfo arithmeticInfo;
    GenericDeepInfo fallbackInfo;

    for (const Case& test : cases) {
        ExpressionProgram canonical;
        ExpressionProgram runtime;
        ExpressionContext fixed;
        fixed.c = test.fixedC;
        ExpressionError error;
        bool okay = canonical.compile(test.source, &error) && canonical.specialize(fixed, test.pixel, runtime, &error);
        mpf_set_d(centerRe, test.centerRe);
        mpf_set_d(centerIm, test.centerIm);
        std::vector<float> backendOutput((size_t)W * H, EMPTYPIXEL);
        std::vector<float> directOutput(backendOutput.size(), EMPTYPIXEL);
        std::vector<float> mpfrOutput(backendOutput.size(), EMPTYPIXEL);
        Mandel engine(W, H, MXIT, 1, backendOutput.data());
        engine.setPrecision(2048);
        std::atomic<float> progress{0.0f};
        std::unique_ptr<IComputeBackend> backend = createComputeBackend("cpu");
        ComputeRequest request;
        request.mode = ComputeMode::Expression;
        request.cpuEngine = &engine;
        request.centerRe = centerRe;
        request.centerIm = centerIm;
        request.scale = scale;
        request.width = W;
        request.height = H;
        request.sub = 1;
        request.maxIterations = MXIT;
        request.iterations = backendOutput.data();
        request.progress = &progress;
        request.expressionSource = &canonical;
        request.expression = &runtime;
        request.expressionFixed = &fixed;
        request.expressionPixel = test.pixel;
        request.expressionBailout = 4.0;
        request.expressionColoring = formula::ExpressionColoring::Raw;
        backend->resetCancellation();
        okay = okay && backend->compute(request);

        ExpressionDeepRenderRequest directRequest;
        directRequest.canonicalProgram = &canonical;
        directRequest.runtimeProgram = &runtime;
        directRequest.center.realMpf = centerRe;
        directRequest.center.imaginaryMpf = centerIm;
        directRequest.scale.mpf = scale;
        directRequest.fixed = fixed;
        directRequest.pixelParameter = test.pixel;
        directRequest.width = W;
        directRequest.height = H;
        directRequest.maxIterations = MXIT;
        directRequest.bailout = 4.0;
        directRequest.output = directOutput.data();
        directRequest.outputCount = directOutput.size();
        directRequest.precision.viewBits = 2048;
        ExpressionDeepRenderResult directResult;
        okay = okay && formula::renderExpressionDeepFrame(directRequest, directResult);

        ExpressionDeepRenderRequest mpfrRequest = directRequest;
        mpfrRequest.output = mpfrOutput.data();
        mpfrRequest.forceMpfrFallbackForVerification = true;
        ExpressionDeepRenderResult mpfrResult;
        okay = okay && formula::renderExpressionDeepFrame(mpfrRequest, mpfrResult);

        const GenericDeepInfo info = backend->lastGenericDeepInfo();
        const bool exact = backendOutput == directOutput && backendOutput == mpfrOutput;
        const bool telemetry = backend->lastComputeUsedGenericDeepPath() && !backend->lastComputeUsedCustomDeepPath() && info.used && info.settled && info.success && info.pixelCount == backendOutput.size() && progress.load(std::memory_order_relaxed) == 1.0f && info.fastPixelCount == directResult.fastPixelCount && info.fallbackPixelCount == directResult.fallbackPixelCount;
        if (!okay || !exact || !telemetry) {
            ++failures;
            printf("  generic backend failed [%s] okay=%d exact=%d "
                   "used/success=%d/%d fast/fallback=%llu/%llu "
                   "status=%s error=%s\n",
                   test.name, okay ? 1 : 0, exact ? 1 : 0, info.used ? 1 : 0, info.success ? 1 : 0, (unsigned long long)info.fastPixelCount, (unsigned long long)info.fallbackPixelCount, info.status.c_str(), info.error.c_str());
        }
        if (strcmp(test.name, "arithmetic-c") == 0) arithmeticInfo = info;
        if (strcmp(test.name, "branch-cut-fallback") == 0) fallbackInfo = info;
    }

    ExpressionProgram genericCanonical;
    ExpressionProgram genericRuntime;
    ExpressionContext genericFixed;
    ExpressionError genericError;
    bool genericReady = genericCanonical.compile("z*z+c+0", &genericError) && genericCanonical.specialize(genericFixed, FormulaParameter::C, genericRuntime, &genericError);
    mpf_set_ui(centerRe, 0);
    mpf_set_ui(centerIm, 0);
    std::vector<float> unsupported((size_t)W * H, EMPTYPIXEL);
    Mandel unsupportedEngine(W, H, MXIT, 1, unsupported.data());
    std::unique_ptr<IComputeBackend> backend = createComputeBackend("cpu");
    ComputeRequest request;
    request.mode = ComputeMode::Expression;
    request.cpuEngine = &unsupportedEngine;
    request.centerRe = centerRe;
    request.centerIm = centerIm;
    request.scale = scale;
    request.width = W;
    request.height = H;
    request.sub = 1;
    request.maxIterations = MXIT;
    request.iterations = unsupported.data();
    request.expressionSource = &genericCanonical;
    request.expression = &genericRuntime;
    request.expressionFixed = &genericFixed;
    request.expressionPixel = FormulaParameter::C;
    request.expressionBailout = 4.0;
    request.expressionColoring = formula::ExpressionColoring::Feather;
    backend->resetCancellation();
    if (!genericReady || backend->compute(request) || backend->lastComputeUsedGenericDeepPath() || std::count(unsupported.begin(), unsupported.end(), EMPTYPIXEL) != (ptrdiff_t)unsupported.size()) {
        ++failures;
        printf("  unsupported deep coloring did not fail empty\n");
    }

    std::vector<float> shallow((size_t)W * H, EMPTYPIXEL);
    std::vector<float> shallowExpected(shallow.size(), EMPTYPIXEL);
    Mandel shallowEngine(W, H, MXIT, 1, shallow.data());
    Mandel shallowOracle(W, H, MXIT, 1, shallowExpected.data());
    mpf_set_d(scale, formula::CUSTOM_DIRECT_ZOOM_LIMIT);
    request.cpuEngine = &shallowEngine;
    request.iterations = shallow.data();
    request.expressionColoring = formula::ExpressionColoring::Raw;
    const bool shallowExpectedOkay = shallowOracle.ComputeExpression(centerRe, centerIm, scale, genericRuntime, genericFixed, FormulaParameter::C, MXIT, 4.0, formula::ExpressionColoring::Raw);
    backend->resetCancellation();
    if (!shallowExpectedOkay || !backend->compute(request) || backend->lastComputeUsedGenericDeepPath() || shallow != shallowExpected) {
        ++failures;
        printf("  scale<=1e12 direct dispatch changed\n");
    }

    ExpressionProgram undefinedCanonical;
    ExpressionProgram undefinedRuntime;
    ExpressionContext undefinedFixed;
    bool undefinedReady = undefinedCanonical.compile("1/(z-z)", &genericError) && undefinedCanonical.specialize(undefinedFixed, FormulaParameter::InitialZ, undefinedRuntime, &genericError);
    mpf_set_str(scale, "1e500", 10);
    std::fill(unsupported.begin(), unsupported.end(), EMPTYPIXEL);
    request.cpuEngine = &unsupportedEngine;
    request.iterations = unsupported.data();
    request.expressionSource = &undefinedCanonical;
    request.expression = &undefinedRuntime;
    request.expressionFixed = &undefinedFixed;
    request.expressionPixel = FormulaParameter::InitialZ;
    backend->resetCancellation();
    const bool undefinedResult = backend->compute(request);
    const GenericDeepInfo undefinedInfo = backend->lastGenericDeepInfo();
    if (!undefinedReady || undefinedResult || !backend->lastComputeUsedGenericDeepPath() || undefinedInfo.success || undefinedInfo.status != "undefined-pixel" || std::count(unsupported.begin(), unsupported.end(), EMPTYPIXEL) != (ptrdiff_t)unsupported.size()) {
        ++failures;
        printf("  undefined generic frame was success-shaped "
               "status=%s error=%s\n",
               undefinedInfo.status.c_str(), undefinedInfo.error.c_str());
    }

    std::unique_ptr<IComputeBackend> warp = createComputeBackend("warp");
    if (!warp || warp->info().fallback) {
        ++failures;
    } else {
        std::fill(unsupported.begin(), unsupported.end(), EMPTYPIXEL);
        request.cpuEngine = &unsupportedEngine;
        request.iterations = unsupported.data();
        request.expressionSource = &genericCanonical;
        request.expression = &genericRuntime;
        request.expressionFixed = &genericFixed;
        request.expressionPixel = FormulaParameter::C;
        warp->resetCancellation();
        const bool warpOkay = warp->compute(request);
        const GenericDeepInfo warpInfo = warp->lastGenericDeepInfo();
        if (!warpOkay || warp->lastComputeUsedGpuPath() || !warp->lastComputeUsedGenericDeepPath() || !warpInfo.used || !warpInfo.success) {
            ++failures;
            printf("  D3D generic deep CPU delegation failed\n");
        }
    }

    ExpressionProgram cancelCanonical;
    ExpressionProgram cancelRuntime;
    ExpressionContext cancelFixed;
    const bool cancelReady = cancelCanonical.compile("z*z+c+0.000001*n", &genericError) && cancelCanonical.specialize(cancelFixed, FormulaParameter::C, cancelRuntime, &genericError);
    constexpr int CW = 9;
    constexpr int CH = 7;
    constexpr int CMAX = 1000000;
    std::vector<float> cancelled((size_t)CW * CH, EMPTYPIXEL);
    Mandel cancelEngine(CW, CH, CMAX, 1, cancelled.data());
    ComputeRequest cancelRequest;
    cancelRequest.mode = ComputeMode::Expression;
    cancelRequest.cpuEngine = &cancelEngine;
    cancelRequest.centerRe = centerRe;
    cancelRequest.centerIm = centerIm;
    cancelRequest.scale = scale;
    cancelRequest.width = CW;
    cancelRequest.height = CH;
    cancelRequest.sub = 1;
    cancelRequest.maxIterations = CMAX;
    cancelRequest.iterations = cancelled.data();
    cancelRequest.expressionSource = &cancelCanonical;
    cancelRequest.expression = &cancelRuntime;
    cancelRequest.expressionFixed = &cancelFixed;
    cancelRequest.expressionPixel = FormulaParameter::C;
    cancelRequest.expressionBailout = 1e100;
    cancelRequest.expressionColoring = formula::ExpressionColoring::Raw;
    backend->resetCancellation();
    auto future = std::async(std::launch::async, [&] { return backend->compute(cancelRequest); });
    for (int wait = 0; wait < 200; ++wait) {
        const GenericDeepInfo info = backend->lastGenericDeepInfo();
        if (info.used && !info.settled) break;
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    const auto cancelStart = std::chrono::steady_clock::now();
    backend->cancel();
    const bool bounded = future.wait_for(std::chrono::seconds(3)) == std::future_status::ready;
    const bool cancelResult = bounded ? future.get() : true;
    const double cancelSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - cancelStart).count();
    const GenericDeepInfo cancelInfo = backend->lastGenericDeepInfo();
    backend->resetCancellation();
    const size_t empty = std::count(cancelled.begin(), cancelled.end(), EMPTYPIXEL);
    if (!cancelReady || !bounded || cancelResult || cancelSeconds >= 3.0 || !cancelInfo.cancelled || cancelInfo.status != "cancelled" || empty == 0) {
        ++failures;
        printf("  generic cancellation failed bounded/result/time/"
               "empty/status=%d/%d/%.3f/%zu/%s\n",
               bounded ? 1 : 0, cancelResult ? 1 : 0, cancelSeconds, empty, cancelInfo.status.c_str());
    }

    const double fallbackRate = arithmeticInfo.pixelCount ? 100.0 * arithmeticInfo.fallbackPixelCount / arithmeticInfo.pixelCount : 100.0;
    if (!arithmeticInfo.success || arithmeticInfo.fastPixelCount == 0 || arithmeticInfo.fallbackPixelCount == arithmeticInfo.pixelCount) ++failures;
    if (!fallbackInfo.success || fallbackInfo.fallbackPixelCount == 0) ++failures;
    printf("=== generic certified deep production backend\n");
    printf("  e500 c/z0 arithmetic, entire, meromorphic, branch, "
           "real/piecewise, Arg/Polar exact MPFR parity\n");
    printf("  arithmetic total/reference/Taylor/fallback "
           "%.6f/%.6f/%.6f/%.6f s; fallback=%llu/%llu (%.2f%%)\n",
           arithmeticInfo.totalSeconds, arithmeticInfo.referenceSeconds, arithmeticInfo.taylorSeconds, arithmeticInfo.fallbackSeconds, (unsigned long long)arithmeticInfo.fallbackPixelCount, (unsigned long long)arithmeticInfo.pixelCount, fallbackRate);
    printf("  branch-cut fallback=%llu/%llu\n", (unsigned long long)fallbackInfo.fallbackPixelCount, (unsigned long long)fallbackInfo.pixelCount);
    printf("  cancellation %.3f s empty=%zu/%zu\n", cancelSeconds, empty, cancelled.size());
    printf("  => %s\n\n", failures == 0 ? "PASS" : "CHECK (generic deep integration failure)");
    mpf_clears(centerRe, centerIm, scale, (mpf_ptr)0);
    return failures == 0 ? 0 : 1;
}

static int runGpuBenchmarkCase(int width, int height) {
    std::unique_ptr<IComputeBackend> gpu = createComputeBackend("gpu");
    printf("=== D3D11 hardware GPU benchmark\n");
    if (!gpu || gpu->info().fallback || !gpu->info().hardwareAccelerated) {
        printf("  SKIP: %s\n\n", gpu ? gpu->info().detail.c_str() : "backend creation failed");
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
            maxDifference = std::max(maxDifference, std::fabs((double)cpuOutput[i] - (double)gpuOutput[i]));
        }
    }
    double speedup = gpuSeconds > 0.0 ? cpuSeconds / gpuSeconds : 0.0;
    bool accurate = warmupOk && warmupUsedGpu && gpuOk && measuredUsedGpu && empty == 0 && classMismatch == 0 && floorMismatch == 0 && maxDifference <= 0.01;
    bool fastEnough = speedup >= 5.0;
    printf("  backend=%s  %s\n", gpu->info().name.c_str(), gpu->info().detail.c_str());
    printf("  frame=%dx%d mxit=%d CPU/GPU=%.3f/%.3f s speedup=%.2fx\n", width, height, mxit, cpuSeconds, gpuSeconds, speedup);
    printf("  empty=%d class=%d floor=%d max smooth diff=%.6g\n", empty, classMismatch, floorMismatch, maxDifference);
    printf("  => %s\n\n", accurate && fastEnough ? "PASS" : (!accurate ? "CHECK (accuracy failure)" : "CHECK (below 5x acceptance threshold)"));

    mpf_clears(centerRe, centerIm, scale, (mpf_ptr)0);
    return accurate && fastEnough ? 0 : 1;
}

static std::string pow10(int n) {
    std::string s = "1";
    s.append(n, '0');
    return s;
}

struct DoubleDouble {
    double high = 0.0;
    double low = 0.0;
};

static DoubleDouble normalizeDoubleDouble(double high, double low) {
    const double sum = high + low;
    return {sum, low - (sum - high)};
}

static DoubleDouble twoSumDoubleDouble(double left, double right) {
    const double sum = left + right;
    const double virtualRight = sum - left;
    return {sum, (left - (sum - virtualRight)) + (right - virtualRight)};
}

static DoubleDouble addDoubleDouble(DoubleDouble left, DoubleDouble right) {
    DoubleDouble high = twoSumDoubleDouble(left.high, right.high);
    const DoubleDouble low = twoSumDoubleDouble(left.low, right.low);
    high.low += low.high;
    high = normalizeDoubleDouble(high.high, high.low);
    high.low += low.low;
    return normalizeDoubleDouble(high.high, high.low);
}

static DoubleDouble negateDoubleDouble(DoubleDouble value) {
    return {-value.high, -value.low};
}

static DoubleDouble subtractDoubleDouble(DoubleDouble left, DoubleDouble right) {
    return addDoubleDouble(left, negateDoubleDouble(right));
}

static DoubleDouble multiplyDoubleDouble(DoubleDouble left, DoubleDouble right) {
    const double product = left.high * right.high;
    double error = std::fma(left.high, right.high, -product);
    error += left.high * right.low + left.low * right.high;
    DoubleDouble result = normalizeDoubleDouble(product, error);
    result.low += left.low * right.low;
    return normalizeDoubleDouble(result.high, result.low);
}

static DoubleDouble multiplyDoubleDouble(DoubleDouble value, double scalar) {
    const double product = value.high * scalar;
    const double error = std::fma(value.high, scalar, -product) + value.low * scalar;
    return normalizeDoubleDouble(product, error);
}

static DoubleDouble absoluteDoubleDouble(DoubleDouble value) {
    return value.high < 0.0 || (value.high == 0.0 && value.low < 0.0) ? negateDoubleDouble(value) : value;
}

static bool greaterDoubleDouble(DoubleDouble left, double right) {
    return left.high > right || (left.high == right && left.low > 0.0);
}

static DoubleDouble splitMpfr(mpfr_srcptr value, mpfr_ptr temporary) {
    DoubleDouble result;
    result.high = mpfr_get_d(value, MPFR_RNDN);
    mpfr_set_d(temporary, result.high, MPFR_RNDN);
    mpfr_sub(temporary, value, temporary, MPFR_RNDN);
    result.low = mpfr_get_d(temporary, MPFR_RNDN);
    return result;
}

static bool renderBurningShipDoubleDouble(const char* centerReal, const char* centerImaginary, const char* scaleText, int width, int height, int maxIterations, double bailout, std::vector<float>& output, double& seconds) {
    constexpr mpfr_prec_t precision = 512;
    mpfr_t centerRe, centerIm, scale, dxHalf, dyHalf, temporary, coordinate;
    mpfr_inits2(precision, centerRe, centerIm, scale, dxHalf, dyHalf, temporary, coordinate, (mpfr_ptr)0);
    const bool parsed = mpfr_set_str(centerRe, centerReal, 10, MPFR_RNDN) == 0 && mpfr_set_str(centerIm, centerImaginary, 10, MPFR_RNDN) == 0 && mpfr_set_str(scale, scaleText, 10, MPFR_RNDN) == 0 && mpfr_sgn(scale) > 0;
    if (!parsed) {
        mpfr_clears(centerRe, centerIm, scale, dxHalf, dyHalf, temporary, coordinate, (mpfr_ptr)0);
        return false;
    }
    mpfr_mul_ui(temporary, scale, static_cast<unsigned long>(width - 1), MPFR_RNDN);
    mpfr_ui_div(dxHalf, 2, temporary, MPFR_RNDN);
    mpfr_mul_ui(temporary, scale, static_cast<unsigned long>(width), MPFR_RNDN);
    mpfr_mul_ui(temporary, temporary, static_cast<unsigned long>(height - 1), MPFR_RNDN);
    mpfr_ui_div(dyHalf, static_cast<unsigned long>(height), temporary, MPFR_RNDN);
    mpfr_mul_ui(dyHalf, dyHalf, 2, MPFR_RNDN);

    std::vector<DoubleDouble> realCoordinates(static_cast<size_t>(width));
    std::vector<DoubleDouble> imaginaryCoordinates(static_cast<size_t>(height));
    for (int x = 0; x < width; ++x) {
        const long centered = static_cast<long>(2LL * x - (width - 1LL));
        mpfr_mul_si(coordinate, dxHalf, centered, MPFR_RNDN);
        mpfr_add(coordinate, coordinate, centerRe, MPFR_RNDN);
        realCoordinates[static_cast<size_t>(x)] = splitMpfr(coordinate, temporary);
    }
    for (int y = 0; y < height; ++y) {
        const long centered = static_cast<long>(2LL * y - (height - 1LL));
        mpfr_mul_si(coordinate, dyHalf, centered, MPFR_RNDN);
        mpfr_add(coordinate, coordinate, centerIm, MPFR_RNDN);
        imaginaryCoordinates[static_cast<size_t>(y)] = splitMpfr(coordinate, temporary);
    }
    mpfr_clears(centerRe, centerIm, scale, dxHalf, dyHalf, temporary, coordinate, (mpfr_ptr)0);

    output.assign(static_cast<size_t>(width) * height, formula::ExpressionDeepInteriorPixel);
    const double bailoutSquared = bailout * bailout;
    const Clock::time_point start = Clock::now();
#pragma omp parallel for schedule(dynamic, 8)
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            DoubleDouble real{};
            DoubleDouble imaginary{};
            const DoubleDouble constantReal = realCoordinates[static_cast<size_t>(x)];
            const DoubleDouble constantImaginary = imaginaryCoordinates[static_cast<size_t>(y)];
            float result = formula::ExpressionDeepInteriorPixel;
            for (int iteration = 0; iteration < maxIterations; ++iteration) {
                const DoubleDouble absoluteReal = absoluteDoubleDouble(real);
                const DoubleDouble absoluteImaginary = absoluteDoubleDouble(imaginary);
                const DoubleDouble realSquared = multiplyDoubleDouble(absoluteReal, absoluteReal);
                const DoubleDouble imaginarySquared = multiplyDoubleDouble(absoluteImaginary, absoluteImaginary);
                const DoubleDouble product = multiplyDoubleDouble(absoluteReal, absoluteImaginary);
                real = addDoubleDouble(subtractDoubleDouble(realSquared, imaginarySquared), constantReal);
                imaginary = addDoubleDouble(multiplyDoubleDouble(product, 2.0), constantImaginary);
                const DoubleDouble magnitudeSquared = addDoubleDouble(multiplyDoubleDouble(real, real), multiplyDoubleDouble(imaginary, imaginary));
                if (greaterDoubleDouble(magnitudeSquared, bailoutSquared)) {
                    result = static_cast<float>(iteration + 1);
                    break;
                }
            }
            output[static_cast<size_t>(y) * width + x] = result;
        }
    }
    seconds = since(start);
    return true;
}

static bool renderBurningShipMpf(const char* centerReal, const char* centerImaginary, const char* scaleText, int width, int height, int maxIterations, double bailout, mp_bitcnt_t precision, std::vector<float>& output, double& seconds) {
    mpf_t centerRe, centerIm, scale, dxHalf, dyHalf, temporary;
    mpf_init2(centerRe, precision);
    mpf_init2(centerIm, precision);
    mpf_init2(scale, precision);
    mpf_init2(dxHalf, precision);
    mpf_init2(dyHalf, precision);
    mpf_init2(temporary, precision);
    const bool parsed = mpf_set_str(centerRe, centerReal, 10) == 0 && mpf_set_str(centerIm, centerImaginary, 10) == 0 && mpf_set_str(scale, scaleText, 10) == 0 && mpf_sgn(scale) > 0;
    if (!parsed) {
        mpf_clears(centerRe, centerIm, scale, dxHalf, dyHalf, temporary, (mpf_ptr)0);
        return false;
    }
    mpf_mul_ui(temporary, scale, static_cast<unsigned long>(width - 1));
    mpf_ui_div(dxHalf, 2, temporary);
    mpf_mul_ui(temporary, scale, static_cast<unsigned long>(width));
    mpf_mul_ui(temporary, temporary, static_cast<unsigned long>(height - 1));
    mpf_ui_div(dyHalf, static_cast<unsigned long>(height), temporary);
    mpf_mul_ui(dyHalf, dyHalf, 2);

    output.assign(static_cast<size_t>(width) * height, formula::ExpressionDeepInteriorPixel);
    const double bailoutSquared = bailout * bailout;
    const Clock::time_point start = Clock::now();
#pragma omp parallel
    {
        mpf_t constantRe, constantIm, real, imaginary, absoluteReal, absoluteImaginary, realSquared, imaginarySquared, product, nextReal, nextImaginary, magnitudeSquared, coordinate;
        mpf_init2(constantRe, precision);
        mpf_init2(constantIm, precision);
        mpf_init2(real, precision);
        mpf_init2(imaginary, precision);
        mpf_init2(absoluteReal, precision);
        mpf_init2(absoluteImaginary, precision);
        mpf_init2(realSquared, precision);
        mpf_init2(imaginarySquared, precision);
        mpf_init2(product, precision);
        mpf_init2(nextReal, precision);
        mpf_init2(nextImaginary, precision);
        mpf_init2(magnitudeSquared, precision);
        mpf_init2(coordinate, precision);
#pragma omp for schedule(dynamic, 8)
        for (int y = 0; y < height; ++y) {
            const long centeredY = static_cast<long>(2LL * y - (height - 1LL));
            mpf_mul_ui(coordinate, dyHalf, static_cast<unsigned long>(std::labs(centeredY)));
            if (centeredY < 0) mpf_neg(coordinate, coordinate);
            mpf_add(constantIm, coordinate, centerIm);
            for (int x = 0; x < width; ++x) {
                const long centeredX = static_cast<long>(2LL * x - (width - 1LL));
                mpf_mul_ui(coordinate, dxHalf, static_cast<unsigned long>(std::labs(centeredX)));
                if (centeredX < 0) mpf_neg(coordinate, coordinate);
                mpf_add(constantRe, coordinate, centerRe);
                mpf_set_ui(real, 0);
                mpf_set_ui(imaginary, 0);
                float result = formula::ExpressionDeepInteriorPixel;
                for (int iteration = 0; iteration < maxIterations; ++iteration) {
                    mpf_abs(absoluteReal, real);
                    mpf_abs(absoluteImaginary, imaginary);
                    mpf_mul(realSquared, absoluteReal, absoluteReal);
                    mpf_mul(imaginarySquared, absoluteImaginary, absoluteImaginary);
                    mpf_sub(nextReal, realSquared, imaginarySquared);
                    mpf_add(nextReal, nextReal, constantRe);
                    mpf_mul(product, absoluteReal, absoluteImaginary);
                    mpf_mul_ui(nextImaginary, product, 2);
                    mpf_add(nextImaginary, nextImaginary, constantIm);
                    mpf_set(real, nextReal);
                    mpf_set(imaginary, nextImaginary);
                    mpf_mul(realSquared, real, real);
                    mpf_mul(imaginarySquared, imaginary, imaginary);
                    mpf_add(magnitudeSquared, realSquared, imaginarySquared);
                    if (mpf_cmp_d(magnitudeSquared, bailoutSquared) > 0) {
                        result = static_cast<float>(iteration + 1);
                        break;
                    }
                }
                output[static_cast<size_t>(y) * width + x] = result;
            }
        }
        mpf_clears(constantRe, constantIm, real, imaginary, absoluteReal, absoluteImaginary, realSquared, imaginarySquared, product, nextReal, nextImaginary, magnitudeSquared, coordinate, (mpf_ptr)0);
    }
    seconds = since(start);
    mpf_clears(centerRe, centerIm, scale, dxHalf, dyHalf, temporary, (mpf_ptr)0);
    return true;
}

static bool renderBurningShipBigFixed(const char* centerReal, const char* centerImaginary, const char* scaleText, int width, int height, int maxIterations, double bailout, int limbs, std::vector<float>& output, double& seconds) {
    const mp_bitcnt_t conversionPrecision = std::max<mp_bitcnt_t>(512, static_cast<mp_bitcnt_t>(64 * (limbs + 2)));
    mpf_t centerRe, centerIm, scale, dxHalf, dyHalf, temporary, coordinate;
    mpf_init2(centerRe, conversionPrecision);
    mpf_init2(centerIm, conversionPrecision);
    mpf_init2(scale, conversionPrecision);
    mpf_init2(dxHalf, conversionPrecision);
    mpf_init2(dyHalf, conversionPrecision);
    mpf_init2(temporary, conversionPrecision);
    mpf_init2(coordinate, conversionPrecision);
    const bool parsed = mpf_set_str(centerRe, centerReal, 10) == 0 && mpf_set_str(centerIm, centerImaginary, 10) == 0 && mpf_set_str(scale, scaleText, 10) == 0 && mpf_sgn(scale) > 0;
    if (!parsed) {
        mpf_clears(centerRe, centerIm, scale, dxHalf, dyHalf, temporary, coordinate, (mpf_ptr)0);
        return false;
    }
    mpf_mul_ui(temporary, scale, static_cast<unsigned long>(width - 1));
    mpf_ui_div(dxHalf, 2, temporary);
    mpf_mul_ui(temporary, scale, static_cast<unsigned long>(width));
    mpf_mul_ui(temporary, temporary, static_cast<unsigned long>(height - 1));
    mpf_ui_div(dyHalf, static_cast<unsigned long>(height), temporary);
    mpf_mul_ui(dyHalf, dyHalf, 2);

    std::vector<BigFixed> realCoordinates;
    std::vector<BigFixed> imaginaryCoordinates;
    realCoordinates.reserve(static_cast<size_t>(width));
    imaginaryCoordinates.reserve(static_cast<size_t>(height));
    for (int x = 0; x < width; ++x) {
        const long centered = static_cast<long>(2LL * x - (width - 1LL));
        mpf_mul_ui(coordinate, dxHalf, static_cast<unsigned long>(std::labs(centered)));
        if (centered < 0) mpf_neg(coordinate, coordinate);
        mpf_add(coordinate, coordinate, centerRe);
        realCoordinates.emplace_back(limbs);
        bf_from_mpf(realCoordinates.back(), coordinate, limbs);
    }
    for (int y = 0; y < height; ++y) {
        const long centered = static_cast<long>(2LL * y - (height - 1LL));
        mpf_mul_ui(coordinate, dyHalf, static_cast<unsigned long>(std::labs(centered)));
        if (centered < 0) mpf_neg(coordinate, coordinate);
        mpf_add(coordinate, coordinate, centerIm);
        imaginaryCoordinates.emplace_back(limbs);
        bf_from_mpf(imaginaryCoordinates.back(), coordinate, limbs);
    }
    mpf_clears(centerRe, centerIm, scale, dxHalf, dyHalf, temporary, coordinate, (mpf_ptr)0);

    output.assign(static_cast<size_t>(width) * height, formula::ExpressionDeepInteriorPixel);
    const unsigned long long bailoutSquaredInteger = static_cast<unsigned long long>(std::ceil(bailout * bailout));
    const Clock::time_point start = Clock::now();
#pragma omp parallel
    {
        std::vector<uint64_t> scratch(static_cast<size_t>(2 * limbs));
        BigFixed real(limbs), imaginary(limbs), realSquare(limbs), imaginarySquare(limbs), product(limbs), realPart(limbs), imaginaryPart(limbs), nextReal(limbs), nextImaginary(limbs), magnitudeSquared(limbs), threshold(limbs);
        threshold.setInt(static_cast<long long>(bailoutSquaredInteger));
#pragma omp for schedule(dynamic, 8)
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                real.setZero();
                imaginary.setZero();
                realSquare.setZero();
                imaginarySquare.setZero();
                const BigFixed& constantReal = realCoordinates[static_cast<size_t>(x)];
                const BigFixed& constantImaginary = imaginaryCoordinates[static_cast<size_t>(y)];
                float result = formula::ExpressionDeepInteriorPixel;
                for (int iteration = 0; iteration < maxIterations; ++iteration) {
                    bf_sub(realPart, realSquare, imaginarySquare);
                    bf_mul(product, real, imaginary, scratch.data());
                    if (!product.isZero()) product.sign = 1;
                    bf_add(imaginaryPart, product, product);
                    bf_add(nextReal, realPart, constantReal);
                    bf_add(nextImaginary, imaginaryPart, constantImaginary);
                    bf_sqr(realSquare, nextReal, scratch.data());
                    bf_sqr(imaginarySquare, nextImaginary, scratch.data());
                    bf_add(magnitudeSquared, realSquare, imaginarySquare);
                    real.m.swap(nextReal.m);
                    real.sign = nextReal.sign;
                    imaginary.m.swap(nextImaginary.m);
                    imaginary.sign = nextImaginary.sign;
                    if (bf_mag_cmp(magnitudeSquared.m.data(), threshold.m.data(), limbs) > 0) {
                        result = static_cast<float>(iteration + 1);
                        break;
                    }
                }
                output[static_cast<size_t>(y) * width + x] = result;
            }
        }
    }
    seconds = since(start);
    return true;
}

static bool renderExpressionCentered4(const formula::ExpressionProgram& canonical, const formula::ExpressionProgram& runtime, const formula::ExpressionContext& fixed, const char* centerReal, const char* centerImaginary, const char* scaleText, int width, int height, int referenceX, int referenceY, int maxIterations, double bailout, std::vector<float>& output, double& seconds, size_t& failedPixels) {
    formula::ExpressionReferenceBuildRequest referenceRequest;
    referenceRequest.canonicalProgram = &canonical;
    referenceRequest.runtimeProgram = &runtime;
    referenceRequest.pixelParameter = FormulaParameter::C;
    mpf_t referenceReal, referenceImaginary, referenceScale, referenceTemporary, referenceDxHalf, referenceDyHalf;
    mpf_init2(referenceReal, 512);
    mpf_init2(referenceImaginary, 512);
    mpf_init2(referenceScale, 512);
    mpf_init2(referenceTemporary, 512);
    mpf_init2(referenceDxHalf, 512);
    mpf_init2(referenceDyHalf, 512);
    bool referenceParsed = mpf_set_str(referenceReal, centerReal, 10) == 0 && mpf_set_str(referenceImaginary, centerImaginary, 10) == 0 && mpf_set_str(referenceScale, scaleText, 10) == 0 && mpf_sgn(referenceScale) > 0;
    if (referenceParsed) {
        mpf_mul_ui(referenceTemporary, referenceScale, static_cast<unsigned long>(width - 1));
        mpf_ui_div(referenceDxHalf, 2, referenceTemporary);
        mpf_mul_ui(referenceTemporary, referenceScale, static_cast<unsigned long>(width));
        mpf_mul_ui(referenceTemporary, referenceTemporary, static_cast<unsigned long>(height - 1));
        mpf_ui_div(referenceDyHalf, static_cast<unsigned long>(height), referenceTemporary);
        mpf_mul_ui(referenceDyHalf, referenceDyHalf, 2);
        const long centeredX = static_cast<long>(2LL * referenceX - (width - 1LL));
        const long centeredY = static_cast<long>(2LL * referenceY - (height - 1LL));
        mpf_mul_ui(referenceTemporary, referenceDxHalf, static_cast<unsigned long>(std::labs(centeredX)));
        if (centeredX < 0) mpf_neg(referenceTemporary, referenceTemporary);
        mpf_add(referenceReal, referenceReal, referenceTemporary);
        mpf_mul_ui(referenceTemporary, referenceDyHalf, static_cast<unsigned long>(std::labs(centeredY)));
        if (centeredY < 0) mpf_neg(referenceTemporary, referenceTemporary);
        mpf_add(referenceImaginary, referenceImaginary, referenceTemporary);
        referenceRequest.center.realMpf = referenceReal;
        referenceRequest.center.imaginaryMpf = referenceImaginary;
    }
    referenceRequest.fixed = fixed;
    referenceRequest.maxIterations = maxIterations;
    referenceRequest.bailout = bailout;
    referenceRequest.precision.requestedBits = 512;
    referenceRequest.precision.minimumBits = 53;
    referenceRequest.precision.guardBits = 0;
    referenceRequest.precision.maximumBits = 4096;
    formula::ExpressionReferenceOrbitResult reference;
    const bool referenceBuilt = referenceParsed && formula::buildExpressionReferenceOrbit(referenceRequest, reference);
    mpf_clears(referenceReal, referenceImaginary, referenceScale, referenceTemporary, referenceDxHalf, referenceDyHalf, (mpf_ptr)0);
    if (!referenceBuilt || !reference.valid || reference.samples.size() != static_cast<size_t>(maxIterations) || reference.escaped || reference.undefined) {
        printf("  centered VM4 reference invalid samples=%zu escaped=%d undefined=%d error=%s\n", reference.samples.size(), reference.escaped ? 1 : 0, reference.undefined ? 1 : 0, reference.error.c_str());
        return false;
    }

    auto reconstruct = [](const formula::ScaledComplexShadow& primary, const formula::ScaledComplexShadow& defect, formula::Complex& outputValue) {
        formula::ScaledComplexValue scaled;
        return formula::makeScaledComplexValue(primary, defect, scaled) == formula::ScaledArithmeticStatus::Success && formula::scaledValueToDouble(scaled, outputValue);
    };
    formula::Complex referenceC;
    formula::Complex referenceZ0;
    if (!reconstruct(reference.c, reference.cDefect, referenceC) || !reconstruct(reference.z0, reference.z0Defect, referenceZ0)) return false;
    std::vector<formula::Complex> referenceZ(static_cast<size_t>(maxIterations));
    std::vector<formula::Complex> referenceNodes(reference.tape.size());
    for (size_t node = 0; node < reference.tape.size(); ++node)
        if (!reconstruct(reference.tape[node].output, reference.tape[node].outputDefect, referenceNodes[node])) return false;
    for (int iteration = 0; iteration < maxIterations; ++iteration) {
        const formula::ExpressionReferenceSample& sample = reference.samples[static_cast<size_t>(iteration)];
        if (!reconstruct(sample.z, sample.zDefect, referenceZ[static_cast<size_t>(iteration)])) return false;
    }

    constexpr mpfr_prec_t GeometryPrecision = 512;
    mpfr_t scale, temporary, dxHalf, dyHalf, coordinate;
    mpfr_inits2(GeometryPrecision, scale, temporary, dxHalf, dyHalf, coordinate, (mpfr_ptr)0);
    const bool parsed = mpfr_set_str(scale, scaleText, 10, MPFR_RNDN) == 0 && mpfr_sgn(scale) > 0;
    if (!parsed) {
        mpfr_clears(scale, temporary, dxHalf, dyHalf, coordinate, (mpfr_ptr)0);
        return false;
    }
    mpfr_mul_ui(temporary, scale, static_cast<unsigned long>(width - 1), MPFR_RNDN);
    mpfr_ui_div(dxHalf, 2, temporary, MPFR_RNDN);
    mpfr_mul_ui(temporary, scale, static_cast<unsigned long>(width), MPFR_RNDN);
    mpfr_mul_ui(temporary, temporary, static_cast<unsigned long>(height - 1), MPFR_RNDN);
    mpfr_ui_div(dyHalf, static_cast<unsigned long>(height), temporary, MPFR_RNDN);
    mpfr_mul_ui(dyHalf, dyHalf, 2, MPFR_RNDN);
    std::vector<double> xOffsets(static_cast<size_t>(width));
    std::vector<double> yOffsets(static_cast<size_t>(height));
    for (int x = 0; x < width; ++x) {
        const long centered = static_cast<long>(2LL * x - (width - 1LL));
        mpfr_mul_si(coordinate, dxHalf, centered, MPFR_RNDN);
        xOffsets[static_cast<size_t>(x)] = mpfr_get_d(coordinate, MPFR_RNDN);
    }
    for (int y = 0; y < height; ++y) {
        const long centered = static_cast<long>(2LL * y - (height - 1LL));
        mpfr_mul_si(coordinate, dyHalf, centered, MPFR_RNDN);
        yOffsets[static_cast<size_t>(y)] = mpfr_get_d(coordinate, MPFR_RNDN);
    }
    mpfr_clears(scale, temporary, dxHalf, dyHalf, coordinate, (mpfr_ptr)0);

    output.assign(static_cast<size_t>(width) * height, formula::ExpressionDeepInteriorPixel);
    std::atomic<size_t> failures{0};
    const Clock::time_point start = Clock::now();
#pragma omp parallel for schedule(dynamic, 1)
    for (int y = 0; y < height; ++y) {
        formula::ExpressionContext referenceContext = fixed;
        referenceContext.c = referenceC;
        referenceContext.z0 = referenceZ0;
        for (int xBase = 0; xBase < width; xBase += 4) {
            const int lanes = std::min(4, width - xBase);
            std::array<formula::ExpressionDeltaContext, 4> deltas{};
            std::array<formula::Complex, 4> stateDelta{};
            std::array<bool, 4> active{};
            std::array<int, 4> pixelIndices{};
            for (int lane = 0; lane < lanes; ++lane) {
                const int x = xBase + lane;
                pixelIndices[lane] = y * width + x;
                deltas[lane].c = {xOffsets[static_cast<size_t>(x)] - xOffsets[static_cast<size_t>(referenceX)], yOffsets[static_cast<size_t>(y)] - yOffsets[static_cast<size_t>(referenceY)]};
                active[lane] = true;
            }
            int activeCount = lanes;
            for (int iteration = 0; iteration < maxIterations && activeCount > 0; ++iteration) {
                referenceContext.iteration = iteration;
                referenceContext.z = referenceZ[static_cast<size_t>(iteration)];
                for (int lane = 0; lane < 4; ++lane) deltas[lane].z = stateDelta[lane];
                formula::ExpressionCenteredResult evaluated[4];
                const formula::ExpressionReferenceSample& sample = reference.samples[static_cast<size_t>(iteration)];
                const formula::Complex* nodeBases = referenceNodes.data() + static_cast<size_t>(sample.tapeOffset);
                if (!formula::ExpressionCenteredEvaluator::evaluate4WithNodeBases(runtime, referenceContext, nodeBases, deltas.data(), evaluated)) {
                    for (int lane = 0; lane < lanes; ++lane)
                        if (active[lane]) {
                            active[lane] = false;
                            ++failures;
                            --activeCount;
                        }
                    break;
                }
                const formula::Complex nextReference = nodeBases[sample.rootNode];
                for (int lane = 0; lane < lanes; ++lane) {
                    if (!active[lane]) continue;
                    if (!evaluated[lane].success()) {
                        active[lane] = false;
                        ++failures;
                        --activeCount;
                        continue;
                    }
                    stateDelta[lane] = evaluated[lane].delta;
                    const formula::Complex actual = nextReference + stateDelta[lane];
                    if (!std::isfinite(actual.real()) || !std::isfinite(actual.imag()) || std::hypot(actual.real(), actual.imag()) > bailout) {
                        output[static_cast<size_t>(pixelIndices[lane])] = static_cast<float>(iteration + 1);
                        active[lane] = false;
                        --activeCount;
                    }
                }
            }
        }
    }
    seconds = since(start);
    failedPixels = failures.load(std::memory_order_relaxed);
    return true;
}

static int runCustomSlowdownCase(int width, int height) {
    using formula::ExpressionContext;
    using formula::ExpressionDeepFallbackReason;
    using formula::ExpressionDeepRenderRequest;
    using formula::ExpressionDeepRenderResult;
    using formula::ExpressionError;
    using formula::ExpressionJit4;
    using formula::ExpressionOrbitPlan;
    using formula::ExpressionProgram;

    const int maxIterations = [] {
        const char* value = getenv("MANDEL_CUSTOM_SLOW_MXIT");
        return value ? std::clamp(atoi(value), 1, 5000000) : 2000;
    }();
    constexpr double bailout = 100.0;
    const char* source = getenv("MANDEL_CUSTOM_SLOW_SOURCE");
    if (!source) source = "sqr(complex(abs(real(z)),abs(imag(z))))+c";
    const char* centerReal = getenv("MANDEL_CUSTOM_SLOW_CX");
    if (!centerReal) centerReal = "-1.013951002213813310632862698887121834129";
    const char* centerImaginary = getenv("MANDEL_CUSTOM_SLOW_CY");
    if (!centerImaginary) centerImaginary = "-0.7988691125646760914741501921763298573252";
    const char* scaleText = getenv("MANDEL_CUSTOM_SLOW_ZOOM");
    if (!scaleText) scaleText = "3.245427860252859436180864346639097915255e12";

    ExpressionProgram canonical;
    ExpressionProgram runtime;
    ExpressionError error;
    ExpressionContext fixed;
    if (!canonical.compile(source, &error) || !canonical.specialize(fixed, FormulaParameter::C, runtime, &error)) {
        printf("custom slowdown compile failed: %s\n", error.message.c_str());
        return 1;
    }
    printf("  formula source=%s instructions=%zu piecewise=%d\n", source, runtime.instructionCount(), runtime.piecewiseQuadraticKind() == formula::ExpressionPiecewiseQuadraticKind::BurningShip ? 1 : 0);
    const bool burningShipResearchEligible = runtime.piecewiseQuadraticKind() == formula::ExpressionPiecewiseQuadraticKind::BurningShip;
    if (!burningShipResearchEligible && (getenv("MANDEL_CUSTOM_SLOW_DD") || getenv("MANDEL_CUSTOM_SLOW_MPF_BITS") || getenv("MANDEL_CUSTOM_SLOW_BF_LIMBS"))) printf("  DD/mpf/BigFixed research lanes SKIP: configured formula is not exact Burning Ship\n");

    const size_t pixelCount = static_cast<size_t>(width) * static_cast<size_t>(height);
    auto makeRequest = [&](std::vector<float>& output) {
        output.assign(pixelCount, formula::ExpressionDeepEmptyPixel);
        ExpressionDeepRenderRequest request;
        request.canonicalProgram = &canonical;
        request.runtimeProgram = &runtime;
        request.center.realDecimal = centerReal;
        request.center.imaginaryDecimal = centerImaginary;
        request.scale.decimal = scaleText;
        request.fixed = fixed;
        request.pixelParameter = FormulaParameter::C;
        request.width = width;
        request.height = height;
        request.maxIterations = maxIterations;
        request.bailout = bailout;
        request.output = output.data();
        request.outputCount = output.size();
        const char* minimumBits = getenv("MANDEL_CUSTOM_SLOW_BITS");
        request.precision.minimumBits = minimumBits ? std::clamp<mpfr_prec_t>(static_cast<mpfr_prec_t>(strtoll(minimumBits, nullptr, 10)), 53, 4096) : 128;
        const char* guardBits = getenv("MANDEL_CUSTOM_SLOW_GUARD");
        request.precision.guardBits = guardBits ? std::clamp<mpfr_prec_t>(static_cast<mpfr_prec_t>(strtoll(guardBits, nullptr, 10)), 0, 4096) : 64;
        request.precision.maximumBits = 4096;
        const char* fallbackGuardBits = getenv("MANDEL_CUSTOM_SLOW_FALLBACK_GUARD");
        request.memory.fallbackGuardBits = fallbackGuardBits ? std::clamp<mpfr_prec_t>(static_cast<mpfr_prec_t>(strtoll(fallbackGuardBits, nullptr, 10)), 0, 4096) : 128;
        request.threading.tileWidth = 16;
        request.threading.tileHeight = 8;
        return request;
    };
    auto render = [&](const char* name, ExpressionDeepRenderRequest request, ExpressionDeepRenderResult& result, double& seconds) {
        const Clock::time_point start = Clock::now();
        const bool okay = formula::renderExpressionDeepFrame(request, result);
        seconds = since(start);
        printf("  %-20s %.3f s fast/fallback=%llu/%llu preflight=%llu/%llu %s Taylor=%llu/%llu preflight iter/ops/folds=%llu/%llu/%llu fast iter/ops/folds=%llu/%llu/%llu total_iter=%llu precision=%lld/%lld specialized=%d bigfixed=%d mpfr pixels/iter/periodic/generic=%llu/%llu/%llu/%llu bf pixels/iter=%llu/%llu\n", name, seconds, (unsigned long long)result.fastPixelCount, (unsigned long long)result.fallbackPixelCount, (unsigned long long)result.preflightFallbackCount, (unsigned long long)result.preflightSampleCount, formula::expressionDeepPreflightDecisionName(result.preflightDecision), (unsigned long long)result.taylorAcceptedJetCount, (unsigned long long)result.taylorAttemptedJetCount, (unsigned long long)result.preflightIterationCount, (unsigned long long)result.preflightOperationCount, (unsigned long long)result.preflightFoldOperationCount, (unsigned long long)result.fastIterationCount, (unsigned long long)result.fastOperationCount, (unsigned long long)result.fastFoldOperationCount, (unsigned long long)result.totalIterations,
               (long long)result.selectedPrecision, (long long)result.fallbackPrecision, result.usedSpecializedPiecewiseMpfr ? 1 : 0, result.usedPiecewiseBigFixed ? 1 : 0, (unsigned long long)result.specializedPiecewiseMpfrPixelCount, (unsigned long long)result.specializedPiecewiseMpfrIterationCount, (unsigned long long)result.specializedPiecewiseMpfrPeriodicPixelCount, (unsigned long long)result.genericMpfrPeriodicPixelCount, (unsigned long long)result.piecewiseBigFixedPixelCount, (unsigned long long)result.piecewiseBigFixedIterationCount);
        printf("    first uncertain preflight:");
        for (size_t bin = 0; bin < result.preflightFirstUncertainHistogram.size(); ++bin)
            if (result.preflightFirstUncertainHistogram[bin]) printf(" b%zu=%llu", bin, (unsigned long long)result.preflightFirstUncertainHistogram[bin]);
        printf(" reasons:");
        for (size_t reason = 0; reason < result.preflightFallbackReasonCounts.size(); ++reason)
            if (result.preflightFallbackReasonCounts[reason]) printf(" %s=%llu", formula::expressionDeepFallbackReasonName(static_cast<ExpressionDeepFallbackReason>(reason)), (unsigned long long)result.preflightFallbackReasonCounts[reason]);
        printf("\n");
        printf("    full fallback reasons:");
        for (size_t reason = 0; reason < result.fallbackReasonCounts.size(); ++reason)
            if (result.fallbackReasonCounts[reason]) printf(" %s=%llu", formula::expressionDeepFallbackReasonName(static_cast<ExpressionDeepFallbackReason>(reason)), (unsigned long long)result.fallbackReasonCounts[reason]);
        printf(" first uncertain:");
        for (size_t bin = 0; bin < result.fallbackFirstUncertainHistogram.size(); ++bin)
            if (result.fallbackFirstUncertainHistogram[bin]) printf(" b%zu=%llu", bin, (unsigned long long)result.fallbackFirstUncertainHistogram[bin]);
        printf("\n");
        if (!okay) printf("    status=%s error=%s\n", formula::expressionDeepRenderStatusName(result.status), result.error.c_str());
        return okay;
    };

    std::vector<float> production;
    ExpressionDeepRenderRequest productionRequest = makeRequest(production);
    ExpressionDeepRenderResult productionResult;
    double productionSeconds = 0.0;
    bool okay = render("production", productionRequest, productionResult, productionSeconds);

    std::vector<float> specialized;
    ExpressionDeepRenderRequest specializedRequest = makeRequest(specialized);
    specializedRequest.preflight.enable = false;
    specializedRequest.taylor.enableTaylor = false;
    specializedRequest.forceMpfrFallbackForVerification = true;
    specializedRequest.disablePiecewiseBigFixedForVerification = true;
    ExpressionDeepRenderResult specializedResult;
    double specializedSeconds = 0.0;
    okay = render("all-MPFR specialized", specializedRequest, specializedResult, specializedSeconds) && okay;

    std::vector<float> generic;
    ExpressionDeepRenderRequest genericRequest = makeRequest(generic);
    genericRequest.preflight.enable = false;
    genericRequest.taylor.enableTaylor = false;
    genericRequest.forceMpfrFallbackForVerification = true;
    genericRequest.disableSpecializedPiecewiseMpfrForVerification = true;
    ExpressionDeepRenderResult genericResult;
    double genericSeconds = 0.0;
    okay = render("all-MPFR generic", genericRequest, genericResult, genericSeconds) && okay;
    size_t centeredReference = pixelCount;
    if (getenv("MANDEL_CUSTOM_SLOW_CENTERED4")) {
        std::vector<std::pair<long long, size_t>> interiorCandidates;
        long long nearestDistance = std::numeric_limits<long long>::max();
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                const size_t index = static_cast<size_t>(y) * width + x;
                if (generic[index] >= 0.0f) continue;
                const long long dx = 2LL * x - (width - 1LL);
                const long long dy = 2LL * y - (height - 1LL);
                const long long distance = dx * dx + dy * dy;
                interiorCandidates.emplace_back(distance, index);
                if (distance < nearestDistance) {
                    centeredReference = index;
                    nearestDistance = distance;
                }
            }
        }
        std::sort(interiorCandidates.begin(), interiorCandidates.end());
        printf("  interior candidates:");
        for (size_t candidate = 0; candidate < std::min<size_t>(interiorCandidates.size(), 12); ++candidate) printf(" (%zu,%zu)", interiorCandidates[candidate].second % static_cast<size_t>(width), interiorCandidates[candidate].second / static_cast<size_t>(width));
        printf("\n");
        if (centeredReference != pixelCount) printf("  nearest interior pixel=(%zu,%zu)\n", centeredReference % static_cast<size_t>(width), centeredReference / static_cast<size_t>(width));
    }

    if (getenv("MANDEL_CUSTOM_SLOW_GUARD") || getenv("MANDEL_CUSTOM_SLOW_FALLBACK_GUARD")) {
        std::vector<float> highPrecision(pixelCount, formula::ExpressionDeepEmptyPixel);
        ExpressionDeepRenderRequest highRequest = genericRequest;
        highRequest.output = highPrecision.data();
        highRequest.outputCount = highPrecision.size();
        highRequest.precision.requestedBits = 319;
        highRequest.precision.minimumBits = 53;
        highRequest.precision.guardBits = 0;
        highRequest.memory.fallbackGuardBits = 128;
        ExpressionDeepRenderResult highResult;
        double highSeconds = 0.0;
        okay = render("all-MPFR high GT", highRequest, highResult, highSeconds) && okay;
        size_t classMismatch = 0;
        double maxDifference = 0.0;
        double sumDifference = 0.0;
        std::vector<double> differences;
        differences.reserve(pixelCount);
        for (size_t index = 0; index < pixelCount; ++index) {
            const bool lowInterior = generic[index] < 0.0f;
            const bool highInterior = highPrecision[index] < 0.0f;
            if (lowInterior != highInterior) {
                ++classMismatch;
            } else if (!lowInterior) {
                const double difference = std::fabs(static_cast<double>(generic[index]) - static_cast<double>(highPrecision[index]));
                maxDifference = std::max(maxDifference, difference);
                sumDifference += difference;
                differences.push_back(difference);
            }
        }
        const double meanDifference = differences.empty() ? 0.0 : sumDifference / differences.size();
        double p99Difference = 0.0;
        if (!differences.empty()) {
            const size_t percentileIndex = static_cast<size_t>(std::ceil(0.99 * differences.size())) - 1;
            std::nth_element(differences.begin(), differences.begin() + percentileIndex, differences.end());
            p99Difference = differences[percentileIndex];
        }
        printf("  low-vs-high GT class=%zu max/mean/p99=%.3f/%.3f/%.3f low/high=%.3f/%.3f s\n", classMismatch, maxDifference, meanDifference, p99Difference, genericSeconds, highSeconds);
    }

    if (getenv("MANDEL_CUSTOM_SLOW_CENTERED4")) {
        std::vector<float> centered;
        double centeredSeconds = 0.0;
        size_t failedPixels = 0;
        if (centeredReference == pixelCount) {
            printf("  centered VM4 SKIP: frame has no bounded reference candidate\n");
        } else if (!renderExpressionCentered4(canonical, runtime, fixed, centerReal, centerImaginary, scaleText, width, height, static_cast<int>(centeredReference % static_cast<size_t>(width)), static_cast<int>(centeredReference / static_cast<size_t>(width)), maxIterations, bailout, centered, centeredSeconds, failedPixels)) {
            printf("  centered VM4 render failed\n");
            okay = false;
        } else {
            size_t classMismatch = 0;
            double maxDifference = 0.0;
            double sumDifference = 0.0;
            std::vector<double> differences;
            differences.reserve(pixelCount);
            int diagnostics = 0;
            for (size_t index = 0; index < pixelCount; ++index) {
                const bool candidateInterior = centered[index] < 0.0f;
                const bool oracleInterior = generic[index] < 0.0f;
                if (candidateInterior != oracleInterior) {
                    ++classMismatch;
                    if (diagnostics++ < 16) printf("    centered mismatch (%zu,%zu) candidate=%.9g GT=%.9g class\n", index % static_cast<size_t>(width), index / static_cast<size_t>(width), centered[index], generic[index]);
                } else if (!candidateInterior) {
                    const double difference = std::fabs(static_cast<double>(centered[index]) - static_cast<double>(generic[index]));
                    maxDifference = std::max(maxDifference, difference);
                    sumDifference += difference;
                    differences.push_back(difference);
                    if (difference != 0.0 && diagnostics++ < 16) printf("    centered mismatch (%zu,%zu) candidate=%.9g GT=%.9g diff=%.9g\n", index % static_cast<size_t>(width), index / static_cast<size_t>(width), centered[index], generic[index], difference);
                }
            }
            const double meanDifference = differences.empty() ? 0.0 : sumDifference / differences.size();
            double p99Difference = 0.0;
            if (!differences.empty()) {
                const size_t percentileIndex = static_cast<size_t>(std::ceil(0.99 * differences.size())) - 1;
                std::nth_element(differences.begin(), differences.begin() + percentileIndex, differences.end());
                p99Difference = differences[percentileIndex];
            }
            printf("  centered VM4 GT class=%zu max/mean/p99=%.3f/%.3f/%.3f failed=%zu time=%.3f s\n", classMismatch, maxDifference, meanDifference, p99Difference, failedPixels, centeredSeconds);
            for (int threshold : {512, 768, 900, 1024}) {
                uint64_t fallbackIterations = 0;
                size_t fallbackCount = 0;
                size_t selectiveClassMismatch = 0;
                double selectiveMaxDifference = 0.0;
                double selectiveSumDifference = 0.0;
                std::vector<double> selectiveDifferences;
                for (size_t index = 0; index < pixelCount; ++index) {
                    const float candidate = centered[index];
                    const bool fallback = candidate < 0.0f || candidate >= static_cast<float>(threshold);
                    const int oracleIterations = generic[index] < 0.0f ? maxIterations : static_cast<int>(generic[index]);
                    if (fallback) {
                        ++fallbackCount;
                        fallbackIterations += static_cast<uint64_t>(oracleIterations);
                    } else {
                        const bool candidateInterior = candidate < 0.0f;
                        const bool oracleInterior = generic[index] < 0.0f;
                        if (candidateInterior != oracleInterior) {
                            ++selectiveClassMismatch;
                        } else if (!candidateInterior) {
                            const double difference = std::fabs(static_cast<double>(candidate) - generic[index]);
                            selectiveMaxDifference = std::max(selectiveMaxDifference, difference);
                            selectiveSumDifference += difference;
                            selectiveDifferences.push_back(difference);
                        }
                    }
                }
                const double selectiveMeanDifference = selectiveDifferences.empty() ? 0.0 : selectiveSumDifference / selectiveDifferences.size();
                double selectiveP99Difference = 0.0;
                if (!selectiveDifferences.empty()) {
                    const size_t percentile = static_cast<size_t>(std::ceil(0.99 * selectiveDifferences.size())) - 1;
                    std::nth_element(selectiveDifferences.begin(), selectiveDifferences.begin() + percentile, selectiveDifferences.end());
                    selectiveP99Difference = selectiveDifferences[percentile];
                }
                printf("    threshold=%d fallback=%zu/%zu MPFR-iterations=%llu/%llu (%.1f%%) quality=%zu/%.0f/%.3f/%.0f\n", threshold, fallbackCount, pixelCount, (unsigned long long)fallbackIterations, (unsigned long long)genericResult.totalIterations, genericResult.totalIterations ? 100.0 * fallbackIterations / genericResult.totalIterations : 0.0, selectiveClassMismatch, selectiveMaxDifference, selectiveMeanDifference, selectiveP99Difference);
            }
        }
    }

    if (getenv("MANDEL_CUSTOM_SLOW_DD") && burningShipResearchEligible) {
        std::vector<float> doubleDouble;
        double doubleDoubleSeconds = 0.0;
        if (!renderBurningShipDoubleDouble(centerReal, centerImaginary, scaleText, width, height, maxIterations, bailout, doubleDouble, doubleDoubleSeconds)) {
            printf("  double-double render failed\n");
            okay = false;
        } else {
            size_t classMismatch = 0;
            double maxDifference = 0.0;
            double sumDifference = 0.0;
            std::vector<double> differences;
            differences.reserve(pixelCount);
            for (size_t index = 0; index < pixelCount; ++index) {
                const bool candidateInterior = doubleDouble[index] < 0.0f;
                const bool oracleInterior = generic[index] < 0.0f;
                if (candidateInterior != oracleInterior) {
                    ++classMismatch;
                } else if (!candidateInterior) {
                    const double difference = std::fabs(static_cast<double>(doubleDouble[index]) - static_cast<double>(generic[index]));
                    maxDifference = std::max(maxDifference, difference);
                    sumDifference += difference;
                    differences.push_back(difference);
                }
            }
            const double meanDifference = differences.empty() ? 0.0 : sumDifference / differences.size();
            double p99Difference = 0.0;
            if (!differences.empty()) {
                const size_t percentileIndex = static_cast<size_t>(std::ceil(0.99 * differences.size())) - 1;
                std::nth_element(differences.begin(), differences.begin() + percentileIndex, differences.end());
                p99Difference = differences[percentileIndex];
            }
            printf("  double-double GT class=%zu max/mean/p99=%.3f/%.3f/%.3f time=%.3f s\n", classMismatch, maxDifference, meanDifference, p99Difference, doubleDoubleSeconds);
        }
    }

    if (const char* mpfBits = getenv("MANDEL_CUSTOM_SLOW_MPF_BITS"); mpfBits && burningShipResearchEligible) {
        const mp_bitcnt_t precision = static_cast<mp_bitcnt_t>(std::clamp<unsigned long long>(strtoull(mpfBits, nullptr, 10), 64, 4096));
        std::vector<float> mpfOutput;
        double mpfSeconds = 0.0;
        if (!renderBurningShipMpf(centerReal, centerImaginary, scaleText, width, height, maxIterations, bailout, precision, mpfOutput, mpfSeconds)) {
            printf("  fixed-mpf render failed\n");
            okay = false;
        } else {
            size_t classMismatch = 0;
            double maxDifference = 0.0;
            double sumDifference = 0.0;
            std::vector<double> differences;
            differences.reserve(pixelCount);
            for (size_t index = 0; index < pixelCount; ++index) {
                const bool candidateInterior = mpfOutput[index] < 0.0f;
                const bool oracleInterior = generic[index] < 0.0f;
                if (candidateInterior != oracleInterior) {
                    ++classMismatch;
                } else if (!candidateInterior) {
                    const double difference = std::fabs(static_cast<double>(mpfOutput[index]) - static_cast<double>(generic[index]));
                    maxDifference = std::max(maxDifference, difference);
                    sumDifference += difference;
                    differences.push_back(difference);
                }
            }
            const double meanDifference = differences.empty() ? 0.0 : sumDifference / differences.size();
            double p99Difference = 0.0;
            if (!differences.empty()) {
                const size_t percentileIndex = static_cast<size_t>(std::ceil(0.99 * differences.size())) - 1;
                std::nth_element(differences.begin(), differences.begin() + percentileIndex, differences.end());
                p99Difference = differences[percentileIndex];
            }
            printf("  fixed-mpf %llu-bit GT class=%zu max/mean/p99=%.3f/%.3f/%.3f time=%.3f s\n", static_cast<unsigned long long>(precision), classMismatch, maxDifference, meanDifference, p99Difference, mpfSeconds);
        }
    }

    if (const char* limbText = getenv("MANDEL_CUSTOM_SLOW_BF_LIMBS"); limbText && burningShipResearchEligible) {
        const int limbs = std::clamp(atoi(limbText), 3, 32);
        std::vector<float> fixedOutput;
        double fixedSeconds = 0.0;
        if (!renderBurningShipBigFixed(centerReal, centerImaginary, scaleText, width, height, maxIterations, bailout, limbs, fixedOutput, fixedSeconds)) {
            printf("  BigFixed render failed\n");
            okay = false;
        } else {
            size_t classMismatch = 0;
            double maxDifference = 0.0;
            double sumDifference = 0.0;
            std::vector<double> differences;
            differences.reserve(pixelCount);
            for (size_t index = 0; index < pixelCount; ++index) {
                const bool candidateInterior = fixedOutput[index] < 0.0f;
                const bool oracleInterior = generic[index] < 0.0f;
                if (candidateInterior != oracleInterior) {
                    ++classMismatch;
                } else if (!candidateInterior) {
                    const double difference = std::fabs(static_cast<double>(fixedOutput[index]) - static_cast<double>(generic[index]));
                    maxDifference = std::max(maxDifference, difference);
                    sumDifference += difference;
                    differences.push_back(difference);
                }
            }
            const double meanDifference = differences.empty() ? 0.0 : sumDifference / differences.size();
            double p99Difference = 0.0;
            if (!differences.empty()) {
                const size_t percentileIndex = static_cast<size_t>(std::ceil(0.99 * differences.size())) - 1;
                std::nth_element(differences.begin(), differences.begin() + percentileIndex, differences.end());
                p99Difference = differences[percentileIndex];
            }
            printf("  BigFixed L=%d GT class=%zu max/mean/p99=%.3f/%.3f/%.3f time=%.3f s\n", limbs, classMismatch, maxDifference, meanDifference, p99Difference, fixedSeconds);
        }
    }

    if (production != generic || specialized != generic) {
        size_t productionMismatch = 0;
        size_t specializedMismatch = 0;
        for (size_t index = 0; index < pixelCount; ++index) {
            productionMismatch += production[index] != generic[index];
            specializedMismatch += specialized[index] != generic[index];
        }
        printf("  MPFR parity mismatch production/specialized=%zu/%zu\n", productionMismatch, specializedMismatch);
        okay = false;
    }

    if (pixelCount <= 50000) {
        std::vector<float> current;
        ExpressionDeepRenderRequest currentRequest = makeRequest(current);
        currentRequest.preflight.enable = false;
        ExpressionDeepRenderResult currentResult;
        double currentSeconds = 0.0;
        okay = render("certified current", currentRequest, currentResult, currentSeconds) && okay;

        std::vector<float> perStep;
        ExpressionDeepRenderRequest perStepRequest = makeRequest(perStep);
        perStepRequest.preflight.enable = false;
        perStepRequest.taylor.enableTaylor = false;
        ExpressionDeepRenderResult perStepResult;
        double perStepSeconds = 0.0;
        okay = render("scaled no Taylor", perStepRequest, perStepResult, perStepSeconds) && okay;
        if (current != generic || perStep != generic) {
            printf("  certified variant output mismatch\n");
            okay = false;
        }
    }

    mpf_t centerRe, centerIm, scale;
    mpf_init2(centerRe, 256);
    mpf_init2(centerIm, 256);
    mpf_init2(scale, 256);
    const bool parsed = mpf_set_str(centerRe, centerReal, 10) == 0 && mpf_set_str(centerIm, centerImaginary, 10) == 0 && mpf_set_str(scale, scaleText, 10) == 0;
    std::vector<float> direct(pixelCount, formula::ExpressionDeepEmptyPixel);
    Mandel directRenderer(width, height, maxIterations, 1, direct.data());
    ExpressionOrbitPlan orbitPlan;
    ExpressionJit4 jit;
    const bool havePlan = orbitPlan.build(runtime, &error) && orbitPlan.profitable();
    const bool haveJit = VERIFY_JIT && (havePlan ? jit.compile(orbitPlan) : jit.compile(runtime));
    const Clock::time_point directStart = Clock::now();
    const bool directOkay = parsed && directRenderer.ComputeExpression(centerRe, centerIm, scale, runtime, fixed, FormulaParameter::C, maxIterations, bailout, formula::ExpressionColoring::Raw, haveJit ? &jit : nullptr, havePlan ? &orbitPlan : nullptr);
    const double directSeconds = since(directStart);
    mpf_clears(centerRe, centerIm, scale, (mpf_ptr)0);
    size_t directMismatch = 0;
    size_t directClassMismatch = 0;
    double directMaxDifference = 0.0;
    double directSumDifference = 0.0;
    std::vector<double> directDifferences;
    directDifferences.reserve(pixelCount);
    for (size_t index = 0; directOkay && index < pixelCount; ++index) {
        if (direct[index] != generic[index]) ++directMismatch;
        const bool directInterior = direct[index] < 0.0f;
        const bool genericInterior = generic[index] < 0.0f;
        if (directInterior != genericInterior) {
            ++directClassMismatch;
        } else if (!directInterior) {
            const double difference = std::fabs(static_cast<double>(direct[index]) - static_cast<double>(generic[index]));
            directMaxDifference = std::max(directMaxDifference, difference);
            directSumDifference += difference;
            directDifferences.push_back(difference);
        }
    }
    const double directMeanDifference = directDifferences.empty() ? 0.0 : directSumDifference / directDifferences.size();
    double directP99Difference = 0.0;
    if (!directDifferences.empty()) {
        const size_t percentileIndex = static_cast<size_t>(std::ceil(0.99 * directDifferences.size())) - 1;
        std::nth_element(directDifferences.begin(), directDifferences.begin() + percentileIndex, directDifferences.end());
        directP99Difference = directDifferences[percentileIndex];
    }
    const double zoom = std::strtod(scaleText, nullptr);
    const double centerReDouble = std::strtod(centerReal, nullptr);
    const double centerImDouble = std::strtod(centerImaginary, nullptr);
    const double halfWidth = 2.0 / zoom;
    const double halfHeight = halfWidth * static_cast<double>(height) / static_cast<double>(width);
    const double dx = 2.0 * halfWidth / static_cast<double>(width - 1);
    const double dy = 2.0 * halfHeight / static_cast<double>(height - 1);
    auto ulp = [](double value) {
        const double magnitude = std::fabs(value);
        return std::nextafter(magnitude, std::numeric_limits<double>::infinity()) - magnitude;
    };
    const double xUlp = ulp(std::fabs(centerReDouble) + halfWidth);
    const double yUlp = ulp(std::fabs(centerImDouble) + halfHeight);
    const double coordinateUlpMargin = std::min(dx / xUlp, dy / yUlp);
    printf("  %-20s %.3f s mismatch/class=%zu/%zu max/mean/p99=%.3f/%.3f/%.3f coordinate margin=%.2f ulp; recurrence proof rejected\n", "direct binary64", directSeconds, directMismatch, directClassMismatch, directMaxDifference, directMeanDifference, directP99Difference, coordinateUlpMargin);

    const bool preflightOkay = productionResult.preflightAttempted && productionResult.preflightRejectedFast && productionResult.preflightSampleCount >= 8 && productionResult.preflightFallbackCount == productionResult.preflightSampleCount && productionResult.fastPixelCount == 0 && productionResult.fallbackPixelCount == pixelCount;
    const double speedup = productionSeconds > 0.0 ? genericSeconds / productionSeconds : 0.0;
    printf("=== custom diffabs production slowdown %dx%d\n", width, height);
    printf("  production/generic MPFR %.3f/%.3f s speedup %.2fx direct exact=%d\n", productionSeconds, genericSeconds, speedup, directOkay && directMismatch == 0 ? 1 : 0);
    printf("  => %s\n\n", okay && preflightOkay ? "PASS" : "CHECK (custom slowdown failure)");
    return okay && preflightOkay ? 0 : 1;
}

int main(int argc, char** argv) {
    std::string which = (argc > 1) ? argv[1] : "all";
    const int W = (argc > 2) ? atoi(argv[2]) : 120;
    const int H = (argc > 3) ? atoi(argv[3]) : 80;

    // step=1: brute-force the whole image. step>1: sparse oracle (feasible when
    // the location needs a high mxit and deep precision).
    TestCase shallow{"shallow (whole set)", "-0.5", "0", pow10(0), 5000, 1};
    std::string deep51_scale = "3831277";
    deep51_scale.append(45, '0'); // 3.831277e51
    TestCase deep{"deep51 (3.8e51 stress, maxit 2M)", testcases::deep51_x, testcases::deep51_y, deep51_scale, 2000000, 32};
    TestCase ticktock{"ticktock (1e141, glitch stress)", testcases::ticktock_x, testcases::ticktock_y, pow10(141), 200000, 12};
    TestCase flake{"flake (1e157, glitch stress)", testcases::flake_x, testcases::flake_y, pow10(157), 200000, 12};
    // Exact Misiurewicz parameter c=i. Its critical orbit enters a repelling
    // period-2 cycle; 1e1000 pixel deltas amplify and escape after ~2650
    // iterations, making this a stable, genuinely deep exterior stress case.
    TestCase exterior1000{"exterior (Misiurewicz c=i, 1e1000)", "0", "1", pow10(1000), 10000, 1};
    // A maxit-sensitive period-1329 component boundary near c=i. This exposed a
    // deep-reference bug where even and odd image sizes followed different
    // reference pixels and flipped almost the entire image classification.
    const std::string parityRe = "0.87658692204855777255720127240577294219664089579606052252153341247110316040359567346941066530198264333715350076914731265081285605596159889659482671637665433036047479449476861255305601870749958937987160494800193599269036070384980604718254475114301562016620167318503180246283533354808944647785156562902419144968711947754196044352006391776697598867374433836375101918562145980477898443489741920475118804535699922081074341333337897641607372542163474622122834005970824831688364959496394361968461917806403290834529241567373997027330589016e-500";
    std::string parityIm = "0.";
    parityIm.append(498, '9');
    parityIm += "9737732145267474019035193255007042699493697319030874172320234538271496366171723116782567955351384344445642267031664918871592758856681040534525429832129077075898037974437238709557577566760153377050600768882679447728184707424726714911339143303870647306255931485936680095871239708214775864294277313659941154622160496034912392379929833856520170763139643443634966071512688131637522087272155766755998538950242005247687225131897953568014863667043990248406875087728885765117677234374808520874070003498657884677143809794193877711208425737725547735873602";
    TestCase parity1000{"reference parity (period-1329 boundary, 1e1000)", parityRe, parityIm, pow10(1000), 10000, 16};
    // A 1e876 needle point (floatexp path, high-half multiply range). Deep exterior
    // stress with a different reference orbit than the 1e1000 cases.
    TestCase deep876{"deep876 (1e876 needle, floatexp)", testcases::deep1_x, testcases::deep1_y, pow10(876), 300000, 24};

    // Sub-pixel-nucleus stress (regression for the exterior-reference glitch):
    // the deep1 period-14164 minibrot is ~1e-876 wide, so at 1e737 it is ~139
    // orders of magnitude sub-pixel. The whole frame is therefore genuinely
    // EXTERIOR (mpmath confirms the centre escapes at n=41351 at this precision).
    // A full-image oracle catches both failure modes this exposed: a periodic
    // reference that renders the sub-pixel nucleus as a false-interior "black
    // star", and an escaping reference (periodic off / EDE / SAC) that is never
    // re-referenced and returns garbage smooth values. maxMisPct = 0 => any
    // false-interior pixel fails; step = 1 => every pixel is checked.
    TestCase subpixel{"subpixel (deep1 nucleus 1e-876 @ 1e737, all-exterior)", testcases::deep1_x, testcases::deep1_y, pow10(737), 60000, 1};
    subpixel.maxMisPct = 0.0;
    subpixel.maxDiff = 0.5;
    subpixel.w = 48;
    subpixel.h = 36;

    // Acceptance gate for the real, resolved period-14164 minibrot at its natural
    // depth. The saved full-image oracle was generated at 1100 decimal digits and
    // its 81-interior / 351-exterior classification was stable at 1000 and 1100
    // digits. A render-precision oracle falsely reports the whole frame interior,
    // so this case MUST use the saved high-precision golden rather than ComputeDirect.
    TestCase minibrot875{"minibrot875 (resolved period-14164 body + exterior)", testcases::deep1_x, testcases::deep1_y, pow10(875), 300000, 1};
    minibrot875.maxMisPct = 0.0;
    minibrot875.maxDiff = 1.0;
    minibrot875.w = 24;
    minibrot875.h = 18;
    minibrot875.goldenPath = "tests/golden/minibrot875.f32";

    // Shift the 1e875 frame right by 0.55 screen widths. The period-14164 nucleus
    // sits just outside the left edge, so findNucleus returns 0 and the exact view
    // centre is an EXTERIOR full reference (escapes near iteration 161971), while
    // the minibrot body remains visible along the left edge (27/432 interior).
    // This guards the full-reference/deep-zero path independently of the interior
    // periodic-nucleus reference used by minibrot875.
    TestCase extref875{"extref875 (exterior full reference, visible minibrot)", testcases::deep1_x, testcases::deep1_y, pow10(875), 300000, 1};
    extref875.maxMisPct = 0.0;
    extref875.maxDiff = 1.0;
    extref875.w = 24;
    extref875.h = 18;
    extref875.goldenPath = "tests/golden/minibrot875_extref.f32";
    extref875.centerShiftRe = "2.2e-875";

    // Moderate-zoom high-iteration stress: BLA is ineffective and ~39% of pixels
    // are interior, so SIMD perturbation/interior detection dominates. The current
    // double perturbation has known smooth-value drift on a small chaotic tail; lock
    // classification exactly and bound max/mean/p99 so performance work cannot worsen it.
    TestCase slowpoint{"slowpoint (7.826796e8, interior-detector stress)", "-1.1758621450236620370", "-0.2447677973532398022", "782679600", 500000, 1};
    slowpoint.maxMisPct = 0.0;
    slowpoint.maxDiff = 90000.0;
    slowpoint.maxMeanDiff = 300.0;
    slowpoint.maxP99Diff = 3000.0;
    slowpoint.w = 48;
    slowpoint.h = 32;
    slowpoint.goldenPath = "tests/golden/slowpoint_7e8.f32";

    std::string point31Scale = "7354177";
    point31Scale.append(25, '0');
    TestCase point31{"point31 (7.354177e31, SA+BLA+SIMD stress)", "-0.749139567333446841955467474699747367338762518832278501811", "0.040823298514634751035521346975478853963578400940553676068", point31Scale, 500000, 1};
    point31.maxMisPct = 0.0;
    point31.maxDiff = 12000.0;
    point31.maxMeanDiff = 200.0;
    point31.maxP99Diff = 4000.0;
    point31.w = 48;
    point31.h = 32;
    point31.goldenPath = "tests/golden/point_7e31.f32";

    // Optional reference-footprint sweep: override mxit for all cases.
    if (const char* e = getenv("MANDEL_MXIT")) {
        int m = atoi(e);
        shallow.mxit = deep.mxit = ticktock.mxit = flake.mxit = exterior1000.mxit = parity1000.mxit = m;
    }
    if (const char* e = getenv("MANDEL_ORACLE_STEP")) {
        int step = std::max(1, atoi(e));
        shallow.step = deep.step = ticktock.step = flake.step = exterior1000.step = parity1000.step = step;
        deep876.step = subpixel.step = minibrot875.step = extref875.step = slowpoint.step = point31.step = step;
    }
    // Optional scale override (decimal exponent) to find interior-containing frames.
    if (const char* e = getenv("MANDEL_SCALE_EXP")) {
        std::string sc = pow10(atoi(e));
        deep.scale = ticktock.scale = flake.scale = exterior1000.scale = parity1000.scale = sc;
    }

    int rc = 0;
    if (which == "all" || which == "shallow") rc |= runCase(shallow, W, H);
    if (which == "all" || which == "deep") rc |= runCase(deep, W, H);
    if (which == "all" || which == "ticktock") rc |= runCase(ticktock, W, H);
    if (which == "all" || which == "flake") rc |= runCase(flake, W, H);
    if (which == "exterior1000") rc |= runCase(exterior1000, W, H);
    if (which == "parity1000") rc |= runCase(parity1000, W, H);
    if (which == "deep876") rc |= runCase(deep876, W, H);
    if (which == "subpixel") rc |= runCase(subpixel, W, H);
    if (which == "minibrot875") rc |= runCase(minibrot875, W, H);
    if (which == "extref875") rc |= runCase(extref875, W, H);
    if (which == "slowpoint") rc |= runCase(slowpoint, W, H);
    if (which == "point31") rc |= runCase(point31, W, H);
    if (which == "gui875") rc |= runAdaptiveGuiCase();
    if (which == "julia") rc |= runJuliaCase(false);
    if (which == "julia-ede") rc |= runJuliaCase(true);
    if (which == "julia-dendrite") rc |= runJuliaCase(false, true);
    if (which == "julia-critical") rc |= runJuliaCase(false, true, true);
    if (which == "expression") {
        rc |= runExpressionCoreCase();
        rc |= runExpressionCenteredCase();
        rc |= runExpressionColoringCase();
    }
    if (which == "expression-centered") rc |= runExpressionCenteredCase();
    if (which == "expression-reference") rc |= runExpressionReferenceCase();
    if (which == "expression-scaled") rc |= runExpressionScaledCase();
    if (which == "expression-taylor") rc |= runExpressionTaylorCase();
    if (which == "expression-deep-render") rc |= runExpressionDeepRenderCase();
    if (which == "expression-coloring") rc |= runExpressionColoringCase();
    if (which == "expression-orbit") rc |= runExpressionOrbitCase();
    if (which == "expression-oracle" || which == "oracle") rc |= runExpressionOracleCase();
    if (which == "expression-suite" || which == "suite") rc |= runFormulaRegressionSuite();
    if (which == "custom-deep") rc |= runCustomDeepZoomCase();
    if (which == "expression-residual") rc |= runExpressionResidualSuite();
    if (which == "formula-bench") rc |= runGenericFormulaProfile();
    if (which == "formula-periodic") rc |= runGenericPeriodicityCase();
    if (which == "multibrot") rc |= runMultibrotCase();
    if (which == "backend") rc |= runBackendCase();
    if (which == "generic-deep") rc |= runGenericDeepBackendCase();
    if (which == "custom-slowdown") rc |= runCustomSlowdownCase(W, H);
    if (which == "gpu") rc |= runGpuBenchmarkCase(argc > 2 ? W : 1920, argc > 3 ? H : 1080);
    return rc;
}
