// Headless verification + benchmark harness for the Mandelbrot renderer.
//
// For a given view it computes the escape-time buffer two ways:
//   (1) the perturbation / series-approximation / rebasing engine (Compute), and
//   (2) a brute-force full-precision reference (ComputeDirect / accuratePointCompute).
// It then reports how well they agree, plus timings and a checksum for
// regression tracking. No libpng / OpenGL / GLUT dependencies.
//
// Usage: verify [case]
//   case = shallow | deep | ticktock | flake | exterior1000 | parity1000 | all
//   The 1e1000 cases are excluded from "all" because their 3400-bit GMP oracles
//   are intentionally much more expensive than the regular regression set.

#include <gmp.h>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <string>
#include <chrono>

#include "mandel_perturbation.h"
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
};

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

static int runCase(const TestCase& tc, int W, int H) {
    // Precision follows the visible scale; 64 guard bits retain center digits
    // beyond the last visible decimal without making shallow tests use the full
    // textual coordinate length.
    int precision = static_cast<int>(tc.scale.size() * log(10) / log(2)) + 64;
    mpf_set_default_prec(precision);

    mpf_t cre, cim, scale;
    mpf_init_set_str(cre, tc.cx.c_str(), 10);
    mpf_init_set_str(cim, tc.cy.c_str(), 10);
    mpf_init_set_str(scale, tc.scale.c_str(), 10);

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

    bool run_oracle = (getenv("MANDEL_NOORACLE") == nullptr);
    t0 = Clock::now();
    if (run_oracle) mandel.ComputeDirect(tc.mxit, itd, tc.step, c_method);
    double t_direct = since(t0);

    // Compare only pixels the (possibly sparse) oracle actually computed.
    long n_int_both = 0, n_ext_both = 0, n_class_mismatch = 0, n_sampled = 0;
    double max_diff = 0, sum_diff = 0;
    int worst_i = -1, worst_j = -1;
    int first_class_i = -1, first_class_j = -1;
    for (int i = 0; i < H; ++i) {
        for (int j = 0; j < W; ++j) {
            float d = itd[i * W + j];
            if (d == EMPTYPIXEL) continue;      // not sampled by the oracle
            float p = itp[i * W + j];
            ++n_sampled;
            bool p_int = (p < 0);
            bool d_int = (d < 0);
            if (p_int != d_int) {
                if (first_class_i < 0) { first_class_i = i; first_class_j = j; }
                ++n_class_mismatch;
                continue;
            }
            if (p_int) { ++n_int_both; continue; }
            ++n_ext_both;
            double diff = std::fabs((double)p - (double)d);
            sum_diff += diff;
            if (diff > max_diff) { max_diff = diff; worst_i = i; worst_j = j; }
        }
    }

    double mean_diff = n_ext_both ? sum_diff / n_ext_both : 0.0;
    // Interior cost estimate over the FULL image: interior pixels (iter < 0)
    // each ran ~mxit iterations; exterior escape-iteration ~ the smooth value.
    long n_interior = 0; double ext_iters = 0;
    for (int q = 0; q < W * H; ++q) {
        if (itp[q] < 0) ++n_interior;
        else ext_iters += itp[q];
    }
    double int_iters = (double)n_interior * tc.mxit;
    double frac_int = (int_iters + ext_iters) > 0 ? int_iters / (int_iters + ext_iters) : 0.0;
    printf("  perturbation : %8.3f s   (full %dx%d image)\n", t_pert, W, H);
    printf("  interior px  : %ld / %d (%.1f%%)   est. iters interior/exterior = %.3g / %.3g  (interior = %.1f%% of all iters)\n",
           n_interior, W * H, 100.0 * n_interior / (W * H), int_iters, ext_iters, 100.0 * frac_int);
    printf("  direct (ref) : %8.3f s   (%ld pixels sampled, step=%d)\n",
           t_direct, n_sampled, tc.step);
    printf("  interior both: %ld   exterior both: %ld   class mismatch: %ld (%.3f%% of sampled)\n",
           n_int_both, n_ext_both, n_class_mismatch,
           n_sampled ? 100.0 * n_class_mismatch / n_sampled : 0.0);
    if (first_class_i >= 0) {
        int q = first_class_i * W + first_class_j;
        printf("  first class mismatch @ %d,%d: pert=%.4f ref=%.4f\n",
               first_class_i, first_class_j, itp[q], itd[q]);
    }
    printf("  escape-time diff vs reference:  max %.6g  mean %.6g", max_diff, mean_diff);
    if (worst_i >= 0)
        printf("   (worst @ %d,%d: pert=%.4f ref=%.4f)",
               worst_i, worst_j, itp[worst_i * W + worst_j], itd[worst_i * W + worst_j]);
    printf("\n  checksum(pert)=0x%08x\n", checksum(itp, W * H));

    // Boundary pixels can legitimately disagree, but the bulk must match.
    double mismatch_pct = n_sampled ? 100.0 * n_class_mismatch / n_sampled : 0.0;
    bool ok = (max_diff < 0.5) && (mismatch_pct < 1.0);
    printf("  => %s\n\n", ok ? "PASS" : "CHECK (see mismatches above)");

    delete[] itp;
    delete[] itd;
    mpf_clear(cre);
    mpf_clear(cim);
    mpf_clear(scale);
    return ok ? 0 : 1;
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

    // Optional reference-footprint sweep: override mxit for all cases.
    if (const char* e = getenv("MANDEL_MXIT")) {
        int m = atoi(e);
        shallow.mxit = deep.mxit = ticktock.mxit = flake.mxit = exterior1000.mxit = parity1000.mxit = m;
    }
    if (const char* e = getenv("MANDEL_ORACLE_STEP")) {
        int step = std::max(1, atoi(e));
        shallow.step = deep.step = ticktock.step = flake.step = exterior1000.step = parity1000.step = step;
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
    return rc;
}
