// Headless verification + benchmark harness for the Mandelbrot renderer.
//
// For a given view it computes the escape-time buffer two ways:
//   (1) the perturbation / series-approximation / rebasing engine (Compute), and
//   (2) a brute-force full-precision reference (ComputeDirect / accuratePointCompute).
// It then reports how well they agree, plus timings and a checksum for
// regression tracking. No libpng / OpenGL / GLUT dependencies.
//
// Usage: verify [case]
//   case = shallow | deep | all   (default: all)

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
    const char* cx;
    const char* cy;
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
    // Precision: enough bits for the scale magnitude plus guard.
    int precision = static_cast<int>(tc.scale.size() * log(10) / log(2)) + 40;
    mpf_set_default_prec(precision);

    mpf_t cre, cim, scale;
    mpf_init_set_str(cre, tc.cx, 10);
    mpf_init_set_str(cim, tc.cy, 10);
    mpf_init_set_str(scale, tc.scale.c_str(), 10);

    float* itp = new float[W * H];
    float* itd = new float[W * H];
    for (int i = 0; i < W * H; ++i) itp[i] = itd[i] = EMPTYPIXEL;

    Mandel mandel(W, H, tc.mxit, 1, itp);
    mandel.setPrecision(precision);

    printf("=== %s  (%dx%d, scale=1e%zu, mxit=%d, prec=%d bits)\n",
           tc.name, W, H, tc.scale.size() - 1, tc.mxit, precision);

    auto t0 = Clock::now();
    mandel.Compute(cre, cim, scale, tc.mxit);
    double t_pert = since(t0);

    bool run_oracle = (getenv("MANDEL_NOORACLE") == nullptr);
    t0 = Clock::now();
    if (run_oracle) mandel.ComputeDirect(tc.mxit, itd, tc.step);
    double t_direct = since(t0);

    // Compare only pixels the (possibly sparse) oracle actually computed.
    long n_int_both = 0, n_ext_both = 0, n_class_mismatch = 0, n_sampled = 0;
    double max_diff = 0, sum_diff = 0;
    int worst_i = -1, worst_j = -1;
    for (int i = 0; i < H; ++i) {
        for (int j = 0; j < W; ++j) {
            float d = itd[i * W + j];
            if (d == EMPTYPIXEL) continue;      // not sampled by the oracle
            float p = itp[i * W + j];
            ++n_sampled;
            bool p_int = (p < 0);
            bool d_int = (d < 0);
            if (p_int != d_int) {
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
    TestCase deep{ "deep (1e30 into deep1)", testcases::deep1_x, testcases::deep1_y, pow10(30), 30000, 1 };
    TestCase ticktock{ "ticktock (1e141, glitch stress)", testcases::ticktock_x, testcases::ticktock_y, pow10(141), 200000, 12 };
    TestCase flake{ "flake (1e157, glitch stress)", testcases::flake_x, testcases::flake_y, pow10(157), 200000, 12 };

    // Optional reference-footprint sweep: override mxit for all cases.
    if (const char* e = getenv("MANDEL_MXIT")) {
        int m = atoi(e);
        shallow.mxit = deep.mxit = ticktock.mxit = flake.mxit = m;
    }
    // Optional scale override (decimal exponent) to find interior-containing frames.
    if (const char* e = getenv("MANDEL_SCALE_EXP")) {
        std::string sc = pow10(atoi(e));
        deep.scale = ticktock.scale = flake.scale = sc;
    }

    int rc = 0;
    if (which == "all" || which == "shallow")  rc |= runCase(shallow, W, H);
    if (which == "all" || which == "deep")     rc |= runCase(deep, W, H);
    if (which == "all" || which == "ticktock") rc |= runCase(ticktock, W, H);
    if (which == "all" || which == "flake")    rc |= runCase(flake, W, H);
    return rc;
}
