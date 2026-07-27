// Minimal COMPUTE-ONLY timing harness (no coloring / no BMP overhead) so speed
// comparisons reflect the engine only. Times Mandel::Compute for one location.
//
// Usage: ss_bench <cx> <cy> <scale> <mxit> <W> <H> <sub> [ede] [trials] [sac]
//   sub=1 -> base 1spp pass; sub>1 -> base + adaptive supersample (SUPER_SAMPLING).
//   ede=1 adds EXTERIOR_DIST_EST. Prints the best (min) Compute time in seconds.
//   MANDEL_SIMD / MANDEL_SS_K / MANDEL_COARSE env vars apply as usual.

#include <gmp.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <chrono>
#include "mandel_perturbation.h"

using Clock = std::chrono::high_resolution_clock;

int main(int argc, char** argv) {
    if (argc < 8) { fprintf(stderr, "usage: ss_bench cx cy scale mxit W H sub [ede] [trials] [sac]\n"); return 1; }
    const char* cx = argv[1]; const char* cy = argv[2]; const char* scaleStr = argv[3];
    int mxit = atoi(argv[4]); int W = atoi(argv[5]); int H = atoi(argv[6]); int sub = atoi(argv[7]);
    int ede = argc > 8 ? atoi(argv[8]) : 0;
    int trials = argc > 9 ? atoi(argv[9]) : 3;
    int sac = argc > 10 ? atoi(argv[10]) : 0;
    if (sub < 1) sub = 1;

    double sd = atof(scaleStr);
    int scaleExp = (sd > 1.0) ? (int)(log10(sd) + 0.5) : 0;
    int precision = (int)((scaleExp + 20) * log(10) / log(2)) + 64;
    mpf_set_default_prec(precision);

    int cmethod = 0;
    if (sub > 1) cmethod |= ColoringMethod::SUPER_SAMPLING;
    if (ede)     cmethod |= ColoringMethod::EXTERIOR_DIST_EST;
    if (sac)     cmethod |= ColoringMethod::STRIPE_AVERAGE;

    mpf_t mcx, mcy, msc;
    mpf_init_set_str(mcx, cx, 10); mpf_init_set_str(mcy, cy, 10); mpf_init_set_str(msc, scaleStr, 10);

    size_t N = (size_t)W * H * sub * sub;
    float* iter = new float[N];

    // Correctness self-check: compute the SS buffer with BLA off (trusted, no
    // stale-table risk) vs on, and report interior/exterior class mismatches and
    // exterior value diffs. BLA on must match BLA off closely.
    if (getenv("MANDEL_SSCHECK")) {
        float* a = new float[N];
        _putenv_s("MANDEL_BLA", "0");
        { for (size_t i = 0; i < N; ++i) iter[i] = EMPTYPIXEL;
          Mandel m(W, H, mxit, sub, iter); m.setPrecision(precision);
          m.Compute(mcx, mcy, msc, mxit, cmethod); }
        for (size_t i = 0; i < N; ++i) a[i] = iter[i];
        _putenv_s("MANDEL_BLA", "1");
        { for (size_t i = 0; i < N; ++i) iter[i] = EMPTYPIXEL;
          Mandel m(W, H, mxit, sub, iter); m.setPrecision(precision);
          m.Compute(mcx, mcy, msc, mxit, cmethod); }
        long cls = 0, extboth = 0, sampled = 0; double maxd = 0, sumd = 0;
        for (size_t i = 0; i < N; ++i) {
            if (a[i] == EMPTYPIXEL || iter[i] == EMPTYPIXEL) continue;
            ++sampled;
            bool ai = a[i] < 0, bi = iter[i] < 0;
            if (ai != bi) { ++cls; continue; }
            if (!ai) { double d = fabs((double)a[i] - iter[i]); sumd += d; if (d > maxd) maxd = d; ++extboth; }
        }
        printf("SSCHECK sub=%d: sampled=%ld class-mismatch(on vs off)=%ld  ext maxdiff=%.6g meandiff=%.6g\n",
               sub, sampled, cls, maxd, extboth ? sumd / extboth : 0.0);
        delete[] a; delete[] iter; return 0;
    }

    double best = 1e30;
    for (int t = 0; t < trials; ++t) {
        for (size_t i = 0; i < N; ++i) iter[i] = EMPTYPIXEL;
        Mandel m(W, H, mxit, sub, iter);
        m.setPrecision(precision);
        auto t0 = Clock::now();
        m.Compute(mcx, mcy, msc, mxit, cmethod);
        double s = std::chrono::duration<double>(Clock::now() - t0).count();
        if (s < best) best = s;
    }
    printf("%.4f\n", best);
    delete[] iter;
    return 0;
}
