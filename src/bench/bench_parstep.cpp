// Quantifies WHY a single deep-zoom reference computation should not be threaded:
// compares the time of one big-integer multiply at 1e1000 depth (L~52 limbs) and
// one full reference-orbit step against the raw cost of an OpenMP thread hand-off
// (an empty parallel region / barrier). If the multiply is far cheaper than the
// hand-off, splitting it across threads is a guaranteed loss.
#include <gmp.h>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <chrono>
#include <vector>
#include <omp.h>

using Clock = std::chrono::high_resolution_clock;
static double ns_per(Clock::time_point a, Clock::time_point b, long n) {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count() / (double)n;
}
static uint64_t rnd64() {
    uint64_t x = 0;
    for (int i = 0; i < 4; ++i) x = (x << 16) ^ (uint64_t)(rand() & 0xffff);
    return x;
}

int main() {
    srand(7);
    const int L = 52; // ~1e1000 depth (3328-bit reference)
    std::vector<uint64_t> a(L), b(L), full(2 * L);
    for (int i = 0; i < L; ++i) {
        a[i] = rnd64();
        b[i] = rnd64();
    }
    const long N = 3000000;

    // 1) one big multiply
    auto t0 = Clock::now();
    for (long k = 0; k < N; ++k) {
        mpn_mul_n((mp_limb_t*)full.data(), (mp_limb_t*)a.data(), (mp_limb_t*)b.data(), L);
        a[0] ^= full[L];
    }
    auto t1 = Clock::now();
    double mul_ns = ns_per(t0, t1, N);

    // 2) one full complex-square orbit step (2 muls Karatsuba-style + adds) via mpf
    mpf_t zr, zi, cr, ci, t, u, v;
    mp_bitcnt_t p = (mp_bitcnt_t)64 * L;
    mpf_init2(zr, p);
    mpf_init2(zi, p);
    mpf_init2(cr, p);
    mpf_init2(ci, p);
    mpf_init2(t, p);
    mpf_init2(u, p);
    mpf_init2(v, p);
    mpf_set_d(zr, 0.1);
    mpf_set_d(zi, 0.2);
    mpf_set_d(cr, -0.5);
    mpf_set_d(ci, 0.0);
    t0 = Clock::now();
    for (long k = 0; k < N; ++k) {
        mpf_add(t, zr, zi);
        mpf_sub(u, zr, zi);
        mpf_mul(v, t, u);   // zr^2 - zi^2
        mpf_mul(t, zr, zi); // zr*zi
        mpf_add(zr, v, cr);
        mpf_add(t, t, t);
        mpf_add(zi, t, ci);
    }
    t1 = Clock::now();
    double step_ns = ns_per(t0, t1, N);
    mpf_clear(zr);
    mpf_clear(zi);
    mpf_clear(cr);
    mpf_clear(ci);
    mpf_clear(t);
    mpf_clear(u);
    mpf_clear(v);

    // 3) OpenMP thread hand-off: empty parallel region (fork/join + barrier)
    const long M = 300000;
    volatile long sink = 0;
    t0 = Clock::now();
    for (long k = 0; k < M; ++k) {
#pragma omp parallel
        { sink += omp_get_thread_num(); }
    }
    t1 = Clock::now();
    double fork_ns = ns_per(t0, t1, M);

    // 4) just a barrier inside a persistent region
    double bar_ns = 0;
#pragma omp parallel
    {
        auto b0 = Clock::now();
        for (long k = 0; k < M; ++k) {
#pragma omp barrier
        }
        auto b1 = Clock::now();
#pragma omp master
        bar_ns = ns_per(b0, b1, M);
    }

    printf("threads available: %d\n", omp_get_max_threads());
    printf("  one mpn_mul_n  (L=52, 3328-bit)     : %8.1f ns\n", mul_ns);
    printf("  one orbit step (2 mul + adds, mpf)  : %8.1f ns\n", step_ns);
    printf("  OpenMP fork/join (empty region)     : %8.1f ns   = %.1fx one multiply\n", fork_ns, fork_ns / mul_ns);
    printf("  OpenMP barrier   (persistent team)  : %8.1f ns   = %.1fx one multiply\n", bar_ns, bar_ns / mul_ns);
    return 0;
}
