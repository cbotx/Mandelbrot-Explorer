// Standalone A/B for the custom BigFixed backend vs GMP mpf_t:
//   (1) correctness of mul / add / sub over random O(1) values, and
//   (2) throughput of a raw multiply and of a complex-square reference-orbit
//       loop (the actual Mandelbrot hot path), at the limb counts that matter.
//
// Answers directly: is the custom fixed-point arithmetic correct, and is it
// faster than GMP, at what precision.

#include <gmp.h>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <vector>

#include "bigfixed.h"

using Clock = std::chrono::high_resolution_clock;
static double secs(Clock::time_point a, Clock::time_point b) {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count() / 1e9;
}

static uint64_t rnd64() {
    // 64 random bits from rand() (good enough for test vectors).
    uint64_t x = 0;
    for (int i = 0; i < 4; ++i) x = (x << 16) ^ (uint64_t)(rand() & 0xffff);
    return x;
}

// Random BigFixed with |value| < 4 (integer part 0..3), full random fraction.
static void randBF(BigFixed& a, int L) {
    a.setL(L);
    a.sign = (rand() & 1) ? 1 : -1;
    for (int i = 0; i < L - 1; ++i) a.m[i] = rnd64();
    a.m[L - 1] = (uint64_t)(rand() & 3);
    int nz = 0; for (int i = 0; i < L; ++i) nz |= (a.m[i] != 0);
    if (!nz) a.sign = 0;
}

// |mpf a - mpf b| as a double.
static double mpfAbsDiff(const mpf_t a, const mpf_t b) {
    mpf_t d; mpf_init2(d, mpf_get_prec(a));
    mpf_sub(d, a, b);
    mpf_abs(d, d);
    double r = mpf_get_d(d);
    mpf_clear(d);
    return r;
}

static void correctness(int L, int trials) {
    const mp_bitcnt_t hp = (mp_bitcnt_t)64 * (L + 3);
    double ulp = std::ldexp(1.0, -64 * (L - 1));
    double worst_mul = 0, worst_add = 0, worst_sub = 0;

    mpf_t ma, mb, mref, mgot;
    mpf_init2(ma, hp); mpf_init2(mb, hp); mpf_init2(mref, hp); mpf_init2(mgot, hp);
    std::vector<uint64_t> tmp(2 * L);
    BigFixed a(L), b(L), r(L);

    for (int t = 0; t < trials; ++t) {
        randBF(a, L); randBF(b, L);
        bf_to_mpf(ma, a); bf_to_mpf(mb, b);

        bf_mul(r, a, b, tmp.data());
        bf_to_mpf(mgot, r);
        mpf_mul(mref, ma, mb);
        worst_mul = std::max(worst_mul, mpfAbsDiff(mref, mgot) / ulp);

        bf_add(r, a, b);
        bf_to_mpf(mgot, r);
        mpf_add(mref, ma, mb);
        worst_add = std::max(worst_add, mpfAbsDiff(mref, mgot) / ulp);

        bf_sub(r, a, b);
        bf_to_mpf(mgot, r);
        mpf_sub(mref, ma, mb);
        worst_sub = std::max(worst_sub, mpfAbsDiff(mref, mgot) / ulp);
    }
    printf("  correctness (%d trials): mul %.2f ulp  add %.2f ulp  sub %.2f ulp  %s\n",
           trials, worst_mul, worst_add, worst_sub,
           (worst_mul <= 1.001 && worst_add < 0.001 && worst_sub < 0.001) ? "OK" : "FAIL");
    mpf_clear(ma); mpf_clear(mb); mpf_clear(mref); mpf_clear(mgot);
}

static void benchMul(int L, long iters) {
    std::vector<uint64_t> tmp(2 * L);
    BigFixed a(L), b(L), r(L);
    randBF(a, L); randBF(b, L); a.sign = b.sign = 1;

    // BigFixed: chain r=a*b; a=r to force a dependency.
    auto t0 = Clock::now();
    for (long i = 0; i < iters; ++i) { bf_mul(r, a, b, tmp.data()); a.m.swap(r.m); }
    auto t1 = Clock::now();
    double bf_ns = secs(t0, t1) / iters * 1e9;

    mpf_t ma, mb, mr;
    mpf_init2(ma, (mp_bitcnt_t)64 * L); mpf_init2(mb, (mp_bitcnt_t)64 * L); mpf_init2(mr, (mp_bitcnt_t)64 * L);
    bf_to_mpf(ma, a); bf_to_mpf(mb, b);
    t0 = Clock::now();
    for (long i = 0; i < iters; ++i) { mpf_mul(mr, ma, mb); mpf_set(ma, mr); }
    t1 = Clock::now();
    double gmp_ns = secs(t0, t1) / iters * 1e9;
    mpf_clear(ma); mpf_clear(mb); mpf_clear(mr);

    printf("  mul  L=%2d (%4d bit): BigFixed %7.1f ns   GMP %7.1f ns   speedup %.2fx\n",
           L, 64 * L, bf_ns, gmp_ns, gmp_ns / bf_ns);
}

// Complex square z = z^2 + c, the reference-orbit inner step.
static void benchOrbit(int L, long iters) {
    std::vector<uint64_t> tmp(2 * L);
    BigFixed zr(L), zi(L), cr(L), ci(L), t1v(L), t2v(L), t3v(L);
    zr.setInt(0); zi.setInt(0);
    // c = -0.5 + 0i  (stays bounded forever).
    cr.setL(L); cr.sign = -1; cr.m[L - 1] = 0; cr.m[L - 2] = 0x8000000000000000ull; // 0.5
    ci.setInt(0);

    auto t0 = Clock::now();
    for (long i = 0; i < iters; ++i) {
        bf_mul(t1v, zr, zr, tmp.data());   // zr^2
        bf_mul(t2v, zi, zi, tmp.data());   // zi^2
        bf_mul(t3v, zr, zi, tmp.data());   // zr*zi
        bf_sub(t1v, t1v, t2v);             // zr^2 - zi^2
        bf_add(zr, t1v, cr);               // zr' = zr^2 - zi^2 + cr
        bf_add(t3v, t3v, t3v);             // 2 zr zi
        bf_add(zi, t3v, ci);               // zi' = 2 zr zi + ci
    }
    auto t1 = Clock::now();
    double bf_ns = secs(t0, t1) / iters * 1e9;

    mpf_t mzr, mzi, mcr, mci, a, bb, cc;
    mp_bitcnt_t p = (mp_bitcnt_t)64 * L;
    mpf_init2(mzr, p); mpf_init2(mzi, p); mpf_init2(mcr, p); mpf_init2(mci, p);
    mpf_init2(a, p); mpf_init2(bb, p); mpf_init2(cc, p);
    mpf_set_ui(mzr, 0); mpf_set_ui(mzi, 0); mpf_set_d(mcr, -0.5); mpf_set_ui(mci, 0);
    t0 = Clock::now();
    for (long i = 0; i < iters; ++i) {
        mpf_mul(a, mzr, mzr);
        mpf_mul(bb, mzi, mzi);
        mpf_mul(cc, mzr, mzi);
        mpf_sub(a, a, bb);
        mpf_add(mzr, a, mcr);
        mpf_add(cc, cc, cc);
        mpf_add(mzi, cc, mci);
    }
    t1 = Clock::now();
    double gmp_ns = secs(t0, t1) / iters * 1e9;
    mpf_clear(mzr); mpf_clear(mzi); mpf_clear(mcr); mpf_clear(mci);
    mpf_clear(a); mpf_clear(bb); mpf_clear(cc);

    printf("  orbit L=%2d (%4d bit): BigFixed %7.1f ns   GMP %7.1f ns   speedup %.2fx\n",
           L, 64 * L, bf_ns, gmp_ns, gmp_ns / bf_ns);
}

int main() {
    srand(12345);
    int Ls[] = { 3, 9, 24, 48 };   // ~1e30, 1e150, 1e450, 1e900 zoom depths

    printf("== correctness (BigFixed vs GMP mpf) ==\n");
    for (int L : Ls) { printf(" L=%d:\n", L); correctness(L, 20000); }

    printf("\n== raw multiply throughput ==\n");
    for (int L : Ls) benchMul(L, 2000000);

    printf("\n== complex-square reference-orbit step ==\n");
    for (int L : Ls) benchOrbit(L, 1000000);
    return 0;
}
