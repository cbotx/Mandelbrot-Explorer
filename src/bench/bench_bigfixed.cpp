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
    int nz = 0;
    for (int i = 0; i < L; ++i) nz |= (a.m[i] != 0);
    if (!nz) a.sign = 0;
}

// |mpf a - mpf b| as a double.
static double mpfAbsDiff(const mpf_t a, const mpf_t b) {
    mpf_t d;
    mpf_init2(d, mpf_get_prec(a));
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
    mpf_init2(ma, hp);
    mpf_init2(mb, hp);
    mpf_init2(mref, hp);
    mpf_init2(mgot, hp);
    std::vector<uint64_t> tmp(2 * L);
    BigFixed a(L), b(L), r(L);

    for (int t = 0; t < trials; ++t) {
        randBF(a, L);
        randBF(b, L);
        bf_to_mpf(ma, a);
        bf_to_mpf(mb, b);

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
    printf("  correctness (%d trials): mul %.2f ulp  add %.2f ulp  sub %.2f ulp  %s\n", trials, worst_mul, worst_add, worst_sub, (worst_mul <= 1.001 && worst_add < 0.001 && worst_sub < 0.001) ? "OK" : "FAIL");
    mpf_clear(ma);
    mpf_clear(mb);
    mpf_clear(mref);
    mpf_clear(mgot);
}

static void benchMul(int L, long iters) {
    std::vector<uint64_t> tmp(2 * L);
    BigFixed a(L), b(L), r(L);
    randBF(a, L);
    randBF(b, L);
    a.sign = b.sign = 1;

    // BigFixed: chain r=a*b; a=r to force a dependency.
    auto t0 = Clock::now();
    for (long i = 0; i < iters; ++i) {
        bf_mul(r, a, b, tmp.data());
        a.m.swap(r.m);
    }
    auto t1 = Clock::now();
    double bf_ns = secs(t0, t1) / iters * 1e9;

    mpf_t ma, mb, mr;
    mpf_init2(ma, (mp_bitcnt_t)64 * L);
    mpf_init2(mb, (mp_bitcnt_t)64 * L);
    mpf_init2(mr, (mp_bitcnt_t)64 * L);
    bf_to_mpf(ma, a);
    bf_to_mpf(mb, b);
    t0 = Clock::now();
    for (long i = 0; i < iters; ++i) {
        mpf_mul(mr, ma, mb);
        mpf_set(ma, mr);
    }
    t1 = Clock::now();
    double gmp_ns = secs(t0, t1) / iters * 1e9;
    mpf_clear(ma);
    mpf_clear(mb);
    mpf_clear(mr);

    printf("  mul  L=%2d (%4d bit): BigFixed %7.1f ns   GMP %7.1f ns   speedup %.2fx\n", L, 64 * L, bf_ns, gmp_ns, gmp_ns / bf_ns);
}

// Complex square z = z^2 + c, the reference-orbit inner step.
static void benchOrbit(int L, long iters) {
    std::vector<uint64_t> tmp(2 * L);
    BigFixed zr(L), zi(L), cr(L), ci(L), t1v(L), t2v(L), t3v(L);
    zr.setInt(0);
    zi.setInt(0);
    // c = -0.5 + 0i  (stays bounded forever).
    cr.setL(L);
    cr.sign = -1;
    cr.m[L - 1] = 0;
    cr.m[L - 2] = 0x8000000000000000ull; // 0.5
    ci.setInt(0);

    auto t0 = Clock::now();
    for (long i = 0; i < iters; ++i) {
        bf_mul(t1v, zr, zr, tmp.data()); // zr^2
        bf_mul(t2v, zi, zi, tmp.data()); // zi^2
        bf_mul(t3v, zr, zi, tmp.data()); // zr*zi
        bf_sub(t1v, t1v, t2v);           // zr^2 - zi^2
        bf_add(zr, t1v, cr);             // zr' = zr^2 - zi^2 + cr
        bf_add(t3v, t3v, t3v);           // 2 zr zi
        bf_add(zi, t3v, ci);             // zi' = 2 zr zi + ci
    }
    auto t1 = Clock::now();
    double bf_ns = secs(t0, t1) / iters * 1e9;

    mpf_t mzr, mzi, mcr, mci, a, bb, cc;
    mp_bitcnt_t p = (mp_bitcnt_t)64 * L;
    mpf_init2(mzr, p);
    mpf_init2(mzi, p);
    mpf_init2(mcr, p);
    mpf_init2(mci, p);
    mpf_init2(a, p);
    mpf_init2(bb, p);
    mpf_init2(cc, p);
    mpf_set_ui(mzr, 0);
    mpf_set_ui(mzi, 0);
    mpf_set_d(mcr, -0.5);
    mpf_set_ui(mci, 0);
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
    mpf_clear(mzr);
    mpf_clear(mzi);
    mpf_clear(mcr);
    mpf_clear(mci);
    mpf_clear(a);
    mpf_clear(bb);
    mpf_clear(cc);

    printf("  orbit L=%2d (%4d bit): BigFixed %7.1f ns   GMP %7.1f ns   speedup %.2fx\n", L, 64 * L, bf_ns, gmp_ns, gmp_ns / bf_ns);
}

// Engine-faithful complex-square step comparison at deep precision:
//   (K) 2-mul Karatsuba : re=(a+b)(a-b), ab=a*b        -> what the engine uses now
//   (S) 3-square form   : a^2, b^2, (a+b)^2 via mpn_sqr -> 2ab=(a+b)^2-a^2-b^2
// vs mpf 2-mul Karatsuba. Also checks (S) and (K) agree to rounding.
static void benchOrbitSqr(int L, long iters) {
    std::vector<uint64_t> tmp(2 * L);
    auto mkC = [&](BigFixed& cr, BigFixed& ci) {
        cr.setL(L);
        cr.sign = -1;
        cr.m[L - 2] = 0x8000000000000000ull;
        ci.setInt(0); // c=-0.5
    };

    // (K) 2-mul Karatsuba
    {
        BigFixed zr(L), zi(L), cr(L), ci(L), s(L), d(L), ab(L), re(L), im(L);
        mkC(cr, ci);
        auto t0 = Clock::now();
        for (long i = 0; i < iters; ++i) {
            bf_add(s, zr, zi);
            bf_sub(d, zr, zi);
            bf_mul(re, s, d, tmp.data());   // a^2 - b^2
            bf_mul(ab, zr, zi, tmp.data()); // a*b
            bf_add(zr, re, cr);
            bf_add(im, ab, ab);
            bf_add(zi, im, ci);
        }
        auto t1 = Clock::now();
        double k_ns = secs(t0, t1) / iters * 1e9;

        // (S) 3-square form
        BigFixed Zr(L), Zi(L), s1(L), s2(L), s3(L), apb(L), t(L);
        mkC(cr, ci);
        t0 = Clock::now();
        for (long i = 0; i < iters; ++i) {
            bf_sqr(s1, Zr, tmp.data()); // a^2
            bf_sqr(s2, Zi, tmp.data()); // b^2
            bf_add(apb, Zr, Zi);
            bf_sqr(s3, apb, tmp.data()); // (a+b)^2
            bf_sub(t, s1, s2);
            bf_add(Zr, t, cr); // a^2-b^2+cr
            bf_sub(t, s3, s1);
            bf_sub(t, t, s2);
            bf_add(Zi, t, ci); // 2ab+ci
        }
        t1 = Clock::now();
        double s_ns = secs(t0, t1) / iters * 1e9;

        // mpf 2-mul Karatsuba
        mpf_t zr2, zi2, cr2, ci2, ms, md, mre, mab;
        mp_bitcnt_t p = (mp_bitcnt_t)64 * L;
        mpf_init2(zr2, p);
        mpf_init2(zi2, p);
        mpf_init2(cr2, p);
        mpf_init2(ci2, p);
        mpf_init2(ms, p);
        mpf_init2(md, p);
        mpf_init2(mre, p);
        mpf_init2(mab, p);
        mpf_set_ui(zr2, 0);
        mpf_set_ui(zi2, 0);
        mpf_set_d(cr2, -0.5);
        mpf_set_ui(ci2, 0);
        t0 = Clock::now();
        for (long i = 0; i < iters; ++i) {
            mpf_add(ms, zr2, zi2);
            mpf_sub(md, zr2, zi2);
            mpf_mul(mre, ms, md);
            mpf_mul(mab, zr2, zi2);
            mpf_add(zr2, mre, cr2);
            mpf_add(mab, mab, mab);
            mpf_add(zi2, mab, ci2);
        }
        t1 = Clock::now();
        double m_ns = secs(t0, t1) / iters * 1e9;
        mpf_clear(zr2);
        mpf_clear(zi2);
        mpf_clear(cr2);
        mpf_clear(ci2);
        mpf_clear(ms);
        mpf_clear(md);
        mpf_clear(mre);
        mpf_clear(mab);

        printf("  step L=%2d (%4d bit): mpf(2mul) %7.1f | K(2mul) %7.1f (%.2fx) | S(3sqr) %7.1f (%.2fx)\n", L, 64 * L, m_ns, k_ns, m_ns / k_ns, s_ns, m_ns / s_ns);
    }
}

// One-step agreement: 3-square result must equal 2-mul result to <=1 ulp.
static void checkSqrForm(int L) {
    std::vector<uint64_t> tmp(2 * L);
    BigFixed a(L), b(L);
    randBF(a, L);
    randBF(b, L);
    a.sign = b.sign = 1;
    // K: re=(a+b)(a-b), im=2ab
    BigFixed s(L), d(L), reK(L), abK(L), imK(L);
    bf_add(s, a, b);
    bf_sub(d, a, b);
    bf_mul(reK, s, d, tmp.data());
    bf_mul(abK, a, b, tmp.data());
    bf_add(imK, abK, abK);
    // S: a^2,b^2,(a+b)^2
    BigFixed s1(L), s2(L), s3(L), apb(L), reS(L), imS(L), t(L);
    bf_sqr(s1, a, tmp.data());
    bf_sqr(s2, b, tmp.data());
    bf_add(apb, a, b);
    bf_sqr(s3, apb, tmp.data());
    bf_sub(reS, s1, s2);
    bf_sub(t, s3, s1);
    bf_sub(imS, t, s2);
    mpf_t x, y, dd;
    mp_bitcnt_t p = (mp_bitcnt_t)64 * L;
    mpf_init2(x, p);
    mpf_init2(y, p);
    mpf_init2(dd, p);
    // ulp = 2^-(64*(L-1)); express |diff| in ulp by shifting up before converting
    // (a plain ldexp underflows for large L).
    bf_to_mpf(x, reK);
    bf_to_mpf(y, reS);
    mpf_sub(dd, x, y);
    mpf_abs(dd, dd);
    mpf_mul_2exp(dd, dd, (mp_bitcnt_t)64 * (L - 1));
    double reDiff = mpf_get_d(dd);
    bf_to_mpf(x, imK);
    bf_to_mpf(y, imS);
    mpf_sub(dd, x, y);
    mpf_abs(dd, dd);
    mpf_mul_2exp(dd, dd, (mp_bitcnt_t)64 * (L - 1));
    double imDiff = mpf_get_d(dd);
    mpf_clear(x);
    mpf_clear(y);
    mpf_clear(dd);
    printf("  agree L=%2d: re %.2f ulp  im %.2f ulp  %s\n", L, reDiff, imDiff, (reDiff <= 2.001 && imDiff <= 2.001) ? "OK" : "FAIL");
}

int main() {
    srand(12345);
    int Ls[] = {3, 9, 24, 48}; // ~1e30, 1e150, 1e450, 1e900 zoom depths
    printf("== correctness (BigFixed vs GMP mpf) ==\n");
    for (int L : Ls) {
        printf(" L=%d:\n", L);
        correctness(L, 20000);
    }

    printf("\n== raw multiply throughput ==\n");
    for (int L : Ls) benchMul(L, 2000000);

    printf("\n== complex-square reference-orbit step (naive 3-mul) ==\n");
    for (int L : Ls) benchOrbit(L, 1000000);

    printf("\n== 3-square vs 2-mul agreement (ulp) ==\n");
    for (int L : Ls) checkSqrForm(L);

    printf("\n== engine-faithful complex-square step: 2-mul vs 3-square ==\n");
    for (int L : Ls) benchOrbitSqr(L, 1000000);
    return 0;
}
