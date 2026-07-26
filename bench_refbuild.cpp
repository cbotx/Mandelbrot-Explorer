// Microbenchmark: reference-orbit complex squaring z^2+c at fixed precision.
// Compares the current mpf approach against raw mpn (fixed-point) variants to
// estimate the payoff of hand-rolling GMP for the Mandelbrot reference build.
//
//   A: mpf Karatsuba complex square (matches calCoefficient)
//   B: mpn Karatsuba  (a+b)(a-b), a*b  -> 2 mpn_mul_n
//   C: mpn direct     a^2, b^2, a*b    -> 2 mpn_sqr + 1 mpn_mul_n
//
// All at L limbs; the orbit values are irrelevant to schoolbook-mul timing, so we
// iterate on fixed data and accumulate a checksum to defeat dead-code elimination.
#include <gmp.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <chrono>

static double nows() {
    return std::chrono::duration_cast<std::chrono::duration<double>>(
        std::chrono::high_resolution_clock::now().time_since_epoch()).count();
}

int main(int argc, char** argv) {
    int bits = argc > 1 ? atoi(argv[1]) : 1502;   // deep432 ~= 1502 bits
    int N    = argc > 2 ? atoi(argv[2]) : 132398;  // deep432 reflen
    int reps = argc > 3 ? atoi(argv[3]) : 20;      // repeat for stable timing
    int L = (bits + 63) / 64;                       // limbs
    printf("bits=%d  L=%d limbs  N=%d  reps=%d\n", bits, L, N, reps);

    // ---- Method A: mpf ----
    mpf_set_default_prec(bits);
    mpf_t zr, zi, cr, ci, t1, t2, nre, nim;
    mpf_init(zr); mpf_init(zi); mpf_init(cr); mpf_init(ci);
    mpf_init(t1); mpf_init(t2); mpf_init(nre); mpf_init(nim);
    // bounded-ish operands
    mpf_set_d(cr, -0.75); mpf_set_d(ci, 0.06);
    double bestA = 1e30;
    for (int r = 0; r < reps; ++r) {
        mpf_set_d(zr, 0.1); mpf_set_d(zi, 0.2);
        double t0 = nows();
        for (int i = 0; i < N; ++i) {
            // Karatsuba complex square: re=(a+b)(a-b), im=2ab, +c
            mpf_add(t1, zr, zi);
            mpf_sub(t2, zr, zi);
            mpf_mul(nim, zr, zi);
            mpf_mul(nre, t1, t2);
            mpf_mul_ui(nim, nim, 2);
            mpf_add(nre, nre, cr);
            mpf_add(nim, nim, ci);
            // keep bounded: halve (cheap, mimics staying in the disk)
            mpf_div_ui(zr, nre, 2);
            mpf_div_ui(zi, nim, 2);
        }
        double dt = nows() - t0;
        if (dt < bestA) bestA = dt;
    }
    double csA = mpf_get_d(zr) + mpf_get_d(zi);

    // ---- mpn helpers: fixed L-limb signed magnitudes ----
    // value = sign * (mant / 2^(64L-2));  |value| < 2.  Product truncates to high L.
    mp_limb_t *a = new mp_limb_t[L], *b = new mp_limb_t[L];
    mp_limb_t *Cr = new mp_limb_t[L], *Ci = new mp_limb_t[L];
    mp_limb_t *prod = new mp_limb_t[2 * L];
    mp_limb_t *s1 = new mp_limb_t[L + 1], *s2 = new mp_limb_t[L + 1];
    mp_limb_t *re = new mp_limb_t[L], *im = new mp_limb_t[L];
    for (int i = 0; i < L; ++i) { a[i] = 0x9e3779b97f4a7c15ULL * (i + 1); b[i] = 0xC2B2AE3D27D4EB4FULL * (i + 3); }
    a[L - 1] &= 0x3fffffffffffffffULL; b[L - 1] &= 0x3fffffffffffffffULL;   // |.|<2
    memcpy(Cr, a, L * sizeof(mp_limb_t)); memcpy(Ci, b, L * sizeof(mp_limb_t));

    // ---- Method B: mpn Karatsuba (2 mul) ----
    double bestB = 1e30;
    mp_limb_t sink = 0;
    for (int r = 0; r < reps; ++r) {
        mp_limb_t x[64], y[64];
        memcpy(x, a, L * sizeof(mp_limb_t)); memcpy(y, b, L * sizeof(mp_limb_t));
        double t0 = nows();
        for (int i = 0; i < N; ++i) {
            mpn_add_n(s1, x, y, L);              // a+b   (ignore top carry for timing)
            mpn_sub_n(s2, x, y, L);              // a-b
            mpn_mul_n(prod, s1, s2, L);          // (a+b)(a-b) -> 2L
            memcpy(re, prod + L, L * sizeof(mp_limb_t));   // truncate high L
            mpn_mul_n(prod, x, y, L);            // a*b -> 2L
            memcpy(im, prod + L, L * sizeof(mp_limb_t));
            mpn_lshift(im, im, L, 1);            // 2ab
            mpn_add_n(re, re, Cr, L);            // +c
            mpn_add_n(im, im, Ci, L);
            mpn_rshift(x, re, L, 1);             // keep bounded
            mpn_rshift(y, im, L, 1);
        }
        double dt = nows() - t0;
        if (dt < bestB) bestB = dt;
        sink ^= x[0] ^ y[L - 1];
    }

    // ---- Method C: mpn direct (2 sqr + 1 mul) ----
    double bestC = 1e30;
    for (int r = 0; r < reps; ++r) {
        mp_limb_t x[64], y[64], sa[128], sb[128];
        memcpy(x, a, L * sizeof(mp_limb_t)); memcpy(y, b, L * sizeof(mp_limb_t));
        double t0 = nows();
        for (int i = 0; i < N; ++i) {
            mpn_sqr(sa, x, L);                   // a^2 -> 2L
            mpn_sqr(sb, y, L);                   // b^2 -> 2L
            mpn_sub_n(re, sa + L, sb + L, L);    // a^2-b^2 (high L)
            mpn_mul_n(prod, x, y, L);            // a*b
            memcpy(im, prod + L, L * sizeof(mp_limb_t));
            mpn_lshift(im, im, L, 1);
            mpn_add_n(re, re, Cr, L);
            mpn_add_n(im, im, Ci, L);
            mpn_rshift(x, re, L, 1);
            mpn_rshift(y, im, L, 1);
        }
        double dt = nows() - t0;
        if (dt < bestC) bestC = dt;
        sink ^= x[0] ^ y[L - 1];
    }

    // ---- Method D: FUSED high-half short products ----
    // mul_hi(hi, x, y, L, g): high L limbs of the magnitude product x*y, computing
    // only product limbs >= L-1-g (g guard limbs). Drops the low ~half of the
    // schoolbook partial products that we truncate away anyway. Exact for limbs
    // > L-1-g; limb L-1-g may miss a tiny carry (<~L) that reaches the lowest kept
    // output limb only with probability ~2^-64g. g=2 => effectively exact.
    auto mul_hi = [](mp_limb_t* hi, const mp_limb_t* x, const mp_limb_t* y, int L, int g) {
        int base = L - 1 - g;                     // lowest product limb we accumulate
        int accsz = L + g + 1;                     // acc[k] = product limb base+k, k=0..L+g
        mp_limb_t acc[160];
        memset(acc, 0, accsz * sizeof(mp_limb_t));
        for (int i = 0; i < L; ++i) {
            int jlo = base - i; if (jlo < 0) jlo = 0;
            int len = L - jlo;
            int accidx = (i + jlo) - base;         // = max(0, i-base)
            mp_limb_t cy = mpn_addmul_1(acc + accidx, y + jlo, len, x[i]);
            if (accidx + len < accsz) mpn_add_1(acc + accidx + len, acc + accidx + len, accsz - (accidx + len), cy);
        }
        memcpy(hi, acc + (g + 1), L * sizeof(mp_limb_t));   // product limbs L..2L-1
    };
    double bestD = 1e30;
    for (int r = 0; r < reps; ++r) {
        mp_limb_t x[64], y[64], s1_[64], s2_[64], hre[64], him[64];
        memcpy(x, a, L * sizeof(mp_limb_t)); memcpy(y, b, L * sizeof(mp_limb_t));
        double t0 = nows();
        for (int i = 0; i < N; ++i) {
            // re = (a+b)(a-b), im = 2ab  via two high-half products, +c, /2 to bound.
            mpn_add_n(s1_, x, y, L);              // a+b
            mpn_sub_n(s2_, x, y, L);              // a-b (magnitude ignoring sign for timing)
            mul_hi(hre, s1_, s2_, L, 2);          // hi((a+b)(a-b))
            mul_hi(him, x, y, L, 2);              // hi(a*b)
            mpn_lshift(him, him, L, 1);           // 2ab
            mpn_add_n(hre, hre, Cr, L);           // +cr
            mpn_add_n(him, him, Ci, L);           // +ci
            mpn_rshift(x, hre, L, 1);
            mpn_rshift(y, him, L, 1);
        }
        double dt = nows() - t0;
        if (dt < bestD) bestD = dt;
        sink ^= x[0] ^ y[L - 1];
    }

    // Correctness spot-check: mul_hi vs full mpn_mul_n high half.
    {
        mp_limb_t full[128], hi[64];
        mpn_mul_n(full, a, b, L);
        auto mh = [&](mp_limb_t* o, const mp_limb_t* X, const mp_limb_t* Y, int g){
            int base=L-1-g, accsz=L+g+1; mp_limb_t acc[160]; memset(acc,0,accsz*sizeof(mp_limb_t));
            for(int i=0;i<L;++i){int jlo=base-i; if(jlo<0)jlo=0; int len=L-jlo; int ai=(i+jlo)-base;
                mp_limb_t cy=mpn_addmul_1(acc+ai,Y+jlo,len,X[i]); if(ai+len<accsz) mpn_add_1(acc+ai+len,acc+ai+len,accsz-(ai+len),cy);}
            memcpy(o,acc+(g+1),L*sizeof(mp_limb_t)); };
        mh(hi, a, b, 2);
        long long maxerr = 0;
        for (int k = 0; k < L; ++k) { long long e = (long long)(hi[k] - full[L + k]); if (e < 0) e = -e; if (e > maxerr) maxerr = e; }
        printf("[mul_hi vs full high-half: max limb diff = %lld]\n", maxerr);
    }

    printf("A mpf Karatsuba : %.4f s  (%.1f ns/iter)\n", bestA, bestA / N * 1e9);
    printf("B mpn Karatsuba : %.4f s  (%.1f ns/iter)  speedup %.2fx\n", bestB, bestB / N * 1e9, bestA / bestB);
    printf("C mpn 2sqr+1mul : %.4f s  (%.1f ns/iter)  speedup %.2fx\n", bestC, bestC / N * 1e9, bestA / bestC);
    printf("D fused hi-half : %.4f s  (%.1f ns/iter)  speedup %.2fx\n", bestD, bestD / N * 1e9, bestA / bestD);
    printf("[csum A=%.3f sink=%llu]\n", csA, (unsigned long long)sink);
    return 0;
}
