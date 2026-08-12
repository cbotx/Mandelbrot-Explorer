// Extreme-depth reference-orbit multiply study. At deep zoom the orbit step is
// dominated by big-integer multiplies, and BigFixed is forced to call GMP's own
// mpn_mul_n -- the same routine mpf_mul uses -- so per-multiply we cannot beat
// GMP. The ONE Mandelbrot-specific edge: mpf_mul computes the full 2L-limb
// product then discards the low half, but the fixed-point orbit only needs the
// high L limbs (plus one guard limb for rounding). A high-half ("short") product
// does less work. This bench quantifies:
//   (1) mpn_mul_n (full)  vs  mpf_mul  vs  mpn_sqr  vs a schoolbook high-half,
//   (2) the full engine-faithful 2-mul Karatsuba complex-square step,
// at the limb counts that matter for extreme depth (L = 16..64, ~1e300..1e1200).

#include <gmp.h>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <chrono>
#include <vector>
#include <intrin.h>
#include <immintrin.h> // _mulx_u64 (BMI2), _addcarryx_u64 (ADX)

using Clock = std::chrono::high_resolution_clock;
static double ns_per(Clock::time_point a, Clock::time_point b, long n) {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count() / (double)n;
}
static uint64_t rnd64() {
    uint64_t x = 0;
    for (int i = 0; i < 4; ++i) x = (x << 16) ^ (uint64_t)(rand() & 0xffff);
    return x;
}

// Schoolbook HIGH-half product: keep only limbs [L-2 .. 2L-1] of the 2L product,
// i.e. compute columns k = i+j >= L-2. r must hold at least L+2 limbs, indexed so
// r[t] corresponds to product limb (L-2 + t). Uses _umul128 + add-carry.
static inline void mulhigh_school(uint64_t* r, const uint64_t* a, const uint64_t* b, int L) {
    const int lo = L - 2;        // lowest product column we keep
    const int nout = 2 * L - lo; // number of kept columns
    std::memset(r, 0, sizeof(uint64_t) * (size_t)nout);
    for (int i = 0; i < L; ++i) {
        uint64_t ai = a[i];
        if (!ai) continue;
        int jstart = lo - i;
        if (jstart < 0) jstart = 0;
        uint64_t carry = 0;
        for (int j = jstart; j < L; ++j) {
            uint64_t hi, plo = _umul128(ai, b[j], &hi);
            int t = (i + j) - lo;
            unsigned char c = _addcarry_u64(0, plo, carry, &plo);
            hi += c;
            c = _addcarry_u64(0, r[t], plo, &r[t]);
            carry = hi + c;
        }
        int tcarry = (i + L) - lo;
        if (tcarry < nout) {
            unsigned char c = _addcarry_u64(0, r[tcarry], carry, &r[tcarry]);
            for (int t = tcarry + 1; c && t < nout; ++t) c = _addcarry_u64(c, r[t], 0, &r[t]);
        }
    }
}

// Tuned HIGH-half product using MULX (BMI2) + ADCX/ADOX (ADX, two carry chains).
// Keeps product columns [lo .. 2L-1], lo = max(0, L-1-GUARD); r[t] = column lo+t,
// needs nout = 2L-lo limbs (+2 headroom). The GUARD dropped-carry columns absorb
// carries from the discarded low half so the kept high L limbs are exact for
// GUARD >= ~2. Does ~half the partial products of a full multiply.
static inline void mulhigh_adx(uint64_t* r, const uint64_t* a, const uint64_t* b, int L, int GUARD) {
    int lo = L - 1 - GUARD;
    if (lo < 0) lo = 0;
    const int nout = 2 * L - lo;
    std::memset(r, 0, sizeof(uint64_t) * (size_t)(nout + 2));
    for (int i = 0; i < L; ++i) {
        uint64_t ai = a[i];
        int jstart = lo - i;
        if (jstart < 0) jstart = 0;
        unsigned char cf = 0, of = 0;
        int t = i + jstart - lo;
        for (int j = jstart; j < L; ++j, ++t) {
            unsigned long long hi, plo = _mulx_u64(ai, b[j], &hi);
            cf = _addcarryx_u64(cf, r[t], plo, &r[t]);        // CF chain: low halves
            of = _addcarryx_u64(of, r[t + 1], hi, &r[t + 1]); // OF chain: high halves
        }
        // t is now t_last+1; flush the two trailing carries.
        unsigned char c = _addcarry_u64(0, r[t], cf, &r[t]);
        for (int u = t + 1; c && u <= nout; ++u) c = _addcarry_u64(c, r[u], 0, &r[u]);
        c = _addcarry_u64(0, r[t + 1], of, &r[t + 1]);
        for (int u = t + 2; c && u <= nout; ++u) c = _addcarry_u64(c, r[u], 0, &r[u]);
    }
}

// HIGH-half built from GMP's OWN tuned asm primitive mpn_addmul_1, applied per
// row but only over the kept high columns (>= lo). This gets GMP's per-limb
// throughput while doing only ~half the partial products of a full multiply.
static inline void mulhigh_addmul(uint64_t* r, const uint64_t* a, const uint64_t* b, int L, int GUARD) {
    int lo = L - 1 - GUARD;
    if (lo < 0) lo = 0;
    const int nout = 2 * L - lo;
    std::memset(r, 0, sizeof(uint64_t) * (size_t)(nout + 1));
    for (int i = 0; i < L; ++i) {
        int jstart = lo - i;
        if (jstart < 0) jstart = 0;
        int n = L - jstart;
        int t = i + jstart - lo;
        mp_limb_t carry = mpn_addmul_1((mp_limb_t*)(r + t), (const mp_limb_t*)(b + jstart), (mp_size_t)n, (mp_limb_t)a[i]);
        int u = t + n;
        mp_limb_t cc = carry;
        while (cc && u <= nout) {
            unsigned char c = _addcarry_u64(0, r[u], cc, &r[u]);
            cc = c;
            ++u;
        }
    }
}

// HIGH-half SQUARE a^2 exploiting symmetry: off-diagonal products a[i]a[j] (i<j)
// summed once via mpn_addmul_1 over high columns, doubled, then diagonal a[i]^2
// added. ~half the partial products of a high-half general multiply. r indexed as
// mulhigh_addmul (r[t] = column lo+t); exact top L limbs for GUARD>=~2.
static inline void sqrhigh_addmul(uint64_t* r, const uint64_t* a, int L, int GUARD) {
    int lo = L - 1 - GUARD;
    if (lo < 0) lo = 0;
    const int nout = 2 * L - lo;
    std::memset(r, 0, sizeof(uint64_t) * (size_t)(nout + 2));
    // off-diagonal sum_{i<j} a[i]a[j], high columns only
    for (int i = 0; i < L - 1; ++i) {
        int jstart = i + 1;
        if (i + jstart < lo) jstart = lo - i;
        if (jstart < i + 1) jstart = i + 1;
        int n = L - jstart;
        if (n <= 0) continue;
        int t = i + jstart - lo;
        mp_limb_t carry = mpn_addmul_1((mp_limb_t*)(r + t), (const mp_limb_t*)(a + jstart), (mp_size_t)n, (mp_limb_t)a[i]);
        int u = t + n;
        mp_limb_t cc = carry;
        while (cc && u <= nout) {
            unsigned char c = _addcarry_u64(0, r[u], cc, &r[u]);
            cc = c;
            ++u;
        }
    }
    mpn_lshift((mp_limb_t*)r, (const mp_limb_t*)r, (mp_size_t)(nout + 1), 1); // double
    // diagonal a[i]^2 at column 2i (both limbs), high only
    for (int i = 0; i < L; ++i) {
        int col = 2 * i;
        if (col + 1 < lo) continue;
        unsigned long long hi, lo64 = _mulx_u64(a[i], a[i], &hi);
        if (col >= lo) {
            int t = col - lo;
            unsigned char c = _addcarry_u64(0, r[t], lo64, &r[t]);
            c = _addcarry_u64(c, r[t + 1], hi, &r[t + 1]);
            for (int u = t + 2; c && u <= nout; ++u) c = _addcarry_u64(c, r[u], 0, &r[u]);
        } else { // col == lo-1: only the hi limb (column lo) is kept
            unsigned char c = _addcarry_u64(0, r[0], hi, &r[0]);
            for (int u = 1; c && u <= nout; ++u) c = _addcarry_u64(c, r[u], 0, &r[u]);
        }
    }
}

static void study(int L, long iters) {
    std::vector<uint64_t> a(L), b(L), full(2 * L), hi(L + 4), hx(2 * L + 4);
    for (int i = 0; i < L; ++i) {
        a[i] = rnd64();
        b[i] = rnd64();
    }

    auto t0 = Clock::now();
    for (long k = 0; k < iters; ++k) {
        mpn_mul_n((mp_limb_t*)full.data(), (mp_limb_t*)a.data(), (mp_limb_t*)b.data(), L);
        a[0] ^= full[L];
    }
    auto t1 = Clock::now();
    double full_ns = ns_per(t0, t1, iters);

    // (2) mpn_sqr
    t0 = Clock::now();
    for (long k = 0; k < iters; ++k) {
        mpn_sqr((mp_limb_t*)full.data(), (mp_limb_t*)a.data(), L);
        a[0] ^= full[L];
    }
    t1 = Clock::now();
    double sqr_ns = ns_per(t0, t1, iters);

    // (3) schoolbook high-half
    t0 = Clock::now();
    for (long k = 0; k < iters; ++k) {
        mulhigh_school(hi.data(), a.data(), b.data(), L);
        a[0] ^= hi[2];
    }
    t1 = Clock::now();
    double hh_ns = ns_per(t0, t1, iters);

    // (3b) tuned ADX high-half (MULX + ADCX/ADOX), GUARD=2
    t0 = Clock::now();
    for (long k = 0; k < iters; ++k) {
        mulhigh_adx(hx.data(), a.data(), b.data(), L, 2);
        a[0] ^= hx[2];
    }
    t1 = Clock::now();
    double ax_ns = ns_per(t0, t1, iters);

    // (3c) high-half via GMP's own mpn_addmul_1 per row, GUARD=2
    t0 = Clock::now();
    for (long k = 0; k < iters; ++k) {
        mulhigh_addmul(hx.data(), a.data(), b.data(), L, 2);
        a[0] ^= hx[2];
    }
    t1 = Clock::now();
    double am_ns = ns_per(t0, t1, iters);

    // (4) mpf_mul at prec = 64*L bits (what the reference build actually uses)
    mpf_t ma, mb, mr;
    mp_bitcnt_t p = (mp_bitcnt_t)64 * L;
    mpf_init2(ma, p);
    mpf_init2(mb, p);
    mpf_init2(mr, p);
    mpf_set_ui(ma, 3);
    mpf_div_ui(ma, ma, 7);
    mpf_set_ui(mb, 5);
    mpf_div_ui(mb, mb, 11);
    t0 = Clock::now();
    for (long k = 0; k < iters; ++k) { mpf_mul(mr, ma, mb); }
    t1 = Clock::now();
    double mpf_ns = ns_per(t0, t1, iters);
    mpf_clear(ma);
    mpf_clear(mb);
    mpf_clear(mr);

    printf("  L=%2d (%5db): mpf %7.1f | mpn_mul_n %7.1f (%.2fx) | mpn_sqr %6.1f (%.2fx) | ADX_hi %7.1f (%.2fx) | addmul_hi %7.1f (%.2fx)\n", L, 64 * L, mpf_ns, full_ns, mpf_ns / full_ns, sqr_ns, mpf_ns / sqr_ns, ax_ns, mpf_ns / ax_ns, am_ns, mpf_ns / am_ns);
}

// Correctness: high-half top L limbs must match mpn_mul_n's top limbs exactly.
static void checkHigh(int L) {
    std::vector<uint64_t> a(L), b(L), full(2 * L), hi(L + 4), hx(2 * L + 4);
    for (int i = 0; i < L; ++i) {
        a[i] = rnd64();
        b[i] = rnd64();
    }
    mpn_mul_n((mp_limb_t*)full.data(), (mp_limb_t*)a.data(), (mp_limb_t*)b.data(), L);
    mulhigh_school(hi.data(), a.data(), b.data(), L);
    int badS = 0;
    {
        const int lo = L - 2;
        for (int col = lo; col < 2 * L; ++col)
            if (full[col] != hi[col - lo]) ++badS;
    }
    // ADX with GUARD=2: kept high L limbs are columns [L-1 .. 2L-2]; verify those.
    const int GUARD = 2;
    int lo = L - 1 - GUARD;
    if (lo < 0) lo = 0;
    mulhigh_adx(hx.data(), a.data(), b.data(), L, GUARD);
    int badA = 0;
    for (int col = L - 1; col < 2 * L; ++col)
        if (full[col] != hx[col - lo]) ++badA;
    std::vector<uint64_t> hm(2 * L + 4);
    mulhigh_addmul(hm.data(), a.data(), b.data(), L, GUARD);
    int badM = 0;
    for (int col = L - 1; col < 2 * L; ++col)
        if (full[col] != hm[col - lo]) ++badM;
    // high-half square vs mpn_sqr full
    std::vector<uint64_t> fsq(2 * L), hs(2 * L + 4);
    mpn_sqr((mp_limb_t*)fsq.data(), (const mp_limb_t*)a.data(), L);
    sqrhigh_addmul(hs.data(), a.data(), L, GUARD);
    int badQ = 0;
    for (int col = L - 1; col < 2 * L; ++col)
        if (fsq[col] != hs[col - lo]) ++badQ;
    printf("  correctness L=%2d: ADX_hi %s  addmul_hi %s  sqrhigh %s\n", L, badA == 0 ? "EXACT" : "approx", badM == 0 ? "EXACT" : "approx", badQ == 0 ? "EXACT" : "approx");
}

// Complex-square step cost: 2 high-half muls (engine's Karatsuba form) vs 3 high-half
// squares (a^2,b^2,(a+b)^2). Both need the same adds; we time just the multiplies.
static void complexStep(int L, long iters) {
    std::vector<uint64_t> a(L), b(L), s(L), r0(2 * L + 4), r1(2 * L + 4), r2(2 * L + 4);
    for (int i = 0; i < L; ++i) {
        a[i] = rnd64() >> 2;
        b[i] = rnd64() >> 2;
    }
    const int G = 3;
    auto t0 = Clock::now();
    for (long k = 0; k < iters; ++k) {
        for (int i = 0; i < L; ++i) s[i] = a[i] + b[i];      // a+b (approx; timing only)
        mulhigh_addmul(r0.data(), s.data(), a.data(), L, G); // stand-in for (a+b)(a-b)
        mulhigh_addmul(r1.data(), a.data(), b.data(), L, G); // a*b
        a[0] ^= r0[2] ^ r1[2];
    }
    auto t1 = Clock::now();
    double mul2 = ns_per(t0, t1, iters);
    t0 = Clock::now();
    for (long k = 0; k < iters; ++k) {
        for (int i = 0; i < L; ++i) s[i] = a[i] + b[i];
        sqrhigh_addmul(r0.data(), a.data(), L, G); // a^2
        sqrhigh_addmul(r1.data(), b.data(), L, G); // b^2
        sqrhigh_addmul(r2.data(), s.data(), L, G); // (a+b)^2
        a[0] ^= r0[2] ^ r1[2] ^ r2[2];
    }
    t1 = Clock::now();
    double sqr3 = ns_per(t0, t1, iters);
    printf("  L=%2d: 2*hh_mul %7.1f ns | 3*hh_sqr %7.1f ns | sqr-form %.2fx\n", L, mul2, sqr3, mul2 / sqr3);
}

int main() {
    srand(999);
    printf("== high-half correctness ==\n");
    for (int L : {6, 10, 16, 24, 48}) checkHigh(L);
    printf("\n== multiply variants (ns each; speedup vs mpf, higher=faster) ==\n");
    for (int L : {6, 8, 10, 12, 16, 20, 24, 32, 48, 56, 60, 64, 68, 72, 80, 88, 96, 128}) study(L, 1000000);
    printf("\n== complex-square: 2 high-half muls vs 3 high-half squares ==\n");
    for (int L : {16, 20, 24, 32, 48, 64}) complexStep(L, 1000000);
    return 0;
}
