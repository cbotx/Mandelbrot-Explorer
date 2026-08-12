#ifndef __BIGFIXED_H__
#define __BIGFIXED_H__

// BigFixed: fixed-precision, fixed-point big real tailored to Mandelbrot
// perturbation. A value is stored as
//     value = sign * ( sum_{i=0}^{L-1} m[i] * 2^(64*i) ) * 2^(-64*(L-1))
// i.e. little-endian limbs with the binary point just below the top limb, so
// m[L-1] is the integer part and m[L-2..0] are the fraction. Reference-orbit
// values (z, c) are O(1) in magnitude, so a single integer limb is plenty; the
// remaining L-1 limbs carry ~64*(L-1) fractional bits of precision.
//
// The point is fixed for every value of a given L, so add/sub are plain
// limbwise operations (no exponent bookkeeping) and multiply is a fixed-size
// schoolbook product followed by a fixed right shift -- far less overhead than
// the generic mpf_t. Uses MSVC x64 intrinsics.

#include <cstdint>
#include <cstring>
#include <vector>
#include <intrin.h>
#include <gmp.h>
#include "floatexp.h"

// ---- raw magnitude helpers (operate on L-limb little-endian arrays) --------

// Compare magnitudes: -1 if a<b, 0 if a==b, +1 if a>b.
static inline int bf_mag_cmp(const uint64_t* a, const uint64_t* b, int L) {
    for (int i = L - 1; i >= 0; --i)
        if (a[i] != b[i]) return a[i] < b[i] ? -1 : 1;
    return 0;
}

// r = a + b (same length). Returns carry out of the top limb.
static inline unsigned char bf_mag_add(uint64_t* r, const uint64_t* a, const uint64_t* b, int L) {
    unsigned char c = 0;
    for (int i = 0; i < L; ++i) c = _addcarry_u64(c, a[i], b[i], &r[i]);
    return c;
}

// r = a - b, assuming a >= b. Returns borrow (should be 0 when a>=b).
static inline unsigned char bf_mag_sub(uint64_t* r, const uint64_t* a, const uint64_t* b, int L) {
    unsigned char br = 0;
    for (int i = 0; i < L; ++i) br = _subborrow_u64(br, a[i], b[i], &r[i]);
    return br;
}

// Full 2L-limb product of two L-limb magnitudes into r2. Uses GMP's tuned mpn
// limb-multiply (asm + Karatsuba/Toom at larger L), which is far faster than a
// hand-rolled schoolbook loop above ~a few limbs; the schoolbook version is kept
// below (bf_mag_mul_full_school) only as a reference/fallback. mpn_mul_n needs
// r2 disjoint from a and b (callers pass separate scratch), and mp_limb_t is a
// 64-bit unsigned on this target so the uint64_t limbs map directly.
static inline void bf_mag_mul_full_school(uint64_t* r2, const uint64_t* a, const uint64_t* b, int L) {
    std::memset(r2, 0, sizeof(uint64_t) * (size_t)(2 * L));
    for (int i = 0; i < L; ++i) {
        uint64_t ai = a[i];
        if (!ai) continue;
        uint64_t carry = 0;
        for (int j = 0; j < L; ++j) {
            uint64_t hi;
            uint64_t lo = _umul128(ai, b[j], &hi);
            unsigned char c = _addcarry_u64(0, lo, carry, &lo); // lo += carry
            hi += c;                                            // cannot overflow
            c = _addcarry_u64(0, r2[i + j], lo, &r2[i + j]);    // acc += lo
            carry = hi + c;                                     // hi + carry-out
        }
        r2[i + L] = carry;
    }
}

static inline void bf_mag_mul_full(uint64_t* r2, const uint64_t* a, const uint64_t* b, int L) {
    mpn_mul_n(reinterpret_cast<mp_limb_t*>(r2), reinterpret_cast<const mp_limb_t*>(a), reinterpret_cast<const mp_limb_t*>(b), (mp_size_t)L);
}

// HIGH-half product: fills hbuf[t] = product column (lo+t), lo = L-1-GUARD, for
// columns [lo .. 2L-1] (nout = 2L-lo limbs, +1 carry headroom). Built from GMP's
// own tuned asm primitive mpn_addmul_1 applied per row but only over the kept
// high columns -- every partial product with i+j >= lo is included, so the result
// limbs (columns >= L-1) are EXACT; only the dropped low half (< lo) is skipped,
// which is ~half the work. Beats a full mpn_mul_n for L >= ~16 (below that the
// per-row call overhead dominates). Requires L >= GUARD+2.
static inline void bf_mag_mulhigh(uint64_t* hbuf, const uint64_t* a, const uint64_t* b, int L, int GUARD) {
    const int lo = L - 1 - GUARD;
    const int nout = 2 * L - lo;
    std::memset(hbuf, 0, sizeof(uint64_t) * (size_t)(nout + 1));
    for (int i = 0; i < L; ++i) {
        int jstart = lo - i;
        if (jstart < 0) jstart = 0;
        const int n = L - jstart;
        const int t = i + jstart - lo;
        mp_limb_t carry = mpn_addmul_1(reinterpret_cast<mp_limb_t*>(hbuf + t), reinterpret_cast<const mp_limb_t*>(b + jstart), (mp_size_t)n, (mp_limb_t)a[i]);
        int u = t + n;
        mp_limb_t cc = carry;
        while (cc && u <= nout) {
            unsigned char c = _addcarry_u64(0, hbuf[u], cc, &hbuf[u]);
            cc = c;
            ++u;
        }
    }
}

// r[L] = round( (a[L] * b[L]) >> 64*(L-1) ), scratch tmp[2L]. For deep references
// (L large) this only computes the needed HIGH half of the product (mpf/mpn_mul_n
// waste ~half computing the discarded low limbs); for small L the full product is
// faster.
static inline void bf_mag_mulshift(uint64_t* r, const uint64_t* a, const uint64_t* b, int L, uint64_t* tmp) {
    constexpr int GUARD = 3;
    if (L >= 16 && L < 72) {
        bf_mag_mulhigh(tmp, a, b, L, GUARD);               // tmp[t] = column (L-1-GUARD+t)
        for (int i = 0; i < L; ++i) r[i] = tmp[i + GUARD]; // result limb i = column L-1+i
        if (tmp[GUARD - 1] >> 63) {                        // round-nearest from column L-2
            unsigned char c = _addcarry_u64(0, r[0], 1, &r[0]);
            for (int i = 1; i < L && c; ++i) c = _addcarry_u64(c, r[i], 0, &r[i]);
        }
        return;
    }
    bf_mag_mul_full(tmp, a, b, L);
    const int drop = L - 1; // limbs discarded from the bottom
    for (int i = 0; i < L; ++i) r[i] = tmp[drop + i];
    // Round to nearest using the most-significant discarded bit.
    if (drop > 0 && (tmp[drop - 1] >> 63)) {
        unsigned char c = _addcarry_u64(0, r[0], 1, &r[0]);
        for (int i = 1; i < L && c; ++i) c = _addcarry_u64(c, r[i], 0, &r[i]);
    }
}

// r[L] = round( (a[L]^2) >> 64*(L-1) ), scratch tmp[2L]. Uses GMP's dedicated
// squaring (mpn_sqr), which is ~1.5-1.85x faster than a general multiply, so a
// complex square done as three squares (a^2, b^2, (a+b)^2) can beat the two
// general multiplies of the Karatsuba complex-square form at deep precision.
static inline void bf_mag_sqrshift(uint64_t* r, const uint64_t* a, int L, uint64_t* tmp) {
    mpn_sqr(reinterpret_cast<mp_limb_t*>(tmp), reinterpret_cast<const mp_limb_t*>(a), (mp_size_t)L);
    const int drop = L - 1;
    for (int i = 0; i < L; ++i) r[i] = tmp[drop + i];
    if (drop > 0 && (tmp[drop - 1] >> 63)) {
        unsigned char c = _addcarry_u64(0, r[0], 1, &r[0]);
        for (int i = 1; i < L && c; ++i) c = _addcarry_u64(c, r[i], 0, &r[i]);
    }
}

// ---- BigFixed value ---------------------------------------------------------

struct BigFixed {
    int L = 0;
    int sign = 0;            // +1, -1, or 0 for exact zero
    std::vector<uint64_t> m; // L limbs, little-endian

    BigFixed() = default;
    explicit BigFixed(int L_) : L(L_), sign(0), m((size_t)L_, 0ull) {}

    void setL(int L_) {
        L = L_;
        sign = 0;
        m.assign((size_t)L_, 0ull);
    }
    bool isZero() const { return sign == 0; }

    void setZero() {
        sign = 0;
        std::memset(m.data(), 0, sizeof(uint64_t) * (size_t)L);
    }

    // Set from a small integer value v (|v| fits in the integer limb).
    void setInt(long long v) {
        setZero();
        if (v == 0) {
            sign = 0;
            return;
        }
        sign = v < 0 ? -1 : 1;
        unsigned long long a = v < 0 ? (unsigned long long)(-(v + 1)) + 1ull : (unsigned long long)v;
        m[L - 1] = a; // integer part lives in the top limb
    }

    // value as a double (may under/overflow to 0/inf like the old path). Uses the
    // top two nonzero limbs (>=53 bits) and scales by the correct exponent, so it is
    // correct for any L and any magnitude (the old full-Horner form broke for large
    // L: its overflow-guard early-break did not adjust the final scale).
    double toDouble() const {
        if (sign == 0) return 0.0;
        int top = L - 1;
        while (top >= 0 && m[top] == 0) --top;
        if (top < 0) return 0.0;
        double hi = (double)m[top];
        double lo = (top >= 1) ? (double)m[top - 1] : 0.0;
        double val = hi * 18446744073709551616.0 /* 2^64 */ + lo; // place 2^(64*(top-1))
        val = std::ldexp(val, 64 * (top - (L - 1) - 1));          // * 2^(64*(top-L))
        return sign < 0 ? -val : val;
    }
};

// r = signed(a) + signed(b) where the operand signs are passed explicitly, so
// subtraction needs no temporary copy (bf_sub just flips sb). Values are bounded
// (|z|,|c| = O(1)), so magnitude add never overflows the integer limb.
inline void bf_addsub(BigFixed& r, const BigFixed& a, int sa, const BigFixed& b, int sb) {
    const int L = a.L;
    if (sa == 0) {
        r.L = L;
        if ((int)r.m.size() != L) r.m.resize(L);
        std::memcpy(r.m.data(), b.m.data(), sizeof(uint64_t) * (size_t)L);
        r.sign = sb;
        return;
    }
    if (sb == 0) {
        r.L = L;
        if ((int)r.m.size() != L) r.m.resize(L);
        std::memcpy(r.m.data(), a.m.data(), sizeof(uint64_t) * (size_t)L);
        r.sign = sa;
        return;
    }
    r.L = L;
    if ((int)r.m.size() != L) r.m.resize(L);
    if (sa == sb) {
        bf_mag_add(r.m.data(), a.m.data(), b.m.data(), L); // bounded: no overflow
        r.sign = sa;
    } else {
        int c = bf_mag_cmp(a.m.data(), b.m.data(), L);
        if (c == 0) {
            r.setZero();
            return;
        }
        if (c > 0) {
            bf_mag_sub(r.m.data(), a.m.data(), b.m.data(), L);
            r.sign = sa;
        } else {
            bf_mag_sub(r.m.data(), b.m.data(), a.m.data(), L);
            r.sign = sb;
        }
    }
}

// r = a + b  (a, b, r may alias only if same object not required here).
inline void bf_add(BigFixed& r, const BigFixed& a, const BigFixed& b) {
    bf_addsub(r, a, a.sign, b, b.sign);
}

inline void bf_sub(BigFixed& r, const BigFixed& a, const BigFixed& b) {
    bf_addsub(r, a, a.sign, b, -b.sign); // no copy of b
}

// r = a * b, using caller scratch tmp[2L] (avoids allocation in hot loops).
inline void bf_mul(BigFixed& r, const BigFixed& a, const BigFixed& b, uint64_t* tmp) {
    const int L = a.L;
    r.L = L;
    if ((int)r.m.size() != L) r.m.resize(L);
    if (a.sign == 0 || b.sign == 0) {
        r.setZero();
        return;
    }
    bf_mag_mulshift(r.m.data(), a.m.data(), b.m.data(), L, tmp);
    // A rounded-to-zero magnitude keeps a nonzero sign, which is harmless: every
    // consumer (toDouble/bf_to_fe/bf_to_mpf) treats a zero magnitude as value 0.
    r.sign = (a.sign == b.sign) ? 1 : -1;
}

// r = a*a (always >= 0), using GMP's faster dedicated squaring. Scratch tmp[2L].
inline void bf_sqr(BigFixed& r, const BigFixed& a, uint64_t* tmp) {
    const int L = a.L;
    r.L = L;
    if ((int)r.m.size() != L) r.m.resize(L);
    if (a.sign == 0) {
        r.setZero();
        return;
    }
    bf_mag_sqrshift(r.m.data(), a.m.data(), L, tmp);
    r.sign = 1;
}

// ---- GMP interop (for tests / A-B verification) -----------------------------

// Exact BigFixed -> mpf_t conversion.
inline void bf_to_mpf(mpf_t out, const BigFixed& a) {
    mpz_t z;
    mpz_init(z);
    mpz_import(z, (size_t)a.L, -1 /*order: least significant first*/, sizeof(uint64_t), 0 /*native endian*/, 0, a.m.data());
    mpf_set_z(out, z);
    mpf_div_2exp(out, out, (mp_bitcnt_t)64 * (a.L - 1));
    if (a.sign < 0) mpf_neg(out, out);
    mpz_clear(z);
}

// mpf_t -> BigFixed (round toward zero at the fixed-point resolution).
inline void bf_from_mpf(BigFixed& a, const mpf_t in, int L) {
    a.setL(L);
    if (mpf_sgn(in) == 0) return;
    a.sign = mpf_sgn(in);
    mpf_t t;
    mpf_init2(t, (mp_bitcnt_t)64 * L + 8);
    mpf_abs(t, in);
    mpf_mul_2exp(t, t, (mp_bitcnt_t)64 * (L - 1)); // scale so integer part is the value
    mpz_t z;
    mpz_init(z);
    mpz_set_f(z, t); // truncates fraction
    size_t count = 0;
    mpz_export(a.m.data(), &count, -1, sizeof(uint64_t), 0, 0, z);
    for (size_t i = count; i < (size_t)L; ++i) a.m[i] = 0;
    mpz_clear(z);
    mpf_clear(t);
}

// BigFixed -> FloatExp (m * 2^e, 0.5 <= |m| < 1), relative precision from the top
// nonzero limbs. Value = mant * 2^(-64*(L-1)); build the 53-bit mantissa from the
// two highest nonzero limbs so tiny |z| (leading-zero limbs) keep full relative
// precision -- this is what fixed-point needs to feed the deep floatexp shadow.
static inline FloatExp bf_to_fe(const BigFixed& a) {
    if (a.sign == 0) return FloatExp{0.0, 0};
    int top = a.L - 1;
    while (top >= 0 && a.m[top] == 0) --top;
    if (top < 0) return FloatExp{0.0, 0};
    // 128-bit window from limbs top, top-1 -> a double mantissa.
    unsigned long long hi = a.m[top];
    unsigned long long lo = (top >= 1) ? a.m[top - 1] : 0ull;
    double mant = (double)hi * 18446744073709551616.0 + (double)lo; // hi*2^64 + lo
    // (hi*2^64+lo) sits at place 2^(64*(top-1)); value = mant * 2^(64*(top-1)) *
    // 2^(-64*(L-1)) = mant * 2^(64*(top-L)).
    int64_t e = (int64_t)64 * (top - a.L);
    FloatExp r = fe_from(mant); // normalises mant into [0.5,1)*2^k
    r.e += e;
    return a.sign < 0 ? fe_neg(r) : r;
}

#endif
