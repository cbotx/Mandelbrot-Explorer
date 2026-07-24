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

// Full 2L-limb product of two L-limb magnitudes into r2 (schoolbook).
static inline void bf_mag_mul_full(uint64_t* r2, const uint64_t* a, const uint64_t* b, int L) {
    std::memset(r2, 0, sizeof(uint64_t) * (size_t)(2 * L));
    for (int i = 0; i < L; ++i) {
        uint64_t ai = a[i];
        if (!ai) continue;
        uint64_t carry = 0;
        for (int j = 0; j < L; ++j) {
            uint64_t hi;
            uint64_t lo = _umul128(ai, b[j], &hi);
            unsigned char c = _addcarry_u64(0, lo, carry, &lo);   // lo += carry
            hi += c;                                              // cannot overflow
            c = _addcarry_u64(0, r2[i + j], lo, &r2[i + j]);      // acc += lo
            carry = hi + c;                                       // hi + carry-out
        }
        r2[i + L] = carry;
    }
}

// r[L] = round( (a[L] * b[L]) >> 64*(L-1) ), scratch tmp[2L].
static inline void bf_mag_mulshift(uint64_t* r, const uint64_t* a, const uint64_t* b, int L, uint64_t* tmp) {
    bf_mag_mul_full(tmp, a, b, L);
    const int drop = L - 1;                     // limbs discarded from the bottom
    for (int i = 0; i < L; ++i) r[i] = tmp[drop + i];
    // Round to nearest using the most-significant discarded bit.
    if (drop > 0 && (tmp[drop - 1] >> 63)) {
        unsigned char c = _addcarry_u64(0, r[0], 1, &r[0]);
        for (int i = 1; i < L && c; ++i) c = _addcarry_u64(c, r[i], 0, &r[i]);
    }
}

// ---- BigFixed value ---------------------------------------------------------

struct BigFixed {
    int L = 0;
    int sign = 0;                 // +1, -1, or 0 for exact zero
    std::vector<uint64_t> m;      // L limbs, little-endian

    BigFixed() = default;
    explicit BigFixed(int L_) : L(L_), sign(0), m((size_t)L_, 0ull) {}

    void setL(int L_) { L = L_; sign = 0; m.assign((size_t)L_, 0ull); }
    bool isZero() const { return sign == 0; }

    void setZero() { sign = 0; std::memset(m.data(), 0, sizeof(uint64_t) * (size_t)L); }

    // Set from a small integer value v (|v| fits in the integer limb).
    void setInt(long long v) {
        setZero();
        if (v == 0) { sign = 0; return; }
        sign = v < 0 ? -1 : 1;
        unsigned long long a = v < 0 ? (unsigned long long)(-(v + 1)) + 1ull : (unsigned long long)v;
        m[L - 1] = a;             // integer part lives in the top limb
    }

    // value as a double (may under/overflow to 0/inf like the old path).
    double toDouble() const {
        if (sign == 0) return 0.0;
        double acc = 0.0;
        // Sum from the top limb; scale so the top limb is the integer part.
        for (int i = L - 1; i >= 0; --i) {
            acc = acc * 18446744073709551616.0 /* 2^64 */ + (double)m[i];
            if (acc > 1e300) break;   // lower limbs are negligible once large
        }
        // We summed m as an integer of L limbs; divide by 2^(64*(L-1)).
        acc = std::ldexp(acc, -64 * (L - 1));
        return sign < 0 ? -acc : acc;
    }
};

// r = a + b  (a, b, r may alias only if same object not required here).
inline void bf_add(BigFixed& r, const BigFixed& a, const BigFixed& b) {
    const int L = a.L;
    if (a.sign == 0) { r = b; return; }
    if (b.sign == 0) { r = a; return; }
    if (a.sign == b.sign) {
        r.L = L; r.m.resize(L); r.sign = a.sign;
        bf_mag_add(r.m.data(), a.m.data(), b.m.data(), L);   // overflow assumed absent (bounded values)
    } else {
        int c = bf_mag_cmp(a.m.data(), b.m.data(), L);
        r.L = L; r.m.resize(L);
        if (c == 0) { r.setZero(); return; }
        if (c > 0) { bf_mag_sub(r.m.data(), a.m.data(), b.m.data(), L); r.sign = a.sign; }
        else       { bf_mag_sub(r.m.data(), b.m.data(), a.m.data(), L); r.sign = b.sign; }
    }
}

inline void bf_sub(BigFixed& r, const BigFixed& a, const BigFixed& b) {
    BigFixed nb = b; nb.sign = -b.sign; bf_add(r, a, nb);
}

// r = a * b, using caller scratch tmp[2L] (avoids allocation in hot loops).
inline void bf_mul(BigFixed& r, const BigFixed& a, const BigFixed& b, uint64_t* tmp) {
    const int L = a.L;
    r.L = L; r.m.resize(L);
    if (a.sign == 0 || b.sign == 0) { r.setZero(); return; }
    bf_mag_mulshift(r.m.data(), a.m.data(), b.m.data(), L, tmp);
    r.sign = (a.sign == b.sign) ? 1 : -1;
    // A product that rounds to zero magnitude is exact zero.
    int nz = 0; for (int i = 0; i < L; ++i) nz |= (r.m[i] != 0);
    if (!nz) r.sign = 0;
}

// ---- GMP interop (for tests / A-B verification) -----------------------------

// Exact BigFixed -> mpf_t conversion.
inline void bf_to_mpf(mpf_t out, const BigFixed& a) {
    mpz_t z;
    mpz_init(z);
    mpz_import(z, (size_t)a.L, -1 /*order: least significant first*/, sizeof(uint64_t),
               0 /*native endian*/, 0, a.m.data());
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
    mpf_mul_2exp(t, t, (mp_bitcnt_t)64 * (L - 1));   // scale so integer part is the value
    mpz_t z;
    mpz_init(z);
    mpz_set_f(z, t);                                 // truncates fraction
    size_t count = 0;
    mpz_export(a.m.data(), &count, -1, sizeof(uint64_t), 0, 0, z);
    for (size_t i = count; i < (size_t)L; ++i) a.m[i] = 0;
    mpz_clear(z);
    mpf_clear(t);
}

#endif
