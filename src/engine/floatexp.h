#ifndef __FLOATEXP_H__
#define __FLOATEXP_H__

// Minimal floating-point-with-extended-exponent type for deep-zoom rescaling.
//
// A value is  m * 2^e  with the mantissa normalized to 0.5 <= |m| < 1 (or the
// exact zero {0,0}). This extends the exponent range far past IEEE double's
// ~1e-308 underflow, which is what caps naive double perturbation at ~1e320.
//
// It is used only for the *scale* S (and the rare full-range steps) in the
// rescaled iteration z = S w: the per-pixel delta w stays an O(1) double, so the
// hot loop keeps native double speed while S carries the deep exponent.

#include <cmath>
#include <cstdint>

struct FloatExp {
    double m;  // mantissa, 0.5 <= |m| < 1, or 0
    int64_t e; // exponent (base 2)
};

// ---- construction / conversion -------------------------------------------
static inline FloatExp fe_from(double x) {
    if (x == 0.0) return FloatExp{0.0, 0};
    int ex;
    double mm = std::frexp(x, &ex);
    return FloatExp{mm, (int64_t)ex};
}

// To double; underflows to 0 / overflows to +-inf just like the hardware would.
static inline double fe_to_double(FloatExp a) {
    if (a.m == 0.0) return 0.0;
    int64_t e = a.e;
    if (e > 1100) return a.m > 0 ? INFINITY : -INFINITY;
    if (e < -1100) return 0.0;
    return std::ldexp(a.m, (int)e);
}

static inline FloatExp fe_renorm(double m, int64_t e) {
    if (m == 0.0) return FloatExp{0.0, 0};
    int ex;
    double mm = std::frexp(m, &ex);
    return FloatExp{mm, e + ex};
}

// ---- arithmetic ----------------------------------------------------------
static inline FloatExp fe_mul(FloatExp a, FloatExp b) {
    if (a.m == 0.0 || b.m == 0.0) return FloatExp{0.0, 0};
    // |m*m| in [0.25, 1): one frexp normalizes it.
    return fe_renorm(a.m * b.m, a.e + b.e);
}

static inline FloatExp fe_mul_d(FloatExp a, double b) {
    if (a.m == 0.0 || b == 0.0) return FloatExp{0.0, 0};
    return fe_renorm(a.m * b, a.e);
}

// Multiply by 2^k (exact, no rounding).
static inline FloatExp fe_scale2(FloatExp a, int64_t k) {
    return a.m == 0.0 ? a : FloatExp{a.m, a.e + k};
}

static inline FloatExp fe_add(FloatExp a, FloatExp b) {
    if (a.m == 0.0) return b;
    if (b.m == 0.0) return a;
    int64_t d = a.e - b.e;
    if (d > 54) return a;  // b negligible vs a
    if (d < -54) return b; // a negligible vs b
    if (d >= 0)
        return fe_renorm(a.m + std::ldexp(b.m, -(int)d), a.e);
    else
        return fe_renorm(b.m + std::ldexp(a.m, (int)d), b.e);
}

static inline FloatExp fe_neg(FloatExp a) {
    return FloatExp{-a.m, a.e};
}
static inline FloatExp fe_sub(FloatExp a, FloatExp b) {
    return fe_add(a, fe_neg(b));
}

static inline FloatExp fe_div(FloatExp a, FloatExp b) {
    if (a.m == 0.0) return FloatExp{0.0, 0};
    return fe_renorm(a.m / b.m, a.e - b.e);
}

// ---- comparison ----------------------------------------------------------
// |a| < |b|
static inline bool fe_abs_less(FloatExp a, FloatExp b) {
    if (a.m == 0.0) return b.m != 0.0;
    if (b.m == 0.0) return false;
    if (a.e != b.e) return a.e < b.e;
    return std::fabs(a.m) < std::fabs(b.m);
}

// sqrt of a non-negative FloatExp magnitude value.
static inline FloatExp fe_sqrt(FloatExp a) {
    if (a.m <= 0.0) return FloatExp{0.0, 0};
    // value = m * 2^e; split exponent parity so the residual exponent is even.
    int64_t e = a.e;
    double m = a.m;
    if (e & 1) {
        m *= 2.0;
        --e;
    } // now e even, m in [1,2)
    return fe_renorm(std::sqrt(m), e / 2);
}

#endif
