#include <iostream>
#include <iomanip>
#include <complex>
#include <gmp.h>
#include <assert.h>
#include <chrono>

#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <array>
#include <utility>
#include <cstring>
#include <cstdlib>
#include <immintrin.h>
#include <omp.h>

#include "mandel_perturbation.h"
#include "float_math.h"
#include "dll_interface.h"

// BLA profiling: per-thread padded counters (no false sharing / contention).
// [tid][0]=iterations skipped, [1]=BLA applies, [2]=normal steps.
static long long g_bla_stat[64][8];
static int g_bla_noescape = -1;   // env MANDEL_BLA_NOESCAPE: force-skip the escape check (timing only)
static long long g_fe_fallback = 0;
static long long g_fe_stat[64][3];   // [tid] 0=BLA skip-iters 1=BLA applies 2=normal steps (MANDEL_PROFILE)
static long long g_blafe_safe[64] = { 0 };   // [tid] tryBLAfe overflow-safe floatexp fallbacks
static int g_int_rep = -2;       // env MANDEL_INT_REP: exact state-repetition interior detection (default on; -2=unread)
static int g_bigfixed = -2;      // env MANDEL_BIGFIXED: -2 unparsed, -1 auto (on for floatexp), 0 off, >=1 force on

// Monotonic seconds, for the MANDEL_PROFILE phase timers.
static inline double nowSec() {
    return std::chrono::duration_cast<std::chrono::duration<double>>(
        std::chrono::high_resolution_clock::now().time_since_epoch()).count();
}

// Stripe Average Coloring window. Averaging the stripe term over the whole orbit
// washes out at deep zoom (orbits are long and nearly identical, so Feather goes
// flat); the per-pixel structure lives in the escaping tail. Weighting recent
// iterations recovers it, and since it is a pure function of the orbit it stays
// pan-invariant (no reference halo). A rectangular "last W" window steps sharply
// when the escape count shifts a feature past its hard far edge (a visible
// brightness line); an exponential window has infinite support so a BLA-skip
// reset leaks a decayed tail (a faint halo). A triangular (Bartlett) window fixes
// both: its weight fades linearly to 0 at age W, so there is no far-edge cliff
// AND its support is finite (age>=W has weight 0), so it stays reference-clean.
// MANDEL_SACWIN<=0 restores the classic full average.
static int g_sac_win = -2;   // -2 = unread
static inline int sacWindow() {
    if (g_sac_win == -2) { const char* e = getenv("MANDEL_SACWIN"); g_sac_win = e ? atoi(e) : 256; }
    return g_sac_win;
}
namespace {
struct SacAccum {
    static const int MAXW = 1024;
    double ring[MAXW];
    double RS = 0.0, TS = 0.0, full = 0.0, last = 0.0;   // rect sum, triangular-weighted sum
    int W = 0, pos = 0, fill = 0, cnt = 0;
    void init(int w) { W = w > MAXW ? MAXW : (w < 0 ? 0 : w); RS = TS = full = last = 0.0; pos = fill = cnt = 0; }
    inline void push(double zr, double zi) {
        double t = 0.5 + 0.5 * sin(7.0 * atan2(zi, zr));
        full += t; last = t; ++cnt;
        if (W > 0) {
            double evicted = (fill == W) ? ring[pos] : (++fill, 0.0);
            TS += (double)W * t - RS;         // newest gets weight W; every other weight -1
            RS += t - evicted;
            ring[pos] = t; pos = (pos + 1) % W;
        }
    }
    inline void reset_window() { RS = TS = 0.0; pos = fill = 0; }   // a BLA skip breaks the tail
    inline void add_full(double stripeSum, int n) { full += stripeSum; cnt += n; }  // W==0 skip restore
    inline float value(double zrad, double R) const {
        double frac = 1.0 - log(log(zrad) * 0.5 / log(R)) / log(2.0);
        frac -= floor(frac);
        double S, WS;   // weighted stripe sum and total weight
        if (W > 0) { S = TS; WS = (double)fill * W - (double)fill * (fill - 1) / 2.0; }   // triangular weights
        else       { S = full; WS = cnt; }
        double a1 = WS > 0 ? S / WS : 0.0;
        double lw = W > 0 ? (double)W : 1.0;                       // the last term's weight
        double a2 = WS > lw ? (S - lw * last) / (WS - lw) : a1;    // average without the newest term
        return (float)(a2 + (a1 - a2) * frac);
    }
};
// Orbit-trap accumulator: tracks the orbit's closest approach to a composite trap
// (a point at 0, the Pickover cross x=0/y=0, and the ring |z|=0.5) plus the angle
// at the point-trap minimum, and maps them to a palette coordinate. Deterministic
// (no randomness); the trap shape is the design parameter.
struct TrapAccum {
    double minPoint = 1e300, trapAngle = 0.0, minCross = 1e300, minCircle = 1e300;
    inline void push(double zr, double zi) {
        double m2 = zr * zr + zi * zi;
        if (m2 < minPoint) { minPoint = m2; trapAngle = atan2(zi, zr); }
        double cross = std::min(std::fabs(zr), std::fabs(zi));
        if (cross < minCross) minCross = cross;
        double circle = std::fabs(std::sqrt(m2) - 0.5);
        if (circle < minCircle) minCircle = circle;
    }
    inline float value(double mu) const {
        double d = std::min(std::min(std::sqrt(minPoint), minCross * 1.5), minCircle);
        // Clamp so -log10(d) >= 0: beyond ~1 the orbit never came near a trap, so there
        // is no "hit" (trap = 0). This ALSO keeps the returned value non-negative --
        // a negative colour value would collide with the baseU<0 interior sentinel in
        // the AA re-colour path and paint far-field pixels solid black.
        if (d > 1.0) d = 1.0;
        double trap = -std::log10(std::max(d, 1e-14));      // >= 0
        if (mu < 0) mu = 0;                                 // far-field fast escapes
        // NB: no atan2(trapAngle) term -- a wrapping angle added to a linear palette
        // coordinate seams at the +/-pi branch cut (a visible left/right divide when the
        // closest-point direction crosses the negative real axis). Distance + smooth
        // count already give a rich, seamless trap.
        return (float)(0.17 * trap + 0.025 * mu);           // always >= 0, seamless
    }
};
}

// GMP mpf -> floatexp (m * 2^e, 0.5 <= |m| < 1), no underflow.
static inline FloatExp mpf_to_fe(mpf_srcptr x) {
    if (mpf_sgn(x) == 0) return FloatExp{ 0.0, 0 };
    long ex; double d = mpf_get_d_2exp(&ex, x);
    return FloatExp{ d, (int64_t)ex };
}


// ---------------------------------------------------------------------------
// GMP minibrot nucleus finder for the periodic-reference experiment. Ball
// (interval) arithmetic finds the dominant period in a disk; Newton refines the
// critical orbit to return to 0 (the nucleus). All in GMP so it works at any
// depth (radii below double's underflow). Ported from src/bench/periodic_ref_bench.
// ---------------------------------------------------------------------------
namespace {
struct GScratch {
    mpf_t t1, t2, t3;
    explicit GScratch(mp_bitcnt_t p) { mpf_init2(t1, p); mpf_init2(t2, p); mpf_init2(t3, p); }
    ~GScratch() { mpf_clear(t1); mpf_clear(t2); mpf_clear(t3); }
};
// z = z^2 + c (Karatsuba square: 2 muls)
inline void g_sqadd(mpf_t zr, mpf_t zi, mpf_srcptr cr, mpf_srcptr ci, GScratch& s) {
    mpf_add(s.t1, zr, zi); mpf_sub(s.t2, zr, zi); mpf_mul(s.t3, zr, zi);
    mpf_mul(zr, s.t1, s.t2); mpf_add(zr, zr, cr);
    mpf_mul_ui(zi, s.t3, 2); mpf_add(zi, zi, ci);
}
inline void g_abs(mpf_t out, mpf_srcptr zr, mpf_srcptr zi, GScratch& s) {
    mpf_mul(s.t1, zr, zr); mpf_mul(s.t2, zi, zi); mpf_add(s.t1, s.t1, s.t2); mpf_sqrt(out, s.t1);
}
// Ball-arithmetic period: propagate the disk of c-values (centre, radius r) under
// z->z^2+c as a disk (centre z_n, radius R_n): R' = (2|z|+R)R + r (+ a small
// precision-scaled roundoff allowance). First p with |z_p| <= R_p => a period-p
// nucleus lies in the disk. Returns period (0 if none / escaped).
static int g_ballPeriod(mpf_srcptr cr, mpf_srcptr ci, mpf_srcptr r, int maxp, mp_bitcnt_t prec) {
    GScratch s(prec);
    mpf_t zr, zi, R, R2m, zabs2, nb, ru, t;
    mpf_init2(zr, prec); mpf_init2(zi, prec); mpf_init2(R, prec); mpf_init2(R2m, prec);
    mpf_init2(zabs2, prec); mpf_init2(nb, prec); mpf_init2(ru, prec); mpf_init2(t, prec);
    mpf_set_ui(zr, 0); mpf_set_ui(zi, 0); mpf_set_ui(R, 0);
    const double cab = std::hypot(mpf_get_d(cr), mpf_get_d(ci));   // |c| ~ O(1)
    mpf_set_ui(ru, 1); mpf_div_2exp(ru, ru, prec > 64 ? prec - 32 : prec / 2);
    int period = 0;
    for (int p = 1; p <= maxp; ++p) {
        // Orbit magnitudes are O(1), so |z| for the radius recurrence is taken in
        // double (no mpf_sqrt): R' = (2|z| + R) R + r + roundoff. R stays in GMP
        // because r (= 2/scale) can be far below double's underflow.
        double za = std::hypot(mpf_get_d(zr), mpf_get_d(zi));
        mpf_set_d(t, 2.0 * za); mpf_add(nb, t, R); mpf_mul(nb, nb, R); mpf_add(nb, nb, r);
        mpf_set_d(t, za * za + cab + 1.0); mpf_mul(t, t, ru); mpf_add(nb, nb, t);
        mpf_set(R, nb);
        g_sqadd(zr, zi, cr, ci, s);
        // Exact disk-contains-0 test |z|^2 <= R^2 (GMP, no sqrt).
        mpf_mul(s.t1, zr, zr); mpf_mul(s.t2, zi, zi); mpf_add(zabs2, s.t1, s.t2);
        mpf_mul(R2m, R, R);
        if (mpf_cmp(zabs2, R2m) <= 0) { period = p; break; }
        // Escape: the whole disk left the bailout region. R, r are tiny; |c|+2 dominates.
        double zb = std::hypot(mpf_get_d(zr), mpf_get_d(zi));
        if (zb > cab + 2.0 + mpf_get_d(R)) break;
    }
    mpf_clears(zr, zi, R, R2m, zabs2, nb, ru, t, (mpf_ptr)0);
    return period;
}
// Newton-refine (cr,ci) to the period-p nucleus: c <- c - z_p / (dz_p/dc).
static bool g_newton(mpf_t cr, mpf_t ci, int period, mp_bitcnt_t prec, int maxit) {
    GScratch s(prec);
    mpf_t zr, zi, dr, di, den, nr, ni, tol, step2;
    mpf_init2(zr, prec); mpf_init2(zi, prec); mpf_init2(dr, prec); mpf_init2(di, prec);
    mpf_init2(den, prec); mpf_init2(nr, prec); mpf_init2(ni, prec); mpf_init2(tol, prec); mpf_init2(step2, prec);
    mpf_set_ui(tol, 1); mpf_div_2exp(tol, tol, prec > 96 ? prec - 48 : prec / 2); mpf_mul(tol, tol, tol);
    bool conv = false;
    for (int it = 0; it < maxit; ++it) {
        mpf_set_ui(zr, 0); mpf_set_ui(zi, 0); mpf_set_ui(dr, 0); mpf_set_ui(di, 0);
        for (int i = 0; i < period; ++i) {                      // z_p, dz_p from z0=0
            mpf_mul(s.t1, zr, dr); mpf_mul(s.t2, zi, di); mpf_sub(s.t3, s.t1, s.t2);
            mpf_mul_ui(s.t3, s.t3, 2); mpf_add_ui(s.t3, s.t3, 1);
            mpf_mul(s.t1, zr, di); mpf_mul(s.t2, zi, dr); mpf_add(di, s.t1, s.t2); mpf_mul_ui(di, di, 2);
            mpf_set(dr, s.t3);
            g_sqadd(zr, zi, cr, ci, s);
        }
        mpf_mul(s.t1, dr, dr); mpf_mul(s.t2, di, di); mpf_add(den, s.t1, s.t2);
        if (mpf_sgn(den) == 0) break;
        mpf_mul(s.t1, zr, dr); mpf_mul(s.t2, zi, di); mpf_add(s.t3, s.t1, s.t2); mpf_div(nr, s.t3, den);
        mpf_mul(s.t1, zi, dr); mpf_mul(s.t2, zr, di); mpf_sub(s.t3, s.t1, s.t2); mpf_div(ni, s.t3, den);
        mpf_sub(cr, cr, nr); mpf_sub(ci, ci, ni);
        mpf_mul(s.t1, nr, nr); mpf_mul(s.t2, ni, ni); mpf_add(step2, s.t1, s.t2);
        if (mpf_cmp(step2, tol) <= 0) { conv = true; break; }
    }
    mpf_clears(zr, zi, dr, di, den, nr, ni, tol, step2, (mpf_ptr)0);
    return conv;
}
} // namespace

Mandel::Mandel(int width, int height, int max_iteration, int sub, float* iter) : _w(width), _h(height), _mxit(max_iteration), _sub(sub), _iter(iter) {
    assert(width > 0);
    assert(height > 0);
    assert(max_iteration > 0);
    assert(sub % 2);
    { const char* e = getenv("MANDEL_BAILOUT"); if (e && atof(e) > 0) _ESCAPE_RADIUS = (float)atof(e); }
    // The full-precision reference/derivative orbits are only ever read at the
    // current and previous index during the build (the delta loop uses the double
    // and floatexp shadows below), so two rotating buffers replace mxit+1 mpf_t
    // each -- saving the whole serial O(mxit) alloc/init and hundreds of MB.
    _z_re = new mpf_t[2];
    _z_im = new mpf_t[2];
    _zf = new Comp[_mxit + 1];
    _zfr = new Float[_mxit + 1];
    _zfi = new Float[_mxit + 1];
    _zfr_fe = new FloatExp[_mxit + 1];
    _zfi_fe = new FloatExp[_mxit + 1];

    _df = new Comp[_mxit + 1];
    _dfr = new Float[_mxit + 1];
    _dfi = new Float[_mxit + 1];

    _done = new bool[width * height * sub * sub];
    _z_m3 = new Float[_mxit + 1];

    mpf_init(_c0_re);
    mpf_init(_c0_im);
    mpf_init(_dx);
    mpf_init(_dy);
    mpf_init(_scale);
    mpf_init(_ref_z_re);
    mpf_init(_ref_z_im);

    mpf_init(_t1);
    mpf_init(_t2);
    for (int i = 0; i < 2; ++i) {
        mpf_init(_z_re[i]);
        mpf_init(_z_im[i]);
    }
}

Mandel::~Mandel() {
    delete[] _done;
    delete[] _z_m3;

    mpf_clear(_c0_re);
    mpf_clear(_c0_im);
    mpf_clear(_dx);
    mpf_clear(_dy);
    mpf_clear(_scale);
    mpf_clear(_ref_z_re);
    mpf_clear(_ref_z_im);

    mpf_clear(_t1);
    mpf_clear(_t2);
    for (int i = 0; i < 2; ++i) {
        mpf_clear(_z_re[i]);
        mpf_clear(_z_im[i]);
    }
    delete[] _zf;
    delete[] _z_re;
    delete[] _z_im;

    delete[] _df;
    delete[] _zfr_fe;
    delete[] _zfi_fe;
}

bool Mandel::attractor(double zz_re, double zz_im, const double c_re, const double c_im, int period) const {
    const double epsilon = 1e-20;
    double z_re, z_im;
    double dz_re, dz_im;
    double zz1_re, zz1_im;
    double tmp, t2;
    for (int j = 0; j < 64; ++j) {
        z_re = zz_re;
        z_im = zz_im;
        dz_re = 1;
        dz_im = 0;
        for (int i = 0; i < period; ++i) {
            tmp = (z_re * dz_re - z_im * dz_im) * 2;
            dz_im = (z_re * dz_im + z_im * dz_re) * 2;
            dz_re = tmp;

            tmp = (z_re * z_re - z_im * z_im) + c_re;
            z_im = (z_re * z_im) * 2 + c_im;
            z_re = tmp;
        }
        z_re -= zz_re;
        z_im -= zz_im;
        tmp = dz_re - 1;
        t2 = tmp * tmp + dz_im * dz_im;
        zz1_re = (z_re * tmp + z_im * dz_im) / t2;
        zz1_im = (z_im * tmp - z_re * dz_im) / t2;
        if (zz1_re * zz1_re + zz1_im * zz1_im < epsilon) {
            return (dz_re * dz_re + dz_im * dz_im <= 1);
        }
        zz_re -= zz1_re;
        zz_im -= zz1_im;
    }
    return false;
}


double Mandel::floatPointCompute(Float c_re, Float c_im, int mxit, int c_method, double* normalOut) const {
    Float z_re = c_re;
    Float z_im = c_im;
    Float d_re = 2;
    Float d_im = 0;
    Float dc_re = 1;
    Float dc_im = 0;
    Float tmp;
    const bool sac = (c_method & ColoringMethod::STRIPE_AVERAGE) != 0;
    const bool normal = (c_method & ColoringMethod::NORMAL_MAP) != 0;
    const bool trap = (c_method & ColoringMethod::ORBIT_TRAP) != 0;
    const bool de_ovl = (c_method & ColoringMethod::DE_OVERLAY) != 0;
    // The total derivative dz/dc is iterated for EDE (distance), the normal map
    // (its argument) and the DE overlay (distance); track it for any of them.
    const bool deriv = (c_method & ColoringMethod::EXTERIOR_DIST_EST) || normal || de_ovl;
    SacAccum sacc; if (sac) sacc.init(sacWindow());
    TrapAccum trapc;
    int i = 1;
    while (i < mxit) {
        tmp = 2.0 * (d_re * z_re - d_im * z_im);
        d_im = 2.0 * (d_re * z_im + d_im * z_re);
        d_re = tmp;

        if (deriv) {
            tmp = 2.0 * (dc_re * z_re - dc_im * z_im) + 1;
            dc_im = 2.0 * (dc_re * z_im + dc_im * z_re);
            dc_re = tmp;
        }
        
        tmp = z_re * z_re - z_im * z_im + c_re;
        z_im = 2.0 * z_re * z_im + c_im;
        z_re = tmp;

        if (sac) sacc.push((double)z_re, (double)z_im);
        if (trap) trapc.push((double)z_re, (double)z_im);

        tmp = z_re * z_re + z_im * z_im;
        if (tmp > _ESCAPE_RADIUS * _ESCAPE_RADIUS) {
            if (c_method & ColoringMethod::EXTERIOR_DIST_EST) {
                return sqrt(tmp) * log(tmp) / sqrt(dc_re * dc_re + dc_im * dc_im);
            } else if (normal) {
                // normal angle = arg(z) - arg(dz/dc); base colour = smooth value.
                if (normalOut) *normalOut = atan2((double)z_im, (double)z_re) - atan2((double)dc_im, (double)dc_re);
                return (i + 1 - log(log(tmp) / 2 / log(2)) / log(2));
            } else if (de_ovl) {
                // base colour = smooth value; overlay = (raw) distance estimate.
                if (normalOut) *normalOut = sqrt(tmp) * log(tmp) / sqrt(dc_re * dc_re + dc_im * dc_im);
                return (i + 1 - log(log(tmp) / 2 / log(2)) / log(2));
            } else if (trap) {
                return trapc.value(i + 1 - log(log(tmp) / 2 / log(2)) / log(2));
            } else if (sac) {
                return sacc.value((double)tmp, (double)_ESCAPE_RADIUS);
            } else {
                return (i + 1 - log(log(tmp) / 2 / log(2)) / log(2));
            }
        }

        tmp = d_re * d_re + d_im * d_im;
        if (tmp < 0.000000001) return trap ? trapc.value(0.0) : -2;   // interior: trap-colour or sentinel
        ++i;
    }
    return trap ? trapc.value(0.0) : -2;   // interior (hit mxit)

}

// AVX2 shallow Mandelbrot with lane refilling: processes a LIST of `count` pixels
// (arbitrary c = cre[k] + i*cim[k]) 4-wide, and the moment a lane escapes/goes
// interior/hits mxit it is immediately reloaded with the next pending pixel, so no
// SIMD lane idles on the divergent boundary (the reason a naive 4-at-a-time SIMD is
// slower than scalar). Writes floatPointCompute-equivalent values (raw EDE distance
// / smooth iter / -2) into out[0..count). Used for both the base row and the
// adaptive supersample sub-pixels.
void Mandel::solveShallowSimdList(const double* cre, const double* cim, int count,
                                  float* out, int mxit, int c_method) const {
    const bool ede = (c_method & ColoringMethod::EXTERIOR_DIST_EST) != 0;
    const double ESC2 = (double)_ESCAPE_RADIUS * _ESCAPE_RADIUS;
    const double LG2 = log(2.0);
    const __m256d two = _mm256_set1_pd(2.0), one = _mm256_set1_pd(1.0);
    const __m256d ESC2v = _mm256_set1_pd(ESC2), intEps = _mm256_set1_pd(1e-9);
    const __m256d mxitv = _mm256_set1_pd((double)mxit);

    alignas(32) double cr_[4], ci_[4], zr_[4], zi_[4], dr_[4], di_[4], dcr_[4], dci_[4], jv_[4];
    int lanePix[4]; int next = 0, activeCount = 0;
    auto loadLane = [&](int l) {
        if (next < count) {
            cr_[l] = cre[next]; ci_[l] = cim[next];
            zr_[l] = cr_[l]; zi_[l] = ci_[l];
            dr_[l] = 2; di_[l] = 0; dcr_[l] = 1; dci_[l] = 0; jv_[l] = 1;
            lanePix[l] = next++; ++activeCount;
        } else {
            cr_[l] = ci_[l] = zr_[l] = zi_[l] = di_[l] = dci_[l] = 0;
            dr_[l] = 2; dcr_[l] = 1; jv_[l] = 1; lanePix[l] = -1;
        }
    };
    for (int l = 0; l < 4; ++l) loadLane(l);

    __m256d cr = _mm256_load_pd(cr_), ci = _mm256_load_pd(ci_);
    __m256d zr = _mm256_load_pd(zr_), zi = _mm256_load_pd(zi_);
    __m256d dr = _mm256_load_pd(dr_), di = _mm256_load_pd(di_);
    __m256d dcr = _mm256_load_pd(dcr_), dci = _mm256_load_pd(dci_);
    __m256d jv = _mm256_load_pd(jv_);

    while (activeCount > 0) {
        // one iteration i = jv:  d'=2zd, dc'=2z*dc+1, z'=z^2+c   (all from old z)
        __m256d ndr = _mm256_mul_pd(two, _mm256_sub_pd(_mm256_mul_pd(dr, zr), _mm256_mul_pd(di, zi)));
        __m256d ndi = _mm256_mul_pd(two, _mm256_add_pd(_mm256_mul_pd(dr, zi), _mm256_mul_pd(di, zr)));
        __m256d ndcr = _mm256_add_pd(_mm256_mul_pd(two, _mm256_sub_pd(_mm256_mul_pd(dcr, zr), _mm256_mul_pd(dci, zi))), one);
        __m256d ndci = _mm256_mul_pd(two, _mm256_add_pd(_mm256_mul_pd(dcr, zi), _mm256_mul_pd(dci, zr)));
        __m256d nzr = _mm256_add_pd(_mm256_sub_pd(_mm256_mul_pd(zr, zr), _mm256_mul_pd(zi, zi)), cr);
        __m256d nzi = _mm256_add_pd(_mm256_mul_pd(two, _mm256_mul_pd(zr, zi)), ci);
        zr = nzr; zi = nzi; dr = ndr; di = ndi; dcr = ndcr; dci = ndci;

        __m256d zrad = _mm256_add_pd(_mm256_mul_pd(zr, zr), _mm256_mul_pd(zi, zi));
        __m256d dmag = _mm256_add_pd(_mm256_mul_pd(dr, dr), _mm256_mul_pd(di, di));
        __m256d jnext = _mm256_add_pd(jv, one);
        int em = _mm256_movemask_pd(_mm256_cmp_pd(zrad, ESC2v, _CMP_GT_OQ));
        int im = _mm256_movemask_pd(_mm256_cmp_pd(dmag, intEps, _CMP_LT_OQ));
        int jm = _mm256_movemask_pd(_mm256_cmp_pd(jnext, mxitv, _CMP_GE_OQ));
        int actmask = (lanePix[0] >= 0) | ((lanePix[1] >= 0) << 1) | ((lanePix[2] >= 0) << 2) | ((lanePix[3] >= 0) << 3);
        int need = (em | im | jm) & actmask;
        if (need == 0) { jv = jnext; continue; }        // fast path: advance i, no scalar work

        // slow path: finish + refill the lanes that ended this iteration.
        alignas(32) double zrad_[4], dcr2_[4], dci2_[4];
        _mm256_store_pd(zrad_, zrad); _mm256_store_pd(dcr2_, dcr); _mm256_store_pd(dci2_, dci);
        _mm256_store_pd(cr_, cr); _mm256_store_pd(ci_, ci);
        _mm256_store_pd(zr_, zr); _mm256_store_pd(zi_, zi);
        _mm256_store_pd(dr_, dr); _mm256_store_pd(di_, di);
        _mm256_store_pd(dcr_, dcr); _mm256_store_pd(dci_, dci);
        _mm256_store_pd(jv_, jv);
        for (int l = 0; l < 4; ++l) {
            if (lanePix[l] < 0) continue;
            bool fin = false; float res = -2.0f;
            if (em & (1 << l)) {
                if (ede) res = (float)(sqrt(zrad_[l]) * log(zrad_[l]) / sqrt(dcr2_[l]*dcr2_[l] + dci2_[l]*dci2_[l]));
                else     res = (float)(jv_[l] + 1 - log(log(zrad_[l]) / 2 / LG2) / LG2);
                fin = true;
            } else if (im & (1 << l)) { res = -2.0f; fin = true; }
            else if (jm & (1 << l)) { res = -2.0f; fin = true; }
            if (fin) { out[lanePix[l]] = res; --activeCount; loadLane(l); }
            else jv_[l] += 1;
        }
        cr = _mm256_load_pd(cr_); ci = _mm256_load_pd(ci_);
        zr = _mm256_load_pd(zr_); zi = _mm256_load_pd(zi_);
        dr = _mm256_load_pd(dr_); di = _mm256_load_pd(di_);
        dcr = _mm256_load_pd(dcr_); dci = _mm256_load_pd(dci_);
        jv = _mm256_load_pd(jv_);
    }
}

float Mandel::accuratePointCompute(mpf_t c_re, mpf_t c_im, int mxit, int c_method) const {
    mpf_t t1, t2;
    mpf_init(t1);
    mpf_init(t2);

    // z = c;
    mpf_t z_re, z_im;
    mpf_init_set(z_re, c_re);
    mpf_init_set(z_im, c_im);
    
    // d = 2.0;
    mpf_t d_re, d_im;
    mpf_init(d_re);
    mpf_init(d_im);
    mpf_set_str(d_re, "2", 10);
    // dc = 1  (derivative w.r.t. the parameter c; drives the EDE distance
    // estimate). Iterated only when EDE is requested, mirroring floatPointCompute.
    const bool ede = (c_method & ColoringMethod::EXTERIOR_DIST_EST) != 0;
    const bool sac = (c_method & ColoringMethod::STRIPE_AVERAGE) != 0;
    SacAccum sacc; if (sac) sacc.init(sacWindow());
    mpf_t dc_re, dc_im, e1, e2;
    mpf_init_set_ui(dc_re, 1); mpf_init(dc_im); mpf_init(e1); mpf_init(e2);
    int i = 1;
    float res = -2;
    while (i < mxit) {
        // d = 2.0 * d * z;
        mpf_mul(t1, d_re, z_im);
        mpf_mul(d_re, d_re, z_re);
        mpf_mul(t2, d_im, z_im);
        mpf_mul(d_im, d_im, z_re);
        mpf_sub(d_re, d_re, t2);
        mpf_mul_ui(d_re, d_re, 2);
        mpf_add(d_im, d_im, t1);
        mpf_mul_ui(d_im, d_im, 2);

        // dc = 2.0 * dc * z + 1   (uses the pre-update z, like floatPointCompute)
        if (ede) {
            mpf_mul(e1, dc_re, z_im); mpf_mul(e2, dc_im, z_re); mpf_add(e2, e1, e2);
            mpf_mul(e1, dc_re, z_re); mpf_mul(t2, dc_im, z_im); mpf_sub(e1, e1, t2);
            mpf_mul_ui(dc_re, e1, 2); mpf_add_ui(dc_re, dc_re, 1);
            mpf_mul_ui(dc_im, e2, 2);
        }

        // z = z * z + c;
        mpf_mul(t1, z_re, z_im);
        mpf_mul(z_re, z_re, z_re);
        mpf_mul(z_im, z_im, z_im);
        mpf_sub(z_re, z_re, z_im);
        mpf_add(z_re, z_re, c_re);
        mpf_mul_ui(z_im, t1, 2);
        mpf_add(z_im, z_im, c_im);

        if (sac) sacc.push(mpf_get_d(z_re), mpf_get_d(z_im));   // ground-truth stripe average

        // auto re_im = z.get_real_imag();
        // if (re_im.first * re_im.first + re_im.second * re_im.second > _ESCAPE_RADIUS) return (i + 1 - log(log(static_cast<float>(re_im.first * re_im.first + re_im.second * re_im.second))) / log(2));
        mpf_mul(t1, z_re, z_re);
        mpf_mul(t2, z_im, z_im);
        mpf_add(t1, t1, t2);
            
        
        double rad = mpf_get_d(t1);

        if (rad > _ESCAPE_RADIUS * _ESCAPE_RADIUS) {
            if (ede) {
                // EDE distance estimate = |z| * log(|z|^2) / |dc|, computed via
                // mpf to avoid overflowing double when |dc| is huge (deep zoom).
                mpf_mul(e1, dc_re, dc_re); mpf_mul(e2, dc_im, dc_im);
                mpf_add(e1, e1, e2); mpf_sqrt(e1, e1);
                mpf_set_d(e2, sqrt(rad) * log(rad));
                mpf_div(e1, e2, e1);
                res = (float)mpf_get_d(e1);
            } else if (sac) {
                res = sacc.value(rad, (double)_ESCAPE_RADIUS);
            } else {
                res = (i + 1 - log(log(rad) / 2 / log(2)) / log(2));
            }
            break;
        }

        // re_im = d.get_real_imag();
        // if (re_im.first * re_im.first + re_im.second * re_im.second < 0.000000001) return -2;
        mpf_mul(t1, d_re, d_re);
        mpf_mul(t2, d_im, d_im);
        mpf_add(t1, t1, t2);
        rad = mpf_get_d(t1);
        if (rad < 0.000000001) break;
        
        ++i;
    }
    mpf_clear(d_re);
    mpf_clear(d_im);
    mpf_clear(dc_re);
    mpf_clear(dc_im);
    mpf_clear(e1);
    mpf_clear(e2);
    mpf_clear(z_re);
    mpf_clear(z_im);
    mpf_clear(t1);
    mpf_clear(t2);
    return res;
}

// Sensitivity helper for the mxit-boundary check: instead of re-running a full
// accuratePointCompute of the whole reference orbit (mxit+64 GMP iterations) just
// to learn whether the reference escapes just past mxit, CONTINUE the orbit we
// already built for `extra` more iterations from its tail Z_last. Same orbit, same
// escape decision, but ~mxit/extra times cheaper. Precision-consistent (uses the
// reference's own representation: BigFixed tail when active, else the mpf tail).
bool Mandel::refTailEscapes(int last, int extra) const {
    mpf_t zr, zi, t1, t2;
    mpf_init(zr); mpf_init(zi); mpf_init(t1); mpf_init(t2);
    if (_use_bigfixed) { bf_to_mpf(zr, _bz_re[last & 1]); bf_to_mpf(zi, _bz_im[last & 1]); }
    else { mpf_set(zr, _z_re[last & 1]); mpf_set(zi, _z_im[last & 1]); }
    const double R2 = (double)_ESCAPE_RADIUS * (double)_ESCAPE_RADIUS;
    bool esc = false;
    for (int k = 0; k < extra && !esc; ++k) {
        mpf_mul(t1, zr, zi);            // zr*zi
        mpf_mul(zr, zr, zr);            // zr^2
        mpf_mul(t2, zi, zi);            // zi^2
        mpf_sub(zr, zr, t2);            // zr^2 - zi^2
        mpf_add(zr, zr, _ref_z_re);     // + cr
        mpf_mul_ui(zi, t1, 2);          // 2 zr zi
        mpf_add(zi, zi, _ref_z_im);     // + ci
        mpf_mul(t1, zr, zr); mpf_mul(t2, zi, zi); mpf_add(t1, t1, t2);
        if (mpf_get_d(t1) > R2) esc = true;
    }
    mpf_clears(zr, zi, t1, t2, (mpf_ptr)0);
    return esc;
}


void Mandel::ComputeDirect(int mxit, float* out, int step, int c_method) {
    // Brute-force ground truth: iterate every sampled pixel in full mpf_t
    // precision. Uses _c0_re/_c0_im/_dx/_dy from the most recent Compute() call.
    if (step < 1) step = 1;
    // EDE values are normalised by the pixel spacing, matching the engine paths
    // (e.g. `_iter[idx] /= dx_f`), so the oracle and engine EDE are comparable.
    const bool ede = (c_method & ColoringMethod::EXTERIOR_DIST_EST) != 0;
    const double dxf = mpf_get_d(_dx);
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < _h; i += step) {
        if (_flag_halt) continue;
        mpf_t cre, cim, tmp;
        mpf_init(cre);
        mpf_init(cim);
        mpf_init(tmp);
        for (int j = 0; j < _w; j += step) {
            mpf_mul_ui(tmp, _dx, j);
            mpf_add(cre, _c0_re, tmp);
            mpf_mul_ui(tmp, _dy, i);
            mpf_add(cim, _c0_im, tmp);
            float v = accuratePointCompute(cre, cim, mxit, c_method);
            if (ede && v >= 0) v /= dxf;
            out[i * _w + j] = v;
        }
        mpf_clear(cre);
        mpf_clear(cim);
        mpf_clear(tmp);
    }
}

int Mandel::findNucleus(mpf_t nuc_re, mpf_t nuc_im, int maxp) {
    const mp_bitcnt_t prec = mpf_get_prec(_c0_re);
    mpf_t cre, cim, r, dh, dist, dr, di;
    mpf_init2(cre, prec); mpf_init2(cim, prec); mpf_init2(r, prec); mpf_init2(dh, prec);
    mpf_init2(dist, prec); mpf_init2(dr, prec); mpf_init2(di, prec);
    // View centre = c0 + (halfwidth, halfheight); search radius = halfwidth = 2/scale.
    mpf_set_ui(r, 2); mpf_div(r, r, _scale);
    mpf_add(cre, _c0_re, r);
    mpf_mul_ui(dh, r, _h); mpf_div_ui(dh, dh, _w);
    mpf_add(cim, _c0_im, dh);
    int period = g_ballPeriod(cre, cim, r, maxp, prec);
    int result = 0;
    if (period > 0) {
        mpf_set(nuc_re, cre); mpf_set(nuc_im, cim);
        if (g_newton(nuc_re, nuc_im, period, prec, 200)) {
            // Accept only if the refined nucleus is inside the search disk.
            mpf_sub(dr, nuc_re, cre); mpf_sub(di, nuc_im, cim);
            mpf_mul(dr, dr, dr); mpf_mul(di, di, di); mpf_add(dist, dr, di); mpf_sqrt(dist, dist);
            if (mpf_cmp(dist, r) <= 0) result = period;
        }
    }
    mpf_clears(cre, cim, r, dh, dist, dr, di, (mpf_ptr)0);
    return result;
}

int Mandel::createPeriodicRef(int period, int mxit, int c_method) {
    // The reference IS the nucleus (bounded, exactly period-p). Build ONE period in
    // GMP (reusing calCoefficient, which fills every double/floatexp shadow); the
    // delta loop then indexes it modulo p (Z_i = Z_{i-p}), so no mxit-long orbit is
    // materialised.  _ref_z is already the nucleus.
    (void)mxit;
    _ref_z_f = Comp{ mpf_get_ld(_ref_z_re), mpf_get_ld(_ref_z_im) };
    _ref_virtual = true;   // the nucleus is not a rendered pixel (skip centre-pixel erase; every pixel is computed)
    // Fractional pixel coordinate of the nucleus. The floatexp dc = _dxfe*(px-_ref_x)
    // uses it directly; the double path uses the nearest pixel _ref plus the
    // sub-pixel remainder _ref_frac so its integer-pixel dc equals c_pixel - nucleus.
    mpf_sub(_t1, _ref_z_re, _c0_re); mpf_div(_t1, _t1, _dx); _ref_x = mpf_get_d(_t1);
    mpf_sub(_t1, _ref_z_im, _c0_im); mpf_div(_t1, _t1, _dy); _ref_y = mpf_get_d(_t1);
    int px = (int)std::lround(_ref_x), py = (int)std::lround(_ref_y);
    _ref = { py, px, 0, 0 };
    _ref_frac_re = (Float)((_ref_x - px) * mpf_get_ld(_dx));   // nucleus - pixel_c(_ref)
    _ref_frac_im = (Float)((_ref_y - py) * mpf_get_ld(_dy));

    // Z_0 = nucleus (mirror createRef's index-0 setup).
    mpf_set(_z_re[0], _ref_z_re); mpf_set(_z_im[0], _ref_z_im);
    _use_bigfixed = false;   // periodic nucleus orbit is short (one period) and does not init BigFixed state; use mpf
    _zf[0] = _ref_z_f; _zfr[0] = _ref_z_f.real(); _zfi[0] = _ref_z_f.imag();
    _zfr_fe[0] = mpf_to_fe(_ref_z_re); _zfi_fe[0] = mpf_to_fe(_ref_z_im);
    _z_m3[0] = (_zfr[0] * _zfr[0] + _zfi[0] * _zfi[0]) / 1000000;
    _df[0] = Comp{ 1 }; _dfr[0] = 1.0; _dfi[0] = 0.0;
    _dfe_r = FloatExp{ 1.0, 0 }; _dfe_i = FloatExp{ 0.0, 0 };
    // No series approximation: the double path's SA-init loop (order 0, Adf[0]=SA_delta)
    // then produces dz = dc, and each pixel starts at reference index 1 (_SA_it = 0).
    mpf_mul_ui(_t1, _dx, _w); mpf_div_ui(_t1, _t1, 2);
    mpf_mul_ui(_t2, _dy, _h); mpf_div_ui(_t2, _t2, 2);
    _SA_delta = { mpf_get_ld(_t1), mpf_get_ld(_t2) };
    _SA_flag = false; _SA_it = 0; _SA_order = 0;
    for (int i = 0; i < _SA_N; ++i) _Adf_old[i] = _Bdf_old[i] = 0;
    _Adf_old[0] = _SA_delta;

    // One period: Z_1..Z_{period-1} via the exact GMP recurrence (pr_it = 0 so no
    // series-approximation work). The nucleus orbit is bounded (never escapes).
    // Force the floatexp shadow fill so Compute can escalate double->floatexp after
    // inspecting the reference (the double path itself ignores _zfr_fe).
    bool save_fe = _use_floatexp; _use_floatexp = true;
    for (int i = 1; i < period; ++i) calCoefficient(i, 0, c_method);
    _use_floatexp = save_fe;

    _ref_period = period;
    return period;
}


// Build the deep-zoom reference orbit and return its length. Extracted verbatim
// from Compute's method==1 path: try a periodic-nucleus reference when eligible,
// else a full reference from the centre pixel; then apply the content-aware
// floatexp escalation (rebuild in floatexp if the orbit passes too close to zero)
// and the mxit-boundary sensitivity rebuild (from the exact centre). Sets
// _ref_period / _ref_bounded / _use_floatexp and accumulates the pf_ref timer.
int Mandel::buildReferenceOrbit(std::set<std::array<int, 4>>& s, mpf_t scale, int mxit,
                                int c_method, bool profile, double& pf_ref) {
    double tk;
    int ref_it = 0;
    _ref_period = 0;
    // Periodic-reference optimisation: if a minibrot nucleus sits in the view, use
    // it as a bounded period-p reference (build ONE period, index modulo p) instead
    // of a full mxit-long orbit. Auto-enabled for the floatexp deep path
    // (scale > 1e280); MANDEL_PERIODIC=1 forces on, =0 off. Requires BLA on and
    // non-EDE/non-SAC (the mod-p index / skip cap live only in those paths).
    bool periodic_done = false;
    static int periodic_env = -1;
    if (periodic_env < 0) { const char* e = getenv("MANDEL_PERIODIC"); periodic_env = e ? atoi(e) : -1; }
    bool periodic_want = (periodic_env == 1) ||
                         (periodic_env != 0 && mpf_cmp_d(scale, 1e280) > 0);
    if (periodic_want && _use_bla &&
        !(c_method & ColoringMethod::EXTERIOR_DIST_EST) &&
        !(c_method & ColoringMethod::NORMAL_MAP) &&
        !(c_method & ColoringMethod::DE_OVERLAY) &&
        !(c_method & ColoringMethod::ORBIT_TRAP) &&
        !(c_method & ColoringMethod::STRIPE_AVERAGE)) {
        tk = nowSec();
        int p = findNucleus(_ref_z_re, _ref_z_im, mxit);
        if (p > 2 && p <= mxit) {
            ref_it = createPeriodicRef(p, mxit, c_method);
            _ref_bounded = true;
            periodic_done = true;
            // Content-aware escalation: if the nucleus orbit passes close enough to
            // zero that 2*Z*dz underflows, run the exact floatexp path.
            if (!_use_floatexp) {
                double refmin2 = 1e300; bool uf = false;
                for (int q = 1; q < p; ++q) { double m2 = _z_m3[q] * 1e6; if (m2 <= 0.0) { uf = true; break; } if (m2 < refmin2) refmin2 = m2; }
                double dcmax = 2.0 / mpf_get_d(scale);
                if (uf || std::sqrt(refmin2) * dcmax < 1e-300) _use_floatexp = true;
            }
            if (profile) fprintf(stderr, "  [profile] periodic reference: period=%d floatexp=%d\n", p, _use_floatexp ? 1 : 0);
        } else if (profile) {
            fprintf(stderr, "  [profile] periodic reference: no in-view nucleus (period=%d), using full reference\n", p);
        }
        pf_ref += nowSec() - tk;
    }
    if (!periodic_done) {
        s.insert({ _h / 2, _w / 2, 0, 0 });
        tk = nowSec(); ref_it = createRef(s, mxit, mxit, true, c_method); pf_ref += nowSec() - tk;
        // Content-aware floatexp escalation. Even below the hard 1e280 gate the plain
        // double delta loop loses precision when the reference orbit passes close to
        // zero (2*Z*dz falls into the denormal range for far pixels). Detect it from
        // the just-built reference's closest approach to zero and redo it in floatexp;
        // references that stay clear of zero (e.g. flake@1e157) keep the fast double
        // path. Only meaningful once dz^2 can underflow (scale > ~1e150).
        if (!_use_floatexp && mpf_cmp_d(scale, 1e150) > 0) {
            double refmin2 = 1e300; bool underflowed = false;
            for (int q = 1; q <= ref_it; ++q) {
                double m2 = _z_m3[q] * 1e6;                 // |Z_q|^2 (double shadow)
                if (m2 <= 0.0) { underflowed = true; break; }   // |Z| below ~1.5e-154
                if (m2 < refmin2) refmin2 = m2;
            }
            double dcmax = 2.0 / mpf_get_d(scale);          // ~max pixel |dc|
            if (underflowed || std::sqrt(refmin2) * dcmax < 1e-300) {
                _use_floatexp = true;
                s.insert({ _h / 2, _w / 2, 0, 0 });
                tk = nowSec(); ref_it = createRef(s, mxit, mxit, true, c_method, true); pf_ref += nowSec() - tk;
            }
        }
        if (_use_floatexp && ref_it >= mxit - 16) {
            // "Sensitive" = the reference escapes at or just beyond mxit (its
            // classification then depends on accumulated double rounding), so it is
            // rebuilt from the exact geometric centre. Continue the orbit we already
            // built ~80 iterations past its tail instead of re-running a full mpf
            // accuratePointCompute of the whole orbit + derivative (which cost ~3s at
            // 1e876, dwarfing the reference build). refTailEscapes uses the
            // reference's own representation so it is exact and at least as sensitive.
            bool sensitive = ref_it < mxit || refTailEscapes(ref_it, 80);
            if (sensitive) {
                // Rebuild from the exact geometric centre so the result cannot depend
                // on which centre pixel exists for even/odd dimensions.
                s.insert({ _h / 2, _w / 2, 0, 0 });
                tk = nowSec(); ref_it = createRef(s, mxit, mxit, true, c_method, true); pf_ref += nowSec() - tk;
            }
        }
        _ref_bounded = (ref_it >= mxit);
    }
    return ref_it;
}


void Mandel::Compute(mpf_t c_re, mpf_t c_im, mpf_t scale, int mxit, int c_method,
                     int full_h, int row_base) {
    std::cout << "mxit: " << mxit << '\n';
    _progress_total = 0;
    progressSet(0.0);
    _mx_coef = -1;
    _ref_cnt = 0;
    int iteration = 0;

    // _scale = scale;
    mpf_set(_scale, scale);

    // Strip export: this Mandel is only `_h` rows tall but represents rows
    // [row_base, row_base+_h) of a taller `full_h`-row image, so derive the grid
    // from full_h and shift the origin down by row_base rows. full_h==0 => the
    // Mandel is the whole image (interactive path, unchanged).
    int rows = (full_h > 0) ? full_h : _h;
    mpf_t dw, dh;
    mpf_init_set_ui(dw, 2);
    mpf_div(dw, dw, scale);
    mpf_init_set(dh, dw);
    mpf_div_ui(dh, dh, _w);
    mpf_mul_ui(dh, dh, rows);
    mpf_sub(_c0_re, c_re, dw);
    mpf_sub(_c0_im, c_im, dh);

    mpf_div_ui(_dx, dw, _w - 1);
    mpf_mul_ui(_dx, _dx, 2);
    mpf_div_ui(_dy, dh, rows - 1);
    mpf_mul_ui(_dy, _dy, 2);

    if (row_base > 0) {                 // move origin to this strip's first row
        mpf_t off; mpf_init(off);
        mpf_mul_ui(off, _dy, row_base);
        mpf_add(_c0_im, _c0_im, off);
        mpf_clear(off);
    }

    mpf_clear(dw);
    mpf_clear(dh);
    
    std::set<std::array<int, 4>> s;

    int method = 0;

    // Shallow zoom: method = 0, basic algorithm
    // Deep zoom: method = 1, Perturbation + Series Approximation + Rebase
    if (mpf_cmp_d(scale, 1e6) > 0) method = 1;
    if (getenv("MANDEL_FORCE_PERT")) method = 1;   // exercise perturbation at any scale
    
    if (_flag_halt) return;
    int ref_it = 0, pr_it = 0;
    if (method == 0) {
        Float c0_re_f = mpf_get_ld(_c0_re);
        Float c0_im_f = mpf_get_ld(_c0_im);
        Float dx_f = mpf_get_ld(_dx);
        Float dy_f = mpf_get_ld(_dy);
        static int simd0 = -1;
        if (simd0 < 0) { const char* e = getenv("MANDEL_SIMD"); simd0 = e ? atoi(e) : 1; }
        const bool ede = (c_method & ColoringMethod::EXTERIOR_DIST_EST) != 0;
        const bool normalmap = (c_method & ColoringMethod::NORMAL_MAP) != 0;
        const bool trapmap = (c_method & ColoringMethod::ORBIT_TRAP) != 0;
        const bool deovl = (c_method & ColoringMethod::DE_OVERLAY) != 0;
        // SAC, the normal map, orbit traps and the DE overlay take the scalar path.
        const bool simd0on = simd0 && !(c_method & ColoringMethod::STRIPE_AVERAGE) && !normalmap && !trapmap && !deovl;   // shallow SIMD has no SAC
        // Coarse-to-fine: a ~1/16-work strided pass paints the whole frame with a
        // blocky preview (picked up immediately by the async display) before the
        // full-resolution pass below sharpens it. ~instant feedback on every view.
        static int coarse_on = -1;
        if (coarse_on < 0) { const char* e = getenv("MANDEL_COARSE"); coarse_on = e ? atoi(e) : 1; }
        if (coarse_on && !_flag_halt) {
            const int C = 4;
#pragma omp parallel for schedule(dynamic, 2)
            for (int ci = 0; ci < _h; ci += C) {
                if (_flag_halt) continue;
                Float cim = c0_im_f + dy_f * ci;
                for (int cj = 0; cj < _w; cj += C) {
                    double val = floatPointCompute(c0_re_f + dx_f * cj, cim, mxit, c_method);
                    if (ede && val >= 0) val /= dx_f;
                    for (int bi = ci; bi < ci + C && bi < _h; ++bi)
                        for (int bj = cj; bj < cj + C && bj < _w; ++bj)
                            _iter[getIndex(bi, bj, 0, 0)] = (float)val;
                }
            }
        }
        progressBegin(_h, coarse_on ? 0.05 : 0.0, coarse_on ? 0.95 : 1.0);
#pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < _h; ++i) {
            if (_flag_halt) continue;
            Float cim = c0_im_f + dy_f * i;
            if (simd0on) {
                std::vector<double> cre(_w), cimv(_w);
                std::vector<float> row(_w);
                for (int j = 0; j < _w; ++j) { cre[j] = c0_re_f + dx_f * j; cimv[j] = cim; }
                solveShallowSimdList(cre.data(), cimv.data(), _w, row.data(), mxit, c_method);
                for (int j = 0; j < _w; ++j) {
                    double val = row[j];
                    if (ede && val >= 0) val /= dx_f;
                    _iter[getIndex(i, j, 0, 0)] = (float)val;
                }
            } else {
                for (int j = 0; j < _w; ++j) {
                    if (_flag_halt) break;
                    int idx = getIndex(i, j, 0, 0);
                    if (normalmap && _normal) {
                        double nrm = 0.0;
                        _iter[idx] = (float)floatPointCompute(c0_re_f + dx_f * j, cim, mxit, c_method, &nrm);
                        _normal[idx] = (float)nrm;
                    } else if (deovl && _normal) {
                        double de = 0.0;
                        _iter[idx] = (float)floatPointCompute(c0_re_f + dx_f * j, cim, mxit, c_method, &de);
                        _normal[idx] = (_iter[idx] < 0) ? 0.0f : (float)(de / dx_f);   // pixel-normalised DE
                    } else {
                        _iter[idx] = floatPointCompute(c0_re_f + dx_f * j, cim, mxit, c_method);
                        if (ede && _iter[idx] >= 0) _iter[idx] /= dx_f;
                    }
                }
            }
            progressAdvance();
        }
    }
    else if (method == 1) {
        // BLA (bivariate linear approximation) skips runs of reference iterations for
        // a big deep-zoom speedup (~6x on 1e50+ high-iteration views). It is now
        // CLASSIFICATION-safe (periodicity detector restarted after each skip) and
        // carries the EDE derivative (BLA-with-derivative), so it is accuracy-safe
        // under distance estimation too. Ground-truth diff at 3.8e51/mxit=2M shows
        // BLA-on matches the GMP oracle BETTER than BLA-off (0 vs 521/12288 pixel
        // misclassifications; equal escape-time floor) because rebasing avoids the
        // naive-perturbation double-precision drift. Hence ON by default; set
        // MANDEL_BLA=0 to force it off.
        { const char* e = getenv("MANDEL_BLA"); _use_bla = e ? (atoi(e) != 0) : true; }
        { const char* e = getenv("MANDEL_BLA_EPS"); _bla_eps = e ? atof(e) : 0.0; }
        { const char* e = getenv("MANDEL_BLA_MINSKIP"); int ms = e ? atoi(e) : 8;
          _bla_minlevel = 0; while ((1 << (_bla_minlevel + 1)) <= ms) ++_bla_minlevel; }
        if (g_bla_noescape < 0) { const char* e = getenv("MANDEL_BLA_NOESCAPE"); g_bla_noescape = e ? atoi(e) : 0; }
        memset(g_bla_stat, 0, sizeof(g_bla_stat));
        memset(g_fe_stat, 0, sizeof(g_fe_stat));
        _simd_measured = false; _simd_bla_idle = false;   // re-measured per frame from the coarse probe
        { const char* e = getenv("MANDEL_INTERIOR"); _use_interior = !e || atoi(e); }
        { const char* e = getenv("MANDEL_INT_EPS"); double ep = e ? atof(e) : 1e-13; _interior_eps2 = ep * ep; }
        { const char* e = getenv("MANDEL_INT_CONFIRM"); _interior_confirm = e ? atoi(e) : 4; }
        // Deep-zoom: below double's ~1e300 representability wall the delta must be
        // rescaled (floatexp scale + double delta). This hard gate catches the
        // representability limit; a finer precision gate is applied after the
        // reference is built (see the escalation below).
        { const char* e = getenv("MANDEL_FE");
          _use_floatexp = e ? atoi(e) != 0 : (mpf_cmp_d(scale, 1e280) > 0); }
        g_fe_fallback = 0;
        memset(g_blafe_safe, 0, sizeof(g_blafe_safe));
        _fe_cutoff_sensitive = false;
        if (g_int_rep == -2) { const char* e = getenv("MANDEL_INT_REP"); g_int_rep = e ? atoi(e) : 1; }   // default ON
        _dxfe = mpf_to_fe(_dx); _dyfe = mpf_to_fe(_dy);
        double pf_ref = 0, pf_step = 0, pf_bla = 0;
        const bool profile = getenv("MANDEL_PROFILE") != nullptr;
        auto now = [] { return std::chrono::duration_cast<std::chrono::duration<double>>(
                            std::chrono::high_resolution_clock::now().time_since_epoch()).count(); };
        double tk;
        Float c0_re_f = mpf_get_ld(_c0_re);
        Float c0_im_f = mpf_get_ld(_c0_im);
        floatPointCompute(c0_re_f, c0_im_f, mxit, c_method);
        ref_it = buildReferenceOrbit(s, scale, mxit, c_method, profile, pf_ref);
        if (_use_bla) { tk = now(); buildBLA(ref_it, (c_method & (ColoringMethod::EXTERIOR_DIST_EST | ColoringMethod::NORMAL_MAP | ColoringMethod::DE_OVERLAY)) != 0); pf_bla += now() - tk; }
        // Reference-orbit stripe prefix sum, so a BLA skip can restore its omitted
        // stripe-average contributions (z ~ X during a valid skip) instead of
        // dropping them -- otherwise the dropped fraction grows with distance from
        // the reference, giving Feather a pan-dependent radial halo.
        if (_use_bla && _use_floatexp && (c_method & ColoringMethod::STRIPE_AVERAGE) && sacWindow() <= 0) {
            const double* zfp = reinterpret_cast<const double*>(_zf);
            int N = ref_it + 1;
            _sacRefPre.assign((size_t)N + 1, 0.0);
            for (int m = 1; m <= N; ++m) {
                double xr = zfp[2 * (m - 1)], xi = zfp[2 * (m - 1) + 1];
                _sacRefPre[m] = _sacRefPre[m - 1] + (0.5 + 0.5 * sin(7.0 * atan2(xi, xr)));
            }
        } else {
            _sacRefPre.clear();
        }
        // Coarse strided grid computed once and reused for two purposes:
        //   (a) interior auto-gate: if no probed pixel is interior, disable the
        //       (costly) periodicity detection so exterior-only frames pay ~zero.
        //   (b) blocky preview: paint each computed grid point across its C x C
        //       block so the deep path shows a preview that SHARPENS in place
        //       instead of scattered single points ("measles") -- the same
        //       coarse-to-fine feedback the shallow path already gives.
        static int coarse_on = -1;
        if (coarse_on < 0) { const char* e = getenv("MANDEL_COARSE"); coarse_on = e ? atoi(e) : 1; }
        if ((_use_interior || coarse_on) && !_flag_halt) {
            const int C = 8;
            std::set<std::array<int, 4>> probe;
            for (int i = 0; i < _h; i += C)
                for (int j = 0; j < _w; j += C)
                    probe.insert({ i, j, 0, 0 });
            if (!probe.empty() && !_flag_halt) stepParallel(probe, ref_it, mxit, c_method);
            // Measure how much BLA actually skipped on the probe. If it barely skips
            // (compute-bound), the SIMD double kernel wins; if it skips most steps
            // (memory-bound deep zoom), SIMD is ~1x/slightly worse, so stay scalar.
            // This measured gate is more reliable than the rmax<dcmax proxy: a view
            // can have rmax>dcmax yet still never trigger a valid skip (e.g. night
            // city at 1.7e40) -- there SIMD still wins.
            {
                long long sk = 0, no = 0;
                for (int t = 0; t < 64; ++t) { sk += g_bla_stat[t][0]; no += g_bla_stat[t][2]; }
                double frac = (sk + no) ? (double)sk / (double)(sk + no) : 0.0;
                _simd_bla_idle = !_use_bla || (frac < 0.5);
                _simd_measured = true;
                if (profile)
                    fprintf(stderr, "  [profile] SIMD gate: probe BLA skip-frac=%.3f -> %s\n",
                            frac, _simd_bla_idle ? "SIMD (compute-bound)" : "scalar (BLA-effective)");
            }
            if (_use_interior) {
                int nint = 0;
                for (int i = 0; i < _h; i += C)
                    for (int j = 0; j < _w; j += C)
                        if (_iter[getIndex(i, j, 0, 0)] == -2) ++nint;
                if (nint == 0) _use_interior = false;
            }
            if (coarse_on && !_flag_halt) {
                // Replicate each grid sample across its block (centre subpixels).
                // Glitched grid points stay EMPTYPIXEL -> a few holes, filled by the
                // full pass below; far better than a field of isolated dots. The
                // image-centre pixel is the reference and is NOT recomputed by the
                // full pass (it is erased from the set below), so preserve it.
                int cidx = getIndex(_h / 2, _w / 2, 0, 0);
                float saveCenter = _iter[cidx];
                for (int i = 0; i < _h; i += C) {
                    for (int j = 0; j < _w; j += C) {
                        float v = _iter[getIndex(i, j, 0, 0)];
                        if (v == EMPTYPIXEL) continue;
                        for (int bi = i; bi < i + C && bi < _h; ++bi)
                            for (int bj = j; bj < j + C && bj < _w; ++bj)
                                if (bi != i || bj != j)
                                    _iter[getIndex(bi, bj, 0, 0)] = v;
                    }
                }
                _iter[cidx] = saveCenter;
            }
        }
        // Insert all pixels to render. Pixels already resolved by the interior
        // probe get recomputed (a negligible ~1/64 fraction); this avoids relying
        // on the caller having pre-initialised _iter to EMPTYPIXEL.
        for (int i = 0; i < _h; ++i)
            for (int j = 0; j < _w; ++j)
                s.insert({ i, j, 0, 0 });
        if (!_ref_virtual)
            s.erase({ _h / 2, _w / 2, 0, 0 });   // pixel reference was already computed
        progressBegin((int)s.size(), 0.08, 0.92);
        
        int ref_pass = 0;
        while (!s.empty()) {
            if (_flag_halt) return;
            size_t pending_before = s.size();
            tk = now(); stepParallel(s, ref_it, mxit, c_method);
            double pass_time = now() - tk; pf_step += pass_time;
            if (profile)
                fprintf(stderr, "  [profile] ref-pass=%d ref-it=%d pending=%zu resolved=%zu remaining=%zu delta=%.3f s\n",
                        ref_pass, ref_it, pending_before, pending_before - s.size(), s.size(), pass_time);
            ++ref_pass;
            if (_flag_halt) return;
            if (s.empty()) break;
            tk = now(); ref_it = createRef(s, mxit, mxit, false, c_method); pf_ref += now() - tk;
            _ref_bounded = (ref_it >= mxit);
            // createRef removes the chosen reference pixel from s. If it was the
            // final unresolved base pixel AND no adaptive-SS pass follows, there is
            // no next delta pass, so building a potentially mxit-long BLA table is
            // pure waste. With SS enabled, subpixels may still need this new
            // reference; its previous BLA is stale and must be rebuilt.
            const bool ss_follows =
                (_sub > 1) && (c_method & ColoringMethod::SUPER_SAMPLING);
            if ((!s.empty() || ss_follows) && _use_bla) {
                tk = now(); buildBLA(ref_it, (c_method & (ColoringMethod::EXTERIOR_DIST_EST | ColoringMethod::NORMAL_MAP | ColoringMethod::DE_OVERLAY)) != 0);
                pf_bla += now() - tk;
            }
        }
        if (profile) {
            long long sk = 0, ap = 0, no = 0;
            for (int t = 0; t < 64; ++t) { sk += g_bla_stat[t][0]; ap += g_bla_stat[t][1]; no += g_bla_stat[t][2]; }
            fprintf(stderr, "  [profile] reference-orbit (GMP): %.3f s   delta-loop (double): %.3f s   buildBLA: %.3f s   refs=%d\n",
                    pf_ref, pf_step, pf_bla, _ref_cnt);
            if (_use_bla)
                fprintf(stderr, "  [profile] BLA: applies=%lld skipped=%lld normal-steps=%lld  avg-skip=%.1f  skip-frac=%.1f%%\n",
                        ap, sk, no, ap ? (double)sk / ap : 0.0, (sk + no) ? 100.0 * sk / (sk + no) : 0.0);
            if (_use_floatexp) {
                long long fsk = 0, fap = 0, fst = 0;
                for (int t = 0; t < 64; ++t) { fsk += g_fe_stat[t][0]; fap += g_fe_stat[t][1]; fst += g_fe_stat[t][2]; }
                fprintf(stderr, "  [profile] FE-BLA: applies=%lld skip-iters=%lld normal-steps=%lld  avg-skip=%.1f  skip-frac=%.1f%%\n",
                        fap, fsk, fst, fap ? (double)fsk / fap : 0.0, (fsk + fst) ? 100.0 * fsk / (fsk + fst) : 0.0);
                fprintf(stderr, "  [profile] FE: cutoff-sensitive=%d GMP-fallback-pixels=%lld\n",
                        _fe_cutoff_sensitive ? 1 : 0, g_fe_fallback);
                long long safe = 0; for (int t = 0; t < 64; ++t) safe += g_blafe_safe[t];
                fprintf(stderr, "  [profile] FE-BLA overflow-safe fallbacks=%lld\n", safe);
            }
        }
    }
    bool is_super_sampling = (_sub > 1) && (c_method & ColoringMethod::SUPER_SAMPLING);
    if (!is_super_sampling) { progressSet(1.0); return; }
    if (_flag_halt) return;
    
    // Oversampling pixels that differs from neighbours with (_sub x _sub) pixels
    std::vector<std::array<int, 2>> v;
    int mix_cnt = 0;
    // Flag detector threshold: supersample a pixel when its iteration value differs
    // from a neighbour by more than exp(log(mxit)/K). Larger K => lower threshold =>
    // more pixels supersampled (catches colour-complex pixels the base AA can't
    // resolve). Tunable via MANDEL_SS_K (default 8).
    double ss_k = 8.0;
    { const char* e = getenv("MANDEL_SS_K"); if (e) ss_k = atof(e); }
    double ss_thresh = log((double)mxit) / ss_k;
    // Feather (SAC) stores a [0,1] stripe value, not an iteration count, so the
    // log(mxit) threshold below never fires for it. The exterior stripe pattern is
    // high-frequency *everywhere*, so cheap adaptive SS can't clean it (it would
    // have to supersample ~every pixel -> use Export's uniform SS for that). By
    // default we only supersample the set boundary (interior stores -2, exterior
    // is >=0, so a sign change marks the edge) -- cheap and it removes the most
    // visible jaggies. Set MANDEL_SAC_SS_STEPS>0 to also flag stripe edges whose
    // colour changes by more than that many palette steps/pixel (colour index =
    // stripe*density/20*colP, colP=2048); lower = more supersampling (slower).
    const bool sac = (c_method & ColoringMethod::STRIPE_AVERAGE) != 0;
    double sac_steps = 0.0;   // 0 = boundary only
    { const char* e = getenv("MANDEL_SAC_SS_STEPS"); if (e) sac_steps = atof(e); }
    const double sac_gain = (double)_ss_density * (2048.0 / 20.0);
    // Orbit trap stores a small trap VALUE (not an iteration count), so the log(mxit)
    // iteration-difference threshold below never fires -> SS had no effect. Flag the
    // set boundary plus pixels whose trap colour changes by more than TRAP_SS_STEPS
    // palette entries/pixel (colour index = trap*density/60*colP). Trap detail is
    // high-frequency, so keep the default moderate; lower = more SS (slower).
    const bool trap = (c_method & ColoringMethod::ORBIT_TRAP) != 0;
    double trap_steps = 3.0;
    { const char* e = getenv("MANDEL_TRAP_SS_STEPS"); if (e) trap_steps = atof(e); }
    const double trap_gain = (double)_ss_density / 60.0 * 2048.0;
    for (int i = 0; i < _h; ++i) {
        for (int j = 0; j < _w; ++j) {
            bool need_sample = false;
            if (c_method & ColoringMethod::EXTERIOR_DIST_EST) {
                for (int xi = -1; xi <= 1; ++xi) {
                    for (int yi = -1; yi <= 1; ++yi) {
                        if (xi == 0 && yi == 0) continue;
                        int ny = i + yi;
                        int nx = j + xi;
                        if (nx >= 0 && nx < _w && ny >= 0 && ny < _h) {
                            if (_iter[getIndex(i, j, 0, 0)] * _iter[getIndex(ny, nx, 0, 0)] < 0) {
                                need_sample = true;
                                break;
                            }
                        }
                    }
                    if (need_sample) break;
                }
            } else {
                float cv = _iter[getIndex(i, j, 0, 0)];
                float diff = 0; bool boundary = false;
                for (int xi = -1; xi <= 1; ++xi) {
                    for (int yi = -1; yi <= 1; ++yi) {
                        if (xi == 0 && yi == 0) continue;
                        int ny = i + yi;
                        int nx = j + xi;
                        if (nx >= 0 && nx < _w && ny >= 0 && ny < _h) {
                            float nv = _iter[getIndex(ny, nx, 0, 0)];
                            diff = std::max(diff, std::abs(cv - nv));
                            if ((sac || trap) && cv * nv < 0) boundary = true;   // interior/exterior edge
                        }
                    }
                }
                if (sac)       need_sample = boundary || (sac_steps > 0 && diff * sac_gain > sac_steps);
                else if (trap) need_sample = boundary || (diff * trap_gain > trap_steps);
                else           need_sample = (log(diff) > ss_thresh);
            }
            if (need_sample) {
                ++mix_cnt;
                v.push_back({ i, j });
                for (int xi = -_sub / 2; xi <= _sub / 2; ++xi) {
                    for (int yi = -_sub / 2; yi <= _sub / 2; ++yi) {
                        if (xi == 0 && yi == 0) continue;
                        s.insert({ i, j, yi, xi });
                    }
                }
            }
        }
    }
    if (getenv("MANDEL_PROFILE")) {
        FILE* f = nullptr; fopen_s(&f, "build\\ss_profile.txt", "a");
        if (f) { fprintf(f, "[ss] supersample-flagged %d / %d px = %.1f%%  (sub=%d -> +%d subpx each)\n",
                         mix_cnt, _w * _h, 100.0 * mix_cnt / (_w * _h), _sub, _sub * _sub - 1); fclose(f); }
    }    
    if (_flag_halt) return;
    if (method == 0) {
        Float c0_re_f = mpf_get_ld(_c0_re);
        Float c0_im_f = mpf_get_ld(_c0_im);
        Float dx_f = mpf_get_ld(_dx);
        Float dy_f = mpf_get_ld(_dy);
        static int simd0 = -1;
        if (simd0 < 0) { const char* e = getenv("MANDEL_SIMD"); simd0 = e ? atoi(e) : 1; }
        const bool ede = (c_method & ColoringMethod::EXTERIOR_DIST_EST) != 0;
        const bool simd0on = simd0 && !(c_method & ColoringMethod::STRIPE_AVERAGE);   // shallow SIMD has no SAC
        const int nsub = _sub * _sub;
#pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < v.size(); ++i) {
            if (_flag_halt) continue;
            if (simd0on && nsub <= 128) {
                // gather this flagged pixel's sub^2-1 subpixel c-coords, solve 4-wide
                double cre[128], cim[128]; float out[128];
                std::array<int, 4> arrs[128];
                int cnt = 0;
                for (int xi = -_sub / 2; xi <= _sub / 2; ++xi)
                    for (int yi = -_sub / 2; yi <= _sub / 2; ++yi) {
                        if (xi == 0 && yi == 0) continue;
                        cre[cnt] = c0_re_f + dx_f * v[i][1] + dx_f * xi / _sub;
                        cim[cnt] = c0_im_f + dy_f * v[i][0] + dy_f * yi / _sub;
                        arrs[cnt] = { v[i][0], v[i][1], yi, xi };
                        ++cnt;
                    }
                solveShallowSimdList(cre, cim, cnt, out, mxit, c_method);
                for (int k = 0; k < cnt; ++k) {
                    double it = out[k];
                    if (ede && it >= 0) it /= dx_f;
                    setPixel(arrs[k], it);
                }
            } else {
                for (int xi = -_sub / 2; xi <= _sub / 2; ++xi) {
                    for (int yi = -_sub / 2; yi <= _sub / 2; ++yi) {
                        if (xi == 0 && yi == 0) continue;
                        std::array<int, 4> arr = { v[i][0], v[i][1], yi, xi };
                        double iteration = floatPointCompute(c0_re_f + dx_f * v[i][1] + dx_f * xi / _sub, c0_im_f + dy_f * v[i][0] + dy_f * yi / _sub, mxit, c_method);
                        if (c_method & ColoringMethod::EXTERIOR_DIST_EST) {
                            if (iteration >= 0) iteration /= dx_f;
                        }
                        setPixel(arr, iteration);
                    }
                }
            }
        }
    }
    else if (method == 1) {
        // The first stepParallel uses the base reference's BLA table, which is
        // still valid (its dcmax covers the whole image incl. subpixels). But once
        // a glitched subpixel is re-referenced by createRef, that table is stale
        // (coefficients belong to the base reference) and rebuilding it per
        // re-reference costs more than the skips save. So keep BLA for the first
        // pass and disable it for the (few) refinement passes -- correct by
        // construction. BLA stays on for the base pass regardless.
        bool saved_bla = _use_bla;
        while (!s.empty()) {
            // std::cout << ref_it << " === " << _SA_it << ' ' << _SA_order << ' ';
            if (_flag_halt) { _use_bla = saved_bla; return; }
            stepParallel(s, ref_it, mxit, c_method);
            if (_flag_halt) { _use_bla = saved_bla; return; }
            if (s.empty()) break;
            ref_it = createRef(s, mxit, mxit, false, c_method);
            _use_bla = false;   // base BLA table is now stale for the new reference
        }
        _use_bla = saved_bla;
    }
    progressSet(1.0);
}


void Mandel::setPrecision(int precision) {
    mpf_set_prec(_c0_re, precision);
    mpf_set_prec(_c0_im, precision);
    mpf_set_prec(_dx, precision);
    mpf_set_prec(_dy, precision);
    mpf_set_prec(_scale, precision);
    mpf_set_prec(_ref_z_re, precision);
    mpf_set_prec(_ref_z_im, precision);

    mpf_set_prec(_t1, precision);
    mpf_set_prec(_t2, precision);
    for (int i = 0; i < 2; ++i) {
        mpf_set_prec(_z_re[i], precision);
        mpf_set_prec(_z_im[i], precision);
    }
}

inline bool Mandel::escape(Comp& z) const { return (std::abs(z) > _ESCAPE_RADIUS); }
inline bool Mandel::escape(mpf_t& z_re, mpf_t& z_im) const {
    // auto ri = z.get_real_imag();
    mpf_t t1, t2;
    mpf_init(t1);
    mpf_init(t2);
    mpf_mul(t1, z_re, z_re);
    mpf_mul(t2, z_im, z_im);
    mpf_add(t1, t1, t2);
    double rad = mpf_get_d(t1);
    mpf_clear(t1);
    mpf_clear(t2);
    return (rad > _ESCAPE_RADIUS);
}

inline float Mandel::getEscapeTime(Comp& z, int i) const {
    double z_re = (double)z.real();
    double z_im = (double)z.imag();
    double rad = z_re * z_re + z_im * z_im;
    return (i + 1 - log(log(rad) / 2 / log(2)) / log(2));
}

inline float Mandel::getEscapeTime(mpf_t& z_re, mpf_t& z_im, int i) {
    // auto ri = z.get_real_imag();
    mpf_mul(_t1, z_re, z_re);
    mpf_mul(_t2, z_im, z_im);
    mpf_add(_t1, _t1, _t2);
    double rad = mpf_get_d(_t1);
    return (i + 1 - log(log(rad) / 2 / log(2)) / log(2));
}

void Mandel::setPixel(std::array<int, 4> p, float iteration) const {
    _iter[getIndex(p)] = iteration;
}

void Mandel::stepParallel(std::set<std::array<int, 4>>& s, int mx_ref_it, int mxit, int c_method) {
    std::vector<std::array<int, 4>> v;
    int n = s.size();
    for (auto& p : s) v.push_back(p);
    memset(_done, 0, n * sizeof(bool));
    auto glitch_p = std::unique_ptr<Float[]>(new Float[n]);

    const bool ede = (c_method & ColoringMethod::EXTERIOR_DIST_EST) != 0;
    const bool normal = (c_method & ColoringMethod::NORMAL_MAP) != 0;
    const bool de_ovl = (c_method & ColoringMethod::DE_OVERLAY) != 0;
    // The total derivative dz/dc is tracked for EDE, the normal map and DE overlay.
    const bool deriv = ede || normal || de_ovl;
    const bool sac = (c_method & ColoringMethod::STRIPE_AVERAGE) != 0;
    // Orbit traps need every orbit point, so BLA (which skips runs) is disabled.
    const bool trap = (c_method & ColoringMethod::ORBIT_TRAP) != 0;
    const bool use_bla_loop = _use_bla && !trap;
    static int simd_env = -1;
    if (simd_env < 0) { const char* e = getenv("MANDEL_SIMD"); simd_env = e ? atoi(e) : 1; }

    if (_use_floatexp) {
        // ---- deep-zoom rescaled (floatexp) path: correct past ~1e320 ----
        // Each pixel's dc is built in floatexp (dxfe/dyfe * integer offset) so it
        // never underflows, then iterated with the rescaled delta loop.
        static int fesimd = -1;
        if (fesimd < 0) { const char* e = getenv("MANDEL_FESIMD"); fesimd = e ? atoi(e) : 1; }
        // The rescaled step vectorises cleanly (4-wide, byte-identical output) but
        // the deep loop is memory/latency-bound on the per-iteration reference
        // gather, so SIMD only ~matches scalar (same as the double-path solveSimd4).
        // BLA already skips most steps and is checked per-lane scalar every
        // iteration, so keep the default (BLA) path scalar and use SIMD only when
        // BLA is off (where the quadratic step is the whole cost).
        const bool feSimdOn = fesimd && !_use_bla && !sac && !deriv && !trap;

        std::vector<FloatExp> Dcr(n), Dci(n);
        std::vector<float> val(n);
        std::vector<float> nrm((normal || de_ovl) ? n : 0);
#pragma omp parallel for schedule(dynamic, 64)
        for (int i = 0; i < (int)v.size(); ++i) {
            if (_flag_halt) continue;
            auto arr = v[i];
            glitch_p[i] = -1;
            double px = arr[1] + (double)arr[3] / _sub;
            double py = arr[0] + (double)arr[2] / _sub;
            Dcr[i] = fe_mul_d(_dxfe, px - _ref_x);
            Dci[i] = fe_mul_d(_dyfe, py - _ref_y);
        }
        // Write each finished pixel to _iter as soon as it is computed (not in a
        // separate pass) so the GUI shows the deep render filling in progressively.
        // Rare cutoff-sensitive frames resolve uncertain pixels with the GMP oracle.
        auto finalize = [&](int i) {
            auto arr = v[i];
            float value = val[i];
            if (_fe_cutoff_sensitive && (value < 0 || value > mxit - 16)) {
                mpf_t cre, cim, t;
                mpf_init_set(cre, _c0_re); mpf_init_set(cim, _c0_im); mpf_init(t);
                mpf_mul_ui(t, _dx, arr[1]); mpf_add(cre, cre, t);
                mpf_mul_ui(t, _dx, abs(arr[3]));
                if (arr[3] < 0) mpf_neg(t, t);
                mpf_div_ui(t, t, _sub); mpf_add(cre, cre, t);
                mpf_mul_ui(t, _dy, arr[0]); mpf_add(cim, cim, t);
                mpf_mul_ui(t, _dy, abs(arr[2]));
                if (arr[2] < 0) mpf_neg(t, t);
                mpf_div_ui(t, t, _sub); mpf_add(cim, cim, t);
                value = accuratePointCompute(cre, cim, mxit, c_method);
                mpf_clears(cre, cim, t, (mpf_ptr)0);
#pragma omp atomic
                ++g_fe_fallback;
            }
            setPixel(arr, value);
            if ((normal || de_ovl) && _normal) _normal[getIndex(arr)] = nrm[i];
            markDone(i);
        };
        if (feSimdOn) {
#pragma omp parallel for schedule(dynamic, 1)
            for (int gg = 0; gg < n; gg += 4) {
                if (_flag_halt) continue;
                int lanes = n - gg < 4 ? n - gg : 4;
                solveRescaledSimd4(Dcr.data(), Dci.data(), gg, lanes, mx_ref_it, mxit, c_method, val.data());
                for (int l = 0; l < lanes; ++l) finalize(gg + l);
            }
        } else {
#pragma omp parallel for schedule(dynamic, 1)
            for (int i = 0; i < (int)v.size(); ++i) {
                if (_flag_halt) continue;
                val[i] = pixelRescaled(Dcr[i], Dci[i], mx_ref_it, mxit, c_method, (normal || de_ovl) ? &nrm[i] : nullptr);
                finalize(i);
            }
        }
        // Zhuoran rebasing (rebase whenever |z| < |dz|) keeps the floatexp path
        // glitch-free, so every pixel is resolved in one pass; drop them all.
        { int i = 0;
          for (auto it = s.begin(); it != s.end(); ++i)
              if (_done[i]) it = s.erase(it); else ++it; }
        return;
    }

    // solveSimd4 now mirrors the scalar BLA + interior-detection loop per lane, so
    // SIMD is no longer gated off when BLA / interior are on. It stays scalar only
    // for the paths solveSimd4 does not implement (EDE/normal derivative, SAC, trap).
    // The compute-bound gate (`bla_idle`) keeps SIMD off on BLA-effective deep views
    // where BLA already skips most steps (memory-bound; SIMD ~1x or slightly worse).
    // Prefer the coarse-probe MEASUREMENT (actual skip fraction); before it is taken
    // (the probe pass itself) fall back to the rmax<dcmax proxy. `_ref_period == 0`:
    // solveSimd4 has no periodic (mod-p, nucleus-fraction) support, so a periodic
    // reference must stay on the scalar path regardless of the compute-bound gate.
    const bool bla_idle = _simd_measured ? _simd_bla_idle : (!_use_bla || (_bla_rmax2 < _bla_dcmax2));
    static int simd_force = -1;
    if (simd_force < 0) { const char* e = getenv("MANDEL_SIMD_ALL"); simd_force = e ? atoi(e) : 0; }
    const bool use_simd = simd_env && !deriv && !sac && !trap && _ref_period == 0 && (bla_idle || simd_force);

    if (use_simd) {
        // ---- setup: per-pixel dc + series-approximation initial delta ----
        std::vector<double> Adcr(n), Adci(n), Adzr(n), Adzi(n);
#pragma omp parallel for schedule(dynamic, 64)
        for (int i = 0; i < (int)v.size(); ++i) {
            if (_flag_halt) continue;
            auto arr = v[i];
            mpf_t t1, t2;
            mpf_init(t1);
            mpf_init(t2);
            glitch_p[i] = -1;
            int xpix = _sub * (arr[1] - _ref[1]) + (arr[3] - _ref[3]);
            int ypix = _sub * (arr[0] - _ref[0]) + (arr[2] - _ref[2]);
            mpf_mul_ui(t1, _dx, abs(xpix));
            mpf_mul_ui(t2, _dy, abs(ypix));
            if (xpix < 0) mpf_neg(t1, t1);
            if (ypix < 0) mpf_neg(t2, t2);
            if (_sub > 1) {
                mpf_div_ui(t1, t1, _sub);
                mpf_div_ui(t2, t2, _sub);
            }
            Comp dc{ mpf_get_ld(t1), mpf_get_ld(t2) };
            mpf_clear(t1);
            mpf_clear(t2);
            Comp dz = { 0 };
            for (int x = _SA_order; x >= 0; --x) {
                dz += _Adf_old[x];
                dz *= dc;
                dz /= _SA_delta;
            }
            Adcr[i] = dc.real(); Adci[i] = dc.imag();
            Adzr[i] = dz.real(); Adzi[i] = dz.imag();
        }
        // ---- solve: AVX2, 4 pixels per group ----
#pragma omp parallel for schedule(dynamic, 1)
        for (int g = 0; g < n; g += 4) {
            if (_flag_halt) continue;
            int lanes = n - g < 4 ? n - g : 4;
            solveSimd4(v.data(), g, lanes, Adcr.data(), Adci.data(), Adzr.data(), Adzi.data(),
                       mx_ref_it, mxit, glitch_p.get());
        }
    } else {
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < v.size(); ++i) {
        if (_flag_halt) continue;
        auto arr = v[i];
        mpf_t t1, t2;
        mpf_init(t1);
        mpf_init(t2);
        glitch_p[i] = -1;
        // Comp dc = static_cast<Comp>(_dx * (v[i][1] - _ref[1]) + _dx / _sub * (v[i][3] - _ref[3]) + _dy * (v[i][0] - _ref[0]) + _dy / _sub * (v[i][2] - _ref[2]));
        int xpix = _sub * (arr[1] - _ref[1]) + (arr[3] - _ref[3]);
        int ypix = _sub * (arr[0] - _ref[0]) + (arr[2] - _ref[2]);
        mpf_mul_ui(t1, _dx, abs(xpix));
        mpf_mul_ui(t2, _dy, abs(ypix));
        if (xpix < 0) mpf_neg(t1, t1);
        if (ypix < 0) mpf_neg(t2, t2);
        if (_sub > 1) {
            mpf_div_ui(t1, t1, _sub);
            mpf_div_ui(t2, t2, _sub);
        }
        Comp dc{ mpf_get_ld(t1), mpf_get_ld(t2) };
        mpf_clear(t1);
        mpf_clear(t2);
        // Periodic reference: the integer-pixel dc above is relative to the nearest
        // pixel _ref; subtract the nucleus's sub-pixel offset so dc = c_pixel - nucleus.
        if (_ref_period) { dc = Comp{ dc.real() - _ref_frac_re, dc.imag() - _ref_frac_im }; }
        Float dxf = mpf_get_ld(_dx);
        Comp dz = { 0 };
        Comp dd = { 0 };
        for (int x = _SA_order; x >= 0; --x) {
            dz += _Adf_old[x];
            dz *= dc;
            dz /= _SA_delta;
        }
        if (deriv) {
            for (int x = _SA_order; x >= 0; --x) {
                dd += _Bdf_old[x];
                dd *= dc;
                dd /= _SA_delta;
            }
        }
        int j = _SA_it + 1;

        Float dzr = dz.real();
        Float dzi = dz.imag();
        Float dcr = dc.real();
        Float dci = dc.imag();
        Float ddr = dd.real();
        Float ddi = dd.imag();
        Float tmp;
        Float zr, zi, dr, di;
        Float dzr2, dzi2, zrad;


        int k = j;
        // Periodic reference: read the orbit modulo the period. rk mirrors k and
        // rkm1 mirrors k-1, maintained incrementally so the non-periodic path
        // (per == 0, rk == k, rkm1 == k-1 throughout) stays byte-identical.
        const int per = _ref_period;
        int rk = k, rkm1 = k - 1;
        Float m = 1e100;
        // Brent-style periodicity detection (interior): compare the full pixel
        // value z=(zr,zi) to a saved "tortoise" whose position advances at
        // geometric intervals; a close return means the orbit has entered an
        // attracting cycle, so the pixel is interior. Works for any period.
        Float zsr = 1e30, zsi = 1e30;   // sentinel: no match before first save
        int save_j = j, period_win = 1;
        // convergence confirmation state (0 = searching for a candidate cycle)
        int conf_P = 0, conf_next = 0, conf_count = 0, conf_giveup = 0;
        Float conf_D2 = 0, conf_zr = 0, conf_zi = 0;
        // Stripe Average Coloring, averaged over the last sacWindow() iterations
        // (the escaping tail) so it stays rich at deep zoom and pan-invariant. A
        // BLA skip breaks the tail, so it resets the window.
        SacAccum sacc; if (sac) sacc.init(sacWindow());
        TrapAccum trapc;
        // Exact state-repetition interior detector (MANDEL_INT_REP, default on). The
        // full pixel state is (dzr, dzi, k); it is deterministic, so if it EXACTLY
        // repeats, the orbit loops forever -> provably interior (an escaping orbit
        // never repeats its state). Zero false positives by construction;
        // BLA-compatible (BLA is deterministic). Brent-style: save a state, compare at
        // geometric intervals. Gated on _use_interior (auto-off for exterior-only
        // views -> no overhead) and !trap (trap's interior value depends on the orbit
        // accumulator).
        const bool int_rep = g_int_rep > 0 && _use_interior && !trap;
        double sdzr = 1e300, sdzi = 1e300; int sk = -1; int save_j2 = j, twin2 = 1;
        while (j < mxit) {
            if (_flag_halt) break;
            if (int_rep) {
                if (k == sk && dzr == sdzr && dzi == sdzi) {   // exact repeat -> interior
                    setPixel(arr, -2); markDone(i);
                    break;
                }
                if (j - save_j2 >= twin2) { sdzr = dzr; sdzi = dzi; sk = k; save_j2 = j; twin2 += twin2; }
            }
            if (use_bla_loop && dzr * dzr + dzi * dzi < _bla_rmax2) {
                // Only attempt BLA when dz is small enough that some level could be
                // valid -- avoids the per-iteration lookup cost when dz is large.
                // BLA start index s = k-1 (loop-top invariant: dz is at ref index k-1).
                int skip = tryBLA(rkm1, dzr, dzi, ddr, ddi, dcr, dci, deriv,
                                  _ESCAPE_RADIUS * _ESCAPE_RADIUS, mx_ref_it);
                if (skip > 0) {
                    int tid = omp_get_thread_num() & 63;
                    g_bla_stat[tid][0] += skip; ++g_bla_stat[tid][1];
                    k += skip; j += skip;
                    rk += skip; if (per && rk >= per) rk -= per; rkm1 = (per && rk == 0) ? per - 1 : rk - 1;
                    if (sac) sacc.reset_window();   // the jump breaks the tail window
                    // A BLA skip jumps j forward, making the tortoise/hare periodicity
                    // detector stale, so restart it. (The exact-state interior detector
                    // below is unaffected -- it compares the full (dz,k) state, which is
                    // deterministic across skips.)
                    zsr = 1e30; zsi = 1e30; save_j = j; period_win = 1;
                    conf_P = 0; conf_giveup = 0;
                    continue;
                }
                ++g_bla_stat[omp_get_thread_num() & 63][2];
            }
            // dd = 2 * (dd*z + dz*(d+dd))
            if (deriv) {
                if (k == 0) {
                    tmp = (dzr * (1 + ddr) - dzi * (1 + ddi)) * 2;
                    ddi = (dzr * (1 + ddi) + dzi * (1 + ddr)) * 2;
                } else {
                    tmp = (ddr * _zfr[k - 1] - ddi * _zfi[k - 1] + dzr * (_dfr[k - 1] + ddr) - dzi * (_dfi[k - 1] + ddi)) * 2;
                    ddi = (ddr * _zfi[k - 1] + ddi * _zfr[k - 1] + dzr * (_dfi[k - 1] + ddi) + dzi * (_dfr[k - 1] + ddr)) * 2;
                }
                ddr = tmp;
                dr = ddr + _dfr[k];
                di = ddi + _dfi[k];
            }
            // dz = dz^2 + 2*dz*z + dc
            dzr2 = dzr * dzr;
            dzi2 = dzi * dzi;
            if (k == 0) {
                tmp = dzr2 - dzi2 + dcr;
                dzi = dzr * dzi * 2 + dci;
            }
            else {
                tmp = (dzr * _zfr[rkm1] - dzi * _zfi[rkm1]) * 2 + dzr2 - dzi2 + dcr;
                dzi = (dzr * _zfi[rkm1] + dzi * _zfr[rkm1]) * 2 + dzr * dzi * 2 + dci;
            }
            dzr = tmp;
            zr = dzr + _zfr[rk];
            zi = dzi + _zfi[rk];
            zrad = zr * zr + zi * zi;

            if (sac) sacc.push(zr, zi);
            if (trap) trapc.push(zr, zi);

            ++k;
            rkm1 = rk; if (per) { if (++rk == per) rk = 0; } else rk = k;
            if (zrad > _ESCAPE_RADIUS * _ESCAPE_RADIUS) {
                
                if (c_method & ColoringMethod::EXTERIOR_DIST_EST) {
                    setPixel(arr, sqrt(zrad) / dxf * log(zrad) / sqrt(dr * dr + di * di));
                } else if (normal) {
                    // base colour = smooth value; normal angle = arg(z) - arg(dz/dc).
                    setPixel(arr, j + 1 - log(log((double)zrad) / 2 / log(2)) / log(2));
                    if (_normal) _normal[getIndex(arr)] = (float)(atan2((double)zi, (double)zr) - atan2((double)di, (double)dr));
                } else if (de_ovl) {
                    // base colour = smooth value; overlay = pixel-normalised DE.
                    setPixel(arr, j + 1 - log(log((double)zrad) / 2 / log(2)) / log(2));
                    if (_normal) _normal[getIndex(arr)] = (float)(sqrt(zrad) / dxf * log(zrad) / sqrt(dr * dr + di * di));
                } else if (trap) {
                    setPixel(arr, trapc.value(j + 1 - log(log((double)zrad) / 2 / log(2)) / log(2)));
                } else if (sac) {
                    setPixel(arr, sacc.value((double)zrad, (double)_ESCAPE_RADIUS));
                } else {
                    setPixel(arr, j + 1 - log(log((double)zrad) / 2 / log(2)) / log(2));
                }
                markDone(i);
                break;
            }
            if (_use_interior) {
                if (conf_P > 0) {
                    // Confirm a candidate cycle over several periods. Each period must
                    // (a) return close to the (updated) anchor -- a real cycle returns
                    // every period, a coincidental approach does not -- and (b) have
                    // multiplier |dz|^2 product = prod(4*zrad) < 1 (attracting). Only a
                    // genuine attracting cycle sustains both => interior.
                    conf_D2 *= 4 * zrad;
                    if (conf_D2 > 1e18) conf_D2 = 1e18;
                    if (j >= conf_next) {
                        Float pr = zr - conf_zr, pi = zi - conf_zi;
                        bool ret = (pr * pr + pi * pi < _interior_eps2 * zrad);
                        if (ret && conf_D2 < 1.0) {
                            if (++conf_count >= _interior_confirm) { setPixel(arr, trap ? trapc.value(0.0) : -2); markDone(i); break; }
                            conf_zr = zr; conf_zi = zi; conf_D2 = 1; conf_next = j + conf_P;
                        } else {
                            conf_P = 0; ++conf_giveup;            // not a sustained attracting cycle
                        }
                    }
                } else if ((j & 15) == 0 && conf_giveup < 3) {
                    // Search for a candidate cycle only every 16 iterations, and give up
                    // after a few failed confirmations, so escaping (exterior) pixels pay
                    // almost nothing. A close (relative) return starts a confirmation.
                    Float pr = zr - zsr, pi = zi - zsi;
                    if (pr * pr + pi * pi < _interior_eps2 * zrad) {
                        conf_P = j - save_j; if (conf_P < 1) conf_P = 1;
                        conf_zr = zr; conf_zi = zi; conf_D2 = 1; conf_next = j + conf_P; conf_count = 0;
                    }
                    if (j - save_j >= period_win) { zsr = zr; zsi = zi; save_j = j; period_win += period_win; }
                }
            }
            // ** Zhuoran method: rebase to original reference with index = 0.
            // Periodic: wrap instead of the ref-end rebase (rk already wraps).
            if (zrad < dzr * dzr + dzi * dzi || (!per && k == mx_ref_it)) {
                if ((dzr * dzr + dzi * dzi) / zrad > 10000000) {
                    // Significant magnitude change causes precision loss.
                    glitch_p[i] = zrad / _z_m3[per ? j % per : j];
                    break;
                }
                dzr = zr;
                dzi = zi;
                if (deriv) {   // EDE, normal-map AND DE-overlay all track dz/dc; reset on rebase
                    ddr = dr - 1;
                    ddi = di;
                }
                k = 0; rk = 0; rkm1 = per ? per - 1 : -1;
            }
            if (!per && k == mx_ref_it) {
                glitch_p[i] = (double)zrad / _z_m3[j];
                break;
            }
            ++j;
        }
        if (j >= mxit) {
            setPixel(arr, trap ? trapc.value(0.0) : -2);
            markDone(i);
        }
    }
    }   // end else (scalar path)
    if (_flag_halt) return;
    // get new reference with minimum |z + dz|
    Float min_orbit = 1e100;
    for (int i = 0; i < n; ++i) {
        if (glitch_p[i] >= 0 && glitch_p[i] < min_orbit) {
            min_orbit = glitch_p[i];
            _new_ref = v[i];
        }
    }

    int i = 0;
    int remove_cnt = 0;
    for (auto iter = s.begin(); iter != s.end(); ++i) {
        if (_done[i]) {
            iter = s.erase(iter);
            ++remove_cnt;
        }
        else {
            ++iter;
        }
    }

    // printf("%.2f%% %d\n", 100.0 * remove_cnt / (_w * _h), remove_cnt);
}

// Build the BLA table from the (low-precision) reference orbit in _zfr/_zfi.
// dc-independent: uses the largest pixel offset |_SA_delta| as the |Δc| bound so
// the radii are valid for every pixel. Level 0 = single steps at ref index
// 1..reflen-1; higher levels merge neighbour pairs (skip 2^p).
void Mandel::buildBLA(int reflen, bool ede) {
    _bla.clear();
    _blaD.clear();
    _bla_rmax2 = 0.0;
    if (reflen < 3) return;
    double eps = _bla_eps > 0 ? _bla_eps : ldexp(1.0, -53);   // negligible vs double rounding
    // Max |dc| over ALL pixels relative to the ACTUAL reference position. The old
    // bound |_SA_delta| assumes a centre reference; after glitch re-referencing the
    // reference is off-centre, so edge pixels have larger |dc| and a centre-based
    // bound makes the radii too loose -> BLA skips over an escape (flake@1e157
    // misclassified corner pixels). Farthest image corner from (_ref_x,_ref_y),
    // +0.5 px for the sub-pixel grid extent.
    double dxf = std::fabs(mpf_get_ld(_dx)), dyf = std::fabs(mpf_get_ld(_dy));
    double mxx = std::max(_ref_x, (double)(_w - 1) - _ref_x) + 0.5;
    double mxy = std::max(_ref_y, (double)(_h - 1) - _ref_y) + 0.5;
    double dcmax = std::sqrt(dxf * dxf * mxx * mxx + dyf * dyf * mxy * mxy);
    _bla_dcmax2 = dcmax * dcmax;   // frame max |dc|^2, for the SIMD compute-bound gate

    std::vector<BLAEntry> lvl0;
    std::vector<BLADeriv> lvl0D;
    lvl0.reserve(reflen - 1);
    if (ede) lvl0D.reserve(reflen - 1);
    for (int s = 1; s < reflen; ++s) {
        double zr = _zfr[s], zi = _zfi[s];
        double Zmag = sqrt(zr * zr + zi * zi);
        double ar = 2.0 * zr, ai = 2.0 * zi;             // A = J_f(Z) = 2 Z_s, |B| = 1
        double Amag = 2.0 * Zmag;
        // mathr / Zhuoran single-step BLA validity radius:
        //   r = max(0, eps * (|Z| - max|dc|) / (|J_f(Z)| + 1))
        double R = eps * (Zmag - dcmax) / (Amag + 1.0);
        double r2 = R > 0 ? R * R : 0.0;
        if (r2 > _bla_rmax2) _bla_rmax2 = r2;
        lvl0.push_back({ ar, ai, 1.0, 0.0, r2, 1 });
        // Under EDE also carry the derivative-delta coupling C = 2 D_s (dd -> C dz).
        if (ede) lvl0D.push_back({ 2.0 * _dfr[s], 2.0 * _dfi[s], 0.0, 0.0 });
    }
    _bla.push_back(std::move(lvl0));
    if (ede) _blaD.push_back(std::move(lvl0D));

    while (_bla.back().size() > 1) {
        const std::vector<BLAEntry>& prev = _bla.back();
        const std::vector<BLADeriv>* prevD = ede ? &_blaD.back() : nullptr;
        std::vector<BLAEntry> nxt;
        std::vector<BLADeriv> nxtD;
        nxt.reserve(prev.size() / 2);
        if (ede) nxtD.reserve(prev.size() / 2);
        for (size_t i = 0; i + 1 < prev.size(); i += 2) {
            const BLAEntry& x = prev[i];
            const BLAEntry& y = prev[i + 1];
            BLAEntry z;
            z.ar = y.ar * x.ar - y.ai * x.ai;                 // A_z = A_y A_x
            z.ai = y.ar * x.ai + y.ai * x.ar;
            z.br = y.ar * x.br - y.ai * x.bi + y.br;          // B_z = A_y B_x + B_y
            z.bi = y.ar * x.bi + y.ai * x.br + y.bi;
            if (ede) {
                const BLADeriv& xd = (*prevD)[i];
                const BLADeriv& yd = (*prevD)[i + 1];
                BLADeriv zd;
                // Derivative couplings: C_z = C_y A_x + A_y C_x, E_z = C_y B_x + A_y E_x + E_y
                zd.cr = (yd.cr * x.ar - yd.ci * x.ai) + (y.ar * xd.cr - y.ai * xd.ci);
                zd.ci = (yd.cr * x.ai + yd.ci * x.ar) + (y.ar * xd.ci + y.ai * xd.cr);
                zd.er = (yd.cr * x.br - yd.ci * x.bi) + (y.ar * xd.er - y.ai * xd.ei) + yd.er;
                zd.ei = (yd.cr * x.bi + yd.ci * x.br) + (y.ar * xd.ei + y.ai * xd.er) + yd.ei;
                nxtD.push_back(zd);
            }
            double Rx = sqrt(x.r2), Ry = sqrt(y.r2);
            double Axmag = sqrt(x.ar * x.ar + x.ai * x.ai);
            double Bxmag = sqrt(x.br * x.br + x.bi * x.bi);
            double cand = Axmag > 0 ? (Ry - Bxmag * dcmax) / Axmag : 0.0;
            double R = std::min(Rx, cand);
            z.r2 = R > 0 ? R * R : 0.0;
            z.l = x.l + y.l;
            nxt.push_back(z);
        }
        _bla.push_back(std::move(nxt));   // odd leftover dropped; still available lower
        if (ede) _blaD.push_back(std::move(nxtD));
    }
    if (getenv("MANDEL_PROFILE"))
        fprintf(stderr, "  [profile] buildBLA: reflen=%d levels=%zu dcmax=%.3e rmax=%.3e |Z1|=%.3e SA_it=%d SA_order=%d\n",
                reflen, _bla.size(), dcmax, sqrt(_bla_rmax2), sqrt(_zfr[1] * _zfr[1] + _zfi[1] * _zfi[1]), _SA_it, _SA_order);
}

// Largest valid BLA starting at reference index s. Levels have skip 2^p and start
// at ref index 1 + i*2^p, so s is a valid start at level p iff (s-1) is a multiple
// of 2^p. Returns skip applied (0 if none), updating dz on success.
int Mandel::tryBLA(int s, double& dzr, double& dzi, double& ddr, double& ddi,
                   double dcr, double dci, bool ede, double ESC2, int mx_ref_it) const {
    if (s < 1 || _bla.empty()) return 0;
    double zmag2 = dzr * dzr + dzi * dzi;
    int maxp = (int)_bla.size() - 1;
    // Levels start at ref index 1 + i*2^p, so s is aligned at exactly levels
    // 0..tz where tz = trailing zeros of (s-1); start the search there (all these
    // levels are aligned, so no per-level alignment test is needed).
    int startp = (s == 1) ? maxp : (int)_tzcnt_u32((unsigned)(s - 1));
    if (startp > maxp) startp = maxp;
    if (startp < _bla_minlevel) return 0;                   // not aligned enough for a worthwhile skip
    for (int p = startp; p >= _bla_minlevel; --p) {
        int i = (s - 1) >> p;
        if (i >= (int)_bla[p].size()) continue;             // past the right edge at this level
        const BLAEntry& b = _bla[p][i];
        if (b.r2 <= 0 || zmag2 >= b.r2) continue;           // radius too small -> smaller skip
        int land = s + b.l;                                 // dz lands at ref index `land`
        if (land >= mx_ref_it - 1) continue;                // leave one step of headroom before the ref-end rebase
        double nzr = b.ar * dzr - b.ai * dzi + b.br * dcr - b.bi * dci;
        double nzi = b.ar * dzi + b.ai * dzr + b.br * dci + b.bi * dcr;
        if (!g_bla_noescape) {
            // Never let a BLA skip overshoot the pixel's escape. (Previously this
            // was skipped for bounded references assuming dz stays negligibly
            // small, but at a usable eps dz is not tiny and a skip can jump past
            // an escape, producing wrong iteration counts.)
            double fr = _zfr[land] + nzr, fi = _zfi[land] + nzi;
            if (fr * fr + fi * fi > ESC2) return 0;
        }
        if (ede) {
            // dd -> C*dz + A*dd + E*dc  (derivative delta carried through the skip)
            const BLADeriv& bd = _blaD[p][i];
            double nddr = (bd.cr * dzr - bd.ci * dzi) + (b.ar * ddr - b.ai * ddi) + (bd.er * dcr - bd.ei * dci);
            double nddi = (bd.cr * dzi + bd.ci * dzr) + (b.ar * ddi + b.ai * ddr) + (bd.er * dci + bd.ei * dcr);
            ddr = nddr; ddi = nddi;
        }
        dzr = nzr; dzi = nzi;
        return b.l;
    }
    return 0;
}

// Floatexp counterpart of tryBLA for the rescaled kernel. dz is carried as
// dz = S*(wr,wi); the linear map dz -> A*dz + B*dc and the |dz|^2 < r2 test are
// done in floatexp so they stay correct past double underflow, then the result
// is re-rescaled back into (S,wr,wi). Reuses the same BLA table (built in double)
// -- the coefficients A,B are O(1) and representable in double.
int Mandel::tryBLAfe(int s, FloatExp& S, FloatExp S2, double& wr, double& wi,
                     double dr, double di, double ESC2, int mx_ref_it,
                     double* outAB) const {
    if (s < 1 || _bla.empty()) return 0;
    int maxp = (int)_bla.size() - 1;
    int startp = (s == 1) ? maxp : (int)_tzcnt_u32((unsigned)(s - 1));
    if (startp > maxp) startp = maxp;
    if (startp < _bla_minlevel) return 0;   // alignment too low: cheap out before dz2
    // |dz|^2 = S^2 (wr^2 + wi^2); S2 = S*S is cached by the caller (constant
    // between rescales/rebases, so ~1 fe_mul saved per call).
    FloatExp dz2 = fe_mul(S2, fe_from(wr * wr + wi * wi));
    for (int p = startp; p >= _bla_minlevel; --p) {
        int i = (s - 1) >> p;
        if (i >= (int)_bla[p].size()) continue;
        const BLAEntry& b = _bla[p][i];
        if (b.r2 <= 0 || !fe_abs_less(dz2, fe_from(b.r2))) continue;   // radius too small
        int land = s + b.l;
        if (land >= mx_ref_it - 1) continue;
        // The BLA map nz = A*dz + B*dc shares the delta scale S: with dz = S*w and
        // dc = S*d (d = dc/S, kept in double by the caller), nz = S*(A*w + B*d). A,B
        // are the double table coefficients and w,d are O(1), so the new unit delta
        // w' = A*w + B*d is normally computed entirely in double -- no underflow (S is
        // factored out) -- replacing ~8 fe_mul_d + 6 fe_add (each a frexp/ldexp) with
        // plain double on the hot 99.9%-skip path. BUT over a long skip |A| can be
        // enormous (up to the double range), so if |w'| is large enough that w'^2 (the
        // re-rescale) or A*w itself could overflow, fall back to an exact floatexp
        // re-rescale (rare: only very long skips near escape).
        double nwr = (b.ar * wr - b.ai * wi) + (b.br * dr - b.bi * di);
        double nwi = (b.ar * wi + b.ai * wr) + (b.br * di + b.bi * dr);
        const double BIG = 1e150;   // w'^2 must stay < DBL_MAX (~1.8e308)
        const double* zfp = reinterpret_cast<const double*>(_zf);
        if (std::isfinite(nwr) && std::isfinite(nwi) && std::fabs(nwr) < BIG && std::fabs(nwi) < BIG) {
            // fast double path. Never skip past an escape: z = X_land + S*w'.
            double fr = zfp[2 * land]     + fe_to_double(fe_mul_d(S, nwr));
            double fi = zfp[2 * land + 1] + fe_to_double(fe_mul_d(S, nwi));
            if (fr * fr + fi * fi > ESC2) return 0;
            // Accept: re-rescale dz = S*w' -> keep |w| ~ 1, fold |w'| into S.
            double wm2 = nwr * nwr + nwi * nwi;
            if (wm2 == 0.0) { S = FloatExp{ 1.0, 0 }; wr = wi = 0.0; }
            else {
                double wmag = std::sqrt(wm2);
                S = fe_mul(S, fe_from(wmag));
                double inv = 1.0 / wmag;
                wr = nwr * inv; wi = nwi * inv;
            }
        } else {
            // Overflow-safe floatexp re-rescale: recompute nz = A*dz + B*dc directly in
            // floatexp (dz = S*w, dc = S*d), so the huge A never leaves the extended
            // exponent range. Matches the pre-optimization behaviour.
            g_blafe_safe[omp_get_thread_num() & 63]++;
            FloatExp dzr = fe_mul_d(S, wr), dzi = fe_mul_d(S, wi);
            FloatExp dcrfe = fe_mul_d(S, dr), dcife = fe_mul_d(S, di);
            FloatExp nzr = fe_add(fe_sub(fe_mul_d(dzr, b.ar), fe_mul_d(dzi, b.ai)),
                                  fe_sub(fe_mul_d(dcrfe, b.br), fe_mul_d(dcife, b.bi)));
            FloatExp nzi = fe_add(fe_add(fe_mul_d(dzr, b.ai), fe_mul_d(dzi, b.ar)),
                                  fe_add(fe_mul_d(dcrfe, b.bi), fe_mul_d(dcife, b.br)));
            double fr = zfp[2 * land] + fe_to_double(nzr), fi = zfp[2 * land + 1] + fe_to_double(nzi);
            if (fr * fr + fi * fi > ESC2) return 0;
            FloatExp Snew = fe_sqrt(fe_add(fe_mul(nzr, nzr), fe_mul(nzi, nzi)));
            if (Snew.m == 0.0) { S = FloatExp{ 1.0, 0 }; wr = wi = 0.0; }
            else { S = Snew; wr = fe_to_double(fe_div(nzr, S)); wi = fe_to_double(fe_div(nzi, S)); }
        }
        // Export the linear map (A, B) so the caller can carry the EDE total
        // derivative through the same skip: J -> A*J + B (see pixelRescaled).
        if (outAB) { outAB[0] = b.ar; outAB[1] = b.ai; outAB[2] = b.br; outAB[3] = b.bi; }
        return b.l;
    }
    return 0;
}

// scalar loop in the same order (no FMA contraction), so the per-pixel result
// is bit-identical to the scalar path. The reference orbit is read per lane via
// gather; escape / rebase / glitch decisions are done in a short scalar pass
// over the 4 lanes (rare events), keeping the arithmetic vectorized. BLA skips
// and periodicity-based interior detection are mirrored per lane (see below), so
// this path is a faithful vector image of the scalar stepParallel loop.
void Mandel::solveSimd4(const std::array<int, 4>* v, int g, int lanes,
                        const double* Adcr, const double* Adci,
                        const double* Adzr, const double* Adzi,
                        int mx_ref_it, int mxit, Float* glitch_p) {
    const double ESC2 = (double)_ESCAPE_RADIUS * _ESCAPE_RADIUS;
    const double LG2 = log(2.0);

    alignas(32) double dzr_[4] = { 0,0,0,0 }, dzi_[4] = { 0,0,0,0 };
    alignas(32) double dcr_[4] = { 0,0,0,0 }, dci_[4] = { 0,0,0,0 };
    int k_[4] = { 0,0,0,0 }, j_[4] = { 0,0,0,0 };
    bool act[4] = { false,false,false,false };
    // Per-lane Brent periodicity (interior) detector state, mirroring the scalar
    // stepParallel loop exactly. Only exercised when _use_interior.
    double zsr[4], zsi[4], confD2[4], confZr[4], confZi[4];
    int save_j[4], pwin[4], confP[4], confNext[4], confCount[4], confGiveup[4];
    for (int l = 0; l < lanes; ++l) {
        dzr_[l] = Adzr[g + l]; dzi_[l] = Adzi[g + l];
        dcr_[l] = Adcr[g + l]; dci_[l] = Adci[g + l];
        k_[l] = j_[l] = _SA_it + 1;
        act[l] = true;
        zsr[l] = 1e30; zsi[l] = 1e30; save_j[l] = j_[l]; pwin[l] = 1;
        confP[l] = 0; confNext[l] = 0; confCount[l] = 0; confGiveup[l] = 0;
        confD2[l] = 0.0; confZr[l] = 0.0; confZi[l] = 0.0;
        // Mirror the scalar `while (j < mxit)` entry guard: if series approximation
        // stayed valid to the iteration cap (_SA_it + 1 >= mxit), the pixel is
        // interior (-2) and must NOT enter the loop (which would gather _zfr past
        // the reference end). Rare (tiny views hugging the reference) but real.
        if (j_[l] >= mxit) { setPixel(v[g + l], -2.f); markDone(g + l); act[l] = false; }
    }

    bool anyactive = false;
    for (int l = 0; l < lanes; ++l) if (act[l]) { anyactive = true; break; }

    __m256d dzr = _mm256_load_pd(dzr_), dzi = _mm256_load_pd(dzi_);
    const __m256d dcr = _mm256_load_pd(dcr_), dci = _mm256_load_pd(dci_);
    const __m256d two = _mm256_set1_pd(2.0);
    const __m256d ESC2v = _mm256_set1_pd(ESC2);
    const __m256d rmax2v = _mm256_set1_pd(_bla_rmax2);
    const bool use_bla = _use_bla;
    const bool use_int = _use_interior;

    while (anyactive) {
        if (_flag_halt) return;

        // ---- BLA: skip a run of reference iterations for lanes whose |dz| is small
        // enough (deep-into-structure). Vector-gated: if no active lane is eligible
        // (|dz|^2 < rmax2) this is one compare + movemask and the fast path is
        // untouched. Mirrors the scalar loop's BLA block (start index s = k-1). A
        // skipping lane advances itself and sits out this round's vector step (like
        // the scalar `continue`). ----
        int blaSkipped = 0;
        if (use_bla) {
            __m256d dzmag_top = _mm256_add_pd(_mm256_mul_pd(dzr, dzr), _mm256_mul_pd(dzi, dzi));
            int elig = _mm256_movemask_pd(_mm256_cmp_pd(dzmag_top, rmax2v, _CMP_LT_OQ));
            if (elig) {
                _mm256_store_pd(dzr_, dzr); _mm256_store_pd(dzi_, dzi);
                for (int l = 0; l < lanes; ++l) {
                    if (!act[l] || !(elig & (1 << l))) continue;
                    double ddr = 0, ddi = 0;
                    int skip = tryBLA(k_[l] - 1, dzr_[l], dzi_[l], ddr, ddi,
                                      dcr_[l], dci_[l], false, ESC2, mx_ref_it);
                    int tid = omp_get_thread_num() & 63;
                    if (skip > 0) {
                        g_bla_stat[tid][0] += skip; ++g_bla_stat[tid][1];
                        k_[l] += skip; j_[l] += skip; blaSkipped |= 1 << l;
                        // a BLA skip jumps j forward, so the periodicity detector's
                        // tortoise/hare and any in-progress confirmation are stale.
                        zsr[l] = 1e30; zsi[l] = 1e30; save_j[l] = j_[l]; pwin[l] = 1;
                        confP[l] = 0; confGiveup[l] = 0;
                        if (j_[l] >= mxit) {          // skip reached the cap -> interior
                            setPixel(v[g + l], -2.f); markDone(g + l);
                            act[l] = false; blaSkipped &= ~(1 << l);
                        }
                    } else {
                        ++g_bla_stat[tid][2];
                    }
                }
                if (blaSkipped) { dzr = _mm256_load_pd(dzr_); dzi = _mm256_load_pd(dzi_); }
            }
        }

        // ---- manual gather + vector quadratic step (all lanes; skipped/inactive
        // lanes compute unused values, filtered out below). ----
        alignas(32) double zprr[4], zpri[4], zcrr[4], zcri[4];
        for (int l = 0; l < 4; ++l) {
            int kk = k_[l];
            zcrr[l] = _zfr[kk]; zcri[l] = _zfi[kk];
            if (kk) { zprr[l] = _zfr[kk - 1]; zpri[l] = _zfi[kk - 1]; }
            else    { zprr[l] = 0.0;         zpri[l] = 0.0; }
        }
        __m256d Zpr = _mm256_load_pd(zprr), Zpi = _mm256_load_pd(zpri);
        __m256d Zcr = _mm256_load_pd(zcrr), Zci = _mm256_load_pd(zcri);

        __m256d dzr2 = _mm256_mul_pd(dzr, dzr);
        __m256d dzi2 = _mm256_mul_pd(dzi, dzi);
        // dzr' = (dzr*Zpr - dzi*Zpi)*2 + dzr2 - dzi2 + dcr
        __m256d t = _mm256_sub_pd(_mm256_mul_pd(dzr, Zpr), _mm256_mul_pd(dzi, Zpi));
        t = _mm256_mul_pd(t, two);
        t = _mm256_add_pd(t, dzr2);
        t = _mm256_sub_pd(t, dzi2);
        __m256d ndzr = _mm256_add_pd(t, dcr);
        // dzi' = (dzr*Zpi + dzi*Zpr)*2 + (dzr*dzi)*2 + dci
        __m256d u = _mm256_add_pd(_mm256_mul_pd(dzr, Zpi), _mm256_mul_pd(dzi, Zpr));
        u = _mm256_mul_pd(u, two);
        __m256d w = _mm256_mul_pd(_mm256_mul_pd(dzr, dzi), two);
        u = _mm256_add_pd(u, w);
        __m256d ndzi = _mm256_add_pd(u, dci);

        // Preserve BLA-skipped lanes' deltas (they did not step this round).
        if (blaSkipped) {
            alignas(32) double ns_r[4], ns_i[4];
            _mm256_store_pd(ns_r, ndzr); _mm256_store_pd(ns_i, ndzi);
            for (int l = 0; l < lanes; ++l) if (blaSkipped & (1 << l)) { ns_r[l] = dzr_[l]; ns_i[l] = dzi_[l]; }
            ndzr = _mm256_load_pd(ns_r); ndzi = _mm256_load_pd(ns_i);
        }
        dzr = ndzr; dzi = ndzi;
        __m256d zr = _mm256_add_pd(dzr, Zcr);
        __m256d zi = _mm256_add_pd(dzi, Zci);
        __m256d zrad = _mm256_add_pd(_mm256_mul_pd(zr, zr), _mm256_mul_pd(zi, zi));
        __m256d dzmag = _mm256_add_pd(_mm256_mul_pd(dzr, dzr), _mm256_mul_pd(dzi, dzi));

        // Vector escape / rebase masks (raw; per-lane guarded below). Fused pass:
        // advance k, build refit, and fold escape/rebase/interior-due into `need`.
        int em = _mm256_movemask_pd(_mm256_cmp_pd(zrad, ESC2v, _CMP_GT_OQ));
        int rm = _mm256_movemask_pd(_mm256_cmp_pd(zrad, dzmag, _CMP_LT_OQ));
        int refit = 0, need = 0;
        for (int l = 0; l < lanes; ++l) {
            if (!act[l] || (blaSkipped & (1 << l))) continue;   // not a stepping lane
            int bit = 1 << l;
            ++k_[l];
            if (em & bit) need |= bit;
            else if (rm & bit) need |= bit;
            if (k_[l] == mx_ref_it) { refit |= bit; need |= bit; }
            // interior detector needs this iteration's z during a confirmation (every
            // iteration) or at a search probe (every 16th) -> take the per-lane path.
            if (use_int && (confP[l] > 0 || ((j_[l] & 15) == 0 && confGiveup[l] < 3))) need |= bit;
        }
        if (need == 0) {
            anyactive = false;
            for (int l = 0; l < lanes; ++l) {
                if (!act[l]) continue;
                if (blaSkipped & (1 << l)) { anyactive = true; continue; }   // skipped: re-loop
                if (++j_[l] >= mxit) { setPixel(v[g + l], -2.f); markDone(g + l); act[l] = false; }
                else anyactive = true;
            }
            continue;   // dzr/dzi already hold the updated deltas
        }

        // slow path: at least one stepping lane escaped / needs interior / rebased.
        alignas(32) double zr_[4], zi_[4], zrad_[4], dzmag_[4], dzrs[4], dzis[4];
        _mm256_store_pd(zr_, zr);     _mm256_store_pd(zi_, zi);
        _mm256_store_pd(zrad_, zrad); _mm256_store_pd(dzmag_, dzmag);
        _mm256_store_pd(dzrs, dzr);   _mm256_store_pd(dzis, dzi);

        bool changed = false;
        for (int l = 0; l < lanes; ++l) {
            if (!act[l] || (blaSkipped & (1 << l))) continue;
            if (em & (1 << l)) {                         // escape (checked first)
                setPixel(v[g + l], (float)(j_[l] + 1 - log(log(zrad_[l]) / 2 / LG2) / LG2));
                markDone(g + l); act[l] = false; continue;
            }
            // Interior (periodicity) detection, mirroring stepParallel. It only
            // affects whether/when a non-escaping pixel is classified interior
            // (value -2); exterior pixels and their smooth values are untouched.
            if (use_int) {
                double zrl = zrad_[l];
                if (confP[l] > 0) {
                    confD2[l] *= 4.0 * zrl; if (confD2[l] > 1e18) confD2[l] = 1e18;
                    if (j_[l] >= confNext[l]) {
                        double pr = zr_[l] - confZr[l], pi = zi_[l] - confZi[l];
                        if (pr * pr + pi * pi < _interior_eps2 * zrl && confD2[l] < 1.0) {
                            if (++confCount[l] >= _interior_confirm) { setPixel(v[g + l], -2.f); markDone(g + l); act[l] = false; continue; }
                            confZr[l] = zr_[l]; confZi[l] = zi_[l]; confD2[l] = 1; confNext[l] = j_[l] + confP[l];
                        } else { confP[l] = 0; ++confGiveup[l]; }
                    }
                } else if ((j_[l] & 15) == 0 && confGiveup[l] < 3) {
                    double pr = zr_[l] - zsr[l], pi = zi_[l] - zsi[l];
                    if (pr * pr + pi * pi < _interior_eps2 * zrl) {
                        confP[l] = j_[l] - save_j[l]; if (confP[l] < 1) confP[l] = 1;
                        confZr[l] = zr_[l]; confZi[l] = zi_[l]; confD2[l] = 1; confNext[l] = j_[l] + confP[l]; confCount[l] = 0;
                    }
                    if (j_[l] - save_j[l] >= pwin[l]) { zsr[l] = zr_[l]; zsi[l] = zi_[l]; save_j[l] = j_[l]; pwin[l] += pwin[l]; }
                }
            }
            if ((rm & (1 << l)) || (refit & (1 << l))) {  // Zhuoran rebase
                if (dzmag_[l] / zrad_[l] > 10000000.0) {
                    glitch_p[g + l] = zrad_[l] / _z_m3[j_[l]];
                    act[l] = false; continue;
                }
                dzrs[l] = zr_[l]; dzis[l] = zi_[l]; k_[l] = 0; changed = true;
            }
            if (k_[l] == mx_ref_it) {
                glitch_p[g + l] = zrad_[l] / _z_m3[j_[l]];
                act[l] = false; continue;
            }
            if (++j_[l] >= mxit) {
                setPixel(v[g + l], -2.f);
                markDone(g + l); act[l] = false; continue;
            }
        }
        if (changed) { dzr = _mm256_load_pd(dzrs); dzi = _mm256_load_pd(dzis); }
        anyactive = false;
        for (int l = 0; l < lanes; ++l) if (act[l]) { anyactive = true; break; }
    }
}

// Deep-zoom rescaled (z = S w) perturbation for one pixel. See header comment:
// w stays an O(1) double, the floatexp scale S carries the deep exponent, so the
// inner loop is native-double yet correct far past double's ~1e320 underflow.
float Mandel::pixelRescaled(FloatExp dcr, FloatExp dci, int mx_ref_it, int mxit, int c_method, float* normalOut) const {
    const double ESC2 = (double)_ESCAPE_RADIUS * _ESCAPE_RADIUS;
    const double LG2 = log(2.0);
    struct FeStat { long long* s; long long sk = 0, skl = 0, st = 0;
        FeStat(long long* p) : s(p) {} ~FeStat() { s[0] += skl; s[1] += sk; s[2] += st; } }
        g(g_fe_stat[omp_get_thread_num() & 63]);
    // The stored reference is _zfr[k] = X_{k+1} (X_0 = 0 is the implicit critical
    // point). Access the orbit as X_m: X_0 = 0, X_m = _zfr[m-1]. Rebasing resets to
    // the critical point m = 0 (X_0 = 0), matching the double path's k==0 case.
    // Read the orbit through the interleaved _zf (re,im adjacent) so each X_m load
    // touches one cache line instead of two separate _zfr/_zfi streams: this loop
    // is memory-bound, so halving the reference footprint is the real lever.
    const double* zfp = reinterpret_cast<const double*>(_zf);
    const int reflen = mx_ref_it + 1;
    // Periodic reference (MANDEL_PERIODIC): the orbit repeats every `per` entries,
    // so read it at rm = (m-1) mod per and never force-rebase on running off the
    // end (it wraps instead). rm is maintained incrementally so the non-periodic
    // path (per == 0, rm == m-1 throughout) stays byte-identical.
    const int per = _ref_period;
    int rm = 0;                          // reference index for X_m (mirrors m-1)

    FloatExp S = fe_sqrt(fe_add(fe_mul(dcr, dcr), fe_mul(dci, dci)));   // |dc|
    double wr, wi, dr, di;
    if (S.m == 0.0) { S = FloatExp{ 1.0, 0 }; wr = wi = dr = di = 0.0; }
    else {
        wr = fe_to_double(fe_div(dcr, S)); wi = fe_to_double(fe_div(dci, S));   // unit
        dr = wr; di = wi;                                                        // d = dc / S
    }
    double s = fe_to_double(S);
    FloatExp S2 = fe_mul(S, S);            // cached |S|^2 for tryBLAfe's radius test
    int m = 1, iter = 1;                  // dz_1 = dc at reference index m = 1

    double zsr = 1e30, zsi = 1e30; int save_iter = 1, period_win = 1;
    int conf_P = 0, conf_next = 0, conf_count = 0, conf_giveup = 0;
    double conf_D2 = 0, conf_zr = 0, conf_zi = 0;
    // Stripe Average Coloring. Averaged over the last sacWindow() iterations (the
    // escaping tail) so Feather keeps its structure at deep zoom; a BLA skip breaks
    // the tail, so it resets the window. W<=0 uses the classic full average, and
    // then restores a skip's omitted points via the reference-orbit prefix sum.
    const bool sac = (c_method & ColoringMethod::STRIPE_AVERAGE) != 0;
    SacAccum sacc; if (sac) sacc.init(sacWindow());
    const bool trap = (c_method & ColoringMethod::ORBIT_TRAP) != 0;   // BLA disabled below
    TrapAccum trapc;

    // Exterior distance estimation. The double path tracks the derivative in
    // double, but at deep zoom the true derivative |dz/dc| ~ 1/dx overflows it,
    // so carry the total derivative J = dz/dc in its own rescaled floatexp
    // (scale SJ, unit (jr,ji)) exactly like the delta w. J is rebase-invariant
    // (it is defined on the reconstructed z), advances as J' = 2 z J + 1 per
    // step, and as J -> A J + B through a BLA skip (A,B are the same linear-map
    // coefficients tryBLAfe uses for dz -- the derivative of dz = A dz + B dc).
    const bool ede = (c_method & ColoringMethod::EXTERIOR_DIST_EST) != 0;
    const bool normal = (c_method & ColoringMethod::NORMAL_MAP) != 0;
    const bool de_ovl = (c_method & ColoringMethod::DE_OVERLAY) != 0;
    const bool deriv = ede || normal || de_ovl;   // track J = dz/dc for EDE (|J|), normal (arg J), DE overlay
    FloatExp SJ{ 1.0, 0 };            // |J| scale; J_1 = d(dc)/dc = 1
    double jr = 1.0, ji = 0.0, invSJd = 1.0;   // invSJd = 1/SJ carries the "+1" seed

    // Exact state-repetition interior detector (see stepParallel). The floatexp
    // pixel state is (wr, wi, S, m); it is deterministic, so a bit-exact repeat
    // proves the orbit loops -> interior. Zero false positives by construction.
    const bool int_rep = g_int_rep > 0 && _use_interior && !trap;
    double swr = 1e300, swi = 1e300; FloatExp sS{ 1e300, 0 }; int sm = -1, save_iter2 = 1, twin2r = 1;

    while (iter < mxit) {
        if (_flag_halt) break;
        if (int_rep) {
            if (m == sm && wr == swr && wi == swi && S.m == sS.m && S.e == sS.e)
                return -2.f;                       // exact state repeat -> provably interior
            if (iter - save_iter2 >= twin2r) { swr = wr; swi = wi; sS = S; sm = m; save_iter2 = iter; twin2r += twin2r; }
        }
        // BLA: skip a run of reference iterations when the rescaled dz is small
        // enough (deep zoom -> almost always). dz is at reference index m (double
        // path convention: table index s = m-1). On a skip, dz/S/wr/wi are
        // advanced; recompute the scalar scale s and unit d = dc/S, and restart
        // the periodicity detector (a multi-iter skip leaves its state stale).
        if (_use_bla && !trap && m >= 2) {
            double ab[4];
            int skip = tryBLAfe(rm, S, S2, wr, wi, dr, di, ESC2, per ? per : reflen, deriv ? ab : nullptr);
            if (skip > 0) {
                if (sac) {
                    if (sacc.W > 0) sacc.reset_window();      // tail broken by the jump
                    else if (!_sacRefPre.empty() && m + skip < (int)_sacRefPre.size())
                        sacc.add_full(_sacRefPre[m + skip] - _sacRefPre[m], skip);   // full-avg restore
                }
                if (deriv) {
                    // J -> A*J + B, carried in floatexp (A can be large over a run).
                    FloatExp Jr = fe_mul_d(SJ, jr), Ji = fe_mul_d(SJ, ji);
                    FloatExp nJr = fe_add(fe_sub(fe_mul_d(Jr, ab[0]), fe_mul_d(Ji, ab[1])), fe_from(ab[2]));
                    FloatExp nJi = fe_add(fe_add(fe_mul_d(Jr, ab[1]), fe_mul_d(Ji, ab[0])), fe_from(ab[3]));
                    SJ = fe_sqrt(fe_add(fe_mul(nJr, nJr), fe_mul(nJi, nJi)));
                    if (SJ.m == 0.0) { SJ = FloatExp{ 1.0, 0 }; jr = ji = 0.0; invSJd = 1.0; }
                    else { jr = fe_to_double(fe_div(nJr, SJ)); ji = fe_to_double(fe_div(nJi, SJ)); invSJd = 1.0 / fe_to_double(SJ); }
                }
                m += skip; iter += skip; g.sk++; g.skl += skip;
                rm += skip; if (per && rm >= per) rm -= per;   // land < per (skip is capped)
                s = fe_to_double(S); S2 = fe_mul(S, S);
                dr = fe_to_double(fe_div(dcr, S)); di = fe_to_double(fe_div(dci, S));
                zsr = 1e30; zsi = 1e30; save_iter = iter; period_win = 1;
                conf_P = 0; conf_giveup = 0;
                continue;
            }
        }
        // w' = 2 X_m w + s w^2 + d   (dz_m -> dz_{m+1}); X_0 = 0.
        double Xr = m ? zfp[2 * rm] : 0.0, Xi = m ? zfp[2 * rm + 1] : 0.0;
        double zmr = Xr + s * wr, zmi = Xi + s * wi;   // z at index m (pre-step), for J' = 2 z J + 1
        double w2r = wr * wr - wi * wi, w2i = 2.0 * wr * wi;
        double nwr = 2.0 * (Xr * wr - Xi * wi) + s * w2r + dr;
        double nwi = 2.0 * (Xr * wi + Xi * wr) + s * w2i + di;
        wr = nwr; wi = nwi; ++m; ++iter; g.st++;
        rm++; if (per && rm == per) rm = 0;             // advance the (mod-per) ref index
        if (deriv) {
            // J_{m+1} = 2 z_m J_m + 1 (exact; z_m is the reconstructed value).
            double a = 2.0 * (zmr * jr - zmi * ji) + invSJd;   // "+1" is invSJd in the SJ scale
            double b = 2.0 * (zmr * ji + zmi * jr);
            jr = a; ji = b;
            double jm2 = jr * jr + ji * ji;
            if (jm2 > 1e16 || (jm2 > 0.0 && jm2 < 1e-16)) {   // keep |j| ~ O(1)
                FloatExp jmag = fe_sqrt(fe_from(jm2));
                SJ = fe_mul(SJ, jmag);
                double inv = 1.0 / fe_to_double(jmag);
                jr *= inv; ji *= inv;
                invSJd = 1.0 / fe_to_double(SJ);
            }
        }

        Xr = m ? zfp[2 * rm] : 0.0; Xi = m ? zfp[2 * rm + 1] : 0.0;
        double zr = Xr + s * wr, zi = Xi + s * wi;      // z = X_m + S w
        double zrad = zr * zr + zi * zi;
        if (sac) sacc.push(zr, zi);
        if (trap) trapc.push(zr, zi);
        if (zrad > ESC2) {
            if (ede) {
                // dist = |z| log|z|^2 / (dx * |z'|); |z'| = SJ * |j| in floatexp.
                double num = sqrt(zrad) * log(zrad);
                FloatExp Jmag = fe_mul_d(SJ, sqrt(jr * jr + ji * ji));
                return (float)fe_to_double(fe_div(fe_from(num), fe_mul(_dxfe, Jmag)));
            }
            if (normal) {
                // normal angle = arg(z) - arg(J); SJ > 0 so arg(J) = atan2(ji, jr).
                if (normalOut) *normalOut = (float)(atan2(zi, zr) - atan2(ji, jr));
                return (float)((double)iter - log(log(zrad) / 2.0 / LG2) / LG2);
            }
            if (de_ovl) {
                // base colour = smooth value; overlay = pixel-normalised DE.
                double num = sqrt(zrad) * log(zrad);
                FloatExp Jmag = fe_mul_d(SJ, sqrt(jr * jr + ji * ji));
                if (normalOut) *normalOut = (float)fe_to_double(fe_div(fe_from(num), fe_mul(_dxfe, Jmag)));
                return (float)((double)iter - log(log(zrad) / 2.0 / LG2) / LG2);
            }
            if (trap) {
                return trapc.value((double)iter - log(log(zrad) / 2.0 / LG2) / LG2);
            }
            if (!sac)
                return (float)((double)iter - log(log(zrad) / 2.0 / LG2) / LG2);
            return sacc.value(zrad, (double)_ESCAPE_RADIUS);
        }

        if (_use_interior) {
            if (conf_P > 0) {
                conf_D2 *= 4.0 * zrad; if (conf_D2 > 1e18) conf_D2 = 1e18;
                if (iter >= conf_next) {
                    double pr = zr - conf_zr, pi = zi - conf_zi;
                    if (pr * pr + pi * pi < _interior_eps2 * zrad && conf_D2 < 1.0) {
                        if (++conf_count >= _interior_confirm) return trap ? trapc.value(0.0) : -2.f;
                        conf_zr = zr; conf_zi = zi; conf_D2 = 1; conf_next = iter + conf_P;
                    } else { conf_P = 0; ++conf_giveup; }
                }
            } else if ((iter & 15) == 0 && conf_giveup < 3) {
                double pr = zr - zsr, pi = zi - zsi;
                if (pr * pr + pi * pi < _interior_eps2 * zrad) {
                    conf_P = iter - save_iter; if (conf_P < 1) conf_P = 1;
                    conf_zr = zr; conf_zi = zi; conf_D2 = 1; conf_next = iter + conf_P; conf_count = 0;
                }
                if (iter - save_iter >= period_win) { zsr = zr; zsi = zi; save_iter = iter; period_win += period_win; }
            }
        }

        // periodic rescale to keep |w| ~ O(1)
        double wmag2 = wr * wr + wi * wi;
        if (wmag2 > 1e16 || (wmag2 < 1e-16 && wmag2 > 0.0)) {
            FloatExp wmag = fe_sqrt(fe_from(wmag2));
            S = fe_mul(S, wmag); s = fe_to_double(S); S2 = fe_mul(S, S);
            double inv = 1.0 / fe_to_double(wmag);
            wr *= inv; wi *= inv;
            dr = fe_to_double(fe_div(dcr, S)); di = fe_to_double(fe_div(dci, S));
        }

        // Zhuoran rebase to the critical point (m = 0). Rebase whenever the pixel
        // orbit is smaller than the delta (|z|^2 < |dz|^2), matching the double
        // path -- NOT only when |z| is tiny in double. Far-from-reference pixels
        // have a large delta (dz ~ O(1)); gating rebasing behind zrad<1e-8 misses
        // their |z|<|dz| rebases, so the delta grows and its 53-bit mantissa drifts
        // -> large-area glitch. dzmag2 = |S w|^2 in double (0 if S underflows deep,
        // where the zrad<1e-8 branch still handles the near-zero reference passes).
        double dzmag2 = s * s * (wr * wr + wi * wi);
        if (zrad < 1e-8 || zrad < dzmag2 || (!per && m >= reflen)) {
            FloatExp Xmr = m ? _zfr_fe[rm] : FloatExp{ 0.0, 0 };
            FloatExp Xmi = m ? _zfi_fe[rm] : FloatExp{ 0.0, 0 };
            FloatExp Swr = fe_mul_d(S, wr), Swi = fe_mul_d(S, wi);      // S w = dz_true
            FloatExp zrfe = fe_add(Xmr, Swr), zife = fe_add(Xmi, Swi);
            FloatExp zradfe = fe_add(fe_mul(zrfe, zrfe), fe_mul(zife, zife));
            FloatExp dzfe = fe_add(fe_mul(Swr, Swr), fe_mul(Swi, Swi));
            if ((!per && m >= reflen) || fe_abs_less(zradfe, dzfe)) {
                FloatExp Snew = fe_sqrt(zradfe);
                if (Snew.m == 0.0) { wr = wi = 0.0; }
                else {
                    S = Snew; s = fe_to_double(S); S2 = fe_mul(S, S);
                    wr = fe_to_double(fe_div(zrfe, S)); wi = fe_to_double(fe_div(zife, S));
                    dr = fe_to_double(fe_div(dcr, S)); di = fe_to_double(fe_div(dci, S));
                }
                m = 0; rm = -1;      // critical point; next ++ makes rm = 0 (mirrors m-1)
            }
        }
    }
    return trap ? trapc.value(0.0) : -2.f;   // interior (hit maxit): trap-colour or sentinel
}

// AVX2 4-wide rescaled deep-zoom kernel. Each lane carries its own scale S, delta
// (wr,wi), reference index m and iteration count, so lanes stay independent and
// coherent adjacent pixels vectorise well. The heavy quadratic step is done in a
// __m256d; BLA skips, |w| rescales and Zhuoran rebases are rare per-lane scalar
// events. Mirrors pixelRescaled op-for-op (see that function for the math).
void Mandel::solveRescaledSimd4(const FloatExp* Dcr, const FloatExp* Dci, int g, int lanes,
                                int mx_ref_it, int mxit, int c_method, float* out) const {
    const double ESC2 = (double)_ESCAPE_RADIUS * _ESCAPE_RADIUS;
    const double LG2 = log(2.0);
    const int reflen = mx_ref_it + 1;

    alignas(32) double wr[4] = { 0,0,0,0 }, wi[4] = { 0,0,0,0 };
    alignas(32) double sS[4] = { 1,1,1,1 }, dr[4] = { 0,0,0,0 }, di[4] = { 0,0,0,0 };
    FloatExp S[4], dcr[4], dci[4];
    int m[4] = { 0,0,0,0 }, iter[4] = { 0,0,0,0 };
    bool act[4] = { false,false,false,false };

    for (int l = 0; l < lanes; ++l) {
        FloatExp DCr = Dcr[g + l], DCi = Dci[g + l];
        dcr[l] = DCr; dci[l] = DCi;
        FloatExp Sl = fe_sqrt(fe_add(fe_mul(DCr, DCr), fe_mul(DCi, DCi)));
        if (Sl.m == 0.0) { S[l] = FloatExp{ 1.0, 0 }; wr[l] = wi[l] = dr[l] = di[l] = 0.0; }
        else {
            S[l] = Sl;
            wr[l] = fe_to_double(fe_div(DCr, Sl)); wi[l] = fe_to_double(fe_div(DCi, Sl));
            dr[l] = wr[l]; di[l] = wi[l];
        }
        sS[l] = fe_to_double(S[l]);
        m[l] = 1; iter[l] = 1;
        out[g + l] = -2.f;      // interior default
        act[l] = true;
    }

    const __m256d two = _mm256_set1_pd(2.0);
    bool anyactive = (lanes > 0);
    while (anyactive) {
        if (_flag_halt) return;

        // Phase 1: per-lane scalar BLA. A skipping lane advances itself and sits
        // out this round's vector step (stepMask bit clear).
        int stepMask = 0;
        for (int l = 0; l < lanes; ++l) {
            if (!act[l]) continue;
            if (iter[l] >= mxit) { out[g + l] = -2.f; act[l] = false; continue; }
            bool skipped = false;
            if (_use_bla && m[l] >= 2) {
                int skip = tryBLAfe(m[l] - 1, S[l], fe_mul(S[l], S[l]), wr[l], wi[l], dr[l], di[l], ESC2, reflen);
                if (skip > 0) {
                    m[l] += skip; iter[l] += skip;
                    sS[l] = fe_to_double(S[l]);
                    dr[l] = fe_to_double(fe_div(dcr[l], S[l]));
                    di[l] = fe_to_double(fe_div(dci[l], S[l]));
                    skipped = true;
                }
            }
            if (!skipped) stepMask |= 1 << l;
        }

        // Phase 2: gather reference values and run the quadratic step for the
        // lanes that are stepping this round (garbage lanes compute unused zeros).
        alignas(32) double Xpr[4] = { 0,0,0,0 }, Xpi[4] = { 0,0,0,0 };
        alignas(32) double Xqr[4] = { 0,0,0,0 }, Xqi[4] = { 0,0,0,0 };
        for (int l = 0; l < 4; ++l) {
            if (!(stepMask & (1 << l))) continue;
            int mm = m[l];
            if (mm) { Xpr[l] = _zfr[mm - 1]; Xpi[l] = _zfi[mm - 1]; }   // X_m
            int pi = mm > mx_ref_it ? mx_ref_it : mm;                  // X_{m+1} = _zfr[mm]
            Xqr[l] = _zfr[pi]; Xqi[l] = _zfi[pi];
        }
        __m256d vwr = _mm256_load_pd(wr), vwi = _mm256_load_pd(wi);
        __m256d vs = _mm256_load_pd(sS), vdr = _mm256_load_pd(dr), vdi = _mm256_load_pd(di);
        __m256d vXpr = _mm256_load_pd(Xpr), vXpi = _mm256_load_pd(Xpi);
        __m256d w2r = _mm256_sub_pd(_mm256_mul_pd(vwr, vwr), _mm256_mul_pd(vwi, vwi));
        __m256d w2i = _mm256_mul_pd(two, _mm256_mul_pd(vwr, vwi));
        // nwr = 2*(Xpr*wr - Xpi*wi) + s*w2r + dr
        __m256d a = _mm256_sub_pd(_mm256_mul_pd(vXpr, vwr), _mm256_mul_pd(vXpi, vwi));
        __m256d nwr = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(two, a), _mm256_mul_pd(vs, w2r)), vdr);
        // nwi = 2*(Xpr*wi + Xpi*wr) + s*w2i + di
        __m256d b = _mm256_add_pd(_mm256_mul_pd(vXpr, vwi), _mm256_mul_pd(vXpi, vwr));
        __m256d nwi = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(two, b), _mm256_mul_pd(vs, w2i)), vdi);
        __m256d vXqr = _mm256_load_pd(Xqr), vXqi = _mm256_load_pd(Xqi);
        __m256d zr = _mm256_add_pd(vXqr, _mm256_mul_pd(vs, nwr));
        __m256d zi = _mm256_add_pd(vXqi, _mm256_mul_pd(vs, nwi));
        __m256d zrad = _mm256_add_pd(_mm256_mul_pd(zr, zr), _mm256_mul_pd(zi, zi));
        __m256d wmag2 = _mm256_add_pd(_mm256_mul_pd(nwr, nwr), _mm256_mul_pd(nwi, nwi));
        alignas(32) double nwr_[4], nwi_[4], zrad_[4], wmag2_[4];
        _mm256_store_pd(nwr_, nwr); _mm256_store_pd(nwi_, nwi);
        _mm256_store_pd(zrad_, zrad); _mm256_store_pd(wmag2_, wmag2);

        // Phase 3: per-lane scalar commit + rare rescale / rebase / escape events.
        for (int l = 0; l < lanes; ++l) {
            if (!(stepMask & (1 << l))) continue;
            wr[l] = nwr_[l]; wi[l] = nwi_[l]; ++m[l]; ++iter[l];
            double zrl = zrad_[l];
            if (zrl > ESC2) {
                out[g + l] = (float)((double)iter[l] - log(log(zrl) / 2.0 / LG2) / LG2);
                act[l] = false; continue;
            }
            double wm = wmag2_[l];
            if (wm > 1e16 || (wm < 1e-16 && wm > 0.0)) {          // keep |w| ~ O(1)
                FloatExp wmag = fe_sqrt(fe_from(wm));
                S[l] = fe_mul(S[l], wmag); sS[l] = fe_to_double(S[l]);
                double inv = 1.0 / fe_to_double(wmag);
                wr[l] *= inv; wi[l] *= inv;
                dr[l] = fe_to_double(fe_div(dcr[l], S[l])); di[l] = fe_to_double(fe_div(dci[l], S[l]));
            }
            double dzmag2 = sS[l] * sS[l] * (wr[l] * wr[l] + wi[l] * wi[l]);   // |S w|^2
            if (zrl < 1e-8 || zrl < dzmag2 || m[l] >= reflen) {   // Zhuoran rebase (|z|<|dz|)
                int mm = m[l];
                FloatExp Xmr = mm ? _zfr_fe[mm - 1] : FloatExp{ 0.0, 0 };
                FloatExp Xmi = mm ? _zfi_fe[mm - 1] : FloatExp{ 0.0, 0 };
                FloatExp Swr = fe_mul_d(S[l], wr[l]), Swi = fe_mul_d(S[l], wi[l]);
                FloatExp zrfe = fe_add(Xmr, Swr), zife = fe_add(Xmi, Swi);
                FloatExp zradfe = fe_add(fe_mul(zrfe, zrfe), fe_mul(zife, zife));
                FloatExp dzfe = fe_add(fe_mul(Swr, Swr), fe_mul(Swi, Swi));
                if (mm >= reflen || fe_abs_less(zradfe, dzfe)) {
                    FloatExp Snew = fe_sqrt(zradfe);
                    if (Snew.m == 0.0) { wr[l] = wi[l] = 0.0; }
                    else {
                        S[l] = Snew; sS[l] = fe_to_double(S[l]);
                        wr[l] = fe_to_double(fe_div(zrfe, S[l])); wi[l] = fe_to_double(fe_div(zife, S[l]));
                        dr[l] = fe_to_double(fe_div(dcr[l], S[l])); di[l] = fe_to_double(fe_div(dci[l], S[l]));
                    }
                    m[l] = 0;
                }
            }
        }

        anyactive = false;
        for (int l = 0; l < lanes; ++l) if (act[l]) { anyactive = true; break; }
    }
}

int Mandel::createRef(std::set<std::array<int, 4>>& s, int pr_it, int mxit, bool random,
                      int c_method, bool view_center) {
    if (s.empty()) return false;
    _ref_period = 0;   // a normal (non-periodic) reference; delta loop indexes linearly
    ++_ref_cnt;
    std::array<int, 4> p{ _h / 2, _w / 2, 0, 0 };
    if (view_center) {
        _ref_virtual = true;
        _fe_cutoff_sensitive = false;
        _ref_x = (_w - 1) * 0.5;
        _ref_y = (_h - 1) * 0.5;
        // Exact geometric center = c0 + spacing*(size-1)/2. This may be a
        // half-pixel position for even dimensions, but is independent of which
        // actual center pixel would otherwise be selected.
        mpf_set(_ref_z_re, _c0_re);
        mpf_mul_ui(_t1, _dx, _w - 1); mpf_div_ui(_t1, _t1, 2);
        mpf_add(_ref_z_re, _ref_z_re, _t1);
        mpf_set(_ref_z_im, _c0_im);
        mpf_mul_ui(_t1, _dy, _h - 1); mpf_div_ui(_t1, _t1, 2);
        mpf_add(_ref_z_im, _ref_z_im, _t1);
    } else {
        _ref_virtual = false;
        if (random) {
            auto r = rand() % s.size();
            auto it = std::begin(s);
            std::advance(it, r);
            p = *it;
            s.erase(it);
        } else {
            // Glitch-guided: reuse the pixel with the smallest orbit magnitude
            // chosen in stepParallel. Fall back to a random pending pixel.
            auto it = s.find(_new_ref);
            if (it == s.end()) {
                auto r = rand() % s.size();
                it = std::begin(s);
                std::advance(it, r);
            }
            p = *it;
            s.erase(it);
        }
        _ref_x = p[1] + (double)p[3] / _sub;
        _ref_y = p[0] + (double)p[2] / _sub;
        // _ref_z = _c0 + pixel/sub offsets.
        mpf_set(_ref_z_re, _c0_re);
        mpf_set(_ref_z_im, _c0_im);
        mpf_mul_ui(_t1, _dx, p[1]);
        mpf_add(_ref_z_re, _ref_z_re, _t1);
        mpf_mul_ui(_t1, _dx, abs(p[3]));
        if (p[3] < 0) mpf_neg(_t1, _t1);
        mpf_div_ui(_t1, _t1, _sub);
        mpf_add(_ref_z_re, _ref_z_re, _t1);
        mpf_mul_ui(_t1, _dy, p[0]);
        mpf_add(_ref_z_im, _ref_z_im, _t1);
        mpf_mul_ui(_t1, _dy, abs(p[2]));
        if (p[2] < 0) mpf_neg(_t1, _t1);
        mpf_div_ui(_t1, _t1, _sub);
        mpf_add(_ref_z_im, _ref_z_im, _t1);
    }

    // _ref_z_f = static_cast<Comp>(_ref_z);
    _ref_z_f = Comp{ mpf_get_ld(_ref_z_re), mpf_get_ld(_ref_z_im) };
    
    _ref = p;
    
    // _z[0] = _ref_z;
    mpf_set(_z_re[0], _ref_z_re);
    mpf_set(_z_im[0], _ref_z_im);
    if (g_bigfixed == -2) { const char* e = getenv("MANDEL_BIGFIXED"); g_bigfixed = e ? atoi(e) : -1; }
    // Auto (env unset): use BigFixed on the deep floatexp path, where its high-half
    // short product beats mpf; env forces on (>=1) or off (0). mpf elsewhere.
    _use_bigfixed = (g_bigfixed == -1) ? _use_floatexp : (g_bigfixed > 0);
    if (_use_bigfixed) {
        // Fixed-point orbit: L limbs = precision/64 + guard for near-zero passes
        // (fixed-point loses relative precision as |Z|->0; the extra guard limbs
        // keep >=53 significant bits so the deep floatexp shadow stays accurate).
        int prec = (int)mpf_get_prec(_ref_z_re);
        _bfL = (prec + 63) / 64 + 2;                 // ceil(prec/64) fractional + integer + 1 round-guard limb
        _bftmp.assign((size_t)2 * _bfL, 0ull);
        bf_from_mpf(_bc_re, _ref_z_re, _bfL); bf_from_mpf(_bc_im, _ref_z_im, _bfL);
        _bz_re[0] = _bc_re; _bz_im[0] = _bc_im;       // z_0 = c
        _bt1.setL(_bfL); _bt2.setL(_bfL); _bab.setL(_bfL); _bre.setL(_bfL); _bim.setL(_bfL);
    }

    _zf[0] = _ref_z_f;
    _zfr[0] = _ref_z_f.real();
    _zfi[0] = _ref_z_f.imag();
    if (_use_floatexp) { _zfr_fe[0] = mpf_to_fe(_ref_z_re); _zfi_fe[0] = mpf_to_fe(_ref_z_im); }

    _df[0] = Comp{ 1 };
    _dfr[0] = 1.0; _dfi[0] = 0.0;
    _dfe_r = FloatExp{ 1.0, 0 }; _dfe_i = FloatExp{ 0.0, 0 };   // D_0 = 1

    // _SA_delta = static_cast<Comp>(_dx * _w / 2 + _dy * _h / 2);
    mpf_mul_ui(_t1, _dx, _w);
    mpf_div_ui(_t1, _t1, 2);
    mpf_mul_ui(_t2, _dy, _h);
    mpf_div_ui(_t2, _t2, 2);
    _SA_delta = { mpf_get_ld(_t1), mpf_get_ld(_t2) };
    for (int i = 0; i < _SA_N; ++i) _Adf_old[i] = _Bdf_old[i] = 0;
    _Adf_old[0] = _SA_delta;

    _SA_flag = true;
    _SA_it = 0;

    for (int i = 1; i <= mxit; ++i) {
        if (_flag_halt) break;
        if (!calCoefficient(i, pr_it, c_method)) {
            if (_ref_virtual && _use_floatexp && i >= mxit - 16)
                _fe_cutoff_sensitive = true;
            if (!_ref_virtual) {
                // The reference (centre) pixel is coloured here, not in the pixel
                // loop. getEscapeTime only gives the smooth iteration count, so for
                // Feather/EDE compute the reference's own value (one GMP pixel) or
                // it shows a stray dot in the wrong colour space at the view centre.
                if (c_method) {
                    float rv = accuratePointCompute(_ref_z_re, _ref_z_im, mxit, c_method);
                    if ((c_method & ColoringMethod::EXTERIOR_DIST_EST) && rv >= 0) rv /= (float)mpf_get_d(_dx);
                    setPixel(_ref, rv);
                } else {
                    if (_use_bigfixed) { bf_to_mpf(_z_re[i & 1], _bz_re[i & 1]); bf_to_mpf(_z_im[i & 1], _bz_im[i & 1]); }
                    setPixel(_ref, getEscapeTime(_z_re[i & 1], _z_im[i & 1], i));
                }
            }
            return i;
        }
    }
    if (_ref_virtual && _use_floatexp) {
        // A reference that only escapes just beyond maxit makes classification
        // exquisitely sensitive to accumulated double perturbation rounding.
        _fe_cutoff_sensitive =
            accuratePointCompute(_ref_z_re, _ref_z_im, mxit + 64, c_method) >= 0;
    }
    return mxit;
}

bool Mandel::calCoefficient(int i, int pr_it, int c_method) {
    // Track the reference derivative for EDE (|D|), the normal map (arg D) and the
    // DE overlay (|D|).
    const bool deriv = (c_method & ColoringMethod::EXTERIOR_DIST_EST)
                    || (c_method & ColoringMethod::NORMAL_MAP)
                    || (c_method & ColoringMethod::DE_OVERLAY);
    // The full-precision orbit uses two rotating buffers: c = current (i), p =
    // previous (i-1). Opposite parity, so the two are always distinct.
    const int c = i & 1, p = (i - 1) & 1;
    // _z[i] = _z[i - 1] * _z[i - 1] + _ref_z;  (complex square via Karatsuba:
    // re = (a+b)(a-b), im = 2ab -> 2 big multiplies instead of 3)
    // _z[i] = _z[i - 1] * _z[i - 1] + _ref_z;  (complex square via Karatsuba:
    // re = (a+b)(a-b), im = 2ab -> 2 big multiplies instead of 3)
    if (_use_bigfixed) {
        // Fixed-point (BigFixed) orbit step: same Karatsuba complex square, but the
        // fixed-size fixed-point bignum avoids mpf's dynamic-size/rounding overhead.
        const BigFixed& a = _bz_re[p]; const BigFixed& b = _bz_im[p];
        bf_add(_bt1, a, b);                        // a + b
        bf_sub(_bt2, a, b);                        // a - b
        bf_mul(_bab, a, b, _bftmp.data());         // a * b
        bf_mul(_bre, _bt1, _bt2, _bftmp.data());   // a^2 - b^2
        bf_add(_bim, _bab, _bab);                  // 2ab
        bf_add(_bz_re[c], _bre, _bc_re);           // + c
        bf_add(_bz_im[c], _bim, _bc_im);
        if (_use_floatexp) {
            // Deep path needs both shadows. bf_to_fe already extracts the top-2-limb
            // mantissa+exponent; derive the double from it via fe_to_double (a cheap
            // ldexp) instead of a second limb scan in toDouble(). fe_to_double(bf_to_fe(x))
            // == x.toDouble() exactly (same top bits, same rounding), so this is exact.
            FloatExp fr = bf_to_fe(_bz_re[c]), fi = bf_to_fe(_bz_im[c]);
            _zfr_fe[i] = fr; _zfi_fe[i] = fi;
            double zr = fe_to_double(fr), zi = fe_to_double(fi);
            _zf[i] = Comp{ zr, zi };
            _zfr[i] = zr; _zfi[i] = zi;
        } else {
            double zr = _bz_re[c].toDouble(), zi = _bz_im[c].toDouble();
            _zf[i] = Comp{ zr, zi };
            _zfr[i] = zr; _zfi[i] = zi;
        }
    } else {
    mpf_add(_t1, _z_re[p], _z_im[p]);       // a + b
    mpf_sub(_t2, _z_re[p], _z_im[p]);       // a - b
    mpf_mul(_z_im[c], _z_re[p], _z_im[p]);  // a*b   (read a,b before _z_re[c] write)
    mpf_mul(_z_re[c], _t1, _t2);            // a^2 - b^2
    mpf_mul_ui(_z_im[c], _z_im[c], 2);      // 2ab
    mpf_add(_z_re[c], _z_re[c], _ref_z_re);
    mpf_add(_z_im[c], _z_im[c], _ref_z_im);

    // _zf[i] = static_cast<Comp>(_z[i]);
    _zf[i] = Comp{ mpf_get_ld(_z_re[c]), mpf_get_ld(_z_im[c]) };
    _zfr[i] = _zf[i].real();
    _zfi[i] = _zf[i].imag();
    if (_use_floatexp) { _zfr_fe[i] = mpf_to_fe(_z_re[c]); _zfi_fe[i] = mpf_to_fe(_z_im[c]); }
    }
    // _z_m3[i] = tmp_z.abs().get_real_imag().first / 1000; // for Pauldelbrot condition
    _z_m3[i] = { (_zfr[i] * _zfr[i] + _zfi[i] * _zfi[i]) / 1000000 };

    if (deriv) {
        // D[i] = 2 * D[i-1] * Z[i-1] + 1  (only feeds the double shadow _dfr/_dfi).
        // Iterate directly, no mpf: floatexp on the deep path (Z[i-1] underflows the
        // double shadow near a minibrot's zero passes), plain double otherwise.
        if (_use_floatexp) {
            FloatExp a = _dfe_r, b = _dfe_i;
            FloatExp zr = _zfr_fe[i - 1], zi = _zfi_fe[i - 1];
            FloatExp re = fe_add(fe_scale2(fe_sub(fe_mul(a, zr), fe_mul(b, zi)), 1), FloatExp{ 1.0, 0 });
            FloatExp im = fe_scale2(fe_add(fe_mul(a, zi), fe_mul(b, zr)), 1);
            _dfe_r = re; _dfe_i = im;
            _dfr[i] = fe_to_double(re); _dfi[i] = fe_to_double(im);
        } else {
            double a = _dfr[i - 1], b = _dfi[i - 1], zr = _zfr[i - 1], zi = _zfi[i - 1];
            _dfr[i] = 2.0 * (a * zr - b * zi) + 1.0;
            _dfi[i] = 2.0 * (a * zi + b * zr);
        }
        _df[i] = Comp{ _dfr[i], _dfi[i] };
    }

    if (escape(_zf[i])) return false;
    if (i <= pr_it) {
        if (_SA_flag) {
            for (int j = 0; j < _SA_N; ++j) {
                _Adf_new[j] = Comp{ 2 } * _zf[i - 1] * _Adf_old[j];
                for (int k = 0; k < j / 2; ++k) {
                    _Adf_new[j] += Comp{ 2 } * _Adf_old[k] * _Adf_old[j - k - 1];
                }
                if (j % 2) _Adf_new[j] += _Adf_old[j / 2] * _Adf_old[j / 2];
            }
            
            if (deriv) {
                for (int j = 0; j < _SA_N; ++j) {
                    _Bdf_new[j] = _Bdf_old[j] * _zf[i - 1] + _Adf_old[j] * _df[i - 1];
                    for (int k = 0; k < j; ++k) {
                        _Bdf_new[j] += _Adf_old[k] * _Bdf_old[j - k - 1];
                    }
                    _Bdf_new[j] *= Comp{ 2 };
                }
            }
            _Adf_new[0] += _SA_delta;
            int order = SACheckMagnitude();
            
            if (order < 0) {
                _SA_flag = false;
            }
            else {
                _SA_it = i;
                _SA_order = order;
                for (int j = 0; j < _SA_N; ++j) _Adf_old[j] = _Adf_new[j];
                if (c_method & ColoringMethod::EXTERIOR_DIST_EST) {
                    for (int j = 0; j < _SA_N; ++j) _Bdf_old[j] = _Bdf_new[j];
                }
            }
        }
        _mx_coef = i;
    }
    return true;
}

int Mandel::SACheckMagnitude() const {
    int pre_mn[_SA_N], suf_mx[_SA_N];
    int re, im;
    re = get_exp(_Adf_new[0].real());
    im = get_exp(_Adf_new[0].imag());
    pre_mn[0] = std::min(re, im);
    
    re = get_exp(_Adf_new[_SA_N - 1].real());
    im = get_exp(_Adf_new[_SA_N - 1].imag());
    suf_mx[_SA_N - 1] = std::max(re, im);
    for (int i = 1; i < _SA_N; ++i) {
        re = get_exp(_Adf_new[i].real());
        im = get_exp(_Adf_new[i].imag());
        pre_mn[i] = std::max(pre_mn[i - 1], std::min(re, im));
    }
    for (int i = _SA_N - 2; i >= 0; --i) {
        re = get_exp(_Adf_new[i].real());
        im = get_exp(_Adf_new[i].imag());
        suf_mx[i] = std::max(suf_mx[i + 1], std::max(re, im));
    }
    
    for (int i = 0; i < _SA_N - 10; ++i) {
        if (pre_mn[i] - suf_mx[i + 1] >= 80) return std::max(i, 5);
    }
    return -1;
}

inline int Mandel::getIndex(std::array<int, 4>& arr) const {
    return (arr[0] * _sub + _sub / 2 + arr[2]) * _w * _sub + (arr[1] * _sub + _sub / 2 + arr[3]);
}

inline int Mandel::getIndex(int i, int j, int u, int v) const {
    return (i * _sub + _sub / 2 + u) * _w * _sub + (j * _sub + _sub / 2 + v);
}

void Mandel::SetHalt(bool flag) {
    _flag_halt = flag;
}

void Mandel::SetProgress(std::atomic<float>* progress, float offset, float scale) {
    _progress = progress;
    _progress_offset = offset;
    _progress_scale = scale;
}

void Mandel::progressSet(double local) {
    if (_progress)
        _progress->store((float)(_progress_offset + _progress_scale
            * std::clamp(local, 0.0, 1.0)), std::memory_order_relaxed);
}

void Mandel::progressBegin(int total, double begin, double span) {
    _progress_done.store(0, std::memory_order_relaxed);
    _progress_total = total;
    _progress_report_step = std::max(1, total / 200);
    _progress_begin = begin;
    _progress_span = span;
    progressSet(begin);
}

void Mandel::progressAdvance() {
    if (!_progress || _progress_total <= 0) return;
    int done = _progress_done.fetch_add(1, std::memory_order_relaxed) + 1;
    if (done == _progress_total || done % _progress_report_step == 0)
        progressSet(_progress_begin + _progress_span * done / _progress_total);
}

inline void Mandel::markDone(int i) {
    if (_done[i]) return;
    _done[i] = true;
    progressAdvance();
}
