// Step 1-2 of period detection: find the dominant minibrot period in a view via
// ball (interval) arithmetic, then Newton-refine to the nucleus. Standalone /
// verifiable before engine integration.
//
// Ball period: iterate the disk of c-values (center c0, radius r) under z->z^2+c.
// The image is ~a disk with centre z_n and radius R_n:
//     R_{n+1} = (2|z_n| + R_n) R_n + r,   z_{n+1} = z_n^2 + c0.
// The first p with |z_p| < R_p means the image disk contains 0 => a period-p
// nucleus lies in the original disk.
//
// Newton: refine c so that z_p(c) = 0 (the periodic critical orbit returns to 0),
// using c <- c - z_p / (dz_p/dc).
//
// Usage: period cx cy scaleExp [maxp]

#include <gmp.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>

#include "floatexp.h"

struct Cmpf { mpf_t re, im; };
static void cinit(Cmpf& z) { mpf_init(z.re); mpf_init(z.im); }
static void cclear(Cmpf& z) { mpf_clear(z.re); mpf_clear(z.im); }

// scratch
static mpf_t t1, t2, t3;

// z = z^2 + c
static void csqadd(Cmpf& z, const Cmpf& c) {
    mpf_mul(t1, z.re, z.re);     // re^2
    mpf_mul(t2, z.im, z.im);     // im^2
    mpf_mul(t3, z.re, z.im);     // re*im
    mpf_sub(z.re, t1, t2);       // re^2 - im^2
    mpf_add(z.re, z.re, c.re);
    mpf_mul_ui(z.im, t3, 2);     // 2 re im
    mpf_add(z.im, z.im, c.im);
}

// (a) = (a)*(b)  complex, using scratch; result into out (out may not alias a,b)
static void cmul(Cmpf& out, const Cmpf& a, const Cmpf& b) {
    mpf_mul(t1, a.re, b.re);
    mpf_mul(t2, a.im, b.im);
    mpf_mul(t3, a.re, b.im);
    mpf_sub(out.re, t1, t2);
    mpf_mul(t1, a.im, b.re);
    mpf_add(out.im, t3, t1);
}

static double cabs_d(const Cmpf& z) {
    return hypot(mpf_get_d(z.re), mpf_get_d(z.im));
}

static FloatExp mpf_to_fe(mpf_srcptr x) {
    if (mpf_sgn(x) == 0) return FloatExp{ 0.0, 0 };
    long e;
    return FloatExp{ mpf_get_d_2exp(&e, x), (int64_t)e };
}

static FloatExp cabs_fe(const Cmpf& z) {
    FloatExp re = mpf_to_fe(z.re), im = mpf_to_fe(z.im);
    return fe_sqrt(fe_add(fe_mul(re, re), fe_mul(im, im)));
}

// Ball-arithmetic period finder. Returns period (0 if none within maxp).
static int findPeriod(const Cmpf& c0, FloatExp r, int maxp) {
    Cmpf z; cinit(z);
    mpf_set_ui(z.re, 0); mpf_set_ui(z.im, 0);
    FloatExp R{ 0.0, 0 };
    int period = 0;
    for (int p = 1; p <= maxp; ++p) {
        FloatExp za = cabs_fe(z);
        R = fe_add(fe_mul(fe_add(fe_mul_d(za, 2.0), R), R), r);
        csqadd(z, c0);
        FloatExp zb = cabs_fe(z);
        if (fe_abs_less(zb, R)) { period = p; break; }
        if (cabs_d(z) > 4.0) break;    // escaped: no bounded atom here
    }
    cclear(z);
    return period;
}

// Newton refine c0 -> nucleus of the given period. Returns iterations used.
static int newtonNucleus(Cmpf& c, int period, int maxit) {
    Cmpf z, dz, num; cinit(z); cinit(dz); cinit(num);
    int it = 0;
    for (; it < maxit; ++it) {
        mpf_set_ui(z.re, 0); mpf_set_ui(z.im, 0);
        mpf_set_ui(dz.re, 0); mpf_set_ui(dz.im, 0);
        for (int i = 0; i < period; ++i) {
            // dz = 2 z dz + 1
            cmul(num, z, dz);
            mpf_mul_ui(num.re, num.re, 2); mpf_mul_ui(num.im, num.im, 2);
            mpf_add_ui(num.re, num.re, 1);
            mpf_set(dz.re, num.re); mpf_set(dz.im, num.im);
            // z = z^2 + c
            csqadd(z, c);
        }
        // delta = z / dz  (complex divide); c -= delta
        mpf_mul(t1, dz.re, dz.re); mpf_mul(t2, dz.im, dz.im); mpf_add(t1, t1, t2); // |dz|^2
        // num = z * conj(dz)
        mpf_mul(num.re, z.re, dz.re); mpf_mul(t3, z.im, dz.im); mpf_add(num.re, num.re, t3);
        mpf_mul(num.im, z.im, dz.re); mpf_mul(t3, z.re, dz.im); mpf_sub(num.im, num.im, t3);
        mpf_div(num.re, num.re, t1); mpf_div(num.im, num.im, t1);
        mpf_sub(c.re, c.re, num.re); mpf_sub(c.im, c.im, num.im);
        FloatExp step = cabs_fe(num);
        if (step.m == 0.0 || step.e < -(int64_t)mpf_get_prec(c.re) + 32) { ++it; break; }
    }
    // residual |z_period(c)|
    mpf_set_ui(z.re, 0); mpf_set_ui(z.im, 0);
    for (int i = 0; i < period; ++i) csqadd(z, c);
    FloatExp resid = cabs_fe(z);
    printf("  Newton: %d iters, residual |z_p| = %.17g * 2^%lld\n",
           it, resid.m, (long long)resid.e);
    printf("  nucleus re = "); mpf_out_str(stdout, 10, 0, c.re); putchar('\n');
    printf("  nucleus im = "); mpf_out_str(stdout, 10, 0, c.im); putchar('\n');
    cclear(z); cclear(dz); cclear(num);
    return it;
}

int main(int argc, char** argv) {
    if (argc < 4) { fprintf(stderr, "usage: period cx cy scaleExp [maxp]\n"); return 1; }
    const char* cx = argv[1]; const char* cy = argv[2];
    int scaleExp = atoi(argv[3]);
    int maxp = argc > 4 ? atoi(argv[4]) : 4000000;

    std::string scale = "1"; for (int i = 0; i < scaleExp; ++i) scale += "0";
    // Near a Misiurewicz point, a nucleus at parameter distance d can have a
    // component width on the order of d^2. Keep about twice the view depth so
    // Newton can resolve the nucleus well enough to inspect that component.
    int precision = 2 * (int)(scale.size() * log(10) / log(2)) + 128;
    mpf_set_default_prec(precision);
    mpf_init(t1); mpf_init(t2); mpf_init(t3);

    Cmpf c0; cinit(c0);
    mpf_set_str(c0.re, cx, 10); mpf_set_str(c0.im, cy, 10);
    mpf_t radius, scale_mpf;
    mpf_init_set_ui(radius, 2);
    mpf_init_set_str(scale_mpf, scale.c_str(), 10);
    mpf_div(radius, radius, scale_mpf);
    FloatExp r = mpf_to_fe(radius);           // view radius (half-width), no underflow

    printf("scale=1e%d  view-radius=%.17g * 2^%lld  prec=%d bits  maxp=%d\n",
           scaleExp, r.m, (long long)r.e, precision, maxp);
    int p = findPeriod(c0, r, maxp);
    printf("dominant period = %d\n", p);
    if (p > 0) newtonNucleus(c0, p, 200);
    mpf_clear(radius);
    mpf_clear(scale_mpf);
    return 0;
}
