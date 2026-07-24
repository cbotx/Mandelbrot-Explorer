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

// Ball-arithmetic period finder. Returns period (0 if none within maxp).
static int findPeriod(const Cmpf& c0, double r, int maxp) {
    Cmpf z; cinit(z);
    mpf_set_ui(z.re, 0); mpf_set_ui(z.im, 0);
    double R = 0.0;
    int period = 0;
    for (int p = 1; p <= maxp; ++p) {
        double za = cabs_d(z);
        R = (2 * za + R) * R + r;
        csqadd(z, c0);
        double zb = cabs_d(z);
        if (zb < R) { period = p; break; }
        if (zb > 4.0) break;    // escaped: no bounded atom here
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
        double step = cabs_d(num);
        if (step < 1e-300 || (step == 0)) { ++it; break; }
        // convergence: step tiny relative to c
        if (step < ldexp(cabs_d(c) + 1e-300, -80)) { ++it; break; }
    }
    // residual |z_period(c)|
    mpf_set_ui(z.re, 0); mpf_set_ui(z.im, 0);
    for (int i = 0; i < period; ++i) csqadd(z, c);
    double resid = cabs_d(z);
    printf("  Newton: %d iters, residual |z_p| = %.3e\n", it, resid);
    cclear(z); cclear(dz); cclear(num);
    return it;
}

int main(int argc, char** argv) {
    if (argc < 4) { fprintf(stderr, "usage: period cx cy scaleExp [maxp]\n"); return 1; }
    const char* cx = argv[1]; const char* cy = argv[2];
    int scaleExp = atoi(argv[3]);
    int maxp = argc > 4 ? atoi(argv[4]) : 4000000;

    std::string scale = "1"; for (int i = 0; i < scaleExp; ++i) scale += "0";
    int precision = (int)(scale.size() * log(10) / log(2)) + 40;
    mpf_set_default_prec(precision);
    mpf_init(t1); mpf_init(t2); mpf_init(t3);

    Cmpf c0; cinit(c0);
    mpf_set_str(c0.re, cx, 10); mpf_set_str(c0.im, cy, 10);
    double r = 2.0 * pow(10.0, -scaleExp);   // view radius (half-width)

    printf("scale=1e%d  view-radius=%.3e  prec=%d bits  maxp=%d\n", scaleExp, r, precision, maxp);
    int p = findPeriod(c0, r, maxp);
    printf("dominant period = %d\n", p);
    if (p > 0) newtonNucleus(c0, p, 200);
    return 0;
}
