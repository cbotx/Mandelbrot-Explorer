// Clean-room standalone perturbation + BLA on a small tile (fr-p0).
//
// Purpose: isolate BLA correctness and speed from the full engine (no floatexp,
// no series approximation, no glitch-reference hunting, no SIMD). Everything is
// plain double perturbation with Zhuoran rebasing, plus an optional BLA table
// built exactly like the engine's buildBLA/tryBLA. It renders a WxH tile three
// ways and reports how they agree, the BLA skip fraction, and timings:
//
//   (A) GMP per-pixel brute force  -> ground truth iteration count (sparse grid)
//   (B) double perturbation, BLA off
//   (C) double perturbation, BLA on
//
// Metrics:
//   * iter-diff (B vs A) on the sparse oracle grid  -> perturbation correct?
//   * iter-diff (C vs B) on the full tile           -> BLA exact?
//   * render PSNR (C vs B)                           -> BLA visual quality
//   * BLA skip fraction, avg skip, time(B) vs time(C)
//
// The location must be shallow enough that dc and dc^2 fit in a normal double
// (roughly zoom <= 1e150); deeper zooms need floatexp and are out of scope here.
//
// Usage: pbench <location> <scaleExp> [W H maxit oracleStep]
//   location = deep1 | ticktock | flake | seahorse | <cx> <cy> is not supported;
//              pick a named location (coords embedded from test_cases.h).
// Env: PB_EPS (BLA eps, default 2^-53), PB_MINSKIP (min skip, default 8).

#include <gmp.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <chrono>

#include "test_cases.h"

using Clock = std::chrono::high_resolution_clock;
static double secs(Clock::time_point a, Clock::time_point b) {
    return std::chrono::duration_cast<std::chrono::duration<double>>(b - a).count();
}

static const double ESCAPE_R = 100.0;
static const double ESC2 = ESCAPE_R * ESCAPE_R;
static const double LOG2 = log(2.0);

// ---- reference orbit (double), computed in GMP at the center -------------
struct Ref {
    std::vector<double> zr, zi;   // X_n, n = 0..len  (X_0 = 0)
    int len = 0;                  // last valid index; orbit escaped or hit maxit
    bool bounded = false;         // did not escape within maxit
};

static Ref computeRef(mpf_t cr, mpf_t ci, int maxit) {
    Ref R;
    R.zr.assign(maxit + 2, 0.0);
    R.zi.assign(maxit + 2, 0.0);
    mpf_t xr, xi, xr2, xi2, t;
    mpf_inits(xr, xi, xr2, xi2, t, (mpf_ptr)0);
    mpf_set_ui(xr, 0); mpf_set_ui(xi, 0);
    R.zr[0] = 0.0; R.zi[0] = 0.0;
    int n = 0;
    for (; n < maxit; ++n) {
        // x = x^2 + c
        mpf_mul(xr2, xr, xr);
        mpf_mul(xi2, xi, xi);
        mpf_mul(t, xr, xi);       // t = xr*xi
        mpf_sub(xr, xr2, xi2);
        mpf_add(xr, xr, cr);      // new xr = xr^2 - xi^2 + cr
        mpf_mul_2exp(xi, t, 1);   // 2*xr*xi
        mpf_add(xi, xi, ci);      // new xi
        R.zr[n + 1] = mpf_get_d(xr);
        R.zi[n + 1] = mpf_get_d(xi);
        double zr = R.zr[n + 1], zi = R.zi[n + 1];
        if (zr * zr + zi * zi > ESC2) { R.len = n + 1; mpf_clears(xr, xi, xr2, xi2, t, (mpf_ptr)0); return R; }
    }
    R.len = n;            // == maxit
    R.bounded = true;
    mpf_clears(xr, xi, xr2, xi2, t, (mpf_ptr)0);
    return R;
}

// ---- BLA table -----------------------------------------------------------
struct BLAEntry { double ar, ai, br, bi, r2; int l; };

struct BLATable {
    std::vector<std::vector<BLAEntry>> lvl;   // lvl[p] : skip 2^p, start index 1 + i*2^p
    double rmax2 = 0.0;
    int minlevel = 0;
};

static BLATable buildBLA(const Ref& ref, double eps, double dcmax, int minskip) {
    BLATable T;
    T.minlevel = 0; while ((1 << (T.minlevel + 1)) <= minskip) ++T.minlevel;
    int reflen = ref.len;
    if (reflen < 3) return T;
    std::vector<BLAEntry> l0;
    l0.reserve(reflen - 1);
    for (int s = 1; s < reflen; ++s) {
        double ar = 2.0 * ref.zr[s], ai = 2.0 * ref.zi[s];   // A = 2 X_s
        double Amag = sqrt(ar * ar + ai * ai);
        double Rv = Amag > 0 ? eps * Amag - dcmax / Amag : 0.0;   // |B| = 1
        double r2 = Rv > 0 ? Rv * Rv : 0.0;
        if (r2 > T.rmax2) T.rmax2 = r2;
        l0.push_back({ ar, ai, 1.0, 0.0, r2, 1 });
    }
    T.lvl.push_back(std::move(l0));
    while (T.lvl.back().size() > 1) {
        const std::vector<BLAEntry>& prev = T.lvl.back();
        std::vector<BLAEntry> nxt;
        nxt.reserve(prev.size() / 2);
        for (size_t i = 0; i + 1 < prev.size(); i += 2) {
            const BLAEntry& x = prev[i];
            const BLAEntry& y = prev[i + 1];
            BLAEntry z;
            z.ar = y.ar * x.ar - y.ai * x.ai;
            z.ai = y.ar * x.ai + y.ai * x.ar;
            z.br = y.ar * x.br - y.ai * x.bi + y.br;
            z.bi = y.ar * x.bi + y.ai * x.br + y.bi;
            double Rx = sqrt(x.r2), Ry = sqrt(y.r2);
            double Axmag = sqrt(x.ar * x.ar + x.ai * x.ai);
            double Bxmag = sqrt(x.br * x.br + x.bi * x.bi);
            double cand = Axmag > 0 ? (Ry - Bxmag * dcmax) / Axmag : 0.0;
            double Rv = std::min(Rx, cand);
            z.r2 = Rv > 0 ? Rv * Rv : 0.0;
            z.l = x.l + y.l;
            nxt.push_back(z);
        }
        T.lvl.push_back(std::move(nxt));
    }
    return T;
}

static inline int trailing_zeros(unsigned v) {
    if (v == 0) return 31;
    int c = 0; while ((v & 1u) == 0) { v >>= 1; ++c; } return c;
}

// Largest valid BLA starting at reference index s (dz currently at index s).
// On success updates dz and returns skip; else returns 0.
static int tryBLA(const BLATable& T, const Ref& ref, int s,
                  double& dzr, double& dzi, double dcr, double dci, int reflen) {
    if (s < 1 || T.lvl.empty()) return 0;
    double zmag2 = dzr * dzr + dzi * dzi;
    int maxp = (int)T.lvl.size() - 1;
    int startp = (s == 1) ? maxp : trailing_zeros((unsigned)(s - 1));
    if (startp > maxp) startp = maxp;
    if (startp < T.minlevel) return 0;
    for (int p = startp; p >= T.minlevel; --p) {
        int i = (s - 1) >> p;
        if (i >= (int)T.lvl[p].size()) continue;
        const BLAEntry& b = T.lvl[p][i];
        if (b.r2 <= 0 || zmag2 >= b.r2) continue;
        int land = s + b.l;
        if (land >= reflen) continue;                 // keep landing inside stored ref
        double nzr = b.ar * dzr - b.ai * dzi + b.br * dcr - b.bi * dci;
        double nzi = b.ar * dzi + b.ai * dzr + b.br * dci + b.bi * dcr;
        double fr = ref.zr[land] + nzr, fi = ref.zi[land] + nzi;
        if (fr * fr + fi * fi > ESC2) return 0;        // would overshoot escape
        dzr = nzr; dzi = nzi;
        return b.l;
    }
    return 0;
}

// ---- per-pixel perturbation (double) -------------------------------------
// Returns smooth escape iteration, or -1 for interior (did not escape). Fills
// skip/normal counters when a BLA table is supplied.
static double pixelPert(const Ref& ref, const BLATable* T,
                        double dcr, double dci, int maxit,
                        long long& skips, long long& applies, long long& normals) {
    int m = 0;                 // reference index; dz is at index m
    double dzr = 0.0, dzi = 0.0;
    int iter = 0;
    int reflen = ref.len;
    double rmax2 = T ? T->rmax2 : 0.0;
    while (iter < maxit) {
        bool did = false;
        if (T && dzr * dzr + dzi * dzi < rmax2) {
            int sk = tryBLA(*T, ref, m, dzr, dzi, dcr, dci, reflen);
            if (sk > 0) { m += sk; iter += sk; skips += sk; ++applies; did = true; }
            else ++normals;
        }
        if (!did) {
            // dz_{m+1} = 2 X_m dz + dz^2 + dc
            double zr = ref.zr[m], zi = ref.zi[m];
            double ndzr = 2.0 * (zr * dzr - zi * dzi) + (dzr * dzr - dzi * dzi) + dcr;
            double ndzi = 2.0 * (zr * dzi + zi * dzr) + 2.0 * dzr * dzi + dci;
            dzr = ndzr; dzi = ndzi;
            ++m; ++iter;
        }
        double zr = ref.zr[m] + dzr, zi = ref.zi[m] + dzi;
        double zrad = zr * zr + zi * zi;
        if (zrad > ESC2) {
            return (double)iter + 1.0 - log(log(zrad) / 2.0 / LOG2) / LOG2;
        }
        // Zhuoran rebase
        if (zrad < dzr * dzr + dzi * dzi || m >= reflen) {
            dzr = zr; dzi = zi; m = 0;
        }
    }
    return -1.0;   // interior
}

// ---- GMP brute-force ground truth ----------------------------------------
static double pixelBrute(mpf_t cr, mpf_t ci, int maxit) {
    mpf_t xr, xi, xr2, xi2, t;
    mpf_inits(xr, xi, xr2, xi2, t, (mpf_ptr)0);
    mpf_set_ui(xr, 0); mpf_set_ui(xi, 0);
    double res = -1.0;
    for (int n = 0; n < maxit; ++n) {
        mpf_mul(xr2, xr, xr);
        mpf_mul(xi2, xi, xi);
        mpf_mul(t, xr, xi);
        mpf_sub(xr, xr2, xi2); mpf_add(xr, xr, cr);
        mpf_mul_2exp(xi, t, 1); mpf_add(xi, xi, ci);
        double zr = mpf_get_d(xr), zi = mpf_get_d(xi);
        double zrad = zr * zr + zi * zi;
        if (zrad > ESC2) { res = (double)(n + 1) + 1.0 - log(log(zrad) / 2.0 / LOG2) / LOG2; break; }
    }
    mpf_clears(xr, xi, xr2, xi2, t, (mpf_ptr)0);
    return res;
}

// smooth iteration -> grayscale byte (cyclic), interior -> 0
static int shade(double v) {
    if (v < 0) return 0;
    double t = 0.5 + 0.5 * sin(v * 0.15);
    int b = (int)(t * 255.0 + 0.5);
    return b < 0 ? 0 : (b > 255 ? 255 : b);
}

struct Loc { const char* name; const char* x; const char* y; };

int main(int argc, char** argv) {
    if (argc < 3) {
        fprintf(stderr, "usage: %s <location> <scaleExp> [W H maxit oracleStep]\n"
                        "  location = deep1 | ticktock | flake | seahorse\n", argv[0]);
        return 2;
    }
    std::string loc = argv[1];
    int scaleExp = atoi(argv[2]);
    int W = argc > 3 ? atoi(argv[3]) : 128;
    int H = argc > 4 ? atoi(argv[4]) : 128;
    int maxit = argc > 5 ? atoi(argv[5]) : 50000;
    int ostep = argc > 6 ? atoi(argv[6]) : 16;
    double eps = (getenv("PB_EPS") ? atof(getenv("PB_EPS")) : ldexp(1.0, -53));
    int minskip = (getenv("PB_MINSKIP") ? atoi(getenv("PB_MINSKIP")) : 8);

    Loc locs[] = {
        { "deep1", testcases::deep1_x, testcases::deep1_y },
        { "ticktock", testcases::ticktock_x, testcases::ticktock_y },
        { "flake", testcases::flake_x, testcases::flake_y },
    };
    const Loc* L = nullptr;
    for (auto& e : locs) if (loc == e.name) { L = &e; break; }
    if (!L) { fprintf(stderr, "unknown location '%s'\n", loc.c_str()); return 2; }

    int precision = (int)(scaleExp * log(10) / log(2)) + 64;
    if (getenv("PB_PREC")) precision = atoi(getenv("PB_PREC"));
    mpf_set_default_prec(precision);

    mpf_t cr, ci, scale, dw, dx, half;
    mpf_inits(cr, ci, scale, dw, dx, half, (mpf_ptr)0);
    mpf_set_str(cr, L->x, 10);
    mpf_set_str(ci, L->y, 10);
    // scale = 10^scaleExp ; half-width dw = 2/scale ; dx = 2*dw/(W-1)
    { std::string s = "1"; s.append(scaleExp, '0'); mpf_set_str(scale, s.c_str(), 10); }
    mpf_set_ui(dw, 2); mpf_div(dw, dw, scale);
    mpf_mul_ui(dx, dw, 2); mpf_div_ui(dx, dx, W - 1);

    printf("=== pbench %s  zoom=1e%d  tile=%dx%d  maxit=%d  prec=%d bits  eps=%.3g minskip=%d\n",
           loc.c_str(), scaleExp, W, H, maxit, precision, eps, minskip);

    // dc of the tile corner (double) to check double is adequate.
    double dx_d = mpf_get_d(dx);
    double dc_corner = dx_d * (W / 2);
    printf("    dx=%.3e  corner |dc|~%.3e  |dc|^2~%.3e  (need > 1e-300)\n",
           dx_d, dc_corner, dc_corner * dc_corner);

    // reference orbit at center
    auto t0 = Clock::now();
    Ref ref = computeRef(cr, ci, maxit);
    double t_ref = secs(t0, Clock::now());
    printf("    reference orbit: len=%d  %s  (%.3f s)\n",
           ref.len, ref.bounded ? "BOUNDED" : "escaped", t_ref);

    double dcmax = dc_corner;              // |dc| bound for a dc-independent table
    BLATable T = buildBLA(ref, eps, dcmax, minskip);
    printf("    BLA levels=%zu  rmax=%.3e (r^2=%.3e)\n",
           T.lvl.size(), sqrt(T.rmax2), T.rmax2);

    std::vector<double> itB(W * H), itC(W * H);

    // dc helper
    auto dcAt = [&](int i, int j, double& dcr, double& dci) {
        dcr = dx_d * (j - (W - 1) / 2.0);
        dci = dx_d * (i - (H - 1) / 2.0);
    };

    long long dummy = 0;
    // (B) BLA off
    t0 = Clock::now();
    for (int i = 0; i < H; ++i)
        for (int j = 0; j < W; ++j) {
            double dcr, dci; dcAt(i, j, dcr, dci);
            itB[i * W + j] = pixelPert(ref, nullptr, dcr, dci, maxit, dummy, dummy, dummy);
        }
    double t_B = secs(t0, Clock::now());

    // (C) BLA on
    long long skips = 0, applies = 0, normals = 0;
    t0 = Clock::now();
    for (int i = 0; i < H; ++i)
        for (int j = 0; j < W; ++j) {
            double dcr, dci; dcAt(i, j, dcr, dci);
            itC[i * W + j] = pixelPert(ref, &T, dcr, dci, maxit, skips, applies, normals);
        }
    double t_C = secs(t0, Clock::now());

    // ---- BLA exactness (C vs B) over full tile, + render PSNR ----
    double maxd = 0, sumd = 0; long ext = 0, clsmis = 0; int wi = -1, wj = -1;
    double mse = 0;
    for (int i = 0; i < H; ++i)
        for (int j = 0; j < W; ++j) {
            double b = itB[i * W + j], c = itC[i * W + j];
            bool bi = b < 0, ci2 = c < 0;
            if (bi != ci2) ++clsmis;
            else if (!bi) { double d = fabs(b - c); sumd += d; ++ext; if (d > maxd) { maxd = d; wi = i; wj = j; } }
            int sb = shade(b), sc = shade(c); double e = sb - sc; mse += e * e;
        }
    mse /= (W * H);
    double psnr = mse > 0 ? 10.0 * log10(255.0 * 255.0 / mse) : 1e9;

    // ---- perturbation correctness (B vs GMP oracle) over sparse grid ----
    double omaxd = 0, osumd = 0; long osamp = 0, ocls = 0; int owi = -1, owj = -1;
    double owtruth = 0;
    t0 = Clock::now();
    for (int i = 0; i < H; i += ostep)
        for (int j = 0; j < W; j += ostep) {
            mpf_t pcr, pci, off;
            mpf_inits(pcr, pci, off, (mpf_ptr)0);
            // pcr = cr + dx*(j-(W-1)/2)  -- MUST match dcAt's (W-1)/2.0 sub-pixel
            // center exactly, so use the integer (2j-(W-1)) then halve.
            mpf_set_si(off, 2 * j - (W - 1)); mpf_mul(off, off, dx); mpf_div_ui(off, off, 2); mpf_add(pcr, cr, off);
            mpf_set_si(off, 2 * i - (H - 1)); mpf_mul(off, off, dx); mpf_div_ui(off, off, 2); mpf_add(pci, ci, off);
            double truth = pixelBrute(pcr, pci, maxit);
            double b = itB[i * W + j];
            ++osamp;
            bool ti = truth < 0, bi = b < 0;
            if (ti != bi) ++ocls;
            else if (!ti) { double d = fabs(truth - b); osumd += d; if (d > omaxd) { omaxd = d; owi = i; owj = j; owtruth = truth; } }
            if (getenv("PB_DUMP")) printf("      (%3d,%3d) truth=%10.4f  B=%10.4f  diff=%9.4g\n", i, j, truth, b, ti||bi?-1.0:fabs(truth-b));
            mpf_clears(pcr, pci, off, (mpf_ptr)0);
        }
    double t_oracle = secs(t0, Clock::now());

    long long totalsteps = skips + normals;
    printf("\n  -- timing --\n");
    printf("    perturbation BLA off : %.3f s\n", t_B);
    printf("    perturbation BLA on  : %.3f s   speedup x%.2f\n", t_C, t_B / t_C);
    printf("    GMP oracle (%ld px)   : %.3f s\n", osamp, t_oracle);
    printf("\n  -- BLA skip stats --\n");
    printf("    applies=%lld  skipped-iters=%lld  normal-steps=%lld  avg-skip=%.1f  skip-frac=%.1f%%\n",
           applies, skips, normals, applies ? (double)skips / applies : 0.0,
           totalsteps ? 100.0 * skips / totalsteps : 0.0);
    printf("\n  -- BLA exactness (C vs B, full tile) --\n");
    printf("    exterior-both=%ld  class-mismatch=%ld  max-iter-diff=%.4g  mean=%.4g\n",
           ext, clsmis, maxd, ext ? sumd / ext : 0.0);
    if (wi >= 0) printf("    worst @ (%d,%d): B=%.5f C=%.5f\n", wi, wj, itB[wi * W + wj], itC[wi * W + wj]);
    printf("    render PSNR (C vs B) = %.2f dB\n", psnr);
    printf("\n  -- perturbation correctness (B vs GMP oracle, %ldpx grid) --\n", osamp);
    printf("    class-mismatch=%ld  max-iter-diff=%.4g  mean=%.4g\n",
           ocls, omaxd, osamp - ocls > 0 ? osumd / (osamp - ocls) : 0.0);
    if (owi >= 0) printf("    worst @ (%d,%d): truth=%.5f B=%.5f\n", owi, owj, owtruth, itB[owi * W + owj]);

    // write the BLA-on render for eyeballing
    if (getenv("PB_BMP")) {
        std::string fn = std::string("pbench_") + loc + ".pgm";
        FILE* f = fopen(fn.c_str(), "wb");
        if (f) {
            fprintf(f, "P5\n%d %d\n255\n", W, H);
            for (int i = 0; i < H; ++i) for (int j = 0; j < W; ++j) { unsigned char px = (unsigned char)shade(itC[i * W + j]); fwrite(&px, 1, 1, f); }
            fclose(f);
            printf("\n  wrote %s\n", fn.c_str());
        }
    }

    mpf_clears(cr, ci, scale, dw, dx, half, (mpf_ptr)0);
    return 0;
}
