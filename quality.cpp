// Standalone anti-aliasing quality experiment for chaotic Mandelbrot regions.
//
// For a few high-gradient "chaotic" views it compares supersampling patterns at
// equal sample budget against a dense ground truth, reporting PSNR (higher =
// closer to ground truth = better AA). No external deps: a self-contained double
// escape-time + a smooth cyclic palette. This decides which pattern to adopt in
// the real renderer before touching it.

#include <cstdio>
#include <cmath>
#include <cstdint>
#include <vector>
#include <string>
#include <algorithm>

static const double ESCAPE2 = 100.0 * 100.0;
static const double LOG2 = log(2.0);

// Smooth escape time; -1 for interior (does not escape).
static double smoothIter(double cr, double ci, int mxit) {
    double zr = 0, zi = 0;
    for (int i = 0; i < mxit; ++i) {
        double zr2 = zr * zr, zi2 = zi * zi;
        if (zr2 + zi2 > ESCAPE2) {
            double m = zr2 + zi2;
            return (i - log(log(m) / 2 / LOG2) / LOG2);
        }
        zi = 2 * zr * zi + ci;
        zr = zr2 - zi2 + cr;
    }
    return -1;
}

// Deterministic smooth palette (cyclic). Interior -> black. Returns 0..255 RGB.
static void palette(double it, double* rgb) {
    if (it < 0) { rgb[0] = rgb[1] = rgb[2] = 0; return; }
    double t = pow(it, 0.6) * 0.15;   // compress + cycle
    rgb[0] = 128 + 127 * sin(t + 0.0);
    rgb[1] = 128 + 127 * sin(t + 2.094);
    rgb[2] = 128 + 127 * sin(t + 4.188);
}

// sRGB-ish gamma helpers for linear-light averaging.
static inline double toLinear(double c) { return pow(c / 255.0, 2.2); }
static inline double toGamma(double l)  { return 255.0 * pow(std::max(0.0, l), 1.0 / 2.2); }

// Van der Corput / Halton low-discrepancy sequence.
static double halton(int i, int base) {
    double f = 1, r = 0;
    while (i > 0) { f /= base; r += f * (i % base); i /= base; }
    return r;
}

// Hash-based per-sample jitter in [0,1).
static double jitter(int x, int y, int k) {
    uint32_t h = (uint32_t)(x * 73856093) ^ (uint32_t)(y * 19349663) ^ (uint32_t)(k * 83492791);
    h ^= h >> 13; h *= 0x5bd1e995u; h ^= h >> 15;
    return (h & 0xffffff) / (double)0x1000000;
}

enum Pattern { CENTER, REGULAR, JITTERED, ROTGRID, HALTON };
enum Avg { LINEAR, GAMMA2, NAIVE };   // linear-light, gamma-2 RMS (current), naive sRGB mean

// Adaptive/progressive sampler: draws jittered samples in small batches, tracks
// the running linear-light mean + variance (Welford) per channel, and stops once
// the standard error of the mean falls below `tol` (relative to full white) or
// `maxN` is reached. Returns the pixel colour and sets *nUsed to the sample count.
static void adaptivePixel(double px, double py, double dx, double dy, int mxit,
                          double tol, int maxN, int pxi, int pyi, double* out, int* nUsed) {
    double mean[3] = { 0,0,0 }, M2[3] = { 0,0,0 };
    int n = 0;
    const int BATCH = 4, MINN = 4;
    while (n < maxN) {
        for (int k = 0; k < BATCH; ++k) {
            double sx = jitter(pxi, pyi, 2 * n) - 0.5, sy = jitter(pxi, pyi, 2 * n + 1) - 0.5;
            double rgb[3]; palette(smoothIter(px + sx * dx, py + sy * dy, mxit), rgb);
            ++n;
            for (int c = 0; c < 3; ++c) {
                double v = toLinear(rgb[c]);
                double d = v - mean[c]; mean[c] += d / n; M2[c] += d * (v - mean[c]);
            }
        }
        if (n >= MINN) {                       // standard error of the mean, worst channel
            double se = 0;
            for (int c = 0; c < 3; ++c) { double var = M2[c] / (n - 1); se = std::max(se, sqrt(var / n)); }
            if (toGamma(mean[0] + se) - toGamma(mean[0]) < tol * 255 && se < tol) break;
        }
    }
    for (int c = 0; c < 3; ++c) out[c] = toGamma(mean[c]);
    *nUsed = n;
}


// Average one output pixel's color for a given sampling pattern with `n` samples
// per axis (n*n total, except CENTER/HALTON), using averaging mode `avg`.
static void samplePixel(double px, double py, double dx, double dy, int mxit,
                        Pattern pat, int n, int pxi, int pyi, double* out, Avg avg = LINEAR) {
    double s0 = 0, s1 = 0, s2 = 0; int cnt = 0;
    auto add = [&](double sx, double sy) {
        double rgb[3]; palette(smoothIter(px + sx * dx, py + sy * dy, mxit), rgb);
        if (avg == LINEAR)      { s0 += toLinear(rgb[0]); s1 += toLinear(rgb[1]); s2 += toLinear(rgb[2]); }
        else if (avg == GAMMA2) { s0 += rgb[0] * rgb[0];  s1 += rgb[1] * rgb[1];  s2 += rgb[2] * rgb[2]; }
        else                    { s0 += rgb[0];           s1 += rgb[1];           s2 += rgb[2]; }
        ++cnt;
    };
    if (pat == CENTER) { add(0, 0); }
    else if (pat == HALTON) {
        int tot = n * n;
        for (int k = 0; k < tot; ++k) add(halton(k + 1, 2) - 0.5, halton(k + 1, 3) - 0.5);
    } else {
        for (int a = 0; a < n; ++a)
            for (int b = 0; b < n; ++b) {
                double cx = (a + 0.5) / n - 0.5, cy = (b + 0.5) / n - 0.5;
                if (pat == REGULAR) add(cx, cy);
                else if (pat == JITTERED)
                    add((a + jitter(pxi, pyi, a * n + b)) / n - 0.5,
                        (b + jitter(pxi, pyi, a * n + b + 1000)) / n - 0.5);
                else if (pat == ROTGRID) {   // rotated grid (atan(1/2) ~ 26.57 deg)
                    double rx = cx * 0.8944 - cy * 0.4472, ry = cx * 0.4472 + cy * 0.8944;
                    add(rx, ry);
                }
            }
    }
    out[0] = toGamma(s0 / cnt); out[1] = toGamma(s1 / cnt); out[2] = toGamma(s2 / cnt);
    if (avg == GAMMA2) { out[0] = sqrt(s0 / cnt); out[1] = sqrt(s1 / cnt); out[2] = sqrt(s2 / cnt); }
    else if (avg == NAIVE) { out[0] = s0 / cnt; out[1] = s1 / cnt; out[2] = s2 / cnt; }
}

struct View { const char* name; double cx, cy, scale; int mxit; };

static void run(const View& v, int W, int H) {
    double dw = 2.0 / v.scale, dh = dw * H / W;
    double dx = 2 * dw / (W - 1), dy = 2 * dh / (H - 1);
    double x0 = v.cx - dw, y0 = v.cy - dh;

    // Ground truth: unbiased low-discrepancy estimate (Halton, 1024 samples/pixel).
    std::vector<double> gt(W * H * 3);
    for (int i = 0; i < H; ++i)
        for (int j = 0; j < W; ++j)
            samplePixel(x0 + j * dx, y0 + i * dy, dx, dy, v.mxit, HALTON, 32, j, i, &gt[(i * W + j) * 3], LINEAR);

    struct Cand { const char* name; Pattern pat; int n; Avg avg; };
    Cand cands[] = {
        { "regular 3x3 (9)",   REGULAR,  3, LINEAR },
        { "jittered 3x3 (9)",  JITTERED, 3, LINEAR },
        { "rot-grid 3x3 (9)",  ROTGRID,  3, LINEAR },
        { "halton 9 spp",      HALTON,   3, LINEAR },
        { "regular 5x5 (25)",  REGULAR,  5, LINEAR },
        { "jittered 5x5 (25)", JITTERED, 5, LINEAR },
        { "halton 25 spp",     HALTON,   5, LINEAR },
    };
    printf("=== %s  (%dx%d, scale=%.3g, mxit=%d)  [ground truth = Halton 1024 spp] ===\n", v.name, W, H, v.scale, v.mxit);
    for (const Cand& c : cands) {
        double mse = 0; double px[3];
        for (int i = 0; i < H; ++i)
            for (int j = 0; j < W; ++j) {
                samplePixel(x0 + j * dx, y0 + i * dy, dx, dy, v.mxit, c.pat, c.n, j, i, px, c.avg);
                for (int k = 0; k < 3; ++k) { double d = px[k] - gt[(i * W + j) * 3 + k]; mse += d * d; }
            }
        mse /= (W * H * 3);
        double psnr = mse > 0 ? 10 * log10(255.0 * 255.0 / mse) : 99;
        printf("  %-20s  PSNR %6.2f dB   MSE %8.3f\n", c.name, psnr, mse);
    }
    // Adaptive: variable samples per pixel until the estimate converges.
    for (double tol : { 0.02, 0.01, 0.005 }) {
        double mse = 0; double px[3]; long tot = 0;
        for (int i = 0; i < H; ++i)
            for (int j = 0; j < W; ++j) {
                int nUsed;
                adaptivePixel(x0 + j * dx, y0 + i * dy, dx, dy, v.mxit, tol, 64, j, i, px, &nUsed);
                tot += nUsed;
                for (int k = 0; k < 3; ++k) { double d = px[k] - gt[(i * W + j) * 3 + k]; mse += d * d; }
            }
        mse /= (W * H * 3);
        double psnr = mse > 0 ? 10 * log10(255.0 * 255.0 / mse) : 99;
        printf("  adaptive tol=%.3f     PSNR %6.2f dB   MSE %8.3f   avg %.1f spp\n",
               tol, psnr, mse, (double)tot / (W * H));
    }
    printf("\n");
}

int main() {
    int W = 200, H = 150;
    View views[] = {
        { "seahorse valley boundary", -0.745,     0.113,    2000.0,  2000 },
        { "spiral filaments",         -0.761574, -0.0847596, 8000.0,  3000 },
        { "deep boundary dust",       -1.7687789, 0.0017389, 50000.0, 5000 },
    };
    for (const View& v : views) run(v, W, H);
    return 0;
}

