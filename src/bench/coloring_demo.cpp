// Coloring-method demo (brute-force double, shallow views only) to compare
// visual styles without touching the production renderer. Reuses the GUI's LAB
// palette rebuild. Not for deep zoom -- this is an experiment/eyeballing tool.
//
// Usage: coloring_demo out.bmp W H cx cy scale mxit palIdx mode [SS]
//   mode: 0 smooth, 1 stripe average, 2 analytical normal lighting,
//         3 point trap, 4 cross trap, 5 circle trap, 6 composite trap

#include <cstdio>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <chrono>
#include <string>
#include <vector>
#include <array>
#include "interpolate.h"
#include "color.h"

static const int CP = 2048;
static float color_map[3][CP];
struct Stop {
    float pos, r, g, b;
};
static std::vector<std::vector<Stop>> pals = {
    {{0, 0, 70, 100}, {0.16f, 32, 107, 203}, {0.42f, 237, 255, 255}, {0.6425f, 255, 170, 0}, {0.8575f, 0, 2, 0}},                                                                        // sunrise
    {{0.05f, 30, 90, 200}, {0.17f, 242, 248, 255}, {0.30f, 25, 155, 95}, {0.42f, 10, 12, 20}, {0.55f, 125, 70, 185}, {0.67f, 250, 244, 252}, {0.80f, 240, 180, 55}, {0.92f, 18, 12, 8}}, // gems
    {{0, 4, 12, 38}, {0.22f, 10, 64, 120}, {0.46f, 36, 150, 190}, {0.66f, 150, 220, 230}, {0.85f, 244, 252, 255}},                                                                       // ocean
};
static void rebuild(const std::vector<Stop>& st) {
    int n = (int)st.size();
    int m = 3 * n;
    std::vector<float> xs(m), out(CP);
    std::vector<std::array<float, 3>> lab(n);
    for (int i = 0; i < n; ++i) rgb2lab(st[i].r, st[i].g, st[i].b, lab[i][0], lab[i][1], lab[i][2]);
    for (int p = -1; p <= 1; ++p)
        for (int i = 0; i < n; ++i) xs[(p + 1) * n + i] = st[i].pos + p;
    for (int c = 0; c < 3; ++c) {
        std::vector<float> ys(m);
        for (int p = 0; p < 3; ++p)
            for (int i = 0; i < n; ++i) ys[p * n + i] = lab[i][c];
        mono_cubic_interpolate(xs.data(), ys.data(), m, out.data(), CP);
        for (int i = 0; i < CP; ++i) color_map[c][i] = out[i];
    }
    for (int i = 0; i < CP; ++i) {
        float r, g, b;
        lab2rgb(color_map[0][i], color_map[1][i], color_map[2][i], r, g, b);
        color_map[0][i] = r;
        color_map[1][i] = g;
        color_map[2][i] = b;
    }
}
static void palColor(float t, float& r, float& g, float& b) { // t in [0,1)
    int x = (int)(t * CP) % CP;
    if (x < 0) x += CP;
    r = color_map[0][x];
    g = color_map[1][x];
    b = color_map[2][x];
}
static inline float toLin(float c) {
    return powf(c / 255.f, 2.2f);
}
static inline float toGam(float l) {
    return 255.f * powf(l < 0 ? 0 : l, 1 / 2.2f);
}

static void smoothPalette(double mu, float& R, float& G, float& B) {
    float l = logf((float)std::max(mu, 0.0) + 2);
    float t = powf(l, l * l / 60.f);
    palColor(t - floorf(t), R, G, B);
}

static void pixel(double cr, double ci, int mxit, int mode, float& R, float& G, float& B) {
    double zr = 0, zi = 0, dr = 0, di = 0, sum = 0, lastAdd = 0;
    int cnt = 0;
    double minPoint = 1e300, minCross = 1e300, minCircle = 1e300;
    double trapAngle = 0;
    const bool doStripe = mode == 1, doNormal = mode == 2, doTrap = mode >= 3;
    const double freq = 7.0;
    const double R2 = 1e20, logR = log(sqrt(R2));
    for (int i = 0; i < mxit; ++i) {
        double ndr = 0, ndi = 0;
        if (doNormal) {
            // d(z^2+c)/dc = 2*z*dz/dc + 1 (analytical normal lighting).
            ndr = 2 * (zr * dr - zi * di) + 1;
            ndi = 2 * (zr * di + zi * dr);
        }
        double zr2 = zr * zr, zi2 = zi * zi;
        double nzr = zr2 - zi2 + cr, nzi = 2 * zr * zi + ci;
        zr = nzr;
        zi = nzi;
        if (doNormal) {
            dr = ndr;
            di = ndi;
        }
        double m2 = zr * zr + zi * zi;
        if (doTrap) {
            double point = m2, cross = std::min(std::fabs(zr), std::fabs(zi));
            double circle = std::fabs(std::sqrt(m2) - 0.5);
            if (point < minPoint) {
                minPoint = point;
                trapAngle = atan2(zi, zr);
            }
            minCross = std::min(minCross, cross);
            minCircle = std::min(minCircle, circle);
        }
        if (doStripe) {
            double add = 0.5 + 0.5 * sin(freq * atan2(zi, zr));
            lastAdd = add;
            sum += add;
            ++cnt;
        }
        if (m2 > R2) {
            double mu = i + 1 - log(log(m2) / 2 / log(2)) / log(2);
            if (mode == 1) {
                double aw = sum / cnt, awo = cnt > 1 ? (sum - lastAdd) / (cnt - 1) : aw;
                double frac = 1.0 - log(log(sqrt(m2)) / logR) / log(2.0);
                if (frac < 0) frac = 0;
                if (frac > 1) frac = 1;
                double avg = awo + (aw - awo) * frac;
                palColor((float)fmod(avg, 1.0), R, G, B);
                return;
            }
            if (mode == 2) {
                smoothPalette(mu, R, G, B);
                // Complex u=z/(dz/dc) points along the exterior-potential normal.
                double den = dr * dr + di * di;
                if (den > 0) {
                    double ux = (zr * dr + zi * di) / den, uy = (zi * dr - zr * di) / den;
                    double un = hypot(ux, uy);
                    if (un > 0) {
                        ux /= un;
                        uy /= un;
                    }
                    const double nz = 1.0, lx = -0.45, ly = -0.35, lz = 0.82;
                    double nn = sqrt(1 + nz * nz), ln = sqrt(lx * lx + ly * ly + lz * lz);
                    double diffuse = std::max(0.0, (ux * lx + uy * ly + nz * lz) / (nn * ln));
                    double shade = 0.35 + 0.75 * diffuse;
                    R = (float)std::min(255.0, R * shade);
                    G = (float)std::min(255.0, G * shade);
                    B = (float)std::min(255.0, B * shade);
                }
                return;
            }
            if (mode >= 3) {
                double d = mode == 3 ? sqrt(minPoint) : mode == 4 ? minCross
                                                    : mode == 5   ? minCircle
                                                                  : std::min({sqrt(minPoint), minCross * 1.5, minCircle});
                // Log distance exposes nested trap contours; mix a little smooth
                // dwell and the angle at closest approach to avoid flat plateaus.
                double trap = -log10(std::max(d, 1e-14));
                double t = 0.17 * trap + 0.025 * mu + (mode == 6 ? 0.10 * trapAngle : 0.0);
                palColor((float)(t - floor(t)), R, G, B);
                double glow = 0.38 + 0.62 * exp(-3.0 * std::min(d, 1.0));
                R = (float)(R * glow);
                G = (float)(G * glow);
                B = (float)(B * glow);
                return;
            }
            smoothPalette(mu, R, G, B);
            return;
        }
    }
    R = G = B = 0;
}

static void writeBMP(const char* p, const std::vector<uint8_t>& rgb, int W, int H) {
    int row = (W * 3 + 3) & ~3, ds = row * H;
    uint8_t fh[14] = {'B', 'M'};
    uint32_t fs = 14 + 40 + ds;
    memcpy(fh + 2, &fs, 4);
    uint32_t off = 54;
    memcpy(fh + 10, &off, 4);
    uint8_t ih[40] = {0};
    uint32_t v;
    v = 40;
    memcpy(ih, &v, 4);
    v = W;
    memcpy(ih + 4, &v, 4);
    v = H;
    memcpy(ih + 8, &v, 4);
    uint16_t pl = 1, bp = 24;
    memcpy(ih + 12, &pl, 2);
    memcpy(ih + 14, &bp, 2);
    v = ds;
    memcpy(ih + 20, &v, 4);
    FILE* f = fopen(p, "wb");
    fwrite(fh, 1, 14, f);
    fwrite(ih, 1, 40, f);
    std::vector<uint8_t> ln(row, 0);
    for (int y = H - 1; y >= 0; --y) {
        for (int x = 0; x < W; ++x) {
            const uint8_t* q = &rgb[(y * W + x) * 3];
            ln[x * 3] = q[2];
            ln[x * 3 + 1] = q[1];
            ln[x * 3 + 2] = q[0];
        }
        fwrite(ln.data(), 1, row, f);
    }
    fclose(f);
}

int main(int argc, char** argv) {
    if (argc < 10) {
        fprintf(stderr, "usage: coloring_demo out W H cx cy scale mxit palIdx mode [SS]\n");
        return 1;
    }
    const char* out = argv[1];
    int W = atoi(argv[2]), H = atoi(argv[3]);
    double cx = atof(argv[4]), cy = atof(argv[5]), scale = atof(argv[6]);
    int mxit = atoi(argv[7]);
    int palIdx = atoi(argv[8]), mode = atoi(argv[9]);
    int SS = argc > 10 ? atoi(argv[10]) : 3;
    if (palIdx < 0 || palIdx >= (int)pals.size()) palIdx = 0;
    rebuild(pals[palIdx]);
    const char* names[] = {"smooth", "stripe-average", "normal-light", "point-trap", "cross-trap", "circle-trap", "composite-trap"};
    if (mode < 0 || mode > 6) {
        fprintf(stderr, "mode must be 0..6\n");
        return 1;
    }
    int Ws = W * SS, Hs = H * SS;
    double half = 2.0 / scale; // view half-width in complex units (approx)
    std::vector<uint8_t> img(W * H * 3);
    auto t0 = std::chrono::steady_clock::now();
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < H; ++i)
        for (int j = 0; j < W; ++j) {
            float lr = 0, lg = 0, lb = 0;
            for (int a = 0; a < SS; ++a)
                for (int b = 0; b < SS; ++b) {
                    double px = (j * SS + b + 0.5) / Ws, py = (i * SS + a + 0.5) / Hs;
                    double cr = cx + (px * 2 - 1) * half;
                    double ci = cy - (py * 2 - 1) * half * H / W;
                    float R, G, B;
                    pixel(cr, ci, mxit, mode, R, G, B);
                    lr += toLin(R);
                    lg += toLin(G);
                    lb += toLin(B);
                }
            int n = SS * SS;
            uint8_t* p = &img[(i * W + j) * 3];
            p[0] = (uint8_t)(toGam(lr / n) + 0.5f);
            p[1] = (uint8_t)(toGam(lg / n) + 0.5f);
            p[2] = (uint8_t)(toGam(lb / n) + 0.5f);
        }
    double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    writeBMP(out, img, W, H);
    printf("wrote %s mode=%s compute+color=%.3f ms\n", out, names[mode], ms);
    return 0;
}
