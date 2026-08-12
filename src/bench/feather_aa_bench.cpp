// Quantitative Feather (SAC) anti-aliasing validation.
//
// Renders a high-sub uniformly sampled engine grid as ground truth, then compares:
//   point.bmp       - centre sample
//   old1d.bmp       - previous |gradient|-wide 1D palette interval
//   box2d.bmp       - exact 2D affine phase / box-pixel palette integral
//   adaptive5.bmp   - real 5x samples on interior/exterior edges + box2d elsewhere
//
// Usage:
//   feather_aa_bench cx cy scale mxit density [W H GTsub]

#include <gmp.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "mandel_perturbation.h"
#include "interpolate.h"

static constexpr int P = 2048;
static float palette[3][P];
static double pal1[3][P + 1], pal2[3][P + 1], palMean[3];
static double density = 60.0;

static double srgbToLinear(double v) {
    double c = std::clamp(v / 255.0, 0.0, 1.0);
    return c <= 0.04045 ? c / 12.92 : std::pow((c + 0.055) / 1.055, 2.4);
}
static double linearToSrgb(double v) {
    v = std::clamp(v, 0.0, 1.0);
    return (v <= 0.0031308 ? 12.92 * v : 1.055 * std::pow(v, 1.0 / 2.4) - 0.055) * 255.0;
}
static void initPalette() {
    static const float pos[] = {0.0f, 0.16f, 0.42f, 0.6425f, 0.8575f, 1.0f};
    static const float col[3][6] = {{0, 32, 237, 255, 0, 0}, {70, 107, 255, 170, 2, 70}, {100, 203, 255, 0, 0, 100}};
    for (int c = 0; c < 3; ++c) {
        mono_cubic_interpolate(pos, col[c], 6, palette[c], P);
        pal1[c][0] = 0;
        for (int i = 0; i < P; ++i) pal1[c][i + 1] = pal1[c][i] + srgbToLinear(palette[c][i]);
        pal2[c][0] = 0;
        for (int i = 0; i < P; ++i) pal2[c][i + 1] = pal2[c][i] + 0.5 * (pal1[c][i] + pal1[c][i + 1]);
        palMean[c] = pal1[c][P] / P;
    }
}
static double phase(float sac) {
    return sac * (density / 20.0) * P;
}
static int palIndex(double u) {
    int i = (int)std::floor(u) % P;
    return i < 0 ? i + P : i;
}
static void pointColor(float v, double rgb[3]) {
    if (v < 0) {
        rgb[0] = rgb[1] = rgb[2] = 0;
        return;
    }
    int i = palIndex(phase(v));
    for (int c = 0; c < 3; ++c) rgb[c] = palette[c][i];
}
static double prefix1(int c, double x) {
    double q = std::floor(x / P), r = x - q * P;
    int i = (int)r;
    if (i >= P) i = P - 1;
    double t = r - i;
    return q * pal1[c][P] + pal1[c][i] * (1.0 - t) + pal1[c][i + 1] * t;
}
static double prefix2(int c, double x) {
    double A = pal1[c][P], B = pal2[c][P];
    double q = std::floor(x / P), r = x - q * P;
    int i = (int)r;
    if (i >= P) i = P - 1;
    double t = r - i, d = pal1[c][i + 1] - pal1[c][i];
    double within = pal2[c][i] + pal1[c][i] * t + 0.5 * d * t * t;
    return A * P * q * (q - 1.0) * 0.5 + q * B + q * A * r + within;
}
static double box2d(int c, double center, double dx, double dy) {
    double wx = std::abs(dx), wy = std::abs(dy);
    constexpr double eps = 1e-4;
    if (wx < eps && wy < eps) return srgbToLinear(palette[c][palIndex(center)]);
    if (wx < eps) return (prefix1(c, center + wy * 0.5) - prefix1(c, center - wy * 0.5)) / wy;
    if (wy < eps) return (prefix1(c, center + wx * 0.5) - prefix1(c, center - wx * 0.5)) / wx;
    double ax = wx * 0.5, ay = wy * 0.5;
    return std::clamp((prefix2(c, center + ax + ay) - prefix2(c, center - ax + ay) - prefix2(c, center + ax - ay) + prefix2(c, center - ax - ay)) / (wx * wy), 0.0, 1.0);
}
static bool ext(float v) {
    return v != EMPTYPIXEL && v >= 0;
}
static void gradients(float v, float l, float r, float u, float d, double& dx, double& dy) {
    double f = phase(v);
    if (ext(l) && ext(r))
        dx = 0.5 * (phase(r) - phase(l));
    else if (ext(r))
        dx = phase(r) - f;
    else if (ext(l))
        dx = f - phase(l);
    else
        dx = 0;
    if (ext(u) && ext(d))
        dy = 0.5 * (phase(d) - phase(u));
    else if (ext(d))
        dy = phase(d) - f;
    else if (ext(u))
        dy = f - phase(u);
    else
        dy = 0;
}
static void old1dColor(float v, float l, float r, float u, float d, double rgb[3]) {
    if (v < 0) {
        rgb[0] = rgb[1] = rgb[2] = 0;
        return;
    }
    double dx, dy;
    gradients(v, l, r, u, d, dx, dy);
    double width = std::sqrt(dx * dx + dy * dy), center = phase(v);
    for (int c = 0; c < 3; ++c) {
        double lin;
        if (width < 1)
            lin = srgbToLinear(palette[c][palIndex(center)]);
        else if (width >= P)
            lin = palMean[c];
        else
            lin = (prefix1(c, center + width * 0.5) - prefix1(c, center - width * 0.5)) / width;
        rgb[c] = linearToSrgb(lin);
    }
}
static void box2dColor(float v, float l, float r, float u, float d, double rgb[3]) {
    if (v < 0) {
        rgb[0] = rgb[1] = rgb[2] = 0;
        return;
    }
    double dx, dy;
    gradients(v, l, r, u, d, dx, dy);
    double center = phase(v);
    for (int c = 0; c < 3; ++c) rgb[c] = linearToSrgb(box2d(c, center, dx, dy));
}
static void quadColor(float v, float l, float r, float u, float d, float ul, float ur, float dl, float dr, int sub, double rgb[3]) {
    if (v < 0) {
        rgb[0] = rgb[1] = rgb[2] = 0;
        return;
    }
    double f = phase(v), gx, gy;
    gradients(v, l, r, u, d, gx, gy);
    double hxx = ext(l) && ext(r) ? phase(r) - 2.0 * f + phase(l) : 0.0;
    double hyy = ext(u) && ext(d) ? phase(d) - 2.0 * f + phase(u) : 0.0;
    double hxy = ext(ul) && ext(ur) && ext(dl) && ext(dr) ? 0.25 * (phase(dr) - phase(dl) - phase(ur) + phase(ul)) : 0.0;
    double lin[3] = {};
    for (int iy = 0; iy < sub; ++iy)
        for (int ix = 0; ix < sub; ++ix) {
            double x = (ix + 0.5) / sub - 0.5, y = (iy + 0.5) / sub - 0.5;
            double center = f + gx * x + gy * y + 0.5 * hxx * x * x + hxy * x * y + 0.5 * hyy * y * y;
            double dx = (gx + hxx * x + hxy * y) / sub;
            double dy = (gy + hxy * x + hyy * y) / sub;
            for (int c = 0; c < 3; ++c) lin[c] += box2d(c, center, dx, dy);
        }
    for (int c = 0; c < 3; ++c) rgb[c] = linearToSrgb(lin[c] / (sub * sub));
}
static void blockColor(const float* it, int W, int sub, int y, int x, double rgb[3]) {
    int stride = W * sub, cc = sub / 2;
    int center = (y * sub + cc) * stride + x * sub + cc;
    if (it[center - (stride + 1) * cc] == EMPTYPIXEL) {
        pointColor(it[center], rgb);
        return;
    }
    double sum[3] = {};
    for (int a = -cc; a <= cc; ++a)
        for (int b = -cc; b <= cc; ++b) {
            double p[3];
            pointColor(it[center + a * stride + b], p);
            for (int c = 0; c < 3; ++c) sum[c] += srgbToLinear(p[c]);
        }
    for (int c = 0; c < 3; ++c) rgb[c] = linearToSrgb(sum[c] / (sub * sub));
}
static void uniform3FromGt(const float* it, int W, int gtSub, int y, int x, double rgb[3]) {
    int stride = W * gtSub, cc = gtSub / 2, step = gtSub / 3;
    int center = (y * gtSub + cc) * stride + x * gtSub + cc;
    double sum[3] = {};
    for (int a = -1; a <= 1; ++a)
        for (int b = -1; b <= 1; ++b) {
            double p[3];
            pointColor(it[center + a * step * stride + b * step], p);
            for (int c = 0; c < 3; ++c) sum[c] += srgbToLinear(p[c]);
        }
    for (int c = 0; c < 3; ++c) rgb[c] = linearToSrgb(sum[c] / 9.0);
}
static void gauss2FromGt(const float* it, int W, int gtSub, int y, int x, double rgb[3]) {
    int stride = W * gtSub, cc = gtSub / 2;
    int step = std::max(1, (int)std::lround(gtSub / (2.0 * std::sqrt(3.0))));
    int center = (y * gtSub + cc) * stride + x * gtSub + cc;
    double sum[3] = {};
    for (int a : {-step, step})
        for (int b : {-step, step}) {
            double p[3];
            pointColor(it[center + a * stride + b], p);
            for (int c = 0; c < 3; ++c) sum[c] += srgbToLinear(p[c]);
        }
    for (int c = 0; c < 3; ++c) rgb[c] = linearToSrgb(sum[c] * 0.25);
}
static void gauss3FromGt(const float* it, int W, int gtSub, int y, int x, double rgb[3]) {
    int stride = W * gtSub, cc = gtSub / 2;
    int step = std::max(1, (int)std::lround(gtSub * std::sqrt(3.0 / 5.0) * 0.5));
    int center = (y * gtSub + cc) * stride + x * gtSub + cc;
    const int off[3] = {-step, 0, step};
    const double wt[3] = {5.0 / 18.0, 8.0 / 18.0, 5.0 / 18.0};
    double sum[3] = {};
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b) {
            double p[3];
            pointColor(it[center + off[a] * stride + off[b]], p);
            for (int c = 0; c < 3; ++c) sum[c] += wt[a] * wt[b] * srgbToLinear(p[c]);
        }
    for (int c = 0; c < 3; ++c) rgb[c] = linearToSrgb(sum[c]);
}
static void corner3FromGt(const float* it, int W, int gtSub, int y, int x, double rgb[3], double& phaseSpan) {
    int stride = W * gtSub, cc = gtSub / 2, step = gtSub / 3;
    int center = (y * gtSub + cc) * stride + x * gtSub + cc;
    double sum[3] = {}, lo = phase(it[center]), hi = lo;
    for (int a : {-step, step})
        for (int b : {-step, step}) {
            float v = it[center + a * stride + b];
            double p[3];
            pointColor(v, p);
            for (int c = 0; c < 3; ++c) sum[c] += srgbToLinear(p[c]);
            if (v >= 0) {
                lo = std::min(lo, phase(v));
                hi = std::max(hi, phase(v));
            }
        }
    for (int c = 0; c < 3; ++c) rgb[c] = linearToSrgb(sum[c] * 0.25);
    phaseSpan = (hi - lo) / P; // palette cycles
}
static bool flagged(const float* it, int W, int sub, int y, int x) {
    int stride = W * sub, cc = sub / 2;
    int center = (y * sub + cc) * stride + x * sub + cc;
    return it[center - (stride + 1) * cc] != EMPTYPIXEL;
}
static void highres3Color(const float* it, int W, int y, int x, double rgb[3]) {
    int stride = W * 3;
    double sum[3] = {};
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b) {
            double p[3];
            pointColor(it[(y * 3 + a) * stride + x * 3 + b], p);
            for (int c = 0; c < 3; ++c) sum[c] += srgbToLinear(p[c]);
        }
    for (int c = 0; c < 3; ++c) rgb[c] = linearToSrgb(sum[c] / 9.0);
}
static void highres2Color(const float* it, int W, int y, int x, double rgb[3]) {
    int stride = W * 2;
    double sum[3] = {};
    for (int a = 0; a < 2; ++a)
        for (int b = 0; b < 2; ++b) {
            double p[3];
            pointColor(it[(y * 2 + a) * stride + x * 2 + b], p);
            for (int c = 0; c < 3; ++c) sum[c] += srgbToLinear(p[c]);
        }
    for (int c = 0; c < 3; ++c) rgb[c] = linearToSrgb(sum[c] * 0.25);
}
static void writeBmp(const char* path, const std::vector<uint8_t>& rgb, int W, int H) {
    int row = (W * 3 + 3) & ~3, bytes = row * H;
    uint8_t fh[14] = {'B', 'M'}, ih[40] = {};
    uint32_t fs = 54 + bytes, off = 54, hs = 40, w = W, h = H, ds = bytes;
    uint16_t planes = 1, bpp = 24;
    memcpy(fh + 2, &fs, 4);
    memcpy(fh + 10, &off, 4);
    memcpy(ih, &hs, 4);
    memcpy(ih + 4, &w, 4);
    memcpy(ih + 8, &h, 4);
    memcpy(ih + 12, &planes, 2);
    memcpy(ih + 14, &bpp, 2);
    memcpy(ih + 20, &ds, 4);
    FILE* f = fopen(path, "wb");
    fwrite(fh, 1, 14, f);
    fwrite(ih, 1, 40, f);
    std::vector<uint8_t> line(row);
    for (int y = H - 1; y >= 0; --y) {
        for (int x = 0; x < W; ++x) {
            const uint8_t* p = &rgb[((size_t)y * W + x) * 3];
            line[x * 3] = p[2];
            line[x * 3 + 1] = p[1];
            line[x * 3 + 2] = p[0];
        }
        fwrite(line.data(), 1, row, f);
    }
    fclose(f);
}
struct Error {
    double sum = 0, sum2 = 0, max = 0;
    std::vector<double> pixels;
    void add(const double a[3], const double b[3]) {
        double e = 0;
        for (int c = 0; c < 3; ++c) {
            double d = std::abs(a[c] - b[c]);
            sum += d;
            sum2 += d * d;
            max = std::max(max, d);
            e += d;
        }
        pixels.push_back(e / 3.0);
    }
    void print(const char* name) {
        std::sort(pixels.begin(), pixels.end());
        double n = pixels.size() * 3.0;
        size_t p95 = (size_t)(pixels.size() * 0.95);
        printf("%-12s mean=%7.3f  RMSE=%7.3f  p95=%7.3f  max=%7.3f\n", name, sum / n, std::sqrt(sum2 / n), pixels[p95], max);
    }
};

int main(int argc, char** argv) {
    if (argc < 6) {
        fprintf(stderr, "usage: feather_aa_bench cx cy scale mxit density [W H GTsub]\n");
        return 1;
    }
    const char* cx = argv[1];
    const char* cy = argv[2];
    const char* scale = argv[3];
    int mxit = atoi(argv[4]);
    density = atof(argv[5]);
    int W = argc > 6 ? atoi(argv[6]) : 240, H = argc > 7 ? atoi(argv[7]) : 160;
    int gtSub = argc > 8 ? atoi(argv[8]) : 9;
    if (!(gtSub & 1)) ++gtSub;
    int precision = (int)((std::max(0.0, std::log10(atof(scale))) + 30) * std::log2(10.0)) + 64;
    mpf_set_default_prec(precision);
    initPalette();

    mpf_t mcx, mcy, msc;
    mpf_init_set_str(mcx, cx, 10);
    mpf_init_set_str(mcy, cy, 10);
    mpf_init_set_str(msc, scale, 10);
    int method = ColoringMethod::STRIPE_AVERAGE | ColoringMethod::SUPER_SAMPLING;

    size_t ng = (size_t)W * H * gtSub * gtSub;
    float* gt = new float[ng];
    std::fill(gt, gt + ng, (float)EMPTYPIXEL);
    _putenv_s("MANDEL_SAC_SS_STEPS", "0.000001");
    {
        Mandel m(W, H, mxit, gtSub, gt);
        m.setPrecision(precision);
        m.setDensity((float)density);
        m.Compute(mcx, mcy, msc, mxit, method);
    }

    size_t na = (size_t)W * H * 25;
    float* ad = new float[na];
    std::fill(ad, ad + na, (float)EMPTYPIXEL);
    _putenv_s("MANDEL_SAC_SS_STEPS", "0");
    {
        Mandel m(W, H, mxit, 5, ad);
        m.setPrecision(precision);
        m.setDensity((float)density);
        m.Compute(mcx, mcy, msc, mxit, method);
    }

    size_t nh = (size_t)W * H * 9;
    float* hi = new float[nh];
    std::fill(hi, hi + nh, (float)EMPTYPIXEL);
    {
        Mandel m(W * 3, H * 3, mxit, 1, hi);
        m.setPrecision(precision);
        m.Compute(mcx, mcy, msc, mxit, ColoringMethod::STRIPE_AVERAGE);
    }
    size_t nh2 = (size_t)W * H * 4;
    float* hi2 = new float[nh2];
    std::fill(hi2, hi2 + nh2, (float)EMPTYPIXEL);
    {
        Mandel m(W * 2, H * 2, mxit, 1, hi2);
        m.setPrecision(precision);
        m.Compute(mcx, mcy, msc, mxit, ColoringMethod::STRIPE_AVERAGE);
    }

    std::vector<uint8_t> imgGt((size_t)W * H * 3), imgPoint(imgGt.size()), imgOld(imgGt.size()), img2d(imgGt.size()), imgQ2(imgGt.size()), imgQ4(imgGt.size()), imgQ8(imgGt.size()), imgG2(imgGt.size()), imgU3(imgGt.size()), imgG3(imgGt.size()), imgHi2(imgGt.size()), imgHi3(imgGt.size()), imgAd(imgGt.size());
    Error ep, eo, e2, eq2, eq4, eq8, eg2, eu3, eg3, ea;
    Error ehi2, ehi3;
    static const double progT[] = {0.05, 0.1, 0.25, 0.5, 1.0, 2.0};
    Error prog[6];
    long refined[6] = {};
    long adCount = 0;
    int gs = W * gtSub, gc = gtSub / 2, as = W * 5, ac = 2;
    for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x) {
            int gi = (y * gtSub + gc) * gs + x * gtSub + gc;
            float v = gt[gi], l = x ? gt[gi - gtSub] : EMPTYPIXEL, r = x + 1 < W ? gt[gi + gtSub] : EMPTYPIXEL, u = y ? gt[gi - gs * gtSub] : EMPTYPIXEL, d = y + 1 < H ? gt[gi + gs * gtSub] : EMPTYPIXEL;
            float ul = x && y ? gt[gi - gs * gtSub - gtSub] : EMPTYPIXEL, ur = x + 1 < W && y ? gt[gi - gs * gtSub + gtSub] : EMPTYPIXEL, dl = x && y + 1 < H ? gt[gi + gs * gtSub - gtSub] : EMPTYPIXEL, dr = x + 1 < W && y + 1 < H ? gt[gi + gs * gtSub + gtSub] : EMPTYPIXEL;
            double cg[3], cp[3], co[3], c2[3], cq2[3], cq4[3], cq8[3], cg2[3], cu3[3], cg3[3], chi2[3], chi3[3], cc3[3], ca[3], span;
            blockColor(gt, W, gtSub, y, x, cg);
            gauss2FromGt(gt, W, gtSub, y, x, cg2);
            uniform3FromGt(gt, W, gtSub, y, x, cu3);
            gauss3FromGt(gt, W, gtSub, y, x, cg3);
            highres2Color(hi2, W, y, x, chi2);
            highres3Color(hi, W, y, x, chi3);
            corner3FromGt(gt, W, gtSub, y, x, cc3, span);
            pointColor(v, cp);
            old1dColor(v, l, r, u, d, co);
            box2dColor(v, l, r, u, d, c2);
            quadColor(v, l, r, u, d, ul, ur, dl, dr, 2, cq2);
            quadColor(v, l, r, u, d, ul, ur, dl, dr, 4, cq4);
            quadColor(v, l, r, u, d, ul, ur, dl, dr, 8, cq8);
            int ai = (y * 5 + ac) * as + x * 5 + ac;
            if (flagged(ad, W, 5, y, x)) {
                blockColor(ad, W, 5, y, x, ca);
                ++adCount;
            } else {
                float av = ad[ai], al = x ? ad[ai - 5] : EMPTYPIXEL, ar = x + 1 < W ? ad[ai + 5] : EMPTYPIXEL, au = y ? ad[ai - as * 5] : EMPTYPIXEL, avd = y + 1 < H ? ad[ai + as * 5] : EMPTYPIXEL;
                box2dColor(av, al, ar, au, avd, ca);
            }
            ep.add(cp, cg);
            eo.add(co, cg);
            e2.add(c2, cg);
            eq2.add(cq2, cg);
            eq4.add(cq4, cg);
            eq8.add(cq8, cg);
            eg2.add(cg2, cg);
            eu3.add(cu3, cg);
            eg3.add(cg3, cg);
            ehi2.add(chi2, cg);
            ehi3.add(chi3, cg);
            ea.add(ca, cg);
            for (int k = 0; k < 6; ++k) {
                bool refine = span > progT[k];
                prog[k].add(refine ? cu3 : cc3, cg);
                if (refine) ++refined[k];
            }
            double* colors[] = {cg, cp, co, c2, cq2, cq4, cq8, cg2, cu3, cg3, chi2, chi3, ca};
            std::vector<uint8_t>* imgs[] = {&imgGt, &imgPoint, &imgOld, &img2d, &imgQ2, &imgQ4, &imgQ8, &imgG2, &imgU3, &imgG3, &imgHi2, &imgHi3, &imgAd};
            for (int k = 0; k < 13; ++k)
                for (int c = 0; c < 3; ++c) (*imgs[k])[((size_t)y * W + x) * 3 + c] = (uint8_t)std::clamp(colors[k][c] + 0.5, 0.0, 255.0);
        }
    printf("W=%d H=%d GT=%dx density=%.1f  adaptive5 flagged=%.2f%%\n", W, H, gtSub, density, 100.0 * adCount / (W * H));
    ep.print("point");
    eo.print("old1d");
    e2.print("box2d");
    eq2.print("quad2");
    eq4.print("quad4");
    eq8.print("quad8");
    eg2.print("gauss2");
    eu3.print("uniform3");
    eg3.print("gauss3");
    ehi2.print("highres2");
    ehi3.print("highres3");
    ea.print("adaptive5");
    puts("progressive centre+4 corners -> full 3x:");
    for (int k = 0; k < 6; ++k) {
        char name[32];
        snprintf(name, sizeof(name), "phase>%.2f", progT[k]);
        prog[k].print(name);
        printf("             refined=%6.2f%%  evals/pixel=%5.2f\n", 100.0 * refined[k] / (W * H), 5.0 + 4.0 * refined[k] / (W * H));
    }
    auto timeFilter = [&](int mode) {
        volatile double sink = 0;
        auto start = std::chrono::steady_clock::now();
        for (int rep = 0; rep < 5; ++rep)
            for (int y = 1; y + 1 < H; ++y)
                for (int x = 1; x + 1 < W; ++x) {
                    int gi = (y * gtSub + gc) * gs + x * gtSub + gc;
                    float v = gt[gi], l = gt[gi - gtSub], r = gt[gi + gtSub], u = gt[gi - gs * gtSub], d = gt[gi + gs * gtSub], ul = gt[gi - gs * gtSub - gtSub], ur = gt[gi - gs * gtSub + gtSub], dl = gt[gi + gs * gtSub - gtSub], dr = gt[gi + gs * gtSub + gtSub];
                    double out[3];
                    if (mode == 0)
                        old1dColor(v, l, r, u, d, out);
                    else if (mode == 1)
                        box2dColor(v, l, r, u, d, out);
                    else
                        quadColor(v, l, r, u, d, ul, ur, dl, dr, 4, out);
                    sink = sink + out[0];
                }
        double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count() / 5.0;
        return ms;
    };
    printf("color-only time/frame: old1d=%.2fms box2d=%.2fms quad4=%.2fms (single thread)\n", timeFilter(0), timeFilter(1), timeFilter(2));
    writeBmp("feather_gt.bmp", imgGt, W, H);
    writeBmp("feather_point.bmp", imgPoint, W, H);
    writeBmp("feather_old1d.bmp", imgOld, W, H);
    writeBmp("feather_box2d.bmp", img2d, W, H);
    writeBmp("feather_quad2.bmp", imgQ2, W, H);
    writeBmp("feather_quad4.bmp", imgQ4, W, H);
    writeBmp("feather_quad8.bmp", imgQ8, W, H);
    writeBmp("feather_gauss2.bmp", imgG2, W, H);
    writeBmp("feather_uniform3.bmp", imgU3, W, H);
    writeBmp("feather_gauss3.bmp", imgG3, W, H);
    writeBmp("feather_highres2.bmp", imgHi2, W, H);
    writeBmp("feather_highres3.bmp", imgHi3, W, H);
    writeBmp("feather_adaptive5.bmp", imgAd, W, H);
    delete[] gt;
    delete[] ad;
    delete[] hi;
    delete[] hi2;
    mpf_clears(mcx, mcy, msc, (mpf_ptr)0);
    return 0;
}
