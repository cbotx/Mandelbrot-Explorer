// Self-contained image renderer for eyeballing the engine output. Uses the real
// Mandel compute + the repo's palette (Ultra-Fractal control points + monotone
// cubic interpolation + the color_func mapping), renders at 3x supersampling and
// box-downsamples in linear light, writes a 24-bit BMP. No libpng/GLUT needed.
//
// Usage: render [<out.bmp> <W> <H> <cx> <cy> <scaleExp> <mxit> [SS]]
// With no arguments, renders a ready-to-view 1e1000 exterior demo.

#include <gmp.h>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <chrono>

#include "mandel_perturbation.h"
#include "interpolate.h"

static const int colP = 2048;
static float color_map[3][colP];
static float color_density = 60.0f;   // MANDEL_DENS overrides (comparison tooling)

static void colorMapInit() {
    static const int N = 6;
    static const float pos[N] = { 0.0f, 0.16f, 0.42f, 0.6425f, 0.8575f, 1.0f };
    static const float col[3][N] = { {0,32,237,255,0,0}, {70,107,255,170,2,70}, {100,203,255,0,0,100} };
    for (int i = 0; i < 3; ++i)
        mono_cubic_interpolate(pos, col[i], N, color_map[i], colP);
}

static inline float color_func(float it) {
    return powf(logf(it + 2), logf(it + 2) * logf(it + 2) * color_density / 3600.0f);
}

static int g_sac = 0;
static int g_ede = 0;

static void getColor(float it, float& r, float& g, float& b) {
    if (it < 0) { r = g = b = 0; return; }   // interior -> black
    // SAC value is already in [0,1]; map it around the palette a few times for
    // banded feather texture. Iteration counts use the log-power curve.
    float f = g_sac ? it * (color_density / 20.0f)
            : g_ede ? tanhf(it * color_density / 3600.0f * 5.0f)
            : color_func(it);
    int x = (int)(f * colP) % colP;
    if (x < 0) x += colP;
    r = color_map[0][x]; g = color_map[1][x]; b = color_map[2][x];
}

static inline float toLin(float c) { return powf(c / 255.f, 2.2f); }
static inline float toGam(float l) { return 255.f * powf(l < 0 ? 0 : l, 1.f / 2.2f); }

static void writeBMP(const char* path, const std::vector<uint8_t>& rgb, int W, int H) {
    int row = (W * 3 + 3) & ~3;               // 4-byte aligned rows
    int dataSize = row * H;
    uint8_t fh[14] = { 'B','M' };
    uint32_t fsize = 14 + 40 + dataSize;
    memcpy(fh + 2, &fsize, 4);
    uint32_t off = 54; memcpy(fh + 10, &off, 4);
    uint8_t ih[40] = { 0 }; uint32_t v;
    v = 40; memcpy(ih + 0, &v, 4);
    v = W; memcpy(ih + 4, &v, 4);
    v = H; memcpy(ih + 8, &v, 4);
    uint16_t planes = 1, bpp = 24; memcpy(ih + 12, &planes, 2); memcpy(ih + 14, &bpp, 2);
    v = dataSize; memcpy(ih + 20, &v, 4);
    FILE* f = fopen(path, "wb");
    fwrite(fh, 1, 14, f); fwrite(ih, 1, 40, f);
    std::vector<uint8_t> line(row, 0);
    for (int y = H - 1; y >= 0; --y) {          // BMP is bottom-up
        for (int x = 0; x < W; ++x) {
            const uint8_t* p = &rgb[(y * W + x) * 3];
            line[x * 3 + 0] = p[2]; line[x * 3 + 1] = p[1]; line[x * 3 + 2] = p[0]; // BGR
        }
        fwrite(line.data(), 1, row, f);
    }
    fclose(f);
}

int main(int argc, char** argv) {
    if (argc != 1 && argc < 8) {
        fprintf(stderr, "usage: render [out.bmp W H cx cy scaleExp mxit [SS]]\n");
        return 1;
    }
    const bool demo = argc == 1;
    const char* out = demo ? "deep_exterior_1e1000.bmp" : argv[1];
    int W = demo ? 960 : atoi(argv[2]);
    int H = demo ? 540 : atoi(argv[3]);
    const char* cx = demo ? "0" : argv[4];
    const char* cy = demo ? "1" : argv[5];
    const char* scaleArg = demo ? "1000" : argv[6];
    int mxit = demo ? 10000 : atoi(argv[7]);
    const int SS = demo ? 2 : (argc > 8 ? atoi(argv[8]) : 3);
    int Ws = W * SS, Hs = H * SS;

    // scaleArg is either a power-of-10 exponent ("51") or a decimal magnitude
    // ("3.831277e51"); expand either to an integer digit string for GMP.
    std::string sa = scaleArg, scale;
    if (sa.find('e') != std::string::npos || sa.find('E') != std::string::npos || sa.find('.') != std::string::npos) {
        size_t ep = sa.find_first_of("eE");
        std::string mant = (ep == std::string::npos) ? sa : sa.substr(0, ep);
        int e10 = (ep == std::string::npos) ? 0 : atoi(sa.substr(ep + 1).c_str());
        std::string sign;
        if (!mant.empty() && (mant[0] == '-' || mant[0] == '+')) { sign = (mant[0] == '-') ? "-" : ""; mant.erase(0, 1); }
        size_t dp = mant.find('.');
        int frac = (dp == std::string::npos) ? 0 : (int)(mant.size() - dp - 1);
        if (dp != std::string::npos) mant.erase(dp, 1);
        int zeros = e10 - frac;
        scale = sign + mant; if (zeros > 0) scale.append(zeros, '0');
    } else {
        int scaleExp = atoi(sa.c_str());
        scale = "1"; for (int i = 0; i < scaleExp; ++i) scale += "0";
    }
    int precision = (int)(scale.size() * log(10) / log(2)) + 64;
    mpf_set_default_prec(precision);

    colorMapInit();
    if (getenv("MANDEL_DENS")) color_density = (float)atof(getenv("MANDEL_DENS"));

    float* iter = new float[Ws * Hs];
    for (int i = 0; i < Ws * Hs; ++i) iter[i] = EMPTYPIXEL;
    Mandel mandel(Ws, Hs, mxit, 1, iter);
    mandel.setPrecision(precision);
    mpf_t mcx, mcy, msc;
    mpf_init_set_str(mcx, cx, 10);
    mpf_init_set_str(mcy, cy, 10);
    mpf_init_set_str(msc, scale.c_str(), 10);

    printf("Rendering %dx%d (SS=%d -> %dx%d), scale=%s (%zu digits), mxit=%d, prec=%d ...\n",
           W, H, SS, Ws, Hs, scaleArg, scale.size(), mxit, precision);
    int cmethod = (getenv("MANDEL_EDE") && atoi(getenv("MANDEL_EDE"))) ? ColoringMethod::EXTERIOR_DIST_EST : 0;
    if (cmethod & ColoringMethod::EXTERIOR_DIST_EST) g_ede = 1;
    if (getenv("MANDEL_SAC") && atoi(getenv("MANDEL_SAC"))) { cmethod |= ColoringMethod::STRIPE_AVERAGE; g_sac = 1; }
    mandel.Compute(mcx, mcy, msc, mxit, cmethod);
    if (getenv("MANDEL_STRIPS") && atoi(getenv("MANDEL_STRIPS")) > 0) {
        // Validation: re-render in horizontal strips (each a strip-sized Mandel
        // covering rows [base,base+sh) of the full Hs-row image) and overwrite
        // iter. Should match the single Compute above.
        int nstrips = atoi(getenv("MANDEL_STRIPS"));
        int sh = (Hs + nstrips - 1) / nstrips;
        for (int base = 0; base < Hs; base += sh) {
            int h = std::min(sh, Hs - base);
            std::vector<float> sbuf((size_t)Ws * h, EMPTYPIXEL);
            Mandel strip(Ws, h, mxit, 1, sbuf.data());
            strip.setPrecision(precision);
            strip.Compute(mcx, mcy, msc, mxit, cmethod, Hs, base);
            for (int r = 0; r < h; ++r)
                for (int c = 0; c < Ws; ++c) iter[(size_t)(base + r) * Ws + c] = sbuf[(size_t)r * Ws + c];
        }
    }
    if (getenv("MANDEL_ORACLE")) mandel.ComputeDirect(mxit, iter, 1, cmethod);   // brute-force GMP ground truth

    if (getenv("MANDEL_DUMPRAW")) {   // raw per-pixel value (float) for exact diffs
        std::string rp = std::string(out) + ".raw";
        FILE* rf = fopen(rp.c_str(), "wb");
        if (rf) { fwrite(iter, sizeof(float), (size_t)Ws * Hs, rf); fclose(rf); }
    }

    // Colour hi-res, box-downsample SSxSS in linear light.
    std::vector<uint8_t> img(W * H * 3);
    auto _tc0 = std::chrono::high_resolution_clock::now();
#pragma omp parallel for schedule(dynamic, 8)
    for (int i = 0; i < H; ++i)
        for (int j = 0; j < W; ++j) {
            float lr = 0, lg = 0, lb = 0;
            for (int a = 0; a < SS; ++a)
                for (int b = 0; b < SS; ++b) {
                    float r, g, bb; getColor(iter[(i * SS + a) * Ws + (j * SS + b)], r, g, bb);
                    lr += toLin(r); lg += toLin(g); lb += toLin(bb);
                }
            int n = SS * SS;
            uint8_t* p = &img[(i * W + j) * 3];
            p[0] = (uint8_t)(toGam(lr / n) + 0.5f);
            p[1] = (uint8_t)(toGam(lg / n) + 0.5f);
            p[2] = (uint8_t)(toGam(lb / n) + 0.5f);
        }
    writeBMP(out, img, W, H);
    if (getenv("MANDEL_PROFILE")) {
        double _tcol = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - _tc0).count();
        fprintf(stderr, "  [timing] coloring pass: %.4f s  (%dx%d, SS=%d)\n", _tcol, W, H, SS);
    }
    printf("Wrote %s\n", out);
    return 0;
}
