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

#include "mandel_perturbation.h"
#include "interpolate.h"

static const int colP = 2048;
static float color_map[3][colP];
static const float color_density = 60.0f;

static void colorMapInit() {
    static const int N = 6;
    static const float pos[N] = { 0.0f, 0.16f, 0.42f, 0.6425f, 0.8575f, 1.0f };
    static const float col[3][N] = { {0,32,237,255,0,0}, {70,107,255,170,2,70}, {100,203,255,0,0,100} };
    for (int i = 0; i < 3; ++i)
        mono_cubic_interpolate(pos, col[i], N, color_map[i], colP);
}

static inline float color_func(float it) {
    return powf(logf(it + 2), logf(it + 2) * logf(it + 2) / color_density);
}

static void getColor(float it, float& r, float& g, float& b) {
    if (it < 0) { r = g = b = 0; return; }   // interior -> black
    float f = color_func(it);
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
    int scaleExp = demo ? 1000 : atoi(argv[6]);
    int mxit = demo ? 10000 : atoi(argv[7]);
    const int SS = demo ? 2 : (argc > 8 ? atoi(argv[8]) : 3);
    int Ws = W * SS, Hs = H * SS;

    std::string scale = "1"; for (int i = 0; i < scaleExp; ++i) scale += "0";
    int precision = (int)(scale.size() * log(10) / log(2)) + 64;
    mpf_set_default_prec(precision);

    colorMapInit();

    float* iter = new float[Ws * Hs];
    for (int i = 0; i < Ws * Hs; ++i) iter[i] = EMPTYPIXEL;
    Mandel mandel(Ws, Hs, mxit, 1, iter);
    mandel.setPrecision(precision);
    mpf_t mcx, mcy, msc;
    mpf_init_set_str(mcx, cx, 10);
    mpf_init_set_str(mcy, cy, 10);
    mpf_init_set_str(msc, scale.c_str(), 10);

    printf("Rendering %dx%d (SS=%d -> %dx%d), scale=1e%d, mxit=%d, prec=%d ...\n",
           W, H, SS, Ws, Hs, scaleExp, mxit, precision);
    mandel.Compute(mcx, mcy, msc, mxit, 0);

    // Colour hi-res, box-downsample SSxSS in linear light.
    std::vector<uint8_t> img(W * H * 3);
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
    printf("Wrote %s\n", out);
    return 0;
}
