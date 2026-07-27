// Palette preview: render one view with a chosen candidate palette, using the
// SAME LAB-space mono-cubic rebuild + non-EDE colorFunction as the GUI, so the
// output matches what the app would show. For comparing palette options.
//
// Usage: palette_preview out.bmp W H cx cy scaleExp mxit paletteIdx [density] [SS]

#include <gmp.h>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <array>

#include "mandel_perturbation.h"
#include "interpolate.h"
#include "color.h"

static const int CP = 2048;
static float color_map[3][CP];
static float color_density = 60.0f;

struct Stop { float pos, r, g, b; };

// Candidate palettes (name + stops). Index selects one.
struct Pal { const char* name; std::vector<Stop> stops; };
static std::vector<Pal> palettes = {
    // Natural / gemstone colours. Saturated hues bridged by white or near-black.
    { "reef",     { {0,0,45,65},{0.15f,25,155,195},{0.40f,240,255,255},{0.60f,255,190,45},{0.78f,220,60,40},{0.90f,12,3,10} } },
    { "emerald",  { {0,6,45,30},{0.16f,40,160,90},{0.42f,235,255,235},{0.6425f,225,205,120},{0.8575f,6,18,12} } },
    { "amethyst", { {0,32,12,60},{0.18f,120,70,180},{0.40f,200,175,230},{0.58f,245,242,252},{0.8575f,18,8,30} } },
    { "topaz",    { {0,35,20,8},{0.18f,170,95,30},{0.40f,235,175,70},{0.60f,250,235,200},{0.85f,20,12,6} } },
    // Gemstone medleys: 4 jewel hues, alternating white/black bridges (2 white, 2 black).
    { "gems",     { {0.05f,30,90,200},{0.17f,242,248,255},{0.30f,25,155,95},{0.42f,10,12,20},{0.55f,125,70,185},{0.67f,250,244,252},{0.80f,240,180,55},{0.92f,18,12,8} } },
    { "gemsB",    { {0.05f,30,90,200},{0.17f,242,248,255},{0.30f,25,155,95},{0.42f,10,12,20},{0.55f,205,45,75},{0.67f,250,244,252},{0.80f,240,180,55},{0.92f,18,12,8} } },
};

static void rebuild(const std::vector<Stop>& stops) {
    int n = (int)stops.size();
    std::vector<std::array<float, 3>> lab(n);
    for (int i = 0; i < n; ++i)
        rgb2lab(stops[i].r, stops[i].g, stops[i].b, lab[i][0], lab[i][1], lab[i][2]);
    // Tile over three periods so the spline is periodic across the cycle wrap
    // (matches PaletteEditor::rebuild in the GUI).
    int m = 3 * n;
    std::vector<float> xs(m), out(CP);
    for (int p = -1; p <= 1; ++p)
        for (int i = 0; i < n; ++i)
            xs[(p + 1) * n + i] = stops[i].pos + p;
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
        color_map[0][i] = r; color_map[1][i] = g; color_map[2][i] = b;
    }
}

static inline float colorFunction(float it) {   // non-EDE, matches win32_main
    float l = logf(it + 2.0f);
    return powf(l, l * l / color_density);
}
static void getColor(float it, float& r, float& g, float& b) {
    if (it < 0) { r = g = b = 0; return; }
    int x = (int)(colorFunction(it) * CP) % CP; if (x < 0) x += CP;
    r = color_map[0][x]; g = color_map[1][x]; b = color_map[2][x];
}
static inline float toLin(float c) { return powf(c / 255.f, 2.2f); }
static inline float toGam(float l) { return 255.f * powf(l < 0 ? 0 : l, 1.f / 2.2f); }

static void writeBMP(const char* path, const std::vector<uint8_t>& rgb, int W, int H) {
    int row = (W * 3 + 3) & ~3, dataSize = row * H;
    uint8_t fh[14] = { 'B','M' }; uint32_t fsize = 14 + 40 + dataSize;
    memcpy(fh + 2, &fsize, 4); uint32_t off = 54; memcpy(fh + 10, &off, 4);
    uint8_t ih[40] = { 0 }; uint32_t v;
    v = 40; memcpy(ih, &v, 4); v = W; memcpy(ih + 4, &v, 4); v = H; memcpy(ih + 8, &v, 4);
    uint16_t planes = 1, bpp = 24; memcpy(ih + 12, &planes, 2); memcpy(ih + 14, &bpp, 2);
    v = dataSize; memcpy(ih + 20, &v, 4);
    FILE* f = fopen(path, "wb"); fwrite(fh, 1, 14, f); fwrite(ih, 1, 40, f);
    std::vector<uint8_t> line(row, 0);
    for (int y = H - 1; y >= 0; --y) {
        for (int x = 0; x < W; ++x) { const uint8_t* p = &rgb[(y * W + x) * 3]; line[x*3]=p[2]; line[x*3+1]=p[1]; line[x*3+2]=p[0]; }
        fwrite(line.data(), 1, row, f);
    }
    fclose(f);
}

int main(int argc, char** argv) {
    if (argc < 9) { fprintf(stderr, "usage: palette_preview out W H cx cy scaleExp mxit palIdx [density] [SS]\n"); return 1; }
    const char* out = argv[1]; int W = atoi(argv[2]), H = atoi(argv[3]);
    const char* cx = argv[4]; const char* cy = argv[5]; const char* scaleArg = argv[6];
    int mxit = atoi(argv[7]); int palIdx = atoi(argv[8]);
    if (argc > 9) color_density = (float)atof(argv[9]);
    int SS = argc > 10 ? atoi(argv[10]) : 3;
    if (palIdx < 0 || palIdx >= (int)palettes.size()) palIdx = 0;
    int Ws = W * SS, Hs = H * SS;

    // Gradient-strip mode: dump the 1D palette (color_map) as a horizontal strip
    // to inspect the interpolation directly. Trigger with scaleArg == "strip".
    if (std::string(scaleArg) == "strip") {
        rebuild(palettes[palIdx].stops);
        std::vector<uint8_t> img((size_t)W * H * 3);
        for (int j = 0; j < W; ++j) {
            int idx = (int)((float)j / W * CP) % CP; if (idx < 0) idx += CP;
            for (int i = 0; i < H; ++i) {
                uint8_t* p = &img[((size_t)i * W + j) * 3];
                p[0] = (uint8_t)(color_map[0][idx] + 0.5f);
                p[1] = (uint8_t)(color_map[1][idx] + 0.5f);
                p[2] = (uint8_t)(color_map[2][idx] + 0.5f);
            }
        }
        writeBMP(out, img, W, H);
        printf("wrote strip %s  palette=%s\n", out, palettes[palIdx].name);
        return 0;
    }

    std::string sa = scaleArg, scale;
    if (sa.find_first_of("eE.") != std::string::npos) {
        size_t ep = sa.find_first_of("eE");
        std::string mant = (ep == std::string::npos) ? sa : sa.substr(0, ep);
        int e10 = (ep == std::string::npos) ? 0 : atoi(sa.substr(ep + 1).c_str());
        std::string sign; if (!mant.empty() && (mant[0]=='-'||mant[0]=='+')) { sign = mant[0]=='-'?"-":""; mant.erase(0,1); }
        size_t dp = mant.find('.'); int frac = (dp==std::string::npos)?0:(int)(mant.size()-dp-1);
        if (dp != std::string::npos) mant.erase(dp,1);
        int zeros = e10 - frac;
        if (zeros >= 0) { scale = sign + mant; scale.append(zeros, '0'); }
        else {
            int point = (int)mant.size() + zeros;   // decimal point position from the left
            if (point <= 0) scale = sign + "0." + std::string(-point, '0') + mant;
            else scale = sign + mant.substr(0, point) + "." + mant.substr(point);
        }
    } else { int se = atoi(sa.c_str()); scale = "1"; for (int i=0;i<se;++i) scale += "0"; }
    int precision = (int)(scale.size() * log(10) / log(2)) + 64;
    mpf_set_default_prec(precision);

    rebuild(palettes[palIdx].stops);
    float* iter = new float[Ws * Hs];
    for (int i = 0; i < Ws * Hs; ++i) iter[i] = EMPTYPIXEL;
    Mandel mandel(Ws, Hs, mxit, 1, iter);
    mandel.setPrecision(precision);
    mpf_t mcx, mcy, msc;
    mpf_init_set_str(mcx, cx, 10); mpf_init_set_str(mcy, cy, 10); mpf_init_set_str(msc, scale.c_str(), 10);
    mandel.Compute(mcx, mcy, msc, mxit, 0);   // non-EDE

    std::vector<uint8_t> img(W * H * 3);
    for (int i = 0; i < H; ++i) for (int j = 0; j < W; ++j) {
        float lr=0,lg=0,lb=0;
        for (int a=0;a<SS;++a) for (int b=0;b<SS;++b) { float r,g,bb; getColor(iter[(i*SS+a)*Ws+(j*SS+b)], r,g,bb); lr+=toLin(r);lg+=toLin(g);lb+=toLin(bb); }
        int n=SS*SS; uint8_t* p=&img[(i*W+j)*3];
        p[0]=(uint8_t)(toGam(lr/n)+0.5f); p[1]=(uint8_t)(toGam(lg/n)+0.5f); p[2]=(uint8_t)(toGam(lb/n)+0.5f);
    }
    writeBMP(out, img, W, H);
    printf("wrote %s  palette=%s\n", out, palettes[palIdx].name);
    return 0;
}
