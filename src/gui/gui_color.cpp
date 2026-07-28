// GUI colouring: palette mapping, sRGB linear-light averaging, analytic AA,
// and the relief / normal-map post-shading. Implements Image.h.
#include <algorithm>
#include <cmath>
#include <cstdint>

#include "Image.h"
#include "interpolate.h"
#include "mandel_perturbation.h"

float color_map[3][colP];
float color_density = 60.0f;
float color_phase = 0.0f;

// Relief (screen-space slope) lighting. Off by default; enabled by the "Relief"
// coloring mode. Light from the upper-left, moderate slope.
int   relief_on = 0;
int   normal_light_on = 0;
int   de_overlay_on = 0;
float relief_light_az = 2.3f;
float relief_light_el = 0.55f;
float relief_strength = 1.0f;

// INTERIOR_SENTINEL (-2.0f) is declared in Image.h. A far-field point can escape
// with a slightly NEGATIVE smooth count, so "interior" must be the exact sentinel.
static inline bool isInterior(float it) { return it == INTERIOR_SENTINEL; }

static float colorFunction(float it, int method) {
    if (method & ColoringMethod::STRIPE_AVERAGE)
        return it * (color_density / 20.0f);   // SAC value in [0,1] -> banded palette
    if (method & ColoringMethod::ORBIT_TRAP)
        return it;                              // orbit-trap value is already a palette coordinate
    if (method & ColoringMethod::EXTERIOR_DIST_EST)
        return tanhf(it * color_density / 3600.0f * 5.0f);
    if (it < 0.0f) it = 0.0f;                  // far-field fast escapes: small negative count
    float l = logf(it + 2.0f);
    return powf(l, l * l * color_density / 3600.0f);
}

void getColor(float iteration, float& r, float& g, float& b, int method) {
    if (isInterior(iteration)) { r = g = b = 0; return; }
    int x = (int)(colorFunction(iteration, method) * colP + color_phase) % colP;
    if (x < 0) x += colP;
    r = color_map[0][x]; g = color_map[1][x]; b = color_map[2][x];
}

void getColor(float iteration, uint8_t& r, uint8_t& g, uint8_t& b, int method) {
    float fr, fg, fb;
    getColor(iteration, fr, fg, fb, method);
    r = (uint8_t)fr; g = (uint8_t)fg; b = (uint8_t)fb;
}

// ---- gamma-correct (sRGB linear-light) colour averaging -------------------
// Averaging in linear light makes a downsampled / analytically-filtered pixel
// match the physical "view from far" (radiance averages linearly), unlike the
// previous RMS/gamma-2 approximation.
float g_srgb2lin[256];
static bool g_lutReady = false;
// Linear-radiance [0,1] -> sRGB [0,255] lookup, so the hot re-colour paths avoid
// a pow() per channel. 8192 entries + linear interpolation keeps the error below
// ~0.2 of an 8-bit level (steepest near 0, where lerp halves the residual).
static constexpr int ENC_N = 8192;
float g_enc255[ENC_N + 1];
static void initSrgbLut() {
    for (int i = 0; i < 256; ++i) {
        double c = i / 255.0;
        g_srgb2lin[i] = (float)(c <= 0.04045 ? c / 12.92 : std::pow((c + 0.055) / 1.055, 2.4));
    }
    for (int i = 0; i <= ENC_N; ++i) {
        double L = (double)i / ENC_N;
        g_enc255[i] = (float)((L <= 0.0031308 ? 12.92 * L : 1.055 * std::pow(L, 1.0 / 2.4) - 0.055) * 255.0);
    }
    g_lutReady = true;
}
// Fast linear[0,1] -> sRGB[0,255] via the LUT (lerp). Clamps out-of-range input.
float srgbEncode255(double L) {
    if (L <= 0.0) return 0.0f;
    if (L >= 1.0) return 255.0f;
    double f = L * ENC_N; int i = (int)f; double t = f - i;
    return g_enc255[i] * (float)(1.0 - t) + g_enc255[i + 1] * (float)t;
}
static inline double srgb2linF(double v255) {           // fractional sRGB 0..255 -> linear
    double c = v255 / 255.0; if (c < 0) c = 0; else if (c > 1) c = 1;
    return c <= 0.04045 ? c / 12.92 : std::pow((c + 0.055) / 1.055, 2.4);
}
float srgbEncode(double L) {                            // linear 0..1 -> sRGB 0..1
    if (L <= 0.0) return 0.0f;
    if (L >= 1.0) return 1.0f;
    return (float)(L <= 0.0031308 ? 12.92 * L : 1.055 * std::pow(L, 1.0 / 2.4) - 0.055);
}
static inline double enc255(double L) { return srgbEncode(L) * 255.0; }
void ensureSrgbLut() { if (!g_lutReady) initSrgbLut(); }

// Lambert slope (relief) post-shade shared by the live view and the PNG export.
// height: per-pixel smooth value, NaN for interior/empty (left unshaded).
void applyReliefTo(uint8_t* rgb, const float* height, int W, int H) {
    const double lx = std::cos((double)relief_light_el) * std::cos((double)relief_light_az);
    const double ly = std::cos((double)relief_light_el) * std::sin((double)relief_light_az);
    const double lz = std::sin((double)relief_light_el);
    const double str = relief_strength;
#pragma omp parallel for schedule(static)
    for (int i = 0; i < H; ++i)
        for (int j = 0; j < W; ++j) {
            float hc = height[(size_t)i * W + j];
            if (std::isnan(hc)) continue;                 // interior/empty: unshaded
            auto Hh = [&](int yy, int xx) -> float {
                float t = height[(size_t)yy * W + xx];
                return std::isnan(t) ? hc : t;            // clamp across set boundary
            };
            float hl = j > 0     ? Hh(i, j - 1) : hc;
            float hr = j < W - 1 ? Hh(i, j + 1) : hc;
            float hu = i > 0     ? Hh(i - 1, j) : hc;
            float hd = i < H - 1 ? Hh(i + 1, j) : hc;
            double nx = -((double)(hr - hl) * 0.5) * str;
            double ny = -((double)(hd - hu) * 0.5) * str;
            double nz = 1.0;
            double inv = 1.0 / std::sqrt(nx * nx + ny * ny + nz * nz);
            double d = (nx * lx + ny * ly + nz * lz) * inv;
            if (d < 0) d = 0;
            double sh = 0.25 + 0.85 * d;                  // ambient + diffuse
            uint8_t* q = rgb + ((size_t)i * W + j) * 3;
            q[0] = (uint8_t)(srgbEncode(g_srgb2lin[q[0]] * sh) * 255.0 + 0.5);
            q[1] = (uint8_t)(srgbEncode(g_srgb2lin[q[1]] * sh) * 255.0 + 0.5);
            q[2] = (uint8_t)(srgbEncode(g_srgb2lin[q[2]] * sh) * 255.0 + 0.5);
        }
}

// Analytic normal-map post-shade: the engine gives a per-pixel surface-normal
// angle (arg(z)-arg(dz/dc)); tilt the unit (cos,sin) by relief_strength and
// Lambert-shade against the light. Deep-zoom-correct (unlike screen-space relief).
void applyNormalLightTo(uint8_t* rgb, const float* angle, int W, int H) {
    const double lx = std::cos((double)relief_light_el) * std::cos((double)relief_light_az);
    const double ly = std::cos((double)relief_light_el) * std::sin((double)relief_light_az);
    const double lz = std::sin((double)relief_light_el);
    const double str = relief_strength;
#pragma omp parallel for schedule(static)
    for (int i = 0; i < H; ++i)
        for (int j = 0; j < W; ++j) {
            float a = angle[(size_t)i * W + j];
            if (std::isnan(a)) continue;                  // interior/empty: unshaded
            double nx = std::cos((double)a) * str, ny = std::sin((double)a) * str, nz = 1.0;
            double inv = 1.0 / std::sqrt(nx * nx + ny * ny + nz * nz);
            double d = (nx * lx + ny * ly + nz * lz) * inv;
            if (d < 0) d = 0;
            double sh = 0.3 + 0.8 * d;
            uint8_t* q = rgb + ((size_t)i * W + j) * 3;
            q[0] = (uint8_t)(srgbEncode(g_srgb2lin[q[0]] * sh) * 255.0 + 0.5);
            q[1] = (uint8_t)(srgbEncode(g_srgb2lin[q[1]] * sh) * 255.0 + 0.5);
            q[2] = (uint8_t)(srgbEncode(g_srgb2lin[q[2]] * sh) * 255.0 + 0.5);
        }
}

// Distance-estimate B&W overlay: multiply the smooth colour by a shade that goes
// dark near the set boundary (small DE) and to full colour far away, drawing the
// fine filament structure over the base. relief_strength tunes the falloff.
void applyDEOverlayTo(uint8_t* rgb, const float* de, int W, int H) {
    // The distance estimate is drawn as a PURE black&white layer on top of the smooth
    // gradient: the gradient supplies the colour of the smooth exterior, while near the
    // set boundary the delicate DE structure takes over as greyscale (the filament lace
    // in black&white, fading to black on the set). No hue comes from the DE itself.
    const double k = 1.5 + 3.0 * (double)relief_strength;   // strength -> falloff distance
#pragma omp parallel for schedule(static)
    for (int i = 0; i < H; ++i)
        for (int j = 0; j < W; ++j) {
            float dd = de[(size_t)i * W + j];
            if (std::isnan(dd)) continue;                 // interior/empty: leave as-is
            double bw = (double)dd / ((double)dd + k);    // 0 on the boundary -> 1 far away
            double gray = bw;                             // the B&W distance value
            // colour weight: full gradient colour far away, pure grey near the boundary
            double cw = bw * bw;
            uint8_t* q = rgb + ((size_t)i * W + j) * 3;
            double r = g_srgb2lin[q[0]] * cw + gray * (1.0 - cw);
            double g = g_srgb2lin[q[1]] * cw + gray * (1.0 - cw);
            double b = g_srgb2lin[q[2]] * cw + gray * (1.0 - cw);
            q[0] = (uint8_t)(srgbEncode(r) * 255.0 + 0.5);
            q[1] = (uint8_t)(srgbEncode(g) * 255.0 + 0.5);
            q[2] = (uint8_t)(srgbEncode(b) * 255.0 + 0.5);
        }
}
// Integrates the palette in LINEAR light so the analytic average is gamma-correct
// and consistent with the supersampled (SmoothColor) path.
static double g_palInt[3][colP + 1];   // prefix sum of linear(color_map), [colP] = full sum
static double g_palMean[3];            // sRGB(mean linear) of the whole palette (full-cycle limit)
float g_palLin[3][colP];               // per-entry linear radiance of the palette (SS re-colour)

void prepareColorFilter() {
    if (!g_lutReady) initSrgbLut();
    for (int c = 0; c < 3; ++c) {
        double s = 0; g_palInt[c][0] = 0;
        for (int i = 0; i < colP; ++i) {
            double lin = srgb2linF(color_map[c][i]);
            g_palLin[c][i] = (float)lin;
            s += lin; g_palInt[c][i + 1] = s;
        }
        g_palMean[c] = enc255(s / colP);
    }
}

// prefix integral of channel-c squares from 0..x (x in palette-index units, any real)
static double palPrefix(int c, double x) {
    double q = std::floor(x / colP);
    double r = x - q * colP;                 // r in [0, colP)
    int i = (int)r; if (i >= colP) i = colP - 1;
    double t = r - i;
    double part = g_palInt[c][i] * (1.0 - t) + g_palInt[c][i + 1] * t;
    return q * g_palInt[c][colP] + part;
}

// Phase-independent analysis: the palette-index centre (baseU) and the palette
// footprint width. Neither depends on color_phase, so both can be cached and
// re-shaded cheaply as the phase animates. baseU < 0 marks interior/empty.
void colorAnalyzeAA(float v, float vL, float vR, float vU, float vD, int method,
                    float& baseU, float& width) {
    if (isInterior(v)) { baseU = -1.0f; width = 0.0f; return; }   // interior
    // A neighbour counts for the gradient if it is a real exterior sample (not
    // empty, not the interior sentinel). Far-field escapes have a small negative
    // smooth count but are still exterior, so test against the sentinel, not < 0.
    auto ext = [](float x) { return x != EMPTYPIXEL && x != INTERIOR_SENTINEL; };
    float f = colorFunction(v, method);
    auto cf = [&](float nv) { return colorFunction(nv, method); };
    // screen-space gradient of the (monotonic) colour value f, from neighbours
    float fx = 0, fy = 0;
    if (ext(vL) && ext(vR)) fx = 0.5f * (cf(vR) - cf(vL));
    else if (ext(vR)) fx = cf(vR) - f;
    else if (ext(vL)) fx = f - cf(vL);
    if (ext(vU) && ext(vD)) fy = 0.5f * (cf(vD) - cf(vU));
    else if (ext(vD)) fy = cf(vD) - f;
    else if (ext(vU)) fy = f - cf(vU);
    float gradf = std::sqrt(fx * fx + fy * fy);          // colour cycles per pixel
    width = (float)((double)gradf * colP);               // palette entries spanned
    baseU = (float)((double)f * colP);
}

// Phase-dependent shading of an analysed (baseU, width) pixel. Cheap: no
// colorFunction / gradient, just a palette-integral average shifted by phase.
void colorShadeAA(float baseU, float width, float phase, uint8_t& r, uint8_t& g, uint8_t& b) {
    if (baseU < 0) { r = g = b = 0; return; }            // interior
    double centerU = (double)baseU + phase;
    if (width < 1.0) {                                   // sub-cycle: point sample
        int x = ((int)centerU) % colP; if (x < 0) x += colP;
        r = (uint8_t)color_map[0][x]; g = (uint8_t)color_map[1][x]; b = (uint8_t)color_map[2][x];
        return;
    }
    if (width >= colP) {                                 // >= a full cycle: palette mean
        r = (uint8_t)g_palMean[0]; g = (uint8_t)g_palMean[1]; b = (uint8_t)g_palMean[2];
        return;
    }
    double a = centerU - 0.5 * width, e = centerU + 0.5 * width;
    auto avg = [&](int c) {
        double m = (palPrefix(c, e) - palPrefix(c, a)) / width;   // mean linear radiance
        if (m < 0) m = 0;                                          // guard fp cancellation
        return srgbEncode255(m);                                   // LUT encode to sRGB [0,255]
    };
    r = (uint8_t)avg(0); g = (uint8_t)avg(1); b = (uint8_t)avg(2);
}

float colorBaseIndex(float iteration, int method) {
    if (isInterior(iteration)) return -1.0f;
    return (float)((double)colorFunction(iteration, method) * colP);
}

void getColorAA(float v, float vL, float vR, float vU, float vD,
                uint8_t& r, uint8_t& g, uint8_t& b, int method) {
    float baseU, width;
    colorAnalyzeAA(v, vL, vR, vU, vD, method, baseU, width);
    colorShadeAA(baseU, width, color_phase, r, g, b);
}


void rgbRotate(float& r, float& g, float& b, float rad) {
    float c = cosf(rad), s = sinf(rad), t = 1.0f / 3.0f, q = sqrtf(t);
    float m1 = c + t * (1 - c), m2 = t * (1 - c) - q * s, m3 = t * (1 - c) + q * s;
    float nr = r * m1 + g * m2 + b * m3;
    float ng = r * m3 + g * m1 + b * m2;
    float nb = r * m2 + g * m3 + b * m1;
    r = std::clamp(nr, 0.0f, 255.0f);
    g = std::clamp(ng, 0.0f, 255.0f);
    b = std::clamp(nb, 0.0f, 255.0f);
}

void colorMapInitialize() {
    static const float pos[] = { 0.0f, 0.16f, 0.42f, 0.6425f, 0.8575f, 1.0f };
    static const float col[3][6] = {
        {0,32,237,255,0,0}, {70,107,255,170,2,70}, {100,203,255,0,0,100}
    };
    for (int c = 0; c < 3; ++c) mono_cubic_interpolate(pos, col[c], 6, color_map[c], colP);
}
