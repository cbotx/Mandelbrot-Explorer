#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#ifndef _WIN32_WINNT
#define _WIN32_WINNT 0x0A00
#endif
#ifndef WINVER
#define WINVER 0x0A00
#endif
#include <windows.h>
#include <windowsx.h>
#include <commdlg.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

#include "Image.h"
#include "color.h"
#include "interpolate.h"
#include "mandel_navigator.h"

#pragma comment(lib, "comdlg32.lib")
#pragma comment(lib, "gdi32.lib")
#pragma comment(lib, "user32.lib")

float color_map[3][colP];
float color_density = 60.0f;

static float colorFunction(float it, int method) {
    if (method & ColoringMethod::EXTERIOR_DIST_EST)
        return tanhf(it / color_density * 5.0f);
    float l = logf(it + 2.0f);
    return powf(l, l * l / color_density);
}

void getColor(float iteration, float& r, float& g, float& b, int method) {
    if (iteration < 0) { r = g = b = 0; return; }
    int x = (int)(colorFunction(iteration, method) * colP) % colP;
    if (x < 0) x += colP;
    r = color_map[0][x]; g = color_map[1][x]; b = color_map[2][x];
}

void getColor(float iteration, uint8_t& r, uint8_t& g, uint8_t& b, int method) {
    float fr, fg, fb;
    getColor(iteration, fr, fg, fb, method);
    r = (uint8_t)fr; g = (uint8_t)fg; b = (uint8_t)fb;
}

// ---- analytic anti-aliasing (gradient-based palette prefiltering) ----------
// The GUI's supersampling averages colours in squared (RMS / gamma-2 linear)
// space, so the analytic filter integrates palette SQUARES and returns the RMS
// to stay consistent with the supersampled path.
static double g_palInt[3][colP + 1];   // prefix sum of color_map^2, [colP] = full sum
static double g_palMean[3];            // RMS of the whole palette (full-cycle limit)

void prepareColorFilter() {
    for (int c = 0; c < 3; ++c) {
        double s = 0; g_palInt[c][0] = 0;
        for (int i = 0; i < colP; ++i) {
            double v = color_map[c][i];
            s += v * v; g_palInt[c][i + 1] = s;
        }
        g_palMean[c] = std::sqrt(s / colP);
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

void getColorAA(float v, float vL, float vR, float vU, float vD,
                uint8_t& r, uint8_t& g, uint8_t& b, int method) {
    if (v < 0) { r = g = b = 0; return; }    // interior
    auto ext = [](float x) { return x != EMPTYPIXEL && x >= 0.0f; };
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
    double width = (double)gradf * colP;                 // palette entries spanned
    double centerU = (double)f * colP;

    if (width < 1.0) {                                   // sub-cycle: point sample
        int x = (int)centerU % colP; if (x < 0) x += colP;
        r = (uint8_t)color_map[0][x]; g = (uint8_t)color_map[1][x]; b = (uint8_t)color_map[2][x];
        return;
    }
    if (width >= colP) {                                 // >= a full cycle: palette RMS
        r = (uint8_t)g_palMean[0]; g = (uint8_t)g_palMean[1]; b = (uint8_t)g_palMean[2];
        return;
    }
    double a = centerU - 0.5 * width, e = centerU + 0.5 * width;
    auto rms = [&](int c) {
        double m = (palPrefix(c, e) - palPrefix(c, a)) / width;   // mean of squares
        if (m < 0) m = 0;                                          // guard fp cancellation
        double s = std::sqrt(m);
        return s > 255.0 ? 255.0f : (float)s;
    };
    r = (uint8_t)rms(0); g = (uint8_t)rms(1); b = (uint8_t)rms(2);
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

namespace {

constexpr int RENDER_W = 900;
constexpr int RENDER_H = 600;
constexpr int PANEL_W = 330;
constexpr int STATUS_H = 30;
constexpr UINT TIMER_ID = 1;

// --- theme -----------------------------------------------------------------
const COLORREF CLR_BG        = RGB(18, 20, 26);
const COLORREF CLR_PANEL     = RGB(26, 29, 37);
const COLORREF CLR_CARD      = RGB(37, 41, 52);
const COLORREF CLR_CARD_HOV  = RGB(48, 53, 66);
const COLORREF CLR_TRACK     = RGB(52, 57, 70);
const COLORREF CLR_ACCENT    = RGB(94, 148, 245);
const COLORREF CLR_ACCENT_HI = RGB(126, 172, 250);
const COLORREF CLR_GREEN     = RGB(112, 200, 158);
const COLORREF CLR_TEXT      = RGB(224, 229, 238);
const COLORREF CLR_TEXT_DIM  = RGB(140, 149, 166);
const COLORREF CLR_BORDER    = RGB(58, 64, 80);

// hit-test ids for self-drawn widgets
enum Hit {
    H_NONE, H_VIEW, H_GRADIENT,
    H_RESET, H_RENDER, H_SAVE, H_COPY,
    H_MAXFIELD, H_MAXTRACK, H_DENSTRACK,
    H_SS, H_EDE, H_PRESET_SNOWY, H_PRESET_SUNRISE, H_COLOR
};

struct Stop { float pos, r, g, b; };

struct PaletteEditor {
    std::vector<Stop> stops;
    int selected = 0;
    bool dragging = false;

    void sunrise() {
        stops = { {0.0f,0,70,100}, {0.16f,32,107,203}, {0.42f,237,255,255},
                  {0.6425f,255,170,0}, {0.8575f,0,2,0} };
        selected = 0; rebuild();
    }
    void snowy() {
        stops = { {0.0f,255,255,255}, {0.12f,32,107,203}, {0.34f,214,223,255} };
        selected = 0; rebuild();
    }
    void rebuild() {
        if (stops.empty()) return;
        int n = (int)stops.size();
        std::vector<float> xs(n + 2), ys(n + 2), out(colP);
        std::vector<std::array<float, 3>> lab(n);
        for (int i = 0; i < n; ++i)
            rgb2lab(stops[i].r, stops[i].g, stops[i].b, lab[i][0], lab[i][1], lab[i][2]);
        xs[0] = stops.back().pos - 1.0f;
        xs[n + 1] = stops.front().pos + 1.0f;
        for (int i = 0; i < n; ++i) xs[i + 1] = stops[i].pos;
        for (int c = 0; c < 3; ++c) {
            ys[0] = lab.back()[c]; ys[n + 1] = lab.front()[c];
            for (int i = 0; i < n; ++i) ys[i + 1] = lab[i][c];
            mono_cubic_interpolate(xs.data(), ys.data(), n + 2, out.data(), colP);
            for (int i = 0; i < colP; ++i) color_map[c][i] = out[i];
        }
        for (int i = 0; i < colP; ++i) {
            float r, g, b;
            lab2rgb(color_map[0][i], color_map[1][i], color_map[2][i], r, g, b);
            color_map[0][i] = r; color_map[1][i] = g; color_map[2][i] = b;
        }
    }
    std::array<uint8_t, 3> sample(float p) const {
        int i = std::clamp((int)(p * (colP - 1)), 0, colP - 1);
        return { (uint8_t)color_map[0][i], (uint8_t)color_map[1][i], (uint8_t)color_map[2][i] };
    }
};

static std::wstring widen(const std::string& s) {
    if (s.empty()) return {};
    int n = MultiByteToWideChar(CP_UTF8, 0, s.data(), (int)s.size(), nullptr, 0);
    std::wstring w(n, 0);
    MultiByteToWideChar(CP_UTF8, 0, s.data(), (int)s.size(), w.data(), n);
    return w;
}

static bool inRect(const RECT& r, int x, int y) {
    return x >= r.left && x < r.right && y >= r.top && y < r.bottom;
}

// --- GDI drawing helpers ---------------------------------------------------
static void fillRound(HDC dc, RECT r, COLORREF fill, COLORREF border, int rad) {
    HBRUSH b = CreateSolidBrush(fill);
    HPEN p = border == CLR_BG ? (HPEN)GetStockObject(NULL_PEN) : CreatePen(PS_SOLID, 1, border);
    HGDIOBJ ob = SelectObject(dc, b), op = SelectObject(dc, p);
    RoundRect(dc, r.left, r.top, r.right, r.bottom, rad, rad);
    SelectObject(dc, ob); SelectObject(dc, op);
    DeleteObject(b);
    if (border != CLR_BG) DeleteObject(p);
}
static void fillRect(HDC dc, RECT r, COLORREF c) {
    HBRUSH b = CreateSolidBrush(c); FillRect(dc, &r, b); DeleteObject(b);
}
static void drawText(HDC dc, RECT r, const std::wstring& s, COLORREF c, HFONT f, UINT fmt) {
    HGDIOBJ of = SelectObject(dc, f);
    SetTextColor(dc, c); SetBkMode(dc, TRANSPARENT);
    DrawTextW(dc, s.c_str(), -1, &r, fmt);
    SelectObject(dc, of);
}

class App {
public:
    HWND hwnd = nullptr;
    HFONT fUi = nullptr, fBold = nullptr, fSmall = nullptr, fMono = nullptr;

    std::unique_ptr<MandelNavigator> nav;
    std::vector<uint8_t> bitmap = std::vector<uint8_t>((size_t)RENDER_W * RENDER_H * 3, 0);
    std::vector<uint8_t> display = std::vector<uint8_t>((size_t)RENDER_W * RENDER_H * 3, 0); // BGR
    PaletteEditor palette;

    // widget rects (computed in layout())
    RECT rcReset{}, rcRender{}, rcSave{}, rcCopy{};
    RECT rcLocation{}, rcMaxField{}, rcMaxTrack{}, rcDensTrack{};
    RECT rcSS{}, rcEDE{}, rcSnowy{}, rcSunrise{}, rcColor{}, rcGradient{};

    // state
    int maxIter = 500000;
    bool autoMax = true;   // zoom-adaptive iteration cap (disabled on manual edit)
    bool ssOn = false, edeOn = true;
    int presetIdx = 0; // 0 snowy, 1 sunrise
    bool navDragging = false, wasComputing = false;
    Hit hover = H_NONE, pressed = H_NONE;
    bool maxEditing = false; std::wstring maxBuf; int caretTick = 0;
    int liveFrames = 0;   // frames still needing per-tick repaint (animation/compute)
    int dpi = 96;         // display DPI; all metrics scale by dpi/96
    std::chrono::steady_clock::time_point renderStart;
    double lastRenderMs = 0;

    int S(int v) const { return MulDiv(v, dpi, 96); }   // scale a design px to device px
    void keepLive(int frames = 60) { liveFrames = std::max(liveFrames, frames); }

    // ---- geometry ----
    RECT viewRect() const {
        RECT rc; GetClientRect(hwnd, &rc);
        int aw = std::max(1, (int)rc.right - S(PANEL_W));
        int ah = std::max(1, (int)rc.bottom - S(STATUS_H));
        double s = std::min((double)aw / RENDER_W, (double)ah / RENDER_H);
        int w = std::max(1, (int)(RENDER_W * s)), h = std::max(1, (int)(RENDER_H * s));
        return { (aw - w) / 2, (ah - h) / 2, (aw - w) / 2 + w, (ah - h) / 2 + h };
    }

    void layout() {
        RECT rc; GetClientRect(hwnd, &rc);
        int px = rc.right - S(PANEL_W) + S(18), w = S(PANEL_W) - S(36), y = S(18);
        int g = S(8), bh = S(30);
        int bw = (w - 3 * g) / 4;
        rcReset   = { px,              y, px + bw,          y + bh };
        rcRender  = { px + bw + g,     y, px + 2*bw + g,    y + bh };
        rcSave    = { px + 2*(bw+g),   y, px + 3*bw + 2*g,  y + bh };
        rcCopy    = { px + 3*(bw+g),   y, px + 4*bw + 3*g,  y + bh };
        y += bh + S(14);
        rcLocation = { px, y, px + w, y + S(84) }; y += S(84) + S(20);
        rcMaxField = { px + w - S(96), y - S(2), px + w, y + S(24) };
        rcMaxTrack = { px, y + S(32), px + w, y + S(46) }; y += S(60);
        rcDensTrack = { px, y + S(26), px + w, y + S(40) }; y += S(54);
        rcSS  = { px, y, px + w, y + bh }; y += S(38);
        rcEDE = { px, y, px + w, y + bh }; y += S(46);
        int half = (w - g) / 2;
        rcSnowy   = { px, y, px + half, y + bh };
        rcSunrise = { px + half + g, y, px + w, y + bh }; y += S(44);
        rcGradient = { px, y + S(22), px + w, y + S(62) }; y += S(84);
        rcColor = { px, y, px + w, y + bh };
    }

    void createFonts() {
        if (fUi) { DeleteObject(fUi); DeleteObject(fBold); DeleteObject(fSmall); DeleteObject(fMono); }
        fUi   = CreateFontW(-S(15),0,0,0,FW_NORMAL,0,0,0,DEFAULT_CHARSET,0,0,CLEARTYPE_QUALITY,0,L"Segoe UI");
        fBold = CreateFontW(-S(15),0,0,0,FW_SEMIBOLD,0,0,0,DEFAULT_CHARSET,0,0,CLEARTYPE_QUALITY,0,L"Segoe UI");
        fSmall= CreateFontW(-S(12),0,0,0,FW_NORMAL,0,0,0,DEFAULT_CHARSET,0,0,CLEARTYPE_QUALITY,0,L"Segoe UI");
        fMono = CreateFontW(-S(12),0,0,0,FW_NORMAL,0,0,0,DEFAULT_CHARSET,0,0,CLEARTYPE_QUALITY,0,L"Consolas");
    }

    // ---- view compositing ----
    void buildDisplay() {
        auto b = nav->GetBoundary();
        double lx = b[0][0], ly = b[0][1], rx = b[2][0], ry = b[2][1];
        // The reference/bitmap is y-up (row 0 = smallest imaginary = bottom) while
        // the display DIB is top-down; flip vertically so the image is upright and
        // pan/zoom track the cursor consistently on the y axis.
        for (int y = 0; y < RENDER_H; ++y) {
            double fy = ly + (ry - ly) * ((RENDER_H - 1 - y) + 0.5) / RENDER_H;
            int sy = (int)(fy * RENDER_H);
            for (int x = 0; x < RENDER_W; ++x) {
                double fx = lx + (rx - lx) * (x + 0.5) / RENDER_W;
                int sx = (int)(fx * RENDER_W);
                size_t d = ((size_t)y * RENDER_W + x) * 3;
                if (sx < 0 || sx >= RENDER_W || sy < 0 || sy >= RENDER_H) {
                    display[d] = display[d+1] = display[d+2] = 0;
                } else {
                    size_t s = ((size_t)sy * RENDER_W + sx) * 3;
                    display[d] = bitmap[s+2]; display[d+1] = bitmap[s+1]; display[d+2] = bitmap[s];
                }
            }
        }
    }
    void commitView() {
        buildDisplay();
        // display is top-down/upright; bitmap is y-up, so flip rows back when baking
        // the warped preview as the new base (keeps the transition frame upright).
        for (int y = 0; y < RENDER_H; ++y) {
            const uint8_t* src = &display[(size_t)y * RENDER_W * 3];
            uint8_t* dst = &bitmap[(size_t)(RENDER_H - 1 - y) * RENDER_W * 3];
            for (int x = 0; x < RENDER_W; ++x) {
                dst[x*3] = src[x*3+2]; dst[x*3+1] = src[x*3+1]; dst[x*3+2] = src[x*3];
            }
        }
        nav->UpdateCoords();
        startRender();
    }
    void startRender() {
        if (autoMax) { maxIter = nav->SuggestMxit(); nav->SetMxit(maxIter); }
        renderStart = std::chrono::steady_clock::now();
        wasComputing = true;
        nav->StartCompute();
        keepLive();
        InvalidateRect(hwnd, nullptr, FALSE);
    }
    POINT mapToRender(int x, int y) const {
        RECT r = viewRect();
        int px = (int)((x - r.left) * (double)RENDER_W / std::max(1L, r.right - r.left));
        int py = (int)((y - r.top) * (double)RENDER_H / std::max(1L, r.bottom - r.top));
        return { std::clamp(px, 0, RENDER_W - 1), std::clamp(py, 0, RENDER_H - 1) };
    }

    // ---- value mapping ----
    double maxToT() const { return std::clamp((log10((double)maxIter) - 2.0) / 4.0, 0.0, 1.0); }
    void setMaxFromT(double t, bool render) {
        autoMax = false;   // slider is a manual override
        maxIter = (int)std::round(pow(10.0, 2.0 + 4.0 * std::clamp(t, 0.0, 1.0)));
        maxIter = std::clamp(maxIter, 100, 1000000);
        nav->SetMxit(maxIter);
        if (render) startRender();
    }
    void setMethodFlag(int flag, bool on) {
        int m = nav->GetCMethod();
        if (on) m |= flag; else m &= ~flag;
        nav->SetCMethod(m);
        startRender();
    }

    // ---- palette / color ----
    int nearestStop(int x, int tol) const {
        int best = -1, dist = tol + 1;
        for (int i = 0; i < (int)palette.stops.size(); ++i) {
            int sx = rcGradient.left + (int)(palette.stops[i].pos * (rcGradient.right - rcGradient.left));
            int d = abs(x - sx);
            if (d < dist) { dist = d; best = i; }
        }
        return best;
    }
    void gradientDown(int x) {
        int n = nearestStop(x, 9);
        if (n >= 0) palette.selected = n;
        else {
            float p = std::clamp((float)(x - rcGradient.left) / (rcGradient.right - rcGradient.left), 0.f, 1.f);
            auto c = palette.sample(p);
            Stop s{ p, (float)c[0], (float)c[1], (float)c[2] };
            auto it = std::lower_bound(palette.stops.begin(), palette.stops.end(), p,
                [](const Stop& a, float v){ return a.pos < v; });
            palette.selected = (int)(it - palette.stops.begin());
            palette.stops.insert(it, s);
        }
        palette.dragging = true; SetCapture(hwnd);
        InvalidateRect(hwnd, nullptr, FALSE);
    }
    void gradientMove(int x) {
        if (!palette.dragging || palette.stops.empty()) return;
        int i = palette.selected;
        float p = std::clamp((float)(x - rcGradient.left) / (rcGradient.right - rcGradient.left), 0.f, 1.f);
        float lo = i ? palette.stops[i-1].pos + 0.002f : 0.f;
        float hi = i + 1 < (int)palette.stops.size() ? palette.stops[i+1].pos - 0.002f : 1.f;
        palette.stops[i].pos = std::clamp(p, lo, hi);
        palette.rebuild(); nav->SetRedisplay();
        InvalidateRect(hwnd, nullptr, FALSE);
    }
    void gradientDelete(int x) {
        int i = nearestStop(x, 9);
        if (i >= 0 && palette.stops.size() > 2) {
            palette.stops.erase(palette.stops.begin() + i);
            palette.selected = std::min(i, (int)palette.stops.size() - 1);
            palette.rebuild(); nav->SetRedisplay();
            InvalidateRect(hwnd, nullptr, FALSE);
        }
    }
    void chooseSelectedColor() {
        if (palette.stops.empty()) return;
        Stop& s = palette.stops[palette.selected];
        static COLORREF custom[16]{};
        CHOOSECOLORW cc{ sizeof(cc) };
        cc.hwndOwner = hwnd;
        cc.rgbResult = RGB((BYTE)s.r, (BYTE)s.g, (BYTE)s.b);
        cc.lpCustColors = custom; cc.Flags = CC_FULLOPEN | CC_RGBINIT;
        if (ChooseColorW(&cc)) {
            s.r = (float)GetRValue(cc.rgbResult);
            s.g = (float)GetGValue(cc.rgbResult);
            s.b = (float)GetBValue(cc.rgbResult);
            palette.rebuild(); nav->SetRedisplay();
            InvalidateRect(hwnd, nullptr, FALSE);
        }
    }

    void copyLocation() {
        std::wstring t = widen(nav->GetLocationText());
        if (!OpenClipboard(hwnd)) return;
        EmptyClipboard();
        size_t bytes = (t.size() + 1) * sizeof(wchar_t);
        HGLOBAL h = GlobalAlloc(GMEM_MOVEABLE, bytes);
        memcpy(GlobalLock(h), t.c_str(), bytes); GlobalUnlock(h);
        SetClipboardData(CF_UNICODETEXT, h); CloseClipboard();
    }
    static bool writeBMP(const wchar_t* path, const std::vector<uint8_t>& bgr) {
        FILE* f = nullptr;
        if (_wfopen_s(&f, path, L"wb") || !f) return false;
        int row = (RENDER_W * 3 + 3) & ~3, data = row * RENDER_H;
        uint8_t fh[14]{'B','M'}, ih[40]{};
        uint32_t size = 54 + data, off = 54, hs = 40, w = RENDER_W, h = RENDER_H;
        uint16_t planes = 1, bpp = 24;
        memcpy(fh+2,&size,4); memcpy(fh+10,&off,4);
        memcpy(ih,&hs,4); memcpy(ih+4,&w,4); memcpy(ih+8,&h,4);
        memcpy(ih+12,&planes,2); memcpy(ih+14,&bpp,2); memcpy(ih+20,&data,4);
        fwrite(fh,1,14,f); fwrite(ih,1,40,f);
        std::vector<uint8_t> line(row);
        for (int y = RENDER_H - 1; y >= 0; --y) {
            memcpy(line.data(), &bgr[(size_t)y * RENDER_W * 3], RENDER_W * 3);
            fwrite(line.data(), 1, row, f);
        }
        fclose(f); return true;
    }
    void saveImage() {
        wchar_t file[MAX_PATH] = L"mandelbrot.bmp";
        OPENFILENAMEW ofn{ sizeof(ofn) };
        ofn.hwndOwner = hwnd; ofn.lpstrFilter = L"Bitmap (*.bmp)\0*.bmp\0";
        ofn.lpstrFile = file; ofn.nMaxFile = MAX_PATH; ofn.lpstrDefExt = L"bmp";
        ofn.Flags = OFN_OVERWRITEPROMPT | OFN_PATHMUSTEXIST;
        if (GetSaveFileNameW(&ofn)) { buildDisplay(); writeBMP(file, display); }
    }

    void commitMaxEdit() {
        if (!maxEditing) return;
        maxEditing = false;
        if (maxBuf.empty()) { autoMax = true; }          // empty commit -> back to auto
        else { autoMax = false; maxIter = std::clamp(_wtoi(maxBuf.c_str()), 100, 1000000); nav->SetMxit(maxIter); }
        startRender();
    }

    // ---- widget drawing ----
    void drawButton(HDC dc, RECT r, const std::wstring& s, Hit id, bool accent) {
        bool hov = hover == id, prs = pressed == id;
        COLORREF bg = accent ? (prs ? RGB(70,120,220) : hov ? CLR_ACCENT_HI : CLR_ACCENT)
                             : (prs ? CLR_CARD_HOV : hov ? CLR_CARD_HOV : CLR_CARD);
        fillRound(dc, r, bg, accent ? bg : CLR_BORDER, S(8));
        drawText(dc, r, s, accent ? RGB(255,255,255) : CLR_TEXT, fUi,
                 DT_CENTER | DT_VCENTER | DT_SINGLELINE);
    }
    void drawSlider(HDC dc, RECT track, double t, Hit id) {
        int cy = (track.top + track.bottom) / 2, hb = S(3);
        RECT bar = { track.left, cy - hb, track.right, cy + hb };
        fillRound(dc, bar, CLR_TRACK, CLR_BG, S(6));
        int kx = track.left + (int)(t * (track.right - track.left));
        RECT fill = { track.left, cy - hb, kx, cy + hb };
        fillRound(dc, fill, CLR_ACCENT, CLR_BG, S(6));
        int kr = hover == id || pressed == id ? S(9) : S(7);
        HBRUSH kb = CreateSolidBrush(RGB(238,242,250));
        HPEN kp = CreatePen(PS_SOLID, S(2), CLR_ACCENT);
        HGDIOBJ ob = SelectObject(dc, kb), op = SelectObject(dc, kp);
        Ellipse(dc, kx - kr, cy - kr, kx + kr, cy + kr);
        SelectObject(dc, ob); SelectObject(dc, op); DeleteObject(kb); DeleteObject(kp);
    }
    void drawToggle(HDC dc, RECT r, const std::wstring& label, bool on, Hit id) {
        bool hov = hover == id;
        fillRound(dc, r, hov ? CLR_CARD_HOV : CLR_CARD, CLR_BORDER, S(8));
        RECT tr = r; tr.left += S(12);
        drawText(dc, tr, label, CLR_TEXT, fUi, DT_LEFT | DT_VCENTER | DT_SINGLELINE);
        int pw = S(40), ph = S(20), cy = (r.top + r.bottom) / 2, pad = S(12), inset = S(3);
        RECT pill = { r.right - pad - pw, cy - ph/2, r.right - pad, cy + ph/2 };
        fillRound(dc, pill, on ? CLR_ACCENT : CLR_TRACK, CLR_BG, ph);
        int kr = ph/2 - inset, kx = on ? pill.right - kr - inset : pill.left + kr + inset;
        HBRUSH kb = CreateSolidBrush(RGB(245,247,252));
        HGDIOBJ ob = SelectObject(dc, kb), op = SelectObject(dc, GetStockObject(NULL_PEN));
        Ellipse(dc, kx - kr, cy - kr, kx + kr, cy + kr);
        SelectObject(dc, ob); SelectObject(dc, op); DeleteObject(kb);
    }
    void drawSeg(HDC dc, RECT r, const std::wstring& s, bool sel, Hit id) {
        bool hov = hover == id;
        fillRound(dc, r, sel ? CLR_ACCENT : (hov ? CLR_CARD_HOV : CLR_CARD),
                  sel ? CLR_ACCENT : CLR_BORDER, S(8));
        drawText(dc, r, s, sel ? RGB(255,255,255) : CLR_TEXT, fUi,
                 DT_CENTER | DT_VCENTER | DT_SINGLELINE);
    }
    void label(HDC dc, int x, int y, const std::wstring& s) {
        RECT r = { x, y, x + S(300), y + S(18) };
        drawText(dc, r, s, CLR_TEXT_DIM, fSmall, DT_LEFT | DT_TOP | DT_SINGLELINE);
    }

    void paint() {
        PAINTSTRUCT ps; HDC wdc = BeginPaint(hwnd, &ps);
        RECT rc; GetClientRect(hwnd, &rc);
        HDC dc = CreateCompatibleDC(wdc);
        HBITMAP bmp = CreateCompatibleBitmap(wdc, rc.right, rc.bottom);
        HGDIOBJ obmp = SelectObject(dc, bmp);
        fillRect(dc, rc, CLR_BG);

        // fractal view
        RECT vr = viewRect();
        BITMAPINFO bi{}; bi.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
        bi.bmiHeader.biWidth = RENDER_W; bi.bmiHeader.biHeight = -RENDER_H;
        bi.bmiHeader.biPlanes = 1; bi.bmiHeader.biBitCount = 24; bi.bmiHeader.biCompression = BI_RGB;
        SetStretchBltMode(dc, HALFTONE);
        StretchDIBits(dc, vr.left, vr.top, vr.right-vr.left, vr.bottom-vr.top,
                      0,0, RENDER_W, RENDER_H, display.data(), &bi, DIB_RGB_COLORS, SRCCOPY);

        // panel
        RECT panel = { rc.right - S(PANEL_W), 0, rc.right, rc.bottom };
        fillRect(dc, panel, CLR_PANEL);
        HPEN sep = CreatePen(PS_SOLID, 1, CLR_BORDER);
        HGDIOBJ osep = SelectObject(dc, sep);
        MoveToEx(dc, panel.left, 0, nullptr); LineTo(dc, panel.left, rc.bottom);
        SelectObject(dc, osep); DeleteObject(sep);

        drawButton(dc, rcReset, L"Reset", H_RESET, false);
        drawButton(dc, rcRender, L"Render", H_RENDER, true);
        drawButton(dc, rcSave, L"Save", H_SAVE, false);
        drawButton(dc, rcCopy, L"Copy", H_COPY, false);

        // location card
        fillRound(dc, rcLocation, CLR_CARD, CLR_BORDER, S(8));
        RECT lt = rcLocation; lt.left += S(10); lt.top += S(6); lt.right -= S(10); lt.bottom -= S(6);
        drawText(dc, lt, widen(nav->GetLocationText()), CLR_TEXT, fMono,
                 DT_LEFT | DT_TOP | DT_WORDBREAK | DT_NOPREFIX | DT_EDITCONTROL);

        // max iterations
        label(dc, rcMaxField.left - S(140), rcMaxField.top + S(4), autoMax ? L"Max iterations (auto)" : L"Max iterations");
        {
            bool hov = hover == H_MAXFIELD || maxEditing;
            fillRound(dc, rcMaxField, hov ? CLR_CARD_HOV : CLR_CARD,
                      maxEditing ? CLR_ACCENT : CLR_BORDER, S(6));
            std::wstring t = maxEditing ? maxBuf : std::to_wstring(maxIter);
            if (maxEditing && (caretTick / 15) % 2 == 0) t += L"|";
            RECT tr = rcMaxField; tr.right -= S(10);
            drawText(dc, tr, t, CLR_TEXT, fUi, DT_RIGHT | DT_VCENTER | DT_SINGLELINE);
        }
        drawSlider(dc, rcMaxTrack, maxToT(), H_MAXTRACK);

        // density
        wchar_t db[48]; swprintf_s(db, L"Color density: %.0f", color_density);
        label(dc, rcDensTrack.left, rcDensTrack.top - S(22), db);
        drawSlider(dc, rcDensTrack, std::clamp((color_density - 10.0) / 190.0, 0.0, 1.0), H_DENSTRACK);

        drawToggle(dc, rcSS, L"5x supersampling", ssOn, H_SS);
        drawToggle(dc, rcEDE, L"Exterior distance estimation", edeOn, H_EDE);

        label(dc, rcSnowy.left, rcSnowy.top - S(20), L"Palette");
        drawSeg(dc, rcSnowy, L"Snowy", presetIdx == 0, H_PRESET_SNOWY);
        drawSeg(dc, rcSunrise, L"Sunrise", presetIdx == 1, H_PRESET_SUNRISE);

        // gradient bar
        label(dc, rcGradient.left, rcGradient.top - S(20), L"Gradient  (drag / right-click / dbl-click)");
        for (int x = rcGradient.left; x < rcGradient.right; ++x) {
            float p = (float)(x - rcGradient.left) / std::max(1L, rcGradient.right - rcGradient.left - 1);
            auto c = palette.sample(p);
            HPEN pen = CreatePen(PS_SOLID, 1, RGB(c[0], c[1], c[2]));
            HGDIOBJ op = SelectObject(dc, pen);
            MoveToEx(dc, x, rcGradient.top, nullptr); LineTo(dc, x, rcGradient.bottom);
            SelectObject(dc, op); DeleteObject(pen);
        }
        for (int i = 0; i < (int)palette.stops.size(); ++i) {
            int x = rcGradient.left + (int)(palette.stops[i].pos * (rcGradient.right - rcGradient.left));
            int tw = S(6), th = S(10);
            POINT pts[3]{ {x, rcGradient.bottom + 1}, {x - tw, rcGradient.bottom + 1 + th}, {x + tw, rcGradient.bottom + 1 + th} };
            auto& s = palette.stops[i];
            HBRUSH br = CreateSolidBrush(RGB((BYTE)s.r,(BYTE)s.g,(BYTE)s.b));
            HPEN pn = CreatePen(PS_SOLID, i == palette.selected ? S(2) : 1, RGB(255,255,255));
            HGDIOBJ ob2 = SelectObject(dc, br), op2 = SelectObject(dc, pn);
            Polygon(dc, pts, 3);
            SelectObject(dc, ob2); SelectObject(dc, op2); DeleteObject(br); DeleteObject(pn);
        }

        // selected color button (shows swatch)
        {
            bool hov = hover == H_COLOR;
            fillRound(dc, rcColor, hov ? CLR_CARD_HOV : CLR_CARD, CLR_BORDER, S(8));
            RECT sw = { rcColor.left + S(8), rcColor.top + S(6), rcColor.left + S(34), rcColor.bottom - S(6) };
            if (!palette.stops.empty()) {
                auto& s = palette.stops[palette.selected];
                fillRound(dc, sw, RGB((BYTE)s.r,(BYTE)s.g,(BYTE)s.b), CLR_BORDER, S(5));
            }
            RECT tr = rcColor; tr.left += S(44);
            drawText(dc, tr, L"Edit selected stop color", CLR_TEXT, fUi, DT_LEFT | DT_VCENTER | DT_SINGLELINE);
        }

        // status bar
        RECT status = { 0, rc.bottom - S(STATUS_H), rc.right - S(PANEL_W), rc.bottom };
        fillRect(dc, status, CLR_PANEL);
        std::wstring st = nav->IsComputing()
            ? L"  Rendering..."
            : L"  Ready   last render " + std::to_wstring((int)lastRenderMs) +
              L" ms    Drag pan   Wheel zoom   R reset   Space render   S save   C copy";
        RECT str = status; str.left += S(10);
        drawText(dc, str, st, CLR_TEXT_DIM, fSmall, DT_LEFT | DT_VCENTER | DT_SINGLELINE);

        BitBlt(wdc, 0, 0, rc.right, rc.bottom, dc, 0, 0, SRCCOPY);
        SelectObject(dc, obmp); DeleteObject(bmp); DeleteDC(dc);
        EndPaint(hwnd, &ps);
    }

    Hit hitTest(int x, int y) {
        if (inRect(rcReset,x,y)) return H_RESET;
        if (inRect(rcRender,x,y)) return H_RENDER;
        if (inRect(rcSave,x,y)) return H_SAVE;
        if (inRect(rcCopy,x,y)) return H_COPY;
        if (inRect(rcMaxField,x,y)) return H_MAXFIELD;
        RECT mt = rcMaxTrack; mt.top -= S(8); mt.bottom += S(8); if (inRect(mt,x,y)) return H_MAXTRACK;
        RECT dt = rcDensTrack; dt.top -= S(8); dt.bottom += S(8); if (inRect(dt,x,y)) return H_DENSTRACK;
        if (inRect(rcSS,x,y)) return H_SS;
        if (inRect(rcEDE,x,y)) return H_EDE;
        if (inRect(rcSnowy,x,y)) return H_PRESET_SNOWY;
        if (inRect(rcSunrise,x,y)) return H_PRESET_SUNRISE;
        if (inRect(rcColor,x,y)) return H_COLOR;
        RECT gr = rcGradient; gr.left -= S(8); gr.right += S(8); gr.top -= S(6); gr.bottom += S(14);
        if (inRect(gr,x,y)) return H_GRADIENT;
        RECT vr = viewRect(); if (inRect(vr,x,y)) return H_VIEW;
        return H_NONE;
    }

    void applyPreset(int idx) {
        presetIdx = idx;
        if (idx == 0) palette.snowy(); else palette.sunrise();
        nav->SetRedisplay();
        InvalidateRect(hwnd, nullptr, FALSE);
    }

    void timer() {
        bool computing = nav->IsComputing();
        bool active = computing || navDragging || palette.dragging || maxEditing ||
                      pressed == H_MAXTRACK || pressed == H_DENSTRACK || liveFrames > 0;
        if (!active) return;                 // idle: no repaint, no flicker, no CPU spin
        if (liveFrames > 0) --liveFrames;
        nav->Update();
        nav->UpdateBitmap(bitmap.data());
        if (wasComputing && !computing)
            lastRenderMs = std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - renderStart).count();
        wasComputing = computing;
        buildDisplay();
        if (maxEditing) ++caretTick;
        InvalidateRect(hwnd, nullptr, FALSE);
    }

    LRESULT message(UINT msg, WPARAM wp, LPARAM lp) {
        switch (msg) {
        case WM_CREATE: {
            dpi = (int)GetDpiForWindow(hwnd);
            if (dpi <= 0) dpi = 96;
            createFonts();
            palette.snowy();
            nav = std::make_unique<MandelNavigator>(RENDER_W, RENDER_H, 5, 1000000, 1.0, 220.0);
            nav->SetMxit(maxIter);
            nav->BindFixImageCallback(fixCallback);
            layout(); startRender();
            SetTimer(hwnd, TIMER_ID, 16, nullptr);
            return 0;
        }
        case WM_DPICHANGED: {
            dpi = HIWORD(wp);
            createFonts();
            RECT* sug = (RECT*)lp;
            SetWindowPos(hwnd, nullptr, sug->left, sug->top,
                         sug->right - sug->left, sug->bottom - sug->top,
                         SWP_NOZORDER | SWP_NOACTIVATE);
            layout(); InvalidateRect(hwnd, nullptr, FALSE);
            return 0;
        }
        case WM_SIZE: layout(); InvalidateRect(hwnd, nullptr, FALSE); return 0;
        case WM_ERASEBKGND: return 1;
        case WM_PAINT: paint(); return 0;
        case WM_TIMER: if (wp == TIMER_ID) timer(); return 0;
        case WM_MOUSEMOVE: {
            int x = GET_X_LPARAM(lp), y = GET_Y_LPARAM(lp);
            if (palette.dragging) { gradientMove(x); return 0; }
            if (pressed == H_MAXTRACK) { setMaxFromT((double)(x - rcMaxTrack.left)/(rcMaxTrack.right-rcMaxTrack.left), false); InvalidateRect(hwnd,nullptr,FALSE); return 0; }
            if (pressed == H_DENSTRACK) { color_density = (float)std::clamp(10.0 + 190.0*(x-rcDensTrack.left)/(rcDensTrack.right-rcDensTrack.left),10.0,200.0); nav->SetRedisplay(); InvalidateRect(hwnd,nullptr,FALSE); return 0; }
            if (navDragging) { POINT p = mapToRender(x,y); nav->Drag(p.x,p.y); return 0; }
            Hit h = hitTest(x,y);
            if (h != hover) { hover = h; InvalidateRect(hwnd, nullptr, FALSE); }
            return 0;
        }
        case WM_LBUTTONDOWN: {
            int x = GET_X_LPARAM(lp), y = GET_Y_LPARAM(lp);
            Hit h = hitTest(x,y);
            if (maxEditing && h != H_MAXFIELD) commitMaxEdit();
            pressed = h;
            if (h == H_GRADIENT) { gradientDown(x); return 0; }
            if (h == H_MAXTRACK) { SetCapture(hwnd); setMaxFromT((double)(x-rcMaxTrack.left)/(rcMaxTrack.right-rcMaxTrack.left), false); InvalidateRect(hwnd,nullptr,FALSE); return 0; }
            if (h == H_DENSTRACK) { SetCapture(hwnd); color_density=(float)std::clamp(10.0+190.0*(x-rcDensTrack.left)/(rcDensTrack.right-rcDensTrack.left),10.0,200.0); nav->SetRedisplay(); InvalidateRect(hwnd,nullptr,FALSE); return 0; }
            if (h == H_MAXFIELD) { maxEditing = true; maxBuf.clear(); caretTick = 0; InvalidateRect(hwnd,nullptr,FALSE); return 0; }
            if (h == H_VIEW) { POINT p = mapToRender(x,y); navDragging = true; SetCapture(hwnd); nav->DragStart(p.x,p.y); }
            InvalidateRect(hwnd, nullptr, FALSE);
            return 0;
        }
        case WM_LBUTTONUP: {
            int x = GET_X_LPARAM(lp), y = GET_Y_LPARAM(lp);
            Hit h = hitTest(x,y);
            if (palette.dragging) { palette.dragging = false; ReleaseCapture(); }
            else if (navDragging) { navDragging = false; nav->DragEnd(); ReleaseCapture(); }
            else if (pressed == H_MAXTRACK || pressed == H_DENSTRACK) {
                ReleaseCapture();
                if (pressed == H_MAXTRACK) startRender(); else startRender();
            } else if (h == pressed && h != H_NONE) {
                switch (h) {
                case H_RESET: nav->Reset(); renderStart = std::chrono::steady_clock::now(); wasComputing = true; keepLive(); break;
                case H_RENDER: startRender(); break;
                case H_SAVE: saveImage(); break;
                case H_COPY: copyLocation(); break;
                case H_SS: ssOn = !ssOn; setMethodFlag(ColoringMethod::SUPER_SAMPLING, ssOn); break;
                case H_EDE: edeOn = !edeOn; applyPreset(edeOn ? 0 : 1); setMethodFlag(ColoringMethod::EXTERIOR_DIST_EST, edeOn); break;
                case H_PRESET_SNOWY: applyPreset(0); break;
                case H_PRESET_SUNRISE: applyPreset(1); break;
                case H_COLOR: chooseSelectedColor(); break;
                default: break;
                }
            }
            pressed = H_NONE;
            InvalidateRect(hwnd, nullptr, FALSE);
            return 0;
        }
        case WM_LBUTTONDBLCLK: {
            int x = GET_X_LPARAM(lp), y = GET_Y_LPARAM(lp);
            Hit h = hitTest(x,y);
            if (h == H_GRADIENT) chooseSelectedColor();
            else if (h == H_VIEW) { POINT p = mapToRender(x,y); nav->ZoomIn(p.x,p.y); keepLive(); }
            return 0;
        }
        case WM_RBUTTONDOWN: {
            int x = GET_X_LPARAM(lp), y = GET_Y_LPARAM(lp);
            Hit h = hitTest(x,y);
            if (h == H_GRADIENT) gradientDelete(x);
            else if (h == H_VIEW) { POINT p = mapToRender(x,y); nav->ZoomOut(p.x,p.y); keepLive(); }
            return 0;
        }
        case WM_MOUSEWHEEL: {
            POINT q{ GET_X_LPARAM(lp), GET_Y_LPARAM(lp) }; ScreenToClient(hwnd, &q);
            if (inRect(viewRect(), q.x, q.y)) {
                POINT p = mapToRender(q.x, q.y);
                if (GET_WHEEL_DELTA_WPARAM(wp) > 0) nav->ZoomIn(p.x,p.y); else nav->ZoomOut(p.x,p.y);
                keepLive();
            }
            return 0;
        }
        case WM_CHAR:
            if (maxEditing) {
                wchar_t c = (wchar_t)wp;
                if (c >= '0' && c <= '9') { if (maxBuf.size() < 7) maxBuf += c; }
                else if (c == 8 && !maxBuf.empty()) maxBuf.pop_back();
                else if (c == '\r') commitMaxEdit();
                else if (c == 27) maxEditing = false;
                InvalidateRect(hwnd, nullptr, FALSE);
                return 0;
            }
            return 0;
        case WM_KEYDOWN:
            if (maxEditing) return 0;
            if (wp == 'R') { nav->Reset(); renderStart = std::chrono::steady_clock::now(); wasComputing = true; keepLive(); }
            else if (wp == VK_SPACE) startRender();
            else if (wp == 'S') saveImage();
            else if (wp == 'C') copyLocation();
            else if (wp == VK_ADD || wp == VK_OEM_PLUS) { nav->ZoomIn(RENDER_W/2, RENDER_H/2); keepLive(); }
            else if (wp == VK_SUBTRACT || wp == VK_OEM_MINUS) { nav->ZoomOut(RENDER_W/2, RENDER_H/2); keepLive(); }
            return 0;
        case WM_GETMINMAXINFO: {
            auto* m = (MINMAXINFO*)lp;
            m->ptMinTrackSize.x = S(900); m->ptMinTrackSize.y = S(640);
            return 0;
        }
        case WM_DESTROY:
            KillTimer(hwnd, TIMER_ID);
            nav.reset();
            DeleteObject(fUi); DeleteObject(fBold); DeleteObject(fSmall); DeleteObject(fMono);
            PostQuitMessage(0);
            return 0;
        }
        return DefWindowProcW(hwnd, msg, wp, lp);
    }

    static void fixCallback();
};

App* g_app = nullptr;
void App::fixCallback() { if (g_app) g_app->commitView(); }

LRESULT CALLBACK windowProc(HWND hwnd, UINT msg, WPARAM wp, LPARAM lp) {
    App* app = (App*)GetWindowLongPtrW(hwnd, GWLP_USERDATA);
    if (msg == WM_NCCREATE) {
        app = (App*)((CREATESTRUCTW*)lp)->lpCreateParams;
        app->hwnd = hwnd;
        SetWindowLongPtrW(hwnd, GWLP_USERDATA, (LONG_PTR)app);
    }
    return app ? app->message(msg, wp, lp) : DefWindowProcW(hwnd, msg, wp, lp);
}

} // namespace

int WINAPI wWinMain(HINSTANCE instance, HINSTANCE, PWSTR, int show) {
    // Per-Monitor-V2 DPI awareness (crisp text + WM_DPICHANGED); fall back on
    // older Windows to system DPI awareness.
    if (!SetProcessDpiAwarenessContext(DPI_AWARENESS_CONTEXT_PER_MONITOR_AWARE_V2))
        SetProcessDPIAware();
    WNDCLASSEXW wc{ sizeof(wc) };
    wc.style = CS_HREDRAW | CS_VREDRAW | CS_DBLCLKS | CS_OWNDC;
    wc.lpfnWndProc = windowProc;
    wc.hInstance = instance;
    wc.hCursor = LoadCursor(nullptr, IDC_ARROW);
    wc.hIcon = LoadIcon(nullptr, IDI_APPLICATION);
    wc.hbrBackground = nullptr;
    wc.lpszClassName = L"MandelbrotExplorerNative";
    RegisterClassExW(&wc);

    App app; g_app = &app;
    // Size the initial window for the system DPI, but clamp to the monitor work
    // area so it never exceeds the screen (which would clip the side panel).
    UINT sysdpi = GetDpiForSystem();
    int winW = MulDiv(1260, sysdpi, 96), winH = MulDiv(800, sysdpi, 96);
    RECT wa{}; SystemParametersInfoW(SPI_GETWORKAREA, 0, &wa, 0);
    winW = std::min(winW, (int)(wa.right - wa.left));
    winH = std::min(winH, (int)(wa.bottom - wa.top));
    int wx = wa.left + ((wa.right - wa.left) - winW) / 2;
    int wy = wa.top + ((wa.bottom - wa.top) - winH) / 2;
    HWND hwnd = CreateWindowExW(0, wc.lpszClassName, L"Mandelbrot Explorer",
        WS_OVERLAPPEDWINDOW | WS_VISIBLE, wx, wy, winW, winH,
        nullptr, nullptr, instance, &app);
    if (!hwnd) return 1;
    ShowWindow(hwnd, show); UpdateWindow(hwnd);

    MSG msg;
    while (GetMessageW(&msg, nullptr, 0, 0) > 0) {
        TranslateMessage(&msg);
        DispatchMessageW(&msg);
    }
    g_app = nullptr;
    return (int)msg.wParam;
}
