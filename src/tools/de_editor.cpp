// Standalone DE-gradient editor. Renders a few DE-rich preset views ONCE (caching
// the smooth escape buffer + the per-pixel distance-estimate buffer), then lets you
// tune the DE-overlay gradient with live trackbars -- re-colouring is milliseconds,
// so the preview updates in real time as you drag. Use it to dial in the falloff,
// vignette darkness, edge-glow strength/position and glow colour, then read the
// values off the panel and tell me which to bake into the app.
//
// The gradient model (per exterior pixel; interior stays black):
//   bw    = de / (de + k)                         0 on boundary .. 1 far exterior
//   base  = paletteColour * bw^gamma              dark vignette (gamma = darkness)
//   edge  = bw^ea * (1-bw)^eb   (peak at edgePos) filament-edge mask
//   glowC = lerp(glowRGB, paletteColour, tint)    white .. palette-tinted glow
//   out   = base + glowStrength * edge * glowC

#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <commctrl.h>
#include <gmp.h>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "mandel_perturbation.h"
#include "interpolate.h"

#pragma comment(lib, "comctl32.lib")

// ------------------------------------------------------------------ palette ---
static const int colP = 2048;
static float color_map[3][colP];
static void colorMapInit() {
    static const int N = 6;
    static const float pos[N] = { 0.0f, 0.16f, 0.42f, 0.6425f, 0.8575f, 1.0f };
    static const float col[3][N] = { {0,32,237,255,0,0}, {70,107,255,170,2,70}, {100,203,255,0,0,100} };
    for (int i = 0; i < 3; ++i) mono_cubic_interpolate(pos, col[i], N, color_map[i], colP);
}
static inline float color_func(float it, float dens) {
    return powf(logf(it + 2), logf(it + 2) * logf(it + 2) * dens / 3600.0f);
}
// Base smooth colour (0..255) for a pixel's escape value; <0 => interior (black).
static void getColorBase(float it, float dens, float& r, float& g, float& b) {
    if (it < 0) { r = g = b = 0; return; }
    int x = (int)(color_func(it, dens) * colP) % colP; if (x < 0) x += colP;
    r = color_map[0][x]; g = color_map[1][x]; b = color_map[2][x];
}
static inline double toLin(double c) { return pow(c / 255.0, 2.2); }
static inline double toGam(double l) { return 255.0 * pow(l < 0 ? 0 : l, 1.0 / 2.2); }

// ------------------------------------------------------------------ presets ---
struct Preset { const char* name; const char* x; const char* y; const char* scale; int mxit; float dens; };
static Preset g_presets[] = {
    { "Night City (1.69e40)",
      "-1.768628917759850520844734198472848718821423994141176532908",
      "0.001395534274228826510747662517373603005419245032078944176",
      "1.691960e40", 500000, 83.4f },
    { "Seahorse spiral (1e9)",
      "-0.743643887037151", "0.131825904205330", "1e9", 40000, 60.0f },
    { "Minibrot antenna (1e5)",
      "-1.2568853", "0.3803475", "1e5", 20000, 60.0f },
};
static const int NPRESET = (int)(sizeof(g_presets) / sizeof(g_presets[0]));

static const int PW = 640, PH = 480, SS = 2;   // preview output size + supersampling

struct Cache { std::vector<float> iter, de; bool done = false; };
static Cache g_cache[NPRESET];
static int g_cur = 0;
static volatile bool g_rendering = false;

// Expand a scale string like "1.691960e40" or "1e9" to an integer digit string.
static std::string expandScale(const char* sa_) {
    std::string sa = sa_, scale;
    if (sa.find_first_of("eE.") != std::string::npos) {
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
        int e = atoi(sa.c_str()); scale = "1"; for (int i = 0; i < e; ++i) scale += "0";
    }
    return scale;
}

static void renderPreset(int p) {
    if (g_cache[p].done) return;
    g_rendering = true;
    const int Ws = PW * SS, Hs = PH * SS;
    Cache& c = g_cache[p];
    c.iter.assign((size_t)Ws * Hs, EMPTYPIXEL);
    c.de.assign((size_t)Ws * Hs, 0.0f);
    std::string scale = expandScale(g_presets[p].scale);
    int precision = (int)(scale.size() * log(10) / log(2)) + 64;
    mpf_set_default_prec(precision);
    mpf_t mcx, mcy, msc;
    mpf_init_set_str(mcx, g_presets[p].x, 10);
    mpf_init_set_str(mcy, g_presets[p].y, 10);
    mpf_init_set_str(msc, scale.c_str(), 10);
    Mandel mandel(Ws, Hs, g_presets[p].mxit, 1, c.iter.data());
    mandel.setPrecision(precision);
    mandel.setNormalBuffer(c.de.data());
    mandel.setDensity(g_presets[p].dens);
    mandel.Compute(mcx, mcy, msc, g_presets[p].mxit, ColoringMethod::DE_OVERLAY);
    mpf_clear(mcx); mpf_clear(mcy); mpf_clear(msc);
    c.done = true;
    g_rendering = false;
}

// ------------------------------------------------------------- gradient params ---
// The DE overlay as an RGBA gradient over the smooth base. de is proportional to
// the true distance from the set boundary (the light source), so the fade follows
// the optical inverse-square law, blended in linear light:
//   colour: BLACK -> glow over bw in [0, riseW], shaped by riseCurve
//   alpha:  1 within a core radius coreDe, then intensity ~ 1/(1+(r/scaleDe)^fadeExp)
//           fadeExp = 2 is the physical inverse-square law (1/r^2); higher = tighter
struct GParams {
    double k = 2.0;          // bw = de/(de+k)  (colour-ramp scale)
    double riseW = 0.08;     // black -> glow rise width (edge sharpness, in bw)
    double riseCurve = 1.0;  // rise shape: >1 slow-start=more black, <1 fast-start=more white
    double coreDe = 0.2;     // glow stays fully opaque within this DE radius (light core)
    double scaleDe = 1.5;    // optical falloff scale (DE distance where it half-dims)
    double fadeExp = 2.0;    // falloff exponent: 2 = inverse-square (1/r^2) optical law
    double glowR = 160, glowG = 160, glowB = 160;   // glow colour (grey by default)
} g_p;

// level = glow-colour amount (0 black .. 1 full glow); alpha = overlay opacity.
static void deOverlay(double de, double& level, double& alpha) {
    double bw = de / (de + g_p.k);
    level = (bw >= g_p.riseW) ? 1.0 : pow(bw / (g_p.riseW + 1e-9), g_p.riseCurve);
    double r = de - g_p.coreDe;                       // distance beyond the light core
    if (r <= 0.0) alpha = 1.0;
    else alpha = 1.0 / (1.0 + pow(r / (g_p.scaleDe + 1e-9), g_p.fadeExp));  // inverse-power
}

static std::vector<uint8_t> g_img;   // PW*PH*3, BGR top-down (DIB-ready)

static void recolor() {
    Cache& c = g_cache[g_cur];
    if (!c.done) return;
    const int Ws = PW * SS;
    if ((int)g_img.size() != PW * PH * 3) g_img.assign((size_t)PW * PH * 3, 0);
    const double gcR = toLin(g_p.glowR), gcG = toLin(g_p.glowG), gcB = toLin(g_p.glowB);
#pragma omp parallel for schedule(dynamic, 8)
    for (int y = 0; y < PH; ++y)
        for (int x = 0; x < PW; ++x) {
            double lr = 0, lg = 0, lb = 0;
            for (int a = 0; a < SS; ++a)
                for (int b = 0; b < SS; ++b) {
                    int yy = y * SS + a, xx = x * SS + b;
                    float it = c.iter[(size_t)yy * Ws + xx];
                    float br, bg, bb; getColorBase(it, g_presets[g_cur].dens, br, bg, bb);
                    double clinR = toLin(br), clinG = toLin(bg), clinB = toLin(bb);
                    if (it < 0) continue;                     // interior -> black
                    double de = c.de[(size_t)yy * Ws + xx];
                    double level, alpha; deOverlay(de, level, alpha);
                    double glR = gcR * level, glG = gcG * level, glB = gcB * level;
                    // composite the overlay (glow, opacity alpha) over the smooth base
                    lr += glR * alpha + clinR * (1.0 - alpha);
                    lg += glG * alpha + clinG * (1.0 - alpha);
                    lb += glB * alpha + clinB * (1.0 - alpha);
                }
            const double n = SS * SS;
            uint8_t* q = &g_img[((size_t)y * PW + x) * 3];
            q[0] = (uint8_t)(toGam(lb / n) + 0.5);   // BGR for DIB
            q[1] = (uint8_t)(toGam(lg / n) + 0.5);
            q[2] = (uint8_t)(toGam(lr / n) + 0.5);
        }
}

// ------------------------------------------------------------------- Win32 UI ---
enum { PANEL_W = 300, CLIENT_H = 574,
       STRIPX = 12, STRIPY = 496, STRIPW = 616, STRIPH = 46 };
struct Slider { std::wstring label; HWND track, lab; double lo, hi; double* p; };
static std::vector<Slider> g_sliders;
static HWND g_hwnd, g_status;
static const int ID_PRESET0 = 2000;

static HWND mkTrack(HWND parent, int x, int y, int id) {
    return CreateWindowExW(0, TRACKBAR_CLASSW, L"", WS_CHILD | WS_VISIBLE | TBS_HORZ | TBS_NOTICKS,
                           x, y, PANEL_W - 24, 22, parent, (HMENU)(INT_PTR)id, GetModuleHandleW(nullptr), nullptr);
}
static void addSlider(HWND parent, int& y, const char* label, double lo, double hi, double* p) {
    int id = 1000 + (int)g_sliders.size();
    HWND lab = CreateWindowExW(0, L"STATIC", L"", WS_CHILD | WS_VISIBLE, PW + 12, y, PANEL_W - 24, 16, parent, nullptr, GetModuleHandleW(nullptr), nullptr);
    HWND tr = mkTrack(parent, PW + 12, y + 16, id);
    SendMessageW(tr, TBM_SETRANGE, TRUE, MAKELPARAM(0, 1000));
    SendMessageW(tr, TBM_SETPOS, TRUE, (LPARAM)((*p - lo) / (hi - lo) * 1000.0));
    wchar_t wl[128]; MultiByteToWideChar(CP_UTF8, 0, label, -1, wl, 128);
    g_sliders.push_back({ wl, tr, lab, lo, hi, p });
    y += 44;
}
static void refreshLabels() {
    for (auto& s : g_sliders) {
        wchar_t buf[160]; swprintf(buf, 160, L"%s = %.3f", s.label.c_str(), *s.p);
        SetWindowTextW(s.lab, buf);
    }
    // status line: copy-pasteable summary of all params
    wchar_t st[256];
    swprintf(st, 256, L"k=%.3f riseW=%.3f curve=%.3f\ncoreDe=%.3f scaleDe=%.3f exp=%.2f\nrgb=%.0f,%.0f,%.0f",
             g_p.k, g_p.riseW, g_p.riseCurve, g_p.coreDe, g_p.scaleDe, g_p.fadeExp, g_p.glowR, g_p.glowG, g_p.glowB);
    SetWindowTextW(g_status, st);
}

static void blit(HDC dc) {
    BITMAPINFO bi{}; bi.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
    bi.bmiHeader.biWidth = PW; bi.bmiHeader.biHeight = -PH;   // top-down
    bi.bmiHeader.biPlanes = 1; bi.bmiHeader.biBitCount = 24; bi.bmiHeader.biCompression = BI_RGB;
    if ((int)g_img.size() == PW * PH * 3)
        SetDIBitsToDevice(dc, 0, 0, PW, PH, 0, 0, 0, PH, g_img.data(), &bi, DIB_RGB_COLORS);
}

// The gradient strip: the DE overlay (glow colour + alpha) composited over a
// checkerboard so transparency is visible. x maps bw 0 (boundary, left) -> 1
// (far exterior, right); shows the black core, glow, fade and transparent rest.
static void drawStrip(HDC dc) {
    std::vector<uint8_t> strip((size_t)STRIPW * STRIPH * 3);
    const double gcR = toLin(g_p.glowR), gcG = toLin(g_p.glowG), gcB = toLin(g_p.glowB);
    for (int x = 0; x < STRIPW; ++x) {
        double bw = (double)x / (STRIPW - 1);
        double de = g_p.k * bw / (1.0 - bw + 1e-9);
        double level, alpha; deOverlay(de, level, alpha);
        double oR = gcR * level, oG = gcG * level, oB = gcB * level;
        for (int y = 0; y < STRIPH; ++y) {
            double ck = (((x / 8) + (y / 8)) & 1) ? 0.75 : 0.35;   // checker (linear grey)
            double r = oR * alpha + ck * (1 - alpha);
            double g = oG * alpha + ck * (1 - alpha);
            double b = oB * alpha + ck * (1 - alpha);
            uint8_t* q = &strip[((size_t)y * STRIPW + x) * 3];
            q[0] = (uint8_t)(toGam(b) + 0.5); q[1] = (uint8_t)(toGam(g) + 0.5); q[2] = (uint8_t)(toGam(r) + 0.5);
        }
    }
    BITMAPINFO bi{}; bi.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
    bi.bmiHeader.biWidth = STRIPW; bi.bmiHeader.biHeight = -STRIPH;
    bi.bmiHeader.biPlanes = 1; bi.bmiHeader.biBitCount = 24; bi.bmiHeader.biCompression = BI_RGB;
    SetDIBitsToDevice(dc, STRIPX, STRIPY, STRIPW, STRIPH, 0, 0, 0, STRIPH, strip.data(), &bi, DIB_RGB_COLORS);
    // labels
    SetBkMode(dc, TRANSPARENT); SetTextColor(dc, RGB(210, 210, 210));
    TextOutW(dc, STRIPX, STRIPY + STRIPH + 2, L"boundary (de=0)", 15);
    TextOutW(dc, STRIPX + STRIPW - 96, STRIPY + STRIPH + 2, L"far / transparent", 17);
}

static LRESULT CALLBACK WndProc(HWND h, UINT m, WPARAM w, LPARAM l) {
    switch (m) {
    case WM_HSCROLL: {
        for (auto& s : g_sliders)
            if ((HWND)l == s.track) {
                int pos = (int)SendMessageW(s.track, TBM_GETPOS, 0, 0);
                *s.p = s.lo + (s.hi - s.lo) * pos / 1000.0;
            }
        recolor(); refreshLabels();
        InvalidateRect(h, nullptr, FALSE);
        return 0;
    }
    case WM_COMMAND:
        if (LOWORD(w) >= ID_PRESET0 && LOWORD(w) < ID_PRESET0 + NPRESET) {
            g_cur = LOWORD(w) - ID_PRESET0;
            SetWindowTextW(g_status, L"Rendering preset...");
            UpdateWindow(h);
            renderPreset(g_cur);
            recolor(); refreshLabels();
            InvalidateRect(h, nullptr, FALSE);
        }
        return 0;
    case WM_PAINT: {
        PAINTSTRUCT ps; HDC dc = BeginPaint(h, &ps);
        RECT rc; GetClientRect(h, &rc);
        HBRUSH bg = CreateSolidBrush(RGB(30, 30, 34));
        RECT panel{ PW, 0, rc.right, rc.bottom }; FillRect(dc, &panel, bg);
        RECT below{ 0, PH, PW, rc.bottom }; FillRect(dc, &below, bg);
        DeleteObject(bg);
        blit(dc);
        drawStrip(dc);
        EndPaint(h, &ps);
        return 0;
    }
    case WM_DESTROY: PostQuitMessage(0); return 0;
    }
    return DefWindowProcW(h, m, w, l);
}

int WINAPI wWinMain(HINSTANCE hInst, HINSTANCE, PWSTR, int) {
    SetProcessDPIAware();
    INITCOMMONCONTROLSEX ic{ sizeof(ic), ICC_BAR_CLASSES }; InitCommonControlsEx(&ic);
    colorMapInit();
    WNDCLASSW wc{}; wc.lpfnWndProc = WndProc; wc.hInstance = hInst; wc.lpszClassName = L"DEGradEditor";
    wc.hCursor = LoadCursorW(nullptr, IDC_ARROW); wc.hbrBackground = (HBRUSH)(COLOR_WINDOW + 1);
    RegisterClassW(&wc);
    RECT r{ 0, 0, PW + PANEL_W, CLIENT_H }; AdjustWindowRect(&r, WS_OVERLAPPEDWINDOW, FALSE);
    g_hwnd = CreateWindowW(L"DEGradEditor", L"DE Gradient Editor", WS_OVERLAPPEDWINDOW & ~WS_MAXIMIZEBOX & ~WS_THICKFRAME,
                           CW_USEDEFAULT, CW_USEDEFAULT, r.right - r.left, r.bottom - r.top, nullptr, nullptr, hInst, nullptr);
    int y = 10;
    for (int i = 0; i < NPRESET; ++i) {
        wchar_t wl[128]; MultiByteToWideChar(CP_UTF8, 0, g_presets[i].name, -1, wl, 128);
        CreateWindowW(L"BUTTON", wl, WS_CHILD | WS_VISIBLE | BS_PUSHBUTTON,
                      PW + 12, y, PANEL_W - 24, 26, g_hwnd, (HMENU)(INT_PTR)(ID_PRESET0 + i), hInst, nullptr);
        y += 30;
    }
    y += 8;
    addSlider(g_hwnd, y, "Boundary scale k (de->bw)", 0.1, 8.0, &g_p.k);
    addSlider(g_hwnd, y, "Glow rise width (edge)", 0.005, 0.3, &g_p.riseW);
    addSlider(g_hwnd, y, "Rise curve (>1 blacker, <1 whiter)", 0.2, 5.0, &g_p.riseCurve);
    addSlider(g_hwnd, y, "Light core radius (DE)", 0.0, 3.0, &g_p.coreDe);
    addSlider(g_hwnd, y, "Optical falloff scale (DE)", 0.1, 8.0, &g_p.scaleDe);
    addSlider(g_hwnd, y, "Falloff exponent (2 = inverse-square)", 0.5, 4.0, &g_p.fadeExp);
    addSlider(g_hwnd, y, "Glow colour R", 0.0, 255.0, &g_p.glowR);
    addSlider(g_hwnd, y, "Glow colour G", 0.0, 255.0, &g_p.glowG);
    addSlider(g_hwnd, y, "Glow colour B", 0.0, 255.0, &g_p.glowB);
    g_status = CreateWindowExW(0, L"STATIC", L"", WS_CHILD | WS_VISIBLE, PW + 12, y + 6, PANEL_W - 24, 56, g_hwnd, nullptr, hInst, nullptr);

    ShowWindow(g_hwnd, SW_SHOW);
    SetWindowTextW(g_status, L"Rendering preset...");
    UpdateWindow(g_hwnd);
    renderPreset(0); recolor(); refreshLabels();
    InvalidateRect(g_hwnd, nullptr, FALSE);

    MSG msg;
    while (GetMessageW(&msg, nullptr, 0, 0)) { TranslateMessage(&msg); DispatchMessageW(&msg); }
    return 0;
}
