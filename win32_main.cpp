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
#include <commctrl.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "Image.h"
#include "color.h"
#include "interpolate.h"
#include "mandel_navigator.h"
#include "mandel_perturbation.h"

#pragma comment(lib, "comdlg32.lib")
#pragma comment(lib, "gdi32.lib")
#pragma comment(lib, "user32.lib")

float color_map[3][colP];
float color_density = 60.0f;

static float colorFunction(float it, int method) {
    if (method & ColoringMethod::STRIPE_AVERAGE)
        return it * (color_density / 20.0f);   // SAC value in [0,1] -> banded palette
    if (method & ColoringMethod::EXTERIOR_DIST_EST)
        return tanhf(it * color_density / 3600.0f * 5.0f);
    float l = logf(it + 2.0f);
    return powf(l, l * l * color_density / 3600.0f);
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

// ---- gamma-correct (sRGB linear-light) colour averaging -------------------
// Averaging in linear light makes a downsampled / analytically-filtered pixel
// match the physical "view from far" (radiance averages linearly), unlike the
// previous RMS/gamma-2 approximation.
float g_srgb2lin[256];
static bool g_lutReady = false;
static void initSrgbLut() {
    for (int i = 0; i < 256; ++i) {
        double c = i / 255.0;
        g_srgb2lin[i] = (float)(c <= 0.04045 ? c / 12.92 : std::pow((c + 0.055) / 1.055, 2.4));
    }
    g_lutReady = true;
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

// ---- analytic anti-aliasing (gradient-based palette prefiltering) ----------
// Integrates the palette in LINEAR light so the analytic average is gamma-correct
// and consistent with the supersampled (SmoothColor) path.
static double g_palInt[3][colP + 1];   // prefix sum of linear(color_map), [colP] = full sum
static double g_palMean[3];            // sRGB(mean linear) of the whole palette (full-cycle limit)

void prepareColorFilter() {
    if (!g_lutReady) initSrgbLut();
    for (int c = 0; c < 3; ++c) {
        double s = 0; g_palInt[c][0] = 0;
        for (int i = 0; i < colP; ++i) {
            s += srgb2linF(color_map[c][i]); g_palInt[c][i + 1] = s;
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
    if (width >= colP) {                                 // >= a full cycle: palette mean
        r = (uint8_t)g_palMean[0]; g = (uint8_t)g_palMean[1]; b = (uint8_t)g_palMean[2];
        return;
    }
    double a = centerU - 0.5 * width, e = centerU + 0.5 * width;
    auto avg = [&](int c) {
        double m = (palPrefix(c, e) - palPrefix(c, a)) / width;   // mean linear radiance
        if (m < 0) m = 0;                                          // guard fp cancellation
        double s = enc255(m);                                      // encode back to sRGB
        return s > 255.0 ? 255.0f : (float)s;
    };
    r = (uint8_t)avg(0); g = (uint8_t)avg(1); b = (uint8_t)avg(2);
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

// ---- high-resolution export (strip renderer) -----------------------------
// Renders the current view at an arbitrary W x H with ss x ss supersampling,
// in horizontal strips so peak memory is ~one strip + the reference orbit (a
// few hundred MB) regardless of output size. Colours with the live palette /
// density / method and box-downsamples in linear light. Reports progress and
// honours a cancel flag; meant to run on a background thread.
struct ExportState {
    std::atomic<float> progress{ 0.0f };
    std::atomic<bool> cancel{ false };
    std::atomic<bool> done{ false };
    std::atomic<bool> ok{ false };
    std::atomic<Mandel*> cur{ nullptr };   // current strip's engine, for cancel
    std::mutex curMx;                       // guards cur while the UI calls SetHalt
};

// Match the on-screen coordinate rectangle: identical if the export aspect
// equals the view's (RENDER_W:RENDER_H), otherwise keep the centre and expand
// the short side outward. scale_out = scale_view / max(1, aspect / (W0/H0)).
static void exportScale(mpf_t scale_view, int W, int H, mpf_t scale_out) {
    double a = (double)W / (double)H;
    double factor = a / ((double)RENDER_W / (double)RENDER_H);
    if (factor < 1.0) factor = 1.0;
    mpf_set_prec(scale_out, mpf_get_prec(scale_view));
    if (factor <= 1.0) { mpf_set(scale_out, scale_view); return; }
    mpf_t f; mpf_init_set_d(f, factor);
    mpf_div(scale_out, scale_view, f);
    mpf_clear(f);
}

static void exportRender(mpf_t cx, mpf_t cy, mpf_t scale_view,
                         int W, int H, int ss, int mxit, int cmethod,
                         std::vector<uint8_t>& out, ExportState* st) {
    if (!g_lutReady) initSrgbLut();
    out.assign((size_t)W * H * 3, 0);
    unsigned long prec = mpf_get_prec(scale_view);
    mpf_t scale_e; mpf_init2(scale_e, prec);
    exportScale(scale_view, W, H, scale_e);

    const int Wss = W * ss, Hss = H * ss;
    // strip height (output rows), bounded so the strip buffer stays ~<=200 MB.
    long budget = 200L * 1024 * 1024;
    int sh = (int)std::max<long>(1, budget / ((long)Wss * ss * 5));
    sh = std::clamp(sh, 1, H);
    const int nstrips = (H + sh - 1) / sh;

    std::vector<float> ibuf((size_t)Wss * sh * ss, EMPTYPIXEL);
    Mandel mandel(Wss, sh * ss, mxit, 1, ibuf.data());
    mandel.setPrecision(prec);
    st->cur = &mandel;

    for (int s = 0; s < nstrips && !st->cancel; ++s) {
        int obase = s * sh;
        int oh = std::min(sh, H - obase);
        std::fill(ibuf.begin(), ibuf.end(), (float)EMPTYPIXEL);
        mandel.SetHalt(false);
        mandel.Compute(cx, cy, scale_e, mxit, cmethod, Hss, obase * ss);
        if (st->cancel) break;
#pragma omp parallel for schedule(dynamic, 4)
        for (int oi = 0; oi < oh; ++oi) {
            for (int oj = 0; oj < W; ++oj) {
                double lr = 0, lg = 0, lb = 0; int n = ss * ss;
                for (int a = 0; a < ss; ++a)
                    for (int b = 0; b < ss; ++b) {
                        float v = ibuf[(size_t)(oi * ss + a) * Wss + (oj * ss + b)];
                        uint8_t r, g, bl; getColor(v, r, g, bl, cmethod);
                        lr += g_srgb2lin[r]; lg += g_srgb2lin[g]; lb += g_srgb2lin[bl];
                    }
                uint8_t* p = &out[((size_t)(obase + oi) * W + oj) * 3];
                p[0] = (uint8_t)(srgbEncode(lr / n) * 255.0 + 0.5);
                p[1] = (uint8_t)(srgbEncode(lg / n) * 255.0 + 0.5);
                p[2] = (uint8_t)(srgbEncode(lb / n) * 255.0 + 0.5);
            }
        }
        st->progress = (float)(s + 1) / nstrips;
    }
    { std::lock_guard<std::mutex> lk(st->curMx); st->cur = nullptr; }
    mpf_clear(scale_e);
    st->ok = !st->cancel;
    st->done = true;
}

// Write W x H top-down RGB to a 24-bit BMP (stored bottom-up, BGR).
static bool writeExportBMP(const wchar_t* path, const std::vector<uint8_t>& rgb, int W, int H) {
    int row = (W * 3 + 3) & ~3, dataSize = row * H;
    uint8_t fh[14] = { 'B','M' }; uint32_t v;
    v = 14 + 40 + dataSize; memcpy(fh + 2, &v, 4);
    v = 54; memcpy(fh + 10, &v, 4);
    uint8_t ih[40] = { 0 };
    v = 40; memcpy(ih + 0, &v, 4);
    v = (uint32_t)W; memcpy(ih + 4, &v, 4);
    v = (uint32_t)H; memcpy(ih + 8, &v, 4);
    uint16_t planes = 1, bpp = 24; memcpy(ih + 12, &planes, 2); memcpy(ih + 14, &bpp, 2);
    v = (uint32_t)dataSize; memcpy(ih + 20, &v, 4);
    FILE* f = _wfopen(path, L"wb"); if (!f) return false;
    fwrite(fh, 1, 14, f); fwrite(ih, 1, 40, f);
    std::vector<uint8_t> line(row, 0);
    for (int y = H - 1; y >= 0; --y) {
        const uint8_t* src = &rgb[(size_t)y * W * 3];
        for (int x = 0; x < W; ++x) { line[x*3+0]=src[x*3+2]; line[x*3+1]=src[x*3+1]; line[x*3+2]=src[x*3+0]; }
        fwrite(line.data(), 1, row, f);
    }
    fclose(f);
    return true;
}

// ---- Export dialog (native controls popup) -------------------------------
enum {
    IDC_RES = 2001, IDC_SS, IDC_WEDIT, IDC_HEDIT, IDC_EXPORT, IDC_CANCEL, IDC_PROG
};
struct ResPreset { const wchar_t* name; int w, h; };

struct ExportDlg {
    HWND hwnd = nullptr, owner = nullptr;
    MandelNavigator* nav = nullptr;
    int cmethod = 0, mxit = 0;
    mpf_t cx, cy, scale;
    HWND resCombo = 0, ssCombo = 0, wEdit = 0, hEdit = 0, exportBtn = 0, cancelBtn = 0, prog = 0;
    std::vector<ResPreset> presets;               // last entry is Custom (w=h=0)
    RECT pvRect{};                                 // preview area (client coords)
    std::vector<uint8_t> preview; int pvW = 0, pvH = 0;   // BGR top-down for StretchDIBits
    std::thread pvThread; ExportState pvState; std::atomic<int> pvGen{ 0 };
    std::thread exThread; ExportState exState; std::vector<uint8_t> exImg;
    bool exporting = false; int outW = 1920, outH = 1280, ss = 2;
    ExportDlg() { mpf_init(cx); mpf_init(cy); mpf_init(scale); }
    ~ExportDlg() {
        { std::lock_guard<std::mutex> lk(pvState.curMx); pvState.cancel = true; Mandel* mm = pvState.cur.load(); if (mm) mm->SetHalt(true); }
        if (pvThread.joinable()) pvThread.join();
        { std::lock_guard<std::mutex> lk(exState.curMx); exState.cancel = true; Mandel* mm = exState.cur.load(); if (mm) mm->SetHalt(true); }
        if (exThread.joinable()) exThread.join();
        mpf_clear(cx); mpf_clear(cy); mpf_clear(scale);
    }
};

// Standard sizes; only those fitting the monitor are offered (plus Custom).
static void buildResPresets(ExportDlg* d) {
    static const ResPreset all[] = {
        { L"1350 x 900  (view 3:2)", 1350, 900 },
        { L"1920 x 1080  (1080p)", 1920, 1080 },
        { L"2560 x 1440  (1440p)", 2560, 1440 },
        { L"2700 x 1800  (3:2)", 2700, 1800 },
        { L"3840 x 2160  (4K)", 3840, 2160 },
        { L"5120 x 2880  (5K)", 5120, 2880 },
        { L"7680 x 4320  (8K)", 7680, 4320 },
    };
    int mw = GetSystemMetrics(SM_CXSCREEN), mh = GetSystemMetrics(SM_CYSCREEN);
    d->presets.clear();
    for (const auto& r : all) if (r.w <= mw && r.h <= mh) d->presets.push_back(r);
    if (d->presets.empty()) d->presets.push_back({ L"1350 x 900  (view 3:2)", 1350, 900 });
    d->presets.push_back({ L"Custom\u2026", 0, 0 });
}

static void exportGetSize(ExportDlg* d) {
    int sel = (int)SendMessageW(d->resCombo, CB_GETCURSEL, 0, 0);
    if (sel < 0) sel = 0;
    if (d->presets[sel].w > 0) { d->outW = d->presets[sel].w; d->outH = d->presets[sel].h; }
    else {
        wchar_t buf[16];
        GetWindowTextW(d->wEdit, buf, 16); d->outW = std::clamp(_wtoi(buf), 16, 30000);
        GetWindowTextW(d->hEdit, buf, 16); d->outH = std::clamp(_wtoi(buf), 16, 30000);
    }
    int sssel = (int)SendMessageW(d->ssCombo, CB_GETCURSEL, 0, 0);
    d->ss = std::clamp(sssel + 1, 1, 4);
}

// Kick off a small preview render on a background thread (cancels any prior one).
static void startPreview(ExportDlg* d) {
    { std::lock_guard<std::mutex> lk(d->pvState.curMx); d->pvState.cancel = true; Mandel* mm = d->pvState.cur.load(); if (mm) mm->SetHalt(true); }
    if (d->pvThread.joinable()) d->pvThread.join();
    exportGetSize(d);
    int pw = 300, ph = std::max(1, (int)std::lround(300.0 * d->outH / d->outW));
    int gen = ++d->pvGen;
    d->pvState.cancel = false;
    HWND hw = d->hwnd;
    d->pvThread = std::thread([d, pw, ph, gen, hw]() {
        std::vector<uint8_t> rgb;
        exportRender(d->cx, d->cy, d->scale, pw, ph, 1, d->mxit, d->cmethod, rgb, &d->pvState);
        if (gen != d->pvGen || d->pvState.cancel) return;
        // to BGR top-down for StretchDIBits
        std::vector<uint8_t> bgr((size_t)pw * ph * 3);
        for (size_t i = 0; i < (size_t)pw * ph; ++i) { bgr[i*3]=rgb[i*3+2]; bgr[i*3+1]=rgb[i*3+1]; bgr[i*3+2]=rgb[i*3]; }
        d->preview.swap(bgr); d->pvW = pw; d->pvH = ph;
        InvalidateRect(hw, nullptr, FALSE);
    });
}

LRESULT CALLBACK ExportWndProc(HWND h, UINT m, WPARAM wp, LPARAM lp);

static void showExportDialog(HWND owner, MandelNavigator* nav) {
    static bool registered = false;
    HINSTANCE inst = (HINSTANCE)GetWindowLongPtrW(owner, GWLP_HINSTANCE);
    if (!registered) {
        INITCOMMONCONTROLSEX icc{ sizeof(icc), ICC_PROGRESS_CLASS | ICC_STANDARD_CLASSES };
        InitCommonControlsEx(&icc);
        WNDCLASSW wc{}; wc.lpfnWndProc = ExportWndProc; wc.hInstance = inst;
        wc.hCursor = LoadCursor(nullptr, IDC_ARROW);
        wc.hbrBackground = (HBRUSH)(COLOR_BTNFACE + 1);
        wc.lpszClassName = L"MandelExportDlg";
        RegisterClassW(&wc); registered = true;
    }
    ExportDlg* d = new ExportDlg();
    d->owner = owner; d->nav = nav;
    d->cmethod = nav->GetCMethod(); d->mxit = nav->GetMxit();
    nav->GetView(d->cx, d->cy, d->scale);
    buildResPresets(d);
    UINT dpi = GetDpiForWindow(owner); if (!dpi) dpi = 96;
    int W = MulDiv(560, dpi, 96), H = MulDiv(470, dpi, 96);
    RECT orc; GetWindowRect(owner, &orc);
    int x = orc.left + 80, y = orc.top + 80;
    HWND hw = CreateWindowExW(0, L"MandelExportDlg", L"Export image",
        WS_POPUPWINDOW | WS_CAPTION | WS_VISIBLE, x, y, W, H, owner, nullptr, inst, d);
    EnableWindow(owner, FALSE);
    d->hwnd = hw;
}

static void exportDoSave(ExportDlg* d) {
    wchar_t path[MAX_PATH] = L"mandelbrot.bmp";
    OPENFILENAMEW ofn{}; ofn.lStructSize = sizeof(ofn); ofn.hwndOwner = d->hwnd;
    ofn.lpstrFilter = L"Bitmap (*.bmp)\0*.bmp\0"; ofn.lpstrFile = path; ofn.nMaxFile = MAX_PATH;
    ofn.lpstrDefExt = L"bmp"; ofn.Flags = OFN_OVERWRITEPROMPT | OFN_PATHMUSTEXIST;
    if (!GetSaveFileNameW(&ofn)) return;
    exportGetSize(d);
    d->exState.cancel = false; d->exState.done = false; d->exState.progress = 0;
    d->exporting = true;
    EnableWindow(d->exportBtn, FALSE); EnableWindow(d->resCombo, FALSE);
    EnableWindow(d->ssCombo, FALSE); EnableWindow(d->wEdit, FALSE); EnableWindow(d->hEdit, FALSE);
    SetWindowTextW(d->cancelBtn, L"Cancel");
    std::wstring save = path;
    int W = d->outW, H = d->outH, ss = d->ss;
    if (d->exThread.joinable()) d->exThread.join();
    d->exThread = std::thread([d, W, H, ss, save]() {
        exportRender(d->cx, d->cy, d->scale, W, H, ss, d->mxit, d->cmethod, d->exImg, &d->exState);
        if (d->exState.ok) writeExportBMP(save.c_str(), d->exImg, W, H);
    });
}

LRESULT CALLBACK ExportWndProc(HWND h, UINT m, WPARAM wp, LPARAM lp) {
    ExportDlg* d = (ExportDlg*)GetWindowLongPtrW(h, GWLP_USERDATA);
    switch (m) {
    case WM_CREATE: {
        CREATESTRUCTW* cs = (CREATESTRUCTW*)lp;
        d = (ExportDlg*)cs->lpCreateParams;
        SetWindowLongPtrW(h, GWLP_USERDATA, (LONG_PTR)d);
        HINSTANCE inst = cs->hInstance;
        UINT dpi = GetDpiForWindow(h); if (!dpi) dpi = 96;
        auto S = [dpi](int v) { return MulDiv(v, dpi, 96); };
        HFONT font = (HFONT)GetStockObject(DEFAULT_GUI_FONT);
        auto mk = [&](const wchar_t* cls, const wchar_t* txt, DWORD style, int x, int y, int w, int hh, int id) -> HWND {
            HWND c = CreateWindowExW(0, cls, txt, WS_CHILD | WS_VISIBLE | style,
                S(x), S(y), S(w), S(hh), h, (HMENU)(INT_PTR)id, inst, nullptr);
            SendMessageW(c, WM_SETFONT, (WPARAM)font, TRUE); return c;
        };
        mk(L"STATIC", L"Resolution", 0, 16, 14, 100, 18, 0);
        d->resCombo = mk(L"COMBOBOX", L"", CBS_DROPDOWNLIST | WS_VSCROLL, 16, 34, 220, 240, IDC_RES);
        for (auto& p : d->presets) SendMessageW(d->resCombo, CB_ADDSTRING, 0, (LPARAM)p.name);
        SendMessageW(d->resCombo, CB_SETCURSEL, 0, 0);
        mk(L"STATIC", L"Custom  W x H", 0, 250, 14, 120, 18, 0);
        d->wEdit = mk(L"EDIT", L"1920", WS_BORDER | ES_NUMBER | ES_RIGHT, 250, 34, 64, 22, IDC_WEDIT);
        mk(L"STATIC", L"x", SS_CENTER, 316, 36, 12, 18, 0);
        d->hEdit = mk(L"EDIT", L"1280", WS_BORDER | ES_NUMBER | ES_RIGHT, 330, 34, 64, 22, IDC_HEDIT);
        mk(L"STATIC", L"Supersampling", 0, 410, 14, 110, 18, 0);
        d->ssCombo = mk(L"COMBOBOX", L"", CBS_DROPDOWNLIST, 410, 34, 116, 160, IDC_SS);
        for (const wchar_t* s : { L"1x (fast)", L"2x", L"3x", L"4x (best)" })
            SendMessageW(d->ssCombo, CB_ADDSTRING, 0, (LPARAM)s);
        SendMessageW(d->ssCombo, CB_SETCURSEL, 1, 0);
        d->pvRect = { S(16), S(74), S(544), S(360) };
        d->prog = mk(L"msctls_progress32", L"", 0, 16, 372, 528, 16, IDC_PROG);
        d->exportBtn = mk(L"BUTTON", L"Export\u2026", BS_DEFPUSHBUTTON, 344, 398, 92, 30, IDC_EXPORT);
        d->cancelBtn = mk(L"BUTTON", L"Close", 0, 444, 398, 92, 30, IDC_CANCEL);
        EnableWindow(d->wEdit, FALSE); EnableWindow(d->hEdit, FALSE);
        SetTimer(h, 1, 100, nullptr);
        startPreview(d);
        return 0;
    }
    case WM_COMMAND: {
        int id = LOWORD(wp), code = HIWORD(wp);
        if (id == IDC_RES && code == CBN_SELCHANGE) {
            int sel = (int)SendMessageW(d->resCombo, CB_GETCURSEL, 0, 0);
            bool custom = (sel >= 0 && d->presets[sel].w == 0);
            EnableWindow(d->wEdit, custom); EnableWindow(d->hEdit, custom);
            startPreview(d);
        } else if ((id == IDC_WEDIT || id == IDC_HEDIT) && code == EN_KILLFOCUS) {
            startPreview(d);
        } else if (id == IDC_EXPORT && code == BN_CLICKED) {
            if (!d->exporting) exportDoSave(d);
        } else if (id == IDC_CANCEL && code == BN_CLICKED) {
            if (d->exporting) { d->exState.cancel = true; std::lock_guard<std::mutex> lk(d->exState.curMx); Mandel* mm = d->exState.cur.load(); if (mm) mm->SetHalt(true); }
            else DestroyWindow(h);
        }
        return 0;
    }
    case WM_TIMER: {
        if (d->exporting) {
            SendMessageW(d->prog, PBM_SETPOS, (int)(d->exState.progress * 100.0f + 0.5f), 0);
            if (d->exState.done) {
                if (d->exThread.joinable()) d->exThread.join();
                d->exporting = false;
                SendMessageW(d->prog, PBM_SETPOS, d->exState.ok ? 100 : 0, 0);
                EnableWindow(d->exportBtn, TRUE); EnableWindow(d->resCombo, TRUE);
                EnableWindow(d->ssCombo, TRUE);
                int sel = (int)SendMessageW(d->resCombo, CB_GETCURSEL, 0, 0);
                bool custom = (sel >= 0 && d->presets[sel].w == 0);
                EnableWindow(d->wEdit, custom); EnableWindow(d->hEdit, custom);
                SetWindowTextW(d->cancelBtn, L"Close");
                if (d->exState.ok) MessageBoxW(h, L"Saved.", L"Export", MB_OK | MB_ICONINFORMATION);
            }
        }
        return 0;
    }
    case WM_PAINT: {
        PAINTSTRUCT ps; HDC dc = BeginPaint(h, &ps);
        RECT r = d->pvRect;
        FrameRect(dc, &r, (HBRUSH)GetStockObject(GRAY_BRUSH));
        wchar_t info[64]; swprintf_s(info, L"%d x %d,  %dx SS", d->outW, d->outH, d->ss);
        RECT tr = { r.left, r.bottom + 4, r.right, r.bottom + 24 };
        SetBkMode(dc, TRANSPARENT); DrawTextW(dc, info, -1, &tr, DT_CENTER | DT_TOP);
        if (d->pvW > 0 && !d->preview.empty()) {
            // aspect-fit the preview inside pvRect
            int bw = r.right - r.left - 4, bh = r.bottom - r.top - 4;
            double s = std::min((double)bw / d->pvW, (double)bh / d->pvH);
            int dw = (int)(d->pvW * s), dh = (int)(d->pvH * s);
            int dx = r.left + 2 + (bw - dw) / 2, dy = r.top + 2 + (bh - dh) / 2;
            BITMAPINFO bi{}; bi.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
            bi.bmiHeader.biWidth = d->pvW; bi.bmiHeader.biHeight = -d->pvH;
            bi.bmiHeader.biPlanes = 1; bi.bmiHeader.biBitCount = 24; bi.bmiHeader.biCompression = BI_RGB;
            SetStretchBltMode(dc, HALFTONE);
            StretchDIBits(dc, dx, dy, dw, dh, 0, 0, d->pvW, d->pvH, d->preview.data(), &bi, DIB_RGB_COLORS, SRCCOPY);
        } else {
            RECT cr = r; DrawTextW(dc, L"Rendering preview\u2026", -1, &cr, DT_CENTER | DT_VCENTER | DT_SINGLELINE);
        }
        EndPaint(h, &ps);
        return 0;
    }
    case WM_CLOSE:
        if (d && d->exporting) { d->exState.cancel = true; std::lock_guard<std::mutex> lk(d->exState.curMx); Mandel* mm = d->exState.cur.load(); if (mm) mm->SetHalt(true); }
        DestroyWindow(h);
        return 0;
    case WM_DESTROY:
        if (d) {
            KillTimer(h, 1);
            EnableWindow(d->owner, TRUE); SetForegroundWindow(d->owner);
            delete d;
            SetWindowLongPtrW(h, GWLP_USERDATA, 0);
        }
        return 0;
    }
    return DefWindowProcW(h, m, wp, lp);
}

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
    H_RESET, H_RENDER, H_SAVE, H_COPY, H_PASTE,
    H_MAXTRACK, H_DENSTRACK,
    H_SS, H_EDE, H_PALETTE_DD, H_COLOR, H_GALLERY_DD
};

// Gallery demo presets: a saved location plus the exact render settings used to
// produce it, so selecting one reproduces the image. zoom is scientific notation
// (expanded by expandSci); palette is looked up by name in palettePresets().
struct Preset {
    const wchar_t* name;
    const char* x; const char* y; const char* zoom;
    float density; bool ss; int coloring;   // 0 Smooth, 1 Distance(EDE), 2 Feather
    const wchar_t* palette; int maxIter;
};
static const std::vector<Preset>& galleryPresets() {
    static const std::vector<Preset> p = {
        { L"Night City  (1.7e40)",
          "-1.768628917759850520844734198472848718821423994141176532908",
          "0.001395534274228826510747662517373603005419245032078944176",
          "1.691960e40", 83.4f, true, 0, L"Sunrise", 500000 },
        { L"Golden Phoenix  (5.1e292)",
          "-1.74961551043225917132558762203092997406776582824486737043789087410512096670138841060878237427473435515004168931589322189249239986606828774155860523237102635565490176581179889267432026148118400022509532606522827603302653557653285809596137929818429125986636205433675211119215873019002373733544536685782529789935",
          "0.00000033552806488437213922924936314667758682179580368291347248996372066040948346974233057158256293565061238343085946164873110399056182528356343156076648247679264617131919657718521028825787075233133790754370605292961104472660681262156196225123453401177034629580032892768133871611413536923738263673919333840137",
          "5.119695e292", 169.0f, true, 2, L"Sunrise", 500000 },
    };
    return p;
}

struct Stop { float pos, r, g, b; };

struct NamedPal { const wchar_t* name; std::vector<Stop> stops; };

// Finalized palette presets (natural / gemstone; saturated hues bridged by
// white or near-black). Index 0 is the default.
static const std::vector<NamedPal>& palettePresets() {
    static const std::vector<NamedPal> p = {
        { L"Sunrise",  { {0,0,70,100},{0.16f,32,107,203},{0.42f,237,255,255},{0.6425f,255,170,0},{0.8575f,0,2,0} } },
        { L"Reef",     { {0,0,45,65},{0.15f,25,155,195},{0.40f,240,255,255},{0.60f,255,190,45},{0.78f,220,60,40},{0.90f,12,3,10} } },
        { L"Emerald",  { {0,6,45,30},{0.16f,40,160,90},{0.42f,235,255,235},{0.6425f,225,205,120},{0.8575f,6,18,12} } },
        { L"Amethyst", { {0,32,12,60},{0.18f,120,70,180},{0.40f,200,175,230},{0.58f,245,242,252},{0.8575f,18,8,30} } },
        { L"Topaz",    { {0,35,20,8},{0.18f,170,95,30},{0.40f,235,175,70},{0.60f,250,235,200},{0.85f,20,12,6} } },
        { L"Gems",     { {0.05f,30,90,200},{0.17f,242,248,255},{0.30f,25,155,95},{0.42f,10,12,20},{0.55f,125,70,185},{0.67f,250,244,252},{0.80f,240,180,55},{0.92f,18,12,8} } },
        { L"Ocean",    { {0,4,12,38},{0.22f,10,64,120},{0.46f,36,150,190},{0.66f,150,220,230},{0.85f,244,252,255} } },
        { L"Lava",     { {0,10,8,20},{0.16f,140,25,30},{0.36f,240,110,30},{0.56f,250,215,120},{0.72f,180,210,230},{0.88f,50,70,140} } },
        { L"Snowy",    { {0,255,255,255},{0.12f,32,107,203},{0.34f,214,223,255} } },
    };
    return p;
}

struct PaletteEditor {
    std::vector<Stop> stops;
    int selected = 0;
    bool dragging = false;

    void load(int idx) {
        const auto& pr = palettePresets();
        if (idx < 0 || idx >= (int)pr.size()) return;
        stops = pr[idx].stops; selected = 0; rebuild();
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

// Cap the x:/y: lines of the location text to at most 2 mono lines each (cpl
// chars/line), truncating with "..." so the zoom line is always visible.
static std::wstring truncLocation(const std::string& raw, int cpl) {
    if (cpl < 8) cpl = 8;
    std::wstring out; std::string acc;
    auto flush = [&](bool last) {
        std::string L = acc; acc.clear();
        if ((L.rfind("x:", 0) == 0 || L.rfind("y:", 0) == 0) && (int)L.size() > cpl) {
            std::string l1 = L.substr(0, cpl), rest = L.substr(cpl);
            if ((int)rest.size() > cpl) rest = rest.substr(0, std::max(0, cpl - 3)) + "...";
            L = l1 + "\n" + rest;
        }
        for (char c : L) out.push_back(c == '\n' ? L'\n' : (wchar_t)(unsigned char)c);
        if (!last) out.push_back(L'\n');
    };
    for (char c : raw) { if (c == '\r') continue; if (c == '\n') flush(false); else acc.push_back(c); }
    flush(true);
    return out;
}

// Expand a scientific-notation number ("4.361669e+03") to a plain decimal string
// ("4361.669") so GMP mpf_set_str (base 10, no exponent) can parse it.
static std::string expandSci(const std::string& in) {
    std::string s = in;
    size_t ep = s.find_first_of("eE");
    std::string mant = (ep == std::string::npos) ? s : s.substr(0, ep);
    long e10 = (ep == std::string::npos) ? 0 : atol(s.substr(ep + 1).c_str());
    std::string sign;
    if (!mant.empty() && (mant[0] == '+' || mant[0] == '-')) { if (mant[0] == '-') sign = "-"; mant.erase(0, 1); }
    size_t dp = mant.find('.');
    int frac = (dp == std::string::npos) ? 0 : (int)(mant.size() - dp - 1);
    if (dp != std::string::npos) mant.erase(dp, 1);
    if (mant.empty()) return "";
    long shift = e10 - frac;                 // value = mant (integer) * 10^shift
    std::string out;
    if (shift >= 0) { out = mant; out.append((size_t)shift, '0'); }
    else {
        long point = (long)mant.size() + shift;   // decimal point position from the left
        if (point <= 0) out = "0." + std::string((size_t)(-point), '0') + mant;
        else out = mant.substr(0, (size_t)point) + "." + mant.substr((size_t)point);
    }
    return sign + out;
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
    RECT rcReset{}, rcRender{}, rcSave{}, rcCopy{}, rcPaste{};
    RECT rcLocation{}, rcMaxTrack{}, rcDensTrack{};
    RECT rcSS{}, rcColoringDD{}, rcPaletteDD{}, rcColor{}, rcGradient{};
    RECT rcGalleryDD{};

    // state
    int maxIter = 500000;
    bool ssOn = false;
    int coloringIdx = 0;                    // 0 Smooth, 1 Distance (EDE), 2 Feather (SAC)
    bool coloringOpen = false;              // coloring dropdown expanded
    int coloringHover = -1;                 // hovered item while open (-1 none)
    int paletteIdx = 0;                    // index into palettePresets()
    bool paletteOpen = false;              // dropdown expanded
    int paletteHover = -1;                 // hovered item while open (-1 none)
    bool galleryOpen = false;              // gallery dropdown expanded
    int galleryHover = -1;                 // hovered gallery item (-1 none)
    bool navDragging = false, wasComputing = false;
    Hit hover = H_NONE, pressed = H_NONE;
    int liveFrames = 0;   // frames still needing per-tick repaint (animation/compute)
    int dpi = 96;         // display DPI; all metrics scale by dpi/96
    std::chrono::steady_clock::time_point renderStart;
    double lastRenderMs = 0;
    // Cached back buffer + partial-repaint state (avoids reallocating the buffer
    // and rebuilding the whole UI panel on every drag/zoom/compute tick).
    HDC memDC = nullptr; HBITMAP memBmp = nullptr; HGDIOBJ memOld = nullptr;
    int memW = 0, memH = 0;
    bool needFull = true;          // force a full redraw (buffer (re)created / first paint)
    bool fractalOnlyTick = false;  // set by the timer for pure fractal animation frames
    double lastPresentMs = 0;

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
        int bw = (w - 4 * g) / 5;
        rcReset   = { px,              y, px + bw,          y + bh };
        rcRender  = { px + bw + g,     y, px + 2*bw + g,    y + bh };
        rcSave    = { px + 2*(bw+g),   y, px + 3*bw + 2*g,  y + bh };
        rcCopy    = { px + 3*(bw+g),   y, px + 4*bw + 3*g,  y + bh };
        rcPaste   = { px + 4*(bw+g),   y, px + 5*bw + 4*g,  y + bh };
        y += bh + S(14);
        rcLocation = { px, y, px + w, y + S(100) }; y += S(100) + S(16);
        rcMaxTrack = { px, y + S(24), px + w, y + S(38) }; y += S(52);
        rcDensTrack = { px, y + S(24), px + w, y + S(38) }; y += S(52);
        rcSS  = { px, y, px + w, y + bh }; y += S(56);
        rcColoringDD = { px, y, px + w, y + bh }; y += S(56);
        rcPaletteDD = { px, y, px + w, y + bh }; y += S(44);
        rcGradient = { px, y + S(22), px + w, y + S(62) }; y += S(84);
        rcColor = { px, y, px + w, y + bh }; y += bh + S(16);
        rcGalleryDD = { px, y, px + w, y + bh };
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
    // Max-iteration slider: log scale from 100 (10^2) to 5,000,000.
    static double maxItExpSpan() { return log10(5000000.0) - 2.0; }   // ~4.699
    double maxToT() const { return std::clamp((log10((double)maxIter) - 2.0) / maxItExpSpan(), 0.0, 1.0); }
    void setMaxFromT(double t, bool render) {
        maxIter = (int)std::round(pow(10.0, 2.0 + maxItExpSpan() * std::clamp(t, 0.0, 1.0)));
        maxIter = std::clamp(maxIter, 100, 5000000);
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
    void pasteLocation() {
        if (!OpenClipboard(hwnd)) return;
        std::string txt;
        if (HANDLE h = GetClipboardData(CF_UNICODETEXT)) {
            if (wchar_t* p = (wchar_t*)GlobalLock(h)) {
                for (wchar_t* q = p; *q; ++q) txt.push_back(*q < 128 ? (char)*q : ' ');
                GlobalUnlock(h);
            }
        }
        CloseClipboard();
        if (txt.empty()) return;
        auto val = [&](const char* key) -> std::string {
            size_t p = txt.find(key); if (p == std::string::npos) return "";
            p += strlen(key);
            while (p < txt.size() && (txt[p] == ' ' || txt[p] == '\t')) ++p;
            size_t e = p; while (e < txt.size() && txt[e] != '\r' && txt[e] != '\n') ++e;
            std::string v = txt.substr(p, e - p);
            while (!v.empty() && (v.back() == ' ' || v.back() == '\t')) v.pop_back();
            return v;
        };
        std::string xs = val("x:"), ys = val("y:"), zs = val("zoom:");
        if (xs.empty() || ys.empty() || zs.empty()) return;
        std::string scale = expandSci(zs);
        if (scale.empty()) return;
        if (nav->SetLocation(xs, ys, scale)) startRender();
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
    // A slider caption: dim label on the left, brighter value right-aligned, on
    // the line just above the track (used by Max iterations + Color density).
    void labelRow(HDC dc, const RECT& track, const std::wstring& lab, const std::wstring& val) {
        RECT r = { track.left, track.top - S(22), track.right, track.top - S(4) };
        drawText(dc, r, lab, CLR_TEXT_DIM, fSmall, DT_LEFT | DT_TOP | DT_SINGLELINE);
        drawText(dc, r, val, CLR_TEXT, fSmall, DT_RIGHT | DT_TOP | DT_SINGLELINE);
    }

    // ---- palette dropdown ----
    static std::array<int,3> sampleStops(const std::vector<Stop>& st, float t) {
        if (st.empty()) return { 0,0,0 };
        if (t <= st.front().pos) return { (int)st.front().r,(int)st.front().g,(int)st.front().b };
        if (t >= st.back().pos)  return { (int)st.back().r,(int)st.back().g,(int)st.back().b };
        for (size_t i = 0; i + 1 < st.size(); ++i)
            if (t >= st[i].pos && t <= st[i+1].pos) {
                float f = (t - st[i].pos) / (st[i+1].pos - st[i].pos + 1e-6f);
                return { (int)(st[i].r + (st[i+1].r - st[i].r)*f),
                         (int)(st[i].g + (st[i+1].g - st[i].g)*f),
                         (int)(st[i].b + (st[i+1].b - st[i].b)*f) };
            }
        return { (int)st.back().r,(int)st.back().g,(int)st.back().b };
    }
    void drawSwatch(HDC dc, RECT r, const std::vector<Stop>& st) {
        auto cl = [](int v){ return v < 0 ? 0 : (v > 255 ? 255 : v); };
        int w = r.right - r.left;
        for (int x = 0; x < w; ++x) {
            auto c = sampleStops(st, (float)x / (w > 1 ? w - 1 : 1));
            RECT col = { r.left + x, r.top, r.left + x + 1, r.bottom };
            HBRUSH b = CreateSolidBrush(RGB(cl(c[0]), cl(c[1]), cl(c[2])));
            FillRect(dc, &col, b); DeleteObject(b);
        }
    }
    int paletteItemH() const { return S(28); }
    RECT paletteListRect() const {
        int n = (int)palettePresets().size();
        return { rcPaletteDD.left, rcPaletteDD.bottom + S(2), rcPaletteDD.right,
                 rcPaletteDD.bottom + S(2) + n * paletteItemH() };
    }
    int paletteItemAt(int x, int y) const {
        RECT lr = paletteListRect();
        if (!inRect(lr, x, y)) return -1;
        int i = (y - lr.top) / paletteItemH();
        return (i >= 0 && i < (int)palettePresets().size()) ? i : -1;
    }
    // ---- coloring-mode dropdown (Smooth / Distance / Feather) ----
    static const wchar_t* coloringName(int i) {
        static const wchar_t* n[3] = { L"Smooth", L"Distance (EDE)", L"Feather (stripe)" };
        return n[i < 0 ? 0 : (i > 2 ? 2 : i)];
    }
    int coloringItemH() const { return S(28); }
    RECT coloringListRect() const {
        int n = 3;
        return { rcColoringDD.left, rcColoringDD.bottom + S(2), rcColoringDD.right,
                 rcColoringDD.bottom + S(2) + n * coloringItemH() };
    }
    int coloringItemAt(int x, int y) const {
        RECT lr = coloringListRect();
        if (x < lr.left || x > lr.right || y < lr.top || y > lr.bottom) return -1;
        int i = (y - lr.top) / coloringItemH();
        return (i < 0 || i > 2) ? -1 : i;
    }
    void drawColoringDD(HDC dc) {
        bool hov = hover == H_EDE;
        fillRound(dc, rcColoringDD, hov ? CLR_CARD_HOV : CLR_CARD,
                  coloringOpen ? CLR_ACCENT : CLR_BORDER, S(8));
        RECT tr = rcColoringDD; tr.left += S(12); tr.right -= S(28);
        drawText(dc, tr, coloringName(coloringIdx), CLR_TEXT, fUi, DT_LEFT | DT_VCENTER | DT_SINGLELINE);
        RECT cr = rcColoringDD; cr.right -= S(12);
        drawText(dc, cr, coloringOpen ? L"\u25B2" : L"\u25BC", CLR_TEXT_DIM, fSmall, DT_RIGHT | DT_VCENTER | DT_SINGLELINE);
    }
    void drawColoringList(HDC dc) {
        RECT lr = coloringListRect();
        fillRound(dc, lr, CLR_CARD, CLR_ACCENT, S(8));
        int ih = coloringItemH();
        for (int i = 0; i < 3; ++i) {
            RECT ir = { lr.left, lr.top + i * ih, lr.right, lr.top + (i + 1) * ih };
            if (i == coloringHover) fillRect(dc, { ir.left + S(3), ir.top, ir.right - S(3), ir.bottom }, CLR_CARD_HOV);
            RECT tr = ir; tr.left += S(14);
            drawText(dc, tr, coloringName(i), i == coloringIdx ? CLR_ACCENT : CLR_TEXT, fUi, DT_LEFT | DT_VCENTER | DT_SINGLELINE);
        }
    }
    void selectColoring(int idx) {
        if (idx < 0 || idx > 2) return;
        coloringIdx = idx;
        int m = nav->GetCMethod();
        m &= ~(ColoringMethod::EXTERIOR_DIST_EST | ColoringMethod::STRIPE_AVERAGE);
        if (idx == 1) m |= ColoringMethod::EXTERIOR_DIST_EST;
        else if (idx == 2) m |= ColoringMethod::STRIPE_AVERAGE;
        nav->SetCMethod(m);
        startRender();
    }

    // ---- gallery (demo presets) ----
    int galleryItemH() const { return S(28); }
    RECT galleryListRect() const {
        int n = (int)galleryPresets().size();
        int h = n * galleryItemH();
        // Gallery sits at the panel bottom, so the list opens upward.
        return { rcGalleryDD.left, rcGalleryDD.top - S(2) - h, rcGalleryDD.right, rcGalleryDD.top - S(2) };
    }
    int galleryItemAt(int x, int y) const {
        RECT lr = galleryListRect();
        if (x < lr.left || x > lr.right || y < lr.top || y > lr.bottom) return -1;
        int i = (y - lr.top) / galleryItemH();
        return (i >= 0 && i < (int)galleryPresets().size()) ? i : -1;
    }
    void drawGalleryDD(HDC dc) {
        bool hov = hover == H_GALLERY_DD;
        fillRound(dc, rcGalleryDD, hov ? CLR_CARD_HOV : CLR_CARD,
                  galleryOpen ? CLR_ACCENT : CLR_BORDER, S(8));
        RECT tr = rcGalleryDD; tr.left += S(12); tr.right -= S(28);
        drawText(dc, tr, L"Gallery", CLR_TEXT, fUi, DT_LEFT | DT_VCENTER | DT_SINGLELINE);
        RECT cr = rcGalleryDD; cr.right -= S(12);
        drawText(dc, cr, galleryOpen ? L"\u25B2" : L"\u25BC", CLR_TEXT_DIM, fSmall, DT_RIGHT | DT_VCENTER | DT_SINGLELINE);
    }
    void drawGalleryList(HDC dc) {
        RECT lr = galleryListRect();
        fillRound(dc, lr, CLR_CARD, CLR_ACCENT, S(8));
        const auto& gp = galleryPresets();
        int ih = galleryItemH();
        for (int i = 0; i < (int)gp.size(); ++i) {
            RECT ir = { lr.left, lr.top + i * ih, lr.right, lr.top + (i + 1) * ih };
            if (i == galleryHover) fillRect(dc, { ir.left + S(3), ir.top, ir.right - S(3), ir.bottom }, CLR_CARD_HOV);
            RECT tr = ir; tr.left += S(14);
            drawText(dc, tr, gp[i].name, CLR_TEXT, fUi, DT_LEFT | DT_VCENTER | DT_SINGLELINE);
        }
    }
    void applyPreset(int idx) {
        const auto& gp = galleryPresets();
        if (idx < 0 || idx >= (int)gp.size()) return;
        const Preset& p = gp[idx];
        // render settings first, then the location (which kicks off the render)
        maxIter = std::clamp(p.maxIter, 100, 5000000); nav->SetMxit(maxIter);
        color_density = std::clamp(p.density, 10.0f, 200.0f);
        ssOn = p.ss;
        coloringIdx = p.coloring;
        int m = 0;
        if (ssOn) m |= ColoringMethod::SUPER_SAMPLING;
        if (coloringIdx == 1) m |= ColoringMethod::EXTERIOR_DIST_EST;
        else if (coloringIdx == 2) m |= ColoringMethod::STRIPE_AVERAGE;
        nav->SetCMethod(m);
        // palette by name
        const auto& pr = palettePresets();
        for (int i = 0; i < (int)pr.size(); ++i)
            if (wcscmp(pr[i].name, p.palette) == 0) { paletteIdx = i; palette.load(i); break; }
        std::string scale = expandSci(p.zoom);
        if (!scale.empty() && nav->SetLocation(p.x, p.y, scale)) startRender();
        InvalidateRect(hwnd, nullptr, FALSE);
    }

    void drawPaletteDD(HDC dc) {
        bool hov = hover == H_PALETTE_DD;
        fillRound(dc, rcPaletteDD, hov ? CLR_CARD_HOV : CLR_CARD,
                  paletteOpen ? CLR_ACCENT : CLR_BORDER, S(8));
        RECT sw = { rcPaletteDD.left + S(8), rcPaletteDD.top + S(6), rcPaletteDD.left + S(52), rcPaletteDD.bottom - S(6) };
        drawSwatch(dc, sw, palette.stops);
        RECT tr = rcPaletteDD; tr.left += S(62); tr.right -= S(28);
        drawText(dc, tr, palettePresets()[paletteIdx].name, CLR_TEXT, fUi, DT_LEFT | DT_VCENTER | DT_SINGLELINE);
        RECT cr = rcPaletteDD; cr.right -= S(12);
        drawText(dc, cr, paletteOpen ? L"\u25B2" : L"\u25BC", CLR_TEXT_DIM, fSmall, DT_RIGHT | DT_VCENTER | DT_SINGLELINE);
    }
    void drawPaletteList(HDC dc) {
        RECT lr = paletteListRect();
        fillRound(dc, lr, CLR_CARD, CLR_ACCENT, S(8));
        const auto& pr = palettePresets();
        int ih = paletteItemH();
        for (int i = 0; i < (int)pr.size(); ++i) {
            RECT ir = { lr.left, lr.top + i * ih, lr.right, lr.top + (i + 1) * ih };
            if (i == paletteHover) fillRect(dc, { ir.left + S(3), ir.top, ir.right - S(3), ir.bottom }, CLR_CARD_HOV);
            RECT sw = { ir.left + S(10), ir.top + S(5), ir.left + S(46), ir.bottom - S(5) };
            drawSwatch(dc, sw, pr[i].stops);
            RECT tr = ir; tr.left += S(56);
            drawText(dc, tr, pr[i].name, i == paletteIdx ? CLR_ACCENT : CLR_TEXT, fUi, DT_LEFT | DT_VCENTER | DT_SINGLELINE);
        }
    }
    void selectPalette(int idx) {
        const auto& pr = palettePresets();
        if (idx < 0 || idx >= (int)pr.size()) return;
        paletteIdx = idx; palette.load(idx);
        nav->SetRedisplay(); keepLive(); InvalidateRect(hwnd, nullptr, FALSE);
    }

    void paint() {
        PAINTSTRUCT ps; HDC wdc = BeginPaint(hwnd, &ps);
        RECT rc; GetClientRect(hwnd, &rc);
        auto _pt0 = std::chrono::steady_clock::now();
        // Cached back buffer: recreate only when the client size changes, not per
        // frame. Any content change forces a full redraw via needFull.
        if (!memDC || memW != (int)rc.right || memH != (int)rc.bottom) {
            if (memDC) { SelectObject(memDC, memOld); DeleteObject(memBmp); DeleteDC(memDC); }
            memDC = CreateCompatibleDC(wdc);
            memBmp = CreateCompatibleBitmap(wdc, rc.right, rc.bottom);
            memOld = SelectObject(memDC, memBmp);
            memW = rc.right; memH = rc.bottom; needFull = true;
        }
        HDC dc = memDC;
        // Full redraw unless this is a pure fractal-animation tick from the timer.
        bool full = needFull || !fractalOnlyTick;
        needFull = false; fractalOnlyTick = false;
        if (full) fillRect(dc, rc, CLR_BG);   // bg + letterbox (full paints only)

        // fractal view (always -- it changes every tick)
        RECT vr = viewRect();
        BITMAPINFO bi{}; bi.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
        bi.bmiHeader.biWidth = RENDER_W; bi.bmiHeader.biHeight = -RENDER_H;
        bi.bmiHeader.biPlanes = 1; bi.bmiHeader.biBitCount = 24; bi.bmiHeader.biCompression = BI_RGB;
        SetStretchBltMode(dc, HALFTONE);
        StretchDIBits(dc, vr.left, vr.top, vr.right-vr.left, vr.bottom-vr.top,
                      0,0, RENDER_W, RENDER_H, display.data(), &bi, DIB_RGB_COLORS, SRCCOPY);

        // panel + all controls: only rebuilt when UI state changed. On fractal/
        // compute/drag ticks the panel persists in the cached buffer.
        if (full) {
        // panel
        RECT panel = { rc.right - S(PANEL_W), 0, rc.right, rc.bottom };
        fillRect(dc, panel, CLR_PANEL);
        HPEN sep = CreatePen(PS_SOLID, 1, CLR_BORDER);
        HGDIOBJ osep = SelectObject(dc, sep);
        MoveToEx(dc, panel.left, 0, nullptr); LineTo(dc, panel.left, rc.bottom);
        SelectObject(dc, osep); DeleteObject(sep);

        drawButton(dc, rcReset, L"Reset", H_RESET, false);
        drawButton(dc, rcRender, L"Render", H_RENDER, true);
        drawButton(dc, rcSave, L"Export", H_SAVE, false);
        drawButton(dc, rcCopy, L"Copy", H_COPY, false);
        drawButton(dc, rcPaste, L"Paste", H_PASTE, false);

        drawGalleryDD(dc);

        // location card -- cap x/y to 2 lines each so the zoom line always shows
        fillRound(dc, rcLocation, CLR_CARD, CLR_BORDER, S(8));
        RECT lt = rcLocation; lt.left += S(10); lt.top += S(6); lt.right -= S(10); lt.bottom -= S(6);
        {
            HGDIOBJ ofm = SelectObject(dc, fMono);
            SIZE cs{}; GetTextExtentPoint32W(dc, L"0", 1, &cs);
            SelectObject(dc, ofm);
            int cpl = cs.cx > 0 ? (int)((lt.right - lt.left) / cs.cx) : 40;
            drawText(dc, lt, truncLocation(nav->GetLocationText(), cpl), CLR_TEXT, fMono,
                     DT_LEFT | DT_TOP | DT_NOPREFIX);
        }

        // max iterations (read-only value; drag slider + mouse-wheel to change)
        labelRow(dc, rcMaxTrack, L"Max iterations", std::to_wstring(maxIter));
        drawSlider(dc, rcMaxTrack, maxToT(), H_MAXTRACK);

        // color density (read-only; drag slider = integer, mouse-wheel = 0.1)
        wchar_t db[32]; swprintf_s(db, L"%.1f", color_density);
        labelRow(dc, rcDensTrack, L"Color density", db);
        drawSlider(dc, rcDensTrack, std::clamp((color_density - 10.0) / 190.0, 0.0, 1.0), H_DENSTRACK);

        drawToggle(dc, rcSS, L"5x supersampling", ssOn, H_SS);
        label(dc, rcColoringDD.left, rcColoringDD.top - S(20), L"Coloring");
        drawColoringDD(dc);

        label(dc, rcPaletteDD.left, rcPaletteDD.top - S(20), L"Palette");
        drawPaletteDD(dc);

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
        // palette dropdown list -- drawn last so it overlays the widgets below it
        if (paletteOpen) drawPaletteList(dc);
        if (coloringOpen) drawColoringList(dc);
        if (galleryOpen) drawGalleryList(dc);
        } // end panel (full paints only)

        // status bar (always -- text changes with compute state / timings)
        RECT status = { 0, rc.bottom - S(STATUS_H), rc.right - S(PANEL_W), rc.bottom };
        fillRect(dc, status, CLR_PANEL);
        std::wstring st = nav->IsComputing()
            ? L"  Rendering..."
            : L"  Ready   last render " + std::to_wstring((int)lastRenderMs) + L" ms";
        st += L"   present " + std::to_wstring((int)(lastPresentMs + 0.5)) + L" ms";
        if (!nav->IsComputing())
            st += L"    Drag pan   Wheel zoom   R reset   Space render   S save   C copy";
        RECT str = status; str.left += S(10);
        drawText(dc, str, st, CLR_TEXT_DIM, fSmall, DT_LEFT | DT_VCENTER | DT_SINGLELINE);

        BitBlt(wdc, 0, 0, rc.right, rc.bottom, dc, 0, 0, SRCCOPY);
        EndPaint(hwnd, &ps);
        lastPresentMs = std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - _pt0).count();
    }

    Hit hitTest(int x, int y) {
        if (inRect(rcReset,x,y)) return H_RESET;
        if (inRect(rcRender,x,y)) return H_RENDER;
        if (inRect(rcSave,x,y)) return H_SAVE;
        if (inRect(rcCopy,x,y)) return H_COPY;
        if (inRect(rcPaste,x,y)) return H_PASTE;
        if (inRect(rcGalleryDD,x,y)) return H_GALLERY_DD;
        RECT mt = rcMaxTrack; mt.top -= S(8); mt.bottom += S(8); if (inRect(mt,x,y)) return H_MAXTRACK;
        RECT dt = rcDensTrack; dt.top -= S(8); dt.bottom += S(8); if (inRect(dt,x,y)) return H_DENSTRACK;
        if (inRect(rcSS,x,y)) return H_SS;
        if (inRect(rcColoringDD,x,y)) return H_EDE;
        if (inRect(rcPaletteDD,x,y)) return H_PALETTE_DD;
        if (inRect(rcColor,x,y)) return H_COLOR;
        RECT gr = rcGradient; gr.left -= S(8); gr.right += S(8); gr.top -= S(6); gr.bottom += S(14);
        if (inRect(gr,x,y)) return H_GRADIENT;
        RECT vr = viewRect(); if (inRect(vr,x,y)) return H_VIEW;
        return H_NONE;
    }

    void timer() {
        bool computing = nav->IsComputing();
        bool active = computing || wasComputing || navDragging || palette.dragging ||
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
        // Pure fractal-animation frames (pan/zoom/compute) can skip the panel
        // rebuild; UI-control interaction still gets a full paint.
        bool uiInteract = palette.dragging ||
                          pressed == H_MAXTRACK || pressed == H_DENSTRACK;
        fractalOnlyTick = !uiInteract;
        InvalidateRect(hwnd, nullptr, FALSE);
    }

    LRESULT message(UINT msg, WPARAM wp, LPARAM lp) {
        // Any non-timer/paint message may change UI state, so force the next paint
        // to be a full redraw (the timer re-enables fractal-only ticks).
        if (msg != WM_TIMER && msg != WM_PAINT) fractalOnlyTick = false;
        switch (msg) {
        case WM_CREATE: {
            dpi = (int)GetDpiForWindow(hwnd);
            if (dpi <= 0) dpi = 96;
            createFonts();
            palette.load(paletteIdx);
            nav = std::make_unique<MandelNavigator>(RENDER_W, RENDER_H, 5, 5000000, 1.0, 220.0);
            nav->SetMxit(maxIter);
            nav->SetCMethod(coloringIdx == 1 ? ColoringMethod::EXTERIOR_DIST_EST
                          : coloringIdx == 2 ? ColoringMethod::STRIPE_AVERAGE : 0);
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
            if (pressed == H_DENSTRACK) { color_density = (float)std::round(std::clamp(10.0 + 190.0*(x-rcDensTrack.left)/(rcDensTrack.right-rcDensTrack.left),10.0,200.0)); nav->SetRedisplay(); InvalidateRect(hwnd,nullptr,FALSE); return 0; }
            if (navDragging) { POINT p = mapToRender(x,y); nav->Drag(p.x,p.y); return 0; }
            if (paletteOpen) { int it = paletteItemAt(x,y); if (it != paletteHover) { paletteHover = it; InvalidateRect(hwnd,nullptr,FALSE); } }
            if (coloringOpen) { int it = coloringItemAt(x,y); if (it != coloringHover) { coloringHover = it; InvalidateRect(hwnd,nullptr,FALSE); } }
            if (galleryOpen) { int it = galleryItemAt(x,y); if (it != galleryHover) { galleryHover = it; InvalidateRect(hwnd,nullptr,FALSE); } }
            Hit h = hitTest(x,y);
            if (h != hover) { hover = h; InvalidateRect(hwnd, nullptr, FALSE); }
            return 0;
        }
        case WM_LBUTTONDOWN: {
            int x = GET_X_LPARAM(lp), y = GET_Y_LPARAM(lp);
            if (paletteOpen) {   // a click while the list is open selects or dismisses it
                int it = paletteItemAt(x,y);
                if (it >= 0) selectPalette(it);
                paletteOpen = false; paletteHover = -1; pressed = H_NONE;
                InvalidateRect(hwnd, nullptr, FALSE);
                return 0;
            }
            if (coloringOpen) {
                int it = coloringItemAt(x,y);
                if (it >= 0) selectColoring(it);
                coloringOpen = false; coloringHover = -1; pressed = H_NONE;
                InvalidateRect(hwnd, nullptr, FALSE);
                return 0;
            }
            if (galleryOpen) {
                int it = galleryItemAt(x,y);
                if (it >= 0) applyPreset(it);
                galleryOpen = false; galleryHover = -1; pressed = H_NONE;
                InvalidateRect(hwnd, nullptr, FALSE);
                return 0;
            }
            Hit h = hitTest(x,y);
            if (h == H_PALETTE_DD) { paletteOpen = true; paletteHover = -1; pressed = H_NONE; InvalidateRect(hwnd, nullptr, FALSE); return 0; }
            if (h == H_EDE) { coloringOpen = true; coloringHover = -1; pressed = H_NONE; InvalidateRect(hwnd, nullptr, FALSE); return 0; }
            if (h == H_GALLERY_DD) { galleryOpen = true; galleryHover = -1; pressed = H_NONE; InvalidateRect(hwnd, nullptr, FALSE); return 0; }
            pressed = h;
            if (h == H_GRADIENT) { gradientDown(x); return 0; }
            if (h == H_MAXTRACK) { SetCapture(hwnd); setMaxFromT((double)(x-rcMaxTrack.left)/(rcMaxTrack.right-rcMaxTrack.left), false); InvalidateRect(hwnd,nullptr,FALSE); return 0; }
            if (h == H_DENSTRACK) { SetCapture(hwnd); color_density=(float)std::round(std::clamp(10.0+190.0*(x-rcDensTrack.left)/(rcDensTrack.right-rcDensTrack.left),10.0,200.0)); nav->SetRedisplay(); InvalidateRect(hwnd,nullptr,FALSE); return 0; }
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
                case H_SAVE: showExportDialog(hwnd, nav.get()); break;
                case H_COPY: copyLocation(); break;
                case H_PASTE: pasteLocation(); break;
                case H_SS: ssOn = !ssOn; setMethodFlag(ColoringMethod::SUPER_SAMPLING, ssOn); break;
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
            int wd = GET_WHEEL_DELTA_WPARAM(wp);
            // Mouse-wheel over a slider fine-tunes it (color density in 0.1 steps,
            // max iterations ~2% per notch on its log scale).
            RECT dt = rcDensTrack; dt.top -= S(10); dt.bottom += S(10);
            if (inRect(dt, q.x, q.y)) {
                double v = color_density + (wd > 0 ? 0.1 : -0.1);
                color_density = (float)(std::round(std::clamp(v, 10.0, 200.0) * 10.0) / 10.0);
                nav->SetRedisplay(); keepLive(); InvalidateRect(hwnd, nullptr, FALSE);
                return 0;
            }
            RECT mt = rcMaxTrack; mt.top -= S(10); mt.bottom += S(10);
            if (inRect(mt, q.x, q.y)) {
                int nv = (int)std::round(maxIter * (wd > 0 ? 1.02 : 1.0 / 1.02));
                if (nv == maxIter) nv += (wd > 0 ? 1 : -1);   // always move at small values
                maxIter = std::clamp(nv, 100, 5000000);
                nav->SetMxit(maxIter); startRender();
                return 0;
            }
            if (inRect(viewRect(), q.x, q.y)) {
                POINT p = mapToRender(q.x, q.y);
                if (wd > 0) nav->ZoomIn(p.x,p.y); else nav->ZoomOut(p.x,p.y);
                keepLive();
            }
            return 0;
        }
        case WM_KEYDOWN:
            if (wp == 'R') { nav->Reset(); renderStart = std::chrono::steady_clock::now(); wasComputing = true; keepLive(); }
            else if (wp == VK_SPACE) startRender();
            else if (wp == 'S') saveImage();
            else if (wp == 'C') copyLocation();
            else if (wp == 'V') pasteLocation();
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
            if (memDC) { SelectObject(memDC, memOld); DeleteObject(memBmp); DeleteDC(memDC); memDC = nullptr; }
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
