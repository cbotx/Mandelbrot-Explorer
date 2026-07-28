// High-resolution PNG export: the strip renderer + the custom-drawn export
// dialog, split out of win32_main.cpp. Entry point: showExportDialog().
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#ifndef _WIN32_WINNT
#define _WIN32_WINNT 0x0A00
#endif
#include <windows.h>
#include <windowsx.h>
#include <commdlg.h>
#include <wincodec.h>
#include <gmp.h>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "gui_theme.h"
#include "gui_export.h"
#include "Image.h"
#include "mandel_navigator.h"
#include "mandel_perturbation.h"

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
    ensureSrgbLut();
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

    // Relief needs neighbouring output pixels for the slope; strips break vertical
    // neighbours, so accumulate a full-image height field and post-shade at the end.
    const bool relief = relief_on != 0;
    std::vector<float> htbuf;
    if (relief) htbuf.assign((size_t)W * H, std::numeric_limits<float>::quiet_NaN());
    // Analytic normal map / DE overlay: the engine fills a per-subpixel second
    // buffer (normal angle, or the pixel-normalized distance estimate).
    const bool normal = normal_light_on != 0;
    const bool deovl = de_overlay_on != 0;
    std::vector<float> normbuf, nfield;
    if (normal || deovl) {
        normbuf.assign((size_t)Wss * sh * ss, 0.0f);
        mandel.setNormalBuffer(normbuf.data());
        nfield.assign((size_t)W * H, std::numeric_limits<float>::quiet_NaN());
    }

    for (int s = 0; s < nstrips && !st->cancel; ++s) {
        int obase = s * sh;
        int oh = std::min(sh, H - obase);
        std::fill(ibuf.begin(), ibuf.end(), (float)EMPTYPIXEL);
        mandel.SetHalt(false);
        mandel.SetProgress(&st->progress, (float)s / nstrips, 0.95f / nstrips);
        mandel.Compute(cx, cy, scale_e, mxit, cmethod, Hss, obase * ss);
        mandel.SetProgress(nullptr);
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
                // Engine rows are y-up (row 0 = smallest imaginary = bottom); flip
                // vertically so the PNG/preview is upright like the main view.
                size_t orow = (size_t)(H - 1 - (obase + oi));
                uint8_t* p = &out[(orow * W + oj) * 3];
                p[0] = (uint8_t)(srgbEncode(lr / n) * 255.0 + 0.5);
                p[1] = (uint8_t)(srgbEncode(lg / n) * 255.0 + 0.5);
                p[2] = (uint8_t)(srgbEncode(lb / n) * 255.0 + 0.5);
                if (relief) {
                    float hv = ibuf[(size_t)(oi * ss + ss / 2) * Wss + (oj * ss + ss / 2)];
                    htbuf[orow * W + oj] =
                        (hv == EMPTYPIXEL || hv == INTERIOR_SENTINEL)
                            ? std::numeric_limits<float>::quiet_NaN() : (hv < 0 ? 0.0f : hv);
                }
                if (normal || deovl) {
                    size_t ci = (size_t)(oi * ss + ss / 2) * Wss + (oj * ss + ss / 2);
                    float hv = ibuf[ci];
                    nfield[orow * W + oj] =
                        (hv == EMPTYPIXEL || hv == INTERIOR_SENTINEL)
                            ? std::numeric_limits<float>::quiet_NaN() : normbuf[ci];
                }
            }
        }
        st->progress = (float)(s + 1) / nstrips;
    }
    if (relief && !st->cancel) applyReliefTo(out.data(), htbuf.data(), W, H);
    if (normal && !st->cancel) applyNormalLightTo(out.data(), nfield.data(), W, H);
    if (deovl && !st->cancel) applyDEOverlayTo(out.data(), nfield.data(), W, H);
    { std::lock_guard<std::mutex> lk(st->curMx); st->cur = nullptr; }
    mpf_clear(scale_e);
    st->ok = !st->cancel;
    st->done = true;
}

// Collision-free default save name: mandelbrot_<YYYYMMDD>_<HHMMSS>.png (local time).
std::wstring defaultSaveName() {
    SYSTEMTIME t; GetLocalTime(&t);
    wchar_t buf[64];
    swprintf_s(buf, L"mandelbrot_%04d%02d%02d_%02d%02d%02d.png",
        t.wYear, t.wMonth, t.wDay, t.wHour, t.wMinute, t.wSecond);
    return buf;
}

// Write W x H top-down RGB to a PNG via the Windows Imaging Component (a system
// DLL -- nothing to ship). Runs on the export thread, so it inits COM itself.
bool writeExportPNG(const wchar_t* path, const std::vector<uint8_t>& rgb, int W, int H) {
    bool comInit = SUCCEEDED(CoInitializeEx(nullptr, COINIT_APARTMENTTHREADED));
    bool ok = false;
    IWICImagingFactory* fac = nullptr;
    if (SUCCEEDED(CoCreateInstance(CLSID_WICImagingFactory, nullptr, CLSCTX_INPROC_SERVER, IID_PPV_ARGS(&fac)))) {
        IWICStream* stream = nullptr;
        if (SUCCEEDED(fac->CreateStream(&stream)) &&
            SUCCEEDED(stream->InitializeFromFilename(path, GENERIC_WRITE))) {
            IWICBitmapEncoder* enc = nullptr;
            if (SUCCEEDED(fac->CreateEncoder(GUID_ContainerFormatPng, nullptr, &enc)) &&
                SUCCEEDED(enc->Initialize(stream, WICBitmapEncoderNoCache))) {
                IWICBitmapFrameEncode* frame = nullptr; IPropertyBag2* props = nullptr;
                if (SUCCEEDED(enc->CreateNewFrame(&frame, &props)) && SUCCEEDED(frame->Initialize(props))) {
                    frame->SetSize((UINT)W, (UINT)H);
                    WICPixelFormatGUID fmt = GUID_WICPixelFormat24bppBGR;
                    frame->SetPixelFormat(&fmt);
                    std::vector<uint8_t> bgr((size_t)W * H * 3);
                    for (size_t i = 0; i < (size_t)W * H; ++i) { bgr[i*3]=rgb[i*3+2]; bgr[i*3+1]=rgb[i*3+1]; bgr[i*3+2]=rgb[i*3]; }
                    if (IsEqualGUID(fmt, GUID_WICPixelFormat24bppBGR) &&
                        SUCCEEDED(frame->WritePixels((UINT)H, (UINT)(W * 3), (UINT)bgr.size(), bgr.data())) &&
                        SUCCEEDED(frame->Commit()) && SUCCEEDED(enc->Commit())) ok = true;
                    if (props) props->Release();
                    frame->Release();
                }
                enc->Release();
            }
            stream->Release();
        }
        fac->Release();
    }
    if (comInit) CoUninitialize();
    return ok;
}

// ---- Export dialog (custom-drawn, matching the main panel's dark theme) ----
struct ResPreset { const wchar_t* name; int w, h; };

// self-drawn widget ids for the export dialog
enum ExHit { EX_NONE, EX_RES, EX_SS, EX_WFIELD, EX_HFIELD, EX_EXPORT, EX_CANCEL };

struct ExportDlg {
    HWND hwnd = nullptr, owner = nullptr;
    MandelNavigator* nav = nullptr;
    int cmethod = 0, mxit = 0;
    mpf_t cx, cy, scale;

    HFONT fUi = nullptr, fBold = nullptr, fSmall = nullptr;
    int dpi = 96;

    std::vector<ResPreset> presets;               // last entry is Custom (w=h=0)
    int resSel = 0, ssSel = 1;                     // selected resolution / SS index
    std::wstring wText = L"1920", hText = L"1280"; // custom W/H fields

    bool resOpen = false, ssOpen = false;          // dropdown list expanded
    int editField = 0;                             // 0 none, 1 W, 2 H
    int hover = EX_NONE, hoverItem = -1;           // hover highlight
    int caretTick = 0; bool caretOn = true;

    RECT rcRes{}, rcSS{}, rcW{}, rcH{}, rcProg{}, rcCancel{}, rcExport{}, rcPreview{};
    int cliH = 466;                               // client height (logical px), may shrink to fit

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
        if (fUi) DeleteObject(fUi); if (fBold) DeleteObject(fBold); if (fSmall) DeleteObject(fSmall);
        mpf_clear(cx); mpf_clear(cy); mpf_clear(scale);
    }
    int S(int v) const { return MulDiv(v, dpi, 96); }
    bool custom() const { return resSel >= 0 && resSel < (int)presets.size() && presets[resSel].w == 0; }
    void layout() {
        int H = cliH;
        rcRes     = { S(16),  S(32),  S(332), S(66) };
        rcSS      = { S(348), S(32),  S(528), S(66) };
        rcW       = { S(16),  S(98),  S(150), S(132) };
        rcH       = { S(178), S(98),  S(312), S(132) };
        int btnTop = S(H - 52), btnBot = S(H - 16);
        rcCancel  = { S(332), btnTop, S(424), btnBot };
        rcExport  = { S(436), btnTop, S(528), btnBot };
        int progTop = S(H - 68);
        rcProg    = { S(16), progTop, S(528), progTop + S(12) };
        rcPreview = { S(16), S(148),  S(528), progTop - S(14) };
    }
};

// Standard sizes (largest first); only those fitting the monitor are offered
// (plus Custom). Defaults the selection to the view-matching 1350x900 entry.
static void buildResPresets(ExportDlg* d) {
    static const ResPreset all[] = {
        { L"7680 x 4320  (8K)", 7680, 4320 },
        { L"5120 x 2880  (5K)", 5120, 2880 },
        { L"3840 x 2160  (4K)", 3840, 2160 },
        { L"2700 x 1800  (3:2)", 2700, 1800 },
        { L"2560 x 1440  (1440p)", 2560, 1440 },
        { L"1920 x 1080  (1080p)", 1920, 1080 },
        { L"1350 x 900  (view 3:2)", 1350, 900 },
    };
    int mw = GetSystemMetrics(SM_CXSCREEN), mh = GetSystemMetrics(SM_CYSCREEN);
    d->presets.clear();
    for (const auto& r : all) if (r.w <= mw && r.h <= mh) d->presets.push_back(r);
    if (d->presets.empty()) d->presets.push_back({ L"1350 x 900  (view 3:2)", 1350, 900 });
    d->presets.push_back({ L"Custom\u2026", 0, 0 });
    d->resSel = 0;
    for (int i = 0; i < (int)d->presets.size(); ++i)
        if (d->presets[i].w == 1350 && d->presets[i].h == 900) { d->resSel = i; break; }
}

static const int kSSFactors[] = { 1, 2, 3, 4, 6, 8 };
static const wchar_t* const kSSNames[] = { L"1\u00d7  (fast)", L"2\u00d7", L"3\u00d7", L"4\u00d7", L"6\u00d7", L"8\u00d7  (max)" };
static const int kSSCount = 6;

static void exportGetSize(ExportDlg* d) {
    if (d->custom()) {
        d->outW = std::clamp(_wtoi(d->wText.c_str()), 16, 30000);
        d->outH = std::clamp(_wtoi(d->hText.c_str()), 16, 30000);
    } else {
        d->outW = d->presets[d->resSel].w;
        d->outH = d->presets[d->resSel].h;
    }
    d->ss = kSSFactors[std::clamp(d->ssSel, 0, kSSCount - 1)];
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

// Commit the currently-edited custom W/H field (clamp + refresh preview).
static void exCommitEdit(ExportDlg* d) {
    if (!d->editField) return;
    std::wstring& t = (d->editField == 1) ? d->wText : d->hText;
    int v = std::clamp(_wtoi(t.c_str()), 16, 30000);
    t = std::to_wstring(v);
    d->editField = 0;
    startPreview(d);
}

// --- custom widget painters (dark theme, rounded cards) --------------------
static void exDrawDD(ExportDlg* d, HDC dc, RECT r, const std::wstring& text, bool hovered, bool open) {
    fillRound(dc, r, (hovered || open) ? CLR_CARD_HOV : CLR_CARD, open ? CLR_ACCENT : CLR_BORDER, d->S(8));
    RECT tr = r; tr.left += d->S(12); tr.right -= d->S(28);
    drawText(dc, tr, text, CLR_TEXT, d->fUi, DT_LEFT | DT_VCENTER | DT_SINGLELINE);
    int cx = r.right - d->S(18), cy = (r.top + r.bottom) / 2, w = d->S(4);
    POINT pts[3] = { { cx - w, cy - w / 2 }, { cx + w, cy - w / 2 }, { cx, cy + w / 2 } };
    HBRUSH b = CreateSolidBrush(CLR_TEXT_DIM); HPEN pen = CreatePen(PS_SOLID, 1, CLR_TEXT_DIM);
    HGDIOBJ ob = SelectObject(dc, b), op = SelectObject(dc, pen);
    Polygon(dc, pts, 3);
    SelectObject(dc, ob); SelectObject(dc, op); DeleteObject(b); DeleteObject(pen);
}
static void exDrawField(ExportDlg* d, HDC dc, RECT r, const std::wstring& text, bool focused, bool enabled) {
    fillRound(dc, r, enabled ? CLR_CARD : CLR_PANEL, focused ? CLR_ACCENT : CLR_BORDER, d->S(8));
    RECT tr = r; tr.left += d->S(12); tr.right -= d->S(12);
    drawText(dc, tr, text, enabled ? CLR_TEXT : CLR_TEXT_DIM, d->fUi, DT_LEFT | DT_VCENTER | DT_SINGLELINE);
    if (focused && enabled && d->caretOn) {
        HGDIOBJ of = SelectObject(dc, d->fUi);
        SIZE sz{}; GetTextExtentPoint32W(dc, text.c_str(), (int)text.size(), &sz);
        SelectObject(dc, of);
        int cx = tr.left + sz.cx + d->S(1), mid = (r.top + r.bottom) / 2;
        HPEN pen = CreatePen(PS_SOLID, std::max(1, d->S(1)), CLR_ACCENT); HGDIOBJ op = SelectObject(dc, pen);
        MoveToEx(dc, cx, mid - d->S(9), nullptr); LineTo(dc, cx, mid + d->S(9));
        SelectObject(dc, op); DeleteObject(pen);
    }
}
static void exDrawBtn(ExportDlg* d, HDC dc, RECT r, const std::wstring& text, bool hovered, bool accent, bool disabled) {
    COLORREF fill = disabled ? CLR_CARD : accent ? (hovered ? CLR_ACCENT_HI : CLR_ACCENT)
                                                  : (hovered ? CLR_CARD_HOV : CLR_CARD);
    fillRound(dc, r, fill, accent && !disabled ? fill : CLR_BORDER, d->S(8));
    COLORREF tc = disabled ? CLR_TEXT_DIM : (accent ? RGB(10, 14, 22) : CLR_TEXT);
    drawText(dc, r, text, tc, d->fBold, DT_CENTER | DT_VCENTER | DT_SINGLELINE);
}
static void exDrawProg(ExportDlg* d, HDC dc, RECT r, float frac) {
    fillRound(dc, r, CLR_TRACK, CLR_BORDER, d->S(5));
    frac = std::clamp(frac, 0.0f, 1.0f);
    if (frac > 0.0f) {
        RECT f = r; f.right = r.left + (int)((r.right - r.left) * frac);
        if (f.right < f.left + d->S(10)) f.right = std::min(r.right, f.left + d->S(10));
        fillRound(dc, f, CLR_ACCENT, CLR_ACCENT, d->S(5));
    }
}
static bool exPtIn(const RECT& r, POINT p) { return p.x >= r.left && p.x < r.right && p.y >= r.top && p.y < r.bottom; }
static RECT exItemRect(ExportDlg* d, const RECT& dd, int i) {
    int ih = d->S(30);
    return { dd.left, dd.bottom + i * ih, dd.right, dd.bottom + (i + 1) * ih };
}
static int exHitTest(ExportDlg* d, POINT p) {
    if (exPtIn(d->rcRes, p)) return EX_RES;
    if (exPtIn(d->rcSS, p)) return EX_SS;
    if (exPtIn(d->rcW, p)) return EX_WFIELD;
    if (exPtIn(d->rcH, p)) return EX_HFIELD;
    if (exPtIn(d->rcExport, p)) return EX_EXPORT;
    if (exPtIn(d->rcCancel, p)) return EX_CANCEL;
    return EX_NONE;
}

static UINT exWinDpi(HWND h) {
    UINT d = GetDpiForWindow(h); return d ? d : 96;
}

static void exMakeFonts(ExportDlg* d) {
    if (d->fUi) DeleteObject(d->fUi);
    if (d->fBold) DeleteObject(d->fBold);
    if (d->fSmall) DeleteObject(d->fSmall);
    auto mk = [&](int px, int weight) {
        return CreateFontW(-MulDiv(px, d->dpi, 96), 0, 0, 0, weight, FALSE, FALSE, FALSE,
            DEFAULT_CHARSET, OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS, CLEARTYPE_QUALITY,
            DEFAULT_PITCH | FF_DONTCARE, L"Segoe UI");
    };
    d->fUi = mk(15, FW_NORMAL); d->fBold = mk(15, FW_SEMIBOLD); d->fSmall = mk(12, FW_NORMAL);
}

// Resize the window frame so its client area exactly holds the layout at dpi.
// Size the window frame so its client area holds the layout at the current DPI,
// shrinking the (preview) height if the monitor work area can't fit the ideal.
static void exFitWindow(HWND h, ExportDlg* d, const RECT* pos) {
    const DWORD style = WS_POPUP | WS_CAPTION | WS_SYSMENU;
    // non-client frame thickness at this DPI
    RECT fr = { 0, 0, d->S(544), d->S(466) };
    AdjustWindowRectExForDpi(&fr, style, FALSE, 0, d->dpi);
    int frameH = (fr.bottom - fr.top) - d->S(466);
    // clamp client height to the monitor's work area
    HMONITOR mon = MonitorFromWindow(h, MONITOR_DEFAULTTONEAREST);
    MONITORINFO mi{ sizeof(mi) };
    int cliHpx = d->S(466);
    if (GetMonitorInfoW(mon, &mi)) {
        int maxClient = (mi.rcWork.bottom - mi.rcWork.top) - frameH - d->S(24);
        cliHpx = std::min(cliHpx, maxClient);
    }
    d->cliH = std::clamp(MulDiv(cliHpx, 96, d->dpi), 380, 466);
    d->layout();
    RECT rc = { 0, 0, d->S(544), d->S(d->cliH) };
    AdjustWindowRectExForDpi(&rc, style, FALSE, 0, d->dpi);
    int w = rc.right - rc.left, ht = rc.bottom - rc.top;
    UINT flags = SWP_NOZORDER | SWP_NOACTIVATE | (pos ? 0 : SWP_NOMOVE);
    SetWindowPos(h, nullptr, pos ? pos->left : 0, pos ? pos->top : 0, w, ht, flags);
}

LRESULT CALLBACK ExportWndProc(HWND h, UINT m, WPARAM wp, LPARAM lp);

void showExportDialog(HWND owner, MandelNavigator* nav) {
    static bool registered = false;
    HINSTANCE inst = (HINSTANCE)GetWindowLongPtrW(owner, GWLP_HINSTANCE);
    if (!registered) {
        WNDCLASSW wc{}; wc.lpfnWndProc = ExportWndProc; wc.hInstance = inst;
        wc.hCursor = LoadCursor(nullptr, IDC_ARROW);
        wc.hbrBackground = nullptr;
        wc.lpszClassName = L"MandelExportDlg";
        RegisterClassW(&wc); registered = true;
    }
    ExportDlg* d = new ExportDlg();
    d->owner = owner; d->nav = nav;
    d->cmethod = nav->GetCMethod(); d->mxit = nav->GetMxit();
    nav->GetView(d->cx, d->cy, d->scale);
    buildResPresets(d);
    UINT dpi = exWinDpi(owner);
    d->dpi = dpi;
    RECT rc = { 0, 0, MulDiv(544, dpi, 96), MulDiv(466, dpi, 96) };
    DWORD style = WS_POPUP | WS_CAPTION | WS_SYSMENU;
    AdjustWindowRectExForDpi(&rc, style, FALSE, 0, dpi);
    int W = rc.right - rc.left, H = rc.bottom - rc.top;
    RECT orc; GetWindowRect(owner, &orc);
    HWND hw = CreateWindowExW(0, L"MandelExportDlg", L"Export image",
        style | WS_VISIBLE, orc.left + 80, orc.top + 80, W, H, owner, nullptr, inst, d);
    d->hwnd = hw;
    // The window's effective DPI is only reliable once it exists; refit its frame
    // to the real DPI (WM_CREATE's SetWindowPos is unreliable on some systems).
    UINT rdpi = exWinDpi(hw);
    if (rdpi != dpi) { d->dpi = rdpi; exMakeFonts(d); d->layout(); }
    exFitWindow(hw, d, nullptr);
    EnableWindow(owner, FALSE);
}

static void exportDoSave(ExportDlg* d) {
    wchar_t path[MAX_PATH]; wcscpy_s(path, defaultSaveName().c_str());
    OPENFILENAMEW ofn{}; ofn.lStructSize = sizeof(ofn); ofn.hwndOwner = d->hwnd;
    ofn.lpstrFilter = L"PNG image (*.png)\0*.png\0"; ofn.lpstrFile = path; ofn.nMaxFile = MAX_PATH;
    ofn.lpstrDefExt = L"png"; ofn.Flags = OFN_OVERWRITEPROMPT | OFN_PATHMUSTEXIST;
    if (!GetSaveFileNameW(&ofn)) return;
    exportGetSize(d);
    d->exState.cancel = false; d->exState.done = false; d->exState.ok = false; d->exState.progress = 0;
    d->exporting = true;
    d->resOpen = d->ssOpen = false; d->editField = 0;
    std::wstring save = path;
    int W = d->outW, H = d->outH, ss = d->ss;
    if (d->exThread.joinable()) d->exThread.join();
    d->exThread = std::thread([d, W, H, ss, save]() {
        exportRender(d->cx, d->cy, d->scale, W, H, ss, d->mxit, d->cmethod, d->exImg, &d->exState);
        if (d->exState.ok) writeExportPNG(save.c_str(), d->exImg, W, H);
    });
    InvalidateRect(d->hwnd, nullptr, FALSE);
}

// Draws the whole dialog surface into dc (a cw x ch device area at d->dpi).
// Shared by WM_PAINT and the offscreen layout dump.
static void exPaintContent(ExportDlg* d, HDC dc, int cw, int ch) {
    fillRect(dc, { 0, 0, cw, ch }, CLR_BG);
    auto S = [&](int v) { return d->S(v); };
    drawText(dc, { S(16), S(12), S(332), S(30) }, L"Resolution", CLR_TEXT_DIM, d->fSmall, DT_LEFT | DT_VCENTER | DT_SINGLELINE);
    drawText(dc, { S(348), S(12), S(528), S(30) }, L"Supersampling", CLR_TEXT_DIM, d->fSmall, DT_LEFT | DT_VCENTER | DT_SINGLELINE);
    drawText(dc, { S(16), S(78), S(312), S(96) }, L"Custom size  (width \u00d7 height)", CLR_TEXT_DIM, d->fSmall, DT_LEFT | DT_VCENTER | DT_SINGLELINE);

    exDrawDD(d, dc, d->rcRes, d->presets[d->resSel].name, d->hover == EX_RES, d->resOpen);
    exDrawDD(d, dc, d->rcSS, kSSNames[std::clamp(d->ssSel, 0, kSSCount - 1)], d->hover == EX_SS, d->ssOpen);

    bool cust = d->custom();
    exDrawField(d, dc, d->rcW, cust ? d->wText : std::to_wstring(d->outW), d->editField == 1, cust);
    exDrawField(d, dc, d->rcH, cust ? d->hText : std::to_wstring(d->outH), d->editField == 2, cust);
    drawText(dc, { d->rcW.right, d->rcW.top, d->rcH.left, d->rcW.bottom }, L"\u00d7", CLR_TEXT_DIM, d->fUi, DT_CENTER | DT_VCENTER | DT_SINGLELINE);

    fillRound(dc, d->rcPreview, CLR_PANEL, CLR_BORDER, d->S(8));
    RECT pr = d->rcPreview;
    if (d->pvW > 0 && !d->preview.empty()) {
        int bw = (pr.right - pr.left) - d->S(8), bh = (pr.bottom - pr.top) - d->S(8);
        double s = std::min((double)bw / d->pvW, (double)bh / d->pvH);
        int dw = (int)(d->pvW * s), dh = (int)(d->pvH * s);
        int dx = pr.left + d->S(4) + (bw - dw) / 2, dy = pr.top + d->S(4) + (bh - dh) / 2;
        BITMAPINFO bi{}; bi.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
        bi.bmiHeader.biWidth = d->pvW; bi.bmiHeader.biHeight = -d->pvH;
        bi.bmiHeader.biPlanes = 1; bi.bmiHeader.biBitCount = 24; bi.bmiHeader.biCompression = BI_RGB;
        SetStretchBltMode(dc, HALFTONE);
        StretchDIBits(dc, dx, dy, dw, dh, 0, 0, d->pvW, d->pvH, d->preview.data(), &bi, DIB_RGB_COLORS, SRCCOPY);
    } else {
        drawText(dc, pr, L"Rendering preview\u2026", CLR_TEXT_DIM, d->fUi, DT_CENTER | DT_VCENTER | DT_SINGLELINE);
    }

    exDrawProg(d, dc, d->rcProg, d->exporting ? d->exState.progress.load() : 0.0f);

    wchar_t info[96];
    if (d->exporting) swprintf_s(info, L"Rendering\u2026 %d%%", (int)(d->exState.progress.load() * 100 + 0.5f));
    else swprintf_s(info, L"%d \u00d7 %d px    \u00b7    %d\u00d7 supersampling", d->outW, d->outH, d->ss);
    RECT infoR = { S(16), d->rcCancel.top, d->rcCancel.left - S(12), d->rcCancel.bottom };
    drawText(dc, infoR, info, CLR_TEXT_DIM, d->fSmall, DT_LEFT | DT_VCENTER | DT_SINGLELINE);

    exDrawBtn(d, dc, d->rcCancel, d->exporting ? L"Cancel" : L"Close", d->hover == EX_CANCEL, false, false);
    exDrawBtn(d, dc, d->rcExport, L"Export\u2026", d->hover == EX_EXPORT, true, d->exporting);

    auto drawList = [&](const RECT& dd, int count, int sel, auto nameOf) {
        for (int i = 0; i < count; ++i) {
            RECT ir = exItemRect(d, dd, i);
            fillRect(dc, ir, d->hoverItem == i ? CLR_CARD_HOV : CLR_CARD);
            RECT tr = ir; tr.left += d->S(12);
            drawText(dc, tr, nameOf(i), i == sel ? CLR_ACCENT_HI : CLR_TEXT, d->fUi, DT_LEFT | DT_VCENTER | DT_SINGLELINE);
        }
        RECT box = { dd.left, dd.bottom, dd.right, exItemRect(d, dd, count - 1).bottom };
        HPEN pen = CreatePen(PS_SOLID, 1, CLR_BORDER);
        HGDIOBJ op = SelectObject(dc, pen), ob2 = SelectObject(dc, GetStockObject(NULL_BRUSH));
        Rectangle(dc, box.left, box.top, box.right, box.bottom);
        SelectObject(dc, op); SelectObject(dc, ob2); DeleteObject(pen);
    };
    if (d->resOpen) drawList(d->rcRes, (int)d->presets.size(), d->resSel, [&](int i) { return std::wstring(d->presets[i].name); });
    else if (d->ssOpen) drawList(d->rcSS, kSSCount, d->ssSel, [&](int i) { return std::wstring(kSSNames[i]); });
}

LRESULT CALLBACK ExportWndProc(HWND h, UINT m, WPARAM wp, LPARAM lp) {
    ExportDlg* d = (ExportDlg*)GetWindowLongPtrW(h, GWLP_USERDATA);
    switch (m) {
    case WM_CREATE: {
        CREATESTRUCTW* cs = (CREATESTRUCTW*)lp;
        d = (ExportDlg*)cs->lpCreateParams;
        SetWindowLongPtrW(h, GWLP_USERDATA, (LONG_PTR)d);
        d->dpi = exWinDpi(h);
        exMakeFonts(d);
        d->layout();
        SetTimer(h, 1, 100, nullptr);
        startPreview(d);
        return 0;
    }
    case WM_DPICHANGED: {
        d->dpi = HIWORD(wp);
        exMakeFonts(d);
        d->layout();
        RECT* pr = (RECT*)lp;
        exFitWindow(h, d, pr);
        InvalidateRect(h, nullptr, FALSE);
        return 0;
    }
    case WM_ERASEBKGND:
        return 1;
    case WM_MOUSEMOVE: {
        POINT p = { GET_X_LPARAM(lp), GET_Y_LPARAM(lp) };
        int nh = EX_NONE, ni = -1;
        if (d->resOpen) { for (int i = 0; i < (int)d->presets.size(); ++i) if (exPtIn(exItemRect(d, d->rcRes, i), p)) ni = i; }
        else if (d->ssOpen) { for (int i = 0; i < kSSCount; ++i) if (exPtIn(exItemRect(d, d->rcSS, i), p)) ni = i; }
        else nh = exHitTest(d, p);
        if (nh != d->hover || ni != d->hoverItem) { d->hover = nh; d->hoverItem = ni; InvalidateRect(h, nullptr, FALSE); }
        return 0;
    }
    case WM_LBUTTONDOWN: {
        POINT p = { GET_X_LPARAM(lp), GET_Y_LPARAM(lp) };
        SetFocus(h);
        if (d->resOpen) {
            for (int i = 0; i < (int)d->presets.size(); ++i)
                if (exPtIn(exItemRect(d, d->rcRes, i), p)) {
                    d->resSel = i; d->resOpen = false; startPreview(d); InvalidateRect(h, nullptr, FALSE); return 0;
                }
            d->resOpen = false; InvalidateRect(h, nullptr, FALSE); return 0;
        }
        if (d->ssOpen) {
            for (int i = 0; i < kSSCount; ++i)
                if (exPtIn(exItemRect(d, d->rcSS, i), p)) {
                    d->ssSel = i; d->ssOpen = false; startPreview(d); InvalidateRect(h, nullptr, FALSE); return 0;
                }
            d->ssOpen = false; InvalidateRect(h, nullptr, FALSE); return 0;
        }
        int hit = exHitTest(d, p);
        if (d->editField && hit != EX_WFIELD && hit != EX_HFIELD) exCommitEdit(d);
        if (d->exporting && hit != EX_CANCEL) return 0;
        switch (hit) {
        case EX_RES: d->resOpen = true; d->ssOpen = false; d->editField = 0; break;
        case EX_SS:  d->ssOpen = true; d->resOpen = false; d->editField = 0; break;
        case EX_WFIELD: if (d->custom()) { d->editField = 1; d->caretOn = true; } break;
        case EX_HFIELD: if (d->custom()) { d->editField = 2; d->caretOn = true; } break;
        case EX_EXPORT: if (!d->exporting) exportDoSave(d); break;
        case EX_CANCEL:
            if (d->exporting) { d->exState.cancel = true; std::lock_guard<std::mutex> lk(d->exState.curMx); Mandel* mm = d->exState.cur.load(); if (mm) mm->SetHalt(true); }
            else { DestroyWindow(h); return 0; }
            break;
        default: break;
        }
        InvalidateRect(h, nullptr, FALSE);
        return 0;
    }
    case WM_KEYDOWN:
        if (wp == VK_ESCAPE) {
            if (d->resOpen || d->ssOpen) { d->resOpen = d->ssOpen = false; InvalidateRect(h, nullptr, FALSE); }
            else if (d->editField) { d->editField = 0; InvalidateRect(h, nullptr, FALSE); }
            else if (!d->exporting) DestroyWindow(h);
        }
        return 0;
    case WM_CHAR:
        if (d->editField && d->custom()) {
            std::wstring& t = (d->editField == 1) ? d->wText : d->hText;
            wchar_t ch = (wchar_t)wp;
            if (ch >= L'0' && ch <= L'9') { if (t.size() < 5) t.push_back(ch); }
            else if (ch == L'\b') { if (!t.empty()) t.pop_back(); }
            else if (ch == L'\r') exCommitEdit(d);
            InvalidateRect(h, nullptr, FALSE);
        }
        return 0;
    case WM_TIMER: {
        if (++d->caretTick >= 5) {
            d->caretTick = 0; d->caretOn = !d->caretOn;
            if (d->editField) { RECT* fr = (d->editField == 1) ? &d->rcW : &d->rcH; InvalidateRect(h, fr, FALSE); }
        }
        if (d->exporting) {
            InvalidateRect(h, nullptr, FALSE);
            if (d->exState.done) {
                if (d->exThread.joinable()) d->exThread.join();
                bool ok = d->exState.ok;
                d->exporting = false;
                InvalidateRect(h, nullptr, FALSE);
                if (ok) MessageBoxW(h, L"Image saved.", L"Export", MB_OK | MB_ICONINFORMATION);
            }
        }
        return 0;
    }
    case WM_PAINT: {
        PAINTSTRUCT ps; HDC wdc = BeginPaint(h, &ps);
        RECT cr; GetClientRect(h, &cr);
        HDC dc = CreateCompatibleDC(wdc);
        HBITMAP bm = CreateCompatibleBitmap(wdc, cr.right, cr.bottom);
        HGDIOBJ obm = SelectObject(dc, bm);
        exPaintContent(d, dc, cr.right, cr.bottom);
        BitBlt(wdc, 0, 0, cr.right, cr.bottom, dc, 0, 0, SRCCOPY);
        SelectObject(dc, obm); DeleteObject(bm); DeleteDC(dc);
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
