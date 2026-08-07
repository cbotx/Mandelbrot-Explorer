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
#include <wincodec.h>
#include <timeapi.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <limits>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "Image.h"
#include "color.h"
#include "gui_theme.h"
#include "gui_export.h"
#include "formula_dialog.h"
#include "formula_editor_panel.h"
#include "interpolate.h"
#include "mandel_navigator.h"
#include "mandel_perturbation.h"
#include "orbit_overlay.h"
#include "orbit_thumbnail.h"

#pragma comment(lib, "comdlg32.lib")
#pragma comment(lib, "gdi32.lib")
#pragma comment(lib, "user32.lib")
#pragma comment(lib, "winmm.lib")


namespace {

// Fixed render geometry + theme + self-drawn widget helpers live in gui_theme.h.
constexpr UINT WM_OPEN_FORMULA_TEST = WM_APP + 17;



// hit-test ids for self-drawn widgets
enum Hit {
    H_NONE, H_VIEW, H_GRADIENT,
    H_RESET, H_RENDER, H_SAVE, H_COPY, H_PASTE,
    H_MAXTRACK, H_DENSTRACK, H_PHASETRACK, H_SPEEDTRACK,
    H_RELAZ, H_RELEL, H_RELSTR,
    H_SS, H_FORMULA, H_JULIA, H_ORBIT, H_EDE, H_PALETTE_DD, H_COLOR, H_GALLERY_DD, H_PANELSB
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
        // --- rainbow family (Lab-interpolated -> smooth, non-muddy transitions) ---
        // Balanced jewel-toned spectrum; hues spaced evenly, wraps violet -> red.
        { L"Rainbow",  { {0.00f,205,60,72},{0.15f,232,132,52},{0.29f,240,206,92},{0.44f,110,190,110},{0.58f,58,176,184},{0.72f,72,112,202},{0.86f,150,92,192} } },
        // Higher chroma but kept in gamut, so colours read clean rather than neon-dirty.
        { L"Vibrant",  { {0.00f,232,38,58},{0.15f,255,138,20},{0.29f,250,222,40},{0.44f,64,198,82},{0.58f,26,188,200},{0.72f,44,96,228},{0.86f,166,58,214} } },
        // Pastel: each hue lifted toward a bright neutral grey (soft, premium, airy).
        { L"Pastel",   { {0.00f,232,168,172},{0.15f,242,202,166},{0.29f,244,228,182},{0.44f,186,222,190},{0.58f,168,218,222},{0.72f,176,192,236},{0.86f,206,182,230} } },
        // Rainbow with one bright peak + one deep-violet trough per cycle (a b&w value
        // sweep of the same period overlaid on the spectrum) -> strong depth/contrast.
        { L"Prism",    { {0.00f,28,10,44},{0.12f,92,58,168},{0.26f,66,120,212},{0.40f,74,196,196},{0.50f,238,246,236},{0.62f,246,202,86},{0.76f,232,112,58},{0.88f,150,44,96} } },
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
        std::vector<std::array<float, 3>> lab(n);
        for (int i = 0; i < n; ++i)
            rgb2lab(stops[i].r, stops[i].g, stops[i].b, lab[i][0], lab[i][1], lab[i][2]);
        // Tile the stops over three periods [-1, 0, +1] and sample the middle one,
        // so the monotone-cubic spline is truly periodic (continuous + C1) across
        // the color-cycle wrap (index colP-1 -> 0). A single ghost on each side
        // (the previous scheme) left a visible seam for presets whose stops don't
        // span [0,1] -- e.g. Gems (0.05..0.92) jumped ~11 dLab at the wrap.
        int m = 3 * n;
        std::vector<float> xs(m), out(colP);
        for (int p = -1; p <= 1; ++p)
            for (int i = 0; i < n; ++i)
                xs[(p + 1) * n + i] = stops[i].pos + p;
        for (int c = 0; c < 3; ++c) {
            std::vector<float> ys(m);
            for (int p = 0; p < 3; ++p)
                for (int i = 0; i < n; ++i) ys[p * n + i] = lab[i][c];
            mono_cubic_interpolate(xs.data(), ys.data(), m, out.data(), colP);
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

// Keep x/y on one line each so mode, coordinates and zoom always fit the card.
static std::wstring truncLocation(const std::string& raw, int cpl) {
    if (cpl < 8) cpl = 8;
    std::wstring out; std::string acc;
    auto flush = [&](bool last) {
        std::string L = acc; acc.clear();
        if ((L.rfind("x:", 0) == 0 || L.rfind("y:", 0) == 0) && (int)L.size() > cpl) {
            L = L.substr(0, std::max(0, cpl - 3)) + "...";
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
class App {
public:
    HWND hwnd = nullptr;
    HFONT fUi = nullptr, fBold = nullptr, fSmall = nullptr;
    HFONT fMono = nullptr, fMicro = nullptr;

    std::unique_ptr<MandelNavigator> nav;
    std::vector<uint8_t> bitmap = std::vector<uint8_t>((size_t)RENDER_W * RENDER_H * 3, 0);
    std::vector<uint8_t> display = std::vector<uint8_t>((size_t)RENDER_W * RENDER_H * 3, 0); // BGR
    std::vector<uint8_t> scaled;   // view-resolution upscale of `display` (BGR top-down)
    int scaledW = 0, scaledH = 0, scaledStride = 0;   // scaledStride = DWORD-aligned row bytes
    // Native-resolution rendering: the fractal is computed at the current view size
    // (renderW x renderH), so the on-screen blit is 1:1 (no upscale -> no beating).
    // These start at the compile-time default and follow the window on resize.
    int renderW = RENDER_W, renderH = RENDER_H;
    bool inSizeMove = false;       // between WM_ENTERSIZEMOVE/EXITSIZEMOVE (defer retarget)
    PaletteEditor palette;

    // widget rects (computed in layout())
    RECT rcReset{}, rcRender{}, rcSave{}, rcCopy{}, rcPaste{};
    RECT rcLocation{}, rcMaxTrack{}, rcDensTrack{}, rcPhaseTrack{}, rcSpeedTrack{};
    RECT rcReliefAz{}, rcReliefEl{}, rcReliefStr{};
    RECT rcSS{}, rcFormula{}, rcJulia{}, rcColoringDD{}, rcPaletteDD{}, rcColor{}, rcGradient{};
    RECT rcGalleryDD{};
    RECT rcOrbitToggle{}, rcOrbitThumb{};

    // state
    int maxIter = 500000;
    int savedMaxIter = 500000;
    bool ssOn = false;
    int coloringIdx = 0;                    // 0 Smooth, 1 Distance (EDE), 2 Feather (SAC)
    bool coloringOpen = false;              // coloring dropdown expanded
    int coloringHover = -1;                 // hovered item while open (-1 none)
    int paletteIdx = 0;                    // index into palettePresets()
    bool orbitOn = false;
    bool juliaUiEnabled = false;           // internal-only until parameters/perf are production ready
    FormulaDialogConfig formulaConfig;
    std::unique_ptr<FormulaEditorPanel> formulaEditor;
    int formulaEditorExpansionDip = 0;
    int formulaEditorHeightExpansionDip = 0;
    int savedColoringIdx = 0;
    bool savedSsOn = false, savedOrbitOn = false;
    bool hasSavedMandelUi = false;
    std::unique_ptr<OrbitWorker> orbitWorker;
    OrbitResult orbitResult;
    std::unique_ptr<OrbitThumbnailWorker> orbitThumbnailWorker;
    std::vector<uint8_t> orbitThumbnail;
    bool orbitThumbnailBuilding = false;
    uint64_t orbitThumbnailGeneration = 0;
    bool orbitThumbnailSmoke = false;
    int orbitThumbnailSmokeStage = 0;
    uint64_t orbitThumbnailSmokeFirstGeneration = 0;
    uint64_t orbitThumbnailSmokeExpectedGeneration = 0;
    std::chrono::steady_clock::time_point orbitThumbnailSmokeDeadline{};
    std::chrono::steady_clock::time_point lastOrbitRequest{};
    int lastOrbitX = -10000, lastOrbitY = -10000;
    bool paletteOpen = false;              // dropdown expanded
    int paletteHover = -1;                 // hovered item while open (-1 none)
    bool galleryOpen = false;              // gallery dropdown expanded
    int galleryHover = -1;                 // hovered gallery item (-1 none)
    bool navDragging = false, wasComputing = false;
    bool gradientFocused = false;
    Hit hover = H_NONE, pressed = H_NONE;
    int liveFrames = 0;   // frames still needing per-tick repaint (animation/compute)
    // ---- panel vertical scroll (the control stack can exceed the window height) --
    int panelScroll = 0;      // current scroll offset in device px (>=0)
    int panelContentH = 0;    // natural stacked height of all panel controls
    int panelViewH = 0;       // visible panel height (client bottom)
    bool panelSbDrag = false; int panelSbGrabY = 0, panelSbGrabScroll = 0;
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
    double lastRecolorMs = 0;   // last RecolorPhase/UpdateBitmap cost (animation diag)
    double animFps = 0;         // smoothed animation frame rate
    std::chrono::steady_clock::time_point lastFrameTs{};
    std::chrono::steady_clock::time_point lastRenderTs{};  // for the ~60fps redraw cap
    bool benchMode = false;     // MANDEL_GUI_BENCH: headless recolor micro-benchmark
    bool orbitBench = false;    // MANDEL_GUI_ORBIT_BENCH: compute one centre orbit and exit
    // MANDEL_GUI_FPSLOG=<sec>: run the REAL animation/timer loop (not the isolated
    // micro-bench) and log the measured sustained animFps + per-frame costs, then exit.
    bool fpsLogMode = false;
    double fpsLogDur = 0.0;
    std::chrono::steady_clock::time_point fpsLogStart{}, fpsLogLast{};
    int fpsLogSamples = 0; double fpsLogSum = 0, fpsLogMin = 1e9;

    // ---- palette-phase animation --------------------------------------------
    // animSpeed is the signed cycle rate in [-1, 1]; 0 (slider centre) = paused.
    // At |1| the palette cycles once every PHASE_SECONDS_AT_MAX seconds. Phase is
    // advanced by wall-clock time so the speed is frame-rate independent.
    float animSpeed = 0.0f;
    static constexpr double PHASE_SECONDS_AT_MAX = 2.0;
    std::chrono::steady_clock::time_point lastAnimTick = std::chrono::steady_clock::now();

    int S(int v) const { return MulDiv(v, dpi, 96); }   // scale a design px to device px
    void keepLive(int frames = 60) { liveFrames = std::max(liveFrames, frames); }

    static constexpr int ORBIT_W = 280;
    static constexpr int ORBIT_H = 192;
    static constexpr double ORBIT_X0 = -2.5;
    static constexpr double ORBIT_X1 = 1.0;
    static constexpr double ORBIT_Y0 = -1.2;
    static constexpr double ORBIT_Y1 = 1.2;

    void rebuildOrbitThumbnail() {
        orbitThumbnail.clear();
        orbitThumbnailBuilding = false;
        orbitThumbnailGeneration = 0;
        if (!orbitThumbnailWorker) return;
        if (!orbitOn || !nav) {
            orbitThumbnailWorker->cancel();
            return;
        }
        std::shared_ptr<const formula::ExpressionOrbitSnapshot> expression;
        if (nav->IsExpression()) {
            expression = nav->GetExpressionOrbitSnapshot();
            if (!expression) {
                orbitThumbnailWorker->cancel();
                return;
            }
        }
        orbitThumbnailGeneration = orbitThumbnailWorker->request(
            ORBIT_W, ORBIT_H,
            ORBIT_X0, ORBIT_X1, ORBIT_Y0, ORBIT_Y1,
            160, std::move(expression));
        orbitThumbnailBuilding = true;
        needFull = true;
        if (hwnd) {
            RECT dirty = rcOrbitThumb;
            dirty.top -= S(20);
            InvalidateRect(hwnd, &dirty, FALSE);
        }
    }

    POINT orbitMap(double re, double im, const RECT& r) const {
        double x = r.left + (re - ORBIT_X0) / (ORBIT_X1 - ORBIT_X0) *
                              (r.right - r.left - 1);
        double y = r.top + (ORBIT_Y1 - im) / (ORBIT_Y1 - ORBIT_Y0) *
                             (r.bottom - r.top - 1);
        constexpr double coordinateLimit = 10000000.0;
        return {
            (LONG)std::clamp(x, -coordinateLimit, coordinateLimit),
            (LONG)std::clamp(y, -coordinateLimit, coordinateLimit)
        };
    }

    void drawOrbitThumbnail(HDC dc) {
        if (!orbitOn || rcOrbitThumb.right <= rcOrbitThumb.left) return;
        const wchar_t* thumbnailLabel = L"Orbit in the Mandelbrot set";
        if (nav && nav->IsExpression()) {
            auto expression = nav->GetExpressionOrbitSnapshot();
            FormulaParameter pixelParameter = expression
                ? expression->pixelParameter : formulaConfig.pixelParameter;
            thumbnailLabel = pixelParameter == FormulaParameter::InitialZ
                ? L"Custom z0 initial-value plane"
                : L"Custom c parameter plane";
        }
        label(dc, rcOrbitThumb.left, rcOrbitThumb.top - S(20),
              thumbnailLabel);
        fillRound(dc, rcOrbitThumb, CLR_CARD, CLR_BORDER, S(8));
        RECT image = rcOrbitThumb;
        InflateRect(&image, -S(6), -S(6));
        const bool hasBackground =
            orbitThumbnail.size() == (size_t)ORBIT_W * ORBIT_H * 3;
        if (hasBackground) {
            BITMAPINFO bi{}; bi.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
            bi.bmiHeader.biWidth = ORBIT_W; bi.bmiHeader.biHeight = -ORBIT_H;
            bi.bmiHeader.biPlanes = 1; bi.bmiHeader.biBitCount = 24;
            bi.bmiHeader.biCompression = BI_RGB;
            SetStretchBltMode(dc, HALFTONE);
            StretchDIBits(dc, image.left, image.top,
                          image.right - image.left,
                          image.bottom - image.top,
                          0, 0, ORBIT_W, ORBIT_H,
                          orbitThumbnail.data(), &bi,
                          DIB_RGB_COLORS, SRCCOPY);
        } else {
            fillRound(dc, image, RGB(16, 13, 10), RGB(16, 13, 10), S(4));
        }

        int saved = SaveDC(dc);
        IntersectClipRect(dc, image.left, image.top, image.right, image.bottom);
        if (hasBackground && orbitResult.points.size() >= 2) {
            std::vector<POINT> points;
            points.reserve(orbitResult.points.size());
            for (const OrbitPoint& p : orbitResult.points) {
                if (!std::isfinite(p.re) || !std::isfinite(p.im))
                    break;
                points.push_back(orbitMap(p.re, p.im, image));
            }
            HPEN line = CreatePen(PS_SOLID, std::max(1, S(1)), RGB(255, 230, 96));
            HGDIOBJ old = SelectObject(dc, line);
            if (points.size() >= 2)
                Polyline(dc, points.data(), (int)points.size());
            SelectObject(dc, old); DeleteObject(line);
        }
        if (hasBackground && orbitResult.generation != 0) {
            POINT cp = orbitMap(
                orbitResult.pixelRe, orbitResult.pixelIm, image);
            HBRUSH marker = CreateSolidBrush(CLR_GREEN);
            HGDIOBJ oldBrush = SelectObject(dc, marker);
            HGDIOBJ oldPen = SelectObject(dc, GetStockObject(NULL_PEN));
            int mr = std::max(2, S(3));
            Ellipse(dc, cp.x - mr, cp.y - mr,
                    cp.x + mr + 1, cp.y + mr + 1);
            SelectObject(dc, oldPen);
            SelectObject(dc, oldBrush);
            DeleteObject(marker);
        }
        RestoreDC(dc, saved);

        if (!hasBackground) {
            drawText(dc, image,
                     orbitThumbnailBuilding
                         ? L"Building orbit background..."
                         : L"Background unavailable",
                     CLR_TEXT, fSmall,
                     DT_CENTER | DT_VCENTER | DT_SINGLELINE);
        } else if (orbitResult.generation == 0) {
            drawText(dc, image, L"Move over the main image", CLR_TEXT, fSmall,
                     DT_CENTER | DT_VCENTER | DT_SINGLELINE);
        } else {
            wchar_t text[96];
            swprintf_s(text, L"%d iter%s   %.1f ms", orbitResult.iterations,
                       orbitResult.escaped ? L"  escaped" : L"", orbitResult.computeMs);
            RECT info = image; info.left += S(6); info.bottom -= S(4);
            drawText(dc, info, text, CLR_TEXT, fSmall,
                     DT_LEFT | DT_BOTTOM | DT_SINGLELINE);
        }
    }

    void requestOrbit(int x, int y, bool applyBoundary = true) {
        if (!orbitOn || navDragging || !nav || !nav->SupportsOrbitOverlay()) return;
        std::shared_ptr<const formula::ExpressionOrbitSnapshot> expression;
        if (nav->IsExpression()) {
            expression = nav->GetExpressionOrbitSnapshot();
            if (!expression) return;
        }
        if (!orbitWorker) orbitWorker = std::make_unique<OrbitWorker>();
        auto now = std::chrono::steady_clock::now();
        if (lastOrbitRequest.time_since_epoch().count()) {
            double ms = std::chrono::duration<double, std::milli>(
                now - lastOrbitRequest).count();
            if (ms < 12.0 && std::abs(x - lastOrbitX) < 2 && std::abs(y - lastOrbitY) < 2)
                return;
        }
        POINT p = mapToRender(x, y);
        if (applyBoundary) {
            auto boundary = nav->GetBoundary();
            double fx = boundary[0][0] + (boundary[2][0] - boundary[0][0]) *
                        (p.x + 0.5) / renderW;
            double fy = boundary[0][1] + (boundary[2][1] - boundary[0][1]) *
                        ((renderH - 1 - p.y) + 0.5) / renderH;
            int sourceX = std::clamp((int)(fx * renderW), 0, renderW - 1);
            int sourceYUp = std::clamp((int)(fy * renderH), 0, renderH - 1);
            p = { sourceX, renderH - 1 - sourceYUp };
        }
        mp_bitcnt_t precision = nav->GetViewPrecision();
        mpf_t re, im, scale;
        mpf_init2(re, precision); mpf_init2(im, precision); mpf_init2(scale, precision);
        nav->GetView(re, im, scale);
        orbitWorker->request(re, im, scale, p.x, p.y, renderW, renderH, maxIter,
                             nav->GetFormulaContext(), std::move(expression),
                             nav->GetCustomDeepZoomPlan());
        mpf_clear(re); mpf_clear(im); mpf_clear(scale);
        lastOrbitRequest = now; lastOrbitX = x; lastOrbitY = y;
    }

    void requestCurrentOrbit() {
        if (!orbitOn || !nav || !nav->SupportsOrbitOverlay()) return;
        RECT vr = viewRect();
        int x = inRect(vr, lastOrbitX, lastOrbitY)
            ? lastOrbitX : (vr.left + vr.right) / 2;
        int y = inRect(vr, lastOrbitX, lastOrbitY)
            ? lastOrbitY : (vr.top + vr.bottom) / 2;
        lastOrbitRequest = {};
        requestOrbit(x, y);
    }

    void finishOrbitThumbnailSmoke(bool passed, const char* detail) {
        FILE* file = nullptr;
        fopen_s(&file, "build\\orbit_thumbnail_smoke.txt", "wb");
        if (file) {
            fprintf(file, "%s: %s\n", passed ? "PASS" : "FAIL", detail);
            fclose(file);
        }
        orbitThumbnailSmoke = false;
        PostQuitMessage(passed ? 0 : 1);
    }

    void startOrbitThumbnailSmoke() {
        FormulaDialogConfig first = formulaConfig;
        first.source = "sin(z)+c+p0";
        first.pixelParameter = FormulaParameter::C;
        first.parameters[0] = { 0.125, -0.0625 };
        first.bailout = 100.0;
        if (!applyFormulaConfig(first, false)) {
            finishOrbitThumbnailSmoke(false, "initial Custom formula rejected");
            return;
        }
        orbitThumbnailSmokeFirstGeneration = orbitThumbnailGeneration;

        FormulaDialogConfig latest = first;
        latest.source = "z*z+c+p0";
        latest.pixelParameter = FormulaParameter::InitialZ;
        latest.fixedC = { -0.25, 0.0 };
        latest.parameters[0] = { 0.015625, 0.0 };
        latest.bailout = 8.0;
        if (!applyFormulaConfig(latest, false)) {
            finishOrbitThumbnailSmoke(false, "latest z0 formula rejected");
            return;
        }
        orbitThumbnailSmokeExpectedGeneration = orbitThumbnailGeneration;
        if (orbitThumbnailSmokeExpectedGeneration <=
            orbitThumbnailSmokeFirstGeneration) {
            finishOrbitThumbnailSmoke(false, "generation did not advance");
            return;
        }
        orbitThumbnailSmokeStage = 0;
        orbitThumbnailSmokeDeadline =
            std::chrono::steady_clock::now() + std::chrono::seconds(20);
    }

    void checkOrbitThumbnailSmoke(const OrbitThumbnailResult& result) {
        if (!orbitThumbnailSmoke) return;
        const bool complete =
            orbitThumbnail.size() ==
                (size_t)ORBIT_W * ORBIT_H * 3 &&
            result.width == ORBIT_W && result.height == ORBIT_H;
        if (orbitThumbnailSmokeStage == 0) {
            if (!complete || !result.expression ||
                result.pixelParameter != FormulaParameter::InitialZ ||
                result.generation !=
                    orbitThumbnailSmokeExpectedGeneration) {
                finishOrbitThumbnailSmoke(
                    false, "latest Custom z0 background mismatch");
                return;
            }
            orbitThumbnailSmokeStage = 1;
            restoreMandelbrotUi(false);
            orbitThumbnailSmokeExpectedGeneration =
                orbitThumbnailGeneration;
            orbitThumbnailSmokeDeadline =
                std::chrono::steady_clock::now() +
                std::chrono::seconds(20);
            return;
        }
        if (!complete || result.expression ||
            result.generation != orbitThumbnailSmokeExpectedGeneration) {
            finishOrbitThumbnailSmoke(
                false, "Mandelbrot background restore mismatch");
            return;
        }
        finishOrbitThumbnailSmoke(
            true, "Custom latest-generation z0 and Mandel restore completed");
    }

    // Set the palette phase from a phase-slider x (0..colP over the track).
    void setPhaseFromX(int x) {
        double t = (double)(x - rcPhaseTrack.left) / std::max(1L, rcPhaseTrack.right - rcPhaseTrack.left);
        color_phase = (float)(std::clamp(t, 0.0, 1.0) * colP);
        recolorPhaseNow();
    }
    // Set the animation speed from a speed-slider x; snaps to the centre (0).
    void setSpeedFromX(int x) {
        double t = (double)(x - rcSpeedTrack.left) / std::max(1L, rcSpeedTrack.right - rcSpeedTrack.left);
        t = std::clamp(t, 0.0, 1.0);
        const double snaps[] = { 0.5 };
        t = snapT(t, snaps, 1, (double)(rcSpeedTrack.right - rcSpeedTrack.left));
        animSpeed = (float)((t - 0.5) * 2.0);
        lastAnimTick = std::chrono::steady_clock::now();   // avoid a phase jump on (re)start
    }
    // Relief light/strength sliders (post-shade only; no fractal recompute).
    static constexpr double RELIEF_STR_MAX = 3.0;
    // DE-overlay slider ranges: in DE mode the azimuth+strength slots become the two
    // DE tunables (boundary scale de_k and optical falloff scale de_scale).
    static constexpr double DE_K_MIN = 0.1, DE_K_MAX = 8.0;
    static constexpr double DE_SCALE_MIN = 0.1, DE_SCALE_MAX = 8.0;
    void setReliefAzFromX(int x) {
        double t = std::clamp((double)(x - rcReliefAz.left) / std::max(1L, rcReliefAz.right - rcReliefAz.left), 0.0, 1.0);
        if (de_overlay_on) de_k = (float)(DE_K_MIN + t * (DE_K_MAX - DE_K_MIN));
        else relief_light_az = (float)(t * 2.0 * 3.14159265358979);
        recolorPhaseNow();
    }
    void setReliefElFromX(int x) {
        double t = std::clamp((double)(x - rcReliefEl.left) / std::max(1L, rcReliefEl.right - rcReliefEl.left), 0.0, 1.0);
        relief_light_el = (float)(0.05 + t * (1.5 - 0.05));   // ~3deg..86deg above horizon
        recolorPhaseNow();
    }
    void setReliefStrFromX(int x) {
        double t = std::clamp((double)(x - rcReliefStr.left) / std::max(1L, rcReliefStr.right - rcReliefStr.left), 0.0, 1.0);
        if (de_overlay_on) de_scale = (float)(DE_SCALE_MIN + t * (DE_SCALE_MAX - DE_SCALE_MIN));
        else relief_strength = (float)(t * RELIEF_STR_MAX);
        recolorPhaseNow();
    }
    // Immediate cheap re-colour + present for the current phase (paused-drag feedback).
    void recolorPhaseNow() {
        color_phase = std::fmod(color_phase, (float)colP); if (color_phase < 0) color_phase += colP;
        if (nav) { nav->RecolorPhase(bitmap.data()); buildDisplay(); }
        InvalidateRect(hwnd, nullptr, FALSE);
    }

    // ---- geometry ----
    int formulaDockWidth() const {
        return formulaEditor && formulaEditor->visible()
            ? S(FormulaEditorPanel::DESIGN_WIDTH) : 0;
    }
    int panelRight() const {
        RECT rc{}; GetClientRect(hwnd, &rc);
        return rc.right - formulaDockWidth();
    }
    RECT viewRect() const {
        RECT rc; GetClientRect(hwnd, &rc);
        int aw = std::max(1, panelRight() - S(PANEL_W));
        int ah = std::max(1, (int)rc.bottom - S(STATUS_H));
        double s = std::min((double)aw / RENDER_W, (double)ah / RENDER_H);
        int w = std::max(4, (int)(RENDER_W * s)), h = std::max(1, (int)(RENDER_H * s));
        // The render buffer == this view size and is blitted 1:1; GDI DIB rows must be
        // DWORD-aligned, so keep the width a multiple of 4 (w*3 then divisible by 4).
        w &= ~3;
        return { (aw - w) / 2, (ah - h) / 2, (aw - w) / 2 + w, (ah - h) / 2 + h };
    }

    void layout() {
        RECT rc; GetClientRect(hwnd, &rc);
        int dockW = formulaDockWidth();
        int panelR = rc.right - dockW;
        if (formulaEditor && formulaEditor->visible())
            formulaEditor->move({ panelR, 0, rc.right, rc.bottom });
        int px = panelR - S(PANEL_W) + S(18), w = S(PANEL_W) - S(36), y = S(18);
        int g = S(8), topH = S(30), bh = S(34);
        int bw = (w - 4 * g) / 5;
        rcReset   = { px,              y, px + bw,          y + topH };
        rcRender  = { px + bw + g,     y, px + 2*bw + g,    y + topH };
        rcSave    = { px + 2*(bw+g),   y, px + 3*bw + 2*g,  y + topH };
        rcCopy    = { px + 3*(bw+g),   y, px + 4*bw + 3*g,  y + topH };
        rcPaste   = { px + 4*(bw+g),   y, px + 5*bw + 4*g,  y + topH };
        y += topH + S(14);
        int locationH = nav && nav->IsJulia() ? S(122) : S(96);
        rcLocation = { px, y, px + w, y + locationH }; y += locationH + S(16);
        rcMaxTrack = { px, y + S(24), px + w, y + S(38) }; y += S(44);
        rcDensTrack = { px, y + S(24), px + w, y + S(38) }; y += S(44);
        rcPhaseTrack = { px, y + S(24), px + w, y + S(38) }; y += S(44);
        rcSpeedTrack = { px, y + S(24), px + w, y + S(38) }; y += S(44);
        // Unlabeled toggle: box sits at y (no label above), then a uniform 14px gap.
        rcSS  = { px, y, px + w, y + bh }; y += bh + S(14);
        rcFormula = { px, y, px + w, y + bh }; y += bh + S(14);
        if (juliaUiEnabled) {
            rcJulia = { px, y, px + w, y + bh }; y += bh + S(14);
        } else {
            rcJulia = RECT{};
        }
        rcOrbitToggle = { px, y, px + w, y + bh }; y += bh + S(14);
        if (orbitOn) {
            int oh = (int)((double)w * ORBIT_H / ORBIT_W);
            rcOrbitThumb = { px, y + S(20), px + w, y + S(20) + oh };
            y += S(20) + oh + S(14);
        } else {
            rcOrbitThumb = RECT{};
        }
        // Labeled dropdowns reserve S(24) at the top for their caption (drawn at
        // box.top-S(20)) -- the same rhythm as the sliders and the gradient bar -- so
        // the gap to the next control stays uniform instead of double-gapping after a
        // dropdown and overlapping before one.
        rcColoringDD = { px, y + S(30), px + w, y + S(30) + bh }; y += S(30) + bh + S(14);
        // Relief/normal light use azimuth+elevation+strength; DE overlay uses only the
        // falloff strength (there is no light direction), so it shows a single slider.
        if (relief_on || normal_light_on) {
            rcReliefAz  = { px, y + S(24), px + w, y + S(38) }; y += S(52);
            rcReliefEl  = { px, y + S(24), px + w, y + S(38) }; y += S(52);
            rcReliefStr = { px, y + S(24), px + w, y + S(38) }; y += S(52);
        } else if (de_overlay_on) {
            rcReliefAz  = { px, y + S(24), px + w, y + S(38) }; y += S(52);
            rcReliefEl  = RECT{};
            rcReliefStr = { px, y + S(24), px + w, y + S(38) }; y += S(52);
        } else {
            rcReliefAz = rcReliefEl = rcReliefStr = RECT{};
        }
        rcPaletteDD = { px, y + S(30), px + w, y + S(30) + bh }; y += S(30) + bh + S(14);
        rcGradient = { px, y - S(4), px + w, y + S(30) }; y += S(44);
        rcColor = { px, y, px + w, y + bh }; y += bh + S(16);
        rcGalleryDD = { px, y, px + w, y + bh }; y += bh;

        // The stacked controls can exceed the window height (e.g. in Relief mode);
        // make the panel vertically scrollable. Everything above is laid out from an
        // unscrolled origin; shift all rects up by the (clamped) scroll offset.
        panelContentH = y + S(18);
        // An expanded dropdown list is an overlay drawn below its box; its height is not
        // part of the normal stacked content, so extend the scrollable range to reach the
        // bottom items (otherwise the last entry stays clipped even at max scroll). Rects
        // are still unscrolled here, so the list rects are at their natural positions.
        if (coloringOpen) panelContentH = std::max(panelContentH, (int)coloringListRect().bottom + S(18));
        if (paletteOpen)  panelContentH = std::max(panelContentH, (int)paletteListRect().bottom  + S(18));
        if (galleryOpen)  panelContentH = std::max(panelContentH, (int)galleryListRect().bottom  + S(18));
        panelViewH = rc.bottom;
        int maxs = std::max(0, panelContentH - panelViewH);
        panelScroll = std::clamp(panelScroll, 0, maxs);
        if (panelScroll != 0) {
            RECT* all[] = { &rcReset, &rcRender, &rcSave, &rcCopy, &rcPaste,
                            &rcLocation, &rcMaxTrack, &rcDensTrack, &rcPhaseTrack, &rcSpeedTrack,
                            &rcReliefAz, &rcReliefEl, &rcReliefStr,
                            &rcSS, &rcFormula, &rcJulia, &rcOrbitToggle, &rcOrbitThumb,
                            &rcColoringDD, &rcPaletteDD, &rcGradient, &rcColor, &rcGalleryDD };
            for (RECT* r : all) OffsetRect(r, 0, -panelScroll);
        }
    }

    int panelMaxScroll() const { return std::max(0, panelContentH - panelViewH); }
    void panelScrollBy(int dyDevice) {
        int maxs = panelMaxScroll();
        int ns = std::clamp(panelScroll + dyDevice, 0, maxs);
        if (ns == panelScroll) return;
        panelScroll = ns; layout(); needFull = true; InvalidateRect(hwnd, nullptr, FALSE);
    }

    void switchJuliaMode(bool enable, bool render = true) {
        if (!enable) {
            if (!nav->IsMandelbrot()) restoreMandelbrotUi(render);
            return;
        }
        if (nav->IsJulia()) {
            if (render) startRender();
            return;
        }
        if (nav->IsMandelbrot() && !hasSavedMandelUi) {
            savedColoringIdx = coloringIdx;
            savedSsOn = ssOn;
            savedOrbitOn = orbitOn;
            savedMaxIter = maxIter;
            hasSavedMandelUi = true;
        }
        if (orbitWorker) orbitWorker->cancel();
        if (orbitThumbnailWorker) orbitThumbnailWorker->cancel();
        orbitOn = false;
        orbitResult = OrbitResult{};
        orbitThumbnail.clear();
        orbitThumbnailBuilding = false;
        ssOn = false;
        coloringIdx = coloringIdx == 1 ? 1 : 0;
        relief_on = normal_light_on = de_overlay_on = 0;
        nav->SetJuliaMode(true);
        layout();
        needFull = true;
        if (render) startRender();
    }

    void restoreMandelbrotUi(bool render = true) {
        if (nav->IsMandelbrot()) return;
        if (orbitWorker) orbitWorker->cancel();
        orbitResult = OrbitResult{};
        maxIter = savedMaxIter;
        nav->SetMxit(maxIter);
        nav->RestoreMandelbrotMode();
        coloringIdx = savedColoringIdx;
        ssOn = savedSsOn;
        orbitOn = savedOrbitOn;
        relief_on = coloringIdx == 3;
        normal_light_on = coloringIdx == 4;
        de_overlay_on = coloringIdx == 6;
        hasSavedMandelUi = false;
        rebuildOrbitThumbnail();
        layout();
        needFull = true;
        if (render) startRender();
        if (render && orbitOn) requestCurrentOrbit();
    }

    bool applyFormulaConfig(const FormulaDialogConfig& candidate, bool render = true) {
        formula::ExpressionError error;
        bool entering = !nav->IsExpression();
        bool capturedMandelUi = false;
        if (entering && nav->IsMandelbrot() && !hasSavedMandelUi) {
            savedColoringIdx = coloringIdx;
            savedSsOn = ssOn;
            savedOrbitOn = orbitOn;
            savedMaxIter = maxIter;
            hasSavedMandelUi = true;
            capturedMandelUi = true;
        }
        if (!nav->SetExpressionFormula(candidate.source, candidate.pixelParameter,
                                       candidate.fixedZ0, candidate.fixedC,
                                       candidate.parameters, candidate.bailout, &error)) {
            wchar_t message[384];
            swprintf_s(message, L"Formula error at character %zu:\n%hs",
                       error.position + 1, error.message.c_str());
            MessageBoxW(hwnd, message, L"Invalid formula", MB_OK | MB_ICONERROR);
            if (capturedMandelUi) hasSavedMandelUi = false;
            return false;
        }
        if (entering) {
            maxIter = std::min(maxIter, 2000);
            nav->SetMxit(maxIter);
        }
        formulaConfig = candidate;
        if (formulaEditor && formulaEditor->visible())
            formulaEditor->setConfig(formulaConfig);
        if (orbitWorker) orbitWorker->cancel();
        orbitResult = OrbitResult{};
        rebuildOrbitThumbnail();
        ssOn = false;
        if (entering ||
            !expressionColoringIndexSupported(
                coloringIdx, nav->ExpressionSupportsDistance()))
            coloringIdx = 0;
        int expressionMethod = 0;
        if (coloringIdx == 1)
            expressionMethod = ColoringMethod::EXTERIOR_DIST_EST;
        else if (coloringIdx == 2)
            expressionMethod = ColoringMethod::STRIPE_AVERAGE;
        else if (coloringIdx == 5)
            expressionMethod = ColoringMethod::ORBIT_TRAP;
        nav->SetCMethod(expressionMethod);
        relief_on = normal_light_on = de_overlay_on = 0;
        layout();
        needFull = true;
        if (render) startRender();
        if (render && orbitOn) requestCurrentOrbit();
        return true;
    }

    void expandWindowForFormulaEditor() {
        if (formulaEditorExpansionDip > 0 || IsZoomed(hwnd)) return;
        RECT wr{}; GetWindowRect(hwnd, &wr);
        HMONITOR monitor = MonitorFromWindow(hwnd, MONITOR_DEFAULTTONEAREST);
        MONITORINFO mi{ sizeof(mi) };
        if (!GetMonitorInfoW(monitor, &mi)) return;
        int desired = S(FormulaEditorPanel::DESIGN_WIDTH);
        int currentW = wr.right - wr.left;
        int currentH = wr.bottom - wr.top;
        int grow = std::min(desired,
            std::max(0, (int)(mi.rcWork.right - mi.rcWork.left) - currentW));
        int growHeight = std::min(
            std::max(0, S(907) - currentH),
            std::max(0, (int)(mi.rcWork.bottom - mi.rcWork.top) - currentH));
        if (grow <= 0 && growHeight <= 0) return;
        int newW = currentW + grow;
        int newH = currentH + growHeight;
        int newX = std::clamp(
            wr.left - grow / 2, mi.rcWork.left, mi.rcWork.right - newW);
        int newY = std::clamp(
            wr.top - growHeight / 2,
            mi.rcWork.top, mi.rcWork.bottom - newH);
        SetWindowPos(hwnd, nullptr, newX, newY, newW, newH,
                     SWP_NOZORDER | SWP_NOACTIVATE);
        formulaEditorExpansionDip = MulDiv(grow, 96, dpi);
        formulaEditorHeightExpansionDip = MulDiv(growHeight, 96, dpi);
    }

    void closeFormulaEditor() {
        if (!formulaEditor) return;
        formulaEditor->hide();
        if ((formulaEditorExpansionDip > 0 ||
             formulaEditorHeightExpansionDip > 0)) {
            bool maximized = IsZoomed(hwnd) != FALSE;
            WINDOWPLACEMENT placement{ sizeof(placement) };
            RECT wr{};
            if (maximized) {
                if (GetWindowPlacement(hwnd, &placement))
                    wr = placement.rcNormalPosition;
            } else {
                GetWindowRect(hwnd, &wr);
            }
            HMONITOR monitor = MonitorFromRect(
                &wr, MONITOR_DEFAULTTONEAREST);
            MONITORINFO mi{ sizeof(mi) };
            bool restored = false;
            if (wr.right > wr.left && wr.bottom > wr.top &&
                GetMonitorInfoW(monitor, &mi)) {
                int width = wr.right - wr.left;
                int height = wr.bottom - wr.top;
                int shrink = std::min(
                    S(formulaEditorExpansionDip),
                    std::max(0, width - S(900)));
                int shrinkHeight = std::min(
                    S(formulaEditorHeightExpansionDip),
                    std::max(0, height - S(640)));
                if (shrink > 0 || shrinkHeight > 0) {
                    int newWidth = width - shrink;
                    int newHeight = height - shrinkHeight;
                    int newX = std::clamp(
                        wr.left + shrink / 2, mi.rcWork.left,
                        mi.rcWork.right - newWidth);
                    int newY = std::clamp(
                        wr.top + shrinkHeight / 2, mi.rcWork.top,
                        mi.rcWork.bottom - newHeight);
                    if (maximized) {
                        placement.rcNormalPosition = {
                            newX, newY, newX + newWidth, newY + newHeight
                        };
                        restored = SetWindowPlacement(
                            hwnd, &placement) != FALSE;
                    } else {
                        restored = SetWindowPos(
                            hwnd, nullptr, newX, newY,
                            newWidth, newHeight,
                            SWP_NOZORDER | SWP_NOACTIVATE) != FALSE;
                    }
                } else {
                    restored = true;
                }
            }
            if (restored) {
                formulaEditorExpansionDip = 0;
                formulaEditorHeightExpansionDip = 0;
            }
        }
        layout();
        if (!inSizeMove && !orbitBench) retargetToView();
        needFull = true;
        InvalidateRect(hwnd, nullptr, TRUE);
    }

    void showFormulaEditor() {
        if (!formulaEditor) return;
        bool opening = !formulaEditor->visible();
        if (opening) {
            formulaEditor->show(formulaConfig);
            expandWindowForFormulaEditor();
            layout();
            if (!inSizeMove && !orbitBench) retargetToView();
            needFull = true;
            InvalidateRect(hwnd, nullptr, TRUE);
        } else {
            closeFormulaEditor();
        }
    }
    // Scrollbar thumb rect on the panel's right edge (empty if nothing to scroll).
    RECT panelScrollbarRect() const {
        int maxs = panelMaxScroll();
        if (maxs <= 0 || panelContentH <= 0) return RECT{};
        RECT rc; GetClientRect(hwnd, &rc);
        int sbw = S(7);
        int trackX = panelRight() - sbw - S(2);
        int trackTop = S(2), trackBot = rc.bottom - S(2);
        int trackH = std::max(1, trackBot - trackTop);
        int thumbH = std::max(S(28), (int)((double)trackH * panelViewH / panelContentH));
        thumbH = std::min(thumbH, trackH);
        int thumbY = trackTop + (int)((double)(trackH - thumbH) * panelScroll / maxs);
        return { trackX, thumbY, trackX + sbw, thumbY + thumbH };
    }

    void createFonts() {
        if (fUi) {
            DeleteObject(fUi); DeleteObject(fBold); DeleteObject(fSmall);
            DeleteObject(fMono); DeleteObject(fMicro);
        }
        fUi = CreateFontW(-S(12),0,0,0,FW_NORMAL,0,0,0,
            DEFAULT_CHARSET,0,0,CLEARTYPE_QUALITY,0,L"Segoe UI");
        fBold = CreateFontW(-S(12),0,0,0,FW_SEMIBOLD,0,0,0,
            DEFAULT_CHARSET,0,0,CLEARTYPE_QUALITY,0,L"Segoe UI");
        fSmall = CreateFontW(-S(11),0,0,0,FW_NORMAL,0,0,0,
            DEFAULT_CHARSET,0,0,CLEARTYPE_QUALITY,0,L"Segoe UI");
        fMono = CreateFontW(-S(12),0,0,0,FW_MEDIUM,0,0,0,
            DEFAULT_CHARSET,0,0,CLEARTYPE_QUALITY,0,L"Cascadia Mono");
        fMicro = CreateFontW(-S(11),0,0,0,FW_NORMAL,0,0,0,
            DEFAULT_CHARSET,0,0,CLEARTYPE_QUALITY,0,L"Segoe UI");
    }

    // ---- view compositing ----
    void buildDisplay() {
        auto b = nav->GetBoundary();
        double lx = b[0][0], ly = b[0][1], rx = b[2][0], ry = b[2][1];
        // The reference/bitmap is y-up (row 0 = smallest imaginary = bottom) while
        // the display DIB is top-down; flip vertically so the image is upright and
        // pan/zoom track the cursor consistently on the y axis. Parallelised: this
        // runs every animation/pan frame, so it must stay well under the 16 ms budget.
#pragma omp parallel for schedule(static)
        for (int y = 0; y < renderH; ++y) {
            double fy = ly + (ry - ly) * ((renderH - 1 - y) + 0.5) / renderH;
            int sy = (int)(fy * renderH);
            for (int x = 0; x < renderW; ++x) {
                double fx = lx + (rx - lx) * (x + 0.5) / renderW;
                int sx = (int)(fx * renderW);
                size_t d = ((size_t)y * renderW + x) * 3;
                if (sx < 0 || sx >= renderW || sy < 0 || sy >= renderH) {
                    display[d] = display[d+1] = display[d+2] = 0;
                } else {
                    size_t s = ((size_t)sy * renderW + sx) * 3;
                    display[d] = bitmap[s+2]; display[d+1] = bitmap[s+1]; display[d+2] = bitmap[s];
                }
            }
        }
    }
    // Parallel bilinear upscale of `display` (renderW x renderH) to the view size.
    // Only used transiently while a window resize is in progress (before the render
    // buffers retarget); in the native steady state the blit is 1:1 (no scaling).
    void scaleDisplay(int vW, int vH) {
        // DIB scanlines must be DWORD-aligned; pad each row to a 4-byte multiple so
        // StretchDIBits reads it with the stride it derives from biWidth (otherwise
        // rows drift -> shear + channel separation for non-4-aligned widths).
        int stride = ((vW * 3) + 3) & ~3;
        if (scaledW != vW || scaledH != vH || (int)scaled.size() != stride * vH) {
            scaled.assign((size_t)stride * vH, 0); scaledW = vW; scaledH = vH; scaledStride = stride;
        }
        const uint8_t* src = display.data();
        uint8_t* dst = scaled.data();
        const double sxr = (double)renderW / vW, syr = (double)renderH / vH;
#pragma omp parallel for schedule(static)
        for (int y = 0; y < vH; ++y) {
            double fsy = (y + 0.5) * syr - 0.5;
            int y0 = (int)std::floor(fsy); double ty = fsy - y0;
            int y0c = std::clamp(y0, 0, renderH - 1), y1c = std::clamp(y0 + 1, 0, renderH - 1);
            const uint8_t* r0 = src + (size_t)y0c * renderW * 3;
            const uint8_t* r1 = src + (size_t)y1c * renderW * 3;
            uint8_t* d = dst + (size_t)y * stride;
            for (int x = 0; x < vW; ++x) {
                double fsx = (x + 0.5) * sxr - 0.5;
                int x0 = (int)std::floor(fsx); double tx = fsx - x0;
                int x0c = std::clamp(x0, 0, renderW - 1) * 3, x1c = std::clamp(x0 + 1, 0, renderW - 1) * 3;
                double w00 = (1 - tx) * (1 - ty), w10 = tx * (1 - ty), w01 = (1 - tx) * ty, w11 = tx * ty;
                for (int ch = 0; ch < 3; ++ch)
                    d[x * 3 + ch] = (uint8_t)(r0[x0c + ch] * w00 + r0[x1c + ch] * w10 +
                                              r1[x0c + ch] * w01 + r1[x1c + ch] * w11 + 0.5);
            }
        }
    }
    void commitView() {
        buildDisplay();
        // display is top-down/upright; bitmap is y-up, so flip rows back when baking
        // the warped preview as the new base (keeps the transition frame upright).
        for (int y = 0; y < renderH; ++y) {
            const uint8_t* src = &display[(size_t)y * renderW * 3];
            uint8_t* dst = &bitmap[(size_t)(renderH - 1 - y) * renderW * 3];
            for (int x = 0; x < renderW; ++x) {
                dst[x*3] = src[x*3+2]; dst[x*3+1] = src[x*3+1]; dst[x*3+2] = src[x*3];
            }
        }
        nav->UpdateCoords();
        needFull = true;
        startRender();
        if (orbitOn && inRect(viewRect(), lastOrbitX, lastOrbitY)) {
            lastOrbitRequest = {};
            requestOrbit(lastOrbitX, lastOrbitY, false);
        }
    }
    void startRender() {
        renderStart = std::chrono::steady_clock::now();
        wasComputing = true;
        nav->StartCompute();
        keepLive();
        InvalidateRect(hwnd, nullptr, FALSE);
    }
    // Re-point the render buffers at a new view size (native-resolution display).
    // Rebuilds the navigator's sample buffers + our RGB buffers, then re-renders.
    void retargetRender(int vW, int vH) {
        if (!nav || vW < 8 || vH < 8) return;
        if (vW == renderW && vH == renderH) return;
        renderW = vW; renderH = vH;
        bitmap.assign((size_t)renderW * renderH * 3, 0);
        display.assign((size_t)renderW * renderH * 3, 0);
        nav->Resize(renderW, renderH);
        startRender();
    }
    void retargetToView() {
        RECT vr = viewRect();
        int vW = vr.right - vr.left, vH = vr.bottom - vr.top;
        // Cap the native render budget: the engine holds _sub*_sub floats/pixel in two
        // arrays, so an unbounded maximized window on a huge/high-DPI monitor could
        // allocate pathologically (and an OOM would be fatal). 4 Mpx keeps the two
        // sample arrays around ~0.8 GB (_sub=5) and covers typical large windows incl.
        // ~1440p natively; only beyond it (approaching 4K maximized) do we render a bit
        // smaller and let the 1:1 blit fall back to the scaled path.
        const double MAX_PX = 4.0e6;
        double px = (double)vW * (double)vH;
        if (px > MAX_PX) {
            double s = std::sqrt(MAX_PX / px);
            vW = std::max(4, (int)(vW * s)) & ~3;
            vH = std::max(1, (int)(vH * s));
        }
        retargetRender(vW, vH);
    }
    POINT mapToRender(int x, int y) const {
        RECT r = viewRect();
        int px = (int)((x - r.left) * (double)renderW / std::max(1L, r.right - r.left));
        int py = (int)((y - r.top) * (double)renderH / std::max(1L, r.bottom - r.top));
        return { std::clamp(px, 0, renderW - 1), std::clamp(py, 0, renderH - 1) };
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
        gradientFocused = true;
        SetFocus(hwnd);
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
        if (i < 0) return;
        palette.selected = i;
        gradientDeleteSelected();
        InvalidateRect(hwnd, nullptr, FALSE);
    }
    void gradientDeleteSelected() {
        if (palette.stops.size() <= 2 ||
            palette.selected < 0 ||
            palette.selected >= (int)palette.stops.size())
            return;
        int i = palette.selected;
        palette.stops.erase(palette.stops.begin() + i);
        palette.selected = std::min(i, (int)palette.stops.size() - 1);
        palette.rebuild(); nav->SetRedisplay();
        InvalidateRect(hwnd, nullptr, FALSE);
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
        std::string location = nav->GetLocationText();
        if (nav->IsExpression()) {
            char line[160];
            location += "\r\npixel: ";
            location += formulaConfig.pixelParameter == FormulaParameter::InitialZ ? "z0" : "c";
            snprintf(line, sizeof(line),
                     "\r\nfixed_z0_re: %.17g\r\nfixed_z0_im: %.17g"
                     "\r\nfixed_c_re: %.17g\r\nfixed_c_im: %.17g"
                     "\r\nbailout: %.17g",
                     formulaConfig.fixedZ0.real(), formulaConfig.fixedZ0.imag(),
                     formulaConfig.fixedC.real(), formulaConfig.fixedC.imag(),
                     formulaConfig.bailout);
            location += line;
            for (int i = 0; i < 8; ++i) {
                snprintf(line, sizeof(line), "\r\np%d_re: %.17g\r\np%d_im: %.17g",
                         i, formulaConfig.parameters[i].real(),
                         i, formulaConfig.parameters[i].imag());
                location += line;
            }
        }
        std::wstring t = widen(location);
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
        std::string mode = val("mode:");
        std::string xs = val("x:"), ys = val("y:"), zs = val("zoom:");
        if (xs.empty() || ys.empty() || zs.empty()) return;
        std::string scale = expandSci(zs);
        if (scale.empty()) return;
        auto validNumbers = [](const std::vector<std::string>& values,
                               int positiveIndex = -1) {
            size_t chars = 0; for (const std::string& value : values) chars = std::max(chars, value.size());
            mp_bitcnt_t precision = std::max<mp_bitcnt_t>(64, (mp_bitcnt_t)(chars * 3.3219) + 40);
            mpf_t parsed; mpf_init2(parsed, precision);
            bool ok = true;
            for (int i = 0; i < (int)values.size(); ++i) {
                if (mpf_set_str(parsed, values[i].c_str(), 10) != 0 ||
                    (i == positiveIndex && mpf_sgn(parsed) <= 0)) {
                    ok = false; break;
                }
            }
            mpf_clear(parsed);
            return ok;
        };
        if (!validNumbers({ xs, ys, scale }, 2)) return;
        if (mode == "expression") {
            FormulaDialogConfig candidate;
            candidate.source = val("formula:");
            std::string pixel = val("pixel:");
            candidate.pixelParameter = pixel == "z0"
                ? FormulaParameter::InitialZ : FormulaParameter::C;
            std::vector<std::string> values = {
                val("fixed_z0_re:"), val("fixed_z0_im:"),
                val("fixed_c_re:"), val("fixed_c_im:"), val("bailout:")
            };
            for (int i = 0; i < 8; ++i) {
                values.push_back(val(("p" + std::to_string(i) + "_re:").c_str()));
                values.push_back(val(("p" + std::to_string(i) + "_im:").c_str()));
            }
            if (candidate.source.empty() ||
                std::any_of(values.begin(), values.end(),
                            [](const std::string& value) { return value.empty(); }) ||
                !validNumbers(values))
                return;
            candidate.fixedZ0 = { atof(values[0].c_str()), atof(values[1].c_str()) };
            candidate.fixedC = { atof(values[2].c_str()), atof(values[3].c_str()) };
            candidate.bailout = atof(values[4].c_str());
            if (!(candidate.bailout > 0.0)) return;
            for (int i = 0; i < 8; ++i)
                candidate.parameters[i] = {
                    atof(values[5 + 2 * i].c_str()),
                    atof(values[6 + 2 * i].c_str())
                };
            if (!applyFormulaConfig(candidate, false)) return;
        } else if (mode == "julia") {
            if (!juliaUiEnabled) return;
            std::string cr = val("c_re:"), ci = val("c_im:");
            if (cr.empty() || ci.empty()) return;
            if (!validNumbers({ cr, ci })) return;
            bool wasJulia = nav->IsJulia();
            switchJuliaMode(true, false);
            if (!nav->SetJuliaC(cr, ci)) {
                if (!wasJulia) switchJuliaMode(false, false);
                return;
            }
        } else if (!nav->IsMandelbrot()) {
            restoreMandelbrotUi(false);
        }
        if (nav->SetLocation(xs, ys, scale)) {
            needFull = true;
            startRender();
        }
    }
    void saveImage() {
        wchar_t file[MAX_PATH]; wcscpy_s(file, defaultSaveName().c_str());
        OPENFILENAMEW ofn{ sizeof(ofn) };
        ofn.hwndOwner = hwnd; ofn.lpstrFilter = L"PNG image (*.png)\0*.png\0";
        ofn.lpstrFile = file; ofn.nMaxFile = MAX_PATH; ofn.lpstrDefExt = L"png";
        ofn.Flags = OFN_OVERWRITEPROMPT | OFN_PATHMUSTEXIST;
        if (GetSaveFileNameW(&ofn)) {
            buildDisplay();
            std::vector<uint8_t> rgb((size_t)renderW * renderH * 3);   // display is BGR -> RGB
            for (size_t i = 0; i < (size_t)renderW * renderH; ++i) {
                rgb[i * 3] = display[i * 3 + 2]; rgb[i * 3 + 1] = display[i * 3 + 1]; rgb[i * 3 + 2] = display[i * 3];
            }
            writeExportPNG(file, rgb, renderW, renderH);
        }
    }

    // ---- widget drawing ----
    void drawButton(HDC dc, RECT r, const std::wstring& s, Hit id,
                   bool accent, bool activeOutline = false) {
        bool hov = hover == id, prs = pressed == id;
        COLORREF bg = accent
            ? (prs ? RGB(70,120,220) : hov ? CLR_ACCENT_HI : CLR_ACCENT)
            : (prs || hov ? CLR_CARD_HOV : CLR_CARD);
        COLORREF border = activeOutline ? CLR_ACCENT
            : (accent ? bg : CLR_BORDER);
        fillRound(dc, r, bg, border, S(8));
        if (activeOutline) {
            RECT stripe = {r.left + S(1), r.top + S(5),
                          r.left + S(4), r.bottom - S(5)};
            fillRound(dc, stripe, CLR_ACCENT, CLR_ACCENT, S(3));
        }
        HFONT font = id >= H_RESET && id <= H_PASTE ? fMicro : fUi;
        drawText(dc, r, s, accent ? RGB(255,255,255) : CLR_TEXT, font,
                 DT_CENTER | DT_VCENTER | DT_SINGLELINE);
    }
    void drawSlider(HDC dc, RECT track, double t, Hit id,
                    const double* snaps = nullptr, int nsnaps = 0) {
        int cy = (track.top + track.bottom) / 2, hb = S(3);
        RECT bar = { track.left, cy - hb, track.right, cy + hb };
        fillRound(dc, bar, CLR_TRACK, CLR_BG, S(6));
        int kx = track.left + (int)(t * (track.right - track.left));
        RECT fill = { track.left, cy - hb, kx, cy + hb };
        fillRound(dc, fill, CLR_ACCENT, CLR_BG, S(6));
        // Snap detent marker lines: a thin tick straddling the track at each snap.
        if (snaps && nsnaps > 0) {
            HPEN dp = CreatePen(PS_SOLID, S(1), CLR_TEXT_DIM);
            HGDIOBJ odp = SelectObject(dc, dp);
            for (int i = 0; i < nsnaps; ++i) {
                int sx = track.left + (int)(snaps[i] * (track.right - track.left));
                MoveToEx(dc, sx, cy - S(9), nullptr); LineTo(dc, sx, cy + S(9));
            }
            SelectObject(dc, odp); DeleteObject(dp);
        }
        int kr = hover == id || pressed == id ? S(9) : S(7);
        HBRUSH kb = CreateSolidBrush(RGB(238,242,250));
        HPEN kp = CreatePen(PS_SOLID, S(2), CLR_ACCENT);
        HGDIOBJ ob = SelectObject(dc, kb), op = SelectObject(dc, kp);
        Ellipse(dc, kx - kr, cy - kr, kx + kr, cy + kr);
        SelectObject(dc, ob); SelectObject(dc, op); DeleteObject(kb); DeleteObject(kp);
    }
    // Snap a normalized slider value t to the nearest snap point within a
    // distance-aware radius: the radius at each snap is a fraction of the gap to
    // its nearest neighbour (or the track edge), clamped to a comfortable pixel
    // range. Returns the (possibly snapped) value. trackPx = track pixel width.
    double snapT(double t, const double* snaps, int nsnaps, double trackPx) {
        double best = t, bestDpx = 1e18;
        for (int i = 0; i < nsnaps; ++i) {
            double s = snaps[i];
            double lo = s, hi = 1.0 - s;                 // gaps to the track edges
            for (int k = 0; k < nsnaps; ++k) if (k != i) {
                double d = s - snaps[k];
                if (d > 0) lo = std::min(lo, d); else hi = std::min(hi, -d);
            }
            double radiusPx = std::clamp(std::min(lo, hi) * 0.18 * trackPx, 6.0, 28.0);
            double dpx = std::fabs(t - s) * trackPx;
            if (dpx <= radiusPx && dpx < bestDpx) { bestDpx = dpx; best = s; }
        }
        return best;
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
    void drawFormulaButton(HDC dc) {
        bool open = formulaEditor && formulaEditor->visible();
        bool active = open || nav->IsExpression();
        bool hov = hover == H_FORMULA, prs = pressed == H_FORMULA;
        fillRound(dc, rcFormula,
                  prs || hov ? CLR_CARD_HOV : CLR_CARD,
                  active ? CLR_ACCENT : CLR_BORDER, S(8));
        if (active) {
            RECT stripe = {
                rcFormula.left + S(1), rcFormula.top + S(5),
                rcFormula.left + S(4), rcFormula.bottom - S(5)
            };
            fillRound(dc, stripe, CLR_ACCENT, CLR_ACCENT, S(3));
        }
        RECT title = rcFormula;
        title.left += S(12);
        title.right -= S(72);
        drawText(dc, title, L"Formula editor", CLR_TEXT, fUi,
                 DT_LEFT | DT_VCENTER | DT_SINGLELINE);
        const wchar_t* state = open ? L"OPEN"
            : (nav->IsExpression() ? L"CUSTOM" : L"");
        if (*state) {
            RECT status = rcFormula;
            status.left = status.right - S(68);
            status.right -= S(12);
            drawText(dc, status, state, CLR_ACCENT_HI, fMicro,
                     DT_RIGHT | DT_VCENTER | DT_SINGLELINE);
        }
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
    // ---- coloring-mode dropdown (Smooth / Distance / Feather / Relief) ----
    static const wchar_t* coloringName(int i) {
        static const wchar_t* n[7] = { L"Smooth iteration", L"Distance (EDE)", L"Feather (stripe)", L"Relief (3D light)", L"Normal light (3D)", L"Orbit trap", L"DE + smooth" };
        return n[i < 0 ? 0 : (i > 6 ? 6 : i)];
    }
    int coloringItemH() const { return S(28); }
    RECT coloringListRect() const {
        int n = 7;
        return { rcColoringDD.left, rcColoringDD.bottom + S(2), rcColoringDD.right,
                 rcColoringDD.bottom + S(2) + n * coloringItemH() };
    }
    int coloringItemAt(int x, int y) const {
        RECT lr = coloringListRect();
        if (x < lr.left || x > lr.right || y < lr.top || y > lr.bottom) return -1;
        int i = (y - lr.top) / coloringItemH();
        return (i < 0 || i > 6) ? -1 : i;
    }
    void drawColoringDD(HDC dc) {
        bool hov = hover == H_EDE;
        fillRound(dc, rcColoringDD, hov ? CLR_CARD_HOV : CLR_CARD,
                  coloringOpen ? CLR_ACCENT : CLR_BORDER, S(8));
        RECT tr = rcColoringDD; tr.left += S(12); tr.right -= S(28);
        const wchar_t* currentName =
            nav->IsExpression() && coloringIdx == 0 &&
                    !nav->ExpressionSupportsDistance()
                ? L"Iteration (raw)" : coloringName(coloringIdx);
        drawText(dc, tr, currentName, CLR_TEXT, fUi,
                 DT_LEFT | DT_VCENTER | DT_SINGLELINE);
        RECT cr = rcColoringDD; cr.right -= S(12);
        drawText(dc, cr, coloringOpen ? L"\u25B2" : L"\u25BC", CLR_TEXT_DIM, fSmall, DT_RIGHT | DT_VCENTER | DT_SINGLELINE);
    }
    void drawColoringList(HDC dc) {
        RECT lr = coloringListRect();
        fillRound(dc, lr, CLR_CARD, CLR_ACCENT, S(8));
        int ih = coloringItemH();
        for (int i = 0; i < 7; ++i) {
            RECT ir = { lr.left, lr.top + i * ih, lr.right, lr.top + (i + 1) * ih };
            if (i == coloringHover) fillRect(dc, { ir.left + S(3), ir.top, ir.right - S(3), ir.bottom }, CLR_CARD_HOV);
            RECT tr = ir; tr.left += S(14);
            COLORREF color = ((nav->IsJulia() && !nav->IsExpression() &&
                               i > 1) ||
                              (nav->IsExpression() &&
                               !expressionColoringIndexSupported(
                                   i, nav->ExpressionSupportsDistance())))
                ? CLR_TEXT_DIM : (i == coloringIdx ? CLR_ACCENT : CLR_TEXT);
            const wchar_t* name =
                nav->IsExpression() && i == 0 &&
                        !nav->ExpressionSupportsDistance()
                    ? L"Iteration (raw)" : coloringName(i);
            drawText(dc, tr, name, color, fUi,
                     DT_LEFT | DT_VCENTER | DT_SINGLELINE);
        }
    }
    void selectColoring(int idx) {
        if (idx < 0 || idx > 6) return;
        if (nav->IsJulia() && !nav->IsExpression() && idx > 1) return;
        if (nav->IsExpression() &&
            !expressionColoringIndexSupported(
                idx, nav->ExpressionSupportsDistance()))
            return;
        coloringIdx = idx;
        relief_on = (idx == 3) ? 1 : 0;
        normal_light_on = (idx == 4) ? 1 : 0;
        de_overlay_on = (idx == 6) ? 1 : 0;
        int m = nav->GetCMethod();
        m &= ~(ColoringMethod::EXTERIOR_DIST_EST | ColoringMethod::STRIPE_AVERAGE | ColoringMethod::NORMAL_MAP | ColoringMethod::ORBIT_TRAP | ColoringMethod::DE_OVERLAY);
        if (idx == 1) m |= ColoringMethod::EXTERIOR_DIST_EST;
        else if (idx == 2) m |= ColoringMethod::STRIPE_AVERAGE;
        else if (idx == 4) m |= ColoringMethod::NORMAL_MAP;
        else if (idx == 5) m |= ColoringMethod::ORBIT_TRAP;
        else if (idx == 6) m |= ColoringMethod::DE_OVERLAY;
        // idx 3 (Relief) uses the plain smooth field (method 0) + slope post-shade;
        // idx 4 (Normal light) uses the engine's analytic normal + Lambert shade;
        // idx 5 (Orbit trap) outputs a trap palette coordinate (BLA off -> slow deep);
        // idx 6 (DE + smooth) draws the distance-estimate B&W layer over the smooth base.
        nav->SetCMethod(m);
        layout();   // light/strength slider appears/disappears -> reflow the panel
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
        if (!nav->IsMandelbrot()) restoreMandelbrotUi(false);
        const Preset& p = gp[idx];
        // render settings first, then the location (which kicks off the render)
        maxIter = std::clamp(p.maxIter, 100, 5000000); nav->SetMxit(maxIter);
        color_density = std::clamp(p.density, 10.0f, 200.0f);
        ssOn = p.ss;
        coloringIdx = p.coloring;
        // Keep the coloring flags + method mask consistent with the preset's coloring
        // (as selectColoring does), otherwise a leftover relief/normal/DE overlay from
        // the previous mode stays active -> a doubled image and its sliders never leave.
        relief_on       = (coloringIdx == 3) ? 1 : 0;
        normal_light_on = (coloringIdx == 4) ? 1 : 0;
        de_overlay_on   = (coloringIdx == 6) ? 1 : 0;
        int m = 0;
        if (ssOn) m |= ColoringMethod::SUPER_SAMPLING;
        if (coloringIdx == 1) m |= ColoringMethod::EXTERIOR_DIST_EST;
        else if (coloringIdx == 2) m |= ColoringMethod::STRIPE_AVERAGE;
        else if (coloringIdx == 4) m |= ColoringMethod::NORMAL_MAP;
        else if (coloringIdx == 5) m |= ColoringMethod::ORBIT_TRAP;
        else if (coloringIdx == 6) m |= ColoringMethod::DE_OVERLAY;
        nav->SetCMethod(m);
        // palette by name
        const auto& pr = palettePresets();
        for (int i = 0; i < (int)pr.size(); ++i)
            if (wcscmp(pr[i].name, p.palette) == 0) { paletteIdx = i; palette.load(i); break; }
        layout();   // light sliders appear/disappear with the preset's coloring -> reflow
        std::string scale = expandSci(p.zoom);
        if (!scale.empty() && nav->SetLocation(p.x, p.y, scale)) {
            needFull = true;
            startRender();
        }
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
        int vW = vr.right - vr.left, vH = vr.bottom - vr.top;
        BITMAPINFO bi{}; bi.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
        bi.bmiHeader.biPlanes = 1; bi.bmiHeader.biBitCount = 24; bi.bmiHeader.biCompression = BI_RGB;
        // Native steady state: the render buffer matches the view exactly -> 1:1 blit
        // of `display` (no resampling, no beating). The scale/GDI paths are only hit
        // transiently while a window resize is in flight (before the buffers retarget).
        static int stretchEnv = -1;
        if (stretchEnv < 0) { const char* e = getenv("MANDEL_STRETCH"); stretchEnv = e ? atoi(e) : 1; }
        if (vW == renderW && vH == renderH) {
            bi.bmiHeader.biWidth = renderW; bi.bmiHeader.biHeight = -renderH;
            StretchDIBits(dc, vr.left, vr.top, vW, vH, 0, 0, renderW, renderH,
                          display.data(), &bi, DIB_RGB_COLORS, SRCCOPY);   // 1:1 copy
        } else if (stretchEnv == 1 || stretchEnv == 3) {
            scaleDisplay(vW, vH);
            bi.bmiHeader.biWidth = vW; bi.bmiHeader.biHeight = -vH;
            StretchDIBits(dc, vr.left, vr.top, vW, vH, 0, 0, vW, vH,
                          scaled.data(), &bi, DIB_RGB_COLORS, SRCCOPY);   // 1:1 copy
        } else {
            bi.bmiHeader.biWidth = renderW; bi.bmiHeader.biHeight = -renderH;
            SetStretchBltMode(dc, stretchEnv == 2 ? HALFTONE : COLORONCOLOR);
            StretchDIBits(dc, vr.left, vr.top, vW, vH,
                          0, 0, renderW, renderH, display.data(), &bi, DIB_RGB_COLORS, SRCCOPY);
        }

        // panel + all controls: only rebuilt when UI state changed. On fractal/
        // compute/drag ticks the panel persists in the cached buffer.
        if (full) {
        // panel
        int panelR = panelRight();
        RECT panel = { panelR - S(PANEL_W), 0, panelR, rc.bottom };
        fillRect(dc, panel, CLR_PANEL);
        RECT panelRail = {
            panel.right - S(11), panel.top,
            panel.right, panel.bottom
        };
        fillRect(dc, panelRail, RGB(21, 24, 33));
        HPEN sep = CreatePen(PS_SOLID, 1, CLR_BORDER);
        HGDIOBJ osep = SelectObject(dc, sep);
        MoveToEx(dc, panel.left, 0, nullptr); LineTo(dc, panel.left, rc.bottom);
        SelectObject(dc, osep); DeleteObject(sep);

        drawButton(dc, rcReset, L"Reset", H_RESET, false);
        drawButton(dc, rcRender, L"Render", H_RENDER, false);
        drawButton(dc, rcSave, L"Export", H_SAVE, false);
        drawButton(dc, rcCopy, L"Copy", H_COPY, false);
        drawButton(dc, rcPaste, L"Paste", H_PASTE, false);

        drawGalleryDD(dc);

        // location card -- keep mode, x, y and zoom on four stable lines
        fillRound(dc, rcLocation, CLR_CARD, CLR_BORDER, S(8));
        RECT lt = rcLocation; lt.left += S(10); lt.top += S(6); lt.right -= S(10); lt.bottom -= S(6);
        {
            std::string location = nav->GetLocationText();
            if (nav->IsExpression()) {
                size_t formulaLine = location.find("\r\nformula:");
                if (formulaLine != std::string::npos)
                    location.resize(formulaLine);
            }
            HGDIOBJ ofm = SelectObject(dc, fMono);
            SIZE cs{}; GetTextExtentPoint32W(dc, L"0", 1, &cs);
            SelectObject(dc, ofm);
            int cpl = cs.cx > 0 ? (int)((lt.right - lt.left) / cs.cx) : 40;
            drawText(dc, lt, truncLocation(location, cpl), CLR_TEXT, fMono,
                     DT_LEFT | DT_TOP | DT_NOPREFIX);
        }

        // max iterations (read-only value; drag slider + mouse-wheel to change)
        labelRow(dc, rcMaxTrack, L"Max iterations", std::to_wstring(maxIter));
        drawSlider(dc, rcMaxTrack, maxToT(), H_MAXTRACK);

        // color density (read-only; drag slider = integer, mouse-wheel = 0.1)
        wchar_t db[32]; swprintf_s(db, L"%.1f", color_density);
        labelRow(dc, rcDensTrack, L"Color density", db);
        drawSlider(dc, rcDensTrack, std::clamp((color_density - 10.0) / 190.0, 0.0, 1.0), H_DENSTRACK);

        // palette-cycle phase (0..1 = one full cycle) and animation speed
        wchar_t pb[32]; swprintf_s(pb, L"%.2f", color_phase / (float)colP);
        labelRow(dc, rcPhaseTrack, L"Gradient phase", pb);
        drawSlider(dc, rcPhaseTrack, std::clamp(color_phase / (double)colP, 0.0, 1.0), H_PHASETRACK);

        wchar_t sb[32];
        if (std::fabs(animSpeed) < 1e-3f) swprintf_s(sb, L"paused");
        else swprintf_s(sb, L"%+.2f", animSpeed);
        labelRow(dc, rcSpeedTrack, L"Cycle speed", sb);
        const double speedSnaps[] = { 0.5 };
        drawSlider(dc, rcSpeedTrack, std::clamp((animSpeed + 1.0) / 2.0, 0.0, 1.0), H_SPEEDTRACK,
                   speedSnaps, 1);

        drawToggle(dc, rcSS, L"5x supersampling", ssOn, H_SS);
        drawFormulaButton(dc);
        if (juliaUiEnabled)
            drawToggle(dc, rcJulia, L"Julia set (experimental)", nav->IsJulia(), H_JULIA);
        drawToggle(dc, rcOrbitToggle, L"Orbit overlay", orbitOn, H_ORBIT);
        if (orbitOn) drawOrbitThumbnail(dc);
        label(dc, rcColoringDD.left, rcColoringDD.top - S(30), L"Coloring");
        drawColoringDD(dc);

        if (relief_on || normal_light_on) {
            wchar_t rab[32]; swprintf_s(rab, L"%.0f\u00B0", relief_light_az * 180.0f / 3.14159265f);
            labelRow(dc, rcReliefAz, L"Light azimuth", rab);
            drawSlider(dc, rcReliefAz, std::clamp(relief_light_az / (2.0 * 3.14159265358979), 0.0, 1.0), H_RELAZ);
            wchar_t reb[32]; swprintf_s(reb, L"%.0f\u00B0", relief_light_el * 180.0f / 3.14159265f);
            labelRow(dc, rcReliefEl, L"Light elevation", reb);
            drawSlider(dc, rcReliefEl, std::clamp((relief_light_el - 0.05) / (1.5 - 0.05), 0.0, 1.0), H_RELEL);
            wchar_t rsb[32]; swprintf_s(rsb, L"%.2f", relief_strength);
            labelRow(dc, rcReliefStr, L"Relief strength", rsb);
            drawSlider(dc, rcReliefStr, std::clamp(relief_strength / RELIEF_STR_MAX, 0.0, 1.0), H_RELSTR);
        } else if (de_overlay_on) {
            wchar_t kb[32]; swprintf_s(kb, L"%.2f", de_k);
            labelRow(dc, rcReliefAz, L"DE boundary scale", kb);
            drawSlider(dc, rcReliefAz, std::clamp((de_k - DE_K_MIN) / (DE_K_MAX - DE_K_MIN), 0.0, 1.0), H_RELAZ);
            wchar_t scb[32]; swprintf_s(scb, L"%.2f", de_scale);
            labelRow(dc, rcReliefStr, L"DE falloff scale", scb);
            drawSlider(dc, rcReliefStr, std::clamp((de_scale - DE_SCALE_MIN) / (DE_SCALE_MAX - DE_SCALE_MIN), 0.0, 1.0), H_RELSTR);
        }

        label(dc, rcPaletteDD.left, rcPaletteDD.top - S(30), L"Palette");
        drawPaletteDD(dc);

        // Gradient: click empty space to add, drag to move, Delete/Backspace or
        // right-click to remove, and double-click to edit a stop's color.
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
        if (gradientFocused) {
            HBRUSH focus = CreateSolidBrush(CLR_ACCENT);
            FrameRect(dc, &rcGradient, focus);
            DeleteObject(focus);
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
            drawText(dc, tr, L"Edit color  \u00b7  Delete removes stop",
                     CLR_TEXT, fUi,
                     DT_LEFT | DT_VCENTER | DT_SINGLELINE);
        }
        // palette dropdown list -- drawn last so it overlays the widgets below it
        if (paletteOpen) drawPaletteList(dc);
        if (coloringOpen) drawColoringList(dc);
        if (galleryOpen) drawGalleryList(dc);
        // panel scrollbar thumb (only when the control stack overflows)
        RECT sbThumb = panelScrollbarRect();
        if (sbThumb.right > sbThumb.left) {
            bool hov = hover == H_PANELSB || pressed == H_PANELSB;
            fillRound(dc, sbThumb, hov ? CLR_ACCENT : CLR_TRACK, CLR_PANEL, S(3));
        }
        } // end panel (full paints only)

        // status bar (always -- text changes with compute state / timings)
        RECT status = {
            0, rc.bottom - S(STATUS_H),
            panelRight() - S(PANEL_W), rc.bottom
        };
        fillRect(dc, status, CLR_PANEL);
        std::wstring st = nav->IsComputing()
            ? L"  Rendering..."
            : L"  Ready   last render " + std::to_wstring((int)lastRenderMs) + L" ms";
        st += L"   present " + std::to_wstring((int)(lastPresentMs + 0.5)) + L" ms";
        {
            const ComputeBackendInfo& backend = nav->GetBackendInfo();
            st += L"   backend " + widen(backend.name);
            if (backend.fallback) st += L" (fallback)";
            else if (backend.hardwareAccelerated && !nav->IsComputing() &&
                     !nav->LastComputeUsedGpuPath())
                st += L" (CPU frame)";
        }
        if (nav->IsExpression()) {
            std::string acceleration = nav->GetExpressionAccelerationText();
            if (!acceleration.empty()) st += L"   formula " + widen(acceleration);
        }
        if (animFps > 0) {
            wchar_t fb[64];
            swprintf_s(fb, L"   anim %.0f fps (recolor %.1f ms)", animFps, lastRecolorMs);
            st += fb;
        }
        if (!nav->IsComputing())
            st += L"    Drag pan   Wheel zoom   R reset   Space render   S save   C copy";
        RECT str = status; str.left += S(10);
        drawText(dc, str, st, CLR_TEXT_DIM, fSmall, DT_LEFT | DT_VCENTER | DT_SINGLELINE);

        if (full) {
            BitBlt(wdc, 0, 0, rc.right, rc.bottom, dc, 0, 0, SRCCOPY);
        } else {
            if (orbitOn) drawOrbitThumbnail(dc);
            // Fractal-only tick: only the view and the status bar changed; the panel
            // is unchanged on screen, so skip re-blitting it (saves a big copy when
            // the window is large).
            BitBlt(wdc, vr.left, vr.top, vW, vH, dc, vr.left, vr.top, SRCCOPY);
            BitBlt(wdc, status.left, status.top, status.right - status.left,
                   status.bottom - status.top, dc, status.left, status.top, SRCCOPY);
            if (orbitOn) {
                RECT orb = rcOrbitThumb; orb.top -= S(20);
                BitBlt(wdc, orb.left, orb.top, orb.right - orb.left, orb.bottom - orb.top,
                       dc, orb.left, orb.top, SRCCOPY);
            }
        }
        EndPaint(hwnd, &ps);
        lastPresentMs = std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - _pt0).count();
    }

    Hit hitTest(int x, int y) {
        { RECT sb = panelScrollbarRect(); if (sb.right > sb.left && inRect(sb, x, y)) return H_PANELSB; }
        if (inRect(rcReset,x,y)) return H_RESET;
        if (inRect(rcRender,x,y)) return H_RENDER;
        if (inRect(rcSave,x,y)) return H_SAVE;
        if (inRect(rcCopy,x,y)) return H_COPY;
        if (inRect(rcPaste,x,y)) return H_PASTE;
        if (inRect(rcGalleryDD,x,y)) return H_GALLERY_DD;
        RECT mt = rcMaxTrack; mt.top -= S(8); mt.bottom += S(8); if (inRect(mt,x,y)) return H_MAXTRACK;
        RECT dt = rcDensTrack; dt.top -= S(8); dt.bottom += S(8); if (inRect(dt,x,y)) return H_DENSTRACK;
        RECT pt = rcPhaseTrack; pt.top -= S(8); pt.bottom += S(8); if (inRect(pt,x,y)) return H_PHASETRACK;
        RECT sp = rcSpeedTrack; sp.top -= S(8); sp.bottom += S(8); if (inRect(sp,x,y)) return H_SPEEDTRACK;
        if (relief_on || normal_light_on || de_overlay_on) {
            RECT ra = rcReliefAz;  ra.top -= S(8); ra.bottom += S(8); if (inRect(ra,x,y)) return H_RELAZ;
            RECT re = rcReliefEl;  re.top -= S(8); re.bottom += S(8); if (inRect(re,x,y)) return H_RELEL;
            RECT rs = rcReliefStr; rs.top -= S(8); rs.bottom += S(8); if (inRect(rs,x,y)) return H_RELSTR;
        }
        if (inRect(rcSS,x,y)) return H_SS;
        if (inRect(rcFormula,x,y)) return H_FORMULA;
        if (juliaUiEnabled && inRect(rcJulia,x,y)) return H_JULIA;
        if (inRect(rcOrbitToggle,x,y)) return H_ORBIT;
        if (inRect(rcColoringDD,x,y)) return H_EDE;
        if (inRect(rcPaletteDD,x,y)) return H_PALETTE_DD;
        if (inRect(rcColor,x,y)) return H_COLOR;
        RECT gr = rcGradient; gr.left -= S(8); gr.right += S(8); gr.top -= S(6); gr.bottom += S(14);
        if (inRect(gr,x,y)) return H_GRADIENT;
        RECT vr = viewRect(); if (inRect(vr,x,y)) return H_VIEW;
        return H_NONE;
    }

    // Env-gated micro-benchmark: on the first settled frame, time the cached phase
    // re-colour (RecolorPhase) vs a full re-colour (UpdateBitmap) and write the
    // result to build\gui_bench.txt, then quit. MANDEL_GUI_SS selects SS on/off.
    void runRecolorBench() {
        benchMode = false;   // run once
        // Engine compute (render) time: from the last startRender() to this first settle.
        double renderMs = std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - renderStart).count();
        nav->SetRedisplay(); nav->UpdateBitmap(bitmap.data());   // ensure settled cache is built
        const int N = 60;
        auto t0 = std::chrono::steady_clock::now();
        for (int k = 0; k < N; ++k) { color_phase += 3.0f; nav->RecolorPhase(bitmap.data()); }
        double rp = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count() / N;
        auto t1 = std::chrono::steady_clock::now();
        for (int k = 0; k < N; ++k) { color_phase += 3.0f; nav->SetRedisplay(); nav->UpdateBitmap(bitmap.data()); }
        double ub = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t1).count() / N;
        // Present-path timing: buildDisplay (warp/copy) and a full synchronous paint
        // (StretchDIBits + status + BitBlt) with the panel skipped, as in animation.
        auto t2 = std::chrono::steady_clock::now();
        for (int k = 0; k < N; ++k) buildDisplay();
        double bd = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t2).count() / N;
        auto t3 = std::chrono::steady_clock::now();
        for (int k = 0; k < N; ++k) { fractalOnlyTick = true; InvalidateRect(hwnd, nullptr, FALSE); UpdateWindow(hwnd); }
        double pt = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t3).count() / N;
        RECT rcw; GetClientRect(hwnd, &rcw); RECT vr = viewRect();
        // Dump the upscaled `scaled` buffer to a BMP (stride-aware) so the DIB
        // alignment can be verified headlessly.
        if (getenv("MANDEL_GUI_DUMP")) {
            int vW = vr.right - vr.left, vH = vr.bottom - vr.top;
            buildDisplay(); scaleDisplay(vW, vH);
            FILE* bf = nullptr; fopen_s(&bf, "build\\gui_scaled.bmp", "wb");
            if (bf) {
                uint32_t ds = (uint32_t)scaledStride * vH, fs = 54 + ds, off = 54;
                uint8_t fh[14] = { 'B','M' }; memcpy(fh + 2, &fs, 4); memcpy(fh + 10, &off, 4);
                uint8_t ih[40] = { 0 }; uint32_t v; v = 40; memcpy(ih, &v, 4);
                int32_t iw = vW, ihh = vH; memcpy(ih + 4, &iw, 4); memcpy(ih + 8, &ihh, 4);
                uint16_t pl = 1, bp = 24; memcpy(ih + 12, &pl, 2); memcpy(ih + 14, &bp, 2);
                memcpy(ih + 20, &ds, 4);
                fwrite(fh, 1, 14, bf); fwrite(ih, 1, 40, bf);
                for (int y = vH - 1; y >= 0; --y) fwrite(scaled.data() + (size_t)y * scaledStride, 1, scaledStride, bf);
                fclose(bf);
            }
            // Final composited frame (RGB, renderW x renderH) exactly as shown.
            FILE* ff = nullptr; fopen_s(&ff, "build\\gui_frame.bmp", "wb");
            if (ff) {
                int stride = (renderW * 3 + 3) & ~3; uint32_t ds = (uint32_t)stride * renderH, fs = 54 + ds, off = 54;
                uint8_t fh[14] = { 'B','M' }; memcpy(fh + 2, &fs, 4); memcpy(fh + 10, &off, 4);
                uint8_t ih[40] = { 0 }; uint32_t v; v = 40; memcpy(ih, &v, 4);
                int32_t iw = renderW, ihh = renderH; memcpy(ih + 4, &iw, 4); memcpy(ih + 8, &ihh, 4);
                uint16_t pl = 1, bp = 24; memcpy(ih + 12, &pl, 2); memcpy(ih + 14, &bp, 2);
                memcpy(ih + 20, &ds, 4);
                fwrite(fh, 1, 14, ff); fwrite(ih, 1, 40, ff);
                std::vector<uint8_t> row(stride, 0);
                for (int y = renderH - 1; y >= 0; --y) {
                    for (int x = 0; x < renderW; ++x) {
                        const uint8_t* s = bitmap.data() + ((size_t)y * renderW + x) * 3;
                        row[x * 3] = s[2]; row[x * 3 + 1] = s[1]; row[x * 3 + 2] = s[0];   // RGB->BGR
                    }
                    fwrite(row.data(), 1, stride, ff);
                }
                fclose(ff);
            }
            // Actual presented frame (post-StretchDIBits) grabbed from the back buffer,
            // so the on-screen 1:1 blit can be verified headlessly (alignment / no shear).
            if (memDC && memBmp) {
                BITMAPINFO pbi{}; pbi.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
                pbi.bmiHeader.biWidth = memW; pbi.bmiHeader.biHeight = memH; pbi.bmiHeader.biPlanes = 1;
                pbi.bmiHeader.biBitCount = 24; pbi.bmiHeader.biCompression = BI_RGB;
                int pst = (memW * 3 + 3) & ~3; std::vector<uint8_t> pbuf((size_t)pst * memH);
                if (GetDIBits(memDC, memBmp, 0, memH, pbuf.data(), &pbi, DIB_RGB_COLORS)) {
                    FILE* pf = nullptr; fopen_s(&pf, "build\\gui_present.bmp", "wb");
                    if (pf) {
                        uint32_t ds = (uint32_t)pst * memH, fs = 54 + ds, off = 54;
                        uint8_t fh[14] = { 'B','M' }; memcpy(fh + 2, &fs, 4); memcpy(fh + 10, &off, 4);
                        uint8_t ih[40] = { 0 }; uint32_t v = 40; memcpy(ih, &v, 4);
                        int32_t iw = memW, ihh = memH; memcpy(ih + 4, &iw, 4); memcpy(ih + 8, &ihh, 4);
                        uint16_t pl = 1, bp = 24; memcpy(ih + 12, &pl, 2); memcpy(ih + 14, &bp, 2);
                        memcpy(ih + 20, &ds, 4);
                        fwrite(fh, 1, 14, pf); fwrite(ih, 1, 40, pf);
                        fwrite(pbuf.data(), 1, (size_t)pst * memH, pf);
                        fclose(pf);
                    }
                }
            }
        }
        FILE* f = nullptr; fopen_s(&f, "build\\gui_bench.txt", "a");
        if (f) {
            fprintf(f, "SS=%d  render=%.1f ms  recolor=%.3f ms  full=%.3f  buildDisplay=%.3f  present(paint)=%.3f  view=%ldx%ld win=%ldx%ld\n",
                    ssOn ? 1 : 0, renderMs, rp, ub, bd, pt,
                    vr.right - vr.left, vr.bottom - vr.top, rcw.right, rcw.bottom);
            fclose(f);
        }
        PostQuitMessage(0);
    }

    void timer() {
        if (orbitOn && orbitThumbnailWorker) {
            OrbitThumbnailResult latest;
            if (orbitThumbnailWorker->takeLatest(latest)) {
                orbitThumbnail = std::move(latest.pixels);
                orbitThumbnailBuilding = false;
                orbitThumbnailGeneration = latest.generation;
                fractalOnlyTick = true;
                RECT dirty = rcOrbitThumb;
                dirty.top -= S(20);
                InvalidateRect(hwnd, &dirty, FALSE);
                checkOrbitThumbnailSmoke(latest);
            }
        }
        if (orbitThumbnailSmoke &&
            std::chrono::steady_clock::now() >
                orbitThumbnailSmokeDeadline) {
            finishOrbitThumbnailSmoke(
                false, "thumbnail worker timed out");
            return;
        }
        if (orbitOn && orbitWorker) {
            OrbitResult latest;
            if (orbitWorker->takeLatest(latest)) {
                orbitResult = std::move(latest);
                if (orbitBench) {
                    FILE* f = nullptr; fopen_s(&f, "build\\orbit_bench.txt", "a");
                    if (f) {
                        fprintf(f, "orbit=%.3f ms iterations=%d points=%zu escaped=%d pixel=(%.17g,%.17g)\n",
                                orbitResult.computeMs, orbitResult.iterations,
                                orbitResult.points.size(), orbitResult.escaped ? 1 : 0,
                                orbitResult.pixelRe, orbitResult.pixelIm);
                        fclose(f);
                    }
                    PostQuitMessage(0);
                }
                fractalOnlyTick = true;
                RECT dirty = rcOrbitThumb; dirty.top -= S(20);
                InvalidateRect(hwnd, &dirty, FALSE);
            }
        }
        bool computing = nav->IsComputing();
        if (benchMode && !fpsLogMode && !computing) { runRecolorBench(); return; }
        bool anim = std::fabs(animSpeed) > 1e-4f;
        bool active = computing || wasComputing || navDragging || palette.dragging ||
                      pressed == H_MAXTRACK || pressed == H_DENSTRACK ||
                      pressed == H_PHASETRACK || pressed == H_SPEEDTRACK ||
                      pressed == H_RELAZ || pressed == H_RELEL || pressed == H_RELSTR ||
                      liveFrames > 0 || anim;
        auto now = std::chrono::steady_clock::now();
        if (!active) { lastAnimTick = now; return; }   // idle: no repaint, no CPU spin
        // Frame-rate cap: the timer runs fine (see SetTimer) only to avoid the 16 ms
        // WM_TIMER beat -- where a ~15 ms handler crossing the 16 ms tick boundary was
        // coalesced to the next tick (~31 ms => ~32 fps). Skip ticks that arrive sooner
        // than ~60 fps so a cheap frame doesn't over-render (waste CPU).
        if (lastRenderTs.time_since_epoch().count()) {
            static double capMs = -1;
            if (capMs < 0) { const char* c = getenv("MANDEL_GUI_CAPMS"); capMs = c ? atof(c) : 12.0; }
            double sinceRender = std::chrono::duration<double, std::milli>(now - lastRenderTs).count();
            if (sinceRender < capMs) return;
        }
        lastRenderTs = now;
        if (liveFrames > 0) --liveFrames;
        nav->Update();
        // Advance the palette phase by elapsed wall-clock time so the cycle speed
        // is identical at 30 or 60 fps. animSpeed==1 => one cycle / PHASE_SECONDS_AT_MAX.
        double dt = std::chrono::duration<double>(now - lastAnimTick).count();
        lastAnimTick = now;
        if (anim) {
            double rate = (double)animSpeed / PHASE_SECONDS_AT_MAX * colP;   // index units / s
            color_phase = (float)std::fmod((double)color_phase + rate * dt, (double)colP);
            if (color_phase < 0) color_phase += colP;
        }
        // Cheap cached phase re-colour when we are only cycling the palette on a
        // settled frame; full re-colour while computing / after settle / on other
        // control drags (which also (re)builds the phase cache).
        bool phaseOnly = !computing && !wasComputing && !navDragging &&
                         pressed != H_MAXTRACK && pressed != H_DENSTRACK &&
                         (anim || pressed == H_PHASETRACK || pressed == H_SPEEDTRACK ||
                          pressed == H_RELAZ || pressed == H_RELEL || pressed == H_RELSTR);
        auto _rc0 = std::chrono::steady_clock::now();
        if (phaseOnly) nav->RecolorPhase(bitmap.data());
        else nav->UpdateBitmap(bitmap.data());
        lastRecolorMs = std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - _rc0).count();
        // Smoothed frame rate while animating (measured cadence, not the timer set point).
        if (anim || pressed == H_PHASETRACK || pressed == H_RELAZ ||
            pressed == H_RELEL || pressed == H_RELSTR) {
            if (lastFrameTs.time_since_epoch().count()) {
                double dtf = std::chrono::duration<double>(now - lastFrameTs).count();
                if (dtf > 0) { double f = 1.0 / dtf; animFps = animFps ? animFps * 0.85 + f * 0.15 : f; }
            }
            lastFrameTs = now;
        } else { lastFrameTs = {}; animFps = 0; }
        if (wasComputing && !computing)
            lastRenderMs = std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - renderStart).count();
        wasComputing = computing;
        buildDisplay();
        // Pure fractal/phase-animation frames (pan/zoom/compute/cycle) can skip the
        // panel rebuild; UI-control interaction still gets a full paint.
        bool uiInteract = palette.dragging ||
                          pressed == H_MAXTRACK || pressed == H_DENSTRACK ||
                          pressed == H_PHASETRACK || pressed == H_SPEEDTRACK ||
                          pressed == H_RELAZ || pressed == H_RELEL || pressed == H_RELSTR;
        fractalOnlyTick = !uiInteract;
        InvalidateRect(hwnd, nullptr, FALSE);
        // Force the paint NOW (synchronously) so the logged animFps reflects the real
        // recolor+buildDisplay+present cost per frame, then sample the sustained rate.
        if (fpsLogMode) {
            UpdateWindow(hwnd);   // process the WM_PAINT synchronously (real present)
            if (!computing) {
                auto tnow = std::chrono::steady_clock::now();
                double handlerMs = std::chrono::duration<double, std::milli>(tnow - now).count();
                if (fpsLogStart.time_since_epoch().count() == 0) { fpsLogStart = tnow; fpsLogLast = tnow; }
                double sinceLast = std::chrono::duration<double>(tnow - fpsLogLast).count();
                if (sinceLast >= 0.5) {
                    fpsLogLast = tnow;
                    if (animFps > 0) { fpsLogSum += animFps; ++fpsLogSamples; if (animFps < fpsLogMin) fpsLogMin = animFps; }
                    RECT vr = viewRect();
                    FILE* f = nullptr; fopen_s(&f, "build\\gui_fpslog.txt", "a");
                    if (f) {
                        fprintf(f, "[fpslog] anim=%.1f fps (%.1f ms/frame)  handler=%.2f ms  recolor=%.2f  present=%.2f  view=%ldx%ld ss=%d\n",
                                animFps, animFps > 0 ? 1000.0 / animFps : 0.0, handlerMs, lastRecolorMs, lastPresentMs,
                                vr.right - vr.left, vr.bottom - vr.top, ssOn ? 1 : 0);
                        fclose(f);
                    }
                }
                if (std::chrono::duration<double>(tnow - fpsLogStart).count() >= fpsLogDur) {
                    double avg = fpsLogSamples ? fpsLogSum / fpsLogSamples : 0.0;
                    printf("[fpslog] DONE avg=%.1f fps  min=%.1f fps  samples=%d  render=%dx%d ss=%d\n",
                           avg, (fpsLogMin < 1e8 ? fpsLogMin : 0.0), fpsLogSamples, renderW, renderH, ssOn ? 1 : 0);
                    fflush(stdout);
                    FILE* f = nullptr; fopen_s(&f, "build\\gui_fpslog.txt", "a");
                    if (f) { fprintf(f, "avg=%.1f fps  min=%.1f fps  recolor=%.2f present=%.2f render=%dx%d ss=%d\n",
                                     avg, (fpsLogMin < 1e8 ? fpsLogMin : 0.0), lastRecolorMs, lastPresentMs, renderW, renderH, ssOn ? 1 : 0);
                             fclose(f); }
                    PostQuitMessage(0);
                }
            }
        }
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
            formulaEditor = std::make_unique<FormulaEditorPanel>();
            FormulaEditorCallbacks editorCallbacks;
            editorCallbacks.apply = [this](const FormulaDialogConfig& config) {
                return applyFormulaConfig(config);
            };
            editorCallbacks.useMandelbrot = [this]() {
                restoreMandelbrotUi();
            };
            editorCallbacks.close = [this]() {
                closeFormulaEditor();
            };
            if (!formulaEditor->create(hwnd, dpi, std::move(editorCallbacks)))
                formulaEditor.reset();
            juliaUiEnabled = getenv("MANDEL_EXPERIMENTAL_JULIA") != nullptr ||
                             getenv("MANDEL_GUI_JULIA") != nullptr;
            orbitThumbnailWorker =
                std::make_unique<OrbitThumbnailWorker>();
            orbitThumbnailSmoke =
                getenv("MANDEL_GUI_ORBIT_THUMBNAIL_SMOKE") != nullptr;
            if (const char* e = getenv("MANDEL_GUI_ORBIT")) orbitOn = atoi(e) != 0;
            orbitBench = getenv("MANDEL_GUI_ORBIT_BENCH") != nullptr;
            if (orbitBench || orbitThumbnailSmoke) orbitOn = true;
            if (orbitOn) rebuildOrbitThumbnail();
            const bool anyBench =
                getenv("MANDEL_GUI_BENCH") ||
                orbitBench || orbitThumbnailSmoke;
            if (anyBench) {
                const char* gx = getenv("MANDEL_GUI_CX"), *gy = getenv("MANDEL_GUI_CY");
                const char* gz = getenv("MANDEL_GUI_ZOOM");
                if (gx && gy && gz) {
                    std::string sc = expandSci(gz);
                    if (!sc.empty()) nav->SetLocation(gx, gy, sc);
                }
                if (const char* ww = getenv("MANDEL_GUI_WINW")) {
                    int winw = atoi(ww); const char* wh = getenv("MANDEL_GUI_WINH");
                    int winh = wh ? atoi(wh) : (int)(winw * 0.66);
                    if (winw > 200 && winh > 200)
                        MoveWindow(hwnd, 40, 40, winw, winh, TRUE);
                }
            }
            if (const char* e = getenv("MANDEL_GUI_JULIA"); e && atoi(e)) {
                switchJuliaMode(true, false);
            }
            if (const char* expression = getenv("MANDEL_GUI_FORMULA")) {
                FormulaDialogConfig config = formulaConfig;
                config.source = expression;
                const char* plane = getenv("MANDEL_GUI_FORMULA_PIXEL");
                config.pixelParameter = plane && strcmp(plane, "z0") == 0
                    ? FormulaParameter::InitialZ : FormulaParameter::C;
                if (const char* value = getenv("MANDEL_GUI_FIXED_Z0_RE")) config.fixedZ0.real(atof(value));
                if (const char* value = getenv("MANDEL_GUI_FIXED_Z0_IM")) config.fixedZ0.imag(atof(value));
                if (const char* value = getenv("MANDEL_GUI_FIXED_C_RE")) config.fixedC.real(atof(value));
                if (const char* value = getenv("MANDEL_GUI_FIXED_C_IM")) config.fixedC.imag(atof(value));
                if (const char* value = getenv("MANDEL_GUI_FORMULA_BAILOUT")) config.bailout = atof(value);
                applyFormulaConfig(config, false);
                // z0-plane activation intentionally resets its view; apply the
                // requested benchmark location after selecting the binding.
                if (anyBench) {
                    const char* gx = getenv("MANDEL_GUI_CX"), *gy = getenv("MANDEL_GUI_CY");
                    const char* gz = getenv("MANDEL_GUI_ZOOM");
                    if (gx && gy && gz) {
                        std::string sc = expandSci(gz);
                        if (!sc.empty()) nav->SetLocation(gx, gy, sc);
                    }
                }
                if (getenv("MANDEL_GUI_FORMULA_RESTORE"))
                    restoreMandelbrotUi(false);
                if (getenv("MANDEL_GUI_FORMULA_JULIA_RESTORE")) {
                    switchJuliaMode(true, false);
                    switchJuliaMode(false, false);
                }
            }
            if (orbitThumbnailSmoke)
                startOrbitThumbnailSmoke();
            if (getenv("MANDEL_GUI_BENCH")) {
                benchMode = true;
                const char* e = getenv("MANDEL_GUI_SS");
                ssOn = nav->IsMandelbrot() && (!e || atoi(e));
                if (ssOn) nav->SetCMethod(nav->GetCMethod() | ColoringMethod::SUPER_SAMPLING);
                if (const char* cc = getenv("MANDEL_GUI_COLOR")) selectColoring(atoi(cc));
                if (const char* ds = getenv("MANDEL_GUI_DENS")) color_density = (float)atof(ds);
                if (const char* pp = getenv("MANDEL_GUI_PAL")) { paletteIdx = atoi(pp); palette.load(paletteIdx); }
                // Apply a gallery preset (e.g. 0 = Night City) for realistic testing.
                if (const char* gi = getenv("MANDEL_GUI_GALLERY")) applyPreset(atoi(gi));
                // Force the palette dropdown open (verify it scrolls fully into view).
                if (getenv("MANDEL_GUI_OPENPAL")) {
                    paletteOpen = true; layout();
                    int over = (int)paletteListRect().bottom - panelViewH + S(8); if (over > 0) panelScrollBy(over);
                }
                // Real sustained-fps logging mode: run the actual timer/animation loop
                // and log the measured animFps, instead of the one-shot micro-bench.
                if (const char* fl = getenv("MANDEL_GUI_FPSLOG")) {
                    fpsLogMode = true; fpsLogDur = atof(fl); if (fpsLogDur <= 0) fpsLogDur = 5.0;
                    const char* as = getenv("MANDEL_GUI_ANIMSPEED");
                    animSpeed = as ? (float)atof(as) : 0.6f;
                    lastAnimTick = std::chrono::steady_clock::now();
                }
            }
            layout();
            if (orbitBench || orbitThumbnailSmoke) {
                // Worker smoke/bench runs must not compete with a full frame render.
            } else if (benchMode) {
                // A hidden benchmark window receives no initial WM_SIZE, so size the
                // navigator explicitly instead of silently rendering at 900x600.
                int oldW = renderW, oldH = renderH;
                retargetToView();
                if (renderW == oldW && renderH == oldH) startRender();
            } else {
                startRender();
            }
            if (orbitOn && !orbitThumbnailSmoke) {
                RECT vr = viewRect();
                requestOrbit((vr.left + vr.right) / 2, (vr.top + vr.bottom) / 2);
            }
            // Fine timer period (~4 ms) so the ~60 fps redraw cap in timer() paces the
            // animation instead of the coarse 16 ms WM_TIMER, which coalesced a ~15 ms
            // handler to the next tick (~31 ms => ~32 fps). The cap prevents over-render.
            { int tms = 4; if (const char* t = getenv("MANDEL_GUI_TIMERMS")) { tms = atoi(t); if (tms < 1) tms = 1; }
              SetTimer(hwnd, TIMER_ID, tms, nullptr); }
            if (getenv("MANDEL_GUI_OPEN_FORMULA"))
                PostMessageW(hwnd, WM_OPEN_FORMULA_TEST, 0, 0);
            return 0;
        }
        case WM_OPEN_FORMULA_TEST:
            showFormulaEditor();
            return 0;
        case WM_DPICHANGED: {
            dpi = HIWORD(wp);
            createFonts();
            if (formulaEditor) formulaEditor->setDpi(dpi);
            RECT* sug = (RECT*)lp;
            SetWindowPos(hwnd, nullptr, sug->left, sug->top,
                         sug->right - sug->left, sug->bottom - sug->top,
                         SWP_NOZORDER | SWP_NOACTIVATE);
            layout(); InvalidateRect(hwnd, nullptr, FALSE);
            return 0;
        }
        case WM_ENTERSIZEMOVE: inSizeMove = true; return 0;
        case WM_EXITSIZEMOVE:
            // Interactive resize finished -> retarget the render buffers to the new
            // view size (native resolution) and re-render once.
            inSizeMove = false; retargetToView(); layout();
            InvalidateRect(hwnd, nullptr, TRUE); return 0;
        case WM_SIZE:
            layout();
            // Maximize/restore/DPI changes arrive without ENTERSIZEMOVE; retarget them
            // immediately. Interactive border drags defer to WM_EXITSIZEMOVE (above),
            // painting a transient scaled preview until then.
            if (wp != SIZE_MINIMIZED && !inSizeMove && !orbitBench) retargetToView();
            InvalidateRect(hwnd, nullptr, FALSE); return 0;
        case WM_ERASEBKGND: return 1;
        case WM_PAINT: paint(); return 0;
        case WM_TIMER: if (wp == TIMER_ID) timer(); return 0;
        case WM_MOUSEMOVE: {
            int x = GET_X_LPARAM(lp), y = GET_Y_LPARAM(lp);
            if (palette.dragging) { gradientMove(x); return 0; }
            if (pressed == H_PANELSB) {
                int maxs = panelMaxScroll();
                RECT rc; GetClientRect(hwnd, &rc);
                int trackH = std::max(1, (int)(rc.bottom - 2 * S(2)));
                int thumbH = std::max(S(28), (int)((double)trackH * panelViewH / std::max(1, panelContentH)));
                int span = std::max(1, trackH - thumbH);
                int ns = panelSbGrabScroll + (int)((double)(y - panelSbGrabY) * maxs / span);
                panelScroll = std::clamp(ns, 0, maxs);
                layout(); needFull = true; InvalidateRect(hwnd, nullptr, FALSE);
                return 0;
            }
            if (pressed == H_MAXTRACK) { setMaxFromT((double)(x - rcMaxTrack.left)/(rcMaxTrack.right-rcMaxTrack.left), false); InvalidateRect(hwnd,nullptr,FALSE); return 0; }
            if (pressed == H_DENSTRACK) { color_density = (float)std::round(std::clamp(10.0 + 190.0*(x-rcDensTrack.left)/(rcDensTrack.right-rcDensTrack.left),10.0,200.0)); nav->SetRedisplay(); InvalidateRect(hwnd,nullptr,FALSE); return 0; }
            if (pressed == H_PHASETRACK) { setPhaseFromX(x); return 0; }
            if (pressed == H_SPEEDTRACK) { setSpeedFromX(x); InvalidateRect(hwnd,nullptr,FALSE); return 0; }
            if (pressed == H_RELAZ) { setReliefAzFromX(x); return 0; }
            if (pressed == H_RELEL) { setReliefElFromX(x); return 0; }
            if (pressed == H_RELSTR) { setReliefStrFromX(x); return 0; }
            if (navDragging) { POINT p = mapToRender(x,y); nav->Drag(p.x,p.y); return 0; }
            if (orbitOn && inRect(viewRect(), x, y)) requestOrbit(x, y);
            if (paletteOpen) { int it = paletteItemAt(x,y); if (it != paletteHover) { paletteHover = it; InvalidateRect(hwnd,nullptr,FALSE); } }
            if (coloringOpen) { int it = coloringItemAt(x,y); if (it != coloringHover) { coloringHover = it; InvalidateRect(hwnd,nullptr,FALSE); } }
            if (galleryOpen) { int it = galleryItemAt(x,y); if (it != galleryHover) { galleryHover = it; InvalidateRect(hwnd,nullptr,FALSE); } }
            Hit h = hitTest(x,y);
            if (h != hover) { hover = h; InvalidateRect(hwnd, nullptr, FALSE); }
            return 0;
        }
        case WM_LBUTTONDOWN: {
            int x = GET_X_LPARAM(lp), y = GET_Y_LPARAM(lp);
            if (formulaEditor && formulaEditor->visible())
                formulaEditor->dismissPopups();
            gradientFocused = false;
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
            gradientFocused = h == H_GRADIENT;
            if (h == H_PALETTE_DD) { paletteOpen = true; paletteHover = -1; pressed = H_NONE; layout(); int over = (int)paletteListRect().bottom - panelViewH + S(8); if (over > 0) panelScrollBy(over); InvalidateRect(hwnd, nullptr, FALSE); return 0; }
            if (h == H_EDE) { coloringOpen = true; coloringHover = -1; pressed = H_NONE; layout(); int over = (int)coloringListRect().bottom - panelViewH + S(8); if (over > 0) panelScrollBy(over); InvalidateRect(hwnd, nullptr, FALSE); return 0; }
            if (h == H_GALLERY_DD) { galleryOpen = true; galleryHover = -1; pressed = H_NONE; layout(); int over = (int)galleryListRect().bottom - panelViewH + S(8); if (over > 0) panelScrollBy(over); InvalidateRect(hwnd, nullptr, FALSE); return 0; }
            pressed = h;
            if (h == H_PANELSB) { SetCapture(hwnd); panelSbDrag = true; panelSbGrabY = y; panelSbGrabScroll = panelScroll; return 0; }
            if (h == H_GRADIENT) { gradientDown(x); return 0; }
            if (h == H_MAXTRACK) { SetCapture(hwnd); setMaxFromT((double)(x-rcMaxTrack.left)/(rcMaxTrack.right-rcMaxTrack.left), false); InvalidateRect(hwnd,nullptr,FALSE); return 0; }
            if (h == H_DENSTRACK) { SetCapture(hwnd); color_density=(float)std::round(std::clamp(10.0+190.0*(x-rcDensTrack.left)/(rcDensTrack.right-rcDensTrack.left),10.0,200.0)); nav->SetRedisplay(); InvalidateRect(hwnd,nullptr,FALSE); return 0; }
            if (h == H_PHASETRACK) { SetCapture(hwnd); setPhaseFromX(x); return 0; }
            if (h == H_SPEEDTRACK) { SetCapture(hwnd); setSpeedFromX(x); InvalidateRect(hwnd,nullptr,FALSE); return 0; }
            if (h == H_RELAZ) { SetCapture(hwnd); setReliefAzFromX(x); return 0; }
            if (h == H_RELEL) { SetCapture(hwnd); setReliefElFromX(x); return 0; }
            if (h == H_RELSTR) { SetCapture(hwnd); setReliefStrFromX(x); return 0; }
            if (h == H_VIEW) { POINT p = mapToRender(x,y); navDragging = true; SetCapture(hwnd); nav->DragStart(p.x,p.y); }
            InvalidateRect(hwnd, nullptr, FALSE);
            return 0;
        }
        case WM_LBUTTONUP: {
            int x = GET_X_LPARAM(lp), y = GET_Y_LPARAM(lp);
            Hit h = hitTest(x,y);
            if (palette.dragging) { palette.dragging = false; ReleaseCapture(); }
            else if (panelSbDrag) { panelSbDrag = false; ReleaseCapture(); }
            else if (navDragging) { navDragging = false; nav->DragEnd(); ReleaseCapture(); }
            else if (pressed == H_MAXTRACK || pressed == H_DENSTRACK) {
                ReleaseCapture();
                if (pressed == H_MAXTRACK) startRender(); else startRender();
            } else if (pressed == H_PHASETRACK || pressed == H_SPEEDTRACK ||
                       pressed == H_RELAZ || pressed == H_RELEL || pressed == H_RELSTR) {
                ReleaseCapture();   // phase/speed/relief re-colour only; no fractal recompute
            } else if (h == pressed && h != H_NONE) {
                switch (h) {
                case H_RESET: nav->Reset(); needFull = true; renderStart = std::chrono::steady_clock::now(); wasComputing = true; keepLive(); break;
                case H_RENDER: startRender(); break;
                case H_SAVE:
                    if (!nav->IsMandelbrot()) saveImage();
                    else showExportDialog(hwnd, nav.get());
                    break;
                case H_COPY: copyLocation(); break;
                case H_PASTE: pasteLocation(); break;
                case H_SS:
                    if (nav->IsMandelbrot()) {
                        ssOn = !ssOn;
                        setMethodFlag(ColoringMethod::SUPER_SAMPLING, ssOn);
                    }
                    break;
                case H_FORMULA:
                    showFormulaEditor();
                    break;
                case H_JULIA: {
                    bool enable = !nav->IsJulia();
                    switchJuliaMode(enable);
                    break;
                }
                case H_ORBIT:
                    if (nav->SupportsOrbitOverlay()) orbitOn = !orbitOn;
                    if (orbitOn) {
                        rebuildOrbitThumbnail();
                    } else {
                        if (orbitWorker) orbitWorker->cancel();
                        if (orbitThumbnailWorker)
                            orbitThumbnailWorker->cancel();
                        orbitThumbnail.clear();
                        orbitThumbnailBuilding = false;
                    }
                    orbitResult = OrbitResult{};
                    layout(); needFull = true;
                    if (orbitOn) requestCurrentOrbit();
                    break;
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
            if (h == H_GRADIENT) {
                gradientFocused = true;
                gradientDelete(x);
            }
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
            // Phase wheel: nudge the palette phase (~1/256 cycle per notch), no snap.
            RECT pt = rcPhaseTrack; pt.top -= S(10); pt.bottom += S(10);
            if (inRect(pt, q.x, q.y)) {
                color_phase += (wd > 0 ? 1.0f : -1.0f) * (colP / 256.0f);
                recolorPhaseNow();
                return 0;
            }
            // Speed wheel: fine-tune the cycle speed (bypasses the centre snap so a
            // tiny non-zero speed is reachable); clamped to [-1, 1].
            RECT sp = rcSpeedTrack; sp.top -= S(10); sp.bottom += S(10);
            if (inRect(sp, q.x, q.y)) {
                animSpeed = (float)std::clamp(animSpeed + (wd > 0 ? 0.02 : -0.02), -1.0, 1.0);
                lastAnimTick = std::chrono::steady_clock::now();
                InvalidateRect(hwnd, nullptr, FALSE);
                return 0;
            }
            if (relief_on || normal_light_on || de_overlay_on) {
                RECT ra = rcReliefAz; ra.top -= S(10); ra.bottom += S(10);
                if (inRect(ra, q.x, q.y)) {
                    if (de_overlay_on) {
                        de_k = (float)std::clamp((double)de_k + (wd > 0 ? 0.1 : -0.1), DE_K_MIN, DE_K_MAX);
                    } else {
                        const double TWO_PI = 6.28318530717959;
                        double v = relief_light_az + (wd > 0 ? 0.05236 : -0.05236);   // ~3deg/notch
                        v = std::fmod(v, TWO_PI); if (v < 0) v += TWO_PI;
                        relief_light_az = (float)v;
                    }
                    recolorPhaseNow(); return 0;
                }
                RECT re = rcReliefEl; re.top -= S(10); re.bottom += S(10);
                if (inRect(re, q.x, q.y)) {
                    relief_light_el = (float)std::clamp((double)relief_light_el + (wd > 0 ? 0.03 : -0.03), 0.05, 1.5);
                    recolorPhaseNow(); return 0;
                }
                RECT rs = rcReliefStr; rs.top -= S(10); rs.bottom += S(10);
                if (inRect(rs, q.x, q.y)) {
                    if (de_overlay_on) {
                        de_scale = (float)std::clamp((double)de_scale + (wd > 0 ? 0.05 : -0.05), DE_SCALE_MIN, DE_SCALE_MAX);
                    } else {
                        relief_strength = (float)std::clamp((double)relief_strength + (wd > 0 ? 0.05 : -0.05), 0.0, RELIEF_STR_MAX);
                    }
                    recolorPhaseNow(); return 0;
                }
            }
            if (inRect(viewRect(), q.x, q.y)) {
                POINT p = mapToRender(q.x, q.y);
                if (wd > 0) nav->ZoomIn(p.x,p.y); else nav->ZoomOut(p.x,p.y);
                keepLive();
                return 0;
            }
            // Anywhere else over the panel: scroll the control stack (if it overflows).
            RECT rcC; GetClientRect(hwnd, &rcC);
            int panelR = panelRight();
            if (q.x >= panelR - S(PANEL_W) && q.x < panelR &&
                panelMaxScroll() > 0)
                panelScrollBy(wd > 0 ? -S(48) : S(48));
            return 0;
        }
        case WM_KEYDOWN:
            if (wp == 'R') { nav->Reset(); needFull = true; renderStart = std::chrono::steady_clock::now(); wasComputing = true; keepLive(); }
            else if (wp == VK_SPACE) startRender();
            else if (wp == 'S') saveImage();
            else if (wp == 'C') copyLocation();
            else if (wp == 'V') pasteLocation();
            else if ((wp == VK_DELETE || wp == VK_BACK) &&
                     gradientFocused)
                gradientDeleteSelected();
            else if (wp == VK_ADD || wp == VK_OEM_PLUS) { nav->ZoomIn(renderW/2, renderH/2); keepLive(); }
            else if (wp == VK_SUBTRACT || wp == VK_OEM_MINUS) { nav->ZoomOut(renderW/2, renderH/2); keepLive(); }
            return 0;
        case WM_GETMINMAXINFO: {
            auto* m = (MINMAXINFO*)lp;
            m->ptMinTrackSize.x =
                S(formulaEditor && formulaEditor->visible() ? 1100 : 900);
            m->ptMinTrackSize.y = S(640);
            return 0;
        }
        case WM_DESTROY:
            KillTimer(hwnd, TIMER_ID);
            orbitThumbnailWorker.reset();
            orbitWorker.reset();
            formulaEditor.reset();
            nav.reset();
            if (memDC) { SelectObject(memDC, memOld); DeleteObject(memBmp); DeleteDC(memDC); memDC = nullptr; }
            DeleteObject(fUi); DeleteObject(fBold); DeleteObject(fSmall);
            DeleteObject(fMono); DeleteObject(fMicro);
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
    // Raise the system timer resolution to 1 ms so the 16 ms animation WM_TIMER
    // actually fires at ~60 Hz. At the default ~15.6 ms resolution a 16 ms timer
    // is coalesced to the next tick (~31 ms => ~32 fps, with jitter).
    timeBeginPeriod(1);
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
    const bool benchHeadless =
        (getenv("MANDEL_GUI_BENCH") ||
         getenv("MANDEL_GUI_ORBIT_BENCH") ||
         getenv("MANDEL_GUI_ORBIT_THUMBNAIL_SMOKE")) &&
        !getenv("MANDEL_GUI_SHOW");
    HWND hwnd = CreateWindowExW(0, wc.lpszClassName, L"Mandelbrot Explorer",
        WS_OVERLAPPEDWINDOW | WS_CLIPCHILDREN | (benchHeadless ? 0 : WS_VISIBLE),
        wx, wy, winW, winH,
        nullptr, nullptr, instance, &app);
    if (!hwnd) return 1;
    ShowWindow(hwnd, benchHeadless ? SW_HIDE : show); UpdateWindow(hwnd);
    if (!benchHeadless && getenv("MANDEL_GUI_MAX")) ShowWindow(hwnd, SW_MAXIMIZE);

    MSG msg;
    while (GetMessageW(&msg, nullptr, 0, 0) > 0) {
        TranslateMessage(&msg);
        DispatchMessageW(&msg);
    }
    g_app = nullptr;
    timeEndPeriod(1);
    return (int)msg.wParam;
}
