#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include "formula_editor_panel.h"

#include <windowsx.h>

#include <algorithm>
#include <array>
#include <cerrno>
#include <cmath>
#include <complex>
#include <cstring>
#include <cwchar>
#include <cwctype>
#include <string>
#include <utility>
#include <vector>

#include "formula_expression.h"
#include "gui_theme.h"
#include "ui_framework.h"

namespace {

constexpr wchar_t PANEL_CLASS[] = L"MandelFormulaEditorPanel";
constexpr UINT WM_FORMULA_CARET_CHANGED = WM_APP + 73;

constexpr int ID_SOURCE = 1001;
constexpr int ID_PRESET = 1002;
constexpr int ID_BAILOUT = 1003;
constexpr int ID_REAL = 1004;
constexpr int ID_IMAGINARY = 1005;

constexpr int ACTION_HEIGHT = 84;
constexpr int CONTENT_HEIGHT = 674;

constexpr int HIT_NONE = 0;
constexpr int HIT_COPY = 1;
constexpr int HIT_PASTE = 2;
constexpr int HIT_C_PLANE = 3;
constexpr int HIT_Z0_PLANE = 4;
constexpr int HIT_RANGE_MINUS = 5;
constexpr int HIT_RANGE_PLUS = 6;
constexpr int HIT_RANGE_RESET = 7;
constexpr int HIT_REVERT = 8;
constexpr int HIT_MANDELBROT = 9;
constexpr int HIT_APPLY = 10;
constexpr int HIT_CLOSE = 11;
constexpr int HIT_PICKER = 12;
constexpr int HIT_VARIABLE_BASE = 100;
constexpr int HIT_FUNCTION_BASE = 200;

enum class InspectorValue {
    Z,
    C,
    Z0,
    N,
    P0,
    P1,
    P2,
    P3,
    P4,
    P5,
    P6,
    P7
};

struct Preset {
    const wchar_t* label;
    const char* source;
};

const std::array<Preset, 8>& presets() {
    static const std::array<Preset, 8> values = {{
        { L"Custom", nullptr },
        { L"Quadratic z*z+c", "z*z+c" },
        { L"Cubic z*z*z+c", "z*z*z+c" },
        { L"Burning Ship sqr(complex(abs(real(z)),abs(imag(z))))+c",
          "sqr(complex(abs(real(z)),abs(imag(z))))+c" },
        { L"Sine sin(z)+c", "sin(z)+c" },
        { L"Parameter polynomial z*z+c+p0*z", "z*z+c+p0*z" },
        { L"Iteration drift z*z+c+0.0001*n", "z*z+c+0.0001*n" },
        { L"Branch power exp(p0*log(z))+c", "exp(p0*log(z))+c" }
    }};
    return values;
}

const std::array<const wchar_t*, 12>& variableLabels() {
    static const std::array<const wchar_t*, 12> labels = {
        L"z", L"c", L"z0", L"n",
        L"p0", L"p1", L"p2", L"p3", L"p4", L"p5", L"p6", L"p7"
    };
    return labels;
}

const std::array<const wchar_t*, 16>& functionLabels() {
    static const std::array<const wchar_t*, 16> labels = {
        L"sqr", L"^", L"sqrt", L"abs", L"sin", L"cos", L"tan", L"exp",
        L"log", L"conj", L"real", L"imag", L"norm", L"arg", L"polar", L"complex"
    };
    return labels;
}

bool widenUtf8(const std::string& value, std::wstring& result) {
    result.clear();
    if (value.empty()) return true;
    int count = MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS,
                                    value.data(), static_cast<int>(value.size()),
                                    nullptr, 0);
    if (count <= 0) return false;
    result.resize(static_cast<size_t>(count));
    return MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS,
                               value.data(), static_cast<int>(value.size()),
                               result.data(), count) == count;
}

bool narrowUtf8(const std::wstring& value, std::string& result) {
    result.clear();
    if (value.empty()) return true;
    int count = WideCharToMultiByte(CP_UTF8, WC_ERR_INVALID_CHARS,
                                    value.data(), static_cast<int>(value.size()),
                                    nullptr, 0, nullptr, nullptr);
    if (count <= 0) return false;
    result.resize(static_cast<size_t>(count));
    return WideCharToMultiByte(CP_UTF8, WC_ERR_INVALID_CHARS,
                               value.data(), static_cast<int>(value.size()),
                               result.data(), count, nullptr, nullptr) == count;
}

std::wstring getWindowText(HWND control) {
    int count = GetWindowTextLengthW(control);
    if (count <= 0) return {};
    std::wstring result(static_cast<size_t>(count) + 1, L'\0');
    int copied = GetWindowTextW(control, result.data(), count + 1);
    result.resize(static_cast<size_t>(std::max(0, copied)));
    return result;
}

bool parseFiniteNumber(HWND control, double& value) {
    std::wstring text = getWindowText(control);
    const wchar_t* begin = text.c_str();
    while (*begin && std::iswspace(*begin)) ++begin;
    if (!*begin) return false;

    errno = 0;
    wchar_t* end = nullptr;
    double parsed = std::wcstod(begin, &end);
    while (end && *end && std::iswspace(*end)) ++end;
    if (errno == ERANGE || end == begin || !end || *end != L'\0' ||
        !std::isfinite(parsed)) {
        return false;
    }
    value = parsed;
    return true;
}

void setNumber(HWND control, double value) {
    wchar_t buffer[80];
    swprintf_s(buffer, L"%.17g", value);
    SetWindowTextW(control, buffer);
}

bool containsPoint(const RECT& rect, int x, int y) {
    return x >= rect.left && x < rect.right && y >= rect.top && y < rect.bottom;
}

bool finiteComplex(const std::complex<double>& value) {
    return std::isfinite(value.real()) && std::isfinite(value.imag());
}

InspectorValue variableFromButton(size_t index) {
    static const std::array<InspectorValue, 12> values = {
        InspectorValue::Z, InspectorValue::C, InspectorValue::Z0, InspectorValue::N,
        InspectorValue::P0, InspectorValue::P1, InspectorValue::P2, InspectorValue::P3,
        InspectorValue::P4, InspectorValue::P5, InspectorValue::P6, InspectorValue::P7
    };
    return values[index];
}

int parameterIndex(InspectorValue value) {
    int raw = static_cast<int>(value) - static_cast<int>(InspectorValue::P0);
    return raw >= 0 && raw < 8 ? raw : -1;
}

const wchar_t* inspectorName(InspectorValue value) {
    switch (value) {
    case InspectorValue::Z: return L"z";
    case InspectorValue::C: return L"c";
    case InspectorValue::Z0: return L"z0";
    case InspectorValue::N: return L"n";
    case InspectorValue::P0: return L"p0";
    case InspectorValue::P1: return L"p1";
    case InspectorValue::P2: return L"p2";
    case InspectorValue::P3: return L"p3";
    case InspectorValue::P4: return L"p4";
    case InspectorValue::P5: return L"p5";
    case InspectorValue::P6: return L"p6";
    case InspectorValue::P7: return L"p7";
    }
    return L"";
}

struct ButtonSpec {
    RECT rect{};
    int hit = HIT_NONE;
    std::wstring label;
};

} // namespace

struct FormulaEditorPanel::Impl {
    HWND hwnd = nullptr;
    HWND owner = nullptr;
    HWND sourceEdit = nullptr;
    HWND presetCombo = nullptr;
    HWND bailoutEdit = nullptr;
    HWND realEdit = nullptr;
    HWND imaginaryEdit = nullptr;
    WNDPROC sourceProc = nullptr;

    FormulaEditorCallbacks callbacks;
    FormulaDialogConfig working;
    FormulaDialogConfig applied;
    InspectorValue selected = InspectorValue::Z0;

    int dpi = 96;
    ui::Resources resources;
    ui::BackBuffer backBuffer;
    ui::ScrollState scroll;
    double pickerRange = 2.0;
    bool syncing = false;
    bool draggingPoint = false;
    bool trackingMouse = false;
    int hoverHit = HIT_NONE;
    int pressedHit = HIT_NONE;

    std::wstring status;
    bool statusError = false;

    HFONT uiFont = nullptr;
    HFONT boldFont = nullptr;
    HFONT smallFont = nullptr;
    HFONT monoFont = nullptr;
    HBRUSH panelBrush = nullptr;
    HBRUSH cardBrush = nullptr;

    ~Impl() {
        destroy();
        deleteGdiObjects();
    }

    int scale(int value) const {
        return MulDiv(value, dpi > 0 ? dpi : 96, 96);
    }

    int unscale(int value) const {
        return MulDiv(value, 96, dpi > 0 ? dpi : 96);
    }

    int clientWidthDip() const {
        if (!hwnd) return FormulaEditorPanel::DESIGN_WIDTH;
        RECT rect{};
        GetClientRect(hwnd, &rect);
        return std::max(1, unscale(rect.right - rect.left));
    }

    int clientHeightDip() const {
        if (!hwnd) return 1;
        RECT rect{};
        GetClientRect(hwnd, &rect);
        return std::max(1, unscale(rect.bottom - rect.top));
    }

    int actionTopDip() const {
        return std::max(0, clientHeightDip() - ACTION_HEIGHT);
    }

    RECT sourceRectDip() const {
        return { 18, 60, std::max(100, clientWidthDip() - 18), 90 };
    }

    RECT presetRectDip() const {
        int right = std::max(180, clientWidthDip() - 18);
        int copyLeft = right - 72 - 8 - 72;
        return { 82, 96, std::max(154, copyLeft - 10), 126 };
    }

    RECT bailoutRectDip() const {
        int right = std::max(280, clientWidthDip() - 18);
        int left = std::max(280, right - 242);
        return { left, 316, right, 346 };
    }

    RECT planeRectDip() const {
        int available = clientWidthDip() - 18 - 17 - 190 - 18;
        int vertical = actionTopDip() - 395 - 10;
        int size = std::clamp(std::min(available, vertical), 220, 300);
        return { 18, 395, 18 + size, 395 + size };
    }

    int inspectorLeftDip() const {
        return planeRectDip().right + 17;
    }

    RECT realRectDip() const {
        return { inspectorLeftDip(), 503, std::max(inspectorLeftDip() + 80,
                  clientWidthDip() - 18), 533 };
    }

    RECT imaginaryRectDip() const {
        return { inspectorLeftDip(), 563, std::max(inspectorLeftDip() + 80,
                  clientWidthDip() - 18), 593 };
    }

    RECT toPixelRect(RECT rect, bool content) const {
        if (content) {
            rect.top -= scroll.position();
            rect.bottom -= scroll.position();
        }
        return {
            scale(rect.left), scale(rect.top),
            scale(rect.right), scale(rect.bottom)
        };
    }

    std::vector<ButtonSpec> contentButtons() const {
        std::vector<ButtonSpec> buttons;
        int right = std::max(180, clientWidthDip() - 18);

        buttons.push_back({ { right - 152, 96, right - 80, 126 }, HIT_COPY, L"Copy" });
        buttons.push_back({ { right - 72, 96, right, 126 }, HIT_PASTE, L"Paste" });

        const auto& variables = variableLabels();
        int left = 18;
        int gap = 4;
        int available = std::max(12, right - left - gap * 11);
        int x = left;
        for (size_t i = 0; i < variables.size(); ++i) {
            int remaining = static_cast<int>(variables.size() - i);
            int width = available / remaining;
            available -= width;
            buttons.push_back({
                { x, 158, x + width, 186 },
                HIT_VARIABLE_BASE + static_cast<int>(i), variables[i]
            });
            x += width + gap;
        }

        const auto& functions = functionLabels();
        for (int row = 0; row < 2; ++row) {
            gap = 5;
            available = std::max(8, right - left - gap * 7);
            x = left;
            for (int column = 0; column < 8; ++column) {
                int remaining = 8 - column;
                int width = available / remaining;
                available -= width;
                int index = row * 8 + column;
                int top = 216 + row * 34;
                buttons.push_back({
                    { x, top, x + width, top + 28 },
                    HIT_FUNCTION_BASE + index, functions[static_cast<size_t>(index)]
                });
                x += width + gap;
            }
        }

        buttons.push_back({ { 18, 316, 138, 346 }, HIT_C_PLANE, L"c-plane" });
        buttons.push_back({ { 138, 316, 258, 346 }, HIT_Z0_PLANE, L"z0-plane" });

        int inspectorLeft = inspectorLeftDip();
        int inspectorRight = std::max(inspectorLeft + 100, right);
        int rangeWidth = inspectorRight - inspectorLeft;
        int small = std::min(42, std::max(28, (rangeWidth - 12) / 4));
        buttons.push_back({
            { inspectorLeft, 635, inspectorLeft + small, 665 },
            HIT_RANGE_MINUS, L"-"
        });
        buttons.push_back({
            { inspectorLeft + small + 6, 635, inspectorLeft + small * 2 + 6, 665 },
            HIT_RANGE_PLUS, L"+"
        });
        buttons.push_back({
            { inspectorLeft + small * 2 + 12, 635, inspectorRight, 665 },
            HIT_RANGE_RESET, L"Reset"
        });
        return buttons;
    }

    std::vector<ButtonSpec> actionButtons() const {
        std::vector<ButtonSpec> buttons;
        int left = 18;
        int right = std::max(left + 200, clientWidthDip() - 18);
        int gap = 8;
        int total = std::max(4, right - left - gap * 3);
        const std::array<int, 4> weights = { 20, 28, 23, 29 };
        const std::array<int, 4> hits = {
            HIT_REVERT, HIT_MANDELBROT, HIT_APPLY, HIT_CLOSE
        };
        const std::array<const wchar_t*, 4> labels = {
            L"Revert", L"Mandelbrot", L"Apply", L"Close"
        };
        int x = left;
        int remainingWeight = 100;
        int remainingWidth = total;
        int top = actionTopDip() + 32;
        for (size_t i = 0; i < weights.size(); ++i) {
            int width = i + 1 == weights.size()
                ? remainingWidth
                : MulDiv(remainingWidth, weights[i], remainingWeight);
            buttons.push_back({
                { x, top, x + width, top + 38 }, hits[i], labels[i]
            });
            x += width + gap;
            remainingWidth -= width;
            remainingWeight -= weights[i];
        }
        return buttons;
    }

    bool isButtonActive(int hit) const {
        if (hit == HIT_C_PLANE)
            return working.pixelParameter != FormulaParameter::InitialZ;
        if (hit == HIT_Z0_PLANE)
            return working.pixelParameter == FormulaParameter::InitialZ;
        if (hit >= HIT_VARIABLE_BASE &&
            hit < HIT_VARIABLE_BASE + static_cast<int>(variableLabels().size())) {
            return selected == variableFromButton(
                static_cast<size_t>(hit - HIT_VARIABLE_BASE));
        }
        return hit == HIT_APPLY;
    }

    int hitTest(int xDip, int yDip) const {
        if (yDip >= actionTopDip()) {
            for (const ButtonSpec& button : actionButtons()) {
                if (containsPoint(button.rect, xDip, yDip)) return button.hit;
            }
            return HIT_NONE;
        }

        int contentY = yDip + scroll.position();
        for (const ButtonSpec& button : contentButtons()) {
            if (containsPoint(button.rect, xDip, contentY)) return button.hit;
        }
        if (containsPoint(planeRectDip(), xDip, contentY)) return HIT_PICKER;
        return HIT_NONE;
    }

    bool selectedEditable(const FormulaDialogConfig& config) const {
        int parameter = parameterIndex(selected);
        if (parameter >= 0) return true;
        if (selected == InspectorValue::C)
            return config.pixelParameter == FormulaParameter::InitialZ;
        if (selected == InspectorValue::Z0)
            return config.pixelParameter != FormulaParameter::InitialZ;
        return false;
    }

    std::complex<double>* selectedValue(FormulaDialogConfig& config) const {
        int parameter = parameterIndex(selected);
        if (parameter >= 0) return &config.parameters[static_cast<size_t>(parameter)];
        if (selected == InspectorValue::C) return &config.fixedC;
        if (selected == InspectorValue::Z0) return &config.fixedZ0;
        return nullptr;
    }

    const std::complex<double>* selectedValue(const FormulaDialogConfig& config) const {
        int parameter = parameterIndex(selected);
        if (parameter >= 0) return &config.parameters[static_cast<size_t>(parameter)];
        if (selected == InspectorValue::C) return &config.fixedC;
        if (selected == InspectorValue::Z0) return &config.fixedZ0;
        return nullptr;
    }

    std::wstring inspectorDescription() const {
        int parameter = parameterIndex(selected);
        if (parameter >= 0)
            return L"Editable complex formula parameter.";
        if (selected == InspectorValue::C) {
            return selectedEditable(working)
                ? L"Fixed c used while z0 is supplied by each pixel."
                : L"c is supplied by each pixel in c-plane and is read-only here.";
        }
        if (selected == InspectorValue::Z0) {
            return selectedEditable(working)
                ? L"Fixed initial z used while c is supplied by each pixel."
                : L"z0 is supplied by each pixel in z0-plane and is read-only here.";
        }
        if (selected == InspectorValue::Z)
            return L"z is the current iteration value and cannot be edited.";
        return L"n is the iteration index and cannot be edited.";
    }

    void setStatus(std::wstring text, bool error) {
        status = std::move(text);
        statusError = error;
        if (hwnd) InvalidateRect(hwnd, nullptr, FALSE);
    }

    void deleteGdiObjects() {
        resources.reset();
        uiFont = boldFont = smallFont = monoFont = nullptr;
        panelBrush = cardBrush = nullptr;
    }

    bool recreateFonts() {
        if (!resources.create(dpi)) return false;
        uiFont = resources.regular();
        boldFont = resources.semibold();
        smallFont = resources.small();
        monoFont = resources.mono();
        panelBrush = resources.panelBrush();
        cardBrush = resources.cardBrush();

        if (sourceEdit) SendMessageW(sourceEdit, WM_SETFONT, reinterpret_cast<WPARAM>(monoFont), TRUE);
        if (presetCombo) SendMessageW(presetCombo, WM_SETFONT, reinterpret_cast<WPARAM>(uiFont), TRUE);
        if (bailoutEdit) SendMessageW(bailoutEdit, WM_SETFONT, reinterpret_cast<WPARAM>(monoFont), TRUE);
        if (realEdit) SendMessageW(realEdit, WM_SETFONT, reinterpret_cast<WPARAM>(monoFont), TRUE);
        if (imaginaryEdit) SendMessageW(imaginaryEdit, WM_SETFONT, reinterpret_cast<WPARAM>(monoFont), TRUE);
        if (presetCombo) {
            SendMessageW(presetCombo, CB_SETITEMHEIGHT, static_cast<WPARAM>(-1), scale(28));
            SendMessageW(presetCombo, CB_SETITEMHEIGHT, 0, scale(24));
            SendMessageW(presetCombo, CB_SETDROPPEDWIDTH, scale(530), 0);
        }
        return true;
    }

    HWND createControl(const wchar_t* className, const wchar_t* text,
                       DWORD style, DWORD exStyle, int id) const {
        HINSTANCE instance = reinterpret_cast<HINSTANCE>(
            GetWindowLongPtrW(hwnd, GWLP_HINSTANCE));
        return CreateWindowExW(
            exStyle, className, text, WS_CHILD | WS_VISIBLE | style,
            0, 0, 10, 10, hwnd, reinterpret_cast<HMENU>(static_cast<INT_PTR>(id)),
            instance, nullptr);
    }

    bool createControls() {
        sourceEdit = createControl(
            L"EDIT", L"",
            WS_TABSTOP | WS_BORDER | ES_LEFT | ES_AUTOHSCROLL,
            0, ID_SOURCE);
        presetCombo = createControl(
            L"COMBOBOX", L"",
            WS_TABSTOP | WS_BORDER | WS_VSCROLL | CBS_DROPDOWNLIST |
                CBS_OWNERDRAWFIXED | CBS_HASSTRINGS,
            0, ID_PRESET);
        bailoutEdit = createControl(
            L"EDIT", L"",
            WS_TABSTOP | WS_BORDER | ES_RIGHT | ES_AUTOHSCROLL,
            0, ID_BAILOUT);
        realEdit = createControl(
            L"EDIT", L"",
            WS_TABSTOP | WS_BORDER | ES_RIGHT | ES_AUTOHSCROLL,
            0, ID_REAL);
        imaginaryEdit = createControl(
            L"EDIT", L"",
            WS_TABSTOP | WS_BORDER | ES_RIGHT | ES_AUTOHSCROLL,
            0, ID_IMAGINARY);
        if (!sourceEdit || !presetCombo || !bailoutEdit || !realEdit || !imaginaryEdit)
            return false;

        SendMessageW(sourceEdit, EM_SETLIMITTEXT,
                     formula::ExpressionProgram::MAX_SOURCE, 0);
        SendMessageW(bailoutEdit, EM_SETLIMITTEXT, 79, 0);
        SendMessageW(realEdit, EM_SETLIMITTEXT, 79, 0);
        SendMessageW(imaginaryEdit, EM_SETLIMITTEXT, 79, 0);

        for (const Preset& preset : presets())
            SendMessageW(presetCombo, CB_ADDSTRING, 0,
                         reinterpret_cast<LPARAM>(preset.label));

        SetWindowLongPtrW(sourceEdit, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(this));
        SetLastError(ERROR_SUCCESS);
        sourceProc = reinterpret_cast<WNDPROC>(
            SetWindowLongPtrW(sourceEdit, GWLP_WNDPROC,
                              reinterpret_cast<LONG_PTR>(&Impl::sourceWindowProc)));
        if (!sourceProc && GetLastError() != ERROR_SUCCESS) return false;

        return recreateFonts();
    }

    bool onCreate() {
        if (!createControls()) return false;
        syncAllControls();
        updateScrollInfo();
        return true;
    }

    void destroy() {
        if (hwnd && IsWindow(hwnd)) DestroyWindow(hwnd);
        hwnd = nullptr;
        sourceEdit = presetCombo = bailoutEdit = realEdit = imaginaryEdit = nullptr;
        sourceProc = nullptr;
    }

    static bool registerWindowClass(HINSTANCE instance) {
        WNDCLASSEXW existing{};
        existing.cbSize = sizeof(existing);
        if (GetClassInfoExW(instance, PANEL_CLASS, &existing)) return true;

        WNDCLASSEXW windowClass{};
        windowClass.cbSize = sizeof(windowClass);
        windowClass.style = CS_DBLCLKS;
        windowClass.lpfnWndProc = &Impl::windowProc;
        windowClass.hInstance = instance;
        windowClass.hCursor = LoadCursorW(nullptr, IDC_ARROW);
        windowClass.hbrBackground = nullptr;
        windowClass.lpszClassName = PANEL_CLASS;
        if (RegisterClassExW(&windowClass)) return true;
        return GetLastError() == ERROR_CLASS_ALREADY_EXISTS;
    }

    bool create(HWND newOwner, int newDpi, FormulaEditorCallbacks newCallbacks) {
        destroy();
        deleteGdiObjects();
        owner = newOwner;
        callbacks = std::move(newCallbacks);
        dpi = newDpi > 0 ? newDpi : 96;
        if (!owner || !IsWindow(owner)) return false;

        HINSTANCE instance = reinterpret_cast<HINSTANCE>(
            GetWindowLongPtrW(owner, GWLP_HINSTANCE));
        if (!instance) instance = GetModuleHandleW(nullptr);
        if (!registerWindowClass(instance)) return false;

        hwnd = CreateWindowExW(
            WS_EX_CONTROLPARENT, PANEL_CLASS, L"",
            WS_CHILD | WS_CLIPCHILDREN | WS_CLIPSIBLINGS | WS_VSCROLL,
            0, 0, scale(FormulaEditorPanel::DESIGN_WIDTH), scale(600),
            owner, nullptr, instance, this);
        return hwnd != nullptr;
    }

    void setDpi(int newDpi) {
        newDpi = newDpi > 0 ? newDpi : 96;
        if (dpi == newDpi && uiFont) return;
        dpi = newDpi;
        if (!hwnd) return;
        if (!recreateFonts()) return;
        updateScrollInfo();
        layoutControls();
        InvalidateRect(hwnd, nullptr, TRUE);
    }

    void placeControl(HWND control, RECT rect, bool combo = false) {
        if (!control) return;
        int viewportBottom = scale(actionTopDip());
        RECT pixels = toPixelRect(rect, true);
        bool visible = pixels.top >= 0 && pixels.bottom <= viewportBottom &&
                       pixels.right > pixels.left;
        int height = combo ? scale(270) : pixels.bottom - pixels.top;
        SetWindowPos(control, nullptr, pixels.left, pixels.top,
                     std::max(1, static_cast<int>(pixels.right - pixels.left)),
                     std::max(1, height),
                     SWP_NOZORDER | SWP_NOACTIVATE |
                         (visible ? SWP_SHOWWINDOW : SWP_HIDEWINDOW));
    }

    void layoutControls() {
        if (!hwnd) return;
        placeControl(sourceEdit, sourceRectDip());
        placeControl(presetCombo, presetRectDip(), true);
        placeControl(bailoutEdit, bailoutRectDip());
        placeControl(realEdit, realRectDip());
        placeControl(imaginaryEdit, imaginaryRectDip());
    }

    void updateScrollInfo() {
        if (!hwnd) return;
        int page = std::max(1, actionTopDip());
        scroll.configure(CONTENT_HEIGHT, page);
        scroll.apply(hwnd);
        layoutControls();
    }

    void setScrollPosition(int position) {
        if (!scroll.setPosition(position)) return;
        scroll.apply(hwnd);
        layoutControls();
        InvalidateRect(hwnd, nullptr, FALSE);
    }

    void ensureVisible(int top, int height) {
        int page = std::max(1, actionTopDip());
        int position = scroll.position();
        if (top < position) position = top;
        else if (top + height > position + page) position = top + height - page;
        setScrollPosition(position);
    }

    void focusNumeric(HWND edit, int top) {
        ensureVisible(top, 30);
        SetFocus(edit);
        SendMessageW(edit, EM_SETSEL, 0, -1);
    }

    void focusFormula(size_t position = static_cast<size_t>(-1)) {
        ensureVisible(sourceRectDip().top, sourceRectDip().bottom - sourceRectDip().top);
        SetFocus(sourceEdit);
        if (position == static_cast<size_t>(-1)) {
            SendMessageW(sourceEdit, EM_SETSEL, 0, -1);
        } else {
            int length = GetWindowTextLengthW(sourceEdit);
            int caret = std::clamp(static_cast<int>(position), 0, length);
            SendMessageW(sourceEdit, EM_SETSEL, caret, caret);
        }
    }

    int presetIndexForSource(const std::string& source) const {
        const auto& values = presets();
        for (size_t i = 1; i < values.size(); ++i) {
            if (source == values[i].source) return static_cast<int>(i);
        }
        return 0;
    }

    void syncFormulaControls() {
        std::wstring source;
        if (!widenUtf8(working.source, source)) {
            source.clear();
            setStatus(L"Formula source is not valid UTF-8.", true);
        }
        syncing = true;
        SetWindowTextW(sourceEdit, source.c_str());
        SendMessageW(presetCombo, CB_SETCURSEL,
                     presetIndexForSource(working.source), 0);
        syncing = false;
    }

    void ensurePointVisible() {
        const std::complex<double>* value = selectedValue(working);
        if (!value || !finiteComplex(*value)) return;
        double magnitude = std::max(std::abs(value->real()), std::abs(value->imag()));
        while (magnitude > pickerRange * 0.92 && pickerRange < 1.0e12)
            pickerRange *= 2.0;
    }

    void syncInspectorControls() {
        std::complex<double> value{};
        const std::complex<double>* selectedPointer = selectedValue(working);
        if (selectedPointer) value = *selectedPointer;

        syncing = true;
        setNumber(realEdit, value.real());
        setNumber(imaginaryEdit, value.imag());
        syncing = false;

        BOOL enabled = selectedEditable(working) ? TRUE : FALSE;
        EnableWindow(realEdit, enabled);
        EnableWindow(imaginaryEdit, enabled);
        ensurePointVisible();
        InvalidateRect(hwnd, nullptr, FALSE);
    }

    void syncAllControls() {
        if (!sourceEdit) return;
        syncFormulaControls();
        syncing = true;
        setNumber(bailoutEdit, working.bailout);
        syncing = false;
        syncInspectorControls();
        layoutControls();
    }

    void show(const FormulaDialogConfig& config) {
        working = config;
        applied = config;
        selected = config.pixelParameter == FormulaParameter::InitialZ
            ? InspectorValue::C : InspectorValue::Z0;
        pickerRange = 2.0;
        scroll.setPosition(0);
        status.clear();
        statusError = false;
        syncAllControls();
        updateScrollInfo();
        ShowWindow(hwnd, SW_SHOWNOACTIVATE);
        InvalidateRect(hwnd, nullptr, TRUE);
    }

    void setConfig(const FormulaDialogConfig& config) {
        working = config;
        applied = config;
        status.clear();
        statusError = false;
        syncAllControls();
    }

    void hide() {
        if (!hwnd) return;
        if (GetCapture() == hwnd) ReleaseCapture();
        draggingPoint = false;
        pressedHit = HIT_NONE;
        ShowWindow(hwnd, SW_HIDE);
    }

    bool clipboardFailure(const wchar_t* operation) {
        std::wstring message = operation;
        message += L" failed because the clipboard is unavailable.";
        setStatus(std::move(message), true);
        return false;
    }

    bool writeClipboard(const std::wstring& text) {
        if (!OpenClipboard(hwnd)) return clipboardFailure(L"Copy");
        if (!EmptyClipboard()) {
            CloseClipboard();
            return clipboardFailure(L"Copy");
        }

        SIZE_T bytes = (text.size() + 1) * sizeof(wchar_t);
        HGLOBAL memory = GlobalAlloc(GMEM_MOVEABLE, bytes);
        if (!memory) {
            CloseClipboard();
            return clipboardFailure(L"Copy");
        }
        void* data = GlobalLock(memory);
        if (!data) {
            GlobalFree(memory);
            CloseClipboard();
            return clipboardFailure(L"Copy");
        }
        std::memcpy(data, text.c_str(), bytes);
        GlobalUnlock(memory);
        if (!SetClipboardData(CF_UNICODETEXT, memory)) {
            GlobalFree(memory);
            CloseClipboard();
            return clipboardFailure(L"Copy");
        }
        CloseClipboard();
        return true;
    }

    bool readClipboard(std::wstring& text) {
        text.clear();
        if (!OpenClipboard(hwnd)) return clipboardFailure(L"Paste");

        bool ok = false;
        if (IsClipboardFormatAvailable(CF_UNICODETEXT)) {
            HGLOBAL memory = GetClipboardData(CF_UNICODETEXT);
            if (memory) {
                const wchar_t* data = static_cast<const wchar_t*>(GlobalLock(memory));
                SIZE_T size = GlobalSize(memory) / sizeof(wchar_t);
                if (data && size > 0) {
                    size_t length = 0;
                    while (length < size && data[length]) ++length;
                    text.assign(data, length);
                    ok = true;
                }
                if (data) GlobalUnlock(memory);
            }
        } else if (IsClipboardFormatAvailable(CF_TEXT)) {
            HGLOBAL memory = GetClipboardData(CF_TEXT);
            if (memory) {
                const char* data = static_cast<const char*>(GlobalLock(memory));
                SIZE_T size = GlobalSize(memory);
                if (data && size > 0) {
                    size_t length = 0;
                    while (length < size && data[length]) ++length;
                    int count = MultiByteToWideChar(CP_ACP, 0, data,
                                                    static_cast<int>(length),
                                                    nullptr, 0);
                    if (length == 0) {
                        ok = true;
                    } else if (count > 0) {
                        text.resize(static_cast<size_t>(count));
                        if (MultiByteToWideChar(CP_ACP, 0, data,
                                                static_cast<int>(length),
                                                text.data(), count) == count) {
                            ok = true;
                        }
                    }
                }
                if (data) GlobalUnlock(memory);
            }
        }
        CloseClipboard();
        if (!ok) {
            setStatus(L"Paste failed: the clipboard does not contain readable text.", true);
            return false;
        }
        return true;
    }

    bool copyFormulaSelection(bool cut) {
        DWORD start = 0, end = 0;
        SendMessageW(sourceEdit, EM_GETSEL,
                     reinterpret_cast<WPARAM>(&start),
                     reinterpret_cast<LPARAM>(&end));
        if (start == end) return true;
        std::wstring text = getWindowText(sourceEdit);
        size_t first = std::min<size_t>(start, text.size());
        size_t last = std::min<size_t>(end, text.size());
        if (!writeClipboard(text.substr(first, last - first))) return false;
        if (cut) SendMessageW(sourceEdit, EM_REPLACESEL, TRUE,
                              reinterpret_cast<LPARAM>(L""));
        setStatus(cut ? L"Formula selection cut." : L"Formula selection copied.", false);
        return true;
    }

    bool pasteFormulaSelection() {
        std::wstring clipboard;
        if (!readClipboard(clipboard)) return false;
        DWORD start = 0, end = 0;
        SendMessageW(sourceEdit, EM_GETSEL,
                     reinterpret_cast<WPARAM>(&start),
                     reinterpret_cast<LPARAM>(&end));
        size_t currentLength = static_cast<size_t>(GetWindowTextLengthW(sourceEdit));
        size_t replaced = static_cast<size_t>(end - start);
        if (currentLength - replaced + clipboard.size() >
            formula::ExpressionProgram::MAX_SOURCE) {
            setStatus(L"Paste failed: the formula would exceed 4096 characters.", true);
            return false;
        }
        SendMessageW(sourceEdit, EM_REPLACESEL, TRUE,
                     reinterpret_cast<LPARAM>(clipboard.c_str()));
        setStatus(L"Text pasted into the formula.", false);
        return true;
    }

    void copyCompleteFormula() {
        if (writeClipboard(getWindowText(sourceEdit)))
            setStatus(L"Complete formula copied to the clipboard.", false);
    }

    void pasteCompleteFormula() {
        std::wstring clipboard;
        if (!readClipboard(clipboard)) return;
        if (clipboard.size() > formula::ExpressionProgram::MAX_SOURCE) {
            setStatus(L"Paste failed: the formula exceeds 4096 characters.", true);
            return;
        }
        SetWindowTextW(sourceEdit, clipboard.c_str());
        SendMessageW(sourceEdit, EM_SETSEL,
                     static_cast<WPARAM>(clipboard.size()),
                     static_cast<LPARAM>(clipboard.size()));
        setStatus(L"Complete formula replaced from the clipboard.", false);
        focusFormula(clipboard.size());
    }

    bool commitInspector(bool reportErrors) {
        if (!selectedEditable(working)) return true;
        double real = 0.0;
        double imaginary = 0.0;
        if (!parseFiniteNumber(realEdit, real)) {
            if (reportErrors) {
                setStatus(L"Real value must be a finite number.", true);
                focusNumeric(realEdit, realRectDip().top);
            }
            return false;
        }
        if (!parseFiniteNumber(imaginaryEdit, imaginary)) {
            if (reportErrors) {
                setStatus(L"Imaginary value must be a finite number.", true);
                focusNumeric(imaginaryEdit, imaginaryRectDip().top);
            }
            return false;
        }
        std::complex<double>* value = selectedValue(working);
        if (value) *value = { real, imaginary };
        return true;
    }

    void liveSyncInspector() {
        if (syncing || !selectedEditable(working)) return;
        double real = 0.0;
        double imaginary = 0.0;
        if (!parseFiniteNumber(realEdit, real) ||
            !parseFiniteNumber(imaginaryEdit, imaginary)) {
            return;
        }
        std::complex<double>* value = selectedValue(working);
        if (value) {
            *value = { real, imaginary };
            InvalidateRect(hwnd, nullptr, FALSE);
        }
    }

    bool selectInspector(InspectorValue value) {
        if (selected == value) return true;
        if (!commitInspector(true)) return false;
        selected = value;
        syncInspectorControls();
        return true;
    }

    bool insertFormulaText(const std::wstring& token, bool function) {
        DWORD start = 0, end = 0;
        SendMessageW(sourceEdit, EM_GETSEL,
                     reinterpret_cast<WPARAM>(&start),
                     reinterpret_cast<LPARAM>(&end));
        std::wstring insertion = token;
        if (function) insertion += L"()";
        size_t currentLength = static_cast<size_t>(GetWindowTextLengthW(sourceEdit));
        if (currentLength - static_cast<size_t>(end - start) + insertion.size() >
            formula::ExpressionProgram::MAX_SOURCE) {
            setStatus(L"Insertion failed: the formula would exceed 4096 characters.", true);
            return false;
        }
        SendMessageW(sourceEdit, EM_REPLACESEL, TRUE,
                     reinterpret_cast<LPARAM>(insertion.c_str()));
        DWORD caret = start + static_cast<DWORD>(insertion.size());
        if (function) caret = start + static_cast<DWORD>(token.size()) + 1;
        SendMessageW(sourceEdit, EM_SETSEL, caret, caret);
        focusFormula(caret);
        return true;
    }

    void selectIdentifierAtCaret() {
        std::wstring text = getWindowText(sourceEdit);
        if (text.empty()) return;
        DWORD selectionStart = 0, selectionEnd = 0;
        SendMessageW(sourceEdit, EM_GETSEL,
                     reinterpret_cast<WPARAM>(&selectionStart),
                     reinterpret_cast<LPARAM>(&selectionEnd));
        size_t position = std::min<size_t>(selectionEnd, text.size());
        auto identifierCharacter = [](wchar_t ch) {
            return std::iswalnum(ch) || ch == L'_';
        };
        size_t probe = position;
        if (probe >= text.size() || !identifierCharacter(text[probe])) {
            if (probe == 0 || !identifierCharacter(text[probe - 1])) return;
            --probe;
        }
        size_t first = probe;
        size_t last = probe + 1;
        while (first > 0 && identifierCharacter(text[first - 1])) --first;
        while (last < text.size() && identifierCharacter(text[last])) ++last;
        std::wstring identifier = text.substr(first, last - first);
        std::transform(identifier.begin(), identifier.end(), identifier.begin(),
                       [](wchar_t ch) { return std::towlower(ch); });

        InspectorValue value;
        bool found = true;
        if (identifier == L"c") value = InspectorValue::C;
        else if (identifier == L"z0") value = InspectorValue::Z0;
        else if (identifier == L"z") value = InspectorValue::Z;
        else if (identifier == L"n") value = InspectorValue::N;
        else if (identifier.size() == 2 && identifier[0] == L'p' &&
                 identifier[1] >= L'0' && identifier[1] <= L'7') {
            value = static_cast<InspectorValue>(
                static_cast<int>(InspectorValue::P0) + identifier[1] - L'0');
        } else {
            found = false;
        }
        if (found) selectInspector(value);
    }

    void choosePreset() {
        int index = static_cast<int>(
            SendMessageW(presetCombo, CB_GETCURSEL, 0, 0));
        if (index <= 0 || index >= static_cast<int>(presets().size())) return;
        const Preset& preset = presets()[static_cast<size_t>(index)];
        std::wstring source;
        if (!widenUtf8(preset.source, source)) {
            setStatus(L"Preset text is not valid UTF-8.", true);
            return;
        }
        syncing = true;
        SetWindowTextW(sourceEdit, source.c_str());
        SendMessageW(presetCombo, CB_SETCURSEL, index, 0);
        syncing = false;
        working.source = preset.source;
        SendMessageW(sourceEdit, EM_SETSEL, source.size(), source.size());
        setStatus(L"Preset loaded; press Apply to use it.", false);
        selectIdentifierAtCaret();
    }

    void choosePlane(FormulaParameter plane) {
        if (working.pixelParameter == plane) return;
        if (!commitInspector(true)) return;
        working.pixelParameter = plane;
        syncInspectorControls();
        setStatus(plane == FormulaParameter::InitialZ
            ? L"z0-plane selected; fixed c is editable."
            : L"c-plane selected; fixed z0 is editable.", false);
    }

    void changeRange(double factor) {
        pickerRange = std::clamp(pickerRange * factor, 1.0e-12, 1.0e12);
        InvalidateRect(hwnd, nullptr, FALSE);
    }

    void resetRange() {
        pickerRange = 2.0;
        InvalidateRect(hwnd, nullptr, FALSE);
    }

    void updatePointFromMouse(int xDip, int contentYDip) {
        if (!selectedEditable(working)) return;
        RECT plane = planeRectDip();
        double width = static_cast<double>(
            std::max(1L, plane.right - plane.left - 1));
        double height = static_cast<double>(
            std::max(1L, plane.bottom - plane.top - 1));
        double normalizedX = std::clamp(
            (xDip - plane.left) / width, 0.0, 1.0);
        double normalizedY = std::clamp(
            (contentYDip - plane.top) / height, 0.0, 1.0);
        double real = (normalizedX * 2.0 - 1.0) * pickerRange;
        double imaginary = (1.0 - normalizedY * 2.0) * pickerRange;
        std::complex<double>* value = selectedValue(working);
        if (!value) return;
        *value = { real, imaginary };
        syncing = true;
        setNumber(realEdit, real);
        setNumber(imaginaryEdit, imaginary);
        syncing = false;
        setStatus(std::wstring(inspectorName(selected)) +
                  L" updated; press Apply to use it.", false);
        InvalidateRect(hwnd, nullptr, FALSE);
    }

    bool validateStoredComplexes(FormulaDialogConfig& candidate) {
        struct Entry {
            std::complex<double>* value;
            InspectorValue inspector;
            const wchar_t* name;
        };
        std::array<Entry, 9> entries{};
        entries[0] = candidate.pixelParameter == FormulaParameter::InitialZ
            ? Entry{ &candidate.fixedC, InspectorValue::C, L"c" }
            : Entry{ &candidate.fixedZ0, InspectorValue::Z0, L"z0" };
        for (int i = 0; i < 8; ++i) {
            entries[static_cast<size_t>(i + 1)] = {
                &candidate.parameters[static_cast<size_t>(i)],
                static_cast<InspectorValue>(
                    static_cast<int>(InspectorValue::P0) + i),
                inspectorName(static_cast<InspectorValue>(
                    static_cast<int>(InspectorValue::P0) + i))
            };
        }
        for (const Entry& entry : entries) {
            if (finiteComplex(*entry.value)) continue;
            selected = entry.inspector;
            working = candidate;
            syncInspectorControls();
            std::wstring message = entry.name;
            message += L" must have finite real and imaginary values.";
            setStatus(std::move(message), true);
            if (!std::isfinite(entry.value->real()))
                focusNumeric(realEdit, realRectDip().top);
            else
                focusNumeric(imaginaryEdit, imaginaryRectDip().top);
            return false;
        }
        return true;
    }

    bool readCandidate(FormulaDialogConfig& candidate) {
        candidate = working;

        std::wstring sourceText = getWindowText(sourceEdit);
        if (!narrowUtf8(sourceText, candidate.source)) {
            setStatus(L"Formula source cannot be converted to UTF-8.", true);
            focusFormula();
            return false;
        }

        double bailout = 0.0;
        if (!parseFiniteNumber(bailoutEdit, bailout)) {
            setStatus(L"Bailout must be a finite number.", true);
            focusNumeric(bailoutEdit, bailoutRectDip().top);
            return false;
        }
        if (!(bailout > 0.0)) {
            setStatus(L"Bailout must be positive.", true);
            focusNumeric(bailoutEdit, bailoutRectDip().top);
            return false;
        }
        candidate.bailout = bailout;

        double real = 0.0;
        double imaginary = 0.0;
        if (!parseFiniteNumber(realEdit, real)) {
            setStatus(L"Real value must be a finite number.", true);
            focusNumeric(realEdit, realRectDip().top);
            return false;
        }
        if (!parseFiniteNumber(imaginaryEdit, imaginary)) {
            setStatus(L"Imaginary value must be a finite number.", true);
            focusNumeric(imaginaryEdit, imaginaryRectDip().top);
            return false;
        }
        if (selectedEditable(candidate)) {
            std::complex<double>* value = selectedValue(candidate);
            if (value) *value = { real, imaginary };
        }
        if (!validateStoredComplexes(candidate)) return false;

        formula::ExpressionProgram program;
        formula::ExpressionError error;
        if (!program.compile(candidate.source, &error)) {
            std::wstring detail;
            if (!widenUtf8(error.message, detail)) detail = L"Invalid formula.";
            wchar_t prefix[80];
            swprintf_s(prefix, L"Formula error at character %zu: ", error.position + 1);
            setStatus(std::wstring(prefix) + detail, true);
            focusFormula(error.position);
            return false;
        }
        return true;
    }

    void applyChanges() {
        FormulaDialogConfig candidate;
        if (!readCandidate(candidate)) return;
        if (!callbacks.apply) {
            setStatus(L"Apply is unavailable because no callback was provided.", true);
            return;
        }
        if (!callbacks.apply(candidate)) {
            setStatus(L"The formula was not accepted; staged edits were kept.", true);
            return;
        }
        working = candidate;
        applied = candidate;
        syncAllControls();
        setStatus(L"Formula applied.", false);
    }

    void revertChanges() {
        working = applied;
        syncAllControls();
        setStatus(L"Reverted to the last successfully applied configuration.", false);
    }

    void useMandelbrot() {
        if (callbacks.useMandelbrot) callbacks.useMandelbrot();
    }

    void closePanel() {
        hide();
        if (callbacks.close) callbacks.close();
    }

    void invokeHit(int hit) {
        if (hit == HIT_COPY) {
            copyCompleteFormula();
            return;
        }
        if (hit == HIT_PASTE) {
            pasteCompleteFormula();
            return;
        }
        if (hit == HIT_C_PLANE) {
            choosePlane(FormulaParameter::C);
            return;
        }
        if (hit == HIT_Z0_PLANE) {
            choosePlane(FormulaParameter::InitialZ);
            return;
        }
        if (hit == HIT_RANGE_MINUS) {
            changeRange(0.5);
            return;
        }
        if (hit == HIT_RANGE_PLUS) {
            changeRange(2.0);
            return;
        }
        if (hit == HIT_RANGE_RESET) {
            resetRange();
            return;
        }
        if (hit == HIT_REVERT) {
            revertChanges();
            return;
        }
        if (hit == HIT_MANDELBROT) {
            useMandelbrot();
            return;
        }
        if (hit == HIT_APPLY) {
            applyChanges();
            return;
        }
        if (hit == HIT_CLOSE) {
            closePanel();
            return;
        }
        if (hit >= HIT_VARIABLE_BASE &&
            hit < HIT_VARIABLE_BASE + static_cast<int>(variableLabels().size())) {
            size_t index = static_cast<size_t>(hit - HIT_VARIABLE_BASE);
            if (insertFormulaText(variableLabels()[index], false))
                selectInspector(variableFromButton(index));
            return;
        }
        if (hit >= HIT_FUNCTION_BASE &&
            hit < HIT_FUNCTION_BASE + static_cast<int>(functionLabels().size())) {
            size_t index = static_cast<size_t>(hit - HIT_FUNCTION_BASE);
            const wchar_t* token = functionLabels()[index];
            insertFormulaText(token, index != 1);
        }
    }

    void drawButton(HDC dc, const ButtonSpec& button, bool content) {
        RECT rect = toPixelRect(button.rect, content);
        bool active = isButtonActive(button.hit);
        bool hovered = hoverHit == button.hit;
        bool pressed = pressedHit == button.hit;
        ui::ButtonStyle style = active
            ? ui::ButtonStyle::Accent : ui::ButtonStyle::Normal;
        if (button.hit == HIT_MANDELBROT)
            style = ui::ButtonStyle::Positive;
        else if (button.hit == HIT_CLOSE)
            style = ui::ButtonStyle::Subtle;
        ui::drawButton(dc, rect, button.label, uiFont, style,
                       hovered, pressed, true, scale(6));
    }

    void drawContentText(HDC dc, RECT rect, const std::wstring& text,
                         COLORREF color, HFONT font, UINT format) {
        drawText(dc, toPixelRect(rect, true), text, color, font, format);
    }

    void drawPlane(HDC dc) {
        RECT logical = planeRectDip();
        RECT rect = toPixelRect(logical, true);
        fillRect(dc, rect, CLR_CARD);

        HBRUSH borderBrush = CreateSolidBrush(CLR_BORDER);
        FrameRect(dc, &rect, borderBrush);
        DeleteObject(borderBrush);

        HPEN gridPen = CreatePen(PS_SOLID, 1, CLR_TRACK);
        HPEN axisPen = CreatePen(PS_SOLID, 1, CLR_ACCENT);
        HGDIOBJ oldPen = SelectObject(dc, gridPen);
        for (int i = 1; i < 8; ++i) {
            int x = rect.left + MulDiv(rect.right - rect.left, i, 8);
            int y = rect.top + MulDiv(rect.bottom - rect.top, i, 8);
            MoveToEx(dc, x, rect.top, nullptr);
            LineTo(dc, x, rect.bottom);
            MoveToEx(dc, rect.left, y, nullptr);
            LineTo(dc, rect.right, y);
        }
        SelectObject(dc, axisPen);
        int centerX = (rect.left + rect.right) / 2;
        int centerY = (rect.top + rect.bottom) / 2;
        MoveToEx(dc, centerX, rect.top, nullptr);
        LineTo(dc, centerX, rect.bottom);
        MoveToEx(dc, rect.left, centerY, nullptr);
        LineTo(dc, rect.right, centerY);
        SelectObject(dc, oldPen);
        DeleteObject(gridPen);
        DeleteObject(axisPen);

        wchar_t range[64];
        swprintf_s(range, L"\u00b1%.6g", pickerRange);
        RECT rangeRect = rect;
        rangeRect.left += scale(6);
        rangeRect.top += scale(4);
        drawText(dc, rangeRect, range, CLR_TEXT_DIM, smallFont,
                 DT_LEFT | DT_TOP | DT_SINGLELINE);

        const std::complex<double>* value = selectedValue(working);
        if (value && finiteComplex(*value)) {
            double normalizedX = value->real() / pickerRange * 0.5 + 0.5;
            double normalizedY = 0.5 - value->imag() / pickerRange * 0.5;
            normalizedX = std::clamp(normalizedX, 0.0, 1.0);
            normalizedY = std::clamp(normalizedY, 0.0, 1.0);
            int x = rect.left + static_cast<int>(
                normalizedX * std::max(0L, rect.right - rect.left - 1));
            int y = rect.top + static_cast<int>(
                normalizedY * std::max(0L, rect.bottom - rect.top - 1));
            int radius = scale(6);
            COLORREF pointColor = selectedEditable(working) ? CLR_ACCENT_HI : CLR_TEXT_DIM;
            HBRUSH pointBrush = CreateSolidBrush(pointColor);
            HPEN pointPen = CreatePen(PS_SOLID, std::max(1, scale(1)), CLR_BG);
            HGDIOBJ oldBrush = SelectObject(dc, pointBrush);
            oldPen = SelectObject(dc, pointPen);
            Ellipse(dc, x - radius, y - radius, x + radius + 1, y + radius + 1);
            SelectObject(dc, oldBrush);
            SelectObject(dc, oldPen);
            DeleteObject(pointBrush);
            DeleteObject(pointPen);
        }

        if (!selectedEditable(working)) {
            RECT overlay{
                rect.left + scale(64), centerY - scale(18),
                rect.right - scale(64), centerY + scale(18)
            };
            fillRound(dc, overlay, CLR_PANEL, CLR_BORDER, scale(6));
            drawText(dc, overlay, L"Read-only", CLR_TEXT_DIM, uiFont,
                     DT_CENTER | DT_VCENTER | DT_SINGLELINE);
        }
    }

    void paintTo(HDC dc, const RECT& client) {
        fillRect(dc, client, CLR_PANEL);
        int actionTopPixels = scale(actionTopDip());

        int saved = SaveDC(dc);
        IntersectClipRect(dc, client.left, client.top, client.right, actionTopPixels);

        drawContentText(dc, { 18, 12, clientWidthDip() - 18, 36 },
                        L"Formula editor", CLR_TEXT, boldFont,
                        DT_LEFT | DT_VCENTER | DT_SINGLELINE);
        drawContentText(dc, { 18, 38, clientWidthDip() - 18, 58 },
                        L"Formula source (z')", CLR_TEXT_DIM, smallFont,
                        DT_LEFT | DT_VCENTER | DT_SINGLELINE);
        drawContentText(dc, { 18, 98, 78, 124 },
                        L"Preset", CLR_TEXT_DIM, uiFont,
                        DT_LEFT | DT_VCENTER | DT_SINGLELINE);
        drawContentText(dc, { 18, 134, clientWidthDip() - 18, 154 },
                        L"Variables", CLR_TEXT_DIM, smallFont,
                        DT_LEFT | DT_VCENTER | DT_SINGLELINE);
        drawContentText(dc, { 18, 192, clientWidthDip() - 18, 212 },
                        L"Functions", CLR_TEXT_DIM, smallFont,
                        DT_LEFT | DT_VCENTER | DT_SINGLELINE);
        drawContentText(dc, { 18, 294, 258, 313 },
                        L"Pixel-bound axis", CLR_TEXT_DIM, smallFont,
                        DT_LEFT | DT_VCENTER | DT_SINGLELINE);

        RECT bailout = bailoutRectDip();
        drawContentText(dc, { bailout.left, 294, bailout.right, 313 },
                        L"Bailout |z|", CLR_TEXT_DIM, smallFont,
                        DT_LEFT | DT_VCENTER | DT_SINGLELINE);
        drawContentText(dc, { 18, 365, clientWidthDip() - 18, 389 },
                        L"Selected variable and complex value", CLR_TEXT, boldFont,
                        DT_LEFT | DT_VCENTER | DT_SINGLELINE);

        for (const ButtonSpec& button : contentButtons())
            drawButton(dc, button, true);

        drawPlane(dc);

        int inspectorLeft = inspectorLeftDip();
        int inspectorRight = std::max(inspectorLeft + 80, clientWidthDip() - 18);
        drawContentText(dc, { inspectorLeft, 395, inspectorRight, 421 },
                        inspectorName(selected), CLR_ACCENT_HI, boldFont,
                        DT_LEFT | DT_VCENTER | DT_SINGLELINE);
        drawContentText(dc, { inspectorLeft, 424, inspectorRight, 486 },
                        inspectorDescription(), CLR_TEXT_DIM, smallFont,
                        DT_LEFT | DT_TOP | DT_WORDBREAK | DT_END_ELLIPSIS);
        drawContentText(dc, { inspectorLeft, 484, inspectorRight, 502 },
                        L"Real", CLR_TEXT_DIM, smallFont,
                        DT_LEFT | DT_VCENTER | DT_SINGLELINE);
        drawContentText(dc, { inspectorLeft, 544, inspectorRight, 562 },
                        L"Imaginary", CLR_TEXT_DIM, smallFont,
                        DT_LEFT | DT_VCENTER | DT_SINGLELINE);
        wchar_t rangeText[80];
        swprintf_s(rangeText, L"Picker range: \u00b1%.6g", pickerRange);
        drawContentText(dc, { inspectorLeft, 606, inspectorRight, 631 },
                        rangeText, CLR_TEXT_DIM, smallFont,
                        DT_LEFT | DT_VCENTER | DT_SINGLELINE | DT_END_ELLIPSIS);

        RestoreDC(dc, saved);

        RECT actionRect{ 0, actionTopPixels, client.right, client.bottom };
        fillRect(dc, actionRect, CLR_BG);
        HPEN separator = CreatePen(PS_SOLID, 1, CLR_BORDER);
        HGDIOBJ oldPen = SelectObject(dc, separator);
        MoveToEx(dc, 0, actionTopPixels, nullptr);
        LineTo(dc, client.right, actionTopPixels);
        SelectObject(dc, oldPen);
        DeleteObject(separator);

        std::wstring statusText = status.empty()
            ? L"Edits are staged until Apply." : status;
        RECT statusRect{
            scale(18), actionTopPixels + scale(4),
            client.right - scale(18), actionTopPixels + scale(27)
        };
        drawText(dc, statusRect, statusText,
                 statusError ? RGB(245, 132, 132) : CLR_TEXT_DIM,
                 smallFont, DT_LEFT | DT_VCENTER | DT_SINGLELINE | DT_END_ELLIPSIS);

        for (const ButtonSpec& button : actionButtons())
            drawButton(dc, button, false);
    }

    void paint() {
        PAINTSTRUCT paintStruct{};
        HDC target = BeginPaint(hwnd, &paintStruct);
        RECT client{};
        GetClientRect(hwnd, &client);
        int width = std::max(1, static_cast<int>(client.right - client.left));
        int height = std::max(1, static_cast<int>(client.bottom - client.top));

        HDC memory = backBuffer.begin(target, width, height);
        if (memory) {
            paintTo(memory, client);
            backBuffer.present(target, client);
        } else {
            paintTo(target, client);
        }
        EndPaint(hwnd, &paintStruct);
    }

    void drawPresetItem(const DRAWITEMSTRUCT* item) {
        if (!item || item->CtlID != ID_PRESET) return;
        int index = item->itemID == static_cast<UINT>(-1)
            ? static_cast<int>(SendMessageW(presetCombo, CB_GETCURSEL, 0, 0))
            : static_cast<int>(item->itemID);
        if (index < 0 || index >= static_cast<int>(presets().size())) index = 0;
        bool selectedItem = (item->itemState & ODS_SELECTED) != 0;
        fillRect(item->hDC, item->rcItem, selectedItem ? CLR_ACCENT : CLR_CARD);
        RECT textRect = item->rcItem;
        textRect.left += scale(7);
        textRect.right -= scale(5);
        drawText(item->hDC, textRect, presets()[static_cast<size_t>(index)].label,
                 selectedItem ? RGB(255, 255, 255) : CLR_TEXT, uiFont,
                 DT_LEFT | DT_VCENTER | DT_SINGLELINE | DT_END_ELLIPSIS);
        if (item->itemState & ODS_FOCUS) {
            RECT focus = item->rcItem;
            InflateRect(&focus, -scale(2), -scale(2));
            DrawFocusRect(item->hDC, &focus);
        }
    }

    void onCommand(WPARAM wp) {
        int id = LOWORD(wp);
        int notification = HIWORD(wp);
        if (id == ID_PRESET && notification == CBN_SELCHANGE) {
            if (!syncing) choosePreset();
            return;
        }
        if (notification != EN_CHANGE || syncing) return;
        if (id == ID_SOURCE) {
            SendMessageW(presetCombo, CB_SETCURSEL, 0, 0);
            std::wstring text = getWindowText(sourceEdit);
            std::string utf8;
            if (narrowUtf8(text, utf8)) working.source = std::move(utf8);
            setStatus(L"Formula edited; press Apply to use it.", false);
        } else if (id == ID_REAL || id == ID_IMAGINARY) {
            liveSyncInspector();
        } else if (id == ID_BAILOUT) {
            double value = 0.0;
            if (parseFiniteNumber(bailoutEdit, value)) working.bailout = value;
        }
    }

    void onVerticalScroll(WPARAM wp) {
        SCROLLINFO info{};
        info.cbSize = sizeof(info);
        info.fMask = SIF_TRACKPOS;
        GetScrollInfo(hwnd, SB_VERT, &info);
        if (!scroll.handleCommand(LOWORD(wp), info.nTrackPos)) return;
        scroll.apply(hwnd);
        layoutControls();
        InvalidateRect(hwnd, nullptr, FALSE);
    }

    void onMouseMove(int xPixels, int yPixels) {
        int x = unscale(xPixels);
        int y = unscale(yPixels);
        if (draggingPoint) {
            updatePointFromMouse(x, y + scroll.position());
            SetCursor(LoadCursorW(nullptr, IDC_CROSS));
            return;
        }

        int hit = hitTest(x, y);
        if (hoverHit != hit) {
            hoverHit = hit;
            InvalidateRect(hwnd, nullptr, FALSE);
        }
        if (!trackingMouse) {
            TRACKMOUSEEVENT event{};
            event.cbSize = sizeof(event);
            event.dwFlags = TME_LEAVE;
            event.hwndTrack = hwnd;
            trackingMouse = TrackMouseEvent(&event) != FALSE;
        }
        if (hit != HIT_NONE)
            SetCursor(LoadCursorW(nullptr, hit == HIT_PICKER ? IDC_CROSS : IDC_HAND));
    }

    void onLeftButtonDown(int xPixels, int yPixels) {
        int x = unscale(xPixels);
        int y = unscale(yPixels);
        int hit = hitTest(x, y);
        if (hit == HIT_PICKER) {
            if (!selectedEditable(working)) {
                setStatus(inspectorDescription(), true);
                return;
            }
            draggingPoint = true;
            SetCapture(hwnd);
            updatePointFromMouse(x, y + scroll.position());
            return;
        }
        if (hit != HIT_NONE) {
            pressedHit = hit;
            SetCapture(hwnd);
            InvalidateRect(hwnd, nullptr, FALSE);
        }
    }

    void onLeftButtonUp(int xPixels, int yPixels) {
        int x = unscale(xPixels);
        int y = unscale(yPixels);
        if (draggingPoint) {
            updatePointFromMouse(x, y + scroll.position());
            draggingPoint = false;
            if (GetCapture() == hwnd) ReleaseCapture();
            return;
        }
        int hit = hitTest(x, y);
        int invoke = hit == pressedHit ? hit : HIT_NONE;
        pressedHit = HIT_NONE;
        if (GetCapture() == hwnd) ReleaseCapture();
        InvalidateRect(hwnd, nullptr, FALSE);
        if (invoke != HIT_NONE) invokeHit(invoke);
    }

    LRESULT handleMessage(UINT message, WPARAM wp, LPARAM lp) {
        switch (message) {
        case WM_CREATE:
            return onCreate() ? 0 : -1;
        case WM_SIZE:
            updateScrollInfo();
            InvalidateRect(hwnd, nullptr, TRUE);
            return 0;
        case WM_PAINT:
            paint();
            return 0;
        case WM_ERASEBKGND:
            return 1;
        case WM_COMMAND:
            onCommand(wp);
            return 0;
        case WM_DRAWITEM:
            if (wp == ID_PRESET) {
                drawPresetItem(reinterpret_cast<const DRAWITEMSTRUCT*>(lp));
                return TRUE;
            }
            break;
        case WM_MEASUREITEM:
            if (wp == ID_PRESET) {
                MEASUREITEMSTRUCT* item = reinterpret_cast<MEASUREITEMSTRUCT*>(lp);
                item->itemHeight = scale(24);
                return TRUE;
            }
            break;
        case WM_CTLCOLOREDIT: {
            HDC dc = reinterpret_cast<HDC>(wp);
            SetBkColor(dc, CLR_CARD);
            SetTextColor(dc, CLR_TEXT);
            return reinterpret_cast<LRESULT>(cardBrush);
        }
        case WM_CTLCOLORSTATIC:
        case WM_CTLCOLORLISTBOX: {
            HDC dc = reinterpret_cast<HDC>(wp);
            SetBkColor(dc, CLR_CARD);
            SetTextColor(dc, IsWindowEnabled(reinterpret_cast<HWND>(lp))
                ? CLR_TEXT : CLR_TEXT_DIM);
            return reinterpret_cast<LRESULT>(cardBrush);
        }
        case WM_VSCROLL:
            onVerticalScroll(wp);
            return 0;
        case WM_MOUSEWHEEL: {
            int steps = GET_WHEEL_DELTA_WPARAM(wp) / WHEEL_DELTA;
            setScrollPosition(scroll.position() - steps * 48);
            return 0;
        }
        case WM_MOUSEMOVE:
            onMouseMove(GET_X_LPARAM(lp), GET_Y_LPARAM(lp));
            return 0;
        case WM_MOUSELEAVE:
            trackingMouse = false;
            hoverHit = HIT_NONE;
            InvalidateRect(hwnd, nullptr, FALSE);
            return 0;
        case WM_LBUTTONDOWN:
            onLeftButtonDown(GET_X_LPARAM(lp), GET_Y_LPARAM(lp));
            return 0;
        case WM_LBUTTONUP:
            onLeftButtonUp(GET_X_LPARAM(lp), GET_Y_LPARAM(lp));
            return 0;
        case WM_CAPTURECHANGED:
            draggingPoint = false;
            pressedHit = HIT_NONE;
            InvalidateRect(hwnd, nullptr, FALSE);
            return 0;
        case WM_FORMULA_CARET_CHANGED:
            selectIdentifierAtCaret();
            return 0;
        case WM_SETFOCUS:
            if (sourceEdit && IsWindowVisible(sourceEdit)) SetFocus(sourceEdit);
            return 0;
        case WM_CLOSE:
            closePanel();
            return 0;
        case WM_NCDESTROY:
        {
            HWND destroyedWindow = hwnd;
            SetWindowLongPtrW(destroyedWindow, GWLP_USERDATA, 0);
            LRESULT result = DefWindowProcW(destroyedWindow, message, wp, lp);
            hwnd = nullptr;
            return result;
        }
        }
        return DefWindowProcW(hwnd, message, wp, lp);
    }

    static LRESULT CALLBACK windowProc(HWND window, UINT message,
                                       WPARAM wp, LPARAM lp) {
        Impl* self = reinterpret_cast<Impl*>(
            GetWindowLongPtrW(window, GWLP_USERDATA));
        if (message == WM_NCCREATE) {
            CREATESTRUCTW* create = reinterpret_cast<CREATESTRUCTW*>(lp);
            self = static_cast<Impl*>(create->lpCreateParams);
            self->hwnd = window;
            SetWindowLongPtrW(window, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(self));
        }
        if (!self) return DefWindowProcW(window, message, wp, lp);
        return self->handleMessage(message, wp, lp);
    }

    static LRESULT CALLBACK sourceWindowProc(HWND edit, UINT message,
                                             WPARAM wp, LPARAM lp) {
        Impl* self = reinterpret_cast<Impl*>(
            GetWindowLongPtrW(edit, GWLP_USERDATA));
        WNDPROC original = self ? self->sourceProc : nullptr;
        if (!original) return DefWindowProcW(edit, message, wp, lp);

        if (message == WM_COPY) {
            self->copyFormulaSelection(false);
            return 0;
        }
        if (message == WM_CUT) {
            self->copyFormulaSelection(true);
            return 0;
        }
        if (message == WM_PASTE) {
            self->pasteFormulaSelection();
            return 0;
        }
        if (message == WM_CHAR &&
            (wp == 1 || wp == 3 || wp == 22 || wp == 24)) {
            return 0;
        }
        if (message == WM_KEYDOWN && (GetKeyState(VK_CONTROL) & 0x8000)) {
            switch (wp) {
            case 'A':
                SendMessageW(edit, EM_SETSEL, 0, -1);
                return 0;
            case 'C':
                self->copyFormulaSelection(false);
                return 0;
            case 'X':
                self->copyFormulaSelection(true);
                return 0;
            case 'V':
                self->pasteFormulaSelection();
                return 0;
            }
        }

        LRESULT result = CallWindowProcW(original, edit, message, wp, lp);
        if (message == WM_LBUTTONUP || message == WM_KEYUP ||
            message == WM_SETFOCUS) {
            if (self->hwnd) PostMessageW(self->hwnd, WM_FORMULA_CARET_CHANGED, 0, 0);
        }
        if (message == WM_NCDESTROY) {
            self->sourceEdit = nullptr;
            self->sourceProc = nullptr;
        }
        return result;
    }
};

FormulaEditorPanel::FormulaEditorPanel()
    : _impl(std::make_unique<Impl>()) {}

FormulaEditorPanel::~FormulaEditorPanel() = default;

bool FormulaEditorPanel::create(HWND owner, int dpi,
                                FormulaEditorCallbacks callbacks) {
    return _impl->create(owner, dpi, std::move(callbacks));
}

void FormulaEditorPanel::show(const FormulaDialogConfig& config) {
    if (_impl->hwnd) _impl->show(config);
}

void FormulaEditorPanel::hide() {
    _impl->hide();
}

bool FormulaEditorPanel::visible() const {
    return _impl->hwnd && IsWindowVisible(_impl->hwnd) != FALSE;
}

void FormulaEditorPanel::move(const RECT& bounds) {
    if (!_impl->hwnd) return;
    SetWindowPos(_impl->hwnd, nullptr, bounds.left, bounds.top,
                 std::max(0L, bounds.right - bounds.left),
                 std::max(0L, bounds.bottom - bounds.top),
                 SWP_NOZORDER | SWP_NOACTIVATE);
}

void FormulaEditorPanel::setDpi(int dpi) {
    _impl->setDpi(dpi);
}

void FormulaEditorPanel::setConfig(const FormulaDialogConfig& config) {
    if (_impl->hwnd) _impl->setConfig(config);
}

HWND FormulaEditorPanel::hwnd() const {
    return _impl->hwnd;
}
