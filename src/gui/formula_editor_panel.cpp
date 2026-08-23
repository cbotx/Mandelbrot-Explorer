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
#include <cstdint>
#include <cstring>
#include <cwchar>
#include <cwctype>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "formula_expression.h"
#include "formula_editor_accessibility.h"
#include "gui_theme.h"
#include "ui_controls.h"
#include "ui_framework.h"

namespace {

constexpr wchar_t PANEL_CLASS[] = L"MandelFormulaEditorPanelV2";
constexpr UINT CARET_TIMER = 0x7f31;

constexpr int ID_FORMULA_SERVICE = 1801;
constexpr int ID_BAILOUT_SERVICE = 1802;
constexpr int ID_REAL_SERVICE = 1803;
constexpr int ID_IMAGINARY_SERVICE = 1804;

constexpr int HEADER_HEIGHT = 62;
constexpr int FOOTER_HEIGHT = 62;
constexpr int BODY_PADDING = 14;
constexpr int CONTENT_HEIGHT = 710;

constexpr int HIT_NONE = 0;
constexpr int HIT_CLOSE = 1;
constexpr int HIT_COPY = 2;
constexpr int HIT_PASTE = 3;
constexpr int HIT_C_PLANE = 4;
constexpr int HIT_Z0_PLANE = 5;
constexpr int HIT_PICKER = 6;
constexpr int HIT_RANGE_OUT = 7;
constexpr int HIT_RANGE_RESET = 8;
constexpr int HIT_RANGE_IN = 9;
constexpr int HIT_MANDELBROT = 10;
constexpr int HIT_REVERT = 11;
constexpr int HIT_APPLY = 12;
constexpr int HIT_FUNCTION_TAB_BASE = 30;
constexpr int HIT_VARIABLE_BASE = 100;
constexpr int HIT_FUNCTION_BASE = 200;

constexpr int FOCUS_FORMULA = 1000;
constexpr int FOCUS_PRESET = 1001;
constexpr int FOCUS_BAILOUT = 1002;
constexpr int FOCUS_REAL = 1003;
constexpr int FOCUS_IMAGINARY = 1004;
constexpr int FOCUS_SCROLLBAR = 1005;

inline const COLORREF DOCK_BG = RGB(23, 26, 33);
inline const COLORREF HEADER_BG = RGB(27, 31, 40);
inline const COLORREF SOFT_BORDER = RGB(45, 50, 64);
inline const COLORREF EDIT_BG = RGB(17, 20, 26);
inline const COLORREF VALUE_CARD = RGB(29, 33, 42);
inline const COLORREF TOKEN_VARIABLE = RGB(142, 184, 255);
inline const COLORREF TOKEN_FUNCTION = RGB(186, 151, 255);
inline const COLORREF TOKEN_NUMBER = RGB(255, 214, 107);
inline const COLORREF TOKEN_OPERATOR = RGB(135, 147, 168);
inline const COLORREF TOKEN_UNKNOWN = RGB(242, 125, 140);
inline const COLORREF DANGER = RGB(242, 125, 140);

enum class InspectorValue { Z,
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
                            P7 };

enum class TokenKind { Identifier,
                       Number,
                       Operator,
                       Punctuation,
                       Space };

struct Token {
    size_t first = 0;
    size_t last = 0;
    TokenKind kind = TokenKind::Punctuation;
    std::wstring normalized;
    bool function = false;
    bool selectable = false;
    InspectorValue inspector = InspectorValue::Z;
};

struct Preset {
    const wchar_t* label;
    const char* source;
    bool setP0 = false;
    std::complex<double> p0{};
};

enum class ButtonKind { Normal,
                        Chip,
                        Function,
                        Segment,
                        Primary,
                        Icon,
                        Tab,
                        Range };

struct ButtonSpec {
    RECT rect{};
    int hit = HIT_NONE;
    std::wstring label;
    ButtonKind kind = ButtonKind::Normal;
    bool enabled = true;
    bool active = false;
    bool content = false;
};

struct FocusEntry {
    int id = 0;
    RECT rect{};
    RECT contentRect{};
    bool content = false;
    bool enabled = true;
};

struct PanelLayout {
    RECT client{};
    RECT header{};
    RECT body{};
    RECT footer{};
    RECT contentClip{};
    RECT presetCard{};
    RECT preset{};
    RECT formulaCard{};
    RECT formulaField{};
    RECT formulaMeta{};
    RECT variableCard{};
    RECT functionCard{};
    RECT planeCard{};
    RECT capabilityCard{};
    RECT valueCard{};
    RECT bailout{};
    RECT picker{};
    RECT realField{};
    RECT imaginaryField{};
    RECT complexPreview{};
    int contentLeftDip = 0;
    int contentRightDip = 0;
    int contentOriginDip = 0;
};

const std::array<Preset, 8>& presets() {
    static const std::array<Preset, 8> values = {{{L"Custom", nullptr, false, {}}, {L"Quadratic", "z*z+c", false, {}}, {L"Cubic", "z*z*z+c", false, {}}, {L"Burning Ship", "sqr(complex(abs(real(z)),abs(imag(z))))+c", false, {}}, {L"Sine", "sin(z)+c", false, {}}, {L"Parameter polynomial", "z*z+c+p0*z", true, {0.25, 0.0}}, {L"Iteration drift", "z*z+c+0.0001*n", false, {}}, {L"Branch power", "exp(p0*log(z))+c", true, {2.0, 0.0}}}};
    return values;
}

const std::array<const wchar_t*, 12>& variableLabels() {
    static const std::array<const wchar_t*, 12> labels = {L"z", L"c", L"z0", L"n", L"p0", L"p1", L"p2", L"p3", L"p4", L"p5", L"p6", L"p7"};
    return labels;
}

const std::array<std::array<const wchar_t*, 16>, 3>& functionSets() {
    static const std::array<std::array<const wchar_t*, 16>, 3> sets = {{{{L"sqr", L"^", L"sqrt", L"abs", L"sin", L"cos", L"tan", L"exp", L"log", L"conj", L"real", L"imag", L"norm", L"arg", L"polar", L"complex"}}, {{L"conj", L"real", L"imag", L"norm", L"arg", L"abs", L"sqrt", L"complex", L"polar", L"exp", L"log", L"log10", L"sinh", L"cosh", L"tanh", L"pow"}}, {{L"sqr", L"^", L"sqrt", L"abs", L"sin", L"cos", L"tan", L"sinh", L"cosh", L"tanh", L"exp", L"log", L"log10", L"conj", L"complex", L"polar"}}}};
    return sets;
}

bool containsPoint(RECT rect, POINT point) {
    return point.x >= rect.left && point.x < rect.right && point.y >= rect.top && point.y < rect.bottom;
}

RECT intersected(RECT a, RECT b) {
    RECT result{};
    if (!IntersectRect(&result, &a, &b)) return {};
    return result;
}

bool nonempty(RECT rect) {
    return rect.right > rect.left && rect.bottom > rect.top;
}

bool widenUtf8(const std::string& value, std::wstring& result) {
    result.clear();
    if (value.empty()) return true;
    int count = MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, value.data(), static_cast<int>(value.size()), nullptr, 0);
    if (count <= 0) return false;
    result.resize(static_cast<size_t>(count));
    return MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, value.data(), static_cast<int>(value.size()), result.data(), count) == count;
}

bool narrowUtf8(const std::wstring& value, std::string& result) {
    result.clear();
    if (value.empty()) return true;
    int count = WideCharToMultiByte(CP_UTF8, WC_ERR_INVALID_CHARS, value.data(), static_cast<int>(value.size()), nullptr, 0, nullptr, nullptr);
    if (count <= 0) return false;
    result.resize(static_cast<size_t>(count));
    return WideCharToMultiByte(CP_UTF8, WC_ERR_INVALID_CHARS, value.data(), static_cast<int>(value.size()), result.data(), count, nullptr, nullptr) == count;
}

size_t wideIndexForUtf8Offset(const std::string& source, size_t offset) {
    offset = std::min(offset, source.size());
    if (offset == 0) return 0;
    int count = MultiByteToWideChar(CP_UTF8, 0, source.data(), static_cast<int>(offset), nullptr, 0);
    return count > 0 ? static_cast<size_t>(count) : offset;
}

bool parseFiniteNumber(const std::wstring& text, double& value) {
    const wchar_t* begin = text.c_str();
    while (*begin && std::iswspace(*begin)) ++begin;
    if (!*begin) return false;
    errno = 0;
    wchar_t* end = nullptr;
    double parsed = std::wcstod(begin, &end);
    while (end && *end && std::iswspace(*end)) ++end;
    if (errno == ERANGE || !end || end == begin || *end || !std::isfinite(parsed)) { return false; }
    value = parsed;
    return true;
}

std::wstring formatNumber(double value, int precision = 12) {
    wchar_t buffer[96]{};
    swprintf_s(buffer, L"%.*g", precision, value);
    return buffer;
}

std::wstring lowerAscii(std::wstring value) {
    for (wchar_t& ch : value) {
        if (ch >= L'A' && ch <= L'Z') ch = ch - L'A' + L'a';
    }
    return value;
}

bool inspectorForName(const std::wstring& name, InspectorValue& value) {
    if (name == L"z")
        value = InspectorValue::Z;
    else if (name == L"c")
        value = InspectorValue::C;
    else if (name == L"z0")
        value = InspectorValue::Z0;
    else if (name == L"n")
        value = InspectorValue::N;
    else if (name.size() == 2 && name[0] == L'p' && name[1] >= L'0' && name[1] <= L'7') {
        value = static_cast<InspectorValue>(static_cast<int>(InspectorValue::P0) + name[1] - L'0');
    } else {
        return false;
    }
    return true;
}

InspectorValue variableFromButton(size_t index) {
    static const std::array<InspectorValue, 12> values = {InspectorValue::Z, InspectorValue::C, InspectorValue::Z0, InspectorValue::N, InspectorValue::P0, InspectorValue::P1, InspectorValue::P2, InspectorValue::P3, InspectorValue::P4, InspectorValue::P5, InspectorValue::P6, InspectorValue::P7};
    return values[index];
}

int parameterIndex(InspectorValue value) {
    int index = static_cast<int>(value) - static_cast<int>(InspectorValue::P0);
    return index >= 0 && index < 8 ? index : -1;
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

bool finiteComplex(const std::complex<double>& value) {
    return std::isfinite(value.real()) && std::isfinite(value.imag());
}

bool sameConfig(const FormulaDialogConfig& a, const FormulaDialogConfig& b) {
    if (a.source != b.source || a.pixelParameter != b.pixelParameter || a.fixedZ0 != b.fixedZ0 || a.fixedC != b.fixedC || a.bailout != b.bailout) { return false; }
    return a.parameters == b.parameters;
}

std::vector<Token> tokenizeFormula(const std::wstring& text) {
    std::vector<Token> tokens;
    size_t i = 0;
    while (i < text.size()) {
        wchar_t ch = text[i];
        Token token;
        token.first = i;
        if (std::iswspace(ch)) {
            token.kind = TokenKind::Space;
            while (i < text.size() && std::iswspace(text[i])) ++i;
        } else if ((ch >= L'A' && ch <= L'Z') || (ch >= L'a' && ch <= L'z') || ch == L'_') {
            token.kind = TokenKind::Identifier;
            ++i;
            while (i < text.size()) {
                wchar_t next = text[i];
                if (!((next >= L'A' && next <= L'Z') || (next >= L'a' && next <= L'z') || (next >= L'0' && next <= L'9') || next == L'_')) { break; }
                ++i;
            }
            token.normalized = lowerAscii(text.substr(token.first, i - token.first));
            size_t probe = i;
            while (probe < text.size() && std::iswspace(text[probe])) ++probe;
            token.function = probe < text.size() && text[probe] == L'(';
            token.selectable = inspectorForName(token.normalized, token.inspector);
        } else if ((ch >= L'0' && ch <= L'9') || (ch == L'.' && i + 1 < text.size() && text[i + 1] >= L'0' && text[i + 1] <= L'9')) {
            token.kind = TokenKind::Number;
            bool exponent = false;
            ++i;
            while (i < text.size()) {
                wchar_t next = text[i];
                if ((next >= L'0' && next <= L'9') || next == L'.') {
                    ++i;
                } else if (!exponent && (next == L'e' || next == L'E')) {
                    exponent = true;
                    ++i;
                    if (i < text.size() && (text[i] == L'+' || text[i] == L'-')) { ++i; }
                } else {
                    break;
                }
            }
            if (i < text.size() && text[i] == L'i') ++i;
        } else {
            token.kind = (ch == L'+' || ch == L'-' || ch == L'*' || ch == L'/' || ch == L'^' || ch == L'=') ? TokenKind::Operator : TokenKind::Punctuation;
            ++i;
        }
        token.last = i;
        tokens.push_back(std::move(token));
    }
    return tokens;
}

bool writeClipboard(HWND owner, const std::wstring& text) {
    if (!OpenClipboard(owner)) return false;
    if (!EmptyClipboard()) {
        CloseClipboard();
        return false;
    }
    SIZE_T bytes = (text.size() + 1u) * sizeof(wchar_t);
    HGLOBAL memory = GlobalAlloc(GMEM_MOVEABLE, bytes);
    if (!memory) {
        CloseClipboard();
        return false;
    }
    void* data = GlobalLock(memory);
    if (!data) {
        GlobalFree(memory);
        CloseClipboard();
        return false;
    }
    std::memcpy(data, text.c_str(), bytes);
    GlobalUnlock(memory);
    if (!SetClipboardData(CF_UNICODETEXT, memory)) {
        GlobalFree(memory);
        CloseClipboard();
        return false;
    }
    CloseClipboard();
    return true;
}

bool readClipboard(HWND owner, std::wstring& text) {
    text.clear();
    if (!OpenClipboard(owner)) return false;
    HGLOBAL memory = static_cast<HGLOBAL>(GetClipboardData(CF_UNICODETEXT));
    bool ok = false;
    if (memory) {
        const wchar_t* data = static_cast<const wchar_t*>(GlobalLock(memory));
        SIZE_T capacity = GlobalSize(memory) / sizeof(wchar_t);
        if (data && capacity > 0) {
            size_t length = 0;
            while (length < capacity && data[length]) ++length;
            text.assign(data, length);
            ok = true;
        }
        if (data) GlobalUnlock(memory);
    }
    CloseClipboard();
    return ok;
}

std::wstring environmentValue(const wchar_t* name) {
    DWORD length = GetEnvironmentVariableW(name, nullptr, 0);
    if (length == 0) return {};
    std::wstring value(static_cast<size_t>(length), L'\0');
    DWORD copied = GetEnvironmentVariableW(name, value.data(), length);
    if (copied == 0 || copied >= length) return {};
    value.resize(static_cast<size_t>(copied));
    return value;
}

bool writeSelfTestResult(const std::wstring& path, const std::vector<std::wstring>& failures) {
    if (path.empty()) return false;
    std::wstring report = failures.empty() ? L"PASS\r\n" : L"FAIL\r\n";
    for (const std::wstring& failure : failures) {
        report += failure;
        report += L"\r\n";
    }
    int length = WideCharToMultiByte(CP_UTF8, 0, report.data(), static_cast<int>(report.size()), nullptr, 0, nullptr, nullptr);
    if (length <= 0) return false;
    std::string utf8(static_cast<size_t>(length), '\0');
    if (WideCharToMultiByte(CP_UTF8, 0, report.data(), static_cast<int>(report.size()), utf8.data(), length, nullptr, nullptr) != length) { return false; }
    HANDLE file = CreateFileW(path.c_str(), GENERIC_WRITE, FILE_SHARE_READ, nullptr, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, nullptr);
    if (file == INVALID_HANDLE_VALUE) return false;
    DWORD written = 0;
    BOOL okay = WriteFile(file, utf8.data(), static_cast<DWORD>(utf8.size()), &written, nullptr);
    CloseHandle(file);
    return okay && written == utf8.size();
}

HFONT createFont(int dpi, int dipHeight, int weight, const wchar_t* face) {
    return CreateFontW(-std::max(1, MulDiv(dipHeight, dpi, 96)), 0, 0, 0, weight, FALSE, FALSE, FALSE, DEFAULT_CHARSET, OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS, CLEARTYPE_QUALITY, DEFAULT_PITCH | FF_DONTCARE, face);
}

void drawSeparator(HDC dc, int x1, int y1, int x2, int y2, COLORREF color) {
    HPEN pen = CreatePen(PS_SOLID, 1, color);
    if (!pen) return;
    HGDIOBJ old = SelectObject(dc, pen);
    MoveToEx(dc, x1, y1, nullptr);
    LineTo(dc, x2, y2);
    SelectObject(dc, old);
    DeleteObject(pen);
}

void drawCircle(HDC dc, int centerX, int centerY, int radius, COLORREF fill, COLORREF border) {
    HBRUSH brush = CreateSolidBrush(fill);
    HPEN pen = CreatePen(PS_SOLID, 1, border);
    if (!brush || !pen) {
        if (brush) DeleteObject(brush);
        if (pen) DeleteObject(pen);
        return;
    }
    HGDIOBJ oldBrush = SelectObject(dc, brush);
    HGDIOBJ oldPen = SelectObject(dc, pen);
    Ellipse(dc, centerX - radius, centerY - radius, centerX + radius + 1, centerY + radius + 1);
    SelectObject(dc, oldBrush);
    SelectObject(dc, oldPen);
    DeleteObject(brush);
    DeleteObject(pen);
}

} // namespace

struct FormulaEditorPanel::Impl {
    struct AccessibilityGuard {
        bool action = false;
        int snapshotDepth = 0;
    };

    HWND hwnd = nullptr;
    HWND owner = nullptr;
    FormulaEditorAccessibilityProvider* accessibility = nullptr;
    std::shared_ptr<AccessibilityGuard> accessibilityGuard = std::make_shared<AccessibilityGuard>();
    FormulaEditorCallbacks callbacks;
    FormulaDialogConfig working;
    FormulaDialogConfig applied;

    int dpi = 96;
    ui::Resources resources;
    ui::BackBuffer backBuffer;
    ui::TextField formulaField;
    ui::TextField bailoutField;
    ui::TextField realField;
    ui::TextField imaginaryField;
    ui::Dropdown presetDropdown;
    ui::Scrollbar scrollbar;
    ui::HitRouter hitRouter;

    HFONT formulaFont = nullptr;
    HFONT chipFont = nullptr;
    HFONT tinyFont = nullptr;
    HFONT sectionFont = nullptr;
    HFONT functionFont = nullptr;
    HFONT controlFont = nullptr;
    HFONT headerFont = nullptr;
    HFONT actionFont = nullptr;

    PanelLayout layout;
    std::vector<ButtonSpec> buttons;
    std::vector<FocusEntry> focusEntries;
    std::vector<Token> tokens;
    formula::ExpressionProgram formulaProgram;
    formula::ExpressionError formulaError;

    InspectorValue selected = InspectorValue::Z0;
    int functionTab = 0;
    int focusedId = FOCUS_FORMULA;
    int hoverHit = HIT_NONE;
    int pressedHit = HIT_NONE;
    int hoveredFormulaToken = -1;
    double pickerRange = 2.0;
    bool formulaValid = false;
    bool syncing = false;
    bool trackingMouse = false;
    bool draggingPicker = false;
    bool selfTestRan = false;
    ui::TextField* draggingText = nullptr;
    std::wstring status;
    bool statusError = false;
    std::shared_ptr<int> lifetime = std::make_shared<int>(0);

    ~Impl() {
        lifetime.reset();
        destroy();
        destroyFormulaEditorAccessibility(accessibility);
        deleteExtraFonts();
    }

    int scale(int value) const { return MulDiv(value, dpi > 0 ? dpi : 96, 96); }

    int unscale(int value) const { return MulDiv(value, 96, dpi > 0 ? dpi : 96); }

    int clientWidthDip() const {
        if (!hwnd) return FormulaEditorPanel::DESIGN_WIDTH;
        RECT rect{};
        GetClientRect(hwnd, &rect);
        return std::max(1, unscale(rect.right - rect.left));
    }

    int clientHeightDip() const {
        if (!hwnd) return 760;
        RECT rect{};
        GetClientRect(hwnd, &rect);
        return std::max(1, unscale(rect.bottom - rect.top));
    }

    RECT pixelRect(int left, int top, int right, int bottom) const { return {scale(left), scale(top), scale(right), scale(bottom)}; }

    RECT contentPixelRect(int left, int top, int right, int bottom) const {
        int origin = HEADER_HEIGHT + BODY_PADDING - scrollbar.position();
        return pixelRect(left, origin + top, right, origin + bottom);
    }

    RECT contentLocalRect(int left, int top, int right, int bottom) const { return {left, top, right, bottom}; }

    void deleteExtraFonts() {
        if (formulaFont) DeleteObject(formulaFont);
        if (chipFont) DeleteObject(chipFont);
        if (tinyFont) DeleteObject(tinyFont);
        if (sectionFont) DeleteObject(sectionFont);
        if (functionFont) DeleteObject(functionFont);
        if (controlFont) DeleteObject(controlFont);
        if (headerFont) DeleteObject(headerFont);
        if (actionFont) DeleteObject(actionFont);
        formulaFont = chipFont = tinyFont = sectionFont = nullptr;
        functionFont = controlFont = headerFont = actionFont = nullptr;
    }

    bool recreateResources() {
        if (!resources.create(dpi)) return false;
        HFONT newFormula = createFont(dpi, 16, FW_MEDIUM, L"Cascadia Mono");
        HFONT newChip = createFont(dpi, 12, FW_MEDIUM, L"Cascadia Mono");
        HFONT newTiny = createFont(dpi, 10, FW_NORMAL, L"Segoe UI");
        HFONT newSection = createFont(dpi, 10, FW_BOLD, L"Segoe UI");
        HFONT newFunction = createFont(dpi, 10, FW_MEDIUM, L"Cascadia Mono");
        HFONT newControl = createFont(dpi, 11, FW_MEDIUM, L"Cascadia Mono");
        HFONT newHeader = createFont(dpi, 16, FW_BOLD, L"Segoe UI");
        HFONT newAction = createFont(dpi, 14, FW_SEMIBOLD, L"Segoe UI");
        if (!newFormula || !newChip || !newTiny || !newSection || !newFunction || !newControl || !newHeader || !newAction) {
            if (newFormula) DeleteObject(newFormula);
            if (newChip) DeleteObject(newChip);
            if (newTiny) DeleteObject(newTiny);
            if (newSection) DeleteObject(newSection);
            if (newFunction) DeleteObject(newFunction);
            if (newControl) DeleteObject(newControl);
            if (newHeader) DeleteObject(newHeader);
            if (newAction) DeleteObject(newAction);
            return false;
        }
        deleteExtraFonts();
        formulaFont = newFormula;
        chipFont = newChip;
        tinyFont = newTiny;
        sectionFont = newSection;
        functionFont = newFunction;
        controlFont = newControl;
        headerFont = newHeader;
        actionFont = newAction;
        return true;
    }

    bool selectedEditable(const FormulaDialogConfig& config) const {
        int parameter = parameterIndex(selected);
        if (parameter >= 0) return true;
        if (selected == InspectorValue::C) return config.pixelParameter == FormulaParameter::InitialZ;
        if (selected == InspectorValue::Z0) return config.pixelParameter != FormulaParameter::InitialZ;
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

    std::wstring inspectorTitle() const {
        int parameter = parameterIndex(selected);
        if (parameter >= 0) return std::wstring(L"Parameter ") + inspectorName(selected);
        if (selected == InspectorValue::N) return L"Iteration index n";
        if (selected == InspectorValue::Z) return L"Orbit state z";
        return std::wstring(L"Variable ") + inspectorName(selected);
    }

    std::wstring inspectorDescription() const {
        if (parameterIndex(selected) >= 0) return L"Drag in the square or enter exact values";
        if (selected == InspectorValue::Z) return L"Computed orbit value; read-only";
        if (selected == InspectorValue::N) return L"Current iteration index; read-only";
        if (!selectedEditable(working)) return L"Bound to each pixel in the selected plane";
        return L"Fixed complex value for the selected plane";
    }

    int usageCount(InspectorValue value) const {
        int count = 0;
        for (const Token& token : tokens) {
            if (token.selectable && token.inspector == value) ++count;
        }
        return count;
    }

    std::wstring usageLabel() const {
        int count = usageCount(selected);
        if (count == 0) return L"not used";
        if (count == 1) return L"used once";
        wchar_t buffer[48]{};
        swprintf_s(buffer, L"used %d times", count);
        return buffer;
    }

    bool dirty() const { return !sameConfig(working, applied); }

    void setStatus(std::wstring text, bool error) {
        status = std::move(text);
        statusError = error;
        if (hwnd) InvalidateRect(hwnd, nullptr, FALSE);
    }

    void notifyAccessibility(DWORD event, LONG key) { notifyFormulaEditorAccessibility(accessibility, event, key); }

    void refreshFormulaAnalysis() {
        tokens = tokenizeFormula(formulaField.text());
        formulaError = {};
        formulaValid = formulaProgram.compile(working.source, &formulaError);
    }

    int presetIndexForSource(const std::string& source) const {
        const auto& values = presets();
        for (size_t i = 1; i < values.size(); ++i) {
            if (source == values[i].source) return static_cast<int>(i);
        }
        return 0;
    }

    InspectorValue initialInspector() const {
        for (int i = 0; i < 8; ++i) {
            InspectorValue value = static_cast<InspectorValue>(static_cast<int>(InspectorValue::P0) + i);
            if (usageCount(value) > 0) return value;
        }
        return working.pixelParameter == FormulaParameter::InitialZ ? InspectorValue::C : InspectorValue::Z0;
    }

    void ensurePointVisible() {
        const std::complex<double>* value = selectedValue(working);
        if (!value || !finiteComplex(*value)) return;
        double magnitude = std::max(std::abs(value->real()), std::abs(value->imag()));
        while (magnitude > pickerRange * 0.92 && pickerRange < 1.0e12) pickerRange *= 2.0;
    }

    void syncInspectorFields() {
        std::complex<double> value{};
        const std::complex<double>* pointer = selectedValue(working);
        if (pointer) value = *pointer;
        syncing = true;
        realField.setText(formatNumber(value.real(), 17));
        imaginaryField.setText(formatNumber(value.imag(), 17));
        bool editable = selectedEditable(working);
        realField.setEnabled(editable);
        imaginaryField.setEnabled(editable);
        syncing = false;
        ensurePointVisible();
        if (hwnd) InvalidateRect(hwnd, nullptr, FALSE);
    }

    void syncAllControls() {
        std::wstring source;
        if (!widenUtf8(working.source, source)) {
            source.clear();
            setStatus(L"Formula source is not valid UTF-8.", true);
        }
        syncing = true;
        formulaField.setText(source);
        formulaField.setSelection(source.size(), source.size());
        bailoutField.setText(formatNumber(working.bailout, 17));
        presetDropdown.setSelectedIndex(presetIndexForSource(working.source));
        syncing = false;
        refreshFormulaAnalysis();
        syncInspectorFields();
        updateLayout();
    }

    void onFormulaEvent(ui::TextField::Event event, const std::wstring& detail) {
        if (event == ui::TextField::Event::ClipboardError) {
            setStatus(detail, true);
            return;
        }
        if (event == ui::TextField::Event::TabForward) {
            focusNext(false);
            return;
        }
        if (event == ui::TextField::Event::TabBackward) {
            focusNext(true);
            return;
        }
        if (event == ui::TextField::Event::FocusGained) {
            focusedId = FOCUS_FORMULA;
            presetDropdown.setFocused(false);
            notifyAccessibility(EVENT_OBJECT_FOCUS, FORMULA_ACC_EXPRESSION);
        }
        if (event == ui::TextField::Event::Changed && !syncing) {
            std::wstring source = formulaField.text();
            std::string utf8;
            if (!narrowUtf8(source, utf8)) {
                formulaValid = false;
                setStatus(L"Formula contains text that is not valid UTF-8.", true);
            } else {
                working.source = std::move(utf8);
                presetDropdown.setSelectedIndex(0);
                refreshFormulaAnalysis();
                if (formulaValid) {
                    setStatus(L"Formula edited; changes are staged.", false);
                } else {
                    std::wstring error;
                    if (!widenUtf8(formulaError.message, error)) error = L"Invalid expression";
                    setStatus(L"Invalid formula: " + error, true);
                }
                selectInspectorAtCaret(false);
                notifyAccessibility(EVENT_OBJECT_VALUECHANGE, FORMULA_ACC_EXPRESSION);
            }
        } else if (event == ui::TextField::Event::CaretMoved) {
            selectInspectorAtCaret(false);
        }
        if (hwnd) InvalidateRect(hwnd, nullptr, FALSE);
    }

    void onNumericEvent(int focus, ui::TextField::Event event, const std::wstring& detail) {
        if (event == ui::TextField::Event::ClipboardError) {
            setStatus(detail, true);
            return;
        }
        if (event == ui::TextField::Event::TabForward) {
            focusNext(false);
            return;
        }
        if (event == ui::TextField::Event::TabBackward) {
            focusNext(true);
            return;
        }
        if (event == ui::TextField::Event::FocusGained) {
            focusedId = focus;
            presetDropdown.setFocused(false);
            LONG key = focus == FOCUS_BAILOUT ? FORMULA_ACC_BAILOUT : (focus == FOCUS_REAL ? FORMULA_ACC_REAL : FORMULA_ACC_IMAGINARY);
            notifyAccessibility(EVENT_OBJECT_FOCUS, key);
        }
        if (event != ui::TextField::Event::Changed || syncing) {
            if (hwnd) InvalidateRect(hwnd, nullptr, FALSE);
            return;
        }

        if (focus == FOCUS_BAILOUT) {
            double value = 0.0;
            if (parseFiniteNumber(bailoutField.text(), value) && value > 0.0) {
                working.bailout = value;
                setStatus(L"Bailout changed; changes are staged.", false);
                notifyAccessibility(EVENT_OBJECT_VALUECHANGE, FORMULA_ACC_BAILOUT);
            }
        } else if (selectedEditable(working)) {
            double real = 0.0;
            double imaginary = 0.0;
            if (parseFiniteNumber(realField.text(), real) && parseFiniteNumber(imaginaryField.text(), imaginary)) {
                std::complex<double>* value = selectedValue(working);
                if (value) {
                    *value = {real, imaginary};
                    ensurePointVisible();
                    setStatus(std::wstring(inspectorName(selected)) + L" changed; changes are staged.", false);
                    notifyAccessibility(EVENT_OBJECT_VALUECHANGE, focus == FOCUS_REAL ? FORMULA_ACC_REAL : FORMULA_ACC_IMAGINARY);
                    notifyAccessibility(EVENT_OBJECT_VALUECHANGE, FORMULA_ACC_PICKER);
                }
            }
        }
        if (hwnd) InvalidateRect(hwnd, nullptr, FALSE);
    }

    bool createControls() {
        if (!formulaField.create(hwnd, ID_FORMULA_SERVICE, formula::ExpressionProgram::MAX_SOURCE, [this](ui::TextField::Event event, const std::wstring& detail) { onFormulaEvent(event, detail); })) { return false; }
        if (!bailoutField.create(hwnd, ID_BAILOUT_SERVICE, 79, [this](ui::TextField::Event event, const std::wstring& detail) { onNumericEvent(FOCUS_BAILOUT, event, detail); })) { return false; }
        if (!realField.create(hwnd, ID_REAL_SERVICE, 79, [this](ui::TextField::Event event, const std::wstring& detail) { onNumericEvent(FOCUS_REAL, event, detail); })) { return false; }
        if (!imaginaryField.create(hwnd, ID_IMAGINARY_SERVICE, 79, [this](ui::TextField::Event event, const std::wstring& detail) { onNumericEvent(FOCUS_IMAGINARY, event, detail); })) { return false; }
        formulaField.setPlaceholder(L"Enter an orbit expression");
        bailoutField.setPlaceholder(L"100");
        realField.setPlaceholder(L"0");
        imaginaryField.setPlaceholder(L"0");

        std::vector<std::wstring> presetItems;
        for (const Preset& preset : presets()) presetItems.emplace_back(preset.label);
        presetDropdown.setItems(std::move(presetItems));
        presetDropdown.setCallback([this](int index) { choosePreset(index); });
        scrollbar.setCallback([this](int) {
            updateLayout();
            if (hwnd) InvalidateRect(hwnd, nullptr, FALSE);
            notifyAccessibility(EVENT_OBJECT_VALUECHANGE, FORMULA_ACC_SCROLLBAR);
        });
        return true;
    }

    void destroyControls() {
        formulaField.destroy();
        bailoutField.destroy();
        realField.destroy();
        imaginaryField.destroy();
        presetDropdown.close();
    }

    bool onCreate() {
        if (!recreateResources()) return false;
        if (!createControls()) return false;
        accessibility = createFormulaEditorAccessibility(hwnd);
        if (!accessibility) return false;
        SetTimer(hwnd, CARET_TIMER, 250, nullptr);
        refreshFormulaAnalysis();
        updateLayout();
        return true;
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
        lifetime = std::make_shared<int>(0);
        owner = newOwner;
        accessibilityGuard = std::make_shared<AccessibilityGuard>();
        callbacks = std::move(newCallbacks);
        dpi = newDpi > 0 ? newDpi : 96;
        if (!owner || !IsWindow(owner)) return false;

        HINSTANCE instance = reinterpret_cast<HINSTANCE>(GetWindowLongPtrW(owner, GWLP_HINSTANCE));
        if (!instance) instance = GetModuleHandleW(nullptr);
        if (!registerWindowClass(instance)) return false;

        hwnd = CreateWindowExW(WS_EX_CONTROLPARENT, PANEL_CLASS, L"", WS_CHILD | WS_CLIPCHILDREN | WS_CLIPSIBLINGS, 0, 0, scale(FormulaEditorPanel::DESIGN_WIDTH), scale(870), owner, nullptr, instance, this);
        return hwnd != nullptr;
    }

    void destroy() {
        if (hwnd && IsWindow(hwnd)) DestroyWindow(hwnd);
        destroyFormulaEditorAccessibility(accessibility);
        hwnd = nullptr;
        owner = nullptr;
        draggingText = nullptr;
        draggingPicker = false;
        pressedHit = HIT_NONE;
    }

    void setDpi(int newDpi) {
        newDpi = newDpi > 0 ? newDpi : 96;
        if (newDpi == dpi && formulaFont) return;
        dpi = newDpi;
        if (!hwnd) return;
        if (!recreateResources()) return;
        updateLayout();
        InvalidateRect(hwnd, nullptr, TRUE);
    }

    void addButton(RECT rect, int hit, const wchar_t* label, ButtonKind kind, bool enabled, bool active, bool content, RECT local = {}) {
        ButtonSpec button{rect, hit, label, kind, enabled, active, content};
        buttons.push_back(std::move(button));
        RECT hitRect = content ? intersected(rect, layout.contentClip) : rect;
        if (nonempty(hitRect)) hitRouter.add(hit, hitRect, enabled);
        focusEntries.push_back({hit, rect, local, content, enabled});
    }

    void addFocus(int id, RECT rect, RECT local, bool content, bool enabled) { focusEntries.push_back({id, rect, local, content, enabled}); }

    void updateLayout() {
        if (!hwnd) return;
        RECT client{};
        GetClientRect(hwnd, &client);
        layout.client = client;

        int widthDip = clientWidthDip();
        int heightDip = clientHeightDip();
        int bodyTop = HEADER_HEIGHT;
        int bodyBottom = std::max(bodyTop, heightDip - FOOTER_HEIGHT);
        int viewport = std::max(1, bodyBottom - bodyTop - BODY_PADDING * 2);
        scrollbar.configure(CONTENT_HEIGHT, viewport);

        int scrollbarDip = scrollbar.visible() ? 10 : 0;
        int contentLeft = 16;
        int contentRight = std::max(contentLeft + 420, widthDip - 16 - scrollbarDip);
        int available = contentRight - contentLeft;
        int leftWidth = std::clamp(available * 46 / 100, 220, 246);
        int gap = 12;
        int rightLeft = contentLeft + leftWidth + gap;
        int rightWidth = std::max(190, contentRight - rightLeft);

        layout.contentLeftDip = contentLeft;
        layout.contentRightDip = contentRight;
        layout.contentOriginDip = HEADER_HEIGHT + BODY_PADDING - scrollbar.position();
        layout.header = pixelRect(0, 0, widthDip, HEADER_HEIGHT);
        layout.body = pixelRect(0, bodyTop, widthDip, bodyBottom);
        layout.footer = pixelRect(0, bodyBottom, widthDip, heightDip);
        layout.contentClip = layout.body;

        layout.formulaCard = contentPixelRect(contentLeft, 21, contentRight, 112);
        layout.formulaField = contentPixelRect(contentLeft + 63, 27, contentRight - 12, 79);
        layout.formulaMeta = contentPixelRect(contentLeft + 1, 80, contentRight - 1, 111);

        layout.variableCard = contentPixelRect(contentLeft, 124, contentLeft + leftWidth, 268);
        layout.functionCard = contentPixelRect(contentLeft, 278, contentLeft + leftWidth, 477);
        layout.planeCard = contentPixelRect(contentLeft, 489, contentLeft + leftWidth, 604);
        layout.capabilityCard = contentPixelRect(contentLeft, 614, contentLeft + leftWidth, 659);
        layout.presetCard = contentPixelRect(contentLeft, 669, contentLeft + leftWidth, 710);
        layout.valueCard = contentPixelRect(rightLeft, 124, contentRight, 718);

        layout.bailout = contentPixelRect(contentLeft + 82, 566, contentLeft + leftWidth - 10, 593);

        int presetInnerLeft = contentLeft + 53;
        int presetInnerRight = contentLeft + leftWidth - 10;
        int smallActionWidth = 39;
        int pasteLeft = presetInnerRight - smallActionWidth;
        int copyLeft = pasteLeft - 5 - smallActionWidth;
        int dropdownLeft = presetInnerLeft;
        int dropdownRight = std::max(dropdownLeft + 70, copyLeft - 6);
        layout.preset = contentPixelRect(dropdownLeft, 675, dropdownRight, 704);

        int pickerSize = std::clamp(rightWidth - 44, 170, 235);
        int pickerLeft = rightLeft + 31;
        int pickerTop = 224;
        layout.picker = contentPixelRect(pickerLeft, pickerTop, pickerLeft + pickerSize, pickerTop + pickerSize);
        int valueInnerLeft = rightLeft + 13;
        int valueInnerRight = contentRight - 13;
        int fieldGap = 7;
        int fieldWidth = std::max(60, (valueInnerRight - valueInnerLeft - fieldGap) / 2);
        layout.realField = contentPixelRect(valueInnerLeft, 528, valueInnerLeft + fieldWidth, 559);
        layout.imaginaryField = contentPixelRect(valueInnerLeft + fieldWidth + fieldGap, 528, valueInnerRight, 559);
        layout.complexPreview = contentPixelRect(valueInnerLeft, 568, valueInnerRight, 595);

        formulaField.setBounds(layout.formulaField);
        bailoutField.setBounds(layout.bailout);
        realField.setBounds(layout.realField);
        imaginaryField.setBounds(layout.imaginaryField);
        presetDropdown.setBounds(layout.preset);
        presetDropdown.setFocused(focusedId == FOCUS_PRESET);
        scrollbar.setBounds(pixelRect(widthDip - 10, bodyTop, widthDip, bodyBottom));

        buttons.clear();
        focusEntries.clear();
        hitRouter.clear();

        addButton(pixelRect(widthDip - 46, 16, widthDip - 16, 46), HIT_CLOSE, L"\u00d7", ButtonKind::Icon, true, false, false);
        addFocus(FOCUS_FORMULA, layout.formulaField, contentLocalRect(contentLeft + 63, 27, contentRight - 12, 79), true, true);
        addFocus(FOCUS_PRESET, layout.preset, contentLocalRect(dropdownLeft, 675, dropdownRight, 704), true, true);
        addButton(contentPixelRect(copyLeft, 675, copyLeft + smallActionWidth, 704), HIT_COPY, L"Copy", ButtonKind::Normal, true, false, true, contentLocalRect(copyLeft, 675, copyLeft + smallActionWidth, 704));
        addButton(contentPixelRect(pasteLeft, 675, pasteLeft + smallActionWidth, 704), HIT_PASTE, L"Paste", ButtonKind::Normal, true, false, true, contentLocalRect(pasteLeft, 675, pasteLeft + smallActionWidth, 704));

        int chipLeft = contentLeft + 11;
        int chipRight = contentLeft + leftWidth - 11;
        int chipGap = 5;
        int chipWidth = (chipRight - chipLeft - chipGap * 3) / 4;
        for (int i = 0; i < 12; ++i) {
            int row = i / 4;
            int column = i % 4;
            int left = chipLeft + column * (chipWidth + chipGap);
            int top = 157 + row * 35;
            InspectorValue value = variableFromButton(static_cast<size_t>(i));
            addButton(contentPixelRect(left, top, left + chipWidth, top + 29), HIT_VARIABLE_BASE + i, variableLabels()[static_cast<size_t>(i)], ButtonKind::Chip, true, selected == value, true, contentLocalRect(left, top, left + chipWidth, top + 29));
        }

        int tabLeft = contentLeft + 11;
        int tabRight = contentLeft + leftWidth - 11;
        int tabGap = 4;
        int tabWidth = (tabRight - tabLeft - tabGap * 2) / 3;
        static const std::array<const wchar_t*, 3> tabLabels = {L"Common", L"Complex", L"All"};
        for (int i = 0; i < 3; ++i) {
            int left = tabLeft + i * (tabWidth + tabGap);
            int right = i == 2 ? tabRight : left + tabWidth;
            addButton(contentPixelRect(left, 311, right, 336), HIT_FUNCTION_TAB_BASE + i, tabLabels[static_cast<size_t>(i)], ButtonKind::Tab, true, functionTab == i, true, contentLocalRect(left, 311, right, 336));
        }

        int functionLeft = contentLeft + 11;
        int functionRight = contentLeft + leftWidth - 11;
        int functionGap = 4;
        int functionWidth = (functionRight - functionLeft - functionGap * 3) / 4;
        const auto& labels = functionSets()[static_cast<size_t>(functionTab)];
        for (int i = 0; i < 16; ++i) {
            int row = i / 4;
            int column = i % 4;
            int left = functionLeft + column * (functionWidth + functionGap);
            int top = 343 + row * 33;
            addButton(contentPixelRect(left, top, left + functionWidth, top + 27), HIT_FUNCTION_BASE + i, labels[static_cast<size_t>(i)], ButtonKind::Function, true, false, true, contentLocalRect(left, top, left + functionWidth, top + 27));
        }

        int segmentLeft = contentLeft + 11;
        int segmentRight = contentLeft + leftWidth - 11;
        int segmentMiddle = (segmentLeft + segmentRight) / 2;
        addButton(contentPixelRect(segmentLeft, 523, segmentMiddle, 554), HIT_C_PLANE, L"c plane", ButtonKind::Segment, true, working.pixelParameter != FormulaParameter::InitialZ, true, contentLocalRect(segmentLeft, 523, segmentMiddle, 554));
        addButton(contentPixelRect(segmentMiddle, 523, segmentRight, 554), HIT_Z0_PLANE, L"z0 plane", ButtonKind::Segment, true, working.pixelParameter == FormulaParameter::InitialZ, true, contentLocalRect(segmentMiddle, 523, segmentRight, 554));
        addFocus(FOCUS_BAILOUT, layout.bailout, contentLocalRect(contentLeft + 82, 566, contentLeft + leftWidth - 10, 593), true, true);

        addButton(layout.picker, HIT_PICKER, L"", ButtonKind::Normal, true, false, true, contentLocalRect(pickerLeft, pickerTop, pickerLeft + pickerSize, pickerTop + pickerSize));

        int rangeRight = valueInnerRight;
        int rangeButtonWidth = 25;
        int rangeGap = 4;
        int rangeLeft = rangeRight - rangeButtonWidth * 3 - rangeGap * 2;
        addButton(contentPixelRect(rangeLeft, 486, rangeLeft + rangeButtonWidth, 509), HIT_RANGE_OUT, L"-", ButtonKind::Range, true, false, true, contentLocalRect(rangeLeft, 486, rangeLeft + rangeButtonWidth, 509));
        addButton(contentPixelRect(rangeLeft + rangeButtonWidth + rangeGap, 486, rangeLeft + rangeButtonWidth * 2 + rangeGap, 509), HIT_RANGE_RESET, L"0", ButtonKind::Range, true, false, true, contentLocalRect(rangeLeft + rangeButtonWidth + rangeGap, 486, rangeLeft + rangeButtonWidth * 2 + rangeGap, 509));
        addButton(contentPixelRect(rangeLeft + (rangeButtonWidth + rangeGap) * 2, 486, rangeRight, 509), HIT_RANGE_IN, L"+", ButtonKind::Range, true, false, true, contentLocalRect(rangeLeft + (rangeButtonWidth + rangeGap) * 2, 486, rangeRight, 509));
        addFocus(FOCUS_REAL, layout.realField, contentLocalRect(valueInnerLeft, 528, valueInnerLeft + fieldWidth, 559), true, selectedEditable(working));
        addFocus(FOCUS_IMAGINARY, layout.imaginaryField, contentLocalRect(valueInnerLeft + fieldWidth + fieldGap, 528, valueInnerRight, 559), true, selectedEditable(working));
        if (scrollbar.visible()) { addFocus(FOCUS_SCROLLBAR, scrollbar.bounds(), {}, false, true); }

        int footerTop = bodyBottom;
        int footerBottom = heightDip;
        int buttonTop = footerTop + 14;
        int buttonBottom = std::min(footerBottom - 10, buttonTop + 36);
        int applyRight = widthDip - 16;
        int applyLeft = applyRight - 140;
        int revertRight = applyLeft - 9;
        int revertLeft = revertRight - 75;
        int mandelbrotRight = revertLeft - 9;
        int mandelbrotLeft = mandelbrotRight - 113;
        addButton(pixelRect(mandelbrotLeft, buttonTop, mandelbrotRight, buttonBottom), HIT_MANDELBROT, L"Mandelbrot", ButtonKind::Normal, true, false, false);
        addButton(pixelRect(revertLeft, buttonTop, revertRight, buttonBottom), HIT_REVERT, L"Revert", ButtonKind::Normal, true, false, false);
        addButton(pixelRect(applyLeft, buttonTop, applyRight, buttonBottom), HIT_APPLY, L"Apply & render", ButtonKind::Primary, true, false, false);

        auto focusRank = [](int id) {
            if (id == FOCUS_FORMULA) return 0;
            if (id >= HIT_VARIABLE_BASE && id < HIT_VARIABLE_BASE + 12) return 100 + id - HIT_VARIABLE_BASE;
            if (id >= HIT_FUNCTION_TAB_BASE && id < HIT_FUNCTION_TAB_BASE + 3) return 200 + id - HIT_FUNCTION_TAB_BASE;
            if (id >= HIT_FUNCTION_BASE && id < HIT_FUNCTION_BASE + 16) return 300 + id - HIT_FUNCTION_BASE;
            if (id == HIT_C_PLANE) return 400;
            if (id == HIT_Z0_PLANE) return 401;
            if (id == FOCUS_BAILOUT) return 402;
            if (id == FOCUS_PRESET) return 500;
            if (id == HIT_COPY) return 501;
            if (id == HIT_PASTE) return 502;
            if (id == HIT_PICKER) return 600;
            if (id == HIT_RANGE_OUT) return 601;
            if (id == HIT_RANGE_RESET) return 602;
            if (id == HIT_RANGE_IN) return 603;
            if (id == FOCUS_REAL) return 604;
            if (id == FOCUS_IMAGINARY) return 605;
            if (id == FOCUS_SCROLLBAR) return 700;
            if (id == HIT_MANDELBROT) return 800;
            if (id == HIT_REVERT) return 801;
            if (id == HIT_APPLY) return 802;
            if (id == HIT_CLOSE) return 803;
            return 900;
        };
        std::stable_sort(focusEntries.begin(), focusEntries.end(), [&focusRank](const FocusEntry& a, const FocusEntry& b) { return focusRank(a.id) < focusRank(b.id); });
    }

    const ButtonSpec* findButton(int hit) const {
        for (const ButtonSpec& button : buttons)
            if (button.hit == hit) return &button;
        return nullptr;
    }

    RECT accessibleScreenRect(RECT rect, bool content, bool& visible) const {
        RECT result = content ? intersected(rect, layout.contentClip) : rect;
        visible = hwnd && IsWindowVisible(hwnd) && nonempty(result);
        if (!nonempty(result)) result = rect;
        if (!nonempty(result)) return {};
        POINT corners[2] = {{result.left, result.top}, {result.right, result.bottom}};
        MapWindowPoints(hwnd, nullptr, corners, 2);
        return {corners[0].x, corners[0].y, corners[1].x, corners[1].y};
    }

    void appendAccessible(std::vector<FormulaAccessibleItem>& items, LONG key, std::wstring name, LONG role, RECT rect, bool content, bool enabled, bool selectedState, bool focusedState, std::wstring value, std::wstring description, std::wstring defaultAction, bool forceInvisible = false) const {
        bool visible = false;
        RECT screen = accessibleScreenRect(rect, content, visible);
        LONG state = STATE_SYSTEM_FOCUSABLE;
        if (role == ROLE_SYSTEM_RADIOBUTTON || role == ROLE_SYSTEM_PAGETAB || role == ROLE_SYSTEM_LISTITEM) { state |= STATE_SYSTEM_SELECTABLE; }
        if (!enabled) state |= STATE_SYSTEM_UNAVAILABLE;
        if (selectedState) state |= STATE_SYSTEM_SELECTED;
        if (focusedState && visible) state |= STATE_SYSTEM_FOCUSED;
        bool panelVisible = hwnd && IsWindowVisible(hwnd);
        if (!panelVisible || forceInvisible)
            state |= STATE_SYSTEM_INVISIBLE;
        else if (!visible)
            state |= STATE_SYSTEM_OFFSCREEN;
        if (forceInvisible) screen = {};
        items.push_back({key, std::move(name), std::move(value), std::move(description), std::move(defaultAction), role, state, screen});
    }

    void accessibleSnapshot(std::vector<FormulaAccessibleItem>& items) const {
        items.clear();
        auto buttonRect = [this](int hit) {
            const ButtonSpec* button = findButton(hit);
            return button ? button->rect : RECT{};
        };
        HWND keyboardFocus = GetFocus();
        bool panelHasFocus = keyboardFocus == hwnd || (keyboardFocus && IsChild(hwnd, keyboardFocus));
        auto focused = [this, panelHasFocus](int id) { return panelHasFocus && focusedId == id; };

        appendAccessible(items, FORMULA_ACC_EXPRESSION, L"Formula expression", ROLE_SYSTEM_TEXT, layout.formulaField, true, true, false, focused(FOCUS_FORMULA), formulaField.text(), L"Orbit expression. Editing stages the formula until Apply.", L"Edit");
        int presetIndex = presetDropdown.selectedIndex();
        std::wstring presetValue = presetIndex >= 0 && presetIndex < static_cast<int>(presets().size()) ? presets()[static_cast<size_t>(presetIndex)].label : L"";
        appendAccessible(items, FORMULA_ACC_PRESET, L"Formula preset", ROLE_SYSTEM_COMBOBOX, layout.preset, true, true, false, focused(FOCUS_PRESET), std::move(presetValue), presetDropdown.open() ? L"Preset list expanded. Activate to accept the current item." : L"Preset list collapsed. Activate to open it.", presetDropdown.open() ? L"Accept" : L"Open");
        if (presetDropdown.open())
            items.back().state |= STATE_SYSTEM_EXPANDED;
        else
            items.back().state |= STATE_SYSTEM_COLLAPSED;

        appendAccessible(items, FORMULA_ACC_COPY, L"Copy complete formula", ROLE_SYSTEM_PUSHBUTTON, buttonRect(HIT_COPY), true, true, false, focused(HIT_COPY), L"", L"Copies the entire formula expression.", L"Press");
        appendAccessible(items, FORMULA_ACC_PASTE, L"Paste complete formula", ROLE_SYSTEM_PUSHBUTTON, buttonRect(HIT_PASTE), true, true, false, focused(HIT_PASTE), L"", L"Replaces the entire formula from the clipboard.", L"Press");

        for (int i = 0; i < 12; ++i) {
            InspectorValue value = variableFromButton(static_cast<size_t>(i));
            std::wstring label = variableLabels()[static_cast<size_t>(i)];
            std::wstring description = L"Inserts variable " + label + L" and selects its value inspector.";
            appendAccessible(items, FORMULA_ACC_VARIABLE_BASE + i, L"Variable " + label, ROLE_SYSTEM_PUSHBUTTON, buttonRect(HIT_VARIABLE_BASE + i), true, true, selected == value, focused(HIT_VARIABLE_BASE + i), usageCount(value) == 0 ? L"not used" : std::to_wstring(usageCount(value)) + L" uses", std::move(description), L"Insert");
        }

        static const std::array<const wchar_t*, 3> tabLabels = {L"Common", L"Complex", L"All"};
        for (int i = 0; i < 3; ++i) { appendAccessible(items, FORMULA_ACC_TAB_BASE + i, std::wstring(tabLabels[static_cast<size_t>(i)]) + L" functions tab", ROLE_SYSTEM_PAGETAB, buttonRect(HIT_FUNCTION_TAB_BASE + i), true, true, functionTab == i, focused(HIT_FUNCTION_TAB_BASE + i), L"", L"Selects the visible function button set.", L"Select"); }

        const auto& labels = functionSets()[static_cast<size_t>(functionTab)];
        for (int i = 0; i < 16; ++i) {
            std::wstring label = labels[static_cast<size_t>(i)];
            appendAccessible(items, FORMULA_ACC_FUNCTION_BASE + i, L"Function " + label, ROLE_SYSTEM_PUSHBUTTON, buttonRect(HIT_FUNCTION_BASE + i), true, true, false, focused(HIT_FUNCTION_BASE + i), L"", L"Inserts this function into the formula.", L"Insert");
        }

        bool cPlane = working.pixelParameter != FormulaParameter::InitialZ;
        appendAccessible(items, FORMULA_ACC_C_PLANE, L"c-plane", ROLE_SYSTEM_RADIOBUTTON, buttonRect(HIT_C_PLANE), true, true, cPlane, focused(HIT_C_PLANE), L"", L"Each pixel supplies c; z0 is a fixed editable value.", L"Select");
        appendAccessible(items, FORMULA_ACC_Z0_PLANE, L"z0-plane", ROLE_SYSTEM_RADIOBUTTON, buttonRect(HIT_Z0_PLANE), true, true, !cPlane, focused(HIT_Z0_PLANE), L"", L"Each pixel supplies z0; c is a fixed editable value.", L"Select");
        appendAccessible(items, FORMULA_ACC_BAILOUT, L"Bailout value", ROLE_SYSTEM_TEXT, layout.bailout, true, true, false, focused(FOCUS_BAILOUT), bailoutField.text(), L"Positive escape radius.", L"Edit");

        const std::complex<double>* complexValue = selectedValue(working);
        std::wstring pickerValue = complexValue ? formatNumber(complexValue->real()) + L", " + formatNumber(complexValue->imag()) : L"Read only";
        bool pickerEnabled = selectedEditable(working);
        appendAccessible(items, FORMULA_ACC_PICKER, L"Complex value picker", ROLE_SYSTEM_SLIDER, layout.picker, true, pickerEnabled, false, focused(HIT_PICKER), std::move(pickerValue),
                         L"Complex plane picker. Arrow keys move by one fiftieth "
                         L"of the current range.",
                         L"Focus");
        appendAccessible(items, FORMULA_ACC_RANGE_OUT, L"Increase picker range", ROLE_SYSTEM_PUSHBUTTON, buttonRect(HIT_RANGE_OUT), true, true, false, focused(HIT_RANGE_OUT), formatNumber(pickerRange), L"Doubles the complex picker range.", L"Press");
        appendAccessible(items, FORMULA_ACC_RANGE_RESET, L"Reset picker range", ROLE_SYSTEM_PUSHBUTTON, buttonRect(HIT_RANGE_RESET), true, true, false, focused(HIT_RANGE_RESET), formatNumber(pickerRange), L"Resets the picker range to 2.", L"Press");
        appendAccessible(items, FORMULA_ACC_RANGE_IN, L"Decrease picker range", ROLE_SYSTEM_PUSHBUTTON, buttonRect(HIT_RANGE_IN), true, true, false, focused(HIT_RANGE_IN), formatNumber(pickerRange), L"Halves the complex picker range.", L"Press");
        appendAccessible(items, FORMULA_ACC_REAL, L"Real value", ROLE_SYSTEM_TEXT, layout.realField, true, realField.enabled(), false, focused(FOCUS_REAL), realField.text(), L"Real component of the selected fixed complex value.", L"Edit");
        appendAccessible(items, FORMULA_ACC_IMAGINARY, L"Imaginary value", ROLE_SYSTEM_TEXT, layout.imaginaryField, true, imaginaryField.enabled(), false, focused(FOCUS_IMAGINARY), imaginaryField.text(), L"Imaginary component of the selected fixed complex value.", L"Edit");

        wchar_t scrollValue[64]{};
        swprintf_s(scrollValue, L"%d of %d", scrollbar.position(), scrollbar.maximumPosition());
        appendAccessible(items, FORMULA_ACC_SCROLLBAR, L"Formula editor scrollbar", ROLE_SYSTEM_SCROLLBAR, scrollbar.bounds(), false, scrollbar.visible(), false, focused(FOCUS_SCROLLBAR), scrollValue,
                         L"Vertical content position. Arrow keys scroll one line; "
                         L"Page Up and Page Down scroll one page.",
                         L"Page down", !scrollbar.visible());

        appendAccessible(items, FORMULA_ACC_MANDELBROT, L"Mandelbrot", ROLE_SYSTEM_PUSHBUTTON, buttonRect(HIT_MANDELBROT), false, true, false, focused(HIT_MANDELBROT), L"", L"Switches back to the built-in Mandelbrot formula.", L"Press");
        appendAccessible(items, FORMULA_ACC_REVERT, L"Revert", ROLE_SYSTEM_PUSHBUTTON, buttonRect(HIT_REVERT), false, true, false, focused(HIT_REVERT), L"", L"Restores the last applied formula configuration.", L"Press");
        appendAccessible(items, FORMULA_ACC_APPLY, L"Apply and render", ROLE_SYSTEM_PUSHBUTTON, buttonRect(HIT_APPLY), false, true, false, focused(HIT_APPLY), L"", L"Validates and applies staged changes, then renders.", L"Press");
        appendAccessible(items, FORMULA_ACC_CLOSE, L"Close formula editor", ROLE_SYSTEM_PUSHBUTTON, buttonRect(HIT_CLOSE), false, true, false, focused(HIT_CLOSE), L"", L"Closes the Formula editor dock.", L"Press");

        for (int i = 0; i < static_cast<int>(presets().size()); ++i) {
            RECT itemRect = presetDropdown.itemBounds(i, layout.contentClip);
            appendAccessible(items, FORMULA_ACC_PRESET_ITEM_BASE + i, std::wstring(L"Preset ") + presets()[static_cast<size_t>(i)].label, ROLE_SYSTEM_LISTITEM, itemRect, false, true, presetIndex == i, false, presets()[static_cast<size_t>(i)].label, L"Selects this formula preset.", L"Choose", !presetDropdown.open());
        }

        auto accessibleRank = [](LONG key) {
            if (key == FORMULA_ACC_EXPRESSION) return 0L;
            if (key >= FORMULA_ACC_VARIABLE_BASE && key < FORMULA_ACC_VARIABLE_BASE + 12) return 100L + key - FORMULA_ACC_VARIABLE_BASE;
            if (key >= FORMULA_ACC_TAB_BASE && key < FORMULA_ACC_TAB_BASE + 3) return 200L + key - FORMULA_ACC_TAB_BASE;
            if (key >= FORMULA_ACC_FUNCTION_BASE && key < FORMULA_ACC_FUNCTION_BASE + 16) return 300L + key - FORMULA_ACC_FUNCTION_BASE;
            if (key == FORMULA_ACC_C_PLANE) return 400L;
            if (key == FORMULA_ACC_Z0_PLANE) return 401L;
            if (key == FORMULA_ACC_BAILOUT) return 402L;
            if (key == FORMULA_ACC_PRESET) return 500L;
            if (key == FORMULA_ACC_COPY) return 501L;
            if (key == FORMULA_ACC_PASTE) return 502L;
            if (key >= FORMULA_ACC_PRESET_ITEM_BASE) return 510L + key - FORMULA_ACC_PRESET_ITEM_BASE;
            if (key == FORMULA_ACC_PICKER) return 600L;
            if (key == FORMULA_ACC_RANGE_OUT) return 601L;
            if (key == FORMULA_ACC_RANGE_RESET) return 602L;
            if (key == FORMULA_ACC_RANGE_IN) return 603L;
            if (key == FORMULA_ACC_REAL) return 604L;
            if (key == FORMULA_ACC_IMAGINARY) return 605L;
            if (key == FORMULA_ACC_SCROLLBAR) return 700L;
            if (key == FORMULA_ACC_MANDELBROT) return 800L;
            if (key == FORMULA_ACC_REVERT) return 801L;
            if (key == FORMULA_ACC_APPLY) return 802L;
            if (key == FORMULA_ACC_CLOSE) return 803L;
            return 900L;
        };
        std::stable_sort(items.begin(), items.end(), [&accessibleRank](const FormulaAccessibleItem& a, const FormulaAccessibleItem& b) { return accessibleRank(a.key) < accessibleRank(b.key); });
    }

    LONG accessibleKeyForFocus(int id) const {
        if (id == FOCUS_FORMULA) return FORMULA_ACC_EXPRESSION;
        if (id == FOCUS_PRESET) return FORMULA_ACC_PRESET;
        if (id == FOCUS_BAILOUT) return FORMULA_ACC_BAILOUT;
        if (id == FOCUS_REAL) return FORMULA_ACC_REAL;
        if (id == FOCUS_IMAGINARY) return FORMULA_ACC_IMAGINARY;
        if (id == FOCUS_SCROLLBAR) return FORMULA_ACC_SCROLLBAR;
        if (id == HIT_CLOSE) return FORMULA_ACC_CLOSE;
        if (id == HIT_COPY) return FORMULA_ACC_COPY;
        if (id == HIT_PASTE) return FORMULA_ACC_PASTE;
        if (id == HIT_C_PLANE) return FORMULA_ACC_C_PLANE;
        if (id == HIT_Z0_PLANE) return FORMULA_ACC_Z0_PLANE;
        if (id == HIT_PICKER) return FORMULA_ACC_PICKER;
        if (id == HIT_RANGE_OUT) return FORMULA_ACC_RANGE_OUT;
        if (id == HIT_RANGE_RESET) return FORMULA_ACC_RANGE_RESET;
        if (id == HIT_RANGE_IN) return FORMULA_ACC_RANGE_IN;
        if (id == HIT_MANDELBROT) return FORMULA_ACC_MANDELBROT;
        if (id == HIT_REVERT) return FORMULA_ACC_REVERT;
        if (id == HIT_APPLY) return FORMULA_ACC_APPLY;
        if (id >= HIT_VARIABLE_BASE && id < HIT_VARIABLE_BASE + 12) { return FORMULA_ACC_VARIABLE_BASE + id - HIT_VARIABLE_BASE; }
        if (id >= HIT_FUNCTION_TAB_BASE && id < HIT_FUNCTION_TAB_BASE + 3) { return FORMULA_ACC_TAB_BASE + id - HIT_FUNCTION_TAB_BASE; }
        if (id >= HIT_FUNCTION_BASE && id < HIT_FUNCTION_BASE + 16) { return FORMULA_ACC_FUNCTION_BASE + id - HIT_FUNCTION_BASE; }
        return FORMULA_ACC_SELF;
    }

    void notifyDropdownTransition(bool wasOpen) {
        bool isOpen = presetDropdown.open();
        if (wasOpen == isOpen) return;
        notifyAccessibility(isOpen ? EVENT_SYSTEM_MENUPOPUPSTART : EVENT_SYSTEM_MENUPOPUPEND, FORMULA_ACC_PRESET);
        notifyAccessibility(EVENT_OBJECT_STATECHANGE, FORMULA_ACC_PRESET);
    }

    bool accessibilityFocusKey(LONG key, bool selectAction) {
        if (key == FORMULA_ACC_SELF) {
            SetFocus(hwnd);
            return true;
        }
        if (key == FORMULA_ACC_EXPRESSION) {
            setFocused(FOCUS_FORMULA);
            if (selectAction) formulaField.selectAll();
            return true;
        }
        if (key == FORMULA_ACC_PRESET) {
            setFocused(FOCUS_PRESET);
            return true;
        }
        if (key == FORMULA_ACC_BAILOUT) {
            setFocused(FOCUS_BAILOUT);
            if (selectAction) bailoutField.selectAll();
            return true;
        }
        if (key == FORMULA_ACC_REAL) {
            if (!realField.enabled()) return false;
            setFocused(FOCUS_REAL);
            if (selectAction) realField.selectAll();
            return true;
        }
        if (key == FORMULA_ACC_IMAGINARY) {
            if (!imaginaryField.enabled()) return false;
            setFocused(FOCUS_IMAGINARY);
            if (selectAction) imaginaryField.selectAll();
            return true;
        }
        if (key == FORMULA_ACC_SCROLLBAR) {
            if (!scrollbar.visible()) return false;
            setFocused(FOCUS_SCROLLBAR);
            return true;
        }
        if (key >= FORMULA_ACC_PRESET_ITEM_BASE && key < FORMULA_ACC_PRESET_ITEM_BASE + static_cast<LONG>(presets().size())) {
            setFocused(FOCUS_PRESET);
            if (selectAction) {
                int index = static_cast<int>(key - FORMULA_ACC_PRESET_ITEM_BASE);
                bool wasOpen = presetDropdown.open();
                presetDropdown.setSelectedIndex(index);
                choosePreset(index);
                presetDropdown.close();
                notifyDropdownTransition(wasOpen);
                notifyAccessibility(EVENT_OBJECT_SELECTION, FORMULA_ACC_PRESET_ITEM_BASE + index);
                notifyAccessibility(EVENT_OBJECT_VALUECHANGE, FORMULA_ACC_PRESET);
            }
            return true;
        }

        int focus = HIT_NONE;
        if (key == FORMULA_ACC_COPY)
            focus = HIT_COPY;
        else if (key == FORMULA_ACC_PASTE)
            focus = HIT_PASTE;
        else if (key == FORMULA_ACC_C_PLANE)
            focus = HIT_C_PLANE;
        else if (key == FORMULA_ACC_Z0_PLANE)
            focus = HIT_Z0_PLANE;
        else if (key == FORMULA_ACC_PICKER) {
            if (!selectedEditable(working)) return false;
            focus = HIT_PICKER;
        } else if (key == FORMULA_ACC_RANGE_OUT)
            focus = HIT_RANGE_OUT;
        else if (key == FORMULA_ACC_RANGE_RESET)
            focus = HIT_RANGE_RESET;
        else if (key == FORMULA_ACC_RANGE_IN)
            focus = HIT_RANGE_IN;
        else if (key == FORMULA_ACC_MANDELBROT)
            focus = HIT_MANDELBROT;
        else if (key == FORMULA_ACC_REVERT)
            focus = HIT_REVERT;
        else if (key == FORMULA_ACC_APPLY)
            focus = HIT_APPLY;
        else if (key == FORMULA_ACC_CLOSE)
            focus = HIT_CLOSE;
        else if (key >= FORMULA_ACC_VARIABLE_BASE && key < FORMULA_ACC_VARIABLE_BASE + 12) {
            focus = HIT_VARIABLE_BASE + static_cast<int>(key - FORMULA_ACC_VARIABLE_BASE);
            if (selectAction) { selectInspector(variableFromButton(static_cast<size_t>(key - FORMULA_ACC_VARIABLE_BASE))); }
        } else if (key >= FORMULA_ACC_TAB_BASE && key < FORMULA_ACC_TAB_BASE + 3) {
            focus = HIT_FUNCTION_TAB_BASE + static_cast<int>(key - FORMULA_ACC_TAB_BASE);
            if (selectAction) {
                int oldTab = functionTab;
                functionTab = static_cast<int>(key - FORMULA_ACC_TAB_BASE);
                updateLayout();
                InvalidateRect(hwnd, nullptr, FALSE);
                notifyAccessibility(EVENT_OBJECT_SELECTION, FORMULA_ACC_TAB_BASE + functionTab);
                if (oldTab != functionTab) notifyAccessibility(EVENT_OBJECT_REORDER, FORMULA_ACC_SELF);
            }

        } else if (key >= FORMULA_ACC_FUNCTION_BASE && key < FORMULA_ACC_FUNCTION_BASE + 16) {
            focus = HIT_FUNCTION_BASE + static_cast<int>(key - FORMULA_ACC_FUNCTION_BASE);
        }
        if (focus == HIT_NONE) return false;
        setFocused(focus);
        if (selectAction && key == FORMULA_ACC_C_PLANE)
            choosePlane(FormulaParameter::C);
        else if (selectAction && key == FORMULA_ACC_Z0_PLANE)
            choosePlane(FormulaParameter::InitialZ);
        return true;
    }

    bool accessibilityDefaultAction(LONG key) {
        if (key == FORMULA_ACC_EXPRESSION || key == FORMULA_ACC_BAILOUT || key == FORMULA_ACC_REAL || key == FORMULA_ACC_IMAGINARY) { return accessibilityFocusKey(key, true); }
        if (key == FORMULA_ACC_PRESET) {
            setFocused(FOCUS_PRESET);
            bool wasOpen = presetDropdown.open();
            presetDropdown.keyDown(VK_RETURN, false, layout.contentClip);
            notifyDropdownTransition(wasOpen);
            InvalidateRect(hwnd, nullptr, FALSE);
            return true;
        }
        if (key >= FORMULA_ACC_PRESET_ITEM_BASE && key < FORMULA_ACC_PRESET_ITEM_BASE + static_cast<LONG>(presets().size())) { return accessibilityFocusKey(key, true); }
        if (key == FORMULA_ACC_PICKER) return accessibilityFocusKey(key, false);
        if (key == FORMULA_ACC_SCROLLBAR) {
            if (!scrollbar.visible()) return false;
            int page = std::max(48, clientHeightDip() - HEADER_HEIGHT - FOOTER_HEIGHT - BODY_PADDING * 2);
            scrollbar.scrollBy(page);
            return true;
        }

        int hit = HIT_NONE;
        if (key == FORMULA_ACC_COPY)
            hit = HIT_COPY;
        else if (key == FORMULA_ACC_PASTE)
            hit = HIT_PASTE;
        else if (key == FORMULA_ACC_C_PLANE)
            hit = HIT_C_PLANE;
        else if (key == FORMULA_ACC_Z0_PLANE)
            hit = HIT_Z0_PLANE;
        else if (key == FORMULA_ACC_RANGE_OUT)
            hit = HIT_RANGE_OUT;
        else if (key == FORMULA_ACC_RANGE_RESET)
            hit = HIT_RANGE_RESET;
        else if (key == FORMULA_ACC_RANGE_IN)
            hit = HIT_RANGE_IN;
        else if (key == FORMULA_ACC_MANDELBROT)
            hit = HIT_MANDELBROT;
        else if (key == FORMULA_ACC_REVERT)
            hit = HIT_REVERT;
        else if (key == FORMULA_ACC_APPLY)
            hit = HIT_APPLY;
        else if (key == FORMULA_ACC_CLOSE)
            hit = HIT_CLOSE;
        else if (key >= FORMULA_ACC_VARIABLE_BASE && key < FORMULA_ACC_VARIABLE_BASE + 12) {
            hit = HIT_VARIABLE_BASE + static_cast<int>(key - FORMULA_ACC_VARIABLE_BASE);
        } else if (key >= FORMULA_ACC_TAB_BASE && key < FORMULA_ACC_TAB_BASE + 3) {
            hit = HIT_FUNCTION_TAB_BASE + static_cast<int>(key - FORMULA_ACC_TAB_BASE);
        } else if (key >= FORMULA_ACC_FUNCTION_BASE && key < FORMULA_ACC_FUNCTION_BASE + 16) {
            hit = HIT_FUNCTION_BASE + static_cast<int>(key - FORMULA_ACC_FUNCTION_BASE);
        }
        if (hit == HIT_NONE) return false;
        setFocused(hit);
        invokeHit(hit);
        return true;
    }

    bool performAccessibilityAction(LONG key, FormulaAccessibilityAction action) {
        if (action == FormulaAccessibilityAction::Default) return accessibilityDefaultAction(key);
        return accessibilityFocusKey(key, action == FormulaAccessibilityAction::Select);
    }

    void ensureContentVisible(RECT contentRect) {
        int viewport = std::max(1, clientHeightDip() - HEADER_HEIGHT - FOOTER_HEIGHT - BODY_PADDING * 2);
        int position = scrollbar.position();
        if (contentRect.top < position)
            position = contentRect.top;
        else if (contentRect.bottom > position + viewport)
            position = contentRect.bottom - viewport;
        scrollbar.setPosition(position);
    }

    const FocusEntry* findFocusEntry(int id) const {
        for (const FocusEntry& entry : focusEntries)
            if (entry.id == id) return &entry;
        return nullptr;
    }

    ui::TextField* textFieldForFocus(int id) {
        if (id == FOCUS_FORMULA) return &formulaField;
        if (id == FOCUS_BAILOUT) return &bailoutField;
        if (id == FOCUS_REAL) return &realField;
        if (id == FOCUS_IMAGINARY) return &imaginaryField;
        return nullptr;
    }

    void setFocused(int id) {
        bool dropdownWasOpen = presetDropdown.open();
        const FocusEntry* entry = findFocusEntry(id);
        if (entry && entry->content) {
            RECT local = entry->contentRect;
            ensureContentVisible(local);
            updateLayout();
        }
        focusedId = id;
        presetDropdown.setFocused(id == FOCUS_PRESET);
        if (id != FOCUS_PRESET) presetDropdown.close();

        ui::TextField* field = textFieldForFocus(id);
        if (field && field->enabled()) {
            field->focus();
        } else if (hwnd) {
            SetFocus(hwnd);
        }
        notifyDropdownTransition(dropdownWasOpen);
        notifyAccessibility(EVENT_OBJECT_FOCUS, accessibleKeyForFocus(id));
        if (hwnd) InvalidateRect(hwnd, nullptr, FALSE);
    }

    void focusNext(bool reverse) {
        if (focusEntries.empty()) return;
        int current = -1;
        for (size_t i = 0; i < focusEntries.size(); ++i) {
            if (focusEntries[i].id == focusedId) {
                current = static_cast<int>(i);
                break;
            }
        }
        int count = static_cast<int>(focusEntries.size());
        for (int step = 1; step <= count; ++step) {
            int index = current < 0 ? (reverse ? count - 1 : 0) : (current + (reverse ? -step : step)) % count;
            if (index < 0) index += count;
            if (focusEntries[static_cast<size_t>(index)].enabled) {
                setFocused(focusEntries[static_cast<size_t>(index)].id);
                return;
            }
        }
    }

    bool commitInspector(bool reportErrors) {
        if (!selectedEditable(working)) return true;
        double real = 0.0;
        double imaginary = 0.0;
        if (!parseFiniteNumber(realField.text(), real)) {
            if (reportErrors) {
                setStatus(L"Real value must be a finite number.", true);
                setFocused(FOCUS_REAL);
                realField.selectAll();
            }
            return false;
        }
        if (!parseFiniteNumber(imaginaryField.text(), imaginary)) {
            if (reportErrors) {
                setStatus(L"Imaginary value must be a finite number.", true);
                setFocused(FOCUS_IMAGINARY);
                imaginaryField.selectAll();
            }
            return false;
        }
        std::complex<double>* value = selectedValue(working);
        if (value) *value = {real, imaginary};
        return true;
    }

    bool selectInspector(InspectorValue value, bool reportErrors = true) {
        if (selected == value) {
            if (hwnd) InvalidateRect(hwnd, nullptr, FALSE);
            return true;
        }
        if (!commitInspector(reportErrors)) return false;
        InspectorValue previous = selected;
        selected = value;
        syncInspectorFields();
        updateLayout();
        notifyAccessibility(EVENT_OBJECT_SELECTION, FORMULA_ACC_VARIABLE_BASE + static_cast<LONG>(value));
        notifyAccessibility(EVENT_OBJECT_STATECHANGE, FORMULA_ACC_VARIABLE_BASE + static_cast<LONG>(previous));
        notifyAccessibility(EVENT_OBJECT_VALUECHANGE, FORMULA_ACC_REAL);
        notifyAccessibility(EVENT_OBJECT_VALUECHANGE, FORMULA_ACC_IMAGINARY);
        notifyAccessibility(EVENT_OBJECT_VALUECHANGE, FORMULA_ACC_PICKER);
        return true;
    }

    const Token* tokenAt(size_t position) const {
        for (const Token& token : tokens) {
            if (position >= token.first && position < token.last) return &token;
        }
        if (position > 0) {
            for (const Token& token : tokens) {
                if (position == token.last && token.first < token.last) return &token;
            }
        }
        return nullptr;
    }

    int tokenIndexAtPoint(POINT point) const {
        std::vector<ui::TextRangeStyle> styles = formulaStyles();
        for (size_t i = 0; i < tokens.size(); ++i) {
            const Token& token = tokens[i];
            if (!token.selectable) continue;
            int leading = i < styles.size() ? styles[i].paddingBefore : 0;
            RECT bounds = formulaField.textRangeBounds(token.first, token.last, leading);
            if (containsPoint(bounds, point)) return static_cast<int>(i);
        }
        return -1;
    }

    void selectInspectorAtCaret(bool selectToken) {
        auto selection = formulaField.selection();
        size_t position = selection.first != selection.second ? selection.first : selection.second;
        const Token* token = tokenAt(position);
        if (!token || !token->selectable) return;
        if (!selectInspector(token->inspector, false)) return;
        if (selectToken) formulaField.setSelection(token->first, token->last);
    }

    bool insertFormulaText(const std::wstring& token, bool function) {
        auto selection = formulaField.selection();
        std::wstring current = formulaField.text();
        size_t first = std::min(selection.first, current.size());
        size_t last = std::min(selection.second, current.size());
        std::wstring insertion;
        size_t selectionFirst = 0;
        size_t selectionLast = 0;

        if (function) {
            if (last > first) {
                insertion = token + L"(" + current.substr(first, last - first) + L")";
                selectionFirst = first + insertion.size();
                selectionLast = selectionFirst;
            } else {
                insertion = token + L"()";
                selectionFirst = first + token.size() + 1u;
                selectionLast = selectionFirst;
            }
        } else {
            insertion = token;
            selectionFirst = first + insertion.size();
            selectionLast = selectionFirst;
        }

        if (current.size() - (last - first) + insertion.size() > formula::ExpressionProgram::MAX_SOURCE) {
            setStatus(L"Insertion failed: the formula is too long.", true);
            return false;
        }
        formulaField.setSelection(first, last);
        if (!formulaField.replaceSelection(insertion, true)) {
            setStatus(L"Insertion failed.", true);
            return false;
        }
        formulaField.setSelection(selectionFirst, selectionLast);
        setFocused(FOCUS_FORMULA);
        return true;
    }

    void choosePreset(int index) {
        if (index <= 0 || index >= static_cast<int>(presets().size())) {
            presetDropdown.setSelectedIndex(0);
            return;
        }
        const Preset& preset = presets()[static_cast<size_t>(index)];
        std::wstring source;
        if (!widenUtf8(preset.source, source)) {
            setStatus(L"Preset text is not valid UTF-8.", true);
            return;
        }
        working.source = preset.source;
        if (preset.setP0) working.parameters[0] = preset.p0;
        syncing = true;
        formulaField.setText(source);
        formulaField.setSelection(source.size(), source.size());
        presetDropdown.setSelectedIndex(index);
        syncing = false;
        refreshFormulaAnalysis();
        if (selected == InspectorValue::P0 && preset.setP0) syncInspectorFields();
        setStatus(std::wstring(preset.label) + L" preset loaded; changes are staged.", false);
        notifyAccessibility(EVENT_OBJECT_SELECTION, FORMULA_ACC_PRESET_ITEM_BASE + index);
        notifyAccessibility(EVENT_OBJECT_VALUECHANGE, FORMULA_ACC_PRESET);
        notifyAccessibility(EVENT_OBJECT_VALUECHANGE, FORMULA_ACC_EXPRESSION);
        setFocused(FOCUS_FORMULA);
        triggerPreview();
    }

    void copyCompleteFormula() {
        if (!writeClipboard(hwnd, formulaField.text())) {
            setStatus(L"Copy failed because the clipboard is unavailable.", true);
            return;
        }
        setStatus(L"Complete formula copied to the clipboard.", false);
    }

    void pasteCompleteFormula() {
        std::wstring clipboard;
        if (!readClipboard(hwnd, clipboard)) {
            setStatus(L"Paste failed because the clipboard has no readable text.", true);
            return;
        }
        if (clipboard.size() > formula::ExpressionProgram::MAX_SOURCE) {
            setStatus(L"Paste failed: the formula is longer than 4096 characters.", true);
            return;
        }
        formulaField.selectAll();
        if (!formulaField.replaceSelection(clipboard, true)) {
            setStatus(L"Paste failed.", true);
            return;
        }
        formulaField.setSelection(clipboard.size(), clipboard.size());
        setStatus(L"Complete formula replaced from the clipboard.", false);
        setFocused(FOCUS_FORMULA);
    }

    void choosePlane(FormulaParameter plane) {
        if (working.pixelParameter == plane) return;
        if (!commitInspector(true)) return;
        working.pixelParameter = plane;
        syncInspectorFields();
        updateLayout();
        setStatus(plane == FormulaParameter::InitialZ ? L"z0 plane selected; fixed c is editable." : L"c plane selected; fixed z0 is editable.", false);
        notifyAccessibility(EVENT_OBJECT_SELECTION, plane == FormulaParameter::InitialZ ? FORMULA_ACC_Z0_PLANE : FORMULA_ACC_C_PLANE);
        notifyAccessibility(EVENT_OBJECT_VALUECHANGE, FORMULA_ACC_REAL);
        notifyAccessibility(EVENT_OBJECT_VALUECHANGE, FORMULA_ACC_IMAGINARY);
        notifyAccessibility(EVENT_OBJECT_STATECHANGE, FORMULA_ACC_PICKER);
        triggerPreview();
    }

    void changeRange(double factor) {
        pickerRange = std::clamp(pickerRange * factor, 1.0e-12, 1.0e12);
        if (hwnd) InvalidateRect(hwnd, nullptr, FALSE);
        notifyAccessibility(EVENT_OBJECT_VALUECHANGE, FORMULA_ACC_PICKER);
    }

    void updatePickerFromMouse(POINT point) {
        if (!selectedEditable(working)) {
            setStatus(inspectorDescription(), true);
            return;
        }
        RECT rect = layout.picker;
        double width = std::max(1L, rect.right - rect.left - 1);
        double height = std::max(1L, rect.bottom - rect.top - 1);
        double x = std::clamp((point.x - rect.left) / width, 0.0, 1.0);
        double y = std::clamp((point.y - rect.top) / height, 0.0, 1.0);
        double real = (x * 2.0 - 1.0) * pickerRange;
        double imaginary = (1.0 - y * 2.0) * pickerRange;
        std::complex<double>* value = selectedValue(working);
        if (!value) return;
        *value = {real, imaginary};
        syncing = true;
        realField.setText(formatNumber(real, 12));
        imaginaryField.setText(formatNumber(imaginary, 12));
        syncing = false;
        setStatus(std::wstring(inspectorName(selected)) + L" changed; changes are staged.", false);
        notifyAccessibility(EVENT_OBJECT_VALUECHANGE, FORMULA_ACC_PICKER);
        notifyAccessibility(EVENT_OBJECT_VALUECHANGE, FORMULA_ACC_REAL);
        notifyAccessibility(EVENT_OBJECT_VALUECHANGE, FORMULA_ACC_IMAGINARY);
        InvalidateRect(hwnd, nullptr, FALSE);
        triggerPreview();
    }

    bool validateComplexes(FormulaDialogConfig& candidate) {
        std::array<std::pair<std::complex<double>*, InspectorValue>, 9> values{};
        values[0] = candidate.pixelParameter == FormulaParameter::InitialZ ? std::make_pair(&candidate.fixedC, InspectorValue::C) : std::make_pair(&candidate.fixedZ0, InspectorValue::Z0);
        for (int i = 0; i < 8; ++i) { values[static_cast<size_t>(i + 1)] = {&candidate.parameters[static_cast<size_t>(i)], static_cast<InspectorValue>(static_cast<int>(InspectorValue::P0) + i)}; }
        for (const auto& entry : values) {
            if (finiteComplex(*entry.first)) continue;
            selected = entry.second;
            working = candidate;
            syncInspectorFields();
            setStatus(std::wstring(inspectorName(selected)) + L" must contain finite values.", true);
            if (!std::isfinite(entry.first->real())) {
                setFocused(FOCUS_REAL);
                realField.selectAll();
            } else {
                setFocused(FOCUS_IMAGINARY);
                imaginaryField.selectAll();
            }
            return false;
        }
        return true;
    }

    bool readCandidate(FormulaDialogConfig& candidate) {
        candidate = working;
        std::wstring sourceText = formulaField.text();
        if (!narrowUtf8(sourceText, candidate.source)) {
            setStatus(L"Formula source cannot be converted to UTF-8.", true);
            setFocused(FOCUS_FORMULA);
            formulaField.selectAll();
            return false;
        }

        double bailout = 0.0;
        if (!parseFiniteNumber(bailoutField.text(), bailout) || !(bailout > 0.0)) {
            setStatus(L"Bailout must be a positive finite number.", true);
            setFocused(FOCUS_BAILOUT);
            bailoutField.selectAll();
            return false;
        }
        candidate.bailout = bailout;

        if (!commitInspector(true)) return false;
        candidate = working;
        candidate.source.clear();
        if (!narrowUtf8(sourceText, candidate.source)) return false;
        candidate.bailout = bailout;
        if (!validateComplexes(candidate)) return false;

        formula::ExpressionProgram program;
        formula::ExpressionError error;
        if (!program.compile(candidate.source, &error)) {
            std::wstring detail;
            if (!widenUtf8(error.message, detail)) detail = L"Invalid formula";
            wchar_t prefix[96]{};
            swprintf_s(prefix, L"Formula error at character %zu: ", error.position + 1u);
            setStatus(std::wstring(prefix) + detail, true);
            size_t caret = wideIndexForUtf8Offset(candidate.source, error.position);
            size_t end = std::min(caret + 1u, sourceText.size());
            setFocused(FOCUS_FORMULA);
            formulaField.setSelection(caret, end);
            return false;
        }
        return true;
    }

    void triggerPreview() {
        FormulaDialogConfig candidate;
        if (!readCandidate(candidate)) return;
        if (callbacks.stagePreview) {
            callbacks.stagePreview(candidate);
        }
    }

    void applyChanges() {
        FormulaDialogConfig candidate;
        if (!readCandidate(candidate)) return;
        if (!callbacks.apply) {
            setStatus(L"Apply is unavailable.", true);
            return;
        }
        auto apply = callbacks.apply;
        std::weak_ptr<int> guard = lifetime;
        bool accepted = apply(candidate);
        if (guard.expired()) return;
        if (!accepted) {
            setStatus(L"The formula was not accepted; staged edits were kept.", true);
            return;
        }
        working = candidate;
        applied = candidate;
        syncAllControls();
        setStatus(L"Formula applied and render requested.", false);
        notifyAccessibility(EVENT_OBJECT_VALUECHANGE, FORMULA_ACC_EXPRESSION);
        notifyAccessibility(EVENT_OBJECT_VALUECHANGE, FORMULA_ACC_BAILOUT);
    }

    void revertChanges() {
        working = applied;
        pickerRange = 2.0;
        refreshFormulaAnalysis();
        selected = initialInspector();
        syncAllControls();
        setStatus(L"Reverted to the last applied configuration.", false);
        notifyAccessibility(EVENT_OBJECT_VALUECHANGE, FORMULA_ACC_EXPRESSION);
        notifyAccessibility(EVENT_OBJECT_VALUECHANGE, FORMULA_ACC_BAILOUT);
        notifyAccessibility(EVENT_OBJECT_VALUECHANGE, FORMULA_ACC_REAL);
        notifyAccessibility(EVENT_OBJECT_VALUECHANGE, FORMULA_ACC_IMAGINARY);
    }

    void useMandelbrot() {
        auto callback = callbacks.useMandelbrot;
        if (callback) callback();
    }

    void closePanel() {
        auto callback = callbacks.close;
        hide();
        if (callback) callback();
    }

    void invokeHit(int hit) {
        if (hit == HIT_CLOSE) {
            closePanel();
            return;
        }
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
        if (hit == HIT_RANGE_OUT) {
            changeRange(2.0);
            return;
        }
        if (hit == HIT_RANGE_RESET) {
            pickerRange = 2.0;
            InvalidateRect(hwnd, nullptr, FALSE);
            notifyAccessibility(EVENT_OBJECT_VALUECHANGE, FORMULA_ACC_PICKER);
            return;
        }
        if (hit == HIT_RANGE_IN) {
            changeRange(0.5);
            return;
        }
        if (hit == HIT_MANDELBROT) {
            useMandelbrot();
            return;
        }
        if (hit == HIT_REVERT) {
            revertChanges();
            return;
        }
        if (hit == HIT_APPLY) {
            applyChanges();
            return;
        }
        if (hit >= HIT_FUNCTION_TAB_BASE && hit < HIT_FUNCTION_TAB_BASE + 3) {
            functionTab = hit - HIT_FUNCTION_TAB_BASE;
            updateLayout();
            InvalidateRect(hwnd, nullptr, FALSE);
            notifyAccessibility(EVENT_OBJECT_SELECTION, FORMULA_ACC_TAB_BASE + functionTab);
            notifyAccessibility(EVENT_OBJECT_REORDER, FORMULA_ACC_SELF);
            return;
        }
        if (hit >= HIT_VARIABLE_BASE && hit < HIT_VARIABLE_BASE + 12) {
            size_t index = static_cast<size_t>(hit - HIT_VARIABLE_BASE);
            if (insertFormulaText(variableLabels()[index], false)) selectInspector(variableFromButton(index));
            return;
        }
        if (hit >= HIT_FUNCTION_BASE && hit < HIT_FUNCTION_BASE + 16) {
            size_t index = static_cast<size_t>(hit - HIT_FUNCTION_BASE);
            const wchar_t* label = functionSets()[static_cast<size_t>(functionTab)][index];
            bool function = std::wcscmp(label, L"^") != 0;
            insertFormulaText(label, function);
        }
    }

    bool runSelfTest() {
        std::vector<std::wstring> failures;
        auto check = [&failures](bool condition, const wchar_t* name) {
            if (!condition) failures.emplace_back(name);
        };
        check(FormulaDialogConfig{}.bailout == 100.0, L"default bailout is 100");

        FormulaDialogConfig savedWorking = working;
        FormulaDialogConfig savedApplied = applied;
        FormulaEditorCallbacks savedCallbacks = callbacks;
        InspectorValue savedSelected = selected;
        int savedFunctionTab = functionTab;
        int savedFocused = focusedId;
        int savedScroll = scrollbar.position();
        double savedRange = pickerRange;
        std::wstring savedStatus = status;
        bool savedStatusError = statusError;
        int savedDpi = dpi;
        RECT savedClient{};
        GetClientRect(hwnd, &savedClient);

        HWND formulaBackend = GetDlgItem(hwnd, ID_FORMULA_SERVICE);
        check(formulaBackend != nullptr, L"text backend exists");

        formulaField.setText(L"z");
        formulaField.setSelection(1, 1);
        SendMessageW(formulaBackend, WM_CHAR, L'+', 1);
        SendMessageW(formulaBackend, WM_CHAR, L'c', 1);
        check(formulaField.text() == L"z+c", L"keyboard typing");
        SendMessageW(formulaBackend, WM_KEYDOWN, VK_LEFT, 0);
        check(formulaField.selection().second == 2, L"keyboard navigation");
        SendMessageW(formulaBackend, WM_KEYDOWN, VK_DELETE, 0);
        check(formulaField.text() == L"z+", L"delete key");
        formulaField.setText(L"abc");
        formulaField.setSelection(2, 2);
        SendMessageW(formulaBackend, WM_CHAR, VK_BACK, 1);
        check(formulaField.text() == L"ac", L"backspace key");
        formulaField.setText(L"z");
        formulaField.setSelection(1, 1);
        SendMessageW(formulaBackend, WM_CHAR, 0x03bb, 1);
        check(formulaField.text() == L"z\u03bb", L"unicode input");

        formulaField.setText(L"abc");
        formulaField.setSelection(1, 2);
        SendMessageW(formulaBackend, WM_COPY, 0, 0);
        std::wstring clipboard;
        check(readClipboard(hwnd, clipboard) && clipboard == L"b", L"selection copy");
        SendMessageW(formulaBackend, WM_CUT, 0, 0);
        check(formulaField.text() == L"ac", L"selection cut");
        check(formulaField.undo() && formulaField.text() == L"abc", L"text undo");
        check(writeClipboard(hwnd, L"X"), L"clipboard setup");
        formulaField.setSelection(1, 2);
        SendMessageW(formulaBackend, WM_PASTE, 0, 0);
        check(formulaField.text() == L"aXc", L"selection paste");

        formulaField.setText(L"z*z+c");
        copyCompleteFormula();
        check(readClipboard(hwnd, clipboard) && clipboard == L"z*z+c", L"complete copy");
        check(writeClipboard(hwnd, L"sin(z)+c"), L"complete paste setup");
        pasteCompleteFormula();
        check(formulaField.text() == L"sin(z)+c", L"complete paste");

        std::wstring maximum(79, L'1');
        bailoutField.setText(maximum);
        bailoutField.setSelection(maximum.size(), maximum.size());
        HWND bailoutBackend = GetDlgItem(hwnd, ID_BAILOUT_SERVICE);
        SendMessageW(bailoutBackend, WM_CHAR, L'2', 1);
        check(bailoutField.text().size() == maximum.size(), L"maximum text length");

        working.source = "z";
        formulaField.setText(L"z");
        refreshFormulaAnalysis();
        presetDropdown.setSelectedIndex(0);
        check(presetDropdown.keyDown(VK_SPACE, false, layout.contentClip) && presetDropdown.open(), L"dropdown keyboard open");
        presetDropdown.keyDown(VK_DOWN, false, layout.contentClip);
        presetDropdown.keyDown(VK_RETURN, false, layout.contentClip);
        check(presetDropdown.selectedIndex() == 1 && working.source == "z*z+c", L"dropdown keyboard selection");
        presetDropdown.keyDown(VK_SPACE, false, layout.contentClip);
        presetDropdown.mouseDown({layout.body.left + 2, layout.body.bottom - 2}, layout.contentClip);
        check(!presetDropdown.open(), L"dropdown outside dismissal");
        presetDropdown.close();
        POINT presetCenter{(layout.preset.left + layout.preset.right) / 2, (layout.preset.top + layout.preset.bottom) / 2};
        onLeftButtonDown(presetCenter);
        check(presetDropdown.open() && GetCapture() != hwnd, L"dropdown opens without mouse capture");
        dismissPopups();
        check(!presetDropdown.open(), L"parent-driven popup dismissal");

        working.source = "z+p0+P0";
        formulaField.setText(L"z+p0+P0");
        formulaField.setSelection(0, 0);
        refreshFormulaAnalysis();
        selected = InspectorValue::Z;
        syncInspectorFields();
        updateLayout();

        HDC windowDc = GetDC(hwnd);
        HDC memoryDc = windowDc ? CreateCompatibleDC(windowDc) : nullptr;
        HBITMAP bitmap = windowDc ? CreateCompatibleBitmap(windowDc, std::max(1L, layout.client.right), std::max(1L, layout.client.bottom)) : nullptr;
        HGDIOBJ oldBitmap = memoryDc && bitmap ? SelectObject(memoryDc, bitmap) : nullptr;
        if (memoryDc && bitmap) paintTo(memoryDc);
        check(memoryDc && bitmap, L"offscreen validation surface");

        if (memoryDc && bitmap) {
            HGDIOBJ oldFont = SelectObject(memoryDc, formulaFont);
            SIZE prefix{};
            SIZE tokenSize{};
            GetTextExtentPoint32W(memoryDc, L"z+", 2, &prefix);
            GetTextExtentPoint32W(memoryDc, L"p0", 2, &tokenSize);
            SelectObject(memoryDc, oldFont);
            POINT tokenPoint{layout.formulaField.left + prefix.cx + tokenSize.cx / 2, (layout.formulaField.top + layout.formulaField.bottom) / 2};
            onLeftButtonDown(tokenPoint);
            onLeftButtonUp(tokenPoint);
            auto tokenSelection = formulaField.selection();
            check(selected == InspectorValue::P0 && tokenSelection.first == 2 && tokenSelection.second == 4 && usageCount(InspectorValue::P0) == 2, L"token click selection");
            std::vector<ui::TextRangeStyle> styles = formulaStyles();
            RECT p0Bounds = formulaField.textRangeBounds(tokens[2].first, tokens[2].last, styles[2].paddingBefore);
            POINT tokenRight{p0Bounds.right - 1, (p0Bounds.top + p0Bounds.bottom) / 2};
            formulaField.setSelection(0, 0);
            selected = InspectorValue::Z;
            onLeftButtonDown(tokenRight);
            onLeftButtonUp(tokenRight);
            tokenSelection = formulaField.selection();
            check(selected == InspectorValue::P0 && tokenSelection.first == 2 && tokenSelection.second == 4, L"token rendered-boundary click");
        }

        formulaField.setText(L"z+c");
        working.source = "z+c";
        formulaField.setSelection(0, 1);
        check(insertFormulaText(L"p0", false) && formulaField.text() == L"p0+c", L"variable insertion");
        formulaField.setText(L"z+");
        working.source = "z+";
        formulaField.setSelection(2, 2);
        check(insertFormulaText(L"sin", true) && formulaField.text() == L"z+sin()" && formulaField.selection().first == 6, L"function insertion caret");
        formulaField.setText(L"z+c");
        working.source = "z+c";
        formulaField.setSelection(0, 3);
        check(insertFormulaText(L"sqr", true) && formulaField.text() == L"sqr(z+c)", L"function selection wrapping");

        selected = InspectorValue::P0;
        working.parameters[0] = {};
        pickerRange = 2.0;
        syncInspectorFields();
        updateLayout();
        POINT pickerPoint{layout.picker.left + (layout.picker.right - layout.picker.left) * 3 / 4, layout.picker.top + (layout.picker.bottom - layout.picker.top) / 4};
        updatePickerFromMouse(pickerPoint);
        check(std::abs(working.parameters[0].real() - 1.0) < 0.04 && std::abs(working.parameters[0].imag() - 1.0) < 0.04, L"picker synchronization");
        realField.setText(L"0.25");
        imaginaryField.setText(L"-0.5");
        check(commitInspector(true) && working.parameters[0] == std::complex<double>(0.25, -0.5), L"manual complex synchronization");

        working.pixelParameter = FormulaParameter::C;
        selected = InspectorValue::C;
        syncInspectorFields();
        check(!selectedEditable(working), L"pixel-bound value read-only");
        choosePlane(FormulaParameter::InitialZ);
        check(selectedEditable(working), L"fixed plane value editable");

        bool applyCalled = false;
        bool mandelbrotCalled = false;
        bool closeCalled = false;
        FormulaDialogConfig callbackCandidate;
        callbacks.apply = [&applyCalled, &callbackCandidate](const FormulaDialogConfig& candidate) {
            applyCalled = true;
            callbackCandidate = candidate;
            return true;
        };
        callbacks.useMandelbrot = [&mandelbrotCalled]() { mandelbrotCalled = true; };
        callbacks.close = [&closeCalled]() { closeCalled = true; };

        working.source = "z*z+c";
        formulaField.setText(L"z*z+c");
        bailoutField.setText(L"4");
        selected = InspectorValue::P0;
        realField.setEnabled(true);
        imaginaryField.setEnabled(true);
        realField.setText(L"0.25");
        imaginaryField.setText(L"-0.5");
        refreshFormulaAnalysis();
        applyChanges();
        check(applyCalled && callbackCandidate.source == "z*z+c" && applied.source == "z*z+c", L"transactional apply");
        working.source = "z*z*z+c";
        formulaField.setText(L"z*z*z+c");
        revertChanges();
        check(working.source == "z*z+c" && formulaField.text() == L"z*z+c", L"transactional revert");
        useMandelbrot();
        check(mandelbrotCalled, L"Mandelbrot callback");

        applyCalled = false;
        working.source = "z+*c";
        formulaField.setText(L"z+*c");
        refreshFormulaAnalysis();
        applyChanges();
        auto errorSelection = formulaField.selection();
        check(!applyCalled && !formulaValid && errorSelection.first == 2 && errorSelection.second == 3, L"parse error focus selection");

        ShowWindow(hwnd, SW_SHOW);
        closePanel();
        check(closeCalled && IsWindowVisible(hwnd) == FALSE, L"transactional close");
        ShowWindow(hwnd, SW_SHOW);

        SetWindowPos(hwnd, nullptr, 0, 0, scale(FormulaEditorPanel::DESIGN_WIDTH), scale(640), SWP_NOMOVE | SWP_NOZORDER | SWP_NOACTIVATE);
        updateLayout();
        check(scrollbar.visible() && layout.header.bottom == scale(HEADER_HEIGHT) && layout.footer.bottom == scale(640), L"compact resize custom scrolling");
        SetWindowPos(hwnd, nullptr, 0, 0, scale(FormulaEditorPanel::DESIGN_WIDTH), scale(870), SWP_NOMOVE | SWP_NOZORDER | SWP_NOACTIVATE);
        updateLayout();
        check(!scrollbar.visible(), L"default height no scrolling");

        for (int testDpi : {96, 144, 192}) {
            setDpi(testDpi);
            SetWindowPos(hwnd, nullptr, 0, 0, scale(FormulaEditorPanel::DESIGN_WIDTH), scale(870), SWP_NOMOVE | SWP_NOZORDER | SWP_NOACTIVATE);
            updateLayout();
            int pickerWidth = layout.picker.right - layout.picker.left;
            int pickerHeight = layout.picker.bottom - layout.picker.top;
            check(std::abs(clientWidthDip() - FormulaEditorPanel::DESIGN_WIDTH) <= 1 && layout.header.bottom == scale(HEADER_HEIGHT) && layout.footer.bottom == scale(870) && std::abs(pickerWidth - pickerHeight) <= 1 && !scrollbar.visible(), testDpi == 96 ? L"96 DPI deterministic geometry" : (testDpi == 144 ? L"144 DPI deterministic geometry" : L"192 DPI deterministic geometry"));
        }
        setDpi(savedDpi);
        SetWindowPos(hwnd, nullptr, 0, 0, scale(FormulaEditorPanel::DESIGN_WIDTH), scale(870), SWP_NOMOVE | SWP_NOZORDER | SWP_NOACTIVATE);
        updateLayout();

        struct NativeControlCheck {
            int edits = 0;
            bool okay = true;
        } nativeCheck;
        EnumChildWindows(
            hwnd,
            [](HWND child, LPARAM parameter) -> BOOL {
                auto* result = reinterpret_cast<NativeControlCheck*>(parameter);
                wchar_t className[64]{};
                GetClassNameW(child, className, 64);
                if (_wcsicmp(className, L"COMBOBOX") == 0 || _wcsicmp(className, L"SCROLLBAR") == 0) { result->okay = false; }
                if (_wcsicmp(className, L"EDIT") == 0) {
                    ++result->edits;
                    LONG_PTR style = GetWindowLongPtrW(child, GWL_STYLE);
                    RECT rect{};
                    GetWindowRect(child, &rect);
                    if ((style & (WS_BORDER | WS_VSCROLL | WS_HSCROLL)) != 0 || rect.right - rect.left > 1) { result->okay = false; }
                }
                return TRUE;
            },
            reinterpret_cast<LPARAM>(&nativeCheck));
        check(nativeCheck.okay && nativeCheck.edits == 4, L"no visible native control chrome");

        IAccessible* accessible = nullptr;
        HRESULT accessibleResult = AccessibleObjectFromWindow(hwnd, OBJID_CLIENT, IID_IAccessible, reinterpret_cast<void**>(&accessible));
        check(SUCCEEDED(accessibleResult) && accessible, L"OBJID_CLIENT accessibility provider");
        LONG accessibleChildren = 0;
        if (accessible) { check(SUCCEEDED(accessible->get_accChildCount(&accessibleChildren)) && accessibleChildren >= 50, L"accessible child inventory"); }
        auto findAccessibleChild = [accessible, accessibleChildren](const wchar_t* expectedName) -> LONG {
            if (!accessible) return 0;
            for (LONG childId = 1; childId <= accessibleChildren; ++childId) {
                VARIANT child{};
                child.vt = VT_I4;
                child.lVal = childId;
                BSTR name = nullptr;
                if (SUCCEEDED(accessible->get_accName(child, &name)) && name && std::wcscmp(name, expectedName) == 0) {
                    SysFreeString(name);
                    return childId;
                }
                if (name) SysFreeString(name);
            }
            return 0;
        };
        LONG expressionChild = findAccessibleChild(L"Formula expression");
        LONG presetChild = findAccessibleChild(L"Formula preset");
        LONG bailoutChild = findAccessibleChild(L"Bailout value");
        LONG pickerChild = findAccessibleChild(L"Complex value picker");
        LONG variableChild = findAccessibleChild(L"Variable z");
        std::wstring accessibleFunctionName = std::wstring(L"Function ") + functionSets()[static_cast<size_t>(functionTab)][0];
        LONG functionChild = findAccessibleChild(accessibleFunctionName.c_str());
        LONG scrollbarChild = findAccessibleChild(L"Formula editor scrollbar");
        LONG applyChild = findAccessibleChild(L"Apply and render");
        LONG closeChild = findAccessibleChild(L"Close formula editor");
        LONG presetItemChild = findAccessibleChild(L"Preset Custom");
        check(expressionChild && presetChild && bailoutChild && pickerChild && variableChild && functionChild && scrollbarChild && applyChild && closeChild && presetItemChild, L"accessible key control names");

        auto checkAccessibleRoleState = [&check, accessible](LONG childId, LONG expectedRole, LONG requiredState, const wchar_t* failure) {
            if (!accessible || !childId) {
                check(false, failure);
                return;
            }
            VARIANT child{};
            child.vt = VT_I4;
            child.lVal = childId;
            VARIANT role{};
            VARIANT state{};
            bool okay = SUCCEEDED(accessible->get_accRole(child, &role)) && role.vt == VT_I4 && role.lVal == expectedRole && SUCCEEDED(accessible->get_accState(child, &state)) && state.vt == VT_I4 && (state.lVal & requiredState) == requiredState;
            VariantClear(&role);
            VariantClear(&state);
            check(okay, failure);
        };
        checkAccessibleRoleState(expressionChild, ROLE_SYSTEM_TEXT, STATE_SYSTEM_FOCUSABLE, L"accessible formula role and state");
        checkAccessibleRoleState(presetChild, ROLE_SYSTEM_COMBOBOX, STATE_SYSTEM_FOCUSABLE | STATE_SYSTEM_COLLAPSED, L"accessible preset role and state");
        checkAccessibleRoleState(applyChild, ROLE_SYSTEM_PUSHBUTTON, STATE_SYSTEM_FOCUSABLE, L"accessible action role and state");

        if (accessible) {
            VARIANT selection{};
            bool multipleSelection = false;
            if (SUCCEEDED(accessible->get_accSelection(&selection)) && selection.vt == VT_UNKNOWN && selection.punkVal) {
                IEnumVARIANT* enumerator = nullptr;
                if (SUCCEEDED(selection.punkVal->QueryInterface(IID_IEnumVARIANT, reinterpret_cast<void**>(&enumerator))) && enumerator) {
                    VARIANT values[8]{};
                    ULONG fetched = 0;
                    HRESULT next = enumerator->Next(8, values, &fetched);
                    multipleSelection = (next == S_OK || next == S_FALSE) && fetched >= 3;
                    for (ULONG i = 0; i < fetched; ++i) VariantClear(&values[i]);
                    enumerator->Release();
                }
            }
            VariantClear(&selection);
            check(multipleSelection, L"accessible multiple selection enumeration");

            presetDropdown.keyDown(VK_SPACE, false, layout.contentClip);
            VARIANT child{};
            child.vt = VT_I4;
            child.lVal = presetItemChild;
            LONG left = 0, top = 0, width = 0, height = 0;
            bool locationOkay = SUCCEEDED(accessible->accLocation(&left, &top, &width, &height, child));
            RECT expected = presetDropdown.itemBounds(0, layout.contentClip);
            POINT expectedTop{expected.left, expected.top};
            POINT expectedBottom{expected.right, expected.bottom};
            ClientToScreen(hwnd, &expectedTop);
            ClientToScreen(hwnd, &expectedBottom);
            locationOkay = locationOkay && left == expectedTop.x && top == expectedTop.y && width == expectedBottom.x - expectedTop.x && height == expectedBottom.y - expectedTop.y;
            check(locationOkay, L"accessible DPI-scaled preset item location");
            presetDropdown.close();
        }

        if (accessible && expressionChild) {
            VARIANT self{};
            self.vt = VT_I4;
            self.lVal = CHILDID_SELF;
            VARIANT first{};
            VARIANT next{};
            VARIANT expression{};
            expression.vt = VT_I4;
            expression.lVal = expressionChild;
            check(SUCCEEDED(accessible->accNavigate(NAVDIR_FIRSTCHILD, self, &first)) && first.vt == VT_I4 && first.lVal == expressionChild && SUCCEEDED(accessible->accNavigate(NAVDIR_NEXT, expression, &next)) && next.vt == VT_I4 && next.lVal == variableChild, L"accessible sequential navigation");
            VariantClear(&first);
            VariantClear(&next);
            VARIANT child{};
            child.vt = VT_I4;
            child.lVal = expressionChild;
            check(SUCCEEDED(accessible->accSelect(SELFLAG_TAKEFOCUS, child)) && formulaField.focused(), L"accessible formula focus action");
            BSTR value = nullptr;
            check(SUCCEEDED(accessible->get_accValue(child, &value)) && value, L"accessible formula value");
            if (value) SysFreeString(value);
        }

        SetWindowPos(hwnd, nullptr, 0, 0, scale(FormulaEditorPanel::DESIGN_WIDTH), scale(640), SWP_NOMOVE | SWP_NOZORDER | SWP_NOACTIVATE);
        updateLayout();
        if (accessible && scrollbarChild) {
            VARIANT child{};
            child.vt = VT_I4;
            child.lVal = scrollbarChild;
            VARIANT state{};
            check(SUCCEEDED(accessible->get_accState(child, &state)) && state.vt == VT_I4 && (state.lVal & STATE_SYSTEM_INVISIBLE) == 0, L"accessible visible custom scrollbar");
            VariantClear(&state);
        }
        SetWindowPos(hwnd, nullptr, 0, 0, scale(FormulaEditorPanel::DESIGN_WIDTH), scale(870), SWP_NOMOVE | SWP_NOZORDER | SWP_NOACTIVATE);
        updateLayout();

        DWORD detachedGdiBefore = GetGuiResources(GetCurrentProcess(), GR_GDIOBJECTS);
        DWORD detachedUserBefore = GetGuiResources(GetCurrentProcess(), GR_USEROBJECTS);
        bool detachOkay = true;
        for (int i = 0; i < 32; ++i) {
            FormulaEditorAccessibilityProvider* temporary = createFormulaEditorAccessibility(hwnd);
            IAccessible* disconnected = acquireFormulaEditorAccessibility(temporary);
            detachFormulaEditorAccessibility(temporary);
            LONG count = -1;
            if (!disconnected || disconnected->get_accChildCount(&count) != CO_E_OBJNOTCONNECTED) { detachOkay = false; }
            if (disconnected) disconnected->Release();
            destroyFormulaEditorAccessibility(temporary);
        }
        DWORD detachedGdiAfter = GetGuiResources(GetCurrentProcess(), GR_GDIOBJECTS);
        DWORD detachedUserAfter = GetGuiResources(GetCurrentProcess(), GR_USEROBJECTS);
        check(detachOkay && (detachedGdiBefore == 0 || detachedGdiAfter <= detachedGdiBefore + 1) && (detachedUserBefore == 0 || detachedUserAfter <= detachedUserBefore + 1), L"accessible provider detach resource stability");
        if (accessible) accessible->Release();

        if (memoryDc && bitmap) {
            DWORD before = GetGuiResources(GetCurrentProcess(), GR_GDIOBJECTS);
            for (int i = 0; i < 96; ++i) paintTo(memoryDc);
            DWORD after = GetGuiResources(GetCurrentProcess(), GR_GDIOBJECTS);
            check(before == 0 || after <= before + 1, L"GDI resource stability");
        }

        if (oldBitmap && memoryDc) SelectObject(memoryDc, oldBitmap);
        if (bitmap) DeleteObject(bitmap);
        if (memoryDc) DeleteDC(memoryDc);
        if (windowDc) ReleaseDC(hwnd, windowDc);

        callbacks = std::move(savedCallbacks);
        working = savedWorking;
        applied = savedApplied;
        selected = savedSelected;
        functionTab = savedFunctionTab;
        pickerRange = savedRange;
        status = std::move(savedStatus);
        statusError = savedStatusError;
        setDpi(savedDpi);
        SetWindowPos(hwnd, nullptr, 0, 0, savedClient.right - savedClient.left, savedClient.bottom - savedClient.top, SWP_NOMOVE | SWP_NOZORDER | SWP_NOACTIVATE);
        syncAllControls();
        scrollbar.setPosition(savedScroll);
        updateLayout();
        setFocused(savedFocused);
        ShowWindow(hwnd, SW_SHOW);

        std::wstring resultPath = environmentValue(L"MANDEL_FORMULA_EDITOR_SELFTEST_RESULT");
        if (!writeSelfTestResult(resultPath, failures)) { writeSelfTestResult(L"build\\formula_editor_selftest.txt", failures); }
        if (failures.empty()) {
            OutputDebugStringW(L"Formula editor self-test: PASS\n");
        } else {
            OutputDebugStringW(L"Formula editor self-test: FAIL\n");
            for (const std::wstring& failure : failures) {
                OutputDebugStringW(failure.c_str());
                OutputDebugStringW(L"\n");
            }
        }
        return failures.empty();
    }

    void runSelfTestIfRequested() {
        if (selfTestRan || environmentValue(L"MANDEL_FORMULA_EDITOR_SELFTEST").empty()) { return; }
        selfTestRan = true;
        runSelfTest();
        if (!environmentValue(L"MANDEL_FORMULA_EDITOR_SELFTEST_EXIT").empty()) { PostMessageW(owner, WM_CLOSE, 0, 0); }
    }

    void show(const FormulaDialogConfig& config) {
        working = config;
        applied = config;
        pickerRange = 2.0;
        scrollbar.setPosition(0);
        status.clear();
        statusError = false;
        syncing = true;
        std::wstring source;
        if (!widenUtf8(config.source, source)) source.clear();
        formulaField.setText(source);
        formulaField.setSelection(source.size(), source.size());
        syncing = false;
        refreshFormulaAnalysis();
        selected = initialInspector();
        syncAllControls();
        ShowWindow(hwnd, SW_SHOWNOACTIVATE);
        InvalidateRect(hwnd, nullptr, TRUE);
        notifyAccessibility(EVENT_OBJECT_SHOW, FORMULA_ACC_SELF);
        runSelfTestIfRequested();
    }

    void setConfig(const FormulaDialogConfig& config) {
        working = config;
        applied = config;
        status.clear();
        statusError = false;
        syncAllControls();
    }

    void dismissPopups() {
        if (!presetDropdown.open()) return;
        bool wasOpen = presetDropdown.open();
        presetDropdown.close();
        notifyDropdownTransition(wasOpen);
        InvalidateRect(hwnd, nullptr, FALSE);
    }

    void hide() {
        if (!hwnd) return;
        bool dropdownWasOpen = presetDropdown.open();
        presetDropdown.close();
        notifyDropdownTransition(dropdownWasOpen);
        if (GetCapture() == hwnd) ReleaseCapture();
        draggingText = nullptr;
        draggingPicker = false;
        pressedHit = HIT_NONE;
        ShowWindow(hwnd, SW_HIDE);
        notifyAccessibility(EVENT_OBJECT_HIDE, FORMULA_ACC_SELF);
    }

    std::vector<ui::TextRangeStyle> formulaStyles() const {
        std::vector<ui::TextRangeStyle> styles;
        std::wstring source = formulaField.text();
        for (size_t index = 0; index < tokens.size(); ++index) {
            const Token& token = tokens[index];
            ui::TextRangeStyle style;
            style.first = token.first;
            style.last = token.last;
            switch (token.kind) {
            case TokenKind::Identifier:
                if (token.selectable)
                    style.text = TOKEN_VARIABLE;
                else if (token.function)
                    style.text = TOKEN_FUNCTION;
                else
                    style.text = TOKEN_UNKNOWN;
                break;
            case TokenKind::Number: style.text = TOKEN_NUMBER; break;
            case TokenKind::Operator:
                style.text = TOKEN_OPERATOR;
                if (token.first < source.size() && source[token.first] != L'=') {
                    if (token.first == 0 || !std::iswspace(source[token.first - 1])) style.paddingBefore = scale(3);
                    if (token.last >= source.size() || !std::iswspace(source[token.last])) style.paddingAfter = scale(3);
                }
                break;
            case TokenKind::Punctuation: style.text = TOKEN_OPERATOR; break;
            case TokenKind::Space: style.text = CLR_TEXT; break;
            }
            if (token.selectable) {
                style.paddingBefore += scale(3);
                style.paddingAfter += scale(3);
            }
            if (token.selectable && token.inspector == selected) {
                style.text = RGB(255, 255, 255);
                style.background = RGB(39, 56, 88);
                style.border = CLR_ACCENT;
            } else if (token.selectable && static_cast<int>(index) == hoveredFormulaToken) {
                style.background = RGB(35, 43, 56);
                style.border = RGB(70, 80, 100);
            }
            styles.push_back(style);
        }
        return styles;
    }

    void drawSectionTitle(HDC dc, RECT rect, const std::wstring& title, const std::wstring& hint) {
        RECT titleRect = rect;
        int oldExtra = SetTextCharacterExtra(dc, scale(1));
        drawText(dc, titleRect, title, CLR_TEXT_DIM, sectionFont, DT_LEFT | DT_VCENTER | DT_SINGLELINE);
        SetTextCharacterExtra(dc, oldExtra);
        if (!hint.empty()) { drawText(dc, titleRect, hint, RGB(98, 108, 128), tinyFont, DT_RIGHT | DT_VCENTER | DT_SINGLELINE); }
    }

    void drawHeader(HDC dc) {
        fillRect(dc, layout.header, HEADER_BG);
        drawSeparator(dc, layout.header.left, layout.header.bottom - 1, layout.header.right, layout.header.bottom - 1, SOFT_BORDER);
        RECT title = pixelRect(18, 11, 300, 34);
        drawText(dc, title, L"Formula editor", CLR_TEXT, headerFont, DT_LEFT | DT_VCENTER | DT_SINGLELINE);
        RECT subtitle = pixelRect(18, 33, 360, 51);
        drawText(dc, subtitle, L"Build the orbit, then bind values visually", CLR_TEXT_DIM, resources.small(), DT_LEFT | DT_VCENTER | DT_SINGLELINE);

        int widthDip = clientWidthDip();
        RECT badge = pixelRect(widthDip - 97, 20, widthDip - 56, 43);
        fillRound(dc, badge, RGB(24, 47, 43), RGB(45, 99, 80), scale(11));
        drawText(dc, badge, L"DOCKED", CLR_GREEN, tinyFont, DT_CENTER | DT_VCENTER | DT_SINGLELINE);
    }

    void drawFormulaSection(HDC dc) {
        int left = layout.contentLeftDip;
        int right = layout.contentRightDip;
        drawSectionTitle(dc, contentPixelRect(left, 0, right, 21), L"ORBIT EXPRESSION", L"Click a variable to edit its value");
        RECT glow = layout.formulaCard;
        InflateRect(&glow, scale(3), scale(3));
        fillRound(dc, glow, DOCK_BG, RGB(27, 38, 56), scale(11));
        drawCardWithColor(dc, layout.formulaCard, EDIT_BG, RGB(70, 80, 100), scale(9));

        RECT prefix = contentPixelRect(left + 14, 27, left + 60, 79);
        drawText(dc, prefix, L"z' =", RGB(119, 131, 153), formulaFont, DT_LEFT | DT_VCENTER | DT_SINGLELINE);
        ui::TextFieldStyle style;
        style.drawChrome = false;
        style.fill = EDIT_BG;
        style.text = CLR_TEXT;
        style.horizontalPadding = 0;
        style.selection = RGB(56, 91, 150);
        formulaField.draw(dc, formulaFont, style, formulaStyles());

        drawSeparator(dc, layout.formulaMeta.left, layout.formulaMeta.top, layout.formulaMeta.right, layout.formulaMeta.top, RGB(37, 42, 53));
        int dotX = contentPixelRect(left + 14, 0, left + 14, 0).left;
        int centerY = (layout.formulaMeta.top + layout.formulaMeta.bottom) / 2;
        drawCircle(dc, dotX, centerY, scale(2), formulaValid ? CLR_GREEN : DANGER, formulaValid ? CLR_GREEN : DANGER);

        std::wstring meta;
        if (formulaValid) {
            wchar_t buffer[96]{};
            swprintf_s(buffer, L"Valid formula \u00b7 %zu operations", formulaProgram.instructionCount());
            meta = buffer;
        } else {
            std::wstring error;
            if (!widenUtf8(formulaError.message, error)) error = L"Invalid formula";
            meta = error;
        }
        RECT metaText = layout.formulaMeta;
        metaText.left += scale(22);
        metaText.right -= scale(128);
        drawText(dc, metaText, meta, formulaValid ? CLR_TEXT_DIM : DANGER, tinyFont, DT_LEFT | DT_VCENTER | DT_SINGLELINE | DT_END_ELLIPSIS);

        std::wstring badge;
        if (formulaValid && formulaProgram.fastPath() == formula::ExpressionProgram::FastPath::IntegerPowerPlusC) {
            wchar_t buffer[48]{};
            swprintf_s(buffer, L"FAST z^%d+c", formulaProgram.fastIntegerPower());
            badge = buffer;
        } else if (formulaValid && formulaProgram.avx2Compatible()) {
            badge = L"AVX2 READY";
        } else if (formulaValid) {
            badge = L"SCALAR";
        }
        RECT badgeRect = layout.formulaMeta;
        badgeRect.left = badgeRect.right - scale(122);
        badgeRect.right -= scale(10);
        drawText(dc, badgeRect, badge, TOKEN_VARIABLE, tinyFont, DT_RIGHT | DT_VCENTER | DT_SINGLELINE);
    }

    static void drawCardWithColor(HDC dc, RECT rect, COLORREF fill, COLORREF border, int radius) { fillRound(dc, rect, fill, border, radius); }

    void drawButton(HDC dc, const ButtonSpec& button) {
        bool hovered = hoverHit == button.hit;
        bool pressed = pressedHit == button.hit;
        COLORREF fill = CLR_CARD;
        COLORREF border = CLR_BORDER;
        COLORREF text = CLR_TEXT;
        HFONT font = resources.regular();
        int radius = scale(6);

        switch (button.kind) {
        case ButtonKind::Primary:
            fill = pressed ? RGB(75, 119, 198) : (hovered ? CLR_ACCENT_HI : CLR_ACCENT);
            border = fill;
            text = RGB(255, 255, 255);
            font = actionFont;
            radius = scale(8);
            break;
        case ButtonKind::Chip:
            font = chipFont;
            if (button.active) {
                fill = RGB(50, 64, 95);
                border = CLR_ACCENT;
                text = RGB(255, 255, 255);
            } else if (hovered) {
                fill = CLR_CARD_HOV;
                border = RGB(89, 98, 118);
            }
            break;
        case ButtonKind::Function:
            font = functionFont;
            fill = hovered ? CLR_CARD_HOV : RGB(36, 40, 51);
            border = hovered ? CLR_ACCENT : RGB(52, 58, 72);
            text = RGB(191, 200, 215);
            break;
        case ButtonKind::Segment:
            fill = button.active ? CLR_ACCENT : EDIT_BG;
            border = EDIT_BG;
            text = button.active ? RGB(255, 255, 255) : CLR_TEXT_DIM;
            break;
        case ButtonKind::Tab:
            fill = button.active ? CLR_CARD : VALUE_CARD;
            border = fill;
            text = button.active ? CLR_TEXT : CLR_TEXT_DIM;
            font = tinyFont;
            break;
        case ButtonKind::Icon:
            fill = hovered ? CLR_CARD_HOV : CLR_CARD;
            border = CLR_BORDER;
            text = hovered ? CLR_TEXT : CLR_TEXT_DIM;
            font = resources.regular();
            break;
        case ButtonKind::Range:
            font = functionFont;
            if (hovered) fill = CLR_CARD_HOV;
            break;
        case ButtonKind::Normal:
        default:
            if (hovered) fill = CLR_CARD_HOV;
            break;
        }
        if (button.hit == HIT_COPY || button.hit == HIT_PASTE) font = tinyFont;
        if (pressed && button.kind != ButtonKind::Primary) fill = RGB(31, 35, 44);
        if (!button.enabled) {
            fill = VALUE_CARD;
            border = SOFT_BORDER;
            text = RGB(92, 99, 113);
        }
        fillRound(dc, button.rect, fill, border, radius);
        RECT textRect = button.rect;
        if (pressed) OffsetRect(&textRect, 0, 1);
        drawText(dc, textRect, button.label, text, font, DT_CENTER | DT_VCENTER | DT_SINGLELINE | DT_END_ELLIPSIS | DT_NOPREFIX);
        if (button.kind == ButtonKind::Chip) {
            int variable = button.hit - HIT_VARIABLE_BASE;
            if (variable >= 0 && variable < 12 && usageCount(variableFromButton(static_cast<size_t>(variable))) > 0) {
                int dotX = button.rect.right - scale(6);
                int dotY = button.rect.top + scale(6);
                drawCircle(dc, dotX, dotY, scale(2), CLR_GREEN, CLR_GREEN);
            }
        }
        if (focusedId == button.hit) {
            RECT ring = button.rect;
            InflateRect(&ring, scale(1), scale(1));
            ui::drawFocusRing(dc, ring, radius + scale(1));
        }
    }

    void drawLeftColumn(HDC dc) {
        drawCardWithColor(dc, layout.variableCard, VALUE_CARD, SOFT_BORDER, scale(9));
        drawSectionTitle(dc, contentPixelRect(layout.contentLeftDip + 11, 139, layout.contentLeftDip + unscale(layout.variableCard.right - layout.variableCard.left) - 11, 156), L"INSERT VARIABLE", L"at cursor");

        drawCardWithColor(dc, layout.functionCard, VALUE_CARD, SOFT_BORDER, scale(9));
        drawSectionTitle(dc, contentPixelRect(layout.contentLeftDip + 11, 288, layout.contentLeftDip + unscale(layout.functionCard.right - layout.functionCard.left) - 11, 306), L"INSERT FUNCTION", L"wrap selection");

        drawCardWithColor(dc, layout.planeCard, VALUE_CARD, SOFT_BORDER, scale(9));
        drawSectionTitle(dc, contentPixelRect(layout.contentLeftDip + 11, 499, layout.contentLeftDip + unscale(layout.planeCard.right - layout.planeCard.left) - 11, 517), L"PIXEL PLANE", L"what each pixel controls");
        RECT segmentBackground = contentPixelRect(layout.contentLeftDip + 9, 521, layout.contentLeftDip + unscale(layout.planeCard.right - layout.planeCard.left) - 9, 556);
        fillRound(dc, segmentBackground, EDIT_BG, EDIT_BG, scale(7));
        RECT bailoutLabel = contentPixelRect(layout.contentLeftDip + 11, 566, layout.contentLeftDip + 79, 593);
        drawText(dc, bailoutLabel, L"Bailout |z|", CLR_TEXT_DIM, tinyFont, DT_LEFT | DT_VCENTER | DT_SINGLELINE);
        ui::TextFieldStyle numericStyle;
        numericStyle.radius = scale(6);
        numericStyle.horizontalPadding = scale(8);
        bailoutField.draw(dc, controlFont, numericStyle);

        drawCardWithColor(dc, layout.capabilityCard, RGB(25, 45, 41), RGB(43, 83, 68), scale(7));
        std::wstring capability;
        if (!formulaValid) {
            capability = L"Unavailable until the formula is valid.";
        } else {
            capability = L"expression bytecode";
            if (formulaProgram.avx2Compatible()) capability += L", AVX2 batch";
            if (formulaProgram.derivativeCompatible()) capability += L", derivatives";
            capability += L".";
        }
        RECT capabilityText = layout.capabilityCard;
        capabilityText.left += scale(9);
        capabilityText.right -= scale(9);
        capabilityText.top += scale(6);
        capabilityText.bottom -= scale(5);
        if (formulaValid) {
            RECT available = capabilityText;
            available.right = available.left + scale(54);
            drawText(dc, available, L"Available:", CLR_GREEN, sectionFont, DT_LEFT | DT_TOP | DT_SINGLELINE);
            capabilityText.left = available.right;
        }
        drawText(dc, capabilityText, capability, formulaValid ? RGB(167, 183, 173) : CLR_TEXT_DIM, tinyFont, DT_LEFT | DT_TOP | DT_WORDBREAK | DT_END_ELLIPSIS);

        drawCardWithColor(dc, layout.presetCard, VALUE_CARD, SOFT_BORDER, scale(8));
        RECT presetLabel = contentPixelRect(layout.contentLeftDip + 11, 675, layout.contentLeftDip + 50, 704);
        drawText(dc, presetLabel, L"Preset", CLR_TEXT_DIM, tinyFont, DT_LEFT | DT_VCENTER | DT_SINGLELINE);
    }

    void drawPicker(HDC dc) {
        RECT rect = layout.picker;
        drawCardWithColor(dc, rect, RGB(17, 21, 28), RGB(83, 96, 118), scale(7));
        HPEN grid = CreatePen(PS_SOLID, 1, RGB(42, 47, 58));
        HPEN axis = CreatePen(PS_SOLID, 1, RGB(89, 98, 118));
        if (grid && axis) {
            HGDIOBJ old = SelectObject(dc, grid);
            for (int i = 1; i < 4; ++i) {
                int x = rect.left + MulDiv(rect.right - rect.left, i, 4);
                int y = rect.top + MulDiv(rect.bottom - rect.top, i, 4);
                MoveToEx(dc, x, rect.top, nullptr);
                LineTo(dc, x, rect.bottom);
                MoveToEx(dc, rect.left, y, nullptr);
                LineTo(dc, rect.right, y);
            }
            SelectObject(dc, axis);
            int centerX = (rect.left + rect.right) / 2;
            int centerY = (rect.top + rect.bottom) / 2;
            MoveToEx(dc, centerX, rect.top, nullptr);
            LineTo(dc, centerX, rect.bottom);
            MoveToEx(dc, rect.left, centerY, nullptr);
            LineTo(dc, rect.right, centerY);
            SelectObject(dc, old);
        }
        if (grid) DeleteObject(grid);
        if (axis) DeleteObject(axis);

        int centerX = (rect.left + rect.right) / 2;
        int centerY = (rect.top + rect.bottom) / 2;
        drawCircle(dc, centerX, centerY, scale(2), RGB(107, 115, 130), RGB(107, 115, 130));

        const std::complex<double>* value = selectedValue(working);
        std::complex<double> point = value ? *value : std::complex<double>{};
        if (finiteComplex(point)) {
            double normalizedX = point.real() / pickerRange * 0.5 + 0.5;
            double normalizedY = 0.5 - point.imag() / pickerRange * 0.5;
            normalizedX = std::clamp(normalizedX, 0.0, 1.0);
            normalizedY = std::clamp(normalizedY, 0.0, 1.0);
            int x = rect.left + static_cast<int>(normalizedX * (rect.right - rect.left - 1));
            int y = rect.top + static_cast<int>(normalizedY * (rect.bottom - rect.top - 1));
            drawCircle(dc, x, y, scale(10), RGB(40, 59, 92), RGB(40, 59, 92));
            drawCircle(dc, x, y, scale(6), selectedEditable(working) ? CLR_ACCENT : CLR_TEXT_DIM, RGB(255, 255, 255));
        }

        if (!selectedEditable(working)) {
            RECT overlay{rect.left + scale(48), centerY - scale(16), rect.right - scale(48), centerY + scale(16)};
            fillRound(dc, overlay, RGB(24, 27, 35), CLR_BORDER, scale(6));
            drawText(dc, overlay, L"Read-only", CLR_TEXT_DIM, resources.small(), DT_CENTER | DT_VCENTER | DT_SINGLELINE);
        }
        if (focusedId == HIT_PICKER) {
            RECT ring = rect;
            InflateRect(&ring, scale(1), scale(1));
            ui::drawFocusRing(dc, ring, scale(8));
        }
    }

    std::wstring complexPreviewText() const {
        std::complex<double> value{};
        const std::complex<double>* pointer = selectedValue(working);
        if (pointer) value = *pointer;
        std::wstring result = inspectorName(selected);
        result += L" = ";
        result += formatNumber(value.real(), 7);
        result += value.imag() < 0.0 ? L" - " : L" + ";
        result += formatNumber(std::abs(value.imag()), 7);
        result += L"i";
        return result;
    }

    void drawValueColumn(HDC dc) {
        drawCardWithColor(dc, layout.valueCard, VALUE_CARD, SOFT_BORDER, scale(9));
        int rightLeft = unscale(layout.valueCard.left);
        int contentRight = layout.contentRightDip;
        drawSectionTitle(dc, contentPixelRect(rightLeft + 12, 139, contentRight - 12, 156), L"SELECTED VALUE", usageLabel());

        RECT symbol = contentPixelRect(rightLeft + 13, 170, rightLeft + 47, 204);
        fillRound(dc, symbol, RGB(43, 64, 100), CLR_ACCENT, scale(8));
        drawText(dc, symbol, inspectorName(selected), RGB(255, 255, 255), formulaFont, DT_CENTER | DT_VCENTER | DT_SINGLELINE);

        RECT title = contentPixelRect(rightLeft + 56, 164, contentRight - 68, 183);
        drawText(dc, title, inspectorTitle(), CLR_TEXT, chipFont, DT_LEFT | DT_VCENTER | DT_SINGLELINE | DT_END_ELLIPSIS);
        RECT description = contentPixelRect(rightLeft + 56, 181, contentRight - 12, 202);
        drawText(dc, description, inspectorDescription(), CLR_TEXT_DIM, tinyFont, DT_LEFT | DT_VCENTER | DT_SINGLELINE | DT_END_ELLIPSIS);

        RECT typeBadge = contentPixelRect(contentRight - 68, 166, contentRight - 13, 187);
        fillRound(dc, typeBadge, RGB(35, 48, 73), RGB(35, 48, 73), scale(8));
        drawText(dc, typeBadge, selectedEditable(working) ? L"COMPLEX" : L"READ ONLY", RGB(169, 197, 251), tinyFont, DT_CENTER | DT_VCENTER | DT_SINGLELINE);

        drawPicker(dc);
        RECT yAxis = layout.picker;
        yAxis.right = yAxis.left - scale(4);
        yAxis.left -= scale(22);
        drawText(dc, yAxis, L"Im", CLR_TEXT_DIM, tinyFont, DT_CENTER | DT_VCENTER | DT_SINGLELINE);
        RECT xAxis{layout.picker.left, layout.picker.bottom, layout.picker.right, layout.picker.bottom + scale(14)};
        drawText(dc, xAxis, L"Re", CLR_TEXT_DIM, tinyFont, DT_CENTER | DT_VCENTER | DT_SINGLELINE);

        wchar_t range[96]{};
        swprintf_s(range, L"Range: -%.3g \u2026 %.3g", pickerRange, pickerRange);
        RECT rangeText = contentPixelRect(rightLeft + 13, 486, contentRight - 95, 509);
        drawText(dc, rangeText, range, CLR_TEXT_DIM, tinyFont, DT_LEFT | DT_VCENTER | DT_SINGLELINE | DT_END_ELLIPSIS);

        int valueInnerLeft = rightLeft + 13;
        int valueInnerRight = contentRight - 13;
        int gap = 7;
        int width = (valueInnerRight - valueInnerLeft - gap) / 2;
        drawText(dc, contentPixelRect(valueInnerLeft, 513, valueInnerLeft + width, 527), L"Real", CLR_TEXT_DIM, tinyFont, DT_LEFT | DT_VCENTER | DT_SINGLELINE);
        drawText(dc, contentPixelRect(valueInnerLeft + width + gap, 513, valueInnerRight, 527), L"Imaginary", CLR_TEXT_DIM, tinyFont, DT_LEFT | DT_VCENTER | DT_SINGLELINE);

        ui::TextFieldStyle numericStyle;
        numericStyle.radius = scale(6);
        numericStyle.horizontalPadding = scale(8);
        realField.draw(dc, controlFont, numericStyle);
        imaginaryField.draw(dc, controlFont, numericStyle);

        fillRound(dc, layout.complexPreview, RGB(21, 24, 32), RGB(21, 24, 32), scale(6));
        drawText(dc, layout.complexPreview, complexPreviewText(), RGB(223, 232, 247), functionFont, DT_CENTER | DT_VCENTER | DT_SINGLELINE | DT_END_ELLIPSIS);
    }

    void drawContent(HDC dc) {
        int saved = SaveDC(dc);
        IntersectClipRect(dc, layout.contentClip.left, layout.contentClip.top, layout.contentClip.right, layout.contentClip.bottom);
        drawFormulaSection(dc);
        drawLeftColumn(dc);
        drawValueColumn(dc);
        for (const ButtonSpec& button : buttons) {
            if (button.content && button.hit != HIT_PICKER) drawButton(dc, button);
        }
        ui::DropdownStyle dropdownStyle;
        dropdownStyle.radius = scale(6);
        dropdownStyle.horizontalPadding = scale(8);
        presetDropdown.draw(dc, resources.small(), dropdownStyle);
        presetDropdown.drawPopup(dc, resources.small(), dropdownStyle, layout.contentClip);
        RestoreDC(dc, saved);
        scrollbar.draw(dc, {});
    }

    void drawFooter(HDC dc) {
        fillRect(dc, layout.footer, HEADER_BG);
        drawSeparator(dc, layout.footer.left, layout.footer.top, layout.footer.right, layout.footer.top, SOFT_BORDER);
        int footerTopDip = clientHeightDip() - FOOTER_HEIGHT;
        int buttonBoundary = clientWidthDip() - 385;
        RECT headline = pixelRect(16, footerTopDip + 12, std::max(90, buttonBoundary), footerTopDip + 28);
        std::wstring headlineText;
        COLORREF headlineColor = CLR_GREEN;
        if (statusError || !formulaValid) {
            headlineText = L"Invalid";
            headlineColor = DANGER;
        } else if (dirty()) {
            headlineText = L"Staged";
            headlineColor = TOKEN_NUMBER;
        } else {
            headlineText = L"Ready";
        }
        drawText(dc, headline, headlineText, headlineColor, tinyFont, DT_LEFT | DT_VCENTER | DT_SINGLELINE);

        RECT detail = pixelRect(16, footerTopDip + 28, std::max(90, buttonBoundary), footerTopDip + 51);
        std::wstring detailText = status.empty() ? (dirty() ? L"Changes are staged until Apply" : L"Configuration matches the last Apply") : status;
        drawText(dc, detail, detailText, statusError ? DANGER : CLR_TEXT_DIM, tinyFont, DT_LEFT | DT_TOP | DT_SINGLELINE | DT_END_ELLIPSIS);

        for (const ButtonSpec& button : buttons) {
            if (!button.content && button.hit != HIT_CLOSE) { drawButton(dc, button); }
        }
    }

    void paintTo(HDC dc) {
        fillRect(dc, layout.client, DOCK_BG);
        drawHeader(dc);
        drawContent(dc);
        for (const ButtonSpec& button : buttons) {
            if (button.hit == HIT_CLOSE) drawButton(dc, button);
        }
        drawFooter(dc);
    }

    void paint() {
        PAINTSTRUCT paintStruct{};
        HDC target = BeginPaint(hwnd, &paintStruct);
        RECT client{};
        GetClientRect(hwnd, &client);
        int width = std::max(1L, client.right - client.left);
        int height = std::max(1L, client.bottom - client.top);
        HDC memory = backBuffer.begin(target, width, height);
        if (memory) {
            paintTo(memory);
            backBuffer.present(target, client);
        } else {
            paintTo(target);
        }
        EndPaint(hwnd, &paintStruct);
    }

    bool pointInContent(POINT point) const { return containsPoint(layout.contentClip, point); }

    void beginTextMouse(ui::TextField& field, int focus, POINT point, bool selectToken) {
        setFocused(focus);
        size_t position = field.indexAtPoint(point);
        int tokenIndex = &field == &formulaField ? tokenIndexAtPoint(point) : -1;
        bool extend = (GetKeyState(VK_SHIFT) & 0x8000) != 0;
        if (!field.mouseDown(point, extend)) return;
        draggingText = &field;
        SetCapture(hwnd);
        if (&field == &formulaField) {
            const Token* token = tokenIndex >= 0 ? &tokens[static_cast<size_t>(tokenIndex)] : tokenAt(position);
            if (token && token->selectable && selectInspector(token->inspector, false)) {
                if (selectToken) formulaField.setSelection(token->first, token->last);
            }
        }
        InvalidateRect(hwnd, nullptr, FALSE);
    }

    void onLeftButtonDown(POINT point) {
        bool dropdownWasOpen = presetDropdown.open();
        auto dropdownResult = presetDropdown.mouseDown(point, layout.contentClip);
        notifyDropdownTransition(dropdownWasOpen);
        if (dropdownResult != ui::Dropdown::MouseResult::Ignored) {
            setFocused(FOCUS_PRESET);
            if (dropdownResult == ui::Dropdown::MouseResult::SelectionChanged) {
                notifyAccessibility(EVENT_OBJECT_SELECTION, FORMULA_ACC_PRESET);
                notifyAccessibility(EVENT_OBJECT_VALUECHANGE, FORMULA_ACC_PRESET);
            }
            InvalidateRect(hwnd, nullptr, FALSE);
            return;
        }

        if (scrollbar.mouseDown(point)) {
            setFocused(FOCUS_SCROLLBAR);
            SetCapture(hwnd);
            InvalidateRect(hwnd, nullptr, FALSE);
            return;
        }

        if (pointInContent(point)) {
            if (containsPoint(layout.formulaField, point)) {
                beginTextMouse(formulaField, FOCUS_FORMULA, point, true);
                return;
            }
            if (containsPoint(layout.bailout, point)) {
                beginTextMouse(bailoutField, FOCUS_BAILOUT, point, false);
                return;
            }
            if (containsPoint(layout.realField, point) && realField.enabled()) {
                beginTextMouse(realField, FOCUS_REAL, point, false);
                return;
            }
            if (containsPoint(layout.imaginaryField, point) && imaginaryField.enabled()) {
                beginTextMouse(imaginaryField, FOCUS_IMAGINARY, point, false);
                return;
            }
            if (containsPoint(layout.picker, point)) {
                setFocused(HIT_PICKER);
                if (!selectedEditable(working)) {
                    setStatus(inspectorDescription(), true);
                    return;
                }
                draggingPicker = true;
                SetCapture(hwnd);
                updatePickerFromMouse(point);
                return;
            }
        }

        int hit = hitRouter.hit(point.x, point.y);
        if (hit != HIT_NONE) {
            setFocused(hit);
            pressedHit = hit;
            hitRouter.press(point.x, point.y);
            SetCapture(hwnd);
            InvalidateRect(hwnd, nullptr, FALSE);
        } else {
            SetFocus(hwnd);
        }
    }

    void onLeftButtonDoubleClick(POINT point) {
        if (!pointInContent(point)) return;
        if (containsPoint(layout.formulaField, point)) {
            beginTextMouse(formulaField, FOCUS_FORMULA, point, true);
            formulaField.selectWordAt(formulaField.indexAtPoint(point));
        } else if (containsPoint(layout.bailout, point)) {
            setFocused(FOCUS_BAILOUT);
            bailoutField.selectAll();
        } else if (containsPoint(layout.realField, point) && realField.enabled()) {
            setFocused(FOCUS_REAL);
            realField.selectAll();
        } else if (containsPoint(layout.imaginaryField, point) && imaginaryField.enabled()) {
            setFocused(FOCUS_IMAGINARY);
            imaginaryField.selectAll();
        }
        InvalidateRect(hwnd, nullptr, FALSE);
    }

    void onMouseMove(POINT point) {
        if (draggingText) {
            if (draggingText->mouseMove(point)) InvalidateRect(hwnd, nullptr, FALSE);
            SetCursor(LoadCursorW(nullptr, IDC_IBEAM));
            return;
        }
        if (draggingPicker) {
            updatePickerFromMouse(point);
            SetCursor(LoadCursorW(nullptr, IDC_CROSS));
            return;
        }
        if (scrollbar.dragging()) {
            if (scrollbar.mouseMove(point)) InvalidateRect(hwnd, nullptr, FALSE);
            return;
        }

        bool changed = presetDropdown.mouseMove(point, layout.contentClip);
        int newHoveredToken = -1;
        if (containsPoint(layout.formulaField, point)) { newHoveredToken = tokenIndexAtPoint(point); }
        if (newHoveredToken != hoveredFormulaToken) {
            hoveredFormulaToken = newHoveredToken;
            changed = true;
        }
        changed = scrollbar.mouseMove(point) || changed;
        if (hitRouter.move(point.x, point.y)) {
            hoverHit = hitRouter.hovered();
            changed = true;
        }
        if (changed) InvalidateRect(hwnd, nullptr, FALSE);

        if (!trackingMouse) {
            TRACKMOUSEEVENT event{};
            event.cbSize = sizeof(event);
            event.dwFlags = TME_LEAVE;
            event.hwndTrack = hwnd;
            trackingMouse = TrackMouseEvent(&event) != FALSE;
        }

        if (containsPoint(layout.formulaField, point) || containsPoint(layout.bailout, point) || (realField.enabled() && containsPoint(layout.realField, point)) || (imaginaryField.enabled() && containsPoint(layout.imaginaryField, point))) {
            SetCursor(LoadCursorW(nullptr, hoveredFormulaToken >= 0 ? IDC_HAND : IDC_IBEAM));
        } else if (containsPoint(layout.picker, point)) {
            SetCursor(LoadCursorW(nullptr, IDC_CROSS));
        } else if (hoverHit != HIT_NONE || containsPoint(layout.preset, point) || scrollbar.contains(point)) {
            SetCursor(LoadCursorW(nullptr, IDC_HAND));
        }
    }

    void onLeftButtonUp(POINT point) {
        if (draggingText) {
            draggingText->mouseUp();
            draggingText = nullptr;
            if (GetCapture() == hwnd) ReleaseCapture();
            return;
        }
        if (draggingPicker) {
            updatePickerFromMouse(point);
            draggingPicker = false;
            if (GetCapture() == hwnd) ReleaseCapture();
            return;
        }
        if (scrollbar.dragging()) {
            scrollbar.mouseUp();
            if (GetCapture() == hwnd) ReleaseCapture();
            InvalidateRect(hwnd, nullptr, FALSE);
            return;
        }

        int hit = hitRouter.release(point.x, point.y);
        int invoke = hit == pressedHit ? hit : HIT_NONE;
        pressedHit = HIT_NONE;
        if (GetCapture() == hwnd) ReleaseCapture();
        InvalidateRect(hwnd, nullptr, FALSE);
        if (invoke != HIT_NONE) invokeHit(invoke);
    }

    bool keyDown(UINT key) {
        bool shift = (GetKeyState(VK_SHIFT) & 0x8000) != 0;
        bool alt = (GetKeyState(VK_MENU) & 0x8000) != 0;
        if (key == VK_TAB) {
            focusNext(shift);
            return true;
        }
        if (focusedId == FOCUS_PRESET) {
            bool wasOpen = presetDropdown.open();
            if (!presetDropdown.keyDown(key, alt, layout.contentClip)) {
                // Continue with the remaining panel shortcuts.
            } else {
                notifyDropdownTransition(wasOpen);
                notifyAccessibility(EVENT_OBJECT_VALUECHANGE, FORMULA_ACC_PRESET);
                InvalidateRect(hwnd, nullptr, FALSE);
                return true;
            }
        }
        if (key == VK_ESCAPE && presetDropdown.open()) {
            bool wasOpen = presetDropdown.open();
            presetDropdown.close();
            notifyDropdownTransition(wasOpen);
            InvalidateRect(hwnd, nullptr, FALSE);
            return true;
        }
        if (focusedId == FOCUS_SCROLLBAR && scrollbar.visible()) {
            int line = 48;
            int page = std::max(line, clientHeightDip() - HEADER_HEIGHT - FOOTER_HEIGHT - BODY_PADDING * 2);
            if (key == VK_UP)
                scrollbar.scrollBy(-line);
            else if (key == VK_DOWN)
                scrollbar.scrollBy(line);
            else if (key == VK_PRIOR)
                scrollbar.scrollBy(-page);
            else if (key == VK_NEXT)
                scrollbar.scrollBy(page);
            else if (key == VK_HOME)
                scrollbar.setPosition(0);
            else if (key == VK_END)
                scrollbar.setPosition(scrollbar.maximumPosition());
            else
                return false;
            return true;
        }
        if (focusedId == HIT_PICKER && selectedEditable(working) && (key == VK_LEFT || key == VK_RIGHT || key == VK_UP || key == VK_DOWN)) {
            std::complex<double>* value = selectedValue(working);
            if (!value) return true;
            double step = pickerRange / 50.0;
            if (key == VK_LEFT) value->real(value->real() - step);
            if (key == VK_RIGHT) value->real(value->real() + step);
            if (key == VK_UP) value->imag(value->imag() + step);
            if (key == VK_DOWN) value->imag(value->imag() - step);
            syncInspectorFields();
            setStatus(std::wstring(inspectorName(selected)) + L" changed; changes are staged.", false);
            notifyAccessibility(EVENT_OBJECT_VALUECHANGE, FORMULA_ACC_PICKER);
            notifyAccessibility(EVENT_OBJECT_VALUECHANGE, FORMULA_ACC_REAL);
            notifyAccessibility(EVENT_OBJECT_VALUECHANGE, FORMULA_ACC_IMAGINARY);
            triggerPreview();
            return true;
        }
        if (key == VK_RETURN || key == VK_SPACE) {
            if (focusedId != FOCUS_FORMULA && focusedId != FOCUS_BAILOUT && focusedId != FOCUS_REAL && focusedId != FOCUS_IMAGINARY && focusedId != FOCUS_PRESET) {
                invokeHit(focusedId);
                return true;
            }
        }
        return false;
    }

    LRESULT handleMessage(UINT message, WPARAM wp, LPARAM lp) {
        switch (message) {
        case WM_GETOBJECT:
            if (static_cast<LONG>(lp) == OBJID_CLIENT) return formulaEditorAccessibilityGetObject(accessibility, wp, lp);
            break;
        case WM_FORMULA_ACCESSIBILITY_SNAPSHOT: {
            if (reinterpret_cast<FormulaEditorAccessibilityProvider*>(wp) != accessibility || !lp) { return FALSE; }
            auto* request = reinterpret_cast<FormulaAccessibilitySnapshotRequest*>(lp);
            if (!request->items) return FALSE;
            std::shared_ptr<AccessibilityGuard> guard = accessibilityGuard;
            if (!guard || guard->snapshotDepth != 0) return FALSE;
            ++guard->snapshotDepth;
            accessibleSnapshot(*request->items);
            --guard->snapshotDepth;
            return TRUE;
        }
        case WM_FORMULA_ACCESSIBILITY_ACTION: {
            if (reinterpret_cast<FormulaEditorAccessibilityProvider*>(wp) != accessibility || !lp) { return FALSE; }
            auto* request = reinterpret_cast<FormulaAccessibilityActionRequest*>(lp);
            std::shared_ptr<AccessibilityGuard> guard = accessibilityGuard;
            if (!guard || guard->action) return FALSE;
            guard->action = true;
            request->handled = performAccessibilityAction(request->key, request->action);
            guard->action = false;
            return request->handled ? TRUE : FALSE;
        }
        case WM_CREATE: return onCreate() ? 0 : -1;
        case WM_DESTROY:
            KillTimer(hwnd, CARET_TIMER);
            destroyFormulaEditorAccessibility(accessibility);
            destroyControls();
            return 0;
        case WM_SIZE:
            updateLayout();
            InvalidateRect(hwnd, nullptr, TRUE);
            return 0;
        case WM_PAINT: paint(); return 0;
        case WM_ERASEBKGND: return 1;
        case WM_TIMER:
            if (wp == CARET_TIMER) {
                InvalidateRect(hwnd, &layout.body, FALSE);
                return 0;
            }
            break;
        case WM_MOUSEWHEEL:
            if (scrollbar.wheel(GET_WHEEL_DELTA_WPARAM(wp), 48)) { return 0; }
            break;
        case WM_MOUSEMOVE: onMouseMove({GET_X_LPARAM(lp), GET_Y_LPARAM(lp)}); return 0;
        case WM_MOUSELEAVE:
            trackingMouse = false;
            hoverHit = HIT_NONE;
            hoveredFormulaToken = -1;
            hitRouter.cancel();
            presetDropdown.mouseMove({-32000, -32000}, layout.contentClip);
            scrollbar.mouseMove({-32000, -32000});
            InvalidateRect(hwnd, nullptr, FALSE);
            return 0;
        case WM_LBUTTONDOWN: onLeftButtonDown({GET_X_LPARAM(lp), GET_Y_LPARAM(lp)}); return 0;
        case WM_LBUTTONDBLCLK: onLeftButtonDoubleClick({GET_X_LPARAM(lp), GET_Y_LPARAM(lp)}); return 0;
        case WM_LBUTTONUP: onLeftButtonUp({GET_X_LPARAM(lp), GET_Y_LPARAM(lp)}); return 0;
        case WM_CAPTURECHANGED:
            if (draggingText) draggingText->mouseUp();
            draggingText = nullptr;
            draggingPicker = false;
            scrollbar.mouseUp();
            pressedHit = HIT_NONE;
            InvalidateRect(hwnd, nullptr, FALSE);
            return 0;
        case WM_KEYDOWN:
        case WM_SYSKEYDOWN:
            if (keyDown(static_cast<UINT>(wp))) return 0;
            break;
        case WM_GETDLGCODE: return DLGC_WANTTAB | DLGC_WANTARROWS | DLGC_WANTCHARS;
        case WM_SETFOCUS:
            if (ui::TextField* field = textFieldForFocus(focusedId)) {
                if (field->enabled()) field->focus();
            }
            InvalidateRect(hwnd, nullptr, FALSE);
            return 0;
        case WM_KILLFOCUS:
            if (presetDropdown.open()) {
                bool wasOpen = presetDropdown.open();
                presetDropdown.close();
                notifyDropdownTransition(wasOpen);
                if (GetCapture() == hwnd) ReleaseCapture();
            }
            InvalidateRect(hwnd, nullptr, FALSE);
            return 0;
        case WM_CLOSE: closePanel(); return 0;
        case WM_NCDESTROY: {
            HWND destroyed = hwnd;
            destroyFormulaEditorAccessibility(accessibility);
            SetWindowLongPtrW(destroyed, GWLP_USERDATA, 0);
            LRESULT result = DefWindowProcW(destroyed, message, wp, lp);
            hwnd = nullptr;
            return result;
        }
        }
        return DefWindowProcW(hwnd, message, wp, lp);
    }

    static LRESULT CALLBACK windowProc(HWND window, UINT message, WPARAM wp, LPARAM lp) {
        Impl* self = reinterpret_cast<Impl*>(GetWindowLongPtrW(window, GWLP_USERDATA));
        if (message == WM_NCCREATE) {
            CREATESTRUCTW* create = reinterpret_cast<CREATESTRUCTW*>(lp);
            self = static_cast<Impl*>(create->lpCreateParams);
            self->hwnd = window;
            SetWindowLongPtrW(window, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(self));
        }
        if (!self) return DefWindowProcW(window, message, wp, lp);
        return self->handleMessage(message, wp, lp);
    }
};

FormulaEditorPanel::FormulaEditorPanel() : _impl(std::make_unique<Impl>()) {
}

FormulaEditorPanel::~FormulaEditorPanel() = default;

bool FormulaEditorPanel::create(HWND owner, int dpi, FormulaEditorCallbacks callbacks) {
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
    SetWindowPos(_impl->hwnd, nullptr, bounds.left, bounds.top, std::max(0L, bounds.right - bounds.left), std::max(0L, bounds.bottom - bounds.top), SWP_NOZORDER | SWP_NOACTIVATE);
}

void FormulaEditorPanel::setDpi(int dpi) {
    _impl->setDpi(dpi);
}

void FormulaEditorPanel::setConfig(const FormulaDialogConfig& config) {
    if (_impl->hwnd) _impl->setConfig(config);
}

void FormulaEditorPanel::dismissPopups() {
    _impl->dismissPopups();
}

HWND FormulaEditorPanel::hwnd() const {
    return _impl->hwnd;
}
