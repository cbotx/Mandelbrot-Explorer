#include "formula_dialog.h"

#include <algorithm>
#include <array>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <vector>

#include "formula_expression.h"

namespace {

enum : int {
    ID_SOURCE = 100,
    ID_PLANE,
    ID_Z0_RE,
    ID_Z0_IM,
    ID_C_RE,
    ID_C_IM,
    ID_BAILOUT,
    ID_PARAM_BASE = 200,
    ID_APPLY = 400,
    ID_MANDELBROT,
    ID_CANCEL
};

struct DialogState {
    HWND hwnd = nullptr;
    HWND owner = nullptr;
    HWND source = nullptr, plane = nullptr, z0Re = nullptr, z0Im = nullptr;
    HWND cRe = nullptr, cIm = nullptr, bailout = nullptr;
    std::array<HWND, 16> parameters{};
    HFONT font = nullptr;
    FormulaDialogConfig config;
    FormulaDialogResult result = FormulaDialogResult::Cancel;
    bool done = false;
};

std::wstring widen(const std::string& value) {
    if (value.empty()) return {};
    int count = MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS,
                                    value.data(), (int)value.size(), nullptr, 0);
    if (count <= 0) return {};
    std::wstring result((size_t)count, L'\0');
    MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, value.data(),
                        (int)value.size(), result.data(), count);
    return result;
}

std::string narrow(const std::wstring& value) {
    if (value.empty()) return {};
    int count = WideCharToMultiByte(CP_UTF8, WC_ERR_INVALID_CHARS,
                                    value.data(), (int)value.size(),
                                    nullptr, 0, nullptr, nullptr);
    if (count <= 0) return {};
    std::string result((size_t)count, '\0');
    WideCharToMultiByte(CP_UTF8, WC_ERR_INVALID_CHARS, value.data(),
                        (int)value.size(), result.data(), count, nullptr, nullptr);
    return result;
}

std::wstring getText(HWND control) {
    int count = GetWindowTextLengthW(control);
    std::wstring result((size_t)count + 1, L'\0');
    GetWindowTextW(control, result.data(), count + 1);
    result.resize((size_t)count);
    return result;
}

void setNumber(HWND control, double value) {
    wchar_t text[64];
    swprintf_s(text, L"%.17g", value);
    SetWindowTextW(control, text);
}

bool parseNumber(HWND control, double& value) {
    std::wstring text = getText(control);
    std::string utf8 = narrow(text);
    if (utf8.empty()) return false;
    char* end = nullptr;
    errno = 0;
    double parsed = std::strtod(utf8.c_str(), &end);
    while (end && (*end == ' ' || *end == '\t')) ++end;
    if (errno == ERANGE || !end || *end != '\0' || !std::isfinite(parsed))
        return false;
    value = parsed;
    return true;
}

HWND addControl(DialogState* state, const wchar_t* cls, const wchar_t* text,
                DWORD style, int x, int y, int w, int h, int id = 0,
                DWORD exStyle = 0) {
    HWND control = CreateWindowExW(exStyle, cls, text, WS_CHILD | WS_VISIBLE | style,
        x, y, w, h, state->hwnd, (HMENU)(INT_PTR)id,
        (HINSTANCE)GetWindowLongPtrW(state->hwnd, GWLP_HINSTANCE), nullptr);
    SendMessageW(control, WM_SETFONT, (WPARAM)state->font, TRUE);
    return control;
}

void addLabel(DialogState* state, const wchar_t* text, int x, int y, int w, int h = 22) {
    addControl(state, L"STATIC", text, SS_LEFT, x, y, w, h);
}

HWND addEdit(DialogState* state, int id, int x, int y, int w, int h = 26) {
    return addControl(state, L"EDIT", L"", WS_TABSTOP | ES_AUTOHSCROLL,
                      x, y, w, h, id, WS_EX_CLIENTEDGE);
}

bool readConfig(DialogState* state) {
    FormulaDialogConfig candidate = state->config;
    candidate.source = narrow(getText(state->source));
    formula::ExpressionProgram program;
    formula::ExpressionError error;
    if (!program.compile(candidate.source, &error)) {
        wchar_t message[384];
        swprintf_s(message, L"Formula error at character %zu:\n%hs",
                   error.position + 1, error.message.c_str());
        MessageBoxW(state->hwnd, message, L"Invalid formula", MB_OK | MB_ICONERROR);
        SetFocus(state->source);
        return false;
    }
    candidate.pixelParameter =
        SendMessageW(state->plane, CB_GETCURSEL, 0, 0) == 1
            ? FormulaParameter::InitialZ : FormulaParameter::C;
    double z0r, z0i, cr, ci;
    if (!parseNumber(state->z0Re, z0r) ||
        !parseNumber(state->z0Im, z0i) ||
        !parseNumber(state->cRe, cr) ||
        !parseNumber(state->cIm, ci) ||
        !parseNumber(state->bailout, candidate.bailout) ||
        !(candidate.bailout > 0.0)) {
        MessageBoxW(state->hwnd, L"Fixed values and bailout must be finite numbers; bailout must be positive.",
                    L"Invalid numeric value", MB_OK | MB_ICONERROR);
        return false;
    }
    candidate.fixedZ0 = { z0r, z0i };
    candidate.fixedC = { cr, ci };
    for (int i = 0; i < 8; ++i) {
        double re, im;
        if (!parseNumber(state->parameters[2 * i], re) ||
            !parseNumber(state->parameters[2 * i + 1], im)) {
            MessageBoxW(state->hwnd, L"Every parameter component must be a finite number.",
                        L"Invalid parameter", MB_OK | MB_ICONERROR);
            return false;
        }
        candidate.parameters[i] = { re, im };
    }
    state->config = std::move(candidate);
    return true;
}

LRESULT CALLBACK FormulaWndProc(HWND hwnd, UINT message, WPARAM wp, LPARAM lp) {
    DialogState* state = (DialogState*)GetWindowLongPtrW(hwnd, GWLP_USERDATA);
    if (message == WM_NCCREATE) {
        state = (DialogState*)((CREATESTRUCTW*)lp)->lpCreateParams;
        state->hwnd = hwnd;
        SetWindowLongPtrW(hwnd, GWLP_USERDATA, (LONG_PTR)state);
    }
    if (!state) return DefWindowProcW(hwnd, message, wp, lp);
    switch (message) {
    case WM_COMMAND:
        switch (LOWORD(wp)) {
        case ID_APPLY:
            if (readConfig(state)) {
                state->result = FormulaDialogResult::Apply;
                state->done = true;
                DestroyWindow(hwnd);
            }
            return 0;
        case ID_MANDELBROT:
            state->result = FormulaDialogResult::UseMandelbrot;
            state->done = true;
            DestroyWindow(hwnd);
            return 0;
        case ID_CANCEL:
            state->done = true;
            DestroyWindow(hwnd);
            return 0;
        }
        break;
    case WM_CLOSE:
        state->done = true;
        DestroyWindow(hwnd);
        return 0;
    case WM_DESTROY:
        if (state->font) { DeleteObject(state->font); state->font = nullptr; }
        EnableWindow(state->owner, TRUE);
        state->done = true;
        return 0;
    }
    return DefWindowProcW(hwnd, message, wp, lp);
}

} // namespace

FormulaDialogResult showFormulaDialog(HWND owner, FormulaDialogConfig& config) {
    static bool registered = false;
    HINSTANCE instance = (HINSTANCE)GetWindowLongPtrW(owner, GWLP_HINSTANCE);
    if (!registered) {
        WNDCLASSW wc{};
        wc.lpfnWndProc = FormulaWndProc;
        wc.hInstance = instance;
        wc.hCursor = LoadCursor(nullptr, IDC_ARROW);
        wc.hbrBackground = (HBRUSH)(COLOR_WINDOW + 1);
        wc.lpszClassName = L"MandelFormulaDialog";
        if (!RegisterClassW(&wc) && GetLastError() != ERROR_CLASS_ALREADY_EXISTS)
            return FormulaDialogResult::Cancel;
        registered = true;
    }

    DialogState state;
    state.owner = owner;
    state.config = config;
    UINT dpi = GetDpiForWindow(owner); if (!dpi) dpi = 96;
    auto S = [dpi](int value) { return MulDiv(value, dpi, 96); };
    state.font = CreateFontW(-S(14), 0, 0, 0, FW_NORMAL, FALSE, FALSE, FALSE,
        DEFAULT_CHARSET, OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS,
        CLEARTYPE_QUALITY, DEFAULT_PITCH | FF_DONTCARE, L"Segoe UI");
    DWORD style = WS_POPUP | WS_CAPTION | WS_SYSMENU;
    RECT frame{ 0, 0, S(680), S(590) };
    AdjustWindowRectExForDpi(&frame, style, FALSE, 0, dpi);
    RECT ownerRect; GetWindowRect(owner, &ownerRect);
    int width = frame.right - frame.left, height = frame.bottom - frame.top;
    HWND hwnd = CreateWindowExW(WS_EX_DLGMODALFRAME, L"MandelFormulaDialog",
        L"Formula editor", style | WS_VISIBLE,
        ownerRect.left + S(80), ownerRect.top + S(40), width, height,
        owner, nullptr, instance, &state);
    if (!hwnd) { DeleteObject(state.font); return FormulaDialogResult::Cancel; }

    int margin = S(18), labelW = S(92), fieldW = S(540), rowH = S(28), y = S(16);
    addLabel(&state, L"Formula z' =", margin, y + S(4), labelW);
    state.source = addEdit(&state, ID_SOURCE, margin + labelW, y, fieldW, rowH);
    SetWindowTextW(state.source, widen(config.source).c_str());
    y += S(40);
    addLabel(&state, L"Pixel axis", margin, y + S(4), labelW);
    state.plane = addControl(&state, L"COMBOBOX", L"",
        WS_TABSTOP | CBS_DROPDOWNLIST, margin + labelW, y, S(220), S(180), ID_PLANE);
    SendMessageW(state.plane, CB_ADDSTRING, 0, (LPARAM)L"c  (parameter plane)");
    SendMessageW(state.plane, CB_ADDSTRING, 0, (LPARAM)L"z0  (dynamical plane)");
    SendMessageW(state.plane, CB_SETCURSEL,
        config.pixelParameter == FormulaParameter::InitialZ ? 1 : 0, 0);
    y += S(42);

    int reX = margin + labelW, imX = reX + S(230), numW = S(210);
    addLabel(&state, L"Fixed z0", margin, y + S(4), labelW);
    state.z0Re = addEdit(&state, ID_Z0_RE, reX, y, numW);
    state.z0Im = addEdit(&state, ID_Z0_IM, imX, y, numW);
    setNumber(state.z0Re, config.fixedZ0.real()); setNumber(state.z0Im, config.fixedZ0.imag());
    y += S(34);
    addLabel(&state, L"Fixed c", margin, y + S(4), labelW);
    state.cRe = addEdit(&state, ID_C_RE, reX, y, numW);
    state.cIm = addEdit(&state, ID_C_IM, imX, y, numW);
    setNumber(state.cRe, config.fixedC.real()); setNumber(state.cIm, config.fixedC.imag());
    y += S(34);
    addLabel(&state, L"Bailout |z|", margin, y + S(4), labelW);
    state.bailout = addEdit(&state, ID_BAILOUT, reX, y, numW);
    setNumber(state.bailout, config.bailout);
    y += S(44);

    addLabel(&state, L"Parameters", margin, y, labelW);
    addLabel(&state, L"real", reX, y, numW);
    addLabel(&state, L"imaginary", imX, y, numW);
    y += S(24);
    for (int i = 0; i < 8; ++i) {
        wchar_t label[8]; swprintf_s(label, L"p%d", i);
        addLabel(&state, label, margin + S(46), y + S(4), S(40));
        state.parameters[2 * i] = addEdit(&state, ID_PARAM_BASE + 2 * i, reX, y, numW);
        state.parameters[2 * i + 1] = addEdit(&state, ID_PARAM_BASE + 2 * i + 1, imX, y, numW);
        setNumber(state.parameters[2 * i], config.parameters[i].real());
        setNumber(state.parameters[2 * i + 1], config.parameters[i].imag());
        y += S(32);
    }

    int buttonY = S(538), buttonW = S(142), gap = S(12);
    addControl(&state, L"BUTTON", L"Apply formula", WS_TABSTOP | BS_DEFPUSHBUTTON,
               margin, buttonY, buttonW, S(34), ID_APPLY);
    addControl(&state, L"BUTTON", L"Use Mandelbrot", WS_TABSTOP,
               margin + buttonW + gap, buttonY, buttonW, S(34), ID_MANDELBROT);
    addControl(&state, L"BUTTON", L"Cancel", WS_TABSTOP,
               margin + 2 * (buttonW + gap), buttonY, buttonW, S(34), ID_CANCEL);

    EnableWindow(owner, FALSE);
    MSG message{};
    int quitCode = -1;
    while (!state.done) {
        BOOL got = GetMessageW(&message, nullptr, 0, 0);
        if (got <= 0) {
            if (got == 0) quitCode = (int)message.wParam;
            if (IsWindow(hwnd)) DestroyWindow(hwnd);
            break;
        }
        if (!IsDialogMessageW(hwnd, &message)) {
            TranslateMessage(&message);
            DispatchMessageW(&message);
        }
    }
    if (IsWindow(hwnd)) DestroyWindow(hwnd);
    EnableWindow(owner, TRUE);
    SetForegroundWindow(owner);
    if (quitCode >= 0) PostQuitMessage(quitCode);
    if (state.result == FormulaDialogResult::Apply) config = std::move(state.config);
    return state.result;
}
