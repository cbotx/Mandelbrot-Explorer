#include "ui_framework.h"

#include "gui_theme.h"

#include <imm.h>

#include <algorithm>
#include <cstdint>
#include <limits>
#include <new>
#include <string>
#include <utility>

#pragma comment(lib, "Imm32.lib")

namespace {

template <typename T>
class UniqueGdiObject {
public:
    UniqueGdiObject() = default;
    explicit UniqueGdiObject(T object, bool owned = true)
        : _object(object), _owned(owned) {}
    ~UniqueGdiObject() { reset(); }

    UniqueGdiObject(const UniqueGdiObject&) = delete;
    UniqueGdiObject& operator=(const UniqueGdiObject&) = delete;

    T get() const { return _object; }
    bool owned() const { return _owned; }

    T release() {
        T object = _object;
        _object = nullptr;
        _owned = false;
        return object;
    }

    void reset(T object = nullptr, bool owned = true) {
        if (_object && _owned) {
            DeleteObject(_object);
        }
        _object = object;
        _owned = owned;
    }

private:
    T _object = nullptr;
    bool _owned = false;
};

class UniqueDc {
public:
    UniqueDc() = default;
    explicit UniqueDc(HDC dc) : _dc(dc) {}
    ~UniqueDc() {
        if (_dc) {
            DeleteDC(_dc);
        }
    }

    UniqueDc(const UniqueDc&) = delete;
    UniqueDc& operator=(const UniqueDc&) = delete;

    HDC get() const { return _dc; }

    HDC release() {
        HDC dc = _dc;
        _dc = nullptr;
        return dc;
    }

private:
    HDC _dc = nullptr;
};

class ScopedSelection {
public:
    ScopedSelection(HDC dc, HGDIOBJ object) : _dc(dc) {
        if (_dc && object) {
            _old = SelectObject(_dc, object);
            _valid = _old && _old != HGDI_ERROR;
        }
    }

    ~ScopedSelection() {
        if (_valid) {
            SelectObject(_dc, _old);
        }
    }

    ScopedSelection(const ScopedSelection&) = delete;
    ScopedSelection& operator=(const ScopedSelection&) = delete;

    bool valid() const { return _valid; }

    HGDIOBJ detach() {
        if (!_valid) {
            return nullptr;
        }
        _valid = false;
        return _old;
    }

private:
    HDC _dc = nullptr;
    HGDIOBJ _old = nullptr;
    bool _valid = false;
};

void deleteFont(HFONT font, bool owned) {
    if (font && owned) DeleteObject(font);
}

HFONT createUiFont(int dpi, int dipHeight, int weight,
                   const wchar_t* faceName) {
    int height = MulDiv(dipHeight, dpi, 96);
    if (height <= 0) {
        height = 1;
    }

    return CreateFontW(-height, 0, 0, 0, weight, FALSE, FALSE, FALSE,
                       DEFAULT_CHARSET, OUT_DEFAULT_PRECIS,
                       CLIP_DEFAULT_PRECIS, CLEARTYPE_QUALITY,
                       DEFAULT_PITCH | FF_DONTCARE, faceName);
}

UniqueGdiObject<HFONT> createFontOrFallback(int dpi, int dipHeight,
                                            int weight,
                                            const wchar_t* faceName) {
    HFONT font = createUiFont(dpi, dipHeight, weight, faceName);
    if (font) {
        return UniqueGdiObject<HFONT>(font, true);
    }

    font = static_cast<HFONT>(GetStockObject(DEFAULT_GUI_FONT));
    return UniqueGdiObject<HFONT>(font, false);
}

COLORREF blendColor(COLORREF from, COLORREF to, unsigned amount) {
    amount = (std::min)(amount, 255u);
    const unsigned inverse = 255u - amount;
    const auto channel = [inverse, amount](BYTE a, BYTE b) {
        return static_cast<BYTE>((static_cast<unsigned>(a) * inverse +
                                  static_cast<unsigned>(b) * amount + 127u) /
                                 255u);
    };
    return RGB(channel(GetRValue(from), GetRValue(to)),
               channel(GetGValue(from), GetGValue(to)),
               channel(GetBValue(from), GetBValue(to)));
}

constexpr wchar_t kTextServiceProperty[] =
    L"MandelbrotExplorer.Ui.TextService.Context";
constexpr size_t kNativeEditLimit = 0x7ffffffeu;

struct TextServiceContext {
    ui::TextService* owner = nullptr;
    WNDPROC originalProc = nullptr;
    std::wstring lastText;
    DWORD selectionFirst = 0;
    DWORD selectionLast = 0;
    unsigned mutationDepth = 0;
};

TextServiceContext* textContext(HWND window) {
    return static_cast<TextServiceContext*>(
        GetPropW(window, kTextServiceProperty));
}

LRESULT callOriginal(TextServiceContext* context, HWND window, UINT message,
                     WPARAM wp, LPARAM lp) {
    if (context && context->originalProc) {
        return CallWindowProcW(context->originalProc, window, message, wp, lp);
    }
    return DefWindowProcW(window, message, wp, lp);
}

bool readWindowText(HWND window, std::wstring& result) {
    const int length = GetWindowTextLengthW(window);
    if (length < 0) {
        return false;
    }

    try {
        std::wstring buffer(static_cast<size_t>(length) + 1u, L'\0');
        const int copied =
            GetWindowTextW(window, buffer.data(), length + 1);
        if (copied < 0) {
            return false;
        }
        buffer.resize(static_cast<size_t>(copied));
        result.swap(buffer);
        return true;
    } catch (const std::bad_alloc&) {
        return false;
    }
}

std::pair<DWORD, DWORD> readSelection(HWND window) {
    DWORD first = 0;
    DWORD last = 0;
    SendMessageW(window, EM_GETSEL, reinterpret_cast<WPARAM>(&first),
                 reinterpret_cast<LPARAM>(&last));
    if (first > last) {
        std::swap(first, last);
    }
    return {first, last};
}

struct SnapshotChanges {
    bool text = false;
    bool selection = false;
};

SnapshotChanges updateSnapshot(HWND window, TextServiceContext* context) {
    SnapshotChanges changes;
    if (!context) {
        return changes;
    }

    std::wstring currentText;
    if (readWindowText(window, currentText)) {
        changes.text = currentText != context->lastText;
        context->lastText.swap(currentText);
    }

    const auto currentSelection = readSelection(window);
    changes.selection =
        currentSelection.first != context->selectionFirst ||
        currentSelection.second != context->selectionLast;
    context->selectionFirst = currentSelection.first;
    context->selectionLast = currentSelection.second;
    return changes;
}

bool replacementFits(HWND window, size_t replacementLength,
                     size_t maximumLength) {
    const int nativeLength = GetWindowTextLengthW(window);
    const size_t currentLength =
        nativeLength > 0 ? static_cast<size_t>(nativeLength) : 0u;
    const auto selected = readSelection(window);
    const size_t first =
        (std::min)(static_cast<size_t>(selected.first), currentLength);
    const size_t last =
        (std::min)(static_cast<size_t>(selected.second), currentLength);
    const size_t retainedLength = currentLength - (last - first);

    return retainedLength <= maximumLength &&
           replacementLength <= maximumLength - retainedLength;
}

bool clipboardTextLength(HWND owner, size_t& length) {
    length = 0;
    if (!OpenClipboard(owner)) {
        return false;
    }

    HGLOBAL data = static_cast<HGLOBAL>(GetClipboardData(CF_UNICODETEXT));
    if (data) {
        const SIZE_T bytes = GlobalSize(data);
        const wchar_t* text = static_cast<const wchar_t*>(GlobalLock(data));
        if (!text) {
            CloseClipboard();
            return false;
        }
        const size_t capacity = static_cast<size_t>(
            bytes / static_cast<SIZE_T>(sizeof(wchar_t)));
        while (length < capacity && text[length] != L'\0') ++length;
        GlobalUnlock(data);
    } else {
        data = static_cast<HGLOBAL>(GetClipboardData(CF_TEXT));
        if (!data) {
            CloseClipboard();
            return false;
        }
        const SIZE_T bytes = GlobalSize(data);
        const char* text = static_cast<const char*>(GlobalLock(data));
        if (!text) {
            CloseClipboard();
            return false;
        }
        size_t byteLength = 0;
        while (byteLength < bytes && text[byteLength] != '\0') ++byteLength;
        int wideLength = MultiByteToWideChar(
            CP_ACP, 0, text, static_cast<int>(byteLength), nullptr, 0);
        GlobalUnlock(data);
        if (wideLength < 0) {
            CloseClipboard();
            return false;
        }
        length = static_cast<size_t>(wideLength);
    }
    CloseClipboard();
    return true;
}

bool imeResultLength(HWND window, size_t& length) {
    length = 0;
    HIMC context = ImmGetContext(window);
    if (!context) {
        return false;
    }

    const LONG bytes =
        ImmGetCompositionStringW(context, GCS_RESULTSTR, nullptr, 0);
    ImmReleaseContext(window, context);
    if (bytes < 0) {
        return false;
    }

    length = static_cast<size_t>(bytes) / sizeof(wchar_t);
    return true;
}

bool isMutationMessage(UINT message, LPARAM lp) {
    switch (message) {
    case WM_CHAR:
    case WM_PASTE:
    case WM_CUT:
    case WM_CLEAR:
    case WM_UNDO:
    case EM_REPLACESEL:
    case WM_SETTEXT:
        return true;
    case WM_IME_COMPOSITION:
        return (lp & GCS_RESULTSTR) != 0;
    default:
        return false;
    }
}

bool isNavigationMessage(UINT message, WPARAM wp) {
    switch (message) {
    case EM_SETSEL:
    case WM_LBUTTONDOWN:
    case WM_LBUTTONUP:
    case WM_LBUTTONDBLCLK:
        return true;
    case WM_MOUSEMOVE:
        return (wp & MK_LBUTTON) != 0;
    case WM_KEYDOWN:
    case WM_SYSKEYDOWN:
        return true;
    default:
        return false;
    }
}

bool replacementLengthForMessage(HWND window, UINT message, WPARAM wp,
                                 LPARAM lp, size_t& replacementLength,
                                 bool& replacesAll) {
    replacementLength = 0;
    replacesAll = false;

    switch (message) {
    case WM_CHAR:
        if (wp < L' ' || wp == 0x7fu) {
            return false;
        }
        replacementLength = 1;
        return true;

    case WM_PASTE:
        return clipboardTextLength(window, replacementLength);

    case EM_REPLACESEL: {
        const wchar_t* replacement =
            reinterpret_cast<const wchar_t*>(lp);
        replacementLength =
            replacement ? std::char_traits<wchar_t>::length(replacement) : 0u;
        return true;
    }

    case WM_SETTEXT: {
        const wchar_t* replacement =
            reinterpret_cast<const wchar_t*>(lp);
        replacementLength =
            replacement ? std::char_traits<wchar_t>::length(replacement) : 0u;
        replacesAll = true;
        return true;
    }

    case WM_IME_COMPOSITION:
        if ((lp & GCS_RESULTSTR) != 0) {
            return imeResultLength(window, replacementLength);
        }
        return false;

    default:
        return false;
    }
}

bool messageExceedsLimit(HWND window, UINT message, WPARAM wp, LPARAM lp,
                         size_t maximumLength) {
    size_t replacementLength = 0;
    bool replacesAll = false;
    if (!replacementLengthForMessage(window, message, wp, lp,
                                     replacementLength, replacesAll)) {
        return false;
    }

    return replacesAll ? replacementLength > maximumLength
                       : !replacementFits(window, replacementLength,
                                          maximumLength);
}

int clampNonnegative(int value) {
    return value < 0 ? 0 : value;
}

} // namespace

namespace ui {

RECT DpiScale::px(RECT dipRect) const {
    return {px(dipRect.left), px(dipRect.top), px(dipRect.right),
            px(dipRect.bottom)};
}

Resources::~Resources() {
    reset();
}

bool Resources::create(int dpi) {
    const int effectiveDpi = dpi > 0 ? dpi : 96;

    auto regular =
        createFontOrFallback(effectiveDpi, 14, FW_NORMAL, L"Segoe UI");
    auto semibold =
        createFontOrFallback(effectiveDpi, 15, FW_SEMIBOLD, L"Segoe UI");
    auto small =
        createFontOrFallback(effectiveDpi, 11, FW_NORMAL, L"Segoe UI");
    auto mono =
        createFontOrFallback(effectiveDpi, 14, FW_NORMAL, L"Consolas");
    UniqueGdiObject<HBRUSH> panel(CreateSolidBrush(CLR_PANEL));
    UniqueGdiObject<HBRUSH> card(CreateSolidBrush(CLR_CARD));

    if (!regular.get() || !semibold.get() || !small.get() || !mono.get() ||
        !panel.get() || !card.get()) {
        return false;
    }

    const bool regularOwned = regular.owned();
    const bool semiboldOwned = semibold.owned();
    const bool smallOwned = small.owned();
    const bool monoOwned = mono.owned();
    reset();
    _scale.setDpi(effectiveDpi);
    _regular = regular.release();
    _semibold = semibold.release();
    _small = small.release();
    _mono = mono.release();
    _regularOwned = regularOwned;
    _semiboldOwned = semiboldOwned;
    _smallOwned = smallOwned;
    _monoOwned = monoOwned;
    _panelBrush = panel.release();
    _cardBrush = card.release();
    return true;
}

void Resources::reset() {
    deleteFont(_regular, _regularOwned);
    deleteFont(_semibold, _semiboldOwned);
    deleteFont(_small, _smallOwned);
    deleteFont(_mono, _monoOwned);
    if (_panelBrush) {
        DeleteObject(_panelBrush);
    }
    if (_cardBrush) {
        DeleteObject(_cardBrush);
    }

    _regular = nullptr;
    _semibold = nullptr;
    _small = nullptr;
    _mono = nullptr;
    _regularOwned = false;
    _semiboldOwned = false;
    _smallOwned = false;
    _monoOwned = false;
    _panelBrush = nullptr;
    _cardBrush = nullptr;
}

BackBuffer::~BackBuffer() {
    reset();
}

HDC BackBuffer::begin(HDC target, int width, int height) {
    if (!target || width <= 0 || height <= 0) {
        return nullptr;
    }
    if (_dc && _bitmap && width == _width && height == _height) {
        return _dc;
    }

    UniqueDc newDc(CreateCompatibleDC(target));
    if (!newDc.get()) {
        return nullptr;
    }

    UniqueGdiObject<HBITMAP> newBitmap(
        CreateCompatibleBitmap(target, width, height));
    if (!newBitmap.get()) {
        return nullptr;
    }

    ScopedSelection selection(newDc.get(), newBitmap.get());
    if (!selection.valid()) {
        return nullptr;
    }

    const HGDIOBJ oldBitmap = selection.detach();
    HDC committedDc = newDc.release();
    HBITMAP committedBitmap = newBitmap.release();

    reset();
    _dc = committedDc;
    _bitmap = committedBitmap;
    _oldBitmap = oldBitmap;
    _width = width;
    _height = height;
    return _dc;
}

void BackBuffer::present(HDC target, const RECT& area) const {
    if (!target || !_dc) {
        return;
    }

    const int width = area.right - area.left;
    const int height = area.bottom - area.top;
    if (width <= 0 || height <= 0) {
        return;
    }

    BitBlt(target, area.left, area.top, width, height, _dc, area.left,
           area.top, SRCCOPY);
}

void BackBuffer::reset() {
    bool bitmapDeselected = !_dc;
    if (_dc && _oldBitmap && _oldBitmap != HGDI_ERROR) {
        const HGDIOBJ selected = SelectObject(_dc, _oldBitmap);
        bitmapDeselected = selected && selected != HGDI_ERROR;
    }
    if (_dc && _bitmap && !bitmapDeselected) {
        DeleteDC(_dc);
        _dc = nullptr;
    }
    if (_bitmap) {
        DeleteObject(_bitmap);
    }
    if (_dc) {
        DeleteDC(_dc);
    }

    _dc = nullptr;
    _bitmap = nullptr;
    _oldBitmap = nullptr;
    _width = 0;
    _height = 0;
}

bool ScrollState::configure(int contentExtent, int viewportExtent) {
    const int newContentExtent = clampNonnegative(contentExtent);
    const int newViewportExtent = clampNonnegative(viewportExtent);
    const int newMaximum =
        newContentExtent > newViewportExtent
            ? newContentExtent - newViewportExtent
            : 0;
    const int newPosition = (std::min)(_position, newMaximum);

    const bool changed = _contentExtent != newContentExtent ||
                         _viewportExtent != newViewportExtent ||
                         _position != newPosition;
    _contentExtent = newContentExtent;
    _viewportExtent = newViewportExtent;
    _position = newPosition;
    return changed;
}

bool ScrollState::setPosition(int position) {
    const int clamped = (std::max)(0, (std::min)(position, maximumPosition()));
    if (clamped == _position) {
        return false;
    }
    _position = clamped;
    return true;
}

bool ScrollState::scrollBy(int delta) {
    const std::int64_t candidate =
        static_cast<std::int64_t>(_position) +
        static_cast<std::int64_t>(delta);
    const std::int64_t clamped =
        (std::max)(std::int64_t{0},
                   (std::min)(candidate,
                              static_cast<std::int64_t>(maximumPosition())));
    return setPosition(static_cast<int>(clamped));
}

bool ScrollState::handleCommand(int command, int trackPosition) {
    const int lineStep = (std::max)(1, _viewportExtent / 12);
    switch (command) {
    case SB_LINEUP:
        return scrollBy(-lineStep);
    case SB_LINEDOWN:
        return scrollBy(lineStep);
    case SB_PAGEUP:
        return scrollBy(-_viewportExtent);
    case SB_PAGEDOWN:
        return scrollBy(_viewportExtent);
    case SB_THUMBPOSITION:
    case SB_THUMBTRACK:
        return setPosition(trackPosition);
    case SB_TOP:
        return setPosition(0);
    case SB_BOTTOM:
        return setPosition(maximumPosition());
    default:
        return false;
    }
}

void ScrollState::apply(HWND window, int bar) const {
    if (!window) {
        return;
    }

    SCROLLINFO info{};
    info.cbSize = sizeof(info);
    info.fMask = SIF_RANGE | SIF_PAGE | SIF_POS;
    info.nMin = 0;
    info.nMax = _contentExtent > 0 ? _contentExtent - 1 : 0;
    info.nPage = static_cast<UINT>(_viewportExtent);
    info.nPos = _position;
    SetScrollInfo(window, bar, &info, TRUE);
}

int ScrollState::maximumPosition() const {
    return _contentExtent > _viewportExtent
               ? _contentExtent - _viewportExtent
               : 0;
}

void HitRouter::clear() {
    _regions.clear();
    _hovered = 0;
    _pressed = 0;
}

void HitRouter::add(int id, RECT bounds, bool enabled) {
    _regions.push_back(HitRegion{id, bounds, enabled});
}

int HitRouter::hit(int x, int y) const {
    for (auto region = _regions.rbegin(); region != _regions.rend(); ++region) {
        if (region->enabled && x >= region->bounds.left &&
            x < region->bounds.right && y >= region->bounds.top &&
            y < region->bounds.bottom) {
            return region->id;
        }
    }
    return 0;
}

bool HitRouter::move(int x, int y) {
    const int newHovered = hit(x, y);
    if (newHovered == _hovered) {
        return false;
    }
    _hovered = newHovered;
    return true;
}

bool HitRouter::press(int x, int y) {
    const int newHovered = hit(x, y);
    const bool changed =
        newHovered != _hovered || newHovered != _pressed;
    _hovered = newHovered;
    _pressed = newHovered;
    return changed;
}

int HitRouter::release(int x, int y) {
    const int released = hit(x, y);
    const int clicked =
        _pressed != 0 && released == _pressed ? released : 0;
    _hovered = released;
    _pressed = 0;
    return clicked;
}

bool HitRouter::cancel() {
    if (_hovered == 0 && _pressed == 0) {
        return false;
    }
    _hovered = 0;
    _pressed = 0;
    return true;
}

void drawButton(HDC dc, RECT bounds, const std::wstring& text, HFONT font,
                ButtonStyle style, bool hovered, bool pressed, bool enabled,
                int radius) {
    if (!dc || bounds.right <= bounds.left || bounds.bottom <= bounds.top) {
        return;
    }

    COLORREF fill = CLR_CARD;
    COLORREF border = CLR_BORDER;
    COLORREF textColor = CLR_TEXT;

    switch (style) {
    case ButtonStyle::Accent:
        fill = CLR_ACCENT;
        border = CLR_ACCENT;
        textColor = RGB(255, 255, 255);
        if (hovered) {
            fill = CLR_ACCENT_HI;
            border = CLR_ACCENT_HI;
        }
        if (pressed) {
            fill = blendColor(CLR_ACCENT, CLR_BG, 52);
            border = fill;
        }
        break;

    case ButtonStyle::Positive:
        fill = CLR_GREEN;
        border = CLR_GREEN;
        textColor = RGB(255, 255, 255);
        if (hovered) {
            fill = blendColor(CLR_GREEN, RGB(255, 255, 255), 38);
            border = fill;
        }
        if (pressed) {
            fill = blendColor(CLR_GREEN, CLR_BG, 52);
            border = fill;
        }
        break;

    case ButtonStyle::Subtle:
        fill = CLR_PANEL;
        border = CLR_PANEL;
        textColor = CLR_TEXT_DIM;
        if (hovered) {
            fill = CLR_CARD;
            border = CLR_BORDER;
            textColor = CLR_TEXT;
        }
        if (pressed) {
            fill = CLR_CARD_HOV;
            border = CLR_BORDER;
            textColor = CLR_TEXT;
        }
        break;

    case ButtonStyle::Normal:
    default:
        if (hovered) {
            fill = CLR_CARD_HOV;
        }
        if (pressed) {
            fill = blendColor(CLR_CARD_HOV, CLR_BG, 48);
        }
        break;
    }

    if (!enabled) {
        fill = blendColor(fill, CLR_PANEL, 150);
        border = blendColor(border, CLR_PANEL, 160);
        textColor = blendColor(CLR_TEXT_DIM, CLR_PANEL, 55);
    }

    fillRound(dc, bounds, fill, border, (std::max)(0, radius));
    RECT textBounds = bounds;
    if (pressed && enabled) {
        OffsetRect(&textBounds, 0, 1);
    }
    drawText(dc, textBounds, text, textColor, font,
             DT_CENTER | DT_VCENTER | DT_SINGLELINE | DT_END_ELLIPSIS);
}

void drawCard(HDC dc, RECT bounds, int radius) {
    if (!dc || bounds.right <= bounds.left || bounds.bottom <= bounds.top) {
        return;
    }
    fillRound(dc, bounds, CLR_CARD, CLR_BORDER, (std::max)(0, radius));
}

void drawFocusRing(HDC dc, RECT bounds, int radius) {
    if (!dc || bounds.right <= bounds.left || bounds.bottom <= bounds.top) {
        return;
    }

    UniqueGdiObject<HPEN> pen(
        CreatePen(PS_SOLID, 2, CLR_ACCENT_HI));
    HBRUSH hollow =
        static_cast<HBRUSH>(GetStockObject(HOLLOW_BRUSH));
    if (!pen.get() || !hollow) {
        return;
    }

    ScopedSelection selectedPen(dc, pen.get());
    ScopedSelection selectedBrush(dc, hollow);
    if (!selectedPen.valid() || !selectedBrush.valid()) {
        return;
    }

    RoundRect(dc, bounds.left, bounds.top, bounds.right, bounds.bottom,
              (std::max)(0, radius), (std::max)(0, radius));
}

TextService::~TextService() {
    destroy();
}

bool TextService::create(HWND parent, int controlId, size_t maximumLength,
                         Callback callback) {
    destroy();
    if (!parent || !IsWindow(parent)) {
        return false;
    }

    const size_t effectiveMaximum =
        (std::min)(maximumLength, kNativeEditLimit);
    HWND edit = CreateWindowExW(
        0, L"EDIT", L"",
        WS_CHILD | WS_VISIBLE | WS_TABSTOP | ES_LEFT | ES_AUTOHSCROLL,
        -32000, -32000, 1, 1, parent,
        reinterpret_cast<HMENU>(static_cast<INT_PTR>(controlId)),
        GetModuleHandleW(nullptr), nullptr);
    if (!edit) {
        return false;
    }

    TextServiceContext* context =
        new (std::nothrow) TextServiceContext();
    if (!context) {
        DestroyWindow(edit);
        return false;
    }
    context->owner = this;

    if (!SetPropW(edit, kTextServiceProperty,
                  reinterpret_cast<HANDLE>(context))) {
        delete context;
        DestroyWindow(edit);
        return false;
    }

    SetLastError(ERROR_SUCCESS);
    WNDPROC original = reinterpret_cast<WNDPROC>(SetWindowLongPtrW(
        edit, GWLP_WNDPROC, reinterpret_cast<LONG_PTR>(&TextService::editProc)));
    if (!original && GetLastError() != ERROR_SUCCESS) {
        RemovePropW(edit, kTextServiceProperty);
        delete context;
        DestroyWindow(edit);
        return false;
    }

    context->originalProc = original;
    _edit = edit;
    _originalProc = original;
    _callback = std::move(callback);
    _maximumLength = effectiveMaximum;

    const size_t nativeLimit =
        effectiveMaximum == 0 ? 1u : effectiveMaximum;
    SendMessageW(_edit, EM_SETLIMITTEXT,
                 static_cast<WPARAM>(nativeLimit), 0);
    updateSnapshot(_edit, context);
    return true;
}

void TextService::destroy() {
    HWND edit = _edit;
    if (!edit) {
        _originalProc = nullptr;
        _callback = {};
        _maximumLength = 0;
        _syncing = false;
        return;
    }

    _syncing = true;
    TextServiceContext* context =
        IsWindow(edit) ? textContext(edit) : nullptr;
    WNDPROC original =
        context && context->originalProc ? context->originalProc
                                         : _originalProc;
    bool restored = false;

    if (IsWindow(edit) && original) {
        WNDPROC current = reinterpret_cast<WNDPROC>(
            GetWindowLongPtrW(edit, GWLP_WNDPROC));
        if (current == &TextService::editProc) {
            SetLastError(ERROR_SUCCESS);
            const LONG_PTR previous = SetWindowLongPtrW(
                edit, GWLP_WNDPROC, reinterpret_cast<LONG_PTR>(original));
            restored =
                previous != 0 || GetLastError() == ERROR_SUCCESS;
        }
    }

    _edit = nullptr;
    _originalProc = nullptr;

    if (context) {
        context->owner = nullptr;
    }

    if (restored) {
        if (context) {
            if (textContext(edit) == context) {
                RemovePropW(edit, kTextServiceProperty);
            }
            delete context;
            context = nullptr;
        }
    }

    if (IsWindow(edit)) {
        DestroyWindow(edit);
    }

    _callback = {};
    _maximumLength = 0;
    _syncing = false;
}

bool TextService::setText(const std::wstring& value) {
    if (!_edit || !IsWindow(_edit) || value.size() > _maximumLength) {
        if (_edit && value.size() > _maximumLength) {
            MessageBeep(MB_ICONWARNING);
        }
        return false;
    }

    const bool wasSyncing = _syncing;
    _syncing = true;
    const LRESULT result = SendMessageW(
        _edit, WM_SETTEXT, 0, reinterpret_cast<LPARAM>(value.c_str()));
    _syncing = wasSyncing;
    return result != FALSE;
}

std::wstring TextService::text() const {
    std::wstring value;
    if (_edit && IsWindow(_edit)) {
        readWindowText(_edit, value);
    }
    return value;
}

void TextService::setInputAnchor(POINT clientPoint, int lineHeight) {
    if (!_edit || !IsWindow(_edit)) return;
    int height = (std::max)(1, lineHeight);
    SetWindowPos(_edit, nullptr, clientPoint.x, clientPoint.y,
                 1, height, SWP_NOZORDER | SWP_NOACTIVATE);

    HIMC context = ImmGetContext(_edit);
    if (!context) return;
    COMPOSITIONFORM composition{};
    composition.dwStyle = CFS_POINT;
    composition.ptCurrentPos = { 0, height };
    ImmSetCompositionWindow(context, &composition);
    CANDIDATEFORM candidate{};
    candidate.dwIndex = 0;
    candidate.dwStyle = CFS_CANDIDATEPOS;
    candidate.ptCurrentPos = { 0, height };
    ImmSetCandidateWindow(context, &candidate);
    ImmReleaseContext(_edit, context);
}

void TextService::focus() {
    if (_edit && IsWindow(_edit)) {
        SetFocus(_edit);
    }
}

bool TextService::focused() const {
    return _edit && IsWindow(_edit) && GetFocus() == _edit;
}

void TextService::setSelection(size_t first, size_t last) {
    if (!_edit || !IsWindow(_edit)) {
        return;
    }

    const int nativeLength = GetWindowTextLengthW(_edit);
    const size_t length =
        nativeLength > 0 ? static_cast<size_t>(nativeLength) : 0u;
    first = (std::min)(first, length);
    last = (std::min)(last, length);
    const size_t nativeMaximum =
        static_cast<size_t>((std::numeric_limits<LONG>::max)());
    first = (std::min)(first, nativeMaximum);
    last = (std::min)(last, nativeMaximum);

    SendMessageW(_edit, EM_SETSEL, static_cast<WPARAM>(first),
                 static_cast<LPARAM>(last));
}

std::pair<size_t, size_t> TextService::selection() const {
    if (!_edit || !IsWindow(_edit)) {
        return {0u, 0u};
    }

    const auto selected = readSelection(_edit);
    return {static_cast<size_t>(selected.first),
            static_cast<size_t>(selected.second)};
}

bool TextService::replaceSelection(const std::wstring& replacement,
                                   bool allowUndo) {
    if (!_edit || !IsWindow(_edit)) {
        return false;
    }
    if (!replacementFits(_edit, replacement.size(), _maximumLength)) {
        MessageBeep(MB_ICONWARNING);
        return false;
    }

    SendMessageW(_edit, EM_REPLACESEL, allowUndo ? TRUE : FALSE,
                 reinterpret_cast<LPARAM>(replacement.c_str()));
    return true;
}

bool TextService::undo() {
    if (!_edit || !IsWindow(_edit) ||
        SendMessageW(_edit, EM_CANUNDO, 0, 0) == FALSE) {
        return false;
    }
    return SendMessageW(_edit, WM_UNDO, 0, 0) != FALSE;
}

void TextService::selectAll() {
    if (!_edit || !IsWindow(_edit)) {
        return;
    }
    SendMessageW(_edit, EM_SETSEL, 0, -1);
}

LRESULT CALLBACK TextService::editProc(HWND window, UINT message, WPARAM wp,
                                       LPARAM lp) {
    TextServiceContext* context = textContext(window);
    if (!context) {
        return DefWindowProcW(window, message, wp, lp);
    }

    if (message == WM_NCDESTROY) {
        WNDPROC original = context->originalProc;
        const LRESULT result =
            original ? CallWindowProcW(original, window, message, wp, lp)
                     : DefWindowProcW(window, message, wp, lp);

        if (textContext(window) == context) {
            TextService* owner = context->owner;
            if (owner && owner->_edit == window) {
                owner->_edit = nullptr;
                owner->_originalProc = nullptr;
            }
            RemovePropW(window, kTextServiceProperty);
            delete context;
        }
        return result;
    }

    const bool mutation = isMutationMessage(message, lp);
    const bool navigation = isNavigationMessage(message, wp);

    TextService* owner = context->owner;
    if (owner && mutation &&
        messageExceedsLimit(window, message, wp, lp,
                            owner->_maximumLength)) {
        MessageBeep(MB_ICONWARNING);
        return message == WM_SETTEXT ? FALSE : 0;
    }

    const bool outerMutation = mutation && context->mutationDepth == 0;
    if (mutation) {
        ++context->mutationDepth;
    }

    const LRESULT result =
        callOriginal(context, window, message, wp, lp);

    TextServiceContext* current = textContext(window);
    if (current != context) {
        return result;
    }

    if (mutation && current->mutationDepth > 0) {
        --current->mutationDepth;
    }
    if (mutation && !outerMutation) {
        return result;
    }

    if (message == WM_SETFOCUS || message == WM_KILLFOCUS) {
        owner = current->owner;
        if (owner) {
            owner->notify(message == WM_SETFOCUS ? Event::FocusGained
                                                 : Event::FocusLost);
        }
        return result;
    }

    if (mutation || navigation) {
        const SnapshotChanges changes = updateSnapshot(window, current);
        owner = current->owner;
        if (!owner || owner->_syncing) {
            return result;
        }

        TextServiceContext* eventContext = current;
        const DWORD eventSelectionFirst = current->selectionFirst;
        const DWORD eventSelectionLast = current->selectionLast;
        if (changes.text) {
            owner->notify(Event::Changed);
            current = textContext(window);
            if (current != eventContext || current->owner != owner) {
                return result;
            }
        }
        if (changes.selection &&
            current->selectionFirst == eventSelectionFirst &&
            current->selectionLast == eventSelectionLast) {
            owner->notify(Event::CaretMoved);
        }
    }

    return result;
}

void TextService::notify(Event event) {
    if (!_callback) return;
    Callback callback = _callback;
    callback(event);
}

void PopupState::show(int count, int current) {
    itemCount = clampNonnegative(count);
    if (itemCount == 0) {
        close();
        return;
    }

    open = true;
    selected = (std::max)(0, (std::min)(current, itemCount - 1));
    hovered = -1;
}

void PopupState::close() {
    open = false;
    selected = -1;
    hovered = -1;
    itemCount = 0;
}

bool PopupState::move(int delta) {
    if (!open || itemCount <= 0 || delta == 0) {
        return false;
    }

    const int oldSelected = selected;
    const int oldHovered = hovered;
    const std::int64_t count = itemCount;
    std::int64_t next =
        static_cast<std::int64_t>((std::max)(0, selected));
    next = (next + static_cast<std::int64_t>(delta) % count) % count;
    if (next < 0) {
        next += count;
    }

    selected = static_cast<int>(next);
    hovered = -1;
    return selected != oldSelected || hovered != oldHovered;
}

bool PopupState::setHovered(int index) {
    int newHovered = -1;
    if (open && itemCount > 0 && index >= 0) {
        newHovered = (std::min)(index, itemCount - 1);
    }
    if (newHovered == hovered) {
        return false;
    }
    hovered = newHovered;
    return true;
}

int PopupState::accept() {
    if (!open || itemCount <= 0) {
        close();
        return -1;
    }

    const int candidate = hovered >= 0 ? hovered : selected;
    const int accepted =
        (std::max)(0, (std::min)(candidate, itemCount - 1));
    close();
    return accepted;
}

} // namespace ui
