#include "ui_controls.h"

#include <windowsx.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <cwctype>
#include <limits>
#include <utility>

#include "gui_theme.h"

namespace {

constexpr wchar_t TEXT_FIELD_PROPERTY[] = L"MandelbrotExplorer.Ui.TextField.Context";

bool contains(RECT rect, POINT point) {
    return point.x >= rect.left && point.x < rect.right && point.y >= rect.top && point.y < rect.bottom;
}

RECT intersected(RECT a, RECT b) {
    RECT result{};
    if (!IntersectRect(&result, &a, &b)) return {};
    return result;
}

COLORREF blend(COLORREF from, COLORREF to, unsigned amount) {
    amount = std::min(amount, 255u);
    unsigned inverse = 255u - amount;
    auto channel = [inverse, amount](BYTE a, BYTE b) { return static_cast<BYTE>((static_cast<unsigned>(a) * inverse + static_cast<unsigned>(b) * amount + 127u) / 255u); };
    return RGB(channel(GetRValue(from), GetRValue(to)), channel(GetGValue(from), GetGValue(to)), channel(GetBValue(from), GetBValue(to)));
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
    void* target = GlobalLock(memory);
    if (!target) {
        GlobalFree(memory);
        CloseClipboard();
        return false;
    }
    std::memcpy(target, text.c_str(), bytes);
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
    bool ok = false;
    HGLOBAL memory = static_cast<HGLOBAL>(GetClipboardData(CF_UNICODETEXT));
    if (memory) {
        const wchar_t* source = static_cast<const wchar_t*>(GlobalLock(memory));
        SIZE_T capacity = GlobalSize(memory) / sizeof(wchar_t);
        if (source && capacity > 0) {
            size_t length = 0;
            while (length < capacity && source[length]) ++length;
            text.assign(source, length);
            ok = true;
        }
        if (source) GlobalUnlock(memory);
    } else {
        memory = static_cast<HGLOBAL>(GetClipboardData(CF_TEXT));
        if (memory) {
            const char* source = static_cast<const char*>(GlobalLock(memory));
            SIZE_T capacity = GlobalSize(memory);
            if (source && capacity > 0) {
                size_t bytes = 0;
                while (bytes < capacity && source[bytes]) ++bytes;
                int length = MultiByteToWideChar(CP_ACP, 0, source, static_cast<int>(bytes), nullptr, 0);
                if (length >= 0) {
                    text.resize(static_cast<size_t>(length));
                    ok = length == 0 || MultiByteToWideChar(CP_ACP, 0, source, static_cast<int>(bytes), text.data(), length) == length;
                }
            }
            if (source) GlobalUnlock(memory);
        }
    }
    CloseClipboard();
    return ok;
}

void drawLine(HDC dc, int x1, int y1, int x2, int y2, COLORREF color) {
    HPEN pen = CreatePen(PS_SOLID, 1, color);
    if (!pen) return;
    HGDIOBJ old = SelectObject(dc, pen);
    MoveToEx(dc, x1, y1, nullptr);
    LineTo(dc, x2, y2);
    SelectObject(dc, old);
    DeleteObject(pen);
}

} // namespace

namespace ui {

TextField::~TextField() {
    destroy();
}

bool TextField::create(HWND parent, int controlId, size_t maximumLength, Callback callback) {
    destroy();
    _lifetime = std::make_shared<int>(0);
    _parent = parent;
    _callback = std::move(callback);
    _maximumLength = maximumLength;
    _blinkReset = GetTickCount();
    if (!_service.create(parent, controlId, maximumLength, [this](TextService::Event event) { onTextServiceEvent(event); })) {
        _parent = nullptr;
        _callback = {};
        return false;
    }

    HWND edit = _service.hwnd();
    if (!SetPropW(edit, TEXT_FIELD_PROPERTY, reinterpret_cast<HANDLE>(this))) {
        destroy();
        return false;
    }
    SetLastError(ERROR_SUCCESS);
    _editOriginal = reinterpret_cast<WNDPROC>(SetWindowLongPtrW(edit, GWLP_WNDPROC, reinterpret_cast<LONG_PTR>(&TextField::editProc)));
    if (!_editOriginal && GetLastError() != ERROR_SUCCESS) {
        RemovePropW(edit, TEXT_FIELD_PROPERTY);
        destroy();
        return false;
    }
    return true;
}

void TextField::destroy() {
    _lifetime.reset();
    HWND edit = _service.hwnd();
    if (edit && IsWindow(edit)) {
        WNDPROC current = reinterpret_cast<WNDPROC>(GetWindowLongPtrW(edit, GWLP_WNDPROC));
        if (current == &TextField::editProc && _editOriginal) { SetWindowLongPtrW(edit, GWLP_WNDPROC, reinterpret_cast<LONG_PTR>(_editOriginal)); }
        RemovePropW(edit, TEXT_FIELD_PROPERTY);
    }
    _editOriginal = nullptr;
    _service.destroy();
    _parent = nullptr;
    _callback = {};
    _advances.clear();
    _maximumLength = 0;
    _dragging = false;
}

void TextField::setBounds(RECT bounds) {
    _bounds = bounds;
}

void TextField::setEnabled(bool enabled) {
    if (_enabled == enabled) return;
    _enabled = enabled;
    if (!enabled) _dragging = false;
    HWND edit = _service.hwnd();
    if (edit) EnableWindow(edit, enabled ? TRUE : FALSE);
}

void TextField::setPlaceholder(std::wstring placeholder) {
    _placeholder = std::move(placeholder);
}

bool TextField::setText(const std::wstring& value) {
    bool result = _service.setText(value);
    if (result) {
        _activeCaret = std::min(_activeCaret, value.size());
        _scrollX = 0;
        resetBlink();
    }
    return result;
}

std::wstring TextField::text() const {
    return _service.text();
}

void TextField::focus() {
    if (_enabled) _service.focus();
}

bool TextField::focused() const {
    return _service.focused();
}

void TextField::setSelection(size_t first, size_t last) {
    _activeCaret = last;
    _caretPreference = last < first ? -1 : (last > first ? 1 : 0);
    resetBlink();
    _service.setSelection(first, last);
}

std::pair<size_t, size_t> TextField::selection() const {
    return _service.selection();
}

void TextField::selectAll() {
    _activeCaret = text().size();
    resetBlink();
    _service.selectAll();
}

bool TextField::replaceSelection(const std::wstring& replacement, bool allowUndo) {
    resetBlink();
    return _service.replaceSelection(replacement, allowUndo);
}

bool TextField::undo() {
    resetBlink();
    return _service.undo();
}

void TextField::updateAdvances(HDC dc, HFONT font, const std::wstring& value, int availableWidth, const std::vector<TextRangeStyle>& ranges) {
    _advances.assign(value.size() + 1u, 0);
    if (value.empty()) {
        _scrollX = 0;
        return;
    }

    HGDIOBJ oldFont = font ? SelectObject(dc, font) : nullptr;
    std::vector<int> widths(value.size());
    SIZE size{};
    int fit = 0;
    if (GetTextExtentExPointW(dc, value.c_str(), static_cast<int>(value.size()), (std::numeric_limits<int>::max)(), &fit, widths.data(), &size)) {
        for (size_t i = 0; i < value.size(); ++i) _advances[i + 1u] = widths[i];
    } else {
        for (size_t i = 0; i < value.size(); ++i) {
            SIZE character{};
            GetTextExtentPoint32W(dc, value.data() + i, 1, &character);
            _advances[i + 1u] = _advances[i] + character.cx;
        }
    }
    if (oldFont) SelectObject(dc, oldFont);

    std::vector<int> padding(value.size() + 1u, 0);
    for (const TextRangeStyle& range : ranges) {
        size_t first = std::min(range.first, value.size());
        size_t last = std::min(range.last, value.size());
        padding[first] += std::max(0, range.paddingBefore);
        padding[last] += std::max(0, range.paddingAfter);
    }
    int accumulatedPadding = 0;
    for (size_t i = 0; i < _advances.size(); ++i) {
        accumulatedPadding += padding[i];
        _advances[i] += accumulatedPadding;
    }

    size_t caret = std::min(_activeCaret, value.size());
    int caretX = _advances[caret];
    int margin = std::max(3, availableWidth / 16);
    if (caretX - _scrollX > availableWidth - margin) _scrollX = caretX - availableWidth + margin;
    if (caretX - _scrollX < margin) _scrollX = std::max(0, caretX - margin);
    int maximum = std::max(0, _advances.back() - availableWidth + margin);
    _scrollX = std::clamp(_scrollX, 0, maximum);
}

int TextField::xForIndex(size_t index) const {
    if (_advances.empty()) return _inner.left - _scrollX;
    index = std::min(index, _advances.size() - 1u);
    return _inner.left + _advances[index] - _scrollX;
}

void TextField::draw(HDC dc, HFONT font, const TextFieldStyle& style, const std::vector<TextRangeStyle>& ranges) {
    if (!dc || _bounds.right <= _bounds.left || _bounds.bottom <= _bounds.top) { return; }

    if (style.drawChrome) {
        COLORREF border = focused() ? CLR_ACCENT : style.border;
        COLORREF fill = _enabled ? style.fill : blend(style.fill, CLR_PANEL, 96);
        fillRound(dc, _bounds, fill, border, style.radius);
        if (focused()) {
            RECT ring = _bounds;
            InflateRect(&ring, 1, 1);
            drawFocusRing(dc, ring, style.radius + 1);
        }
    }

    _inner = _bounds;
    _inner.left += style.horizontalPadding;
    _inner.right -= style.horizontalPadding;
    if (_inner.right <= _inner.left) return;

    std::wstring value = text();
    updateAdvances(dc, font, value, _inner.right - _inner.left, ranges);

    int saved = SaveDC(dc);
    IntersectClipRect(dc, _inner.left, _inner.top, _inner.right, _inner.bottom);
    SetBkMode(dc, TRANSPARENT);
    HGDIOBJ oldFont = font ? SelectObject(dc, font) : nullptr;

    TEXTMETRICW metrics{};
    GetTextMetricsW(dc, &metrics);
    int lineHeight = std::max(1, static_cast<int>(_inner.bottom - _inner.top));
    int textHeight = std::max(1, static_cast<int>(metrics.tmHeight));
    int textY = _inner.top + (lineHeight - textHeight) / 2;

    if (value.empty()) {
        if (!_placeholder.empty()) {
            SetTextColor(dc, style.placeholder);
            ExtTextOutW(dc, _inner.left, textY, 0, nullptr, _placeholder.c_str(), static_cast<UINT>(_placeholder.size()), nullptr);
        }
    } else {
        for (const TextRangeStyle& range : ranges) {
            size_t first = std::min(range.first, value.size());
            size_t last = std::min(range.last, value.size());
            if (last <= first || range.background == CLR_INVALID) continue;
            RECT highlight{xForIndex(first) - std::max(0, range.paddingBefore), textY - 1, xForIndex(last), textY + textHeight + 1};
            fillRound(dc, highlight, range.background, range.border == CLR_INVALID ? range.background : range.border, 4);
        }

        auto selected = selection();
        size_t selectionFirst = std::min(selected.first, value.size());
        size_t selectionLast = std::min(selected.second, value.size());
        if (selectionLast > selectionFirst) {
            RECT selectionRect{xForIndex(selectionFirst), textY - 1, xForIndex(selectionLast), textY + textHeight + 1};
            fillRect(dc, selectionRect, focused() ? style.selection : blend(style.selection, style.fill, 105));
        }

        std::vector<size_t> boundaries{0u, value.size(), selectionFirst, selectionLast};
        for (const TextRangeStyle& range : ranges) {
            boundaries.push_back(std::min(range.first, value.size()));
            boundaries.push_back(std::min(range.last, value.size()));
        }
        std::sort(boundaries.begin(), boundaries.end());
        boundaries.erase(std::unique(boundaries.begin(), boundaries.end()), boundaries.end());

        for (size_t i = 0; i + 1u < boundaries.size(); ++i) {
            size_t first = boundaries[i];
            size_t last = boundaries[i + 1u];
            if (last <= first) continue;
            COLORREF color = _enabled ? style.text : style.disabledText;
            for (const TextRangeStyle& range : ranges) {
                if (first >= range.first && first < range.last) { color = _enabled ? range.text : style.disabledText; }
            }
            if (first >= selectionFirst && first < selectionLast) color = style.selectionText;
            SetTextColor(dc, color);
            ExtTextOutW(dc, xForIndex(first), textY, 0, nullptr, value.data() + first, static_cast<UINT>(last - first), nullptr);
        }
    }

    if (focused() && _enabled) {
        auto selected = selection();
        if (selected.first == selected.second) {
            DWORD blinkTime = GetCaretBlinkTime();
            if (blinkTime == 0 || blinkTime == INFINITE) blinkTime = 530;
            DWORD elapsed = GetTickCount() - _blinkReset;
            if ((elapsed / blinkTime) % 2u == 0u) {
                size_t caret = std::min(_activeCaret, value.size());
                int x = xForIndex(caret);
                drawLine(dc, x, textY - 1, x, textY + textHeight + 1, style.caret);
            }
        }
        size_t caret = std::min(_activeCaret, value.size());
        POINT anchor{xForIndex(caret), textY};
        _service.setInputAnchor(anchor, textHeight + 2);
    }

    if (oldFont) SelectObject(dc, oldFont);
    RestoreDC(dc, saved);
}

size_t TextField::indexAtPoint(POINT point) const {
    if (_advances.empty()) return 0;
    int target = point.x - _inner.left + _scrollX;
    if (target <= 0) return 0;
    auto found = std::lower_bound(_advances.begin(), _advances.end(), target);
    if (found == _advances.end()) return _advances.size() - 1u;
    size_t index = static_cast<size_t>(found - _advances.begin());
    if (index > 0u) {
        int previous = _advances[index - 1u];
        int next = _advances[index];
        if (target - previous < next - target) --index;
    }
    return index;
}

RECT TextField::textRangeBounds(size_t first, size_t last, int leadingPadding) const {
    if (_advances.empty()) return {};
    size_t maximum = _advances.size() - 1u;
    first = std::min(first, maximum);
    last = std::min(last, maximum);
    if (last < first) std::swap(first, last);
    return {xForIndex(first) - std::max(0, leadingPadding), _inner.top, xForIndex(last), _inner.bottom};
}

bool TextField::mouseDown(POINT point, bool extendSelection) {
    if (!_enabled || !contains(_bounds, point)) return false;
    size_t position = indexAtPoint(point);
    auto selected = selection();
    if (extendSelection) {
        if (selected.first == selected.second)
            _dragAnchor = selected.first;
        else
            _dragAnchor = _activeCaret == selected.first ? selected.second : selected.first;
    } else {
        _dragAnchor = position;
    }
    _activeCaret = position;
    _caretPreference = position < _dragAnchor ? -1 : (position > _dragAnchor ? 1 : 0);
    _dragging = true;
    resetBlink();
    std::weak_ptr<int> lifetime = _lifetime;
    _service.focus();
    if (lifetime.expired()) return true;
    _service.setSelection(_dragAnchor, position);
    return true;
}

bool TextField::mouseMove(POINT point) {
    if (!_dragging || !_enabled) return false;
    size_t position = indexAtPoint(point);
    if (position == _activeCaret) return false;
    _activeCaret = position;
    _caretPreference = position < _dragAnchor ? -1 : (position > _dragAnchor ? 1 : 0);
    resetBlink();
    _service.setSelection(_dragAnchor, position);
    return true;
}

void TextField::mouseUp() {
    _dragging = false;
}

void TextField::selectWordAt(size_t position) {
    std::wstring value = text();
    position = std::min(position, value.size());
    auto word = [](wchar_t ch) { return std::iswalnum(ch) || ch == L'_'; };
    size_t probe = position;
    if (probe == value.size() || (probe < value.size() && !word(value[probe]))) {
        if (probe == 0 || !word(value[probe - 1u])) return;
        --probe;
    }
    size_t first = probe;
    size_t last = probe + 1u;
    while (first > 0u && word(value[first - 1u])) --first;
    while (last < value.size() && word(value[last])) ++last;
    setSelection(first, last);
}

void TextField::onTextServiceEvent(TextService::Event event) {
    resetBlink();
    if (event == TextService::Event::Changed) {
        auto selected = _service.selection();
        if (selected.first == selected.second)
            _activeCaret = selected.second;
        else if (_caretPreference < 0)
            _activeCaret = selected.first;
        else if (_caretPreference > 0)
            _activeCaret = selected.second;
        else if (_activeCaret != selected.first && _activeCaret != selected.second)
            _activeCaret = selected.second;
        _caretPreference = 0;
        notify(Event::Changed);
    } else if (event == TextService::Event::CaretMoved) {
        auto selected = _service.selection();
        if (selected.first == selected.second)
            _activeCaret = selected.second;
        else if (_caretPreference < 0)
            _activeCaret = selected.first;
        else if (_caretPreference > 0)
            _activeCaret = selected.second;
        else if (_activeCaret != selected.first && _activeCaret != selected.second)
            _activeCaret = selected.second;
        _caretPreference = 0;
        notify(Event::CaretMoved);
    } else if (event == TextService::Event::FocusGained) {
        notify(Event::FocusGained);
    } else if (event == TextService::Event::FocusLost) {
        _dragging = false;
        notify(Event::FocusLost);
    }
}

void TextField::notify(Event event, const std::wstring& detail) {
    if (!_callback) return;
    Callback callback = _callback;
    callback(event, detail);
}

bool TextField::copySelection(bool cut) {
    auto selected = selection();
    std::wstring value = text();
    size_t first = std::min(selected.first, value.size());
    size_t last = std::min(selected.second, value.size());
    if (last <= first) return true;
    if (!writeClipboard(_parent, value.substr(first, last - first))) {
        notify(Event::ClipboardError, L"Copy failed because the clipboard is unavailable.");
        return false;
    }
    if (cut) {
        setSelection(first, last);
        replaceSelection(L"", true);
    }
    return true;
}

bool TextField::pasteSelection() {
    std::wstring clipboard;
    if (!readClipboard(_parent, clipboard)) {
        notify(Event::ClipboardError, L"Paste failed because the clipboard has no readable text.");
        return false;
    }
    if (!replaceSelection(clipboard, true)) {
        notify(Event::ClipboardError, L"Paste failed because the text is too long.");
        return false;
    }
    return true;
}

void TextField::resetBlink() {
    _blinkReset = GetTickCount();
}

LRESULT CALLBACK TextField::editProc(HWND window, UINT message, WPARAM wp, LPARAM lp) {
    TextField* self = static_cast<TextField*>(GetPropW(window, TEXT_FIELD_PROPERTY));
    if (!self) return DefWindowProcW(window, message, wp, lp);
    return self->handleEditMessage(window, message, wp, lp);
}

LRESULT TextField::handleEditMessage(HWND window, UINT message, WPARAM wp, LPARAM lp) {
    WNDPROC original = _editOriginal;
    if (!original) return DefWindowProcW(window, message, wp, lp);

    if (message == WM_NCDESTROY) {
        RemovePropW(window, TEXT_FIELD_PROPERTY);
        _editOriginal = nullptr;
        return CallWindowProcW(original, window, message, wp, lp);
    }
    if (message == WM_COPY) {
        copySelection(false);
        return 0;
    }
    if (message == WM_CUT) {
        copySelection(true);
        return 0;
    }
    if (message == WM_PASTE) {
        pasteSelection();
        return 0;
    }
    if (message == WM_KEYDOWN || message == WM_SYSKEYDOWN) {
        bool control = (GetKeyState(VK_CONTROL) & 0x8000) != 0;
        bool shift = (GetKeyState(VK_SHIFT) & 0x8000) != 0;
        if (shift) {
            if (wp == VK_LEFT || wp == VK_HOME || wp == VK_UP)
                _caretPreference = -1;
            else if (wp == VK_RIGHT || wp == VK_END || wp == VK_DOWN)
                _caretPreference = 1;
        } else if (wp == VK_LEFT || wp == VK_RIGHT || wp == VK_HOME || wp == VK_END || wp == VK_UP || wp == VK_DOWN) {
            _caretPreference = 0;
        }
        if (wp == VK_TAB) {
            notify(shift ? Event::TabBackward : Event::TabForward);
            return 0;
        }
        if (wp == VK_RETURN) {
            notify(Event::Enter);
            return 0;
        }
        if (wp == VK_ESCAPE) {
            notify(Event::Escape);
            return 0;
        }
        if (control) {
            switch (wp) {
            case 'A':
                selectAll();
                notify(Event::CaretMoved);
                return 0;
            case 'C': copySelection(false); return 0;
            case 'X': copySelection(true); return 0;
            case 'V': pasteSelection(); return 0;
            case 'Z': undo(); return 0;
            default: break;
            }
        }
    }
    if (message == WM_CHAR && (wp == 1 || wp == 3 || wp == 22 || wp == 24 || wp == 26 || wp == L'\t' || wp == L'\r')) { return 0; }
    return CallWindowProcW(original, window, message, wp, lp);
}

void Dropdown::setBounds(RECT bounds) {
    _bounds = bounds;
    _itemHeight = std::max(20L, bounds.bottom - bounds.top);
}

void Dropdown::setItems(std::vector<std::wstring> items) {
    _items = std::move(items);
    if (_items.empty()) {
        _selected = -1;
        close();
    } else {
        _selected = std::clamp(_selected, 0, static_cast<int>(_items.size()) - 1);
    }
}

void Dropdown::setSelectedIndex(int index) {
    if (_items.empty()) {
        _selected = -1;
        return;
    }
    _selected = std::clamp(index, 0, static_cast<int>(_items.size()) - 1);
}

void Dropdown::close() {
    _popup.close();
    _hovered = false;
}

RECT Dropdown::popupBounds(RECT clip) const {
    if (!_popup.open || _items.empty()) return {};
    int rows = std::min(_visibleRows, static_cast<int>(_items.size()));
    int height = rows * _itemHeight + 2;
    RECT result{_bounds.left, _bounds.bottom + 3, _bounds.right, _bounds.bottom + 3 + height};
    if (result.bottom > clip.bottom && _bounds.top - 3 - height >= clip.top) {
        result.bottom = _bounds.top - 3;
        result.top = result.bottom - height;
    }
    if (result.bottom > clip.bottom) result.bottom = clip.bottom;
    if (result.top < clip.top) result.top = clip.top;
    return result;
}

RECT Dropdown::itemBounds(int index, RECT clip) const {
    if (!_popup.open || index < 0 || index >= static_cast<int>(_items.size())) { return {}; }
    RECT popup = popupBounds(clip);
    RECT item{popup.left + 1, popup.top + 1 + index * _itemHeight, popup.right - 1, popup.top + 1 + (index + 1) * _itemHeight};
    return intersected(item, popup);
}

void Dropdown::openPopup(RECT clip) {
    if (_items.empty()) return;
    _popup.show(static_cast<int>(_items.size()), std::max(0, _selected));
    RECT popup = popupBounds(clip);
    int rows = std::max(1L, (popup.bottom - popup.top - 2) / _itemHeight);
    _visibleRows = std::min(rows, static_cast<int>(_items.size()));
}

int Dropdown::itemAtPoint(POINT point, RECT clip) const {
    RECT popup = popupBounds(clip);
    if (!contains(popup, point)) return -1;
    int index = (point.y - popup.top - 1) / std::max(1, _itemHeight);
    return index >= 0 && index < static_cast<int>(_items.size()) ? index : -1;
}

bool Dropdown::accept(int index) {
    if (index < 0 || index >= static_cast<int>(_items.size())) {
        close();
        return false;
    }
    bool changed = index != _selected;
    _selected = index;
    close();
    if (changed && _callback) _callback(index);
    return changed;
}

void Dropdown::draw(HDC dc, HFONT font, const DropdownStyle& style) const {
    COLORREF fill = _hovered ? style.hoverFill : style.fill;
    fillRound(dc, _bounds, fill, _focused ? CLR_ACCENT : style.border, style.radius);
    RECT text = _bounds;
    text.left += style.horizontalPadding;
    text.right -= 24;
    std::wstring label = _selected >= 0 && _selected < static_cast<int>(_items.size()) ? _items[static_cast<size_t>(_selected)] : L"";
    drawText(dc, text, label, style.text, font, DT_LEFT | DT_VCENTER | DT_SINGLELINE | DT_END_ELLIPSIS);

    int centerX = _bounds.right - 12;
    int centerY = (_bounds.top + _bounds.bottom) / 2;
    HPEN pen = CreatePen(PS_SOLID, 1, style.dimText);
    if (pen) {
        HGDIOBJ old = SelectObject(dc, pen);
        MoveToEx(dc, centerX - 4, centerY - 2, nullptr);
        LineTo(dc, centerX, centerY + 2);
        LineTo(dc, centerX + 4, centerY - 2);
        SelectObject(dc, old);
        DeleteObject(pen);
    }
    if (_focused) {
        RECT ring = _bounds;
        InflateRect(&ring, 1, 1);
        drawFocusRing(dc, ring, style.radius + 1);
    }
}

void Dropdown::drawPopup(HDC dc, HFONT font, const DropdownStyle& style, RECT clip) const {
    if (!_popup.open) return;
    RECT popup = popupBounds(clip);
    if (popup.right <= popup.left || popup.bottom <= popup.top) return;

    int saved = SaveDC(dc);
    IntersectClipRect(dc, clip.left, clip.top, clip.right, clip.bottom);
    fillRound(dc, popup, style.popupFill, style.border, style.radius);
    int rows = std::min(static_cast<int>(_items.size()), std::max(0, static_cast<int>((popup.bottom - popup.top - 2) / _itemHeight)));
    for (int i = 0; i < rows; ++i) {
        RECT item{popup.left + 1, popup.top + 1 + i * _itemHeight, popup.right - 1, popup.top + 1 + (i + 1) * _itemHeight};
        bool selected = i == _popup.selected;
        bool hovered = i == _popup.hovered;
        if (selected || hovered) {
            COLORREF fill = selected ? style.selectedFill : style.hoverFill;
            fillRound(dc, item, fill, fill, std::max(2, style.radius - 2));
        }
        RECT text = item;
        text.left += style.horizontalPadding;
        text.right -= style.horizontalPadding;
        drawText(dc, text, _items[static_cast<size_t>(i)], style.text, font, DT_LEFT | DT_VCENTER | DT_SINGLELINE | DT_END_ELLIPSIS);
    }
    RestoreDC(dc, saved);
}

bool Dropdown::mouseMove(POINT point, RECT clip) {
    bool changed = false;
    bool hovered = contains(_bounds, point);
    if (hovered != _hovered) {
        _hovered = hovered;
        changed = true;
    }
    if (_popup.open) { changed = _popup.setHovered(itemAtPoint(point, clip)) || changed; }
    return changed;
}

Dropdown::MouseResult Dropdown::mouseDown(POINT point, RECT clip) {
    if (_popup.open) {
        int index = itemAtPoint(point, clip);
        if (index >= 0) { return accept(index) ? MouseResult::SelectionChanged : MouseResult::Consumed; }
        if (contains(_bounds, point)) {
            close();
            return MouseResult::Consumed;
        }
        close();
        return MouseResult::Ignored;
    }
    if (!contains(_bounds, point)) return MouseResult::Ignored;
    openPopup(clip);
    return MouseResult::Consumed;
}

bool Dropdown::keyDown(UINT virtualKey, bool altDown, RECT clip) {
    if ((virtualKey == VK_DOWN && altDown) || virtualKey == VK_F4 || virtualKey == VK_SPACE) {
        if (_popup.open)
            close();
        else
            openPopup(clip);
        return true;
    }
    if (virtualKey == VK_ESCAPE && _popup.open) {
        close();
        return true;
    }
    if (virtualKey == VK_RETURN) {
        if (_popup.open) {
            int accepted = _popup.accept();
            accept(accepted);
        } else {
            openPopup(clip);
        }
        return true;
    }
    if (virtualKey == VK_UP || virtualKey == VK_DOWN) {
        int delta = virtualKey == VK_UP ? -1 : 1;
        if (!_popup.open) {
            openPopup(clip);
        } else {
            _popup.move(delta);
        }
        return true;
    }
    return false;
}

bool Scrollbar::configure(int contentExtent, int viewportExtent) {
    bool changed = _state.configure(contentExtent, viewportExtent);
    if (!visible()) {
        _dragging = false;
        _hovered = false;
    }
    return changed;
}

bool Scrollbar::setPosition(int position) {
    if (!_state.setPosition(position)) return false;
    changed();
    return true;
}

bool Scrollbar::scrollBy(int delta) {
    if (!_state.scrollBy(delta)) return false;
    changed();
    return true;
}

bool Scrollbar::wheel(int wheelDelta, int step) {
    _wheelRemainder += wheelDelta;
    int notches = _wheelRemainder / WHEEL_DELTA;
    _wheelRemainder %= WHEEL_DELTA;
    if (notches == 0) return false;
    return scrollBy(-notches * std::max(1, step));
}

bool Scrollbar::contains(POINT point) const {
    return visible() && ::contains(_bounds, point);
}

RECT Scrollbar::thumbBounds() const {
    if (!visible()) return {};
    int trackHeight = std::max(1L, _bounds.bottom - _bounds.top - 4);
    int content = std::max(1, _state.contentExtent());
    int viewport = std::max(1, _state.viewportExtent());
    int minimumThumb = std::max(8, MulDiv(std::max(1L, _bounds.right - _bounds.left), 24, 10));
    int thumbHeight = std::clamp(MulDiv(trackHeight, viewport, content), std::min(minimumThumb, trackHeight), trackHeight);
    int travel = std::max(0, trackHeight - thumbHeight);
    int top = _bounds.top + 2;
    if (_state.maximumPosition() > 0) top += MulDiv(travel, _state.position(), _state.maximumPosition());
    return {_bounds.left + 2, top, _bounds.right - 2, top + thumbHeight};
}

void Scrollbar::draw(HDC dc, const ScrollbarStyle& style) const {
    if (!visible()) return;
    fillRect(dc, _bounds, style.rail);
    drawLine(dc, _bounds.left, _bounds.top, _bounds.left, _bounds.bottom, style.border);
    RECT thumb = thumbBounds();
    COLORREF color = (_hovered || _dragging) ? style.thumbHover : style.thumb;
    fillRound(dc, thumb, color, color, std::max(2L, (thumb.right - thumb.left) / 2));
}

bool Scrollbar::mouseDown(POINT point) {
    if (!contains(point)) return false;
    RECT thumb = thumbBounds();
    if (::contains(thumb, point)) {
        _dragging = true;
        _dragOffset = point.y - thumb.top;
        return true;
    }
    int page = std::max(1, _state.viewportExtent() - 24);
    scrollBy(point.y < thumb.top ? -page : page);
    return true;
}

bool Scrollbar::mouseMove(POINT point) {
    bool hovered = contains(point);
    bool changedHover = hovered != _hovered;
    _hovered = hovered;
    if (!_dragging) return changedHover;

    RECT thumb = thumbBounds();
    int trackTop = _bounds.top + 2;
    int trackHeight = std::max(1, static_cast<int>(_bounds.bottom - _bounds.top - 4));
    int thumbHeight = std::max(1, static_cast<int>(thumb.bottom - thumb.top));
    int travel = std::max(1, trackHeight - thumbHeight);
    int top = std::clamp(static_cast<int>(point.y - _dragOffset), trackTop, trackTop + travel);
    int position = MulDiv(top - trackTop, _state.maximumPosition(), travel);
    return setPosition(position) || changedHover;
}

bool Scrollbar::mouseUp() {
    bool wasDragging = _dragging;
    _dragging = false;
    return wasDragging;
}

void Scrollbar::changed() {
    if (_callback) _callback(_state.position());
}

} // namespace ui
