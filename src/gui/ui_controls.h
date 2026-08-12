#ifndef MANDEL_UI_CONTROLS_H
#define MANDEL_UI_CONTROLS_H

#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>

#include <cstddef>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "ui_framework.h"

namespace ui {

struct TextRangeStyle {
    size_t first = 0;
    size_t last = 0;
    COLORREF text = RGB(224, 229, 238);
    COLORREF background = CLR_INVALID;
    COLORREF border = CLR_INVALID;
    int paddingBefore = 0;
    int paddingAfter = 0;
};

struct TextFieldStyle {
    COLORREF fill = RGB(17, 20, 26);
    COLORREF border = RGB(58, 64, 80);
    COLORREF text = RGB(224, 229, 238);
    COLORREF disabledText = RGB(140, 149, 166);
    COLORREF selection = RGB(58, 96, 158);
    COLORREF selectionText = RGB(255, 255, 255);
    COLORREF caret = RGB(224, 229, 238);
    COLORREF placeholder = RGB(103, 112, 130);
    int radius = 6;
    int horizontalPadding = 8;
    bool drawChrome = true;
};

class TextField {
  public:
    enum class Event { Changed,
                       CaretMoved,
                       FocusGained,
                       FocusLost,
                       TabForward,
                       TabBackward,
                       Enter,
                       Escape,
                       ClipboardError };
    using Callback = std::function<void(Event, const std::wstring&)>;

    TextField() = default;
    ~TextField();
    TextField(const TextField&) = delete;
    TextField& operator=(const TextField&) = delete;

    bool create(HWND parent, int controlId, size_t maximumLength, Callback callback);
    void destroy();

    void setBounds(RECT bounds);
    const RECT& bounds() const { return _bounds; }
    void setEnabled(bool enabled);
    bool enabled() const { return _enabled; }
    void setPlaceholder(std::wstring placeholder);

    bool setText(const std::wstring& text);
    std::wstring text() const;
    void focus();
    bool focused() const;
    void setSelection(size_t first, size_t last);
    std::pair<size_t, size_t> selection() const;
    void selectAll();
    bool replaceSelection(const std::wstring& text, bool allowUndo = true);
    bool undo();

    void draw(HDC dc, HFONT font, const TextFieldStyle& style, const std::vector<TextRangeStyle>& ranges = {});
    RECT textRangeBounds(size_t first, size_t last, int leadingPadding = 0) const;
    size_t indexAtPoint(POINT point) const;
    bool mouseDown(POINT point, bool extendSelection);
    bool mouseMove(POINT point);
    void mouseUp();
    bool dragging() const { return _dragging; }
    void selectWordAt(size_t position);

  private:
    static LRESULT CALLBACK editProc(HWND window, UINT message, WPARAM wp, LPARAM lp);
    LRESULT handleEditMessage(HWND window, UINT message, WPARAM wp, LPARAM lp);
    void onTextServiceEvent(TextService::Event event);
    void notify(Event event, const std::wstring& detail = {});
    bool copySelection(bool cut);
    bool pasteSelection();
    void updateAdvances(HDC dc, HFONT font, const std::wstring& value, int availableWidth, const std::vector<TextRangeStyle>& ranges);
    int xForIndex(size_t index) const;
    void resetBlink();

    TextService _service;
    HWND _parent = nullptr;
    WNDPROC _editOriginal = nullptr;
    Callback _callback;
    RECT _bounds{};
    RECT _inner{};
    std::wstring _placeholder;
    std::vector<int> _advances;
    size_t _maximumLength = 0;
    size_t _dragAnchor = 0;
    size_t _activeCaret = 0;
    int _scrollX = 0;
    DWORD _blinkReset = 0;
    int _caretPreference = 0;
    bool _enabled = true;
    bool _dragging = false;
    std::shared_ptr<int> _lifetime = std::make_shared<int>(0);
};

struct DropdownStyle {
    COLORREF fill = RGB(37, 41, 52);
    COLORREF hoverFill = RGB(48, 53, 66);
    COLORREF popupFill = RGB(29, 33, 42);
    COLORREF selectedFill = RGB(58, 96, 158);
    COLORREF border = RGB(58, 64, 80);
    COLORREF text = RGB(224, 229, 238);
    COLORREF dimText = RGB(140, 149, 166);
    int radius = 6;
    int horizontalPadding = 8;
};

class Dropdown {
  public:
    enum class MouseResult { Ignored,
                             Consumed,
                             SelectionChanged };
    using Callback = std::function<void(int)>;

    void setBounds(RECT bounds);
    const RECT& bounds() const { return _bounds; }
    void setItems(std::vector<std::wstring> items);
    const std::vector<std::wstring>& items() const { return _items; }
    void setSelectedIndex(int index);
    int selectedIndex() const { return _selected; }
    void setCallback(Callback callback) { _callback = std::move(callback); }
    void setFocused(bool focused) { _focused = focused; }
    bool focused() const { return _focused; }
    bool open() const { return _popup.open; }
    void close();

    void draw(HDC dc, HFONT font, const DropdownStyle& style) const;
    void drawPopup(HDC dc, HFONT font, const DropdownStyle& style, RECT clip) const;
    bool mouseMove(POINT point, RECT clip);
    MouseResult mouseDown(POINT point, RECT clip);
    bool keyDown(UINT virtualKey, bool altDown, RECT clip);
    RECT popupBounds(RECT clip) const;
    RECT itemBounds(int index, RECT clip) const;

  private:
    void openPopup(RECT clip);
    int itemAtPoint(POINT point, RECT clip) const;
    bool accept(int index);

    RECT _bounds{};
    std::vector<std::wstring> _items;
    PopupState _popup;
    Callback _callback;
    int _selected = -1;
    int _itemHeight = 26;
    int _visibleRows = 8;
    bool _focused = false;
    bool _hovered = false;
};

struct ScrollbarStyle {
    COLORREF rail = RGB(21, 24, 33);
    COLORREF thumb = RGB(75, 82, 98);
    COLORREF thumbHover = RGB(92, 101, 121);
    COLORREF border = RGB(43, 48, 59);
};

class Scrollbar {
  public:
    using Callback = std::function<void(int)>;

    void setBounds(RECT bounds) { _bounds = bounds; }
    const RECT& bounds() const { return _bounds; }
    void setCallback(Callback callback) { _callback = std::move(callback); }
    bool configure(int contentExtent, int viewportExtent);
    bool visible() const { return _state.maximumPosition() > 0; }
    int position() const { return _state.position(); }
    int maximumPosition() const { return _state.maximumPosition(); }
    bool setPosition(int position);
    bool scrollBy(int delta);
    bool wheel(int wheelDelta, int step);

    void draw(HDC dc, const ScrollbarStyle& style) const;
    bool mouseDown(POINT point);
    bool mouseMove(POINT point);
    bool mouseUp();
    bool dragging() const { return _dragging; }
    bool contains(POINT point) const;

  private:
    RECT thumbBounds() const;
    void changed();

    RECT _bounds{};
    ScrollState _state;
    Callback _callback;
    int _dragOffset = 0;
    int _wheelRemainder = 0;
    bool _hovered = false;
    bool _dragging = false;
};

} // namespace ui

#endif
