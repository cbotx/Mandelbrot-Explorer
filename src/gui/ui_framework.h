#ifndef MANDEL_UI_FRAMEWORK_H
#define MANDEL_UI_FRAMEWORK_H

#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>

#include <functional>
#include <string>
#include <utility>
#include <vector>

namespace ui {

class DpiScale {
  public:
    explicit DpiScale(int dpi = 96) : _dpi(dpi > 0 ? dpi : 96) {}
    int dpi() const { return _dpi; }
    void setDpi(int dpi) { _dpi = dpi > 0 ? dpi : 96; }
    int px(int dip) const { return MulDiv(dip, _dpi, 96); }
    int dip(int pixels) const { return MulDiv(pixels, 96, _dpi); }
    RECT px(RECT dipRect) const;

  private:
    int _dpi;
};

class Resources {
  public:
    Resources() = default;
    ~Resources();
    Resources(const Resources&) = delete;
    Resources& operator=(const Resources&) = delete;

    bool create(int dpi);
    void reset();
    int dpi() const { return _scale.dpi(); }
    const DpiScale& scale() const { return _scale; }

    HFONT regular() const { return _regular; }
    HFONT semibold() const { return _semibold; }
    HFONT small() const { return _small; }
    HFONT mono() const { return _mono; }
    HBRUSH panelBrush() const { return _panelBrush; }
    HBRUSH cardBrush() const { return _cardBrush; }

  private:
    DpiScale _scale;
    HFONT _regular = nullptr;
    HFONT _semibold = nullptr;
    HFONT _small = nullptr;
    HFONT _mono = nullptr;
    bool _regularOwned = false;
    bool _semiboldOwned = false;
    bool _smallOwned = false;
    bool _monoOwned = false;
    HBRUSH _panelBrush = nullptr;
    HBRUSH _cardBrush = nullptr;
};

class BackBuffer {
  public:
    BackBuffer() = default;
    ~BackBuffer();
    BackBuffer(const BackBuffer&) = delete;
    BackBuffer& operator=(const BackBuffer&) = delete;

    HDC begin(HDC target, int width, int height);
    void present(HDC target, const RECT& area) const;
    void reset();
    int width() const { return _width; }
    int height() const { return _height; }

  private:
    HDC _dc = nullptr;
    HBITMAP _bitmap = nullptr;
    HGDIOBJ _oldBitmap = nullptr;
    int _width = 0;
    int _height = 0;
};

class ScrollState {
  public:
    bool configure(int contentExtent, int viewportExtent);
    bool setPosition(int position);
    bool scrollBy(int delta);
    bool handleCommand(int command, int trackPosition);
    void apply(HWND window, int bar = SB_VERT) const;

    int position() const { return _position; }
    int contentExtent() const { return _contentExtent; }
    int viewportExtent() const { return _viewportExtent; }
    int maximumPosition() const;

  private:
    int _contentExtent = 0;
    int _viewportExtent = 0;
    int _position = 0;
};

struct HitRegion {
    int id = 0;
    RECT bounds{};
    bool enabled = true;
};

class HitRouter {
  public:
    void clear();
    void add(int id, RECT bounds, bool enabled = true);
    int hit(int x, int y) const;
    bool move(int x, int y);
    bool press(int x, int y);
    int release(int x, int y);
    bool cancel();

    int hovered() const { return _hovered; }
    int pressed() const { return _pressed; }

  private:
    std::vector<HitRegion> _regions;
    int _hovered = 0;
    int _pressed = 0;
};

enum class ButtonStyle { Normal,
                         Accent,
                         Positive,
                         Subtle };

void drawButton(HDC dc, RECT bounds, const std::wstring& text, HFONT font, ButtonStyle style, bool hovered, bool pressed, bool enabled = true, int radius = 7);
void drawCard(HDC dc, RECT bounds, int radius = 8);
void drawFocusRing(HDC dc, RECT bounds, int radius = 7);

class TextService {
  public:
    enum class Event { Changed,
                       CaretMoved,
                       FocusGained,
                       FocusLost };
    using Callback = std::function<void(Event)>;

    TextService() = default;
    ~TextService();
    TextService(const TextService&) = delete;
    TextService& operator=(const TextService&) = delete;

    bool create(HWND parent, int controlId, size_t maximumLength, Callback callback);
    void destroy();
    HWND hwnd() const { return _edit; }

    bool setText(const std::wstring& text);
    std::wstring text() const;
    void setInputAnchor(POINT clientPoint, int lineHeight);
    void focus();
    bool focused() const;
    void setSelection(size_t first, size_t last);
    std::pair<size_t, size_t> selection() const;
    bool replaceSelection(const std::wstring& text, bool allowUndo = true);
    bool undo();
    void selectAll();

  private:
    static LRESULT CALLBACK editProc(HWND window, UINT message, WPARAM wp, LPARAM lp);
    void notify(Event event);

    HWND _edit = nullptr;
    WNDPROC _originalProc = nullptr;
    Callback _callback;
    size_t _maximumLength = 0;
    bool _syncing = false;
};

struct PopupState {
    bool open = false;
    int selected = -1;
    int hovered = -1;
    int itemCount = 0;

    void show(int count, int current);
    void close();
    bool move(int delta);
    bool setHovered(int index);
    int accept();
};

} // namespace ui

#endif
