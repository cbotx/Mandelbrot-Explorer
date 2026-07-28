#ifndef GUI_THEME_H
#define GUI_THEME_H

#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include <string>

// Fixed fractal render buffer + side-panel / status-bar geometry.
constexpr int RENDER_W = 900;
constexpr int RENDER_H = 600;
constexpr int PANEL_W = 330;
constexpr int STATUS_H = 30;
constexpr UINT TIMER_ID = 1;

// Dark theme shared by the main panel and the custom-drawn export dialog.
inline const COLORREF CLR_BG        = RGB(18, 20, 26);
inline const COLORREF CLR_PANEL     = RGB(26, 29, 37);
inline const COLORREF CLR_CARD      = RGB(37, 41, 52);
inline const COLORREF CLR_CARD_HOV  = RGB(48, 53, 66);
inline const COLORREF CLR_TRACK     = RGB(52, 57, 70);
inline const COLORREF CLR_ACCENT    = RGB(94, 148, 245);
inline const COLORREF CLR_ACCENT_HI = RGB(126, 172, 250);
inline const COLORREF CLR_GREEN     = RGB(112, 200, 158);
inline const COLORREF CLR_TEXT      = RGB(224, 229, 238);
inline const COLORREF CLR_TEXT_DIM  = RGB(140, 149, 166);
inline const COLORREF CLR_BORDER    = RGB(58, 64, 80);

// Shared self-drawn widgets.
void fillRound(HDC dc, RECT r, COLORREF fill, COLORREF border, int rad);
void fillRect(HDC dc, RECT r, COLORREF c);
void drawText(HDC dc, RECT r, const std::wstring& s, COLORREF c, HFONT f, UINT fmt);

#endif
