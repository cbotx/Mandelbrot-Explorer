#include "gui_theme.h"

void fillRound(HDC dc, RECT r, COLORREF fill, COLORREF border, int rad) {
    HBRUSH b = CreateSolidBrush(fill);
    HPEN p = border == CLR_BG ? (HPEN)GetStockObject(NULL_PEN) : CreatePen(PS_SOLID, 1, border);
    HGDIOBJ ob = SelectObject(dc, b), op = SelectObject(dc, p);
    RoundRect(dc, r.left, r.top, r.right, r.bottom, rad, rad);
    SelectObject(dc, ob); SelectObject(dc, op);
    DeleteObject(b);
    if (border != CLR_BG) DeleteObject(p);
}
void fillRect(HDC dc, RECT r, COLORREF c) {
    HBRUSH b = CreateSolidBrush(c); FillRect(dc, &r, b); DeleteObject(b);
}
void drawText(HDC dc, RECT r, const std::wstring& s, COLORREF c, HFONT f, UINT fmt) {
    HGDIOBJ of = SelectObject(dc, f);
    SetTextColor(dc, c); SetBkMode(dc, TRANSPARENT);
    DrawTextW(dc, s.c_str(), -1, &r, fmt);
    SelectObject(dc, of);
}
