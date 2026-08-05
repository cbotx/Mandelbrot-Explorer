#ifndef MANDEL_FORMULA_EDITOR_PANEL_H
#define MANDEL_FORMULA_EDITOR_PANEL_H

#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>

#include <functional>
#include <memory>

#include "formula_dialog.h"

struct FormulaEditorCallbacks {
    std::function<bool(const FormulaDialogConfig&)> apply;
    std::function<void()> useMandelbrot;
    std::function<void()> close;
};

class FormulaEditorPanel {
public:
    static constexpr int DESIGN_WIDTH = 570;

    FormulaEditorPanel();
    ~FormulaEditorPanel();
    FormulaEditorPanel(const FormulaEditorPanel&) = delete;
    FormulaEditorPanel& operator=(const FormulaEditorPanel&) = delete;

    bool create(HWND owner, int dpi, FormulaEditorCallbacks callbacks);
    void show(const FormulaDialogConfig& config);
    void hide();
    bool visible() const;
    void move(const RECT& bounds);
    void setDpi(int dpi);
    void setConfig(const FormulaDialogConfig& config);
    void dismissPopups();
    HWND hwnd() const;

private:
    struct Impl;
    std::unique_ptr<Impl> _impl;
};

#endif
