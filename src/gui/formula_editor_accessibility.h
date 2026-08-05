#ifndef MANDEL_FORMULA_EDITOR_ACCESSIBILITY_H
#define MANDEL_FORMULA_EDITOR_ACCESSIBILITY_H

#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include <oleacc.h>
#ifdef small
#undef small
#endif

#include <string>
#include <vector>

constexpr UINT WM_FORMULA_ACCESSIBILITY_SNAPSHOT = WM_APP + 0x5a0;
constexpr UINT WM_FORMULA_ACCESSIBILITY_ACTION = WM_APP + 0x5a1;

enum FormulaAccessibleKey : LONG {
    FORMULA_ACC_SELF = 0,
    FORMULA_ACC_EXPRESSION = 1,
    FORMULA_ACC_PRESET = 2,
    FORMULA_ACC_COPY = 3,
    FORMULA_ACC_PASTE = 4,
    FORMULA_ACC_C_PLANE = 5,
    FORMULA_ACC_Z0_PLANE = 6,
    FORMULA_ACC_BAILOUT = 7,
    FORMULA_ACC_PICKER = 8,
    FORMULA_ACC_REAL = 9,
    FORMULA_ACC_IMAGINARY = 10,
    FORMULA_ACC_SCROLLBAR = 11,
    FORMULA_ACC_RANGE_OUT = 12,
    FORMULA_ACC_RANGE_RESET = 13,
    FORMULA_ACC_RANGE_IN = 14,
    FORMULA_ACC_MANDELBROT = 15,
    FORMULA_ACC_REVERT = 16,
    FORMULA_ACC_APPLY = 17,
    FORMULA_ACC_CLOSE = 18,
    FORMULA_ACC_VARIABLE_BASE = 100,
    FORMULA_ACC_TAB_BASE = 200,
    FORMULA_ACC_FUNCTION_BASE = 300,
    FORMULA_ACC_PRESET_ITEM_BASE = 400
};

struct FormulaAccessibleItem {
    LONG key = 0;
    std::wstring name;
    std::wstring value;
    std::wstring description;
    std::wstring defaultAction;
    LONG role = ROLE_SYSTEM_CLIENT;
    LONG state = 0;
    RECT screenRect{};
};

struct FormulaAccessibilitySnapshotRequest {
    std::vector<FormulaAccessibleItem>* items = nullptr;
};

enum class FormulaAccessibilityAction {
    Focus,
    Select,
    Default
};

struct FormulaAccessibilityActionRequest {
    LONG key = 0;
    FormulaAccessibilityAction action = FormulaAccessibilityAction::Default;
    bool handled = false;
};

class FormulaEditorAccessibilityProvider;

FormulaEditorAccessibilityProvider* createFormulaEditorAccessibility(HWND hwnd);
void destroyFormulaEditorAccessibility(
    FormulaEditorAccessibilityProvider*& provider);
void detachFormulaEditorAccessibility(
    FormulaEditorAccessibilityProvider* provider);
LRESULT formulaEditorAccessibilityGetObject(
    FormulaEditorAccessibilityProvider* provider, WPARAM wp, LPARAM lp);
IAccessible* acquireFormulaEditorAccessibility(
    FormulaEditorAccessibilityProvider* provider);
void notifyFormulaEditorAccessibility(
    FormulaEditorAccessibilityProvider* provider, DWORD event, LONG key);

#endif
