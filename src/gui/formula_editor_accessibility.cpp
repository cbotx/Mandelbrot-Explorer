#include "formula_editor_accessibility.h"

#include <algorithm>
#include <new>

namespace {

bool validChildVariant(const VARIANT& child) {
    return child.vt == VT_I4 && child.lVal >= CHILDID_SELF;
}

HRESULT copyString(const std::wstring& value, BSTR* result) {
    if (!result) return E_INVALIDARG;
    *result = SysAllocStringLen(
        value.data(), static_cast<UINT>(value.size()));
    return *result || value.empty() ? S_OK : E_OUTOFMEMORY;
}

void setEmpty(VARIANT* value) {
    if (value) VariantInit(value);
}

class SelectedChildrenEnumerator final : public IEnumVARIANT {
public:
    explicit SelectedChildrenEnumerator(std::vector<LONG> children,
                                        ULONG position = 0)
        : _children(std::move(children)), _position(position) {}

    HRESULT STDMETHODCALLTYPE QueryInterface(
        REFIID iid, void** object) override {
        if (!object) return E_INVALIDARG;
        *object = nullptr;
        if (iid == IID_IUnknown || iid == IID_IEnumVARIANT) {
            *object = static_cast<IEnumVARIANT*>(this);
            AddRef();
            return S_OK;
        }
        return E_NOINTERFACE;
    }

    ULONG STDMETHODCALLTYPE AddRef() override {
        return static_cast<ULONG>(InterlockedIncrement(&_references));
    }

    ULONG STDMETHODCALLTYPE Release() override {
        LONG references = InterlockedDecrement(&_references);
        if (references == 0) delete this;
        return static_cast<ULONG>(references);
    }

    HRESULT STDMETHODCALLTYPE Next(
        ULONG count, VARIANT* values, ULONG* fetched) override {
        if (!values || (count > 1 && !fetched)) return E_INVALIDARG;
        ULONG produced = 0;
        while (produced < count && _position < _children.size()) {
            VariantInit(&values[produced]);
            values[produced].vt = VT_I4;
            values[produced].lVal = _children[_position];
            ++produced;
            ++_position;
        }
        if (fetched) *fetched = produced;
        return produced == count ? S_OK : S_FALSE;
    }

    HRESULT STDMETHODCALLTYPE Skip(ULONG count) override {
        size_t remaining = _children.size() - _position;
        size_t skipped = (std::min)(remaining, static_cast<size_t>(count));
        _position += skipped;
        return skipped == count ? S_OK : S_FALSE;
    }

    HRESULT STDMETHODCALLTYPE Reset() override {
        _position = 0;
        return S_OK;
    }

    HRESULT STDMETHODCALLTYPE Clone(IEnumVARIANT** clone) override {
        if (!clone) return E_INVALIDARG;
        *clone = new (std::nothrow)
            SelectedChildrenEnumerator(_children, _position);
        return *clone ? S_OK : E_OUTOFMEMORY;
    }

private:
    volatile LONG _references = 1;
    std::vector<LONG> _children;
    size_t _position = 0;
};

} // namespace

class FormulaEditorAccessibilityProvider final : public IAccessible {
public:
    explicit FormulaEditorAccessibilityProvider(HWND hwnd)
        : _hwnd(hwnd), _comThread(GetCurrentThreadId()) {
        HRESULT initialized =
            CoInitializeEx(nullptr, COINIT_APARTMENTTHREADED);
        _uninitializeCom = SUCCEEDED(initialized) ? 1 : 0;
    }

    ~FormulaEditorAccessibilityProvider() {
        uninitializeCom();
    }

    void detach() {
        InterlockedExchangePointer(
            reinterpret_cast<PVOID volatile*>(&_hwnd), nullptr);
        uninitializeCom();
    }

    HWND hwnd() const {
        return static_cast<HWND>(InterlockedCompareExchangePointer(
            reinterpret_cast<PVOID volatile*>(
                const_cast<HWND*>(&_hwnd)),
            nullptr, nullptr));
    }

    bool snapshot(std::vector<FormulaAccessibleItem>& items) const {
        items.clear();
        HWND window = hwnd();
        if (!window || !IsWindow(window)) return false;
        FormulaAccessibilitySnapshotRequest request{&items};
        return SendMessageW(
                   window, WM_FORMULA_ACCESSIBILITY_SNAPSHOT,
                   reinterpret_cast<WPARAM>(this),
                   reinterpret_cast<LPARAM>(&request)) != 0;
    }

    bool action(LONG key, FormulaAccessibilityAction actionKind) const {
        HWND window = hwnd();
        if (!window || !IsWindow(window)) return false;
        FormulaAccessibilityActionRequest request{
            key, actionKind, false};
        SendMessageW(
            window, WM_FORMULA_ACCESSIBILITY_ACTION,
            reinterpret_cast<WPARAM>(this),
            reinterpret_cast<LPARAM>(&request));
        return request.handled;
    }

    void notify(DWORD event, LONG key) const {
        HWND window = hwnd();
        if (!window || !IsWindow(window)) return;
        LONG childId = CHILDID_SELF;
        if (key != FORMULA_ACC_SELF) {
            std::vector<FormulaAccessibleItem> items;
            if (!snapshot(items)) return;
            auto found = std::find_if(
                items.begin(), items.end(),
                [key](const FormulaAccessibleItem& item) {
                    return item.key == key;
                });
            if (found == items.end()) return;
            childId = static_cast<LONG>(
                std::distance(items.begin(), found) + 1);
        }
        NotifyWinEvent(event, window, OBJID_CLIENT, childId);
    }

    HRESULT STDMETHODCALLTYPE QueryInterface(
        REFIID iid, void** object) override {
        if (!object) return E_INVALIDARG;
        *object = nullptr;
        if (iid == IID_IUnknown || iid == IID_IDispatch ||
            iid == IID_IAccessible) {
            *object = static_cast<IAccessible*>(this);
            AddRef();
            return S_OK;
        }
        return E_NOINTERFACE;
    }

    ULONG STDMETHODCALLTYPE AddRef() override {
        return static_cast<ULONG>(InterlockedIncrement(&_references));
    }

    ULONG STDMETHODCALLTYPE Release() override {
        LONG references = InterlockedDecrement(&_references);
        if (references == 0) delete this;
        return static_cast<ULONG>(references);
    }

    HRESULT STDMETHODCALLTYPE GetTypeInfoCount(UINT* count) override {
        if (!count) return E_INVALIDARG;
        *count = 0;
        return S_OK;
    }

    HRESULT STDMETHODCALLTYPE GetTypeInfo(
        UINT, LCID, ITypeInfo**) override {
        return E_NOTIMPL;
    }

    HRESULT STDMETHODCALLTYPE GetIDsOfNames(
        REFIID, LPOLESTR*, UINT, LCID, DISPID*) override {
        return DISP_E_UNKNOWNNAME;
    }

    HRESULT STDMETHODCALLTYPE Invoke(
        DISPID, REFIID, LCID, WORD, DISPPARAMS*,
        VARIANT*, EXCEPINFO*, UINT*) override {
        return DISP_E_MEMBERNOTFOUND;
    }

    HRESULT STDMETHODCALLTYPE get_accParent(IDispatch** parent) override {
        if (!parent) return E_INVALIDARG;
        *parent = nullptr;
        HWND window = hwnd();
        if (!window) return CO_E_OBJNOTCONNECTED;
        HWND parentWindow = GetParent(window);
        if (!parentWindow) return S_FALSE;
        IAccessible* accessible = nullptr;
        HRESULT result = AccessibleObjectFromWindow(
            parentWindow, OBJID_CLIENT, IID_IAccessible,
            reinterpret_cast<void**>(&accessible));
        if (FAILED(result)) return result;
        *parent = static_cast<IDispatch*>(accessible);
        return S_OK;
    }

    HRESULT STDMETHODCALLTYPE get_accChildCount(LONG* count) override {
        if (!count) return E_INVALIDARG;
        *count = 0;
        std::vector<FormulaAccessibleItem> items;
        if (!snapshot(items)) return CO_E_OBJNOTCONNECTED;
        *count = static_cast<LONG>(items.size());
        return S_OK;
    }

    HRESULT STDMETHODCALLTYPE get_accChild(
        VARIANT child, IDispatch** result) override {
        if (!result) return E_INVALIDARG;
        *result = nullptr;
        if (!validChildVariant(child)) return E_INVALIDARG;
        return S_FALSE;
    }

    HRESULT STDMETHODCALLTYPE get_accName(
        VARIANT child, BSTR* name) override {
        if (!validChildVariant(child)) return E_INVALIDARG;
        if (child.lVal == CHILDID_SELF)
            return copyString(L"Formula editor", name);
        FormulaAccessibleItem item;
        HRESULT result = resolve(child.lVal, item);
        return SUCCEEDED(result) ? copyString(item.name, name) : result;
    }

    HRESULT STDMETHODCALLTYPE get_accValue(
        VARIANT child, BSTR* value) override {
        if (!validChildVariant(child)) return E_INVALIDARG;
        if (child.lVal == CHILDID_SELF) return S_FALSE;
        FormulaAccessibleItem item;
        HRESULT result = resolve(child.lVal, item);
        if (FAILED(result)) return result;
        if (item.value.empty()) return S_FALSE;
        return copyString(item.value, value);
    }

    HRESULT STDMETHODCALLTYPE get_accDescription(
        VARIANT child, BSTR* description) override {
        if (!validChildVariant(child)) return E_INVALIDARG;
        if (child.lVal == CHILDID_SELF) {
            return copyString(
                L"Build an orbit formula and configure its complex values.",
                description);
        }
        FormulaAccessibleItem item;
        HRESULT result = resolve(child.lVal, item);
        if (FAILED(result)) return result;
        if (item.description.empty()) return S_FALSE;
        return copyString(item.description, description);
    }

    HRESULT STDMETHODCALLTYPE get_accRole(
        VARIANT child, VARIANT* role) override {
        if (!role || !validChildVariant(child)) return E_INVALIDARG;
        VariantInit(role);
        role->vt = VT_I4;
        if (child.lVal == CHILDID_SELF) {
            role->lVal = ROLE_SYSTEM_PANE;
            return S_OK;
        }
        FormulaAccessibleItem item;
        HRESULT result = resolve(child.lVal, item);
        if (FAILED(result)) return result;
        role->lVal = item.role;
        return S_OK;
    }

    HRESULT STDMETHODCALLTYPE get_accState(
        VARIANT child, VARIANT* state) override {
        if (!state || !validChildVariant(child)) return E_INVALIDARG;
        VariantInit(state);
        state->vt = VT_I4;
        if (child.lVal == CHILDID_SELF) {
            HWND window = hwnd();
            if (!window) return CO_E_OBJNOTCONNECTED;
            state->lVal = STATE_SYSTEM_FOCUSABLE;
            if (!IsWindowVisible(window))
                state->lVal |= STATE_SYSTEM_INVISIBLE;
            HWND focus = GetFocus();
            if (focus == window || (focus && IsChild(window, focus)))
                state->lVal |= STATE_SYSTEM_FOCUSED;
            return S_OK;
        }
        FormulaAccessibleItem item;
        HRESULT result = resolve(child.lVal, item);
        if (FAILED(result)) return result;
        state->lVal = item.state;
        return S_OK;
    }

    HRESULT STDMETHODCALLTYPE get_accHelp(
        VARIANT, BSTR* help) override {
        if (help) *help = nullptr;
        return S_FALSE;
    }

    HRESULT STDMETHODCALLTYPE get_accHelpTopic(
        BSTR* helpFile, VARIANT, LONG* topic) override {
        if (helpFile) *helpFile = nullptr;
        if (topic) *topic = -1;
        return S_FALSE;
    }

    HRESULT STDMETHODCALLTYPE get_accKeyboardShortcut(
        VARIANT, BSTR* shortcut) override {
        if (shortcut) *shortcut = nullptr;
        return S_FALSE;
    }

    HRESULT STDMETHODCALLTYPE get_accFocus(VARIANT* focus) override {
        if (!focus) return E_INVALIDARG;
        VariantInit(focus);
        std::vector<FormulaAccessibleItem> items;
        if (!snapshot(items)) return CO_E_OBJNOTCONNECTED;
        auto found = std::find_if(
            items.begin(), items.end(),
            [](const FormulaAccessibleItem& item) {
                return (item.state & STATE_SYSTEM_FOCUSED) != 0;
            });
        if (found != items.end()) {
            focus->vt = VT_I4;
            focus->lVal = static_cast<LONG>(
                std::distance(items.begin(), found) + 1);
            return S_OK;
        }
        HWND window = hwnd();
        HWND keyboardFocus = GetFocus();
        if (window &&
            (keyboardFocus == window ||
             (keyboardFocus && IsChild(window, keyboardFocus)))) {
            focus->vt = VT_I4;
            focus->lVal = CHILDID_SELF;
            return S_OK;
        }
        return S_FALSE;
    }

    HRESULT STDMETHODCALLTYPE get_accSelection(
        VARIANT* selection) override {
        if (!selection) return E_INVALIDARG;
        VariantInit(selection);
        std::vector<FormulaAccessibleItem> items;
        if (!snapshot(items)) return CO_E_OBJNOTCONNECTED;
        std::vector<LONG> selected;
        for (size_t i = 0; i < items.size(); ++i) {
            if ((items[i].state & STATE_SYSTEM_SELECTED) != 0)
                selected.push_back(static_cast<LONG>(i + 1));
        }
        if (selected.empty()) return S_FALSE;
        if (selected.size() == 1) {
            selection->vt = VT_I4;
            selection->lVal = selected.front();
            return S_OK;
        }
        auto* enumerator = new (std::nothrow)
            SelectedChildrenEnumerator(std::move(selected));
        if (!enumerator) return E_OUTOFMEMORY;
        selection->vt = VT_UNKNOWN;
        selection->punkVal = enumerator;
        return S_OK;
    }

    HRESULT STDMETHODCALLTYPE get_accDefaultAction(
        VARIANT child, BSTR* actionName) override {
        if (!validChildVariant(child)) return E_INVALIDARG;
        if (child.lVal == CHILDID_SELF) return S_FALSE;
        FormulaAccessibleItem item;
        HRESULT result = resolve(child.lVal, item);
        if (FAILED(result)) return result;
        if (item.defaultAction.empty()) return S_FALSE;
        return copyString(item.defaultAction, actionName);
    }

    HRESULT STDMETHODCALLTYPE accSelect(
        LONG flags, VARIANT child) override {
        if (!validChildVariant(child)) return E_INVALIDARG;
        const LONG supported =
            SELFLAG_TAKEFOCUS | SELFLAG_TAKESELECTION;
        if ((flags & supported) == 0 ||
            (flags & ~supported) != 0) {
            return E_INVALIDARG;
        }
        if (child.lVal == CHILDID_SELF) {
            return action(FORMULA_ACC_SELF,
                          FormulaAccessibilityAction::Focus)
                ? S_OK : E_FAIL;
        }
        FormulaAccessibleItem item;
        HRESULT result = resolve(child.lVal, item);
        if (FAILED(result)) return result;
        FormulaAccessibilityAction operation =
            (flags & SELFLAG_TAKESELECTION)
                ? FormulaAccessibilityAction::Select
                : FormulaAccessibilityAction::Focus;
        return action(item.key, operation) ? S_OK : E_FAIL;
    }

    HRESULT STDMETHODCALLTYPE accLocation(
        LONG* left, LONG* top, LONG* width, LONG* height,
        VARIANT child) override {
        if (!left || !top || !width || !height ||
            !validChildVariant(child)) {
            return E_INVALIDARG;
        }
        RECT rect{};
        if (child.lVal == CHILDID_SELF) {
            HWND window = hwnd();
            if (!window) return CO_E_OBJNOTCONNECTED;
            GetWindowRect(window, &rect);
        } else {
            FormulaAccessibleItem item;
            HRESULT result = resolve(child.lVal, item);
            if (FAILED(result)) return result;
            rect = item.screenRect;
        }
        *left = rect.left;
        *top = rect.top;
        *width = std::max(0L, rect.right - rect.left);
        *height = std::max(0L, rect.bottom - rect.top);
        return S_OK;
    }

    HRESULT STDMETHODCALLTYPE accNavigate(
        LONG direction, VARIANT start, VARIANT* destination) override {
        if (!destination || !validChildVariant(start))
            return E_INVALIDARG;
        VariantInit(destination);
        std::vector<FormulaAccessibleItem> items;
        if (!snapshot(items)) return CO_E_OBJNOTCONNECTED;
        LONG count = static_cast<LONG>(items.size());
        LONG target = 0;
        auto navigable = [&items](LONG childId) {
            return childId >= 1 &&
                   childId <= static_cast<LONG>(items.size()) &&
                   (items[static_cast<size_t>(childId - 1)].state &
                    STATE_SYSTEM_INVISIBLE) == 0;
        };
        if (start.lVal == CHILDID_SELF) {
            if (direction == NAVDIR_FIRSTCHILD) {
                for (LONG candidate = 1; candidate <= count; ++candidate) {
                    if (navigable(candidate)) {
                        target = candidate;
                        break;
                    }
                }
            } else if (direction == NAVDIR_LASTCHILD) {
                for (LONG candidate = count; candidate >= 1; --candidate) {
                    if (navigable(candidate)) {
                        target = candidate;
                        break;
                    }
                }
            }
        } else if (start.lVal >= 1 && start.lVal <= count) {
            if (direction == NAVDIR_NEXT || direction == NAVDIR_DOWN) {
                for (LONG candidate = start.lVal + 1;
                     candidate <= count; ++candidate) {
                    if (navigable(candidate)) {
                        target = candidate;
                        break;
                    }
                }
            } else if (direction == NAVDIR_PREVIOUS ||
                       direction == NAVDIR_UP) {
                for (LONG candidate = start.lVal - 1;
                     candidate >= 1; --candidate) {
                    if (navigable(candidate)) {
                        target = candidate;
                        break;
                    }
                }
            }
        }
        if (target == 0) return S_FALSE;
        destination->vt = VT_I4;
        destination->lVal = target;
        return S_OK;
    }

    HRESULT STDMETHODCALLTYPE accHitTest(
        LONG x, LONG y, VARIANT* child) override {
        if (!child) return E_INVALIDARG;
        VariantInit(child);
        std::vector<FormulaAccessibleItem> items;
        if (!snapshot(items)) return CO_E_OBJNOTCONNECTED;
        for (size_t i = items.size(); i-- > 0;) {
            const FormulaAccessibleItem& item = items[i];
            if ((item.state & (STATE_SYSTEM_INVISIBLE |
                               STATE_SYSTEM_OFFSCREEN)) != 0) {
                continue;
            }
            if (x >= item.screenRect.left && x < item.screenRect.right &&
                y >= item.screenRect.top && y < item.screenRect.bottom) {
                child->vt = VT_I4;
                child->lVal = static_cast<LONG>(i + 1);
                return S_OK;
            }
        }
        HWND window = hwnd();
        RECT rect{};
        if (window && GetWindowRect(window, &rect) &&
            x >= rect.left && x < rect.right &&
            y >= rect.top && y < rect.bottom) {
            child->vt = VT_I4;
            child->lVal = CHILDID_SELF;
            return S_OK;
        }
        return S_FALSE;
    }

    HRESULT STDMETHODCALLTYPE accDoDefaultAction(
        VARIANT child) override {
        if (!validChildVariant(child) ||
            child.lVal == CHILDID_SELF) {
            return E_INVALIDARG;
        }
        FormulaAccessibleItem item;
        HRESULT result = resolve(child.lVal, item);
        if (FAILED(result)) return result;
        if ((item.state & (STATE_SYSTEM_UNAVAILABLE |
                           STATE_SYSTEM_INVISIBLE)) != 0) {
            return E_FAIL;
        }
        return action(item.key, FormulaAccessibilityAction::Default)
            ? S_OK : E_FAIL;
    }

    HRESULT STDMETHODCALLTYPE put_accName(
        VARIANT, BSTR) override {
        return E_ACCESSDENIED;
    }

    HRESULT STDMETHODCALLTYPE put_accValue(
        VARIANT, BSTR) override {
        return E_ACCESSDENIED;
    }

private:
    void uninitializeCom() {
        if (GetCurrentThreadId() != _comThread) return;
        if (InterlockedExchange(&_uninitializeCom, 0) != 0)
            CoUninitialize();
    }

    HRESULT resolve(LONG childId, FormulaAccessibleItem& item) const {
        std::vector<FormulaAccessibleItem> items;
        if (!snapshot(items)) return CO_E_OBJNOTCONNECTED;
        if (childId < 1 ||
            childId > static_cast<LONG>(items.size())) {
            return E_INVALIDARG;
        }
        item = std::move(items[static_cast<size_t>(childId - 1)]);
        return S_OK;
    }

    volatile LONG _references = 1;
    HWND volatile _hwnd = nullptr;
    DWORD _comThread = 0;
    volatile LONG _uninitializeCom = 0;
};

FormulaEditorAccessibilityProvider* createFormulaEditorAccessibility(
    HWND hwnd) {
    if (!hwnd) return nullptr;
    return new (std::nothrow) FormulaEditorAccessibilityProvider(hwnd);
}

void destroyFormulaEditorAccessibility(
    FormulaEditorAccessibilityProvider*& provider) {
    if (!provider) return;
    provider->detach();
    provider->Release();
    provider = nullptr;
}

void detachFormulaEditorAccessibility(
    FormulaEditorAccessibilityProvider* provider) {
    if (provider) provider->detach();
}

LRESULT formulaEditorAccessibilityGetObject(
    FormulaEditorAccessibilityProvider* provider, WPARAM wp, LPARAM lp) {
    if (!provider || static_cast<LONG>(lp) != OBJID_CLIENT ||
        !provider->hwnd()) {
        return 0;
    }
    return LresultFromObject(
        IID_IAccessible, wp, static_cast<IAccessible*>(provider));
}

IAccessible* acquireFormulaEditorAccessibility(
    FormulaEditorAccessibilityProvider* provider) {
    if (!provider) return nullptr;
    provider->AddRef();
    return static_cast<IAccessible*>(provider);
}

void notifyFormulaEditorAccessibility(
    FormulaEditorAccessibilityProvider* provider,
    DWORD event, LONG key) {
    if (provider) provider->notify(event, key);
}
