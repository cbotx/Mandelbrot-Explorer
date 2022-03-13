#ifndef __WIDGET_H__
#define __WIDGET_H__

#include <vector>
#include <array>

#include "glut_canvas.h"
#include "observer.h"
#include "Image.h"
#include "mandel_navigator.h"


class Widget {
protected:
    GlutCanvas* _canvas;
    int _w;
    int _h;
    int _x;
    int _y;
    const int MIN_W = 16;
    bool _require_redraw = true;

public:
    explicit Widget(GlutCanvas* canvas, int width, int height, int x, int y);

    virtual ~Widget() = default;

    bool RequireRedraw();

    void SetRedraw(bool redraw);

    virtual void Update();

    virtual void Render() = 0;

    virtual void MouseDownLC(int x, int y);

    virtual void MouseUpLC(int x, int y);

    virtual void MouseDownRC(int x, int y);

    virtual void MouseUpRC(int x, int y);

    virtual void MouseMove(int x, int y);

    bool HitTest(int x, int y);


};

class Label : public Widget {
private:
    const char* _text;
public:
    explicit Label(GlutCanvas* canvas, int width, int height, int x, int y, const char* text);
    
    virtual ~Label() = default;

    void Render() override;
};

class GradientBar : public Widget, virtual public ColorObserver, virtual public ColorSubject {
private:
    GlutCanvas::Texture _texture;
    uint8_t* _color;
    std::vector<std::array<float, 4>> _p;
    float* _mp, * _mp_full;
    static constexpr float DMSIZE = 5;
    int _drag_idx = -1;
    int _selected = 0;
    MandelNavigator* _nav = nullptr;

public:
    explicit GradientBar(GlutCanvas* canvas, int width, int height, int x, int y);
    
    virtual ~GradientBar() = default;

    void SetSceneSunRise();

    void SetSceneSnowy();

    void Update() override;
    
    void Render() override;

    void MouseDownLC(int x, int y) override;

    void MouseUpLC(int x, int y) override;

    void MouseDownRC(int x, int y) override;

    void MouseUpRC(int x, int y) override;

    void MouseMove(int x, int y) override;

    void UpdateColor(float R, float G, float B) override;

    void BindNavigator(MandelNavigator* nav);

};

class ColorPickerBarV;

class ColorPickerDisk : public Widget, virtual public ColorObserver, virtual public ColorSubject {
    friend class ColorPickerBarV;
private:
    GlutCanvas::Texture _texture;
    uint8_t* _color;
    std::array<int, 2> _p;
    static const int RADIUS = 100;
    static constexpr float DMSIZE = 8;

    float _fh, _fs, _fv;
    ColorPickerBarV* _cpb;

public:
    ColorPickerDisk(GlutCanvas* canvas, int width, int height, int x, int y);

    virtual ~ColorPickerDisk() = default;

    void Bind(ColorPickerBarV* cpb);

    void Update() override;

    void Render() override;

    void MouseDownLC(int x, int y) override;

    void MouseUpLC(int x, int y) override;

    void MouseDownRC(int x, int y) override;

    void MouseUpRC(int x, int y) override;

    void MouseMove(int x, int y) override;

    void UpdateColor(float R, float G, float B) override;

};

class ColorPickerBarV : public Widget {
private:
    GlutCanvas::Texture _texture;
    uint8_t* _color;
    int _p;
    static const int RADIUS = 100;
    static constexpr float DMSIZE = 2;

    ColorPickerDisk* _cpd;

public:
    ColorPickerBarV(GlutCanvas* canvas, int width, int height, int x, int y);

    virtual ~ColorPickerBarV() = default;

    void Bind(ColorPickerDisk* cpd);

    void Update() override;

    void Render() override;

    void MouseDownLC(int x, int y) override;

    void MouseUpLC(int x, int y) override;

    void MouseDownRC(int x, int y) override;

    void MouseUpRC(int x, int y) override;

    void MouseMove(int x, int y) override;

};

class ScrollHBase : public Widget {
protected:
    float _min, _max, _val, _ratio;
    static constexpr float DMSIZE = 6;

public:
    ScrollHBase(GlutCanvas* canvas, int width, int height, int x, int y, int mn, int mx, int val);

    virtual ~ScrollHBase() = default;

    void Update() override;

    void Render() override;

    void MouseDownLC(int x, int y) override;

    void MouseUpLC(int x, int y) override;

    void MouseDownRC(int x, int y) override;

    void MouseUpRC(int x, int y) override;

    void MouseMove(int x, int y) override;

protected:
    virtual void valueOnChange() = 0;

};

class ScrollHColor : public ScrollHBase {
private:
    MandelNavigator* _nav = nullptr;

public:
    ScrollHColor(GlutCanvas* canvas, int width, int height, int x, int y, int mn, int mx, int val);

    virtual ~ScrollHColor() = default;

    void BindNavigator(MandelNavigator* nav);

private:
    void valueOnChange() override;
};

class ScrollHMxit : public ScrollHBase {
private:
    MandelNavigator* _nav = nullptr;

public:
    ScrollHMxit(GlutCanvas* canvas, int width, int height, int x, int y, int mn, int mx, int val);

    virtual ~ScrollHMxit() = default;

    void BindNavigator(MandelNavigator* nav);
    
private:
    void valueOnChange() override;

};


class CheckerBase : public Widget {
protected:
    bool _state;
    static constexpr float SIZE_RATIO = 0.8;

public:
    CheckerBase(GlutCanvas* canvas, int width, int height, int x, int y, bool state);

    virtual ~CheckerBase() = default;

    void Update() override;

    void Render() override;

    void MouseDownLC(int x, int y) override;

protected:
    virtual void valueOnChange() = 0;

};


class CheckerSuperSampling : public CheckerBase {
private:
    MandelNavigator* _nav = nullptr;

public:
    CheckerSuperSampling(GlutCanvas* canvas, int width, int height, int x, int y, bool state);

    virtual ~CheckerSuperSampling() = default;

    void BindNavigator(MandelNavigator* nav);

private:
    void valueOnChange() override;

};


class CheckerExteriorDE : public CheckerBase {
private:
    MandelNavigator* _nav = nullptr;
    GradientBar* _gb = nullptr;
    
public:
    CheckerExteriorDE(GlutCanvas* canvas, int width, int height, int x, int y, bool state);

    virtual ~CheckerExteriorDE() = default;

    void BindNavigator(MandelNavigator* nav);

    void BindGradientBar(GradientBar* gb);

private:
    void valueOnChange() override;

};

#endif