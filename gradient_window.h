#ifndef __GRADIENT_WINDOW_H__
#define __GRADIENT_WINDOW_H__


#include <vector>

#include "glut_canvas.h"
#include "widget.h"
#include "mandel_navigator.h"

class GradientWindow {
private:
    int _w, _h;
    bool _mouse_down;
    Widget* _wclick;

    GlutCanvas* _canvas;
    GradientBar* _gb;
    ColorPickerDisk* _cpd;
    ColorPickerBarV* _cpb;
    ScrollHColor* _sc;
    ScrollHMxit* _sc_mxit;
    Label* _l_mxit, * _l_super_sampling, * _l_ext_de;
    CheckerSuperSampling* _ch_super_sampling;
    CheckerExteriorDE* _ch_ext_de;
    MandelNavigator* _nav = nullptr;
    
    std::vector<Widget* > _widgets;
public:
    std::vector<int> _windows;

public:
    GradientWindow() = delete;

    explicit GradientWindow(int width, int height, const char* name);

    void ResetWindowSize();

    void Render();

    void MouseDownLC(int x, int y);

    void MouseUpLC(int x, int y);

    void MouseDownRC(int x, int y);

    void MouseUpRC(int x, int y);

    void MouseMove(int x, int y);

    void BindNavigator(MandelNavigator* nav);
};

extern GradientWindow* gradient_window;

void gradient_window_display();

void gradient_mouse_move(int x, int y);

void gradient_mouse_click(int button, int state, int x, int y);

#endif