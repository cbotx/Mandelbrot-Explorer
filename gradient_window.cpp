#include <vector>

#include "gradient_window.h"



GradientWindow* gradient_window;

void gradient_window_display() {
    gradient_window->Render();
}

void gradient_mouse_move(int x, int y) {
    gradient_window->MouseMove(x, y);
}

void gradient_reshape(int x, int y) {
    gradient_window->ResetWindowSize();
}

void gradient_mouse_click(int button, int state, int x, int y) {
    if (state == GLUT_DOWN) {
        if (button == GLUT_LEFT_BUTTON) {
            gradient_window->MouseDownLC(x, y);
        } else if (button == GLUT_RIGHT_BUTTON) {
            gradient_window->MouseDownRC(x, y);
        }
    } else if (state == GLUT_UP) {
        if (button == GLUT_LEFT_BUTTON) {
            gradient_window->MouseUpLC(x, y);
        }  else if (button == GLUT_RIGHT_BUTTON) {
            gradient_window->MouseUpRC(x, y);
        }
    }
}

GradientWindow::GradientWindow(int width, int height, const char* name) : _w(width), _h(height) {
    _canvas = new GlutCanvas();
    _canvas->GlutCreateWindow(width, height, 710, 100, name);
    _canvas->GlutBindDisplayFunc(gradient_window_display);
    _canvas->GlutBindReshapeFunc(gradient_reshape);
    _canvas->GlutBindMouseFunc(gradient_mouse_click);
    _canvas->GlutBindMotionFunc(gradient_mouse_move);
    

    _cpd = new ColorPickerDisk(_canvas, 200, 200, 10, height - 10 - 200);
    _cpb = new ColorPickerBarV(_canvas, 20, 200, 220, height - 10 - 200);
    _gb = new GradientBar(_canvas, 750, 40, 10, height - 220 - 40);
    _sc = new ScrollHColor(_canvas, 750, 20, 10, height - 270 - 20, 10, 200, 60);
    _sc_mxit = new ScrollHMxit(_canvas, 370, 20, 370, height - 10 - 20, 10000, 1000000, 500000);
    _l_mxit = new Label(_canvas, 40, 20, 250, height - 8 - 20, "Max Iteration");
    _l_super_sampling = new Label(_canvas, 40, 20, 280, height - 38 - 20, "Supersampling");
    _l_ext_de = new Label(_canvas, 40, 20, 280, height - 70 - 18, "Exterior Distance Estimation");
    _ch_super_sampling = new CheckerSuperSampling(_canvas, 20, 20, 250, height - 40 - 20, false);
    _ch_ext_de = new CheckerExteriorDE(_canvas, 20, 20, 250, height - 70 - 20, true);


    _cpd->Bind(_cpb);
    _cpb->Bind(_cpd);

    _widgets.push_back(_cpd);
    _widgets.push_back(_cpb);
    _widgets.push_back(_gb);
    _widgets.push_back(_sc);
    _widgets.push_back(_sc_mxit);
    _widgets.push_back(_l_mxit);
    _widgets.push_back(_l_super_sampling);
    _widgets.push_back(_l_ext_de);
    _widgets.push_back(_ch_super_sampling);
    _widgets.push_back(_ch_ext_de);

    _gb->Attach(_cpd);
    _cpd->Attach(_gb);

    _ch_ext_de->BindGradientBar(_gb);

    _wclick = nullptr;
}

void GradientWindow::ResetWindowSize() {
    _canvas->GlutResetWindowSize();
    
    for (auto widget : _widgets) {
        widget->SetRedraw(true);
    }
    _canvas->Redisplay();
}

void GradientWindow::Render() {
    // _canvas->BeginDraw();
    bool flag = false;
    for (auto widget : _widgets) {
        if (widget->RequireRedraw()) {
            widget->Render();
            widget->SetRedraw(false);
            flag = true;
        }
    }
    if (flag) _canvas->EndDraw();
}

void GradientWindow::MouseDownLC(int x, int y) {
    y = _h - y;
    _mouse_down = true;
    for (auto& widget : _widgets) {
        if (widget->HitTest(x, y)) {
            _wclick = widget;
            widget->MouseDownLC(x, y);
            break;
        }
    }
}

void GradientWindow::MouseUpLC(int x, int y) {
    y = _h - y;
    _mouse_down = false;
    if (_wclick) _wclick->MouseUpLC(x, y);
    _wclick = nullptr;
    for (int window : _windows) {
        _canvas->Redisplay(window);
    }
}

void GradientWindow::MouseDownRC(int x, int y) {
    y = _h - y;
    for (auto& widget : _widgets) {
        if (widget->HitTest(x, y)) {
            widget->MouseDownRC(x, y);
            break;
        }
    }
}

void GradientWindow::MouseUpRC(int x, int y) {
    y = _h - y;

}

void GradientWindow::MouseMove(int x, int y) {
    y = _h - y;
    if (_wclick) _wclick->MouseMove(x, y);
}

void GradientWindow::BindNavigator(MandelNavigator* nav) {
    _nav = nav;
    _gb->BindNavigator(_nav);
    _sc->BindNavigator(_nav);
    _sc_mxit->BindNavigator(_nav);
    _ch_super_sampling->BindNavigator(_nav);
    _ch_ext_de->BindNavigator(_nav);
}