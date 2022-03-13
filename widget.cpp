#include <iostream>
#include <algorithm>
#include <math.h>

#include "widget.h"
#include "interpolate.h"
#include "color.h"

Widget::Widget(GlutCanvas* canvas, int width, int height, int x, int y) : _canvas(canvas), _w(width), _h(height), _x(x), _y(y) {
    
}

bool Widget::RequireRedraw() { return _require_redraw; }

void Widget::SetRedraw(bool redraw) { _require_redraw = redraw; }

void Widget::Update() {}

void Widget::MouseDownLC(int x, int y) {}

void Widget::MouseUpLC(int x, int y) {}

void Widget::MouseDownRC(int x, int y) {}

void Widget::MouseUpRC(int x, int y) {}

void Widget::MouseMove(int x, int y) {}

bool Widget::HitTest(int x, int y) {
    return (x >= _x && x <= _x + _w && y >= _y && y <= _y + _h);
}


Label::Label(GlutCanvas* canvas, int width, int height, int x, int y, const char* text)
            : Widget(canvas, width, height, x, y), _text(text) {

}

void Label::Render() {
    _canvas->RenderText(_x, _y, 1, 1, 1, reinterpret_cast<const unsigned char*>(_text));
}


GradientBar::GradientBar(GlutCanvas* canvas, int width, int height, int x, int y) : Widget(canvas, width, height, x, y) {
    _texture = _canvas->CreateTexture(false);
    _color = new uint8_t[width * 3];
    _mp = new float[_w * 3];
    _mp_full = new float[colP * 3];
    SetSceneSnowy();
}

void GradientBar::SetSceneSunRise() {
    static const int N = 6;
    static const float pos[N] = { 0.0, 0.16, 0.42, 0.6425, 0.8575, 1.0 };
    static const float col[3][N] = { {0, 32, 237, 255, 0, 0}, {70, 107, 255, 170, 2, 70}, {100, 203, 255, 0, 0, 100} };
    _p.clear();
    for (int i = 0; i < N - 1; ++i) {
        float l_s, a_s, b_s;
        rgb2lab(col[0][i], col[1][i], col[2][i], l_s, a_s, b_s);
        _p.push_back({pos[i], l_s, a_s, b_s});
    }
    Update();

    int idx = (int)(_p[0][0] * _w) * 3;
    Notify(_color[idx], _color[idx + 1], _color[idx + 2]);
}


void GradientBar::SetSceneSnowy() {
    static const int N = 4;
    static const float pos[N] = { 0.0, 0.12, 0.34, 1.0 };
    static const float col[3][N] = { {255, 32, 214, 255}, {255, 107, 223, 255}, {255, 203, 255, 255} };
    _p.clear();
    for (int i = 0; i < N - 1; ++i) {
        float l_s, a_s, b_s;
        rgb2lab(col[0][i], col[1][i], col[2][i], l_s, a_s, b_s);
        _p.push_back({pos[i], l_s, a_s, b_s});
    }
    Update();

    int idx = (int)(_p[0][0] * _w) * 3;
    Notify(_color[idx], _color[idx + 1], _color[idx + 2]);
}

void GradientBar::BindNavigator(MandelNavigator* nav) {
    _nav = nav;
}

void GradientBar::Update() {
    auto p = _p;
    sort(p.begin(), p.end());
    float* xs = new float[p.size() + 2];
    float* ys = new float[p.size() + 2];
    for (int c = 0; c < 3; ++c) {
        for (int i = 0; i < p.size(); ++i) {
            xs[i + 1] = p[i][0];
            ys[i + 1] = p[i][1 + c];
        }
        xs[0] = xs[p.size()] - 1;
        ys[0] = ys[p.size()];
        xs[p.size() + 1] = xs[1] + 1;
        ys[p.size() + 1] = ys[1];
        mono_cubic_interpolate(xs, ys, p.size() + 2, _mp + c * _w, _w);
        // cos_interpolate(xs, ys, p.size(), _mp + c * _w, _w);
        mono_cubic_interpolate(xs, ys, p.size() + 2, _mp_full + c * colP, colP);
        // cos_interpolate(xs, ys, p.size(), _mp_full + c * colP, colP);
    }
    for (int i = 0; i < colP; ++i) {
        float r, g, b;
        lab2rgb(_mp_full[i], _mp_full[i + colP], _mp_full[i + colP * 2], r, g, b);
        color_map[0][i] = r;
        color_map[1][i] = g;
        color_map[2][i] = b;
    }
    for (int i = 0; i < _w; ++i) {
        float r, g, b;
        lab2rgb(_mp[i], _mp[i + _w], _mp[i + _w * 2], r, g, b);
        _color[i * 3 + 0] = r;
        _color[i * 3 + 1] = g;
        _color[i * 3 + 2] = b;
    }
    delete[] xs;
    delete[] ys;
    if (_nav) _nav->SetRedisplay();
    SetRedraw(true);
    _canvas->Redisplay();
}

void GradientBar::Render() {
    std::vector<std::array<float, 2>> v;
    v.push_back({0.f + _x - DMSIZE, 0.f + _y});
    v.push_back({0.f + _x, 0.f + _y});
    v.push_back({0.f + _x, 0.f + _y + _h});
    v.push_back({0.f + _x - DMSIZE, 0.f + _y + _h});
    _canvas->DrawPolygon(v, 0, 0, 0, true);
    
    v.clear();
    v.push_back({0.f + _x + _w, 0.f + _y});
    v.push_back({0.f + _x + _w + DMSIZE, 0.f + _y});
    v.push_back({0.f + _x + _w + DMSIZE, 0.f + _y + _h});
    v.push_back({0.f + _x + _w, 0.f + _y + _h});
    _canvas->DrawPolygon(v, 0, 0, 0, true);

    _canvas->BitmapToTexture(_texture, _color, _w, 1);
    _canvas->RenderTextureP(_texture, {0, 0, 0, 1, 1, 1, 1, 0}, {_x, _y, _x, _y + _h, _x + _w, _y + _h, _x + _w, _y});
    
    for (int i = 0; i < _p.size(); ++i) {
        auto& p = _p[i];
        int idx = (int)(p[0] * _w) * 3;
        v.clear();
        v.push_back({p[0] * _w + _x - DMSIZE, _y + _h / 2.f});
        v.push_back({p[0] * _w + _x, -DMSIZE + _y + _h / 2.f});
        v.push_back({p[0] * _w + _x + DMSIZE, _y + _h / 2.f});
        v.push_back({p[0] * _w + _x, DMSIZE + _y + _h / 2.f});
        _canvas->DrawPolygon(v, 1 - _color[idx] / 255.f, 1 - _color[idx + 1] / 255.f, 1 - _color[idx + 2] / 255.f, i == _selected);
    }
}


void GradientBar::UpdateColor(float R, float G, float B) {
    float l_s, a_s, b_s;
    rgb2lab(R, G, B, l_s, a_s, b_s);
    _p[_selected][1] = l_s;
    _p[_selected][2] = a_s;
    _p[_selected][3] = b_s;
    Update();
}

void GradientBar::MouseDownLC(int x, int y) {
    for (int i = 0; i < _p.size(); ++i) {
        if (std::abs(x - _p[i][0] * _w - _x) <= DMSIZE) {
            _selected = _drag_idx = i;
            _p[_drag_idx][0] = std::min(std::max((x - _x) * 1.0 / _w, 0.0), 1.0);
            Update();
            int idx = (int)(_p[_selected][0] * _w) * 3;
            Notify(_color[idx], _color[idx + 1], _color[idx + 2]);
            return;
        }
    }
    _p.push_back({1.f * (x - _x) / _w, _mp[x - _x], _mp[_w + x - _x], _mp[_w * 2 + x - _x]});
    _selected = _drag_idx = _p.size() - 1;
    Update();
    int idx = (int)(_p[_selected][0] * _w) * 3;
    Notify(_color[idx], _color[idx + 1], _color[idx + 2]);
}

void GradientBar::MouseUpLC(int x, int y) {
    _drag_idx = -1;
}

void GradientBar::MouseDownRC(int x, int y) {
    for (int i = 0; i < _p.size(); ++i) {
        if (std::abs(x - _p[i][0] * _w - _x) <= DMSIZE) {
            _selected = 0;
            _p.erase(_p.begin() + i);
            if (_p.empty()) _p.push_back({0, 0, 0, 0});
            Update();
            
            int idx = (int)(_p[_selected][0] * _w) * 3;
            Notify(_color[idx], _color[idx + 1], _color[idx + 2]);
            break;
        }
    }
}

void GradientBar::MouseUpRC(int x, int y) {}

void GradientBar::MouseMove(int x, int y) {
    if (_drag_idx >= 0) {
        _p[_drag_idx][0] = std::min(std::max((x - _x) * 1.0 / _w, 0.0), 1.0);
        Update();
    }
}


ColorPickerDisk::ColorPickerDisk(GlutCanvas* canvas, int width, int height, int x, int y) : Widget(canvas, width, height, x, y) {
    _texture = _canvas->CreateTexture(false);
    _color = new uint8_t[width * height * 3];
    _fh = _fs = 0;
    _fv = 1;
    Update();
}

void ColorPickerDisk::Bind(ColorPickerBarV* cpb) {
    _cpb = cpb;
}

void ColorPickerDisk::Update() {
    for (int i = 0; i < _h; ++i) {
        for (int j = 0; j < _w; ++j) {
            int idx = (i * _w + j) * 3;
            float r, g, b;
            float rad = std::sqrt(std::pow(i - _h / 2.0, 2) + std::pow(j - _w / 2.0, 2));
            if (rad <= RADIUS) {
                hsv2rgb(std::atan2(i - _h / 2.0, j - _w / 2.0) / 3.14159 * 180, rad / RADIUS, _fv, r, g, b);
            } else {
                r = g = b = 0;
            }
            _color[idx] = r;
            _color[idx + 1] = g;
            _color[idx + 2] = b;
        }
    }
    _canvas->BitmapToTexture(_texture, _color, _w, _h);
    _p = { (int)(std::cos(_fh / 180 * 3.14159) * _fs * RADIUS + _w / 2.0), (int)(std::sin(_fh / 180 * 3.14159) * _fs * RADIUS + _h / 2.0) };

    SetRedraw(true);
    _canvas->Redisplay();
}

void ColorPickerDisk::Render() {
    
    std::vector<std::array<float, 2>> v;
    v.push_back({_x - DMSIZE, _y - DMSIZE});
    v.push_back({_x + _w + DMSIZE, _y - DMSIZE});
    v.push_back({_x + _w + DMSIZE, _y + _h + DMSIZE});
    v.push_back({_x - DMSIZE, _y + _h + DMSIZE});
    _canvas->DrawPolygon(v, 0, 0, 0, true);

    _canvas->RenderTextureP(_texture, {0, 0, 0, 1, 1, 1, 1, 0}, {_x, _y, _x, _y + _h, _x + _w, _y + _h, _x + _w, _y});
    int idx = (_p[1] * _w + _p[0]) * 3;
    v.clear();
    v.push_back({_p[0] + _x - DMSIZE, 0.f + _p[1] + _y});
    v.push_back({0.f + _p[0] + _x, -DMSIZE + _p[1] + _y});
    v.push_back({_p[0] + _x + DMSIZE, 0.f + _p[1] + _y});
    v.push_back({0.f + _p[0] + _x, DMSIZE + _p[1] + _y});
    _canvas->DrawPolygon(v, 1 - _color[idx] / 255.f, 1 - _color[idx + 1] / 255.f, 1 - _color[idx + 2] / 255.f);
}

void ColorPickerDisk::UpdateColor(float R, float G, float B) {
    rgb2hsv(R, G, B, _fh, _fs, _fv);
    Update();
    _cpb->Update();
}

void ColorPickerDisk::MouseDownLC(int x, int y) {
    MouseMove(x, y);
}

void ColorPickerDisk::MouseUpLC(int x, int y) {}

void ColorPickerDisk::MouseDownRC(int x, int y) {}

void ColorPickerDisk::MouseUpRC(int x, int y) {}

void ColorPickerDisk::MouseMove(int x, int y) {
    float rad = std::sqrt(std::pow(x - _x - _w / 2.0, 2) + std::pow(y - _y - _h / 2.0, 2));
    if (rad < (RADIUS - 1)) {
        _p = {x - _x, y - _y};
    } else {
        _p = {(int)((x - _x - _w / 2.0) / rad * (RADIUS - 1) + _w / 2.0), (int)((y - _y - _h / 2.0) / rad * (RADIUS - 1) + _h / 2.0)};
    }
    int idx = (_p[1] * _w + _p[0]) * 3;
    float _;
    rgb2hsv(_color[idx], _color[idx + 1], _color[idx + 2], _fh, _fs, _);
    Notify(_color[idx], _color[idx + 1], _color[idx + 2]);
    
    SetRedraw(true);
    _canvas->Redisplay();
}



ColorPickerBarV::ColorPickerBarV(GlutCanvas* canvas, int width, int height, int x, int y) : Widget(canvas, width, height, x, y) {
    _texture = _canvas->CreateTexture(false);
    
    _color = new uint8_t[MIN_W * height * 3];
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < MIN_W; ++j) {
            int idx = (i * MIN_W + j) * 3;
            _color[idx + 0] = i * 255.0 / height;
            _color[idx + 1] = i * 255.0 / height;
            _color[idx + 2] = i * 255.0 / height;
        }
    }
    _p = _h;
    _canvas->BitmapToTexture(_texture, _color, MIN_W, _h);
}

void ColorPickerBarV::Bind(ColorPickerDisk* cpd) {
    _cpd = cpd;
}

void ColorPickerBarV::Update() {
    _p = _cpd->_fv * _h;
    SetRedraw(true);
    _canvas->Redisplay();
}

void ColorPickerBarV::Render() {
    std::vector<std::array<float, 2>> v;
    v.push_back({_x - DMSIZE, _y - DMSIZE});
    v.push_back({_x + _w + DMSIZE, _y - DMSIZE});
    v.push_back({_x + _w + DMSIZE, _y + _h + DMSIZE});
    v.push_back({_x - DMSIZE, _y + _h + DMSIZE});
    _canvas->DrawPolygon(v, 0, 0, 0, true);

    _canvas->RenderTextureP(_texture, {0, 0, 0, 1, 1, 1, 1, 0}, {_x, _y, _x, _y + _h, _x + _w, _y + _h, _x + _w, _y});
    int idx = (_p * MIN_W) * 3;
    v.clear();
    v.push_back({0.f + _x + 1, 0.f + _p + _y - DMSIZE});
    v.push_back({0.f + _x + _w, 0.f + _p + _y - DMSIZE});
    v.push_back({0.f + _x + _w, 0.f + _p + _y + DMSIZE});
    v.push_back({0.f + _x + 1, 0.f + _p + _y + DMSIZE});
    _canvas->DrawPolygon(v, 1 - _color[idx] / 255.f, 1 - _color[idx + 1] / 255.f, 1 - _color[idx + 2] / 255.f, true);

}

void ColorPickerBarV::MouseDownLC(int x, int y) {
    MouseMove(x, y);
}

void ColorPickerBarV::MouseUpLC(int x, int y) {}

void ColorPickerBarV::MouseDownRC(int x, int y) {}

void ColorPickerBarV::MouseUpRC(int x, int y) {}

void ColorPickerBarV::MouseMove(int x, int y) {
    _p = std::min(std::max(y - _y, 0), _h);
    _cpd->_fv = 1.0 * _p / _h;
    _cpd->Update();
    
    int idx = (_cpd->_p[1] * _cpd->_w + _cpd->_p[0]) * 3;
    _cpd->Notify(_cpd->_color[idx], _cpd->_color[idx + 1], _cpd->_color[idx + 2]);
    
    SetRedraw(true);
}


ScrollHBase::ScrollHBase(GlutCanvas* canvas, int width, int height, int x, int y, int mn, int mx, int val)
        : Widget(canvas, width, height, x, y), _min(mn), _max(mx), _val(val) {
    _ratio = (_val - _min) / (_max - _min);
}

void ScrollHBase::Update() {
    SetRedraw(true);
    _canvas->Redisplay();
}

void ScrollHBase::Render() {
    std::vector<std::array<float, 2>> v;
    v.push_back({-0.5f + _x, -0.5f + _y});
    v.push_back({0.f + _x + _w, -0.5f + _y});
    v.push_back({0.f + _x + _w, 0.f + _y + _h});
    v.push_back({-0.5f + _x, 0.f + _y + _h});
    _canvas->DrawPolygon(v, 0, 0, 0, true);

    if ((_w - 2 * DMSIZE) * _ratio > DMSIZE) {
        v.clear();
        v.push_back({0.f + _x + DMSIZE, 0.f + _y + _h / 2.f});
        v.push_back({0.f + _x + (_w - 2 * DMSIZE) * _ratio, 0.f + _y + _h / 2.f});
        _canvas->DrawPolygon(v, 1, 1, 1, false);
    }

    if ((_w - 2 * DMSIZE) * _ratio + 2 * DMSIZE < _w - DMSIZE) {
        v.clear();
        v.push_back({0.f + _x + (_w - 2 * DMSIZE) * _ratio + 2 * DMSIZE, 0.f + _y + _h / 2.f});
        v.push_back({0.f + _x + _w - DMSIZE, 0.f + _y + _h / 2.f});
        _canvas->DrawPolygon(v, 1, 1, 1, false);
    }

    v.clear();
    v.push_back({0.f + _x + (_w - 2 * DMSIZE) * _ratio, 0.f + _y});
    v.push_back({0.f + _x + (_w - 2 * DMSIZE) * _ratio + 2 * DMSIZE, 0.f + _y});
    v.push_back({0.f + _x + (_w - 2 * DMSIZE) * _ratio + 2 * DMSIZE, 0.f + _y + _h});
    v.push_back({0.f + _x + (_w - 2 * DMSIZE) * _ratio, 0.f + _y + _h});
    _canvas->DrawPolygon(v, 1, 1, 1, false);

}

void ScrollHBase::MouseDownLC(int x, int y) {
    MouseMove(x, y);
}

void ScrollHBase::MouseUpLC(int x, int y) {}

void ScrollHBase::MouseDownRC(int x, int y) {}

void ScrollHBase::MouseUpRC(int x, int y) {}

void ScrollHBase::MouseMove(int x, int y) {
    _ratio = std::min(std::max((x - _x - DMSIZE) / (_w - 2.f * DMSIZE), 0.f), 1.f);
    _val = _ratio * (_max - _min) + _min;
    SetRedraw(true);
    valueOnChange();
}



ScrollHColor::ScrollHColor(GlutCanvas* canvas, int width, int height, int x, int y, int mn, int mx, int val)
        : ScrollHBase(canvas, width, height, x, y, mn, mx, val) {
}

void ScrollHColor::BindNavigator(MandelNavigator* nav) {
    _nav = nav;
}

void ScrollHColor::valueOnChange() {
    color_density = _val;
    if (_nav) _nav->SetRedisplay();
    _canvas->Redisplay();
}



ScrollHMxit::ScrollHMxit(GlutCanvas* canvas, int width, int height, int x, int y, int mn, int mx, int val)
        : ScrollHBase(canvas, width, height, x, y, mn, mx, val) {
}

void ScrollHMxit::BindNavigator(MandelNavigator* nav) {
    _nav = nav;
    _nav->SetMxit(_val);
}

void ScrollHMxit::valueOnChange() {
    if (_nav) _nav->SetMxit(_val);
}




CheckerBase::CheckerBase(GlutCanvas* canvas, int width, int height, int x, int y, bool state)
        : Widget(canvas, width, height, x, y), _state(state) {
}

void CheckerBase::Update() {   
    SetRedraw(true);
    _canvas->Redisplay();
}

void CheckerBase::Render() {
    std::vector<std::array<float, 2>> v;
    v.push_back({0.f + _x, 0.f + _y});
    v.push_back({0.f + _x + _w, 0.f + _y});
    v.push_back({0.f + _x + _w, 0.f + _y + _h});
    v.push_back({0.f + _x, 0.f + _y + _h});
    _canvas->DrawPolygon(v, 0, 0, 0, true);
    _canvas->DrawPolygon(v, 1, 1, 1, false);

    if (_state) {
        v.clear();
        float margin_size = (_w - _w * SIZE_RATIO) * 0.5f;
        v.push_back({0.f + _x + margin_size, 0.f + _y + margin_size});
        v.push_back({0.f + _x + _w - margin_size, 0.f + _y + margin_size});
        v.push_back({0.f + _x + _w - margin_size, 0.f + _y + _h - margin_size});
        v.push_back({0.f + _x + margin_size, 0.f + _y + _h - margin_size});
        _canvas->DrawPolygon(v, 1, 1, 1, true);
    }
}

void CheckerBase::MouseDownLC(int x, int y) {
    _state = !_state;
    SetRedraw(true);
    valueOnChange();
}


CheckerSuperSampling::CheckerSuperSampling(GlutCanvas* canvas, int width, int height, int x, int y, bool state)
        : CheckerBase(canvas, width, height, x, y, state) {
}

void CheckerSuperSampling::BindNavigator(MandelNavigator* nav) {
    _nav = nav;
}

void CheckerSuperSampling::valueOnChange() {
    if (_nav) {
        int c_method = _nav->GetCMethod();
        c_method |= ColoringMethod::SUPER_SAMPLING;
        if (!_state) c_method ^= ColoringMethod::SUPER_SAMPLING;
        _nav->SetCMethod(c_method);
        _nav->StartCompute();
        _canvas->Redisplay();
    }
}



CheckerExteriorDE::CheckerExteriorDE(GlutCanvas* canvas, int width, int height, int x, int y, bool state)
        : CheckerBase(canvas, width, height, x, y, state) {
}

void CheckerExteriorDE::BindNavigator(MandelNavigator* nav) {
    _nav = nav;
}

void CheckerExteriorDE::BindGradientBar(GradientBar* gb) {
    _gb = gb;
}

void CheckerExteriorDE::valueOnChange() {
    if (_gb) {
        if (_state) {
            _gb->SetSceneSnowy();   
        } else {
            _gb->SetSceneSunRise();
        }
    }

    if (_nav) {
        int c_method = _nav->GetCMethod();
        c_method |= ColoringMethod::EXTERIOR_DIST_EST;
        if (!_state) c_method ^= ColoringMethod::EXTERIOR_DIST_EST;
        _nav->SetCMethod(c_method);
        _nav->StartCompute();
        _canvas->Redisplay();
    }
}
