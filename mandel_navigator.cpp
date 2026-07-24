#include <gmp.h>
#include <cassert>
#include <future>
#include <vector>

#include "navigator.h"
#include "mandel_perturbation.h"
#include "float_math.h"
#include "Image.h"

#include "mandel_navigator.h"

MandelNavigator::MandelNavigator(int width, int height, int sub, int max_iteration, double zoom_step, double zoom_time) 
        : Navigator(width, height, zoom_step, zoom_time), _sub(sub) {
    assert(sub % 2);
    _iter = new float[width * height * sub * sub];
    _mandel = new Mandel(width, height, max_iteration, sub, _iter);
    _shift_idx = (_w * _sub) * (_sub / 2);
    _require_update = true;
    mpf_init_set_d(_z_re, -0.5);
    mpf_init_set_d(_z_im, 0.0);
    mpf_init_set_d(_scale, 1.0);
    mpf_init(_t);
}

MandelNavigator::~MandelNavigator() {
    InterruptCompute();
    delete[] _iter;
    delete _mandel;
    mpf_clear(_z_re);
    mpf_clear(_z_im);
    mpf_clear(_scale);
    mpf_clear(_t);
}

void MandelNavigator::Reset() {
    mpf_set_d(_z_re, -0.5);
    mpf_set_d(_z_im, 0.0);
    mpf_set_d(_scale, 1.0);
    StartCompute();
}

void MandelNavigator::StartCompute() {
    InterruptCompute();
    for (int i = 0; i < _w * _h * _sub * _sub; ++i) _iter[i] = EMPTYPIXEL;
    auto compute_task = [this]() {
        this->_mandel->Compute(this->_z_re, this->_z_im, this->_scale, this->_mxit, this->_c_method);
        this->_require_update = true;
    };
    // _task = std::async(&Mandel::Compute, _mandel, _z_re, _z_im, _scale, _mxit, _c_method);
    _task = std::async(std::launch::async, compute_task);
}

void MandelNavigator::InterruptCompute() {
    if (_task.valid()) {
        _mandel->SetHalt(true);
        _task.wait();
        _mandel->SetHalt(false);
    }
}

bool MandelNavigator::IsComputing() {
    if (!_task.valid()) return false;
    std::chrono::milliseconds span(0);
    if (_task.wait_for(span) == std::future_status::ready) return false;
    return true;
}

std::string MandelNavigator::GetLocationText() const {
    int digits = std::max(18, (int)(mpf_get_prec(_z_re) * log(2.0) / log(10.0)));
    std::vector<char> buf((size_t)digits * 2 + 256);
    gmp_snprintf(buf.data(), buf.size(), "x: %.*Ff\r\ny: %.*Ff\r\nzoom: %.6Fe",
                 digits, _z_re, digits, _z_im, _scale);
    return std::string(buf.data());
}

void MandelNavigator::UpdateCoords() {
    mpf_set_d(_t, _k);
    mpf_mul(_scale, _scale, _t);
    int precision = std::abs(get_exp(_scale)) + 30;
    mpf_set_prec(_scale, precision);
    mpf_set_prec(_z_re, precision);
    mpf_set_prec(_z_im, precision);
    mpf_set_prec(_t, precision);
    _mandel->setPrecision(precision);
    mpf_set_d(_t, 2.0 * (_k - 1.0 + 2.0 * _display_dx / _w));
    mpf_div(_t, _t, _scale);
    mpf_sub(_z_re, _z_re, _t);

    mpf_set_d(_t, 2.0 * _h / _w * (_k - 1.0 + 2.0 * _display_dy / _h));
    mpf_div(_t, _t, _scale);
    mpf_sub(_z_im, _z_im, _t);
    gmp_printf("\nx: %.*Ff\ny: %.*Ff\nzoom: %.2Fe\n", (int)(precision * log(2) / log(10)), _z_re, (int)(precision * log(2) / log(10)), _z_im, _scale);
}

void MandelNavigator::UpdateBitmap(uint8_t* bitmap) {
    if (!_require_update && !IsComputing()) return;
    _require_update = false;
    const int stride = _w * _sub;
    const int c = _sub / 2;
    prepareColorFilter();   // rebuild palette box-filter integral once per frame
    const bool ss = (_c_method & ColoringMethod::SUPER_SAMPLING) != 0;
    for (int i = 0; i < _h; ++i) {
        for (int j = 0; j < _w; ++j) {
            int idx_bmp = (i * _w + j) * 3;
            int idx = (i * _sub + c) * stride + (j * _sub + c);
            if (_iter[idx] == EMPTYPIXEL) continue;
            // In SS mode a flagged pixel has its full sub-block computed (top-left
            // corner subpixel filled): average it. Otherwise this is the 1-sample
            // base layer -> analytic AA from the centre + its 4 pixel neighbours.
            if (ss && _iter[idx - (stride + 1) * c] != EMPTYPIXEL) {
                SmoothColor(bitmap + idx_bmp, idx, _c_method);
            } else {
                float vL = j > 0        ? _iter[idx - _sub]         : EMPTYPIXEL;
                float vR = j < _w - 1   ? _iter[idx + _sub]         : EMPTYPIXEL;
                float vU = i > 0        ? _iter[idx - stride * _sub] : EMPTYPIXEL;
                float vD = i < _h - 1   ? _iter[idx + stride * _sub] : EMPTYPIXEL;
                getColorAA(_iter[idx], vL, vR, vU, vD,
                           bitmap[idx_bmp], bitmap[idx_bmp + 1], bitmap[idx_bmp + 2], _c_method);
            }
        }
    }
}

void MandelNavigator::SmoothColor(uint8_t* bitmap_pixel, int idx, int _c_method) {
    double rs = 0, gs = 0, bs = 0;
    double ws = 0;
    // Average sub-samples in LINEAR light (sRGB-correct), matching getColorAA.
    for (int i = idx - _shift_idx; i <= idx + _shift_idx; i += _w * _sub) {
        for (int j = -_sub / 2; j <= _sub / 2; ++j) {
            float r, g, b;
            getColor(_iter[i + j], r, g, b, _c_method);
            int ri = (int)(r + 0.5f); if (ri < 0) ri = 0; else if (ri > 255) ri = 255;
            int gi = (int)(g + 0.5f); if (gi < 0) gi = 0; else if (gi > 255) gi = 255;
            int bi = (int)(b + 0.5f); if (bi < 0) bi = 0; else if (bi > 255) bi = 255;
            rs += g_srgb2lin[ri];
            gs += g_srgb2lin[gi];
            bs += g_srgb2lin[bi];
            ws += 1;
        }
    }
    rs /= ws;
    gs /= ws;
    bs /= ws;
    bitmap_pixel[0] = (uint8_t)(srgbEncode(rs) * 255.0 + 0.5);
    bitmap_pixel[1] = (uint8_t)(srgbEncode(gs) * 255.0 + 0.5);
    bitmap_pixel[2] = (uint8_t)(srgbEncode(bs) * 255.0 + 0.5);
}

void MandelNavigator::SetMxit(int mxit) {
    _mxit = mxit;
}

int MandelNavigator::GetCMethod() {
    return _c_method;
}

void MandelNavigator::SetCMethod(int c_method) {
    _c_method = c_method;
}

void MandelNavigator::SetRedisplay() {
    _require_update = true;
}