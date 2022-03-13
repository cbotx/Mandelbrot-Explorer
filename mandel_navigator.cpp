#include <mpir.h>
#include <future>

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
    _task = std::async(compute_task);
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
    std::chrono::milliseconds span(50);
    if (_task.wait_for(span) == std::future_status::ready) return false;
    return true;
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
    for (int i = 0; i < _h; ++i) {
        for (int j = 0; j < _w; ++j) {
            int idx_bmp = (i * _w + j) * 3;
            int idx = (i * _sub + _sub / 2) * _w * _sub + (j * _sub + _sub / 2);
            if (_c_method & ColoringMethod::SUPER_SAMPLING) {
                if (_iter[idx] != EMPTYPIXEL) {
                    if (_iter[idx - (_w * _sub + 1) * (_sub / 2)] == EMPTYPIXEL) {
                        getColor(_iter[idx], bitmap[idx_bmp], bitmap[idx_bmp + 1], bitmap[idx_bmp + 2], _c_method);
                    } else {
                        SmoothColor(bitmap + idx_bmp, idx, _c_method);
                    }
                }
            } else {
                if (_iter[idx] != EMPTYPIXEL) {
                    getColor(_iter[idx], bitmap[idx_bmp], bitmap[idx_bmp + 1], bitmap[idx_bmp + 2], _c_method);
                }
            }

        }
    }
}

void MandelNavigator::SmoothColor(uint8_t* bitmap_pixel, int idx, int _c_method) {
    float rs, gs, bs;
    float ws = 0;
    
    for (int i = idx - _shift_idx; i <= idx + _shift_idx; i += _w * _sub) {
        for (int j = -_sub / 2; j <= _sub / 2; ++j) {
            float r, g, b;
            getColor(_iter[i + j], r, g, b, _c_method);
            rs += r * r;
            gs += g * g;
            bs += b * b;
            ws += 1;
        }
    }
    rs /= ws;
    gs /= ws;
    bs /= ws;
    // bitmap_pixel[0] = 255;
    // bitmap_pixel[1] = 0;
    // bitmap_pixel[2] = 0;
    bitmap_pixel[0] = sqrt(rs);
    bitmap_pixel[1] = sqrt(gs);
    bitmap_pixel[2] = sqrt(bs);
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