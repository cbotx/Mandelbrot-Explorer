#include <gmp.h>
#include <cassert>
#include <future>
#include <vector>
#include <string>
#include <algorithm>

#include "navigator.h"
#include "mandel_perturbation.h"
#include "float_math.h"
#include "Image.h"

#include "mandel_navigator.h"

MandelNavigator::MandelNavigator(int width, int height, int sub, int max_iteration, double zoom_step, double zoom_time) 
        : Navigator(width, height, zoom_step, zoom_time),
          _mxit(max_iteration), _sub(sub), _adaptive_sub(sub) {
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

void MandelNavigator::ConfigureSampling() {
    bool want_uniform = (_c_method & ColoringMethod::SUPER_SAMPLING)
                     && (_c_method & ColoringMethod::STRIPE_AVERAGE);
    if (want_uniform == _uniform_feather) return;

    delete[] _iter;
    delete _mandel;
    _uniform_feather = want_uniform;
    _sub = want_uniform ? 2 : _adaptive_sub;
    size_t count = (size_t)_w * _h * _sub * _sub;
    _iter = new float[count];
    if (want_uniform)
        _mandel = new Mandel(_w * _sub, _h * _sub, _mxit, 1, _iter);
    else
        _mandel = new Mandel(_w, _h, _mxit, _sub, _iter);
    _mandel->setPrecision((int)mpf_get_prec(_scale));
    _shift_idx = (_w * _sub) * (_sub / 2);
    _cache_valid = false;                       // sampling mode changed -> cache stale
}

void MandelNavigator::StartCompute() {
    InterruptCompute();
    ConfigureSampling();
    _cache_valid = false;                       // fractal changing -> phase cache stale
    for (int i = 0; i < _w * _h * _sub * _sub; ++i) _iter[i] = EMPTYPIXEL;
    auto compute_task = [this]() {
        this->_mandel->setDensity(color_density);   // for the SAC adaptive-SS detector
        int method = this->_uniform_feather
            ? (this->_c_method & ~ColoringMethod::SUPER_SAMPLING) : this->_c_method;
        this->_mandel->Compute(this->_z_re, this->_z_im, this->_scale, this->_mxit, method);
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

bool MandelNavigator::SetLocation(const std::string& xs, const std::string& ys, const std::string& ss) {
    int prec = (int)(std::max(xs.size(), ys.size()) * 3.3219) + 40;
    if (prec < 64) prec = 64;
    mpf_set_prec(_scale, prec); mpf_set_prec(_z_re, prec);
    mpf_set_prec(_z_im, prec); mpf_set_prec(_t, prec);
    _mandel->setPrecision(prec);
    if (mpf_set_str(_scale, ss.c_str(), 10) != 0 || mpf_sgn(_scale) <= 0) return false;
    if (mpf_set_str(_z_re, xs.c_str(), 10) != 0) return false;
    if (mpf_set_str(_z_im, ys.c_str(), 10) != 0) return false;
    JumpReset();
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
    bool computing = IsComputing();
    if (!_require_update && !computing && !_need_settle) return;
    _require_update = false;
    // While a render is in progress the frame is a coarse/partial preview: point-
    // sample it so the coarse blocks stay clean. Analytic AA (getColorAA) reads
    // neighbours, and across the artificial coarse-block boundaries it averages a
    // wide palette range into a muddy/grey edge -- the "grey block borders". AA is
    // correct only on the settled (fully computed) frame, so defer it until then.
    // _need_settle guarantees exactly one AA pass after computing finishes even if
    // the settle tick would otherwise be skipped.
    const bool settled = !computing;
    _need_settle = !settled;
    prepareColorFilter();   // also initializes the linear-light LUT
    if (_uniform_feather) {
        const int stride = _w * _sub;
        if (settled) _baseUsub.assign((size_t)_w * _h * _sub * _sub, EMPTYPIXEL);
#pragma omp parallel for schedule(dynamic, 8)
        for (int i = 0; i < _h; ++i) {
            for (int j = 0; j < _w; ++j) {
                double rs = 0, gs = 0, bs = 0; int n = 0;
                for (int a = 0; a < _sub; ++a) for (int b = 0; b < _sub; ++b) {
                    int sidx = (i * _sub + a) * stride + j * _sub + b;
                    float v = _iter[sidx];
                    if (settled) _baseUsub[sidx] = (v == EMPTYPIXEL) ? EMPTYPIXEL
                                                                     : colorBaseIndex(v, _c_method);
                    if (v == EMPTYPIXEL) continue;
                    float r, g, bl; getColor(v, r, g, bl, _c_method);
                    int ri = std::clamp((int)(r + 0.5f), 0, 255);
                    int gi = std::clamp((int)(g + 0.5f), 0, 255);
                    int bi = std::clamp((int)(bl + 0.5f), 0, 255);
                    rs += g_srgb2lin[ri]; gs += g_srgb2lin[gi]; bs += g_srgb2lin[bi]; ++n;
                }
                if (!n) continue;
                uint8_t* p = bitmap + ((size_t)i * _w + j) * 3;
                p[0] = (uint8_t)(srgbEncode(rs / n) * 255.0 + 0.5);
                p[1] = (uint8_t)(srgbEncode(gs / n) * 255.0 + 0.5);
                p[2] = (uint8_t)(srgbEncode(bs / n) * 255.0 + 0.5);
            }
        }
        if (settled) { _cache_valid = true; _cache_density = color_density; _cache_method = _c_method; }
        return;
    }
    const int stride = _w * _sub;
    const int c = _sub / 2;
    const bool ss = (_c_method & ColoringMethod::SUPER_SAMPLING) != 0;
    // The phase cache is only usable when every pixel took the analytic-AA path
    // (SS off): SS-flagged pixels are averaged sub-blocks that this cache does
    // not store, so those modes fall back to a full re-colour during animation.
    const bool cacheable = !ss && settled;
    if (cacheable) { _baseU.assign((size_t)_w * _h, -1.0f); _widthC.assign((size_t)_w * _h, 0.0f); }
#pragma omp parallel for schedule(dynamic, 8)
    for (int i = 0; i < _h; ++i) {
        for (int j = 0; j < _w; ++j) {
            int idx_bmp = (i * _w + j) * 3;
            int idx = (i * _sub + c) * stride + (j * _sub + c);
            if (_iter[idx] == EMPTYPIXEL) continue;
            // In SS mode a flagged pixel has its full sub-block computed (top-left
            // corner subpixel filled): average it. Otherwise this is the 1-sample
            // base layer -> analytic AA from the centre + its 4 pixel neighbours,
            // but ONLY once settled (see above) -- while computing, point-sample.
            if (ss && _iter[idx - (stride + 1) * c] != EMPTYPIXEL) {
                SmoothColor(bitmap + idx_bmp, idx, _c_method);
            } else if (settled) {
                float vL = j > 0        ? _iter[idx - _sub]         : EMPTYPIXEL;
                float vR = j < _w - 1   ? _iter[idx + _sub]         : EMPTYPIXEL;
                float vU = i > 0        ? _iter[idx - stride * _sub] : EMPTYPIXEL;
                float vD = i < _h - 1   ? _iter[idx + stride * _sub] : EMPTYPIXEL;
                float baseU, width;
                colorAnalyzeAA(_iter[idx], vL, vR, vU, vD, _c_method, baseU, width);
                colorShadeAA(baseU, width, color_phase,
                             bitmap[idx_bmp], bitmap[idx_bmp + 1], bitmap[idx_bmp + 2]);
                if (cacheable) { _baseU[(size_t)i * _w + j] = baseU; _widthC[(size_t)i * _w + j] = width; }
            } else {
                float r, g, b; getColor(_iter[idx], r, g, b, _c_method);
                bitmap[idx_bmp]     = (uint8_t)std::clamp((int)(r + 0.5f), 0, 255);
                bitmap[idx_bmp + 1] = (uint8_t)std::clamp((int)(g + 0.5f), 0, 255);
                bitmap[idx_bmp + 2] = (uint8_t)std::clamp((int)(b + 0.5f), 0, 255);
            }
        }
    }
    if (cacheable) { _cache_valid = true; _cache_density = color_density; _cache_method = _c_method; }
}

void MandelNavigator::RecolorPhase(uint8_t* bitmap) {
    // Fast animation re-colour: re-shade the cached settled frame with the live
    // color_phase. If the cache is stale (density/method changed) or unusable
    // (SS-averaged mode, or a render in progress), fall back to a full re-colour.
    bool computing = IsComputing();
    if (!_cache_valid || computing ||
        _cache_density != color_density || _cache_method != _c_method) {
        _require_update = true;
        UpdateBitmap(bitmap);
        return;
    }
    // The palette (and thus its box-filter integral) may have been edited since
    // the cache was built; refresh it so an animated phase uses current colours.
    // Cheap (colP*3) relative to the per-pixel re-colour. The per-pixel analysis
    // (baseU/width) is palette-independent, so it stays cached.
    prepareColorFilter();
    if (_uniform_feather) {
        const int stride = _w * _sub;
#pragma omp parallel for schedule(dynamic, 8)
        for (int i = 0; i < _h; ++i) {
            for (int j = 0; j < _w; ++j) {
                double rs = 0, gs = 0, bs = 0; int n = 0;
                for (int a = 0; a < _sub; ++a) for (int b = 0; b < _sub; ++b) {
                    float bu = _baseUsub[(size_t)(i * _sub + a) * stride + j * _sub + b];
                    if (bu == EMPTYPIXEL) continue;
                    int ri, gi, bi;
                    if (bu < 0) { ri = gi = bi = 0; }
                    else {
                        int x = ((int)(bu + color_phase)) % colP; if (x < 0) x += colP;
                        ri = std::clamp((int)(color_map[0][x] + 0.5f), 0, 255);
                        gi = std::clamp((int)(color_map[1][x] + 0.5f), 0, 255);
                        bi = std::clamp((int)(color_map[2][x] + 0.5f), 0, 255);
                    }
                    rs += g_srgb2lin[ri]; gs += g_srgb2lin[gi]; bs += g_srgb2lin[bi]; ++n;
                }
                if (!n) continue;
                uint8_t* p = bitmap + ((size_t)i * _w + j) * 3;
                p[0] = (uint8_t)(srgbEncode(rs / n) * 255.0 + 0.5);
                p[1] = (uint8_t)(srgbEncode(gs / n) * 255.0 + 0.5);
                p[2] = (uint8_t)(srgbEncode(bs / n) * 255.0 + 0.5);
            }
        }
        return;
    }
#pragma omp parallel for schedule(dynamic, 8)
    for (int i = 0; i < _h; ++i) {
        for (int j = 0; j < _w; ++j) {
            size_t p = (size_t)i * _w + j;
            uint8_t* q = bitmap + p * 3;
            colorShadeAA(_baseU[p], _widthC[p], color_phase, q[0], q[1], q[2]);
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

void MandelNavigator::GetView(mpf_t re, mpf_t im, mpf_t scale) const {
    mpf_set_prec(re, mpf_get_prec(_z_re)); mpf_set(re, _z_re);
    mpf_set_prec(im, mpf_get_prec(_z_im)); mpf_set(im, _z_im);
    mpf_set_prec(scale, mpf_get_prec(_scale)); mpf_set(scale, _scale);
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