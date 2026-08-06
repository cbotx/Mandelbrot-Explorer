#include <gmp.h>
#include <cassert>
#include <cmath>
#include <future>
#include <limits>
#include <vector>
#include <string>
#include <algorithm>
#include <cstdlib>

#include "navigator.h"
#include "mandel_perturbation.h"
#include "float_math.h"
#include "Image.h"

#include "mandel_navigator.h"

// The engine marks set-interior pixels with the sentinel -2.0f (EMPTYPIXEL=-10
// marks not-yet-computed). Relief skips both so they carry no slope.
// (INTERIOR_SENTINEL is declared in Image.h.)
// Cache marker in _baseU: this pixel is an SS SmoothColor sub-block; its phase
// re-colour reads per-subpixel base indices from _baseUsub instead of colorShadeAA.
static constexpr float SS_SENTINEL = -3.0f;

MandelNavigator::MandelNavigator(int width, int height, int sub, int max_iteration, double zoom_step, double zoom_time) 
        : Navigator(width, height, zoom_step, zoom_time),
          _mxit(max_iteration), _sub(sub), _adaptive_sub(sub) {
    assert(sub % 2);
    _iter = new float[width * height * sub * sub];
    _normal = new float[width * height * sub * sub];
    _mandel = new Mandel(width, height, max_iteration, sub, _iter);
    _mandel->setNormalBuffer(_normal);
    _shift_idx = (_w * _sub) * (_sub / 2);
    _require_update = true;
    mpf_init_set_d(_z_re, -0.5);
    mpf_init_set_d(_z_im, 0.0);
    mpf_init_set_d(_scale, 1.0);
    mpf_init_set_d(_julia_c_re, 0.0);
    mpf_init_set_d(_julia_c_im, 0.0);
    mpf_init(_saved_mandel_re);
    mpf_init(_saved_mandel_im);
    mpf_init(_saved_mandel_scale);
    mpf_init(_t);
    _backend = createComputeBackend(getenv("MANDEL_BACKEND"));
}

MandelNavigator::~MandelNavigator() {
    InterruptCompute();
    delete[] _iter;
    delete[] _normal;
    delete _mandel;
    mpf_clear(_z_re);
    mpf_clear(_z_im);
    mpf_clear(_scale);
    mpf_clear(_julia_c_re);
    mpf_clear(_julia_c_im);
    mpf_clear(_saved_mandel_re);
    mpf_clear(_saved_mandel_im);
    mpf_clear(_saved_mandel_scale);
    mpf_clear(_t);
}

void MandelNavigator::Reset() {
    mpf_set_d(_z_re, IsJulia() ? 0.0 : -0.5);
    mpf_set_d(_z_im, 0.0);
    mpf_set_d(_scale, 1.0);
    StartCompute();
}

void MandelNavigator::Resize(int width, int height) {
    if (width < 8) width = 8;
    if (height < 8) height = 8;
    if (width == _w && height == _h) return;
    InterruptCompute();
    _w = width; _h = height;
    // Rebuild the sample buffers + engine at the new size (mirrors ConfigureSampling).
    delete[] _iter;
    delete[] _normal;
    delete _mandel;
    size_t count = (size_t)_w * _h * _sub * _sub;
    _iter = new float[count];
    _normal = new float[count];
    if (_uniform_feather)
        _mandel = new Mandel(_w * _sub, _h * _sub, _mxit, 1, _iter);
    else
        _mandel = new Mandel(_w, _h, _mxit, _sub, _iter);
    _mandel->setNormalBuffer(_normal);
    _mandel->setPrecision((int)mpf_get_prec(_scale));
    _shift_idx = (_w * _sub) * (_sub / 2);
    _cache_valid = false;      // size changed -> phase/relief/normal caches stale
    _require_update = true;
}

void MandelNavigator::ConfigureSampling() {
    bool want_uniform = IsMandelbrot() && (_c_method & ColoringMethod::SUPER_SAMPLING)
                     && (_c_method & ColoringMethod::STRIPE_AVERAGE);
    int wanted_sub = IsMandelbrot() ? (want_uniform ? 2 : _adaptive_sub) : 1;
    if (want_uniform == _uniform_feather && wanted_sub == _sub) return;

    delete[] _iter;
    delete _mandel;
    _uniform_feather = want_uniform;
    _sub = wanted_sub;
    size_t count = (size_t)_w * _h * _sub * _sub;
    _iter = new float[count];
    delete[] _normal;
    _normal = new float[count];
    if (want_uniform)
        _mandel = new Mandel(_w * _sub, _h * _sub, _mxit, 1, _iter);
    else
        _mandel = new Mandel(_w, _h, _mxit, _sub, _iter);
    _mandel->setNormalBuffer(_normal);
    _mandel->setPrecision((int)mpf_get_prec(_scale));
    _shift_idx = (_w * _sub) * (_sub / 2);
    _cache_valid = false;                       // sampling mode changed -> cache stale
}

void MandelNavigator::StartCompute() {
    InterruptCompute();
    _backend->resetCancellation();
    ConfigureSampling();
    _cache_valid = false;                       // fractal changing -> phase cache stale
    // Clear BOTH buffers. _iter marks pixels unresolved; _normal (the DE / normal-map
    // field) must be reset to NaN too, otherwise buildNormalField reads the PREVIOUS
    // view's DE for a pixel whose _iter is freshly filled but whose _normal has not yet
    // been rewritten -- the DE overlay then renders at the old (pre-pan) screen
    // positions for one transient (a ghost), while the smooth base looks correct.
    const float NaN = std::numeric_limits<float>::quiet_NaN();
    const size_t cnt = (size_t)_w * _h * _sub * _sub;
    for (size_t i = 0; i < cnt; ++i) { _iter[i] = EMPTYPIXEL; _normal[i] = NaN; }
    auto compute_task = [this]() {
        this->_mandel->setDensity(color_density);   // for the SAC adaptive-SS detector
        int method = this->_uniform_feather
            ? (this->_c_method & ~ColoringMethod::SUPER_SAMPLING) : this->_c_method;
        ComputeRequest request;
        request.cpuEngine = this->_mandel;
        request.centerRe = this->_z_re;
        request.centerIm = this->_z_im;
        request.scale = this->_scale;
        request.width = this->_uniform_feather ? this->_w * this->_sub : this->_w;
        request.height = this->_uniform_feather ? this->_h * this->_sub : this->_h;
        request.sub = this->_uniform_feather ? 1 : this->_sub;
        request.maxIterations = this->_mxit;
        request.coloringMethod = method;
        request.iterations = this->_iter;
        request.normal = this->_normal;
        if (this->IsExpression()) {
            request.mode = ComputeMode::Expression;
            request.expression = &this->_expressionRuntimeProgram;
            request.expressionFixed = &this->_expressionFixed;
#if defined(MANDEL_ENABLE_ASMJIT)
            request.expressionJit = this->_expressionUseJit
                ? &this->_expressionJit : nullptr;
#endif
            request.expressionPixel = this->_expressionPixel;
            request.expressionBailout = this->_expressionBailout;
            request.expressionColoring =
                this->_expressionProgram.fastIntegerPower() >= 2
                    ? ((method & ColoringMethod::EXTERIOR_DIST_EST)
                        ? formula::ExpressionColoring::Distance
                        : formula::ExpressionColoring::Smooth)
                    : formula::ExpressionColoring::Raw;
        } else if (this->IsJulia()) {
            request.mode = ComputeMode::Julia;
            request.fixedCRe = this->_julia_c_re;
            request.fixedCIm = this->_julia_c_im;
            request.coloringMethod =
                method & ColoringMethod::EXTERIOR_DIST_EST;
        }
        this->_backend->compute(request);
        this->_require_update = true;
    };
    // _task = std::async(&Mandel::Compute, _mandel, _z_re, _z_im, _scale, _mxit, _c_method);
    _task = std::async(std::launch::async, compute_task);
}

void MandelNavigator::InterruptCompute() {
    if (_task.valid()) {
        _backend->cancel();
        _task.wait();
        _backend->resetCancellation();
    }
}

bool MandelNavigator::IsComputing() {
    if (!_task.valid()) return false;
    std::chrono::milliseconds span(0);
    if (_task.wait_for(span) == std::future_status::ready) return false;
    return true;
}

std::string MandelNavigator::GetLocationText() const {
    int precision = (int)mpf_get_prec(_z_re);
    if (IsJulia())
        precision = std::max(precision, (int)std::max(mpf_get_prec(_julia_c_re),
                                                       mpf_get_prec(_julia_c_im)));
    int digits = std::max(18, (int)(precision * log(2.0) / log(10.0)));
    std::vector<char> buf((size_t)digits * 5 + 512);
    if (IsExpression()) {
        gmp_snprintf(buf.data(), buf.size(),
                     "mode: expression\r\nx: %.*Ff\r\ny: %.*Ff\r\nzoom: %.6Fe",
                     digits, _z_re, digits, _z_im, _scale);
        return std::string(buf.data()) + "\r\nformula: " + _expressionProgram.source();
    } else if (IsJulia())
        gmp_snprintf(buf.data(), buf.size(),
                     "mode: julia\r\nx: %.*Ff\r\ny: %.*Ff\r\nzoom: %.6Fe\r\n"
                     "c_re: %.*Ff\r\nc_im: %.*Ff",
                     digits, _z_re, digits, _z_im, _scale,
                     digits, _julia_c_re, digits, _julia_c_im);
    else
        gmp_snprintf(buf.data(), buf.size(), "x: %.*Ff\r\ny: %.*Ff\r\nzoom: %.6Fe",
                     digits, _z_re, digits, _z_im, _scale);
    return std::string(buf.data());
}

bool MandelNavigator::SetLocation(const std::string& xs, const std::string& ys, const std::string& ss) {
    int prec = (int)(std::max(xs.size(), ys.size()) * 3.3219) + 40;
    if (prec < 64) prec = 64;
    mpf_t newScale, newRe, newIm;
    mpf_init2(newScale, prec); mpf_init2(newRe, prec); mpf_init2(newIm, prec);
    bool ok = mpf_set_str(newScale, ss.c_str(), 10) == 0 && mpf_sgn(newScale) > 0 &&
              mpf_set_str(newRe, xs.c_str(), 10) == 0 &&
              mpf_set_str(newIm, ys.c_str(), 10) == 0;
    if (!ok) {
        mpf_clears(newScale, newRe, newIm, (mpf_ptr)0);
        return false;
    }
    if (!IsMandelbrot() && mpf_cmp_d(newScale, 1e12) > 0) mpf_set_d(newScale, 1e12);
    mpf_set_prec(_scale, prec); mpf_set(_scale, newScale);
    mpf_set_prec(_z_re, prec); mpf_set(_z_re, newRe);
    mpf_set_prec(_z_im, prec); mpf_set(_z_im, newIm);
    mpf_set_prec(_t, prec);
    _mandel->setPrecision(prec);
    mpf_clears(newScale, newRe, newIm, (mpf_ptr)0);
    JumpReset();
    return true;
}

void MandelNavigator::UpdateCoords() {
    double effectiveK = _k;
    if (!IsMandelbrot() && _k > 1.0) {
        double scale = mpf_get_d(_scale);
        if (scale > 0.0) effectiveK = std::min(_k, 1e12 / scale);
    }
    double effectiveDisplayDx = _display_dx, effectiveDisplayDy = _display_dy;
    if (effectiveK != _k && _k != 1.0) {
        double zoomX = _dx - (_display_dx - _dx) / (_k - 1.0);
        double zoomY = _dy - (_display_dy - _dy) / (_k - 1.0);
        effectiveDisplayDx = _dx - (zoomX - _dx) * (effectiveK - 1.0);
        effectiveDisplayDy = _dy - (zoomY - _dy) * (effectiveK - 1.0);
    }
    mpf_set_d(_t, effectiveK);
    mpf_mul(_scale, _scale, _t);
    if (IsJulia() && mpf_cmp_d(_scale, 1e12) > 0) mpf_set_d(_scale, 1e12);
    int precision = std::abs(get_exp(_scale)) + 30;
    mpf_set_prec(_scale, precision);
    mpf_set_prec(_z_re, precision);
    mpf_set_prec(_z_im, precision);
    mpf_set_prec(_t, precision);
    _mandel->setPrecision(precision);
    mpf_set_d(_t, 2.0 * (effectiveK - 1.0 + 2.0 * effectiveDisplayDx / _w));
    mpf_div(_t, _t, _scale);
    mpf_sub(_z_re, _z_re, _t);

    mpf_set_d(_t, 2.0 * _h / _w *
                  (effectiveK - 1.0 + 2.0 * effectiveDisplayDy / _h));
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
    // Rebuild the relief / normal-DE field EVERY frame (not only when settled) from the
    // current _iter, so the post-shade layer always tracks the base it is drawn over.
    // Building only on settle left the overlay from the previous location while the base
    // showed the new one -> the two layers visibly slid apart during a pan/zoom redraw.
    if (relief_on) buildReliefHeight();
    if (normal_light_on || de_overlay_on) buildNormalField();
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
                    float bu = (v == EMPTYPIXEL) ? EMPTYPIXEL : colorBaseIndex(v, _c_method);
                    if (settled) _baseUsub[sidx] = bu;
                    if (v == EMPTYPIXEL) continue;
                    if (bu < 0) { ++n; continue; }          // interior -> black (linear 0)
                    int x = ((int)(bu + color_phase)) & (colP - 1);
                    rs += g_palLin[0][x]; gs += g_palLin[1][x]; bs += g_palLin[2][x]; ++n;
                }
                if (!n) continue;
                uint8_t* p = bitmap + ((size_t)i * _w + j) * 3;
                double inv = 1.0 / n;
                p[0] = (uint8_t)(srgbEncode255(rs * inv) + 0.5f);
                p[1] = (uint8_t)(srgbEncode255(gs * inv) + 0.5f);
                p[2] = (uint8_t)(srgbEncode255(bs * inv) + 0.5f);
            }
        }
        if (settled) { _cache_valid = true; _cache_density = color_density; _cache_method = _c_method; }
        applyRelief(bitmap);
        applyNormalLight(bitmap);        applyDEOverlay(bitmap);
        return;
    }
    const int stride = _w * _sub;
    const int c = _sub / 2;
    const bool ss = (_c_method & ColoringMethod::SUPER_SAMPLING) != 0;
    // Cache the settled frame for cheap phase re-colour during animation. Non-SS
    // pixels cache (baseU,width) for colorShadeAA; SS SmoothColor pixels cache
    // per-subpixel base indices in _baseUsub and are marked with SS_SENTINEL.
    const bool cacheable = settled;
    if (cacheable) {
        _baseU.assign((size_t)_w * _h, -1.0f); _widthC.assign((size_t)_w * _h, 0.0f);
        if (ss) { size_t cnt = (size_t)_w * _h * _sub * _sub;
                  if (_baseUsub.size() != cnt) _baseUsub.assign(cnt, EMPTYPIXEL); }
    }
#pragma omp parallel for schedule(dynamic, 8)
    for (int i = 0; i < _h; ++i) {
        for (int j = 0; j < _w; ++j) {
            int idx_bmp = (i * _w + j) * 3;
            int idx = (i * _sub + c) * stride + (j * _sub + c);
            if (_iter[idx] == EMPTYPIXEL) continue;
            // Post-shade overlays (relief / normal-light / DE) were snapshotted from
            // _iter just above; the async compute keeps filling pixels while this
            // base-colour loop runs, so a pixel finished *after* that snapshot would be
            // base-coloured here yet get no overlay -> for one frame it flashes as the
            // bare base (e.g. the smooth gradient with the DE lace missing) before the
            // next frame's snapshot restores the overlay. Keep such a pixel on the warped
            // preview (which carries a consistent base+overlay) until the overlay field
            // includes it. Interior pixels carry no overlay, so colour them as usual.
            if (!settled && _iter[idx] != INTERIOR_SENTINEL) {
                size_t pf = (size_t)i * _w + j;
                if ((de_overlay_on || normal_light_on) && std::isnan(_normalField[pf])) continue;
                if (relief_on && std::isnan(_reliefHt[pf])) continue;
            }
            // In SS mode a flagged pixel has its full sub-block computed (top-left
            // corner subpixel filled): average it. Otherwise this is the 1-sample
            // base layer -> analytic AA from the centre + its 4 pixel neighbours,
            // but ONLY once settled (see above) -- while computing, point-sample.
            if (ss && _iter[idx - (stride + 1) * c] != EMPTYPIXEL) {
                if (cacheable) {
                    cacheSmoothBlock(idx, _c_method);
                    shadeSmoothBlockCached(bitmap + idx_bmp, idx);
                    _baseU[(size_t)i * _w + j] = SS_SENTINEL;
                } else {
                    SmoothColor(bitmap + idx_bmp, idx, _c_method);
                }
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
    applyRelief(bitmap);
    applyNormalLight(bitmap);    applyDEOverlay(bitmap);
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
                    if (bu < 0) { ++n; continue; }         // interior -> black (linear 0)
                    int x = ((int)(bu + color_phase)) & (colP - 1);
                    rs += g_palLin[0][x]; gs += g_palLin[1][x]; bs += g_palLin[2][x]; ++n;
                }
                if (!n) continue;
                uint8_t* p = bitmap + ((size_t)i * _w + j) * 3;
                double inv = 1.0 / n;
                p[0] = (uint8_t)(srgbEncode255(rs * inv) + 0.5f);
                p[1] = (uint8_t)(srgbEncode255(gs * inv) + 0.5f);
                p[2] = (uint8_t)(srgbEncode255(bs * inv) + 0.5f);
            }
        }
        applyRelief(bitmap);
        applyNormalLight(bitmap);        applyDEOverlay(bitmap);
        return;
    }
#pragma omp parallel for schedule(dynamic, 8)
    for (int i = 0; i < _h; ++i) {
        const int stride = _w * _sub, c = _sub / 2;
        for (int j = 0; j < _w; ++j) {
            size_t p = (size_t)i * _w + j;
            uint8_t* q = bitmap + p * 3;
            if (_baseU[p] == SS_SENTINEL) {
                int idx = (i * _sub + c) * stride + (j * _sub + c);
                shadeSmoothBlockCached(q, idx);
            } else {
                colorShadeAA(_baseU[p], _widthC[p], color_phase, q[0], q[1], q[2]);
            }
        }
    }
    applyRelief(bitmap);
    applyNormalLight(bitmap);    applyDEOverlay(bitmap);
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

// Fill _baseUsub for one SmoothColor sub-block from _iter (phase-independent
// analysis, done once when the frame settles).
void MandelNavigator::cacheSmoothBlock(int idx, int method) {
    for (int i = idx - _shift_idx; i <= idx + _shift_idx; i += _w * _sub)
        for (int j = -_sub / 2; j <= _sub / 2; ++j)
            _baseUsub[i + j] = colorBaseIndex(_iter[i + j], method);
}

// Re-shade a SmoothColor sub-block from the cached per-subpixel base indices with
// the live color_phase, averaging in linear light -- matches SmoothColor exactly
// but without any colorFunction (pow/log) work, so it is cheap enough for 60 fps.
void MandelNavigator::shadeSmoothBlockCached(uint8_t* out, int idx) const {
    double rs = 0, gs = 0, bs = 0, ws = 0;
    for (int i = idx - _shift_idx; i <= idx + _shift_idx; i += _w * _sub)
        for (int j = -_sub / 2; j <= _sub / 2; ++j) {
            float bu = _baseUsub[i + j];
            if (bu < 0) { ws += 1; continue; }             // interior/empty -> black (linear 0)
            int x = ((int)(bu + color_phase)) & (colP - 1);  // colP is 2^11 -> mask == %colP (bu,phase >= 0)
            rs += g_palLin[0][x]; gs += g_palLin[1][x]; bs += g_palLin[2][x]; ws += 1;
        }
    double inv = ws > 0.0 ? 1.0 / ws : 0.0;
    out[0] = (uint8_t)(srgbEncode255(rs * inv) + 0.5f);
    out[1] = (uint8_t)(srgbEncode255(gs * inv) + 0.5f);
    out[2] = (uint8_t)(srgbEncode255(bs * inv) + 0.5f);
}

// Build a per-pixel height field (centre subpixel smooth value); interior/empty
// pixels are NaN so they carry no slope and stay unshaded.
void MandelNavigator::buildReliefHeight() {
    const int stride = _w * _sub;
    const int c = _sub / 2;
    const float NaN = std::numeric_limits<float>::quiet_NaN();
    _reliefHt.assign((size_t)_w * _h, 0.0f);
#pragma omp parallel for schedule(static)
    for (int i = 0; i < _h; ++i)
        for (int j = 0; j < _w; ++j) {
            float v = _iter[(size_t)(i * _sub + c) * stride + (j * _sub + c)];
            _reliefHt[(size_t)i * _w + j] =
                (v == EMPTYPIXEL || v == INTERIOR_SENTINEL) ? NaN : (v < 0 ? 0.0f : v);
        }
}

// Post-shade the finished bitmap by Lambert slope lighting derived from the
// screen-space gradient of the height field. Composes with any coloring and the
// phase animation (recomputed each frame from the cached height).
void MandelNavigator::applyRelief(uint8_t* bitmap) {
    if (!relief_on || _reliefHt.size() != (size_t)_w * _h) return;
    applyReliefTo(bitmap, _reliefHt.data(), _w, _h);
}

// Per-pixel analytic normal angle from the engine's _normal (centre subpixel);
// NaN for interior/empty so those stay unshaded.
void MandelNavigator::buildNormalField() {
    const int stride = _w * _sub;
    const int c = _sub / 2;
    const float NaN = std::numeric_limits<float>::quiet_NaN();
    _normalField.assign((size_t)_w * _h, 0.0f);
#pragma omp parallel for schedule(static)
    for (int i = 0; i < _h; ++i)
        for (int j = 0; j < _w; ++j) {
            size_t sidx = (size_t)(i * _sub + c) * stride + (j * _sub + c);
            float v = _iter[sidx];
            _normalField[(size_t)i * _w + j] =
                (v == EMPTYPIXEL || v == INTERIOR_SENTINEL) ? NaN : _normal[sidx];
        }
}

void MandelNavigator::applyNormalLight(uint8_t* bitmap) {
    if (!normal_light_on || _normalField.size() != (size_t)_w * _h) return;
    applyNormalLightTo(bitmap, _normalField.data(), _w, _h);
}

void MandelNavigator::applyDEOverlay(uint8_t* bitmap) {
    if (!de_overlay_on || _normalField.size() != (size_t)_w * _h) return;
    applyDEOverlayTo(bitmap, _normalField.data(), _w, _h);
}

void MandelNavigator::SetMxit(int mxit) {
    _mxit = mxit;
}

void MandelNavigator::GetView(mpf_t re, mpf_t im, mpf_t scale) const {
    mpf_set_prec(re, mpf_get_prec(_z_re)); mpf_set(re, _z_re);
    mpf_set_prec(im, mpf_get_prec(_z_im)); mpf_set(im, _z_im);
    mpf_set_prec(scale, mpf_get_prec(_scale)); mpf_set(scale, _scale);
}

void MandelNavigator::GetJuliaC(mpf_t re, mpf_t im) const {
    mpf_set_prec(re, mpf_get_prec(_julia_c_re)); mpf_set(re, _julia_c_re);
    mpf_set_prec(im, mpf_get_prec(_julia_c_im)); mpf_set(im, _julia_c_im);
}

bool MandelNavigator::SetJuliaC(const std::string& re, const std::string& im) {
    int precision = (int)(std::max(re.size(), im.size()) * 3.3219) + 40;
    precision = std::max(64, precision);
    mpf_t parsedRe, parsedIm;
    mpf_init2(parsedRe, precision); mpf_init2(parsedIm, precision);
    bool ok = mpf_set_str(parsedRe, re.c_str(), 10) == 0 &&
              mpf_set_str(parsedIm, im.c_str(), 10) == 0;
    if (ok) {
        mpf_set_prec(_julia_c_re, precision); mpf_set(_julia_c_re, parsedRe);
        mpf_set_prec(_julia_c_im, precision); mpf_set(_julia_c_im, parsedIm);
    }
    mpf_clear(parsedRe); mpf_clear(parsedIm);
    return ok;
}

void MandelNavigator::SetJuliaMode(bool enabled) {
    if (enabled == IsJulia()) return;
    InterruptCompute();
    if (enabled) {
        if (IsExpression()) RestoreMandelbrotMode();
        SaveMandelbrotState();
        mpf_set(_julia_c_re, _z_re);
        mpf_set(_julia_c_im, _z_im);
        _formula = quadraticJulia();
        _c_method &= ColoringMethod::EXTERIOR_DIST_EST;
        mpf_set_ui(_z_re, 0); mpf_set_ui(_z_im, 0); mpf_set_ui(_scale, 1);
    } else {
        RestoreMandelbrotMode();
        return;
    }
    JumpReset();
    mpf_set_prec(_t, mpf_get_prec(_scale));
    ConfigureSampling();
    _cache_valid = false;
    _require_update = true;
}

void MandelNavigator::SaveMandelbrotState() {
    if (_has_saved_mandel_view) return;
    mp_bitcnt_t precision = std::max({ mpf_get_prec(_z_re), mpf_get_prec(_z_im),
                                       mpf_get_prec(_scale) });
    for (mpf_ptr value : { _saved_mandel_re, _saved_mandel_im, _saved_mandel_scale,
                           _julia_c_re, _julia_c_im })
        mpf_set_prec(value, precision);
    mpf_set(_saved_mandel_re, _z_re);
    mpf_set(_saved_mandel_im, _z_im);
    mpf_set(_saved_mandel_scale, _scale);
    _has_saved_mandel_view = true;
    _saved_mandel_method = _c_method;
}

void MandelNavigator::RestoreMandelbrotMode() {
    InterruptCompute();
#if defined(MANDEL_ENABLE_ASMJIT)
    _expressionUseJit = false;
    _expressionJit.reset();
#endif
    _formula = quadraticMandelbrot();
    if (_has_saved_mandel_view) {
        mpf_set_prec(_z_re, mpf_get_prec(_saved_mandel_re));
        mpf_set_prec(_z_im, mpf_get_prec(_saved_mandel_im));
        mpf_set_prec(_scale, mpf_get_prec(_saved_mandel_scale));
        mpf_set(_z_re, _saved_mandel_re);
        mpf_set(_z_im, _saved_mandel_im);
        mpf_set(_scale, _saved_mandel_scale);
        _c_method = _saved_mandel_method;
    } else {
        mpf_set_d(_z_re, -0.5); mpf_set_ui(_z_im, 0); mpf_set_ui(_scale, 1);
    }
    _has_saved_mandel_view = false;
    JumpReset();
    mpf_set_prec(_t, mpf_get_prec(_scale));
    ConfigureSampling();
    _cache_valid = false;
    _require_update = true;
}

bool MandelNavigator::SetExpressionFormula(
        const std::string& source, FormulaParameter pixel,
        std::complex<double> fixedZ0, std::complex<double> fixedC,
        const std::array<std::complex<double>, 8>& parameters,
        double bailout, formula::ExpressionError* error) {
    formula::ExpressionProgram compiled;
    if (!compiled.compile(source, error))
        return false;
    if ((pixel != FormulaParameter::C &&
         pixel != FormulaParameter::InitialZ) ||
        !(bailout > 0.0) || !std::isfinite(bailout)) {
        if (error) {
            error->position = 0;
            error->message = pixel != FormulaParameter::C &&
                             pixel != FormulaParameter::InitialZ
                ? "unsupported pixel parameter"
                : "bailout must be finite and positive";
        }
        return false;
    }
    formula::ExpressionContext fixed;
    fixed.z0 = fixedZ0;
    fixed.c = fixedC;
    fixed.parameters = parameters;
    formula::ExpressionProgram runtime;
    if (!compiled.specialize(fixed, pixel, runtime, error))
        return false;
#if defined(MANDEL_ENABLE_ASMJIT)
    formula::ExpressionJit4 runtimeJit;
    bool useRuntimeJit =
        runtime.fastPath() == formula::ExpressionProgram::FastPath::None &&
        runtimeJit.compile(runtime);
#endif

    InterruptCompute();
    if (IsJulia()) RestoreMandelbrotMode();
    if (!IsExpression()) SaveMandelbrotState();
    bool planeChanged = IsExpression() && _expressionPixel != pixel;
    _expressionProgram = std::move(compiled);
    _expressionRuntimeProgram = std::move(runtime);
#if defined(MANDEL_ENABLE_ASMJIT)
    _expressionUseJit = useRuntimeJit;
    _expressionJit = std::move(runtimeJit);
#endif
    _expressionFixed = fixed;
    _expressionPixel = pixel;
    _expressionBailout = bailout;
    _formula = expressionFormula();
    _formula.slice.pixel = pixel;
    _c_method = 0;
    if (planeChanged || pixel == FormulaParameter::InitialZ) {
        mpf_set_ui(_z_re, 0); mpf_set_ui(_z_im, 0); mpf_set_ui(_scale, 1);
    } else if (mpf_cmp_d(_scale, 1e12) > 0) {
        mpf_set_d(_scale, 1e12);
    }
    JumpReset();
    ConfigureSampling();
    _cache_valid = false;
    _require_update = true;
    return true;
}

std::string MandelNavigator::GetExpressionAccelerationText() const {
    if (!IsExpression()) return {};
    if (_expressionProgram.fastPath() !=
        formula::ExpressionProgram::FastPath::None) {
        int power = _expressionProgram.fastIntegerPower();
        if (power == 3 && _expressionPixel == FormulaParameter::C) {
            const char* enabled = std::getenv("MANDEL_CUBIC_RESIDUAL");
            const char* residual =
                std::getenv("MANDEL_EXPR_RESIDUAL_POWER");
            const char* series =
                std::getenv("MANDEL_EXPR_CUBIC_SA");
            const char* configured =
                std::getenv("MANDEL_CUBIC_RESIDUAL_SCALE");
            double threshold = configured
                ? std::atof(configured) : 1e8;
            if ((!enabled || std::atoi(enabled) != 0) &&
                (!residual || std::atoi(residual) != 0) &&
                (!series || std::atoi(series) != 0) &&
                std::isfinite(threshold) && threshold > 0.0 &&
                mpf_cmp_d(_scale, threshold) >= 0)
                return "cubic SA / AVX2 adaptive";
        }
        return power <= 8
            ? "integer-power AVX2" : "integer-power scalar";
    }
#if defined(MANDEL_ENABLE_ASMJIT)
    if (_expressionUseJit)
        return "W^X JIT";
    return _expressionRuntimeProgram.avx2Compatible()
        ? "AVX2"
        : (_expressionRuntimeProgram.batchCompatible()
            ? "hybrid AVX2" : "scalar");
#else
    return _expressionRuntimeProgram.avx2Compatible()
        ? "AVX2"
        : (_expressionRuntimeProgram.batchCompatible()
            ? "hybrid AVX2" : "scalar");
#endif
}

int MandelNavigator::GetCMethod() {
    return _c_method;
}

void MandelNavigator::SetCMethod(int c_method) {
    _c_method = !IsMandelbrot()
        ? (c_method & ColoringMethod::EXTERIOR_DIST_EST)
        : c_method;
    if (IsExpression())
        _c_method = ExpressionSupportsDistance()
            ? (c_method & ColoringMethod::EXTERIOR_DIST_EST) : 0;
}

void MandelNavigator::SetRedisplay() {
    _require_update = true;
}
