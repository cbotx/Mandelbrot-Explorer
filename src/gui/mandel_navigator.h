#ifndef __MANDEL_NAVIGATOR_H__
#define __MANDEL_NAVIGATOR_H__


#include <gmp.h>
#include <atomic>
#include <future>
#include <string>
#include <vector>

#include "navigator.h"
#include "mandel_perturbation.h"
#include "formula_spec.h"
#include "float_math.h"

class MandelNavigator : public Navigator {
private:
    Mandel* _mandel;
    float* _iter;
    mpf_t _z_re, _z_im, _scale;
    mpf_t _t;
    int _mxit;
    int _sub;
    int _adaptive_sub;
    bool _uniform_feather = false;
    int _c_method = ColoringMethod::EXTERIOR_DIST_EST;
    int _shift_idx;
    std::atomic_bool _require_update{true};
    bool _need_settle = false;   // a computing frame was point-sampled; force one AA pass once settled

    // ---- palette-phase animation cache -------------------------------------
    // Phase-independent per-pixel colouring analysis, rebuilt whenever the frame
    // settles. RecolorPhase re-shades these with the live color_phase each frame
    // without re-running colorFunction / the neighbourhood gradient.
    std::vector<float> _baseU;      // AA path: palette-index centre (_w*_h)
    std::vector<float> _widthC;     // AA path: palette footprint width (_w*_h)
    std::vector<float> _baseUsub;   // Feather path: per-subpixel base index (sub grid)
    bool _cache_valid = false;
    float _cache_density = -1.0f;   // color_density the cache was built with
    int _cache_method = -1;         // _c_method the cache was built with

    // ---- relief (screen-space slope) height cache --------------------------
    // Per-pixel smooth height (centre subpixel), rebuilt on each settled frame;
    // applyRelief() re-shades from it each frame so the light stays animatable.
    std::vector<float> _reliefHt;   // (_w*_h), NaN for interior/empty pixels
    // Analytic normal-map: the engine fills _normal (sub-grid) with the per-pixel
    // surface-normal angle; _normalField is the per-output-pixel angle (centre
    // subpixel), NaN for interior/empty. applyNormalLight() shades from it.
    float* _normal = nullptr;
    std::vector<float> _normalField;

    std::future<void> _task;

public:
    MandelNavigator(int width, int height, int sub, int max_iteration, double zoom_step, double zoom_time);

    virtual ~MandelNavigator();

    void Reset();

    // Re-target the render buffers to a new pixel size (native-resolution display:
    // the fractal is computed at the window's view size, no upscaling). Rebuilds
    // _iter/_normal/_mandel and invalidates the phase cache; the caller re-renders.
    void Resize(int width, int height);

    void StartCompute();

    void InterruptCompute();

    void UpdateCoords();

    void UpdateBitmap(uint8_t* bitmap);

    // Cheap palette-phase re-colour for animation: re-shades the cached settled
    // frame with the live color_phase (no fractal recompute, no colorFunction).
    // Falls back to a full UpdateBitmap while a render is in progress or if the
    // cache is stale (density/method changed).
    void RecolorPhase(uint8_t* bitmap);

    void SetMxit(int mxit);
    int GetMxit() const { return _mxit; }
    // Copy the current view: center (re/im) and scale into caller-owned mpf_t.
    void GetView(mpf_t re, mpf_t im, mpf_t scale) const;
    mp_bitcnt_t GetViewPrecision() const { return mpf_get_prec(_scale); }
    static constexpr FormulaContext GetFormulaContext() { return quadraticMandelbrot(); }
    
    int GetCMethod();

    void SetCMethod(int c_method);

    void SetRedisplay();

    bool IsComputing();

    std::string GetLocationText() const;

    // Jump to an absolute location: x/y as decimal strings, scale as a plain
    // decimal string (scientific already expanded). Sets precision from the digit
    // count, resets the preview transform. Returns false on a parse error.
    bool SetLocation(const std::string& x, const std::string& y, const std::string& scale);

private:
    void ConfigureSampling();
    void SmoothColor(uint8_t* bitmap_pixel, int idx, int _c_method);
    // SS phase cache: fill per-subpixel base palette indices for a SmoothColor
    // sub-block (settle), and re-shade a sub-block from those cached indices with
    // the live phase (animation) -- byte-identical to SmoothColor but pow/log-free.
    void cacheSmoothBlock(int idx, int method);
    void shadeSmoothBlockCached(uint8_t* out, int idx) const;
    // Build _reliefHt from the settled _iter field (centre subpixel per pixel).
    void buildReliefHeight();
    // Multiply the finished bitmap by per-pixel Lambert slope shading.
    void applyRelief(uint8_t* bitmap);
    // Build _normalField from the engine's _normal (centre subpixel per pixel) and
    // apply the analytic normal-map Lambert shade to the finished bitmap.
    void buildNormalField();
    void applyNormalLight(uint8_t* bitmap);
    void applyDEOverlay(uint8_t* bitmap);
};


#endif