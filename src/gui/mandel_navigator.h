#ifndef __MANDEL_NAVIGATOR_H__
#define __MANDEL_NAVIGATOR_H__


#include <gmp.h>
#include <atomic>
#include <future>
#include <string>

#include "navigator.h"
#include "mandel_perturbation.h"
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

    std::future<void> _task;

public:
    MandelNavigator(int width, int height, int sub, int max_iteration, double zoom_step, double zoom_time);

    virtual ~MandelNavigator();

    void Reset();

    void StartCompute();

    void InterruptCompute();

    void UpdateCoords();

    void UpdateBitmap(uint8_t* bitmap);

    void SetMxit(int mxit);
    int GetMxit() const { return _mxit; }
    // Copy the current view: center (re/im) and scale into caller-owned mpf_t.
    void GetView(mpf_t re, mpf_t im, mpf_t scale) const;
    
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
};


#endif