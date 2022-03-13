#ifndef __MANDEL_NAVIGATOR_H__
#define __MANDEL_NAVIGATOR_H__


#include <mpir.h>
#include <future>

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
    int _c_method = ColoringMethod::EXTERIOR_DIST_EST;
    int _shift_idx;
    bool _require_update;

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
    
    int GetCMethod();

    void SetCMethod(int c_method);

    void SetRedisplay();

    bool IsComputing();

private:
    void SmoothColor(uint8_t* bitmap_pixel, int idx, int _c_method);
};


#endif