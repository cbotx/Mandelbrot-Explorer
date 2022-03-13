#include <mpir.h>

#include "mandel_perturbation.h"
#include "dll_interface.h"

static Mandel* mandel;

void __stdcall __halt() {
    mandel->SetHalt(true);
}

void __stdcall __core_init(int width, int height, int mxit, bool hq, float* const iter) {
    int sub = hq ? 3 : 1;
    mandel = new Mandel(width, height, mxit, sub, iter);
}

void __stdcall __mandel_compute(int nx, const char xc[],
                                int ny, const char yc[],
                                int nz, const char zc[]) {
    mandel->SetHalt(false);
    int precision = static_cast<int>(nz * log(10) / log(2)) + 30;

    mpf_t c_re, c_im, scale;
    mpf_set_default_prec(precision);
    mandel->setPrecision(precision);
    mpf_init_set_str(c_re, xc, 10);
    mpf_init_set_str(c_im, yc, 10);
    mpf_init_set_str(scale, zc, 10);
    mandel->Compute(c_re, c_im, scale, 0);
    mpf_clear(c_re);
    mpf_clear(c_im);
    mpf_clear(scale);
}