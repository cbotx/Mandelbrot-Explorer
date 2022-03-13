#ifndef __INTERPOLATE_H__
#define __INTERPOLATE_H__

void mono_cubic_interpolate(const float* xs, const float* ys, const int length, float* mp, const int mp_sz);

void cos_interpolate(const float* xs, const float* ys, const int length, float* mp, const int mp_sz);

#endif