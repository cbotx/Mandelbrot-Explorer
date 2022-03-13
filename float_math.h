#ifndef __FLOAT_MATH_H__
#define __FLOAT_MATH_H__

#include <mpir.h>


union idouble {
    double d;
    struct {
        unsigned long long man : 52;
        unsigned long long exp : 11;
        unsigned long long sgn : 1;
    } i;
};

union ildouble {
    long double d;
    struct {
        unsigned long long man : 63;
        unsigned long long exp : 15;
        unsigned long long sgn : 1;
    } i;
};

long double mpf_get_ld(mpf_t a);

int get_exp(long double ld);

int get_exp(mpf_t a);

void bit_print(long double ld);

#endif