#ifndef __FLOAT_MATH_H__
#define __FLOAT_MATH_H__

#include <gmp.h>

double mpf_get_ld(mpf_t a);

int get_exp(long double ld);

int get_exp(mpf_t a);

void bit_print(long double ld);

#endif