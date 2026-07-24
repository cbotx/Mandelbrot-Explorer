#include <iostream>
#include <cmath>
#include <gmp.h>

#include "float_math.h"

// Portable conversion from an mpf_t to the low-precision Float type used by the
// perturbation kernel. Uses only the public GMP API (mpf_get_d_2exp) instead of
// reaching into GMP''s limb layout, so it is correct under MSVC (where
// long double == double) as well as toolchains with 80-bit long double.
//
// Note: the mantissa carried through is that of a double (53 bits). On targets
// where Float is a true 80-bit long double this drops ~11 bits of the seed
// value; the perturbation deltas dominate and rebasing absorbs the difference.
double mpf_get_ld(mpf_t a) {
    signed long int e;
    double d = mpf_get_d_2exp(&e, a);   // a = d * 2^e, |d| in [0.5, 1)
    return ldexp(d, (int)e);
}

// Unbiased base-2 exponent of a long double: floor(log2|ld|). Zero maps to a
// large negative sentinel so series-approximation magnitude checks treat it as
// negligible; the sentinel keeps exponent differences within int range.
int get_exp(long double ld) {
    if (ld == 0.0L) return -0x40000000;
    int e;
    (void)frexpl(ld, &e);               // ld = m * 2^e, |m| in [0.5, 1)
    return e - 1;                       // hardware unbiased exponent
}

// Base-2 exponent of an mpf_t: matches the original routine (floor(log2|a|)+1).
int get_exp(mpf_t a) {
    if (mpf_sgn(a) == 0) return 0;
    signed long int e;
    (void)mpf_get_d_2exp(&e, a);        // a = d * 2^e, |d| in [0.5, 1)
    return (int)e;
}

void bit_print(long double ld) {
    int e;
    long double m = frexpl(ld, &e);
    std::cout << m << " x 2^" << e << '\n';
}
