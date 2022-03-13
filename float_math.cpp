#include <iostream>
#include <bitset>
#include <intrin.h>
#include <mpir.h>

#include "float_math.h"

long double mpf_get_ld(mpf_t a) {
    int exp = a->_mp_exp * 64;
    ildouble idb{ 1 };
    if (a->_mp_size < 0) {
        idb.i.sgn = 1;
    }
    for (int i = abs(a->_mp_size) - 1; i >= 0; --i) {
        int leading_zeros = __lzcnt64(*(a->_mp_d + i));
        exp -= leading_zeros;
        if (leading_zeros < 64) {
            idb.i.man = (*(a->_mp_d + i) << leading_zeros + 1 >> 1);
            if (leading_zeros > 0 && i > 0) {
                idb.i.man |= (*(a->_mp_d + i - 1) >> 64 - leading_zeros);
            }
            --exp;
            if (exp + idb.i.exp <= 0) return 0;
            idb.i.exp += exp;
            return idb.d;
        }
    }
    return 0;
}

int get_exp(long double ld) {
    ildouble idb{ ld };
    return idb.i.exp - 16383;
}

int get_exp(mpf_t a) {
    int exp = a->_mp_exp * 64;
    for (int i = abs(a->_mp_size) - 1; i >= 0; --i) {
        int leading_zeros = __lzcnt64(*(a->_mp_d + i));
        exp -= leading_zeros;
        if (leading_zeros < 64) break;
    }
    return exp;
}

void bit_print(long double ld) {
    std::bitset<64> bs = *(reinterpret_cast<unsigned long long*>(&ld) + 1);
    std::cout << bs << ' ';
    bs = *reinterpret_cast<unsigned long long*>(&ld);
    std::cout << bs << '\n';
}