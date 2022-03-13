#include <cstdint>

#include "utils.h"


uint32_t burtle_hash(uint32_t a) {
    a = (a+0x7ed55d16) + (a<<12);
    a = (a^0xc761c23c) ^ (a>>19);
    a = (a+0x165667b1) + (a<<5);
    a = (a+0xd3a2646c) ^ (a<<9);
    a = (a+0xfd7046c5) + (a<<3);
    a = (a^0xb55a4f09) ^ (a>>16);
    return a;
}

// [-0.5, 0.5)
double jitter(uint32_t x, uint32_t y, uint32_t c) {
  return burtle_hash(x + burtle_hash(y + burtle_hash(c))) / (double) (0x100000000LL) - 0.5;
}

