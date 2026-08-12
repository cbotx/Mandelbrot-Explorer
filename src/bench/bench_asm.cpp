// De-risk: how close is a hand-written mulx addmul_1 to GMP's tuned mpn_addmul_1?
// If within ~1.5x, a high-half kernel built from it (half the work) can beat GMP.
#include <gmp.h>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <chrono>
#include <vector>

extern "C" uint64_t mm_addmul_1(uint64_t* rp, const uint64_t* up, int64_t n, uint64_t v);

using Clock = std::chrono::high_resolution_clock;
static double ns_per(Clock::time_point a, Clock::time_point b, long n) {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count() / (double)n;
}
static uint64_t rnd64() {
    uint64_t x = 0;
    for (int i = 0; i < 4; ++i) x = (x << 16) ^ (uint64_t)(rand() & 0xffff);
    return x;
}

int main() {
    srand(2024);
    printf("== mm_addmul_1 (hand asm) vs mpn_addmul_1 (GMP) ==\n");
    for (int n : {8, 16, 24, 32, 48, 64}) {
        std::vector<uint64_t> up(n), r1(n + 1), r2(n + 1);
        for (int i = 0; i < n; ++i) {
            up[i] = rnd64();
            uint64_t x = rnd64();
            r1[i] = x;
            r2[i] = x;
        }
        uint64_t v = rnd64();
        // correctness
        uint64_t c1 = mpn_addmul_1((mp_limb_t*)r1.data(), (const mp_limb_t*)up.data(), n, v);
        uint64_t c2 = mm_addmul_1(r2.data(), up.data(), n, v);
        bool ok = (c1 == c2) && (memcmp(r1.data(), r2.data(), n * 8) == 0);

        const long IT = 5000000;
        // reset each time is costly; just chain in-place (values grow but timing is representative)
        auto t0 = Clock::now();
        for (long k = 0; k < IT; ++k) { v ^= mpn_addmul_1((mp_limb_t*)r1.data(), (const mp_limb_t*)up.data(), n, v); }
        auto t1 = Clock::now();
        double gmp = ns_per(t0, t1, IT);
        t0 = Clock::now();
        for (long k = 0; k < IT; ++k) { v ^= mm_addmul_1(r2.data(), up.data(), n, v); }
        t1 = Clock::now();
        double mine = ns_per(t0, t1, IT);

        printf("  n=%2d: GMP %6.1f ns | mine %6.1f ns | ratio %.2fx  %s\n", n, gmp, mine, gmp / mine, ok ? "correct" : "WRONG");
    }
    return 0;
}
