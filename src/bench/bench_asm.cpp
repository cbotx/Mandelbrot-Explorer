// De-risk: how close is a hand-written mulx addmul_1 to GMP's tuned mpn_addmul_1?
// If within ~1.5x, a high-half kernel built from it (half the work) can beat GMP.
#include <gmp.h>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <chrono>
#include <immintrin.h>
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

static __forceinline void rescaledStep(__m256d& wr, __m256d& wi, __m256d xr, __m256d xi, __m256d scale, __m256d dr, __m256d di) {
    const __m256d factorReal = _mm256_add_pd(_mm256_add_pd(xr, xr), _mm256_mul_pd(scale, wr));
    const __m256d factorImaginary = _mm256_add_pd(_mm256_add_pd(xi, xi), _mm256_mul_pd(scale, wi));
    const __m256d nextReal = _mm256_add_pd(_mm256_sub_pd(_mm256_mul_pd(factorReal, wr), _mm256_mul_pd(factorImaginary, wi)), dr);
    const __m256d nextImaginary = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(factorReal, wi), _mm256_mul_pd(factorImaginary, wr)), di);
    wr = nextReal;
    wi = nextImaginary;
}

static void benchRescaledPipelines() {
    constexpr long iterations = 50000000;
    const __m256d xr = _mm256_set_pd(0.12, -0.08, 0.05, -0.11);
    const __m256d xi = _mm256_set_pd(-0.09, 0.07, 0.1, 0.04);
    const __m256d scale = _mm256_set1_pd(0.0001);
    const __m256d dr = _mm256_set_pd(1e-5, -2e-5, 3e-5, -4e-5);
    const __m256d di = _mm256_set_pd(-3e-5, 4e-5, -1e-5, 2e-5);
    __m256d wr0 = _mm256_set_pd(0.3, -0.2, 0.1, -0.4);
    __m256d wi0 = _mm256_set_pd(-0.15, 0.25, -0.35, 0.05);

    auto start = Clock::now();
    for (long iteration = 0; iteration < iterations; ++iteration) rescaledStep(wr0, wi0, xr, xi, scale, dr, di);
    auto stop = Clock::now();
    const double singleNs = ns_per(start, stop, iterations);

    __m256d wr1 = _mm256_set_pd(-0.22, 0.17, -0.31, 0.28);
    __m256d wi1 = _mm256_set_pd(0.19, -0.27, 0.13, -0.16);
    start = Clock::now();
    for (long iteration = 0; iteration < iterations; ++iteration) {
        rescaledStep(wr0, wi0, xr, xi, scale, dr, di);
        rescaledStep(wr1, wi1, xr, xi, scale, dr, di);
    }
    stop = Clock::now();
    const double dualNs = ns_per(start, stop, iterations * 2);

    alignas(32) double sink[4];
    _mm256_store_pd(sink, _mm256_add_pd(_mm256_add_pd(wr0, wi0), _mm256_add_pd(wr1, wi1)));
    printf("\n== AVX2 rescaled-step software pipeline ==\n");
    printf("  single-chain %.3f ns/group | dual-chain %.3f ns/group | throughput %.2fx | sink=%.17g\n", singleNs, dualNs, singleNs / dualNs, sink[0] + sink[1] + sink[2] + sink[3]);
}

struct LookupEntry {
    double radiusSquared;
    int skip;
};

struct LookupCoefficient {
    double ar, ai, br, bi;
};

struct FlatLookup {
    LookupEntry entry;
    LookupCoefficient coefficient;
};

static void benchBlaLayouts() {
    constexpr int referenceLength = 16384;
    constexpr int levels = 14;
    std::vector<std::vector<LookupEntry>> entries(levels);
    std::vector<std::vector<LookupCoefficient>> coefficients(levels);
    for (int level = 0; level < levels; ++level) {
        const int count = (referenceLength - 1) >> level;
        entries[level].resize(count);
        coefficients[level].resize(count);
        for (int index = 0; index < count; ++index) {
            entries[level][index] = {std::ldexp(1.0 + (index & 7) * 0.03125, -2 * level), 1 << level};
            coefficients[level][index] = {1.0 + level * 0.001, level * 0.0001, 0.5 + index * 1e-8, -0.25};
        }
    }

    std::vector<uint32_t> offsets(referenceLength + 1);
    std::vector<FlatLookup> flat;
    flat.reserve(referenceLength * 2);
    for (int reference = 1; reference < referenceLength; ++reference) {
        offsets[reference] = static_cast<uint32_t>(flat.size());
        int startLevel = reference == 1 ? levels - 1 : static_cast<int>(_tzcnt_u32(static_cast<unsigned>(reference - 1)));
        if (startLevel >= levels) startLevel = levels - 1;
        for (int level = startLevel; level >= 0; --level) {
            const int index = (reference - 1) >> level;
            if (index < static_cast<int>(entries[level].size())) flat.push_back({entries[level][index], coefficients[level][index]});
        }
    }
    offsets[referenceLength] = static_cast<uint32_t>(flat.size());

    constexpr long iterations = 20000000;
    volatile double sink = 0.0;
    auto start = Clock::now();
    for (long iteration = 0; iteration < iterations; ++iteration) {
        const int reference = 1 + static_cast<int>((static_cast<uint64_t>(iteration) * 11400714819323198485ull) % (referenceLength - 2));
        const double magnitude = std::ldexp(0.75 + (iteration & 7) * 0.03125, -2 * (iteration % 8));
        int startLevel = reference == 1 ? levels - 1 : static_cast<int>(_tzcnt_u32(static_cast<unsigned>(reference - 1)));
        if (startLevel >= levels) startLevel = levels - 1;
        for (int level = startLevel; level >= 0; --level) {
            const int index = (reference - 1) >> level;
            if (index >= static_cast<int>(entries[level].size())) continue;
            const LookupEntry& entry = entries[level][index];
            if (magnitude >= entry.radiusSquared) continue;
            const LookupCoefficient& coefficient = coefficients[level][index];
            sink = sink + coefficient.ar + coefficient.br + entry.skip;
            break;
        }
    }
    auto stop = Clock::now();
    const double hierarchyNs = ns_per(start, stop, iterations);

    start = Clock::now();
    for (long iteration = 0; iteration < iterations; ++iteration) {
        const int reference = 1 + static_cast<int>((static_cast<uint64_t>(iteration) * 11400714819323198485ull) % (referenceLength - 2));
        const double magnitude = std::ldexp(0.75 + (iteration & 7) * 0.03125, -2 * (iteration % 8));
        for (uint32_t index = offsets[reference]; index < offsets[reference + 1]; ++index) {
            const FlatLookup& candidate = flat[index];
            if (magnitude >= candidate.entry.radiusSquared) continue;
            sink = sink + candidate.coefficient.ar + candidate.coefficient.br + candidate.entry.skip;
            break;
        }
    }
    stop = Clock::now();
    const double flatNs = ns_per(start, stop, iterations);

    printf("\n== BLA candidate layout lookup ==\n");
    printf("  hierarchy %.3f ns/query | flat %.3f ns/query | speedup %.2fx | candidates=%zu | sink=%.17g\n", hierarchyNs, flatNs, hierarchyNs / flatNs, flat.size(), sink);
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
    benchRescaledPipelines();
    benchBlaLayouts();
    return 0;
}
