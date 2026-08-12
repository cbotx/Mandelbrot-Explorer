// Standalone A/B for the custom BigFixed backend vs GMP mpf_t:
//   (1) correctness of mul / add / sub over random O(1) values, and
//   (2) throughput of a raw multiply and of a complex-square reference-orbit
//       loop (the actual Mandelbrot hot path), at the limb counts that matter.
//
// Answers directly: is the custom fixed-point arithmetic correct, and is it
// faster than GMP, at what precision.

#include <gmp.h>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <vector>

#include "bigfixed.h"

using Clock = std::chrono::high_resolution_clock;
static double secs(Clock::time_point a, Clock::time_point b) {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count() / 1e9;
}

static uint64_t rnd64() {
    // 64 random bits from rand() (good enough for test vectors).
    uint64_t x = 0;
    for (int i = 0; i < 4; ++i) x = (x << 16) ^ (uint64_t)(rand() & 0xffff);
    return x;
}

static __declspec(noinline) FloatExp bfToFeBaseline(const BigFixed& a) {
    return bf_to_fe(a);
}

static __declspec(noinline) FloatExp bfToFeNoFrexp(const BigFixed& a) {
    if (a.sign == 0) return FloatExp{0.0, 0};
    int top = a.L - 1;
    while (top >= 0 && a.m[top] == 0) --top;
    if (top < 0) return FloatExp{0.0, 0};
    const uint64_t hi = a.m[top];
    const uint64_t lo = top >= 1 ? a.m[top - 1] : 0;
    const double value = static_cast<double>(hi) * 18446744073709551616.0 + static_cast<double>(lo);
    uint64_t bits = 0;
    std::memcpy(&bits, &value, sizeof(bits));
    const int exponent = static_cast<int>((bits >> 52) & 0x7ff) - 1022;
    bits = (bits & ((uint64_t{1} << 52) - 1)) | (uint64_t{1022} << 52);
    double mantissa = 0.0;
    std::memcpy(&mantissa, &bits, sizeof(mantissa));
    FloatExp result{mantissa, static_cast<int64_t>(exponent) + static_cast<int64_t>(64) * (top - a.L)};
    return a.sign < 0 ? fe_neg(result) : result;
}

static void benchToFe(int L, long iterations) {
    constexpr size_t count = 256;
    std::vector<BigFixed> values;
    values.reserve(count);
    for (size_t index = 0; index < count; ++index) {
        values.emplace_back(L);
        BigFixed& value = values.back();
        value.sign = (index & 1) ? 1 : -1;
        const int leadingZeroLimbs = static_cast<int>(index % static_cast<size_t>(std::min(L, 8)));
        const int top = L - 1 - leadingZeroLimbs;
        for (int limb = 0; limb <= top; ++limb) value.m[limb] = rnd64();
        if (value.m[top] == 0) value.m[top] = 1;
    }

    int mismatches = 0;
    for (const BigFixed& value : values) {
        const FloatExp baseline = bfToFeBaseline(value);
        const FloatExp candidate = bfToFeNoFrexp(value);
        uint64_t baselineBits = 0;
        uint64_t candidateBits = 0;
        std::memcpy(&baselineBits, &baseline.m, sizeof(baselineBits));
        std::memcpy(&candidateBits, &candidate.m, sizeof(candidateBits));
        mismatches += baselineBits != candidateBits || baseline.e != candidate.e;
    }

    volatile double mantissaSink = 0.0;
    volatile int64_t exponentSink = 0;
    auto start = Clock::now();
    for (long iteration = 0; iteration < iterations; ++iteration) {
        const FloatExp result = bfToFeBaseline(values[static_cast<size_t>(iteration) & (count - 1)]);
        mantissaSink = mantissaSink + result.m;
        exponentSink = exponentSink ^ result.e;
    }
    auto stop = Clock::now();
    const double baselineNs = secs(start, stop) / iterations * 1e9;

    start = Clock::now();
    for (long iteration = 0; iteration < iterations; ++iteration) {
        const FloatExp result = bfToFeNoFrexp(values[static_cast<size_t>(iteration) & (count - 1)]);
        mantissaSink = mantissaSink + result.m;
        exponentSink = exponentSink ^ result.e;
    }
    stop = Clock::now();
    const double candidateNs = secs(start, stop) / iterations * 1e9;

    printf("  to_fe L=%2d: frexp %6.1f ns | bits %6.1f ns | speedup %.2fx | mismatch=%d | sink=%.3g/%lld\n", L, baselineNs, candidateNs, baselineNs / candidateNs, mismatches, mantissaSink, (long long)exponentSink);
}

// Random BigFixed with |value| < 4 (integer part 0..3), full random fraction.
static void randBF(BigFixed& a, int L) {
    a.setL(L);
    a.sign = (rand() & 1) ? 1 : -1;
    for (int i = 0; i < L - 1; ++i) a.m[i] = rnd64();
    a.m[L - 1] = (uint64_t)(rand() & 3);
    int nz = 0;
    for (int i = 0; i < L; ++i) nz |= (a.m[i] != 0);
    if (!nz) a.sign = 0;
}

static __forceinline void bfAddSubReady(BigFixed& result, const BigFixed& left, int leftSign, const BigFixed& right, int rightSign) {
    const int L = left.L;
    if (leftSign == 0) {
        std::memcpy(result.m.data(), right.m.data(), sizeof(uint64_t) * static_cast<size_t>(L));
        result.sign = rightSign;
        return;
    }
    if (rightSign == 0) {
        std::memcpy(result.m.data(), left.m.data(), sizeof(uint64_t) * static_cast<size_t>(L));
        result.sign = leftSign;
        return;
    }
    if (leftSign == rightSign) {
        bf_mag_add(result.m.data(), left.m.data(), right.m.data(), L);
        result.sign = leftSign;
        return;
    }
    const int comparison = bf_mag_cmp(left.m.data(), right.m.data(), L);
    if (comparison == 0) {
        result.sign = 0;
        std::memset(result.m.data(), 0, sizeof(uint64_t) * static_cast<size_t>(L));
    } else if (comparison > 0) {
        bf_mag_sub(result.m.data(), left.m.data(), right.m.data(), L);
        result.sign = leftSign;
    } else {
        bf_mag_sub(result.m.data(), right.m.data(), left.m.data(), L);
        result.sign = rightSign;
    }
}

static __forceinline void bfAddReady(BigFixed& result, const BigFixed& left, const BigFixed& right) {
    bfAddSubReady(result, left, left.sign, right, right.sign);
}

static __forceinline void bfSubReady(BigFixed& result, const BigFixed& left, const BigFixed& right) {
    bfAddSubReady(result, left, left.sign, right, -right.sign);
}

static __forceinline void bfSumDiffReady(BigFixed& sum, BigFixed& difference, const BigFixed& left, const BigFixed& right) {
    const int L = left.L;
    if (left.sign == 0) {
        std::memcpy(sum.m.data(), right.m.data(), sizeof(uint64_t) * static_cast<size_t>(L));
        std::memcpy(difference.m.data(), right.m.data(), sizeof(uint64_t) * static_cast<size_t>(L));
        sum.sign = right.sign;
        difference.sign = -right.sign;
        return;
    }
    if (right.sign == 0) {
        std::memcpy(sum.m.data(), left.m.data(), sizeof(uint64_t) * static_cast<size_t>(L));
        std::memcpy(difference.m.data(), left.m.data(), sizeof(uint64_t) * static_cast<size_t>(L));
        sum.sign = left.sign;
        difference.sign = left.sign;
        return;
    }

    const int comparison = bf_mag_cmp(left.m.data(), right.m.data(), L);
    const uint64_t* larger = comparison >= 0 ? left.m.data() : right.m.data();
    const uint64_t* smaller = comparison >= 0 ? right.m.data() : left.m.data();
    uint64_t* addOutput = left.sign == right.sign ? sum.m.data() : difference.m.data();
    uint64_t* subtractOutput = left.sign == right.sign ? difference.m.data() : sum.m.data();
    unsigned char carry = 0;
    unsigned char borrow = 0;
    for (int limb = 0; limb < L; ++limb) {
        carry = _addcarry_u64(carry, left.m[limb], right.m[limb], &addOutput[limb]);
        borrow = _subborrow_u64(borrow, larger[limb], smaller[limb], &subtractOutput[limb]);
    }

    if (left.sign == right.sign) {
        sum.sign = left.sign;
        difference.sign = comparison == 0 ? 0 : (comparison > 0 ? left.sign : -left.sign);
    } else {
        difference.sign = left.sign;
        sum.sign = comparison == 0 ? 0 : (comparison > 0 ? left.sign : right.sign);
    }
}

static __forceinline void bfMulReady(BigFixed& result, const BigFixed& left, const BigFixed& right, uint64_t* scratch) {
    if (left.sign == 0 || right.sign == 0) {
        result.sign = 0;
        std::memset(result.m.data(), 0, sizeof(uint64_t) * static_cast<size_t>(left.L));
        return;
    }
    bf_mag_mulshift(result.m.data(), left.m.data(), right.m.data(), left.L, scratch);
    result.sign = left.sign == right.sign ? 1 : -1;
}

static __declspec(noinline) void bfComplexSquareAddPrepared(BigFixed& nextReal, BigFixed& nextImaginary, const BigFixed& real, const BigFixed& imaginary, const BigFixed& constantReal, const BigFixed& constantImaginary, BigFixed& sum, BigFixed& difference, BigFixed& product, BigFixed& squareReal, BigFixed& squareImaginary, uint64_t* scratch) {
    bfAddReady(sum, real, imaginary);
    bfSubReady(difference, real, imaginary);
    bfMulReady(product, real, imaginary, scratch);
    bfMulReady(squareReal, sum, difference, scratch);
    bfAddReady(squareImaginary, product, product);
    bfAddReady(nextReal, squareReal, constantReal);
    bfAddReady(nextImaginary, squareImaginary, constantImaginary);
}

static void bfMagnitudeSquareHigh(uint64_t* output, const uint64_t* input, int L, int guard) {
    const int lowColumn = L - 1 - guard;
    const int outputLimbs = 2 * L - lowColumn;
    std::memset(output, 0, sizeof(uint64_t) * static_cast<size_t>(outputLimbs + 2));
    for (int left = 0; left < L - 1; ++left) {
        int rightStart = left + 1;
        if (left + rightStart < lowColumn) rightStart = lowColumn - left;
        if (rightStart < left + 1) rightStart = left + 1;
        const int count = L - rightStart;
        if (count <= 0) continue;
        const int target = left + rightStart - lowColumn;
        mp_limb_t carry = mpn_addmul_1(reinterpret_cast<mp_limb_t*>(output + target), reinterpret_cast<const mp_limb_t*>(input + rightStart), static_cast<mp_size_t>(count), static_cast<mp_limb_t>(input[left]));
        int carryIndex = target + count;
        while (carry && carryIndex <= outputLimbs) {
            const unsigned char nextCarry = _addcarry_u64(0, output[carryIndex], carry, &output[carryIndex]);
            carry = nextCarry;
            ++carryIndex;
        }
    }
    mpn_lshift(reinterpret_cast<mp_limb_t*>(output), reinterpret_cast<const mp_limb_t*>(output), static_cast<mp_size_t>(outputLimbs + 1), 1);
    for (int index = 0; index < L; ++index) {
        const int column = 2 * index;
        if (column + 1 < lowColumn) continue;
        uint64_t high = 0;
        const uint64_t low = _mulx_u64(input[index], input[index], &high);
        if (column >= lowColumn) {
            const int target = column - lowColumn;
            unsigned char carry = _addcarry_u64(0, output[target], low, &output[target]);
            carry = _addcarry_u64(carry, output[target + 1], high, &output[target + 1]);
            for (int carryIndex = target + 2; carry && carryIndex <= outputLimbs; ++carryIndex) carry = _addcarry_u64(carry, output[carryIndex], 0, &output[carryIndex]);
        } else {
            unsigned char carry = _addcarry_u64(0, output[0], high, &output[0]);
            for (int carryIndex = 1; carry && carryIndex <= outputLimbs; ++carryIndex) carry = _addcarry_u64(carry, output[carryIndex], 0, &output[carryIndex]);
        }
    }
}

static __forceinline void bfSquareHighReady(BigFixed& result, const BigFixed& input, uint64_t* scratch) {
    constexpr int guard = 3;
    if (input.sign == 0) {
        result.sign = 0;
        std::memset(result.m.data(), 0, sizeof(uint64_t) * static_cast<size_t>(input.L));
        return;
    }
    bfMagnitudeSquareHigh(scratch, input.m.data(), input.L, guard);
    std::memcpy(result.m.data(), scratch + guard, sizeof(uint64_t) * static_cast<size_t>(input.L));
    if (scratch[guard - 1] >> 63) {
        unsigned char carry = _addcarry_u64(0, result.m[0], 1, &result.m[0]);
        for (int limb = 1; limb < input.L && carry; ++limb) carry = _addcarry_u64(carry, result.m[limb], 0, &result.m[limb]);
    }
    result.sign = 1;
}

static void benchReadyStep(int L, long iterations) {
    constexpr size_t inputCount = 16;
    std::vector<BigFixed> realInputs;
    std::vector<BigFixed> imaginaryInputs;
    realInputs.reserve(inputCount);
    imaginaryInputs.reserve(inputCount);
    for (size_t index = 0; index < inputCount; ++index) {
        realInputs.emplace_back(L);
        imaginaryInputs.emplace_back(L);
        randBF(realInputs.back(), L);
        randBF(imaginaryInputs.back(), L);
        realInputs.back().m[L - 1] &= 1;
        imaginaryInputs.back().m[L - 1] &= 1;
    }

    BigFixed cReal(L), cImaginary(L);
    randBF(cReal, L);
    randBF(cImaginary, L);
    cReal.m[L - 1] = 0;
    cImaginary.m[L - 1] = 0;

    std::vector<uint64_t> scratch(2 * L);
    BigFixed sum(L), difference(L), product(L), realPart(L), imaginaryPart(L), nextReal(L), nextImaginary(L);
    BigFixed readySum(L), readyDifference(L), readyProduct(L), readyRealPart(L), readyImaginaryPart(L), readyNextReal(L), readyNextImaginary(L);
    BigFixed fusedSum(L), fusedDifference(L), fusedProduct(L), fusedRealPart(L), fusedImaginaryPart(L), fusedNextReal(L), fusedNextImaginary(L);
    BigFixed integratedSum(L), integratedDifference(L), integratedProduct(L), integratedRealPart(L), integratedImaginaryPart(L), integratedNextReal(L), integratedNextImaginary(L);

    int mismatches = 0;
    for (size_t index = 0; index < inputCount; ++index) {
        const BigFixed& real = realInputs[index];
        const BigFixed& imaginary = imaginaryInputs[index];
        bf_add(sum, real, imaginary);
        bf_sub(difference, real, imaginary);
        bf_mul(product, real, imaginary, scratch.data());
        bf_mul(realPart, sum, difference, scratch.data());
        bf_add(imaginaryPart, product, product);
        bf_add(nextReal, realPart, cReal);
        bf_add(nextImaginary, imaginaryPart, cImaginary);

        bfAddReady(readySum, real, imaginary);
        bfSubReady(readyDifference, real, imaginary);
        bfMulReady(readyProduct, real, imaginary, scratch.data());
        bfMulReady(readyRealPart, readySum, readyDifference, scratch.data());
        bfAddReady(readyImaginaryPart, readyProduct, readyProduct);
        bfAddReady(readyNextReal, readyRealPart, cReal);
        bfAddReady(readyNextImaginary, readyImaginaryPart, cImaginary);

        bfSumDiffReady(fusedSum, fusedDifference, real, imaginary);
        bfMulReady(fusedProduct, real, imaginary, scratch.data());
        bfMulReady(fusedRealPart, fusedSum, fusedDifference, scratch.data());
        bfAddReady(fusedImaginaryPart, fusedProduct, fusedProduct);
        bfAddReady(fusedNextReal, fusedRealPart, cReal);
        bfAddReady(fusedNextImaginary, fusedImaginaryPart, cImaginary);

        bfComplexSquareAddPrepared(integratedNextReal, integratedNextImaginary, real, imaginary, cReal, cImaginary, integratedSum, integratedDifference, integratedProduct, integratedRealPart, integratedImaginaryPart, scratch.data());

        mismatches += nextReal.sign != readyNextReal.sign || nextImaginary.sign != readyNextImaginary.sign || nextReal.m != readyNextReal.m || nextImaginary.m != readyNextImaginary.m || nextReal.sign != fusedNextReal.sign || nextImaginary.sign != fusedNextImaginary.sign || nextReal.m != fusedNextReal.m || nextImaginary.m != fusedNextImaginary.m || nextReal.sign != integratedNextReal.sign || nextImaginary.sign != integratedNextImaginary.sign || nextReal.m != integratedNextReal.m || nextImaginary.m != integratedNextImaginary.m;
    }

    volatile uint64_t sink = 0;
    auto start = Clock::now();
    for (long iteration = 0; iteration < iterations; ++iteration) {
        const size_t index = static_cast<size_t>(iteration) & (inputCount - 1);
        const BigFixed& real = realInputs[index];
        const BigFixed& imaginary = imaginaryInputs[index];
        bf_add(sum, real, imaginary);
        bf_sub(difference, real, imaginary);
        bf_mul(product, real, imaginary, scratch.data());
        bf_mul(realPart, sum, difference, scratch.data());
        bf_add(imaginaryPart, product, product);
        bf_add(nextReal, realPart, cReal);
        bf_add(nextImaginary, imaginaryPart, cImaginary);
        sink ^= nextReal.m[0] ^ nextImaginary.m[0];
    }
    auto stop = Clock::now();
    const double baselineNs = secs(start, stop) / iterations * 1e9;

    start = Clock::now();
    for (long iteration = 0; iteration < iterations; ++iteration) {
        const size_t index = static_cast<size_t>(iteration) & (inputCount - 1);
        const BigFixed& real = realInputs[index];
        const BigFixed& imaginary = imaginaryInputs[index];
        bfAddReady(readySum, real, imaginary);
        bfSubReady(readyDifference, real, imaginary);
        bfMulReady(readyProduct, real, imaginary, scratch.data());
        bfMulReady(readyRealPart, readySum, readyDifference, scratch.data());
        bfAddReady(readyImaginaryPart, readyProduct, readyProduct);
        bfAddReady(readyNextReal, readyRealPart, cReal);
        bfAddReady(readyNextImaginary, readyImaginaryPart, cImaginary);
        sink ^= readyNextReal.m[0] ^ readyNextImaginary.m[0];
    }
    stop = Clock::now();
    const double readyNs = secs(start, stop) / iterations * 1e9;

    start = Clock::now();
    for (long iteration = 0; iteration < iterations; ++iteration) {
        const size_t index = static_cast<size_t>(iteration) & (inputCount - 1);
        const BigFixed& real = realInputs[index];
        const BigFixed& imaginary = imaginaryInputs[index];
        bfSumDiffReady(fusedSum, fusedDifference, real, imaginary);
        bfMulReady(fusedProduct, real, imaginary, scratch.data());
        bfMulReady(fusedRealPart, fusedSum, fusedDifference, scratch.data());
        bfAddReady(fusedImaginaryPart, fusedProduct, fusedProduct);
        bfAddReady(fusedNextReal, fusedRealPart, cReal);
        bfAddReady(fusedNextImaginary, fusedImaginaryPart, cImaginary);
        sink ^= fusedNextReal.m[0] ^ fusedNextImaginary.m[0];
    }
    stop = Clock::now();
    const double fusedNs = secs(start, stop) / iterations * 1e9;

    start = Clock::now();
    for (long iteration = 0; iteration < iterations; ++iteration) {
        const size_t index = static_cast<size_t>(iteration) & (inputCount - 1);
        const BigFixed& real = realInputs[index];
        const BigFixed& imaginary = imaginaryInputs[index];
        bfComplexSquareAddPrepared(integratedNextReal, integratedNextImaginary, real, imaginary, cReal, cImaginary, integratedSum, integratedDifference, integratedProduct, integratedRealPart, integratedImaginaryPart, scratch.data());
        sink ^= integratedNextReal.m[0] ^ integratedNextImaginary.m[0];
    }
    stop = Clock::now();
    const double integratedNs = secs(start, stop) / iterations * 1e9;

    printf("  ready L=%2d: baseline %7.1f ns | ready %7.1f ns (%.2fx) | sumdiff %7.1f ns (%.2fx) | integrated %7.1f ns (%.2fx) | mismatch=%d | sink=%llx\n", L, baselineNs, readyNs, baselineNs / readyNs, fusedNs, baselineNs / fusedNs, integratedNs, baselineNs / integratedNs, mismatches, (unsigned long long)sink);
}

static void benchHighSquareStep(int L, long iterations) {
    constexpr size_t inputCount = 16;
    std::vector<BigFixed> realInputs;
    std::vector<BigFixed> imaginaryInputs;
    realInputs.reserve(inputCount);
    imaginaryInputs.reserve(inputCount);
    for (size_t index = 0; index < inputCount; ++index) {
        realInputs.emplace_back(L);
        imaginaryInputs.emplace_back(L);
        randBF(realInputs.back(), L);
        randBF(imaginaryInputs.back(), L);
        realInputs.back().m[L - 1] &= 1;
        imaginaryInputs.back().m[L - 1] &= 1;
    }

    BigFixed cReal(L), cImaginary(L);
    randBF(cReal, L);
    randBF(cImaginary, L);
    cReal.m[L - 1] = 0;
    cImaginary.m[L - 1] = 0;

    std::vector<uint64_t> scratch(2 * L + 2);
    BigFixed sum(L), difference(L), product(L), realPart(L), imaginaryPart(L), nextReal(L), nextImaginary(L);
    BigFixed realSquare(L), imaginarySquare(L), sumSquare(L), squareTemporary(L), squareRealPart(L), squareImaginaryPart(L), squareNextReal(L), squareNextImaginary(L);

    int mismatches = 0;
    for (size_t index = 0; index < inputCount; ++index) {
        const BigFixed& real = realInputs[index];
        const BigFixed& imaginary = imaginaryInputs[index];
        bfAddReady(sum, real, imaginary);
        bfSubReady(difference, real, imaginary);
        bfMulReady(product, real, imaginary, scratch.data());
        bfMulReady(realPart, sum, difference, scratch.data());
        bfAddReady(imaginaryPart, product, product);
        bfAddReady(nextReal, realPart, cReal);
        bfAddReady(nextImaginary, imaginaryPart, cImaginary);

        bfSquareHighReady(realSquare, real, scratch.data());
        bfSquareHighReady(imaginarySquare, imaginary, scratch.data());
        bfAddReady(sum, real, imaginary);
        bfSquareHighReady(sumSquare, sum, scratch.data());
        bfSubReady(squareRealPart, realSquare, imaginarySquare);
        bfSubReady(squareTemporary, sumSquare, realSquare);
        bfSubReady(squareImaginaryPart, squareTemporary, imaginarySquare);
        bfAddReady(squareNextReal, squareRealPart, cReal);
        bfAddReady(squareNextImaginary, squareImaginaryPart, cImaginary);

        mismatches += nextReal.sign != squareNextReal.sign || nextImaginary.sign != squareNextImaginary.sign || nextReal.m != squareNextReal.m || nextImaginary.m != squareNextImaginary.m;
    }

    volatile uint64_t sink = 0;
    auto start = Clock::now();
    for (long iteration = 0; iteration < iterations; ++iteration) {
        const size_t index = static_cast<size_t>(iteration) & (inputCount - 1);
        const BigFixed& real = realInputs[index];
        const BigFixed& imaginary = imaginaryInputs[index];
        bfAddReady(sum, real, imaginary);
        bfSubReady(difference, real, imaginary);
        bfMulReady(product, real, imaginary, scratch.data());
        bfMulReady(realPart, sum, difference, scratch.data());
        bfAddReady(imaginaryPart, product, product);
        bfAddReady(nextReal, realPart, cReal);
        bfAddReady(nextImaginary, imaginaryPart, cImaginary);
        sink ^= nextReal.m[0] ^ nextImaginary.m[0];
    }
    auto stop = Clock::now();
    const double multiplyNs = secs(start, stop) / iterations * 1e9;

    start = Clock::now();
    for (long iteration = 0; iteration < iterations; ++iteration) {
        const size_t index = static_cast<size_t>(iteration) & (inputCount - 1);
        const BigFixed& real = realInputs[index];
        const BigFixed& imaginary = imaginaryInputs[index];
        bfSquareHighReady(realSquare, real, scratch.data());
        bfSquareHighReady(imaginarySquare, imaginary, scratch.data());
        bfAddReady(sum, real, imaginary);
        bfSquareHighReady(sumSquare, sum, scratch.data());
        bfSubReady(squareRealPart, realSquare, imaginarySquare);
        bfSubReady(squareTemporary, sumSquare, realSquare);
        bfSubReady(squareImaginaryPart, squareTemporary, imaginarySquare);
        bfAddReady(squareNextReal, squareRealPart, cReal);
        bfAddReady(squareNextImaginary, squareImaginaryPart, cImaginary);
        sink ^= squareNextReal.m[0] ^ squareNextImaginary.m[0];
    }
    stop = Clock::now();
    const double squareNs = secs(start, stop) / iterations * 1e9;

    printf("  high-square L=%3d: 2mul %8.1f ns | 3sqr %8.1f ns | speedup %.2fx | exact-mismatch=%d | sink=%llx\n", L, multiplyNs, squareNs, multiplyNs / squareNs, mismatches, (unsigned long long)sink);
}

// |mpf a - mpf b| as a double.
static double mpfAbsDiff(const mpf_t a, const mpf_t b) {
    mpf_t d;
    mpf_init2(d, mpf_get_prec(a));
    mpf_sub(d, a, b);
    mpf_abs(d, d);
    double r = mpf_get_d(d);
    mpf_clear(d);
    return r;
}

static void correctness(int L, int trials) {
    const mp_bitcnt_t hp = (mp_bitcnt_t)64 * (L + 3);
    double ulp = std::ldexp(1.0, -64 * (L - 1));
    double worst_mul = 0, worst_add = 0, worst_sub = 0;

    mpf_t ma, mb, mref, mgot;
    mpf_init2(ma, hp);
    mpf_init2(mb, hp);
    mpf_init2(mref, hp);
    mpf_init2(mgot, hp);
    std::vector<uint64_t> tmp(2 * L);
    BigFixed a(L), b(L), r(L);

    for (int t = 0; t < trials; ++t) {
        randBF(a, L);
        randBF(b, L);
        bf_to_mpf(ma, a);
        bf_to_mpf(mb, b);

        bf_mul(r, a, b, tmp.data());
        bf_to_mpf(mgot, r);
        mpf_mul(mref, ma, mb);
        worst_mul = std::max(worst_mul, mpfAbsDiff(mref, mgot) / ulp);

        bf_add(r, a, b);
        bf_to_mpf(mgot, r);
        mpf_add(mref, ma, mb);
        worst_add = std::max(worst_add, mpfAbsDiff(mref, mgot) / ulp);

        bf_sub(r, a, b);
        bf_to_mpf(mgot, r);
        mpf_sub(mref, ma, mb);
        worst_sub = std::max(worst_sub, mpfAbsDiff(mref, mgot) / ulp);
    }
    printf("  correctness (%d trials): mul %.2f ulp  add %.2f ulp  sub %.2f ulp  %s\n", trials, worst_mul, worst_add, worst_sub, (worst_mul <= 1.001 && worst_add < 0.001 && worst_sub < 0.001) ? "OK" : "FAIL");
    mpf_clear(ma);
    mpf_clear(mb);
    mpf_clear(mref);
    mpf_clear(mgot);
}

static void benchMul(int L, long iters) {
    std::vector<uint64_t> tmp(2 * L);
    BigFixed a(L), b(L), r(L);
    randBF(a, L);
    randBF(b, L);
    a.sign = b.sign = 1;

    // BigFixed: chain r=a*b; a=r to force a dependency.
    auto t0 = Clock::now();
    for (long i = 0; i < iters; ++i) {
        bf_mul(r, a, b, tmp.data());
        a.m.swap(r.m);
    }
    auto t1 = Clock::now();
    double bf_ns = secs(t0, t1) / iters * 1e9;

    mpf_t ma, mb, mr;
    mpf_init2(ma, (mp_bitcnt_t)64 * L);
    mpf_init2(mb, (mp_bitcnt_t)64 * L);
    mpf_init2(mr, (mp_bitcnt_t)64 * L);
    bf_to_mpf(ma, a);
    bf_to_mpf(mb, b);
    t0 = Clock::now();
    for (long i = 0; i < iters; ++i) {
        mpf_mul(mr, ma, mb);
        mpf_set(ma, mr);
    }
    t1 = Clock::now();
    double gmp_ns = secs(t0, t1) / iters * 1e9;
    mpf_clear(ma);
    mpf_clear(mb);
    mpf_clear(mr);

    printf("  mul  L=%2d (%4d bit): BigFixed %7.1f ns   GMP %7.1f ns   speedup %.2fx\n", L, 64 * L, bf_ns, gmp_ns, gmp_ns / bf_ns);
}

// Complex square z = z^2 + c, the reference-orbit inner step.
static void benchOrbit(int L, long iters) {
    std::vector<uint64_t> tmp(2 * L);
    BigFixed zr(L), zi(L), cr(L), ci(L), t1v(L), t2v(L), t3v(L);
    zr.setInt(0);
    zi.setInt(0);
    // c = -0.5 + 0i  (stays bounded forever).
    cr.setL(L);
    cr.sign = -1;
    cr.m[L - 1] = 0;
    cr.m[L - 2] = 0x8000000000000000ull; // 0.5
    ci.setInt(0);

    auto t0 = Clock::now();
    for (long i = 0; i < iters; ++i) {
        bf_mul(t1v, zr, zr, tmp.data()); // zr^2
        bf_mul(t2v, zi, zi, tmp.data()); // zi^2
        bf_mul(t3v, zr, zi, tmp.data()); // zr*zi
        bf_sub(t1v, t1v, t2v);           // zr^2 - zi^2
        bf_add(zr, t1v, cr);             // zr' = zr^2 - zi^2 + cr
        bf_add(t3v, t3v, t3v);           // 2 zr zi
        bf_add(zi, t3v, ci);             // zi' = 2 zr zi + ci
    }
    auto t1 = Clock::now();
    double bf_ns = secs(t0, t1) / iters * 1e9;

    mpf_t mzr, mzi, mcr, mci, a, bb, cc;
    mp_bitcnt_t p = (mp_bitcnt_t)64 * L;
    mpf_init2(mzr, p);
    mpf_init2(mzi, p);
    mpf_init2(mcr, p);
    mpf_init2(mci, p);
    mpf_init2(a, p);
    mpf_init2(bb, p);
    mpf_init2(cc, p);
    mpf_set_ui(mzr, 0);
    mpf_set_ui(mzi, 0);
    mpf_set_d(mcr, -0.5);
    mpf_set_ui(mci, 0);
    t0 = Clock::now();
    for (long i = 0; i < iters; ++i) {
        mpf_mul(a, mzr, mzr);
        mpf_mul(bb, mzi, mzi);
        mpf_mul(cc, mzr, mzi);
        mpf_sub(a, a, bb);
        mpf_add(mzr, a, mcr);
        mpf_add(cc, cc, cc);
        mpf_add(mzi, cc, mci);
    }
    t1 = Clock::now();
    double gmp_ns = secs(t0, t1) / iters * 1e9;
    mpf_clear(mzr);
    mpf_clear(mzi);
    mpf_clear(mcr);
    mpf_clear(mci);
    mpf_clear(a);
    mpf_clear(bb);
    mpf_clear(cc);

    printf("  orbit L=%2d (%4d bit): BigFixed %7.1f ns   GMP %7.1f ns   speedup %.2fx\n", L, 64 * L, bf_ns, gmp_ns, gmp_ns / bf_ns);
}

// Engine-faithful complex-square step comparison at deep precision:
//   (K) 2-mul Karatsuba : re=(a+b)(a-b), ab=a*b        -> what the engine uses now
//   (S) 3-square form   : a^2, b^2, (a+b)^2 via mpn_sqr -> 2ab=(a+b)^2-a^2-b^2
// vs mpf 2-mul Karatsuba. Also checks (S) and (K) agree to rounding.
static void benchOrbitSqr(int L, long iters) {
    std::vector<uint64_t> tmp(2 * L);
    auto mkC = [&](BigFixed& cr, BigFixed& ci) {
        cr.setL(L);
        cr.sign = -1;
        cr.m[L - 2] = 0x8000000000000000ull;
        ci.setInt(0); // c=-0.5
    };

    // (K) 2-mul Karatsuba
    {
        BigFixed zr(L), zi(L), cr(L), ci(L), s(L), d(L), ab(L), re(L), im(L);
        mkC(cr, ci);
        auto t0 = Clock::now();
        for (long i = 0; i < iters; ++i) {
            bf_add(s, zr, zi);
            bf_sub(d, zr, zi);
            bf_mul(re, s, d, tmp.data());   // a^2 - b^2
            bf_mul(ab, zr, zi, tmp.data()); // a*b
            bf_add(zr, re, cr);
            bf_add(im, ab, ab);
            bf_add(zi, im, ci);
        }
        auto t1 = Clock::now();
        double k_ns = secs(t0, t1) / iters * 1e9;

        // (S) 3-square form
        BigFixed Zr(L), Zi(L), s1(L), s2(L), s3(L), apb(L), t(L);
        mkC(cr, ci);
        t0 = Clock::now();
        for (long i = 0; i < iters; ++i) {
            bf_sqr(s1, Zr, tmp.data()); // a^2
            bf_sqr(s2, Zi, tmp.data()); // b^2
            bf_add(apb, Zr, Zi);
            bf_sqr(s3, apb, tmp.data()); // (a+b)^2
            bf_sub(t, s1, s2);
            bf_add(Zr, t, cr); // a^2-b^2+cr
            bf_sub(t, s3, s1);
            bf_sub(t, t, s2);
            bf_add(Zi, t, ci); // 2ab+ci
        }
        t1 = Clock::now();
        double s_ns = secs(t0, t1) / iters * 1e9;

        // mpf 2-mul Karatsuba
        mpf_t zr2, zi2, cr2, ci2, ms, md, mre, mab;
        mp_bitcnt_t p = (mp_bitcnt_t)64 * L;
        mpf_init2(zr2, p);
        mpf_init2(zi2, p);
        mpf_init2(cr2, p);
        mpf_init2(ci2, p);
        mpf_init2(ms, p);
        mpf_init2(md, p);
        mpf_init2(mre, p);
        mpf_init2(mab, p);
        mpf_set_ui(zr2, 0);
        mpf_set_ui(zi2, 0);
        mpf_set_d(cr2, -0.5);
        mpf_set_ui(ci2, 0);
        t0 = Clock::now();
        for (long i = 0; i < iters; ++i) {
            mpf_add(ms, zr2, zi2);
            mpf_sub(md, zr2, zi2);
            mpf_mul(mre, ms, md);
            mpf_mul(mab, zr2, zi2);
            mpf_add(zr2, mre, cr2);
            mpf_add(mab, mab, mab);
            mpf_add(zi2, mab, ci2);
        }
        t1 = Clock::now();
        double m_ns = secs(t0, t1) / iters * 1e9;
        mpf_clear(zr2);
        mpf_clear(zi2);
        mpf_clear(cr2);
        mpf_clear(ci2);
        mpf_clear(ms);
        mpf_clear(md);
        mpf_clear(mre);
        mpf_clear(mab);

        printf("  step L=%2d (%4d bit): mpf(2mul) %7.1f | K(2mul) %7.1f (%.2fx) | S(3sqr) %7.1f (%.2fx)\n", L, 64 * L, m_ns, k_ns, m_ns / k_ns, s_ns, m_ns / s_ns);
    }
}

// One-step agreement: 3-square result must equal 2-mul result to <=1 ulp.
static void checkSqrForm(int L) {
    std::vector<uint64_t> tmp(2 * L);
    BigFixed a(L), b(L);
    randBF(a, L);
    randBF(b, L);
    a.sign = b.sign = 1;
    // K: re=(a+b)(a-b), im=2ab
    BigFixed s(L), d(L), reK(L), abK(L), imK(L);
    bf_add(s, a, b);
    bf_sub(d, a, b);
    bf_mul(reK, s, d, tmp.data());
    bf_mul(abK, a, b, tmp.data());
    bf_add(imK, abK, abK);
    // S: a^2,b^2,(a+b)^2
    BigFixed s1(L), s2(L), s3(L), apb(L), reS(L), imS(L), t(L);
    bf_sqr(s1, a, tmp.data());
    bf_sqr(s2, b, tmp.data());
    bf_add(apb, a, b);
    bf_sqr(s3, apb, tmp.data());
    bf_sub(reS, s1, s2);
    bf_sub(t, s3, s1);
    bf_sub(imS, t, s2);
    mpf_t x, y, dd;
    mp_bitcnt_t p = (mp_bitcnt_t)64 * L;
    mpf_init2(x, p);
    mpf_init2(y, p);
    mpf_init2(dd, p);
    // ulp = 2^-(64*(L-1)); express |diff| in ulp by shifting up before converting
    // (a plain ldexp underflows for large L).
    bf_to_mpf(x, reK);
    bf_to_mpf(y, reS);
    mpf_sub(dd, x, y);
    mpf_abs(dd, dd);
    mpf_mul_2exp(dd, dd, (mp_bitcnt_t)64 * (L - 1));
    double reDiff = mpf_get_d(dd);
    bf_to_mpf(x, imK);
    bf_to_mpf(y, imS);
    mpf_sub(dd, x, y);
    mpf_abs(dd, dd);
    mpf_mul_2exp(dd, dd, (mp_bitcnt_t)64 * (L - 1));
    double imDiff = mpf_get_d(dd);
    mpf_clear(x);
    mpf_clear(y);
    mpf_clear(dd);
    printf("  agree L=%2d: re %.2f ulp  im %.2f ulp  %s\n", L, reDiff, imDiff, (reDiff <= 2.001 && imDiff <= 2.001) ? "OK" : "FAIL");
}

int main() {
    srand(12345);
    int Ls[] = {3, 9, 24, 48}; // ~1e30, 1e150, 1e450, 1e900 zoom depths
    printf("== correctness (BigFixed vs GMP mpf) ==\n");
    for (int L : Ls) {
        printf(" L=%d:\n", L);
        correctness(L, 20000);
    }

    printf("\n== raw multiply throughput ==\n");
    for (int L : Ls) benchMul(L, 2000000);

    printf("\n== complex-square reference-orbit step (naive 3-mul) ==\n");
    for (int L : Ls) benchOrbit(L, 1000000);

    printf("\n== 3-square vs 2-mul agreement (ulp) ==\n");
    for (int L : Ls) checkSqrForm(L);

    printf("\n== engine-faithful complex-square step: 2-mul vs 3-square ==\n");
    for (int L : Ls) benchOrbitSqr(L, 1000000);

    printf("\n== BigFixed to FloatExp conversion ==\n");
    for (int L : Ls) benchToFe(L, 5000000);

    printf("\n== pre-sized reference step wrappers ==\n");
    for (int L : Ls) benchReadyStep(L, 1000000);

    printf("\n== production-step high-half square crossover ==\n");
    for (int L : {48, 64, 96, 128}) benchHighSquareStep(L, 200000);
    return 0;
}
