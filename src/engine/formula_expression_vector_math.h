#ifndef MANDEL_FORMULA_EXPRESSION_VECTOR_MATH_H
#define MANDEL_FORMULA_EXPRESSION_VECTOR_MATH_H

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <immintrin.h>

namespace formula {
namespace detail {

template <size_t N>
inline __m256d vectorMathEvaluatePolynomial(__m256d value, const double (&coefficients)[N]) {
    __m256d result = _mm256_set1_pd(coefficients[N - 1]);
    for (size_t index = N - 1; index-- > 0;) result = _mm256_add_pd(_mm256_set1_pd(coefficients[index]), _mm256_mul_pd(value, result));
    return result;
}

inline void vectorMathRealSinCos(__m256d value, __m256d& sine, __m256d& cosine, int (&quadrants)[4]) {
    constexpr double InverseHalfPi = 0.63661977236758134308;
    constexpr double HalfPiHigh = 1.57079632673412561417;
    constexpr double HalfPiLow = 6.07710050650619224932e-11;
    constexpr double SineCoefficients[] = {
        -1.66666666666666657415e-1,
        8.33333333333333321769e-3,
        -1.98412698412698412530e-4,
        2.75573192239858925112e-6,
        -2.50521083854417187751e-8,
        1.60590438368216133409e-10,
        -7.64716373181981640551e-13,
        2.81145725434552059811e-15,
        -8.22063524662432949554e-18,
    };
    constexpr double CosineCoefficients[] = {
        -5.00000000000000000000e-1,
        4.16666666666666643537e-2,
        -1.38888888888888894189e-3,
        2.48015873015873015658e-5,
        -2.75573192239858925112e-7,
        2.08767569878680989792e-9,
        -1.14707455977297245073e-11,
        4.77947733238738525345e-14,
        -1.56192069685862252711e-16,
        4.11031762331216484407e-19,
    };

    const __m256d absolute = _mm256_andnot_pd(_mm256_set1_pd(-0.0), value);
    __m256d reduced;
    if (_mm256_movemask_pd(_mm256_cmp_pd(absolute, _mm256_set1_pd(0.78539816339744830962), _CMP_LE_OQ)) == 15) {
        quadrants[0] = quadrants[1] = quadrants[2] = quadrants[3] = 0;
        reduced = value;
    } else {
        const __m256d small = _mm256_cmp_pd(absolute, _mm256_set1_pd(0.78539816339744830962), _CMP_LE_OQ);
        __m256d multiple = _mm256_round_pd(_mm256_mul_pd(value, _mm256_set1_pd(InverseHalfPi)), _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
        multiple = _mm256_blendv_pd(multiple, _mm256_setzero_pd(), small);
        _mm_storeu_si128(reinterpret_cast<__m128i*>(quadrants), _mm256_cvttpd_epi32(multiple));
        reduced = _mm256_sub_pd(value, _mm256_mul_pd(multiple, _mm256_set1_pd(HalfPiHigh)));
        reduced = _mm256_sub_pd(reduced, _mm256_mul_pd(multiple, _mm256_set1_pd(HalfPiLow)));
    }
    const __m256d squared = _mm256_mul_pd(reduced, reduced);
    sine = _mm256_add_pd(reduced, _mm256_mul_pd(_mm256_mul_pd(reduced, squared), vectorMathEvaluatePolynomial(squared, SineCoefficients)));
    cosine = _mm256_add_pd(_mm256_set1_pd(1.0), _mm256_mul_pd(squared, vectorMathEvaluatePolynomial(squared, CosineCoefficients)));
}

inline __m256d vectorMathPositiveExp(__m256d value) {
    constexpr double InverseLn2 = 1.44269504088896338700;
    constexpr double Ln2High = 6.93147180369123816490e-1;
    constexpr double Ln2Low = 1.90821492927058770002e-10;
    constexpr double ExpCoefficients[] = {
        1.66666666666666019037e-1,
        -2.77777777770155933842e-3,
        6.61375632143793436117e-5,
        -1.65339022054652515390e-6,
        4.13813679705723846039e-8,
    };

    const __m256d multiple = _mm256_round_pd(_mm256_mul_pd(value, _mm256_set1_pd(InverseLn2)), _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
    const __m256d high = _mm256_sub_pd(value, _mm256_mul_pd(multiple, _mm256_set1_pd(Ln2High)));
    const __m256d low = _mm256_mul_pd(multiple, _mm256_set1_pd(Ln2Low));
    const __m256d reduced = _mm256_sub_pd(high, low);
    const __m256d squared = _mm256_mul_pd(reduced, reduced);
    const __m256d correction = _mm256_sub_pd(reduced, _mm256_mul_pd(squared, vectorMathEvaluatePolynomial(squared, ExpCoefficients)));
    const __m256d fraction = _mm256_div_pd(_mm256_mul_pd(reduced, correction), _mm256_sub_pd(_mm256_set1_pd(2.0), correction));
    const __m256d polynomial = _mm256_sub_pd(_mm256_set1_pd(1.0), _mm256_sub_pd(_mm256_sub_pd(low, fraction), high));

    alignas(16) int exponents[4];
    _mm_store_si128(reinterpret_cast<__m128i*>(exponents), _mm256_cvttpd_epi32(multiple));
    alignas(32) double scales[4];
    for (int lane = 0; lane < 4; ++lane) {
        const uint64_t bits = static_cast<uint64_t>(exponents[lane] + 1023) << 52;
        std::memcpy(&scales[lane], &bits, sizeof(bits));
    }
    return _mm256_mul_pd(polynomial, _mm256_load_pd(scales));
}

inline void vectorMathRealSinhCosh(__m256d value, __m256d& hyperbolicSine, __m256d& hyperbolicCosine) {
    constexpr double SinhCoefficients[] = {
        1.66666666666666657415e-1,
        8.33333333333333321769e-3,
        1.98412698412698412530e-4,
        2.75573192239858925112e-6,
        2.50521083854417187751e-8,
        1.60590438368216133409e-10,
        7.64716373181981640551e-13,
        2.81145725434552059811e-15,
    };
    constexpr double CoshCoefficients[] = {
        5.00000000000000000000e-1,
        4.16666666666666643537e-2,
        1.38888888888888894189e-3,
        2.48015873015873015658e-5,
        2.75573192239858925112e-7,
        2.08767569878680989792e-9,
        1.14707455977297245073e-11,
        4.77947733238738525345e-14,
        1.56192069685862252711e-16,
    };

    const __m256d signMask = _mm256_set1_pd(-0.0);
    const __m256d absolute = _mm256_andnot_pd(signMask, value);
    const __m256d squared = _mm256_mul_pd(value, value);
    const __m256d smallSine = _mm256_add_pd(value, _mm256_mul_pd(_mm256_mul_pd(value, squared), vectorMathEvaluatePolynomial(squared, SinhCoefficients)));
    const __m256d smallCosine = _mm256_add_pd(_mm256_set1_pd(1.0), _mm256_mul_pd(squared, vectorMathEvaluatePolynomial(squared, CoshCoefficients)));
    const __m256d small = _mm256_cmp_pd(absolute, _mm256_set1_pd(0.5), _CMP_LE_OQ);
    if (_mm256_movemask_pd(small) == 15) {
        hyperbolicSine = smallSine;
        hyperbolicCosine = smallCosine;
        return;
    }

    const __m256d exponential = vectorMathPositiveExp(absolute);
    const __m256d reciprocal = _mm256_div_pd(_mm256_set1_pd(1.0), exponential);
    const __m256d half = _mm256_set1_pd(0.5);
    __m256d largeSine = _mm256_mul_pd(half, _mm256_sub_pd(exponential, reciprocal));
    largeSine = _mm256_xor_pd(largeSine, _mm256_and_pd(value, signMask));
    const __m256d largeCosine = _mm256_mul_pd(half, _mm256_add_pd(exponential, reciprocal));
    hyperbolicSine = _mm256_blendv_pd(largeSine, smallSine, small);
    hyperbolicCosine = _mm256_blendv_pd(largeCosine, smallCosine, small);
}

inline uint8_t evaluateVectorComplexSinCosPairFast(__m256d inputReal, __m256d inputImaginary, __m256d& sineReal, __m256d& sineImaginary, __m256d& cosineReal, __m256d& cosineImaginary) {
    const __m256d signMask = _mm256_set1_pd(-0.0);
    const __m256d absoluteReal = _mm256_andnot_pd(signMask, inputReal);
    const __m256d absoluteImaginary = _mm256_andnot_pd(signMask, inputImaginary);
    __m256d safe = _mm256_cmp_pd(absoluteReal, _mm256_set1_pd(128.0), _CMP_LE_OQ);
    safe = _mm256_and_pd(safe, _mm256_cmp_pd(absoluteImaginary, _mm256_set1_pd(20.0), _CMP_LE_OQ));
    safe = _mm256_and_pd(safe, _mm256_cmp_pd(inputReal, _mm256_setzero_pd(), _CMP_NEQ_OQ));
    safe = _mm256_and_pd(safe, _mm256_cmp_pd(inputImaginary, _mm256_setzero_pd(), _CMP_NEQ_OQ));
    const uint8_t safeMask = static_cast<uint8_t>(_mm256_movemask_pd(safe));
    if (safeMask == 0) {
        sineReal = sineImaginary = cosineReal = cosineImaginary = _mm256_setzero_pd();
        return 0;
    }

    const __m256d evaluationReal = _mm256_and_pd(inputReal, safe);
    const __m256d evaluationImaginary = _mm256_and_pd(inputImaginary, safe);
    __m256d realSine, realCosine, imaginarySine, imaginaryCosine;
    int quadrants[4];
    alignas(32) double reducedSine[4], reducedCosine[4];
    vectorMathRealSinCos(evaluationReal, realSine, realCosine, quadrants);
    _mm256_store_pd(reducedSine, realSine);
    _mm256_store_pd(reducedCosine, realCosine);
    vectorMathRealSinhCosh(evaluationImaginary, imaginarySine, imaginaryCosine);
    for (int lane = 0; lane < 4; ++lane) {
        double laneSine = reducedSine[lane];
        double laneCosine = reducedCosine[lane];
        switch (quadrants[lane] & 3) {
        case 1:
            std::swap(laneSine, laneCosine);
            laneCosine = -laneCosine;
            break;
        case 2:
            laneSine = -laneSine;
            laneCosine = -laneCosine;
            break;
        case 3:
            std::swap(laneSine, laneCosine);
            laneSine = -laneSine;
            break;
        }
        reducedSine[lane] = laneSine;
        reducedCosine[lane] = laneCosine;
    }
    realSine = _mm256_load_pd(reducedSine);
    realCosine = _mm256_load_pd(reducedCosine);
    sineReal = _mm256_mul_pd(realSine, imaginaryCosine);
    sineImaginary = _mm256_mul_pd(realCosine, imaginarySine);
    cosineReal = _mm256_mul_pd(realCosine, imaginaryCosine);
    cosineImaginary = _mm256_xor_pd(_mm256_mul_pd(realSine, imaginarySine), signMask);
    return safeMask;
}

inline uint8_t evaluateVectorComplexSinCosFast(__m256d inputReal, __m256d inputImaginary, bool cosineOperation, __m256d& outputReal, __m256d& outputImaginary) {
    __m256d sineReal, sineImaginary, cosineReal, cosineImaginary;
    const uint8_t safeMask = evaluateVectorComplexSinCosPairFast(inputReal, inputImaginary, sineReal, sineImaginary, cosineReal, cosineImaginary);
    outputReal = cosineOperation ? cosineReal : sineReal;
    outputImaginary = cosineOperation ? cosineImaginary : sineImaginary;
    return safeMask;
}

uint8_t evaluateVectorComplexExpFast(__m256d inputReal, __m256d inputImaginary, __m256d& outputReal, __m256d& outputImaginary);

} // namespace detail
} // namespace formula

#endif
