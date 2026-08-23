#include "formula_expression_centered.h"
#include "formula_expression_vector_math.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

#ifdef _MSC_VER
#pragma float_control(precise, on, push)
#pragma fp_contract(off)
#endif

namespace formula {

namespace {

struct CenteredValue {
    Complex base{};
    Complex delta{};
    Complex endpoint{};
    ExpressionCenteredStatus status = ExpressionCenteredStatus::Success;
    bool denominatorInstability = false;
};

enum class FastEntireOperation : uint8_t { Sin,
                                           Cos,
                                           Sinh,
                                           Cosh,
                                           Exp };

bool finite(Complex value) {
    return std::isfinite(value.real()) && std::isfinite(value.imag());
}

bool zero(Complex value) {
    return value.real() == 0.0 && value.imag() == 0.0;
}

Complex nanComplex() {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    return {nan, nan};
}

CenteredValue finish(Complex base, Complex delta, Complex endpoint, ExpressionCenteredStatus status = ExpressionCenteredStatus::Success, bool denominatorInstability = false) {
    if (status == ExpressionCenteredStatus::Success && (!finite(base) || !finite(delta) || !finite(endpoint))) status = ExpressionCenteredStatus::NonFinite;
    return {base, delta, endpoint, status, denominatorInstability};
}

ExpressionCenteredStatus combine(ExpressionCenteredStatus left, ExpressionCenteredStatus right) {
    return left != ExpressionCenteredStatus::Success ? left : right;
}

bool segmentMayContainZero(Complex base, Complex endpoint) {
    if (zero(base) || zero(endpoint)) return true;
    const double minReal = std::min(base.real(), endpoint.real());
    const double maxReal = std::max(base.real(), endpoint.real());
    const double minImaginary = std::min(base.imag(), endpoint.imag());
    const double maxImaginary = std::max(base.imag(), endpoint.imag());
    const long double dx = (long double)endpoint.real() - (long double)base.real();
    const long double dy = (long double)endpoint.imag() - (long double)base.imag();
    if (dx != 0.0L || dy != 0.0L) {
        const long double x0 = base.real();
        const long double y0 = base.imag();
        if (dx != 0.0L) {
            const long double t = -x0 / dx;
            if (t >= 0.0L && t <= 1.0L && y0 + t * dy == 0.0L) return true;
        } else if (x0 == 0.0L && dy != 0.0L) {
            const long double t = -y0 / dy;
            if (t >= 0.0L && t <= 1.0L) return true;
        }
    }
    return minReal <= 0.0 && maxReal >= 0.0 && minImaginary <= 0.0 && maxImaginary >= 0.0;
}

struct CutCheck {
    bool singular = false;
    bool cut = false;
};

CutCheck checkNegativeRealCut(Complex base, Complex endpoint) {
    CutCheck result;
    if (zero(base) || zero(endpoint)) {
        result.singular = true;
        return result;
    }

    const long double x0 = base.real();
    const long double y0 = base.imag();
    const long double dx = (long double)endpoint.real() - x0;
    const long double dy = (long double)endpoint.imag() - y0;
    if (dx != 0.0L) {
        const long double t = -x0 / dx;
        if (t >= 0.0L && t <= 1.0L && y0 + t * dy == 0.0L) {
            result.singular = true;
            return result;
        }
    } else if (x0 == 0.0L && dy != 0.0L) {
        const long double t = -y0 / dy;
        if (t >= 0.0L && t <= 1.0L) {
            result.singular = true;
            return result;
        }
    }

    if ((base.imag() == 0.0 && base.real() < 0.0) || (endpoint.imag() == 0.0 && endpoint.real() < 0.0)) {
        result.cut = true;
        return result;
    }
    if (dy != 0.0L) {
        const long double t = -y0 / dy;
        if (t >= 0.0L && t <= 1.0L) {
            const long double crossingReal = x0 + t * dx;
            if (crossingReal <= 0.0L) result.cut = true;
        }
    } else if (y0 == 0.0L) {
        const long double x1 = endpoint.real();
        if (std::min(x0, x1) < 0.0L) result.cut = true;
    }
    return result;
}

bool signedZeroCutSideChanges(Complex base, Complex endpoint) {
    return base.real() < 0.0 && endpoint.real() < 0.0 && base.imag() == 0.0 && endpoint.imag() == 0.0 && std::signbit(base.imag()) != std::signbit(endpoint.imag());
}

bool tangentSegmentMayHitPole(Complex base, Complex endpoint) {
    if (!finite(base) || !finite(endpoint)) return true;
    const double minIm = std::min(base.imag(), endpoint.imag());
    const double maxIm = std::max(base.imag(), endpoint.imag());
    if (minIm > 0.0 || maxIm < 0.0) return false;

    double minX = std::min(base.real(), endpoint.real());
    double maxX = std::max(base.real(), endpoint.real());
    if (base.imag() != endpoint.imag()) {
        const double t = base.imag() / (base.imag() - endpoint.imag());
        const double xCross = base.real() + t * (endpoint.real() - base.real());
        minX = std::min(minX, xCross);
        maxX = std::max(maxX, xCross);
    }

    constexpr double Pi = 3.14159265358979323846;
    constexpr double HalfPi = 1.57079632679489661923;
    constexpr double Margin = 1e-4;

    const double kMin = std::ceil((minX - Margin - HalfPi) / Pi);
    const double kMax = std::floor((maxX + Margin - HalfPi) / Pi);
    return kMin <= kMax;
}

bool hyperbolicTangentSegmentMayHitPole(Complex base, Complex endpoint) {
    Complex swappedBase{base.imag(), base.real()};
    Complex swappedEndpoint{endpoint.imag(), endpoint.real()};
    return tangentSegmentMayHitPole(swappedBase, swappedEndpoint);
}

Complex centeredMultiply(Complex leftBase, Complex leftDelta, Complex rightBase, Complex rightDelta) {
    // Bytecode v1 fixes this order: a*db, then b*da, then da*db.
    Complex result = leftBase * rightDelta;
    result += rightBase * leftDelta;
    result += leftDelta * rightDelta;
    return result;
}

Complex centeredSquare(Complex base, Complex delta) {
    Complex result = base * delta;
    result += base * delta;
    result += delta * delta;
    return result;
}

Complex complexExpm1(Complex value) {
    if (zero(value)) return {};
    const double x = value.real();
    const double y = value.imag();
    const double sine = std::sin(y);
    const double cosine = std::cos(y);
    if (x < -0.5) {
        const double exponential = std::exp(x);
        return {exponential * cosine - 1.0, exponential * sine};
    }
    const double realExpm1 = std::expm1(x);
    const double halfSine = std::sin(y * 0.5);
    const double cosineMinusOne = -2.0 * halfSine * halfSine;
    double real = realExpm1 * cosine;
    real += cosineMinusOne;
    double imaginary = realExpm1 * sine;
    imaginary += sine;
    return {real, imaginary};
}

Complex relativeExponentialDelta(Complex baseOutput, Complex perturbedOutput, Complex exponentDelta) {
    if (zero(exponentDelta)) return {};
    if (exponentDelta.real() > 0.0) {
        Complex relative = complexExpm1(-exponentDelta);
        return -(perturbedOutput * relative);
    }
    return baseOutput * complexExpm1(exponentDelta);
}

Complex twiceSinHalf(Complex value) {
    const double magnitude = std::max(std::abs(value.real()), std::abs(value.imag()));
    if (magnitude < 1e-8) return value;
    return 2.0 * std::sin(value * 0.5);
}

Complex twiceSinhHalf(Complex value) {
    const double magnitude = std::max(std::abs(value.real()), std::abs(value.imag()));
    if (magnitude < 1e-8) return value;
    return 2.0 * std::sinh(value * 0.5);
}

Complex centeredSine(Complex base, Complex delta) {
    Complex midpoint = base + delta * 0.5;
    return std::cos(midpoint) * twiceSinHalf(delta);
}

Complex centeredCosine(Complex base, Complex delta) {
    Complex midpoint = base + delta * 0.5;
    return -std::sin(midpoint) * twiceSinHalf(delta);
}

Complex centeredHyperbolicSine(Complex base, Complex delta) {
    Complex midpoint = base + delta * 0.5;
    return std::cosh(midpoint) * twiceSinhHalf(delta);
}

Complex centeredHyperbolicCosine(Complex base, Complex delta) {
    Complex midpoint = base + delta * 0.5;
    return std::sinh(midpoint) * twiceSinhHalf(delta);
}

uint8_t vectorComplexSinhCoshFast(__m256d inputReal, __m256d inputImaginary, bool cosineOperation, __m256d& outputReal, __m256d& outputImaginary) {
    __m256d transformedReal = _mm256_xor_pd(inputImaginary, _mm256_set1_pd(-0.0));
    __m256d transformedRealOutput, transformedImaginaryOutput;
    const uint8_t safeMask = detail::evaluateVectorComplexSinCosFast(transformedReal, inputReal, cosineOperation, transformedRealOutput, transformedImaginaryOutput);
    if (cosineOperation) {
        outputReal = transformedRealOutput;
        outputImaginary = transformedImaginaryOutput;
    } else {
        outputReal = transformedImaginaryOutput;
        outputImaginary = _mm256_xor_pd(transformedRealOutput, _mm256_set1_pd(-0.0));
    }
    return safeMask;
}

uint8_t vectorComplexSinhCoshPairFast(__m256d inputReal, __m256d inputImaginary, __m256d& sineReal, __m256d& sineImaginary, __m256d& cosineReal, __m256d& cosineImaginary) {
    const __m256d transformedReal = _mm256_xor_pd(inputImaginary, _mm256_set1_pd(-0.0));
    __m256d transformedSineReal, transformedSineImaginary, transformedCosineReal, transformedCosineImaginary;
    const uint8_t safeMask = detail::evaluateVectorComplexSinCosPairFast(transformedReal, inputReal, transformedSineReal, transformedSineImaginary, transformedCosineReal, transformedCosineImaginary);
    sineReal = transformedSineImaginary;
    sineImaginary = _mm256_xor_pd(transformedSineReal, _mm256_set1_pd(-0.0));
    cosineReal = transformedCosineReal;
    cosineImaginary = transformedCosineImaginary;
    return safeMask;
}

uint8_t fusedCenteredEntireVectors(FastEntireOperation operation, __m256d inputBaseRe, __m256d inputBaseIm, __m256d outputBaseRe, __m256d outputBaseIm, bool haveCompanion, __m256d companionRe, __m256d companionIm, __m256d deltaRe, __m256d deltaIm, uint8_t validMask, __m256d& outputDeltaRe, __m256d& outputDeltaIm, __m256d* derivativeFactorRe, __m256d* derivativeFactorIm) {
    outputDeltaRe = _mm256_setzero_pd();
    outputDeltaIm = _mm256_setzero_pd();
    const bool needDerivative = derivativeFactorRe && derivativeFactorIm;
    if (needDerivative) {
        *derivativeFactorRe = _mm256_setzero_pd();
        *derivativeFactorIm = _mm256_setzero_pd();
    }
    const __m256d maximum = _mm256_set1_pd(std::numeric_limits<double>::max());
    const __m256d absoluteMask = _mm256_castsi256_pd(_mm256_set1_epi64x(0x7fffffffffffffffLL));
    __m256d usableBase = _mm256_cmp_pd(_mm256_and_pd(inputBaseRe, absoluteMask), maximum, _CMP_LE_OQ);
    usableBase = _mm256_and_pd(usableBase, _mm256_cmp_pd(_mm256_and_pd(inputBaseIm, absoluteMask), maximum, _CMP_LE_OQ));
    usableBase = _mm256_and_pd(usableBase, _mm256_cmp_pd(_mm256_and_pd(outputBaseRe, absoluteMask), maximum, _CMP_LE_OQ));
    usableBase = _mm256_and_pd(usableBase, _mm256_cmp_pd(_mm256_and_pd(outputBaseIm, absoluteMask), maximum, _CMP_LE_OQ));
    usableBase = _mm256_and_pd(usableBase, _mm256_cmp_pd(inputBaseRe, _mm256_setzero_pd(), _CMP_NEQ_OQ));
    usableBase = _mm256_and_pd(usableBase, _mm256_cmp_pd(inputBaseIm, _mm256_setzero_pd(), _CMP_NEQ_OQ));
    usableBase = _mm256_and_pd(usableBase, _mm256_cmp_pd(outputBaseRe, _mm256_setzero_pd(), _CMP_NEQ_OQ));
    usableBase = _mm256_and_pd(usableBase, _mm256_cmp_pd(outputBaseIm, _mm256_setzero_pd(), _CMP_NEQ_OQ));
    validMask &= static_cast<uint8_t>(_mm256_movemask_pd(usableBase));
    if (validMask == 0) return 0;

    const __m256d half = _mm256_set1_pd(0.5);
    const __m256d halfDeltaRe = _mm256_mul_pd(deltaRe, half);
    const __m256d halfDeltaIm = _mm256_mul_pd(deltaIm, half);
    const __m256d midpointRe = _mm256_add_pd(inputBaseRe, halfDeltaRe);
    const __m256d midpointIm = _mm256_add_pd(inputBaseIm, halfDeltaIm);
    __m256d factorRe, factorIm, halfFunctionRe, halfFunctionIm;
    uint8_t factorMask = 0;
    uint8_t halfFunctionMask = 0;
    double scale = 2.0;
    auto multiply = [](__m256d leftRe, __m256d leftIm, __m256d rightRe, __m256d rightIm, __m256d& productRe, __m256d& productIm) {
        productRe = _mm256_sub_pd(_mm256_mul_pd(leftRe, rightRe), _mm256_mul_pd(leftIm, rightIm));
        productIm = _mm256_add_pd(_mm256_mul_pd(leftRe, rightIm), _mm256_mul_pd(leftIm, rightRe));
    };
    if (haveCompanion && operation != FastEntireOperation::Exp) {
        __m256d halfSineRe, halfSineIm, halfCosineRe, halfCosineIm;
        const __m256d signMask = _mm256_set1_pd(-0.0);
        const __m256d absoluteDeltaRe = _mm256_andnot_pd(signMask, deltaRe);
        const __m256d absoluteDeltaIm = _mm256_andnot_pd(signMask, deltaIm);
        __m256d tiny = _mm256_cmp_pd(absoluteDeltaRe, _mm256_set1_pd(1e-8), _CMP_LT_OQ);
        tiny = _mm256_and_pd(tiny, _mm256_cmp_pd(absoluteDeltaIm, _mm256_set1_pd(1e-8), _CMP_LT_OQ));
        tiny = _mm256_and_pd(tiny, _mm256_cmp_pd(deltaRe, _mm256_setzero_pd(), _CMP_NEQ_OQ));
        tiny = _mm256_and_pd(tiny, _mm256_cmp_pd(deltaIm, _mm256_setzero_pd(), _CMP_NEQ_OQ));
        const uint8_t tinyMask = validMask & static_cast<uint8_t>(_mm256_movemask_pd(tiny));
        if (tinyMask == validMask) {
            __m256d squareRe, squareIm;
            squareRe = _mm256_sub_pd(_mm256_mul_pd(halfDeltaRe, halfDeltaRe), _mm256_mul_pd(halfDeltaIm, halfDeltaIm));
            squareIm = _mm256_add_pd(_mm256_mul_pd(halfDeltaRe, halfDeltaIm), _mm256_mul_pd(halfDeltaIm, halfDeltaRe));
            const __m256d halfSquareRe = _mm256_mul_pd(_mm256_set1_pd(0.5), squareRe);
            const __m256d halfSquareIm = _mm256_mul_pd(_mm256_set1_pd(0.5), squareIm);
            halfSineRe = halfDeltaRe;
            halfSineIm = halfDeltaIm;
            if (operation == FastEntireOperation::Sin || operation == FastEntireOperation::Cos) {
                halfCosineRe = _mm256_sub_pd(_mm256_set1_pd(1.0), halfSquareRe);
                halfCosineIm = _mm256_xor_pd(halfSquareIm, signMask);
            } else {
                halfCosineRe = _mm256_add_pd(_mm256_set1_pd(1.0), halfSquareRe);
                halfCosineIm = halfSquareIm;
            }
            halfFunctionMask = tinyMask;
        } else {
            __m256d vectorSineRe, vectorSineIm, vectorCosineRe, vectorCosineIm;
            const uint8_t vectorMask = operation == FastEntireOperation::Sin || operation == FastEntireOperation::Cos ? detail::evaluateVectorComplexSinCosPairFast(halfDeltaRe, halfDeltaIm, vectorSineRe, vectorSineIm, vectorCosineRe, vectorCosineIm) : vectorComplexSinhCoshPairFast(halfDeltaRe, halfDeltaIm, vectorSineRe, vectorSineIm, vectorCosineRe, vectorCosineIm);
            if (tinyMask == 0) {
                halfSineRe = vectorSineRe;
                halfSineIm = vectorSineIm;
                halfCosineRe = vectorCosineRe;
                halfCosineIm = vectorCosineIm;
            } else {
                __m256d squareRe = _mm256_sub_pd(_mm256_mul_pd(halfDeltaRe, halfDeltaRe), _mm256_mul_pd(halfDeltaIm, halfDeltaIm));
                __m256d squareIm = _mm256_add_pd(_mm256_mul_pd(halfDeltaRe, halfDeltaIm), _mm256_mul_pd(halfDeltaIm, halfDeltaRe));
                const __m256d halfSquareRe = _mm256_mul_pd(_mm256_set1_pd(0.5), squareRe);
                const __m256d halfSquareIm = _mm256_mul_pd(_mm256_set1_pd(0.5), squareIm);
                const __m256d tinySineRe = halfDeltaRe;
                const __m256d tinySineIm = halfDeltaIm;
                const __m256d tinyCosineRe = operation == FastEntireOperation::Sin || operation == FastEntireOperation::Cos ? _mm256_sub_pd(_mm256_set1_pd(1.0), halfSquareRe) : _mm256_add_pd(_mm256_set1_pd(1.0), halfSquareRe);
                const __m256d tinyCosineIm = operation == FastEntireOperation::Sin || operation == FastEntireOperation::Cos ? _mm256_xor_pd(halfSquareIm, signMask) : halfSquareIm;
                const __m256d tinyVectorMask = _mm256_castsi256_pd(_mm256_set_epi64x((tinyMask & 8) ? -1LL : 0LL, (tinyMask & 4) ? -1LL : 0LL, (tinyMask & 2) ? -1LL : 0LL, (tinyMask & 1) ? -1LL : 0LL));
                halfSineRe = _mm256_blendv_pd(vectorSineRe, tinySineRe, tinyVectorMask);
                halfSineIm = _mm256_blendv_pd(vectorSineIm, tinySineIm, tinyVectorMask);
                halfCosineRe = _mm256_blendv_pd(vectorCosineRe, tinyCosineRe, tinyVectorMask);
                halfCosineIm = _mm256_blendv_pd(vectorCosineIm, tinyCosineIm, tinyVectorMask);
            }
            halfFunctionMask = tinyMask | vectorMask;
        }
        const __m256d firstRe = companionRe;
        const __m256d firstIm = companionIm;
        const __m256d secondRe = outputBaseRe;
        const __m256d secondIm = outputBaseIm;
        const __m256d firstProductRe = _mm256_sub_pd(_mm256_mul_pd(firstRe, halfCosineRe), _mm256_mul_pd(firstIm, halfCosineIm));
        const __m256d firstProductIm = _mm256_add_pd(_mm256_mul_pd(firstRe, halfCosineIm), _mm256_mul_pd(firstIm, halfCosineRe));
        const __m256d secondProductRe = _mm256_sub_pd(_mm256_mul_pd(secondRe, halfSineRe), _mm256_mul_pd(secondIm, halfSineIm));
        const __m256d secondProductIm = _mm256_add_pd(_mm256_mul_pd(secondRe, halfSineIm), _mm256_mul_pd(secondIm, halfSineRe));
        if (operation == FastEntireOperation::Sin) {
            factorRe = _mm256_sub_pd(firstProductRe, secondProductRe);
            factorIm = _mm256_sub_pd(firstProductIm, secondProductIm);
        } else {
            factorRe = _mm256_add_pd(firstProductRe, secondProductRe);
            factorIm = _mm256_add_pd(firstProductIm, secondProductIm);
        }
        halfFunctionRe = halfSineRe;
        halfFunctionIm = halfSineIm;
        factorMask = halfFunctionMask;
        if (operation == FastEntireOperation::Cos) scale = -2.0;
        if (needDerivative) {
            __m256d firstSineRe, firstSineIm, secondCosineRe, secondCosineIm;
            multiply(firstRe, firstIm, halfSineRe, halfSineIm, firstSineRe, firstSineIm);
            multiply(secondRe, secondIm, halfCosineRe, halfCosineIm, secondCosineRe, secondCosineIm);
            __m256d midpointValueRe;
            __m256d midpointValueIm;
            if (operation == FastEntireOperation::Cos) {
                midpointValueRe = _mm256_sub_pd(secondCosineRe, firstSineRe);
                midpointValueIm = _mm256_sub_pd(secondCosineIm, firstSineIm);
            } else {
                midpointValueRe = _mm256_add_pd(secondCosineRe, firstSineRe);
                midpointValueIm = _mm256_add_pd(secondCosineIm, firstSineIm);
            }
            __m256d midpointSineRe, midpointSineIm;
            multiply(midpointValueRe, midpointValueIm, halfSineRe, halfSineIm, midpointSineRe, midpointSineIm);
            midpointSineRe = _mm256_add_pd(midpointSineRe, midpointSineRe);
            midpointSineIm = _mm256_add_pd(midpointSineIm, midpointSineIm);
            if (operation == FastEntireOperation::Sin) {
                *derivativeFactorRe = _mm256_sub_pd(firstRe, midpointSineRe);
                *derivativeFactorIm = _mm256_sub_pd(firstIm, midpointSineIm);
            } else if (operation == FastEntireOperation::Cos) {
                *derivativeFactorRe = _mm256_xor_pd(_mm256_add_pd(firstRe, midpointSineRe), signMask);
                *derivativeFactorIm = _mm256_xor_pd(_mm256_add_pd(firstIm, midpointSineIm), signMask);
            } else {
                *derivativeFactorRe = _mm256_add_pd(firstRe, midpointSineRe);
                *derivativeFactorIm = _mm256_add_pd(firstIm, midpointSineIm);
            }
        }
    } else {
        switch (operation) {
        case FastEntireOperation::Sin:
            factorMask = detail::evaluateVectorComplexSinCosFast(midpointRe, midpointIm, true, factorRe, factorIm);
            halfFunctionMask = detail::evaluateVectorComplexSinCosFast(halfDeltaRe, halfDeltaIm, false, halfFunctionRe, halfFunctionIm);
            break;
        case FastEntireOperation::Cos:
            factorMask = detail::evaluateVectorComplexSinCosFast(midpointRe, midpointIm, false, factorRe, factorIm);
            halfFunctionMask = detail::evaluateVectorComplexSinCosFast(halfDeltaRe, halfDeltaIm, false, halfFunctionRe, halfFunctionIm);
            scale = -2.0;
            break;
        case FastEntireOperation::Sinh:
            factorMask = vectorComplexSinhCoshFast(midpointRe, midpointIm, true, factorRe, factorIm);
            halfFunctionMask = vectorComplexSinhCoshFast(halfDeltaRe, halfDeltaIm, false, halfFunctionRe, halfFunctionIm);
            break;
        case FastEntireOperation::Cosh:
            factorMask = vectorComplexSinhCoshFast(midpointRe, midpointIm, false, factorRe, factorIm);
            halfFunctionMask = vectorComplexSinhCoshFast(halfDeltaRe, halfDeltaIm, false, halfFunctionRe, halfFunctionIm);
            break;
        case FastEntireOperation::Exp: {
            if (haveCompanion) {
                __m256d halfExpRe, halfExpIm;
                factorMask = detail::evaluateVectorComplexExpFast(halfDeltaRe, halfDeltaIm, halfExpRe, halfExpIm);
                factorRe = _mm256_sub_pd(_mm256_mul_pd(outputBaseRe, halfExpRe), _mm256_mul_pd(outputBaseIm, halfExpIm));
                factorIm = _mm256_add_pd(_mm256_mul_pd(outputBaseRe, halfExpIm), _mm256_mul_pd(outputBaseIm, halfExpRe));
            } else {
                factorMask = detail::evaluateVectorComplexExpFast(midpointRe, midpointIm, factorRe, factorIm);
            }
            halfFunctionMask = vectorComplexSinhCoshFast(halfDeltaRe, halfDeltaIm, false, halfFunctionRe, halfFunctionIm);
            break;
        }
        }
    }
    uint8_t safeMask = validMask & factorMask & halfFunctionMask;
    outputDeltaRe = _mm256_mul_pd(_mm256_set1_pd(scale), _mm256_sub_pd(_mm256_mul_pd(factorRe, halfFunctionRe), _mm256_mul_pd(factorIm, halfFunctionIm)));
    outputDeltaIm = _mm256_mul_pd(_mm256_set1_pd(scale), _mm256_add_pd(_mm256_mul_pd(factorRe, halfFunctionIm), _mm256_mul_pd(factorIm, halfFunctionRe)));
    const __m256d endpointRe = _mm256_add_pd(outputBaseRe, outputDeltaRe);
    const __m256d endpointIm = _mm256_add_pd(outputBaseIm, outputDeltaIm);
    if (needDerivative && operation == FastEntireOperation::Exp) {
        *derivativeFactorRe = endpointRe;
        *derivativeFactorIm = endpointIm;
    }
    __m256d usable = _mm256_cmp_pd(_mm256_and_pd(outputDeltaRe, absoluteMask), maximum, _CMP_LE_OQ);
    usable = _mm256_and_pd(usable, _mm256_cmp_pd(_mm256_and_pd(outputDeltaIm, absoluteMask), maximum, _CMP_LE_OQ));
    usable = _mm256_and_pd(usable, _mm256_cmp_pd(outputDeltaRe, _mm256_setzero_pd(), _CMP_NEQ_OQ));
    usable = _mm256_and_pd(usable, _mm256_cmp_pd(outputDeltaIm, _mm256_setzero_pd(), _CMP_NEQ_OQ));
    usable = _mm256_and_pd(usable, _mm256_cmp_pd(_mm256_and_pd(endpointRe, absoluteMask), maximum, _CMP_LE_OQ));
    usable = _mm256_and_pd(usable, _mm256_cmp_pd(_mm256_and_pd(endpointIm, absoluteMask), maximum, _CMP_LE_OQ));
    usable = _mm256_and_pd(usable, _mm256_cmp_pd(endpointRe, _mm256_setzero_pd(), _CMP_NEQ_OQ));
    usable = _mm256_and_pd(usable, _mm256_cmp_pd(endpointIm, _mm256_setzero_pd(), _CMP_NEQ_OQ));
    if (needDerivative) {
        usable = _mm256_and_pd(usable, _mm256_cmp_pd(_mm256_and_pd(*derivativeFactorRe, absoluteMask), maximum, _CMP_LE_OQ));
        usable = _mm256_and_pd(usable, _mm256_cmp_pd(_mm256_and_pd(*derivativeFactorIm, absoluteMask), maximum, _CMP_LE_OQ));
        if (!haveCompanion) usable = _mm256_setzero_pd();
    }
    safeMask &= static_cast<uint8_t>(_mm256_movemask_pd(usable));
    return safeMask;
}

uint8_t fusedCenteredEntireVectors(FastEntireOperation operation, Complex inputBase, Complex outputBase, const Complex* companion, __m256d deltaRe, __m256d deltaIm, uint8_t validMask, __m256d& outputDeltaRe, __m256d& outputDeltaIm, __m256d* derivativeFactorRe, __m256d* derivativeFactorIm) {
    const Complex companionValue = companion ? *companion : Complex{};
    return fusedCenteredEntireVectors(operation, _mm256_set1_pd(inputBase.real()), _mm256_set1_pd(inputBase.imag()), _mm256_set1_pd(outputBase.real()), _mm256_set1_pd(outputBase.imag()), companion != nullptr, _mm256_set1_pd(companionValue.real()), _mm256_set1_pd(companionValue.imag()), deltaRe, deltaIm, validMask, outputDeltaRe, outputDeltaIm, derivativeFactorRe, derivativeFactorIm);
}

ExpressionCenteredStatus centeredDivisionDelta(Complex numeratorDelta, Complex denominatorBase, Complex denominatorDelta, Complex denominatorEndpoint, Complex quotientBase, Complex& output, bool guardLinearSegment) {
    if (zero(denominatorBase) || zero(denominatorEndpoint) || (guardLinearSegment && segmentMayContainZero(denominatorBase, denominatorEndpoint))) return ExpressionCenteredStatus::Singular;
    Complex numerator = numeratorDelta - quotientBase * denominatorDelta;
    output = numerator / denominatorEndpoint;
    return finite(output) ? ExpressionCenteredStatus::Success : ExpressionCenteredStatus::NonFinite;
}

double centeredRealProduct(double leftBase, double leftDelta, double rightBase, double rightDelta) {
    double result = leftBase * rightDelta;
    result += rightBase * leftDelta;
    result += leftDelta * rightDelta;
    return result;
}

ExpressionCenteredStatus centeredTangent(Complex base, Complex delta, Complex endpoint, Complex baseOutput, Complex endpointOutput, Complex& output) {
    if (zero(delta)) {
        output = {};
        return ExpressionCenteredStatus::Success;
    }
    if (tangentSegmentMayHitPole(base, endpoint)) return ExpressionCenteredStatus::Singular;
    if (!finite(baseOutput) || !finite(endpointOutput)) return ExpressionCenteredStatus::NonFinite;

    const bool sameHighHalfPlane = std::abs(base.imag()) > 16.0 && std::abs(endpoint.imag()) > 16.0 && std::signbit(base.imag()) == std::signbit(endpoint.imag());
    if (!sameHighHalfPlane) {
        Complex sineBase = std::sin(base);
        Complex cosineBase = std::cos(base);
        Complex cosineEndpoint = std::cos(endpoint);
        Complex sineDelta = centeredSine(base, delta);
        Complex cosineDelta = centeredCosine(base, delta);
        if (!finite(sineBase) || !finite(cosineBase) || !finite(cosineEndpoint) || !finite(sineDelta) || !finite(cosineDelta)) return ExpressionCenteredStatus::NonFinite;
        return centeredDivisionDelta(sineDelta, cosineBase, cosineDelta, cosineEndpoint, baseOutput, output, true);
    }

    const double sign = std::signbit(base.imag()) ? -1.0 : 1.0;
    const double q = std::exp(-2.0 * std::abs(base.imag()));
    const double q1 = std::exp(-2.0 * std::abs(endpoint.imag()));
    const double qDelta = relativeExponentialDelta({q, 0.0}, {q1, 0.0}, {-2.0 * sign * delta.imag(), 0.0}).real();
    const double angle = 2.0 * base.real();
    const double angleDelta = 2.0 * delta.real();
    const double sine = std::sin(angle);
    const double cosine = std::cos(angle);
    const double sineDelta = centeredSine({angle, 0.0}, {angleDelta, 0.0}).real();
    const double cosineDelta = centeredCosine({angle, 0.0}, {angleDelta, 0.0}).real();
    double qSquaredDelta = q * qDelta;
    qSquaredDelta += q * qDelta;
    qSquaredDelta += qDelta * qDelta;
    const double numeratorRealDelta = 2.0 * centeredRealProduct(q, qDelta, sine, sineDelta);
    const double numeratorImaginaryDelta = -sign * qSquaredDelta;
    double denominatorDelta = qSquaredDelta;
    denominatorDelta += 2.0 * centeredRealProduct(q, qDelta, cosine, cosineDelta);
    double denominator = 1.0 + q * q;
    denominator += 2.0 * q * cosine;
    double endpointDenominator = 1.0 + q1 * q1;
    endpointDenominator += 2.0 * q1 * std::cos(2.0 * endpoint.real());
    return centeredDivisionDelta({numeratorRealDelta, numeratorImaginaryDelta}, {denominator, 0.0}, {denominatorDelta, 0.0}, {endpointDenominator, 0.0}, baseOutput, output, true);
}

ExpressionCenteredStatus centeredHyperbolicTangent(Complex base, Complex delta, Complex endpoint, Complex baseOutput, Complex endpointOutput, Complex& output) {
    if (zero(delta)) {
        output = {};
        return ExpressionCenteredStatus::Success;
    }
    if (hyperbolicTangentSegmentMayHitPole(base, endpoint)) return ExpressionCenteredStatus::Singular;
    if (!finite(baseOutput) || !finite(endpointOutput)) return ExpressionCenteredStatus::NonFinite;

    const bool sameHighHalfPlane = std::abs(base.real()) > 16.0 && std::abs(endpoint.real()) > 16.0 && std::signbit(base.real()) == std::signbit(endpoint.real());
    if (!sameHighHalfPlane) {
        Complex sineBase = std::sinh(base);
        Complex cosineBase = std::cosh(base);
        Complex cosineEndpoint = std::cosh(endpoint);
        Complex sineDelta = centeredHyperbolicSine(base, delta);
        Complex cosineDelta = centeredHyperbolicCosine(base, delta);
        if (!finite(sineBase) || !finite(cosineBase) || !finite(cosineEndpoint) || !finite(sineDelta) || !finite(cosineDelta)) return ExpressionCenteredStatus::NonFinite;
        return centeredDivisionDelta(sineDelta, cosineBase, cosineDelta, cosineEndpoint, baseOutput, output, true);
    }

    const double sign = std::signbit(base.real()) ? -1.0 : 1.0;
    const double q = std::exp(-2.0 * std::abs(base.real()));
    const double q1 = std::exp(-2.0 * std::abs(endpoint.real()));
    const double qDelta = relativeExponentialDelta({q, 0.0}, {q1, 0.0}, {-2.0 * sign * delta.real(), 0.0}).real();
    const double angle = 2.0 * base.imag();
    const double angleDelta = 2.0 * delta.imag();
    const double sine = std::sin(angle);
    const double cosine = std::cos(angle);
    const double sineDelta = centeredSine({angle, 0.0}, {angleDelta, 0.0}).real();
    const double cosineDelta = centeredCosine({angle, 0.0}, {angleDelta, 0.0}).real();
    double qSquaredDelta = q * qDelta;
    qSquaredDelta += q * qDelta;
    qSquaredDelta += qDelta * qDelta;
    const double numeratorRealDelta = -sign * qSquaredDelta;
    const double numeratorImaginaryDelta = 2.0 * centeredRealProduct(q, qDelta, sine, sineDelta);
    double denominatorDelta = qSquaredDelta;
    denominatorDelta += 2.0 * centeredRealProduct(q, qDelta, cosine, cosineDelta);
    double denominator = 1.0 + q * q;
    denominator += 2.0 * q * cosine;
    double endpointDenominator = 1.0 + q1 * q1;
    endpointDenominator += 2.0 * q1 * std::cos(2.0 * endpoint.imag());
    return centeredDivisionDelta({numeratorRealDelta, numeratorImaginaryDelta}, {denominator, 0.0}, {denominatorDelta, 0.0}, {endpointDenominator, 0.0}, baseOutput, output, true);
}

ExpressionCenteredStatus centeredLogarithm(Complex base, Complex delta, Complex endpoint, Complex& output) {
    if (zero(base) || zero(endpoint)) return ExpressionCenteredStatus::Singular;
    if (signedZeroCutSideChanges(base, endpoint)) return ExpressionCenteredStatus::BranchUncertain;
    if (zero(delta)) {
        output = {};
        return ExpressionCenteredStatus::Success;
    }
    CutCheck cut = checkNegativeRealCut(base, endpoint);
    if (cut.singular) return ExpressionCenteredStatus::Singular;
    if (cut.cut) return ExpressionCenteredStatus::BranchUncertain;
    Complex ratio = delta / base;
    if (!finite(ratio)) return ExpressionCenteredStatus::NonFinite;
    const double x = ratio.real();
    const double y = ratio.imag();
    const double onePlusX = 1.0 + x;
    const double magnitude = std::hypot(onePlusX, y);
    if (magnitude == 0.0) return ExpressionCenteredStatus::Singular;

    double real;
    if (std::max(std::abs(x), std::abs(y)) <= 0.5) {
        double squared = x * x;
        squared += y * y;
        double logArgument = x + x;
        logArgument += squared;
        if (logArgument <= -1.0) {
            real = std::log(magnitude);
        } else {
            real = 0.5 * std::log1p(logArgument);
        }
    } else {
        real = std::log(magnitude);
    }
    output = {real, std::atan2(y, onePlusX)};
    return finite(output) ? ExpressionCenteredStatus::Success : ExpressionCenteredStatus::NonFinite;
}

Complex centeredNormDelta(Complex base, Complex delta) {
    Complex cross = std::conj(base) * delta;
    double real = cross.real() + cross.real();
    real += std::norm(delta);
    return {real, 0.0};
}

ExpressionCenteredStatus centeredAbsolute(Complex base, Complex delta, Complex endpoint, Complex& output) {
    const double baseRadius = std::abs(base);
    const double perturbedRadius = std::abs(endpoint);
    const double scale = std::max(baseRadius, perturbedRadius);
    if (!std::isfinite(scale)) return ExpressionCenteredStatus::NonFinite;
    if (scale == 0.0) {
        output = {};
        return ExpressionCenteredStatus::Success;
    }
    Complex scaledBase = base / scale;
    Complex scaledDelta = delta / scale;
    double numerator = centeredNormDelta(scaledBase, scaledDelta).real();
    double denominator = baseRadius / scale + perturbedRadius / scale;
    if (denominator == 0.0) return ExpressionCenteredStatus::Singular;
    double residual = numerator / denominator;
    residual *= scale;
    output = {residual, 0.0};
    return finite(output) ? ExpressionCenteredStatus::Success : ExpressionCenteredStatus::NonFinite;
}

ExpressionCenteredStatus centeredSquareRoot(Complex base, Complex delta, Complex endpoint, Complex baseRoot, Complex endpointRoot, Complex& output) {
    if (zero(base) && zero(endpoint)) return ExpressionCenteredStatus::Singular;
    if (signedZeroCutSideChanges(base, endpoint)) return ExpressionCenteredStatus::BranchUncertain;
    if (zero(delta)) {
        output = {};
        return ExpressionCenteredStatus::Success;
    }
    CutCheck cut = checkNegativeRealCut(base, endpoint);
    if (cut.cut) return ExpressionCenteredStatus::BranchUncertain;
    Complex denominator = endpointRoot + baseRoot;
    if (zero(denominator) || cut.singular) return ExpressionCenteredStatus::Singular;
    output = delta / denominator;
    return finite(output) ? ExpressionCenteredStatus::Success : ExpressionCenteredStatus::NonFinite;
}

ExpressionCenteredStatus centeredPower(Complex base, Complex baseDelta, Complex baseEndpoint, Complex exponent, Complex exponentDelta, Complex exponentEndpoint, Complex baseOutput, Complex endpointOutput, Complex& output) {
    if (zero(base)) {
        if (zero(baseEndpoint)) {
            output = {};
            return ExpressionCenteredStatus::Success;
        }
        if (exponent.real() > 0.0 && exponentEndpoint.real() > 0.0) {
            output = std::pow(baseEndpoint, exponentEndpoint);
            return finite(output) ? ExpressionCenteredStatus::Success : ExpressionCenteredStatus::NonFinite;
        }
        return ExpressionCenteredStatus::Singular;
    }
    if (zero(baseEndpoint)) return ExpressionCenteredStatus::Singular;
    if (signedZeroCutSideChanges(base, baseEndpoint)) return ExpressionCenteredStatus::BranchUncertain;
    if (zero(baseDelta) && zero(exponentDelta)) {
        output = {};
        return ExpressionCenteredStatus::Success;
    }
    CutCheck cut = checkNegativeRealCut(base, baseEndpoint);
    if (cut.singular) return ExpressionCenteredStatus::Singular;
    if (cut.cut) return ExpressionCenteredStatus::BranchUncertain;

    Complex logBase = std::log(base);
    Complex logDelta;
    ExpressionCenteredStatus status = centeredLogarithm(base, baseDelta, baseEndpoint, logDelta);
    if (status != ExpressionCenteredStatus::Success) return status;
    Complex productDelta = centeredMultiply(exponent, exponentDelta, logBase, logDelta);
    if (!finite(exponentEndpoint) || !finite(baseOutput) || !finite(endpointOutput)) return ExpressionCenteredStatus::NonFinite;
    output = relativeExponentialDelta(baseOutput, endpointOutput, productDelta);
    return finite(output) ? ExpressionCenteredStatus::Success : ExpressionCenteredStatus::NonFinite;
}

ExpressionCenteredStatus centeredPolar(const CenteredValue& radius, const CenteredValue& angle, Complex baseOutput, Complex endpointOutput, Complex& output) {
    if (radius.base.imag() != 0.0 || radius.delta.imag() != 0.0 || radius.endpoint.imag() != 0.0 || angle.base.imag() != 0.0 || angle.delta.imag() != 0.0 || angle.endpoint.imag() != 0.0) return ExpressionCenteredStatus::Undefined;
    const double baseRadius = radius.base.real();
    const double radiusDelta = radius.delta.real();
    const double baseAngle = angle.base.real();
    const double angleDelta = angle.delta.real();
    const double perturbedRadius = radius.endpoint.real();
    const double perturbedAngle = angle.endpoint.real();
    if (!std::isfinite(baseRadius) || !std::isfinite(radiusDelta) || !std::isfinite(baseAngle) || !std::isfinite(angleDelta) || !std::isfinite(perturbedRadius) || !std::isfinite(perturbedAngle)) return ExpressionCenteredStatus::NonFinite;
    if (baseRadius < 0.0 || perturbedRadius < 0.0) return ExpressionCenteredStatus::Undefined;
    if (!finite(baseOutput) || !finite(endpointOutput)) return ExpressionCenteredStatus::NonFinite;
    Complex angleResidual = complexExpm1({0.0, angleDelta});
    Complex result = baseOutput * angleResidual;
    result += radiusDelta * std::polar(1.0, perturbedAngle);
    output = result;
    return finite(output) ? ExpressionCenteredStatus::Success : ExpressionCenteredStatus::NonFinite;
}

} // namespace

const char* expressionCenteredStatusName(ExpressionCenteredStatus status) {
    switch (status) {
    case ExpressionCenteredStatus::Success: return "success";
    case ExpressionCenteredStatus::BranchUncertain: return "branch-uncertain";
    case ExpressionCenteredStatus::Singular: return "singular";
    case ExpressionCenteredStatus::Undefined: return "undefined";
    case ExpressionCenteredStatus::Unsupported: return "unsupported";
    case ExpressionCenteredStatus::NonFinite: return "nonfinite";
    case ExpressionCenteredStatus::InvalidProgram: return "invalid-program";
    }
    return "invalid-program";
}

ExpressionCenteredResult ExpressionCenteredEvaluator::evaluate(const ExpressionProgram& program, const ExpressionContext& reference, const ExpressionDeltaContext& inputDelta) {
    ExpressionCenteredResult invalid;
    invalid.base = nanComplex();
    invalid.delta = nanComplex();
    if (!program._valid) return invalid;

    std::array<CenteredValue, ExpressionProgram::MAX_STACK> stack;
    size_t top = 0;
    auto pushUnchanged = [&](Complex value) { stack[top++] = finish(value, {}, value); };
    auto pushPerturbed = [&](Complex base, Complex delta) {
        Complex endpoint = base;
        endpoint += delta;
        stack[top++] = finish(base, delta, endpoint);
    };
    auto unary = [&](const ExpressionProgram::Instruction& instruction) {
        CenteredValue input = stack[top - 1];
        Complex base = ExpressionProgram::evaluateUnary(instruction.op, input.base);
        Complex endpoint = ExpressionProgram::evaluateUnary(instruction.op, input.endpoint);
        if (input.status != ExpressionCenteredStatus::Success) {
            stack[top - 1] = finish(base, nanComplex(), endpoint, input.status, input.denominatorInstability);
            return;
        }

        Complex delta;
        ExpressionCenteredStatus status = ExpressionCenteredStatus::Success;
        bool denominatorInstability = input.denominatorInstability;
        switch (instruction.op) {
        case ExpressionProgram::Op::Negate: delta = -input.delta; break;
        case ExpressionProgram::Op::Square: delta = centeredSquare(input.base, input.delta); break;
        case ExpressionProgram::Op::Sin: delta = centeredSine(input.base, input.delta); break;
        case ExpressionProgram::Op::Cos: delta = centeredCosine(input.base, input.delta); break;
        case ExpressionProgram::Op::Tan:
            status = centeredTangent(input.base, input.delta, input.endpoint, base, endpoint, delta);
            denominatorInstability = denominatorInstability || status != ExpressionCenteredStatus::Success;
            break;
        case ExpressionProgram::Op::Sinh: delta = centeredHyperbolicSine(input.base, input.delta); break;
        case ExpressionProgram::Op::Cosh: delta = centeredHyperbolicCosine(input.base, input.delta); break;
        case ExpressionProgram::Op::Tanh:
            status = centeredHyperbolicTangent(input.base, input.delta, input.endpoint, base, endpoint, delta);
            denominatorInstability = denominatorInstability || status != ExpressionCenteredStatus::Success;
            break;
        case ExpressionProgram::Op::Exp: {
            if (!finite(base) || !finite(endpoint)) {
                status = ExpressionCenteredStatus::NonFinite;
            } else {
                delta = relativeExponentialDelta(base, endpoint, input.delta);
            }
            break;
        }
        case ExpressionProgram::Op::Log: status = centeredLogarithm(input.base, input.delta, input.endpoint, delta); break;
        case ExpressionProgram::Op::Log10:
            status = centeredLogarithm(input.base, input.delta, input.endpoint, delta);
            if (status == ExpressionCenteredStatus::Success) delta /= std::log(10.0);
            break;
        case ExpressionProgram::Op::Sqrt: status = centeredSquareRoot(input.base, input.delta, input.endpoint, base, endpoint, delta); break;
        case ExpressionProgram::Op::Abs: status = centeredAbsolute(input.base, input.delta, input.endpoint, delta); break;
        case ExpressionProgram::Op::Norm: delta = centeredNormDelta(input.base, input.delta); break;
        case ExpressionProgram::Op::Arg:
            status = centeredLogarithm(input.base, input.delta, input.endpoint, delta);
            if (status == ExpressionCenteredStatus::Success) delta = {delta.imag(), 0.0};
            break;
        case ExpressionProgram::Op::Conjugate: delta = std::conj(input.delta); break;
        case ExpressionProgram::Op::Real: delta = {input.delta.real(), 0.0}; break;
        case ExpressionProgram::Op::Imaginary: delta = {input.delta.imag(), 0.0}; break;
        default:
            status = ExpressionCenteredStatus::Unsupported;
            delta = nanComplex();
            break;
        }
        stack[top - 1] = finish(base, delta, endpoint, status, denominatorInstability);
    };
    auto binary = [&](const ExpressionProgram::Instruction& instruction) {
        CenteredValue right = stack[--top];
        CenteredValue left = stack[top - 1];
        Complex base = ExpressionProgram::evaluateBinary(instruction.op, left.base, right.base);
        Complex endpoint = ExpressionProgram::evaluateBinary(instruction.op, left.endpoint, right.endpoint);
        ExpressionCenteredStatus status = combine(left.status, right.status);
        bool denominatorInstability = left.denominatorInstability || right.denominatorInstability;
        if (status != ExpressionCenteredStatus::Success) {
            stack[top - 1] = finish(base, nanComplex(), endpoint, status, denominatorInstability);
            return;
        }

        Complex delta;
        switch (instruction.op) {
        case ExpressionProgram::Op::Add: delta = left.delta + right.delta; break;
        case ExpressionProgram::Op::Subtract: delta = left.delta - right.delta; break;
        case ExpressionProgram::Op::Multiply: delta = centeredMultiply(left.base, left.delta, right.base, right.delta); break;
        case ExpressionProgram::Op::Divide: {
            Complex quotientDelta;
            status = centeredDivisionDelta(left.delta, right.base, right.delta, right.endpoint, base, quotientDelta, true);
            denominatorInstability = denominatorInstability || status != ExpressionCenteredStatus::Success;
            delta = quotientDelta;
            break;
        }
        case ExpressionProgram::Op::Power: status = centeredPower(left.base, left.delta, left.endpoint, right.base, right.delta, right.endpoint, base, endpoint, delta); break;
        case ExpressionProgram::Op::MakeComplex: delta = {left.delta.real(), right.delta.real()}; break;
        case ExpressionProgram::Op::Polar: status = centeredPolar(left, right, base, endpoint, delta); break;
        default:
            status = ExpressionCenteredStatus::Unsupported;
            delta = nanComplex();
            break;
        }
        stack[top - 1] = finish(base, delta, endpoint, status, denominatorInstability);
    };

    for (const ExpressionProgram::Instruction& instruction : program._code) {
        switch (instruction.op) {
        case ExpressionProgram::Op::Constant: pushUnchanged(instruction.value); break;
        case ExpressionProgram::Op::Z: pushPerturbed(reference.z, inputDelta.z); break;
        case ExpressionProgram::Op::C: pushPerturbed(reference.c, inputDelta.c); break;
        case ExpressionProgram::Op::Z0: pushPerturbed(reference.z0, inputDelta.z0); break;
        case ExpressionProgram::Op::Iteration: pushUnchanged({(double)reference.iteration, 0.0}); break;
        case ExpressionProgram::Op::Parameter: pushPerturbed(reference.parameters[instruction.argument], inputDelta.parameters[instruction.argument]); break;
        case ExpressionProgram::Op::OrbitInvariant: stack[top++] = finish(nanComplex(), nanComplex(), nanComplex(), ExpressionCenteredStatus::Unsupported); break;
        case ExpressionProgram::Op::Negate:
        case ExpressionProgram::Op::Square:
        case ExpressionProgram::Op::Sin:
        case ExpressionProgram::Op::Cos:
        case ExpressionProgram::Op::Tan:
        case ExpressionProgram::Op::Sinh:
        case ExpressionProgram::Op::Cosh:
        case ExpressionProgram::Op::Tanh:
        case ExpressionProgram::Op::Exp:
        case ExpressionProgram::Op::Log:
        case ExpressionProgram::Op::Log10:
        case ExpressionProgram::Op::Sqrt:
        case ExpressionProgram::Op::Abs:
        case ExpressionProgram::Op::Norm:
        case ExpressionProgram::Op::Arg:
        case ExpressionProgram::Op::Conjugate:
        case ExpressionProgram::Op::Real:
        case ExpressionProgram::Op::Imaginary: unary(instruction); break;
        case ExpressionProgram::Op::Add:
        case ExpressionProgram::Op::Subtract:
        case ExpressionProgram::Op::Multiply:
        case ExpressionProgram::Op::Divide:
        case ExpressionProgram::Op::Power:
        case ExpressionProgram::Op::MakeComplex:
        case ExpressionProgram::Op::Polar: binary(instruction); break;
        }
    }
    return {stack[0].base, stack[0].delta, stack[0].status, stack[0].denominatorInstability};
}

bool ExpressionCenteredEvaluator::evaluate4(const ExpressionProgram& program, const ExpressionContext& reference, const ExpressionDeltaContext* inputDeltas, ExpressionCenteredResult* results) {
    return evaluate4WithNodeBasesImpl(program, reference, nullptr, inputDeltas, results);
}

bool ExpressionCenteredEvaluator::evaluate4WithNodeBases(const ExpressionProgram& program, const ExpressionContext& reference, const Complex* nodeBases, const ExpressionDeltaContext* inputDeltas, ExpressionCenteredResult* results) {
    return evaluate4WithNodeBasesImpl(program, reference, nodeBases, inputDeltas, results);
}

bool ExpressionCenteredEvaluator::evaluate4FastEntire(const ExpressionProgram& program, const ExpressionContext& reference, const ExpressionDeltaContext* inputDeltas, ExpressionCenteredResult* results) {
    return evaluate4FastEntireImpl(program, reference, nullptr, nullptr, inputDeltas, results, nullptr);
}

bool ExpressionCenteredEvaluator::evaluate4WithNodeBasesFastEntire(const ExpressionProgram& program, const ExpressionContext& reference, const Complex* nodeBases, const Complex* nodeAuxiliaries, const ExpressionDeltaContext* inputDeltas, ExpressionCenteredResult* results) {
    return evaluate4FastEntireImpl(program, reference, nodeBases, nodeAuxiliaries, inputDeltas, results, nullptr);
}

bool ExpressionCenteredEvaluator::evaluate4WithNodeBasesFastEntireDerivative(const ExpressionProgram& program, const ExpressionContext& reference, const Complex* nodeBases, const Complex* nodeAuxiliaries, const ExpressionDeltaContext* inputDeltas, ExpressionCenteredResult* results, double* jacobianNorms) {
    return evaluate4FastEntireImpl(program, reference, nodeBases, nodeAuxiliaries, inputDeltas, results, jacobianNorms);
}

bool ExpressionCenteredEvaluator::evaluate4WithLaneNodeBasesFastEntire(const ExpressionProgram& program, const ExpressionContext* const* references, const Complex* const* nodeBases, const Complex* const* nodeAuxiliaries, const ExpressionDeltaContext* const* inputDeltas, uint8_t activeMask, ExpressionCenteredResult* results) {
    return evaluate4LaneFastEntireImpl(program, references, nodeBases, nodeAuxiliaries, inputDeltas, activeMask, results, nullptr);
}

bool ExpressionCenteredEvaluator::evaluate4WithLaneNodeBasesFastEntireDerivative(const ExpressionProgram& program, const ExpressionContext* const* references, const Complex* const* nodeBases, const Complex* const* nodeAuxiliaries, const ExpressionDeltaContext* const* inputDeltas, uint8_t activeMask, ExpressionCenteredResult* results, double* jacobianNorms) {
    return evaluate4LaneFastEntireImpl(program, references, nodeBases, nodeAuxiliaries, inputDeltas, activeMask, results, jacobianNorms);
}

bool ExpressionCenteredEvaluator::evaluate4FastEntireImpl(const ExpressionProgram& program, const ExpressionContext& reference, const Complex* nodeBases, const Complex* nodeAuxiliaries, const ExpressionDeltaContext* inputDeltas, ExpressionCenteredResult* results, double* derivatives) {
    if (!program._valid || !inputDeltas || !results || (derivatives && !nodeAuxiliaries)) return false;
    for (const ExpressionProgram::Instruction& instruction : program._code) {
        switch (instruction.op) {
        case ExpressionProgram::Op::Constant:
        case ExpressionProgram::Op::Z:
        case ExpressionProgram::Op::C:
        case ExpressionProgram::Op::Z0:
        case ExpressionProgram::Op::Iteration:
        case ExpressionProgram::Op::Parameter:
        case ExpressionProgram::Op::Negate:
        case ExpressionProgram::Op::Add:
        case ExpressionProgram::Op::Subtract:
        case ExpressionProgram::Op::Multiply:
        case ExpressionProgram::Op::Divide:
        case ExpressionProgram::Op::Square:
        case ExpressionProgram::Op::Sin:
        case ExpressionProgram::Op::Cos:
        case ExpressionProgram::Op::Tan:
        case ExpressionProgram::Op::Sinh:
        case ExpressionProgram::Op::Cosh:
        case ExpressionProgram::Op::Tanh:
        case ExpressionProgram::Op::Exp:
        case ExpressionProgram::Op::Log:
        case ExpressionProgram::Op::Log10:
        case ExpressionProgram::Op::Sqrt:
        case ExpressionProgram::Op::Power:
        case ExpressionProgram::Op::Abs:
        case ExpressionProgram::Op::Norm:
        case ExpressionProgram::Op::Conjugate:
        case ExpressionProgram::Op::Real:
        case ExpressionProgram::Op::Imaginary:
        case ExpressionProgram::Op::MakeComplex: break;
        default: return derivatives ? false : evaluate4WithNodeBasesImpl(program, reference, nodeBases, inputDeltas, results);
        }
    }

    struct FastBatchValue {
        Complex base{};
        __m256d deltaReal = _mm256_setzero_pd();
        __m256d deltaImaginary = _mm256_setzero_pd();
        __m256d derivativeReal = _mm256_setzero_pd();
        __m256d derivativeImaginary = _mm256_setzero_pd();
        __m256d derivativeYReal = _mm256_setzero_pd();
        __m256d derivativeYImaginary = _mm256_setzero_pd();
        std::array<ExpressionCenteredStatus, 4> status{};
        std::array<bool, 4> denominatorInstability{};
    };
    alignas(32) std::array<FastBatchValue, ExpressionProgram::MAX_STACK> stack;
    size_t top = 0;
    const __m256d signMask = _mm256_set1_pd(-0.0);
    const __m256d absoluteMask = _mm256_castsi256_pd(_mm256_set1_epi64x(0x7fffffffffffffffLL));
    const __m256d maximum = _mm256_set1_pd(std::numeric_limits<double>::max());
    auto loadReal = [&](auto getter) { return _mm256_set_pd(getter(inputDeltas[3]).real(), getter(inputDeltas[2]).real(), getter(inputDeltas[1]).real(), getter(inputDeltas[0]).real()); };
    auto loadImaginary = [&](auto getter) { return _mm256_set_pd(getter(inputDeltas[3]).imag(), getter(inputDeltas[2]).imag(), getter(inputDeltas[1]).imag(), getter(inputDeltas[0]).imag()); };
    auto multiply = [](__m256d leftReal, __m256d leftImaginary, __m256d rightReal, __m256d rightImaginary, __m256d& outputReal, __m256d& outputImaginary) {
        outputReal = _mm256_sub_pd(_mm256_mul_pd(leftReal, rightReal), _mm256_mul_pd(leftImaginary, rightImaginary));
        outputImaginary = _mm256_add_pd(_mm256_mul_pd(leftReal, rightImaginary), _mm256_mul_pd(leftImaginary, rightReal));
    };
    auto normalize = [&](FastBatchValue& value) {
        bool allSuccessful = true;
        for (ExpressionCenteredStatus status : value.status) allSuccessful = allSuccessful && status == ExpressionCenteredStatus::Success;
        if (allSuccessful && finite(value.base)) {
            __m256d usable = _mm256_cmp_pd(_mm256_and_pd(value.deltaReal, absoluteMask), maximum, _CMP_LE_OQ);
            usable = _mm256_and_pd(usable, _mm256_cmp_pd(_mm256_and_pd(value.deltaImaginary, absoluteMask), maximum, _CMP_LE_OQ));
            const __m256d endpointReal = _mm256_add_pd(_mm256_set1_pd(value.base.real()), value.deltaReal);
            const __m256d endpointImaginary = _mm256_add_pd(_mm256_set1_pd(value.base.imag()), value.deltaImaginary);
            usable = _mm256_and_pd(usable, _mm256_cmp_pd(_mm256_and_pd(endpointReal, absoluteMask), maximum, _CMP_LE_OQ));
            usable = _mm256_and_pd(usable, _mm256_cmp_pd(_mm256_and_pd(endpointImaginary, absoluteMask), maximum, _CMP_LE_OQ));
            if (_mm256_movemask_pd(usable) == 15) return;
        }
        alignas(32) double real[4], imaginary[4];
        _mm256_store_pd(real, value.deltaReal);
        _mm256_store_pd(imaginary, value.deltaImaginary);
        for (int lane = 0; lane < 4; ++lane) {
            Complex delta{real[lane], imaginary[lane]};
            if (value.status[lane] == ExpressionCenteredStatus::Success && (!finite(value.base) || !finite(delta) || !finite(value.base + delta))) value.status[lane] = ExpressionCenteredStatus::NonFinite;
            if (value.status[lane] != ExpressionCenteredStatus::Success) {
                delta = nanComplex();
                real[lane] = delta.real();
                imaginary[lane] = delta.imag();
            }
        }
        value.deltaReal = _mm256_load_pd(real);
        value.deltaImaginary = _mm256_load_pd(imaginary);
    };
    auto pushUnchanged = [&](Complex value) {
        FastBatchValue& output = stack[top++];
        output.base = value;
        output.deltaReal = _mm256_setzero_pd();
        output.deltaImaginary = _mm256_setzero_pd();
        output.derivativeReal = _mm256_setzero_pd();
        output.derivativeImaginary = _mm256_setzero_pd();
        output.derivativeYReal = _mm256_setzero_pd();
        output.derivativeYImaginary = _mm256_setzero_pd();
        output.status.fill(finite(value) ? ExpressionCenteredStatus::Success : ExpressionCenteredStatus::NonFinite);
        output.denominatorInstability.fill(false);
    };
    auto pushPerturbed = [&](Complex base, auto getter, bool derivativeSeed) {
        FastBatchValue& output = stack[top++];
        output.base = base;
        output.deltaReal = loadReal(getter);
        output.deltaImaginary = loadImaginary(getter);
        output.derivativeReal = derivativeSeed ? _mm256_set1_pd(1.0) : _mm256_setzero_pd();
        output.derivativeImaginary = _mm256_setzero_pd();
        output.derivativeYReal = _mm256_setzero_pd();
        output.derivativeYImaginary = derivativeSeed ? _mm256_set1_pd(1.0) : _mm256_setzero_pd();
        output.status.fill(ExpressionCenteredStatus::Success);
        output.denominatorInstability.fill(false);
        normalize(output);
    };
    auto scalarEntireLane = [&](ExpressionProgram::Op operation, Complex inputBase, Complex inputDelta, Complex outputBase, Complex& outputDelta, Complex& derivativeFactor) {
        const Complex inputEndpoint = inputBase + inputDelta;
        const Complex outputEndpoint = ExpressionProgram::evaluateUnary(operation, inputEndpoint);
        outputDelta = nanComplex();
        derivativeFactor = nanComplex();
        ExpressionCenteredStatus status = ExpressionCenteredStatus::Success;
        switch (operation) {
        case ExpressionProgram::Op::Sin:
            outputDelta = centeredSine(inputBase, inputDelta);
            derivativeFactor = std::cos(inputEndpoint);
            break;
        case ExpressionProgram::Op::Cos:
            outputDelta = centeredCosine(inputBase, inputDelta);
            derivativeFactor = -std::sin(inputEndpoint);
            break;
        case ExpressionProgram::Op::Tan:
            status = centeredTangent(inputBase, inputDelta, inputEndpoint, outputBase, outputEndpoint, outputDelta);
            if (status == ExpressionCenteredStatus::Success) derivativeFactor = Complex{1.0, 0.0} + outputEndpoint * outputEndpoint;
            break;
        case ExpressionProgram::Op::Sinh:
            outputDelta = centeredHyperbolicSine(inputBase, inputDelta);
            derivativeFactor = std::cosh(inputEndpoint);
            break;
        case ExpressionProgram::Op::Cosh:
            outputDelta = centeredHyperbolicCosine(inputBase, inputDelta);
            derivativeFactor = std::sinh(inputEndpoint);
            break;
        case ExpressionProgram::Op::Tanh:
            status = centeredHyperbolicTangent(inputBase, inputDelta, inputEndpoint, outputBase, outputEndpoint, outputDelta);
            if (status == ExpressionCenteredStatus::Success) derivativeFactor = Complex{1.0, 0.0} - outputEndpoint * outputEndpoint;
            break;
        case ExpressionProgram::Op::Exp:
            if (!finite(outputBase) || !finite(outputEndpoint))
                status = ExpressionCenteredStatus::NonFinite;
            else {
                outputDelta = relativeExponentialDelta(outputBase, outputEndpoint, inputDelta);
                derivativeFactor = outputEndpoint;
            }
            break;
        case ExpressionProgram::Op::Log:
            status = centeredLogarithm(inputBase, inputDelta, inputEndpoint, outputDelta);
            if (status == ExpressionCenteredStatus::Success) derivativeFactor = Complex{1.0, 0.0} / inputEndpoint;
            break;
        case ExpressionProgram::Op::Log10:
            status = centeredLogarithm(inputBase, inputDelta, inputEndpoint, outputDelta);
            if (status == ExpressionCenteredStatus::Success) {
                outputDelta /= std::log(10.0);
                derivativeFactor = (Complex{1.0, 0.0} / inputEndpoint) / std::log(10.0);
            }
            break;
        case ExpressionProgram::Op::Sqrt:
            status = centeredSquareRoot(inputBase, inputDelta, inputEndpoint, outputBase, outputEndpoint, outputDelta);
            if (status == ExpressionCenteredStatus::Success) derivativeFactor = Complex{0.5, 0.0} / outputEndpoint;
            break;
        case ExpressionProgram::Op::Abs:
            status = centeredAbsolute(inputBase, inputDelta, inputEndpoint, outputDelta);
            if (status == ExpressionCenteredStatus::Success) {
                const double mag = std::abs(inputEndpoint);
                derivativeFactor = mag > 0.0 ? inputEndpoint / mag : Complex{0.0, 0.0};
            }
            break;
        default: return ExpressionCenteredStatus::Unsupported;
        }
        return finish(outputBase, outputDelta, outputEndpoint, status).status;
    };

    for (size_t instructionIndex = 0; instructionIndex < program._code.size(); ++instructionIndex) {
        const ExpressionProgram::Instruction& instruction = program._code[instructionIndex];
        const Complex sharedBase = nodeBases ? nodeBases[instructionIndex] : Complex{};
        switch (instruction.op) {
        case ExpressionProgram::Op::Constant: pushUnchanged(nodeBases ? sharedBase : instruction.value); break;
        case ExpressionProgram::Op::Z: pushPerturbed(nodeBases ? sharedBase : reference.z, [](const ExpressionDeltaContext& delta) { return delta.z; }, true); break;
        case ExpressionProgram::Op::C: pushPerturbed(nodeBases ? sharedBase : reference.c, [](const ExpressionDeltaContext& delta) { return delta.c; }, false); break;
        case ExpressionProgram::Op::Z0: pushPerturbed(nodeBases ? sharedBase : reference.z0, [](const ExpressionDeltaContext& delta) { return delta.z0; }, false); break;
        case ExpressionProgram::Op::Iteration: pushUnchanged(nodeBases ? sharedBase : Complex{(double)reference.iteration, 0.0}); break;
        case ExpressionProgram::Op::Parameter: {
            const int parameter = instruction.argument;
            pushPerturbed(nodeBases ? sharedBase : reference.parameters[parameter], [parameter](const ExpressionDeltaContext& delta) { return delta.parameters[parameter]; }, false);
            break;
        }
        case ExpressionProgram::Op::Negate: {
            FastBatchValue& value = stack[top - 1];
            value.base = nodeBases ? sharedBase : -value.base;
            value.deltaReal = _mm256_xor_pd(value.deltaReal, signMask);
            value.deltaImaginary = _mm256_xor_pd(value.deltaImaginary, signMask);
            if (derivatives) {
                value.derivativeReal = _mm256_xor_pd(value.derivativeReal, signMask);
                value.derivativeImaginary = _mm256_xor_pd(value.derivativeImaginary, signMask);
                value.derivativeYReal = _mm256_xor_pd(value.derivativeYReal, signMask);
                value.derivativeYImaginary = _mm256_xor_pd(value.derivativeYImaginary, signMask);
            }
            normalize(value);
            break;
        }
        case ExpressionProgram::Op::Square: {
            FastBatchValue input = stack[top - 1];
            FastBatchValue& output = stack[top - 1];
            output.base = nodeBases ? sharedBase : input.base * input.base;
            __m256d firstReal, firstImaginary, squaredReal, squaredImaginary;
            multiply(_mm256_set1_pd(input.base.real()), _mm256_set1_pd(input.base.imag()), input.deltaReal, input.deltaImaginary, firstReal, firstImaginary);
            multiply(input.deltaReal, input.deltaImaginary, input.deltaReal, input.deltaImaginary, squaredReal, squaredImaginary);
            output.deltaReal = _mm256_add_pd(_mm256_add_pd(firstReal, firstReal), squaredReal);
            output.deltaImaginary = _mm256_add_pd(_mm256_add_pd(firstImaginary, firstImaginary), squaredImaginary);
            if (derivatives) {
                const __m256d endpointReal = _mm256_add_pd(_mm256_set1_pd(input.base.real()), input.deltaReal);
                const __m256d endpointImaginary = _mm256_add_pd(_mm256_set1_pd(input.base.imag()), input.deltaImaginary);
                multiply(endpointReal, endpointImaginary, input.derivativeReal, input.derivativeImaginary, output.derivativeReal, output.derivativeImaginary);
                output.derivativeReal = _mm256_add_pd(output.derivativeReal, output.derivativeReal);
                output.derivativeImaginary = _mm256_add_pd(output.derivativeImaginary, output.derivativeImaginary);
                multiply(endpointReal, endpointImaginary, input.derivativeYReal, input.derivativeYImaginary, output.derivativeYReal, output.derivativeYImaginary);
                output.derivativeYReal = _mm256_add_pd(output.derivativeYReal, output.derivativeYReal);
                output.derivativeYImaginary = _mm256_add_pd(output.derivativeYImaginary, output.derivativeYImaginary);
            }
            normalize(output);
            break;
        }
        case ExpressionProgram::Op::Sin:
        case ExpressionProgram::Op::Cos:
        case ExpressionProgram::Op::Tan:
        case ExpressionProgram::Op::Sinh:
        case ExpressionProgram::Op::Cosh:
        case ExpressionProgram::Op::Tanh:
        case ExpressionProgram::Op::Exp:
        case ExpressionProgram::Op::Log:
        case ExpressionProgram::Op::Log10:
        case ExpressionProgram::Op::Sqrt:
        case ExpressionProgram::Op::Abs: {
            FastBatchValue input = stack[top - 1];
            FastBatchValue& output = stack[top - 1];
            output.base = nodeBases ? sharedBase : ExpressionProgram::evaluateUnary(instruction.op, input.base);
            FastEntireOperation operation = FastEntireOperation::Sin;
            switch (instruction.op) {
            case ExpressionProgram::Op::Cos: operation = FastEntireOperation::Cos; break;
            case ExpressionProgram::Op::Sinh: operation = FastEntireOperation::Sinh; break;
            case ExpressionProgram::Op::Cosh: operation = FastEntireOperation::Cosh; break;
            case ExpressionProgram::Op::Exp: operation = FastEntireOperation::Exp; break;
            default: break;
            }
            const Complex* companion = nodeAuxiliaries ? nodeAuxiliaries + instructionIndex : nullptr;
            uint8_t validMask = 0;
            for (int lane = 0; lane < 4; ++lane)
                if (input.status[lane] == ExpressionCenteredStatus::Success) validMask |= uint8_t{1} << lane;
            __m256d derivativeFactorReal = _mm256_setzero_pd();
            __m256d derivativeFactorImaginary = _mm256_setzero_pd();
            const bool meromorphicOrBranchOperation = instruction.op == ExpressionProgram::Op::Tan || instruction.op == ExpressionProgram::Op::Tanh || instruction.op == ExpressionProgram::Op::Log || instruction.op == ExpressionProgram::Op::Log10 || instruction.op == ExpressionProgram::Op::Sqrt || instruction.op == ExpressionProgram::Op::Abs;
            const uint8_t fastMask = meromorphicOrBranchOperation ? 0 : fusedCenteredEntireVectors(operation, input.base, output.base, companion, input.deltaReal, input.deltaImaginary, validMask, output.deltaReal, output.deltaImaginary, derivatives ? &derivativeFactorReal : nullptr, derivatives ? &derivativeFactorImaginary : nullptr);
            output.status = input.status;
            output.denominatorInstability = input.denominatorInstability;
            if (validMask != 0 && fastMask == validMask) {
                if (derivatives) {
                    multiply(derivativeFactorReal, derivativeFactorImaginary, input.derivativeReal, input.derivativeImaginary, output.derivativeReal, output.derivativeImaginary);
                    multiply(derivativeFactorReal, derivativeFactorImaginary, input.derivativeYReal, input.derivativeYImaginary, output.derivativeYReal, output.derivativeYImaginary);
                }
                normalize(output);
                break;
            }
            alignas(32) double inputReal[4], inputImaginary[4], outputReal[4], outputImaginary[4], factorReal[4], factorImaginary[4];
            _mm256_store_pd(inputReal, input.deltaReal);
            _mm256_store_pd(inputImaginary, input.deltaImaginary);
            _mm256_store_pd(outputReal, output.deltaReal);
            _mm256_store_pd(outputImaginary, output.deltaImaginary);
            if (derivatives) {
                _mm256_store_pd(factorReal, derivativeFactorReal);
                _mm256_store_pd(factorImaginary, derivativeFactorImaginary);
            }
            for (int lane = 0; lane < 4; ++lane) {
                if (input.status[lane] != ExpressionCenteredStatus::Success) {
                    const Complex nan = nanComplex();
                    outputReal[lane] = nan.real();
                    outputImaginary[lane] = nan.imag();
                    if (derivatives) {
                        factorReal[lane] = nan.real();
                        factorImaginary[lane] = nan.imag();
                    }
                } else if ((fastMask & (uint8_t{1} << lane)) == 0) {
                    Complex outputDelta;
                    Complex derivativeFactor;
                    output.status[lane] = scalarEntireLane(instruction.op, input.base, {inputReal[lane], inputImaginary[lane]}, output.base, outputDelta, derivativeFactor);
                    if (meromorphicOrBranchOperation && output.status[lane] != ExpressionCenteredStatus::Success) output.denominatorInstability[lane] = true;
                    outputReal[lane] = outputDelta.real();
                    outputImaginary[lane] = outputDelta.imag();
                    if (derivatives) {
                        factorReal[lane] = derivativeFactor.real();
                        factorImaginary[lane] = derivativeFactor.imag();
                    }
                }
            }
            output.deltaReal = _mm256_load_pd(outputReal);
            output.deltaImaginary = _mm256_load_pd(outputImaginary);
            if (derivatives) {
                derivativeFactorReal = _mm256_load_pd(factorReal);
                derivativeFactorImaginary = _mm256_load_pd(factorImaginary);
                multiply(derivativeFactorReal, derivativeFactorImaginary, input.derivativeReal, input.derivativeImaginary, output.derivativeReal, output.derivativeImaginary);
                multiply(derivativeFactorReal, derivativeFactorImaginary, input.derivativeYReal, input.derivativeYImaginary, output.derivativeYReal, output.derivativeYImaginary);
            }
            normalize(output);
            break;
        }
        case ExpressionProgram::Op::Norm:
        case ExpressionProgram::Op::Conjugate:
        case ExpressionProgram::Op::Real:
        case ExpressionProgram::Op::Imaginary: {
            FastBatchValue input = stack[top - 1];
            FastBatchValue& output = stack[top - 1];
            output.base = nodeBases ? sharedBase : ExpressionProgram::evaluateUnary(instruction.op, input.base);
            if (instruction.op == ExpressionProgram::Op::Conjugate) {
                output.deltaReal = input.deltaReal;
                output.deltaImaginary = _mm256_xor_pd(input.deltaImaginary, signMask);
                if (derivatives) {
                    output.derivativeReal = input.derivativeReal;
                    output.derivativeImaginary = _mm256_xor_pd(input.derivativeImaginary, signMask);
                    output.derivativeYReal = input.derivativeYReal;
                    output.derivativeYImaginary = _mm256_xor_pd(input.derivativeYImaginary, signMask);
                }
            } else if (instruction.op == ExpressionProgram::Op::Real) {
                output.deltaReal = input.deltaReal;
                output.deltaImaginary = _mm256_setzero_pd();
                if (derivatives) {
                    output.derivativeReal = input.derivativeReal;
                    output.derivativeImaginary = _mm256_setzero_pd();
                    output.derivativeYReal = input.derivativeYReal;
                    output.derivativeYImaginary = _mm256_setzero_pd();
                }
            } else if (instruction.op == ExpressionProgram::Op::Imaginary) {
                output.deltaReal = input.deltaImaginary;
                output.deltaImaginary = _mm256_setzero_pd();
                if (derivatives) {
                    output.derivativeReal = input.derivativeImaginary;
                    output.derivativeImaginary = _mm256_setzero_pd();
                    output.derivativeYReal = input.derivativeYImaginary;
                    output.derivativeYImaginary = _mm256_setzero_pd();
                }
            } else {
                const __m256d baseReal = _mm256_set1_pd(input.base.real());
                const __m256d baseImaginary = _mm256_set1_pd(input.base.imag());
                const __m256d endpointReal = _mm256_add_pd(baseReal, input.deltaReal);
                const __m256d endpointImaginary = _mm256_add_pd(baseImaginary, input.deltaImaginary);
                const __m256d cross = _mm256_add_pd(_mm256_mul_pd(baseReal, input.deltaReal), _mm256_mul_pd(baseImaginary, input.deltaImaginary));
                const __m256d deltaNorm = _mm256_add_pd(_mm256_mul_pd(input.deltaReal, input.deltaReal), _mm256_mul_pd(input.deltaImaginary, input.deltaImaginary));
                output.deltaReal = _mm256_add_pd(_mm256_add_pd(cross, cross), deltaNorm);
                output.deltaImaginary = _mm256_setzero_pd();
                if (derivatives) {
                    const __m256d twiceEndpointReal = _mm256_add_pd(endpointReal, endpointReal);
                    const __m256d twiceEndpointImaginary = _mm256_add_pd(endpointImaginary, endpointImaginary);
                    output.derivativeReal = _mm256_add_pd(_mm256_mul_pd(twiceEndpointReal, input.derivativeReal), _mm256_mul_pd(twiceEndpointImaginary, input.derivativeImaginary));
                    output.derivativeImaginary = _mm256_setzero_pd();
                    output.derivativeYReal = _mm256_add_pd(_mm256_mul_pd(twiceEndpointReal, input.derivativeYReal), _mm256_mul_pd(twiceEndpointImaginary, input.derivativeYImaginary));
                    output.derivativeYImaginary = _mm256_setzero_pd();
                }
            }
            output.status = input.status;
            normalize(output);
            break;
        }
        case ExpressionProgram::Op::Add:
        case ExpressionProgram::Op::Subtract:
        case ExpressionProgram::Op::Multiply:
        case ExpressionProgram::Op::Divide:
        case ExpressionProgram::Op::Power:
        case ExpressionProgram::Op::MakeComplex: {
            FastBatchValue right = stack[--top];
            FastBatchValue left = stack[top - 1];
            FastBatchValue& output = stack[top - 1];
            output.base = nodeBases ? sharedBase : ExpressionProgram::evaluateBinary(instruction.op, left.base, right.base);
            for (int lane = 0; lane < 4; ++lane) {
                output.status[lane] = combine(left.status[lane], right.status[lane]);
                output.denominatorInstability[lane] = left.denominatorInstability[lane] || right.denominatorInstability[lane];
            }
            if (instruction.op == ExpressionProgram::Op::MakeComplex) {
                output.deltaReal = left.deltaReal;
                output.deltaImaginary = right.deltaReal;
                if (derivatives) {
                    output.derivativeReal = left.derivativeReal;
                    output.derivativeImaginary = right.derivativeReal;
                    output.derivativeYReal = left.derivativeYReal;
                    output.derivativeYImaginary = right.derivativeYReal;
                }
            } else if (instruction.op == ExpressionProgram::Op::Add) {
                output.deltaReal = _mm256_add_pd(left.deltaReal, right.deltaReal);
                output.deltaImaginary = _mm256_add_pd(left.deltaImaginary, right.deltaImaginary);
                if (derivatives) {
                    output.derivativeReal = _mm256_add_pd(left.derivativeReal, right.derivativeReal);
                    output.derivativeImaginary = _mm256_add_pd(left.derivativeImaginary, right.derivativeImaginary);
                    output.derivativeYReal = _mm256_add_pd(left.derivativeYReal, right.derivativeYReal);
                    output.derivativeYImaginary = _mm256_add_pd(left.derivativeYImaginary, right.derivativeYImaginary);
                }
            } else if (instruction.op == ExpressionProgram::Op::Subtract) {
                output.deltaReal = _mm256_sub_pd(left.deltaReal, right.deltaReal);
                output.deltaImaginary = _mm256_sub_pd(left.deltaImaginary, right.deltaImaginary);
                if (derivatives) {
                    output.derivativeReal = _mm256_sub_pd(left.derivativeReal, right.derivativeReal);
                    output.derivativeImaginary = _mm256_sub_pd(left.derivativeImaginary, right.derivativeImaginary);
                    output.derivativeYReal = _mm256_sub_pd(left.derivativeYReal, right.derivativeYReal);
                    output.derivativeYImaginary = _mm256_sub_pd(left.derivativeYImaginary, right.derivativeYImaginary);
                }
            } else if (instruction.op == ExpressionProgram::Op::Divide) {
                alignas(32) double leftDeltaReal[4], leftDeltaImaginary[4], rightDeltaReal[4], rightDeltaImaginary[4], outputReal[4], outputImaginary[4];
                alignas(32) double leftDerivativeReal[4], leftDerivativeImaginary[4], leftDerivativeYReal[4], leftDerivativeYImaginary[4];
                alignas(32) double rightDerivativeReal[4], rightDerivativeImaginary[4], rightDerivativeYReal[4], rightDerivativeYImaginary[4];
                alignas(32) double outputDerivativeReal[4]{}, outputDerivativeImaginary[4]{}, outputDerivativeYReal[4]{}, outputDerivativeYImaginary[4]{};
                _mm256_store_pd(leftDeltaReal, left.deltaReal);
                _mm256_store_pd(leftDeltaImaginary, left.deltaImaginary);
                _mm256_store_pd(rightDeltaReal, right.deltaReal);
                _mm256_store_pd(rightDeltaImaginary, right.deltaImaginary);
                if (derivatives) {
                    _mm256_store_pd(leftDerivativeReal, left.derivativeReal);
                    _mm256_store_pd(leftDerivativeImaginary, left.derivativeImaginary);
                    _mm256_store_pd(leftDerivativeYReal, left.derivativeYReal);
                    _mm256_store_pd(leftDerivativeYImaginary, left.derivativeYImaginary);
                    _mm256_store_pd(rightDerivativeReal, right.derivativeReal);
                    _mm256_store_pd(rightDerivativeImaginary, right.derivativeImaginary);
                    _mm256_store_pd(rightDerivativeYReal, right.derivativeYReal);
                    _mm256_store_pd(rightDerivativeYImaginary, right.derivativeYImaginary);
                }
                for (int lane = 0; lane < 4; ++lane) {
                    Complex delta = nanComplex();
                    if (output.status[lane] == ExpressionCenteredStatus::Success) {
                        const Complex leftDelta{leftDeltaReal[lane], leftDeltaImaginary[lane]};
                        const Complex rightDelta{rightDeltaReal[lane], rightDeltaImaginary[lane]};
                        const Complex leftEndpoint = left.base + leftDelta;
                        const Complex rightEndpoint = right.base + rightDelta;
                        output.status[lane] = centeredDivisionDelta(leftDelta, right.base, rightDelta, rightEndpoint, output.base, delta, true);
                        if (output.status[lane] != ExpressionCenteredStatus::Success) output.denominatorInstability[lane] = true;
                        if (derivatives && output.status[lane] == ExpressionCenteredStatus::Success) {
                            const Complex denominatorSquared = rightEndpoint * rightEndpoint;
                            const Complex derivative = (Complex{leftDerivativeReal[lane], leftDerivativeImaginary[lane]} * rightEndpoint - leftEndpoint * Complex{rightDerivativeReal[lane], rightDerivativeImaginary[lane]}) / denominatorSquared;
                            const Complex derivativeY = (Complex{leftDerivativeYReal[lane], leftDerivativeYImaginary[lane]} * rightEndpoint - leftEndpoint * Complex{rightDerivativeYReal[lane], rightDerivativeYImaginary[lane]}) / denominatorSquared;
                            if (!finite(derivative) || !finite(derivativeY)) {
                                output.status[lane] = ExpressionCenteredStatus::NonFinite;
                                output.denominatorInstability[lane] = true;
                            } else {
                                outputDerivativeReal[lane] = derivative.real();
                                outputDerivativeImaginary[lane] = derivative.imag();
                                outputDerivativeYReal[lane] = derivativeY.real();
                                outputDerivativeYImaginary[lane] = derivativeY.imag();
                            }
                        }
                    }
                    outputReal[lane] = delta.real();
                    outputImaginary[lane] = delta.imag();
                }
                output.deltaReal = _mm256_load_pd(outputReal);
                output.deltaImaginary = _mm256_load_pd(outputImaginary);
                if (derivatives) {
                    output.derivativeReal = _mm256_load_pd(outputDerivativeReal);
                    output.derivativeImaginary = _mm256_load_pd(outputDerivativeImaginary);
                    output.derivativeYReal = _mm256_load_pd(outputDerivativeYReal);
                    output.derivativeYImaginary = _mm256_load_pd(outputDerivativeYImaginary);
                }
            } else if (instruction.op == ExpressionProgram::Op::Power) {
                alignas(32) double leftDeltaReal[4], leftDeltaImaginary[4], rightDeltaReal[4], rightDeltaImaginary[4], outputReal[4], outputImaginary[4];
                alignas(32) double leftDerivativeReal[4], leftDerivativeImaginary[4], leftDerivativeYReal[4], leftDerivativeYImaginary[4];
                alignas(32) double rightDerivativeReal[4], rightDerivativeImaginary[4], rightDerivativeYReal[4], rightDerivativeYImaginary[4];
                alignas(32) double outputDerivativeReal[4]{}, outputDerivativeImaginary[4]{}, outputDerivativeYReal[4]{}, outputDerivativeYImaginary[4]{};
                _mm256_store_pd(leftDeltaReal, left.deltaReal);
                _mm256_store_pd(leftDeltaImaginary, left.deltaImaginary);
                _mm256_store_pd(rightDeltaReal, right.deltaReal);
                _mm256_store_pd(rightDeltaImaginary, right.deltaImaginary);
                if (derivatives) {
                    _mm256_store_pd(leftDerivativeReal, left.derivativeReal);
                    _mm256_store_pd(leftDerivativeImaginary, left.derivativeImaginary);
                    _mm256_store_pd(leftDerivativeYReal, left.derivativeYReal);
                    _mm256_store_pd(leftDerivativeYImaginary, left.derivativeYImaginary);
                    _mm256_store_pd(rightDerivativeReal, right.derivativeReal);
                    _mm256_store_pd(rightDerivativeImaginary, right.derivativeImaginary);
                    _mm256_store_pd(rightDerivativeYReal, right.derivativeYReal);
                    _mm256_store_pd(rightDerivativeYImaginary, right.derivativeYImaginary);
                }
                for (int lane = 0; lane < 4; ++lane) {
                    Complex delta = nanComplex();
                    if (output.status[lane] == ExpressionCenteredStatus::Success) {
                        const Complex leftDelta{leftDeltaReal[lane], leftDeltaImaginary[lane]};
                        const Complex rightDelta{rightDeltaReal[lane], rightDeltaImaginary[lane]};
                        const Complex leftEndpoint = left.base + leftDelta;
                        const Complex rightEndpoint = right.base + rightDelta;
                        const Complex outputEndpoint = ExpressionProgram::evaluateBinary(ExpressionProgram::Op::Power, leftEndpoint, rightEndpoint);
                        output.status[lane] = centeredPower(left.base, leftDelta, leftEndpoint, right.base, rightDelta, rightEndpoint, output.base, outputEndpoint, delta);
                        if (derivatives && output.status[lane] == ExpressionCenteredStatus::Success) {
                            if (zero(leftEndpoint)) {
                                outputDerivativeReal[lane] = 0.0;
                                outputDerivativeImaginary[lane] = 0.0;
                                outputDerivativeYReal[lane] = 0.0;
                                outputDerivativeYImaginary[lane] = 0.0;
                            } else {
                                const Complex termU = (rightEndpoint / leftEndpoint) * Complex{leftDerivativeReal[lane], leftDerivativeImaginary[lane]};
                                const Complex termV = std::log(leftEndpoint) * Complex{rightDerivativeReal[lane], rightDerivativeImaginary[lane]};
                                const Complex derivative = outputEndpoint * (termU + termV);

                                const Complex termUY = (rightEndpoint / leftEndpoint) * Complex{leftDerivativeYReal[lane], leftDerivativeYImaginary[lane]};
                                const Complex termVY = std::log(leftEndpoint) * Complex{rightDerivativeYReal[lane], rightDerivativeYImaginary[lane]};
                                const Complex derivativeY = outputEndpoint * (termUY + termVY);

                                if (!finite(derivative) || !finite(derivativeY)) {
                                    output.status[lane] = ExpressionCenteredStatus::NonFinite;
                                } else {
                                    outputDerivativeReal[lane] = derivative.real();
                                    outputDerivativeImaginary[lane] = derivative.imag();
                                    outputDerivativeYReal[lane] = derivativeY.real();
                                    outputDerivativeYImaginary[lane] = derivativeY.imag();
                                }
                            }
                        }
                    }
                    outputReal[lane] = delta.real();
                    outputImaginary[lane] = delta.imag();
                }
                output.deltaReal = _mm256_load_pd(outputReal);
                output.deltaImaginary = _mm256_load_pd(outputImaginary);
                if (derivatives) {
                    output.derivativeReal = _mm256_load_pd(outputDerivativeReal);
                    output.derivativeImaginary = _mm256_load_pd(outputDerivativeImaginary);
                    output.derivativeYReal = _mm256_load_pd(outputDerivativeYReal);
                    output.derivativeYImaginary = _mm256_load_pd(outputDerivativeYImaginary);
                }
            } else if (instruction.op == ExpressionProgram::Op::Multiply) {
                __m256d firstReal, firstImaginary, secondReal, secondImaginary, thirdReal, thirdImaginary;
                multiply(_mm256_set1_pd(left.base.real()), _mm256_set1_pd(left.base.imag()), right.deltaReal, right.deltaImaginary, firstReal, firstImaginary);
                multiply(_mm256_set1_pd(right.base.real()), _mm256_set1_pd(right.base.imag()), left.deltaReal, left.deltaImaginary, secondReal, secondImaginary);
                multiply(left.deltaReal, left.deltaImaginary, right.deltaReal, right.deltaImaginary, thirdReal, thirdImaginary);
                output.deltaReal = _mm256_add_pd(_mm256_add_pd(firstReal, secondReal), thirdReal);
                output.deltaImaginary = _mm256_add_pd(_mm256_add_pd(firstImaginary, secondImaginary), thirdImaginary);
                if (derivatives) {
                    const __m256d leftEndpointReal = _mm256_add_pd(_mm256_set1_pd(left.base.real()), left.deltaReal);
                    const __m256d leftEndpointImaginary = _mm256_add_pd(_mm256_set1_pd(left.base.imag()), left.deltaImaginary);
                    const __m256d rightEndpointReal = _mm256_add_pd(_mm256_set1_pd(right.base.real()), right.deltaReal);
                    const __m256d rightEndpointImaginary = _mm256_add_pd(_mm256_set1_pd(right.base.imag()), right.deltaImaginary);
                    multiply(left.derivativeReal, left.derivativeImaginary, rightEndpointReal, rightEndpointImaginary, firstReal, firstImaginary);
                    multiply(leftEndpointReal, leftEndpointImaginary, right.derivativeReal, right.derivativeImaginary, secondReal, secondImaginary);
                    output.derivativeReal = _mm256_add_pd(firstReal, secondReal);
                    output.derivativeImaginary = _mm256_add_pd(firstImaginary, secondImaginary);
                    multiply(left.derivativeYReal, left.derivativeYImaginary, rightEndpointReal, rightEndpointImaginary, firstReal, firstImaginary);
                    multiply(leftEndpointReal, leftEndpointImaginary, right.derivativeYReal, right.derivativeYImaginary, secondReal, secondImaginary);
                    output.derivativeYReal = _mm256_add_pd(firstReal, secondReal);
                    output.derivativeYImaginary = _mm256_add_pd(firstImaginary, secondImaginary);
                }
            }
            normalize(output);
            break;
        }
        default: return false;
        }
    }
    alignas(32) double resultReal[4], resultImaginary[4], derivativeReal[4], derivativeImaginary[4], derivativeYReal[4], derivativeYImaginary[4];
    _mm256_store_pd(resultReal, stack[0].deltaReal);
    _mm256_store_pd(resultImaginary, stack[0].deltaImaginary);
    if (derivatives) {
        _mm256_store_pd(derivativeReal, stack[0].derivativeReal);
        _mm256_store_pd(derivativeImaginary, stack[0].derivativeImaginary);
        _mm256_store_pd(derivativeYReal, stack[0].derivativeYReal);
        _mm256_store_pd(derivativeYImaginary, stack[0].derivativeYImaginary);
    }
    for (int lane = 0; lane < 4; ++lane) {
        results[lane] = {stack[0].base, {resultReal[lane], resultImaginary[lane]}, stack[0].status[lane], stack[0].denominatorInstability[lane]};
        if (derivatives) {
            const double a = derivativeReal[lane];
            const double b = derivativeYReal[lane];
            const double c = derivativeImaginary[lane];
            const double d = derivativeYImaginary[lane];
            const double scale = std::max({std::abs(a), std::abs(b), std::abs(c), std::abs(d)});
            double norm = 0.0;
            if (scale != 0.0 && std::isfinite(scale)) {
                const double as = a / scale;
                const double bs = b / scale;
                const double cs = c / scale;
                const double ds = d / scale;
                const double sum = as * as + bs * bs + cs * cs + ds * ds;
                const double determinant = as * ds - bs * cs;
                const double discriminant = std::max(0.0, sum * sum - 4.0 * determinant * determinant);
                norm = scale * std::sqrt(0.5 * (sum + std::sqrt(discriminant)));
            } else if (!std::isfinite(scale)) {
                norm = std::numeric_limits<double>::quiet_NaN();
            }
            derivatives[lane] = results[lane].success() && std::isfinite(norm) ? norm : std::numeric_limits<double>::quiet_NaN();
        }
    }
    return true;
}

bool ExpressionCenteredEvaluator::evaluate4LaneFastEntireImpl(const ExpressionProgram& program, const ExpressionContext* const* references, const Complex* const* inputNodeBases, const Complex* const* inputNodeAuxiliaries, const ExpressionDeltaContext* const* inputDeltas, uint8_t activeMask, ExpressionCenteredResult* results, double* derivatives) {
    activeMask &= 15;
    if (!program._valid || !references || !inputNodeBases || !inputDeltas || !results) return false;
    if (activeMask == 0) {
        for (int lane = 0; lane < 4; ++lane) {
            results[lane] = {nanComplex(), nanComplex(), ExpressionCenteredStatus::InvalidProgram};
            if (derivatives) derivatives[lane] = std::numeric_limits<double>::quiet_NaN();
        }
        return true;
    }
    for (const ExpressionProgram::Instruction& instruction : program._code) {
        switch (instruction.op) {
        case ExpressionProgram::Op::Constant:
        case ExpressionProgram::Op::Z:
        case ExpressionProgram::Op::C:
        case ExpressionProgram::Op::Z0:
        case ExpressionProgram::Op::Iteration:
        case ExpressionProgram::Op::Parameter:
        case ExpressionProgram::Op::Negate:
        case ExpressionProgram::Op::Add:
        case ExpressionProgram::Op::Subtract:
        case ExpressionProgram::Op::Multiply:
        case ExpressionProgram::Op::Divide:
        case ExpressionProgram::Op::Square:
        case ExpressionProgram::Op::Sin:
        case ExpressionProgram::Op::Cos:
        case ExpressionProgram::Op::Tan:
        case ExpressionProgram::Op::Sinh:
        case ExpressionProgram::Op::Cosh:
        case ExpressionProgram::Op::Tanh:
        case ExpressionProgram::Op::Exp:
        case ExpressionProgram::Op::Log:
        case ExpressionProgram::Op::Log10:
        case ExpressionProgram::Op::Sqrt:
        case ExpressionProgram::Op::Power:
        case ExpressionProgram::Op::Abs:
        case ExpressionProgram::Op::Norm:
        case ExpressionProgram::Op::Conjugate:
        case ExpressionProgram::Op::Real:
        case ExpressionProgram::Op::Imaginary:
        case ExpressionProgram::Op::MakeComplex: break;
        default: return false;
        }
    }

    int firstActive = 0;
    while ((activeMask & (uint8_t{1} << firstActive)) == 0) ++firstActive;
    std::array<const ExpressionContext*, 4> laneReferences;
    std::array<const ExpressionDeltaContext*, 4> laneDeltas;
    std::array<const Complex*, 4> nodeBases;
    std::array<const Complex*, 4> nodeAuxiliaries;
    const bool haveAuxiliaries = inputNodeAuxiliaries != nullptr;
    if (derivatives && !haveAuxiliaries) return false;
    for (int lane = 0; lane < 4; ++lane) {
        const bool active = (activeMask & (uint8_t{1} << lane)) != 0;
        if (active && (!references[lane] || !inputNodeBases[lane] || !inputDeltas[lane] || (haveAuxiliaries && !inputNodeAuxiliaries[lane]))) return false;
        laneReferences[lane] = references[active ? lane : firstActive];
        laneDeltas[lane] = inputDeltas[active ? lane : firstActive];
        nodeBases[lane] = inputNodeBases[active ? lane : firstActive];
        nodeAuxiliaries[lane] = haveAuxiliaries ? inputNodeAuxiliaries[active ? lane : firstActive] : nullptr;
    }

    auto gatherComponent = [](const std::array<const Complex*, 4>& pointers, size_t index, size_t component) {
        static_assert(sizeof(Complex) == 2 * sizeof(double), "AVX2 node gather requires interleaved complex doubles");
        const double* anchor = reinterpret_cast<const double*>(pointers[0] + index) + component;
        const intptr_t anchorAddress = reinterpret_cast<intptr_t>(anchor);
        const __m256i offsets = _mm256_set_epi64x(reinterpret_cast<intptr_t>(reinterpret_cast<const double*>(pointers[3] + index) + component) - anchorAddress, reinterpret_cast<intptr_t>(reinterpret_cast<const double*>(pointers[2] + index) + component) - anchorAddress, reinterpret_cast<intptr_t>(reinterpret_cast<const double*>(pointers[1] + index) + component) - anchorAddress, 0);
        return _mm256_i64gather_pd(anchor, offsets, 1);
    };
    auto gatherReal = [&](const std::array<const Complex*, 4>& pointers, size_t index) { return gatherComponent(pointers, index, 0); };
    auto gatherImaginary = [&](const std::array<const Complex*, 4>& pointers, size_t index) { return gatherComponent(pointers, index, 1); };

    struct FastLaneBatchValue {
        __m256d baseReal = _mm256_setzero_pd();
        __m256d baseImaginary = _mm256_setzero_pd();
        __m256d deltaReal = _mm256_setzero_pd();
        __m256d deltaImaginary = _mm256_setzero_pd();
        __m256d derivativeReal = _mm256_setzero_pd();
        __m256d derivativeImaginary = _mm256_setzero_pd();
        __m256d derivativeYReal = _mm256_setzero_pd();
        __m256d derivativeYImaginary = _mm256_setzero_pd();
        std::array<ExpressionCenteredStatus, 4> status{};
        std::array<bool, 4> denominatorInstability{};
    };
    alignas(32) std::array<FastLaneBatchValue, ExpressionProgram::MAX_STACK> stack;
    size_t top = 0;
    const __m256d signMask = _mm256_set1_pd(-0.0);
    const __m256d absoluteMask = _mm256_castsi256_pd(_mm256_set1_epi64x(0x7fffffffffffffffLL));
    const __m256d maximum = _mm256_set1_pd(std::numeric_limits<double>::max());
    auto loadReferenceReal = [&](auto getter) { return _mm256_set_pd(getter(*laneReferences[3]).real(), getter(*laneReferences[2]).real(), getter(*laneReferences[1]).real(), getter(*laneReferences[0]).real()); };
    auto loadReferenceImaginary = [&](auto getter) { return _mm256_set_pd(getter(*laneReferences[3]).imag(), getter(*laneReferences[2]).imag(), getter(*laneReferences[1]).imag(), getter(*laneReferences[0]).imag()); };
    auto loadReal = [&](auto getter) { return _mm256_set_pd(getter(*laneDeltas[3]).real(), getter(*laneDeltas[2]).real(), getter(*laneDeltas[1]).real(), getter(*laneDeltas[0]).real()); };
    auto loadImaginary = [&](auto getter) { return _mm256_set_pd(getter(*laneDeltas[3]).imag(), getter(*laneDeltas[2]).imag(), getter(*laneDeltas[1]).imag(), getter(*laneDeltas[0]).imag()); };
    auto multiply = [](__m256d leftReal, __m256d leftImaginary, __m256d rightReal, __m256d rightImaginary, __m256d& outputReal, __m256d& outputImaginary) {
        outputReal = _mm256_sub_pd(_mm256_mul_pd(leftReal, rightReal), _mm256_mul_pd(leftImaginary, rightImaginary));
        outputImaginary = _mm256_add_pd(_mm256_mul_pd(leftReal, rightImaginary), _mm256_mul_pd(leftImaginary, rightReal));
    };
    auto successfulMask = [](const FastLaneBatchValue& value) {
        uint8_t mask = 0;
        for (int lane = 0; lane < 4; ++lane)
            if (value.status[lane] == ExpressionCenteredStatus::Success) mask |= uint8_t{1} << lane;
        return mask;
    };
    auto normalize = [&](FastLaneBatchValue& value) {
        const uint8_t successMask = successfulMask(value);
        __m256d usable = _mm256_cmp_pd(_mm256_and_pd(value.baseReal, absoluteMask), maximum, _CMP_LE_OQ);
        usable = _mm256_and_pd(usable, _mm256_cmp_pd(_mm256_and_pd(value.baseImaginary, absoluteMask), maximum, _CMP_LE_OQ));
        usable = _mm256_and_pd(usable, _mm256_cmp_pd(_mm256_and_pd(value.deltaReal, absoluteMask), maximum, _CMP_LE_OQ));
        usable = _mm256_and_pd(usable, _mm256_cmp_pd(_mm256_and_pd(value.deltaImaginary, absoluteMask), maximum, _CMP_LE_OQ));
        const __m256d endpointReal = _mm256_add_pd(value.baseReal, value.deltaReal);
        const __m256d endpointImaginary = _mm256_add_pd(value.baseImaginary, value.deltaImaginary);
        usable = _mm256_and_pd(usable, _mm256_cmp_pd(_mm256_and_pd(endpointReal, absoluteMask), maximum, _CMP_LE_OQ));
        usable = _mm256_and_pd(usable, _mm256_cmp_pd(_mm256_and_pd(endpointImaginary, absoluteMask), maximum, _CMP_LE_OQ));
        if ((static_cast<uint8_t>(_mm256_movemask_pd(usable)) & successMask) == successMask) return;
        alignas(32) double baseReal[4], baseImaginary[4], deltaReal[4], deltaImaginary[4];
        _mm256_store_pd(baseReal, value.baseReal);
        _mm256_store_pd(baseImaginary, value.baseImaginary);
        _mm256_store_pd(deltaReal, value.deltaReal);
        _mm256_store_pd(deltaImaginary, value.deltaImaginary);
        for (int lane = 0; lane < 4; ++lane) {
            if (value.status[lane] != ExpressionCenteredStatus::Success) continue;
            const Complex base{baseReal[lane], baseImaginary[lane]};
            const Complex delta{deltaReal[lane], deltaImaginary[lane]};
            if (!finite(base) || !finite(delta) || !finite(base + delta)) {
                value.status[lane] = ExpressionCenteredStatus::NonFinite;
                const Complex nan = nanComplex();
                deltaReal[lane] = nan.real();
                deltaImaginary[lane] = nan.imag();
            }
        }
        value.deltaReal = _mm256_load_pd(deltaReal);
        value.deltaImaginary = _mm256_load_pd(deltaImaginary);
    };
    auto setInitialStatus = [&](FastLaneBatchValue& output) {
        for (int lane = 0; lane < 4; ++lane) output.status[lane] = (activeMask & (uint8_t{1} << lane)) != 0 ? ExpressionCenteredStatus::Success : ExpressionCenteredStatus::InvalidProgram;
    };
    auto pushUnchanged = [&](Complex value) {
        FastLaneBatchValue& output = stack[top++];
        output.baseReal = _mm256_set1_pd(value.real());
        output.baseImaginary = _mm256_set1_pd(value.imag());
        output.deltaReal = _mm256_setzero_pd();
        output.deltaImaginary = _mm256_setzero_pd();
        output.derivativeReal = _mm256_setzero_pd();
        output.derivativeImaginary = _mm256_setzero_pd();
        output.derivativeYReal = _mm256_setzero_pd();
        output.derivativeYImaginary = _mm256_setzero_pd();
        setInitialStatus(output);
        output.denominatorInstability.fill(false);
    };
    auto pushIteration = [&] {
        FastLaneBatchValue& output = stack[top++];
        output.baseReal = _mm256_set_pd(static_cast<double>(laneReferences[3]->iteration), static_cast<double>(laneReferences[2]->iteration), static_cast<double>(laneReferences[1]->iteration), static_cast<double>(laneReferences[0]->iteration));
        output.baseImaginary = _mm256_setzero_pd();
        output.deltaReal = _mm256_setzero_pd();
        output.deltaImaginary = _mm256_setzero_pd();
        output.derivativeReal = _mm256_setzero_pd();
        output.derivativeImaginary = _mm256_setzero_pd();
        output.derivativeYReal = _mm256_setzero_pd();
        output.derivativeYImaginary = _mm256_setzero_pd();
        setInitialStatus(output);
        output.denominatorInstability.fill(false);
    };
    auto pushPerturbed = [&](auto referenceGetter, auto deltaGetter, bool derivativeSeed) {
        FastLaneBatchValue& output = stack[top++];
        output.baseReal = loadReferenceReal(referenceGetter);
        output.baseImaginary = loadReferenceImaginary(referenceGetter);
        output.deltaReal = loadReal(deltaGetter);
        output.deltaImaginary = loadImaginary(deltaGetter);
        output.derivativeReal = derivativeSeed ? _mm256_set1_pd(1.0) : _mm256_setzero_pd();
        output.derivativeImaginary = _mm256_setzero_pd();
        output.derivativeYReal = _mm256_setzero_pd();
        output.derivativeYImaginary = derivativeSeed ? _mm256_set1_pd(1.0) : _mm256_setzero_pd();
        setInitialStatus(output);
        output.denominatorInstability.fill(false);
    };
    auto scalarEntireLane = [&](ExpressionProgram::Op operation, Complex inputBase, Complex inputDelta, Complex outputBase, Complex& outputDelta, Complex& derivativeFactor) {
        const Complex inputEndpoint = inputBase + inputDelta;
        const Complex outputEndpoint = ExpressionProgram::evaluateUnary(operation, inputEndpoint);
        outputDelta = nanComplex();
        derivativeFactor = nanComplex();
        ExpressionCenteredStatus status = ExpressionCenteredStatus::Success;
        switch (operation) {
        case ExpressionProgram::Op::Sin:
            outputDelta = centeredSine(inputBase, inputDelta);
            derivativeFactor = std::cos(inputEndpoint);
            break;
        case ExpressionProgram::Op::Cos:
            outputDelta = centeredCosine(inputBase, inputDelta);
            derivativeFactor = -std::sin(inputEndpoint);
            break;
        case ExpressionProgram::Op::Tan:
            status = centeredTangent(inputBase, inputDelta, inputEndpoint, outputBase, outputEndpoint, outputDelta);
            if (status == ExpressionCenteredStatus::Success) derivativeFactor = Complex{1.0, 0.0} + outputEndpoint * outputEndpoint;
            break;
        case ExpressionProgram::Op::Sinh:
            outputDelta = centeredHyperbolicSine(inputBase, inputDelta);
            derivativeFactor = std::cosh(inputEndpoint);
            break;
        case ExpressionProgram::Op::Cosh:
            outputDelta = centeredHyperbolicCosine(inputBase, inputDelta);
            derivativeFactor = std::sinh(inputEndpoint);
            break;
        case ExpressionProgram::Op::Tanh:
            status = centeredHyperbolicTangent(inputBase, inputDelta, inputEndpoint, outputBase, outputEndpoint, outputDelta);
            if (status == ExpressionCenteredStatus::Success) derivativeFactor = Complex{1.0, 0.0} - outputEndpoint * outputEndpoint;
            break;
        case ExpressionProgram::Op::Exp:
            if (!finite(outputBase) || !finite(outputEndpoint))
                status = ExpressionCenteredStatus::NonFinite;
            else {
                outputDelta = relativeExponentialDelta(outputBase, outputEndpoint, inputDelta);
                derivativeFactor = outputEndpoint;
            }
            break;
        case ExpressionProgram::Op::Log:
            status = centeredLogarithm(inputBase, inputDelta, inputEndpoint, outputDelta);
            if (status == ExpressionCenteredStatus::Success) derivativeFactor = Complex{1.0, 0.0} / inputEndpoint;
            break;
        case ExpressionProgram::Op::Log10:
            status = centeredLogarithm(inputBase, inputDelta, inputEndpoint, outputDelta);
            if (status == ExpressionCenteredStatus::Success) {
                outputDelta /= std::log(10.0);
                derivativeFactor = (Complex{1.0, 0.0} / inputEndpoint) / std::log(10.0);
            }
            break;
        case ExpressionProgram::Op::Sqrt:
            status = centeredSquareRoot(inputBase, inputDelta, inputEndpoint, outputBase, outputEndpoint, outputDelta);
            if (status == ExpressionCenteredStatus::Success) derivativeFactor = Complex{0.5, 0.0} / outputEndpoint;
            break;
        case ExpressionProgram::Op::Abs:
            status = centeredAbsolute(inputBase, inputDelta, inputEndpoint, outputDelta);
            if (status == ExpressionCenteredStatus::Success) {
                const double mag = std::abs(inputEndpoint);
                derivativeFactor = mag > 0.0 ? inputEndpoint / mag : Complex{0.0, 0.0};
            }
            break;
        default: return ExpressionCenteredStatus::Unsupported;
        }
        return finish(outputBase, outputDelta, outputEndpoint, status).status;
    };

    for (size_t instructionIndex = 0; instructionIndex < program._code.size(); ++instructionIndex) {
        const ExpressionProgram::Instruction& instruction = program._code[instructionIndex];
        switch (instruction.op) {
        case ExpressionProgram::Op::Constant: pushUnchanged(instruction.value); break;
        case ExpressionProgram::Op::Iteration: pushIteration(); break;
        case ExpressionProgram::Op::Z: pushPerturbed([](const ExpressionContext& reference) { return reference.z; }, [](const ExpressionDeltaContext& delta) { return delta.z; }, true); break;
        case ExpressionProgram::Op::C: pushPerturbed([](const ExpressionContext& reference) { return reference.c; }, [](const ExpressionDeltaContext& delta) { return delta.c; }, false); break;
        case ExpressionProgram::Op::Z0: pushPerturbed([](const ExpressionContext& reference) { return reference.z0; }, [](const ExpressionDeltaContext& delta) { return delta.z0; }, false); break;
        case ExpressionProgram::Op::Parameter: {
            const int parameter = instruction.argument;
            pushPerturbed([parameter](const ExpressionContext& reference) { return reference.parameters[parameter]; }, [parameter](const ExpressionDeltaContext& delta) { return delta.parameters[parameter]; }, false);
            break;
        }
        case ExpressionProgram::Op::Negate: {
            FastLaneBatchValue& value = stack[top - 1];
            value.baseReal = gatherReal(nodeBases, instructionIndex);
            value.baseImaginary = gatherImaginary(nodeBases, instructionIndex);
            value.deltaReal = _mm256_xor_pd(value.deltaReal, signMask);
            value.deltaImaginary = _mm256_xor_pd(value.deltaImaginary, signMask);
            if (derivatives) {
                value.derivativeReal = _mm256_xor_pd(value.derivativeReal, signMask);
                value.derivativeImaginary = _mm256_xor_pd(value.derivativeImaginary, signMask);
                value.derivativeYReal = _mm256_xor_pd(value.derivativeYReal, signMask);
                value.derivativeYImaginary = _mm256_xor_pd(value.derivativeYImaginary, signMask);
            }
            normalize(value);
            break;
        }
        case ExpressionProgram::Op::Square: {
            FastLaneBatchValue input = stack[top - 1];
            FastLaneBatchValue& output = stack[top - 1];
            output.baseReal = gatherReal(nodeBases, instructionIndex);
            output.baseImaginary = gatherImaginary(nodeBases, instructionIndex);
            __m256d firstReal, firstImaginary, squaredReal, squaredImaginary;
            multiply(input.baseReal, input.baseImaginary, input.deltaReal, input.deltaImaginary, firstReal, firstImaginary);
            multiply(input.deltaReal, input.deltaImaginary, input.deltaReal, input.deltaImaginary, squaredReal, squaredImaginary);
            output.deltaReal = _mm256_add_pd(_mm256_add_pd(firstReal, firstReal), squaredReal);
            output.deltaImaginary = _mm256_add_pd(_mm256_add_pd(firstImaginary, firstImaginary), squaredImaginary);
            if (derivatives) {
                const __m256d endpointReal = _mm256_add_pd(input.baseReal, input.deltaReal);
                const __m256d endpointImaginary = _mm256_add_pd(input.baseImaginary, input.deltaImaginary);
                multiply(endpointReal, endpointImaginary, input.derivativeReal, input.derivativeImaginary, output.derivativeReal, output.derivativeImaginary);
                output.derivativeReal = _mm256_add_pd(output.derivativeReal, output.derivativeReal);
                output.derivativeImaginary = _mm256_add_pd(output.derivativeImaginary, output.derivativeImaginary);
                multiply(endpointReal, endpointImaginary, input.derivativeYReal, input.derivativeYImaginary, output.derivativeYReal, output.derivativeYImaginary);
                output.derivativeYReal = _mm256_add_pd(output.derivativeYReal, output.derivativeYReal);
                output.derivativeYImaginary = _mm256_add_pd(output.derivativeYImaginary, output.derivativeYImaginary);
            }
            normalize(output);
            break;
        }
        case ExpressionProgram::Op::Sin:
        case ExpressionProgram::Op::Cos:
        case ExpressionProgram::Op::Tan:
        case ExpressionProgram::Op::Sinh:
        case ExpressionProgram::Op::Cosh:
        case ExpressionProgram::Op::Tanh:
        case ExpressionProgram::Op::Exp:
        case ExpressionProgram::Op::Log:
        case ExpressionProgram::Op::Log10:
        case ExpressionProgram::Op::Sqrt:
        case ExpressionProgram::Op::Abs: {
            FastLaneBatchValue input = stack[top - 1];
            FastLaneBatchValue& output = stack[top - 1];
            output.baseReal = gatherReal(nodeBases, instructionIndex);
            output.baseImaginary = gatherImaginary(nodeBases, instructionIndex);
            FastEntireOperation operation = FastEntireOperation::Sin;
            switch (instruction.op) {
            case ExpressionProgram::Op::Cos: operation = FastEntireOperation::Cos; break;
            case ExpressionProgram::Op::Sinh: operation = FastEntireOperation::Sinh; break;
            case ExpressionProgram::Op::Cosh: operation = FastEntireOperation::Cosh; break;
            case ExpressionProgram::Op::Exp: operation = FastEntireOperation::Exp; break;
            default: break;
            }
            const __m256d companionReal = haveAuxiliaries ? gatherReal(nodeAuxiliaries, instructionIndex) : _mm256_setzero_pd();
            const __m256d companionImaginary = haveAuxiliaries ? gatherImaginary(nodeAuxiliaries, instructionIndex) : _mm256_setzero_pd();
            const uint8_t validMask = successfulMask(input);
            __m256d derivativeFactorReal = _mm256_setzero_pd();
            __m256d derivativeFactorImaginary = _mm256_setzero_pd();
            const bool meromorphicOrBranchOperation = instruction.op == ExpressionProgram::Op::Tan || instruction.op == ExpressionProgram::Op::Tanh || instruction.op == ExpressionProgram::Op::Log || instruction.op == ExpressionProgram::Op::Log10 || instruction.op == ExpressionProgram::Op::Sqrt || instruction.op == ExpressionProgram::Op::Abs;
            const uint8_t fastMask = meromorphicOrBranchOperation ? 0 : fusedCenteredEntireVectors(operation, input.baseReal, input.baseImaginary, output.baseReal, output.baseImaginary, haveAuxiliaries, companionReal, companionImaginary, input.deltaReal, input.deltaImaginary, validMask, output.deltaReal, output.deltaImaginary, derivatives ? &derivativeFactorReal : nullptr, derivatives ? &derivativeFactorImaginary : nullptr);
            output.status = input.status;
            output.denominatorInstability = input.denominatorInstability;
            if (validMask != 0 && fastMask == validMask) {
                if (derivatives) {
                    multiply(derivativeFactorReal, derivativeFactorImaginary, input.derivativeReal, input.derivativeImaginary, output.derivativeReal, output.derivativeImaginary);
                    multiply(derivativeFactorReal, derivativeFactorImaginary, input.derivativeYReal, input.derivativeYImaginary, output.derivativeYReal, output.derivativeYImaginary);
                }
                normalize(output);
                break;
            }
            alignas(32) double inputBaseReal[4], inputBaseImaginary[4], inputReal[4], inputImaginary[4], outputBaseReal[4], outputBaseImaginary[4], outputReal[4], outputImaginary[4], factorReal[4], factorImaginary[4];
            _mm256_store_pd(inputBaseReal, input.baseReal);
            _mm256_store_pd(inputBaseImaginary, input.baseImaginary);
            _mm256_store_pd(inputReal, input.deltaReal);
            _mm256_store_pd(inputImaginary, input.deltaImaginary);
            _mm256_store_pd(outputBaseReal, output.baseReal);
            _mm256_store_pd(outputBaseImaginary, output.baseImaginary);
            _mm256_store_pd(outputReal, output.deltaReal);
            _mm256_store_pd(outputImaginary, output.deltaImaginary);
            if (derivatives) {
                _mm256_store_pd(factorReal, derivativeFactorReal);
                _mm256_store_pd(factorImaginary, derivativeFactorImaginary);
            }
            for (int lane = 0; lane < 4; ++lane) {
                if (input.status[lane] != ExpressionCenteredStatus::Success) continue;
                if ((fastMask & (uint8_t{1} << lane)) != 0) continue;
                Complex outputDelta;
                Complex derivativeFactor;
                output.status[lane] = scalarEntireLane(instruction.op, {inputBaseReal[lane], inputBaseImaginary[lane]}, {inputReal[lane], inputImaginary[lane]}, {outputBaseReal[lane], outputBaseImaginary[lane]}, outputDelta, derivativeFactor);
                if (meromorphicOrBranchOperation && output.status[lane] != ExpressionCenteredStatus::Success) output.denominatorInstability[lane] = true;
                outputReal[lane] = outputDelta.real();
                outputImaginary[lane] = outputDelta.imag();
                if (derivatives) {
                    factorReal[lane] = derivativeFactor.real();
                    factorImaginary[lane] = derivativeFactor.imag();
                }
            }
            output.deltaReal = _mm256_load_pd(outputReal);
            output.deltaImaginary = _mm256_load_pd(outputImaginary);
            if (derivatives) {
                derivativeFactorReal = _mm256_load_pd(factorReal);
                derivativeFactorImaginary = _mm256_load_pd(factorImaginary);
                multiply(derivativeFactorReal, derivativeFactorImaginary, input.derivativeReal, input.derivativeImaginary, output.derivativeReal, output.derivativeImaginary);
                multiply(derivativeFactorReal, derivativeFactorImaginary, input.derivativeYReal, input.derivativeYImaginary, output.derivativeYReal, output.derivativeYImaginary);
            }
            normalize(output);
            break;
        }
        case ExpressionProgram::Op::Norm:
        case ExpressionProgram::Op::Conjugate:
        case ExpressionProgram::Op::Real:
        case ExpressionProgram::Op::Imaginary: {
            FastLaneBatchValue input = stack[top - 1];
            FastLaneBatchValue& output = stack[top - 1];
            output.baseReal = gatherReal(nodeBases, instructionIndex);
            output.baseImaginary = gatherImaginary(nodeBases, instructionIndex);
            if (instruction.op == ExpressionProgram::Op::Conjugate) {
                output.deltaReal = input.deltaReal;
                output.deltaImaginary = _mm256_xor_pd(input.deltaImaginary, signMask);
                if (derivatives) {
                    output.derivativeReal = input.derivativeReal;
                    output.derivativeImaginary = _mm256_xor_pd(input.derivativeImaginary, signMask);
                    output.derivativeYReal = input.derivativeYReal;
                    output.derivativeYImaginary = _mm256_xor_pd(input.derivativeYImaginary, signMask);
                }
            } else if (instruction.op == ExpressionProgram::Op::Real) {
                output.deltaReal = input.deltaReal;
                output.deltaImaginary = _mm256_setzero_pd();
                if (derivatives) {
                    output.derivativeReal = input.derivativeReal;
                    output.derivativeImaginary = _mm256_setzero_pd();
                    output.derivativeYReal = input.derivativeYReal;
                    output.derivativeYImaginary = _mm256_setzero_pd();
                }
            } else if (instruction.op == ExpressionProgram::Op::Imaginary) {
                output.deltaReal = input.deltaImaginary;
                output.deltaImaginary = _mm256_setzero_pd();
                if (derivatives) {
                    output.derivativeReal = input.derivativeImaginary;
                    output.derivativeImaginary = _mm256_setzero_pd();
                    output.derivativeYReal = input.derivativeYImaginary;
                    output.derivativeYImaginary = _mm256_setzero_pd();
                }
            } else {
                const __m256d endpointReal = _mm256_add_pd(input.baseReal, input.deltaReal);
                const __m256d endpointImaginary = _mm256_add_pd(input.baseImaginary, input.deltaImaginary);
                const __m256d cross = _mm256_add_pd(_mm256_mul_pd(input.baseReal, input.deltaReal), _mm256_mul_pd(input.baseImaginary, input.deltaImaginary));
                const __m256d deltaNorm = _mm256_add_pd(_mm256_mul_pd(input.deltaReal, input.deltaReal), _mm256_mul_pd(input.deltaImaginary, input.deltaImaginary));
                output.deltaReal = _mm256_add_pd(_mm256_add_pd(cross, cross), deltaNorm);
                output.deltaImaginary = _mm256_setzero_pd();
                if (derivatives) {
                    const __m256d twiceEndpointReal = _mm256_add_pd(endpointReal, endpointReal);
                    const __m256d twiceEndpointImaginary = _mm256_add_pd(endpointImaginary, endpointImaginary);
                    output.derivativeReal = _mm256_add_pd(_mm256_mul_pd(twiceEndpointReal, input.derivativeReal), _mm256_mul_pd(twiceEndpointImaginary, input.derivativeImaginary));
                    output.derivativeImaginary = _mm256_setzero_pd();
                    output.derivativeYReal = _mm256_add_pd(_mm256_mul_pd(twiceEndpointReal, input.derivativeYReal), _mm256_mul_pd(twiceEndpointImaginary, input.derivativeYImaginary));
                    output.derivativeYImaginary = _mm256_setzero_pd();
                }
            }
            output.status = input.status;
            normalize(output);
            break;
        }
        case ExpressionProgram::Op::Add:
        case ExpressionProgram::Op::Subtract:
        case ExpressionProgram::Op::Multiply:
        case ExpressionProgram::Op::Divide:
        case ExpressionProgram::Op::Power:
        case ExpressionProgram::Op::MakeComplex: {
            FastLaneBatchValue right = stack[--top];
            FastLaneBatchValue left = stack[top - 1];
            FastLaneBatchValue& output = stack[top - 1];
            output.baseReal = gatherReal(nodeBases, instructionIndex);
            output.baseImaginary = gatherImaginary(nodeBases, instructionIndex);
            for (int lane = 0; lane < 4; ++lane) {
                output.status[lane] = combine(left.status[lane], right.status[lane]);
                output.denominatorInstability[lane] = left.denominatorInstability[lane] || right.denominatorInstability[lane];
            }
            if (instruction.op == ExpressionProgram::Op::MakeComplex) {
                output.deltaReal = left.deltaReal;
                output.deltaImaginary = right.deltaReal;
                if (derivatives) {
                    output.derivativeReal = left.derivativeReal;
                    output.derivativeImaginary = right.derivativeReal;
                    output.derivativeYReal = left.derivativeYReal;
                    output.derivativeYImaginary = right.derivativeYReal;
                }
            } else if (instruction.op == ExpressionProgram::Op::Add) {
                output.deltaReal = _mm256_add_pd(left.deltaReal, right.deltaReal);
                output.deltaImaginary = _mm256_add_pd(left.deltaImaginary, right.deltaImaginary);
                if (derivatives) {
                    output.derivativeReal = _mm256_add_pd(left.derivativeReal, right.derivativeReal);
                    output.derivativeImaginary = _mm256_add_pd(left.derivativeImaginary, right.derivativeImaginary);
                    output.derivativeYReal = _mm256_add_pd(left.derivativeYReal, right.derivativeYReal);
                    output.derivativeYImaginary = _mm256_add_pd(left.derivativeYImaginary, right.derivativeYImaginary);
                }
            } else if (instruction.op == ExpressionProgram::Op::Subtract) {
                output.deltaReal = _mm256_sub_pd(left.deltaReal, right.deltaReal);
                output.deltaImaginary = _mm256_sub_pd(left.deltaImaginary, right.deltaImaginary);
                if (derivatives) {
                    output.derivativeReal = _mm256_sub_pd(left.derivativeReal, right.derivativeReal);
                    output.derivativeImaginary = _mm256_sub_pd(left.derivativeImaginary, right.derivativeImaginary);
                    output.derivativeYReal = _mm256_sub_pd(left.derivativeYReal, right.derivativeYReal);
                    output.derivativeYImaginary = _mm256_sub_pd(left.derivativeYImaginary, right.derivativeYImaginary);
                }
            } else if (instruction.op == ExpressionProgram::Op::Divide) {
                alignas(32) double leftBaseReal[4], leftBaseImaginary[4], rightBaseReal[4], rightBaseImaginary[4], quotientBaseReal[4], quotientBaseImaginary[4];
                alignas(32) double leftDeltaReal[4], leftDeltaImaginary[4], rightDeltaReal[4], rightDeltaImaginary[4], outputReal[4], outputImaginary[4];
                alignas(32) double leftDerivativeReal[4], leftDerivativeImaginary[4], leftDerivativeYReal[4], leftDerivativeYImaginary[4];
                alignas(32) double rightDerivativeReal[4], rightDerivativeImaginary[4], rightDerivativeYReal[4], rightDerivativeYImaginary[4];
                alignas(32) double outputDerivativeReal[4]{}, outputDerivativeImaginary[4]{}, outputDerivativeYReal[4]{}, outputDerivativeYImaginary[4]{};
                _mm256_store_pd(leftBaseReal, left.baseReal);
                _mm256_store_pd(leftBaseImaginary, left.baseImaginary);
                _mm256_store_pd(rightBaseReal, right.baseReal);
                _mm256_store_pd(rightBaseImaginary, right.baseImaginary);
                _mm256_store_pd(quotientBaseReal, output.baseReal);
                _mm256_store_pd(quotientBaseImaginary, output.baseImaginary);
                _mm256_store_pd(leftDeltaReal, left.deltaReal);
                _mm256_store_pd(leftDeltaImaginary, left.deltaImaginary);
                _mm256_store_pd(rightDeltaReal, right.deltaReal);
                _mm256_store_pd(rightDeltaImaginary, right.deltaImaginary);
                if (derivatives) {
                    _mm256_store_pd(leftDerivativeReal, left.derivativeReal);
                    _mm256_store_pd(leftDerivativeImaginary, left.derivativeImaginary);
                    _mm256_store_pd(leftDerivativeYReal, left.derivativeYReal);
                    _mm256_store_pd(leftDerivativeYImaginary, left.derivativeYImaginary);
                    _mm256_store_pd(rightDerivativeReal, right.derivativeReal);
                    _mm256_store_pd(rightDerivativeImaginary, right.derivativeImaginary);
                    _mm256_store_pd(rightDerivativeYReal, right.derivativeYReal);
                    _mm256_store_pd(rightDerivativeYImaginary, right.derivativeYImaginary);
                }
                for (int lane = 0; lane < 4; ++lane) {
                    Complex delta = nanComplex();
                    if (output.status[lane] == ExpressionCenteredStatus::Success) {
                        const Complex leftBase{leftBaseReal[lane], leftBaseImaginary[lane]};
                        const Complex rightBase{rightBaseReal[lane], rightBaseImaginary[lane]};
                        const Complex leftDelta{leftDeltaReal[lane], leftDeltaImaginary[lane]};
                        const Complex rightDelta{rightDeltaReal[lane], rightDeltaImaginary[lane]};
                        const Complex leftEndpoint = leftBase + leftDelta;
                        const Complex rightEndpoint = rightBase + rightDelta;
                        output.status[lane] = centeredDivisionDelta(leftDelta, rightBase, rightDelta, rightEndpoint, {quotientBaseReal[lane], quotientBaseImaginary[lane]}, delta, true);
                        if (output.status[lane] != ExpressionCenteredStatus::Success) output.denominatorInstability[lane] = true;
                        if (derivatives && output.status[lane] == ExpressionCenteredStatus::Success) {
                            const Complex denominatorSquared = rightEndpoint * rightEndpoint;
                            const Complex derivative = (Complex{leftDerivativeReal[lane], leftDerivativeImaginary[lane]} * rightEndpoint - leftEndpoint * Complex{rightDerivativeReal[lane], rightDerivativeImaginary[lane]}) / denominatorSquared;
                            const Complex derivativeY = (Complex{leftDerivativeYReal[lane], leftDerivativeYImaginary[lane]} * rightEndpoint - leftEndpoint * Complex{rightDerivativeYReal[lane], rightDerivativeYImaginary[lane]}) / denominatorSquared;
                            if (!finite(derivative) || !finite(derivativeY)) {
                                output.status[lane] = ExpressionCenteredStatus::NonFinite;
                                output.denominatorInstability[lane] = true;
                            } else {
                                outputDerivativeReal[lane] = derivative.real();
                                outputDerivativeImaginary[lane] = derivative.imag();
                                outputDerivativeYReal[lane] = derivativeY.real();
                                outputDerivativeYImaginary[lane] = derivativeY.imag();
                            }
                        }
                    }
                    outputReal[lane] = delta.real();
                    outputImaginary[lane] = delta.imag();
                }
                output.deltaReal = _mm256_load_pd(outputReal);
                output.deltaImaginary = _mm256_load_pd(outputImaginary);
                if (derivatives) {
                    output.derivativeReal = _mm256_load_pd(outputDerivativeReal);
                    output.derivativeImaginary = _mm256_load_pd(outputDerivativeImaginary);
                    output.derivativeYReal = _mm256_load_pd(outputDerivativeYReal);
                    output.derivativeYImaginary = _mm256_load_pd(outputDerivativeYImaginary);
                }
            } else if (instruction.op == ExpressionProgram::Op::Power) {
                alignas(32) double leftBaseReal[4], leftBaseImaginary[4], rightBaseReal[4], rightBaseImaginary[4], outputBaseReal[4], outputBaseImaginary[4];
                alignas(32) double leftDeltaReal[4], leftDeltaImaginary[4], rightDeltaReal[4], rightDeltaImaginary[4], outputReal[4], outputImaginary[4];
                alignas(32) double leftDerivativeReal[4], leftDerivativeImaginary[4], leftDerivativeYReal[4], leftDerivativeYImaginary[4];
                alignas(32) double rightDerivativeReal[4], rightDerivativeImaginary[4], rightDerivativeYReal[4], rightDerivativeYImaginary[4];
                alignas(32) double outputDerivativeReal[4]{}, outputDerivativeImaginary[4]{}, outputDerivativeYReal[4]{}, outputDerivativeYImaginary[4]{};
                _mm256_store_pd(leftBaseReal, left.baseReal);
                _mm256_store_pd(leftBaseImaginary, left.baseImaginary);
                _mm256_store_pd(rightBaseReal, right.baseReal);
                _mm256_store_pd(rightBaseImaginary, right.baseImaginary);
                _mm256_store_pd(outputBaseReal, output.baseReal);
                _mm256_store_pd(outputBaseImaginary, output.baseImaginary);
                _mm256_store_pd(leftDeltaReal, left.deltaReal);
                _mm256_store_pd(leftDeltaImaginary, left.deltaImaginary);
                _mm256_store_pd(rightDeltaReal, right.deltaReal);
                _mm256_store_pd(rightDeltaImaginary, right.deltaImaginary);
                if (derivatives) {
                    _mm256_store_pd(leftDerivativeReal, left.derivativeReal);
                    _mm256_store_pd(leftDerivativeImaginary, left.derivativeImaginary);
                    _mm256_store_pd(leftDerivativeYReal, left.derivativeYReal);
                    _mm256_store_pd(leftDerivativeYImaginary, left.derivativeYImaginary);
                    _mm256_store_pd(rightDerivativeReal, right.derivativeReal);
                    _mm256_store_pd(rightDerivativeImaginary, right.derivativeImaginary);
                    _mm256_store_pd(rightDerivativeYReal, right.derivativeYReal);
                    _mm256_store_pd(rightDerivativeYImaginary, right.derivativeYImaginary);
                }
                for (int lane = 0; lane < 4; ++lane) {
                    Complex delta = nanComplex();
                    if (output.status[lane] == ExpressionCenteredStatus::Success) {
                        const Complex leftBase{leftBaseReal[lane], leftBaseImaginary[lane]};
                        const Complex rightBase{rightBaseReal[lane], rightBaseImaginary[lane]};
                        const Complex outputBase{outputBaseReal[lane], outputBaseImaginary[lane]};
                        const Complex leftDelta{leftDeltaReal[lane], leftDeltaImaginary[lane]};
                        const Complex rightDelta{rightDeltaReal[lane], rightDeltaImaginary[lane]};
                        const Complex leftEndpoint = leftBase + leftDelta;
                        const Complex rightEndpoint = rightBase + rightDelta;
                        const Complex outputEndpoint = ExpressionProgram::evaluateBinary(ExpressionProgram::Op::Power, leftEndpoint, rightEndpoint);
                        output.status[lane] = centeredPower(leftBase, leftDelta, leftEndpoint, rightBase, rightDelta, rightEndpoint, outputBase, outputEndpoint, delta);
                        if (derivatives && output.status[lane] == ExpressionCenteredStatus::Success) {
                            if (zero(leftEndpoint)) {
                                outputDerivativeReal[lane] = 0.0;
                                outputDerivativeImaginary[lane] = 0.0;
                                outputDerivativeYReal[lane] = 0.0;
                                outputDerivativeYImaginary[lane] = 0.0;
                            } else {
                                const Complex termU = (rightEndpoint / leftEndpoint) * Complex{leftDerivativeReal[lane], leftDerivativeImaginary[lane]};
                                const Complex termV = std::log(leftEndpoint) * Complex{rightDerivativeReal[lane], rightDerivativeImaginary[lane]};
                                const Complex derivative = outputEndpoint * (termU + termV);

                                const Complex termUY = (rightEndpoint / leftEndpoint) * Complex{leftDerivativeYReal[lane], leftDerivativeYImaginary[lane]};
                                const Complex termVY = std::log(leftEndpoint) * Complex{rightDerivativeYReal[lane], rightDerivativeYImaginary[lane]};
                                const Complex derivativeY = outputEndpoint * (termUY + termVY);

                                if (!finite(derivative) || !finite(derivativeY)) {
                                    output.status[lane] = ExpressionCenteredStatus::NonFinite;
                                } else {
                                    outputDerivativeReal[lane] = derivative.real();
                                    outputDerivativeImaginary[lane] = derivative.imag();
                                    outputDerivativeYReal[lane] = derivativeY.real();
                                    outputDerivativeYImaginary[lane] = derivativeY.imag();
                                }
                            }
                        }
                    }
                    outputReal[lane] = delta.real();
                    outputImaginary[lane] = delta.imag();
                }
                output.deltaReal = _mm256_load_pd(outputReal);
                output.deltaImaginary = _mm256_load_pd(outputImaginary);
                if (derivatives) {
                    output.derivativeReal = _mm256_load_pd(outputDerivativeReal);
                    output.derivativeImaginary = _mm256_load_pd(outputDerivativeImaginary);
                    output.derivativeYReal = _mm256_load_pd(outputDerivativeYReal);
                    output.derivativeYImaginary = _mm256_load_pd(outputDerivativeYImaginary);
                }
            } else {
                __m256d firstReal, firstImaginary, secondReal, secondImaginary, thirdReal, thirdImaginary;
                multiply(left.baseReal, left.baseImaginary, right.deltaReal, right.deltaImaginary, firstReal, firstImaginary);
                multiply(right.baseReal, right.baseImaginary, left.deltaReal, left.deltaImaginary, secondReal, secondImaginary);
                multiply(left.deltaReal, left.deltaImaginary, right.deltaReal, right.deltaImaginary, thirdReal, thirdImaginary);
                output.deltaReal = _mm256_add_pd(_mm256_add_pd(firstReal, secondReal), thirdReal);
                output.deltaImaginary = _mm256_add_pd(_mm256_add_pd(firstImaginary, secondImaginary), thirdImaginary);
                if (derivatives) {
                    const __m256d leftEndpointReal = _mm256_add_pd(left.baseReal, left.deltaReal);
                    const __m256d leftEndpointImaginary = _mm256_add_pd(left.baseImaginary, left.deltaImaginary);
                    const __m256d rightEndpointReal = _mm256_add_pd(right.baseReal, right.deltaReal);
                    const __m256d rightEndpointImaginary = _mm256_add_pd(right.baseImaginary, right.deltaImaginary);
                    multiply(left.derivativeReal, left.derivativeImaginary, rightEndpointReal, rightEndpointImaginary, firstReal, firstImaginary);
                    multiply(leftEndpointReal, leftEndpointImaginary, right.derivativeReal, right.derivativeImaginary, secondReal, secondImaginary);
                    output.derivativeReal = _mm256_add_pd(firstReal, secondReal);
                    output.derivativeImaginary = _mm256_add_pd(firstImaginary, secondImaginary);
                    multiply(left.derivativeYReal, left.derivativeYImaginary, rightEndpointReal, rightEndpointImaginary, firstReal, firstImaginary);
                    multiply(leftEndpointReal, leftEndpointImaginary, right.derivativeYReal, right.derivativeYImaginary, secondReal, secondImaginary);
                    output.derivativeYReal = _mm256_add_pd(firstReal, secondReal);
                    output.derivativeYImaginary = _mm256_add_pd(firstImaginary, secondImaginary);
                }
            }
            normalize(output);
            break;
        }
        default: return false;
        }
    }
    if (top != 1) return false;
    normalize(stack[0]);
    alignas(32) double baseReal[4], baseImaginary[4], resultReal[4], resultImaginary[4], derivativeReal[4], derivativeImaginary[4], derivativeYReal[4], derivativeYImaginary[4];
    _mm256_store_pd(baseReal, stack[0].baseReal);
    _mm256_store_pd(baseImaginary, stack[0].baseImaginary);
    _mm256_store_pd(resultReal, stack[0].deltaReal);
    _mm256_store_pd(resultImaginary, stack[0].deltaImaginary);
    if (derivatives) {
        _mm256_store_pd(derivativeReal, stack[0].derivativeReal);
        _mm256_store_pd(derivativeImaginary, stack[0].derivativeImaginary);
        _mm256_store_pd(derivativeYReal, stack[0].derivativeYReal);
        _mm256_store_pd(derivativeYImaginary, stack[0].derivativeYImaginary);
    }
    for (int lane = 0; lane < 4; ++lane) {
        if ((activeMask & (uint8_t{1} << lane)) == 0) {
            results[lane] = {nanComplex(), nanComplex(), ExpressionCenteredStatus::InvalidProgram};
            if (derivatives) derivatives[lane] = std::numeric_limits<double>::quiet_NaN();
            continue;
        }
        results[lane] = {{baseReal[lane], baseImaginary[lane]}, {resultReal[lane], resultImaginary[lane]}, stack[0].status[lane], stack[0].denominatorInstability[lane]};
        if (derivatives) {
            const double a = derivativeReal[lane];
            const double b = derivativeYReal[lane];
            const double c = derivativeImaginary[lane];
            const double d = derivativeYImaginary[lane];
            const double scale = std::max({std::abs(a), std::abs(b), std::abs(c), std::abs(d)});
            double norm = 0.0;
            if (scale != 0.0 && std::isfinite(scale)) {
                const double as = a / scale;
                const double bs = b / scale;
                const double cs = c / scale;
                const double ds = d / scale;
                const double sum = as * as + bs * bs + cs * cs + ds * ds;
                const double determinant = as * ds - bs * cs;
                const double discriminant = std::max(0.0, sum * sum - 4.0 * determinant * determinant);
                norm = scale * std::sqrt(0.5 * (sum + std::sqrt(discriminant)));
            } else if (!std::isfinite(scale)) {
                norm = std::numeric_limits<double>::quiet_NaN();
            }
            derivatives[lane] = results[lane].success() && std::isfinite(norm) ? norm : std::numeric_limits<double>::quiet_NaN();
        }
    }
    return true;
}

bool ExpressionCenteredEvaluator::evaluate4WithNodeBasesImpl(const ExpressionProgram& program, const ExpressionContext& reference, const Complex* nodeBases, const ExpressionDeltaContext* inputDeltas, ExpressionCenteredResult* results) {
    if (!program._valid || !inputDeltas || !results) return false;
    struct BatchValue {
        Complex base{};
        std::array<Complex, 4> delta{};
        std::array<Complex, 4> endpoint{};
        std::array<ExpressionCenteredStatus, 4> status{};
        std::array<bool, 4> denominatorInstability{};
    };
    std::array<BatchValue, ExpressionProgram::MAX_STACK> stack;
    size_t top = 0;
    auto pushUnchanged = [&](Complex value) {
        BatchValue& output = stack[top++];
        output.base = value;
        for (int lane = 0; lane < 4; ++lane) {
            output.delta[lane] = {};
            output.endpoint[lane] = value;
            output.status[lane] = finite(value) ? ExpressionCenteredStatus::Success : ExpressionCenteredStatus::NonFinite;
            output.denominatorInstability[lane] = false;
        }
    };
    auto pushPerturbed = [&](Complex base, auto deltaGetter) {
        BatchValue& output = stack[top++];
        output.base = base;
        for (int lane = 0; lane < 4; ++lane) {
            output.delta[lane] = deltaGetter(inputDeltas[lane]);
            output.endpoint[lane] = base + output.delta[lane];
            output.status[lane] = finite(base) && finite(output.delta[lane]) && finite(output.endpoint[lane]) ? ExpressionCenteredStatus::Success : ExpressionCenteredStatus::NonFinite;
            output.denominatorInstability[lane] = false;
        }
    };
    auto unary = [&](const ExpressionProgram::Instruction& instruction, Complex sharedBase) {
        BatchValue input = stack[top - 1];
        BatchValue& output = stack[top - 1];
        output.base = sharedBase;
        for (int lane = 0; lane < 4; ++lane) {
            output.endpoint[lane] = ExpressionProgram::evaluateUnary(instruction.op, input.endpoint[lane]);
            ExpressionCenteredStatus status = input.status[lane];
            bool denominatorInstability = input.denominatorInstability[lane];
            Complex delta;
            if (status == ExpressionCenteredStatus::Success) {
                switch (instruction.op) {
                case ExpressionProgram::Op::Negate: delta = -input.delta[lane]; break;
                case ExpressionProgram::Op::Square: delta = centeredSquare(input.base, input.delta[lane]); break;
                case ExpressionProgram::Op::Sin: delta = centeredSine(input.base, input.delta[lane]); break;
                case ExpressionProgram::Op::Cos: delta = centeredCosine(input.base, input.delta[lane]); break;
                case ExpressionProgram::Op::Tan:
                    status = centeredTangent(input.base, input.delta[lane], input.endpoint[lane], output.base, output.endpoint[lane], delta);
                    denominatorInstability = denominatorInstability || status != ExpressionCenteredStatus::Success;
                    break;
                case ExpressionProgram::Op::Sinh: delta = centeredHyperbolicSine(input.base, input.delta[lane]); break;
                case ExpressionProgram::Op::Cosh: delta = centeredHyperbolicCosine(input.base, input.delta[lane]); break;
                case ExpressionProgram::Op::Tanh:
                    status = centeredHyperbolicTangent(input.base, input.delta[lane], input.endpoint[lane], output.base, output.endpoint[lane], delta);
                    denominatorInstability = denominatorInstability || status != ExpressionCenteredStatus::Success;
                    break;
                case ExpressionProgram::Op::Exp:
                    if (!finite(output.base) || !finite(output.endpoint[lane]))
                        status = ExpressionCenteredStatus::NonFinite;
                    else
                        delta = relativeExponentialDelta(output.base, output.endpoint[lane], input.delta[lane]);
                    break;
                case ExpressionProgram::Op::Log: status = centeredLogarithm(input.base, input.delta[lane], input.endpoint[lane], delta); break;
                case ExpressionProgram::Op::Log10:
                    status = centeredLogarithm(input.base, input.delta[lane], input.endpoint[lane], delta);
                    if (status == ExpressionCenteredStatus::Success) delta /= std::log(10.0);
                    break;
                case ExpressionProgram::Op::Sqrt: status = centeredSquareRoot(input.base, input.delta[lane], input.endpoint[lane], output.base, output.endpoint[lane], delta); break;
                case ExpressionProgram::Op::Abs: status = centeredAbsolute(input.base, input.delta[lane], input.endpoint[lane], delta); break;
                case ExpressionProgram::Op::Norm: delta = centeredNormDelta(input.base, input.delta[lane]); break;
                case ExpressionProgram::Op::Arg:
                    status = centeredLogarithm(input.base, input.delta[lane], input.endpoint[lane], delta);
                    if (status == ExpressionCenteredStatus::Success) delta = {delta.imag(), 0.0};
                    break;
                case ExpressionProgram::Op::Conjugate: delta = std::conj(input.delta[lane]); break;
                case ExpressionProgram::Op::Real: delta = {input.delta[lane].real(), 0.0}; break;
                case ExpressionProgram::Op::Imaginary: delta = {input.delta[lane].imag(), 0.0}; break;
                default:
                    status = ExpressionCenteredStatus::Unsupported;
                    delta = nanComplex();
                    break;
                }
            } else {
                delta = nanComplex();
            }
            CenteredValue finished = finish(output.base, delta, output.endpoint[lane], status, denominatorInstability);
            output.delta[lane] = finished.delta;
            output.endpoint[lane] = finished.endpoint;
            output.status[lane] = finished.status;
            output.denominatorInstability[lane] = finished.denominatorInstability;
        }
    };
    auto binary = [&](const ExpressionProgram::Instruction& instruction, Complex sharedBase) {
        BatchValue right = stack[--top];
        BatchValue left = stack[top - 1];
        BatchValue& output = stack[top - 1];
        output.base = sharedBase;
        for (int lane = 0; lane < 4; ++lane) {
            output.endpoint[lane] = ExpressionProgram::evaluateBinary(instruction.op, left.endpoint[lane], right.endpoint[lane]);
            ExpressionCenteredStatus status = combine(left.status[lane], right.status[lane]);
            bool denominatorInstability = left.denominatorInstability[lane] || right.denominatorInstability[lane];
            Complex delta;
            if (status == ExpressionCenteredStatus::Success) {
                switch (instruction.op) {
                case ExpressionProgram::Op::Add: delta = left.delta[lane] + right.delta[lane]; break;
                case ExpressionProgram::Op::Subtract: delta = left.delta[lane] - right.delta[lane]; break;
                case ExpressionProgram::Op::Multiply: delta = centeredMultiply(left.base, left.delta[lane], right.base, right.delta[lane]); break;
                case ExpressionProgram::Op::Divide:
                    status = centeredDivisionDelta(left.delta[lane], right.base, right.delta[lane], right.endpoint[lane], output.base, delta, true);
                    denominatorInstability = denominatorInstability || status != ExpressionCenteredStatus::Success;
                    break;
                case ExpressionProgram::Op::Power: status = centeredPower(left.base, left.delta[lane], left.endpoint[lane], right.base, right.delta[lane], right.endpoint[lane], output.base, output.endpoint[lane], delta); break;
                case ExpressionProgram::Op::MakeComplex: delta = {left.delta[lane].real(), right.delta[lane].real()}; break;
                case ExpressionProgram::Op::Polar: {
                    CenteredValue leftValue = finish(left.base, left.delta[lane], left.endpoint[lane], left.status[lane]);
                    CenteredValue rightValue = finish(right.base, right.delta[lane], right.endpoint[lane], right.status[lane]);
                    status = centeredPolar(leftValue, rightValue, output.base, output.endpoint[lane], delta);
                    break;
                }
                default:
                    status = ExpressionCenteredStatus::Unsupported;
                    delta = nanComplex();
                    break;
                }
            } else {
                delta = nanComplex();
            }
            CenteredValue finished = finish(output.base, delta, output.endpoint[lane], status, denominatorInstability);
            output.delta[lane] = finished.delta;
            output.endpoint[lane] = finished.endpoint;
            output.status[lane] = finished.status;
            output.denominatorInstability[lane] = finished.denominatorInstability;
        }
    };

    for (size_t instructionIndex = 0; instructionIndex < program._code.size(); ++instructionIndex) {
        const ExpressionProgram::Instruction& instruction = program._code[instructionIndex];
        const Complex sharedBase = nodeBases ? nodeBases[instructionIndex] : Complex{};
        switch (instruction.op) {
        case ExpressionProgram::Op::Constant: pushUnchanged(nodeBases ? sharedBase : instruction.value); break;
        case ExpressionProgram::Op::Z: pushPerturbed(nodeBases ? sharedBase : reference.z, [](const ExpressionDeltaContext& delta) { return delta.z; }); break;
        case ExpressionProgram::Op::C: pushPerturbed(nodeBases ? sharedBase : reference.c, [](const ExpressionDeltaContext& delta) { return delta.c; }); break;
        case ExpressionProgram::Op::Z0: pushPerturbed(nodeBases ? sharedBase : reference.z0, [](const ExpressionDeltaContext& delta) { return delta.z0; }); break;
        case ExpressionProgram::Op::Iteration: pushUnchanged(nodeBases ? sharedBase : Complex{(double)reference.iteration, 0.0}); break;
        case ExpressionProgram::Op::Parameter: {
            const int parameter = instruction.argument;
            pushPerturbed(nodeBases ? sharedBase : reference.parameters[parameter], [parameter](const ExpressionDeltaContext& delta) { return delta.parameters[parameter]; });
            break;
        }
        case ExpressionProgram::Op::OrbitInvariant: {
            BatchValue& output = stack[top++];
            output.base = nodeBases ? sharedBase : nanComplex();
            for (int lane = 0; lane < 4; ++lane) {
                output.delta[lane] = nanComplex();
                output.endpoint[lane] = output.base;
                output.status[lane] = ExpressionCenteredStatus::Unsupported;
                output.denominatorInstability[lane] = false;
            }
            break;
        }
        case ExpressionProgram::Op::Negate:
        case ExpressionProgram::Op::Square:
        case ExpressionProgram::Op::Sin:
        case ExpressionProgram::Op::Cos:
        case ExpressionProgram::Op::Tan:
        case ExpressionProgram::Op::Sinh:
        case ExpressionProgram::Op::Cosh:
        case ExpressionProgram::Op::Tanh:
        case ExpressionProgram::Op::Exp:
        case ExpressionProgram::Op::Log:
        case ExpressionProgram::Op::Log10:
        case ExpressionProgram::Op::Sqrt:
        case ExpressionProgram::Op::Abs:
        case ExpressionProgram::Op::Norm:
        case ExpressionProgram::Op::Arg:
        case ExpressionProgram::Op::Conjugate:
        case ExpressionProgram::Op::Real:
        case ExpressionProgram::Op::Imaginary: unary(instruction, nodeBases ? sharedBase : ExpressionProgram::evaluateUnary(instruction.op, stack[top - 1].base)); break;
        case ExpressionProgram::Op::Add:
        case ExpressionProgram::Op::Subtract:
        case ExpressionProgram::Op::Multiply:
        case ExpressionProgram::Op::Divide:
        case ExpressionProgram::Op::Power:
        case ExpressionProgram::Op::MakeComplex:
        case ExpressionProgram::Op::Polar: binary(instruction, nodeBases ? sharedBase : ExpressionProgram::evaluateBinary(instruction.op, stack[top - 2].base, stack[top - 1].base)); break;
        }
    }
    if (top != 1) return false;
    for (int lane = 0; lane < 4; ++lane) results[lane] = {stack[0].base, stack[0].delta[lane], stack[0].status[lane], stack[0].denominatorInstability[lane]};
    return true;
}

} // namespace formula

#ifdef _MSC_VER
#pragma float_control(pop)
#endif
