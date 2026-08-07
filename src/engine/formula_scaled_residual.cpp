#include "formula_scaled_residual.h"

#include <algorithm>
#include <cfloat>
#include <cstring>
#include <limits>

namespace formula {

namespace {

bool checkedAddExponent(
        int64_t left, int64_t right, int64_t& output) {
    if ((right > 0 &&
         left > std::numeric_limits<int64_t>::max() - right) ||
        (right < 0 &&
         left < std::numeric_limits<int64_t>::min() - right))
        return false;
    output = left + right;
    return true;
}

bool checkedAddThreeExponents(
        int64_t first, int64_t second, int64_t third,
        int64_t& output) {
    int64_t partial = 0;
    if (checkedAddExponent(first, third, partial) &&
        checkedAddExponent(partial, second, output))
        return true;
    if (checkedAddExponent(second, third, partial) &&
        checkedAddExponent(partial, first, output))
        return true;
    return checkedAddExponent(first, second, partial) &&
           checkedAddExponent(partial, third, output);
}

ScaledArithmeticStatus normalize(
        double mantissa, int64_t exponent,
        ScaledRealValue& output) {
    ScaledRealValue result{};
    if (!std::isfinite(mantissa)) {
        result.mantissa = mantissa;
        output = result;
        return ScaledArithmeticStatus::Nonfinite;
    }
    if (mantissa == 0.0) {
        result.mantissa = mantissa;
        output = result;
        return ScaledArithmeticStatus::Success;
    }
    int adjustment = 0;
    result.mantissa = std::frexp(mantissa, &adjustment);
    if (!checkedAddExponent(
            exponent, static_cast<int64_t>(adjustment),
            result.exponent))
        return ScaledArithmeticStatus::ExponentRange;
    output = result;
    return ScaledArithmeticStatus::Success;
}

ScaledArithmeticStatus validate(
        const ScaledRealValue& value) {
    if (!value.isFinite())
        return ScaledArithmeticStatus::Nonfinite;
    return value.isNormalized()
        ? ScaledArithmeticStatus::Success
        : ScaledArithmeticStatus::ExponentRange;
}

ScaledArithmeticStatus validate(
        const ScaledComplexValue& value) {
    ScaledArithmeticStatus status = validate(value.re);
    if (status != ScaledArithmeticStatus::Success)
        return status;
    return validate(value.im);
}

double signedZero(bool negative) {
    return negative ? -0.0 : 0.0;
}

ScaledRealValue absoluteValue(ScaledRealValue value) {
    value.mantissa = std::abs(value.mantissa);
    return value;
}

ScaledArithmeticStatus maxComponentMagnitude(
        const ScaledComplexValue& value,
        ScaledRealValue& output) {
    ScaledArithmeticStatus status = validate(value);
    if (status != ScaledArithmeticStatus::Success)
        return status;
    ScaledRealValue real = absoluteValue(value.re);
    ScaledRealValue imaginary = absoluteValue(value.im);
    output = compareScaledNonnegative(real, imaginary) >= 0
        ? real : imaginary;
    return ScaledArithmeticStatus::Success;
}

ScaledArithmeticStatus magnitudeBound(
        const ScaledComplexValue& value,
        ScaledRealValue& output) {
    ScaledRealValue maximum;
    ScaledArithmeticStatus status =
        maxComponentMagnitude(value, maximum);
    if (status != ScaledArithmeticStatus::Success)
        return status;
    ScaledRealValue two;
    makeScaledRealValue(2.0, two);
    return scaledMultiply(maximum, two, output);
}

ExpressionScaledResidualStatus residualStatus(
        ScaledArithmeticStatus status) {
    switch (status) {
    case ScaledArithmeticStatus::Success:
        return ExpressionScaledResidualStatus::Success;
    case ScaledArithmeticStatus::Nonfinite:
        return ExpressionScaledResidualStatus::Nonfinite;
    case ScaledArithmeticStatus::Singular:
        return ExpressionScaledResidualStatus::Singular;
    case ScaledArithmeticStatus::ExponentRange:
        return ExpressionScaledResidualStatus::ExponentRange;
    }
    return ExpressionScaledResidualStatus::ExponentRange;
}

bool validShadowEncoding(const ScaledRealShadow& shadow) {
    if (!std::isfinite(shadow.mantissa))
        return shadow.exponent == 0;
    if (shadow.mantissa == 0.0)
        return shadow.exponent == 0;
    return std::abs(shadow.mantissa) >= 0.5 &&
           std::abs(shadow.mantissa) < 1.0 &&
           shadow.exponent >=
               static_cast<int64_t>(mpfr_get_emin()) &&
           shadow.exponent <=
               static_cast<int64_t>(mpfr_get_emax());
}

bool validShadowEncoding(const ScaledComplexShadow& shadow) {
    return validShadowEncoding(shadow.re) &&
           validShadowEncoding(shadow.im);
}

bool finiteShadow(const ScaledComplexShadow& shadow) {
    return std::isfinite(shadow.re.mantissa) &&
           std::isfinite(shadow.im.mantissa);
}

bool operationNeedsCompanion(ExpressionOracleOperation operation) {
    switch (operation) {
    case ExpressionOracleOperation::Sin:
    case ExpressionOracleOperation::Cos:
    case ExpressionOracleOperation::Sinh:
    case ExpressionOracleOperation::Cosh:
        return true;
    default:
        return false;
    }
}

struct SeriesValue {
    ScaledComplexValue value;
    ScaledRealValue remainder;
};

ScaledArithmeticStatus doubledMagnitude(
        const ScaledComplexValue& value,
        ScaledRealValue& output) {
    ScaledRealValue magnitude;
    ScaledArithmeticStatus status =
        magnitudeBound(value, magnitude);
    if (status != ScaledArithmeticStatus::Success)
        return status;
    ScaledRealValue two;
    makeScaledRealValue(2.0, two);
    return scaledMultiply(magnitude, two, output);
}

bool negligibleRelative(
        const ScaledComplexValue& term,
        const ScaledComplexValue& sum,
        ScaledArithmeticStatus& status) {
    ScaledRealValue termMagnitude, sumMagnitude, threshold;
    status = magnitudeBound(term, termMagnitude);
    if (status != ScaledArithmeticStatus::Success)
        return false;
    status = magnitudeBound(sum, sumMagnitude);
    if (status != ScaledArithmeticStatus::Success)
        return false;
    if (termMagnitude.isZero()) return true;
    if (sumMagnitude.isZero()) return false;
    ScaledRealValue factor;
    status = makeScaledRealValue(
        std::ldexp(1.0, -56), factor);
    if (status != ScaledArithmeticStatus::Success)
        return false;
    status = scaledMultiply(
        sumMagnitude, factor, threshold);
    if (status != ScaledArithmeticStatus::Success)
        return false;
    return compareScaledNonnegative(
        termMagnitude, threshold) <= 0;
}

ScaledArithmeticStatus checkSeriesRange(
        const ScaledComplexValue& delta) {
    ScaledArithmeticStatus status = validate(delta);
    if (status != ScaledArithmeticStatus::Success)
        return status;
    ScaledRealValue magnitude, cutoff;
    status = maxComponentMagnitude(delta, magnitude);
    if (status != ScaledArithmeticStatus::Success)
        return status;
    makeScaledRealValue(0.0625, cutoff);
    return compareScaledNonnegative(magnitude, cutoff) <= 0
        ? ScaledArithmeticStatus::Success
        : ScaledArithmeticStatus::Singular;
}

ScaledArithmeticStatus expMinusOneSeries(
        const ScaledComplexValue& delta,
        SeriesValue& output) {
    output = {};
    ScaledArithmeticStatus status = checkSeriesRange(delta);
    if (status != ScaledArithmeticStatus::Success)
        return status;
    if (delta.isZero()) return ScaledArithmeticStatus::Success;

    ScaledComplexValue sum = delta;
    ScaledComplexValue term = delta;
    for (int order = 2; order <= 20; ++order) {
        ScaledComplexValue next;
        status = scaledMultiply(term, delta, next);
        if (status != ScaledArithmeticStatus::Success)
            return status;
        status = scaledDivideByDouble(
            next, static_cast<double>(order), next);
        if (status != ScaledArithmeticStatus::Success)
            return status;
        ScaledArithmeticStatus comparisonStatus =
            ScaledArithmeticStatus::Success;
        if (negligibleRelative(
                next, sum, comparisonStatus)) {
            status = doubledMagnitude(
                next, output.remainder);
            if (status != ScaledArithmeticStatus::Success)
                return status;
            output.value = sum;
            return ScaledArithmeticStatus::Success;
        }
        if (comparisonStatus !=
                ScaledArithmeticStatus::Success)
            return comparisonStatus;
        status = scaledAdd(sum, next, sum);
        if (status != ScaledArithmeticStatus::Success)
            return status;
        term = next;
    }
    ScaledComplexValue next;
    status = scaledMultiply(term, delta, next);
    if (status != ScaledArithmeticStatus::Success)
        return status;
    status = scaledDivideByDouble(next, 21.0, next);
    if (status != ScaledArithmeticStatus::Success)
        return status;
    status = doubledMagnitude(next, output.remainder);
    if (status == ScaledArithmeticStatus::Success)
        output.value = sum;
    return status;
}

ScaledArithmeticStatus oddSeries(
        const ScaledComplexValue& delta,
        bool hyperbolic,
        SeriesValue& output) {
    output = {};
    ScaledArithmeticStatus status = checkSeriesRange(delta);
    if (status != ScaledArithmeticStatus::Success)
        return status;
    if (delta.isZero()) return ScaledArithmeticStatus::Success;

    ScaledComplexValue deltaSquared;
    status = scaledMultiply(
        delta, delta, deltaSquared);
    if (status != ScaledArithmeticStatus::Success)
        return status;
    ScaledComplexValue sum = delta;
    ScaledComplexValue term = delta;
    for (int power = 3; power <= 25; power += 2) {
        ScaledComplexValue next;
        status = scaledMultiply(
            term, deltaSquared, next);
        if (status != ScaledArithmeticStatus::Success)
            return status;
        const double divisor =
            static_cast<double>((power - 1) * power);
        status = scaledDivideByDouble(
            next, divisor, next);
        if (status != ScaledArithmeticStatus::Success)
            return status;
        if (!hyperbolic) {
            status = scaledNegate(next, next);
            if (status != ScaledArithmeticStatus::Success)
                return status;
        }
        ScaledArithmeticStatus comparisonStatus =
            ScaledArithmeticStatus::Success;
        if (negligibleRelative(
                next, sum, comparisonStatus)) {
            status = doubledMagnitude(
                next, output.remainder);
            if (status != ScaledArithmeticStatus::Success)
                return status;
            output.value = sum;
            return ScaledArithmeticStatus::Success;
        }
        if (comparisonStatus !=
                ScaledArithmeticStatus::Success)
            return comparisonStatus;
        status = scaledAdd(sum, next, sum);
        if (status != ScaledArithmeticStatus::Success)
            return status;
        term = next;
    }
    output.value = sum;
    return doubledMagnitude(term, output.remainder);
}

ScaledArithmeticStatus evenMinusOneSeries(
        const ScaledComplexValue& delta,
        bool hyperbolic,
        SeriesValue& output) {
    output = {};
    ScaledArithmeticStatus status = checkSeriesRange(delta);
    if (status != ScaledArithmeticStatus::Success)
        return status;
    if (delta.isZero()) return ScaledArithmeticStatus::Success;

    ScaledComplexValue deltaSquared;
    status = scaledMultiply(
        delta, delta, deltaSquared);
    if (status != ScaledArithmeticStatus::Success)
        return status;
    ScaledComplexValue term;
    status = scaledDivideByDouble(
        deltaSquared, 2.0, term);
    if (status != ScaledArithmeticStatus::Success)
        return status;
    if (!hyperbolic) {
        status = scaledNegate(term, term);
        if (status != ScaledArithmeticStatus::Success)
            return status;
    }
    ScaledComplexValue sum = term;
    for (int power = 4; power <= 24; power += 2) {
        ScaledComplexValue next;
        status = scaledMultiply(
            term, deltaSquared, next);
        if (status != ScaledArithmeticStatus::Success)
            return status;
        const double divisor =
            static_cast<double>((power - 1) * power);
        status = scaledDivideByDouble(
            next, divisor, next);
        if (status != ScaledArithmeticStatus::Success)
            return status;
        if (!hyperbolic) {
            status = scaledNegate(next, next);
            if (status != ScaledArithmeticStatus::Success)
                return status;
        }
        ScaledArithmeticStatus comparisonStatus =
            ScaledArithmeticStatus::Success;
        if (negligibleRelative(
                next, sum, comparisonStatus)) {
            status = doubledMagnitude(
                next, output.remainder);
            if (status != ScaledArithmeticStatus::Success)
                return status;
            output.value = sum;
            return ScaledArithmeticStatus::Success;
        }
        if (comparisonStatus !=
                ScaledArithmeticStatus::Success)
            return comparisonStatus;
        status = scaledAdd(sum, next, sum);
        if (status != ScaledArithmeticStatus::Success)
            return status;
        term = next;
    }
    output.value = sum;
    return doubledMagnitude(term, output.remainder);
}

ScaledArithmeticStatus combineRemainder(
        const ScaledComplexValue& coefficient,
        const ScaledRealValue& seriesRemainder,
        ScaledRealValue& output) {
    if (seriesRemainder.isZero()) {
        output = {};
        return ScaledArithmeticStatus::Success;
    }
    ScaledRealValue coefficientMagnitude;
    ScaledArithmeticStatus status =
        magnitudeBound(coefficient, coefficientMagnitude);
    if (status != ScaledArithmeticStatus::Success)
        return status;
    return scaledMultiply(
        coefficientMagnitude, seriesRemainder, output);
}

ScaledArithmeticStatus addRemainders(
        const ScaledRealValue& left,
        const ScaledRealValue& right,
        ScaledRealValue& output) {
    return scaledAdd(left, right, output);
}

} // namespace

ScaledArithmeticStatus makeScaledRealValue(
        double value, ScaledRealValue& output) {
    return normalize(value, 0, output);
}

ScaledArithmeticStatus makeScaledComplexValue(
        Complex value, ScaledComplexValue& output) {
    ScaledComplexValue result;
    ScaledArithmeticStatus status =
        makeScaledRealValue(value.real(), result.re);
    if (status != ScaledArithmeticStatus::Success) {
        output = result;
        return status;
    }
    status = makeScaledRealValue(
        value.imag(), result.im);
    output = result;
    return status;
}

ScaledArithmeticStatus makeScaledRealValue(
        mpfr_srcptr value, ScaledRealValue& output) {
    ScaledRealShadow shadow;
    if (!makeScaledRealShadow(value, shadow))
        return ScaledArithmeticStatus::ExponentRange;
    return makeScaledRealValue(shadow, output);
}

ScaledArithmeticStatus makeScaledComplexValue(
        const MpfrComplex& value,
        ScaledComplexValue& output) {
    ScaledComplexShadow shadow;
    if (!makeScaledComplexShadow(value, shadow))
        return ScaledArithmeticStatus::ExponentRange;
    return makeScaledComplexValue(shadow, output);
}

ScaledArithmeticStatus makeScaledRealValue(
        const ScaledRealShadow& shadow,
        ScaledRealValue& output) {
    return normalize(
        shadow.mantissa, shadow.exponent, output);
}

ScaledArithmeticStatus makeScaledComplexValue(
        const ScaledComplexShadow& shadow,
        ScaledComplexValue& output) {
    ScaledComplexValue result;
    ScaledArithmeticStatus status =
        makeScaledRealValue(shadow.re, result.re);
    if (status != ScaledArithmeticStatus::Success) {
        output = result;
        return status;
    }
    status = makeScaledRealValue(
        shadow.im, result.im);
    output = result;
    return status;
}

ScaledArithmeticStatus makeScaledComplexValue(
        const ScaledComplexShadow& primary,
        const ScaledComplexShadow& defect,
        ScaledComplexValue& output) {
    ScaledComplexValue left, right;
    ScaledArithmeticStatus status =
        makeScaledComplexValue(primary, left);
    if (status != ScaledArithmeticStatus::Success)
        return status;
    status = makeScaledComplexValue(defect, right);
    if (status != ScaledArithmeticStatus::Success)
        return status;
    return scaledAdd(left, right, output);
}

bool setMpfrFromScaledValue(
        mpfr_ptr output, const ScaledRealValue& value,
        mpfr_rnd_t rounding) {
    if (validate(value) != ScaledArithmeticStatus::Success)
        return false;
    ScaledRealShadow shadow;
    shadow.mantissa = value.mantissa;
    shadow.exponent = value.exponent;
    return setMpfrFromScaledShadow(
        output, shadow, rounding);
}

bool setMpfrFromScaledValue(
        MpfrComplex& output,
        const ScaledComplexValue& value,
        mpfr_rnd_t rounding) {
    return setMpfrFromScaledValue(
               output.re, value.re, rounding) &&
           setMpfrFromScaledValue(
               output.im, value.im, rounding);
}

bool scaledValueToDouble(
        const ScaledRealValue& value, double& output) {
    if (validate(value) != ScaledArithmeticStatus::Success)
        return false;
    if (value.isZero()) {
        output = value.mantissa;
        return true;
    }
    if (value.exponent < -1074 ||
        value.exponent > 1024)
        return false;
    double converted = std::ldexp(
        value.mantissa,
        static_cast<int>(value.exponent));
    if (!std::isfinite(converted) ||
        converted == 0.0)
        return false;
    output = converted;
    return true;
}

bool scaledValueToDouble(
        const ScaledComplexValue& value,
        Complex& output) {
    double real = 0.0, imaginary = 0.0;
    if (!scaledValueToDouble(value.re, real) ||
        !scaledValueToDouble(value.im, imaginary))
        return false;
    output = { real, imaginary };
    return true;
}

ScaledArithmeticStatus scaledNegate(
        const ScaledRealValue& value,
        ScaledRealValue& output) {
    ScaledArithmeticStatus status = validate(value);
    if (status != ScaledArithmeticStatus::Success)
        return status;
    output = value;
    output.mantissa = -output.mantissa;
    return ScaledArithmeticStatus::Success;
}

ScaledArithmeticStatus scaledNegate(
        const ScaledComplexValue& value,
        ScaledComplexValue& output) {
    ScaledComplexValue result;
    ScaledArithmeticStatus status =
        scaledNegate(value.re, result.re);
    if (status != ScaledArithmeticStatus::Success)
        return status;
    status = scaledNegate(value.im, result.im);
    if (status == ScaledArithmeticStatus::Success)
        output = result;
    return status;
}

ScaledArithmeticStatus scaledAdd(
        const ScaledRealValue& left,
        const ScaledRealValue& right,
        ScaledRealValue& output) {
    ScaledArithmeticStatus status = validate(left);
    if (status != ScaledArithmeticStatus::Success)
        return status;
    status = validate(right);
    if (status != ScaledArithmeticStatus::Success)
        return status;
    if (left.isZero() && right.isZero())
        return normalize(
            left.mantissa + right.mantissa, 0, output);
    if (left.isZero()) {
        output = right;
        return ScaledArithmeticStatus::Success;
    }
    if (right.isZero()) {
        output = left;
        return ScaledArithmeticStatus::Success;
    }

    const ScaledRealValue* large = &left;
    const ScaledRealValue* small = &right;
    if (right.exponent > left.exponent)
        std::swap(large, small);
    if (large->exponent >
            std::numeric_limits<int64_t>::min() + 1075 &&
        small->exponent < large->exponent - 1075) {
        output = *large;
        return ScaledArithmeticStatus::Success;
    }
    const int64_t difference =
        small->exponent - large->exponent;
    if (difference < -1075) {
        output = *large;
        return ScaledArithmeticStatus::Success;
    }
    const double sum =
        large->mantissa +
        std::ldexp(
            small->mantissa,
            static_cast<int>(difference));
    return normalize(sum, large->exponent, output);
}

ScaledArithmeticStatus scaledSubtract(
        const ScaledRealValue& left,
        const ScaledRealValue& right,
        ScaledRealValue& output) {
    ScaledRealValue negative;
    ScaledArithmeticStatus status =
        scaledNegate(right, negative);
    if (status != ScaledArithmeticStatus::Success)
        return status;
    return scaledAdd(left, negative, output);
}

ScaledArithmeticStatus scaledMultiply(
        const ScaledRealValue& left,
        const ScaledRealValue& right,
        ScaledRealValue& output) {
    ScaledArithmeticStatus status = validate(left);
    if (status != ScaledArithmeticStatus::Success)
        return status;
    status = validate(right);
    if (status != ScaledArithmeticStatus::Success)
        return status;
    if (left.isZero() || right.isZero()) {
        output = {};
        output.mantissa = signedZero(
            std::signbit(left.mantissa) !=
            std::signbit(right.mantissa));
        return ScaledArithmeticStatus::Success;
    }
    ScaledRealValue factor;
    ScaledArithmeticStatus factorStatus = normalize(
        left.mantissa * right.mantissa, 0, factor);
    if (factorStatus != ScaledArithmeticStatus::Success)
        return factorStatus;
    int64_t exponent = 0;
    if (!checkedAddThreeExponents(
            left.exponent, right.exponent,
            factor.exponent, exponent))
        return ScaledArithmeticStatus::ExponentRange;
    output = factor;
    output.exponent = exponent;
    return ScaledArithmeticStatus::Success;
}

ScaledArithmeticStatus scaledDivideByDouble(
        const ScaledRealValue& value, double divisor,
        ScaledRealValue& output) {
    ScaledArithmeticStatus status = validate(value);
    if (status != ScaledArithmeticStatus::Success)
        return status;
    if (!std::isfinite(divisor))
        return ScaledArithmeticStatus::Nonfinite;
    if (divisor == 0.0)
        return ScaledArithmeticStatus::Singular;
    if (value.isZero()) {
        output = {};
        output.mantissa = signedZero(
            std::signbit(value.mantissa) !=
            std::signbit(divisor));
        return ScaledArithmeticStatus::Success;
    }
    int divisorExponent = 0;
    double divisorMantissa =
        std::frexp(divisor, &divisorExponent);
    ScaledRealValue quotient;
    ScaledArithmeticStatus quotientStatus = normalize(
        value.mantissa / divisorMantissa, 0, quotient);
    if (quotientStatus != ScaledArithmeticStatus::Success)
        return quotientStatus;
    int64_t exponent = 0;
    if (!checkedAddThreeExponents(
            value.exponent,
            -static_cast<int64_t>(divisorExponent),
            quotient.exponent, exponent))
        return ScaledArithmeticStatus::ExponentRange;
    output = quotient;
    output.exponent = exponent;
    return ScaledArithmeticStatus::Success;
}

ScaledArithmeticStatus scaledAdd(
        const ScaledComplexValue& left,
        const ScaledComplexValue& right,
        ScaledComplexValue& output) {
    ScaledComplexValue result;
    ScaledArithmeticStatus status =
        scaledAdd(left.re, right.re, result.re);
    if (status != ScaledArithmeticStatus::Success)
        return status;
    status = scaledAdd(
        left.im, right.im, result.im);
    if (status == ScaledArithmeticStatus::Success)
        output = result;
    return status;
}

ScaledArithmeticStatus scaledSubtract(
        const ScaledComplexValue& left,
        const ScaledComplexValue& right,
        ScaledComplexValue& output) {
    ScaledComplexValue result;
    ScaledArithmeticStatus status =
        scaledSubtract(left.re, right.re, result.re);
    if (status != ScaledArithmeticStatus::Success)
        return status;
    status = scaledSubtract(
        left.im, right.im, result.im);
    if (status == ScaledArithmeticStatus::Success)
        output = result;
    return status;
}

ScaledArithmeticStatus scaledMultiply(
        const ScaledComplexValue& left,
        const ScaledComplexValue& right,
        ScaledComplexValue& output) {
    ScaledRealValue realReal, imaginaryImaginary;
    ScaledRealValue realImaginary, imaginaryReal;
    ScaledComplexValue result;
    ScaledArithmeticStatus status =
        scaledMultiply(left.re, right.re, realReal);
    if (status != ScaledArithmeticStatus::Success)
        return status;
    status = scaledMultiply(
        left.im, right.im, imaginaryImaginary);
    if (status != ScaledArithmeticStatus::Success)
        return status;
    status = scaledSubtract(
        realReal, imaginaryImaginary, result.re);
    if (status != ScaledArithmeticStatus::Success)
        return status;
    status = scaledMultiply(
        left.re, right.im, realImaginary);
    if (status != ScaledArithmeticStatus::Success)
        return status;
    status = scaledMultiply(
        left.im, right.re, imaginaryReal);
    if (status != ScaledArithmeticStatus::Success)
        return status;
    status = scaledAdd(
        realImaginary, imaginaryReal, result.im);
    if (status == ScaledArithmeticStatus::Success)
        output = result;
    return status;
}

ScaledArithmeticStatus scaledMultiplyByDouble(
        const ScaledComplexValue& value,
        double multiplier,
        ScaledComplexValue& output) {
    ScaledRealValue scalar;
    ScaledArithmeticStatus status =
        makeScaledRealValue(multiplier, scalar);
    if (status != ScaledArithmeticStatus::Success)
        return status;
    ScaledComplexValue result;
    status = scaledMultiply(
        value.re, scalar, result.re);
    if (status != ScaledArithmeticStatus::Success)
        return status;
    status = scaledMultiply(
        value.im, scalar, result.im);
    if (status == ScaledArithmeticStatus::Success)
        output = result;
    return status;
}

ScaledArithmeticStatus scaledDivideByDouble(
        const ScaledComplexValue& value, double divisor,
        ScaledComplexValue& output) {
    ScaledComplexValue result;
    ScaledArithmeticStatus status =
        scaledDivideByDouble(
            value.re, divisor, result.re);
    if (status != ScaledArithmeticStatus::Success)
        return status;
    status = scaledDivideByDouble(
        value.im, divisor, result.im);
    if (status == ScaledArithmeticStatus::Success)
        output = result;
    return status;
}

ScaledArithmeticStatus scaledNormSquared(
        const ScaledComplexValue& value,
        ScaledRealValue& output) {
    ScaledRealValue realSquared, imaginarySquared;
    ScaledArithmeticStatus status =
        scaledMultiply(value.re, value.re, realSquared);
    if (status != ScaledArithmeticStatus::Success)
        return status;
    status = scaledMultiply(
        value.im, value.im, imaginarySquared);
    if (status != ScaledArithmeticStatus::Success)
        return status;
    return scaledAdd(
        realSquared, imaginarySquared, output);
}

int compareScaledNonnegative(
        const ScaledRealValue& left,
        const ScaledRealValue& right) {
    if (left.isZero())
        return right.isZero() ? 0 : -1;
    if (right.isZero()) return 1;
    if (left.exponent != right.exponent)
        return left.exponent < right.exponent ? -1 : 1;
    const double leftMagnitude = std::abs(left.mantissa);
    const double rightMagnitude = std::abs(right.mantissa);
    if (leftMagnitude == rightMagnitude) return 0;
    return leftMagnitude < rightMagnitude ? -1 : 1;
}

const char* expressionScaledResidualStatusName(
        ExpressionScaledResidualStatus status) {
    switch (status) {
    case ExpressionScaledResidualStatus::Success:
        return "success";
    case ExpressionScaledResidualStatus::BranchUncertain:
        return "branch-uncertain";
    case ExpressionScaledResidualStatus::Singular:
        return "singular";
    case ExpressionScaledResidualStatus::Unsupported:
        return "unsupported";
    case ExpressionScaledResidualStatus::Nonfinite:
        return "nonfinite";
    case ExpressionScaledResidualStatus::ExponentRange:
        return "exponent-range";
    case ExpressionScaledResidualStatus::InvalidTape:
        return "invalid-tape";
    case ExpressionScaledResidualStatus::InvalidInput:
        return "invalid-input";
    }
    return "invalid-tape";
}

void ExpressionScaledResidualEvaluator::reset() {
    _program = nullptr;
    _reference = nullptr;
    _states.clear();
    _layoutStack.clear();
    _preparationStatus =
        ExpressionScaledResidualStatus::InvalidTape;
    _error.clear();
    _ready = false;
}

bool ExpressionScaledResidualEvaluator::prepare(
        const ExpressionProgram& program,
        const ExpressionReferenceOrbitResult& reference) {
    reset();
    auto fail = [&](const char* message) {
        _error = message;
        return false;
    };
    if (!program.valid() || !reference.valid ||
        reference.status !=
            ExpressionReferenceBuildStatus::Success)
        return fail("program or reference is invalid");
    if (reference.programSemanticHash !=
            program.semanticHash() ||
        reference.programSource != program.source())
        return fail("reference semantic identity mismatch");
    if (reference.sampleCount !=
            reference.samples.size())
        return fail("reference sample count mismatch");
    if (program.instructionCount() == 0 ||
        program.instructionCount() >
            std::numeric_limits<uint16_t>::max())
        return fail("program cannot be represented by the tape");
    const ScaledComplexShadow* referenceShadows[] = {
        &reference.c, &reference.cDefect,
        &reference.z0, &reference.z0Defect,
        &reference.pixel, &reference.pixelDefect,
        &reference.initialZ, &reference.initialZDefect
    };
    for (const ScaledComplexShadow* shadow :
         referenceShadows)
        if (!validShadowEncoding(*shadow))
            return fail("malformed compact reference input");

    auto operationOf = [](ExpressionProgram::Op op) {
        using Op = ExpressionProgram::Op;
        switch (op) {
        case Op::Constant:
            return ExpressionOracleOperation::Constant;
        case Op::Z: return ExpressionOracleOperation::Z;
        case Op::C: return ExpressionOracleOperation::C;
        case Op::Z0: return ExpressionOracleOperation::Z0;
        case Op::Iteration:
            return ExpressionOracleOperation::Iteration;
        case Op::Parameter:
            return ExpressionOracleOperation::Parameter;
        case Op::OrbitInvariant:
            return ExpressionOracleOperation::OrbitInvariant;
        case Op::Negate:
            return ExpressionOracleOperation::Negate;
        case Op::Add: return ExpressionOracleOperation::Add;
        case Op::Subtract:
            return ExpressionOracleOperation::Subtract;
        case Op::Multiply:
            return ExpressionOracleOperation::Multiply;
        case Op::Divide:
            return ExpressionOracleOperation::Divide;
        case Op::Power:
            return ExpressionOracleOperation::Power;
        case Op::Square:
            return ExpressionOracleOperation::Square;
        case Op::Sin: return ExpressionOracleOperation::Sin;
        case Op::Cos: return ExpressionOracleOperation::Cos;
        case Op::Tan: return ExpressionOracleOperation::Tan;
        case Op::Sinh:
            return ExpressionOracleOperation::Sinh;
        case Op::Cosh:
            return ExpressionOracleOperation::Cosh;
        case Op::Tanh:
            return ExpressionOracleOperation::Tanh;
        case Op::Exp: return ExpressionOracleOperation::Exp;
        case Op::Log: return ExpressionOracleOperation::Log;
        case Op::Log10:
            return ExpressionOracleOperation::Log10;
        case Op::Sqrt:
            return ExpressionOracleOperation::Sqrt;
        case Op::Abs: return ExpressionOracleOperation::Abs;
        case Op::Norm:
            return ExpressionOracleOperation::Norm;
        case Op::Arg: return ExpressionOracleOperation::Arg;
        case Op::Conjugate:
            return ExpressionOracleOperation::Conjugate;
        case Op::Real: return ExpressionOracleOperation::Real;
        case Op::Imaginary:
            return ExpressionOracleOperation::Imaginary;
        case Op::MakeComplex:
            return ExpressionOracleOperation::MakeComplex;
        case Op::Polar:
            return ExpressionOracleOperation::Polar;
        }
        return ExpressionOracleOperation::OrbitInvariant;
    };

    _layoutStack.reserve(program.stackDepth());
    for (size_t sampleIndex = 0;
         sampleIndex < reference.samples.size();
         ++sampleIndex) {
        const ExpressionReferenceSample& sample =
            reference.samples[sampleIndex];
        if (sample.iteration !=
                static_cast<int>(sampleIndex) ||
            sample.tapeCount !=
                program.instructionCount() ||
            sample.tapeOffset > reference.tape.size() ||
            sample.tapeCount >
                reference.tape.size() -
                    static_cast<size_t>(sample.tapeOffset))
            return fail("reference sample range mismatch");
        if (!validShadowEncoding(sample.z) ||
            !validShadowEncoding(sample.zDefect) ||
            !validShadowEncoding(sample.next) ||
            !validShadowEncoding(sample.rootDefect))
            return fail("malformed compact reference sample");
        _layoutStack.clear();
        for (size_t instructionIndex = 0;
             instructionIndex < program._code.size();
             ++instructionIndex) {
            const ExpressionProgram::Instruction& instruction =
                program._code[instructionIndex];
            const ExpressionReferenceTapeNode& node =
                reference.tape[
                    static_cast<size_t>(sample.tapeOffset) +
                    instructionIndex];
            const int operands =
                ExpressionProgram::operandCount(instruction.op);
            if (operands < 0 ||
                _layoutStack.size() <
                    static_cast<size_t>(operands))
                return fail("malformed expression layout");
            uint16_t expectedLeft = UINT16_MAX;
            uint16_t expectedRight = UINT16_MAX;
            if (operands == 1) {
                expectedLeft = _layoutStack.back();
                _layoutStack.pop_back();
            } else if (operands == 2) {
                expectedRight = _layoutStack.back();
                _layoutStack.pop_back();
                expectedLeft = _layoutStack.back();
                _layoutStack.pop_back();
            }
            if (node.operation !=
                    operationOf(instruction.op) ||
                node.argument != instruction.argument ||
                node.leftNode != expectedLeft ||
                node.rightNode != expectedRight)
                return fail("reference tape layout mismatch");
            if (!validShadowEncoding(node.output) ||
                !validShadowEncoding(node.outputDefect))
                return fail("malformed reference output");
            const bool hasAuxiliary =
                (node.flags &
                 (OracleTraceHasCompanion |
                  OracleTraceHasDenominator |
                  OracleTraceHasLogarithmBase)) != 0;
            if (operationNeedsCompanion(node.operation) &&
                !(node.flags & OracleTraceHasCompanion))
                return fail("required companion is missing");
            if (hasAuxiliary &&
                (!validShadowEncoding(node.auxiliary) ||
                 !validShadowEncoding(
                     node.auxiliaryDefect)))
                return fail("malformed reference companion");
            if (operationNeedsCompanion(node.operation) &&
                !(node.flags & OracleTraceUndefined) &&
                (!finiteShadow(node.auxiliary) ||
                 !finiteShadow(node.auxiliaryDefect)))
                return fail("required companion is nonfinite");
            _layoutStack.push_back(
                static_cast<uint16_t>(instructionIndex));
        }
        if (_layoutStack.size() != 1 ||
            sample.rootNode != _layoutStack.back())
            return fail("reference tape root mismatch");
    }

    _program = &program;
    _reference = &reference;
    _states.resize(program.instructionCount());
    _preparationStatus =
        ExpressionScaledResidualStatus::Success;
    _ready = true;
    return true;
}

ExpressionScaledResidualResult
ExpressionScaledResidualEvaluator::evaluate(
        size_t sampleIndex,
        const ExpressionScaledResidualInput& input) {
    ExpressionScaledResidualResult result;
    result.status = _preparationStatus;
    if (!_ready || !_program || !_reference)
        return result;
    if (sampleIndex >= _reference->samples.size()) {
        result.status =
            ExpressionScaledResidualStatus::InvalidTape;
        return result;
    }
    const ExpressionReferenceSample& sample =
        _reference->samples[sampleIndex];
    if (input.iteration != sample.iteration) {
        result.status =
            ExpressionScaledResidualStatus::InvalidInput;
        return result;
    }
    auto validateInput = [&](const ScaledComplexValue& value) {
        ScaledArithmeticStatus status = validate(value);
        if (status == ScaledArithmeticStatus::Nonfinite)
            result.status =
                ExpressionScaledResidualStatus::Nonfinite;
        else if (status != ScaledArithmeticStatus::Success)
            result.status =
                ExpressionScaledResidualStatus::InvalidInput;
        return status == ScaledArithmeticStatus::Success;
    };
    if (!validateInput(input.z) ||
        !validateInput(input.c) ||
        !validateInput(input.z0))
        return result;
    for (const ScaledComplexValue& parameter :
         input.parameters)
        if (!validateInput(parameter))
            return result;

    const size_t offset =
        static_cast<size_t>(sample.tapeOffset);
    auto failArithmetic = [&](
            ScaledArithmeticStatus status) {
        result.status = residualStatus(status);
        return result;
    };
    auto nodeBase = [&](uint16_t localNode,
                        ScaledComplexValue& output) {
        const ExpressionReferenceTapeNode& node =
            _reference->tape[offset + localNode];
        return makeScaledComplexValue(
            node.output, node.outputDefect, output);
    };
    auto updateRemainder = [&](
            const ScaledRealValue& estimate) {
        if (compareScaledNonnegative(
                estimate, result.remainderEstimate) > 0)
            result.remainderEstimate = estimate;
    };
    auto nonlinearResult = [&](
            const ExpressionReferenceTapeNode& node,
            const ScaledComplexValue& delta,
            ScaledComplexValue& output) {
        ScaledComplexValue base, companion;
        SeriesValue first, second;
        ScaledArithmeticStatus status =
            makeScaledComplexValue(
                node.output, node.outputDefect, base);
        if (status != ScaledArithmeticStatus::Success)
            return status;
        if (node.operation !=
                ExpressionOracleOperation::Exp) {
            status = makeScaledComplexValue(
                node.auxiliary,
                node.auxiliaryDefect,
                companion);
            if (status != ScaledArithmeticStatus::Success)
                return status;
        }

        ScaledComplexValue left, right;
        ScaledRealValue leftRemainder, rightRemainder;
        switch (node.operation) {
        case ExpressionOracleOperation::Exp:
            status = expMinusOneSeries(delta, first);
            if (status == ScaledArithmeticStatus::Singular)
                return status;
            if (status != ScaledArithmeticStatus::Success)
                return status;
            status = scaledMultiply(
                base, first.value, output);
            if (status != ScaledArithmeticStatus::Success)
                return status;
            status = combineRemainder(
                base, first.remainder, leftRemainder);
            if (status == ScaledArithmeticStatus::Success)
                updateRemainder(leftRemainder);
            break;
        case ExpressionOracleOperation::Sin:
        case ExpressionOracleOperation::Cos:
            status = evenMinusOneSeries(
                delta, false, first);
            if (status != ScaledArithmeticStatus::Success)
                return status;
            status = oddSeries(delta, false, second);
            if (status != ScaledArithmeticStatus::Success)
                return status;
            status = scaledMultiply(
                base, first.value, left);
            if (status != ScaledArithmeticStatus::Success)
                return status;
            status = scaledMultiply(
                companion, second.value, right);
            if (status != ScaledArithmeticStatus::Success)
                return status;
            if (node.operation ==
                    ExpressionOracleOperation::Cos) {
                status = scaledNegate(right, right);
                if (status != ScaledArithmeticStatus::Success)
                    return status;
            }
            status = scaledAdd(left, right, output);
            if (status != ScaledArithmeticStatus::Success)
                return status;
            status = combineRemainder(
                base, first.remainder, leftRemainder);
            if (status != ScaledArithmeticStatus::Success)
                return status;
            status = combineRemainder(
                companion, second.remainder,
                rightRemainder);
            if (status != ScaledArithmeticStatus::Success)
                return status;
            status = addRemainders(
                leftRemainder, rightRemainder,
                leftRemainder);
            if (status == ScaledArithmeticStatus::Success)
                updateRemainder(leftRemainder);
            break;
        case ExpressionOracleOperation::Sinh:
        case ExpressionOracleOperation::Cosh:
            status = evenMinusOneSeries(
                delta, true, first);
            if (status != ScaledArithmeticStatus::Success)
                return status;
            status = oddSeries(delta, true, second);
            if (status != ScaledArithmeticStatus::Success)
                return status;
            status = scaledMultiply(
                base, first.value, left);
            if (status != ScaledArithmeticStatus::Success)
                return status;
            status = scaledMultiply(
                companion, second.value, right);
            if (status != ScaledArithmeticStatus::Success)
                return status;
            status = scaledAdd(left, right, output);
            if (status != ScaledArithmeticStatus::Success)
                return status;
            status = combineRemainder(
                base, first.remainder, leftRemainder);
            if (status != ScaledArithmeticStatus::Success)
                return status;
            status = combineRemainder(
                companion, second.remainder,
                rightRemainder);
            if (status != ScaledArithmeticStatus::Success)
                return status;
            status = addRemainders(
                leftRemainder, rightRemainder,
                leftRemainder);
            if (status == ScaledArithmeticStatus::Success)
                updateRemainder(leftRemainder);
            break;
        default:
            return ScaledArithmeticStatus::Singular;
        }
        if (status == ScaledArithmeticStatus::Success &&
            !delta.isZero()) {
            result.uncertified = true;
            ++result.seriesOperationCount;
        }
        return status;
    };

    for (size_t local = 0;
         local < sample.tapeCount; ++local) {
        const ExpressionReferenceTapeNode& node =
            _reference->tape[offset + local];
        ++result.operationCount;
        if (node.flags & OracleTraceUndefined) {
            if (node.flags & OracleTraceSingularPoint)
                result.status =
                    ExpressionScaledResidualStatus::Singular;
            else if (node.flags &
                     OracleTraceBranchSensitive)
                result.status =
                    ExpressionScaledResidualStatus::BranchUncertain;
            else
                result.status =
                    ExpressionScaledResidualStatus::Nonfinite;
            return result;
        }
        if (!finiteShadow(node.output) ||
            !finiteShadow(node.outputDefect)) {
            result.status =
                ExpressionScaledResidualStatus::Nonfinite;
            return result;
        }
        ScaledComplexValue& output =
            _states[local].residual;
        output = {};
        ScaledArithmeticStatus status =
            ScaledArithmeticStatus::Success;
        const ScaledComplexValue* left =
            node.leftNode == UINT16_MAX
                ? nullptr
                : &_states[node.leftNode].residual;
        const ScaledComplexValue* right =
            node.rightNode == UINT16_MAX
                ? nullptr
                : &_states[node.rightNode].residual;
        switch (node.operation) {
        case ExpressionOracleOperation::Constant:
        case ExpressionOracleOperation::Iteration:
            break;
        case ExpressionOracleOperation::Z:
            output = input.z;
            break;
        case ExpressionOracleOperation::C:
            output = input.c;
            break;
        case ExpressionOracleOperation::Z0:
            output = input.z0;
            break;
        case ExpressionOracleOperation::Parameter:
            if (node.argument >= input.parameters.size()) {
                result.status =
                    ExpressionScaledResidualStatus::InvalidTape;
                return result;
            }
            output = input.parameters[node.argument];
            break;
        case ExpressionOracleOperation::Negate:
            status = scaledNegate(*left, output);
            break;
        case ExpressionOracleOperation::Add:
            status = scaledAdd(*left, *right, output);
            break;
        case ExpressionOracleOperation::Subtract:
            status = scaledSubtract(
                *left, *right, output);
            break;
        case ExpressionOracleOperation::Multiply: {
            ScaledComplexValue leftBase, rightBase;
            ScaledComplexValue first, second, cross;
            status = nodeBase(
                node.leftNode, leftBase);
            if (status != ScaledArithmeticStatus::Success)
                break;
            status = nodeBase(
                node.rightNode, rightBase);
            if (status != ScaledArithmeticStatus::Success)
                break;
            status = scaledMultiply(
                leftBase, *right, first);
            if (status != ScaledArithmeticStatus::Success)
                break;
            status = scaledMultiply(
                rightBase, *left, second);
            if (status != ScaledArithmeticStatus::Success)
                break;
            status = scaledMultiply(
                *left, *right, cross);
            if (status != ScaledArithmeticStatus::Success)
                break;
            status = scaledAdd(first, second, output);
            if (status != ScaledArithmeticStatus::Success)
                break;
            status = scaledAdd(output, cross, output);
            break;
        }
        case ExpressionOracleOperation::Square: {
            ScaledComplexValue base, linear, quadratic;
            status = nodeBase(node.leftNode, base);
            if (status != ScaledArithmeticStatus::Success)
                break;
            status = scaledMultiply(
                base, *left, linear);
            if (status != ScaledArithmeticStatus::Success)
                break;
            status = scaledMultiplyByDouble(
                linear, 2.0, linear);
            if (status != ScaledArithmeticStatus::Success)
                break;
            status = scaledMultiply(
                *left, *left, quadratic);
            if (status != ScaledArithmeticStatus::Success)
                break;
            status = scaledAdd(
                linear, quadratic, output);
            break;
        }
        case ExpressionOracleOperation::Exp:
        case ExpressionOracleOperation::Sin:
        case ExpressionOracleOperation::Cos:
        case ExpressionOracleOperation::Sinh:
        case ExpressionOracleOperation::Cosh:
            status = nonlinearResult(
                node, *left, output);
            if (status == ScaledArithmeticStatus::Singular) {
                result.status =
                    ExpressionScaledResidualStatus::Unsupported;
                return result;
            }
            break;
        case ExpressionOracleOperation::Conjugate:
            output.re = left->re;
            status = scaledNegate(
                left->im, output.im);
            break;
        case ExpressionOracleOperation::Real:
            output.re = left->re;
            output.im = {};
            break;
        case ExpressionOracleOperation::Imaginary:
            output.re = left->im;
            output.im = {};
            break;
        case ExpressionOracleOperation::MakeComplex:
            output.re = left->re;
            output.im = right->re;
            break;
        case ExpressionOracleOperation::Norm: {
            ScaledComplexValue base, product;
            ScaledRealValue deltaNorm, linear;
            status = nodeBase(node.leftNode, base);
            if (status != ScaledArithmeticStatus::Success)
                break;
            ScaledComplexValue conjugateBase = base;
            status = scaledNegate(
                conjugateBase.im,
                conjugateBase.im);
            if (status != ScaledArithmeticStatus::Success)
                break;
            status = scaledMultiply(
                conjugateBase, *left, product);
            if (status != ScaledArithmeticStatus::Success)
                break;
            status = scaledMultiplyByDouble(
                product, 2.0, product);
            if (status != ScaledArithmeticStatus::Success)
                break;
            linear = product.re;
            status = scaledNormSquared(
                *left, deltaNorm);
            if (status != ScaledArithmeticStatus::Success)
                break;
            status = scaledAdd(
                linear, deltaNorm, output.re);
            output.im = {};
            break;
        }
        case ExpressionOracleOperation::Divide:
            result.status =
                (node.flags & OracleTraceSingularPoint)
                    ? ExpressionScaledResidualStatus::Singular
                    : ExpressionScaledResidualStatus::Unsupported;
            return result;
        case ExpressionOracleOperation::Log:
        case ExpressionOracleOperation::Log10:
        case ExpressionOracleOperation::Sqrt:
        case ExpressionOracleOperation::Power:
        case ExpressionOracleOperation::Arg:
            result.status =
                (node.flags & OracleTraceSingularPoint)
                    ? ExpressionScaledResidualStatus::Singular
                    : ExpressionScaledResidualStatus::BranchUncertain;
            return result;
        case ExpressionOracleOperation::Tan:
        case ExpressionOracleOperation::Tanh:
            result.status =
                (node.flags & OracleTraceSingularPoint)
                    ? ExpressionScaledResidualStatus::Singular
                    : ExpressionScaledResidualStatus::BranchUncertain;
            return result;
        case ExpressionOracleOperation::Abs:
        case ExpressionOracleOperation::Polar:
        case ExpressionOracleOperation::OrbitInvariant:
            result.status =
                ExpressionScaledResidualStatus::Unsupported;
            return result;
        }
        if (status != ScaledArithmeticStatus::Success)
            return failArithmetic(status);
    }
    result.residual =
        _states[sample.rootNode].residual;
    result.status =
        ExpressionScaledResidualStatus::Success;
    return result;
}

ScaledArithmeticStatus makeExpressionResidualNextPrimaryState(
        const ExpressionReferenceSample& sample,
        const ScaledComplexValue& exactResidual,
        ExpressionScaledPrimaryRelativeState& primaryRelative) {
    ScaledArithmeticStatus status =
        makeScaledComplexValue(
            sample.rootDefect,
            primaryRelative.referenceDefect);
    if (status != ScaledArithmeticStatus::Success)
        return status;
    primaryRelative.residual = exactResidual;
    return ScaledArithmeticStatus::Success;
}

ScaledArithmeticStatus resetExpressionResidualToExactSample(
        const ExpressionReferenceSample& sample,
        const ExpressionScaledPrimaryRelativeState& primaryRelative,
        ScaledComplexValue& exactResidual) {
    ScaledComplexValue defect;
    ScaledArithmeticStatus status =
        makeScaledComplexValue(
            sample.zDefect, defect);
    if (status != ScaledArithmeticStatus::Success)
        return status;
    ScaledComplexValue defectCorrection;
    status = scaledSubtract(
        primaryRelative.referenceDefect,
        defect, defectCorrection);
    if (status != ScaledArithmeticStatus::Success)
        return status;
    return scaledAdd(
        defectCorrection,
        primaryRelative.residual,
        exactResidual);
}

} // namespace formula
