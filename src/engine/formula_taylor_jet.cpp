#include "formula_taylor_jet.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <new>
#include <utility>

namespace formula {

namespace {

using Clock = std::chrono::steady_clock;

bool validRadius(const ScaledRealValue& value) {
    return value.isNormalized() && value.mantissa >= 0.0;
}

ScaledRealValue absoluteValue(ScaledRealValue value) {
    value.mantissa = std::abs(value.mantissa);
    return value;
}

ExpressionTaylorJetStatus statusForArithmetic(
        ScaledArithmeticStatus status) {
    switch (status) {
    case ScaledArithmeticStatus::Success:
        return ExpressionTaylorJetStatus::Success;
    case ScaledArithmeticStatus::Nonfinite:
        return ExpressionTaylorJetStatus::Nonfinite;
    case ScaledArithmeticStatus::Singular:
        return ExpressionTaylorJetStatus::UnsupportedProgram;
    case ScaledArithmeticStatus::ExponentRange:
        return ExpressionTaylorJetStatus::ExponentRange;
    }
    return ExpressionTaylorJetStatus::ExponentRange;
}

bool checkedAddSize(
        size_t& total, size_t count, size_t bytes) {
    if (bytes != 0 &&
        count >
            (std::numeric_limits<size_t>::max() - total) /
                bytes)
        return false;
    total += count * bytes;
    return true;
}

bool checkedShiftExponent(
        ScaledRealValue& value, int64_t shift) {
    if (value.isZero()) return true;
    if ((shift > 0 &&
         value.exponent >
             std::numeric_limits<int64_t>::max() - shift) ||
        (shift < 0 &&
         value.exponent <
             std::numeric_limits<int64_t>::min() - shift))
        return false;
    value.exponent += shift;
    return true;
}

ScaledArithmeticStatus upperMagnitude(
        const ScaledComplexBall& value,
        ScaledRealValue& output) {
    if (!validRadius(value.radius))
        return ScaledArithmeticStatus::ExponentRange;
    const ScaledRealValue real = absoluteValue(value.value.re);
    const ScaledRealValue imaginary =
        absoluteValue(value.value.im);
    ScaledArithmeticStatus status =
        scaledAddUp(real, imaginary, output);
    if (status != ScaledArithmeticStatus::Success)
        return status;
    ScaledRealValue sum;
    status = scaledAddUp(output, value.radius, sum);
    if (status != ScaledArithmeticStatus::Success)
        return status;
    output = sum;
    status = scaledAddUp(output, value.radius, sum);
    if (status == ScaledArithmeticStatus::Success)
        output = sum;
    return status;
}

ScaledArithmeticStatus upperMagnitude(
        const ScaledComplexValue& value,
        ScaledRealValue& output) {
    const ScaledRealValue real = absoluteValue(value.re);
    const ScaledRealValue imaginary =
        absoluteValue(value.im);
    return scaledAddUp(real, imaginary, output);
}

ScaledArithmeticStatus upperMaxComponent(
        const ScaledComplexBall& value,
        ScaledRealValue& output) {
    if (!validRadius(value.radius))
        return ScaledArithmeticStatus::ExponentRange;
    const ScaledRealValue real =
        absoluteValue(value.value.re);
    const ScaledRealValue imaginary =
        absoluteValue(value.value.im);
    output =
        compareScaledNonnegative(real, imaginary) >= 0
        ? real : imaginary;
    ScaledRealValue sum;
    const ScaledArithmeticStatus status =
        scaledAddUp(output, value.radius, sum);
    if (status == ScaledArithmeticStatus::Success)
        output = sum;
    return status;
}

ScaledArithmeticStatus addUp(
        ScaledRealValue& target,
        const ScaledRealValue& value) {
    ScaledRealValue sum;
    const ScaledArithmeticStatus status =
        scaledAddUp(target, value, sum);
    if (status == ScaledArithmeticStatus::Success)
        target = sum;
    return status;
}

ScaledArithmeticStatus multiplyUp(
        const ScaledRealValue& left,
        const ScaledRealValue& right,
        ScaledRealValue& output) {
    return scaledMultiplyUp(left, right, output);
}

ScaledArithmeticStatus makeScaledNonnegativeDownward(
        mpfr_srcptr value, ScaledRealValue& output) {
    if (!mpfr_number_p(value) || mpfr_sgn(value) < 0)
        return ScaledArithmeticStatus::Nonfinite;
    if (mpfr_zero_p(value)) {
        output = {};
        return ScaledArithmeticStatus::Success;
    }
    long exponent = 0;
    double mantissa =
        mpfr_get_d_2exp(&exponent, value, MPFR_RNDZ);
    if (!(mantissa > 0.0) || !std::isfinite(mantissa))
        return ScaledArithmeticStatus::ExponentRange;
    if (mantissa >= 1.0) {
        mantissa *= 0.5;
        if (exponent == LONG_MAX)
            return ScaledArithmeticStatus::ExponentRange;
        ++exponent;
    }
    output.mantissa = mantissa;
    output.exponent = static_cast<int64_t>(exponent);
    return output.isNormalized()
        ? ScaledArithmeticStatus::Success
        : ScaledArithmeticStatus::ExponentRange;
}

enum class PrincipalFunction : uint8_t {
    Log,
    Log10,
    Sqrt,
    Exp
};

enum class PrincipalClearanceFailure : uint8_t {
    None,
    Zero,
    Cut
};

ScaledArithmeticStatus encodeMpfrComplexBall(
        const MpfrComplex& value,
        mpfr_prec_t precision,
        ScaledComplexBall& output) {
    output = {};
    if (!mpfr_number_p(value.re) ||
        !mpfr_number_p(value.im) ||
        precision < MPFR_PREC_MIN)
        return ScaledArithmeticStatus::Nonfinite;
    ScaledArithmeticStatus status =
        makeScaledComplexValue(value, output.value);
    if (status != ScaledArithmeticStatus::Success)
        return status;

    MpfrComplex encoded(precision);
    mpfr_t firstError, secondError, radius, scale;
    mpfr_inits2(
        precision, firstError, secondError, radius, scale,
        (mpfr_ptr)0);
    const mpfr_flags_t saved = mpfr_flags_save();
    mpfr_flags_clear(MPFR_FLAGS_ALL);
    bool okay =
        setMpfrFromScaledValue(
            encoded, output.value, MPFR_RNDN);
    if (okay) {
        mpfr_set_zero(radius, 0);
        mpfr_sub(
            firstError, value.re, encoded.re, MPFR_RNDD);
        mpfr_sub(
            secondError, value.re, encoded.re, MPFR_RNDU);
        mpfr_abs(firstError, firstError, MPFR_RNDU);
        mpfr_abs(secondError, secondError, MPFR_RNDU);
        mpfr_max(firstError, firstError, secondError, MPFR_RNDU);
        mpfr_max(radius, radius, firstError, MPFR_RNDU);
        mpfr_sub(
            firstError, value.im, encoded.im, MPFR_RNDD);
        mpfr_sub(
            secondError, value.im, encoded.im, MPFR_RNDU);
        mpfr_abs(firstError, firstError, MPFR_RNDU);
        mpfr_abs(secondError, secondError, MPFR_RNDU);
        mpfr_max(firstError, firstError, secondError, MPFR_RNDU);
        mpfr_max(radius, radius, firstError, MPFR_RNDU);

        mpfr_set_ui(scale, 1, MPFR_RNDU);
        mpfr_abs(firstError, value.re, MPFR_RNDU);
        mpfr_max(scale, scale, firstError, MPFR_RNDU);
        mpfr_abs(firstError, value.im, MPFR_RNDU);
        mpfr_max(scale, scale, firstError, MPFR_RNDU);
        // Every local principal value uses a short stable sequence of
        // correctly rounded MPFR operations. Eight guard bits at the working
        // precision cover that sequence in addition to encoding error.
        if (precision >
                static_cast<mpfr_prec_t>(LONG_MAX - 8)) {
            okay = false;
        } else {
            mpfr_mul_2si(
                scale, scale,
                8 - static_cast<long>(precision),
                MPFR_RNDU);
            mpfr_max(radius, radius, scale, MPFR_RNDU);
            okay =
                mpfr_number_p(radius) &&
                mpfr_sgn(radius) >= 0 &&
                mpfr_flags_test(
                    MPFR_FLAGS_UNDERFLOW |
                    MPFR_FLAGS_OVERFLOW |
                    MPFR_FLAGS_ERANGE) == 0;
        }
    }
    if (okay)
        status = makeScaledNonnegativeUpward(
            radius, output.radius);
    else
        status = ScaledArithmeticStatus::ExponentRange;
    mpfr_flags_restore(saved, MPFR_FLAGS_ALL);
    mpfr_clears(
        firstError, secondError, radius, scale,
        (mpfr_ptr)0);
    if (status != ScaledArithmeticStatus::Success)
        return status;
    return certifyScaledMpfrExponentRange(output);
}

ScaledArithmeticStatus makePrincipalFunctionBall(
        const ScaledComplexValue& input,
        PrincipalFunction function,
        mpfr_prec_t precision,
        ScaledComplexBall& output) {
    output = {};
    if (precision < 256)
        precision = 256;
    if (precision > MPFR_PREC_MAX - 64)
        return ScaledArithmeticStatus::ExponentRange;
    precision += 64;

    MpfrComplex argument(precision);
    MpfrComplex value(precision);
    mpfr_t magnitude, first, second, third;
    mpfr_inits2(
        precision, magnitude, first, second, third,
        (mpfr_ptr)0);
    const mpfr_flags_t saved = mpfr_flags_save();
    mpfr_flags_clear(MPFR_FLAGS_ALL);
    bool okay =
        setMpfrFromScaledValue(
            argument, input, MPFR_RNDN);
    if (okay) {
        switch (function) {
        case PrincipalFunction::Log:
        case PrincipalFunction::Log10:
            mpfr_hypot(
                magnitude, argument.re, argument.im,
                MPFR_RNDN);
            okay = !mpfr_zero_p(magnitude);
            if (okay) {
                mpfr_log(value.re, magnitude, MPFR_RNDN);
                mpfr_atan2(
                    value.im, argument.im, argument.re,
                    MPFR_RNDN);
                if (function == PrincipalFunction::Log10) {
                    mpfr_const_log2(first, MPFR_RNDN);
                    mpfr_log_ui(second, 5, MPFR_RNDN);
                    mpfr_add(first, first, second, MPFR_RNDN);
                    mpfr_div(
                        value.re, value.re, first,
                        MPFR_RNDN);
                    mpfr_div(
                        value.im, value.im, first,
                        MPFR_RNDN);
                }
            }
            break;
        case PrincipalFunction::Sqrt:
            okay = !(mpfr_zero_p(argument.re) &&
                     mpfr_zero_p(argument.im));
            if (okay) {
                mpfr_hypot(
                    magnitude, argument.re, argument.im,
                    MPFR_RNDN);
                if (mpfr_sgn(argument.re) >= 0) {
                    mpfr_add(
                        first, magnitude, argument.re,
                        MPFR_RNDN);
                    mpfr_div_2ui(
                        first, first, 1, MPFR_RNDN);
                    mpfr_sqrt(value.re, first, MPFR_RNDN);
                    mpfr_mul_2ui(
                        second, value.re, 1, MPFR_RNDN);
                    if (mpfr_zero_p(second)) {
                        okay = false;
                    } else {
                        mpfr_div(
                            value.im, argument.im, second,
                            MPFR_RNDN);
                    }
                } else {
                    mpfr_sub(
                        first, magnitude, argument.re,
                        MPFR_RNDN);
                    mpfr_div_2ui(
                        first, first, 1, MPFR_RNDN);
                    mpfr_sqrt(first, first, MPFR_RNDN);
                    mpfr_copysign(
                        value.im, first, argument.im,
                        MPFR_RNDN);
                    mpfr_abs(second, argument.im, MPFR_RNDN);
                    mpfr_abs(third, value.im, MPFR_RNDN);
                    mpfr_mul_2ui(
                        third, third, 1, MPFR_RNDN);
                    if (mpfr_zero_p(third)) {
                        okay = false;
                    } else {
                        mpfr_div(
                            value.re, second, third,
                            MPFR_RNDN);
                    }
                }
            }
            break;
        case PrincipalFunction::Exp:
            mpfr_exp(magnitude, argument.re, MPFR_RNDN);
            mpfr_cos(first, argument.im, MPFR_RNDN);
            mpfr_sin(second, argument.im, MPFR_RNDN);
            mpfr_mul(
                value.re, magnitude, first, MPFR_RNDN);
            mpfr_mul(
                value.im, magnitude, second, MPFR_RNDN);
            break;
        }
    }
    okay =
        okay &&
        mpfr_number_p(value.re) &&
        mpfr_number_p(value.im) &&
        mpfr_flags_test(
            MPFR_FLAGS_UNDERFLOW |
            MPFR_FLAGS_OVERFLOW |
            MPFR_FLAGS_ERANGE) == 0;
    ScaledArithmeticStatus status =
        okay
        ? encodeMpfrComplexBall(value, precision, output)
        : ScaledArithmeticStatus::ExponentRange;
    mpfr_flags_restore(saved, MPFR_FLAGS_ALL);
    mpfr_clears(
        magnitude, first, second, third,
        (mpfr_ptr)0);
    return status;
}

ScaledArithmeticStatus makeReciprocalLn10Ball(
        mpfr_prec_t precision,
        ScaledComplexBall& output) {
    if (precision < 256)
        precision = 256;
    if (precision > MPFR_PREC_MAX - 64)
        return ScaledArithmeticStatus::ExponentRange;
    precision += 64;
    MpfrComplex value(precision);
    mpfr_t logarithmFive;
    mpfr_init2(logarithmFive, precision);
    const mpfr_flags_t saved = mpfr_flags_save();
    mpfr_flags_clear(MPFR_FLAGS_ALL);
    mpfr_const_log2(value.re, MPFR_RNDN);
    mpfr_log_ui(logarithmFive, 5, MPFR_RNDN);
    mpfr_add(
        value.re, value.re, logarithmFive, MPFR_RNDN);
    mpfr_ui_div(value.re, 1, value.re, MPFR_RNDN);
    mpfr_set_zero(value.im, 0);
    const bool okay =
        mpfr_number_p(value.re) &&
        mpfr_sgn(value.re) > 0 &&
        mpfr_flags_test(
            MPFR_FLAGS_UNDERFLOW |
            MPFR_FLAGS_OVERFLOW |
            MPFR_FLAGS_ERANGE) == 0;
    const ScaledArithmeticStatus status =
        okay
        ? encodeMpfrComplexBall(value, precision, output)
        : ScaledArithmeticStatus::ExponentRange;
    mpfr_flags_restore(saved, MPFR_FLAGS_ALL);
    mpfr_clear(logarithmFive);
    return status;
}

ScaledArithmeticStatus logarithmSeriesTailBound(
        const ScaledRealValue& radius,
        int retainedOrder,
        ScaledRealValue& output) {
    output = {};
    if (!validRadius(radius) || retainedOrder < 1)
        return ScaledArithmeticStatus::ExponentRange;
    if (radius.isZero())
        return ScaledArithmeticStatus::Success;
    mpfr_t r, oneMinus, power, denominator, bound;
    mpfr_inits2(
        256, r, oneMinus, power, denominator, bound,
        (mpfr_ptr)0);
    const mpfr_flags_t saved = mpfr_flags_save();
    mpfr_flags_clear(MPFR_FLAGS_ALL);
    bool okay =
        setMpfrFromScaledValue(r, radius, MPFR_RNDU);
    if (okay) {
        mpfr_ui_sub(oneMinus, 1, r, MPFR_RNDD);
        okay = mpfr_sgn(oneMinus) > 0;
    }
    if (okay) {
        const unsigned long next =
            static_cast<unsigned long>(retainedOrder + 1);
        mpfr_pow_ui(power, r, next, MPFR_RNDU);
        mpfr_mul_ui(
            denominator, oneMinus, next, MPFR_RNDD);
        mpfr_div(bound, power, denominator, MPFR_RNDU);
        okay =
            mpfr_number_p(bound) &&
            mpfr_sgn(bound) >= 0 &&
            mpfr_flags_test(
                MPFR_FLAGS_UNDERFLOW |
                MPFR_FLAGS_OVERFLOW |
                MPFR_FLAGS_ERANGE) == 0;
    }
    ScaledArithmeticStatus status =
        okay
        ? makeScaledNonnegativeUpward(bound, output)
        : ScaledArithmeticStatus::Singular;
    mpfr_flags_restore(saved, MPFR_FLAGS_ALL);
    mpfr_clears(
        r, oneMinus, power, denominator, bound,
        (mpfr_ptr)0);
    return status;
}

ScaledArithmeticStatus principalBranchClearance(
        const ScaledComplexValue& center,
        const ScaledRealValue& variation,
        ScaledRealValue& zeroClearance,
        ScaledRealValue& cutClearance,
        PrincipalClearanceFailure& failure) {
    zeroClearance = {};
    cutClearance = {};
    failure = PrincipalClearanceFailure::None;
    if (!validRadius(variation) ||
        certifyScaledMpfrExponentRange(center) !=
            ScaledArithmeticStatus::Success)
        return ScaledArithmeticStatus::ExponentRange;
    MpfrComplex value(256);
    mpfr_t uncertainty, real, imaginary;
    mpfr_t zeroDistance, cutDistance;
    mpfr_t zeroLower, cutLower;
    mpfr_inits2(
        256, uncertainty, real, imaginary,
        zeroDistance, cutDistance, zeroLower, cutLower,
        (mpfr_ptr)0);
    const mpfr_flags_t saved = mpfr_flags_save();
    mpfr_flags_clear(MPFR_FLAGS_ALL);
    bool okay =
        setMpfrFromScaledValue(value, center, MPFR_RNDN) &&
        setMpfrFromScaledValue(
            uncertainty, variation, MPFR_RNDU);
    if (okay) {
        mpfr_abs(real, value.re, MPFR_RNDD);
        mpfr_abs(imaginary, value.im, MPFR_RNDD);
        mpfr_max(
            zeroDistance, real, imaginary, MPFR_RNDD);
        // This is the exact L-infinity distance to the closed negative-real
        // ray: the endpoint when Re(center)>0, otherwise vertical distance.
        if (mpfr_sgn(value.re) > 0) {
            mpfr_max(
                cutDistance, value.re, imaginary,
                MPFR_RNDD);
        } else {
            mpfr_set(
                cutDistance, imaginary, MPFR_RNDD);
        }
        mpfr_sub(
            zeroLower, zeroDistance, uncertainty,
            MPFR_RNDD);
        mpfr_sub(
            cutLower, cutDistance, uncertainty,
            MPFR_RNDD);
        okay =
            mpfr_number_p(zeroLower) &&
            mpfr_number_p(cutLower) &&
            mpfr_flags_test(
                MPFR_FLAGS_UNDERFLOW |
                MPFR_FLAGS_OVERFLOW |
                MPFR_FLAGS_ERANGE) == 0;
    }
    ScaledArithmeticStatus status =
        ScaledArithmeticStatus::ExponentRange;
    if (okay) {
        status = ScaledArithmeticStatus::Success;
        const bool zeroPositive = mpfr_sgn(zeroLower) > 0;
        const bool cutPositive = mpfr_sgn(cutLower) > 0;
        if (zeroPositive) {
            status = makeScaledNonnegativeDownward(
                zeroLower, zeroClearance);
        } else {
            failure = PrincipalClearanceFailure::Zero;
        }
        if (status == ScaledArithmeticStatus::Success) {
            const ScaledArithmeticStatus cutStatus =
                cutPositive
                ? makeScaledNonnegativeDownward(
                      cutLower, cutClearance)
                : ScaledArithmeticStatus::Success;
            if (cutStatus != ScaledArithmeticStatus::Success)
                status = cutStatus;
        }
        if (!cutPositive &&
            failure == PrincipalClearanceFailure::None)
            failure = PrincipalClearanceFailure::Cut;
        if (status == ScaledArithmeticStatus::Success &&
            failure != PrincipalClearanceFailure::None)
            status = ScaledArithmeticStatus::Singular;
    }
    mpfr_flags_restore(saved, MPFR_FLAGS_ALL);
    mpfr_clears(
        uncertainty, real, imaginary,
        zeroDistance, cutDistance, zeroLower, cutLower,
        (mpfr_ptr)0);
    return status;
}

ScaledArithmeticStatus lowerMaxComponentAfterVariation(
        const ScaledComplexValue& center,
        const ScaledRealValue& variation,
        ScaledRealValue& output) {
    if (!validRadius(variation))
        return ScaledArithmeticStatus::ExponentRange;
    MpfrComplex value(256);
    mpfr_t radius, real, imaginary, magnitude, lower;
    mpfr_inits2(
        256, radius, real, imaginary, magnitude, lower,
        (mpfr_ptr)0);
    const mpfr_flags_t saved = mpfr_flags_save();
    mpfr_flags_clear(MPFR_FLAGS_ALL);
    bool reconstructed =
        setMpfrFromScaledValue(value, center, MPFR_RNDN) &&
        setMpfrFromScaledValue(
            radius, variation, MPFR_RNDU);
    bool positive = false;
    if (reconstructed) {
        mpfr_abs(real, value.re, MPFR_RNDD);
        mpfr_abs(imaginary, value.im, MPFR_RNDD);
        mpfr_max(
            magnitude, real, imaginary, MPFR_RNDD);
        mpfr_sub(lower, magnitude, radius, MPFR_RNDD);
        reconstructed =
            mpfr_number_p(lower) &&
            mpfr_flags_test(
                MPFR_FLAGS_UNDERFLOW |
                MPFR_FLAGS_OVERFLOW |
                MPFR_FLAGS_ERANGE) == 0;
        positive =
            reconstructed && mpfr_sgn(lower) > 0;
    }
    ScaledArithmeticStatus status =
        reconstructed
        ? ScaledArithmeticStatus::Singular
        : ScaledArithmeticStatus::ExponentRange;
    if (positive)
        status = makeScaledNonnegativeDownward(
            lower, output);
    mpfr_flags_restore(saved, MPFR_FLAGS_ALL);
    mpfr_clears(
        radius, real, imaginary, magnitude, lower,
        (mpfr_ptr)0);
    return status;
}

ScaledArithmeticStatus multiplyNonnegativeDown(
        const ScaledRealValue& left,
        const ScaledRealValue& right,
        ScaledRealValue& output) {
    if (!validRadius(left) || !validRadius(right))
        return ScaledArithmeticStatus::ExponentRange;
    if (left.isZero() || right.isZero()) {
        output = {};
        return ScaledArithmeticStatus::Success;
    }
    mpfr_t first, second, product;
    mpfr_inits2(256, first, second, product, (mpfr_ptr)0);
    const mpfr_flags_t saved = mpfr_flags_save();
    mpfr_flags_clear(MPFR_FLAGS_ALL);
    bool okay =
        setMpfrFromScaledValue(
            first, left, MPFR_RNDD) &&
        setMpfrFromScaledValue(
            second, right, MPFR_RNDD);
    if (okay) {
        mpfr_mul(product, first, second, MPFR_RNDD);
        okay =
            mpfr_number_p(product) &&
            mpfr_sgn(product) > 0 &&
            mpfr_flags_test(
                MPFR_FLAGS_UNDERFLOW |
                MPFR_FLAGS_OVERFLOW |
                MPFR_FLAGS_ERANGE) == 0;
    }
    ScaledArithmeticStatus status =
        ScaledArithmeticStatus::ExponentRange;
    if (okay)
        status = makeScaledNonnegativeDownward(
            product, output);
    mpfr_flags_restore(saved, MPFR_FLAGS_ALL);
    mpfr_clears(first, second, product, (mpfr_ptr)0);
    return status;
}

ScaledArithmeticStatus divideNonnegativeUp(
        const ScaledRealValue& numerator,
        const ScaledRealValue& denominator,
        ScaledRealValue& output) {
    if (!validRadius(numerator) ||
        !validRadius(denominator) ||
        denominator.isZero())
        return ScaledArithmeticStatus::Singular;
    mpfr_t left, right, quotient;
    mpfr_inits2(256, left, right, quotient, (mpfr_ptr)0);
    const mpfr_flags_t saved = mpfr_flags_save();
    mpfr_flags_clear(MPFR_FLAGS_ALL);
    bool okay =
        setMpfrFromScaledValue(
            left, numerator, MPFR_RNDU) &&
        setMpfrFromScaledValue(
            right, denominator, MPFR_RNDD);
    if (okay) {
        mpfr_div(quotient, left, right, MPFR_RNDU);
        okay =
            mpfr_number_p(quotient) &&
            mpfr_sgn(quotient) >= 0 &&
            mpfr_flags_test(
                MPFR_FLAGS_UNDERFLOW |
                MPFR_FLAGS_OVERFLOW |
                MPFR_FLAGS_ERANGE) == 0;
    }
    ScaledArithmeticStatus status =
        ScaledArithmeticStatus::ExponentRange;
    if (okay)
        status = makeScaledNonnegativeUpward(
            quotient, output);
    mpfr_flags_restore(saved, MPFR_FLAGS_ALL);
    mpfr_clears(left, right, quotient, (mpfr_ptr)0);
    return status;
}

ScaledArithmeticStatus reciprocalScaledValue(
        const ScaledComplexValue& input,
        ScaledComplexBall& output) {
    output = {};
    if (input.isZero())
        return ScaledArithmeticStatus::Singular;
    MpfrComplex value(256);
    if (!setMpfrFromScaledValue(
            value, input, MPFR_RNDN))
        return ScaledArithmeticStatus::ExponentRange;

    mpfr_t scaledReal, scaledImaginary;
    mpfr_t denominatorLow, denominatorHigh;
    mpfr_t numerator, lower, upper, midpoint;
    mpfr_t encoded, firstError, secondError, radius;
    mpfr_inits2(
        256, scaledReal, scaledImaginary,
        denominatorLow, denominatorHigh, numerator,
        lower, upper, midpoint, encoded,
        firstError, secondError, radius,
        (mpfr_ptr)0);
    const mpfr_flags_t saved = mpfr_flags_save();
    mpfr_flags_clear(MPFR_FLAGS_ALL);

    const mpfr_exp_t realExponent =
        mpfr_zero_p(value.re)
        ? mpfr_get_emin()
        : mpfr_get_exp(value.re);
    const mpfr_exp_t imaginaryExponent =
        mpfr_zero_p(value.im)
        ? mpfr_get_emin()
        : mpfr_get_exp(value.im);
    const mpfr_exp_t scaleExponent =
        std::max(realExponent, imaginaryExponent);
    bool okay =
        scaleExponent >= LONG_MIN &&
        scaleExponent <= LONG_MAX &&
        scaleExponent != LONG_MIN;
    if (okay) {
        mpfr_div_2si(
            scaledReal, value.re,
            static_cast<long>(scaleExponent),
            MPFR_RNDN);
        mpfr_div_2si(
            scaledImaginary, value.im,
            static_cast<long>(scaleExponent),
            MPFR_RNDN);
        mpfr_sqr(
            denominatorLow, scaledReal, MPFR_RNDD);
        mpfr_sqr(
            firstError, scaledImaginary, MPFR_RNDD);
        mpfr_add(
            denominatorLow, denominatorLow,
            firstError, MPFR_RNDD);
        mpfr_sqr(
            denominatorHigh, scaledReal, MPFR_RNDU);
        mpfr_sqr(
            firstError, scaledImaginary, MPFR_RNDU);
        mpfr_add(
            denominatorHigh, denominatorHigh,
            firstError, MPFR_RNDU);
        okay =
            mpfr_number_p(denominatorLow) &&
            mpfr_sgn(denominatorLow) > 0 &&
            mpfr_number_p(denominatorHigh);
    }

    ScaledRealValue componentRadius;
    auto makeComponent = [&](
            mpfr_srcptr source,
            ScaledRealValue& component) {
        mpfr_set(numerator, source, MPFR_RNDN);
        if (mpfr_sgn(numerator) >= 0) {
            mpfr_div(
                lower, numerator,
                denominatorHigh, MPFR_RNDD);
            mpfr_div(
                upper, numerator,
                denominatorLow, MPFR_RNDU);
        } else {
            mpfr_div(
                lower, numerator,
                denominatorLow, MPFR_RNDD);
            mpfr_div(
                upper, numerator,
                denominatorHigh, MPFR_RNDU);
        }
        mpfr_mul_2si(
            lower, lower,
            static_cast<long>(-scaleExponent),
            MPFR_RNDD);
        mpfr_mul_2si(
            upper, upper,
            static_cast<long>(-scaleExponent),
            MPFR_RNDU);
        mpfr_add(midpoint, lower, upper, MPFR_RNDN);
        mpfr_div_2ui(midpoint, midpoint, 1, MPFR_RNDN);
        ScaledArithmeticStatus status =
            makeScaledRealValue(midpoint, component);
        if (status != ScaledArithmeticStatus::Success)
            return status;
        if (!setMpfrFromScaledValue(
                encoded, component, MPFR_RNDN))
            return ScaledArithmeticStatus::ExponentRange;
        mpfr_sub(
            firstError, encoded, lower, MPFR_RNDU);
        mpfr_abs(
            firstError, firstError, MPFR_RNDU);
        mpfr_sub(
            secondError, upper, encoded, MPFR_RNDU);
        mpfr_abs(
            secondError, secondError, MPFR_RNDU);
        mpfr_max(radius, firstError, secondError, MPFR_RNDU);
        ScaledRealValue localRadius;
        status = makeScaledNonnegativeUpward(
            radius, localRadius);
        if (status != ScaledArithmeticStatus::Success)
            return status;
        if (compareScaledNonnegative(
                localRadius, componentRadius) > 0)
            componentRadius = localRadius;
        return ScaledArithmeticStatus::Success;
    };

    ScaledArithmeticStatus status =
        ScaledArithmeticStatus::ExponentRange;
    if (okay) {
        status = makeComponent(
            scaledReal, output.value.re);
        if (status == ScaledArithmeticStatus::Success) {
            mpfr_neg(
                scaledImaginary,
                scaledImaginary, MPFR_RNDN);
            status = makeComponent(
                scaledImaginary, output.value.im);
        }
    }
    if (status == ScaledArithmeticStatus::Success) {
        output.radius = componentRadius;
        if (mpfr_flags_test(
                MPFR_FLAGS_UNDERFLOW |
                MPFR_FLAGS_OVERFLOW |
                MPFR_FLAGS_ERANGE) != 0)
            status = ScaledArithmeticStatus::ExponentRange;
    }
    mpfr_flags_restore(saved, MPFR_FLAGS_ALL);
    mpfr_clears(
        scaledReal, scaledImaginary,
        denominatorLow, denominatorHigh, numerator,
        lower, upper, midpoint, encoded,
        firstError, secondError, radius,
        (mpfr_ptr)0);
    if (status != ScaledArithmeticStatus::Success)
        return status;
    return certifyScaledMpfrExponentRange(output);
}

ScaledArithmeticStatus finiteSeriesTailBound(
        const ScaledRealValue& radius,
        int retainedOrder,
        ScaledRealValue& output) {
    output = {};
    if (!validRadius(radius) || retainedOrder < 0)
        return ScaledArithmeticStatus::ExponentRange;
    if (radius.isZero())
        return ScaledArithmeticStatus::Success;
    if (certifyScaledMpfrExponentRange(radius) !=
            ScaledArithmeticStatus::Success)
        return ScaledArithmeticStatus::ExponentRange;

    mpfr_t r, exponential, power, factorial, bound;
    mpfr_inits2(
        256, r, exponential, power, factorial, bound,
        (mpfr_ptr)0);
    const mpfr_flags_t saved = mpfr_flags_save();
    mpfr_flags_clear(MPFR_FLAGS_ALL);
    bool okay =
        setMpfrFromScaledValue(r, radius, MPFR_RNDU);
    if (okay) {
        mpfr_exp(exponential, r, MPFR_RNDU);
        mpfr_pow_ui(
            power, r,
            static_cast<unsigned long>(retainedOrder + 1),
            MPFR_RNDU);
        mpfr_mul(bound, exponential, power, MPFR_RNDU);
        mpfr_fac_ui(
            factorial,
            static_cast<unsigned long>(retainedOrder + 1),
            MPFR_RNDN);
        mpfr_div(bound, bound, factorial, MPFR_RNDU);
        const mpfr_flags_t rangeFlags = mpfr_flags_test(
            MPFR_FLAGS_UNDERFLOW |
            MPFR_FLAGS_OVERFLOW |
            MPFR_FLAGS_ERANGE);
        okay =
            rangeFlags == 0 &&
            mpfr_number_p(bound) &&
            mpfr_sgn(bound) >= 0 &&
            !mpfr_zero_p(bound);
    }
    ScaledArithmeticStatus status =
        ScaledArithmeticStatus::ExponentRange;
    if (okay)
        status = makeScaledNonnegativeUpward(
            bound, output);
    mpfr_flags_restore(saved, MPFR_FLAGS_ALL);
    mpfr_clears(
        r, exponential, power, factorial, bound,
        (mpfr_ptr)0);
    if (status != ScaledArithmeticStatus::Success)
        return status;
    return certifyScaledMpfrExponentRange(output);
}

struct TaylorThreshold {
    ScaledRealValue midpoint;
    ScaledRealValue error;
};

bool makeThreshold(double bailout, TaylorThreshold& threshold) {
    ScaledComplexValue radius;
    return makeScaledRealValue(bailout, radius.re) ==
               ScaledArithmeticStatus::Success &&
           certifiedScaledNormSquared(
               radius, {}, threshold.midpoint,
               threshold.error) ==
               ScaledArithmeticStatus::Success;
}

enum class InsideStatus : uint8_t {
    Inside,
    Uncertain,
    Error
};

InsideStatus certifyInside(
        const ScaledComplexValue& base,
        const ScaledComplexBall* coefficients,
        size_t coefficientCount,
        const ScaledRealValue& remainder,
        const TaylorThreshold& threshold,
        ScaledRealValue& margin,
        ScaledArithmeticStatus& arithmetic) {
    if (!coefficients || coefficientCount == 0 ||
        !validRadius(remainder)) {
        arithmetic = ScaledArithmeticStatus::ExponentRange;
        return InsideStatus::Error;
    }
    ScaledComplexBall center;
    center.value = base;
    arithmetic = certifiedScaledAdd(
        center, coefficients[0], center);
    if (arithmetic != ScaledArithmeticStatus::Success)
        return InsideStatus::Error;
    ScaledRealValue tail = remainder;
    for (size_t coefficient = 1;
         coefficient < coefficientCount; ++coefficient) {
        ScaledRealValue magnitude;
        arithmetic = upperMagnitude(
            coefficients[coefficient], magnitude);
        if (arithmetic != ScaledArithmeticStatus::Success)
            return InsideStatus::Error;
        arithmetic = addUp(tail, magnitude);
        if (arithmetic != ScaledArithmeticStatus::Success)
            return InsideStatus::Error;
    }
    arithmetic = addUp(center.radius, tail);
    if (arithmetic != ScaledArithmeticStatus::Success)
        return InsideStatus::Error;
    arithmetic = certifyScaledMpfrExponentRange(center);
    if (arithmetic != ScaledArithmeticStatus::Success)
        return InsideStatus::Error;

    ScaledRealValue norm, normError;
    arithmetic = certifiedScaledNormSquared(
        center.value, center.radius, norm, normError);
    if (arithmetic != ScaledArithmeticStatus::Success)
        return InsideStatus::Error;
    ScaledRealValue normUpper;
    arithmetic = scaledAddUp(norm, normError, normUpper);
    if (arithmetic != ScaledArithmeticStatus::Success)
        return InsideStatus::Error;
    ScaledRealValue insideLimit;
    arithmetic = scaledAddUp(
        normUpper, threshold.error, insideLimit);
    if (arithmetic != ScaledArithmeticStatus::Success)
        return InsideStatus::Error;
    if (compareScaledNonnegative(
            insideLimit, threshold.midpoint) > 0)
        return InsideStatus::Uncertain;
    arithmetic = scaledSubtract(
        threshold.midpoint, insideLimit, margin);
    if (arithmetic != ScaledArithmeticStatus::Success)
        return InsideStatus::Error;
    ScaledRealValue conservativeMargin;
    arithmetic = scaledDivideByDouble(
        margin, 2.0, conservativeMargin);
    if (arithmetic != ScaledArithmeticStatus::Success)
        return InsideStatus::Error;
    margin = conservativeMargin;
    return InsideStatus::Inside;
}

InsideStatus certifyInside(
        const ScaledComplexValue& base,
        const ScaledComplexBall* coefficients,
        int order,
        const ScaledRealValue& remainder,
        const TaylorThreshold& threshold,
        ScaledRealValue& margin,
        ScaledArithmeticStatus& arithmetic) {
    if (order < 0) {
        arithmetic = ScaledArithmeticStatus::ExponentRange;
        return InsideStatus::Error;
    }
    return certifyInside(
        base, coefficients,
        static_cast<size_t>(order) + 1,
        remainder, threshold, margin, arithmetic);
}

struct Candidate {
    ExpressionTaylorJetStatus status =
        ExpressionTaylorJetStatus::NoCoverage;
    std::string reason;
    int order = 0;
    int landing = 0;
    std::vector<ScaledComplexBall> coefficients;
    ScaledRealValue remainder;
    std::vector<ScaledRealValue> margins;
    uint64_t operations = 0;
    uint64_t bivariateConvolutionOperations = 0;
    uint64_t functionSeriesCount = 0;
    uint64_t functionSeriesOperations = 0;
    int maximumFunctionSeriesOrder = 0;
    ScaledRealValue maximumFunctionSeriesTail;
    uint64_t reciprocalCount = 0;
    uint64_t reciprocalOperations = 0;
    int maximumReciprocalOrder = 0;
    ScaledRealValue minimumDenominatorClearance;
    ScaledRealValue maximumReciprocalTail;
    bool hasDenominatorClearance = false;
    bool poleRejected = false;
    uint64_t branchCompositionCount = 0;
    uint64_t branchCompositionOperations = 0;
    int maximumBranchSeriesOrder = 0;
    ScaledRealValue maximumBranchSeriesTail;
    ScaledRealValue minimumBranchCutClearance;
    ScaledRealValue minimumBranchZeroClearance;
    bool hasBranchCutClearance = false;
    bool hasBranchZeroClearance = false;
    bool branchRejected = false;
    uint64_t absBranchCount = 0;
    uint64_t absPositiveCellCount = 0;
    uint64_t absNegativeCellCount = 0;
    ScaledRealValue minimumFoldClearance;
    bool hasFoldClearance = false;
    bool foldRejected = false;
    int foldRejectionIteration = -1;
    std::string foldRejectionReason;
    bool landingUsesSampleOutput = false;
    size_t memoryBytes = 0;
};

} // namespace

const char* expressionTaylorJetStatusName(
        ExpressionTaylorJetStatus status) {
    switch (status) {
    case ExpressionTaylorJetStatus::Success:
        return "success";
    case ExpressionTaylorJetStatus::NoCoverage:
        return "no-coverage";
    case ExpressionTaylorJetStatus::InvalidRequest:
        return "invalid-request";
    case ExpressionTaylorJetStatus::UnsupportedProgram:
        return "unsupported-program";
    case ExpressionTaylorJetStatus::InvalidTape:
        return "invalid-tape";
    case ExpressionTaylorJetStatus::ResourceLimit:
        return "resource-limit";
    case ExpressionTaylorJetStatus::Cancelled:
        return "cancelled";
    case ExpressionTaylorJetStatus::ExponentRange:
        return "exponent-range";
    case ExpressionTaylorJetStatus::Nonfinite:
        return "nonfinite";
    case ExpressionTaylorJetStatus::AccuracyBudget:
        return "accuracy-budget";
    case ExpressionTaylorJetStatus::BailoutUncertain:
        return "bailout-uncertain";
    case ExpressionTaylorJetStatus::PoleRejected:
        return "pole-rejected";
    case ExpressionTaylorJetStatus::BranchRejected:
        return "branch-rejected";
    }
    return "invalid-request";
}

const char* expressionTaylorJetLayoutName(
        ExpressionTaylorJetLayout layout) {
    switch (layout) {
    case ExpressionTaylorJetLayout::ComplexUnivariate:
        return "complex-univariate";
    case ExpressionTaylorJetLayout::RealBivariate:
        return "real-bivariate";
    }
    return "complex-univariate";
}

namespace {

constexpr size_t MaximumBivariateMonomials =
    (ExpressionTaylorMaximumBivariateOrder + 1) *
    (ExpressionTaylorMaximumBivariateOrder + 2) / 2;

bool triangularStart(size_t degree, size_t& start) {
    size_t left = degree;
    size_t right = degree + 1;
    if (right == 0)
        return false;
    if ((left & 1) == 0)
        left /= 2;
    else
        right /= 2;
    if (right != 0 &&
        left > std::numeric_limits<size_t>::max() / right)
        return false;
    start = left * right;
    return true;
}

} // namespace

bool expressionTaylorBivariateMonomialCount(
        int order, size_t& count) {
    count = 0;
    if (order < 0)
        return false;
    const size_t degree = static_cast<size_t>(order);
    size_t start = 0;
    if (!triangularStart(degree + 1, start))
        return false;
    count = start;
    return true;
}

bool expressionTaylorBivariateIndex(
        int order, int qDegree, int conjugateDegree,
        size_t& index) {
    index = 0;
    if (order < 0 || qDegree < 0 ||
        conjugateDegree < 0 ||
        qDegree > order - conjugateDegree)
        return false;
    size_t count = 0;
    if (!expressionTaylorBivariateMonomialCount(
            order, count))
        return false;
    const size_t degree =
        static_cast<size_t>(qDegree) +
        static_cast<size_t>(conjugateDegree);
    size_t start = 0;
    if (!triangularStart(degree, start) ||
        static_cast<size_t>(conjugateDegree) >
            std::numeric_limits<size_t>::max() - start)
        return false;
    index =
        start + static_cast<size_t>(conjugateDegree);
    return index < count;
}

bool expressionTaylorBivariateExponents(
        int order, size_t index, int& qDegree,
        int& conjugateDegree) {
    qDegree = 0;
    conjugateDegree = 0;
    size_t count = 0;
    if (!expressionTaylorBivariateMonomialCount(
            order, count) ||
        index >= count)
        return false;
    size_t low = 0;
    size_t high = static_cast<size_t>(order) + 1;
    while (low + 1 < high) {
        const size_t middle = low + (high - low) / 2;
        size_t start = 0;
        if (!triangularStart(middle, start))
            return false;
        if (start <= index)
            low = middle;
        else
            high = middle;
    }
    size_t start = 0;
    if (!triangularStart(low, start))
        return false;
    const size_t conjugate = index - start;
    if (conjugate > low ||
        low > static_cast<size_t>(
            std::numeric_limits<int>::max()))
        return false;
    conjugateDegree = static_cast<int>(conjugate);
    qDegree = static_cast<int>(low - conjugate);
    return true;
}

bool makeExpressionTaylorFrameScale(
        const ScaledRealValue& maximumRealMagnitude,
        const ScaledRealValue& maximumImaginaryMagnitude,
        const ScaledRealValue& maximumComponentError,
        ScaledComplexValue& parameterScale) {
    if (!validRadius(maximumRealMagnitude) ||
        !validRadius(maximumImaginaryMagnitude) ||
        !validRadius(maximumComponentError))
        return false;
    ScaledRealValue bound;
    if (scaledAddUp(
            maximumRealMagnitude,
            maximumImaginaryMagnitude,
            bound) != ScaledArithmeticStatus::Success ||
        addUp(bound, maximumComponentError) !=
            ScaledArithmeticStatus::Success ||
        addUp(bound, maximumComponentError) !=
            ScaledArithmeticStatus::Success ||
        bound.isZero() ||
        bound.exponent ==
            std::numeric_limits<int64_t>::max())
        return false;
    parameterScale = {};
    parameterScale.re.mantissa = 0.5;
    parameterScale.re.exponent = bound.exponent + 1;
    return certifyScaledMpfrExponentRange(parameterScale) ==
           ScaledArithmeticStatus::Success;
}

bool makeExpressionTaylorNormalizedQ(
        const ScaledComplexBall& offset,
        const ScaledComplexValue& parameterScale,
        ScaledComplexBall& q) {
    if (!validRadius(offset.radius) ||
        parameterScale.re.mantissa != 0.5 ||
        parameterScale.im.isZero() == false ||
        !parameterScale.re.isNormalized())
        return false;
    const int64_t exponent = parameterScale.re.exponent;
    int64_t shift = 0;
    if (exponent >= 1) {
        shift = -(exponent - 1);
    } else {
        if (exponent <
                std::numeric_limits<int64_t>::min() + 2)
            return false;
        shift = 1 - exponent;
    }
    q = offset;
    if (!checkedShiftExponent(q.value.re, shift) ||
        !checkedShiftExponent(q.value.im, shift) ||
        !checkedShiftExponent(q.radius, shift))
        return false;
    return certifyScaledMpfrExponentRange(q) ==
           ScaledArithmeticStatus::Success;
}

bool expressionTaylorQInsideUnitDisk(
        const ScaledComplexBall& q) {
    if (certifyScaledMpfrExponentRange(q) !=
            ScaledArithmeticStatus::Success)
        return false;
    ScaledRealValue bound;
    if (scaledAddUp(
            absoluteValue(q.value.re),
            absoluteValue(q.value.im),
            bound) != ScaledArithmeticStatus::Success ||
        addUp(bound, q.radius) !=
            ScaledArithmeticStatus::Success ||
        addUp(bound, q.radius) !=
            ScaledArithmeticStatus::Success)
        return false;
    ScaledRealValue one;
    one.mantissa = 0.5;
    one.exponent = 1;
    return compareScaledNonnegative(bound, one) <= 0;
}

class ExpressionRealTaylorJetBuilder {
public:
    static bool build(
            const ExpressionTaylorJetRequest& request,
            ExpressionTaylorJetResult& result) {
        const Clock::time_point start = Clock::now();
        result = {};
        result.layout =
            ExpressionTaylorJetLayout::RealBivariate;
        auto finish = [&](ExpressionTaylorJetStatus status,
                          const std::string& reason,
                          bool success) {
            result.status = status;
            result.failureReason = reason;
            result.buildSeconds =
                std::chrono::duration<double>(
                    Clock::now() - start).count();
            return success;
        };

        if (!request.program || !request.reference ||
            !request.program->valid() ||
            !request.reference->valid ||
            request.reference->status !=
                ExpressionReferenceBuildStatus::Success ||
            !request.reference->
                certifiedAgainstHigherPrecision ||
            (request.program->scaledResidualCapability() !=
                 ExpressionScaledResidualCapability::
                     CertifiedRealCandidate &&
             request.program->scaledResidualCapability() !=
                 ExpressionScaledResidualCapability::
                     CertifiedPiecewiseCandidate) ||
            (request.pixelParameter != FormulaParameter::C &&
             request.pixelParameter !=
                 FormulaParameter::InitialZ) ||
            request.minimumOrder < 8 ||
            request.maximumOrder > 20 ||
            request.maximumBivariateOrder < 8 ||
            request.maximumBivariateOrder >
                ExpressionTaylorMaximumBivariateOrder ||
            request.minimumOrder > request.preferredOrder ||
            request.preferredOrder > request.maximumOrder ||
            request.minimumOrder >
                request.maximumBivariateOrder ||
            request.minimumLanding < 1 ||
            !(request.bailout > 0.0) ||
            !std::isfinite(request.bailout) ||
            !(request.accuracyBudget > 0.0) ||
            !std::isfinite(request.accuracyBudget) ||
            !request.parameterScale.isNormalized() ||
            request.parameterScale.isZero())
            return finish(
                ExpressionTaylorJetStatus::InvalidRequest,
                "real-bivariate Taylor jet request is invalid",
                false);
        if (request.reference->programSemanticHash !=
                request.program->semanticHash() ||
            request.reference->programSource !=
                request.program->source())
            return finish(
                ExpressionTaylorJetStatus::InvalidTape,
                "real-bivariate Taylor reference semantic identity mismatch",
                false);
        if (request.reference->samples.size() < 2)
            return finish(
                ExpressionTaylorJetStatus::NoCoverage,
                "real-bivariate Taylor reference has no skippable prefix",
                false);
        if (certifyScaledMpfrExponentRange(
                request.parameterScale) !=
                ScaledArithmeticStatus::Success)
            return finish(
                ExpressionTaylorJetStatus::ExponentRange,
                "real-bivariate Taylor parameter scale is outside MPFR guards",
                false);

        using Op = ExpressionProgram::Op;
        for (const ExpressionProgram::Instruction& instruction :
             request.program->_code) {
            switch (instruction.op) {
            case Op::Constant:
            case Op::Z:
            case Op::C:
            case Op::Z0:
            case Op::Iteration:
            case Op::Parameter:
            case Op::Negate:
            case Op::Add:
            case Op::Subtract:
            case Op::Multiply:
            case Op::Square:
            case Op::Norm:
            case Op::Conjugate:
            case Op::Real:
            case Op::Imaginary:
            case Op::MakeComplex:
            case Op::Abs:
                break;
            default:
                return finish(
                    ExpressionTaylorJetStatus::
                        UnsupportedProgram,
                    "operation has no certified real-bivariate Taylor semantics",
                    false);
            }
        }

        size_t validationBytes = 0;
        if (!ExpressionScaledResidualEvaluator::
                estimateWorkspaceBytes(
                    *request.program, *request.reference,
                    validationBytes))
            return finish(
                ExpressionTaylorJetStatus::ResourceLimit,
                "real-bivariate tape validation workspace size overflow",
                false);
        if (request.memoryLimitBytes != 0 &&
            validationBytes > request.memoryLimitBytes)
            return finish(
                ExpressionTaylorJetStatus::ResourceLimit,
                "real-bivariate tape validation exceeds memory policy",
                false);
        ExpressionScaledResidualEvaluator validator;
        try {
            if (!validator.prepare(
                    *request.program, *request.reference))
                return finish(
                    ExpressionTaylorJetStatus::InvalidTape,
                    validator.error().empty()
                        ? "real-bivariate tape validation failed"
                        : validator.error(),
                    false);
        } catch (const std::bad_alloc&) {
            return finish(
                ExpressionTaylorJetStatus::ResourceLimit,
                "real-bivariate tape validation allocation failed",
                false);
        } catch (const std::length_error&) {
            return finish(
                ExpressionTaylorJetStatus::ResourceLimit,
                "real-bivariate tape validation length overflow",
                false);
        }
        if (validator.workspaceBytes() > validationBytes)
            return finish(
                ExpressionTaylorJetStatus::ResourceLimit,
                "real-bivariate tape validation exceeded its preflight bound",
                false);

        ScaledRealValue accuracyBudget;
        if (makeScaledRealValue(
                request.accuracyBudget, accuracyBudget) !=
                ScaledArithmeticStatus::Success)
            return finish(
                ExpressionTaylorJetStatus::InvalidRequest,
                "real-bivariate Taylor accuracy budget is not representable",
                false);
        TaylorThreshold threshold;
        if (!makeThreshold(request.bailout, threshold))
            return finish(
                ExpressionTaylorJetStatus::InvalidRequest,
                "real-bivariate Taylor bailout is not representable",
                false);

        const int maximumLanding = std::min<int>(
            request.maximumCandidateIteration > 0
                ? request.maximumCandidateIteration
                : static_cast<int>(
                      request.reference->samples.size()),
            static_cast<int>(
                request.reference->samples.size()));
        if (maximumLanding < request.minimumLanding)
            return finish(
                ExpressionTaylorJetStatus::NoCoverage,
                "real-bivariate reference is shorter than minimum landing",
                false);

        Candidate best;
        bool haveBest = false;
        uint64_t totalOperations = 0;
        uint64_t totalConvolutionOperations = 0;
        size_t peakMemory = 0;
        ExpressionTaylorJetStatus lastStatus =
            ExpressionTaylorJetStatus::NoCoverage;
        std::string lastReason =
            "real-bivariate Taylor prefix has no coverage";

        const int maximumOrder = std::min(
            request.maximumOrder,
            request.maximumBivariateOrder);
        for (int order = request.minimumOrder;
             order <= maximumOrder; ++order) {
            if (request.shouldCancel && request.shouldCancel())
                return finish(
                    ExpressionTaylorJetStatus::Cancelled,
                    "real-bivariate Taylor build cancelled",
                    false);

            size_t monomialCount = 0;
            if (!expressionTaylorBivariateMonomialCount(
                    order, monomialCount) ||
                monomialCount == 0 ||
                monomialCount >
                    MaximumBivariateMonomials) {
                lastStatus =
                    ExpressionTaylorJetStatus::ResourceLimit;
                lastReason =
                    "real-bivariate monomial count overflow";
                break;
            }
            const size_t nodeCount =
                request.program->instructionCount();
            if (nodeCount != 0 &&
                monomialCount >
                    std::numeric_limits<size_t>::max() /
                        nodeCount) {
                lastStatus =
                    ExpressionTaylorJetStatus::ResourceLimit;
                lastReason =
                    "real-bivariate node count overflow";
                break;
            }
            const size_t nodeEntries =
                nodeCount * monomialCount;
            size_t candidateBytes = validationBytes;
            if (haveBest &&
                (!checkedAddSize(
                     candidateBytes,
                     best.coefficients.capacity(),
                     sizeof(ScaledComplexBall)) ||
                 !checkedAddSize(
                     candidateBytes,
                     best.margins.capacity(),
                     sizeof(ScaledRealValue)))) {
                lastStatus =
                    ExpressionTaylorJetStatus::ResourceLimit;
                lastReason =
                    "retained real-bivariate candidate size overflow";
                break;
            }
            if (!checkedAddSize(
                    candidateBytes, nodeEntries,
                    sizeof(ScaledComplexBall)) ||
                !checkedAddSize(
                    candidateBytes, nodeCount,
                    sizeof(ScaledRealValue)) ||
                !checkedAddSize(
                    candidateBytes, nodeCount,
                    sizeof(uint8_t)) ||
                !checkedAddSize(
                    candidateBytes, nodeCount,
                    sizeof(uint8_t)) ||
                !checkedAddSize(
                    candidateBytes, monomialCount * 2,
                    sizeof(ScaledComplexBall)) ||
                !checkedAddSize(
                    candidateBytes,
                    static_cast<size_t>(maximumLanding) + 1,
                    sizeof(ScaledRealValue))) {
                lastStatus =
                    ExpressionTaylorJetStatus::ResourceLimit;
                lastReason =
                    "real-bivariate Taylor workspace size overflow";
                break;
            }
            peakMemory = std::max(
                peakMemory, candidateBytes);
            if (request.memoryLimitBytes != 0 &&
                candidateBytes >
                    request.memoryLimitBytes) {
                lastStatus =
                    ExpressionTaylorJetStatus::ResourceLimit;
                lastReason =
                    "real-bivariate Taylor workspace exceeds memory policy";
                break;
            }

            Candidate candidate;
            candidate.order = order;
            candidate.memoryBytes = candidateBytes;
            try {
                std::array<
                    uint8_t, MaximumBivariateMonomials>
                    qDegrees{};
                std::array<
                    uint8_t, MaximumBivariateMonomials>
                    conjugateDegrees{};
                std::array<
                    size_t, MaximumBivariateMonomials>
                    conjugateIndices{};
                for (size_t index = 0;
                     index < monomialCount; ++index) {
                    int qDegree = 0;
                    int conjugateDegree = 0;
                    size_t conjugateIndex = 0;
                    if (!expressionTaylorBivariateExponents(
                            order, index, qDegree,
                            conjugateDegree) ||
                        !expressionTaylorBivariateIndex(
                            order, conjugateDegree,
                            qDegree, conjugateIndex)) {
                        candidate.status =
                            ExpressionTaylorJetStatus::
                                ResourceLimit;
                        candidate.reason =
                            "real-bivariate index mapping failed";
                        break;
                    }
                    qDegrees[index] =
                        static_cast<uint8_t>(qDegree);
                    conjugateDegrees[index] =
                        static_cast<uint8_t>(
                            conjugateDegree);
                    conjugateIndices[index] =
                        conjugateIndex;
                }
                if (!candidate.reason.empty())
                    throw std::length_error(
                        candidate.reason);

                std::vector<ScaledComplexBall> nodes(
                    nodeEntries);
                std::vector<ScaledRealValue> remainders(
                    nodeCount);
                std::vector<uint8_t> exactReal(
                    nodeCount);
                std::vector<uint8_t> exactZero(
                    nodeCount);
                std::vector<ScaledComplexBall> state(
                    monomialCount);
                std::vector<ScaledComplexBall> nextState(
                    monomialCount);
                candidate.margins.reserve(
                    static_cast<size_t>(maximumLanding) + 1);
                size_t actualCandidateBytes =
                    validationBytes;
                if ((haveBest &&
                     (!checkedAddSize(
                          actualCandidateBytes,
                          best.coefficients.capacity(),
                          sizeof(ScaledComplexBall)) ||
                      !checkedAddSize(
                          actualCandidateBytes,
                          best.margins.capacity(),
                          sizeof(ScaledRealValue)))) ||
                    !checkedAddSize(
                        actualCandidateBytes,
                        nodes.capacity(),
                        sizeof(ScaledComplexBall)) ||
                    !checkedAddSize(
                        actualCandidateBytes,
                        remainders.capacity(),
                        sizeof(ScaledRealValue)) ||
                    !checkedAddSize(
                        actualCandidateBytes,
                        exactReal.capacity(),
                        sizeof(uint8_t)) ||
                    !checkedAddSize(
                        actualCandidateBytes,
                        exactZero.capacity(),
                        sizeof(uint8_t)) ||
                    !checkedAddSize(
                        actualCandidateBytes,
                        state.capacity(),
                        sizeof(ScaledComplexBall)) ||
                    !checkedAddSize(
                        actualCandidateBytes,
                        nextState.capacity(),
                        sizeof(ScaledComplexBall)) ||
                    !checkedAddSize(
                        actualCandidateBytes,
                        candidate.margins.capacity(),
                        sizeof(ScaledRealValue))) {
                    candidate.status =
                        ExpressionTaylorJetStatus::
                            ResourceLimit;
                    candidate.reason =
                        "real-bivariate Taylor vector capacity calculation overflow";
                    throw std::length_error(
                        candidate.reason);
                }
                candidate.memoryBytes =
                    actualCandidateBytes;
                peakMemory = std::max(
                    peakMemory, actualCandidateBytes);
                if (request.memoryLimitBytes != 0 &&
                    actualCandidateBytes >
                        request.memoryLimitBytes) {
                    candidate.status =
                        ExpressionTaylorJetStatus::
                            ResourceLimit;
                    candidate.reason =
                        "real-bivariate Taylor vector capacities exceed memory policy";
                    throw std::length_error(
                        candidate.reason);
                }

                auto count = [](uint64_t& target,
                                uint64_t amount = 1) {
                    target =
                        amount >
                            std::numeric_limits<uint64_t>::
                                max() - target
                        ? std::numeric_limits<uint64_t>::max()
                        : target + amount;
                };
                auto coefficient = [&](
                        size_t node, size_t index)
                        -> ScaledComplexBall& {
                    return nodes[
                        node * monomialCount + index];
                };
                auto nodeBase = [&](
                        size_t sampleOffset, uint16_t node,
                        ScaledComplexValue& value) {
                    const ExpressionReferenceTapeNode& tape =
                        request.reference->tape[
                            sampleOffset + node];
                    ScaledArithmeticStatus status =
                        makeScaledComplexValue(
                            tape.output,
                            tape.outputDefect, value);
                    if (status !=
                            ScaledArithmeticStatus::Success)
                        return status;
                    return certifyScaledMpfrExponentRange(
                        value);
                };
                auto addRemainder = [&](
                        ScaledRealValue& target,
                        const ScaledRealValue& value) {
                    count(candidate.operations);
                    return addUp(target, value);
                };
                auto addBall = [&](
                        const ScaledComplexBall& left,
                        const ScaledComplexBall& right,
                        ScaledComplexBall& output) {
                    count(candidate.operations);
                    return certifiedScaledAdd(
                        left, right, output);
                };
                auto subtractBall = [&](
                        const ScaledComplexBall& left,
                        const ScaledComplexBall& right,
                        ScaledComplexBall& output) {
                    count(candidate.operations);
                    return certifiedScaledSubtract(
                        left, right, output);
                };
                auto multiplyBall = [&](
                        const ScaledComplexBall& left,
                        const ScaledComplexBall& right,
                        ScaledComplexBall& output) {
                    count(candidate.operations);
                    return certifiedScaledMultiply(
                        left, right, output);
                };
                auto multiplyBound = [&](
                        const ScaledRealValue& left,
                        const ScaledRealValue& right,
                        ScaledRealValue& output) {
                    count(candidate.operations);
                    return multiplyUp(left, right, output);
                };
                auto conjugateBall = [&](
                        const ScaledComplexBall& input,
                        ScaledComplexBall& output) {
                    output.value.re = input.value.re;
                    output.radius = input.radius;
                    count(candidate.operations);
                    return scaledNegate(
                        input.value.im, output.value.im);
                };
                auto halfBall = [&](
                        const ScaledComplexBall& input,
                        ScaledComplexBall& output) {
                    count(candidate.operations, 2);
                    ScaledArithmeticStatus status =
                        scaledDivideByDouble(
                            input.value, 2.0,
                            output.value);
                    if (status ==
                            ScaledArithmeticStatus::Success)
                        status = scaledDivideByDouble(
                            input.radius, 2.0,
                            output.radius);
                    return status;
                };
                auto rotateMinusI = [&](
                        const ScaledComplexBall& input,
                        ScaledComplexBall& output) {
                    output.value.re = input.value.im;
                    output.radius = input.radius;
                    count(candidate.operations);
                    return scaledNegate(
                        input.value.re, output.value.im);
                };
                auto rotatePlusI = [&](
                        const ScaledComplexBall& input,
                        ScaledComplexBall& output) {
                    output.value.im = input.value.re;
                    output.radius = input.radius;
                    count(candidate.operations);
                    return scaledNegate(
                        input.value.im, output.value.re);
                };
                auto realCoefficient = [&](
                        uint16_t source, size_t index,
                        ScaledComplexBall& output) {
                    ScaledComplexBall swapped;
                    ScaledArithmeticStatus status =
                        conjugateBall(
                            coefficient(
                                source,
                                conjugateIndices[index]),
                            swapped);
                    ScaledComplexBall sum;
                    if (status ==
                            ScaledArithmeticStatus::Success)
                        status = addBall(
                            coefficient(source, index),
                            swapped, sum);
                    if (status ==
                            ScaledArithmeticStatus::Success)
                        status = halfBall(sum, output);
                    return status;
                };
                auto imaginaryCoefficient = [&](
                        uint16_t source, size_t index,
                        ScaledComplexBall& output) {
                    ScaledComplexBall swapped;
                    ScaledArithmeticStatus status =
                        conjugateBall(
                            coefficient(
                                source,
                                conjugateIndices[index]),
                            swapped);
                    ScaledComplexBall difference;
                    if (status ==
                            ScaledArithmeticStatus::Success)
                        status = subtractBall(
                            coefficient(source, index),
                            swapped, difference);
                    ScaledComplexBall rotated;
                    if (status ==
                            ScaledArithmeticStatus::Success)
                        status = rotateMinusI(
                            difference, rotated);
                    if (status ==
                            ScaledArithmeticStatus::Success)
                        status = halfBall(
                            rotated, output);
                    return status;
                };
                auto polynomialSup = [&](
                        uint16_t node,
                        ScaledRealValue& output) {
                    output = {};
                    for (size_t index = 0;
                         index < monomialCount; ++index) {
                        ScaledRealValue magnitude;
                        ScaledArithmeticStatus status =
                            upperMagnitude(
                                coefficient(node, index),
                                magnitude);
                        count(candidate.operations);
                        if (status !=
                                ScaledArithmeticStatus::
                                    Success)
                            return status;
                        status = addRemainder(
                            output, magnitude);
                        if (status !=
                                ScaledArithmeticStatus::
                                    Success)
                            return status;
                    }
                    return ScaledArithmeticStatus::Success;
                };
                auto certifyAbsBranch = [&](
                        size_t sampleOffset,
                        uint16_t source,
                        int& sign,
                        ScaledRealValue& clearance) {
                    sign = 0;
                    clearance = {};
                    if (!exactReal[source])
                        return ScaledArithmeticStatus::Singular;
                    const ExpressionReferenceTapeNode& tape =
                        request.reference->tape[
                            sampleOffset + source];
                    ScaledComplexBall primary;
                    ScaledComplexBall defect;
                    ScaledComplexBall center;
                    ScaledArithmeticStatus status =
                        makeScaledComplexValue(
                            tape.output, primary.value);
                    if (status ==
                            ScaledArithmeticStatus::Success)
                        status = makeScaledComplexValue(
                            tape.outputDefect,
                            defect.value);
                    if (status ==
                            ScaledArithmeticStatus::Success)
                        status = addBall(
                            primary, defect, center);
                    if (status ==
                            ScaledArithmeticStatus::Success)
                        status = addRemainder(
                            center.radius,
                            tape.outputError);
                    if (status !=
                            ScaledArithmeticStatus::Success)
                        return status;
                    if (!center.value.im.isZero())
                        return ScaledArithmeticStatus::Singular;

                    status = addBall(
                        center, coefficient(source, 0),
                        center);
                    if (status !=
                            ScaledArithmeticStatus::Success)
                        return status;
                    if (!center.value.im.isZero())
                        return ScaledArithmeticStatus::Singular;

                    ScaledRealValue radius =
                        remainders[source];
                    status = addRemainder(
                        radius, center.radius);
                    if (status !=
                            ScaledArithmeticStatus::Success)
                        return status;
                    for (size_t index = 1;
                         index < monomialCount; ++index) {
                        ScaledRealValue magnitude;
                        status = upperMagnitude(
                            coefficient(source, index),
                            magnitude);
                        count(candidate.operations);
                        if (status !=
                                ScaledArithmeticStatus::Success)
                            return status;
                        status = addRemainder(
                            radius, magnitude);
                        if (status !=
                                ScaledArithmeticStatus::Success)
                            return status;
                    }

                    const ScaledRealValue magnitude =
                        absoluteValue(center.value.re);
                    if (center.value.re.mantissa > 0.0 &&
                        compareScaledNonnegative(
                            magnitude, radius) > 0)
                        sign = 1;
                    else if (
                        center.value.re.mantissa < 0.0 &&
                        compareScaledNonnegative(
                            magnitude, radius) > 0)
                        sign = -1;
                    else
                        return ScaledArithmeticStatus::Success;

                    status = scaledSubtract(
                        magnitude, radius, clearance);
                    count(candidate.operations);
                    if (status !=
                            ScaledArithmeticStatus::Success)
                        return status;
                    ScaledRealValue conservative;
                    status = scaledDivideByDouble(
                        clearance, 2.0, conservative);
                    count(candidate.operations);
                    if (status ==
                            ScaledArithmeticStatus::Success)
                        clearance = conservative;
                    return status;
                };
                auto multiplyPolynomials = [&](
                        size_t sampleOffset,
                        uint16_t left, uint16_t right,
                        bool conjugateRight,
                        uint16_t outputNode) {
                    ScaledComplexValue leftBaseValue;
                    ScaledComplexValue rightBaseValue;
                    ScaledArithmeticStatus status =
                        nodeBase(
                            sampleOffset, left,
                            leftBaseValue);
                    if (status !=
                            ScaledArithmeticStatus::Success)
                        return status;
                    status = nodeBase(
                        sampleOffset, right,
                        rightBaseValue);
                    if (status !=
                            ScaledArithmeticStatus::Success)
                        return status;
                    ScaledComplexBall leftBase;
                    ScaledComplexBall rightBase;
                    leftBase.value = leftBaseValue;
                    rightBase.value = rightBaseValue;
                    if (conjugateRight) {
                        ScaledComplexBall conjugated;
                        status = conjugateBall(
                            rightBase, conjugated);
                        if (status !=
                                ScaledArithmeticStatus::
                                    Success)
                            return status;
                        rightBase = conjugated;
                        rightBaseValue = rightBase.value;
                    }

                    std::array<
                        ScaledRealValue,
                        MaximumBivariateMonomials>
                            leftMagnitudes{};
                    std::array<
                        ScaledRealValue,
                        MaximumBivariateMonomials>
                            rightMagnitudes{};
                    ScaledRealValue leftSup;
                    ScaledRealValue rightSup;
                    for (size_t index = 0;
                         index < monomialCount; ++index) {
                        ScaledComplexBall rightValue;
                        if (conjugateRight) {
                            status = conjugateBall(
                                coefficient(
                                    right,
                                    conjugateIndices[index]),
                                rightValue);
                        } else {
                            rightValue =
                                coefficient(right, index);
                        }
                        if (status !=
                                ScaledArithmeticStatus::
                                    Success)
                            return status;

                        ScaledComplexBall sum;
                        ScaledComplexBall term;
                        status = multiplyBall(
                            leftBase, rightValue, sum);
                        if (status ==
                                ScaledArithmeticStatus::
                                    Success)
                            status = multiplyBall(
                                rightBase,
                                coefficient(left, index),
                                term);
                        if (status ==
                                ScaledArithmeticStatus::
                                    Success)
                            status = addBall(
                                sum, term, sum);
                        if (status !=
                                ScaledArithmeticStatus::
                                    Success)
                            return status;
                        coefficient(outputNode, index) =
                            sum;

                        status = upperMagnitude(
                            coefficient(left, index),
                            leftMagnitudes[index]);
                        count(candidate.operations);
                        if (status ==
                                ScaledArithmeticStatus::
                                    Success)
                            status = upperMagnitude(
                                rightValue,
                                rightMagnitudes[index]);
                        count(candidate.operations);
                        if (status ==
                                ScaledArithmeticStatus::
                                    Success)
                            status = addRemainder(
                                leftSup,
                                leftMagnitudes[index]);
                        if (status ==
                                ScaledArithmeticStatus::
                                    Success)
                            status = addRemainder(
                                rightSup,
                                rightMagnitudes[index]);
                        if (status !=
                                ScaledArithmeticStatus::
                                    Success)
                            return status;
                    }

                    ScaledRealValue leftBaseMagnitude;
                    ScaledRealValue rightBaseMagnitude;
                    status = upperMagnitude(
                        leftBaseValue,
                        leftBaseMagnitude);
                    count(candidate.operations);
                    if (status ==
                            ScaledArithmeticStatus::Success)
                        status = upperMagnitude(
                            rightBaseValue,
                            rightBaseMagnitude);
                    count(candidate.operations);
                    if (status !=
                            ScaledArithmeticStatus::Success)
                        return status;

                    ScaledRealValue tail;
                    auto addProduct = [&](
                            const ScaledRealValue& first,
                            const ScaledRealValue& second,
                            bool doubleProduct = false,
                            bool convolution = false) {
                        ScaledRealValue product;
                        ScaledArithmeticStatus local =
                            multiplyBound(
                                first, second, product);
                        if (convolution)
                            count(
                                candidate.
                                    bivariateConvolutionOperations);
                        if (local !=
                                ScaledArithmeticStatus::
                                    Success)
                            return local;
                        if (doubleProduct) {
                            ScaledRealValue doubled;
                            local = scaledAddUp(
                                product, product,
                                doubled);
                            count(candidate.operations);
                            if (local !=
                                    ScaledArithmeticStatus::
                                        Success)
                                return local;
                            product = doubled;
                        }
                        local = addRemainder(
                            tail, product);
                        if (convolution)
                            count(
                                candidate.
                                    bivariateConvolutionOperations);
                        return local;
                    };
                    status = addProduct(
                        leftBaseMagnitude,
                        remainders[right]);
                    if (status ==
                            ScaledArithmeticStatus::Success)
                        status = addProduct(
                            rightBaseMagnitude,
                            remainders[left]);
                    if (status ==
                            ScaledArithmeticStatus::Success)
                        status = addProduct(
                            leftSup, remainders[right]);
                    if (status ==
                            ScaledArithmeticStatus::Success)
                        status = addProduct(
                            rightSup, remainders[left]);
                    if (status ==
                            ScaledArithmeticStatus::Success)
                        status = addProduct(
                            remainders[left],
                            remainders[right], true);
                    if (status !=
                            ScaledArithmeticStatus::Success)
                        return status;

                    for (size_t first = 0;
                         first < monomialCount; ++first) {
                        for (size_t second = 0;
                             second < monomialCount;
                             ++second) {
                            const int qDegree =
                                qDegrees[first] +
                                qDegrees[second];
                            const int conjugateDegree =
                                conjugateDegrees[first] +
                                conjugateDegrees[second];
                            if (qDegree + conjugateDegree >
                                    order) {
                                status = addProduct(
                                    leftMagnitudes[first],
                                    rightMagnitudes[second],
                                    false, true);
                                if (status !=
                                        ScaledArithmeticStatus::
                                            Success)
                                    return status;
                                continue;
                            }
                            size_t outputIndex = 0;
                            if (!expressionTaylorBivariateIndex(
                                    order, qDegree,
                                    conjugateDegree,
                                    outputIndex))
                                return
                                    ScaledArithmeticStatus::
                                        ExponentRange;
                            ScaledComplexBall rightValue;
                            if (conjugateRight) {
                                status = conjugateBall(
                                    coefficient(
                                        right,
                                        conjugateIndices[
                                            second]),
                                    rightValue);
                            } else {
                                rightValue =
                                    coefficient(
                                        right, second);
                            }
                            if (status !=
                                    ScaledArithmeticStatus::
                                        Success)
                                return status;
                            ScaledComplexBall term;
                            status = multiplyBall(
                                coefficient(left, first),
                                rightValue, term);
                            count(
                                candidate.
                                    bivariateConvolutionOperations);
                            if (status ==
                                    ScaledArithmeticStatus::
                                        Success)
                                status = addBall(
                                    coefficient(
                                        outputNode,
                                        outputIndex),
                                    term,
                                    coefficient(
                                        outputNode,
                                        outputIndex));
                            count(
                                candidate.
                                    bivariateConvolutionOperations);
                            if (status !=
                                    ScaledArithmeticStatus::
                                        Success)
                                return status;
                        }
                    }
                    remainders[outputNode] = tail;
                    return ScaledArithmeticStatus::Success;
                };
                auto addNodeRoundoff = [&](
                        const ExpressionReferenceTapeNode& node,
                        size_t sampleOffset,
                        ScaledRealValue& remainder) {
                    if (request.reference->
                            certificationPrecision <= 16)
                        return
                            ScaledArithmeticStatus::Success;
                    const bool arithmeticOperation =
                        node.operation ==
                            ExpressionOracleOperation::Add ||
                        node.operation ==
                            ExpressionOracleOperation::
                                Subtract ||
                        node.operation ==
                            ExpressionOracleOperation::
                                Multiply ||
                        node.operation ==
                            ExpressionOracleOperation::
                                Square ||
                        node.operation ==
                            ExpressionOracleOperation::Norm ||
                        node.operation ==
                            ExpressionOracleOperation::Abs;
                    if (!arithmeticOperation)
                        return
                            ScaledArithmeticStatus::Success;
                    auto fullBound = [&](
                            uint16_t child,
                            ScaledRealValue& output) {
                        ScaledComplexValue base;
                        ScaledArithmeticStatus localStatus =
                            nodeBase(
                                sampleOffset, child, base);
                        if (localStatus !=
                                ScaledArithmeticStatus::
                                    Success)
                            return localStatus;
                        ScaledRealValue baseMagnitude;
                        localStatus = upperMagnitude(
                            base, baseMagnitude);
                        count(candidate.operations);
                        if (localStatus !=
                                ScaledArithmeticStatus::
                                    Success)
                            return localStatus;
                        ScaledRealValue polynomial;
                        localStatus = polynomialSup(
                            child, polynomial);
                        if (localStatus !=
                                ScaledArithmeticStatus::
                                    Success)
                            return localStatus;
                        output = baseMagnitude;
                        localStatus = addRemainder(
                            output, polynomial);
                        if (localStatus ==
                                ScaledArithmeticStatus::
                                    Success)
                            localStatus = addRemainder(
                                output,
                                remainders[child]);
                        if (localStatus ==
                                ScaledArithmeticStatus::
                                    Success)
                            localStatus = addRemainder(
                                output,
                                remainders[child]);
                        return localStatus;
                    };

                    ScaledRealValue left;
                    ScaledRealValue right;
                    ScaledRealValue operationMagnitude;
                    ScaledArithmeticStatus status =
                        fullBound(node.leftNode, left);
                    if (status !=
                            ScaledArithmeticStatus::Success)
                        return status;
                    if (node.operation ==
                            ExpressionOracleOperation::
                                Square ||
                        node.operation ==
                            ExpressionOracleOperation::Norm) {
                        status = multiplyBound(
                            left, left,
                            operationMagnitude);
                    } else if (
                        node.operation ==
                            ExpressionOracleOperation::Abs) {
                        operationMagnitude = left;
                    } else {
                        status = fullBound(
                            node.rightNode, right);
                        if (status !=
                                ScaledArithmeticStatus::
                                    Success)
                            return status;
                        if (node.operation ==
                                ExpressionOracleOperation::
                                    Add ||
                            node.operation ==
                                ExpressionOracleOperation::
                                    Subtract) {
                            operationMagnitude = left;
                            status = addRemainder(
                                operationMagnitude, right);
                        } else {
                            status = multiplyBound(
                                left, right,
                                operationMagnitude);
                        }
                    }
                    if (status !=
                            ScaledArithmeticStatus::Success)
                        return status;
                    ScaledRealValue scale;
                    scale.mantissa = 0.5;
                    scale.exponent =
                        9 - static_cast<int64_t>(
                            request.reference->
                                certificationPrecision);
                    ScaledRealValue roundoff;
                    status = multiplyBound(
                        operationMagnitude, scale,
                        roundoff);
                    if (status !=
                            ScaledArithmeticStatus::Success)
                        return status;
                    return addRemainder(
                        remainder, roundoff);
                };

                size_t qIndex = 0;
                if (!expressionTaylorBivariateIndex(
                        order, 1, 0, qIndex)) {
                    candidate.status =
                        ExpressionTaylorJetStatus::
                            ResourceLimit;
                    candidate.reason =
                        "real-bivariate q index is unavailable";
                } else if (request.pixelParameter ==
                               FormulaParameter::InitialZ) {
                    state[qIndex].value =
                        request.parameterScale;
                }
                ScaledRealValue stateRemainder =
                    request.reference->samples[0].zError;
                if (candidate.reason.empty() &&
                    request.pixelParameter ==
                        FormulaParameter::InitialZ &&
                    addRemainder(
                        stateRemainder,
                        request.reference->z0Error) !=
                        ScaledArithmeticStatus::Success) {
                    candidate.status =
                        ExpressionTaylorJetStatus::
                            ExponentRange;
                    candidate.reason =
                        "initial real-bivariate Taylor radius overflow";
                }

                if (candidate.reason.empty()) {
                    ScaledComplexValue initialBase;
                    ScaledArithmeticStatus arithmetic =
                        makeScaledComplexValue(
                            request.reference->samples[0].z,
                            request.reference->samples[0].
                                zDefect,
                            initialBase);
                    ScaledRealValue margin;
                    const InsideStatus inside =
                        arithmetic ==
                                ScaledArithmeticStatus::
                                    Success
                        ? certifyInside(
                              initialBase, state.data(),
                              monomialCount, stateRemainder,
                              threshold, margin, arithmetic)
                        : InsideStatus::Error;
                    if (inside == InsideStatus::Inside) {
                        candidate.margins.push_back(
                            margin);
                    } else {
                        candidate.status =
                            inside ==
                                InsideStatus::Uncertain
                            ? ExpressionTaylorJetStatus::
                                  BailoutUncertain
                            : statusForArithmetic(arithmetic);
                        candidate.reason =
                            "initial real-bivariate Taylor frame overlaps bailout";
                    }
                }

                int landing = 0;
                while (candidate.reason.empty() &&
                       landing < maximumLanding) {
                    if (request.shouldCancel &&
                        request.shouldCancel()) {
                        candidate.status =
                            ExpressionTaylorJetStatus::
                                Cancelled;
                        candidate.reason =
                            "real-bivariate Taylor build cancelled";
                        break;
                    }
                    const ExpressionReferenceSample& sample =
                        request.reference->samples[
                            static_cast<size_t>(landing)];
                    const size_t sampleOffset =
                        static_cast<size_t>(
                            sample.tapeOffset);
                    std::fill(
                        nodes.begin(), nodes.end(),
                        ScaledComplexBall{});
                    std::fill(
                        remainders.begin(),
                        remainders.end(),
                        ScaledRealValue{});
                    std::fill(
                        exactReal.begin(),
                        exactReal.end(), uint8_t{ 0 });
                    std::fill(
                        exactZero.begin(),
                        exactZero.end(), uint8_t{ 0 });
                    ScaledArithmeticStatus arithmetic =
                        ScaledArithmeticStatus::Success;

                    for (size_t local = 0;
                         local < sample.tapeCount;
                         ++local) {
                        const ExpressionReferenceTapeNode&
                            node =
                                request.reference->tape[
                                    sampleOffset + local];
                        const Op op =
                            request.program->_code[local].
                                op;
                        auto copyNode = [&](
                                uint16_t source) {
                            for (size_t index = 0;
                                 index < monomialCount;
                                 ++index)
                                coefficient(local, index) =
                                    coefficient(
                                        source, index);
                            remainders[local] =
                                remainders[source];
                        };

                        switch (op) {
                        case Op::Constant:
                            exactReal[local] =
                                request.program->_code[local].
                                    value.imag() == 0.0;
                            exactZero[local] =
                                request.program->_code[local].
                                    value.real() == 0.0 &&
                                request.program->_code[local].
                                    value.imag() == 0.0;
                            break;
                        case Op::Iteration:
                            exactReal[local] = 1;
                            break;
                        case Op::Parameter:
                            break;
                        case Op::Z:
                            for (size_t index = 0;
                                 index < monomialCount;
                                 ++index)
                                coefficient(local, index) =
                                    state[index];
                            remainders[local] =
                                stateRemainder;
                            break;
                        case Op::C:
                            if (request.pixelParameter ==
                                    FormulaParameter::C)
                                coefficient(
                                    local, qIndex).value =
                                        request.
                                            parameterScale;
                            remainders[local] =
                                request.reference->cError;
                            break;
                        case Op::Z0:
                            if (request.pixelParameter ==
                                    FormulaParameter::
                                        InitialZ)
                                coefficient(
                                    local, qIndex).value =
                                        request.
                                            parameterScale;
                            remainders[local] =
                                request.reference->z0Error;
                            break;
                        case Op::Negate:
                            copyNode(node.leftNode);
                            exactReal[local] =
                                exactReal[node.leftNode];
                            exactZero[local] =
                                exactZero[node.leftNode];
                            for (size_t index = 0;
                                 index < monomialCount;
                                 ++index) {
                                arithmetic = scaledNegate(
                                    coefficient(
                                        node.leftNode,
                                        index).value,
                                    coefficient(
                                        local,
                                        index).value);
                                count(candidate.operations);
                                if (arithmetic !=
                                        ScaledArithmeticStatus::
                                            Success)
                                    break;
                            }
                            break;
                        case Op::Add:
                        case Op::Subtract:
                            exactReal[local] =
                                exactReal[node.leftNode] &&
                                exactReal[node.rightNode];
                            exactZero[local] =
                                exactZero[node.leftNode] &&
                                exactZero[node.rightNode];
                            for (size_t index = 0;
                                 index < monomialCount;
                                 ++index) {
                                arithmetic =
                                    op == Op::Add
                                    ? addBall(
                                          coefficient(
                                              node.leftNode,
                                              index),
                                          coefficient(
                                              node.rightNode,
                                              index),
                                          coefficient(
                                              local, index))
                                    : subtractBall(
                                          coefficient(
                                              node.leftNode,
                                              index),
                                          coefficient(
                                              node.rightNode,
                                              index),
                                          coefficient(
                                              local, index));
                                if (arithmetic !=
                                        ScaledArithmeticStatus::
                                            Success)
                                    break;
                            }
                            if (arithmetic ==
                                    ScaledArithmeticStatus::
                                        Success) {
                                remainders[local] =
                                    remainders[
                                        node.leftNode];
                                arithmetic = addRemainder(
                                    remainders[local],
                                    remainders[
                                        node.rightNode]);
                            }
                            break;
                        case Op::Multiply:
                            exactReal[local] =
                                exactReal[node.leftNode] &&
                                exactReal[node.rightNode];
                            exactZero[local] =
                                exactZero[node.leftNode] ||
                                exactZero[node.rightNode];
                            arithmetic =
                                multiplyPolynomials(
                                    sampleOffset,
                                    node.leftNode,
                                    node.rightNode, false,
                                    static_cast<uint16_t>(
                                        local));
                            break;
                        case Op::Square:
                            exactReal[local] =
                                exactReal[node.leftNode];
                            exactZero[local] =
                                exactZero[node.leftNode];
                            arithmetic =
                                multiplyPolynomials(
                                    sampleOffset,
                                    node.leftNode,
                                    node.leftNode, false,
                                    static_cast<uint16_t>(
                                        local));
                            break;
                        case Op::Conjugate:
                            exactReal[local] =
                                exactReal[node.leftNode];
                            exactZero[local] =
                                exactZero[node.leftNode];
                            for (size_t index = 0;
                                 index < monomialCount;
                                 ++index) {
                                arithmetic =
                                    conjugateBall(
                                        coefficient(
                                            node.leftNode,
                                            conjugateIndices[
                                                index]),
                                        coefficient(
                                            local, index));
                                if (arithmetic !=
                                        ScaledArithmeticStatus::
                                            Success)
                                    break;
                            }
                            remainders[local] =
                                remainders[node.leftNode];
                            break;
                        case Op::Real:
                            exactReal[local] = 1;
                            exactZero[local] =
                                exactZero[node.leftNode];
                            for (size_t index = 0;
                                 index < monomialCount;
                                 ++index) {
                                arithmetic =
                                    realCoefficient(
                                        node.leftNode,
                                        index,
                                        coefficient(
                                            local, index));
                                if (arithmetic !=
                                        ScaledArithmeticStatus::
                                            Success)
                                    break;
                            }
                            remainders[local] =
                                remainders[node.leftNode];
                            break;
                        case Op::Imaginary:
                            exactReal[local] = 1;
                            exactZero[local] =
                                exactReal[node.leftNode] ||
                                exactZero[node.leftNode];
                            for (size_t index = 0;
                                 index < monomialCount;
                                 ++index) {
                                arithmetic =
                                    imaginaryCoefficient(
                                        node.leftNode,
                                        index,
                                        coefficient(
                                            local, index));
                                if (arithmetic !=
                                        ScaledArithmeticStatus::
                                            Success)
                                    break;
                            }
                            remainders[local] =
                                remainders[node.leftNode];
                            break;
                        case Op::MakeComplex:
                            exactReal[local] =
                                exactZero[node.rightNode];
                            exactZero[local] =
                                exactZero[node.leftNode] &&
                                exactZero[node.rightNode];
                            for (size_t index = 0;
                                 index < monomialCount;
                                 ++index) {
                                ScaledComplexBall realLeft;
                                ScaledComplexBall realRight;
                                ScaledComplexBall imaginary;
                                arithmetic =
                                    realCoefficient(
                                        node.leftNode,
                                        index, realLeft);
                                if (arithmetic ==
                                        ScaledArithmeticStatus::
                                            Success)
                                    arithmetic =
                                        realCoefficient(
                                            node.rightNode,
                                            index,
                                            realRight);
                                if (arithmetic ==
                                        ScaledArithmeticStatus::
                                            Success)
                                    arithmetic =
                                        rotatePlusI(
                                            realRight,
                                            imaginary);
                                if (arithmetic ==
                                        ScaledArithmeticStatus::
                                            Success)
                                    arithmetic = addBall(
                                        realLeft,
                                        imaginary,
                                        coefficient(
                                            local, index));
                                if (arithmetic !=
                                        ScaledArithmeticStatus::
                                            Success)
                                    break;
                            }
                            remainders[local] =
                                remainders[node.leftNode];
                            if (arithmetic ==
                                    ScaledArithmeticStatus::
                                        Success)
                                arithmetic = addRemainder(
                                    remainders[local],
                                    remainders[
                                        node.rightNode]);
                            break;
                        case Op::Norm:
                            exactReal[local] = 1;
                            exactZero[local] =
                                exactZero[node.leftNode];
                            arithmetic =
                                multiplyPolynomials(
                                    sampleOffset,
                                    node.leftNode,
                                    node.leftNode, true,
                                    static_cast<uint16_t>(
                                        local));
                            break;
                        case Op::Abs: {
                            ++candidate.absBranchCount;
                            int sign = 0;
                            ScaledRealValue clearance;
                            arithmetic = certifyAbsBranch(
                                sampleOffset,
                                node.leftNode, sign,
                                clearance);
                            if (arithmetic !=
                                    ScaledArithmeticStatus::
                                        Success)
                                break;
                            if (sign == 0) {
                                candidate.status =
                                    ExpressionTaylorJetStatus::
                                        BranchRejected;
                                candidate.reason =
                                    "absolute-value input enclosure touches or crosses zero";
                                candidate.foldRejected = true;
                                candidate.
                                    foldRejectionIteration =
                                        landing;
                                candidate.
                                    foldRejectionReason =
                                        candidate.reason;
                                break;
                            }
                            if (!candidate.
                                    hasFoldClearance ||
                                compareScaledNonnegative(
                                    clearance,
                                    candidate.
                                        minimumFoldClearance) <
                                    0) {
                                candidate.
                                    minimumFoldClearance =
                                        clearance;
                                candidate.
                                    hasFoldClearance = true;
                            }
                            copyNode(node.leftNode);
                            exactReal[local] = 1;
                            exactZero[local] =
                                exactZero[node.leftNode];
                            if (sign > 0) {
                                ++candidate.
                                    absPositiveCellCount;
                            } else {
                                ++candidate.
                                    absNegativeCellCount;
                                for (size_t index = 0;
                                     index <
                                         monomialCount;
                                     ++index) {
                                    arithmetic =
                                        scaledNegate(
                                            coefficient(
                                                local,
                                                index).
                                                value,
                                            coefficient(
                                                local,
                                                index).
                                                value);
                                    count(
                                        candidate.operations);
                                    if (arithmetic !=
                                            ScaledArithmeticStatus::
                                                Success)
                                        break;
                                    if (coefficient(
                                            local, index).
                                            value.isZero()) {
                                        coefficient(
                                            local, index).
                                            value.re.
                                                mantissa =
                                                    0.0;
                                        coefficient(
                                            local, index).
                                            value.im.
                                                mantissa =
                                                    0.0;
                                    }
                                }
                            }
                            break;
                        }
                        default:
                            arithmetic =
                                ScaledArithmeticStatus::
                                    Singular;
                            break;
                        }
                        if (!candidate.reason.empty())
                            break;
                        if (arithmetic !=
                                ScaledArithmeticStatus::
                                    Success)
                            break;
                        arithmetic = addRemainder(
                            remainders[local],
                            node.outputError);
                        if (arithmetic ==
                                ScaledArithmeticStatus::
                                    Success)
                            arithmetic = addNodeRoundoff(
                                    node, sampleOffset,
                                remainders[local]);
                        if (arithmetic !=
                                ScaledArithmeticStatus::
                                    Success)
                            break;
                        for (size_t index = 0;
                             index < monomialCount;
                             ++index) {
                            if (certifyScaledMpfrExponentRange(
                                    coefficient(
                                        local, index)) !=
                                    ScaledArithmeticStatus::
                                        Success) {
                                arithmetic =
                                    ScaledArithmeticStatus::
                                        ExponentRange;
                                break;
                            }
                        }
                        if (arithmetic !=
                                ScaledArithmeticStatus::
                                    Success ||
                            certifyScaledMpfrExponentRange(
                                remainders[local]) !=
                                ScaledArithmeticStatus::
                                    Success) {
                            arithmetic =
                                ScaledArithmeticStatus::
                                    ExponentRange;
                            break;
                        }
                    }

                    if (!candidate.reason.empty())
                        break;
                    if (arithmetic !=
                            ScaledArithmeticStatus::Success) {
                        candidate.status =
                            statusForArithmetic(arithmetic);
                        candidate.reason =
                            "real-bivariate Taylor bytecode propagation failed";
                        break;
                    }
                    const uint16_t root = sample.rootNode;
                    if (compareScaledNonnegative(
                            remainders[root],
                            accuracyBudget) > 0) {
                        candidate.status =
                            ExpressionTaylorJetStatus::
                                AccuracyBudget;
                        candidate.reason =
                            "real-bivariate Taylor truncation radius exceeds accuracy budget";
                        break;
                    }

                    ScaledComplexValue outputBase;
                    arithmetic = makeScaledComplexValue(
                        sample.next, sample.rootDefect,
                        outputBase);
                    ScaledRealValue margin;
                    const InsideStatus inside =
                        arithmetic ==
                                ScaledArithmeticStatus::
                                    Success
                        ? certifyInside(
                              outputBase,
                              &coefficient(root, 0),
                              monomialCount,
                              remainders[root],
                              threshold, margin,
                              arithmetic)
                        : InsideStatus::Error;
                    if (inside != InsideStatus::Inside) {
                        candidate.status =
                            inside ==
                                InsideStatus::Uncertain
                            ? ExpressionTaylorJetStatus::
                                  BailoutUncertain
                            : statusForArithmetic(arithmetic);
                        candidate.reason =
                            "real-bivariate Taylor prefix reaches an uncertain or escaping state";
                        break;
                    }

                    if (static_cast<size_t>(landing + 1) ==
                            request.reference->
                                samples.size()) {
                        for (size_t index = 0;
                             index < monomialCount;
                             ++index)
                            state[index] =
                                coefficient(root, index);
                        stateRemainder =
                            remainders[root];
                        ++landing;
                        candidate.margins.push_back(
                            margin);
                        candidate.
                            landingUsesSampleOutput = true;
                        break;
                    }

                    const ExpressionReferenceSample&
                        nextSample =
                            request.reference->samples[
                                static_cast<size_t>(
                                    landing + 1)];
                    ScaledComplexValue nextBaseValue;
                    arithmetic = makeScaledComplexValue(
                        nextSample.z,
                        nextSample.zDefect,
                        nextBaseValue);
                    ScaledComplexBall rootBase;
                    ScaledComplexBall nextBase;
                    rootBase.value = outputBase;
                    nextBase.value = nextBaseValue;
                    ScaledComplexBall rebase;
                    if (arithmetic ==
                            ScaledArithmeticStatus::Success)
                        arithmetic = subtractBall(
                            rootBase, nextBase, rebase);
                    for (size_t index = 0;
                         arithmetic ==
                             ScaledArithmeticStatus::
                                 Success &&
                         index < monomialCount;
                         ++index)
                        nextState[index] =
                            coefficient(root, index);
                    if (arithmetic ==
                            ScaledArithmeticStatus::Success)
                        arithmetic = addBall(
                            nextState[0], rebase,
                            nextState[0]);
                    if (arithmetic !=
                            ScaledArithmeticStatus::Success) {
                        candidate.status =
                            statusForArithmetic(arithmetic);
                        candidate.reason =
                            "real-bivariate Taylor compact rebase failed";
                        break;
                    }
                    state.swap(nextState);
                    stateRemainder = remainders[root];
                    ++landing;
                    candidate.margins.push_back(margin);
                }

                candidate.landing = landing;
                candidate.coefficients =
                    std::move(state);
                candidate.remainder = stateRemainder;
                if (candidate.reason.empty())
                    candidate.status =
                        ExpressionTaylorJetStatus::Success;
            } catch (const std::bad_alloc&) {
                candidate.status =
                    ExpressionTaylorJetStatus::
                        ResourceLimit;
                candidate.reason =
                    "real-bivariate Taylor workspace allocation failed";
            } catch (const std::length_error&) {
                if (candidate.reason.empty()) {
                    candidate.status =
                        ExpressionTaylorJetStatus::
                            ResourceLimit;
                    candidate.reason =
                        "real-bivariate Taylor workspace length overflow";
                }
            }

            totalOperations =
                candidate.operations >
                    std::numeric_limits<uint64_t>::max() -
                        totalOperations
                ? std::numeric_limits<uint64_t>::max()
                : totalOperations + candidate.operations;
            totalConvolutionOperations =
                candidate.bivariateConvolutionOperations >
                    std::numeric_limits<uint64_t>::max() -
                        totalConvolutionOperations
                ? std::numeric_limits<uint64_t>::max()
                : totalConvolutionOperations +
                      candidate.
                          bivariateConvolutionOperations;
            lastStatus = candidate.status;
            lastReason = candidate.reason.empty()
                ? expressionTaylorJetStatusName(
                      candidate.status)
                : candidate.reason;
            if (!haveBest ||
                candidate.landing > best.landing ||
                (candidate.landing == best.landing &&
                 candidate.order < best.order)) {
                best = std::move(candidate);
                haveBest = true;
            }
            if (haveBest &&
                best.landing >= maximumLanding)
                break;
            if (lastStatus ==
                    ExpressionTaylorJetStatus::Cancelled ||
                lastStatus ==
                    ExpressionTaylorJetStatus::
                        ResourceLimit ||
                lastStatus ==
                    ExpressionTaylorJetStatus::
                        ExponentRange ||
                lastStatus ==
                    ExpressionTaylorJetStatus::Nonfinite)
                break;
        }

        result.operationCount = totalOperations;
        result.bivariateConvolutionOperationCount =
            totalConvolutionOperations;
        result.memoryBytes = peakMemory;
        result.pixelParameter = request.pixelParameter;
        result.parameterScale = request.parameterScale;
        result.programSemanticHash =
            request.program->semanticHash();
        if (lastStatus ==
                ExpressionTaylorJetStatus::Cancelled)
            return finish(
                ExpressionTaylorJetStatus::Cancelled,
                lastReason, false);
        if (!haveBest)
            return finish(
                lastStatus, lastReason, false);

        try {
            result.coefficients.reserve(
                best.coefficients.size());
            result.coefficientRadii.reserve(
                best.coefficients.size());
        } catch (const std::bad_alloc&) {
            return finish(
                ExpressionTaylorJetStatus::ResourceLimit,
                "real-bivariate Taylor result allocation failed",
                false);
        } catch (const std::length_error&) {
            return finish(
                ExpressionTaylorJetStatus::ResourceLimit,
                "real-bivariate Taylor result length overflow",
                false);
        }
        size_t resultPeak = validationBytes;
        if (!checkedAddSize(
                resultPeak,
                best.coefficients.capacity(),
                sizeof(ScaledComplexBall)) ||
            !checkedAddSize(
                resultPeak, best.margins.capacity(),
                sizeof(ScaledRealValue)) ||
            !checkedAddSize(
                resultPeak, result.coefficients.capacity(),
                sizeof(ScaledComplexValue)) ||
            !checkedAddSize(
                resultPeak,
                result.coefficientRadii.capacity(),
                sizeof(ScaledRealValue)))
            return finish(
                ExpressionTaylorJetStatus::ResourceLimit,
                "real-bivariate Taylor result memory calculation overflow",
                false);
        peakMemory = std::max(peakMemory, resultPeak);
        result.memoryBytes = peakMemory;
        if (request.memoryLimitBytes != 0 &&
            peakMemory > request.memoryLimitBytes)
            return finish(
                ExpressionTaylorJetStatus::ResourceLimit,
                "real-bivariate Taylor retained result exceeds memory policy",
                false);

        result.order = best.order;
        result.monomialCount =
            best.coefficients.size();
        result.landingIteration = best.landing;
        result.landingSample =
            best.landingUsesSampleOutput
            ? static_cast<size_t>(best.landing - 1)
            : static_cast<size_t>(best.landing);
        result.landingUsesSampleOutput =
            best.landingUsesSampleOutput;
        result.absBranchCount =
            best.absBranchCount;
        result.absPositiveCellCount =
            best.absPositiveCellCount;
        result.absNegativeCellCount =
            best.absNegativeCellCount;
        result.minimumFoldClearance =
            best.minimumFoldClearance;
        result.foldRejected =
            best.foldRejected;
        result.foldRejectionIteration =
            best.foldRejectionIteration;
        result.foldRejectionReason =
            best.foldRejectionReason;
        result.remainderRadius = best.remainder;
        result.intermediateEscapeMargins =
            std::move(best.margins);
        for (const ScaledComplexBall& coefficient :
             best.coefficients) {
            result.coefficients.push_back(
                coefficient.value);
            result.coefficientRadii.push_back(
                coefficient.radius);
        }
        if (best.landing < request.minimumLanding)
            return finish(
                best.status ==
                    ExpressionTaylorJetStatus::Success
                    ? ExpressionTaylorJetStatus::NoCoverage
                    : best.status,
                best.reason.empty()
                    ? "real-bivariate Taylor prefix is shorter than minimum landing"
                    : best.reason,
                false);

        result.valid = true;
        result.certified = true;
        return finish(
            ExpressionTaylorJetStatus::Success,
            best.reason, true);
    }
};

bool ExpressionTaylorJetBuilder::build(
        const ExpressionTaylorJetRequest& request,
        ExpressionTaylorJetResult& result) {
    const Clock::time_point start = Clock::now();
    result = {};
    if (request.program && request.program->valid() &&
        (request.program->scaledResidualCapability() ==
             ExpressionScaledResidualCapability::
                 CertifiedRealCandidate ||
         request.program->scaledResidualCapability() ==
             ExpressionScaledResidualCapability::
                 CertifiedPiecewiseCandidate))
        return ExpressionRealTaylorJetBuilder::build(
            request, result);
    auto finish = [&](ExpressionTaylorJetStatus status,
                      const std::string& reason,
                      bool success) {
        result.status = status;
        result.failureReason = reason;
        result.buildSeconds =
            std::chrono::duration<double>(
                Clock::now() - start).count();
        return success;
    };

    if (!request.program || !request.reference ||
        !request.program->valid() ||
        !request.reference->valid ||
        request.reference->status !=
            ExpressionReferenceBuildStatus::Success ||
        !request.reference->certifiedAgainstHigherPrecision ||
        (request.program->scaledResidualCapability() !=
             ExpressionScaledResidualCapability::
                 ExactCenteredArithmetic &&
         request.program->scaledResidualCapability() !=
             ExpressionScaledResidualCapability::
                 CertifiedEntireCandidate &&
         request.program->scaledResidualCapability() !=
             ExpressionScaledResidualCapability::
                 CertifiedMeromorphicCandidate &&
         request.program->scaledResidualCapability() !=
             ExpressionScaledResidualCapability::
                 CertifiedBranchCandidate) ||
        (request.pixelParameter != FormulaParameter::C &&
         request.pixelParameter !=
             FormulaParameter::InitialZ) ||
        request.minimumOrder < 8 ||
        request.maximumOrder > 20 ||
        request.maximumCompositionOrder < 1 ||
        request.maximumCompositionOrder > 24 ||
        request.minimumOrder > request.preferredOrder ||
        request.preferredOrder > request.maximumOrder ||
        request.minimumLanding < 1 ||
        !(request.bailout > 0.0) ||
        !std::isfinite(request.bailout) ||
        !(request.accuracyBudget > 0.0) ||
        !std::isfinite(request.accuracyBudget) ||
        !request.parameterScale.isNormalized() ||
        request.parameterScale.isZero())
        return finish(
            ExpressionTaylorJetStatus::InvalidRequest,
            "Taylor jet request is invalid", false);
    if (request.reference->programSemanticHash !=
            request.program->semanticHash() ||
        request.reference->programSource !=
            request.program->source())
        return finish(
            ExpressionTaylorJetStatus::InvalidTape,
            "Taylor reference semantic identity mismatch",
            false);
    if (request.reference->samples.size() < 2)
        return finish(
            ExpressionTaylorJetStatus::NoCoverage,
            "Taylor reference has no skippable prefix",
            false);
    if (certifyScaledMpfrExponentRange(
            request.parameterScale) !=
            ScaledArithmeticStatus::Success)
        return finish(
            ExpressionTaylorJetStatus::ExponentRange,
            "Taylor parameter scale is outside MPFR guards",
            false);

    using Op = ExpressionProgram::Op;
    for (const ExpressionProgram::Instruction& instruction :
         request.program->_code) {
        switch (instruction.op) {
        case Op::Constant:
        case Op::Z:
        case Op::C:
        case Op::Z0:
        case Op::Iteration:
        case Op::Parameter:
        case Op::Negate:
        case Op::Add:
        case Op::Subtract:
        case Op::Multiply:
        case Op::Divide:
        case Op::Square:
        case Op::Exp:
        case Op::Sin:
        case Op::Cos:
        case Op::Tan:
        case Op::Sinh:
        case Op::Cosh:
        case Op::Tanh:
        case Op::Log:
        case Op::Log10:
        case Op::Sqrt:
        case Op::Power:
            break;
        default:
            return finish(
                ExpressionTaylorJetStatus::
                    UnsupportedProgram,
                "operation has no univariate complex-q Taylor semantics",
                false);
        }
    }

    size_t validationBytes = 0;
    if (!ExpressionScaledResidualEvaluator::
            estimateWorkspaceBytes(
                *request.program, *request.reference,
                validationBytes))
        return finish(
            ExpressionTaylorJetStatus::ResourceLimit,
            "Taylor tape validation workspace size overflow",
            false);
    if (request.memoryLimitBytes != 0 &&
        validationBytes > request.memoryLimitBytes)
        return finish(
            ExpressionTaylorJetStatus::ResourceLimit,
            "Taylor tape validation exceeds memory policy",
            false);

    ExpressionScaledResidualEvaluator validator;
    try {
        if (!validator.prepare(
                *request.program, *request.reference))
            return finish(
                ExpressionTaylorJetStatus::InvalidTape,
                validator.error().empty()
                    ? "Taylor reference tape validation failed"
                    : validator.error(),
                false);
    } catch (const std::bad_alloc&) {
        return finish(
            ExpressionTaylorJetStatus::ResourceLimit,
            "Taylor tape validation allocation failed",
            false);
    } catch (const std::length_error&) {
        return finish(
            ExpressionTaylorJetStatus::ResourceLimit,
            "Taylor tape validation length overflow",
            false);
    }
    if (validator.workspaceBytes() > validationBytes)
        return finish(
            ExpressionTaylorJetStatus::ResourceLimit,
            "Taylor tape validation exceeded its preflight bound",
            false);

    ScaledRealValue accuracyBudget;
    if (makeScaledRealValue(
            request.accuracyBudget, accuracyBudget) !=
            ScaledArithmeticStatus::Success)
        return finish(
            ExpressionTaylorJetStatus::InvalidRequest,
            "Taylor accuracy budget is not representable",
            false);
    TaylorThreshold threshold;
    if (!makeThreshold(request.bailout, threshold))
        return finish(
            ExpressionTaylorJetStatus::InvalidRequest,
            "Taylor bailout is not representable", false);

    const int maximumLanding = std::min<int>(
        request.maximumCandidateIteration > 0
            ? request.maximumCandidateIteration
            : static_cast<int>(
                  request.reference->samples.size()),
        static_cast<int>(
            request.reference->samples.size()));
    if (maximumLanding < request.minimumLanding)
        return finish(
            ExpressionTaylorJetStatus::NoCoverage,
            "Taylor reference is shorter than minimum landing",
            false);

    Candidate best;
    bool haveBest = false;
    const bool meromorphicProgram =
        request.program->scaledResidualCapability() ==
            ExpressionScaledResidualCapability::
                CertifiedMeromorphicCandidate;
    const bool branchProgram =
        request.program->scaledResidualCapability() ==
            ExpressionScaledResidualCapability::
                CertifiedBranchCandidate;
    uint64_t totalOperations = 0;
    size_t peakMemory = 0;
    ExpressionTaylorJetStatus lastStatus =
        ExpressionTaylorJetStatus::NoCoverage;
    std::string lastReason = "Taylor prefix has no coverage";
    size_t branchMpfrScratchBytes = 0;
    if (branchProgram) {
        const size_t precisionBits =
            static_cast<size_t>(
                request.reference->certificationPrecision) + 64;
        const size_t limbs =
            (precisionBits + GMP_NUMB_BITS - 1) /
            GMP_NUMB_BITS;
        size_t mpfrValueBytes = sizeof(__mpfr_struct);
        if (!checkedAddSize(
                mpfrValueBytes, limbs,
                sizeof(mp_limb_t)) ||
            !checkedAddSize(
                branchMpfrScratchBytes, 14,
                mpfrValueBytes))
            return finish(
                ExpressionTaylorJetStatus::ResourceLimit,
                "branch Taylor MPFR scratch size overflow",
                false);
    }

    for (int order = request.minimumOrder;
         order <= request.maximumOrder; ++order) {
        if (request.shouldCancel && request.shouldCancel())
            return finish(
                ExpressionTaylorJetStatus::Cancelled,
                "Taylor build cancelled", false);
        const size_t stride =
            static_cast<size_t>(order) + 1;
        const size_t nodeCount =
            request.program->instructionCount();
        size_t candidateBytes = validationBytes;
        if (!checkedAddSize(
                candidateBytes, 1,
                branchMpfrScratchBytes)) {
            lastStatus =
                ExpressionTaylorJetStatus::ResourceLimit;
            lastReason =
                "branch Taylor MPFR scratch size overflow";
            break;
        }
        if (haveBest &&
            (!checkedAddSize(
                 candidateBytes,
                 best.coefficients.capacity(),
                 sizeof(ScaledComplexBall)) ||
             !checkedAddSize(
                 candidateBytes,
                 best.margins.capacity(),
                 sizeof(ScaledRealValue)))) {
            lastStatus =
                ExpressionTaylorJetStatus::ResourceLimit;
            lastReason =
                "retained Taylor candidate size overflow";
            break;
        }
        if (!checkedAddSize(
                candidateBytes, nodeCount * stride,
                sizeof(ScaledComplexBall)) ||
            !checkedAddSize(
                candidateBytes, nodeCount,
                sizeof(ScaledRealValue)) ||
            !checkedAddSize(
                candidateBytes,
                stride *
                    ((meromorphicProgram || branchProgram)
                         ? 11 : 8),
                sizeof(ScaledComplexBall)) ||
            !checkedAddSize(
                candidateBytes,
                (static_cast<size_t>(maximumLanding) + 1) *
                    2,
                sizeof(ScaledRealValue))) {
            lastStatus =
                ExpressionTaylorJetStatus::ResourceLimit;
            lastReason = "Taylor workspace size overflow";
            break;
        }
        peakMemory = std::max(peakMemory, candidateBytes);
        if (request.memoryLimitBytes != 0 &&
            candidateBytes > request.memoryLimitBytes) {
            lastStatus =
                ExpressionTaylorJetStatus::ResourceLimit;
            lastReason = "Taylor workspace exceeds memory policy";
            break;
        }

        Candidate candidate;
        candidate.order = order;
        candidate.memoryBytes = candidateBytes;
        try {
            std::vector<ScaledComplexBall> nodes(
                nodeCount * stride);
            std::vector<ScaledRealValue> remainders(
                nodeCount);
            std::vector<ScaledComplexBall> state(stride);
            std::vector<ScaledComplexBall> nextState(stride);
            std::vector<ScaledComplexBall> functionTerm(stride);
            std::vector<ScaledComplexBall> functionScratch(stride);
            std::vector<ScaledComplexBall> functionFirst(stride);
            std::vector<ScaledComplexBall> functionSecond(stride);
            std::vector<ScaledComplexBall>
                meromorphicNumerator;
            std::vector<ScaledComplexBall>
                meromorphicDenominator;
            std::vector<ScaledComplexBall>
                reciprocalCoefficients;
            std::vector<ScaledComplexBall>
                quotientCoefficients;
            if (meromorphicProgram || branchProgram) {
                meromorphicNumerator.resize(stride);
                meromorphicDenominator.resize(stride);
                reciprocalCoefficients.resize(stride);
                quotientCoefficients.resize(stride);
            }
            candidate.margins.reserve(
                static_cast<size_t>(maximumLanding) + 1);

            auto coefficient = [&](
                    size_t node, int degree)
                    -> ScaledComplexBall& {
                return nodes[
                    node * stride +
                    static_cast<size_t>(degree)];
            };
            auto nodeBase = [&](
                    size_t sampleOffset, uint16_t node,
                    ScaledComplexValue& value) {
                const ExpressionReferenceTapeNode& tape =
                    request.reference->tape[
                        sampleOffset + node];
                ScaledArithmeticStatus status =
                    makeScaledComplexValue(
                        tape.output,
                        tape.outputDefect, value);
                if (status !=
                        ScaledArithmeticStatus::Success)
                    return status;
                return certifyScaledMpfrExponentRange(value);
            };
            auto polynomialSup = [&](
                    size_t node,
                    ScaledRealValue& output) {
                output = {};
                for (int degree = 0;
                     degree <= order; ++degree) {
                    ScaledRealValue magnitude;
                    ScaledArithmeticStatus status =
                        upperMagnitude(
                            coefficient(node, degree),
                            magnitude);
                    ++candidate.operations;
                    if (status !=
                            ScaledArithmeticStatus::Success)
                        return status;
                    status = addUp(output, magnitude);
                    ++candidate.operations;
                    if (status !=
                            ScaledArithmeticStatus::Success)
                        return status;
                }
                return ScaledArithmeticStatus::Success;
            };
            auto addRemainder = [&](
                    ScaledRealValue& target,
                    const ScaledRealValue& value) {
                ++candidate.operations;
                return addUp(target, value);
            };
            auto addBall = [&](
                    const ScaledComplexBall& left,
                    const ScaledComplexBall& right,
                    ScaledComplexBall& output) {
                ++candidate.operations;
                return certifiedScaledAdd(
                    left, right, output);
            };
            auto subtractBall = [&](
                    const ScaledComplexBall& left,
                    const ScaledComplexBall& right,
                    ScaledComplexBall& output) {
                ++candidate.operations;
                return certifiedScaledSubtract(
                    left, right, output);
            };
            auto multiplyBall = [&](
                    const ScaledComplexBall& left,
                    const ScaledComplexBall& right,
                    ScaledComplexBall& output) {
                ++candidate.operations;
                return certifiedScaledMultiply(
                    left, right, output);
            };
            auto multiplyBound = [&](
                    const ScaledRealValue& left,
                    const ScaledRealValue& right,
                    ScaledRealValue& output) {
                ++candidate.operations;
                return multiplyUp(left, right, output);
            };
            auto functionAddBall = [&](
                    const ScaledComplexBall& left,
                    const ScaledComplexBall& right,
                    ScaledComplexBall& output) {
                ++candidate.functionSeriesOperations;
                return addBall(left, right, output);
            };
            auto functionSubtractBall = [&](
                    const ScaledComplexBall& left,
                    const ScaledComplexBall& right,
                    ScaledComplexBall& output) {
                ++candidate.functionSeriesOperations;
                return subtractBall(left, right, output);
            };
            auto functionMultiplyBall = [&](
                    const ScaledComplexBall& left,
                    const ScaledComplexBall& right,
                    ScaledComplexBall& output) {
                ++candidate.functionSeriesOperations;
                return multiplyBall(left, right, output);
            };
            auto functionAddRemainder = [&](
                    ScaledRealValue& target,
                    const ScaledRealValue& value) {
                ++candidate.functionSeriesOperations;
                return addRemainder(target, value);
            };
            auto functionMultiplyBound = [&](
                    const ScaledRealValue& left,
                    const ScaledRealValue& right,
                    ScaledRealValue& output) {
                ++candidate.functionSeriesOperations;
                return multiplyBound(left, right, output);
            };
            auto polynomialArraySup = [&](
                    const ScaledComplexBall* coefficients,
                    ScaledRealValue& output) {
                output = {};
                for (int degree = 0;
                     degree <= order; ++degree) {
                    ScaledRealValue magnitude;
                    ScaledArithmeticStatus status =
                        upperMagnitude(
                            coefficients[
                                static_cast<size_t>(degree)],
                            magnitude);
                    ++candidate.operations;
                    ++candidate.functionSeriesOperations;
                    if (status !=
                            ScaledArithmeticStatus::Success)
                        return status;
                    status = functionAddRemainder(
                        output, magnitude);
                    if (status !=
                            ScaledArithmeticStatus::Success)
                        return status;
                }
                return ScaledArithmeticStatus::Success;
            };
            auto multiplyPolynomials = [&](
                    const ScaledComplexBall* left,
                    const ScaledRealValue& leftRemainder,
                    const ScaledComplexBall* right,
                    const ScaledRealValue& rightRemainder,
                    std::vector<ScaledComplexBall>& output,
                    ScaledRealValue& outputRemainder) {
                std::fill(
                    output.begin(), output.end(),
                    ScaledComplexBall{});
                ScaledArithmeticStatus status =
                    ScaledArithmeticStatus::Success;
                for (int degree = 0;
                     degree <= order; ++degree) {
                    bool haveSum = false;
                    ScaledComplexBall sum;
                    for (int first = 0;
                         first <= degree; ++first) {
                        ScaledComplexBall term;
                        status = functionMultiplyBall(
                            left[
                                static_cast<size_t>(first)],
                            right[
                                static_cast<size_t>(
                                    degree - first)],
                            term);
                        if (status !=
                                ScaledArithmeticStatus::Success)
                            return status;
                        if (!haveSum) {
                            sum = term;
                            haveSum = true;
                        } else {
                            status = functionAddBall(
                                sum, term, sum);
                            if (status !=
                                    ScaledArithmeticStatus::
                                        Success)
                                return status;
                        }
                    }
                    output[static_cast<size_t>(degree)] =
                        sum;
                }

                ScaledRealValue leftSup, rightSup;
                status = polynomialArraySup(
                    left, leftSup);
                if (status !=
                        ScaledArithmeticStatus::Success)
                    return status;
                status = polynomialArraySup(
                    right, rightSup);
                if (status !=
                        ScaledArithmeticStatus::Success)
                    return status;

                outputRemainder = {};
                auto addProduct = [&](
                        const ScaledRealValue& first,
                        const ScaledRealValue& second,
                        bool doubleProduct = false) {
                    ScaledRealValue product;
                    ScaledArithmeticStatus local =
                        functionMultiplyBound(
                            first, second, product);
                    if (local !=
                            ScaledArithmeticStatus::Success)
                        return local;
                    if (doubleProduct) {
                        ScaledRealValue doubled;
                        local = scaledAddUp(
                            product, product, doubled);
                        ++candidate.operations;
                        ++candidate.functionSeriesOperations;
                        if (local !=
                                ScaledArithmeticStatus::
                                    Success)
                            return local;
                        product = doubled;
                    }
                    return functionAddRemainder(
                        outputRemainder, product);
                };
                status = addProduct(
                    leftSup, rightRemainder);
                if (status ==
                        ScaledArithmeticStatus::Success)
                    status = addProduct(
                        rightSup, leftRemainder);
                if (status ==
                        ScaledArithmeticStatus::Success)
                    status = addProduct(
                        leftRemainder, rightRemainder,
                        true);
                for (int first = 0;
                     status ==
                         ScaledArithmeticStatus::Success &&
                     first <= order; ++first) {
                    ScaledRealValue firstMagnitude;
                    status = upperMagnitude(
                        left[static_cast<size_t>(first)],
                        firstMagnitude);
                    ++candidate.operations;
                    ++candidate.functionSeriesOperations;
                    for (int second =
                             order + 1 - first;
                         status ==
                             ScaledArithmeticStatus::Success &&
                         second <= order; ++second) {
                        if (second < 0) continue;
                        ScaledRealValue secondMagnitude;
                        status = upperMagnitude(
                            right[
                                static_cast<size_t>(second)],
                            secondMagnitude);
                        ++candidate.operations;
                        ++candidate.functionSeriesOperations;
                        if (status ==
                                ScaledArithmeticStatus::Success)
                            status = addProduct(
                                firstMagnitude,
                                secondMagnitude);
                    }
                }
                return status;
            };
            auto scalePolynomial = [&](
                    const ScaledComplexBall* input,
                    const ScaledRealValue& inputRemainder,
                    const ScaledComplexBall& factor,
                    std::vector<ScaledComplexBall>& output,
                    ScaledRealValue& outputRemainder) {
                ScaledArithmeticStatus status =
                    ScaledArithmeticStatus::Success;
                for (int degree = 0;
                     degree <= order; ++degree) {
                    status = functionMultiplyBall(
                        input[static_cast<size_t>(degree)],
                        factor,
                        output[static_cast<size_t>(degree)]);
                    if (status !=
                            ScaledArithmeticStatus::Success)
                        return status;
                }
                ScaledRealValue factorMagnitude;
                status = upperMagnitude(
                    factor, factorMagnitude);
                ++candidate.operations;
                ++candidate.functionSeriesOperations;
                if (status !=
                        ScaledArithmeticStatus::Success)
                    return status;
                return functionMultiplyBound(
                    factorMagnitude, inputRemainder,
                    outputRemainder);
            };
            auto addPolynomialTerm = [&](
                    std::vector<ScaledComplexBall>& target,
                    ScaledRealValue& targetRemainder,
                    const std::vector<ScaledComplexBall>& term,
                    const ScaledRealValue& termRemainder,
                    bool negate) {
                for (int degree = 0;
                     degree <= order; ++degree) {
                    ScaledComplexBall contribution =
                        term[static_cast<size_t>(degree)];
                    if (negate) {
                        ScaledArithmeticStatus local =
                            scaledNegate(
                                contribution.value,
                                contribution.value);
                        ++candidate.operations;
                        ++candidate.functionSeriesOperations;
                        if (local !=
                                ScaledArithmeticStatus::Success)
                            return local;
                    }
                    ScaledArithmeticStatus local =
                        functionAddBall(
                            target[
                                static_cast<size_t>(degree)],
                            contribution,
                            target[
                                static_cast<size_t>(degree)]);
                    if (local !=
                            ScaledArithmeticStatus::Success)
                        return local;
                }
                return functionAddRemainder(
                    targetRemainder, termRemainder);
            };
            auto reciprocalBall = [&](
                    int divisor,
                    ScaledComplexBall& output) {
                output = {};
                if (divisor <= 0)
                    return ScaledArithmeticStatus::
                        ExponentRange;
                const double reciprocal =
                    1.0 / static_cast<double>(divisor);
                ScaledArithmeticStatus status =
                    makeScaledRealValue(
                        reciprocal, output.value.re);
                if (status !=
                        ScaledArithmeticStatus::Success)
                    return status;
                if (divisor == 1)
                    return status;
                return makeScaledRealValue(
                    std::ldexp(
                        std::abs(reciprocal), -50),
                    output.radius);
            };
            auto addNodeRoundoff = [&](
                    size_t local,
                    const ExpressionReferenceTapeNode& node,
                    size_t sampleOffset,
                    ScaledRealValue& remainder) {
                if (request.reference->
                        certificationPrecision <= 16)
                    return ScaledArithmeticStatus::Success;
                const bool arithmeticOperation =
                    node.operation ==
                        ExpressionOracleOperation::Add ||
                    node.operation ==
                        ExpressionOracleOperation::Subtract ||
                    node.operation ==
                        ExpressionOracleOperation::Multiply ||
                    node.operation ==
                        ExpressionOracleOperation::Square;
                const bool nonlinearOutputOperation =
                    node.operation ==
                        ExpressionOracleOperation::Exp ||
                    node.operation ==
                        ExpressionOracleOperation::Sin ||
                    node.operation ==
                        ExpressionOracleOperation::Cos ||
                    node.operation ==
                        ExpressionOracleOperation::Sinh ||
                    node.operation ==
                        ExpressionOracleOperation::Cosh ||
                    node.operation ==
                        ExpressionOracleOperation::Divide ||
                    node.operation ==
                        ExpressionOracleOperation::Tan ||
                    node.operation ==
                        ExpressionOracleOperation::Tanh;
                const bool branchOutputOperation =
                    node.operation ==
                        ExpressionOracleOperation::Log ||
                    node.operation ==
                        ExpressionOracleOperation::Log10 ||
                    node.operation ==
                        ExpressionOracleOperation::Sqrt ||
                    node.operation ==
                        ExpressionOracleOperation::Power;
                if (!arithmeticOperation &&
                    !nonlinearOutputOperation &&
                    !branchOutputOperation)
                    return ScaledArithmeticStatus::Success;

                auto fullBound = [&](
                        uint16_t child,
                        ScaledRealValue& output) {
                    ScaledComplexValue base;
                    ScaledArithmeticStatus status =
                        nodeBase(
                            sampleOffset, child, base);
                    if (status !=
                            ScaledArithmeticStatus::Success)
                        return status;
                    ScaledRealValue baseMagnitude;
                    status = upperMagnitude(
                        base, baseMagnitude);
                    if (status !=
                            ScaledArithmeticStatus::Success)
                        return status;
                    ScaledRealValue polynomial;
                    status = polynomialSup(
                        child, polynomial);
                    if (status !=
                            ScaledArithmeticStatus::Success)
                        return status;
                    output = baseMagnitude;
                    status = addRemainder(
                        output, polynomial);
                    if (status !=
                            ScaledArithmeticStatus::Success)
                        return status;
                    status = addRemainder(
                        output, remainders[child]);
                    if (status !=
                            ScaledArithmeticStatus::Success)
                        return status;
                    return addRemainder(
                        output, remainders[child]);
                };

                ScaledRealValue left, right;
                ScaledRealValue operationMagnitude;
                ScaledArithmeticStatus status =
                    ScaledArithmeticStatus::Success;
                const bool directOutputOperation =
                    nonlinearOutputOperation ||
                    branchOutputOperation;
                if (directOutputOperation) {
                    status = fullBound(
                        static_cast<uint16_t>(local),
                        operationMagnitude);
                } else {
                    status = fullBound(node.leftNode, left);
                    if (status !=
                            ScaledArithmeticStatus::Success)
                        return status;
                }
                if (!directOutputOperation &&
                    node.operation ==
                        ExpressionOracleOperation::Square) {
                    status = multiplyBound(
                        left, left, operationMagnitude);
                } else if (!directOutputOperation) {
                    status = fullBound(
                        node.rightNode, right);
                    if (status !=
                            ScaledArithmeticStatus::Success)
                        return status;
                    if (node.operation ==
                            ExpressionOracleOperation::Add ||
                        node.operation ==
                            ExpressionOracleOperation::Subtract) {
                        operationMagnitude = left;
                        status = addRemainder(
                            operationMagnitude, right);
                    } else {
                        status = multiplyBound(
                            left, right,
                            operationMagnitude);
                    }
                }
                if (status !=
                        ScaledArithmeticStatus::Success)
                    return status;
                ScaledRealValue scale;
                scale.mantissa = 0.5;
                scale.exponent =
                    9 - static_cast<int64_t>(
                        request.reference->
                            certificationPrecision);
                ScaledRealValue roundoff;
                status = multiplyBound(
                    operationMagnitude, scale, roundoff);
                if (status !=
                        ScaledArithmeticStatus::Success)
                    return status;
                return addRemainder(
                    remainder, roundoff);
            };
            auto composeEntire = [&](
                    size_t local,
                    const ExpressionReferenceTapeNode& node,
                    ScaledRealValue& outputRemainder,
                    ScaledComplexBall* destination = nullptr,
                    ExpressionOracleOperation
                        operationOverride =
                            ExpressionOracleOperation::
                                Constant,
                    const ScaledComplexBall*
                        baseOverride = nullptr,
                    const ScaledComplexBall*
                        companionOverride = nullptr,
                    const ScaledComplexBall*
                        inputOverride = nullptr,
                    const ScaledRealValue*
                        inputRemainderOverride = nullptr) {
                auto failComposition = [&](
                        ExpressionTaylorJetStatus status,
                        const char* reason) {
                    candidate.status = status;
                    candidate.reason = reason;
                    return false;
                };
                if ((node.flags & OracleTraceUndefined) != 0)
                    return failComposition(
                        ExpressionTaylorJetStatus::Nonfinite,
                        "entire-function reference node is undefined");

                const ExpressionOracleOperation operation =
                    operationOverride ==
                        ExpressionOracleOperation::Constant
                    ? node.operation : operationOverride;
                const bool needsCompanion =
                    operation ==
                        ExpressionOracleOperation::Sin ||
                    operation ==
                        ExpressionOracleOperation::Cos ||
                    operation ==
                        ExpressionOracleOperation::Sinh ||
                    operation ==
                        ExpressionOracleOperation::Cosh;
                if (!companionOverride && needsCompanion &&
                    (node.flags & OracleTraceHasCompanion) == 0)
                    return failComposition(
                        ExpressionTaylorJetStatus::InvalidTape,
                        "entire-function reference companion is missing");

                ScaledComplexBall base;
                ScaledArithmeticStatus arithmetic =
                    ScaledArithmeticStatus::Success;
                if (baseOverride) {
                    base = *baseOverride;
                } else {
                    arithmetic = makeScaledComplexValue(
                        node.output, node.outputDefect,
                        base.value);
                    base.radius = node.outputError;
                }
                if (arithmetic ==
                        ScaledArithmeticStatus::Success)
                    arithmetic =
                        certifyScaledMpfrExponentRange(base);
                if (arithmetic !=
                        ScaledArithmeticStatus::Success)
                    return failComposition(
                        statusForArithmetic(arithmetic),
                        "entire-function reference output is outside certified range");

                ScaledComplexBall companion;
                if (needsCompanion) {
                    if (companionOverride) {
                        companion = *companionOverride;
                    } else {
                        arithmetic = makeScaledComplexValue(
                            node.auxiliary,
                            node.auxiliaryDefect,
                            companion.value);
                        companion.radius =
                            node.auxiliaryError;
                    }
                    if (arithmetic ==
                            ScaledArithmeticStatus::Success)
                        arithmetic =
                            certifyScaledMpfrExponentRange(
                                companion);
                    if (arithmetic !=
                            ScaledArithmeticStatus::Success)
                        return failComposition(
                            statusForArithmetic(arithmetic),
                            "entire-function reference companion is outside certified range");
                }

                const uint16_t inputNode = node.leftNode;
                if (!inputOverride &&
                    inputNode == UINT16_MAX)
                    return failComposition(
                        ExpressionTaylorJetStatus::InvalidTape,
                        "entire-function reference input is missing");
                const ScaledComplexBall* inputCoefficients =
                    inputOverride
                    ? inputOverride
                    : &coefficient(inputNode, 0);
                const ScaledRealValue inputRemainder =
                    inputRemainderOverride
                    ? *inputRemainderOverride
                    : remainders[inputNode];
                ScaledRealValue inputRadius;
                arithmetic = polynomialArraySup(
                    inputCoefficients,
                    inputRadius);
                if (arithmetic ==
                        ScaledArithmeticStatus::Success)
                    arithmetic = functionAddRemainder(
                        inputRadius,
                        inputRemainder);
                if (arithmetic ==
                        ScaledArithmeticStatus::Success)
                    arithmetic = functionAddRemainder(
                        inputRadius,
                        inputRemainder);
                if (arithmetic ==
                        ScaledArithmeticStatus::Success)
                    arithmetic =
                        certifyScaledMpfrExponentRange(
                            inputRadius);
                if (arithmetic !=
                        ScaledArithmeticStatus::Success)
                    return failComposition(
                        statusForArithmetic(arithmetic),
                        "entire-function input radius is outside certified range");

                std::fill(
                    functionFirst.begin(),
                    functionFirst.end(),
                    ScaledComplexBall{});
                std::fill(
                    functionSecond.begin(),
                    functionSecond.end(),
                    ScaledComplexBall{});
                for (int degree = 0;
                     degree <= order; ++degree)
                    functionTerm[
                        static_cast<size_t>(degree)] =
                            inputCoefficients[
                                static_cast<size_t>(degree)];
                ScaledRealValue termRemainder =
                    inputRemainder;
                ScaledRealValue firstRemainder;
                ScaledRealValue secondRemainder;
                ScaledRealValue certifiedTail;
                int retainedOrder = 0;

                for (int seriesOrder = 1;
                     seriesOrder <=
                         request.maximumCompositionOrder;
                     ++seriesOrder) {
                    if (seriesOrder > 1) {
                        ScaledRealValue productRemainder;
                        arithmetic = multiplyPolynomials(
                            functionTerm.data(),
                            termRemainder,
                            inputCoefficients,
                            inputRemainder,
                            functionScratch,
                            productRemainder);
                        if (arithmetic !=
                                ScaledArithmeticStatus::
                                    Success)
                            return failComposition(
                                statusForArithmetic(arithmetic),
                                "entire-function polynomial power failed");
                        ScaledComplexBall reciprocal;
                        arithmetic = reciprocalBall(
                            seriesOrder, reciprocal);
                        if (arithmetic !=
                                ScaledArithmeticStatus::
                                    Success)
                            return failComposition(
                                statusForArithmetic(arithmetic),
                                "entire-function reciprocal is not representable");
                        arithmetic = scalePolynomial(
                            functionScratch.data(),
                            productRemainder,
                            reciprocal,
                            functionTerm,
                            termRemainder);
                        if (arithmetic !=
                                ScaledArithmeticStatus::
                                    Success)
                            return failComposition(
                                statusForArithmetic(arithmetic),
                                "entire-function Taylor term scaling failed");
                    }

                    const bool exponential =
                        operation ==
                            ExpressionOracleOperation::Exp;
                    const bool hyperbolic =
                        operation ==
                            ExpressionOracleOperation::Sinh ||
                        operation ==
                            ExpressionOracleOperation::Cosh;
                    const bool even =
                        (seriesOrder & 1) == 0;
                    const bool useFirst =
                        exponential || even;
                    bool negate = false;
                    if (!exponential && !hyperbolic) {
                        if (even)
                            negate =
                                ((seriesOrder / 2) & 1) != 0;
                        else
                            negate =
                                (((seriesOrder - 1) / 2) &
                                 1) != 0;
                    }
                    arithmetic = addPolynomialTerm(
                        useFirst
                            ? functionFirst
                            : functionSecond,
                        useFirst
                            ? firstRemainder
                            : secondRemainder,
                        functionTerm,
                        termRemainder,
                        negate);
                    if (arithmetic !=
                            ScaledArithmeticStatus::Success)
                        return failComposition(
                            statusForArithmetic(arithmetic),
                            "entire-function Taylor accumulation failed");

                    ScaledRealValue scalarTail;
                    arithmetic = finiteSeriesTailBound(
                        inputRadius, seriesOrder,
                        scalarTail);
                    ++candidate.functionSeriesOperations;
                    ++candidate.operations;
                    if (arithmetic !=
                            ScaledArithmeticStatus::Success)
                        return failComposition(
                            statusForArithmetic(arithmetic),
                            "entire-function finite-series tail is not certifiable");

                    ScaledRealValue weight;
                    arithmetic = upperMagnitude(
                        base, weight);
                    ++candidate.operations;
                    ++candidate.functionSeriesOperations;
                    if (arithmetic ==
                            ScaledArithmeticStatus::Success &&
                        needsCompanion) {
                        ScaledRealValue companionMagnitude;
                        arithmetic = upperMagnitude(
                            companion,
                            companionMagnitude);
                        ++candidate.operations;
                        ++candidate.functionSeriesOperations;
                        if (arithmetic ==
                                ScaledArithmeticStatus::Success)
                            arithmetic =
                                functionAddRemainder(
                                    weight,
                                    companionMagnitude);
                    }
                    if (arithmetic ==
                            ScaledArithmeticStatus::Success)
                        arithmetic = functionMultiplyBound(
                            weight, scalarTail,
                            certifiedTail);
                    if (arithmetic !=
                            ScaledArithmeticStatus::Success)
                        return failComposition(
                            statusForArithmetic(arithmetic),
                            "entire-function weighted tail overflowed");
                    if (compareScaledNonnegative(
                            certifiedTail,
                            accuracyBudget) <= 0) {
                        retainedOrder = seriesOrder;
                        break;
                    }
                }
                if (retainedOrder == 0)
                    return failComposition(
                        ExpressionTaylorJetStatus::
                            AccuracyBudget,
                        "entire-function finite-series tail exceeds accuracy budget");

                ScaledRealValue firstOutputRemainder;
                arithmetic = scalePolynomial(
                    functionFirst.data(),
                    firstRemainder,
                    base,
                    functionScratch,
                    firstOutputRemainder);
                if (arithmetic !=
                        ScaledArithmeticStatus::Success)
                    return failComposition(
                        statusForArithmetic(arithmetic),
                        "entire-function reference output multiplication failed");

                if (!needsCompanion) {
                    for (int degree = 0;
                         degree <= order; ++degree)
                        (destination
                            ? destination[
                                  static_cast<size_t>(degree)]
                            : coefficient(local, degree)) =
                            functionScratch[
                                static_cast<size_t>(degree)];
                    outputRemainder =
                        firstOutputRemainder;
                } else {
                    ScaledRealValue secondOutputRemainder;
                    arithmetic = scalePolynomial(
                        functionSecond.data(),
                        secondRemainder,
                        companion,
                        functionTerm,
                        secondOutputRemainder);
                    if (arithmetic !=
                            ScaledArithmeticStatus::Success)
                        return failComposition(
                            statusForArithmetic(arithmetic),
                            "entire-function companion multiplication failed");
                    const bool subtractCompanion =
                        operation ==
                            ExpressionOracleOperation::Cos;
                    for (int degree = 0;
                         degree <= order; ++degree) {
                        arithmetic = subtractCompanion
                            ? functionSubtractBall(
                                  functionScratch[
                                      static_cast<size_t>(
                                          degree)],
                                  functionTerm[
                                      static_cast<size_t>(
                                          degree)],
                                  destination
                                      ? destination[
                                            static_cast<size_t>(
                                                degree)]
                                      : coefficient(
                                            local, degree))
                            : functionAddBall(
                                  functionScratch[
                                      static_cast<size_t>(
                                          degree)],
                                  functionTerm[
                                      static_cast<size_t>(
                                          degree)],
                                  destination
                                      ? destination[
                                            static_cast<size_t>(
                                                degree)]
                                      : coefficient(
                                            local, degree));
                        if (arithmetic !=
                                ScaledArithmeticStatus::
                                    Success)
                            return failComposition(
                                statusForArithmetic(arithmetic),
                                "entire-function reference combination failed");
                    }
                    outputRemainder =
                        firstOutputRemainder;
                    arithmetic = functionAddRemainder(
                        outputRemainder,
                        secondOutputRemainder);
                    if (arithmetic !=
                            ScaledArithmeticStatus::Success)
                        return failComposition(
                            statusForArithmetic(arithmetic),
                            "entire-function remainder combination failed");
                }
                arithmetic = functionAddRemainder(
                    outputRemainder, certifiedTail);
                if (arithmetic !=
                        ScaledArithmeticStatus::Success)
                    return failComposition(
                        statusForArithmetic(arithmetic),
                        "entire-function tail accumulation failed");

                ++candidate.functionSeriesCount;
                candidate.maximumFunctionSeriesOrder =
                    std::max(
                        candidate.maximumFunctionSeriesOrder,
                        retainedOrder);
                if (compareScaledNonnegative(
                        certifiedTail,
                        candidate.maximumFunctionSeriesTail) >
                    0)
                    candidate.maximumFunctionSeriesTail =
                        certifiedTail;
                return true;
            };
            auto composePrincipalLogVariation = [&](
                    const ExpressionReferenceTapeNode& node,
                    size_t sampleOffset,
                    uint16_t inputNode,
                    bool base10,
                    ScaledComplexBall* destination,
                    ScaledRealValue& outputRemainder,
                    ScaledComplexBall& baseFunction,
                    ScaledComplexValue*
                        expansionCenter = nullptr) {
                auto failBranch = [&](
                        ExpressionTaylorJetStatus status,
                        const char* reason,
                        bool rejected = false) {
                    candidate.status = status;
                    candidate.reason = reason;
                    candidate.branchRejected =
                        candidate.branchRejected || rejected;
                    return false;
                };
                if (!destination || inputNode == UINT16_MAX)
                    return failBranch(
                        ExpressionTaylorJetStatus::InvalidTape,
                        "branch-function input tape is missing");
                ScaledComplexBall inputCenter;
                ScaledArithmeticStatus arithmetic =
                    nodeBase(
                        sampleOffset, inputNode,
                        inputCenter.value);
                if (arithmetic !=
                        ScaledArithmeticStatus::Success)
                    return failBranch(
                        statusForArithmetic(arithmetic),
                        "branch-function reference input is outside certified range");
                arithmetic = addBall(
                    inputCenter,
                    coefficient(inputNode, 0),
                    inputCenter);
                if (arithmetic !=
                        ScaledArithmeticStatus::Success)
                    return failBranch(
                        statusForArithmetic(arithmetic),
                        "branch-function constant center failed");
                if (expansionCenter)
                    *expansionCenter = inputCenter.value;

                std::fill(
                    meromorphicNumerator.begin(),
                    meromorphicNumerator.end(),
                    ScaledComplexBall{});
                for (int degree = 1;
                     degree <= order; ++degree)
                    meromorphicNumerator[
                        static_cast<size_t>(degree)] =
                            coefficient(inputNode, degree);
                ScaledRealValue deltaRemainder =
                    remainders[inputNode];
                arithmetic = functionAddRemainder(
                    deltaRemainder, inputCenter.radius);
                inputCenter.radius = {};
                if (arithmetic !=
                        ScaledArithmeticStatus::Success)
                    return failBranch(
                        statusForArithmetic(arithmetic),
                        "branch-function input uncertainty overflowed");

                ScaledRealValue enclosureVariation =
                    deltaRemainder;
                for (int degree = 1;
                     degree <= order; ++degree) {
                    ScaledRealValue magnitude;
                    arithmetic = upperMagnitude(
                        meromorphicNumerator[
                            static_cast<size_t>(degree)],
                        magnitude);
                    ++candidate.operations;
                    if (arithmetic ==
                            ScaledArithmeticStatus::Success)
                        arithmetic = functionAddRemainder(
                            enclosureVariation, magnitude);
                    if (arithmetic !=
                            ScaledArithmeticStatus::Success)
                        return failBranch(
                            statusForArithmetic(arithmetic),
                            "branch-function enclosure overflowed");
                }

                ScaledRealValue zeroClearance;
                ScaledRealValue cutClearance;
                PrincipalClearanceFailure clearanceFailure =
                    PrincipalClearanceFailure::None;
                arithmetic = principalBranchClearance(
                    inputCenter.value,
                    enclosureVariation,
                    zeroClearance,
                    cutClearance,
                    clearanceFailure);
                ++candidate.operations;
                auto recordClearance = [](
                        const ScaledRealValue& clearance,
                        ScaledRealValue& minimum,
                        bool& haveMinimum) {
                    if (!haveMinimum ||
                        compareScaledNonnegative(
                            clearance, minimum) < 0) {
                        minimum = clearance;
                        haveMinimum = true;
                    }
                };
                recordClearance(
                    zeroClearance,
                    candidate.minimumBranchZeroClearance,
                    candidate.hasBranchZeroClearance);
                recordClearance(
                    cutClearance,
                    candidate.minimumBranchCutClearance,
                    candidate.hasBranchCutClearance);
                if (arithmetic ==
                        ScaledArithmeticStatus::Singular) {
                    return failBranch(
                        ExpressionTaylorJetStatus::
                            BranchRejected,
                        clearanceFailure ==
                                PrincipalClearanceFailure::Zero
                            ? "branch input enclosure touches or contains zero"
                            : "branch input enclosure touches the negative-real cut or has an ambiguous lip",
                        true);
                }
                if (arithmetic !=
                        ScaledArithmeticStatus::Success)
                    return failBranch(
                        statusForArithmetic(arithmetic),
                        "branch-function clearance is outside certified range");
                if ((node.flags & OracleTraceUndefined) != 0)
                    return failBranch(
                        ExpressionTaylorJetStatus::Nonfinite,
                        "branch-function reference node is undefined");

                ScaledComplexBall reciprocal;
                arithmetic = reciprocalScaledValue(
                    inputCenter.value, reciprocal);
                ++candidate.operations;
                if (arithmetic !=
                        ScaledArithmeticStatus::Success)
                    return failBranch(
                        ExpressionTaylorJetStatus::
                            BranchRejected,
                        "branch expansion center is zero or not representable",
                        true);
                ScaledRealValue normalizedRemainder;
                arithmetic = scalePolynomial(
                    meromorphicNumerator.data(),
                    deltaRemainder,
                    reciprocal,
                    meromorphicDenominator,
                    normalizedRemainder);
                if (arithmetic !=
                        ScaledArithmeticStatus::Success)
                    return failBranch(
                        statusForArithmetic(arithmetic),
                        "branch normalized polynomial failed");

                ScaledRealValue rho;
                arithmetic = polynomialArraySup(
                    meromorphicDenominator.data(), rho);
                if (arithmetic ==
                        ScaledArithmeticStatus::Success)
                    arithmetic = functionAddRemainder(
                        rho, normalizedRemainder);
                if (arithmetic ==
                        ScaledArithmeticStatus::Success)
                    arithmetic = functionAddRemainder(
                        rho, normalizedRemainder);
                ScaledRealValue one;
                if (arithmetic ==
                        ScaledArithmeticStatus::Success)
                    arithmetic = makeScaledRealValue(
                        1.0, one);
                if (arithmetic !=
                        ScaledArithmeticStatus::Success)
                    return failBranch(
                        statusForArithmetic(arithmetic),
                        "branch normalized radius overflowed");
                if (compareScaledNonnegative(rho, one) >= 0)
                    return failBranch(
                        ExpressionTaylorJetStatus::
                            AccuracyBudget,
                        "branch normalized radius is not strictly below one");

                std::fill(
                    functionFirst.begin(),
                    functionFirst.end(),
                    ScaledComplexBall{});
                for (int degree = 0;
                     degree <= order; ++degree)
                    functionTerm[
                        static_cast<size_t>(degree)] =
                            meromorphicDenominator[
                                static_cast<size_t>(degree)];
                ScaledRealValue termRemainder =
                    normalizedRemainder;
                ScaledRealValue seriesRemainder;
                ScaledRealValue certifiedTail;
                ScaledRealValue branchTailBudget;
                arithmetic = scaledDivideByDouble(
                    accuracyBudget, 16.0,
                    branchTailBudget);
                ++candidate.operations;
                if (arithmetic !=
                        ScaledArithmeticStatus::Success)
                    return failBranch(
                        statusForArithmetic(arithmetic),
                        "branch tail budget is not representable");
                int retainedOrder = 0;
                for (int seriesOrder = 1;
                     seriesOrder <=
                         request.maximumCompositionOrder;
                     ++seriesOrder) {
                    ScaledComplexBall reciprocalOrder;
                    arithmetic = reciprocalBall(
                        seriesOrder, reciprocalOrder);
                    if (arithmetic !=
                            ScaledArithmeticStatus::Success)
                        return failBranch(
                            statusForArithmetic(arithmetic),
                            "branch-series reciprocal is not representable");
                    ScaledRealValue scaledTermRemainder;
                    arithmetic = scalePolynomial(
                        functionTerm.data(),
                        termRemainder,
                        reciprocalOrder,
                        functionScratch,
                        scaledTermRemainder);
                    if (arithmetic !=
                            ScaledArithmeticStatus::Success)
                        return failBranch(
                            statusForArithmetic(arithmetic),
                            "branch-series term scaling failed");
                    arithmetic = addPolynomialTerm(
                        functionFirst,
                        seriesRemainder,
                        functionScratch,
                        scaledTermRemainder,
                        (seriesOrder & 1) == 0);
                    if (arithmetic !=
                            ScaledArithmeticStatus::Success)
                        return failBranch(
                            statusForArithmetic(arithmetic),
                            "branch-series accumulation failed");

                    arithmetic = logarithmSeriesTailBound(
                        rho, seriesOrder, certifiedTail);
                    ++candidate.operations;
                    ++candidate.functionSeriesOperations;
                    if (arithmetic !=
                            ScaledArithmeticStatus::Success)
                        return failBranch(
                            ExpressionTaylorJetStatus::
                                AccuracyBudget,
                            "branch logarithm tail is not certifiable");
                    if (compareScaledNonnegative(
                            certifiedTail,
                            branchTailBudget) <= 0) {
                        retainedOrder = seriesOrder;
                        break;
                    }
                    ScaledRealValue nextRemainder;
                    arithmetic = multiplyPolynomials(
                        functionTerm.data(),
                        termRemainder,
                        meromorphicDenominator.data(),
                        normalizedRemainder,
                        functionScratch,
                        nextRemainder);
                    if (arithmetic !=
                            ScaledArithmeticStatus::Success)
                        return failBranch(
                            statusForArithmetic(arithmetic),
                            "branch polynomial power failed");
                    functionTerm.swap(functionScratch);
                    termRemainder = nextRemainder;
                }
                if (retainedOrder == 0)
                    return failBranch(
                        ExpressionTaylorJetStatus::
                            AccuracyBudget,
                        "branch logarithm tail exceeds accuracy budget");

                const ScaledComplexBall* series =
                    functionFirst.data();
                if (base10) {
                    ScaledComplexBall reciprocalLn10;
                    arithmetic = makeReciprocalLn10Ball(
                        request.reference->
                            certificationPrecision,
                        reciprocalLn10);
                    ++candidate.operations;
                    if (arithmetic !=
                            ScaledArithmeticStatus::Success)
                        return failBranch(
                            statusForArithmetic(arithmetic),
                            "certified reciprocal ln(10) is unavailable");
                    ScaledRealValue scaledRemainder;
                    arithmetic = scalePolynomial(
                        functionFirst.data(),
                        seriesRemainder,
                        reciprocalLn10,
                        functionSecond,
                        scaledRemainder);
                    if (arithmetic !=
                            ScaledArithmeticStatus::Success)
                        return failBranch(
                            statusForArithmetic(arithmetic),
                            "log10 polynomial scaling failed");
                    ScaledRealValue factorMagnitude;
                    arithmetic = upperMagnitude(
                        reciprocalLn10, factorMagnitude);
                    ++candidate.operations;
                    if (arithmetic ==
                            ScaledArithmeticStatus::Success)
                        arithmetic = functionMultiplyBound(
                            factorMagnitude,
                            certifiedTail,
                            certifiedTail);
                    if (arithmetic !=
                            ScaledArithmeticStatus::Success)
                        return failBranch(
                            statusForArithmetic(arithmetic),
                            "log10 tail scaling failed");
                    seriesRemainder = scaledRemainder;
                    series = functionSecond.data();
                }
                for (int degree = 0;
                     degree <= order; ++degree)
                    destination[
                        static_cast<size_t>(degree)] =
                            series[
                                static_cast<size_t>(degree)];
                outputRemainder = seriesRemainder;
                arithmetic = functionAddRemainder(
                    outputRemainder, certifiedTail);
                if (arithmetic !=
                        ScaledArithmeticStatus::Success)
                    return failBranch(
                        statusForArithmetic(arithmetic),
                        "branch-series tail accumulation failed");

                arithmetic = makePrincipalFunctionBall(
                    inputCenter.value,
                    base10
                        ? PrincipalFunction::Log10
                        : PrincipalFunction::Log,
                    request.reference->
                        certificationPrecision,
                    baseFunction);
                ++candidate.operations;
                if (arithmetic !=
                        ScaledArithmeticStatus::Success)
                    return failBranch(
                        statusForArithmetic(arithmetic),
                        "principal logarithm base is outside certified range");

                ++candidate.branchCompositionCount;
                candidate.maximumBranchSeriesOrder =
                    std::max(
                        candidate.maximumBranchSeriesOrder,
                        retainedOrder);
                if (compareScaledNonnegative(
                        certifiedTail,
                        candidate.maximumBranchSeriesTail) >
                    0)
                    candidate.maximumBranchSeriesTail =
                        certifiedTail;
                return true;
            };
            auto inflateReferenceValue = [&](
                    ScaledComplexBall& value) {
                if (request.reference->
                        certificationPrecision <= 16)
                    return ScaledArithmeticStatus::Success;
                ScaledRealValue magnitude;
                ScaledArithmeticStatus status =
                    upperMagnitude(value, magnitude);
                ++candidate.operations;
                ++candidate.reciprocalOperations;
                if (status !=
                        ScaledArithmeticStatus::Success)
                    return status;
                ScaledRealValue scale;
                scale.mantissa = 0.5;
                scale.exponent =
                    14 - static_cast<int64_t>(
                        request.reference->
                            certificationPrecision);
                ScaledRealValue roundoff;
                status = multiplyBound(
                    magnitude, scale, roundoff);
                ++candidate.reciprocalOperations;
                if (status !=
                        ScaledArithmeticStatus::Success)
                    return status;
                status = addRemainder(
                    value.radius, roundoff);
                ++candidate.reciprocalOperations;
                return status;
            };
            auto composeReciprocal = [&](
                    const ScaledComplexBall* denominator,
                    const ScaledRealValue&
                        denominatorRemainder,
                    std::vector<ScaledComplexBall>& output,
                    ScaledRealValue& outputRemainder,
                    const ScaledRealValue*
                        clearanceScale = nullptr) {
                auto failReciprocal = [&](
                        ExpressionTaylorJetStatus status,
                        const char* reason,
                        bool pole = false) {
                    candidate.status = status;
                    candidate.reason = reason;
                    candidate.poleRejected =
                        candidate.poleRejected || pole;
                    return false;
                };
                if (!denominator)
                    return failReciprocal(
                        ExpressionTaylorJetStatus::
                            InvalidTape,
                        "reciprocal denominator is missing");

                ScaledComplexBall constant = denominator[0];
                ScaledRealValue inputRemainder =
                    denominatorRemainder;
                ScaledArithmeticStatus arithmetic =
                    addRemainder(
                        inputRemainder, constant.radius);
                ++candidate.reciprocalOperations;
                constant.radius = {};
                if (arithmetic !=
                        ScaledArithmeticStatus::Success)
                    return failReciprocal(
                        statusForArithmetic(arithmetic),
                        "reciprocal constant uncertainty overflowed");

                ScaledRealValue polynomialVariation;
                // Reciprocal clearance and tails use the maximum-component
                // norm throughout. Complex products therefore carry the
                // sharp-enough factor two used below.
                for (int degree = 1;
                     degree <= order; ++degree) {
                    ScaledRealValue magnitude;
                    arithmetic = upperMaxComponent(
                        denominator[
                            static_cast<size_t>(degree)],
                        magnitude);
                    ++candidate.operations;
                    ++candidate.reciprocalOperations;
                    if (arithmetic !=
                            ScaledArithmeticStatus::Success)
                        return failReciprocal(
                            statusForArithmetic(arithmetic),
                            "reciprocal denominator coefficient bound failed");
                    arithmetic = addRemainder(
                        polynomialVariation, magnitude);
                    ++candidate.reciprocalOperations;
                    if (arithmetic !=
                            ScaledArithmeticStatus::Success)
                        return failReciprocal(
                            statusForArithmetic(arithmetic),
                            "reciprocal denominator variation overflowed");
                }
                arithmetic = addRemainder(
                    polynomialVariation, inputRemainder);
                ++candidate.reciprocalOperations;
                if (arithmetic !=
                        ScaledArithmeticStatus::Success)
                    return failReciprocal(
                        statusForArithmetic(arithmetic),
                        "reciprocal denominator remainder overflowed");

                ScaledRealValue clearance;
                arithmetic =
                    lowerMaxComponentAfterVariation(
                    constant.value,
                    polynomialVariation, clearance);
                ++candidate.operations;
                ++candidate.reciprocalOperations;
                if (arithmetic ==
                        ScaledArithmeticStatus::Singular)
                    return failReciprocal(
                        ExpressionTaylorJetStatus::
                            PoleRejected,
                        "denominator ball overlaps zero",
                        true);
                if (arithmetic !=
                        ScaledArithmeticStatus::Success)
                    return failReciprocal(
                        statusForArithmetic(arithmetic),
                        "denominator clearance is outside certified range");
                ScaledRealValue conservativeClearance;
                arithmetic = scaledDivideByDouble(
                    clearance, 2.0,
                    conservativeClearance);
                ++candidate.operations;
                ++candidate.reciprocalOperations;
                if (arithmetic !=
                        ScaledArithmeticStatus::Success ||
                    conservativeClearance.isZero())
                    return failReciprocal(
                        ExpressionTaylorJetStatus::
                            PoleRejected,
                        "denominator clearance has no conservative margin",
                        true);
                ScaledRealValue reportedClearance =
                    conservativeClearance;
                if (clearanceScale) {
                    // If both factors have max-component lower bounds a,b,
                    // their product has max-component magnitude at least
                    // a*b/sqrt(2); division by two is conservative.
                    arithmetic = multiplyNonnegativeDown(
                        clearance, *clearanceScale,
                        reportedClearance);
                    ++candidate.operations;
                    ++candidate.reciprocalOperations;
                    if (arithmetic ==
                            ScaledArithmeticStatus::Success)
                        arithmetic = scaledDivideByDouble(
                            reportedClearance, 2.0,
                            reportedClearance);
                    ++candidate.operations;
                    ++candidate.reciprocalOperations;
                    if (arithmetic !=
                            ScaledArithmeticStatus::Success)
                        return failReciprocal(
                            statusForArithmetic(arithmetic),
                            "denominator clearance product overflowed");
                }
                if (!candidate.hasDenominatorClearance ||
                    compareScaledNonnegative(
                        reportedClearance,
                        candidate.
                            minimumDenominatorClearance) < 0) {
                    candidate.minimumDenominatorClearance =
                        reportedClearance;
                    candidate.hasDenominatorClearance = true;
                }

                std::fill(
                    output.begin(), output.end(),
                    ScaledComplexBall{});
                arithmetic = reciprocalScaledValue(
                    constant.value, output[0]);
                ++candidate.operations;
                ++candidate.reciprocalOperations;
                if (arithmetic !=
                        ScaledArithmeticStatus::Success)
                    return failReciprocal(
                        statusForArithmetic(arithmetic),
                        "reciprocal constant is not representable");

                for (int degree = 1;
                     degree <= order; ++degree) {
                    bool haveSum = false;
                    ScaledComplexBall sum;
                    for (int first = 1;
                         first <= degree; ++first) {
                        ScaledComplexBall term;
                        arithmetic = multiplyBall(
                            denominator[
                                static_cast<size_t>(first)],
                            output[
                                static_cast<size_t>(
                                    degree - first)],
                            term);
                        ++candidate.reciprocalOperations;
                        if (arithmetic !=
                                ScaledArithmeticStatus::
                                    Success)
                            break;
                        if (!haveSum) {
                            sum = term;
                            haveSum = true;
                        } else {
                            arithmetic = addBall(
                                sum, term, sum);
                            ++candidate.reciprocalOperations;
                            if (arithmetic !=
                                    ScaledArithmeticStatus::
                                        Success)
                                break;
                        }
                    }
                    if (arithmetic !=
                            ScaledArithmeticStatus::Success)
                        break;
                    ScaledComplexBall coefficient;
                    arithmetic = multiplyBall(
                        output[0], sum, coefficient);
                    ++candidate.reciprocalOperations;
                    if (arithmetic !=
                            ScaledArithmeticStatus::Success)
                        break;
                    arithmetic = scaledNegate(
                        coefficient.value,
                        coefficient.value);
                    ++candidate.operations;
                    ++candidate.reciprocalOperations;
                    if (arithmetic !=
                            ScaledArithmeticStatus::Success)
                        break;
                    output[
                        static_cast<size_t>(degree)] =
                            coefficient;
                }
                if (arithmetic !=
                        ScaledArithmeticStatus::Success)
                    return failReciprocal(
                        statusForArithmetic(arithmetic),
                        "formal reciprocal coefficient division failed");

                ScaledRealValue identityResidual;
                for (int degree = 0;
                     degree <= 2 * order; ++degree) {
                    bool haveSum = false;
                    ScaledComplexBall sum;
                    const int firstMinimum =
                        std::max(0, degree - order);
                    const int firstMaximum =
                        std::min(order, degree);
                    for (int first = firstMinimum;
                         first <= firstMaximum; ++first) {
                        const ScaledComplexBall&
                            denominatorCoefficient =
                                first == 0
                                ? constant
                                : denominator[
                                      static_cast<size_t>(
                                          first)];
                        ScaledComplexBall term;
                        arithmetic = multiplyBall(
                            denominatorCoefficient,
                            output[
                                static_cast<size_t>(
                                    degree - first)],
                            term);
                        ++candidate.reciprocalOperations;
                        if (arithmetic !=
                                ScaledArithmeticStatus::
                                    Success)
                            break;
                        if (!haveSum) {
                            sum = term;
                            haveSum = true;
                        } else {
                            arithmetic = addBall(
                                sum, term, sum);
                            ++candidate.reciprocalOperations;
                            if (arithmetic !=
                                    ScaledArithmeticStatus::
                                        Success)
                                break;
                        }
                    }
                    if (arithmetic !=
                            ScaledArithmeticStatus::Success)
                        break;
                    if (degree == 0) {
                        ScaledComplexBall one;
                        arithmetic = makeScaledRealValue(
                            1.0, one.value.re);
                        ++candidate.operations;
                        ++candidate.reciprocalOperations;
                        if (arithmetic ==
                                ScaledArithmeticStatus::
                                    Success)
                            arithmetic = subtractBall(
                                one, sum, sum);
                        ++candidate.reciprocalOperations;
                    } else {
                        arithmetic = scaledNegate(
                            sum.value, sum.value);
                        ++candidate.operations;
                        ++candidate.reciprocalOperations;
                    }
                    if (arithmetic !=
                            ScaledArithmeticStatus::Success)
                        break;
                    ScaledRealValue magnitude;
                    arithmetic = upperMaxComponent(
                        sum, magnitude);
                    ++candidate.operations;
                    ++candidate.reciprocalOperations;
                    if (arithmetic ==
                            ScaledArithmeticStatus::Success)
                        arithmetic = addRemainder(
                            identityResidual, magnitude);
                    ++candidate.reciprocalOperations;
                    if (arithmetic !=
                            ScaledArithmeticStatus::Success)
                        break;
                }
                if (arithmetic !=
                        ScaledArithmeticStatus::Success)
                    return failReciprocal(
                        statusForArithmetic(arithmetic),
                        "reciprocal identity residual bound failed");

                ScaledRealValue reciprocalSup;
                for (int degree = 0;
                     degree <= order; ++degree) {
                    ScaledRealValue magnitude;
                    arithmetic = upperMaxComponent(
                        output[
                            static_cast<size_t>(degree)],
                        magnitude);
                    ++candidate.operations;
                    ++candidate.reciprocalOperations;
                    if (arithmetic ==
                            ScaledArithmeticStatus::Success)
                        arithmetic = addRemainder(
                            reciprocalSup, magnitude);
                    ++candidate.reciprocalOperations;
                    if (arithmetic !=
                            ScaledArithmeticStatus::Success)
                        break;
                }
                ScaledRealValue inputAmplification;
                if (arithmetic ==
                        ScaledArithmeticStatus::Success)
                    arithmetic = multiplyBound(
                        inputRemainder, reciprocalSup,
                        inputAmplification);
                ++candidate.reciprocalOperations;
                if (arithmetic ==
                        ScaledArithmeticStatus::Success) {
                    ScaledRealValue doubled;
                    arithmetic = scaledAddUp(
                        inputAmplification,
                        inputAmplification, doubled);
                    ++candidate.operations;
                    ++candidate.reciprocalOperations;
                    if (arithmetic ==
                            ScaledArithmeticStatus::Success)
                        inputAmplification = doubled;
                }
                if (arithmetic ==
                        ScaledArithmeticStatus::Success)
                    arithmetic = addRemainder(
                        identityResidual,
                        inputAmplification);
                ++candidate.reciprocalOperations;
                if (arithmetic !=
                        ScaledArithmeticStatus::Success)
                    return failReciprocal(
                        statusForArithmetic(arithmetic),
                        "reciprocal input-remainder amplification failed");

                ScaledRealValue doubledResidual;
                if (arithmetic ==
                        ScaledArithmeticStatus::Success)
                    arithmetic = scaledAddUp(
                        identityResidual,
                        identityResidual,
                        doubledResidual);
                ++candidate.operations;
                ++candidate.reciprocalOperations;
                if (arithmetic ==
                        ScaledArithmeticStatus::Success)
                    arithmetic = divideNonnegativeUp(
                    doubledResidual, clearance,
                    outputRemainder);
                ++candidate.operations;
                ++candidate.reciprocalOperations;
                if (arithmetic !=
                        ScaledArithmeticStatus::Success)
                    return failReciprocal(
                        statusForArithmetic(arithmetic),
                        "reciprocal omitted tail is not representable");
                ++candidate.reciprocalCount;
                candidate.maximumReciprocalOrder =
                    std::max(
                        candidate.maximumReciprocalOrder,
                        order);
                if (compareScaledNonnegative(
                        outputRemainder,
                        candidate.maximumReciprocalTail) >
                    0)
                    candidate.maximumReciprocalTail =
                        outputRemainder;
                return true;
            };
            auto dividePolynomials = [&](
                    const ScaledComplexBall* numerator,
                    const ScaledRealValue&
                        numeratorRemainder,
                    const ScaledComplexBall* denominator,
                    const ScaledRealValue&
                        denominatorRemainder,
                    const ScaledComplexBall&
                        referenceOutput,
                    ScaledComplexBall* output,
                    ScaledRealValue& outputRemainder,
                    const ScaledRealValue*
                        clearanceScale = nullptr) {
                ScaledRealValue reciprocalRemainder;
                if (!composeReciprocal(
                        denominator,
                        denominatorRemainder,
                        reciprocalCoefficients,
                        reciprocalRemainder,
                        clearanceScale))
                    return false;
                const uint64_t operationsBefore =
                    candidate.operations;
                ScaledArithmeticStatus arithmetic =
                    multiplyPolynomials(
                        numerator, numeratorRemainder,
                        reciprocalCoefficients.data(),
                        reciprocalRemainder,
                        quotientCoefficients,
                        outputRemainder);
                candidate.reciprocalOperations +=
                    candidate.operations - operationsBefore;
                if (arithmetic !=
                        ScaledArithmeticStatus::Success) {
                    candidate.status =
                        statusForArithmetic(arithmetic);
                    candidate.reason =
                        "division polynomial multiplication failed";
                    return false;
                }
                for (int degree = 0;
                     degree <= order; ++degree)
                    output[static_cast<size_t>(degree)] =
                        quotientCoefficients[
                            static_cast<size_t>(degree)];
                arithmetic = subtractBall(
                    output[0], referenceOutput,
                    output[0]);
                ++candidate.reciprocalOperations;
                if (arithmetic !=
                        ScaledArithmeticStatus::Success) {
                    candidate.status =
                        statusForArithmetic(arithmetic);
                    candidate.reason =
                        "division reference rebase failed";
                    return false;
                }
                return true;
            };

            if (request.pixelParameter ==
                    FormulaParameter::InitialZ) {
                if (order < 1) {
                    candidate.status =
                        ExpressionTaylorJetStatus::
                            InvalidRequest;
                    candidate.reason =
                        "Taylor order cannot represent the parameter";
                } else {
                    state[1].value =
                        request.parameterScale;
                }
            }
            ScaledRealValue stateRemainder =
                request.reference->samples[0].zError;
            if (request.pixelParameter ==
                    FormulaParameter::InitialZ &&
                addRemainder(
                    stateRemainder,
                    request.reference->z0Error) !=
                    ScaledArithmeticStatus::Success) {
                candidate.status =
                    ExpressionTaylorJetStatus::
                        ExponentRange;
                candidate.reason =
                    "initial Taylor radius overflow";
            }

            if (candidate.reason.empty()) {
                ScaledComplexValue initialBase;
                ScaledArithmeticStatus arithmetic =
                    makeScaledComplexValue(
                        request.reference->samples[0].z,
                        request.reference->samples[0].
                            zDefect,
                        initialBase);
                ScaledRealValue margin;
                const InsideStatus inside =
                    arithmetic ==
                            ScaledArithmeticStatus::Success
                    ? certifyInside(
                          initialBase, state.data(), order,
                          stateRemainder, threshold, margin,
                          arithmetic)
                    : InsideStatus::Error;
                if (inside == InsideStatus::Inside) {
                    candidate.margins.push_back(margin);
                } else {
                    candidate.status =
                        inside == InsideStatus::Uncertain
                        ? ExpressionTaylorJetStatus::
                              BailoutUncertain
                        : statusForArithmetic(arithmetic);
                    candidate.reason =
                        "initial Taylor frame overlaps bailout";
                }
            }

            int landing = 0;
            while (candidate.reason.empty() &&
                   landing < maximumLanding) {
                if (request.shouldCancel &&
                    request.shouldCancel()) {
                    candidate.status =
                        ExpressionTaylorJetStatus::Cancelled;
                    candidate.reason =
                        "Taylor build cancelled";
                    break;
                }
                const ExpressionReferenceSample& sample =
                    request.reference->samples[
                        static_cast<size_t>(landing)];
                const size_t sampleOffset =
                    static_cast<size_t>(sample.tapeOffset);
                std::fill(nodes.begin(), nodes.end(),
                          ScaledComplexBall{});
                std::fill(remainders.begin(),
                          remainders.end(),
                          ScaledRealValue{});
                ScaledArithmeticStatus arithmetic =
                    ScaledArithmeticStatus::Success;

                for (size_t local = 0;
                     local < sample.tapeCount; ++local) {
                    const ExpressionReferenceTapeNode& node =
                        request.reference->tape[
                            sampleOffset + local];
                    const Op op =
                        request.program->_code[local].op;
                    auto copyNode = [&](uint16_t source) {
                        for (int degree = 0;
                             degree <= order; ++degree)
                            coefficient(local, degree) =
                                coefficient(source, degree);
                        remainders[local] =
                            remainders[source];
                    };

                    switch (op) {
                    case Op::Constant:
                    case Op::Iteration:
                        break;
                    case Op::Z:
                        for (int degree = 0;
                             degree <= order; ++degree)
                            coefficient(local, degree) =
                                state[
                                    static_cast<size_t>(
                                        degree)];
                        remainders[local] =
                            stateRemainder;
                        break;
                    case Op::C:
                        if (request.pixelParameter ==
                                FormulaParameter::C)
                            coefficient(local, 1).value =
                                request.parameterScale;
                        remainders[local] =
                            request.reference->cError;
                        break;
                    case Op::Z0:
                        if (request.pixelParameter ==
                                FormulaParameter::
                                    InitialZ)
                            coefficient(local, 1).value =
                                request.parameterScale;
                        remainders[local] =
                            request.reference->z0Error;
                        break;
                    case Op::Parameter:
                        break;
                    case Op::Negate:
                        copyNode(node.leftNode);
                        for (int degree = 0;
                             degree <= order; ++degree) {
                            arithmetic = scaledNegate(
                                coefficient(
                                    node.leftNode,
                                    degree).value,
                                coefficient(
                                    local, degree).value);
                            ++candidate.operations;
                            if (arithmetic !=
                                    ScaledArithmeticStatus::
                                        Success)
                                break;
                        }
                        break;
                    case Op::Add:
                    case Op::Subtract:
                        for (int degree = 0;
                             degree <= order; ++degree) {
                            arithmetic =
                                op == Op::Add
                                ? addBall(
                                      coefficient(
                                          node.leftNode,
                                          degree),
                                      coefficient(
                                          node.rightNode,
                                          degree),
                                      coefficient(
                                          local, degree))
                                : subtractBall(
                                      coefficient(
                                          node.leftNode,
                                          degree),
                                      coefficient(
                                          node.rightNode,
                                          degree),
                                      coefficient(
                                          local, degree));
                            if (arithmetic !=
                                    ScaledArithmeticStatus::
                                        Success)
                                break;
                        }
                        if (arithmetic ==
                                ScaledArithmeticStatus::
                                    Success) {
                            remainders[local] =
                                remainders[node.leftNode];
                            arithmetic = addRemainder(
                                remainders[local],
                                remainders[node.rightNode]);
                        }
                        break;
                    case Op::Multiply:
                    case Op::Square: {
                        const uint16_t left = node.leftNode;
                        const uint16_t right =
                            op == Op::Square
                            ? node.leftNode
                            : node.rightNode;
                        ScaledComplexValue leftBaseValue;
                        ScaledComplexValue rightBaseValue;
                        arithmetic = nodeBase(
                            sampleOffset, left,
                            leftBaseValue);
                        if (arithmetic !=
                                ScaledArithmeticStatus::
                                    Success)
                            break;
                        arithmetic = nodeBase(
                            sampleOffset, right,
                            rightBaseValue);
                        if (arithmetic !=
                                ScaledArithmeticStatus::
                                    Success)
                            break;
                        ScaledComplexBall leftBase;
                        ScaledComplexBall rightBase;
                        leftBase.value = leftBaseValue;
                        rightBase.value = rightBaseValue;
                        for (int degree = 0;
                             degree <= order; ++degree) {
                            ScaledComplexBall sum;
                            ScaledComplexBall term;
                            arithmetic = multiplyBall(
                                leftBase,
                                coefficient(right, degree),
                                sum);
                            if (arithmetic !=
                                    ScaledArithmeticStatus::
                                        Success)
                                break;
                            arithmetic = multiplyBall(
                                rightBase,
                                coefficient(left, degree),
                                term);
                            if (arithmetic !=
                                    ScaledArithmeticStatus::
                                        Success)
                                break;
                            arithmetic = addBall(
                                sum, term, sum);
                            if (arithmetic !=
                                    ScaledArithmeticStatus::
                                        Success)
                                break;
                            for (int first = 0;
                                 first <= degree; ++first) {
                                arithmetic = multiplyBall(
                                    coefficient(
                                        left, first),
                                    coefficient(
                                        right,
                                        degree - first),
                                    term);
                                if (arithmetic !=
                                        ScaledArithmeticStatus::
                                            Success)
                                    break;
                                arithmetic = addBall(
                                    sum, term, sum);
                                if (arithmetic !=
                                        ScaledArithmeticStatus::
                                            Success)
                                    break;
                            }
                            if (arithmetic !=
                                    ScaledArithmeticStatus::
                                        Success)
                                break;
                            coefficient(local, degree) =
                                sum;
                        }
                        if (arithmetic !=
                                ScaledArithmeticStatus::
                                    Success)
                            break;

                        ScaledRealValue leftBaseMagnitude;
                        ScaledRealValue rightBaseMagnitude;
                        ScaledRealValue leftSup;
                        ScaledRealValue rightSup;
                        arithmetic = upperMagnitude(
                            leftBaseValue,
                            leftBaseMagnitude);
                        if (arithmetic ==
                                ScaledArithmeticStatus::
                                    Success)
                            arithmetic = upperMagnitude(
                                rightBaseValue,
                                rightBaseMagnitude);
                        if (arithmetic ==
                                ScaledArithmeticStatus::
                                    Success)
                            arithmetic = polynomialSup(
                                left, leftSup);
                        if (arithmetic ==
                                ScaledArithmeticStatus::
                                    Success)
                            arithmetic = polynomialSup(
                                right, rightSup);
                        if (arithmetic !=
                                ScaledArithmeticStatus::
                                    Success)
                            break;

                        ScaledRealValue tail;
                        auto addProduct = [&](
                                const ScaledRealValue& first,
                                const ScaledRealValue& second,
                                bool bothAreComponentRadii =
                                    false) {
                            ScaledRealValue product;
                            ScaledArithmeticStatus status =
                                multiplyBound(
                                    first, second, product);
                            if (status !=
                                    ScaledArithmeticStatus::
                                        Success)
                                return status;
                            if (bothAreComponentRadii) {
                                ScaledRealValue doubled;
                                status = scaledAddUp(
                                    product, product,
                                    doubled);
                                ++candidate.operations;
                                if (status !=
                                        ScaledArithmeticStatus::
                                            Success)
                                    return status;
                                product = doubled;
                            }
                            return addRemainder(
                                tail, product);
                        };
                        arithmetic = addProduct(
                            leftBaseMagnitude,
                            remainders[right]);
                        if (arithmetic ==
                                ScaledArithmeticStatus::
                                    Success)
                            arithmetic = addProduct(
                                rightBaseMagnitude,
                                remainders[left]);
                        if (arithmetic ==
                                ScaledArithmeticStatus::
                                    Success)
                            arithmetic = addProduct(
                                leftSup,
                                remainders[right]);
                        if (arithmetic ==
                                ScaledArithmeticStatus::
                                    Success)
                            arithmetic = addProduct(
                                rightSup,
                                remainders[left]);
                        if (arithmetic ==
                                ScaledArithmeticStatus::
                                    Success)
                            arithmetic = addProduct(
                                remainders[left],
                                remainders[right], true);
                        for (int first = 0;
                             arithmetic ==
                                 ScaledArithmeticStatus::
                                     Success &&
                             first <= order; ++first) {
                            ScaledRealValue firstMagnitude;
                            arithmetic = upperMagnitude(
                                coefficient(left, first),
                                firstMagnitude);
                            for (int second =
                                     order + 1 - first;
                                 arithmetic ==
                                     ScaledArithmeticStatus::
                                         Success &&
                                 second <= order; ++second) {
                                if (second < 0) continue;
                                ScaledRealValue
                                    secondMagnitude;
                                arithmetic =
                                    upperMagnitude(
                                        coefficient(
                                            right, second),
                                        secondMagnitude);
                                if (arithmetic ==
                                        ScaledArithmeticStatus::
                                            Success)
                                    arithmetic = addProduct(
                                        firstMagnitude,
                                        secondMagnitude);
                            }
                        }
                        if (arithmetic ==
                                ScaledArithmeticStatus::
                                    Success)
                            remainders[local] = tail;
                        break;
                    }
                    case Op::Divide: {
                        if ((node.flags &
                             OracleTraceSingularPoint) != 0) {
                            candidate.status =
                                ExpressionTaylorJetStatus::
                                    PoleRejected;
                            candidate.reason =
                                "division reference is an exact pole";
                            candidate.poleRejected = true;
                            break;
                        }
                        if ((node.flags &
                             OracleTraceUndefined) != 0) {
                            candidate.status =
                                ExpressionTaylorJetStatus::
                                    Nonfinite;
                            candidate.reason =
                                "division reference is undefined";
                            break;
                        }
                        if ((node.flags &
                             OracleTraceHasDenominator) == 0 ||
                            node.leftNode == UINT16_MAX ||
                            node.rightNode == UINT16_MAX) {
                            candidate.status =
                                ExpressionTaylorJetStatus::
                                    InvalidTape;
                            candidate.reason =
                                "division denominator tape is missing";
                            break;
                        }
                        ScaledComplexBall numeratorBase;
                        ScaledComplexBall denominatorBase;
                        arithmetic = nodeBase(
                            sampleOffset, node.leftNode,
                            numeratorBase.value);
                        if (arithmetic ==
                                ScaledArithmeticStatus::
                                    Success)
                            arithmetic = nodeBase(
                                sampleOffset, node.rightNode,
                                denominatorBase.value);
                        if (arithmetic !=
                                ScaledArithmeticStatus::
                                    Success)
                            break;
                        arithmetic = addBall(
                            numeratorBase,
                            coefficient(
                                node.leftNode, 0),
                            meromorphicNumerator[0]);
                        if (arithmetic ==
                                ScaledArithmeticStatus::
                                    Success)
                            arithmetic = addBall(
                                denominatorBase,
                                coefficient(
                                    node.rightNode, 0),
                                meromorphicDenominator[0]);
                        for (int degree = 1;
                             arithmetic ==
                                 ScaledArithmeticStatus::
                                     Success &&
                             degree <= order; ++degree) {
                            meromorphicNumerator[
                                static_cast<size_t>(degree)] =
                                    coefficient(
                                        node.leftNode,
                                        degree);
                            meromorphicDenominator[
                                static_cast<size_t>(degree)] =
                                    coefficient(
                                        node.rightNode,
                                        degree);
                        }
                        ScaledComplexBall referenceOutput;
                        if (arithmetic ==
                                ScaledArithmeticStatus::
                                    Success)
                            arithmetic =
                                makeScaledComplexValue(
                                    node.output,
                                    node.outputDefect,
                                    referenceOutput.value);
                        if (arithmetic !=
                                ScaledArithmeticStatus::
                                    Success)
                            break;
                        if (!dividePolynomials(
                                meromorphicNumerator.data(),
                                remainders[node.leftNode],
                                meromorphicDenominator.data(),
                                remainders[node.rightNode],
                                referenceOutput,
                                &coefficient(local, 0),
                                remainders[local]))
                            arithmetic =
                                ScaledArithmeticStatus::
                                    Success;
                        break;
                    }
                    case Op::Exp:
                    case Op::Sin:
                    case Op::Cos:
                    case Op::Sinh:
                    case Op::Cosh:
                        if (!composeEntire(
                                local, node,
                                remainders[local]))
                            arithmetic =
                                ScaledArithmeticStatus::
                                    Success;
                        break;
                    case Op::Log:
                    case Op::Log10: {
                        const uint64_t operationsBefore =
                            candidate.operations;
                        ScaledComplexBall baseFunction;
                        if (!composePrincipalLogVariation(
                                node, sampleOffset,
                                node.leftNode,
                                op == Op::Log10,
                                &coefficient(local, 0),
                                remainders[local],
                                baseFunction)) {
                            arithmetic =
                                ScaledArithmeticStatus::
                                    Success;
                            break;
                        }
                        ScaledComplexBall referenceOutput;
                        arithmetic = makeScaledComplexValue(
                            node.output, node.outputDefect,
                            referenceOutput.value);
                        ScaledComplexBall baseDelta;
                        if (arithmetic ==
                                ScaledArithmeticStatus::
                                    Success)
                            arithmetic = subtractBall(
                                baseFunction,
                                referenceOutput,
                                baseDelta);
                        if (arithmetic ==
                                ScaledArithmeticStatus::
                                    Success)
                            arithmetic = addBall(
                                coefficient(local, 0),
                                baseDelta,
                                coefficient(local, 0));
                        candidate.branchCompositionOperations +=
                            candidate.operations -
                            operationsBefore;
                        break;
                    }
                    case Op::Sqrt: {
                        const uint64_t operationsBefore =
                            candidate.operations;
                        ScaledComplexBall logBase;
                        ScaledComplexValue inputCenter;
                        ScaledRealValue logRemainder;
                        if (!composePrincipalLogVariation(
                                node, sampleOffset,
                                node.leftNode, false,
                                meromorphicNumerator.data(),
                                logRemainder, logBase,
                                &inputCenter)) {
                            arithmetic =
                                ScaledArithmeticStatus::
                                    Success;
                            break;
                        }
                        ScaledComplexBall oneHalf;
                        arithmetic = makeScaledRealValue(
                            0.5, oneHalf.value.re);
                        ScaledRealValue halfLogRemainder;
                        if (arithmetic ==
                                ScaledArithmeticStatus::
                                    Success)
                            arithmetic = scalePolynomial(
                                meromorphicNumerator.data(),
                                logRemainder,
                                oneHalf,
                                meromorphicDenominator,
                                halfLogRemainder);
                        ScaledComplexBall sqrtBase;
                        if (arithmetic ==
                                ScaledArithmeticStatus::
                                    Success)
                            arithmetic =
                                makePrincipalFunctionBall(
                                    inputCenter,
                                    PrincipalFunction::Sqrt,
                                    request.reference->
                                        certificationPrecision,
                                    sqrtBase);
                        if (arithmetic !=
                                ScaledArithmeticStatus::
                                    Success)
                            break;
                        if (!composeEntire(
                                local, node,
                                remainders[local],
                                &coefficient(local, 0),
                                ExpressionOracleOperation::
                                    Exp,
                                &sqrtBase, nullptr,
                                meromorphicDenominator.data(),
                                &halfLogRemainder)) {
                            arithmetic =
                                ScaledArithmeticStatus::
                                    Success;
                            break;
                        }
                        ScaledComplexBall referenceOutput;
                        arithmetic = makeScaledComplexValue(
                            node.output, node.outputDefect,
                            referenceOutput.value);
                        ScaledComplexBall baseDelta;
                        if (arithmetic ==
                                ScaledArithmeticStatus::
                                    Success)
                            arithmetic = subtractBall(
                                sqrtBase,
                                referenceOutput,
                                baseDelta);
                        if (arithmetic ==
                                ScaledArithmeticStatus::
                                    Success)
                            arithmetic = addBall(
                                coefficient(local, 0),
                                baseDelta,
                                coefficient(local, 0));
                        candidate.branchCompositionOperations +=
                            candidate.operations -
                            operationsBefore;
                        break;
                    }
                    case Op::Power: {
                        const uint64_t operationsBefore =
                            candidate.operations;
                        if ((node.flags &
                             OracleTraceHasLogarithmBase) == 0 ||
                            node.leftNode == UINT16_MAX ||
                            node.rightNode == UINT16_MAX) {
                            candidate.status =
                                ExpressionTaylorJetStatus::
                                    InvalidTape;
                            candidate.reason =
                                "power logarithm companion is missing";
                            break;
                        }
                        ScaledComplexBall logBase;
                        ScaledRealValue logRemainder;
                        if (!composePrincipalLogVariation(
                                node, sampleOffset,
                                node.leftNode, false,
                                meromorphicNumerator.data(),
                                logRemainder, logBase)) {
                            arithmetic =
                                ScaledArithmeticStatus::
                                    Success;
                            break;
                        }

                        ScaledComplexBall storedLogBase;
                        arithmetic = makeScaledComplexValue(
                            node.auxiliary,
                            node.auxiliaryDefect,
                            storedLogBase.value);
                        storedLogBase.radius =
                            node.auxiliaryError;
                        ScaledComplexBall logBaseDelta;
                        if (arithmetic ==
                                ScaledArithmeticStatus::
                                    Success)
                            arithmetic = subtractBall(
                                logBase,
                                storedLogBase,
                                logBaseDelta);
                        if (arithmetic ==
                                ScaledArithmeticStatus::
                                    Success)
                            arithmetic = addBall(
                                meromorphicNumerator[0],
                                logBaseDelta,
                                meromorphicNumerator[0]);
                        if (arithmetic ==
                                ScaledArithmeticStatus::
                                    Success)
                            arithmetic = addBall(
                                meromorphicNumerator[0],
                                storedLogBase,
                                meromorphicNumerator[0]);
                        if (arithmetic !=
                                ScaledArithmeticStatus::
                                    Success)
                            break;

                        ScaledComplexBall exponentCenter;
                        arithmetic = nodeBase(
                            sampleOffset, node.rightNode,
                            exponentCenter.value);
                        if (arithmetic ==
                                ScaledArithmeticStatus::
                                    Success)
                            arithmetic = addBall(
                                exponentCenter,
                                coefficient(
                                    node.rightNode, 0),
                                exponentCenter);
                        ScaledRealValue exponentRemainder =
                            remainders[node.rightNode];
                        if (arithmetic ==
                                ScaledArithmeticStatus::
                                    Success)
                            arithmetic = functionAddRemainder(
                                exponentRemainder,
                                exponentCenter.radius);
                        exponentCenter.radius = {};
                        std::fill(
                            meromorphicDenominator.begin(),
                            meromorphicDenominator.end(),
                            ScaledComplexBall{});
                        meromorphicDenominator[0].value =
                            exponentCenter.value;
                        for (int degree = 1;
                             degree <= order; ++degree)
                            meromorphicDenominator[
                                static_cast<size_t>(degree)] =
                                    coefficient(
                                        node.rightNode,
                                        degree);
                        ScaledRealValue productRemainder;
                        if (arithmetic ==
                                ScaledArithmeticStatus::
                                    Success)
                            arithmetic = multiplyPolynomials(
                                meromorphicNumerator.data(),
                                logRemainder,
                                meromorphicDenominator.data(),
                                exponentRemainder,
                                functionScratch,
                                productRemainder);

                        ScaledComplexValue productCenter;
                        if (arithmetic ==
                                ScaledArithmeticStatus::
                                    Success)
                            arithmetic = scaledMultiply(
                                exponentCenter.value,
                                logBase.value,
                                productCenter);
                        ScaledComplexBall productCenterPoint;
                        productCenterPoint.value =
                            productCenter;
                        if (arithmetic ==
                                ScaledArithmeticStatus::
                                    Success)
                            arithmetic = subtractBall(
                                functionScratch[0],
                                productCenterPoint,
                                functionScratch[0]);
                        for (int degree = 0;
                             arithmetic ==
                                 ScaledArithmeticStatus::
                                     Success &&
                             degree <= order; ++degree)
                            quotientCoefficients[
                                static_cast<size_t>(degree)] =
                                    functionScratch[
                                        static_cast<size_t>(
                                            degree)];
                        ScaledComplexBall exponentialBase;
                        if (arithmetic ==
                                ScaledArithmeticStatus::
                                    Success)
                            arithmetic =
                                makePrincipalFunctionBall(
                                    productCenter,
                                    PrincipalFunction::Exp,
                                    request.reference->
                                        certificationPrecision,
                                    exponentialBase);
                        if (arithmetic !=
                                ScaledArithmeticStatus::
                                    Success)
                            break;
                        if (!composeEntire(
                                local, node,
                                remainders[local],
                                &coefficient(local, 0),
                                ExpressionOracleOperation::
                                    Exp,
                                &exponentialBase, nullptr,
                                quotientCoefficients.data(),
                                &productRemainder)) {
                            arithmetic =
                                ScaledArithmeticStatus::
                                    Success;
                            break;
                        }
                        ScaledComplexBall referenceOutput;
                        arithmetic = makeScaledComplexValue(
                            node.output, node.outputDefect,
                            referenceOutput.value);
                        ScaledComplexBall baseDelta;
                        if (arithmetic ==
                                ScaledArithmeticStatus::
                                    Success)
                            arithmetic = subtractBall(
                                exponentialBase,
                                referenceOutput,
                                baseDelta);
                        if (arithmetic ==
                                ScaledArithmeticStatus::
                                    Success)
                            arithmetic = addBall(
                                coefficient(local, 0),
                                baseDelta,
                                coefficient(local, 0));
                        candidate.branchCompositionOperations +=
                            candidate.operations -
                            operationsBefore;
                        break;
                    }
                    case Op::Tan:
                    case Op::Tanh: {
                        if ((node.flags &
                             OracleTraceSingularPoint) != 0) {
                            candidate.status =
                                ExpressionTaylorJetStatus::
                                    PoleRejected;
                            candidate.reason =
                                "tangent reference is an exact pole";
                            candidate.poleRejected = true;
                            break;
                        }
                        if ((node.flags &
                             OracleTraceUndefined) != 0) {
                            candidate.status =
                                ExpressionTaylorJetStatus::
                                    Nonfinite;
                            candidate.reason =
                                "tangent reference is undefined";
                            break;
                        }
                        if ((node.flags &
                             OracleTraceHasDenominator) == 0 ||
                            (node.flags &
                             OracleTraceHasCompanion) == 0 ||
                            node.leftNode == UINT16_MAX) {
                            candidate.status =
                                ExpressionTaylorJetStatus::
                                    InvalidTape;
                            candidate.reason =
                                "tangent denominator companion is missing";
                            break;
                        }
                        if ((node.flags &
                             OracleTraceHasAsymptoticLogPhase) !=
                                0) {
                            candidate.status =
                                ExpressionTaylorJetStatus::
                                    AccuracyBudget;
                            candidate.reason =
                                "asymptotic tangent reference has no direct denominator jet";
                            break;
                        }
                        ScaledComplexBall referenceOutput;
                        ScaledComplexBall tangentBase;
                        ScaledComplexBall denominatorBase;
                        arithmetic = makeScaledComplexValue(
                            node.output, node.outputDefect,
                            referenceOutput.value);
                        tangentBase = referenceOutput;
                        tangentBase.radius =
                            node.outputError;
                        if (arithmetic ==
                                ScaledArithmeticStatus::
                                    Success)
                            arithmetic =
                                makeScaledComplexValue(
                                    node.auxiliary,
                                    node.auxiliaryDefect,
                                    denominatorBase.value);
                        denominatorBase.radius =
                            node.auxiliaryError;
                        if (arithmetic ==
                                ScaledArithmeticStatus::
                                    Success)
                            arithmetic =
                                inflateReferenceValue(
                                    tangentBase);
                        if (arithmetic ==
                                ScaledArithmeticStatus::
                                    Success)
                            arithmetic =
                                inflateReferenceValue(
                                    denominatorBase);
                        if (arithmetic !=
                                ScaledArithmeticStatus::
                                    Success)
                            break;

                        ScaledRealValue pointVariation =
                            denominatorBase.radius;
                        ScaledRealValue pointClearance;
                        if (arithmetic ==
                                ScaledArithmeticStatus::
                                    Success)
                            arithmetic =
                                lowerMaxComponentAfterVariation(
                                    denominatorBase.value,
                                    pointVariation,
                                    pointClearance);
                        ++candidate.operations;
                        ++candidate.reciprocalOperations;
                        if (arithmetic ==
                                ScaledArithmeticStatus::
                                    Singular) {
                            candidate.status =
                                ExpressionTaylorJetStatus::
                                    PoleRejected;
                            candidate.reason =
                                "tangent point denominator ball overlaps zero";
                            candidate.poleRejected = true;
                            arithmetic =
                                ScaledArithmeticStatus::
                                    Success;
                            break;
                        }
                        if (arithmetic !=
                                ScaledArithmeticStatus::
                                    Success)
                            break;
                        ScaledRealValue conservativePoint;
                        arithmetic = scaledDivideByDouble(
                            pointClearance, 2.0,
                            conservativePoint);
                        ++candidate.operations;
                        ++candidate.reciprocalOperations;
                        if (arithmetic !=
                                ScaledArithmeticStatus::
                                    Success ||
                            conservativePoint.isZero()) {
                            candidate.status =
                                ExpressionTaylorJetStatus::
                                    PoleRejected;
                            candidate.reason =
                                "tangent point denominator has no conservative margin";
                            candidate.poleRejected = true;
                            arithmetic =
                                ScaledArithmeticStatus::
                                    Success;
                            break;
                        }
                        ScaledComplexBall zero;
                        ScaledComplexBall one;
                        arithmetic = makeScaledRealValue(
                            1.0, one.value.re);
                        ++candidate.operations;
                        ++candidate.
                            functionSeriesOperations;
                        if (arithmetic !=
                                ScaledArithmeticStatus::
                                    Success)
                            break;
                        ScaledRealValue sineRemainder;
                        ScaledRealValue cosineRemainder;
                        const bool hyperbolic =
                            op == Op::Tanh;
                        if (!composeEntire(
                                local, node,
                                sineRemainder,
                                meromorphicNumerator.data(),
                                hyperbolic
                                    ? ExpressionOracleOperation::
                                          Sinh
                                    : ExpressionOracleOperation::
                                          Sin,
                                &zero,
                                &one) ||
                            !composeEntire(
                                local, node,
                                cosineRemainder,
                                meromorphicDenominator.data(),
                                hyperbolic
                                    ? ExpressionOracleOperation::
                                          Cosh
                                    : ExpressionOracleOperation::
                                          Cos,
                                &one,
                                &zero)) {
                            arithmetic =
                                ScaledArithmeticStatus::
                                    Success;
                            break;
                        }
                        arithmetic = addBall(
                            one,
                            meromorphicDenominator[0],
                            meromorphicDenominator[0]);
                        ScaledRealValue
                            numeratorRemainder;
                        if (arithmetic ==
                                ScaledArithmeticStatus::
                                    Success)
                            arithmetic = scalePolynomial(
                                meromorphicDenominator.
                                    data(),
                                cosineRemainder,
                                tangentBase,
                                functionScratch,
                                numeratorRemainder);
                        if (arithmetic ==
                                ScaledArithmeticStatus::
                                    Success)
                            arithmetic =
                                addPolynomialTerm(
                                    functionScratch,
                                    numeratorRemainder,
                                    meromorphicNumerator,
                                    sineRemainder,
                                    false);
                        ScaledRealValue scaledSineRemainder;
                        if (arithmetic ==
                                ScaledArithmeticStatus::
                                    Success)
                            arithmetic = scalePolynomial(
                                meromorphicNumerator.data(),
                                sineRemainder,
                                tangentBase,
                                functionFirst,
                                scaledSineRemainder);
                        for (int degree = 0;
                             arithmetic ==
                                 ScaledArithmeticStatus::
                                     Success &&
                             degree <= order; ++degree)
                            functionTerm[
                                static_cast<size_t>(degree)] =
                                    meromorphicDenominator[
                                        static_cast<size_t>(
                                            degree)];
                        ScaledRealValue denominatorRemainder =
                            cosineRemainder;
                        if (arithmetic ==
                                ScaledArithmeticStatus::
                                    Success)
                            arithmetic =
                                addPolynomialTerm(
                                    functionTerm,
                                    denominatorRemainder,
                                    functionFirst,
                                    scaledSineRemainder,
                                    !hyperbolic);
                        if (arithmetic !=
                                ScaledArithmeticStatus::
                                    Success)
                            break;
                        if (!dividePolynomials(
                                functionScratch.data(),
                                numeratorRemainder,
                                functionTerm.data(),
                                denominatorRemainder,
                                referenceOutput,
                                &coefficient(local, 0),
                                remainders[local],
                                &conservativePoint))
                            arithmetic =
                                ScaledArithmeticStatus::
                                    Success;
                        break;
                    }
                    default:
                        arithmetic =
                            ScaledArithmeticStatus::Singular;
                        break;
                    }
                    if (!candidate.reason.empty())
                        break;
                    if (arithmetic !=
                            ScaledArithmeticStatus::Success)
                        break;
                    arithmetic = addRemainder(
                        remainders[local],
                        node.outputError);
                    if (arithmetic ==
                            ScaledArithmeticStatus::Success)
                        arithmetic = addNodeRoundoff(
                            local, node, sampleOffset,
                            remainders[local]);
                    if (arithmetic !=
                            ScaledArithmeticStatus::Success)
                        break;
                    for (int degree = 0;
                         degree <= order; ++degree) {
                        if (certifyScaledMpfrExponentRange(
                                coefficient(
                                    local, degree)) !=
                                ScaledArithmeticStatus::
                                    Success) {
                            arithmetic =
                                ScaledArithmeticStatus::
                                    ExponentRange;
                            break;
                        }
                    }
                    if (arithmetic !=
                            ScaledArithmeticStatus::Success ||
                        certifyScaledMpfrExponentRange(
                            remainders[local]) !=
                            ScaledArithmeticStatus::Success) {
                        arithmetic =
                            ScaledArithmeticStatus::
                                ExponentRange;
                        break;
                    }
                }

                if (!candidate.reason.empty())
                    break;
                if (arithmetic !=
                        ScaledArithmeticStatus::Success) {
                    candidate.status =
                        statusForArithmetic(arithmetic);
                    candidate.reason =
                        "Taylor bytecode propagation failed";
                    break;
                }
                const uint16_t root = sample.rootNode;
                if (compareScaledNonnegative(
                        remainders[root],
                        accuracyBudget) > 0) {
                    candidate.status =
                        ExpressionTaylorJetStatus::
                            AccuracyBudget;
                    candidate.reason =
                        "Taylor truncation radius exceeds accuracy budget";
                    break;
                }

                ScaledComplexValue outputBase;
                arithmetic = makeScaledComplexValue(
                    sample.next, sample.rootDefect,
                    outputBase);
                ScaledRealValue margin;
                const InsideStatus inside =
                    arithmetic ==
                            ScaledArithmeticStatus::Success
                    ? certifyInside(
                          outputBase,
                          &coefficient(root, 0), order,
                          remainders[root], threshold,
                          margin, arithmetic)
                    : InsideStatus::Error;
                if (inside != InsideStatus::Inside) {
                    candidate.status =
                        inside == InsideStatus::Uncertain
                        ? ExpressionTaylorJetStatus::
                              BailoutUncertain
                        : statusForArithmetic(arithmetic);
                    candidate.reason =
                        "Taylor prefix reaches an uncertain or escaping state";
                    break;
                }

                if (static_cast<size_t>(landing + 1) ==
                        request.reference->samples.size()) {
                    for (int degree = 0;
                         degree <= order; ++degree)
                        state[
                            static_cast<size_t>(degree)] =
                                coefficient(root, degree);
                    stateRemainder = remainders[root];
                    ++landing;
                    candidate.margins.push_back(margin);
                    candidate.landingUsesSampleOutput = true;
                    break;
                }

                const ExpressionReferenceSample& nextSample =
                    request.reference->samples[
                        static_cast<size_t>(landing + 1)];
                ScaledComplexValue nextBaseValue;
                arithmetic = makeScaledComplexValue(
                    nextSample.z, nextSample.zDefect,
                    nextBaseValue);
                ScaledComplexBall rootBase;
                ScaledComplexBall nextBase;
                rootBase.value = outputBase;
                nextBase.value = nextBaseValue;
                ScaledComplexBall rebase;
                if (arithmetic ==
                        ScaledArithmeticStatus::Success)
                    arithmetic = subtractBall(
                        rootBase, nextBase, rebase);
                for (int degree = 0;
                     arithmetic ==
                         ScaledArithmeticStatus::Success &&
                     degree <= order; ++degree)
                    nextState[
                        static_cast<size_t>(degree)] =
                            coefficient(root, degree);
                if (arithmetic ==
                        ScaledArithmeticStatus::Success)
                    arithmetic = addBall(
                        nextState[0], rebase,
                        nextState[0]);
                if (arithmetic !=
                        ScaledArithmeticStatus::Success) {
                    candidate.status =
                        statusForArithmetic(arithmetic);
                    candidate.reason =
                        "Taylor compact rebase failed";
                    break;
                }
                state.swap(nextState);
                stateRemainder = remainders[root];
                ++landing;
                candidate.margins.push_back(margin);
            }

            candidate.landing = landing;
            candidate.coefficients = state;
            candidate.remainder = stateRemainder;
            if (candidate.reason.empty()) {
                candidate.status =
                    ExpressionTaylorJetStatus::Success;
            }
        } catch (const std::bad_alloc&) {
            candidate.status =
                ExpressionTaylorJetStatus::ResourceLimit;
            candidate.reason =
                "Taylor workspace allocation failed";
        } catch (const std::length_error&) {
            candidate.status =
                ExpressionTaylorJetStatus::ResourceLimit;
            candidate.reason =
                "Taylor workspace length overflow";
        }

        totalOperations += candidate.operations;
        lastStatus = candidate.status;
        lastReason = candidate.reason.empty()
            ? expressionTaylorJetStatusName(candidate.status)
            : candidate.reason;
        if (!haveBest ||
            candidate.landing > best.landing ||
            (candidate.landing == best.landing &&
             candidate.order < best.order)) {
            best = std::move(candidate);
            haveBest = true;
        }
        if (haveBest && best.landing >= maximumLanding)
            break;
        if (lastStatus ==
                ExpressionTaylorJetStatus::Cancelled ||
            lastStatus ==
                ExpressionTaylorJetStatus::ResourceLimit ||
            lastStatus ==
                ExpressionTaylorJetStatus::ExponentRange ||
            lastStatus ==
                ExpressionTaylorJetStatus::Nonfinite)
            break;
    }

    result.operationCount = totalOperations;
    result.memoryBytes = peakMemory;
    result.pixelParameter = request.pixelParameter;
    result.parameterScale = request.parameterScale;
    result.programSemanticHash =
        request.program->semanticHash();
    if (lastStatus ==
            ExpressionTaylorJetStatus::Cancelled)
        return finish(
            ExpressionTaylorJetStatus::Cancelled,
            lastReason, false);
    if (!haveBest)
        return finish(lastStatus, lastReason, false);

    size_t resultPeak = validationBytes;
    if (!checkedAddSize(
            resultPeak, best.coefficients.capacity(),
            sizeof(ScaledComplexBall)) ||
        !checkedAddSize(
            resultPeak, best.margins.capacity(),
            sizeof(ScaledRealValue)) ||
        !checkedAddSize(
            resultPeak, best.coefficients.size(),
            sizeof(ScaledComplexValue)) ||
        !checkedAddSize(
            resultPeak, best.coefficients.size(),
            sizeof(ScaledRealValue)))
        return finish(
            ExpressionTaylorJetStatus::ResourceLimit,
            "Taylor result memory calculation overflow",
            false);
    peakMemory = std::max(peakMemory, resultPeak);
    result.memoryBytes = peakMemory;
    if (request.memoryLimitBytes != 0 &&
        peakMemory > request.memoryLimitBytes)
        return finish(
            ExpressionTaylorJetStatus::ResourceLimit,
            "Taylor retained result exceeds memory policy",
            false);

    result.order = best.order;
    result.layout =
        ExpressionTaylorJetLayout::ComplexUnivariate;
    result.monomialCount =
        static_cast<size_t>(best.order) + 1;
    result.landingIteration = best.landing;
    result.landingSample =
        best.landingUsesSampleOutput
        ? static_cast<size_t>(best.landing - 1)
        : static_cast<size_t>(best.landing);
    result.landingUsesSampleOutput =
        best.landingUsesSampleOutput;
    result.maximumFunctionSeriesOrder =
        best.maximumFunctionSeriesOrder;
    result.functionSeriesCount =
        best.functionSeriesCount;
    result.functionSeriesOperationCount =
        best.functionSeriesOperations;
    result.maximumFunctionSeriesTail =
        best.maximumFunctionSeriesTail;
    result.maximumReciprocalOrder =
        best.maximumReciprocalOrder;
    result.reciprocalCount =
        best.reciprocalCount;
    result.reciprocalOperationCount =
        best.reciprocalOperations;
    result.minimumDenominatorClearance =
        best.minimumDenominatorClearance;
    result.maximumReciprocalTail =
        best.maximumReciprocalTail;
    result.poleRejected = best.poleRejected;
    result.maximumBranchSeriesOrder =
        best.maximumBranchSeriesOrder;
    result.branchCompositionCount =
        best.branchCompositionCount;
    result.branchCompositionOperationCount =
        best.branchCompositionOperations;
    result.maximumBranchSeriesTail =
        best.maximumBranchSeriesTail;
    result.minimumBranchCutClearance =
        best.minimumBranchCutClearance;
    result.minimumBranchZeroClearance =
        best.minimumBranchZeroClearance;
    result.branchRejected = best.branchRejected;
    result.remainderRadius = best.remainder;
    result.intermediateEscapeMargins =
        std::move(best.margins);
    result.coefficients.reserve(
        best.coefficients.size());
    result.coefficientRadii.reserve(
        best.coefficients.size());
    for (const ScaledComplexBall& coefficient :
         best.coefficients) {
        result.coefficients.push_back(
            coefficient.value);
        result.coefficientRadii.push_back(
            coefficient.radius);
    }
    if (best.landing < request.minimumLanding)
        return finish(
            best.status == ExpressionTaylorJetStatus::Success
                ? ExpressionTaylorJetStatus::NoCoverage
                : best.status,
            best.reason.empty()
                ? "Taylor prefix is shorter than minimum landing"
                : best.reason,
            false);

    result.valid = true;
    result.certified = true;
    return finish(
        ExpressionTaylorJetStatus::Success,
        best.reason, true);
}

bool ExpressionTaylorJetEvaluator::evaluate(
        const ExpressionTaylorJetResult& jet,
        const ScaledComplexBall& q,
        ExpressionTaylorJetEvaluation& result) {
    result = {};
    size_t expectedCoefficientCount = 0;
    if (jet.layout ==
            ExpressionTaylorJetLayout::RealBivariate) {
        if (!expressionTaylorBivariateMonomialCount(
                jet.order, expectedCoefficientCount))
            expectedCoefficientCount = 0;
    } else if (jet.order >= 0) {
        expectedCoefficientCount =
            static_cast<size_t>(jet.order) + 1;
    }
    if (!jet.valid || !jet.certified ||
        jet.order < 1 ||
        jet.coefficients.size() !=
            expectedCoefficientCount ||
        jet.monomialCount != expectedCoefficientCount ||
        jet.coefficientRadii.size() !=
            jet.coefficients.size() ||
        !validRadius(jet.remainderRadius) ||
        !expressionTaylorQInsideUnitDisk(q)) {
        result.status =
            ExpressionTaylorJetStatus::InvalidRequest;
        return false;
    }
    ScaledComplexBall value;
    if (jet.layout ==
            ExpressionTaylorJetLayout::RealBivariate) {
        ScaledComplexBall conjugateQ = q;
        ScaledArithmeticStatus arithmetic =
            scaledNegate(
                q.value.im,
                conjugateQ.value.im);
        ++result.operationCount;
        if (arithmetic !=
                ScaledArithmeticStatus::Success) {
            result.status =
                statusForArithmetic(arithmetic);
            return false;
        }
        bool haveValue = false;
        for (int conjugateDegree = jet.order;
             conjugateDegree >= 0;
             --conjugateDegree) {
            const int maximumQDegree =
                jet.order - conjugateDegree;
            size_t coefficientIndex = 0;
            if (!expressionTaylorBivariateIndex(
                    jet.order, maximumQDegree,
                    conjugateDegree,
                    coefficientIndex)) {
                result.status =
                    ExpressionTaylorJetStatus::
                        InvalidRequest;
                return false;
            }
            ScaledComplexBall row;
            row.value =
                jet.coefficients[coefficientIndex];
            row.radius =
                jet.coefficientRadii[coefficientIndex];
            for (int qDegree = maximumQDegree - 1;
                 qDegree >= 0; --qDegree) {
                ScaledComplexBall product;
                arithmetic = certifiedScaledMultiply(
                    row, q, product);
                ++result.operationCount;
                if (arithmetic !=
                        ScaledArithmeticStatus::
                            Success) {
                    result.status =
                        statusForArithmetic(arithmetic);
                    return false;
                }
                if (!expressionTaylorBivariateIndex(
                        jet.order, qDegree,
                        conjugateDegree,
                        coefficientIndex)) {
                    result.status =
                        ExpressionTaylorJetStatus::
                            InvalidRequest;
                    return false;
                }
                ScaledComplexBall coefficient;
                coefficient.value =
                    jet.coefficients[
                        coefficientIndex];
                coefficient.radius =
                    jet.coefficientRadii[
                        coefficientIndex];
                const bool productZero =
                    product.value.isZero() &&
                    product.radius.isZero();
                const bool coefficientZero =
                    coefficient.value.isZero() &&
                    coefficient.radius.isZero();
                if (productZero && coefficientZero) {
                    row = product;
                    row.value.re.mantissa =
                        std::signbit(product.value.re.mantissa) ||
                        std::signbit(coefficient.value.re.mantissa)
                        ? -0.0 : 0.0;
                    row.value.im.mantissa =
                        std::signbit(product.value.im.mantissa) ||
                        std::signbit(coefficient.value.im.mantissa)
                        ? -0.0 : 0.0;
                } else if (productZero) {
                    row = coefficient;
                } else if (coefficientZero) {
                    row = product;
                } else {
                    arithmetic = certifiedScaledAdd(
                        product, coefficient, row);
                    ++result.operationCount;
                    if (arithmetic !=
                            ScaledArithmeticStatus::
                                Success) {
                        result.status =
                            statusForArithmetic(arithmetic);
                        return false;
                    }
                }
            }
            if (!haveValue) {
                value = row;
                haveValue = true;
            } else {
                ScaledComplexBall product;
                arithmetic = certifiedScaledMultiply(
                    value, conjugateQ, product);
                ++result.operationCount;
                if (arithmetic !=
                        ScaledArithmeticStatus::
                            Success) {
                    result.status =
                        statusForArithmetic(arithmetic);
                    return false;
                }
                const bool productZero =
                    product.value.isZero() &&
                    product.radius.isZero();
                const bool rowZero =
                    row.value.isZero() &&
                    row.radius.isZero();
                if (productZero && rowZero) {
                    value = product;
                    value.value.re.mantissa =
                        std::signbit(product.value.re.mantissa) ||
                        std::signbit(row.value.re.mantissa)
                        ? -0.0 : 0.0;
                    value.value.im.mantissa =
                        std::signbit(product.value.im.mantissa) ||
                        std::signbit(row.value.im.mantissa)
                        ? -0.0 : 0.0;
                } else if (productZero) {
                    value = row;
                } else if (rowZero) {
                    value = product;
                } else {
                    arithmetic = certifiedScaledAdd(
                        product, row, value);
                    ++result.operationCount;
                    if (arithmetic !=
                            ScaledArithmeticStatus::
                                Success) {
                        result.status =
                            statusForArithmetic(arithmetic);
                        return false;
                    }
                }
            }
        }
    } else {
        value.value = jet.coefficients.back();
        value.radius = jet.coefficientRadii.back();
        for (int degree = jet.order - 1;
             degree >= 0; --degree) {
            ScaledComplexBall product;
            ScaledArithmeticStatus arithmetic =
                certifiedScaledMultiply(
                    value, q, product);
            ++result.operationCount;
            if (arithmetic !=
                    ScaledArithmeticStatus::Success) {
                result.status =
                    statusForArithmetic(arithmetic);
                return false;
            }
            ScaledComplexBall coefficient;
            coefficient.value =
                jet.coefficients[
                    static_cast<size_t>(degree)];
            coefficient.radius =
                jet.coefficientRadii[
                    static_cast<size_t>(degree)];
            arithmetic = certifiedScaledAdd(
                product, coefficient, value);
            ++result.operationCount;
            if (arithmetic !=
                    ScaledArithmeticStatus::Success) {
                result.status =
                    statusForArithmetic(arithmetic);
                return false;
            }
        }
    }
    if (addUp(
            value.radius, jet.remainderRadius) !=
            ScaledArithmeticStatus::Success ||
        certifyScaledMpfrExponentRange(value) !=
            ScaledArithmeticStatus::Success) {
        result.status =
            ExpressionTaylorJetStatus::ExponentRange;
        return false;
    }
    result.residual = value;
    result.valid = true;
    result.status = ExpressionTaylorJetStatus::Success;
    return true;
}

} // namespace formula
