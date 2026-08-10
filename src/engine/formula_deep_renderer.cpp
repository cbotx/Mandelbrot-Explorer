#include "formula_deep_renderer.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <climits>
#include <cmath>
#include <limits>
#include <memory>
#include <mutex>
#include <new>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace formula {

namespace {

using Clock = std::chrono::steady_clock;
constexpr mpfr_rnd_t RND = MPFR_RNDN;
constexpr uint8_t NoFallbackReason =
    static_cast<uint8_t>(ExpressionDeepFallbackReason::Count);

double secondsSince(const Clock::time_point& start) {
    return std::chrono::duration<double>(
        Clock::now() - start).count();
}

bool finiteComplex(Complex value) {
    return std::isfinite(value.real()) &&
           std::isfinite(value.imag());
}

bool decimalSignificandIsNonzero(const std::string& text) {
    size_t end = text.find_first_of("eE");
    if (end == std::string::npos) end = text.size();
    for (size_t i = 0; i < end; ++i)
        if (text[i] >= '1' && text[i] <= '9')
            return true;
    return false;
}

size_t decimalDigits(const std::string& text) {
    size_t digits = 0;
    size_t end = text.find_first_of("eE");
    if (end == std::string::npos) end = text.size();
    for (size_t i = 0; i < end; ++i)
        if (text[i] >= '0' && text[i] <= '9')
            ++digits;
    return digits;
}

uint64_t decimalBits(const std::string& text) {
    return static_cast<uint64_t>(
        std::ceil(decimalDigits(text) *
                  3.32192809488736234787)) + 2;
}

uint64_t mpfBits(mpf_srcptr value) {
    if (!value) return 0;
    const uint64_t limbs =
        static_cast<uint64_t>(mpf_size(value));
    if (limbs >
            std::numeric_limits<uint64_t>::max() /
                static_cast<uint64_t>(GMP_NUMB_BITS))
        return std::numeric_limits<uint64_t>::max();
    return limbs * static_cast<uint64_t>(GMP_NUMB_BITS);
}

bool validCenterRepresentation(
        const ExpressionReferenceExactInput& center) {
    if (center.usesMpf())
        return center.realMpf && center.imaginaryMpf &&
               !center.usesDecimal();
    return !center.realDecimal.empty() &&
           !center.usesMpf();
}

bool validScaleRepresentation(
        const ExpressionDeepExactReal& scale) {
    return scale.usesMpf()
        ? !scale.usesDecimal()
        : scale.usesDecimal();
}

bool setDecimal(mpfr_ptr output, const std::string& text) {
    const mpfr_flags_t saved = mpfr_flags_save();
    mpfr_flags_clear(MPFR_FLAGS_ALL);
    const int status =
        mpfr_set_str(output, text.c_str(), 10, RND);
    const mpfr_flags_t rangeFlags = mpfr_flags_test(
        MPFR_FLAGS_UNDERFLOW |
        MPFR_FLAGS_OVERFLOW |
        MPFR_FLAGS_ERANGE);
    const bool okay =
        status == 0 && rangeFlags == 0 &&
        mpfr_number_p(output) &&
        !(mpfr_zero_p(output) &&
          decimalSignificandIsNonzero(text));
    mpfr_flags_restore(saved, MPFR_FLAGS_ALL);
    return okay;
}

bool setExactMpf(mpfr_ptr output, mpf_srcptr input) {
    if (!input) return false;
    const mpfr_flags_t saved = mpfr_flags_save();
    mpfr_flags_clear(MPFR_FLAGS_ALL);
    mpfr_set_f(output, input, RND);
    const mpfr_flags_t rangeFlags = mpfr_flags_test(
        MPFR_FLAGS_UNDERFLOW |
        MPFR_FLAGS_OVERFLOW |
        MPFR_FLAGS_ERANGE);
    const bool okay =
        rangeFlags == 0 &&
        mpfr_number_p(output) &&
        !(mpf_sgn(input) != 0 && mpfr_zero_p(output));
    mpfr_flags_restore(saved, MPFR_FLAGS_ALL);
    return okay;
}

bool setExactReal(
        mpfr_ptr output, const ExpressionDeepExactReal& input) {
    if (input.usesMpf()) {
        if (input.usesDecimal()) return false;
        return setExactMpf(output, input.mpf);
    }
    return input.usesDecimal() &&
           setDecimal(output, input.decimal);
}

bool setExactCenter(
        MpfrComplex& output,
        const ExpressionReferenceExactInput& input) {
    if (input.usesMpf()) {
        if (!input.realMpf || !input.imaginaryMpf ||
            input.usesDecimal())
            return false;
        if (!setExactMpf(output.re, input.realMpf) ||
            !setExactMpf(output.im, input.imaginaryMpf))
            return false;
    } else {
        if (input.realDecimal.empty()) return false;
        const std::string imaginary =
            input.imaginaryDecimal.empty()
                ? "0" : input.imaginaryDecimal;
        if (!setDecimal(output.re, input.realDecimal) ||
            !setDecimal(output.im, imaginary))
            return false;
    }
    return mpfr_number_p(output.re) &&
           mpfr_number_p(output.im);
}

class ExactGeometry {
public:
    explicit ExactGeometry(mpfr_prec_t precision)
        : center(precision) {
        mpfr_inits2(
            precision, scale, dxHalf, dyHalf, temporary,
            (mpfr_ptr)0);
    }

    ~ExactGeometry() {
        mpfr_clears(
            scale, dxHalf, dyHalf, temporary,
            (mpfr_ptr)0);
    }

    bool initialize(const ExpressionDeepRenderRequest& request) {
        if (!setExactCenter(center, request.center) ||
            !setExactReal(scale, request.scale) ||
            mpfr_sgn(scale) <= 0)
            return false;

        // dx/2 = 2/(scale*(width-1)).
        mpfr_mul_ui(
            temporary, scale,
            static_cast<unsigned long>(request.width - 1),
            RND);
        mpfr_ui_div(dxHalf, 2, temporary, RND);

        // dy/2 = 2*height/(scale*width*(height-1)).
        mpfr_mul_ui(
            temporary, scale,
            static_cast<unsigned long>(request.width),
            RND);
        mpfr_mul_ui(
            temporary, temporary,
            static_cast<unsigned long>(request.height - 1),
            RND);
        mpfr_ui_div(
            dyHalf,
            static_cast<unsigned long>(request.height),
            temporary, RND);
        mpfr_mul_ui(dyHalf, dyHalf, 2, RND);
        return mpfr_number_p(dxHalf) &&
               mpfr_number_p(dyHalf) &&
               !mpfr_zero_p(dxHalf) &&
               !mpfr_zero_p(dyHalf);
    }

    bool coordinate(
            int x, int y,
            const ExpressionDeepRenderRequest& request,
            MpfrComplex& output) {
        const long centeredX =
            static_cast<long>(2LL * x -
                              (request.width - 1LL));
        const long centeredY =
            static_cast<long>(2LL * y -
                              (request.height - 1LL));
        if (centeredX == 0) {
            mpfr_set(output.re, center.re, RND);
        } else {
            mpfr_mul_si(output.re, dxHalf, centeredX, RND);
            mpfr_add(output.re, output.re, center.re, RND);
        }
        if (centeredY == 0) {
            mpfr_set(output.im, center.im, RND);
        } else {
            mpfr_mul_si(output.im, dyHalf, centeredY, RND);
            mpfr_add(output.im, output.im, center.im, RND);
        }
        return mpfr_number_p(output.re) &&
               mpfr_number_p(output.im);
    }

    MpfrComplex center;
    mpfr_t scale;
    mpfr_t dxHalf;
    mpfr_t dyHalf;
    mpfr_t temporary;
};

struct FallbackWorkerWorkspace {
    FallbackWorkerWorkspace(
            mpfr_prec_t precision,
            const ExpressionDeepRenderRequest& request)
        : geometry(precision),
          geometryReady(geometry.initialize(request)),
          context(precision),
          pixel(precision),
          next(precision),
          magnitudeStorage(precision) {}

    ExactGeometry geometry;
    bool geometryReady = false;
    ExpressionOracleContext context;
    MpfrComplex pixel;
    MpfrComplex next;
    MpfrComplex magnitudeStorage;
};

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

uint64_t ceilLog2(uint64_t value) {
    uint64_t result = 0;
    uint64_t power = 1;
    while (power < value && result < 63) {
        power <<= 1;
        ++result;
    }
    return result;
}

bool selectAutomaticViewBits(
        const ExpressionDeepRenderRequest& request,
        ExpressionReferencePrecisionPolicy& policy,
        std::string& error) {
    const uint64_t centerBits = request.center.usesMpf()
        ? std::max(
            mpfBits(request.center.realMpf),
            mpfBits(request.center.imaginaryMpf))
        : std::max(
            decimalBits(request.center.realDecimal),
            decimalBits(request.center.imaginaryDecimal));
    const uint64_t scaleBits = request.scale.usesMpf()
        ? mpfBits(request.scale.mpf)
        : decimalBits(request.scale.decimal);
    const uint64_t inputBits =
        std::max(centerBits, scaleBits);
    if (inputBits >
            static_cast<uint64_t>(MPFR_PREC_MAX)) {
        error = "exact input precision exceeds MPFR";
        return false;
    }
    policy.requestedBits = std::max(
        policy.requestedBits,
        static_cast<mpfr_prec_t>(inputBits));

    const mpfr_prec_t probePrecision = 256;
    ExactGeometry geometry(probePrecision);
    if (!geometry.initialize(request)) {
        error = "failed to parse exact center or positive scale";
        return false;
    }

    uint64_t required = 53;
    const mpfr_exp_t scaleExponent =
        mpfr_get_exp(geometry.scale);
    if (scaleExponent > 0) {
        required = std::max<uint64_t>(
            required,
            static_cast<uint64_t>(scaleExponent) +
                ceilLog2(static_cast<uint64_t>(
                    std::max(request.width, request.height))) +
                8);
    }
    auto coverAddition = [&](mpfr_srcptr center,
                             mpfr_srcptr step) {
        if (mpfr_zero_p(center) || mpfr_zero_p(step))
            return;
        const mpfr_exp_t difference =
            mpfr_get_exp(center) - mpfr_get_exp(step);
        if (difference > 0)
            required = std::max<uint64_t>(
                required,
                static_cast<uint64_t>(difference) + 8);
    };
    coverAddition(geometry.center.re, geometry.dxHalf);
    coverAddition(geometry.center.im, geometry.dyHalf);
    if (required >
            static_cast<uint64_t>(MPFR_PREC_MAX)) {
        error = "view precision exceeds MPFR";
        return false;
    }
    policy.viewBits = std::max(
        policy.viewBits,
        static_cast<mpfr_prec_t>(required));
    return true;
}

bool selectDeepPrecision(
        const ExpressionReferencePrecisionPolicy& policy,
        mpfr_prec_t& precision,
        std::string& error) {
    const uint64_t base = std::max<uint64_t>({
        53,
        static_cast<uint64_t>(policy.requestedBits),
        static_cast<uint64_t>(policy.viewBits),
        static_cast<uint64_t>(policy.minimumBits)
    });
    const uint64_t guard =
        static_cast<uint64_t>(policy.guardBits);
    if (base >
            std::numeric_limits<uint64_t>::max() - guard) {
        error = "MPFR precision calculation overflow";
        return false;
    }
    const uint64_t selected = base + guard;
    if (selected >
            static_cast<uint64_t>(policy.maximumBits) ||
        selected >
            static_cast<uint64_t>(
                ExpressionReferencePrecisionPolicy::
                    ApplicationMaximumBits) ||
        selected >
            static_cast<uint64_t>(MPFR_PREC_MAX)) {
        error = "required MPFR precision exceeds policy maximum";
        return false;
    }
    precision = static_cast<mpfr_prec_t>(selected);
    return true;
}

void configureFixed(
        const ExpressionContext& fixed,
        ExpressionOracleContext& context) {
    context.z.set(fixed.z.real(), fixed.z.imag());
    context.c.set(fixed.c.real(), fixed.c.imag());
    context.z0.set(fixed.z0.real(), fixed.z0.imag());
    for (size_t i = 0; i < fixed.parameters.size(); ++i) {
        context.parameters[i].set(
            fixed.parameters[i].real(),
            fixed.parameters[i].imag());
    }
    context.iteration = fixed.iteration;
}

ExpressionDeepFallbackReason reasonForCapability(
        ExpressionScaledResidualCapability capability) {
    switch (capability) {
    case ExpressionScaledResidualCapability::
            CertifiedEntireCandidate:
    case ExpressionScaledResidualCapability::
            CertifiedMeromorphicCandidate:
    case ExpressionScaledResidualCapability::
            CertifiedBranchCandidate:
    case ExpressionScaledResidualCapability::
            CertifiedRealCandidate:
    case ExpressionScaledResidualCapability::
            CertifiedPiecewiseCandidate:
        return ExpressionDeepFallbackReason::
            CertificationFailure;
    case ExpressionScaledResidualCapability::UncertifiedSeries:
        return ExpressionDeepFallbackReason::UncertifiedSeries;
    case ExpressionScaledResidualCapability::BranchSensitive:
        return ExpressionDeepFallbackReason::BranchSensitive;
    case ExpressionScaledResidualCapability::Unsupported:
        return ExpressionDeepFallbackReason::UnsupportedOperation;
    case ExpressionScaledResidualCapability::ExactCenteredArithmetic:
        break;
    }
    return ExpressionDeepFallbackReason::InvalidTape;
}

bool certifiedTaylorCapability(
        ExpressionScaledResidualCapability capability) {
    return capability ==
               ExpressionScaledResidualCapability::
                   CertifiedEntireCandidate ||
           capability ==
               ExpressionScaledResidualCapability::
                   CertifiedMeromorphicCandidate ||
           capability ==
               ExpressionScaledResidualCapability::
                   CertifiedBranchCandidate ||
           capability ==
               ExpressionScaledResidualCapability::
                   CertifiedRealCandidate ||
           capability ==
               ExpressionScaledResidualCapability::
                   CertifiedPiecewiseCandidate;
}

bool certifiedReferenceCapability(
        ExpressionScaledResidualCapability capability) {
    return capability ==
               ExpressionScaledResidualCapability::
                   ExactCenteredArithmetic ||
           certifiedTaylorCapability(capability);
}

ExpressionDeepFallbackReason reasonForResidualStatus(
        ExpressionScaledResidualStatus status) {
    switch (status) {
    case ExpressionScaledResidualStatus::BranchUncertain:
        return ExpressionDeepFallbackReason::BranchSensitive;
    case ExpressionScaledResidualStatus::Singular:
        return ExpressionDeepFallbackReason::Singular;
    case ExpressionScaledResidualStatus::Unsupported:
        return ExpressionDeepFallbackReason::UnsupportedOperation;
    case ExpressionScaledResidualStatus::Nonfinite:
        return ExpressionDeepFallbackReason::Nonfinite;
    case ExpressionScaledResidualStatus::ExponentRange:
        return ExpressionDeepFallbackReason::ExponentRange;
    case ExpressionScaledResidualStatus::InvalidTape:
    case ExpressionScaledResidualStatus::InvalidInput:
        return ExpressionDeepFallbackReason::InvalidTape;
    case ExpressionScaledResidualStatus::Success:
        break;
    }
    return ExpressionDeepFallbackReason::InvalidTape;
}

ExpressionDeepFallbackReason reasonForArithmeticStatus(
        ScaledArithmeticStatus status) {
    switch (status) {
    case ScaledArithmeticStatus::Nonfinite:
        return ExpressionDeepFallbackReason::Nonfinite;
    case ScaledArithmeticStatus::ExponentRange:
        return ExpressionDeepFallbackReason::ExponentRange;
    case ScaledArithmeticStatus::Singular:
        return ExpressionDeepFallbackReason::Singular;
    case ScaledArithmeticStatus::Success:
        break;
    }
    return ExpressionDeepFallbackReason::ReconstructionFailure;
}

enum class BailoutDecision : uint8_t {
    Inside,
    Outside,
    Uncertain,
    Error
};

struct BailoutThreshold {
    ScaledRealValue midpoint;
    ScaledRealValue error;
};

bool makeBailoutThreshold(
        double bailout, BailoutThreshold& threshold) {
    ScaledComplexValue radius;
    if (makeScaledRealValue(bailout, radius.re) !=
            ScaledArithmeticStatus::Success ||
        certifiedScaledNormSquared(
            radius, {}, threshold.midpoint,
            threshold.error) !=
                ScaledArithmeticStatus::Success)
        return false;
    return true;
}

BailoutDecision decideBailout(
        const ScaledComplexValue& value,
        const ScaledRealValue& radius,
        const BailoutThreshold& threshold,
        ScaledArithmeticStatus& arithmeticStatus) {
    ScaledComplexBall state{ value, radius };
    arithmeticStatus =
        certifyScaledMpfrExponentRange(state);
    if (arithmeticStatus != ScaledArithmeticStatus::Success)
        return BailoutDecision::Error;
    ScaledRealValue norm, error;
    arithmeticStatus = certifiedScaledNormSquared(
        value, radius, norm, error);
    if (arithmeticStatus != ScaledArithmeticStatus::Success)
        return BailoutDecision::Error;
    arithmeticStatus =
        certifyScaledMpfrExponentRange(norm);
    if (arithmeticStatus != ScaledArithmeticStatus::Success)
        return BailoutDecision::Error;
    arithmeticStatus =
        certifyScaledMpfrExponentRange(error);
    if (arithmeticStatus != ScaledArithmeticStatus::Success)
        return BailoutDecision::Error;
    ScaledRealValue thresholdUpper, outsideMargin;
    arithmeticStatus = scaledAddUp(
        threshold.midpoint, threshold.error,
        thresholdUpper);
    if (arithmeticStatus != ScaledArithmeticStatus::Success)
        return BailoutDecision::Error;
    arithmeticStatus =
        certifyScaledMpfrExponentRange(thresholdUpper);
    if (arithmeticStatus != ScaledArithmeticStatus::Success)
        return BailoutDecision::Error;
    arithmeticStatus = scaledAddUp(
        thresholdUpper, error, outsideMargin);
    if (arithmeticStatus != ScaledArithmeticStatus::Success)
        return BailoutDecision::Error;
    arithmeticStatus =
        certifyScaledMpfrExponentRange(outsideMargin);
    if (arithmeticStatus != ScaledArithmeticStatus::Success)
        return BailoutDecision::Error;
    if (compareScaledNonnegative(norm, outsideMargin) > 0)
        return BailoutDecision::Outside;
    ScaledRealValue normUpper, insideMargin;
    arithmeticStatus = scaledAddUp(
        norm, error, normUpper);
    if (arithmeticStatus != ScaledArithmeticStatus::Success)
        return BailoutDecision::Error;
    arithmeticStatus =
        certifyScaledMpfrExponentRange(normUpper);
    if (arithmeticStatus != ScaledArithmeticStatus::Success)
        return BailoutDecision::Error;
    arithmeticStatus = scaledAddUp(
        normUpper, threshold.error, insideMargin);
    if (arithmeticStatus != ScaledArithmeticStatus::Success)
        return BailoutDecision::Error;
    arithmeticStatus =
        certifyScaledMpfrExponentRange(insideMargin);
    if (arithmeticStatus != ScaledArithmeticStatus::Success)
        return BailoutDecision::Error;
    if (compareScaledNonnegative(
            insideMargin, threshold.midpoint) <= 0)
        return BailoutDecision::Inside;
    return BailoutDecision::Uncertain;
}

struct ScaledOffset {
    ScaledRealValue value;
    ScaledRealValue error;
};

ScaledRealValue absoluteScaled(ScaledRealValue value) {
    value.mantissa = std::abs(value.mantissa);
    return value;
}

bool componentDiscrepancy(
        mpfr_srcptr exact,
        const ScaledRealValue& base,
        const ScaledRealValue& offset,
        mpfr_prec_t precision,
        ScaledRealValue& error) {
    mpfr_t baseValue, offsetValue, lower, upper;
    mpfr_t differenceLower, differenceUpper, maximum;
    mpfr_inits2(
        precision, baseValue, offsetValue, lower, upper,
        differenceLower, differenceUpper, maximum,
        (mpfr_ptr)0);
    const bool reconstructed =
        setMpfrFromScaledValue(baseValue, base) &&
        setMpfrFromScaledValue(offsetValue, offset);
    if (reconstructed) {
        mpfr_add(
            lower, baseValue, offsetValue, MPFR_RNDD);
        mpfr_add(
            upper, baseValue, offsetValue, MPFR_RNDU);
        mpfr_sub(
            differenceLower, exact, upper, MPFR_RNDD);
        mpfr_sub(
            differenceUpper, exact, lower, MPFR_RNDU);
        mpfr_abs(
            differenceLower, differenceLower, MPFR_RNDU);
        mpfr_abs(
            differenceUpper, differenceUpper, MPFR_RNDU);
        if (mpfr_cmp(
                differenceLower, differenceUpper) >= 0)
            mpfr_set(maximum, differenceLower, MPFR_RNDU);
        else
            mpfr_set(maximum, differenceUpper, MPFR_RNDU);
    }
    const bool okay = reconstructed &&
        makeScaledNonnegativeUpward(maximum, error) ==
            ScaledArithmeticStatus::Success;
    mpfr_clears(
        baseValue, offsetValue, lower, upper,
        differenceLower, differenceUpper, maximum,
        (mpfr_ptr)0);
    return okay;
}

bool inflateRadius(
        ScaledRealValue& radius, int bits) {
    if (bits == 0 || radius.isZero())
        return true;
    if (radius.exponent >
            std::numeric_limits<int64_t>::max() - bits)
        return false;
    radius.exponent += bits;
    return true;
}

bool inflateReferenceErrors(
        ExpressionReferenceOrbitResult& reference,
        int bits) {
    if (!inflateRadius(reference.cError, bits) ||
        !inflateRadius(reference.z0Error, bits) ||
        !inflateRadius(reference.pixelError, bits) ||
        !inflateRadius(reference.initialZError, bits))
        return false;
    for (ExpressionReferenceSample& sample :
         reference.samples)
        if (!inflateRadius(sample.zError, bits) ||
            !inflateRadius(sample.nextError, bits) ||
            !inflateRadius(sample.rootError, bits))
            return false;
    for (ExpressionReferenceTapeNode& node :
         reference.tape)
        if (!inflateRadius(node.outputError, bits) ||
            !inflateRadius(node.auxiliaryError, bits))
            return false;
    return true;
}

struct FastPixelResult {
    bool decided = false;
    float output = ExpressionDeepEmptyPixel;
    ExpressionDeepFallbackReason reason =
        ExpressionDeepFallbackReason::InvalidTape;
    uint64_t iterations = 0;
};

ExpressionDeepRenderStatus mapReferenceStatus(
        ExpressionReferenceBuildStatus status) {
    switch (status) {
    case ExpressionReferenceBuildStatus::ProgramMismatch:
        return ExpressionDeepRenderStatus::ProgramMismatch;
    case ExpressionReferenceBuildStatus::PrecisionOutOfRange:
        return ExpressionDeepRenderStatus::PrecisionOutOfRange;
    case ExpressionReferenceBuildStatus::ResourceLimit:
        return ExpressionDeepRenderStatus::ResourceLimit;
    case ExpressionReferenceBuildStatus::InvalidRequest:
    case ExpressionReferenceBuildStatus::InputParseError:
        return ExpressionDeepRenderStatus::InvalidRequest;
    case ExpressionReferenceBuildStatus::Success:
    case ExpressionReferenceBuildStatus::UnsupportedProgram:
    case ExpressionReferenceBuildStatus::CompactionOutOfRange:
        return ExpressionDeepRenderStatus::ReferenceFailure;
    }
    return ExpressionDeepRenderStatus::ReferenceFailure;
}

} // namespace

const char* expressionDeepRenderStatusName(
        ExpressionDeepRenderStatus status) {
    switch (status) {
    case ExpressionDeepRenderStatus::Success: return "success";
    case ExpressionDeepRenderStatus::Cancelled: return "cancelled";
    case ExpressionDeepRenderStatus::InvalidRequest:
        return "invalid-request";
    case ExpressionDeepRenderStatus::ProgramMismatch:
        return "program-mismatch";
    case ExpressionDeepRenderStatus::PrecisionOutOfRange:
        return "precision-out-of-range";
    case ExpressionDeepRenderStatus::ResourceLimit:
        return "resource-limit";
    case ExpressionDeepRenderStatus::ReferenceFailure:
        return "reference-failure";
    case ExpressionDeepRenderStatus::UndefinedPixel:
        return "undefined-pixel";
    case ExpressionDeepRenderStatus::InternalError:
        return "internal-error";
    }
    return "internal-error";
}

const char* expressionDeepFallbackReasonName(
        ExpressionDeepFallbackReason reason) {
    switch (reason) {
    case ExpressionDeepFallbackReason::UncertifiedSeries:
        return "uncertified-series";
    case ExpressionDeepFallbackReason::BranchSensitive:
        return "branch-sensitive";
    case ExpressionDeepFallbackReason::UnsupportedOperation:
        return "unsupported-operation";
    case ExpressionDeepFallbackReason::Singular: return "singular";
    case ExpressionDeepFallbackReason::Nonfinite: return "nonfinite";
    case ExpressionDeepFallbackReason::ExponentRange:
        return "exponent-range";
    case ExpressionDeepFallbackReason::InvalidTape:
        return "invalid-tape";
    case ExpressionDeepFallbackReason::ReferenceExhausted:
        return "reference-exhausted";
    case ExpressionDeepFallbackReason::CertificationFailure:
        return "certification-failure";
    case ExpressionDeepFallbackReason::BailoutUncertain:
        return "bailout-uncertain";
    case ExpressionDeepFallbackReason::ReconstructionFailure:
        return "reconstruction-failure";
    case ExpressionDeepFallbackReason::Count: break;
    }
    return "invalid";
}

bool renderExpressionDeepFrame(
        const ExpressionDeepRenderRequest& request,
        ExpressionDeepRenderResult& result) {
    result = {};
    auto fail = [&](ExpressionDeepRenderStatus status,
                    const std::string& error) {
        result.status = status;
        result.error = error;
        result.success = false;
        result.cancelled =
            status == ExpressionDeepRenderStatus::Cancelled;
        return false;
    };

    try {
        if (!request.runtimeProgram ||
            !request.runtimeProgram->valid())
            return fail(
                ExpressionDeepRenderStatus::InvalidRequest,
                "runtime program is missing or invalid");
        if (request.canonicalProgram &&
            !request.canonicalProgram->valid())
            return fail(
                ExpressionDeepRenderStatus::InvalidRequest,
                "canonical program is invalid");
        if (!validCenterRepresentation(request.center) ||
            !validScaleRepresentation(request.scale))
            return fail(
                ExpressionDeepRenderStatus::InvalidRequest,
                "exact center or scale representation is invalid");
        if (request.pixelParameter != FormulaParameter::C &&
            request.pixelParameter != FormulaParameter::InitialZ)
            return fail(
                ExpressionDeepRenderStatus::InvalidRequest,
                "pixel binding must be c or z0");
        if (request.width < 2 || request.height < 2 ||
            request.width > (LONG_MAX + 1LL) / 2 ||
            request.height > (LONG_MAX + 1LL) / 2)
            return fail(
                ExpressionDeepRenderStatus::InvalidRequest,
                "frame dimensions are invalid");
        if (request.maxIterations < 1 ||
            request.maxIterations > (1 << 24))
            return fail(
                ExpressionDeepRenderStatus::InvalidRequest,
                "iteration count cannot be represented exactly in float output");
        if (!(request.bailout > 0.0) ||
            !std::isfinite(request.bailout))
            return fail(
                ExpressionDeepRenderStatus::InvalidRequest,
                "bailout must be finite and positive");
        if (!request.output)
            return fail(
                ExpressionDeepRenderStatus::InvalidRequest,
                "output buffer is null");
        if (request.threading.tileWidth < 1 ||
            request.threading.tileHeight < 1 ||
            request.threading.threads < 0 ||
            request.threading.threads > 1024)
            return fail(
                ExpressionDeepRenderStatus::InvalidRequest,
                "thread or tile policy is invalid");
        if (request.memory.fallbackGuardBits < 0)
            return fail(
                ExpressionDeepRenderStatus::InvalidRequest,
                "fallback guard precision is invalid");
        if (request.taylor.minimumLanding < 1 ||
            request.taylor.minimumOrder < 8 ||
            request.taylor.maximumOrder > 20 ||
            request.taylor.maximumBivariateOrder < 8 ||
            request.taylor.maximumBivariateOrder >
                ExpressionTaylorMaximumBivariateOrder ||
            request.taylor.minimumOrder >
                request.taylor.order ||
            request.taylor.order >
                request.taylor.maximumOrder ||
            request.taylor.maximumCompositionOrder < 1 ||
            request.taylor.maximumCompositionOrder > 24 ||
            request.taylor.maximumCandidateIteration < 0 ||
            !(request.taylor.accuracyBudget > 0.0) ||
            !std::isfinite(
                request.taylor.accuracyBudget))
            return fail(
                ExpressionDeepRenderStatus::InvalidRequest,
                "Taylor policy is invalid");
        if (request.verificationErrorInflationBits < 0 ||
            request.verificationErrorInflationBits > 4096)
            return fail(
                ExpressionDeepRenderStatus::InvalidRequest,
                "verification error inflation is invalid");
        if (request.verificationFault <
                ExpressionDeepVerificationFault::None ||
            request.verificationFault >
                ExpressionDeepVerificationFault::
                    FallbackIterationAllocation)
            return fail(
                ExpressionDeepRenderStatus::InvalidRequest,
                "verification fault selection is invalid");
        if (request.precision.requestedBits < 0 ||
            request.precision.viewBits < 0 ||
            request.precision.guardBits < 0 ||
            request.precision.minimumBits < 0 ||
            request.precision.maximumBits < MPFR_PREC_MIN ||
            request.precision.maximumBits > MPFR_PREC_MAX)
            return fail(
                ExpressionDeepRenderStatus::InvalidRequest,
                "reference precision policy is invalid");
        if (!finiteComplex(request.fixed.z) ||
            !finiteComplex(request.fixed.c) ||
            !finiteComplex(request.fixed.z0))
            return fail(
                ExpressionDeepRenderStatus::InvalidRequest,
                "fixed context is nonfinite");
        for (Complex parameter : request.fixed.parameters)
            if (!finiteComplex(parameter))
                return fail(
                    ExpressionDeepRenderStatus::InvalidRequest,
                    "fixed parameter is nonfinite");

        ExpressionProgram independentlySpecialized;
        if (request.canonicalProgram) {
            ExpressionError expressionError;
            if (!request.canonicalProgram->specialize(
                    request.fixed, request.pixelParameter,
                    independentlySpecialized,
                    &expressionError))
                return fail(
                    ExpressionDeepRenderStatus::InvalidRequest,
                    expressionError.message.empty()
                        ? "canonical specialization failed"
                        : expressionError.message);
            if (request.runtimeProgram->source() !=
                    request.canonicalProgram->source() ||
                !request.runtimeProgram->
                    semanticallyEquivalent(
                        independentlySpecialized))
                return fail(
                    ExpressionDeepRenderStatus::ProgramMismatch,
                    "runtime expression does not match canonical specialization");
        }

        const uint64_t width =
            static_cast<uint64_t>(request.width);
        const uint64_t height =
            static_cast<uint64_t>(request.height);
        if (height != 0 &&
            width >
                std::numeric_limits<uint64_t>::max() /
                    height)
            return fail(
                ExpressionDeepRenderStatus::ResourceLimit,
                "pixel count overflow");
        const uint64_t pixelCount64 = width * height;
        if (pixelCount64 >
                static_cast<uint64_t>(
                    std::numeric_limits<size_t>::max()) ||
            request.outputCount <
                static_cast<size_t>(pixelCount64))
            return fail(
                ExpressionDeepRenderStatus::InvalidRequest,
                "output buffer is smaller than the frame");
        if (pixelCount64 >
                std::numeric_limits<uint64_t>::max() /
                    static_cast<uint64_t>(
                        request.maxIterations))
            return fail(
                ExpressionDeepRenderStatus::ResourceLimit,
                "iteration resource multiplication overflow");
        const size_t pixelCount =
            static_cast<size_t>(pixelCount64);
        const uint64_t tilesX =
            (width +
             static_cast<uint64_t>(
                 request.threading.tileWidth) - 1) /
            static_cast<uint64_t>(
                request.threading.tileWidth);
        const uint64_t tilesY =
            (height +
             static_cast<uint64_t>(
                 request.threading.tileHeight) - 1) /
            static_cast<uint64_t>(
                request.threading.tileHeight);
        if (tilesY != 0 &&
            tilesX >
                std::numeric_limits<uint64_t>::max() /
                    tilesY)
            return fail(
                ExpressionDeepRenderStatus::ResourceLimit,
                "tile count overflow");
        const uint64_t tileCount = tilesX * tilesY;
        if (tileCount >
                static_cast<uint64_t>(LLONG_MAX) ||
            pixelCount64 >
                static_cast<uint64_t>(LLONG_MAX))
            return fail(
                ExpressionDeepRenderStatus::ResourceLimit,
                "OpenMP loop range is too large");

        int threadCount = request.threading.threads;
#ifdef _OPENMP
        if (threadCount == 0)
            threadCount = omp_get_max_threads();
#else
        threadCount = 1;
#endif
        threadCount = std::max(1, threadCount);

        result.capability =
            request.runtimeProgram->scaledResidualCapability();
        const bool certifiedTaylorCandidate =
            certifiedTaylorCapability(result.capability);
        const bool certifiedPiecewiseCandidate =
            result.capability ==
                ExpressionScaledResidualCapability::
                    CertifiedPiecewiseCandidate;
        const bool piecewisePerStepEligible =
            certifiedPiecewiseCandidate &&
            !request.runtimeProgram->
                scaledResidualRequiresTaylor();
        const bool certifiedTaylorEligible =
            !certifiedTaylorCandidate ||
            piecewisePerStepEligible ||
            (request.taylor.enableTaylor &&
             request.maxIterations >
                 request.taylor.minimumLanding);
        bool runFast =
            !request.forceMpfrFallbackForVerification &&
            certifiedTaylorEligible &&
            (certifiedReferenceCapability(
                 result.capability) ||
            (result.capability ==
                 ExpressionScaledResidualCapability::
                    UncertifiedSeries &&
             request.allowUncertifiedForBenchmark));
        bool certificationUnavailable =
            request.forceMpfrFallbackForVerification ||
            (certifiedTaylorCandidate &&
             !certifiedTaylorEligible);

        size_t rendererBaseBytes = 0;
        if ((runFast &&
             (!checkedAddSize(
                  rendererBaseBytes, request.width,
                  sizeof(ScaledOffset)) ||
              !checkedAddSize(
                  rendererBaseBytes, request.height,
                  sizeof(ScaledOffset)))) ||
            !checkedAddSize(
                rendererBaseBytes, pixelCount,
                sizeof(uint8_t)) ||
            !checkedAddSize(
                rendererBaseBytes, pixelCount,
                sizeof(size_t)))
            return fail(
                ExpressionDeepRenderStatus::ResourceLimit,
                "renderer memory calculation overflow");
        size_t fastThreadBytes = 2048;
        if (runFast && !checkedAddSize(
                fastThreadBytes,
                request.runtimeProgram->instructionCount(),
                sizeof(ScaledComplexValue) +
                    sizeof(ScaledRealValue) +
                    sizeof(uint16_t)))
            return fail(
                ExpressionDeepRenderStatus::ResourceLimit,
                "thread workspace calculation overflow");
        size_t fastThreadBytesTotal = 0;
        if (runFast && !checkedAddSize(
                fastThreadBytesTotal,
                static_cast<size_t>(threadCount),
                fastThreadBytes))
            return fail(
                ExpressionDeepRenderStatus::ResourceLimit,
                "thread workspace multiplication overflow");
        size_t rendererBytes = rendererBaseBytes;
        if (!checkedAddSize(
                rendererBytes, 1,
                fastThreadBytesTotal))
            return fail(
                ExpressionDeepRenderStatus::ResourceLimit,
                "renderer memory calculation overflow");
        result.rendererBytes = rendererBytes;
        if (request.memory.memoryLimitBytes != 0 &&
            rendererBytes >
                request.memory.memoryLimitBytes) {
            if (!certifiedTaylorCandidate)
                return fail(
                    ExpressionDeepRenderStatus::ResourceLimit,
                    "renderer exceeds memory limit");
            runFast = false;
            certificationUnavailable = true;
            rendererBaseBytes = 0;
            if (!checkedAddSize(
                    rendererBaseBytes, pixelCount,
                    sizeof(uint8_t)) ||
                !checkedAddSize(
                    rendererBaseBytes, pixelCount,
                    sizeof(size_t)))
                return fail(
                    ExpressionDeepRenderStatus::ResourceLimit,
                    "fallback renderer memory calculation overflow");
            fastThreadBytesTotal = 0;
            rendererBytes = rendererBaseBytes;
            result.rendererBytes = rendererBytes;
            if (rendererBytes >
                    request.memory.memoryLimitBytes)
                return fail(
                    ExpressionDeepRenderStatus::ResourceLimit,
                    "fallback renderer exceeds memory limit");
        }

        ExpressionReferencePrecisionPolicy precision =
            request.precision;
        std::string precisionError;
        if (!selectAutomaticViewBits(
                request, precision, precisionError))
            return fail(
                ExpressionDeepRenderStatus::InvalidRequest,
                precisionError);
        if (!selectDeepPrecision(
                precision, result.selectedPrecision,
                precisionError))
            return fail(
                ExpressionDeepRenderStatus::
                    PrecisionOutOfRange,
                precisionError);
        if (result.selectedPrecision <=
                MPFR_PREC_MAX -
                    request.memory.fallbackGuardBits)
            result.fallbackPrecision =
                result.selectedPrecision +
                request.memory.fallbackGuardBits;

        if (runFast &&
            certifiedReferenceCapability(
                result.capability)) {
            const mpfr_prec_t certificationGuard =
                std::max<mpfr_prec_t>(
                    128, precision.guardBits);
            if (certificationGuard >
                    precision.maximumBits ||
                certificationGuard >
                    ExpressionReferencePrecisionPolicy::
                        ApplicationMaximumBits ||
                result.selectedPrecision >
                    MPFR_PREC_MAX - certificationGuard ||
                result.selectedPrecision >
                    precision.maximumBits -
                        certificationGuard ||
                result.selectedPrecision >
                    ExpressionReferencePrecisionPolicy::
                        ApplicationMaximumBits -
                        certificationGuard) {
                runFast = false;
                certificationUnavailable = true;
            } else {
                result.certificationPrecision =
                    result.selectedPrecision +
                    certificationGuard;
            }
        }

        std::fill(
            request.output,
            request.output + pixelCount,
            ExpressionDeepEmptyPixel);

        std::atomic_bool cancelled{ false };
        auto pollCancellation = [&]() {
            if (cancelled.load(std::memory_order_acquire))
                return true;
            if (request.shouldCancel &&
                request.shouldCancel()) {
                cancelled.store(true);
                return true;
            }
            return false;
        };
        std::mutex progressMutex;
        auto notifyProgress = [&](
                ExpressionDeepRenderPhase phase,
                uint64_t completed, uint64_t total) {
            if (!request.progress) return;
            std::lock_guard<std::mutex> lock(progressMutex);
            request.progress(phase, completed, total);
        };

        ExpressionReferenceOrbitResult reference;
        if (runFast) {
            notifyProgress(
                ExpressionDeepRenderPhase::Reference,
                0, pixelCount64);
            ExpressionReferenceBuildRequest referenceRequest;
            referenceRequest.canonicalProgram =
                request.canonicalProgram;
            referenceRequest.runtimeProgram =
                request.runtimeProgram;
            referenceRequest.pixelParameter =
                request.pixelParameter;
            referenceRequest.center = request.center;
            referenceRequest.fixed = request.fixed;
            referenceRequest.bailout = request.bailout;
            referenceRequest.maxIterations =
                request.maxIterations;
            referenceRequest.precision = precision;
            referenceRequest.certificationPrecision =
                certifiedReferenceCapability(
                    result.capability)
                ? result.certificationPrecision : 0;
            referenceRequest.memoryLimitBytes =
                request.memory.memoryLimitBytes;
            referenceRequest.shouldCancel = pollCancellation;

            const Clock::time_point referenceStart =
                Clock::now();
            const bool referenceBuilt =
                buildExpressionReferenceOrbit(
                    referenceRequest, reference);
            result.referenceSeconds =
                secondsSince(referenceStart);
            result.referenceBytes = reference.memoryBytes;
            if (cancelled.load(std::memory_order_acquire) ||
                reference.cancelled)
                return fail(
                    ExpressionDeepRenderStatus::Cancelled,
                    "render cancelled while building the reference");
            if (!referenceBuilt) {
                const bool mayFallback =
                    certifiedReferenceCapability(
                        result.capability) &&
                    reference.status !=
                        ExpressionReferenceBuildStatus::
                            ProgramMismatch &&
                    reference.status !=
                        ExpressionReferenceBuildStatus::
                            InvalidRequest &&
                    reference.status !=
                        ExpressionReferenceBuildStatus::
                            InputParseError;
                if (!mayFallback)
                    return fail(
                        mapReferenceStatus(reference.status),
                        reference.error.empty()
                            ? "reference build failed"
                            : reference.error);
                runFast = false;
                certificationUnavailable = true;
                reference = {};
                result.referenceBytes = 0;
            } else if (!reference.valid) {
                runFast = false;
                certificationUnavailable = true;
                reference = {};
                result.referenceBytes = 0;
            } else if (
                certifiedReferenceCapability(
                    result.capability) &&
                (!reference.
                     certifiedAgainstHigherPrecision ||
                 reference.certificationPrecision !=
                     result.certificationPrecision)) {
                runFast = false;
                certificationUnavailable = true;
                reference = {};
                result.referenceBytes = 0;
            }
            if (runFast &&
                !inflateReferenceErrors(
                    reference,
                    request.verificationErrorInflationBits))
                return fail(
                    ExpressionDeepRenderStatus::InternalError,
                    "certification radius inflation overflow");
            if (runFast &&
                request.memory.memoryLimitBytes != 0 &&
                (reference.memoryBytes >
                     request.memory.memoryLimitBytes ||
                 rendererBytes >
                     request.memory.memoryLimitBytes -
                         reference.memoryBytes)) {
                if (!certifiedTaylorCandidate)
                    return fail(
                        ExpressionDeepRenderStatus::ResourceLimit,
                        "reference plus renderer exceeds memory limit");
                runFast = false;
                certificationUnavailable = true;
                reference = {};
                result.referenceBytes = 0;
                rendererBaseBytes = 0;
                if (!checkedAddSize(
                        rendererBaseBytes, pixelCount,
                        sizeof(uint8_t)) ||
                    !checkedAddSize(
                        rendererBaseBytes, pixelCount,
                        sizeof(size_t)))
                    return fail(
                        ExpressionDeepRenderStatus::ResourceLimit,
                        "fallback renderer memory calculation overflow");
                fastThreadBytesTotal = 0;
                rendererBytes = rendererBaseBytes;
                result.rendererBytes = rendererBytes;
            }
        }

        std::vector<uint8_t> fallbackReason(
            pixelCount, NoFallbackReason);
        std::vector<ScaledOffset> xOffsets;
        std::vector<ScaledOffset> yOffsets;
        ScaledComplexBall initialReference;
        bool initialReferenceExponentUnsafe = false;
        BailoutThreshold bailoutThreshold;
        std::unique_ptr<ExpressionScaledResidualEvaluator>
            validationEvaluator;
        ExpressionTaylorJetResult taylorJet;
        bool useTaylor = false;
        size_t retainedTaylorBytes = 0;
        if (runFast) {
            ExactGeometry geometry(reference.precision);
            ExactGeometry certificationGeometry(
                reference.certificationPrecision != 0
                ? reference.certificationPrecision
                : reference.precision);
            if (!geometry.initialize(request) ||
                !certificationGeometry.initialize(request))
                return fail(
                    ExpressionDeepRenderStatus::InternalError,
                    "failed to reconstruct reference geometry");
            if (!makeBailoutThreshold(
                    request.bailout, bailoutThreshold))
                return fail(
                    ExpressionDeepRenderStatus::InternalError,
                    "scaled geometry or bailout is not representable");

            xOffsets.resize(
                static_cast<size_t>(request.width));
            yOffsets.resize(
                static_cast<size_t>(request.height));
            ScaledComplexValue pixelBase;
            const ScaledComplexShadow& pixelPrimary =
                request.pixelParameter == FormulaParameter::C
                ? reference.c : reference.z0;
            const ScaledComplexShadow& pixelDefect =
                request.pixelParameter == FormulaParameter::C
                ? reference.cDefect : reference.z0Defect;
            if (makeScaledComplexValue(
                    pixelPrimary, pixelDefect,
                    pixelBase) !=
                        ScaledArithmeticStatus::Success)
                return fail(
                    ExpressionDeepRenderStatus::InternalError,
                    "compact pixel reference is not representable");
            ScaledComplexBall initialPrimary;
            ScaledComplexBall initialDefect;
            if (makeScaledComplexValue(
                    reference.initialZ,
                    initialPrimary.value) !=
                        ScaledArithmeticStatus::Success ||
                makeScaledComplexValue(
                    reference.initialZDefect,
                    initialDefect.value) !=
                        ScaledArithmeticStatus::Success ||
                certifiedScaledAdd(
                    initialPrimary, initialDefect,
                    initialReference) !=
                        ScaledArithmeticStatus::Success ||
                scaledAddUp(
                    initialReference.radius,
                    reference.initialZError,
                    initialReference.radius) !=
                        ScaledArithmeticStatus::Success)
                return fail(
                    ExpressionDeepRenderStatus::InternalError,
                    "compact initial reference is not representable");
            initialReferenceExponentUnsafe =
                certifyScaledMpfrExponentRange(
                    initialPrimary) !=
                        ScaledArithmeticStatus::Success ||
                certifyScaledMpfrExponentRange(
                    initialDefect) !=
                        ScaledArithmeticStatus::Success ||
                certifyScaledMpfrExponentRange(
                    initialReference) !=
                        ScaledArithmeticStatus::Success;
            MpfrComplex lowOffset(reference.precision);
            MpfrComplex exactCoordinate(
                certificationGeometry.center.precision());
            for (int x = 0; x < request.width; ++x) {
                const long centered =
                    static_cast<long>(
                        2LL * x - (request.width - 1LL));
                ScaledOffset& offset =
                    xOffsets[static_cast<size_t>(x)];
                mpfr_mul_si(
                    lowOffset.re, geometry.dxHalf,
                    centered, RND);
                if (makeScaledRealValue(
                        lowOffset.re, offset.value) !=
                        ScaledArithmeticStatus::Success ||
                    !certificationGeometry.coordinate(
                        x, 0, request, exactCoordinate) ||
                    !componentDiscrepancy(
                        exactCoordinate.re, pixelBase.re,
                        offset.value,
                        certificationGeometry.center.precision(),
                        offset.error) ||
                    !inflateRadius(
                        offset.error,
                        request.
                            verificationErrorInflationBits))
                    return fail(
                        ExpressionDeepRenderStatus::InternalError,
                        "certified x coordinate construction failed");
            }
            for (int y = 0; y < request.height; ++y) {
                const long centered =
                    static_cast<long>(
                        2LL * y - (request.height - 1LL));
                ScaledOffset& offset =
                    yOffsets[static_cast<size_t>(y)];
                mpfr_mul_si(
                    lowOffset.im, geometry.dyHalf,
                    centered, RND);
                if (makeScaledRealValue(
                        lowOffset.im, offset.value) !=
                        ScaledArithmeticStatus::Success ||
                    !certificationGeometry.coordinate(
                        0, y, request, exactCoordinate) ||
                    !componentDiscrepancy(
                        exactCoordinate.im, pixelBase.im,
                        offset.value,
                        certificationGeometry.center.precision(),
                        offset.error) ||
                    !inflateRadius(
                        offset.error,
                        request.
                            verificationErrorInflationBits))
                    return fail(
                        ExpressionDeepRenderStatus::InternalError,
                        "certified y coordinate construction failed");
            }
            validationEvaluator =
                std::make_unique<
                    ExpressionScaledResidualEvaluator>();
            if (reference.defectPending ||
                !validationEvaluator->prepare(
                    *request.runtimeProgram, reference)) {
                runFast = false;
                std::fill(
                    fallbackReason.begin(),
                    fallbackReason.end(),
                    static_cast<uint8_t>(
                        ExpressionDeepFallbackReason::
                            InvalidTape));
            } else {
                size_t validationThreadBytesTotal = 0;
                if (!checkedAddSize(
                        validationThreadBytesTotal,
                        static_cast<size_t>(threadCount),
                        validationEvaluator->workspaceBytes()))
                    return fail(
                        ExpressionDeepRenderStatus::ResourceLimit,
                        "validation workspace multiplication overflow");
                size_t validationPeak = rendererBaseBytes;
                if (!checkedAddSize(
                        validationPeak, 1,
                        validationThreadBytesTotal))
                    return fail(
                        ExpressionDeepRenderStatus::ResourceLimit,
                        "validation workspace calculation overflow");
                rendererBytes = std::max(
                    rendererBytes, validationPeak);
                result.rendererBytes = rendererBytes;
                if (request.memory.memoryLimitBytes != 0 &&
                    (reference.memoryBytes >
                         request.memory.memoryLimitBytes ||
                     rendererBytes >
                         request.memory.memoryLimitBytes -
                             reference.memoryBytes)) {
                    if (!certifiedTaylorCandidate)
                        return fail(
                            ExpressionDeepRenderStatus::ResourceLimit,
                            "validation workspace exceeds memory limit");
                    runFast = false;
                    certificationUnavailable = true;
                    reference = {};
                    result.referenceBytes = 0;
                    validationEvaluator.reset();
                    std::vector<ScaledOffset>().swap(xOffsets);
                    std::vector<ScaledOffset>().swap(yOffsets);
                    rendererBaseBytes = 0;
                    if (!checkedAddSize(
                            rendererBaseBytes, pixelCount,
                            sizeof(uint8_t)) ||
                        !checkedAddSize(
                            rendererBaseBytes, pixelCount,
                            sizeof(size_t)))
                        return fail(
                            ExpressionDeepRenderStatus::ResourceLimit,
                            "fallback renderer memory calculation overflow");
                    fastThreadBytesTotal = 0;
                    rendererBytes = rendererBaseBytes;
                    result.rendererBytes = rendererBytes;
                }
            }
            // Validation is complete; destroy its retained vector capacities
            // before Taylor construction and the parallel render peak.
            validationEvaluator.reset();

            if (runFast && request.taylor.enableTaylor &&
                certifiedReferenceCapability(
                    result.capability) &&
                request.maxIterations >
                    request.taylor.minimumLanding) {
                result.taylorAttempted = true;
                ScaledRealValue maximumRealMagnitude;
                ScaledRealValue maximumImaginaryMagnitude;
                ScaledRealValue maximumComponentError;
                for (const ScaledOffset& offset : xOffsets) {
                    const ScaledRealValue magnitude =
                        absoluteScaled(offset.value);
                    if (compareScaledNonnegative(
                            magnitude,
                            maximumRealMagnitude) > 0)
                        maximumRealMagnitude = magnitude;
                    if (compareScaledNonnegative(
                            offset.error,
                            maximumComponentError) > 0)
                        maximumComponentError =
                            offset.error;
                }
                for (const ScaledOffset& offset : yOffsets) {
                    const ScaledRealValue magnitude =
                        absoluteScaled(offset.value);
                    if (compareScaledNonnegative(
                            magnitude,
                            maximumImaginaryMagnitude) > 0)
                        maximumImaginaryMagnitude = magnitude;
                    if (compareScaledNonnegative(
                            offset.error,
                            maximumComponentError) > 0)
                        maximumComponentError =
                            offset.error;
                }
                ScaledComplexValue parameterScale;
                if (makeExpressionTaylorFrameScale(
                        maximumRealMagnitude,
                        maximumImaginaryMagnitude,
                        maximumComponentError,
                        parameterScale)) {
                    ExpressionTaylorJetRequest jetRequest;
                    jetRequest.program =
                        request.runtimeProgram;
                    jetRequest.reference = &reference;
                    jetRequest.pixelParameter =
                        request.pixelParameter;
                    jetRequest.parameterScale =
                        parameterScale;
                    jetRequest.minimumOrder =
                        request.taylor.minimumOrder;
                    jetRequest.preferredOrder =
                        request.taylor.order;
                    jetRequest.maximumOrder =
                        request.taylor.maximumOrder;
                    jetRequest.maximumBivariateOrder =
                        request.taylor.
                            maximumBivariateOrder;
                    jetRequest.maximumCompositionOrder =
                        request.taylor.
                            maximumCompositionOrder;
                    jetRequest.minimumLanding =
                        request.taylor.minimumLanding;
                    jetRequest.maximumCandidateIteration =
                        request.taylor.
                            maximumCandidateIteration;
                    jetRequest.bailout = request.bailout;
                    jetRequest.accuracyBudget =
                        request.taylor.accuracyBudget;
                    if (request.memory.memoryLimitBytes != 0 &&
                        reference.memoryBytes <=
                            request.memory.memoryLimitBytes &&
                        rendererBaseBytes <=
                            request.memory.memoryLimitBytes -
                                reference.memoryBytes)
                        jetRequest.memoryLimitBytes =
                            request.memory.memoryLimitBytes -
                            reference.memoryBytes -
                            rendererBaseBytes;
                    else
                        jetRequest.memoryLimitBytes =
                            request.memory.memoryLimitBytes;
                    jetRequest.shouldCancel =
                        pollCancellation;
                    useTaylor =
                        ExpressionTaylorJetBuilder::build(
                            jetRequest, taylorJet);
                } else {
                    taylorJet.status =
                        ExpressionTaylorJetStatus::
                            ExponentRange;
                    taylorJet.failureReason =
                        "Taylor frame normalization failed";
                }
                result.taylorBuildSeconds =
                    taylorJet.buildSeconds;
                result.taylorMemoryBytes =
                    taylorJet.memoryBytes;
                result.taylorStatus =
                    taylorJet.status;
                result.taylorFailureReason =
                    taylorJet.failureReason;
                result.taylorOrder = taylorJet.order;
                result.taylorLayout =
                    taylorJet.layout;
                result.taylorMonomialCount =
                    taylorJet.monomialCount;
                result.
                    taylorBivariateConvolutionOperationCount =
                        taylorJet.
                            bivariateConvolutionOperationCount;
                result.taylorCoveredIterations =
                    taylorJet.landingIteration;
                result.taylorMaximumFunctionSeriesOrder =
                    taylorJet.maximumFunctionSeriesOrder;
                result.taylorFunctionSeriesCount =
                    taylorJet.functionSeriesCount;
                result.taylorFunctionSeriesOperationCount =
                    taylorJet.functionSeriesOperationCount;
                result.taylorMaximumFunctionSeriesTail =
                    taylorJet.maximumFunctionSeriesTail;
                result.taylorMaximumReciprocalOrder =
                    taylorJet.maximumReciprocalOrder;
                result.taylorReciprocalCount =
                    taylorJet.reciprocalCount;
                result.taylorReciprocalOperationCount =
                    taylorJet.reciprocalOperationCount;
                result.taylorMinimumDenominatorClearance =
                    taylorJet.minimumDenominatorClearance;
                result.taylorMaximumReciprocalTail =
                    taylorJet.maximumReciprocalTail;
                result.taylorPoleRejected =
                    taylorJet.poleRejected;
                result.taylorMaximumBranchSeriesOrder =
                    taylorJet.maximumBranchSeriesOrder;
                result.taylorBranchCompositionCount =
                    taylorJet.branchCompositionCount;
                result.taylorBranchCompositionOperationCount =
                    taylorJet.branchCompositionOperationCount;
                result.taylorMaximumBranchSeriesTail =
                    taylorJet.maximumBranchSeriesTail;
                result.taylorMinimumBranchCutClearance =
                    taylorJet.minimumBranchCutClearance;
                result.taylorMinimumBranchZeroClearance =
                    taylorJet.minimumBranchZeroClearance;
                result.taylorBranchRejected =
                    taylorJet.branchRejected;
                result.taylorArgCompositionCount =
                    taylorJet.argCompositionCount;
                result.taylorArgRejectionReason =
                    taylorJet.argRejectionReason;
                result.taylorPolarCompositionCount =
                    taylorJet.polarCompositionCount;
                result.taylorMinimumPolarRadiusClearance =
                    taylorJet.minimumPolarRadiusClearance;
                result.taylorPolarRejected =
                    taylorJet.polarRejected;
                result.taylorPolarRejectionReason =
                    taylorJet.polarRejectionReason;
                result.taylorAbsBranchCount =
                    taylorJet.absBranchCount;
                result.taylorAbsPositiveCellCount =
                    taylorJet.absPositiveCellCount;
                result.taylorAbsNegativeCellCount =
                    taylorJet.absNegativeCellCount;
                result.taylorMinimumFoldClearance =
                    taylorJet.minimumFoldClearance;
                result.taylorFoldRejected =
                    taylorJet.foldRejected;
                result.taylorFoldRejectionIteration =
                    taylorJet.foldRejectionIteration;
                result.taylorFoldRejectionReason =
                    taylorJet.foldRejectionReason;
                if (cancelled.load(
                        std::memory_order_acquire) ||
                    taylorJet.status ==
                        ExpressionTaylorJetStatus::
                            Cancelled)
                    return fail(
                        ExpressionDeepRenderStatus::
                            Cancelled,
                        "render cancelled while building Taylor jet");
                if (useTaylor &&
                    certifiedTaylorCapability(
                        result.capability) &&
                    (!certifiedPiecewiseCandidate ||
                     taylorJet.argCompositionCount > 0 ||
                     taylorJet.polarCompositionCount > 0) &&
                    !request.allowUncertifiedForBenchmark &&
                    taylorJet.landingIteration <
                        request.maxIterations) {
                    useTaylor = false;
                    result.taylorFailureReason =
                        "certified Taylor jet does not cover the full iteration horizon";
                }
                if (useTaylor &&
                    request.taylor.
                        requirePredictedBenefit) {
                    const long double saved =
                        static_cast<long double>(
                            pixelCount64) *
                        taylorJet.landingIteration * 16.0L;
                    const long double evaluation =
                        static_cast<long double>(
                            pixelCount64) *
                        (taylorJet.layout ==
                                 ExpressionTaylorJetLayout::
                                     RealBivariate
                             ? 2.0L *
                                   static_cast<long double>(
                                       taylorJet.
                                           monomialCount) +
                                   taylorJet.order + 2.0L
                             : 2.0L * taylorJet.order +
                                   2.0L);
                    const long double predictedCost =
                        static_cast<long double>(
                            taylorJet.operationCount) +
                        evaluation;
                    if (saved <= predictedCost) {
                        useTaylor = false;
                        result.taylorFailureReason =
                            "Taylor cost predicts no frame benefit";
                    }
                }
                if (useTaylor) {
                    size_t persistentTaylorBytes = 0;
                    if (!checkedAddSize(
                            persistentTaylorBytes,
                            taylorJet.coefficients.capacity(),
                            sizeof(ScaledComplexValue)) ||
                        !checkedAddSize(
                            persistentTaylorBytes,
                            taylorJet.coefficientRadii.capacity(),
                            sizeof(ScaledRealValue)) ||
                        !checkedAddSize(
                            persistentTaylorBytes,
                            taylorJet.
                                intermediateEscapeMargins.
                                    capacity(),
                            sizeof(ScaledRealValue))) {
                        useTaylor = false;
                        result.taylorFailureReason =
                            "Taylor coefficient memory calculation overflow";
                    }
                    size_t fastWithTaylor = rendererBytes;
                    size_t buildWithTaylor =
                        rendererBaseBytes;
                    if (useTaylor &&
                        (!checkedAddSize(
                             fastWithTaylor, 1,
                             persistentTaylorBytes) ||
                         !checkedAddSize(
                             buildWithTaylor, 1,
                             taylorJet.memoryBytes))) {
                        useTaylor = false;
                        result.taylorFailureReason =
                            "Taylor peak memory calculation overflow";
                    }
                    const size_t withTaylor = std::max(
                        fastWithTaylor, buildWithTaylor);
                    if (useTaylor &&
                        (request.memory.memoryLimitBytes != 0 &&
                         (reference.memoryBytes >
                              request.memory.memoryLimitBytes ||
                          withTaylor >
                              request.memory.memoryLimitBytes -
                                  reference.memoryBytes))) {
                        useTaylor = false;
                        result.taylorFailureReason =
                            "Taylor coefficients exceed renderer memory policy";
                    } else if (useTaylor) {
                        retainedTaylorBytes =
                            persistentTaylorBytes;
                        rendererBytes = withTaylor;
                        result.rendererBytes = rendererBytes;
                    }
                }
                if (!useTaylor) {
                    std::vector<ScaledComplexValue>().swap(
                        taylorJet.coefficients);
                    std::vector<ScaledRealValue>().swap(
                        taylorJet.coefficientRadii);
                    std::vector<ScaledRealValue>().swap(
                        taylorJet.intermediateEscapeMargins);
                }
                result.taylorAccepted = useTaylor;
            }
            if (runFast &&
                certifiedTaylorCapability(
                    result.capability) &&
                !piecewisePerStepEligible &&
                !request.allowUncertifiedForBenchmark &&
                !useTaylor) {
                runFast = false;
                certificationUnavailable = true;
                reference = {};
                result.referenceBytes = 0;
                std::vector<ScaledOffset>().swap(xOffsets);
                std::vector<ScaledOffset>().swap(yOffsets);
                retainedTaylorBytes = 0;
                rendererBaseBytes = 0;
                if (!checkedAddSize(
                        rendererBaseBytes, pixelCount,
                        sizeof(uint8_t)) ||
                    !checkedAddSize(
                        rendererBaseBytes, pixelCount,
                        sizeof(size_t)))
                    return fail(
                        ExpressionDeepRenderStatus::ResourceLimit,
                        "fallback renderer memory calculation overflow");
                fastThreadBytesTotal = 0;
                rendererBytes = rendererBaseBytes;
                result.rendererBytes = rendererBytes;
                std::fill(
                    fallbackReason.begin(),
                    fallbackReason.end(),
                    static_cast<uint8_t>(
                        ExpressionDeepFallbackReason::
                            CertificationFailure));
            }
        }
        if (!runFast &&
            fallbackReason[0] == NoFallbackReason) {
            const uint8_t reason = static_cast<uint8_t>(
                certificationUnavailable
                ? ExpressionDeepFallbackReason::
                    CertificationFailure
                : reasonForCapability(result.capability));
            std::fill(
                fallbackReason.begin(),
                fallbackReason.end(), reason);
        }

        notifyProgress(
            ExpressionDeepRenderPhase::Fast,
            0, pixelCount64);
        const Clock::time_point fastStart = Clock::now();
        std::atomic<uint64_t> fastPixels{ 0 };
        std::atomic<uint64_t> totalIterations{ 0 };
        std::atomic<uint64_t> taylorAcceptedPixels{ 0 };
        std::atomic<uint64_t> taylorFallbackPixels{ 0 };
        std::atomic<uint64_t> taylorEvaluationNanoseconds{ 0 };
        std::atomic<uint64_t> taylorResidualNanoseconds{ 0 };
        std::atomic<uint64_t> fastCompleted{ 0 };
        std::atomic_bool fastResourceError{ false };
        std::atomic_bool fastInternalError{ false };
        std::atomic_bool fastIterationFaultInjected{ false };
        if (runFast) {
#pragma omp parallel num_threads(threadCount)
            {
                ExpressionScaledResidualEvaluator evaluator;
                bool prepared = false;
                bool workerReady = true;
                try {
                    if (request.verificationFault ==
                            ExpressionDeepVerificationFault::
                                FastWorkerAllocation)
                        throw std::bad_alloc();
                    prepared = evaluator.prepare(
                        *request.runtimeProgram, reference);
                } catch (const std::bad_alloc&) {
                    fastResourceError.store(
                        true, std::memory_order_release);
                    workerReady = false;
                } catch (const std::length_error&) {
                    fastResourceError.store(
                        true, std::memory_order_release);
                    workerReady = false;
                } catch (...) {
                    fastInternalError.store(
                        true, std::memory_order_release);
                    workerReady = false;
                }
                uint64_t localIterations = 0;
#pragma omp for schedule(dynamic, 1)
                for (long long tile = 0;
                     tile < static_cast<long long>(tileCount);
                     ++tile) {
                  try {
                    if (!workerReady ||
                        fastResourceError.load(
                            std::memory_order_acquire) ||
                        fastInternalError.load(
                            std::memory_order_acquire))
                        continue;
                    if (request.verificationFault ==
                            ExpressionDeepVerificationFault::
                                FastIterationAllocation &&
                        !fastIterationFaultInjected.exchange(
                            true, std::memory_order_acq_rel))
                        throw std::bad_alloc();
                    if (pollCancellation()) continue;
                    const int tileX =
                        static_cast<int>(
                            tile % static_cast<long long>(
                                tilesX));
                    const int tileY =
                        static_cast<int>(
                            tile / static_cast<long long>(
                                tilesX));
                    const int xBegin =
                        tileX * request.threading.tileWidth;
                    const int yBegin =
                        tileY * request.threading.tileHeight;
                    const int xEnd = std::min(
                        request.width,
                        xBegin + request.threading.tileWidth);
                    const int yEnd = std::min(
                        request.height,
                        yBegin + request.threading.tileHeight);
                    uint64_t completedInTile = 0;
                    for (int y = yBegin; y < yEnd; ++y) {
                        for (int x = xBegin; x < xEnd; ++x) {
                            if (pollCancellation()) break;
                            const size_t index =
                                static_cast<size_t>(y) *
                                    request.width + x;
                            FastPixelResult pixel;
                            if (!prepared) {
                                pixel.reason =
                                    ExpressionDeepFallbackReason::
                                        InvalidTape;
                            } else {
                                ScaledComplexBall offset;
                                const ScaledOffset& xOffset =
                                    xOffsets[
                                        static_cast<size_t>(x)];
                                const ScaledOffset& yOffset =
                                    yOffsets[
                                        static_cast<size_t>(y)];
                                offset.value.re = xOffset.value;
                                offset.value.im = yOffset.value;
                                offset.radius =
                                    compareScaledNonnegative(
                                        xOffset.error,
                                        yOffset.error) >= 0
                                    ? xOffset.error
                                    : yOffset.error;
                                ExpressionScaledResidualInput input;
                                ScaledComplexBall stateDelta;
                                ScaledArithmeticStatus arithmetic =
                                    initialReferenceExponentUnsafe
                                    ? ScaledArithmeticStatus::
                                        ExponentRange
                                    : certifyScaledMpfrExponentRange(
                                        offset);
                                if (arithmetic !=
                                        ScaledArithmeticStatus::
                                            Success) {
                                    pixel.reason =
                                        reasonForArithmeticStatus(
                                            arithmetic);
                                } else if (
                                    request.pixelParameter ==
                                        FormulaParameter::C) {
                                    input.c = offset.value;
                                    input.cError =
                                        offset.radius;
                                    stateDelta.radius =
                                        initialReference.radius;
                                } else {
                                    input.z0 = offset.value;
                                    input.z0Error =
                                        offset.radius;
                                    stateDelta = offset;
                                }

                                ScaledComplexBall initialBase;
                                initialBase.value =
                                    initialReference.value;
                                ScaledComplexBall initialValue;
                                if (pixel.reason ==
                                        ExpressionDeepFallbackReason::
                                            InvalidTape)
                                    arithmetic =
                                        certifyScaledMpfrExponentRange(
                                            initialBase);
                                if (arithmetic ==
                                        ScaledArithmeticStatus::
                                            Success)
                                    arithmetic =
                                        certifyScaledMpfrExponentRange(
                                            stateDelta);
                                if (arithmetic ==
                                        ScaledArithmeticStatus::
                                            Success)
                                    arithmetic = certifiedScaledAdd(
                                            initialBase, stateDelta,
                                            initialValue);
                                if (arithmetic ==
                                        ScaledArithmeticStatus::
                                            Success)
                                    arithmetic =
                                        certifyScaledMpfrExponentRange(
                                            initialValue);
                                if (arithmetic !=
                                        ScaledArithmeticStatus::
                                            Success) {
                                    pixel.reason =
                                        reasonForArithmeticStatus(
                                            arithmetic);
                                } else {
                                    ScaledArithmeticStatus gateStatus =
                                        ScaledArithmeticStatus::
                                            Success;
                                    BailoutDecision decision =
                                        decideBailout(
                                            initialValue.value,
                                            initialValue.radius,
                                            bailoutThreshold,
                                            gateStatus);
                                    if (decision ==
                                            BailoutDecision::Error) {
                                        pixel.reason =
                                            reasonForArithmeticStatus(
                                                gateStatus);
                                    } else if (decision ==
                                            BailoutDecision::Uncertain) {
                                        pixel.reason =
                                            ExpressionDeepFallbackReason::
                                                BailoutUncertain;
                                    } else if (decision ==
                                            BailoutDecision::Outside) {
                                        pixel.decided = true;
                                        pixel.output = 0.0f;
                                    } else {
                                        int firstIteration = 0;
                                        if (useTaylor) {
                                            const Clock::time_point
                                                evaluationStart =
                                                    Clock::now();
                                            ScaledComplexBall q;
                                            ScaledComplexBall
                                                landingDelta;
                                            bool landed = false;
                                            if (makeExpressionTaylorNormalizedQ(
                                                    offset,
                                                    taylorJet.
                                                        parameterScale,
                                                    q) &&
                                                expressionTaylorQInsideUnitDisk(
                                                    q)) {
                                                ExpressionTaylorJetEvaluation
                                                    landing;
                                                if (ExpressionTaylorJetEvaluator::
                                                        evaluate(
                                                            taylorJet,
                                                            q,
                                                            landing)) {
                                                    landingDelta =
                                                        landing.residual;
                                                    const ExpressionReferenceSample&
                                                        landingSample =
                                                            reference.samples[
                                                                taylorJet.
                                                                    landingSample];
                                                    ScaledComplexBall
                                                        landingBase;
                                                    ScaledComplexBall
                                                        landingValue;
                                                    arithmetic =
                                                        makeScaledComplexValue(
                                                            taylorJet.
                                                                landingUsesSampleOutput
                                                            ? landingSample.next
                                                            : landingSample.z,
                                                            taylorJet.
                                                                landingUsesSampleOutput
                                                            ? landingSample.
                                                                  rootDefect
                                                            : landingSample.
                                                                  zDefect,
                                                            landingBase.value);
                                                    if (arithmetic ==
                                                            ScaledArithmeticStatus::
                                                                Success)
                                                        arithmetic =
                                                            certifiedScaledAdd(
                                                                landingBase,
                                                                landingDelta,
                                                                landingValue);
                                                    if (arithmetic ==
                                                            ScaledArithmeticStatus::
                                                                Success)
                                                        arithmetic =
                                                            certifyScaledMpfrExponentRange(
                                                                landingValue);
                                                    if (arithmetic ==
                                                            ScaledArithmeticStatus::
                                                                Success) {
                                                        decision =
                                                            decideBailout(
                                                                landingValue.
                                                                    value,
                                                                landingValue.
                                                                    radius,
                                                                bailoutThreshold,
                                                                gateStatus);
                                                        landed =
                                                            decision ==
                                                                BailoutDecision::
                                                                    Inside;
                                                    }
                                                }
                                            }
                                            taylorEvaluationNanoseconds.
                                                fetch_add(
                                                    static_cast<uint64_t>(
                                                        std::chrono::
                                                            duration_cast<
                                                                std::chrono::
                                                                    nanoseconds>(
                                                                Clock::now() -
                                                                evaluationStart).
                                                            count()),
                                                    std::memory_order_relaxed);
                                            if (landed) {
                                                stateDelta =
                                                    landingDelta;
                                                firstIteration =
                                                    taylorJet.
                                                        landingIteration;
                                                taylorAcceptedPixels.
                                                    fetch_add(
                                                        1,
                                                        std::memory_order_relaxed);
                                                if (firstIteration ==
                                                        request.
                                                            maxIterations) {
                                                    pixel.decided = true;
                                                    pixel.output =
                                                        ExpressionDeepInteriorPixel;
                                                }
                                            } else {
                                                taylorFallbackPixels.
                                                    fetch_add(
                                                        1,
                                                        std::memory_order_relaxed);
                                                if (certifiedTaylorCapability(
                                                        result.capability) &&
                                                    !certifiedPiecewiseCandidate &&
                                                    !request.
                                                        allowUncertifiedForBenchmark) {
                                                    pixel.reason =
                                                        ExpressionDeepFallbackReason::
                                                            CertificationFailure;
                                                    firstIteration =
                                                        request.maxIterations;
                                                }
                                            }
                                        }
                                        const Clock::time_point
                                            residualStart =
                                                Clock::now();
                                        for (int iteration =
                                                 firstIteration;
                                             iteration <
                                                 request.maxIterations;
                                             ++iteration) {
                                            if (pollCancellation())
                                                break;
                                            if (static_cast<size_t>(
                                                    iteration) >=
                                                    reference.samples.size()) {
                                                pixel.reason =
                                                    ExpressionDeepFallbackReason::
                                                        ReferenceExhausted;
                                                break;
                                            }
                                            input.z =
                                                stateDelta.value;
                                            input.zError =
                                                stateDelta.radius;
                                            input.iteration = iteration;
                                            ExpressionScaledResidualResult
                                                evaluated =
                                                    evaluator.evaluate(
                                                        static_cast<size_t>(
                                                            iteration),
                                                        input);
                                            ++pixel.iterations;
                                            if (evaluated.status !=
                                                    ExpressionScaledResidualStatus::
                                                        Success) {
                                                pixel.reason =
                                                    reasonForResidualStatus(
                                                        evaluated.status);
                                                break;
                                            }
                                            if (evaluated.uncertified &&
                                                !request.
                                                    allowUncertifiedForBenchmark) {
                                                pixel.reason =
                                                    ExpressionDeepFallbackReason::
                                                        UncertifiedSeries;
                                                break;
                                            }
                                            if ((result.capability ==
                                                     ExpressionScaledResidualCapability::
                                                         ExactCenteredArithmetic ||
                                                 certifiedPiecewiseCandidate) &&
                                                !evaluated.certified) {
                                                pixel.reason =
                                                    ExpressionDeepFallbackReason::
                                                        CertificationFailure;
                                                break;
                                            }
                                            const ExpressionReferenceSample&
                                                sample =
                                                    reference.samples[
                                                        static_cast<size_t>(
                                                            iteration)];
                                            ScaledComplexBall outputBase;
                                            ScaledComplexBall residualBall;
                                            ScaledComplexBall actualOutput;
                                            arithmetic =
                                                makeScaledComplexValue(
                                                    sample.next,
                                                    sample.rootDefect,
                                                    outputBase.value);
                                            if (arithmetic ==
                                                    ScaledArithmeticStatus::
                                                        Success)
                                                arithmetic =
                                                    certifyScaledMpfrExponentRange(
                                                        outputBase);
                                            residualBall.value =
                                                evaluated.residual;
                                            residualBall.radius =
                                                evaluated.radius;
                                            if (arithmetic ==
                                                    ScaledArithmeticStatus::
                                                        Success)
                                                arithmetic =
                                                    certifiedScaledAdd(
                                                    outputBase,
                                                    residualBall,
                                                    actualOutput);
                                            if (arithmetic ==
                                                    ScaledArithmeticStatus::
                                                        Success)
                                                arithmetic =
                                                    certifyScaledMpfrExponentRange(
                                                        actualOutput);
                                            if (arithmetic !=
                                                    ScaledArithmeticStatus::
                                                        Success) {
                                                pixel.reason =
                                                    reasonForArithmeticStatus(
                                                        arithmetic);
                                                break;
                                            }
                                            decision = decideBailout(
                                                actualOutput.value,
                                                actualOutput.radius,
                                                bailoutThreshold,
                                                gateStatus);
                                            if (decision ==
                                                    BailoutDecision::Error) {
                                                pixel.reason =
                                                    reasonForArithmeticStatus(
                                                        gateStatus);
                                                break;
                                            }
                                            if (decision ==
                                                    BailoutDecision::Uncertain) {
                                                pixel.reason =
                                                    ExpressionDeepFallbackReason::
                                                        BailoutUncertain;
                                                break;
                                            }
                                            if (decision ==
                                                    BailoutDecision::Outside) {
                                                pixel.decided = true;
                                                pixel.output =
                                                    static_cast<float>(
                                                        iteration + 1);
                                                break;
                                            }
                                            if (iteration + 1 ==
                                                    request.maxIterations) {
                                                pixel.decided = true;
                                                pixel.output =
                                                    ExpressionDeepInteriorPixel;
                                                break;
                                            }
                                            if (static_cast<size_t>(
                                                    iteration + 1) >=
                                                    reference.samples.size()) {
                                                pixel.reason =
                                                    ExpressionDeepFallbackReason::
                                                        ReferenceExhausted;
                                                break;
                                            }
                                            ScaledComplexBall nextBase;
                                            arithmetic =
                                                makeScaledComplexValue(
                                                    reference.samples[
                                                        static_cast<size_t>(
                                                            iteration + 1)].
                                                        z,
                                                    reference.samples[
                                                        static_cast<size_t>(
                                                            iteration + 1)].
                                                        zDefect,
                                                    nextBase.value);
                                            if (arithmetic ==
                                                    ScaledArithmeticStatus::
                                                        Success)
                                                arithmetic =
                                                    certifyScaledMpfrExponentRange(
                                                        nextBase);
                                            if (arithmetic ==
                                                    ScaledArithmeticStatus::
                                                        Success)
                                                arithmetic =
                                                    certifiedScaledSubtract(
                                                        actualOutput,
                                                        nextBase,
                                                        stateDelta);
                                            if (arithmetic ==
                                                    ScaledArithmeticStatus::
                                                        Success)
                                                arithmetic =
                                                    certifyScaledMpfrExponentRange(
                                                        stateDelta);
                                            if (arithmetic !=
                                                    ScaledArithmeticStatus::
                                                        Success) {
                                                pixel.reason =
                                                    reasonForArithmeticStatus(
                                                        arithmetic);
                                                break;
                                            }
                                        }
                                        if (useTaylor)
                                            taylorResidualNanoseconds.
                                                fetch_add(
                                                    static_cast<uint64_t>(
                                                        std::chrono::
                                                            duration_cast<
                                                                std::chrono::
                                                                    nanoseconds>(
                                                                Clock::now() -
                                                                residualStart).
                                                            count()),
                                                    std::memory_order_relaxed);
                                    }
                                }
                            }
                            localIterations += pixel.iterations;
                            if (cancelled.load(
                                    std::memory_order_acquire))
                                break;
                            if (pixel.decided) {
                                if (!cancelled.load()) {
                                    request.output[index] =
                                        pixel.output;
                                    if (cancelled.load()) {
                                        request.output[index] =
                                            ExpressionDeepEmptyPixel;
                                    } else {
                                        fastPixels.fetch_add(
                                            1,
                                            std::memory_order_relaxed);
                                    }
                                }
                            } else {
                                fallbackReason[index] =
                                    static_cast<uint8_t>(
                                        pixel.reason);
                            }
                            ++completedInTile;
                        }
                        if (cancelled.load(
                                std::memory_order_acquire))
                            break;
                    }
                    const uint64_t completed =
                        fastCompleted.fetch_add(
                            completedInTile,
                            std::memory_order_relaxed) +
                        completedInTile;
                    notifyProgress(
                        ExpressionDeepRenderPhase::Fast,
                        completed, pixelCount64);
                  } catch (const std::bad_alloc&) {
                    fastResourceError.store(
                        true, std::memory_order_release);
                  } catch (const std::length_error&) {
                    fastResourceError.store(
                        true, std::memory_order_release);
                  } catch (...) {
                    fastInternalError.store(
                        true, std::memory_order_release);
                  }
                }
                totalIterations.fetch_add(
                    localIterations,
                    std::memory_order_relaxed);
            }
        }
        result.fastSeconds = secondsSince(fastStart);
        result.fastPixelCount =
            fastPixels.load(std::memory_order_relaxed);
        result.totalIterations =
            totalIterations.load(std::memory_order_relaxed);
        result.taylorAcceptedPixelCount =
            taylorAcceptedPixels.load(
                std::memory_order_relaxed);
        result.taylorFallbackPixelCount =
            taylorFallbackPixels.load(
                std::memory_order_relaxed);
        result.taylorEvaluationSeconds =
            static_cast<double>(
                taylorEvaluationNanoseconds.load(
                    std::memory_order_relaxed)) /
            1.0e9;
        result.taylorResidualSeconds =
            static_cast<double>(
                taylorResidualNanoseconds.load(
                    std::memory_order_relaxed)) /
            1.0e9;
        if (fastResourceError.load(
                std::memory_order_acquire))
            return fail(
                ExpressionDeepRenderStatus::ResourceLimit,
                "fast worker allocation failed");
        if (fastInternalError.load(
                std::memory_order_acquire))
            return fail(
                ExpressionDeepRenderStatus::InternalError,
                "fast worker failed");

        std::vector<size_t> fallbackQueue;
        fallbackQueue.reserve(pixelCount);
        for (size_t index = 0; index < pixelCount; ++index) {
            const uint8_t encoded = fallbackReason[index];
            if (encoded == NoFallbackReason) continue;
            fallbackQueue.push_back(index);
            const auto reason =
                static_cast<ExpressionDeepFallbackReason>(
                    encoded);
            const size_t reasonIndex =
                static_cast<size_t>(reason);
            if (reasonIndex <
                    result.fallbackReasonCounts.size())
                ++result.fallbackReasonCounts[reasonIndex];
            if (reason ==
                    ExpressionDeepFallbackReason::
                        UncertifiedSeries ||
                reason ==
                    ExpressionDeepFallbackReason::
                        BranchSensitive ||
                reason ==
                    ExpressionDeepFallbackReason::
                        CertificationFailure ||
                reason ==
                    ExpressionDeepFallbackReason::
                        BailoutUncertain)
                ++result.uncertainPixelCount;
        }
        result.fallbackPixelCount =
            static_cast<uint64_t>(fallbackQueue.size());
        for (uint64_t count : result.fallbackReasonCounts)
            result.maxFallbackReasonCount =
                std::max(
                    result.maxFallbackReasonCount, count);

        for (uint64_t tile = 0; tile < tileCount; ++tile) {
            const int tileX =
                static_cast<int>(tile % tilesX);
            const int tileY =
                static_cast<int>(tile / tilesX);
            const int xBegin =
                tileX * request.threading.tileWidth;
            const int yBegin =
                tileY * request.threading.tileHeight;
            const int xEnd = std::min(
                request.width,
                xBegin + request.threading.tileWidth);
            const int yEnd = std::min(
                request.height,
                yBegin + request.threading.tileHeight);
            uint64_t fallbackInTile = 0;
            for (int y = yBegin; y < yEnd; ++y)
                for (int x = xBegin; x < xEnd; ++x)
                    fallbackInTile +=
                        fallbackReason[
                            static_cast<size_t>(y) *
                                request.width + x] !=
                        NoFallbackReason;
            if (fallbackInTile == 0) continue;
            ++result.fallbackTileCount;
            const uint64_t pixelsInTile =
                static_cast<uint64_t>(xEnd - xBegin) *
                static_cast<uint64_t>(yEnd - yBegin);
            result.maxTileFallbackRate = std::max(
                result.maxTileFallbackRate,
                static_cast<double>(fallbackInTile) /
                    static_cast<double>(pixelsInTile));
        }

        if (cancelled.load(std::memory_order_acquire))
            return fail(
                ExpressionDeepRenderStatus::Cancelled,
                "render cancelled during the fast pass");

        if (!fallbackQueue.empty()) {
            if (result.selectedPrecision >
                    MPFR_PREC_MAX -
                        request.memory.fallbackGuardBits)
                return fail(
                    ExpressionDeepRenderStatus::
                        PrecisionOutOfRange,
                    "fallback precision overflow");
            result.fallbackPrecision =
                result.selectedPrecision +
                request.memory.fallbackGuardBits;
            if (result.fallbackPrecision >
                    ExpressionReferencePrecisionPolicy::
                        ApplicationMaximumBits ||
                result.fallbackPrecision >
                    request.precision.maximumBits)
                return fail(
                    ExpressionDeepRenderStatus::
                        PrecisionOutOfRange,
                    "fallback precision exceeds policy maximum");
            const size_t limbs =
                (static_cast<size_t>(
                     result.fallbackPrecision) +
                 GMP_NUMB_BITS - 1) /
                GMP_NUMB_BITS;
            size_t mpfrValueBytes =
                sizeof(__mpfr_struct);
            if (!checkedAddSize(
                    mpfrValueBytes, limbs,
                    sizeof(mp_limb_t)))
                return fail(
                    ExpressionDeepRenderStatus::ResourceLimit,
                    "MPFR workspace calculation overflow");
            // Worker state (34), oracle scratch (12), and up to six
            // simultaneous operation temporaries, plus the bytecode stack.
            size_t mpfrValuesPerThread = 52;
            if (!checkedAddSize(
                    mpfrValuesPerThread,
                    request.runtimeProgram->stackDepth(), 2))
                return fail(
                    ExpressionDeepRenderStatus::ResourceLimit,
                    "oracle stack calculation overflow");
            size_t fallbackThreadBytes = 0;
            if (!checkedAddSize(
                    fallbackThreadBytes,
                    mpfrValuesPerThread,
                    mpfrValueBytes))
                return fail(
                    ExpressionDeepRenderStatus::ResourceLimit,
                    "fallback workspace calculation overflow");
            size_t fallbackThreadBytesTotal = 0;
            if (!checkedAddSize(
                    fallbackThreadBytesTotal,
                    static_cast<size_t>(threadCount),
                    fallbackThreadBytes))
                return fail(
                    ExpressionDeepRenderStatus::ResourceLimit,
                    "fallback workspace multiplication overflow");
            size_t fallbackRendererBytes =
                rendererBaseBytes;
            if (!checkedAddSize(
                    fallbackRendererBytes, 1,
                    retainedTaylorBytes) ||
                !checkedAddSize(
                    fallbackRendererBytes, 1,
                    std::max(
                        fastThreadBytesTotal,
                        fallbackThreadBytesTotal)))
                return fail(
                    ExpressionDeepRenderStatus::ResourceLimit,
                    "renderer peak memory calculation overflow");
            rendererBytes = std::max(
                rendererBytes, fallbackRendererBytes);
            result.rendererBytes = rendererBytes;
            if (request.memory.memoryLimitBytes != 0 &&
                (reference.memoryBytes >
                     request.memory.memoryLimitBytes ||
                 rendererBytes >
                     request.memory.memoryLimitBytes -
                         reference.memoryBytes))
                return fail(
                    ExpressionDeepRenderStatus::ResourceLimit,
                    "MPFR fallback exceeds memory limit");
        }

        notifyProgress(
            ExpressionDeepRenderPhase::Fallback,
            0,
            static_cast<uint64_t>(fallbackQueue.size()));
        const Clock::time_point fallbackStart = Clock::now();
        std::atomic<uint64_t> fallbackCompleted{ 0 };
        std::atomic<uint64_t> undefinedPixels{ 0 };
        std::atomic_bool fallbackResourceError{ false };
        std::atomic_bool fallbackInternalError{ false };
        std::atomic_bool fallbackIterationFaultInjected{
            false
        };
        if (!fallbackQueue.empty()) {
#pragma omp parallel num_threads(threadCount)
            {
                std::unique_ptr<FallbackWorkerWorkspace>
                    workspace;
                try {
                    if (request.verificationFault ==
                            ExpressionDeepVerificationFault::
                                FallbackWorkerAllocation)
                        throw std::bad_alloc();
                    workspace =
                        std::make_unique<
                            FallbackWorkerWorkspace>(
                                result.fallbackPrecision,
                                request);
                    if (!workspace->geometryReady)
                        fallbackInternalError.store(
                            true,
                            std::memory_order_release);
                } catch (const std::bad_alloc&) {
                    fallbackResourceError.store(
                        true, std::memory_order_release);
                } catch (const std::length_error&) {
                    fallbackResourceError.store(
                        true, std::memory_order_release);
                } catch (...) {
                    fallbackInternalError.store(
                        true, std::memory_order_release);
                }
                uint64_t localIterations = 0;
#pragma omp for schedule(dynamic, 8)
                for (long long queueIndex = 0;
                     queueIndex <
                         static_cast<long long>(
                             fallbackQueue.size());
                     ++queueIndex) {
                  try {
                    if (!workspace ||
                        fallbackResourceError.load(
                            std::memory_order_acquire) ||
                        fallbackInternalError.load(
                            std::memory_order_acquire))
                        continue;
                    if (request.verificationFault ==
                            ExpressionDeepVerificationFault::
                                FallbackIterationAllocation &&
                        !fallbackIterationFaultInjected.exchange(
                            true, std::memory_order_acq_rel))
                        throw std::bad_alloc();
                    if (pollCancellation()) continue;
                    ExactGeometry& geometry =
                        workspace->geometry;
                    ExpressionOracleContext& context =
                        workspace->context;
                    MpfrComplex& pixel = workspace->pixel;
                    MpfrComplex& next = workspace->next;
                    mpfr_ptr magnitude =
                        workspace->magnitudeStorage.re;
                    const size_t outputIndex =
                        fallbackQueue[
                            static_cast<size_t>(queueIndex)];
                    const int y = static_cast<int>(
                        outputIndex /
                            static_cast<size_t>(
                                request.width));
                    const int x = static_cast<int>(
                        outputIndex %
                            static_cast<size_t>(
                                request.width));
                    configureFixed(request.fixed, context);
                    bool undefined =
                        !geometry.coordinate(
                            x, y, request, pixel);
                    if (!undefined) {
                        if (request.pixelParameter ==
                                FormulaParameter::C) {
                            context.c.set(pixel);
                        } else {
                            context.z0.set(pixel);
                        }
                        context.z.set(context.z0);
                    }

                    float output =
                        ExpressionDeepInteriorPixel;
                    bool decided = false;
                    if (!undefined) {
                        if (mpfr_nan_p(context.z.re) ||
                            mpfr_nan_p(context.z.im)) {
                            undefined = true;
                        } else if (
                            mpfr_inf_p(context.z.re) ||
                            mpfr_inf_p(context.z.im)) {
                            output = 0.0f;
                            decided = true;
                        } else {
                            mpfr_hypot(
                                magnitude,
                                context.z.re, context.z.im,
                                RND);
                            if (mpfr_nan_p(magnitude)) {
                                undefined = true;
                            } else if (
                                mpfr_inf_p(magnitude) ||
                                mpfr_cmp_d(
                                    magnitude,
                                    request.bailout) > 0) {
                                output = 0.0f;
                                decided = true;
                            }
                        }
                    }
                    for (int iteration = 0;
                         !undefined && !decided &&
                         iteration < request.maxIterations;
                         ++iteration) {
                        if (pollCancellation()) break;
                        context.iteration = iteration;
                        std::string oracleError;
                        const bool oracleDefined =
                            ExpressionOracle::evaluate(
                                *request.runtimeProgram,
                                context, next,
                                &oracleError);
                        ++localIterations;
                        if (mpfr_nan_p(next.re) ||
                            mpfr_nan_p(next.im)) {
                            undefined = true;
                            break;
                        }
                        if (mpfr_inf_p(next.re) ||
                            mpfr_inf_p(next.im)) {
                            output = static_cast<float>(
                                iteration + 1);
                            decided = true;
                            break;
                        }
                        if (!oracleDefined ||
                            !mpfr_number_p(next.re) ||
                            !mpfr_number_p(next.im)) {
                            undefined = true;
                            break;
                        }
                        context.z.set(next);
                        mpfr_hypot(
                            magnitude,
                            context.z.re, context.z.im,
                            RND);
                        if (mpfr_nan_p(magnitude)) {
                            undefined = true;
                            break;
                        }
                        if (mpfr_inf_p(magnitude) ||
                            mpfr_cmp_d(
                                magnitude,
                                request.bailout) > 0) {
                            output = static_cast<float>(
                                iteration + 1);
                            decided = true;
                        } else if (iteration + 1 ==
                                       request.maxIterations) {
                            output =
                                ExpressionDeepInteriorPixel;
                            decided = true;
                        }
                    }
                    if (cancelled.load(
                            std::memory_order_acquire))
                        continue;
                    if (undefined) {
                        undefinedPixels.fetch_add(
                            1,
                            std::memory_order_relaxed);
                    } else if (decided) {
                        if (!cancelled.load()) {
                            request.output[outputIndex] =
                                output;
                            if (cancelled.load())
                                request.output[outputIndex] =
                                    ExpressionDeepEmptyPixel;
                        }
                    }
                    const uint64_t completed =
                        fallbackCompleted.fetch_add(
                            1,
                            std::memory_order_relaxed) + 1;
                    if ((completed & 31) == 0)
                        notifyProgress(
                            ExpressionDeepRenderPhase::Fallback,
                            completed,
                            static_cast<uint64_t>(
                                fallbackQueue.size()));
                  } catch (const std::bad_alloc&) {
                    fallbackResourceError.store(
                        true, std::memory_order_release);
                  } catch (const std::length_error&) {
                    fallbackResourceError.store(
                        true, std::memory_order_release);
                  } catch (...) {
                    fallbackInternalError.store(
                        true, std::memory_order_release);
                  }
                }
                totalIterations.fetch_add(
                    localIterations,
                    std::memory_order_relaxed);
            }
        }
        result.fallbackSeconds =
            secondsSince(fallbackStart);
        result.totalIterations =
            totalIterations.load(std::memory_order_relaxed);
        result.undefinedPixelCount =
            undefinedPixels.load(std::memory_order_relaxed);
        if (fallbackResourceError.load(
                std::memory_order_acquire))
            return fail(
                ExpressionDeepRenderStatus::ResourceLimit,
                "fallback worker allocation failed");
        if (fallbackInternalError.load(
                std::memory_order_acquire))
            return fail(
                ExpressionDeepRenderStatus::InternalError,
                "fallback geometry initialization failed");
        if (cancelled.load(std::memory_order_acquire))
            return fail(
                ExpressionDeepRenderStatus::Cancelled,
                "render cancelled during MPFR fallback");
        if (result.undefinedPixelCount != 0)
            return fail(
                ExpressionDeepRenderStatus::UndefinedPixel,
                "one or more pixel orbits are undefined; their output remains empty");

        notifyProgress(
            ExpressionDeepRenderPhase::Complete,
            pixelCount64, pixelCount64);
        result.status = ExpressionDeepRenderStatus::Success;
        result.success = true;
        result.cancelled = false;
        result.error.clear();
        return true;
    } catch (const std::bad_alloc&) {
        return fail(
            ExpressionDeepRenderStatus::ResourceLimit,
            "renderer allocation failed");
    } catch (const std::length_error&) {
        return fail(
            ExpressionDeepRenderStatus::ResourceLimit,
            "renderer allocation exceeds container limits");
    } catch (...) {
        return fail(
            ExpressionDeepRenderStatus::InternalError,
            "renderer failed with an unexpected exception");
    }
}

} // namespace formula
