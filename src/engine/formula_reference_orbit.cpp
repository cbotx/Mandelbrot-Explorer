#include "formula_reference_orbit.h"
#include "formula_scaled_residual.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <new>
#include <stdexcept>

namespace formula {

namespace {

constexpr mpfr_rnd_t RND = MPFR_RNDN;
constexpr mpfr_rnd_t SHADOW_RND = MPFR_RNDZ;

bool finiteComplex(Complex value) {
    return std::isfinite(value.real()) && std::isfinite(value.imag());
}

size_t decimalDigits(const std::string& text) {
    size_t digits = 0;
    size_t end = text.find_first_of("eE");
    if (end == std::string::npos) end = text.size();
    for (size_t i = 0; i < end; ++i) {
        if (text[i] >= '0' && text[i] <= '9') ++digits;
    }
    return digits;
}

bool choosePrecision(const ExpressionReferenceBuildRequest& request, mpfr_prec_t& precision, std::string& error) {
    const ExpressionReferencePrecisionPolicy& policy = request.precision;
    if (policy.requestedBits < 0 || policy.viewBits < 0 || policy.guardBits < 0 || policy.minimumBits < 0 || policy.maximumBits < MPFR_PREC_MIN || policy.maximumBits > MPFR_PREC_MAX) {
        error = "invalid MPFR precision policy";
        return false;
    }

    uint64_t inputBits = 0;
    if (request.center.usesMpf()) {
        auto significantBits = [](mpf_srcptr value) {
            const mp_size_t limbCount = mpf_size(value);
            if (limbCount == 0) return uint64_t{0};
            const size_t highBit = mpn_sizeinbase(value->_mp_d, limbCount, 2);
            const mp_bitcnt_t lowBit = mpn_scan1(value->_mp_d, 0);
            return highBit > lowBit ? static_cast<uint64_t>(highBit - lowBit) : uint64_t{1};
        };
        // GMP may retain a carry limb beyond mpf_get_prec(), while
        // mpf_set_prec_raw() may leave more live bits than the nominal
        // precision. Count the actual nonzero significand span so mpfr_set_f
        // remains exact without charging zero padding as precision.
        inputBits = std::max(significantBits(request.center.realMpf), significantBits(request.center.imaginaryMpf));
    } else {
        size_t digits = std::max(decimalDigits(request.center.realDecimal), decimalDigits(request.center.imaginaryDecimal));
        inputBits = static_cast<uint64_t>(std::ceil(digits * 3.32192809488736234787)) + 2;
    }
    uint64_t base = std::max<uint64_t>({53, static_cast<uint64_t>(policy.requestedBits), static_cast<uint64_t>(policy.viewBits), static_cast<uint64_t>(policy.minimumBits), inputBits});
    uint64_t guard = static_cast<uint64_t>(policy.guardBits);
    if (base > std::numeric_limits<uint64_t>::max() - guard) {
        error = "MPFR precision calculation overflow";
        return false;
    }
    uint64_t selected = base + guard;
    if (selected > static_cast<uint64_t>(policy.maximumBits) || selected > static_cast<uint64_t>(MPFR_PREC_MAX)) {
        error = "required MPFR precision exceeds policy maximum";
        return false;
    }
    precision = static_cast<mpfr_prec_t>(selected);
    return true;
}

bool setExactInput(const ExpressionReferenceExactInput& input, MpfrComplex& value) {
    if (input.usesMpf()) {
        if (!input.realMpf || !input.imaginaryMpf || input.usesDecimal()) return false;
        if (mpfr_set_f(value.re, input.realMpf, RND) != 0 || mpfr_set_f(value.im, input.imaginaryMpf, RND) != 0) return false;
    } else {
        if (!input.usesDecimal()) return false;
        const std::string imaginary = input.imaginaryDecimal.empty() ? "0" : input.imaginaryDecimal;
        if (!value.set(input.realDecimal, imaginary)) return false;
    }
    return mpfr_number_p(value.re) && mpfr_number_p(value.im);
}

void setDoubleContext(const ExpressionContext& input, ExpressionOracleContext& output) {
    output.z.set(input.z.real(), input.z.imag());
    output.c.set(input.c.real(), input.c.imag());
    output.z0.set(input.z0.real(), input.z0.imag());
    for (size_t i = 0; i < input.parameters.size(); ++i) { output.parameters[i].set(input.parameters[i].real(), input.parameters[i].imag()); }
    output.iteration = input.iteration;
}

bool reconstructExpansion(MpfrComplex& output, MpfrComplex& term, const ScaledComplexShadow& primary, const ScaledComplexShadow& defect, const ScaledComplexExpansionTail* tail, size_t termCount, mpfr_rnd_t rounding) {
    if (termCount < 1 || termCount > 4 || !setMpfrFromScaledShadow(output, primary, rounding)) return false;
    const ScaledComplexShadow* residuals[] = {&defect, tail ? &tail->residual2 : nullptr, tail ? &tail->residual3 : nullptr};
    for (size_t index = 1; index < termCount; ++index) {
        if (!residuals[index - 1] || !setMpfrFromScaledShadow(term, *residuals[index - 1], rounding)) return false;
        if (!mpfr_zero_p(term.re)) mpfr_add(output.re, output.re, term.re, rounding);
        if (!mpfr_zero_p(term.im)) mpfr_add(output.im, output.im, term.im, rounding);
    }
    return true;
}

bool compactExpansion(const MpfrComplex& value, size_t termCount, ScaledComplexShadow& primary, ScaledComplexShadow& defect, ScaledComplexExpansionTail* tail, MpfrComplex& reconstructed, MpfrComplex& term, MpfrComplex& difference) {
    primary = {};
    defect = {};
    if (tail) *tail = {};
    if (!makeScaledComplexShadow(value, primary)) return false;
    if (!mpfr_number_p(value.re) || !mpfr_number_p(value.im)) return true;
    ScaledComplexShadow* residuals[] = {&defect, tail ? &tail->residual2 : nullptr, tail ? &tail->residual3 : nullptr};
    for (size_t index = 1; index < termCount; ++index) {
        if (!residuals[index - 1] || !reconstructExpansion(reconstructed, term, primary, defect, tail, index, RND)) return false;
        const mpfr_flags_t saved = mpfr_flags_save();
        mpfr_flags_clear(MPFR_FLAGS_ALL);
        mpfr_sub(difference.re, value.re, reconstructed.re, RND);
        mpfr_sub(difference.im, value.im, reconstructed.im, RND);
        const mpfr_flags_t rangeFlags = mpfr_flags_test(MPFR_FLAGS_UNDERFLOW | MPFR_FLAGS_OVERFLOW | MPFR_FLAGS_ERANGE);
        mpfr_flags_restore(saved, MPFR_FLAGS_ALL);
        if (rangeFlags != 0 || !mpfr_number_p(difference.re) || !mpfr_number_p(difference.im) || !makeScaledComplexShadow(difference, *residuals[index - 1])) return false;
    }
    return true;
}

bool subtractOutward(mpfr_srcptr exactComponent, mpfr_srcptr reconstructedComponent, mpfr_ptr lower, mpfr_ptr upper, mpfr_ptr output) {
    const mpfr_flags_t saved = mpfr_flags_save();
    mpfr_flags_clear(MPFR_FLAGS_ALL);
    mpfr_sub(lower, exactComponent, reconstructedComponent, MPFR_RNDD);
    mpfr_sub(upper, exactComponent, reconstructedComponent, MPFR_RNDU);
    mpfr_abs(lower, lower, MPFR_RNDU);
    mpfr_abs(upper, upper, MPFR_RNDU);
    const mpfr_flags_t rangeFlags = mpfr_flags_test(MPFR_FLAGS_UNDERFLOW | MPFR_FLAGS_OVERFLOW | MPFR_FLAGS_ERANGE);
    mpfr_flags_restore(saved, MPFR_FLAGS_ALL);
    if (rangeFlags != 0 || !mpfr_number_p(lower) || !mpfr_number_p(upper)) return false;
    mpfr_set(output, mpfr_cmp(lower, upper) >= 0 ? lower : upper, MPFR_RNDU);
    return true;
}

bool makeUpwardScaledRadius(mpfr_srcptr value, ScaledRealValue& output) {
    output = {};
    if (!mpfr_number_p(value) || mpfr_sgn(value) < 0) return false;
    if (mpfr_zero_p(value)) return true;
    long exponent = 0;
    double mantissa = mpfr_get_d_2exp(&exponent, value, MPFR_RNDU);
    if (!std::isfinite(mantissa) || mantissa <= 0.0) return false;
    int adjustment = 0;
    mantissa = std::frexp(mantissa, &adjustment);
    if ((adjustment > 0 && exponent > std::numeric_limits<int64_t>::max() - adjustment) || (adjustment < 0 && exponent < std::numeric_limits<int64_t>::min() - adjustment)) return false;
    output.mantissa = mantissa;
    output.exponent = static_cast<int64_t>(exponent) + adjustment;
    return output.isNormalized();
}

bool compactError(const MpfrComplex& exact, const ScaledComplexShadow& primary, const ScaledComplexShadow& defect, const ScaledComplexExpansionTail* tail, MpfrComplex& reconstructedValue, MpfrComplex& reconstructionTerm, MpfrComplex& differenceLower, MpfrComplex& differenceUpper, mpfr_ptr componentMaximum, mpfr_ptr maximum, ScaledRealValue& output) {
    if (tail) {
        if (!reconstructExpansion(reconstructedValue, reconstructionTerm, primary, defect, tail, 4, RND)) return false;
    } else {
        ScaledComplexValue reconstructed;
        if (makeScaledComplexValue(primary, defect, reconstructed) != ScaledArithmeticStatus::Success || !setMpfrFromScaledValue(reconstructedValue, reconstructed)) return false;
    }
    if (!subtractOutward(exact.re, reconstructedValue.re, differenceLower.re, differenceUpper.re, maximum) || !subtractOutward(exact.im, reconstructedValue.im, differenceLower.im, differenceUpper.im, componentMaximum)) return false;
    if (mpfr_cmp(componentMaximum, maximum) > 0) mpfr_set(maximum, componentMaximum, MPFR_RNDU);
    return makeUpwardScaledRadius(maximum, output);
}

bool outsideBailout(const MpfrComplex& value, double bailout, mpfr_ptr magnitude) {
    if (!mpfr_number_p(value.re) || !mpfr_number_p(value.im)) return true;
    mpfr_hypot(magnitude, value.re, value.im, RND);
    return !mpfr_number_p(magnitude) || mpfr_cmp_d(magnitude, bailout) > 0;
}

size_t retainedBytes(const ExpressionReferenceOrbitResult& result) {
    size_t bytes = sizeof(result) + result.samples.capacity() * sizeof(ExpressionReferenceSample) + result.tape.capacity() * sizeof(ExpressionReferenceTapeNode) + result.error.capacity() + result.programSource.capacity();
    if (result.fourTerm) bytes += sizeof(ExpressionReferenceFourTermData) + result.fourTerm->samples.capacity() * sizeof(ExpressionReferenceSampleFourTerm) + result.fourTerm->tape.capacity() * sizeof(ExpressionReferenceTapeNodeFourTerm);
    return bytes;
}

bool addEstimate(uint64_t& total, uint64_t count, uint64_t bytes) {
    if (bytes != 0 && count > (std::numeric_limits<uint64_t>::max() - total) / bytes) return false;
    total += count * bytes;
    return true;
}

bool estimatePeakBytes(const ExpressionProgram& program, uint64_t iterations, uint64_t tapeNodes, mpfr_prec_t precision, ExpressionReferenceCompaction compaction, uint64_t& peakBytes) {
    uint64_t persistent = sizeof(ExpressionReferenceOrbitResult) + 512;
    if (!addEstimate(persistent, iterations, sizeof(ExpressionReferenceSample)) || !addEstimate(persistent, tapeNodes, sizeof(ExpressionReferenceTapeNode)) || !addEstimate(persistent, static_cast<uint64_t>(program.source().size()) + 1, sizeof(char))) return false;
    if (compaction == ExpressionReferenceCompaction::FourTermCertifiedTransfer && (!addEstimate(persistent, 1, sizeof(ExpressionReferenceFourTermData)) || !addEstimate(persistent, iterations, sizeof(ExpressionReferenceSampleFourTerm)) || !addEstimate(persistent, tapeNodes, sizeof(ExpressionReferenceTapeNodeFourTerm)))) return false;

    const uint64_t instructionCount = static_cast<uint64_t>(program.instructionCount());
    const uint64_t stackDepth = static_cast<uint64_t>(program.stackDepth());
    uint64_t transient = sizeof(ExpressionOracleContext) + sizeof(ExpressionOracleTrace) + 6 * sizeof(MpfrComplex);
    if (!addEstimate(transient, stackDepth, sizeof(MpfrComplex)) || !addEstimate(transient, stackDepth, sizeof(uint16_t)) || !addEstimate(transient, instructionCount, sizeof(ExpressionOracleTraceNode))) return false;

    const uint64_t bits = static_cast<uint64_t>(precision);
    const uint64_t limbBits = static_cast<uint64_t>(GMP_NUMB_BITS);
    const uint64_t limbs = (bits + limbBits - 1) / limbBits;
    const uint64_t mpfrValueBytes = (limbs + 2) * sizeof(mp_limb_t) + 32;

    // Base context/compaction values plus oracle stack, scratch, trace
    // output+auxiliary values, operand copies, operation temporaries, and a
    // conservative second working set for MPFR's internal transient storage.
    uint64_t mpfrValues = 64;
    if (!addEstimate(mpfrValues, stackDepth, 2) || !addEstimate(mpfrValues, instructionCount, 4) || !addEstimate(transient, mpfrValues, mpfrValueBytes) || !addEstimate(transient, mpfrValues, mpfrValueBytes)) return false;

    peakBytes = persistent;
    return addEstimate(peakBytes, 1, transient);
}

} // namespace

bool makeScaledRealShadow(mpfr_srcptr value, ScaledRealShadow& output) {
    ScaledRealShadow result{};
    if (mpfr_nan_p(value)) {
        result.mantissa = std::numeric_limits<double>::quiet_NaN();
        output = result;
        return true;
    }
    if (mpfr_inf_p(value)) {
        result.mantissa = mpfr_signbit(value) ? -std::numeric_limits<double>::infinity() : std::numeric_limits<double>::infinity();
        output = result;
        return true;
    }
    if (mpfr_zero_p(value)) {
        result.mantissa = mpfr_signbit(value) ? -0.0 : 0.0;
        output = result;
        return true;
    }
    long exponent = 0;
    result.mantissa = mpfr_get_d_2exp(&exponent, value, SHADOW_RND);
    result.exponent = static_cast<int64_t>(exponent);
    if (!result.isFinite() || result.exponent < static_cast<int64_t>(mpfr_get_emin()) || result.exponent > static_cast<int64_t>(mpfr_get_emax())) return false;
    output = result;
    return true;
}

bool makeScaledComplexShadow(const MpfrComplex& value, ScaledComplexShadow& output) {
    ScaledComplexShadow result{};
    if (!makeScaledRealShadow(value.re, result.re) || !makeScaledRealShadow(value.im, result.im)) return false;
    output = result;
    return true;
}

bool setMpfrFromScaledShadow(mpfr_ptr output, const ScaledRealShadow& shadow, mpfr_rnd_t rounding) {
    if (shadow.isZero()) {
        mpfr_set_zero(output, std::signbit(shadow.mantissa) ? -1 : 1);
        return true;
    }
    if (shadow.isNan()) {
        mpfr_set_nan(output);
        return true;
    }
    if (shadow.isInfinity()) {
        mpfr_set_inf(output, shadow.mantissa < 0.0 ? -1 : 1);
        return true;
    }
    if (!shadow.isFinite() || shadow.exponent < static_cast<int64_t>(std::numeric_limits<mpfr_exp_t>::min()) || shadow.exponent > static_cast<int64_t>(std::numeric_limits<mpfr_exp_t>::max())) return false;
    mpfr_set_d(output, shadow.mantissa, rounding);
    mpfr_exp_t current = mpfr_get_exp(output);
    mpfr_exp_t exponent = static_cast<mpfr_exp_t>(shadow.exponent);
    if ((exponent > 0 && current > std::numeric_limits<mpfr_exp_t>::max() - exponent) || (exponent < 0 && current < std::numeric_limits<mpfr_exp_t>::min() - exponent)) return false;
    const mpfr_exp_t target = current + exponent;
    if (target < mpfr_get_emin() || target > mpfr_get_emax()) return false;
    return mpfr_set_exp(output, target) == 0;
}

bool setMpfrFromScaledShadow(MpfrComplex& output, const ScaledComplexShadow& shadow, mpfr_rnd_t rounding) {
    return setMpfrFromScaledShadow(output.re, shadow.re, rounding) && setMpfrFromScaledShadow(output.im, shadow.im, rounding);
}

bool reconstructMpfrFromShadows(MpfrComplex& output, const ScaledComplexShadow& primary, const ScaledComplexShadow& defect, mpfr_rnd_t rounding) {
    return reconstructMpfrFromExpansion(output, primary, defect, nullptr, rounding);
}

bool reconstructMpfrFromExpansion(MpfrComplex& output, const ScaledComplexShadow& primary, const ScaledComplexShadow& defect, const ScaledComplexExpansionTail* tail, mpfr_rnd_t rounding) {
    MpfrComplex term(output.precision());
    return reconstructExpansion(output, term, primary, defect, tail, tail ? 4 : 2, rounding);
}

static bool buildExpressionReferenceOrbitImpl(const ExpressionReferenceBuildRequest& request, ExpressionReferenceOrbitResult& result) {
    result = {};
    result.compaction = request.compaction;
    auto fail = [&](ExpressionReferenceBuildStatus status, const char* message) {
        result.status = status;
        result.error = message;
        result.memoryBytes = retainedBytes(result);
        return false;
    };
    auto succeed = [&] {
        result.status = ExpressionReferenceBuildStatus::Success;
        result.valid = true;
        result.sampleCount = result.samples.size();
        result.memoryBytes = retainedBytes(result);
        return true;
    };

    if ((!request.canonicalProgram && !request.runtimeProgram) || (request.pixelParameter != FormulaParameter::C && request.pixelParameter != FormulaParameter::InitialZ) || request.maxIterations < 0 || !(request.bailout > 0.0) || !std::isfinite(request.bailout) || (request.center.usesMpf() == request.center.usesDecimal()) || (request.center.usesDecimal() && request.center.realDecimal.empty()) || (request.compaction != ExpressionReferenceCompaction::TwoTerm && request.compaction != ExpressionReferenceCompaction::FourTermCertifiedTransfer) || (request.compaction == ExpressionReferenceCompaction::FourTermCertifiedTransfer && request.certificationPrecision == 0)) { return fail(ExpressionReferenceBuildStatus::InvalidRequest, "invalid expression reference request"); }
    if ((request.center.usesMpf() && (!request.center.realMpf || !request.center.imaginaryMpf)) || !finiteComplex(request.fixed.z0) || !finiteComplex(request.fixed.c) || !std::all_of(request.fixed.parameters.begin(), request.fixed.parameters.end(), finiteComplex)) { return fail(ExpressionReferenceBuildStatus::InvalidRequest, "reference inputs must be finite and complete"); }

    ExpressionProgram specialized;
    const ExpressionProgram* runtime = request.runtimeProgram;
    if (request.canonicalProgram) {
        if (!request.canonicalProgram->valid()) { return fail(ExpressionReferenceBuildStatus::InvalidRequest, "canonical expression program is invalid"); }
        ExpressionError expressionError;
        if (!request.canonicalProgram->specialize(request.fixed, request.pixelParameter, specialized, &expressionError)) {
            result.status = ExpressionReferenceBuildStatus::InvalidRequest;
            result.error = expressionError.message;
            result.memoryBytes = retainedBytes(result);
            return false;
        }
        if (runtime) {
            if (!runtime->valid() || runtime->source() != request.canonicalProgram->source() || !runtime->semanticallyEquivalent(specialized)) { return fail(ExpressionReferenceBuildStatus::ProgramMismatch, "runtime expression does not match canonical specialization"); }
        } else {
            runtime = &specialized;
        }
    }
    if (!runtime || !runtime->valid()) { return fail(ExpressionReferenceBuildStatus::InvalidRequest, "runtime expression program is invalid"); }
    if (runtime->containsOrbitInvariant()) { return fail(ExpressionReferenceBuildStatus::UnsupportedProgram, "orbit-plan bytecode is not a reference-orbit program"); }

    mpfr_prec_t precision = 0;
    std::string precisionError;
    if (!choosePrecision(request, precision, precisionError)) {
        result.status = ExpressionReferenceBuildStatus::PrecisionOutOfRange;
        result.error = precisionError;
        result.memoryBytes = retainedBytes(result);
        return false;
    }

    result.precision = precision;
    result.programSemanticHash = runtime->semanticHash();
    result.canonicalSemanticHash = request.canonicalProgram ? request.canonicalProgram->semanticHash() : 0;
    result.programSource = runtime->source();

    if (request.shouldCancel && request.shouldCancel()) {
        result.cancelled = true;
        return succeed();
    }

    const uint64_t instructionCount = static_cast<uint64_t>(runtime->instructionCount());
    const uint64_t iterations = static_cast<uint64_t>(request.maxIterations);
    if (instructionCount != 0 && iterations > std::numeric_limits<uint64_t>::max() / instructionCount) { return fail(ExpressionReferenceBuildStatus::ResourceLimit, "reference tape size overflow"); }
    const uint64_t tapeNodes = instructionCount * iterations;
    if (iterations > static_cast<uint64_t>(result.samples.max_size()) || tapeNodes > static_cast<uint64_t>(result.tape.max_size())) { return fail(ExpressionReferenceBuildStatus::ResourceLimit, "reference orbit exceeds vector capacity"); }
    if (precision > ExpressionReferencePrecisionPolicy::ApplicationMaximumBits) { return fail(ExpressionReferenceBuildStatus::ResourceLimit, "selected precision exceeds application resource cap"); }
    uint64_t estimatedPeak = 0;
    if (!estimatePeakBytes(*runtime, iterations, tapeNodes, precision, request.compaction, estimatedPeak) || (request.memoryLimitBytes != 0 && estimatedPeak > static_cast<uint64_t>(request.memoryLimitBytes))) { return fail(ExpressionReferenceBuildStatus::ResourceLimit, "reference orbit exceeds peak memory limit"); }
    if (request.certificationPrecision != 0) {
        if (request.certificationPrecision <= precision || request.certificationPrecision > request.precision.maximumBits || request.certificationPrecision > ExpressionReferencePrecisionPolicy::ApplicationMaximumBits) { return fail(ExpressionReferenceBuildStatus::PrecisionOutOfRange, "higher reference certification precision is invalid"); }
        uint64_t certificationPeak = 0;
        if (!estimatePeakBytes(*runtime, iterations, tapeNodes, request.certificationPrecision, request.compaction, certificationPeak) || estimatedPeak > std::numeric_limits<uint64_t>::max() - certificationPeak || (request.memoryLimitBytes != 0 && estimatedPeak + certificationPeak > static_cast<uint64_t>(request.memoryLimitBytes))) { return fail(ExpressionReferenceBuildStatus::ResourceLimit, "certified reference exceeds peak memory limit"); }
        result.certificationPrecision = request.certificationPrecision;
    }

    ExpressionOracleContext context(precision);
    setDoubleContext(request.fixed, context);
    MpfrComplex center(precision);
    if (!setExactInput(request.center, center)) { return fail(ExpressionReferenceBuildStatus::InputParseError, "failed to parse exact reference center"); }
    if (request.pixelParameter == FormulaParameter::C)
        context.c.set(center);
    else
        context.z0.set(center);
    context.z.set(context.z0);

    const bool certify = request.certificationPrecision != 0;
    const mpfr_prec_t certificationPrecision = certify ? request.certificationPrecision : precision;
    ExpressionOracleContext certificationContext(certificationPrecision);
    setDoubleContext(request.fixed, certificationContext);
    MpfrComplex certificationCenter(certificationPrecision);
    if (!setExactInput(request.center, certificationCenter)) { return fail(ExpressionReferenceBuildStatus::InputParseError, "failed to parse higher-precision reference center"); }
    if (request.pixelParameter == FormulaParameter::C)
        certificationContext.c.set(certificationCenter);
    else
        certificationContext.z0.set(certificationCenter);
    certificationContext.z.set(certificationContext.z0);

    std::shared_ptr<ExpressionReferenceFourTermData> fourTerm;
    if (request.compaction == ExpressionReferenceCompaction::FourTermCertifiedTransfer) {
        fourTerm = std::make_shared<ExpressionReferenceFourTermData>();
        result.fourTerm = fourTerm;
    }

    MpfrComplex reconstructed(precision);
    MpfrComplex reconstructionTerm(precision);
    MpfrComplex difference(precision);
    MpfrComplex radiusReconstructed(certificationPrecision);
    MpfrComplex radiusReconstructionTerm(certificationPrecision);
    MpfrComplex radiusDifferenceLower(certificationPrecision);
    MpfrComplex radiusDifferenceUpper(certificationPrecision);
    MpfrComplex radiusMaximumStorage(certificationPrecision);
    mpfr_ptr radiusComponentMaximum = radiusMaximumStorage.re;
    mpfr_ptr radiusMaximum = radiusMaximumStorage.im;
    auto compact = [&](const MpfrComplex& value, const MpfrComplex& exact, ScaledComplexShadow& shadow, ScaledComplexShadow& defect, ScaledComplexExpansionTail* tail, ScaledRealValue& error) {
        shadow = {};
        defect = {};
        if (tail) *tail = {};
        error = {};
        const size_t termCount = tail ? 4 : 2;
        if (!compactExpansion(value, termCount, shadow, defect, tail, reconstructed, reconstructionTerm, difference)) return false;
        return !certify || compactError(exact, shadow, defect, tail, radiusReconstructed, radiusReconstructionTerm, radiusDifferenceLower, radiusDifferenceUpper, radiusComponentMaximum, radiusMaximum, error);
    };
    if (!compact(context.c, certificationContext.c, result.c, result.cDefect, fourTerm ? &fourTerm->c : nullptr, result.cError) || !compact(context.z0, certificationContext.z0, result.z0, result.z0Defect, fourTerm ? &fourTerm->z0 : nullptr, result.z0Error) || !compact(center, certificationCenter, result.pixel, result.pixelDefect, fourTerm ? &fourTerm->pixel : nullptr, result.pixelError)) { return fail(ExpressionReferenceBuildStatus::CompactionOutOfRange, "reference input cannot be represented by compact shadows"); }
    result.initialZ = result.z0;
    result.initialZDefect = result.z0Defect;
    result.initialZError = result.z0Error;
    if (fourTerm) fourTerm->initialZ = fourTerm->z0;

    result.samples.reserve(static_cast<size_t>(iterations));
    result.tape.reserve(static_cast<size_t>(tapeNodes));
    if (fourTerm) {
        fourTerm->samples.reserve(static_cast<size_t>(iterations));
        fourTerm->tape.reserve(static_cast<size_t>(tapeNodes));
    }

    MpfrComplex magnitudeStorage(precision);
    mpfr_ptr magnitude = magnitudeStorage.re;
    if (outsideBailout(context.z, request.bailout, magnitude)) {
        result.escaped = true;
        result.escapeIteration = 0;
    } else {
        MpfrComplex next(precision);
        MpfrComplex certificationNext(certificationPrecision);
        for (int iteration = 0; iteration < request.maxIterations; ++iteration) {
            if (request.shouldCancel && request.shouldCancel()) {
                result.cancelled = true;
                break;
            }

            context.iteration = iteration;
            ExpressionOracleTrace trace;
            std::string oracleError;
            bool defined = ExpressionOracle::evaluateTrace(*runtime, context, next, trace, &oracleError);
            ExpressionOracleTrace certificationTrace;
            std::string certificationError;
            bool certificationDefined = true;
            if (certify) {
                certificationContext.iteration = iteration;
                certificationDefined = ExpressionOracle::evaluateTrace(*runtime, certificationContext, certificationNext, certificationTrace, &certificationError);
            } else {
                certificationNext.set(next);
            }
            if (trace.nodes.empty() || trace.nodes.size() > std::numeric_limits<uint16_t>::max()) { return fail(ExpressionReferenceBuildStatus::UnsupportedProgram, "oracle trace did not produce a compact program tape"); }
            if (certify && (certificationTrace.nodes.size() != trace.nodes.size() || certificationDefined != defined)) { return fail(ExpressionReferenceBuildStatus::CompactionOutOfRange, "higher-precision reference did not converge to the same finite trace"); }

            ExpressionReferenceSample sample;
            ExpressionReferenceSampleFourTerm sampleFourTerm;
            sample.iteration = iteration;
            if (!compact(context.z, certificationContext.z, sample.z, sample.zDefect, fourTerm ? &sampleFourTerm.z : nullptr, sample.zError) || !compact(next, certificationNext, sample.next, sample.rootDefect, fourTerm ? &sampleFourTerm.next : nullptr, sample.nextError)) { return fail(ExpressionReferenceBuildStatus::CompactionOutOfRange, "reference sample cannot be represented by compact shadows"); }
            sample.tapeOffset = result.tape.size();
            sample.tapeCount = static_cast<uint16_t>(trace.nodes.size());
            sample.rootNode = static_cast<uint16_t>(trace.nodes.size() - 1);
            for (size_t traceIndex = 0; traceIndex < trace.nodes.size(); ++traceIndex) {
                const ExpressionOracleTraceNode& traced = trace.nodes[traceIndex];
                const ExpressionOracleTraceNode& exact = certify ? certificationTrace.nodes[traceIndex] : traced;
                if (certify && (traced.operation != exact.operation || traced.argument != exact.argument || traced.leftNode != exact.leftNode || traced.rightNode != exact.rightNode)) { return fail(ExpressionReferenceBuildStatus::CompactionOutOfRange, "higher-precision reference tape layout mismatch"); }
                ExpressionReferenceTapeNode node;
                ExpressionReferenceTapeNodeFourTerm nodeFourTerm;
                if (!compact(traced.output, exact.output, node.output, node.outputDefect, fourTerm ? &nodeFourTerm.output : nullptr, node.outputError)) { return fail(ExpressionReferenceBuildStatus::CompactionOutOfRange, "reference tape output cannot be represented by compact shadows"); }
                if (traced.flags & (OracleTraceHasCompanion | OracleTraceHasDenominator | OracleTraceHasLogarithmBase)) {
                    if (!compact(traced.auxiliary, exact.auxiliary, node.auxiliary, node.auxiliaryDefect, fourTerm ? &nodeFourTerm.auxiliary : nullptr, node.auxiliaryError)) { return fail(ExpressionReferenceBuildStatus::CompactionOutOfRange, "reference tape auxiliary cannot be represented by compact shadows"); }
                }
                node.leftNode = traced.leftNode;
                node.rightNode = traced.rightNode;
                node.flags = traced.flags;
                node.argument = traced.argument;
                node.operation = traced.operation;
                node.cut = traced.cut;
                node.clearance = traced.clearance;
                node.certification = traced.certification;
                result.tape.push_back(node);
                if (fourTerm) fourTerm->tape.push_back(nodeFourTerm);
            }
            sample.rootError = result.tape[static_cast<size_t>(sample.tapeOffset) + sample.rootNode].outputError;
            result.samples.push_back(sample);
            if (fourTerm) fourTerm->samples.push_back(sampleFourTerm);

            if (!defined) {
                result.undefined = true;
                result.undefinedIteration = iteration + 1;
                result.error = oracleError;
                break;
            }
            context.z.set(next);
            if (certify) certificationContext.z.set(certificationNext);
            if (outsideBailout(context.z, request.bailout, magnitude)) {
                result.escaped = true;
                result.escapeIteration = iteration + 1;
                break;
            }
        }
    }
    result.certifiedAgainstHigherPrecision = certify;
    return succeed();
}

bool buildExpressionReferenceOrbit(const ExpressionReferenceBuildRequest& request, ExpressionReferenceOrbitResult& result) {
    try {
        return buildExpressionReferenceOrbitImpl(request, result);
    } catch (const std::bad_alloc&) {
        std::vector<ExpressionReferenceSample>().swap(result.samples);
        std::vector<ExpressionReferenceTapeNode>().swap(result.tape);
        result.fourTerm.reset();
        result.status = ExpressionReferenceBuildStatus::ResourceLimit;
        result.error.clear();
        result.valid = false;
        result.sampleCount = 0;
        result.memoryBytes = retainedBytes(result);
        return false;
    } catch (const std::length_error&) {
        std::vector<ExpressionReferenceSample>().swap(result.samples);
        std::vector<ExpressionReferenceTapeNode>().swap(result.tape);
        result.fourTerm.reset();
        result.status = ExpressionReferenceBuildStatus::ResourceLimit;
        result.error.clear();
        result.valid = false;
        result.sampleCount = 0;
        result.memoryBytes = retainedBytes(result);
        return false;
    }
}

} // namespace formula
