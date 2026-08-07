#include "formula_reference_orbit.h"

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
    return std::isfinite(value.real()) &&
           std::isfinite(value.imag());
}

size_t decimalDigits(const std::string& text) {
    size_t digits = 0;
    size_t end = text.find_first_of("eE");
    if (end == std::string::npos) end = text.size();
    for (size_t i = 0; i < end; ++i) {
        if (text[i] >= '0' && text[i] <= '9')
            ++digits;
    }
    return digits;
}

bool choosePrecision(
        const ExpressionReferenceBuildRequest& request,
        mpfr_prec_t& precision, std::string& error) {
    const ExpressionReferencePrecisionPolicy& policy =
        request.precision;
    if (policy.requestedBits < 0 || policy.viewBits < 0 ||
        policy.guardBits < 0 || policy.minimumBits < 0 ||
        policy.maximumBits < MPFR_PREC_MIN ||
        policy.maximumBits > MPFR_PREC_MAX) {
        error = "invalid MPFR precision policy";
        return false;
    }

    uint64_t inputBits = 0;
    if (request.center.usesMpf()) {
        const uint64_t limbBits =
            static_cast<uint64_t>(GMP_NUMB_BITS);
        const uint64_t realLimbs =
            static_cast<uint64_t>(
                mpf_size(request.center.realMpf));
        const uint64_t imaginaryLimbs =
            static_cast<uint64_t>(
                mpf_size(request.center.imaginaryMpf));
        const uint64_t usedLimbs =
            std::max(realLimbs, imaginaryLimbs);
        if (usedLimbs >
            std::numeric_limits<uint64_t>::max() /
                limbBits) {
            error = "GMP input precision calculation overflow";
            return false;
        }
        // mpf can retain a carry limb beyond mpf_get_prec(), and
        // mpf_set_prec_raw() can make the nominal precision smaller than the
        // live significand. Cover every used limb so mpfr_set_f is exact.
        inputBits = usedLimbs * limbBits;
    } else {
        size_t digits = std::max(
            decimalDigits(request.center.realDecimal),
            decimalDigits(request.center.imaginaryDecimal));
        inputBits = static_cast<uint64_t>(
            std::ceil(digits * 3.32192809488736234787)) + 2;
    }
    uint64_t base = std::max<uint64_t>({
        53,
        static_cast<uint64_t>(policy.requestedBits),
        static_cast<uint64_t>(policy.viewBits),
        static_cast<uint64_t>(policy.minimumBits),
        inputBits
    });
    uint64_t guard = static_cast<uint64_t>(policy.guardBits);
    if (base > std::numeric_limits<uint64_t>::max() - guard) {
        error = "MPFR precision calculation overflow";
        return false;
    }
    uint64_t selected = base + guard;
    if (selected >
            static_cast<uint64_t>(policy.maximumBits) ||
        selected > static_cast<uint64_t>(MPFR_PREC_MAX)) {
        error = "required MPFR precision exceeds policy maximum";
        return false;
    }
    precision = static_cast<mpfr_prec_t>(selected);
    return true;
}

bool setExactInput(
        const ExpressionReferenceExactInput& input,
        MpfrComplex& value) {
    if (input.usesMpf()) {
        if (!input.realMpf || !input.imaginaryMpf ||
            input.usesDecimal())
            return false;
        if (mpfr_set_f(value.re, input.realMpf, RND) != 0 ||
            mpfr_set_f(value.im, input.imaginaryMpf, RND) != 0)
            return false;
    } else {
        if (!input.usesDecimal())
            return false;
        const std::string imaginary =
            input.imaginaryDecimal.empty()
                ? "0" : input.imaginaryDecimal;
        if (!value.set(input.realDecimal, imaginary))
            return false;
    }
    return mpfr_number_p(value.re) &&
           mpfr_number_p(value.im);
}

void setDoubleContext(
        const ExpressionContext& input,
        ExpressionOracleContext& output) {
    output.z.set(input.z.real(), input.z.imag());
    output.c.set(input.c.real(), input.c.imag());
    output.z0.set(input.z0.real(), input.z0.imag());
    for (size_t i = 0; i < input.parameters.size(); ++i) {
        output.parameters[i].set(
            input.parameters[i].real(),
            input.parameters[i].imag());
    }
    output.iteration = input.iteration;
}

bool quantizationDefect(
        const MpfrComplex& value,
        const ScaledComplexShadow& shadow,
        MpfrComplex& reconstructed,
        MpfrComplex& difference,
        ScaledComplexShadow& defect) {
    defect = {};
    if (!mpfr_number_p(value.re) ||
        !mpfr_number_p(value.im))
        return true;
    if (!setMpfrFromScaledShadow(reconstructed, shadow))
        return false;

    const mpfr_flags_t saved = mpfr_flags_save();
    mpfr_flags_clear(MPFR_FLAGS_ALL);
    mpfr_sub(
        difference.re, value.re, reconstructed.re, RND);
    mpfr_sub(
        difference.im, value.im, reconstructed.im, RND);
    const mpfr_flags_t rangeFlags = mpfr_flags_test(
        MPFR_FLAGS_UNDERFLOW |
        MPFR_FLAGS_OVERFLOW |
        MPFR_FLAGS_ERANGE);
    mpfr_flags_restore(saved, MPFR_FLAGS_ALL);
    if (rangeFlags != 0 ||
        !mpfr_number_p(difference.re) ||
        !mpfr_number_p(difference.im))
        return false;
    return makeScaledComplexShadow(difference, defect);
}

bool outsideBailout(
        const MpfrComplex& value, double bailout,
        mpfr_ptr magnitude) {
    if (!mpfr_number_p(value.re) ||
        !mpfr_number_p(value.im))
        return true;
    mpfr_hypot(magnitude, value.re, value.im, RND);
    return !mpfr_number_p(magnitude) ||
           mpfr_cmp_d(magnitude, bailout) > 0;
}

size_t retainedBytes(
        const ExpressionReferenceOrbitResult& result) {
    return sizeof(result) +
           result.samples.capacity() *
               sizeof(ExpressionReferenceSample) +
           result.tape.capacity() *
               sizeof(ExpressionReferenceTapeNode) +
           result.error.capacity() +
           result.programSource.capacity();
}

bool addEstimate(
        uint64_t& total, uint64_t count,
        uint64_t bytes) {
    if (bytes != 0 &&
        count >
            (std::numeric_limits<uint64_t>::max() -
             total) / bytes)
        return false;
    total += count * bytes;
    return true;
}

bool estimatePeakBytes(
        const ExpressionProgram& program,
        uint64_t iterations, uint64_t tapeNodes,
        mpfr_prec_t precision, uint64_t& peakBytes) {
    uint64_t persistent =
        sizeof(ExpressionReferenceOrbitResult) + 512;
    if (!addEstimate(
            persistent, iterations,
            sizeof(ExpressionReferenceSample)) ||
        !addEstimate(
            persistent, tapeNodes,
            sizeof(ExpressionReferenceTapeNode)) ||
        !addEstimate(
            persistent,
            static_cast<uint64_t>(program.source().size()) + 1,
            sizeof(char)))
        return false;

    const uint64_t instructionCount =
        static_cast<uint64_t>(program.instructionCount());
    const uint64_t stackDepth =
        static_cast<uint64_t>(program.stackDepth());
    uint64_t transient =
        sizeof(ExpressionOracleContext) +
        sizeof(ExpressionOracleTrace) +
        5 * sizeof(MpfrComplex);
    if (!addEstimate(
            transient, stackDepth,
            sizeof(MpfrComplex)) ||
        !addEstimate(
            transient, stackDepth,
            sizeof(uint16_t)) ||
        !addEstimate(
            transient, instructionCount,
            sizeof(ExpressionOracleTraceNode)))
        return false;

    const uint64_t bits =
        static_cast<uint64_t>(precision);
    const uint64_t limbBits =
        static_cast<uint64_t>(GMP_NUMB_BITS);
    const uint64_t limbs =
        (bits + limbBits - 1) / limbBits;
    const uint64_t mpfrValueBytes =
        (limbs + 2) * sizeof(mp_limb_t) + 32;

    // Base context/compaction values plus oracle stack, scratch, trace
    // output+auxiliary values, operand copies, operation temporaries, and a
    // conservative second working set for MPFR's internal transient storage.
    uint64_t mpfrValues = 64;
    if (!addEstimate(mpfrValues, stackDepth, 2) ||
        !addEstimate(mpfrValues, instructionCount, 4) ||
        !addEstimate(
            transient, mpfrValues,
            mpfrValueBytes) ||
        !addEstimate(
            transient, mpfrValues,
            mpfrValueBytes))
        return false;

    peakBytes = persistent;
    return addEstimate(peakBytes, 1, transient);
}

} // namespace

bool makeScaledRealShadow(
        mpfr_srcptr value, ScaledRealShadow& output) {
    ScaledRealShadow result{};
    if (mpfr_nan_p(value)) {
        result.mantissa =
            std::numeric_limits<double>::quiet_NaN();
        output = result;
        return true;
    }
    if (mpfr_inf_p(value)) {
        result.mantissa = mpfr_signbit(value)
            ? -std::numeric_limits<double>::infinity()
            : std::numeric_limits<double>::infinity();
        output = result;
        return true;
    }
    if (mpfr_zero_p(value)) {
        result.mantissa =
            mpfr_signbit(value) ? -0.0 : 0.0;
        output = result;
        return true;
    }
    long exponent = 0;
    result.mantissa =
        mpfr_get_d_2exp(&exponent, value, SHADOW_RND);
    result.exponent = static_cast<int64_t>(exponent);
    if (!result.isFinite() ||
        result.exponent <
            static_cast<int64_t>(mpfr_get_emin()) ||
        result.exponent >
            static_cast<int64_t>(mpfr_get_emax()))
        return false;
    output = result;
    return true;
}

bool makeScaledComplexShadow(
        const MpfrComplex& value,
        ScaledComplexShadow& output) {
    ScaledComplexShadow result{};
    if (!makeScaledRealShadow(value.re, result.re) ||
        !makeScaledRealShadow(value.im, result.im))
        return false;
    output = result;
    return true;
}

bool setMpfrFromScaledShadow(
        mpfr_ptr output, const ScaledRealShadow& shadow,
        mpfr_rnd_t rounding) {
    if (shadow.isZero()) {
        mpfr_set_zero(
            output, std::signbit(shadow.mantissa) ? -1 : 1);
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
    if (!shadow.isFinite() ||
        shadow.exponent <
            static_cast<int64_t>(
                std::numeric_limits<mpfr_exp_t>::min()) ||
        shadow.exponent >
            static_cast<int64_t>(
                std::numeric_limits<mpfr_exp_t>::max()))
        return false;
    mpfr_set_d(output, shadow.mantissa, rounding);
    mpfr_exp_t current = mpfr_get_exp(output);
    mpfr_exp_t exponent =
        static_cast<mpfr_exp_t>(shadow.exponent);
    if ((exponent > 0 &&
         current >
             std::numeric_limits<mpfr_exp_t>::max() - exponent) ||
        (exponent < 0 &&
         current <
             std::numeric_limits<mpfr_exp_t>::min() - exponent))
        return false;
    const mpfr_exp_t target = current + exponent;
    if (target < mpfr_get_emin() ||
        target > mpfr_get_emax())
        return false;
    return mpfr_set_exp(output, target) == 0;
}

bool setMpfrFromScaledShadow(
        MpfrComplex& output,
        const ScaledComplexShadow& shadow,
        mpfr_rnd_t rounding) {
    return setMpfrFromScaledShadow(
               output.re, shadow.re, rounding) &&
           setMpfrFromScaledShadow(
               output.im, shadow.im, rounding);
}

bool reconstructMpfrFromShadows(
        MpfrComplex& output,
        const ScaledComplexShadow& primary,
        const ScaledComplexShadow& defect,
        mpfr_rnd_t rounding) {
    MpfrComplex remainder(output.precision());
    if (!setMpfrFromScaledShadow(
            output, primary, rounding) ||
        !setMpfrFromScaledShadow(
            remainder, defect, rounding))
        return false;
    mpfr_add(output.re, output.re, remainder.re, rounding);
    mpfr_add(output.im, output.im, remainder.im, rounding);
    return true;
}

static bool buildExpressionReferenceOrbitImpl(
        const ExpressionReferenceBuildRequest& request,
        ExpressionReferenceOrbitResult& result) {
    result = {};
    auto fail = [&](ExpressionReferenceBuildStatus status,
                    const char* message) {
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

    if ((!request.canonicalProgram &&
         !request.runtimeProgram) ||
        (request.pixelParameter != FormulaParameter::C &&
         request.pixelParameter != FormulaParameter::InitialZ) ||
        request.maxIterations < 0 ||
        !(request.bailout > 0.0) ||
        !std::isfinite(request.bailout) ||
        (request.center.usesMpf() ==
         request.center.usesDecimal()) ||
        (request.center.usesDecimal() &&
         request.center.realDecimal.empty())) {
        return fail(
            ExpressionReferenceBuildStatus::InvalidRequest,
            "invalid expression reference request");
    }
    if ((request.center.usesMpf() &&
         (!request.center.realMpf ||
          !request.center.imaginaryMpf)) ||
        !finiteComplex(request.fixed.z0) ||
        !finiteComplex(request.fixed.c) ||
        !std::all_of(
            request.fixed.parameters.begin(),
            request.fixed.parameters.end(),
            finiteComplex)) {
        return fail(
            ExpressionReferenceBuildStatus::InvalidRequest,
            "reference inputs must be finite and complete");
    }

    ExpressionProgram specialized;
    const ExpressionProgram* runtime =
        request.runtimeProgram;
    if (request.canonicalProgram) {
        if (!request.canonicalProgram->valid()) {
            return fail(
                ExpressionReferenceBuildStatus::InvalidRequest,
                "canonical expression program is invalid");
        }
        ExpressionError expressionError;
        if (!request.canonicalProgram->specialize(
                request.fixed, request.pixelParameter,
                specialized, &expressionError)) {
            result.status =
                ExpressionReferenceBuildStatus::InvalidRequest;
            result.error = expressionError.message;
            result.memoryBytes = retainedBytes(result);
            return false;
        }
        if (runtime) {
            if (!runtime->valid() ||
                runtime->source() !=
                    request.canonicalProgram->source() ||
                !runtime->semanticallyEquivalent(specialized)) {
                return fail(
                    ExpressionReferenceBuildStatus::ProgramMismatch,
                    "runtime expression does not match canonical specialization");
            }
        } else {
            runtime = &specialized;
        }
    }
    if (!runtime || !runtime->valid()) {
        return fail(
            ExpressionReferenceBuildStatus::InvalidRequest,
            "runtime expression program is invalid");
    }
    if (runtime->containsOrbitInvariant()) {
        return fail(
            ExpressionReferenceBuildStatus::UnsupportedProgram,
            "orbit-plan bytecode is not a reference-orbit program");
    }

    mpfr_prec_t precision = 0;
    std::string precisionError;
    if (!choosePrecision(
            request, precision, precisionError)) {
        result.status =
            ExpressionReferenceBuildStatus::PrecisionOutOfRange;
        result.error = precisionError;
        result.memoryBytes = retainedBytes(result);
        return false;
    }

    result.precision = precision;
    result.programSemanticHash = runtime->semanticHash();
    result.canonicalSemanticHash =
        request.canonicalProgram
            ? request.canonicalProgram->semanticHash()
            : 0;
    result.programSource = runtime->source();

    if (request.shouldCancel &&
        request.shouldCancel()) {
        result.cancelled = true;
        return succeed();
    }

    const uint64_t instructionCount =
        static_cast<uint64_t>(runtime->instructionCount());
    const uint64_t iterations =
        static_cast<uint64_t>(request.maxIterations);
    if (instructionCount != 0 &&
        iterations >
            std::numeric_limits<uint64_t>::max() /
                instructionCount) {
        return fail(
            ExpressionReferenceBuildStatus::ResourceLimit,
            "reference tape size overflow");
    }
    const uint64_t tapeNodes =
        instructionCount * iterations;
    if (iterations >
            static_cast<uint64_t>(
                result.samples.max_size()) ||
        tapeNodes >
            static_cast<uint64_t>(
                result.tape.max_size())) {
        return fail(
            ExpressionReferenceBuildStatus::ResourceLimit,
            "reference orbit exceeds vector capacity");
    }
    if (precision >
            ExpressionReferencePrecisionPolicy::
                ApplicationMaximumBits) {
        return fail(
            ExpressionReferenceBuildStatus::ResourceLimit,
            "selected precision exceeds application resource cap");
    }
    uint64_t estimatedPeak = 0;
    if (!estimatePeakBytes(
            *runtime, iterations, tapeNodes,
            precision, estimatedPeak) ||
        (request.memoryLimitBytes != 0 &&
         estimatedPeak >
            static_cast<uint64_t>(
                request.memoryLimitBytes))) {
        return fail(
            ExpressionReferenceBuildStatus::ResourceLimit,
            "reference orbit exceeds peak memory limit");
    }

    ExpressionOracleContext context(precision);
    setDoubleContext(request.fixed, context);
    MpfrComplex center(precision);
    if (!setExactInput(request.center, center)) {
        return fail(
            ExpressionReferenceBuildStatus::InputParseError,
            "failed to parse exact reference center");
    }
    if (request.pixelParameter == FormulaParameter::C)
        context.c.set(center);
    else
        context.z0.set(center);
    context.z.set(context.z0);

    MpfrComplex reconstructed(precision);
    MpfrComplex difference(precision);
    auto compact = [&](const MpfrComplex& value,
                       ScaledComplexShadow& shadow,
                       ScaledComplexShadow& defect) {
        shadow = {};
        defect = {};
        return makeScaledComplexShadow(value, shadow) &&
               quantizationDefect(
                   value, shadow, reconstructed,
                   difference, defect);
    };
    if (!compact(context.c, result.c, result.cDefect) ||
        !compact(context.z0, result.z0, result.z0Defect) ||
        !compact(center, result.pixel, result.pixelDefect)) {
        return fail(
            ExpressionReferenceBuildStatus::
                CompactionOutOfRange,
            "reference input cannot be represented by compact shadows");
    }
    result.initialZ = result.z0;
    result.initialZDefect = result.z0Defect;

    result.samples.reserve(
        static_cast<size_t>(iterations));
    result.tape.reserve(
        static_cast<size_t>(tapeNodes));

    MpfrComplex magnitudeStorage(precision);
    mpfr_ptr magnitude = magnitudeStorage.re;
    if (outsideBailout(
            context.z, request.bailout, magnitude)) {
        result.escaped = true;
        result.escapeIteration = 0;
    } else {
        MpfrComplex next(precision);
        for (int iteration = 0;
             iteration < request.maxIterations; ++iteration) {
            if (request.shouldCancel &&
                request.shouldCancel()) {
                result.cancelled = true;
                break;
            }

            context.iteration = iteration;
            ExpressionOracleTrace trace;
            std::string oracleError;
            bool defined = ExpressionOracle::evaluateTrace(
                *runtime, context, next, trace, &oracleError);
            if (trace.nodes.empty() ||
                trace.nodes.size() >
                    std::numeric_limits<uint16_t>::max()) {
                return fail(
                    ExpressionReferenceBuildStatus::UnsupportedProgram,
                    "oracle trace did not produce a compact program tape");
            }

            ExpressionReferenceSample sample;
            sample.iteration = iteration;
            if (!compact(
                    context.z, sample.z,
                    sample.zDefect) ||
                !compact(
                    next, sample.next,
                    sample.rootDefect)) {
                return fail(
                    ExpressionReferenceBuildStatus::
                        CompactionOutOfRange,
                    "reference sample cannot be represented by compact shadows");
            }
            sample.tapeOffset = result.tape.size();
            sample.tapeCount =
                static_cast<uint16_t>(trace.nodes.size());
            sample.rootNode =
                static_cast<uint16_t>(trace.nodes.size() - 1);
            for (const ExpressionOracleTraceNode& traced :
                 trace.nodes) {
                ExpressionReferenceTapeNode node;
                if (!compact(
                        traced.output, node.output,
                        node.outputDefect)) {
                    return fail(
                        ExpressionReferenceBuildStatus::
                            CompactionOutOfRange,
                        "reference tape output cannot be represented by compact shadows");
                }
                if (traced.flags &
                    (OracleTraceHasCompanion |
                     OracleTraceHasDenominator |
                     OracleTraceHasLogarithmBase)) {
                    if (!compact(
                            traced.auxiliary,
                            node.auxiliary,
                            node.auxiliaryDefect)) {
                        return fail(
                            ExpressionReferenceBuildStatus::
                                CompactionOutOfRange,
                            "reference tape auxiliary cannot be represented by compact shadows");
                    }
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
            }
            result.samples.push_back(sample);

            if (!defined) {
                result.undefined = true;
                result.undefinedIteration = iteration + 1;
                result.error = oracleError;
                break;
            }
            context.z.set(next);
            if (outsideBailout(
                    context.z, request.bailout, magnitude)) {
                result.escaped = true;
                result.escapeIteration = iteration + 1;
                break;
            }
        }
    }
    return succeed();
}

bool buildExpressionReferenceOrbit(
        const ExpressionReferenceBuildRequest& request,
        ExpressionReferenceOrbitResult& result) {
    try {
        return buildExpressionReferenceOrbitImpl(
            request, result);
    } catch (const std::bad_alloc&) {
        std::vector<ExpressionReferenceSample>().swap(
            result.samples);
        std::vector<ExpressionReferenceTapeNode>().swap(
            result.tape);
        result.status =
            ExpressionReferenceBuildStatus::ResourceLimit;
        result.error.clear();
        result.valid = false;
        result.sampleCount = 0;
        result.memoryBytes = retainedBytes(result);
        return false;
    } catch (const std::length_error&) {
        std::vector<ExpressionReferenceSample>().swap(
            result.samples);
        std::vector<ExpressionReferenceTapeNode>().swap(
            result.tape);
        result.status =
            ExpressionReferenceBuildStatus::ResourceLimit;
        result.error.clear();
        result.valid = false;
        result.sampleCount = 0;
        result.memoryBytes = retainedBytes(result);
        return false;
    }
}

} // namespace formula
