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
        int order,
        const ScaledRealValue& remainder,
        const TaylorThreshold& threshold,
        ScaledRealValue& margin,
        ScaledArithmeticStatus& arithmetic) {
    if (!coefficients || order < 0 ||
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
    for (int coefficient = 1;
         coefficient <= order; ++coefficient) {
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
    }
    return "invalid-request";
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

bool ExpressionTaylorJetBuilder::build(
        const ExpressionTaylorJetRequest& request,
        ExpressionTaylorJetResult& result) {
    const Clock::time_point start = Clock::now();
    result = {};
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
        request.program->scaledResidualCapability() !=
            ExpressionScaledResidualCapability::
                ExactCenteredArithmetic ||
        (request.pixelParameter != FormulaParameter::C &&
         request.pixelParameter !=
             FormulaParameter::InitialZ) ||
        request.minimumOrder < 8 ||
        request.maximumOrder > 20 ||
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
        case Op::Square:
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
                  request.reference->samples.size()) - 1,
        static_cast<int>(
            request.reference->samples.size()) - 1);
    if (maximumLanding < request.minimumLanding)
        return finish(
            ExpressionTaylorJetStatus::NoCoverage,
            "Taylor reference is shorter than minimum landing",
            false);

    Candidate best;
    bool haveBest = false;
    uint64_t totalOperations = 0;
    size_t peakMemory = 0;
    ExpressionTaylorJetStatus lastStatus =
        ExpressionTaylorJetStatus::NoCoverage;
    std::string lastReason = "Taylor prefix has no coverage";

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
                candidateBytes, nodeCount * stride,
                sizeof(ScaledComplexBall)) ||
            !checkedAddSize(
                candidateBytes, nodeCount,
                sizeof(ScaledRealValue)) ||
            !checkedAddSize(
                candidateBytes, stride * 4,
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
            auto addNodeRoundoff = [&](
                    const ExpressionReferenceTapeNode& node,
                    size_t sampleOffset,
                    ScaledRealValue& remainder) {
                if (request.reference->
                        certificationPrecision <= 16)
                    return ScaledArithmeticStatus::Success;
                if (node.operation !=
                        ExpressionOracleOperation::Add &&
                    node.operation !=
                        ExpressionOracleOperation::Subtract &&
                    node.operation !=
                        ExpressionOracleOperation::Multiply &&
                    node.operation !=
                        ExpressionOracleOperation::Square)
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
                ScaledArithmeticStatus status =
                    fullBound(node.leftNode, left);
                if (status !=
                        ScaledArithmeticStatus::Success)
                    return status;
                ScaledRealValue operationMagnitude;
                if (node.operation ==
                        ExpressionOracleOperation::Square) {
                    status = multiplyBound(
                        left, left, operationMagnitude);
                } else {
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
                    default:
                        arithmetic =
                            ScaledArithmeticStatus::Singular;
                        break;
                    }
                    if (arithmetic !=
                            ScaledArithmeticStatus::Success)
                        break;
                    arithmetic = addRemainder(
                        remainders[local],
                        node.outputError);
                    if (arithmetic ==
                            ScaledArithmeticStatus::Success)
                        arithmetic = addNodeRoundoff(
                            node, sampleOffset,
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

    result.order = best.order;
    result.landingIteration = best.landing;
    result.landingSample =
        static_cast<size_t>(best.landing);
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
    if (!jet.valid || !jet.certified ||
        jet.order < 1 ||
        jet.coefficients.size() !=
            static_cast<size_t>(jet.order) + 1 ||
        jet.coefficientRadii.size() !=
            jet.coefficients.size() ||
        !validRadius(jet.remainderRadius) ||
        !expressionTaylorQInsideUnitDisk(q)) {
        result.status =
            ExpressionTaylorJetStatus::InvalidRequest;
        return false;
    }
    ScaledComplexBall value;
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
