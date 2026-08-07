#ifndef MANDEL_FORMULA_SCALED_RESIDUAL_H
#define MANDEL_FORMULA_SCALED_RESIDUAL_H

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "formula_reference_orbit.h"

namespace formula {

// A finite nonzero value is mantissa*2^exponent with
// 0.5 <= abs(mantissa) < 1. Unlike a reference shadow, this is an arithmetic
// value: every operation normalizes its result and reports exponent failure.
struct ScaledRealValue {
    double mantissa = 0.0;
    int64_t exponent = 0;

    bool isZero() const { return mantissa == 0.0; }
    bool isFinite() const { return std::isfinite(mantissa); }
    bool isNormalized() const {
        return isZero()
            ? exponent == 0
            : isFinite() &&
              std::abs(mantissa) >= 0.5 &&
              std::abs(mantissa) < 1.0;
    }
};

struct ScaledComplexValue {
    ScaledRealValue re;
    ScaledRealValue im;

    bool isFinite() const {
        return re.isFinite() && im.isFinite();
    }
    bool isNormalized() const {
        return re.isNormalized() && im.isNormalized();
    }
    bool isZero() const {
        return re.isZero() && im.isZero();
    }
};

enum class ScaledArithmeticStatus : uint8_t {
    Success,
    Nonfinite,
    Singular,
    ExponentRange
};

ScaledArithmeticStatus makeScaledRealValue(
    double value, ScaledRealValue& output);
ScaledArithmeticStatus makeScaledComplexValue(
    Complex value, ScaledComplexValue& output);
ScaledArithmeticStatus makeScaledRealValue(
    mpfr_srcptr value, ScaledRealValue& output);
ScaledArithmeticStatus makeScaledComplexValue(
    const MpfrComplex& value, ScaledComplexValue& output);
ScaledArithmeticStatus makeScaledRealValue(
    const ScaledRealShadow& shadow, ScaledRealValue& output);
ScaledArithmeticStatus makeScaledComplexValue(
    const ScaledComplexShadow& shadow, ScaledComplexValue& output);
ScaledArithmeticStatus makeScaledComplexValue(
    const ScaledComplexShadow& primary,
    const ScaledComplexShadow& defect,
    ScaledComplexValue& output);

bool setMpfrFromScaledValue(
    mpfr_ptr output, const ScaledRealValue& value,
    mpfr_rnd_t rounding = MPFR_RNDN);
bool setMpfrFromScaledValue(
    MpfrComplex& output, const ScaledComplexValue& value,
    mpfr_rnd_t rounding = MPFR_RNDN);
bool scaledValueToDouble(
    const ScaledRealValue& value, double& output);
bool scaledValueToDouble(
    const ScaledComplexValue& value, Complex& output);

ScaledArithmeticStatus scaledNegate(
    const ScaledRealValue& value, ScaledRealValue& output);
ScaledArithmeticStatus scaledNegate(
    const ScaledComplexValue& value, ScaledComplexValue& output);
ScaledArithmeticStatus scaledAdd(
    const ScaledRealValue& left, const ScaledRealValue& right,
    ScaledRealValue& output);
ScaledArithmeticStatus scaledSubtract(
    const ScaledRealValue& left, const ScaledRealValue& right,
    ScaledRealValue& output);
ScaledArithmeticStatus scaledMultiply(
    const ScaledRealValue& left, const ScaledRealValue& right,
    ScaledRealValue& output);
ScaledArithmeticStatus scaledDivideByDouble(
    const ScaledRealValue& value, double divisor,
    ScaledRealValue& output);
ScaledArithmeticStatus scaledAdd(
    const ScaledComplexValue& left,
    const ScaledComplexValue& right,
    ScaledComplexValue& output);
ScaledArithmeticStatus scaledSubtract(
    const ScaledComplexValue& left,
    const ScaledComplexValue& right,
    ScaledComplexValue& output);
ScaledArithmeticStatus scaledMultiply(
    const ScaledComplexValue& left,
    const ScaledComplexValue& right,
    ScaledComplexValue& output);
ScaledArithmeticStatus scaledMultiplyByDouble(
    const ScaledComplexValue& value, double multiplier,
    ScaledComplexValue& output);
ScaledArithmeticStatus scaledDivideByDouble(
    const ScaledComplexValue& value, double divisor,
    ScaledComplexValue& output);
ScaledArithmeticStatus scaledNormSquared(
    const ScaledComplexValue& value, ScaledRealValue& output);

// Inputs are finite nonnegative values. The return value has strcmp-like
// ordering. Zero is supported.
int compareScaledNonnegative(
    const ScaledRealValue& left,
    const ScaledRealValue& right);

enum class ExpressionScaledResidualStatus : uint8_t {
    Success,
    BranchUncertain,
    Singular,
    Unsupported,
    Nonfinite,
    ExponentRange,
    InvalidTape,
    InvalidInput
};

const char* expressionScaledResidualStatusName(
    ExpressionScaledResidualStatus status);

struct ExpressionScaledResidualInput {
    // Every delta is relative to the exact two-term reference value
    // primary+defect for the selected sample and fixed context.
    ScaledComplexValue z;
    ScaledComplexValue c;
    ScaledComplexValue z0;
    std::array<ScaledComplexValue, 8> parameters{};
    int iteration = 0;
};

struct ExpressionScaledResidualResult {
    ExpressionScaledResidualStatus status =
        ExpressionScaledResidualStatus::InvalidTape;
    ScaledComplexValue residual;
    // Maximum coefficient-scaled first omitted local-series term. This is
    // diagnostic only, not an interval bound or branch certificate.
    ScaledRealValue remainderEstimate;
    bool uncertified = false;
    size_t operationCount = 0;
    size_t seriesOperationCount = 0;
};

class ExpressionScaledResidualEvaluator {
public:
    // The program and reference must remain alive and unchanged. Evaluation
    // reuses internal node storage and is intentionally not thread-safe.
    bool prepare(
        const ExpressionProgram& program,
        const ExpressionReferenceOrbitResult& reference);
    void reset();

    bool ready() const { return _ready; }
    ExpressionScaledResidualStatus preparationStatus() const {
        return _preparationStatus;
    }
    const std::string& error() const { return _error; }

    ExpressionScaledResidualResult evaluate(
        size_t sampleIndex,
        const ExpressionScaledResidualInput& input);

private:
    struct NodeState {
        ScaledComplexValue residual;
    };

    const ExpressionProgram* _program = nullptr;
    const ExpressionReferenceOrbitResult* _reference = nullptr;
    std::vector<NodeState> _states;
    std::vector<uint16_t> _layoutStack;
    ExpressionScaledResidualStatus _preparationStatus =
        ExpressionScaledResidualStatus::InvalidTape;
    std::string _error;
    bool _ready = false;
};

// The evaluator result is relative to sample.next+sample.rootDefect. Relative
// to sample.next primary, the state is the exact unevaluated two-term sum
// rootDefect+residual. Keeping the terms separate prevents a large compact
// defect from rounding away an e500 residual before the next reset.
struct ExpressionScaledPrimaryRelativeState {
    ScaledComplexValue referenceDefect;
    ScaledComplexValue residual;
};

ScaledArithmeticStatus makeExpressionResidualNextPrimaryState(
    const ExpressionReferenceSample& sample,
    const ScaledComplexValue& exactResidual,
    ExpressionScaledPrimaryRelativeState& primaryRelative);
ScaledArithmeticStatus resetExpressionResidualToExactSample(
    const ExpressionReferenceSample& sample,
    const ExpressionScaledPrimaryRelativeState& primaryRelative,
    ScaledComplexValue& exactResidual);

} // namespace formula

#endif
