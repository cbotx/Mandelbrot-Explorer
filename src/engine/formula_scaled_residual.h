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

struct ScaledComplexValue {
    ScaledRealValue re;
    ScaledRealValue im;

    bool isFinite() const { return re.isFinite() && im.isFinite(); }
    bool isNormalized() const { return re.isNormalized() && im.isNormalized(); }
    bool isZero() const { return re.isZero() && im.isZero(); }
};

struct ScaledComplexBall {
    ScaledComplexValue value;
    // Maximum absolute component error.
    ScaledRealValue radius;
};

struct ScaledRealBall {
    ScaledRealValue value;
    ScaledRealValue radius;
};

enum class ScaledArithmeticStatus : uint8_t { Success,
                                              Nonfinite,
                                              Singular,
                                              ExponentRange };

ScaledArithmeticStatus makeScaledRealValue(double value, ScaledRealValue& output);
ScaledArithmeticStatus makeScaledComplexValue(Complex value, ScaledComplexValue& output);
ScaledArithmeticStatus makeScaledRealValue(mpfr_srcptr value, ScaledRealValue& output);
ScaledArithmeticStatus makeScaledNonnegativeUpward(mpfr_srcptr value, ScaledRealValue& output);
ScaledArithmeticStatus makeScaledComplexValue(const MpfrComplex& value, ScaledComplexValue& output);
ScaledArithmeticStatus makeScaledRealValue(const ScaledRealShadow& shadow, ScaledRealValue& output);
ScaledArithmeticStatus makeScaledComplexValue(const ScaledComplexShadow& shadow, ScaledComplexValue& output);
ScaledArithmeticStatus makeScaledComplexValue(const ScaledComplexShadow& primary, const ScaledComplexShadow& defect, ScaledComplexValue& output);

bool setMpfrFromScaledValue(mpfr_ptr output, const ScaledRealValue& value, mpfr_rnd_t rounding = MPFR_RNDN);
bool setMpfrFromScaledValue(MpfrComplex& output, const ScaledComplexValue& value, mpfr_rnd_t rounding = MPFR_RNDN);
bool scaledValueToDouble(const ScaledRealValue& value, double& output);
bool scaledValueToDouble(const ScaledComplexValue& value, Complex& output);

ScaledArithmeticStatus scaledNegate(const ScaledRealValue& value, ScaledRealValue& output);
ScaledArithmeticStatus scaledNegate(const ScaledComplexValue& value, ScaledComplexValue& output);
ScaledArithmeticStatus scaledAdd(const ScaledRealValue& left, const ScaledRealValue& right, ScaledRealValue& output);
ScaledArithmeticStatus scaledSubtract(const ScaledRealValue& left, const ScaledRealValue& right, ScaledRealValue& output);
ScaledArithmeticStatus scaledMultiply(const ScaledRealValue& left, const ScaledRealValue& right, ScaledRealValue& output);
ScaledArithmeticStatus scaledDivideByDouble(const ScaledRealValue& value, double divisor, ScaledRealValue& output);
ScaledArithmeticStatus scaledAdd(const ScaledComplexValue& left, const ScaledComplexValue& right, ScaledComplexValue& output);
ScaledArithmeticStatus scaledSubtract(const ScaledComplexValue& left, const ScaledComplexValue& right, ScaledComplexValue& output);
ScaledArithmeticStatus scaledMultiply(const ScaledComplexValue& left, const ScaledComplexValue& right, ScaledComplexValue& output);
ScaledArithmeticStatus scaledMultiplyByDouble(const ScaledComplexValue& value, double multiplier, ScaledComplexValue& output);
ScaledArithmeticStatus scaledDivideByDouble(const ScaledComplexValue& value, double divisor, ScaledComplexValue& output);
ScaledArithmeticStatus scaledNormSquared(const ScaledComplexValue& value, ScaledRealValue& output);
ScaledArithmeticStatus scaledAddUp(const ScaledRealValue& left, const ScaledRealValue& right, ScaledRealValue& output);
ScaledArithmeticStatus scaledMultiplyUp(const ScaledRealValue& left, const ScaledRealValue& right, ScaledRealValue& output);
ScaledArithmeticStatus certifiedScaledAdd(const ScaledComplexBall& left, const ScaledComplexBall& right, ScaledComplexBall& output);
ScaledArithmeticStatus certifiedScaledSubtract(const ScaledComplexBall& left, const ScaledComplexBall& right, ScaledComplexBall& output);
ScaledArithmeticStatus certifiedScaledMultiply(const ScaledComplexBall& left, const ScaledComplexBall& right, ScaledComplexBall& output);
ScaledArithmeticStatus certifiedScaledNormSquared(const ScaledComplexValue& value, const ScaledRealValue& radius, ScaledRealValue& midpoint, ScaledRealValue& error);
// Certification is relative to finite MPFR arithmetic. These guards reject
// values close enough to the runtime exponent limits that MPFR overflow or
// underflow could differ from int64-exponent scaled arithmetic.
ScaledArithmeticStatus certifyScaledMpfrExponentRange(const ScaledRealValue& value);
ScaledArithmeticStatus certifyScaledMpfrExponentRange(const ScaledComplexValue& value);
ScaledArithmeticStatus certifyScaledMpfrExponentRange(const ScaledComplexBall& value);

// Inputs are finite nonnegative values. The return value has strcmp-like
// ordering. Zero is supported.
int compareScaledNonnegative(const ScaledRealValue& left, const ScaledRealValue& right);

enum class ExpressionScaledResidualStatus : uint8_t { Success,
                                                      BranchUncertain,
                                                      Singular,
                                                      Unsupported,
                                                      Nonfinite,
                                                      ExponentRange,
                                                      InvalidTape,
                                                      InvalidInput };

const char* expressionScaledResidualStatusName(ExpressionScaledResidualStatus status);

// Computes |reference+residual|-|reference| without subtracting the two
// magnitudes. Both arguments are correlated real balls: residual bounds the
// endpoint displacement from the same reference value. A ball that contains
// zero cannot certify an absolute-value branch.
ExpressionScaledResidualStatus certifiedScaledDiffAbsReal(const ScaledRealBall& reference, const ScaledRealBall& residual, ScaledRealBall& output);

struct ExpressionScaledResidualInput {
    // Every delta is relative to the fast reconstruction of primary+defect;
    // the corresponding input error includes compact/reference discrepancy.
    ScaledComplexValue z;
    ScaledRealValue zError;
    ScaledComplexValue c;
    ScaledRealValue cError;
    ScaledComplexValue z0;
    ScaledRealValue z0Error;
    std::array<ScaledComplexValue, 8> parameters{};
    std::array<ScaledRealValue, 8> parameterErrors{};
    int iteration = 0;
};

struct ExpressionScaledResidualResult {
    ExpressionScaledResidualStatus status = ExpressionScaledResidualStatus::InvalidTape;
    ScaledComplexValue residual;
    // Certified maximum-component radius for arithmetic-only V1 programs,
    // relative to the reference's higher-precision finite MPFR oracle.
    ScaledRealValue radius;
    // Maximum coefficient-scaled first omitted local-series term. This is
    // diagnostic only, not an interval bound or branch certificate.
    ScaledRealValue remainderEstimate;
    bool uncertified = false;
    bool certified = false;
    size_t operationCount = 0;
    size_t seriesOperationCount = 0;
    size_t foldOperationCount = 0;
    size_t uncertainFoldCount = 0;
};

class ExpressionScaledResidualEvaluator {
  public:
    // The program and reference must remain alive and unchanged. Evaluation
    // reuses internal node storage and is intentionally not thread-safe.
    bool prepare(const ExpressionProgram& program, const ExpressionReferenceOrbitResult& reference);
    static bool estimateWorkspaceBytes(const ExpressionProgram& program, const ExpressionReferenceOrbitResult& reference, size_t& bytes);
    static bool estimateReferenceCacheBytes(size_t tapeNodeCount, size_t& bytes);
    void reset();

    bool ready() const { return _ready; }
    size_t workspaceBytes() const { return _states.capacity() * sizeof(NodeState) + _layoutStack.capacity() * sizeof(uint16_t) + _sampleExponentRangeUnsafe.capacity() * sizeof(uint8_t) + _sampleUndefinedStatus.capacity() * sizeof(uint8_t) + _nodeExponentRangeUnsafe.capacity() * sizeof(uint8_t) + (_nodeBases.capacity() + _nodeAuxiliaries.capacity()) * sizeof(ScaledComplexValue) + _nodeBaseMagnitudes.capacity() * sizeof(ScaledRealValue); }
    ExpressionScaledResidualStatus preparationStatus() const { return _preparationStatus; }
    const std::string& error() const { return _error; }

    ExpressionScaledResidualResult evaluate(size_t sampleIndex, const ExpressionScaledResidualInput& input);

  private:
    struct NodeState {
        ScaledComplexValue residual;
        ScaledRealValue radius;
        bool realValued = false;
        bool zeroValued = false;
    };

    const ExpressionProgram* _program = nullptr;
    const ExpressionReferenceOrbitResult* _reference = nullptr;
    std::vector<NodeState> _states;
    std::vector<uint16_t> _layoutStack;
    std::vector<uint8_t> _sampleExponentRangeUnsafe;
    std::vector<uint8_t> _sampleUndefinedStatus;
    std::vector<uint8_t> _nodeExponentRangeUnsafe;
    std::vector<ScaledComplexValue> _nodeBases;
    std::vector<ScaledComplexValue> _nodeAuxiliaries;
    std::vector<ScaledRealValue> _nodeBaseMagnitudes;
    ExpressionScaledResidualStatus _preparationStatus = ExpressionScaledResidualStatus::InvalidTape;
    std::string _error;
    int64_t _mpfrSafeMinimum = 0;
    int64_t _mpfrSafeMaximum = 0;
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

ScaledArithmeticStatus makeExpressionResidualNextPrimaryState(const ExpressionReferenceSample& sample, const ScaledComplexValue& exactResidual, ExpressionScaledPrimaryRelativeState& primaryRelative);
ScaledArithmeticStatus resetExpressionResidualToExactSample(const ExpressionReferenceSample& sample, const ExpressionScaledPrimaryRelativeState& primaryRelative, ScaledComplexValue& exactResidual);

} // namespace formula

#endif
