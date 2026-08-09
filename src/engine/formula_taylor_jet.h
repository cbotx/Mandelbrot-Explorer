#ifndef MANDEL_FORMULA_TAYLOR_JET_H
#define MANDEL_FORMULA_TAYLOR_JET_H

#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

#include "formula_scaled_residual.h"

namespace formula {

enum class ExpressionTaylorJetStatus : uint8_t {
    Success,
    NoCoverage,
    InvalidRequest,
    UnsupportedProgram,
    InvalidTape,
    ResourceLimit,
    Cancelled,
    ExponentRange,
    Nonfinite,
    AccuracyBudget,
    BailoutUncertain
};

const char* expressionTaylorJetStatusName(
    ExpressionTaylorJetStatus status);

struct ExpressionTaylorJetRequest {
    const ExpressionProgram* program = nullptr;
    const ExpressionReferenceOrbitResult* reference = nullptr;
    FormulaParameter pixelParameter = FormulaParameter::C;
    // The exact parameterization is parameter = reference + D*q. D has no
    // uncertainty and may use the full scaled exponent range.
    ScaledComplexValue parameterScale;
    int minimumOrder = 8;
    int preferredOrder = 12;
    int maximumOrder = 20;
    int minimumLanding = 8;
    // Zero selects the last reference sample that can be used as a landing.
    int maximumCandidateIteration = 0;
    double bailout = 4.0;
    // Uniform maximum-component tail budget on |q| <= 1.
    double accuracyBudget = 0x1p-40;
    size_t memoryLimitBytes = size_t{ 1 } << 30;
    std::function<bool()> shouldCancel;
};

struct ExpressionTaylorJetResult {
    ExpressionTaylorJetStatus status =
        ExpressionTaylorJetStatus::InvalidRequest;
    std::string failureReason;
    bool valid = false;
    bool certified = false;
    FormulaParameter pixelParameter = FormulaParameter::C;
    ScaledComplexValue parameterScale;
    uint64_t programSemanticHash = 0;
    int landingIteration = 0;
    size_t landingSample = 0;
    int order = 0;
    // Index zero is the constant residual caused by compact rebasing.
    std::vector<ScaledComplexValue> coefficients;
    std::vector<ScaledRealValue> coefficientRadii;
    ScaledRealValue remainderRadius;
    std::vector<ScaledRealValue> intermediateEscapeMargins;
    double buildSeconds = 0.0;
    uint64_t operationCount = 0;
    size_t memoryBytes = 0;
};

struct ExpressionTaylorJetEvaluation {
    ExpressionTaylorJetStatus status =
        ExpressionTaylorJetStatus::InvalidRequest;
    bool valid = false;
    ScaledComplexBall residual;
    uint64_t operationCount = 0;
};

class ExpressionTaylorJetBuilder {
public:
    static bool build(
        const ExpressionTaylorJetRequest& request,
        ExpressionTaylorJetResult& result);
};

class ExpressionTaylorJetEvaluator {
public:
    // q may carry a certified construction radius. Evaluation is allocation
    // free and the immutable jet can be shared by renderer threads.
    static bool evaluate(
        const ExpressionTaylorJetResult& jet,
        const ScaledComplexBall& q,
        ExpressionTaylorJetEvaluation& result);
};

// Chooses an exact positive power-of-two D that encloses
// |Re(delta)| + |Im(delta)| plus both component errors.
bool makeExpressionTaylorFrameScale(
    const ScaledRealValue& maximumRealMagnitude,
    const ScaledRealValue& maximumImaginaryMagnitude,
    const ScaledRealValue& maximumComponentError,
    ScaledComplexValue& parameterScale);

// Divides an arbitrary-exponent offset by a power-of-two frame scale without
// losing exponent information.
bool makeExpressionTaylorNormalizedQ(
    const ScaledComplexBall& offset,
    const ScaledComplexValue& parameterScale,
    ScaledComplexBall& q);

// Uses the L1 enclosure, so a true result rigorously implies |q| <= 1.
bool expressionTaylorQInsideUnitDisk(
    const ScaledComplexBall& q);

} // namespace formula

#endif
