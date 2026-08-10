#ifndef MANDEL_FORMULA_DEEP_RENDERER_H
#define MANDEL_FORMULA_DEEP_RENDERER_H

#include <gmp.h>
#include <mpfr.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>

#include "formula_reference_orbit.h"
#include "formula_scaled_residual.h"
#include "formula_taylor_jet.h"

namespace formula {

constexpr float ExpressionDeepEmptyPixel = -10.0f;
constexpr float ExpressionDeepInteriorPixel = -2.0f;

struct ExpressionDeepExactReal {
    mpf_srcptr mpf = nullptr;
    std::string decimal;

    bool usesMpf() const { return mpf != nullptr; }
    bool usesDecimal() const { return !decimal.empty(); }
};

struct ExpressionDeepThreadPolicy {
    int tileWidth = 16;
    int tileHeight = 8;
    // Zero selects the OpenMP runtime default.
    int threads = 0;
};

struct ExpressionDeepMemoryPolicy {
    size_t memoryLimitBytes = size_t{ 1 } << 30;
    mpfr_prec_t fallbackGuardBits = 128;
};

struct ExpressionDeepTaylorPolicy {
    bool enableTaylor = true;
    int minimumLanding = 8;
    int order = 12;
    int minimumOrder = 8;
    int maximumOrder = 20;
    int maximumBivariateOrder =
        ExpressionTaylorMaximumBivariateOrder;
    int maximumCompositionOrder = 24;
    int maximumCandidateIteration = 0;
    double accuracyBudget = 0x1p-40;
    // Reject a built jet when its certified work estimate cannot amortize
    // across the frame.
    bool requirePredictedBenefit = true;
};

enum class ExpressionDeepRenderPhase : uint8_t {
    Reference,
    Fast,
    Fallback,
    Complete
};

enum class ExpressionDeepFallbackReason : uint8_t {
    UncertifiedSeries,
    BranchSensitive,
    UnsupportedOperation,
    Singular,
    Nonfinite,
    ExponentRange,
    InvalidTape,
    ReferenceExhausted,
    CertificationFailure,
    BailoutUncertain,
    ReconstructionFailure,
    Count
};

enum class ExpressionDeepRenderStatus : uint8_t {
    Success,
    Cancelled,
    InvalidRequest,
    ProgramMismatch,
    PrecisionOutOfRange,
    ResourceLimit,
    ReferenceFailure,
    UndefinedPixel,
    InternalError
};

enum class ExpressionDeepVerificationFault : uint8_t {
    None,
    FastWorkerAllocation,
    FastIterationAllocation,
    FallbackWorkerAllocation,
    FallbackIterationAllocation
};

struct ExpressionDeepRenderRequest {
    const ExpressionProgram* canonicalProgram = nullptr;
    const ExpressionProgram* runtimeProgram = nullptr;
    ExpressionReferenceExactInput center;
    ExpressionDeepExactReal scale;
    ExpressionContext fixed;
    FormulaParameter pixelParameter = FormulaParameter::C;
    int width = 0;
    int height = 0;
    int maxIterations = 0;
    double bailout = 4.0;
    float* output = nullptr;
    size_t outputCount = 0;
    ExpressionReferencePrecisionPolicy precision;
    ExpressionDeepThreadPolicy threading;
    ExpressionDeepMemoryPolicy memory;
    ExpressionDeepTaylorPolicy taylor;
    // This is a verification/benchmark switch. Production callers must leave
    // it false because local transcendental series are not interval-certified.
    bool allowUncertifiedForBenchmark = false;
    // Verification-only all-MPFR baseline for the identical formula.
    bool forceMpfrFallbackForVerification = false;
    // Verification-only outward inflation of every certification radius.
    int verificationErrorInflationBits = 0;
    // Verification-only deterministic exception injection. Production callers
    // must leave this at None.
    ExpressionDeepVerificationFault verificationFault =
        ExpressionDeepVerificationFault::None;
    // Cancellation may be polled concurrently by worker threads. Progress
    // callbacks are serialized by the renderer.
    std::function<bool()> shouldCancel;
    std::function<void(
        ExpressionDeepRenderPhase, uint64_t, uint64_t)> progress;
};

struct ExpressionDeepRenderResult {
    ExpressionDeepRenderStatus status =
        ExpressionDeepRenderStatus::InvalidRequest;
    std::string error;
    bool success = false;
    bool cancelled = false;
    double referenceSeconds = 0.0;
    double fastSeconds = 0.0;
    double fallbackSeconds = 0.0;
    uint64_t fastPixelCount = 0;
    uint64_t fallbackPixelCount = 0;
    uint64_t uncertainPixelCount = 0;
    uint64_t undefinedPixelCount = 0;
    std::array<uint64_t,
        static_cast<size_t>(ExpressionDeepFallbackReason::Count)>
        fallbackReasonCounts{};
    uint64_t maxFallbackReasonCount = 0;
    uint64_t fallbackTileCount = 0;
    double maxTileFallbackRate = 0.0;
    size_t referenceBytes = 0;
    size_t rendererBytes = 0;
    mpfr_prec_t selectedPrecision = 0;
    mpfr_prec_t certificationPrecision = 0;
    mpfr_prec_t fallbackPrecision = 0;
    uint64_t totalIterations = 0;
    bool taylorAttempted = false;
    bool taylorAccepted = false;
    int taylorOrder = 0;
    ExpressionTaylorJetLayout taylorLayout =
        ExpressionTaylorJetLayout::ComplexUnivariate;
    size_t taylorMonomialCount = 0;
    uint64_t taylorBivariateConvolutionOperationCount = 0;
    int taylorCoveredIterations = 0;
    int taylorMaximumFunctionSeriesOrder = 0;
    uint64_t taylorFunctionSeriesCount = 0;
    uint64_t taylorFunctionSeriesOperationCount = 0;
    ScaledRealValue taylorMaximumFunctionSeriesTail;
    int taylorMaximumReciprocalOrder = 0;
    uint64_t taylorReciprocalCount = 0;
    uint64_t taylorReciprocalOperationCount = 0;
    ScaledRealValue taylorMinimumDenominatorClearance;
    ScaledRealValue taylorMaximumReciprocalTail;
    bool taylorPoleRejected = false;
    int taylorMaximumBranchSeriesOrder = 0;
    uint64_t taylorBranchCompositionCount = 0;
    uint64_t taylorBranchCompositionOperationCount = 0;
    ScaledRealValue taylorMaximumBranchSeriesTail;
    ScaledRealValue taylorMinimumBranchCutClearance;
    ScaledRealValue taylorMinimumBranchZeroClearance;
    bool taylorBranchRejected = false;
    double taylorBuildSeconds = 0.0;
    // Aggregate worker time; it may exceed fastSeconds under parallelism.
    double taylorEvaluationSeconds = 0.0;
    double taylorResidualSeconds = 0.0;
    uint64_t taylorAcceptedPixelCount = 0;
    uint64_t taylorFallbackPixelCount = 0;
    size_t taylorMemoryBytes = 0;
    ExpressionTaylorJetStatus taylorStatus =
        ExpressionTaylorJetStatus::NoCoverage;
    std::string taylorFailureReason;
    ExpressionScaledResidualCapability capability =
        ExpressionScaledResidualCapability::Unsupported;
};

const char* expressionDeepRenderStatusName(
    ExpressionDeepRenderStatus status);
const char* expressionDeepFallbackReasonName(
    ExpressionDeepFallbackReason reason);

// Renders finite-iteration escape classifications. The arithmetic fast path is
// certified relative to an independently iterated higher-precision MPFR oracle;
// ambiguous bailout intervals fall back per pixel. Exp/sin/cos/sinh/cosh may
// use a certified whole-prefix Taylor jet; if that jet cannot cover the full
// requested horizon profitably, the formula remains all-MPFR. Divide/tan/tanh
// use a distinct certified-meromorphic Taylor tier and reject any full-frame
// denominator neighborhood that is not proven pole-free. Log/log10/sqrt/power
// use a certified principal-branch tier and reject any full-frame input
// neighborhood that is not proven clear of zero and the negative-real cut.
// Conjugate/real/imaginary/norm/complex formulas use a separate certified
// real-bivariate q/conjugate(q) tier. Unsupported mixed transcendental real
// formulas release fast resources and render through MPFR.
// Interior output means only that no escape was observed before maxIterations,
// not mathematical membership.
bool renderExpressionDeepFrame(
    const ExpressionDeepRenderRequest& request,
    ExpressionDeepRenderResult& result);

} // namespace formula

#endif
