#ifndef MANDEL_FORMULA_DEEP_RENDERER_H
#define MANDEL_FORMULA_DEEP_RENDERER_H

#include <gmp.h>
#include <mpfr.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

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
    size_t memoryLimitBytes = size_t{1} << 30;
    mpfr_prec_t fallbackGuardBits = 128;
};

struct ExpressionDeepTaylorPolicy {
    bool enableTaylor = true;
    bool enableTileTaylor = true;
    // Coherently traverses one accepted univariate jet for up to four pixels.
    // Unsupported layouts or failed lane preparation use scalar evaluation.
    bool enableBatchEvaluation = true;
    int minimumLanding = 8;
    int order = 12;
    int minimumOrder = 8;
    int maximumOrder = 20;
    int maximumBivariateOrder = ExpressionTaylorMaximumBivariateOrder;
    int maximumCompositionOrder = 24;
    int maximumCandidateIteration = 0;
    double accuracyBudget = 0x1p-40;
    int maximumDepth = 5;
    int minimumTileWidth = 4;
    int minimumTileHeight = 4;
    size_t maximumJetCount = 256;
    size_t maximumRejectedBeforeFirstAcceptance = 2;
    // Zero uses the remaining renderer memory budget.
    size_t maximumJetMemoryBytes = 0;
    // Reject a built jet when its certified work estimate cannot amortize
    // across the frame or tile.
    bool requirePredictedBenefit = true;
};

struct ExpressionDeepPreflightPolicy {
    bool enable = true;
    size_t maximumSamples = 16;
    size_t minimumSamples = 8;
    uint32_t earlyRejectMinimumFirstUncertainIteration = 8;
};

struct ExpressionDeepTransferPolicy {
    // Opt-in prototype for ExactCenteredArithmetic. It builds one certified
    // terminal segment over the full frame and otherwise routes directly to
    // the existing exact MPFR fallback.
    bool enableCertifiedSegments = false;
    bool requirePredictedBenefit = true;
    size_t minimumPixelCount = 256;
    int minimumIterations = 32;
};

struct ExpressionDeepCenteredPolicy {
    // Default-off exact/entire/real-smooth c- or z0-plane candidate. It is
    // never exposed as a render mode: failed admission or sparse verification
    // continues through exact MPFR.
    bool enableAdaptiveCandidate = false;
    bool requirePredictedBenefit = true;
    size_t minimumPixelCount = 256;
    int minimumIterations = 256;
    size_t maximumReferences = 5;
    size_t preflightSamples = 256;
    double maximumPredictedFallbackRate = 0.125;
};

enum class ExpressionDeepRenderPhase : uint8_t { Reference,
                                                 Preflight,
                                                 Fast,
                                                 Fallback,
                                                 Complete };

enum class ExpressionDeepPreflightDecision : uint8_t { NotRun,
                                                       ContinueCertifiedFast,
                                                       DirectMpfr };

constexpr size_t ExpressionDeepUncertaintyHistogramBins = 16;

enum class ExpressionDeepFallbackReason : uint8_t { UncertifiedSeries,
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
                                                    Count };

enum class ExpressionDeepCenteredFallbackReason : uint8_t { CandidateFailure,
                                                            Nonfinite,
                                                            ReferenceOutputDisagreement,
                                                            StateDisagreement,
                                                            InteriorConservatism,
                                                            DoubleDoubleRejected,
                                                            StructuralInconsistency,
                                                            DenominatorInstability,
                                                            Count };

constexpr uint8_t expressionDeepCenteredFallbackReasonMask(ExpressionDeepCenteredFallbackReason reason) {
    return uint8_t{1} << static_cast<uint8_t>(reason);
}

constexpr uint8_t ExpressionDeepCenteredFallbackMandatoryMask = expressionDeepCenteredFallbackReasonMask(ExpressionDeepCenteredFallbackReason::CandidateFailure) | expressionDeepCenteredFallbackReasonMask(ExpressionDeepCenteredFallbackReason::Nonfinite) | expressionDeepCenteredFallbackReasonMask(ExpressionDeepCenteredFallbackReason::ReferenceOutputDisagreement) | expressionDeepCenteredFallbackReasonMask(ExpressionDeepCenteredFallbackReason::DoubleDoubleRejected) | expressionDeepCenteredFallbackReasonMask(ExpressionDeepCenteredFallbackReason::StructuralInconsistency) | expressionDeepCenteredFallbackReasonMask(ExpressionDeepCenteredFallbackReason::DenominatorInstability);

enum class ExpressionDeepRenderStatus : uint8_t { Success,
                                                  Cancelled,
                                                  InvalidRequest,
                                                  ProgramMismatch,
                                                  PrecisionOutOfRange,
                                                  ResourceLimit,
                                                  ReferenceFailure,
                                                  UndefinedPixel,
                                                  InternalError };

enum class ExpressionDeepVerificationFault : uint8_t { None,
                                                       FastWorkerAllocation,
                                                       FastIterationAllocation,
                                                       FallbackWorkerAllocation,
                                                       FallbackIterationAllocation };

struct ExpressionDeepRenderRequest {
    const ExpressionProgram* canonicalProgram = nullptr;
    const ExpressionProgram* runtimeProgram = nullptr;
    ExpressionReferenceExactInput center;
    ExpressionDeepExactReal scale;
    ExpressionContext fixed;
    FormulaParameter pixelParameter = FormulaParameter::C;
    int width = 0;
    int height = 0;
    int fullHeight = 0;
    int rowBase = 0;
    int maxIterations = 0;
    double bailout = 4.0;
    float* output = nullptr;
    size_t outputCount = 0;
    ExpressionReferencePrecisionPolicy precision;
    ExpressionDeepThreadPolicy threading;
    ExpressionDeepMemoryPolicy memory;
    ExpressionDeepTaylorPolicy taylor;
    ExpressionDeepPreflightPolicy preflight;
    ExpressionDeepTransferPolicy transfer;
    ExpressionDeepCenteredPolicy centered;
    // This is a verification/benchmark switch. Production callers must leave
    // it false because local transcendental series are not interval-certified.
    bool allowUncertifiedForBenchmark = false;
    // Verification-only all-MPFR baseline for the identical formula.
    bool forceMpfrFallbackForVerification = false;
    // Verification-only generic bytecode-oracle baseline.
    bool disableSpecializedPiecewiseMpfrForVerification = false;
    // Verification-only specialized-MPFR baseline without BigFixed.
    bool disablePiecewiseBigFixedForVerification = false;
    // Verification-only outward inflation of every certification radius.
    int verificationErrorInflationBits = 0;
    // Verification-only selector calibration overrides. Production callers
    // must leave these at zero to use the fixed production thresholds.
    double verificationCenteredPrimaryErrorThreshold = 0.0;
    double verificationCenteredStateDisagreementThreshold = 0.0;
    // Verification-only adaptive-fallback precision/sample overrides.
    mpfr_prec_t verificationCenteredFallbackCandidatePrecision = 0;
    size_t verificationCenteredFallbackValidationSamples = 0;
    mpfr_prec_t verificationCenteredMandatoryFallbackCandidatePrecision = 0;
    size_t verificationCenteredMandatoryFallbackValidationSamples = 0;
    // Verification-only deterministic exception injection. Production callers
    // must leave this at None.
    ExpressionDeepVerificationFault verificationFault = ExpressionDeepVerificationFault::None;
    // Cancellation may be polled concurrently by worker threads. Progress
    // callbacks are serialized by the renderer.
    std::function<bool()> shouldCancel;
    std::function<void(ExpressionDeepRenderPhase, uint64_t, uint64_t)> progress;
};

struct ExpressionDeepRenderResult {
    ExpressionDeepRenderStatus status = ExpressionDeepRenderStatus::InvalidRequest;
    std::string error;
    bool success = false;
    bool cancelled = false;
    double referenceSeconds = 0.0;
    double preflightSeconds = 0.0;
    double fastSeconds = 0.0;
    double fallbackSeconds = 0.0;
    bool preflightAttempted = false;
    bool preflightRejectedFast = false;
    ExpressionDeepPreflightDecision preflightDecision = ExpressionDeepPreflightDecision::NotRun;
    uint64_t preflightSampleCount = 0;
    uint64_t preflightFallbackCount = 0;
    uint64_t preflightFastCount = 0;
    uint64_t preflightIterationCount = 0;
    uint64_t preflightOperationCount = 0;
    uint64_t preflightFoldOperationCount = 0;
    uint64_t preflightUncertainFoldCount = 0;
    uint32_t preflightMinimumFirstUncertainIteration = UINT32_MAX;
    uint32_t preflightMaximumFirstUncertainIteration = 0;
    uint64_t preflightAvoidedFastPixelCount = 0;
    uint64_t preflightAvoidedIterationEstimate = 0;
    uint64_t preflightAvoidedOperationEstimate = 0;
    std::array<uint64_t, ExpressionDeepUncertaintyHistogramBins> preflightFirstUncertainHistogram{};
    std::array<uint64_t, static_cast<size_t>(ExpressionDeepFallbackReason::Count)> preflightFallbackReasonCounts{};
    uint64_t fastPixelCount = 0;
    uint64_t fallbackPixelCount = 0;
    uint64_t uncertainPixelCount = 0;
    uint64_t undefinedPixelCount = 0;
    std::array<uint64_t, static_cast<size_t>(ExpressionDeepFallbackReason::Count)> fallbackReasonCounts{};
    uint64_t maxFallbackReasonCount = 0;
    uint64_t fallbackTileCount = 0;
    double maxTileFallbackRate = 0.0;
    size_t referenceBytes = 0;
    size_t rendererBytes = 0;
    mpfr_prec_t selectedPrecision = 0;
    mpfr_prec_t certificationPrecision = 0;
    mpfr_prec_t fallbackPrecision = 0;
    uint64_t totalIterations = 0;
    uint64_t fastIterationCount = 0;
    uint64_t fastOperationCount = 0;
    uint64_t fastSeriesOperationCount = 0;
    uint64_t fastFoldOperationCount = 0;
    uint64_t fastUncertainFoldCount = 0;
    std::array<uint64_t, ExpressionDeepUncertaintyHistogramBins> fallbackFirstUncertainHistogram{};
    bool transferAttempted = false;
    bool transferAccepted = false;
    uint64_t transferAcceptedSegmentCount = 0;
    int transferCoveredIterations = 0;
    uint64_t transferSkippedIterationCount = 0;
    double transferBuildSeconds = 0.0;
    double transferApplySeconds = 0.0;
    size_t transferMemoryBytes = 0;
    ScaledRealValue transferFinalRadius;
    bool centeredAttempted = false;
    bool centeredAccepted = false;
    bool centeredPreflightRejected = false;
    uint64_t centeredReferenceCount = 0;
    uint64_t centeredReferenceAttemptCount = 0;
    uint64_t centeredPreflightSampleCount = 0;
    uint64_t centeredPreflightPrimaryRiskFlagCount = 0;
    uint64_t centeredPreflightSecondaryEvaluationCount = 0;
    uint64_t centeredPreflightInitialFlagCount = 0;
    uint64_t centeredPreflightAdditionalReferenceEvaluationCount = 0;
    uint64_t centeredPreflightFlagCount = 0;
    uint64_t centeredPrimaryRiskFlagCount = 0;
    uint64_t centeredSecondaryEvaluationCount = 0;
    uint64_t centeredSelectorFlagCount = 0;
    uint64_t centeredAdditionalReferenceEvaluationCount = 0;
    uint64_t centeredHierarchyEvaluationCount = 0;
    uint64_t centeredFinalFallbackFlagCount = 0;
    uint64_t centeredDoubleDoubleVerifiedPixelCount = 0;
    uint64_t centeredDoubleDoubleAgreementCount = 0;
    uint64_t centeredDoubleDoubleRejectedCount = 0;
    uint64_t centeredDoubleDoubleIterationCount = 0;
    uint64_t centeredVectorStepCount = 0;
    uint64_t centeredVectorActiveLaneCount = 0;
    double centeredPreflightSeconds = 0.0;
    double centeredPrimarySeconds = 0.0;
    double centeredSecondarySeconds = 0.0;
    double centeredCandidateSeconds = 0.0;
    double centeredSelectorSeconds = 0.0;
    double centeredAdditionalReferenceSeconds = 0.0;
    double centeredFinalSelectorSeconds = 0.0;
    double centeredDoubleDoubleSeconds = 0.0;
    mpfr_prec_t centeredFallbackCandidatePrecision = 0;
    mpfr_prec_t centeredFallbackFullPrecision = 0;
    bool centeredFallbackValidationIsEmpirical = false;
    uint64_t centeredFallbackValidationSampleCount = 0;
    uint64_t centeredFallbackValidationMismatchCount = 0;
    bool centeredFallbackUpgraded = false;
    uint64_t centeredFallbackMandatoryFullPixelCount = 0;
    mpfr_prec_t centeredFallbackMandatoryCandidatePrecision = 0;
    uint64_t centeredFallbackMandatoryCandidatePixelCount = 0;
    uint64_t centeredFallbackMandatoryCandidateIterationCount = 0;
    double centeredFallbackMandatoryCandidateSeconds = 0.0;
    uint64_t centeredFallbackMandatoryFullPrecisionPixelCount = 0;
    uint64_t centeredFallbackMandatoryFullPrecisionIterationCount = 0;
    double centeredFallbackMandatoryFullPrecisionSeconds = 0.0;
    uint64_t centeredFallbackMandatoryValidationSampleCount = 0;
    uint64_t centeredFallbackMandatoryValidationMismatchCount = 0;
    bool centeredFallbackMandatoryUpgraded = false;
    uint64_t centeredFallbackMandatoryUpgradedPixelCount = 0;
    uint64_t centeredFallbackLowEligiblePixelCount = 0;
    uint64_t centeredFallbackUpgradedPixelCount = 0;
    std::array<uint64_t, static_cast<size_t>(ExpressionDeepCenteredFallbackReason::Count)> centeredFallbackPrecisionReasonCounts{};
    std::vector<uint8_t> centeredFallbackPrecisionReasonMask;
    uint64_t centeredFallbackLowPrecisionPixelCount = 0;
    uint64_t centeredFallbackLowPrecisionIterationCount = 0;
    double centeredFallbackLowPrecisionSeconds = 0.0;
    uint64_t centeredFallbackFullPrecisionPixelCount = 0;
    uint64_t centeredFallbackFullPrecisionIterationCount = 0;
    double centeredFallbackFullPrecisionSeconds = 0.0;
    uint64_t centeredFallbackMandatoryFullIterationCount = 0;
    double centeredFallbackMandatoryFullSeconds = 0.0;
    bool usedSpecializedPiecewiseMpfr = false;
    uint64_t specializedPiecewiseMpfrPixelCount = 0;
    uint64_t specializedPiecewiseMpfrIterationCount = 0;
    uint64_t specializedPiecewiseMpfrPeriodicPixelCount = 0;
    uint64_t genericMpfrPeriodicPixelCount = 0;
    bool usedPiecewiseBigFixed = false;
    uint64_t piecewiseBigFixedPixelCount = 0;
    uint64_t piecewiseBigFixedIterationCount = 0;
    bool taylorAttempted = false;
    bool taylorAccepted = false;
    int taylorOrder = 0;
    ExpressionTaylorJetLayout taylorLayout = ExpressionTaylorJetLayout::ComplexUnivariate;
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
    uint64_t taylorArgCompositionCount = 0;
    std::string taylorArgRejectionReason;
    uint64_t taylorPolarCompositionCount = 0;
    ScaledRealValue taylorMinimumPolarRadiusClearance;
    bool taylorPolarRejected = false;
    std::string taylorPolarRejectionReason;
    uint64_t taylorAbsBranchCount = 0;
    uint64_t taylorAbsPositiveCellCount = 0;
    uint64_t taylorAbsNegativeCellCount = 0;
    ScaledRealValue taylorMinimumFoldClearance;
    bool taylorFoldRejected = false;
    int taylorFoldRejectionIteration = -1;
    std::string taylorFoldRejectionReason;
    double taylorBuildSeconds = 0.0;
    // Planner/setup time excluding time reported by Taylor builders.
    double taylorPlanningSeconds = 0.0;
    // Aggregate worker time; it may exceed fastSeconds under parallelism.
    double taylorEvaluationSeconds = 0.0;
    double taylorResidualSeconds = 0.0;
    uint64_t taylorAcceptedPixelCount = 0;
    uint64_t taylorFallbackPixelCount = 0;
    uint64_t taylorAttemptedJetCount = 0;
    uint64_t taylorAcceptedJetCount = 0;
    uint64_t taylorRejectedJetCount = 0;
    uint64_t taylorTileSplitCount = 0;
    int taylorMaximumTileDepth = 0;
    uint64_t taylorAcceptedPixelCoverage = 0;
    double taylorWeightedLanding = 0.0;
    uint64_t taylorFoldRejectedJetCount = 0;
    uint64_t taylorCutRejectedJetCount = 0;
    uint64_t taylorPoleRejectedJetCount = 0;
    uint64_t taylorTileMapHash = 0;
    size_t taylorMemoryBytes = 0;
    size_t taylorRetainedBytes = 0;
    ExpressionTaylorJetStatus taylorStatus = ExpressionTaylorJetStatus::NoCoverage;
    std::string taylorFailureReason;
    ExpressionScaledResidualCapability capability = ExpressionScaledResidualCapability::Unsupported;
};

const char* expressionDeepRenderStatusName(ExpressionDeepRenderStatus status);
const char* expressionDeepFallbackReasonName(ExpressionDeepFallbackReason reason);
const char* expressionDeepCenteredFallbackReasonName(ExpressionDeepCenteredFallbackReason reason);
const char* expressionDeepPreflightDecisionName(ExpressionDeepPreflightDecision decision);

// Renders finite-iteration escape classifications. The arithmetic fast path is
// certified relative to an independently iterated higher-precision MPFR oracle;
// ambiguous bailout intervals fall back per pixel. Taylor first attempts one
// full-frame jet. When enabled, a rejected fold/cut/pole neighborhood is then
// subdivided deterministically into certified tile-local P+D*q jets, all
// relative to the same higher-precision reference tape. Accepted immutable
// jets land into the existing certified residual tail; uncovered leaves retain
// the per-step scaled path when supported and otherwise use MPFR. Existing
// piecewise formulas first run a bounded deterministic certified sample; when
// every representative pixel becomes uncertain, it can only reject fast work
// and routes the whole frame to MPFR. Recognized componentwise-abs quadratics
// use an operation-equivalent reusable MPFR kernel after that rejection.
// Existing entire, meromorphic, principal-branch, Arg, Polar, real-bivariate,
// bailout,
// first-escape, and profitability gates remain in force for every tile.
// Interior output means only that no escape was observed before maxIterations,
// not mathematical membership.
bool renderExpressionDeepFrame(const ExpressionDeepRenderRequest& request, ExpressionDeepRenderResult& result);

} // namespace formula

#endif
