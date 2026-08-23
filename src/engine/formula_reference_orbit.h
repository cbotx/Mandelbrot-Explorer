#ifndef MANDEL_FORMULA_REFERENCE_ORBIT_H
#define MANDEL_FORMULA_REFERENCE_ORBIT_H

#include <gmp.h>
#include <mpfr.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "formula_expression.h"
#include "formula_expression_oracle.h"

namespace formula {

// A finite value is mantissa*2^exponent with a nonzero binary64 mantissa.
// Zero/NaN/infinity are encoded by the mantissa itself. Finite shadows avoid
// binary64's exponent range, but reconstruction is intentionally limited to
// the current MPFR runtime range [mpfr_get_emin(), mpfr_get_emax()]. The
// int64 field is a stable storage type, not a promise that every int64
// exponent is reconstructible by the linked MPFR build.
struct ScaledRealShadow {
    double mantissa = 0.0;
    int64_t exponent = 0;

    bool isZero() const { return mantissa == 0.0; }
    bool isNan() const { return std::isnan(mantissa); }
    bool isInfinity() const { return std::isinf(mantissa); }
    bool isFinite() const { return std::isfinite(mantissa) && mantissa != 0.0; }
};

struct ScaledComplexShadow {
    ScaledRealShadow re;
    ScaledRealShadow im;
};

struct ScaledComplexExpansionTail {
    // Third and fourth terms. The existing defect is the second term.
    ScaledComplexShadow residual2;
    ScaledComplexShadow residual3;
};

// A finite arithmetic value is mantissa*2^exponent with
// 0.5 <= abs(mantissa) < 1. Error radii use the same representation and are
// always finite and nonnegative.
struct ScaledRealValue {
    double mantissa = 0.0;
    int64_t exponent = 0;

    bool isZero() const { return mantissa == 0.0; }
    bool isFinite() const { return std::isfinite(mantissa); }
    bool isNormalized() const { return isZero() ? exponent == 0 : isFinite() && std::abs(mantissa) >= 0.5 && std::abs(mantissa) < 1.0; }
};

// Conversion is fallible: every finite result is guaranteed to be
// reconstructible in the linked MPFR runtime exponent range.
bool makeScaledRealShadow(mpfr_srcptr value, ScaledRealShadow& output);
bool makeScaledComplexShadow(const MpfrComplex& value, ScaledComplexShadow& output);
bool setMpfrFromScaledShadow(mpfr_ptr output, const ScaledRealShadow& shadow, mpfr_rnd_t rounding = MPFR_RNDN);
bool setMpfrFromScaledShadow(MpfrComplex& output, const ScaledComplexShadow& shadow, mpfr_rnd_t rounding = MPFR_RNDN);
bool reconstructMpfrFromShadows(MpfrComplex& output, const ScaledComplexShadow& primary, const ScaledComplexShadow& defect, mpfr_rnd_t rounding = MPFR_RNDN);
bool reconstructMpfrFromExpansion(MpfrComplex& output, const ScaledComplexShadow& primary, const ScaledComplexShadow& defect, const ScaledComplexExpansionTail* tail, mpfr_rnd_t rounding = MPFR_RNDN);

struct ExpressionReferenceExactInput {
    // Select exactly one representation. Both mpf pointers must be non-null,
    // or realDecimal must be non-empty (an empty imaginaryDecimal means zero).
    // Any nonempty decimal field selects the decimal representation, so decimal
    // text can never be silently ignored when MPF pointers are also present.
    mpf_srcptr realMpf = nullptr;
    mpf_srcptr imaginaryMpf = nullptr;
    std::string realDecimal;
    std::string imaginaryDecimal;

    bool usesMpf() const { return realMpf != nullptr || imaginaryMpf != nullptr; }
    bool usesDecimal() const { return !realDecimal.empty() || !imaginaryDecimal.empty(); }
};

struct ExpressionReferencePrecisionPolicy {
    // Hard builder cap, independent of a caller's narrower maximumBits policy.
    static constexpr mpfr_prec_t ApplicationMaximumBits = 1 << 16;

    mpfr_prec_t requestedBits = 0;
    mpfr_prec_t viewBits = 0;
    mpfr_prec_t guardBits = 64;
    mpfr_prec_t minimumBits = 128;
    mpfr_prec_t maximumBits = 1 << 20;
};

enum class ExpressionReferenceBuildStatus : uint8_t { Success,
                                                      InvalidRequest,
                                                      ProgramMismatch,
                                                      UnsupportedProgram,
                                                      PrecisionOutOfRange,
                                                      InputParseError,
                                                      CompactionOutOfRange,
                                                      ResourceLimit };

enum class ExpressionReferenceCompaction : uint8_t {
    TwoTerm,
    FourTermCertifiedTransfer
};

struct ExpressionReferenceBuildRequest {
    // canonicalProgram is the user formula before specialization.
    // runtimeProgram is the fixed-parameter specialization used by the
    // renderer. If both are supplied, the builder independently specializes
    // canonicalProgram and rejects a semantic/source mismatch.
    const ExpressionProgram* canonicalProgram = nullptr;
    const ExpressionProgram* runtimeProgram = nullptr;
    FormulaParameter pixelParameter = FormulaParameter::C;
    ExpressionReferenceExactInput center;
    ExpressionContext fixed;
    double bailout = 4.0;
    int maxIterations = 0;
    ExpressionReferencePrecisionPolicy precision;
    // When nonzero, compact values are additionally compared against an
    // independently iterated oracle at this higher precision. The resulting
    // radii certify the compact data relative to that finite-precision oracle.
    mpfr_prec_t certificationPrecision = 0;
    // Four terms are opt-in and reserved for consumers that explicitly prove
    // they consume every retained residual. Existing reference users retain
    // the two-term representation and behavior.
    ExpressionReferenceCompaction compaction = ExpressionReferenceCompaction::TwoTerm;
    size_t memoryLimitBytes = size_t{1} << 30;
    std::function<bool()> shouldCancel;
};

struct ExpressionReferenceTapeNode {
    // outputError/auxiliaryError are upward maximum-component radii from the
    // fast reconstructed value to the independently iterated certification
    // oracle when certifiedAgainstHigherPrecision is true.
    ScaledComplexShadow output;
    ScaledComplexShadow outputDefect;
    ScaledRealValue outputError;
    ScaledComplexShadow auxiliary;
    ScaledComplexShadow auxiliaryDefect;
    ScaledRealValue auxiliaryError;
    uint16_t leftNode = UINT16_MAX;
    uint16_t rightNode = UINT16_MAX;
    uint16_t flags = OracleTraceNone;
    uint8_t argument = 0;
    ExpressionOracleOperation operation = ExpressionOracleOperation::Constant;
    ExpressionOracleCutLocation cut = ExpressionOracleCutLocation::NotApplicable;
    ExpressionOraclePointClearance clearance = ExpressionOraclePointClearance::NotApplicable;
    ExpressionOracleCertification certification = ExpressionOracleCertification::PointOnlyNotCertified;
};

struct ExpressionReferenceSample {
    int iteration = 0;
    ScaledComplexShadow z;
    ScaledComplexShadow zDefect;
    ScaledRealValue zError;
    ScaledComplexShadow next;
    // rootDefect is only the compact quantization correction from the MPFR
    // trace root to next. It is not a full formula-recurrence defect. The
    // centered evaluator includes its reconstruction in both the state and
    // the propagated radius rather than treating primary shadows as exact.
    ScaledComplexShadow rootDefect;
    ScaledRealValue nextError;
    ScaledRealValue rootError;
    uint64_t tapeOffset = 0;
    uint16_t tapeCount = 0;
    uint16_t rootNode = UINT16_MAX;
};

struct ExpressionReferenceTapeNodeFourTerm {
    ScaledComplexExpansionTail output;
    ScaledComplexExpansionTail auxiliary;
};

struct ExpressionReferenceSampleFourTerm {
    ScaledComplexExpansionTail z;
    ScaledComplexExpansionTail next;
};

struct ExpressionReferenceFourTermData {
    ScaledComplexExpansionTail c;
    ScaledComplexExpansionTail z0;
    ScaledComplexExpansionTail pixel;
    ScaledComplexExpansionTail initialZ;
    std::vector<ExpressionReferenceSampleFourTerm> samples;
    std::vector<ExpressionReferenceTapeNodeFourTerm> tape;
};

struct ExpressionReferenceOrbitResult {
    ExpressionReferenceBuildStatus status = ExpressionReferenceBuildStatus::InvalidRequest;
    std::string error;
    bool valid = false;
    bool escaped = false;
    bool undefined = false;
    bool cancelled = false;
    // false means the compact root quantization correction is present for
    // every retained sample. This does not certify a recurrence neighborhood.
    bool defectPending = false;
    int escapeIteration = -1;
    int undefinedIteration = -1;
    mpfr_prec_t precision = 0;
    mpfr_prec_t certificationPrecision = 0;
    bool certifiedAgainstHigherPrecision = false;
    uint64_t programSemanticHash = 0;
    uint64_t canonicalSemanticHash = 0;
    std::string programSource;
    size_t sampleCount = 0;
    size_t memoryBytes = 0;
    ExpressionOracleCertification branchCertification = ExpressionOracleCertification::PointOnlyNotCertified;
    ExpressionReferenceCompaction compaction = ExpressionReferenceCompaction::TwoTerm;

    ScaledComplexShadow c;
    ScaledComplexShadow cDefect;
    ScaledRealValue cError;
    ScaledComplexShadow z0;
    ScaledComplexShadow z0Defect;
    ScaledRealValue z0Error;
    ScaledComplexShadow pixel;
    ScaledComplexShadow pixelDefect;
    ScaledRealValue pixelError;
    ScaledComplexShadow initialZ;
    ScaledComplexShadow initialZDefect;
    ScaledRealValue initialZError;
    std::vector<ExpressionReferenceSample> samples;
    std::vector<ExpressionReferenceTapeNode> tape;
    // Immutable sidecar keeps default sample/tape layouts unchanged. It is
    // present exactly when compaction is FourTermCertifiedTransfer.
    std::shared_ptr<const ExpressionReferenceFourTermData> fourTerm;
};

// This builds reference/tape data only. It is intentionally not connected to
// production GUI dispatch. Arithmetic references can be certified relative to
// a requested higher-precision finite MPFR iteration. Branch and pole metadata
// still classify only the MPFR point. Certified Taylor construction must
// independently enclose the full parameter frame; this builder never promotes
// point cut metadata into a neighborhood certificate.
//
// The intended low-precision base is the two-term compact value
// output+outputDefect stored at every tape node, not a fresh evaluation of the
// formula from rounded z shadows. The centered evaluator propagates the
// perturbation and a conservative maximum-component radius, reconstructing as
// next+rootDefect+perturbation. Auxiliary companions use the same two-term
// convention. This avoids assuming that primary shadows form an exact orbit.
bool buildExpressionReferenceOrbit(const ExpressionReferenceBuildRequest& request, ExpressionReferenceOrbitResult& result);

} // namespace formula

#endif
