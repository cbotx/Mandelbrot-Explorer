#ifndef MANDEL_FORMULA_REFERENCE_ORBIT_H
#define MANDEL_FORMULA_REFERENCE_ORBIT_H

#include <gmp.h>
#include <mpfr.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
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
    bool isInfinity() const {
        return std::isinf(mantissa);
    }
    bool isFinite() const {
        return std::isfinite(mantissa) &&
               mantissa != 0.0;
    }
};

struct ScaledComplexShadow {
    ScaledRealShadow re;
    ScaledRealShadow im;
};

// Conversion is fallible: every finite result is guaranteed to be
// reconstructible in the linked MPFR runtime exponent range.
bool makeScaledRealShadow(
    mpfr_srcptr value, ScaledRealShadow& output);
bool makeScaledComplexShadow(
    const MpfrComplex& value, ScaledComplexShadow& output);
bool setMpfrFromScaledShadow(
    mpfr_ptr output, const ScaledRealShadow& shadow,
    mpfr_rnd_t rounding = MPFR_RNDN);
bool setMpfrFromScaledShadow(
    MpfrComplex& output, const ScaledComplexShadow& shadow,
    mpfr_rnd_t rounding = MPFR_RNDN);
bool reconstructMpfrFromShadows(
    MpfrComplex& output, const ScaledComplexShadow& primary,
    const ScaledComplexShadow& defect,
    mpfr_rnd_t rounding = MPFR_RNDN);

struct ExpressionReferenceExactInput {
    // Select exactly one representation. Both mpf pointers must be non-null,
    // or realDecimal must be non-empty (an empty imaginaryDecimal means zero).
    // Any nonempty decimal field selects the decimal representation, so decimal
    // text can never be silently ignored when MPF pointers are also present.
    mpf_srcptr realMpf = nullptr;
    mpf_srcptr imaginaryMpf = nullptr;
    std::string realDecimal;
    std::string imaginaryDecimal;

    bool usesMpf() const {
        return realMpf != nullptr || imaginaryMpf != nullptr;
    }
    bool usesDecimal() const {
        return !realDecimal.empty() ||
               !imaginaryDecimal.empty();
    }
};

struct ExpressionReferencePrecisionPolicy {
    // Hard builder cap, independent of a caller's narrower maximumBits policy.
    static constexpr mpfr_prec_t ApplicationMaximumBits =
        1 << 16;

    mpfr_prec_t requestedBits = 0;
    mpfr_prec_t viewBits = 0;
    mpfr_prec_t guardBits = 64;
    mpfr_prec_t minimumBits = 128;
    mpfr_prec_t maximumBits = 1 << 20;
};

enum class ExpressionReferenceBuildStatus : uint8_t {
    Success,
    InvalidRequest,
    ProgramMismatch,
    UnsupportedProgram,
    PrecisionOutOfRange,
    InputParseError,
    CompactionOutOfRange,
    ResourceLimit
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
    size_t memoryLimitBytes = size_t{ 1 } << 30;
    std::function<bool()> shouldCancel;
};

struct ExpressionReferenceTapeNode {
    ScaledComplexShadow output;
    ScaledComplexShadow outputDefect;
    ScaledComplexShadow auxiliary;
    ScaledComplexShadow auxiliaryDefect;
    uint16_t leftNode = UINT16_MAX;
    uint16_t rightNode = UINT16_MAX;
    uint16_t flags = OracleTraceNone;
    uint8_t argument = 0;
    ExpressionOracleOperation operation =
        ExpressionOracleOperation::Constant;
    ExpressionOracleCutLocation cut =
        ExpressionOracleCutLocation::NotApplicable;
    ExpressionOraclePointClearance clearance =
        ExpressionOraclePointClearance::NotApplicable;
    ExpressionOracleCertification certification =
        ExpressionOracleCertification::PointOnlyNotCertified;
};

struct ExpressionReferenceSample {
    int iteration = 0;
    ScaledComplexShadow z;
    ScaledComplexShadow zDefect;
    ScaledComplexShadow next;
    // rootDefect is only the compact quantization correction from the MPFR
    // trace root to next. It is not a full formula-recurrence defect. A future
    // centered evaluator must include it when rebasing onto next rather than
    // treating the rounded primary shadows as an exact orbit.
    ScaledComplexShadow rootDefect;
    uint64_t tapeOffset = 0;
    uint16_t tapeCount = 0;
    uint16_t rootNode = UINT16_MAX;
};

struct ExpressionReferenceOrbitResult {
    ExpressionReferenceBuildStatus status =
        ExpressionReferenceBuildStatus::InvalidRequest;
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
    uint64_t programSemanticHash = 0;
    uint64_t canonicalSemanticHash = 0;
    std::string programSource;
    size_t sampleCount = 0;
    size_t memoryBytes = 0;
    ExpressionOracleCertification branchCertification =
        ExpressionOracleCertification::PointOnlyNotCertified;

    ScaledComplexShadow c;
    ScaledComplexShadow cDefect;
    ScaledComplexShadow z0;
    ScaledComplexShadow z0Defect;
    ScaledComplexShadow pixel;
    ScaledComplexShadow pixelDefect;
    ScaledComplexShadow initialZ;
    ScaledComplexShadow initialZDefect;
    std::vector<ExpressionReferenceSample> samples;
    std::vector<ExpressionReferenceTapeNode> tape;
};

// This builds reference/tape data only. It is intentionally not connected to
// production GUI dispatch. Branch and pole metadata classify the MPFR point;
// interval-certified neighborhood clearance remains future work.
//
// The intended low-precision base is the two-term compact value
// output+outputDefect stored at every tape node, not a fresh evaluation of the
// formula from rounded z shadows. A future centered evaluator propagates only
// the perturbation from that stored reference and reconstructs the root as
// next+rootDefect+perturbation. Auxiliary companions use the same two-term
// convention. This avoids assuming that primary shadows form an exact orbit.
bool buildExpressionReferenceOrbit(
    const ExpressionReferenceBuildRequest& request,
    ExpressionReferenceOrbitResult& result);

} // namespace formula

#endif
