#ifndef MANDEL_FORMULA_EXPRESSION_ORACLE_H
#define MANDEL_FORMULA_EXPRESSION_ORACLE_H

#include <mpfr.h>

#include <array>
#include <complex>
#include <cstdint>
#include <string>
#include <vector>

#include "formula_expression.h"

namespace formula {

class MpfrComplex {
public:
    explicit MpfrComplex(mpfr_prec_t precision = 256);
    MpfrComplex(const MpfrComplex& other);
    MpfrComplex& operator=(const MpfrComplex& other);
    ~MpfrComplex();

    void set(double real, double imaginary = 0.0);
    bool set(const std::string& real, const std::string& imaginary = "0");
    void set(const MpfrComplex& other);
    std::complex<double> toDouble() const;
    mpfr_prec_t precision() const;

    mpfr_t re;
    mpfr_t im;
};

struct ExpressionOracleContext {
    explicit ExpressionOracleContext(mpfr_prec_t precision = 256);

    MpfrComplex z;
    MpfrComplex c;
    MpfrComplex z0;
    std::array<MpfrComplex, 8> parameters;
    int iteration = 0;
};

enum class ExpressionOracleOperation : uint8_t {
    Constant,
    Z, C, Z0, Iteration, Parameter,
    OrbitInvariant,
    Negate, Add, Subtract, Multiply, Divide, Power,
    Square,
    Sin, Cos, Tan, Sinh, Cosh, Tanh,
    Exp, Log, Log10, Sqrt,
    Abs, Norm, Arg, Conjugate, Real, Imaginary,
    MakeComplex, Polar
};

enum ExpressionOracleTraceFlag : uint16_t {
    OracleTraceNone = 0,
    OracleTraceNonlinear = 1 << 0,
    OracleTraceBranchSensitive = 1 << 1,
    OracleTraceHasCompanion = 1 << 2,
    OracleTraceHasDenominator = 1 << 3,
    OracleTraceHasLogarithmBase = 1 << 4,
    OracleTraceUndefined = 1 << 5,
    OracleTraceSingularPoint = 1 << 6,
    // For the large-argument tan/tanh branch, auxiliary stores
    // (logScale, reducedPhase), where q=exp(2*logScale). reducedPhase is
    // 2*Re(input) for tan or 2*Im(input) for tanh, reduced by atan2 to
    // [-pi,pi]. This preserves asymptotically tiny sensitivity without
    // requiring q to fit MPFR's runtime exponent range.
    OracleTraceHasAsymptoticLogPhase = 1 << 7
};

enum class ExpressionOracleCutLocation : uint8_t {
    NotApplicable,
    UpperHalfPlane,
    LowerHalfPlane,
    PositiveRealAxis,
    NegativeRealUpperLip,
    NegativeRealLowerLip,
    Origin
};

enum class ExpressionOraclePointClearance : uint8_t {
    NotApplicable,
    ClearAtPoint,
    ZeroAtPoint,
    NonFiniteAtPoint
};

enum class ExpressionOracleCertification : uint8_t {
    // This is a classification of the exact reference point only. It is not an
    // interval proof that a future residual neighborhood avoids a cut or pole.
    PointOnlyNotCertified
};

struct ExpressionOracleTraceNode {
    explicit ExpressionOracleTraceNode(mpfr_prec_t precision);

    size_t instructionIndex = 0;
    uint8_t argument = 0;
    ExpressionOracleOperation operation =
        ExpressionOracleOperation::Constant;
    uint16_t flags = OracleTraceNone;
    uint16_t leftNode = UINT16_MAX;
    uint16_t rightNode = UINT16_MAX;
    ExpressionOracleCutLocation cut =
        ExpressionOracleCutLocation::NotApplicable;
    ExpressionOraclePointClearance clearance =
        ExpressionOraclePointClearance::NotApplicable;
    ExpressionOracleCertification certification =
        ExpressionOracleCertification::PointOnlyNotCertified;
    MpfrComplex output;
    // sin/cos and sinh/cosh companions, tan/tanh or divide denominators, or
    // log(base) for power according to flags. When
    // OracleTraceHasAsymptoticLogPhase is set for tan/tanh, auxiliary instead
    // stores (logScale, reducedPhase); OracleTraceHasDenominator then refers
    // to the stable scaled denominator 1+q^2+2q*cos(reducedPhase).
    MpfrComplex auxiliary;
};

struct ExpressionOracleTrace {
    mpfr_prec_t precision = 0;
    bool exactDomain = false;
    std::vector<ExpressionOracleTraceNode> nodes;
};

class ExpressionOracle {
public:
    static bool evaluate(const ExpressionProgram& program,
                         const ExpressionOracleContext& context,
                         MpfrComplex& output,
                         std::string* error = nullptr);
    static bool evaluateTrace(const ExpressionProgram& program,
                              const ExpressionOracleContext& context,
                              MpfrComplex& output,
                              ExpressionOracleTrace& trace,
                              std::string* error = nullptr);
    // Releases the reusable evaluator owned by the calling thread.
    static void releaseThreadWorkspace();

private:
    static bool evaluateInternal(const ExpressionProgram& program,
                                 const ExpressionOracleContext& context,
                                 MpfrComplex& output,
                                 ExpressionOracleTrace* trace,
                                 std::string* error);
};

} // namespace formula

#endif
