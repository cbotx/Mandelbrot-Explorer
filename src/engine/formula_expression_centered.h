#ifndef MANDEL_FORMULA_EXPRESSION_CENTERED_H
#define MANDEL_FORMULA_EXPRESSION_CENTERED_H

#include <array>
#include <cstdint>

#include "formula_expression.h"

namespace formula {

struct ExpressionDeltaContext {
    Complex z{};
    Complex c{};
    Complex z0{};
    std::array<Complex, 8> parameters{};
};

enum class ExpressionCenteredStatus : uint8_t { Success,
                                                BranchUncertain,
                                                Singular,
                                                Undefined,
                                                Unsupported,
                                                NonFinite,
                                                InvalidProgram };

struct ExpressionCenteredResult {
    Complex base{};
    Complex delta{};
    ExpressionCenteredStatus status = ExpressionCenteredStatus::InvalidProgram;
    bool denominatorInstability = false;

    bool success() const { return status == ExpressionCenteredStatus::Success; }
};

const char* expressionCenteredStatusName(ExpressionCenteredStatus status);

class ExpressionCenteredEvaluator {
  public:
    // A non-success status is an explicit request to re-reference or use a
    // higher-precision fallback; this evaluator never substitutes one.
    static ExpressionCenteredResult evaluate(const ExpressionProgram& program, const ExpressionContext& reference, const ExpressionDeltaContext& delta);
    static bool evaluate4(const ExpressionProgram& program, const ExpressionContext& reference, const ExpressionDeltaContext* deltas, ExpressionCenteredResult* results);
    static bool evaluate4WithNodeBases(const ExpressionProgram& program, const ExpressionContext& reference, const Complex* nodeBases, const ExpressionDeltaContext* deltas, ExpressionCenteredResult* results);
    // Approximate AVX2 analytic execution for the adaptive candidate.
    // Strict centered/scaled/Taylor consumers must use the APIs above.
    static bool evaluate4FastEntire(const ExpressionProgram& program, const ExpressionContext& reference, const ExpressionDeltaContext* deltas, ExpressionCenteredResult* results);
    static bool evaluate4WithNodeBasesFastEntire(const ExpressionProgram& program, const ExpressionContext& reference, const Complex* nodeBases, const Complex* nodeAuxiliaries, const ExpressionDeltaContext* deltas, ExpressionCenteredResult* results);
    // The final array contains real 2x2 state-Jacobian operator norms. For
    // holomorphic programs this equals the traditional complex derivative
    // magnitude; real-smooth operations use their full bivariate Jacobian.
    static bool evaluate4WithNodeBasesFastEntireDerivative(const ExpressionProgram& program, const ExpressionContext& reference, const Complex* nodeBases, const Complex* nodeAuxiliaries, const ExpressionDeltaContext* deltas, ExpressionCenteredResult* results, double* jacobianNorms);
    static bool evaluate4WithLaneNodeBasesFastEntire(const ExpressionProgram& program, const ExpressionContext* const* references, const Complex* const* nodeBases, const Complex* const* nodeAuxiliaries, const ExpressionDeltaContext* const* deltas, uint8_t activeMask, ExpressionCenteredResult* results);
    static bool evaluate4WithLaneNodeBasesFastEntireDerivative(const ExpressionProgram& program, const ExpressionContext* const* references, const Complex* const* nodeBases, const Complex* const* nodeAuxiliaries, const ExpressionDeltaContext* const* deltas, uint8_t activeMask, ExpressionCenteredResult* results, double* jacobianNorms);

  private:
    static bool evaluate4WithNodeBasesImpl(const ExpressionProgram& program, const ExpressionContext& reference, const Complex* nodeBases, const ExpressionDeltaContext* deltas, ExpressionCenteredResult* results);
    static bool evaluate4FastEntireImpl(const ExpressionProgram& program, const ExpressionContext& reference, const Complex* nodeBases, const Complex* nodeAuxiliaries, const ExpressionDeltaContext* deltas, ExpressionCenteredResult* results, double* derivativeMagnitudes);
    static bool evaluate4LaneFastEntireImpl(const ExpressionProgram& program, const ExpressionContext* const* references, const Complex* const* nodeBases, const Complex* const* nodeAuxiliaries, const ExpressionDeltaContext* const* deltas, uint8_t activeMask, ExpressionCenteredResult* results, double* derivativeMagnitudes);
};

} // namespace formula

#endif
