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

    bool success() const { return status == ExpressionCenteredStatus::Success; }
};

const char* expressionCenteredStatusName(ExpressionCenteredStatus status);

class ExpressionCenteredEvaluator {
  public:
    // A non-success status is an explicit request to re-reference or use a
    // higher-precision fallback; this evaluator never substitutes one.
    static ExpressionCenteredResult evaluate(const ExpressionProgram& program, const ExpressionContext& reference, const ExpressionDeltaContext& delta);
};

} // namespace formula

#endif
