#ifndef MANDEL_FORMULA_EXPRESSION_ORBIT_H
#define MANDEL_FORMULA_EXPRESSION_ORBIT_H

#include <functional>
#include <vector>

#include "formula_expression.h"

namespace formula {

struct ExpressionOrbitSnapshot {
    ExpressionProgram program;
    ExpressionContext fixed;
    FormulaParameter pixelParameter = FormulaParameter::C;
    double bailout = 4.0;

    bool valid() const;
};

struct ExpressionOrbitEvaluation {
    std::vector<Complex> points;
    int iterations = 0;
    bool escaped = false;
    bool cancelled = false;
};

struct ExpressionOrbitClassification {
    int iterations = 0;
    bool escaped = false;
    bool cancelled = false;
};

bool evaluateExpressionOrbit(const ExpressionOrbitSnapshot& snapshot, Complex pixel, int maxIterations, ExpressionOrbitEvaluation& result, const std::function<bool()>& shouldCancel = {});

bool classifyExpressionOrbit(const ExpressionOrbitSnapshot& snapshot, Complex pixel, int maxIterations, ExpressionOrbitClassification& result, const std::function<bool()>& shouldCancel = {});

} // namespace formula

#endif
