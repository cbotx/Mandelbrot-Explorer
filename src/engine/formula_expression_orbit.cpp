#include "formula_expression_orbit.h"

#include <array>
#include <cmath>

namespace formula {

bool ExpressionOrbitSnapshot::valid() const {
    return program.valid() &&
           (pixelParameter == FormulaParameter::C ||
            pixelParameter == FormulaParameter::InitialZ) &&
           bailout > 0.0 && std::isfinite(bailout);
}

bool evaluateExpressionOrbit(
        const ExpressionOrbitSnapshot& snapshot, Complex pixel,
        int maxIterations, ExpressionOrbitEvaluation& result,
        const std::function<bool()>& shouldCancel) {
    result = ExpressionOrbitEvaluation{};
    if (!snapshot.valid() || maxIterations < 0)
        return false;

    ExpressionContext context = snapshot.fixed;
    if (snapshot.pixelParameter == FormulaParameter::C)
        context.c = pixel;
    else
        context.z0 = pixel;
    context.z = context.z0;

    result.points.reserve((size_t)maxIterations + 1);
    result.points.push_back(context.z);
    if (!std::isfinite(context.z.real()) ||
        !std::isfinite(context.z.imag()) ||
        std::hypot(context.z.real(), context.z.imag()) > snapshot.bailout) {
        result.escaped = true;
        return true;
    }

    std::array<Complex, ExpressionProgram::MAX_STACK> stack;
    for (int n = 0; n < maxIterations; ++n) {
        if (shouldCancel && shouldCancel()) {
            result.cancelled = true;
            break;
        }
        context.iteration = n;
        context.z = snapshot.program.evaluate(
            context, stack.data(), snapshot.program.stackDepth());
        result.points.push_back(context.z);
        result.iterations = n + 1;
        if (!std::isfinite(context.z.real()) ||
            !std::isfinite(context.z.imag()) ||
            std::hypot(context.z.real(), context.z.imag()) > snapshot.bailout) {
            result.escaped = true;
            break;
        }
    }
    return true;
}

} // namespace formula
