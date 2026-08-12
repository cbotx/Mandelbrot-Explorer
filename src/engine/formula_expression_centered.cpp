#include "formula_expression_centered.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

#ifdef _MSC_VER
#pragma float_control(precise, on, push)
#pragma fp_contract(off)
#endif

namespace formula {

namespace {

struct CenteredValue {
    Complex base{};
    Complex delta{};
    Complex endpoint{};
    ExpressionCenteredStatus status = ExpressionCenteredStatus::Success;
};

constexpr double PI = 3.141592653589793238462643383279502884;

bool finite(Complex value) {
    return std::isfinite(value.real()) && std::isfinite(value.imag());
}

bool zero(Complex value) {
    return value.real() == 0.0 && value.imag() == 0.0;
}

Complex nanComplex() {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    return {nan, nan};
}

CenteredValue finish(Complex base, Complex delta, Complex endpoint, ExpressionCenteredStatus status = ExpressionCenteredStatus::Success) {
    if (status == ExpressionCenteredStatus::Success && (!finite(base) || !finite(delta) || !finite(endpoint))) status = ExpressionCenteredStatus::NonFinite;
    return {base, delta, endpoint, status};
}

ExpressionCenteredStatus combine(ExpressionCenteredStatus left, ExpressionCenteredStatus right) {
    return left != ExpressionCenteredStatus::Success ? left : right;
}

bool segmentMayContainZero(Complex base, Complex endpoint) {
    if (zero(base) || zero(endpoint)) return true;
    const double minReal = std::min(base.real(), endpoint.real());
    const double maxReal = std::max(base.real(), endpoint.real());
    const double minImaginary = std::min(base.imag(), endpoint.imag());
    const double maxImaginary = std::max(base.imag(), endpoint.imag());
    const long double dx = (long double)endpoint.real() - (long double)base.real();
    const long double dy = (long double)endpoint.imag() - (long double)base.imag();
    if (dx != 0.0L || dy != 0.0L) {
        const long double x0 = base.real();
        const long double y0 = base.imag();
        if (dx != 0.0L) {
            const long double t = -x0 / dx;
            if (t >= 0.0L && t <= 1.0L && y0 + t * dy == 0.0L) return true;
        } else if (x0 == 0.0L && dy != 0.0L) {
            const long double t = -y0 / dy;
            if (t >= 0.0L && t <= 1.0L) return true;
        }
    }
    return minReal <= 0.0 && maxReal >= 0.0 && minImaginary <= 0.0 && maxImaginary >= 0.0;
}

struct CutCheck {
    bool singular = false;
    bool cut = false;
};

CutCheck checkNegativeRealCut(Complex base, Complex endpoint) {
    CutCheck result;
    if (zero(base) || zero(endpoint)) {
        result.singular = true;
        return result;
    }

    const long double x0 = base.real();
    const long double y0 = base.imag();
    const long double dx = (long double)endpoint.real() - x0;
    const long double dy = (long double)endpoint.imag() - y0;
    if (dx != 0.0L) {
        const long double t = -x0 / dx;
        if (t >= 0.0L && t <= 1.0L && y0 + t * dy == 0.0L) {
            result.singular = true;
            return result;
        }
    } else if (x0 == 0.0L && dy != 0.0L) {
        const long double t = -y0 / dy;
        if (t >= 0.0L && t <= 1.0L) {
            result.singular = true;
            return result;
        }
    }

    if ((base.imag() == 0.0 && base.real() < 0.0) || (endpoint.imag() == 0.0 && endpoint.real() < 0.0)) {
        result.cut = true;
        return result;
    }
    if (dy != 0.0L) {
        const long double t = -y0 / dy;
        if (t >= 0.0L && t <= 1.0L) {
            const long double crossingReal = x0 + t * dx;
            if (crossingReal <= 0.0L) result.cut = true;
        }
    } else if (y0 == 0.0L) {
        const long double x1 = endpoint.real();
        if (std::min(x0, x1) < 0.0L) result.cut = true;
    }
    return result;
}

bool signedZeroCutSideChanges(Complex base, Complex endpoint) {
    return base.real() < 0.0 && endpoint.real() < 0.0 && base.imag() == 0.0 && endpoint.imag() == 0.0 && std::signbit(base.imag()) != std::signbit(endpoint.imag());
}

bool nearPeriodicPole(double value) {
    if (!std::isfinite(value)) return true;
    const double reduced = std::remainder(value - PI * 0.5, PI);
    const double tolerance = 64.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, std::abs(value)) + 1e-14;
    return std::abs(reduced) <= tolerance;
}

bool intervalContainsPeriodicPole(double first, double second) {
    if (!std::isfinite(first) || !std::isfinite(second)) return true;
    double low = std::min(first, second);
    double high = std::max(first, second);
    if (high - low >= PI) return true;
    double k = std::ceil((low - PI * 0.5) / PI);
    double pole = PI * 0.5 + k * PI;
    const double tolerance = 64.0 * std::numeric_limits<double>::epsilon() * std::max({1.0, std::abs(low), std::abs(high)}) + 1e-14;
    return pole <= high + tolerance;
}

bool tangentSegmentMayHitPole(Complex base, Complex endpoint) {
    const double y0 = base.imag();
    const double y1 = endpoint.imag();
    if (!(std::min(y0, y1) <= 0.0 && std::max(y0, y1) >= 0.0)) return false;
    const long double dx = (long double)endpoint.real() - (long double)base.real();
    const long double dy = (long double)endpoint.imag() - (long double)base.imag();
    if (dy == 0.0L) { return y0 == 0.0 && intervalContainsPeriodicPole(base.real(), endpoint.real()); }
    long double t = -(long double)y0 / dy;
    if (t < 0.0L || t > 1.0L) return false;
    long double crossing = (long double)base.real() + t * dx;
    return nearPeriodicPole((double)crossing);
}

bool hyperbolicTangentSegmentMayHitPole(Complex base, Complex endpoint) {
    Complex swappedBase{base.imag(), base.real()};
    Complex swappedEndpoint{endpoint.imag(), endpoint.real()};
    return tangentSegmentMayHitPole(swappedBase, swappedEndpoint);
}

Complex centeredMultiply(Complex leftBase, Complex leftDelta, Complex rightBase, Complex rightDelta) {
    // Bytecode v1 fixes this order: a*db, then b*da, then da*db.
    Complex result = leftBase * rightDelta;
    result += rightBase * leftDelta;
    result += leftDelta * rightDelta;
    return result;
}

Complex centeredSquare(Complex base, Complex delta) {
    Complex result = base * delta;
    result += base * delta;
    result += delta * delta;
    return result;
}

Complex complexExpm1(Complex value) {
    if (zero(value)) return {};
    const double x = value.real();
    const double y = value.imag();
    const double sine = std::sin(y);
    const double cosine = std::cos(y);
    if (x < -0.5) {
        const double exponential = std::exp(x);
        return {exponential * cosine - 1.0, exponential * sine};
    }
    const double realExpm1 = std::expm1(x);
    const double halfSine = std::sin(y * 0.5);
    const double cosineMinusOne = -2.0 * halfSine * halfSine;
    double real = realExpm1 * cosine;
    real += cosineMinusOne;
    double imaginary = realExpm1 * sine;
    imaginary += sine;
    return {real, imaginary};
}

Complex relativeExponentialDelta(Complex baseOutput, Complex perturbedOutput, Complex exponentDelta) {
    if (zero(exponentDelta)) return {};
    if (exponentDelta.real() > 0.0) {
        Complex relative = complexExpm1(-exponentDelta);
        return -(perturbedOutput * relative);
    }
    return baseOutput * complexExpm1(exponentDelta);
}

Complex twiceSinHalf(Complex value) {
    const double magnitude = std::max(std::abs(value.real()), std::abs(value.imag()));
    if (magnitude < 1e-8) return value;
    return 2.0 * std::sin(value * 0.5);
}

Complex twiceSinhHalf(Complex value) {
    const double magnitude = std::max(std::abs(value.real()), std::abs(value.imag()));
    if (magnitude < 1e-8) return value;
    return 2.0 * std::sinh(value * 0.5);
}

Complex centeredSine(Complex base, Complex delta) {
    Complex midpoint = base + delta * 0.5;
    return std::cos(midpoint) * twiceSinHalf(delta);
}

Complex centeredCosine(Complex base, Complex delta) {
    Complex midpoint = base + delta * 0.5;
    return -std::sin(midpoint) * twiceSinHalf(delta);
}

Complex centeredHyperbolicSine(Complex base, Complex delta) {
    Complex midpoint = base + delta * 0.5;
    return std::cosh(midpoint) * twiceSinhHalf(delta);
}

Complex centeredHyperbolicCosine(Complex base, Complex delta) {
    Complex midpoint = base + delta * 0.5;
    return std::sinh(midpoint) * twiceSinhHalf(delta);
}

ExpressionCenteredStatus centeredDivisionDelta(Complex numeratorDelta, Complex denominatorBase, Complex denominatorDelta, Complex denominatorEndpoint, Complex quotientBase, Complex& output, bool guardLinearSegment) {
    if (zero(denominatorBase) || zero(denominatorEndpoint) || (guardLinearSegment && segmentMayContainZero(denominatorBase, denominatorEndpoint))) return ExpressionCenteredStatus::Singular;
    Complex numerator = numeratorDelta - quotientBase * denominatorDelta;
    output = numerator / denominatorEndpoint;
    return finite(output) ? ExpressionCenteredStatus::Success : ExpressionCenteredStatus::NonFinite;
}

double centeredRealProduct(double leftBase, double leftDelta, double rightBase, double rightDelta) {
    double result = leftBase * rightDelta;
    result += rightBase * leftDelta;
    result += leftDelta * rightDelta;
    return result;
}

ExpressionCenteredStatus centeredTangent(Complex base, Complex delta, Complex endpoint, Complex baseOutput, Complex endpointOutput, Complex& output) {
    if (zero(delta)) {
        output = {};
        return ExpressionCenteredStatus::Success;
    }
    if (tangentSegmentMayHitPole(base, endpoint)) return ExpressionCenteredStatus::Singular;
    if (!finite(baseOutput) || !finite(endpointOutput)) return ExpressionCenteredStatus::NonFinite;

    const bool sameHighHalfPlane = std::abs(base.imag()) > 16.0 && std::abs(endpoint.imag()) > 16.0 && std::signbit(base.imag()) == std::signbit(endpoint.imag());
    if (!sameHighHalfPlane) {
        Complex sineBase = std::sin(base);
        Complex cosineBase = std::cos(base);
        Complex cosineEndpoint = std::cos(endpoint);
        Complex sineDelta = centeredSine(base, delta);
        Complex cosineDelta = centeredCosine(base, delta);
        if (!finite(sineBase) || !finite(cosineBase) || !finite(cosineEndpoint) || !finite(sineDelta) || !finite(cosineDelta)) return ExpressionCenteredStatus::NonFinite;
        return centeredDivisionDelta(sineDelta, cosineBase, cosineDelta, cosineEndpoint, baseOutput, output, false);
    }

    const double sign = std::signbit(base.imag()) ? -1.0 : 1.0;
    const double q = std::exp(-2.0 * std::abs(base.imag()));
    const double q1 = std::exp(-2.0 * std::abs(endpoint.imag()));
    const double qDelta = relativeExponentialDelta({q, 0.0}, {q1, 0.0}, {-2.0 * sign * delta.imag(), 0.0}).real();
    const double angle = 2.0 * base.real();
    const double angleDelta = 2.0 * delta.real();
    const double sine = std::sin(angle);
    const double cosine = std::cos(angle);
    const double sineDelta = centeredSine({angle, 0.0}, {angleDelta, 0.0}).real();
    const double cosineDelta = centeredCosine({angle, 0.0}, {angleDelta, 0.0}).real();
    double qSquaredDelta = q * qDelta;
    qSquaredDelta += q * qDelta;
    qSquaredDelta += qDelta * qDelta;
    const double numeratorRealDelta = 2.0 * centeredRealProduct(q, qDelta, sine, sineDelta);
    const double numeratorImaginaryDelta = -sign * qSquaredDelta;
    double denominatorDelta = qSquaredDelta;
    denominatorDelta += 2.0 * centeredRealProduct(q, qDelta, cosine, cosineDelta);
    double denominator = 1.0 + q * q;
    denominator += 2.0 * q * cosine;
    double endpointDenominator = 1.0 + q1 * q1;
    endpointDenominator += 2.0 * q1 * std::cos(2.0 * endpoint.real());
    return centeredDivisionDelta({numeratorRealDelta, numeratorImaginaryDelta}, {denominator, 0.0}, {denominatorDelta, 0.0}, {endpointDenominator, 0.0}, baseOutput, output, false);
}

ExpressionCenteredStatus centeredHyperbolicTangent(Complex base, Complex delta, Complex endpoint, Complex baseOutput, Complex endpointOutput, Complex& output) {
    if (zero(delta)) {
        output = {};
        return ExpressionCenteredStatus::Success;
    }
    if (hyperbolicTangentSegmentMayHitPole(base, endpoint)) return ExpressionCenteredStatus::Singular;
    if (!finite(baseOutput) || !finite(endpointOutput)) return ExpressionCenteredStatus::NonFinite;

    const bool sameHighHalfPlane = std::abs(base.real()) > 16.0 && std::abs(endpoint.real()) > 16.0 && std::signbit(base.real()) == std::signbit(endpoint.real());
    if (!sameHighHalfPlane) {
        Complex sineBase = std::sinh(base);
        Complex cosineBase = std::cosh(base);
        Complex cosineEndpoint = std::cosh(endpoint);
        Complex sineDelta = centeredHyperbolicSine(base, delta);
        Complex cosineDelta = centeredHyperbolicCosine(base, delta);
        if (!finite(sineBase) || !finite(cosineBase) || !finite(cosineEndpoint) || !finite(sineDelta) || !finite(cosineDelta)) return ExpressionCenteredStatus::NonFinite;
        return centeredDivisionDelta(sineDelta, cosineBase, cosineDelta, cosineEndpoint, baseOutput, output, false);
    }

    const double sign = std::signbit(base.real()) ? -1.0 : 1.0;
    const double q = std::exp(-2.0 * std::abs(base.real()));
    const double q1 = std::exp(-2.0 * std::abs(endpoint.real()));
    const double qDelta = relativeExponentialDelta({q, 0.0}, {q1, 0.0}, {-2.0 * sign * delta.real(), 0.0}).real();
    const double angle = 2.0 * base.imag();
    const double angleDelta = 2.0 * delta.imag();
    const double sine = std::sin(angle);
    const double cosine = std::cos(angle);
    const double sineDelta = centeredSine({angle, 0.0}, {angleDelta, 0.0}).real();
    const double cosineDelta = centeredCosine({angle, 0.0}, {angleDelta, 0.0}).real();
    double qSquaredDelta = q * qDelta;
    qSquaredDelta += q * qDelta;
    qSquaredDelta += qDelta * qDelta;
    const double numeratorRealDelta = -sign * qSquaredDelta;
    const double numeratorImaginaryDelta = 2.0 * centeredRealProduct(q, qDelta, sine, sineDelta);
    double denominatorDelta = qSquaredDelta;
    denominatorDelta += 2.0 * centeredRealProduct(q, qDelta, cosine, cosineDelta);
    double denominator = 1.0 + q * q;
    denominator += 2.0 * q * cosine;
    double endpointDenominator = 1.0 + q1 * q1;
    endpointDenominator += 2.0 * q1 * std::cos(2.0 * endpoint.imag());
    return centeredDivisionDelta({numeratorRealDelta, numeratorImaginaryDelta}, {denominator, 0.0}, {denominatorDelta, 0.0}, {endpointDenominator, 0.0}, baseOutput, output, false);
}

ExpressionCenteredStatus centeredLogarithm(Complex base, Complex delta, Complex endpoint, Complex& output) {
    if (zero(base) || zero(endpoint)) return ExpressionCenteredStatus::Singular;
    if (signedZeroCutSideChanges(base, endpoint)) return ExpressionCenteredStatus::BranchUncertain;
    if (zero(delta)) {
        output = {};
        return ExpressionCenteredStatus::Success;
    }
    CutCheck cut = checkNegativeRealCut(base, endpoint);
    if (cut.singular) return ExpressionCenteredStatus::Singular;
    if (cut.cut) return ExpressionCenteredStatus::BranchUncertain;
    Complex ratio = delta / base;
    if (!finite(ratio)) return ExpressionCenteredStatus::NonFinite;
    const double x = ratio.real();
    const double y = ratio.imag();
    const double onePlusX = 1.0 + x;
    const double magnitude = std::hypot(onePlusX, y);
    if (magnitude == 0.0) return ExpressionCenteredStatus::Singular;

    double real;
    if (std::max(std::abs(x), std::abs(y)) <= 0.5) {
        double squared = x * x;
        squared += y * y;
        double logArgument = x + x;
        logArgument += squared;
        if (logArgument <= -1.0) {
            real = std::log(magnitude);
        } else {
            real = 0.5 * std::log1p(logArgument);
        }
    } else {
        real = std::log(magnitude);
    }
    output = {real, std::atan2(y, onePlusX)};
    return finite(output) ? ExpressionCenteredStatus::Success : ExpressionCenteredStatus::NonFinite;
}

Complex centeredNormDelta(Complex base, Complex delta) {
    Complex cross = std::conj(base) * delta;
    double real = cross.real() + cross.real();
    real += std::norm(delta);
    return {real, 0.0};
}

ExpressionCenteredStatus centeredAbsolute(Complex base, Complex delta, Complex endpoint, Complex& output) {
    const double baseRadius = std::abs(base);
    const double perturbedRadius = std::abs(endpoint);
    const double scale = std::max(baseRadius, perturbedRadius);
    if (!std::isfinite(scale)) return ExpressionCenteredStatus::NonFinite;
    if (scale == 0.0) {
        output = {};
        return ExpressionCenteredStatus::Success;
    }
    Complex scaledBase = base / scale;
    Complex scaledDelta = delta / scale;
    double numerator = centeredNormDelta(scaledBase, scaledDelta).real();
    double denominator = baseRadius / scale + perturbedRadius / scale;
    if (denominator == 0.0) return ExpressionCenteredStatus::Singular;
    double residual = numerator / denominator;
    residual *= scale;
    output = {residual, 0.0};
    return finite(output) ? ExpressionCenteredStatus::Success : ExpressionCenteredStatus::NonFinite;
}

ExpressionCenteredStatus centeredSquareRoot(Complex base, Complex delta, Complex endpoint, Complex baseRoot, Complex endpointRoot, Complex& output) {
    if (zero(base) && zero(endpoint)) return ExpressionCenteredStatus::Singular;
    if (signedZeroCutSideChanges(base, endpoint)) return ExpressionCenteredStatus::BranchUncertain;
    if (zero(delta)) {
        output = {};
        return ExpressionCenteredStatus::Success;
    }
    CutCheck cut = checkNegativeRealCut(base, endpoint);
    if (cut.cut) return ExpressionCenteredStatus::BranchUncertain;
    Complex denominator = endpointRoot + baseRoot;
    if (zero(denominator) || cut.singular) return ExpressionCenteredStatus::Singular;
    output = delta / denominator;
    return finite(output) ? ExpressionCenteredStatus::Success : ExpressionCenteredStatus::NonFinite;
}

ExpressionCenteredStatus centeredPower(Complex base, Complex baseDelta, Complex baseEndpoint, Complex exponent, Complex exponentDelta, Complex exponentEndpoint, Complex baseOutput, Complex endpointOutput, Complex& output) {
    if (zero(base) || zero(baseEndpoint)) return ExpressionCenteredStatus::Singular;
    if (signedZeroCutSideChanges(base, baseEndpoint)) return ExpressionCenteredStatus::BranchUncertain;
    if (zero(baseDelta) && zero(exponentDelta)) {
        output = {};
        return ExpressionCenteredStatus::Success;
    }
    CutCheck cut = checkNegativeRealCut(base, baseEndpoint);
    if (cut.singular) return ExpressionCenteredStatus::Singular;
    if (cut.cut) return ExpressionCenteredStatus::BranchUncertain;

    Complex logBase = std::log(base);
    Complex logDelta;
    ExpressionCenteredStatus status = centeredLogarithm(base, baseDelta, baseEndpoint, logDelta);
    if (status != ExpressionCenteredStatus::Success) return status;
    Complex productDelta = centeredMultiply(exponent, exponentDelta, logBase, logDelta);
    if (!finite(exponentEndpoint) || !finite(baseOutput) || !finite(endpointOutput)) return ExpressionCenteredStatus::NonFinite;
    output = relativeExponentialDelta(baseOutput, endpointOutput, productDelta);
    return finite(output) ? ExpressionCenteredStatus::Success : ExpressionCenteredStatus::NonFinite;
}

ExpressionCenteredStatus centeredPolar(const CenteredValue& radius, const CenteredValue& angle, Complex baseOutput, Complex endpointOutput, Complex& output) {
    if (radius.base.imag() != 0.0 || radius.delta.imag() != 0.0 || radius.endpoint.imag() != 0.0 || angle.base.imag() != 0.0 || angle.delta.imag() != 0.0 || angle.endpoint.imag() != 0.0) return ExpressionCenteredStatus::Undefined;
    const double baseRadius = radius.base.real();
    const double radiusDelta = radius.delta.real();
    const double baseAngle = angle.base.real();
    const double angleDelta = angle.delta.real();
    const double perturbedRadius = radius.endpoint.real();
    const double perturbedAngle = angle.endpoint.real();
    if (!std::isfinite(baseRadius) || !std::isfinite(radiusDelta) || !std::isfinite(baseAngle) || !std::isfinite(angleDelta) || !std::isfinite(perturbedRadius) || !std::isfinite(perturbedAngle)) return ExpressionCenteredStatus::NonFinite;
    if (baseRadius < 0.0 || perturbedRadius < 0.0) return ExpressionCenteredStatus::Undefined;
    if (!finite(baseOutput) || !finite(endpointOutput)) return ExpressionCenteredStatus::NonFinite;
    Complex angleResidual = complexExpm1({0.0, angleDelta});
    Complex result = baseOutput * angleResidual;
    result += radiusDelta * std::polar(1.0, perturbedAngle);
    output = result;
    return finite(output) ? ExpressionCenteredStatus::Success : ExpressionCenteredStatus::NonFinite;
}

} // namespace

const char* expressionCenteredStatusName(ExpressionCenteredStatus status) {
    switch (status) {
    case ExpressionCenteredStatus::Success: return "success";
    case ExpressionCenteredStatus::BranchUncertain: return "branch-uncertain";
    case ExpressionCenteredStatus::Singular: return "singular";
    case ExpressionCenteredStatus::Undefined: return "undefined";
    case ExpressionCenteredStatus::Unsupported: return "unsupported";
    case ExpressionCenteredStatus::NonFinite: return "nonfinite";
    case ExpressionCenteredStatus::InvalidProgram: return "invalid-program";
    }
    return "invalid-program";
}

ExpressionCenteredResult ExpressionCenteredEvaluator::evaluate(const ExpressionProgram& program, const ExpressionContext& reference, const ExpressionDeltaContext& inputDelta) {
    ExpressionCenteredResult invalid;
    invalid.base = nanComplex();
    invalid.delta = nanComplex();
    if (!program._valid) return invalid;

    std::array<CenteredValue, ExpressionProgram::MAX_STACK> stack;
    size_t top = 0;
    auto pushUnchanged = [&](Complex value) { stack[top++] = finish(value, {}, value); };
    auto pushPerturbed = [&](Complex base, Complex delta) {
        Complex endpoint = base;
        endpoint += delta;
        stack[top++] = finish(base, delta, endpoint);
    };
    auto unary = [&](const ExpressionProgram::Instruction& instruction) {
        CenteredValue input = stack[top - 1];
        Complex base = ExpressionProgram::evaluateUnary(instruction.op, input.base);
        Complex endpoint = ExpressionProgram::evaluateUnary(instruction.op, input.endpoint);
        if (input.status != ExpressionCenteredStatus::Success) {
            stack[top - 1] = finish(base, nanComplex(), endpoint, input.status);
            return;
        }

        Complex delta;
        ExpressionCenteredStatus status = ExpressionCenteredStatus::Success;
        switch (instruction.op) {
        case ExpressionProgram::Op::Negate: delta = -input.delta; break;
        case ExpressionProgram::Op::Square: delta = centeredSquare(input.base, input.delta); break;
        case ExpressionProgram::Op::Sin: delta = centeredSine(input.base, input.delta); break;
        case ExpressionProgram::Op::Cos: delta = centeredCosine(input.base, input.delta); break;
        case ExpressionProgram::Op::Tan: status = centeredTangent(input.base, input.delta, input.endpoint, base, endpoint, delta); break;
        case ExpressionProgram::Op::Sinh: delta = centeredHyperbolicSine(input.base, input.delta); break;
        case ExpressionProgram::Op::Cosh: delta = centeredHyperbolicCosine(input.base, input.delta); break;
        case ExpressionProgram::Op::Tanh: status = centeredHyperbolicTangent(input.base, input.delta, input.endpoint, base, endpoint, delta); break;
        case ExpressionProgram::Op::Exp: {
            if (!finite(base) || !finite(endpoint)) {
                status = ExpressionCenteredStatus::NonFinite;
            } else {
                delta = relativeExponentialDelta(base, endpoint, input.delta);
            }
            break;
        }
        case ExpressionProgram::Op::Log: status = centeredLogarithm(input.base, input.delta, input.endpoint, delta); break;
        case ExpressionProgram::Op::Log10:
            status = centeredLogarithm(input.base, input.delta, input.endpoint, delta);
            if (status == ExpressionCenteredStatus::Success) delta /= std::log(10.0);
            break;
        case ExpressionProgram::Op::Sqrt: status = centeredSquareRoot(input.base, input.delta, input.endpoint, base, endpoint, delta); break;
        case ExpressionProgram::Op::Abs: status = centeredAbsolute(input.base, input.delta, input.endpoint, delta); break;
        case ExpressionProgram::Op::Norm: delta = centeredNormDelta(input.base, input.delta); break;
        case ExpressionProgram::Op::Arg:
            status = centeredLogarithm(input.base, input.delta, input.endpoint, delta);
            if (status == ExpressionCenteredStatus::Success) delta = {delta.imag(), 0.0};
            break;
        case ExpressionProgram::Op::Conjugate: delta = std::conj(input.delta); break;
        case ExpressionProgram::Op::Real: delta = {input.delta.real(), 0.0}; break;
        case ExpressionProgram::Op::Imaginary: delta = {input.delta.imag(), 0.0}; break;
        default:
            status = ExpressionCenteredStatus::Unsupported;
            delta = nanComplex();
            break;
        }
        stack[top - 1] = finish(base, delta, endpoint, status);
    };
    auto binary = [&](const ExpressionProgram::Instruction& instruction) {
        CenteredValue right = stack[--top];
        CenteredValue left = stack[top - 1];
        Complex base = ExpressionProgram::evaluateBinary(instruction.op, left.base, right.base);
        Complex endpoint = ExpressionProgram::evaluateBinary(instruction.op, left.endpoint, right.endpoint);
        ExpressionCenteredStatus status = combine(left.status, right.status);
        if (status != ExpressionCenteredStatus::Success) {
            stack[top - 1] = finish(base, nanComplex(), endpoint, status);
            return;
        }

        Complex delta;
        switch (instruction.op) {
        case ExpressionProgram::Op::Add: delta = left.delta + right.delta; break;
        case ExpressionProgram::Op::Subtract: delta = left.delta - right.delta; break;
        case ExpressionProgram::Op::Multiply: delta = centeredMultiply(left.base, left.delta, right.base, right.delta); break;
        case ExpressionProgram::Op::Divide: {
            Complex quotientDelta;
            status = centeredDivisionDelta(left.delta, right.base, right.delta, right.endpoint, base, quotientDelta, true);
            delta = quotientDelta;
            break;
        }
        case ExpressionProgram::Op::Power: status = centeredPower(left.base, left.delta, left.endpoint, right.base, right.delta, right.endpoint, base, endpoint, delta); break;
        case ExpressionProgram::Op::MakeComplex: delta = {left.delta.real(), right.delta.real()}; break;
        case ExpressionProgram::Op::Polar: status = centeredPolar(left, right, base, endpoint, delta); break;
        default:
            status = ExpressionCenteredStatus::Unsupported;
            delta = nanComplex();
            break;
        }
        stack[top - 1] = finish(base, delta, endpoint, status);
    };

    for (const ExpressionProgram::Instruction& instruction : program._code) {
        switch (instruction.op) {
        case ExpressionProgram::Op::Constant: pushUnchanged(instruction.value); break;
        case ExpressionProgram::Op::Z: pushPerturbed(reference.z, inputDelta.z); break;
        case ExpressionProgram::Op::C: pushPerturbed(reference.c, inputDelta.c); break;
        case ExpressionProgram::Op::Z0: pushPerturbed(reference.z0, inputDelta.z0); break;
        case ExpressionProgram::Op::Iteration: pushUnchanged({(double)reference.iteration, 0.0}); break;
        case ExpressionProgram::Op::Parameter: pushPerturbed(reference.parameters[instruction.argument], inputDelta.parameters[instruction.argument]); break;
        case ExpressionProgram::Op::OrbitInvariant: stack[top++] = finish(nanComplex(), nanComplex(), nanComplex(), ExpressionCenteredStatus::Unsupported); break;
        case ExpressionProgram::Op::Negate:
        case ExpressionProgram::Op::Square:
        case ExpressionProgram::Op::Sin:
        case ExpressionProgram::Op::Cos:
        case ExpressionProgram::Op::Tan:
        case ExpressionProgram::Op::Sinh:
        case ExpressionProgram::Op::Cosh:
        case ExpressionProgram::Op::Tanh:
        case ExpressionProgram::Op::Exp:
        case ExpressionProgram::Op::Log:
        case ExpressionProgram::Op::Log10:
        case ExpressionProgram::Op::Sqrt:
        case ExpressionProgram::Op::Abs:
        case ExpressionProgram::Op::Norm:
        case ExpressionProgram::Op::Arg:
        case ExpressionProgram::Op::Conjugate:
        case ExpressionProgram::Op::Real:
        case ExpressionProgram::Op::Imaginary: unary(instruction); break;
        case ExpressionProgram::Op::Add:
        case ExpressionProgram::Op::Subtract:
        case ExpressionProgram::Op::Multiply:
        case ExpressionProgram::Op::Divide:
        case ExpressionProgram::Op::Power:
        case ExpressionProgram::Op::MakeComplex:
        case ExpressionProgram::Op::Polar: binary(instruction); break;
        }
    }
    return {stack[0].base, stack[0].delta, stack[0].status};
}

} // namespace formula

#ifdef _MSC_VER
#pragma float_control(pop)
#endif
