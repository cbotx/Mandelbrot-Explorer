#ifndef MANDEL_CUSTOM_DEEP_ZOOM_H
#define MANDEL_CUSTOM_DEEP_ZOOM_H

#include <algorithm>
#include <cmath>
#include <limits>

#include <gmp.h>

#include "formula_expression.h"
#include "formula_reference_orbit.h"

namespace formula {

constexpr double CUSTOM_DIRECT_ZOOM_LIMIT = 1.0e12;
constexpr double CUSTOM_DEEP_MIN_ESCAPE_RADIUS = 2.0;

enum class CustomDeepZoomPath : uint8_t {
    DirectExpression,
    QuadraticPerturbation
};

enum class ExpressionProductionPath : uint8_t {
    DirectExpression,
    QuadraticPerturbation,
    GenericCertifiedDeep,
    UnsupportedDeep
};

enum class CustomDeepZoomOutputAdapter : uint8_t {
    None,
    SmoothExpression,
    DistanceExpression,
    FeatherExpression,
    OrbitTrapExpression
};

struct CustomDeepZoomPlan {
    CustomDeepZoomPath path = CustomDeepZoomPath::DirectExpression;
    CustomDeepZoomOutputAdapter outputAdapter =
        CustomDeepZoomOutputAdapter::None;
    double escapeRadius = 0.0;
    bool canonicalQuadratic = false;
    bool cPlane = false;
    bool standardInitialZ = false;
    bool safeEscapeRadius = false;
    bool coloringCompatible = false;
    bool initialIterationCompatible = false;

    bool recurrenceCompatible() const {
        return canonicalQuadratic && cPlane && standardInitialZ &&
               safeEscapeRadius;
    }

    bool canZoomBeyondDirectLimit() const {
        return recurrenceCompatible() && coloringCompatible &&
               initialIterationCompatible;
    }

    bool usesQuadraticPerturbation() const {
        return path == CustomDeepZoomPath::QuadraticPerturbation;
    }
};

struct ExpressionProductionPlan {
    ExpressionProductionPath path =
        ExpressionProductionPath::DirectExpression;
    CustomDeepZoomPlan quadratic;
    bool aboveDirectLimit = false;
    bool rawOutput = false;

    bool canRender() const {
        return path != ExpressionProductionPath::UnsupportedDeep;
    }

    bool canZoomBeyondDirectLimit() const {
        return path == ExpressionProductionPath::QuadraticPerturbation ||
               path == ExpressionProductionPath::GenericCertifiedDeep;
    }

    bool usesQuadraticPerturbation() const {
        return path == ExpressionProductionPath::QuadraticPerturbation;
    }

    bool usesGenericCertifiedDeep() const {
        return path == ExpressionProductionPath::GenericCertifiedDeep;
    }
};

inline bool genericDeepPrecisionSupported(
        mpf_srcptr scale, mpf_srcptr centerRe,
        mpf_srcptr centerIm, int width, int height) {
    if (!scale || !centerRe || !centerIm ||
        width < 2 || height < 2 || mpf_sgn(scale) <= 0)
        return false;
    uint64_t required = std::max<uint64_t>({
        mpf_get_prec(scale),
        mpf_get_prec(centerRe),
        mpf_get_prec(centerIm),
        128
    });
    auto exponentOf = [](mpf_srcptr value) {
        if (mpf_sgn(value) == 0) return 0L;
        long exponent = 0;
        mpf_get_d_2exp(&exponent, value);
        return exponent;
    };
    const long scaleExponent = exponentOf(scale);
    const long centerExponent = std::max(
        exponentOf(centerRe), exponentOf(centerIm));
    uint64_t dimensionBits = 0;
    uint64_t dimension = static_cast<uint64_t>(
        std::max(width, height));
    while ((uint64_t{1} << dimensionBits) < dimension &&
           dimensionBits < 63)
        ++dimensionBits;
    if (scaleExponent > 0) {
        const uint64_t scaleBits =
            static_cast<uint64_t>(scaleExponent);
        const uint64_t centerBits = centerExponent > 0
            ? static_cast<uint64_t>(centerExponent) : 0;
        if (scaleBits >
                std::numeric_limits<uint64_t>::max() -
                    centerBits - dimensionBits - 64)
            return false;
        required = std::max(
            required,
            scaleBits + centerBits + dimensionBits + 8);
    }
    constexpr uint64_t rendererGuardBits = 64;
    constexpr uint64_t fallbackGuardBits = 128;
    constexpr uint64_t reservedBits =
        rendererGuardBits + fallbackGuardBits;
    const uint64_t maximum =
        static_cast<uint64_t>(
            ExpressionReferencePrecisionPolicy::
                ApplicationMaximumBits);
    return required <= maximum - reservedBits;
}

inline bool isStandardMandelbrotInitialZ(Complex z0) {
    return z0.real() == 0.0 && z0.imag() == 0.0 &&
           !std::signbit(z0.real()) && !std::signbit(z0.imag());
}

inline double customDeepMaxEscapeRadius() {
    // The escaped value can be approximately radius^2 and coloring then forms
    // its squared magnitude. Keep enough headroom for (radius^2 + radius)^2.
    static const double maxRadius =
        0.5 * std::sqrt(std::sqrt(
            std::numeric_limits<double>::max()));
    return maxRadius;
}

inline bool isCustomDeepEscapeRadiusSupported(double radius) {
    return std::isfinite(radius) &&
           radius >= CUSTOM_DEEP_MIN_ESCAPE_RADIUS &&
           radius <= customDeepMaxEscapeRadius();
}

inline CustomDeepZoomPlan makeCustomDeepZoomPlan(
        const ExpressionProgram& sourceProgram,
        const ExpressionProgram& runtimeProgram,
        const ExpressionContext& fixed,
        FormulaParameter pixelParameter,
        double bailout,
        ExpressionColoring coloring,
        mpf_srcptr scale = nullptr,
        mpf_srcptr centerRe = nullptr,
        mpf_srcptr centerIm = nullptr,
        int width = 0,
        int height = 0) {
    CustomDeepZoomPlan plan;
    plan.escapeRadius = bailout;
    plan.canonicalQuadratic =
        sourceProgram.isCanonicalQuadraticPlusC() &&
        runtimeProgram.isCanonicalQuadraticPlusC();
    plan.cPlane = pixelParameter == FormulaParameter::C;
    plan.standardInitialZ = isStandardMandelbrotInitialZ(fixed.z0);
    plan.safeEscapeRadius =
        isCustomDeepEscapeRadiusSupported(bailout);
    if (!scale ||
        mpf_cmp_d(scale, CUSTOM_DIRECT_ZOOM_LIMIT) <= 0) {
        plan.initialIterationCompatible = true;
    } else if (centerRe && centerIm && width >= 2 && height >= 2) {
        const double re = mpf_get_d(centerRe);
        const double im = mpf_get_d(centerIm);
        const double zoom = mpf_get_d(scale);
        if (std::isfinite(re) && std::isfinite(im) && zoom > 0.0) {
            const double halfWidth = std::isfinite(zoom)
                ? 2.0 / zoom : 0.0;
            const double halfHeight =
                halfWidth * (double)height / (double)width;
            const double farthest = std::hypot(
                std::fabs(re) + halfWidth,
                std::fabs(im) + halfHeight);
            // The production engine starts its optimized loop at z1=c and does
            // not test z1. Routing is therefore exact only when every pixel's
            // z1 is proven inside the bailout disk.
            const double guard = std::max(
                1.0e-9, bailout * 1.0e-12);
            plan.initialIterationCompatible =
                farthest <= bailout - guard;
        }
    }
    if (coloring == ExpressionColoring::Smooth) {
        // The direct expression renderer uses n-log(log|z|)/log(2), while the
        // production quadratic engine's historical convention subtracts
        // log(log|z|/log(2))/log(2). The deep writer applies the constant offset.
        plan.coloringCompatible = true;
        plan.outputAdapter = CustomDeepZoomOutputAdapter::SmoothExpression;
    } else if (coloring == ExpressionColoring::Distance) {
        // The quadratic engine uses log(|z|^2); the expression renderer uses
        // log(|z|), so its pixel-normalised distance is exactly one half.
        plan.coloringCompatible = true;
        plan.outputAdapter = CustomDeepZoomOutputAdapter::DistanceExpression;
    } else if (coloring == ExpressionColoring::Feather) {
        plan.coloringCompatible = true;
        plan.outputAdapter = CustomDeepZoomOutputAdapter::FeatherExpression;
    } else if (coloring == ExpressionColoring::OrbitTrap) {
        plan.coloringCompatible = true;
        plan.outputAdapter =
            CustomDeepZoomOutputAdapter::OrbitTrapExpression;
    }
    // Raw has a different iteration convention and remains direct/capped.
    if (plan.canZoomBeyondDirectLimit() && scale &&
        mpf_cmp_d(scale, CUSTOM_DIRECT_ZOOM_LIMIT) > 0) {
        plan.path = CustomDeepZoomPath::QuadraticPerturbation;
    }
    return plan;
}

inline ExpressionProductionPlan makeExpressionProductionPlan(
        const ExpressionProgram& sourceProgram,
        const ExpressionProgram& runtimeProgram,
        const ExpressionContext& fixed,
        FormulaParameter pixelParameter,
        double bailout,
        ExpressionColoring coloring,
        mpf_srcptr scale,
        mpf_srcptr centerRe,
        mpf_srcptr centerIm,
        int width,
        int height) {
    ExpressionProductionPlan plan;
    plan.quadratic = makeCustomDeepZoomPlan(
        sourceProgram, runtimeProgram, fixed, pixelParameter,
        bailout, coloring, scale, centerRe, centerIm,
        width, height);
    plan.rawOutput = coloring == ExpressionColoring::Raw;
    plan.aboveDirectLimit =
        scale && mpf_cmp_d(scale, CUSTOM_DIRECT_ZOOM_LIMIT) > 0;
    if (!plan.aboveDirectLimit) return plan;
    if (plan.quadratic.usesQuadraticPerturbation()) {
        plan.path = ExpressionProductionPath::QuadraticPerturbation;
    } else if (plan.rawOutput &&
               genericDeepPrecisionSupported(
                   scale, centerRe, centerIm,
                   width, height)) {
        plan.path = ExpressionProductionPath::GenericCertifiedDeep;
    } else {
        plan.path = ExpressionProductionPath::UnsupportedDeep;
    }
    return plan;
}

} // namespace formula

#endif
