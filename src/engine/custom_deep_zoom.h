#ifndef MANDEL_CUSTOM_DEEP_ZOOM_H
#define MANDEL_CUSTOM_DEEP_ZOOM_H

#include <algorithm>
#include <cmath>
#include <limits>

#include <gmp.h>

#include "formula_expression.h"

namespace formula {

constexpr double CUSTOM_DIRECT_ZOOM_LIMIT = 1.0e12;
constexpr double CUSTOM_DEEP_MIN_ESCAPE_RADIUS = 2.0;

enum class CustomDeepZoomPath : uint8_t {
    DirectExpression,
    QuadraticPerturbation
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

} // namespace formula

#endif
