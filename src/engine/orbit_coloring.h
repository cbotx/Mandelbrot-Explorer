#ifndef MANDEL_ORBIT_COLORING_H
#define MANDEL_ORBIT_COLORING_H

#include <algorithm>
#include <cmath>
#include <cstdlib>

namespace orbitcolor {

// A finite triangular tail window avoids both the hard edge of a rectangular
// window and the reference leakage of infinite-support weighting. A nonpositive
// setting preserves the classic full-orbit average.
inline int stripeWindow() {
    static int window = [] {
        const char* value = std::getenv("MANDEL_SACWIN");
        return value ? std::atoi(value) : 256;
    }();
    return window;
}

struct SacAccum {
    static constexpr int MAXW = 1024;
    double ring[MAXW];
    double RS = 0.0, TS = 0.0, full = 0.0, last = 0.0;
    int W = 0, pos = 0, fill = 0, cnt = 0;

    void init(int window) {
        W = window > MAXW ? MAXW : (window < 0 ? 0 : window);
        RS = TS = full = last = 0.0;
        pos = fill = cnt = 0;
    }

    inline void push(double zr, double zi) {
        double term = 0.5 + 0.5 * std::sin(7.0 * std::atan2(zi, zr));
        full += term;
        last = term;
        ++cnt;
        if (W > 0) {
            double evicted = fill == W ? ring[pos] : (++fill, 0.0);
            TS += (double)W * term - RS;
            RS += term - evicted;
            ring[pos] = term;
            pos = (pos + 1) % W;
        }
    }

    inline void reset_window() {
        RS = TS = 0.0;
        pos = fill = 0;
    }

    inline void add_full(double stripeSum, int count) {
        full += stripeSum;
        cnt += count;
    }

    inline void weighted(double& sum, double& weight) const {
        if (W > 0) {
            sum = TS;
            weight =
                (double)fill * W - (double)fill * (fill - 1) / 2.0;
        } else {
            sum = full;
            weight = cnt;
        }
    }

    inline float exactValue() const {
        double sum, weight;
        weighted(sum, weight);
        return (float)(weight > 0.0 ? sum / weight : 0.0);
    }

    inline float interpolatedValue(double fraction) const {
        double sum, weight;
        weighted(sum, weight);
        double including = weight > 0.0 ? sum / weight : 0.0;
        double lastWeight = W > 0 ? (double)W : 1.0;
        double excluding = weight > lastWeight
            ? (sum - lastWeight * last) / (weight - lastWeight)
            : including;
        fraction -= std::floor(fraction);
        return (float)(excluding + (including - excluding) * fraction);
    }

    inline float value(double zrad, double radius) const {
        double fraction =
            1.0 - std::log(std::log(zrad) * 0.5 / std::log(radius)) /
                      std::log(2.0);
        return interpolatedValue(fraction);
    }

    inline float powerValue(double magnitude, double radius, int power) const {
        if (power < 2 || !(radius > 1.0) || !(magnitude > 1.0) ||
            !std::isfinite(magnitude))
            return exactValue();
        double fraction =
            1.0 - std::log(std::log(magnitude) / std::log(radius)) /
                      std::log((double)power);
        return std::isfinite(fraction)
            ? interpolatedValue(fraction) : exactValue();
    }
};

struct TrapAccum {
    double minPoint = 1e300;
    double trapAngle = 0.0;
    double minCross = 1e300;
    double minCircle = 1e300;

    inline void push(double zr, double zi) {
        double magnitudeSquared = zr * zr + zi * zi;
        if (magnitudeSquared < minPoint) {
            minPoint = magnitudeSquared;
            trapAngle = std::atan2(zi, zr);
        }
        double cross = std::min(std::fabs(zr), std::fabs(zi));
        if (cross < minCross) minCross = cross;
        double circle = std::fabs(std::sqrt(magnitudeSquared) - 0.5);
        if (circle < minCircle) minCircle = circle;
    }

    inline float value(double mu) const {
        double distance = std::min(
            std::min(std::sqrt(minPoint), minCross * 1.5), minCircle);
        if (distance > 1.0) distance = 1.0;
        double trap =
            -std::log10(std::max(distance, 1e-14));
        if (mu < 0.0) mu = 0.0;
        // Distance plus escape count is seamless across atan2's branch cut.
        return (float)(0.17 * trap + 0.025 * mu);
    }
};

enum class FormulaColorMode {
    Feather,
    OrbitTrap
};

struct FormulaColorAccum {
    FormulaColorMode mode = FormulaColorMode::Feather;
    SacAccum stripe;
    TrapAccum trap;

    void init(FormulaColorMode value) {
        mode = value;
        stripe.init(stripeWindow());
        trap = {};
    }

    inline void push(double re, double im) {
        if (!std::isfinite(re) || !std::isfinite(im)) return;
        if (mode == FormulaColorMode::Feather)
            stripe.push(re, im);
        else
            trap.push(re, im);
    }

    inline float escaped(int iteration, double magnitude,
                         int power, double bailout) const {
        double result;
        if (mode == FormulaColorMode::Feather) {
            result = power > 1
                ? stripe.powerValue(magnitude, bailout, power)
                : stripe.exactValue();
        } else {
            double mu = (double)iteration;
            if (power > 1 && magnitude > 1.0 &&
                std::isfinite(magnitude)) {
                double smooth = iteration -
                    std::log(std::log(magnitude)) /
                    std::log((double)power);
                if (std::isfinite(smooth)) mu = smooth;
            }
            result = trap.value(mu);
        }
        return std::isfinite(result) && result >= 0.0
            ? (float)result : 0.0f;
    }

    inline float interior() const {
        if (mode == FormulaColorMode::Feather) return -2.0f;
        double result = trap.value(0.0);
        return std::isfinite(result) && result >= 0.0
            ? (float)result : 0.0f;
    }
};

} // namespace orbitcolor

#endif
