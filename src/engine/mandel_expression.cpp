#include "mandel_perturbation.h"

#include <algorithm>
#include <array>
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <immintrin.h>
#include <limits>

#include "float_math.h"
#if defined(MANDEL_ENABLE_ASMJIT)
#include "formula_expression_jit.h"
#endif

namespace {

constexpr int CUBIC_SA_TERMS = 30;

struct CubicSeries {
    std::array<formula::Complex, CUBIC_SA_TERMS> terms{};
    formula::Complex scale{};
    int iteration = 0;
    int order = 0;
};

int cubicSeriesOrder(
        const std::array<formula::Complex, CUBIC_SA_TERMS>& terms) {
    double prefix = 0.0;
    std::array<double, CUBIC_SA_TERMS> suffix{};
    for (int i = CUBIC_SA_TERMS - 1; i >= 0; --i) {
        double magnitude = std::abs(terms[(size_t)i]);
        if (!std::isfinite(magnitude)) return -1;
        suffix[(size_t)i] = i + 1 < CUBIC_SA_TERMS
            ? std::max(magnitude, suffix[(size_t)i + 1])
            : magnitude;
    }
    for (int i = 0; i < CUBIC_SA_TERMS - 10; ++i) {
        prefix = std::max(prefix, std::abs(terms[(size_t)i]));
        double tolerance =
            std::max(prefix, std::numeric_limits<double>::min()) *
            std::ldexp(1.0, -80);
        if (suffix[(size_t)i + 1] <= tolerance)
            return std::max(i, 5);
    }
    return -1;
}

CubicSeries buildCubicSeries(
        const std::vector<formula::Complex>& orbit,
        FormulaParameter pixelParameter, formula::Complex scale,
        int mxit, double bailout, volatile bool* halt) {
    CubicSeries best;
    best.scale = scale;
    if (scale == formula::Complex{} ||
        !std::isfinite(scale.real()) ||
        !std::isfinite(scale.imag()))
        return best;

    std::array<formula::Complex, CUBIC_SA_TERMS> old{};
    std::array<formula::Complex, CUBIC_SA_TERMS> next{};
    std::array<formula::Complex, CUBIC_SA_TERMS - 1> square{};
    if (pixelParameter == FormulaParameter::InitialZ)
        old[0] = scale;

    int limit = std::min(mxit, (int)orbit.size() - 1);
    for (int iteration = 0; iteration < limit; ++iteration) {
        if (*halt) break;
        square.fill({});
        for (int a = 0; a < CUBIC_SA_TERMS; ++a) {
            for (int b = 0;
                 a + b < CUBIC_SA_TERMS - 1; ++b) {
                square[(size_t)(a + b)] +=
                    old[(size_t)a] * old[(size_t)b];
            }
        }

        formula::Complex z = orbit[(size_t)iteration];
        formula::Complex zSquared = z * z;
        for (int coefficient = 0;
             coefficient < CUBIC_SA_TERMS; ++coefficient) {
            formula::Complex value =
                (3.0 * zSquared) * old[(size_t)coefficient];
            if (coefficient >= 1) {
                value += (3.0 * z) *
                    square[(size_t)(coefficient - 1)];
            }
            if (coefficient >= 2) {
                formula::Complex cubic{};
                for (int squareIndex = 0;
                     squareIndex <= coefficient - 2; ++squareIndex) {
                    cubic += square[(size_t)squareIndex] *
                        old[(size_t)(coefficient - 2 - squareIndex)];
                }
                value += cubic;
            }
            next[(size_t)coefficient] = value;
        }
        if (pixelParameter == FormulaParameter::C)
            next[0] += scale;

        int order = cubicSeriesOrder(next);
        if (order < 0) break;
        double deltaBound = 0.0;
        for (int coefficient = 0; coefficient <= order; ++coefficient)
            deltaBound += std::abs(next[(size_t)coefficient]);
        double absoluteBound =
            std::abs(orbit[(size_t)iteration + 1]) + deltaBound;
        if (!std::isfinite(absoluteBound) ||
            absoluteBound > bailout)
            break;
        best.terms = next;
        best.iteration = iteration + 1;
        best.order = order;
        old = next;
    }
    return best;
}

void evaluateCubicSeries(
        const CubicSeries& series, formula::Complex pixelDelta,
        formula::Complex& delta, formula::Complex& derivative) {
    delta = {};
    derivative = {};
    if (series.iteration <= 0 || series.scale == formula::Complex{})
        return;
    formula::Complex q = pixelDelta / series.scale;
    for (int coefficient = series.order; coefficient >= 0; --coefficient) {
        delta += series.terms[(size_t)coefficient];
        delta *= q;
        derivative *= q;
        derivative +=
            (double)(coefficient + 1) *
            series.terms[(size_t)coefficient] / series.scale;
    }
}

float initialPowerResult(formula::Complex z, FormulaParameter pixelParameter,
                         formula::ExpressionColoring coloring,
                         int power, double bailout, double dx) {
    double magnitude = std::hypot(z.real(), z.imag());
    if (std::isfinite(z.real()) && std::isfinite(z.imag()) &&
        magnitude <= bailout)
        return EMPTYPIXEL;
    if (coloring == formula::ExpressionColoring::Distance &&
        pixelParameter == FormulaParameter::InitialZ) {
        return (float)(magnitude * std::log(magnitude) / std::fabs(dx));
    }
    if (coloring == formula::ExpressionColoring::Smooth &&
        std::isfinite(magnitude) && magnitude > 1.0) {
        return (float)(-std::log(std::log(magnitude)) /
                       std::log((double)power));
    }
    return 0.0f;
}

bool solveIntegerPowerSimdRow(
        double startRe, double dx, double pixelIm, int count,
        int power, FormulaParameter pixelParameter,
        const formula::ExpressionContext& fixed, int mxit,
        double bailout, formula::ExpressionColoring coloring,
        float* output, volatile bool* halt) {
    const bool distance = coloring == formula::ExpressionColoring::Distance;
    const double bailoutSquared = bailout * bailout;
    const double logPower = std::log((double)power);
    const __m256d one = _mm256_set1_pd(1.0);
    const __m256d powerVector = _mm256_set1_pd((double)power);
    const __m256d bailoutVector = _mm256_set1_pd(bailoutSquared);
    const __m256d maxDouble = _mm256_set1_pd(DBL_MAX);
    const __m256d signMask = _mm256_set1_pd(-0.0);
    const __m256d mxitVector = _mm256_set1_pd((double)mxit);

    alignas(32) double cr_[4], ci_[4], zr_[4], zi_[4];
    alignas(32) double dr_[4], di_[4], iteration_[4];
    int lanePixel[4] = { -1, -1, -1, -1 };
    int nextPixel = 0;
    int activeCount = 0;

    auto loadLane = [&](int lane) {
        lanePixel[lane] = -1;
        while (nextPixel < count) {
            int pixelIndex = nextPixel++;
            formula::Complex pixel{
                startRe + dx * pixelIndex, pixelIm
            };
            formula::Complex z = pixelParameter == FormulaParameter::InitialZ
                ? pixel : fixed.z0;
            formula::Complex c = pixelParameter == FormulaParameter::C
                ? pixel : fixed.c;
            float initial = initialPowerResult(
                z, pixelParameter, coloring, power, bailout, dx);
            if (initial != EMPTYPIXEL) {
                output[pixelIndex] = initial;
                continue;
            }
            cr_[lane] = c.real();
            ci_[lane] = c.imag();
            zr_[lane] = z.real();
            zi_[lane] = z.imag();
            dr_[lane] =
                pixelParameter == FormulaParameter::InitialZ ? 1.0 : 0.0;
            di_[lane] = 0.0;
            iteration_[lane] = 0.0;
            lanePixel[lane] = pixelIndex;
            ++activeCount;
            return;
        }
        cr_[lane] = ci_[lane] = zr_[lane] = zi_[lane] = 0.0;
        dr_[lane] = di_[lane] = iteration_[lane] = 0.0;
    };

    for (int lane = 0; lane < 4; ++lane) loadLane(lane);

    __m256d cr = _mm256_load_pd(cr_);
    __m256d ci = _mm256_load_pd(ci_);
    __m256d zr = _mm256_load_pd(zr_);
    __m256d zi = _mm256_load_pd(zi_);
    __m256d dr = _mm256_load_pd(dr_);
    __m256d di = _mm256_load_pd(di_);
    __m256d iteration = _mm256_load_pd(iteration_);
    unsigned steps = 0;

    while (activeCount > 0) {
        if ((steps++ & 255u) == 0u && *halt) return false;

        // Match the scalar std::complex multiplication order: form z^(d-1)
        // by repeated complex products, then multiply once more for z^d.
        __m256d powerMinusOneRe = zr;
        __m256d powerMinusOneIm = zi;
        for (int exponent = 2; exponent < power; ++exponent) {
            __m256d nextRe = _mm256_sub_pd(
                _mm256_mul_pd(powerMinusOneRe, zr),
                _mm256_mul_pd(powerMinusOneIm, zi));
            __m256d nextIm = _mm256_add_pd(
                _mm256_mul_pd(powerMinusOneRe, zi),
                _mm256_mul_pd(powerMinusOneIm, zr));
            powerMinusOneRe = nextRe;
            powerMinusOneIm = nextIm;
        }
        __m256d nextZr = _mm256_add_pd(
            _mm256_sub_pd(
                _mm256_mul_pd(powerMinusOneRe, zr),
                _mm256_mul_pd(powerMinusOneIm, zi)),
            cr);
        __m256d nextZi = _mm256_add_pd(
            _mm256_add_pd(
                _mm256_mul_pd(powerMinusOneRe, zi),
                _mm256_mul_pd(powerMinusOneIm, zr)),
            ci);

        __m256d nextDr = dr;
        __m256d nextDi = di;
        if (distance) {
            nextDr = _mm256_mul_pd(
                powerVector,
                _mm256_sub_pd(
                    _mm256_mul_pd(powerMinusOneRe, dr),
                    _mm256_mul_pd(powerMinusOneIm, di)));
            nextDi = _mm256_mul_pd(
                powerVector,
                _mm256_add_pd(
                    _mm256_mul_pd(powerMinusOneRe, di),
                    _mm256_mul_pd(powerMinusOneIm, dr)));
            if (pixelParameter == FormulaParameter::C)
                nextDr = _mm256_add_pd(nextDr, one);
        }
        zr = nextZr;
        zi = nextZi;
        dr = nextDr;
        di = nextDi;

        __m256d magnitudeSquared = _mm256_add_pd(
            _mm256_mul_pd(zr, zr), _mm256_mul_pd(zi, zi));
        __m256d nextIteration = _mm256_add_pd(iteration, one);
        int escaped = _mm256_movemask_pd(_mm256_cmp_pd(
            magnitudeSquared, bailoutVector, _CMP_GT_OQ));
        __m256d absZr = _mm256_andnot_pd(signMask, zr);
        __m256d absZi = _mm256_andnot_pd(signMask, zi);
        int nonfinite = _mm256_movemask_pd(_mm256_or_pd(
            _mm256_or_pd(
                _mm256_cmp_pd(zr, zr, _CMP_UNORD_Q),
                _mm256_cmp_pd(zi, zi, _CMP_UNORD_Q)),
            _mm256_or_pd(
                _mm256_cmp_pd(absZr, maxDouble, _CMP_GT_OQ),
                _mm256_cmp_pd(absZi, maxDouble, _CMP_GT_OQ))));
        int reachedMaximum = _mm256_movemask_pd(_mm256_cmp_pd(
            nextIteration, mxitVector, _CMP_GE_OQ));
        int activeMask =
            (lanePixel[0] >= 0) |
            ((lanePixel[1] >= 0) << 1) |
            ((lanePixel[2] >= 0) << 2) |
            ((lanePixel[3] >= 0) << 3);
        int finishMask =
            (escaped | nonfinite | reachedMaximum) & activeMask;
        if (finishMask == 0) {
            iteration = nextIteration;
            continue;
        }

        _mm256_store_pd(cr_, cr);
        _mm256_store_pd(ci_, ci);
        _mm256_store_pd(zr_, zr);
        _mm256_store_pd(zi_, zi);
        _mm256_store_pd(dr_, dr);
        _mm256_store_pd(di_, di);
        _mm256_store_pd(iteration_, iteration);
        for (int lane = 0; lane < 4; ++lane) {
            if (lanePixel[lane] < 0) continue;
            int bit = 1 << lane;
            if ((finishMask & bit) == 0) {
                iteration_[lane] += 1.0;
                continue;
            }

            float result = -2.0f;
            if ((escaped | nonfinite) & bit) {
                double magnitude = std::hypot(zr_[lane], zi_[lane]);
                if (coloring == formula::ExpressionColoring::Smooth &&
                    std::isfinite(magnitude) && magnitude > 1.0) {
                    result = (float)(iteration_[lane] + 1.0 -
                        std::log(std::log(magnitude)) / logPower);
                } else if (distance && std::isfinite(magnitude)) {
                    double derivative =
                        std::hypot(dr_[lane], di_[lane]);
                    double denominator =
                        derivative * std::fabs(dx);
                    result = denominator > 0.0 &&
                             std::isfinite(denominator)
                        ? (float)(magnitude * std::log(magnitude) /
                                  denominator)
                        : 0.0f;
                } else {
                    result = (float)(iteration_[lane] + 1.0);
                }
            }
            output[lanePixel[lane]] = result;
            --activeCount;
            loadLane(lane);
        }

        cr = _mm256_load_pd(cr_);
        ci = _mm256_load_pd(ci_);
        zr = _mm256_load_pd(zr_);
        zi = _mm256_load_pd(zi_);
        dr = _mm256_load_pd(dr_);
        di = _mm256_load_pd(di_);
        iteration = _mm256_load_pd(iteration_);
    }
    return true;
}

} // namespace

bool Mandel::ComputeExpression(mpf_t center_re, mpf_t center_im, mpf_t scale,
                               const formula::ExpressionProgram& program,
                               const formula::ExpressionContext& fixed,
                               FormulaParameter pixelParameter,
                               int mxit, double bailout,
                               formula::ExpressionColoring coloring,
                               const formula::ExpressionJit4* jit) {
    if (_sub != 1 || !program.valid() || mxit < 1 || !(bailout > 0.0) ||
        !std::isfinite(bailout) ||
        (pixelParameter != FormulaParameter::C &&
         pixelParameter != FormulaParameter::InitialZ))
        return false;

    if (_flag_halt) return false;
    std::fill(_iter, _iter + (size_t)_w * _h, EMPTYPIXEL);
    mpf_set(_scale, scale);

    mpf_t dw, dh;
    mpf_init_set_ui(dw, 2);
    mpf_div(dw, dw, scale);
    mpf_init_set(dh, dw);
    mpf_mul_ui(dh, dh, _h);
    mpf_div_ui(dh, dh, _w);
    mpf_sub(_c0_re, center_re, dw);
    mpf_sub(_c0_im, center_im, dh);
    mpf_mul_ui(_dx, dw, 2); mpf_div_ui(_dx, _dx, _w - 1);
    mpf_mul_ui(_dy, dh, 2); mpf_div_ui(_dy, _dy, _h - 1);
    mpf_clear(dw); mpf_clear(dh);

    const double startRe = mpf_get_ld(_c0_re), startIm = mpf_get_ld(_c0_im);
    const double dx = mpf_get_ld(_dx), dy = mpf_get_ld(_dy);
    const int integerPower =
        program.fastPath() == formula::ExpressionProgram::FastPath::IntegerPowerPlusC
            ? program.fastIntegerPower() : 0;
    const char* powerSimdSetting = std::getenv("MANDEL_EXPR_POWER_SIMD");
    const char* vectorSetting = std::getenv("MANDEL_EXPR_VECTOR");
    const bool vectorEnabled =
        !vectorSetting || std::atoi(vectorSetting) != 0;
    const double bailoutSquared = bailout * bailout;
    const bool powerSimd =
        integerPower >= 2 && integerPower <= 8 &&
        bailoutSquared >= DBL_MIN && std::isfinite(bailoutSquared) &&
        (!powerSimdSetting || std::atoi(powerSimdSetting) != 0);
    if (integerPower < 2 || bailout < 1.0)
        coloring = formula::ExpressionColoring::Raw;
    // Reserve the final progress slot for the successful completion commit.
    // Completed rows alone can therefore never publish exactly 100%.
    progressBegin(_h + 1, 0.0, 1.0);
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < _h; ++i) {
        if (_flag_halt) continue;
        bool rowCompleted = true;
        const double pixelIm = startIm + dy * i;
        if (powerSimd) {
            rowCompleted = solveIntegerPowerSimdRow(
                startRe, dx, pixelIm, _w, integerPower,
                pixelParameter, fixed, mxit, bailout, coloring,
                _iter + (size_t)i * _w, &_flag_halt);
            if (rowCompleted) progressAdvance();
            continue;
        }
        if (integerPower == 0 && program.avx2Compatible() &&
            vectorEnabled) {
            for (int j = 0; j < _w; j += 4) {
                if (_flag_halt) { rowCompleted = false; break; }
                int lanes = std::min(4, _w - j);
                formula::ExpressionContext contexts[4] = { fixed, fixed, fixed, fixed };
                formula::Complex outputs[4]{};
                bool active[4] = { false, false, false, false };
                float results[4] = { -2.0f, -2.0f, -2.0f, -2.0f };
                int activeCount = 0;
                for (int lane = 0; lane < lanes; ++lane) {
                    formula::Complex pixel{ startRe + dx * (j + lane), pixelIm };
                    if (pixelParameter == FormulaParameter::C) contexts[lane].c = pixel;
                    else contexts[lane].z0 = pixel;
                    contexts[lane].z = contexts[lane].z0;
                    if (!std::isfinite(contexts[lane].z.real()) ||
                        !std::isfinite(contexts[lane].z.imag()) ||
                        std::hypot(contexts[lane].z.real(),
                                   contexts[lane].z.imag()) > bailout) {
                        results[lane] = 0.0f;
                    } else {
                        active[lane] = true;
                        ++activeCount;
                    }
                }
#if defined(MANDEL_ENABLE_ASMJIT)
                const bool useJit = jit && jit->valid();
                if (useJit) {
                    rowCompleted = jit->evaluateOrbit(
                        contexts, lanes, mxit, bailout,
                        results, &_flag_halt);
                    if (!rowCompleted) break;
                    for (int lane = 0; lane < lanes; ++lane)
                        _iter[i * _w + j + lane] = results[lane];
                    continue;
                }
#endif
                for (int n = 0; n < mxit && activeCount > 0; ++n) {
                    if ((n & 255) == 0 && _flag_halt) {
                        rowCompleted = false;
                        break;
                    }
                    for (int lane = 0; lane < 4; ++lane)
                        contexts[lane].iteration = n;
                    if (!program.evaluate4(contexts, outputs)) {
                        rowCompleted = false;
                        break;
                    }
                    for (int lane = 0; lane < lanes; ++lane) {
                        if (!active[lane]) continue;
                        contexts[lane].z = outputs[lane];
                        double re = outputs[lane].real(), im = outputs[lane].imag();
                        if (!std::isfinite(re) || !std::isfinite(im) ||
                            std::hypot(re, im) > bailout) {
                            results[lane] = (float)(n + 1);
                            active[lane] = false;
                            --activeCount;
                        }
                    }
                }
                if (!rowCompleted) break;
                for (int lane = 0; lane < lanes; ++lane)
                    _iter[i * _w + j + lane] = results[lane];
            }
            if (rowCompleted) progressAdvance();
            continue;
        }
        for (int j = 0; j < _w; ++j) {
            if (_flag_halt) { rowCompleted = false; break; }
            formula::ExpressionContext context = fixed;
            formula::Complex pixel{ startRe + dx * j, pixelIm };
            if (pixelParameter == FormulaParameter::C) context.c = pixel;
            else context.z0 = pixel;
            context.z = context.z0;
            std::array<formula::Complex, formula::ExpressionProgram::MAX_STACK> stack;

            float result = -2.0f;
            if (!std::isfinite(context.z.real()) || !std::isfinite(context.z.imag()) ||
                std::hypot(context.z.real(), context.z.imag()) > bailout) {
                double magnitude = std::hypot(context.z.real(), context.z.imag());
                if (coloring == formula::ExpressionColoring::Distance &&
                    pixelParameter == FormulaParameter::InitialZ) {
                    result = (float)(magnitude * std::log(magnitude) /
                                     std::fabs(dx));
                } else if (coloring == formula::ExpressionColoring::Smooth &&
                           std::isfinite(magnitude) && magnitude > 1.0) {
                    result = (float)(-std::log(std::log(magnitude)) /
                                     std::log((double)integerPower));
                } else {
                    result = 0.0f;
                }
            } else {
                if (integerPower >= 2) {
                    formula::Complex z = context.z;
                    const formula::Complex c = context.c;
                    formula::Complex derivative =
                        pixelParameter == FormulaParameter::InitialZ
                            ? formula::Complex{ 1.0, 0.0 } : formula::Complex{};
                    for (int n = 0; n < mxit; ++n) {
                        if ((n & 255) == 0 && _flag_halt) {
                            rowCompleted = false;
                            break;
                        }
                        formula::Complex powerMinusOne = z;
                        for (int power = 1; power < integerPower - 1; ++power)
                            powerMinusOne *= z;
                        formula::Complex next = powerMinusOne * z;
                        formula::Complex nextDerivative =
                            (double)integerPower * powerMinusOne * derivative;
                        if (pixelParameter == FormulaParameter::C)
                            nextDerivative += 1.0;
                        z = next + c;
                        derivative = nextDerivative;
                        bool escaped = !std::isfinite(z.real()) ||
                                       !std::isfinite(z.imag()) ||
                                       std::hypot(z.real(), z.imag()) > bailout;
                        if (escaped) {
                            double magnitude = std::hypot(z.real(), z.imag());
                            if (coloring == formula::ExpressionColoring::Smooth &&
                                std::isfinite(magnitude) && magnitude > 1.0) {
                                result = (float)(n + 1 -
                                    std::log(std::log(magnitude)) /
                                    std::log((double)integerPower));
                            } else if (coloring == formula::ExpressionColoring::Distance &&
                                       std::isfinite(magnitude)) {
                                double denominator = std::abs(derivative) *
                                                     std::fabs(dx);
                                result = denominator > 0.0 &&
                                         std::isfinite(denominator)
                                    ? (float)(magnitude * std::log(magnitude) /
                                              denominator)
                                    : 0.0f;
                            } else {
                                result = (float)(n + 1);
                            }
                            break;
                        }
                    }
                } else {
                    for (int n = 0; n < mxit; ++n) {
                        if ((n & 255) == 0 && _flag_halt) {
                            rowCompleted = false;
                            break;
                        }
                        context.iteration = n;
                        context.z = program.evaluate(context, stack.data(),
                                                     program.stackDepth());
                        double re = context.z.real(), im = context.z.imag();
                        if (!std::isfinite(re) || !std::isfinite(im) ||
                            std::hypot(re, im) > bailout) {
                            result = (float)(n + 1);
                            break;
                        }
                    }
                }
            }
            if (!rowCompleted) break;
            _iter[i * _w + j] = result;
        }
        if (rowCompleted) progressAdvance();
    }
    bool completed = !_flag_halt;
    if (completed) progressSet(1.0);
    return completed;
}

bool Mandel::ComputeExpressionResidual(
        mpf_t center_re, mpf_t center_im, mpf_t scale,
        const formula::ExpressionProgram& program,
        const formula::ExpressionContext& fixed,
        FormulaParameter pixelParameter, int mxit, double bailout,
        bool* usedPerturbation,
        formula::ExpressionColoring coloring,
        int* seriesIterations,
        int minimumSeriesIterations) {
    if (usedPerturbation) *usedPerturbation = false;
    if (seriesIterations) *seriesIterations = 0;
    if (_sub != 1 || !program.valid() || mxit < 1 || !(bailout > 0.0) ||
        !std::isfinite(bailout) ||
        (pixelParameter != FormulaParameter::C &&
         pixelParameter != FormulaParameter::InitialZ))
        return false;
    if (_flag_halt) return false;
    if (bailout < 1.0)
        coloring = formula::ExpressionColoring::Raw;

    mpf_t dw, dh;
    mpf_init_set_ui(dw, 2);
    mpf_div(dw, dw, scale);
    mpf_init_set(dh, dw);
    mpf_mul_ui(dh, dh, _h);
    mpf_div_ui(dh, dh, _w);
    mpf_sub(_c0_re, center_re, dw);
    mpf_sub(_c0_im, center_im, dh);
    mpf_mul_ui(_dx, dw, 2); mpf_div_ui(_dx, _dx, _w - 1);
    mpf_mul_ui(_dy, dh, 2); mpf_div_ui(_dy, _dy, _h - 1);
    mpf_set(_scale, scale);
    mpf_clear(dw); mpf_clear(dh);

    formula::ExpressionContext reference = fixed;
    const int integerPower =
        program.fastPath() ==
            formula::ExpressionProgram::FastPath::IntegerPowerPlusC
        ? program.fastIntegerPower() : 0;
    const char* residualPowerSetting =
        std::getenv("MANDEL_EXPR_RESIDUAL_POWER");
    const bool cubicResidual =
        integerPower == 3 &&
        (!residualPowerSetting || std::atoi(residualPowerSetting) != 0);
    if (coloring != formula::ExpressionColoring::Raw &&
        !cubicResidual)
        return false;
    formula::Complex center{ mpf_get_ld(center_re), mpf_get_ld(center_im) };
    if (pixelParameter == FormulaParameter::C) reference.c = center;
    else reference.z0 = center;
    reference.z = reference.z0;
    std::vector<formula::Complex> orbit((size_t)mxit + 1);
    orbit[0] = reference.z;
    std::array<formula::Complex, formula::ExpressionProgram::MAX_STACK> refStack;
    progressBegin(mxit + _h + 1, 0.0, 1.0);
    bool boundedReference =
        std::isfinite(reference.z.real()) && std::isfinite(reference.z.imag()) &&
        std::hypot(reference.z.real(), reference.z.imag()) <= bailout;
    for (int n = 0; n < mxit && boundedReference; ++n) {
        if (_flag_halt) return false;
        reference.iteration = n;
        reference.z = orbit[n];
        if (cubicResidual) {
            formula::Complex square = orbit[n] * orbit[n];
            orbit[n + 1] = square * orbit[n] + reference.c;
        } else {
            orbit[n + 1] = program.evaluate(
                reference, refStack.data(), program.stackDepth());
        }
        boundedReference =
            std::isfinite(orbit[n + 1].real()) &&
            std::isfinite(orbit[n + 1].imag()) &&
            std::hypot(orbit[n + 1].real(), orbit[n + 1].imag()) <= bailout;
        progressAdvance();
    }
    if (!boundedReference) {
        return ComputeExpression(center_re, center_im, scale, program, fixed,
                                 pixelParameter, mxit, bailout,
                                 coloring, nullptr);
    }
    if (usedPerturbation) *usedPerturbation = true;

    progressSet(0.0);
    std::fill(_iter, _iter + (size_t)_w * _h, EMPTYPIXEL);
    const double startRe = mpf_get_ld(_c0_re), startIm = mpf_get_ld(_c0_im);
    const double dx = mpf_get_ld(_dx), dy = mpf_get_ld(_dy);
    const double halfWidth = (_w - 1) * 0.5;
    const double halfHeight = (_h - 1) * 0.5;
    CubicSeries series;
    const char* seriesSetting = std::getenv("MANDEL_EXPR_CUBIC_SA");
    if (cubicResidual &&
        (!seriesSetting || std::atoi(seriesSetting) != 0)) {
        series = buildCubicSeries(
            orbit, pixelParameter,
            { dx * halfWidth, dy * halfHeight },
            mxit, bailout, &_flag_halt);
        if (series.iteration < 8)
            series = {};
        if (seriesIterations) *seriesIterations = series.iteration;
        if (minimumSeriesIterations > 0 &&
            series.iteration < minimumSeriesIterations) {
            if (usedPerturbation) *usedPerturbation = false;
            return ComputeExpression(
                center_re, center_im, scale, program, fixed,
                pixelParameter, mxit, bailout, coloring, nullptr);
        }
    }
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < _h; ++i) {
        if (_flag_halt) continue;
        bool rowCompleted = true;
        double pixelIm = startIm + dy * i;
        for (int j = 0; j < _w; ++j) {
            if (_flag_halt) { rowCompleted = false; break; }
            formula::ExpressionContext context = fixed;
            formula::Complex pixel{ startRe + dx * j, pixelIm };
            if (pixelParameter == FormulaParameter::C) context.c = pixel;
            else context.z0 = pixel;
            formula::Complex parameterDelta{};
            formula::Complex delta;
            formula::Complex pixelGridDelta{
                dx * (j - halfWidth), dy * (i - halfHeight)
            };
            if (cubicResidual) {
                if (pixelParameter == FormulaParameter::C) {
                    parameterDelta = pixelGridDelta;
                    context.c = reference.c + parameterDelta;
                    delta = context.z0 - orbit[0];
                } else {
                    context.z0 = reference.z0 + pixelGridDelta;
                    delta = pixelGridDelta;
                }
            } else {
                delta = context.z0 - orbit[0];
            }
            float result = -2.0f;
            formula::Complex absolute = orbit[0] + delta;
            formula::Complex initialDelta = delta;
            formula::Complex initialAbsolute = absolute;
            formula::Complex derivative =
                pixelParameter == FormulaParameter::InitialZ
                    ? formula::Complex{ 1.0, 0.0 }
                    : formula::Complex{};
            int firstIteration = 0;
            bool initialEscaped =
                !std::isfinite(initialAbsolute.real()) ||
                !std::isfinite(initialAbsolute.imag()) ||
                std::hypot(
                    initialAbsolute.real(),
                    initialAbsolute.imag()) > bailout;
            if (!initialEscaped &&
                cubicResidual && series.iteration > 0) {
                evaluateCubicSeries(
                    series, pixelGridDelta, delta, derivative);
                firstIteration = series.iteration;
                absolute = orbit[(size_t)firstIteration] + delta;
                bool landingEscaped =
                    !std::isfinite(absolute.real()) ||
                    !std::isfinite(absolute.imag()) ||
                    std::hypot(
                        absolute.real(), absolute.imag()) > bailout;
                if (landingEscaped) {
                    delta = initialDelta;
                    derivative =
                        pixelParameter == FormulaParameter::InitialZ
                            ? formula::Complex{ 1.0, 0.0 }
                            : formula::Complex{};
                    firstIteration = 0;
                    absolute = initialAbsolute;
                }
            }
            if (!std::isfinite(absolute.real()) ||
                !std::isfinite(absolute.imag()) ||
                std::hypot(absolute.real(), absolute.imag()) > bailout) {
                result = cubicResidual
                    ? initialPowerResult(
                        absolute, pixelParameter, coloring,
                        3, bailout, dx)
                    : 0.0f;
            } else {
                std::array<formula::Complex,
                           formula::ExpressionProgram::MAX_STACK> stack;
                for (int n = firstIteration; n < mxit; ++n) {
                    if ((n & 255) == 0 && _flag_halt) {
                        rowCompleted = false;
                        break;
                    }
                    formula::Complex next;
                    if (cubicResidual) {
                        formula::Complex absoluteCurrent =
                            orbit[n] + delta;
                        formula::Complex absoluteSquared =
                            absoluteCurrent * absoluteCurrent;
                        formula::Complex nextDerivative =
                            (3.0 * absoluteSquared) * derivative;
                        if (pixelParameter == FormulaParameter::C)
                            nextDerivative += 1.0;
                        formula::Complex deltaSquared = delta * delta;
                        formula::Complex referenceSquared =
                            orbit[n] * orbit[n];
                        delta =
                            (3.0 * referenceSquared) * delta +
                            (3.0 * orbit[n]) * deltaSquared +
                            deltaSquared * delta + parameterDelta;
                        next = orbit[n + 1] + delta;
                        derivative = nextDerivative;
                    } else {
                        context.iteration = n;
                        context.z = orbit[n] + delta;
                        next = program.evaluate(
                            context, stack.data(), program.stackDepth());
                        delta = next - orbit[n + 1];
                    }
                    if (!std::isfinite(next.real()) ||
                        !std::isfinite(next.imag()) ||
                        std::hypot(next.real(), next.imag()) > bailout) {
                        double magnitude =
                            std::hypot(next.real(), next.imag());
                        if (coloring ==
                                formula::ExpressionColoring::Smooth &&
                            std::isfinite(magnitude) &&
                            magnitude > 1.0) {
                            result = (float)(n + 1 -
                                std::log(std::log(magnitude)) /
                                std::log(3.0));
                        } else if (coloring ==
                                       formula::ExpressionColoring::Distance &&
                                   std::isfinite(magnitude)) {
                            double denominator =
                                std::abs(derivative) * std::fabs(dx);
                            result = denominator > 0.0 &&
                                     std::isfinite(denominator)
                                ? (float)(
                                    magnitude * std::log(magnitude) /
                                    denominator)
                                : 0.0f;
                        } else {
                            result = (float)(n + 1);
                        }
                        break;
                    }
                }
            }
            if (!rowCompleted) break;
            _iter[i * _w + j] = result;
        }
        if (rowCompleted) progressAdvance();
    }
    bool completed = !_flag_halt;
    if (completed) progressSet(1.0);
    return completed;
}
