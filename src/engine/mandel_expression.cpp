#include "mandel_perturbation.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

#include "float_math.h"
#include "formula_expression_jit.h"

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
    progressSet(0.0);
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
        if (integerPower == 0 && program.avx2Compatible()) {
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
                formula::ExpressionJitInput4 jitInput;
                formula::ExpressionJitOutput4 jitOutput;
                const bool useJit = jit && jit->valid();
                if (useJit) jitInput.setContexts(contexts);
                for (int n = 0; n < mxit && activeCount > 0; ++n) {
                    for (int lane = 0; lane < 4; ++lane)
                        contexts[lane].iteration = n;
                    if (useJit) {
                        for (int lane = 0; lane < 4; ++lane) {
                            jitInput.vectors[formula::ExpressionJitInput4::Z_RE][lane] =
                                contexts[lane].z.real();
                            jitInput.vectors[formula::ExpressionJitInput4::Z_IM][lane] =
                                contexts[lane].z.imag();
                            jitInput.vectors[formula::ExpressionJitInput4::N_RE][lane] =
                                (double)n;
                        }
                        jit->evaluate(jitInput, jitOutput);
                        for (int lane = 0; lane < 4; ++lane)
                            outputs[lane] = { jitOutput.re[lane], jitOutput.im[lane] };
                    } else {
                        if (!program.evaluate4(contexts, outputs)) {
                            rowCompleted = false;
                            break;
                        }
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
                for (int lane = 0; lane < lanes; ++lane)
                    _iter[i * _w + j + lane] = results[lane];
                if (!rowCompleted) break;
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
            _iter[i * _w + j] = result;
        }
        if (rowCompleted) progressAdvance();
    }
    bool completed = !_flag_halt;
    if (completed) progressSet(1.0);
    return completed;
}
