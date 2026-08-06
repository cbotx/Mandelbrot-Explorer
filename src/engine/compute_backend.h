#ifndef MANDEL_COMPUTE_BACKEND_H
#define MANDEL_COMPUTE_BACKEND_H

#include <gmp.h>

#include <memory>
#include <string>

#include "formula_expression.h"
#include "formula_spec.h"

class Mandel;

enum class ComputeMode {
    Mandelbrot,
    Julia,
    Expression
};

struct ComputeRequest {
    ComputeMode mode = ComputeMode::Mandelbrot;
    Mandel* cpuEngine = nullptr;

    mpf_ptr centerRe = nullptr;
    mpf_ptr centerIm = nullptr;
    mpf_ptr scale = nullptr;
    mpf_ptr fixedCRe = nullptr;
    mpf_ptr fixedCIm = nullptr;

    int width = 0;
    int height = 0;
    int sub = 1;
    int maxIterations = 0;
    int coloringMethod = 0;
    float* iterations = nullptr;
    float* normal = nullptr;

    const formula::ExpressionProgram* expression = nullptr;
    const formula::ExpressionContext* expressionFixed = nullptr;
    const formula::ExpressionOrbitPlan* expressionPlan = nullptr;
    const formula::ExpressionJit4* expressionJit = nullptr;
    FormulaParameter expressionPixel = FormulaParameter::C;
    formula::ExpressionColoring expressionColoring =
        formula::ExpressionColoring::Raw;
    double expressionBailout = 4.0;
};

struct ComputeBackendInfo {
    std::string name;
    std::string detail;
    bool hardwareAccelerated = false;
    bool fallback = false;
};

class IComputeBackend {
public:
    virtual ~IComputeBackend() = default;
    virtual const ComputeBackendInfo& info() const = 0;
    virtual bool lastComputeUsedGpuPath() const = 0;
    virtual bool compute(const ComputeRequest& request) = 0;
    virtual void cancel() = 0;
    virtual void resetCancellation() = 0;
};

std::unique_ptr<IComputeBackend> createComputeBackend(const char* requested);

#endif
