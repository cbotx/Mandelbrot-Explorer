#ifndef MANDEL_COMPUTE_BACKEND_H
#define MANDEL_COMPUTE_BACKEND_H

#include <gmp.h>

#include <array>
#include <atomic>
#include <cstdint>
#include <memory>
#include <string>

#include "formula_expression.h"
#include "formula_spec.h"

class Mandel;

enum class ComputeMode { Mandelbrot,
                         Julia,
                         Expression };

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
    std::atomic<float>* progress = nullptr;

    const formula::ExpressionProgram* expressionSource = nullptr;
    const formula::ExpressionProgram* expression = nullptr;
    const formula::ExpressionContext* expressionFixed = nullptr;
    const formula::ExpressionOrbitPlan* expressionPlan = nullptr;
    const formula::ExpressionJit4* expressionJit = nullptr;
    FormulaParameter expressionPixel = FormulaParameter::C;
    formula::ExpressionColoring expressionColoring = formula::ExpressionColoring::Raw;
    double expressionBailout = 4.0;
};

struct ComputeBackendInfo {
    std::string name;
    std::string detail;
    bool hardwareAccelerated = false;
    bool fallback = false;
};

struct GenericDeepInfo {
    bool used = false;
    bool settled = false;
    bool success = false;
    bool cancelled = false;
    bool taylorAccepted = false;
    bool preflightAttempted = false;
    bool preflightRejectedFast = false;
    bool specializedPiecewiseMpfr = false;
    bool piecewiseBigFixed = false;
    std::string status;
    std::string error;
    std::string phase;
    std::string predictedPath;
    float progress = 0.0f;
    uint64_t pixelCount = 0;
    uint64_t preflightSampleCount = 0;
    uint64_t preflightFallbackCount = 0;
    uint64_t preflightAvoidedFastPixelCount = 0;
    uint64_t preflightIterationCount = 0;
    uint64_t preflightOperationCount = 0;
    uint64_t preflightFoldOperationCount = 0;
    std::array<uint64_t, 16> preflightFirstUncertainHistogram{};
    uint64_t fastPixelCount = 0;
    uint64_t fallbackPixelCount = 0;
    uint64_t taylorPixelCoverage = 0;
    uint64_t totalIterationCount = 0;
    uint64_t specializedPiecewiseMpfrPixelCount = 0;
    uint64_t specializedPiecewiseMpfrIterationCount = 0;
    uint64_t specializedPiecewiseMpfrPeriodicPixelCount = 0;
    uint64_t genericMpfrPeriodicPixelCount = 0;
    uint64_t piecewiseBigFixedPixelCount = 0;
    uint64_t piecewiseBigFixedIterationCount = 0;
    uint64_t selectedPrecision = 0;
    uint64_t fallbackPrecision = 0;
    double totalSeconds = 0.0;
    double referenceSeconds = 0.0;
    double preflightSeconds = 0.0;
    double taylorSeconds = 0.0;
    double fastSeconds = 0.0;
    double fallbackSeconds = 0.0;
};

class IComputeBackend {
  public:
    virtual ~IComputeBackend() = default;
    virtual const ComputeBackendInfo& info() const = 0;
    virtual bool lastComputeUsedGpuPath() const = 0;
    virtual bool lastComputeUsedCustomDeepPath() const { return false; }
    virtual bool lastComputeUsedGenericDeepPath() const { return false; }
    virtual GenericDeepInfo lastGenericDeepInfo() const { return {}; }
    virtual bool compute(const ComputeRequest& request) = 0;
    virtual void cancel() = 0;
    virtual void resetCancellation() = 0;
};

std::unique_ptr<IComputeBackend> createComputeBackend(const char* requested);

#endif
