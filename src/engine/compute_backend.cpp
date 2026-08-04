#include "compute_backend.h"

#include <cmath>
#include <cstring>
#include <mutex>

#include "formula_expression_jit.h"
#include "mandel_perturbation.h"

namespace {

class CpuComputeBackend final : public IComputeBackend {
public:
    explicit CpuComputeBackend(bool fallback, const char* requested) {
        _info.name = "CPU";
        _info.detail = "OpenMP + AVX2";
        _info.fallback = fallback;
        if (fallback) {
            _info.detail += " (requested ";
            _info.detail += requested ? requested : "unknown";
            _info.detail += "; unavailable)";
        }
    }

    const ComputeBackendInfo& info() const override { return _info; }

    bool compute(const ComputeRequest& request) override {
        if (!valid(request)) return false;
        {
            std::lock_guard<std::mutex> lock(_mutex);
            if (_active) return false;
            _active = request.cpuEngine;
            if (_cancelRequested) {
                _active->SetHalt(true);
                _active = nullptr;
                return false;
            }
            _active->SetHalt(false);
        }
        ActiveGuard guard{ this };

        bool result = true;
        switch (request.mode) {
        case ComputeMode::Mandelbrot:
            request.cpuEngine->Compute(
                request.centerRe, request.centerIm, request.scale,
                request.maxIterations, request.coloringMethod);
            break;
        case ComputeMode::Julia:
            request.cpuEngine->ComputeJulia(
                request.centerRe, request.centerIm, request.scale,
                request.fixedCRe, request.fixedCIm,
                request.maxIterations, request.coloringMethod);
            break;
        case ComputeMode::Expression:
            result = request.cpuEngine->ComputeExpression(
                request.centerRe, request.centerIm, request.scale,
                *request.expression, *request.expressionFixed,
                request.expressionPixel, request.maxIterations,
                request.expressionBailout, request.expressionColoring,
                request.expressionJit);
            break;
        }
        std::lock_guard<std::mutex> lock(_mutex);
        return result && !_cancelRequested;
    }

    void cancel() override {
        std::lock_guard<std::mutex> lock(_mutex);
        _cancelRequested = true;
        if (_active) _active->SetHalt(true);
    }

    void resetCancellation() override {
        std::lock_guard<std::mutex> lock(_mutex);
        if (!_active) _cancelRequested = false;
    }

private:
    ComputeBackendInfo _info;
    struct ActiveGuard {
        CpuComputeBackend* backend;
        ~ActiveGuard() { backend->clearActive(); }
    };

    std::mutex _mutex;
    Mandel* _active = nullptr;
    bool _cancelRequested = false;

    void clearActive() {
        std::lock_guard<std::mutex> lock(_mutex);
        _active = nullptr;
    }

    static bool valid(const ComputeRequest& request) {
        if (!request.cpuEngine || !request.centerRe || !request.centerIm ||
            !request.scale || mpf_sgn(request.scale) <= 0 ||
            request.width < 2 || request.height < 2 ||
            request.sub < 1 || request.maxIterations < 1 ||
            !request.iterations)
            return false;
        switch (request.mode) {
        case ComputeMode::Mandelbrot:
            return true;
        case ComputeMode::Julia:
            return request.sub == 1 && request.fixedCRe && request.fixedCIm &&
                   (request.coloringMethod &
                    ~ColoringMethod::EXTERIOR_DIST_EST) == 0;
        case ComputeMode::Expression:
            return request.expression && request.expression->valid() &&
                   request.expressionFixed &&
                   request.sub == 1 &&
                   std::isfinite(request.expressionBailout) &&
                   request.expressionBailout > 0.0 &&
                   (request.expressionPixel == FormulaParameter::C ||
                    request.expressionPixel == FormulaParameter::InitialZ);
        default:
            return false;
        }
    }
};

} // namespace

std::unique_ptr<IComputeBackend> createComputeBackend(const char* requested) {
    bool cpu = !requested || !*requested || _stricmp(requested, "cpu") == 0 ||
               _stricmp(requested, "auto") == 0;
    return std::make_unique<CpuComputeBackend>(!cpu, requested);
}
