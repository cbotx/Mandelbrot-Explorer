#include "compute_backend.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <mutex>

#include "compute_backend_d3d11.h"
#include "formula_expression_jit.h"
#include "mandel_perturbation.h"

namespace {

class CpuComputeBackend final : public IComputeBackend {
public:
    explicit CpuComputeBackend(bool fallback, const char* requested,
                               const std::string& reason = {}) {
        _info.name = "CPU";
        _info.detail = "OpenMP + AVX2";
        _info.fallback = fallback;
        if (fallback) {
            _info.detail += " (requested ";
            _info.detail += requested ? requested : "unknown";
            _info.detail += "; unavailable";
            if (!reason.empty()) {
                _info.detail += ": ";
                _info.detail += reason;
            }
            _info.detail += ")";
        }
    }

    const ComputeBackendInfo& info() const override { return _info; }
    bool lastComputeUsedGpuPath() const override { return false; }

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
        {
            const char* cubicSetting =
                std::getenv("MANDEL_CUBIC_RESIDUAL");
            const char* residualPowerSetting =
                std::getenv("MANDEL_EXPR_RESIDUAL_POWER");
            const char* seriesSetting =
                std::getenv("MANDEL_EXPR_CUBIC_SA");
            const char* thresholdSetting =
                std::getenv("MANDEL_CUBIC_RESIDUAL_SCALE");
            double threshold = thresholdSetting
                ? std::atof(thresholdSetting) : 1e8;
            bool useCubicResidual =
                (!cubicSetting || std::atoi(cubicSetting) != 0) &&
                (!residualPowerSetting ||
                 std::atoi(residualPowerSetting) != 0) &&
                (!seriesSetting || std::atoi(seriesSetting) != 0) &&
                std::isfinite(threshold) && threshold > 0.0 &&
                request.expression->fastIntegerPower() == 3 &&
                request.expressionPixel == FormulaParameter::C &&
                mpf_cmp_d(request.scale, threshold) >= 0;
            if (useCubicResidual) {
                result = request.cpuEngine->ComputeExpressionResidual(
                    request.centerRe, request.centerIm, request.scale,
                    *request.expression, *request.expressionFixed,
                    request.expressionPixel, request.maxIterations,
                    request.expressionBailout, nullptr,
                    request.expressionColoring, nullptr,
                    std::max(8, request.maxIterations * 9 / 10));
            } else {
                result = request.cpuEngine->ComputeExpression(
                    request.centerRe, request.centerIm, request.scale,
                    *request.expression, *request.expressionFixed,
                    request.expressionPixel, request.maxIterations,
                    request.expressionBailout, request.expressionColoring,
                    request.expressionJit, request.expressionPlan);
            }
            break;
        }
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
    if (cpu) return std::make_unique<CpuComputeBackend>(false, requested);

    const bool gpu = _stricmp(requested, "gpu") == 0 ||
                     _stricmp(requested, "d3d11") == 0;
    const bool warp = _stricmp(requested, "warp") == 0 ||
                      _stricmp(requested, "d3d11-warp") == 0;
    if (gpu || warp) {
        std::string error;
        auto backend = createD3D11ComputeBackend(
            warp, std::make_unique<CpuComputeBackend>(false, "cpu"), &error);
        if (backend) return backend;
        return std::make_unique<CpuComputeBackend>(true, requested, error);
    }
    return std::make_unique<CpuComputeBackend>(
        true, requested, "unknown backend name");
}
