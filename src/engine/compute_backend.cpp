#include "compute_backend.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <mutex>

#include "compute_backend_d3d11.h"
#include "formula_deep_renderer.h"
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
    bool lastComputeUsedCustomDeepPath() const override {
        return _lastCustomDeep.load(std::memory_order_acquire);
    }
    bool lastComputeUsedGenericDeepPath() const override {
        return _lastGenericDeep.load(std::memory_order_acquire);
    }
    GenericDeepInfo lastGenericDeepInfo() const override {
        std::lock_guard<std::mutex> lock(_genericMutex);
        return _genericInfo;
    }

    bool compute(const ComputeRequest& request) override {
        if (!valid(request)) return false;
        {
            std::lock_guard<std::mutex> lock(_mutex);
            if (_active) return false;
            _active = request.cpuEngine;
            if (_cancelRequested.load(std::memory_order_acquire)) {
                _active->SetHalt(true);
                _active = nullptr;
                return false;
            }
            _active->SetHalt(false);
        }
        ActiveGuard guard{ this };
        _lastCustomDeep.store(false, std::memory_order_release);
        _lastGenericDeep.store(false, std::memory_order_release);
        {
            std::lock_guard<std::mutex> lock(_genericMutex);
            _genericInfo = {};
        }

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
            formula::ExpressionProductionPlan productionPlan;
            if (request.expressionSource) {
                productionPlan = formula::makeExpressionProductionPlan(
                    *request.expressionSource, *request.expression,
                    *request.expressionFixed, request.expressionPixel,
                    request.expressionBailout, request.expressionColoring,
                    request.scale, request.centerRe, request.centerIm,
                    request.width, request.height);
            }
            const bool aboveDirectLimit =
                mpf_cmp_d(
                    request.scale,
                    formula::CUSTOM_DIRECT_ZOOM_LIMIT) > 0;
            if (productionPlan.usesQuadraticPerturbation()) {
                const formula::CustomDeepZoomPlan& deepZoom =
                    productionPlan.quadratic;
                const int method =
                    request.coloringMethod &
                    ~ColoringMethod::SUPER_SAMPLING;
                const bool methodMatches =
                    (deepZoom.outputAdapter ==
                         formula::CustomDeepZoomOutputAdapter::
                             SmoothExpression &&
                     method == 0) ||
                    (deepZoom.outputAdapter ==
                         formula::CustomDeepZoomOutputAdapter::
                             DistanceExpression &&
                     method ==
                         ColoringMethod::EXTERIOR_DIST_EST) ||
                    (deepZoom.outputAdapter ==
                         formula::CustomDeepZoomOutputAdapter::
                             FeatherExpression &&
                     method ==
                         ColoringMethod::STRIPE_AVERAGE) ||
                    (deepZoom.outputAdapter ==
                         formula::CustomDeepZoomOutputAdapter::
                             OrbitTrapExpression &&
                     method ==
                         ColoringMethod::ORBIT_TRAP);
                if (!methodMatches) {
                    std::fill(
                        request.iterations,
                        request.iterations +
                            static_cast<size_t>(
                                request.width) * request.height,
                        formula::ExpressionDeepEmptyPixel);
                    result = false;
                    break;
                }
                Mandel::ScopedCustomCompute customCompute(
                    *request.cpuEngine, deepZoom.escapeRadius,
                    deepZoom.outputAdapter);
                if (!customCompute.active()) {
                    std::fill(
                        request.iterations,
                        request.iterations +
                            static_cast<size_t>(
                                request.width) * request.height,
                        formula::ExpressionDeepEmptyPixel);
                    result = false;
                    break;
                }
                request.cpuEngine->Compute(
                    request.centerRe, request.centerIm, request.scale,
                    request.maxIterations, request.coloringMethod);
                _lastCustomDeep.store(true, std::memory_order_release);
                break;
            }
            if (aboveDirectLimit) {
                if (!request.expressionSource ||
                    !productionPlan.usesGenericCertifiedDeep()) {
                    std::fill(
                        request.iterations,
                        request.iterations +
                            static_cast<size_t>(
                                request.width) * request.height,
                        formula::ExpressionDeepEmptyPixel);
                    result = false;
                    break;
                }
                result = computeGenericDeep(request);
                break;
            }
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
        return result &&
               !_cancelRequested.load(std::memory_order_acquire);
    }

    void cancel() override {
        std::lock_guard<std::mutex> lock(_mutex);
        _cancelRequested.store(true, std::memory_order_release);
        if (_active) _active->SetHalt(true);
    }

    void resetCancellation() override {
        std::lock_guard<std::mutex> lock(_mutex);
        if (!_active)
            _cancelRequested.store(false, std::memory_order_release);
    }

private:
    ComputeBackendInfo _info;
    std::atomic_bool _lastCustomDeep{false};
    std::atomic_bool _lastGenericDeep{false};
    struct ActiveGuard {
        CpuComputeBackend* backend;
        ~ActiveGuard() { backend->clearActive(); }
    };

    std::mutex _mutex;
    Mandel* _active = nullptr;
    std::atomic_bool _cancelRequested{false};
    mutable std::mutex _genericMutex;
    GenericDeepInfo _genericInfo;

    void clearActive() {
        std::lock_guard<std::mutex> lock(_mutex);
        _active = nullptr;
    }

    static const char* phaseName(
            formula::ExpressionDeepRenderPhase phase) {
        switch (phase) {
        case formula::ExpressionDeepRenderPhase::Reference:
            return "reference/Taylor planning";
        case formula::ExpressionDeepRenderPhase::Fast:
            return "certified Taylor";
        case formula::ExpressionDeepRenderPhase::Fallback:
            return "MPFR fallback";
        case formula::ExpressionDeepRenderPhase::Complete:
            return "complete";
        }
        return "unknown";
    }

    static float phaseProgress(
            formula::ExpressionDeepRenderPhase phase,
            uint64_t completed, uint64_t total) {
        const float fraction = total
            ? static_cast<float>(
                static_cast<double>(completed) /
                static_cast<double>(total))
            : 1.0f;
        switch (phase) {
        case formula::ExpressionDeepRenderPhase::Reference:
            return 0.03f;
        case formula::ExpressionDeepRenderPhase::Fast:
            return 0.12f + 0.68f * fraction;
        case formula::ExpressionDeepRenderPhase::Fallback:
            return 0.80f + 0.19f * fraction;
        case formula::ExpressionDeepRenderPhase::Complete:
            return 1.0f;
        }
        return 0.0f;
    }

    bool computeGenericDeep(const ComputeRequest& request) {
        const uint64_t pixelCount =
            static_cast<uint64_t>(request.width) *
            static_cast<uint64_t>(request.height);
        std::fill(
            request.iterations,
            request.iterations + static_cast<size_t>(pixelCount),
            formula::ExpressionDeepEmptyPixel);
        _lastGenericDeep.store(true, std::memory_order_release);
        {
            std::lock_guard<std::mutex> lock(_genericMutex);
            _genericInfo = {};
            _genericInfo.used = true;
            _genericInfo.status = "running";
            _genericInfo.phase = "reference/Taylor planning";
            _genericInfo.progress = 0.0f;
            _genericInfo.pixelCount = pixelCount;
        }
        if (request.progress)
            request.progress->store(0.0f, std::memory_order_relaxed);

        formula::ExpressionDeepRenderRequest deepRequest;
        deepRequest.canonicalProgram = request.expressionSource;
        deepRequest.runtimeProgram = request.expression;
        deepRequest.center.realMpf = request.centerRe;
        deepRequest.center.imaginaryMpf = request.centerIm;
        deepRequest.scale.mpf = request.scale;
        deepRequest.fixed = *request.expressionFixed;
        deepRequest.pixelParameter = request.expressionPixel;
        deepRequest.width = request.width;
        deepRequest.height = request.height;
        deepRequest.maxIterations = request.maxIterations;
        deepRequest.bailout = request.expressionBailout;
        deepRequest.output = request.iterations;
        deepRequest.outputCount =
            static_cast<size_t>(pixelCount);
        deepRequest.precision.viewBits =
            static_cast<mpfr_prec_t>(std::max({
                mpf_get_prec(request.centerRe),
                mpf_get_prec(request.centerIm),
                mpf_get_prec(request.scale)
            }));
        deepRequest.shouldCancel = [this, engine = request.cpuEngine] {
            return _cancelRequested.load(std::memory_order_acquire) ||
                   engine->HaltRequested();
        };
        deepRequest.progress =
            [this, progress = request.progress](
                formula::ExpressionDeepRenderPhase phase,
                uint64_t completed, uint64_t total) {
                const float value =
                    phaseProgress(phase, completed, total);
                if (progress)
                    progress->store(
                        value, std::memory_order_relaxed);
                std::lock_guard<std::mutex> lock(_genericMutex);
                _genericInfo.phase = phaseName(phase);
                _genericInfo.progress = value;
            };

        const auto started = std::chrono::steady_clock::now();
        formula::ExpressionDeepRenderResult deepResult;
        const bool success =
            formula::renderExpressionDeepFrame(
                deepRequest, deepResult);
        const double totalSeconds =
            std::chrono::duration<double>(
                std::chrono::steady_clock::now() - started).count();
        {
            std::lock_guard<std::mutex> lock(_genericMutex);
            _genericInfo.settled = true;
            _genericInfo.success = success;
            _genericInfo.cancelled = deepResult.cancelled;
            _genericInfo.taylorAccepted =
                deepResult.taylorAccepted;
            _genericInfo.status =
                formula::expressionDeepRenderStatusName(
                    deepResult.status);
            _genericInfo.error = deepResult.error;
            _genericInfo.phase =
                success ? "complete" :
                (deepResult.cancelled ? "cancelled" : "failed");
            _genericInfo.progress =
                success ? 1.0f : _genericInfo.progress;
            _genericInfo.fastPixelCount =
                deepResult.fastPixelCount;
            _genericInfo.fallbackPixelCount =
                deepResult.fallbackPixelCount;
            _genericInfo.taylorPixelCoverage =
                deepResult.taylorAcceptedPixelCoverage;
            _genericInfo.totalSeconds = totalSeconds;
            _genericInfo.referenceSeconds =
                deepResult.referenceSeconds;
            _genericInfo.taylorSeconds =
                deepResult.taylorBuildSeconds +
                deepResult.fastSeconds;
            _genericInfo.fastSeconds = deepResult.fastSeconds;
            _genericInfo.fallbackSeconds =
                deepResult.fallbackSeconds;
        }
        if (success && request.progress)
            request.progress->store(1.0f, std::memory_order_relaxed);
        return success;
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
