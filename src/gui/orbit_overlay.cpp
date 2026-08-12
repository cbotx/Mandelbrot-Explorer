#include "orbit_overlay.h"

#include <algorithm>
#include <chrono>
#include <cmath>

using Clock = std::chrono::steady_clock;

OrbitWorker::OrbitWorker() {
    mpf_init(_centerRe);
    mpf_init(_centerIm);
    mpf_init(_scale);
    _thread = std::thread(&OrbitWorker::run, this);
}

OrbitWorker::~OrbitWorker() {
    {
        std::lock_guard<std::mutex> lock(_mutex);
        _requestedGeneration.fetch_add(1);
        _stop = true;
    }
    _wake.notify_one();
    if (_thread.joinable()) _thread.join();
    mpf_clear(_centerRe);
    mpf_clear(_centerIm);
    mpf_clear(_scale);
}

void OrbitWorker::request(mpf_srcptr centerRe, mpf_srcptr centerIm, mpf_srcptr scale, int pixelX, int pixelY, int width, int height, int maxIterations, FormulaContext formulaContext, std::shared_ptr<const formula::ExpressionOrbitSnapshot> expression, formula::CustomDeepZoomPlan customDeepZoom) {
    const mp_bitcnt_t precision = std::max({mpf_get_prec(centerRe), mpf_get_prec(centerIm), mpf_get_prec(scale)});
    uint64_t generation;
    OrbitResult staleResult;
    std::shared_ptr<const formula::ExpressionOrbitSnapshot> staleExpression;
    {
        std::lock_guard<std::mutex> lock(_mutex);
        generation = _requestedGeneration.fetch_add(1) + 1;
        staleResult = std::move(_latest);
        _latest = OrbitResult{};
        staleExpression = std::move(_expression);
        mpf_set_prec(_centerRe, precision);
        mpf_set_prec(_centerIm, precision);
        mpf_set_prec(_scale, precision);
        mpf_set(_centerRe, centerRe);
        mpf_set(_centerIm, centerIm);
        mpf_set(_scale, scale);
        _pixelX = pixelX;
        _pixelY = pixelY;
        _width = std::max(2, width);
        _height = std::max(2, height);
        _maxIterations = std::max(1, std::min(maxIterations, 2048));
        _formula = formulaContext;
        _expression = std::move(expression);
        _customDeepZoom = customDeepZoom;
        _pendingGeneration = generation;
        _hasRequest = true;
    }
    _wake.notify_one();
}

bool OrbitWorker::takeLatest(OrbitResult& result) {
    std::lock_guard<std::mutex> lock(_mutex);
    if (_latest.generation == 0 || _latest.generation == _deliveredGeneration) return false;
    result = _latest;
    _deliveredGeneration = _latest.generation;
    return true;
}

void OrbitWorker::cancel() {
    OrbitResult staleResult;
    std::shared_ptr<const formula::ExpressionOrbitSnapshot> staleExpression;
    {
        std::lock_guard<std::mutex> lock(_mutex);
        _requestedGeneration.fetch_add(1);
        _hasRequest = false;
        staleResult = std::move(_latest);
        _latest = OrbitResult{};
        staleExpression = std::move(_expression);
        _deliveredGeneration = 0;
    }
}

void OrbitWorker::run() {
    for (;;) {
        mp_bitcnt_t precision;
        int pixelX, pixelY, width, height, maxIterations;
        uint64_t generation;
        FormulaContext formulaContext;
        std::shared_ptr<const formula::ExpressionOrbitSnapshot> expression;
        formula::CustomDeepZoomPlan customDeepZoom;
        mpf_t centerRe, centerIm, scale;

        {
            std::unique_lock<std::mutex> lock(_mutex);
            _wake.wait(lock, [&] { return _stop || _hasRequest; });
            if (_stop) return;
            precision = std::max({mpf_get_prec(_centerRe), mpf_get_prec(_centerIm), mpf_get_prec(_scale)});
            mpf_init2(centerRe, precision);
            mpf_init2(centerIm, precision);
            mpf_init2(scale, precision);
            mpf_set(centerRe, _centerRe);
            mpf_set(centerIm, _centerIm);
            mpf_set(scale, _scale);
            pixelX = _pixelX;
            pixelY = _pixelY;
            width = _width;
            height = _height;
            maxIterations = _maxIterations;
            generation = _pendingGeneration;
            formulaContext = _formula;
            expression = _expression;
            customDeepZoom = _customDeepZoom;
            _hasRequest = false;
        }

        OrbitResult result;
        result.generation = generation;
        auto begin = Clock::now();

        mpf_t dw, dh, dx, dy, cRe, cIm, t, zr, zi, nr, ni, zr2, zi2, escapeSquared;
        mpf_inits(dw, dh, dx, dy, cRe, cIm, t, zr, zi, nr, ni, zr2, zi2, escapeSquared, (mpf_ptr)0);
        for (mpf_ptr value : {dw, dh, dx, dy, cRe, cIm, t, zr, zi, nr, ni, zr2, zi2, escapeSquared}) mpf_set_prec(value, precision);

        mpf_set_ui(dw, 2);
        mpf_div(dw, dw, scale);
        mpf_set(dh, dw);
        mpf_mul_ui(dh, dh, (unsigned long)height);
        mpf_div_ui(dh, dh, (unsigned long)width);

        mpf_mul_ui(dx, dw, 2);
        mpf_div_ui(dx, dx, (unsigned long)(width - 1));
        mpf_mul_ui(dy, dh, 2);
        mpf_div_ui(dy, dy, (unsigned long)(height - 1));

        pixelX = std::clamp(pixelX, 0, width - 1);
        pixelY = std::clamp(pixelY, 0, height - 1);
        mpf_sub(cRe, centerRe, dw);
        mpf_mul_ui(t, dx, (unsigned long)pixelX);
        mpf_add(cRe, cRe, t);
        mpf_sub(cIm, centerIm, dh);
        mpf_mul_ui(t, dy, (unsigned long)(height - 1 - pixelY));
        mpf_add(cIm, cIm, t);
        result.pixelRe = mpf_get_d(cRe);
        result.pixelIm = mpf_get_d(cIm);
        result.pixelParameter = expression ? expression->pixelParameter : formulaContext.slice.pixel;

        const bool regularQuadratic = formulaContext.formula.id == FormulaId::PowerPlusC && formulaContext.formula.power == 2 && formulaContext.slice.pixel == FormulaParameter::C;
        const bool customDeepQuadratic = formulaContext.formula.id == FormulaId::Expression && customDeepZoom.usesQuadraticPerturbation();
        if (regularQuadratic || customDeepQuadratic) {
            const double escapeRadius = customDeepQuadratic ? customDeepZoom.escapeRadius : 4.0;
            mpf_set_d(escapeSquared, escapeRadius);
            mpf_mul(escapeSquared, escapeSquared, escapeSquared);
            result.usedGmpQuadratic = true;
            mpf_set_ui(zr, 0);
            mpf_set_ui(zi, 0);
            result.points.reserve((size_t)maxIterations + 1);
            result.points.push_back({0.0f, 0.0f});
            for (int i = 0; i < maxIterations; ++i) {
                if ((i & 63) == 0 && _requestedGeneration.load() != generation) break;
                mpf_mul(nr, zr, zi);
                mpf_mul(zr2, zr, zr);
                mpf_mul(zi2, zi, zi);
                mpf_sub(zr, zr2, zi2);
                mpf_add(zr, zr, cRe);
                mpf_mul_ui(zi, nr, 2);
                mpf_add(zi, zi, cIm);
                result.points.push_back({(float)mpf_get_d(zr), (float)mpf_get_d(zi)});
                result.iterations = i + 1;

                mpf_mul(zr2, zr, zr);
                mpf_mul(zi2, zi, zi);
                mpf_add(t, zr2, zi2);
                if (mpf_cmp(t, escapeSquared) > 0) {
                    result.escaped = true;
                    break;
                }
            }
        } else if (formulaContext.formula.id == FormulaId::Expression && expression && expression->valid()) {
            formula::ExpressionOrbitEvaluation evaluation;
            formula::evaluateExpressionOrbit(*expression, {result.pixelRe, result.pixelIm}, maxIterations, evaluation, [this, generation] { return _requestedGeneration.load() != generation; });
            result.points.reserve(evaluation.points.size());
            for (const formula::Complex& point : evaluation.points) { result.points.push_back({(float)point.real(), (float)point.imag()}); }
            result.iterations = evaluation.iterations;
            result.escaped = evaluation.escaped;
        }
        result.computeMs = std::chrono::duration<double, std::milli>(Clock::now() - begin).count();
        mpf_clears(dw, dh, dx, dy, cRe, cIm, t, zr, zi, nr, ni, zr2, zi2, escapeSquared, (mpf_ptr)0);
        mpf_clears(centerRe, centerIm, scale, (mpf_ptr)0);

        {
            std::lock_guard<std::mutex> lock(_mutex);
            if (_requestedGeneration.load() != generation) continue;
            _latest = std::move(result);
        }
    }
}
