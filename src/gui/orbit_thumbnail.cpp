#include "orbit_thumbnail.h"

#include <algorithm>
#include <chrono>
#include <cmath>

namespace {

using Clock = std::chrono::steady_clock;

} // namespace

bool classifyOrbitThumbnailPixel(const formula::ExpressionOrbitSnapshot* expression, formula::Complex pixel, int maxIterations, OrbitThumbnailPixel& result, const std::function<bool()>& shouldCancel) {
    result = OrbitThumbnailPixel{};
    if (maxIterations < 0) return false;
    if (expression) {
        formula::ExpressionOrbitClassification classification;
        if (!formula::classifyExpressionOrbit(*expression, pixel, maxIterations, classification, shouldCancel)) { return false; }
        result.iterations = classification.iterations;
        result.escaped = classification.escaped;
        result.cancelled = classification.cancelled;
        return true;
    }

    double zr = 0.0;
    double zi = 0.0;
    for (int n = 0; n < maxIterations; ++n) {
        if (shouldCancel && shouldCancel()) {
            result.cancelled = true;
            break;
        }
        double nextReal = zr * zr - zi * zi + pixel.real();
        zi = 2.0 * zr * zi + pixel.imag();
        zr = nextReal;
        result.iterations = n + 1;
        if (!std::isfinite(zr) || !std::isfinite(zi) || zr * zr + zi * zi > 16.0) {
            result.escaped = true;
            break;
        }
    }
    return true;
}

OrbitThumbnailWorker::OrbitThumbnailWorker() : _thread(&OrbitThumbnailWorker::run, this) {
}

OrbitThumbnailWorker::~OrbitThumbnailWorker() {
    {
        std::lock_guard<std::mutex> lock(_mutex);
        _requestedGeneration.fetch_add(1);
        _stop = true;
    }
    _wake.notify_one();
    if (_thread.joinable()) _thread.join();
}

uint64_t OrbitThumbnailWorker::request(int width, int height, double x0, double x1, double y0, double y1, int maxIterations, std::shared_ptr<const formula::ExpressionOrbitSnapshot> expression) {
    uint64_t generation;
    OrbitThumbnailResult staleResult;
    std::shared_ptr<const formula::ExpressionOrbitSnapshot> staleExpression;
    {
        std::lock_guard<std::mutex> lock(_mutex);
        generation = _requestedGeneration.fetch_add(1) + 1;
        staleResult = std::move(_latest);
        staleExpression = std::move(_expression);
        _width = std::max(2, width);
        _height = std::max(2, height);
        _x0 = x0;
        _x1 = x1;
        _y0 = y0;
        _y1 = y1;
        _maxIterations = std::clamp(maxIterations, 1, 160);
        _expression = std::move(expression);
        _pendingGeneration = generation;
        _hasRequest = true;
        _latestReady = false;
    }
    _wake.notify_one();
    return generation;
}

bool OrbitThumbnailWorker::takeLatest(OrbitThumbnailResult& result) {
    std::lock_guard<std::mutex> lock(_mutex);
    if (!_latestReady) return false;
    result = std::move(_latest);
    _latestReady = false;
    return true;
}

void OrbitThumbnailWorker::cancel() {
    OrbitThumbnailResult staleResult;
    std::shared_ptr<const formula::ExpressionOrbitSnapshot> staleExpression;
    {
        std::lock_guard<std::mutex> lock(_mutex);
        _requestedGeneration.fetch_add(1);
        _hasRequest = false;
        _latestReady = false;
        staleResult = std::move(_latest);
        staleExpression = std::move(_expression);
    }
}

void OrbitThumbnailWorker::run() {
    for (;;) {
        int width, height, maxIterations;
        double x0, x1, y0, y1;
        uint64_t generation;
        std::shared_ptr<const formula::ExpressionOrbitSnapshot> expression;
        {
            std::unique_lock<std::mutex> lock(_mutex);
            _wake.wait(lock, [this] { return _stop || _hasRequest; });
            if (_stop) return;
            width = _width;
            height = _height;
            maxIterations = _maxIterations;
            x0 = _x0;
            x1 = _x1;
            y0 = _y0;
            y1 = _y1;
            generation = _pendingGeneration;
            expression = _expression;
            _hasRequest = false;
        }

        OrbitThumbnailResult result;
        result.width = width;
        result.height = height;
        result.maxIterations = maxIterations;
        result.expression = expression != nullptr;
        if (expression) result.pixelParameter = expression->pixelParameter;
        result.generation = generation;
        result.pixels.assign((size_t)width * height * 3, 0);
        auto begin = Clock::now();
        auto cancelled = [this, generation] { return _requestedGeneration.load() != generation; };

        bool complete = true;
        for (int y = 0; y < height && complete; ++y) {
            double imaginary = y1 - (y1 - y0) * y / std::max(1, height - 1);
            for (int x = 0; x < width; ++x) {
                if (cancelled()) {
                    complete = false;
                    break;
                }
                double real = x0 + (x1 - x0) * x / std::max(1, width - 1);
                OrbitThumbnailPixel classification;
                if (!classifyOrbitThumbnailPixel(expression.get(), {real, imaginary}, maxIterations, classification, cancelled) || classification.cancelled) {
                    complete = false;
                    break;
                }
                uint8_t* pixel = &result.pixels[((size_t)y * width + x) * 3];
                if (!classification.escaped) {
                    pixel[0] = 16;
                    pixel[1] = 13;
                    pixel[2] = 10;
                } else {
                    double t = std::sqrt((double)classification.iterations / maxIterations);
                    pixel[0] = (uint8_t)(30 + 95 * t);
                    pixel[1] = (uint8_t)(24 + 135 * t);
                    pixel[2] = (uint8_t)(20 + 215 * t);
                }
            }
        }
        if (!complete || cancelled()) continue;
        result.computeMs = std::chrono::duration<double, std::milli>(Clock::now() - begin).count();

        {
            std::lock_guard<std::mutex> lock(_mutex);
            if (_requestedGeneration.load() != generation) continue;
            _latest = std::move(result);
            _latestReady = true;
        }
    }
}
