#ifndef MANDEL_ORBIT_THUMBNAIL_H
#define MANDEL_ORBIT_THUMBNAIL_H

#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <functional>
#include <memory>
#include <mutex>
#include <thread>
#include <vector>

#include "formula_expression_orbit.h"

struct OrbitThumbnailPixel {
    int iterations = 0;
    bool escaped = false;
    bool cancelled = false;
};

bool classifyOrbitThumbnailPixel(
    const formula::ExpressionOrbitSnapshot* expression,
    formula::Complex pixel, int maxIterations,
    OrbitThumbnailPixel& result,
    const std::function<bool()>& shouldCancel = {});

struct OrbitThumbnailResult {
    std::vector<uint8_t> pixels;
    int width = 0;
    int height = 0;
    int maxIterations = 0;
    bool expression = false;
    FormulaParameter pixelParameter = FormulaParameter::C;
    double computeMs = 0.0;
    uint64_t generation = 0;
};

class OrbitThumbnailWorker {
public:
    OrbitThumbnailWorker();
    ~OrbitThumbnailWorker();

    OrbitThumbnailWorker(const OrbitThumbnailWorker&) = delete;
    OrbitThumbnailWorker& operator=(const OrbitThumbnailWorker&) = delete;

    uint64_t request(
        int width, int height,
        double x0, double x1, double y0, double y1,
        int maxIterations,
        std::shared_ptr<const formula::ExpressionOrbitSnapshot> expression);
    bool takeLatest(OrbitThumbnailResult& result);
    void cancel();

private:
    void run();

    std::thread _thread;
    std::mutex _mutex;
    std::condition_variable _wake;
    bool _stop = false;
    bool _hasRequest = false;
    bool _latestReady = false;

    int _width = 1;
    int _height = 1;
    int _maxIterations = 1;
    double _x0 = -2.5;
    double _x1 = 1.0;
    double _y0 = -1.2;
    double _y1 = 1.2;
    std::shared_ptr<const formula::ExpressionOrbitSnapshot> _expression;
    uint64_t _pendingGeneration = 0;
    std::atomic<uint64_t> _requestedGeneration{0};
    OrbitThumbnailResult _latest;
};

#endif
