#ifndef MANDEL_ORBIT_OVERLAY_H
#define MANDEL_ORBIT_OVERLAY_H

#include <gmp.h>

#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <memory>
#include <mutex>
#include <thread>
#include <vector>

#include "formula_expression_orbit.h"
#include "formula_spec.h"

struct OrbitPoint {
    float re;
    float im;
};

struct OrbitResult {
    std::vector<OrbitPoint> points;
    double pixelRe = 0.0;
    double pixelIm = 0.0;
    FormulaParameter pixelParameter = FormulaParameter::C;
    double computeMs = 0.0;
    int iterations = 0;
    bool escaped = false;
    uint64_t generation = 0;
};

class OrbitWorker {
public:
    OrbitWorker();
    ~OrbitWorker();

    OrbitWorker(const OrbitWorker&) = delete;
    OrbitWorker& operator=(const OrbitWorker&) = delete;

    void request(mpf_srcptr centerRe, mpf_srcptr centerIm, mpf_srcptr scale,
                 int pixelX, int pixelY, int width, int height, int maxIterations,
                 FormulaContext formulaContext,
                 std::shared_ptr<const formula::ExpressionOrbitSnapshot>
                     expression = nullptr);
    bool takeLatest(OrbitResult& result);
    void cancel();

private:
    void run();

    std::thread _thread;
    std::mutex _mutex;
    std::condition_variable _wake;
    bool _stop = false;
    bool _hasRequest = false;

    mpf_t _centerRe, _centerIm, _scale;
    int _pixelX = 0, _pixelY = 0, _width = 1, _height = 1, _maxIterations = 1;
    FormulaContext _formula = quadraticMandelbrot();
    std::shared_ptr<const formula::ExpressionOrbitSnapshot> _expression;
    uint64_t _pendingGeneration = 0;
    std::atomic<uint64_t> _requestedGeneration{0};

    OrbitResult _latest;
    uint64_t _deliveredGeneration = 0;
};

#endif
