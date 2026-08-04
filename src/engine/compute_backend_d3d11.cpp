#include "compute_backend_d3d11.h"

#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include <d3d11.h>
#include <d3dcompiler.h>
#include <dxgi1_2.h>
#include <wrl/client.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "compute_backend.h"
#include "mandel_perturbation.h"

using Microsoft::WRL::ComPtr;

namespace {

constexpr UINT kThreadsPerGroup = 64;
constexpr UINT kPixelsPerChunk = 524288;
constexpr int kMaxGpuIterations = 100000;
constexpr UINT kGpuLongPrefix = 256;
constexpr int kCpuRefineBlock = 128;
constexpr float kAnalyticInterior = -3.0f;
constexpr float kCpuRefineIteration = 64.0f;

const char kShaderSource[] = R"(
cbuffer Params : register(b0) {
    float2 c0Re;
    float2 c0Im;
    float2 dx;
    float2 dy;
    uint width;
    uint height;
    uint maxIterations;
    uint basePixel;
    uint pixelCount;
    float2 escapeSquared;
    uint padding;
};

RWStructuredBuffer<float> Output : register(u0);

float2 dsAdd(float2 a, float2 b) {
    precise float s = a.x + b.x;
    precise float v = s - a.x;
    precise float e = (a.x - (s - v)) + (b.x - v);
    e = e + a.y + b.y;
    precise float hi = s + e;
    precise float lo = e - (hi - s);
    return float2(hi, lo);
}

float2 dsSub(float2 a, float2 b) {
    return dsAdd(a, float2(-b.x, -b.y));
}

float2 dsMul(float2 a, float2 b) {
    const float splitter = 8193.0f;
    precise float p = a.x * b.x;

    precise float ca = splitter * a.x;
    precise float ah = ca - (ca - a.x);
    precise float al = a.x - ah;
    precise float cb = splitter * b.x;
    precise float bh = cb - (cb - b.x);
    precise float bl = b.x - bh;

    precise float e = ((ah * bh - p) + ah * bl + al * bh) + al * bl;
    e = e + a.x * b.y + a.y * b.x + a.y * b.y;
    precise float hi = p + e;
    precise float lo = e - (hi - p);
    return float2(hi, lo);
}

float2 dsScale(float2 a, float b) {
    return dsMul(a, float2(b, 0.0f));
}

bool dsGreater(float2 a, float2 b) {
    return a.x > b.x || (a.x == b.x && a.y > b.y);
}

bool inMainCardioidOrBulb(float2 x, float2 y) {
    float2 y2 = dsMul(y, y);
    float2 xm = dsAdd(x, float2(-0.25f, 0.0f));
    float2 q = dsAdd(dsMul(xm, xm), y2);
    float2 cardioidLeft = dsMul(q, dsAdd(q, xm));
    float2 cardioidRight = dsScale(y2, 0.25f);
    bool cardioid = dsGreater(cardioidRight, cardioidLeft);

    float2 xp = dsAdd(x, float2(1.0f, 0.0f));
    float2 bulbRadius = dsAdd(dsMul(xp, xp), y2);
    bool bulb = dsGreater(float2(0.0625f, 0.0f), bulbRadius);
    return cardioid || bulb;
}

float solvePixelDs(uint pixel, uint prefix) {
    uint y = pixel / width;
    uint x = pixel - y * width;
    float2 cr = dsAdd(c0Re, dsScale(dx, (float)x));
    float2 ci = dsAdd(c0Im, dsScale(dy, (float)y));
    if (inMainCardioidOrBulb(cr, ci)) return -3.0f;

    float2 zr = cr;
    float2 zi = ci;
    float dr = 2.0f;
    float di = 0.0f;
    uint i = 1;
    uint iterationLimit = min(maxIterations, prefix);
    while (i < iterationLimit) {
        precise float nextDr = 2.0f * (dr * zr.x - di * zi.x);
        precise float nextDi = 2.0f * (dr * zi.x + di * zr.x);
        float2 zr2 = dsMul(zr, zr);
        float2 zi2 = dsMul(zi, zi);
        float2 zri = dsMul(zr, zi);
        zr = dsAdd(dsSub(zr2, zi2), cr);
        zi = dsAdd(dsScale(zri, 2.0f), ci);
        dr = nextDr;
        di = nextDi;

        float2 radiusSquared = dsAdd(dsMul(zr, zr), dsMul(zi, zi));
        if (dsGreater(radiusSquared, escapeSquared)) {
            float radius = radiusSquared.x + radiusSquared.y;
            const float invLog2 = 1.4426950408889634f;
            return (float)(i + 1) -
                   log(log(radius) * (0.5f * invLog2)) * invLog2;
        }
        if (dr * dr + di * di < 1.0e-9f) return -4.0f;
        ++i;
    }
    return -5.0f;
}

bool inMainCardioidOrBulbFast(float x, float y) {
    precise float y2 = y * y;
    precise float xm = x - 0.25f;
    precise float q = xm * xm + y2;
    bool cardioid = q * (q + xm) < 0.25f * y2;
    precise float xp = x + 1.0f;
    bool bulb = xp * xp + y2 < 0.0625f;
    return cardioid || bulb;
}

float solvePixelFast(uint pixel) {
    uint y = pixel / width;
    uint x = pixel - y * width;
    precise float fx = (float)x;
    precise float fy = (float)y;
    precise float cr = c0Re.x + dx.x * fx;
    cr = cr + (c0Re.y + dx.y * fx);
    precise float ci = c0Im.x + dy.x * fy;
    ci = ci + (c0Im.y + dy.y * fy);
    if (inMainCardioidOrBulbFast(cr, ci)) return -3.0f;

    precise float zr = cr;
    precise float zi = ci;
    precise float dr = 2.0f;
    precise float di = 0.0f;
    float threshold = escapeSquared.x + escapeSquared.y;
    uint i = 1;
    uint iterationLimit = min(maxIterations, 256u);
    while (i < iterationLimit) {
        precise float nextDr = 2.0f * (dr * zr - di * zi);
        precise float nextDi = 2.0f * (dr * zi + di * zr);
        precise float nextZr = zr * zr - zi * zi + cr;
        precise float nextZi = 2.0f * zr * zi + ci;
        zr = nextZr;
        zi = nextZi;
        dr = nextDr;
        di = nextDi;

        precise float radiusSquared = zr * zr + zi * zi;
        if (radiusSquared > threshold) {
            const float invLog2 = 1.4426950408889634f;
            return (float)(i + 1) -
                   log(log(radiusSquared) * (0.5f * invLog2)) * invLog2;
        }
        if (dr * dr + di * di < 1.0e-9f) return -4.0f;
        ++i;
    }
    return -5.0f;
}

[numthreads(64, 1, 1)]
void main(uint3 dispatchThreadId : SV_DispatchThreadID) {
    uint localPixel = dispatchThreadId.x;
    if (localPixel >= pixelCount) return;
    uint pixel = basePixel + localPixel;
    Output[pixel] = solvePixelDs(pixel, 64u);
}

[numthreads(64, 1, 1)]
void mainFast(uint3 dispatchThreadId : SV_DispatchThreadID) {
    uint localPixel = dispatchThreadId.x;
    if (localPixel >= pixelCount) return;
    uint pixel = basePixel + localPixel;
    Output[pixel] = solvePixelFast(pixel);
}

[numthreads(64, 1, 1)]
void mainLong(uint3 dispatchThreadId : SV_DispatchThreadID) {
    uint localPixel = dispatchThreadId.x;
    if (localPixel >= pixelCount) return;
    uint pixel = basePixel + localPixel;
    float previous = Output[pixel];
    if (previous == -3.0f ||
        (previous >= 0.0f && previous < 64.0f))
        return;
    float candidate = solvePixelDs(pixel, 256u);
    if (previous == -4.0f && candidate == -4.0f) {
        Output[pixel] = -7.0f;
    } else if (previous >= 0.0f && candidate >= 0.0f &&
               (uint)previous == (uint)candidate &&
               abs(previous - candidate) <= 0.125f) {
        Output[pixel] = candidate;
    } else if (previous == -5.0f && candidate == -5.0f) {
        Output[pixel] = -5.0f;
    } else {
        Output[pixel] = -6.0f;
    }
}
)";

std::string hresultText(HRESULT hr) {
    char message[64];
    sprintf_s(message, "HRESULT 0x%08X", static_cast<unsigned>(hr));
    return message;
}

std::string narrow(const wchar_t* text) {
    if (!text || !*text) return {};
    int size = WideCharToMultiByte(CP_UTF8, 0, text, -1, nullptr, 0, nullptr, nullptr);
    if (size <= 1) return {};
    std::string result(static_cast<size_t>(size), '\0');
    WideCharToMultiByte(CP_UTF8, 0, text, -1, result.data(), size, nullptr, nullptr);
    result.pop_back();
    return result;
}

struct FloatPair {
    float hi;
    float lo;
};

FloatPair splitDouble(double value) {
    FloatPair result{};
    result.hi = static_cast<float>(value);
    result.lo = static_cast<float>(value - static_cast<double>(result.hi));
    return result;
}

bool inMainCardioidOrBulb(double x, double y) {
    const double y2 = y * y;
    const double xm = x - 0.25;
    const double q = xm * xm + y2;
    return q * (q + xm) < 0.25 * y2 ||
           (x + 1.0) * (x + 1.0) + y2 < 0.0625;
}

struct alignas(16) ShaderParams {
    FloatPair c0Re;
    FloatPair c0Im;
    FloatPair dx;
    FloatPair dy;
    UINT width;
    UINT height;
    UINT maxIterations;
    UINT basePixel;
    UINT pixelCount;
    FloatPair escapeSquared;
    UINT padding;
};
static_assert(sizeof(ShaderParams) == 64, "HLSL constant layout mismatch");

class D3D11ComputeBackend final : public IComputeBackend {
public:
    D3D11ComputeBackend(bool warp, std::unique_ptr<IComputeBackend> cpuFallback)
        : _warp(warp), _cpu(std::move(cpuFallback)) {}

    bool initialize(std::string* error) {
        UINT flags = 0;
        const D3D_FEATURE_LEVEL requested[] = {
            D3D_FEATURE_LEVEL_11_1, D3D_FEATURE_LEVEL_11_0
        };
        D3D_DRIVER_TYPE driver = _warp
            ? D3D_DRIVER_TYPE_WARP : D3D_DRIVER_TYPE_HARDWARE;
        HRESULT hr = D3D11CreateDevice(
            nullptr, driver, nullptr, flags, requested,
            static_cast<UINT>(sizeof(requested) / sizeof(requested[0])),
            D3D11_SDK_VERSION,
            &_device, &_featureLevel, &_context);
        if (hr == E_INVALIDARG) {
            hr = D3D11CreateDevice(
                nullptr, driver, nullptr, flags, requested + 1, 1,
                D3D11_SDK_VERSION, &_device, &_featureLevel, &_context);
        }
        if (FAILED(hr)) {
            if (error) *error = "D3D11 device creation failed (" + hresultText(hr) + ")";
            return false;
        }
        if (_featureLevel < D3D_FEATURE_LEVEL_11_0) {
            if (error) *error = "D3D feature level 11.0 is required";
            return false;
        }
        DXGI_ADAPTER_DESC1 adapter{};
        bool haveAdapter = adapterDescription(&adapter);
        bool softwareAdapter =
            haveAdapter &&
            ((adapter.Flags & DXGI_ADAPTER_FLAG_SOFTWARE) != 0 ||
             (adapter.VendorId == 0x1414 && adapter.DeviceId == 0x008c));
        if (!_warp && softwareAdapter) {
            if (error) {
                *error = "no physical D3D11 adapter (found " +
                         narrow(adapter.Description) + ")";
            }
            return false;
        }

        const UINT compileFlags =
            D3DCOMPILE_ENABLE_STRICTNESS |
            D3DCOMPILE_IEEE_STRICTNESS |
            D3DCOMPILE_OPTIMIZATION_LEVEL3;
        auto compileShader = [&](const char* entry,
                                 ComPtr<ID3D11ComputeShader>& shader) {
            ComPtr<ID3DBlob> shaderBlob;
            ComPtr<ID3DBlob> errors;
            HRESULT result = D3DCompile(
                kShaderSource, sizeof(kShaderSource) - 1,
                "mandelbrot_shallow.hlsl", nullptr, nullptr, entry, "cs_5_0",
                compileFlags, 0, &shaderBlob, &errors);
            if (FAILED(result)) {
                if (error) {
                    if (errors && errors->GetBufferPointer())
                        *error = static_cast<const char*>(errors->GetBufferPointer());
                    else
                        *error = std::string(entry) +
                                 " shader compilation failed (" +
                                 hresultText(result) + ")";
                }
                return false;
            }
            result = _device->CreateComputeShader(
                shaderBlob->GetBufferPointer(), shaderBlob->GetBufferSize(),
                nullptr, &shader);
            if (FAILED(result) && error)
                *error = std::string(entry) + " shader creation failed (" +
                         hresultText(result) + ")";
            return SUCCEEDED(result);
        };
        if (!compileShader("main", _shader) ||
            !compileShader("mainFast", _shaderFast) ||
            !compileShader("mainLong", _shaderLong))
            return false;

        D3D11_BUFFER_DESC constants{};
        constants.ByteWidth = sizeof(ShaderParams);
        constants.Usage = D3D11_USAGE_DEFAULT;
        constants.BindFlags = D3D11_BIND_CONSTANT_BUFFER;
        hr = _device->CreateBuffer(&constants, nullptr, &_constantBuffer);
        if (FAILED(hr)) {
            if (error) *error = "constant buffer creation failed (" + hresultText(hr) + ")";
            return false;
        }

        D3D11_QUERY_DESC query{};
        query.Query = D3D11_QUERY_EVENT;
        hr = _device->CreateQuery(&query, &_completionQuery);
        if (FAILED(hr)) {
            if (error) *error = "completion query creation failed (" + hresultText(hr) + ")";
            return false;
        }

        std::string adapterName =
            haveAdapter ? narrow(adapter.Description) : std::string();
        _info.name = _warp ? "D3D11 WARP" : "D3D11 GPU";
        _info.detail = adapterName.empty()
            ? "FP32/2xFP32 bounded-prefix compute"
            : adapterName + "; FP32/2xFP32 bounded-prefix compute";
        _info.detail += "; max 100k iterations";
        _info.hardwareAccelerated = !_warp;
        _info.fallback = false;
        return true;
    }

    const ComputeBackendInfo& info() const override { return _info; }
    bool lastComputeUsedGpuPath() const override {
        return _lastComputeUsedGpuPath.load(std::memory_order_acquire);
    }

    bool compute(const ComputeRequest& request) override {
        bool expected = false;
        if (!_computing.compare_exchange_strong(expected, true))
            return false;
        ComputeGuard guard{_computing};
        _lastComputeUsedGpuPath.store(false, std::memory_order_release);
        if (_cancelRequested.load(std::memory_order_acquire))
            return false;

        if (!_gpuHealthy || !eligible(request))
            return computeCpu(request);

        if (computeGpu(request)) {
            _lastComputeUsedGpuPath.store(true, std::memory_order_release);
            return !_cancelRequested.load(std::memory_order_acquire);
        }
        if (_cancelRequested.load(std::memory_order_acquire))
            return false;

        if (FAILED(_device->GetDeviceRemovedReason())) _gpuHealthy = false;
        return computeCpu(request);
    }

    void cancel() override {
        _cancelRequested.store(true, std::memory_order_release);
        _cpu->cancel();
    }

    void resetCancellation() override {
        if (!_computing.load(std::memory_order_acquire)) {
            _cancelRequested.store(false, std::memory_order_release);
            _cpu->resetCancellation();
        }
    }

private:
    struct ComputeGuard {
        std::atomic_bool& flag;
        ~ComputeGuard() { flag.store(false, std::memory_order_release); }
    };

    bool _warp = false;
    bool _gpuHealthy = true;
    std::unique_ptr<IComputeBackend> _cpu;
    ComputeBackendInfo _info;
    std::atomic_bool _computing{false};
    std::atomic_bool _cancelRequested{false};
    std::atomic_bool _lastComputeUsedGpuPath{false};

    D3D_FEATURE_LEVEL _featureLevel = D3D_FEATURE_LEVEL_9_1;
    ComPtr<ID3D11Device> _device;
    ComPtr<ID3D11DeviceContext> _context;
    ComPtr<ID3D11ComputeShader> _shader;
    ComPtr<ID3D11ComputeShader> _shaderFast;
    ComPtr<ID3D11ComputeShader> _shaderLong;
    ComPtr<ID3D11Buffer> _constantBuffer;
    ComPtr<ID3D11Query> _completionQuery;
    ComPtr<ID3D11Buffer> _outputBuffer;
    ComPtr<ID3D11UnorderedAccessView> _outputUav;
    ComPtr<ID3D11Buffer> _stagingBuffer;
    UINT _outputCapacity = 0;
    std::vector<float> _readback;
    std::vector<int> _refinePixels;
    std::vector<double> _refineRe;
    std::vector<double> _refineIm;
    std::vector<float> _refineOutput;

    bool computeCpu(const ComputeRequest& request) {
        if (_cancelRequested.load(std::memory_order_acquire))
            _cpu->cancel();
        else
            _cpu->resetCancellation();
        return _cpu->compute(request);
    }

    static bool eligible(const ComputeRequest& request) {
        return request.mode == ComputeMode::Mandelbrot &&
               request.cpuEngine && request.centerRe && request.centerIm &&
               request.scale && mpf_sgn(request.scale) > 0 &&
               mpf_cmp_d(request.scale, 1.0e6) <= 0 &&
               request.width >= 2 && request.height >= 2 &&
               request.sub >= 1 && (request.sub & 1) != 0 &&
               request.maxIterations >= 2 &&
               request.maxIterations <= kMaxGpuIterations &&
               request.iterations &&
               request.coloringMethod == 0;
    }

    static bool fastFloatSafe(const double* values, int width, int height) {
        const double lastRe = values[0] + values[2] * (width - 1);
        const double lastIm = values[1] + values[3] * (height - 1);
        const double maxCoordinate = std::max({
            std::fabs(values[0]), std::fabs(values[1]),
            std::fabs(lastRe), std::fabs(lastIm), 1.0
        });
        float coordinate = static_cast<float>(maxCoordinate);
        float next = std::nextafter(
            coordinate, std::numeric_limits<float>::infinity());
        double ulp = static_cast<double>(next) - coordinate;
        double pixelStep =
            std::min(std::fabs(values[2]), std::fabs(values[3]));
        return std::isfinite(next) && pixelStep >= ulp * 16.0;
    }

    bool ensureOutput(UINT count) {
        if (_outputBuffer && _outputUav && _stagingBuffer &&
            count <= _outputCapacity)
            return true;
        _outputCapacity = 0;
        _outputBuffer.Reset();
        _outputUav.Reset();
        _stagingBuffer.Reset();

        D3D11_BUFFER_DESC output{};
        output.ByteWidth = count * sizeof(float);
        output.Usage = D3D11_USAGE_DEFAULT;
        output.BindFlags = D3D11_BIND_UNORDERED_ACCESS;
        output.MiscFlags = D3D11_RESOURCE_MISC_BUFFER_STRUCTURED;
        output.StructureByteStride = sizeof(float);
        HRESULT hr = _device->CreateBuffer(&output, nullptr, &_outputBuffer);
        if (FAILED(hr)) return false;

        D3D11_UNORDERED_ACCESS_VIEW_DESC outputView{};
        outputView.Format = DXGI_FORMAT_UNKNOWN;
        outputView.ViewDimension = D3D11_UAV_DIMENSION_BUFFER;
        outputView.Buffer.NumElements = count;
        hr = _device->CreateUnorderedAccessView(
            _outputBuffer.Get(), &outputView, &_outputUav);
        if (FAILED(hr)) return false;

        D3D11_BUFFER_DESC staging = output;
        staging.Usage = D3D11_USAGE_STAGING;
        staging.BindFlags = 0;
        staging.CPUAccessFlags = D3D11_CPU_ACCESS_READ;
        staging.MiscFlags = 0;
        staging.StructureByteStride = 0;
        hr = _device->CreateBuffer(&staging, nullptr, &_stagingBuffer);
        if (FAILED(hr)) return false;
        _outputCapacity = count;
        return true;
    }

    bool computeGpu(const ComputeRequest& request) {
        using ProfileClock = std::chrono::steady_clock;
        const auto totalStart = ProfileClock::now();
        const bool profile = getenv("MANDEL_GPU_PROFILE") != nullptr;
        const uint64_t count64 =
            static_cast<uint64_t>(request.width) * request.height;
        if (count64 == 0 ||
            count64 > std::numeric_limits<UINT>::max() / sizeof(float))
            return false;
        const UINT count = static_cast<UINT>(count64);
        if (!ensureOutput(count)) return false;

        const mp_bitcnt_t precision = std::max({
            mpf_get_prec(request.centerRe), mpf_get_prec(request.centerIm),
            mpf_get_prec(request.scale)
        });
        mpf_t dw, dh, c0Re, c0Im, dx, dy;
        mpf_init2(dw, precision); mpf_init2(dh, precision);
        mpf_init2(c0Re, precision); mpf_init2(c0Im, precision);
        mpf_init2(dx, precision); mpf_init2(dy, precision);
        mpf_set_ui(dw, 2); mpf_div(dw, dw, request.scale);
        mpf_set(dh, dw); mpf_div_ui(dh, dh, request.width);
        mpf_mul_ui(dh, dh, request.height);
        mpf_sub(c0Re, request.centerRe, dw);
        mpf_sub(c0Im, request.centerIm, dh);
        mpf_mul_ui(dx, dw, 2); mpf_div_ui(dx, dx, request.width - 1);
        mpf_mul_ui(dy, dh, 2); mpf_div_ui(dy, dy, request.height - 1);

        const double values[] = {
            mpf_get_d(c0Re), mpf_get_d(c0Im), mpf_get_d(dx), mpf_get_d(dy)
        };
        mpf_clears(dw, dh, c0Re, c0Im, dx, dy, (mpf_ptr)0);
        for (double value : values) {
            if (!std::isfinite(value) || std::fabs(value) > 1.0e8)
                return false;
        }
        const double lastRe =
            values[0] + values[2] * (request.width - 1);
        const double lastIm =
            values[1] + values[3] * (request.height - 1);
        if (!std::isfinite(lastRe) || !std::isfinite(lastIm) ||
            std::fabs(lastRe) > 1.0e8 || std::fabs(lastIm) > 1.0e8)
            return false;

        ShaderParams params{};
        params.c0Re = splitDouble(values[0]);
        params.c0Im = splitDouble(values[1]);
        params.dx = splitDouble(values[2]);
        params.dy = splitDouble(values[3]);
        const FloatPair coordinatePairs[] = {
            params.c0Re, params.c0Im, params.dx, params.dy
        };
        for (const FloatPair pair : coordinatePairs)
            if (!std::isfinite(pair.hi) || !std::isfinite(pair.lo))
                return false;
        params.width = static_cast<UINT>(request.width);
        params.height = static_cast<UINT>(request.height);
        params.maxIterations = static_cast<UINT>(request.maxIterations);
        double escapeRadius = request.cpuEngine->escapeRadius();
        if (!std::isfinite(escapeRadius) || escapeRadius < 2.0 ||
            escapeRadius > 1.0e8)
            return false;
        double escapeSquared = request.cpuEngine->escapeRadiusSquared();
        if (!std::isfinite(escapeSquared)) return false;
        params.escapeSquared = splitDouble(escapeSquared);
        if (!std::isfinite(params.escapeSquared.hi) ||
            !std::isfinite(params.escapeSquared.lo))
            return false;

        const bool fastFloat =
            fastFloatSafe(values, request.width, request.height);
        const float cpuRefineIteration =
            static_cast<float>(kGpuLongPrefix);
        const UINT gpuIterationPrefix = kGpuLongPrefix;

        ID3D11Buffer* constantBuffers[] = {_constantBuffer.Get()};
        ID3D11UnorderedAccessView* uavs[] = {_outputUav.Get()};
        _context->CSSetShader(
            fastFloat ? _shaderFast.Get() : _shader.Get(), nullptr, 0);
        _context->CSSetConstantBuffers(0, 1, constantBuffers);
        _context->CSSetUnorderedAccessViews(0, 1, uavs, nullptr);

        const auto dispatchStart = ProfileClock::now();
        for (UINT base = 0; base < count; base += kPixelsPerChunk) {
            if (_cancelRequested.load(std::memory_order_acquire)) {
                unbind();
                return false;
            }
            params.basePixel = base;
            params.pixelCount = std::min(kPixelsPerChunk, count - base);
            _context->UpdateSubresource(
                _constantBuffer.Get(), 0, nullptr, &params, 0, 0);
            UINT groups =
                (params.pixelCount + kThreadsPerGroup - 1) / kThreadsPerGroup;
            _context->Dispatch(groups, 1, 1);
            _context->CSSetShader(_shaderLong.Get(), nullptr, 0);
            _context->Dispatch(groups, 1, 1);
            _context->End(_completionQuery.Get());
            _context->Flush();
            HRESULT completed = S_FALSE;
            unsigned spins = 0;
            do {
                completed = _context->GetData(
                    _completionQuery.Get(), nullptr, 0, 0);
                if (completed == S_FALSE) {
                    if ((++spins & 255u) == 0)
                        SwitchToThread();
                    else
                        YieldProcessor();
                }
            } while (completed == S_FALSE);
            if (FAILED(completed)) {
                unbind();
                return false;
            }
            if (_cancelRequested.load(std::memory_order_acquire)) {
                unbind();
                return false;
            }
            _context->CSSetShader(
                fastFloat ? _shaderFast.Get() : _shader.Get(), nullptr, 0);
        }
        const auto dispatchEnd = ProfileClock::now();
        unbind();

        const auto readbackStart = ProfileClock::now();
        _context->CopyResource(_stagingBuffer.Get(), _outputBuffer.Get());
        D3D11_MAPPED_SUBRESOURCE mapped{};
        HRESULT hr = _context->Map(
            _stagingBuffer.Get(), 0, D3D11_MAP_READ, 0, &mapped);
        if (FAILED(hr)) return false;
        _readback.resize(count);
        std::memcpy(_readback.data(), mapped.pData, count * sizeof(float));
        _context->Unmap(_stagingBuffer.Get(), 0);
        if (_cancelRequested.load(std::memory_order_acquire)) return false;
        const auto readbackEnd = ProfileClock::now();

        // The GPU resolves robust short escapes only. Long or bounded trajectories
        // are sparse and numerically sensitive, so finish those with the engine's
        // exact FP64 scalar semantics instead of letting divergent GPU lanes run
        // to the full iteration limit.
        const auto refineStart = ProfileClock::now();
        _refinePixels.resize(count);
        std::atomic<int> refineCount{0};
        long long analytic = 0, trustedInterior = 0;
        long long derivativeInterior = 0, prefixTail = 0;
        long long unstable = 0, longEscape = 0;
#pragma omp parallel for schedule(static) reduction(+:analytic,trustedInterior,derivativeInterior,prefixTail,unstable,longEscape)
        for (int pixel = 0; pixel < static_cast<int>(count); ++pixel) {
            if (_cancelRequested.load(std::memory_order_relaxed)) continue;
            float value = _readback[pixel];
            int y = pixel / request.width;
            int x = pixel - y * request.width;
            double cr = values[0] + values[2] * x;
            double ci = values[1] + values[3] * y;
            if (value == kAnalyticInterior) {
                if (inMainCardioidOrBulb(cr, ci)) {
                    _readback[pixel] = -2.0f;
                    ++analytic;
                    continue;
                }
            } else if (value == -7.0f) {
                _readback[pixel] = -2.0f;
                ++trustedInterior;
                continue;
            } else if (value >= 0.0f && value < cpuRefineIteration) {
                continue;
            }
            if (value == -4.0f) ++derivativeInterior;
            else if (value == -5.0f) ++prefixTail;
            else if (value == -6.0f) ++unstable;
            else ++longEscape;
            int slot = refineCount.fetch_add(1, std::memory_order_relaxed);
            _refinePixels[slot] = pixel;
        }
        const int refined = refineCount.load(std::memory_order_relaxed);
        _refineRe.resize(refined);
        _refineIm.resize(refined);
        _refineOutput.resize(refined);
#pragma omp parallel for schedule(static)
        for (int slot = 0; slot < refined; ++slot) {
            int pixel = _refinePixels[slot];
            int y = pixel / request.width;
            int x = pixel - y * request.width;
            _refineRe[slot] = values[0] + values[2] * x;
            _refineIm[slot] = values[1] + values[3] * y;
        }
        const int refineBlocks =
            (refined + kCpuRefineBlock - 1) / kCpuRefineBlock;
#pragma omp parallel for schedule(dynamic, 1)
        for (int block = 0; block < refineBlocks; ++block) {
            if (_cancelRequested.load(std::memory_order_relaxed)) continue;
            int base = block * kCpuRefineBlock;
            int blockCount = std::min(kCpuRefineBlock, refined - base);
            request.cpuEngine->ComputeShallowPoints(
                _refineRe.data() + base, _refineIm.data() + base,
                blockCount, _refineOutput.data() + base,
                request.maxIterations);
        }
#pragma omp parallel for schedule(static)
        for (int slot = 0; slot < refined; ++slot) {
            _readback[_refinePixels[slot]] = _refineOutput[slot];
        }
        if (_cancelRequested.load(std::memory_order_acquire)) return false;
        const auto refineEnd = ProfileClock::now();

        if (request.sub == 1) {
            std::memcpy(
                request.iterations, _readback.data(), count * sizeof(float));
        } else {
            const size_t stride =
                static_cast<size_t>(request.width) * request.sub;
            const int center = request.sub / 2;
            for (int y = 0; y < request.height; ++y) {
                size_t destination =
                    (static_cast<size_t>(y) * request.sub + center) * stride +
                    center;
                const float* source =
                    _readback.data() + static_cast<size_t>(y) * request.width;
                for (int x = 0; x < request.width; ++x)
                    request.iterations[destination +
                        static_cast<size_t>(x) * request.sub] = source[x];
            }
        }
        if (profile) {
            auto milliseconds = [](auto begin, auto end) {
                return std::chrono::duration<double, std::milli>(
                    end - begin).count();
            };
            const auto totalEnd = ProfileClock::now();
            fprintf(stderr,
                    "  gpu phases: dispatch=%.3f ms readback=%.3f ms "
                    "refine=%.3f ms total=%.3f ms "
                    "gpu-only=%lld analytic=%lld trusted=%lld refined=%lld "
                    "(deriv=%lld tail=%lld unstable=%lld escape=%lld) "
                    "prefix=%u mode=%s\n",
                    milliseconds(dispatchStart, dispatchEnd),
                    milliseconds(readbackStart, readbackEnd),
                    milliseconds(refineStart, refineEnd),
                    milliseconds(totalStart, totalEnd),
                    static_cast<long long>(count) - analytic -
                        trustedInterior - refined,
                    analytic, trustedInterior,
                    static_cast<long long>(refined),
                    derivativeInterior, prefixTail, unstable, longEscape,
                    gpuIterationPrefix,
                    fastFloat ? "fp32+2xfp32" : "2xfp32");
            fflush(stderr);
        }
        return true;
    }

    void unbind() {
        ID3D11UnorderedAccessView* nullUav = nullptr;
        ID3D11Buffer* nullBuffer = nullptr;
        _context->CSSetUnorderedAccessViews(0, 1, &nullUav, nullptr);
        _context->CSSetConstantBuffers(0, 1, &nullBuffer);
        _context->CSSetShader(nullptr, nullptr, 0);
    }

    bool adapterDescription(DXGI_ADAPTER_DESC1* description) const {
        if (!description) return false;
        ComPtr<IDXGIDevice> dxgiDevice;
        ComPtr<IDXGIAdapter> adapter;
        ComPtr<IDXGIAdapter1> adapter1;
        if (FAILED(_device.As(&dxgiDevice)) ||
            FAILED(dxgiDevice->GetAdapter(&adapter)) ||
            FAILED(adapter.As(&adapter1)) ||
            FAILED(adapter1->GetDesc1(description)))
            return false;
        return true;
    }
};

} // namespace

std::unique_ptr<IComputeBackend> createD3D11ComputeBackend(
    bool warp, std::unique_ptr<IComputeBackend> cpuFallback,
    std::string* error) {
    auto backend =
        std::make_unique<D3D11ComputeBackend>(warp, std::move(cpuFallback));
    if (!backend->initialize(error)) return nullptr;
    return backend;
}
