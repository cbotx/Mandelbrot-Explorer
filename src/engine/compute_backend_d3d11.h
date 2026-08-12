#ifndef MANDEL_COMPUTE_BACKEND_D3D11_H
#define MANDEL_COMPUTE_BACKEND_D3D11_H

#include <memory>
#include <string>

class IComputeBackend;

std::unique_ptr<IComputeBackend> createD3D11ComputeBackend(bool warp, std::unique_ptr<IComputeBackend> cpuFallback, std::string* error);

#endif
