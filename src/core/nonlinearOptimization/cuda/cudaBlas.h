#pragma once

#include "cudaDeviceBuffer.h"

#include <cublas_v2.h>

namespace pgo::NonlinearOptimization::detail
{

inline void checkBlas(cublasStatus_t status, const char *operation)
{
  if (status != CUBLAS_STATUS_SUCCESS)
    throw std::runtime_error(std::string(operation) + ": cuBLAS status " + std::to_string(static_cast<int>(status)));
}

class BlasHandle
{
public:
  BlasHandle() { checkBlas(cublasCreate(&handle_), "cublasCreate"); }
  ~BlasHandle()
  {
    if (handle_)
      cublasDestroy(handle_);
  }
  BlasHandle(const BlasHandle &) = delete;
  BlasHandle &operator=(const BlasHandle &) = delete;
  operator cublasHandle_t() const { return handle_; }

  double dot(int n, const double *a, const double *b) const
  {
    double value = 0;
    checkBlas(cublasDdot(handle_, n, a, 1, b, 1, &value), "cublasDdot");
    return value;
  }
  double norm(int n, const double *a) const
  {
    double value = 0;
    checkBlas(cublasDnrm2(handle_, n, a, 1, &value), "cublasDnrm2");
    return value;
  }
  double maxAbs(int n, const double *a) const
  {
    const int index = maxAbsIndex(n, a);
    double value = 0;
    if (index >= 0)
      CudaUtilities::checkCuda(cudaMemcpy(&value, a + index, sizeof(double), cudaMemcpyDeviceToHost), "CUDA infinity norm");
    return value < 0 ? -value : value;
  }
  int maxAbsIndex(int n, const double *a) const
  {
    int index = 0;
    checkBlas(cublasIdamax(handle_, n, a, 1, &index), "cublasIdamax");
    return index - 1;
  }

private:
  cublasHandle_t handle_ = nullptr;
};

// A separate handle prevents device-result reductions from changing the
// pointer-mode assumptions of ordinary energy and Newton scalar evaluations.
class DeviceBlasHandle
{
public:
  DeviceBlasHandle()
  {
    checkBlas(cublasCreate(&handle_), "cublasCreate device reductions");
    const auto status = cublasSetPointerMode(handle_, CUBLAS_POINTER_MODE_DEVICE);
    if (status != CUBLAS_STATUS_SUCCESS) {
      cublasDestroy(handle_);
      handle_ = nullptr;
      checkBlas(status, "cublasSetPointerMode device reductions");
    }
  }
  ~DeviceBlasHandle()
  {
    if (handle_)
      cublasDestroy(handle_);
  }
  DeviceBlasHandle(const DeviceBlasHandle &) = delete;
  DeviceBlasHandle &operator=(const DeviceBlasHandle &) = delete;
  void dot(int n, const double *a, const double *b, double *deviceResult) const
  {
    checkBlas(cublasDdot(handle_, n, a, 1, b, 1, deviceResult), "cublasDdot device result");
  }
  void norm(int n, const double *a, double *deviceResult) const
  {
    checkBlas(cublasDnrm2(handle_, n, a, 1, deviceResult), "cublasDnrm2 device result");
  }

private:
  cublasHandle_t handle_ = nullptr;
};

}  // namespace pgo::NonlinearOptimization::detail
