#include "cudaInit.h"
#include "cudaUtils.h"

#include <cuda_runtime.h>

int pgo::CudaUtilities::initDevice(int deviceID)
{
  int deviceCount = 0;
  cudaError_t err = cudaGetDeviceCount(&deviceCount);
  if (err != cudaSuccess || deviceCount == 0) {
    SPDLOG_LOGGER_ERROR(pgo::Logging::lgr(), "No CUDA devices found: {}", cudaGetErrorString(err));
    return 1;
  }

  if (deviceID < 0 || deviceID >= deviceCount) {
    SPDLOG_LOGGER_ERROR(pgo::Logging::lgr(), "CUDA device {} is outside [0, {}).", deviceID, deviceCount);
    return 1;
  }
  const int device = deviceID;
  err = cudaSetDevice(device);
  if (err != cudaSuccess) {
    SPDLOG_LOGGER_ERROR(pgo::Logging::lgr(), "Failed to set CUDA device {}: {}", device, cudaGetErrorString(err));
    return 1;
  }
  cudaDeviceProp deviceProp{};
  err = cudaGetDeviceProperties(&deviceProp, device);
  if (err != cudaSuccess) {
    SPDLOG_LOGGER_ERROR(pgo::Logging::lgr(), "Failed to query CUDA device {}: {}", device, cudaGetErrorString(err));
    return 1;
  }
  SPDLOG_LOGGER_INFO(pgo::Logging::lgr(), "Using CUDA device {}: {}", device, deviceProp.name);
  SPDLOG_LOGGER_INFO(pgo::Logging::lgr(), "  Total global memory: {} bytes ({} MB)",
    deviceProp.totalGlobalMem, deviceProp.totalGlobalMem / (1024 * 1024));
  SPDLOG_LOGGER_INFO(pgo::Logging::lgr(), "  Compute capability: {}.{}", deviceProp.major, deviceProp.minor);
  SPDLOG_LOGGER_INFO(pgo::Logging::lgr(), "  Multiprocessor count: {}", deviceProp.multiProcessorCount);
  SPDLOG_LOGGER_INFO(pgo::Logging::lgr(), "  Max shared memory per block: {} bytes ({} KB)",
    deviceProp.sharedMemPerBlockOptin, deviceProp.sharedMemPerBlockOptin / 1024);

  return 0;  // Success
}
