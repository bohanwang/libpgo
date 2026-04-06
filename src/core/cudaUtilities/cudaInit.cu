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

  int device = deviceID;  // Select the first device (you can modify this to select a specific device)
  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, device);
  SPDLOG_LOGGER_INFO(pgo::Logging::lgr(), "Using CUDA device {}: {}", device, deviceProp.name);
  SPDLOG_LOGGER_INFO(pgo::Logging::lgr(), "  Total global memory: {} bytes ({} MB)",
    deviceProp.totalGlobalMem, deviceProp.totalGlobalMem / (1024 * 1024));
  SPDLOG_LOGGER_INFO(pgo::Logging::lgr(), "  Compute capability: {}.{}", deviceProp.major, deviceProp.minor);
  SPDLOG_LOGGER_INFO(pgo::Logging::lgr(), "  Multiprocessor count: {}", deviceProp.multiProcessorCount);
  SPDLOG_LOGGER_INFO(pgo::Logging::lgr(), "  Max shared memory per block: {} bytes ({} KB)",
    deviceProp.sharedMemPerBlockOptin, deviceProp.sharedMemPerBlockOptin / 1024);

  return 0;

  err = cudaSetDevice(device);
  if (err != cudaSuccess) {
    SPDLOG_LOGGER_ERROR(pgo::Logging::lgr(), "Failed to set CUDA device {}: {}", device, cudaGetErrorString(err));
    return 1;
  }

  return 0;  // Success
}