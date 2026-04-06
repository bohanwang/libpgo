#pragma once

#include "cudaUtils.h"

#include "pgoLogging.h"

#include <atomic>

#include <stdlib.h>

namespace pgo
{
namespace CudaUtilities
{

template<typename T>
struct LinearMappedMemory
{
  T *deviceData, *hostData;
  int n;

  LinearMappedMemory(int n_);
  ~LinearMappedMemory();

protected:
  int allocateMemory();
  int deallocateMemory();
};

template<typename T>
int LinearMappedMemory<T>::allocateMemory()
{
  CUDA_CHECK(cudaHostAlloc(&hostData, n * sizeof(T), cudaHostAllocDefault | cudaHostAllocMapped), return -1);
  CUDA_CHECK(cudaHostGetDevicePointer(&deviceData, hostData, 0), return -1);
  CUDA_CHECK(cudaDeviceSynchronize(), return -1);

  unsigned int flags = 0;
  CUDA_CHECK(cudaHostGetFlags(&flags, hostData), return -1);
  if (flags & cudaHostAllocMapped) {
    SPDLOG_LOGGER_DEBUG(pgo::Logging::lgr(), "Memory is mapped");
  }
  else {
    SPDLOG_LOGGER_DEBUG(pgo::Logging::lgr(), "Memory is not mapped");
  }

  return 0;
}

template<typename T>
int LinearMappedMemory<T>::deallocateMemory()
{
  CUDA_CHECK(cudaFreeHost(hostData), return -1);
  hostData = nullptr;
  deviceData = nullptr;
  n = 0;

  return 0;
}

template<typename T>
LinearMappedMemory<T>::LinearMappedMemory(int n_):
  n(n_)
{
  int ret = allocateMemory();
  if (ret != 0) {
    exit(EXIT_FAILURE);
  }
}

template<typename T>
LinearMappedMemory<T>::~LinearMappedMemory()
{
  int ret = deallocateMemory();
  if (ret != 0) {
    exit(EXIT_FAILURE);
  }
}
}  // namespace CudaUtilities
}  // namespace pgo