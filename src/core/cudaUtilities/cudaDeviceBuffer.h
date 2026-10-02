#pragma once

#include <cuda_runtime_api.h>

#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

namespace pgo::CudaUtilities
{

inline void checkCuda(cudaError_t status, const char *operation)
{
  if (status != cudaSuccess)
    throw std::runtime_error(std::string(operation) + ": " + cudaGetErrorString(status));
}

// Owning device allocation. Transfers are synchronous; kernels and zero() use
// the default CUDA stream. resize() retains capacity but does not preserve data
// when it has to grow. Device objects never escape through Eigen host views.
template<typename T>
class DeviceBuffer
{
  static_assert(std::is_trivially_copyable_v<T>);

public:
  DeviceBuffer() = default;
  explicit DeviceBuffer(std::size_t size) { resize(size); }
  ~DeviceBuffer()
  {
    if (data_)
      cudaFree(data_);
  }
  DeviceBuffer(const DeviceBuffer &) = delete;
  DeviceBuffer &operator=(const DeviceBuffer &) = delete;
  DeviceBuffer(DeviceBuffer &&other) noexcept { swap(other); }
  DeviceBuffer &operator=(DeviceBuffer &&other) noexcept
  {
    if (this != &other) {
      DeviceBuffer tmp(std::move(other));
      swap(tmp);
    }
    return *this;
  }

  void resize(std::size_t size)
  {
    if (size > capacity_) {
      if (size > std::numeric_limits<std::size_t>::max() / sizeof(T))
        throw std::length_error("CUDA device buffer allocation overflows size_t.");
      T *replacement = nullptr;
      checkCuda(cudaMalloc(reinterpret_cast<void **>(&replacement), size * sizeof(T)), "cudaMalloc");
      if (data_)
        cudaFree(data_);
      data_ = replacement;
      capacity_ = size;
    }
    size_ = size;
  }

  void upload(const T *host, std::size_t size)
  {
    resize(size);
    if (size)
      checkCuda(cudaMemcpy(data_, host, size * sizeof(T), cudaMemcpyHostToDevice), "CUDA upload");
  }
  void download(T *host, std::size_t size) const
  {
    if (size > size_)
      throw std::out_of_range("CUDA download exceeds buffer size.");
    if (size)
      checkCuda(cudaMemcpy(host, data_, size * sizeof(T), cudaMemcpyDeviceToHost), "CUDA download");
  }
  void zero()
  {
    if (size_)
      checkCuda(cudaMemsetAsync(data_, 0, size_ * sizeof(T)), "CUDA zero buffer");
  }
  T *data() { return data_; }
  const T *data() const { return data_; }
  std::size_t size() const { return size_; }
  void swap(DeviceBuffer &other) noexcept
  {
    std::swap(data_, other.data_);
    std::swap(size_, other.size_);
    std::swap(capacity_, other.capacity_);
  }

private:
  T *data_ = nullptr;
  std::size_t size_ = 0;
  std::size_t capacity_ = 0;
};

}  // namespace pgo::CudaUtilities
