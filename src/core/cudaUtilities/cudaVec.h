#pragma once

#include "cudaUtils.h"

#include <cmath>

namespace pgo
{
namespace CudaUtilities
{
template<typename T, int N>
class
#if !defined(__CUDACC__)
  alignas(32)
#endif
    CudaVec
{
public:
  __host__ __device__ CudaVec() {}
  __host__ __device__ explicit CudaVec(const T &value)
  {
#pragma unroll
    for (int i = 0; i < N; ++i) {
      dataPtr[i] = value;
    }
  }

  __host__ __device__ explicit CudaVec(const T &v0, const T &v1, const T &v2)
  {
    static_assert(N == 3, "CudaVec constructor only supports 3 components.");
    dataPtr[0] = v0;
    dataPtr[1] = v1;
    dataPtr[2] = v2;
  }

  __host__ __device__ explicit CudaVec(const T &v0, const T &v1)
  {
    static_assert(N == 2, "CudaVec constructor only supports 2 components.");
    dataPtr[0] = v0;
    dataPtr[1] = v1;
  }

  __host__ __device__ explicit CudaVec(const T v[])
  {
#pragma unroll
    for (int i = 0; i < N; ++i) {
      dataPtr[i] = v[i];
    }
  }

  __host__ __device__ T &x()
  {
    static_assert(N >= 1, "CudaVec has no x component.");
    return dataPtr[0];
  }

  __host__ __device__ const T &x() const
  {
    static_assert(N >= 1, "CudaVec has no x component.");
    return dataPtr[0];
  }

  __host__ __device__ T &y()
  {
    static_assert(N >= 2, "CudaVec has no y component.");
    return dataPtr[1];
  }

  __host__ __device__ const T &y() const
  {
    static_assert(N >= 2, "CudaVec has no y component.");
    return dataPtr[1];
  }

  __host__ __device__ T &z()
  {
    static_assert(N >= 3, "CudaVec has no z component.");
    return dataPtr[2];
  }

  __host__ __device__ const T &z() const
  {
    static_assert(N >= 3, "CudaVec has no z component.");
    return dataPtr[2];
  }

  __host__ __device__ T &operator[](int index)
  {
    return dataPtr[index];
  }

  __host__ __device__ const T &operator[](int index) const
  {
    return dataPtr[index];
  }

  __host__ __device__ CudaVec<T, N> operator+(const CudaVec<T, N> &other) const
  {
    CudaVec<T, N> result;

#pragma unroll
    for (int i = 0; i < N; ++i) {
      result[i] = dataPtr[i] + other[i];
    }
    return result;
  }

  __host__ __device__ CudaVec<T, N> operator-(const CudaVec<T, N> &other) const
  {
    CudaVec<T, N> result;

#pragma unroll
    for (int i = 0; i < N; ++i) {
      result[i] = dataPtr[i] - other[i];
    }
    return result;
  }

  __host__ __device__ CudaVec<T, N> operator*(const T &scalar) const
  {
    CudaVec<T, N> result;

#pragma unroll
    for (int i = 0; i < N; ++i) {
      result[i] = dataPtr[i] * scalar;
    }
    return result;
  }

  __host__ __device__ CudaVec<T, N> operator/(const T &scalar) const
  {
    CudaVec<T, N> result;

#pragma unroll
    for (int i = 0; i < N; ++i) {
      result[i] = dataPtr[i] / scalar;
    }
    return result;
  }

  __host__ __device__ CudaVec<T, N> &operator+=(const CudaVec<T, N> &other)
  {
#pragma unroll
    for (int i = 0; i < N; ++i) {
      dataPtr[i] += other[i];
    }
    return *this;
  }

  __host__ __device__ CudaVec<T, N> &operator-=(const CudaVec<T, N> &other)
  {
#pragma unroll
    for (int i = 0; i < N; ++i) {
      dataPtr[i] -= other[i];
    }
    return *this;
  }

  __host__ __device__ CudaVec<T, N> &operator*=(const T &scalar)
  {
#pragma unroll
    for (int i = 0; i < N; ++i) {
      dataPtr[i] *= scalar;
    }
    return *this;
  }

  __host__ __device__ CudaVec<T, N> &operator/=(const T &scalar)
  {
#pragma unroll
    for (int i = 0; i < N; ++i) {
      dataPtr[i] /= scalar;
    }
    return *this;
  }

  template<typename T2>
  __host__ __device__ CudaVec<T2, N> cast() const
  {
    CudaVec<T2, N> result;

#pragma unroll
    for (int i = 0; i < N; ++i) {
      result[i] = static_cast<T2>(dataPtr[i]);
    }
    return result;
  }

  __host__ __device__ T norm() const
  {
    T s = 0;

#pragma unroll
    for (int i = 0; i < N; ++i) {
      s += dataPtr[i] * dataPtr[i];
    }
    return std::sqrt(s);
  }

  __host__ __device__ T length() const
  {
    return norm();
  }

  __host__ __device__ T squaredNorm() const
  {
    T s = 0;

#pragma unroll
    for (int i = 0; i < N; ++i) {
      s += dataPtr[i] * dataPtr[i];
    }
    return s;
  }

  __host__ __device__ T dot(const CudaVec<T, N> &other) const
  {
    T s = 0;

#pragma unroll
    for (int i = 0; i < N; ++i) {
      s += dataPtr[i] * other.dataPtr[i];
    }
    return s;
  }

  __host__ __device__ T *data()
  {
    return dataPtr;
  }

  __host__ __device__ const T *data() const
  {
    return dataPtr;
  }

private:
  T dataPtr[N];
};

using CudaV3i = CudaVec<int, 3>;
using CudaV2i = CudaVec<int, 2>;
using CudaV3f = CudaVec<float, 3>;
}  // namespace CudaUtilities
}  // namespace pgo
