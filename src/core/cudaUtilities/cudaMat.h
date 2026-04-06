#pragma once

#include "cudaUtils.h"
#include "cudaVec.h"

#include <cassert>
namespace pgo
{
namespace CudaUtilities
{
template<typename T, int Rows, int Cols>
class
#if !defined(__CUDACC__)
  alignas(32)
#endif
    CudaMat
{
public:
  __host__ __device__ CudaMat() noexcept
  {
  }

  __host__ __device__ T &operator()(int row, int col)
  {
    assert(row >= 0 && row < Rows);
    assert(col >= 0 && col < Cols);

    return dataPtr[col * Rows + row];
  }

  __host__ __device__ const T &operator()(int row, int col) const
  {
    assert(row >= 0 && row < Rows);
    assert(col >= 0 && col < Cols);

    return dataPtr[col * Rows + row];
  }

  __host__ __device__ CudaMat<T, Rows, Cols> operator+(const CudaMat<T, Rows, Cols> &other) const
  {
    CudaMat<T, Rows, Cols> result;

#pragma unroll
    for (int i = 0; i < Rows * Cols; ++i) {
      result.dataPtr[i] = dataPtr[i] + other.dataPtr[i];
    }
    return result;
  }

  __host__ __device__ CudaMat<T, Rows, Cols> operator-(const CudaMat<T, Rows, Cols> &other) const
  {
    CudaMat<T, Rows, Cols> result;

#pragma unroll
    for (int i = 0; i < Rows * Cols; ++i) {
      result.dataPtr[i] = dataPtr[i] - other.dataPtr[i];
    }
    return result;
  }

  __host__ __device__ CudaMat<T, Rows, Cols> operator*(T scalar) const
  {
    CudaMat<T, Rows, Cols> result;
#pragma unroll
    for (int i = 0; i < Rows * Cols; ++i) {
      result.dataPtr[i] = dataPtr[i] * scalar;
    }
    return result;
  }

  __host__ __device__ CudaVec<T, Rows> operator*(const CudaVec<T, Cols> &vec) const
  {
    CudaVec<T, Rows> result;

#pragma unroll
    for (int i = 0; i < Rows; ++i) {
      T sum = T(0);

#pragma unroll
      for (int j = 0; j < Cols; ++j) {
        sum += (*this)(i, j) * vec[j];
      }
      result[i] = sum;
    }
    return result;
  }

  __host__ __device__ CudaMat<T, Rows, Cols> operator/(T scalar) const
  {
    CudaMat<T, Rows, Cols> result;
#pragma unroll
    for (int i = 0; i < Rows * Cols; ++i) {
      result.dataPtr[i] = dataPtr[i] / scalar;
    }
    return result;
  }

  __host__ __device__ CudaMat<T, Rows, Cols> &operator+=(const CudaMat<T, Rows, Cols> &other)
  {
#pragma unroll
    for (int i = 0; i < Rows * Cols; ++i) {
      dataPtr[i] += other.dataPtr[i];
    }
    return *this;
  }

  __host__ __device__ CudaMat<T, Rows, Cols> &operator-=(const CudaMat<T, Rows, Cols> &other)
  {
#pragma unroll
    for (int i = 0; i < Rows * Cols; ++i) {
      dataPtr[i] -= other.dataPtr[i];
    }
    return *this;
  }

  __host__ __device__ CudaMat<T, Rows, Cols> &operator*=(T scalar)
  {
#pragma unroll
    for (int i = 0; i < Rows * Cols; ++i) {
      dataPtr[i] *= scalar;
    }
    return *this;
  }

  __host__ __device__ CudaMat<T, Rows, Cols> &operator/=(T scalar)
  {
#pragma unroll
    for (int i = 0; i < Rows * Cols; ++i) {
      dataPtr[i] /= scalar;
    }
    return *this;
  }

  __host__ __device__ void setZero()
  {
#pragma unroll
    for (int i = 0; i < Rows * Cols; ++i) {
      dataPtr[i] = T(0);
    }
  }

  __host__ __device__ void setIdentity()
  {
    setZero();
#pragma unroll
    for (int i = 0; i < Rows; ++i) {
      for (int j = 0; j < Cols; ++j) {
        if (i == j) {
          (*this)(i, j) = T(1);
        }
      }
    }
  }

  __host__ __device__ T *data()
  {
    return dataPtr;
  }

  __host__ __device__ const T *data() const
  {
    return dataPtr;
  }

  template<int Rows2, int Cols2>
  __host__ __device__ CudaMat<T, Rows2, Cols2> block(int rowOffset, int colOffset) const
  {
    static_assert(Rows2 <= Rows && Cols2 <= Cols, "Block size exceeds matrix dimensions");

    CudaMat<T, Rows2, Cols2> block;

#pragma unroll
    for (int i = 0; i < Rows2; ++i) {
#pragma unroll
      for (int j = 0; j < Cols2; ++j) {
        int srcRow = rowOffset + i;
        int srcCol = colOffset + j;
        block(i, j) = (*this)(srcRow, srcCol);
      }
    }

    return block;
  }

private:
  T dataPtr[Rows * Cols];
};

}  // namespace CudaUtilities
}  // namespace pgo