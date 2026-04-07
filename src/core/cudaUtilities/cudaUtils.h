#pragma once

#include "pgoLogging.h"

#include <cuda_runtime.h>

#include <array>
#include <cstring>

#define IF_CUDA_ERR(err, jump)                                                                                        \
  do {                                                                                                                \
    if (err != cudaSuccess) {                                                                                         \
      spdlog::error("CUDA error at {}:{}: code: {}; msg: {}", __FILE__, __LINE__, int(err), cudaGetErrorString(err)); \
      jump;                                                                                                           \
    }                                                                                                                 \
  } while (0)

#define CUBLAS_CHECK(func, jump)                                                                                     \
  {                                                                                                                  \
    cublasStatus_t status = (func);                                                                                  \
    if (status != CUBLAS_STATUS_SUCCESS) {                                                                           \
      spdlog::error("CUBLAS API failed at {}:{} with error: {}", __FILE__, __LINE__, cublasGetStatusString(status)); \
      jump;                                                                                                          \
    }                                                                                                                \
  }

#define CUSPARSE_CHECK(func, jump)                                                                             \
  do {                                                                                                         \
    cusparseStatus_t err = (func);                                                                             \
    if (err != CUSPARSE_STATUS_SUCCESS) {                                                                      \
      spdlog::error("cuSPARSE error at {}:{} with error {}", __FILE__, __LINE__, cusparseGetErrorString(err)); \
      jump;                                                                                                    \
    }                                                                                                          \
  } while (0)

#define CUSOLVER_CHECK(func, jump)                                                          \
  do {                                                                                      \
    cusolverStatus_t err = (func);                                                          \
    if (err != CUSOLVER_STATUS_SUCCESS) {                                                   \
      spdlog::error("cuSPARSE error at {}:{} with error {}", __FILE__, __LINE__, int(err)); \
      jump;                                                                                 \
    }                                                                                       \
  } while (0)

#define CUDSS_CHECK(func, jump)                                                          \
  do {                                                                                   \
    cudssStatus_t err = (func);                                                          \
    if (err != CUDSS_STATUS_SUCCESS) {                                                   \
      spdlog::error("CUDSS error at {}:{} with error {}", __FILE__, __LINE__, int(err)); \
      jump;                                                                              \
    }                                                                                    \
  } while (0)

#define IF_CHECK(call, jump)                                            \
  do {                                                                  \
    int ret = call;                                                     \
    if (ret != 0) {                                                     \
      spdlog::error("Error at {}:{}: {}\n", __FILE__, __LINE__, #call); \
      jump;                                                             \
    }                                                                   \
  } while (0)

template<typename T, int N>
struct ArrayCmp
{
  bool operator()(const std::array<T, N> &f1, const std::array<T, N> &f2) const
  {
    return std::memcmp(f1.data(), f2.data(), sizeof(f1)) < 0;
  }
};

__host__ __device__ inline int elt(int row, int col, int numRows)
{
  return col * numRows + row;
}

#if defined(_MSC_VER)
__host__ __device__ inline float arg(float a)
{
  return a;
}
#endif
