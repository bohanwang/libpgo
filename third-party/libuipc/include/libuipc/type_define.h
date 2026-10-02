#pragma once
// libpgo compatibility bridge for the selected libuipc CUDA math routines.
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cassert>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <limits>
#ifdef __CUDACC__
#define UIPC_GENERIC __host__ __device__
#else
#define UIPC_GENERIC
#endif
#define UIPC_INLINE inline
#define UIPC_KERNEL_ASSERT(condition, ...) assert(condition)
#define UIPC_KERNEL_ERROR_WITH_LOCATION(...) assert(false)
namespace uipc {
using Float = double;
using IndexT = int32_t;
using Eigen::Matrix;
using Eigen::Vector;
using Vector2i = Vector<IndexT, 2>;
using Vector3i = Vector<IndexT, 3>;
using Vector4i = Vector<IndexT, 4>;
using Vector2 = Vector<double, 2>;
using Vector3 = Vector<double, 3>;
using Vector4 = Vector<double, 4>;
using Matrix3x3 = Matrix<double, 3, 3>;
}
