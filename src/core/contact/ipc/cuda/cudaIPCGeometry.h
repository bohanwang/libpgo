/*
copyright to Bohan Wang
*/

#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>

namespace pgo::Contact::CIPC::CudaIPCGeometry
{

// Interior line-line distance derivatives through the closest-point envelope.
// Orthogonalizing the two line directions avoids forming ac-b^2 or subtracting
// large nearly equal terms in the inverse of the closest-point normal matrix.
template<typename T>
EIGEN_DEVICE_FUNC void interiorEdgeEdge(const Eigen::Matrix<T, 3, 1> &a,
  const Eigen::Matrix<T, 3, 1> &b, const Eigen::Matrix<T, 3, 1> &c,
  const Eigen::Matrix<T, 3, 1> &d, T &distance,
  Eigen::Matrix<T, 12, 1> *gradient, Eigen::Matrix<T, 12, 12> *hessian)
{
  using V3 = Eigen::Matrix<T, 3, 1>;
  const V3 u = b - a, v = d - c, w = a - c, n = u.cross(v);
  const T nn = n.squaredNorm(), uu = u.squaredNorm();
  const T s = v.cross(w).dot(n) / nn, t = u.cross(w).dot(n) / nn;
  const V3 residual = (w.dot(n) / nn) * n;
  distance = residual.squaredNorm();
  const T weights[4] = { T(1) - s, s, t - T(1), -t };
  if (gradient)
    for (int i = 0; i < 4; ++i)
      gradient->template segment<3>(3 * i) = T(2) * weights[i] * residual;
  if (!hessian)
    return;
  const T ratio = u.dot(v) / uu;
  const V3 perpendicular = n.cross(u) / uu;
  const T ds[4] = { -1, 1, 0, 0 }, dt[4] = { 0, 0, 1, -1 };
  Eigen::Matrix<T, 12, 1> first, second;
  using std::sqrt;
  const T firstScale = sqrt(T(2) / uu), secondScale = sqrt(T(2) * uu / nn);
  for (int i = 0; i < 4; ++i) {
    first.template segment<3>(3 * i) = firstScale * (ds[i] * residual + weights[i] * u);
    second.template segment<3>(3 * i) = secondScale * ((dt[i] + ratio * ds[i]) * residual - weights[i] * perpendicular);
    for (int j = 0; j < 4; ++j)
      hessian->template block<3, 3>(3 * i, 3 * j) = T(2) * weights[i] * weights[j] * Eigen::Matrix<T, 3, 3>::Identity();
  }
  *hessian -= first * first.transpose() + second * second.transpose();
}

}  // namespace pgo::Contact::CIPC::CudaIPCGeometry
