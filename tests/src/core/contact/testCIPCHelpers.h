#pragma once

#include "EigenDef.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace pgo::Contact::CIPCTest
{

namespace ES = pgo::EigenSupport;

inline std::tuple<ES::MXd, ES::MXi> makeTwoTriangleMesh()
{
  ES::MXd V(6, 3);
  V <<
    0.00, 0.00, 0.00,
    1.00, 0.00, 0.00,
    0.00, 1.00, 0.00,
    0.23, 0.17, 0.045,
    1.15, 0.28, 0.068,
    0.19, 1.11, 0.026;

  ES::MXi F(2, 3);
  F <<
    0, 1, 2,
    3, 4, 5;

  return { V, F };
}

inline ES::VXd flattenPositions(const ES::MXd &V)
{
  ES::VXd x(V.rows() * 3);
  for (int i = 0; i < V.rows(); ++i)
    x.segment<3>(3 * i) = V.row(i).transpose();
  return x;
}

inline ES::VXd finiteDifferenceGradient(
  const std::function<double(const ES::VXd &)> &energyFn,
  const ES::VXd &x,
  double h)
{
  ES::VXd grad(x.size());
  for (int i = 0; i < x.size(); ++i) {
    ES::VXd xp = x;
    ES::VXd xm = x;
    xp[i] += h;
    xm[i] -= h;
    grad[i] = (energyFn(xp) - energyFn(xm)) / (2.0 * h);
  }
  return grad;
}

inline ES::MXd finiteDifferenceHessian(
  const std::function<ES::VXd(const ES::VXd &)> &gradientFn,
  const ES::VXd &x,
  double h)
{
  ES::MXd H(x.size(), x.size());
  for (int i = 0; i < x.size(); ++i) {
    ES::VXd xp = x;
    ES::VXd xm = x;
    xp[i] += h;
    xm[i] -= h;
    const ES::VXd gp = gradientFn(xp);
    const ES::VXd gm = gradientFn(xm);
    H.col(i) = (gp - gm) / (2.0 * h);
  }
  return H;
}

inline double relativeError(const ES::VXd &a, const ES::VXd &b)
{
  const double denom = std::max({1.0, a.norm(), b.norm()});
  return (a - b).norm() / denom;
}

inline double relativeError(const ES::MXd &A, const ES::MXd &B)
{
  const double denom = std::max({1.0, A.norm(), B.norm()});
  return (A - B).norm() / denom;
}

inline double computeFloorEnergy(const ES::VXd &x, double floorHeight, double floorKappa, int floorAxis = 2)
{
  if (floorAxis < 0 || floorAxis > 2)
    throw std::invalid_argument("floorAxis must be 0, 1, or 2.");
  double energy = 0.0;
  for (int vi = 0; vi < x.size() / 3; ++vi) {
    const double dz = x[3 * vi + floorAxis] - floorHeight;
    if (dz < 0.0)
      energy += 0.5 * floorKappa * dz * dz;
  }
  return energy;
}

inline ES::VXd computeFloorGradient(const ES::VXd &x, double floorHeight, double floorKappa, int floorAxis = 2)
{
  if (floorAxis < 0 || floorAxis > 2)
    throw std::invalid_argument("floorAxis must be 0, 1, or 2.");
  ES::VXd g = ES::VXd::Zero(x.size());
  for (int vi = 0; vi < x.size() / 3; ++vi) {
    const double dz = x[3 * vi + floorAxis] - floorHeight;
    if (dz < 0.0)
      g[3 * vi + floorAxis] = floorKappa * dz;
  }
  return g;
}

inline ES::MXd computeFloorHessian(const ES::VXd &x, double floorHeight, double floorKappa, int floorAxis = 2)
{
  if (floorAxis < 0 || floorAxis > 2)
    throw std::invalid_argument("floorAxis must be 0, 1, or 2.");
  ES::MXd H = ES::MXd::Zero(x.size(), x.size());
  for (int vi = 0; vi < x.size() / 3; ++vi) {
    const double dz = x[3 * vi + floorAxis] - floorHeight;
    if (dz < 0.0)
      H(3 * vi + floorAxis, 3 * vi + floorAxis) = floorKappa;
  }
  return H;
}

inline ES::MXd sparseToDense(const ES::SpMatD &H)
{
  return ES::MXd(H);
}

}  // namespace pgo::Contact::CIPCTest
