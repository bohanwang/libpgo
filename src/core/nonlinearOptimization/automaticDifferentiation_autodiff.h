/*
author: Bohan Wang
copyright to USC, MIT
*/

#pragma once

#include "EigenSupport.h"

#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>

namespace pgo
{
namespace NonlinearOptimization
{
namespace AutomaticDifferentiation_autodiff
{
using dual2nd = autodiff::dual2nd;
using dualV3d = Eigen::Matrix<dual2nd, 3, 1>;

template<int n, typename Func, typename Derived, typename... Args>
double computeEnergy(Func func, const Eigen::MatrixBase<Derived> &x, Args &&...args)
{
  dual2nd dualx[n];
  for (int i = 0; i < n; i++) {
    dualx[i] = x(i);
  }

  dual2nd eng = func(dualx, std::forward<Args>(args)...);
  return eng.val.val;
}

template<typename Func, typename Derived, typename... Args>
double computeEnergy(int n, Func func, const Eigen::MatrixBase<Derived> &x, Args &&...args)
{
  if (n != (int)x.size()) {
    return -1.0;
  }

  std::vector<dual2nd> dualx(x.size());
  for (typename Eigen::MatrixBase<Derived>::Index i = 0; i < x.size(); i++) {
    dualx[i] = x(i);
  }

  dual2nd eng = func(dualx.data(), std::forward<Args>(args)...);
  return eng.val.val;
}

template<int n, typename Func, typename Derived, typename... Args>
void computeGradient(Func func, const Eigen::MatrixBase<Derived> &x, Eigen::Matrix<double, n, 1> &grad, Args &&...args)
{
  dual2nd dualx[n];
  for (int i = 0; i < n; i++) {
    dualx[i] = x(i);
  }

  for (int i = 0; i < n; i++) {
    autodiff::dual1st val = autodiff::derivative(func, autodiff::wrt(dualx[i]), autodiff::at(dualx, std::forward<Args>(args)...));
    grad[i] = val.val;
  }
}

template<typename Func, typename Derived, typename... Args>
void computeGradient(Func func, const Eigen::MatrixBase<Derived> &x, EigenSupport::VXd &grad, Args &&...args)
{
  std::vector<dual2nd> dualx(x.size());
  for (typename Eigen::MatrixBase<Derived>::Index i = 0; i < x.size(); i++) {
    dualx[i] = x(i);
  }

  for (typename Eigen::MatrixBase<Derived>::Index i = 0; i < x.size(); i++) {
    grad[i] = autodiff::derivative(func, autodiff::wrt(dualx[i]), autodiff::at(dualx.data(), std::forward<Args>(args)...));
  }
}

template<int n, typename Func, typename Derived, typename... Args>
void computeHessian(Func func, const Eigen::MatrixBase<Derived> &x, Eigen::Matrix<double, n, n> &h, Args &&...args)
{
  dual2nd dualx[n];
  for (int i = 0; i < n; i++) {
    dualx[i] = x(i);
  }

  for (int i = 0; i < n; i++) {
    for (int j = i; j < n; j++) {
      h(i, j) = autodiff::derivative<2>(func, autodiff::wrt(dualx[i], dualx[j]), autodiff::at(dualx, std::forward<Args>(args)...));

      if (i != j)
        h(j, i) = h(i, j);
    }
  }
}

template<typename Func, typename Derived, typename... Args>
void computeHessian(Func func, const Eigen::MatrixBase<Derived> &x, EigenSupport::MXd &h, Args &&...args)
{
  std::vector<dual2nd> dualx(x.size());
  int n = (int)x.size();
  for (int i = 0; i < n; i++) {
    dualx[i] = x(i);
  }

  for (int i = 0; i < n; i++) {
    for (int j = i; j < n; j++) {
      h(i, j) = autodiff::derivative<2>(func, autodiff::wrt(dualx[i], dualx[j]), autodiff::at(dualx.data(), std::forward<Args>(args)...));

      if (j != i)
        h(j, i) = h(i, j);
    }
  }
}

}  // namespace AutomaticDifferentiation_autodiff
}  // namespace NonlinearOptimization
}  // namespace pgo
