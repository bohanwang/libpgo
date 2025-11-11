/*
author: Bohan Wang
copyright to USC
*/

#pragma once

#include <memory>
#include <functional>

namespace pgo
{
namespace NonlinearOptimization
{
class PotentialEnergy;
class ConstraintFunctions;

class FiniteDifference
{
public:
  enum METHOD
  {
    M_THREE_POINT,
    M_FIVE_POINT
  };
  FiniteDifference(METHOD m, double eps);

  void testEnergy(std::shared_ptr<const PotentialEnergy> energy, bool testGradient, bool testHessian, double range = -1.0,
    const double *x = nullptr, int dofRange = -1, double *gradRelError = nullptr, double *hessRelError = nullptr) const;

  typedef std::function<void(const double *, double *, double *, double *)> EvalFunc;
  void testEnergy(EvalFunc evalFunc, int n, bool testGradient, bool testHessian, double range = -1.0, const double *x = nullptr,
    double *gradRelError = nullptr, double *hessRelError = nullptr) const;

  void testConstraints(std::shared_ptr<const ConstraintFunctions> c, double range = -1.0, const double *x = nullptr, const double *lambda = nullptr) const;

  typedef std::function<void(const double *, double *, double *)> EvalVecFunc;
  void testVecFunc(EvalVecFunc evalFunc, int m, int n, double range = -1.0, const double *x = nullptr,
  double *err = nullptr) const;

  static void randomSeq(double range, double *x, int n);

  void gradient(const double *x, double *grad, int n, std::vector<int> *dofs, std::function<double(const double *)> energyFunc);
  void hessian(const double *x, double *hess, int n, std::vector<int> *dofs, std::function<void(const double *, double *)> gradFunc);

protected:
  METHOD method;
  double eps;

  int left, right;
  const double *coeffs;
};
}  // namespace NonlinearOptimization
}  // namespace pgo