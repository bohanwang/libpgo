/*
author: Bohan Wang
copyright to USC, MIT
*/

#pragma once

#include "potentialEnergy.h"

#include <numeric>
#include <vector>

namespace pgo
{
namespace PredefinedPotentialEnergies
{
class LinearPotentialEnergy : public NonlinearOptimization::PotentialEnergy
{
public:
  // b^T x
  LinearPotentialEnergy(const EigenSupport::VXd &b_);

  void setDOFs(const std::vector<int> &dofs);

  virtual double func(EigenSupport::ConstRefVecXd x) const override { return x.dot(b); }
  virtual void gradient(EigenSupport::ConstRefVecXd, EigenSupport::RefVecXd grad) const override { grad = b; }
  virtual void hessian(EigenSupport::ConstRefVecXd, EigenSupport::SpMatD &) const override {}

  virtual void createHessian(EigenSupport::SpMatD &hess) const override { hess = EigenSupport::SpMatD(); }
  virtual void getDOFs(std::vector<int> &dofs) const override { dofs = allDOFs; }
  virtual int getNumDOFs() const override { return (int)allDOFs.size(); }

  virtual int isQuadratic() const override { return 0; }
  virtual int hasHessianVector() const override { return 0; }

protected:
  std::vector<int> allDOFs;
  // const EigenSupport::VXd &b;
  EigenSupport::VXd b;
};
}  // namespace NonlinearOptimization
}  // namespace pgo