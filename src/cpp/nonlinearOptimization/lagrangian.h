/*
author: Bohan Wang
copyright to USC
*/

#pragma once

#include "potentialEnergy.h"
#include "constraintFunctions.h"

namespace pgo
{
namespace NonlinearOptimization
{
class Lagrangian : public PotentialEnergy
{
public:
  Lagrangian(PotentialEnergy_const_p energy, ConstraintFunctions_const_p constraints = nullptr);
  virtual ~Lagrangian() {}

  virtual double func(EigenSupport::ConstRefVecXd x) const override;
  virtual void gradient(EigenSupport::ConstRefVecXd x, EigenSupport::RefVecXd grad) const override;
  virtual void hessian(EigenSupport::ConstRefVecXd x, EigenSupport::SpMatD &hess) const override;
  virtual void createHessian(EigenSupport::SpMatD &h) const override { h = hess; }

  virtual void getDOFs(std::vector<int> &dofs) const override { dofs = allDOFs; }
  virtual int getNumDOFs() const override { return (int)allDOFs.size(); }

protected:
  std::vector<int> allDOFs;
  PotentialEnergy_const_p energy;
  ConstraintFunctions_const_p constraints;

  mutable EigenSupport::VXd grad, g;
  mutable EigenSupport::SpMatD hess, energyHessian, constraintHessian, constraintJac, constraintJacT;
  EigenSupport::SpMatI cJacMappingBL, cJacMappingTR, cHessMapping, cJacTMapping, eMapping;
};

}  // namespace NonlinearOptimization
}  // namespace pgo