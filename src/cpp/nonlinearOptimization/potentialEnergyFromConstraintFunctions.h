/*
author: Bohan Wang
copyright to USC
*/

#pragma once

#include "potentialEnergy.h"
#include "constraintFunctions.h"

#include <memory>

namespace pgo
{
namespace ES = EigenSupport;

namespace NonlinearOptimization
{
class PotentialEnergyConstraintFunctions : public PotentialEnergy
{
public:
  PotentialEnergyConstraintFunctions(int nAll, std::shared_ptr<const ConstraintFunctions> cnstt);
  virtual ~PotentialEnergyConstraintFunctions();

  virtual double func(EigenSupport::ConstRefVecXd x) const override;
  virtual void gradient(EigenSupport::ConstRefVecXd x, EigenSupport::RefVecXd grad) const override;
  virtual void hessian(EigenSupport::ConstRefVecXd x, EigenSupport::SpMatD &hess) const override;

  virtual void createHessian(EigenSupport::SpMatD &hess) const override { hess = hessAll; }
  virtual void getDOFs(std::vector<int> &dofs) const override { dofs = allDOFs; }
  virtual int getNumDOFs() const override { return (int)allDOFs.size(); }

protected:
  std::shared_ptr<const ConstraintFunctions> cnstt;
  std::vector<int> allDOFs;

  mutable EigenSupport::VXd g;
  mutable EigenSupport::SpMatD jac, lambdaHessian, JTJ;
  EigenSupport::SpMatD hessAll;
  EigenSupport::SymbolicMmData *mmData;
  EigenSupport::SpMatI JTJMapping, lambdaHessMapping;
};

}  // namespace NonlinearOptimization
}  // namespace pgo
