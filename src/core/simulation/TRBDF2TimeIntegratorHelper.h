#pragma once

#include "potentialEnergy.h"

#include <memory>
#include <vector>

namespace pgo
{
namespace Simulation
{
class TRBDF2TimeIntegrator;

class TRBDF2TimeIntegratorEnergy : public NonlinearOptimization::PotentialEnergy
{
public:
  TRBDF2TimeIntegratorEnergy(TRBDF2TimeIntegrator *integrator_, const EigenSupport::SpMatD &A, const EigenSupport::VXd &b);

  virtual double func(EigenSupport::ConstRefVecXd x) const override;
  virtual void gradient(EigenSupport::ConstRefVecXd x, EigenSupport::RefVecXd grad) const override;
  virtual void hessian(EigenSupport::ConstRefVecXd x, EigenSupport::SpMatD &hess) const override;
  virtual void gradient_hessian(EigenSupport::ConstRefVecXd x, EigenSupport::RefVecXd grad, EigenSupport::SpMatD &hess) const override;
  virtual void createHessian(EigenSupport::SpMatD &hess) const override;
  virtual void getDOFs(std::vector<int> &dofs) const override;
  virtual int getNumDOFs() const override;

  virtual NonlinearOptimization::MaxStepResult computeMaxStepLimit(EigenSupport::ConstRefVecXd x, EigenSupport::ConstRefVecXd dx) const override;

  virtual int isHessianTopologyFixed() const override;
  virtual void hessianDirect(EigenSupport::ConstRefVecXd x, EigenSupport::SpMatD &hess) const override;

  void printImplicitEnergy(EigenSupport::ConstRefVecXd x) const;

protected:
  TRBDF2TimeIntegrator *intg;
  const EigenSupport::SpMatD &A;
  const EigenSupport::VXd &b;
};
}  // namespace OptimizationBasedIntegrator
}  // namespace VegaFEM
