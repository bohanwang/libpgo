#pragma once

#include "EigenSupport.h"
#include "potentialEnergy.h"

#include <memory>
#include <vector>

namespace pgo
{
namespace Simulation
{
class ImplicitBackwardEulerTimeIntegrator;

class ImplicitBackwardEulerEnergy : public NonlinearOptimization::PotentialEnergy
{
public:
  ImplicitBackwardEulerEnergy(ImplicitBackwardEulerTimeIntegrator *integrator_);

  virtual double func(EigenSupport::ConstRefVecXd x) const override;
  virtual void gradient(EigenSupport::ConstRefVecXd x, EigenSupport::RefVecXd grad) const override;
  virtual void hessian(EigenSupport::ConstRefVecXd x, EigenSupport::SpMatD &hess) const override;
  virtual void createHessian(EigenSupport::SpMatD &hess) const override;
  virtual void getDOFs(std::vector<int> &dofs) const override;
  virtual int getNumDOFs() const override;

  void printImplicitEnergy(EigenSupport::ConstRefVecXd x) const;

protected:
  ImplicitBackwardEulerTimeIntegrator *intg;
  bool approxHessian = true;
};

}  // namespace OptimizationBasedIntegrator
}  // namespace VegaFEM