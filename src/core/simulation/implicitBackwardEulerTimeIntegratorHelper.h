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
  virtual double func_grad_hessian(EigenSupport::ConstRefVecXd x, EigenSupport::RefVecXd grad, EigenSupport::SpMatD &hess) const override;
  virtual void gradient_hessian(EigenSupport::ConstRefVecXd x, EigenSupport::RefVecXd grad, EigenSupport::SpMatD &hess) const override;
  virtual void createHessian(EigenSupport::SpMatD &hess) const override;
  virtual void getDOFs(std::vector<int> &dofs) const override;
  virtual int getNumDOFs() const override;

  virtual NonlinearOptimization::MaxStepResult computeMaxStepLimit(EigenSupport::ConstRefVecXd x, EigenSupport::ConstRefVecXd dx) const override;

  virtual int isHessianTopologyFixed() const override;
  virtual void hessianDirect(EigenSupport::ConstRefVecXd x, EigenSupport::SpMatD &hess) const override;

  void clearCachedImplicitEnergyComponents() const;
  void printImplicitEnergy(EigenSupport::ConstRefVecXd x, bool allowCachedComponents = false) const;

protected:
  void cacheImplicitEnergyComponents(double mainEnergy, const std::vector<double> &componentEnergies) const;
  bool hasCachedImplicitEnergyComponents() const;

  ImplicitBackwardEulerTimeIntegrator *intg;
  bool approxHessian = true;
  mutable bool hasCachedEnergyComponents = false;
  mutable double cachedMainEnergy = 0.0;
  mutable std::vector<double> cachedImplicitModelEnergies;
};

}  // namespace OptimizationBasedIntegrator
}  // namespace VegaFEM
