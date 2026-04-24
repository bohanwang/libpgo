#pragma once

#include "EigenSupport.h"
#include "potentialEnergy.h"

#include <memory>
#include <atomic>
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

  virtual double computeMaxStepSize(EigenSupport::ConstRefVecXd x, EigenSupport::ConstRefVecXd dx) const override;
  virtual void resetSolveMaxStepStats() const override;
  virtual void recordLineSearchStepDiagnostics(
    double feasibleAlpha,
    double lineSearchAlpha,
    double effectiveAlpha) const override;
  virtual void getFeasibleAlphaClampBreakdown(
    double &materialAlpha,
    double &contactAlpha) const override;
  double getMinFeasibleAlphaThisSolve() const { return minFeasibleAlphaThisSolve_.load(std::memory_order_relaxed); }
  double getMinLineSearchAlphaThisSolve() const { return minLineSearchAlphaThisSolve_.load(std::memory_order_relaxed); }
  double getMinEffectiveAlphaThisSolve() const { return minEffectiveAlphaThisSolve_.load(std::memory_order_relaxed); }

  virtual int isHessianTopologyFixed() const override;
  virtual void hessianDirect(EigenSupport::ConstRefVecXd x, EigenSupport::SpMatD &hess) const override;

  void printImplicitEnergy(EigenSupport::ConstRefVecXd x) const;

protected:
  ImplicitBackwardEulerTimeIntegrator *intg;
  bool approxHessian = true;
  mutable std::atomic<double> currentMaterialFeasibleAlpha_{1.0};
  mutable std::atomic<double> currentContactFeasibleAlpha_{1.0};
  mutable std::atomic<double> minFeasibleAlphaThisSolve_{1.0};
  mutable std::atomic<double> minLineSearchAlphaThisSolve_{1.0};
  mutable std::atomic<double> minEffectiveAlphaThisSolve_{1.0};
};

}  // namespace OptimizationBasedIntegrator
}  // namespace VegaFEM
