#pragma once

#include "timeIntegrator.h"

namespace pgo
{
namespace Simulation
{
class ImplicitBackwardEulerEnergy;
class TimeIntegratorSolver;

class ImplicitBackwardEulerTimeIntegrator : public TimeIntegrator
{
public:
  ImplicitBackwardEulerTimeIntegrator(const EigenSupport::SpMatD &massMatrix, 
    std::shared_ptr<const NonlinearOptimization::PotentialEnergy> elasticPotential,
    double massDampingCoeff, double stiffnessDampingCoeff, double timestep, int numIter = 0, double eps = 1e-5);

  virtual void doTimestep(int updateq = 1, int verbose = 0, int printResidual = 0) override;
  virtual std::shared_ptr<const NonlinearOptimization::PotentialEnergy> getInternalEnergy() const override;

  const EigenSupport::VXd &getLastSolution() const { return z; }
  void setSolution(EigenSupport::ConstRefVecXd newz);

protected:
  virtual void finiteDifferenceTestIntegratorEnergy(EigenSupport::ConstRefVecXd x) const override;

  void updateA();
  void updateD();
  void updateb();

  // internal buffers
  EigenSupport::SpMatD A, D;
  EigenSupport::VXd b;

  EigenSupport::VXd z;
  EigenSupport::VXd qz, qz1, qz2, temp0, zero;

  friend class ImplicitBackwardEulerEnergy;

  std::shared_ptr<ImplicitBackwardEulerEnergy> eulerEnergy;
  std::shared_ptr<TimeIntegratorSolver> solver;
};

}  // namespace OptimizationBasedIntegrator
}  // namespace VegaFEM
