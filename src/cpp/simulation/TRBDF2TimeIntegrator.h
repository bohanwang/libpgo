#pragma once

#include "timeIntegrator.h"
#include "timeIntegratorSolver.h"

namespace pgo
{

namespace NonlinearOptimization
{
class PotentialEnergy;
}  // namespace NonlinearOptimization

namespace Simulation
{
class TRBDF2TimeIntegratorEnergy;

class TRBDF2TimeIntegrator : public TimeIntegrator
{
public:
  TRBDF2TimeIntegrator(const EigenSupport::SpMatD &massMatrix, 
    std::shared_ptr<const NonlinearOptimization::PotentialEnergy> elasticPotential,
    double massDampingCoeff, double stiffnessDampingCoeff, double gamma, double timestep, int numIter = 0, double eps = 1e-6);

  virtual void doTimestep(int updateq = 1, int verbose = 0, int printResidual = 0) override;
  virtual std::shared_ptr<const NonlinearOptimization::PotentialEnergy> getInternalEnergy() const override;

  const EigenSupport::VXd &getLastSolutionTR() const { return z1; }
  const EigenSupport::VXd &getLastSolutionBDF2() const { return z2; }

  void setGamma(double gamma);
  virtual void setTimestep(double t) override;

protected:
  virtual void finiteDifferenceTestIntegratorEnergy(EigenSupport::ConstRefVecXd x) const override;

  void updateA1();
  void updateA2();

  void updateD();

  void updateb1();
  void updateb2();

  void updateCoeffs();

  void solve(EigenSupport::VXd &x, std::shared_ptr<TRBDF2TimeIntegratorEnergy> eng, int verbose, int printResidual);

  // internal buffers
  EigenSupport::SpMatD A1, A2, D;
  EigenSupport::VXd b1, b2;

  EigenSupport::VXd z1, z2;
  EigenSupport::VXd qz, temp0;
  EigenSupport::VXd qy, qvely, qaccy;
  EigenSupport::VXd zero;

  double alpha;
  double beta[10];
  double y;

  int stage = 0;

  friend class TRBDF2TimeIntegratorEnergy;

  std::shared_ptr<TRBDF2TimeIntegratorEnergy> trEnergy;
  std::shared_ptr<TRBDF2TimeIntegratorEnergy> bdf2Energy;
  std::shared_ptr<TimeIntegratorSolver> solver[2];
};

}  // namespace OptimizationBasedIntegrator
}  // namespace VegaFEM
