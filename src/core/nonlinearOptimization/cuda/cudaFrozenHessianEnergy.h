#pragma once

#include "cudaPotentialEnergy.h"
#include "EigenDef.h"

#include <memory>

namespace pgo::NonlinearOptimization
{

// E(x) = scale/2 * (x-reference)^T H(reference) (x-reference).
// A separate model owns the frozen device Hessian. Do not prepare that model
// elsewhere while this energy is in use. This represents Rayleigh stiffness
// damping with scale = beta / timestep; Newton never updates H(reference).
class CudaFrozenHessianEnergy final : public CudaPotentialEnergy
{
public:
  explicit CudaFrozenHessianEnergy(std::shared_ptr<CudaPotentialEnergy> model);
  ~CudaFrozenHessianEnergy() override;
  CudaFrozenHessianEnergy(const CudaFrozenHessianEnergy &) = delete;
  CudaFrozenHessianEnergy &operator=(const CudaFrozenHessianEnergy &) = delete;

  void setReference(EigenSupport::ConstRefVecXd reference, double scale);
  int getNumDOFs() const override;
  double func(const double *deviceX) override;
  EnergyChange energyChange(const double *deviceBase, const double *deviceTrial) override;
  void addGradient(const double *deviceX, double *deviceGradient) override;
  void prepareHessian(const double *deviceX) override;
  void addHessianVector(const double *deviceV, double *deviceOut) const override;
  void addDiagonal(double *deviceDiagonal) const override;
  double computeMaxStepSize(const double *, const double *) const override { return 1.0; }

private:
  class Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace pgo::NonlinearOptimization
