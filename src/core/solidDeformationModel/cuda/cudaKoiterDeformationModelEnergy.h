#pragma once

#include "cudaPotentialEnergy.h"
#include "EigenDef.h"

#include <memory>

namespace pgo::SolidDeformationModel
{
class SimulationMesh;

// CPU Koiter-StVK six-vertex patches with fixed rest fundamental forms.
// Parameters are [membrane E, membrane nu, bending E, bending nu, thickness]
// for every element, in the same order as DeformationModelEnergy::setElasticParams.
// Input is displacement. Projected Hessians clamp the full 18x18 element matrix.
class CudaKoiterDeformationModelEnergy final : public NonlinearOptimization::CudaPotentialEnergy
{
public:
  CudaKoiterDeformationModelEnergy(const SimulationMesh &mesh,
    EigenSupport::ConstRefVecXd elasticParameters, bool zeroRestCurvature = false);
  ~CudaKoiterDeformationModelEnergy() override;
  CudaKoiterDeformationModelEnergy(const CudaKoiterDeformationModelEnergy &) = delete;
  CudaKoiterDeformationModelEnergy &operator=(const CudaKoiterDeformationModelEnergy &) = delete;

  int getNumDOFs() const override;
  double func(const double *deviceX) override;
  EnergyChange energyChange(const double *deviceBase, const double *deviceTrial) override;
  void addGradient(const double *deviceX, double *deviceGradient) override;
  void prepareHessian(const double *deviceX) override;
  void prepareExactHessian(const double *deviceX) override;
  void addHessianVector(const double *deviceV, double *deviceOut) const override;
  void addDiagonal(double *deviceDiagonal) const override;
  // CPU Koiter currently supplies no material feasibility bound.
  double computeMaxStepSize(const double *, const double *) const override { return 1.0; }

private:
  void prepareHessian(const double *deviceX, bool project);
  class Impl;
  std::unique_ptr<Impl> data;
};
}  // namespace pgo::SolidDeformationModel
