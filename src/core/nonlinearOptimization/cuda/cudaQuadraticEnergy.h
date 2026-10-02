#pragma once

#include "cudaPotentialEnergy.h"
#include "EigenDef.h"

#include <memory>

namespace pgo::NonlinearOptimization
{

// E(x) = 0.5 (x-reference)^T A (x-reference)
//        + linear^T (x-reference) + constant.
// The symmetric sparse A and all vectors are uploaded once. This supports
// consistent mass matrices, fixed damping, external forces and pulling terms.
class CudaQuadraticEnergy : public CudaPotentialEnergy
{
public:
  CudaQuadraticEnergy(const EigenSupport::SpMatD &A,
    EigenSupport::ConstRefVecXd reference, EigenSupport::ConstRefVecXd linear, double constant = 0);
  ~CudaQuadraticEnergy();
  CudaQuadraticEnergy(const CudaQuadraticEnergy &) = delete;
  CudaQuadraticEnergy &operator=(const CudaQuadraticEnergy &) = delete;

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
