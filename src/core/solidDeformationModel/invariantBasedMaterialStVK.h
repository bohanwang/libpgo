/*
author: Bohan Wang
copyright to USC,MIT
*/

#pragma once

#include "invariantBasedMaterial.h"

namespace pgo
{
namespace SolidDeformationModel
{
class InvariantBasedMaterialStVK : public InvariantBasedMaterial
{
public:
  InvariantBasedMaterialStVK(double E, double nu, double compressionRatio);
  virtual ~InvariantBasedMaterialStVK() {}

  virtual double compute_psi(const double invariants[3]) const override;
  virtual void compute_dpsi_dI(const double invariants[3], double gradient[3]) const override;
  virtual void compute_d2psi_dI2(const double invariants[3], double hessian[6]) const override;

  void setMaterial(double mu_, double lambda_) { this->mu = mu_, this->lambda = lambda_; }

protected:
  double mu, lambda, coeffJ;
};

}  // namespace SolidDeformationModel
}  // namespace pgo