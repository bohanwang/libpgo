/*
author: Bohan Wang
copyright to USC,MIT,NUS
*/

#pragma once

#include "plasticModel3DDeformationGradient.h"

namespace pgo
{
namespace SolidDeformationModel
{
class PlasticModel3DConstant : public PlasticModel3DDeformationGradient
{
public:
  PlasticModel3DConstant(const double Fp_[9]);
  ~PlasticModel3DConstant() {}

  virtual int getNumParameters() const override { return 0; }
  virtual void computeA(const double *param, double A[9]) const override;
  virtual void computeAInv(const double *param, double AInv[9]) const override;
  virtual double compute_detA(const double *param) const override { return detFp; }

protected:
  double Fp[9], FpInv[9], detFp;
};

}  // namespace SolidDeformationModel
}  // namespace pgo
