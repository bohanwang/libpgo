/*
author: Bohan Wang
copyright to USC, MIT
*/

#pragma once

#include "potentialEnergyAligningMeshConnectivity.h"

namespace pgo
{
namespace ConstraintPotentialEnergies
{
class MultipleVertexPulling : public PotentialEnergyAligningMeshConnectivity
{
public:
  MultipleVertexPulling(const EigenSupport::SpMatD &Koff, const double *restPositionsAll,
    int numPts, const int *vertexIndices, const double *tgt, const double *bcCoeff, int isDisp);
  virtual double func(EigenSupport::ConstRefVecXd u) const override;
  virtual void gradient(EigenSupport::ConstRefVecXd u, EigenSupport::RefVecXd grad) const override;
  virtual void hessian(EigenSupport::ConstRefVecXd, EigenSupport::SpMatD &hess) const override;

  void setCoeff(double v) { coeffAll = v; }
  void setCoeff(const double *v);
  void setTargetPos(const double *tgt);

  void printErrorInfo(EigenSupport::ConstRefVecXd u) const;

protected:
  typedef Eigen::Matrix<EigenSupport::IDX, 3, 3> M3i;
  std::vector<M3i, Eigen::aligned_allocator<M3i>> KIndices;

  EigenSupport::VXd tgtp, restpAll;
  double coeffAll = 1.0;
  std::vector<double> coeffs;
  std::vector<int> vertexIndices;

  int isDisp;
};
}  // namespace ConstraintPotentialEnergies
}  // namespace pgo