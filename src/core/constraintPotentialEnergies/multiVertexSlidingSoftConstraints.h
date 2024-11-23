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
struct MultipleVertexSlidingBuf;

class MultipleVertexSliding : public PotentialEnergyAligningMeshConnectivity
{
public:
  MultipleVertexSliding(const EigenSupport::SpMatD &Koff, const EigenSupport::VXd &restPositionsAll,
    int numPoints, const double *constraintCoeffs, const double *constraintTargetPositions, const double *constraintNormals,
    const int *idx, double bcCoeff, int useN, int checkPenetration, int isInputDisplacement);

  virtual double func(EigenSupport::ConstRefVecXd u) const override;
  virtual void gradient(EigenSupport::ConstRefVecXd u, EigenSupport::RefVecXd grad) const override;
  virtual void hessian(EigenSupport::ConstRefVecXd u, EigenSupport::SpMatD &hess) const override;

  void setCoeff(double v) { coeffAll = v; }
  void setCoeff(const double *v);
  void setTargetPos(const double *tgt);
  void setNormals(const double *n);
  void setCheckContact(int check) { checkPenetration = check; }
  void computeHessian();

  void printErrorInfo(EigenSupport::ConstRefVecXd u) const;

protected:
  const EigenSupport::VXd &restpAll;
  EigenSupport::VXd normals;
  EigenSupport::VXd tgtp;
  EigenSupport::VXd coeffs;

  std::vector<int> vertexIndices;

  EigenSupport::EigenArray<EigenSupport::M3i> hessianIndices;
  EigenSupport::SpMatD hessianConstant;

  double coeffAll;
  int useNormal;
  int checkPenetration;
  int isDisp = 0;

  std::shared_ptr<MultipleVertexSlidingBuf> buf;
};
}  // namespace ConstraintPotentialEnergies
}  // namespace pgo