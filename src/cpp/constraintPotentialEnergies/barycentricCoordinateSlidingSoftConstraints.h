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
struct BarycentricCoordinateSlidingBuf;

class BarycentricCoordinateSliding : public PotentialEnergyAligningMeshConnectivity
{
public:
  BarycentricCoordinateSliding(const EigenSupport::SpMatD &Koff, const EigenSupport::VXd &restPositionsAll,
    int numPoints, const double *constraintCoeffs, const double *constraintTargetPositions, const double *constraintNormals,
    const int *idx, const double *barycentricWeights, int numElementVertices, double bcCoeff, int useN, int checkPenetration, int isInputDisplacement);

  virtual double func(EigenSupport::ConstRefVecXd u) const override;
  virtual void gradient(EigenSupport::ConstRefVecXd u, EigenSupport::RefVecXd grad) const override;
  virtual void hessian(EigenSupport::ConstRefVecXd u, EigenSupport::SpMatD &hess) const override;

  void setCoeff(double v) { coeffAll = v; }
  void setCoeff(const double *v);
  void setTargetPos(const double *tgt);
  void setNormals(const double *n);
  void computeHessian();

  void printErrorInfo(EigenSupport::ConstRefVecXd u) const;

protected:
  const EigenSupport::VXd &restpAll;
  EigenSupport::VXd normals;
  EigenSupport::VXd tgtp;
  EigenSupport::VXd coeffs;

  std::vector<int> vertexIndices;
  std::vector<double> vertexBarycentricWeights;
  EigenSupport::EigenArray<EigenSupport::MXi> hessianIndices;
  EigenSupport::SpMatD hessianConstant;

  double coeffAll;
  int numElementVertices;
  int useNormal;
  int checkPenetration;
  int isInputDisplacement;

  std::shared_ptr<BarycentricCoordinateSlidingBuf> buf;
};

}  // namespace ConstraintPotentialEnergies
}  // namespace pgo