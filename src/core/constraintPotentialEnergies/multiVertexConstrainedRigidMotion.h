/*
author: Bohan Wang
copyright to USC, MIT
*/

#pragma once

#include "potentialEnergy.h"

namespace pgo
{
namespace ConstraintPotentialEnergies
{
class MultipleVertexConstrainedRigidMotion : public NonlinearOptimization::PotentialEnergy
{
public:
  MultipleVertexConstrainedRigidMotion(const EigenSupport::VXd &restPositions, const std::vector<int> &vertexIndices, const double *vertexWeights,
    double coeffR, double coefft, bool flexibleRigidMotion, bool isDisp);

  void setDOFs(const std::vector<int> &dofs);

  virtual double func(EigenSupport::ConstRefVecXd u) const override;
  virtual void gradient(EigenSupport::ConstRefVecXd u, EigenSupport::RefVecXd grad) const override;
  virtual void hessian(EigenSupport::ConstRefVecXd u, EigenSupport::SpMatD &hess) const override;

  void setCoeff(double r, double t) { coeff[0] = r, coeff[1] = t; }
  void setCoeff(const double *v) { vertexWeights.assign(v, v + vertexIndices.size()); }

  void setCenter(const double t[3]) { targetCenter = EigenSupport::Mp<const EigenSupport::V3d>(t); }
  void setRotation(const double R[9]) { targetRotation = EigenSupport::Mp<const EigenSupport::M3d>(R); }

  void setRestCenter(const double origin[3] = nullptr);
  EigenSupport::V3d getRestCenter() const { return centerRest; }

  virtual void createHessian(EigenSupport::SpMatD &hess) const override { hess = hessTemplate; }
  virtual void getDOFs(std::vector<int> &dofs) const override { dofs = allDOFs; }
  virtual int getNumDOFs() const override { return (int)allDOFs.size(); }

protected:
  EigenSupport::V3d computePosition(EigenSupport::ConstRefVecXd x, int vid) const;

  const EigenSupport::VXd &restPositions;
  std::vector<int> vertexIndices;
  std::unordered_map<int, int> vertexIDtoIndices;
  std::vector<double> vertexWeights;

  double coeff[2];
  bool flexibleRigidMotion = false;
  bool isDisp = true;

  EigenSupport::V3d targetCenter;
  EigenSupport::M3d targetRotation;

  EigenSupport::V3d centerRest;
  EigenSupport::EigenArray<EigenSupport::V3d> q;

  EigenSupport::SpMatD Z, ZTZ;
  std::vector<EigenSupport::SpMatD> Mi;

  EigenSupport::SpMatD hessTemplate;
  std::vector<int> allDOFs;
  EigenSupport::EntryMap entryMap;
};

}  // namespace ConstraintPotentialEnergies
}  // namespace pgo