/*
author: Bohan Wang
copyright to USC,MIT
*/

#pragma once

#include "constraintFunctions.h"
#include "deformationModel.h"

#include <tbb/spin_mutex.h>

#include <memory>
#include <vector>
#include <functional>

namespace pgo
{
namespace SolidDeformationModel
{
class DeformationModelManager;

class PrescribedPrincipleStressConstraintFunctions : public NonlinearOptimization::ConstraintFunctions
{
public:
  PrescribedPrincipleStressConstraintFunctions(int nAll, int dofOffset, int numElements, const int *elementIDs, const DeformationModelManager *tetMeshDMM);
  virtual ~PrescribedPrincipleStressConstraintFunctions() {}

  using XToPosFunc = std::function<void(const EigenSupport::V3d &, int offset, EigenSupport::V3d &)>;
  void setXToPosFunc(XToPosFunc func) { xToPosFunc = func; }
  void setTargetPHat(const double *phat) { targetPrincipleStress = EigenSupport::Mp<const EigenSupport::VXd>(phat, elements.size() * 3); }

  void computeForceFromTargetPHat(EigenSupport::ConstRefVecXd x, EigenSupport::RefVecXd fext) const;
  double computeSurfaceNormalTractionFromElement(EigenSupport::ConstRefVecXd x, const EigenSupport::V3d &n, int eleID) const;

  virtual void func(EigenSupport::ConstRefVecXd x, EigenSupport::RefVecXd g) const override;
  virtual void jacobian(EigenSupport::ConstRefVecXd x, EigenSupport::SpMatD &jac) const override;
  virtual void hessian(EigenSupport::ConstRefVecXd x, EigenSupport::ConstRefVecXd lambda, EigenSupport::SpMatD &hess) const override;

  virtual bool isLinear() const override { return false; }
  virtual bool isQuadratic() const override { return false; }
  virtual bool hasHessianVector() const override { return false; }

protected:
  int dofStart;
  const DeformationModelManager *tetMeshDMM;
  std::vector<int> elements;
  XToPosFunc xToPosFunc;
  EigenSupport::VXd targetPrincipleStress;
  std::vector<DeformationModel::CacheData *> elementCacheData;

  EigenSupport::EntryMap jacEntries, hessEntries;

  mutable std::vector<tbb::spin_mutex> hessLocks;
};

}  // namespace SolidDeformationModel
}  // namespace pgo