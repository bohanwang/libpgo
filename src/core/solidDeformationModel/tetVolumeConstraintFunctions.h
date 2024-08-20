/*
author: Bohan Wang
copyright to USC,MIT
*/

#pragma once

#include "constraintFunctions.h"
#include "simulationMesh.h"

#include <tbb/spin_mutex.h>

#include <memory>
#include <vector>

class TetMesh;

namespace pgo
{
namespace SolidDeformationModel
{
class TetVolumeConstraintFunctions : public NonlinearOptimization::ConstraintFunctions
{
public:
  TetVolumeConstraintFunctions(const SimulationMesh *tetMesh, int nAll, const EigenSupport::VXd *restPosition = nullptr, const EigenSupport::M3Xd *DmInv = nullptr);
  virtual ~TetVolumeConstraintFunctions();

  void setDmInv(const EigenSupport::M3Xd &DmInv_);
  void setElementFlags(const int *flags) { elementFlags.assign(flags, flags + nele); }

  virtual void func(EigenSupport::ConstRefVecXd x, EigenSupport::RefVecXd g) const override;
  virtual void jacobian(EigenSupport::ConstRefVecXd x, EigenSupport::SpMatD &jac) const override;
  virtual void hessian(EigenSupport::ConstRefVecXd x, EigenSupport::ConstRefVecXd lambda, EigenSupport::SpMatD &hess) const override;

  virtual bool isLinear() const override { return false; }
  virtual bool isQuadratic() const override { return false; }
  virtual bool hasHessianVector() const override { return false; }

protected:
  const SimulationMesh *tetMesh;
  const EigenSupport::VXd *restPosition;

  EigenSupport::M3Xd DmInv;
  EigenSupport::EigenArray<EigenSupport::M9x12d> dFdx;

  typedef Eigen::Matrix<EigenSupport::IDX, 4, 1> JacIndex;
  typedef Eigen::Matrix<EigenSupport::IDX, 12, 12> HessIndex;

  EigenSupport::EigenArray<JacIndex> jacobianIndices;
  EigenSupport::EigenArray<HessIndex> hessIndices;
  std::vector<int> elementFlags;

  mutable std::vector<tbb::spin_mutex> hessLocks;

  int nele;
};

}  // namespace SolidDeformationModel
}  // namespace pgo