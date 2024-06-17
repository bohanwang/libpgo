/*
author: Bohan Wang
copyright to USC,MIT
*/

#pragma once

#include "constraintFunctions.h"

#include <tbb/spin_mutex.h>

#include <memory>
#include <vector>

class TetMesh;

namespace pgo
{
namespace SolidDeformationModel
{
class SegmentChainConstraintFunctions : public NonlinearOptimization::ConstraintFunctions
{
public:
  SegmentChainConstraintFunctions(int nAll, int numPoints, const double *const points, int isCircle = 0);
  virtual ~SegmentChainConstraintFunctions();

  virtual void func(EigenSupport::ConstRefVecXd x, EigenSupport::RefVecXd g) const override;
  virtual void jacobian(EigenSupport::ConstRefVecXd x, EigenSupport::SpMatD &jac) const override;
  virtual void hessian(EigenSupport::ConstRefVecXd x, EigenSupport::ConstRefVecXd lambda, EigenSupport::SpMatD &hess) const override;

  virtual bool isLinear() const override { return false; }
  virtual bool isQuadratic() const override { return false; }
  virtual bool hasHessianVector() const override { return false; }

protected:
  int nele;
  int n;

  EigenSupport::VXd restPositions;
  EigenSupport::VXd restLengths;

  typedef Eigen::Matrix<EigenSupport::IDX, 2, 1> JacIndex;
  typedef Eigen::Matrix<EigenSupport::IDX, 6, 6> HessIndex;

  EigenSupport::EigenArray<JacIndex> jacobianIndices;
  EigenSupport::EigenArray<HessIndex> hessIndices;

  mutable std::vector<tbb::spin_mutex> hessLocks;
};

}  // namespace SolidDeformationModel
}  // namespace pgo