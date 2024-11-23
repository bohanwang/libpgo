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
class SegmentBinormalConstraintFunctions : public NonlinearOptimization::ConstraintFunctions
{
public:
  SegmentBinormalConstraintFunctions(int nAll, int positionDOFStart, int segmentDOFStart, int numSegments, const double *positions = nullptr);
  virtual ~SegmentBinormalConstraintFunctions();

  virtual void func(EigenSupport::ConstRefVecXd x, EigenSupport::RefVecXd g) const override;
  virtual void jacobian(EigenSupport::ConstRefVecXd x, EigenSupport::SpMatD &jac) const override;
  virtual void hessian(EigenSupport::ConstRefVecXd x, EigenSupport::ConstRefVecXd lambda, EigenSupport::SpMatD &hess) const override;

  virtual bool isLinear() const override { return false; }
  virtual bool isQuadratic() const override { return false; }
  virtual bool hasHessianVector() const override { return false; }

protected:
  int nSeg;
  int positionDOFStart;
  int segmentDOFStart;

  EigenSupport::VXd restPositions;

  typedef Eigen::Matrix<EigenSupport::IDX, 3, 1> JacIndex;
  typedef Eigen::Matrix<EigenSupport::IDX, 9, 9> HessIndex;

  EigenSupport::EigenArray<JacIndex> jacobianIndices;
  EigenSupport::EigenArray<HessIndex> hessIndices;

  mutable std::vector<tbb::spin_mutex> hessLocks;
};

}  // namespace SolidDeformationModel
}  // namespace pgo