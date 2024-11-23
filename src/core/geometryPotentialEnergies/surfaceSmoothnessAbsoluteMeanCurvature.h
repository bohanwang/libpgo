#pragma once

#include "potentialEnergy.h"
#include "triMeshGeo.h"

#include <array>

namespace pgo
{
namespace PredefinedPotentialEnergies
{
struct SurfaceSmoothnessAbsoluteMeanCurvatureBuf;

class SurfaceSmoothnessAbsoluteMeanCurvature : public NonlinearOptimization::PotentialEnergy
{
public:
  SurfaceSmoothnessAbsoluteMeanCurvature(const EigenSupport::VXd &restPositions, const Mesh::TriMeshGeo &mesh, double eps = 1e-8);
  virtual ~SurfaceSmoothnessAbsoluteMeanCurvature();

  void setDOFs(const std::vector<int> &dofs);
  void updateRestInfo();

  virtual double func(EigenSupport::ConstRefVecXd x) const override;
  virtual void gradient(EigenSupport::ConstRefVecXd x, EigenSupport::RefVecXd grad) const override;
  virtual void hessian(EigenSupport::ConstRefVecXd x, EigenSupport::SpMatD &hess) const override;

  virtual void createHessian(EigenSupport::SpMatD &hess) const override { hess = hessTemplate; }
  virtual void getDOFs(std::vector<int> &dofs) const override { dofs = allDOFs; }
  virtual int getNumDOFs() const override { return (int)allDOFs.size(); }

protected:
  const EigenSupport::VXd &restPositions;
  const Mesh::TriMeshGeo &mesh;

  std::vector<std::array<int, 4>> surfaceQuads;
  std::vector<Eigen::Matrix<double, 3, 12>> elementVertexWeights;
  EigenSupport::VXd elementWeights;
  EigenSupport::VXd restValues;

  EigenSupport::SpMatD hessTemplate;
  std::vector<int> allDOFs;

  std::vector<Eigen::Matrix<std::ptrdiff_t, 12, 12>> elementOffsets;

  double eps = 1e-8;

  std::shared_ptr<SurfaceSmoothnessAbsoluteMeanCurvatureBuf> buf;
};

}  // namespace PredefinedPotentialEnergies
}  // namespace pgo