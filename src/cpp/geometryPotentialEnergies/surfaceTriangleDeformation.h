#pragma once

#include "potentialEnergy.h"
#include "triMeshGeo.h"

#include <array>

namespace pgo
{
namespace PredefinedPotentialEnergies
{
struct SurfaceTriangleDeformationBuf;

class SurfaceTriangleDeformation : public NonlinearOptimization::PotentialEnergy
{
public:
  SurfaceTriangleDeformation(const EigenSupport::VXd &restPositions, const Mesh::TriMeshGeo &mesh, double eps = 1e-8);
  virtual ~SurfaceTriangleDeformation();

  void setDOFs(const std::vector<int> &dofs);
  void updateRestInfo();

  virtual double func(EigenSupport::ConstRefVecXd x) const override;
  virtual void gradient(EigenSupport::ConstRefVecXd x, EigenSupport::RefVecXd grad) const override;
  virtual void hessian(EigenSupport::ConstRefVecXd x, EigenSupport::SpMatD &hess) const override;

  virtual void createHessian(EigenSupport::SpMatD &hess) const { hess = hessTemplate; }
  virtual void getDOFs(std::vector<int> &dofs) const { dofs = allDOFs; }
  virtual int getNumDOFs() const { return (int)allDOFs.size(); }

protected:
  const EigenSupport::VXd &restPositions;
  const Mesh::TriMeshGeo &mesh;
  double coeffs[2] = { 1.0, 0.0 };

  Eigen::Matrix<double, 6, 9> dFdx;
  EigenSupport::VXd elementWeights;

  EigenSupport::SpMatD hessTemplate;
  std::vector<int> allDOFs;

  std::vector<Eigen::Matrix<std::ptrdiff_t, 9, 9>> elementOffsets;
  std::vector<EigenSupport::M2d> elementAbarInv;

  double eps = 1e-8;

  std::shared_ptr<SurfaceTriangleDeformationBuf> buf;
};

}  // namespace PredefinedPotentialEnergies
}  // namespace pgo