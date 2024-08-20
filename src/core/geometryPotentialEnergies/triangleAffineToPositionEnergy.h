#pragma once

#include "potentialEnergy.h"
#include "triMeshGeo.h"

namespace pgo
{
namespace PredefinedPotentialEnergies
{
struct TriangleAffineToPositionEnergyBuf;

class TriangleAffineToPositionEnergy : public pgo::NonlinearOptimization::PotentialEnergy
{
public:
  TriangleAffineToPositionEnergy(int totalNumDOFs, const Mesh::TriMeshGeo &mesh, std::shared_ptr<const PotentialEnergy> originalEnergy);
  virtual ~TriangleAffineToPositionEnergy();

  void computeWA();
  void compute_A_orig(pgo::EigenSupport::ConstRefVecXd A_tri, pgo::EigenSupport::RefVecXd A_orig) const;

  virtual double func(pgo::EigenSupport::ConstRefVecXd x) const override;
  virtual void gradient(pgo::EigenSupport::ConstRefVecXd x, pgo::EigenSupport::RefVecXd grad) const override;
  virtual void hessian(pgo::EigenSupport::ConstRefVecXd x, pgo::EigenSupport::SpMatD &hess) const override;

  virtual void createHessian(pgo::EigenSupport::SpMatD &hess) const override { hess = hessTemplate; }
  virtual void getDOFs(std::vector<int> &dofs) const override { dofs = allDOFs; }
  virtual int getNumDOFs() const override { return (int)allDOFs.size(); }
  virtual int isQuadratic() const override { return 0; }
  virtual int hasHessianVector() const override { return 0; }

private:
  int nAll, nRestDOFs = 0;
  const Mesh::TriMeshGeo &inputMesh;

  std::vector<std::vector<std::tuple<int, double>>> vertexNeighboringTriangleWeights;
  std::vector<pgo::EigenSupport::TripletD> entries;
  pgo::EigenSupport::SpMatD WA;

  std::shared_ptr<const PotentialEnergy> originalEnergy;
  pgo::EigenSupport::SpMatD hessTemplate;
  std::vector<int> allDOFs;

  std::shared_ptr<TriangleAffineToPositionEnergyBuf> buf;
};
}  // namespace PredefinedPotentialEnergies
}  // namespace pgo
