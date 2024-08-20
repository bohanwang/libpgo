#pragma once

#include "potentialEnergy.h"
#include "triMeshGeo.h"

namespace pgo
{
namespace PredefinedPotentialEnergies
{
struct VertexAffineToPositionEnergyBuf;

class VertexAffineToPositionEnergy : public pgo::NonlinearOptimization::PotentialEnergy
{
public:
  VertexAffineToPositionEnergy(int totalNumDOFs, const pgo::EigenSupport::VXd &restPosition, std::shared_ptr<const PotentialEnergy> originalEnergy);
  virtual ~VertexAffineToPositionEnergy();

  void computeAToX();
  void compute_x_orig(pgo::EigenSupport::ConstRefVecXd x, pgo::EigenSupport::RefVecXd x_orig) const;

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
  const pgo::EigenSupport::VXd &restPosition;
  pgo::EigenSupport::SpMatD AToX;

  std::shared_ptr<const PotentialEnergy> originalEnergy;
  pgo::EigenSupport::SpMatD hessTemplate;
  std::vector<int> allDOFs;

  std::shared_ptr<VertexAffineToPositionEnergyBuf> buf;
};
}  // namespace PredefinedPotentialEnergies
}  // namespace pgo
