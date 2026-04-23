/*
copyright to Bohan Wang
*/

#pragma once

#include "surfaceIPCCore.h"
#include "potentialEnergy.h"

#include <vector>

namespace pgo
{
namespace Contact
{
namespace CIPC
{

using namespace pgo::EigenSupport;

class CIPCPotentialEnergy : public NonlinearOptimization::PotentialEnergy
{
public:
  CIPCPotentialEnergy(double dhat_, double kappa_, bool isInputDisp_,
    double eps_ee_ = 0.0,
    bool useFloor_ = false, double floorHeight_ = -1e-4, double floorKappa_ = 0.1):
    dhat(dhat_), kappa(kappa_), eps_ee(eps_ee_),
    useFloor(useFloor_), floorHeight(floorHeight_), floorKappa(floorKappa_),
    isInputDisp(isInputDisp_) {}

  double dhat = 1e-1;
  double kappa = 0.1;
  double eps_ee = 0.0;

  bool useFloor = false;
  double floorHeight = -1e-4;
  double floorKappa = 0.1;

  double slackness = 1.0;

  void setMesh(const MXd &V, const MXi &F);

  virtual double func(EigenSupport::ConstRefVecXd x) const override;
  virtual void gradient(EigenSupport::ConstRefVecXd x, EigenSupport::RefVecXd grad) const override;
  virtual void hessian(EigenSupport::ConstRefVecXd x, EigenSupport::SpMatD &hess) const override;
  virtual void createHessian(EigenSupport::SpMatD &hess) const override;
  virtual void getDOFs(std::vector<int> &dofs) const override { dofs = allDOFs_; }
  virtual int getNumDOFs() const override { return static_cast<int>(allDOFs_.size()); }
  virtual double computeMaxStepSize(EigenSupport::ConstRefVecXd x, EigenSupport::ConstRefVecXd dx) const override;
  std::int64_t getContactClampCount() const { return core.getContactClampCount(); }
  double getMinContactFeasibleAlphaThisSolve() const { return core.getMinContactFeasibleAlphaThisSolve(); }
  void resetContactMaxStepStats() const { core.resetContactMaxStepStats(); }

  virtual int isHessianTopologyFixed() const override { return 0; }
  virtual void hessianDirect(EigenSupport::ConstRefVecXd x, EigenSupport::SpMatD &hess) const override;

  const std::vector<PTPair> &getPTPairs() const { return core.getPTPairs(); }
  const std::vector<EEPair> &getEEPairs() const { return core.getEEPairs(); }

private:
  void syncCoreParametersFromWrapper() const;
  VXd toSurfacePositions(EigenSupport::ConstRefVecXd x) const;
  VXd toSurfaceDisplacements(EigenSupport::ConstRefVecXd dx) const;
  double computeFloorEnergy(const VXd &x_surf) const;
  void addFloorGradient(const VXd &x_surf, EigenSupport::RefVecXd grad) const;
  void addFloorHessian(const VXd &x_surf, EigenSupport::SpMatD &hess) const;

  mutable SurfaceIPCCore core;
  bool isInputDisp;
  VXd restPosition;
  std::vector<int> allDOFs_;
};

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
