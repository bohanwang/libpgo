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

class EmbeddedSurfaceIPCPotentialEnergy : public NonlinearOptimization::PotentialEnergy
{
public:
  EmbeddedSurfaceIPCPotentialEnergy(
    const EigenSupport::MXd &surfaceRestVertices,
    const EigenSupport::MXi &surfaceTriangles,
    const EigenSupport::SpMatD &surfaceFromSimulationDispMap,
    const SurfaceIPCCore::Parameters &ipcParams = {});

  virtual double func(EigenSupport::ConstRefVecXd simulationDisplacements) const override;
  virtual void gradient(
    EigenSupport::ConstRefVecXd simulationDisplacements,
    EigenSupport::RefVecXd simulationGradient) const override;
  virtual void hessian(EigenSupport::ConstRefVecXd simulationDisplacements, EigenSupport::SpMatD &simulationHessian) const override;
  virtual void createHessian(EigenSupport::SpMatD &simulationHessian) const override;
  virtual void hessianDirect(
    EigenSupport::ConstRefVecXd simulationDisplacements,
    EigenSupport::SpMatD &simulationHessian) const override;
  virtual double computeMaxStepSize(
    EigenSupport::ConstRefVecXd simulationDisplacements,
    EigenSupport::ConstRefVecXd trialSimulationDisplacements) const override;
  std::int64_t getContactClampCount() const { return surfaceIPCCore_.getContactClampCount(); }
  double getMinContactFeasibleAlphaThisSolve() const { return surfaceIPCCore_.getMinContactFeasibleAlphaThisSolve(); }
  void resetContactMaxStepStats() const { surfaceIPCCore_.resetContactMaxStepStats(); }

  virtual void getDOFs(std::vector<int> &dofs) const override { dofs = simulationDOFs_; }
  virtual int getNumDOFs() const override { return static_cast<int>(simulationDOFs_.size()); }
  virtual int isHessianTopologyFixed() const override { return 0; }

private:
  VXd computeSurfaceDisplacementsFromSimulationDisplacements(EigenSupport::ConstRefVecXd simulationDisplacements) const;
  VXd computeSurfacePositionsFromSimulationDisplacements(EigenSupport::ConstRefVecXd simulationDisplacements) const;
  void validateSimulationDisplacementSize(EigenSupport::ConstRefVecXd simulationDisplacements) const;

  SurfaceIPCCore surfaceIPCCore_;
  VXd surfaceRestPositions_;
  SpMatD surfaceFromSimulationDispMap_;
  std::vector<int> simulationDOFs_;
};

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
