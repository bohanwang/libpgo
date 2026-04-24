/*
copyright to Bohan Wang
*/

#pragma once

#include "potentialEnergy.h"

#include <vector>

namespace pgo
{
namespace Contact
{
namespace CIPC
{

using namespace pgo::EigenSupport;

class MappedSurfacePotentialEnergy : public NonlinearOptimization::PotentialEnergy
{
public:
  MappedSurfacePotentialEnergy(
    const EigenSupport::MXd &surfaceRestVertices,
    const EigenSupport::SpMatD &surfaceFromSimulationDispMap);

  virtual double func(EigenSupport::ConstRefVecXd simulationDisplacements) const override;
  virtual void gradient(
    EigenSupport::ConstRefVecXd simulationDisplacements,
    EigenSupport::RefVecXd simulationGradient) const override;
  virtual void hessian(EigenSupport::ConstRefVecXd, EigenSupport::SpMatD &) const override;
  virtual void createHessian(EigenSupport::SpMatD &) const override;
  virtual void hessianDirect(
    EigenSupport::ConstRefVecXd simulationDisplacements,
    EigenSupport::SpMatD &simulationHessian) const override;
  virtual double computeMaxStepSize(
    EigenSupport::ConstRefVecXd simulationDisplacements,
    EigenSupport::ConstRefVecXd trialSimulationDisplacements) const override;

  virtual void getDOFs(std::vector<int> &dofs) const override { dofs = simulationDOFs_; }
  virtual int getNumDOFs() const override { return static_cast<int>(simulationDOFs_.size()); }
  virtual int isHessianTopologyFixed() const override { return 0; }

protected:
  void validateSimulationDisplacementSize(EigenSupport::ConstRefVecXd simulationDisplacements) const;
  VXd computeSurfaceDisplacementsFromSimulationDisplacements(EigenSupport::ConstRefVecXd simulationDisplacements) const;
  VXd computeSurfacePositionsFromSimulationDisplacements(EigenSupport::ConstRefVecXd simulationDisplacements) const;

  virtual double computeSurfaceEnergy(EigenSupport::ConstRefVecXd surfacePositions) const = 0;
  virtual void computeSurfaceGradient(
    EigenSupport::ConstRefVecXd surfacePositions,
    EigenSupport::RefVecXd surfaceGradient) const = 0;
  virtual void computeSurfaceHessian(
    EigenSupport::ConstRefVecXd surfacePositions,
    EigenSupport::SpMatD &surfaceHessian) const = 0;
  virtual double computeSurfaceMaxStepSize(
    EigenSupport::ConstRefVecXd surfacePositions,
    EigenSupport::ConstRefVecXd surfaceDisplacements) const;

private:
  VXd surfaceRestPositions_;
  SpMatD surfaceFromSimulationDispMap_;
  std::vector<int> simulationDOFs_;
};

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
