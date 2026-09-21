/*
copyright to Bohan Wang
*/

#pragma once

#include "mappedSurfacePotentialEnergy.h"
#include "ipc/core/surfaceIPCCore.h"

#include <cstdint>
#include <vector>

namespace pgo
{
namespace Contact
{
namespace CIPC
{

using namespace pgo::EigenSupport;

class EmbeddedSurfaceIPCPotentialEnergy : public MappedSurfacePotentialEnergy
{
public:
  EmbeddedSurfaceIPCPotentialEnergy(
    const EigenSupport::MXd &surfaceRestVertices,
    const EigenSupport::MXi &surfaceTriangles,
    const EigenSupport::SpMatD &surfaceFromSimulationDispMap,
    const SurfaceIPCCore::Parameters &ipcParams = {});
  EmbeddedSurfaceIPCPotentialEnergy(
    const EigenSupport::MXd &surfaceRestVertices,
    const EigenSupport::MXi &surfaceTriangles,
    const EigenSupport::SpMatD &surfaceFromSimulationDispMap,
    const std::vector<uint8_t> &vertexIsDeformableMask,
    const SurfaceIPCCore::Parameters &ipcParams = {});

  int getNumSurfaceVertices() const { return surfaceIPCCore_.getNumSurfaceVertices(); }
  int getNumSurfaceTriangles() const { return surfaceIPCCore_.getNumSurfaceTriangles(); }
  bool isSurfaceVertexDeformable(int vi) const { return surfaceIPCCore_.isSurfaceVertexDeformable(vi); }
  void validateCollisionFreeState(EigenSupport::ConstRefVecXd simulationDisplacements) const;
  void updateFriction(EigenSupport::ConstRefVecXd referenceDisplacements,
    EigenSupport::ConstRefVecXd laggedDisplacements, double timestep);
  double getFrictionCoeff() const { return surfaceIPCCore_.getParameters().frictionCoeff; }
  std::size_t getNumFrictionPairs() const { return surfaceIPCCore_.getNumFrictionPairs(); }

private:
  virtual double computeSurfaceEnergy(EigenSupport::ConstRefVecXd surfacePositions) const override;
  virtual void computeSurfaceGradient(
    EigenSupport::ConstRefVecXd surfacePositions,
    EigenSupport::RefVecXd surfaceGradient) const override;
  virtual void computeSurfaceHessian(
    EigenSupport::ConstRefVecXd surfacePositions,
    EigenSupport::SpMatD &surfaceHessian) const override;
  virtual double computeSurfaceMaxStepSize(
    EigenSupport::ConstRefVecXd surfacePositions,
    EigenSupport::ConstRefVecXd surfaceDisplacements) const override;

  void ensurePreparedForSurfacePositions(EigenSupport::ConstRefVecXd surfacePositions) const;

  SurfaceIPCCore surfaceIPCCore_;
};

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
