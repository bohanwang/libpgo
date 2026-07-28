/*
copyright to Bohan Wang
*/

#pragma once

#include "mappedSurfacePotentialEnergy.h"
#include "ipc/core/surfaceIPCCore.h"

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
