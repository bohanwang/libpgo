/*
copyright to Bohan Wang
*/

#pragma once

#include "mappedSurfacePotentialEnergy.h"
#include "ipc/core/surfaceIPCCore.h"
#include "ipc/external/obstacleSurface.h"

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
    const SurfaceIPCCore::Parameters &ipcParams = {},
    std::vector<ObstacleSurface> obstacleSurfaces = {});

  void setObstacleTime(double t);
  void markObstacleStatic(int32_t objectId);

private:
  void cacheEnergyActiveSet(SurfaceIPCActiveSet activeSet) const;
  const SurfaceIPCActiveSet *cachedEnergyActiveSetFor(EigenSupport::ConstRefVecXd surfacePositions) const;
  void clearCachedEnergyActiveSet() const;

  virtual double computeSurfaceEnergy(EigenSupport::ConstRefVecXd surfacePositions) const override;
  virtual void computeSurfaceGradient(
    EigenSupport::ConstRefVecXd surfacePositions,
    EigenSupport::RefVecXd surfaceGradient) const override;
  virtual void computeSurfaceHessian(
    EigenSupport::ConstRefVecXd surfacePositions,
    EigenSupport::SpMatD &surfaceHessian) const override;
  virtual void computeSurfaceGradHessian(
    EigenSupport::ConstRefVecXd surfacePositions,
    EigenSupport::RefVecXd surfaceGradient,
    EigenSupport::SpMatD &surfaceHessian) const override;
  virtual void computeSurfaceFuncGrad(
    EigenSupport::ConstRefVecXd surfacePositions,
    double &surfaceEnergy,
    EigenSupport::RefVecXd surfaceGradient) const override;
  virtual void computeSurfaceAll(
    EigenSupport::ConstRefVecXd surfacePositions,
    double &surfaceEnergy,
    EigenSupport::RefVecXd surfaceGradient,
    EigenSupport::SpMatD &surfaceHessian) const override;
  virtual NonlinearOptimization::MaxStepResult computeSurfaceMaxStepLimit(
    EigenSupport::ConstRefVecXd surfacePositions,
    EigenSupport::ConstRefVecXd surfaceDisplacements) const override;
  virtual void beginSurfaceLineSearch(
    EigenSupport::ConstRefVecXd surfacePositions,
    EigenSupport::ConstRefVecXd surfaceDisplacements) const override;
  virtual void endSurfaceLineSearch() const override;

  SurfaceIPCCore surfaceIPCCore_;
  mutable bool hasCachedEnergyActiveSet_ = false;
  mutable SurfaceIPCActiveSet cachedEnergyActiveSet_;
  mutable bool hasLineSearchActiveSet_ = false;
  mutable bool hasLineSearchEnergyState_ = false;
  mutable SurfaceIPCActiveSet lineSearchActiveSet_;
};

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
