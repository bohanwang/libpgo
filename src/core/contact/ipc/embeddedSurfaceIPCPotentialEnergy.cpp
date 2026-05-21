/*
copyright to Bohan Wang
*/

#include "ipc/embeddedSurfaceIPCPotentialEnergy.h"

#include "ipc/profiling/surfaceIPCProfiling.h"
#include "scopedProfileSection.h"

#include <stdexcept>
#include <utility>

namespace pgo
{
namespace Contact
{
namespace CIPC
{

EmbeddedSurfaceIPCPotentialEnergy::EmbeddedSurfaceIPCPotentialEnergy(
  const EigenSupport::MXd &surfaceRestVertices,
  const EigenSupport::MXi &surfaceTriangles,
  const EigenSupport::SpMatD &surfaceFromSimulationDispMap,
  const SurfaceIPCCore::Parameters &ipcParams,
  std::vector<ObstacleSurface> obstacleSurfaces):
  MappedSurfacePotentialEnergy(surfaceRestVertices, surfaceFromSimulationDispMap),
  surfaceIPCCore_(ipcParams, std::move(obstacleSurfaces))
{
  if (surfaceTriangles.cols() != 3)
    throw std::invalid_argument("surfaceTriangles must be an M x 3 triangle index matrix.");
  if (surfaceTriangles.size() > 0) {
    if (surfaceTriangles.minCoeff() < 0 || surfaceTriangles.maxCoeff() >= surfaceRestVertices.rows()) {
      throw std::invalid_argument("surfaceTriangles contains an out-of-range vertex index.");
    }
  }

  surfaceIPCCore_.setMesh(surfaceRestVertices, surfaceTriangles);
}

void EmbeddedSurfaceIPCPotentialEnergy::cacheEnergyActiveSet(SurfaceIPCActiveSet activeSet) const
{
  cachedEnergyActiveSet_ = std::move(activeSet);
  hasCachedEnergyActiveSet_ = true;
}

const SurfaceIPCActiveSet *EmbeddedSurfaceIPCPotentialEnergy::cachedEnergyActiveSetFor(
  EigenSupport::ConstRefVecXd surfacePositions) const
{
  if (!hasCachedEnergyActiveSet_)
    return nullptr;
  if (cachedEnergyActiveSet_.positions.size() != surfacePositions.size())
    return nullptr;
  if (!(cachedEnergyActiveSet_.positions.array() == surfacePositions.array()).all())
    return nullptr;
  return &cachedEnergyActiveSet_;
}

void EmbeddedSurfaceIPCPotentialEnergy::clearCachedEnergyActiveSet() const
{
  cachedEnergyActiveSet_.clear();
  hasCachedEnergyActiveSet_ = false;
}

double EmbeddedSurfaceIPCPotentialEnergy::computeSurfaceEnergy(EigenSupport::ConstRefVecXd surfacePositions) const
{
  Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kEnergy);
  clearCachedEnergyActiveSet();
  if (hasLineSearchActiveSet_) {
    lineSearchActiveSet_.positions = surfacePositions;
    hasLineSearchEnergyState_ = true;
    return surfaceIPCCore_.computeEnergy(lineSearchActiveSet_);
  }

  SurfaceIPCActiveSet activeSet = surfaceIPCCore_.buildActiveSet(surfacePositions);
  const double energy = surfaceIPCCore_.computeEnergy(activeSet);
  cacheEnergyActiveSet(std::move(activeSet));
  return energy;
}

void EmbeddedSurfaceIPCPotentialEnergy::computeSurfaceGradient(
  EigenSupport::ConstRefVecXd surfacePositions,
  EigenSupport::RefVecXd surfaceGradient) const
{
  surfaceIPCCore_.computeGradient(surfacePositions, surfaceGradient);
}

void EmbeddedSurfaceIPCPotentialEnergy::computeSurfaceHessian(
  EigenSupport::ConstRefVecXd surfacePositions,
  EigenSupport::SpMatD &surfaceHessian) const
{
  surfaceIPCCore_.computeHessian(surfacePositions, surfaceHessian);
}

void EmbeddedSurfaceIPCPotentialEnergy::computeSurfaceGradHessian(
  EigenSupport::ConstRefVecXd surfacePositions,
  EigenSupport::RefVecXd surfaceGradient,
  EigenSupport::SpMatD &surfaceHessian) const
{
  if (const SurfaceIPCActiveSet *cachedActiveSet = cachedEnergyActiveSetFor(surfacePositions)) {
    surfaceIPCCore_.computeGradient(*cachedActiveSet, surfaceGradient);
    surfaceIPCCore_.computeHessian(*cachedActiveSet, surfaceHessian);
    clearCachedEnergyActiveSet();
    return;
  }

  const SurfaceIPCActiveSet activeSet = surfaceIPCCore_.buildActiveSet(surfacePositions);
  surfaceIPCCore_.computeGradient(activeSet, surfaceGradient);
  surfaceIPCCore_.computeHessian(activeSet, surfaceHessian);
}

void EmbeddedSurfaceIPCPotentialEnergy::computeSurfaceFuncGrad(
  EigenSupport::ConstRefVecXd surfacePositions,
  double &surfaceEnergy,
  EigenSupport::RefVecXd surfaceGradient) const
{
  if (const SurfaceIPCActiveSet *cachedActiveSet = cachedEnergyActiveSetFor(surfacePositions)) {
    surfaceEnergy = surfaceIPCCore_.computeEnergy(*cachedActiveSet);
    surfaceIPCCore_.computeGradient(*cachedActiveSet, surfaceGradient);
    clearCachedEnergyActiveSet();
    return;
  }

  SurfaceIPCActiveSet activeSet = surfaceIPCCore_.buildActiveSet(surfacePositions);
  surfaceEnergy = surfaceIPCCore_.computeEnergy(activeSet);
  surfaceIPCCore_.computeGradient(activeSet, surfaceGradient);
  cacheEnergyActiveSet(std::move(activeSet));
}

void EmbeddedSurfaceIPCPotentialEnergy::computeSurfaceAll(
  EigenSupport::ConstRefVecXd surfacePositions,
  double &surfaceEnergy,
  EigenSupport::RefVecXd surfaceGradient,
  EigenSupport::SpMatD &surfaceHessian) const
{
  if (const SurfaceIPCActiveSet *cachedActiveSet = cachedEnergyActiveSetFor(surfacePositions)) {
    EigenSupport::VXd localGradient = EigenSupport::VXd::Zero(surfaceGradient.size());
    surfaceIPCCore_.computeAll(*cachedActiveSet, surfaceEnergy, localGradient, surfaceHessian);
    surfaceGradient = localGradient;
    clearCachedEnergyActiveSet();
    return;
  }

  const SurfaceIPCActiveSet activeSet = surfaceIPCCore_.buildActiveSet(surfacePositions);
  EigenSupport::VXd localGradient = EigenSupport::VXd::Zero(surfaceGradient.size());
  surfaceIPCCore_.computeAll(activeSet, surfaceEnergy, localGradient, surfaceHessian);
  surfaceGradient = localGradient;
}

NonlinearOptimization::MaxStepResult EmbeddedSurfaceIPCPotentialEnergy::computeSurfaceMaxStepLimit(
  EigenSupport::ConstRefVecXd surfacePositions,
  EigenSupport::ConstRefVecXd surfaceDisplacements) const
{
  return surfaceIPCCore_.computeMaxStepLimit(surfacePositions, surfaceDisplacements);
}

void EmbeddedSurfaceIPCPotentialEnergy::beginSurfaceLineSearch(
  EigenSupport::ConstRefVecXd surfacePositions,
  EigenSupport::ConstRefVecXd surfaceDisplacements) const
{
  clearCachedEnergyActiveSet();
  lineSearchActiveSet_ = surfaceIPCCore_.buildLineSearchActiveSetSuperset(surfacePositions, surfaceDisplacements);
  hasLineSearchActiveSet_ = true;
  hasLineSearchEnergyState_ = false;
}

void EmbeddedSurfaceIPCPotentialEnergy::endSurfaceLineSearch() const
{
  if (hasLineSearchActiveSet_ && hasLineSearchEnergyState_)
    cacheEnergyActiveSet(lineSearchActiveSet_);
  lineSearchActiveSet_.clear();
  hasLineSearchActiveSet_ = false;
  hasLineSearchEnergyState_ = false;
}

void EmbeddedSurfaceIPCPotentialEnergy::setObstacleTime(double t)
{
  clearCachedEnergyActiveSet();
  lineSearchActiveSet_.clear();
  hasLineSearchActiveSet_ = false;
  hasLineSearchEnergyState_ = false;
  surfaceIPCCore_.setObstacleTime(t);
}

void EmbeddedSurfaceIPCPotentialEnergy::markObstacleStatic(int32_t objectId)
{
  clearCachedEnergyActiveSet();
  lineSearchActiveSet_.clear();
  hasLineSearchActiveSet_ = false;
  hasLineSearchEnergyState_ = false;
  surfaceIPCCore_.markObstacleStatic(objectId);
}

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
