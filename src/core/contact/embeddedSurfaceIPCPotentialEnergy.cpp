/*
copyright to Bohan Wang
*/

#include "embeddedSurfaceIPCPotentialEnergy.h"

#include <stdexcept>

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

void EmbeddedSurfaceIPCPotentialEnergy::ensurePreparedForSurfacePositions(
  EigenSupport::ConstRefVecXd surfacePositions) const
{
  if (!surfaceIPCCore_.preparedState().isPreparedFor(surfacePositions))
    surfaceIPCCore_.prepareForSurfacePositions(surfacePositions);
}

double EmbeddedSurfaceIPCPotentialEnergy::computeSurfaceEnergy(EigenSupport::ConstRefVecXd surfacePositions) const
{
  ensurePreparedForSurfacePositions(surfacePositions);
  return surfaceIPCCore_.computeEnergyWithPreparedPairs();
}

void EmbeddedSurfaceIPCPotentialEnergy::computeSurfaceGradient(
  EigenSupport::ConstRefVecXd surfacePositions,
  EigenSupport::RefVecXd surfaceGradient) const
{
  ensurePreparedForSurfacePositions(surfacePositions);
  surfaceIPCCore_.computeGradientWithPreparedPairs(surfaceGradient);
}

void EmbeddedSurfaceIPCPotentialEnergy::computeSurfaceHessian(
  EigenSupport::ConstRefVecXd surfacePositions,
  EigenSupport::SpMatD &surfaceHessian) const
{
  ensurePreparedForSurfacePositions(surfacePositions);
  surfaceIPCCore_.computeHessianWithPreparedPairs(surfaceHessian);
}

NonlinearOptimization::MaxStepResult EmbeddedSurfaceIPCPotentialEnergy::computeSurfaceMaxStepLimit(
  EigenSupport::ConstRefVecXd surfacePositions,
  EigenSupport::ConstRefVecXd surfaceDisplacements) const
{
  return surfaceIPCCore_.computeMaxStepLimit(surfacePositions, surfaceDisplacements);
}

void EmbeddedSurfaceIPCPotentialEnergy::updateObstacleStage(double tStart, double tEnd)
{
  surfaceIPCCore_.updateObstacleStage(tStart, tEnd);
}

void EmbeddedSurfaceIPCPotentialEnergy::invalidatePreparedState()
{
  surfaceIPCCore_.invalidatePreparedState();
}

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
