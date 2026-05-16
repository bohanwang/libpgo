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

double EmbeddedSurfaceIPCPotentialEnergy::computeSurfaceEnergy(EigenSupport::ConstRefVecXd surfacePositions) const
{
  return surfaceIPCCore_.computeEnergy(surfacePositions);
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

void EmbeddedSurfaceIPCPotentialEnergy::computeSurfaceFuncGrad(
  EigenSupport::ConstRefVecXd surfacePositions,
  double &surfaceEnergy,
  EigenSupport::RefVecXd surfaceGradient) const
{
  const SurfaceIPCActiveSet activeSet = surfaceIPCCore_.buildActiveSet(surfacePositions);
  surfaceEnergy = surfaceIPCCore_.computeEnergy(activeSet);
  surfaceIPCCore_.computeGradient(activeSet, surfaceGradient);
}

void EmbeddedSurfaceIPCPotentialEnergy::computeSurfaceAll(
  EigenSupport::ConstRefVecXd surfacePositions,
  double &surfaceEnergy,
  EigenSupport::RefVecXd surfaceGradient,
  EigenSupport::SpMatD &surfaceHessian) const
{
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

void EmbeddedSurfaceIPCPotentialEnergy::updateObstacleStage(double tStart, double tEnd)
{
  surfaceIPCCore_.updateObstacleStage(tStart, tEnd);
}

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
