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
  const SurfaceIPCCore::Parameters &ipcParams):
  MappedSurfacePotentialEnergy(surfaceRestVertices, surfaceFromSimulationDispMap),
  surfaceIPCCore_(ipcParams)
{
  surfaceIPCCore_.setMesh(surfaceRestVertices, surfaceTriangles);
}

EmbeddedSurfaceIPCPotentialEnergy::EmbeddedSurfaceIPCPotentialEnergy(
  const EigenSupport::MXd &surfaceRestVertices,
  const EigenSupport::MXi &surfaceTriangles,
  const EigenSupport::SpMatD &surfaceFromSimulationDispMap,
  const std::vector<uint8_t> &vertexIsDeformableMask,
  const SurfaceIPCCore::Parameters &ipcParams):
  MappedSurfacePotentialEnergy(surfaceRestVertices, surfaceFromSimulationDispMap),
  surfaceIPCCore_(ipcParams)
{
  if (vertexIsDeformableMask.size() != static_cast<std::size_t>(surfaceRestVertices.rows()))
    throw std::invalid_argument("EmbeddedSurfaceIPCPotentialEnergy deformable-role mask length must equal the number of surface vertices.");

  for (Eigen::Index outer = 0; outer < surfaceFromSimulationDispMap.outerSize(); ++outer) {
    for (EigenSupport::SpMatD::InnerIterator it(surfaceFromSimulationDispMap, outer); it; ++it) {
      const Eigen::Index vertex = it.row() / 3;
      if (vertexIsDeformableMask[static_cast<std::size_t>(vertex)] == 0 && it.value() != 0.0)
        throw std::invalid_argument("EmbeddedSurfaceIPCPotentialEnergy external vertices must have zero displacement-map rows.");
    }
  }

  surfaceIPCCore_.setMesh(surfaceRestVertices, surfaceTriangles, vertexIsDeformableMask);
}

void EmbeddedSurfaceIPCPotentialEnergy::validateCollisionFreeState(
  EigenSupport::ConstRefVecXd simulationDisplacements) const
{
  surfaceIPCCore_.validateCollisionFreeSurfacePositions(
    computeSurfacePositionsFromSimulationDisplacements(simulationDisplacements));
}

void EmbeddedSurfaceIPCPotentialEnergy::ensurePreparedForSurfacePositions(
  EigenSupport::ConstRefVecXd surfacePositions) const
{
  if (!surfaceIPCCore_.isPreparedFor(surfacePositions))
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

double EmbeddedSurfaceIPCPotentialEnergy::computeSurfaceMaxStepSize(
  EigenSupport::ConstRefVecXd surfacePositions,
  EigenSupport::ConstRefVecXd surfaceDisplacements) const
{
  return surfaceIPCCore_.computeMaxStepSize(surfacePositions, surfaceDisplacements);
}

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
