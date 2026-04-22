/*
copyright to Bohan Wang
*/

#include "embeddedSurfaceIPCPotentialEnergy.h"

#include "scopedProfileSection.h"
#include "surfaceIPCProfiling.h"

#include <numeric>
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
  surfaceIPCCore_(ipcParams),
  surfaceFromSimulationDispMap_(surfaceFromSimulationDispMap)
{
  if (surfaceRestVertices.cols() != 3)
    throw std::invalid_argument("surfaceRestVertices must be an N x 3 matrix.");
  if (surfaceRestVertices.rows() <= 0)
    throw std::invalid_argument("surfaceRestVertices must contain at least one vertex.");
  if (surfaceFromSimulationDispMap_.rows() != surfaceRestVertices.rows() * 3)
    throw std::invalid_argument("surfaceFromSimulationDispMap row count must equal 3 * numSurfaceVertices.");
  if (surfaceFromSimulationDispMap_.cols() <= 0)
    throw std::invalid_argument("surfaceFromSimulationDispMap must contain at least one simulation DOF.");
  if (surfaceTriangles.cols() != 3)
    throw std::invalid_argument("surfaceTriangles must be an M x 3 triangle index matrix.");
  if (surfaceTriangles.size() > 0) {
    if (surfaceTriangles.minCoeff() < 0 || surfaceTriangles.maxCoeff() >= surfaceRestVertices.rows()) {
      throw std::invalid_argument("surfaceTriangles contains an out-of-range vertex index.");
    }
  }

  surfaceRestPositions_.resize(surfaceRestVertices.rows() * 3);
  for (int vi = 0; vi < surfaceRestVertices.rows(); ++vi)
    surfaceRestPositions_.segment<3>(3 * vi) = surfaceRestVertices.row(vi).transpose();

  simulationDOFs_.resize(static_cast<std::size_t>(surfaceFromSimulationDispMap_.cols()));
  std::iota(simulationDOFs_.begin(), simulationDOFs_.end(), 0);

  surfaceIPCCore_.setMesh(surfaceRestVertices, surfaceTriangles);
}

void EmbeddedSurfaceIPCPotentialEnergy::validateSimulationDisplacementSize(EigenSupport::ConstRefVecXd simulationDisplacements) const
{
  if (simulationDisplacements.size() != surfaceFromSimulationDispMap_.cols())
    throw std::invalid_argument("Simulation displacement vector size does not match surfaceFromSimulationDispMap column count.");
}

VXd EmbeddedSurfaceIPCPotentialEnergy::computeSurfaceDisplacementsFromSimulationDisplacements(
  EigenSupport::ConstRefVecXd simulationDisplacements) const
{
  validateSimulationDisplacementSize(simulationDisplacements);
  Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kAdapterMapToSurface);
  return VXd(surfaceFromSimulationDispMap_ * simulationDisplacements);
}

VXd EmbeddedSurfaceIPCPotentialEnergy::computeSurfacePositionsFromSimulationDisplacements(
  EigenSupport::ConstRefVecXd simulationDisplacements) const
{
  return surfaceRestPositions_ + computeSurfaceDisplacementsFromSimulationDisplacements(simulationDisplacements);
}

double EmbeddedSurfaceIPCPotentialEnergy::func(EigenSupport::ConstRefVecXd simulationDisplacements) const
{
  Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kAdapterFunc);
  const VXd surfacePositions = computeSurfacePositionsFromSimulationDisplacements(simulationDisplacements);
  return surfaceIPCCore_.computeEnergy(surfacePositions);
}

void EmbeddedSurfaceIPCPotentialEnergy::gradient(
  EigenSupport::ConstRefVecXd simulationDisplacements,
  EigenSupport::RefVecXd simulationGradient) const
{
  Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kAdapterGradient);
  const VXd surfacePositions = computeSurfacePositionsFromSimulationDisplacements(simulationDisplacements);

  VXd surfaceGradient(surfaceRestPositions_.size());
  surfaceIPCCore_.computeGradient(surfacePositions, surfaceGradient);

  {
    Profiling::ScopedProfileSection pullbackProfile(SurfaceIPCProfileSections::kAdapterPullbackGradient);
    simulationGradient = surfaceFromSimulationDispMap_.transpose() * surfaceGradient;
  }
}

void EmbeddedSurfaceIPCPotentialEnergy::hessian(EigenSupport::ConstRefVecXd, EigenSupport::SpMatD &) const
{
  throw std::runtime_error("EmbeddedSurfaceIPCPotentialEnergy::hessian() should not be called directly. Use hessianDirect() instead.");
}

void EmbeddedSurfaceIPCPotentialEnergy::createHessian(EigenSupport::SpMatD &) const
{
  throw std::runtime_error("EmbeddedSurfaceIPCPotentialEnergy::createHessian() should not be called directly. Use hessianDirect() instead.");
}

void EmbeddedSurfaceIPCPotentialEnergy::hessianDirect(
  EigenSupport::ConstRefVecXd simulationDisplacements,
  EigenSupport::SpMatD &simulationHessian) const
{
  Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kAdapterHessianDirect);
  const VXd surfacePositions = computeSurfacePositionsFromSimulationDisplacements(simulationDisplacements);

  SpMatD surfaceHessian;
  surfaceIPCCore_.computeHessian(surfacePositions, surfaceHessian);

  {
    Profiling::ScopedProfileSection pullbackProfile(SurfaceIPCProfileSections::kAdapterPullbackHessian);
    simulationHessian = surfaceFromSimulationDispMap_.transpose() * surfaceHessian * surfaceFromSimulationDispMap_;
  }
}

double EmbeddedSurfaceIPCPotentialEnergy::computeMaxStepSize(
  EigenSupport::ConstRefVecXd simulationDisplacements,
  EigenSupport::ConstRefVecXd trialSimulationDisplacements) const
{
  Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kAdapterMaxStep);
  validateSimulationDisplacementSize(trialSimulationDisplacements);
  const VXd surfacePositions = computeSurfacePositionsFromSimulationDisplacements(simulationDisplacements);
  const VXd trialSurfaceDisplacements = computeSurfaceDisplacementsFromSimulationDisplacements(trialSimulationDisplacements);
  return surfaceIPCCore_.computeMaxStepSize(surfacePositions, trialSurfaceDisplacements);
}

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
