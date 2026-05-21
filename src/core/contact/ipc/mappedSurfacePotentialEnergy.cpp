/*
copyright to Bohan Wang
*/

#include "mappedSurfacePotentialEnergy.h"

#include "scopedProfileSection.h"
#include "ipc/profiling/surfaceIPCProfiling.h"

#include <numeric>
#include <stdexcept>

namespace pgo
{
namespace Contact
{
namespace CIPC
{

MappedSurfacePotentialEnergy::MappedSurfacePotentialEnergy(
  const EigenSupport::MXd &surfaceRestVertices,
  const EigenSupport::SpMatD &surfaceFromSimulationDispMap):
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

  surfaceRestPositions_.resize(surfaceRestVertices.rows() * 3);
  for (int vi = 0; vi < surfaceRestVertices.rows(); ++vi)
    surfaceRestPositions_.segment<3>(3 * vi) = surfaceRestVertices.row(vi).transpose();

  simulationDOFs_.resize(static_cast<std::size_t>(surfaceFromSimulationDispMap_.cols()));
  std::iota(simulationDOFs_.begin(), simulationDOFs_.end(), 0);
}

void MappedSurfacePotentialEnergy::validateSimulationDisplacementSize(EigenSupport::ConstRefVecXd simulationDisplacements) const
{
  if (simulationDisplacements.size() != surfaceFromSimulationDispMap_.cols())
    throw std::invalid_argument("Simulation displacement vector size does not match surfaceFromSimulationDispMap column count.");
}

VXd MappedSurfacePotentialEnergy::computeSurfaceDisplacementsFromSimulationDisplacements(
  EigenSupport::ConstRefVecXd simulationDisplacements) const
{
  validateSimulationDisplacementSize(simulationDisplacements);
  Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kAdapterMapToSurface);
  return VXd(surfaceFromSimulationDispMap_ * simulationDisplacements);
}

VXd MappedSurfacePotentialEnergy::computeSurfacePositionsFromSimulationDisplacements(
  EigenSupport::ConstRefVecXd simulationDisplacements) const
{
  return surfaceRestPositions_ + computeSurfaceDisplacementsFromSimulationDisplacements(simulationDisplacements);
}

double MappedSurfacePotentialEnergy::func(EigenSupport::ConstRefVecXd simulationDisplacements) const
{
  Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kAdapterFunc);
  const VXd surfacePositions = computeSurfacePositionsFromSimulationDisplacements(simulationDisplacements);
  return computeSurfaceEnergy(surfacePositions);
}

void MappedSurfacePotentialEnergy::gradient(
  EigenSupport::ConstRefVecXd simulationDisplacements,
  EigenSupport::RefVecXd simulationGradient) const
{
  Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kAdapterGradient);
  const VXd surfacePositions = computeSurfacePositionsFromSimulationDisplacements(simulationDisplacements);

  VXd surfaceGradient = VXd::Zero(surfaceRestPositions_.size());
  computeSurfaceGradient(surfacePositions, surfaceGradient);

  {
    Profiling::ScopedProfileSection pullbackProfile(SurfaceIPCProfileSections::kAdapterPullbackGradient);
    simulationGradient = surfaceFromSimulationDispMap_.transpose() * surfaceGradient;
  }
}

void MappedSurfacePotentialEnergy::computeSurfaceFuncGrad(
  EigenSupport::ConstRefVecXd surfacePositions,
  double &surfaceEnergy,
  EigenSupport::RefVecXd surfaceGradient) const
{
  surfaceEnergy = computeSurfaceEnergy(surfacePositions);
  computeSurfaceGradient(surfacePositions, surfaceGradient);
}

void MappedSurfacePotentialEnergy::computeSurfaceAll(
  EigenSupport::ConstRefVecXd surfacePositions,
  double &surfaceEnergy,
  EigenSupport::RefVecXd surfaceGradient,
  EigenSupport::SpMatD &surfaceHessian) const
{
  surfaceEnergy = computeSurfaceEnergy(surfacePositions);
  computeSurfaceGradient(surfacePositions, surfaceGradient);
  computeSurfaceHessian(surfacePositions, surfaceHessian);
}

void MappedSurfacePotentialEnergy::computeSurfaceGradHessian(
  EigenSupport::ConstRefVecXd surfacePositions,
  EigenSupport::RefVecXd surfaceGradient,
  EigenSupport::SpMatD &surfaceHessian) const
{
  computeSurfaceGradient(surfacePositions, surfaceGradient);
  computeSurfaceHessian(surfacePositions, surfaceHessian);
}

double MappedSurfacePotentialEnergy::func_grad(
  EigenSupport::ConstRefVecXd simulationDisplacements,
  EigenSupport::RefVecXd simulationGradient) const
{
  Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kAdapterGradient);
  const VXd surfacePositions = computeSurfacePositionsFromSimulationDisplacements(simulationDisplacements);

  double surfaceEnergy = 0.0;
  VXd surfaceGradient = VXd::Zero(surfaceRestPositions_.size());
  computeSurfaceFuncGrad(surfacePositions, surfaceEnergy, surfaceGradient);

  {
    Profiling::ScopedProfileSection pullbackProfile(SurfaceIPCProfileSections::kAdapterPullbackGradient);
    simulationGradient = surfaceFromSimulationDispMap_.transpose() * surfaceGradient;
  }
  return surfaceEnergy;
}

double MappedSurfacePotentialEnergy::func_grad_hessian(
  EigenSupport::ConstRefVecXd simulationDisplacements,
  EigenSupport::RefVecXd simulationGradient,
  EigenSupport::SpMatD &simulationHessian) const
{
  Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kAdapterHessianDirect);
  const VXd surfacePositions = computeSurfacePositionsFromSimulationDisplacements(simulationDisplacements);

  double surfaceEnergy = 0.0;
  VXd surfaceGradient = VXd::Zero(surfaceRestPositions_.size());
  SpMatD surfaceHessian(surfaceRestPositions_.size(), surfaceRestPositions_.size());
  computeSurfaceAll(surfacePositions, surfaceEnergy, surfaceGradient, surfaceHessian);

  {
    Profiling::ScopedProfileSection pullbackProfile(SurfaceIPCProfileSections::kAdapterPullbackGradient);
    simulationGradient = surfaceFromSimulationDispMap_.transpose() * surfaceGradient;
  }
  {
    Profiling::ScopedProfileSection pullbackProfile(SurfaceIPCProfileSections::kAdapterPullbackHessian);
    simulationHessian = surfaceFromSimulationDispMap_.transpose() * surfaceHessian * surfaceFromSimulationDispMap_;
  }
  return surfaceEnergy;
}

void MappedSurfacePotentialEnergy::gradient_hessian(
  EigenSupport::ConstRefVecXd simulationDisplacements,
  EigenSupport::RefVecXd simulationGradient,
  EigenSupport::SpMatD &simulationHessian) const
{
  Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kAdapterHessianDirect);
  const VXd surfacePositions = computeSurfacePositionsFromSimulationDisplacements(simulationDisplacements);

  VXd surfaceGradient = VXd::Zero(surfaceRestPositions_.size());
  SpMatD surfaceHessian(surfaceRestPositions_.size(), surfaceRestPositions_.size());
  computeSurfaceGradHessian(surfacePositions, surfaceGradient, surfaceHessian);

  {
    Profiling::ScopedProfileSection pullbackProfile(SurfaceIPCProfileSections::kAdapterPullbackGradient);
    simulationGradient = surfaceFromSimulationDispMap_.transpose() * surfaceGradient;
  }
  {
    Profiling::ScopedProfileSection pullbackProfile(SurfaceIPCProfileSections::kAdapterPullbackHessian);
    simulationHessian = surfaceFromSimulationDispMap_.transpose() * surfaceHessian * surfaceFromSimulationDispMap_;
  }
}

void MappedSurfacePotentialEnergy::hessian(EigenSupport::ConstRefVecXd, EigenSupport::SpMatD &) const
{
  throw std::runtime_error("MappedSurfacePotentialEnergy::hessian() should not be called directly. Use hessianDirect() instead.");
}

void MappedSurfacePotentialEnergy::createHessian(EigenSupport::SpMatD &) const
{
  throw std::runtime_error("MappedSurfacePotentialEnergy::createHessian() should not be called directly. Use hessianDirect() instead.");
}

void MappedSurfacePotentialEnergy::hessianDirect(
  EigenSupport::ConstRefVecXd simulationDisplacements,
  EigenSupport::SpMatD &simulationHessian) const
{
  Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kAdapterHessianDirect);
  const VXd surfacePositions = computeSurfacePositionsFromSimulationDisplacements(simulationDisplacements);

  SpMatD surfaceHessian(surfaceRestPositions_.size(), surfaceRestPositions_.size());
  computeSurfaceHessian(surfacePositions, surfaceHessian);

  {
    Profiling::ScopedProfileSection pullbackProfile(SurfaceIPCProfileSections::kAdapterPullbackHessian);
    simulationHessian = surfaceFromSimulationDispMap_.transpose() * surfaceHessian * surfaceFromSimulationDispMap_;
  }
}

NonlinearOptimization::MaxStepResult MappedSurfacePotentialEnergy::computeMaxStepLimit(
  EigenSupport::ConstRefVecXd simulationDisplacements,
  EigenSupport::ConstRefVecXd trialSimulationDisplacements) const
{
  Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kAdapterMaxStep);
  const VXd surfacePositions = computeSurfacePositionsFromSimulationDisplacements(simulationDisplacements);
  const VXd trialSurfaceDisplacements = computeSurfaceDisplacementsFromSimulationDisplacements(trialSimulationDisplacements);
  return computeSurfaceMaxStepLimit(surfacePositions, trialSurfaceDisplacements);
}

void MappedSurfacePotentialEnergy::beginLineSearch(
  EigenSupport::ConstRefVecXd simulationDisplacements,
  EigenSupport::ConstRefVecXd trialSimulationDisplacements) const
{
  const VXd surfacePositions = computeSurfacePositionsFromSimulationDisplacements(simulationDisplacements);
  const VXd trialSurfaceDisplacements = computeSurfaceDisplacementsFromSimulationDisplacements(trialSimulationDisplacements);
  beginSurfaceLineSearch(surfacePositions, trialSurfaceDisplacements);
}

void MappedSurfacePotentialEnergy::endLineSearch() const
{
  endSurfaceLineSearch();
}

NonlinearOptimization::MaxStepResult MappedSurfacePotentialEnergy::computeSurfaceMaxStepLimit(
  EigenSupport::ConstRefVecXd,
  EigenSupport::ConstRefVecXd) const
{
  return NonlinearOptimization::MaxStepResult::unconstrained();
}

void MappedSurfacePotentialEnergy::beginSurfaceLineSearch(
  EigenSupport::ConstRefVecXd,
  EigenSupport::ConstRefVecXd) const
{
}

void MappedSurfacePotentialEnergy::endSurfaceLineSearch() const
{
}

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
