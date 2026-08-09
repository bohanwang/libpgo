/*
copyright to Bohan Wang
*/

#include "mappedSurfacePotentialEnergy.h"

#include <cmath>
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
  if (!surfaceRestVertices.allFinite())
    throw std::invalid_argument("surfaceRestVertices must be finite.");
  if (surfaceFromSimulationDispMap_.rows() != surfaceRestVertices.rows() * 3)
    throw std::invalid_argument("surfaceFromSimulationDispMap row count must equal 3 * numSurfaceVertices.");
  if (surfaceFromSimulationDispMap_.cols() <= 0)
    throw std::invalid_argument("surfaceFromSimulationDispMap must contain at least one simulation DOF.");
  for (Eigen::Index i = 0; i < surfaceFromSimulationDispMap_.nonZeros(); ++i) {
    if (!std::isfinite(surfaceFromSimulationDispMap_.valuePtr()[i]))
      throw std::invalid_argument("surfaceFromSimulationDispMap must contain only finite coefficients.");
  }

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
  return VXd(surfaceFromSimulationDispMap_ * simulationDisplacements);
}

VXd MappedSurfacePotentialEnergy::computeSurfacePositionsFromSimulationDisplacements(
  EigenSupport::ConstRefVecXd simulationDisplacements) const
{
  return surfaceRestPositions_ + computeSurfaceDisplacementsFromSimulationDisplacements(simulationDisplacements);
}

double MappedSurfacePotentialEnergy::func(EigenSupport::ConstRefVecXd simulationDisplacements) const
{
  const VXd surfacePositions = computeSurfacePositionsFromSimulationDisplacements(simulationDisplacements);
  return computeSurfaceEnergy(surfacePositions);
}

void MappedSurfacePotentialEnergy::gradient(
  EigenSupport::ConstRefVecXd simulationDisplacements,
  EigenSupport::RefVecXd simulationGradient) const
{
  const VXd surfacePositions = computeSurfacePositionsFromSimulationDisplacements(simulationDisplacements);

  VXd surfaceGradient = VXd::Zero(surfaceRestPositions_.size());
  computeSurfaceGradient(surfacePositions, surfaceGradient);

  {
    simulationGradient = surfaceFromSimulationDispMap_.transpose() * surfaceGradient;
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
  const VXd surfacePositions = computeSurfacePositionsFromSimulationDisplacements(simulationDisplacements);

  SpMatD surfaceHessian(surfaceRestPositions_.size(), surfaceRestPositions_.size());
  computeSurfaceHessian(surfacePositions, surfaceHessian);

  {
    simulationHessian = surfaceFromSimulationDispMap_.transpose() * surfaceHessian * surfaceFromSimulationDispMap_;
  }
}

double MappedSurfacePotentialEnergy::computeMaxStepSize(
  EigenSupport::ConstRefVecXd simulationDisplacements,
  EigenSupport::ConstRefVecXd trialSimulationDisplacements) const
{
  const VXd surfacePositions = computeSurfacePositionsFromSimulationDisplacements(simulationDisplacements);
  const VXd trialSurfaceDisplacements = computeSurfaceDisplacementsFromSimulationDisplacements(trialSimulationDisplacements);
  return computeSurfaceMaxStepSize(surfacePositions, trialSurfaceDisplacements);
}

double MappedSurfacePotentialEnergy::computeSurfaceMaxStepSize(
  EigenSupport::ConstRefVecXd,
  EigenSupport::ConstRefVecXd) const
{
  return 1.0;
}

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
