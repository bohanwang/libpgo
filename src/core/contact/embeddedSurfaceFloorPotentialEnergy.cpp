/*
copyright to Bohan Wang
*/

#include "embeddedSurfaceFloorPotentialEnergy.h"

#include <cmath>
#include <stdexcept>
#include <vector>

namespace pgo
{
namespace Contact
{
namespace CIPC
{

namespace
{
int floorAxisToIndex(FloorAxis axis)
{
  switch (axis) {
    case FloorAxis::X:
      return 0;
    case FloorAxis::Y:
      return 1;
    case FloorAxis::Z:
      return 2;
    default:
      throw std::invalid_argument("FloorPenaltyParameters.floorAxis must be X, Y, or Z.");
  }
}
}  // namespace

EmbeddedSurfaceFloorPotentialEnergy::EmbeddedSurfaceFloorPotentialEnergy(
  const EigenSupport::MXd &surfaceRestVertices,
  const EigenSupport::SpMatD &surfaceFromSimulationDispMap,
  const FloorPenaltyParameters &params):
  MappedSurfacePotentialEnergy(surfaceRestVertices, surfaceFromSimulationDispMap),
  params_(params)
{
  (void)floorAxisToIndex(params_.floorAxis);
  if (!std::isfinite(params_.floorHeight))
    throw std::invalid_argument("FloorPenaltyParameters.floorHeight must be finite.");
  if (!std::isfinite(params_.floorKappa))
    throw std::invalid_argument("FloorPenaltyParameters.floorKappa must be finite.");
}

double EmbeddedSurfaceFloorPotentialEnergy::computeSurfaceEnergy(EigenSupport::ConstRefVecXd surfacePositions) const
{
  const int axis = floorAxisToIndex(params_.floorAxis);
  double energy = 0.0;
  for (int vi = 0; vi < surfacePositions.size() / 3; ++vi) {
    const double dz = surfacePositions[3 * vi + axis] - params_.floorHeight;
    if (dz < 0.0)
      energy += 0.5 * params_.floorKappa * dz * dz;
  }
  return energy;
}

void EmbeddedSurfaceFloorPotentialEnergy::computeSurfaceGradient(
  EigenSupport::ConstRefVecXd surfacePositions,
  EigenSupport::RefVecXd surfaceGradient) const
{
  const int axis = floorAxisToIndex(params_.floorAxis);
  surfaceGradient.setZero();
  for (int vi = 0; vi < surfacePositions.size() / 3; ++vi) {
    const double dz = surfacePositions[3 * vi + axis] - params_.floorHeight;
    if (dz < 0.0)
      surfaceGradient[3 * vi + axis] = params_.floorKappa * dz;
  }
}

void EmbeddedSurfaceFloorPotentialEnergy::computeSurfaceHessian(
  EigenSupport::ConstRefVecXd surfacePositions,
  EigenSupport::SpMatD &surfaceHessian) const
{
  const int axis = floorAxisToIndex(params_.floorAxis);
  std::vector<EigenSupport::TripletD> triplets;
  triplets.reserve(static_cast<std::size_t>(surfacePositions.size() / 3));
  for (int vi = 0; vi < surfacePositions.size() / 3; ++vi) {
    const double dz = surfacePositions[3 * vi + axis] - params_.floorHeight;
    if (dz < 0.0) {
      const int row = 3 * vi + axis;
      triplets.emplace_back(row, row, params_.floorKappa);
    }
  }

  surfaceHessian.resize(surfacePositions.size(), surfacePositions.size());
  surfaceHessian.setFromTriplets(triplets.begin(), triplets.end());
}

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
