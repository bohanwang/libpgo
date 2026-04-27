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

double floorSideToSign(FloorSide side)
{
  switch (side) {
    case FloorSide::LOWER:
      return 1.0;
    case FloorSide::UPPER:
      return -1.0;
    default:
      throw std::invalid_argument("FloorPenaltyParameters.floorSide must be LOWER or UPPER.");
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
  (void)floorSideToSign(params_.floorSide);
  if (!std::isfinite(params_.floorHeight))
    throw std::invalid_argument("FloorPenaltyParameters.floorHeight must be finite.");
  if (!std::isfinite(params_.floorKappa))
    throw std::invalid_argument("FloorPenaltyParameters.floorKappa must be finite.");
}

void EmbeddedSurfaceFloorPotentialEnergy::setFloorHeight(double h)
{
  if (!std::isfinite(h))
    throw std::invalid_argument("FloorPenaltyParameters.floorHeight must be finite.");
  params_.floorHeight = h;
}

double EmbeddedSurfaceFloorPotentialEnergy::floorHeight() const
{
  return params_.floorHeight;
}

double EmbeddedSurfaceFloorPotentialEnergy::computeSurfaceEnergy(EigenSupport::ConstRefVecXd surfacePositions) const
{
  const int axis = floorAxisToIndex(params_.floorAxis);
  const double sideSign = floorSideToSign(params_.floorSide);
  double energy = 0.0;
  for (int vi = 0; vi < surfacePositions.size() / 3; ++vi) {
    const double dzEff = sideSign * (surfacePositions[3 * vi + axis] - params_.floorHeight);
    if (dzEff < 0.0)
      energy += 0.5 * params_.floorKappa * dzEff * dzEff;
  }
  return energy;
}

void EmbeddedSurfaceFloorPotentialEnergy::computeSurfaceGradient(
  EigenSupport::ConstRefVecXd surfacePositions,
  EigenSupport::RefVecXd surfaceGradient) const
{
  const int axis = floorAxisToIndex(params_.floorAxis);
  const double sideSign = floorSideToSign(params_.floorSide);
  surfaceGradient.setZero();
  for (int vi = 0; vi < surfacePositions.size() / 3; ++vi) {
    const double dzEff = sideSign * (surfacePositions[3 * vi + axis] - params_.floorHeight);
    if (dzEff < 0.0)
      surfaceGradient[3 * vi + axis] = params_.floorKappa * dzEff * sideSign;
  }
}

void EmbeddedSurfaceFloorPotentialEnergy::computeSurfaceHessian(
  EigenSupport::ConstRefVecXd surfacePositions,
  EigenSupport::SpMatD &surfaceHessian) const
{
  const int axis = floorAxisToIndex(params_.floorAxis);
  const double sideSign = floorSideToSign(params_.floorSide);
  std::vector<EigenSupport::TripletD> triplets;
  triplets.reserve(static_cast<std::size_t>(surfacePositions.size() / 3));
  for (int vi = 0; vi < surfacePositions.size() / 3; ++vi) {
    const double dzEff = sideSign * (surfacePositions[3 * vi + axis] - params_.floorHeight);
    if (dzEff < 0.0) {
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
