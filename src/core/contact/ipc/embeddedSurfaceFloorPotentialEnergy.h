/*
copyright to Bohan Wang
*/

#pragma once

#include "mappedSurfacePotentialEnergy.h"

#include <limits>

namespace pgo
{
namespace Contact
{
namespace CIPC
{

enum class FloorAxis : int
{
  INVALID = -1,
  X = 0,
  Y = 1,
  Z = 2
};

enum class FloorSide : int
{
  LOWER = 1,
  UPPER = -1
};

struct FloorPenaltyParameters
{
  FloorAxis floorAxis = FloorAxis::INVALID;
  FloorSide floorSide = FloorSide::LOWER;
  double floorHeight = std::numeric_limits<double>::quiet_NaN();
  double floorKappa = std::numeric_limits<double>::quiet_NaN();
};

class EmbeddedSurfaceFloorPotentialEnergy : public MappedSurfacePotentialEnergy
{
public:
  EmbeddedSurfaceFloorPotentialEnergy(
    const EigenSupport::MXd &surfaceRestVertices,
    const EigenSupport::SpMatD &surfaceFromSimulationDispMap,
    const FloorPenaltyParameters &params);

  void setFloorHeight(double h);
  double floorHeight() const;

protected:
  virtual double computeSurfaceEnergy(EigenSupport::ConstRefVecXd surfacePositions) const override;
  virtual void computeSurfaceGradient(
    EigenSupport::ConstRefVecXd surfacePositions,
    EigenSupport::RefVecXd surfaceGradient) const override;
  virtual void computeSurfaceHessian(
    EigenSupport::ConstRefVecXd surfacePositions,
    EigenSupport::SpMatD &surfaceHessian) const override;

private:
  FloorPenaltyParameters params_;
};

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
