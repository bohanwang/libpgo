#pragma once

#include "EigenDef.h"

#include <cstdint>
#include <vector>

namespace pgo
{
namespace Contact
{
namespace CIPC
{

struct ObstacleSurfaceView
{
  int32_t objectId = -1;
  const EigenSupport::VXd *previousPositions = nullptr;
  const EigenSupport::VXd *currentPositions = nullptr;
  const EigenSupport::MXi *triangles = nullptr;
  const EigenSupport::MXi *uniqueEdges = nullptr;
  const std::vector<double> *triAreas = nullptr;
  const std::vector<double> *edgeLengths = nullptr;
};

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
