#pragma once

#include "EigenDef.h"
#include "ipc/core/surfaceIPCPairs.h"
#include "ipc/external/obstacleSurface.h"
#include "ipc/topology/surfaceIPCTopology.h"

#include <memory>
#include <vector>

namespace pgo
{
namespace Contact
{
namespace CIPC
{

void buildSelfPairs(
  const SurfaceIPCTopology &topology,
  EigenSupport::ConstRefVecXd positions,
  double dhat,
  SelfPairSet &pairs);

void buildExternalPairs(
  const SurfaceIPCTopology &topology,
  EigenSupport::ConstRefVecXd positions,
  const std::vector<std::shared_ptr<ObstacleSurface>> &obstacles,
  double dhatExternal,
  ExternalPairSet &pairs);

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
