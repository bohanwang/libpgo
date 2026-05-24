#pragma once

#include "EigenDef.h"
#include "ipc/core/surfaceIPCPairs.h"
#include "ipc/external/obstacleSurface.h"
#include "ipc/topology/surfaceIPCTopology.h"

#include <vector>

namespace pgo
{
namespace Contact
{
namespace IPC
{

void buildSelfPairs(
  const SurfaceIPCTopology &topology,
  EigenSupport::ConstRefVecXd positions,
  double dhat,
  SelfPairSet &pairs);

void buildSelfPairsLineSearchSuperset(
  const SurfaceIPCTopology &topology,
  EigenSupport::ConstRefVecXd positions,
  EigenSupport::ConstRefVecXd displacements,
  double dhat,
  SelfPairSet &pairs);

void buildExternalPairs(
  const SurfaceIPCTopology &topology,
  EigenSupport::ConstRefVecXd positions,
  const std::vector<ObstacleSurface> &obstacles,
  double dhatExternal,
  ExternalPairSet &pairs);

void buildExternalPairsLineSearchSuperset(
  const SurfaceIPCTopology &topology,
  EigenSupport::ConstRefVecXd positions,
  EigenSupport::ConstRefVecXd displacements,
  const std::vector<ObstacleSurface> &obstacles,
  double dhatExternal,
  ExternalPairSet &pairs);

}  // namespace IPC
}  // namespace Contact
}  // namespace pgo
