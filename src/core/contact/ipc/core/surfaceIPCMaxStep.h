#pragma once

#include "EigenDef.h"
#include "ipc/external/obstacleSurface.h"
#include "ipc/topology/surfaceIPCTopology.h"

#include <vector>

namespace pgo
{
namespace Contact
{
namespace CIPC
{

double computeSelfMaxStep(
  const SurfaceIPCTopology &topology,
  EigenSupport::ConstRefVecXd x,
  EigenSupport::ConstRefVecXd dx,
  double dhat,
  double slackness);

double computeExternalMaxStep(
  const SurfaceIPCTopology &topology,
  EigenSupport::ConstRefVecXd x,
  EigenSupport::ConstRefVecXd dx,
  const std::vector<ObstacleSurface> &obstacles,
  double dhatExternal,
  double slackness);

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
