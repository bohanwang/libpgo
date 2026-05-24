#pragma once

#include "EigenDef.h"
#include "ipc/external/obstacleSurface.h"
#include "ipc/topology/surfaceIPCTopology.h"

#include <vector>

namespace pgo
{
namespace Contact
{
namespace IPC
{

// `thickness` is the CCD minimum-separation distance (xi). The returned alpha
// is the largest step for which no entity pair comes within `thickness` of
// another. With thickness = 0, this is classic contact-at-zero-distance CCD.
// Broad-phase swept AABBs are inflated by `thickness` so the prune stays
// sound under min-separation contact semantics.
double computeSelfMaxStep(
  const SurfaceIPCTopology &topology,
  EigenSupport::ConstRefVecXd x,
  EigenSupport::ConstRefVecXd dx,
  double dhat,
  double slackness,
  double thickness = 0.0);

double computeExternalMaxStep(
  const SurfaceIPCTopology &topology,
  EigenSupport::ConstRefVecXd x,
  EigenSupport::ConstRefVecXd dx,
  const std::vector<ObstacleSurface> &obstacles,
  double dhatExternal,
  double slackness,
  double thickness = 0.0);

}  // namespace IPC
}  // namespace Contact
}  // namespace pgo
