/*
copyright to Bohan Wang
*/

#pragma once

#include "EigenDef.h"
#include "ipc/core/surfaceIPCPairs.h"
#include "../topology/surfaceIPCTopology.h"

#include <array>
#include <vector>

namespace pgo
{
namespace Contact
{
namespace CIPC
{

class SurfaceIPCSelfBroadPhase
{
public:
  void buildPairs(
    const SurfaceIPCTopology &topology,
    EigenSupport::ConstRefVecXd positions,
    double dhat,
    std::vector<PTPair> &ptPairs,
    std::vector<EEPair> &eePairs) const;

  // Exact edge-triangle intersection test over primitives that share no
  // vertex and involve at least one deformable vertex. The distance-based
  // checks cannot see two triangles that pass through each other: every
  // point-triangle and edge-edge distance of such a pair is strictly positive,
  // so the barrier is zero there and CCD then keeps the surfaces from
  // separating. Each result is {edge index, triangle index}, sorted ascending.
  void findIntersectingEdgeTrianglePairs(
    const SurfaceIPCTopology &topology,
    EigenSupport::ConstRefVecXd positions,
    std::vector<std::array<int, 2>> &intersections) const;
};

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
