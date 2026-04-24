/*
copyright to Bohan Wang
*/

#pragma once

#include "EigenDef.h"
#include "ipc/core/surfaceIPCPairs.h"
#include "../topology/surfaceIPCTopology.h"

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
};

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
