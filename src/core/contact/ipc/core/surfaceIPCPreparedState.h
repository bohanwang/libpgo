#pragma once

#include "EigenDef.h"
#include "ipc/core/surfaceIPCPairs.h"

namespace pgo
{
namespace Contact
{
namespace CIPC
{

struct SurfaceIPCPreparedState
{
  bool hasState = false;
  EigenSupport::VXd positions;
  SelfPairSet selfPairs;
  ExternalPairSet externalPairs;

  void clear()
  {
    hasState = false;
    positions.resize(0);
    selfPairs.clear();
    externalPairs.clear();
  }

  bool isPreparedFor(EigenSupport::ConstRefVecXd x) const
  {
    return hasState &&
      positions.size() == x.size() &&
      (positions.array() == x.array()).all();
  }
};

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
