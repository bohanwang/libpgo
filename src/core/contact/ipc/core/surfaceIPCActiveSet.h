#pragma once

#include "EigenDef.h"
#include "ipc/core/surfaceIPCPairs.h"

#include <cstddef>

namespace pgo
{
namespace Contact
{
namespace CIPC
{

struct SurfaceIPCActiveSet
{
  EigenSupport::VXd positions;
  SelfPairSet selfPairs;
  ExternalPairSet externalPairs;

  void clear()
  {
    positions.resize(0);
    selfPairs.clear();
    externalPairs.clear();
  }

  std::size_t size() const
  {
    return selfPairs.size() + externalPairs.size();
  }
};

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
