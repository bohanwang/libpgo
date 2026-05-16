#pragma once

#include "EigenDef.h"
#include "ipc/core/surfaceIPCPairs.h"

#include <vector>

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
  std::vector<PTPair> ptPairs;
  std::vector<EEPair> eePairs;
  std::vector<ExternalPTPair> externalPTPairs;
  std::vector<ExternalTPPair> externalTPPairs;
  std::vector<ExternalEEPair> externalEEPairs;

  void clear()
  {
    hasState = false;
    positions.resize(0);
    ptPairs.clear();
    eePairs.clear();
    externalPTPairs.clear();
    externalTPPairs.clear();
    externalEEPairs.clear();
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
