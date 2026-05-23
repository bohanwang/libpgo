#pragma once

#include "ipc/embeddedSurfaceFloorPotentialEnergy.h"
#include "setup/setup.h"

#include <vector>

namespace pgo
{
class ConfigFileJSON;
}

namespace pgo::RunIPCSim
{
struct ParsedFloorConfig
{
  Contact::CIPC::FloorPenaltyParameters params;
  IpcFloorMotionState motionState;
};

const char *floorAxisToString(Contact::CIPC::FloorAxis axis);
const char *floorSideToString(Contact::CIPC::FloorSide side);
std::vector<ParsedFloorConfig> parseFloorsConfig(const pgo::ConfigFileJSON &jconfig);
}  // namespace pgo::RunIPCSim
