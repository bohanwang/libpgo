#pragma once

#include "embeddedSurfaceFloorPotentialEnergy.h"
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
  Contact::IPC::FloorPenaltyParameters params;
  IpcFloorMotionState motionState;
};

const char *floorAxisToString(Contact::IPC::FloorAxis axis);
const char *floorSideToString(Contact::IPC::FloorSide side);
std::vector<ParsedFloorConfig> parseFloorsConfig(const pgo::ConfigFileJSON &jconfig);
}  // namespace pgo::RunIPCSim
