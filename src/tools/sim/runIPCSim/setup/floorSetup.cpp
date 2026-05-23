#include "setup/floorSetup.h"

#include "configFileJSON.h"
#include "setup/setupCommon.h"

#include <cmath>
#include <string>

namespace pgo::RunIPCSim
{
using pgo::Contact::CIPC::FloorAxis;
using pgo::Contact::CIPC::FloorPenaltyParameters;
using pgo::Contact::CIPC::FloorSide;

FloorAxis parseFloorAxis(const std::string &axis)
{
  if (axis == "x")
    return FloorAxis::X;
  if (axis == "y")
    return FloorAxis::Y;
  if (axis == "z")
    return FloorAxis::Z;
  throwConfigError("floor `axis` must be one of: x, y, z.");
}
const char *floorAxisToString(FloorAxis axis)
{
  switch (axis) {
    case FloorAxis::X:
      return "x";
    case FloorAxis::Y:
      return "y";
    case FloorAxis::Z:
      return "z";
    default:
      return "invalid";
  }
}
FloorSide parseFloorSide(const std::string &side)
{
  if (side == "lower")
    return FloorSide::LOWER;
  if (side == "upper")
    return FloorSide::UPPER;
  throwConfigError("floor `side` must be `lower` or `upper`.");
}
const char *floorSideToString(FloorSide side)
{
  switch (side) {
    case FloorSide::LOWER:
      return "lower";
    case FloorSide::UPPER:
      return "upper";
    default:
      return "invalid";
  }
}
std::vector<ParsedFloorConfig> parseFloorsConfig(const pgo::ConfigFileJSON &jconfig)
{
  for (const char *legacyField : { "use-floor", "floor-axis", "floor-height", "floor-kappa" }) {
    if (jconfig.exist(legacyField))
      throwConfigError(std::string("runIPCSim floor config field `") + legacyField + "` has been replaced by `floors[]`.");
  }

  std::vector<ParsedFloorConfig> floors;
  if (!jconfig.exist("floors"))
    return floors;

  const auto &floorsJson = jconfig.handle()["floors"];
  if (!floorsJson.is_array())
    throwConfigError("`floors` must be a JSON array.");

  floors.reserve(floorsJson.size());
  for (std::size_t floorIndex = 0; floorIndex < floorsJson.size(); ++floorIndex) {
    const auto &floorJson = floorsJson.at(floorIndex);
    if (!floorJson.is_object())
      throwConfigError("Each `floors[]` entry must be a JSON object.");

    if (!floorJson.contains("axis"))
      throwConfigError("Missing required field `floors[].axis`.");
    if (!floorJson.contains("kappa"))
      throwConfigError("Missing required field `floors[].kappa`.");

    const bool hasHeight = floorJson.contains("height");
    const bool hasMotion = floorJson.contains("motion");
    if (hasHeight == hasMotion)
      throwConfigError("Each `floors[]` entry must provide exactly one of `height` or `motion`.");

    ParsedFloorConfig floorConfig;
    floorConfig.params.floorAxis = parseFloorAxis(floorJson.at("axis").get<std::string>());
    floorConfig.params.floorSide = floorJson.contains("side") ? parseFloorSide(floorJson.at("side").get<std::string>()) : FloorSide::LOWER;
    floorConfig.params.floorKappa = floorJson.at("kappa").get<double>();
    if (!std::isfinite(floorConfig.params.floorKappa))
      throwConfigError("`floors[].kappa` must be finite.");

    if (hasHeight) {
      floorConfig.params.floorHeight = floorJson.at("height").get<double>();
      if (!std::isfinite(floorConfig.params.floorHeight))
        throwConfigError("`floors[].height` must be finite.");
      floorConfig.motionState.hasMotion = false;
      floorConfig.motionState.heightStart = floorConfig.params.floorHeight;
      floorConfig.motionState.heightEnd = floorConfig.params.floorHeight;
    }
    else {
      const auto &motionJson = floorJson.at("motion");
      if (!motionJson.is_object())
        throwConfigError("`floors[].motion` must be a JSON object.");
      for (const char *field : { "height-start", "height-end", "frame-start", "frame-end" }) {
        if (!motionJson.contains(field))
          throwConfigError(std::string("Missing required field `floors[].motion.") + field + "`.");
      }
      floorConfig.motionState.hasMotion = true;
      floorConfig.motionState.heightStart = motionJson.at("height-start").get<double>();
      floorConfig.motionState.heightEnd = motionJson.at("height-end").get<double>();
      floorConfig.motionState.frameStart = motionJson.at("frame-start").get<int>();
      floorConfig.motionState.frameEnd = motionJson.at("frame-end").get<int>();
      if (!std::isfinite(floorConfig.motionState.heightStart) || !std::isfinite(floorConfig.motionState.heightEnd))
        throwConfigError("`floors[].motion` heights must be finite.");
      if (floorConfig.motionState.frameEnd < floorConfig.motionState.frameStart)
        throwConfigError("`floors[].motion.frame-end` must be greater than or equal to `frame-start`.");
      floorConfig.params.floorHeight = floorHeightAtFrame(floorConfig.motionState, 0);
    }

    floors.push_back(floorConfig);
  }

  return floors;
}
double floorHeightAtFrame(const IpcFloorMotionState &motion, int frame)
{
  if (!motion.hasMotion)
    return motion.heightStart;

  if (frame <= motion.frameStart)
    return motion.heightStart;
  if (frame >= motion.frameEnd)
    return motion.heightEnd;

  const double denom = static_cast<double>(motion.frameEnd - motion.frameStart);
  const double alpha = denom > 0.0 ? static_cast<double>(frame - motion.frameStart) / denom : 1.0;
  return motion.heightStart * (1.0 - alpha) + motion.heightEnd * alpha;
}
}  // namespace pgo::RunIPCSim
