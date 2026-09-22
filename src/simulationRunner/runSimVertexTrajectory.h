#pragma once

// Prescribed vertex trajectories ("pgo-vertex-trajectory" JSON files) used by
// the `fixed-vertices` entries with a `trajectory` field. A track groups the
// vertices that share one set of keyframes; every key stores the simulation
// time and one position per track vertex. Sampling is linear between keys and
// holds the first/last key outside their range.
#include "EigenSupport.h"

#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace pgo
{
namespace RunSim
{

struct VertexTrajectoryTrack
{
  std::string name;
  std::vector<int> vertices;                     // simulation-mesh vertex indices, ascending
  std::vector<double> times;                     // seconds, strictly increasing
  std::vector<EigenSupport::VXd> keyPositions;   // one 3 * numVertices() vector per key
  std::optional<double> releaseTime;             // the track is free for time >= releaseTime

  int numVertices() const { return static_cast<int>(vertices.size()); }
  int numKeys() const { return static_cast<int>(times.size()); }
  bool isActiveAt(double time) const { return !releaseTime.has_value() || time < *releaseTime; }
  // Clamped piecewise-linear interpolation of the key positions.
  EigenSupport::VXd sampleAt(double time) const;
  // Largest vertex speed implied by consecutive keys, in length units per second.
  double peakSpeed() const;
};

struct VertexTrajectory
{
  enum class Space
  {
    Position,      // absolute coordinates in the mesh file's frame (before `scale`)
    Displacement,  // offsets from the rest position (before `scale`)
  };

  Space space = Space::Position;
  std::vector<VertexTrajectoryTrack> tracks;

  double lastKeyTime() const;
};

// Parses and validates the JSON text of a trajectory file. `label` names the
// source in error messages.
VertexTrajectory parseVertexTrajectory(const std::string &jsonText, const std::string &label);
VertexTrajectory loadVertexTrajectory(const std::filesystem::path &path);

}  // namespace RunSim
}  // namespace pgo
