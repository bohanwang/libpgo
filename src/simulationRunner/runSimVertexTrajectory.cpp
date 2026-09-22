#include "runSimVertexTrajectory.h"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>

namespace pgo::RunSim
{
namespace
{
namespace ES = pgo::EigenSupport;

[[noreturn]] void fail(const std::string &label, const std::string &message)
{
  throw std::invalid_argument("Vertex trajectory " + label + ": " + message);
}

double requireFiniteNumber(const nlohmann::json &value, const std::string &label, const std::string &what)
{
  if (!value.is_number())
    fail(label, what + " must be a number.");
  const double number = value.get<double>();
  if (!std::isfinite(number))
    fail(label, what + " must be finite.");
  return number;
}

VertexTrajectoryTrack parseTrack(const nlohmann::json &track, std::size_t trackIndex, const std::string &label)
{
  const std::string trackLabel = "tracks[" + std::to_string(trackIndex) + "]";
  if (!track.is_object())
    fail(label, trackLabel + " must be an object.");

  VertexTrajectoryTrack result;
  if (track.contains("name")) {
    if (!track["name"].is_string())
      fail(label, trackLabel + ".name must be a string.");
    result.name = track["name"].get<std::string>();
  }
  else {
    result.name = trackLabel;
  }

  if (!track.contains("vertices") || !track["vertices"].is_array() || track["vertices"].empty())
    fail(label, trackLabel + " requires a non-empty `vertices` array.");
  std::vector<int> vertices;
  for (const auto &v : track["vertices"]) {
    if (!v.is_number_integer() || v.get<long long>() < 0)
      fail(label, trackLabel + ".vertices must hold non-negative integers.");
    vertices.push_back(v.get<int>());
  }
  if (std::set<int>(vertices.begin(), vertices.end()).size() != vertices.size())
    fail(label, trackLabel + ".vertices contains a repeated vertex.");
  const int numVertices = static_cast<int>(vertices.size());

  if (!track.contains("keys") || !track["keys"].is_array() || track["keys"].empty())
    fail(label, trackLabel + " requires a non-empty `keys` array.");

  std::vector<double> times;
  std::vector<ES::VXd> positions;
  for (std::size_t ki = 0; ki < track["keys"].size(); ++ki) {
    const auto &key = track["keys"][ki];
    const std::string keyLabel = trackLabel + ".keys[" + std::to_string(ki) + "]";
    if (!key.is_object() || !key.contains("t") || !key.contains("p"))
      fail(label, keyLabel + " must be an object with `t` and `p`.");
    const double t = requireFiniteNumber(key["t"], label, keyLabel + ".t");
    if (!times.empty() && t <= times.back())
      fail(label, keyLabel + ".t must be greater than the previous key time.");
    const auto &p = key["p"];
    if (!p.is_array() || static_cast<int>(p.size()) != numVertices)
      fail(label, keyLabel + ".p must list exactly one position per track vertex (" + std::to_string(numVertices) + ").");
    ES::VXd flat(3 * numVertices);
    for (int vi = 0; vi < numVertices; ++vi) {
      const auto &xyz = p[vi];
      if (!xyz.is_array() || xyz.size() != 3)
        fail(label, keyLabel + ".p[" + std::to_string(vi) + "] must be [x, y, z].");
      for (int c = 0; c < 3; ++c)
        flat[3 * vi + c] = requireFiniteNumber(xyz[c], label, keyLabel + ".p[" + std::to_string(vi) + "]");
    }
    times.push_back(t);
    positions.push_back(std::move(flat));
  }

  if (track.contains("release-time")) {
    result.releaseTime = requireFiniteNumber(track["release-time"], label, trackLabel + ".release-time");
  }

  // Store vertices in ascending order so downstream index handling is uniform;
  // permute every key's positions the same way.
  std::vector<int> order(vertices.size());
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](int a, int b) { return vertices[a] < vertices[b]; });
  result.vertices.resize(vertices.size());
  for (std::size_t i = 0; i < order.size(); ++i)
    result.vertices[i] = vertices[order[i]];
  result.times = std::move(times);
  result.keyPositions.reserve(positions.size());
  for (const ES::VXd &flat : positions) {
    ES::VXd permuted(flat.size());
    for (std::size_t i = 0; i < order.size(); ++i)
      permuted.segment<3>(3 * static_cast<Eigen::Index>(i)) = flat.segment<3>(3 * order[i]);
    result.keyPositions.push_back(std::move(permuted));
  }
  return result;
}
}  // namespace

ES::VXd VertexTrajectoryTrack::sampleAt(double time) const
{
  if (times.empty())
    throw std::logic_error("Vertex trajectory track has no keys.");
  if (time <= times.front())
    return keyPositions.front();
  if (time >= times.back())
    return keyPositions.back();
  const auto upper = std::upper_bound(times.begin(), times.end(), time);
  const std::size_t hi = static_cast<std::size_t>(upper - times.begin());
  const std::size_t lo = hi - 1;
  const double w = (time - times[lo]) / (times[hi] - times[lo]);
  return (1.0 - w) * keyPositions[lo] + w * keyPositions[hi];
}

double VertexTrajectoryTrack::peakSpeed() const
{
  double peak = 0.0;
  for (std::size_t k = 1; k < times.size(); ++k) {
    const double dt = times[k] - times[k - 1];
    for (int vi = 0; vi < numVertices(); ++vi) {
      const double step = (keyPositions[k].segment<3>(3 * vi) - keyPositions[k - 1].segment<3>(3 * vi)).norm();
      peak = std::max(peak, step / dt);
    }
  }
  return peak;
}

double VertexTrajectory::lastKeyTime() const
{
  double last = 0.0;
  for (const VertexTrajectoryTrack &track : tracks)
    last = std::max(last, track.times.back());
  return last;
}

VertexTrajectory parseVertexTrajectory(const std::string &jsonText, const std::string &label)
{
  nlohmann::json document;
  try {
    document = nlohmann::json::parse(jsonText);
  }
  catch (const std::exception &err) {
    fail(label, std::string("invalid JSON: ") + err.what());
  }
  if (!document.is_object())
    fail(label, "top level must be an object.");
  if (!document.contains("format") || !document["format"].is_string() ||
    document["format"].get<std::string>() != "pgo-vertex-trajectory")
    fail(label, "`format` must be \"pgo-vertex-trajectory\".");
  if (!document.contains("version") || !document["version"].is_number_integer() || document["version"].get<int>() != 1)
    fail(label, "`version` must be 1.");

  VertexTrajectory trajectory;
  if (document.contains("space")) {
    const auto &space = document["space"];
    if (!space.is_string())
      fail(label, "`space` must be a string.");
    const std::string value = space.get<std::string>();
    if (value == "position")
      trajectory.space = VertexTrajectory::Space::Position;
    else if (value == "displacement")
      trajectory.space = VertexTrajectory::Space::Displacement;
    else
      fail(label, "`space` must be \"position\" or \"displacement\".");
  }
  if (document.contains("interpolation")) {
    const auto &interp = document["interpolation"];
    if (!interp.is_string() || interp.get<std::string>() != "linear")
      fail(label, "`interpolation` must be \"linear\".");
  }

  if (!document.contains("tracks") || !document["tracks"].is_array() || document["tracks"].empty())
    fail(label, "requires a non-empty `tracks` array.");

  std::set<int> usedVertices;
  for (std::size_t ti = 0; ti < document["tracks"].size(); ++ti) {
    VertexTrajectoryTrack track = parseTrack(document["tracks"][ti], ti, label);
    for (int v : track.vertices) {
      if (!usedVertices.insert(v).second)
        fail(label, "vertex " + std::to_string(v) + " appears in more than one track.");
    }
    trajectory.tracks.push_back(std::move(track));
  }
  return trajectory;
}

VertexTrajectory loadVertexTrajectory(const std::filesystem::path &path)
{
  std::ifstream in(path);
  if (!in.is_open())
    throw std::invalid_argument("Vertex trajectory file cannot be opened: " + path.string());
  std::stringstream buffer;
  buffer << in.rdbuf();
  return parseVertexTrajectory(buffer.str(), path.string());
}

}  // namespace pgo::RunSim
