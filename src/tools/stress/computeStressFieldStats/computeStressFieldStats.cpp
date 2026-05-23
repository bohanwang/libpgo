#include "pgoLogging.h"

#include <argparse/argparse.hpp>
#include <fmt/format.h>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

struct Options
{
  std::string stressDir;
  std::string outputPath;
  std::string prefix = "von_mises";
  int frameStart = 0;
  int frameEnd = -1;
};

Options parseOptions(int argc, char *argv[])
{
  argparse::ArgumentParser program("computeStressFieldStats");
  program.add_description(
    "Aggregate per-frame statistics (mean, median, p99, min, max) over a sequence "
    "of stress JSON files like {prefix}{frame:04d}.json (the format written by runIPCSim).");

  program.add_argument("--stress-dir")
    .required()
    .help("Folder containing the per-frame stress JSON files.");
  program.add_argument("--output")
    .required()
    .help("Output JSON file with per-frame stats.");
  program.add_argument("--prefix")
    .default_value(std::string("von_mises"))
    .help("Per-frame filename prefix; files are read as {prefix}{frame:04d}.json.");
  program.add_argument("--frame-start")
    .scan<'i', int>()
    .default_value(0)
    .help("First frame index (inclusive).");
  program.add_argument("--frame-end")
    .scan<'i', int>()
    .default_value(-1)
    .help("Last frame index (exclusive). -1 = auto-detect from files on disk.");

  program.parse_args(argc, argv);

  Options options;
  options.stressDir = program.get<std::string>("--stress-dir");
  options.outputPath = program.get<std::string>("--output");
  options.prefix = program.get<std::string>("--prefix");
  options.frameStart = program.get<int>("--frame-start");
  options.frameEnd = program.get<int>("--frame-end");
  return options;
}

std::filesystem::path framePath(const std::filesystem::path &dir, const std::string &prefix, int frame)
{
  return dir / fmt::format("{}{:04d}.json", prefix, frame);
}

struct FrameStats
{
  int frame = 0;
  double time = 0.0;
  std::size_t count = 0;
  double min = 0.0;
  double mean = 0.0;
  double stddev = 0.0;
  double median = 0.0;
  double p99 = 0.0;
  double max = 0.0;
};

// Linear-interpolated quantile on an already-sorted vector. Matches numpy's
// default 'linear' method: idx = q * (n-1), then lerp between neighbors.
double quantileSorted(const std::vector<double> &sorted, double q)
{
  const std::size_t n = sorted.size();
  if (n == 0)
    return 0.0;
  if (n == 1)
    return sorted[0];

  const double pos = q * static_cast<double>(n - 1);
  const std::size_t lo = static_cast<std::size_t>(std::floor(pos));
  const std::size_t hi = static_cast<std::size_t>(std::ceil(pos));
  if (lo == hi)
    return sorted[lo];

  const double frac = pos - static_cast<double>(lo);
  return sorted[lo] * (1.0 - frac) + sorted[hi] * frac;
}

FrameStats computeFrameStats(const nlohmann::json &doc, const std::filesystem::path &source)
{
  if (!doc.contains("values") || !doc.at("values").is_array())
    throw std::runtime_error("Missing or non-array 'values' field in " + source.string());

  std::vector<double> values = doc.at("values").get<std::vector<double>>();
  if (values.empty())
    throw std::runtime_error("Empty 'values' array in " + source.string());

  std::sort(values.begin(), values.end());

  const double n = static_cast<double>(values.size());

  double sum = 0.0;
  for (double v : values)
    sum += v;
  const double mean = sum / n;

  double sqSum = 0.0;
  for (double v : values) {
    const double d = v - mean;
    sqSum += d * d;
  }

  FrameStats stats;
  stats.frame = doc.value("frame", 0);
  stats.time = doc.value("time", 0.0);
  stats.count = values.size();
  stats.min = values.front();
  stats.max = values.back();
  stats.mean = mean;
  // Population standard deviation: the values are a complete field, not a sample.
  stats.stddev = std::sqrt(sqSum / n);
  stats.median = quantileSorted(values, 0.5);
  stats.p99 = quantileSorted(values, 0.99);
  return stats;
}

int detectFrameEnd(const std::filesystem::path &stressDir, const std::string &prefix, int frameStart)
{
  int frame = frameStart;
  while (std::filesystem::exists(framePath(stressDir, prefix, frame)))
    ++frame;
  return frame;
}

}  // namespace

int main(int argc, char *argv[])
{
  pgo::Logging::init();

  Options options;
  try {
    options = parseOptions(argc, argv);
  }
  catch (const std::exception &err) {
    std::cerr << "Error parsing arguments: " << err.what() << '\n';
    return 1;
  }

  try {
    const std::filesystem::path stressDir(options.stressDir);
    if (!std::filesystem::is_directory(stressDir))
      throw std::runtime_error("stress directory does not exist: " + stressDir.string());

    int frameEnd = options.frameEnd;
    if (frameEnd < 0) {
      frameEnd = detectFrameEnd(stressDir, options.prefix, options.frameStart);
      std::cout << "Auto-detected frame range: [" << options.frameStart << ", " << frameEnd << ")\n";
    }
    if (frameEnd <= options.frameStart)
      throw std::runtime_error(fmt::format(
        "No frames found in {} starting at {}", stressDir.string(), options.frameStart));

    std::string stressType;
    std::string location;
    std::vector<FrameStats> frameStats;
    frameStats.reserve(static_cast<std::size_t>(frameEnd - options.frameStart));

    for (int frame = options.frameStart; frame < frameEnd; ++frame) {
      const std::filesystem::path path = framePath(stressDir, options.prefix, frame);
      std::ifstream in(path);
      if (!in)
        throw std::runtime_error("Failed to open frame file: " + path.string());

      nlohmann::json doc;
      in >> doc;

      if (frame == options.frameStart) {
        stressType = doc.value("stress_type", std::string());
        location = doc.value("location", std::string());
      }

      frameStats.push_back(computeFrameStats(doc, path));
    }

    nlohmann::json out;
    out["stress_type"] = stressType;
    out["location"] = location;
    out["source_dir"] = stressDir.lexically_normal().string();
    out["prefix"] = options.prefix;
    out["frame_start"] = options.frameStart;
    out["frame_end"] = frameEnd;
    out["num_frames"] = static_cast<int>(frameStats.size());

    nlohmann::json framesArr = nlohmann::json::array();
    for (const FrameStats &s : frameStats) {
      framesArr.push_back({
        { "frame", s.frame },
        { "time", s.time },
        { "count", s.count },
        { "min", s.min },
        { "mean", s.mean },
        { "stddev", s.stddev },
        { "median", s.median },
        { "p99", s.p99 },
        { "max", s.max },
      });
    }
    out["frames"] = std::move(framesArr);

    const std::filesystem::path outputPath(options.outputPath);
    if (outputPath.has_parent_path()) {
      std::error_code ec;
      std::filesystem::create_directories(outputPath.parent_path(), ec);
      if (ec)
        throw std::runtime_error("Failed to create output directory: " + outputPath.parent_path().string());
    }

    std::ofstream outFile(outputPath);
    if (!outFile)
      throw std::runtime_error("Failed to open output file: " + outputPath.string());
    outFile << out.dump(2) << '\n';

    std::cout << "Wrote stats for " << frameStats.size() << " frames to "
              << outputPath << '\n';
    return 0;
  }
  catch (const std::exception &err) {
    std::cerr << "Error: " << err.what() << '\n';
    return 1;
  }
}
