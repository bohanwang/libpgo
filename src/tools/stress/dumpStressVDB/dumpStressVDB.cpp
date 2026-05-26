#include "pgoLogging.h"
#include "stressFieldVDBExporter.h"

#include <argparse/argparse.hpp>

#include <filesystem>
#include <iostream>
#include <string>

int main(int argc, char *argv[])
{
  pgo::Logging::init();

  argparse::ArgumentParser program("dumpStressVDB");
  program.add_description(
    "Splat per-tet von Mises stress into a per-frame OpenVDB sequence.\n"
    "Inputs are the .veg mesh and the runIPCSim output folder containing\n"
    "  states/deform####.u and stress/von_mises####.json.");
  program.add_argument("--veg")
    .required()
    .help("Path to the .veg volumetric tet mesh.");
  program.add_argument("--sim-output")
    .required()
    .help("Path to the runIPCSim output folder (contains states/ and stress/).");
  program.add_argument("--out")
    .required()
    .help("Output folder for the .vdb sequence.");
  program.add_argument("--frame-start")
    .scan<'i', int>()
    .default_value(0)
    .help("First frame index (inclusive).");
  program.add_argument("--frame-end")
    .scan<'i', int>()
    .default_value(-1)
    .help("Last frame index (exclusive). -1 = auto-detect from files on disk.");
  program.add_argument("--prefix")
    .default_value(std::string("vonMises"))
    .help("Output file prefix; frames written as {prefix}{frame:04d}.vdb.");
  program.add_argument("--voxel-size")
    .scan<'g', double>()
    .default_value(0.0)
    .help("VDB voxel size in world units. <= 0 auto-derives ~1/2 rest edge length.");

  try {
    program.parse_args(argc, argv);
  }
  catch (const std::exception &err) {
    std::cerr << err.what() << '\n'
              << program;
    return 1;
  }

  const std::filesystem::path vegPath = program.get<std::string>("--veg");
  const std::filesystem::path simOutput = program.get<std::string>("--sim-output");
  const std::filesystem::path outDir = program.get<std::string>("--out");
  const std::string prefix = program.get<std::string>("--prefix");
  const double voxelSize = program.get<double>("--voxel-size");
  const int frameStart = program.get<int>("--frame-start");
  int frameEnd = program.get<int>("--frame-end");

  const std::filesystem::path statesDir = simOutput / "states";
  const std::filesystem::path stressDir = simOutput / "stress";

  if (!std::filesystem::is_directory(statesDir)) {
    std::cerr << "states folder not found: " << statesDir << '\n';
    return 1;
  }
  if (!std::filesystem::is_directory(stressDir)) {
    std::cerr << "stress folder not found: " << stressDir << '\n';
    return 1;
  }

  if (frameEnd < 0) {
    int detected = frameStart;
    while (true) {
      const auto u = statesDir / fmt::format("deform{:04d}.u", detected);
      const auto s = stressDir / fmt::format("von_mises{:04d}.json", detected);
      if (!std::filesystem::exists(u) || !std::filesystem::exists(s))
        break;
      ++detected;
    }
    frameEnd = detected;
    std::cout << "Auto-detected frame range: [" << frameStart << ", " << frameEnd << ")\n";
  }

  if (frameEnd <= frameStart) {
    std::cerr << "No frames found in " << statesDir << " starting at " << frameStart << "\n";
    return 1;
  }

  pgo::AnimationIO::StressFieldVDBExporter exporter;
  if (exporter.loadTetMesh(vegPath.string().c_str()) != 0)
    return 1;
  if (exporter.loadDeformationSequence(statesDir.string().c_str(), "deform{:04d}.u",
        frameStart, frameEnd) != 0)
    return 1;
  if (exporter.loadVonMisesSequence(stressDir.string().c_str(), "von_mises{:04d}.json",
        frameStart, frameEnd) != 0)
    return 1;
  if (exporter.exportAnimationVDB(outDir.string().c_str(), prefix.c_str(), voxelSize) != 0)
    return 1;

  std::cout << "Wrote " << exporter.numFrames() << " VDB frames to " << outDir << '\n';
  return 0;
}
