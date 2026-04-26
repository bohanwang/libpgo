#include "animationLoader.h"
#include "configFileJSON.h"
#include "pgoLogging.h"
#include "initPredicates.h"

#include <filesystem>
#include <iostream>

int main(int argc, char *argv[])
{
  pgo::Mesh::initPredicates();
  pgo::Logging::init();

  if (argc != 2 && argc != 3) {
    std::cerr << "Usage: convertAnimation <anim.json> [output-folder]\n";
    return 1;
  }

  pgo::ConfigFileJSON config;
  if (!config.open(argv[1])) {
    return 1;
  }

  const std::filesystem::path configPath(argv[1]);
  const std::filesystem::path outputFolder =
    (argc == 3) ? std::filesystem::path(argv[2]) :
    (config.exist("output-folder") ? std::filesystem::path(config.getResolvedPath("output-folder", 1)) :
                                     (configPath.has_parent_path() ? configPath.parent_path() : std::filesystem::path(".")));

  std::error_code ec;
  std::filesystem::create_directories(outputFolder, ec);
  if (ec) {
    std::cerr << "Failed to create output folder " << outputFolder << ": " << ec.message() << "\n";
    return 1;
  }

  pgo::AnimationIO::AnimationLoader loader;
  if (loader.load(argv[1]) != 0) {
    return 1;
  }

  return loader.saveABC(outputFolder.string().c_str());
}
