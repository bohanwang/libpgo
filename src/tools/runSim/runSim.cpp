#include "simulationRunner.h"

#include <argparse/argparse.hpp>

#include <filesystem>
#include <iostream>

int main(int argc, char *argv[])
{
  argparse::ArgumentParser program("Run Simulation");
  program.add_argument("config").help("Config File").required();
  program.add_argument("--log")
    .help("Write command-line output to a .log file next to the config file")
    .default_value(false)
    .implicit_value(true);

  try {
    program.parse_args(argc, argv);
  }
  catch (const std::exception &err) {
    std::cerr << err.what() << std::endl;
    std::cerr << program;
    return 1;
  }

  const std::filesystem::path configFilename = program.get<std::string>("config");
  try {
    return pgo::SimulationRunner::runSampledSimulationFromConfig(
      configFilename, program.get<bool>("--log"));
  }
  catch (const std::exception &err) {
    std::cerr << err.what() << std::endl;
    return 1;
  }
}
