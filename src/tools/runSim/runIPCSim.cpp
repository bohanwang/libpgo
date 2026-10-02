#include "simulationRunner.h"
#include "pgoLogging.h"

#include <argparse/argparse.hpp>

#include <filesystem>
#include <iostream>

int main(int argc, char *argv[])
{
  pgo::Logging::init();

  argparse::ArgumentParser program("Run IPC Simulation");
  program.add_argument("config").help("Config File").required();
  program.add_argument("--log")
    .help("Write command-line output to a .log file next to the config file")
    .default_value(false)
    .implicit_value(true);
  program.add_argument("--max-frames")
    .help("Execute at most this many frames without changing the configured motion schedule (-1: all remaining)")
    .default_value(-1)
    .scan<'i', int>();
  program.add_argument("--timing-csv")
    .help("Write per-frame compute/output wall times and terminal solve status to this CSV")
    .default_value(std::string());

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
    pgo::SimulationRunner::IPCSimulationRunOptions options;
    options.maxFrames = program.get<int>("--max-frames");
    options.timingOutput = program.get<std::string>("--timing-csv");
    return pgo::SimulationRunner::runIPCSimulationFromConfig(
      configFilename, program.get<bool>("--log"), options);
  }
  catch (const std::exception &err) {
    std::cerr << err.what() << std::endl;
    return 1;
  }
}
