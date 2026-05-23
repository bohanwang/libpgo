#include "cli/cli.h"

#include <argparse/argparse.hpp>

#include <iostream>
#include <string>

namespace pgo::RunIPCSim
{
namespace
{
void configureRunIPCSimArgumentParser(argparse::ArgumentParser &program)
{
  program.add_argument("config")
    .help("Config File")
    .required();
  program.add_argument("--log")
    .help("Write command-line output to a .log file next to the config file")
    .default_value(false)
    .implicit_value(true);
  program.add_argument("--legacy")
    .help("Run volume legacy penalty-contact configs through the runIPCSim entry")
    .default_value(false)
    .implicit_value(true);
}

RunIPCSimCliOptions readRunIPCSimCliOptions(const argparse::ArgumentParser &program)
{
  RunIPCSimCliOptions options;
  options.configPath = program.get<std::string>("config");
  options.runOptions.enableCliLog = program.get<bool>("--log");
  options.runOptions.contactBackendKind = program.get<bool>("--legacy")
    ? ContactBackendKind::LegacyPenalty
    : ContactBackendKind::Ipc;
  return options;
}
}  // namespace

RunIPCSimCliOptions parseRunIPCSimCli(int argc, char *argv[])
{
  argparse::ArgumentParser program("Run IPC Simulation");
  configureRunIPCSimArgumentParser(program);
  program.parse_args(argc, argv);
  return readRunIPCSimCliOptions(program);
}

int runCli(int argc, char *argv[])
{
  argparse::ArgumentParser program("Run IPC Simulation");
  configureRunIPCSimArgumentParser(program);
  try {
    program.parse_args(argc, argv);
  }
  catch (const std::exception &err) {
    std::cerr << err.what() << std::endl;
    std::cerr << program;
    return 1;
  }

  const RunIPCSimCliOptions options = readRunIPCSimCliOptions(program);
  return runFromConfig(options.configPath, options.runOptions);
}
}  // namespace pgo::RunIPCSim
