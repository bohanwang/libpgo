#include "argparse/argparse.hpp"
#include "configFileJSON.h"
#include "fileService.h"
#include "initPredicates.h"
#include "pgoLogging.h"
#include "runSimCore.h"

#include <cstdio>
#include <exception>
#include <filesystem>
#include <iostream>
#include <string>

namespace {

std::filesystem::path configLogPath(const std::string& configPath) {
    std::filesystem::path logPath(configPath);
    logPath.replace_extension(".log");
    return logPath;
}

bool redirectStdStreamsToLog(const std::filesystem::path& logPath) {
    if (logPath.has_parent_path() && !logPath.parent_path().empty()) {
        std::error_code ec;
        std::filesystem::create_directories(logPath.parent_path(), ec);
        if (ec) {
            std::cerr << "Cannot create log directory: " << logPath.parent_path().string() << std::endl;
            return false;
        }
    }

    std::fflush(stdout);
    std::fflush(stderr);

    const std::string logPathString = logPath.string();
    if (std::freopen(logPathString.c_str(), "w", stdout) == nullptr) {
        std::cerr << "Cannot redirect stdout to log file: " << logPathString << std::endl;
        return false;
    }
    if (std::freopen(logPathString.c_str(), "a", stderr) == nullptr) {
        std::cerr << "Cannot redirect stderr to log file: " << logPathString << std::endl;
        return false;
    }

    return true;
}

}  // namespace

int main(int argc, char* argv[]) {
    argparse::ArgumentParser program("Run Simulation");
    program.add_argument("config").help("Config File").default_value(std::string("examples/box/box-ipc.json"));
    program.add_argument("--deterministic")
        .help("Force deterministic mode by running the simulation single-threaded")
        .default_value(false)
        .implicit_value(true);
    program.add_argument("--save-log")
        .help("Save run stdout and stderr to a .log file next to the config file")
        .default_value(false)
        .implicit_value(true);

    try {
        program.parse_args(argc, argv);
    } catch (const std::exception& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        return 1;
    }

    if (program.get<bool>("--save-log")) {
        const std::filesystem::path logPath = configLogPath(program.get<std::string>("config"));
        if (!redirectStdStreamsToLog(logPath)) {
            return 1;
        }
    }

    pgo::Logging::init();
    pgo::Mesh::initPredicates();

    const pgo::FileService::ConfigPathResolver configPathResolver(program.get<std::string>("config"));

    pgo::ConfigFileJSON jconfig;
    if (jconfig.open(configPathResolver.filePath().c_str()) != true) {
        return 0;
    }

    pgo::api::RunSimConfig config;
    try {
        config = pgo::api::parseRunSimConfig(jconfig, configPathResolver.filePath());
    } catch (const std::exception& err) {
        SPDLOG_LOGGER_ERROR(pgo::Logging::lgr(), "{}", err.what());
        return 1;
    }

    if (program.get<bool>("--deterministic")) {
        config.runtime.deterministicMode = true;
    }

    // config.runtime.restartEnabled = true;
    // config.runtime.restartPause = true;

    try {
        return pgo::api::runSimFromConfig(config);
    } catch (const std::exception& err) {
        SPDLOG_LOGGER_CRITICAL(pgo::Logging::lgr(), "Simulation aborted: {}", err.what());
        return 1;
    }
}
