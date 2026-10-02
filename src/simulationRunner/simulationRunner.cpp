#include "simulationRunner.h"

#include "configFileJSON.h"
#include "pgoLogging.h"

#include <stdexcept>
#include <string>

int pgo::SimulationRunner::runSimulationFromConfig(const std::filesystem::path &configFilename)
{
  try {
    ConfigFileJSON config;
    if (!config.open(configFilename.string().c_str())) {
      SPDLOG_LOGGER_ERROR(Logging::lgr(), "Unable to open simulation config: {}", configFilename.string());
      return 1;
    }

    const std::string contactModel =
      config.exist("contact-model") ? config.getString("contact-model") : "sampled";

    if (contactModel == "sampled") {
      if (config.exist("simulation-backend") && config.getString("simulation-backend") != "cpu")
        throw std::invalid_argument("The CUDA simulation backend requires `contact-model = ipc`.");
      return runSampledSimulationFromConfig(configFilename);
    }
    if (contactModel == "ipc")
      return runIPCSimulationFromConfig(configFilename);

    SPDLOG_LOGGER_ERROR(
      Logging::lgr(), "Unsupported contact model `{}`; expected `sampled` or `ipc`.", contactModel);
    return 1;
  }
  catch (const std::exception &err) {
    SPDLOG_LOGGER_ERROR(Logging::lgr(), "{}", err.what());
    return 1;
  }
}
