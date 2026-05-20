#include "runIPCSimApp.h"

#include "configFileJSON.h"
#include "initPredicates.h"
#include "pgoLogging.h"
#include "runIPCSimConfig.h"
#include "runIPCSimLogging.h"
#include "runIPCSimLoop.h"
#include "runIPCSimOutput.h"
#include "runIPCSimSession.h"
#include "runIPCSimStaticSolve.h"

#include <tbb/global_control.h>

#include <iostream>
#include <stdexcept>

namespace pgo::RunIPCSim
{
bool openRunIPCSimConfig(const std::filesystem::path &configPath, pgo::ConfigFileJSON &config)
{
  if (config.open(configPath.string().c_str()) != true)
    return false;
  return true;
}

IpcSimulationContext buildIpcSimulation(const pgo::ConfigFileJSON &config)
{
  const bool hasTetMesh = config.exist("tet-mesh");
  const bool hasCubicMesh = config.exist("cubic-mesh");
  const bool useVolumePath = hasTetMesh || hasCubicMesh;
  return useVolumePath ? buildVolumeIpcSimulation(config) : buildShellIpcSimulation(config);
}

IpcSimulationContext buildRunIPCSimSimulation(const pgo::ConfigFileJSON &config, const RunIPCSimOptions &options)
{
  if (options.contactBackendKind == ContactBackendKind::LegacyPenalty)
    return buildVolumeLegacyPenaltySimulation(config);

  return buildIpcSimulation(config);
}

int runFromConfig(const std::filesystem::path &configPath, const RunIPCSimOptions &options)
{
  try {
    tbb::global_control c(tbb::global_control::max_allowed_parallelism, 64);

    pgo::ConfigFileJSON config;
    if (!openRunIPCSimConfig(configPath, config))
      return 0;

    const RunIPCSimRuntimeConfig runtimeConfig = parseRunIPCSimRuntimeConfig(config);
    if (runtimeConfig.outputVonMises && !config.exist("tet-mesh") && !config.exist("cubic-mesh"))
      throw std::invalid_argument("runIPCSim `output-von-mises` requires `tet-mesh` or `cubic-mesh`.");

    RunIPCSimOutput output(runtimeConfig.outputFolder);
    output.prepare(runtimeConfig.restartFromU);
    RunIPCSimRunScope runScope(config, runtimeConfig, options, output);
    if (!runtimeConfig.restartFromU)
      std::cout << "restart-from-u=false; clearing output folder " << runtimeConfig.outputFolder << "." << std::endl;
    pgo::Mesh::initPredicates();

    IpcSimulationContext context = buildRunIPCSimSimulation(config, options);
    if (runtimeConfig.simulationMode == RunIPCSimSimulationMode::Static) {
      if (runtimeConfig.restartFromU)
        throw std::invalid_argument("runIPCSim static mode does not support `restart-from-u`.");
      runIPCSimStaticSolve(runtimeConfig, context, output);
    }
    else {
      RunIPCSimSession session = createRunIPCSimSession(runtimeConfig, context);
      restoreRestartStateIfRequested(runtimeConfig, output, session);
      context.contactBackend->initializeAfterRestart(runtimeConfig, context, session);
      runIPCSimLoop(runtimeConfig, context, session, output);
    }
    runScope.logProfileSummaryIfEnabled();
    return 0;
  }
  catch (const std::exception &err) {
    SPDLOG_LOGGER_ERROR(Logging::lgr(), "{}", err.what());
    return 1;
  }
}
}  // namespace pgo::RunIPCSim
