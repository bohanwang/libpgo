#include "runIPCSimConfig.h"

#include "configFileJSON.h"

#include <stdexcept>
#include <string>

namespace pgo::RunIPCSim
{
namespace ES = pgo::EigenSupport;

RunIPCSimRuntimeConfig parseRunIPCSimRuntimeConfig(const pgo::ConfigFileJSON &config)
{
  RunIPCSimRuntimeConfig runtimeConfig;
  runtimeConfig.gravity = ES::Mp<ES::V3d>(config.getValue<std::array<double, 3>>("g", 1).data());
  runtimeConfig.initialVelocity = ES::Mp<ES::V3d>(config.getValue<std::array<double, 3>>("init-vel", 1).data());
  runtimeConfig.timestep = config.getDouble("timestep", 1);
  runtimeConfig.scale = config.getDouble("scale", 1);
  runtimeConfig.solverEps = config.getDouble("solver-eps", 1);
  runtimeConfig.solverMaxIter = config.getInt("solver-max-iter", 1);
  runtimeConfig.dampingParams = config.getValue<std::array<double, 2>>("damping-params", 1);
  runtimeConfig.numSimSteps = config.getInt("num-timestep", 1);
  runtimeConfig.frameGap = config.getInt("dump-interval", 1);

  const std::string simType = config.getString("sim-type");
  if (simType == "dynamic") {
    runtimeConfig.simulationMode = RunIPCSimSimulationMode::Dynamic;
  }
  else if (simType == "static") {
    runtimeConfig.simulationMode = RunIPCSimSimulationMode::Static;
  }
  else {
    throw std::invalid_argument("runIPCSim only supports `sim-type = dynamic` or `sim-type = static`.");
  }

  runtimeConfig.outputFolder = config.getResolvedPath("output", 1);
  runtimeConfig.restartFromU = config.exist("restart-from-u") ? config.getValue<bool>("restart-from-u", 1) : false;
  runtimeConfig.dumpDeformEveryFrame = config.exist("dump_deform_every_frame")
    ? config.getValue<bool>("dump_deform_every_frame", 1)
    : false;
  runtimeConfig.outputVonMises = config.exist("output-von-mises") ? config.getValue<bool>("output-von-mises", 1) : false;
  runtimeConfig.enableProfiling = config.exist("profiling") ? config.getValue<bool>("profiling", 1) : false;
  return runtimeConfig;
}
}  // namespace pgo::RunIPCSim
