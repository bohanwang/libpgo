#pragma once

#include "EigenSupport.h"

#include <array>
#include <filesystem>

namespace pgo
{
class ConfigFileJSON;
}

namespace pgo::RunIPCSim
{
enum class RunIPCSimSimulationMode
{
  Dynamic,
  Static,
};

struct RunIPCSimRuntimeConfig
{
  RunIPCSimSimulationMode simulationMode = RunIPCSimSimulationMode::Dynamic;
  EigenSupport::V3d gravity = EigenSupport::V3d::Zero();
  EigenSupport::V3d initialVelocity = EigenSupport::V3d::Zero();
  double timestep = 0.0;
  double scale = 1.0;
  double solverEps = 0.0;
  int solverMaxIter = 0;
  std::array<double, 2> dampingParams = { 0.0, 0.0 };
  int numSimSteps = 0;
  int frameGap = 1;
  bool restartFromU = false;
  bool dumpDeformEveryFrame = false;
  bool outputVonMises = false;
  bool enableProfiling = false;
  std::filesystem::path outputFolder;
};

RunIPCSimRuntimeConfig parseRunIPCSimRuntimeConfig(const pgo::ConfigFileJSON &config);
}  // namespace pgo::RunIPCSim
