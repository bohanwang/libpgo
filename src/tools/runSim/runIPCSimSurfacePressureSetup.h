#pragma once

#include "EigenSupport.h"
#include "triMeshGeo.h"

namespace pgo
{
class ConfigFileJSON;
}

namespace pgo::RunIPCSim
{
struct ParsedSurfacePressureForceConfig
{
  bool enabled = false;
  double pressure = 0.0;
  int rampSteps = 1;
};

ParsedSurfacePressureForceConfig parseSurfacePressureForceConfig(const pgo::ConfigFileJSON &jconfig);
EigenSupport::VXd computeSurfacePressureSimulationForce(const pgo::Mesh::TriMeshGeo &surfaceMesh,
  const EigenSupport::SpMatD &surfaceFromSimulationDispMap,
  const ParsedSurfacePressureForceConfig &pressureConfig, int simulationDofCount);
}  // namespace pgo::RunIPCSim
