#include "runIPCSimSurfacePressureSetup.h"

#include "configFileJSON.h"
#include "runIPCSimSetupCommon.h"
#include "triMeshPseudoNormal.h"

#include <cmath>

namespace pgo::RunIPCSim
{
namespace ES = pgo::EigenSupport;

ParsedSurfacePressureForceConfig parseSurfacePressureForceConfig(const pgo::ConfigFileJSON &jconfig)
{
  ParsedSurfacePressureForceConfig pressureConfig;
  if (!jconfig.exist("surface-pressure-force"))
    return pressureConfig;

  const auto &pressureJson = jconfig.handle()["surface-pressure-force"];
  if (!pressureJson.is_object())
    throwConfigError("`surface-pressure-force` must be a JSON object.");

  pressureConfig.enabled = pressureJson.contains("enabled")
    ? pressureJson.at("enabled").get<bool>()
    : false;
  if (!pressureConfig.enabled)
    return pressureConfig;

  if (!pressureJson.contains("pressure"))
    throwConfigError("Missing required field `surface-pressure-force.pressure` when enabled.");

  pressureConfig.pressure = pressureJson.at("pressure").get<double>();
  if (pressureJson.contains("ramp-steps"))
    pressureConfig.rampSteps = pressureJson.at("ramp-steps").get<int>();

  if (!std::isfinite(pressureConfig.pressure))
    throwConfigError("`surface-pressure-force.pressure` must be finite.");
  if (pressureConfig.rampSteps <= 0)
    throwConfigError("`surface-pressure-force.ramp-steps` must be positive.");

  return pressureConfig;
}
ES::VXd computeSurfacePressureSimulationForce(const pgo::Mesh::TriMeshGeo &surfaceMesh,
  const ES::SpMatD &surfaceFromSimulationDispMap,
  const ParsedSurfacePressureForceConfig &pressureConfig, int simulationDofCount)
{
  if (!pressureConfig.enabled)
    return {};

  std::vector<double> vertexAreas(surfaceMesh.numVertices(), 0.0);
  surfaceMesh.ref().computeVertexSurfaceAreas(vertexAreas.data());

  pgo::Mesh::TriMeshPseudoNormal meshNormal(surfaceMesh);

  ES::VXd surfaceForce = ES::VXd::Zero(surfaceMesh.numVertices() * 3);
  for (int vi = 0; vi < surfaceMesh.numVertices(); ++vi) {
    surfaceForce.segment<3>(vi * 3) = meshNormal.vtxNormal(vi) * pressureConfig.pressure * vertexAreas[vi] * -1.0;
  }

  ES::VXd simulationForce(surfaceFromSimulationDispMap.cols());
  ES::mv(surfaceFromSimulationDispMap, surfaceForce, simulationForce, 1);
  if (simulationForce.size() != simulationDofCount)
    throwConfigError("surface-pressure-force projected force has an unexpected simulation DOF count.");

  return simulationForce;
}
}  // namespace pgo::RunIPCSim
