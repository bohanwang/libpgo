#include "io/volumeMeshIO.h"

#include "configFileJSON.h"
#include "cubicMesh.h"
#include "tetMesh.h"
#include "volumetricMesh.h"

#include <fmt/format.h>

#include <stdexcept>

namespace pgo::RunSim
{
namespace
{
using VolumetricMesh = pgo::VolumetricMeshes::VolumetricMesh;

std::string elementTypeName(VolumetricMesh::elementType type)
{
  switch (type) {
  case VolumetricMesh::TET:
    return "tet";
  case VolumetricMesh::CUBIC:
    return "cubic";
  default:
    return "invalid";
  }
}
}

VolumeMeshInputConfig parseVolumeMeshInputConfig(const ConfigFileJSON &jconfig)
{
  const bool hasTetMesh = jconfig.exist("tet-mesh");
  const bool hasCubicMesh = jconfig.exist("cubic-mesh");

  if (hasTetMesh == hasCubicMesh) {
    throw std::invalid_argument("runSim expects exactly one of \"tet-mesh\" or \"cubic-mesh\".");
  }

  VolumeMeshInputConfig config;
  if (hasTetMesh) {
    config.configKey = "tet-mesh";
    config.meshFilename = jconfig.getResolvedPath("tet-mesh", 1);
    config.expectedElementType = VolumetricMesh::TET;
  }
  else {
    config.configKey = "cubic-mesh";
    config.meshFilename = jconfig.getResolvedPath("cubic-mesh", 1);
    config.expectedElementType = VolumetricMesh::CUBIC;
  }

  return config;
}

ResolvedRunSimPaths resolveRunSimPaths(const ConfigFileJSON &jconfig)
{
  ResolvedRunSimPaths paths;
  paths.surfaceMeshFilename = jconfig.getResolvedPath("surface-mesh", 1);
  paths.outputPath = jconfig.getResolvedPath("output", 1);

  if (jconfig.exist("fixed-vertices")) {
    for (const auto &fv : jconfig.handle()["fixed-vertices"]) {
      paths.fixedVertexFilenames.push_back(jconfig.resolvePath(fv["filename"].get<std::string>()));
    }
  }

  if (jconfig.exist("external-objects")) {
    for (const auto &jko : jconfig.handle()["external-objects"]) {
      paths.externalObjectFilenames.push_back(jconfig.resolvePath(jko["filename"].get<std::string>()));
    }
  }

  return paths;
}

std::unique_ptr<VolumetricMeshes::VolumetricMesh> loadValidatedVolumeMesh(const VolumeMeshInputConfig &config, double scale)
{
  using VolumetricMesh = pgo::VolumetricMeshes::VolumetricMesh;

  const auto actualType = VolumetricMesh::getElementType(config.meshFilename.c_str());
  if (actualType == VolumetricMesh::INVALID) {
    throw std::invalid_argument(fmt::format("Failed to determine volumetric mesh type for file: {}", config.meshFilename));
  }

  if (actualType != config.expectedElementType) {
    throw std::invalid_argument(fmt::format(
      "Config key \"{}\" expects a {} volumetric mesh, but file {} is {}.",
      config.configKey, elementTypeName((VolumetricMesh::elementType)config.expectedElementType), config.meshFilename, elementTypeName(actualType)));
  }

  std::unique_ptr<VolumetricMeshes::VolumetricMesh> volumetricMesh;
  switch (actualType) {
  case VolumetricMesh::TET:
    volumetricMesh = std::make_unique<VolumetricMeshes::TetMesh>(config.meshFilename.c_str());
    break;
  case VolumetricMesh::CUBIC:
    volumetricMesh = std::make_unique<VolumetricMeshes::CubicMesh>(config.meshFilename.c_str());
    break;
  default:
    throw std::invalid_argument(fmt::format("Unsupported volumetric mesh type for file: {}", config.meshFilename));
  }

  if (scale != 1.0) {
    for (int vi = 0; vi < volumetricMesh->getNumVertices(); vi++) {
      volumetricMesh->setVertex(vi, volumetricMesh->getVertex(vi) * scale);
    }
  }

  return volumetricMesh;
}
}  // namespace pgo::RunSim
