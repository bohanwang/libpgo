#pragma once

#include <memory>
#include <string>
#include <vector>

namespace pgo
{
class ConfigFileJSON;

namespace VolumetricMeshes 
{
class VolumetricMesh;
}

namespace RunSim
{
struct VolumeMeshInputConfig
{
  std::string configKey;
  std::string meshFilename;
  int expectedElementType = -1;
};

struct ResolvedRunSimPaths
{
  std::string surfaceMeshFilename;
  std::string outputPath;
  std::vector<std::string> fixedVertexFilenames;
  std::vector<std::string> externalObjectFilenames;
};

VolumeMeshInputConfig parseVolumeMeshInputConfig(const ConfigFileJSON &jconfig);
ResolvedRunSimPaths resolveRunSimPaths(const ConfigFileJSON &jconfig);
std::unique_ptr<VolumetricMeshes::VolumetricMesh> loadValidatedVolumeMesh(const VolumeMeshInputConfig &config, double scale);
}  // namespace RunSim
}  // namespace pgo
