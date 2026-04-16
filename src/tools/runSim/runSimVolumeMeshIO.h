#pragma once

#include <memory>
#include <string>

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

VolumeMeshInputConfig parseVolumeMeshInputConfig(const ConfigFileJSON &jconfig);
std::unique_ptr<VolumetricMeshes::VolumetricMesh> loadValidatedVolumeMesh(const VolumeMeshInputConfig &config, double scale);
}  // namespace RunSim
}  // namespace pgo
