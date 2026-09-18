#pragma once

#include "EigenSupport.h"
#include "ipc/assembly/ipcCollisionSurfaceAssembler.h"

#include <filesystem>
#include <vector>

namespace pgo
{
class ConfigFileJSON;

namespace RunIPCSim
{

struct IPCExternalObjectSpec
{
  std::filesystem::path path;
  double scale = 1.0;
  EigenSupport::V3d initialTranslation = EigenSupport::V3d::Zero();
  EigenSupport::V3d movement = EigenSupport::V3d::Zero();
};

std::vector<IPCExternalObjectSpec> parseIPCExternalObjectSpecs(
  const ConfigFileJSON &config);

std::vector<Contact::CIPC::IPCExternalSurface> loadIPCExternalSurfaces(
  const std::vector<IPCExternalObjectSpec> &specs);

Contact::CIPC::IPCCollisionSurfaceData prepareIPCCollisionSurface(
  const EigenSupport::MXd &deformableVertices,
  const EigenSupport::MXi &deformableTriangles,
  const EigenSupport::SpMatD &deformableDisplacementMap,
  const ConfigFileJSON &config);

}  // namespace RunIPCSim
}  // namespace pgo
