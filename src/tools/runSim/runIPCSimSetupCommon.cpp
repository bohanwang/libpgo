#include "runIPCSimSetupCommon.h"

#include "configFileJSON.h"

#include <algorithm>
#include <array>
#include <stdexcept>
#include <vector>

namespace pgo::RunIPCSim
{
namespace ES = pgo::EigenSupport;

[[noreturn]] void throwConfigError(const std::string &message)
{
  throw std::invalid_argument(message);
}
void rejectIfPresent(const pgo::ConfigFileJSON &config, const char *field, const char *reason)
{
  if (config.exist(field))
    throwConfigError(std::string("runIPCSim phase2does not support `") + field + "`: " + reason);
}
ES::SpMatD makeIdentityEmbedding(int n3)
{
  std::vector<ES::TripletD> triplets;
  triplets.reserve(n3);
  for (int i = 0; i < n3; ++i)
    triplets.emplace_back(i, i, 1.0);

  ES::SpMatD W(n3, n3);
  W.setFromTriplets(triplets.begin(), triplets.end());
  return W;
}
void validateZeroInitialDisplacement(const pgo::ConfigFileJSON &jconfig)
{
  const ES::V3d initialDisp = ES::Mp<ES::V3d>(jconfig.getValue<std::array<double, 3>>("init-disp", 1).data());
  if (initialDisp.cwiseAbs().maxCoeff() > 0.0)
    throwConfigError("runIPCSim phase2only accepts zero `init-disp`.");
}
void loadSurfaceMeshAndRestPositions(const std::string &surfaceMeshFilename, double scale,
  pgo::Mesh::TriMeshGeo &surfaceMesh, ES::VXd &surfaceRestPositions)
{
  if (!surfaceMesh.load(surfaceMeshFilename))
    throw std::runtime_error("Failed to load surface mesh: " + surfaceMeshFilename);

  for (int vi = 0; vi < surfaceMesh.numVertices(); ++vi)
    surfaceMesh.pos(vi) *= scale;

  surfaceRestPositions.resize(surfaceMesh.numVertices() * 3);
  for (int vi = 0; vi < surfaceMesh.numVertices(); ++vi)
    surfaceRestPositions.segment<3>(vi * 3) = surfaceMesh.pos(vi);
}
SolidDeformationModel::DeformationModelElasticMaterial parseVolumeElasticMaterial(const pgo::ConfigFileJSON &jconfig)
{
  const std::string material = jconfig.getString("elastic-material");
  if (material == "stable-neo")
    return SolidDeformationModel::DeformationModelElasticMaterial::STABLE_NEO;
  if (material == "stvk-vol")
    return SolidDeformationModel::DeformationModelElasticMaterial::STVK_VOL;
  if (material == "koiter-stvk") {
    throwConfigError(
      "runIPCSim phase2tet/cubic only supports `elastic-material = stable-neo` or `stvk-vol`; `koiter-stvk` is shell-only.");
  }

  throwConfigError(
    "runIPCSim phase2tet/cubic only supports `elastic-material = stable-neo` or `stvk-vol`.");
}
bool parseEnableMaterialMaxStep(const pgo::ConfigFileJSON &jconfig)
{
  return jconfig.exist("enable-material-max-step")
    ? jconfig.getValue<bool>("enable-material-max-step", 1)
    : true;
}
}  // namespace pgo::RunIPCSim
