#pragma once

#include "EigenSupport.h"
#include "deformationModelManager.h"
#include "triMeshGeo.h"

#include <string>

namespace pgo
{
class ConfigFileJSON;
}

namespace pgo::RunIPCSim
{
[[noreturn]] void throwConfigError(const std::string &message);
void rejectIfPresent(const pgo::ConfigFileJSON &config, const char *field, const char *reason);
EigenSupport::SpMatD makeIdentityEmbedding(int n3);
void validateZeroInitialDisplacement(const pgo::ConfigFileJSON &jconfig);
void loadSurfaceMeshAndRestPositions(const std::string &surfaceMeshFilename, double scale,
  pgo::Mesh::TriMeshGeo &surfaceMesh, EigenSupport::VXd &surfaceRestPositions);
SolidDeformationModel::DeformationModelElasticMaterial parseVolumeElasticMaterial(const pgo::ConfigFileJSON &jconfig);
bool parseEnableMaterialMaxStep(const pgo::ConfigFileJSON &jconfig);
}  // namespace pgo::RunIPCSim
