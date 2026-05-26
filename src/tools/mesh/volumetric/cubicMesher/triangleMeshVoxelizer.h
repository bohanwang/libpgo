#pragma once

#include "cubicMesh.h"

#include <memory>
#include <string>

namespace cubic_mesher
{

struct TriangleMeshVoxelizerOptions
{
  std::string inputMesh;
  int resolution = 0;
  double E = 1e6;
  double nu = 0.45;
  double density = 1000.0;
};

std::unique_ptr<pgo::VolumetricMeshes::CubicMesh> createTriangleMeshCubicMesh(
  const TriangleMeshVoxelizerOptions &options);

}  // namespace cubic_mesher
