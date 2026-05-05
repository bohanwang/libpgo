#pragma once

#include "EigenSupport.h"
#include "generateFBMSUnionSurfaceSphere.h"
#include "triMeshGeo.h"

#include <string>

namespace pgo::Tools::FBMSUnionSurface
{

struct ThicknessSearchResult;

void extractSurfaceMarchingCubes(const EigenSupport::V3d &bmin, const EigenSupport::V3d &bmax,
  int resolution, const EigenSupport::VXd &field, double isoOffset, Mesh::TriMeshGeo &outMesh);

void extractSurfaceOpenVDB(const Mesh::TriMeshGeo &fbmsMesh, const SphereParameters &sphere,
  double fbmsThickness, double sphereThickness, bool enableTruncating, const std::string &debugFieldMode,
  double voxelSize, double halfWidth, double adaptivity, int smoothSteps,
  Mesh::TriMeshGeo &outMesh);

ThicknessSearchResult findThicknessForVolumeBudgetOpenVDB(
  const Mesh::TriMeshGeo &fbmsMesh, const SphereParameters &sphere,
  double sphereThickness, bool enableTruncating,
  double voxelSize, double halfWidth, double adaptivity, int smoothSteps,
  double targetVolume, double tLo, double tHi, int maxIterations, double relTol);

}  // namespace pgo::Tools::FBMSUnionSurface
