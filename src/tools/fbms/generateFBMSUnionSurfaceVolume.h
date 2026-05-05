#pragma once

#include "EigenSupport.h"
#include "generateFBMSUnionSurfaceSphere.h"
#include "triMeshGeo.h"

#include <string>

namespace pgo::Tools::FBMSUnionSurface
{

double computeMeshVolume(const Mesh::TriMeshGeo &mesh);

struct ThicknessSearchResult
{
  double thickness = 0.0;
  double volume = 0.0;
  Mesh::TriMeshGeo mesh;
  bool budgetExceededAtMin = false;
  bool budgetUnreachedAtMax = false;
  int iterations = 0;
};

ThicknessSearchResult findThicknessForVolumeBudget(
  const EigenSupport::VXd &fbmsDistance, const SphereParameters &sphere,
  const EigenSupport::V3d &bmin, const EigenSupport::V3d &bmax, int resolution,
  double sphereThickness, bool enableTruncating, const std::string &isoOffsetMode, double fixedIsoOffset,
  double targetVolume,
  double tLo, double tHi, int maxIterations, double relTol);

}  // namespace pgo::Tools::FBMSUnionSurface
