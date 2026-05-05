#pragma once

#include "EigenSupport.h"
#include "generateFBMSUnionSurfaceSphere.h"
#include "triMeshGeo.h"

#include <string>

namespace pgo::Tools::FBMSUnionSurface
{

void computeUnionBBox(const Mesh::TriMeshGeo &fbmsMesh, const Mesh::TriMeshGeo &sphereMesh,
  double fbmsThickness, double sphereThickness, double paddingRatio,
  EigenSupport::V3d &bmin, EigenSupport::V3d &bmax);

void computeFBMSDistance(const Mesh::TriMeshGeo &fbmsMesh,
  const EigenSupport::V3d &bmin, const EigenSupport::V3d &bmax, int resolution,
  EigenSupport::VXd &fbmsDistance);

void assembleUnionField(const EigenSupport::VXd &fbmsDistance, const SphereParameters &sphere,
  const EigenSupport::V3d &bmin, const EigenSupport::V3d &bmax, int resolution,
  double fbmsThickness, double sphereThickness, bool enableTruncating, const std::string &debugFieldMode,
  EigenSupport::VXd &unionField, double &fieldMin, double &fieldMax);

double selectIsoOffset(const std::string &isoOffsetMode, double fixedIsoOffset,
  const EigenSupport::VXd &field, double fbmsThickness, double sphereThickness,
  const EigenSupport::V3d &bmin, const EigenSupport::V3d &bmax, int resolution);

}  // namespace pgo::Tools::FBMSUnionSurface
