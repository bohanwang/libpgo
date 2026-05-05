#pragma once

#include "EigenSupport.h"
#include "triMeshGeo.h"

#include <string>

namespace pgo::Tools::FBMSUnionSurface
{

struct SphereParameters
{
  EigenSupport::V3d center = EigenSupport::V3d::Zero();
  double radius = 0.0;
  bool fromHeader = false;
};

SphereParameters loadSphereParameters(const std::string &spherePath, const Mesh::TriMeshGeo &sphereMesh);
int projectBoundaryVerticesToSphere(Mesh::TriMeshGeo &mesh, const SphereParameters &sphere);

}  // namespace pgo::Tools::FBMSUnionSurface
