#pragma once

#include "EigenSupport.h"
#include "generateFBMSUnionSurfaceOptions.h"
#include "triMeshGeo.h"

#include <string>

namespace pgo::Tools::FBMSUnionSurface
{

void loadMeshOrThrow(const std::string &path, Mesh::TriMeshGeo &mesh, const std::string &label);
void printVector(const char *label, const EigenSupport::V3d &v);
void saveSurfaceOrThrow(const Mesh::TriMeshGeo &mesh, const std::string &path);
void applySmallComponentFilterIfRequested(Mesh::TriMeshGeo &rawSurface, const Options &options);

}  // namespace pgo::Tools::FBMSUnionSurface
