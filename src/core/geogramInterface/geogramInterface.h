#pragma once

#include "triMeshGeo.h"

namespace pgo
{
namespace GeogramInterface
{
void initGEO();
Mesh::TriMeshGeo remesh(const std::string &meshFilename, int targetNumPoints, double sizeFactor, double anisotropy);
Mesh::TriMeshGeo remesh(const Mesh::TriMeshGeo &inputMesh, int targetNumPoints, double sizeFactor, double anisotropy);

}  // namespace GeogramInterface
}  // namespace pgo