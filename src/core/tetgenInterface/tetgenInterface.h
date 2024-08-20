#pragma once

#include "EigenDef.h"
#include "triMeshGeo.h"

#include <string>

namespace pgo
{
namespace TetgenInterface
{
struct VoronoiDiagram
{
  // vertices
  EigenSupport::M3Xd vertexPositions;
  // edges
  std::vector<std::tuple<int, int>> edges;
  EigenSupport::M3Xd edgeDirs;
  // facets
  std::vector<std::vector<int>> facets;
  std::vector<std::tuple<int, int>> facetCells;

  std::vector<std::vector<int>> cells;
};

int computeVoronoiDiagram(const std::vector<EigenSupport::V3d> &points, VoronoiDiagram &vd);
int computeTetMesh(const Mesh::TriMeshGeo &mesh, const std::string &switcher, EigenSupport::MXd &vtx, EigenSupport::MXi &tet);
}  // namespace TetgenInterface
}  // namespace pgo