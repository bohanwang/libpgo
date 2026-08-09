/*
copyright to Bohan Wang
*/

#include "surfaceIPCTopology.h"

#include <set>
#include <stdexcept>

namespace pgo
{
namespace Contact
{
namespace CIPC
{

void SurfaceIPCTopology::setMesh(const EigenSupport::MXd &V, const EigenSupport::MXi &F)
{
  setMesh(V, F, std::vector<uint8_t>(static_cast<std::size_t>(V.rows()), uint8_t{ 1 }));
}

void SurfaceIPCTopology::setMesh(const EigenSupport::MXd &V, const EigenSupport::MXi &F,
  const std::vector<uint8_t> &vertexIsDeformableMask)
{
  numVerts = static_cast<int>(V.rows());
  if (vertexIsDeformableMask.size() != static_cast<std::size_t>(numVerts))
    throw std::invalid_argument("SurfaceIPCTopology deformable-role mask length must equal the number of mesh vertices.");
  for (uint8_t value : vertexIsDeformableMask) {
    if (value > 1)
      throw std::invalid_argument("SurfaceIPCTopology deformable-role mask entries must be 0 or 1.");
  }
  vertexIsDeformable = vertexIsDeformableMask;

  triangles.resize(F.rows());
  for (int i = 0; i < static_cast<int>(F.rows()); ++i)
    triangles[i] = { F(i, 0), F(i, 1), F(i, 2) };

  std::set<std::pair<int, int>> edgeSet;
  for (const auto &tri : triangles) {
    for (int j = 0; j < 3; ++j) {
      int a = tri[j], b = tri[(j + 1) % 3];
      if (a > b)
        std::swap(a, b);
      edgeSet.insert({ a, b });
    }
  }
  edges.clear();
  edges.reserve(edgeSet.size());
  for (const auto &e : edgeSet)
    edges.push_back({ e.first, e.second });

  const int nTri = static_cast<int>(triangles.size());
  const int nEdge = static_cast<int>(edges.size());

  triangleHasDeformableVertex.resize(nTri);
  for (int fi = 0; fi < nTri; ++fi) {
    const auto &tri = triangles[fi];
    triangleHasDeformableVertex[fi] = static_cast<uint8_t>(
      isVertexDeformable(tri[0]) || isVertexDeformable(tri[1]) || isVertexDeformable(tri[2]));
  }

  edgeHasDeformableVertex.resize(nEdge);
  for (int ei = 0; ei < nEdge; ++ei) {
    const auto &edge = edges[ei];
    edgeHasDeformableVertex[ei] = static_cast<uint8_t>(
      isVertexDeformable(edge[0]) || isVertexDeformable(edge[1]));
  }

  triArea.resize(nTri);
  for (int fi = 0; fi < nTri; ++fi) {
    EigenSupport::V3d v0 = V.row(triangles[fi][0]).transpose();
    EigenSupport::V3d v1 = V.row(triangles[fi][1]).transpose();
    EigenSupport::V3d v2 = V.row(triangles[fi][2]).transpose();
    triArea[fi] = 0.5 * (v1 - v0).cross(v2 - v0).norm();
  }

  vertexArea.assign(numVerts, 0.0);
  for (int fi = 0; fi < nTri; ++fi) {
    const double areaThird = triArea[fi] / 3.0;
    for (int j = 0; j < 3; ++j)
      vertexArea[triangles[fi][j]] += areaThird;
  }

  edgeLength.resize(nEdge);
  for (int ei = 0; ei < nEdge; ++ei) {
    EigenSupport::V3d e0 = V.row(edges[ei][0]).transpose();
    EigenSupport::V3d e1 = V.row(edges[ei][1]).transpose();
    edgeLength[ei] = (e1 - e0).norm();
  }
}

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
