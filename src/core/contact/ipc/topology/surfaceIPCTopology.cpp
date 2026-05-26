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
namespace IPC
{

void SurfaceIPCTopology::setMesh(const EigenSupport::MXd &V, const EigenSupport::MXi &F)
{
  if (V.cols() != 3)
    throw std::invalid_argument("SurfaceIPCTopology: V must be an N x 3 vertex matrix.");
  if (F.cols() != 3)
    throw std::invalid_argument("SurfaceIPCTopology: F must be an N x 3 triangle index matrix.");
  if (F.size() > 0 && (F.minCoeff() < 0 || F.maxCoeff() >= V.rows()))
    throw std::invalid_argument("SurfaceIPCTopology: F contains an out-of-range vertex index.");

  numVerts = static_cast<int>(V.rows());
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

}  // namespace IPC
}  // namespace Contact
}  // namespace pgo
