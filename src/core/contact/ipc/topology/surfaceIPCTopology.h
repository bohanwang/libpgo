/*
copyright to Bohan Wang
*/

#pragma once

#include "EigenDef.h"

#include <array>
#include <cstdint>
#include <vector>

namespace pgo
{
namespace Contact
{
namespace CIPC
{

struct SurfaceIPCTopology
{
  int numVerts = 0;
  std::vector<std::array<int, 3>> triangles;
  std::vector<std::array<int, 2>> edges;
  std::vector<double> vertexArea;
  std::vector<double> triArea;
  std::vector<double> edgeLength;
  std::vector<uint8_t> vertexIsDeformable;
  std::vector<uint8_t> triangleHasDeformableVertex;
  std::vector<uint8_t> edgeHasDeformableVertex;

  void setMesh(const EigenSupport::MXd &V, const EigenSupport::MXi &F);
  void setMesh(const EigenSupport::MXd &V, const EigenSupport::MXi &F,
    const std::vector<uint8_t> &vertexIsDeformableMask);
  // Recompute triangle areas, vertex areas and edge lengths at the given
  // flattened positions (3 * numVerts). Used when kinematic external vertices
  // are moved or deformed between frames so the barrier and friction weights
  // follow the current geometry.
  void updateMeasures(EigenSupport::ConstRefVecXd positions);
  int numSurfaceDOFs() const { return 3 * numVerts; }
  bool isVertexDeformable(int vi) const { return vertexIsDeformable[vi] != 0; }
  bool triangleContainsDeformableVertex(int fi) const { return triangleHasDeformableVertex[fi] != 0; }
  bool edgeContainsDeformableVertex(int ei) const { return edgeHasDeformableVertex[ei] != 0; }
};

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
