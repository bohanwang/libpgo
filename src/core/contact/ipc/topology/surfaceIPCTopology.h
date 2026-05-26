/*
copyright to Bohan Wang
*/

#pragma once

#include "EigenDef.h"

#include <array>
#include <vector>

namespace pgo
{
namespace Contact
{
namespace IPC
{

struct SurfaceIPCTopology
{
  int numVerts = 0;
  std::vector<std::array<int, 3>> triangles;
  std::vector<std::array<int, 2>> edges;
  std::vector<double> vertexArea;
  std::vector<double> triArea;
  std::vector<double> edgeLength;

  void setMesh(const EigenSupport::MXd &V, const EigenSupport::MXi &F);
  int numSurfaceDOFs() const { return 3 * numVerts; }
};

}  // namespace IPC
}  // namespace Contact
}  // namespace pgo
