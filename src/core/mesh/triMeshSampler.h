#pragma once

#include "triMeshGeo.h"

#include <vector>
#include <array>

namespace pgo
{
namespace Mesh
{
int sampleTriangle(const TriMeshGeo &mesh, int subdivideTriangleCount,
  std::vector<int> &sampleTriangleIDs, std::vector<double> &sampleWeights, std::vector<std::array<double, 3>> &sampleBarycentricWeights);

}
}  // namespace pgo
