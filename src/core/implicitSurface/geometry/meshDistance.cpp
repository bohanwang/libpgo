#include "geometry/meshDistance.h"

#include "libiglInterface.h"

#include <stdexcept>

namespace pgo::ImplicitSurface {

void computeMeshUnsignedDistance(const Mesh::TriMeshGeo &surfaceMesh,
  const GridSpec &grid, DenseGrid &outDistance)
{
  if (outDistance.gridSpec() != grid)
    throw std::runtime_error("computeMeshUnsignedDistance: outDistance grid spec does not match input grid");

  EigenSupport::VXd dist;
  libiglInterface::computeDistanceField(
    surfaceMesh, grid.bmin, grid.bmax, grid.resolution,
    /*robust=*/1, /*sign=*/0, dist);

  const int total = outDistance.size();
  for (int i = 0; i < total; ++i)
    outDistance[i] = dist[i];
}

}  // namespace pgo::ImplicitSurface
