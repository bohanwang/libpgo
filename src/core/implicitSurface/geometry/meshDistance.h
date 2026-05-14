#pragma once

#include "EigenSupport.h"
#include "triMeshGeo.h"
#include "field/gridSpec.h"
#include "field/denseGrid.h"

namespace pgo::ImplicitSurface {

// Compute unsigned distance field from a surface mesh onto a dense grid.
void computeMeshUnsignedDistance(const Mesh::TriMeshGeo &surfaceMesh,
  const GridSpec &grid, DenseGrid &outDistance);

}  // namespace pgo::ImplicitSurface
