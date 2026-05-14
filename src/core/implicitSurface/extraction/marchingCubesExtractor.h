#pragma once

#include "triMeshGeo.h"
#include "field/denseGrid.h"

namespace pgo::ImplicitSurface {

struct MarchingCubesOptions {
  double isoOffset = 0.0;
};

void extractMarchingCubes(const DenseGrid &field, const MarchingCubesOptions &options,
  Mesh::TriMeshGeo &outMesh);

}  // namespace pgo::ImplicitSurface
