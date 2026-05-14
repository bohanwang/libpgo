#pragma once

#include "EigenSupport.h"
#include "triMeshGeo.h"
#include "field/denseGrid.h"

namespace pgo::ImplicitSurface {

struct SphereField {
  EigenSupport::V3d center = EigenSupport::V3d::Zero();
  double radius = 0.0;
};

SphereField computeSphereFieldFromBBox(const Mesh::TriMeshGeo &sphereMesh);
int projectOpenBoundaryToSphere(Mesh::TriMeshGeo &mesh, const SphereField &sphere);

// Evaluate sphere shell field: | |p - center| - radius | - 0.5 * thickness
void thickenSphereShell(const SphereField &sphere, double thickness,
  const GridSpec &grid, DenseGrid &out);

// Evaluate solid ball SDF: |p - center| - radius
void evaluateBallSDF(const SphereField &sphere, const GridSpec &grid, DenseGrid &out);

}  // namespace pgo::ImplicitSurface
