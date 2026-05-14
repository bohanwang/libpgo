#pragma once

#include "EigenSupport.h"
#include "triMeshGeo.h"
#include "geometry/sphereField.h"
#include "operations/booleanOps.h"

#include <memory>
#include <stdexcept>

#ifdef PGO_HAS_OPENVDB
#include <openvdb/openvdb.h>
#include <openvdb/Grid.h>
#endif

namespace pgo::ImplicitSurface {

#ifdef PGO_HAS_OPENVDB
struct OpenVDBLevelSet {
  openvdb::FloatGrid::Ptr grid = nullptr;
  OpenVDBLevelSet() = default;
  explicit OpenVDBLevelSet(openvdb::FloatGrid::Ptr g)
    : grid(std::move(g)) {}
};
#else
struct OpenVDBLevelSet {
  OpenVDBLevelSet() = default;
};
#endif

struct OpenVDBOptions {
  double voxelSize = 0.0;
  double halfWidth = 3.0;
  double adaptivity = 0.0;
  int smoothSteps = 0;
};

void validateOpenVDBOptions(const OpenVDBOptions &options);

// Build shell level-set from an arbitrary surface mesh.
std::unique_ptr<OpenVDBLevelSet> buildOpenVDBShellFromMesh(
  const Mesh::TriMeshGeo &surfaceMesh, double shellThickness,
  const OpenVDBOptions &options);

// Build sphere shell level-set analytically.
std::unique_ptr<OpenVDBLevelSet> buildOpenVDBSphereShell(
  const SphereField &sphere, double sphereShellThickness,
  const OpenVDBOptions &options);

// Build solid ball SDF (for truncation).
std::unique_ptr<OpenVDBLevelSet> buildOpenVDBBallLevelSet(
  const SphereField &sphere, const OpenVDBOptions &options);

// CSG combine two level sets.
std::unique_ptr<OpenVDBLevelSet> combineOpenVDBLevelSets(
  const OpenVDBLevelSet &a, const OpenVDBLevelSet &b, BooleanOp op);

// Extract mesh from level set.
void extractOpenVDBLevelSet(const OpenVDBLevelSet &levelSet,
  const OpenVDBOptions &options, Mesh::TriMeshGeo &outMesh);

}  // namespace pgo::ImplicitSurface
