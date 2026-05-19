#pragma once

#include "EigenDef.h"
#include "ipc/broadPhase/spatialHashGrid.h"

#include <vector>

namespace pgo
{
namespace Contact
{
namespace CIPC
{

// All quantities that are a deterministic function of an obstacle's sampled
// pose (`current_` positions on `ObstacleSurface`). Refreshed by
// `buildObstaclePoseCache` whenever the pose changes.
//
// Spatial AABBs and hashes are stored un-inflated; broad-phase / max-step
// callers are responsible for inflating the dyn-side query box by their own
// margin (Minkowski-sum equivalence), so a single cache serves both pipelines
// regardless of `dhatExternal` / `ccd_thickness`.
struct ObstaclePoseCache
{
  // Barrier-weight derived quantities, consumed by the external barrier
  // assembler.
  std::vector<double> triAreas;
  std::vector<double> edgeLengths;

  // Per-primitive un-inflated AABBs.
  SpatialHashGrid::AABB surfaceBox;
  bool hasSurfaceBox = false;
  std::vector<SpatialHashGrid::AABB> vertBoxes;
  std::vector<SpatialHashGrid::AABB> triBoxes;
  std::vector<SpatialHashGrid::AABB> edgeBoxes;

  // Spatial hashes indexed over the obstacle's triangles / edges using the
  // un-inflated AABBs above.
  SpatialHashGrid triHash{ 0 };
  SpatialHashGrid edgeHash{ 0 };

  // Average obstacle-triangle AABB diagonal (lower-bounded by 1e-6), used as
  // the spatial-hash cell size both for the obstacle hashes above and as a
  // default for any caller-built dyn-side hash that wants to match.
  double cellSize = 0.0;
};

// Refresh every field in `cache` from the supplied pose. Idempotent: vectors
// are resized to match the topology, hashes are cleared and rebuilt, cellSize
// is recomputed. Safe to call on a default-constructed cache.
//
//   positions   : 3*N_verts entries, row-major xyz per vertex
//   triangles   : N_tri x 3 local-index matrix
//   uniqueEdges : N_edge x 2 local-index matrix (derived from triangles)
void buildObstaclePoseCache(
  const EigenSupport::VXd &positions,
  const EigenSupport::MXi &triangles,
  const EigenSupport::MXi &uniqueEdges,
  ObstaclePoseCache &cache);

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
