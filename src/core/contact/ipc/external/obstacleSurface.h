#pragma once

#include "EigenDef.h"
#include "ipc/external/obstaclePoseCache.h"

#include <cstdint>
#include <functional>

namespace pgo
{
namespace Contact
{
namespace CIPC
{

class ObstacleSurface
{
public:
  using TrajectorySampler = std::function<void(double t, EigenSupport::RefVecXd out)>;

  ObstacleSurface(
    EigenSupport::MXd restVertices,    // num_obstacle_vertices x 3
    EigenSupport::MXi triangles,       // num_obstacle_tris   x 3, local index
    TrajectorySampler sampler);        // sampler(t, out) writes 3*num_vertices

  // Sample the obstacle pose at absolute time t and rebuild the pose-derived
  // cache (areas, lengths, AABBs, spatial hashes, cell size). Treat the
  // obstacle as fixed at this pose for the subsequent solve; broad-phase /
  // max-step callers should read derived state via `cache()` rather than
  // rebuilding it themselves.
  void update(double t);

  int32_t                  objectId()         const { return objectId_; }
  const EigenSupport::VXd &restPositions()    const { return rest_; }
  const EigenSupport::VXd &currentPositions() const { return current_; }
  const EigenSupport::MXi &triangles()        const { return triangles_; }
  const EigenSupport::MXi &uniqueEdges()      const { return uniqueEdges_; }
  const ObstaclePoseCache &cache()            const { return cache_; }

  void setObjectId(int32_t id) { objectId_ = id; }

private:
  int32_t              objectId_ = -1;
  EigenSupport::VXd    rest_;
  EigenSupport::VXd    current_;
  EigenSupport::MXi    triangles_;       // local 0-based indices
  EigenSupport::MXi    uniqueEdges_;     // derived from triangles_
  TrajectorySampler    sampler_;
  ObstaclePoseCache    cache_;
};

ObstacleSurface::TrajectorySampler makeLinearTrajectorySampler(
  const EigenSupport::VXd &restPositions,
  const EigenSupport::V3d &velocity,
  double t0 = 0.0);

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
