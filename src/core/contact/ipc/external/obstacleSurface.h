#pragma once

#include "EigenDef.h"
#include "ipc/external/obstacleSurfaceView.h"

#include <cstdint>
#include <functional>
#include <vector>

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

  void update(double tStart, double tEnd);

  int32_t                          objectId()        const { return objectId_; }
  const EigenSupport::VXd &        restPositions()   const { return rest_; }
  const EigenSupport::VXd &        currentPositions()const { return current_; }
  const EigenSupport::VXd &        previousPositions()const{ return previous_; }
  const EigenSupport::MXi &        triangles()       const { return triangles_; }
  const EigenSupport::MXi &        uniqueEdges()     const { return uniqueEdges_; }
  const std::vector<double> &       triAreas()        const { return triAreas_; }
  const std::vector<double> &       edgeLengths()     const { return edgeLengths_; }

  void setObjectId(int32_t id) { objectId_ = id; }
  ObstacleSurfaceView view() const;

private:
  int32_t              objectId_ = -1;
  EigenSupport::VXd    rest_;
  EigenSupport::VXd    previous_;
  EigenSupport::VXd    current_;
  EigenSupport::MXi    triangles_;       // local 0-based indices
  EigenSupport::MXi    uniqueEdges_;     // derived from triangles_
  TrajectorySampler    sampler_;
  std::vector<double>  triAreas_;        // cached from current_ positions
  std::vector<double>  edgeLengths_;     // cached from current_ positions
};

ObstacleSurface::TrajectorySampler makeLinearTrajectorySampler(
  const EigenSupport::VXd &restPositions,
  const EigenSupport::V3d &velocity,
  double t0 = 0.0);

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
