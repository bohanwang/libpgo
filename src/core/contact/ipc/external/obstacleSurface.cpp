#include "obstacleSurface.h"

#include <set>
#include <stdexcept>

namespace pgo
{
namespace Contact
{
namespace CIPC
{

ObstacleSurface::ObstacleSurface(
  EigenSupport::MXd restVertices,
  EigenSupport::MXi triangles,
  TrajectorySampler sampler):
  rest_(restVertices.reshaped()),
  triangles_(std::move(triangles)),
  sampler_(std::move(sampler))
{
  if (triangles_.cols() != 3)
    throw std::invalid_argument("ObstacleSurface: triangles must be an N x 3 matrix.");
  if (!sampler_)
    throw std::invalid_argument("ObstacleSurface: trajectory sampler must not be empty.");

  if (triangles_.size() > 0) {
    if (triangles_.minCoeff() < 0 || triangles_.maxCoeff() >= restVertices.rows())
      throw std::invalid_argument("ObstacleSurface: triangle index out of range.");
  }

  const int numVerts = static_cast<int>(restVertices.rows());

  // Derive unique edges from triangles
  std::set<std::pair<int, int>> edgeSet;
  for (int fi = 0; fi < triangles_.rows(); ++fi) {
    for (int j = 0; j < 3; ++j) {
      int a = triangles_(fi, j);
      int b = triangles_(fi, (j + 1) % 3);
      if (a > b)
        std::swap(a, b);
      edgeSet.insert({ a, b });
    }
  }
  uniqueEdges_.resize(static_cast<int>(edgeSet.size()), 2);
  int row = 0;
  for (const auto &e : edgeSet) {
    uniqueEdges_(row, 0) = e.first;
    uniqueEdges_(row, 1) = e.second;
    ++row;
  }

  previous_.resize(numVerts * 3);
  current_.resize(numVerts * 3);
  previous_.setZero();
  current_.setZero();
}

void ObstacleSurface::update(double tStart, double tEnd)
{
  const int n3 = static_cast<int>(rest_.size());
  sampler_(tStart, previous_);
  sampler_(tEnd, current_);
}

ObstacleSurface::TrajectorySampler makeLinearTrajectorySampler(
  const EigenSupport::VXd &restPositions,
  const EigenSupport::V3d &velocity,
  double t0)
{
  const int n3 = static_cast<int>(restPositions.size());
  return [restPositions, velocity, t0, n3](double t, EigenSupport::RefVecXd out) {
    if (out.size() != n3)
      throw std::runtime_error("makeLinearTrajectorySampler: output vector has wrong size.");
    for (int i = 0; i < n3 / 3; ++i)
      out.segment<3>(3 * i) = restPositions.segment<3>(3 * i) + velocity * (t - t0);
  };
}

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
