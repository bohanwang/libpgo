#include "obstacleSurface.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

namespace pgo
{
namespace Contact
{
namespace IPC
{
namespace
{

constexpr double kCoplanarNormalTolerance = 1e-12;

std::pair<int, int> canonicalEdge(int a, int b)
{
  if (a > b)
    std::swap(a, b);
  return { a, b };
}

EigenSupport::MXi edgeMatrixFromSortedEdges(const std::vector<std::pair<int, int>> &edges)
{
  EigenSupport::MXi edgeMatrix(static_cast<int>(edges.size()), 2);
  for (int ei = 0; ei < static_cast<int>(edges.size()); ++ei) {
    edgeMatrix(ei, 0) = edges[ei].first;
    edgeMatrix(ei, 1) = edges[ei].second;
  }
  return edgeMatrix;
}

EigenSupport::MXi buildUniqueEdgeMatrix(const EigenSupport::MXi &triangles)
{
  std::set<std::pair<int, int>> edgeSet;
  for (int fi = 0; fi < triangles.rows(); ++fi) {
    for (int j = 0; j < 3; ++j)
      edgeSet.insert(canonicalEdge(triangles(fi, j), triangles(fi, (j + 1) % 3)));
  }

  std::vector<std::pair<int, int>> edges(edgeSet.begin(), edgeSet.end());
  return edgeMatrixFromSortedEdges(edges);
}

EigenSupport::MXi buildContactEdgeMatrix(
  const EigenSupport::VXd &positions,
  const EigenSupport::MXi &triangles)
{
  std::map<std::pair<int, int>, std::vector<int>> incidentFaces;
  for (int fi = 0; fi < triangles.rows(); ++fi) {
    for (int j = 0; j < 3; ++j)
      incidentFaces[canonicalEdge(triangles(fi, j), triangles(fi, (j + 1) % 3))].push_back(fi);
  }

  std::vector<EigenSupport::V3d> normals(static_cast<std::size_t>(triangles.rows()), EigenSupport::V3d::Zero());
  std::vector<bool> degenerate(static_cast<std::size_t>(triangles.rows()), false);
  for (int fi = 0; fi < triangles.rows(); ++fi) {
    const EigenSupport::V3d v0 = positions.segment<3>(3 * triangles(fi, 0));
    const EigenSupport::V3d v1 = positions.segment<3>(3 * triangles(fi, 1));
    const EigenSupport::V3d v2 = positions.segment<3>(3 * triangles(fi, 2));
    const EigenSupport::V3d normal = (v1 - v0).cross(v2 - v0);
    const double normalLength = normal.norm();
    if (normalLength <= 0.0) {
      degenerate[static_cast<std::size_t>(fi)] = true;
      continue;
    }
    normals[static_cast<std::size_t>(fi)] = normal / normalLength;
  }

  std::vector<std::pair<int, int>> contactEdges;
  contactEdges.reserve(incidentFaces.size());
  for (const auto &[edge, faces] : incidentFaces) {
    bool keep = faces.size() != 2;
    if (!keep) {
      const int f0 = faces[0];
      const int f1 = faces[1];
      if (degenerate[static_cast<std::size_t>(f0)] || degenerate[static_cast<std::size_t>(f1)]) {
        keep = true;
      }
      else {
        const double absDot = std::min(1.0, std::abs(normals[static_cast<std::size_t>(f0)].dot(normals[static_cast<std::size_t>(f1)])));
        keep = (1.0 - absDot) > kCoplanarNormalTolerance;
      }
    }

    if (keep)
      contactEdges.push_back(edge);
  }

  return edgeMatrixFromSortedEdges(contactEdges);
}

}  // namespace

ObstacleSurface::ObstacleSurface(
  EigenSupport::MXd restVertices,
  EigenSupport::MXi triangles,
  TrajectorySampler sampler):
  triangles_(std::move(triangles)),
  sampler_(std::move(sampler))
{
  if (restVertices.cols() != 3)
    throw std::invalid_argument("ObstacleSurface: restVertices must be an N x 3 vertex matrix.");
  if (triangles_.cols() != 3)
    throw std::invalid_argument("ObstacleSurface: triangles must be an N x 3 matrix.");
  if (!sampler_)
    throw std::invalid_argument("ObstacleSurface: trajectory sampler must not be empty.");

  if (triangles_.size() > 0) {
    if (triangles_.minCoeff() < 0 || triangles_.maxCoeff() >= restVertices.rows())
      throw std::invalid_argument("ObstacleSurface: triangle index out of range.");
  }

  const int numVerts = static_cast<int>(restVertices.rows());
  rest_.resize(numVerts * 3);
  for (int vi = 0; vi < numVerts; ++vi)
    rest_.segment<3>(3 * vi) = restVertices.row(vi).transpose();

  uniqueEdges_ = buildUniqueEdgeMatrix(triangles_);
  contactEdges_ = buildContactEdgeMatrix(rest_, triangles_);

  current_.resize(numVerts * 3);
  current_.setZero();
}

void ObstacleSurface::update(double t)
{
  sampler_(t, current_);
  contactEdges_ = buildContactEdgeMatrix(current_, triangles_);
  buildObstaclePoseCache(current_, triangles_, contactEdges_, cache_);
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

}  // namespace IPC
}  // namespace Contact
}  // namespace pgo
