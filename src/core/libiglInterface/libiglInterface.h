#pragma once

#include "triMeshGeo.h"
#include "EigenDef.h"

namespace pgo
{
namespace libiglInterface
{
void computeDistanceField(const Mesh::TriMeshGeo &mesh,
  const EigenSupport::V3d &bmin, const EigenSupport::V3d &bmax,
  int res, int robust, int sign,
  EigenSupport::VXd &dists,
  EigenSupport::VXi *triIDsOut = nullptr, EigenSupport::MXd *closestPtsOut = nullptr, EigenSupport::MXd *closestPtNormalsOut = nullptr);

void computeMarchingCubes(const EigenSupport::V3d &bmin, const EigenSupport::V3d &bmax, int res, const EigenSupport::VXd &dists, Mesh::TriMeshGeo &outMesh);

void diag(const EigenSupport::VXd &v, EigenSupport::SpMatD &M);
void diag(const Eigen::SparseVector<double> &v, EigenSupport::SpMatD &M);

void graphLaplacianMatrix(const EigenSupport::MXi &F, EigenSupport::SpMatD &K);
void graphLaplacianMatrix(const EigenSupport::SpMatD &Adj, EigenSupport::SpMatD &K);

void boundaryLoop(const Mesh::TriMeshGeo &meshIn, std::vector<int> &boundary);
void boundaryLoops(const Mesh::TriMeshGeo &meshIn, std::vector<std::vector<int>> &boundary);

void meanCuravtures(const Mesh::TriMeshGeo &meshIn, EigenSupport::VXd &h);

void computeParameterization(const Mesh::TriMeshGeo &meshIn, Mesh::TriMeshGeo &meshOut, int mode);

void computeMassMatrix(const Mesh::TriMeshGeo &meshIn, EigenSupport::SpMatD &M, int lumped, int expand3);

}  // namespace libiglInterface
}  // namespace pgo