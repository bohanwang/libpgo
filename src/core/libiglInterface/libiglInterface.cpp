#include "libiglInterface.h"

#include "pgoLogging.h"
#include "EigenSupport.h"

#if defined(__GNUC__) || defined(__GNUG__)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wsign-compare"
#  pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#  pragma GCC diagnostic ignored "-Wunknown-pragmas"
#  pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#  pragma GCC diagnostic ignored "-Wunused-parameter"
#elif defined(_MSC_VER)
#  pragma warning(disable : 4244)
#endif

#include <igl/signed_distance.h>
#include <igl/marching_cubes.h>
#include <igl/adjacency_matrix.h>
#include <igl/sum.h>
#include <igl/boundary_loop.h>
#include <igl/massmatrix.h>
#include <igl/cotmatrix.h>
#include <igl/invert_diag.h>

#include <igl/lscm.h>
#include <igl/arap.h>
#include <igl/boundary_loop.h>
#include <igl/harmonic.h>
#include <igl/map_vertices_to_circle.h>

namespace pgo::libiglInterface
{
void triMeshToMatrix(const Mesh::TriMeshGeo &mesh, EigenSupport::MXd &V, EigenSupport::MXi &F)
{
  V.resize(mesh.numVertices(), 3);
  F.resize(mesh.numTriangles(), 3);

  for (int i = 0; i < mesh.numVertices(); i++) {
    V.row(i) = mesh.pos(i);
  }

  for (int i = 0; i < mesh.numTriangles(); i++) {
    F.row(i) = mesh.tri(i);
  }
}

void matrixToTriMesh(const EigenSupport::MXd &V, const EigenSupport::MXi &F, Mesh::TriMeshGeo &mesh)
{
  mesh.clear();
  for (int vi = 0; vi < (int)V.rows(); vi++) {
    mesh.addPos(V.row(vi));
  }

  for (int fi = 0; fi < (int)F.rows(); fi++) {
    mesh.addTri(F.row(fi));
  }
}

}  // namespace pgo::libiglInterface

void pgo::libiglInterface::computeDistanceField(const Mesh::TriMeshGeo &mesh,
  const EigenSupport::V3d &bmin, const EigenSupport::V3d &bmax,
  int res, int robust, int sign, EigenSupport::VXd &dists,
  EigenSupport::VXi *triIDsOut, EigenSupport::MXd *closestPtsOut, EigenSupport::MXd *closestPtNormalsOut)
{
  // Choose type of signing to use
  igl::SignedDistanceType sign_type = igl::SIGNED_DISTANCE_TYPE_PSEUDONORMAL;
  if (robust) {
    sign_type = igl::SIGNED_DISTANCE_TYPE_FAST_WINDING_NUMBER;
  }

  if (sign == 0) {
    sign_type = igl::SIGNED_DISTANCE_TYPE_UNSIGNED;
  }

  EigenSupport::MXd triMeshV;
  EigenSupport::MXi triMeshF;
  triMeshToMatrix(mesh, triMeshV, triMeshF);

  EigenSupport::MX3d gridPoints(res * res * res, 3);
  EigenSupport::V3d diff = bmax - bmin;
  EigenSupport::V3d delta = diff / (res - 1);
  for (int z = 0; z < res; z++) {
    for (int y = 0; y < res; y++) {
      for (int x = 0; x < res; x++) {
        gridPoints.row(z * res * res + y * res + x) = bmin + delta.cwiseProduct(EigenSupport::V3d(x, y, z));
      }
    }
  }

  dists.resize(res * res * res);
  EigenSupport::VXi triIDs(res * res * res);
  EigenSupport::MXd closestPts(res * res * res, 3);
  EigenSupport::MXd closestPtNormals(res * res * res, 3);
  igl::signed_distance(gridPoints, triMeshV, triMeshF, sign_type, dists, triIDs, closestPts, closestPtNormals);

  if (triIDsOut) {
    *triIDsOut = triIDs;
  }

  if (closestPtsOut) {
    *closestPtsOut = closestPts;
  }

  if (closestPtNormalsOut) {
    *closestPtNormalsOut = closestPtNormals;
  }
}

void pgo::libiglInterface::computeMarchingCubes(const EigenSupport::V3d &bmin, const EigenSupport::V3d &bmax, int res, const EigenSupport::VXd &dists, Mesh::TriMeshGeo &outMesh)
{
  EigenSupport::MX3d gridPoints(res * res * res, 3);
  EigenSupport::V3d diff = bmax - bmin;
  EigenSupport::V3d delta = diff / (res - 1);
  for (int z = 0; z < res; z++) {
    for (int y = 0; y < res; y++) {
      for (int x = 0; x < res; x++) {
        gridPoints.row(z * res * res + y * res + x) = bmin + delta.cwiseProduct(EigenSupport::V3d(x, y, z));
      }
    }
  }

  EigenSupport::MXd V;
  EigenSupport::MXi F;
  igl::marching_cubes(dists, gridPoints, res, res, res, 0.0, V, F);

  matrixToTriMesh(V, F, outMesh);
}

namespace pgo::libiglInterface
{
}  // namespace pgo::libiglInterface

void pgo::libiglInterface::graphLaplacianMatrix(const EigenSupport::MXi &F, EigenSupport::SpMatD &K)
{
  Eigen::SparseMatrix<double> A_col;
  igl::adjacency_matrix(F, A_col);
  EigenSupport::SpMatD A = A_col;

  Eigen::SparseVector<double> Asum;
  igl::sum(A_col, 1, Asum);

  EigenSupport::SpMatD Adiag;
  diag(Asum, Adiag);

  K = A - Adiag;
}

void pgo::libiglInterface::graphLaplacianMatrix(const EigenSupport::SpMatD &Adj, EigenSupport::SpMatD &K)
{
  Eigen::SparseMatrix<double> A_col = Adj;
  Eigen::SparseVector<double> Asum;
  igl::sum(A_col, 1, Asum);

  EigenSupport::SpMatD Adiag;
  diag(Asum, Adiag);

  K = Adj - Adiag;
}

void pgo::libiglInterface::diag(const EigenSupport::VXd &V, EigenSupport::SpMatD &X)
{
  PGO_ALOG(V.rows() == 1 || V.cols() == 1);

  // clear and resize output
  X.setZero();
  X.resize(V.size(), V.size());

  std::vector<EigenSupport::TripletD> entries;
  entries.reserve(V.size());

  // loop over non-zeros
  for (int i = 0; i < V.size(); i++) {
    entries.emplace_back(i, i, V[i]);
  }
  X.setFromTriplets(entries.begin(), entries.end());
}

void pgo::libiglInterface::diag(const Eigen::SparseVector<double> &V, EigenSupport::SpMatD &X)
{
  // clear and resize output
  X.setZero();
  X.resize(V.size(), V.size());

  std::vector<EigenSupport::TripletD> entries;
  entries.reserve(V.size());

  // loop over non-zeros
  for (Eigen::SparseVector<double>::InnerIterator it(V); it; ++it) {
    entries.emplace_back(it.index(), it.index(), it.value());
  }

  X.setFromTriplets(entries.begin(), entries.end());
}

void pgo::libiglInterface::boundaryLoop(const Mesh::TriMeshGeo &meshIn, std::vector<int> &boundary)
{
  EigenSupport::MXd vtx;
  EigenSupport::MXi tri;
  Mesh::triMeshGeoToMatrices(meshIn, vtx, tri);

  igl::boundary_loop(tri, boundary);
}

void pgo::libiglInterface::boundaryLoops(const Mesh::TriMeshGeo &meshIn, std::vector<std::vector<int>> &boundary)
{
  EigenSupport::MXd vtx;
  EigenSupport::MXi tri;

  Mesh::triMeshGeoToMatrices(meshIn, vtx, tri);

  igl::boundary_loop(tri, boundary);
}

void pgo::libiglInterface::meanCuravtures(const Mesh::TriMeshGeo &meshIn, EigenSupport::VXd &h)
{
  EigenSupport::MXd V;
  EigenSupport::MXi F;
  Mesh::triMeshGeoToMatrices(meshIn, V, F);

  EigenSupport::MXd HN;
  Eigen::SparseMatrix<double> L, M, Minv;

  igl::cotmatrix(V, F, L);
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M);
  igl::invert_diag(M, Minv);
  HN = -Minv * (L * V);
  h = HN.rowwise().norm();  // up to sign
}

void pgo::libiglInterface::computeParameterization(const Mesh::TriMeshGeo &meshIn, Mesh::TriMeshGeo &meshOut, int mode)
{
  EigenSupport::MXd V;
  EigenSupport::MXi F;
  Mesh::triMeshGeoToMatrices(meshIn, V, F);

  Eigen::MatrixXd V_uv;

  if (mode == 0) {
    Eigen::MatrixXd initial_guess;

    // Compute the initial solution for ARAP (harmonic parametrization)
    Eigen::VectorXi bnd;
    igl::boundary_loop(F, bnd);

    Eigen::MatrixXd bnd_uv;
    igl::map_vertices_to_circle(V, bnd, bnd_uv);

    igl::harmonic(V, F, bnd, bnd_uv, 1, initial_guess);

    std::cout << 1 << std::endl;

    // Add dynamic regularization to avoid to specify boundary conditions
    igl::ARAPData arap_data;
    arap_data.with_dynamics = true;

    Eigen::VectorXi b = Eigen::VectorXi::Zero(0);
    Eigen::MatrixXd bc = Eigen::MatrixXd::Zero(0, 0);

    // Initialize ARAP
    arap_data.max_iter = 100;

    std::cout << 2 << std::endl;

    // 2 means that we're going to *solve* in 2d
    arap_precomputation(V, F, 2, b, arap_data);

    std::cout << 3 << std::endl;

    // Solve arap using the harmonic map as initial guess
    V_uv = initial_guess;

    arap_solve(bc, arap_data, V_uv);

    std::cout << 4 << std::endl;
  }
  else if (mode == 1) {
    igl::lscm(V, F, V_uv);
  }

  meshOut = meshIn;
  for (int vi = 0; vi < meshOut.numVertices(); vi++) {
    meshOut.pos(vi) = EigenSupport::V3d(V_uv(vi, 0), V_uv(vi, 1), 0);
  }
}

void pgo::libiglInterface::computeMassMatrix(const Mesh::TriMeshGeo &meshIn, EigenSupport::SpMatD &M, int lumped, int expand3)
{
  EigenSupport::MXd V;
  EigenSupport::MXi F;
  Mesh::triMeshGeoToMatrices(meshIn, V, F);

  Eigen::SparseMatrix<double> MTemp;
  igl::massmatrix(V, F, lumped ? igl::MASSMATRIX_TYPE_BARYCENTRIC : igl::MASSMATRIX_TYPE_VORONOI, MTemp);

  PGO_ALOG(MTemp.rows() == V.rows());

  if (expand3) {
    EigenSupport::SpMatD Mrow = MTemp;
    EigenSupport::expand3(Mrow, M);
  }
  else {
    M = MTemp;
  }
}