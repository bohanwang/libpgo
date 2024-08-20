#include "surfaceTriangleDeformation.h"
#include "pgoLogging.h"
#include "triMeshNeighbor.h"
#include "basicAlgorithms.h"
#include "geometryQuery.h"
#include "determinantDerivatives.h"
#include "EigenSupport.h"

#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/spin_mutex.h>

#include <numeric>

using namespace pgo::PredefinedPotentialEnergies;
namespace ES = pgo::EigenSupport;

namespace pgo
{
namespace PredefinedPotentialEnergies
{
struct SurfaceTriangleDeformationBuf
{
  SurfaceTriangleDeformationBuf(int numV):
    locks(numV) {}
  std::vector<tbb::spin_mutex> locks;
};
}  // namespace PredefinedPotentialEnergies
}  // namespace pgo

SurfaceTriangleDeformation::SurfaceTriangleDeformation(const EigenSupport::VXd &restp, const Mesh::TriMeshGeo &meshIn, double eps):
  restPositions(restp), mesh(meshIn)
{
  elementWeights.resize(mesh.numTriangles());
  elementAbarInv.resize(mesh.numTriangles());
  updateRestInfo();

  std::vector<ES::TripletD> entries;
  for (int ei = 0; ei < mesh.numTriangles(); ei++) {
    for (int vi = 0; vi < 3; vi++) {
      int vidx_i = mesh.triVtxID(ei, vi);

      for (int vj = 0; vj < 3; vj++) {
        int vidx_j = mesh.triVtxID(ei, vj);

        for (int r = 0; r < 3; r++) {
          for (int c = 0; c < 3; c++) {
            entries.emplace_back(vidx_i * 3 + r, vidx_j * 3 + c, 1.0);
          }
        }
      }
    }
  }

  hessTemplate.resize(restPositions.size(), restPositions.size());
  hessTemplate.setFromTriplets(entries.begin(), entries.end());

  allDOFs.resize(restPositions.size());
  std::iota(allDOFs.begin(), allDOFs.end(), 0);

  elementOffsets.resize(mesh.numTriangles());
  for (int ei = 0; ei < mesh.numTriangles(); ei++) {
    for (int vi = 0; vi < 3; vi++) {
      int vidx_i = mesh.triVtxID(ei, vi);

      for (int vj = 0; vj < 3; vj++) {
        int vidx_j = mesh.triVtxID(ei, vj);

        for (int r = 0; r < 3; r++) {
          for (int c = 0; c < 3; c++) {
            elementOffsets[ei](vi * 3 + r, vj * 3 + c) = ES::findEntryOffset(hessTemplate, vidx_i * 3 + r, vidx_j * 3 + c);
            PGO_ALOG(elementOffsets[ei](vi * 3 + r, vj * 3 + c) >= 0);
          }
        }
      }
    }
  }

  // F = [x1 - x0, x2 - x0]
  dFdx.setZero();
  dFdx.block<3, 3>(0, 0) = ES::M3d::Identity() * -1.0;
  dFdx.block<3, 3>(3, 0) = ES::M3d::Identity() * -1.0;

  dFdx.block<3, 3>(0, 3) = ES::M3d::Identity();
  dFdx.block<3, 3>(3, 6) = ES::M3d::Identity();

  buf = std::make_shared<SurfaceTriangleDeformationBuf>(mesh.numVertices());
}

SurfaceTriangleDeformation::~SurfaceTriangleDeformation()
{
}

void SurfaceTriangleDeformation::setDOFs(const std::vector<int> &dofs)
{
  PGO_ALOG((int)dofs.size() == (int)allDOFs.size());
  allDOFs = dofs;
}

void SurfaceTriangleDeformation::updateRestInfo()
{
  tbb::parallel_for(0, mesh.numTriangles(), [&](int ei) {
    ES::V3d p[3] = {
      restPositions.segment<3>(mesh.triVtxID(ei, 0) * 3),
      restPositions.segment<3>(mesh.triVtxID(ei, 1) * 3),
      restPositions.segment<3>(mesh.triVtxID(ei, 2) * 3),
    };

    elementWeights[ei] = Mesh::getTriangleArea(p[0], p[1], p[2]);

    ES::M3x2d F;
    F.col(0) = restPositions.segment<3>(mesh.triVtxID(ei, 1) * 3) - restPositions.segment<3>(mesh.triVtxID(ei, 0) * 3);
    F.col(1) = restPositions.segment<3>(mesh.triVtxID(ei, 2) * 3) - restPositions.segment<3>(mesh.triVtxID(ei, 0) * 3);

    // elementAbarInv[ei].setIdentity();
    elementAbarInv[ei] = (F.transpose() * F).fullPivHouseholderQr().inverse();

    Eigen::SelfAdjointEigenSolver<ES::M2d> eigSolver(elementAbarInv[ei], Eigen::ComputeEigenvectors);
    ES::M2d U = eigSolver.eigenvectors();
    ES::V2d S = eigSolver.eigenvalues().array().sqrt();

    elementAbarInv[ei] = U * S.asDiagonal() * U.transpose();
  });
}

// F = [x1 - x0, x2 - x0]
// A = FTF
// S = eig(A)
// energy = 1/2 (|| S(0) - S(1) ||^2 + ||S - I||^2)
double SurfaceTriangleDeformation::func(EigenSupport::ConstRefVecXd x) const
{
  double energyAll = tbb::parallel_reduce(
    tbb::blocked_range<int>(0, mesh.numTriangles()), 0.0,
    [&](const tbb::blocked_range<int> &r, double init) -> double {
      for (int ei = r.begin(); ei != r.end(); ++ei) {
        ES::M3x2d F;
        F.col(0) = x.segment<3>(mesh.triVtxID(ei, 1) * 3) - x.segment<3>(mesh.triVtxID(ei, 0) * 3);
        F.col(1) = x.segment<3>(mesh.triVtxID(ei, 2) * 3) - x.segment<3>(mesh.triVtxID(ei, 0) * 3);

        ES::M2d A = elementAbarInv[ei] * (F.transpose() * F) * elementAbarInv[ei];
        // Eigen::SelfAdjointEigenSolver<ES::M2d> eigSolver(A);
        // ES::V2d S = eigSolver.eigenvalues();

        // double energy = (S(0) - S(1)) * (S(0) - S(1)) * coeffs[0] + (S - ES::V2d::Ones()).squaredNorm() * coeffs[1];
        double energy = A.trace() * A.trace() - 4 * A.determinant();

        // double v1 = (S(0) - S(1)) * (S(0) - S(1));
        // double v2 = A.trace() * A.trace() - 4 * A.determinant();
        // std::cout << v1 << ',' << v2 << std::endl;

        init += energy * elementWeights[ei];
      }

      return init;
    },
    [](double x, double y) -> double {
      return x + y;
    });

  return energyAll * 0.5;
}

// F = [x1 - x0, x2 - x0]
// A = M^-1 FTF
// S = eig(A)
// E = 1/2 (|| S(0) - S(1) ||^2 + ||S - I||^2)
// E = 1/2 ( (trA)^2 - 4 detA ) + 1/2 ((trA)^2 - 2 detA - 2 trA  + 2)
// dE/dA = trA * dtrA/dA - 2 d detA / dA
// dA/dF =
// dA/dFi =
void SurfaceTriangleDeformation::gradient(EigenSupport::ConstRefVecXd x, EigenSupport::RefVecXd grad) const
{
  grad.setZero();

  // for (int ei = 0; ei < mesh.numTriangles(); ei++) {
  tbb::parallel_for(0, mesh.numTriangles(), [&](int ei) {
    ES::M3x2d F;
    F.col(0) = x.segment<3>(mesh.triVtxID(ei, 1) * 3) - x.segment<3>(mesh.triVtxID(ei, 0) * 3);
    F.col(1) = x.segment<3>(mesh.triVtxID(ei, 2) * 3) - x.segment<3>(mesh.triVtxID(ei, 0) * 3);

    ES::M2d A = elementAbarInv[ei] * (F.transpose() * F) * elementAbarInv[ei];
    double trA = A.trace(), detA = A.determinant();

    ES::M2d dtrA_dA = ES::M2d::Identity();
    ES::M2d ddetA_dA;
    pgo::NonlinearOptimization::Determinant::Dim2::ddetA_dA(A.data(), ddetA_dA.data());

    ES::M2d dEdA = trA * dtrA_dA - 2 * ddetA_dA;
    // ES::M2d dEdA = ddetA_dA;

    Eigen::Matrix<double, 4, 6> dAdF;

    for (int c = 0; c < 2; c++) {
      for (int r = 0; r < 3; r++) {
        ES::M3x2d dFi;
        dFi.setZero();
        dFi(r, c) = 1.0;

        ES::M2d dA = elementAbarInv[ei] * (dFi.transpose() * F + F.transpose() * dFi) * elementAbarInv[ei];
        dAdF.col(c * 3 + r) = ES::Mp<ES::V4d>(dA.data());
      }
    }

    ES::V9d dEdx = dFdx.transpose() * dAdF.transpose() * ES::Mp<ES::V4d>(dEdA.data());
    dEdx *= elementWeights[ei];

    for (int i = 0; i < 3; i++) {
      buf->locks[mesh.triVtxID(ei, i)].lock();

      grad.segment<3>(mesh.triVtxID(ei, i) * 3) += dEdx.segment<3>(i * 3) * coeffs[0];

      buf->locks[mesh.triVtxID(ei, i)].unlock();
    }
  });
}

// F = [x1 - x0, x2 - x0]
// A = FTF
// S = eig(A)
// E = 1/2 (|| S(0) - S(1) ||^2 + ||S - I||^2)
// E = 1/2 ( (trA)^2 - 4 detA ) + 1/2 ((trA)^2 - 2 detA - 2 trA  + 2)
// dE/dA = trA * dtrA/dA - 2 d detA / dA
void SurfaceTriangleDeformation::hessian(EigenSupport::ConstRefVecXd x, EigenSupport::SpMatD &hess) const
{
  memset(hess.valuePtr(), 0, sizeof(double) * hess.nonZeros());

  // for (int ei = 0; ei < mesh.numTriangles(); ei++) {
  tbb::parallel_for(0, mesh.numTriangles(), [&](int ei) {
    ES::M3x2d F;
    F.col(0) = x.segment<3>(mesh.triVtxID(ei, 1) * 3) - x.segment<3>(mesh.triVtxID(ei, 0) * 3);
    F.col(1) = x.segment<3>(mesh.triVtxID(ei, 2) * 3) - x.segment<3>(mesh.triVtxID(ei, 0) * 3);

    ES::M2d A = elementAbarInv[ei] * (F.transpose() * F) * elementAbarInv[ei];
    double trA = A.trace(), detA = A.determinant();

    ES::M2d dtrA_dA = ES::M2d::Identity();
    ES::M2d ddetA_dA;
    pgo::NonlinearOptimization::Determinant::Dim2::ddetA_dA(A.data(), ddetA_dA.data());
    ES::M2d dEdA = trA * dtrA_dA - 2 * ddetA_dA;
    // ES::M2d dEdA = ddetA_dA;

    // dE/dA = trA * dtrA/dA - 2 ddetA/dA
    // d2E/dA2 = dtrA/dA * dtrA/dA - 2 d2detA/dA2
    ES::M4d d2detA_dA2;
    pgo::NonlinearOptimization::Determinant::Dim2::d2detA_dA2(A.data(), d2detA_dA2.data());
    ES::M4d d2EdA2 = ES::Mp<ES::V4d>(dtrA_dA.data()) * ES::Mp<ES::V4d>(dtrA_dA.data()).transpose() - 2 * d2detA_dA2;
    // ES::M4d d2EdA2 = d2detA_dA2;

    // dE/dx = dE/dA dA/dx
    // d2E/dx2 = dA/dx d2E/dA2 dA/dx + dE/dA d2A/dx2
    Eigen::Matrix<double, 4, 6> dAdF;
    for (int c = 0; c < 2; c++) {
      for (int r = 0; r < 3; r++) {
        ES::M3x2d dFi;
        dFi.setZero();
        dFi(r, c) = 1.0;

        ES::M2d dA = elementAbarInv[ei] * (dFi.transpose() * F + F.transpose() * dFi) * elementAbarInv[ei];
        dAdF.col(c * 3 + r) = ES::Mp<ES::V4d>(dA.data());
      }
    }

    ES::M6d d2AidF2[4];
    for (int c0 = 0; c0 < 2; c0++) {
      for (int r0 = 0; r0 < 3; r0++) {
        ES::M3x2d dFi;
        dFi.setZero();
        dFi(r0, c0) = 1.0;

        for (int c1 = 0; c1 < 2; c1++) {
          for (int r1 = 0; r1 < 3; r1++) {
            ES::M3x2d dFj;
            dFj.setZero();
            dFj(r1, c1) = 1.0;

            ES::M2d d2AdFidFj = elementAbarInv[ei] * (dFi.transpose() * dFj + dFj.transpose() * dFi) * elementAbarInv[ei];

            for (int k = 0; k < 4; k++) {
              d2AidF2[k](c0 * 3 + r0, c1 * 3 + r1) = d2AdFidFj.data()[k];
            }
          }
        }
      }
    }

    Eigen::Matrix<double, 4, 9> dAdx = dAdF * dFdx;
    Eigen::Matrix<double, 9, 9> d2Aidx2[4];
    for (int k = 0; k < 4; k++) {
      d2Aidx2[k] = dFdx.transpose() * d2AidF2[k] * dFdx;
    }

    // d2E/dx2 = dA/dx d2E/dA2 dA/dx + dE/dA d2A/dx2
    ES::M9d hessLocal;
    hessLocal = dAdx.transpose() * d2EdA2 * dAdx;
    for (int k = 0; k < 4; k++) {
      hessLocal += dEdA.data()[k] * d2Aidx2[k];
    }
    hessLocal *= elementWeights[ei];

    for (int vi = 0; vi < 3; vi++) {
      int vidx_i = mesh.triVtxID(ei, vi);

      for (int vj = 0; vj < 3; vj++) {
        int vidx_j = mesh.triVtxID(ei, vj);

        buf->locks[vidx_i].lock();

        for (int r = 0; r < 3; r++) {
          for (int c = 0; c < 3; c++) {
            std::ptrdiff_t offset = elementOffsets[ei](vi * 3 + r, vj * 3 + c);
            hess.valuePtr()[offset] += hessLocal(vi * 3 + r, vj * 3 + c) * coeffs[0];
          }
        }

        buf->locks[vidx_i].unlock();
      }
    }
  });
}