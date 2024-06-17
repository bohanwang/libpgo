/*
author: Bohan Wang
copyright to USC, MIT
*/

#include "barycentricCoordinateSlidingSoftConstraints.h"

#include "pgoLogging.h"
#include "EigenSupport.h"

#include <tbb/parallel_for.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/spin_mutex.h>

#include <numeric>
#include <iostream>

namespace pgo
{
namespace ConstraintPotentialEnergies
{
struct BarycentricCoordinateSlidingBuf
{
  BarycentricCoordinateSlidingBuf(int n):
    energyTLS(0.0), locks(n) {}

  tbb::enumerable_thread_specific<double> energyTLS;
  std::vector<tbb::spin_mutex> locks;
};
}  // namespace ConstraintPotentialEnergies
}  // namespace pgo

using namespace pgo::ConstraintPotentialEnergies;
namespace ES = pgo::EigenSupport;

BarycentricCoordinateSliding::BarycentricCoordinateSliding(const ES::SpMatD &Koff, const ES::VXd &restPositionsAll,
  int numPoints, const double *constraintCoeffs, const double *constraintTargetPositions, const double *constraintNormals,
  const int *idx, const double *barycentricWeights, int nEleVtx, double bcCoeff, int useN, int cp, int isDisp):
  PotentialEnergyAligningMeshConnectivity(Koff),
  restpAll(restPositionsAll), coeffAll(bcCoeff),
  numElementVertices(nEleVtx), useNormal(useN), checkPenetration(cp), isInputDisplacement(isDisp)
{
  tgtp = ES::Mp<const ES::VXd>(constraintTargetPositions, numPoints * 3);
  normals = ES::Mp<const ES::VXd>(constraintNormals, numPoints * 3);
  coeffs = ES::Mp<const ES::VXd>(constraintCoeffs, numPoints);

  vertexIndices.assign(idx, idx + numPoints * numElementVertices);
  vertexBarycentricWeights.assign(barycentricWeights, barycentricWeights + numPoints * numElementVertices);

  hessianIndices.resize(numPoints);
  for (int pi = 0; pi < numPoints; pi++) {
    hessianIndices[pi].setConstant(numElementVertices * 3, numElementVertices * 3, -1);
    ES::MXi &idx = hessianIndices[pi];
    for (int vi = 0; vi < numElementVertices; vi++) {
      int vidi = vertexIndices[pi * numElementVertices + vi];

      for (int vj = 0; vj < numElementVertices; vj++) {
        int vidj = vertexIndices[pi * numElementVertices + vj];

        for (int dofi = 0; dofi < 3; dofi++) {
          for (int dofj = 0; dofj < 3; dofj++) {
            int localRow = vi * 3 + dofi;
            int localCol = vj * 3 + dofj;

            int globalRow = vidi * 3 + dofi;
            int globalCol = vidj * 3 + dofj;

            idx(localRow, localCol) = (int)ES::findEntryOffset(Koff, globalRow, globalCol);
            PGO_ALOG(idx(localRow, localCol) >= 0);
          }
        }
      }
    }
  }

  hessianConstant = Koff;
  buf = std::make_shared<BarycentricCoordinateSlidingBuf>(Koff.rows());

  if (checkPenetration == 0)
    computeHessian();
}

// E = (n^T (sum w_i (restp + u) - tgtp))^2 * 0.5
double BarycentricCoordinateSliding::func(ES::ConstRefVecXd q) const
{
  for (auto &val : buf->energyTLS) {
    val = 0.0;
  }

  tbb::parallel_for(
    0, (int)coeffs.size(), [&](int ci) {
      if (std::abs(coeffs[ci]) < 1e-10)
        return;

      ES::V3d n = normals.segment<3>(ci * 3);
      ES::V3d p0 = tgtp.segment<3>(ci * 3);
      ES::V3d p;
      p.setZero();
      for (int vi = 0; vi < numElementVertices; vi++) {
        int vid = vertexIndices[ci * numElementVertices + vi];

        ES::V3d pp;
        if (isInputDisplacement) {
          pp = q.segment<3>(vid * 3) + restpAll.segment<3>(vid * 3);
        }
        else {
          pp = q.segment<3>(vid * 3);
        }

        p += pp * vertexBarycentricWeights[ci * numElementVertices + vi];
      }

      double &energyLocal = buf->energyTLS.local();

      if (useNormal) {
        double proj = (p - p0).dot(n);

        if ((checkPenetration && proj < 0) || checkPenetration == 0) {
          energyLocal += proj * proj * 0.5 * coeffs[ci];
        }
      }
      else {
        energyLocal += (p - p0).squaredNorm() * 0.5 * coeffs[ci];
      }
    },
    tbb::static_partitioner());

  double energyAll = std::accumulate(buf->energyTLS.begin(), buf->energyTLS.end(), 0.0) * coeffAll;

  return energyAll;
}

// dE = (n^T (sum w_i (restp + u) - tgtp)) d((n^T (sum w_i (restp + u) - tgtp)))
// dE = (sum w_i (restp + u) - tgtp) d(sum w_i (restp + u) - tgtp)
void BarycentricCoordinateSliding::gradient(ES::ConstRefVecXd q, ES::RefVecXd grad) const
{
  grad.setZero();

  tbb::parallel_for(
    0, (int)coeffs.size(), [&](int ci) {
      if (std::abs(coeffs[ci]) < 1e-10)
        return;

      ES::V3d n = normals.segment<3>(ci * 3);
      ES::V3d p0 = tgtp.segment<3>(ci * 3);
      ES::V3d p;
      p.setZero();
      for (int vi = 0; vi < numElementVertices; vi++) {
        int vid = vertexIndices[ci * numElementVertices + vi];
        ES::V3d pp;
        if (isInputDisplacement) {
          pp = q.segment<3>(vid * 3) + restpAll.segment<3>(vid * 3);
        }
        else {
          pp = q.segment<3>(vid * 3);
        }
        p += pp * vertexBarycentricWeights[ci * numElementVertices + vi];
      }

      ES::V3d diff = p - p0;

      if (useNormal) {
        double proj = diff.dot(n);

        if ((checkPenetration && proj < 0) || checkPenetration == 0) {
          for (int vi = 0; vi < numElementVertices; vi++) {
            ES::V3d gradLocal = n * vertexBarycentricWeights[ci * numElementVertices + vi] * proj * coeffs[ci] * coeffAll;
            int vid = vertexIndices[ci * numElementVertices + vi];

            buf->locks[vid].lock();

            grad.segment<3>(vid * 3) += gradLocal;

            buf->locks[vid].unlock();
          }
        }
      }
      else {
        for (int vi = 0; vi < numElementVertices; vi++) {
          ES::V3d gradLocal = diff * vertexBarycentricWeights[ci * numElementVertices + vi] * coeffs[ci] * coeffAll;
          int vid = vertexIndices[ci * numElementVertices + vi];

          buf->locks[vid].lock();

          grad.segment<3>(vid * 3) += gradLocal;

          buf->locks[vid].unlock();
        }
      }
    },
    tbb::static_partitioner());
}
// use n
// E = (n^T (sum w_i (restp + u) - tgtp)) ((n^T (sum w_i (restp + u) - tgtp)))
// W = [ w1I w2I w3I w4I]
// d2E = W^T (nnT) W

// not use n
// E = (sum w_i (restp + u) - tgtp) (sum w_i (restp + u) - tgtp)
// W = [ w1I w2I w3I w4I]
// d2E = W^T W

void BarycentricCoordinateSliding::computeHessian()
{
  memset(hessianConstant.valuePtr(), 0, sizeof(double) * hessianConstant.nonZeros());

  tbb::parallel_for(
    0, (int)coeffs.size(), [&](int ci) {
      if (std::abs(coeffs[ci]) < 1e-10)
        return;

      ES::V3d n = normals.segment<3>(ci * 3);
      ES::M3d nnT = ES::tensorProduct(n, n) * coeffs[ci];

      if (useNormal == 0) {
        nnT = ES::M3d::Identity() * coeffs[ci];
      }

      const double *weights = vertexBarycentricWeights.data() + (ci * numElementVertices);

      for (int vi = 0; vi < numElementVertices; vi++) {
        int vid = vertexIndices[ci * numElementVertices + vi];

        for (int vj = 0; vj < numElementVertices; vj++) {
          ES::M3d hLocal = nnT * weights[vi] * weights[vj];

          for (int dofi = 0; dofi < 3; dofi++) {
            buf->locks[vid * 3 + dofi].lock();

            for (int dofj = 0; dofj < 3; dofj++) {
              int localRow = vi * 3 + dofi;
              int localCol = vj * 3 + dofj;
              int offset = hessianIndices[ci](localRow, localCol);

              hessianConstant.valuePtr()[offset] += hLocal(dofi, dofj);
            }

            buf->locks[vid * 3 + dofi].unlock();
          }
        }
      }
    },
    tbb::static_partitioner());
}

void BarycentricCoordinateSliding::hessian(ES::ConstRefVecXd q, ES::SpMatD &hess) const
{
  if (checkPenetration) {
    memset(hess.valuePtr(), 0, sizeof(double) * hess.nonZeros());

    tbb::parallel_for(
      0, (int)coeffs.size(), [&](int ci) {
        if (std::abs(coeffs[ci]) < 1e-10)
          return;

        ES::V3d n = normals.segment<3>(ci * 3);
        ES::V3d p0 = tgtp.segment<3>(ci * 3);
        ES::V3d p;
        p.setZero();
        for (int vi = 0; vi < numElementVertices; vi++) {
          int vid = vertexIndices[ci * numElementVertices + vi];
          ES::V3d pp;
          if (isInputDisplacement) {
            pp = q.segment<3>(vid * 3) + restpAll.segment<3>(vid * 3);
          }
          else {
            pp = q.segment<3>(vid * 3);
          }
          p += pp * vertexBarycentricWeights[ci * numElementVertices + vi];
        }

        ES::V3d diff = p - p0;
        if (diff.dot(n) > 0)
          return;

        ES::M3d nnT;
        if (useNormal == 0) {
          nnT = ES::M3d::Identity() * coeffs[ci];
        }
        else {
          nnT = ES::tensorProduct(n, n) * coeffs[ci];
        }

        const double *weights = vertexBarycentricWeights.data() + (ci * numElementVertices);

        for (int vi = 0; vi < numElementVertices; vi++) {
          int vid = vertexIndices[ci * numElementVertices + vi];

          for (int vj = 0; vj < numElementVertices; vj++) {
            ES::M3d hLocal = nnT * weights[vi] * weights[vj];

            for (int dofi = 0; dofi < 3; dofi++) {
              buf->locks[vid * 3 + dofi].lock();

              for (int dofj = 0; dofj < 3; dofj++) {
                int localRow = vi * 3 + dofi;
                int localCol = vj * 3 + dofj;
                int offset = hessianIndices[ci](localRow, localCol);

                hess.valuePtr()[offset] += hLocal(dofi, dofj);
              }

              buf->locks[vid * 3 + dofi].unlock();
            }
          }
        }
      },
      tbb::static_partitioner());
  }
  else {
    memcpy(hess.valuePtr(), hessianConstant.valuePtr(), sizeof(double) * hess.nonZeros());
  }

  // cblas_dscal((int)hess.nonZeros(), coeffAll, hess.valuePtr(), 1);
  (ES::Mp<ES::VXd>(hess.valuePtr(), hess.nonZeros())) *= coeffAll;
}

void BarycentricCoordinateSliding::setCoeff(const double *v)
{
  memcpy(coeffs.data(), v, coeffs.size() * sizeof(double));
}

void BarycentricCoordinateSliding::setTargetPos(const double *tgt)
{
  memcpy(tgtp.data(), tgt, tgtp.size() * sizeof(double));
}

void BarycentricCoordinateSliding::setNormals(const double *n)
{
  memcpy(normals.data(), n, normals.size() * sizeof(double));
}

void BarycentricCoordinateSliding::printErrorInfo(ES::ConstRefVecXd q) const
{
  ES::VXd dist(coeffs.size());
  ES::VXd dist2(coeffs.size());

  tbb::parallel_for(
    0, (int)coeffs.size(), [&](int ci) {
      ES::V3d n = normals.segment<3>(ci * 3);
      ES::V3d p0 = tgtp.segment<3>(ci * 3);
      ES::V3d p;
      p.setZero();
      for (int vi = 0; vi < numElementVertices; vi++) {
        int vid = vertexIndices[ci * numElementVertices + vi];
        ES::V3d pp;
        if (isInputDisplacement) {
          pp = q.segment<3>(vid * 3) + restpAll.segment<3>(vid * 3);
        }
        else {
          pp = q.segment<3>(vid * 3);
        }
        p += pp * vertexBarycentricWeights[ci * numElementVertices + vi];
      }

      ES::V3d diff = p - p0;
      dist[ci] = std::abs(diff.dot(n));
      dist2[ci] = diff.norm();
    },
    tbb::static_partitioner());

  double mind = 1e100, maxd = 0, avgd = 0;
  for (ES::IDX i = 0; i < dist.size(); i++) {
    if (std::abs(coeffs[i]) < 1e-9)
      continue;

    double d = dist[i];
    mind = std::min(mind, d);
    maxd = std::max(maxd, d);
    avgd += d;
  }
  avgd /= (double)(dist.size());
  std::cout << "distance error info (proj): " << mind << '/' << avgd << '/' << maxd << std::endl;

  mind = 1e100, maxd = 0, avgd = 0;
  for (ES::IDX i = 0; i < dist.size(); i++) {
    if (std::abs(coeffs[i]) < 1e-9)
      continue;

    double d = dist2[i];
    mind = std::min(mind, d);
    maxd = std::max(maxd, d);
    avgd += d;
  }
  avgd /= (double)(dist2.size());
  std::cout << "distance error info: " << mind << '/' << avgd << '/' << maxd << std::endl;
}
