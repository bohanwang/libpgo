/*
author: Bohan Wang
copyright to USC, MIT
*/

#include "multiVertexSlidingSoftConstraints.h"
#include "pgoLogging.h"

#include <tbb/parallel_for.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/spin_mutex.h>

#include <numeric>
#include <iostream>

namespace pgo::ConstraintPotentialEnergies
{
struct MultipleVertexSlidingBuf
{
  MultipleVertexSlidingBuf():
    energyTLS(0.0) {}

  tbb::enumerable_thread_specific<double> energyTLS;
  std::vector<tbb::spin_mutex> locks;
};
}  // namespace pgo::ConstraintPotentialEnergies

using namespace pgo::ConstraintPotentialEnergies;

namespace ES = pgo::EigenSupport;

MultipleVertexSliding::MultipleVertexSliding(const ES::SpMatD &Koff, const EigenSupport::VXd &restPositionsAll,
  int numPoints, const double *constraintCoeffs, const double *constraintTargetPositions, const double *constraintNormals, const int *idx,
  double bcCoeff, int useN, int cp, int isd):
  PotentialEnergyAligningMeshConnectivity(Koff),
  restpAll(restPositionsAll),
  coeffAll(bcCoeff), useNormal(useN), checkPenetration(cp), isDisp(isd)
{
  tgtp = ES::Mp<const ES::VXd>(constraintTargetPositions, numPoints * 3);
  normals = ES::Mp<const ES::VXd>(constraintNormals, numPoints * 3);
  coeffs = ES::Mp<const ES::VXd>(constraintCoeffs, numPoints);

  vertexIndices.assign(idx, idx + numPoints);

  hessianIndices.resize(numPoints);
  for (int pi = 0; pi < numPoints; pi++) {
    hessianIndices[pi].setConstant(-1);

    ES::M3i &idx = hessianIndices[pi];
    int vid = vertexIndices[pi];

    for (int dofi = 0; dofi < 3; dofi++) {
      for (int dofj = 0; dofj < 3; dofj++) {
        int localRow = dofi;
        int localCol = dofj;

        int globalRow = vid * 3 + dofi;
        int globalCol = vid * 3 + dofj;

        idx(localRow, localCol) = (int)ES::findEntryOffset(Koff, globalRow, globalCol);
        PGO_ALOG(idx(localRow, localCol) >= 0);
      }
    }
  }
  hessianConstant = Koff;

  buf = std::make_shared<MultipleVertexSlidingBuf>();
  buf->locks = std::vector<tbb::spin_mutex>(std::max(Koff.nonZeros(), Koff.rows()));

  if (checkPenetration == 0)
    computeHessian();
}

double MultipleVertexSliding::func(ES::ConstRefVecXd u) const
{
  for (auto &val : buf->energyTLS) {
    val = 0.0;
  }

  tbb::parallel_for(
    0, (int)coeffs.size(), [&](int ci) {
  // for (int ci = 0; ci < (int)coeffs.size(); ci++) {
    if (std::abs(coeffs[ci]) < 1e-10)
      return;
      // continue;

    int vid = vertexIndices[ci];
    ES::V3d n = normals.segment<3>(ci * 3);
    ES::V3d p0 = tgtp.segment<3>(ci * 3);
    ES::V3d p;
    if (isDisp) {
      p = u.segment<3>(vid * 3) + restpAll.segment<3>(vid * 3);
    }
    else {
      p = u.segment<3>(vid * 3);
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
  }  ,    tbb::static_partitioner());

  double energyAll = std::accumulate(buf->energyTLS.begin(), buf->energyTLS.end(), 0.0) * coeffAll;

  return energyAll;
}

void MultipleVertexSliding::gradient(ES::ConstRefVecXd u, ES::RefVecXd grad) const
{
  grad.setZero();

  tbb::parallel_for(
    0, (int)coeffs.size(), [&](int ci) {
      if (std::abs(coeffs[ci]) < 1e-10)
        return;

      int vid = vertexIndices[ci];
      ES::V3d n = normals.segment<3>(ci * 3);
      ES::V3d p0 = tgtp.segment<3>(ci * 3);
      ES::V3d p;
      if (isDisp) {
        p = u.segment<3>(vid * 3) + restpAll.segment<3>(vid * 3);
      }
      else {
        p = u.segment<3>(vid * 3);
      }
      ES::V3d diff = p - p0;

      if (useNormal) {
        double proj = diff.dot(n);

        if ((checkPenetration && proj < 0) || checkPenetration == 0) {
          ES::V3d gradLocal = n * proj * coeffs[ci] * coeffAll;

          buf->locks[vid].lock();

          grad.segment<3>(vid * 3) += gradLocal;

          buf->locks[vid].unlock();
        }
      }
      else {
        ES::V3d gradLocal = diff * coeffs[ci] * coeffAll;

        buf->locks[vid].lock();

        grad.segment<3>(vid * 3) += gradLocal;

        buf->locks[vid].unlock();
      }
    },
    tbb::static_partitioner());
}

void MultipleVertexSliding::computeHessian()
{
  memset(hessianConstant.valuePtr(), 0, sizeof(double) * hessianConstant.nonZeros());

  tbb::parallel_for(
    0, (int)coeffs.size(), [&](int ci) {
      ES::V3d n = normals.segment<3>(ci * 3);
      ES::M3d nnT = ES::tensorProduct(n, n) * coeffs[ci];

      if (useNormal == 0) {
        nnT = ES::M3d::Identity() * coeffs[ci];
      }

      for (int dofi = 0; dofi < 3; dofi++) {
        for (int dofj = 0; dofj < 3; dofj++) {
          int localRow = dofi;
          int localCol = dofj;
          int offset = hessianIndices[ci](localRow, localCol);

          buf->locks[offset].lock();

          hessianConstant.valuePtr()[offset] += nnT(dofi, dofj);

          buf->locks[offset].unlock();
        }
      }
    },
    tbb::static_partitioner());
}

void MultipleVertexSliding::hessian(ES::ConstRefVecXd u, ES::SpMatD &hess) const
{
  if (checkPenetration) {
    memset(hess.valuePtr(), 0, sizeof(double) * hess.nonZeros());

    tbb::parallel_for(
      0, (int)coeffs.size(), [&](int ci) {
        if (std::abs(coeffs[ci]) < 1e-10)
          return;

        int vid = vertexIndices[ci];
        ES::V3d n = normals.segment<3>(ci * 3);
        ES::V3d p0 = tgtp.segment<3>(ci * 3);
        ES::V3d p;
        if (isDisp) {
          p = u.segment<3>(vid * 3) + restpAll.segment<3>(vid * 3);
        }
        else {
          p = u.segment<3>(vid * 3);
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

        for (int dofi = 0; dofi < 3; dofi++) {
          for (int dofj = 0; dofj < 3; dofj++) {
            int localRow = dofi;
            int localCol = dofj;
            int offset = hessianIndices[ci](localRow, localCol);

            buf->locks[offset].lock();

            hess.valuePtr()[offset] += nnT(dofi, dofj);

            buf->locks[offset].unlock();
          }
        }
      },
      tbb::static_partitioner());
  }
  else {
    memcpy(hess.valuePtr(), hessianConstant.valuePtr(), sizeof(double) * hess.nonZeros());
  }

  (ES::Mp<ES::VXd>(hess.valuePtr(), hess.nonZeros())) *= coeffAll;
}

void MultipleVertexSliding::setCoeff(const double *v)
{
  memcpy(coeffs.data(), v, coeffs.size() * sizeof(double));
}

void MultipleVertexSliding::setTargetPos(const double *tgt)
{
  memcpy(tgtp.data(), tgt, tgtp.size() * sizeof(double));
}

void MultipleVertexSliding::setNormals(const double *n)
{
  memcpy(normals.data(), n, normals.size() * sizeof(double));
}

void MultipleVertexSliding::printErrorInfo(ES::ConstRefVecXd u) const
{
  ES::VXd dist(coeffs.size());
  ES::VXd dist2(coeffs.size());

  tbb::parallel_for(
    0, (int)coeffs.size(), [&](int ci) {
      int vid = vertexIndices[ci];
      ES::V3d n = normals.segment<3>(ci * 3);
      ES::V3d p0 = tgtp.segment<3>(ci * 3);
      ES::V3d p;
      if (isDisp) {
        p = u.segment<3>(vid * 3) + restpAll.segment<3>(vid * 3);
      }
      else {
        p = u.segment<3>(vid * 3);
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
