/*
author: Bohan Wang
copyright to USC, MIT
*/

#include "multiVertexPullingSoftConstraintsPOrder.h"
#include "pgoLogging.h"

#include <tbb/parallel_for.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/spin_mutex.h>

#include <numeric>
#include <iostream>

namespace pgo::ConstraintPotentialEnergies
{
struct MultipleVertexPullingSoftConstraintsPOrderBuf
{
  MultipleVertexPullingSoftConstraintsPOrderBuf(int numLocks):
    energyTLS(0.0), locks(numLocks) {}

  tbb::enumerable_thread_specific<double> energyTLS;
  std::vector<tbb::spin_mutex> locks;
};
}  // namespace pgo::ConstraintPotentialEnergies

using namespace pgo::ConstraintPotentialEnergies;

namespace ES = pgo::EigenSupport;

MultipleVertexPullingSoftConstraintsPOrder::MultipleVertexPullingSoftConstraintsPOrder(const ES::SpMatD &Koff, const EigenSupport::VXd &restPositionsAll,
  int numPoints, const double *constraintCoeffs, const double *constraintTargetPositions, const double *constraintNormals, const int *idx,
  double bcCoeff, int useN, int isd, double pval):
  PotentialEnergyAligningMeshConnectivity(Koff),
  restpAll(restPositionsAll), coeffAll(bcCoeff), useNormal(useN), isDisp(isd),
  pvalue(pval)
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

  buf = std::make_shared<MultipleVertexPullingSoftConstraintsPOrderBuf>(int(Koff.rows()));

  if (std::abs(pvalue - 2) < 1e-10) {
    isL2 = true;
  }
  else {
    isL2 = false;
  }
}

double MultipleVertexPullingSoftConstraintsPOrder::func(ES::ConstRefVecXd u) const
{
  for (auto &val : buf->energyTLS) {
    val = 0.0;
  }

  for (int ci = 0; ci < (int)coeffs.size(); ci++) {
    // tbb::parallel_for(    0, (int)coeffs.size(), [&](int ci) {
    if (std::abs(coeffs[ci]) < 1e-10)
      // return;
      continue;

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

    double proj = 0;
    if (useNormal) {
      proj = (p - p0).dot(n);
    }
    else {
      proj = (p - p0).norm();
    }

    if (isL2) {
      energyLocal += proj * proj * 0.5 * coeffs[ci];
    }
    else {
      energyLocal += std::pow(proj, pvalue) * 0.5 * coeffs[ci];
    }
  }  // ,    tbb::static_partitioner());

  double energyAll = std::accumulate(buf->energyTLS.begin(), buf->energyTLS.end(), 0.0) * coeffAll;

  return energyAll;
}

void MultipleVertexPullingSoftConstraintsPOrder::gradient(ES::ConstRefVecXd u, ES::RefVecXd grad) const
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
      ES::V3d gradLocal;
      if (isL2) {
        if (useNormal) {
          ES::M3d nnT = n * n.transpose();
          gradLocal = nnT * diff * coeffs[ci];
        }
        else {
          gradLocal = diff * coeffs[ci];
        }
      }
      else {
        double proj = 0;
        if (useNormal) {
          proj = (p - p0).dot(n);
        }
        else {
          proj = (p - p0).norm();
        }

        double halfPower = pvalue * 0.5;
        double proj2 = proj * proj;

        if (halfPower - 1 < 0) {
          if (proj2 < 1e-20) {
            proj2 = 1e-20;
          }
        }

        // ||p - p0||^r
        // d ((p - p0)^2)^(r/2)
        // = (r/2) * ((p - p0)^2)^(r / 2 - 1) * d (p - p0)^2 / dp
        // = (r/2) * ((p - p0)^2)^(r / 2 - 1) * 2(p - p0)
        gradLocal = halfPower * std::pow(proj2, halfPower - 1) * diff * coeffs[ci];
      }

      buf->locks[vid].lock();

      grad.segment<3>(vid * 3) += gradLocal;

      buf->locks[vid].unlock();
    },
    tbb::static_partitioner());

  grad *= coeffAll;
}

void MultipleVertexPullingSoftConstraintsPOrder::hessian(ES::ConstRefVecXd u, ES::SpMatD &hess) const
{
  memset(hess.valuePtr(), 0, sizeof(double) * hess.nonZeros());

  tbb::parallel_for(
    0, (int)coeffs.size(), [&](int ci) {
      if (std::abs(coeffs[ci]) < 1e-9)
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
      ES::M3d hessianLocal;
      if (isL2) {
        if (useNormal == 0) {
          hessianLocal = ES::M3d::Identity() * coeffs[ci];
        }
        else {
          hessianLocal = ES::tensorProduct(n, n) * coeffs[ci];
        }
      }
      else {
        double proj = 0;
        if (useNormal) {
          proj = (p - p0).dot(n);
        }
        else {
          proj = (p - p0).norm();
        }

        double halfPower = pvalue * 0.5;
        double proj2 = proj * proj;

        if (halfPower - 1 < 0) {
          if (proj2 < 1e-20) {
            proj2 = 1e-20;
          }
        }

        // ||p - p0||^r
        // d ((p - p0)^2)^(r/2)
        // = (r/2) * ((p - p0)^2)^(r / 2 - 1) * d (p - p0)^2 / dp
        // = (r/2) * ((p - p0)^2)^(r / 2 - 1) * 2(p - p0)
        // = (r/2) * z^(r/2 - 1) * dz/dp
        // dd =
        // = (r/2) * (r/2 - 1) * z^(r/2 - 2) * dz/dp * dzdp + (r/2) * z^(r/2 - 1) * d2z/dp2
        hessianLocal = halfPower * (halfPower - 1) * std::pow(proj2, halfPower - 2) * diff * diff.transpose() * 2;
        hessianLocal += halfPower * std::pow(proj2, halfPower - 1) * ES::M3d::Identity();
        hessianLocal *= coeffs[ci];
      }

      for (int dofi = 0; dofi < 3; dofi++) {
        for (int dofj = 0; dofj < 3; dofj++) {
          int localRow = dofi;
          int localCol = dofj;
          int offset = hessianIndices[ci](localRow, localCol);

          buf->locks[vid].lock();

          hess.valuePtr()[offset] += hessianLocal(dofi, dofj);

          buf->locks[vid].unlock();
        }
      }
    },
    tbb::static_partitioner());

  // cblas_dscal((int)hess.nonZeros(), coeffAll, hess.valuePtr(), 1);
  (ES::Mp<ES::VXd>(hess.valuePtr(), hess.nonZeros())) *= coeffAll;
}

void MultipleVertexPullingSoftConstraintsPOrder::setCoeff(const double *v)
{
  memcpy(coeffs.data(), v, coeffs.size() * sizeof(double));
}

void MultipleVertexPullingSoftConstraintsPOrder::setTargetPos(const double *tgt)
{
  memcpy(tgtp.data(), tgt, tgtp.size() * sizeof(double));
}

void MultipleVertexPullingSoftConstraintsPOrder::setNormals(const double *n)
{
  memcpy(normals.data(), n, normals.size() * sizeof(double));
}

void MultipleVertexPullingSoftConstraintsPOrder::printErrorInfo(ES::ConstRefVecXd u) const
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
