/*
author: Bohan Wang
copyright to USC, MIT
*/

#include "multiVertexPullingSoftConstraints.h"

#include "pgoLogging.h"

#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>

#include <iostream>

using namespace pgo::ConstraintPotentialEnergies;

namespace ES = pgo::EigenSupport;

MultipleVertexPulling::MultipleVertexPulling(const EigenSupport::SpMatD &Koff, const double *restPositionsAll,
  int numPts, const int *vi, const double *tgt, const double *bcCoeff, int isd):
  PotentialEnergyAligningMeshConnectivity(Koff),
  isDisp(isd)
{
  vertexIndices.assign(vi, vi + numPts);
  KIndices.assign(numPts, M3i::Constant(-1));

  for (size_t vi = 0; vi < vertexIndices.size(); vi++) {
    int vid = vertexIndices[vi];

    KIndices[vi] = Eigen::Matrix<ES::IDX, 3, 3>::Constant(-1);
    for (int i = 0; i < 3; i++) {
      KIndices[vi](i, i) = ES::findEntryOffset(Koff, vid * 3 + i, vid * 3 + i);
      PGO_ALOG(KIndices[vi](i, i) >= 0);
    }
  }

  tgtp = ES::Mp<const ES::VXd>(tgt, vertexIndices.size() * 3);
  restpAll = ES::Mp<const ES::VXd>(restPositionsAll, Koff.rows());

  if (bcCoeff) {
    coeffs = ES::Mp<const ES::VXd>(bcCoeff, numPts);
  }
  else {
    coeffs.setConstant(vertexIndices.size(), 1.0);
  }

  masks.setOnes(vertexIndices.size() * 3);
}

void MultipleVertexPulling::setCoeff(const double *v)
{
  coeffs = ES::Mp<const ES::VXd>(v, vertexIndices.size());
}

void MultipleVertexPulling::setTargetPos(const double *tgt)
{
  tgtp = ES::Mp<const ES::VXd>(tgt, vertexIndices.size() * 3);
}

void MultipleVertexPulling::setMasks(const double *v)
{
  masks = ES::Mp<const ES::VXd>(v, vertexIndices.size() * 3);
}

double MultipleVertexPulling::func(ES::ConstRefVecXd u) const
{
  // double eng = 0;
  // for (size_t i = 0; i < vertexIndices.size(); i++) {
  auto computeEnergyForVertex = [&](size_t i) -> double {
    ES::V3d p;
    if (isDisp) {
      p = u.segment<3>(vertexIndices[i] * 3) + restpAll.segment<3>(vertexIndices[i] * 3);
    }
    else {
      p = u.segment<3>(vertexIndices[i] * 3);
    }

    ES::V3d diff = p - tgtp.segment<3>(i * 3);
    diff = diff.cwiseProduct(masks.segment<3>(i * 3));

    return diff.dot(diff) * 0.5 * coeffs[i];
  };

  double eng = tbb::parallel_reduce(tbb::blocked_range<size_t>(0, vertexIndices.size()), 0.0,  //
    [&](const tbb::blocked_range<size_t> &r, double init) -> double {
      for (size_t i = r.begin(); i != r.end(); ++i) {
        init += computeEnergyForVertex(i);
      }
      return init; }, std::plus<double>());

  return eng * coeffAll;
}

void MultipleVertexPulling::gradient(ES::ConstRefVecXd u, ES::RefVecXd grad) const
{
  grad.setZero();

  // for (size_t i = 0; i < vertexIndices.size(); i++) {
  auto computeGradientForVertex = [&](size_t i) {
    int vtx = vertexIndices[i];

    ES::V3d p;
    if (isDisp) {
      p = u.segment<3>(vertexIndices[i] * 3) + restpAll.segment<3>(vertexIndices[i] * 3);
    }
    else {
      p = u.segment<3>(vertexIndices[i] * 3);
    }

    ES::V3d diff = p - tgtp.segment<3>(i * 3);
    diff = diff.cwiseProduct(masks.segment<3>(i * 3)).cwiseProduct(masks.segment<3>(i * 3));

    // E = 1/2  (W(p - pbar))^2
    // dE/dp = W^2 (p - pbar)

    grad.segment<3>(vtx * 3) = diff;
    grad.segment<3>(vtx * 3) *= coeffs[i];
  };

  tbb::parallel_for(tbb::blocked_range<size_t>(0, vertexIndices.size()),
    [&](const tbb::blocked_range<size_t> &r) {
      for (size_t i = r.begin(); i != r.end(); ++i) {
        computeGradientForVertex(i);
      }
    });

  grad *= coeffAll;
}

void MultipleVertexPulling::hessian(ES::ConstRefVecXd, ES::SpMatD &hess) const
{
  memset(hess.valuePtr(), 0, sizeof(double) * hess.nonZeros());

  // for (size_t vi = 0; vi < vertexIndices.size(); vi++) {
  auto computeHessianForVertex = [&](size_t vi) {
    ES::V3d w = masks.segment<3>(vi * 3).cwiseProduct(masks.segment<3>(vi * 3));

    for (int i = 0; i < 3; i++) {
      hess.valuePtr()[KIndices[vi](i, i)] = w[i] * coeffs[vi] * coeffAll;
    }
  };

  tbb::parallel_for(tbb::blocked_range<size_t>(0, vertexIndices.size()),
    [&](const tbb::blocked_range<size_t> &r) {
      for (size_t vi = r.begin(); vi != r.end(); ++vi) {
        computeHessianForVertex(vi);
      }
    });
}

void MultipleVertexPulling::printErrorInfo(ES::ConstRefVecXd u) const
{
  ES::VXd diff(vertexIndices.size() * 3);
  ES::VXd diff1(vertexIndices.size() * 3);
  for (size_t i = 0; i < vertexIndices.size(); i++) {
    ES::V3d p;
    if (isDisp) {
      p = u.segment<3>(vertexIndices[i] * 3) + restpAll.segment<3>(vertexIndices[i] * 3);
    }
    else {
      p = u.segment<3>(vertexIndices[i] * 3);
    }

    diff.segment<3>(i * 3) = p - tgtp.segment<3>(i * 3);
    diff.segment<3>(i * 3) = diff.segment<3>(i * 3).cwiseProduct(masks.segment<3>(i * 3));
    diff1.segment<3>(i * 3) = diff.segment<3>(i * 3) * coeffs[i];
  }

  std::cout << "  ||Wu - bcu||=" << diff.squaredNorm() << std::endl;
  std::cout << "  ||Wu - bcu||_Z=" << diff1.dot(diff) << std::endl;
  std::cout << "  func(u)=" << func(u) << std::endl;

  double mind = 1e100, maxd = 0, avgd = 0;
  for (ES::IDX i = 0; i < diff.size() / 3; i++) {
    if (coeffs[i] < 1e-9)
      continue;

    double d = diff.segment<3>(i * 3).norm();
    mind = std::min(mind, d);
    maxd = std::max(maxd, d);
    avgd += d;
  }
  avgd /= (double)(diff.size() / 3);
  std::cout << "distance error info: " << mind << '/' << avgd << '/' << maxd << std::endl;
}
