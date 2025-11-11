/*
author: Bohan Wang
copyright to USC, MIT
*/

#include "multiVertexPullingSoftConstraints.h"

#include "pgoLogging.h"

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
    coeffs.assign(bcCoeff, bcCoeff + numPts);
  }
  else {
    coeffs.assign(vertexIndices.size(), 1.0);
  }
}

void MultipleVertexPulling::setCoeff(const double *v)
{
  coeffs.assign(v, v + vertexIndices.size());
}

void MultipleVertexPulling::setTargetPos(const double *tgt)
{
  tgtp = ES::Mp<const ES::VXd>(tgt, vertexIndices.size() * 3);
}

double MultipleVertexPulling::func(ES::ConstRefVecXd u) const
{
  double eng = 0;
  for (size_t i = 0; i < vertexIndices.size(); i++) {
    ES::V3d p;
    if (isDisp) {
      p = u.segment<3>(vertexIndices[i] * 3) + restpAll.segment<3>(vertexIndices[i] * 3);
    }
    else {
      p = u.segment<3>(vertexIndices[i] * 3);
    }

    ES::V3d diff = p - tgtp.segment<3>(i * 3);
    eng += diff.dot(diff) * 0.5 * coeffs[i];
  }

  return eng * coeffAll;
}

void MultipleVertexPulling::gradient(ES::ConstRefVecXd u, ES::RefVecXd grad) const
{
  memset(grad.data(), 0, sizeof(double) * grad.size());

  for (size_t i = 0; i < vertexIndices.size(); i++) {
    int vtx = vertexIndices[i];
    ES::V3d p;
    if (isDisp) {
      p = u.segment<3>(vertexIndices[i] * 3) + restpAll.segment<3>(vertexIndices[i] * 3);
    }
    else {
      p = u.segment<3>(vertexIndices[i] * 3);
    }

    ES::V3d diff = p - tgtp.segment<3>(i * 3);

    grad.segment<3>(vtx * 3) = diff;
    grad.segment<3>(vtx * 3) *= coeffs[i];
  }

  grad *= coeffAll;
}

void MultipleVertexPulling::hessian(ES::ConstRefVecXd, ES::SpMatD &hess) const
{
  memset(hess.valuePtr(), 0, sizeof(double) * hess.nonZeros());

  for (size_t vi = 0; vi < vertexIndices.size(); vi++) {
    for (int i = 0; i < 3; i++) {
      hess.valuePtr()[KIndices[vi](i, i)] = coeffs[vi] * coeffAll;
    }
  }
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
