/*
author: Bohan Wang
copyright to USC, MIT
*/

#include "multiVertexConstrainedRigidMotion.h"

#include "polarDecompositionDerivatives.h"
#include "pgoLogging.h"

#include <tbb/parallel_for.h>

#include <numeric>

using namespace pgo::NonlinearOptimization;
using namespace pgo::ConstraintPotentialEnergies;

namespace ES = pgo::EigenSupport;

MultipleVertexConstrainedRigidMotion::MultipleVertexConstrainedRigidMotion(const ES::VXd &rp, const std::vector<int> &vi, const double *vw,
  double coeffR, double coefft, bool frm, bool id):
  restPositions(rp),
  vertexIndices(vi), flexibleRigidMotion(frm), isDisp(id)
{
  if (vw) {
    vertexWeights.assign(vw, vw + vertexIndices.size());
  }
  else {
    vertexWeights.assign(vertexIndices.size(), 1.0);
  }

  coeff[0] = coeffR;
  coeff[1] = coefft;

  for (int i = 0; i < (int)vertexIndices.size(); i++) {
    vertexIDtoIndices[vertexIndices[i]] = i;
  }

  // compute rest info
  double wAll = 0;
  centerRest.setZero();
  for (size_t i = 0; i < vertexIndices.size(); i++) {
    centerRest += restPositions.segment<3>(vertexIndices[i] * 3) * vertexWeights[i];
    wAll += vertexWeights[i];
  }
  centerRest /= wAll;

  for (size_t i = 0; i < vertexIndices.size(); i++) {
    ES::V3d qi = restPositions.segment<3>(vertexIndices[i] * 3) - centerRest;
    q.emplace_back(qi);
  }

  std::vector<ES::TripletD> entries;
  for (size_t i = 0; i < vertexIndices.size(); i++) {
    for (int dof = 0; dof < 3; dof++)
      entries.emplace_back(dof, vertexIndices[i] * 3 + dof, vertexWeights[i] / wAll);
  }
  Z.resize(3, restPositions.size());
  Z.setFromTriplets(entries.begin(), entries.end());

  for (size_t i = 0; i < vertexIndices.size(); i++) {
    entries.clear();
    for (int dof = 0; dof < 3; dof++)
      entries.emplace_back(dof, vertexIndices[i] * 3 + dof, 1.0);

    ES::SpMatD Si(3, restPositions.size());
    Si.setFromTriplets(entries.begin(), entries.end());

    Mi.emplace_back(Si - Z);
  }

  ES::mm(Z, Z, ZTZ, 1);

  entries.clear();
  for (int vi = 0; vi < (int)vertexIndices.size(); vi++) {
    for (int dofi = 0; dofi < 3; dofi++) {
      for (int vj = 0; vj < (int)vertexIndices.size(); vj++) {
        for (int dofj = 0; dofj < 3; dofj++) {
          int globalRow = vertexIndices[vi] * 3 + dofi;
          int globalCol = vertexIndices[vj] * 3 + dofj;

          entries.emplace_back(globalRow, globalCol, 1.0);
        }
      }
    }
  }

  if (flexibleRigidMotion) {
    for (int vi = 0; vi < (int)vertexIndices.size(); vi++) {
      for (int dofi = 0; dofi < 3; dofi++) {
        int globalRow = vertexIndices[vi] * 3 + dofi;
        for (int j = 0; j < 12; j++) {
          entries.emplace_back(globalRow, (int)restPositions.size() + j, 1.0);
          entries.emplace_back((int)restPositions.size() + j, globalRow, 1.0);
        }
      }
    }

    for (int i = 0; i < 12; i++) {
      for (int j = 0; j < 12; j++) {
        entries.emplace_back((int)restPositions.size() + i, (int)restPositions.size() + j, 1.0);
      }
    }

    hessTemplate.resize(restPositions.size() + 12, restPositions.size() + 12);
    hessTemplate.setFromTriplets(entries.begin(), entries.end());
    ES::buildEntryMap(hessTemplate, entryMap);

    allDOFs.resize(restPositions.size() + 12);
    std::iota(allDOFs.begin(), allDOFs.end(), 0);
  }
  else {
    hessTemplate.resize(restPositions.size(), restPositions.size());
    hessTemplate.setFromTriplets(entries.begin(), entries.end());
    ES::buildEntryMap(hessTemplate, entryMap);

    allDOFs.resize(restPositions.size());
    std::iota(allDOFs.begin(), allDOFs.end(), 0);
  }

  targetCenter = centerRest;
  targetRotation.setIdentity();
}

void MultipleVertexConstrainedRigidMotion::setDOFs(const std::vector<int> &dofs)
{
  PGO_ALOG((int)dofs.size() == (int)allDOFs.size());
  allDOFs = dofs;
}

void MultipleVertexConstrainedRigidMotion::setRestCenter(const double origin[3])
{
  if (origin) {
    centerRest = ES::V3d(origin[0], origin[1], origin[2]);
  }
  else {
    double wAll = 0;
    centerRest.setZero();
    for (size_t i = 0; i < vertexIndices.size(); i++) {
      centerRest += restPositions.segment<3>(vertexIndices[i] * 3) * vertexWeights[i];
      wAll += vertexWeights[i];
    }
    centerRest /= wAll;
  }

  for (size_t i = 0; i < vertexIndices.size(); i++) {
    ES::V3d qi = restPositions.segment<3>(vertexIndices[i] * 3) - centerRest;
    q[i] = qi;
  }
}

ES::V3d MultipleVertexConstrainedRigidMotion::computePosition(ES::ConstRefVecXd x, int vid) const
{
  if (isDisp) {
    return restPositions.segment<3>(vid * 3) + x.segment<3>(vid * 3);
  }
  else {
    return x.segment<3>(vid * 3);
  }
}

// tcur = sum wi/wall xi
// F = sum wi/wall * (xi - t) qi^T
// R = polar(F)
double MultipleVertexConstrainedRigidMotion::func(ES::ConstRefVecXd u) const
{
  // compute center
  ES::V3d center = ES::V3d::Zero();
  double wAll = 0;
  for (size_t i = 0; i < vertexIndices.size(); i++) {
    center += computePosition(u, vertexIndices[i]) * vertexWeights[i];
    wAll += vertexWeights[i];
  }
  center /= wAll;

  double eng = 0.0;

  if (flexibleRigidMotion) {
    eng += (center - u.segment<3>(restPositions.size())).squaredNorm() * coeff[1];
  }
  else {
    eng += (center - targetCenter).squaredNorm() * coeff[1];
  }

  if (coeff[0] > 0) {
    // compute rotation
    ES::M3d F = ES::M3d::Zero();
    for (size_t i = 0; i < vertexIndices.size(); i++) {
      ES::V3d qi = q[i];
      ES::V3d pi = computePosition(u, vertexIndices[i]) - center;
      F += ES::tensorProduct(pi, qi) * vertexWeights[i];
    }
    F /= wAll;

    ES::M3d S = PolarDecompositionDerivatives::computeS(F);
    ES::M3d R = PolarDecompositionDerivatives::computeR(F, S);

    if (flexibleRigidMotion) {
      ES::M3d tgtR = ES::Mp<const ES::M3d>(u.data() + restPositions.size() + 3);
      eng += (R - tgtR).squaredNorm() * coeff[0];
    }
    else {
      eng += (R - targetRotation).squaredNorm() * coeff[0];
    }
  }

  return eng * 0.5;
}

// tcur = sum wi/wall xi
// tcur = Z * x
//
// F = sum wi/wall * (xi - t) qi^T
//   = sum wi/wall * (xi - Zx) qi^T
//   = sum wi/wall * (Si x - Z x) qi^T
//   = sum wi/wall * (Si - Z) x qi^T
//   = sum wi/wall * Mi x qi^T
// Fkl = sum wi/wall * Mi_k x qi_l
// dFkl/dx = sum wi/wall * Mi_k I qi_l
//         = sum wi/wall * Mi_k qi_l

// R = polar(F)
// dR/dx = dR/dF * dF/dx

// E = 1/2 (R - Rbar)^2
// dx = (R - Rbar) : dR/dF : dF/dx

// 1/2 (t-tbar)^2
// dx = (Zx - tbar)^T Z

void MultipleVertexConstrainedRigidMotion::gradient(ES::ConstRefVecXd u, ES::RefVecXd grad) const
{
  // compute center
  ES::V3d center = ES::V3d::Zero();
  double wAll = 0;
  for (size_t i = 0; i < vertexIndices.size(); i++) {
    center += computePosition(u, vertexIndices[i]) * vertexWeights[i];
    wAll += vertexWeights[i];
  }
  center /= wAll;

  // compute translation gradient
  if (flexibleRigidMotion) {
    // x
    ES::V3d diff = center - u.segment<3>(restPositions.size());
    ES::mv(Z, diff, grad.head(restPositions.size()), 1);
    // t
    grad.segment<3>(restPositions.size()) = -diff;

    grad.head(restPositions.size() + 3) *= coeff[1];
  }
  else {
    ES::V3d diff = center - targetCenter;
    ES::mv(Z, diff, grad, 1);
    grad *= coeff[1];
  }

  // rotation
  if (coeff[0] > 0) {
    // compute deformation gradient
    ES::M3d F = ES::M3d::Zero();
    for (size_t i = 0; i < vertexIndices.size(); i++) {
      ES::V3d qi = q[i];
      ES::V3d pi = computePosition(u, vertexIndices[i]) - center;
      F += ES::tensorProduct(pi, qi) * vertexWeights[i];
    }
    F /= wAll;

    // compute polar decomposition and derivatives
    ES::M3d R, S, dRdF[9], dSdF[9];
    PolarDecompositionDerivatives::compute(F, R, S, dSdF, dRdF, nullptr, nullptr);

    // E = 1/2 (R - Rbar)^2
    // dx = (R - Rbar) : dR/dF : dF/dx
    ES::M3d Rdiff;
    if (flexibleRigidMotion) {
      Rdiff = R - ES::Mp<const ES::M3d>(u.data() + restPositions.size() + 3);
    }
    else {
      Rdiff = R - targetRotation;
    }

    // compute dE/dF
    double dEdFi[9];
    for (int i = 0; i < 9; i++) {
      dEdFi[i] = Rdiff.cwiseProduct(dRdF[i]).sum();
    }

    ES::VXd dFidx;
    for (int row = 0; row < 3; row++) {
      for (int col = 0; col < 3; col++) {
        int Fdofi = col * 3 + row;

        dFidx.setZero(Z.cols());
        for (size_t i = 0; i < vertexIndices.size(); i++) {
          ES::V3d qi = q[i];

          for (ES::SpMatD::InnerIterator it(Mi[i], row); it; ++it) {
            dFidx(it.col()) += it.value() * qi[col] * vertexWeights[i] / wAll;
          }
        }

        grad.head(restPositions.size()) += dFidx * dEdFi[Fdofi] * coeff[0];
      }
    }

    if (flexibleRigidMotion) {
      // compute dE/dRbar
      grad.tail<9>() = ES::Mp<const ES::V9d>(Rdiff.data()) * -coeff[0];
    }
  }
}

// tcur = sum wi/wall xi
// tcur = Z * x
//
// F = sum wi/wall * (xi - t) qi^T
//   = sum wi/wall * (xi - Zx) qi^T
//   = sum wi/wall * (Si x - Z x) qi^T
//   = sum wi/wall * (Si - Z) x qi^T
//   = sum wi/wall * Mi x qi^T
//
// dF/dxk = sum wi/wall * Mi.col(k) qi^T

// R = polar(F)
// dR/dx = dR/dF * dF/dx

// E = 1/2 (R - Rbar)^2
// dx = (R - Rbar) : dR/dF : dF/dx
// dx2 = dFdx^T : dRdF^T : dR/dF : dF/dx + (R - Rbar) : d2R/dF2 dFdx : dF/dx
// dxidxj = dR/dxj : dR/dxi + (R - Rbar) : (d2R/dF2 dFdxj) : dF/dxi

void MultipleVertexConstrainedRigidMotion::hessian(ES::ConstRefVecXd u, ES::SpMatD &hess) const
{
  memset(hess.valuePtr(), 0, sizeof(double) * hess.nonZeros());

  // 1/2 (t-tbar)^2
  // dx = (Zx - tbar)^T Z
  // dx2 = ZTZ
  // dx dtbar = -Z
  // dtbar2 = I
  // for (ES::IDX rowi = 0; rowi < ZTZ.rows(); rowi++) {
  tbb::parallel_for(0, (int)ZTZ.rows(), [&](int rowi) {
    for (ES::SpMatD::InnerIterator it(ZTZ, rowi); it; ++it) {
      auto iter = entryMap.find(std::pair<int, int>((int)it.row(), (int)it.col()));
      PGO_ALOG(iter != entryMap.end());

      hess.valuePtr()[iter->second] = it.value() * coeff[1];
    }
  });

  if (flexibleRigidMotion) {
    for (ES::IDX rowi = 0; rowi < Z.rows(); rowi++) {
      for (ES::SpMatD::InnerIterator it(Z, rowi); it; ++it) {
        auto iter = entryMap.find(std::pair<int, int>((int)it.row() + (int)restPositions.size(), (int)it.col()));
        PGO_ALOG(iter != entryMap.end());
        hess.valuePtr()[iter->second] += -it.value() * coeff[1];

        iter = entryMap.find(std::pair<int, int>((int)it.col(), (int)it.row() + (int)restPositions.size()));
        PGO_ALOG(iter != entryMap.end());
        hess.valuePtr()[iter->second] = -it.value() * coeff[1];
      }
    }

    for (int i = 0; i < 3; i++) {
      auto iter = entryMap.find(std::pair<int, int>(i + (int)restPositions.size(), i + (int)restPositions.size()));
      PGO_ALOG(iter != entryMap.end());
      hess.valuePtr()[iter->second] = coeff[1];
    }
  }

  if (coeff[0] > 0) {
    ES::V3d center = ES::V3d::Zero();
    double wAll = 0;
    for (size_t i = 0; i < vertexIndices.size(); i++) {
      center += computePosition(u, vertexIndices[i]) * vertexWeights[i];
      wAll += vertexWeights[i];
    }
    center /= wAll;

    ES::M3d F = ES::M3d::Zero();
    for (size_t i = 0; i < vertexIndices.size(); i++) {
      ES::V3d qi = q[i];
      ES::V3d pi = computePosition(u, vertexIndices[i]) - center;
      F += ES::tensorProduct(pi, qi) * vertexWeights[i];
    }
    F /= wAll;

    ES::M3d R, S, dRdF[9], dSdF[9], d2RdF2[9][9], d2SdF2[9][9];
    PolarDecompositionDerivatives::compute(F, R, S, dSdF, dRdF, d2SdF2, d2RdF2);

    // E = 1/2 (R - Rbar)^2
    // dx = (R - Rbar) : dR/dF : dF/dx
    // dx2 = dFdx^T : dRdF^T : dR/dF : dF/dx + (R - Rbar) : d2R/dF2 dF/dx : dF/dx
    // dxidxj = dR/dxj : dR/dxi + (R - Rbar) : (d2R/dF2 dFdxj) : dF/dxi

    // compute dF/dxi
    ES::EigenArray<ES::M3d> dFdxi(vertexIndices.size() * 3, ES::M3d::Zero());
    for (size_t i = 0; i < vertexIndices.size(); i++) {
      ES::V3d qi = q[i];

      for (ES::SpMatD::InnerIterator it(Mi[i], 0); it; ++it) {
        for (int dof = 0; dof < 3; dof++) {
          ES::V3d Mvec = Mi[i].block(0, it.col() + dof, 3, 1);
          ES::M3d Mvecq = ES::tensorProduct(Mvec, qi);

          auto iter = vertexIDtoIndices.find((int)it.col() / 3);
          PGO_ALOG(iter != vertexIDtoIndices.end());
          int idx = iter->second;

          dFdxi[idx * 3 + dof] += Mvecq * (vertexWeights[i] / wAll);
        }
      }
    }

    // compute dR/dxi
    ES::EigenArray<ES::M3d> dRdxi(vertexIndices.size() * 3, ES::M3d::Zero());
    for (int vi = 0; vi < (int)vertexIndices.size(); vi++) {
      for (int dofi = 0; dofi < 3; dofi++) {
        for (int i = 0; i < 9; i++) {
          dRdxi[vi * 3 + dofi] += dRdF[i] * dFdxi[vi * 3 + dofi].data()[i];
        }
      }
    }

    ES::M3d Rdiff;
    if (flexibleRigidMotion) {
      Rdiff = R - ES::Mp<const ES::M3d>(u.data() + restPositions.size() + 3);
    }
    else {
      Rdiff = R - targetRotation;
    }

    // for (int vi = 0; vi < (int)vertexIndices.size(); vi++) {
    tbb::parallel_for(0, (int)vertexIndices.size(), [&](int vi) {
      // for (int vj = 0; vj < (int)vertexIndices.size(); vj++) {
      tbb::parallel_for(0, (int)vertexIndices.size(), [&](int vj) {
        for (int dofi = 0; dofi < 3; dofi++) {
          for (int dofj = 0; dofj < 3; dofj++) {
            double val = dRdxi[vi * 3 + dofi].cwiseProduct(dRdxi[vj * 3 + dofj]).sum();

            // (d2R / dF2 dFdxj): dF/dxi
            ES::M3d d2Rdxidxj;
            d2Rdxidxj.setZero();

            for (int k = 0; k < 9; k++) {
              for (int l = 0; l < 9; l++) {
                d2Rdxidxj += d2RdF2[k][l] * dFdxi[vj * 3 + dofj].data()[k] * dFdxi[vi * 3 + dofi].data()[l];
              }
            }

            val += Rdiff.cwiseProduct(d2Rdxidxj).sum();

            int globalRow = vertexIndices[vi] * 3 + dofi;
            int globalCol = vertexIndices[vj] * 3 + dofj;

            auto iter = entryMap.find(std::pair<int, int>(globalRow, globalCol));
            PGO_ALOG(iter != entryMap.end());

            hess.valuePtr()[iter->second] += val * coeff[0];
          }
        }
      });
    });

    if (flexibleRigidMotion) {
      // dxi dRbar = -dRbar/dRbar : dR/dxi
      // for (int vi = 0; vi < (int)vertexIndices.size(); vi++) {
      tbb::parallel_for(0, (int)vertexIndices.size(), [&](int vi) {
        for (int dofi = 0; dofi < 3; dofi++) {
          for (int dofj = 0; dofj < 9; dofj++) {
            int globalRow = vertexIndices[vi] * 3 + dofi;
            int globalCol = int(restPositions.size()) + 3 + dofj;
            double v = -dRdxi[vi * 3 + dofi].data()[dofj];

            auto iter = entryMap.find(std::pair<int, int>(globalRow, globalCol));
            PGO_ALOG(iter != entryMap.end());
            hess.valuePtr()[iter->second] += v * coeff[0];

            iter = entryMap.find(std::pair<int, int>(globalCol, globalRow));
            PGO_ALOG(iter != entryMap.end());
            hess.valuePtr()[iter->second] += v * coeff[0];
          }
        }
      });
      // dRbar dRbar = I
      for (int i = 0; i < 9; i++) {
        auto iter = entryMap.find(std::pair<int, int>(i + (int)restPositions.size() + 3, i + (int)restPositions.size() + 3));
        PGO_ALOG(iter != entryMap.end());
        hess.valuePtr()[iter->second] = coeff[0];
      }
    }
  }
}