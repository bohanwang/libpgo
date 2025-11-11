/*
author: Bohan Wang
copyright to USC,MIT,NUS
*/

#include "tetVolumeConstraintFunctions.h"
#include "tetMeshDeformationModel.h"

#include "determinantDerivatives.h"
#include "tetMesh.h"
#include "pgoLogging.h"
#include "EigenSupport.h"

#include <tbb/parallel_for.h>

using namespace pgo;
using namespace pgo::SolidDeformationModel;
namespace ES = pgo::EigenSupport;

TetVolumeConstraintFunctions::TetVolumeConstraintFunctions(const SimulationMesh *tm, int nAll, const ES::VXd *rp, const ES::M3Xd *DmInvIn):
  NonlinearOptimization::ConstraintFunctions(nAll)
{
  tetMesh = tm;
  nele = tetMesh->getNumElements();

  if (DmInvIn) {
    setDmInv(*DmInvIn);
  }
  else {
    DmInv.resize(3, nele * 3);
    dFdx.assign(nele, ES::M9x12d::Zero());

    for (int ei = 0; ei < nele; ei++) {
      ES::M3d Dm;

      ES::V12d xlocal;
      for (int i = 0; i < 4; i++) {
        tetMesh->getVertex(ei, i, xlocal.data() + i * 3);
      }
      TetMeshDeformationModel::computeDm(xlocal.data(), Dm.data());

      DmInv.block<3, 3>(0, ei * 3) = Dm.fullPivHouseholderQr().inverse();

      TetMeshDeformationModel::compute_dF_dx(DmInv.data() + ei * 9, dFdx[ei].data());
    }
  }

  if (rp) {
    restPosition = rp;
  }

  std::vector<ES::TripletD> entries;
  for (int ei = 0; ei < nele; ei++) {
    for (int i = 0; i < 4; i++) {
      entries.emplace_back(ei, tetMesh->getVertexIndex(ei, i) * 3, 1.0);
      entries.emplace_back(ei, tetMesh->getVertexIndex(ei, i) * 3 + 1, 1.0);
      entries.emplace_back(ei, tetMesh->getVertexIndex(ei, i) * 3 + 2, 1.0);
    }
  }
  jacobianTemplate.resize(nele, nAll);
  jacobianTemplate.setFromTriplets(entries.begin(), entries.end());

  entries.clear();
  for (int ei = 0; ei < nele; ei++) {
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        for (int k = 0; k < 3; k++) {
          for (int l = 0; l < 3; l++) {
            entries.emplace_back(tetMesh->getVertexIndex(ei, i) * 3 + k, tetMesh->getVertexIndex(ei, j) * 3 + l, 1.0);
          }
        }
      }
    }
  }
  lambdahTemplate.resize(nAll, nAll);
  lambdahTemplate.setFromTriplets(entries.begin(), entries.end());

  jacobianIndices.assign(nele, JacIndex::Constant(-1));
  for (int ei = 0; ei < nele; ei++) {
    JacIndex idx;
    for (int i = 0; i < 4; i++) {
      idx[i] = ES::findEntryOffset(jacobianTemplate, ei, tetMesh->getVertexIndex(ei, i) * 3);
      PGO_ALOG(idx[i] >= 0);
    }
    jacobianIndices[ei] = idx;
  }

  hessIndices.assign(nele, HessIndex::Constant(-1));
  for (int ei = 0; ei < nele; ei++) {
    HessIndex idx;

    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        for (int k = 0; k < 3; k++) {
          for (int l = 0; l < 3; l++) {
            int row = tetMesh->getVertexIndex(ei, i) * 3 + k;
            int col = tetMesh->getVertexIndex(ei, j) * 3 + l;

            int localRow = i * 3 + k;
            int localCol = j * 3 + l;

            idx(localRow, localCol) = ES::findEntryOffset(lambdahTemplate, row, col);
            PGO_ALOG(idx(localRow, localCol) >= 0);
          }
        }
      }
    }
    hessIndices[ei] = idx;
  }
  hessLocks = std::vector<tbb::spin_mutex>(lambdahTemplate.nonZeros());
  elementFlags.assign(nele, 1);
}

TetVolumeConstraintFunctions::~TetVolumeConstraintFunctions()
{
}

void TetVolumeConstraintFunctions::setDmInv(const ES::M3Xd &DmInv_)
{
  DmInv = DmInv_;
  for (int ei = 0; ei < nele; ei++) {
    TetMeshDeformationModel::compute_dF_dx(DmInv.data() + ei * 9, dFdx[ei].data());
  }
}

void TetVolumeConstraintFunctions::func(ES::ConstRefVecXd x, ES::RefVecXd g) const
{
  // for (int ei = 0; ei < nele; ei++) {
  tbb::parallel_for(0, nele, [&](int ei) {
    if (elementFlags[ei] == 0) {
      g[ei] = 1;
      return;
    }

    const int *vertexIndices = tetMesh->getVertexIndices(ei);
    ES::V12d xlocal;
    if (restPosition) {
      for (int i = 0; i < 4; i++) {
        xlocal.segment<3>(i * 3) = restPosition->segment<3>(vertexIndices[i] * 3) + x.segment<3>(vertexIndices[i] * 3);
      }
    }
    else {
      for (int i = 0; i < 4; i++) {
        xlocal.segment<3>(i * 3) = x.segment<3>(vertexIndices[i] * 3);
      }
    }

    ES::M3d Ds;
    TetMeshDeformationModel::computeDs(xlocal.data(), Ds.data());

    ES::M3d DmInvLocal = DmInv.block<3, 3>(0, ei * 3);
    ES::M3d F = Ds * DmInvLocal;

    g[ei] = F.determinant();
  });
}

void TetVolumeConstraintFunctions::jacobian(ES::ConstRefVecXd x, ES::SpMatD &jac) const
{
  tbb::parallel_for(0, nele, [&](int ei) {
    if (elementFlags[ei] == 0) {
      for (int i = 0; i < 4; i++) {
        ES::IDX offsetStart = jacobianIndices[ei][i];
        jac.valuePtr()[offsetStart + 0] = 0;
        jac.valuePtr()[offsetStart + 1] = 0;
        jac.valuePtr()[offsetStart + 2] = 0;
      }

      return;
    }

    const int *vertexIndices = tetMesh->getVertexIndices(ei);

    ES::V12d xlocal;
    if (restPosition) {
      for (int i = 0; i < 4; i++) {
        xlocal.segment<3>(i * 3) = restPosition->segment<3>(vertexIndices[i] * 3) + x.segment<3>(vertexIndices[i] * 3);
      }
    }
    else {
      for (int i = 0; i < 4; i++) {
        xlocal.segment<3>(i * 3) = x.segment<3>(vertexIndices[i] * 3);
      }
    }

    ES::M3d Ds;
    TetMeshDeformationModel::computeDs(xlocal.data(), Ds.data());

    ES::M3d DmInvLocal = DmInv.block<3, 3>(0, ei * 3);

    ES::M3d F = Ds * DmInvLocal;

    // ddetF/dF * dF/dx
    ES::M3d ddetA_dA;
    NonlinearOptimization::Determinant::Dim3::ddetA_dA(F.data(), ddetA_dA.data());

    ES::V9d ddetA_dA_vec = Eigen::Map<const ES::V9d>(ddetA_dA.data());
    ES::V12d ddetA_dx = dFdx[ei].transpose() * ddetA_dA_vec;

    for (int i = 0; i < 4; i++) {
      ES::IDX offsetStart = jacobianIndices[ei][i];
      jac.valuePtr()[offsetStart + 0] = ddetA_dx[i * 3];
      jac.valuePtr()[offsetStart + 1] = ddetA_dx[i * 3 + 1];
      jac.valuePtr()[offsetStart + 2] = ddetA_dx[i * 3 + 2];
    }
  });
}

void TetVolumeConstraintFunctions::hessian(ES::ConstRefVecXd x, ES::ConstRefVecXd lambda, ES::SpMatD &hess) const
{
  memset(hess.valuePtr(), 0, sizeof(double) * hess.nonZeros());

  tbb::parallel_for(0, nele, [&](int ei) {
    if (elementFlags[ei] == 0)
      return;

    const int *vertexIndices = tetMesh->getVertexIndices(ei);

    ES::V12d xlocal;
    if (restPosition) {
      for (int i = 0; i < 4; i++) {
        xlocal.segment<3>(i * 3) = restPosition->segment<3>(vertexIndices[i] * 3) + x.segment<3>(vertexIndices[i] * 3);
      }
    }
    else {
      for (int i = 0; i < 4; i++) {
        xlocal.segment<3>(i * 3) = x.segment<3>(vertexIndices[i] * 3);
      }
    }

    ES::M3d Ds;
    TetMeshDeformationModel::computeDs(xlocal.data(), Ds.data());

    ES::M3d DmInvLocal = DmInv.block<3, 3>(0, ei * 3);

    ES::M3d F = Ds * DmInvLocal;

    // d detF / dx = d detF / dF * dF/dx
    // d2 detF / dx2 = d (d detF / dF * dF/dx) / dx
    //               = (dF/dx)^T * d2 detF/ dF2 * dF/dx

    ES::M9d d2detA_dA2;
    NonlinearOptimization::Determinant::Dim3::d2detA_dA2(F.data(), d2detA_dA2.data());

    ES::M12d d2detA_dx2 = dFdx[ei].transpose() * d2detA_dA2 * dFdx[ei] * lambda[ei];

    for (int col = 0; col < 4; col++) {
      for (int dofj = 0; dofj < 3; dofj++) {
        for (int row = 0; row < 4; row++) {
          for (int dofi = 0; dofi < 3; dofi++) {
            int localRow = row * 3 + dofi;
            int localCol = col * 3 + dofj;

            ES::IDX offset = hessIndices[ei](localRow, localCol);

            hessLocks[offset].lock();

            hess.valuePtr()[offset] += d2detA_dx2(localRow, localCol);

            hessLocks[offset].unlock();
          }
        }
      }
    }
  });
}
