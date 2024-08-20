/*
author: Bohan Wang
copyright to USC,MIT
*/

#include "prescribedPrincipleStressConstraintFunctions.h"
#include "tetMeshDeformationModel.h"
#include "deformationModelManager.h"
#include "simulationMesh.h"

#include "svdDerivatives.h"
#include "pgoLogging.h"

using namespace pgo;
using namespace pgo::SolidDeformationModel;

PrescribedPrincipleStressConstraintFunctions::PrescribedPrincipleStressConstraintFunctions(int nAll, int doff, int numElements, const int *elementIDs, const DeformationModelManager *tmdmm):
  ConstraintFunctions(nAll), dofStart(doff), tetMeshDMM(tmdmm)
{
  elements.assign(elementIDs, elementIDs + numElements);
  targetPrincipleStress.resize(numElements * 3);

  elementCacheData.assign(numElements, nullptr);
  for (int ei = 0; ei < numElements; ei++) {
    elementCacheData[ei] = tetMeshDMM->getDeformationModel(elements[ei])->allocateCacheData();
  }

  std::vector<ES::TripletD> entries;
  for (int i = 0; i < numElements; i++) {
    for (int r = 0; r < 3; r++) {
      for (int j = 0; j < 4; j++) {
        for (int dof = 0; dof < 3; dof++) {
          entries.emplace_back(i * 3 + r, dofStart + tetMeshDMM->getMesh()->getVertexIndex(elements[i], j) * 3 + dof, 1.0);
        }
      }
    }
  }
  jacobianTemplate.resize(numElements * 3, nAll);
  jacobianTemplate.setFromTriplets(entries.begin(), entries.end());
  ES::buildEntryMap(jacobianTemplate, jacEntries);

  entries.clear();
  for (int i = 0; i < numElements; i++) {
    for (int vi = 0; vi < 4; vi++) {
      for (int dofi = 0; dofi < 3; dofi++) {
        int row = tetMeshDMM->getMesh()->getVertexIndex(elements[i], vi) * 3 + dofi;
        for (int vj = 0; vj < 4; vj++) {
          for (int dofj = 0; dofj < 3; dofj++) {
            int col = tetMeshDMM->getMesh()->getVertexIndex(elements[i], vj) * 3 + dofj;
            entries.emplace_back(dofStart + row, dofStart + col, 1.0);
          }
        }
      }
    }
  }

  lambdahTemplate.resize(nAll, nAll);
  lambdahTemplate.setFromTriplets(entries.begin(), entries.end());
  ES::buildEntryMap(lambdahTemplate, hessEntries);

  hessLocks = std::vector<tbb::spin_mutex>(nAll);
}

// g = S(P) - Pbar
void PrescribedPrincipleStressConstraintFunctions::func(ES::ConstRefVecXd x, ES::RefVecXd g) const
{
  for (int i = 0; i < (int)elements.size(); i++) {
    ES::V3d Phat = targetPrincipleStress.segment<3>(i * 3);
    int eleID = elements[i];
    ES::V18d localp;
    for (int j = 0; j < tetMeshDMM->getMesh()->getNumElementVertices(); j++) {
      int vid = tetMeshDMM->getMesh()->getVertexIndex(eleID, j);
      ES::V3d vtxp;
      xToPosFunc(x.segment<3>(dofStart + vid * 3), dofStart + vid * 3, vtxp);
      localp.segment<3>(j * 3) = vtxp;
    }

    const DeformationModel *fem = tetMeshDMM->getDeformationModel(eleID);
    ES::V12d plasticParam, elasticParam;
    plasticParam.setZero();
    elasticParam.setZero();

    fem->prepareData(localp.data(), plasticParam.data(), elasticParam.data(), elementCacheData[i]);

    ES::M3d P;
    fem->computeP(elementCacheData[i], 0, P.data());

    ES::V3d S;
    NonlinearOptimization::SVDDerivatives::unorderedSquareMatrixSVD3(P, S);

    g.segment<3>(i * 3) = S - Phat;
  }
}

void PrescribedPrincipleStressConstraintFunctions::computeForceFromTargetPHat(ES::ConstRefVecXd x, ES::RefVecXd fext) const
{
  fext.setZero();

  for (int i = 0; i < (int)elements.size(); i++) {
    ES::V3d Phat = targetPrincipleStress.segment<3>(i * 3);
    int eleID = elements[i];
    ES::V18d localp;
    for (int j = 0; j < tetMeshDMM->getMesh()->getNumElementVertices(); j++) {
      int vid = tetMeshDMM->getMesh()->getVertexIndex(eleID, j);
      ES::V3d vtxp;
      xToPosFunc(x.segment<3>(dofStart + vid * 3), dofStart + vid * 3, vtxp);
      localp.segment<3>(j * 3) = vtxp;
    }

    const DeformationModel *fem = tetMeshDMM->getDeformationModel(eleID);
    ES::V12d plasticParam, elasticParam;
    plasticParam.setZero();
    elasticParam.setZero();

    fem->prepareData(localp.data(), plasticParam.data(), elasticParam.data(), elementCacheData[i]);

    ES::M3d P;
    fem->computeP(elementCacheData[i], 0, P.data());

    ES::V3d S;
    ES::M3d U, V;
    NonlinearOptimization::SVDDerivatives::unorderedSquareMatrixSVD3(P, S, &U, &V);

    ES::M3d P1 = U * Phat.asDiagonal() * V.transpose();
    // ES::M3d P1 = Phat.asDiagonal();
    ES::V12d f;
    f.setZero();

    const TetMeshDeformationModel *tetFEM = dynamic_cast<const TetMeshDeformationModel *>(fem);
    tetFEM->computeForceFromP(elementCacheData[i], P1.data(), f.data());

    for (int j = 0; j < tetMeshDMM->getMesh()->getNumElementVertices(); j++) {
      int vid = tetMeshDMM->getMesh()->getVertexIndex(eleID, j);
      fext.segment<3>(vid * 3) = f.segment<3>(j * 3);
    }
  }
}

double PrescribedPrincipleStressConstraintFunctions::computeSurfaceNormalTractionFromElement(ES::ConstRefVecXd x, const ES::V3d &n, int eleID) const
{
  const DeformationModel *fem = tetMeshDMM->getDeformationModel(eleID);
  DeformationModel::CacheData *cache = fem->allocateCacheData();
  ES::V12d plasticParam, elasticParam;
  plasticParam.setZero();
  elasticParam.setZero();

  ES::V18d localp;
  for (int j = 0; j < tetMeshDMM->getMesh()->getNumElementVertices(); j++) {
    int vid = tetMeshDMM->getMesh()->getVertexIndex(eleID, j);
    ES::V3d vtxp;
    xToPosFunc(x.segment<3>(dofStart + vid * 3), dofStart + vid * 3, vtxp);
    localp.segment<3>(j * 3) = vtxp;
  }
  fem->prepareData(localp.data(), plasticParam.data(), elasticParam.data(), cache);

  ES::M3d P;
  fem->computeP(cache, 0, P.data());

  fem->freeCacheData(cache);

  return (P * n).dot(n);
}

// dg/dx = dS/dP dP/dF dF/dx
void PrescribedPrincipleStressConstraintFunctions::jacobian(ES::ConstRefVecXd x, ES::SpMatD &jac) const
{
  for (int i = 0; i < (int)elements.size(); i++) {
    ES::V3d Phat = targetPrincipleStress.segment<3>(i * 3);
    int eleID = elements[i];
    ES::V18d localp;
    for (int j = 0; j < tetMeshDMM->getMesh()->getNumElementVertices(); j++) {
      int vid = tetMeshDMM->getMesh()->getVertexIndex(eleID, j);
      ES::V3d vtxp;
      xToPosFunc(x.segment<3>(dofStart + vid * 3), dofStart + vid * 3, vtxp);
      localp.segment<3>(j * 3) = vtxp;
    }

    const DeformationModel *fem = tetMeshDMM->getDeformationModel(eleID);
    ES::V12d plasticParam, elasticParam;
    plasticParam.setZero();
    elasticParam.setZero();

    fem->prepareData(localp.data(), plasticParam.data(), elasticParam.data(), elementCacheData[i]);

    ES::M3d P;
    fem->computeP(elementCacheData[i], 0, P.data());

    ES::V3d S;
    ES::M3d U, V;
    NonlinearOptimization::SVDDerivatives::unorderedSquareMatrixSVD3(P, S, &U, &V);

    ES::V3d dSdPi[9];
    for (int i = 0; i < 9; i++) {
      ES::M3d dP;
      dP.setZero();
      dP.data()[i] = 1.0;
      NonlinearOptimization::SVDDerivatives::unorderedSquareMatrixSVD3Derivatices(P, S, U, V, dP, dSdPi[i]);
    }

    ES::M9d dPdF;
    ES::M9x12d dFdx;
    fem->computedPdF(elementCacheData[i], 0, dPdF.data());
    fem->computedFdx(elementCacheData[i], 0, dFdx.data());

    ES::M9x12d dPdx = dPdF * dFdx;

    Eigen::Matrix<double, 3, 12> dSdx;
    dSdx.setZero();

    // dS/dPi dPi/dx
    for (int i = 0; i < 9; i++) {
      Eigen::Matrix<double, 1, 12> v = dPdx.row(i);
      dSdx += dSdPi[i] * v;
    }

    for (int ci = 0; ci < 12; ci++) {
      for (int ri = 0; ri < 3; ri++) {
        int vid = ci / 3;
        int dof = ci % 3;

        auto it = jacEntries.find(std::make_pair(i * 3 + ri, tetMeshDMM->getMesh()->getVertexIndex(eleID, vid) * 3 + dof));
        PGO_ALOG(it != jacEntries.end());

        jac.valuePtr()[it->second] = dSdx(ri, ci);
      }
    }
  }
}

// dg/dx = dSdP dPdx
// d2g/dx2 = dPdx d2SdP2 dPdx + dSdP d2Pdx2
void PrescribedPrincipleStressConstraintFunctions::hessian(ES::ConstRefVecXd x, ES::ConstRefVecXd lambda, ES::SpMatD &hess) const
{
  for (int ei = 0; ei < (int)elements.size(); ei++) {
    int eleID = elements[ei];
    ES::V18d localp;
    for (int j = 0; j < tetMeshDMM->getMesh()->getNumElementVertices(); j++) {
      int vid = tetMeshDMM->getMesh()->getVertexIndex(eleID, j);
      ES::V3d vtxp;
      xToPosFunc(x.segment<3>(dofStart + vid * 3), dofStart + vid * 3, vtxp);
      localp.segment<3>(j * 3) = vtxp;
    }

    const DeformationModel *fem = tetMeshDMM->getDeformationModel(eleID);
    ES::V12d plasticParam, elasticParam;
    plasticParam.setZero();
    elasticParam.setZero();

    fem->prepareData(localp.data(), plasticParam.data(), elasticParam.data(), elementCacheData[ei]);

    ES::M3d P;
    fem->computeP(elementCacheData[ei], 0, P.data());

    ES::V3d S;
    ES::M3d U, V;
    NonlinearOptimization::SVDDerivatives::unorderedSquareMatrixSVD3(P, S, &U, &V);

    ES::V3d d2SdPidPj[9][9], dSdPi[9];
    for (int i = 0; i < 9; i++) {
      for (int j = 0; j < 9; j++) {
        ES::M3d dP;
        dP.setZero();
        dP.data()[i] = 1.0;

        ES::M3d d2P;
        d2P.setZero();

        NonlinearOptimization::SVDDerivatives::unorderedSquareMatrixSVD3Derivatices(P, S, U, V,
          dP, dSdPi[i], nullptr, nullptr,
          &d2P, &d2SdPidPj[i][j]);
      }
    }

    ES::M9d dPdF;
    ES::M9x12d dFdx;
    fem->computedPdF(elementCacheData[ei], 0, dPdF.data());
    fem->computedFdx(elementCacheData[ei], 0, dFdx.data());
    ES::M9x12d dPdx = dPdF * dFdx;

    ES::M12d hessLocal;
    hessLocal.setZero();

    ES::V3d lambdaLocal = lambda.segment<3>(ei * 3);

    // d2g/dx2 = dPdx d2SdP2 dPdx + dSdP d2Pdx2
    // second term is omit
    for (int si = 0; si < 3; si++) {
      for (int i = 0; i < 9; i++) {
        ES::V12d dPidx = dPdx.row(i);

        for (int j = 0; j < 9; j++) {
          ES::V12d dPjdx = dPdx.row(j);

          ES::M12d h;
          ES::tensorProduct(h, dPjdx, dPidx);
          hessLocal += d2SdPidPj[i][j][si] * lambdaLocal[si] * h;
        }
      }
    }

    for (int ci = 0; ci < 12; ci++) {
      for (int ri = 0; ri < 12; ri++) {
        int vi = ci / 3;
        int dofi = ci % 3;

        int vj = ri / 3;
        int dofj = ri % 3;

        int globalRow = tetMeshDMM->getMesh()->getVertexIndex(eleID, vj) * 3 + dofj;
        int globalCol = tetMeshDMM->getMesh()->getVertexIndex(eleID, vi) * 3 + dofi;

        auto it = hessEntries.find(std::make_pair(globalRow, globalCol));
        PGO_ALOG(it != hessEntries.end());

        hessLocks[globalRow].lock();

        hess.valuePtr()[it->second] += hessLocal(ri, ci);

        hessLocks[globalRow].unlock();
      }
    }
  }
}
