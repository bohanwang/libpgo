/*
author: Bohan Wang
copyright to USC,MIT
*/

#include "deformationModelAssemblerElastic.h"
#include "deformationModelAssemblerElasticImpl.h"
#include "deformationModel.h"
#include "deformationModelManager.h"
#include "simulationMesh.h"
#include "plasticModel.h"

#include "pgoLogging.h"
#include "fmtEigen.h"

#include <tbb/parallel_for.h>

using namespace pgo;
using namespace pgo::SolidDeformationModel;

DeformationModelAssemblerElastic::DeformationModelAssemblerElastic(const DeformationModelManager *dm, const double *basisMatrix, int _numHandles, const double *elementFlags):
  DeformationModelAssembler(dm, DeformationModelAssemblingType::AT_ELASTIC, elementFlags)
{
  numPlasticParams = deformationModelManager->getNumPlasticParameters();
  numElasticParams = deformationModelManager->getNumElasticParameters();
  numHandles = _numHandles;

  SPDLOG_LOGGER_INFO(Logging::lgr(), "Assembler parameter:{},{},{}", numElasticParams, numPlasticParams, numHandles);

  eimpl = new DeformationModelAssemblerElasticImpl(*impl);
  delete impl;
  impl = eimpl;

  if (basisMatrix && numHandles) {
    ES::MXd W = Eigen::Map<const ES::MXd>(basisMatrix, nele, numHandles);
    eimpl->elementWeightMatricesCompressed.assign(nele, ES::VXd::Zero(numHandles));
    eimpl->elementWeightMatricesPlastic.resize(nele, ES::MXd::Zero(numPlasticParams, numPlasticParams * numHandles));

    if (numElasticParams > 0)
      eimpl->elementWeightMatricesElastic.resize(nele, ES::MXd::Zero(numElasticParams, numElasticParams * numHandles));

    for (int ele = 0; ele < nele; ele++) {
      for (int hi = 0; hi < numHandles; hi++) {
        eimpl->elementWeightMatricesPlastic[ele].block(0, hi * numPlasticParams, numPlasticParams, numPlasticParams) =
          ES::MXd::Identity(numPlasticParams, numPlasticParams) * W(ele, hi);

        if (numElasticParams > 0)
          eimpl->elementWeightMatricesElastic[ele].block(0, hi * numElasticParams, numElasticParams, numElasticParams) =
            ES::MXd::Identity(numElasticParams, numElasticParams) * W(ele, hi);
      }

      eimpl->elementWeightMatricesCompressed[ele] = W.row(ele);
    }
  }

  std::vector<ES::TripletD> entries;
  for (int ele = 0; ele < nele; ele++) {
    const int *vertexIndices = deformationModelManager->getMesh()->getVertexIndices(ele);

    // add element to the stiffness matrix
    for (int vi = 0; vi < neleVtx; vi++) {
      for (int vj = 0; vj < neleVtx; vj++) {
        for (int dofi = 0; dofi < 3; dofi++) {
          for (int dofj = 0; dofj < 3; dofj++) {
            if (vertexIndices[vi] >= 0 && vertexIndices[vj] >= 0)
              entries.emplace_back(vertexIndices[vi] * 3 + dofi, vertexIndices[vj] * 3 + dofj, 1.0);
          }
        }
      }
    }
  }

  eimpl->KTemplateHandle.hess.resize(n3, n3);
  eimpl->KTemplateHandle.hess.setFromTriplets(entries.begin(), entries.end());

  eimpl->inverseIndices.resize(nele);
  for (int ele = 0; ele < nele; ele++) {
    const int *vertexIndices = deformationModelManager->getMesh()->getVertexIndices(ele);
    DeformationModelAssemblerElasticImpl::IndexMatrix idxM = DeformationModelAssemblerElasticImpl::IndexMatrix::Constant(-1);

    // upper-left block
    for (int vi = 0; vi < neleVtx; vi++) {
      for (int vj = 0; vj < neleVtx; vj++) {
        for (int dofi = 0; dofi < 3; dofi++) {
          for (int dofj = 0; dofj < 3; dofj++) {
            int localRow = vi * 3 + dofi;
            int localCol = vj * 3 + dofj;

            if (vertexIndices[vi] >= 0 && vertexIndices[vj] >= 0) {
              int globalRow = vertexIndices[vi] * 3 + dofi;
              int globalCol = vertexIndices[vj] * 3 + dofj;

              idxM(localRow, localCol) = ES::findEntryOffset(eimpl->KTemplateHandle.hess, globalRow, globalCol);
            }
            else {
              idxM(localRow, localCol) = -1;
            }
          }
        }
      }
    }

    eimpl->inverseIndices[ele] = idxM;
  }
}

DeformationModelAssemblerCacheData *DeformationModelAssemblerElastic::allocateCache() const
{
  DeformationModelAssemblerCacheData *data = new DeformationModelAssemblerCacheData;
  data->internalForceVertexLocks = std::vector<tbb::spin_mutex>(n3 / 3);
  data->stiffnessMatrixVertexRowLocks = std::vector<tbb::spin_mutex>(n3 / 3);

  data->elementCacheData.assign(nele, nullptr);
  for (int i = 0; i < nele; i++) {
    data->elementCacheData[i] = eimpl->femModels[i]->allocateCacheData();
  }

  return data;
}

void DeformationModelAssemblerElastic::freeCache(DeformationModelAssemblerCacheData *data) const
{
  for (int i = 0; i < nele; i++) {
    eimpl->femModels[i]->freeCacheData(data->elementCacheData[i]);
  }
  delete data;
}

void DeformationModelAssemblerElastic::getFixedPlasticParameters(int ele, double *param) const
{
  if (numPlasticParams == 0)
    return;

  if (enableFullspaceScaling) {
    for (int dof = 0; dof < numPlasticParams; dof++) {
      param[dof] = eimpl->fixedParameters[ele * numPlasticParams + dof];
    }
  }
  else {
    if (numHandles == 0)
      return;

    Eigen::Map<ES::VXd>(param, numPlasticParams) = eimpl->elementWeightMatricesPlastic[ele] * eimpl->fixedParameters.head(eimpl->elementWeightMatricesPlastic[ele].cols());
    eimpl->femModels[ele]->getPlasticModel()->compute_paramfull(param);
  }
}

void DeformationModelAssemblerElastic::getFixedElasticParameters(int ele, double *param) const
{
  if (numElasticParams == 0)
    return;

  int offset = nele * numPlasticParams;

  if (enableFullspaceScaling) {
    for (int dof = 0; dof < numElasticParams; dof++) {
      param[dof] = eimpl->fixedParameters[offset + ele * numElasticParams + dof];
    }
  }
  else {
    if (numHandles == 0)
      return;

    Eigen::Map<ES::VXd>(param, numElasticParams) = eimpl->elementWeightMatricesElastic[ele] * eimpl->fixedParameters.segment(offset, eimpl->elementWeightMatricesElastic[ele].cols());
  }
}

double DeformationModelAssemblerElastic::computeEnergy(const double *x, DeformationModelAssemblerCacheData *data) const
{
  for (auto it = data->energyLocalBuffer.begin(); it != data->energyLocalBuffer.end(); ++it)
    *it = 0.0;

  auto localEnergyFunc = [this, x, data](int ele) {
    if (eimpl->elementFlags[ele] == 0)
      return;

    ES::V18d localp;
    for (int j = 0; j < neleVtx; j++) {
      int vid = deformationModelManager->getMesh()->getVertexIndex(ele, j);
      if (vid >= 0)
        localp.segment<3>(j * 3) = ES::V3d(x[vid * 3], x[vid * 3 + 1], x[vid * 3 + 2]);
      else
        localp.segment<3>(j * 3).setZero();
    }

    ES::V12d plasticParam = ES::V12d::Zero();
    double elasticParam[100];

    getFixedPlasticParameters(ele, plasticParam.data());
    getFixedElasticParameters(ele, elasticParam);

    const DeformationModel *fem = eimpl->femModels[ele];
    fem->prepareData(localp.data(), plasticParam.data(), elasticParam, data->elementCacheData[ele]);
    double energy = fem->computeEnergy(data->elementCacheData[ele]);
    /*
    if (energy < 0) {
      LGI << "Encounter energy " << energy << " < 0\n"
          << "ele: " << ele << ',' << eleType << '\n'
          << "p:\n"
          << localp << '\n'
          << "plastic:\n"
          << plasticParam << '\n'
          << "act: " << activationParam;
    }*/
    data->energyLocalBuffer.local() += energy * eimpl->elementFlags[ele];
  };

  // for (int ele = 0; ele < nele; ele++)
  //   localEnergyFunc(ele);
  tbb::parallel_for(0, nele, localEnergyFunc, data->partitioners[0]);

  double energyAll = 0;
  for (auto it = data->energyLocalBuffer.begin(); it != data->energyLocalBuffer.end(); ++it)
    energyAll += *it;

  return energyAll;
}

void DeformationModelAssemblerElastic::computeVonMisesStresses(const double *x, DeformationModelAssemblerCacheData *data, double *elementStresses) const
{
  auto localStressFunc = [this, x, data, elementStresses](int ele) {
    if (eimpl->elementFlags[ele] == 0)
      return;

    ES::V18d localp;
    for (int j = 0; j < neleVtx; j++) {
      int vid = deformationModelManager->getMesh()->getVertexIndex(ele, j);
      if (vid >= 0)
        localp.segment<3>(j * 3) = ES::V3d(x[vid * 3], x[vid * 3 + 1], x[vid * 3 + 2]);
      else
        localp.segment<3>(j * 3).setZero();
    }

    ES::V12d plasticParam = ES::V12d::Zero();
    double elasticParam[100];

    getFixedPlasticParameters(ele, plasticParam.data());
    getFixedElasticParameters(ele, elasticParam);

    const DeformationModel *fem = eimpl->femModels[ele];
    fem->prepareData(localp.data(), plasticParam.data(), elasticParam, data->elementCacheData[ele]);

    double stress[12];
    int nPt;
    fem->vonMisesStress(data->elementCacheData[ele], nPt, stress);

    for (int j = 0; j < nPt; j++)
      elementStresses[ele * nPt + j] = stress[j];
  };

  // for (int ele = 0; ele < nele; ele++)
  //   localEnergyFunc(ele);
  tbb::parallel_for(0, nele, localStressFunc);
}

void DeformationModelAssemblerElastic::computeMaxStrains(const double *x, DeformationModelAssemblerCacheData *data, double *elementStrains) const
{
  auto localStrainFunc = [this, x, data, elementStrains](int ele) {
    if (eimpl->elementFlags[ele] == 0)
      return;

    ES::V18d localp;
    for (int j = 0; j < neleVtx; j++) {
      int vid = deformationModelManager->getMesh()->getVertexIndex(ele, j);
      if (vid >= 0)
        localp.segment<3>(j * 3) = ES::V3d(x[vid * 3], x[vid * 3 + 1], x[vid * 3 + 2]);
      else
        localp.segment<3>(j * 3).setZero();
    }

    ES::V12d plasticParam = ES::V12d::Zero();
    double elasticParam[100];

    getFixedPlasticParameters(ele, plasticParam.data());
    getFixedElasticParameters(ele, elasticParam);

    const DeformationModel *fem = eimpl->femModels[ele];
    fem->prepareData(localp.data(), plasticParam.data(), elasticParam, data->elementCacheData[ele]);

    double stress[12];
    int nPt;
    fem->maxStrain(data->elementCacheData[ele], nPt, stress);

    for (int j = 0; j < nPt; j++)
      elementStrains[ele * nPt + j] = stress[j];
  };

  // for (int ele = 0; ele < nele; ele++)
  //   localEnergyFunc(ele);
  tbb::parallel_for(0, nele, localStrainFunc);
}

void DeformationModelAssemblerElastic::computeGradient(const double *x, double *grad, DeformationModelAssemblerCacheData *data) const
{
  memset(grad, 0, sizeof(double) * n3);
  auto localGradFunc = [this, x, grad, data](int ele) {
    if (eimpl->elementFlags[ele] == 0)
      return;

    ES::V18d localp;
    for (int j = 0; j < neleVtx; j++) {
      int vid = deformationModelManager->getMesh()->getVertexIndex(ele, j);

      if (vid >= 0)
        localp.segment<3>(j * 3) = ES::V3d(x[vid * 3], x[vid * 3 + 1], x[vid * 3 + 2]);
      else
        localp.segment<3>(j * 3).setZero();
    }

    ES::V12d plasticParam = ES::V12d::Zero();
    double elasticParam[100];

    getFixedPlasticParameters(ele, plasticParam.data());
    getFixedElasticParameters(ele, elasticParam);

    // std::cout << ele << " Fp:";
    // for (int i = 0; i < numPlasticParams; i++) {
    //   std::cout << plasticParam[i] << ',';
    // }
    // std::cout << std::endl;

    const DeformationModel *fem = eimpl->femModels[ele];
    fem->prepareData(localp.data(), plasticParam.data(), elasticParam, data->elementCacheData[ele]);

    ES::V18d localGradx;
    fem->compute_dE_dx(data->elementCacheData[ele], localGradx.data());
    localGradx *= eimpl->elementFlags[ele];

    if (enableSanityCheck) {
      for (int i = 0; i < neleVtx * 3; i++) {
        if (std::isfinite(localGradx[i]) == false) {
          SPDLOG_LOGGER_ERROR(Logging::lgr(), "Encounter weird numbers.\nGrad:\n{}\n;x:{}\n", localGradx, localp);
          SPDLOG_LOGGER_ERROR(Logging::lgr(), "Plastic param: {}, {}, {}\n", plasticParam[0], plasticParam[1], plasticParam[2]);
          exit(EXIT_FAILURE);
        }
      }
    }

    const int *vertexIndices = deformationModelManager->getMesh()->getVertexIndices(ele);
    for (int v = 0; v < neleVtx; v++) {
      if (vertexIndices[v] >= 0) {
        data->internalForceVertexLocks[vertexIndices[v]].lock();

        grad[vertexIndices[v] * 3] += localGradx[v * 3];
        grad[vertexIndices[v] * 3 + 1] += localGradx[v * 3 + 1];
        grad[vertexIndices[v] * 3 + 2] += localGradx[v * 3 + 2];

        data->internalForceVertexLocks[vertexIndices[v]].unlock();
      }
    }
  };

  // for (int ele = 0; ele < nele; ele++) {
  //   localGradFunc(ele);
  // }
  tbb::parallel_for(0, nele, localGradFunc, data->partitioners[1]);

  // clean the numbers
  if (enableSanityCheck) {
    for (int i = 0; i < n3; i++) {
      int fpclass = std::fpclassify(grad[i]);
      if (fpclass == FP_INFINITE || fpclass == FP_NAN) {
        SPDLOG_LOGGER_ERROR(Logging::lgr(), "Encounter weird numbers at {}: {}", i, grad[i]);
        exit(EXIT_FAILURE);
      }
      else if (fpclass == FP_SUBNORMAL) {
        grad[i] = 0;
      }
    }
  }
}

void DeformationModelAssemblerElastic::computeHessian(const double *x, HessianMatrixHandle *hessHandle, DeformationModelAssemblerCacheData *data) const
{
  memset(hessHandle->hess.valuePtr(), 0, sizeof(double) * hessHandle->hess.nonZeros());

  auto localHessFunc = [this, x, &hessHandle, data](int ele) {
    if (eimpl->elementFlags[ele] == 0)
      return;

    ES::V18d localp;
    for (int j = 0; j < neleVtx; j++) {
      int vid = deformationModelManager->getMesh()->getVertexIndex(ele, j);

      if (vid >= 0)
        localp.segment<3>(j * 3) = ES::V3d(x[vid * 3], x[vid * 3 + 1], x[vid * 3 + 2]);
      else
        localp.segment<3>(j * 3).setZero();
    }

    ES::V12d plasticParam = ES::V12d::Zero();
    double elasticParam[100];

    getFixedPlasticParameters(ele, plasticParam.data());
    getFixedElasticParameters(ele, elasticParam);

    const DeformationModel *fem = eimpl->femModels[ele];
    fem->prepareData(localp.data(), plasticParam.data(), elasticParam, data->elementCacheData[ele]);

    double localKData[24 * 24];
    fem->compute_d2E_dx2(data->elementCacheData[ele], localKData);

    ES::Mp<ES::MXd> localK(localKData, neleVtx * 3, neleVtx * 3);
    localK *= eimpl->elementFlags[ele];

    const int *vertexIndices = deformationModelManager->getMesh()->getVertexIndices(ele);
    const auto &idxM = eimpl->inverseIndices[ele];

    // write matrices in place
    for (int vb = 0; vb < neleVtx; vb++) {
      int vIdxB = vertexIndices[vb];
      if (vIdxB < 0)
        continue;

      data->stiffnessMatrixVertexRowLocks[vIdxB].lock();

      for (int va = 0; va < neleVtx; va++) {
        int vIdxA = vertexIndices[va];
        if (vIdxA < 0)
          continue;

        for (int j = 0; j < 3; j++) {
          for (int i = 0; i < 3; i++) {
            int local_col = 3 * vb + j;
            int local_row = 3 * va + i;

            std::ptrdiff_t offset = idxM(local_row, local_col);
            if (offset >= 0)
              hessHandle->hess.valuePtr()[offset] += localK(local_row, local_col);
          }  // i
        }    // j
      }      // vb

      data->stiffnessMatrixVertexRowLocks[vIdxB].unlock();
    }  // va
  };

  // for (int ele = 0; ele < nele; ele++) {
  //   localHessFunc(ele);
  // }

  tbb::parallel_for(0, nele, localHessFunc, data->partitioners[2]);

  // clean the numbers
  if (enableSanityCheck) {
    for (Eigen::Index i = 0; i < hessHandle->hess.nonZeros(); i++) {
      int fpclass = std::fpclassify(hessHandle->hess.valuePtr()[i]);
      if (fpclass == FP_INFINITE || fpclass == FP_NAN) {
        SPDLOG_LOGGER_ERROR(Logging::lgr(), "Encounter weird numbers at {}: {}", i, hessHandle->hess.valuePtr()[i]);
        throw std::logic_error("Encounter weird numbers.");
      }
      else if (fpclass == FP_SUBNORMAL) {
        hessHandle->hess.valuePtr()[i] = 0;
      }
    }
  }
}
