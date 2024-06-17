/*
author: Bohan Wang
copyright to USC,MIT
*/

#include "deformationModelAssemblerElasticMultiPlastic.h"
#include "deformationModelAssemblerElasticMultiPlasticImpl.h"
#include "deformationModel.h"
#include "deformationModelManager.h"
#include "simulationMesh.h"
#include "plasticModel.h"

#include "tetMeshManifold.h"
#include "tetMesh.h"
#include "pgoLogging.h"
#include "fmtEigen.h"

#include <tbb/parallel_for.h>
#include <tbb/task_arena.h>

#include <vector>

using namespace pgo;
using namespace pgo::SolidDeformationModel;

DeformationModelAssemblerElasticMultiPlastic::DeformationModelAssemblerElasticMultiPlastic(const DeformationModelManager *dm,
  const double *basisMatrix, int nh, const double *elementFlags):
  DeformationModelAssembler(dm, DeformationModelAssemblingType::AT_ELASTIC_MULTI_PLASTIC, elementFlags),
  numHandles(nh)
{
  PGO_ALOG(basisMatrix != nullptr && numHandles != 0);

  SPDLOG_LOGGER_INFO(Logging::lgr(), "Assembler parameter:{},{},{}", numElasticParams, numPlasticParams, numHandles);

  epimpl = new DeformationModelAssemblerElasticMultiPlasticImpl(*impl);
  delete impl;
  impl = epimpl;

  numPlasticParams = deformationModelManager->getNumPlasticParameters();
  numElasticParams = deformationModelManager->getNumElasticParameters();

  ES::MXd W = Eigen::Map<const ES::MXd>(basisMatrix, nele, numHandles);
  epimpl->elementWeightMatricesPlastic.assign(nele, ES::MXd::Zero(numPlasticParams, numPlasticParams * numHandles));
  epimpl->elementWeightMatricesCompressed.assign(nele, ES::VXd::Zero(numHandles));

  if (numElasticParams > 0)
    epimpl->elementWeightMatricesElastic.resize(nele, ES::MXd::Zero(numElasticParams, numElasticParams * numHandles));

  tbb::parallel_for(
    0, nele, [&](int ele) {
      // for (int ele = 0; ele < nele; ele++) {
      for (int hi = 0; hi < numHandles; hi++) {
        epimpl->elementWeightMatricesPlastic[ele].block(0, hi * numPlasticParams, numPlasticParams, numPlasticParams) =
          ES::MXd::Identity(numPlasticParams, numPlasticParams) * W(ele, hi);

        if (numElasticParams > 0)
          epimpl->elementWeightMatricesElastic[ele].block(0, hi * numElasticParams, numElasticParams, numElasticParams) =
            ES::MXd::Identity(numElasticParams, numElasticParams) * W(ele, hi);
      }
      epimpl->elementWeightMatricesCompressed[ele] = W.row(ele);
    },
    tbb::static_partitioner());

  std::vector<ES::TripletD> entries;
  tbb::enumerable_thread_specific<std::vector<ES::TripletD>> entriesTLS;

  // for (int ele = 0; ele < nele; ele++) {
  tbb::parallel_for(
    0, nele, [&](int ele) {
      const int *vertexIndices = deformationModelManager->getMesh()->getVertexIndices(ele);
      auto &entriesLocal = entriesTLS.local();

      // add element to the stiffness matrix
      for (int vi = 0; vi < 4; vi++) {
        for (int vj = 0; vj < 4; vj++) {
          for (int dofi = 0; dofi < 3; dofi++) {
            for (int dofj = 0; dofj < 3; dofj++) {
              entriesLocal.emplace_back(vertexIndices[vi] * 3 + dofi, vertexIndices[vj] * 3 + dofj, 1.0);
            }
          }
        }
      }

      // upper-right block
      for (int vi = 0; vi < 4; vi++) {
        for (int dofi = 0; dofi < 3; dofi++) {
          for (int parami = 0; parami < numPlasticParams * numHandles; parami++) {
            entriesLocal.emplace_back(vertexIndices[vi] * 3 + dofi, n3 + parami, 1.0);
          }
        }
      }

      // lower-left block
      for (int parami = 0; parami < numPlasticParams * numHandles; parami++) {
        for (int vi = 0; vi < 4; vi++) {
          for (int dofi = 0; dofi < 3; dofi++) {
            entriesLocal.emplace_back(n3 + parami, vertexIndices[vi] * 3 + dofi, 1.0);
          }
        }
      }
    },
    tbb::static_partitioner());

  for (int parami = 0; parami < numPlasticParams * numHandles; parami++) {
    for (int paramj = 0; paramj < numPlasticParams * numHandles; paramj++) {
      entries.emplace_back(n3 + parami, n3 + paramj, 1.0);
    }
  }

  for (auto it = entriesTLS.begin(); it != entriesTLS.end(); ++it) {
    entries.insert(entries.end(), it->begin(), it->end());
  }

  epimpl->KTemplateHandle.hess.resize(n3 + numPlasticParams * numHandles, n3 + numPlasticParams * numHandles);
  epimpl->KTemplateHandle.hess.setFromTriplets(entries.begin(), entries.end());

  epimpl->inverseIndices.resize(nele);
  // for (int ele = 0; ele < nele; ele++) {
  tbb::parallel_for(
    0, nele, [&](int ele) {
      const int *vertexIndices = deformationModelManager->getMesh()->getVertexIndices(ele);
      DeformationModelAssemblerImpl::IndexMatrix idxM = DeformationModelAssemblerImpl::IndexMatrix::Ones() * -1;

      // upper-left block
      for (int vi = 0; vi < 4; vi++) {
        for (int vj = 0; vj < 4; vj++) {
          for (int dofi = 0; dofi < 3; dofi++) {
            for (int dofj = 0; dofj < 3; dofj++) {
              int globalRow = vertexIndices[vi] * 3 + dofi;
              int globalCol = vertexIndices[vj] * 3 + dofj;

              int localRow = vi * 3 + dofi;
              int localCol = vj * 3 + dofj;

              idxM(localRow, localCol) = ES::findEntryOffset(epimpl->KTemplateHandle.hess, globalRow, globalCol);
            }
          }
        }
      }
      epimpl->inverseIndices[ele] = idxM;
    },
    tbb::static_partitioner());

  epimpl->upperRightBlockInverseIndices = DeformationModelAssemblerImpl::DynamicIndexMatrix::Ones(n3, numPlasticParams * numHandles) * -1;
  // for (ES::IDX outeri = 0; outeri < n3; outeri++) {
  tbb::parallel_for((ES::IDX)0, (ES::IDX)n3,
    [&](ES::IDX outeri) {
      const ES::SpMatD::StorageIndex *rowBegin = epimpl->KTemplateHandle.hess.innerIndexPtr() + epimpl->KTemplateHandle.hess.outerIndexPtr()[outeri];
      const ES::SpMatD::StorageIndex *rowEnd = epimpl->KTemplateHandle.hess.innerIndexPtr() + epimpl->KTemplateHandle.hess.outerIndexPtr()[outeri + 1];

      int i = 0;
      while (rowBegin + i < rowEnd && rowBegin[i] < n3)
        i++;

      while (rowBegin + i < rowEnd) {
        epimpl->upperRightBlockInverseIndices(outeri, rowBegin[i] - n3) = rowBegin + i - epimpl->KTemplateHandle.hess.innerIndexPtr();
        i++;
      }
    },
    tbb::static_partitioner());

  epimpl->lowerLeftBlockInverseIndices = DeformationModelAssemblerImpl::DynamicIndexMatrix::Ones(numPlasticParams * numHandles, n3) * -1;
  for (ES::IDX outeri = n3; outeri < epimpl->KTemplateHandle.hess.outerSize(); outeri++) {
    const ES::SpMatD::StorageIndex *rowBegin = epimpl->KTemplateHandle.hess.innerIndexPtr() + epimpl->KTemplateHandle.hess.outerIndexPtr()[outeri];
    const ES::SpMatD::StorageIndex *rowEnd = epimpl->KTemplateHandle.hess.innerIndexPtr() + epimpl->KTemplateHandle.hess.outerIndexPtr()[outeri + 1];

    for (int i = 0; rowBegin + i < rowEnd && rowBegin[i] < n3; i++) {
      epimpl->lowerLeftBlockInverseIndices(outeri - n3, rowBegin[i]) = rowBegin + i - epimpl->KTemplateHandle.hess.innerIndexPtr();
    }
  }

  epimpl->lowerRightBlockInverseIndices = DeformationModelAssemblerImpl::DynamicIndexMatrix::Ones(numPlasticParams * numHandles, numPlasticParams * numHandles) * -1;
  for (ES::IDX outeri = n3; outeri < epimpl->KTemplateHandle.hess.outerSize(); outeri++) {
    const ES::SpMatD::StorageIndex *rowBegin = epimpl->KTemplateHandle.hess.innerIndexPtr() + epimpl->KTemplateHandle.hess.outerIndexPtr()[outeri];
    const ES::SpMatD::StorageIndex *rowEnd = epimpl->KTemplateHandle.hess.innerIndexPtr() + epimpl->KTemplateHandle.hess.outerIndexPtr()[outeri + 1];

    int i = 0;
    while (rowBegin + i < rowEnd && rowBegin[i] < n3)
      i++;

    while (rowBegin + i < rowEnd) {
      epimpl->lowerRightBlockInverseIndices(outeri - n3, rowBegin[i] - n3) = rowBegin + i - epimpl->KTemplateHandle.hess.innerIndexPtr();
      i++;
    }
  }
}

DeformationModelAssemblerCacheData *DeformationModelAssemblerElasticMultiPlastic::allocateCache() const
{
  DeformationModelAssemblerElasticMultiPlasticCacheData *cache = new DeformationModelAssemblerElasticMultiPlasticCacheData;

  cache->internalForceVertexLocks = std::vector<tbb::spin_mutex>(n3 / 3);
  cache->stiffnessMatrixVertexRowLocks = std::vector<tbb::spin_mutex>(n3 / 3);

  cache->elementCacheData.assign(nele, nullptr);
  for (int i = 0; i < nele; i++) {
    cache->elementCacheData[i] = epimpl->femModels[i]->allocateCacheData();
  }

  cache->da_dzTLS = std::make_shared<tbb::enumerable_thread_specific<ES::MXd>>(ES::MXd::Zero(numPlasticParams, numHandles * numPlasticParams));
  cache->d2ai_dz2TLS = std::make_shared<tbb::enumerable_thread_specific<ES::MXd>>(ES::MXd::Zero(numHandles * numPlasticParams, numHandles * numPlasticParams));
  cache->d2E_dadzTLS = std::make_shared<tbb::enumerable_thread_specific<ES::MXd>>(ES::MXd::Zero(numPlasticParams, numHandles * numPlasticParams));
  cache->d2E_dudzTLS = std::make_shared<tbb::enumerable_thread_specific<ES::MXd>>(ES::MXd::Zero(12, numHandles * numPlasticParams));

  return cache;
}

void DeformationModelAssemblerElasticMultiPlastic::freeCache(DeformationModelAssemblerCacheData *cache) const
{
  DeformationModelAssemblerElasticMultiPlasticCacheData *data = static_cast<DeformationModelAssemblerElasticMultiPlasticCacheData *>(cache);
  for (int i = 0; i < nele; i++) {
    epimpl->femModels[i]->freeCacheData(data->elementCacheData[i]);
  }
  delete data;
}

void DeformationModelAssemblerElasticMultiPlastic::getFixedElasticParameters(int ele, double *param) const
{
  if (numElasticParams == 0)
    return;
  else {
    SPDLOG_LOGGER_CRITICAL(Logging::lgr(), "Not done yet");
    abort();
    // todo
  }
}

void DeformationModelAssemblerElasticMultiPlastic::getPlasticParameters(const double *x_param, int ele, double *params) const
{
  for (int parami = 0; parami < numPlasticParams; parami++) {
    params[parami] = 0.0;
  }

  (Eigen::Map<ES::VXd>(params, numPlasticParams)) = epimpl->elementWeightMatricesPlastic[ele] * Eigen::Map<const ES::VXd>(x_param + n3, numHandles * numPlasticParams);
  deformationModelManager->getDeformationModel(ele)->getPlasticModel()->compute_paramfull(params);
}

double DeformationModelAssemblerElasticMultiPlastic::computeEnergy(const double *x_param, DeformationModelAssemblerCacheData *cache) const
{
  DeformationModelAssemblerElasticMultiPlasticCacheData *data = static_cast<DeformationModelAssemblerElasticMultiPlasticCacheData *>(cache);

  for (auto it = data->energyLocalBuffer.begin(); it != data->energyLocalBuffer.end(); ++it)
    *it = 0.0;

  auto localEnergyFunc = [&](int ele) {
    if (epimpl->elementFlags[ele] == 0)
      return;

    ES::V12d localp;
    for (int j = 0; j < 4; j++) {
      int vid = deformationModelManager->getMesh()->getVertexIndex(ele, j);
      localp.segment<3>(j * 3) = ES::V3d(x_param[vid * 3], x_param[vid * 3 + 1], x_param[vid * 3 + 2]);
    }

    ES::V12d plasticParam = ES::V12d::Zero();
    ES::V12d elasticParam = ES::V12d::Zero();
    getFixedElasticParameters(ele, elasticParam.data());
    getPlasticParameters(x_param, ele, plasticParam.data());

    const DeformationModel *fem = deformationModelManager->getDeformationModel(ele);
    fem->prepareData(localp.data(), plasticParam.data(), elasticParam.data(), data->elementCacheData[ele]);
    data->energyLocalBuffer.local() += fem->computeEnergy(data->elementCacheData[ele]) * epimpl->elementFlags[ele];
  };

  // for (int ele = 0; ele < deformationModelManager->getMesh()->getNumElements(); ele++) {
  tbb::parallel_for(0, deformationModelManager->getMesh()->getNumElements(), localEnergyFunc, data->partitioners[0]);

  double energyAll = 0;
  for (auto it = data->energyLocalBuffer.begin(); it != data->energyLocalBuffer.end(); ++it)
    energyAll += *it;

  return energyAll;
}

void DeformationModelAssemblerElasticMultiPlastic::computeGradient(const double *x_param, double *grad, DeformationModelAssemblerCacheData *cache) const
{
  DeformationModelAssemblerElasticMultiPlasticCacheData *data = static_cast<DeformationModelAssemblerElasticMultiPlasticCacheData *>(cache);
  memset(grad, 0, sizeof(double) * (n3 + numHandles * numPlasticParams));

  for (auto it = data->gradBlockTLS.begin(); it != data->gradBlockTLS.end(); ++it)
    memset(it->data(), 0, it->size() * sizeof(double));

  auto localGradFunc = [&](int ele) {
    if (epimpl->elementFlags[ele] == 0)
      return;

    ES::V12d localp;
    for (int j = 0; j < 4; j++) {
      int vid = deformationModelManager->getMesh()->getVertexIndex(ele, j);
      localp.segment<3>(j * 3) = ES::V3d(x_param[vid * 3], x_param[vid * 3 + 1], x_param[vid * 3 + 2]);
    }

    ES::V12d plasticParam = ES::V12d::Zero();
    ES::V12d elasticParam = ES::V12d::Zero();
    getFixedElasticParameters(ele, elasticParam.data());
    getPlasticParameters(x_param, ele, plasticParam.data());

    const DeformationModel *fem = deformationModelManager->getDeformationModel(ele);
    fem->prepareData(localp.data(), plasticParam.data(), elasticParam.data(), data->elementCacheData[ele]);

    ES::V12d localGradx;
    fem->compute_dE_dx(data->elementCacheData[ele], localGradx.data());
    localGradx *= epimpl->elementFlags[ele];

    for (int i = 0; i < 12; i++) {
      if (std::isfinite(localGradx[i]) == false) {
        SPDLOG_LOGGER_ERROR(Logging::lgr(), "Encounter weird numbers.\nGrad:\n{}\n;x:\n{}", localGradx, localp);
        exit(EXIT_FAILURE);
      }
    }

    const int *vertexIndices = deformationModelManager->getMesh()->getVertexIndices(ele);
    for (int v = 0; v < 4; v++) {
      data->internalForceVertexLocks[vertexIndices[v]].lock();

      grad[vertexIndices[v] * 3] += localGradx[v * 3];
      grad[vertexIndices[v] * 3 + 1] += localGradx[v * 3 + 1];
      grad[vertexIndices[v] * 3 + 2] += localGradx[v * 3 + 2];

      data->internalForceVertexLocks[vertexIndices[v]].unlock();
    }

    ES::V12d localGrada = ES::V12d::Zero();
    fem->compute_dE_da(data->elementCacheData[ele], localGrada.data());
    localGrada *= epimpl->elementFlags[ele];

    for (int i = 0; i < numPlasticParams; i++) {
      if (std::isfinite(localGrada[i]) == false) {
        SPDLOG_LOGGER_ERROR(Logging::lgr(), "Encounter weird numbers.\nGrad:\n{}\n;x:\n{}", localGrada, localp);
        exit(EXIT_FAILURE);
      }
    }

    ES::VXd &gradBlock = data->gradBlockTLS.local();
    ES::MXd &dadz = data->da_dzTLS->local();

    tbb::this_task_arena::isolate([&]() {
      if (gradBlock.rows() == 0)
        gradBlock = ES::VXd::Zero(numHandles * numPlasticParams);

      fem->getPlasticModel()->compute_dparamfull_dparamsub(x_param + n3, epimpl->elementWeightMatricesCompressed[ele].data(), numHandles, dadz.data());
      gradBlock += dadz.transpose() * localGrada.block(0, 0, numPlasticParams, 1);
    });
  };

  // for (int ele = 0; ele < deformationModelManager->getMesh()->getNumElements(); ele++) {
  tbb::parallel_for(0, nele, localGradFunc, data->partitioners[1]);

  for (auto it = data->gradBlockTLS.begin(); it != data->gradBlockTLS.end(); ++it) {
    Eigen::Map<ES::VXd>(grad + n3, numHandles * numPlasticParams) += *it;
  }

  // clean the numbers
  // static int debugid = 0;
  // std::ofstream outfile(fmt::format("f:/grad-{}.txt", debugid++).c_str());
  for (int i = 0; i < n3 + numHandles * numPlasticParams; i++) {
    int fpclass = std::fpclassify(grad[i]);
    if (fpclass == FP_INFINITE || fpclass == FP_NAN) {
      SPDLOG_LOGGER_ERROR(Logging::lgr(), "Encounter weird numbers at {}: {}", i, grad[i]);
      exit(EXIT_FAILURE);
    }
    else if (fpclass == FP_SUBNORMAL) {
      grad[i] = 0;
    }
    // outfile << grad[i] << ' ' << fpclass << std::endl;
  }
}

void DeformationModelAssemblerElasticMultiPlastic::computeHessian(const double *x_param, HessianMatrixHandle *hessHandle, DeformationModelAssemblerCacheData *cache) const
{
  DeformationModelAssemblerElasticMultiPlasticCacheData *data = static_cast<DeformationModelAssemblerElasticMultiPlasticCacheData *>(cache);

  memset(hessHandle->hess.valuePtr(), 0, sizeof(double) * hessHandle->hess.nonZeros());

  for (auto it = data->upperRightBlockTLS.begin(); it != data->upperRightBlockTLS.end(); ++it)
    memset(it->data(), 0, it->size() * sizeof(double));

  for (auto it = data->lowerRightBlockTLS.begin(); it != data->lowerRightBlockTLS.end(); ++it)
    memset(it->data(), 0, it->size() * sizeof(double));

  auto localHessFunc = [&](int ele) {
    if (epimpl->elementFlags[ele] == 0)
      return;

    Eigen::Matrix<double, 12, 1> localp;
    for (int j = 0; j < 4; j++) {
      int vid = deformationModelManager->getMesh()->getVertexIndex(ele, j);
      localp.segment<3>(j * 3) = ES::V3d(x_param[vid * 3], x_param[vid * 3 + 1], x_param[vid * 3 + 2]);
    }

    ES::V12d plasticParam = ES::V12d::Zero();
    ES::V12d elasticParam = ES::V12d::Zero();
    getFixedElasticParameters(ele, elasticParam.data());
    getPlasticParameters(x_param, ele, plasticParam.data());

    const DeformationModel *fem = deformationModelManager->getDeformationModel(ele);
    fem->prepareData(localp.data(), plasticParam.data(), elasticParam.data(), data->elementCacheData[ele]);

    ES::M12d localK;
    fem->compute_d2E_dx2(data->elementCacheData[ele], localK.data());
    localK *= epimpl->elementFlags[ele];

    const int *vertexIndices = deformationModelManager->getMesh()->getVertexIndices(ele);
    const auto &idxM = epimpl->inverseIndices[ele];

    // write matrices in place
    for (int vb = 0; vb < 4; vb++) {
      int vIdxB = vertexIndices[vb];
      data->stiffnessMatrixVertexRowLocks[vIdxB].lock();

      for (int va = 0; va < 4; va++) {
        // int vIdxA = vertexIndices[va];

        for (int j = 0; j < 3; j++) {
          for (int i = 0; i < 3; i++) {
            int local_col = 3 * vb + j;
            int local_row = 3 * va + i;

            std::ptrdiff_t offset = idxM(local_row, local_col);
            hessHandle->hess.valuePtr()[offset] += localK(local_row, local_col);
          }  // i
        }    // j
      }      // vb

      data->stiffnessMatrixVertexRowLocks[vIdxB].unlock();
    }  // va

    double A[12];
    double B[12 * 12];
    double C[12 * 12];

    fem->compute_dE_da(data->elementCacheData[ele], A);
    fem->compute_d2E_dxda(data->elementCacheData[ele], B);
    fem->compute_d2E_da2(data->elementCacheData[ele], C);

    Eigen::Map<ES::MXd> BMap(B, 12, numPlasticParams);
    Eigen::Map<ES::MXd> CMap(C, numPlasticParams, numPlasticParams);
    BMap *= epimpl->elementFlags[ele];
    CMap *= epimpl->elementFlags[ele];

    ES::MXd &d2E_dudz = data->d2E_dudzTLS->local();
    ES::MXd &d2E_dadz = data->d2E_dadzTLS->local();
    ES::MXd &dadz = data->da_dzTLS->local();
    ES::MXd &d2aidz2 = data->d2ai_dz2TLS->local();
    ES::MXd &localUpperRightBlock = data->upperRightBlockTLS.local();
    ES::MXd &localLowerRightBlock = data->lowerRightBlockTLS.local();

    tbb::this_task_arena::isolate([&] {
      fem->getPlasticModel()->compute_dparamfull_dparamsub(x_param + n3, epimpl->elementWeightMatricesCompressed[ele].data(), numHandles, dadz.data());
      d2E_dudz.noalias() = BMap * dadz;

      if (localUpperRightBlock.rows() == 0)
        localUpperRightBlock = ES::MXd::Zero(n3, numHandles * numPlasticParams);

      for (int j = 0; j < 4; j++) {
        for (int dof = 0; dof < 3; dof++) {
          localUpperRightBlock.row(vertexIndices[j] * 3 + dof) += d2E_dudz.row(j * 3 + dof);
        }
      }

      d2E_dadz.noalias() = CMap * dadz;

      if (localLowerRightBlock.rows() == 0)
        localLowerRightBlock = ES::MXd::Zero(numHandles * numPlasticParams, numHandles * numPlasticParams);

      localLowerRightBlock.noalias() += dadz.transpose() * d2E_dadz;

      if (fem->getPlasticModel()->has2OrderDeriv()) {
        for (int pi = 0; pi < numPlasticParams; pi++) {
          fem->getPlasticModel()->compute_d2paramfull_dparamsub2(x_param + n3, epimpl->elementWeightMatricesCompressed[ele].data(),
            numHandles, pi, d2aidz2.data());
          localLowerRightBlock.noalias() += d2aidz2 * A[pi];
        }
      }
    });
  };

  // for (int ele = 0; ele < deformationModelManager->getMesh()->getNumElements(); ele++) {
  //   localHessFunc(ele);
  // }
  tbb::parallel_for(0, nele, localHessFunc, data->partitioners[2]);

  tbb::parallel_for(
    0, n3, [&](int dofi) {
      // for (int dofi = 0; dofi < n3; dofi++) {
      for (auto it = data->upperRightBlockTLS.begin(); it != data->upperRightBlockTLS.end(); ++it) {
        for (int dofj = 0; dofj < numHandles * numPlasticParams; dofj++) {
          hessHandle->hess.valuePtr()[epimpl->upperRightBlockInverseIndices(dofi, dofj)] += (*it)(dofi, dofj);
          hessHandle->hess.valuePtr()[epimpl->lowerLeftBlockInverseIndices(dofj, dofi)] += (*it)(dofi, dofj);
        }
      }
    },
    tbb::static_partitioner());

  for (auto it = data->lowerRightBlockTLS.begin(); it != data->lowerRightBlockTLS.end(); ++it) {
    for (int col = 0; col < numHandles * numPlasticParams; col++) {
      for (int row = 0; row < numHandles * numPlasticParams; row++) {
        hessHandle->hess.valuePtr()[epimpl->lowerRightBlockInverseIndices(row, col)] += (*it)(row, col);
      }
    }
  }

  // clean the numbers
  for (Eigen::Index i = 0; i < hessHandle->hess.nonZeros(); i++) {
    int fpclass = std::fpclassify(hessHandle->hess.valuePtr()[i]);
    if (fpclass == FP_INFINITE || fpclass == FP_NAN) {
      SPDLOG_LOGGER_ERROR(Logging::lgr(), "Encounter weird numbers at {}: {}", i, hessHandle->hess.valuePtr()[i]);
      exit(EXIT_FAILURE);
    }
    else if (fpclass == FP_SUBNORMAL) {
      hessHandle->hess.valuePtr()[i] = 0;
    }
  }
}

void DeformationModelAssemblerElasticMultiPlastic::createWeightMatrix(const SimulationMesh *tetMesh, const int *elementIDs, int numElementIDs, double *WOut)
{
  // std::vector<Vec3d> vertices;
  // std::vector<Vec4i> tets;
  // tetMesh->exportMeshGeometry(vertices, tets);

#if 0
  std::vector<int> fixedVertices;
  for (int ele : elementIDs) {
    for (int j = 0; j < 4; j++) {
      fixedVertices.push_back(deformationModelManager->getMesh()->getVertexIndex(ele, j));
    }
  }

  ES::MXd W(deformationModelManager->getMesh()->getNumVertices(), fixedVertices.size());
  // my_bbw(deformationModelManager->getMesh(), (int)fixedVertices.size(), fixedVertices.data(), W.data(), "../config.opt");
  // WriteMatrixToDisk("a.dat", deformationModelManager->getMesh()->getNumVertices(), fixedVertices.size(), W.data());

  igl_bbw(deformationModelManager->getMesh()->getNumVertices(), vertices.data(),
    deformationModelManager->getMesh()->getNumElements(), tets.data(),
    (int)fixedVertices.size(), fixedVertices.data(),
    W.data());
  WriteMatrixToDisk("b.dat", deformationModelManager->getMesh()->getNumVertices(), fixedVertices.size(), W.data());

  //LGI << W.row(0).sum();

  ES::MXd Wele(deformationModelManager->getMesh()->getNumElements(), elementIDs.size());
  for (int ele = 0; ele < deformationModelManager->getMesh()->getNumElements(); ele++) {
    Vec3d center(0.0);
    for (int j = 0; j < 4; j++) {
      Vec3d p = deformationModelManager->getMesh()->getVertex(ele, j);
      center += p;
    }
    center *= 0.25;

    double w[4];
    deformationModelManager->getMesh()->computeBarycentricWeights(ele, center, w);

    for (size_t elej = 0; elej < elementIDs.size(); elej++) {
      Wele(ele, elej) = 0;

      for (int j = 0; j < 4; j++) {
        for (int k = 0; k < 4; k++) {
          Wele(ele, elej) += W(deformationModelManager->getMesh()->getVertexIndex(ele, k), elej * 4 + j) * w[k];
        }
      }
    }
  }

  //LGI << Wele.row(0).sum();

  return Wele;
#elif 1
  throw std::runtime_error("todo");

  // std::vector<int> fixedVertices;
  // for (int i = 0; i < numElementIDs; i++)
  //   fixedVertices.push_back(tetMesh->getVertexIndex(elementIDs[i], 0));

  // ES::MXd W(tetMesh->getNumVertices(), fixedVertices.size());
  //// my_bbw(deformationModelManager->getMesh(), (int)fixedVertices.size(), fixedVertices.data(), W.data(), "../config.opt");
  //// WriteMatrixToDisk("a.dat", deformationModelManager->getMesh()->getNumVertices(), fixedVertices.size(), W.data());

  //// libiglInterface::igl_bbw(tetMesh->getNumVertices(), vertices.data(),
  ////   tetMesh->getNumElements(), tets.data(),
  ////   (int)fixedVertices.size(), fixedVertices.data(),
  ////   W.data());
  //// WriteMatrixToDisk("b.dat", tetMesh->getNumVertices(), (int)fixedVertices.size(), W.data());

  //// LGI << W.row(0).sum();

  // ES::MXd Wele(tetMesh->getNumElements(), numElementIDs);
  // for (int ele = 0; ele < tetMesh->getNumElements(); ele++) {
  //   Vec3d center(0.0);
  //   for (int j = 0; j < 4; j++) {
  //     Vec3d p = tetMesh->getVertex(ele, j);
  //     center += p;
  //   }
  //   center *= 0.25;

  //  double w[4];
  //  tetMesh->computeBarycentricWeights(ele, center, w);

  //  for (size_t elej = 0; elej < numElementIDs; elej++) {
  //    Wele(ele, elej) = 0;

  //    for (int j = 0; j < 4; j++) {
  //      Wele(ele, elej) += W(tetMesh->getVertexIndex(ele, j), elej) * w[j];
  //    }
  //  }
  //}

  // for (int ele = 0; ele < tetMesh->getNumElements(); ele++) {
  //   Wele.row(ele) /= Wele.row(ele).sum();
  // }

  //// LGI << Wele.row(0).sum();
  // Eigen::Map<ES::MXd>(WOut, Wele.rows(), Wele.cols()) = Wele;

#else
  ES::MXd Wele(deformationModelManager->getMesh()->getNumElements(), elementIDs.size());
  ES::MXd Ws[4] = {
    ES::MXd(deformationModelManager->getMesh()->getNumVertices(), elementIDs.size()),
    ES::MXd(deformationModelManager->getMesh()->getNumVertices(), elementIDs.size()),
    ES::MXd(deformationModelManager->getMesh()->getNumVertices(), elementIDs.size()),
    ES::MXd(deformationModelManager->getMesh()->getNumVertices(), elementIDs.size())
  };

  for (int j = 0; j < 4; j++) {
    std::vector<int> fixedVertices;
    for (int ele : elementIDs) {
      fixedVertices.push_back(deformationModelManager->getMesh()->getVertexIndex(ele, j));
    }
    ES::MXd W(deformationModelManager->getMesh()->getNumVertices(), fixedVertices.size());

    igl_bbw(deformationModelManager->getMesh()->getNumVertices(), vertices.data(),
      deformationModelManager->getMesh()->getNumElements(), tets.data(),
      (int)fixedVertices.size(), fixedVertices.data(),
      Ws[j].data());
  }

  ES::MXd WAll = Ws[0] + Ws[1] + Ws[2] + Ws[3];

  // WriteMatrixToDisk("b.dat", deformationModelManager->getMesh()->getNumVertices(), fixedVertices.size(), W.data());
  for (int ele = 0; ele < deformationModelManager->getMesh()->getNumElements(); ele++) {
    Vec3d center(0.0);
    for (int j = 0; j < 4; j++) {
      Vec3d p = deformationModelManager->getMesh()->getVertex(ele, j);
      center += p;
    }
    center *= 0.25;

    double w[4];
    deformationModelManager->getMesh()->computeBarycentricWeights(ele, center, w);

    for (size_t elej = 0; elej < elementIDs.size(); elej++) {
      Wele(ele, elej) = 0;

      for (int j = 0; j < 4; j++) {
        Wele(ele, elej) += Ws[j](deformationModelManager->getMesh()->getVertexIndex(ele, j), elej) * w[j];
      }
    }
  }

  for (int ele = 0; ele < deformationModelManager->getMesh()->getNumElements(); ele++) {
    Wele.row(ele) /= Wele.row(ele).sum();
  }

  // my_bbw(deformationModelManager->getMesh(), (int)fixedVertices.size(), fixedVertices.data(), W.data(), "../config.opt");
  // WriteMatrixToDisk("a.dat", deformationModelManager->getMesh()->getNumVertices(), fixedVertices.size(), W.data());

  // LGI << W.row(0).sum();

  // LGI << Wele.row(0).sum();

  return Wele;
#endif
}

#if 0

void DeformationModelAssemblerElasticMultiPlastic::computedf_ds(const double *x_param, double *hess) const
{
  memset(hess, 0, sizeof(double) * (nm * 3 * n3));

  for (auto it = data->upperRightBlockTLS.begin(); it != data->upperRightBlockTLS.end(); ++it)
    memset(it->data(), 0, it->size() * sizeof(double));

  tbb::enumerable_thread_specific<ES::MXd> elementLocalBlockTLS;

  auto localHessFunc = [&](int ele) {
    if (elementFlags[ele] == 0)
      return;

    Eigen::Matrix<double, 12, 1> localp;
    for (int j = 0; j < 4; j++) {
      int vid = deformationModelManager->getMesh()->getVertexIndex(ele, j);
      localp.segment<3>(j * 3) = ES::V3d(x_param[vid * 3], x_param[vid * 3 + 1], x_param[vid * 3 + 2]);
    }

    int eleType = muscleModel->getElementTypes()[ele];
    if (eleType >= 0) {
      ES::V3d plasticParam = ES::V3d::Ones();
      double activationParam = 0;
      const int *vertexIndices = deformationModelManager->getMesh()->getVertexIndices(ele);

      getPlasticParameters(x_param, ele, plasticParam.data());
      activationParam = getFixedActivationParamers(eleType, ele);

      const MuscleElasticMaterialModelWithParameterForTetMesh *fem = muscleModel->getMuscleFEM(ele);
      fem->prepareData(localp.data(), plasticParam.data(), &activationParam, data->elementCacheData[ele]);

      Eigen::Matrix<double, 12, 3> B;
      fem->compute_d2E_duda(data->elementCacheData[ele], B.data());

      ES::MXd &elementLocalBlock = elementLocalBlockTLS.local();
      ES::MXd &localUpperRightBlock = data->upperRightBlockTLS.local();

      tbb::this_task_arena::isolate([&] {
        if (elementLocalBlock.rows() == 0)
          elementLocalBlock = ES::MXd::Zero(12, numHandles * numPlasticParams);
        elementLocalBlock.noalias() = B * elementWeightMatrices[ele];

        if (localUpperRightBlock.rows() == 0)
          localUpperRightBlock = ES::MXd::Zero(n3, numHandles * numPlasticParams);

        for (int j = 0; j < 4; j++) {
          for (int dof = 0; dof < 3; dof++) {
            localUpperRightBlock.row(vertexIndices[j] * 3 + dof) += elementLocalBlock.row(j * 3 + dof);
          }
        }
      });
    }
  };

  //for (int ele = 0; ele < deformationModelManager->getMesh()->getNumElements(); ele++) {,
  tbb::parallel_for(0, nele, localHessFunc, data->partitioners[2]);

  for (auto it = data->upperRightBlockTLS.begin(); it != data->upperRightBlockTLS.end(); ++it) {
    Eigen::Map<ES::MXd>(hess, n3, numHandles * numPlasticParams) += *it;
  }

  // clean the numbers
  for (Eigen::Index i = 0; i < n3 * numHandles * numPlasticParams; i++) {
    int fpclass = std::fpclassify(hess[i]);
    if (fpclass == FP_INFINITE || fpclass == FP_NAN) {
      LGI << "Encounter weird numbers.";
      exit(EXIT_FAILURE);
    }
    else if (fpclass == FP_SUBNORMAL) {
      hess[i] = 0;
    }
  }
}



struct MatrixBuffer
{
  MatrixBuffer(int r)
  {
    d3EdudZdu.resize(12 * r * 12);
    d3EdudZda.resize(12 * r * 3);
    d3EdudZdZ.resize(12 * r * r);
  }

  std::vector<double> d3EdudZdu, d3EdudZda, d3EdudZdZ;
};

template<typename T>
inline T &ijk(T *dataPtr, int i, int j, int k, int n1, int n2)
{
  return dataPtr[k * n1 * n2 + j * n1 + k];
}

void DeformationModelAssemblerElasticMultiPlastic::compute3rdOrderTensor(const double *x_param, const double *lambda, ES::SpMatD &hess) const
{
  memset(hess.valuePtr(), 0, sizeof(double) * hess.nonZeros());

  for (auto it = data->upperRightBlockTLS.begin(); it != data->upperRightBlockTLS.end(); ++it)
    memset(it->data(), 0, it->size() * sizeof(double));

  for (auto it = data->lowerRightBlockTLS.begin(); it != data->lowerRightBlockTLS.end(); ++it)
    memset(it->data(), 0, it->size() * sizeof(double));

  int r = numHandles * numPlasticParams;
  MatrixBuffer exampler(r);
  tbb::enumerable_thread_specific<MatrixBuffer> elementMatrixBufferTLS(exampler);

  auto localHessFunc = [&](int ele) {
    if (elementFlags[ele] == 0)
      return;

    Eigen::Matrix<double, 12, 1> localp;
    for (int j = 0; j < 4; j++) {
      int vid = deformationModelManager->getMesh()->getVertexIndex(ele, j);
      localp.segment<3>(j * 3) = ES::V3d(x_param[vid * 3], x_param[vid * 3 + 1], x_param[vid * 3 + 2]);
    }

    int eleType = muscleModel->getElementTypes()[ele];
    ES::V3d plasticParam = ES::V3d::Ones();
    double activationParam = 0;
    if (eleType >= 0) {
      getPlasticParameters(x_param, ele, plasticParam.data());
      activationParam = getFixedActivationParamers(eleType, ele);
    }

    const MuscleElasticMaterialModelWithParameterForTetMesh *fem = muscleModel->getMuscleFEM(ele);
    fem->prepareData(localp.data(), plasticParam.data(), &activationParam, data->elementCacheData[ele]);

    double tensor[144 * 12];
    fem->compute_d3E_du3(data->elementCacheData[ele], tensor);

    const int *vertexIndices = deformationModelManager->getMesh()->getVertexIndices(ele);
    const auto &idxM = inverseIndices[ele];

    for (int j = 0; j < 12; j++) {
      int vidj = j / 3;
      int dofj = j % 3;

      for (int k = 0; k < 12; k++) {
        int vidk = k / 3;
        int dofk = k % 3;

        double val = 0;
        for (int i = 0; i < 12; i++) {
          int vidi = i / 3;
          int dofi = i % 3;
          val += lambda[vertexIndices[vidi] + dofi] * tensor[i * 144 + j * 12 + k];
        }
        std::ptrdiff_t offset = idxM(j, k);

        data->stiffnessMatrixVertexRowLocks[vertexIndices[vidj]].lock();

        hess.valuePtr()[offset] += val;
        data->stiffnessMatrixVertexRowLocks[vertexIndices[vidj]].unlock();
      }
    }

    if (eleType >= 0) {
      double tensor0[36 * 12];
      fem->compute_d3E_dudadu(data->elementCacheData[ele], tensor0);

      double tensor1[36 * 3];
      fem->compute_d3E_dudada(data->elementCacheData[ele], tensor1);

      auto &matrixBuffer = elementMatrixBufferTLS.local();
      ES::MXd &localUpperRightBlock = data->upperRightBlockTLS.local();
      ES::MXd &localLowerRightBlock = data->lowerRightBlockTLS.local();

      tbb::this_task_arena::isolate([&] {
        if (localUpperRightBlock.rows() == 0)
          localUpperRightBlock = ES::MXd::Zero(n3, r);

        if (localLowerRightBlock.rows() == 0)
          localLowerRightBlock = ES::MXd::Zero(r, r);

        // compute d3E / (du dZ du)
        for (int i = 0; i < 12; i++) {
          Eigen::Map<ES::MXd>(matrixBuffer.d3EdudZdu.data() + i * 36, 12, r) =
            Eigen::Map<const Eigen::Matrix<double, 12, 3>>(tensor0 + i * 36) * elementWeightMatrices[ele];
        }
        // compute d3E / (du dZ da)
        for (int i = 0; i < 3; i++) {
          Eigen::Map<ES::MXd>(matrixBuffer.d3EdudZda.data() + i * 36, 12, r) =
            Eigen::Map<const Eigen::Matrix<double, 12, 3>>(tensor1 + i * 36) * elementWeightMatrices[ele];
        }
        // compute d3E / (du dZ dZ)
        for (int i = 0; i < r; i++) {
          Eigen::Map<ES::MXd> d3EdudZdZi(matrixBuffer.d3EdudZdZ.data() + i * r * 12, 12, r);
          d3EdudZdZi.noalias() = ES::MXd::Zero(12, r);

          for (int j = 0; j < 3; j++) {
            d3EdudZdZi += Eigen::Map<ES::MXd>(matrixBuffer.d3EdudZda.data() + j * 36, 12, r) * elementWeightMatrices[ele](j, r);
          }
        }

        for (int j = 0; j < 12; j++) {
          int vidj = j / 3;
          int dofj = j % 3;

          for (int k = 0; k < r; k++) {
            double val = 0;
            for (int i = 0; i < 12; i++) {
              int vidi = i / 3;
              int dofi = i % 3;
              val += lambda[vertexIndices[vidi] + dofi] * matrixBuffer.d3EdudZdu[j * 12 * r + k * 12 + i];
            }

            localUpperRightBlock(vertexIndices[vidj] + dofj, k) += val;
          }
        }

        for (int j = 0; j < r; j++) {
          for (int k = 0; k < r; k++) {
            for (int i = 0; i < 12; i++) {
              int vidi = i / 3;
              int dofi = i % 3;
              localLowerRightBlock(j, k) += lambda[vertexIndices[vidi] * 3 + dofi] * matrixBuffer.d3EdudZdZ[k * 12 * r + j * 12 + i];
            }
          }
        }
      });
    }
  };

  //for (int ele = 0; ele < deformationModelManager->getMesh()->getNumElements(); ele++) {
  //  localHessFunc(ele);
  //}
  tbb::parallel_for(0, nele, localHessFunc, data->partitioners[3]);

  tbb::parallel_for(
    0, n3, [&](int dofi) {
    // for (int dofi = 0; dofi < n3; dofi++) {
    for (auto it = data->upperRightBlockTLS.begin(); it != data->upperRightBlockTLS.end(); ++it) {
      for (int dofj = 0; dofj < r; dofj++) {
        hess.valuePtr()[upperRightBlockInverseIndices(dofi, dofj)] += (*it)(dofi, dofj);
        hess.valuePtr()[lowerLeftBlockInverseIndices(dofj, dofi)] += (*it)(dofi, dofj);
      }
    }
  },
    tbb::static_partitioner());

  for (auto it = data->lowerRightBlockTLS.begin(); it != data->lowerRightBlockTLS.end(); ++it) {
    for (int col = 0; col < r; col++) {
      for (int row = 0; row < r; row++) {
        hess.valuePtr()[lowerRightBlockInverseIndices(row, col)] += (*it)(row, col);
      }
    }
  }

  // clean the numbers
  for (Eigen::Index i = 0; i < hess.nonZeros(); i++) {
    int fpclass = std::fpclassify(hess.valuePtr()[i]);
    if (fpclass == FP_INFINITE || fpclass == FP_NAN) {
      LGI << "Encounter weird numbers.";
      exit(EXIT_FAILURE);
    }
    else if (fpclass == FP_SUBNORMAL) {
      hess.valuePtr()[i] = 0;
    }
  }
}

#endif