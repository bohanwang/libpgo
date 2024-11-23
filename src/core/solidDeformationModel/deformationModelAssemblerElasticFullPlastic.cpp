/*
author: Bohan Wang
copyright to USC,MIT
*/

#include "deformationModelAssemblerElasticFullPlastic.h"
#include "deformationModelAssemblerElasticFullPlasticImpl.h"
#include "deformationModel.h"
#include "deformationModelManager.h"
#include "simulationMesh.h"

#include "pgoLogging.h"
#include "fmtEigen.h"

#include <tbb/parallel_for.h>

using namespace pgo;
using namespace pgo::SolidDeformationModel;

DeformationModelAssemblerElasticFullPlastic::DeformationModelAssemblerElasticFullPlastic(const DeformationModelManager *dm,
  const double *basisMatrix, int nh, const double *elementFlags):
  DeformationModelAssembler(dm, DeformationModelAssemblingType::AT_ELASTIC_FULL_PLASTIC, elementFlags),
  numHandles(nh)
{
  std::vector<ES::TripletD> entries;
  tbb::enumerable_thread_specific<std::vector<ES::TripletD>> entriesTLS;

  numPlasticParams = deformationModelManager->getNumPlasticParameters();
  numElasticParams = deformationModelManager->getNumElasticParameters();

  SPDLOG_LOGGER_INFO(Logging::lgr(), "Assembler parameter:{},{},{}", numElasticParams, numPlasticParams, numHandles);

  epimpl = new DeformationModelAssemblerElasticFullPlasticImpl(*impl);
  delete impl;
  impl = epimpl;

  if (basisMatrix && numHandles) {
    ES::MXd W = Eigen::Map<const ES::MXd>(basisMatrix, nele, numHandles);
    epimpl->elementWeightMatricesCompressed.assign(nele, ES::VXd::Zero(numHandles));

    if (numElasticParams > 0)
      epimpl->elementWeightMatricesElastic.resize(nele, ES::MXd::Zero(numElasticParams, numElasticParams * numHandles));

    for (int ele = 0; ele < nele; ele++) {
      for (int hi = 0; hi < numHandles; hi++) {
        if (numElasticParams > 0)
          epimpl->elementWeightMatricesElastic[ele].block(0, hi * numElasticParams, numElasticParams, numElasticParams) =
            ES::MXd::Identity(numElasticParams, numElasticParams) * W(ele, hi);
      }

      epimpl->elementWeightMatricesCompressed[ele] = W.row(ele);
    }
  }

  // for (int ele = 0; ele < nele; ele++) {
  tbb::parallel_for(
    0, nele, [&](int ele) {
      const int *vertexIndices = deformationModelManager->getMesh()->getVertexIndices(ele);
      auto &entriesLocal = entriesTLS.local();

      // add element to the stiffness matrix
      for (int vi = 0; vi < neleVtx; vi++) {
        for (int vj = 0; vj < neleVtx; vj++) {
          for (int dofi = 0; dofi < 3; dofi++) {
            for (int dofj = 0; dofj < 3; dofj++) {
              entriesLocal.emplace_back(vertexIndices[vi] * 3 + dofi, vertexIndices[vj] * 3 + dofj, 1.0);
            }
          }
        }
      }

      // upper-right block
      for (int vi = 0; vi < neleVtx; vi++) {
        for (int dofi = 0; dofi < 3; dofi++) {
          for (int parami = 0; parami < numPlasticParams; parami++) {
            entriesLocal.emplace_back(vertexIndices[vi] * 3 + dofi, n3 + ele * numPlasticParams + parami, 1.0);
            entriesLocal.emplace_back(n3 + ele * numPlasticParams + parami, vertexIndices[vi] * 3 + dofi, 1.0);
          }
        }
      }

      for (int parami = 0; parami < numPlasticParams; parami++) {
        for (int paramj = 0; paramj < numPlasticParams; paramj++) {
          entriesLocal.emplace_back(n3 + ele * numPlasticParams + parami, n3 + ele * numPlasticParams + paramj, 1.0);
        }
      }
    },
    tbb::static_partitioner());

  for (auto it = entriesTLS.begin(); it != entriesTLS.end(); ++it) {
    entries.insert(entries.end(), it->begin(), it->end());
  }

  epimpl->KTemplateHandle.hess.resize(n3 + numPlasticParams * nele, n3 + numPlasticParams * nele);
  epimpl->KTemplateHandle.hess.setFromTriplets(entries.begin(), entries.end());
  epimpl->inverseIndices.resize(nele);

  // for (int ele = 0; ele < nele; ele++) {
  tbb::parallel_for(
    0, nele, [&](int ele) {
      const int *vertexIndices = deformationModelManager->getMesh()->getVertexIndices(ele);
      DeformationModelAssemblerImpl::IndexMatrix idxM = DeformationModelAssemblerImpl::IndexMatrix::Ones() * -1;

      // upper-left block
      for (int vi = 0; vi < neleVtx; vi++) {
        for (int vj = 0; vj < neleVtx; vj++) {
          for (int dofi = 0; dofi < 3; dofi++) {
            for (int dofj = 0; dofj < 3; dofj++) {
              int globalRow = vertexIndices[vi] * 3 + dofi;
              int globalCol = vertexIndices[vj] * 3 + dofj;

              int localRow = vi * 3 + dofi;
              int localCol = vj * 3 + dofj;

              idxM(localRow, localCol) = ES::findEntryOffset(epimpl->KTemplateHandle.hess, globalRow, globalCol);
              PGO_ALOG(idxM(localRow, localCol) >= 0);
            }
          }
        }
      }

      for (int vi = 0; vi < neleVtx; vi++) {
        for (int dofi = 0; dofi < 3; dofi++) {
          for (int parami = 0; parami < numPlasticParams; parami++) {
            int globalRow = vertexIndices[vi] * 3 + dofi;
            int globalCol = ele * numPlasticParams + n3 + parami;

            int localRow = vi * 3 + dofi;
            int localCol = neleVtx * 3 + parami;

            idxM(localRow, localCol) = ES::findEntryOffset(epimpl->KTemplateHandle.hess, globalRow, globalCol);
            PGO_ALOG(idxM(localRow, localCol) >= 0);

            idxM(localCol, localRow) = ES::findEntryOffset(epimpl->KTemplateHandle.hess, globalCol, globalRow);
            PGO_ALOG(idxM(localCol, localRow) >= 0);
          }
        }
      }

      for (int parami = 0; parami < numPlasticParams; parami++) {
        for (int paramj = 0; paramj < numPlasticParams; paramj++) {
          int globalRow = n3 + ele * numPlasticParams + parami;
          int globalCol = n3 + ele * numPlasticParams + paramj;

          int localRow = neleVtx * 3 + parami;
          int localCol = neleVtx * 3 + paramj;

          idxM(localRow, localCol) = ES::findEntryOffset(epimpl->KTemplateHandle.hess, globalRow, globalCol);
          PGO_ALOG(idxM(localRow, localCol) >= 0);
        }
      }

      epimpl->inverseIndices[ele] = idxM;
    },
    tbb::static_partitioner());
}

DeformationModelAssemblerCacheData *DeformationModelAssemblerElasticFullPlastic::allocateCache() const
{
  DeformationModelAssemblerCacheData *data = new DeformationModelAssemblerCacheData;
  data->internalForceVertexLocks = std::vector<tbb::spin_mutex>(n3 / 3 + nele);
  data->stiffnessMatrixVertexRowLocks = std::vector<tbb::spin_mutex>(n3 / 3 + nele);

  data->elementCacheData.assign(nele, nullptr);
  for (int i = 0; i < nele; i++) {
    data->elementCacheData[i] = epimpl->femModels[i]->allocateCacheData();
  }

  return data;
}

void DeformationModelAssemblerElasticFullPlastic::freeCache(DeformationModelAssemblerCacheData *data) const
{
  for (int i = 0; i < nele; i++) {
    epimpl->femModels[i]->freeCacheData(data->elementCacheData[i]);
  }
  delete data;
}

void DeformationModelAssemblerElasticFullPlastic::getFixedElasticParameters(int ele, double *param) const
{
  if (numElasticParams == 0)
    return;
  for (int dof = 0; dof < numElasticParams; dof++) {
    param[dof] = impl->fixedParameters[ele * numElasticParams + dof];
  }
}

// Consider elastic dofs in the future...
void DeformationModelAssemblerElasticFullPlastic::getElasticParameters(const double *x_param, int ele, double *params) const
{
  for (int parami = 0; parami < numElasticParams; parami++) {
    params[parami] = x_param[n3 + nele * numPlasticParams + ele * numElasticParams + parami];
  }
}

void DeformationModelAssemblerElasticFullPlastic::getPlasticParameters(const double *x_param, int ele, double *params) const
{
  for (int parami = 0; parami < numPlasticParams; parami++) {
    params[parami] = x_param[n3 + ele * numPlasticParams + parami];
  }
}

double DeformationModelAssemblerElasticFullPlastic::computeEnergy(const double *x_param, DeformationModelAssemblerCacheData *data) const
{
  for (auto it = data->energyLocalBuffer.begin(); it != data->energyLocalBuffer.end(); ++it)
    *it = 0.0;

  auto localEnergyFunc = [&](int ele) {
    if (epimpl->elementFlags[ele] == 0)
      return;

    ES::V18d localp;
    for (int j = 0; j < neleVtx; j++) {
      int vid = deformationModelManager->getMesh()->getVertexIndex(ele, j);
      localp.segment<3>(j * 3) = ES::V3d(x_param[vid * 3], x_param[vid * 3 + 1], x_param[vid * 3 + 2]);
    }

    ES::V12d plasticParam = ES::V12d::Zero();
    ES::V12d elasticParam = ES::V12d::Zero();

    getFixedElasticParameters(ele, elasticParam.data());
    getPlasticParameters(x_param, ele, plasticParam.data());

    const DeformationModel *fem = epimpl->femModels[ele];
    fem->prepareData(localp.data(), plasticParam.data(), elasticParam.data(), data->elementCacheData[ele]);

    data->energyLocalBuffer.local() += fem->computeEnergy(data->elementCacheData[ele]) * epimpl->elementFlags[ele];
  };

  // for (int ele = 0; ele < tetMesh->getNumElements(); ele++) {
  tbb::parallel_for(0, nele, localEnergyFunc, data->partitioners[0]);

  double energyAll = 0;
  for (auto it = data->energyLocalBuffer.begin(); it != data->energyLocalBuffer.end(); ++it)
    energyAll += *it;

  return energyAll;
}

void DeformationModelAssemblerElasticFullPlastic::computeGradient(const double *x_param, double *grad, DeformationModelAssemblerCacheData *data) const
{
  memset(grad, 0, sizeof(double) * (n3 + nele * numPlasticParams));

  auto localGradFunc = [&](int ele) {
    if (epimpl->elementFlags[ele] == 0)
      return;

    ES::V18d localp;
    for (int j = 0; j < neleVtx; j++) {
      int vid = deformationModelManager->getMesh()->getVertexIndex(ele, j);
      localp.segment<3>(j * 3) = ES::V3d(x_param[vid * 3], x_param[vid * 3 + 1], x_param[vid * 3 + 2]);
    }

    ES::V12d plasticParam = ES::V12d::Zero();
    ES::V12d elasticParam = ES::V12d::Zero();

    getFixedElasticParameters(ele, elasticParam.data());
    getPlasticParameters(x_param, ele, plasticParam.data());

    const DeformationModel *fem = epimpl->femModels[ele];
    fem->prepareData(localp.data(), plasticParam.data(), elasticParam.data(), data->elementCacheData[ele]);

    ES::V18d localGradx = ES::V18d::Zero();
    fem->compute_dE_dx(data->elementCacheData[ele], localGradx.data());
    localGradx *= epimpl->elementFlags[ele];

    if (enableSanityCheck) {
      for (int i = 0; i < neleVtx * 3; i++) {
        if (std::isfinite(localGradx[i]) == false) {
          SPDLOG_LOGGER_ERROR(Logging::lgr(), "Encounter weird numbers.\nGrad:\n{}\nx:\n{}\n", localGradx, localp);
          SPDLOG_LOGGER_ERROR(Logging::lgr(), "Plastic param:");
          for (int j = 0; j < numPlasticParams; j++) {
            std::cout << plasticParam[j] << ',';
          }
          std::cout << std::endl;
          throw std::logic_error("Weird number");
        }
      }
    }

    const int *vertexIndices = deformationModelManager->getMesh()->getVertexIndices(ele);
    for (int v = 0; v < neleVtx; v++) {
      data->internalForceVertexLocks[vertexIndices[v]].lock();

      grad[vertexIndices[v] * 3] += localGradx[v * 3];
      grad[vertexIndices[v] * 3 + 1] += localGradx[v * 3 + 1];
      grad[vertexIndices[v] * 3 + 2] += localGradx[v * 3 + 2];

      data->internalForceVertexLocks[vertexIndices[v]].unlock();
    }

    ES::V12d localGrada;
    fem->compute_dE_da(data->elementCacheData[ele], localGrada.data());
    localGrada *= epimpl->elementFlags[ele];

    if (enableSanityCheck) {
      for (int i = 0; i < numPlasticParams; i++) {
        if (std::isfinite(localGrada[i]) == false) {
          SPDLOG_LOGGER_ERROR(Logging::lgr(), "Encounter weird numbers.\nGrad:\n{}\nx:\n{}\n", localGrada, localp);
          SPDLOG_LOGGER_ERROR(Logging::lgr(), "Plastic param:");
          for (int j = 0; j < numPlasticParams; j++) {
            std::cout << localGrada[j] << ',';
          }
          std::cout << std::endl;
          throw std::logic_error("Weird number");
        }
      }
    }

    data->internalForceVertexLocks[n3 / 3 + ele].lock();

    for (int j = 0; j < numPlasticParams; j++)
      grad[n3 + ele * numPlasticParams + j] += localGrada[j];

    data->internalForceVertexLocks[n3 / 3 + ele].unlock();
  };

  tbb::parallel_for(0, nele, localGradFunc, data->partitioners[1]);

  // LGI << Eigen::Map<ES::VXd>(grad, n3).norm();
  // LGI << Eigen::Map<ES::VXd>(grad + n3, nele * numPlasticParams).norm();

  // clean the numbers
  if (enableSanityCheck) {
    for (int i = 0; i < n3 + nele * numPlasticParams; i++) {
      int fpclass = std::fpclassify(grad[i]);
      if (fpclass == FP_INFINITE || fpclass == FP_NAN) {
        SPDLOG_LOGGER_ERROR(Logging::lgr(), "Encounter weird numbers.");
        throw std::logic_error("Encounter weird numbers.");
      }
      else if (fpclass == FP_SUBNORMAL) {
        grad[i] = 0;
      }
    }
  }
}

void DeformationModelAssemblerElasticFullPlastic::computeHessian(const double *x_param, HessianMatrixHandle *hessHandle, DeformationModelAssemblerCacheData *data) const
{
  memset(hessHandle->hess.valuePtr(), 0, sizeof(double) * hessHandle->hess.nonZeros());

  auto localHessFunc = [&](int ele) {
    if (epimpl->elementFlags[ele] == 0)
      return;

    ES::V18d localp = ES::V18d::Zero();
    for (int j = 0; j < neleVtx; j++) {
      int vid = deformationModelManager->getMesh()->getVertexIndex(ele, j);
      localp.segment<3>(j * 3) = ES::V3d(x_param[vid * 3], x_param[vid * 3 + 1], x_param[vid * 3 + 2]);
    }

    ES::V12d plasticParam = ES::V12d::Zero();
    ES::V12d elasticParam = ES::V12d::Zero();

    getFixedElasticParameters(ele, elasticParam.data());
    getPlasticParameters(x_param, ele, plasticParam.data());

    const DeformationModel *fem = epimpl->femModels[ele];
    fem->prepareData(localp.data(), plasticParam.data(), elasticParam.data(), data->elementCacheData[ele]);

    double localKData[24 * 24];
    fem->compute_d2E_dx2(data->elementCacheData[ele], localKData);

    ES::Mp<ES::MXd> localK(localKData, neleVtx * 3, neleVtx * 3);
    localK *= epimpl->elementFlags[ele];

    const int *vertexIndices = deformationModelManager->getMesh()->getVertexIndices(ele);
    const auto &idxM = epimpl->inverseIndices[ele];

    // write matrices in place
    for (int vb = 0; vb < neleVtx; vb++) {
      int vIdxB = vertexIndices[vb];
      data->stiffnessMatrixVertexRowLocks[vIdxB].lock();

      for (int va = 0; va < neleVtx; va++) {
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
    }                   // va

    double B[24 * 12];  // = { 0 };
    double C[12 * 12];  // = { 0 };

    fem->compute_d2E_dxda(data->elementCacheData[ele], B);
    fem->compute_d2E_da2(data->elementCacheData[ele], C);

    Eigen::Map<ES::MXd> BMap(B, 3 * neleVtx, numPlasticParams);
    Eigen::Map<ES::MXd> CMap(C, numPlasticParams, numPlasticParams);
    BMap *= epimpl->elementFlags[ele];
    CMap *= epimpl->elementFlags[ele];

    data->stiffnessMatrixVertexRowLocks[n3 / 3 + ele].lock();

    for (int va = 0; va < neleVtx; va++) {
      for (int i = 0; i < 3; i++) {
        for (int parami = 0; parami < numPlasticParams; parami++) {
          int local_row = va * 3 + i;
          int local_col = neleVtx * 3 + parami;

          std::ptrdiff_t offset = idxM(local_row, local_col);
          hessHandle->hess.valuePtr()[offset] += BMap(local_row, parami);

          offset = idxM(local_col, local_row);
          hessHandle->hess.valuePtr()[offset] += BMap(local_row, parami);
        }
      }
    }

    for (int parami = 0; parami < numPlasticParams; parami++) {
      for (int paramj = 0; paramj < numPlasticParams; paramj++) {
        int local_row = neleVtx * 3 + parami;
        int local_col = neleVtx * 3 + paramj;

        std::ptrdiff_t offset = idxM(local_row, local_col);
        hessHandle->hess.valuePtr()[offset] += CMap(parami, paramj);
      }
    }

    data->stiffnessMatrixVertexRowLocks[n3 / 3 + ele].unlock();
  };

  // for (int ele = 0; ele < tetMesh->getNumElements(); ele++) {
  //   localHessFunc(ele);
  // }
  tbb::parallel_for(0, nele, localHessFunc, data->partitioners[2]);

  // clean the numbers
  if (enableSanityCheck) {
    for (Eigen::Index i = 0; i < hessHandle->hess.nonZeros(); i++) {
      int fpclass = std::fpclassify(hessHandle->hess.valuePtr()[i]);
      if (fpclass == FP_INFINITE || fpclass == FP_NAN) {
        SPDLOG_LOGGER_ERROR(Logging::lgr(), "Encounter weird numbers.");
        throw std::logic_error("Encounter weird numbers.");
      }
      else if (fpclass == FP_SUBNORMAL) {
        hessHandle->hess.valuePtr()[i] = 0;
      }
    }
  }
}
