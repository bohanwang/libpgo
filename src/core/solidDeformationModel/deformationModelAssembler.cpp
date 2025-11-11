/*
author: Bohan Wang
copyright to USC,MIT,NUS
*/

#include "deformationModelAssembler.h"
#include "deformationModelManager.h"
#include "simulationMesh.h"
#include "deformationModel.h"
#include "elasticModel.h"
#include "plasticModel.h"

#include "pgoLogging.h"
#include "EigenSupport.h"
#include "fmtEigen.h"

#include <tbb/parallel_for.h>
#include <tbb/enumerable_thread_specific.h>

#include <atomic>

using namespace pgo::SolidDeformationModel;

namespace ES = pgo::EigenSupport;

namespace pgo::SolidDeformationModel
{
class DeformationModelAssemblerCacheData
{
public:
  tbb::enumerable_thread_specific<double> energyLocalBuffer;
  tbb::enumerable_thread_specific<ES::MXd> upperRightBlockTLS, lowerRightBlockTLS;
  tbb::enumerable_thread_specific<ES::VXd> gradBlockTLS;

  tbb::affinity_partitioner partitioners[5];
  std::vector<DeformationModel::CacheData *> elementCacheData;
};
}  // namespace pgo::SolidDeformationModel

DeformationModelAssembler::DeformationModelAssembler(std::shared_ptr<const DeformationModelManager> dm, const double *elementFlags_):
  deformationModelManager(dm)
{
  nele = deformationModelManager->getMesh()->getNumElements();
  nvtx = deformationModelManager->getMesh()->getNumVertices();
  neleVtx = deformationModelManager->getMesh()->getNumElementVertices();
  n3 = 3 * nvtx;

  numElasticParams = deformationModelManager->getDeformationModel(0)->getElasticModel()->getNumParameters();
  numPlasticParams = deformationModelManager->getDeformationModel(0)->getPlasticModel()->getNumParameters();

  if (elementFlags_) {
    elementFlags.assign(elementFlags_, elementFlags_ + nele);
  }
  else {
    elementFlags.assign(nele, 1);
  }

  restPositions.resize(n3);
  for (int vi = 0; vi < deformationModelManager->getMesh()->getNumVertices(); vi++) {
    ES::V3d p;
    deformationModelManager->getMesh()->getVertex(vi, p.data());
    restPositions.segment<3>(vi * 3) = p;
  }

  data = new DeformationModelAssemblerCacheData;

  for (int i = 0; i < nele; i++) {
    femModels.push_back(deformationModelManager->getDeformationModel(i));
    data->elementCacheData.push_back(femModels.back()->allocateCacheData());
  }

  SPDLOG_LOGGER_INFO(Logging::lgr(), "Assembler parameter:{},{}", numElasticParams, numPlasticParams);

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

  KTemplate.resize(n3, n3);
  KTemplate.setFromTriplets(entries.begin(), entries.end());

  elementKInverseIndices.resize(nele);
  for (int ele = 0; ele < nele; ele++) {
    const int *vertexIndices = deformationModelManager->getMesh()->getVertexIndices(ele);
    IndexMatrix idxM;
    idxM.setConstant(-1);

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

              idxM(localRow, localCol) = ES::findEntryOffset(KTemplate, globalRow, globalCol);
            }
            else {
              idxM(localRow, localCol) = -1;
            }
          }
        }
      }
    }

    elementKInverseIndices[ele] = idxM;
  }

  entries.clear();
  for (int ele = 0; ele < nele; ele++) {
    for (int vi = 0; vi < neleVtx; vi++) {
      int vidx = deformationModelManager->getMesh()->getVertexIndex(ele, vi);
      if (vidx < 0) {
        continue;
      }

      for (int dof = 0; dof < 3; dof++) {
        int globalRow = vidx * 3 + dof;
        for (int ep = 0; ep < numElasticParams; ep++)
          entries.emplace_back(globalRow, ele * numElasticParams + ep, 1.0);
      }
    }
  }
  dfdbTemplate.resize(n3, nele * numElasticParams);
  dfdbTemplate.setFromTriplets(entries.begin(), entries.end());

  element_dfdb_InverseIndices.resize(nele);
  for (int ele = 0; ele < nele; ele++) {
    const int *vertexIndices = deformationModelManager->getMesh()->getVertexIndices(ele);
    IndexMatrix idxM;
    idxM.setConstant(-1);

    // upper-left block
    for (int vi = 0; vi < neleVtx; vi++) {
      for (int dofi = 0; dofi < 3; dofi++) {
        for (int ep = 0; ep < numElasticParams; ep++) {
          int localRow = vi * 3 + dofi;
          int localCol = ep;

          if (vertexIndices[vi] >= 0) {
            int globalRow = vertexIndices[vi] * 3 + dofi;
            int globalCol = ele * numElasticParams + ep;

            idxM(localRow, localCol) = ES::findEntryOffset(dfdbTemplate, globalRow, globalCol);
          }
          else {
            idxM(localRow, localCol) = -1;
          }
        }
      }
    }

    element_dfdb_InverseIndices[ele] = idxM;
  }

  entries.clear();
  for (int ele = 0; ele < nele; ele++) {
    for (int vi = 0; vi < neleVtx; vi++) {
      int vidx = deformationModelManager->getMesh()->getVertexIndex(ele, vi);
      if (vidx < 0) {
        continue;
      }
      for (int dof = 0; dof < 3; dof++) {
        int globalRow = vidx * 3 + dof;
        for (int pp = 0; pp < numPlasticParams; pp++)
          entries.emplace_back(globalRow, ele * numPlasticParams + pp, 1.0);
      }
    }
  }
  dfdaTemplate.resize(n3, nele * numPlasticParams);
  dfdaTemplate.setFromTriplets(entries.begin(), entries.end());

  element_dfda_InverseIndices.resize(nele);
  for (int ele = 0; ele < nele; ele++) {
    const int *vertexIndices = deformationModelManager->getMesh()->getVertexIndices(ele);
    IndexMatrix idxM;
    idxM.setConstant(-1);

    // upper-left block
    for (int vi = 0; vi < neleVtx; vi++) {
      for (int dofi = 0; dofi < 3; dofi++) {
        for (int pp = 0; pp < numPlasticParams; pp++) {
          int localRow = vi * 3 + dofi;
          int localCol = pp;

          if (vertexIndices[vi] >= 0) {
            int globalRow = vertexIndices[vi] * 3 + dofi;
            int globalCol = ele * numPlasticParams + pp;

            idxM(localRow, localCol) = ES::findEntryOffset(dfdaTemplate, globalRow, globalCol);
          }
          else {
            idxM(localRow, localCol) = -1;
          }
        }
      }
    }

    element_dfda_InverseIndices[ele] = idxM;
  }
}

DeformationModelAssembler::~DeformationModelAssembler()
{
  if (data)
    delete data;
}

void DeformationModelAssembler::getPlasticParameters(int ele, const double *paramsAll, double *param) const
{
  for (int j = 0; j < numPlasticParams; j++) {
    param[j] = paramsAll[ele * numPlasticParams + j];
  }
}

void DeformationModelAssembler::getElasticParameters(int ele, const double *paramsAll, double *param) const
{
  for (int j = 0; j < numElasticParams; j++) {
    param[j] = paramsAll[ele * numElasticParams + j];
  }
}

double DeformationModelAssembler::computeEnergy(const double *x, const double *plasticParams, const double *elasticParams) const
{
  for (auto it = data->energyLocalBuffer.begin(); it != data->energyLocalBuffer.end(); ++it)
    *it = 0.0;

  auto localEnergyFunc = [this, x, plasticParams, elasticParams](int ele) {
    if (elementFlags[ele] == 0)
      return;

    ES::V24d localp;
    for (int j = 0; j < neleVtx; j++) {
      int vid = deformationModelManager->getMesh()->getVertexIndex(ele, j);
      if (vid >= 0)
        localp.segment<3>(j * 3) = ES::V3d(x[vid * 3], x[vid * 3 + 1], x[vid * 3 + 2]);
      else
        localp.segment<3>(j * 3).setZero();
    }

    double plasticParam[20];
    double elasticParam[20];
    getPlasticParameters(ele, plasticParams, plasticParam);
    getElasticParameters(ele, elasticParams, elasticParam);

    const DeformationModel *fem = femModels[ele];
    fem->prepareData(localp.data(), plasticParam, elasticParam, data->elementCacheData[ele]);
    double energy = fem->computeEnergy(data->elementCacheData[ele]);

    data->energyLocalBuffer.local() += energy * elementFlags[ele];
  };

  // for (int ele = 0; ele < nele; ele++)
  //   localEnergyFunc(ele);
  tbb::parallel_for(0, nele, localEnergyFunc, data->partitioners[0]);

  double energyAll = 0;
  for (auto it = data->energyLocalBuffer.begin(); it != data->energyLocalBuffer.end(); ++it)
    energyAll += *it;

  return energyAll;
}

void DeformationModelAssembler::computeGradient(const double *x, const double *plasticParams, const double *elasticParams, double *grad) const
{
  memset(grad, 0, sizeof(double) * n3);
  auto localGradFunc = [this, x, plasticParams, elasticParams, grad](int ele) {
    if (elementFlags[ele] == 0)
      return;

    ES::V24d localp;
    for (int j = 0; j < neleVtx; j++) {
      int vid = deformationModelManager->getMesh()->getVertexIndex(ele, j);

      if (vid >= 0)
        localp.segment<3>(j * 3) = ES::V3d(x[vid * 3], x[vid * 3 + 1], x[vid * 3 + 2]);
      else
        localp.segment<3>(j * 3).setZero();
    }

    double plasticParam[20];
    double elasticParam[20];
    getPlasticParameters(ele, plasticParams, plasticParam);
    getElasticParameters(ele, elasticParams, elasticParam);

    // std::cout << ele << " Fp:";
    // for (int i = 0; i < numPlasticParams; i++) {
    //   std::cout << plasticParam[i] << ',';
    // }
    // std::cout << std::endl;

    const DeformationModel *fem = femModels[ele];
    fem->prepareData(localp.data(), plasticParam, elasticParam, data->elementCacheData[ele]);

    ES::V18d localGradx;
    fem->compute_dE_dx(data->elementCacheData[ele], localGradx.data());
    localGradx *= elementFlags[ele];

    if (enableSanityCheck) {
      for (int i = 0; i < neleVtx * 3; i++) {
        if (std::isfinite(localGradx[i]) == false) {
          SPDLOG_LOGGER_ERROR(Logging::lgr(), "Ele: {}", ele);
          SPDLOG_LOGGER_ERROR(Logging::lgr(), "Encounter weird numbers.\nGrad:\n{}\n;x:{}\n", localGradx, localp);
          SPDLOG_LOGGER_ERROR(Logging::lgr(), "Plastic param: {}, {}, {}\n", plasticParam[0], plasticParam[1], plasticParam[2]);
          // exit(EXIT_FAILURE);
        }
      }
    }

    const int *vertexIndices = deformationModelManager->getMesh()->getVertexIndices(ele);
    for (int v = 0; v < neleVtx; v++) {
      if (vertexIndices[v] >= 0) {
        // data->internalForceVertexLocks[vertexIndices[v]].lock();

        for (int dof = 0; dof < 3; dof++) {
          std::atomic_ref<double> atomicGrad(grad[vertexIndices[v] * 3 + dof]);
          atomicGrad.fetch_add(localGradx[v * 3 + dof]);
        }

        // data->internalForceVertexLocks[vertexIndices[v]].unlock();
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
        throw std::logic_error("Encounter weird numbers.");
      }
      else if (fpclass == FP_SUBNORMAL) {
        grad[i] = 0;
      }
    }
  }
}

void DeformationModelAssembler::computeHessian(const double *x, const double *plasticParams, const double *elasticParams, EigenSupport::SpMatD &hess) const
{
  memset(hess.valuePtr(), 0, sizeof(double) * hess.nonZeros());

  auto localHessFunc = [this, x, plasticParams, elasticParams, &hess](int ele) {
    if (elementFlags[ele] == 0)
      return;

    ES::V24d localp;
    for (int j = 0; j < neleVtx; j++) {
      int vid = deformationModelManager->getMesh()->getVertexIndex(ele, j);

      if (vid >= 0)
        localp.segment<3>(j * 3) = ES::V3d(x[vid * 3], x[vid * 3 + 1], x[vid * 3 + 2]);
      else
        localp.segment<3>(j * 3).setZero();
    }

    double plasticParam[20];
    double elasticParam[20];
    getPlasticParameters(ele, plasticParams, plasticParam);
    getElasticParameters(ele, elasticParams, elasticParam);

    const DeformationModel *fem = femModels[ele];
    fem->prepareData(localp.data(), plasticParam, elasticParam, data->elementCacheData[ele]);

    double localKData[24 * 24];
    fem->compute_d2E_dx2(data->elementCacheData[ele], localKData);

    ES::Mp<ES::MXd> localK(localKData, neleVtx * 3, neleVtx * 3);
    localK *= elementFlags[ele];

    const int *vertexIndices = deformationModelManager->getMesh()->getVertexIndices(ele);
    const auto &idxM = elementKInverseIndices[ele];

    // write matrices in place
    for (int vb = 0; vb < neleVtx; vb++) {
      int vIdxB = vertexIndices[vb];
      if (vIdxB < 0)
        continue;

      // data->stiffnessMatrixVertexRowLocks[vIdxB].lock();

      for (int va = 0; va < neleVtx; va++) {
        int vIdxA = vertexIndices[va];
        if (vIdxA < 0)
          continue;

        for (int j = 0; j < 3; j++) {
          for (int i = 0; i < 3; i++) {
            int local_col = 3 * vb + j;
            int local_row = 3 * va + i;

            std::ptrdiff_t offset = idxM(local_row, local_col);
            if (offset >= 0) {
              std::atomic_ref<double> hessRef(hess.valuePtr()[offset]);
              hessRef.fetch_add(localK(local_row, local_col));
            }

          }  // i
        }  // j
      }  // vb

      // data->stiffnessMatrixVertexRowLocks[vIdxB].unlock();
    }  // va
  };

  // for (int ele = 0; ele < nele; ele++) {
  //   localHessFunc(ele);
  // }

  tbb::parallel_for(0, nele, localHessFunc, data->partitioners[2]);

  // clean the numbers
  if (enableSanityCheck) {
    for (Eigen::Index i = 0; i < hess.nonZeros(); i++) {
      int fpclass = std::fpclassify(hess.valuePtr()[i]);
      if (fpclass == FP_INFINITE || fpclass == FP_NAN) {
        SPDLOG_LOGGER_ERROR(Logging::lgr(), "Encounter weird numbers at {}: {}", i, hess.valuePtr()[i]);
        throw std::logic_error("Encounter weird numbers.");
      }
      else if (fpclass == FP_SUBNORMAL) {
        hess.valuePtr()[i] = 0;
      }
    }
  }
}

void DeformationModelAssembler::compute_df_da(const double *x, const double *plasticParams, const double *elasticParams, EigenSupport::SpMatD &hess) const
{
  memset(hess.valuePtr(), 0, sizeof(double) * hess.nonZeros());

  auto localHessFunc = [this, x, plasticParams, elasticParams, &hess](int ele) {
    if (elementFlags[ele] == 0)
      return;

    ES::V24d localp;
    for (int j = 0; j < neleVtx; j++) {
      int vid = deformationModelManager->getMesh()->getVertexIndex(ele, j);

      if (vid >= 0)
        localp.segment<3>(j * 3) = ES::V3d(x[vid * 3], x[vid * 3 + 1], x[vid * 3 + 2]);
      else
        localp.segment<3>(j * 3).setZero();
    }

    double plasticParam[20];
    double elasticParam[20];
    getPlasticParameters(ele, plasticParams, plasticParam);
    getElasticParameters(ele, elasticParams, elasticParam);

    const DeformationModel *fem = femModels[ele];
    fem->prepareData(localp.data(), plasticParam, elasticParam, data->elementCacheData[ele]);

    double localKData[24 * 24];
    fem->compute_d2E_dxda(data->elementCacheData[ele], localKData);

    ES::Mp<ES::MXd> localK(localKData, neleVtx * 3, numPlasticParams);
    localK *= elementFlags[ele];

    const int *vertexIndices = deformationModelManager->getMesh()->getVertexIndices(ele);
    const auto &idxM = element_dfda_InverseIndices[ele];

    // write matrices in place
    for (int vb = 0; vb < neleVtx; vb++) {
      int vIdxB = vertexIndices[vb];
      if (vIdxB < 0)
        continue;

      for (int va = 0; va < numPlasticParams; va++) {
        for (int j = 0; j < 3; j++) {
          int local_row = 3 * vb + j;
          int local_col = va;

          std::ptrdiff_t offset = idxM(local_row, local_col);
          if (offset >= 0) {
            std::atomic_ref<double> hessRef(hess.valuePtr()[offset]);
            hessRef.fetch_add(localK(local_row, local_col));
          }
        }  // j
      }  // va
    }  // vb
  };

  for (int ele = 0; ele < nele; ele++) {
    localHessFunc(ele);
  }

  // tbb::parallel_for(0, nele, localHessFunc, data->partitioners[2]);

  // clean the numbers
  if (enableSanityCheck) {
    for (Eigen::Index i = 0; i < hess.nonZeros(); i++) {
      int fpclass = std::fpclassify(hess.valuePtr()[i]);
      if (fpclass == FP_INFINITE || fpclass == FP_NAN) {
        SPDLOG_LOGGER_ERROR(Logging::lgr(), "Encounter weird numbers at {}: {}", i, hess.valuePtr()[i]);
        throw std::logic_error("Encounter weird numbers.");
      }
      else if (fpclass == FP_SUBNORMAL) {
        hess.valuePtr()[i] = 0;
      }
    }
  }
}

void DeformationModelAssembler::compute_df_db(const double *x, const double *plasticParams, const double *elasticParams, EigenSupport::SpMatD &hess) const
{
  memset(hess.valuePtr(), 0, sizeof(double) * hess.nonZeros());

  auto localHessFunc = [this, x, plasticParams, elasticParams, &hess](int ele) {
    if (elementFlags[ele] == 0)
      return;

    ES::V24d localp;
    for (int j = 0; j < neleVtx; j++) {
      int vid = deformationModelManager->getMesh()->getVertexIndex(ele, j);

      if (vid >= 0)
        localp.segment<3>(j * 3) = ES::V3d(x[vid * 3], x[vid * 3 + 1], x[vid * 3 + 2]);
      else
        localp.segment<3>(j * 3).setZero();
    }

    double plasticParam[20];
    double elasticParam[20];
    getPlasticParameters(ele, plasticParams, plasticParam);
    getElasticParameters(ele, elasticParams, elasticParam);

    const DeformationModel *fem = femModels[ele];
    fem->prepareData(localp.data(), plasticParam, elasticParam, data->elementCacheData[ele]);

    double localKData[24 * 24];
    fem->compute_d2E_dxdb(data->elementCacheData[ele], localKData);

    ES::Mp<ES::MXd> localK(localKData, neleVtx * 3, numElasticParams);
    localK *= elementFlags[ele];

    const int *vertexIndices = deformationModelManager->getMesh()->getVertexIndices(ele);
    const auto &idxM = element_dfdb_InverseIndices[ele];

    // write matrices in place
    for (int vb = 0; vb < neleVtx; vb++) {
      int vIdxB = vertexIndices[vb];
      if (vIdxB < 0)
        continue;

      for (int va = 0; va < numElasticParams; va++) {
        for (int j = 0; j < 3; j++) {
          int local_row = 3 * vb + j;
          int local_col = va;

          std::ptrdiff_t offset = idxM(local_row, local_col);
          if (offset >= 0) {
            std::atomic_ref<double> hessRef(hess.valuePtr()[offset]);
            hessRef.fetch_add(localK(local_row, local_col));
          }
        }  // j
      }  // va
    }  // vb
  };

  for (int ele = 0; ele < nele; ele++) {
    localHessFunc(ele);
  }

  // tbb::parallel_for(0, nele, localHessFunc, data->partitioners[2]);

  // clean the numbers
  if (enableSanityCheck) {
    for (Eigen::Index i = 0; i < hess.nonZeros(); i++) {
      int fpclass = std::fpclassify(hess.valuePtr()[i]);
      if (fpclass == FP_INFINITE || fpclass == FP_NAN) {
        SPDLOG_LOGGER_ERROR(Logging::lgr(), "Encounter weird numbers at {}: {}", i, hess.valuePtr()[i]);
        throw std::logic_error("Encounter weird numbers.");
      }
      else if (fpclass == FP_SUBNORMAL) {
        hess.valuePtr()[i] = 0;
      }
    }
  }
}

void DeformationModelAssembler::computeVonMisesStresses(const double *x, const double *plasticParams, const double *elasticParams, double *elementStresses) const
{
}

void DeformationModelAssembler::computeMaxStrains(const double *x, const double *plasticParams, const double *elasticParams, double *elementStrain) const
{
}