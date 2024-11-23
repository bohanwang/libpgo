/*
author: Bohan Wang
copyright to USC,MIT
*/

#pragma once

#include "deformationModel.h"
#include "hessianMatrixHandle.h"

#include "EigenSupport.h"

#include <tbb/spin_mutex.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/partitioner.h>

#include <vector>

namespace pgo
{
namespace ES = EigenSupport;

namespace SolidDeformationModel
{
class DeformationModelAssemblerImpl
{
public:
  DeformationModelAssemblerImpl():
    KTemplateHandle(KTemplate) {}

  DeformationModelAssemblerImpl(const DeformationModelAssemblerImpl &other):
    KTemplate(other.KTemplate), fixedParameters(other.fixedParameters), restPositions(other.restPositions),
    KTemplateHandle(KTemplate), inverseIndices(other.inverseIndices),
    upperRightBlockInverseIndices(other.upperRightBlockInverseIndices),
    lowerLeftBlockInverseIndices(other.lowerLeftBlockInverseIndices),
    lowerRightBlockInverseIndices(other.lowerRightBlockInverseIndices),
    femModels(other.femModels), elementFlags(other.elementFlags) {}

  virtual ~DeformationModelAssemblerImpl() {}

  ES::SpMatD KTemplate;
  ES::VXd fixedParameters;
  ES::VXd restPositions;

  HessianMatrixHandle KTemplateHandle;

  typedef Eigen::Matrix<std::ptrdiff_t, 24, 24> IndexMatrix;
  typedef Eigen::Matrix<std::ptrdiff_t, Eigen::Dynamic, Eigen::Dynamic> DynamicIndexMatrix;

  std::vector<IndexMatrix, Eigen::aligned_allocator<IndexMatrix>> inverseIndices;
  DynamicIndexMatrix upperRightBlockInverseIndices, lowerLeftBlockInverseIndices, lowerRightBlockInverseIndices;

  std::vector<const DeformationModel *> femModels;
  std::vector<double> elementFlags;
};

class DeformationModelAssemblerCacheData
{
public:
  virtual ~DeformationModelAssemblerCacheData() {}

  tbb::enumerable_thread_specific<double> energyLocalBuffer;
  tbb::enumerable_thread_specific<EigenSupport::MXd> upperRightBlockTLS, lowerRightBlockTLS;
  tbb::enumerable_thread_specific<EigenSupport::VXd> gradBlockTLS;

  tbb::affinity_partitioner partitioners[5];
  std::vector<tbb::spin_mutex> internalForceVertexLocks, stiffnessMatrixVertexRowLocks;
  std::vector<DeformationModelCacheData *> elementCacheData;
};

}  // namespace SolidDeformationModel
}  // namespace pgo
