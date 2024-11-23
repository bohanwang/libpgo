/*
author: Bohan Wang
copyright to USC,MIT
*/

#pragma once

#include "deformationModelAssemblerImpl.h"

namespace pgo
{
namespace SolidDeformationModel
{
class DeformationModelAssemblerElasticMultiPlasticImpl : public DeformationModelAssemblerImpl
{
public:
  DeformationModelAssemblerElasticMultiPlasticImpl(const DeformationModelAssemblerImpl &parent):
    DeformationModelAssemblerImpl(parent) {}

  std::vector<ES::MXd> elementWeightMatricesElastic;
  std::vector<ES::MXd> elementWeightMatricesPlastic;
  std::vector<ES::VXd> elementWeightMatricesCompressed;
};

class DeformationModelAssemblerElasticMultiPlasticCacheData : public DeformationModelAssemblerCacheData
{
public:
  std::shared_ptr<tbb::enumerable_thread_specific<ES::MXd>> da_dzTLS;
  std::shared_ptr<tbb::enumerable_thread_specific<ES::MXd>> d2ai_dz2TLS;
  std::shared_ptr<tbb::enumerable_thread_specific<ES::MXd>> d2E_dudzTLS;
  std::shared_ptr<tbb::enumerable_thread_specific<ES::MXd>> d2E_dadzTLS;
};

}  // namespace SolidDeformationModel
}  // namespace pgo
