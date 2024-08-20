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
class DeformationModelAssemblerElasticFullPlasticImpl : public DeformationModelAssemblerImpl
{
public:
  DeformationModelAssemblerElasticFullPlasticImpl(const DeformationModelAssemblerImpl &parent):
    DeformationModelAssemblerImpl(parent) {}

  std::vector<ES::MXd> elementWeightMatricesElastic;
  std::vector<ES::VXd> elementWeightMatricesCompressed;
};
}  // namespace SolidDeformationModel
}  // namespace pgo