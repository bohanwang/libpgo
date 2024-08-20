/*
author: Bohan Wang
copyright to USC,MIT
*/

#pragma once

#include "deformationModelAssemblerImpl.h"

#include <vector>

namespace pgo
{
namespace ES = pgo::EigenSupport;

namespace SolidDeformationModel
{
class DeformationModelAssemblerElasticImpl : public DeformationModelAssemblerImpl
{
public:
  DeformationModelAssemblerElasticImpl(const DeformationModelAssemblerImpl &parent):
    DeformationModelAssemblerImpl(parent) {}

  std::vector<ES::MXd> elementWeightMatricesPlastic, elementWeightMatricesElastic;
  std::vector<ES::VXd> elementWeightMatricesCompressed;
};

}  // namespace SolidDeformationModel
}  // namespace pgo
