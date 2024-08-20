/*
author: Bohan Wang
copyright to USC,MIT
*/

#pragma once

#include "EigenSupport.h"

namespace pgo
{
namespace ES = EigenSupport;

namespace SolidDeformationModel
{
class HessianMatrixHandle
{
public:
  HessianMatrixHandle(ES::SpMatD &h):
    hess(h) {}
  ES::SpMatD &hess;
};
}  // namespace SolidDeformationModel
}  // namespace pgo