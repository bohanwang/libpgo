/*
author: Bohan Wang
copyright to USC,MIT,NUS
*/

#include "plasticModel3DConstant.h"

#include "EigenSupport.h"

#include <cstring>

namespace ES = pgo::EigenSupport;
using namespace pgo::SolidDeformationModel;

PlasticModel3DConstant::PlasticModel3DConstant(const double Fp_[9]):
  PlasticModel3DDeformationGradient(0)
{
  std::memcpy(Fp, Fp_, sizeof(double) * 9);
  (Eigen::Map<ES::M3d>(FpInv)) = (Eigen::Map<ES::M3d>(Fp)).fullPivLu().inverse();
  detFp = (Eigen::Map<ES::M3d>(Fp)).determinant();
}

void PlasticModel3DConstant::computeA(const double *, double A[9]) const
{
  std::memcpy(A, Fp, sizeof(double) * 9);
}

void PlasticModel3DConstant::computeAInv(const double *, double AInv[9]) const
{
  std::memcpy(AInv, FpInv, sizeof(double) * 9);
}
