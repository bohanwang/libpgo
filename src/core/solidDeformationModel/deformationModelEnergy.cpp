/*
author: Bohan Wang
copyright to USC,MIT
*/

#include "deformationModelEnergy.h"
#include "deformationModelAssembler.h"
#include "deformationModelAssemblerElastic.h"
#include "hessianMatrixHandle.h"

#include <numeric>

using namespace pgo;
using namespace pgo::NonlinearOptimization;
using namespace pgo::SolidDeformationModel;

DeformationModelEnergy::DeformationModelEnergy(std::shared_ptr<DeformationModelAssembler> fma, const ES::VXd *restp, int offset):
  forceModelAssembler(fma)
{
  allDOFs.resize(forceModelAssembler->getNumDOFs());
  std::iota(allDOFs.begin(), allDOFs.end(), offset);

  if (restp) {
    restPosition = *restp;
  }

  data = forceModelAssembler->allocateCache();
}

DeformationModelEnergy::~DeformationModelEnergy()
{
  forceModelAssembler->freeCache(data);
}

void DeformationModelEnergy::computeVonMisesStresses(EigenSupport::ConstRefVecXd x, EigenSupport::RefVecXd elementStresses) const
{
  if (restPosition.size()) {
    ES::VXd p = restPosition + x.segment(allDOFs[0], restPosition.size());
    auto ptr = std::dynamic_pointer_cast<DeformationModelAssemblerElastic>(forceModelAssembler);
    if (ptr) {
      ptr->computeVonMisesStresses(p.data(), data, elementStresses.data());
    }
  }
  else {
    auto ptr = std::dynamic_pointer_cast<DeformationModelAssemblerElastic>(forceModelAssembler);
    if (ptr) {
      ptr->computeVonMisesStresses(x.data() + allDOFs[0], data, elementStresses.data());
    }
  }
}

void DeformationModelEnergy::computeMaxStrains(EigenSupport::ConstRefVecXd x, EigenSupport::RefVecXd elementStrains) const
{
  if (restPosition.size()) {
    ES::VXd p = restPosition + x.segment(allDOFs[0], restPosition.size());
    auto ptr = std::dynamic_pointer_cast<DeformationModelAssemblerElastic>(forceModelAssembler);
    if (ptr) {
      ptr->computeMaxStrains(p.data(), data, elementStrains.data());
    }
  }
  else {
    auto ptr = std::dynamic_pointer_cast<DeformationModelAssemblerElastic>(forceModelAssembler);
    if (ptr) {
      ptr->computeMaxStrains(x.data() + allDOFs[0], data, elementStrains.data());
    }
  }
}

double DeformationModelEnergy::func(EigenSupport::ConstRefVecXd x) const
{
  if (restPosition.size()) {
    ES::VXd p = restPosition + x.segment(allDOFs[0], restPosition.size());
    return forceModelAssembler->computeEnergy(p.data(), data);
  }
  else
    return forceModelAssembler->computeEnergy(x.data() + allDOFs[0], data);
}

void DeformationModelEnergy::gradient(EigenSupport::ConstRefVecXd x, EigenSupport::RefVecXd grad) const
{
  if (restPosition.size()) {
    ES::VXd p = restPosition + x.segment(allDOFs[0], restPosition.size());
    forceModelAssembler->computeGradient(p.data(), grad.data(), data);
  }
  else {
    forceModelAssembler->computeGradient(x.data() + allDOFs[0], grad.data(), data);
  }
}

void DeformationModelEnergy::hessian(EigenSupport::ConstRefVecXd x, EigenSupport::SpMatD &hess) const
{
  HessianMatrixHandle handle(hess);

  if (restPosition.size()) {
    ES::VXd p = restPosition + x.segment(allDOFs[0], restPosition.size());
    forceModelAssembler->computeHessian(p.data(), &handle, data);
  }
  else {
    forceModelAssembler->computeHessian(x.data() + allDOFs[0], &handle, data);
  }
}

void DeformationModelEnergy::createHessian(EigenSupport::SpMatD &hess) const
{
  hess = forceModelAssembler->getHessianTemplate()->hess;
}