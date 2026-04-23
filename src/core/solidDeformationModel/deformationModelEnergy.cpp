/*
author: Bohan Wang
copyright to USC,MIT,NUS
*/

#include "deformationModelEnergy.h"

#include "deformationModelAssembler.h"

#include <numeric>

using namespace pgo;
using namespace pgo::NonlinearOptimization;
using namespace pgo::SolidDeformationModel;
namespace ES = pgo::EigenSupport;

namespace
{
ES::VXd assembleAbsolutePositions(ES::ConstRefVecXd x, const ES::VXd &restPosition, int offset, int numDOFs)
{
  if (restPosition.size())
    return restPosition + x.segment(offset, restPosition.size());

  return ES::VXd(Eigen::Map<const ES::VXd>(x.data() + offset, numDOFs));
}

ES::VXd assembleDirectionSlice(ES::ConstRefVecXd dx, int offset, int numDOFs)
{
  return ES::VXd(Eigen::Map<const ES::VXd>(dx.data() + offset, numDOFs));
}
}  // namespace

DeformationModelEnergy::DeformationModelEnergy(std::shared_ptr<DeformationModelAssembler> fma, const ES::VXd *restp, int offset):
  forceModelAssembler(fma)
{
  allDOFs.resize(forceModelAssembler->getNumDOFs());
  std::iota(allDOFs.begin(), allDOFs.end(), offset);

  if (restp) {
    restPosition = *restp;
  }
}

DeformationModelEnergy::~DeformationModelEnergy()
{
}

double DeformationModelEnergy::func(EigenSupport::ConstRefVecXd x) const
{
  if (restPosition.size()) {
    ES::VXd p = restPosition + x.segment(allDOFs[0], restPosition.size());
    return forceModelAssembler->computeEnergy(p.data(), plasticParams.data(), elasticParams.data());
  }
  else
    return forceModelAssembler->computeEnergy(x.data() + allDOFs[0], plasticParams.data(), elasticParams.data());
}

void DeformationModelEnergy::gradient(EigenSupport::ConstRefVecXd x, EigenSupport::RefVecXd grad) const
{
  if (restPosition.size()) {
    ES::VXd p = restPosition + x.segment(allDOFs[0], restPosition.size());
    forceModelAssembler->computeGradient(p.data(), plasticParams.data(), elasticParams.data(), grad.data());
  }
  else {
    forceModelAssembler->computeGradient(x.data() + allDOFs[0], plasticParams.data(), elasticParams.data(), grad.data());
  }
}

void DeformationModelEnergy::hessian(EigenSupport::ConstRefVecXd x, EigenSupport::SpMatD &hess) const
{
  if (restPosition.size()) {
    ES::VXd p = restPosition + x.segment(allDOFs[0], restPosition.size());
    forceModelAssembler->computeHessian(p.data(), plasticParams.data(), elasticParams.data(), hess);
  }
  else {
    forceModelAssembler->computeHessian(x.data() + allDOFs[0], plasticParams.data(), elasticParams.data(), hess);
  }
}

void DeformationModelEnergy::createHessian(EigenSupport::SpMatD &hess) const
{
  hess = forceModelAssembler->getHessianTemplate();
}

double DeformationModelEnergy::computeMaxStepSize(EigenSupport::ConstRefVecXd x, EigenSupport::ConstRefVecXd dx) const
{
  if (!enableMaterialMaxStep_)
    return 1.0;

  const int offset = allDOFs.empty() ? 0 : allDOFs[0];
  const int numDOFs = getNumDOFs();

  const ES::VXd dxLocal = assembleDirectionSlice(dx, offset, numDOFs);
  if (dxLocal.size() == 0 || dxLocal.squaredNorm() == 0.0)
    return 1.0;

  const ES::VXd absolutePositions = assembleAbsolutePositions(x, restPosition, offset, numDOFs);
  const double maxStepSize = forceModelAssembler->computeMaxStepSize(absolutePositions.data(), dxLocal.data());

  if (maxStepSize < 1.0)
    materialClampCount_.fetch_add(1, std::memory_order_relaxed);

  return maxStepSize;
}
