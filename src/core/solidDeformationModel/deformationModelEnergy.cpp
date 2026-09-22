/*
author: Bohan Wang
copyright to USC,MIT,NUS
*/

#include "deformationModelEnergy.h"

#include "deformationModelAssembler.h"
#include "simulationMesh.h"
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

// x, grad and hess are expressed in this energy's own DOFs (size getNumDOFs()).
// The constructor offset only places those DOFs inside an aggregate such as
// PotentialEnergies, which gathers them before calling these methods.
double DeformationModelEnergy::func(EigenSupport::ConstRefVecXd x) const
{
  if (restPosition.size()) {
    ES::VXd p = restPosition + x.head(restPosition.size());
    return forceModelAssembler->computeEnergy(p.data(), plasticParams.data(), elasticParams.data());
  }
  else
    return forceModelAssembler->computeEnergy(x.data(), plasticParams.data(), elasticParams.data());
}

void DeformationModelEnergy::gradient(EigenSupport::ConstRefVecXd x, EigenSupport::RefVecXd grad) const
{
  if (restPosition.size()) {
    ES::VXd p = restPosition + x.head(restPosition.size());
    forceModelAssembler->computeGradient(p.data(), plasticParams.data(), elasticParams.data(), grad.data());
  }
  else {
    forceModelAssembler->computeGradient(x.data(), plasticParams.data(), elasticParams.data(), grad.data());
  }
}

void DeformationModelEnergy::hessian(EigenSupport::ConstRefVecXd x, EigenSupport::SpMatD &hess) const
{
  if (restPosition.size()) {
    ES::VXd p = restPosition + x.head(restPosition.size());
    forceModelAssembler->computeHessian(p.data(), plasticParams.data(), elasticParams.data(), hess);
  }
  else {
    forceModelAssembler->computeHessian(x.data(), plasticParams.data(), elasticParams.data(), hess);
  }
}

void DeformationModelEnergy::createHessian(EigenSupport::SpMatD &hess) const
{
  hess = forceModelAssembler->getHessianTemplate();
}

double DeformationModelEnergy::computeMaxStepSize(EigenSupport::ConstRefVecXd x, EigenSupport::ConstRefVecXd dx) const
{
  if (!enableMaterialMaxStep_) {
    return 1.0;
  }

  const int numDOFs = getNumDOFs();

  const ES::VXd dxLocal = assembleDirectionSlice(dx, 0, numDOFs);
  if (dxLocal.size() == 0 || dxLocal.squaredNorm() == 0.0) {
    return 1.0;
  }

  const ES::VXd absolutePositions = assembleAbsolutePositions(x, restPosition, 0, numDOFs);
  return forceModelAssembler->computeMaxStepObservation(absolutePositions.data(), dxLocal.data()).alpha;
}
