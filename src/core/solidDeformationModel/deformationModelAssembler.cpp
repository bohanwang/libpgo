/*
author: Bohan Wang
copyright to USC,MIT
*/

#include "deformationModelAssembler.h"
#include "deformationModelAssemblerImpl.h"

#include "deformationModelManager.h"
#include "simulationMesh.h"

using namespace pgo::SolidDeformationModel;

DeformationModelAssembler::DeformationModelAssembler(const DeformationModelManager *dm, DeformationModelAssemblingType at, const double *elementFlags):
  deformationModelManager(dm), assemblingType(at)
{
  impl = new DeformationModelAssemblerImpl;
  nele = deformationModelManager->getMesh()->getNumElements();
  nvtx = deformationModelManager->getMesh()->getNumVertices();
  neleVtx = deformationModelManager->getMesh()->getNumElementVertices();
  n3 = 3 * nvtx;

  if (elementFlags) {
    impl->elementFlags.assign(elementFlags, elementFlags + nele);
  }
  else {
    impl->elementFlags.assign(nele, 1);
  }

  impl->restPositions.resize(n3);
  for (int vi = 0; vi < deformationModelManager->getMesh()->getNumVertices(); vi++) {
    ES::V3d p;
    deformationModelManager->getMesh()->getVertex(vi, p.data());
    impl->restPositions.segment<3>(vi * 3) = p;
  }

  for (int i = 0; i < nele; i++) {
    impl->femModels.push_back(deformationModelManager->getDeformationModel(i));
  }
}

DeformationModelAssembler::~DeformationModelAssembler()
{
  if (impl)
    delete impl;
}

void DeformationModelAssembler::setFixedParameters(const double *param, int size)
{
  impl->fixedParameters = Eigen::Map<const ES::VXd>(param, size);
}

int DeformationModelAssembler::getFixedParameters(double *param) const
{
  memcpy(param, impl->fixedParameters.data(), sizeof(double) * impl->fixedParameters.size());
  return (int)impl->fixedParameters.size();
}

const HessianMatrixHandle *DeformationModelAssembler::getHessianTemplate() const
{
  return &impl->KTemplateHandle;
}

void DeformationModelAssembler::compute3rdOrderTensor(const double *, const double *, HessianMatrixHandle *, DeformationModelAssemblerCacheData *) const
{
}

void DeformationModelAssembler::computeHessian(const double *, double *, DeformationModelAssemblerCacheData *) const
{
}