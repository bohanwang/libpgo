/*
author: Bohan Wang
copyright to USC
*/

#include "deformationModelFactory.h"

#include "simulationMesh.h"
#include "deformationModel.h"
#include "deformationModelAssembler.h"
#include "deformationModelEnergy.h"
#include "plasticModel3DDeformationGradient.h"
#include "cubicMesh.h"
#include "tetMesh.h"
#include "volumetricMesh.h"

#include <stdexcept>

namespace pgo::SolidDeformationModel
{
namespace ES = pgo::EigenSupport;

std::shared_ptr<SimulationMesh> makeSimulationMesh(const VolumetricMeshes::VolumetricMesh &mesh)
{
  std::shared_ptr<SimulationMesh> result;
  switch (mesh.getElementType()) {
  case VolumetricMeshes::VolumetricMesh::TET: {
    const auto *tet = dynamic_cast<const VolumetricMeshes::TetMesh *>(&mesh);
    if (!tet)
      throw std::invalid_argument("makeSimulationMesh: element type is TET but object is not a TetMesh.");
    result = loadTetMesh(tet);
    break;
  }
  case VolumetricMeshes::VolumetricMesh::CUBIC: {
    const auto *cubic = dynamic_cast<const VolumetricMeshes::CubicMesh *>(&mesh);
    if (!cubic)
      throw std::invalid_argument("makeSimulationMesh: element type is CUBIC but object is not a CubicMesh.");
    result = loadCubicMesh(cubic);
    break;
  }
  default:
    throw std::invalid_argument("makeSimulationMesh: unsupported volumetric element type.");
  }

  if (!result)
    throw std::runtime_error("makeSimulationMesh: failed to create SimulationMesh from volumetric mesh.");

  return result;
}

DeformationModelBundle makeDeformationModel(
  std::shared_ptr<SimulationMesh> mesh,
  DeformationModelElasticMaterial elastic,
  DeformationModelPlasticMaterial plastic,
  const DeformationModelOptions &opts)
{
  if (!mesh)
    throw std::invalid_argument("makeDeformationModel: mesh is null.");

  DeformationModelBundle bundle;
  bundle.mesh = std::move(mesh);

  bundle.manager = std::make_shared<DeformationModelManager>();
  bundle.manager->setMesh(bundle.mesh.get(), nullptr, nullptr);
  bundle.manager->init(plastic, elastic);
  bundle.manager->setEnforceSPD(opts.enforceSPD ? 1 : 0);

  const int nele = bundle.mesh->getNumElements();
  const int n3 = bundle.mesh->getNumVertices() * 3;

  ES::VXd elementWeights = opts.elementWeights;
  if (elementWeights.size() == 0)
    elementWeights = ES::VXd::Ones(nele);
  else if (static_cast<int>(elementWeights.size()) != nele)
    throw std::invalid_argument("makeDeformationModel: elementWeights size does not match the element count.");

  bundle.assembler = std::make_shared<DeformationModelAssembler>(bundle.manager, elementWeights.data());

  const int numPlasticParams = bundle.manager->getNumPlasticParameters();
  bundle.plasticParams.resize(static_cast<Eigen::Index>(nele) * numPlasticParams);
  if (numPlasticParams > 0) {
    bundle.plasticParams.setZero();
    const ES::M3d identity = ES::M3d::Identity();
    for (int ei = 0; ei < nele; ei++) {
      const auto *pm = dynamic_cast<const PlasticModel3DDeformationGradient *>(
        bundle.manager->getDeformationModel(ei)->getPlasticModel());
      if (!pm)
        throw std::runtime_error("makeDeformationModel: plastic model is not a PlasticModel3DDeformationGradient.");
      pm->toParam(identity.data(), bundle.plasticParams.data() + ei * numPlasticParams);
    }
  }

  bundle.restPosition.resize(n3);
  for (int vi = 0; vi < bundle.mesh->getNumVertices(); vi++) {
    double p[3];
    bundle.mesh->getVertex(vi, p);
    bundle.restPosition.segment<3>(vi * 3) = ES::V3d(p[0], p[1], p[2]);
  }

  bundle.energy = std::make_shared<DeformationModelEnergy>(bundle.assembler, &bundle.restPosition, 0);
  bundle.energy->setEnableMaterialMaxStep(opts.enableMaterialMaxStep);
  bundle.energy->setPlasticParams(bundle.plasticParams);

  return bundle;
}

DeformationModelBundle makeDeformationModel(
  const VolumetricMeshes::VolumetricMesh &mesh,
  DeformationModelElasticMaterial elastic,
  DeformationModelPlasticMaterial plastic,
  const DeformationModelOptions &opts)
{
  return makeDeformationModel(makeSimulationMesh(mesh), elastic, plastic, opts);
}
}  // namespace pgo::SolidDeformationModel
