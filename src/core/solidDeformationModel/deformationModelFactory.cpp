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

std::unique_ptr<SimulationMesh> makeSimulationMesh(const VolumetricMeshes::VolumetricMesh &mesh)
{
  std::unique_ptr<SimulationMesh> result;
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
  std::unique_ptr<SimulationMesh> mesh,
  DeformationModelElasticMaterial elastic,
  DeformationModelPlasticMaterial plastic,
  const DeformationModelOptions &opts)
{
  if (!mesh)
    throw std::invalid_argument("makeDeformationModel: mesh is null.");

  const int nele = mesh->getNumElements();
  const int n3 = mesh->getNumVertices() * 3;

  // Capture the rest pose from the mesh before it is moved into the manager.
  ES::VXd restPosition(n3);
  for (int vi = 0; vi < mesh->getNumVertices(); vi++) {
    double p[3];
    mesh->getVertex(vi, p);
    restPosition.segment<3>(vi * 3) = ES::V3d(p[0], p[1], p[2]);
  }

  auto manager = std::make_unique<DeformationModelManager>(
    std::move(mesh), plastic, elastic,
    opts.enforceSPD ? 1 : 0);

  ES::VXd elementWeights = opts.elementWeights;
  if (elementWeights.size() == 0)
    elementWeights = ES::VXd::Ones(nele);
  else if (static_cast<int>(elementWeights.size()) != nele)
    throw std::invalid_argument("makeDeformationModel: elementWeights size does not match the element count.");

  const int numPlasticParams = manager->getNumPlasticParameters();
  ES::VXd plasticParams(static_cast<Eigen::Index>(nele) * numPlasticParams);
  if (numPlasticParams > 0) {
    plasticParams.setZero();
    const ES::M3d identity = ES::M3d::Identity();
    for (int ei = 0; ei < nele; ei++) {
      const auto *pm = dynamic_cast<const PlasticModel3DDeformationGradient *>(
        manager->getDeformationModel(ei)->getPlasticModel());
      if (!pm)
        throw std::runtime_error("makeDeformationModel: plastic model is not a PlasticModel3DDeformationGradient.");
      pm->toParam(identity.data(), plasticParams.data() + ei * numPlasticParams);
    }
  }

  auto assembler = std::make_unique<DeformationModelAssembler>(std::move(manager), elementWeights.data());

  DeformationModelBundle bundle;
  bundle.restPosition = std::move(restPosition);
  bundle.plasticParams = std::move(plasticParams);
  bundle.energy = std::make_shared<DeformationModelEnergy>(std::move(assembler), &bundle.restPosition, 0);
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
