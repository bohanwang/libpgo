#include "runSimFEMSetup.h"

// Shared FEM setup used by the sampled and IPC runners.
#include "cubicMesh.h"
#include "deformationModel.h"
#include "deformationModelAssembler.h"
#include "deformationModelEnergy.h"
#include "deformationModelManager.h"
#include "plasticModel3DDeformationGradient.h"
#include "simulationMesh.h"
#include "tetMesh.h"
#include "volumetricMesh.h"

#include <stdexcept>
#include <vector>

namespace pgo::RunSim
{
InitializedVolumetricSimulation initializeVolumetricSimulation(
  const VolumetricMeshes::VolumetricMesh &volumetricMesh,
  SolidDeformationModel::DeformationModelElasticMaterial elasticMat,
  SolidDeformationModel::DeformationModelPlasticMaterial plasticMat,
  bool enableMaterialMaxStep)
{
  using namespace pgo::SolidDeformationModel;
  namespace ES = pgo::EigenSupport;

  InitializedVolumetricSimulation initialized;

  switch (volumetricMesh.getElementType()) {
  case VolumetricMeshes::VolumetricMesh::TET: {
    const auto *tetMesh = dynamic_cast<const VolumetricMeshes::TetMesh *>(&volumetricMesh);
    if (!tetMesh) {
      throw std::invalid_argument("Expected a TetMesh after volumetric mesh type validation.");
    }
    initialized.simMesh.reset(loadTetMesh(tetMesh));
    break;
  }
  case VolumetricMeshes::VolumetricMesh::CUBIC: {
    const auto *cubicMesh = dynamic_cast<const VolumetricMeshes::CubicMesh *>(&volumetricMesh);
    if (!cubicMesh) {
      throw std::invalid_argument("Expected a CubicMesh after volumetric mesh type validation.");
    }
    initialized.simMesh.reset(loadCubicMesh(cubicMesh));
    break;
  }
  default:
    throw std::invalid_argument("Unsupported volumetric mesh element type in runSim FEM initialization.");
  }

  if (!initialized.simMesh) {
    throw std::runtime_error("Failed to create SimulationMesh from volumetric mesh.");
  }

  initialized.dmm = std::make_shared<DeformationModelManager>();
  initialized.dmm->setMesh(initialized.simMesh.get(), nullptr, nullptr);
  initialized.dmm->init(plasticMat, elasticMat);
  initialized.dmm->setEnforceSPD(1);

  std::vector<double> elementWeights(initialized.simMesh->getNumElements(), 1.0);
  initialized.assembler = std::make_shared<DeformationModelAssembler>(initialized.dmm, elementWeights.data());

  const int nele = initialized.simMesh->getNumElements();
  const int n3 = initialized.simMesh->getNumVertices() * 3;
  const int numPlasticParams = initialized.dmm->getNumPlasticParameters();

  initialized.plasticity.resize(nele * numPlasticParams);
  if (numPlasticParams > 0) {
    initialized.plasticity.setZero();
    const ES::M3d identity = ES::M3d::Identity();
    for (int ei = 0; ei < nele; ei++) {
      const auto *pm = dynamic_cast<const PlasticModel3DDeformationGradient *>(
        initialized.dmm->getDeformationModel(ei)->getPlasticModel());
      if (!pm) {
        throw std::runtime_error("Plastic model is not of type PlasticModel3DDeformationGradient.");
      }
      pm->toParam(identity.data(), initialized.plasticity.data() + ei * numPlasticParams);
    }
  }

  initialized.restPosition.resize(n3);
  for (int vi = 0; vi < initialized.simMesh->getNumVertices(); vi++) {
    double p[3];
    initialized.simMesh->getVertex(vi, p);
    initialized.restPosition.segment<3>(vi * 3) = ES::V3d(p[0], p[1], p[2]);
  }

  initialized.elasticEnergy = std::make_shared<DeformationModelEnergy>(initialized.assembler, &initialized.restPosition, 0);
  initialized.elasticEnergy->setEnableMaterialMaxStep(enableMaterialMaxStep);
  initialized.elasticEnergy->setPlasticParams(initialized.plasticity);

  return initialized;
}
}  // namespace pgo::RunSim
