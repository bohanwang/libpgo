#pragma once

#include "deformationModelManager.h"
#include "EigenSupport.h"

#include <memory>

namespace pgo
{
namespace SolidDeformationModel
{
class SimulationMesh;
class DeformationModelAssembler;
class DeformationModelEnergy;
}

namespace VolumetricMeshes
{
class VolumetricMesh;
}

namespace RunSim
{
struct InitializedVolumetricSimulation
{
  // The energy is the single owning root: energy -> assembler -> manager -> mesh.
  std::shared_ptr<SolidDeformationModel::DeformationModelEnergy> elasticEnergy;
  EigenSupport::VXd plasticity;
  EigenSupport::VXd restPosition;
};

InitializedVolumetricSimulation initializeVolumetricSimulation(
  const VolumetricMeshes::VolumetricMesh &volumetricMesh,
  SolidDeformationModel::DeformationModelElasticMaterial elasticMat,
  SolidDeformationModel::DeformationModelPlasticMaterial plasticMat =
    SolidDeformationModel::DeformationModelPlasticMaterial::VOLUMETRIC_DOF6,
  bool enableMaterialMaxStep = true);
}  // namespace RunSim
}  // namespace pgo
