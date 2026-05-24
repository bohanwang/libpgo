#include "setup/femSetup.h"

#include "deformationModelFactory.h"

#include <utility>

namespace pgo::RunSim
{
InitializedVolumetricSimulation initializeVolumetricSimulation(
  const VolumetricMeshes::VolumetricMesh &volumetricMesh,
  SolidDeformationModel::DeformationModelElasticMaterial elasticMat,
  SolidDeformationModel::DeformationModelPlasticMaterial plasticMat,
  bool enableMaterialMaxStep)
{
  using namespace pgo::SolidDeformationModel;

  DeformationModelOptions opts;
  opts.enforceSPD = true;
  opts.enableMaterialMaxStep = enableMaterialMaxStep;

  DeformationModelBundle bundle = makeDeformationModel(volumetricMesh, elasticMat, plasticMat, opts);

  InitializedVolumetricSimulation initialized;
  initialized.simMesh = std::move(bundle.mesh);
  initialized.dmm = std::move(bundle.manager);
  initialized.assembler = std::move(bundle.assembler);
  initialized.elasticEnergy = std::move(bundle.energy);
  initialized.plasticity = std::move(bundle.plasticParams);
  initialized.restPosition = std::move(bundle.restPosition);
  return initialized;
}
}  // namespace pgo::RunSim
