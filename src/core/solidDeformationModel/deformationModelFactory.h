/*
author: Bohan Wang
copyright to USC
*/

#pragma once

#include "deformationModelManager.h"  // DeformationModelElasticMaterial / DeformationModelPlasticMaterial
#include "EigenSupport.h"

#include <memory>

namespace pgo
{
namespace VolumetricMeshes
{
class VolumetricMesh;
}

namespace SolidDeformationModel
{
class SimulationMesh;
class DeformationModelAssembler;
class DeformationModelEnergy;

struct DeformationModelOptions
{
  bool enforceSPD = true;
  bool enableMaterialMaxStep = true;
  // Per-element assembler weights; empty means all ones.
  EigenSupport::VXd elementWeights;
};

// One ready-to-use FEM deformation model. The energy is the single owning root of
// the unique_ptr spine (energy -> assembler -> manager -> mesh); borrow the inner
// objects via energy->assembler() / .getDeformationModelManager() / .getMesh(). The
// energy already has its rest pose and (identity) plastic parameters applied.
struct DeformationModelBundle
{
  std::shared_ptr<DeformationModelEnergy> energy;
  EigenSupport::VXd restPosition;
  EigenSupport::VXd plasticParams;
};

// Build a SimulationMesh from a volumetric mesh, dispatching on element type.
std::unique_ptr<SimulationMesh> makeSimulationMesh(const VolumetricMeshes::VolumetricMesh &mesh);

// Collapse the SimulationMesh -> manager -> assembler -> energy construction chain
// (with its init-order and ownership requirements) into a single call.
DeformationModelBundle makeDeformationModel(
  const VolumetricMeshes::VolumetricMesh &mesh,
  DeformationModelElasticMaterial elastic,
  DeformationModelPlasticMaterial plastic = DeformationModelPlasticMaterial::VOLUMETRIC_DOF6,
  const DeformationModelOptions &opts = {});

// Same, taking ownership of an already-built SimulationMesh.
DeformationModelBundle makeDeformationModel(
  std::unique_ptr<SimulationMesh> mesh,
  DeformationModelElasticMaterial elastic,
  DeformationModelPlasticMaterial plastic = DeformationModelPlasticMaterial::VOLUMETRIC_DOF6,
  const DeformationModelOptions &opts = {});
}  // namespace SolidDeformationModel
}  // namespace pgo
