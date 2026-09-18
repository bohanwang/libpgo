#pragma once

// Reusable IPC simulation setup extracted from the command-line tool.
#include "EigenSupport.h"
#include "triMeshGeo.h"

#include <memory>
#include <vector>

namespace pgo
{
class ConfigFileJSON;

namespace SolidDeformationModel
{
class DeformationModelAssembler;
class DeformationModelEnergy;
class DeformationModelManager;
class SimulationMesh;
}  // namespace SolidDeformationModel

namespace ConstraintPotentialEnergies
{
class MultipleVertexPulling;
}

namespace Contact
{
namespace CIPC
{
class EmbeddedSurfaceIPCPotentialEnergy;
}
}  // namespace Contact

namespace RunIPCSim
{
struct IpcSimulationContext
{
  EigenSupport::SpMatD M;
  EigenSupport::VXd simulationRestPosition;
  EigenSupport::VXd initialDisplacement;
  EigenSupport::VXd surfaceRestPositions;
  EigenSupport::SpMatD surfaceFromSimulationDispMap;
  std::shared_ptr<SolidDeformationModel::SimulationMesh> simulationMeshOwner;
  std::shared_ptr<SolidDeformationModel::DeformationModelManager> deformationModelManagerOwner;
  std::shared_ptr<SolidDeformationModel::DeformationModelAssembler> deformationModelAssemblerOwner;
  std::shared_ptr<SolidDeformationModel::DeformationModelEnergy> elasticEnergy;
  std::shared_ptr<EigenSupport::SpMatD> pullingHessianBaseOwner;
  std::vector<std::shared_ptr<ConstraintPotentialEnergies::MultipleVertexPulling>> pullingEnergies;
  std::vector<EigenSupport::VXd> pullingTargets;
  std::vector<EigenSupport::VXd> pullingTargetRests;
  pgo::Mesh::TriMeshGeo surfaceMesh;
  std::shared_ptr<Contact::CIPC::EmbeddedSurfaceIPCPotentialEnergy> collisionHandler;
};

IpcSimulationContext buildShellIpcSimulation(const ConfigFileJSON &jconfig);
IpcSimulationContext buildVolumeIpcSimulation(const ConfigFileJSON &jconfig);
}  // namespace RunIPCSim
}  // namespace pgo
