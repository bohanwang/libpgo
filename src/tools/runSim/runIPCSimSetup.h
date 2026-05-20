#pragma once

#include "EigenSupport.h"
#include "runIPCSimContactBackend.h"
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
}

namespace ConstraintPotentialEnergies
{
class MultipleVertexPulling;
}

namespace NonlinearOptimization
{
class PotentialEnergy;
}

namespace Contact
{
namespace CIPC
{
class EmbeddedSurfaceFloorPotentialEnergy;
class EmbeddedSurfaceIPCPotentialEnergy;
}
}  // namespace Contact

namespace RunIPCSim
{
struct IpcFloorMotionState
{
  bool hasMotion = false;
  double heightStart = 0.0;
  double heightEnd = 0.0;
  int frameStart = 0;
  int frameEnd = 0;
};

struct IpcSimulationContext
{
  EigenSupport::SpMatD M;
  EigenSupport::VXd simulationRestPosition;
  EigenSupport::VXd surfaceRestPositions;
  EigenSupport::VXd plasticParams;
  EigenSupport::VXd elasticParams;
  EigenSupport::SpMatD surfaceFromSimulationDispMap;
  std::shared_ptr<SolidDeformationModel::SimulationMesh> simulationMeshOwner;
  std::shared_ptr<SolidDeformationModel::DeformationModelManager> deformationModelManagerOwner;
  std::shared_ptr<SolidDeformationModel::DeformationModelAssembler> deformationModelAssemblerOwner;
  std::shared_ptr<SolidDeformationModel::DeformationModelEnergy> elasticEnergy;
  std::vector<std::shared_ptr<ConstraintPotentialEnergies::MultipleVertexPulling>> pullingEnergies;
  std::vector<EigenSupport::VXd> pullingTargets;
  std::vector<EigenSupport::VXd> pullingTargetRests;
  pgo::Mesh::TriMeshGeo surfaceMesh;
  std::shared_ptr<RunIPCSimContactBackend> contactBackend;
  std::shared_ptr<Contact::CIPC::EmbeddedSurfaceIPCPotentialEnergy> collisionHandler;
  std::vector<std::shared_ptr<NonlinearOptimization::PotentialEnergy>> extraGeneralImplicitForceModels;
  std::vector<std::shared_ptr<Contact::CIPC::EmbeddedSurfaceFloorPotentialEnergy>> floorPotentialEnergies;
  std::vector<IpcFloorMotionState> floorMotionStates;
  bool surfacePressureForceEnabled = false;
  int surfacePressureRampSteps = 1;
  EigenSupport::VXd surfacePressureSimulationForce;
};

IpcSimulationContext buildShellIpcSimulation(const ConfigFileJSON &jconfig);
IpcSimulationContext buildVolumeIpcSimulation(const ConfigFileJSON &jconfig);
IpcSimulationContext buildVolumeLegacyPenaltySimulation(const ConfigFileJSON &jconfig);
std::shared_ptr<RunIPCSimContactBackend> makeIpcContactBackend();
double floorHeightAtFrame(const IpcFloorMotionState &motion, int frame);
}  // namespace RunIPCSim
}  // namespace pgo
