#pragma once

// Reusable IPC simulation setup extracted from the command-line tool.
#include "EigenSupport.h"
#include "runIPCExternalObjects.h"
#include "runSimVertexTrajectory.h"
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
// A `fixed-vertices` entry driven by a prescribed trajectory: one pulling
// energy per track, retargeted by the runner at every frame and switched off
// at the track's release time.
struct TrajectoryPulling
{
  std::shared_ptr<ConstraintPotentialEnergies::MultipleVertexPulling> energy;
  RunSim::VertexTrajectoryTrack track;
  RunSim::VertexTrajectory::Space space = RunSim::VertexTrajectory::Space::Position;
  EigenSupport::VXd restPositions;  // scaled rest positions of the track vertices
  double scale = 1.0;
  double coeff = 1.0;
  bool released = false;

  // Target positions in simulation coordinates for the given simulation time.
  EigenSupport::VXd targetAt(double time) const;
};

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
  std::vector<TrajectoryPulling> trajectoryPullings;
  pgo::Mesh::TriMeshGeo surfaceMesh;
  std::shared_ptr<Contact::CIPC::EmbeddedSurfaceIPCPotentialEnergy> collisionHandler;
  // Kinematic `ipc-external-objects`, placed per frame by the runner.
  std::vector<IPCKinematicExternalObject> kinematicExternalObjects;
};

IpcSimulationContext buildShellIpcSimulation(const ConfigFileJSON &jconfig);
IpcSimulationContext buildVolumeIpcSimulation(const ConfigFileJSON &jconfig);
}  // namespace RunIPCSim
}  // namespace pgo
