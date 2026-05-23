#include "setup/setup.h"

#include "barycentricCoordinates.h"
#include "configFileJSON.h"
#include "deformationModelAssembler.h"
#include "deformationModelEnergy.h"
#include "deformationModelManager.h"
#include "ipc/embeddedSurfaceFloorPotentialEnergy.h"
#include "ipc/embeddedSurfaceIPCPotentialEnergy.h"
#include "generateMassMatrix.h"
#include "libiglInterface.h"
#include "multiVertexPullingSoftConstraints.h"
#include "pgoLogging.h"
#include "setup/attachmentSetup.h"
#include "setup/floorSetup.h"
#include "setup/obstacleSetup.h"
#include "setup/setupCommon.h"
#include "setup/surfacePressureSetup.h"
#include "setup/femSetup.h"
#include "io/volumeMeshIO.h"
#include "simulationMesh.h"
#include "ipc/core/surfaceIPCCore.h"
#include "volumetricMesh.h"

#include <cmath>
#include <iostream>
#include <memory>

namespace pgo::RunIPCSim
{
namespace ES = pgo::EigenSupport;

Contact::CIPC::SurfaceIPCCore::Parameters makeVolumeIPCParams(const pgo::ConfigFileJSON &jconfig)
{
  const bool ipcHeuristic = jconfig.exist("ipc-heuristic") ? jconfig.getValue<bool>("ipc-heuristic", 1) : false;
  if (ipcHeuristic) {
    throwConfigError(
      "runIPCSim phase2tet/cubic IPC currently requires explicit `ipc-dhat` and `ipc-kappa`; `ipc-heuristic=true` is shell-only.");
  }
  if (!jconfig.exist("ipc-dhat"))
    throwConfigError("Missing required field `ipc-dhat`.");
  if (!jconfig.exist("ipc-kappa"))
    throwConfigError("Missing required field `ipc-kappa`.");

  Contact::CIPC::SurfaceIPCCore::Parameters ipcParams;
  ipcParams.dhat = jconfig.getDouble("ipc-dhat", 1);
  ipcParams.kappa = jconfig.getDouble("ipc-kappa", 1);
  ipcParams.eps_ee = 0.0;
  ipcParams.slackness = 1.0;

  if (jconfig.exist("ipc-dhat-external"))
    ipcParams.dhat_external = jconfig.getDouble("ipc-dhat-external", 1);

  return ipcParams;
}
IpcSimulationContext buildVolumeIpcSimulation(const pgo::ConfigFileJSON &jconfig)
{
  validateZeroInitialDisplacement(jconfig);

  if (!jconfig.exist("fixed-vertices"))
    throwConfigError("Missing required field `fixed-vertices`.");

  const double scale = jconfig.getDouble("scale", 1);
  const Contact::CIPC::SurfaceIPCCore::Parameters ipcParams = makeVolumeIPCParams(jconfig);
  const SolidDeformationModel::DeformationModelElasticMaterial elasticMat = parseVolumeElasticMaterial(jconfig);
  const bool enableMaterialMaxStep = parseEnableMaterialMaxStep(jconfig);
  const std::vector<ParsedFloorConfig> floorConfigs = parseFloorsConfig(jconfig);
  ParsedSurfacePressureForceConfig pressureConfig = parseSurfacePressureForceConfig(jconfig);

  const RunSim::ResolvedRunSimPaths resolvedPaths = RunSim::resolveRunSimPaths(jconfig);
  std::unique_ptr<VolumetricMeshes::VolumetricMesh> volumetricMesh =
    RunSim::loadValidatedVolumeMesh(RunSim::parseVolumeMeshInputConfig(jconfig), scale);

  pgo::Mesh::TriMeshGeo surfaceMesh;
  ES::VXd surfaceRestPositions;
  loadSurfaceMeshAndRestPositions(resolvedPaths.surfaceMeshFilename, scale, surfaceMesh, surfaceRestPositions);

  std::cout << "runIPCSim phase2 volume IPC parameters: "
            << "ipc-heuristic=false, "
            << "source=config, "
            << "enable-material-max-step=" << (enableMaterialMaxStep ? "true" : "false") << ", "
            << "ipc-dhat=" << ipcParams.dhat << ", "
            << "ipc-dhat-external=" << ipcParams.dhat_external << ", "
            << "ipc-kappa=" << ipcParams.kappa << ", "
            << "eps_ee=" << ipcParams.eps_ee << ", "
            << "slackness=" << ipcParams.slackness << ", "
            << "floors=" << floorConfigs.size();
  for (std::size_t floorIndex = 0; floorIndex < floorConfigs.size(); ++floorIndex) {
    const ParsedFloorConfig &floorConfig = floorConfigs[floorIndex];
    std::cout << ", floor[" << floorIndex << "].axis=" << floorAxisToString(floorConfig.params.floorAxis)
              << ", floor[" << floorIndex << "].side=" << floorSideToString(floorConfig.params.floorSide)
              << ", floor[" << floorIndex << "].height=" << floorConfig.params.floorHeight
              << ", floor[" << floorIndex << "].kappa=" << floorConfig.params.floorKappa;
    if (floorConfig.motionState.hasMotion) {
      std::cout << ", floor[" << floorIndex << "].motion=[" << floorConfig.motionState.heightStart
                << "->" << floorConfig.motionState.heightEnd
                << ", frames " << floorConfig.motionState.frameStart
                << "->" << floorConfig.motionState.frameEnd << "]";
    }
  }
  std::cout << ", surface-pressure-force=" << (pressureConfig.enabled ? "true" : "false");
  if (pressureConfig.enabled) {
    std::cout << ", pressure=" << pressureConfig.pressure
              << ", ramp-steps=" << pressureConfig.rampSteps
              << ", direction=opposite-surface-normal";
  }
  std::cout << std::endl;

  pgo::InterpolationCoordinates::BarycentricCoordinates bc(
    surfaceMesh.numVertices(), surfaceRestPositions.data(), volumetricMesh.get());
  ES::SpMatD W = bc.generateInterpolationMatrix();
  if (W.rows() != surfaceMesh.numVertices() * 3)
    throwConfigError("runIPCSim phase2volume setup produced an embedding matrix with unexpected row count.");
  if (W.cols() != volumetricMesh->getNumVertices() * 3)
    throwConfigError("runIPCSim phase2volume setup produced an embedding matrix with unexpected column count.");
  if (W.nonZeros() <= 0)
    throwConfigError("runIPCSim phase2volume setup produced an empty embedding matrix.");

  ES::SpMatD M;
  VolumetricMeshes::GenerateMassMatrix::computeMassMatrix(volumetricMesh.get(), M, true);

  RunSim::InitializedVolumetricSimulation initialized =
    RunSim::initializeVolumetricSimulation(*volumetricMesh, elasticMat,
      SolidDeformationModel::DeformationModelPlasticMaterial::VOLUMETRIC_DOF6,
      enableMaterialMaxStep);
  if (W.cols() != initialized.restPosition.size())
    throwConfigError("runIPCSim phase2volume setup produced an embedding matrix incompatible with simulation DOFs.");

  ES::VXd zero = ES::VXd::Zero(initialized.restPosition.size());
  ES::SpMatD K;
  initialized.elasticEnergy->createHessian(K);
  initialized.elasticEnergy->hessian(zero, K);

  std::vector<std::shared_ptr<ConstraintPotentialEnergies::MultipleVertexPulling>> pullingEnergies;
  std::vector<ES::VXd> pullingTargets;
  std::vector<ES::VXd> pullingTargetRests;
  buildPullingConstraints(jconfig, resolvedPaths.fixedVertexFilenames, initialized.restPosition, K,
    pullingEnergies, pullingTargets, pullingTargetRests);

  ES::MXd V;
  ES::MXi F;
  Mesh::triMeshGeoToMatrices(surfaceMesh, V, F);

  IpcSimulationContext context;
  context.M = std::move(M);
  context.simulationRestPosition = std::move(initialized.restPosition);
  context.surfaceRestPositions = std::move(surfaceRestPositions);
  context.plasticParams = std::move(initialized.plasticity);
  context.surfaceFromSimulationDispMap = std::move(W);
  context.simulationMeshOwner = initialized.simMesh;
  context.deformationModelManagerOwner = initialized.dmm;
  context.deformationModelAssemblerOwner = initialized.assembler;
  context.elasticEnergy = initialized.elasticEnergy;
  context.pullingEnergies = std::move(pullingEnergies);
  context.pullingTargets = std::move(pullingTargets);
  context.pullingTargetRests = std::move(pullingTargetRests);
  context.surfaceMesh = std::move(surfaceMesh);
  std::vector<bool> staticFlags;
  auto obstacles = parseExternalObjects(jconfig, scale, &staticFlags);
  const std::size_t obstacleCount = obstacles.size();
  context.collisionHandler =
    std::make_shared<Contact::CIPC::EmbeddedSurfaceIPCPotentialEnergy>(
      V, F, context.surfaceFromSimulationDispMap, ipcParams, std::move(obstacles));
  context.contactBackend = makeIpcContactBackend();
  for (std::size_t i = 0; i < staticFlags.size(); ++i)
    if (staticFlags[i])
      context.collisionHandler->markObstacleStatic(static_cast<int32_t>(i));
  for (const ParsedFloorConfig &floorConfig : floorConfigs) {
    auto floorEnergy =
      std::make_shared<Contact::CIPC::EmbeddedSurfaceFloorPotentialEnergy>(V, context.surfaceFromSimulationDispMap, floorConfig.params);
    context.floorPotentialEnergies.push_back(floorEnergy);
    context.floorMotionStates.push_back(floorConfig.motionState);
    context.extraGeneralImplicitForceModels.push_back(
      std::move(floorEnergy));
  }
  context.surfacePressureForceEnabled = pressureConfig.enabled;
  context.surfacePressureRampSteps = pressureConfig.rampSteps;
  context.surfacePressureSimulationForce = computeSurfacePressureSimulationForce(
    context.surfaceMesh, context.surfaceFromSimulationDispMap,
    pressureConfig, static_cast<int>(context.simulationRestPosition.size()));
  if (obstacleCount > 0)
    std::cout << ", obstacles=" << obstacleCount;
  return context;
}
}  // namespace pgo::RunIPCSim
