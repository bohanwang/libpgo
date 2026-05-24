#include "setup/setup.h"

#include "configFileJSON.h"
#include "deformationModelAssembler.h"
#include "deformationModelEnergy.h"
#include "deformationModelManager.h"
#include "embeddedSurfaceFloorPotentialEnergy.h"
#include "ipc/embeddedSurfaceIPCPotentialEnergy.h"
#include "libiglInterface.h"
#include "multiVertexPullingSoftConstraints.h"
#include "pgoLogging.h"
#include "setup/attachmentSetup.h"
#include "setup/floorSetup.h"
#include "setup/obstacleSetup.h"
#include "setup/setupCommon.h"
#include "setup/surfacePressureSetup.h"
#include "simulationMesh.h"
#include "ipc/core/surfaceIPCCore.h"

#include <cmath>
#include <iostream>
#include <memory>
#include <string>

namespace pgo::RunIPCSim
{
namespace ES = pgo::EigenSupport;

Contact::IPC::SurfaceIPCCore::Parameters makeShellIPCParams(
  const pgo::ConfigFileJSON &jconfig, const pgo::Mesh::BoundingBox &surfaceBox)
{
  constexpr double kShellYoungsModulus = 1000000.0;
  constexpr double kShellThickness = 3e-3;

  Contact::IPC::SurfaceIPCCore::Parameters ipcParams;
  const bool ipcHeuristic = jconfig.exist("ipc-heuristic") ? jconfig.getValue<bool>("ipc-heuristic", 1) : false;
  if (ipcHeuristic) {
    ipcParams.dhat = surfaceBox.sides().norm() * 1e-3;
    ipcParams.kappa = kShellYoungsModulus * kShellThickness;
  }
  else {
    if (!jconfig.exist("ipc-dhat"))
      throwConfigError("Missing required field `ipc-dhat`.");
    if (!jconfig.exist("ipc-kappa"))
      throwConfigError("Missing required field `ipc-kappa`.");
    ipcParams.dhat = jconfig.getDouble("ipc-dhat", 1);
    ipcParams.kappa = jconfig.getDouble("ipc-kappa", 1);
  }
  ipcParams.eps_ee = 0.0;
  ipcParams.slackness = 1.0;

  if (jconfig.exist("ipc-dhat-external"))
    ipcParams.dhat_external = jconfig.getDouble("ipc-dhat-external", 1);

  return ipcParams;
}
IpcSimulationContext buildShellIpcSimulation(const pgo::ConfigFileJSON &jconfig)
{
  validateZeroInitialDisplacement(jconfig);
  const ParsedSurfacePressureForceConfig pressureConfig = parseSurfacePressureForceConfig(jconfig);
  if (pressureConfig.enabled)
    throwConfigError("`surface-pressure-force` is only supported for runIPCSim volume simulations.");

  if (jconfig.exist("tet-mesh") || jconfig.exist("cubic-mesh")) {
    throwConfigError("runIPCSim phase2shell setup cannot consume tet/cubic mesh inputs.");
  }
  if (!jconfig.exist("fixed-vertices"))
    throwConfigError("Missing required field `fixed-vertices`.");

  const double scale = jconfig.getDouble("scale", 1);
  PGO_ALOG(std::abs(scale - 1.0) < 1e-6);

  const std::string material = jconfig.getString("elastic-material");
  if (material != "koiter-stvk")
    throwConfigError("runIPCSim phase2shell path only supports `elastic-material = koiter-stvk`.");

  const std::string surfaceMeshFilename = jconfig.getResolvedPath("surface-mesh", 1);

  pgo::Mesh::TriMeshGeo surfaceMesh;
  ES::VXd surfaceRestPositions;
  loadSurfaceMeshAndRestPositions(surfaceMeshFilename, 1.0, surfaceMesh, surfaceRestPositions);
  const pgo::Mesh::BoundingBox surfaceBox(surfaceMesh.positions());
  const bool ipcHeuristic = jconfig.exist("ipc-heuristic") ? jconfig.getValue<bool>("ipc-heuristic", 1) : false;
  const bool enableMaterialMaxStep = parseEnableMaterialMaxStep(jconfig);
  const Contact::IPC::SurfaceIPCCore::Parameters ipcParams = makeShellIPCParams(jconfig, surfaceBox);
  const std::vector<ParsedFloorConfig> floorConfigs = parseFloorsConfig(jconfig);

  std::cout << "runIPCSim phase2 shell IPC parameters: "
            << "ipc-heuristic=" << (ipcHeuristic ? "true" : "false") << ", "
            << "source=" << (ipcHeuristic ? "heuristic" : "config") << ", "
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
  std::cout << std::endl;

  SolidDeformationModel::SimulationMeshENuhMaterial matParam(10000, 0.4, 0.001);
  std::shared_ptr<SolidDeformationModel::SimulationMesh> simMesh(
    SolidDeformationModel::loadShellMesh(surfaceMesh, &matParam));
  std::shared_ptr<SolidDeformationModel::DeformationModelManager> dmm =
    std::make_shared<SolidDeformationModel::DeformationModelManager>();

  dmm->setMesh(simMesh.get(), nullptr, nullptr);
  dmm->init(pgo::SolidDeformationModel::DeformationModelPlasticMaterial::SHELL_FF_DOF0,
    pgo::SolidDeformationModel::DeformationModelElasticMaterial::KOITER_STVK);
  dmm->setEnforceSPD(1);

  std::vector<double> elementWeights(simMesh->getNumElements(), 1.0);
  std::shared_ptr<SolidDeformationModel::DeformationModelAssembler> assembler =
    std::make_shared<SolidDeformationModel::DeformationModelAssembler>(dmm, elementWeights.data());

  const int n = simMesh->getNumVertices();
  const int n3 = n * 3;
  const int nele = simMesh->getNumElements();
  constexpr double kShellYoungsModulus = 1000000.0;
  constexpr double kShellThickness = 3e-3;

  ES::VXd elasticParams(5 * nele);
  for (int ei = 0; ei < nele; ++ei) {
    const double E_bend = kShellYoungsModulus;
    const double nu = 0.4;

    elasticParams[ei * 5 + 0] = kShellYoungsModulus;
    elasticParams[ei * 5 + 1] = nu;
    elasticParams[ei * 5 + 2] = E_bend;
    elasticParams[ei * 5 + 3] = nu;
    elasticParams[ei * 5 + 4] = kShellThickness;
  }

  ES::VXd simulationRestPosition(n3);
  for (int vi = 0; vi < n; ++vi) {
    double p[3];
    simMesh->getVertex(vi, p);
    simulationRestPosition.segment<3>(vi * 3) = ES::V3d(p[0], p[1], p[2]);
  }

  std::shared_ptr<SolidDeformationModel::DeformationModelEnergy> elasticEnergy =
    std::make_shared<SolidDeformationModel::DeformationModelEnergy>(assembler, &simulationRestPosition, 0);
  elasticEnergy->setEnableMaterialMaxStep(enableMaterialMaxStep);
  elasticEnergy->setElasticParams(elasticParams);

  ES::VXd zero = ES::VXd::Zero(n3);
  ES::SpMatD K;
  elasticEnergy->createHessian(K);
  elasticEnergy->hessian(zero, K);

  std::vector<std::string> fixedVertexFilenames;
  for (const auto &fv : jconfig.handle()["fixed-vertices"])
    fixedVertexFilenames.push_back(jconfig.resolvePath(fv["filename"].get<std::string>()));

  std::vector<std::shared_ptr<ConstraintPotentialEnergies::MultipleVertexPulling>> pullingEnergies;
  std::vector<ES::VXd> pullingTargets;
  std::vector<ES::VXd> pullingTargetRests;
  buildPullingConstraints(jconfig, fixedVertexFilenames, simulationRestPosition, K, pullingEnergies, pullingTargets, pullingTargetRests);

  ES::SpMatD M;
  libiglInterface::computeMassMatrix(surfaceMesh, M, 1, 1);
  M *= 100;

  ES::MXd V;
  ES::MXi F;
  Mesh::triMeshGeoToMatrices(surfaceMesh, V, F);
  const ES::SpMatD W = makeIdentityEmbedding(n3);

  IpcSimulationContext context;
  context.M = std::move(M);
  context.simulationRestPosition = std::move(simulationRestPosition);
  context.surfaceRestPositions = std::move(surfaceRestPositions);
  context.elasticParams = std::move(elasticParams);
  context.surfaceFromSimulationDispMap = W;
  context.simulationMeshOwner = simMesh;
  context.deformationModelManagerOwner = dmm;
  context.deformationModelAssemblerOwner = assembler;
  context.elasticEnergy = elasticEnergy;
  context.pullingEnergies = std::move(pullingEnergies);
  context.pullingTargets = std::move(pullingTargets);
  context.pullingTargetRests = std::move(pullingTargetRests);
  context.surfaceMesh = std::move(surfaceMesh);
  std::vector<bool> staticFlags;
  auto obstacles = parseExternalObjects(jconfig, 1.0, &staticFlags);
  const std::size_t obstacleCount = obstacles.size();
  context.collisionHandler =
    std::make_shared<Contact::IPC::EmbeddedSurfaceIPCPotentialEnergy>(
      V, F, context.surfaceFromSimulationDispMap, ipcParams, std::move(obstacles));
  context.contactBackend = makeIpcContactBackend();
  for (std::size_t i = 0; i < staticFlags.size(); ++i)
    if (staticFlags[i])
      context.collisionHandler->markObstacleStatic(static_cast<int32_t>(i));
  for (const ParsedFloorConfig &floorConfig : floorConfigs) {
    auto floorEnergy =
      std::make_shared<Contact::IPC::EmbeddedSurfaceFloorPotentialEnergy>(V, context.surfaceFromSimulationDispMap, floorConfig.params);
    context.floorPotentialEnergies.push_back(floorEnergy);
    context.floorMotionStates.push_back(floorConfig.motionState);
    context.extraGeneralImplicitForceModels.push_back(
      std::move(floorEnergy));
  }
  if (obstacleCount > 0)
    std::cout << ", obstacles=" << obstacleCount;
  return context;
}
}  // namespace pgo::RunIPCSim
