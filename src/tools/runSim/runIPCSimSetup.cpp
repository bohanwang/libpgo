#include "runIPCSimSetup.h"

#include "EigenSupport.h"
#include "basicIO.h"
#include "barycentricCoordinates.h"
#include "configFileJSON.h"
#include "deformationModelAssembler.h"
#include "deformationModelEnergy.h"
#include "deformationModelManager.h"
#include "embeddedSurfaceIPCPotentialEnergy.h"
#include "generateMassMatrix.h"
#include "geometryQuery.h"
#include "libiglInterface.h"
#include "multiVertexPullingSoftConstraints.h"
#include "pgoLogging.h"
#include "runSimFEMSetup.h"
#include "runSimVolumeMeshIO.h"
#include "simulationMesh.h"
#include "surfaceIPCCore.h"
#include "volumetricMesh.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace pgo::RunIPCSim
{
namespace
{
namespace ES = pgo::EigenSupport;

[[noreturn]] void throwConfigError(const std::string &message)
{
  throw std::invalid_argument(message);
}

void rejectIfPresent(const pgo::ConfigFileJSON &config, const char *field, const char *reason)
{
  if (config.exist(field))
    throwConfigError(std::string("runIPCSim phase1D does not support `") + field + "`: " + reason);
}

ES::SpMatD makeIdentityEmbedding(int n3)
{
  std::vector<ES::TripletD> triplets;
  triplets.reserve(n3);
  for (int i = 0; i < n3; ++i)
    triplets.emplace_back(i, i, 1.0);

  ES::SpMatD W(n3, n3);
  W.setFromTriplets(triplets.begin(), triplets.end());
  return W;
}

void validateZeroInitialDisplacement(const pgo::ConfigFileJSON &jconfig)
{
  const ES::V3d initialDisp = ES::Mp<ES::V3d>(jconfig.getValue<std::array<double, 3>>("init-disp", 1).data());
  if (initialDisp.cwiseAbs().maxCoeff() > 0.0)
    throwConfigError("runIPCSim phase1D only accepts zero `init-disp`.");
}

Contact::CIPC::SurfaceIPCCore::Parameters makeShellIPCParams(
  const pgo::ConfigFileJSON &jconfig, const pgo::Mesh::BoundingBox &surfaceBox)
{
  constexpr double kShellYoungsModulus = 1000000.0;
  constexpr double kShellThickness = 3e-3;

  Contact::CIPC::SurfaceIPCCore::Parameters ipcParams;
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

  return ipcParams;
}

Contact::CIPC::SurfaceIPCCore::Parameters makeVolumeIPCParams(const pgo::ConfigFileJSON &jconfig)
{
  const bool ipcHeuristic = jconfig.exist("ipc-heuristic") ? jconfig.getValue<bool>("ipc-heuristic", 1) : false;
  if (ipcHeuristic) {
    throwConfigError(
      "runIPCSim phase1D tet/cubic IPC currently requires explicit `ipc-dhat` and `ipc-kappa`; `ipc-heuristic=true` is shell-only.");
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
  return ipcParams;
}

void loadSurfaceMeshAndRestPositions(const std::string &surfaceMeshFilename, double scale,
  pgo::Mesh::TriMeshGeo &surfaceMesh, ES::VXd &surfaceRestPositions)
{
  if (!surfaceMesh.load(surfaceMeshFilename))
    throw std::runtime_error("Failed to load surface mesh: " + surfaceMeshFilename);

  for (int vi = 0; vi < surfaceMesh.numVertices(); ++vi)
    surfaceMesh.pos(vi) *= scale;

  surfaceRestPositions.resize(surfaceMesh.numVertices() * 3);
  for (int vi = 0; vi < surfaceMesh.numVertices(); ++vi)
    surfaceRestPositions.segment<3>(vi * 3) = surfaceMesh.pos(vi);
}

void buildPullingConstraints(const pgo::ConfigFileJSON &jconfig, const std::vector<std::string> &fixedVertexFilenames,
  const ES::VXd &simulationRestPosition, const ES::SpMatD &K,
  std::vector<std::shared_ptr<ConstraintPotentialEnergies::MultipleVertexPulling>> &pullingEnergies,
  std::vector<ES::VXd> &pullingTargets, std::vector<ES::VXd> &pullingTargetRests)
{
  int fixedVertexFileIndex = 0;
  for (const auto &fv : jconfig.handle()["fixed-vertices"]) {
    const std::string &filename = fixedVertexFilenames.at(fixedVertexFileIndex++);
    const std::array<double, 3> movement = fv["movement"].get<std::array<double, 3>>();
    const double attachmentCoeff = fv["coeff"].get<double>();

    std::vector<int> fixedVertices;
    if (BasicIO::read1DText(filename.c_str(), std::back_inserter(fixedVertices)) != 0) {
      throw std::runtime_error("Failed to read fixed vertex file: " + filename);
    }
    std::sort(fixedVertices.begin(), fixedVertices.end());

    ES::VXd tgtVertexPositions(fixedVertices.size() * 3);
    ES::VXd tgtVertexRests(fixedVertices.size() * 3);
    const ES::V3d movementVec(movement[0], movement[1], movement[2]);
    for (int vi = 0; vi < static_cast<int>(fixedVertices.size()); ++vi) {
      tgtVertexPositions.segment<3>(vi * 3) =
        simulationRestPosition.segment<3>(fixedVertices[vi] * 3) + movementVec;
      tgtVertexRests.segment<3>(vi * 3) =
        simulationRestPosition.segment<3>(fixedVertices[vi] * 3);
    }

    auto pullingEnergy = std::make_shared<ConstraintPotentialEnergies::MultipleVertexPulling>(
      K, simulationRestPosition.data(), static_cast<int>(fixedVertices.size()), fixedVertices.data(), tgtVertexPositions.data(), nullptr, 1);
    pullingEnergy->setCoeff(attachmentCoeff);
    pullingEnergies.push_back(pullingEnergy);
    pullingTargets.push_back(tgtVertexPositions);
    pullingTargetRests.push_back(tgtVertexRests);
  }
}

SolidDeformationModel::DeformationModelElasticMaterial parseVolumeElasticMaterial(const pgo::ConfigFileJSON &jconfig)
{
  const std::string material = jconfig.getString("elastic-material");
  if (material == "stable-neo")
    return SolidDeformationModel::DeformationModelElasticMaterial::STABLE_NEO;
  if (material == "stvk-vol")
    return SolidDeformationModel::DeformationModelElasticMaterial::STVK_VOL;
  if (material == "koiter-stvk") {
    throwConfigError(
      "runIPCSim phase1D tet/cubic only supports `elastic-material = stable-neo` or `stvk-vol`; `koiter-stvk` is shell-only.");
  }

  throwConfigError(
    "runIPCSim phase1D tet/cubic only supports `elastic-material = stable-neo` or `stvk-vol`.");
}
}  // namespace

IpcSimulationContext buildShellIpcSimulation(const pgo::ConfigFileJSON &jconfig)
{
  validateZeroInitialDisplacement(jconfig);
  rejectIfPresent(jconfig, "external-objects", "external contact is out of scope for phase1D.");

  if (jconfig.exist("tet-mesh") || jconfig.exist("cubic-mesh")) {
    throwConfigError("runIPCSim phase1D shell setup cannot consume tet/cubic mesh inputs.");
  }
  if (!jconfig.exist("fixed-vertices"))
    throwConfigError("Missing required field `fixed-vertices`.");

  const double scale = jconfig.getDouble("scale", 1);
  PGO_ALOG(std::abs(scale - 1.0) < 1e-6);

  const std::string material = jconfig.getString("elastic-material");
  if (material != "koiter-stvk")
    throwConfigError("runIPCSim phase1D shell path only supports `elastic-material = koiter-stvk`.");

  const std::string surfaceMeshFilename = jconfig.getResolvedPath("surface-mesh", 1);

  pgo::Mesh::TriMeshGeo surfaceMesh;
  ES::VXd surfaceRestPositions;
  loadSurfaceMeshAndRestPositions(surfaceMeshFilename, 1.0, surfaceMesh, surfaceRestPositions);
  const pgo::Mesh::BoundingBox surfaceBox(surfaceMesh.positions());
  const bool ipcHeuristic = jconfig.exist("ipc-heuristic") ? jconfig.getValue<bool>("ipc-heuristic", 1) : false;
  const Contact::CIPC::SurfaceIPCCore::Parameters ipcParams = makeShellIPCParams(jconfig, surfaceBox);

  std::cout << "runIPCSim phase1D shell IPC parameters: "
            << "ipc-heuristic=" << (ipcHeuristic ? "true" : "false") << ", "
            << "source=" << (ipcHeuristic ? "heuristic" : "config") << ", "
            << "ipc-dhat=" << ipcParams.dhat << ", "
            << "ipc-kappa=" << ipcParams.kappa << ", "
            << "eps_ee=" << ipcParams.eps_ee << ", "
            << "slackness=" << ipcParams.slackness << std::endl;

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
  context.surfaceFromSimulationDispMap = W;
  context.simulationMeshOwner = simMesh;
  context.deformationModelManagerOwner = dmm;
  context.deformationModelAssemblerOwner = assembler;
  context.elasticEnergy = elasticEnergy;
  context.pullingEnergies = std::move(pullingEnergies);
  context.pullingTargets = std::move(pullingTargets);
  context.pullingTargetRests = std::move(pullingTargetRests);
  context.surfaceMesh = std::move(surfaceMesh);
  context.collisionHandler =
    std::make_shared<Contact::CIPC::EmbeddedSurfaceIPCPotentialEnergy>(V, F, context.surfaceFromSimulationDispMap, ipcParams);
  return context;
}

IpcSimulationContext buildVolumeIpcSimulation(const pgo::ConfigFileJSON &jconfig)
{
  validateZeroInitialDisplacement(jconfig);
  rejectIfPresent(jconfig, "external-objects", "external contact is out of scope for phase1D.");

  if (!jconfig.exist("fixed-vertices"))
    throwConfigError("Missing required field `fixed-vertices`.");

  const double scale = jconfig.getDouble("scale", 1);
  const Contact::CIPC::SurfaceIPCCore::Parameters ipcParams = makeVolumeIPCParams(jconfig);
  const SolidDeformationModel::DeformationModelElasticMaterial elasticMat = parseVolumeElasticMaterial(jconfig);

  std::cout << "runIPCSim phase1D volume IPC parameters: "
            << "ipc-heuristic=false, "
            << "source=config, "
            << "ipc-dhat=" << ipcParams.dhat << ", "
            << "ipc-kappa=" << ipcParams.kappa << ", "
            << "eps_ee=" << ipcParams.eps_ee << ", "
            << "slackness=" << ipcParams.slackness << std::endl;

  const RunSim::ResolvedRunSimPaths resolvedPaths = RunSim::resolveRunSimPaths(jconfig);
  std::unique_ptr<VolumetricMeshes::VolumetricMesh> volumetricMesh =
    RunSim::loadValidatedVolumeMesh(RunSim::parseVolumeMeshInputConfig(jconfig), scale);

  pgo::Mesh::TriMeshGeo surfaceMesh;
  ES::VXd surfaceRestPositions;
  loadSurfaceMeshAndRestPositions(resolvedPaths.surfaceMeshFilename, scale, surfaceMesh, surfaceRestPositions);

  pgo::InterpolationCoordinates::BarycentricCoordinates bc(
    surfaceMesh.numVertices(), surfaceRestPositions.data(), volumetricMesh.get());
  ES::SpMatD W = bc.generateInterpolationMatrix();
  if (W.rows() != surfaceMesh.numVertices() * 3)
    throwConfigError("runIPCSim phase1D volume setup produced an embedding matrix with unexpected row count.");
  if (W.cols() != volumetricMesh->getNumVertices() * 3)
    throwConfigError("runIPCSim phase1D volume setup produced an embedding matrix with unexpected column count.");
  if (W.nonZeros() <= 0)
    throwConfigError("runIPCSim phase1D volume setup produced an empty embedding matrix.");

  ES::SpMatD M;
  VolumetricMeshes::GenerateMassMatrix::computeMassMatrix(volumetricMesh.get(), M, true);

  RunSim::InitializedVolumetricSimulation initialized =
    RunSim::initializeVolumetricSimulation(*volumetricMesh, elasticMat);
  if (W.cols() != initialized.restPosition.size())
    throwConfigError("runIPCSim phase1D volume setup produced an embedding matrix incompatible with simulation DOFs.");

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
  context.surfaceFromSimulationDispMap = std::move(W);
  context.simulationMeshOwner = initialized.simMesh;
  context.deformationModelManagerOwner = initialized.dmm;
  context.deformationModelAssemblerOwner = initialized.assembler;
  context.elasticEnergy = initialized.elasticEnergy;
  context.pullingEnergies = std::move(pullingEnergies);
  context.pullingTargets = std::move(pullingTargets);
  context.pullingTargetRests = std::move(pullingTargetRests);
  context.surfaceMesh = std::move(surfaceMesh);
  context.collisionHandler =
    std::make_shared<Contact::CIPC::EmbeddedSurfaceIPCPotentialEnergy>(V, F, context.surfaceFromSimulationDispMap, ipcParams);
  return context;
}
}  // namespace pgo::RunIPCSim
