#include "runIPCSimSetup.h"

#include "EigenSupport.h"
#include "basicIO.h"
#include "barycentricCoordinates.h"
#include "configFileJSON.h"
#include "deformationModelAssembler.h"
#include "deformationModelEnergy.h"
#include "deformationModelManager.h"
#include "embeddedSurfaceFloorPotentialEnergy.h"
#include "embeddedSurfaceIPCPotentialEnergy.h"
#include "generateMassMatrix.h"
#include "geometryQuery.h"
#include "libiglInterface.h"
#include "multiVertexPullingSoftConstraints.h"
#include "pgoLogging.h"
#include "runSimFEMSetup.h"
#include "runSimVolumeMeshIO.h"
#include "simulationMesh.h"
#include "ipc/core/surfaceIPCCore.h"
#include "ipc/external/obstacleSurface.h"
#include "volumetricMesh.h"
#include "triMeshPseudoNormal.h"

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
using pgo::Contact::CIPC::FloorAxis;
using pgo::Contact::CIPC::FloorPenaltyParameters;
using pgo::Contact::CIPC::FloorSide;

[[noreturn]] void throwConfigError(const std::string &message)
{
  throw std::invalid_argument(message);
}

void rejectIfPresent(const pgo::ConfigFileJSON &config, const char *field, const char *reason)
{
  if (config.exist(field))
    throwConfigError(std::string("runIPCSim phase2does not support `") + field + "`: " + reason);
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
    throwConfigError("runIPCSim phase2only accepts zero `init-disp`.");
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

  if (jconfig.exist("ipc-dhat-external"))
    ipcParams.dhat_external = jconfig.getDouble("ipc-dhat-external", 1);

  return ipcParams;
}

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
      "runIPCSim phase2tet/cubic only supports `elastic-material = stable-neo` or `stvk-vol`; `koiter-stvk` is shell-only.");
  }

  throwConfigError(
    "runIPCSim phase2tet/cubic only supports `elastic-material = stable-neo` or `stvk-vol`.");
}

bool parseEnableMaterialMaxStep(const pgo::ConfigFileJSON &jconfig)
{
  return jconfig.exist("enable-material-max-step")
    ? jconfig.getValue<bool>("enable-material-max-step", 1)
    : true;
}

struct ParsedFloorConfig
{
  FloorPenaltyParameters params;
  IpcFloorMotionState motionState;
};

struct ParsedSurfacePressureForceConfig
{
  bool enabled = false;
  double pressure = 0.0;
  int rampSteps = 1;
};

FloorAxis parseFloorAxis(const std::string &axis)
{
  if (axis == "x")
    return FloorAxis::X;
  if (axis == "y")
    return FloorAxis::Y;
  if (axis == "z")
    return FloorAxis::Z;
  throwConfigError("floor `axis` must be one of: x, y, z.");
}

const char *floorAxisToString(FloorAxis axis)
{
  switch (axis) {
    case FloorAxis::X:
      return "x";
    case FloorAxis::Y:
      return "y";
    case FloorAxis::Z:
      return "z";
    default:
      return "invalid";
  }
}

FloorSide parseFloorSide(const std::string &side)
{
  if (side == "lower")
    return FloorSide::LOWER;
  if (side == "upper")
    return FloorSide::UPPER;
  throwConfigError("floor `side` must be `lower` or `upper`.");
}

const char *floorSideToString(FloorSide side)
{
  switch (side) {
    case FloorSide::LOWER:
      return "lower";
    case FloorSide::UPPER:
      return "upper";
    default:
      return "invalid";
  }
}

double floorHeightAtFrame(const IpcFloorMotionState &motion, int frame)
{
  if (!motion.hasMotion)
    return motion.heightStart;

  if (frame <= motion.frameStart)
    return motion.heightStart;
  if (frame >= motion.frameEnd)
    return motion.heightEnd;

  const double denom = static_cast<double>(motion.frameEnd - motion.frameStart);
  const double alpha = denom > 0.0 ? static_cast<double>(frame - motion.frameStart) / denom : 1.0;
  return motion.heightStart * (1.0 - alpha) + motion.heightEnd * alpha;
}

std::vector<ParsedFloorConfig> parseFloorsConfig(const pgo::ConfigFileJSON &jconfig)
{
  for (const char *legacyField : { "use-floor", "floor-axis", "floor-height", "floor-kappa" }) {
    if (jconfig.exist(legacyField))
      throwConfigError(std::string("runIPCSim floor config field `") + legacyField + "` has been replaced by `floors[]`.");
  }

  std::vector<ParsedFloorConfig> floors;
  if (!jconfig.exist("floors"))
    return floors;

  const auto &floorsJson = jconfig.handle()["floors"];
  if (!floorsJson.is_array())
    throwConfigError("`floors` must be a JSON array.");

  floors.reserve(floorsJson.size());
  for (std::size_t floorIndex = 0; floorIndex < floorsJson.size(); ++floorIndex) {
    const auto &floorJson = floorsJson.at(floorIndex);
    if (!floorJson.is_object())
      throwConfigError("Each `floors[]` entry must be a JSON object.");

    if (!floorJson.contains("axis"))
      throwConfigError("Missing required field `floors[].axis`.");
    if (!floorJson.contains("kappa"))
      throwConfigError("Missing required field `floors[].kappa`.");

    const bool hasHeight = floorJson.contains("height");
    const bool hasMotion = floorJson.contains("motion");
    if (hasHeight == hasMotion)
      throwConfigError("Each `floors[]` entry must provide exactly one of `height` or `motion`.");

    ParsedFloorConfig floorConfig;
    floorConfig.params.floorAxis = parseFloorAxis(floorJson.at("axis").get<std::string>());
    floorConfig.params.floorSide = floorJson.contains("side") ? parseFloorSide(floorJson.at("side").get<std::string>()) : FloorSide::LOWER;
    floorConfig.params.floorKappa = floorJson.at("kappa").get<double>();
    if (!std::isfinite(floorConfig.params.floorKappa))
      throwConfigError("`floors[].kappa` must be finite.");

    if (hasHeight) {
      floorConfig.params.floorHeight = floorJson.at("height").get<double>();
      if (!std::isfinite(floorConfig.params.floorHeight))
        throwConfigError("`floors[].height` must be finite.");
      floorConfig.motionState.hasMotion = false;
      floorConfig.motionState.heightStart = floorConfig.params.floorHeight;
      floorConfig.motionState.heightEnd = floorConfig.params.floorHeight;
    }
    else {
      const auto &motionJson = floorJson.at("motion");
      if (!motionJson.is_object())
        throwConfigError("`floors[].motion` must be a JSON object.");
      for (const char *field : { "height-start", "height-end", "frame-start", "frame-end" }) {
        if (!motionJson.contains(field))
          throwConfigError(std::string("Missing required field `floors[].motion.") + field + "`.");
      }
      floorConfig.motionState.hasMotion = true;
      floorConfig.motionState.heightStart = motionJson.at("height-start").get<double>();
      floorConfig.motionState.heightEnd = motionJson.at("height-end").get<double>();
      floorConfig.motionState.frameStart = motionJson.at("frame-start").get<int>();
      floorConfig.motionState.frameEnd = motionJson.at("frame-end").get<int>();
      if (!std::isfinite(floorConfig.motionState.heightStart) || !std::isfinite(floorConfig.motionState.heightEnd))
        throwConfigError("`floors[].motion` heights must be finite.");
      if (floorConfig.motionState.frameEnd < floorConfig.motionState.frameStart)
        throwConfigError("`floors[].motion.frame-end` must be greater than or equal to `frame-start`.");
      floorConfig.params.floorHeight = floorHeightAtFrame(floorConfig.motionState, 0);
    }

    floors.push_back(floorConfig);
  }

  return floors;
}

ParsedSurfacePressureForceConfig parseSurfacePressureForceConfig(const pgo::ConfigFileJSON &jconfig)
{
  ParsedSurfacePressureForceConfig pressureConfig;
  if (!jconfig.exist("surface-pressure-force"))
    return pressureConfig;

  const auto &pressureJson = jconfig.handle()["surface-pressure-force"];
  if (!pressureJson.is_object())
    throwConfigError("`surface-pressure-force` must be a JSON object.");

  pressureConfig.enabled = pressureJson.contains("enabled")
    ? pressureJson.at("enabled").get<bool>()
    : false;
  if (!pressureConfig.enabled)
    return pressureConfig;

  if (!pressureJson.contains("pressure"))
    throwConfigError("Missing required field `surface-pressure-force.pressure` when enabled.");

  pressureConfig.pressure = pressureJson.at("pressure").get<double>();
  if (pressureJson.contains("ramp-steps"))
    pressureConfig.rampSteps = pressureJson.at("ramp-steps").get<int>();

  if (!std::isfinite(pressureConfig.pressure))
    throwConfigError("`surface-pressure-force.pressure` must be finite.");
  if (pressureConfig.rampSteps <= 0)
    throwConfigError("`surface-pressure-force.ramp-steps` must be positive.");

  return pressureConfig;
}

ES::VXd computeSurfacePressureSimulationForce(const pgo::Mesh::TriMeshGeo &surfaceMesh,
  const ES::SpMatD &surfaceFromSimulationDispMap,
  const ParsedSurfacePressureForceConfig &pressureConfig, int simulationDofCount)
{
  if (!pressureConfig.enabled)
    return {};

  std::vector<double> vertexAreas(surfaceMesh.numVertices(), 0.0);
  surfaceMesh.ref().computeVertexSurfaceAreas(vertexAreas.data());

  pgo::Mesh::TriMeshPseudoNormal meshNormal(surfaceMesh);

  ES::VXd surfaceForce = ES::VXd::Zero(surfaceMesh.numVertices() * 3);
  for (int vi = 0; vi < surfaceMesh.numVertices(); ++vi) {
    surfaceForce.segment<3>(vi * 3) = meshNormal.vtxNormal(vi) * pressureConfig.pressure * vertexAreas[vi] * -1.0;
  }

  ES::VXd simulationForce(surfaceFromSimulationDispMap.cols());
  ES::mv(surfaceFromSimulationDispMap, surfaceForce, simulationForce, 1);
  if (simulationForce.size() != simulationDofCount)
    throwConfigError("surface-pressure-force projected force has an unexpected simulation DOF count.");

  return simulationForce;
}
std::vector<Contact::CIPC::ObstacleSurface> parseExternalObjects(
  const pgo::ConfigFileJSON &jconfig, double scale)
{
  std::vector<Contact::CIPC::ObstacleSurface> obstacles;
  if (!jconfig.exist("external-objects"))
    return obstacles;

  const auto &extObjs = jconfig.handle()["external-objects"];
  if (!extObjs.is_array())
    throwConfigError("`external-objects` must be a JSON array.");

  obstacles.reserve(extObjs.size());
  for (std::size_t oi = 0; oi < extObjs.size(); ++oi) {
    const auto &objJson = extObjs.at(oi);
    if (!objJson.is_object())
      throwConfigError("Each `external-objects[]` entry must be a JSON object.");

    if (!objJson.contains("filename"))
      throwConfigError("Missing required field `external-objects[].filename`.");
    if (!objJson.contains("movement"))
      throwConfigError("Missing required field `external-objects[].movement`.");

    const std::string filename = jconfig.resolvePath(objJson["filename"].get<std::string>());
    const std::array<double, 3> movementArr = objJson["movement"].get<std::array<double, 3>>();
    const double objScale = objJson.contains("scale") ? objJson["scale"].get<double>() : scale;

    // Load obstacle mesh
    pgo::Mesh::TriMeshGeo obsMesh;
    if (!obsMesh.load(filename))
      throw std::runtime_error("Failed to load external object mesh: " + filename);
    for (int vi = 0; vi < obsMesh.numVertices(); ++vi)
      obsMesh.pos(vi) *= objScale;

    ES::MXd V(obsMesh.numVertices(), 3);
    ES::MXi F(obsMesh.numTriangles(), 3);
    for (int vi = 0; vi < obsMesh.numVertices(); ++vi)
      V.row(vi) = obsMesh.pos(vi).transpose();
    for (int fi = 0; fi < obsMesh.numTriangles(); ++fi)
      F.row(fi) = obsMesh.tri(fi).transpose();

    ES::VXd restFlat(V.rows() * 3);
    for (int vi = 0; vi < V.rows(); ++vi)
      restFlat.segment<3>(vi * 3) = V.row(vi).transpose();

    ES::V3d velocity(movementArr[0], movementArr[1], movementArr[2]);
    auto sampler = Contact::CIPC::makeLinearTrajectorySampler(restFlat, velocity);

    obstacles.emplace_back(std::move(V), std::move(F), std::move(sampler));
  }

  return obstacles;
}
}  // namespace

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
  const Contact::CIPC::SurfaceIPCCore::Parameters ipcParams = makeShellIPCParams(jconfig, surfaceBox);
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
  auto obstacles = parseExternalObjects(jconfig, 1.0);
  const std::size_t obstacleCount = obstacles.size();
  context.collisionHandler =
    std::make_shared<Contact::CIPC::EmbeddedSurfaceIPCPotentialEnergy>(
      V, F, context.surfaceFromSimulationDispMap, ipcParams, std::move(obstacles));
  for (const ParsedFloorConfig &floorConfig : floorConfigs) {
    auto floorEnergy =
      std::make_shared<Contact::CIPC::EmbeddedSurfaceFloorPotentialEnergy>(V, context.surfaceFromSimulationDispMap, floorConfig.params);
    context.floorPotentialEnergies.push_back(floorEnergy);
    context.floorMotionStates.push_back(floorConfig.motionState);
    context.extraGeneralImplicitForceModels.push_back(
      std::move(floorEnergy));
  }
  if (obstacleCount > 0)
    std::cout << ", obstacles=" << obstacleCount;
  return context;
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
  auto obstacles = parseExternalObjects(jconfig, scale);
  const std::size_t obstacleCount = obstacles.size();
  context.collisionHandler =
    std::make_shared<Contact::CIPC::EmbeddedSurfaceIPCPotentialEnergy>(
      V, F, context.surfaceFromSimulationDispMap, ipcParams, std::move(obstacles));
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
