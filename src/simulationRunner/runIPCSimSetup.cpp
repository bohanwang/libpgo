#include "runIPCSimSetup.h"

// IPC-specific construction remains separate from sampled contact.
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
#include "runIPCExternalObjects.h"
#include "runSimVolumeMeshIO.h"
#include "simulationMesh.h"
#include "ipc/core/surfaceIPCCore.h"
#include "volumetricMesh.h"

#include <fmt/format.h>

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

ES::V3d parseInitialDisplacement(const pgo::ConfigFileJSON &jconfig)
{
  const ES::V3d initialDisp = ES::Mp<ES::V3d>(jconfig.getValue<std::array<double, 3>>("init-disp", 1).data());
  if (!initialDisp.allFinite())
    throwConfigError("`init-disp` must contain exactly three finite values.");
  return initialDisp;
}

ES::VXd makeInitialDisplacement(int numSimulationVertices, const ES::V3d &initialDisp)
{
  ES::VXd initialU(numSimulationVertices * 3);
  for (int vi = 0; vi < numSimulationVertices; ++vi)
    initialU.segment<3>(3 * vi) = initialDisp;
  return initialU;
}

double parseScale(const pgo::ConfigFileJSON &jconfig)
{
  const double scale = jconfig.getDouble("scale", 1);
  if (!std::isfinite(scale) || scale <= 0.0)
    throwConfigError("`scale` must be finite and strictly positive.");
  return scale;
}

struct ShellMaterial
{
  std::string name;
  SolidDeformationModel::DeformationModelElasticMaterial type;
  std::vector<double> elementParams;  // elastic parameters shared by every element
  double membraneStiffness;           // scales the heuristic IPC kappa
  double thickness;
  double arealDensity;
};

double parseOptionalPositive(const pgo::ConfigFileJSON &jconfig, const char *key, double defaultValue)
{
  const double value = jconfig.getValue<double>(key, 0, defaultValue);
  if (!std::isfinite(value) || value <= 0.0)
    throwConfigError(std::string("`") + key + "` must be finite and strictly positive.");
  return value;
}

ShellMaterial parseShellMaterial(const pgo::ConfigFileJSON &jconfig)
{
  ShellMaterial material;
  material.name = jconfig.getString("elastic-material");
  material.thickness = parseOptionalPositive(jconfig, "shell-thickness", 3e-3);
  material.arealDensity = parseOptionalPositive(jconfig, "shell-areal-density", 100.0);

  if (material.name != "koiter-stvk")
    throwConfigError("runIPCSim phase1D shell path only supports `elastic-material = koiter-stvk`.");
  const double E = parseOptionalPositive(jconfig, "shell-youngs-modulus", 1000000.0);
  const double bendingE = parseOptionalPositive(jconfig, "shell-bending-youngs-modulus", E);
  material.type = SolidDeformationModel::DeformationModelElasticMaterial::KOITER_STVK;
  material.elementParams = { E, 0.4, bendingE, 0.4, material.thickness };
  material.membraneStiffness = E;
  return material;
}

// The barrier distance comes from exactly one of two fields: `ipc-dhat`, an
// absolute distance, or `ipc-dhat-relative`, a fraction of the deformable
// surface's bounding-box diagonal (measured after `scale`). The relative form
// keeps dhat below the mesh's own feature size on meshes of any size, which
// avoids spurious barrier pairs between neighboring surface primitives.
double parseConfiguredDhat(const pgo::ConfigFileJSON &jconfig, const pgo::Mesh::BoundingBox &surfaceBox, std::string &source)
{
  const bool hasAbsolute = jconfig.exist("ipc-dhat");
  const bool hasRelative = jconfig.exist("ipc-dhat-relative");
  if (hasAbsolute && hasRelative)
    throwConfigError("`ipc-dhat` and `ipc-dhat-relative` are mutually exclusive; specify one of them.");

  if (hasRelative) {
    const double fraction = jconfig.getDouble("ipc-dhat-relative", 1);
    if (!std::isfinite(fraction) || fraction <= 0.0)
      throwConfigError("`ipc-dhat-relative` must be finite and strictly positive.");
    const double diagonal = surfaceBox.sides().norm();
    if (!std::isfinite(diagonal) || diagonal <= 0.0)
      throwConfigError("`ipc-dhat-relative` requires a surface mesh with a non-degenerate bounding box.");
    source = fmt::format("relative ({} x bounding-box diagonal {})", fraction, diagonal);
    return fraction * diagonal;
  }

  if (!hasAbsolute)
    throwConfigError("Missing required field `ipc-dhat` (or `ipc-dhat-relative`).");
  source = "config";
  const double dhat = jconfig.getDouble("ipc-dhat", 1);
  if (!std::isfinite(dhat) || dhat <= 0.0)
    throwConfigError("`ipc-dhat` must be finite and strictly positive.");
  return dhat;
}

Contact::CIPC::SurfaceIPCCore::Parameters makeShellIPCParams(
  const pgo::ConfigFileJSON &jconfig, const pgo::Mesh::BoundingBox &surfaceBox, const ShellMaterial &material,
  std::string &dhatSource)
{
  Contact::CIPC::SurfaceIPCCore::Parameters ipcParams;
  const bool ipcHeuristic = jconfig.exist("ipc-heuristic") ? jconfig.getValue<bool>("ipc-heuristic", 1) : false;
  if (ipcHeuristic) {
    ipcParams.dhat = surfaceBox.sides().norm() * 1e-3;
    ipcParams.kappa = material.membraneStiffness * material.thickness;
    dhatSource = "heuristic";
  }
  else {
    ipcParams.dhat = parseConfiguredDhat(jconfig, surfaceBox, dhatSource);
    if (!jconfig.exist("ipc-kappa"))
      throwConfigError("Missing required field `ipc-kappa`.");
    ipcParams.kappa = jconfig.getDouble("ipc-kappa", 1);
  }
  ipcParams.eps_ee = 0.0;
  ipcParams.slackness = 1.0;
  ipcParams.frictionCoeff = jconfig.getValue<double>("contact-friction-coeff", 0, 0.0);
  ipcParams.frictionEpsV = jconfig.getValue<double>("ipc-friction-epsv", 0, 1e-3);

  return ipcParams;
}

Contact::CIPC::SurfaceIPCCore::Parameters makeVolumeIPCParams(
  const pgo::ConfigFileJSON &jconfig, const pgo::Mesh::BoundingBox &surfaceBox, std::string &dhatSource)
{
  const bool ipcHeuristic = jconfig.exist("ipc-heuristic") ? jconfig.getValue<bool>("ipc-heuristic", 1) : false;
  if (ipcHeuristic) {
    throwConfigError(
      "runIPCSim phase1D tet/cubic IPC currently requires explicit `ipc-dhat` (or `ipc-dhat-relative`) and `ipc-kappa`; `ipc-heuristic=true` is shell-only.");
  }

  Contact::CIPC::SurfaceIPCCore::Parameters ipcParams;
  ipcParams.dhat = parseConfiguredDhat(jconfig, surfaceBox, dhatSource);
  if (!jconfig.exist("ipc-kappa"))
    throwConfigError("Missing required field `ipc-kappa`.");
  ipcParams.kappa = jconfig.getDouble("ipc-kappa", 1);
  ipcParams.eps_ee = 0.0;
  ipcParams.slackness = 1.0;
  ipcParams.frictionCoeff = jconfig.getValue<double>("contact-friction-coeff", 0, 0.0);
  ipcParams.frictionEpsV = jconfig.getValue<double>("ipc-friction-epsv", 0, 1e-3);
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

// Entries without a `filename` (trajectory entries) get an empty string so the
// entry indices stay aligned with the config array.
std::vector<std::string> resolveFixedVertexFilenames(const pgo::ConfigFileJSON &jconfig)
{
  std::vector<std::string> filenames;
  for (const auto &fv : jconfig.handle()["fixed-vertices"]) {
    if (fv.is_object() && fv.contains("filename") && fv["filename"].is_string())
      filenames.push_back(jconfig.resolvePath(fv["filename"].get<std::string>()));
    else
      filenames.push_back(std::string());
  }
  return filenames;
}

void buildPullingConstraints(const pgo::ConfigFileJSON &jconfig, const std::vector<std::string> &fixedVertexFilenames,
  const ES::VXd &simulationRestPosition, const ES::SpMatD &K, double scale,
  std::vector<std::shared_ptr<ConstraintPotentialEnergies::MultipleVertexPulling>> &pullingEnergies,
  std::vector<ES::VXd> &pullingTargets, std::vector<ES::VXd> &pullingTargetRests,
  std::vector<TrajectoryPulling> &trajectoryPullings)
{
  const int numSimulationVertices = static_cast<int>(simulationRestPosition.size() / 3);
  int entryIndex = 0;
  for (const auto &fv : jconfig.handle()["fixed-vertices"]) {
    const std::string entryLabel = "fixed-vertices[" + std::to_string(entryIndex) + "]";
    const std::string &filename = fixedVertexFilenames.at(entryIndex++);
    if (!fv.is_object())
      throwConfigError(entryLabel + " must be an object.");
    const bool hasTrajectory = fv.contains("trajectory");
    const bool hasFile = fv.contains("filename") || fv.contains("movement");
    if (hasTrajectory && hasFile)
      throwConfigError(entryLabel + " must use either `filename` + `movement` or `trajectory`, not both.");
    if (!hasTrajectory && !hasFile)
      throwConfigError(entryLabel + " requires `filename` + `movement` or `trajectory`.");
    if (!fv.contains("coeff") || !fv["coeff"].is_number())
      throwConfigError(entryLabel + " requires a numeric `coeff`.");
    const double attachmentCoeff = fv["coeff"].get<double>();
    if (!std::isfinite(attachmentCoeff) || attachmentCoeff < 0.0)
      throwConfigError(entryLabel + ".coeff must be finite and non-negative.");

    if (hasTrajectory) {
      if (!fv["trajectory"].is_string())
        throwConfigError(entryLabel + ".trajectory must be a file path.");
      const std::string trajectoryPath = jconfig.resolvePath(fv["trajectory"].get<std::string>());
      RunSim::VertexTrajectory trajectory = RunSim::loadVertexTrajectory(trajectoryPath);
      for (RunSim::VertexTrajectoryTrack &track : trajectory.tracks) {
        for (int v : track.vertices) {
          if (v >= numSimulationVertices)
            throwConfigError(entryLabel + " trajectory track `" + track.name + "` references vertex " + std::to_string(v) +
              ", but the simulation mesh has " + std::to_string(numSimulationVertices) + " vertices.");
        }
        TrajectoryPulling pulling;
        pulling.space = trajectory.space;
        pulling.scale = scale;
        pulling.coeff = attachmentCoeff;
        pulling.restPositions.resize(3 * track.numVertices());
        for (int i = 0; i < track.numVertices(); ++i)
          pulling.restPositions.segment<3>(3 * i) = simulationRestPosition.segment<3>(3 * track.vertices[i]);
        pulling.track = std::move(track);
        const ES::VXd initialTarget = pulling.targetAt(0.0);
        pulling.energy = std::make_shared<ConstraintPotentialEnergies::MultipleVertexPulling>(
          K, simulationRestPosition.data(), pulling.track.numVertices(), pulling.track.vertices.data(),
          initialTarget.data(), nullptr, 1);
        pulling.energy->setCoeff(attachmentCoeff);
        std::cout << "runIPCSim trajectory attachment `" << pulling.track.name << "`: "
                  << pulling.track.numVertices() << " vertices, " << pulling.track.numKeys() << " keys over ["
                  << pulling.track.times.front() << ", " << pulling.track.times.back() << "] s, peak speed "
                  << pulling.track.peakSpeed() * scale << " per s, release "
                  << (pulling.track.releaseTime ? fmt::format("at {} s", *pulling.track.releaseTime) : std::string("never"))
                  << ", coeff=" << attachmentCoeff << std::endl;
        trajectoryPullings.push_back(std::move(pulling));
      }
      continue;
    }

    const std::array<double, 3> movement = fv["movement"].get<std::array<double, 3>>();

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

bool parseEnableMaterialMaxStep(const pgo::ConfigFileJSON &jconfig)
{
  return jconfig.exist("enable-material-max-step") ? jconfig.getValue<bool>("enable-material-max-step", 1) : true;
}

void rejectRemovedFloorConfig(const pgo::ConfigFileJSON &jconfig)
{
  if (jconfig.exist("use-floor") || jconfig.exist("floor-axis") ||
    jconfig.exist("floor-height") || jconfig.exist("floor-kappa"))
    throwConfigError("Analytic IPC floor fields were removed; use `ipc-external-objects` floor geometry instead.");
}
}  // namespace

ES::VXd TrajectoryPulling::targetAt(double time) const
{
  ES::VXd target = track.sampleAt(time) * scale;
  if (space == RunSim::VertexTrajectory::Space::Displacement)
    target += restPositions;
  return target;
}

IpcSimulationContext buildShellIpcSimulation(const pgo::ConfigFileJSON &jconfig)
{
  rejectRemovedFloorConfig(jconfig);
  const ES::V3d initialDisp = parseInitialDisplacement(jconfig);

  if (jconfig.exist("tet-mesh") || jconfig.exist("cubic-mesh")) {
    throwConfigError("runIPCSim phase1D shell setup cannot consume tet/cubic mesh inputs.");
  }
  if (!jconfig.exist("fixed-vertices"))
    throwConfigError("Missing required field `fixed-vertices`.");

  const double scale = parseScale(jconfig);

  const ShellMaterial shellMaterial = parseShellMaterial(jconfig);

  const std::string surfaceMeshFilename = jconfig.getResolvedPath("surface-mesh", 1);

  pgo::Mesh::TriMeshGeo surfaceMesh;
  ES::VXd surfaceRestPositions;
  loadSurfaceMeshAndRestPositions(surfaceMeshFilename, scale, surfaceMesh, surfaceRestPositions);
  const pgo::Mesh::BoundingBox surfaceBox(surfaceMesh.positions());
  const bool ipcHeuristic = jconfig.exist("ipc-heuristic") ? jconfig.getValue<bool>("ipc-heuristic", 1) : false;
  const bool enableMaterialMaxStep = parseEnableMaterialMaxStep(jconfig);
  const bool zeroRestCurvature = jconfig.exist("shell-zero-rest-curvature") ? jconfig.getValue<bool>("shell-zero-rest-curvature", 1) : false;
  std::string dhatSource;
  const Contact::CIPC::SurfaceIPCCore::Parameters ipcParams = makeShellIPCParams(jconfig, surfaceBox, shellMaterial, dhatSource);

  std::cout << "runIPCSim phase1D shell material: elastic-material=" << shellMaterial.name << ", element-params=[";
  for (std::size_t i = 0; i < shellMaterial.elementParams.size(); ++i)
    std::cout << (i ? ", " : "") << shellMaterial.elementParams[i];
  std::cout << "], shell-thickness=" << shellMaterial.thickness
            << ", shell-areal-density=" << shellMaterial.arealDensity
            << ", shell-zero-rest-curvature=" << (zeroRestCurvature ? "true" : "false") << std::endl;
  std::cout << "runIPCSim phase1D shell IPC parameters: "
            << "ipc-heuristic=" << (ipcHeuristic ? "true" : "false") << ", "
            << "source=" << dhatSource << ", "
            << "enable-material-max-step=" << (enableMaterialMaxStep ? "true" : "false") << ", "
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
  dmm->init(pgo::SolidDeformationModel::DeformationModelPlasticMaterial::SHELL_FF_DOF0, shellMaterial.type);
  dmm->setEnforceSPD(1);
  if (zeroRestCurvature)
    dmm->zeroShellRestCurvature();

  std::vector<double> elementWeights(simMesh->getNumElements(), 1.0);
  std::shared_ptr<SolidDeformationModel::DeformationModelAssembler> assembler =
    std::make_shared<SolidDeformationModel::DeformationModelAssembler>(dmm, elementWeights.data());

  const int n = simMesh->getNumVertices();
  const int n3 = n * 3;
  const int nele = simMesh->getNumElements();

  const int numElementParams = static_cast<int>(shellMaterial.elementParams.size());
  ES::VXd elasticParams(numElementParams * nele);
  for (int ei = 0; ei < nele; ++ei)
    elasticParams.segment(ei * numElementParams, numElementParams) =
      Eigen::Map<const ES::VXd>(shellMaterial.elementParams.data(), numElementParams);

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
  auto pullingHessianBase = std::make_shared<ES::SpMatD>();
  elasticEnergy->createHessian(*pullingHessianBase);
  elasticEnergy->hessian(zero, *pullingHessianBase);

  const std::vector<std::string> fixedVertexFilenames = resolveFixedVertexFilenames(jconfig);

  std::vector<std::shared_ptr<ConstraintPotentialEnergies::MultipleVertexPulling>> pullingEnergies;
  std::vector<ES::VXd> pullingTargets;
  std::vector<ES::VXd> pullingTargetRests;
  std::vector<TrajectoryPulling> trajectoryPullings;
  buildPullingConstraints(jconfig, fixedVertexFilenames, simulationRestPosition,
    *pullingHessianBase, scale, pullingEnergies, pullingTargets, pullingTargetRests, trajectoryPullings);

  // `initial-surface-mesh`: the same mesh in its state at t = 0, which becomes a
  // per-vertex initial displacement on top of the uniform `init-disp`.
  ES::VXd initialDisplacement = makeInitialDisplacement(n, initialDisp);
  if (jconfig.exist("initial-surface-mesh")) {
    const std::string initialMeshFilename = jconfig.getResolvedPath("initial-surface-mesh", 1);
    pgo::Mesh::TriMeshGeo initialMesh;
    if (!initialMesh.load(initialMeshFilename))
      throw std::runtime_error("Failed to load initial surface mesh: " + initialMeshFilename);
    if (initialMesh.numVertices() != surfaceMesh.numVertices() || initialMesh.numTriangles() != surfaceMesh.numTriangles()) {
      throwConfigError(fmt::format(
        "`initial-surface-mesh` must have the vertex and triangle counts of `surface-mesh` ({} vertices, {} triangles), got {} and {}.",
        surfaceMesh.numVertices(), surfaceMesh.numTriangles(), initialMesh.numVertices(), initialMesh.numTriangles()));
    }
    for (int fi = 0; fi < surfaceMesh.numTriangles(); ++fi) {
      for (int j = 0; j < 3; ++j) {
        if (initialMesh.tri(fi)[j] != surfaceMesh.tri(fi)[j])
          throwConfigError(fmt::format("`initial-surface-mesh` triangle {} differs from `surface-mesh`; the files must share topology.", fi));
      }
    }
    double maxInitialOffset = 0.0;
    for (int vi = 0; vi < n; ++vi) {
      if (!initialMesh.pos(vi).allFinite())
        throwConfigError("`initial-surface-mesh` contains a non-finite vertex.");
      const ES::V3d offset = scale * initialMesh.pos(vi) - surfaceRestPositions.segment<3>(3 * vi);
      initialDisplacement.segment<3>(3 * vi) += offset;
      maxInitialOffset = std::max(maxInitialOffset, offset.norm());
    }
    std::cout << "runIPCSim initial-surface-mesh: " << initialMeshFilename
              << ", max vertex offset from rest " << maxInitialOffset << std::endl;
  }

  ES::SpMatD M;
  libiglInterface::computeMassMatrix(surfaceMesh, M, 1, 1);
  M *= shellMaterial.arealDensity;

  ES::MXd V;
  ES::MXi F;
  Mesh::triMeshGeoToMatrices(surfaceMesh, V, F);
  const ES::SpMatD W = makeIdentityEmbedding(n3);
  IPCPreparedCollisionSurface preparedSurface = prepareIPCCollisionSurface(V, F, W, jconfig);
  const Contact::CIPC::IPCCollisionSurfaceData &collisionSurface = preparedSurface.surface;

  IpcSimulationContext context;
  context.M = std::move(M);
  context.simulationRestPosition = std::move(simulationRestPosition);
  context.initialDisplacement = std::move(initialDisplacement);
  context.surfaceRestPositions = std::move(surfaceRestPositions);
  context.surfaceFromSimulationDispMap = W;
  context.simulationMeshOwner = simMesh;
  context.deformationModelManagerOwner = dmm;
  context.deformationModelAssemblerOwner = assembler;
  context.elasticEnergy = elasticEnergy;
  context.pullingHessianBaseOwner = std::move(pullingHessianBase);
  context.pullingEnergies = std::move(pullingEnergies);
  context.pullingTargets = std::move(pullingTargets);
  context.pullingTargetRests = std::move(pullingTargetRests);
  context.trajectoryPullings = std::move(trajectoryPullings);
  context.surfaceMesh = std::move(surfaceMesh);
  context.collisionHandler =
    std::make_shared<Contact::CIPC::EmbeddedSurfaceIPCPotentialEnergy>(
      collisionSurface.vertices, collisionSurface.triangles, collisionSurface.displacementMap,
      collisionSurface.vertexIsDeformable, ipcParams);
  context.kinematicExternalObjects = std::move(preparedSurface.kinematicObjects);
  return context;
}

IpcSimulationContext buildVolumeIpcSimulation(const pgo::ConfigFileJSON &jconfig)
{
  rejectRemovedFloorConfig(jconfig);
  const ES::V3d initialDisp = parseInitialDisplacement(jconfig);

  if (!jconfig.exist("fixed-vertices"))
    throwConfigError("Missing required field `fixed-vertices`.");

  const double scale = parseScale(jconfig);
  const SolidDeformationModel::DeformationModelElasticMaterial elasticMat = parseVolumeElasticMaterial(jconfig);
  const bool enableMaterialMaxStep = parseEnableMaterialMaxStep(jconfig);

  const RunSim::ResolvedRunSimPaths resolvedPaths = RunSim::resolveRunSimPaths(jconfig);
  std::unique_ptr<VolumetricMeshes::VolumetricMesh> volumetricMesh =
    RunSim::loadValidatedVolumeMesh(RunSim::parseVolumeMeshInputConfig(jconfig), scale);

  pgo::Mesh::TriMeshGeo surfaceMesh;
  ES::VXd surfaceRestPositions;
  loadSurfaceMeshAndRestPositions(resolvedPaths.surfaceMeshFilename, scale, surfaceMesh, surfaceRestPositions);

  // dhat may be relative to the (scaled) surface bounding box, so the surface
  // mesh is loaded first.
  const pgo::Mesh::BoundingBox surfaceBox(surfaceMesh.positions());
  std::string dhatSource;
  const Contact::CIPC::SurfaceIPCCore::Parameters ipcParams = makeVolumeIPCParams(jconfig, surfaceBox, dhatSource);

  std::cout << "runIPCSim phase1D volume IPC parameters: "
            << "ipc-heuristic=false, "
            << "source=" << dhatSource << ", "
            << "enable-material-max-step=" << (enableMaterialMaxStep ? "true" : "false") << ", "
            << "ipc-dhat=" << ipcParams.dhat << ", "
            << "ipc-kappa=" << ipcParams.kappa << ", "
            << "eps_ee=" << ipcParams.eps_ee << ", "
            << "slackness=" << ipcParams.slackness << std::endl;

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
    RunSim::initializeVolumetricSimulation(*volumetricMesh, elasticMat,
      SolidDeformationModel::DeformationModelPlasticMaterial::VOLUMETRIC_DOF6,
      enableMaterialMaxStep);
  if (W.cols() != initialized.restPosition.size())
    throwConfigError("runIPCSim phase1D volume setup produced an embedding matrix incompatible with simulation DOFs.");

  ES::VXd zero = ES::VXd::Zero(initialized.restPosition.size());
  auto pullingHessianBase = std::make_shared<ES::SpMatD>();
  initialized.elasticEnergy->createHessian(*pullingHessianBase);
  initialized.elasticEnergy->hessian(zero, *pullingHessianBase);

  std::vector<std::shared_ptr<ConstraintPotentialEnergies::MultipleVertexPulling>> pullingEnergies;
  std::vector<ES::VXd> pullingTargets;
  std::vector<ES::VXd> pullingTargetRests;
  if (jconfig.exist("initial-surface-mesh"))
    throwConfigError("`initial-surface-mesh` is supported by shell simulations only; tet/cubic simulations use `init-disp` or a restart state.");

  std::vector<TrajectoryPulling> trajectoryPullings;
  buildPullingConstraints(jconfig, resolvedPaths.fixedVertexFilenames, initialized.restPosition,
    *pullingHessianBase, scale,
    pullingEnergies, pullingTargets, pullingTargetRests, trajectoryPullings);

  ES::MXd V;
  ES::MXi F;
  Mesh::triMeshGeoToMatrices(surfaceMesh, V, F);
  IPCPreparedCollisionSurface preparedSurface = prepareIPCCollisionSurface(V, F, W, jconfig);
  const Contact::CIPC::IPCCollisionSurfaceData &collisionSurface = preparedSurface.surface;

  IpcSimulationContext context;
  context.M = std::move(M);
  context.simulationRestPosition = std::move(initialized.restPosition);
  context.initialDisplacement = makeInitialDisplacement(
    static_cast<int>(context.simulationRestPosition.size() / 3), initialDisp);
  context.surfaceRestPositions = std::move(surfaceRestPositions);
  context.surfaceFromSimulationDispMap = std::move(W);
  context.simulationMeshOwner = initialized.simMesh;
  context.deformationModelManagerOwner = initialized.dmm;
  context.deformationModelAssemblerOwner = initialized.assembler;
  context.elasticEnergy = initialized.elasticEnergy;
  context.pullingHessianBaseOwner = std::move(pullingHessianBase);
  context.pullingEnergies = std::move(pullingEnergies);
  context.pullingTargets = std::move(pullingTargets);
  context.pullingTargetRests = std::move(pullingTargetRests);
  context.trajectoryPullings = std::move(trajectoryPullings);
  context.surfaceMesh = std::move(surfaceMesh);
  context.collisionHandler =
    std::make_shared<Contact::CIPC::EmbeddedSurfaceIPCPotentialEnergy>(
      collisionSurface.vertices, collisionSurface.triangles, collisionSurface.displacementMap,
      collisionSurface.vertexIsDeformable, ipcParams);
  context.kinematicExternalObjects = std::move(preparedSurface.kinematicObjects);
  return context;
}
}  // namespace pgo::RunIPCSim
