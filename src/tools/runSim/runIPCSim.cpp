#include "configFileJSON.h"
#include "EigenSupport.h"
#include "deformationModelAssembler.h"
#include "deformationModelEnergy.h"
#include "embeddedSurfaceIPCPotentialEnergy.h"
#include "embeddedSurfaceFloorPotentialEnergy.h"
#include "implicitBackwardEulerTimeIntegrator.h"
#include "implicitBackwardEulerTimeIntegratorHelper.h"
#include "initPredicates.h"
#include "multiVertexPullingSoftConstraints.h"
#include "pgoLogging.h"
#include "runIPCSimSetup.h"
#include "runSimCliLogging.h"
#include "scopedProfileSection.h"
#include "simulationMesh.h"

#include <argparse/argparse.hpp>
#include <fmt/format.h>
#include <nlohmann/json.hpp>
#include <tbb/global_control.h>

#include <array>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

namespace
{
namespace ES = pgo::EigenSupport;

struct OutputDirectories
{
  std::filesystem::path root;
  std::filesystem::path states;
  std::filesystem::path surface;
  std::filesystem::path stress;
};

const double *dataOrNull(const ES::VXd &values)
{
  return values.size() > 0 ? values.data() : nullptr;
}

OutputDirectories makeOutputDirectories(const std::filesystem::path &outputFolder)
{
  return {
    outputFolder,
    outputFolder / "states",
    outputFolder / "surface",
    outputFolder / "stress",
  };
}

std::filesystem::path framePath(const std::filesystem::path &dir, const char *prefix, int frame, const char *extension)
{
  return dir / fmt::format("{}{:04d}{}", prefix, frame, extension);
}

void createOutputSubdirectories(const OutputDirectories &outputDirs)
{
  std::error_code ec;
  std::filesystem::create_directories(outputDirs.root, ec);
  if (ec)
    throw std::runtime_error("Failed to create output folder `" + outputDirs.root.string() + "`: " + ec.message());

  for (const auto &dir : { outputDirs.states, outputDirs.surface, outputDirs.stress }) {
    std::filesystem::create_directories(dir, ec);
    if (ec)
      throw std::runtime_error("Failed to create output subfolder `" + dir.string() + "`: " + ec.message());
  }
}

const char *vonMisesStressLocation(const pgo::SolidDeformationModel::SimulationMesh &mesh)
{
  using pgo::SolidDeformationModel::SimulationMeshType;
  switch (mesh.getElementType()) {
    case SimulationMeshType::TET:
      return "tet_element";
    case SimulationMeshType::CUBIC:
      return "cubic_element";
    default:
      return "element";
  }
}

void logRunIPCSimMaxStepSummary(
  const std::shared_ptr<pgo::SolidDeformationModel::DeformationModelEnergy> &elasticEnergy,
  const std::shared_ptr<pgo::Contact::CIPC::EmbeddedSurfaceIPCPotentialEnergy> &collisionHandler,
  const std::shared_ptr<pgo::Simulation::ImplicitBackwardEulerTimeIntegrator> &integrator)
{
  auto logger = pgo::Logging::lgr();
  if (!logger)
    return;

  (void)elasticEnergy;
  (void)collisionHandler;
  const pgo::NonlinearOptimization::SolveDiagnostics &summary = integrator->getLastSolveDiagnostics();

  if (logger->should_log(spdlog::level::info)) {
    SPDLOG_LOGGER_INFO(logger,
      "runIPCSim max-step summary: materialClampCount={} contactClampCount={} minMaterialFeasibleAlphaThisSolve={} minContactFeasibleAlphaThisSolve={} minFeasibleAlphaThisSolve={} minLineSearchAlphaThisSolve={} minEffectiveAlphaThisSolve={}",
      summary.materialClampCount, summary.contactClampCount,
      summary.minMaterialFeasibleAlpha,
      summary.minContactFeasibleAlpha,
      summary.minFeasibleAlpha,
      summary.minLineSearchAlpha,
      summary.minEffectiveAlpha);
  }
}

std::filesystem::path resolveRunIPCSimLogPath(const std::filesystem::path &outputFolder)
{
  return outputFolder / "runIPCSim.log";
}

void clearOutputDirectory(const std::filesystem::path &outputFolder)
{
  std::error_code ec;
  std::filesystem::remove_all(outputFolder, ec);
  if (ec)
    throw std::runtime_error("Failed to clear output folder `" + outputFolder.string() + "`: " + ec.message());

  std::filesystem::create_directories(outputFolder, ec);
  if (ec)
    throw std::runtime_error("Failed to create output folder `" + outputFolder.string() + "`: " + ec.message());
}

void writeVonMisesStressJson(const OutputDirectories &outputDirs, int frame, double timestep,
  const pgo::RunIPCSim::IpcSimulationContext &context, const ES::VXd &displacement)
{
  if (!context.deformationModelAssemblerOwner || !context.simulationMeshOwner)
    throw std::runtime_error("runIPCSim cannot output von Mises stresses without a simulation mesh and assembler.");

  const int elementCount = context.simulationMeshOwner->getNumElements();
  std::vector<double> elementStresses(elementCount, 0.0);
  const ES::VXd absolutePositions = context.simulationRestPosition + displacement;
  context.deformationModelAssemblerOwner->computeVonMisesStresses(
    absolutePositions.data(),
    dataOrNull(context.plasticParams),
    dataOrNull(context.elasticParams),
    elementStresses.data());

  nlohmann::json stressJson;
  stressJson["frame"] = frame;
  stressJson["time"] = static_cast<double>(frame) * timestep;
  stressJson["stress_type"] = "von_mises";
  stressJson["location"] = vonMisesStressLocation(*context.simulationMeshOwner);
  stressJson["values"] = elementStresses;

  const std::filesystem::path outputPath = framePath(outputDirs.stress, "von_mises", frame, ".json");
  std::ofstream out(outputPath);
  if (!out.is_open())
    throw std::runtime_error("Failed to write von Mises stress JSON: " + outputPath.string());
  out << stressJson.dump(2) << '\n';
}

void logProfileSummary()
{
  auto logger = pgo::Logging::lgr();
  if (!logger || !logger->should_log(spdlog::level::info))
    return;

  const std::vector<pgo::Profiling::ProfileStat> stats = pgo::Profiling::snapshotProfileStatistics();
  SPDLOG_LOGGER_INFO(logger, "runIPCSim profiling summary:");
  for (const pgo::Profiling::ProfileStat &stat : stats) {
    SPDLOG_LOGGER_INFO(logger,
      "profile name={} callCount={} totalSeconds={} maxSeconds={}",
      stat.name, stat.callCount, stat.totalSeconds, stat.maxSeconds);
  }
}

double floorHeightAtFrame(const pgo::RunIPCSim::IpcFloorMotionState &motion, int frame)
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
}

int main(int argc, char *argv[])
{
  using namespace pgo;
  namespace ES = EigenSupport;

  argparse::ArgumentParser program("Run IPC Simulation");
  program.add_argument("config")
    .help("Config File")
    .required();
  program.add_argument("--log")
    .help("Write command-line output to a .log file next to the config file")
    .default_value(false)
    .implicit_value(true);

  try {
    program.parse_args(argc, argv);
  }
  catch (const std::exception &err) {
    std::cerr << err.what() << std::endl;
    std::cerr << program;
    return 1;
  }

  tbb::global_control c(tbb::global_control::max_allowed_parallelism, 64);

  const std::string configFilename = program.get<std::string>("config");
  const bool enableCliLog = program.get<bool>("--log");
  std::unique_ptr<RunSim::ScopedRunSimCliLogRedirect> logRedirect;
  bool enableProfiling = false;

  try {
    ConfigFileJSON jconfig;
    if (jconfig.open(configFilename.c_str()) != true)
      return 0;

    const ES::V3d extAcc = ES::Mp<ES::V3d>(jconfig.getValue<std::array<double, 3>>("g", 1).data());
    const ES::V3d initialVel = ES::Mp<ES::V3d>(jconfig.getValue<std::array<double, 3>>("init-vel", 1).data());
    const double timestep = jconfig.getDouble("timestep", 1);
    const double scale = jconfig.getDouble("scale", 1);
    const double solverEps = jconfig.getDouble("solver-eps", 1);
    const int solverMaxIter = jconfig.getInt("solver-max-iter", 1);
    const std::array<double, 2> dampingParams = jconfig.getValue<std::array<double, 2>>("damping-params", 1);
    const int numSimSteps = jconfig.getInt("num-timestep", 1);
    const int frameGap = jconfig.getInt("dump-interval", 1);
    const std::string simType = jconfig.getString("sim-type");
    if (simType != "dynamic")
      throw std::invalid_argument("runIPCSim phase1D only supports `sim-type = dynamic`.");
    const std::filesystem::path outputFolder = jconfig.getResolvedPath("output", 1);
    const OutputDirectories outputDirs = makeOutputDirectories(outputFolder);
    const bool restartFromU = jconfig.exist("restart-from-u") ? jconfig.getValue<bool>("restart-from-u", 1) : false;
    const bool dumpDeformEveryFrame = jconfig.exist("dump_deform_every_frame")
      ? jconfig.getValue<bool>("dump_deform_every_frame", 1)
      : false;
    const bool outputVonMises = jconfig.exist("output-von-mises") ? jconfig.getValue<bool>("output-von-mises", 1) : false;
    enableProfiling = jconfig.exist("profiling") ? jconfig.getValue<bool>("profiling", 1) : false;

    if (restartFromU) {
      createOutputSubdirectories(outputDirs);
    }
    else {
      clearOutputDirectory(outputFolder);
      createOutputSubdirectories(outputDirs);
    }

    if (enableCliLog) {
      const std::filesystem::path logPath = resolveRunIPCSimLogPath(outputFolder);
      logRedirect = std::make_unique<RunSim::ScopedRunSimCliLogRedirect>(logPath.string());
    }
    pgo::Logging::init(nullptr, RunSim::resolveConfiguredLogLevel(jconfig));

    if (!restartFromU) {
      std::cout << "restart-from-u=false; clearing output folder " << outputFolder << "." << std::endl;
    }
    if (enableProfiling) {
      pgo::Profiling::setProfilingEnabled(true);
      pgo::Profiling::resetProfileStatistics();
    }
    pgo::Mesh::initPredicates();

    const bool hasTetMesh = jconfig.exist("tet-mesh");
    const bool hasCubicMesh = jconfig.exist("cubic-mesh");
    const bool useVolumePath = hasTetMesh || hasCubicMesh;
    if (outputVonMises && !useVolumePath)
      throw std::invalid_argument("runIPCSim `output-von-mises` requires `tet-mesh` or `cubic-mesh`.");

    RunIPCSim::IpcSimulationContext context = useVolumePath
      ? RunIPCSim::buildVolumeIpcSimulation(jconfig)
      : RunIPCSim::buildShellIpcSimulation(jconfig);

    const int n3 = static_cast<int>(context.simulationRestPosition.size());
    const int surfn3 = static_cast<int>(context.surfaceRestPositions.size());
    const int n = n3 / 3;

    ES::VXd g(n3);
    for (int vi = 0; vi < n; ++vi)
      g.segment<3>(vi * 3) = extAcc;

    ES::VXd gravityForce(n3);
    ES::mv(context.M, g, gravityForce);
    ES::VXd fext = gravityForce;

    std::shared_ptr<Simulation::ImplicitBackwardEulerTimeIntegrator> intg =
      std::make_shared<Simulation::ImplicitBackwardEulerTimeIntegrator>(context.M, context.elasticEnergy,
        dampingParams[0], dampingParams[1], timestep, solverMaxIter, solverEps);

    for (auto &pullingEnergy : context.pullingEnergies)
      intg->addImplicitForceModel(pullingEnergy, 0, 0);

    intg->setExternalForce(fext.data());

    ES::VXd u = ES::VXd::Zero(n3);
    ES::VXd uvel = ES::VXd::Zero(n3);
    ES::VXd uacc = ES::VXd::Zero(n3);
    ES::VXd usurf = ES::VXd::Zero(surfn3);

    for (int i = 0; i < n; ++i)
      uvel.segment<3>(i * 3) = initialVel;

    int frameStart = -1;
    if (restartFromU) {
      for (int framei = numSimSteps - 1; framei >= 0; --framei) {
        const std::filesystem::path deformFilename = framePath(outputDirs.states, "deform", framei, ".u");
        if (!std::filesystem::exists(deformFilename))
          continue;

        ES::MXd uMat(n3, 3);
        if (ES::readMatrix(deformFilename.string().c_str(), uMat) == 0) {
          frameStart = framei;
          u.noalias() = uMat.col(0);
          uvel.noalias() = uMat.col(1);
          uacc.noalias() = uMat.col(2);
          std::cout << "Restarting from frame " << framei << std::endl;
          break;
        }
      }
    }

    if (frameStart < 0 && restartFromU)
      std::cout << "No restart state found in " << outputDirs.states << ". Starting from frame 0." << std::endl;
    else if (frameStart < 0)
      std::cout << "Starting from frame 0." << std::endl;

    const double ratioDenom = numSimSteps > 1 ? static_cast<double>(numSimSteps - 1) : 1.0;

    bool executedStep = false;
    for (int framei = frameStart + 1; framei < numSimSteps; ++framei) {
      intg->clearGeneralImplicitForceModel();

      const double ratio = numSimSteps > 1 ? static_cast<double>(framei) / ratioDenom : 1.0;
      for (std::size_t pi = 0; pi < context.pullingEnergies.size(); ++pi) {
        const ES::VXd curTgt = context.pullingTargetRests[pi] * (1.0 - ratio) + context.pullingTargets[pi] * ratio;
        context.pullingEnergies[pi]->setTargetPos(curTgt.data());
        std::cout << "Frame " << framei << ", attachment " << pi << " target: " << curTgt.transpose().head(3) << std::endl;
      }

      intg->addGeneralImplicitForceModel(context.collisionHandler, 0, 0);
      for (std::size_t fi = 0; fi < context.floorPotentialEnergies.size(); ++fi) {
        context.floorPotentialEnergies[fi]->setFloorHeight(floorHeightAtFrame(context.floorMotionStates[fi], framei));
      }
      for (const auto &forceModel : context.extraGeneralImplicitForceModels)
        intg->addGeneralImplicitForceModel(forceModel, 0, 0);
      if (context.surfacePressureForceEnabled) {
        const double ramp = std::min(1.0, static_cast<double>(framei + 1) / static_cast<double>(context.surfacePressureRampSteps));
        fext.noalias() = gravityForce + ramp * context.surfacePressureSimulationForce;
        intg->setExternalForce(fext.data());
      }
      intg->setqState(u, uvel, uacc);
      intg->doTimestep(1, 3, 1);
      executedStep = true;
      intg->getq(u);
      intg->getqvel(uvel);
      intg->getqacc(uacc);
      logRunIPCSimMaxStepSummary(context.elasticEnergy, context.collisionHandler, intg);

      const bool dumpDeformThisFrame = dumpDeformEveryFrame || (framei % frameGap == 0);
      if (dumpDeformThisFrame) {
        ES::MXd uMat(n3, 3);
        uMat.col(0) = u;
        uMat.col(1) = uvel;
        uMat.col(2) = uacc;
        ES::writeMatrix(framePath(outputDirs.states, "deform", framei, ".u").string().c_str(), uMat);
      }

      if (outputVonMises) {
        writeVonMisesStressJson(outputDirs, framei, timestep, context, u);
      }

      if (framei % frameGap == 0) {
        Mesh::TriMeshGeo mesh = context.surfaceMesh;
        ES::mv(context.surfaceFromSimulationDispMap, u, usurf);
        const ES::VXd psurf = context.surfaceRestPositions + usurf;
        for (int vi = 0; vi < mesh.numVertices(); ++vi)
          mesh.pos(vi) = psurf.segment<3>(vi * 3) / scale;
        mesh.save(framePath(outputDirs.surface, "ret", framei / frameGap, ".obj").string());
      }
    }

    if (!executedStep)
      logRunIPCSimMaxStepSummary(context.elasticEnergy, context.collisionHandler, intg);
    if (enableProfiling) {
      logProfileSummary();
      pgo::Profiling::setProfilingEnabled(false);
      pgo::Profiling::resetProfileStatistics();
    }
  }
  catch (const std::exception &err) {
    SPDLOG_LOGGER_ERROR(Logging::lgr(), "{}", err.what());
    if (enableProfiling) {
      pgo::Profiling::setProfilingEnabled(false);
      pgo::Profiling::resetProfileStatistics();
    }
    return 1;
  }

  return 0;
}
