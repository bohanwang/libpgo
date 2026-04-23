#include "configFileJSON.h"
#include "EigenSupport.h"
#include "deformationModelEnergy.h"
#include "embeddedSurfaceIPCPotentialEnergy.h"
#include "implicitBackwardEulerTimeIntegrator.h"
#include "implicitBackwardEulerTimeIntegratorHelper.h"
#include "initPredicates.h"
#include "multiVertexPullingSoftConstraints.h"
#include "pgoLogging.h"
#include "runIPCSimSetup.h"
#include "runSimCliLogging.h"

#include <argparse/argparse.hpp>
#include <fmt/format.h>

#include <array>
#include <filesystem>
#include <iostream>
#include <memory>
#include <stdexcept>

namespace
{
namespace ES = pgo::EigenSupport;

struct SolveMaxStepSummary
{
  double minFeasibleAlphaThisSolve = 1.0;
  double minLineSearchAlphaThisSolve = 1.0;
  double minEffectiveAlphaThisSolve = 1.0;
};

SolveMaxStepSummary currentSolveMaxStepSummary(const std::shared_ptr<pgo::Simulation::ImplicitBackwardEulerTimeIntegrator> &integrator)
{
  const auto internalEnergy = std::dynamic_pointer_cast<const pgo::Simulation::ImplicitBackwardEulerEnergy>(integrator->getInternalEnergy());
  if (!internalEnergy)
    return {};

  return {
    internalEnergy->getMinFeasibleAlphaThisSolve(),
    internalEnergy->getMinLineSearchAlphaThisSolve(),
    internalEnergy->getMinEffectiveAlphaThisSolve(),
  };
}

void resetSolveMaxStepSummary(const std::shared_ptr<pgo::Simulation::ImplicitBackwardEulerTimeIntegrator> &integrator)
{
  const auto internalEnergy = std::dynamic_pointer_cast<const pgo::Simulation::ImplicitBackwardEulerEnergy>(integrator->getInternalEnergy());
  if (internalEnergy)
    internalEnergy->resetSolveMaxStepStats();
}

void logRunIPCSimMaxStepSummary(
  const std::shared_ptr<pgo::SolidDeformationModel::DeformationModelEnergy> &elasticEnergy,
  const std::shared_ptr<pgo::Contact::CIPC::EmbeddedSurfaceIPCPotentialEnergy> &collisionHandler,
  const std::shared_ptr<pgo::Simulation::ImplicitBackwardEulerTimeIntegrator> &integrator)
{
  auto logger = pgo::Logging::lgr();
  if (!logger)
    return;

  const SolveMaxStepSummary summary = currentSolveMaxStepSummary(integrator);
  const auto materialClampCount = elasticEnergy->getMaterialClampCount();
  const auto contactClampCount = collisionHandler->getContactClampCount();

  if (logger->should_log(spdlog::level::info)) {
    SPDLOG_LOGGER_INFO(logger,
      "runIPCSim max-step summary: materialClampCount={} contactClampCount={} minMaterialFeasibleAlphaThisSolve={} minContactFeasibleAlphaThisSolve={} minFeasibleAlphaThisSolve={} minLineSearchAlphaThisSolve={} minEffectiveAlphaThisSolve={}",
      materialClampCount, contactClampCount,
      elasticEnergy->getMinMaterialFeasibleAlphaThisSolve(),
      collisionHandler->getMinContactFeasibleAlphaThisSolve(),
      summary.minFeasibleAlphaThisSolve,
      summary.minLineSearchAlphaThisSolve,
      summary.minEffectiveAlphaThisSolve);
  }
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

  const std::string configFilename = program.get<std::string>("config");
  const bool enableCliLog = program.get<bool>("--log");
  std::unique_ptr<RunSim::ScopedRunSimCliLogRedirect> logRedirect;

  try {
    ConfigFileJSON jconfig;
    if (jconfig.open(configFilename.c_str()) != true)
      return 0;

    if (enableCliLog) {
      const std::filesystem::path logPath = RunSim::deriveDefaultLogPathFromConfig(configFilename);
      logRedirect = std::make_unique<RunSim::ScopedRunSimCliLogRedirect>(logPath.string());
      pgo::Logging::init(nullptr, RunSim::resolveConfiguredLogLevel(jconfig));
    }
    else {
      pgo::Logging::init(nullptr, RunSim::resolveConfiguredLogLevel(jconfig));
    }
    pgo::Mesh::initPredicates();

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
    const std::string outputFolder = jconfig.getResolvedPath("output", 1);

    const bool hasTetMesh = jconfig.exist("tet-mesh");
    const bool hasCubicMesh = jconfig.exist("cubic-mesh");
    const bool useVolumePath = hasTetMesh || hasCubicMesh;

    RunIPCSim::IpcSimulationContext context = useVolumePath
      ? RunIPCSim::buildVolumeIpcSimulation(jconfig)
      : RunIPCSim::buildShellIpcSimulation(jconfig);

    const int n3 = static_cast<int>(context.simulationRestPosition.size());
    const int surfn3 = static_cast<int>(context.surfaceRestPositions.size());
    const int n = n3 / 3;

    ES::VXd g(n3);
    for (int vi = 0; vi < n; ++vi)
      g.segment<3>(vi * 3) = extAcc;

    ES::VXd fext(n3);
    ES::mv(context.M, g, fext);

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

    if (!std::filesystem::exists(outputFolder))
      std::filesystem::create_directories(outputFolder);

    int frameStart = -1;
    for (int framei = numSimSteps - 1; framei >= 0; --framei) {
      const std::string deformFilename = fmt::format("{}/deform{:04d}.u", outputFolder, framei);
      if (!std::filesystem::exists(deformFilename))
        continue;

      ES::MXd uMat(n3, 3);
      if (ES::readMatrix(deformFilename.c_str(), uMat) == 0) {
        frameStart = framei;
        u.noalias() = uMat.col(0);
        uvel.noalias() = uMat.col(1);
        uacc.noalias() = uMat.col(2);
        std::cout << "Restarting from frame " << framei << std::endl;
        break;
      }
    }

    if (frameStart < 0)
      std::cout << "No restart state found in " << outputFolder << ". Starting from frame 0." << std::endl;

    const double ratioDenom = numSimSteps > 1 ? static_cast<double>(numSimSteps - 1) : 1.0;

    bool executedStep = false;
    resetSolveMaxStepSummary(intg);
    context.elasticEnergy->resetMaterialMaxStepStats();
    context.collisionHandler->resetContactMaxStepStats();
    for (int framei = frameStart + 1; framei < numSimSteps; ++framei) {
      intg->clearGeneralImplicitForceModel();

      const double ratio = numSimSteps > 1 ? static_cast<double>(framei) / ratioDenom : 1.0;
      for (std::size_t pi = 0; pi < context.pullingEnergies.size(); ++pi) {
        const ES::VXd curTgt = context.pullingTargetRests[pi] * (1.0 - ratio) + context.pullingTargets[pi] * ratio;
        context.pullingEnergies[pi]->setTargetPos(curTgt.data());
        std::cout << "Frame " << framei << ", attachment " << pi << " target: " << curTgt.transpose().head(3) << std::endl;
      }

      context.elasticEnergy->resetMaterialMaxStepStats();
      context.collisionHandler->resetContactMaxStepStats();
      intg->addGeneralImplicitForceModel(context.collisionHandler, 0, 0);
      intg->setqState(u, uvel, uacc);
      intg->doTimestep(1, 3, 1);
      executedStep = true;
      intg->getq(u);
      intg->getqvel(uvel);
      intg->getqacc(uacc);
      logRunIPCSimMaxStepSummary(context.elasticEnergy, context.collisionHandler, intg);

      ES::MXd uMat(n3, 3);
      uMat.col(0) = u;
      uMat.col(1) = uvel;
      uMat.col(2) = uacc;
      ES::writeMatrix(fmt::format("{}/deform{:04d}.u", outputFolder, framei).c_str(), uMat);

      if (framei % frameGap == 0) {
        Mesh::TriMeshGeo mesh = context.surfaceMesh;
        ES::mv(context.surfaceFromSimulationDispMap, u, usurf);
        const ES::VXd psurf = context.surfaceRestPositions + usurf;
        for (int vi = 0; vi < mesh.numVertices(); ++vi)
          mesh.pos(vi) = psurf.segment<3>(vi * 3) / scale;
        mesh.save(fmt::format("{}/ret{:04d}.obj", outputFolder, framei / frameGap));
      }
    }

    if (!executedStep)
      logRunIPCSimMaxStepSummary(context.elasticEnergy, context.collisionHandler, intg);
  }
  catch (const std::exception &err) {
    SPDLOG_LOGGER_ERROR(Logging::lgr(), "{}", err.what());
    return 1;
  }

  return 0;
}
