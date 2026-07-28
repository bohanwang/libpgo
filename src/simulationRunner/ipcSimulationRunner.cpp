#include "simulationRunner.h"

#include "configFileJSON.h"
#include "EigenSupport.h"
#include "deformationModelEnergy.h"
#include "embeddedSurfaceIPCPotentialEnergy.h"
#include "implicitBackwardEulerTimeIntegrator.h"
#include "initPredicates.h"
#include "multiVertexPullingSoftConstraints.h"
#include "pgoLogging.h"
#include "runIPCSimSetup.h"
#include "runSimCliLogging.h"
#include <fmt/format.h>

#include <array>
#include <filesystem>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

namespace
{
namespace ES = pgo::EigenSupport;

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

}

int pgo::SimulationRunner::runIPCSimulationFromConfig(
  const std::filesystem::path &configFilename, bool enableCliLog)
{
  using namespace pgo;
  namespace ES = EigenSupport;

  try {
    ConfigFileJSON jconfig;
    if (jconfig.open(configFilename.string().c_str()) != true)
      return 1;

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
    const bool restartFromU = jconfig.exist("restart-from-u") ? jconfig.getValue<bool>("restart-from-u", 1) : false;

    if (restartFromU) {
      std::filesystem::create_directories(outputFolder);
    }
    else {
      clearOutputDirectory(outputFolder);
    }

    std::unique_ptr<RunSim::ScopedRunSimCliLogRedirect> logRedirect;
    if (enableCliLog) {
      const std::filesystem::path logPath = outputFolder / "runIPCSim.log";
      logRedirect = std::make_unique<RunSim::ScopedRunSimCliLogRedirect>(logPath.string());
    }
    pgo::Logging::init(nullptr, RunSim::resolveConfiguredLogLevel(jconfig));

    if (!restartFromU) {
      std::cout << "restart-from-u=false; clearing output folder " << outputFolder << "." << std::endl;
    }
    pgo::Mesh::initPredicates();

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

    int frameStart = -1;
    if (restartFromU) {
      for (int framei = numSimSteps - 1; framei >= 0; --framei) {
        const std::string deformFilename = fmt::format("{}/deform{:04d}.u", outputFolder.string(), framei);
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
    }

    if (frameStart < 0 && restartFromU)
      std::cout << "No restart state found in " << outputFolder << ". Starting from frame 0." << std::endl;
    else if (frameStart < 0)
      std::cout << "Starting from frame 0." << std::endl;

    const double ratioDenom = numSimSteps > 1 ? static_cast<double>(numSimSteps - 1) : 1.0;

    for (int framei = frameStart + 1; framei < numSimSteps; ++framei) {
      intg->clearGeneralImplicitForceModel();

      const double ratio = numSimSteps > 1 ? static_cast<double>(framei) / ratioDenom : 1.0;
      for (std::size_t pi = 0; pi < context.pullingEnergies.size(); ++pi) {
        const ES::VXd curTgt = context.pullingTargetRests[pi] * (1.0 - ratio) + context.pullingTargets[pi] * ratio;
        context.pullingEnergies[pi]->setTargetPos(curTgt.data());
        std::cout << "Frame " << framei << ", attachment " << pi << " target: " << curTgt.transpose().head(3) << std::endl;
      }

      intg->addGeneralImplicitForceModel(context.collisionHandler, 0, 0);
      for (const auto &forceModel : context.extraGeneralImplicitForceModels)
        intg->addGeneralImplicitForceModel(forceModel, 0, 0);
      intg->setqState(u, uvel, uacc);
      intg->doTimestep(1, 3, 1);
      intg->getq(u);
      intg->getqvel(uvel);
      intg->getqacc(uacc);

      ES::MXd uMat(n3, 3);
      uMat.col(0) = u;
      uMat.col(1) = uvel;
      uMat.col(2) = uacc;
      ES::writeMatrix(fmt::format("{}/deform{:04d}.u", outputFolder.string(), framei).c_str(), uMat);

      if (framei % frameGap == 0) {
        Mesh::TriMeshGeo mesh = context.surfaceMesh;
        ES::mv(context.surfaceFromSimulationDispMap, u, usurf);
        const ES::VXd psurf = context.surfaceRestPositions + usurf;
        for (int vi = 0; vi < mesh.numVertices(); ++vi)
          mesh.pos(vi) = psurf.segment<3>(vi * 3) / scale;
        mesh.save(fmt::format("{}/ret{:04d}.obj", outputFolder.string(), framei / frameGap));
      }
    }

  }
  catch (const std::exception &err) {
    SPDLOG_LOGGER_ERROR(Logging::lgr(), "{}", err.what());
    return 1;
  }

  return 0;
}
