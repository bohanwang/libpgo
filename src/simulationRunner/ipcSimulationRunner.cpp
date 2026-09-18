#include "simulationRunner.h"

#include "configFileJSON.h"
#include "EigenSupport.h"
#include "deformationModelEnergy.h"
#include "embeddedSurfaceIPCPotentialEnergy.h"
#include "implicitBackwardEulerTimeIntegrator.h"
#include "initPredicates.h"
#include "linearPotentialEnergy.h"
#include "multiVertexPullingSoftConstraints.h"
#include "NewtonSolver.h"
#include "pgoLogging.h"
#include "potentialEnergies.h"
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

void saveSurface(const pgo::RunIPCSim::IpcSimulationContext &context,
  ES::ConstRefVecXd displacement, double scale, const std::filesystem::path &filename)
{
  ES::VXd surfaceDisplacement(context.surfaceRestPositions.size());
  ES::mv(context.surfaceFromSimulationDispMap, displacement, surfaceDisplacement);
  const ES::VXd surfacePositions = context.surfaceRestPositions + surfaceDisplacement;

  pgo::Mesh::TriMeshGeo mesh = context.surfaceMesh;
  for (int vi = 0; vi < mesh.numVertices(); ++vi)
    mesh.pos(vi) = surfacePositions.segment<3>(vi * 3) / scale;

  if (!mesh.save(filename.string()))
    throw std::runtime_error("Failed to save IPC surface output: " + filename.string());
}

}  // namespace

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
    const double scale = jconfig.getDouble("scale", 1);
    const double solverEps = jconfig.getDouble("solver-eps", 1);
    const int solverMaxIter = jconfig.getInt("solver-max-iter", 1);
    const std::string simType = jconfig.getString("sim-type");
    if (simType != "dynamic" && simType != "static")
      throw std::invalid_argument("runIPCSim supports `sim-type` values `dynamic` and `static`.");
    const std::filesystem::path outputFolder = jconfig.getResolvedPath("output", 1);
    const bool restartFromU = jconfig.exist("restart-from-u") ? jconfig.getValue<bool>("restart-from-u", 1) : false;
    if (simType == "static" && restartFromU)
      throw std::invalid_argument("`restart-from-u=true` is only supported for dynamic IPC simulations.");

    std::filesystem::create_directories(outputFolder);

    std::unique_ptr<RunSim::ScopedRunSimCliLogRedirect> logRedirect;
    if (enableCliLog) {
      const std::filesystem::path logPath = outputFolder / "runIPCSim.log";
      logRedirect = std::make_unique<RunSim::ScopedRunSimCliLogRedirect>(logPath.string());
    }
    pgo::Logging::init(nullptr, RunSim::resolveConfiguredLogLevel(jconfig));

    if (!restartFromU)
      std::cout << "restart-from-u=false; existing output files will be overwritten as frames are written." << std::endl;
    pgo::Mesh::initPredicates();

    const bool hasTetMesh = jconfig.exist("tet-mesh");
    const bool hasCubicMesh = jconfig.exist("cubic-mesh");
    const bool useVolumePath = hasTetMesh || hasCubicMesh;

    RunIPCSim::IpcSimulationContext context = useVolumePath ? RunIPCSim::buildVolumeIpcSimulation(jconfig) : RunIPCSim::buildShellIpcSimulation(jconfig);

    const int n3 = static_cast<int>(context.simulationRestPosition.size());
    const int n = n3 / 3;
    if (context.initialDisplacement.size() != n3 || !context.initialDisplacement.allFinite())
      throw std::runtime_error("IPC setup produced an invalid initial displacement.");

    ES::VXd g(n3);
    for (int vi = 0; vi < n; ++vi)
      g.segment<3>(vi * 3) = extAcc;

    ES::VXd fext(n3);
    ES::mv(context.M, g, fext);

    if (simType == "static") {
      ES::VXd u = context.initialDisplacement;
      const std::filesystem::path statePath = outputFolder / "deform0000.u";
      context.collisionHandler->validateCollisionFreeState(u);

      auto externalForcesEnergy =
        std::make_shared<PredefinedPotentialEnergies::LinearPotentialEnergy>(fext);
      auto energyAll =
        std::make_shared<NonlinearOptimization::PotentialEnergies>(n3);
      energyAll->addPotentialEnergy(context.elasticEnergy);
      for (const auto &pullingEnergy : context.pullingEnergies)
        energyAll->addPotentialEnergy(pullingEnergy);
      energyAll->addPotentialEnergy(externalForcesEnergy, -1.0);
      energyAll->addPotentialEnergy(context.collisionHandler);
      energyAll->init();

      NonlinearOptimization::NewtonSolver::SolverParam solverParam;
      NonlinearOptimization::NewtonSolver solver(
        u.data(), solverParam, energyAll, std::vector<int>(), nullptr);
      if (solver.solve(u.data(), solverMaxIter, solverEps, 2) != 0)
        return 1;
      if (!u.allFinite())
        throw std::runtime_error("Static IPC solve produced non-finite displacement.");

      ES::MXd state = ES::MXd::Zero(n3, 3);
      state.col(0) = u;
      if (ES::writeMatrix(statePath.string().c_str(), state) != 0)
        throw std::runtime_error("Failed to save static IPC state: " + statePath.string());
      saveSurface(context, u, scale, outputFolder / "ret0000.obj");
      return 0;
    }

    const ES::V3d initialVel = ES::Mp<ES::V3d>(jconfig.getValue<std::array<double, 3>>("init-vel", 1).data());
    if (!initialVel.allFinite())
      throw std::invalid_argument("`init-vel` must contain exactly three finite values.");
    const double timestep = jconfig.getDouble("timestep", 1);
    const std::array<double, 2> dampingParams = jconfig.getValue<std::array<double, 2>>("damping-params", 1);
    const int numSimSteps = jconfig.getInt("num-timestep", 1);
    const int frameGap = jconfig.getInt("dump-interval", 1);

    std::shared_ptr<Simulation::ImplicitBackwardEulerTimeIntegrator> intg =
      std::make_shared<Simulation::ImplicitBackwardEulerTimeIntegrator>(context.M, context.elasticEnergy,
        dampingParams[0], dampingParams[1], timestep, solverMaxIter, solverEps);

    for (auto &pullingEnergy : context.pullingEnergies)
      intg->addImplicitForceModel(pullingEnergy, 0, 0);

    intg->setExternalForce(fext.data());

    ES::VXd u = context.initialDisplacement;
    ES::VXd uvel = ES::VXd::Zero(n3);
    ES::VXd uacc = ES::VXd::Zero(n3);

    for (int i = 0; i < n; ++i)
      uvel.segment<3>(i * 3) = initialVel;

    int frameStart = -1;
    if (restartFromU) {
      for (int framei = numSimSteps - 1; framei >= 0; --framei) {
        const std::filesystem::path deformPath = outputFolder / fmt::format("deform{:04d}.u", framei);
        const std::string deformFilename = deformPath.string();
        if (!std::filesystem::exists(deformPath))
          continue;

        ES::MXd uMat;
        if (ES::readMatrix(deformFilename.c_str(), uMat) != 0)
          throw std::runtime_error("Failed to read dynamic IPC restart state: " + deformFilename);
        if (uMat.rows() != n3 || uMat.cols() < 3 || !uMat.leftCols(3).allFinite())
          throw std::runtime_error("Invalid dynamic IPC restart state: " + deformFilename);
        frameStart = framei;
        u.noalias() = uMat.col(0);
        uvel.noalias() = uMat.col(1);
        uacc.noalias() = uMat.col(2);
        std::cout << "Restarting from frame " << framei << std::endl;
        break;
      }
    }

    if (frameStart < 0 && restartFromU)
      std::cout << "No restart state found in " << outputFolder << ". Starting from frame 0." << std::endl;
    else if (frameStart < 0)
      std::cout << "Starting from frame 0." << std::endl;

    context.collisionHandler->validateCollisionFreeState(u);

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
      intg->setqState(u, uvel, uacc);
      intg->doTimestep(1, 3, 1);
      intg->getq(u);
      intg->getqvel(uvel);
      intg->getqacc(uacc);

      ES::MXd uMat(n3, 3);
      uMat.col(0) = u;
      uMat.col(1) = uvel;
      uMat.col(2) = uacc;
      const std::filesystem::path statePath = outputFolder / fmt::format("deform{:04d}.u", framei);
      if (ES::writeMatrix(statePath.string().c_str(), uMat) != 0)
        throw std::runtime_error("Failed to save dynamic IPC state: " + statePath.string());
      if (framei % frameGap == 0)
        saveSurface(context, u, scale,
          outputFolder / fmt::format("ret{:04d}.obj", framei / frameGap));
    }
  }
  catch (const std::exception &err) {
    SPDLOG_LOGGER_ERROR(Logging::lgr(), "{}", err.what());
    return 1;
  }

  return 0;
}
