#include "configFileJSON.h"
#include "EigenSupport.h"
#include "initPredicates.h"
#include "triMeshGeo.h"
#include "pgoLogging.h"
#include "simulationMesh.h"
#include "deformationModelManager.h"
#include "basicIO.h"
#include "deformationModelAssembler.h"
#include "deformationModelEnergy.h"
#include "multiVertexPullingSoftConstraints.h"
#include "implicitBackwardEulerTimeIntegrator.h"
#include "generateMassMatrix.h"
#include "libiglInterface.h"
#include "embeddedSurfaceIPCPotentialEnergy.h"
#include "runSimCliLogging.h"

#include <argparse/argparse.hpp>

#include <filesystem>
#include <iostream>
#include <memory>
#include <numeric>
#include <stdexcept>

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
    throwConfigError(std::string("runIPCSim phase1BC does not support `") + field + "`: " + reason);
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

void validateZeroInitialDisplacement(const ES::V3d &initialDisp)
{
  if (initialDisp.cwiseAbs().maxCoeff() > 0.0)
    throwConfigError("runIPCSim phase1BC only accepts zero `init-disp`.");
}
}  // namespace

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

  std::string configFilename = program.get<std::string>("config");
  const bool enableCliLog = program.get<bool>("--log");
  std::unique_ptr<RunSim::ScopedRunSimCliLogRedirect> logRedirect;

  if (enableCliLog) {
    try {
      const std::filesystem::path logPath = RunSim::deriveDefaultLogPathFromConfig(configFilename);
      logRedirect = std::make_unique<RunSim::ScopedRunSimCliLogRedirect>(logPath.string());
      pgo::Logging::init();
    }
    catch (const std::exception &err) {
      std::cerr << err.what() << std::endl;
      return 1;
    }
  }
  else {
    pgo::Logging::init();
  }
  pgo::Mesh::initPredicates();

  try {
    ConfigFileJSON jconfig;
    if (jconfig.open(configFilename.c_str()) != true)
      return 0;

    if (jconfig.exist("tet-mesh") || jconfig.exist("cubic-mesh"))
      throwConfigError("runIPCSim phase1BC is shell-only. tet/cubic inputs are not supported yet.");

    rejectIfPresent(jconfig, "external-objects", "external contact is out of scope for phase1BC.");
    if (!jconfig.exist("fixed-vertices"))
      throwConfigError("Missing required field `fixed-vertices`.");
    if (!jconfig.exist("ipc-dhat"))
      throwConfigError("Missing required field `ipc-dhat`.");
    if (!jconfig.exist("ipc-kappa"))
      throwConfigError("Missing required field `ipc-kappa`.");

    const std::string surfaceMeshFilename = jconfig.getResolvedPath("surface-mesh", 1);
    const ES::V3d extAcc = ES::Mp<ES::V3d>(jconfig.getValue<std::array<double, 3>>("g", 1).data());
    const ES::V3d initialVel = ES::Mp<ES::V3d>(jconfig.getValue<std::array<double, 3>>("init-vel", 1).data());
    const ES::V3d initialDisp = ES::Mp<ES::V3d>(jconfig.getValue<std::array<double, 3>>("init-disp", 1).data());
    validateZeroInitialDisplacement(initialDisp);

    const double timestep = jconfig.getDouble("timestep", 1);
    const double scale = jconfig.getDouble("scale", 1);
    PGO_ALOG(std::abs(scale - 1.0) < 1e-6);

    const double solverEps = jconfig.getDouble("solver-eps", 1);
    const int solverMaxIter = jconfig.getInt("solver-max-iter", 1);
    const std::array<double, 2> dampingParams = jconfig.getValue<std::array<double, 2>>("damping-params", 1);
    const std::string material = jconfig.getString("elastic-material");
    if (material != "koiter-stvk")
      throwConfigError("runIPCSim phase1BC only supports `elastic-material = koiter-stvk`.");

    const int numSimSteps = jconfig.getInt("num-timestep", 1);
    const int frameGap = jconfig.getInt("dump-interval", 1);
    const std::string simType = jconfig.getString("sim-type");
    if (simType != "dynamic")
      throwConfigError("runIPCSim phase1BC only supports `sim-type = dynamic`.");

    const std::string outputFolder = jconfig.getResolvedPath("output", 1);

    Contact::CIPC::SurfaceIPCCore::Parameters ipcParams;
    ipcParams.dhat = jconfig.getDouble("ipc-dhat", 1);
    ipcParams.kappa = jconfig.getDouble("ipc-kappa", 1);
    ipcParams.eps_ee = 0.0;
    ipcParams.slackness = 1.0;

    std::cout << "runIPCSim phase1BC shell IPC parameters: "
              << "ipc-dhat=" << ipcParams.dhat << ", "
              << "ipc-kappa=" << ipcParams.kappa << ", "
              << "eps_ee=" << ipcParams.eps_ee << ", "
              << "slackness=" << ipcParams.slackness << std::endl;

    Mesh::TriMeshGeo surfaceMesh;
    if (surfaceMesh.load(surfaceMeshFilename) != true)
      return 1;

    int surfn = surfaceMesh.numVertices();
    int surfn3 = surfn * 3;
    ES::VXd surfaceRestPositions(surfn3);
    for (int vi = 0; vi < surfn; vi++)
      surfaceRestPositions.segment<3>(vi * 3) = surfaceMesh.pos(vi);

    SolidDeformationModel::SimulationMeshENuhMaterial matParam(10000, 0.4, 0.001);
    std::shared_ptr<SolidDeformationModel::SimulationMesh> simMesh(
      SolidDeformationModel::loadShellMesh(surfaceMesh, &matParam));
    std::shared_ptr<SolidDeformationModel::DeformationModelManager> dmm = std::make_shared<SolidDeformationModel::DeformationModelManager>();

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

    double E = 1000000;
    double h = 3e-3;
    ES::VXd elasticParams(5 * nele);
    for (int ei = 0; ei < nele; ei++) {
      double E_bend = E;
      double nu = 0.4;

      elasticParams[ei * 5 + 0] = E;
      elasticParams[ei * 5 + 1] = nu;
      elasticParams[ei * 5 + 2] = E_bend;
      elasticParams[ei * 5 + 3] = nu;
      elasticParams[ei * 5 + 4] = h;
    }

    ES::VXd restPosition(n3);
    for (int vi = 0; vi < n; vi++) {
      double p[3];
      simMesh->getVertex(vi, p);
      restPosition.segment<3>(vi * 3) = ES::V3d(p[0], p[1], p[2]);
    }

    std::shared_ptr<SolidDeformationModel::DeformationModelEnergy> elasticEnergy =
      std::make_shared<SolidDeformationModel::DeformationModelEnergy>(assembler, &restPosition, 0);
    elasticEnergy->setElasticParams(elasticParams);

    ES::VXd zero(n3);
    zero.setZero();
    ES::SpMatD K;
    elasticEnergy->createHessian(K);
    elasticEnergy->hessian(zero, K);

    std::vector<std::shared_ptr<ConstraintPotentialEnergies::MultipleVertexPulling>> pullingEnergies;
    std::vector<ES::VXd> pullingTargets, pullingTargetRests;
    for (const auto &fv : jconfig.handle()["fixed-vertices"]) {
      std::string filename = jconfig.resolvePath(fv["filename"].get<std::string>());
      std::array<double, 3> movement = fv["movement"].get<std::array<double, 3>>();
      double attachmentCoeff = fv["coeff"].get<double>();

      std::vector<int> fixedVertices;
      if (BasicIO::read1DText(filename.c_str(), std::back_inserter(fixedVertices)) != 0)
        return 1;
      std::sort(fixedVertices.begin(), fixedVertices.end());

      ES::VXd tgtVertexPositions(fixedVertices.size() * 3);
      ES::VXd tgtVertexRests(fixedVertices.size() * 3);
      for (int vi = 0; vi < static_cast<int>(fixedVertices.size()); vi++) {
        tgtVertexPositions.segment<3>(vi * 3) = restPosition.segment<3>(fixedVertices[vi] * 3) + ES::Mp<ES::V3d>(movement.data());
        tgtVertexRests.segment<3>(vi * 3) = restPosition.segment<3>(fixedVertices[vi] * 3);
      }

      auto pullingEnergy = std::make_shared<ConstraintPotentialEnergies::MultipleVertexPulling>(K, restPosition.data(),
        static_cast<int>(fixedVertices.size()), fixedVertices.data(), tgtVertexPositions.data(), nullptr, 1);
      pullingEnergy->setCoeff(attachmentCoeff);
      pullingEnergies.push_back(pullingEnergy);
      pullingTargets.push_back(tgtVertexPositions);
      pullingTargetRests.push_back(tgtVertexRests);
    }

    ES::SpMatD M;
    libiglInterface::computeMassMatrix(surfaceMesh, M, 1, 1);
    M *= 100;

    ES::VXd g(n3);
    for (int vi = 0; vi < n; vi++)
      g.segment<3>(vi * 3) = extAcc;

    ES::VXd fext(n3);
    ES::mv(M, g, fext);

    ES::MXd V;
    ES::MXi F;
    Mesh::triMeshGeoToMatrices(surfaceMesh, V, F);
    std::shared_ptr<Contact::CIPC::EmbeddedSurfaceIPCPotentialEnergy> collisionHandler =
      std::make_shared<Contact::CIPC::EmbeddedSurfaceIPCPotentialEnergy>(V, F, makeIdentityEmbedding(n3), ipcParams);

    std::shared_ptr<Simulation::ImplicitBackwardEulerTimeIntegrator> intg =
      std::make_shared<Simulation::ImplicitBackwardEulerTimeIntegrator>(M, elasticEnergy,
        dampingParams[0], dampingParams[1], timestep, solverMaxIter, solverEps);

    for (auto &pullingEnergy : pullingEnergies)
      intg->addImplicitForceModel(pullingEnergy, 0, 0);

    intg->setExternalForce(fext.data());

    ES::VXd u(n3), uvel(n3), uacc(n3);
    u.setZero();
    uvel.setZero();
    uacc.setZero();

    for (int i = 0; i < n; i++)
      uvel.segment<3>(i * 3) = initialVel;

    if (!std::filesystem::exists(outputFolder))
      std::filesystem::create_directories(outputFolder);

    int frameStart = -1;
    for (int framei = numSimSteps - 1; framei >= 0; framei--) {
      if (!std::filesystem::exists(fmt::format("{}/deform{:04d}.u", outputFolder, framei)))
        continue;

      ES::MXd uMat(n3, 3);
      if (ES::readMatrix(fmt::format("{}/deform{:04d}.u", outputFolder, framei).c_str(), uMat) == 0) {
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

    for (int framei = frameStart + 1; framei < numSimSteps; framei++) {
      intg->clearGeneralImplicitForceModel();

      const double ratio = numSimSteps > 1 ? static_cast<double>(framei) / ratioDenom : 1.0;
      for (std::size_t pi = 0; pi < pullingEnergies.size(); pi++) {
        ES::VXd restTgt = pullingTargetRests[pi];
        ES::VXd curTgt = restTgt * (1 - ratio) + pullingTargets[pi] * ratio;
        pullingEnergies[pi]->setTargetPos(curTgt.data());
        std::cout << "Frame " << framei << ", attachment " << pi << " target: " << curTgt.transpose().head(3) << std::endl;
      }

      intg->addGeneralImplicitForceModel(collisionHandler, 0, 0);
      intg->setqState(u, uvel, uacc);
      intg->doTimestep(1, 3, 1);
      intg->getq(u);
      intg->getqvel(uvel);
      intg->getqacc(uacc);

      if (framei % frameGap == 0) {
        Mesh::TriMeshGeo mesh = surfaceMesh;
        const ES::VXd psurf = surfaceRestPositions + u;
        for (int vi = 0; vi < mesh.numVertices(); vi++)
          mesh.pos(vi) = psurf.segment<3>(vi * 3) / scale;
        mesh.save(fmt::format("{}/ret{:04d}.obj", outputFolder, framei / frameGap));

        ES::MXd uMat(n3, 3);
        uMat.col(0) = u;
        uMat.col(1) = uvel;
        uMat.col(2) = uacc;
        ES::writeMatrix(fmt::format("{}/deform{:04d}.u", outputFolder, framei).c_str(), uMat);
      }
    }
  }
  catch (const std::exception &err) {
    SPDLOG_LOGGER_ERROR(Logging::lgr(), "{}", err.what());
    return 1;
  }

  return 0;
}
