#include "configFileJSON.h"
#include "EigenSupport.h"
#include "initPredicates.h"
#include "triMeshGeo.h"
#include "pgoLogging.h"
#include "geometryQuery.h"
#include "simulationMesh.h"
#include "deformationModelManager.h"
#include "basicIO.h"
#include "deformationModelAssembler.h"
#include "deformationModelEnergy.h"
#include "plasticModel.h"
#include "multiVertexPullingSoftConstraints.h"
#include "implicitBackwardEulerTimeIntegrator.h"
#include "TRBDF2TimeIntegrator.h"
#include "generateMassMatrix.h"
#include "generateSurfaceMesh.h"
#include "linearPotentialEnergy.h"
#include "NewtonRaphsonSolver.h"
#include "createTriMesh.h"
#include "libiglInterface.h"

#include <argparse/argparse.hpp>

#include <iostream>

int main(int argc, char *argv[])
{
  using namespace pgo;
  namespace ES = EigenSupport;

  // git add subparser
  argparse::ArgumentParser program("Run Simulation");
  program.add_argument("config")
    .help("Config File")
    .required();

  try {
    program.parse_args(argc, argv);  // Example: ./main --color orange
  }
  catch (const std::exception &err) {
    std::cerr << err.what() << std::endl;
    std::cerr << program;
    return 1;
  }

  pgo::Logging::init();
  pgo::Mesh::initPredicates();

  std::string configFilename = program.get<std::string>("config");

  ConfigFileJSON jconfig;
  if (jconfig.open(configFilename.c_str()) != true) {
    return 0;
  }

  // surface mesh filename
  std::string surfaceMeshFilename = jconfig.getString("surface-mesh", 1);

  // external acceleration
  ES::V3d extAcc = ES::Mp<ES::V3d>(jconfig.getValue<std::array<double, 3>>("g", 1).data());

  // initial velocity
  ES::V3d initialVel = ES::Mp<ES::V3d>(jconfig.getValue<std::array<double, 3>>("init-vel", 1).data());

  // timestep
  double timestep = jconfig.getDouble("timestep", 1);

  // solver param
  double solverEps = jconfig.getDouble("solver-eps", 1);
  int solverMaxIter = jconfig.getInt("solver-max-iter", 1);

  // damping
  std::array<double, 2> dampingParams = jconfig.getValue<std::array<double, 2>>("damping-params", 1);

  // material
  std::string material = jconfig.getString("elastic-material");
  SolidDeformationModel::DeformationModelElasticMaterial elasticMat;
  SolidDeformationModel::SimulationMeshENuhMaterial matParam(1000000, 0.4, 0.01);

  if (material == "kl-stvk") {
    elasticMat = pgo::SolidDeformationModel::DeformationModelElasticMaterial::STVK;
  }
  else if (material == "kl-linear") {
    elasticMat = pgo::SolidDeformationModel::DeformationModelElasticMaterial::LINEAR;
  }
  else {
    SPDLOG_LOGGER_ERROR(Logging::lgr(), "Unsupported elastic material: {}", material);
    return 1;
  }

  int numSimSteps = jconfig.getInt("num-timestep", 1);
  int frameGap = jconfig.getInt("dump-interval", 1);

  // sim type
  std::string simType = jconfig.getString("sim-type");

  // output
  std::string outputFolder = jconfig.getString("output", 1);

  Mesh::TriMeshGeo surfaceMesh;
  if (surfaceMesh.load(surfaceMeshFilename) != true)
    return 1;

  int surfn = surfaceMesh.numVertices();
  int surfn3 = surfn * 3;
  ES::VXd surfaceRestPositions(surfn3);
  for (int vi = 0; vi < surfn; vi++) {
    surfaceRestPositions.segment<3>(vi * 3) = surfaceMesh.pos(vi);
  }

  // initialize fem
  std::shared_ptr<SolidDeformationModel::SimulationMesh> simMesh(SolidDeformationModel::loadShellMesh(surfaceMesh, &matParam));
  std::shared_ptr<SolidDeformationModel::DeformationModelManager> dmm = std::make_shared<SolidDeformationModel::DeformationModelManager>();

  dmm->setMesh(simMesh.get(), nullptr, nullptr);
  dmm->init(pgo::SolidDeformationModel::DeformationModelPlasticMaterial::SHELL_FF_DOF0, elasticMat, 1);

  std::vector<double> elementWeights(simMesh->getNumElements(), 1.0);
  std::shared_ptr<SolidDeformationModel::DeformationModelAssembler> assembler =
    std::make_shared<SolidDeformationModel::DeformationModelAssembler>(dmm, elementWeights.data());

  int n = simMesh->getNumVertices();
  int n3 = n * 3;
  int nele = simMesh->getNumElements();

  ES::VXd restPosition(n3);
  for (int vi = 0; vi < n; vi++) {
    double p[3];
    simMesh->getVertex(vi, p);
    restPosition.segment<3>(vi * 3) = ES::V3d(p[0], p[1], p[2]);
  }

  std::shared_ptr<SolidDeformationModel::DeformationModelEnergy> elasticEnergy = std::make_shared<SolidDeformationModel::DeformationModelEnergy>(assembler, &restPosition, 0);

  ES::VXd zero(n3);
  zero.setZero();

  ES::SpMatD K;
  elasticEnergy->createHessian(K);
  elasticEnergy->hessian(zero, K);

  // attachments
  std::vector<std::shared_ptr<ConstraintPotentialEnergies::MultipleVertexPulling>> pullingEnergies;
  std::vector<ES::VXd> pullingTargets, pullingTargetRests;
  for (const auto &fv : jconfig.handle()["fixed-vertices"]) {
    std::string filename = fv["filename"].get<std::string>();
    std::array<double, 3> movement = fv["movement"].get<std::array<double, 3>>();
    double attachmentCoeff = fv["coeff"].get<double>();

    std::vector<int> fixedVertices;
    if (BasicIO::read1DText(filename.c_str(), std::back_inserter(fixedVertices)) != 0) {
      return 1;
    }
    std::sort(fixedVertices.begin(), fixedVertices.end());

    ES::VXd tgtVertexPositions(fixedVertices.size() * 3);
    ES::VXd tgtVertexRests(fixedVertices.size() * 3);
    for (int vi = 0; vi < (int)fixedVertices.size(); vi++) {
      tgtVertexPositions.segment<3>(vi * 3) = restPosition.segment<3>(fixedVertices[vi] * 3) + ES::Mp<ES::V3d>(movement.data());
      tgtVertexRests.segment<3>(vi * 3) = restPosition.segment<3>(fixedVertices[vi] * 3);
    }

    // initialize fixed constraints
    auto pullingEnergy = std::make_shared<ConstraintPotentialEnergies::MultipleVertexPulling>(K, restPosition.data(),
      (int)fixedVertices.size(), fixedVertices.data(), tgtVertexPositions.data(), nullptr, 1);
    pullingEnergy->setCoeff(attachmentCoeff);
    pullingEnergies.push_back(pullingEnergy);
    pullingTargets.push_back(tgtVertexPositions);
    pullingTargetRests.push_back(tgtVertexRests);
  }

  ES::SpMatD M;
  libiglInterface::computeMassMatrix(surfaceMesh, M, 1, 1);

  // initialize gravity
  ES::VXd g(n3);
  for (int vi = 0; vi < n; vi++) {
    g.segment<3>(vi * 3) = extAcc;
  }

  ES::VXd fext(n3);
  ES::mv(M, g, fext);

  if (simType == "dynamic") {
    // initialize contact
    std::shared_ptr<Simulation::ImplicitBackwardEulerTimeIntegrator> intg =
      std::make_shared<Simulation::ImplicitBackwardEulerTimeIntegrator>(M, elasticEnergy,
        dampingParams[0], dampingParams[1], timestep, solverMaxIter, solverEps);

#if defined(PGO_HAS_KNITRO)
    intg->setSolverOption(Simulation::TimeIntegratorSolverOption::SO_KNITRO);
    intg->setSolverConfigFile("config.opt");
#endif

    for (auto pullingEnergy : pullingEnergies)
      intg->addImplicitForceModel(pullingEnergy, 0, 0);

    intg->setExternalForce(fext.data());

    ES::VXd curSurfacePos = surfaceRestPositions;
    ES::VXd x = restPosition, u(n3);
    ES::VXd uvel(n3), uacc(n3), usurf(surfn3);

    usurf.setZero();
    u.setZero();
    uvel.setZero();
    uacc.setZero();

    for (int i = 0; i < n; i++) {
      uvel.segment<3>(i * 3) = ES::V3d(initialVel[0], initialVel[1], initialVel[2]);
    }

    if (!std::filesystem::exists(outputFolder)) {
      std::filesystem::create_directories(outputFolder);
    }

    for (int framei = 0; framei < numSimSteps; framei++) {
      intg->clearGeneralImplicitForceModel();

      intg->setqState(u, uvel, uacc);
      intg->doTimestep(1, 2, 1);

      intg->getq(u);
      intg->getq(uvel);
      intg->getq(uacc);

      if (framei % frameGap == 0) {
        ES::VXd psurf = surfaceRestPositions + u;

        Mesh::TriMeshGeo mesh = surfaceMesh;
        for (int vi = 0; vi < mesh.numVertices(); vi++) {
          mesh.pos(vi) = psurf.segment<3>(vi * 3);
        }
        mesh.save(fmt::format("{}/ret{:04d}.obj", outputFolder, framei / frameGap));
      }
    }
  }
  else if (simType == "static") {
    std::shared_ptr<PredefinedPotentialEnergies::LinearPotentialEnergy> externalForcesEnergy = std::make_shared<PredefinedPotentialEnergies::LinearPotentialEnergy>(fext);

    std::shared_ptr<NonlinearOptimization::PotentialEnergies> energyAll = std::make_shared<NonlinearOptimization::PotentialEnergies>(n3);
    energyAll->addPotentialEnergy(elasticEnergy);
    for (auto eng : pullingEnergies)
      energyAll->addPotentialEnergy(eng, 1.0);
    
    energyAll->addPotentialEnergy(externalForcesEnergy, -1.0);
    energyAll->init();

    NonlinearOptimization::NewtonRaphsonSolver::SolverParam solverParam;

    ES::VXd u(n3);
    u.setZero();

    energyAll->printEnergy(u);

    NonlinearOptimization::NewtonRaphsonSolver solver(u.data(), solverParam, energyAll, std::vector<int>(), nullptr);
    solver.solve(u.data(), solverMaxIter, solverEps, 2);

    ES::VXd x = restPosition + u;

    Mesh::TriMeshGeo meshOut = surfaceMesh;
    for (int vi = 0; vi < meshOut.numVertices(); vi++) {
      meshOut.pos(vi) = x.segment<3>(vi * 3);
    }
    meshOut.save(outputFolder);
  }

  return 0;
}