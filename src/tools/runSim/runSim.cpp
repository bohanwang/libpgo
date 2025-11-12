#include "configFileJSON.h"
#include "EigenSupport.h"
#include "initPredicates.h"
#include "triMeshGeo.h"
#include "generateTetMeshMatrix.h"
#include "tetMesh.h"
#include "pgoLogging.h"
#include "geometryQuery.h"
#include "simulationMesh.h"
#include "deformationModelManager.h"
#include "tetMeshDeformationModel.h"
#include "plasticModel3DDeformationGradient.h"
#include "basicIO.h"
#include "deformationModelAssembler.h"
#include "deformationModelEnergy.h"
#include "multiVertexPullingSoftConstraints.h"
#include "implicitBackwardEulerTimeIntegrator.h"
#include "TRBDF2TimeIntegrator.h"
#include "generateMassMatrix.h"
#include "generateSurfaceMesh.h"
#include "barycentricCoordinates.h"
#include "triangleMeshExternalContactHandler.h"
#include "pointPenetrationEnergy.h"
#include "triangleMeshSelfContactHandler.h"
#include "pointTrianglePairCouplingEnergyWithCollision.h"
#include "linearPotentialEnergy.h"
#include "NewtonRaphsonSolver.h"
#include "createTriMesh.h"
#include "finiteDifference.h"

#include <argparse/argparse.hpp>

#include <tbb/global_control.h>

#include <thread>
#include <iostream>

int main(int argc, char *argv[])
{
  using namespace pgo;
  using namespace pgo::EigenSupport;
  namespace ES = EigenSupport;

  tbb::global_control c(tbb::global_control::max_allowed_parallelism, std::min(64, (int)std::thread::hardware_concurrency()));
  tbb::global_control global_limit(tbb::global_control::thread_stack_size, 16 * 1024 * 1024);

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

  // tet mesh filename
  std::string tetMeshFilename = jconfig.getString("tet-mesh", 1);

  // surface mesh filename
  std::string surfaceMeshFilename = jconfig.getString("surface-mesh", 1);

  // external acceleration
  ES::V3d extAcc = ES::Mp<ES::V3d>(jconfig.getValue<std::array<double, 3>>("g", 1).data());

  // initial velocity
  ES::V3d initialVel = ES::Mp<ES::V3d>(jconfig.getValue<std::array<double, 3>>("init-vel", 1).data());
  ES::V3d initialDisp = ES::Mp<ES::V3d>(jconfig.getValue<std::array<double, 3>>("init-disp", 1).data());

  // scene scale
  double scale = jconfig.getDouble("scale", 1);

  // timestep
  double timestep = jconfig.getDouble("timestep", 1);

  // contact param
  double contactK = jconfig.getDouble("contact-stiffness", 1);
  int contactSamples = jconfig.getInt("contact-sample", 1);
  double fricCoeff = jconfig.getDouble("contact-friction-coeff", 1);
  double velEps = jconfig.getDouble("contact-vel-eps", 1);

  // solver param
  double solverEps = jconfig.getDouble("solver-eps", 1);
  int solverMaxIter = jconfig.getInt("solver-max-iter", 1);

  // damping
  std::array<double, 2> dampingParams = jconfig.getValue<std::array<double, 2>>("damping-params", 1);

  // material
  std::string material = jconfig.getString("elastic-material");
  pgo::SolidDeformationModel::DeformationModelElasticMaterial elasticMat;
  if (material == "stable-neo") {
    elasticMat = pgo::SolidDeformationModel::DeformationModelElasticMaterial::STABLE_NEO;
  }
  else if (material == "stvk-vol") {
    elasticMat = pgo::SolidDeformationModel::DeformationModelElasticMaterial::STVK_VOL;
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

  VolumetricMeshes::TetMesh tetMesh(tetMeshFilename.c_str());
  for (int vi = 0; vi < tetMesh.getNumVertices(); vi++) {
    tetMesh.setVertex(vi, tetMesh.getVertex(vi) * scale);
  }

  Mesh::TriMeshGeo surfaceMesh;
  if (surfaceMesh.load(surfaceMeshFilename) != true)
    return 1;

  for (int vi = 0; vi < surfaceMesh.numVertices(); vi++) {
    surfaceMesh.pos(vi) *= scale;
  }

  int surfn = surfaceMesh.numVertices();
  int surfn3 = surfn * 3;
  ES::VXd surfaceRestPositions(surfn3);
  for (int vi = 0; vi < surfn; vi++) {
    surfaceRestPositions.segment<3>(vi * 3) = surfaceMesh.pos(vi);
  }

  // initialize interpolation weights
  InterpolationCoordinates::BarycentricCoordinates bc(surfaceMesh.numVertices(), surfaceRestPositions.data(), &tetMesh);
  ES::SpMatD W = bc.generateInterpolationMatrix();

  // initialize fem
  std::shared_ptr<SolidDeformationModel::SimulationMesh> simMesh(SolidDeformationModel::loadTetMesh(&tetMesh));
  std::shared_ptr<SolidDeformationModel::DeformationModelManager> dmm = std::make_shared<SolidDeformationModel::DeformationModelManager>();

  dmm->setMesh(simMesh.get(), nullptr, nullptr);
  dmm->init(pgo::SolidDeformationModel::DeformationModelPlasticMaterial::VOLUMETRIC_DOF6, elasticMat, 1);

  std::vector<double> elementWeights(simMesh->getNumElements(), 1.0);
  std::shared_ptr<SolidDeformationModel::DeformationModelAssembler> assembler =
    std::make_shared<SolidDeformationModel::DeformationModelAssembler>(dmm, elementWeights.data());

  int n = simMesh->getNumVertices();
  int n3 = n * 3;
  int nele = simMesh->getNumElements();

  ES::VXd plasticity(nele * 6);
  ES::M3d I = ES::M3d::Identity();
  for (int ei = 0; ei < nele; ei++) {
    const SolidDeformationModel::PlasticModel3DDeformationGradient *pm =
      dynamic_cast<const SolidDeformationModel::PlasticModel3DDeformationGradient *>(dmm->getDeformationModel(ei)->getPlasticModel());
    if (!pm) {
      SPDLOG_LOGGER_ERROR(Logging::lgr(), "Plastic model is not of type PlasticModel3DDeformationGradient.");
      return 1;
    }
    pm->toParam(I.data(), plasticity.data() + ei * dmm->getNumPlasticParameters());
  }

  ES::VXd restPosition(n3);
  for (int vi = 0; vi < n; vi++) {
    double p[3];
    simMesh->getVertex(vi, p);
    restPosition.segment<3>(vi * 3) = ES::V3d(p[0], p[1], p[2]);
  }

  std::shared_ptr<SolidDeformationModel::DeformationModelEnergy> elasticEnergy =
    std::make_shared<SolidDeformationModel::DeformationModelEnergy>(assembler, &restPosition, 0);
  elasticEnergy->setPlasticParams(plasticity);

  ES::VXd zero(n3);
  zero.setZero();

  ES::SpMatD K;
  elasticEnergy->createHessian(K);
  elasticEnergy->hessian(zero, K);

  // attachments
  std::vector<std::shared_ptr<ConstraintPotentialEnergies::MultipleVertexPulling>> pullingEnergies;
  std::vector<ES::VXd> pullingTargets, pullingTargetRests;
  Mesh::TriMeshGeo tempMesh;
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
      tempMesh.addPos(tgtVertexPositions.segment<3>(vi * 3));
    }

    // initialize fixed constraints
    auto pullingEnergy = std::make_shared<ConstraintPotentialEnergies::MultipleVertexPulling>(K, restPosition.data(),
      (int)fixedVertices.size(), fixedVertices.data(), tgtVertexPositions.data(), nullptr, 1);
    pullingEnergy->setCoeff(attachmentCoeff);
    pullingEnergies.push_back(pullingEnergy);
    pullingTargets.push_back(tgtVertexPositions);
    pullingTargetRests.push_back(tgtVertexRests);
  }

  tempMesh.save("fv.obj");
  std::vector<std::string> kinematicObjectFilenames;
  std::vector<ES::V3d> kinematicObjectMovements;
  if (jconfig.exist("external-objects")) {
    auto jkinObjects = jconfig.handle()["external-objects"];
    for (const auto &jko : jkinObjects) {
      std::string koFilename = jko["filename"].get<std::string>();
      kinematicObjectFilenames.push_back(koFilename);
      kinematicObjectMovements.push_back(ES::Mp<ES::V3d>(jko["movement"].get<std::array<double, 3>>().data()));
    }
  }

  ES::SpMatD M;
  VolumetricMeshes::GenerateMassMatrix::computeMassMatrix(&tetMesh, M, true);

  // initialize gravity
  ES::VXd g(n3);
  for (int vi = 0; vi < n; vi++) {
    g.segment<3>(vi * 3) = extAcc;
  }

  ES::VXd fext(n3);
  ES::mv(M, g, fext);

  if (simType == "dynamic") {
    // initialize contact
    std::vector<Mesh::TriMeshGeo> kinematicObjects;
    for (const auto &filename : kinematicObjectFilenames) {
      kinematicObjects.emplace_back();
      if (kinematicObjects.back().load(filename) != true) {
        return 1;
      }

      for (int vi = 0; vi < kinematicObjects.back().numVertices(); vi++) {
        kinematicObjects.back().pos(vi) *= scale;
      }
    }

    std::vector<Mesh::TriMeshRef> kinematicObjectsRef;
    for (int i = 0; i < (int)kinematicObjects.size(); i++) {
      kinematicObjectsRef.emplace_back(kinematicObjects[i]);
    }

    std::shared_ptr<Contact::TriangleMeshExternalContactHandler> externalContactHandler;
    std::shared_ptr<Contact::TriangleMeshSelfContactHandler> selfCD;
    if (kinematicObjectsRef.size() && contactK > 0) {
      externalContactHandler = std::make_shared<Contact::TriangleMeshExternalContactHandler>(surfaceMesh.positions(), surfaceMesh.triangles(), n3,
        kinematicObjectsRef, contactSamples, &bc.getEmbeddingVertexIndices(), &bc.getEmbeddingWeights());
    }

    if (contactK > 0) {
      selfCD = std::make_shared<Contact::TriangleMeshSelfContactHandler>(
        surfaceMesh.positions(), surfaceMesh.triangles(), n3, contactSamples, &bc.getEmbeddingVertexIndices(), &bc.getEmbeddingWeights());
    }

    std::shared_ptr<Simulation::ImplicitBackwardEulerTimeIntegrator> intg =
      std::make_shared<Simulation::ImplicitBackwardEulerTimeIntegrator>(M, elasticEnergy,
        dampingParams[0], dampingParams[1], timestep, solverMaxIter, solverEps);

    // std::shared_ptr<Simulation::TRBDF2TimeIntegrator> intg =
    //   std::make_shared<Simulation::TRBDF2TimeIntegrator>(M, elasticEnergy,
    //     dampingParams[0], dampingParams[1], 1, timestep, solverMaxIter, solverEps);

    // #if defined(PGO_HAS_KNITRO)
    //     intg->setSolverOption(Simulation::TimeIntegratorSolverOption::SO_KNITRO);
    //     intg->setSolverConfigFile("config.opt");
    // #endif

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
      u.segment<3>(i * 3) = ES::V3d(initialDisp[0], initialDisp[1], initialDisp[2]);
    }

    if (!std::filesystem::exists(outputFolder)) {
      std::filesystem::create_directories(outputFolder);
    }

    int frameStart = 0;
    for (int framei = numSimSteps - 1; framei >= 0; framei--) {
      if (!std::filesystem::exists(fmt::format("{}/deform{:04d}.u", outputFolder, framei))) {
        std::cerr << "Frame " << framei << " not found." << std::endl;
        continue;
      }

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

    ES::mv(W, u, usurf);

    std::cout << frameStart << std::endl;
    std::cin.get();

    for (size_t eobji = 0; eobji < kinematicObjects.size(); eobji++) {
      ES::V3d movement = kinematicObjectMovements[eobji] / (numSimSteps - 1) * (frameStart + 1);
      for (int vi = 0; vi < kinematicObjects[eobji].numVertices(); vi++) {
        kinematicObjects[eobji].pos(vi) += movement;
      }
      externalContactHandler->updateExternalSurface(eobji, kinematicObjectsRef[eobji]);
    }

    for (int framei = frameStart + 1; framei < numSimSteps; framei++) {
      intg->clearGeneralImplicitForceModel();

      double ratio = (double)framei / (numSimSteps - 1);
      for (size_t pi = 0; pi < pullingEnergies.size(); pi++) {
        ES::VXd restTgt = pullingTargetRests[pi];
        ES::VXd curTgt = restTgt * (1 - ratio) + pullingTargets[pi] * ratio;
        pullingEnergies[pi]->setTargetPos(curTgt.data());

        std::cout << "Frame " << framei << ", attachment " << pi << " target: " << curTgt.transpose().head(3) << std::endl;
      }

      std::shared_ptr<Contact::PointPenetrationEnergy> extContactEnergy;
      Contact::PointPenetrationEnergyBuffer *extContactBuffer = nullptr;
      if (externalContactHandler) {
        externalContactHandler->execute(usurf.data());

        if (externalContactHandler->getNumCollidingSamples()) {
          extContactEnergy = externalContactHandler->buildContactEnergy();
          extContactBuffer = extContactEnergy->allocateBuffer();

          auto posFunc = [&restPosition](const EigenSupport::V3d &u, EigenSupport::V3d &p, int dofStart) {
            p = u + restPosition.segment<3>(dofStart);
          };

          auto lastPosFunc = [&restPosition, &u](const EigenSupport::V3d &x, EigenSupport::V3d &p, int dofStart) {
            p = u.segment<3>(dofStart) + restPosition.segment<3>(dofStart);
          };

          extContactEnergy->setComputePosFunction(posFunc);
          extContactEnergy->setBuffer(extContactBuffer);
          extContactEnergy->setCoeff(contactK);

          extContactEnergy->setFrictionCoeff(fricCoeff);
          extContactEnergy->setComputeLastPosFunction(lastPosFunc);
          extContactEnergy->setVelEps(velEps);
          extContactEnergy->setTimestep(timestep);

          // if (externalContactHandler->getNumCollidingSamples() > 500) {
          //   NonlinearOptimization::FiniteDifference fd(NonlinearOptimization::FiniteDifference::M_FIVE_POINT, 1e-6);
          //   double err[2];

          //  fd.testEnergy(extContactEnergy, true, true, -1, u.data(), 500, err, err + 1);

          //  std::cout << "External contact energy FD test, max rel error: grad " << err[0] << ", hess " << err[1] << std::endl;
          //  exit(1);
          //}

          intg->addGeneralImplicitForceModel(extContactEnergy, 0, 0);
        }
      }

      std::shared_ptr<Contact::PointTrianglePairCouplingEnergyWithCollision> selfContactEnergy;
      Contact::PointTrianglePairCouplingEnergyWithCollisionBuffer *selfContactEnergyBuf = nullptr;

      if (selfCD) {
        selfCD->execute(usurf.data());

        if (selfCD->getCollidingTrianglePair().size() > 0) {
          /*Mesh::TriMeshGeo m;
          for (int i = 0; i < (int)selfCD->getCollidingTrianglePair().size(); i++) {
            const auto &ccdData = selfCD->getCollidingTrianglePair()[i];

            m.addMesh(Mesh::createSingleTriangleMesh(
              surfaceMesh.pos(ccdData.triA, 0),
              surfaceMesh.pos(ccdData.triA, 1),
              surfaceMesh.pos(ccdData.triA, 2)));

            m.addMesh(Mesh::createSingleTriangleMesh(
              surfaceMesh.pos(ccdData.triB, 0),
              surfaceMesh.pos(ccdData.triB, 1),
              surfaceMesh.pos(ccdData.triB, 2)));
          }
          m.save("cd.obj");
          exit(1);*/

          selfCD->handleContactDCD(0, 100);

          selfContactEnergy = selfCD->buildContactEnergy();
          selfContactEnergy->setToPosFunction([&restPosition](const ES::V3d &x, ES::V3d &p, int offset) {
            p = x + restPosition.segment<3>(offset);
          });

          selfContactEnergy->setToLastPosFunction([&restPosition, &u](const ES::V3d &x, ES::V3d &p, int offset) {
            p = restPosition.segment<3>(offset) + u.segment<3>(offset);
          });

          selfContactEnergyBuf = selfContactEnergy->allocateBuffer();
          selfContactEnergy->setBuffer(selfContactEnergyBuf);
          selfContactEnergy->setCoeff(contactK);
          selfContactEnergy->computeClosestPosition(u.data());

          selfContactEnergy->setFrictionCoeff(fricCoeff);
          selfContactEnergy->setTimestep(timestep);
          selfContactEnergy->setVelEps(velEps);

          intg->addGeneralImplicitForceModel(selfContactEnergy, 0, 0);
        }
      }

      intg->setqState(u, uvel, uacc);

      intg->doTimestep(1, 2, 1);

      intg->getq(u);
      intg->getqvel(uvel);
      intg->getqacc(uacc);

      //double Ec = extContactEnergy ? extContactEnergy->func(u) : 0;
      //double Eelastic = elasticEnergy->func(u);
      //double Emass = 0.5 * uvel.transpose() * M * uvel;
      //double Etotal = Ec + Eelastic + Emass;
      //std::cout << "Frame " << framei << ": Eelastic = " << Eelastic << ", Ec = " << Ec << ", Emass = " << Emass << ", Etotal = " << Etotal << std::endl;

      if (extContactBuffer && extContactEnergy) {
        extContactEnergy->freeBuffer(extContactBuffer);
      }

      if (selfContactEnergyBuf && selfContactEnergy) {
        selfContactEnergy->freeBuffer(selfContactEnergyBuf);
      }

      ES::mv(W, u, usurf);

      if (framei % frameGap == 0) {
        ES::VXd psurf = surfaceRestPositions + usurf;

        Mesh::TriMeshGeo mesh = surfaceMesh;
        for (int vi = 0; vi < mesh.numVertices(); vi++) {
          mesh.pos(vi) = psurf.segment<3>(vi * 3) / scale;
        }
        mesh.save(fmt::format("{}/ret{:04d}.obj", outputFolder, framei / frameGap));
      }

      for (size_t eobji = 0; eobji < kinematicObjects.size(); eobji++) {
        ES::V3d movement = kinematicObjectMovements[eobji] / (numSimSteps - 1);
        for (int vi = 0; vi < kinematicObjects[eobji].numVertices(); vi++) {
          kinematicObjects[eobji].pos(vi) += movement;
        }
        externalContactHandler->updateExternalSurface(eobji, kinematicObjectsRef[eobji]);
      }

      ES::MXd uMat(n3, 3);
      uMat.col(0) = u;
      uMat.col(1) = uvel;
      uMat.col(2) = uacc;

      ES::writeMatrix(fmt::format("{}/deform{:04d}.u", outputFolder, framei).c_str(), uMat);
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

    ES::VXd xsurf(surfn3);
    ES::mv(W, x, xsurf);

    Mesh::TriMeshGeo meshOut = surfaceMesh;
    for (int vi = 0; vi < meshOut.numVertices(); vi++) {
      meshOut.pos(vi) = xsurf.segment<3>(vi * 3);
    }
    meshOut.save(outputFolder);
  }

  return 0;
}
