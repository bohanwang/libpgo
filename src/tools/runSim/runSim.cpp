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
#include "basicIO.h"
#include "deformationModelAssemblerElastic.h"
#include "deformationModelEnergy.h"
#include "plasticModel.h"
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

#include <argparse/argparse.hpp>

#include <iostream>

int main(int argc, char *argv[])
{
  using namespace pgo;
  using namespace pgo::EigenSupport;
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

  // tet mesh filename
  std::string tetMeshFilename = jconfig.getString("tet-mesh", 1);

  // surface mesh filename
  std::string surfaceMeshFilename = jconfig.getString("surface-mesh", 1);

  // fixed vertices filename
  std::string fixedVerticesFilename = jconfig.getString("fixed-vertices", 1);

  // external objects filenames
  std::vector<std::string> kinematicObjectFilenames = jconfig.getVectorString("external-objects", 1);

  // external acceleration
  ES::V3d extAcc = ES::Mp<ES::V3d>(jconfig.getValue<std::array<double, 3>>("g", 1).data());

  // initial velocity
  ES::V3d initialVel = ES::Mp<ES::V3d>(jconfig.getValue<std::array<double, 3>>("init-vel", 1).data());

  // scene scale
  double scale = jconfig.getDouble("scale", 1);

  // timestep
  double timestep = jconfig.getDouble("timestep", 1);

  // contact param
  double contactK = jconfig.getDouble("contact-stiffness", 1);
  int contactSamples = jconfig.getInt("contact-sample", 1);
  double fricCoeff = jconfig.getDouble("contact-friction-coeff", 1);
  double velEps = jconfig.getDouble("contact-vel-eps", 1);

  double attachmentCoeff = jconfig.getDouble("attachment-stiffness", 1);

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
  dmm->init(pgo::SolidDeformationModel::DeformationModelPlasticMaterial::VOLUMETRIC_DOF6,
    elasticMat, nullptr);

  std::vector<double> elementWeights(simMesh->getNumElements(), 1.0);
  std::shared_ptr<SolidDeformationModel::DeformationModelAssemblerElastic> assembler =
    std::make_shared<SolidDeformationModel::DeformationModelAssemblerElastic>(dmm.get(), nullptr, 0, elementWeights.data());

  int n = simMesh->getNumVertices();
  int n3 = n * 3;
  int nele = simMesh->getNumElements();
  ES::VXd plasticity(nele * 6);
  ES::M3d I = ES::M3d::Identity();
  for (int ei = 0; ei < nele; ei++) {
    dmm->getDeformationModel(ei)->getPlasticModel()->toParam(I.data(), plasticity.data() + ei * dmm->getNumPlasticParameters());
  }
  assembler->enableFullspaceScale(1);
  assembler->setFixedParameters(plasticity.data(), (int)plasticity.size());

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
  std::shared_ptr<ConstraintPotentialEnergies::MultipleVertexPulling> pullingEnergy;

  if (fixedVerticesFilename.length()) {
    std::vector<int> fixedVertices;
    if (BasicIO::read1DText(fixedVerticesFilename.c_str(), std::back_inserter(fixedVertices)) != 0) {
      return 1;
    }
    std::sort(fixedVertices.begin(), fixedVertices.end());

    ES::VXd tgtVertexPositions(fixedVertices.size() * 3);
    for (int vi = 0; vi < (int)fixedVertices.size(); vi++) {
      tgtVertexPositions.segment<3>(vi * 3) = restPosition.segment<3>(fixedVertices[vi] * 3);
    }

    // initialize fixed constraints
    pullingEnergy = std::make_shared<ConstraintPotentialEnergies::MultipleVertexPulling>(K, restPosition.data(),
      (int)fixedVertices.size(), fixedVertices.data(), tgtVertexPositions.data(), nullptr, 1);
    pullingEnergy->setCoeff(attachmentCoeff);

    Mesh::TriMeshGeo fv;
    for (int vi = 0; vi < (int)fixedVertices.size(); vi++) {
      fv.addPos(restPosition.segment<3>(fixedVertices[vi] * 3));
    }
    fv.save("fv.obj");
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

#if defined(PGO_HAS_KNITRO)
    intg->setSolverOption(Simulation::TimeIntegratorSolverOption::SO_KNITRO);
    intg->setSolverConfigFile("config.opt");
#endif

    if (pullingEnergy)
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
      intg->getq(uvel);
      intg->getq(uacc);

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
    }
  }
  else if (simType == "static") {
    std::shared_ptr<PredefinedPotentialEnergies::LinearPotentialEnergy> externalForcesEnergy = std::make_shared<PredefinedPotentialEnergies::LinearPotentialEnergy>(fext);

    std::shared_ptr<NonlinearOptimization::PotentialEnergies> energyAll = std::make_shared<NonlinearOptimization::PotentialEnergies>(n3);
    energyAll->addPotentialEnergy(elasticEnergy);
    energyAll->addPotentialEnergy(pullingEnergy, 1.0);
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