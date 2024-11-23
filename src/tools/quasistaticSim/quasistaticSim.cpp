#include "tetMesh.h"
#include "EigenSupport.h"
#include "generateTetMeshMatrix.h"
#include "simulationMesh.h"
#include "deformationModelManager.h"
#include "tetMeshDeformationModel.h"
#include "pgoLogging.h"
#include "basicIO.h"
#include "deformationModelAssemblerElastic.h"
#include "deformationModelEnergy.h"
#include "plasticModel.h"
#include "potentialEnergies.h"
#include "multiVertexPullingSoftConstraints.h"
#include "NewtonRaphsonSolver.h"
#include "linearPotentialEnergy.h"
#include "triMeshGeo.h"
#include "barycentricCoordinates.h"

int main(int argc, char* argv[])
{
  using namespace pgo;
  namespace ES = EigenSupport;

  if (argc != 5) {
    std::cerr << argv[0] << " <tet mesh> <obj mesh> <fixed vertices> <out obj mesh>" << std::endl;
    return 1;
  }

  pgo::Logging::init();

  SPDLOG_LOGGER_INFO(Logging::lgr(), "Loading {}...", argv[1]);
  VolumetricMeshes::TetMesh tetMesh(argv[1]);

  SPDLOG_LOGGER_INFO(Logging::lgr(), "Loading {}...", argv[2]);
  Mesh::TriMeshGeo surfaceMesh;

  if (surfaceMesh.load(argv[2]) == false) {
    SPDLOG_LOGGER_ERROR(Logging::lgr(), "Cannot load {}", argv[2]);
    return 1;
  }

  SPDLOG_LOGGER_INFO(Logging::lgr(), "Loading FEM...");

  std::shared_ptr<SolidDeformationModel::SimulationMesh> simMesh(SolidDeformationModel::loadTetMesh(&tetMesh));
  std::shared_ptr<SolidDeformationModel::DeformationModelManager> dmm = std::make_shared<SolidDeformationModel::DeformationModelManager>();

  dmm->setMesh(simMesh.get(), nullptr, nullptr);
  dmm->init(pgo::SolidDeformationModel::DeformationModelPlasticMaterial::VOLUMETRIC_DOF6,
    pgo::SolidDeformationModel::DeformationModelElasticMaterial::STABLE_NEO, nullptr);

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
  std::shared_ptr<SolidDeformationModel::DeformationModelEnergy> elasticEnergy = std::make_shared<SolidDeformationModel::DeformationModelEnergy>(assembler, nullptr, 0);

  ES::VXd restPosition(n3);
  for (int vi = 0; vi < n; vi++) {
    double p[3];
    simMesh->getVertex(vi, p);
    restPosition.segment<3>(vi * 3) = ES::V3d(p[0], p[1], p[2]);
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

  // attachments
  std::vector<int> fixedVertices;
  if (BasicIO::read1DText(argv[3], std::back_inserter(fixedVertices)) != 0) {
    return 1;
  }
  std::sort(fixedVertices.begin(), fixedVertices.end());

  ES::VXd tgtVertexPositions(fixedVertices.size() * 3);
  for (int vi = 0; vi < (int)fixedVertices.size(); vi++) {
    tgtVertexPositions.segment<3>(vi * 3) = restPosition.segment<3>(fixedVertices[vi] * 3);
  }

  ES::SpMatD K;
  elasticEnergy->createHessian(K);
  elasticEnergy->hessian(restPosition, K);

  std::shared_ptr<ConstraintPotentialEnergies::MultipleVertexPulling> pullingEnergy =
    std::make_shared<ConstraintPotentialEnergies::MultipleVertexPulling>(K, restPosition.data(),
      (int)fixedVertices.size(), fixedVertices.data(), tgtVertexPositions.data(), nullptr, 0);

  ES::VXd vertexMasses(n);
  vertexMasses.setZero();
  for (int ele = 0; ele < nele; ele++) {
    double mass = tetMesh.getElementVolume(ele) * 1000.0;
    for (int j = 0; j < 4; j++) {
      vertexMasses[simMesh->getVertexIndex(ele, j)] += mass * 0.25;
    }
  }

  // gravity
  ES::VXd fext(n3);
  for (int vi = 0; vi < n; vi++) {
    fext.segment<3>(vi * 3) = ES::V3d(0, 9.8, 0) * vertexMasses[vi];
  }

  std::shared_ptr<PredefinedPotentialEnergies::LinearPotentialEnergy> externalForcesEnergy = std::make_shared<PredefinedPotentialEnergies::LinearPotentialEnergy>(fext);

  std::shared_ptr<NonlinearOptimization::PotentialEnergies> energyAll = std::make_shared<NonlinearOptimization::PotentialEnergies>(n3);
  energyAll->addPotentialEnergy(elasticEnergy);
  energyAll->addPotentialEnergy(pullingEnergy, 1e5);
  energyAll->addPotentialEnergy(externalForcesEnergy);
  energyAll->init();

  NonlinearOptimization::NewtonRaphsonSolver::SolverParam solverParam;

  ES::VXd x = restPosition;
  NonlinearOptimization::NewtonRaphsonSolver solver(x.data(), solverParam, energyAll, std::vector<int>(), nullptr);
  solver.solve(x.data(), 200, 1e-6, 3);

  ES::VXd xsurf(surfn3);
  ES::mv(W, x, xsurf);

  Mesh::TriMeshGeo meshOut = surfaceMesh;
  for (int vi = 0; vi < meshOut.numVertices(); vi++) {
    meshOut.pos(vi) = xsurf.segment<3>(vi * 3);
  }
  meshOut.save(argv[4]);



  return 0;
}