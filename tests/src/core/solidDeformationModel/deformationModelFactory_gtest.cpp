#include <gtest/gtest.h>

#include "deformationModelFactory.h"

#include "deformationModel.h"
#include "deformationModelAssembler.h"
#include "deformationModelEnergy.h"
#include "deformationModelManager.h"
#include "plasticModel3DDeformationGradient.h"
#include "simulationMesh.h"
#include "tetMesh.h"
#include "pgoLogging.h"

#include <memory>
#include <vector>

namespace
{
namespace ES = pgo::EigenSupport;
using namespace pgo::SolidDeformationModel;

constexpr const char *kTorusVegPath = LIBPGO_TEST_TORUS_VEG;

// Rebuild the FEM energy the long way (the exact chain makeDeformationModel collapses),
// so the test fails if the facade ever drifts from the documented construction protocol.
// The mesh must live as long as the energy (manager stores a raw borrow pointer), so we
// keep it alive outside and use the borrow setMesh.
std::shared_ptr<DeformationModelEnergy> buildEnergyManually(
  const pgo::VolumetricMeshes::TetMesh &tetMesh, ES::VXd &restPositionOut,
  std::unique_ptr<SimulationMesh> &meshOut)
{
  std::unique_ptr<SimulationMesh> mesh = loadTetMesh(&tetMesh);

  const int nele = mesh->getNumElements();
  const int nvtx = mesh->getNumVertices();

  restPositionOut.resize(nvtx * 3);
  for (int vi = 0; vi < nvtx; vi++) {
    double p[3];
    mesh->getVertex(vi, p);
    restPositionOut.segment<3>(vi * 3) = ES::V3d(p[0], p[1], p[2]);
  }

  auto dmm = std::make_unique<DeformationModelManager>();
  dmm->setMesh(std::move(mesh), nullptr, nullptr);  // dmm now owns the mesh
  dmm->init(DeformationModelPlasticMaterial::VOLUMETRIC_DOF6, DeformationModelElasticMaterial::STABLE_NEO);
  dmm->setEnforceSPD(1);

  const int numPlasticParams = dmm->getNumPlasticParameters();
  ES::VXd plasticParams(static_cast<Eigen::Index>(nele) * numPlasticParams);
  plasticParams.setZero();
  const ES::M3d identity = ES::M3d::Identity();
  for (int ei = 0; ei < nele; ei++) {
    const auto *pm = dynamic_cast<const PlasticModel3DDeformationGradient *>(
      dmm->getDeformationModel(ei)->getPlasticModel());
    pm->toParam(identity.data(), plasticParams.data() + ei * numPlasticParams);
  }

  std::vector<double> elementWeights(nele, 1.0);
  auto assembler = std::make_unique<DeformationModelAssembler>(std::move(dmm), elementWeights.data());

  auto energy = std::make_shared<DeformationModelEnergy>(std::move(assembler), &restPositionOut, 0);
  energy->setEnableMaterialMaxStep(true);
  energy->setPlasticParams(plasticParams);

  // Return ownership of the mesh to the caller (the energy's chain borrowed it).
  // With the unique spine, the mesh lives inside the manager, so we don't need to
  // keep it separately. The meshOut is kept for backward compatibility with the
  // old test signature; the test now only needs the energy.
  meshOut = nullptr;
  return energy;
}

ES::VXd perturbed(const ES::VXd &rest)
{
  ES::VXd x = rest;
  for (Eigen::Index i = 0; i < x.size(); i++)
    x[i] += 1e-3 * std::sin(0.7 * static_cast<double>(i));
  return x;
}
}  // namespace

TEST(DeformationModelFactoryGTest, MakeDeformationModelMatchesManualChain)
{
  pgo::Logging::init();

  pgo::VolumetricMeshes::TetMesh tetMesh(kTorusVegPath);

  ES::VXd manualRest;
  std::unique_ptr<SimulationMesh> manualMesh;
  auto manual = buildEnergyManually(tetMesh, manualRest, manualMesh);

  DeformationModelBundle bundle = makeDeformationModel(
    tetMesh, DeformationModelElasticMaterial::STABLE_NEO, DeformationModelPlasticMaterial::VOLUMETRIC_DOF6);

  ASSERT_NE(bundle.energy, nullptr);
  ASSERT_EQ(bundle.energy->getNumDOFs(), manual->getNumDOFs());
  EXPECT_LT((bundle.restPosition - manualRest).cwiseAbs().maxCoeff(), 1e-15);

  for (const ES::VXd &x : { ES::VXd(manualRest), perturbed(manualRest) }) {
    const double fFacade = bundle.energy->func(x);
    const double fManual = manual->func(x);
    EXPECT_NEAR(fFacade, fManual, 1e-9 * (1.0 + std::abs(fManual)));

    ES::VXd gFacade = ES::VXd::Zero(x.size());
    ES::VXd gManual = ES::VXd::Zero(x.size());
    bundle.energy->gradient(x, gFacade);
    manual->gradient(x, gManual);
    EXPECT_LT((gFacade - gManual).cwiseAbs().maxCoeff(), 1e-9 * (1.0 + gManual.cwiseAbs().maxCoeff()));

    ES::SpMatD hFacade, hManual;
    bundle.energy->createHessian(hFacade);
    bundle.energy->hessian(x, hFacade);
    manual->createHessian(hManual);
    manual->hessian(x, hManual);
    EXPECT_EQ(hFacade.nonZeros(), hManual.nonZeros());
    const double hDiff = (ES::MXd(hFacade) - ES::MXd(hManual)).cwiseAbs().maxCoeff();
    EXPECT_LT(hDiff, 1e-9 * (1.0 + ES::MXd(hManual).cwiseAbs().maxCoeff()));
  }
}

TEST(DeformationModelFactoryGTest, MakeSimulationMeshReturnsSharedOwnership)
{
  pgo::Logging::init();

  pgo::VolumetricMeshes::TetMesh tetMesh(kTorusVegPath);
  std::shared_ptr<SimulationMesh> mesh = makeSimulationMesh(tetMesh);
  ASSERT_NE(mesh, nullptr);
  EXPECT_GT(mesh->getNumVertices(), 0);
  EXPECT_GT(mesh->getNumElements(), 0);
}
