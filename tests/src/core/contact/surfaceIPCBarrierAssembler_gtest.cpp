#include <gtest/gtest.h>

#include "ipc/broadPhase/surfaceIPCSelfBroadPhase.h"
#include "ipc/core/surfaceIPCBarrierAssembler.h"
#include "ipc/topology/surfaceIPCTopology.h"
#include "ipc/core/surfaceIPCCore.h"
#include "ipc/geometry/ipcDistancePrimitives.h"

#include "testCIPCHelpers.h"

namespace
{
namespace ES = pgo::EigenSupport;
using pgo::Contact::CIPC::EEPair;
using pgo::Contact::CIPC::PTPair;
using pgo::Contact::CIPC::SurfaceIPCBarrierAssembler;
using pgo::Contact::CIPC::SurfaceIPCCore;
using pgo::Contact::CIPC::SurfaceIPCSelfBroadPhase;
using pgo::Contact::CIPC::SurfaceIPCTopology;
using pgo::Contact::CIPCTest::flattenPositions;
using pgo::Contact::CIPCTest::finiteDifferenceGradient;
using pgo::Contact::CIPCTest::finiteDifferenceHessian;
using pgo::Contact::CIPCTest::makeTwoTriangleMesh;
using pgo::Contact::CIPCTest::relativeError;
using pgo::Contact::CIPCTest::sparseToDense;
namespace distance = pgo::Contact::CIPC::distance;

constexpr double kDhat = 0.1;
constexpr double kKappa = 1.0;
constexpr double kFDStep = 1e-6;

ES::VXd makePointTriangleTransitionState(double pointY)
{
  ES::VXd x(12);
  x << 0.5, pointY, 0.03,
    0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0;
  return x;
}

ES::VXd makeEdgeEdgeTransitionState(double edgeBX)
{
  ES::VXd x(12);
  x << 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    edgeBX, -0.5, 0.03,
    edgeBX, 0.5, 0.03;
  return x;
}

ES::VXd makeMollifiedEdgeEdgeState()
{
  ES::VXd x(12);
  x << 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    0.0, -0.02, 0.03,
    1.0, 0.02, 0.03;
  return x;
}

const std::vector<PTPair> kSinglePointTrianglePair = { { 0, 1, 2, 3, 1.0 } };
const std::vector<EEPair> kSingleEdgeEdgePair = { { 0, 1, 2, 3, 1.0 } };
}  // namespace

TEST(SurfaceIPCBarrierAssemblerGTest, HelperMatchesSurfaceIPCCoreAssembly)
{
  const auto [V, F] = makeTwoTriangleMesh();
  const ES::VXd x = flattenPositions(V);

  SurfaceIPCTopology topology;
  topology.setMesh(V, F);

  constexpr double dhat = 0.1;
  constexpr double kappa = 1.0;
  constexpr double epsEE = 0.0;
  std::vector<PTPair> ptPairs;
  std::vector<EEPair> eePairs;
  SurfaceIPCSelfBroadPhase().buildPairs(topology, x, dhat, ptPairs, eePairs);

  const SurfaceIPCBarrierAssembler assembler;
  const double helperEnergy = assembler.computeEnergy(x, ptPairs, eePairs, topology.numVerts, dhat, kappa, epsEE);
  ES::VXd helperGradient(x.size());
  assembler.computeGradient(x, ptPairs, eePairs, topology.numVerts, dhat, kappa, epsEE, helperGradient);
  ES::SpMatD helperHessian;
  assembler.computeHessian(x, ptPairs, eePairs, topology.numVerts, dhat, kappa, epsEE, true, helperHessian);

  SurfaceIPCCore core;
  SurfaceIPCCore::Parameters params;
  params.dhat = dhat;
  params.kappa = kappa;
  params.eps_ee = epsEE;
  params.slackness = 0.9;
  core.setParameters(params);
  core.setMesh(V, F);

  ES::VXd coreGradient(x.size());
  core.computeGradient(x, coreGradient);
  ES::SpMatD coreHessian;
  core.computeHessian(x, coreHessian);

  EXPECT_NEAR(helperEnergy, core.computeEnergy(x), 1e-12);
  EXPECT_LT(relativeError(helperGradient, coreGradient), 1e-12);
  EXPECT_LT(relativeError(sparseToDense(helperHessian), sparseToDense(coreHessian)), 1e-12);
}

TEST(SurfaceIPCBarrierAssemblerGTest, PointTriangleFeatureTransitionPreservesGradientAndOneSidedHessians)
{
  const SurfaceIPCBarrierAssembler assembler;
  const std::vector<EEPair> noEdgePairs;
  const auto energy = [&](const ES::VXd &x) {
    return assembler.computeEnergy(x, kSinglePointTrianglePair, noEdgePairs, 4, kDhat, kKappa, 0.0);
  };
  const auto gradient = [&](const ES::VXd &x) {
    ES::VXd g = ES::VXd::Zero(x.size());
    assembler.computeGradient(x, kSinglePointTrianglePair, noEdgePairs, 4, kDhat, kKappa, 0.0, g);
    return g;
  };

  const ES::VXd transition = makePointTriangleTransitionState(0.0);
  EXPECT_EQ(distance::classifyPT(transition.segment<3>(0), transition.segment<3>(3),
              transition.segment<3>(6), transition.segment<3>(9)),
    pgo::Contact::CIPC::PTDistType::PT);
  const ES::VXd gradientAtTransition = gradient(transition);
  const ES::VXd fdGradientAtTransition = finiteDifferenceGradient(energy, transition, kFDStep);
  EXPECT_LT(relativeError(gradientAtTransition, fdGradientAtTransition), 2e-5);

  const ES::VXd faceSide = makePointTriangleTransitionState(1e-3);
  const ES::VXd edgeSide = makePointTriangleTransitionState(-1e-3);
  EXPECT_EQ(distance::classifyPT(faceSide.segment<3>(0), faceSide.segment<3>(3),
              faceSide.segment<3>(6), faceSide.segment<3>(9)),
    pgo::Contact::CIPC::PTDistType::PT);
  EXPECT_EQ(distance::classifyPT(edgeSide.segment<3>(0), edgeSide.segment<3>(3),
              edgeSide.segment<3>(6), edgeSide.segment<3>(9)),
    pgo::Contact::CIPC::PTDistType::PE_PT0T1);

  for (const ES::VXd *x : { &faceSide, &edgeSide }) {
    const ES::VXd analyticGradient = gradient(*x);
    const ES::VXd fdGradient = finiteDifferenceGradient(energy, *x, kFDStep);
    EXPECT_LT(relativeError(analyticGradient, fdGradient), 2e-5);

    ES::SpMatD hessian;
    assembler.computeHessian(*x, kSinglePointTrianglePair, noEdgePairs, 4,
      kDhat, kKappa, 0.0, false, hessian);
    const ES::MXd fdHessian = finiteDifferenceHessian(gradient, *x, kFDStep);
    const ES::MXd analyticHessian = sparseToDense(hessian);
    EXPECT_LT(relativeError(analyticHessian, analyticHessian.transpose()), 1e-12);
    EXPECT_LT(relativeError(analyticHessian, fdHessian), 2e-5);
  }
}

TEST(SurfaceIPCBarrierAssemblerGTest, EdgeEdgeFeatureTransitionPreservesGradientAndOneSidedHessians)
{
  const SurfaceIPCBarrierAssembler assembler;
  const std::vector<PTPair> noPointTrianglePairs;
  const auto energy = [&](const ES::VXd &x) {
    return assembler.computeEnergy(x, noPointTrianglePairs, kSingleEdgeEdgePair, 4, kDhat, kKappa, 0.0);
  };
  const auto gradient = [&](const ES::VXd &x) {
    ES::VXd g = ES::VXd::Zero(x.size());
    assembler.computeGradient(x, noPointTrianglePairs, kSingleEdgeEdgePair, 4, kDhat, kKappa, 0.0, g);
    return g;
  };

  const ES::VXd transition = makeEdgeEdgeTransitionState(0.0);
  EXPECT_EQ(distance::classifyEE(transition.segment<3>(0), transition.segment<3>(3),
              transition.segment<3>(6), transition.segment<3>(9)),
    pgo::Contact::CIPC::EEDistType::PE_Ea0_Eb);
  const ES::VXd gradientAtTransition = gradient(transition);
  const ES::VXd fdGradientAtTransition = finiteDifferenceGradient(energy, transition, kFDStep);
  EXPECT_LT(relativeError(gradientAtTransition, fdGradientAtTransition), 2e-5);

  const ES::VXd edgeEdgeSide = makeEdgeEdgeTransitionState(1e-3);
  const ES::VXd pointEdgeSide = makeEdgeEdgeTransitionState(-1e-3);
  EXPECT_EQ(distance::classifyEE(edgeEdgeSide.segment<3>(0), edgeEdgeSide.segment<3>(3),
              edgeEdgeSide.segment<3>(6), edgeEdgeSide.segment<3>(9)),
    pgo::Contact::CIPC::EEDistType::EE);
  EXPECT_EQ(distance::classifyEE(pointEdgeSide.segment<3>(0), pointEdgeSide.segment<3>(3),
              pointEdgeSide.segment<3>(6), pointEdgeSide.segment<3>(9)),
    pgo::Contact::CIPC::EEDistType::PE_Ea0_Eb);

  for (const ES::VXd *x : { &edgeEdgeSide, &pointEdgeSide }) {
    const ES::VXd analyticGradient = gradient(*x);
    const ES::VXd fdGradient = finiteDifferenceGradient(energy, *x, kFDStep);
    EXPECT_LT(relativeError(analyticGradient, fdGradient), 2e-5);

    ES::SpMatD hessian;
    assembler.computeHessian(*x, noPointTrianglePairs, kSingleEdgeEdgePair, 4,
      kDhat, kKappa, 0.0, false, hessian);
    const ES::MXd fdHessian = finiteDifferenceHessian(gradient, *x, kFDStep);
    const ES::MXd analyticHessian = sparseToDense(hessian);
    EXPECT_LT(relativeError(analyticHessian, analyticHessian.transpose()), 1e-12);
    EXPECT_LT(relativeError(analyticHessian, fdHessian), 2e-5);
  }
}

TEST(SurfaceIPCBarrierAssemblerGTest, MollifiedEdgeEdgeBarrierMatchesFiniteDifferences)
{
  const SurfaceIPCBarrierAssembler assembler;
  const std::vector<PTPair> noPointTrianglePairs;
  constexpr double kMollifierEpsilon = 1e-2;
  const ES::VXd x = makeMollifiedEdgeEdgeState();
  const double mollifier = distance::eeMollifier(x.segment<3>(0), x.segment<3>(3),
    x.segment<3>(6), x.segment<3>(9), kMollifierEpsilon);
  EXPECT_GT(mollifier, 0.0);
  EXPECT_LT(mollifier, 1.0);

  const auto energy = [&](const ES::VXd &state) {
    return assembler.computeEnergy(state, noPointTrianglePairs, kSingleEdgeEdgePair, 4,
      kDhat, kKappa, kMollifierEpsilon);
  };
  const auto gradient = [&](const ES::VXd &state) {
    ES::VXd g = ES::VXd::Zero(state.size());
    assembler.computeGradient(state, noPointTrianglePairs, kSingleEdgeEdgePair, 4,
      kDhat, kKappa, kMollifierEpsilon, g);
    return g;
  };
  const ES::VXd analyticGradient = gradient(x);
  const ES::VXd fdGradient = finiteDifferenceGradient(energy, x, kFDStep);
  EXPECT_LT(relativeError(analyticGradient, fdGradient), 2e-5);

  ES::SpMatD hessian;
  assembler.computeHessian(x, noPointTrianglePairs, kSingleEdgeEdgePair, 4,
    kDhat, kKappa, kMollifierEpsilon, false, hessian);
  const ES::MXd analyticHessian = sparseToDense(hessian);
  const ES::MXd fdHessian = finiteDifferenceHessian(gradient, x, kFDStep);
  EXPECT_LT(relativeError(analyticHessian, analyticHessian.transpose()), 1e-12);
  EXPECT_LT(relativeError(analyticHessian, fdHessian), 2e-5);
}
