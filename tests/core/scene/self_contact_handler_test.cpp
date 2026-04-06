#include "pointTrianglePairBarrierEnergy.h"
#include "triangleMeshSelfContactHandler.h"
#include "triMeshGeo.h"
#include "pgoLogging.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace pgo::Contact::test {

namespace {

Mesh::TriMeshGeo makeSeparatedParallelTriangleMesh(double gap) {
    return Mesh::TriMeshGeo(
        std::vector<EigenSupport::V3d>{
            EigenSupport::V3d(0.0, 0.0, 0.0),
            EigenSupport::V3d(1.0, 0.0, 0.0),
            EigenSupport::V3d(0.0, 1.0, 0.0),
            EigenSupport::V3d(0.0, 0.0, gap),
            EigenSupport::V3d(1.0, 0.0, gap),
            EigenSupport::V3d(0.0, 1.0, gap),
        },
        std::vector<EigenSupport::V3i>{
            EigenSupport::V3i(0, 1, 2),
            EigenSupport::V3i(3, 4, 5),
        });
}

Mesh::TriMeshGeo makeIntersectingTriangleMesh() {
    return Mesh::TriMeshGeo(
        std::vector<EigenSupport::V3d>{
            EigenSupport::V3d(0.0, 0.0, 0.0),
            EigenSupport::V3d(1.0, 0.0, 0.0),
            EigenSupport::V3d(0.0, 1.0, 0.0),
            EigenSupport::V3d(0.2, 0.2, -0.2),
            EigenSupport::V3d(0.8, 0.2, 0.2),
            EigenSupport::V3d(0.2, 0.8, 0.2),
        },
        std::vector<EigenSupport::V3i>{
            EigenSupport::V3i(0, 1, 2),
            EigenSupport::V3i(3, 4, 5),
        });
}

Mesh::TriMeshGeo makeThreeLayerTriangleMesh(double gapLowerMiddle, double gapMiddleUpper) {
    const double middleZ = gapLowerMiddle;
    const double upperZ  = gapLowerMiddle + gapMiddleUpper;

    return Mesh::TriMeshGeo(
        std::vector<EigenSupport::V3d>{
            EigenSupport::V3d(0.0, 0.0, 0.0),
            EigenSupport::V3d(1.0, 0.0, 0.0),
            EigenSupport::V3d(0.0, 1.0, 0.0),
            EigenSupport::V3d(0.0, 0.0, middleZ),
            EigenSupport::V3d(1.0, 0.0, middleZ),
            EigenSupport::V3d(0.0, 1.0, middleZ),
            EigenSupport::V3d(0.0, 0.0, upperZ),
            EigenSupport::V3d(1.0, 0.0, upperZ),
            EigenSupport::V3d(0.0, 1.0, upperZ),
        },
        std::vector<EigenSupport::V3i>{
            EigenSupport::V3i(0, 1, 2),
            EigenSupport::V3i(3, 4, 5),
            EigenSupport::V3i(6, 7, 8),
        });
}

EigenSupport::VXd makeZeroDisplacement(const Mesh::TriMeshGeo& mesh) {
    return EigenSupport::VXd::Zero(mesh.numVertices() * 3);
}

EigenSupport::VXd makeSurfacePositions(const Mesh::TriMeshGeo& mesh) {
    EigenSupport::VXd positions(mesh.numVertices() * 3);
    for (int vi = 0; vi < mesh.numVertices(); ++vi) {
        positions.segment<3>(vi * 3) = mesh.pos(vi);
    }

    return positions;
}

std::array<int, 3> sortedTargetSampleIds(const std::array<int, 4>& pair) {
    std::array<int, 3> target{pair[1], pair[2], pair[3]};
    std::sort(target.begin(), target.end());
    return target;
}

TriangleMeshSelfContactHandler makeHandler(const Mesh::TriMeshGeo& mesh) {
    return TriangleMeshSelfContactHandler(mesh.positions(), mesh.triangles(), mesh.numVertices() * 3, 1);
}

}  // namespace

TEST(SelfContactHandlerTest, NearContactActivationUsesUnsignedDistance) {
    Logging::init();

    const Mesh::TriMeshGeo mesh = makeSeparatedParallelTriangleMesh(0.05);
    auto                   handler = makeHandler(mesh);
    const EigenSupport::VXd u = makeZeroDisplacement(mesh);

    handler.execute(u.data(), 0.1);

    EXPECT_GT(handler.getNumActivePairs(), 0);
    EXPECT_EQ(handler.getActiveVtxTriPairs().size(), static_cast<size_t>(handler.getNumActivePairs()));
    EXPECT_EQ(handler.getLastActiveClosestPoints().size(), handler.getNumActivePairs() * 3);
    EXPECT_EQ(handler.getLastActiveNormals().size(), handler.getNumActivePairs() * 3);

    const EigenSupport::VXd& distances = handler.getLastActiveDistances();
    ASSERT_EQ(distances.size(), handler.getNumActivePairs());
    for (int i = 0; i < distances.size(); ++i) {
        EXPECT_LT(distances[i], 0.1);
        EXPECT_NEAR(distances[i], 0.05, 1e-8);
    }
}

TEST(SelfContactHandlerTest, ActivationDistanceThresholdSuppressesPairs) {
    Logging::init();

    const Mesh::TriMeshGeo mesh = makeSeparatedParallelTriangleMesh(0.05);
    auto                   handler = makeHandler(mesh);
    const EigenSupport::VXd u = makeZeroDisplacement(mesh);

    handler.execute(u.data(), 0.01);

    EXPECT_EQ(handler.getNumActivePairs(), 0);
    EXPECT_EQ(handler.getLastActiveDistances().size(), 0);
}

TEST(SelfContactHandlerTest, ZeroBandFallsBackToCollisionOnlyWhenNotColliding) {
    Logging::init();

    const Mesh::TriMeshGeo mesh = makeSeparatedParallelTriangleMesh(0.05);
    auto                   handler = makeHandler(mesh);
    const EigenSupport::VXd u = makeZeroDisplacement(mesh);

    handler.execute(u.data(), 0.0);

    EXPECT_EQ(handler.getNumActivePairs(), 0);
    EXPECT_TRUE(handler.getCollidingTrianglePair().empty());
}

TEST(SelfContactHandlerTest, ZeroBandFallbackHandlesCollidingPairs) {
    Logging::init();

    const Mesh::TriMeshGeo mesh = makeIntersectingTriangleMesh();
    auto                   handler = makeHandler(mesh);
    const EigenSupport::VXd u = makeZeroDisplacement(mesh);

    handler.execute(u.data(), 0.0);

    EXPECT_GT(handler.getNumActivePairs(), 0);
    EXPECT_FALSE(handler.getCollidingTrianglePair().empty());
}

TEST(SelfContactHandlerTest, ClosestTargetSelectionKeepsOneTargetPerSample) {
    Logging::init();

    const Mesh::TriMeshGeo mesh = makeThreeLayerTriangleMesh(0.05, 0.08);
    auto                   handler = makeHandler(mesh);
    const EigenSupport::VXd u = makeZeroDisplacement(mesh);

    handler.execute(u.data(), 0.1);

    const auto& pairs = handler.getActiveVtxTriPairs();
    const auto& distances = handler.getLastActiveDistances();

    ASSERT_EQ(handler.getNumActivePairs(), 9);
    ASSERT_EQ(pairs.size(), 9u);
    ASSERT_EQ(distances.size(), 9);

    int middleLayerPairCount = 0;
    const std::array<int, 3> lowerTarget{0, 1, 2};
    for (size_t pairIdx = 0; pairIdx < pairs.size(); ++pairIdx) {
        if (pairs[pairIdx][0] < 3 || pairs[pairIdx][0] > 5) {
            continue;
        }

        middleLayerPairCount++;
        EXPECT_EQ(sortedTargetSampleIds(pairs[pairIdx]), lowerTarget);
        EXPECT_NEAR(distances[static_cast<int>(pairIdx)], 0.05, 1e-8);
    }

    EXPECT_EQ(middleLayerPairCount, 3);
}

TEST(SelfContactHandlerTest, LegacyPenaltyPathStillWorks) {
    Logging::init();

    const Mesh::TriMeshGeo mesh = makeIntersectingTriangleMesh();
    auto                   handler = makeHandler(mesh);
    const EigenSupport::VXd u = makeZeroDisplacement(mesh);

    handler.execute(u.data());
    EXPECT_FALSE(handler.getCollidingTrianglePair().empty());

    handler.handleContactDCD(0.0, 100);
    EXPECT_FALSE(handler.getCollidingVtxTriPairs().empty());
}

TEST(SelfContactHandlerTest, BarrierBuilderProducesFiniteEnergyForNearContactPairs) {
    Logging::init();

    const Mesh::TriMeshGeo mesh = makeSeparatedParallelTriangleMesh(0.05);
    auto                   handler = makeHandler(mesh);
    const EigenSupport::VXd u = makeZeroDisplacement(mesh);
    const EigenSupport::VXd x = makeSurfacePositions(mesh);

    handler.execute(u.data(), 0.1);
    ASSERT_GT(handler.getNumActivePairs(), 0);

    auto barrier = handler.buildBarrierEnergy(0.1);
    auto* buffer = barrier->allocateBuffer();
    barrier->setBuffer(buffer);
    barrier->setCoeff(5.0);
    barrier->setToPosFunction([](const EigenSupport::V3d& xLocal, EigenSupport::V3d& p, int) { p = xLocal; });

    const double value = barrier->func(x);
    EXPECT_TRUE(std::isfinite(value));
    EXPECT_GT(value, 0.0);

    barrier->freeBuffer(buffer);
}

TEST(SelfContactHandlerTest, AlphaUpperBoundShrinksForInwardDirections) {
    Logging::init();

    const Mesh::TriMeshGeo mesh = makeSeparatedParallelTriangleMesh(0.05);
    auto                   handler = makeHandler(mesh);
    const EigenSupport::VXd currentU = makeZeroDisplacement(mesh);

    handler.execute(currentU.data(), 0.1);
    ASSERT_GT(handler.getNumActivePairs(), 0);

    const auto&             pair = handler.getActiveVtxTriPairs().front();
    const EigenSupport::V3d n = handler.getLastActiveNormals().segment<3>(0);
    EigenSupport::VXd       du = EigenSupport::VXd::Zero(currentU.size());
    du.segment<3>(pair[0] * 3) = -0.1 * n;

    const double dSafe =
        PointTrianglePairBarrierEnergy::normalizePositiveZero(0.1, -1.0);
    const double alphaUpper = handler.computeEmbeddedAlphaUpperBound(currentU, du, 1.0, dSafe);

    EXPECT_GT(alphaUpper, 0.0);
    EXPECT_LT(alphaUpper, 1.0);
}

TEST(SelfContactHandlerTest, AlphaUpperBoundScalesWithAlphaSafety) {
    Logging::init();

    const Mesh::TriMeshGeo mesh = makeSeparatedParallelTriangleMesh(0.05);
    auto                   handler = makeHandler(mesh);
    const EigenSupport::VXd currentU = makeZeroDisplacement(mesh);

    handler.execute(currentU.data(), 0.1);
    ASSERT_GT(handler.getNumActivePairs(), 0);

    const auto&             pair = handler.getActiveVtxTriPairs().front();
    const EigenSupport::V3d n = handler.getLastActiveNormals().segment<3>(0);
    EigenSupport::VXd       du = EigenSupport::VXd::Zero(currentU.size());
    du.segment<3>(pair[0] * 3) = -0.1 * n;

    const double dSafe =
        PointTrianglePairBarrierEnergy::normalizePositiveZero(0.1, -1.0);
    const double alphaSafeOne  = handler.computeEmbeddedAlphaUpperBound(currentU, du, 1.0, dSafe);
    const double alphaSafeHalf = handler.computeEmbeddedAlphaUpperBound(currentU, du, 0.5, dSafe);

    EXPECT_GT(alphaSafeOne, 0.0);
    EXPECT_LT(alphaSafeOne, 1.0);
    EXPECT_NEAR(alphaSafeHalf, alphaSafeOne * 0.5, 1e-12);
}

TEST(SelfContactHandlerTest, AlphaUpperBoundAllowsOutwardRecovery) {
    Logging::init();

    const Mesh::TriMeshGeo mesh = makeSeparatedParallelTriangleMesh(0.05);
    auto                   handler = makeHandler(mesh);
    const EigenSupport::VXd currentU = makeZeroDisplacement(mesh);

    handler.execute(currentU.data(), 0.1);
    ASSERT_GT(handler.getNumActivePairs(), 0);

    const auto&             pair = handler.getActiveVtxTriPairs().front();
    const EigenSupport::V3d n = handler.getLastActiveNormals().segment<3>(0);
    EigenSupport::VXd       du = EigenSupport::VXd::Zero(currentU.size());
    du.segment<3>(pair[0] * 3) = 0.1 * n;

    const double dSafe =
        PointTrianglePairBarrierEnergy::normalizePositiveZero(0.1, -1.0);
    const double alphaUpper = handler.computeEmbeddedAlphaUpperBound(currentU, du, 0.9, dSafe);

    EXPECT_DOUBLE_EQ(alphaUpper, 1.0);
}

TEST(SelfContactHandlerTest, AlphaUpperBoundBlocksFurtherInwardMotionInsideSafeBand) {
    Logging::init();

    const double nearZeroGap = 5e-8;
    const Mesh::TriMeshGeo mesh = makeSeparatedParallelTriangleMesh(nearZeroGap);
    auto                   handler = makeHandler(mesh);
    const EigenSupport::VXd currentU = makeZeroDisplacement(mesh);

    handler.execute(currentU.data(), 0.1);
    ASSERT_GT(handler.getNumActivePairs(), 0);

    EigenSupport::VXd       inwardDu = EigenSupport::VXd::Zero(currentU.size());
    EigenSupport::VXd       outwardDu = EigenSupport::VXd::Zero(currentU.size());

    const auto& pairs = handler.getActiveVtxTriPairs();
    const auto& normals = handler.getLastActiveNormals();
    for (int pi = 0; pi < handler.getNumActivePairs(); ++pi) {
        const int sampleID = pairs[pi][0];
        const EigenSupport::V3d n = normals.segment<3>(pi * 3);
        inwardDu.segment<3>(sampleID * 3) = -0.01 * n;
        outwardDu.segment<3>(sampleID * 3) = 0.01 * n;
    }

    const double dSafe =
        PointTrianglePairBarrierEnergy::normalizePositiveZero(0.1, -1.0);

    const double inwardAlpha = handler.computeEmbeddedAlphaUpperBound(currentU, inwardDu, 1.0, dSafe);
    const double outwardAlpha = handler.computeEmbeddedAlphaUpperBound(currentU, outwardDu, 1.0, dSafe);

    EXPECT_DOUBLE_EQ(inwardAlpha, 0.0);
    EXPECT_DOUBLE_EQ(outwardAlpha, 1.0);
}

TEST(SelfContactHandlerTest, AlphaUpperBoundReturnsOneWithNoActivePairs) {
    Logging::init();

    const Mesh::TriMeshGeo mesh = makeSeparatedParallelTriangleMesh(0.05);
    auto                   handler = makeHandler(mesh);
    const EigenSupport::VXd currentU = makeZeroDisplacement(mesh);
    EigenSupport::VXd       du = EigenSupport::VXd::Zero(currentU.size());

    du.segment<3>(0) = EigenSupport::V3d(0.0, 0.0, -0.1);

    handler.execute(currentU.data(), 0.01);
    ASSERT_EQ(handler.getNumActivePairs(), 0);

    const double dSafe =
        PointTrianglePairBarrierEnergy::normalizePositiveZero(0.1, -1.0);
    EXPECT_DOUBLE_EQ(handler.computeEmbeddedAlphaUpperBound(currentU, du, 0.9, dSafe), 1.0);
}

}  // namespace pgo::Contact::test
