#include "triangleMeshSelfContactDetection.h"
#include "triMeshGeo.h"

#include <gtest/gtest.h>

#include <vector>

namespace pgo::Contact::test {

namespace {

std::vector<double> flattenPositions(const std::vector<EigenSupport::V3d>& positions) {
    std::vector<double> data(static_cast<size_t>(positions.size()) * 3);
    for (size_t i = 0; i < positions.size(); ++i) {
        data[i * 3 + 0] = positions[i][0];
        data[i * 3 + 1] = positions[i][1];
        data[i * 3 + 2] = positions[i][2];
    }

    return data;
}

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

Mesh::TriMeshGeo makeAdjacentTriangleMesh() {
    return Mesh::TriMeshGeo(
        std::vector<EigenSupport::V3d>{
            EigenSupport::V3d(0.0, 0.0, 0.0),
            EigenSupport::V3d(1.0, 0.0, 0.0),
            EigenSupport::V3d(0.0, 1.0, 0.0),
            EigenSupport::V3d(1.0, 1.0, 0.0),
        },
        std::vector<EigenSupport::V3i>{
            EigenSupport::V3i(0, 1, 2),
            EigenSupport::V3i(1, 3, 2),
        });
}

void expectSinglePair(const std::vector<std::pair<int, int>>& pairs, int triA = 0, int triB = 1) {
    ASSERT_EQ(pairs.size(), 1u);
    EXPECT_EQ(pairs[0].first, triA);
    EXPECT_EQ(pairs[0].second, triB);
}

}  // namespace

TEST(SelfContactDetectionTest, DcdNearContactRequiresActivationBand) {
    Mesh::TriMeshGeo mesh      = makeSeparatedParallelTriangleMesh(0.05);
    auto             positions = flattenPositions(mesh.positions());
    Contact::TriangleMeshSelfContactDetection detection(mesh.ref());

    detection.execute(positions.data(), static_cast<const double*>(nullptr));
    EXPECT_TRUE(detection.getCandidateTrianglePairs().empty());

    detection.execute(positions.data(), 0.0);
    EXPECT_TRUE(detection.getCandidateTrianglePairs().empty());

    detection.execute(positions.data(), 0.01);
    EXPECT_TRUE(detection.getCandidateTrianglePairs().empty());

    detection.execute(positions.data(), 0.1);
    expectSinglePair(detection.getCandidateTrianglePairs());
    expectSinglePair(detection.getPotentialCollidingTrianglePairs());
}

TEST(SelfContactDetectionTest, CollisionOnlyOverloadMatchesZeroBandOverload) {
    Mesh::TriMeshGeo mesh      = makeIntersectingTriangleMesh();
    auto             positions = flattenPositions(mesh.positions());
    Contact::TriangleMeshSelfContactDetection detection(mesh.ref());

    detection.execute(positions.data(), static_cast<const double*>(nullptr));
    const auto collisionOnlyPairs = detection.getCandidateTrianglePairs();
    expectSinglePair(collisionOnlyPairs);

    detection.execute(positions.data(), 0.0);
    const auto zeroBandPairs = detection.getCandidateTrianglePairs();
    ASSERT_EQ(collisionOnlyPairs.size(), zeroBandPairs.size());
    ASSERT_EQ(collisionOnlyPairs.size(), 1u);
    EXPECT_EQ(collisionOnlyPairs[0], zeroBandPairs[0]);
}

TEST(SelfContactDetectionTest, SharedVertexPairsRemainFiltered) {
    Mesh::TriMeshGeo mesh      = makeAdjacentTriangleMesh();
    auto             positions = flattenPositions(mesh.positions());
    Contact::TriangleMeshSelfContactDetection detection(mesh.ref());

    detection.execute(positions.data(), 0.1);
    EXPECT_TRUE(detection.getCandidateTrianglePairs().empty());
}

TEST(SelfContactDetectionTest, CcdNearContactUsesSweptBand) {
    Mesh::TriMeshGeo mesh = makeSeparatedParallelTriangleMesh(0.2);

    std::vector<EigenSupport::V3d> endPositions = mesh.positions();
    for (int vi = 3; vi < 6; ++vi) {
        endPositions[vi][2] = 0.05;
    }

    auto start = flattenPositions(mesh.positions());
    auto end   = flattenPositions(endPositions);

    Contact::TriangleMeshSelfContactDetection detection(mesh.ref());

    detection.execute(start.data(), end.data());
    EXPECT_TRUE(detection.getCandidateTrianglePairs().empty());

    detection.execute(start.data(), end.data(), 0.0);
    EXPECT_TRUE(detection.getCandidateTrianglePairs().empty());

    detection.execute(start.data(), end.data(), 0.1);
    expectSinglePair(detection.getCandidateTrianglePairs());
}

}  // namespace pgo::Contact::test
