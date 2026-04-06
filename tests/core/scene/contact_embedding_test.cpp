#include "barycentricCoordinates.h"
#include "cubicMesh.h"
#include "tetMesh.h"
#include "triangleMeshExternalContactHandler.h"
#include "triangleMeshSelfContactHandler.h"
#include "triMeshGeo.h"
#include "pgoLogging.h"

#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <vector>

namespace pgo::Contact::test {

namespace {

std::filesystem::path repoRootPath() {
    return std::filesystem::absolute(std::filesystem::path(__FILE__))
        .parent_path()
        .parent_path()
        .parent_path()
        .parent_path();
}

EigenSupport::VXd makeSurfaceRestPositions(const Mesh::TriMeshGeo& surfaceMesh) {
    EigenSupport::VXd surfaceRestPositions(surfaceMesh.numVertices() * 3);
    for (int vi = 0; vi < surfaceMesh.numVertices(); ++vi) {
        surfaceRestPositions.segment<3>(vi * 3) = surfaceMesh.pos(vi);
    }

    return surfaceRestPositions;
}

double maxAbsDenseDifference(const EigenSupport::SpMatD& lhs, const EigenSupport::SpMatD& rhs) {
    EXPECT_EQ(lhs.rows(), rhs.rows());
    EXPECT_EQ(lhs.cols(), rhs.cols());

    const Eigen::MatrixXd lhsDense(lhs);
    const Eigen::MatrixXd rhsDense(rhs);
    return (lhsDense - rhsDense).cwiseAbs().maxCoeff();
}

template <class MeshType>
InterpolationCoordinates::BarycentricCoordinates makeBarycentricCoordinates(const Mesh::TriMeshGeo& surfaceMesh,
                                                                            const MeshType& mesh) {
    const EigenSupport::VXd surfaceRestPositions = makeSurfaceRestPositions(surfaceMesh);
    return {surfaceMesh.numVertices(), surfaceRestPositions.data(), &mesh};
}

template <class MeshType>
void expectExternalEmbeddingMatchesBarycentric(const std::filesystem::path& meshPath,
                                               const std::filesystem::path& surfacePath) {
    Logging::init();

    MeshType         mesh(meshPath.string().c_str());
    Mesh::TriMeshGeo surfaceMesh;
    ASSERT_TRUE(surfaceMesh.load(surfacePath.string()));

    auto bc = makeBarycentricCoordinates(surfaceMesh, mesh);
    const EigenSupport::SpMatD expected = bc.generateInterpolationMatrix();

    const std::vector<Mesh::TriMeshRef> externalSurfaces;
    Contact::TriangleMeshExternalContactHandler handler(surfaceMesh.positions(), surfaceMesh.triangles(),
                                                        mesh.getNumVertices() * 3, externalSurfaces, 1,
                                                        &bc.getEmbeddingVertexIndices(), &bc.getEmbeddingWeights());

    const EigenSupport::SpMatD& actual = handler.getSampleEmbeddingMatrix();
    ASSERT_EQ(actual.rows(), expected.rows());
    ASSERT_EQ(actual.cols(), expected.cols());
    EXPECT_NEAR(maxAbsDenseDifference(actual, expected), 0.0, 1e-12);
}

void expectSelfEmbeddingMatchesBarycentric(const std::filesystem::path& meshPath,
                                           const std::filesystem::path& surfacePath) {
    Logging::init();

    VolumetricMeshes::CubicMesh mesh(meshPath.string().c_str());
    Mesh::TriMeshGeo            surfaceMesh;
    ASSERT_TRUE(surfaceMesh.load(surfacePath.string()));

    auto bc = makeBarycentricCoordinates(surfaceMesh, mesh);
    const EigenSupport::SpMatD expected = bc.generateInterpolationMatrix();

    Contact::TriangleMeshSelfContactHandler handler(surfaceMesh.positions(), surfaceMesh.triangles(),
                                                    mesh.getNumVertices() * 3, 1, &bc.getEmbeddingVertexIndices(),
                                                    &bc.getEmbeddingWeights());

    const EigenSupport::SpMatD& actual = handler.getSampleEmbeddingMatrix();
    ASSERT_EQ(actual.rows(), expected.rows());
    ASSERT_EQ(actual.cols(), expected.cols());
    EXPECT_NEAR(maxAbsDenseDifference(actual, expected), 0.0, 1e-12);
}

Contact::TriangleMeshExternalContactHandler makePlanarExternalHandler() {
    const std::vector<EigenSupport::V3d> surfaceVertices = {
        EigenSupport::V3d(0.2, 0.2, 0.05),
        EigenSupport::V3d(0.2, 0.8, 0.05),
        EigenSupport::V3d(0.8, 0.2, 0.05),
    };
    const std::vector<EigenSupport::V3i> surfaceTriangles = {
        EigenSupport::V3i(0, 1, 2),
    };

    Mesh::TriMeshGeo externalMesh(
        std::vector<EigenSupport::V3d>{
            EigenSupport::V3d(0.0, 0.0, 0.0),
            EigenSupport::V3d(2.0, 0.0, 0.0),
            EigenSupport::V3d(0.0, 2.0, 0.0),
            EigenSupport::V3d(2.0, 2.0, 0.0),
        },
        std::vector<EigenSupport::V3i>{
            EigenSupport::V3i(0, 1, 2),
            EigenSupport::V3i(1, 3, 2),
        });
    const std::vector<Mesh::TriMeshRef> externalSurfaces = {externalMesh.ref()};

    return Contact::TriangleMeshExternalContactHandler(surfaceVertices, surfaceTriangles,
                                                       static_cast<int>(surfaceVertices.size()) * 3, externalSurfaces,
                                                       1);
}

}  // namespace

TEST(ContactEmbeddingTest, CubicExternalMatchesBarycentricEmbedding) {
    const std::filesystem::path repoRoot = repoRootPath();
    expectExternalEmbeddingMatchesBarycentric<VolumetricMeshes::CubicMesh>(
        repoRoot / "examples" / "cubic-box" / "cubic-box.veg",
        repoRoot / "examples" / "cubic-box" / "cubic-box.obj");
}

TEST(ContactEmbeddingTest, CubicSelfMatchesBarycentricEmbedding) {
    const std::filesystem::path repoRoot = repoRootPath();
    expectSelfEmbeddingMatchesBarycentric(repoRoot / "examples" / "cubic-box" / "cubic-box.veg",
                                          repoRoot / "examples" / "cubic-box" / "cubic-box.obj");
}

TEST(ContactEmbeddingTest, TetExternalMatchesBarycentricEmbedding) {
    const std::filesystem::path repoRoot = repoRootPath();
    expectExternalEmbeddingMatchesBarycentric<VolumetricMeshes::TetMesh>(
        repoRoot / "examples" / "box" / "box.veg", repoRoot / "examples" / "box" / "box.obj");
}

TEST(ContactEmbeddingTest, RejectsMalformedEmbeddingLayout) {
    Logging::init();

    const std::filesystem::path repoRoot = repoRootPath();
    VolumetricMeshes::CubicMesh mesh((repoRoot / "examples" / "cubic-box" / "cubic-box.veg").string().c_str());
    Mesh::TriMeshGeo            surfaceMesh;
    ASSERT_TRUE(surfaceMesh.load((repoRoot / "examples" / "cubic-box" / "cubic-box.obj").string()));

    auto bc = makeBarycentricCoordinates(surfaceMesh, mesh);

    std::vector<int>    badIndices = bc.getEmbeddingVertexIndices();
    std::vector<double> weights    = bc.getEmbeddingWeights();
    badIndices.pop_back();

    const std::vector<Mesh::TriMeshRef> externalSurfaces;
    EXPECT_THROW((void)Contact::TriangleMeshExternalContactHandler(surfaceMesh.positions(), surfaceMesh.triangles(),
                                                                   mesh.getNumVertices() * 3, externalSurfaces, 1,
                                                                   &badIndices, &weights),
                 std::runtime_error);
}

TEST(ContactEmbeddingTest, ExternalBarrierActivationIncludesNearContactSamples) {
    Logging::init();

    Contact::TriangleMeshExternalContactHandler handler = makePlanarExternalHandler();

    const std::vector<EigenSupport::V3d> surfaceVertices = {
        EigenSupport::V3d(0.2, 0.2, 0.05),
        EigenSupport::V3d(0.2, 0.8, 0.05),
        EigenSupport::V3d(0.8, 0.2, 0.05),
    };

    handler.execute(surfaceVertices);
    EXPECT_EQ(handler.getNumActiveSamples(), 0);
    EXPECT_EQ(handler.getNumCollidingSamples(), 0);

    handler.execute(surfaceVertices, 0.1);
    EXPECT_EQ(handler.getNumActiveSamples(), 3);
    EXPECT_EQ(handler.getNumCollidingSamples(), 3);

    const EigenSupport::VXd& signedDistances = handler.getLastSignedDistances();
    ASSERT_EQ(signedDistances.size(), 3);
    for (int i = 0; i < signedDistances.size(); ++i) {
        EXPECT_GT(signedDistances[i], 0.0);
        EXPECT_LT(signedDistances[i], 0.1);
        EXPECT_NEAR(signedDistances[i], 0.05, 1e-8);
    }
}

TEST(ContactEmbeddingTest, FeasibleStepUpperBoundMatchesAnalyticPlaneBound) {
    Logging::init();

    Contact::TriangleMeshExternalContactHandler handler = makePlanarExternalHandler();
    const EigenSupport::VXd currentU = EigenSupport::VXd::Zero(9);
    EigenSupport::VXd       du       = EigenSupport::VXd::Zero(9);
    const double            dSafe    = 0.01;
    const double            alphaSafety = 0.9;

    handler.execute(currentU.data(), 0.1);
    ASSERT_EQ(handler.getNumActiveSamples(), 3);

    for (int vi = 0; vi < 3; ++vi) {
        du[vi * 3 + 2] = -0.08;
    }

    const double alphaUpper = handler.computeSurfaceAlphaUpperBound(currentU, du, alphaSafety, dSafe);
    const double expected   = ((0.05 - dSafe) / 0.08) * alphaSafety;

    EXPECT_NEAR(alphaUpper, expected, 1e-12);
}

TEST(ContactEmbeddingTest, FeasibleStepUpperBoundScalesWithAlphaSafety) {
    Logging::init();

    Contact::TriangleMeshExternalContactHandler handler = makePlanarExternalHandler();
    const EigenSupport::VXd currentU = EigenSupport::VXd::Zero(9);
    EigenSupport::VXd       du       = EigenSupport::VXd::Zero(9);
    const double            dSafe    = 0.01;

    handler.execute(currentU.data(), 0.1);
    ASSERT_EQ(handler.getNumActiveSamples(), 3);

    for (int vi = 0; vi < 3; ++vi) {
        du[vi * 3 + 2] = -0.08;
    }

    const double alphaSafeOne  = handler.computeSurfaceAlphaUpperBound(currentU, du, 1.0, dSafe);
    const double alphaSafeHalf = handler.computeSurfaceAlphaUpperBound(currentU, du, 0.5, dSafe);

    EXPECT_GT(alphaSafeOne, 0.0);
    EXPECT_LT(alphaSafeOne, 1.0);
    EXPECT_NEAR(alphaSafeHalf, alphaSafeOne * 0.5, 1e-12);
}

TEST(ContactEmbeddingTest, FeasibleStepUpperBoundIsOneForMotionAwayFromContact) {
    Logging::init();

    Contact::TriangleMeshExternalContactHandler handler = makePlanarExternalHandler();
    const EigenSupport::VXd currentU = EigenSupport::VXd::Zero(9);
    EigenSupport::VXd       du       = EigenSupport::VXd::Zero(9);

    handler.execute(currentU.data(), 0.1);
    ASSERT_EQ(handler.getNumActiveSamples(), 3);

    for (int vi = 0; vi < 3; ++vi) {
        du[vi * 3 + 2] = 0.02;
    }

    EXPECT_DOUBLE_EQ(handler.computeSurfaceAlphaUpperBound(currentU, du, 0.9, 0.01), 1.0);
}

TEST(ContactEmbeddingTest, FeasibleStepUpperBoundReturnsOneWithNoActiveSamples) {
    Logging::init();

    Contact::TriangleMeshExternalContactHandler handler = makePlanarExternalHandler();
    const EigenSupport::VXd currentU = EigenSupport::VXd::Zero(9);
    EigenSupport::VXd       du       = EigenSupport::VXd::Zero(9);

    for (int vi = 0; vi < 3; ++vi) {
        du[vi * 3 + 2] = -0.04;
    }

    handler.execute(currentU.data(), 0.01);
    ASSERT_EQ(handler.getNumActiveSamples(), 0);

    EXPECT_DOUBLE_EQ(handler.computeSurfaceAlphaUpperBound(currentU, du, 0.9, 0.01), 1.0);
}

TEST(ContactEmbeddingTest, FeasibleStepUpperBoundReturnsZeroWhenAlreadyInsideSafeMargin) {
    Logging::init();

    Contact::TriangleMeshExternalContactHandler handler = makePlanarExternalHandler();
    EigenSupport::VXd       currentU = EigenSupport::VXd::Zero(9);
    const EigenSupport::VXd du       = EigenSupport::VXd::Zero(9);

    for (int vi = 0; vi < 3; ++vi) {
        currentU[vi * 3 + 2] = -0.045;
    }

    handler.execute(currentU.data(), 0.1);
    ASSERT_EQ(handler.getNumActiveSamples(), 3);

    EXPECT_DOUBLE_EQ(handler.computeSurfaceAlphaUpperBound(currentU, du, 0.9, 0.01), 0.0);
}

TEST(ContactEmbeddingTest, FeasibleStepUpperBoundAllowsRecoveryWhenInsideSafeMarginButMovingAway) {
    Logging::init();

    Contact::TriangleMeshExternalContactHandler handler = makePlanarExternalHandler();
    EigenSupport::VXd       currentU = EigenSupport::VXd::Zero(9);
    EigenSupport::VXd       du       = EigenSupport::VXd::Zero(9);

    for (int vi = 0; vi < 3; ++vi) {
        currentU[vi * 3 + 2] = -0.045;
        du[vi * 3 + 2]       = 0.02;
    }

    handler.execute(currentU.data(), 0.1);
    ASSERT_EQ(handler.getNumActiveSamples(), 3);

    EXPECT_DOUBLE_EQ(handler.computeSurfaceAlphaUpperBound(currentU, du, 0.9, 0.01), 1.0);
}

}  // namespace pgo::Contact::test
