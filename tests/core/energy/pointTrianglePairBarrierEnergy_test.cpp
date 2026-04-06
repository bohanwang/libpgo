#include "pointTrianglePairBarrierEnergy.h"

#include "barrierFunction.h"
#include "finiteDifference.h"

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <memory>
#include <vector>

namespace pgo::Contact::test {

namespace {

namespace ES = pgo::EigenSupport;
namespace BF = pgo::NonlinearOptimization::BarrierFunctions;

struct SinglePairBarrierCase {
    double            pointSampleWeight = 1.0;
    ES::V3d           normal            = ES::V3d::UnitZ();
    ES::V3d           baryWeights       = ES::V3d(1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0);
    std::array<int, 4> objectIDs{{0, 0, 0, 0}};
    std::array<int, 4> pointTrianglePairs{{0, 1, 2, 3}};
    std::array<int, 2> objectDOFOffsets{{0, 12}};
    ES::SpMatD         sampleEmbedding;
    std::vector<double> sampleWeights;
    std::shared_ptr<PointTrianglePairBarrierEnergy> energy;

    SinglePairBarrierCase(double pointSampleWeight_, const ES::V3d& normal_, const ES::V3d& baryWeights_, double dhat,
                          double positiveZero)
        : pointSampleWeight(pointSampleWeight_), normal(normal_), baryWeights(baryWeights_) {
        std::vector<ES::TripletD> entries;
        entries.reserve(12);
        for (int dofi = 0; dofi < 12; ++dofi) {
            entries.emplace_back(dofi, dofi, 1.0);
        }
        sampleEmbedding.resize(12, 12);
        sampleEmbedding.setFromTriplets(entries.begin(), entries.end());

        sampleWeights = {pointSampleWeight, 1.0, 1.0, 1.0};

        energy = std::make_shared<PointTrianglePairBarrierEnergy>(
            1, 1, &objectIDs, &pointTrianglePairs, objectDOFOffsets.data(), &sampleEmbedding, &sampleWeights,
            normal.data(), baryWeights.data(), dhat, positiveZero);
        energy->setToPosFunction([](const ES::V3d& x, ES::V3d& p, int) { p = x; });
    }
};

ES::V12d makeBarrierDirection(const ES::V3d& normal, const ES::V3d& baryWeights) {
    ES::V12d direction = ES::V12d::Zero();
    direction.segment<3>(0) = normal;
    direction.segment<3>(3) = -normal * baryWeights[0];
    direction.segment<3>(6) = -normal * baryWeights[1];
    direction.segment<3>(9) = -normal * baryWeights[2];
    return direction;
}

void expectFiniteSparseEntries(const ES::SpMatD& hess) {
    for (Eigen::Index i = 0; i < hess.nonZeros(); ++i) {
        EXPECT_TRUE(std::isfinite(hess.valuePtr()[i]));
    }
}

}  // namespace

TEST(PointTrianglePairBarrierEnergyTest, InactiveRegionReturnsZero) {
    SinglePairBarrierCase barrierCase(2.0, ES::V3d::UnitZ(), ES::V3d(0.2, 0.3, 0.5), 1.0, 1e-6);
    auto                   energy = barrierCase.energy;
    auto*                  buffer = energy->allocateBuffer();
    energy->setBuffer(buffer);
    energy->setCoeff(3.0);

    const ES::VXd x = (ES::VXd(12) << 0.2, -0.1, 1.25, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0).finished();

    EXPECT_DOUBLE_EQ(energy->func(x), 0.0);

    ES::VXd grad(12);
    energy->gradient(x, grad);
    EXPECT_LT(grad.norm(), 1e-12);

    ES::SpMatD hess;
    energy->createHessian(hess);
    energy->hessian(x, hess);
    EXPECT_EQ(hess.rows(), 12);
    EXPECT_EQ(hess.cols(), 12);
    EXPECT_GT(hess.nonZeros(), 0);
    EXPECT_NEAR(Eigen::MatrixXd(hess).norm(), 0.0, 1e-12);

    energy->freeBuffer(buffer);
}

TEST(PointTrianglePairBarrierEnergyTest, ActiveRegionMatchesFormulaAndFiniteDifference) {
    const double pointSampleWeight = 2.0;
    const double dhat              = 1.0;
    const double positiveZero      = 1e-6;
    const double ipcKappa          = 3.0;
    const ES::V3d normal           = ES::V3d::UnitZ();
    const ES::V3d baryWeights(0.2, 0.3, 0.5);

    SinglePairBarrierCase barrierCase(pointSampleWeight, normal, baryWeights, dhat, positiveZero);
    auto                   energy = barrierCase.energy;
    auto*                  buffer = energy->allocateBuffer();
    energy->setBuffer(buffer);
    energy->setCoeff(ipcKappa);

    const ES::VXd x =
        (ES::VXd(12) << 0.1, -0.2, 0.4, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0).finished();
    const double distance = 0.4;
    const double scale    = pointSampleWeight * ipcKappa;
    const ES::V12d direction = makeBarrierDirection(normal, baryWeights);

    EXPECT_NEAR(energy->func(x), scale * BF::logBarrierEnergy(distance, dhat, positiveZero), 1e-12);

    ES::VXd grad(12);
    energy->gradient(x, grad);
    const ES::V12d expectedGrad = direction * (scale * BF::logBarrierGradient(distance, dhat, positiveZero));
    EXPECT_NEAR((grad - expectedGrad).norm(), 0.0, 1e-10);

    ES::SpMatD hess;
    energy->createHessian(hess);
    energy->hessian(x, hess);
    const Eigen::MatrixXd denseHess(hess);
    const Eigen::Matrix<double, 12, 12> expectedHess =
        (direction * direction.transpose()) * (scale * BF::logBarrierHessian(distance, dhat, positiveZero));
    EXPECT_NEAR((denseHess - expectedHess).norm(), 0.0, 1e-10);

    NonlinearOptimization::FiniteDifference fd(NonlinearOptimization::FiniteDifference::M_FIVE_POINT, 1e-6);
    double gradRelError = 0.0;
    double hessRelError = 0.0;
    std::shared_ptr<const NonlinearOptimization::PotentialEnergy> potential = energy;
    fd.testEnergy(potential, true, true, -1.0, x.data(), -1, &gradRelError, &hessRelError);

    EXPECT_LT(gradRelError, 1e-6);
    EXPECT_LT(hessRelError, 1e-5);

    energy->freeBuffer(buffer);
}

TEST(PointTrianglePairBarrierEnergyTest, ClampAvoidsNaNForNonPositiveDistances) {
    SinglePairBarrierCase barrierCase(1.5, ES::V3d::UnitZ(), ES::V3d(0.2, 0.3, 0.5), 1.0, 1e-6);
    auto                   energy = barrierCase.energy;
    auto*                  buffer = energy->allocateBuffer();
    energy->setBuffer(buffer);
    energy->setCoeff(4.0);

    const ES::VXd x =
        (ES::VXd(12) << 0.0, 0.0, -0.2, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0).finished();

    const double value = energy->func(x);
    EXPECT_TRUE(std::isfinite(value));

    ES::VXd grad(12);
    energy->gradient(x, grad);
    for (int i = 0; i < grad.size(); ++i) {
        EXPECT_TRUE(std::isfinite(grad[i]));
    }

    ES::SpMatD hess;
    energy->createHessian(hess);
    energy->hessian(x, hess);
    expectFiniteSparseEntries(hess);

    energy->freeBuffer(buffer);
}

}  // namespace pgo::Contact::test
