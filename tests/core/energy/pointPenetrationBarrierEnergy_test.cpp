#include "pointPenetrationBarrierEnergy.h"

#include "barrierFunction.h"
#include "finiteDifference.h"

#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <vector>

namespace pgo::Contact::test {

namespace {

namespace ES = pgo::EigenSupport;
namespace BF = pgo::NonlinearOptimization::BarrierFunctions;

struct SinglePointBarrierCase {
    double constraintCoeff = 0.0;
    ES::V3d targetPosition = ES::V3d::Zero();
    ES::V3d normal         = ES::V3d::UnitZ();
    std::shared_ptr<PointPenetrationBarrierEnergy> energy;

    SinglePointBarrierCase(double constraintCoeff_, double dhat, double positiveZero)
        : constraintCoeff(constraintCoeff_) {
        static const std::vector<std::vector<int>> barycentricIdx = {{0}};
        static const std::vector<std::vector<double>> barycentricWeights = {{1.0}};

        energy = std::make_shared<PointPenetrationBarrierEnergy>(
            1, 3, &constraintCoeff, targetPosition.data(), normal.data(), barycentricIdx, barycentricWeights, dhat,
            positiveZero);
        energy->setComputePosFunction([](const ES::V3d& x, ES::V3d& p, int) { p = x; });
    }
};

void expectFiniteSparseEntries(const ES::SpMatD& hess) {
    for (Eigen::Index i = 0; i < hess.nonZeros(); ++i) {
        EXPECT_TRUE(std::isfinite(hess.valuePtr()[i]));
    }
}

}  // namespace

TEST(PointPenetrationBarrierEnergyTest, InactiveRegionReturnsZero) {
    SinglePointBarrierCase barrierCase(2.0, 1.0, 1e-6);
    auto                   energy = barrierCase.energy;
    auto*                  buffer = energy->allocateBuffer();
    energy->setBuffer(buffer);
    energy->setCoeff(3.0);

    const ES::VXd u = (ES::VXd(3) << 0.2, -0.1, 1.25).finished();

    EXPECT_DOUBLE_EQ(energy->func(u), 0.0);

    ES::VXd grad(3);
    energy->gradient(u, grad);
    EXPECT_LT(grad.norm(), 1e-12);

    ES::SpMatD hess;
    energy->createHessian(hess);
    energy->hessian(u, hess);
    EXPECT_EQ(hess.rows(), 3);
    EXPECT_EQ(hess.cols(), 3);
    EXPECT_GT(hess.nonZeros(), 0);
    EXPECT_NEAR(Eigen::MatrixXd(hess).norm(), 0.0, 1e-12);

    energy->freeBuffer(buffer);
}

TEST(PointPenetrationBarrierEnergyTest, ActiveRegionMatchesFormulaAndFiniteDifference) {
    const double constraintCoeff = 2.0;
    const double dhat            = 1.0;
    const double positiveZero    = 1e-6;
    const double ipcKappa        = 3.0;

    SinglePointBarrierCase barrierCase(constraintCoeff, dhat, positiveZero);
    auto                   energy = barrierCase.energy;
    auto* buffer = energy->allocateBuffer();
    energy->setBuffer(buffer);
    energy->setCoeff(ipcKappa);

    const ES::VXd u = (ES::VXd(3) << 0.1, -0.2, 0.4).finished();
    const double signedDistance = u[2];
    const double scale          = constraintCoeff * ipcKappa;

    EXPECT_NEAR(energy->func(u), scale * BF::logBarrierEnergy(signedDistance, dhat, positiveZero), 1e-12);

    ES::VXd grad(3);
    energy->gradient(u, grad);
    EXPECT_NEAR(grad[0], 0.0, 1e-12);
    EXPECT_NEAR(grad[1], 0.0, 1e-12);
    EXPECT_NEAR(grad[2], scale * BF::logBarrierGradient(signedDistance, dhat, positiveZero), 1e-10);

    ES::SpMatD hess;
    energy->createHessian(hess);
    energy->hessian(u, hess);
    const Eigen::MatrixXd denseHess(hess);
    const double tangentialBlockNorm = denseHess.topLeftCorner<2, 2>().norm();
    EXPECT_NEAR(denseHess(0, 0), 0.0, 1e-12);
    EXPECT_NEAR(denseHess(1, 1), 0.0, 1e-12);
    EXPECT_NEAR(denseHess(2, 2), scale * BF::logBarrierHessian(signedDistance, dhat, positiveZero), 1e-10);
    EXPECT_NEAR(tangentialBlockNorm, 0.0, 1e-12);
    EXPECT_NEAR(denseHess(0, 2), 0.0, 1e-12);
    EXPECT_NEAR(denseHess(2, 0), 0.0, 1e-12);
    EXPECT_NEAR(denseHess(1, 2), 0.0, 1e-12);
    EXPECT_NEAR(denseHess(2, 1), 0.0, 1e-12);

    NonlinearOptimization::FiniteDifference fd(NonlinearOptimization::FiniteDifference::M_FIVE_POINT, 1e-6);
    double gradRelError = 0.0;
    double hessRelError = 0.0;
    std::shared_ptr<const NonlinearOptimization::PotentialEnergy> potential = energy;
    fd.testEnergy(potential, true, true, -1.0, u.data(), -1, &gradRelError, &hessRelError);

    EXPECT_LT(gradRelError, 1e-6);
    EXPECT_LT(hessRelError, 1e-5);

    energy->freeBuffer(buffer);
}

TEST(PointPenetrationBarrierEnergyTest, ClampAvoidsNaNForNonPositiveDistances) {
    SinglePointBarrierCase barrierCase(1.5, 1.0, 1e-6);
    auto                   energy = barrierCase.energy;
    auto* buffer = energy->allocateBuffer();
    energy->setBuffer(buffer);
    energy->setCoeff(4.0);

    const ES::VXd u = (ES::VXd(3) << 0.0, 0.0, -0.2).finished();

    const double value = energy->func(u);
    EXPECT_TRUE(std::isfinite(value));

    ES::VXd grad(3);
    energy->gradient(u, grad);
    for (int i = 0; i < grad.size(); ++i) {
        EXPECT_TRUE(std::isfinite(grad[i]));
    }

    ES::SpMatD hess;
    energy->createHessian(hess);
    energy->hessian(u, hess);
    expectFiniteSparseEntries(hess);

    energy->freeBuffer(buffer);
}

}  // namespace pgo::Contact::test
