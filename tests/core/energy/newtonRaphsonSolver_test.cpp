#include "NewtonRaphsonSolver.h"

#include <gtest/gtest.h>

#include <memory>
#include <vector>

namespace pgo::NonlinearOptimization::test {

namespace {

namespace ES = pgo::EigenSupport;

class QuadraticWellEnergy final : public PotentialEnergy {
public:
    double func(ES::ConstRefVecXd x) const override { return 0.5 * x.squaredNorm(); }

    void gradient(ES::ConstRefVecXd x, ES::RefVecXd grad) const override { grad = x; }

    void hessian(ES::ConstRefVecXd, ES::SpMatD& hess) const override { hess.coeffRef(0, 0) = 1.0; }

    void createHessian(ES::SpMatD& hess) const override {
        hess.resize(1, 1);
        hess.reserve(1);
        hess.insert(0, 0) = 1.0;
        hess.makeCompressed();
    }

    void getDOFs(std::vector<int>& dofs) const override { dofs = {0}; }

    int getNumDOFs() const override { return 1; }
};

class ConstantEnergyWithIdentityModel final : public PotentialEnergy {
public:
    double func(ES::ConstRefVecXd) const override { return 1.0; }

    void gradient(ES::ConstRefVecXd x, ES::RefVecXd grad) const override { grad = x; }

    void hessian(ES::ConstRefVecXd, ES::SpMatD& hess) const override { hess.coeffRef(0, 0) = 1.0; }

    void createHessian(ES::SpMatD& hess) const override {
        hess.resize(1, 1);
        hess.reserve(1);
        hess.insert(0, 0) = 1.0;
        hess.makeCompressed();
    }

    void getDOFs(std::vector<int>& dofs) const override { dofs = {0}; }

    int getNumDOFs() const override { return 1; }
};

NewtonRaphsonSolver makeDefaultSolver(const std::shared_ptr<const PotentialEnergy>& energy, double x0) {
    NewtonRaphsonSolver::SolverParam solverParam;
    ES::VXd                         x(1);
    x[0] = x0;
    return NewtonRaphsonSolver(x.data(), solverParam, energy, std::vector<int>(), nullptr);
}

}  // namespace

TEST(NewtonRaphsonSolverTest, ReturnsFeasibleStepZeroWhenAlphaUpperBoundIsZero) {
    auto energy = std::make_shared<QuadraticWellEnergy>();

    ES::VXd x(1);
    x[0] = 1.0;

    NewtonRaphsonSolver solver = makeDefaultSolver(energy, x[0]);
    solver.setAlphaTestFunc([](const ES::VXd&, const ES::VXd&) { return 0.0; });

    const int ret = solver.solve(x.data(), 10, 1e-12, 0);

    EXPECT_EQ(ret, NewtonRaphsonSolver::SR_FEASIBLE_STEP_ZERO);
    EXPECT_DOUBLE_EQ(x[0], 1.0);
}

TEST(NewtonRaphsonSolverTest, ReturnsStalledSmallStepWhenLineSearchFailsNearStationaryPoint) {
    auto energy = std::make_shared<ConstantEnergyWithIdentityModel>();

    ES::VXd x(1);
    x[0] = 1e-16;

    NewtonRaphsonSolver solver = makeDefaultSolver(energy, x[0]);

    const int ret = solver.solve(x.data(), 10, 1e-20, 0);

    EXPECT_EQ(ret, NewtonRaphsonSolver::SR_STALLED_SMALL_STEP);
}

TEST(NewtonRaphsonSolverTest, ReturnsMaxIterReachedWhenBudgetExpiresBeforeConvergenceCheck) {
    auto energy = std::make_shared<QuadraticWellEnergy>();

    ES::VXd x(1);
    x[0] = 1.0;

    NewtonRaphsonSolver solver = makeDefaultSolver(energy, x[0]);

    const int ret = solver.solve(x.data(), 1, 1e-20, 0);

    EXPECT_EQ(ret, NewtonRaphsonSolver::SR_MAX_ITER_REACHED);
    EXPECT_NEAR(x[0], 0.0, 1e-12);
}

}  // namespace pgo::NonlinearOptimization::test
