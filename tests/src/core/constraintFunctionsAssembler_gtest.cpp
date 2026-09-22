#include <gtest/gtest.h>

#include "constraintFunctionsAssember.h"
#include "pgoLogging.h"

#include <memory>
#include <vector>

namespace
{
namespace ES = pgo::EigenSupport;
namespace NO = pgo::NonlinearOptimization;

// c(x) = 0.5 (x_a^2 + x_b^2) - r on two global DOFs.
class SphereConstraint1D final : public NO::ConstraintFunction
{
public:
  SphereConstraint1D(int a, int b, double r): r_(r) { dofs = { a, b }; }
  double func(ES::ConstRefVecXd x) const override { return 0.5 * (x[dofs[0]] * x[dofs[0]] + x[dofs[1]] * x[dofs[1]]) - r_; }
  void jacobian(ES::ConstRefVecXd x, ES::SpMatD &jac) const override
  {
    jac.coeffRef(0, 0) = x[dofs[0]];
    jac.coeffRef(0, 1) = x[dofs[1]];
  }
  void hessian(ES::ConstRefVecXd, ES::SpMatD &hess) const override
  {
    hess.coeffRef(0, 0) = 1.0;
    hess.coeffRef(1, 1) = 1.0;
  }
  void createJacobian(ES::SpMatD &jac) const override
  {
    jac.resize(1, 2);
    jac.insert(0, 0) = 0.0;
    jac.insert(0, 1) = 0.0;
    jac.makeCompressed();
  }
  void createHessian(ES::SpMatD &hess) const override
  {
    hess.resize(2, 2);
    hess.insert(0, 0) = 0.0;
    hess.insert(1, 1) = 0.0;
    hess.makeCompressed();
  }

private:
  double r_;
};

// Two rows over all DOFs: c0 = 0.5 |x|^2 - 1 (curved), c1 = sum(x) - 2 (linear).
class SphereAndPlaneConstraintsND final : public NO::ConstraintFunctions
{
public:
  explicit SphereAndPlaneConstraintsND(int n): NO::ConstraintFunctions(n), n_(n)
  {
    std::vector<ES::TripletD> jac, hess;
    for (int i = 0; i < n; ++i) {
      jac.emplace_back(0, i, 1.0);
      jac.emplace_back(1, i, 1.0);
      hess.emplace_back(i, i, 1.0);
    }
    jacobianTemplate.resize(2, n);
    jacobianTemplate.setFromTriplets(jac.begin(), jac.end());
    lambdahTemplate.resize(n, n);
    lambdahTemplate.setFromTriplets(hess.begin(), hess.end());
  }
  void func(ES::ConstRefVecXd x, ES::RefVecXd g) const override
  {
    g[0] = 0.5 * x.squaredNorm() - 1.0;
    g[1] = x.sum() - 2.0;
  }
  void jacobian(ES::ConstRefVecXd x, ES::SpMatD &jac) const override
  {
    for (int i = 0; i < n_; ++i) {
      jac.coeffRef(0, i) = x[i];
      jac.coeffRef(1, i) = 1.0;
    }
  }
  void hessian(ES::ConstRefVecXd, ES::ConstRefVecXd lambda, ES::SpMatD &hess) const override
  {
    for (int i = 0; i < n_; ++i)
      hess.coeffRef(i, i) = lambda[0];
  }
  void hessianVector(ES::ConstRefVecXd, ES::ConstRefVecXd lambda, ES::ConstRefVecXd vec, ES::RefVecXd hessVec) const override
  {
    hessVec = lambda[0] * vec;
  }
  bool hasHessianVector() const override { return true; }

private:
  int n_;
};

std::shared_ptr<NO::ConstraintFunctionsAssembler> makeAssembler()
{
  auto assembler = std::make_shared<NO::ConstraintFunctionsAssembler>(4);
  assembler->addConstraint(std::make_shared<SphereConstraint1D>(0, 1, 1.0));
  assembler->addConstraint(std::make_shared<SphereConstraint1D>(2, 3, 1.0));
  assembler->addConstraint(std::make_shared<SphereAndPlaneConstraintsND>(4));
  assembler->init();
  return assembler;
}
}  // namespace

TEST(ConstraintFunctionsAssemblerGTest, HessianIsWrittenNotAccumulated)
{
  pgo::Logging::init();
  const auto assembler = makeAssembler();
  ASSERT_EQ(assembler->getNumConstraints(), 4);

  ES::VXd x(4), lambda(4);
  x << 0.3, -0.2, 0.7, 0.1;
  lambda << 2.0, 3.0, 5.0, 7.0;  // two 1-D multipliers, then the 2-row N-D block

  ES::SpMatD hess;
  assembler->createHessian(hess);
  assembler->hessian(x, lambda, hess);
  const ES::MXd first = ES::MXd(hess);
  assembler->hessian(x, lambda, hess);
  const ES::MXd second = ES::MXd(hess);

  // sphere_a on DOFs 0,1 weighted 2, sphere_b on DOFs 2,3 weighted 3, N-D row 0 on all DOFs weighted 5
  ES::MXd expected = ES::MXd::Zero(4, 4);
  expected.diagonal() << 2.0 + 5.0, 2.0 + 5.0, 3.0 + 5.0, 3.0 + 5.0;
  EXPECT_LT((first - expected).norm(), 1e-12);
  EXPECT_LT((second - expected).norm(), 1e-12);  // a second evaluation must not double the values
}

TEST(ConstraintFunctionsAssemblerGTest, HessianVectorVisitsEveryNDConstraintWithItsOwnMultipliers)
{
  pgo::Logging::init();
  const auto assembler = makeAssembler();
  ASSERT_TRUE(assembler->hasHessianVector());

  ES::VXd x(4), lambda(4), vec(4);
  x << 0.3, -0.2, 0.7, 0.1;
  lambda << 2.0, 3.0, 5.0, 7.0;
  vec << 1.0, -1.0, 0.5, 2.0;

  // Two 1-D constraints and one N-D constraint: the old loop bound would index
  // a second N-D constraint that does not exist and pass the unsliced lambda.
  ES::VXd hessVec = ES::VXd::Zero(4);
  assembler->hessianVector(x, lambda, vec, hessVec);

  ES::SpMatD hess;
  assembler->createHessian(hess);
  assembler->hessian(x, lambda, hess);
  const ES::VXd expected = hess * vec;
  EXPECT_LT((hessVec - expected).norm(), 1e-12);
  EXPECT_NEAR(hessVec[0], 7.0 * 1.0, 1e-12);
}
