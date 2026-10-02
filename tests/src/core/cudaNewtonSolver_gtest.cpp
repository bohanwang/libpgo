#include <gtest/gtest.h>

#include "cudaDeviceBuffer.h"
#include "cudaNewtonSolver.h"
#include "cudaQuadraticEnergy.h"

#include <cmath>
#include <limits>
#include <memory>

namespace
{
namespace ES = pgo::EigenSupport;
namespace NO = pgo::NonlinearOptimization;
using pgo::CudaUtilities::DeviceBuffer;
using pgo::CudaUtilities::checkCuda;

class CudaNewtonTest : public testing::Test
{
protected:
  void SetUp() override
  {
    int count = 0;
    if (cudaGetDeviceCount(&count) != cudaSuccess || count == 0)
      GTEST_SKIP() << "No CUDA device is accessible.";
    checkCuda(cudaSetDevice(0), "CUDA test device selection");
  }
};

ES::SpMatD scalarMatrix(double diagonal)
{
  ES::SpMatD A(1, 1);
  A.insert(0, 0) = diagonal;
  A.makeCompressed();
  return A;
}

class BoundedQuadratic : public NO::CudaQuadraticEnergy
{
public:
  explicit BoundedQuadratic(double maxStep, double constant = 0):
    NO::CudaQuadraticEnergy(scalarMatrix(1), ES::VXd::Constant(1, 10), ES::VXd::Zero(1), constant), bound(maxStep)
  {
  }
  double computeMaxStepSize(const double *, const double *direction) const override
  {
    double dx;
    checkCuda(cudaMemcpy(&dx, direction, sizeof(double), cudaMemcpyDeviceToHost), "test bound direction");
    return dx == 0 ? 1 : bound;
  }
  double bound;
};

class NonFiniteTrialEnergy : public BoundedQuadratic
{
public:
  NonFiniteTrialEnergy():
    BoundedQuadratic(1) {}
  double func(const double *x) override
  {
    double value = 0;
    checkCuda(cudaMemcpy(&value, x, sizeof(double), cudaMemcpyDeviceToHost), "test scalar read");
    return value == 0 ? BoundedQuadratic::func(x) : std::numeric_limits<double>::infinity();
  }
};

// Emulate one ULP of reduction error on an energy with a large constant. The
// exact quadratic decrease near x=10 is smaller than the represented energy's
// spacing, although the gradient still resolves the Newton correction.
class RoundedQuadratic : public BoundedQuadratic
{
public:
  RoundedQuadratic():
    BoundedQuadratic(1) {}
  double func(const double *x) override
  {
    double coordinate = 0;
    checkCuda(cudaMemcpy(&coordinate, x, sizeof(double), cudaMemcpyDeviceToHost), "test scalar read");
    const double value = BoundedQuadratic::func(x) + 1e12;
    return coordinate > 10 - 5e-6 ? std::nextafter(value, std::numeric_limits<double>::infinity()) : value;
  }
  EnergyChange energyChange(const double *base, const double *trial) override
  {
    // This test deliberately represents an energy without a stable local
    // change implementation, so exercise the interface's scalar fallback.
    return NO::CudaPotentialEnergy::energyChange(base, trial);
  }
};

class RoundedEnergyWithoutResidualImprovement : public RoundedQuadratic
{
public:
  void addGradient(const double *, double *out) override
  {
    double gradient = 0;
    checkCuda(cudaMemcpy(&gradient, out, sizeof(double), cudaMemcpyDeviceToHost), "test gradient read");
    gradient -= 1e-5;
    checkCuda(cudaMemcpy(out, &gradient, sizeof(double), cudaMemcpyHostToDevice), "test gradient write");
  }
};

class ProjectedScalarQuadratic : public BoundedQuadratic
{
public:
  enum class ExactMode
  {
    Valid,
    NegativeCurvature,
    NonFiniteDiagonal,
    NonFiniteProduct
  };
  explicit ProjectedScalarQuadratic(ExactMode mode = ExactMode::Valid):
    BoundedQuadratic(1), mode(mode),
    projected(scalarMatrix(100), ES::VXd::Zero(1), ES::VXd::Zero(1)),
    indefinite(scalarMatrix(-1), ES::VXd::Zero(1), ES::VXd::Zero(1))
  {
  }
  void prepareHessian(const double *) override { exact = false; }
  void prepareExactHessian(const double *) override { exact = true; }
  void addHessianVector(const double *v, double *out) const override
  {
    if (!exact)
      projected.addHessianVector(v, out);
    else if (mode == ExactMode::NegativeCurvature)
      indefinite.addHessianVector(v, out);
    else if (mode == ExactMode::NonFiniteProduct)
      writeNonFinite(out);
    else
      BoundedQuadratic::addHessianVector(v, out);
  }
  void addDiagonal(double *out) const override
  {
    if (!exact)
      projected.addDiagonal(out);
    else if (mode == ExactMode::NegativeCurvature)
      indefinite.addDiagonal(out);
    else if (mode == ExactMode::NonFiniteDiagonal)
      writeNonFinite(out);
    else
      BoundedQuadratic::addDiagonal(out);
  }

private:
  static void writeNonFinite(double *out)
  {
    const double value = std::numeric_limits<double>::quiet_NaN();
    checkCuda(cudaMemcpy(out, &value, sizeof(value), cudaMemcpyHostToDevice), "test invalid exact operator");
  }
  ExactMode mode;
  bool exact = false;
  NO::CudaQuadraticEnergy projected, indefinite;
};

class SmallDisplacementQuadratic : public NO::CudaQuadraticEnergy
{
public:
  SmallDisplacementQuadratic():
    NO::CudaQuadraticEnergy(scalarMatrix(1e15), ES::VXd::Constant(1, 1e-15), ES::VXd::Zero(1)) {}
  double computeMaxStepSize(const double *, const double *direction) const override
  {
    double dx;
    checkCuda(cudaMemcpy(&dx, direction, sizeof(double), cudaMemcpyDeviceToHost), "test small direction");
    return dx == 0 ? 1 : 0.5;
  }
};

// Poison only the unused queued products after one-step convergence. The
// actual operator and final true-residual product are the identity. This
// exercises sticky guards without relying on zero times NaN being harmless.
class PoisonedQueuedProducts : public BoundedQuadratic
{
public:
  explicit PoisonedQueuedProducts(int batch):
    BoundedQuadratic(1), batch(batch) {}
  void addHessianVector(const double *v, double *out) const override
  {
    ++calls;
    if (calls > 1 && calls <= batch) {
      const double value = std::numeric_limits<double>::quiet_NaN();
      checkCuda(cudaMemcpy(out, &value, sizeof(value), cudaMemcpyHostToDevice), "test unused product poison");
    }
    else
      BoundedQuadratic::addHessianVector(v, out);
  }
  mutable int calls = 0;
  int batch;
};

class InvalidProductQuadratic : public BoundedQuadratic
{
public:
  explicit InvalidProductQuadratic(double value):
    BoundedQuadratic(1), value(value) {}
  void addHessianVector(const double *, double *out) const override
  {
    checkCuda(cudaMemcpy(out, &value, sizeof(value), cudaMemcpyHostToDevice), "test invalid product");
  }
  double value;
};

class ConservativeApproachQuadratic : public BoundedQuadratic
{
public:
  ConservativeApproachQuadratic():
    BoundedQuadratic(1) {}
  double computeMaxStepSize(const double *, const double *dx) const override
  {
    double direction;
    checkCuda(cudaMemcpy(&direction, dx, sizeof(double), cudaMemcpyDeviceToHost), "test conservative advancement");
    return std::abs(direction) > 1 ? 0.5 : 1;
  }
};

// An SPD surrogate bends Newton's path around a circular obstacle even though
// the minimizer is hidden behind it from the input. Every individual Newton
// segment is feasible; their endpoint chord need not be.
class CurvedPathQuadratic : public NO::CudaQuadraticEnergy
{
public:
  explicit CurvedPathQuadratic(double scale = 1):
    NO::CudaQuadraticEnergy((ES::MXd::Identity(2, 2) / (scale * scale)).sparseView(),
      (ES::VXd(2) << scale, 0).finished(), ES::VXd::Zero(2)),
    surrogate(((ES::MXd(2, 2) << 2, -1, -1, 1).finished() / (scale * scale)).sparseView(), ES::VXd::Zero(2), ES::VXd::Zero(2)),
    scale(scale) {}
  void addHessianVector(const double *v, double *out) const override { surrogate.addHessianVector(v, out); }
  void addDiagonal(double *out) const override { surrogate.addDiagonal(out); }
  double func(const double *x) override
  {
    ES::V2d value;
    checkCuda(cudaMemcpy(value.data(), x, 2 * sizeof(double), cudaMemcpyDeviceToHost), "test trial read");
    evaluatedHiddenPoint |= bound(ES::V2d::Zero(), value / scale) < 1;
    return NO::CudaQuadraticEnergy::func(x);
  }
  double computeMaxStepSize(const double *x, const double *dx) const override
  {
    ES::V2d value, direction;
    checkCuda(cudaMemcpy(value.data(), x, 2 * sizeof(double), cudaMemcpyDeviceToHost), "test input read");
    checkCuda(cudaMemcpy(direction.data(), dx, 2 * sizeof(double), cudaMemcpyDeviceToHost), "test direction read");
    return bound(value / scale, direction / scale);
  }
  static double bound(const ES::V2d &x, const ES::V2d &dx)
  {
    // Use extended precision for the independent geometric oracle, including
    // near-tangent chords reached by the visibility-constrained solve.
    const long double rx = static_cast<long double>(x[0]) - 0.5L, ry = x[1];
    const long double vx = dx[0], vy = dx[1];
    const long double radius = 0.2L;
    const long double a = vx * vx + vy * vy, b = 2 * (rx * vx + ry * vy), c = rx * rx + ry * ry - radius * radius;
    if (c <= 0)
      return 0;
    if (a == 0)
      return 1;
    const long double discriminant = b * b - 4 * a * c;
    if (discriminant < 0)
      return 1;
    const long double collisionTime = b < 0 ? 2 * c / (-b + std::sqrt(discriminant)) :
                                              (-b - std::sqrt(discriminant)) / (2 * a);
    return collisionTime >= 0 && collisionTime <= 1 ? static_cast<double>(0.9L * collisionTime) : 1;
  }
  bool evaluatedHiddenPoint = false;

private:
  NO::CudaQuadraticEnergy surrogate;
  double scale;
};

class QuantizedResidualQuadratic : public NO::CudaQuadraticEnergy
{
public:
  explicit QuantizedResidualQuadratic(bool blockCorrection = false):
    NO::CudaQuadraticEnergy((ES::M2d() << 1e13, 1e9, 1e9, 1e7).finished().sparseView(),
      (ES::V2d() << -0.016, 0).finished(), (ES::V2d() << 8e-6, -8e-8).finished(), 1e8),
    blockCorrection(blockCorrection) {}
  EnergyChange energyChange(const double *base, const double *trial) override
  {
    return NO::CudaPotentialEnergy::energyChange(base, trial);
  }
  double computeMaxStepSize(const double *, const double *direction) const override
  {
    ES::V2d d;
    checkCuda(cudaMemcpy(d.data(), direction, sizeof(double) * 2, cudaMemcpyDeviceToHost), "test coordinate bound");
    return blockCorrection && d[1] < 0 ? 0 : 1;
  }

private:
  bool blockCorrection;
};

// The energy is finite and stationary at an invalid input. Feasibility must
// be checked before either convergence or the accepted-iterate observer.
class InvalidInitialQuadratic : public BoundedQuadratic
{
public:
  explicit InvalidInitialQuadratic(double initialBound):
    BoundedQuadratic(1), initialBound(initialBound) {}
  double func(const double *x) override
  {
    ++energyEvaluations;
    return BoundedQuadratic::func(x);
  }
  double computeMaxStepSize(const double *, const double *direction) const override
  {
    double dx;
    checkCuda(cudaMemcpy(&dx, direction, sizeof(double), cudaMemcpyDeviceToHost), "test initial direction");
    EXPECT_DOUBLE_EQ(dx, 0);
    ++initialChecks;
    return initialBound;
  }
  double initialBound;
  int energyEvaluations = 0;
  mutable int initialChecks = 0;
};

// A clipped real-valued direction can round across a boundary between two
// doubles. Rechecking the actual endpoint must reject that rounded trial.
class RoundedBoundaryQuadratic : public NO::CudaQuadraticEnergy
{
public:
  RoundedBoundaryQuadratic():
    NO::CudaQuadraticEnergy(scalarMatrix(1e16),
      ES::VXd::Constant(1, 1 + 4 * std::numeric_limits<double>::epsilon()), ES::VXd::Zero(1)) {}
  static bool feasible(double x)
  {
    return static_cast<long double>(x) < 1 + 1.75L * std::numeric_limits<double>::epsilon();
  }
  double computeMaxStepSize(const double *x, const double *dx) const override
  {
    double value, direction;
    checkCuda(cudaMemcpy(&value, x, sizeof(double), cudaMemcpyDeviceToHost), "test rounded boundary input");
    checkCuda(cudaMemcpy(&direction, dx, sizeof(double), cudaMemcpyDeviceToHost), "test rounded boundary direction");
    if (!feasible(value))
      return 0;
    const long double gap = 1 + 1.75L * std::numeric_limits<double>::epsilon() - value;
    return direction > 0 && direction >= gap ? static_cast<double>(0.9L * gap / direction) : 1;
  }
};

// Near a visibility boundary, rounding can leave the blocked coordinate fixed
// while another coordinate advances a few ULPs. Emulate relative reduction
// variation of 2e-13 in the dominant residual: the objective is unchanged and
// the measured norm remains far above the requested convergence tolerance.
class VaryingResidualBoundaryQuadratic : public NO::CudaQuadraticEnergy
{
public:
  VaryingResidualBoundaryQuadratic():
    NO::CudaQuadraticEnergy(ES::M2d::Identity().sparseView(),
      (ES::V2d() << 2, 0.5).finished(), ES::V2d::Zero()) {}
  void addGradient(const double *x, double *out) override
  {
    NO::CudaQuadraticEnergy::addGradient(x, out);
    ES::V2d value;
    checkCuda(cudaMemcpy(value.data(), out, 2 * sizeof(double), cudaMemcpyDeviceToHost), "test varying residual read");
    value[0] *= 1 + (++evaluations % 2 ? 2e-13 : -2e-13);
    checkCuda(cudaMemcpy(out, value.data(), 2 * sizeof(double), cudaMemcpyHostToDevice), "test varying residual write");
  }
  double computeMaxStepSize(const double *x, const double *dx) const override
  {
    ES::V2d value, direction;
    checkCuda(cudaMemcpy(value.data(), x, 2 * sizeof(double), cudaMemcpyDeviceToHost), "test varying boundary input");
    checkCuda(cudaMemcpy(direction.data(), dx, 2 * sizeof(double), cudaMemcpyDeviceToHost), "test varying boundary direction");
    if (!RoundedBoundaryQuadratic::feasible(value[0]))
      return 0;
    const long double gap = 1 + 1.75L * std::numeric_limits<double>::epsilon() - value[0];
    return direction[0] > 0 && direction[0] >= gap ? static_cast<double>(0.9L * gap / direction[0]) : 1;
  }

private:
  int evaluations = 0;
};

TEST_F(CudaNewtonTest, SparseQuadraticMatchesDenseEnergyAndAdditiveOperators)
{
  ES::MXd dense(4, 4);
  dense << 5, 1, 0, -0.5,
    1, 4, 0.5, 0,
    0, 0.5, 3, 1,
    -0.5, 0, 1, 2;
  const ES::SpMatD A = dense.sparseView();
  const ES::VXd reference = (ES::VXd(4) << 1, -2, 3, 4).finished();
  const ES::VXd linear = (ES::VXd(4) << -2, 3, 1, -4).finished();
  const ES::VXd x = (ES::VXd(4) << 0.5, -1, 4, 3).finished();
  const ES::VXd v = (ES::VXd(4) << 0.2, 0.5, -0.3, 1).finished();
  NO::CudaQuadraticEnergy energy(A, reference, linear, 7);
  DeviceBuffer<double> deviceX, deviceV, out;
  deviceX.upload(x.data(), 4);
  deviceV.upload(v.data(), 4);
  const ES::VXd du = x - reference;
  EXPECT_NEAR(energy.func(deviceX.data()), 0.5 * du.dot(dense * du) + linear.dot(du) + 7, 1e-12);

  ES::VXd actual = ES::VXd::Constant(4, 2);
  out.upload(actual.data(), 4);
  energy.addGradient(deviceX.data(), out.data());
  out.download(actual.data(), 4);
  EXPECT_LT((actual - (dense * du + linear + ES::VXd::Constant(4, 2))).norm(), 1e-12);

  actual.setConstant(-3);
  out.upload(actual.data(), 4);
  energy.prepareHessian(deviceX.data());
  energy.addHessianVector(deviceV.data(), out.data());
  out.download(actual.data(), 4);
  EXPECT_LT((actual - (dense * v - ES::VXd::Constant(4, 3))).norm(), 1e-12);

  actual.setOnes();
  out.upload(actual.data(), 4);
  energy.addDiagonal(out.data());
  out.download(actual.data(), 4);
  EXPECT_LT((actual - dense.diagonal() - ES::VXd::Ones(4)).norm(), 1e-12);
}

TEST_F(CudaNewtonTest, ComposedSparseQuadraticConvergesInOneNewtonStep)
{
  ES::MXd dense(3, 3);
  dense << 4, 1, 0, 1, 3, 0.5, 0, 0.5, 2;
  const ES::SpMatD A = dense.sparseView();
  const ES::VXd reference = (ES::VXd(3) << 3, -4, 1).finished();
  const ES::VXd linear = (ES::VXd(3) << 1, 2, -3).finished();
  ES::SpMatD identity(3, 3);
  identity.setIdentity();
  auto first = std::make_shared<NO::CudaQuadraticEnergy>(A, reference, linear);
  auto second = std::make_shared<NO::CudaQuadraticEnergy>(identity, ES::VXd::Zero(3), ES::VXd::Zero(3));
  NO::CudaNewtonSolver::SolverParam parameters;
  parameters.pcgRelativeTolerance = 1e-12;
  NO::CudaNewtonSolver solver({ first, second }, parameters);
  ES::VXd x = ES::VXd::Zero(3);
  const ES::VXd expected = (dense + ES::MXd::Identity(3, 3)).ldlt().solve(dense * reference - linear);
  EXPECT_EQ(solver.solve(x.data(), 1, 1e-10), NO::CudaNewtonSolver::SOLVE_CONVERGED);
  EXPECT_LT((x - expected).norm(), 1e-10);
  EXPECT_EQ(solver.getDiagnostics().newtonIterations, 1);
  EXPECT_GT(solver.getDiagnostics().pcgIterations, 0);
  EXPECT_LT(solver.getDiagnostics().finalGradientInfinityNorm, 1e-10);
}

TEST_F(CudaNewtonTest, FeasibleBoundIsNeverIncreased)
{
  auto bounded = std::make_shared<BoundedQuadratic>(0.25);
  NO::CudaNewtonSolver solver({ bounded });
  double x = 0;
  EXPECT_EQ(solver.solve(&x, 1, 1e-10), NO::CudaNewtonSolver::SOLVE_NOT_CONVERGED);
  EXPECT_NEAR(x, 2.5, 1e-12);
  EXPECT_DOUBLE_EQ(solver.getDiagnostics().lastFeasibleAlpha, 0.25);
  EXPECT_DOUBLE_EQ(solver.getDiagnostics().lastAcceptedAlpha, 0.25);

  bounded->bound = 0;
  EXPECT_EQ(solver.solve(&x, 10, 1e-10), NO::CudaNewtonSolver::SOLVE_NOT_CONVERGED);
  EXPECT_NEAR(x, 2.5, 1e-12);
  EXPECT_EQ(solver.getDiagnostics().stopReason, "zero_feasible_step");
}

TEST_F(CudaNewtonTest, CurvatureFailureRegularizesAndCanBeDisabled)
{
  auto concave = std::make_shared<NO::CudaQuadraticEnergy>(scalarMatrix(-1), ES::VXd::Zero(1), ES::VXd::Zero(1));
  NO::CudaNewtonSolver::SolverParam parameters;
  parameters.maxRegularizationAttempts = 0;
  NO::CudaNewtonSolver unregularized({ concave }, parameters);
  double x = 1;
  EXPECT_EQ(unregularized.solve(&x, 1, 1e-10), NO::CudaNewtonSolver::SOLVE_NOT_CONVERGED);
  EXPECT_DOUBLE_EQ(x, 1);
  EXPECT_EQ(unregularized.getDiagnostics().stopReason, "linear_curvature_or_descent_failure");

  NO::CudaNewtonSolver regularized({ concave });
  EXPECT_EQ(regularized.solve(&x, 1, 1e-10), NO::CudaNewtonSolver::SOLVE_NOT_CONVERGED);
  EXPECT_GT(x, 1);
  EXPECT_TRUE(std::isfinite(x));
  EXPECT_GT(regularized.getDiagnostics().regularizationAttempts, 0);
  EXPECT_GT(regularized.getDiagnostics().maxRegularization, 1);
}

TEST_F(CudaNewtonTest, FailedLineSearchKeepsLastFiniteState)
{
  auto wall = std::make_shared<NonFiniteTrialEnergy>();
  NO::CudaNewtonSolver solver({ wall });
  double x = 0;
  EXPECT_EQ(solver.solve(&x, 4, 1e-10), NO::CudaNewtonSolver::SOLVE_NOT_CONVERGED);
  EXPECT_DOUBLE_EQ(x, 0);
  EXPECT_EQ(solver.getDiagnostics().stopReason, "line_search_failed");
  EXPECT_TRUE(std::isfinite(solver.getDiagnostics().finalEnergy));
}

TEST_F(CudaNewtonTest, ZeroIterationStatusDescribesReturnedState)
{
  auto energy = std::make_shared<BoundedQuadratic>(1);
  NO::CudaNewtonSolver solver({ energy });
  int observedInputs = 0;
  solver.setAcceptedIterateObserver([&](int iteration, const double *, double, double) {
    EXPECT_EQ(iteration, 0);
    ++observedInputs;
  });
  double x = 0;
  EXPECT_EQ(solver.solve(&x, 0, 1e-10), NO::CudaNewtonSolver::SOLVE_NOT_CONVERGED);
  EXPECT_DOUBLE_EQ(x, 0);
  x = 10;
  EXPECT_EQ(solver.solve(&x, 0, 1e-10), NO::CudaNewtonSolver::SOLVE_CONVERGED);
  EXPECT_DOUBLE_EQ(x, 10);
  EXPECT_EQ(observedInputs, 2);
  x = std::numeric_limits<double>::quiet_NaN();
  EXPECT_EQ(solver.solve(&x, 10, 1e-10), NO::CudaNewtonSolver::SOLVE_NUMERICAL_FAILURE);
  EXPECT_TRUE(std::isnan(x));
  EXPECT_EQ(observedInputs, 2);
}

TEST_F(CudaNewtonTest, StrictNewtonToleranceIsNotBlockedByPCGAbsoluteTolerance)
{
  auto energy = std::make_shared<NO::CudaQuadraticEnergy>(scalarMatrix(1),
    ES::VXd::Zero(1), ES::VXd::Constant(1, 1e-13));
  NO::CudaNewtonSolver solver({ energy });
  double x = 0;
  EXPECT_EQ(solver.solve(&x, 1, 1e-16), NO::CudaNewtonSolver::SOLVE_CONVERGED);
  EXPECT_NEAR(x, -1e-13, 1e-28);
  EXPECT_LT(solver.getDiagnostics().finalGradientInfinityNorm, 1e-16);
}

TEST_F(CudaNewtonTest, RoundoffLimitedEnergyAcceptsOnlyMeasuredResidualImprovement)
{
  auto energy = std::make_shared<RoundedQuadratic>();
  NO::CudaNewtonSolver solver({ energy });
  double x = 10 - 1e-5;
  EXPECT_EQ(solver.solve(&x, 1, 1e-10), NO::CudaNewtonSolver::SOLVE_CONVERGED);
  EXPECT_DOUBLE_EQ(x, 10);
  EXPECT_EQ(solver.getDiagnostics().roundoffAcceptedSteps, 1);
  EXPECT_LT(solver.getDiagnostics().finalGradientInfinityNorm, 1e-10);

  energy->bound = 0.25;
  x = 10 - 1e-5;
  EXPECT_EQ(solver.solve(&x, 1, 1e-10), NO::CudaNewtonSolver::SOLVE_NOT_CONVERGED);
  EXPECT_NEAR(x, 10 - 0.75e-5, 2e-15);
  EXPECT_DOUBLE_EQ(solver.getDiagnostics().lastAcceptedAlpha, 0.25);
  EXPECT_EQ(solver.getDiagnostics().roundoffAcceptedSteps, 1);
  EXPECT_GT(solver.getDiagnostics().finalGradientInfinityNorm, 1e-10);
  energy->bound = 0;
  const double unchanged = x;
  EXPECT_EQ(solver.solve(&x, 1, 1e-10), NO::CudaNewtonSolver::SOLVE_NOT_CONVERGED);
  EXPECT_DOUBLE_EQ(x, unchanged);
  EXPECT_EQ(solver.getDiagnostics().stopReason, "zero_feasible_step");
  EXPECT_EQ(solver.getDiagnostics().roundoffAcceptedSteps, 0);

  auto stalledEnergy = std::make_shared<RoundedEnergyWithoutResidualImprovement>();
  NO::CudaNewtonSolver stalledSolver({ stalledEnergy });
  x = 10 - 1e-5;
  EXPECT_EQ(stalledSolver.solve(&x, 1, 1e-10), NO::CudaNewtonSolver::SOLVE_NOT_CONVERGED);
  EXPECT_DOUBLE_EQ(x, 10 - 1e-5);
  EXPECT_EQ(stalledSolver.getDiagnostics().roundoffAcceptedSteps, 0);
  EXPECT_EQ(stalledSolver.getDiagnostics().stopReason, "line_search_failed");
  EXPECT_GT(stalledSolver.getDiagnostics().finalGradientInfinityNorm, 1e-10);
}

TEST_F(CudaNewtonTest, StableEnergyChangeAllowsSmallDescentBelowAbsoluteEnergySpacing)
{
  // A conservative bound allows only 0.1% residual improvement. The exact
  // decrease is meaningful, even though the constant makes both absolute
  // energies identical and a roundoff residual fallback would reject it.
  auto energy = std::make_shared<BoundedQuadratic>(0.001, 1e12);
  NO::CudaNewtonSolver solver({ energy });
  double x = 10 - 1e-5;
  const double before = x;
  EXPECT_EQ(solver.solve(&x, 1, 1e-10), NO::CudaNewtonSolver::SOLVE_NOT_CONVERGED);
  EXPECT_GT(x, before);
  EXPECT_NEAR(x, before + 0.001 * (10 - before), 2e-15);
  EXPECT_EQ(solver.getDiagnostics().newtonIterations, 1);
  EXPECT_EQ(solver.getDiagnostics().roundoffAcceptedSteps, 0);
  EXPECT_DOUBLE_EQ(solver.getDiagnostics().lastAcceptedAlpha, 0.001);
  EXPECT_GT(solver.getDiagnostics().finalGradientInfinityNorm, 0.99 * (10 - before));
  EXPECT_LT(solver.getDiagnostics().finalGradientInfinityNorm, 10 - before);
}

TEST_F(CudaNewtonTest, ExactTangentAcceleratesSlowProjectedConvergence)
{
  auto energy = std::make_shared<ProjectedScalarQuadratic>();
  NO::CudaNewtonSolver::SolverParam params;
  params.exactHessianStartIteration = -1;
  NO::CudaNewtonSolver projected({ energy }, params);
  double x = 0;
  EXPECT_EQ(projected.solve(&x, 3, 1e-10), NO::CudaNewtonSolver::SOLVE_NOT_CONVERGED);
  EXPECT_EQ(projected.getDiagnostics().exactHessianAttempts, 0);
  EXPECT_GT(projected.getDiagnostics().finalGradientInfinityNorm, 9);

  params.exactHessianStartIteration = 2;
  NO::CudaNewtonSolver exact({ energy }, params);
  x = 0;
  EXPECT_EQ(exact.solve(&x, 3, 1e-10), NO::CudaNewtonSolver::SOLVE_CONVERGED);
  EXPECT_DOUBLE_EQ(x, 10);
  EXPECT_EQ(exact.getDiagnostics().newtonIterations, 3);
  EXPECT_EQ(exact.getDiagnostics().exactHessianAttempts, 1);
  EXPECT_EQ(exact.getDiagnostics().exactHessianDirections, 1);
  EXPECT_EQ(exact.getDiagnostics().exactHessianFallbacks, 0);
}

TEST_F(CudaNewtonTest, UnusableExactTangentRestoresProjectedOperator)
{
  for (auto mode : { ProjectedScalarQuadratic::ExactMode::NegativeCurvature,
         ProjectedScalarQuadratic::ExactMode::NonFiniteDiagonal,
         ProjectedScalarQuadratic::ExactMode::NonFiniteProduct }) {
    auto energy = std::make_shared<ProjectedScalarQuadratic>(mode);
    NO::CudaNewtonSolver::SolverParam params;
    params.exactHessianStartIteration = 0;
    params.maxRegularizationAttempts = 0;
    NO::CudaNewtonSolver solver({ energy }, params);
    double x = 0;
    EXPECT_EQ(solver.solve(&x, 1, 1e-10), NO::CudaNewtonSolver::SOLVE_NOT_CONVERGED);
    EXPECT_NEAR(x, 0.1, 1e-15);
    EXPECT_EQ(solver.getDiagnostics().exactHessianAttempts, 1);
    EXPECT_EQ(solver.getDiagnostics().exactHessianDirections, 0);
    EXPECT_EQ(solver.getDiagnostics().exactHessianFallbacks, 1);
    EXPECT_EQ(solver.getDiagnostics().regularizationAttempts, 0);
    EXPECT_TRUE(std::isfinite(solver.getDiagnostics().finalGradientInfinityNorm));
  }
}

TEST_F(CudaNewtonTest, RepresentableSmallDisplacementsRetainStrictConvergence)
{
  auto energy = std::make_shared<SmallDisplacementQuadratic>();
  NO::CudaNewtonSolver solver({ energy });
  double x = 0;
  EXPECT_EQ(solver.solve(&x, 40, 1e-9), NO::CudaNewtonSolver::SOLVE_CONVERGED);
  EXPECT_NEAR(x, 1e-15, 1e-24);
  EXPECT_GT(solver.getDiagnostics().newtonIterations, 1);
  EXPECT_LT(solver.getDiagnostics().finalGradientInfinityNorm, 1e-9);
  EXPECT_DOUBLE_EQ(solver.getDiagnostics().lastAcceptedAlpha, 0.5);
}

TEST_F(CudaNewtonTest, InputChordCheckDistinguishesPhysicalAndNewtonPaths)
{
  NO::CudaNewtonSolver::SolverParam params;
  params.exactHessianStartIteration = -1;
  auto unguardedEnergy = std::make_shared<CurvedPathQuadratic>();
  NO::CudaNewtonSolver unguarded({ unguardedEnergy }, params);
  ES::V2d x = ES::V2d::Zero();
  EXPECT_EQ(unguarded.solve(x.data(), 100, 1e-8), NO::CudaNewtonSolver::SOLVE_CONVERGED);
  EXPECT_LT(CurvedPathQuadratic::bound(ES::V2d::Zero(), x), 1);
  EXPECT_TRUE(unguardedEnergy->evaluatedHiddenPoint);

  auto guardedEnergy = std::make_shared<CurvedPathQuadratic>();
  params.checkInputChord = true;
  NO::CudaNewtonSolver guarded({ guardedEnergy }, params);
  x.setZero();
  EXPECT_EQ(guarded.solve(x.data(), 100, 1e-8), NO::CudaNewtonSolver::SOLVE_NOT_CONVERGED);
  EXPECT_DOUBLE_EQ(CurvedPathQuadratic::bound(ES::V2d::Zero(), x), 1) << x.transpose();
  EXPECT_FALSE(guardedEnergy->evaluatedHiddenPoint);
  EXPECT_GT(guarded.getDiagnostics().inputChordRejections, 0);
  EXPECT_GT(guarded.getDiagnostics().finalGradientInfinityNorm, 1e-8);
  EXPECT_TRUE(x.allFinite());
}

TEST_F(CudaNewtonTest, InputChordCertificationRejectsIncompleteCertificates)
{
  auto energy = std::make_shared<ConservativeApproachQuadratic>();
  NO::CudaNewtonSolver::SolverParam params;
  params.checkInputChord = true;
  NO::CudaNewtonSolver solver({ energy }, params);
  double x = 0;
  EXPECT_EQ(solver.solve(&x, 10, 1e-10), NO::CudaNewtonSolver::SOLVE_NOT_CONVERGED);
  EXPECT_GT(x, 0);
  EXPECT_LE(x, 1);
  EXPECT_GT(solver.getDiagnostics().inputChordRejections, 0);
  EXPECT_GT(solver.getDiagnostics().newtonIterations, 1);
}

TEST_F(CudaNewtonTest, ResidualCorrectionRecoversRepresentableQuantizedRoot)
{
  // The continuous Newton correction in x[0] is below half a displacement
  // ULP. Rounding that component makes its ordinary direction increase the
  // residual. The smaller-magnitude x[1] coordinate can still resolve it.
  NO::CudaNewtonSolver::SolverParam disabled;
  disabled.maxResidualCorrectionCandidates = 0;
  auto energy = std::make_shared<QuantizedResidualQuadratic>();
  NO::CudaNewtonSolver ordinary({ energy }, disabled);
  ES::V2d x(-0.016, 0);
  EXPECT_EQ(ordinary.solve(x.data(), 2, 1e-6), NO::CudaNewtonSolver::SOLVE_NOT_CONVERGED);
  EXPECT_GE(ordinary.getDiagnostics().finalGradientInfinityNorm, 1e-6);

  NO::CudaNewtonSolver::SolverParam enabled;
  enabled.checkInputChord = true;
  enabled.checkTrialChord = true;
  NO::CudaNewtonSolver recovered({ energy }, enabled);
  std::vector<ES::V2d> observed;
  recovered.setAcceptedIterateObserver([&](int iteration, const double *deviceX, double, double) {
    EXPECT_EQ(iteration, static_cast<int>(observed.size()));
    ES::V2d value;
    checkCuda(cudaMemcpy(value.data(), deviceX, sizeof(double) * 2, cudaMemcpyDeviceToHost), "test recovery observer");
    observed.push_back(value);
  });
  x << -0.016, 0;
  EXPECT_EQ(recovered.solve(x.data(), 2, 1e-6), NO::CudaNewtonSolver::SOLVE_CONVERGED);
  EXPECT_DOUBLE_EQ(x[0], -0.016);
  EXPECT_LT(x[1], 0);
  EXPECT_LT(recovered.getDiagnostics().finalGradientInfinityNorm, 1e-6);
  EXPECT_EQ(recovered.getDiagnostics().residualCorrectionAttempts, 1);
  EXPECT_EQ(recovered.getDiagnostics().residualCorrectionSteps, 1);
  EXPECT_LE(recovered.getDiagnostics().residualCorrectionCandidates, 2);
  ASSERT_EQ(observed.size(), recovered.getDiagnostics().newtonIterations + 1);
  EXPECT_EQ(observed.back(), x);
}

TEST_F(CudaNewtonTest, ResidualCorrectionRetainsFeasibilityCertificate)
{
  NO::CudaNewtonSolver::SolverParam params;
  params.checkInputChord = true;
  params.checkTrialChord = true;
  NO::CudaNewtonSolver solver({ std::make_shared<QuantizedResidualQuadratic>(true) }, params);
  ES::V2d x(-0.016, 0);
  EXPECT_EQ(solver.solve(x.data(), 2, 1e-6), NO::CudaNewtonSolver::SOLVE_NOT_CONVERGED);
  EXPECT_DOUBLE_EQ(x[0], -0.016);
  EXPECT_GE(x[1], 0);
  EXPECT_GE(solver.getDiagnostics().finalGradientInfinityNorm, 1e-6);
  EXPECT_EQ(solver.getDiagnostics().residualCorrectionAttempts, 1);
  EXPECT_EQ(solver.getDiagnostics().residualCorrectionSteps, 0);
  EXPECT_GT(solver.getDiagnostics().trialChordRejections, 0);
}

TEST_F(CudaNewtonTest, ObserverRecordsCertifiedCurvedOptimizationPath)
{
  auto energy = std::make_shared<CurvedPathQuadratic>();
  NO::CudaNewtonSolver::SolverParam params;
  params.checkTrialChord = true;
  params.exactHessianStartIteration = -1;
  NO::CudaNewtonSolver solver({ energy }, params);
  std::vector<ES::V2d> path;
  solver.setAcceptedIterateObserver([&](int iteration, const double *deviceX, double energyValue, double gradientInf) {
    EXPECT_EQ(iteration, static_cast<int>(path.size()));
    ES::V2d value;
    checkCuda(cudaMemcpy(value.data(), deviceX, 2 * sizeof(double), cudaMemcpyDeviceToHost), "test optimization path observer");
    const ES::V2d gradient = value - ES::V2d(1, 0);
    EXPECT_NEAR(energyValue, 0.5 * gradient.squaredNorm(), 1e-12);
    EXPECT_NEAR(gradientInf, gradient.cwiseAbs().maxCoeff(), 1e-12);
    if (!path.empty())
      EXPECT_DOUBLE_EQ(CurvedPathQuadratic::bound(path.back(), value - path.back()), 1);
    path.push_back(value);
  });
  ES::V2d x = ES::V2d::Zero();
  EXPECT_EQ(solver.solve(x.data(), 100, 1e-8), NO::CudaNewtonSolver::SOLVE_CONVERGED);
  ASSERT_EQ(path.size(), solver.getDiagnostics().newtonIterations + 1);
  EXPECT_EQ(path.front(), ES::V2d::Zero());
  EXPECT_EQ(path.back(), x);
  EXPECT_LT(CurvedPathQuadratic::bound(path.front(), path.back() - path.front()), 1);
  EXPECT_EQ(solver.getDiagnostics().inputChordRejections, 0);
}

TEST_F(CudaNewtonTest, TrialChordCheckRejectsRoundedBoundaryCrossing)
{
  auto energy = std::make_shared<RoundedBoundaryQuadratic>();
  NO::CudaNewtonSolver ordinary({ energy });
  double x = 1;
  EXPECT_EQ(ordinary.solve(&x, 1, 1e-8), NO::CudaNewtonSolver::SOLVE_NOT_CONVERGED);
  EXPECT_FALSE(RoundedBoundaryQuadratic::feasible(x));

  NO::CudaNewtonSolver::SolverParam params;
  params.checkTrialChord = true;
  NO::CudaNewtonSolver guarded({ energy }, params);
  int observed = 0;
  guarded.setAcceptedIterateObserver([&](int, const double *deviceX, double, double) {
    double value;
    checkCuda(cudaMemcpy(&value, deviceX, sizeof(double), cudaMemcpyDeviceToHost), "test rounded trial observer");
    EXPECT_TRUE(RoundedBoundaryQuadratic::feasible(value));
    ++observed;
  });
  x = 1;
  EXPECT_EQ(guarded.solve(&x, 1, 1e-8), NO::CudaNewtonSolver::SOLVE_NOT_CONVERGED);
  EXPECT_DOUBLE_EQ(x, std::nextafter(1.0, 2.0));
  EXPECT_GT(guarded.getDiagnostics().trialChordRejections, 0);
  EXPECT_EQ(guarded.getDiagnostics().inputChordRejections, 0);
  EXPECT_EQ(observed, guarded.getDiagnostics().newtonIterations + 1);
}

TEST_F(CudaNewtonTest, ObserverFailurePropagatesInsteadOfOmittingAcceptedState)
{
  NO::CudaNewtonSolver solver({ std::make_shared<BoundedQuadratic>(1) });
  int observations = 0;
  solver.setAcceptedIterateObserver([&](int iteration, const double *, double, double) {
    ++observations;
    if (iteration == 1)
      throw std::runtime_error("recording failed");
  });
  double x = 0;
  EXPECT_THROW(solver.solve(&x, 1, 1e-10), std::runtime_error);
  EXPECT_EQ(observations, 2);
}

TEST_F(CudaNewtonTest, BatchedPCGConvergesInsideBatchWithoutExtraUpdates)
{
  for (int size : { 1, 3, 11 }) {
    ES::MXd matrix = ES::MXd::Zero(size, size);
    for (int i = 0; i < size; ++i) {
      matrix(i, i) = 2 + 0.1 * i;
      if (i)
        matrix(i - 1, i) = matrix(i, i - 1) = 0.25;
    }
    ES::VXd reference = ES::VXd::LinSpaced(size, 1, 2);
    int expectedIterations = -1;
    for (int batch : { 0, 1, 8, 16 }) {
      SCOPED_TRACE("size=" + std::to_string(size) + " batch=" + std::to_string(batch));
      auto energy = std::make_shared<NO::CudaQuadraticEnergy>(matrix.sparseView(), reference, ES::VXd::Zero(size));
      NO::CudaNewtonSolver::SolverParam params;
      params.pcgReductionBatchSize = batch;
      params.maxPCGIterations = 64;
      params.pcgRelativeTolerance = 1e-13;
      NO::CudaNewtonSolver solver({ energy }, params);
      ES::VXd x = ES::VXd::Zero(size);
      EXPECT_EQ(solver.solve(x.data(), 1, 1e-10), NO::CudaNewtonSolver::SOLVE_CONVERGED);
      EXPECT_LT((x - reference).norm(), 1e-10);
      EXPECT_EQ(solver.getDiagnostics().newtonIterations, 1);
      if (!batch)
        expectedIterations = solver.getDiagnostics().pcgIterations;
      EXPECT_EQ(solver.getDiagnostics().pcgIterations, expectedIterations);
      if (size == 1)
        EXPECT_EQ(expectedIterations, 1);
    }
  }
}

TEST_F(CudaNewtonTest, BatchedPCGFreezesBeforePoisonedQueuedScratch)
{
  for (int batch : { 8, 16 }) {
    auto energy = std::make_shared<PoisonedQueuedProducts>(batch);
    NO::CudaNewtonSolver::SolverParam params;
    params.pcgReductionBatchSize = batch;
    params.maxPCGIterations = batch;
    NO::CudaNewtonSolver solver({ energy }, params);
    double x = 0;
    EXPECT_EQ(solver.solve(&x, 1, 1e-10), NO::CudaNewtonSolver::SOLVE_CONVERGED);
    EXPECT_DOUBLE_EQ(x, 10);
    EXPECT_EQ(solver.getDiagnostics().pcgIterations, 1);
    EXPECT_DOUBLE_EQ(solver.getDiagnostics().lastLinearResidualNorm, 0);
    EXPECT_EQ(energy->calls, batch + 1);  // Includes the final true-residual product.
  }
}

TEST_F(CudaNewtonTest, BatchedPCGHonorsPartialBatchIterationCaps)
{
  constexpr int size = 16;
  ES::MXd matrix = ES::MXd::Zero(size, size);
  for (int i = 0; i < size; ++i) {
    matrix(i, i) = 2 + 0.1 * i;
    if (i)
      matrix(i - 1, i) = matrix(i, i - 1) = 0.5;
  }
  const ES::VXd reference = ES::VXd::LinSpaced(size, 1, 3);
  for (int cap : { 1, 3, 9 }) {
    ES::VXd expected;
    for (int batch : { 0, 8, 16 }) {
      SCOPED_TRACE("cap=" + std::to_string(cap) + " batch=" + std::to_string(batch));
      auto energy = std::make_shared<NO::CudaQuadraticEnergy>(matrix.sparseView(), reference, ES::VXd::Zero(size));
      NO::CudaNewtonSolver::SolverParam params;
      params.pcgReductionBatchSize = batch;
      params.maxPCGIterations = cap;
      params.pcgRelativeTolerance = 1e-15;
      params.pcgAbsoluteTolerance = 0;
      params.maxRegularizationAttempts = 0;
      NO::CudaNewtonSolver solver({ energy }, params);
      ES::VXd x = ES::VXd::Zero(size);
      EXPECT_EQ(solver.solve(x.data(), 1, 1e-14), NO::CudaNewtonSolver::SOLVE_NOT_CONVERGED);
      EXPECT_EQ(solver.getDiagnostics().pcgIterations, cap);
      EXPECT_EQ(solver.getDiagnostics().newtonIterations, 1);
      EXPECT_DOUBLE_EQ(solver.getDiagnostics().lastAcceptedAlpha, 1);
      EXPECT_NEAR(solver.getDiagnostics().lastLinearResidualNorm, (matrix * (x - reference)).norm(), 1e-13);
      if (!batch)
        expected = x;
      EXPECT_LT((x - expected).norm(), 1e-12);
    }
  }
}

TEST_F(CudaNewtonTest, BatchedPCGCurvatureFailuresAndRetriesResetState)
{
  ES::MXd indefinite(2, 2);
  indefinite << 1, 2, 2, 1;
  for (int batch : { 0, 8, 16 }) {
    NO::CudaNewtonSolver::SolverParam params;
    params.pcgReductionBatchSize = batch;
    params.maxRegularizationAttempts = 0;
    auto energy = std::make_shared<NO::CudaQuadraticEnergy>(indefinite.sparseView(), ES::V2d::Zero(), ES::V2d(1, 0));
    NO::CudaNewtonSolver solver({ energy }, params);
    ES::V2d x = ES::V2d::Zero();
    EXPECT_EQ(solver.solve(x.data(), 1, 1e-10), NO::CudaNewtonSolver::SOLVE_NOT_CONVERGED);
    EXPECT_TRUE(x.isZero());
    EXPECT_EQ(solver.getDiagnostics().pcgIterations, 1);  // The second curvature is negative.
    EXPECT_EQ(solver.getDiagnostics().stopReason, "linear_curvature_or_descent_failure");

    auto concave = std::make_shared<NO::CudaQuadraticEnergy>(scalarMatrix(-1), ES::VXd::Zero(1), ES::VXd::Zero(1));
    NO::CudaNewtonSolver firstFailure({ concave }, params);
    double scalar = 1;
    EXPECT_EQ(firstFailure.solve(&scalar, 1, 1e-10), NO::CudaNewtonSolver::SOLVE_NOT_CONVERGED);
    EXPECT_DOUBLE_EQ(scalar, 1);
    EXPECT_EQ(firstFailure.getDiagnostics().pcgIterations, 0);
    params.maxRegularizationAttempts = 12;
    NO::CudaNewtonSolver retried({ concave }, params);
    EXPECT_EQ(retried.solve(&scalar, 1, 1e-10), NO::CudaNewtonSolver::SOLVE_NOT_CONVERGED);
    EXPECT_GT(scalar, 1);
    EXPECT_TRUE(std::isfinite(scalar));
    EXPECT_GT(retried.getDiagnostics().regularizationAttempts, 0);
    EXPECT_EQ(retried.getDiagnostics().pcgIterations, 1);

    for (auto mode : { ProjectedScalarQuadratic::ExactMode::NegativeCurvature, ProjectedScalarQuadratic::ExactMode::NonFiniteProduct }) {
      auto fallback = std::make_shared<ProjectedScalarQuadratic>(mode);
      params.exactHessianStartIteration = 0;
      params.maxRegularizationAttempts = 0;
      NO::CudaNewtonSolver restored({ fallback }, params);
      scalar = 0;
      EXPECT_EQ(restored.solve(&scalar, 1, 1e-10), NO::CudaNewtonSolver::SOLVE_NOT_CONVERGED);
      EXPECT_NEAR(scalar, 0.1, 1e-15);
      EXPECT_EQ(restored.getDiagnostics().exactHessianFallbacks, 1);
      EXPECT_EQ(restored.getDiagnostics().pcgIterations, 1);
    }
  }
}

TEST_F(CudaNewtonTest, BatchedPCGRejectsNonfiniteProductsAndOverflowedSolution)
{
  for (int batch : { 0, 8, 16 }) {
    NO::CudaNewtonSolver::SolverParam params;
    params.pcgReductionBatchSize = batch;
    params.maxRegularizationAttempts = 0;
    for (double bad : { std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity() }) {
      auto energy = std::make_shared<InvalidProductQuadratic>(bad);
      NO::CudaNewtonSolver solver({ energy }, params);
      double x = 0;
      EXPECT_EQ(solver.solve(&x, 1, 1e-10), NO::CudaNewtonSolver::SOLVE_NUMERICAL_FAILURE);
      EXPECT_DOUBLE_EQ(x, 0);
      EXPECT_EQ(solver.getDiagnostics().pcgIterations, 0);
      EXPECT_EQ(solver.getDiagnostics().newtonIterations, 0);
      EXPECT_EQ(solver.getDiagnostics().stopReason, "non_finite_linear_solve");
    }
    // Finite alpha and recurrence convergence can still overflow dx. The
    // unchanged final true-residual/solution check must reject that direction.
    auto overflow = std::make_shared<NO::CudaQuadraticEnergy>(scalarMatrix(5e-309), ES::VXd::Zero(1), ES::VXd::Ones(1));
    NO::CudaNewtonSolver solver({ overflow }, params);
    double x = 0;
    EXPECT_EQ(solver.solve(&x, 1, 1e-10), NO::CudaNewtonSolver::SOLVE_NUMERICAL_FAILURE);
    EXPECT_DOUBLE_EQ(x, 0);
    EXPECT_EQ(solver.getDiagnostics().pcgIterations, 1);
    EXPECT_EQ(solver.getDiagnostics().newtonIterations, 0);
    EXPECT_EQ(solver.getDiagnostics().stopReason, "non_finite_linear_solve");
  }
}

TEST_F(CudaNewtonTest, BatchedPCGReportsTrueResidualOnIllConditionedSystem)
{
  constexpr int size = 8;
  ES::MXd hilbert(size, size);
  for (int i = 0; i < size; ++i)
    for (int j = 0; j < size; ++j)
      hilbert(i, j) = 1.0 / (i + j + 1);
  const ES::VXd reference = ES::VXd::LinSpaced(size, 1, 2);
  for (int batch : { 0, 8, 16 }) {
    auto energy = std::make_shared<NO::CudaQuadraticEnergy>(hilbert.sparseView(), reference, ES::VXd::Zero(size));
    NO::CudaNewtonSolver::SolverParam params;
    params.pcgReductionBatchSize = batch;
    params.pcgRelativeTolerance = 1e-24;
    params.pcgAbsoluteTolerance = 0;
    params.maxPCGIterations = 128;
    params.maxRegularizationAttempts = 0;
    NO::CudaNewtonSolver solver({ energy }, params);
    ES::VXd x = ES::VXd::Zero(size), initialGradient(size), product(size);
    DeviceBuffer<double> input, output(size);
    input.upload(x.data(), size);
    output.zero();
    energy->addGradient(input.data(), output.data());
    output.download(initialGradient.data(), size);
    EXPECT_EQ(solver.solve(x.data(), 1, 1e-24), NO::CudaNewtonSolver::SOLVE_NOT_CONVERGED);
    ASSERT_EQ(solver.getDiagnostics().newtonIterations, 1);
    ASSERT_DOUBLE_EQ(solver.getDiagnostics().lastAcceptedAlpha, 1);
    input.upload(x.data(), size);
    output.zero();
    energy->addHessianVector(input.data(), output.data());
    output.download(product.data(), size);
    const double trueNorm = (initialGradient + product).norm();
    EXPECT_GT(trueNorm, 1e-24 * initialGradient.norm());
    EXPECT_NEAR(solver.getDiagnostics().lastLinearResidualNorm, trueNorm, 1e-12 * trueNorm);
  }
}

TEST_F(CudaNewtonTest, BatchedPCGRejectsInvalidBatchSize)
{
  auto energy = std::make_shared<BoundedQuadratic>(1);
  for (int size : { -1, 65 }) {
    NO::CudaNewtonSolver::SolverParam params;
    params.pcgReductionBatchSize = size;
    EXPECT_THROW(NO::CudaNewtonSolver({ energy }, params), std::invalid_argument);
  }
}

TEST_F(CudaNewtonTest, StagnationWindowExitsBlockedPlateauWithoutClaimingConvergence)
{
  for (double scale : { 1e-12, 1e12 }) {
    SCOPED_TRACE("length scale=" + std::to_string(scale));
    auto energy = std::make_shared<CurvedPathQuadratic>(scale);
    NO::CudaNewtonSolver::SolverParam params;
    params.checkInputChord = true;
    params.exactHessianStartIteration = -1;
    params.stagnationWindow = 4;
    NO::CudaNewtonSolver solver({ energy }, params);
    int observations = 0;
    ES::V2d lastObserved = ES::V2d::Zero();
    solver.setAcceptedIterateObserver([&](int iteration, const double *deviceX, double, double) {
      EXPECT_EQ(iteration, observations++);
      checkCuda(cudaMemcpy(lastObserved.data(), deviceX, 2 * sizeof(double), cudaMemcpyDeviceToHost), "test accepted plateau state");
    });
    ES::V2d x = ES::V2d::Zero();
    const double tolerance = 1e-8 / scale;
    EXPECT_EQ(solver.solve(x.data(), 200, tolerance), NO::CudaNewtonSolver::SOLVE_NOT_CONVERGED);
    const auto &d = solver.getDiagnostics();
    EXPECT_EQ(d.stopReason, "representational_stagnation");
    EXPECT_EQ(d.stagnationTriggers, 1);
    EXPECT_GT(d.residualCorrectionAttempts, 0);
    EXPECT_EQ(d.residualCorrectionSteps, 0);
    EXPECT_LT(d.newtonIterations, 200);
    EXPECT_EQ(observations, d.newtonIterations + 1);
    EXPECT_EQ(x, lastObserved);
    EXPECT_GT(d.finalGradientInfinityNorm, tolerance);
    EXPECT_LE(d.lastStagnationStepRatio, 16 * std::numeric_limits<double>::epsilon());
    EXPECT_LE(d.lastStagnationNorm2RelativeChange, 64 * std::numeric_limits<double>::epsilon());
    EXPECT_LE(d.lastStagnationInfinityRelativeChange, 64 * std::numeric_limits<double>::epsilon());
    EXPECT_DOUBLE_EQ(CurvedPathQuadratic::bound(ES::V2d::Zero(), x / scale), 1);
    EXPECT_FALSE(energy->evaluatedHiddenPoint);
  }
}

TEST_F(CudaNewtonTest, StagnationWindowPreservesTinyPhysicalProgressAndRecovery)
{
  NO::CudaNewtonSolver::SolverParam params;
  params.stagnationWindow = 4;
  // This fixture intentionally returns only a half-step certificate. Its
  // ordinary Newton path remains valid with the optional input guard off.
  auto small = std::make_shared<SmallDisplacementQuadratic>();
  NO::CudaNewtonSolver tiny({ small }, params);
  double scalar = 0;
  EXPECT_EQ(tiny.solve(&scalar, 40, 1e-9), NO::CudaNewtonSolver::SOLVE_CONVERGED);
  EXPECT_NEAR(scalar, 1e-15, 1e-24);
  EXPECT_EQ(tiny.getDiagnostics().stagnationTriggers, 0);

  params.checkInputChord = true;
  params.checkTrialChord = true;
  auto quantized = std::make_shared<QuantizedResidualQuadratic>();
  NO::CudaNewtonSolver recovery({ quantized }, params);
  ES::V2d x(-0.016, 0);
  EXPECT_EQ(recovery.solve(x.data(), 2, 1e-6), NO::CudaNewtonSolver::SOLVE_CONVERGED);
  EXPECT_GT(recovery.getDiagnostics().residualCorrectionSteps, 0);
  EXPECT_LT(recovery.getDiagnostics().finalGradientInfinityNorm, 1e-6);
  EXPECT_EQ(recovery.getDiagnostics().stagnationTriggers, 0);
}

TEST_F(CudaNewtonTest, StagnationWindowRetainsClippedOneUlpResidualProgress)
{
  // The represented step is only one ULP after an input-chord rejection, but
  // reduces the actual gradient by 25%. Even the opt-in residual threshold
  // with a window of one must preserve it and the original convergence test.
  for (double relativeTolerance : { 64 * std::numeric_limits<double>::epsilon(), 1e-10 }) {
    SCOPED_TRACE("relative tolerance=" + std::to_string(relativeTolerance));
    for (double tolerance : { 1e-8, 7.0 }) {
      NO::CudaNewtonSolver::SolverParam params;
      params.checkInputChord = true;
      params.stagnationWindow = 1;
      params.stagnationRelativeTolerance = relativeTolerance;
      auto energy = std::make_shared<RoundedBoundaryQuadratic>();
      NO::CudaNewtonSolver solver({ energy }, params);
      double x = 1;
      const auto expected = tolerance == 7 ? NO::CudaNewtonSolver::SOLVE_CONVERGED : NO::CudaNewtonSolver::SOLVE_NOT_CONVERGED;
      EXPECT_EQ(solver.solve(&x, 1, tolerance), expected);
      const auto &d = solver.getDiagnostics();
      EXPECT_GT(d.inputChordRejections, 0);
      EXPECT_EQ(d.newtonIterations, 1);
      EXPECT_EQ(d.stagnationTriggers, 0);
      EXPECT_DOUBLE_EQ(x, std::nextafter(1.0, 2.0));
      EXPECT_LT(d.finalGradientInfinityNorm, 7);
      EXPECT_GT(d.finalGradientInfinityNorm, 6);
    }
  }
}

TEST_F(CudaNewtonTest, StagnationRelativeToleranceHandlesResidualVariationWithoutConvergence)
{
  constexpr int cap = 12;
  constexpr double tolerance = 1e-8;
  const double originalThreshold = 64 * std::numeric_limits<double>::epsilon();
  for (bool optIn : { false, true }) {
    SCOPED_TRACE(optIn ? "opt-in residual variation" : "original strict threshold");
    NO::CudaNewtonSolver::SolverParam params;
    EXPECT_DOUBLE_EQ(params.stagnationRelativeTolerance, originalThreshold);
    params.checkInputChord = true;
    params.stagnationWindow = 4;
    if (optIn)
      params.stagnationRelativeTolerance = 1e-10;
    auto energy = std::make_shared<VaryingResidualBoundaryQuadratic>();
    NO::CudaNewtonSolver solver({ energy }, params);
    int observations = 0;
    ES::V2d lastObserved = ES::V2d::Zero();
    solver.setAcceptedIterateObserver([&](int iteration, const double *deviceX, double, double) {
      EXPECT_EQ(iteration, observations++);
      checkCuda(cudaMemcpy(lastObserved.data(), deviceX, 2 * sizeof(double), cudaMemcpyDeviceToHost), "test varying accepted state");
      EXPECT_TRUE(RoundedBoundaryQuadratic::feasible(lastObserved[0]));
    });
    ES::V2d x(1, 0.0625);
    EXPECT_EQ(solver.solve(x.data(), cap, tolerance), NO::CudaNewtonSolver::SOLVE_NOT_CONVERGED);
    const auto &d = solver.getDiagnostics();
    EXPECT_GT(d.inputChordRejections, 0);
    EXPECT_GT(d.finalGradientInfinityNorm, 0.9);
    EXPECT_EQ(observations, d.newtonIterations + 1);
    EXPECT_EQ(x, lastObserved);
    EXPECT_LE(d.lastStagnationStepRatio, 16 * std::numeric_limits<double>::epsilon());
    if (optIn) {
      EXPECT_EQ(d.stopReason, "representational_stagnation");
      EXPECT_EQ(d.stagnationTriggers, 1);
      EXPECT_LT(d.newtonIterations, cap);
      EXPECT_GT(d.residualCorrectionAttempts, 0);
      EXPECT_EQ(d.residualCorrectionSteps, 0);
      EXPECT_LE(d.lastStagnationNorm2RelativeChange, params.stagnationRelativeTolerance);
      EXPECT_LE(d.lastStagnationInfinityRelativeChange, params.stagnationRelativeTolerance);
    }
    else {
      EXPECT_EQ(d.stopReason, "iteration_limit");
      EXPECT_EQ(d.stagnationTriggers, 0);
      EXPECT_EQ(d.newtonIterations, cap);
      EXPECT_GT(d.lastStagnationNorm2RelativeChange, originalThreshold);
      EXPECT_GT(d.lastStagnationInfinityRelativeChange, originalThreshold);
    }
  }
}

TEST_F(CudaNewtonTest, InitialFeasibilityPrecedesConvergenceAndObservation)
{
  for (double bound : { 0.0, 0.5, -1.0, 2.0,
         std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN() }) {
    for (int maxIterations : { 0, 5 }) {
      auto energy = std::make_shared<InvalidInitialQuadratic>(bound);
      NO::CudaNewtonSolver solver({ energy });
      int observations = 0;
      solver.setAcceptedIterateObserver([&](int, const double *, double, double) { ++observations; });
      double x = 10;  // Exact unconstrained stationary point.
      EXPECT_EQ(solver.solve(&x, maxIterations, 1e-6), NO::CudaNewtonSolver::SOLVE_NUMERICAL_FAILURE);
      EXPECT_EQ(solver.getDiagnostics().stopReason, "infeasible_initial_state");
      EXPECT_EQ(solver.getDiagnostics().newtonIterations, 0);
      EXPECT_EQ(energy->initialChecks, 1);
      EXPECT_EQ(energy->energyEvaluations, 0);
      EXPECT_EQ(observations, 0);
      EXPECT_DOUBLE_EQ(x, 10);
    }
  }
  auto disabledGuard = std::make_shared<InvalidInitialQuadratic>(1);
  NO::CudaNewtonSolver solver({ disabledGuard });
  int observations = 0;
  solver.setAcceptedIterateObserver([&](int iteration, const double *, double, double) {
    EXPECT_EQ(iteration, 0);
    ++observations;
  });
  double x = 10;
  EXPECT_EQ(solver.solve(&x, 0, 1e-6), NO::CudaNewtonSolver::SOLVE_CONVERGED);
  EXPECT_EQ(disabledGuard->initialChecks, 1);
  EXPECT_EQ(disabledGuard->energyEvaluations, 1);
  EXPECT_EQ(observations, 1);
}

TEST_F(CudaNewtonTest, StagnationWindowRejectsInvalidLength)
{
  auto energy = std::make_shared<BoundedQuadratic>(1);
  for (int window : { -1, 65 }) {
    NO::CudaNewtonSolver::SolverParam params;
    params.stagnationWindow = window;
    EXPECT_THROW(NO::CudaNewtonSolver({ energy }, params), std::invalid_argument);
  }
}

TEST_F(CudaNewtonTest, StagnationRelativeToleranceRejectsInvalidValues)
{
  auto energy = std::make_shared<BoundedQuadratic>(1);
  for (double value : { -1.0, 0.0, 1.0, 2.0,
         std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN() }) {
    NO::CudaNewtonSolver::SolverParam params;
    params.stagnationRelativeTolerance = value;
    EXPECT_THROW(NO::CudaNewtonSolver({ energy }, params), std::invalid_argument);
  }
}

}  // namespace
