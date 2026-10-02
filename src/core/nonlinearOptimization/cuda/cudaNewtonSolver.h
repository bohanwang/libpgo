#pragma once

#include "cudaPotentialEnergy.h"

#include <functional>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace pgo::NonlinearOptimization
{

class CudaNewtonSolver
{
public:
  // Values intentionally match NewtonSolver::SolveStatus.
  enum SolveStatus
  {
    SOLVE_CONVERGED = 0,
    SOLVE_NOT_CONVERGED = 1,
    SOLVE_NUMERICAL_FAILURE = 2,
  };

  struct SolverParam
  {
    int maxPCGIterations = 0;  // 0: twice the number of DOFs, capped at 10000
    // Zero retains host scalar checks. A positive batch size (at most 64)
    // performs the same per-iteration checks on the device, reading their
    // sticky terminal status once per batch. No convergence checks are skipped.
    int pcgReductionBatchSize = 0;
    // Zero disables the optional stagnation exit. After this many consecutive
    // input-chord-clipped, representationally tiny candidates without measured
    // residual progress, try bounded residual recovery or return not converged.
    int stagnationWindow = 0;
    // Relative change allowed in both residual norms over that window. This
    // controls retry detection only; the requested gradient tolerance is fixed.
    double stagnationRelativeTolerance = 64 * std::numeric_limits<double>::epsilon();
    double pcgRelativeTolerance = 1e-6;
    double pcgAbsoluteTolerance = 1e-12;
    int maxLineSearchIterations = 60;
    double armijoCoeff = 1e-4;
    // After nonpositive curvature, retry (H + shift I) dx = -g. Each retry
    // increases shift tenfold, starting at initialRegularization * max|diag H|.
    int maxRegularizationAttempts = 12;
    double initialRegularization = 1e-8;
    // Once this many steps have been accepted, try the unprojected tangent
    // before the projected operator. Negative values disable these attempts.
    int exactHessianStartIteration = 20;
    // Bounded coordinate recovery when rounding defeats the Newton direction.
    // Zero disables recovery. Accepted trials retain the normal energy and
    // feasibility safeguards and must improve the complete measured residual.
    int maxResidualCorrectionCandidates = 32;
    // Require every represented trial to have a feasible straight path from
    // the solve's input, in addition to each energy's Newton step bound. The
    // caller may retry a shorter physical timestep if this visibility constraint
    // prevents convergence. Only complete original-chord certificates pass.
    bool checkInputChord = false;
    // Certify the current iterate to the represented trial explicitly. This
    // checks the rounded endpoints even after clipping the search direction.
    bool checkTrialChord = false;
  };

  struct Diagnostics
  {
    int newtonIterations = 0;
    int pcgIterations = 0;
    int regularizationAttempts = 0;
    int roundoffAcceptedSteps = 0;
    int exactHessianAttempts = 0;
    int exactHessianDirections = 0;
    int exactHessianFallbacks = 0;
    int inputChordRejections = 0;
    int trialChordRejections = 0;
    int residualCorrectionAttempts = 0;
    int residualCorrectionCandidates = 0;
    int residualCorrectionSteps = 0;
    int stagnationTriggers = 0;
    double lastStagnationStepRatio = 0;
    double lastStagnationNorm2RelativeChange = 0;
    double lastStagnationInfinityRelativeChange = 0;
    double stagnationAnchorNorm2 = 0;
    double stagnationAnchorInfinity = 0;
    double maxRegularization = 0;
    double finalEnergy = 0;
    double finalGradientInfinityNorm = 0;
    double lastLinearResidualNorm = 0;
    double lastFeasibleAlpha = 1;
    double lastAcceptedAlpha = 0;
    std::string stopReason = "not_run";
  };

  explicit CudaNewtonSolver(std::vector<std::shared_ptr<CudaPotentialEnergy>> energies);
  CudaNewtonSolver(std::vector<std::shared_ptr<CudaPotentialEnergy>> energies, const SolverParam &parameters);
  ~CudaNewtonSolver();
  CudaNewtonSolver(const CudaNewtonSolver &) = delete;
  CudaNewtonSolver &operator=(const CudaNewtonSolver &) = delete;

  // Host ingress/egress only. All Newton/PCG vectors remain on the device
  // between these transfers. On failure, returns the last finite accepted x.
  int solve(double *hostX, int maxIterations, double gradientTolerance, int verbose = 0);
  const Diagnostics &getDiagnostics() const;
  // Called synchronously for the finite feasible input (iteration zero) and every
  // accepted update, including residual recovery. The device pointer is valid
  // only during the callback. No observer means no extra state transfers.
  void setAcceptedIterateObserver(std::function<void(int iteration, const double *deviceX, double energy, double gradientInf)> observer);

private:
  class Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace pgo::NonlinearOptimization
