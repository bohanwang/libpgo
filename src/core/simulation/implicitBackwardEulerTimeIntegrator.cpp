#include "implicitBackwardEulerTimeIntegrator.h"
#include "implicitBackwardEulerTimeIntegratorHelper.h"
#include "NewtonSolver.h"
#include "timeIntegratorSolver.h"

#include "potentialEnergies.h"
#include "finiteDifference.h"

#include <tbb/parallel_for.h>

#include <numeric>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace pgo;
using namespace pgo::NonlinearOptimization;
using namespace pgo::Simulation;

namespace ES = pgo::EigenSupport;

ImplicitBackwardEulerTimeIntegrator::ImplicitBackwardEulerTimeIntegrator(
  const EigenSupport::SpMatD &massMatrix,
  std::shared_ptr<const NonlinearOptimization::PotentialEnergy> elasticPotential,
  double massDampingCoeff, double stiffnessDampingCoeff, double t, int numIter, double eps_):
  TimeIntegrator(massMatrix, elasticPotential, massDampingCoeff, stiffnessDampingCoeff, t, numIter, eps_)
{
  A = K;
  D = K;

  b.setZero(n3);
  z.setZero(n3);
  zero.setZero(n3);

  temp0 = q;
  qz = q;
  qz1 = qz;
  qz2 = qz;

  eulerEnergy = std::make_shared<ImplicitBackwardEulerEnergy>(this);
  solver = std::make_shared<TimeIntegratorSolver>();
}

std::shared_ptr<const NonlinearOptimization::PotentialEnergy> ImplicitBackwardEulerTimeIntegrator::getInternalEnergy() const
{
  return eulerEnergy;
}

const NonlinearOptimization::SolveDiagnostics &ImplicitBackwardEulerTimeIntegrator::getLastSolveDiagnostics() const
{
  return solver->getLastSolveDiagnostics();
}

void ImplicitBackwardEulerTimeIntegrator::doTimestep(int updateq, int verbose, int printResidual)
{
  const int ret = tryTimestep(updateq, verbose, printResidual);
  if (ret != 0) {
    throw std::runtime_error(std::string("ImplicitBackwardEulerTimeIntegrator timestep failed with solverRet=") +
      NewtonSolver::solveStatusToString(ret));
  }
}

int ImplicitBackwardEulerTimeIntegrator::tryTimestep(int updateq, int verbose, int printResidual)
{
  assembleImplicitModels();

  updateD();
  updateA();
  updateb();

  // z is the optimization variable u (full position)
  // initialize based on initial guess mode
  if (initialGuessMode == InitialGuessMode::LAST_U_PLUS_VH) {
    for (int i = 0; i < n3; i++) {
      z[i] = q[i] + qvel[i] * timestep;
    }
  }
  else {
    // InitialGuessMode::LAST_U
    for (int i = 0; i < n3; i++) {
      z[i] = q[i];
    }
  }

  if (generalForceModelChanged) {
    // do nothing
  }

  bool needRenew = (constraintsChanged || generalForceModelChanged);
  if (verbose) {
    std::cout << "ImplicitBackwardEuler timestep begin: T" << timestepID
              << " dt=" << timestep << std::endl;
  }

  eulerEnergy->clearCachedImplicitEnergyComponents();
  solverRet = solver->solve(needRenew, z, g, lambda, uRangeLow, uRangeHi,
    constraintsRangeLow, constraintsRangeHi, eulerEnergy, constraints,
    nIter, eps, verbose, solverConfigFilename.length() ? solverConfigFilename.c_str() : nullptr,
    solverOption);

  const SolveDiagnostics &diagnostics = solver->getLastSolveDiagnostics();
  double residualNorm = 0.0;
  double residualMaxNorm = 0.0;
  if (printResidual || verbose || solverRet != 0) {
    if (diagnostics.hasFinalGradientStats) {
      residualNorm = diagnostics.finalGradientNorm;
      residualMaxNorm = diagnostics.finalGradientMaxNorm;
    }
    else {
      ES::VXd residual(n3), rhs = ES::VXd::Zero(n3 - fixedDOFs.size());
      residual.setZero();
      eulerEnergy->gradient(z, residual);
      ES::transferBigToSmall(residual, rhs, rhsb2s);
      residualNorm = rhs.norm();
      residualMaxNorm = rhs.cwiseAbs().maxCoeff();
    }
  }

  const bool acceptedTimestep = solverRet == 0 ||
    solverRet == static_cast<int>(NewtonSolver::SolveStatus::MaxIterations) ||
    solverRet == static_cast<int>(NewtonSolver::SolveStatus::StepTooSmall);

  if (solverRet != 0) {
    std::cout << "Warning: solverRet = " << NewtonSolver::solveStatusToString(solverRet) << "\n";
  }

  if (verbose) {
    std::cout << "ImplicitBackwardEuler timestep end: T" << timestepID
              << " solverRet=" << NewtonSolver::solveStatusToString(solverRet)
              << " residual=" << residualNorm
              << " residualMax=" << residualMaxNorm
              << " accepted=" << (acceptedTimestep ? "true" : "false")
              << std::endl;
  }

  if (printResidual) {
    std::cout << "    T" << timestepID << ": ||g||=" << residualNorm
              << "; Solver Ret: " << solverRet
              << " (" << NewtonSolver::solveStatusToString(solverRet) << ")" << std::endl;

    std::cout << "    Energy components:\n";
    eulerEnergy->printImplicitEnergy(z, diagnostics.hasFinalGradientStats);
  }

  if (!acceptedTimestep) {
    TimeIntegrator::doTimestep(0, verbose, printResidual);
    return solverRet;
  }

  if (finiteDifferenceTestFlag)
    finiteDifferenceTest(z, q);

  // z is now u (full position)
  // q1 = u
  // qvel1 = (u - q) / h
  // qacc1 = (qvel1 - qvel) / h
  tbb::parallel_for(
    0, n3, [&](int i) {
      q1[i] = z[i];
      qvel1[i] = (z[i] - q[i]) / timestep;
      qacc1[i] = (qvel1[i] - qvel[i]) / timestep;
    },
    tbb::static_partitioner());

  if (updateq) {
    proceedTimestep();
  }

  TimeIntegrator::doTimestep(updateq, verbose, printResidual);
  return 0;
}

void ImplicitBackwardEulerTimeIntegrator::setSolution(ES::ConstRefVecXd newz)
{
  z = newz;
  // z is u (full position)
  // q1 = u
  // qvel1 = (u - q) / h
  // qacc1 = (qvel1 - qvel) / h
  tbb::parallel_for(
    0, n3, [&](int i) {
      q1[i] = z[i];
      qvel1[i] = (z[i] - q[i]) / timestep;
      qacc1[i] = (qvel1[i] - qvel[i]) / timestep;
    },
    tbb::static_partitioner());
}

void ImplicitBackwardEulerTimeIntegrator::updateD()
{
  if (generalForceModelChanged)
    D = hessianAll;

  // D = 0;
  memset(D.valuePtr(), 0, sizeof(double) * D.nonZeros());

  // Damping is only supported for energies with fixed hessian topology.
  // Non-fixed-topology energies (e.g. CIPC contact) are skipped here.
  for (size_t i = 0; i < implicitModelsAll.size(); i++) {
    if (!implicitModelsAll[i]->isHessianTopologyFixed())
      continue;

    const ES::SpMatD &M = *implicitModelsAll_M[i];
    const ES::SpMatI &mapping = *implicitModelsAll_Kmaping[i];

    // D += dM * massM
    if (massDampingParamsAll[i] > 0)
      ES::addSmallToBig(massDampingParamsAll[i], M, D, 1.0, mapping, 1);
  }

  for (size_t i = 0; i < implicitModelsAll.size(); i++) {
    if (!implicitModelsAll[i]->isHessianTopologyFixed())
      continue;

    ES::SpMatD &curK = *implicitModelsAll_K[i];
    const ES::SpMatI &mapping = *implicitModelsAll_Kmaping[i];

    // D += dK * K
    if (dampingParamsAll[i] > 0) {
      implicitModelsAll[i]->hessian(q, curK);
      ES::addSmallToBig(dampingParamsAll[i], curK, D, 1.0, mapping, 1);
    }
  }

  std::cout << "D:" << D.coeffs().cwiseAbs().maxCoeff() << std::endl;

  // std::cout << D.norm() << std::endl;
}

void ImplicitBackwardEulerTimeIntegrator::updateA()
{
  if (generalForceModelChanged)
    A = hessianAll;

  // A = 0;
  memset(A.valuePtr(), 0, sizeof(double) * A.nonZeros());

  // A += (1/h)^2 M
  double s = 1.0 / (timestep * timestep);
  ES::addSmallToBig(s, MasK, A, 1.0, Kmapping, 1);

  // A += 1/h D
  s = 1.0 / timestep;

  // cblas_daxpy((int)A.nonZeros(), s, D.valuePtr(), 1, A.valuePtr(), 1);
  (ES::Mp<ES::VXd>(A.valuePtr(), A.nonZeros())) += ES::Mp<const ES::VXd>(D.valuePtr(), D.nonZeros()) * s;
}

void ImplicitBackwardEulerTimeIntegrator::updateb()
{
  // b = fext + 1/h M qvel + A q
  // (the A*q term comes from changing the optimization variable from du to u)
  // 1/h M qvel
  ES::mv(MasK, qvel, b);
  // cblas_dscal(n3, 1.0 / timestep, b.data(), 1);
  b *= 1.0 / timestep;

  // += fext
  // cblas_daxpy(n3, 1.0, fext.data(), 1, b.data(), 1);
  b += f_ext;

  // += A q (shift from du to u variable)
  ES::mv(A, q, b, 1.0, 1.0);
}

void ImplicitBackwardEulerTimeIntegrator::finiteDifferenceTestIntegratorEnergy(ES::ConstRefVecXd x) const
{
  FiniteDifference fd(FiniteDifference::M_FIVE_POINT, 1e-7);
  std::cout << "Test integrator energy " << std::endl;
  fd.testEnergy(eulerEnergy, true, true, -1, x.data());
}
