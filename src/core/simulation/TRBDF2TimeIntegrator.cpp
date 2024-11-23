#include "TRBDF2TimeIntegrator.h"
#include "TRBDF2TimeIntegratorHelper.h"

#include "finiteDifference.h"

#include <tbb/parallel_for.h>

#include <numeric>
#include <iostream>

using namespace pgo;
using namespace pgo::NonlinearOptimization;
using namespace pgo::Simulation;

namespace ES = pgo::EigenSupport;

TRBDF2TimeIntegrator::TRBDF2TimeIntegrator(const EigenSupport::SpMatD &massMatrix,
  std::shared_ptr<const NonlinearOptimization::PotentialEnergy> elasticPotential,
  double massDampingCoeff, double stiffnessDampingCoeff, double gamma, double t, int numIter, double eps_):
  TimeIntegrator(massMatrix, elasticPotential, massDampingCoeff, stiffnessDampingCoeff, t, numIter, eps_),
  y(gamma)
{
  A1 = K;
  A2 = K;
  D = K;

  b1.setZero(n3);
  b2.setZero(n3);

  z1.setZero(n3);
  z2.setZero(n3);

  zero.setZero(n3);

  temp0 = q;
  qz = q;

  qy = q;
  qvely = q;
  qaccy = q;

  trEnergy = std::make_shared<TRBDF2TimeIntegratorEnergy>(this, A1, b1);
  bdf2Energy = std::make_shared<TRBDF2TimeIntegratorEnergy>(this, A2, b2);

  updateCoeffs();

  solver[0] = std::make_shared<TimeIntegratorSolver>();
  solver[1] = std::make_shared<TimeIntegratorSolver>();
}

void TRBDF2TimeIntegrator::setGamma(double gamma)
{
  y = gamma;
  updateCoeffs();
}

void TRBDF2TimeIntegrator::setTimestep(double t)
{
  TimeIntegrator::setTimestep(t);

  updateCoeffs();
}

std::shared_ptr<const NonlinearOptimization::PotentialEnergy> TRBDF2TimeIntegrator::getInternalEnergy() const
{
  return nullptr;
}

void TRBDF2TimeIntegrator::updateCoeffs()
{
  alpha = 2.0 / (y * timestep);

  beta[0] = (2.0 - y) / (y * (1 - y) * (1 - y) * timestep * timestep);
  beta[1] = -beta[0];

  beta[2] = (1.0 - y) / (y * timestep);
  beta[3] = -1.0 / (y * (1.0 - y) * timestep);

  beta[4] = (2.0 - y) * (2.0 - y) / ((1.0 - y) * (1.0 - y) * timestep * timestep);

  beta[5] = 1.0 / (y * (1.0 - y) * timestep);
  beta[6] = -beta[5];

  beta[7] = (2.0 - y) / ((1.0 - y) * timestep);
}

void TRBDF2TimeIntegrator::updateD()
{
  if (generalForceModelChanged)
    D = hessianAll;

  // D = 0;
  memset(D.valuePtr(), 0, sizeof(double) * D.nonZeros());

  for (size_t i = 0; i < implicitModelsAll.size(); i++) {
    const ES::SpMatD &M = *implicitModelsAll_M[i];
    const ES::SpMatI &mapping = *implicitModelsAll_Kmaping[i];

    // std::cout << massDampingParamsAll[i] << ',';

    // D += dM * massM
    if (massDampingParamsAll[i] > 0)
      ES::addSmallToBig(massDampingParamsAll[i], M, D, 1.0, mapping, 1);
  }
  // std::cout << D.norm() << ',';

  for (size_t i = 0; i < implicitModelsAll.size(); i++) {
    ES::SpMatD &curK = *implicitModelsAll_K[i];
    const ES::SpMatI &mapping = *implicitModelsAll_Kmaping[i];

    // D += dK * K
    if (dampingParamsAll[i] > 0) {
      implicitModelsAll[i]->hessian(q, curK);
      ES::addSmallToBig(dampingParamsAll[i], curK, D, 1.0, mapping, 1);
    }
  }

  // std::cout << D.norm() << std::endl;
}

void TRBDF2TimeIntegrator::updateA1()
{
  if (generalForceModelChanged)
    A1 = hessianAll;

  // A = 0;
  memset(A1.valuePtr(), 0, sizeof(double) * A1.nonZeros());

  // A = alpha * alpha * M + alpha * D

  // A += alpha * alpha * M
  ES::addSmallToBig(alpha * alpha, MasK, A1, 1.0, Kmapping, 1);
  // A += alpha * D
  // cblas_daxpy((int)A1.nonZeros(), alpha, D.valuePtr(), 1, A1.valuePtr(), 1);
  (ES::Mp<ES::VXd>(A1.valuePtr(), A1.nonZeros())) += ES::Mp<ES::VXd>(D.valuePtr(), D.nonZeros()) * alpha;
}

void TRBDF2TimeIntegrator::updateA2()
{
  if (generalForceModelChanged)
    A2 = hessianAll;

  // A = 0;
  memset(A2.valuePtr(), 0, sizeof(double) * A2.nonZeros());

  // A =  beta5 * M + beta8 * D

  // A += beta5 * M
  ES::addSmallToBig(beta[4], MasK, A2, 1.0, Kmapping, 1);
  // A += beta8 * D
  // cblas_daxpy((int)A2.nonZeros(), beta[7], D.valuePtr(), 1, A2.valuePtr(), 1);
  (ES::Mp<ES::VXd>(A2.valuePtr(), A2.nonZeros())) += ES::Mp<ES::VXd>(D.valuePtr(), D.nonZeros()) * beta[7];
}

void TRBDF2TimeIntegrator::updateb1()
{
  // b1 = 2 alpha M qvel + M qacc + D qvel + fext

  // b1 = 2 alpha M qvel
  ES::mv(MasK, qvel, b1);
  // cblas_dscal(n3, 2.0 * alpha, b1.data(), 1);
  b1 *= 2.0 * alpha;

  // b1 += M qacc
  ES::mv(MasK, qacc, b1, 1.0, 1.0);

  // b1 += D qvel
  ES::mv(D, qvel, b1, 1.0, 1.0);

  // += fext
  // cblas_daxpy(n3, 1.0, fext.data(), 1, b1.data(), 1);
  b1 += f_ext;

  // cblas_dscal(n3, -1.0, b1.data(), 1);
  b1 *= -1;
}

void TRBDF2TimeIntegrator::updateb2()
{
  // b2 = M (b1 q + b2 qy + b3 qvel + b4 qvely)
  temp0.noalias() = beta[0] * q;
  temp0 += beta[1] * qy;

  temp0 += beta[2] * qvel;
  temp0 += beta[3] * qvely;

  ES::mv(MasK, temp0, b2);

  // b2 += D (b6 q + b7 qy)
  temp0.noalias() = beta[5] * q;
  temp0 += beta[6] * qy;

  ES::mv(D, temp0, b2, 1.0, 1.0);

  // b2 -= fext
  // cblas_daxpy(n3, -1.0, fext.data(), 1, b2.data(), 1);
  b2 -= f_ext;
}

void TRBDF2TimeIntegrator::solve(ES::VXd &x, std::shared_ptr<TRBDF2TimeIntegratorEnergy> eng, int verbose, int printResidual)
{
  bool needRenew = (constraintsChanged || generalForceModelChanged);
  solverRet = solver[stage]->solve(needRenew, x, g, lambda, deltauRangeLow, deltauRangeHi,
    constraintsRangeLow, constraintsRangeHi, eng, constraints,
    nIter, eps, verbose, solverConfigFilename.length() ? solverConfigFilename.c_str() : nullptr,
    solverOption);

  if (printResidual) {
    ES::VXd residual(n3), rhs = ES::VXd::Zero(n3 - fixedDOFs.size());
    residual.setZero();
    eng->gradient(x, residual);
    ES::transferBigToSmall(residual, rhs, rhsb2s);
    std::cout << "    T" << timestepID << ": ||g||=" << rhs.norm() << "; Solver Ret: " << solverRet << std::endl;

    std::cout << "    Energy components:\n";
    eng->printImplicitEnergy(x);

    if (constraints) {
      std::cout << "    Constraints norm:" << lambda.norm() << std::endl;
    }
  }

  if (finiteDifferenceTestFlag)
    finiteDifferenceTest(x, q);
}

void TRBDF2TimeIntegrator::doTimestep(int updateq, int verbose, int printResidual)
{
  assembleImplicitModels();

  updateD();
  updateA1();
  updateb1();

  stage = 0;

  for (int i = 0; i < n3; i++) {
    z1[i] = deltauInitial[i];
  }

  if (generalForceModelChanged) {
    // do nothing
  }

  solve(z1, trEnergy, verbose, printResidual);

  // qy = q + z
  // qvely = alpha z - qvel
  // qaccy = alpha^2 z - 2 alpha qvel - qacc
  tbb::parallel_for(
    0, n3, [&](int i) {
      qy[i] = q[i] + z1[i];
      qvely[i] = alpha * z1[i] - qvel[i];
      qaccy[i] = alpha * alpha * z1[i] - 2.0 * alpha * qvel[i] - qacc[i];
    },
    tbb::static_partitioner());

  // BDF 2
  if (y < 1 - 1e-9) {
    stage = 1;

    updateA2();
    updateb2();

    for (int i = 0; i < n3; i++) {
      // z2 can be z1 too.
      z2[i] = deltauInitial[i];
    }

    if (generalForceModelChanged) {
      // do nothing
    }

    solve(z2, bdf2Energy, verbose, printResidual);

    // q1 = q + z
    // qvel1 = b6 q + b7 qy + b8 z
    // qacc1 = b1 q + b2 qy + b3 qvel + b4 qvely + b5 z
    tbb::parallel_for(
      0, n3, [&](int i) {
        q1[i] = q[i] + z2[i];
        qvel1[i] = beta[5] * q[i] + beta[6] * qy[i] + beta[7] * z2[i];
        qacc1[i] = beta[0] * q[i] + beta[1] * qy[i] + beta[2] * qvel[i] + beta[3] * qvely[i] + beta[4] * z2[i];
      },
      tbb::static_partitioner());

    if (printResidual) {
      ES::VXd residual = MasK * qacc1 + D * qvel1;
      for (size_t i = 0; i < implicitModelsAll.size(); i++) {
        ES::VXd grad(n3);
        grad.setZero();
        implicitModelsAll[i]->gradient(q1, grad);
        residual += grad;
      }
      ES::VXd res1(n3 - fixedDOFs.size());
      res1.setZero();
      ES::transferBigToSmall(residual, res1, rhsb2s);
      std::cout << "residual:" << res1.norm() << std::endl;
    }
  }
  else {
    tbb::parallel_for(
      0, n3, [&](int i) {
        q1[i] = qy[i];
        qvel1[i] = qvely[i];
        qacc1[i] = qaccy[i];
      },
      tbb::static_partitioner());
  }

  if (updateq) {
    proceedTimestep();
  }

  TimeIntegrator::doTimestep(updateq, verbose, printResidual);
}

void TRBDF2TimeIntegrator::finiteDifferenceTestIntegratorEnergy(ES::ConstRefVecXd x) const
{
  FiniteDifference fd(FiniteDifference::M_FIVE_POINT, 1e-7);
  std::cout << "Test integrator energy " << std::endl;
  if (stage == 0)
    fd.testEnergy(trEnergy, true, true, -1, x.data());
  else if (stage == 1)
    fd.testEnergy(bdf2Energy, true, true, -1, x.data());
}
