#include "implicitBackwardEulerTimeIntegrator.h"
#include "implicitBackwardEulerTimeIntegratorHelper.h"
#include "timeIntegratorSolver.h"

#include "potentialEnergies.h"
#include "finiteDifference.h"

#include <tbb/parallel_for.h>

#include <numeric>
#include <iostream>

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

void ImplicitBackwardEulerTimeIntegrator::doTimestep(int updateq, int verbose, int printResidual)
{
  assembleImplicitModels();

  updateD();
  updateA();
  updateb();

  for (int i = 0; i < n3; i++) {
    z[i] = deltauInitial[i];
  }

  if (generalForceModelChanged) {
    // do nothing
  }

  bool needRenew = (constraintsChanged || generalForceModelChanged);
  solverRet = solver->solve(needRenew, z, g, lambda, deltauRangeLow, deltauRangeHi,
    constraintsRangeLow, constraintsRangeHi, eulerEnergy, constraints,
    nIter, eps, verbose, solverConfigFilename.length() ? solverConfigFilename.c_str() : nullptr,
    solverOption);

  if (printResidual) {
    ES::VXd residual(n3), rhs = ES::VXd::Zero(n3 - fixedDOFs.size());
    eulerEnergy->gradient(z, residual);

    ES::transferBigToSmall(residual, rhs, rhsb2s);
    std::cout << "    T" << timestepID << ": ||g||=" << rhs.norm() << "; Solver Ret: " << solverRet << std::endl;

    std::cout << "    Energy components:\n";
    eulerEnergy->printImplicitEnergy(z);
  }

  if (finiteDifferenceTestFlag)
    finiteDifferenceTest(z, q);

  // q1 = q + z
  // qvel1 = 1/h * z
  // qacc1 = 1/h(1/h * z - qvel)
  tbb::parallel_for(
    0, n3, [&](int i) {
      q1[i] = q[i] + z[i];
      qvel1[i] = z[i] / timestep;
      qacc1[i] = (qvel1[i] - qvel[i]) / timestep;
    },
    tbb::static_partitioner());

  if (updateq) {
    proceedTimestep();
  }

  TimeIntegrator::doTimestep(updateq, verbose, printResidual);
}

void ImplicitBackwardEulerTimeIntegrator::setSolution(ES::ConstRefVecXd newz)
{
  z = newz;
  // q1 = q + z
  // qvel1 = 1/h * z
  // qacc1 = 1/h(1/h * z - qvel)
  tbb::parallel_for(
    0, n3, [&](int i) {
      q1[i] = q[i] + z[i];
      qvel1[i] = z[i] / timestep;
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
  // b = fext + 1/h M qvel
  // 1/h M qvel
  ES::mv(MasK, qvel, b);
  // cblas_dscal(n3, 1.0 / timestep, b.data(), 1);
  b *= 1.0 / timestep;

  // += fext
  // cblas_daxpy(n3, 1.0, fext.data(), 1, b.data(), 1);
  b += f_ext;
}

void ImplicitBackwardEulerTimeIntegrator::finiteDifferenceTestIntegratorEnergy(ES::ConstRefVecXd x) const
{
  FiniteDifference fd(FiniteDifference::M_FIVE_POINT, 1e-7);
  std::cout << "Test integrator energy " << std::endl;
  fd.testEnergy(eulerEnergy, true, true, -1, x.data());
}
