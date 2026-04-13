#include "TRBDF2TimeIntegratorHelper.h"
#include "TRBDF2TimeIntegrator.h"

#include <tbb/parallel_for.h>
#include <tbb/partitioner.h>

#include <iostream>

using namespace pgo;
using namespace pgo::Simulation;

namespace ES = pgo::EigenSupport;

TRBDF2TimeIntegratorEnergy::TRBDF2TimeIntegratorEnergy(TRBDF2TimeIntegrator *integrator_, const ES::SpMatD &A_, const ES::VXd &b_):
  intg(integrator_), A(A_), b(b_)
{
}

double TRBDF2TimeIntegratorEnergy::func(ES::ConstRefVecXd x) const
{
  // x is u (the full position), not du
  //std::cout << "Energy: ";
  // 0.5 A u^2
  double energy = ES::vTMv(A, x, intg->temp0, 0) * 0.5;

  //std::cout << energy << ',';
  for (size_t i = 0; i < intg->implicitModelsAll.size(); i++) {
    // elastic(u)
    double tempEnergy = intg->implicitModelsAll[i]->func(x);
    energy += tempEnergy;
  }

  // b^T u
  // double last_term = cblas_ddot(intg->n3, x.data(), 1, b.data(), 1);
  double last_term = x.dot(b);

  energy += last_term;
  // std::cout << last_term << std::endl;

  return energy;
}

void TRBDF2TimeIntegratorEnergy::gradient(ES::ConstRefVecXd x, ES::RefVecXd grad) const
{
  // x is u (the full position)
  // A u
  ES::mv(A, x, grad, 0);

  for (size_t i = 0; i < intg->implicitModelsAll.size(); i++) {
    ES::VXd &fint = *intg->implicitModelsAll_fint[i];

    // fint(u)
    intg->implicitModelsAll[i]->gradient(x, fint);
    // cblas_daxpy(intg->n3, 1.0, fint.data(), 1, grad.data(), 1);
    grad += fint;
  }

  // b
  // cblas_daxpy(intg->n3, 1.0, b.data(), 1, grad.data(), 1);
  grad += b;
}

void TRBDF2TimeIntegratorEnergy::hessian(ES::ConstRefVecXd x, ES::SpMatD &hess) const
{
  // x is u (the full position)
  memset(hess.valuePtr(), 0, sizeof(double) * hess.nonZeros());

  for (size_t i = 0; i < intg->implicitModelsAll.size(); i++) {
    ES::SpMatD &K = *intg->implicitModelsAll_K[i];
    const ES::SpMatI &mapping = *intg->implicitModelsAll_Kmaping[i];

    // beta/h K + K
    intg->implicitModelsAll[i]->hessian(x, K);
    double scale = 1.0;

    ES::addSmallToBig(scale, K, hess, 1.0, mapping, 1);
  }

  // cblas_daxpy((int)A.nonZeros(), 1.0, A.valuePtr(), 1, hess.valuePtr(), 1);
  (ES::Mp<ES::VXd>(hess.valuePtr(), hess.nonZeros())) += ES::Mp<const ES::VXd>(A.valuePtr(), A.nonZeros());
}

void TRBDF2TimeIntegratorEnergy::getDOFs(std::vector<int> &dofs) const
{
  dofs = intg->allDOFs;
}

void TRBDF2TimeIntegratorEnergy::createHessian(ES::SpMatD &hess) const
{
  hess = A;
}

int TRBDF2TimeIntegratorEnergy::getNumDOFs() const
{
  return intg->n3;
}

void TRBDF2TimeIntegratorEnergy::printImplicitEnergy(ES::ConstRefVecXd x) const
{
  // x is u (the full position)
  //std::cout << "Energy: ";
  // 0.5 Au^2
  double energy = ES::vTMv(A, x, intg->temp0, 0) * 0.5;

  // + b^T u
  // double last_term = cblas_ddot(intg->n3, x.data(), 1, b.data(), 1);
  double last_term = x.dot(b);

  energy += last_term;

  std::cout << "  main: " << energy << '\n';

  //std::cout << energy << ',';
  for (size_t i = 0; i < intg->implicitModelsAll.size(); i++) {
    // elastic(u)
    double tempEnergy = intg->implicitModelsAll[i]->func(x);
    std::cout << "  sub " << i << ": " << tempEnergy << '\n';
  }
}

double TRBDF2TimeIntegratorEnergy::computeMaxStepSize(ES::ConstRefVecXd x, ES::ConstRefVecXd dx) const
{
  double maxStepSize = 1.0;
  for (size_t i = 0; i < intg->implicitModelsAll.size(); i++) {
    double s = intg->implicitModelsAll[i]->computeMaxStepSize(x, dx);
    if (s < maxStepSize)
      maxStepSize = s;
  }
  return maxStepSize;
}
