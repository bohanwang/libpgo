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
  //std::cout << "Energy: ";
  // 0.5 A z^2
  double energy = ES::vTMv(A, x, intg->temp0, 0) * 0.5;

  //std::cout << energy << ',';
  //// compute current u
  tbb::parallel_for(
    0, intg->n3, [&](int i) {
      intg->qz[i] = intg->q[i] + x[i];
    },
    tbb::static_partitioner());

  for (size_t i = 0; i < intg->implicitModelsAll.size(); i++) {
    // elastic(q+z)
    double tempEnergy = intg->implicitModelsAll[i]->func(intg->qz);
    energy += tempEnergy;
  }

  // b^T z
  // double last_term = cblas_ddot(intg->n3, x.data(), 1, b.data(), 1);
  double last_term = x.dot(b);

  energy += last_term;
  // std::cout << last_term << std::endl;

  return energy;
}

void TRBDF2TimeIntegratorEnergy::gradient(ES::ConstRefVecXd x, ES::RefVecXd grad) const
{
  // A z
  ES::mv(A, x, grad, 0);

  tbb::parallel_for(
    0, intg->n3, [&](int i) {
      intg->qz[i] = intg->q[i] + x[i];
    },
    tbb::static_partitioner());

  for (size_t i = 0; i < intg->implicitModelsAll.size(); i++) {
    ES::VXd &fint = *intg->implicitModelsAll_fint[i];

    // fint(q+z)
    intg->implicitModelsAll[i]->gradient(intg->qz, fint);
    // cblas_daxpy(intg->n3, 1.0, fint.data(), 1, grad.data(), 1);
    grad += fint;
  }

  // b
  // cblas_daxpy(intg->n3, 1.0, b.data(), 1, grad.data(), 1);
  grad += b;
}

void TRBDF2TimeIntegratorEnergy::hessian(ES::ConstRefVecXd x, ES::SpMatD &hess) const
{
  memset(hess.valuePtr(), 0, sizeof(double) * hess.nonZeros());

  tbb::parallel_for(
    0, intg->n3, [&](int i) {
      intg->qz[i] = intg->q[i] + x[i];
    },
    tbb::static_partitioner());

  for (size_t i = 0; i < intg->implicitModelsAll.size(); i++) {
    ES::SpMatD &K = *intg->implicitModelsAll_K[i];
    const ES::SpMatI &mapping = *intg->implicitModelsAll_Kmaping[i];

    // beta/h K + K
    intg->implicitModelsAll[i]->hessian(intg->qz, K);
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
  //std::cout << "Energy: ";
  // 0.5 Az^2
  double energy = ES::vTMv(A, x, intg->temp0, 0) * 0.5;

  // + b_y^T z
  // double last_term = cblas_ddot(intg->n3, x.data(), 1, b.data(), 1);
  double last_term = x.dot(b);

  energy += last_term;

  std::cout << "  main: " << energy << '\n';

  //std::cout << energy << ',';
  //// compute current u
  tbb::parallel_for(
    0, intg->n3, [&](int i) {
      intg->qz[i] = intg->q[i] + x[i];
    },
    tbb::static_partitioner());

  for (size_t i = 0; i < intg->implicitModelsAll.size(); i++) {
    // elastic(q+z)
    double tempEnergy = intg->implicitModelsAll[i]->func(intg->qz);
    std::cout << "  sub " << i << ": " << tempEnergy << '\n';
  }
}
