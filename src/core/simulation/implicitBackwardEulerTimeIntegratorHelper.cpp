#include "implicitBackwardEulerTimeIntegratorHelper.h"
#include "implicitBackwardEulerTimeIntegrator.h"

#include <tbb/parallel_for.h>

#include <iostream>

using namespace pgo;
using namespace pgo::Simulation;

namespace ES = pgo::EigenSupport;

ImplicitBackwardEulerEnergy::ImplicitBackwardEulerEnergy(ImplicitBackwardEulerTimeIntegrator *integrator_):
  intg(integrator_)
{
}

double ImplicitBackwardEulerEnergy::func(ES::ConstRefVecXd x) const
{
  //std::cout << "Energy: ";
  // 0.5 Az^2
  double energy = ES::vTMv(intg->A, x, intg->temp0, 0) * 0.5;
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
    // std::cout << tempEnergy << ',';
      energy += tempEnergy;
  }

  // - b_y^T z
  // double last_term = cblas_ddot(intg->n3, x.data(), 1, intg->b.data(), 1);
  double last_term = x.dot(intg->b);

  energy -= last_term;
  // std::cout << last_term << std::endl;

  return energy;
}

void ImplicitBackwardEulerEnergy::gradient(ES::ConstRefVecXd x, ES::RefVecXd grad) const
{
  // Az
  ES::mv(intg->A, x, grad, 0);

  tbb::parallel_for(
    0, intg->n3, [&](int i) {
      intg->qz[i] = intg->q[i] + x[i];
    },
    tbb::static_partitioner());

  for (size_t i = 0; i < intg->implicitModelsAll.size(); i++) {
    ES::VXd &fint = *intg->implicitModelsAll_fint[i];
    // ES::SpMatD &K = *intg->implicitModelsAll_K[i];

    // fint(q+z)
    intg->implicitModelsAll[i]->gradient(intg->qz, fint);
    grad += fint;
  }

  //cblas_daxpy(intg->n3, -1.0, intg->b.data(), 1, grad.data(), 1);
  grad -= intg->b;
}

void ImplicitBackwardEulerEnergy::hessian(ES::ConstRefVecXd x, ES::SpMatD &hess) const
{
  // constexpr double eps = 1e-8;
  memset(hess.valuePtr(), 0, sizeof(double) * hess.nonZeros());

  tbb::parallel_for(
    0, intg->n3, [&](int i) {
      intg->qz[i] = intg->q[i] + x[i];
    },
    tbb::static_partitioner());

  for (size_t i = 0; i < intg->implicitModelsAll.size(); i++) {
    ES::SpMatD &K = *intg->implicitModelsAll_K[i];
    // ES::SpMatD &K1 = *intg->implicitModelsAll_K1[i];
    const ES::SpMatI &mapping = *intg->implicitModelsAll_Kmaping[i];

    // beta/h K + K
    intg->implicitModelsAll[i]->hessian(intg->qz, K);
    double scale = 1.0;

    ES::addSmallToBig(scale, K, hess, 1.0, mapping, 1);
  }

  (ES::Mp<ES::VXd>(hess.valuePtr(), hess.nonZeros())) +=ES::Mp<const ES::VXd>(intg->A.valuePtr(), intg->A.nonZeros());
}

void ImplicitBackwardEulerEnergy::getDOFs(std::vector<int> &dofs) const
{
  dofs = intg->allDOFs;
}

void ImplicitBackwardEulerEnergy::createHessian(ES::SpMatD &hess) const
{
  hess = intg->hessianAll;
}

int ImplicitBackwardEulerEnergy::getNumDOFs() const
{
  return intg->n3;
}

void ImplicitBackwardEulerEnergy::printImplicitEnergy(ES::ConstRefVecXd x) const
{
  //std::cout << "Energy: ";
  // 0.5 Az^2
  double energy = ES::vTMv(intg->A, x, intg->temp0, 0) * 0.5;
  // - b_y^T z
  // double last_term = cblas_ddot(intg->n3, x.data(), 1, intg->b.data(), 1);
  double last_term = x.dot(intg->b);
  energy -= last_term;

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
