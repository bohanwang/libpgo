#include "TRBDF2TimeIntegratorHelper.h"
#include "TRBDF2TimeIntegrator.h"
#include "deformationModelEnergy.h"
#include "embeddedSurfaceIPCPotentialEnergy.h"
#include "CIPC.h"

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

void TRBDF2TimeIntegratorEnergy::gradient_hessian(ES::ConstRefVecXd x, ES::RefVecXd grad, ES::SpMatD &hess) const
{
  ES::mv(A, x, grad, 0);

  hess = intg->hessianAll;
  memset(hess.valuePtr(), 0, sizeof(double) * hess.nonZeros());
  (ES::Mp<ES::VXd>(hess.valuePtr(), hess.nonZeros())) += ES::Mp<const ES::VXd>(A.valuePtr(), A.nonZeros());

  for (size_t i = 0; i < intg->implicitModelsAll.size(); i++) {
    ES::VXd &fint = *intg->implicitModelsAll_fint[i];
    fint.setZero();

    if (intg->implicitModelsAll[i]->isHessianTopologyFixed()) {
      ES::SpMatD &K = *intg->implicitModelsAll_K[i];
      const ES::SpMatI &mapping = *intg->implicitModelsAll_Kmaping[i];

      intg->implicitModelsAll[i]->gradient(x, fint);
      intg->implicitModelsAll[i]->hessian(x, K);
      ES::addSmallToBig(1.0, K, hess, 1.0, mapping, 1);
    }
    else {
      ES::SpMatD Ki;
      intg->implicitModelsAll[i]->gradient_hessian(x, fint, Ki);
      if (Ki.nonZeros())
        hess = hess + Ki;
    }

    grad += fint;
  }

  grad += b;
}

double TRBDF2TimeIntegratorEnergy::func_grad_hessian(ES::ConstRefVecXd x, ES::RefVecXd grad, ES::SpMatD &hess) const
{
  ES::mv(A, x, grad, 0);
  double energy = ES::vTMv(A, x, intg->temp0, 0) * 0.5;

  hess = intg->hessianAll;
  memset(hess.valuePtr(), 0, sizeof(double) * hess.nonZeros());
  (ES::Mp<ES::VXd>(hess.valuePtr(), hess.nonZeros())) += ES::Mp<const ES::VXd>(A.valuePtr(), A.nonZeros());

  for (size_t i = 0; i < intg->implicitModelsAll.size(); i++) {
    ES::VXd &fint = *intg->implicitModelsAll_fint[i];
    fint.setZero();

    if (intg->implicitModelsAll[i]->isHessianTopologyFixed()) {
      energy += intg->implicitModelsAll[i]->func(x);

      ES::SpMatD &K = *intg->implicitModelsAll_K[i];
      const ES::SpMatI &mapping = *intg->implicitModelsAll_Kmaping[i];
      intg->implicitModelsAll[i]->gradient(x, fint);
      intg->implicitModelsAll[i]->hessian(x, K);
      ES::addSmallToBig(1.0, K, hess, 1.0, mapping, 1);
    }
    else {
      ES::SpMatD Ki;
      energy += intg->implicitModelsAll[i]->func_grad_hessian(x, fint, Ki);
      if (Ki.nonZeros())
        hess = hess + Ki;
    }

    grad += fint;
  }

  energy += x.dot(b);
  grad += b;
  return energy;
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

int TRBDF2TimeIntegratorEnergy::isHessianTopologyFixed() const
{
  for (size_t i = 0; i < intg->implicitModelsAll.size(); i++) {
    if (!intg->implicitModelsAll[i]->isHessianTopologyFixed())
      return 0;
  }
  return 1;
}

void TRBDF2TimeIntegratorEnergy::hessianDirect(ES::ConstRefVecXd x, ES::SpMatD &hess) const
{
  // Start with hessianAll pattern (covers A + all fixed-topology models)
  hess = intg->hessianAll;
  memset(hess.valuePtr(), 0, sizeof(double) * hess.nonZeros());

  // Add A (mass/damping/timestep contribution)
  (ES::Mp<ES::VXd>(hess.valuePtr(), hess.nonZeros())) += ES::Mp<const ES::VXd>(A.valuePtr(), A.nonZeros());

  // Add fixed-topology models using efficient mapping
  for (size_t i = 0; i < intg->implicitModelsAll.size(); i++) {
    if (intg->implicitModelsAll[i]->isHessianTopologyFixed()) {
      ES::SpMatD &K = *intg->implicitModelsAll_K[i];
      if (K.nonZeros() == 0)
        continue;

      const ES::SpMatI &mapping = *intg->implicitModelsAll_Kmaping[i];
      intg->implicitModelsAll[i]->hessian(x, K);
      ES::addSmallToBig(1.0, K, hess, 1.0, mapping, 1);
    }
  }

  // Add non-fixed-topology models
  for (size_t i = 0; i < intg->implicitModelsAll.size(); i++) {
    if (!intg->implicitModelsAll[i]->isHessianTopologyFixed()) {
      ES::SpMatD Ki;
      intg->implicitModelsAll[i]->hessianDirect(x, Ki);
      if (Ki.nonZeros() == 0)
        continue;

      hess = hess + Ki;
    }
  }
}

NonlinearOptimization::MaxStepResult TRBDF2TimeIntegratorEnergy::computeMaxStepLimit(ES::ConstRefVecXd x, ES::ConstRefVecXd dx) const
{
  NonlinearOptimization::MaxStepResult result = NonlinearOptimization::MaxStepResult::unconstrained();
  for (size_t i = 0; i < intg->implicitModelsAll.size(); i++) {
    const auto &model = intg->implicitModelsAll[i];
    result = NonlinearOptimization::mergeMaxStepResults(result, model->computeMaxStepLimit(x, dx));
  }
  return result;
}
