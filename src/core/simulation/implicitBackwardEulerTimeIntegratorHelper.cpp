#include "implicitBackwardEulerTimeIntegratorHelper.h"
#include "implicitBackwardEulerTimeIntegrator.h"
#include "deformationModelEnergy.h"
#include "ipc/embeddedSurfaceIPCPotentialEnergy.h"

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
  // x is u (the full position), not du
  // 0.5 Au^2
  double energy = ES::vTMv(intg->A, x, intg->temp0, 0) * 0.5;

  for (size_t i = 0; i < intg->implicitModelsAll.size(); i++) {
    // elastic(u)
    double tempEnergy = intg->implicitModelsAll[i]->func(x);
    energy += tempEnergy;
  }

  // - b^T u
  double last_term = x.dot(intg->b);

  energy -= last_term;

  return energy;
}

void ImplicitBackwardEulerEnergy::gradient(ES::ConstRefVecXd x, ES::RefVecXd grad) const
{
  // x is u (the full position)
  // Au
  ES::mv(intg->A, x, grad, 0);

  for (size_t i = 0; i < intg->implicitModelsAll.size(); i++) {
    ES::VXd &fint = *intg->implicitModelsAll_fint[i];
    // ES::SpMatD &K = *intg->implicitModelsAll_K[i];

    // fint(u)
    intg->implicitModelsAll[i]->gradient(x, fint);
    grad += fint;
  }

  //cblas_daxpy(intg->n3, -1.0, intg->b.data(), 1, grad.data(), 1);
  grad -= intg->b;
}

void ImplicitBackwardEulerEnergy::hessian(ES::ConstRefVecXd x, ES::SpMatD &hess) const
{
  // x is u (the full position)
  // constexpr double eps = 1e-8;
  memset(hess.valuePtr(), 0, sizeof(double) * hess.nonZeros());

  for (size_t i = 0; i < intg->implicitModelsAll.size(); i++) {
    ES::SpMatD &K = *intg->implicitModelsAll_K[i];
    // ES::SpMatD &K1 = *intg->implicitModelsAll_K1[i];
    const ES::SpMatI &mapping = *intg->implicitModelsAll_Kmaping[i];

    // beta/h K + K
    intg->implicitModelsAll[i]->hessian(x, K);
    double scale = 1.0;

    ES::addSmallToBig(scale, K, hess, 1.0, mapping, 1);
  }

  (ES::Mp<ES::VXd>(hess.valuePtr(), hess.nonZeros())) +=ES::Mp<const ES::VXd>(intg->A.valuePtr(), intg->A.nonZeros());
}

void ImplicitBackwardEulerEnergy::gradient_hessian(ES::ConstRefVecXd x, ES::RefVecXd grad, ES::SpMatD &hess) const
{
  ES::mv(intg->A, x, grad, 0);

  hess = intg->hessianAll;
  memset(hess.valuePtr(), 0, sizeof(double) * hess.nonZeros());
  (ES::Mp<ES::VXd>(hess.valuePtr(), hess.nonZeros())) += ES::Mp<const ES::VXd>(intg->A.valuePtr(), intg->A.nonZeros());

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

  grad -= intg->b;
}

double ImplicitBackwardEulerEnergy::func_grad_hessian(ES::ConstRefVecXd x, ES::RefVecXd grad, ES::SpMatD &hess) const
{
  // x is u (the full position)
  // A·u
  ES::mv(intg->A, x, grad, 0);

  // 0.5·uᵀ·A·u
  double mainEnergy = ES::vTMv(intg->A, x, intg->temp0, 0) * 0.5;
  double energy = mainEnergy;
  std::vector<double> componentEnergies;
  componentEnergies.reserve(intg->implicitModelsAll.size());

  hess = intg->hessianAll;
  memset(hess.valuePtr(), 0, sizeof(double) * hess.nonZeros());
  (ES::Mp<ES::VXd>(hess.valuePtr(), hess.nonZeros())) += ES::Mp<const ES::VXd>(intg->A.valuePtr(), intg->A.nonZeros());

  for (size_t i = 0; i < intg->implicitModelsAll.size(); i++) {
    ES::VXd &fint = *intg->implicitModelsAll_fint[i];
    fint.setZero();

    if (intg->implicitModelsAll[i]->isHessianTopologyFixed()) {
      // Material: energy + gradient + hessian separately
      const double modelEnergy = intg->implicitModelsAll[i]->func(x);
      energy += modelEnergy;
      componentEnergies.push_back(modelEnergy);

      ES::SpMatD &K = *intg->implicitModelsAll_K[i];
      const ES::SpMatI &mapping = *intg->implicitModelsAll_Kmaping[i];
      intg->implicitModelsAll[i]->gradient(x, fint);
      intg->implicitModelsAll[i]->hessian(x, K);
      ES::addSmallToBig(1.0, K, hess, 1.0, mapping, 1);
    }
    else {
      // IPC: combined energy + gradient + hessian (1 buildActiveSet)
      ES::SpMatD Ki;
      const double modelEnergy = intg->implicitModelsAll[i]->func_grad_hessian(x, fint, Ki);
      energy += modelEnergy;
      componentEnergies.push_back(modelEnergy);
      if (Ki.nonZeros())
        hess = hess + Ki;
    }
    grad += fint;
  }

  mainEnergy -= x.dot(intg->b);
  energy -= x.dot(intg->b);
  grad -= intg->b;
  cacheImplicitEnergyComponents(mainEnergy, componentEnergies);
  return energy;
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

void ImplicitBackwardEulerEnergy::printImplicitEnergy(ES::ConstRefVecXd x, bool allowCachedComponents) const
{
  if (allowCachedComponents && hasCachedImplicitEnergyComponents()) {
    std::cout << "  main: " << cachedMainEnergy << '\n';
    for (size_t i = 0; i < cachedImplicitModelEnergies.size(); i++) {
      std::cout << "  sub " << i << ": " << cachedImplicitModelEnergies[i] << '\n';
    }
    return;
  }

  // x is u (the full position)
  //std::cout << "Energy: ";
  // 0.5 Au^2
  double energy = ES::vTMv(intg->A, x, intg->temp0, 0) * 0.5;
  // - b^T u
  // double last_term = cblas_ddot(intg->n3, x.data(), 1, intg->b.data(), 1);
  double last_term = x.dot(intg->b);
  energy -= last_term;

  std::cout << "  main: " << energy << '\n';

  //std::cout << energy << ',';
  for (size_t i = 0; i < intg->implicitModelsAll.size(); i++) {
    // elastic(u)
    double tempEnergy = intg->implicitModelsAll[i]->func(x);

    std::cout << "  sub " << i << ": " << tempEnergy << '\n';
  }
}

void ImplicitBackwardEulerEnergy::clearCachedImplicitEnergyComponents() const
{
  hasCachedEnergyComponents = false;
  cachedImplicitModelEnergies.clear();
}

void ImplicitBackwardEulerEnergy::cacheImplicitEnergyComponents(
  double mainEnergy, const std::vector<double> &componentEnergies) const
{
  hasCachedEnergyComponents = true;
  cachedMainEnergy = mainEnergy;
  cachedImplicitModelEnergies = componentEnergies;
}

bool ImplicitBackwardEulerEnergy::hasCachedImplicitEnergyComponents() const
{
  if (!hasCachedEnergyComponents)
    return false;
  if (cachedImplicitModelEnergies.size() != intg->implicitModelsAll.size())
    return false;

  return true;
}

int ImplicitBackwardEulerEnergy::isHessianTopologyFixed() const
{
  for (size_t i = 0; i < intg->implicitModelsAll.size(); i++) {
    if (!intg->implicitModelsAll[i]->isHessianTopologyFixed())
      return 0;
  }
  return 1;
}

void ImplicitBackwardEulerEnergy::hessianDirect(ES::ConstRefVecXd x, ES::SpMatD &hess) const
{
  // Start with hessianAll pattern (covers A + all fixed-topology models)
  hess = intg->hessianAll;
  memset(hess.valuePtr(), 0, sizeof(double) * hess.nonZeros());

  // Add A (mass/damping/timestep contribution)
  (ES::Mp<ES::VXd>(hess.valuePtr(), hess.nonZeros())) += ES::Mp<const ES::VXd>(intg->A.valuePtr(), intg->A.nonZeros());

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

NonlinearOptimization::MaxStepResult ImplicitBackwardEulerEnergy::computeMaxStepLimit(ES::ConstRefVecXd x, ES::ConstRefVecXd dx) const
{
  NonlinearOptimization::MaxStepResult result = NonlinearOptimization::MaxStepResult::unconstrained();
  for (size_t i = 0; i < intg->implicitModelsAll.size(); i++) {
    const auto &model = intg->implicitModelsAll[i];
    result = NonlinearOptimization::mergeMaxStepResults(result, model->computeMaxStepLimit(x, dx));
  }
  return result;
}
