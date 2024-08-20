/*
author: Bohan Wang
copyright to USC
*/

#include "IpoptProblem.h"

#include "pgoLogging.h"

#include <coin-or/IpIpoptCalculatedQuantities.hpp>
#include <coin-or/IpIpoptData.hpp>
#include <coin-or/IpTNLPAdapter.hpp>
#include <coin-or/IpOrigIpoptNLP.hpp>

using namespace pgo::NonlinearOptimization;
namespace ES = pgo::EigenSupport;

IpoptProblem::IpoptProblem(PotentialEnergy_const_p p):
  NonlinearProblem(p)
{
  for (int i = 0; i < problem->getNumDOFs(); i++) {
    xlow[i] = -inf();
    xhi[i] = inf();
  }

  additionalEnergy = 0;
  additionalGrad = ES::VXd::Zero(xinit.size());
  additionalH = hessianAll;
  ES::Zero(additionalH);

  xcur = ES::VXd::Zero(problem->getNumDOFs());
}

IpoptProblem::IpoptProblem(PotentialEnergy_const_p p, ConstraintFunctions_const_p c):
  NonlinearProblem(p, c)
{
  for (int i = 0; i < problem->getNumDOFs(); i++) {
    xlow[i] = -inf();
    xhi[i] = inf();
  }

  additionalEnergy = 0;
  additionalGrad = ES::VXd::Zero(xinit.size());
  additionalH = hessianAll;
  ES::Zero(additionalH);

  xcur = ES::VXd::Zero(problem->getNumDOFs());
  lambdafinal = ES::VXd::Zero(constraints->getNumConstraints());
}

IpoptProblem::~IpoptProblem()
{
}

bool IpoptProblem::get_nlp_info(Ipopt::Index &n, Ipopt::Index &m,
  Ipopt::Index &nnz_jac_g, Ipopt::Index &nnz_h_lag, IndexStyleEnum &index_style)
{
  // # theta = numRotations + w vector
  n = problem->getNumDOFs();

  if (constraints) {
    m = constraints->getNumConstraints();
    nnz_jac_g = (int)jac.nonZeros();
  }
  else {
    nnz_jac_g = 0;
    m = 0;
  }

  // nonzeros in the hessian of the lagrangian
  // n * n in the hessian of the objective
  //  and 3 in the hessian of the constraints for x1)
  nnz_h_lag = 0;
  for (Eigen::Index outeri = 0; outeri < hessianAll.outerSize(); outeri++) {
    for (ES::SpMatD::InnerIterator it(hessianAll, outeri); it; ++it) {
      Eigen::Index row = it.row();
      Eigen::Index col = it.col();

      if (row >= col) {
        nnz_h_lag++;
      }
    }
  }

  // We use the standard fortran index style for row/col entries
  index_style = C_STYLE;

  return true;
}

bool IpoptProblem::get_bounds_info(Ipopt::Index n, Ipopt::Number *x_l, Ipopt::Number *x_u, Ipopt::Index m, Ipopt::Number *g_l, Ipopt::Number *g_u)
{
  for (Ipopt::Index i = 0; i < n; i++) {
    x_l[i] = xlow[i];
    x_u[i] = xhi[i];
  }

  for (Ipopt::Index i = 0; i < m; i++) {
    g_l[i] = clow[i];
    g_u[i] = chi[i];
  }

  return true;
}

bool IpoptProblem::get_starting_point(Ipopt::Index, bool init_x, Ipopt::Number *x, bool, Ipopt::Number *, Ipopt::Number *,
  Ipopt::Index m, bool init_lambda, Ipopt::Number *lambda)
{
  if (init_x) {
    memcpy(x, xinit.data(), sizeof(double) * problem->getNumDOFs());
  }

  if (init_lambda) {
    std::cerr << "Error: Require lambda but it is not set." << std::endl;
    (void)m;
    (void)lambda;
  }

  return true;
}

bool IpoptProblem::eval_f(Ipopt::Index n, const Ipopt::Number *x, bool new_x, Ipopt::Number &obj_value)
{
  // we ignore new_x
  (void)new_x;

  obj_value = problem->func(Eigen::Map<const Eigen::VectorXd>(x, n)) + additionalEnergy;

  return true;
}

bool IpoptProblem::eval_grad_f(Ipopt::Index n, const Ipopt::Number *x, bool new_x, Ipopt::Number *grad_f)
{
  memset(grad_f, 0, sizeof(double) * n);

  problem->gradient(Eigen::Map<const Eigen::VectorXd>(x, n), Eigen::Map<Eigen::VectorXd>(grad_f, n));

  Eigen::Map<Eigen::VectorXd>(grad_f, n) += additionalGrad;

  // we ignore new_x
  (void)new_x;

  return true;
}

bool IpoptProblem::eval_g(Ipopt::Index n, const Ipopt::Number *x, bool new_x, Ipopt::Index m, Ipopt::Number *g)
{
  if (!constraints)
    return true;

  constraints->func(Eigen::Map<const Eigen::VectorXd>(x, n), Eigen::Map<Eigen::VectorXd>(g, m));

  // we ignore new_x
  (void)new_x;

  return true;
}

bool IpoptProblem::eval_jac_g(Ipopt::Index n, const Ipopt::Number *x, bool new_x, Ipopt::Index m, Ipopt::Index nele_jac,
  Ipopt::Index *iRow, Ipopt::Index *jCol, Ipopt::Number *values)
{
  // we ignore new_x
  (void)new_x;

  if (!constraints)
    return true;

  PGO_ALOG((Ipopt::Index)jac.outerSize() == m);

  if (values == nullptr) {
    int idx = 0;
    for (Eigen::Index outeri = 0; outeri < jac.outerSize(); outeri++) {
      for (ES::SpMatD::InnerIterator it(jac, outeri); it; ++it) {
        Ipopt::Index row = (Ipopt::Index)it.row();
        Ipopt::Index col = (Ipopt::Index)it.col();

        iRow[idx] = row;
        jCol[idx] = col;

        idx++;
      }
    }
    PGO_ALOG(idx == nele_jac);
  }
  else {
    memset(jac.valuePtr(), 0, sizeof(double) * jac.nonZeros());
    constraints->jacobian(Eigen::Map<const Eigen::VectorXd>(x, n), jac);

    int idx = 0;
    for (Eigen::Index outeri = 0; outeri < jac.outerSize(); outeri++) {
      for (ES::SpMatD::InnerIterator it(jac, outeri); it; ++it) {
        values[idx] = it.value();
        idx++;
      }
    }

    PGO_ALOG(idx == nele_jac);
  }

  return true;
}

bool IpoptProblem::eval_h(Ipopt::Index n, const Ipopt::Number *x, bool new_x, Ipopt::Number obj_factor, Ipopt::Index m, const Ipopt::Number *lambda,
  bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index *iRow, Ipopt::Index *jCol, Ipopt::Number *values)
{
  // we ignore new_x, new_lambda
  (void)new_x;
  (void)new_lambda;

  if (values == nullptr) {
    int idx = 0;
    for (Eigen::Index outeri = 0; outeri < hessianAll.outerSize(); outeri++) {
      for (ES::SpMatD::InnerIterator it(hessianAll, outeri); it; ++it) {
        Ipopt::Index row = (Ipopt::Index)it.row();
        Ipopt::Index col = (Ipopt::Index)it.col();

        if (row >= col) {
          iRow[idx] = row;
          jCol[idx] = col;

          idx++;
        }
      }
    }

    PGO_ALOG(idx == nele_hess);
  }
  else {
    Eigen::Index nnz = hessianAll.nonZeros();
    memset(hessianAll.valuePtr(), 0, sizeof(double) * hessianAll.nonZeros());

    memset(energyHessian.valuePtr(), 0, sizeof(double) * energyHessian.nonZeros());
    problem->hessian(Eigen::Map<const Eigen::VectorXd>(x, n), energyHessian);
    ES::addSmallToBig(obj_factor, energyHessian, hessianAll, 0.0, energyHessianMapping);

    if (constraints) {
      memset(lambdaHessian.valuePtr(), 0, sizeof(double) * lambdaHessian.nonZeros());
      constraints->hessian(Eigen::Map<const Eigen::VectorXd>(x, n), Eigen::Map<const Eigen::VectorXd>(lambda, m), lambdaHessian);
      ES::addSmallToBig(1.0, lambdaHessian, hessianAll, 1.0, lambdaHessianMapping);
    }

    for (Eigen::Index i = 0; i < nnz; i++)
      hessianAll.valuePtr()[i] += additionalH.valuePtr()[i];

    int idx = 0;
    for (Eigen::Index outeri = 0; outeri < hessianAll.outerSize(); outeri++) {
      for (ES::SpMatD::InnerIterator it(hessianAll, outeri); it; ++it) {
        Ipopt::Index row = (Ipopt::Index)it.row();
        Ipopt::Index col = (Ipopt::Index)it.col();

        if (row >= col) {
          values[idx] = it.value();
          idx++;
        }
      }
    }

    PGO_ALOG(idx == nele_hess);
  }

  return true;
}

void IpoptProblem::finalize_solution(Ipopt::SolverReturn status, Ipopt::Index n, const Ipopt::Number *x,
  const Ipopt::Number *z_L, const Ipopt::Number *z_U, Ipopt::Index m,
  const Ipopt::Number *g, const Ipopt::Number *lambda,
  Ipopt::Number obj_value, const Ipopt::IpoptData *ip_data,
  Ipopt::IpoptCalculatedQuantities *ip_cq)
{
  // here is where we would store the solution to variables, or write to a file, etc
  // so we could use the solution. Since the solution is displayed to the console,
  // we currently do nothing here.
  (void)status;
  (void)z_L;
  (void)z_U;
  (void)m;
  (void)g;
  (void)obj_value;
  (void)ip_data;
  (void)ip_cq;

  if (constraints) {
    lambdafinal = Eigen::Map<const Eigen::VectorXd>(lambda, constraints->getNumConstraints());
    gfinal = Eigen::Map<const Eigen::VectorXd>(g, constraints->getNumConstraints());
  }

  xfinal = Eigen::Map<const Eigen::VectorXd>(x, n);
}

bool IpoptProblem::intermediate_callback(Ipopt::AlgorithmMode mode, Ipopt::Index iter, Ipopt::Number obj_value,
  Ipopt::Number inf_pr, Ipopt::Number inf_du, Ipopt::Number mu, Ipopt::Number d_norm, Ipopt::Number regularization_size,
  Ipopt::Number alpha_du, Ipopt::Number alpha_pr, Ipopt::Index ls_trials, const Ipopt::IpoptData *ip_data,
  Ipopt::IpoptCalculatedQuantities *ip_cq)
{
  Ipopt::TNLPAdapter *tnlp_adapter = nullptr;
  if (ip_cq != nullptr) {
    Ipopt::OrigIpoptNLP *orignlp;
    orignlp = dynamic_cast<Ipopt::OrigIpoptNLP *>(GetRawPtr(ip_cq->GetIpoptNLP()));
    if (orignlp != nullptr)
      tnlp_adapter = dynamic_cast<Ipopt::TNLPAdapter *>(GetRawPtr(orignlp->nlp()));
  }

  if (tnlp_adapter == nullptr) {
    SPDLOG_LOGGER_WARN(pgo::Logging::lgr(), "Warning: cannot get current variables.");
    return true;
  }

  tnlp_adapter->ResortX(*ip_data->curr()->x(), xcur.data());

  if (iterationCB)
    iterationCB(xcur, additionalEnergy, additionalGrad, additionalH);

  (void)mode;
  (void)iter;
  (void)obj_value;
  (void)inf_pr;
  (void)inf_du;
  (void)mu;
  (void)d_norm;
  (void)regularization_size;
  (void)alpha_du;
  (void)alpha_pr;
  (void)ls_trials;

  return true;
}
