/*
author: Bohan Wang
copyright to USC
*/

#include "knitroOptimizer.h"
#include "knitroProblem.h"

#include "EigenSupport.h"
#include "pgoLogging.h"

#include <knitro.h>

#include <mutex>
#include <thread>
#include <numeric>

#define KNITRO_ERROR(func, jump)                                                      \
  do {                                                                                \
    int error;                                                                        \
    error = func;                                                                     \
    if (error) {                                                                      \
      SPDLOG_LOGGER_ERROR(Logging::lgr(), "Knitro error happened. Error: {}", error); \
      jump;                                                                           \
    }                                                                                 \
  } while (0)

namespace pgo::NonlinearOptimization
{
class KnitroHandles
{
public:
  KnitroProblem *problem;
  KN_context *kc;
  CB_context *cb;
  int enableMultiEval = 0;

  std::mutex problemLock;
};

}  // namespace pgo::NonlinearOptimization

using namespace pgo;
using namespace pgo::NonlinearOptimization;
namespace ES = pgo::EigenSupport;

static int callbackEvalFC(KN_context_ptr /*kc*/, CB_context_ptr /*cb*/,
  KN_eval_request_ptr const evalRequest, KN_eval_result_ptr const evalResult, void *const userParams)
{
  if (evalRequest->type != KN_RC_EVALFC) {
    SPDLOG_LOGGER_ERROR(Logging::lgr(), "*** callbackEvalFC incorrectly called with eval type {}", evalRequest->type);
    return -1;
  }

  KnitroHandles *handles = reinterpret_cast<KnitroHandles *>(userParams);

  // if (handles->enableMultiEval == 0)
  handles->problemLock.lock();

  KnitroProblem *problem = handles->problem;
  Eigen::Map<const ES::VXd> xmap(evalRequest->x, problem->getn());

  if (problem->isDense()) {
    PotentialEnergyDense_const_p objFunc = problem->getDenseObjectiveFunction();
    ConstraintFunctionsDense_const_p cnstt = problem->getDenseConstraintFunctions();

    if (!objFunc->isQuadratic())
      *evalResult->obj = objFunc->func(xmap);
    else
      *evalResult->obj = 0;

    if (cnstt &&
      !cnstt->isLinear() &&
      !cnstt->isQuadratic()) {
      cnstt->func(xmap, Eigen::Map<ES::VXd>(evalResult->c, problem->getm()));
    }
  }
  else {
    PotentialEnergy_const_p objFunc = problem->getObjectiveFunction();
    ConstraintFunctions_const_p cnstt = problem->getConstraintFunctions();

    if (!objFunc->isQuadratic())
      *evalResult->obj = objFunc->func(xmap);
    else
      *evalResult->obj = 0;

    if (cnstt &&
      !cnstt->isLinear() &&
      !cnstt->isQuadratic()) {
      cnstt->func(xmap, Eigen::Map<ES::VXd>(evalResult->c, problem->getm()));
    }
  }

  // if (handles->enableMultiEval == 0)
  handles->problemLock.unlock();

  return 0;
}

static int callbackEvalGA(KN_context_ptr /*kc*/, CB_context_ptr /*cb*/,
  KN_eval_request_ptr const evalRequest, KN_eval_result_ptr const evalResult, void *const userParams)
{
  if (evalRequest->type != KN_RC_EVALGA) {
    SPDLOG_LOGGER_ERROR(Logging::lgr(), "*** callbackEvalGA incorrectly called with eval type {}", evalRequest->type);
    return -1;
  }

  KnitroHandles *handles = reinterpret_cast<KnitroHandles *>(userParams);

  // if (handles->enableMultiEval == 0)
  handles->problemLock.lock();

  KnitroProblem *problem = handles->problem;
  Eigen::Map<const ES::VXd> xmap(evalRequest->x, problem->getn());

  if (problem->isDense()) {
    PotentialEnergyDense_const_p objFunc = problem->getDenseObjectiveFunction();
    ConstraintFunctionsDense_const_p cnstt = problem->getDenseConstraintFunctions();

    if (!objFunc->isQuadratic())
      objFunc->gradient(xmap, Eigen::Map<ES::VXd>(evalResult->objGrad, problem->getn()));
    // else memset(evalResult->objGrad, 0, sizeof(double) * problem->getn());

    if (cnstt && !cnstt->isLinear() && !cnstt->isQuadratic()) {
      ES::Mp<ES::MXd> jacMap(evalResult->jac, problem->getm(), problem->getn());
      cnstt->jacobian(xmap, jacMap);
    }
  }
  else {
    PotentialEnergy_const_p objFunc = problem->getObjectiveFunction();
    ConstraintFunctions_const_p cnstt = problem->getConstraintFunctions();

    if (!objFunc->isQuadratic())
      objFunc->gradient(xmap, Eigen::Map<ES::VXd>(evalResult->objGrad, problem->getn()));
    // else memset(evalResult->objGrad, 0, sizeof(double) * problem->getn());

    if (cnstt && !cnstt->isLinear() && !cnstt->isQuadratic()) {
      ES::SpMatD &jac = problem->getJacBuffer();
      memset(jac.valuePtr(), 0, jac.nonZeros() * sizeof(double));
      cnstt->jacobian(xmap, problem->getJacBuffer());

      ES::IDX inc = 0;
      for (ES::IDX i = 0; i < jac.outerSize(); i++) {
        for (ES::SpMatD::InnerIterator it(jac, i); it; ++it)
          evalResult->jac[inc++] = it.value();
      }
    }
  }

  // if (handles->enableMultiEval == 0)
  handles->problemLock.unlock();

  return 0;
}

int callbackEvalH(KN_context_ptr /*kc*/, CB_context_ptr /*cb*/,
  KN_eval_request_ptr const evalRequest, KN_eval_result_ptr const evalResult, void *const userParams)
{
  if (evalRequest->type != KN_RC_EVALH &&
    evalRequest->type != KN_RC_EVALH_NO_F &&
    evalRequest->type != KN_RC_EVALHV &&
    evalRequest->type != KN_RC_EVALHV_NO_F) {
    printf("*** callbackEvalH incorrectly called with eval type %d\n",
      evalRequest->type);
    return (-1);
  }

  KnitroHandles *handles = reinterpret_cast<KnitroHandles *>(userParams);

  // if (handles->enableMultiEval == 0)
  handles->problemLock.lock();

  KnitroProblem *problem = handles->problem;
  Eigen::Map<const ES::VXd> xmap(evalRequest->x, problem->getn());
  ES::VXd emptyLambda;

  // for (int i = 0; i < problem->getm(); i++)
  //   LG_ << evalRequest->lambda[i] << ',';
  // LG_ << '\n';

  int retCode = 0;
  switch (evalRequest->type) {
  case KN_RC_EVALH: {
    double sigma = *evalRequest->sigma;
    if (problem->isDense()) {
      if (problem->getDenseObjectiveFunction()->isQuadratic())
        sigma = 0;

      ES::MXd &h = problem->hessianDense(xmap,
        problem->getm() > 0 ? Eigen::Map<const ES::VXd>(evalRequest->lambda, problem->getm()) : emptyLambda,
        sigma);

      int inc = 0;
      for (ES::IDX col = 0; col < h.cols(); col++) {
        for (ES::IDX row = 0; row <= col; row++) {
          evalResult->hess[inc++] = h(row, col);
        }
      }
    }
    else {
      if (problem->getObjectiveFunction()->isQuadratic())
        sigma = 0;

      ES::SpMatD &h = problem->hessian(xmap,
        problem->getm() > 0 ? Eigen::Map<const ES::VXd>(evalRequest->lambda, problem->getm()) : emptyLambda,
        sigma);

      for (ES::IDX i = 0, inc = 0; i < h.outerSize(); i++) {
        for (ES::SpMatD::InnerIterator it(h, i); it; ++it) {
          if (it.row() > it.col())
            continue;

          evalResult->hess[inc++] = it.value();
        }
      }
    }

    retCode = 0;
    break;
  }
  case KN_RC_EVALH_NO_F: {
    if (problem->isDense()) {
      ES::MXd &h = problem->hessianDense(xmap,
        problem->getm() > 0 ? Eigen::Map<const ES::VXd>(evalRequest->lambda, problem->getm()) : emptyLambda,
        0.0);

      int inc = 0;
      for (ES::IDX col = 0; col < h.cols(); col++) {
        for (ES::IDX row = 0; row <= col; row++) {
          evalResult->hess[inc++] = h(row, col);
        }
      }
    }
    else {
      ES::SpMatD &h = problem->hessian(xmap,
        problem->getm() > 0 ? Eigen::Map<const ES::VXd>(evalRequest->lambda, problem->getm()) : emptyLambda,
        0.0);

      for (ES::IDX i = 0, inc = 0; i < h.outerSize(); i++) {
        for (ES::SpMatD::InnerIterator it(h, i); it; ++it) {
          if (it.row() > it.col())
            continue;

          evalResult->hess[inc++] = it.value();
        }
      }
    }
    retCode = 0;
    break;
  }
  case KN_RC_EVALHV: {
    // LGI << "Hessian vector routine";
    PGO_ALOG(problem->isDense() == false);

    double sigma = *evalRequest->sigma;
    if (problem->getObjectiveFunction()->isQuadratic())
      sigma = 0;

    ES::VXd &hVec = problem->hessianVector(xmap,
      problem->getm() > 0 ? Eigen::Map<const ES::VXd>(evalRequest->lambda, problem->getm()) : emptyLambda,
      Eigen::Map<const ES::VXd>(evalRequest->vec, problem->getn()),
      sigma);

    memcpy(evalResult->hessVec, hVec.data(), sizeof(double) * problem->getn());

    retCode = 0;

    break;
  }
  case KN_RC_EVALHV_NO_F: {
    // LGI << "Hessian vector routine";
    PGO_ALOG(problem->isDense() == false);

    ES::VXd &hVec = problem->hessianVector(xmap,
      problem->getm() > 0 ? Eigen::Map<const ES::VXd>(evalRequest->lambda, problem->getm()) : emptyLambda,
      Eigen::Map<const ES::VXd>(evalRequest->vec, problem->getn()),
      0.0);

    memcpy(evalResult->hessVec, hVec.data(), sizeof(double) * problem->getn());

    // LGE << "Not supported yet";
    retCode = 0;

    break;
  }

  default:
    SPDLOG_LOGGER_ERROR(Logging::lgr(), "*** callbackEvalHess incorrectly called with eval code {}", evalRequest->type);
    retCode = 1;
    break;
  }

  // if (handles->enableMultiEval == 0)
  handles->problemLock.unlock();

  return retCode;
}

static int callbackNewPt(KN_context_ptr /*kc*/,
  const double *const x, const double *const lambda, void *const userParams)
{
  KnitroHandles *handles = reinterpret_cast<KnitroHandles *>(userParams);
  KnitroProblem *problem = handles->problem;

  if (problem->getConstraintFunctions() || problem->getDenseConstraintFunctions()) {
    problem->setCurrentSolution(Eigen::Map<const ES::VXd>(x, problem->getn()), Eigen::Map<const ES::VXd>(lambda, problem->getm()));
  }
  else {
    problem->setCurrentSolution(Eigen::Map<const ES::VXd>(x, problem->getn()), ES::VXd());
  }

  problem->iterationCallback();

  return 0;
}

KnitroOptimizer::KnitroOptimizer(KnitroProblem *p)
{
  handles = new KnitroHandles;
  handles->problem = p;

  int error = KN_new(&handles->kc);
  if (error) {
    SPDLOG_LOGGER_ERROR(Logging::lgr(), "Error in create Knitro context.");
    exit(1);
  }

  KNITRO_ERROR(KN_set_int_param(handles->kc, KN_PARAM_ALGORITHM, 1),
    throw std::domain_error("set algorithm"));

  x = new double[p->getn()];
  lambda = new double[p->getn() + p->getm()];
  g = new double[p->getm()];
}

KnitroOptimizer::~KnitroOptimizer()
{
  KN_free(&handles->kc);

  delete handles;
  delete[] x;
  delete[] lambda;
  delete[] g;
}

void KnitroOptimizer::setConfigFile(const char *filename)
{
  if (filename && strlen(filename)) {
    KNITRO_ERROR(KN_load_param_file(handles->kc, filename),
      throw std::domain_error("set config file"));

    KNITRO_ERROR(KN_set_int_param(handles->kc, KN_PARAM_PAR_NUMTHREADS, (int)std::thread::hardware_concurrency() / 2),
      throw std::domain_error("set #threads"));
  }
  else {
    SPDLOG_LOGGER_WARN(Logging::lgr(), "No config file specified. Use default settings.");

    KNITRO_ERROR(KN_set_int_param(handles->kc, KN_PARAM_ALGORITHM, KN_ALG_BAR_DIRECT),
      throw std::domain_error("set algorithm"));

    KNITRO_ERROR(KN_set_int_param(handles->kc, KN_PARAM_HESSOPT, KN_HESSOPT_EXACT),
      throw std::domain_error("set hessopt"));

    KNITRO_ERROR(KN_set_int_param(handles->kc, KN_PARAM_PAR_CONCURRENT_EVALS, KN_PAR_CONCURRENT_EVALS_NO),
      throw std::domain_error("set par_concurrent_evals"));

    KNITRO_ERROR(KN_set_int_param(handles->kc, KN_PARAM_HESSIAN_NO_F, KN_HESSIAN_NO_F_ALLOW),
      throw std::domain_error("set hessian_no_f"));

    KNITRO_ERROR(KN_set_int_param(handles->kc, KN_PARAM_OUTLEV, KN_OUTLEV_ITER_VERBOSE),
      throw std::domain_error("set outlev"));

    setMaxIter(200);
    setFeasTol(1e-8);
    setOptTol(1e-8);

    KNITRO_ERROR(KN_set_int_param(handles->kc, KN_PARAM_DERIVCHECK_TYPE, KN_DERIVCHECK_CENTRAL),
      throw std::domain_error("set derivcheck_type"));

    KNITRO_ERROR(KN_set_int_param(handles->kc, KN_PARAM_DERIVCHECK, KN_DERIVCHECK_NONE),
      throw std::domain_error("set derivcheck"));

    KNITRO_ERROR(KN_set_double_param(handles->kc, KN_PARAM_DERIVCHECK_TOL, 1e-6),
      throw std::domain_error("set derivcheck_tol"));

    KNITRO_ERROR(KN_set_double_param(handles->kc, KN_PARAM_BNDRANGE, 1e20),
      throw std::domain_error("set bndrange"));

    KNITRO_ERROR(KN_set_int_param(handles->kc, KN_PARAM_HONORBNDS, KN_HONORBNDS_ALWAYS),
      throw std::domain_error("set honorbnds"));

    KNITRO_ERROR(KN_set_int_param(handles->kc, KN_PARAM_DATACHECK, KN_DATACHECK_NO),
      throw std::domain_error("set datacheck"));
  }
}

void KnitroOptimizer::setMaxIter(int maxIter)
{
  KNITRO_ERROR(KN_set_int_param(handles->kc, KN_PARAM_MAXIT, maxIter),
    throw std::domain_error("set max iter"));
}

void KnitroOptimizer::setFeasTol(double eps)
{
  KNITRO_ERROR(KN_set_double_param(handles->kc, KN_PARAM_FEASTOL, eps),
    throw std::domain_error("set feasible tol"));
}

void KnitroOptimizer::setOptTol(double eps)
{
  KNITRO_ERROR(KN_set_double_param(handles->kc, KN_PARAM_OPTTOL, eps),
    throw std::domain_error("set optimality tol"));
}

void KnitroOptimizer::setVerbose(int verbose)
{
  KNITRO_ERROR(KN_set_int_param(handles->kc, KN_PARAM_OUTLEV, verbose),
    throw std::domain_error("set verbose level"));
}

void KnitroOptimizer::enableWarmStart(int enable)
{
  KNITRO_ERROR(KN_set_int_param(handles->kc, KN_PARAM_STRAT_WARM_START, enable ? 1 : 0),
    throw std::domain_error("set warm start"));
}

void KnitroOptimizer::enableMultiEvaluation(int enable)
{
  handles->enableMultiEval = enable ? 1 : 0;
  KNITRO_ERROR(KN_set_int_param(handles->kc, KN_PARAM_PAR_CONCURRENT_EVALS, handles->enableMultiEval),
    throw std::domain_error("set multi eval"));
}
void KnitroOptimizer::honorBoundary(int enable)
{
  KNITRO_ERROR(KN_set_int_param(handles->kc, KN_PARAM_HONORBNDS, enable),
    throw std::domain_error("set honor boundary"));
}

void KnitroOptimizer::emphasisFeasibility(int opt)
{
  if (opt == 0) {
    KNITRO_ERROR(KN_set_int_param(handles->kc, KN_PARAM_BAR_FEASIBLE, KN_BAR_FEASIBLE_NO),
      throw std::domain_error("set bar feasible"));
  }
  else if (opt == 1) {
    KNITRO_ERROR(KN_set_int_param(handles->kc, KN_PARAM_BAR_FEASIBLE, KN_BAR_FEASIBLE_STAY),
      throw std::domain_error("set bar feasible"));
  }
  else if (opt == 2) {
    KNITRO_ERROR(KN_set_int_param(handles->kc, KN_PARAM_BAR_FEASIBLE, KN_BAR_FEASIBLE_GET),
      throw std::domain_error("set bar feasible"));
  }
  else if (opt == 3) {
    KNITRO_ERROR(KN_set_int_param(handles->kc, KN_PARAM_BAR_FEASIBLE, KN_BAR_FEASIBLE_GET_STAY),
      throw std::domain_error("set bar feasible"));
  }
  else {
    throw std::domain_error("unknown opt");
  }
}

void KnitroOptimizer::initQuadraticProblem()
{
  std::vector<KNINT> rows;
  std::vector<KNINT> cols;
  std::vector<double> values;

  ES::SpMatD h;
  handles->problem->getObjectiveFunction()->createHessian(h);

  ES::VXd x = ES::VXd::Zero(h.rows()), grad = ES::VXd::Zero(h.rows());
  handles->problem->getObjectiveFunction()->hessian(x, h);
  handles->problem->getObjectiveFunction()->gradient(x, grad);

  for (ES::IDX i = 0; i < h.outerSize(); i++) {
    for (ES::SpMatD::InnerIterator it(h, i); it; ++it) {
      rows.push_back((KNINT)it.row());
      cols.push_back((KNINT)it.col());
      values.push_back(it.value() * 0.5);
    }
  }

  KNITRO_ERROR(KN_add_obj_quadratic_struct(handles->kc, (KNLONG)rows.size(), rows.data(), cols.data(), values.data()),
    throw std::domain_error("set quadratic components"));

  std::vector<int> dofs;
  std::vector<double> linearCoeffs;
  for (ES::IDX i = 0; i < grad.size(); i++) {
    if (fabs(grad[i]) > 0) {
      dofs.push_back((int)i);
      linearCoeffs.push_back(grad[i]);
    }
  }

  if (dofs.size()) {
    KNITRO_ERROR(KN_add_obj_linear_struct(handles->kc, (KNINT)dofs.size(), dofs.data(), linearCoeffs.data()),
      throw std::domain_error("set linear components"));
  }
}

void KnitroOptimizer::init()
{
  KNITRO_ERROR(KN_set_int_param(handles->kc, KN_PARAM_PAR_NUMTHREADS, (int)std::thread::hardware_concurrency()),
    throw std::domain_error("set num threads"));

  /** Initialize Knitro with the problem definition. */
  KNITRO_ERROR(KN_add_vars(handles->kc, handles->problem->getn(), nullptr),
    throw std::domain_error("add variable"));

  KNITRO_ERROR(KN_set_var_lobnds_all(handles->kc, handles->problem->getXLow().data()),
    throw std::domain_error("set x lower bound"));

  KNITRO_ERROR(KN_set_var_upbnds_all(handles->kc, handles->problem->getXHi().data()),
    throw std::domain_error("set x hi bound"));

  /** Define an initial point.  If not set, Knitro will generate one. */
  KNITRO_ERROR(KN_set_var_primal_init_values_all(handles->kc, handles->problem->getXInit().data()),
    throw std::domain_error("set x init"));

  /** Add the constraints and set their lower bounds */
  if (handles->problem->getm() > 0) {
    KNITRO_ERROR(KN_add_cons(handles->kc, handles->problem->getm(), nullptr),
      throw std::domain_error("set # constraints"));

    KNITRO_ERROR(KN_set_con_lobnds_all(handles->kc, handles->problem->getCLow().data()),
      throw std::domain_error("set constraints low bound"));

    KNITRO_ERROR(KN_set_con_upbnds_all(handles->kc, handles->problem->getCHi().data()),
      throw std::domain_error("set constraints hi bound"));

    if (handles->problem->isDense()) {
      PotentialEnergyDense_const_p objFunc = handles->problem->getDenseObjectiveFunction();
      ConstraintFunctionsDense_const_p cnstt = handles->problem->getDenseConstraintFunctions();
      // get indices of the hessian
      std::vector<KNINT> constraintsIndices;

      // if the constraints are linear,
      // we directly use their built-in routine
      if (cnstt->isLinear()) {
        std::vector<std::tuple<int, int, int, double>> qTerms;
        std::vector<std::tuple<int, int, double>> lTerms;
        ES::VXd c(handles->problem->getm());

        cnstt->getLinearAndQuadraticTerms(qTerms, lTerms, c);

        std::vector<KNINT> cIndices, lIndices;
        std::vector<double> coeffs;

        cIndices.reserve(lTerms.size());
        lIndices.reserve(lTerms.size());
        coeffs.reserve(lTerms.size());
        for (const auto &tp : lTerms) {
          cIndices.emplace_back(std::get<0>(tp));
          lIndices.emplace_back(std::get<1>(tp));
          coeffs.emplace_back(std::get<2>(tp));
        }

        KNITRO_ERROR(KN_add_con_linear_struct(handles->kc, (KNLONG)lTerms.size(), cIndices.data(), lIndices.data(), coeffs.data()),
          throw std::domain_error("set linear constraints"));

        ES::VXd clow = handles->problem->getCLow() - c;
        ES::VXd chi = handles->problem->getCHi() - c;

        KNITRO_ERROR(KN_set_con_lobnds_all(handles->kc, clow.data()),
          throw std::domain_error("set constraints low bound"));

        KNITRO_ERROR(KN_set_con_upbnds_all(handles->kc, chi.data()),
          throw std::domain_error("set constraints hi bound"));
      }
      // if the constraints are quadratic,
      // we directly use their built-in routine
      // untested
      else if (cnstt->isQuadratic()) {
        std::vector<std::tuple<int, int, int, double>> qTerms;
        std::vector<std::tuple<int, int, double>> lTerms;
        ES::VXd c(handles->problem->getm());

        cnstt->getLinearAndQuadraticTerms(qTerms, lTerms, c);

        std::vector<KNINT> qcIndices, qIndices1, qIndices2;
        std::vector<KNINT> lcIndices, lIndices;
        std::vector<double> qcoeffs, lcoeffs;

        qcIndices.reserve(qTerms.size());
        qIndices1.reserve(qTerms.size());
        qIndices2.reserve(qTerms.size());
        qcoeffs.reserve(qTerms.size());
        for (const auto &tp : qTerms) {
          qcIndices.emplace_back(std::get<0>(tp));
          qIndices1.emplace_back(std::get<1>(tp));
          qIndices2.emplace_back(std::get<2>(tp));
          qcoeffs.emplace_back(std::get<3>(tp));
        }

        lcIndices.reserve(lTerms.size());
        lIndices.reserve(lTerms.size());
        lcoeffs.reserve(lTerms.size());
        for (const auto &tp : lTerms) {
          lcIndices.emplace_back(std::get<0>(tp));
          lIndices.emplace_back(std::get<1>(tp));
          lcoeffs.emplace_back(std::get<2>(tp));
        }

        KNITRO_ERROR(KN_add_con_quadratic_struct(handles->kc, (KNLONG)qTerms.size(), qcIndices.data(), qIndices1.data(), qIndices2.data(), qcoeffs.data()),
          throw std::domain_error("set quadratic constraints"));

        KNITRO_ERROR(KN_add_con_linear_struct(handles->kc, (KNLONG)lTerms.size(), lcIndices.data(), lIndices.data(), lcoeffs.data()),
          throw std::domain_error("set linear constraints"));

        ES::VXd clow = handles->problem->getCLow() - c;
        ES::VXd chi = handles->problem->getCHi() - c;

        KNITRO_ERROR(KN_set_con_lobnds_all(handles->kc, clow.data()),
          throw std::domain_error("set constraints low bound"));

        KNITRO_ERROR(KN_set_con_upbnds_all(handles->kc, chi.data()),
          throw std::domain_error("set constraints hi bound"));
      }
      // if the constraints are general
      else {
        for (int i = 0; i < handles->problem->getm(); i++) {
          constraintsIndices.push_back(i);
        }
      }

      // if it is a general constraints
      if (constraintsIndices.size()) {
        // if the objective is quadratic
        if (objFunc->isQuadratic()) {
          initQuadraticProblem();

          KNITRO_ERROR(KN_add_eval_callback(handles->kc, KNFALSE,
                         (KNINT)constraintsIndices.size(), constraintsIndices.data(),
                         callbackEvalFC, &handles->cb),
            throw std::domain_error("set fc callback"));

          KNITRO_ERROR(KN_set_cb_grad(handles->kc, handles->cb, 0, nullptr,
                         KN_DENSE_COLMAJOR, nullptr, nullptr,
                         callbackEvalGA),
            throw std::domain_error("set ga callback"));
        }
        // otherwise we address it in the general way
        else {
          KNITRO_ERROR(KN_add_eval_callback(handles->kc, KNTRUE,
                         (KNINT)constraintsIndices.size(), constraintsIndices.data(),
                         callbackEvalFC, &handles->cb),
            throw std::domain_error("set fc callback"));

          KNITRO_ERROR(KN_set_cb_grad(handles->kc, handles->cb, KN_DENSE, nullptr,
                         KN_DENSE_COLMAJOR, nullptr, nullptr,
                         callbackEvalGA),
            throw std::domain_error("set ga callback"));
        }
      }
      // if there is no general constraints
      else {
        // if the objective is quadratic
        if (objFunc->isQuadratic()) {
          initQuadraticProblem();

          KNITRO_ERROR(KN_add_eval_callback(handles->kc, KNFALSE, 0, nullptr,
                         callbackEvalFC, &handles->cb),
            throw std::domain_error("set fc callback"));

          KNITRO_ERROR(KN_set_cb_grad(handles->kc, handles->cb, 0, nullptr, 0, nullptr, nullptr,
                         callbackEvalGA),
            throw std::domain_error("set ga callback"));
        }
        else {
          KNITRO_ERROR(KN_add_eval_callback(handles->kc, KNTRUE, 0, nullptr,
                         callbackEvalFC, &handles->cb),
            throw std::domain_error("set fc callback"));

          KNITRO_ERROR(KN_set_cb_grad(handles->kc, handles->cb, KN_DENSE, nullptr, 0, nullptr, nullptr,
                         callbackEvalGA),
            throw std::domain_error("set ga callback"));
        }
      }
    }
    // if it is a sparse problem
    else {
      // get indices of the hessian
      std::vector<KNINT> constraintsIndices;
      ES::SpMatD &jac = handles->problem->getJacBuffer();
      handles->problem->getConstraintFunctions()->jacobian(handles->problem->getXInit(), jac);

      std::vector<KNINT> jacRows(jac.nonZeros());
      std::vector<KNINT> jacCols(jac.nonZeros());
      std::vector<double> jacCoeffs(jac.nonZeros());

      ES::IDX inc = 0;
      for (ES::IDX r = 0; r < jac.outerSize(); r++) {
        for (ES::SpMatD::InnerIterator it(jac, r); it; ++it) {
          jacRows[inc] = (KNINT)it.row();
          jacCols[inc] = (KNINT)it.col();
          jacCoeffs[inc] = it.value();
          inc++;
        }
      }

      // if the constraints are linear,
      // we directly use their built-in routine
      if (handles->problem->getConstraintFunctions()->isLinear()) {
        KNITRO_ERROR(KN_add_con_linear_struct(handles->kc, (KNLONG)jac.nonZeros(), jacRows.data(), jacCols.data(), jacCoeffs.data()),
          throw std::domain_error("set linear constraints"));

        ES::VXd zero = ES::VXd::Zero(handles->problem->getn());
        ES::VXd d = ES::VXd::Zero(handles->problem->getm());
        handles->problem->getConstraintFunctions()->func(zero, d);

        ES::VXd clow = handles->problem->getCLow() - d;
        ES::VXd chi = handles->problem->getCHi() - d;

        KNITRO_ERROR(KN_set_con_lobnds_all(handles->kc, clow.data()),
          throw std::domain_error("set constraints low bound"));

        KNITRO_ERROR(KN_set_con_upbnds_all(handles->kc, chi.data()),
          throw std::domain_error("set constraints hi bound"));
      }
      // if the constraints are quadratic,
      // we directly use their built-in routine
      // untested
      else if (handles->problem->getConstraintFunctions()->isQuadratic()) {
        ES::SpMatD &lambdaHessian = handles->problem->getLambdaHessianBuffer();
        ES::VXd lambda(handles->problem->getm());

        std::vector<KNINT> rows;
        std::vector<KNINT> var0, var1;
        std::vector<double> coeffs;

        rows.reserve(lambdaHessian.nonZeros() * handles->problem->getm());
        var0.reserve(lambdaHessian.nonZeros() * handles->problem->getm());
        var1.reserve(lambdaHessian.nonZeros() * handles->problem->getm());
        coeffs.reserve(lambdaHessian.nonZeros() * handles->problem->getm());

        for (int i = 0; i < handles->problem->getm(); i++) {
          lambda.noalias() = ES::VXd::Zero(handles->problem->getm());
          lambda[i] = 1.0;

          memset(lambdaHessian.valuePtr(), 0, lambdaHessian.nonZeros() * sizeof(double));
          handles->problem->getConstraintFunctions()->hessian(handles->problem->getXInit(), lambda, lambdaHessian);

          for (ES::IDX r = 0; r < lambdaHessian.outerSize(); r++) {
            for (ES::SpMatD::InnerIterator it(lambdaHessian, r); it; ++it) {
              if (fabs(it.value()) > 0) {
                rows.push_back((KNINT)i);
                var0.push_back((KNINT)it.row());
                var1.push_back((KNINT)it.col());
                coeffs.push_back(it.value());
              }
            }
          }
        }

        KNITRO_ERROR(KN_add_con_quadratic_struct(handles->kc, (KNLONG)rows.size(), rows.data(), var0.data(), var1.data(), coeffs.data()),
          throw std::domain_error("set quadratic constraints"));
      }
      // if the constraints are general
      else {
        for (int i = 0; i < handles->problem->getm(); i++) {
          constraintsIndices.push_back(i);
        }
      }

      // if it is a general constraints
      if (constraintsIndices.size()) {
        // if the objective is quadratic
        if (handles->problem->getObjectiveFunction()->isQuadratic()) {
          initQuadraticProblem();

          KNITRO_ERROR(KN_add_eval_callback(handles->kc, KNFALSE,
                         (KNINT)constraintsIndices.size(), constraintsIndices.data(),
                         callbackEvalFC, &handles->cb),
            throw std::domain_error("set fc callback"));

          KNITRO_ERROR(KN_set_cb_grad(handles->kc, handles->cb, 0, nullptr,
                         (KNLONG)jacRows.size(), jacRows.data(), jacCols.data(),
                         callbackEvalGA),
            throw std::domain_error("set ga callback"));
        }
        // otherwise we address it in the general way
        else {
          KNITRO_ERROR(KN_add_eval_callback(handles->kc, KNTRUE,
                         (KNINT)constraintsIndices.size(), constraintsIndices.data(),
                         callbackEvalFC, &handles->cb),
            throw std::domain_error("set fc callback"));

          KNITRO_ERROR(KN_set_cb_grad(handles->kc, handles->cb, KN_DENSE, nullptr,
                         (KNLONG)jacRows.size(), jacRows.data(), jacCols.data(),
                         callbackEvalGA),
            throw std::domain_error("set ga callback"));
        }
      }
      // if there is no general constraints
      else {
        // if the objective is quadratic
        if (handles->problem->getObjectiveFunction()->isQuadratic()) {
          initQuadraticProblem();

          KNITRO_ERROR(KN_add_eval_callback(handles->kc, KNFALSE, 0, nullptr,
                         callbackEvalFC, &handles->cb),
            throw std::domain_error("set fc callback"));

          KNITRO_ERROR(KN_set_cb_grad(handles->kc, handles->cb, 0, nullptr, 0, nullptr, nullptr,
                         callbackEvalGA),
            throw std::domain_error("set ga callback"));
        }
        else {
          KNITRO_ERROR(KN_add_eval_callback(handles->kc, KNTRUE, 0, nullptr,
                         callbackEvalFC, &handles->cb),
            throw std::domain_error("set fc callback"));

          KNITRO_ERROR(KN_set_cb_grad(handles->kc, handles->cb, KN_DENSE, nullptr, 0, nullptr, nullptr,
                         callbackEvalGA),
            throw std::domain_error("set ga callback"));
        }
      }
    }
  }
  else {
    if (handles->problem->isDense()) {
      if (handles->problem->getDenseObjectiveFunction()->isQuadratic()) {
        initQuadraticProblem();

        KNITRO_ERROR(KN_add_eval_callback(handles->kc, KNFALSE, 0, nullptr,
                       callbackEvalFC, &handles->cb),
          throw std::domain_error("set fc callback"));
      }
      else {
        KNITRO_ERROR(KN_add_eval_callback(handles->kc, KNTRUE, 0, nullptr,
                       callbackEvalFC, &handles->cb),
          throw std::domain_error("set fc callback"));

        KNITRO_ERROR(KN_set_cb_grad(handles->kc, handles->cb, KN_DENSE, nullptr, 0, nullptr, nullptr,
                       callbackEvalGA),
          throw std::domain_error("set ga callback"));
      }
    }
    else {
      if (handles->problem->getObjectiveFunction()->isQuadratic()) {
        initQuadraticProblem();

        KNITRO_ERROR(KN_add_eval_callback(handles->kc, KNFALSE, 0, nullptr,
                       callbackEvalFC, &handles->cb),
          throw std::domain_error("set fc callback"));
      }
      else {
        KNITRO_ERROR(KN_add_eval_callback(handles->kc, KNTRUE, 0, nullptr,
                       callbackEvalFC, &handles->cb),
          throw std::domain_error("set fc callback"));

        KNITRO_ERROR(KN_set_cb_grad(handles->kc, handles->cb, KN_DENSE, nullptr, 0, nullptr, nullptr,
                       callbackEvalGA),
          throw std::domain_error("set ga callback"));
      }
    }
  }

  KNITRO_ERROR(KN_set_cb_user_params(handles->kc, handles->cb, handles),
    throw std::domain_error("set user params in cb"));

  KNITRO_ERROR(KN_set_newpt_callback(handles->kc, callbackNewPt, handles),
    throw std::domain_error("set new pt callback"));

  if (handles->problem->isDense()) {
    KNITRO_ERROR(KN_set_cb_hess(handles->kc, handles->cb, KN_DENSE_COLMAJOR, nullptr, nullptr, callbackEvalH),
      throw std::domain_error("set hessian callback"));

    KNITRO_ERROR(KN_set_int_param(handles->kc, KN_PARAM_LINSOLVER, 0),
      throw std::domain_error("set linear solver"));
  }
  else {
    std::vector<KNINT> hessianRows(handles->problem->getFinalHessianBuffer().nonZeros());
    std::vector<KNINT> hessianCols(handles->problem->getFinalHessianBuffer().nonZeros());

    ES::SpMatD &h = handles->problem->getFinalHessianBuffer();
    ES::IDX nonZeros = 0;
    for (ES::IDX i = 0; i < h.outerSize(); i++) {
      for (ES::SpMatD::InnerIterator it(h, i); it; ++it) {
        if (it.row() > it.col())
          continue;

        hessianRows[nonZeros] = (KNINT)it.row();
        hessianCols[nonZeros] = (KNINT)it.col();
        nonZeros++;
      }
    }

    KNITRO_ERROR(KN_set_cb_hess(handles->kc, handles->cb, nonZeros, hessianRows.data(), hessianCols.data(), callbackEvalH),
      throw std::domain_error("set hessian callback"));
  }

  /** Specify that the user is able to provide evaluations
   *  of the hessian matrix without the objective component.
   *  turned off by default but should be enabled if possible. */
  KNITRO_ERROR(KN_set_int_param(handles->kc, KN_PARAM_HESSIAN_NO_F, KN_HESSIAN_NO_F_ALLOW),
    throw std::domain_error("allow hess without f"));

  /** Set minimize or maximize (if not set, assumed minimize) */
  KNITRO_ERROR(KN_set_obj_goal(handles->kc, KN_OBJGOAL_MINIMIZE),
    throw std::domain_error("goal of the obj is minimize"));
}

int KnitroOptimizer::solve()
{
  /** Solve the problem.
   *
   *  Return status codes are defined in "knitro.h" and described
   *  in the Knitro manual.
   */
  int retCode = KN_solve(handles->kc);
  double energy = 0;
  KNITRO_ERROR(KN_get_solution(handles->kc, &retCode, &energy, x, lambda),
    throw std::domain_error("get solution"));

  if (handles->problem->hasConstraints()) {
    KNITRO_ERROR(KN_get_con_values_all(handles->kc, g),
      throw std::domain_error("get constraints"));
  }

  KNITRO_ERROR(KN_get_abs_feas_error(handles->kc, &feasError),
    throw std::domain_error("get feasibility error"));

  KNITRO_ERROR(KN_get_abs_opt_error(handles->kc, &optError),
    throw std::domain_error("get optimality error"));

  return retCode;
}

void KnitroOptimizer::setxInit(const double *x)
{
  /** Define an initial point.  If not set, Knitro will generate one. */
  KNITRO_ERROR(KN_set_var_primal_init_values_all(handles->kc, x),
    throw std::domain_error("set x init"));
}

void KnitroOptimizer::setxRange(const double *xlow, const double *xhi)
{
  KNITRO_ERROR(KN_set_var_lobnds_all(handles->kc, xlow),
    throw std::domain_error("set x lower bound"));

  KNITRO_ERROR(KN_set_var_upbnds_all(handles->kc, xhi),
    throw std::domain_error("set x hi bound"));
}

void KnitroOptimizer::setcRange(const double *lo, const double *hi)
{
  if (!handles->problem->hasConstraints())
    return;

  ES::VXd clow, chi;
  if (handles->problem->isDense()) {
    if (handles->problem->getDenseConstraintFunctions()->isLinear() ||
      handles->problem->getDenseConstraintFunctions()->isQuadratic()) {
      ES::VXd zero = ES::VXd::Zero(handles->problem->getn());
      ES::VXd d = ES::VXd::Zero(handles->problem->getm());
      handles->problem->getDenseConstraintFunctions()->func(zero, d);

      clow = Eigen::Map<const ES::VXd>(lo, d.size()) - d;
      chi = Eigen::Map<const ES::VXd>(hi, d.size()) - d;
    }
    else {
      clow = Eigen::Map<const ES::VXd>(lo, handles->problem->getm());
      chi = Eigen::Map<const ES::VXd>(hi, handles->problem->getm());
    }
  }
  else {
    if (handles->problem->getConstraintFunctions()->isLinear() ||
      handles->problem->getConstraintFunctions()->isQuadratic()) {
      ES::VXd zero = ES::VXd::Zero(handles->problem->getn());
      ES::VXd d = ES::VXd::Zero(handles->problem->getm());
      handles->problem->getConstraintFunctions()->func(zero, d);

      clow = Eigen::Map<const ES::VXd>(lo, d.size()) - d;
      chi = Eigen::Map<const ES::VXd>(hi, d.size()) - d;
    }
    else {
      clow = Eigen::Map<const ES::VXd>(lo, handles->problem->getm());
      chi = Eigen::Map<const ES::VXd>(hi, handles->problem->getm());
    }
  }

  KNITRO_ERROR(KN_set_con_lobnds_all(handles->kc, clow.data()),
    throw std::domain_error("set constraints low bound"));

  KNITRO_ERROR(KN_set_con_upbnds_all(handles->kc, chi.data()),
    throw std::domain_error("set constraints hi bound"));
}

void KnitroOptimizer::printInfo() const
{
  int niter = 0;
  KNITRO_ERROR(KN_get_number_iters(handles->kc, &niter),
    throw std::domain_error("get number iters"));

  double timeCost = 0;
  KNITRO_ERROR(KN_get_solve_time_real(handles->kc, &timeCost),
    throw std::domain_error("get solve time real"));

  SPDLOG_LOGGER_INFO(Logging::lgr(), "Solver #iterations: {}; Time cost: {}s", niter, timeCost);
}