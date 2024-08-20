#include "timeIntegratorSolver.h"

#if defined(PGO_HAS_KNITRO)
#  include "knitroOptimizer.h"
#  include "knitroProblem.h"
#endif

#if defined(PGO_HAS_IPOPT)
#  include "IpoptOptimizer.h"
#  include "IpoptProblem.h"
#endif

#include "NewtonRaphsonSolver.h"

#include "potentialEnergies.h"
#include "constraintFunctions.h"
#include "lagrangian.h"
#include "minimizeEnergy.h"

using namespace pgo;
using namespace pgo::NonlinearOptimization;
using namespace pgo::Simulation;

namespace ES = EigenSupport;

namespace pgo::Simulation
{
class TimeIntegratorSolverData
{
public:
#if defined(PGO_HAS_KNITRO)
  std::shared_ptr<KnitroOptimizer> opt;
  std::shared_ptr<KnitroProblem> problem;
#endif

  std::shared_ptr<NewtonRaphsonSolver> newtonSolver;
  std::vector<int> fixedDOFs;
  ES::VXd fixedValues;
};

}  // namespace pgo::Simulation

TimeIntegratorSolver::TimeIntegratorSolver()
{
  da = std::make_shared<TimeIntegratorSolverData>();
}

int TimeIntegratorSolver::solve(bool needRenew, ES::VXd &x,
  [[maybe_unused]] ES::VXd &g, [[maybe_unused]] ES::VXd &lambda,
  const ES::VXd &xlow, const ES::VXd &xhi,
  [[maybe_unused]] const ES::VXd &clow, [[maybe_unused]] const ES::VXd &chi,
  std::shared_ptr<const PotentialEnergy> energy,
  [[maybe_unused]] std::shared_ptr<const ConstraintFunctions> constraints,
  int nIter, double eps, int verbose, [[maybe_unused]] const char *solverConfigFilename,
  TimeIntegratorSolverOption op)
{
  if (op == TimeIntegratorSolverOption::SO_KNITRO) {
#if defined(PGO_HAS_KNITRO)
    if (needRenew || !da->problem) {
      if (constraints) {
        da->problem = std::make_unique<KnitroProblem>(energy, constraints);
        da->problem->setInit(x);
        da->problem->setRange(xlow, xhi);
        da->problem->setConstraintsRange(clow, chi);
      }
      else {
        da->problem = std::make_shared<KnitroProblem>(energy);
        da->problem->setInit(x);
        da->problem->setRange(xlow, xhi);
      }

      da->opt = std::make_shared<KnitroOptimizer>(da->problem.get());
      if (solverConfigFilename && strlen(solverConfigFilename))
        da->opt->setConfigFile(solverConfigFilename);

      if (nIter > 0)
        da->opt->setMaxIter(nIter);

      if (eps > 0) {
        // solver->opt[stage]->setFeasTol(eps);
        da->opt->setOptTol(eps);
      }

      if (verbose >= 0) {
        da->opt->setVerbose(verbose);
      }

      da->opt->init();
    }
    else {
      if (constraints) {
        da->opt->setxInit(x.data());
        da->opt->setxRange(xlow.data(), xhi.data());
        da->opt->setcRange(clow.data(), chi.data());
      }
      else {
        da->opt->setxInit(x.data());
        da->opt->setxRange(xlow.data(), xhi.data());
      }
    }

    int solverRet = da->opt->solve();

    /*if (verbose == 0)
    da->opt->printInfo();*/

    x.noalias() = Eigen::Map<const ES::VXd>(da->opt->getx(), da->problem->getn());

    if (constraints) {
      lambda = Eigen::Map<const ES::VXd>(da->opt->getlambda(), da->problem->getm());
      g = Eigen::Map<const ES::VXd>(da->opt->getg(), da->problem->getm());
    }

    // if (solverRet != 0) {
    //   while (1) {
    //     ES::VXd grad(x.size());
    //     grad.setZero();
    //     energy->gradient(x, grad);

    //    for (ES::IDX i = 0; i < x.size(); i++) {
    //      if (xlow[i] == xhi[i]) {
    //        grad[i] = 0.0;
    //      }
    //    }

    //    double f0 = energy->func(x);
    //    double alpha = 1.0;
    //    ES::VXd x1 = x;
    //    LGI << "f0: " << f0;
    //    while (alpha > 1e-4) {
    //      x1 = x + alpha * grad;
    //      double f1 = energy->func(x1);
    //      LGI << "a: " << alpha << "; f:" << f1;
    //      if (f1 < f0) {
    //        break;
    //      }

    //      alpha *= 0.5;
    //    }

    //    if (alpha > 1e-4) {
    //      x = x1;
    //      grad.setZero();
    //      energy->gradient(x, grad);
    //      LGI << "Grad: " << grad.norm();
    //    }
    //    else {
    //      break;
    //    }
    //  }
    //  exit(1);
    //}

    return solverRet;
#else
    throw std::invalid_argument("No available selected solver");
    /*


    if (needRenew) {
      if (constraints) {
        da->L = std::make_shared<Lagrangian>(energy, constraints);
      }
      else {
        da->L = energy;
      }

      da->xinit.resize(da->L->getNumDOFs());

      NewtonRaphsonSolver::SolverParam solverParam;
      solverParam.lsm = NewtonRaphsonSolver::LSM_THUENTEMORE;
      // solverParam.lsm = NewtonRaphsonSolver::LSM_BRENTS;
      // solverParam.sst = NewtonRaphsonSolver::SST_SUBITERATION_ONE;

      da->solver = std::make_shared<NewtonRaphsonSolver>(da->xinit.data(), solverParam, da->L, da->fixedDOFs[0]);
    }

    if (needRenew == false && sameFixedDOFs == false) {
      da->solver->setFixedDOFs(da->fixedDOFs[0]);
    }

    da->xinit.head(energy->getNumDOFs()) = x;
    if (constraints) {
      if (lambda.size())
        da->xinit.tail(constraints->getNumConstraints()) = lambda;
      else
        da->xinit.tail(constraints->getNumConstraints()) = ES::VXd::Zero(constraints->getNumConstraints());
    }

    int ret = da->solver->solve(da->xinit.data(), nIter, eps, verbose);

    x.noalias() = da->xinit.head(energy->getNumDOFs());

    if (lambda.size())
      lambda.noalias() = da->xinit.tail(constraints->getNumConstraints());

    if (g.size() && constraints)
      constraints->func(x, g);

    return 0;
    */
#endif
  }
  else if (op == TimeIntegratorSolverOption::SO_IPOPT) {
#if defined(USE_IPOPT)
    return solveDirect(x, g, lambda, xlow, xhi, clow, chi, energy, constraints, nIter, eps, verbose, nullptr, op);
#else
    throw std::invalid_argument("No IPOPT solver");
#endif
  }
  else if (op == TimeIntegratorSolverOption::SO_NEWTON) {
    if (constraints) {
      throw std::invalid_argument("The solver does not support constraints");
    }
    else {
      da->fixedDOFs.reserve(energy->getNumDOFs());
      da->fixedDOFs.clear();

      for (int i = 0; i < energy->getNumDOFs(); i++)
        if (xlow[i] == xhi[i])
          da->fixedDOFs.push_back((int)i);

      da->fixedValues.resize(da->fixedDOFs.size());
      for (int i = 0; i < (int)da->fixedDOFs.size(); i++) {
        da->fixedValues[i] = xlow[da->fixedDOFs[i]];
      }

      if (needRenew || !da->newtonSolver) {
        NewtonRaphsonSolver::SolverParam sp;

        da->newtonSolver = std::make_shared<NewtonRaphsonSolver>(x.data(), sp, energy, da->fixedDOFs, da->fixedValues.data());
      }
      else {
        da->newtonSolver->setFixedDOFs(da->fixedDOFs, da->fixedValues.data());
      }

      int ret = da->newtonSolver->solve(x.data(), nIter, eps, verbose);

      return ret;
    }
  }
  else {
    throw std::invalid_argument("None of the solver is selected.");
  }
}

int TimeIntegratorSolver::solveDirect(ES::VXd &x, ES::VXd &g, ES::VXd &lambda,
  const ES::VXd &xlow, const ES::VXd &xhi, const ES::VXd &clow, const ES::VXd &chi,
  std::shared_ptr<const PotentialEnergy> energy,
  std::shared_ptr<const ConstraintFunctions> constraints,
  int nIter, double eps, int verbose, const char *solverConfigFilename,
  TimeIntegratorSolverOption solverOption)
{
  int solverRet = 0;

  if (constraints) {
    if (solverOption == TimeIntegratorSolverOption::SO_IPOPT)
      solverRet = EnergyOptimizer::minimize(x, energy, xlow, xhi,
        lambda, g, constraints, clow, chi,
        EnergyOptimizer::SolverType::ST_IPOPT, nIter, eps, verbose);
    else if (solverOption == TimeIntegratorSolverOption::SO_NEWTON)
      solverRet = EnergyOptimizer::minimize(x, energy, xlow, xhi,
        lambda, g, constraints, clow, chi,
        EnergyOptimizer::SolverType::ST_NEWTON, nIter, eps, verbose);
    else if (solverOption == TimeIntegratorSolverOption::SO_KNITRO)
      solverRet = EnergyOptimizer::minimizeUsingKnitro(x, energy, xlow, xhi,
        lambda, g, constraints, clow, chi,
        nIter, eps, verbose, solverConfigFilename);
  }
  else {
    if (solverOption == TimeIntegratorSolverOption::SO_IPOPT)
      solverRet = EnergyOptimizer::minimize(x, energy, xlow, xhi,
        EnergyOptimizer::SolverType::ST_IPOPT, nIter, eps, verbose);
    else if (solverOption == TimeIntegratorSolverOption::SO_NEWTON)
      solverRet = EnergyOptimizer::minimize(x, energy, xlow, xhi,
        EnergyOptimizer::SolverType::ST_NEWTON, nIter, eps, verbose);
    else if (solverOption == TimeIntegratorSolverOption::SO_KNITRO)
      solverRet = EnergyOptimizer::minimizeUsingKnitro(x, energy, xlow, xhi,
        lambda, g, nullptr, ES::VXd(), ES::VXd(), nIter, eps, verbose,
        solverConfigFilename);
  }

  return solverRet;
}
