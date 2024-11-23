/*
author: Bohan Wang
copyright to USC
*/

#include "minimizeEnergy.h"
#include "lagrangian.h"
#include "NewtonRaphsonSolver.h"
#include "potentialEnergy.h"

#ifdef USE_IPOPT
#  include "IpoptProblem.h"
#  include "IpoptOptimizer.h"
#endif  // USE_IPOPT

#ifdef USE_KNITRO
#  include "knitroProblem.h"
#  include "knitroOptimizer.h"

namespace pgo
{
namespace NonlinearOptimization
{
namespace EnergyOptimizer
{
struct KnitroData
{
  std::shared_ptr<KnitroProblem> problem;
  std::shared_ptr<KnitroOptimizer> optimizer;
};
}  // namespace EnergyOptimizer
}  // namespace NonlinearOptimization
}  // namespace pgo
#endif

#include <iostream>

using namespace pgo;
using namespace pgo::NonlinearOptimization;
using namespace pgo::NonlinearOptimization::EnergyOptimizer;

namespace ES = pgo::EigenSupport;

namespace pgo
{
namespace NonlinearOptimization
{
namespace EnergyOptimizer
{
int minimizeUsingIpopt(EigenSupport::RefVecXd x, PotentialEnergy_const_p energy, EigenSupport::ConstRefVecXd xlow, EigenSupport::ConstRefVecXd xhi,
  EigenSupport::RefVecXd lambda, EigenSupport::RefVecXd g, ConstraintFunctions_const_p constraints, EigenSupport::ConstRefVecXd clow, EigenSupport::ConstRefVecXd chi,
  int maxIter, double eps, int verbose);

int minimizeUsingNewton(EigenSupport::RefVecXd x, PotentialEnergy_const_p energy, EigenSupport::ConstRefVecXd xlow, EigenSupport::ConstRefVecXd xhi,
  EigenSupport::RefVecXd lambda, EigenSupport::RefVecXd g, ConstraintFunctions_const_p constraints, EigenSupport::ConstRefVecXd clow, EigenSupport::ConstRefVecXd chi,
  int maxIter, double eps, int verbose);

}  // namespace EnergyOptimizer
}  // namespace NonlinearOptimization
}  // namespace pgo

int EnergyOptimizer::minimize(EigenSupport::RefVecXd x, PotentialEnergy_const_p energy, EigenSupport::ConstRefVecXd xlow, EigenSupport::ConstRefVecXd xhi,
  EigenSupport::RefVecXd lambda, EigenSupport::RefVecXd g, ConstraintFunctions_const_p constraints, EigenSupport::ConstRefVecXd clow, EigenSupport::ConstRefVecXd chi,
  SolverType solverType, int maxIter, double eps, int verbose)
{
  if (solverType == SolverType::ST_IPOPT) {
    return minimizeUsingIpopt(x, energy, xlow, xhi, lambda, g, constraints, clow, chi, maxIter, eps, verbose);
  }
  else if (solverType == SolverType::ST_NEWTON) {
    return minimizeUsingNewton(x, energy, xlow, xhi, lambda, g, constraints, clow, chi, maxIter, eps, verbose);
  }
  else if (solverType == SolverType::ST_KNITRO) {
    return minimizeUsingKnitro(x, energy, xlow, xhi, lambda, g, constraints, clow, chi, maxIter, eps, verbose, nullptr);
  }

  return -1;
}

int EnergyOptimizer::minimize(EigenSupport::RefVecXd x, PotentialEnergy_const_p energy, EigenSupport::ConstRefVecXd xlow, EigenSupport::ConstRefVecXd xhi,
  SolverType solverType, int maxIter, double eps, int verbose)
{
  ES::VXd lambda, g, clow, chi;
  if (solverType == SolverType::ST_IPOPT) {
    return minimizeUsingIpopt(x, energy, xlow, xhi, lambda, g, nullptr, clow, chi, maxIter, eps, verbose);
  }
  else if (solverType == SolverType::ST_NEWTON) {
    return minimizeUsingNewton(x, energy, xlow, xhi, lambda, g, nullptr, clow, chi, maxIter, eps, verbose);
  }
  else if (solverType == SolverType::ST_KNITRO) {
    return minimizeUsingKnitro(x, energy, xlow, xhi, lambda, g, nullptr, clow, chi, maxIter, eps, verbose, nullptr);
  }

  return -1;
}

int EnergyOptimizer::minimizeUsingIpopt(EigenSupport::RefVecXd x, PotentialEnergy_const_p energy, EigenSupport::ConstRefVecXd xlow, EigenSupport::ConstRefVecXd xhi,
  EigenSupport::RefVecXd lambda, EigenSupport::RefVecXd g, ConstraintFunctions_const_p constraints, EigenSupport::ConstRefVecXd clow, EigenSupport::ConstRefVecXd chi,
  int maxIter, double eps, int verbose)
{
#if defined(USE_IPOPT)
  if (constraints) {
    Ipopt::SmartPtr<IpoptProblem> problem = new IpoptProblem(energy, constraints);
    problem->setInit(x);
    problem->setRange(xlow, xhi);
    problem->setConstraintsRange(clow, chi);

    IpoptOptimizer solver(problem);
    solver.setMaxIter(maxIter);
    solver.setTol(eps);
    solver.setVerbose(verbose);

    if (constraints->isLinear())
      solver.setLinearConstraints(true);

    solver.init();

    int solverRet = solver.solve();
    x = problem->getFinalx();

    if (g.size() > 0)
      g = problem->getFinalg();

    if (lambda.size() > 0)
      lambda = problem->getFinalLambda();

    return solverRet;
  }
  else {
    Ipopt::SmartPtr<IpoptProblem> problem = new IpoptProblem(energy);
    problem->setInit(x);
    problem->setRange(xlow, xhi);

    IpoptOptimizer solver(problem);
    solver.setMaxIter(maxIter);
    solver.setTol(eps);
    solver.setVerbose(verbose);

    solver.init();

    int solverRet = solver.solve();
    x = problem->getFinalx();

    return solverRet;
  }
#else
  (void)x, (void)energy;
  (void)xlow, (void)xhi;
  (void)lambda, (void)g, (void)constraints, (void)clow, (void)chi;
  (void)maxIter, (void)eps, (void)verbose;
  throw std::runtime_error("No Ipopt Module");
  return 1;
#endif
}

int EnergyOptimizer::minimizeUsingKnitro(EigenSupport::RefVecXd x, PotentialEnergy_const_p energy, EigenSupport::ConstRefVecXd xlow, EigenSupport::ConstRefVecXd xhi,
  EigenSupport::RefVecXd lambda, EigenSupport::RefVecXd g, ConstraintFunctions_const_p constraints, EigenSupport::ConstRefVecXd clow, EigenSupport::ConstRefVecXd chi,
  int maxIter, double eps, int verbose, const char *filename, int parallelEval, CallbackFunc callback, KnitroData ** /*data*/, double feasTol)
{
#if defined(USE_KNITRO)
  if (constraints) {
    std::unique_ptr<KnitroProblem> problem = std::make_unique<KnitroProblem>(energy, constraints);
    problem->setInit(x);
    problem->setRange(xlow, xhi);
    problem->setConstraintsRange(clow, chi);

    KnitroProblem::CallbackFunc func = [&](KnitroProblem *const ptr) {
      callback(ptr->getCurrentx().data(), (int)ptr->getCurrentx().size(),
        ptr->getCurrentlambda().data(), (int)ptr->getCurrentlambda().size());
    };

    if (callback)
      problem->setIterationCallback(func);

    KnitroOptimizer solver(problem.get());

    if (filename) {
      if (verbose >= 0)
        std::cout << "Solver config filename " << filename << std::endl;

      solver.setConfigFile(filename);
    }

    if (maxIter > 0)
      solver.setMaxIter(maxIter);

    if (eps > 0) {
      // solver.setFeasTol(eps);
      solver.setOptTol(eps);
    }

    if (feasTol > 0)
      solver.setFeasTol(feasTol);

    if (verbose >= 0) {
      solver.setVerbose(verbose);
    }
    else {
      solver.setVerbose(0);
    }

    solver.enableMultiEvaluation(parallelEval);

    solver.init();

    int solverRet = solver.solve();

    if (verbose == 0)
      solver.printInfo();

    x = Eigen::Map<const ES::VXd>(solver.getx(), problem->getn());

    if (lambda.size() > 0)
      lambda = Eigen::Map<const ES::VXd>(solver.getlambda(), problem->getm());

    if (g.size() > 0)
      g = Eigen::Map<const ES::VXd>(solver.getg(), problem->getm());

    return solverRet;
  }
  else {
    std::shared_ptr<KnitroProblem> problem = std::make_shared<KnitroProblem>(energy);
    problem->setInit(x);
    problem->setRange(xlow, xhi);

    KnitroProblem::CallbackFunc func = [&](KnitroProblem *const ptr) {
      callback(ptr->getCurrentx().data(), (int)ptr->getCurrentx().size(),
        nullptr, 0);
    };

    if (callback)
      problem->setIterationCallback(func);

    KnitroOptimizer solver(problem.get());

    if (filename) {
      if (verbose >= 0)
        std::cout << "Solver config filename " << filename << std::endl;
      solver.setConfigFile(filename);
    }

    if (maxIter > 0)
      solver.setMaxIter(maxIter);

    if (eps > 0) {
      solver.setFeasTol(eps);
      solver.setOptTol(eps);
    }

    if (verbose >= 0) {
      solver.setVerbose(verbose);
    }
    else {
      solver.setVerbose(0);
    }

    solver.init();

    int solverRet = solver.solve();

    if (verbose == 0)
      solver.printInfo();

    x = Eigen::Map<const ES::VXd>(solver.getx(), problem->getn());

    return solverRet;
  }
#else
  std::cout << "Knitro not found. Using IPOPT" << std::endl;
  return minimizeUsingIpopt(x, energy, xlow, xhi, lambda, g, constraints, clow, chi, maxIter, eps, verbose);
#endif
}

int EnergyOptimizer::minimizeUsingKnitroDense(EigenSupport::RefVecXd x, PotentialEnergyDense_const_p energy, EigenSupport::ConstRefVecXd xlow, EigenSupport::ConstRefVecXd xhi,
  EigenSupport::RefVecXd lambda, EigenSupport::RefVecXd g, ConstraintFunctionsDense_const_p constraints, EigenSupport::ConstRefVecXd clow, EigenSupport::ConstRefVecXd chi,
  int maxIter, double eps, int verbose, const char *filename, int parallelEval, CallbackFunc callback, KnitroData ** /*data*/, double feasTol)
{
#if defined(USE_KNITRO)
  if (constraints) {
    std::unique_ptr<KnitroProblem> problem = std::make_unique<KnitroProblem>(energy, constraints);
    problem->setInit(x);
    problem->setRange(xlow, xhi);
    problem->setConstraintsRange(clow, chi);

    KnitroProblem::CallbackFunc func = [&](KnitroProblem *const ptr) {
      callback(ptr->getCurrentx().data(), (int)ptr->getCurrentx().size(),
        ptr->getCurrentlambda().data(), (int)ptr->getCurrentlambda().size());
    };

    if (callback)
      problem->setIterationCallback(func);

    KnitroOptimizer solver(problem.get());

    if (filename) {
      if (verbose >= 0)
        std::cout << "Solver config filename " << filename << std::endl;

      solver.setConfigFile(filename);
    }

    if (maxIter > 0)
      solver.setMaxIter(maxIter);

    if (eps > 0) {
      // solver.setFeasTol(eps);
      solver.setOptTol(eps);
    }

    if (feasTol > 0)
      solver.setFeasTol(feasTol);

    if (verbose >= 0) {
      solver.setVerbose(verbose);
    }
    else {
      solver.setVerbose(0);
    }

    solver.enableMultiEvaluation(parallelEval);

    solver.init();

    int solverRet = solver.solve();

    if (verbose == 0)
      solver.printInfo();

    x = Eigen::Map<const ES::VXd>(solver.getx(), problem->getn());

    if (lambda.size() > 0)
      lambda = Eigen::Map<const ES::VXd>(solver.getlambda(), problem->getm());

    if (g.size() > 0)
      g = Eigen::Map<const ES::VXd>(solver.getg(), problem->getm());

    return solverRet;
  }
  else {
    std::shared_ptr<KnitroProblem> problem = std::make_shared<KnitroProblem>(energy);
    problem->setInit(x);
    problem->setRange(xlow, xhi);

    KnitroProblem::CallbackFunc func = [&](KnitroProblem *const ptr) {
      callback(ptr->getCurrentx().data(), (int)ptr->getCurrentx().size(),
        nullptr, 0);
    };

    if (callback)
      problem->setIterationCallback(func);

    KnitroOptimizer solver(problem.get());

    if (filename) {
      if (verbose >= 0)
        std::cout << "Solver config filename " << filename << std::endl;
      solver.setConfigFile(filename);
    }

    if (maxIter > 0)
      solver.setMaxIter(maxIter);

    if (eps > 0) {
      solver.setFeasTol(eps);
      solver.setOptTol(eps);
    }

    if (verbose >= 0) {
      solver.setVerbose(verbose);
    }
    else {
      solver.setVerbose(0);
    }

    solver.init();

    int solverRet = solver.solve();

    if (verbose == 0)
      solver.printInfo();

    x = Eigen::Map<const ES::VXd>(solver.getx(), problem->getn());

    return solverRet;
  }
#else
  std::cout << "Knitro not found. exit(1)" << std::endl;
  exit(1);

  return 1;
#endif
}

int EnergyOptimizer::minimizeUsingApproximateActiveSet(EigenSupport::RefVecXd x, PotentialEnergy_const_p energy, EigenSupport::ConstRefVecXd xlow, EigenSupport::ConstRefVecXd xhi,
  EigenSupport::RefVecXd lambda, EigenSupport::RefVecXd g, ConstraintFunctions_const_p constraints, EigenSupport::ConstRefVecXd clow, EigenSupport::ConstRefVecXd chi,
  int maxIter, double eps, double equalityThreshold, int verbose)
{
#if defined(USE_IPOPT)
  if (constraints) {
    Ipopt::SmartPtr<IpoptProblem> problem = new IpoptProblem(energy, constraints);
    problem->setRange(xlow, xhi);

    ES::VXd clow1 = clow, chi1 = chi;

    auto minf = [&]() {
      problem->setInit(x);
      problem->setConstraintsRange(clow1, chi1);

      IpoptOptimizer solver(problem);
      solver.setMaxIter(maxIter);
      solver.setTol(eps);
      solver.setVerbose(verbose);

      if (constraints->isLinear())
        solver.setLinearConstraints(true);
      solver.init();

      int solverRet = solver.solve();
      x = problem->getFinalx();
      g = problem->getFinalg();
      lambda = problem->getFinalLambda();

      if (solverRet != 0)
        return solverRet;
      else
        return 0;
    };

    // solve first time
    int solverRet = minf();
    if (solverRet < 0)
      return solverRet;

    for (Eigen::Index i = 0; i < g.size(); i++) {
      if (g[i] < 1e-10) {
        clow1[i] = 0.0;
        chi1[i] = 0.0;
      }
    }

    while (1) {
      solverRet = minf();
      if (solverRet < 0)
        return solverRet;

      // deactivate all pulling constraints
      bool isGood = true;
      int counter0 = 0;
      for (Eigen::Index i = 0; i < lambda.size(); i++) {
        if (lambda[i] > 1e-10) {
          clow1[i] = -10;
          chi1[i] = 10;

          isGood = false;
          counter0++;
        }
      }

      // activate all pushing constraints
      int counter1 = 0;
      for (Eigen::Index i = 0; i < g.size(); i++) {
        if (lambda[i] < -1e-10 && g[i] < -1e-10) {
          clow1[i] = 0;
          chi1[i] = 0;

          isGood = false;
          counter1++;
        }
      }

      std::cout << "# pulling forces: " << counter0 << '/' << lambda.size() << std::endl;
      std::cout << "# violated: " << counter1 << '/' << lambda.size() << std::endl;
      std::cout << "# eq: " << std::count_if(chi1.data(), chi1.data() + chi1.size(), [](double val) { return val == 0; }) << std::endl;

      if (isGood)
        break;
    }

    std::cout << std::endl;

    return solverRet;
  }
  else {
    Ipopt::SmartPtr<IpoptProblem> problem = new IpoptProblem(energy);
    problem->setInit(x);
    problem->setRange(xlow, xhi);

    IpoptOptimizer solver(problem);
    solver.setMaxIter(maxIter);
    solver.setTol(eps);
    solver.setVerbose(verbose);

    solver.init();

    int solverRet = solver.solve();
    x = problem->getFinalx();

    return solverRet;
  }
#else
  (void)x, (void)energy;
  (void)xlow, (void)xhi;
  (void)lambda, (void)g, (void)constraints, (void)clow, (void)chi;
  (void)maxIter, (void)eps, (void)verbose;
  (void)equalityThreshold;
  throw std::runtime_error("No Ipopt Module");

  return 1;
#endif
}

int EnergyOptimizer::minimizeUsingNewton(EigenSupport::RefVecXd x, PotentialEnergy_const_p energy, EigenSupport::ConstRefVecXd xlow, EigenSupport::ConstRefVecXd xhi,
  EigenSupport::RefVecXd lambda, EigenSupport::RefVecXd g, ConstraintFunctions_const_p constraints, EigenSupport::ConstRefVecXd clow, EigenSupport::ConstRefVecXd chi,
  int maxIter, double eps, int verbose)
{
  if (constraints) {
    std::cout << "Assume equality constraints C=0. Ignore clow, chi" << std::endl;
    (void)clow, (void)chi;
  }

  std::vector<int> fixedDOFs;
  fixedDOFs.reserve(energy->getNumDOFs());
  for (int i = 0; i < energy->getNumDOFs(); i++)
    if (xlow[i] == xhi[i])
      fixedDOFs.push_back((int)i);

  std::shared_ptr<const PotentialEnergy> L;
  if (constraints) {
    L = std::make_shared<Lagrangian>(energy, constraints);
  }
  else {
    L = energy;
  }

  NewtonRaphsonSolver::SolverParam solverParam;
  solverParam.lsm = NewtonRaphsonSolver::LSM_SIMPLE;
  // solverParam.lsm = NewtonRaphsonSolver::LSM_BRENTS;
  // solverParam.sst = NewtonRaphsonSolver::SST_SUBITERATION_ONE;
  NewtonRaphsonSolver solver(x.data(), solverParam, L, fixedDOFs);
  ES::VXd xinit(L->getNumDOFs());
  xinit.head(energy->getNumDOFs()) = x;

  if (constraints) {
    if (lambda.size())
      xinit.tail(constraints->getNumConstraints()) = lambda;
    else
      xinit.tail(constraints->getNumConstraints()) = ES::VXd::Zero(constraints->getNumConstraints());
  }

  int ret = solver.solve(xinit.data(), maxIter, eps, verbose);

  x = xinit.head(energy->getNumDOFs());

  if (lambda.size())
    lambda = xinit.tail(constraints->getNumConstraints());

  if (g.size() && constraints)
    constraints->func(x, g);

  return ret;
}
