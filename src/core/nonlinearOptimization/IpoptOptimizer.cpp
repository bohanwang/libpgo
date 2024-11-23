/*
author: Bohan Wang
copyright to USC
*/

#include "IpoptOptimizer.h"
#include "IpoptProblem.h"

#include "pgoLogging.h"

#include <coin-or/IpIpoptApplication.hpp>
#include <coin-or/IpSolveStatistics.hpp>

using namespace pgo;
using namespace pgo::NonlinearOptimization;

IpoptOptimizer::IpoptOptimizer(IpoptProblem_ptr prob):
  problem(prob)
{
}

IpoptOptimizer::~IpoptOptimizer()
{
}

int IpoptOptimizer::init()
{
  // Create an instance of the IpoptApplication
  //
  // We are using the factory, since this allows us to compile this
  // example with an Ipopt Windows DLL
  app = IpoptApplicationFactory();

  app->Options()->SetIntegerValue("print_level", printLevel);
  app->Options()->SetNumericValue("tol", tol);
  app->Options()->SetIntegerValue("max_iter", maxIter);
  app->Options()->SetNumericValue("nlp_lower_bound_inf", -IpoptProblem::inf());
  app->Options()->SetNumericValue("nlp_upper_bound_inf", IpoptProblem::inf());
  app->Options()->SetStringValue("linear_solver", "pardisomkl");
  // app->Options()->SetStringValue("nlp_scaling_method", "none");

  if (linearConstraints) {
    app->Options()->SetStringValue("jac_c_constant", "yes");
    app->Options()->SetStringValue("jac_d_constant", "yes");
  }

  // Initialize the IpoptApplication and process the options
  Ipopt::ApplicationReturnStatus status;
  status = app->Initialize();
  if (status != Ipopt::Solve_Succeeded) {
    SPDLOG_LOGGER_ERROR(pgo::Logging::lgr(), "*** Error during initialization!");
    return 1;
  }

  return 0;
}

int IpoptOptimizer::solve()
{
  Ipopt::ApplicationReturnStatus status;

  if (firstTime) {
    status = app->OptimizeTNLP(problem);
    firstTime = 0;
  }
  else {
    status = app->ReOptimizeTNLP(problem);
  }

  if (status == Ipopt::Solve_Succeeded) {
    if (printLevel > 0) {
      // Retrieve some statistics about the solve
      Ipopt::Index iter_count = app->Statistics()->IterationCount();
      SPDLOG_LOGGER_INFO(pgo::Logging::lgr(), "*** The problem solved in {} iterations!", iter_count);

      Ipopt::Number final_obj = app->Statistics()->FinalObjective();
      SPDLOG_LOGGER_INFO(pgo::Logging::lgr(), "*** The final value of the objective function is {}", final_obj);
    }

    return 0;
  }

  return status;
}