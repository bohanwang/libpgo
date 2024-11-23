/*
author: Bohan Wang
copyright to USC
*/

#pragma once

#include <Eigen/Dense>

#include <coin-or/IpIpoptApplication.hpp>
#include <coin-or/IpSmartPtr.hpp>

namespace pgo
{
namespace NonlinearOptimization
{
class IpoptProblem;
typedef Ipopt::SmartPtr<IpoptProblem> IpoptProblem_ptr;

class IpoptOptimizer
{
public:
  IpoptOptimizer(IpoptProblem_ptr problem);
  virtual ~IpoptOptimizer();

  void setVerbose(int verbose) { printLevel = verbose; }
  void setTol(double val) { tol = val; }
  void setMaxIter(int numIter) { maxIter = numIter; }
  void setLinearConstraints(bool isLinear) { linearConstraints = isLinear; }
  int init();

  int solve();

protected:
  IpoptProblem_ptr problem;
  Ipopt::SmartPtr<Ipopt::IpoptApplication> app;
  int firstTime = 1;

  int printLevel = 5;
  double tol = 1e-8;
  int maxIter = 1000;

  bool linearConstraints = false;
};

}  // namespace NonlinearOptimization
}  // namespace pgo