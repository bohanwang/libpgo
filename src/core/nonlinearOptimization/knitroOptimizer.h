/*
author: Bohan Wang
copyright to USC
*/

#pragma once

namespace pgo
{
namespace NonlinearOptimization
{

class KnitroProblem;
class KnitroHandles;

class KnitroOptimizer
{
public:
  KnitroOptimizer(KnitroProblem *problem);
  virtual ~KnitroOptimizer();

  void setConfigFile(const char *filename);
  void setMaxIter(int maxIter);
  void setFeasTol(double eps);
  void setOptTol(double eps);
  void setVerbose(int verbose);

  void setxInit(const double *x);

  void setxRange(const double *xlow, const double *xhi);
  void setcRange(const double *clow, const double *chi);
  void enableWarmStart(int enable);

  void enableMultiEvaluation(int enable);
  void honorBoundary(int enable);
  void emphasisFeasibility(int opt); // 0: no emphasis, 1: emphasis ineq; 2: emphasis feas; 3; all

  void init();

  int solve();

  const double *getx() const { return x; }
  const double *getlambda() const { return lambda; }
  const double *getg() const { return g; }

  void printInfo() const;

  double getFeasibilityError() const { return feasError; }
  double getOptimizationError() const { return optError; }

protected:
  void initQuadraticProblem();

  KnitroHandles *handles;
  double *x, *lambda, *g;
  double feasError, optError;
};

}
}
