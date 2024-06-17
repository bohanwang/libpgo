// copyright: vegaFEM

#ifndef LINESEARCH_H
#define LINESEARCH_H

#include <vector>
#include <functional>
#include <cfloat>

// Implements several line search algorithm.
// In general, line search solves: min f(x) where x = x0 + alpha p
//
// 1) golden search:
// use golden partitioning to narrow down the search range and find a local minima
//
// 2) Brent's method:
// use a combination of parabola fit and golden search to find a local minima
//
// Some additional terms for following line search methods:
// - Sufficient decrease:
// f(x0 + a p) <= f(x0) + c1 a grad(f0)^T p, where grad(f0) = grad(f(x))|x=x0
// and c1 is in (0, 1). In practice, c1 = 1e-4
// - Curvature condition:
// grad f(x0 + a p)^T p >= c2 grad f0^T p, where c2 is in (c1, 1).
// Typically c2 is 0.9 when Newton or quasi-Newton, and 0.1 when nonlinear
// conjugate gradient method
//
// - The Wolfe conditions:
// Combination of the sufficient decrease and curvature condition.
//
// - The strong Wolfe conditions:
// Combination of the sufficient decrease and
// |grad f(x0 + a p)^T p| <= c2 |grad f0^T p|
//
// For smooth and lower-bounded functions it is guaranteed that there are alphas
// satisfying the Wolfe conditions and the Wolfe conditions
//
// - The Goldstein conditions
// f(x0) + (1-c) a grad(f0) p <= f(x0 + a p) <= f(x0) + c a grad(f0) p,
// where c is in (0, 0.5).
// The Goldstein conditions are often used in Newton-type methods but
// are not well suited for quasi-Newton methods.
//
// 3) Backtracking
// Reduce step size until it meets the sufficient decrease condition
// It is welled suited for Newton methods but less appropriated for
// quasi-Newton and conjugate gradient methods.

class LineSearch
{
public:
  
  // if f != nullptr, return the function objective
  // if gradient != nullptr, return the function gradient
  using EvaluateFunction = std::function<int(const double * x, double * f, double * gradient)>;

  LineSearch(int n, EvaluateFunction f, int verbose = 0);
  void setBounds(double lower, double upper) { alphaLowerBound = lower; alphaUpperBound = upper; }
  void setMaxIterations(int mi) { maxIter = mi; }

  void setAbsolutePrecision(double p) { absPrecision = p; }

  void setToForwardSearch(int numIter, double rho) { forwardIter = numIter; backtrack_rho = rho; }

  double f1D(double alpha) const; // evaluate f(x0+alpha p)

  struct Result
  {
    double f = 0.0;
    double alpha = 0.0;
    int numIter = 0;
  };

  Result golden(const double * x0, const double * p, double fx0, double initialAlpha = 1.0) const;

  Result BrentsMethod(const double * x0, const double * p, double fx0, double initialAlpha = 1.0) const;

  // until f(x0 + a p) <= fx0 + c a gx0^T p, do a = rho a
  // assert c \in (0, 1) and rho \in (0, 1)
  // In practice, c is a small value like 1e-4
  Result backtracking(const double * x0, const double * p, double fx0, const double * gx0, double c, double rho,
      double initialAlpha = 1.0) const;

  void setVerbose(int v) { verbose = v; }

protected:

  // input: a and b and fa and fb
  // output: a triplet [a, b, c] such that: f1D(b) <= min( f1D(a), f1D(b) ) and also return f1D(a), f1D(b), f1D(c)
  // function return 0 on success and 1 on failure to find the triplet
  int initializeTriplet(double & a, double & b, double & c, double & fa, double & fb, double & fc) const;

  int n = 0;
  int maxIter = 100;
  EvaluateFunction func;

  double alphaLowerBound = -DBL_MAX;
  double alphaUpperBound = DBL_MAX;
  double absPrecision = 0; // absolute precision theshold to control the precision of the result
  int forwardIter = -1; // control whether to maintain the initial direction p
                        // if >=0, forwardIter is the max number of iterations used to check initial direction
  double backtrack_rho = 0.8;

  int verbose = 0;
  // helper buffers
  mutable std::vector<double> f1DBuffer;
  mutable int f1DCount = 0; // used to count how many f1D called

  // x0 and p used in the optimization: min f(x) where x = x0 + alpha p
  mutable const double * x0 = nullptr;
  mutable const double * p = nullptr;
};



#endif
