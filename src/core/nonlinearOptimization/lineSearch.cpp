/*
copyright: vegafem
*/
#include "lineSearch.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <algorithm>

using namespace std;

using Result = LineSearch::Result;

static const double goldenRatio = 1.61803398875;
static const double goldenPercentage = goldenRatio - 1.0; // 0.618
static const double goldenSmallerPercentage = 2.0 - goldenRatio; // 0.382

// get a new value starting from a and b
static inline double takeGoldenStep(double a, double b)
{
  return b + (b-a) * goldenRatio;
}

// clamp a if fabs(a) is below th
// assert(th >= 0)
static inline double clampFromZero(double a, double th)
{
  if (a > th || a < -th)
    return a;

  if (a >= 0) return th;
  return -th;
}

// check whether a is between b and c: (b, c) or (c, b)
// assert b != c
static inline bool inBetween(double a, double b, double c)
{
  return (a-b) * (c-a) > 0.0;
}

static inline double clamp(double a, double lowerBound, double upperBound)
{
  return (a < lowerBound ? lowerBound : (a > upperBound ? upperBound : a));
}

static inline double vectorDotProduct(int n, const double * va, const double * vb)
{
  double ret = 0.0;
  for(int i = 0; i < n; i++)
    ret += va[i] * vb[i];
  return ret;
}

LineSearch::LineSearch(int n, EvaluateFunction f, int v) : n(n), func(f), verbose(v)
{
  f1DBuffer.resize(n);
}

Result LineSearch::golden(const double * x0, const double * p, double fx0, double initialAlpha) const
{
  this->x0 = x0;
  this->p = p;
  f1DCount = 0;

  double a = 0.0;
  double b = initialAlpha;
  double c = 0.0;
  double fa = fx0, fb = f1D(b), fc = 0.0;

  // number of iterations to control the direction
  for(int i = 0; i < forwardIter; i++)
  {
    if (fb <= fa)
      break;
    b *= backtrack_rho;
    fb = f1D(b);
  }
  if (fb > fa && forwardIter >= 0)
  {
    Result result;
    result.alpha = 0.0;
    result.f = fa;
    result.numIter = f1DCount;
    return result;
  }

  int ret = initializeTriplet(a, b, c, fa, fb, fc);
  if (ret != 0) // fail to find a valid triplet, use the smaller one between fa and fb
  {
    if (verbose)
      cout << "failed to find a triplet: " << a << " " << b << " " << c << " | " << fa << " " << fb << " "  << fc << endl;
    // in this case, c stores the lowest value
    Result result;
    result.alpha = c;
    result.f = fc;
    result.numIter = f1DCount;
    return result;
  }
  // now b is inBetween a and c, and fb <= min(fa, fc)
  assert(fb <= fa && fb <= fc);
  if (verbose)
    cout << "triplet: " << a << " " << b << " " << c << " | " << fa << " " << fb << " "  << fc << endl;
  if (fa == fb && fb == fc && fc == DBL_MAX)
  {
    if (verbose)
      cout << "error: triplet are all invalid values" << endl;
    Result result;
    result.alpha = 0.0;
    result.f = fx0;
    result.numIter = f1DCount;
    return result;
  }

  // now run golden search on [a, b, c]
  // we first make sure dist(a,b) < dist(b,c), i.e. b closer to a than to c
  // then we find a d between b and c and iterate
  if (fabs(a-b) > fabs(b-c)) // dist(b, c) < dist(a, b), b is closer to c
  { // then we swap a and c
    swap(a, c);
    swap(fa, fc);
  }
  // now get d
  double d = 0.0, fd = 0.0;
  d = b + (c - b) * goldenSmallerPercentage;
  fd = f1D(d);

  const double tol = 1e-4;
  for(int iter = 0; iter < maxIter && fabs(a-c) > max(absPrecision*2, tol * (fabs(b) + fabs(d))); iter++)
  {
    if (verbose)
      cout << "iter: " << iter << "/" << maxIter << " (" << a << " " << b << " " << c << ") (" <<
        fa << " " << fb << " " << fc << ")" << endl;
    if (fb < fd)
    {
      //f ^               ----            * END OF GRAPH
      //  | \            /                *
      //  |  ----     ---                 *
      //  |      \  /                     *
      //  |       --                      *
      //  |                               *
      //  ----------------------> x or -x *
      //     a    b  d      c

      // use [d, b, a] as the new [a, b, c]
      c = a;
      a = d;
      fc = fa;
      fa = fd;
    }
    else // fd < fb
    {
      //f ^                               * END OF GRAPH
      //  | \                             *
      //  |  ----          -----          *
      //  |      \        /               *
      //  |       --     /                *
      //  |         \   /                 *
      //  |          \_/                  *
      //  ----------------------> x or -x *
      //     a    b  d      c

      // use [b, d, c] as the new [a, b, c]
      a = b;
      b = d;
      fa = fb;
      fb = fd;
    }

    d = b + (c - b) * goldenSmallerPercentage;
    fd = f1D(d);
  } // end while loop

  // now fb or fd has the smallest value
  if (fb < fd)
  {
    swap(d, b);
    swap(fd, fb);
  }
  // now we make sure fd is the smallest value

  Result result;
  result.alpha = d;
  result.f = fd;
  result.numIter = f1DCount;
  return result;
}

Result LineSearch::BrentsMethod(const double * x0, const double * p, double fx0, double initialAlpha) const
{
  // TODO:
  // add control of absPrecisionDiameter
  this->x0 = x0;
  this->p = p;
  f1DCount = 0;

  double a = 0.0;
  double b = 0.5 * initialAlpha; // XXX 0.5 here is to be consistent with numerical recipe
  double c = 0.0;
  double fa = fx0, fb = f1D(b), fc = 0.0;
  if (isfinite(fx0) == false)
    cout << "Error, input function value fx0 is not finite: " << fx0 << endl;
  assert(isfinite(fx0));
  if (isfinite(fb) == false)
    cout << "Error, function value is not finite: " << fb << endl;
  assert(isfinite(fb));

  int ret = initializeTriplet(a, b, c, fa, fb, fc);


  if (isfinite(fa) == false)
    cout << "Error, function value is not finite: " << fa << endl;
  assert(isfinite(fa));
  if (isfinite(fb) == false)
    cout << "Error, function value is not finite: " << fb << endl;
  assert(isfinite(fb));
  if (isfinite(fc) == false)
    cout << "Error, function value is not finite: " << fc << endl;
  assert(isfinite(fc));


  if (ret != 0) // fail to find a valid triplet, use the smaller one between fa and fb
  {
    if (verbose)
      cout << "failed to find a triplet: " << a << " " << b << " " << c << " | " << fa << " " << fb << " "  << fc << endl;
    // in this case, c stores the lowest value
    Result result;
    result.alpha = c;
    result.f = fc;
    result.numIter = f1DCount;
    return result;
  }
  // now b is inBetween a and c, and fb <= min(fa, fc)
  if (!(fb <= fa && fb <= fc))
  {
    cout << "Error: " << fa << " " << fb << " " << fc <<endl;
  }
  assert(fb <= fa && fb <= fc);
  if (verbose)
    cout << "triplet: " << a << " " << b << " " << c << " | " << fa << " " << fb << " "  << fc << endl;
  if (fa == fb && fb == fc && fc == DBL_MAX)
  {
    if (verbose)
      cout << "error: triplet are all invalid values" << endl;
    Result result;
    result.alpha = 0.0;
    result.f = fx0;
    result.numIter = f1DCount;
    return result;
  }

  // now b is inBetween a and c, and fb <= min(fa, fc)

  // now let's run Brent's method:

  // create a bound from a and c
  double brentLowerBound = a, brentUpperBound = c;
  if (a > c)
    swap(brentLowerBound, brentUpperBound);
  if (verbose)
    cout << "Search alppha bounds: [" << alphaLowerBound << ", " << alphaUpperBound << "]" << endl;
  brentLowerBound = max(brentLowerBound, alphaLowerBound);
  brentUpperBound = min(brentUpperBound, alphaUpperBound);
  double clamped_b = clamp(b, brentLowerBound, brentUpperBound);
  if (clamped_b != b)
  {
    if (verbose)
      cout << "Clamp b into alpha bounds: " << b << "->" << clamped_b << endl;
    b = clamped_b;
    fb = f1D(b);
  }

  // x the point with the very least function valued found so far (or the most recent one in case of a tie)
  // x2: the point with the second least function value
  // x3: the previous value of x3
  double x = b, x2 = b, x3 = b;
  double fx = fb, fx2 = fb, fx3 = fb;

  double step1 = 0.0;  // the distance moved on the last step
  double step2 = 0.0; // the distance moved on the step before last

  for(int iter = 0; iter < maxIter; iter++)
  {
    if (verbose)
      cout << "iter: " << iter << "/" << maxIter << " [" << brentLowerBound << " " << brentUpperBound << "] x = " << x << ", f(x) = " <<
        fx << endl;
    if (absPrecision * 2 > brentUpperBound - brentLowerBound)
    {
      if (verbose)
        cout << "absPrecision " << absPrecision << " reached." << endl;
      break;
    }
    double middle = 0.5 * (brentLowerBound + brentUpperBound);
    double tol = 2 * 0.0001;
    double eps = 1e-10;
    double tol1 = tol * fabs(x) + eps;
    double tol2 = 2.0 * tol1;
    if (fabs(x-middle) + 0.5*(brentUpperBound - brentLowerBound) <= tol2) // the distance between the bounds are small enough
    {
      if (verbose)
        cout << "Relative tolerance of " << tol << " reached." << endl;
      break;
    }
    if (fabs(step2) > tol1)
    { // use x2, x, x3 to fit a parabola (x2=a, x=b, x3=c as in the function for parabola)
      // fit a parabola on a, b, c and compute the vertex as d
      // the equation is:
      //     (b^2-c^2)(fb-fa) - (b^2-a^2)(fb-fc)         (b-c)^2 (fb-fa) - (b-a)^2 (fb-fc)
      // d = ------------------------------------ = b - ------------------------------------- = b + stepSize
      //       2[(b-c)(fb-fa) - (b-a)(fb-fc)]               2[(b-c)(fb-fa) - (b-a)(fb-fc)]
      double bcfbfa = (x-x3) * (fx-fx2);
      double bafbfc = (x-x2) * (fx-fx3);
      double numerator = (x-x3) * bcfbfa - (x-x2) * bafbfc;
      double denominator = (-2) * (bcfbfa - bafbfc);
      // now the stepSize to go from x to the vertex of the parabola is: d - x = stepSize = numerator / denominator
      if (denominator < 0.0) // make sure denominator is >= 0.0
      {
        numerator *= (-1);
        denominator *= (-1);
      }

      double stepBeforeLast = step2;
      step2 = step1;
      // check whether the parabola fit is bad
      // fabs(numerator) >= fabs(0.5*denominator*stepBeforeLast): test whether the movement from the best current value x
      // that is larger or equal to half the movement of the step before last.
      if (fabs(numerator) >= fabs(0.5*denominator*stepBeforeLast) ||
          // also check whether the stepSize will move the value outside of the bounds
          numerator <= denominator*(brentLowerBound-x) || numerator >= denominator*(brentUpperBound-x))
      {
        // parabola is bad, use golden search
        step2 = (x >= middle ? brentLowerBound - x : brentUpperBound - x);
        step1 = goldenSmallerPercentage * step2;
      }
      else // parabola is good
      {
        const double threshold = 1e-20;
        step1 = numerator / clampFromZero(denominator, threshold);
        double xnew = x+step1;
        if (xnew - brentLowerBound < tol2 || brentUpperBound - xnew < tol2) // if xnew is too close to the bounds
        {
          if (middle - x >= 0.0)
            step1 = tol1;
          else
            step1 = -tol1;
        }
      }
    }
    else // use golden search
    {
      step2 = (x >= middle ? brentLowerBound - x : brentUpperBound - x);
      step1 = goldenSmallerPercentage * step2;
    }

    // if fabs(step1) < tol1, clamp step to be tol1
    double xnew = x + clampFromZero(step1, tol1);
    xnew = clamp(xnew, alphaLowerBound, alphaUpperBound);
    double fxnew = f1D(xnew);
    if (verbose)
      cout << "try new x: " << xnew << " f(x) = " << fxnew << endl;
    if (fxnew <= fx)
    {
      if (xnew >= x)
        brentLowerBound = x;
      else
        brentUpperBound = x;
      if (brentLowerBound > brentUpperBound)
      {
        if (verbose)
          cout << "Error: A " << brentLowerBound << " " << brentUpperBound << " " << x << " " << xnew << endl;
      }
      assert(brentLowerBound <= brentUpperBound);
      x3 = x2;
      x2 = x;
      x = xnew;
      fx3 = fx2;
      fx2 = fx;
      fx = fxnew;
    }
    else // fxnew > fx
    {
      if (xnew < x)
        brentLowerBound = xnew;
      else
        brentUpperBound = xnew;

      if (brentLowerBound > brentUpperBound)
      {
        if (verbose)
          cout << "Error: B " << brentLowerBound << " " << brentUpperBound << " " << x << " " << xnew << endl;
      }
      assert(brentLowerBound <= brentUpperBound);

      // w is the second least function value
      if (fxnew <= fx2 || x2 == x)
      {
        x3=x2;
        x2=xnew;
        fx3=fx2;
        fx2=fxnew;
      }
      else if (fxnew <= fx3 || x3 == x || x3 == x2)
      {
        x3=xnew;
        fx3=fxnew;
      }
    }
  } // end iter

  Result result;
  result.alpha = x;
  result.f = fx;
  result.numIter = f1DCount;
//  cout << "assert the returned value is correct:" << endl;
//  double tmpf = f1D(x);
//  cout << tmpf << " " << fx << " " << fabs(fx - tmpf) << endl;
  return result;
}

Result LineSearch::backtracking(const double * x0, const double * p, double fx0, const double * gx0, double c, double rho,
      double initialAlpha) const
{
  this->x0 = x0;
  this->p = p;
  f1DCount = 0;

  double cgp = c * vectorDotProduct(n, gx0, p);
//  double dir = 1.0;
  if (cgp > 0)
  {
    if (verbose)
      cout << "Warning: the direction is not downhill when backtracking: " << cgp << endl;
//    cgp *= -1.0;
//    dir = -1.0; // the direction of p is not correct, we have to reverse it
  }

  double alpha = initialAlpha;
  double fa = f1D(alpha);
  double minf = fa, mina = initialAlpha;
  int iter = 0;
  for(; (fa > fx0 + alpha * cgp) && iter < maxIter; iter++)
  {
    alpha *= rho;
    fa = f1D(alpha);
    if (fa < minf)
    {
      minf = fa;
      mina = alpha;
    }
  }
  if (minf > fx0 && verbose)
    cout << "Warning: the values found are larger than the intial function value: " << minf << " > " << fx0 << endl;

  Result result;
  result.alpha = mina;
  result.f = minf;
  result.numIter = f1DCount;
  return result;
}

int LineSearch::initializeTriplet(double & a, double & b, double & c, double & fa, double & fb, double & fc) const
{
  assert(a != b);

  if (isfinite(fa) == false)
      cout << "Error, in initializeTriplet, function value is not finite: " << fa << endl;
    assert(isfinite(fa));
  if (isfinite(fb) == false)
    cout << "Error, in initializeTriplet, function value is not finite: " << fb << endl;
  assert(isfinite(fb));

  if (fa < fb) // make sure fa >= fb
  {
    swap(fa, fb);
    swap(a, b);
  }

  c = takeGoldenStep(a, b);
  c = clamp(c, alphaLowerBound, alphaUpperBound);
  if (c == b)
    fc = fb;
  else
    fc = f1D(c);

  if (verbose)
    cout << "after golden step, we have a = " << a << " b = " << b << " c = " << c <<
      ", fa = " << fa << " fb = " << fb << " fc = " << fc << endl;

  if (c == b) // b alreay hit bounds, so c cannot be beyond b,
  {           // we have to find another value as the new b between a and c
    if (verbose)
      cout << "b already hit bounds, testing other situation" << endl;
    int numFallbackTryIter = 10;
    for(int iterID = 0; iterID < numFallbackTryIter; iterID++)
    {
      b = (a + c) / 2.0;
      fb = f1D(b);
      if (fb <= fa && fb <= fc)
        return 0;

      swap(a, b);
      swap(fa, fb);
    }
    return 1;
  }

  int iterID = -1;
  while(fc < fb) // while the condition f1D(b) <= min( f1D(a), f1D(c) ) on [a, b, c] are not satisfied
  {
    iterID++;
    if (iterID >= maxIter)
    {
      if (verbose)
        cout << "Max iteration " << maxIter << " on initializing triplet" << endl;
      return 1;
    }
    // fit a parabola on a, b, c and compute the vertex as d
    // the equation is:
    //     (b^2-c^2)(fb-fa) - (b^2-a^2)(fb-fc)
    // d = ------------------------------------
    //       2[(b-c)(fb-fa) - (b-a)(fb-fc)]
    double bcfbfa = (b-c) * (fb-fa);
    double bafbfc = (b-a) * (fb-fc);
    double numerator = (b+c) * bcfbfa - (b+a) * bafbfc;
    const double threshold = 1e-20;
    double denominator = clampFromZero(2 * (bcfbfa - bafbfc), threshold);
    double d = numerator / denominator;
    const double maxStepSizeRatio = 100.0;
    double dMax = b + (c-b) * maxStepSizeRatio; // furtherst d allowed away from b and c
    d = clamp(d, alphaLowerBound, alphaUpperBound);
    dMax = clamp(dMax, alphaLowerBound, alphaUpperBound);
    if (verbose)
      cout << "triplet iter " << iterID <<": " << a << " " << b << " " << c << " "
        << fa << " " << fb << " " << fc << ", d = " << d << " dMax = " << dMax << endl;

    // now we have a function graph like this:
    //f ^                               * END OF GRAPH
    //  | \                             *
    //  |  ----                         *
    //  |      \                        *
    //  |       ----                    *
    //  |           \                   *
    //  ----------------------> x or -x *
    //     a    b    c
    if (inBetween(d, b, c))
    {
      double fd = f1D(d);
      if (fd <= fc) // then [b, d, c] becomes one valid triplet
      {
        a = b;
        b = d;
        fa = fb;
        fb = fd;
        return 0;
      }
      else if (fd <= fb) // now fd is between fb and fc
      {
        // no valid triplet, we have to go to the next iteration
        d = takeGoldenStep(b, c); // create a new d
        // use [b, c, d] as the new [a, b, c] and continue
        a = b;
        b = c;
        c = d;
        fa = fb;
        fb = fc;
        fc = fd;
        continue;
      }
      else // fd > fb, then [a, b, d] becomes one valid triplet
      {
        c = d;
        fc = fd;
        return 0;
      }
    } // end if inBetween(d, b, c)
    else if (inBetween(d, c, dMax))
    {
      double fd = f1D(d);
      if (fd >= fc) // valid triplet: [b, c, d]
      {
        a = b;
        b = c;
        c = d;
        fa = fb;
        fb = fc;
        fc = fd;
        return 0;
      }
      else // fd < fc: no valid triplet, we have to go to the next iteration
      {
        // let's take [c, d, goldenStep(c,d)] as the new [a, b, c] in the next iteration
        a = c;
        b = d;
        fa = fc;
        fb = fd;
        c = takeGoldenStep(a, b);
        fc = f1D(c);
        continue;
      }
    } // end if inBetween(d, c, dMax)
    else if (inBetween(dMax, c, d)) // if d beyond dMax
    {
      // use dMax as d, continue
      a = b;
      b = c;
      c = dMax;
      fa = fb;
      fb = fc;
      fc = f1D(c);
      continue;
    }
    else // d is before b, then the fitted parabola is concave
    {    // it is a bad parabola, let's find a new value from goldenStep and
         // go to next iteration
      a = b;
      b = c;
      c = takeGoldenStep(a, b);
      fa = fb;
      fb = fc;
      fc = f1D(c);
      continue;
    }
    assert(0); // show never reach here
  } // end while loop
  return 0;
}

double LineSearch::f1D(double alpha) const
{
  f1DCount++;
  for(int i = 0; i < n; i++)
  {
    f1DBuffer[i] = x0[i] + p[i] * alpha;
  }
  double f = 0.0;
  int ret = func(f1DBuffer.data(), &f, nullptr);
  // cout << "compute f1D at alpha = " << setprecision(numeric_limits<long double>::digits10 + 1) << 
      // alpha << endl;
  if (ret != 0)
    f = DBL_MAX;
  return f;
}
