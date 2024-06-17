/*
author: Bohan Wang
copyright to USC, MIT
*/

#include "barrierFunction.h"

#include <algorithm>
#include <cmath>

namespace pgo
{
namespace NonlinearOptimization
{
namespace BarrierFunctions
{
double logBarrierEnergy(double v, double vhat, double positiveZero)
{
  v = std::max(v, positiveZero);

  double ratio = v / vhat;

  //-(z - c0)^2 ln(z/c0)
  return -(v - vhat) * (v - vhat) * std::log(ratio);
}

double logBarrierGradient(double v, double vhat, double positiveZero)
{
  v = std::max(v, positiveZero);

  //-(z - c0)^2 ln(z/c0)
  // dx = -2(z - c0) ln(z/c0) - (z - c0)^2 1/(z/c0) 1/c0
  double diff = v - vhat;
  double ratio = v / vhat;
  return -2.0 * diff * std::log(ratio) - diff * diff / v;
}

double logBarrierHessian(double v, double vhat, double positiveZero)
{
  v = std::max(v, positiveZero);

  //-(z - c0)^2 ln(z/c0)
  // dx = -2(z - c0) ln(z/c0) - (z - c0)^2 1/(z/c0) 1/c0
  // d2x = -2 ln (z/c0) - 2(z - c0) 1/z - 2(z- c0) 1/z - (z - c0)^2 (-1) z^-2
  double diff = v - vhat;
  double ratio = v / vhat;

  return -2.0 * std::log(ratio) - 4.0 * diff / v + diff * diff / (v * v);
}

double polyBarrierEnergy(double v, double vhat, int power)
{
  power = power * 2;
  return std::pow(v - vhat, power);
}

double polyBarrierGradient(double v, double vhat, int power)
{
  power = power * 2;
  return std::pow(v - vhat, power - 1) * power;
}

double polyBarrierHessian(double v, double vhat, int power)
{
  power = power * 2;
  return std::pow(v - vhat, power - 2) * power * (power - 1);
}
}  // namespace BarrierFunctions
}  // namespace NonlinearOptimization
}  // namespace pgo
