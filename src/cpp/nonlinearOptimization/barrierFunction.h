/*
author: Bohan Wang
copyright to USC, MIT
*/

#pragma once

namespace pgo
{
namespace NonlinearOptimization
{
namespace BarrierFunctions
{
double logBarrierEnergy(double v, double vhat, double positiveZero = 1e-10);
double logBarrierGradient(double v, double vhat, double positiveZero = 1e-10);
double logBarrierHessian(double v, double vhat, double positiveZero = 1e-10);

double polyBarrierEnergy(double v, double vhat, int power);
double polyBarrierGradient(double v, double vhat, int power);
double polyBarrierHessian(double v, double vhat, int power);
}  // namespace BarrierFunctions
}  // namespace NonlinearOptimization
}  // namespace pgo