/*
author: Bohan Wang
copyright to USC
*/

#pragma once

#include "EigenSupport.h"

#include <memory>
#include <functional>

namespace pgo
{
namespace NonlinearOptimization
{
class PotentialEnergy;
typedef std::shared_ptr<PotentialEnergy> PotentialEnergy_p;
typedef std::shared_ptr<const PotentialEnergy> PotentialEnergy_const_p;

class ConstraintFunctions;
typedef std::shared_ptr<ConstraintFunctions> ConstraintFunctions_p;
typedef std::shared_ptr<const ConstraintFunctions> ConstraintFunctions_const_p;

class PotentialEnergyDense;
typedef std::shared_ptr<PotentialEnergyDense> PotentialEnergyDense_p;
typedef std::shared_ptr<const PotentialEnergyDense> PotentialEnergyDense_const_p;

class ConstraintFunctionsDense;
typedef std::shared_ptr<ConstraintFunctionsDense> ConstraintFunctionsDense_p;
typedef std::shared_ptr<const ConstraintFunctionsDense> ConstraintFunctionsDense_const_p;

namespace EnergyOptimizer
{
enum class SolverType
{
  ST_IPOPT,
  ST_KNITRO,
  ST_NEWTON
};

int minimize(EigenSupport::RefVecXd x, PotentialEnergy_const_p energy, EigenSupport::ConstRefVecXd xlow, EigenSupport::ConstRefVecXd xhi,
  EigenSupport::RefVecXd lambda, EigenSupport::RefVecXd g, ConstraintFunctions_const_p constraints,
  EigenSupport::ConstRefVecXd clow, EigenSupport::ConstRefVecXd chi,
  SolverType solverType, int maxIter, double eps, int verbose);

int minimize(EigenSupport::RefVecXd x, PotentialEnergy_const_p energy, EigenSupport::ConstRefVecXd xlow, EigenSupport::ConstRefVecXd xhi,
  SolverType solverType, int maxIter, double eps, int verbose);

typedef std::function<void(const double *, int, const double *, int)> CallbackFunc;

int minimizeUsingIpopt(EigenSupport::RefVecXd x, PotentialEnergy_const_p energy, EigenSupport::ConstRefVecXd xlow, EigenSupport::ConstRefVecXd xhi,
  EigenSupport::RefVecXd lambda, EigenSupport::RefVecXd g, ConstraintFunctions_const_p constraints, EigenSupport::ConstRefVecXd clow, EigenSupport::ConstRefVecXd chi,
  int maxIter, double eps, int verbose);

int minimizeUsingNewton(EigenSupport::RefVecXd x, PotentialEnergy_const_p energy, EigenSupport::ConstRefVecXd xlow, EigenSupport::ConstRefVecXd xhi,
  EigenSupport::RefVecXd lambda, EigenSupport::RefVecXd g, ConstraintFunctions_const_p constraints, EigenSupport::ConstRefVecXd clow, EigenSupport::ConstRefVecXd chi,
  int maxIter, double eps, int verbose);

struct KnitroData;
int minimizeUsingKnitro(EigenSupport::RefVecXd x, PotentialEnergy_const_p energy, EigenSupport::ConstRefVecXd xlow, EigenSupport::ConstRefVecXd xhi,
  EigenSupport::RefVecXd lambda, EigenSupport::RefVecXd g, ConstraintFunctions_const_p constraints, EigenSupport::ConstRefVecXd clow, EigenSupport::ConstRefVecXd chi,
  int maxIter, double eps, int verbose, const char *configFilename, int parallelEval = 0, CallbackFunc callback = nullptr, KnitroData **data = nullptr, double feasTol = -1.0);

int minimizeUsingKnitroDense(EigenSupport::RefVecXd x, PotentialEnergyDense_const_p energy, EigenSupport::ConstRefVecXd xlow, EigenSupport::ConstRefVecXd xhi,
  EigenSupport::RefVecXd lambda, EigenSupport::RefVecXd g, ConstraintFunctionsDense_const_p constraints, EigenSupport::ConstRefVecXd clow, EigenSupport::ConstRefVecXd chi,
  int maxIter, double eps, int verbose, const char *configFilename, int parallelEval = 0, CallbackFunc callback = nullptr, KnitroData **data = nullptr, double feasTol = -1.0);

int minimizeUsingApproximateActiveSet(EigenSupport::RefVecXd x, PotentialEnergy_const_p energy, EigenSupport::ConstRefVecXd xlow, EigenSupport::ConstRefVecXd xhi,
  EigenSupport::RefVecXd lambda, EigenSupport::RefVecXd g, ConstraintFunctions_const_p constraints, EigenSupport::ConstRefVecXd clow, EigenSupport::ConstRefVecXd chi,
  int maxIter, double eps, double equalityThreshold, int verbose);
}  // namespace EnergyOptimizer

}  // namespace NonlinearOptimization
}  // namespace pgo