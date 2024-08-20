#pragma once

#include "timeIntegratorOptions.h"

#include "EigenSupport.h"
#include "potentialEnergy.h"
#include "potentialEnergies.h"
#include "potentialEnergyAligningMeshConnectivity.h"
#include "constraintFunctions.h"

namespace pgo
{
namespace Simulation
{
class TimeIntegrator
{
public:
  TimeIntegrator(const EigenSupport::SpMatD &massMatrix,
    std::shared_ptr<const NonlinearOptimization::PotentialEnergy> elasticPotential,
    double massDampingCoeff, double stiffnessDampingCoeff, double timestep, int numIter = 0, double eps = 1e-5);
  virtual ~TimeIntegrator();

  void addExternalForce(const double *fext);
  void addExternalForce(int vid, const double fext[3]);
  void setExternalForce(const double *fext);
  void getExternalForce(double *fext) const;
  void clearExternalForce();

  virtual void doTimestep(int updateq = 1, int verbose = 0, int printResidual = 0);
  void proceedTimestep();

  virtual void setTimestep(double t) { timestep = t; }
  void setMaxIter(int niter) { nIter = niter; }
  void setTol(double tol) { eps = tol; }
  double getTimestep() const { return timestep; }

  void setDampingParam(double a, double b) { dM = a, dK = b; }

  void getq(EigenSupport::RefVecXd u) const { u.noalias() = q; }
  void getqvel(EigenSupport::RefVecXd uvel) const { uvel.noalias() = qvel; }
  void getqacc(EigenSupport::RefVecXd uacc) const { uacc.noalias() = qacc; }

  void getq1(EigenSupport::RefVecXd u) const { u.noalias() = q1; }
  void getqvel1(EigenSupport::RefVecXd uvel) const { uvel.noalias() = qvel1; }
  void getqacc1(EigenSupport::RefVecXd uacc) const { uacc.noalias() = qacc1; }

  void setqState(EigenSupport::ConstRefVecXd u, EigenSupport::ConstRefVecXd uvel, EigenSupport::ConstRefVecXd uacc)
  {
    q.noalias() = u;
    qvel.noalias() = uvel;
    qacc.noalias() = uacc;
  }

  void setSolverOption(TimeIntegratorSolverOption so) { solverOption = so; }
  TimeIntegratorSolverOption getSolverOption() const { return solverOption; }
  int getSolverReturn() const { return solverRet; }

  void addImplicitForceModel(std::shared_ptr<ConstraintPotentialEnergies::PotentialEnergyAligningMeshConnectivity> fm, double kd = -1, double md = 0);
  void clearImplicitForceModel();

  void addGeneralImplicitForceModel(std::shared_ptr<NonlinearOptimization::PotentialEnergy> fm, double kd = -1, double md = 0);
  void clearGeneralImplicitForceModel();

  virtual void addConstraints(std::shared_ptr<NonlinearOptimization::ConstraintFunctions> constt, int sign = 0);
  virtual void addConstraints(std::shared_ptr<NonlinearOptimization::ConstraintFunctions> constt, const int *equalities = nullptr);
  virtual void addConstraints(std::shared_ptr<NonlinearOptimization::ConstraintFunctions> constt, EigenSupport::ConstRefVecXd clow, EigenSupport::ConstRefVecXd chi);
  virtual void clearConstraints();

  void getFinalContraint(EigenSupport::RefVecXd g_) const { g_ = g; }
  void getFinalContraintLambda(EigenSupport::RefVecXd l_) const { l_ = lambda; }

  void setDeltauInitial(EigenSupport::ConstRefVecXd deltau) { deltauInitial.noalias() = deltau; }
  void clearDeltauInitial();

  void setDeltauRange(double delta);
  void setDeltauRange(EigenSupport::ConstRefVecXd low, EigenSupport::ConstRefVecXd hi);
  void clearDeltauRange();

  void setFixedVertices(const std::vector<int> &fixedVertices, EigenSupport::ConstRefVecXd fixedRestPosition, EigenSupport::ConstRefVecXd fixedPosition, int isIndexSorted = 1);
  void setFixedVertices(const std::vector<int> &fixedVertices);

  void enableFiniteDifferenceTest(bool testIntegratorEnergy, bool testImplicitEnergy, bool testConstraints);
  void finiteDifferenceTest(double range);

  void setSolverConfigFile(const char *filename) { solverConfigFilename = filename; }

  void resetTimestepID() { timestepID = 0; }
  void setTimestepID(int id) { timestepID = id; }

  virtual std::shared_ptr<const NonlinearOptimization::PotentialEnergy> getInternalEnergy() const = 0;
  virtual void computeInternalForces(const double *u, double *fint) const;

protected:
  void updateFixedDOFs();
  void assembleImplicitModels();

  void finiteDifferenceTestImplicitEnergy(EigenSupport::ConstRefVecXd x) const;
  void finiteDifferenceTestConstraints(EigenSupport::ConstRefVecXd x) const;
  virtual void finiteDifferenceTestIntegratorEnergy(EigenSupport::ConstRefVecXd x) const = 0;
  void finiteDifferenceTest(EigenSupport::ConstRefVecXd x, EigenSupport::ConstRefVecXd u) const;

  EigenSupport::SpMatD MasK, K, K1;
  EigenSupport::SpMatI Kmapping;
  EigenSupport::SpMatD hessianAll;
  std::shared_ptr<const NonlinearOptimization::PotentialEnergy> elasticEnergy;
  double timestep;
  double dK, dM;
  int nIter;
  double eps;

  std::vector<int> fixedDOFs, allDOFs;

  std::vector<std::shared_ptr<ConstraintPotentialEnergies::PotentialEnergyAligningMeshConnectivity>> additionalForceModels;
  std::vector<double> additionalForceModelsDampingParams;
  std::vector<double> additionalForceModelsMassDampingParams;
  std::vector<EigenSupport::SpMatD> additionalForceModels_K;
  std::vector<EigenSupport::SpMatD> additionalForceModels_K1;
  std::vector<EigenSupport::VXd> additionalForceModels_fint;

  std::vector<std::shared_ptr<NonlinearOptimization::PotentialEnergy>> generalAdditionalForceModels;
  std::vector<double> generalForceModelsDampingParams;
  std::vector<double> generalForceModelsMassDampingParams;
  std::vector<EigenSupport::SpMatD> generalAdditionalForceModels_K;
  std::vector<EigenSupport::SpMatD> generalAdditionalForceModels_K1;
  std::vector<EigenSupport::SpMatD> generalAdditionalForceModels_M;
  std::vector<EigenSupport::SpMatI> generalAdditionalForceModels_Kmapping;
  std::vector<EigenSupport::VXd> generalAdditionalForceModels_fint;
  bool generalForceModelChanged = false;

  std::vector<std::shared_ptr<const NonlinearOptimization::PotentialEnergy>> implicitModelsAll;
  std::vector<double> dampingParamsAll;
  std::vector<double> massDampingParamsAll;
  std::vector<EigenSupport::SpMatD *> implicitModelsAll_K;
  std::vector<EigenSupport::SpMatD *> implicitModelsAll_K1;
  std::vector<EigenSupport::SpMatD *> implicitModelsAll_M;
  std::vector<EigenSupport::SpMatI *> implicitModelsAll_Kmaping;
  std::vector<EigenSupport::VXd *> implicitModelsAll_fint;

  enum ImplicitModelType
  {
    IMT_ELASTIC,
    IMT_SAME_TOPOLOGY,
    IMT_GENERAL
  };
  std::vector<ImplicitModelType> implicitModelTypesAll;

  std::shared_ptr<NonlinearOptimization::ConstraintFunctions> constraints;
  bool constraintsChanged = false;

  EigenSupport::VXd q, qvel, qacc;
  EigenSupport::VXd q1, qvel1, qacc1;
  EigenSupport::VXd g, lambda;

  EigenSupport::VXd f_int, f_ext;

  EigenSupport::VXd deltauInitial;
  EigenSupport::VXd deltauRangeLow, deltauRangeHi;
  EigenSupport::VXd constraintsRangeLow, constraintsRangeHi;

  std::vector<int> rhss2b, rhsb2s;

  int n3;
  int nnz;

  int solverRet = 0;
  uint64_t timestepID = 0;
  TimeIntegratorSolverOption solverOption = TimeIntegratorSolverOption::SO_NEWTON;
  EigenSupport::VXd fixedPosition, fixedRestPosition;

  int finiteDifferenceTestFlag = 0;
  std::string solverConfigFilename;
  static double defaultUnknownBoundary;
};
}  // namespace Simulation
}  // namespace pgo
