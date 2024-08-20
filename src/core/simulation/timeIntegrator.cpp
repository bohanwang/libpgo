#include "timeIntegrator.h"
#include "finiteDifference.h"

#include <numeric>
#include <cstring>
#include <iostream>

using namespace pgo;
using namespace pgo::NonlinearOptimization;
using namespace pgo::ConstraintPotentialEnergies;
using namespace pgo::Simulation;

namespace ES = EigenSupport;

TimeIntegrator::TimeIntegrator(const ES::SpMatD &massMatrix, std::shared_ptr<const PotentialEnergy> elasticPotential,
  double massDampingCoeff, double stiffnessDampingCoeff, double t, int numIter, double eps_):
  elasticEnergy(elasticPotential),
  timestep(t),
  dK(stiffnessDampingCoeff), dM(massDampingCoeff), nIter(numIter), eps(eps_)
{
  n3 = elasticEnergy->getNumDOFs();

  elasticEnergy->createHessian(K);
  elasticEnergy->createHessian(K1);
  elasticEnergy->createHessian(MasK);

  memset(MasK.valuePtr(), 0, sizeof(double) * MasK.nonZeros());

  ES::SpMatI M2K;
  ES::small2Big(massMatrix, K, 0, 0, M2K);
  ES::transferSmallToBig(massMatrix, MasK, M2K, 1);

  nnz = (int)K.nonZeros();

  q.setZero(n3);
  qvel.setZero(n3);
  qacc.setZero(n3);

  q1.setZero(n3);
  qvel1.setZero(n3);
  qacc1.setZero(n3);

  f_int.setZero(n3);
  f_ext.setZero(n3);

  allDOFs.resize(n3);
  std::iota(allDOFs.begin(), allDOFs.end(), 0);

  deltauInitial.setZero(n3);

  deltauRangeLow.setZero(n3);
  deltauRangeHi.setZero(n3);

  for (int i = 0; i < n3; i++) {
    deltauRangeLow[i] = -defaultUnknownBoundary;
    deltauRangeHi[i] = defaultUnknownBoundary;
  }

  hessianAll = K;
  ES::small2Big(K, hessianAll, allDOFs, Kmapping);
  rhss2b = allDOFs;
  rhsb2s = allDOFs;
}

TimeIntegrator::TimeIntegrator::~TimeIntegrator()
{
}

void TimeIntegrator::TimeIntegrator::addExternalForce(int vid, const double f[3])
{
  f_ext.segment<3>(vid * 3) += ES::V3d(f[0], f[1], f[2]);
}

void TimeIntegrator::TimeIntegrator::addExternalForce(const double *f)
{
  f_ext += ES::Mp<const ES::VXd>(f, n3);
}

void TimeIntegrator::TimeIntegrator::setExternalForce(const double *f)
{
  f_ext.noalias() = ES::Mp<const ES::VXd>(f, n3);
}

void TimeIntegrator::clearExternalForce()
{
  f_ext.setZero();
}

void TimeIntegrator::getExternalForce(double *f) const
{
  (ES::Mp<ES::VXd>(f, n3)) = f_ext;
}

void TimeIntegrator::proceedTimestep()
{
  q.noalias() = q1;
  qvel.noalias() = qvel1;
  qacc.noalias() = qacc1;

  timestepID++;
}

void TimeIntegrator::addImplicitForceModel(std::shared_ptr<PotentialEnergyAligningMeshConnectivity> fm, double kd, double md)
{
  additionalForceModels.push_back(fm);
  additionalForceModels_K.push_back(K);
  additionalForceModels_K1.push_back(K);
  additionalForceModels_fint.emplace_back(n3);

  if (kd < 0)
    additionalForceModelsDampingParams.push_back(dK);
  else
    additionalForceModelsDampingParams.push_back(kd);

  additionalForceModelsMassDampingParams.push_back(md);
}

void TimeIntegrator::clearImplicitForceModel()
{
  additionalForceModels.clear();
  additionalForceModels_fint.clear();
  additionalForceModels_K.clear();
  additionalForceModels_K1.clear();
  additionalForceModelsDampingParams.clear();
  additionalForceModelsMassDampingParams.clear();
}

void TimeIntegrator::addGeneralImplicitForceModel(std::shared_ptr<PotentialEnergy> fm, double kd, double md)
{
  generalAdditionalForceModels.push_back(fm);

  generalAdditionalForceModels_K.push_back(ES::SpMatD());
  fm->createHessian(generalAdditionalForceModels_K.back());
  generalAdditionalForceModels_K1.push_back(generalAdditionalForceModels_K.back());
  generalAdditionalForceModels_M.push_back(generalAdditionalForceModels_K.back());

  for (Eigen::Index outeri = 0; outeri < generalAdditionalForceModels_M.back().outerSize(); outeri++) {
    for (ES::SpMatD::InnerIterator it(generalAdditionalForceModels_M.back(), outeri); it; ++it) {
      Eigen::Index row = it.row();
      Eigen::Index col = it.col();

      if (MasK.nonZeros()) {
        std::ptrdiff_t offset = ES::findEntryOffset(MasK, (int)row, (int)col);
        if (offset >= 0)
          it.valueRef() = MasK.valuePtr()[offset];
        else
          it.valueRef() = 0;
      }
      else
        it.valueRef() = 0;

      // it.valueRef() = MasK.coeff(row, col);
    }
  }

  generalAdditionalForceModels_fint.emplace_back(n3);

  if (kd < 0)
    generalForceModelsDampingParams.push_back(dK);
  else
    generalForceModelsDampingParams.push_back(kd);
  generalForceModelsMassDampingParams.push_back(md);

  generalForceModelChanged = true;
}

void TimeIntegrator::clearGeneralImplicitForceModel()
{
  if (generalAdditionalForceModels.size()) {
    generalForceModelChanged = true;
  }

  generalAdditionalForceModels.clear();
  generalForceModelsDampingParams.clear();
  generalForceModelsMassDampingParams.clear();
  generalAdditionalForceModels_K.clear();
  generalAdditionalForceModels_K1.clear();
  generalAdditionalForceModels_M.clear();
  generalAdditionalForceModels_fint.clear();
}

void TimeIntegrator::assembleImplicitModels()
{
  implicitModelsAll.clear();
  implicitModelTypesAll.clear();
  dampingParamsAll.clear();
  massDampingParamsAll.clear();

  implicitModelsAll_K.clear();
  implicitModelsAll_K1.clear();
  implicitModelsAll_M.clear();
  implicitModelsAll_fint.clear();
  implicitModelsAll_Kmaping.clear();

  implicitModelsAll.push_back(elasticEnergy);
  implicitModelTypesAll.push_back(IMT_ELASTIC);
  dampingParamsAll.push_back(dK);
  massDampingParamsAll.push_back(dM);

  implicitModelsAll_K.push_back(&K);
  implicitModelsAll_K1.push_back(&K1);
  implicitModelsAll_M.push_back(&MasK);
  implicitModelsAll_fint.push_back(&f_int);
  implicitModelsAll_Kmaping.push_back(&Kmapping);

  int inc = 0;
  for (auto implicit : additionalForceModels) {
    implicitModelsAll.push_back(implicit);
    implicitModelTypesAll.push_back(IMT_SAME_TOPOLOGY);
    dampingParamsAll.push_back(additionalForceModelsDampingParams[inc]);
    massDampingParamsAll.push_back(additionalForceModelsMassDampingParams[inc]);

    implicitModelsAll_K.push_back(&additionalForceModels_K[inc]);
    implicitModelsAll_K1.push_back(&additionalForceModels_K1[inc]);
    implicitModelsAll_M.push_back(&MasK);
    implicitModelsAll_fint.push_back(&additionalForceModels_fint[inc]);
    implicitModelsAll_Kmaping.push_back(&Kmapping);

    inc++;
  }

  if (generalForceModelChanged) {
    generalAdditionalForceModels_Kmapping.resize(generalAdditionalForceModels.size());

    inc = 0;
    for (auto implicit : generalAdditionalForceModels) {
      implicitModelsAll.push_back(implicit);
      implicitModelTypesAll.push_back(IMT_GENERAL);
      dampingParamsAll.push_back(generalForceModelsDampingParams[inc]);
      massDampingParamsAll.push_back(generalForceModelsMassDampingParams[inc]);

      implicitModelsAll_K.push_back(&generalAdditionalForceModels_K[inc]);
      implicitModelsAll_K1.push_back(&generalAdditionalForceModels_K1[inc]);
      implicitModelsAll_M.push_back(&generalAdditionalForceModels_M[inc]);

      implicitModelsAll_fint.push_back(&generalAdditionalForceModels_fint[inc]);
      implicitModelsAll_Kmaping.push_back(&generalAdditionalForceModels_Kmapping[inc]);

      inc++;
    }

    // generate mapping
    if (inc) {
      std::vector<ES::TripletD> entries;
      entries.reserve(K.nonZeros() * 2);
      for (Eigen::Index outeri = 0; outeri < K.outerSize(); outeri++) {
        for (ES::SpMatD::InnerIterator it(K, outeri); it; ++it) {
          entries.emplace_back(
            (ES::SpMatD::StorageIndex)it.row(),
            (ES::SpMatD::StorageIndex)it.col(),
            1.0);
        }
      }

      for (size_t i = 0; i < generalAdditionalForceModels.size(); i++) {
        for (Eigen::Index outeri = 0; outeri < generalAdditionalForceModels_K[i].outerSize(); outeri++) {
          for (ES::SpMatD::InnerIterator it(generalAdditionalForceModels_K[i], outeri); it; ++it) {
            entries.emplace_back(
              (ES::SpMatD::StorageIndex)it.row(),
              (ES::SpMatD::StorageIndex)it.col(),
              1.0);
          }
        }
      }

      hessianAll.resize(n3, n3);
      hessianAll.setFromTriplets(entries.begin(), entries.end());

      for (size_t i = 0; i < generalAdditionalForceModels.size(); i++) {
        ES::small2Big(generalAdditionalForceModels_K[i], hessianAll, allDOFs, generalAdditionalForceModels_Kmapping[i]);
      }
    }
    else {
      hessianAll = K;
    }

    ES::small2Big(K, hessianAll, allDOFs, Kmapping);
  }
  else if (generalAdditionalForceModels.size()) {
    inc = 0;
    for (auto implicit : generalAdditionalForceModels) {
      implicitModelsAll.push_back(implicit);
      implicitModelTypesAll.push_back(IMT_GENERAL);
      dampingParamsAll.push_back(generalForceModelsDampingParams[inc]);
      massDampingParamsAll.push_back(generalForceModelsMassDampingParams[inc]);

      implicitModelsAll_K.push_back(&generalAdditionalForceModels_K[inc]);
      implicitModelsAll_K1.push_back(&generalAdditionalForceModels_K1[inc]);
      implicitModelsAll_M.push_back(&generalAdditionalForceModels_M[inc]);

      implicitModelsAll_fint.push_back(&generalAdditionalForceModels_fint[inc]);
      implicitModelsAll_Kmaping.push_back(&generalAdditionalForceModels_Kmapping[inc]);

      inc++;
    }
  }
}

void TimeIntegrator::clearDeltauInitial()
{
  memset(deltauInitial.data(), 0, sizeof(double) * n3);
}

void TimeIntegrator::setDeltauRange(ES::ConstRefVecXd low, ES::ConstRefVecXd hi)
{
  deltauRangeLow.noalias() = low;
  deltauRangeHi.noalias() = hi;

  updateFixedDOFs();
}

void TimeIntegrator::setDeltauRange(double delta)
{
  for (int i = 0; i < n3; i++) {
    deltauRangeLow[i] = -delta;
    deltauRangeHi[i] = delta;
  }

  updateFixedDOFs();

  std::cout << "delta u range:" << deltauRangeLow[0] << ',' << deltauRangeHi[1] << std::endl;
}

void TimeIntegrator::clearDeltauRange()
{
  for (int i = 0; i < n3; i++) {
    deltauRangeLow[i] = -defaultUnknownBoundary;
    deltauRangeHi[i] = defaultUnknownBoundary;
  }

  updateFixedDOFs();
}

void TimeIntegrator::setFixedVertices(const std::vector<int> &fixedVertices, ES::ConstRefVecXd fixedRestP, ES::ConstRefVecXd fixedP, int isIndexSorted)
{
  fixedDOFs.clear();

  if (isIndexSorted) {
    for (int v : fixedVertices) {
      fixedDOFs.push_back(v * 3);
      fixedDOFs.push_back(v * 3 + 1);
      fixedDOFs.push_back(v * 3 + 2);
    }

    fixedPosition = fixedP;
    fixedRestPosition = fixedRestP;
  }
  else {
    typedef std::pair<int, std::array<double, 6>> FixedInfo;
    std::vector<FixedInfo> fixedInfo;
    fixedInfo.reserve(fixedVertices.size());

    for (size_t i = 0; i < fixedVertices.size(); i++) {
      std::array<double, 6> pos = {
        fixedRestP[i * 3], fixedRestP[i * 3 + 1], fixedRestP[i * 3 + 2],
        fixedP[i * 3], fixedP[i * 3 + 1], fixedP[i * 3 + 2]
      };

      fixedInfo.emplace_back(fixedVertices[i], pos);
    }

    std::sort(fixedInfo.begin(), fixedInfo.end(), [](const FixedInfo &fi1, const FixedInfo &fi2) {
      return fi1.first < fi2.first;
    });

    fixedPosition.resize(fixedInfo.size() * 3);
    fixedRestPosition.resize(fixedInfo.size() * 3);

    for (size_t i = 0; i < fixedInfo.size(); i++) {
      int v = fixedInfo[i].first;
      fixedDOFs.push_back(v * 3);
      fixedDOFs.push_back(v * 3 + 1);
      fixedDOFs.push_back(v * 3 + 2);

      const auto &vec = fixedInfo[i].second;
      fixedRestPosition.segment<3>(i * 3) = ES::V3d(vec[0], vec[1], vec[2]);
      fixedPosition.segment<3>(i * 3) = ES::V3d(vec[3], vec[4], vec[5]);
    }
  }

  rhsb2s.clear();
  rhss2b.clear();
  ES::removeRows(n3, fixedDOFs, rhsb2s, rhss2b);

  updateFixedDOFs();
}

void TimeIntegrator::setFixedVertices(const std::vector<int> &fixedVertices)
{
  fixedDOFs.clear();

  for (int v : fixedVertices) {
    fixedDOFs.push_back(v * 3);
    fixedDOFs.push_back(v * 3 + 1);
    fixedDOFs.push_back(v * 3 + 2);
  }

  fixedPosition = ES::VXd::Zero(fixedDOFs.size());
  fixedRestPosition = ES::VXd::Zero(fixedDOFs.size());

  rhsb2s.clear();
  rhss2b.clear();
  ES::removeRows(n3, fixedDOFs, rhsb2s, rhss2b);

  updateFixedDOFs();
}

void TimeIntegrator::updateFixedDOFs()
{
  for (size_t i = 0; i < fixedDOFs.size(); i++) {
    double targetp = fixedPosition[i];
    double restp = fixedRestPosition[i];
    double curu = q[fixedDOFs[i]];

    double deltau = targetp - curu - restp;

    deltauRangeLow[fixedDOFs[i]] = deltau;
    deltauRangeHi[fixedDOFs[i]] = deltau;
  }
}

void TimeIntegrator::addConstraints(std::shared_ptr<ConstraintFunctions> constt, const int *equalities)
{
  constraints = constt;
  constraintsRangeLow = ES::VXd::Zero(constt->getNumConstraints());
  constraintsRangeHi = ES::VXd::Zero(constt->getNumConstraints());
  g = ES::VXd::Zero(constt->getNumConstraints());
  lambda = ES::VXd::Zero(constt->getNumConstraints());

  if (equalities) {
    for (int i = 0; i < constraints->getNumConstraints(); i++) {
      if (equalities[i] > 0) {
        constraintsRangeLow[i] = 0;
        constraintsRangeHi[i] = defaultUnknownBoundary;
      }
      else if (equalities[i] < 0) {
        constraintsRangeLow[i] = -defaultUnknownBoundary;
        constraintsRangeHi[i] = 0;
      }
    }
  }

  constraintsChanged = true;
}

void TimeIntegrator::addConstraints(std::shared_ptr<ConstraintFunctions> constt, int eq)
{
  constraints = constt;
  constraintsRangeLow = ES::VXd::Zero(constt->getNumConstraints());
  constraintsRangeHi = ES::VXd::Zero(constt->getNumConstraints());
  g = ES::VXd::Zero(constt->getNumConstraints());
  lambda = ES::VXd::Zero(constt->getNumConstraints());

  if (eq != 0) {
    for (int i = 0; i < constraints->getNumConstraints(); i++) {
      if (eq > 0) {
        constraintsRangeLow[i] = 0;
        constraintsRangeHi[i] = defaultUnknownBoundary;
      }
      else if (eq < 0) {
        constraintsRangeLow[i] = -defaultUnknownBoundary;
        constraintsRangeHi[i] = 0;
      }
    }
  }
  constraintsChanged = true;
}

void TimeIntegrator::addConstraints(std::shared_ptr<ConstraintFunctions> constt, ES::ConstRefVecXd clow, ES::ConstRefVecXd chi)
{
  constraints = constt;

  constraintsRangeLow = clow;
  constraintsRangeHi = chi;

  g = ES::VXd::Zero(constt->getNumConstraints());
  lambda = ES::VXd::Zero(constt->getNumConstraints());

  constraintsChanged = true;
}

void TimeIntegrator::clearConstraints()
{
  constraints.reset();
  constraintsChanged = true;
}

void TimeIntegrator::doTimestep(int /*updateq*/, int /*verbose*/, int /*printResidual*/)
{
  generalForceModelChanged = false;
  constraintsChanged = false;
}

void TimeIntegrator::enableFiniteDifferenceTest(bool testIntegratorEnergy, bool testImplicitEnergy, bool testConstraints)
{
  finiteDifferenceTestFlag = 0;
  finiteDifferenceTestFlag |= (testIntegratorEnergy ? 1 : 0);
  finiteDifferenceTestFlag |= (testImplicitEnergy ? 2 : 0);
  finiteDifferenceTestFlag |= (testConstraints ? 4 : 0);
}

void TimeIntegrator::finiteDifferenceTestImplicitEnergy(ES::ConstRefVecXd x) const
{
  FiniteDifference fd(FiniteDifference::M_FIVE_POINT, 1e-7);

  for (size_t i = 0; i < implicitModelsAll.size(); i++) {
    std::cout << "Test energy " << i << std::endl;
    fd.testEnergy(implicitModelsAll[i], true, true, -1, x.data());
  }
}

void TimeIntegrator::finiteDifferenceTestConstraints(ES::ConstRefVecXd x) const
{
  std::cout << "Test current constraints:\n";
  FiniteDifference fd(FiniteDifference::M_FIVE_POINT, 1e-7);

  fd.testConstraints(constraints, -1, x.data());
}

void TimeIntegrator::finiteDifferenceTest(ES::ConstRefVecXd x, ES::ConstRefVecXd u) const
{
  if (finiteDifferenceTestFlag & 1)
    finiteDifferenceTestIntegratorEnergy(x);

  if (finiteDifferenceTestFlag & 2)
    finiteDifferenceTestImplicitEnergy(u);

  if (finiteDifferenceTestFlag & 4)
    finiteDifferenceTestConstraints(x);
}

void TimeIntegrator::finiteDifferenceTest(double range)
{
  ES::VXd u(n3), x(n3);
  FiniteDifference::randomSeq(range, u.data(), n3);
  FiniteDifference::randomSeq(range, x.data(), n3);

  std::cout << "||u||=" << u.norm() << std::endl;
  std::cout << "||x||=" << x.norm() << std::endl;

  finiteDifferenceTest(x, u);
}

void TimeIntegrator::computeInternalForces(const double *u, double *fint) const
{
  ES::Mp<ES::VXd> fintMp(fint, n3);
  ES::Mp<const ES::VXd> uMp(u, n3);

  fintMp.setZero();
  for (const auto &model : implicitModelsAll) {
    ES::VXd grad(n3);
    grad.setZero();

    model->gradient(uMp, grad);
    fintMp += grad;
  }
}

double TimeIntegrator::defaultUnknownBoundary = 1e20;