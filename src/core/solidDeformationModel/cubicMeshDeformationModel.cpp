#include "cubicMeshDeformationModel.h"

#include "elasticModel3DDeformationGradient.h"
#include "plasticModel3DDeformationGradient.h"

#include "EigenSupport.h"

#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace ES = pgo::EigenSupport;

namespace pgo
{
namespace SolidDeformationModel
{
namespace
{
using M3x8d = Eigen::Matrix<double, 3, 8>;
using M9x24d = Eigen::Matrix<double, 9, 24>;

constexpr int kNumVertices = 8;
constexpr int kNumIntegrationPoints = 8;

const int kVertexAlpha[kNumVertices] = { 0, 1, 1, 0, 0, 1, 1, 0 };
const int kVertexBeta[kNumVertices] = { 0, 0, 1, 1, 0, 0, 1, 1 };
const int kVertexGamma[kNumVertices] = { 0, 0, 0, 0, 1, 1, 1, 1 };

const double *optionalParamData(const ES::VXd &param)
{
  return param.size() ? param.data() : nullptr;
}
}  // namespace

class CubicMeshDeformationModelInternal
{
public:
  std::array<ES::V3d, kNumVertices> restX;
  std::array<M3x8d, kNumIntegrationPoints> dN_dabc;
  std::array<ES::M3d, kNumIntegrationPoints> restDmInv;
  std::array<M3x8d, kNumIntegrationPoints> restBm;
  std::array<M9x24d, kNumIntegrationPoints> rest_dFdx;
  std::array<double, kNumIntegrationPoints> weightDetJ;

  const ElasticModel3DDeformationGradient *elasticModel = nullptr;
  const PlasticModel3DDeformationGradient *plasticModel = nullptr;

  static void computeShapeDerivative(double alpha, double beta, double gamma, M3x8d &dN);
  static void computeSVD(const ES::M3d &Fe, ES::M3d &U, ES::M3d &V, ES::V3d &S);
  static void compute_dF_dx(const M3x8d &dN_dX, M9x24d &dFdx);
  static void computeCurrent_dFdx(const M9x24d &rest_dFdx, const ES::M3d &FpInv, M9x24d &dFdx);
  static void compute_d2Fe_dx_dai(const ES::M3d &dAInvdai, const M9x24d &rest_dFdx, M9x24d &d2Fdudai);
  static void compute_dP_dai(const ES::M9d &dPdF, const ES::M3d &dFdai, ES::M3d &dPdai);

  double compute_dV_dai(double weight, double ddetA_dai) const;
  double compute_d2V_daidaj(double weight, double d2detA_daidaj) const;
  void compute_dFe_dai(const ES::M3d &Fref, const ES::M3d &dAInvdai, ES::M3d &dFdai) const;
  void compute_d2Fe_dai_daj(const ES::M3d &Fref, const ES::M3d &dAInvdaidaj, ES::M3d &d2Fdaidaj) const;
  double compute_dpsi_dai(const ES::M3d &Fref, const ES::M3d &dAInv_dai, const ES::M3d &P) const;
  double compute_d2psi_dai_daj(const ES::M3d &Fref, const ES::M3d &dAInv_dai, const ES::M3d &dAInv_daj,
    const ES::M3d &d2AInv_dai_daj, const ES::M3d &P, const ES::M9d &dPdF) const;
};

class CubicMeshDeformationModelCacheData : public DeformationModelCacheData
{
public:
  CubicMeshDeformationModelCacheData(int numPlasticParams_, int numElasticParams_):
    numPlasticParams(numPlasticParams_),
    plasticParam(ES::VXd::Zero(numPlasticParams_)),
    ddetA_da(ES::VXd::Zero(numPlasticParams_)),
    d2detA_da2(ES::MXd::Zero(numPlasticParams_, numPlasticParams_)),
    dAInv_dai(numPlasticParams_, ES::M3d::Zero()),
    d2AInv_dai_daj(numPlasticParams_ * numPlasticParams_, ES::M3d::Zero()),
    materialParam(ES::VXd::Zero(numElasticParams_))
  {
  }

  ES::M3d &d2AInv(int i, int j)
  {
    return d2AInv_dai_daj[i * numPlasticParams + j];
  }

  const ES::M3d &d2AInv(int i, int j) const
  {
    return d2AInv_dai_daj[i * numPlasticParams + j];
  }

  int numPlasticParams = 0;

  ES::VXd plasticParam;
  ES::M3d Fp = ES::M3d::Identity();
  ES::M3d FpInv = ES::M3d::Identity();
  double detFp = 1.0;

  ES::VXd ddetA_da;
  ES::MXd d2detA_da2;

  std::vector<ES::M3d> dAInv_dai;
  std::vector<ES::M3d> d2AInv_dai_daj;

  M3x8d x = M3x8d::Zero();
  std::array<ES::M3d, kNumIntegrationPoints> Fref;
  std::array<ES::M3d, kNumIntegrationPoints> Fe;
  std::array<ES::M3d, kNumIntegrationPoints> U;
  std::array<ES::M3d, kNumIntegrationPoints> V;
  std::array<ES::V3d, kNumIntegrationPoints> S;
  std::array<M9x24d, kNumIntegrationPoints> dFdx;
  std::array<M3x8d, kNumIntegrationPoints> Bm;

  ES::VXd materialParam;
};

void CubicMeshDeformationModel::enableSPD(int)
{
}

CubicMeshDeformationModel::CubicMeshDeformationModel(const double restPositions[24], ElasticModel *elasticModel, PlasticModel *plasticModel):
  DeformationModel(elasticModel, plasticModel)
{
  ind = new CubicMeshDeformationModelInternal;

  for (int vi = 0; vi < kNumVertices; vi++) {
    ind->restX[vi] = ES::V3d(restPositions[vi * 3 + 0], restPositions[vi * 3 + 1], restPositions[vi * 3 + 2]);
  }

  const double offset = 0.5 / std::sqrt(3.0);
  const double gp[2] = { 0.5 - offset, 0.5 + offset };
  const double weight = 1.0 / 8.0;

  M3x8d X;
  for (int vi = 0; vi < kNumVertices; vi++) {
    X.col(vi) = ind->restX[vi];
  }

  int q = 0;
  for (int ia = 0; ia < 2; ia++) {
    for (int ib = 0; ib < 2; ib++) {
      for (int ig = 0; ig < 2; ig++, q++) {
        CubicMeshDeformationModelInternal::computeShapeDerivative(gp[ia], gp[ib], gp[ig], ind->dN_dabc[q]);

        ES::M3d Dm = X * ind->dN_dabc[q].transpose();
        double detDm = Dm.determinant();
        ind->restDmInv[q] = Dm.fullPivLu().inverse();
        ind->weightDetJ[q] = std::abs(detDm) * weight;

        M3x8d dN_dX = ind->restDmInv[q].transpose() * ind->dN_dabc[q];
        ind->restBm[q] = ind->weightDetJ[q] * dN_dX;
        CubicMeshDeformationModelInternal::compute_dF_dx(dN_dX, ind->rest_dFdx[q]);
      }
    }
  }

  ind->elasticModel = dynamic_cast<const ElasticModel3DDeformationGradient *>(elasticModel);
  ind->plasticModel = dynamic_cast<const PlasticModel3DDeformationGradient *>(plasticModel);
  if (ind->elasticModel == nullptr) {
    throw std::invalid_argument("CubicMeshDeformationModel requires ElasticModel3DDeformationGradient.");
  }
  if (ind->plasticModel == nullptr) {
    throw std::invalid_argument("CubicMeshDeformationModel requires PlasticModel3DDeformationGradient.");
  }

  numPlasticParams_ = ind->plasticModel->getNumParameters();
  numElasticParams_ = ind->elasticModel->getNumParameters();
}

CubicMeshDeformationModel::~CubicMeshDeformationModel()
{
  delete ind;
}

DeformationModelCacheData *CubicMeshDeformationModel::allocateCacheData() const
{
  return new CubicMeshDeformationModelCacheData(numPlasticParams_, numElasticParams_);
}

void CubicMeshDeformationModel::freeCacheData(DeformationModelCacheData *data) const
{
  delete data;
}

void CubicMeshDeformationModel::prepareData(const double *x, const double *param, const double *materialParam, CacheData *cacheDataBase) const
{
  CubicMeshDeformationModelCacheData *cacheData = static_cast<CubicMeshDeformationModelCacheData *>(cacheDataBase);

  for (int vi = 0; vi < kNumVertices; vi++) {
    cacheData->x.col(vi) = ES::V3d(x[vi * 3 + 0], x[vi * 3 + 1], x[vi * 3 + 2]);
  }

  ind->plasticModel->computeA(param, cacheData->Fp.data());
  ind->plasticModel->computeAInv(param, cacheData->FpInv.data());
  cacheData->detFp = ind->plasticModel->compute_detA(param);

  for (int i = 0; i < numPlasticParams_; i++) {
    cacheData->plasticParam[i] = param[i];
  }

  if (numPlasticParams_ > 0) {
    ind->plasticModel->compute_ddetA_da(cacheData->plasticParam.data(), cacheData->ddetA_da.data(), numPlasticParams_);
    ind->plasticModel->compute_d2detA_da2(cacheData->plasticParam.data(), cacheData->d2detA_da2.data(), numPlasticParams_);

    for (int i = 0; i < numPlasticParams_; i++) {
      ind->plasticModel->compute_dAInv_da(cacheData->plasticParam.data(), i, cacheData->dAInv_dai[i].data());

      for (int j = 0; j < numPlasticParams_; j++) {
        ind->plasticModel->compute_d2AInv_da2(cacheData->plasticParam.data(), i, j, cacheData->d2AInv(i, j).data());
      }
    }
  }

  for (int i = 0; i < numElasticParams_; i++) {
    cacheData->materialParam[i] = materialParam[i];
  }

  for (int q = 0; q < kNumIntegrationPoints; q++) {
    cacheData->Fref[q] = cacheData->x * ind->dN_dabc[q].transpose() * ind->restDmInv[q];
    cacheData->Fe[q] = cacheData->Fref[q] * cacheData->FpInv;
    CubicMeshDeformationModelInternal::computeSVD(cacheData->Fe[q], cacheData->U[q], cacheData->V[q], cacheData->S[q]);
    CubicMeshDeformationModelInternal::computeCurrent_dFdx(ind->rest_dFdx[q], cacheData->FpInv, cacheData->dFdx[q]);
    cacheData->Bm[q] = cacheData->detFp * cacheData->FpInv.transpose() * ind->restBm[q];
  }
}

double CubicMeshDeformationModel::computeEnergy(const CacheData *cacheDataBase) const
{
  const CubicMeshDeformationModelCacheData *cacheData = static_cast<const CubicMeshDeformationModelCacheData *>(cacheDataBase);
  const double *materialParam = optionalParamData(cacheData->materialParam);

  double energy = 0.0;
  for (int q = 0; q < kNumIntegrationPoints; q++) {
    energy += ind->elasticModel->compute_psi(materialParam, cacheData->Fe[q].data(),
      cacheData->U[q].data(), cacheData->V[q].data(), cacheData->S[q].data()) *
      ind->weightDetJ[q] * cacheData->detFp;
  }
  return energy;
}

void CubicMeshDeformationModel::compute_dE_dx(const CacheData *cacheDataBase, double *grad) const
{
  const CubicMeshDeformationModelCacheData *cacheData = static_cast<const CubicMeshDeformationModelCacheData *>(cacheDataBase);
  const double *materialParam = optionalParamData(cacheData->materialParam);

  Eigen::Map<ES::V24d> gradMap(grad);
  gradMap.setZero();

  for (int q = 0; q < kNumIntegrationPoints; q++) {
    ES::M3d P;
    ind->elasticModel->compute_P(materialParam, cacheData->Fe[q].data(),
      cacheData->U[q].data(), cacheData->V[q].data(), cacheData->S[q].data(), P.data());
    const Eigen::Matrix<double, 3, 8> localForce = P * cacheData->Bm[q];
    gradMap += Eigen::Map<const ES::V24d>(localForce.data());
  }
}

void CubicMeshDeformationModel::compute_d2E_dx2(const CacheData *cacheDataBase, double *hess) const
{
  const CubicMeshDeformationModelCacheData *cacheData = static_cast<const CubicMeshDeformationModelCacheData *>(cacheDataBase);
  const double *materialParam = optionalParamData(cacheData->materialParam);

  ES::M24d hessMat = ES::M24d::Zero();
  for (int q = 0; q < kNumIntegrationPoints; q++) {
    ES::M9d dPdF;
    ind->elasticModel->compute_dPdF(materialParam, cacheData->Fe[q].data(),
      cacheData->U[q].data(), cacheData->V[q].data(), cacheData->S[q].data(), dPdF.data());
    dPdF *= ind->weightDetJ[q] * cacheData->detFp;
    hessMat += cacheData->dFdx[q].transpose() * dPdF * cacheData->dFdx[q];
  }

  Eigen::Map<ES::M24d> hessMap(hess);
  hessMap = hessMat;
}

void CubicMeshDeformationModel::compute_dE_da(const CacheData *cacheDataBase, double *grad) const
{
  if (ind->plasticModel->getNumParameters() == 0)
    return;

  const CubicMeshDeformationModelCacheData *cacheData = static_cast<const CubicMeshDeformationModelCacheData *>(cacheDataBase);
  const double *materialParam = optionalParamData(cacheData->materialParam);

  ES::VXd gradVec = ES::VXd::Zero(numPlasticParams_);
  for (int q = 0; q < kNumIntegrationPoints; q++) {
    double psi = ind->elasticModel->compute_psi(materialParam, cacheData->Fe[q].data(),
      cacheData->U[q].data(), cacheData->V[q].data(), cacheData->S[q].data());

    ES::M3d P;
    ind->elasticModel->compute_P(materialParam, cacheData->Fe[q].data(),
      cacheData->U[q].data(), cacheData->V[q].data(), cacheData->S[q].data(), P.data());

    for (int i = 0; i < numPlasticParams_; i++) {
      double dVda = ind->compute_dV_dai(ind->weightDetJ[q], cacheData->ddetA_da[i]);
      double dpsi_da = ind->compute_dpsi_dai(cacheData->Fref[q], cacheData->dAInv_dai[i], P);
      gradVec[i] += dVda * psi + ind->weightDetJ[q] * cacheData->detFp * dpsi_da;
    }
  }

  for (int i = 0; i < numPlasticParams_; i++) {
    grad[i] = gradVec[i];
  }
}

void CubicMeshDeformationModel::compute_d2E_da2(const CacheData *cacheDataBase, double *hess) const
{
  if (ind->plasticModel->getNumParameters() == 0)
    return;

  const CubicMeshDeformationModelCacheData *cacheData = static_cast<const CubicMeshDeformationModelCacheData *>(cacheDataBase);
  const double *materialParam = optionalParamData(cacheData->materialParam);

  ES::MXd hessMat = ES::MXd::Zero(numPlasticParams_, numPlasticParams_);
  for (int q = 0; q < kNumIntegrationPoints; q++) {
    const double vol = ind->weightDetJ[q] * cacheData->detFp;
    double psi = ind->elasticModel->compute_psi(materialParam, cacheData->Fe[q].data(),
      cacheData->U[q].data(), cacheData->V[q].data(), cacheData->S[q].data());

    ES::M3d P;
    ind->elasticModel->compute_P(materialParam, cacheData->Fe[q].data(),
      cacheData->U[q].data(), cacheData->V[q].data(), cacheData->S[q].data(), P.data());

    ES::M9d dPdF;
    ind->elasticModel->compute_dPdF(materialParam, cacheData->Fe[q].data(),
      cacheData->U[q].data(), cacheData->V[q].data(), cacheData->S[q].data(), dPdF.data());

    for (int i = 0; i < numPlasticParams_; i++) {
      double dVda_i = ind->compute_dV_dai(ind->weightDetJ[q], cacheData->ddetA_da[i]);
      double dpsi_da_i = ind->compute_dpsi_dai(cacheData->Fref[q], cacheData->dAInv_dai[i], P);

      for (int j = 0; j < numPlasticParams_; j++) {
        double dVda_j = ind->compute_dV_dai(ind->weightDetJ[q], cacheData->ddetA_da[j]);
        double d2V = ind->compute_d2V_daidaj(ind->weightDetJ[q], cacheData->d2detA_da2(i, j));
        double dpsi_da_j = ind->compute_dpsi_dai(cacheData->Fref[q], cacheData->dAInv_dai[j], P);
        double d2psi = ind->compute_d2psi_dai_daj(cacheData->Fref[q], cacheData->dAInv_dai[i], cacheData->dAInv_dai[j],
          cacheData->d2AInv(i, j), P, dPdF);
        hessMat(i, j) += d2V * psi + dVda_i * dpsi_da_j + dVda_j * dpsi_da_i + vol * d2psi;
      }
    }
  }

  Eigen::Map<ES::MXd>(hess, numPlasticParams_, numPlasticParams_) = hessMat;
}

void CubicMeshDeformationModel::compute_d2E_dxda(const CacheData *cacheDataBase, double *hess) const
{
  if (ind->plasticModel->getNumParameters() == 0)
    return;

  const CubicMeshDeformationModelCacheData *cacheData = static_cast<const CubicMeshDeformationModelCacheData *>(cacheDataBase);
  const double *materialParam = optionalParamData(cacheData->materialParam);

  ES::MXd mixed = ES::MXd::Zero(24, numPlasticParams_);
  for (int q = 0; q < kNumIntegrationPoints; q++) {
    ES::M3d P;
    ind->elasticModel->compute_P(materialParam, cacheData->Fe[q].data(),
      cacheData->U[q].data(), cacheData->V[q].data(), cacheData->S[q].data(), P.data());

    ES::V24d dpsi_dx = cacheData->dFdx[q].transpose() * Eigen::Map<const ES::V9d>(P.data());

    ES::M9d dPdF;
    ind->elasticModel->compute_dPdF(materialParam, cacheData->Fe[q].data(),
      cacheData->U[q].data(), cacheData->V[q].data(), cacheData->S[q].data(), dPdF.data());

    const double vol = ind->weightDetJ[q] * cacheData->detFp;
    for (int i = 0; i < numPlasticParams_; i++) {
      const double dVda = ind->compute_dV_dai(ind->weightDetJ[q], cacheData->ddetA_da[i]);

      ES::M3d dFda;
      ind->compute_dFe_dai(cacheData->Fref[q], cacheData->dAInv_dai[i], dFda);

      ES::M3d dPda;
      CubicMeshDeformationModelInternal::compute_dP_dai(dPdF, dFda, dPda);

      ES::V24d temp = cacheData->dFdx[q].transpose() * Eigen::Map<const ES::V9d>(dPda.data());

      M9x24d d2F_duda;
      CubicMeshDeformationModelInternal::compute_d2Fe_dx_dai(cacheData->dAInv_dai[i], ind->rest_dFdx[q], d2F_duda);

      mixed.col(i) += dVda * dpsi_dx + vol * (temp + d2F_duda.transpose() * Eigen::Map<const ES::V9d>(P.data()));
    }
  }

  Eigen::Map<ES::MXd>(hess, 24, numPlasticParams_) = mixed;
}

void CubicMeshDeformationModel::compute_dE_db(const CacheData *cacheDataBase, double *grad) const
{
  if (ind->elasticModel->getNumParameters() == 0)
    return;

  const CubicMeshDeformationModelCacheData *cacheData = static_cast<const CubicMeshDeformationModelCacheData *>(cacheDataBase);
  const double *materialParam = optionalParamData(cacheData->materialParam);

  ES::VXd gradVec = ES::VXd::Zero(numElasticParams_);
  for (int q = 0; q < kNumIntegrationPoints; q++) {
    const double vol = ind->weightDetJ[q] * cacheData->detFp;
    for (int i = 0; i < numElasticParams_; i++) {
      gradVec[i] += vol * ind->elasticModel->compute_dpsi_dparam(materialParam, i,
        cacheData->Fe[q].data(), cacheData->U[q].data(), cacheData->V[q].data(), cacheData->S[q].data());
    }
  }

  for (int i = 0; i < numElasticParams_; i++) {
    grad[i] = gradVec[i];
  }
}

void CubicMeshDeformationModel::compute_d2E_db2(const CacheData *cacheDataBase, double *hess) const
{
  if (ind->elasticModel->getNumParameters() == 0)
    return;

  const CubicMeshDeformationModelCacheData *cacheData = static_cast<const CubicMeshDeformationModelCacheData *>(cacheDataBase);
  const double *materialParam = optionalParamData(cacheData->materialParam);

  ES::MXd hessMat = ES::MXd::Zero(numElasticParams_, numElasticParams_);
  for (int q = 0; q < kNumIntegrationPoints; q++) {
    const double vol = ind->weightDetJ[q] * cacheData->detFp;
    for (int i = 0; i < numElasticParams_; i++) {
      for (int j = 0; j < numElasticParams_; j++) {
        hessMat(i, j) += vol * ind->elasticModel->compute_d2psi_dparam2(materialParam, i, j,
          cacheData->Fe[q].data(), cacheData->U[q].data(), cacheData->V[q].data(), cacheData->S[q].data());
      }
    }
  }

  Eigen::Map<ES::MXd>(hess, numElasticParams_, numElasticParams_) = hessMat;
}

void CubicMeshDeformationModel::compute_d2E_dxdb(const CacheData *cacheDataBase, double *hess) const
{
  if (ind->elasticModel->getNumParameters() == 0)
    return;

  const CubicMeshDeformationModelCacheData *cacheData = static_cast<const CubicMeshDeformationModelCacheData *>(cacheDataBase);
  const double *materialParam = optionalParamData(cacheData->materialParam);

  ES::MXd mixed = ES::MXd::Zero(24, numElasticParams_);
  for (int q = 0; q < kNumIntegrationPoints; q++) {
    const double vol = ind->weightDetJ[q] * cacheData->detFp;
    for (int i = 0; i < numElasticParams_; i++) {
      ES::M3d dPdb;
      ind->elasticModel->compute_dP_dparam(materialParam, i, cacheData->Fe[q].data(),
        cacheData->U[q].data(), cacheData->V[q].data(), cacheData->S[q].data(), dPdb.data());
      mixed.col(i) += vol * (cacheData->dFdx[q].transpose() * Eigen::Map<const ES::V9d>(dPdb.data()));
    }
  }

  for (int col = 0; col < numElasticParams_; col++) {
    for (int row = 0; row < 24; row++) {
      hess[col * 24 + row] = mixed(row, col);
    }
  }
}

void CubicMeshDeformationModel::compute_d2E_dadb(const CacheData *cacheDataBase, double *hess) const
{
  if (ind->elasticModel->getNumParameters() == 0 || ind->plasticModel->getNumParameters() == 0)
    return;

  const CubicMeshDeformationModelCacheData *cacheData = static_cast<const CubicMeshDeformationModelCacheData *>(cacheDataBase);
  const double *materialParam = optionalParamData(cacheData->materialParam);

  Eigen::Map<ES::MXd> mixed(hess, numPlasticParams_, numElasticParams_);
  mixed.setZero();

  for (int q = 0; q < kNumIntegrationPoints; q++) {
    const double vol = ind->weightDetJ[q] * cacheData->detFp;

    std::vector<ES::M3d> dPdb(numElasticParams_, ES::M3d::Zero());
    for (int i = 0; i < numElasticParams_; i++) {
      ind->elasticModel->compute_dP_dparam(materialParam, i, cacheData->Fe[q].data(),
        cacheData->U[q].data(), cacheData->V[q].data(), cacheData->S[q].data(), dPdb[i].data());
    }

    std::vector<ES::M3d> dFda(numPlasticParams_, ES::M3d::Zero());
    for (int i = 0; i < numPlasticParams_; i++) {
      ind->compute_dFe_dai(cacheData->Fref[q], cacheData->dAInv_dai[i], dFda[i]);
      double dVda = ind->compute_dV_dai(ind->weightDetJ[q], cacheData->ddetA_da[i]);
      for (int j = 0; j < numElasticParams_; j++) {
        double dpsi_db = ind->elasticModel->compute_dpsi_dparam(materialParam, j,
          cacheData->Fe[q].data(), cacheData->U[q].data(), cacheData->V[q].data(), cacheData->S[q].data());
        mixed(i, j) += dVda * dpsi_db;
      }
    }

    for (int i = 0; i < numPlasticParams_; i++) {
      for (int j = 0; j < numElasticParams_; j++) {
        mixed(i, j) += vol * Eigen::Map<const ES::V9d>(dPdb[j].data()).dot(Eigen::Map<const ES::V9d>(dFda[i].data()));
      }
    }
  }
}

void CubicMeshDeformationModel::vonMisesStress(const CacheData *cacheDataBase, int &nPt, double *stresses) const
{
  const CubicMeshDeformationModelCacheData *cacheData = static_cast<const CubicMeshDeformationModelCacheData *>(cacheDataBase);
  const double *materialParam = optionalParamData(cacheData->materialParam);
  nPt = kNumIntegrationPoints;

  for (int q = 0; q < kNumIntegrationPoints; q++) {
    ES::M3d P;
    ind->elasticModel->compute_P(materialParam, cacheData->Fe[q].data(),
      cacheData->U[q].data(), cacheData->V[q].data(), cacheData->S[q].data(), P.data());

    const double detF = cacheData->Fe[q].determinant();
    ES::M3d cauchyStress = P * cacheData->Fe[q].transpose() / detF;

    const double t1 = std::pow(cauchyStress(0, 0) - cauchyStress(1, 1), 2.0);
    const double t2 = std::pow(cauchyStress(1, 1) - cauchyStress(2, 2), 2.0);
    const double t3 = std::pow(cauchyStress(2, 2) - cauchyStress(0, 0), 2.0);
    const double t4 = 6.0 * (std::pow(cauchyStress(1, 2), 2.0) + std::pow(cauchyStress(2, 0), 2.0) + std::pow(cauchyStress(0, 1), 2.0));
    stresses[q] = std::sqrt((t1 + t2 + t3 + t4) * 0.5);
  }
}

void CubicMeshDeformationModel::maxStrain(const CacheData *cacheDataBase, int &nPt, double *stresses) const
{
  const CubicMeshDeformationModelCacheData *cacheData = static_cast<const CubicMeshDeformationModelCacheData *>(cacheDataBase);
  nPt = kNumIntegrationPoints;

  for (int q = 0; q < kNumIntegrationPoints; q++) {
    ES::M3d E = 0.5 * (cacheData->Fe[q].transpose() * cacheData->Fe[q] - ES::M3d::Identity());
    Eigen::SelfAdjointEigenSolver<ES::M3d> eigSolver(E);
    stresses[q] = eigSolver.eigenvalues().maxCoeff();
  }
}

void CubicMeshDeformationModel::computeF(const double *x, int materialLocationID, double F[9]) const
{
  M3x8d xMat;
  for (int vi = 0; vi < kNumVertices; vi++) {
    xMat.col(vi) = ES::V3d(x[vi * 3 + 0], x[vi * 3 + 1], x[vi * 3 + 2]);
  }

  Eigen::Map<ES::M3d> FMap(F);
  FMap = xMat * ind->dN_dabc[materialLocationID].transpose() * ind->restDmInv[materialLocationID];
}

void CubicMeshDeformationModel::computeFe(const CacheData *cacheDataBase, int materialLocationID, double F[9]) const
{
  const CubicMeshDeformationModelCacheData *cacheData = static_cast<const CubicMeshDeformationModelCacheData *>(cacheDataBase);
  Eigen::Map<ES::M3d> FMap(F);
  FMap = cacheData->Fe[materialLocationID];
}

void CubicMeshDeformationModel::computeP(const CacheData *cacheDataBase, int materialLocationID, double POut[9]) const
{
  const CubicMeshDeformationModelCacheData *cacheData = static_cast<const CubicMeshDeformationModelCacheData *>(cacheDataBase);
  const double *materialParam = optionalParamData(cacheData->materialParam);

  ES::M3d P;
  ind->elasticModel->compute_P(materialParam, cacheData->Fe[materialLocationID].data(),
    cacheData->U[materialLocationID].data(), cacheData->V[materialLocationID].data(), cacheData->S[materialLocationID].data(), P.data());
  Eigen::Map<ES::M3d> PMap(POut);
  PMap = P;
}

void CubicMeshDeformationModel::computedPdF(const CacheData *cacheDataBase, int materialLocationID, double dPdFOut[81]) const
{
  const CubicMeshDeformationModelCacheData *cacheData = static_cast<const CubicMeshDeformationModelCacheData *>(cacheDataBase);
  const double *materialParam = optionalParamData(cacheData->materialParam);

  ES::M9d dPdF;
  ind->elasticModel->compute_dPdF(materialParam, cacheData->Fe[materialLocationID].data(),
    cacheData->U[materialLocationID].data(), cacheData->V[materialLocationID].data(), cacheData->S[materialLocationID].data(), dPdF.data());
  Eigen::Map<ES::M9d> dPdFMap(dPdFOut);
  dPdFMap = dPdF;
}

void CubicMeshDeformationModel::computedFdx(const CacheData *cacheDataBase, int materialLocationID, double *dFdxOut) const
{
  const CubicMeshDeformationModelCacheData *cacheData = static_cast<const CubicMeshDeformationModelCacheData *>(cacheDataBase);
  Eigen::Map<M9x24d> dFdxMap(dFdxOut);
  dFdxMap = cacheData->dFdx[materialLocationID];
}

void CubicMeshDeformationModel::computeForceFromP(const CacheData *cacheDataBase, int materialLocationID, const double P[9], double f[24]) const
{
  const CubicMeshDeformationModelCacheData *cacheData = static_cast<const CubicMeshDeformationModelCacheData *>(cacheDataBase);
  const Eigen::Map<const ES::M3d> PMap(P);
  const Eigen::Matrix<double, 3, 8> localForce = PMap * cacheData->Bm[materialLocationID];
  Eigen::Map<ES::V24d> fMap(f);
  fMap = Eigen::Map<const ES::V24d>(localForce.data());
}

double CubicMeshDeformationModel::getWeightDetJ(int materialLocationID) const
{
  return ind->weightDetJ[materialLocationID];
}

void CubicMeshDeformationModelInternal::computeShapeDerivative(double alpha, double beta, double gamma, M3x8d &dN)
{
  for (int i = 0; i < kNumVertices; i++) {
    const double alphaFactor = kVertexAlpha[i] ? alpha : (1.0 - alpha);
    const double betaFactor = kVertexBeta[i] ? beta : (1.0 - beta);
    const double gammaFactor = kVertexGamma[i] ? gamma : (1.0 - gamma);

    dN(0, i) = (kVertexAlpha[i] ? 1.0 : -1.0) * betaFactor * gammaFactor;
    dN(1, i) = (kVertexBeta[i] ? 1.0 : -1.0) * alphaFactor * gammaFactor;
    dN(2, i) = (kVertexGamma[i] ? 1.0 : -1.0) * alphaFactor * betaFactor;
  }
}

void CubicMeshDeformationModelInternal::computeSVD(const ES::M3d &Fe, ES::M3d &U, ES::M3d &V, ES::V3d &S)
{
  Eigen::JacobiSVD<ES::M3d, Eigen::NoQRPreconditioner> svd(Fe, Eigen::ComputeFullU | Eigen::ComputeFullV);
  U = svd.matrixU();
  V = svd.matrixV();
  S = svd.singularValues();

  if (U.determinant() < 0.0) {
    U.col(2) *= -1.0;
    S(2) *= -1.0;
  }
  if (V.determinant() < 0.0) {
    V.col(2) *= -1.0;
    S(2) *= -1.0;
  }
}

void CubicMeshDeformationModelInternal::compute_dF_dx(const M3x8d &dN_dX, M9x24d &dFdx)
{
  dFdx.setZero();
  for (int vi = 0; vi < kNumVertices; vi++) {
    for (int coord = 0; coord < 3; coord++) {
      const int dof = vi * 3 + coord;
      for (int deriv = 0; deriv < 3; deriv++) {
        dFdx(deriv * 3 + coord, dof) = dN_dX(deriv, vi);
      }
    }
  }
}

void CubicMeshDeformationModelInternal::computeCurrent_dFdx(const M9x24d &rest_dFdx, const ES::M3d &FpInv, M9x24d &dFdx)
{
  for (int col = 0; col < 24; col++) {
    const Eigen::Map<const ES::M3d> dFref(rest_dFdx.col(col).data());
    const ES::M3d dFe = dFref * FpInv;
    dFdx.col(col) = Eigen::Map<const ES::V9d>(dFe.data());
  }
}

void CubicMeshDeformationModelInternal::compute_d2Fe_dx_dai(const ES::M3d &dAInvdai, const M9x24d &rest_dFdx, M9x24d &d2Fdudai)
{
  for (int col = 0; col < 24; col++) {
    const Eigen::Map<const ES::M3d> dFref(rest_dFdx.col(col).data());
    const ES::M3d d2Fe = dFref * dAInvdai;
    d2Fdudai.col(col) = Eigen::Map<const ES::V9d>(d2Fe.data());
  }
}

void CubicMeshDeformationModelInternal::compute_dP_dai(const ES::M9d &dPdF, const ES::M3d &dFdai, ES::M3d &dPdai)
{
  Eigen::Map<ES::V9d>(dPdai.data()) = dPdF * Eigen::Map<const ES::V9d>(dFdai.data());
}

double CubicMeshDeformationModelInternal::compute_dV_dai(double weight, double ddetA_dai) const
{
  return weight * ddetA_dai;
}

double CubicMeshDeformationModelInternal::compute_d2V_daidaj(double weight, double d2detA_daidaj) const
{
  return weight * d2detA_daidaj;
}

void CubicMeshDeformationModelInternal::compute_dFe_dai(const ES::M3d &Fref, const ES::M3d &dAInvdai, ES::M3d &dFdai) const
{
  dFdai = Fref * dAInvdai;
}

void CubicMeshDeformationModelInternal::compute_d2Fe_dai_daj(const ES::M3d &Fref, const ES::M3d &dAInvdaidaj, ES::M3d &d2Fdaidaj) const
{
  d2Fdaidaj = Fref * dAInvdaidaj;
}

double CubicMeshDeformationModelInternal::compute_dpsi_dai(const ES::M3d &Fref, const ES::M3d &dAInv_dai, const ES::M3d &P) const
{
  ES::M3d dFe_dai;
  compute_dFe_dai(Fref, dAInv_dai, dFe_dai);
  return P.cwiseProduct(dFe_dai).sum();
}

double CubicMeshDeformationModelInternal::compute_d2psi_dai_daj(const ES::M3d &Fref, const ES::M3d &dAInv_dai,
  const ES::M3d &dAInv_daj, const ES::M3d &d2AInv_dai_daj, const ES::M3d &P, const ES::M9d &dPdF) const
{
  ES::M3d dFe_dai;
  ES::M3d dFe_daj;
  compute_dFe_dai(Fref, dAInv_dai, dFe_dai);
  compute_dFe_dai(Fref, dAInv_daj, dFe_daj);

  ES::M3d dP_daj;
  compute_dP_dai(dPdF, dFe_daj, dP_daj);

  ES::M3d d2Fe_daidaj;
  compute_d2Fe_dai_daj(Fref, d2AInv_dai_daj, d2Fe_daidaj);

  return dP_daj.cwiseProduct(dFe_dai).sum() + P.cwiseProduct(d2Fe_daidaj).sum();
}

}  // namespace SolidDeformationModel
}  // namespace pgo
