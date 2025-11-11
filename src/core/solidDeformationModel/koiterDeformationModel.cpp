#include "koiterDeformationModel.h"
#include "plasticModel2DFundamentalForms.h"
#include "elasticModel2DFundamentalForms.h"

#include "geometryQuery.h"
#include "EigenSupport.h"
#include "pgoLogging.h"

#include <mutex>

namespace pgo
{
namespace ES = pgo::EigenSupport;
namespace SolidDeformationModel
{
class KoiterDeformationModelInternal
{
public:
  ES::V3d restX[6];
  int hasVtx[6];
  int oppVtx[3] = { 4, 5, 3 };

  ES::M2d restI, restII;
  int enforceSPD;

  ElasticModel2DFundamentalForms *elasticModel;
  PlasticModel2DFundamentalForms *plasticModel;

  ES::M2d compute_a_and_derivatives(
    const ES::V3d x[3],
    Eigen::Matrix<double, 4, 9> *da_dx,
    ES::M9d ahess[4]);

  ES::M2d compute_b_and_derivatives(
    const ES::V3d x[6],
    int hasVtx[6],
    Eigen::Matrix<double, 4, 18> *db_dx,
    ES::M18d bhess[4]);

  ES::V3d secondFundamentalFormEntries(
    const ES::V3d x[6],
    int hasVtx[6],
    Eigen::Matrix<double, 3, 18> *derivative,
    ES::M18d hessian[3]);

  ES::M3d crossMatrix(Eigen::Vector3d v)
  {
    ES::M3d ret;
    ret << 0, -v[2], v[1],
      v[2], 0, -v[0],
      -v[1], v[0], 0;
    return ret;
  }

  ES::V3d faceNormal(const ES::V3d x0,
    const ES::V3d x1,
    const ES::V3d x2,
    Eigen::Matrix<double, 3, 9> *derivative,
    ES::M9d hessian[3]);
};

class KoiterDeformationModelCacheData : public DeformationModelCacheData
{
public:
  ES::V3d x[6];
  ES::M2d a, abar, b, bbar;
  ES::V18d elasticParams, plasticParams;

  double area;

  ElasticModel2DFundamentalForms *elasticModel;
  PlasticModel2DFundamentalForms *plasticModel;
};

}  // namespace SolidDeformationModel
}  // namespace pgo

using namespace pgo;
using namespace pgo::SolidDeformationModel;

KoiterDeformationModel::KoiterDeformationModel(const double X0[3], const double X1[3], const double X2[3],
  const double X3[3], const double X4[3], const double X5[3],
  ElasticModel *elasticModel, PlasticModel *plasticModel, int enforceSPD):
  SolidDeformationModel::DeformationModel(elasticModel, plasticModel)
{
  ind = new KoiterDeformationModelInternal;

  ind->enforceSPD = enforceSPD;

  ind->restX[0] = ES::V3d(X0[0], X0[1], X0[2]);
  ind->restX[1] = ES::V3d(X1[0], X1[1], X1[2]);
  ind->restX[2] = ES::V3d(X2[0], X2[1], X2[2]);
  ind->restX[3] = ES::V3d(X3[0], X3[1], X3[2]);
  ind->restX[4] = ES::V3d(X4[0], X4[1], X4[2]);
  ind->restX[5] = ES::V3d(X5[0], X5[1], X5[2]);

  ind->hasVtx[0] = 1;
  ind->hasVtx[1] = 1;
  ind->hasVtx[2] = 1;

  if (X3[0] == -10496)
    ind->hasVtx[3] = 0;
  else
    ind->hasVtx[3] = 1;

  if (X4[0] == -10496)
    ind->hasVtx[4] = 0;
  else
    ind->hasVtx[4] = 1;

  if (X5[0] == -10496)
    ind->hasVtx[5] = 0;
  else
    ind->hasVtx[5] = 1;

  ind->elasticModel = dynamic_cast<ElasticModel2DFundamentalForms *>(em);
  ind->plasticModel = dynamic_cast<PlasticModel2DFundamentalForms *>(pm);

  ES::M2d abar = ind->compute_a_and_derivatives(ind->restX, nullptr, nullptr);
  ES::M2d bbar = ind->compute_b_and_derivatives(ind->restX, ind->hasVtx, nullptr, nullptr);

  double area = (ind->restX[1] - ind->restX[0]).cross(ind->restX[2] - ind->restX[0]).norm() * 0.5;

  ind->plasticModel->set_abar(abar);
  ind->plasticModel->set_bbar(bbar);
  ind->plasticModel->setArea(area);
}

KoiterDeformationModel::~KoiterDeformationModel()
{
  delete ind;
}

DeformationModelCacheData *KoiterDeformationModel::allocateCacheData() const
{
  return new KoiterDeformationModelCacheData;
}

void KoiterDeformationModel::freeCacheData(CacheData *d) const
{
  delete d;
}

void KoiterDeformationModel::prepareData(const double *x, const double *param, const double *materialParam, CacheData *cacheDataBase) const
{
  KoiterDeformationModelCacheData *cacheData = dynamic_cast<KoiterDeformationModelCacheData *>(cacheDataBase);
  cacheData->x[0] = ES::V3d(x[0], x[1], x[2]);
  cacheData->x[1] = ES::V3d(x[3], x[4], x[5]);
  cacheData->x[2] = ES::V3d(x[6], x[7], x[8]);
  cacheData->x[3] = ES::V3d(x[9], x[10], x[11]);
  cacheData->x[4] = ES::V3d(x[12], x[13], x[14]);
  cacheData->x[5] = ES::V3d(x[15], x[16], x[17]);

  ind->plasticModel->compute_abar(param, cacheData->abar.data());
  ind->plasticModel->compute_bbar(param, cacheData->bbar.data());
  cacheData->area = ind->plasticModel->computeArea(param);

  cacheData->a = ind->compute_a_and_derivatives(cacheData->x, nullptr, nullptr);
  cacheData->b = ind->compute_b_and_derivatives(cacheData->x, ind->hasVtx, nullptr, nullptr);

  for (int i = 0; i < ind->elasticModel->getNumParameters(); i++) {
    cacheData->elasticParams[i] = materialParam[i];
  }

  for (int i = 0; i < ind->plasticModel->getNumParameters(); i++) {
    cacheData->plasticParams[i] = param[i];
  }
}

double KoiterDeformationModel::computeEnergy(const CacheData *cacheDataBase) const
{
  const KoiterDeformationModelCacheData *cacheData = dynamic_cast<const KoiterDeformationModelCacheData *>(cacheDataBase);

  double E1 = ind->elasticModel->compute_psi_a(cacheData->elasticParams.data(), cacheData->a.data(), cacheData->abar.data());
  double E2 = ind->elasticModel->compute_psi_b(cacheData->elasticParams.data(), cacheData->b.data(), cacheData->abar.data(), cacheData->bbar.data());

  return (E1 + E2) * cacheData->area;
}

void KoiterDeformationModel::compute_dE_dx(const CacheData *cacheDataBase, double *grad) const
{
  const KoiterDeformationModelCacheData *cacheData = dynamic_cast<const KoiterDeformationModelCacheData *>(cacheDataBase);

  ES::M4x9d dadx;
  ES::M4x18d dbdx;

  ind->compute_a_and_derivatives(cacheData->x, &dadx, nullptr);
  ind->compute_b_and_derivatives(cacheData->x, ind->hasVtx, &dbdx, nullptr);

  // dE1/dx = dE1/da da/dx
  ES::M2d dEda;
  ES::M2d dEdb;
  ind->elasticModel->compute_dpsi_da(cacheData->elasticParams.data(), cacheData->a.data(), cacheData->abar.data(), dEda.data());
  ind->elasticModel->compute_dpsi_db(cacheData->elasticParams.data(), cacheData->b.data(), cacheData->abar.data(), cacheData->bbar.data(), dEdb.data());

  (ES::Mp<ES::V18d>(grad)).setZero();
  (ES::Mp<ES::V9d>(grad)) = dadx.transpose() * ES::Mp<const ES::V4d>(dEda.data()) * cacheData->area;
  (ES::Mp<ES::V18d>(grad)) += dbdx.transpose() * ES::Mp<const ES::V4d>(dEdb.data()) * cacheData->area;
}

void KoiterDeformationModel::compute_d2E_dx2(const CacheData *cacheDataBase, double *hess) const
{
  const KoiterDeformationModelCacheData *cacheData = dynamic_cast<const KoiterDeformationModelCacheData *>(cacheDataBase);

  ES::M4x9d dadx;
  ES::M4x18d dbdx;

  ES::M9d d2adx2[4];
  ES::M18d d2bdx2[4];

  ind->compute_a_and_derivatives(cacheData->x, &dadx, d2adx2);
  ind->compute_b_and_derivatives(cacheData->x, ind->hasVtx, &dbdx, d2bdx2);

  // dE/dx = dE/da da/dx
  // d2E/dx2 = d(dE/da da/dx) = da/dx d2E/da2 da/dx + dE/da d2a/dx2

  ES::M2d dEda;
  ES::M2d dEdb;
  ind->elasticModel->compute_dpsi_da(cacheData->elasticParams.data(), cacheData->a.data(), cacheData->abar.data(), dEda.data());
  ind->elasticModel->compute_dpsi_db(cacheData->elasticParams.data(), cacheData->b.data(), cacheData->abar.data(), cacheData->bbar.data(), dEdb.data());

  ES::M4d d2Eda2;
  ES::M4d d2Edb2;
  ind->elasticModel->compute_d2psi_da2(cacheData->elasticParams.data(), cacheData->a.data(), cacheData->abar.data(), d2Eda2.data());
  ind->elasticModel->compute_d2psi_db2(cacheData->elasticParams.data(), cacheData->b.data(), cacheData->abar.data(), cacheData->bbar.data(), d2Edb2.data());

  ES::Mp<ES::M18d> hessMap(hess);
  hessMap.setZero();

  hessMap.block<9, 9>(0, 0) += dadx.transpose() * d2Eda2 * dadx * cacheData->area;
  for (int j = 0; j < 4; j++) {
    hessMap.block<9, 9>(0, 0) += dEda.data()[j] * d2adx2[j] * cacheData->area;
  }

  hessMap += dbdx.transpose() * d2Edb2 * dbdx * cacheData->area;
  for (int j = 0; j < 4; j++) {
    hessMap += dEdb.data()[j] * d2bdx2[j] * cacheData->area;
  }
}

void KoiterDeformationModel::compute_d2E_dxda(const CacheData *cacheDataBase, double *hess) const
{
  // E = psi * area
  // dE/dx = dpsi/da da/dx * area
  // d2E/dxdparam = d2psi/dadabar dabar/dparam da/dx * area + dpsi/da darea/dparam da/dx
  const KoiterDeformationModelCacheData *cacheData = dynamic_cast<const KoiterDeformationModelCacheData *>(cacheDataBase);
  ES::V4d dpsi_da, dpsi_db;
  ind->elasticModel->compute_dpsi_da(cacheData->elasticParams.data(), cacheData->a.data(), cacheData->abar.data(), dpsi_da.data());
  ind->elasticModel->compute_dpsi_db(cacheData->elasticParams.data(), cacheData->b.data(), cacheData->abar.data(), cacheData->bbar.data(), dpsi_db.data());

  ES::M4d d2psi_da_dabar, d2psi_db_dabar, d2psi_db_dbbar;
  ind->elasticModel->compute_d2psi_dadabar(cacheData->elasticParams.data(), cacheData->a.data(), cacheData->abar.data(), d2psi_da_dabar.data());
  ind->elasticModel->compute_d2psi_db_dabar(cacheData->elasticParams.data(), cacheData->b.data(), cacheData->abar.data(), cacheData->bbar.data(), d2psi_db_dabar.data());
  ind->elasticModel->compute_d2psi_db_dbbar(cacheData->elasticParams.data(), cacheData->b.data(), cacheData->abar.data(), cacheData->bbar.data(), d2psi_db_dbbar.data());

  ES::M4x9d dadx;
  ES::M4x18d dbdx;
  ind->compute_a_and_derivatives(cacheData->x, &dadx, nullptr);
  ind->compute_b_and_derivatives(cacheData->x, ind->hasVtx, &dbdx, nullptr);

  ES::M4x18d dabar_dF;
  ES::M4x18d dbbar_dF;
  ES::V18d darea_dF;
  ind->plasticModel->compute_dabar_dparam(cacheData->plasticParams.data(), dabar_dF.data());
  ind->plasticModel->compute_dbbar_dparam(cacheData->plasticParams.data(), dbbar_dF.data());
  ind->plasticModel->compute_darea_dparam(cacheData->plasticParams.data(), darea_dF.data());

  ES::V9d dpsi_a_dx = dpsi_da.transpose() * dadx;
  ES::V18d dpsi_b_dx = dpsi_db.transpose() * dbdx;

  ES::Mp<ES::MXd> hessMap(hess, 18, pm->getNumParameters());
  hessMap.setZero();

  hessMap.block(0, 0, 9, pm->getNumParameters()) += dadx.transpose() * d2psi_da_dabar * dabar_dF.leftCols(pm->getNumParameters()) * cacheData->area;
  hessMap.block(0, 0, 9, pm->getNumParameters()) += dpsi_a_dx * darea_dF.head(pm->getNumParameters()).transpose();

  hessMap.block(0, 0, 18, pm->getNumParameters()) += dbdx.transpose() * d2psi_db_dabar * dabar_dF.leftCols(pm->getNumParameters()) * cacheData->area;
  hessMap.block(0, 0, 18, pm->getNumParameters()) += dbdx.transpose() * d2psi_db_dbbar * dbbar_dF.leftCols(pm->getNumParameters()) * cacheData->area;
  hessMap.block(0, 0, 18, pm->getNumParameters()) += dpsi_b_dx * darea_dF.head(pm->getNumParameters()).transpose();
}

void KoiterDeformationModel::compute_d2E_dxdb(const CacheData *cacheDataBase, double *hess) const
{
  // E = psi * area
  // dE/dx = (dpsi/da da/dx + dpsi/db db/dx) * area
  // d2E/dxdparam = d2psi/dadparam da/dx * area + d2psi/dbdparam db/dx * area
  const KoiterDeformationModelCacheData *cacheData = dynamic_cast<const KoiterDeformationModelCacheData *>(cacheDataBase);

  ES::M4x18d d2psi_da_dparam, d2psi_db_dparam;
  ind->elasticModel->compute_d2psi_da_dparam(cacheData->elasticParams.data(), cacheData->a.data(), cacheData->abar.data(), d2psi_da_dparam.data());
  ind->elasticModel->compute_d2psi_db_dparam(cacheData->elasticParams.data(), cacheData->b.data(), cacheData->abar.data(), cacheData->bbar.data(), d2psi_db_dparam.data());

  ES::M4x9d dadx;
  ES::M4x18d dbdx;
  ind->compute_a_and_derivatives(cacheData->x, &dadx, nullptr);
  ind->compute_b_and_derivatives(cacheData->x, ind->hasVtx, &dbdx, nullptr);

  ES::Mp<ES::MXd> hessMap(hess, 18, em->getNumParameters());
  hessMap.setZero();

  hessMap.block(0, 0, 9, em->getNumParameters()) += dadx.transpose() * d2psi_da_dparam.leftCols(em->getNumParameters()) * cacheData->area;
  hessMap.block(0, 0, 18, em->getNumParameters()) += dbdx.transpose() * d2psi_db_dparam.leftCols(em->getNumParameters()) * cacheData->area;
}

ES::M2d KoiterDeformationModelInternal::compute_a_and_derivatives(const ES::V3d x[3],
  Eigen::Matrix<double, 4, 9> *da_dx, ES::M9d ahess[4])
{
  ES::V3d e1 = x[1] - x[0];
  ES::V3d e2 = x[2] - x[0];
  ES::M2d a;
  a(0, 0) = e1.squaredNorm();
  a(0, 1) = e1.dot(e2);
  a(1, 0) = a(0, 1);
  a(1, 1) = e2.squaredNorm();

  if (da_dx) {
    da_dx->setZero();
    da_dx->block<1, 3>(0, 0) = -2.0 * e1.transpose();
    da_dx->block<1, 3>(0, 3) = 2.0 * e1.transpose();

    Eigen::RowVector3d da01_dv0 = (-e1 - e2).transpose();
    da_dx->block<1, 3>(1, 0) = da01_dv0;
    da_dx->block<1, 3>(1, 3) = e2.transpose();
    da_dx->block<1, 3>(1, 6) = e1.transpose();
    da_dx->row(2) = da_dx->row(1);
    da_dx->block<1, 3>(3, 0) = -2.0 * e2.transpose();
    da_dx->block<1, 3>(3, 6) = 2.0 * e2.transpose();
  }

  if (ahess) {
    for (int i = 0; i < 4; ++i) {
      ahess[i].setZero();
    }

    ES::M3d I3 = ES::M3d::Identity();
    ahess[0].block<3, 3>(0, 0) = 2.0 * I3;
    ahess[0].block<3, 3>(0, 3) = -2.0 * I3;
    ahess[0].block<3, 3>(3, 0) = -2.0 * I3;
    ahess[0].block<3, 3>(3, 3) = 2.0 * I3;

    ahess[1].block<3, 3>(0, 0) = 2.0 * I3;
    ahess[1].block<3, 3>(0, 3) = -I3;
    ahess[1].block<3, 3>(0, 6) = -I3;
    ahess[1].block<3, 3>(3, 0) = -I3;
    ahess[1].block<3, 3>(3, 6) = I3;
    ahess[1].block<3, 3>(6, 0) = -I3;
    ahess[1].block<3, 3>(6, 3) = I3;

    ahess[2] = ahess[1];

    ahess[3].block<3, 3>(0, 0) = 2.0 * I3;
    ahess[3].block<3, 3>(0, 6) = -2.0 * I3;
    ahess[3].block<3, 3>(6, 0) = -2.0 * I3;
    ahess[3].block<3, 3>(6, 6) = 2.0 * I3;
  }
  return a;
}

ES::V3d KoiterDeformationModelInternal::secondFundamentalFormEntries(
  const ES::V3d x[6],
  int hasVtx[6],
  Eigen::Matrix<double, 3, 18> *derivative,
  ES::M18d hessian[3])
{
  if (derivative)
    derivative->setZero();

  if (hessian) {
    for (int i = 0; i < 3; i++)
      hessian[i].setZero();
  }
  ES::V3d II;

  ES::V3d oppNormals[3];
  Eigen::Matrix<double, 3, 9> dn[3];
  ES::M9d hn[3][3];

  Eigen::Matrix<double, 3, 9> dcn;
  ES::M9d hcn[3];
  ES::V3d cNormal = faceNormal(x[0], x[1], x[2],
    (derivative || hessian) ? &dcn : nullptr,
    hessian ? hcn : nullptr);

  for (int i = 0; i < 3; i++) {
    if (hasVtx[oppVtx[i]] == 0) {
      oppNormals[i].setZero();
      dn[i].setZero();
      for (int j = 0; j < 3; j++)
        hn[i][j].setZero();
    }
    else {
      ES::V3d x0 = x[oppVtx[i]];
      ES::V3d x1 = x[(i + 2) % 3];
      ES::V3d x2 = x[(i + 1) % 3];

      oppNormals[i] = faceNormal(x0, x1, x2,
        (derivative || hessian) ? &dn[i] : nullptr,
        hessian ? hn[i] : nullptr);
    }
  }

  ES::V3d qs[3];
  ES::V3d mvec[3];
  ES::V3d qvec[3];
  double mnorms[3];
  for (int i = 0; i < 3; i++) {
    qs[i] = x[i];
    mvec[i] = oppNormals[i] + cNormal;
    mnorms[i] = mvec[i].norm();
  }

  for (int i = 0; i < 3; i++) {
    int ip1 = (i + 1) % 3;
    int ip2 = (i + 2) % 3;
    qvec[i] = qs[ip1] + qs[ip2] - 2.0 * qs[i];
  }

  for (int i = 0; i < 3; i++) {
    int ip1 = (i + 1) % 3;
    int ip2 = (i + 2) % 3;
    II[i] = (qs[ip1] + qs[ip2] - 2.0 * qs[i]).dot(oppNormals[i]) / mnorms[i];

    if (derivative) {
      derivative->block<1, 3>(i, 3 * i) += -2.0 * oppNormals[i].transpose() / mnorms[i];
      derivative->block<1, 3>(i, 3 * ip1) += 1.0 * oppNormals[i].transpose() / mnorms[i];
      derivative->block<1, 3>(i, 3 * ip2) += 1.0 * oppNormals[i].transpose() / mnorms[i];

      derivative->block<1, 3>(i, 9 + 3 * (oppVtx[i] - 3)) += qvec[i].transpose() / mnorms[i] * dn[i].block<3, 3>(0, 0);
      derivative->block<1, 3>(i, 3 * ip2) += qvec[i].transpose() / mnorms[i] * dn[i].block<3, 3>(0, 3);
      derivative->block<1, 3>(i, 3 * ip1) += qvec[i].transpose() / mnorms[i] * dn[i].block<3, 3>(0, 6);

      derivative->block<1, 3>(i, 9 + 3 * (oppVtx[i] - 3)) += -qvec[i].dot(oppNormals[i]) / mnorms[i] / mnorms[i] / mnorms[i] * mvec[i].transpose() * dn[i].block<3, 3>(0, 0);
      derivative->block<1, 3>(i, 3 * ip2) += -qvec[i].dot(oppNormals[i]) / mnorms[i] / mnorms[i] / mnorms[i] * mvec[i].transpose() * dn[i].block<3, 3>(0, 3);
      derivative->block<1, 3>(i, 3 * ip1) += -qvec[i].dot(oppNormals[i]) / mnorms[i] / mnorms[i] / mnorms[i] * mvec[i].transpose() * dn[i].block<3, 3>(0, 6);

      derivative->block<1, 3>(i, 0) += -qvec[i].dot(oppNormals[i]) / mnorms[i] / mnorms[i] / mnorms[i] * mvec[i].transpose() * dcn.block<3, 3>(0, 0);
      derivative->block<1, 3>(i, 3) += -qvec[i].dot(oppNormals[i]) / mnorms[i] / mnorms[i] / mnorms[i] * mvec[i].transpose() * dcn.block<3, 3>(0, 3);
      derivative->block<1, 3>(i, 6) += -qvec[i].dot(oppNormals[i]) / mnorms[i] / mnorms[i] / mnorms[i] * mvec[i].transpose() * dcn.block<3, 3>(0, 6);
    }

    if (hessian) {
      int ip1 = (i + 1) % 3;
      int ip2 = (i + 2) % 3;

      int miidx[3];
      miidx[0] = 9 + 3 * (oppVtx[i] - 3);
      miidx[1] = 3 * ip2;
      miidx[2] = 3 * ip1;

      ES::M3d dnij[3];
      for (int j = 0; j < 3; j++)
        dnij[j] = dn[i].block<3, 3>(0, 3 * j);

      for (int j = 0; j < 3; j++) {
        hessian[i].block<3, 3>(miidx[j], 3 * ip1) += (1.0 / mnorms[i]) * dnij[j].transpose();
        hessian[i].block<3, 3>(miidx[j], 3 * ip2) += (1.0 / mnorms[i]) * dnij[j].transpose();
        hessian[i].block<3, 3>(miidx[j], 3 * i) += (-2.0 / mnorms[i]) * dnij[j].transpose();

        ES::V3d dnijTm = dnij[j].transpose() * mvec[i];
        ES::V3d dcnjTm = (dcn.block<3, 3>(0, 3 * j).transpose() * mvec[i]);
        ES::M3d term3 = dnijTm * oppNormals[i].transpose();

        hessian[i].block<3, 3>(miidx[j], 3 * ip1) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * term3;
        hessian[i].block<3, 3>(miidx[j], 3 * ip2) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * term3;
        hessian[i].block<3, 3>(miidx[j], 3 * i) += (2.0 / mnorms[i] / mnorms[i] / mnorms[i]) * term3;

        ES::M3d term4 = dcnjTm * oppNormals[i].transpose();

        hessian[i].block<3, 3>(3 * j, 3 * ip1) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * term4;
        hessian[i].block<3, 3>(3 * j, 3 * ip2) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * term4;
        hessian[i].block<3, 3>(3 * j, 3 * i) += (2.0 / mnorms[i] / mnorms[i] / mnorms[i]) * term4;

        hessian[i].block<3, 3>(3 * ip1, miidx[j]) += (1.0 / mnorms[i]) * dnij[j];
        hessian[i].block<3, 3>(3 * ip2, miidx[j]) += (1.0 / mnorms[i]) * dnij[j];
        hessian[i].block<3, 3>(3 * i, miidx[j]) += (-2.0 / mnorms[i]) * dnij[j];

        for (int k = 0; k < 3; k++) {
          hessian[i].block<3, 3>(miidx[j], miidx[k]) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * (dnij[j].transpose() * mvec[i]) * (qvec[i].transpose() * dnij[k]);
          hessian[i].block<3, 3>(3 * j, miidx[k]) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * (dcn.block<3, 3>(0, 3 * j).transpose() * mvec[i]) * (qvec[i].transpose() * dnij[k]);
        }

        ES::M3d term1 = oppNormals[i] * dnijTm.transpose();
        hessian[i].block<3, 3>(3 * ip1, miidx[j]) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * term1;
        hessian[i].block<3, 3>(3 * ip2, miidx[j]) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * term1;
        hessian[i].block<3, 3>(3 * i, miidx[j]) += (2.0 / mnorms[i] / mnorms[i] / mnorms[i]) * term1;

        ES::M3d term2 = oppNormals[i] * dcnjTm.transpose();
        hessian[i].block<3, 3>(3 * ip1, 3 * j) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * term2;
        hessian[i].block<3, 3>(3 * ip2, 3 * j) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * term2;
        hessian[i].block<3, 3>(3 * i, 3 * j) += (2.0 / mnorms[i] / mnorms[i] / mnorms[i]) * term2;

        ES::V3d dnijTq = dnij[j].transpose() * qvec[i];

        for (int k = 0; k < 3; k++) {
          hessian[i].block<3, 3>(miidx[j], miidx[k]) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * dnijTq * (mvec[i].transpose() * dnij[k]);
          hessian[i].block<3, 3>(miidx[j], 3 * k) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * dnijTq * (mvec[i].transpose() * dcn.block<3, 3>(0, 3 * k));
        }

        double qdoto = qvec[i].dot(oppNormals[i]);

        for (int k = 0; k < 3; k++) {
          hessian[i].block<3, 3>(miidx[j], miidx[k]) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * qdoto * dnij[j].transpose() * dnij[k];
          hessian[i].block<3, 3>(miidx[j], 3 * k) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * qdoto * dnij[j].transpose() * dcn.block<3, 3>(0, 3 * k);
          hessian[i].block<3, 3>(3 * j, miidx[k]) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * qdoto * dcn.block<3, 3>(0, 3 * j).transpose() * dnij[k];
          hessian[i].block<3, 3>(3 * j, 3 * k) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * qdoto * dcn.block<3, 3>(0, 3 * j).transpose() * dcn.block<3, 3>(0, 3 * k);

          hessian[i].block<3, 3>(miidx[j], miidx[k]) += (3.0 / mnorms[i] / mnorms[i] / mnorms[i] / mnorms[i] / mnorms[i]) * qdoto * dnijTm * (mvec[i].transpose() * dnij[k]);
          hessian[i].block<3, 3>(miidx[j], 3 * k) += (3.0 / mnorms[i] / mnorms[i] / mnorms[i] / mnorms[i] / mnorms[i]) * qdoto * dnijTm * (mvec[i].transpose() * dcn.block<3, 3>(0, 3 * k));
          hessian[i].block<3, 3>(3 * j, miidx[k]) += (3.0 / mnorms[i] / mnorms[i] / mnorms[i] / mnorms[i] / mnorms[i]) * qdoto * dcnjTm * (mvec[i].transpose() * dnij[k]);
          hessian[i].block<3, 3>(3 * j, 3 * k) += (3.0 / mnorms[i] / mnorms[i] / mnorms[i] / mnorms[i] / mnorms[i]) * qdoto * dcnjTm * (mvec[i].transpose() * dcn.block<3, 3>(0, 3 * k));

          for (int l = 0; l < 3; l++) {
            hessian[i].block<3, 3>(miidx[j], miidx[k]) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * qdoto * mvec[i][l] * hn[i][l].block<3, 3>(3 * j, 3 * k);
            hessian[i].block<3, 3>(3 * j, 3 * k) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * qdoto * mvec[i][l] * hcn[l].block<3, 3>(3 * j, 3 * k);
            hessian[i].block<3, 3>(miidx[j], miidx[k]) += (1.0 / mnorms[i]) * qvec[i][l] * hn[i][l].block<3, 3>(3 * j, 3 * k);
          }
        }
      }
    }
  }

  return II;
}

ES::M2d KoiterDeformationModelInternal::compute_b_and_derivatives(
  const ES::V3d x[6],
  int hasVtx[6],
  Eigen::Matrix<double, 4, 18> *db_dx,
  ES::M18d bhess[4])
{
  if (db_dx) {
    db_dx->setZero();
  }
  if (bhess) {
    for (int i = 0; i < 4; i++) {
      bhess[i].setZero();
    }
  }
  Eigen::Matrix<double, 3, 18> IIderiv;
  ES::M18d IIhess[3];

  ES::V3d II = secondFundamentalFormEntries(x, hasVtx,
    db_dx ? &IIderiv : nullptr,
    bhess ? IIhess : nullptr);

  ES::M2d result;
  result << II[0] + II[1], II[0], II[0], II[0] + II[2];

  if (db_dx) {
    db_dx->row(0) += IIderiv.row(0);
    db_dx->row(0) += IIderiv.row(1);

    db_dx->row(1) += IIderiv.row(0);
    db_dx->row(2) += IIderiv.row(0);

    db_dx->row(3) += IIderiv.row(0);
    db_dx->row(3) += IIderiv.row(2);
  }

  if (bhess) {
    bhess[0] += IIhess[0];
    bhess[0] += IIhess[1];

    bhess[1] += IIhess[0];
    bhess[2] += IIhess[0];

    bhess[3] += IIhess[0];
    bhess[3] += IIhess[2];
  }

  return result;
}

ES::V3d KoiterDeformationModelInternal::faceNormal(
  const ES::V3d x0, const ES::V3d x1, const ES::V3d x2,
  Eigen::Matrix<double, 3, 9> *derivative,
  ES::M9d hessian[3])
{
  if (derivative)
    derivative->setZero();

  if (hessian) {
    for (int i = 0; i < 3; i++)
      hessian[i].setZero();
  }

  ES::V3d qi0(x0);
  ES::V3d qi1(x1);
  ES::V3d qi2(x2);

  ES::V3d n = (qi1 - qi0).cross(qi2 - qi0);

  if (derivative) {
    derivative->block(0, 0, 3, 3) += crossMatrix(qi2 - qi1);
    derivative->block(0, 3, 3, 3) += crossMatrix(qi0 - qi2);
    derivative->block(0, 6, 3, 3) += crossMatrix(qi1 - qi0);
  }

  if (hessian) {
    for (int j = 0; j < 3; j++) {
      Eigen::Vector3d ej(0, 0, 0);
      ej[j] = 1.0;
      ES::M3d ejc = crossMatrix(ej);
      hessian[j].block(0, 3, 3, 3) -= ejc;
      hessian[j].block(0, 6, 3, 3) += ejc;
      hessian[j].block(3, 6, 3, 3) -= ejc;
      hessian[j].block(3, 0, 3, 3) += ejc;
      hessian[j].block(6, 0, 3, 3) -= ejc;
      hessian[j].block(6, 3, 3, 3) += ejc;
    }
  }

  return n;
}