/*
author: Bohan Wang
copyright to USC,MIT
*/

#include "tetMeshDeformationModel.h"
#include "elasticModel.h"
#include "plasticModel.h"

// #define MKL_DIRECT_CALL_SEQ_JIT

#include "geometryQuery.h"
#include "EigenSupport.h"
#include "pgoLogging.h"

#include <mutex>

namespace ES = pgo::EigenSupport;

typedef Eigen::Matrix<double, 3, 4> M3x4d;
typedef Eigen::Matrix<double, 3, 12> M3x12d;

namespace pgo

{
namespace SolidDeformationModel
{
class TetMeshDeformationModelInternal
{
public:
  ES::V3d restX[4];
  ES::M3d restDmInv;
  M3x4d restBm;
  ES::M9x12d rest_dFdx;
  double restVolume;

  // volume related derivatives
  double compute_dV_dai(double ddetA_dai) const;
  double compute_d2V_daidaj(double d2detA_daidaj) const;
  // F related derivatives
  void compute_dFe_dai(const ES::M3d &Ds, const ES::M3d &dAInvdai, ES::M3d &dFdai) const;
  void compute_d2Fe_dai_daj(const ES::M3d &Ds, const ES::M3d &dAInvdaidaj, ES::M3d &d2Fdaidaj) const;
  void compute_d2Fe_dx_dai(const ES::M3d &dAInvdai, const ES::M9x12d &dFdx, ES::M9x12d &d2Fdudai) const;
  // P related derivatives
  void compute_dP_dai(const ES::M9d &dPdF, const ES::M3d &dFdai, ES::M3d &dPdai) const;
  // psi related derivatives
  double compute_dpsi_dai(const ES::M3d &Ds, const ES::M3d &dAInv_dai, const ES::M3d &P) const;
  double compute_d2psi_dai_daj(const ES::M3d &Ds, const ES::M3d &dAInv_dai, const ES::M3d &dAInv_daj,
    const ES::M3d &d2AInv_dai_daj, const ES::M3d &P, const ES::M9d &dPdF) const;
  // basic
  void computeDs(const ES::V3d &x0, const ES::V3d &x1, const ES::V3d &x2, const ES::V3d &x3, ES::M3d &Ds) const;
  void computeDs(const ES::V12d &x, ES::M3d &Ds) const;
  void computeDmInv(const ES::V3d &X0, const ES::V3d &X1, const ES::V3d &X2, const ES::V3d &X3, ES::M3d &DmInv) const;
  double computeVolume(const double X0[3], const double X1[3], const double X2[3], const double X3[3]) const;
  void compute_dF_dx(const ES::M3d &DmInv, ES::M9x12d &dFdx) const;
  void compute_Bm(const ES::M3d &DmInv, double volume, M3x4d &Bm) const;

  void computeF(const ES::M3d &Ds, const ES::M3d &FpInv, ES::M3d &Fe) const;
  void computeSVD(const ES::M3d &Fe, ES::M3d &U, ES::M3d &V, ES::V3d &S) const;
};

class TetMeshDeformationModelCacheData : public DeformationModelCacheData
{
public:
  const static int maxNumPlasticParam = 12;
  const static int maxNumElasticParam = 12;

  typedef Eigen::Matrix<double, maxNumPlasticParam, 1> Vp;
  typedef Eigen::Matrix<double, maxNumPlasticParam, maxNumPlasticParam> Mp;

  Vp plasticParam;
  ES::M3d Fp, FpInv;
  double detFp;

  Vp ddetA_da;
  Mp d2detA_da2;

  ES::M3d dAInv_dai[maxNumPlasticParam];
  ES::M3d d2AInv_dai_daj[maxNumPlasticParam][maxNumPlasticParam];

  ES::V3d x[4];
  ES::M3d Ds;
  ES::M3d Fe;
  ES::M3d U, V;
  ES::V3d S;

  Eigen::Matrix<double, 100, 1> materialParam;
};

}  // namespace SolidDeformationModel
}  // namespace pgo

using namespace pgo::SolidDeformationModel;

TetMeshDeformationModel::TetMeshDeformationModel(
  const double X0[3], const double X1[3], const double X2[3], const double X3[3],
  const ElasticModel *elasticModel, const PlasticModel *plasticModel):
  DeformationModel(elasticModel, plasticModel)
{
  ind = new TetMeshDeformationModelInternal;
  ind->restX[0] = ES::V3d(X0[0], X0[1], X0[2]);
  ind->restX[1] = ES::V3d(X1[0], X1[1], X1[2]);
  ind->restX[2] = ES::V3d(X2[0], X2[1], X2[2]);
  ind->restX[3] = ES::V3d(X3[0], X3[1], X3[2]);

  ind->restVolume = ind->computeVolume(X0, X1, X2, X3);
  ind->computeDmInv(ind->restX[0], ind->restX[1], ind->restX[2], ind->restX[3], ind->restDmInv);
  ind->compute_Bm(ind->restDmInv, ind->restVolume, ind->restBm);
  ind->compute_dF_dx(ind->restDmInv, ind->rest_dFdx);
}

TetMeshDeformationModel::~TetMeshDeformationModel()
{
  delete ind;
}

DeformationModelCacheData *TetMeshDeformationModel::allocateCacheData() const
{
  return new TetMeshDeformationModelCacheData;
}

void TetMeshDeformationModel::freeCacheData(DeformationModelCacheData *d) const
{
  delete d;
}

void TetMeshDeformationModel::computeF(const double *x, int materialLocationIDs, double F[9]) const
{
  ES::V3d xs[4] = {
    ES::V3d(x[0], x[1], x[2]),
    ES::V3d(x[3], x[4], x[5]),
    ES::V3d(x[6], x[7], x[8]),
    ES::V3d(x[9], x[10], x[11])
  };

  ES::M3d Ds;
  ind->computeDs(xs[0], xs[1], xs[2], xs[3], Ds);
  (Eigen::Map<ES::M3d>(F)) = Ds * ind->restDmInv;
}

void TetMeshDeformationModel::computeP(const CacheData *cacheDataBase, int materialLocationIDs, double POut[9]) const
{
  const TetMeshDeformationModelCacheData *cacheData = static_cast<const TetMeshDeformationModelCacheData *>(cacheDataBase);

  ES::M3d P;
  em->compute_P(cacheData->materialParam.data(), cacheData->Fe.data(),
    cacheData->U.data(), cacheData->V.data(), cacheData->S.data(), P.data());

  (ES::Mp<ES::M3d>(POut)) = P;
}

void TetMeshDeformationModel::computedPdF(const CacheData *cacheDataBase, int materialLocationIDs, double dPdFOut[81]) const
{
  const TetMeshDeformationModelCacheData *cacheData = static_cast<const TetMeshDeformationModelCacheData *>(cacheDataBase);

  ES::M9d dPdF;
  em->compute_dPdF(cacheData->materialParam.data(), cacheData->Fe.data(),
    cacheData->U.data(), cacheData->V.data(), cacheData->S.data(), dPdF.data());

  (Eigen::Map<ES::M9d>(dPdFOut)) = dPdF;
}

void TetMeshDeformationModel::computedFdx(const CacheData *cacheDataBase, int materialLocationIDs, double *dFdxOut) const
{
  const TetMeshDeformationModelCacheData *cacheData = static_cast<const TetMeshDeformationModelCacheData *>(cacheDataBase);

  ES::M9x12d dFdx;
  ES::M3d DmInv = ind->restDmInv * cacheData->FpInv;
  ind->compute_dF_dx(DmInv, dFdx);

  (ES::Mp<ES::M9x12d>(dFdxOut)) = dFdx;
}

void TetMeshDeformationModel::prepareData(const double *x, const double *param, const double *materialParam, CacheData *cacheDataBase) const
{
  TetMeshDeformationModelCacheData *cacheData = static_cast<TetMeshDeformationModelCacheData *>(cacheDataBase);

  cacheData->x[0] = ES::V3d(x[0], x[1], x[2]);
  cacheData->x[1] = ES::V3d(x[3], x[4], x[5]);
  cacheData->x[2] = ES::V3d(x[6], x[7], x[8]);
  cacheData->x[3] = ES::V3d(x[9], x[10], x[11]);

  ind->computeDs(cacheData->x[0], cacheData->x[1], cacheData->x[2], cacheData->x[3], cacheData->Ds);

  pm->computeA(param, cacheData->Fp.data());
  // if (fabs(cacheData->Fp.determinant() - 1) > 1e-8) {
  //   std::cout << cacheData->Fp << std::endl;
  // }

  pm->computeAInv(param, cacheData->FpInv.data());
  cacheData->detFp = pm->compute_detA(param);

  for (int i = 0; i < pm->getNumParameters(); i++)
    cacheData->plasticParam[i] = param[i];

  pm->compute_ddetA_da(cacheData->plasticParam.data(), cacheData->ddetA_da.data(), TetMeshDeformationModelCacheData::maxNumPlasticParam);
  pm->compute_d2detA_da2(cacheData->plasticParam.data(), cacheData->d2detA_da2.data(), TetMeshDeformationModelCacheData::maxNumPlasticParam);

  for (int i = 0; i < pm->getNumParameters(); i++) {
    pm->compute_dAInv_da(cacheData->plasticParam.data(), i, cacheData->dAInv_dai[i].data());

    for (int j = 0; j < pm->getNumParameters(); j++) {
      pm->compute_d2AInv_da2(cacheData->plasticParam.data(), i, j, cacheData->d2AInv_dai_daj[i][j].data());
    }
  }

  ind->computeF(cacheData->Ds, cacheData->FpInv, cacheData->Fe);
  ind->computeSVD(cacheData->Fe, cacheData->U, cacheData->V, cacheData->S);

  for (int i = 0; em && i < em->getNumParameters(); i++)
    cacheData->materialParam[i] = materialParam[i];
}

double TetMeshDeformationModel::computeEnergy(const CacheData *cacheDataBase) const
{
  const TetMeshDeformationModelCacheData *cacheData = static_cast<const TetMeshDeformationModelCacheData *>(cacheDataBase);

  return em->compute_psi(cacheData->materialParam.data(), cacheData->Fe.data(),
           cacheData->U.data(), cacheData->V.data(), cacheData->S.data()) *
    ind->restVolume * cacheData->detFp;
}

void TetMeshDeformationModel::vonMisesStress(const CacheData *cacheDataBase, int &nPt, double *stresses) const
{
  const TetMeshDeformationModelCacheData *cacheData = static_cast<const TetMeshDeformationModelCacheData *>(cacheDataBase);

  ES::M3d P;
  em->compute_P(cacheData->materialParam.data(), cacheData->Fe.data(),
    cacheData->U.data(), cacheData->V.data(), cacheData->S.data(), P.data());

  double detF = cacheData->Fe.determinant();
  ES::M3d cauchyStress = P * cacheData->Fe.transpose() / detF;
  // P = J sigma F^-T
  // J^-1 P FT = sigma

  double t1 = (cauchyStress(0, 0) - cauchyStress(1, 1)) * (cauchyStress(0, 0) - cauchyStress(1, 1));
  double t2 = (cauchyStress(1, 1) - cauchyStress(2, 2)) * (cauchyStress(1, 1) - cauchyStress(2, 2));
  double t3 = (cauchyStress(2, 2) - cauchyStress(0, 0)) * (cauchyStress(2, 2) - cauchyStress(0, 0));
  double t4 = 6 * (cauchyStress(1, 2) * cauchyStress(1, 2) + cauchyStress(2, 0) * cauchyStress(2, 0) + cauchyStress(0, 1) * cauchyStress(0, 1));
  double v = std::sqrt((t1 + t2 + t3 + t4) * 0.5);

  nPt = 1;
  stresses[0] = v;
}

void TetMeshDeformationModel::maxStrain(const CacheData *cacheDataBase, int &nPt, double *stresses) const
{
  const TetMeshDeformationModelCacheData *cacheData = static_cast<const TetMeshDeformationModelCacheData *>(cacheDataBase);

  ES::M3d E = 0.5 * (cacheData->Fe.transpose() * cacheData->Fe - ES::M3d::Identity());
  Eigen::SelfAdjointEigenSolver<ES::M3d> eigSolver(E);
  nPt = 1;
  stresses[0] = eigSolver.eigenvalues().maxCoeff();
}

void TetMeshDeformationModel::computeForceFromP(const CacheData *cacheDataBase, const double P[9], double f[12]) const
{
  const TetMeshDeformationModelCacheData *cacheData = static_cast<const TetMeshDeformationModelCacheData *>(cacheDataBase);

  M3x4d Bm = cacheData->detFp * cacheData->FpInv.transpose() * ind->restBm;
  (Eigen::Map<M3x4d>(f)) = ES::Mp<const ES::M3d>(P) * Bm;
}

void TetMeshDeformationModel::compute_dE_dx(const CacheData *cacheDataBase, double *grad) const
{
  const TetMeshDeformationModelCacheData *cacheData = static_cast<const TetMeshDeformationModelCacheData *>(cacheDataBase);

  ES::M3d P;
  em->compute_P(cacheData->materialParam.data(), cacheData->Fe.data(),
    cacheData->U.data(), cacheData->V.data(), cacheData->S.data(), P.data());

  M3x4d Bm = cacheData->detFp * cacheData->FpInv.transpose() * ind->restBm;
  (Eigen::Map<M3x4d>(grad)) = P * Bm;
}

void TetMeshDeformationModel::compute_d2E_dx2(const CacheData *cacheDataBase, double *hess) const
{
  const TetMeshDeformationModelCacheData *cacheData = static_cast<const TetMeshDeformationModelCacheData *>(cacheDataBase);

  ES::M9d dPdF;
  em->compute_dPdF(cacheData->materialParam.data(), cacheData->Fe.data(),
    cacheData->U.data(), cacheData->V.data(), cacheData->S.data(), dPdF.data());

  dPdF *= ind->restVolume * cacheData->detFp;

  // we don't compute this through dF/du * Fp^-1,
  // because the follow should involve less operations
  ES::M9x12d dFdx;
  ES::M3d DmInv = ind->restDmInv * cacheData->FpInv;
  ind->compute_dF_dx(DmInv, dFdx);

  (Eigen::Map<ES::M12d>(hess)) = dFdx.transpose() * dPdF * dFdx;
}

void TetMeshDeformationModel::compute_d3E_dx3(const CacheData *cacheDataBase, double *tensor) const
{
  const TetMeshDeformationModelCacheData *cacheData = static_cast<const TetMeshDeformationModelCacheData *>(cacheDataBase);

  ES::V12d x0;
  x0.segment<3>(0) = cacheData->x[0];
  x0.segment<3>(3) = cacheData->x[1];
  x0.segment<3>(6) = cacheData->x[2];
  x0.segment<3>(9) = cacheData->x[3];

#if defined(USE_FINITE_DIFF)
  double eps = 1e-7;
  double step[2] = { -eps, eps };

  for (int dof = 0; dof < 12; dof++) {
    ES::V12d x1 = x0;
    ES::M12d hess[2];
    ES::M3d Ds;
    ES::M3d Fe;
    ES::M3d U, V;
    ES::V3d S;

    for (int i = 0; i < 2; i++) {
      x1[dof] = x0[dof] + step[i];

      ind->computeDs(x1, Ds);
      ind->computeF(Ds, cacheData->FpInv, Fe);
      ind->computeSVD(Fe, U, V, S);

      ES::M9d dPdF;
      em->compute_dPdF(cacheData->materialParam.data(),
        Fe.data(), U.data(), V.data(), S.data(), dPdF.data());
      dPdF *= ind->restVolume * cacheData->detFp;

      ES::M9x12d dFdx;
      ES::M3d DmInv = ind->restDmInv * cacheData->FpInv;
      ind->compute_dF_dx(DmInv, dFdx);

      hess[i] = dFdx.transpose() * dPdF * dFdx;
    }

    (Eigen::Map<ES::M12d>(tensor + dof * 144)) = (hess[1] - hess[0]) * 0.5 / eps;
  }
#else
  double d3psidF3[729];
  em->compute_d2PdF2(cacheData->materialParam.data(), cacheData->Fe.data(),
    cacheData->U.data(), cacheData->V.data(), cacheData->S.data(), d3psidF3);

  double volume = ind->restVolume * cacheData->detFp;

  // we don't compute this through dF/du * Fp^-1,
  // because the follow should involve less operations
  ES::M9x12d dFdx;
  ES::M3d DmInv = ind->restDmInv * cacheData->FpInv;
  ind->compute_dF_dx(DmInv, dFdx);

  auto getEntry = [](double tensor[729], int i, int j, int k) -> double & {
    return tensor[k * 81 + j * 9 + i];
  };

  // dE/dxi = V * dpsi/dFj * dFj/dxi
  // d2E/(dxi dxj) = V * sum_r (d(dpsi/dFr dFr/dxi) / dxj)
  //               = V * sum_r (d2psi/(dFrdxj) * dFr/dxi)
  //               = V * sum_r (sum_s (d2psi/(dFr dFs) dFs/dxj) * dFr/dxi)

  // d3E/(dxi dxj dxl) = V * d(sum_r (sum_s (d2psi/(dFr dFs) dFs/dxj) * dFr/dxi)) dxk
  //                   = V * sum_r (sum_s (d3psi/(dFr dFs dxk) dFs/dxj) * dFr/dxi)
  //                   = V * sum_r (sum_s (sum_t (d3psi/(dFr dFs dFt) * dFt/dxk) dFs/dxj) dFr/dxi)

  for (int i = 0; i < 12; i++) {
    for (int j = 0; j < 12; j++) {
      for (int k = 0; k < 12; k++) {
        double entry = 0.0;

        // for each entry
        for (int r = 0; r < 9; r++) {
          double sum_s = 0;

          for (int s = 0; s < 9; s++) {
            double sum_t = 0;

            for (int t = 0; t < 9; t++) {
              sum_t += getEntry(d3psidF3, r, s, t) * dFdx(t, k);
            }  // end for t

            sum_s += sum_t * dFdx(s, j);
          }  // end for s

          entry += sum_s * dFdx(r, j);
        }  // end for r

        d3E_dx3_ijk(tensor, i, j, k) = entry * volume;
      }
    }
  }

#endif
}

void TetMeshDeformationModel::compute_d3E_dxdadx(const CacheData *cacheDataBase, double *tensor) const
{
}

void TetMeshDeformationModel::compute_d3E_dxdada(const CacheData *cacheDataBase, double *tensor) const
{
}

void TetMeshDeformationModel::compute_dE_da(const CacheData *cacheDataBase, double *grad) const
{
  if (pm->getNumParameters() == 0)
    return;

  const TetMeshDeformationModelCacheData *cacheData = static_cast<const TetMeshDeformationModelCacheData *>(cacheDataBase);

  TetMeshDeformationModelCacheData::Vp dvolda;
  for (int i = 0; i < pm->getNumParameters(); i++) {
    dvolda[i] = ind->compute_dV_dai(cacheData->ddetA_da[i]);
  }

  double psi = em->compute_psi(cacheData->materialParam.data(), cacheData->Fe.data(),
    cacheData->U.data(), cacheData->V.data(), cacheData->S.data());
  double vol = ind->restVolume * cacheData->detFp;

  ES::M3d P;
  em->compute_P(cacheData->materialParam.data(), cacheData->Fe.data(),
    cacheData->U.data(), cacheData->V.data(), cacheData->S.data(), P.data());

  TetMeshDeformationModelCacheData::Vp dpsi_da;
  for (int i = 0; i < pm->getNumParameters(); i++) {
    dpsi_da[i] = ind->compute_dpsi_dai(cacheData->Ds, cacheData->dAInv_dai[i], P);
  }

  for (int i = 0; i < pm->getNumParameters(); i++) {
    grad[i] = dvolda[i] * psi + vol * dpsi_da[i];
  }
}

void TetMeshDeformationModel::compute_d2E_da2(const CacheData *cacheDataBase, double *hess) const
{
  if (pm->getNumParameters() == 0)
    return;

  const TetMeshDeformationModelCacheData *cacheData = static_cast<const TetMeshDeformationModelCacheData *>(cacheDataBase);

  // V
  double vol = ind->restVolume * cacheData->detFp;

  // dV da
  TetMeshDeformationModelCacheData::Vp dvolda;
  for (int i = 0; i < pm->getNumParameters(); i++) {
    dvolda[i] = ind->compute_dV_dai(cacheData->ddetA_da[i]);
  }

  // d2V / da2
  TetMeshDeformationModelCacheData::Mp d2volda2;
  for (int i = 0; i < pm->getNumParameters(); i++) {
    for (int j = i; j < pm->getNumParameters(); j++) {
      d2volda2(i, j) = ind->compute_d2V_daidaj(cacheData->d2detA_da2(i, j));
    }
  }

  for (int i = 0; i < pm->getNumParameters(); i++) {
    for (int j = 0; j < i; j++) {
      d2volda2(i, j) = d2volda2(j, i);
    }
  }

  // psi
  double psi = em->compute_psi(cacheData->materialParam.data(),
    cacheData->Fe.data(), cacheData->U.data(), cacheData->V.data(), cacheData->S.data());

  // dpsi / da
  ES::M3d P;
  em->compute_P(cacheData->materialParam.data(), cacheData->Fe.data(),
    cacheData->U.data(), cacheData->V.data(), cacheData->S.data(), P.data());

  TetMeshDeformationModelCacheData::Vp dpsi_da;
  for (int i = 0; i < pm->getNumParameters(); i++) {
    dpsi_da[i] = ind->compute_dpsi_dai(cacheData->Ds, cacheData->dAInv_dai[i], P);
  }

  // d2psi / da2
  ES::M9d dPdF;
  em->compute_dPdF(cacheData->materialParam.data(), cacheData->Fe.data(),
    cacheData->U.data(), cacheData->V.data(), cacheData->S.data(), dPdF.data());

  TetMeshDeformationModelCacheData::Mp d2psi_da2;
  for (int i = 0; i < pm->getNumParameters(); i++) {
    for (int j = i; j < pm->getNumParameters(); j++) {
      d2psi_da2(i, j) = ind->compute_d2psi_dai_daj(cacheData->Ds, cacheData->dAInv_dai[i], cacheData->dAInv_dai[j],
        cacheData->d2AInv_dai_daj[i][j], P, dPdF);
    }
  }

  for (int i = 0; i < pm->getNumParameters(); i++) {
    for (int j = 0; j < i; j++) {
      d2psi_da2(i, j) = d2psi_da2(j, i);
    }
  }

  Eigen::Map<ES::MXd> hessMap(hess, pm->getNumParameters(), pm->getNumParameters());
  for (int i = 0; i < pm->getNumParameters(); i++) {
    for (int j = 0; j < pm->getNumParameters(); j++) {
      hessMap(i, j) = d2volda2(i, j) * psi + dvolda[i] * dpsi_da[j] + dvolda[j] * dpsi_da[i] + vol * d2psi_da2(i, j);
    }
  }
}

void TetMeshDeformationModel::compute_d2E_dxda(const CacheData *cacheDataBase, double *hess) const
{
  if (pm->getNumParameters() == 0)
    return;

  const TetMeshDeformationModelCacheData *cacheData = static_cast<const TetMeshDeformationModelCacheData *>(cacheDataBase);

  Eigen::Matrix<double, 12, TetMeshDeformationModelCacheData::maxNumPlasticParam> d2E_duda;

  ES::M9x12d dFdx;
  ES::M3d DmInv = ind->restDmInv * cacheData->FpInv;
  ind->compute_dF_dx(DmInv, dFdx);

  // P = dpsi/dF
  ES::M3d P;
  em->compute_P(cacheData->materialParam.data(), cacheData->Fe.data(),
    cacheData->U.data(), cacheData->V.data(), cacheData->S.data(), P.data());

  // dpsi/dx = P * dF/dx
  ES::V12d dpsi_dx = dFdx.transpose() * Eigen::Map<const ES::V9d>(P.data());

  // dP/dF
  ES::M9d dPdF;
  em->compute_dPdF(cacheData->materialParam.data(), cacheData->Fe.data(),
    cacheData->U.data(), cacheData->V.data(), cacheData->S.data(), dPdF.data());

  // V
  double vol = ind->restVolume * cacheData->detFp;

  for (int i = 0; i < pm->getNumParameters(); i++) {
    double dVda = ind->compute_dV_dai(cacheData->ddetA_da[i]);

    // dP/da
    ES::M3d dFda, dPda;
    ind->compute_dFe_dai(cacheData->Ds, cacheData->dAInv_dai[i], dFda);
    ind->compute_dP_dai(dPdF, dFda, dPda);

    ES::V12d temp = dFdx.transpose() * Eigen::Map<const ES::V9d>(dPda.data());

    ES::M9x12d d2F_duda;
    ind->compute_d2Fe_dx_dai(cacheData->dAInv_dai[i], ind->rest_dFdx, d2F_duda);

    d2E_duda.col(i) = dVda * dpsi_dx + vol * (temp + d2F_duda.transpose() * Eigen::Map<const ES::V9d>(P.data()));
  }

  (Eigen::Map<ES::MXd>(hess, 12, pm->getNumParameters())) = d2E_duda.block(0, 0, 12, pm->getNumParameters());
}

void TetMeshDeformationModel::compute_dE_db(const CacheData *cacheDataBase, double *grad) const
{
  if (em->getNumParameters() == 0)
    return;

  const TetMeshDeformationModelCacheData *cacheData = static_cast<const TetMeshDeformationModelCacheData *>(cacheDataBase);

  // V
  double V = ind->restVolume * cacheData->detFp;
  for (int i = 0; i < em->getNumParameters(); i++) {
    double val = em->compute_dpsi_dparam(cacheData->materialParam.data(), i,
      cacheData->Fe.data(), cacheData->U.data(), cacheData->V.data(), cacheData->S.data());
    grad[i] = val * V;
  }
}

void TetMeshDeformationModel::compute_d2E_db2(const CacheData *cacheDataBase, double *hess) const
{
  if (em->getNumParameters() == 0)
    return;

  const TetMeshDeformationModelCacheData *cacheData = static_cast<const TetMeshDeformationModelCacheData *>(cacheDataBase);

  // V
  double V = ind->restVolume * cacheData->detFp;
  for (int i = 0; i < em->getNumParameters(); i++) {
    for (int j = 0; j < em->getNumParameters(); j++) {
      double val = em->compute_d2psi_dparam2(cacheData->materialParam.data(), i, j,
        cacheData->Fe.data(), cacheData->U.data(), cacheData->V.data(), cacheData->S.data());
      hess[j * em->getNumParameters() + i] = val * V;
    }
  }
}

void TetMeshDeformationModel::compute_d2E_dxdb(const CacheData *cacheDataBase, double *hess) const
{
  if (em->getNumParameters() == 0)
    return;

  const TetMeshDeformationModelCacheData *cacheData = static_cast<const TetMeshDeformationModelCacheData *>(cacheDataBase);

  ES::M9x12d dFdx;
  ES::M3d DmInv = ind->restDmInv * cacheData->FpInv;
  ind->compute_dF_dx(DmInv, dFdx);

  // V
  double V = ind->restVolume * cacheData->detFp;

  Eigen::Matrix<double, 12, 100> d2E_dudb;
  for (int i = 0; i < em->getNumParameters(); i++) {
    ES::M3d dPdb;
    em->compute_dP_dparam(cacheData->materialParam.data(), i, cacheData->Fe.data(),
      cacheData->U.data(), cacheData->V.data(), cacheData->S.data(), dPdb.data());

    ES::V12d temp = dFdx.transpose() * Eigen::Map<const ES::V9d>(dPdb.data());
    d2E_dudb.col(i) = V * temp;
  }

  for (int col = 0; col < em->getNumParameters(); col++) {
    for (int row = 0; row < 12; row++) {
      hess[col * 12 + row] = d2E_dudb(row, col);
    }
  }
}

void TetMeshDeformationModel::compute_d2E_dadb(const CacheData *cacheDataBase, double *hess) const
{
  if (em->getNumParameters() == 0)
    return;

  if (pm->getNumParameters() == 0)
    return;

  const TetMeshDeformationModelCacheData *cacheData = static_cast<const TetMeshDeformationModelCacheData *>(cacheDataBase);

  Eigen::Map<ES::MXd> d2E_dadb(hess, pm->getNumParameters(), em->getNumParameters());

  // V
  double V = ind->restVolume * cacheData->detFp;

  TetMeshDeformationModelCacheData::Vp dVda;
  for (int i = 0; i < pm->getNumParameters(); i++) {
    dVda[i] = ind->compute_dV_dai(cacheData->ddetA_da[i]);
  }

  for (int i = 0; i < em->getNumParameters(); i++) {
    double dpsi_dbi = em->compute_dpsi_dparam(cacheData->materialParam.data(), i,
      cacheData->Fe.data(), cacheData->U.data(), cacheData->V.data(), cacheData->S.data());

    d2E_dadb.col(i) = dVda.head(pm->getNumParameters()) * dpsi_dbi;
  }

  ES::M3d dPdb[TetMeshDeformationModelCacheData::maxNumElasticParam];
  for (int i = 0; i < em->getNumParameters(); i++) {
    em->compute_dP_dparam(cacheData->materialParam.data(), i, cacheData->Fe.data(),
      cacheData->U.data(), cacheData->V.data(), cacheData->S.data(), dPdb[i].data());
  }

  ES::M3d dFda[TetMeshDeformationModelCacheData::maxNumPlasticParam];
  for (int i = 0; i < pm->getNumParameters(); i++) {
    ind->compute_dFe_dai(cacheData->Ds, cacheData->dAInv_dai[i], dFda[i]);
  }

  for (int i = 0; i < pm->getNumParameters(); i++) {
    for (int j = 0; j < em->getNumParameters(); j++) {
      d2E_dadb(i, j) += Eigen::Map<const ES::V9d>(dPdb[j].data()).dot(Eigen::Map<const ES::V9d>(dFda[i].data())) * V;
    }
  }
}

#if 0
void MuscleElasticMaterialModelWithParameterForTetMesh::compute_d3E_dudadu(MuscleElasticMaterialModelWithParameterForTetMeshCacheData *cacheData, double *tensor) const
{
  ES::V12d x0;
  x0.segment<3>(0) = cacheData->x[0];
  x0.segment<3>(3) = cacheData->x[1];
  x0.segment<3>(6) = cacheData->x[2];
  x0.segment<3>(9) = cacheData->x[3];

  double eps = 1e-7;
  double step[2] = { -eps, eps };

  // we don't compute this through dF/du * RT S^-1 R,
  // because the follow should involve less operations
  ES::M9x12d dFduPrime;
  ES::M3d DmInvPrime = data->restDmInv * data->R.transpose() * cacheData->AInv * data->R;
  compute_dF_du(DmInvPrime.data(), dFduPrime.data());
  // V
  double vol = data->restVolume * cacheData->detA;

  for (int dof = 0; dof < 12; dof++) {
    ES::V12d x1 = x0;
    Eigen::Matrix<double, 12, 3> d2E_duda[2];
    ES::M3d Ds;
    ES::M3d F;
    ES::M3d U, V;
    ES::V3d S;

    for (int si = 0; si < 2; si++) {
      x1[dof] = x0[dof] + step[si];

      computeDs(x1.data(), x1.data() + 3, x1.data() + 6, x1.data() + 9, Ds.data());
      computeF(Ds.data(), cacheData->AInv.data(), F.data());
      computeSVD(F.data(), U.data(), V.data(), S.data());

      // P = dpsi/dF
      ES::M3d P;
      femModel->compute_P(cacheData->materialParam.data(),
        F.data(), U.data(), V.data(), S.data(), P.data());

      // dpsi/du = P * dF/du
      ES::V12d dpsi_du = dFduPrime.transpose() * Eigen::Map<Eigen::Matrix<double, 9, 1>>(P.data());

      // dP/dF
      ES::M9d dPdF;
      femModel->compute_dPdF(cacheData->materialParam.data(),
        F.data(), U.data(), V.data(), S.data(), dPdF.data());

      for (int i = 0; i < 3; i++) {
        double dVda = compute_dV_da(cacheData->ddetA_da[i], i);

        // dP/da
        ES::M3d dPda;
        ES::M3d dFda;
        compute_dF_da(Ds.data(), cacheData->dAInv_da[i].data(), i, dFda.data());
        compute_dP_da(dPdF.data(), dFda.data(), dPda.data());
        Eigen::Matrix<double, 12, 1> temp = dFduPrime.transpose() * Eigen::Map<Eigen::Matrix<double, 9, 1>>(dPda.data());

        ES::M9x12d d2F_duda;
        compute_d2F_duda(cacheData->dAInv_da[i].data(), data->dFdu.data(), i, d2F_duda.data());

        d2E_duda[si].col(i) = dVda * dpsi_du + vol * (temp + d2F_duda.transpose() * Eigen::Map<Eigen::Matrix<double, 9, 1>>(P.data()));
      }
    }

    (Eigen::Map<Eigen::Matrix<double, 12, 3>>(tensor + 36 * dof)) = (d2E_duda[1] - d2E_duda[0]) * 0.5 / eps;
  }
}

void MuscleElasticMaterialModelWithParameterForTetMesh::compute_d3E_dudada(MuscleElasticMaterialModelWithParameterForTetMeshCacheData *cacheData, double *tensor) const
{
  ES::V3d alpha0 = cacheData->alpha;
  double eps = 1e-7;
  double step[2] = { -eps, eps };

  // P = dpsi/dF
  ES::M3d P;
  femModel->compute_P(cacheData->materialParam.data(),
    cacheData->F.data(), cacheData->U.data(), cacheData->V.data(), cacheData->S.data(), P.data());

  // dP/dF
  ES::M9d dPdF;
  femModel->compute_dPdF(cacheData->materialParam.data(),
    cacheData->F.data(), cacheData->U.data(), cacheData->V.data(), cacheData->S.data(), dPdF.data());

  for (int dof = 0; dof < 3; dof++) {
    ES::V3d alpha = alpha0;
    ES::M3d A, AInv;
    ES::M9x12d dFduPrime;
    Eigen::Matrix<double, 12, 3> d2E_duda[2];

    for (int si = 0; si < 2; si++) {
      alpha[dof] = alpha0[dof] + step[si];
      muscleModel->computeA(alpha.data(), A.data());
      AInv = A.fullPivHouseholderQr().inverse();

      ES::M3d DmInvPrime = data->restDmInv * data->R.transpose() * AInv * data->R;
      compute_dF_du(DmInvPrime.data(), dFduPrime.data());

      double detA = A.determinant();
      double vol = data->restVolume * detA;

      // dpsi/du = P * dF/du
      ES::V12d dpsi_du = dFduPrime.transpose() * Eigen::Map<Eigen::Matrix<double, 9, 1>>(P.data());

      ES::V3d ddetA_da;
      muscleModel->compute_ddetA_da(alpha.data(), ddetA_da.data());

      ES::M3d dAInv_da[3];
      muscleModel->compute_dAInv_da(alpha.data(), 0, dAInv_da[0].data());
      muscleModel->compute_dAInv_da(alpha.data(), 1, dAInv_da[1].data());
      muscleModel->compute_dAInv_da(alpha.data(), 2, dAInv_da[2].data());

      for (int i = 0; i < 3; i++) {
        double dVda = compute_dV_da(ddetA_da[i], i);

        // dP/da
        ES::M3d dPda;
        ES::M3d dFda;
        compute_dF_da(cacheData->Ds.data(), dAInv_da[i].data(), i, dFda.data());
        compute_dP_da(dPdF.data(), dFda.data(), dPda.data());
        Eigen::Matrix<double, 12, 1> temp = dFduPrime.transpose() * Eigen::Map<Eigen::Matrix<double, 9, 1>>(dPda.data());

        ES::M9x12d d2F_duda;
        compute_d2F_duda(dAInv_da[i].data(), data->dFdu.data(), i, d2F_duda.data());

        d2E_duda[si].col(i) = dVda * dpsi_du + vol * (temp + d2F_duda.transpose() * Eigen::Map<Eigen::Matrix<double, 9, 1>>(P.data()));
      }
    }

    (Eigen::Map<Eigen::Matrix<double, 12, 3>>(tensor + 36 * dof)) = (d2E_duda[1] - d2E_duda[0]) * 0.5 / eps;
  }
}

#endif

double TetMeshDeformationModelInternal::compute_dV_dai(double ddetA_dai) const
{
  return restVolume * ddetA_dai;
}

double TetMeshDeformationModelInternal::compute_d2V_daidaj(double d2detA_daidaj) const
{
  return restVolume * d2detA_daidaj;
}

void TetMeshDeformationModelInternal::compute_dFe_dai(const ES::M3d &Ds, const ES::M3d &dAInvdai, ES::M3d &dFdai) const
{
  // F * d(Fp^-1)_dai
  dFdai = Ds * restDmInv * dAInvdai;
}

void TetMeshDeformationModelInternal::compute_d2Fe_dai_daj(const ES::M3d &Ds, const ES::M3d &dAInvdaidaj, ES::M3d &d2Fdaidaj) const
{
  // F * d2(Fp^-1)_dai_daj
  d2Fdaidaj = Ds * restDmInv * dAInvdaidaj;
}

void TetMeshDeformationModelInternal::compute_d2Fe_dx_dai(const ES::M3d &dAInvdai, const ES::M9x12d &dFdx, ES::M9x12d &d2Fdudai) const
{
  d2Fdudai = ES::M9x12d::Zero();
  for (int k = 0; k < 12; k++) {
    // dFdx *
    ES::M3d d2Fedukdai = Eigen::Map<const ES::M3d>(dFdx.data() + k * 9) * dAInvdai;
    d2Fdudai.col(k) = Eigen::Map<const ES::V9d>(d2Fedukdai.data());
  }
}

void TetMeshDeformationModelInternal::compute_dP_dai(const ES::M9d &dPdF, const ES::M3d &dFdai, ES::M3d &dPdai) const
{
  (Eigen::Map<ES::V9d>(dPdai.data())) = dPdF * Eigen::Map<const ES::V9d>(dFdai.data());
}

double TetMeshDeformationModelInternal::compute_dpsi_dai(const ES::M3d &Ds, const ES::M3d &dAInv_dai, const ES::M3d &P) const
{
  ES::M3d dFe_dai;
  compute_dFe_dai(Ds, dAInv_dai, dFe_dai);
  return P.cwiseProduct(dFe_dai).sum();
}

double TetMeshDeformationModelInternal::compute_d2psi_dai_daj(const ES::M3d &Ds,
  const ES::M3d &dAInv_dai, const ES::M3d &dAInv_daj, const ES::M3d &d2AInv_dai_daj, const ES::M3d &P, const ES::M9d &dPdF) const
{
  ES::M3d dFe_dai, dFe_daj;
  compute_dFe_dai(Ds, dAInv_dai, dFe_dai);
  compute_dFe_dai(Ds, dAInv_daj, dFe_daj);

  ES::M3d dP_daj;
  compute_dP_dai(dPdF, dFe_daj, dP_daj);

  ES::M3d d2F_dai_daj;
  compute_d2Fe_dai_daj(Ds, d2AInv_dai_daj, d2F_dai_daj);

  return dP_daj.cwiseProduct(dFe_dai).sum() + P.cwiseProduct(d2F_dai_daj).sum();
}

void TetMeshDeformationModelInternal::computeDs(const ES::V3d &x0, const ES::V3d &x1,
  const ES::V3d &x2, const ES::V3d &x3, ES::M3d &Ds) const
{
  Ds.col(0) = x1 - x0;
  Ds.col(1) = x2 - x0;
  Ds.col(2) = x3 - x0;
}

void TetMeshDeformationModelInternal::computeDs(const ES::V12d &x, ES::M3d &Ds) const
{
  Ds.col(0) = x.segment<3>(3) - x.segment<3>(0);
  Ds.col(1) = x.segment<3>(6) - x.segment<3>(0);
  Ds.col(2) = x.segment<3>(9) - x.segment<3>(0);
}

void TetMeshDeformationModelInternal::computeDmInv(const ES::V3d &X0, const ES::V3d &X1,
  const ES::V3d &X2, const ES::V3d &X3, ES::M3d &DmInv) const
{
  ES::M3d Dm;
  Dm.col(0) = X1 - X0;
  Dm.col(1) = X2 - X0;
  Dm.col(2) = X3 - X0;

  DmInv = Dm.fullPivLu().inverse();
}

double TetMeshDeformationModelInternal::computeVolume(const double X0[3], const double X1[3],
  const double X2[3], const double X3[3]) const
{
  return Mesh::getTetVolume(asVec3d(X0), asVec3d(X1), asVec3d(X2), asVec3d(X3));
}

void TetMeshDeformationModelInternal::compute_dF_dx(const ES::M3d &DmInv, ES::M9x12d &dFdx) const
{
  dFdx = ES::M9x12d::Zero();

  double v0 = -(DmInv(0, 0) + DmInv(1, 0) + DmInv(2, 0));
  double v1 = -(DmInv(0, 1) + DmInv(1, 1) + DmInv(2, 1));
  double v2 = -(DmInv(0, 2) + DmInv(1, 2) + DmInv(2, 2));

  // F is in column major
  // dF / dx0
  dFdx(0, 0) = v0;
  dFdx(3, 0) = v1;
  dFdx(6, 0) = v2;

  // dF / dx1
  dFdx(1, 1) = v0;
  dFdx(4, 1) = v1;
  dFdx(7, 1) = v2;

  // dF / dx2
  dFdx(2, 2) = v0;
  dFdx(5, 2) = v1;
  dFdx(8, 2) = v2;

  // dF / dx3
  dFdx(0, 3) = DmInv(0, 0);
  dFdx(3, 3) = DmInv(0, 1);
  dFdx(6, 3) = DmInv(0, 2);

  // dF / dx4
  dFdx(1, 4) = DmInv(0, 0);
  dFdx(4, 4) = DmInv(0, 1);
  dFdx(7, 4) = DmInv(0, 2);

  // dF / dx5
  dFdx(2, 5) = DmInv(0, 0);
  dFdx(5, 5) = DmInv(0, 1);
  dFdx(8, 5) = DmInv(0, 2);

  // dF / dx6
  dFdx(0, 6) = DmInv(1, 0);
  dFdx(3, 6) = DmInv(1, 1);
  dFdx(6, 6) = DmInv(1, 2);

  // dF / du7
  dFdx(1, 7) = DmInv(1, 0);
  dFdx(4, 7) = DmInv(1, 1);
  dFdx(7, 7) = DmInv(1, 2);

  // dF / du8
  dFdx(2, 8) = DmInv(1, 0);
  dFdx(5, 8) = DmInv(1, 1);
  dFdx(8, 8) = DmInv(1, 2);

  // dF / du9
  dFdx(0, 9) = DmInv(2, 0);
  dFdx(3, 9) = DmInv(2, 1);
  dFdx(6, 9) = DmInv(2, 2);

  // dF / du10
  dFdx(1, 10) = DmInv(2, 0);
  dFdx(4, 10) = DmInv(2, 1);
  dFdx(7, 10) = DmInv(2, 2);

  // dF / du11
  dFdx(2, 11) = DmInv(2, 0);
  dFdx(5, 11) = DmInv(2, 1);
  dFdx(8, 11) = DmInv(2, 2);
}

void TetMeshDeformationModelInternal::compute_Bm(const ES::M3d &DmInv, double volume, M3x4d &Bm) const
{
  Bm.block<3, 3>(0, 1) = DmInv.transpose() * volume;
  Bm.col(0) = (Bm.col(1) + Bm.col(2) + Bm.col(3)) * -1.0;

  // LGI << "v1:\n"
  //   << Bm << '\n'
  //   << "v2:\n"
  //   << data->Bm << '\n'
  //  LGI << (Bm - data->Bm).norm();
}

void TetMeshDeformationModelInternal::computeF(const ES::M3d &Ds, const ES::M3d &FpInv, ES::M3d &Fe) const
{
  Fe = Ds * restDmInv * FpInv;
}

void TetMeshDeformationModelInternal::computeSVD(const ES::M3d &Fe, ES::M3d &U, ES::M3d &V, ES::V3d &S) const
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

void TetMeshDeformationModel::computeDs(const double xIn[12], double DsOut[9])
{
  ES::V12d x = Eigen::Map<const ES::V12d>(xIn);
  ES::M3d Ds;

  Ds.col(0) = x.segment<3>(3) - x.segment<3>(0);
  Ds.col(1) = x.segment<3>(6) - x.segment<3>(0);
  Ds.col(2) = x.segment<3>(9) - x.segment<3>(0);

  (Eigen::Map<ES::M3d>(DsOut)) = Ds;
}

void TetMeshDeformationModel::computeDm(const double XIn[12], double DmOut[9])
{
  Eigen::Map<const ES::V3d> X0(XIn), X1(XIn + 3), X2(XIn + 6), X3(XIn + 9);
  ES::M3d Dm;
  Dm.col(0) = X1 - X0;
  Dm.col(1) = X2 - X0;
  Dm.col(2) = X3 - X0;

  (Eigen::Map<ES::M3d>(DmOut)) = Dm;
}

double TetMeshDeformationModel::computeVolume(const double X0[3], const double X1[3],
  const double X2[3], const double X3[3])
{
  return Mesh::getTetVolume(asVec3d(X0), asVec3d(X1), asVec3d(X2), asVec3d(X3));
}

void TetMeshDeformationModel::compute_dF_dx(const double DmInvIn[9], double dFdxOut[9 * 12])
{
  ES::M9x12d dFdx = ES::M9x12d::Zero();

  ES::M3d DmInv = Eigen::Map<const ES::M3d>(DmInvIn);

  double v0 = -(DmInv(0, 0) + DmInv(1, 0) + DmInv(2, 0));
  double v1 = -(DmInv(0, 1) + DmInv(1, 1) + DmInv(2, 1));
  double v2 = -(DmInv(0, 2) + DmInv(1, 2) + DmInv(2, 2));

  // F is in column major
  // dF / dx0
  dFdx(0, 0) = v0;
  dFdx(3, 0) = v1;
  dFdx(6, 0) = v2;

  // dF / dx1
  dFdx(1, 1) = v0;
  dFdx(4, 1) = v1;
  dFdx(7, 1) = v2;

  // dF / dx2
  dFdx(2, 2) = v0;
  dFdx(5, 2) = v1;
  dFdx(8, 2) = v2;

  // dF / dx3
  dFdx(0, 3) = DmInv(0, 0);
  dFdx(3, 3) = DmInv(0, 1);
  dFdx(6, 3) = DmInv(0, 2);

  // dF / dx4
  dFdx(1, 4) = DmInv(0, 0);
  dFdx(4, 4) = DmInv(0, 1);
  dFdx(7, 4) = DmInv(0, 2);

  // dF / dx5
  dFdx(2, 5) = DmInv(0, 0);
  dFdx(5, 5) = DmInv(0, 1);
  dFdx(8, 5) = DmInv(0, 2);

  // dF / dx6
  dFdx(0, 6) = DmInv(1, 0);
  dFdx(3, 6) = DmInv(1, 1);
  dFdx(6, 6) = DmInv(1, 2);

  // dF / du7
  dFdx(1, 7) = DmInv(1, 0);
  dFdx(4, 7) = DmInv(1, 1);
  dFdx(7, 7) = DmInv(1, 2);

  // dF / du8
  dFdx(2, 8) = DmInv(1, 0);
  dFdx(5, 8) = DmInv(1, 1);
  dFdx(8, 8) = DmInv(1, 2);

  // dF / du9
  dFdx(0, 9) = DmInv(2, 0);
  dFdx(3, 9) = DmInv(2, 1);
  dFdx(6, 9) = DmInv(2, 2);

  // dF / du10
  dFdx(1, 10) = DmInv(2, 0);
  dFdx(4, 10) = DmInv(2, 1);
  dFdx(7, 10) = DmInv(2, 2);

  // dF / du11
  dFdx(2, 11) = DmInv(2, 0);
  dFdx(5, 11) = DmInv(2, 1);
  dFdx(8, 11) = DmInv(2, 2);

  (Eigen::Map<ES::M9x12d>(dFdxOut)) = dFdx;
}