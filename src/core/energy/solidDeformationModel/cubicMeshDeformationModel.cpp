/*
author: OpenAI Codex
copyright to USC,MIT,NUS
*/

#include "cubicMeshDeformationModel.h"

#include "elasticModel3DDeformationGradient.h"
#include "plasticModel3DDeformationGradient.h"

#include "EigenSupport.h"
#include "pgoLogging.h"

#include <array>
#include <cmath>

namespace ES = pgo::EigenSupport;
using M3x8d  = Eigen::Matrix<double, 3, 8>;
using M9x24d = Eigen::Matrix<double, 9, 24>;

namespace pgo {
namespace SolidDeformationModel {
namespace {
constexpr int kNumVertices         = 8;
constexpr int kNumQuadraturePoints = 8;
constexpr double kQuadratureWeight = 0.125;
}  // namespace

class CubicMeshDeformationModelInternal {
public:
    struct QuadratureData {
        M3x8d  dN_dabc;  // transpose of the jacobian of shape functions w.r.t. reference coordinates
        ES::M3d DmInv; // Dm = dX/dabc = X * dN_dabc^T, where X is the 3x8 matrix of vertex positions. DmInv is the inverse of Dm.
        M3x8d  restBm; // rest configuration Bm matrix = weightDetJ * DmInv^T * dN_dabc = weightDetJ * J_N(X)^T
        M9x24d rest_dFdx; 
        double weightDetJ; // w * det(Dm)
    };

    ES::V3d restX[kNumVertices];
    QuadratureData quad[kNumQuadraturePoints];

    const ElasticModel3DDeformationGradient* elasticModel;
    const PlasticModel3DDeformationGradient* plasticModel;

    void fillShapeGradients(double alpha, double beta, double gamma, M3x8d& dN_dabc) const;
    void computeDm(const M3x8d& X, const M3x8d& dN_dabc, ES::M3d& Dm) const;
    void computeF(const M3x8d& x, const M3x8d& dN_dabc, const ES::M3d& DmInv, ES::M3d& F) const;
    void compute_dF_dx(const M3x8d& dN_dabc, const ES::M3d& DmInv, M9x24d& dFdx) const;
    void computeSVD(const ES::M3d& Fe, ES::M3d& U, ES::M3d& V, ES::V3d& S) const;

    inline double compute_dV_dai(double weightDetJ, double ddetA_dai) const { return weightDetJ * ddetA_dai; }
    inline double compute_d2V_daidaj(double weightDetJ, double d2detA_daidaj) const {
        return weightDetJ * d2detA_daidaj;
    }
    void compute_dFe_dai(const ES::M3d& Fref, const ES::M3d& dAInvdai, ES::M3d& dFdai) const;
    void compute_d2Fe_dai_daj(const ES::M3d& Fref, const ES::M3d& dAInvdaidaj, ES::M3d& d2Fdaidaj) const;
    void compute_d2Fe_dx_dai(const ES::M3d& dAInvdai, const M9x24d& rest_dFdx, M9x24d& d2Fdudai) const;
    void compute_dP_dai(const ES::M9d& dPdF, const ES::M3d& dFdai, ES::M3d& dPdai) const;
    double compute_dpsi_dai(const ES::M3d& Fref, const ES::M3d& dAInv_dai, const ES::M3d& P) const;
    double compute_d2psi_dai_daj(const ES::M3d& Fref, const ES::M3d& dAInv_dai, const ES::M3d& dAInv_daj,
                                 const ES::M3d& d2AInv_dai_daj, const ES::M3d& P, const ES::M9d& dPdF) const;
};

class CubicMeshDeformationModelCacheData : public DeformationModelCacheData {
public:
    static const int maxNumPlasticParam = 12;
    static const int maxNumElasticParam = 100;

    typedef Eigen::Matrix<double, maxNumPlasticParam, 1>                  Vp;
    typedef Eigen::Matrix<double, maxNumPlasticParam, maxNumPlasticParam> Mp;

    Vp      plasticParam;
    ES::M3d Fp, FpInv;
    double  detFp;

    Vp ddetA_da;
    Mp d2detA_da2;

    ES::M3d dAInv_dai[maxNumPlasticParam];
    ES::M3d d2AInv_dai_daj[maxNumPlasticParam][maxNumPlasticParam];

    ES::V3d x[kNumVertices];
    ES::M3d Fref[kNumQuadraturePoints];
    ES::M3d Fe[kNumQuadraturePoints];
    ES::M3d U[kNumQuadraturePoints], V[kNumQuadraturePoints];
    ES::V3d S[kNumQuadraturePoints];
    M9x24d  dFdx[kNumQuadraturePoints];
    M3x8d   Bm[kNumQuadraturePoints];

    Eigen::Matrix<double, maxNumElasticParam, 1> materialParam;
};

}  // namespace SolidDeformationModel
}  // namespace pgo

using namespace pgo::SolidDeformationModel;

CubicMeshDeformationModel::CubicMeshDeformationModel(const double restPositions[24], ElasticModel* elasticModel,
                                                     PlasticModel* plasticModel)
    : DeformationModel(elasticModel, plasticModel) {
    ind = new CubicMeshDeformationModelInternal;

    for (int vi = 0; vi < kNumVertices; vi++) {
        ind->restX[vi] = ES::V3d(restPositions[vi * 3 + 0], restPositions[vi * 3 + 1], restPositions[vi * 3 + 2]);
    }

    ind->elasticModel = dynamic_cast<const ElasticModel3DDeformationGradient*>(elasticModel);
    ind->plasticModel = dynamic_cast<const PlasticModel3DDeformationGradient*>(plasticModel);
    PGO_ALOG(ind->elasticModel != nullptr);
    PGO_ALOG(ind->plasticModel != nullptr);

    const Eigen::Map<const M3x8d> restX(restPositions);
    const double quadratureCoord[2] = {
        0.5 - 0.5 / std::sqrt(3.0),
        0.5 + 0.5 / std::sqrt(3.0),
    };

    int qid = 0;
    for (int ia = 0; ia < 2; ia++) {
        for (int ib = 0; ib < 2; ib++) {
            for (int ig = 0; ig < 2; ig++) {
                auto& quad = ind->quad[qid++];
                // get the transpose of the jacobian of shape functions w.r.t. reference coordinates
                ind->fillShapeGradients(quadratureCoord[ia], quadratureCoord[ib], quadratureCoord[ig], quad.dN_dabc);

                ES::M3d Dm;
                ind->computeDm(restX, quad.dN_dabc, Dm);
                
                // compute the rest configuration Bm matrix and the deformation gradient jacobian
                quad.DmInv      = Dm.fullPivLu().inverse();
                quad.weightDetJ = kQuadratureWeight * std::abs(Dm.determinant());
                quad.restBm     = quad.weightDetJ * quad.DmInv.transpose() * quad.dN_dabc;
                ind->compute_dF_dx(quad.dN_dabc, quad.DmInv, quad.rest_dFdx);
            }
        }
    }

    numMaterialLocations = kNumQuadraturePoints;
}

CubicMeshDeformationModel::~CubicMeshDeformationModel() {
    delete ind;
}

DeformationModelCacheData* CubicMeshDeformationModel::allocateCacheData() const {
    return new CubicMeshDeformationModelCacheData;
}

void CubicMeshDeformationModel::freeCacheData(DeformationModelCacheData* data) const {
    delete data;
}

void CubicMeshDeformationModel::prepareData(const double* x, const double* param, const double* materialParam,
                                            CacheData* cacheDataBase) const {
    CubicMeshDeformationModelCacheData* cacheData = static_cast<CubicMeshDeformationModelCacheData*>(cacheDataBase);

    for (int vi = 0; vi < kNumVertices; vi++) {
        cacheData->x[vi] = ES::V3d(x[vi * 3 + 0], x[vi * 3 + 1], x[vi * 3 + 2]);
    }

    ind->plasticModel->computeA(param, cacheData->Fp.data());
    ind->plasticModel->computeAInv(param, cacheData->FpInv.data());
    cacheData->detFp = ind->plasticModel->compute_detA(param);

    for (int i = 0; i < ind->plasticModel->getNumParameters(); i++) {
        cacheData->plasticParam[i] = param[i];
    }

    ind->plasticModel->compute_ddetA_da(cacheData->plasticParam.data(), cacheData->ddetA_da.data(),
                                        CubicMeshDeformationModelCacheData::maxNumPlasticParam);
    ind->plasticModel->compute_d2detA_da2(cacheData->plasticParam.data(), cacheData->d2detA_da2.data(),
                                          CubicMeshDeformationModelCacheData::maxNumPlasticParam);

    for (int i = 0; i < ind->plasticModel->getNumParameters(); i++) {
        ind->plasticModel->compute_dAInv_da(cacheData->plasticParam.data(), i, cacheData->dAInv_dai[i].data());
        for (int j = 0; j < ind->plasticModel->getNumParameters(); j++) {
            ind->plasticModel->compute_d2AInv_da2(cacheData->plasticParam.data(), i, j,
                                                  cacheData->d2AInv_dai_daj[i][j].data());
        }
    }

    for (int i = 0; ind->elasticModel && i < ind->elasticModel->getNumParameters(); i++) {
        cacheData->materialParam[i] = materialParam[i];
    }

    const Eigen::Map<const M3x8d> xMat(x);
    for (int qi = 0; qi < kNumQuadraturePoints; qi++) {
        const auto& quad = ind->quad[qi];

        ind->computeF(xMat, quad.dN_dabc, quad.DmInv, cacheData->Fref[qi]);
        cacheData->Fe[qi] = cacheData->Fref[qi] * cacheData->FpInv;
        ind->computeSVD(cacheData->Fe[qi], cacheData->U[qi], cacheData->V[qi], cacheData->S[qi]);

        ind->compute_dF_dx(quad.dN_dabc, quad.DmInv * cacheData->FpInv, cacheData->dFdx[qi]);
        cacheData->Bm[qi] = cacheData->detFp * cacheData->FpInv.transpose() * quad.restBm;
    }
}

double CubicMeshDeformationModel::computeEnergy(const CacheData* cacheDataBase) const {
    const CubicMeshDeformationModelCacheData* cacheData =
        static_cast<const CubicMeshDeformationModelCacheData*>(cacheDataBase);

    double energy = 0.0;
    for (int qi = 0; qi < kNumQuadraturePoints; qi++) {
        energy += ind->elasticModel->compute_psi(cacheData->materialParam.data(), cacheData->Fe[qi].data(),
                                                 cacheData->U[qi].data(), cacheData->V[qi].data(),
                                                 cacheData->S[qi].data()) *
                  ind->quad[qi].weightDetJ * cacheData->detFp;
    }

    return energy;
}

void CubicMeshDeformationModel::compute_dE_dx(const CacheData* cacheDataBase, double* grad) const {
    const CubicMeshDeformationModelCacheData* cacheData =
        static_cast<const CubicMeshDeformationModelCacheData*>(cacheDataBase);

    Eigen::Map<M3x8d> gradMap(grad);
    gradMap.setZero();

    for (int qi = 0; qi < kNumQuadraturePoints; qi++) {
        ES::M3d P;
        ind->elasticModel->compute_P(cacheData->materialParam.data(), cacheData->Fe[qi].data(), cacheData->U[qi].data(),
                                     cacheData->V[qi].data(), cacheData->S[qi].data(), P.data());
        gradMap.noalias() += P * cacheData->Bm[qi];
    }
}

void CubicMeshDeformationModel::compute_d2E_dx2(const CacheData* cacheDataBase, double* hess) const {
    const CubicMeshDeformationModelCacheData* cacheData =
        static_cast<const CubicMeshDeformationModelCacheData*>(cacheDataBase);

    Eigen::Map<ES::M24d> hessMap(hess);
    hessMap.setZero();

    for (int qi = 0; qi < kNumQuadraturePoints; qi++) {
        ES::M9d dPdF;
        ind->elasticModel->compute_dPdF(cacheData->materialParam.data(), cacheData->Fe[qi].data(),
                                        cacheData->U[qi].data(), cacheData->V[qi].data(), cacheData->S[qi].data(),
                                        dPdF.data());
        dPdF *= ind->quad[qi].weightDetJ * cacheData->detFp;
        hessMap.noalias() += cacheData->dFdx[qi].transpose() * dPdF * cacheData->dFdx[qi];
    }
}

void CubicMeshDeformationModel::compute_dE_da(const CacheData* cacheDataBase, double* grad) const {
    if (ind->plasticModel->getNumParameters() == 0) {
        return;
    }

    const CubicMeshDeformationModelCacheData* cacheData =
        static_cast<const CubicMeshDeformationModelCacheData*>(cacheDataBase);

    Eigen::Map<ES::VXd> gradMap(grad, ind->plasticModel->getNumParameters());
    gradMap.setZero();

    for (int qi = 0; qi < kNumQuadraturePoints; qi++) {
        const auto& quad = ind->quad[qi];
        const double psi = ind->elasticModel->compute_psi(cacheData->materialParam.data(), cacheData->Fe[qi].data(),
                                                          cacheData->U[qi].data(), cacheData->V[qi].data(),
                                                          cacheData->S[qi].data());
        ES::M3d P;
        ind->elasticModel->compute_P(cacheData->materialParam.data(), cacheData->Fe[qi].data(), cacheData->U[qi].data(),
                                     cacheData->V[qi].data(), cacheData->S[qi].data(), P.data());

        for (int i = 0; i < ind->plasticModel->getNumParameters(); i++) {
            const double dVda = ind->compute_dV_dai(quad.weightDetJ, cacheData->ddetA_da[i]);
            const double dpsiDa = ind->compute_dpsi_dai(cacheData->Fref[qi], cacheData->dAInv_dai[i], P);
            gradMap[i] += dVda * psi + quad.weightDetJ * cacheData->detFp * dpsiDa;
        }
    }
}

void CubicMeshDeformationModel::compute_d2E_da2(const CacheData* cacheDataBase, double* hess) const {
    if (ind->plasticModel->getNumParameters() == 0) {
        return;
    }

    const CubicMeshDeformationModelCacheData* cacheData =
        static_cast<const CubicMeshDeformationModelCacheData*>(cacheDataBase);

    Eigen::Map<ES::MXd> hessMap(hess, ind->plasticModel->getNumParameters(), ind->plasticModel->getNumParameters());
    hessMap.setZero();

    for (int qi = 0; qi < kNumQuadraturePoints; qi++) {
        const auto& quad = ind->quad[qi];
        const double vol = quad.weightDetJ * cacheData->detFp;
        const double psi = ind->elasticModel->compute_psi(cacheData->materialParam.data(), cacheData->Fe[qi].data(),
                                                          cacheData->U[qi].data(), cacheData->V[qi].data(),
                                                          cacheData->S[qi].data());

        ES::M3d P;
        ind->elasticModel->compute_P(cacheData->materialParam.data(), cacheData->Fe[qi].data(), cacheData->U[qi].data(),
                                     cacheData->V[qi].data(), cacheData->S[qi].data(), P.data());

        ES::M9d dPdF;
        ind->elasticModel->compute_dPdF(cacheData->materialParam.data(), cacheData->Fe[qi].data(),
                                        cacheData->U[qi].data(), cacheData->V[qi].data(), cacheData->S[qi].data(),
                                        dPdF.data());

        for (int i = 0; i < ind->plasticModel->getNumParameters(); i++) {
            const double dVdaI = ind->compute_dV_dai(quad.weightDetJ, cacheData->ddetA_da[i]);
            const double dpsiDaI = ind->compute_dpsi_dai(cacheData->Fref[qi], cacheData->dAInv_dai[i], P);
            for (int j = 0; j < ind->plasticModel->getNumParameters(); j++) {
                const double dVdaJ = ind->compute_dV_dai(quad.weightDetJ, cacheData->ddetA_da[j]);
                const double dpsiDaJ = ind->compute_dpsi_dai(cacheData->Fref[qi], cacheData->dAInv_dai[j], P);
                const double d2V =
                    ind->compute_d2V_daidaj(quad.weightDetJ, cacheData->d2detA_da2(i, j));
                const double d2Psi = ind->compute_d2psi_dai_daj(cacheData->Fref[qi], cacheData->dAInv_dai[i],
                                                                cacheData->dAInv_dai[j],
                                                                cacheData->d2AInv_dai_daj[i][j], P, dPdF);
                hessMap(i, j) += d2V * psi + dVdaI * dpsiDaJ + dVdaJ * dpsiDaI + vol * d2Psi;
            }
        }
    }
}

void CubicMeshDeformationModel::compute_d2E_dxda(const CacheData* cacheDataBase, double* hess) const {
    if (ind->plasticModel->getNumParameters() == 0) {
        return;
    }

    const CubicMeshDeformationModelCacheData* cacheData =
        static_cast<const CubicMeshDeformationModelCacheData*>(cacheDataBase);

    Eigen::Map<ES::MXd> hessMap(hess, 24, ind->plasticModel->getNumParameters());
    hessMap.setZero();

    for (int qi = 0; qi < kNumQuadraturePoints; qi++) {
        const auto& quad = ind->quad[qi];

        ES::M3d P;
        ind->elasticModel->compute_P(cacheData->materialParam.data(), cacheData->Fe[qi].data(), cacheData->U[qi].data(),
                                     cacheData->V[qi].data(), cacheData->S[qi].data(), P.data());

        const ES::V24d dpsi_dx = cacheData->dFdx[qi].transpose() * Eigen::Map<const ES::V9d>(P.data());

        ES::M9d dPdF;
        ind->elasticModel->compute_dPdF(cacheData->materialParam.data(), cacheData->Fe[qi].data(),
                                        cacheData->U[qi].data(), cacheData->V[qi].data(), cacheData->S[qi].data(),
                                        dPdF.data());

        const double vol = quad.weightDetJ * cacheData->detFp;

        for (int i = 0; i < ind->plasticModel->getNumParameters(); i++) {
            const double dVda = ind->compute_dV_dai(quad.weightDetJ, cacheData->ddetA_da[i]);

            ES::M3d dFda, dPda;
            ind->compute_dFe_dai(cacheData->Fref[qi], cacheData->dAInv_dai[i], dFda);
            ind->compute_dP_dai(dPdF, dFda, dPda);

            const ES::V24d temp = cacheData->dFdx[qi].transpose() * Eigen::Map<const ES::V9d>(dPda.data());

            M9x24d d2F_duda;
            ind->compute_d2Fe_dx_dai(cacheData->dAInv_dai[i], quad.rest_dFdx, d2F_duda);

            hessMap.col(i).noalias() +=
                dVda * dpsi_dx + vol * (temp + d2F_duda.transpose() * Eigen::Map<const ES::V9d>(P.data()));
        }
    }
}

void CubicMeshDeformationModel::compute_dE_db(const CacheData* cacheDataBase, double* grad) const {
    if (ind->elasticModel->getNumParameters() == 0) {
        return;
    }

    const CubicMeshDeformationModelCacheData* cacheData =
        static_cast<const CubicMeshDeformationModelCacheData*>(cacheDataBase);

    Eigen::Map<ES::VXd> gradMap(grad, ind->elasticModel->getNumParameters());
    gradMap.setZero();

    for (int qi = 0; qi < kNumQuadraturePoints; qi++) {
        const double vol = ind->quad[qi].weightDetJ * cacheData->detFp;
        for (int i = 0; i < ind->elasticModel->getNumParameters(); i++) {
            gradMap[i] +=
                vol * ind->elasticModel->compute_dpsi_dparam(cacheData->materialParam.data(), i, cacheData->Fe[qi].data(),
                                                             cacheData->U[qi].data(), cacheData->V[qi].data(),
                                                             cacheData->S[qi].data());
        }
    }
}

void CubicMeshDeformationModel::compute_d2E_db2(const CacheData* cacheDataBase, double* hess) const {
    if (ind->elasticModel->getNumParameters() == 0) {
        return;
    }

    const CubicMeshDeformationModelCacheData* cacheData =
        static_cast<const CubicMeshDeformationModelCacheData*>(cacheDataBase);

    Eigen::Map<ES::MXd> hessMap(hess, ind->elasticModel->getNumParameters(), ind->elasticModel->getNumParameters());
    hessMap.setZero();

    for (int qi = 0; qi < kNumQuadraturePoints; qi++) {
        const double vol = ind->quad[qi].weightDetJ * cacheData->detFp;
        for (int i = 0; i < ind->elasticModel->getNumParameters(); i++) {
            for (int j = 0; j < ind->elasticModel->getNumParameters(); j++) {
                hessMap(i, j) += vol * ind->elasticModel->compute_d2psi_dparam2(
                                           cacheData->materialParam.data(), i, j, cacheData->Fe[qi].data(),
                                           cacheData->U[qi].data(), cacheData->V[qi].data(), cacheData->S[qi].data());
            }
        }
    }
}

void CubicMeshDeformationModel::compute_d2E_dxdb(const CacheData* cacheDataBase, double* hess) const {
    if (ind->elasticModel->getNumParameters() == 0) {
        return;
    }

    const CubicMeshDeformationModelCacheData* cacheData =
        static_cast<const CubicMeshDeformationModelCacheData*>(cacheDataBase);

    Eigen::Map<ES::MXd> hessMap(hess, 24, ind->elasticModel->getNumParameters());
    hessMap.setZero();

    for (int qi = 0; qi < kNumQuadraturePoints; qi++) {
        const double vol = ind->quad[qi].weightDetJ * cacheData->detFp;
        for (int i = 0; i < ind->elasticModel->getNumParameters(); i++) {
            ES::M3d dPdb;
            ind->elasticModel->compute_dP_dparam(cacheData->materialParam.data(), i, cacheData->Fe[qi].data(),
                                                 cacheData->U[qi].data(), cacheData->V[qi].data(),
                                                 cacheData->S[qi].data(), dPdb.data());
            hessMap.col(i).noalias() +=
                vol * (cacheData->dFdx[qi].transpose() * Eigen::Map<const ES::V9d>(dPdb.data()));
        }
    }
}

void CubicMeshDeformationModel::compute_d2E_dadb(const CacheData* cacheDataBase, double* hess) const {
    if (ind->plasticModel->getNumParameters() == 0 || ind->elasticModel->getNumParameters() == 0) {
        return;
    }

    const CubicMeshDeformationModelCacheData* cacheData =
        static_cast<const CubicMeshDeformationModelCacheData*>(cacheDataBase);

    Eigen::Map<ES::MXd> hessMap(hess, ind->plasticModel->getNumParameters(), ind->elasticModel->getNumParameters());
    hessMap.setZero();

    for (int qi = 0; qi < kNumQuadraturePoints; qi++) {
        const auto& quad = ind->quad[qi];
        const double vol = quad.weightDetJ * cacheData->detFp;

        std::array<ES::M3d, CubicMeshDeformationModelCacheData::maxNumElasticParam> dPdb;
        for (int j = 0; j < ind->elasticModel->getNumParameters(); j++) {
            ind->elasticModel->compute_dP_dparam(cacheData->materialParam.data(), j, cacheData->Fe[qi].data(),
                                                 cacheData->U[qi].data(), cacheData->V[qi].data(),
                                                 cacheData->S[qi].data(), dPdb[j].data());
        }

        for (int i = 0; i < ind->plasticModel->getNumParameters(); i++) {
            const double dVda = ind->compute_dV_dai(quad.weightDetJ, cacheData->ddetA_da[i]);
            ES::M3d dFda;
            ind->compute_dFe_dai(cacheData->Fref[qi], cacheData->dAInv_dai[i], dFda);

            for (int j = 0; j < ind->elasticModel->getNumParameters(); j++) {
                const double dPsiDb =
                    ind->elasticModel->compute_dpsi_dparam(cacheData->materialParam.data(), j, cacheData->Fe[qi].data(),
                                                           cacheData->U[qi].data(), cacheData->V[qi].data(),
                                                           cacheData->S[qi].data());
                hessMap(i, j) += dVda * dPsiDb +
                                 vol * Eigen::Map<const ES::V9d>(dPdb[j].data()).dot(Eigen::Map<const ES::V9d>(dFda.data()));
            }
        }
    }
}

void CubicMeshDeformationModel::vonMisesStress(const CacheData* cacheDataBase, int& nPt, double* stresses) const {
    const CubicMeshDeformationModelCacheData* cacheData =
        static_cast<const CubicMeshDeformationModelCacheData*>(cacheDataBase);

    nPt = kNumQuadraturePoints;
    for (int qi = 0; qi < kNumQuadraturePoints; qi++) {
        ES::M3d P;
        ind->elasticModel->compute_P(cacheData->materialParam.data(), cacheData->Fe[qi].data(), cacheData->U[qi].data(),
                                     cacheData->V[qi].data(), cacheData->S[qi].data(), P.data());

        const double detF = cacheData->Fe[qi].determinant();
        ES::M3d cauchyStress = P * cacheData->Fe[qi].transpose() / detF;

        const double t1 = (cauchyStress(0, 0) - cauchyStress(1, 1)) *
                          (cauchyStress(0, 0) - cauchyStress(1, 1));
        const double t2 = (cauchyStress(1, 1) - cauchyStress(2, 2)) *
                          (cauchyStress(1, 1) - cauchyStress(2, 2));
        const double t3 = (cauchyStress(2, 2) - cauchyStress(0, 0)) *
                          (cauchyStress(2, 2) - cauchyStress(0, 0));
        const double t4 = 6.0 * (cauchyStress(1, 2) * cauchyStress(1, 2) + cauchyStress(2, 0) * cauchyStress(2, 0) +
                                  cauchyStress(0, 1) * cauchyStress(0, 1));
        stresses[qi] = std::sqrt((t1 + t2 + t3 + t4) * 0.5);
    }
}

void CubicMeshDeformationModel::maxStrain(const CacheData* cacheDataBase, int& nPt, double* stresses) const {
    const CubicMeshDeformationModelCacheData* cacheData =
        static_cast<const CubicMeshDeformationModelCacheData*>(cacheDataBase);

    nPt = kNumQuadraturePoints;
    for (int qi = 0; qi < kNumQuadraturePoints; qi++) {
        ES::M3d E = 0.5 * (cacheData->Fe[qi].transpose() * cacheData->Fe[qi] - ES::M3d::Identity());
        Eigen::SelfAdjointEigenSolver<ES::M3d> eigSolver(E);
        stresses[qi] = eigSolver.eigenvalues().maxCoeff();
    }
}

void CubicMeshDeformationModelInternal::fillShapeGradients(double alpha, double beta, double gamma,
                                                           M3x8d& dN_dabc) const {
    const double a0 = 1.0 - alpha;
    const double b0 = 1.0 - beta;
    const double g0 = 1.0 - gamma;

    dN_dabc.row(0) << -b0 * g0, b0 * g0, beta * g0, -beta * g0, -b0 * gamma, b0 * gamma, beta * gamma,
        -beta * gamma;
    dN_dabc.row(1) << -a0 * g0, -alpha * g0, alpha * g0, a0 * g0, -a0 * gamma, -alpha * gamma, alpha * gamma,
        a0 * gamma;
    dN_dabc.row(2) << -a0 * b0, -alpha * b0, -alpha * beta, -a0 * beta, a0 * b0, alpha * b0, alpha * beta,
        a0 * beta;
}

void CubicMeshDeformationModelInternal::computeDm(const M3x8d& X, const M3x8d& dN_dabc, ES::M3d& Dm) const {
    Dm.noalias() = X * dN_dabc.transpose();
}

void CubicMeshDeformationModelInternal::computeF(const M3x8d& x, const M3x8d& dN_dabc, const ES::M3d& DmInv,
                                                 ES::M3d& F) const {
    F.noalias() = x * dN_dabc.transpose() * DmInv;
}

void CubicMeshDeformationModelInternal::compute_dF_dx(const M3x8d& dN_dabc, const ES::M3d& Dm_inverse, M9x24d& dFdx) const {
    dFdx.setZero();

    const Eigen::Matrix<double, 8, 3> G = dN_dabc.transpose() * Dm_inverse;
    for (int vi = 0; vi < kNumVertices; vi++) {
        for (int dim = 0; dim < 3; dim++) {
            ES::M3d dF = ES::M3d::Zero();
            dF.row(dim) = G.row(vi);
            dFdx.col(vi * 3 + dim) = Eigen::Map<const ES::V9d>(dF.data());
        }
    }
}

void CubicMeshDeformationModelInternal::computeSVD(const ES::M3d& Fe, ES::M3d& U, ES::M3d& V, ES::V3d& S) const {
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

void CubicMeshDeformationModelInternal::compute_dFe_dai(const ES::M3d& Fref, const ES::M3d& dAInvdai,
                                                        ES::M3d& dFdai) const {
    dFdai = Fref * dAInvdai;
}

void CubicMeshDeformationModelInternal::compute_d2Fe_dai_daj(const ES::M3d& Fref, const ES::M3d& dAInvdaidaj,
                                                             ES::M3d& d2Fdaidaj) const {
    d2Fdaidaj = Fref * dAInvdaidaj;
}

void CubicMeshDeformationModelInternal::compute_d2Fe_dx_dai(const ES::M3d& dAInvdai, const M9x24d& rest_dFdx,
                                                            M9x24d& d2Fdudai) const {
    d2Fdudai.setZero();
    for (int k = 0; k < 24; k++) {
        ES::M3d d2Fedukdai = Eigen::Map<const ES::M3d>(rest_dFdx.data() + k * 9) * dAInvdai;
        d2Fdudai.col(k) = Eigen::Map<const ES::V9d>(d2Fedukdai.data());
    }
}

void CubicMeshDeformationModelInternal::compute_dP_dai(const ES::M9d& dPdF, const ES::M3d& dFdai,
                                                       ES::M3d& dPdai) const {
    Eigen::Map<ES::V9d>(dPdai.data()) = dPdF * Eigen::Map<const ES::V9d>(dFdai.data());
}

double CubicMeshDeformationModelInternal::compute_dpsi_dai(const ES::M3d& Fref, const ES::M3d& dAInv_dai,
                                                           const ES::M3d& P) const {
    ES::M3d dFe_dai;
    compute_dFe_dai(Fref, dAInv_dai, dFe_dai);
    return P.cwiseProduct(dFe_dai).sum();
}

double CubicMeshDeformationModelInternal::compute_d2psi_dai_daj(const ES::M3d& Fref, const ES::M3d& dAInv_dai,
                                                                const ES::M3d& dAInv_daj,
                                                                const ES::M3d& d2AInv_dai_daj, const ES::M3d& P,
                                                                const ES::M9d& dPdF) const {
    ES::M3d dFe_dai, dFe_daj;
    compute_dFe_dai(Fref, dAInv_dai, dFe_dai);
    compute_dFe_dai(Fref, dAInv_daj, dFe_daj);

    ES::M3d dP_daj;
    compute_dP_dai(dPdF, dFe_daj, dP_daj);

    ES::M3d d2F_dai_daj;
    compute_d2Fe_dai_daj(Fref, d2AInv_dai_daj, d2F_dai_daj);

    return dP_daj.cwiseProduct(dFe_dai).sum() + P.cwiseProduct(d2F_dai_daj).sum();
}
