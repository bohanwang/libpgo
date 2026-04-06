#include "pointPenetrationBarrierEnergy.h"

#include "EigenSupport.h"
#include "barrierFunction.h"
#include "pgoLogging.h"

#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <tbb/spin_mutex.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <numeric>

using namespace pgo;
using namespace pgo::Contact;

namespace ES = pgo::EigenSupport;
namespace BF = pgo::NonlinearOptimization::BarrierFunctions;

namespace pgo::Contact {
struct PointPenetrationBarrierEnergyBuffer {
public:
    tbb::enumerable_thread_specific<double> energyTLS;
    std::vector<tbb::spin_mutex>            entryLocks;
};

}  // namespace pgo::Contact

namespace {

double normalizePositiveZeroImpl(double dhat, double positiveZero) {
    const double recommended = std::max(1e-10, 1e-6 * dhat);
    const double safeUpper   = dhat * 0.5;

    if (positiveZero <= 0.0) {
        return std::min(recommended, safeUpper);
    }

    return std::min(positiveZero, safeUpper);
}

}  // namespace

double PointPenetrationBarrierEnergy::normalizePositiveZero(double dhat, double positiveZero) {
    return normalizePositiveZeroImpl(dhat, positiveZero);
}

PointPenetrationBarrierEnergy::PointPenetrationBarrierEnergy(
    int np, int na, const double* coeffs, const double* tgtPos, const double* nrms,
    const std::vector<std::vector<int>>& bIdx, const std::vector<std::vector<double>>& bWeights, double dhat_,
    double positiveZero_)
    : nAll(na),
      numPoints(np),
      constraintCoeffs(coeffs),
      constraintTargetPositions(tgtPos),
      constraintNormals(nrms),
      barycentricIdx(bIdx),
      barycentricWeights(bWeights),
      dhat(dhat_),
      positiveZero(normalizePositiveZeroImpl(dhat_, positiveZero_)) {
    std::vector<ES::TripletD> entries;
    for (int pi = 0; pi < numPoints; pi++) {
        for (int vi = 0; vi < static_cast<int>(barycentricIdx[pi].size()); vi++) {
            const int vidi = barycentricIdx[pi][vi];

            for (int vj = 0; vj < static_cast<int>(barycentricIdx[pi].size()); vj++) {
                const int vidj = barycentricIdx[pi][vj];

                for (int dofi = 0; dofi < 3; dofi++) {
                    for (int dofj = 0; dofj < 3; dofj++) {
                        entries.emplace_back(vidi * 3 + dofi, vidj * 3 + dofj, 1.0);
                    }
                }
            }
        }
    }

    hessianTemplate.resize(nAll, nAll);
    hessianTemplate.setFromTriplets(entries.begin(), entries.end());
    ES::buildEntryMap(hessianTemplate, hessianEntryMap);

    allDOFs.resize(nAll);
    std::iota(allDOFs.begin(), allDOFs.end(), 0);
}

PointPenetrationBarrierEnergyBuffer* PointPenetrationBarrierEnergy::allocateBuffer() const {
    auto* newBuffer      = new PointPenetrationBarrierEnergyBuffer;
    newBuffer->entryLocks = std::vector<tbb::spin_mutex>(nAll);
    return newBuffer;
}

void PointPenetrationBarrierEnergy::freeBuffer(PointPenetrationBarrierEnergyBuffer* buffer) const {
    delete buffer;
}

double PointPenetrationBarrierEnergy::func(ES::ConstRefVecXd u) const {
    for (auto& value : buf->energyTLS) {
        value = 0.0;
    }

    tbb::parallel_for(
        0, numPoints,
        [&](int ci) {
            if (std::abs(constraintCoeffs[ci]) < 1e-9) {
                return;
            }

            ES::V3d p = ES::V3d::Zero();
            for (int vi = 0; vi < static_cast<int>(barycentricIdx[ci].size()); vi++) {
                const int vid = barycentricIdx[ci][vi];
                ES::V3d    pp;
                posFunc(u.segment<3>(vid * 3), pp, vid * 3);
                p += pp * barycentricWeights[ci][vi];
            }

            const ES::V3d n        = ES::Mp<const ES::V3d>(constraintNormals + ci * 3);
            const ES::V3d p0       = ES::Mp<const ES::V3d>(constraintTargetPositions + ci * 3);
            const double  distance = (p - p0).dot(n);
            if (distance >= dhat) {
                return;
            }

            double& energyLocal = buf->energyTLS.local();
            energyLocal += constraintCoeffs[ci] * BF::logBarrierEnergy(distance, dhat, positiveZero);
        },
        tbb::static_partitioner());

    return std::accumulate(buf->energyTLS.begin(), buf->energyTLS.end(), 0.0) * coeffAll;
}

void PointPenetrationBarrierEnergy::gradient(ES::ConstRefVecXd u, ES::RefVecXd grad) const {
    grad.setZero();

    tbb::parallel_for(
        0, numPoints,
        [&](int ci) {
            if (std::abs(constraintCoeffs[ci]) < 1e-9) {
                return;
            }

            ES::V3d p = ES::V3d::Zero();
            for (int vi = 0; vi < static_cast<int>(barycentricIdx[ci].size()); vi++) {
                const int vid = barycentricIdx[ci][vi];
                ES::V3d    pp;
                posFunc(u.segment<3>(vid * 3), pp, vid * 3);
                p += pp * barycentricWeights[ci][vi];
            }

            const ES::V3d n        = ES::Mp<const ES::V3d>(constraintNormals + ci * 3);
            const ES::V3d p0       = ES::Mp<const ES::V3d>(constraintTargetPositions + ci * 3);
            const double  distance = (p - p0).dot(n);
            if (distance >= dhat) {
                return;
            }

            const double barrierGradient = constraintCoeffs[ci] * BF::logBarrierGradient(distance, dhat, positiveZero);
            const ES::V3d gradSample = n * barrierGradient;

            for (int vi = 0; vi < static_cast<int>(barycentricIdx[ci].size()); vi++) {
                const int    vid       = barycentricIdx[ci][vi];
                const double w         = barycentricWeights[ci][vi];
                const ES::V3d gradLocal = gradSample * w;

                buf->entryLocks[vid].lock();
                grad.segment<3>(vid * 3) += gradLocal;
                buf->entryLocks[vid].unlock();
            }
        },
        tbb::static_partitioner());

    grad *= coeffAll;
}

void PointPenetrationBarrierEnergy::hessian(ES::ConstRefVecXd u, ES::SpMatD& hess) const {
    std::memset(hess.valuePtr(), 0, sizeof(double) * hess.nonZeros());

    tbb::parallel_for(
        0, numPoints,
        [&](int ci) {
            if (std::abs(constraintCoeffs[ci]) < 1e-9) {
                return;
            }

            ES::V3d p = ES::V3d::Zero();
            for (int vi = 0; vi < static_cast<int>(barycentricIdx[ci].size()); vi++) {
                const int vid = barycentricIdx[ci][vi];
                ES::V3d    pp;
                posFunc(u.segment<3>(vid * 3), pp, vid * 3);
                p += pp * barycentricWeights[ci][vi];
            }

            const ES::V3d n        = ES::Mp<const ES::V3d>(constraintNormals + ci * 3);
            const ES::V3d p0       = ES::Mp<const ES::V3d>(constraintTargetPositions + ci * 3);
            const double  distance = (p - p0).dot(n);
            if (distance >= dhat) {
                return;
            }

            const double barrierHessian = constraintCoeffs[ci] * BF::logBarrierHessian(distance, dhat, positiveZero);
            const ES::M3d nnT           = ES::tensorProduct(n, n) * barrierHessian;

            for (int vi = 0; vi < static_cast<int>(barycentricIdx[ci].size()); vi++) {
                const double wi = barycentricWeights[ci][vi];

                for (int vj = 0; vj < static_cast<int>(barycentricIdx[ci].size()); vj++) {
                    const double wj     = barycentricWeights[ci][vj];
                    const ES::M3d hLocal = nnT * wi * wj;

                    for (int dofi = 0; dofi < 3; dofi++) {
                        for (int dofj = 0; dofj < 3; dofj++) {
                            const int grow = barycentricIdx[ci][vi] * 3 + dofi;
                            const int gcol = barycentricIdx[ci][vj] * 3 + dofj;

                            const auto entry = hessianEntryMap.find(std::make_pair(grow, gcol));
                            PGO_ALOG(entry != hessianEntryMap.end());

                            buf->entryLocks[grow].lock();
                            hess.valuePtr()[entry->second] += hLocal(dofi, dofj);
                            buf->entryLocks[grow].unlock();
                        }
                    }
                }
            }
        },
        tbb::static_partitioner());

    hess *= coeffAll;
}
