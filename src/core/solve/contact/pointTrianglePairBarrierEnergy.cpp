#include "pointTrianglePairBarrierEnergy.h"

#include "EigenSupport.h"
#include "barrierFunction.h"

#include <tbb/concurrent_vector.h>
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
struct PointTrianglePairBarrierEnergyBuffer {
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

double PointTrianglePairBarrierEnergy::normalizePositiveZero(double dhat, double positiveZero) {
    return normalizePositiveZeroImpl(dhat, positiveZero);
}

PointTrianglePairBarrierEnergy::PointTrianglePairBarrierEnergy(
    int numPairs_, int numObjects_, const std::array<int, 4>* objectIDs_,
    const std::array<int, 4>* pointTrianglePairs_, const int* objectDOFOffsets_,
    const ES::SpMatD* sampleEmbeddingWeights_, const std::vector<double>* sampleWeights_,
    const double* constraintNormals_, const double* constraintBarycentricWeights_, double dhat_, double positiveZero_)
    : numPairs(numPairs_),
      numObjects(numObjects_),
      objectIDs(objectIDs_),
      pointTrianglePairs(pointTrianglePairs_),
      objectDOFOffsets(objectDOFOffsets_),
      sampleEmbeddingWeights(sampleEmbeddingWeights_),
      sampleWeights(sampleWeights_),
      constraintNormals(constraintNormals_),
      constraintBarycentricWeights(constraintBarycentricWeights_),
      dhat(dhat_),
      positiveZero(normalizePositiveZeroImpl(dhat_, positiveZero_)) {
    tbb::concurrent_vector<ES::TripletD> entries;

    tbb::parallel_for(0, numPairs, [&](int pi) {
        for (int vi = 0; vi < 4; ++vi) {
            const int sampleIDi = pointTrianglePairs[pi][vi];
            const int objIDi    = objectIDs[pi][vi];

            for (ES::SpMatD::InnerIterator it_i(sampleEmbeddingWeights[objIDi], sampleIDi * 3); it_i; ++it_i) {
                const int actualVtxi = static_cast<int>(it_i.col()) / 3;
                const int rowBase    = objectDOFOffsets[objIDi] + actualVtxi * 3;

                for (int vj = 0; vj < 4; ++vj) {
                    const int sampleIDj = pointTrianglePairs[pi][vj];
                    const int objIDj    = objectIDs[pi][vj];

                    for (ES::SpMatD::InnerIterator it_j(sampleEmbeddingWeights[objIDj], sampleIDj * 3); it_j; ++it_j) {
                        const int actualVtxj = static_cast<int>(it_j.col()) / 3;
                        const int colBase    = objectDOFOffsets[objIDj] + actualVtxj * 3;

                        for (int dofi = 0; dofi < 3; ++dofi) {
                            for (int dofj = 0; dofj < 3; ++dofj) {
                                entries.emplace_back(rowBase + dofi, colBase + dofj, 1.0);
                            }
                        }
                    }
                }
            }
        }
    });

    hessTemplate.resize(objectDOFOffsets[numObjects], objectDOFOffsets[numObjects]);
    hessTemplate.setFromTriplets(entries.begin(), entries.end());
    ES::buildEntryMap(hessTemplate, entryMap);

    allDOFs.resize(objectDOFOffsets[numObjects]);
    std::iota(allDOFs.begin(), allDOFs.end(), 0);
}

PointTrianglePairBarrierEnergyBuffer* PointTrianglePairBarrierEnergy::allocateBuffer() const {
    auto* newBuffer      = new PointTrianglePairBarrierEnergyBuffer;
    newBuffer->entryLocks = std::vector<tbb::spin_mutex>(static_cast<size_t>(hessTemplate.rows()));
    return newBuffer;
}

void PointTrianglePairBarrierEnergy::freeBuffer(PointTrianglePairBarrierEnergyBuffer* buffer) const {
    delete buffer;
}

ES::V3d PointTrianglePairBarrierEnergy::computePosition(ES::ConstRefVecXd x, int objID, int sampleID) const {
    ES::V3d position = ES::V3d::Zero();

    for (ES::SpMatD::InnerIterator it(sampleEmbeddingWeights[objID], sampleID * 3); it; ++it) {
        const int offset = objectDOFOffsets[objID] + static_cast<int>(it.col());
        const ES::V3d localX = x.segment<3>(offset);
        ES::V3d       mappedPosition;
        toPosFunc(localX, mappedPosition, offset);
        position += mappedPosition * it.value();
    }

    return position;
}

double PointTrianglePairBarrierEnergy::getPairWeight(int pairID) const {
    if (!sampleWeights) {
        return 1.0;
    }

    const int pointObjectID = objectIDs[pairID][0];
    const int pointSampleID = pointTrianglePairs[pairID][0];
    if (pointSampleID < 0 || pointSampleID >= static_cast<int>(sampleWeights[pointObjectID].size())) {
        return 1.0;
    }

    return sampleWeights[pointObjectID][pointSampleID];
}

double PointTrianglePairBarrierEnergy::func(ES::ConstRefVecXd x) const {
    for (double& value : buf->energyTLS) {
        value = 0.0;
    }

    tbb::parallel_for(
        0, numPairs,
        [&](int pi) {
            const double pairWeight = getPairWeight(pi);
            if (std::abs(pairWeight) < 1e-12) {
                return;
            }

            const ES::V3d normal = ES::Mp<const ES::V3d>(constraintNormals + pi * 3);
            const ES::V3d baryW  = ES::Mp<const ES::V3d>(constraintBarycentricWeights + pi * 3);

            ES::V3d positions[4];
            for (int localIdx = 0; localIdx < 4; ++localIdx) {
                positions[localIdx] = computePosition(x, objectIDs[pi][localIdx], pointTrianglePairs[pi][localIdx]);
            }

            const double distance =
                normal.dot(positions[0] - positions[1] * baryW[0] - positions[2] * baryW[1] - positions[3] * baryW[2]);
            if (distance >= dhat) {
                return;
            }

            buf->energyTLS.local() += pairWeight * BF::logBarrierEnergy(distance, dhat, positiveZero);
        },
        tbb::static_partitioner());

    return std::accumulate(buf->energyTLS.begin(), buf->energyTLS.end(), 0.0) * coeffAll;
}

void PointTrianglePairBarrierEnergy::gradient(ES::ConstRefVecXd x, ES::RefVecXd grad) const {
    grad.setZero();

    tbb::parallel_for(
        0, numPairs,
        [&](int pi) {
            const double pairWeight = getPairWeight(pi);
            if (std::abs(pairWeight) < 1e-12) {
                return;
            }

            const ES::V3d normal = ES::Mp<const ES::V3d>(constraintNormals + pi * 3);
            const ES::V3d baryW  = ES::Mp<const ES::V3d>(constraintBarycentricWeights + pi * 3);

            ES::V3d positions[4];
            for (int localIdx = 0; localIdx < 4; ++localIdx) {
                positions[localIdx] = computePosition(x, objectIDs[pi][localIdx], pointTrianglePairs[pi][localIdx]);
            }

            const double distance =
                normal.dot(positions[0] - positions[1] * baryW[0] - positions[2] * baryW[1] - positions[3] * baryW[2]);
            if (distance >= dhat) {
                return;
            }

            const double barrierGradient = pairWeight * BF::logBarrierGradient(distance, dhat, positiveZero);
            const double localScale[4]   = {1.0, -baryW[0], -baryW[1], -baryW[2]};

            for (int localIdx = 0; localIdx < 4; ++localIdx) {
                const ES::V3d localGrad = normal * (barrierGradient * localScale[localIdx]);
                const int     objID     = objectIDs[pi][localIdx];
                const int     sampleID  = pointTrianglePairs[pi][localIdx];

                for (ES::SpMatD::InnerIterator it(sampleEmbeddingWeights[objID], sampleID * 3); it; ++it) {
                    const int offset = objectDOFOffsets[objID] + static_cast<int>(it.col());
                    buf->entryLocks[offset].lock();
                    grad.segment<3>(offset) += localGrad * it.value();
                    buf->entryLocks[offset].unlock();
                }
            }
        },
        tbb::static_partitioner());

    grad *= coeffAll;
}

void PointTrianglePairBarrierEnergy::hessian(ES::ConstRefVecXd x, ES::SpMatD& hess) const {
    std::memset(hess.valuePtr(), 0, sizeof(double) * hess.nonZeros());

    tbb::parallel_for(
        0, numPairs,
        [&](int pi) {
            const double pairWeight = getPairWeight(pi);
            if (std::abs(pairWeight) < 1e-12) {
                return;
            }

            const ES::V3d normal = ES::Mp<const ES::V3d>(constraintNormals + pi * 3);
            const ES::V3d baryW  = ES::Mp<const ES::V3d>(constraintBarycentricWeights + pi * 3);

            ES::V3d positions[4];
            for (int localIdx = 0; localIdx < 4; ++localIdx) {
                positions[localIdx] = computePosition(x, objectIDs[pi][localIdx], pointTrianglePairs[pi][localIdx]);
            }

            const double distance =
                normal.dot(positions[0] - positions[1] * baryW[0] - positions[2] * baryW[1] - positions[3] * baryW[2]);
            if (distance >= dhat) {
                return;
            }

            const double barrierHessian = pairWeight * BF::logBarrierHessian(distance, dhat, positiveZero);
            const ES::M3d nnT           = ES::tensorProduct(normal, normal) * barrierHessian;
            const double localScale[4]  = {1.0, -baryW[0], -baryW[1], -baryW[2]};

            for (int vi = 0; vi < 4; ++vi) {
                const int objIDi    = objectIDs[pi][vi];
                const int sampleIDi = pointTrianglePairs[pi][vi];

                for (int vj = 0; vj < 4; ++vj) {
                    const int objIDj    = objectIDs[pi][vj];
                    const int sampleIDj = pointTrianglePairs[pi][vj];
                    const ES::M3d Kij   = nnT * (localScale[vi] * localScale[vj]);

                    for (ES::SpMatD::InnerIterator it_i(sampleEmbeddingWeights[objIDi], sampleIDi * 3); it_i; ++it_i) {
                        const int offsetI = objectDOFOffsets[objIDi] + static_cast<int>(it_i.col());

                        for (ES::SpMatD::InnerIterator it_j(sampleEmbeddingWeights[objIDj], sampleIDj * 3); it_j;
                             ++it_j) {
                            const int offsetJ = objectDOFOffsets[objIDj] + static_cast<int>(it_j.col());
                            const ES::M3d block = Kij * (it_i.value() * it_j.value());

                            for (int dofi = 0; dofi < 3; ++dofi) {
                                for (int dofj = 0; dofj < 3; ++dofj) {
                                    const int grow = offsetI + dofi;
                                    const int gcol = offsetJ + dofj;
                                    const auto entry = entryMap.find(std::make_pair(grow, gcol));
                                    if (entry == entryMap.end()) {
                                        continue;
                                    }

                                    buf->entryLocks[grow].lock();
                                    hess.valuePtr()[entry->second] += block(dofi, dofj);
                                    buf->entryLocks[grow].unlock();
                                }
                            }
                        }
                    }
                }
            }
        },
        tbb::static_partitioner());

    hess *= coeffAll;
}
