/*
author: Bohan Wang
copyright to USC
*/

#include "triangleMeshSelfContactHandler.h"
#include "triangleMeshSelfContactDetection.h"
#include "pointTrianglePairBarrierEnergy.h"
#include "pointTrianglePairCouplingEnergyWithCollision.h"

#include "contactEnergyUtilities.h"
#include "pgoLogging.h"
#include "geometryQuery.h"
#include "triangleSampler.h"
#include "basicAlgorithms.h"
#include "predicates.h"

#include <tbb/parallel_for.h>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/concurrent_vector.h>
#include <tbb/spin_mutex.h>

#include <atomic>
#include <queue>
#include <limits>
#include <stdexcept>
#include <unordered_set>
#include <chrono>

using namespace pgo;
using namespace pgo::Contact;
using namespace pgo::Mesh;
using namespace pgo::BasicAlgorithms;

namespace ES = pgo::EigenSupport;
using hclock = std::chrono::high_resolution_clock;

inline double dura(const hclock::time_point& t1, const hclock::time_point& t2) {
    return std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1e6;
}

namespace {

ES::V3d computeClosestPointFromBarycentric(const ES::V3d& va, const ES::V3d& vb, const ES::V3d& vc,
                                           const ES::V3d& barycentricWeight) {
    return va * barycentricWeight[0] + vb * barycentricWeight[1] + vc * barycentricWeight[2];
}

ES::V3d computeTriangleUnitNormal(const ES::V3d& va, const ES::V3d& vb, const ES::V3d& vc) {
    ES::V3d n = (vb - va).cross(vc - va);
    const double nNorm = n.norm();
    if (nNorm > 1e-12) {
        n /= nNorm;
    } else {
        n.setZero();
    }

    return n;
}

}  // namespace

namespace pgo::Contact {
class TriangleMeshSelfContactHandlerRuntimeData {
public:
    using CCDD = CCDKernel::TriangleCCD::CCDData;
    tbb::enumerable_thread_specific<std::vector<CCDD, tbb::cache_aligned_allocator<CCDD>>> colliingTrianglePairsTLS;

    struct SampleSeedRecord {
        int    seedTriangleID = -1;
        double bestDist2      = std::numeric_limits<double>::infinity();
    };

    struct ActiveSeedRecord {
        int    sampleID       = -1;
        int    seedTriangleID = -1;
        double bestDist2      = std::numeric_limits<double>::infinity();
    };

    std::vector<SampleSeedRecord> sampleSeedRecords;
    std::vector<int>              sampleVisited;
    std::vector<ActiveSeedRecord> activeSeedRecords;

    std::vector<std::array<int, 2>> contactedPointTrianglePairs;
    std::vector<int>                contactedPointTrianglePairsMask;

    struct LocalSearchBuffer {
        int                                         sampleID;
        std::priority_queue<std::pair<double, int>> Q;
        std::unordered_set<int>                     searchFilter;
    };

    tbb::enumerable_thread_specific<LocalSearchBuffer> localSearchBuf;
};
}  // namespace pgo::Contact

TriangleMeshSelfContactHandler::TriangleMeshSelfContactHandler(const std::vector<ES::V3d>& V,
                                                               const std::vector<ES::V3i>& T, int nDOFs,
                                                               int                        subdivideTriangle,
                                                               const std::vector<int>*    vertexEmbeddingIndices,
                                                               const std::vector<double>* vertexEmbeddingWeights)
    : vertices(V), triangles(T), surfaceMeshRef(vertices, triangles), totaln3(nDOFs) {
    rd = std::make_shared<TriangleMeshSelfContactHandlerRuntimeData>();

    SPDLOG_LOGGER_INFO(Logging::lgr(), "Computing triangle area...");

    std::vector<double> triangleAreas(triangles.size(), 0);
    tbb::parallel_for(
        0, (int)triangles.size(),
        [&](int trii) {
            Vec3i   tri    = triangles[trii];
            ES::V3d vtx[3] = {
                vertices[tri[0]],
                vertices[tri[1]],
                vertices[tri[2]],
            };

            triangleAreas[trii] = getTriangleArea(vtx[0], vtx[1], vtx[2]);
        },
        tbb::static_partitioner());

    SPDLOG_LOGGER_INFO(Logging::lgr(), "Sampling surface mesh...");

    triangleSampler = std::make_shared<TriangleSampler>(std::max(subdivideTriangle, 1));
    triangleSamples.assign(triangles.size(), std::vector<SampleInfo>());

    auto sortThree = [](const SampleID& idIn) -> SampleID {
        SampleID id = idIn;

        if (id[0] > id[1])
            std::swap(id[0], id[1]);

        if (id[0] > id[2])
            std::swap(id[0], id[2]);

        if (id[1] > id[2])
            std::swap(id[1], id[2]);

        return id;
    };

    // Sample surface
    tbb::parallel_for(
        0, (int)triangles.size(),
        [&](int trii) {
            ES::V3i tri    = triangles[trii];
            ES::V3d vtx[3] = {
                vertices[tri[0]],
                vertices[tri[1]],
                vertices[tri[2]],
            };

            triangleSampler->visitSample(vtx[0], vtx[1], vtx[2], triangleAreas[trii],
                                         [&](int i, int j, double area, const ES::V3d& baryWeight, const ES::V3d& pos) {
                                             SampleInfo sampleInfo;

                                             UTriKey key         = triangleSampler->getSampleKey(tri, i, j);
                                             sampleInfo.id[0]    = key[0];
                                             sampleInfo.id[1]    = key[1];
                                             sampleInfo.id[2]    = key[2];
                                             sampleInfo.idSorted = sortThree(sampleInfo.id);

                                             sampleInfo.pos          = pos;
                                             sampleInfo.w            = baryWeight;
                                             sampleInfo.sampleWeight = area;
                                             sampleInfo.coord        = Vec2i(i, j);
                                             sampleInfo.triangleID   = trii;
                                             triangleSamples[trii].push_back(sampleInfo);
                                         });
        },
        tbb::static_partitioner());

    // std::size_t count = 0;
    // sampleIDQueryTable.reserve(triangles.size());
    std::atomic<int>                                                                count(0);
    tbb::concurrent_unordered_map<SampleInfo, int, SampleInfoHash, SampleInfoEqual> sampleIDQueryTableCC;

    if (subdivideTriangle > 1) {
        tbb::parallel_for((int)0, (int)triangleSamples.size(), [&](int trii) {
            count.fetch_add((int)triangleSamples[trii].size());

            for (const auto& sample : triangleSamples[trii]) {
                sampleIDQueryTableCC.emplace(sample, 0);
            }

            if (trii % 1000 == 0)
                std::cout << trii << ' ' << std::flush;
        });
        std::cout << std::endl;

        int inc = 0;
        for (auto& pr : sampleIDQueryTableCC) {
            pr.second = inc++;
        }
    }
    // we want to keep the original vertex index
    else {
        for (auto it = triangleSamples.begin(); it != triangleSamples.end(); ++it) {
            count += (int)it->size();
            for (const auto& sample : *it) {
                // there must be two negative and 1 positive IDs
                int maxID = std::max({sample.id[0], sample.id[1], sample.id[2]});
                PGO_ALOG((int64_t)sample.id[0] + (int64_t)sample.id[1] + (int64_t)sample.id[2] == (int64_t)maxID - 2);

                sampleIDQueryTableCC.emplace(sample, maxID / 2);
            }
        }
    }

    sampleIDQueryTable.reserve(sampleIDQueryTableCC.size());
    for (const auto& pr : sampleIDQueryTableCC) {
        sampleIDQueryTable.emplace(pr.first, pr.second);
    }

    SPDLOG_LOGGER_INFO(Logging::lgr(), "Sampling done.");

    std::vector<tbb::spin_mutex> sampleLocks(sampleIDQueryTable.size());

    sampleTriangleIDs.assign(sampleIDQueryTable.size(), std::vector<int>());
    tbb::parallel_for(0, (int)triangleSamples.size(), [&](int tri) {
        for (const auto& sinfo : triangleSamples[tri]) {
            int sid = sampleIDQueryTable.at(sinfo);
            sampleLocks[sid].lock();
            sampleTriangleIDs[sid].push_back(tri);
            sampleLocks[sid].unlock();
        }
    });

    // sampleInfoAndIDs.assign(sampleIDQueryTable.begin(), sampleIDQueryTable.end());
    int numSamples = (int)sampleIDQueryTableCC.size();
    sampleInfoAndIDs.resize(numSamples);
    for (const auto& pr : sampleIDQueryTable) {
        if (pr.second >= 0 && pr.second < numSamples) {
            sampleInfoAndIDs[pr.second] = pr.first;
        } else {
            fprintf(stderr, "pr.second %d out of bounds! (numSamples: %d)\n", pr.second, numSamples);
        }
    }

    vertexID2SampleIDsLinear.resize(vertices.size());
    // double maxDist = 0;
    for (auto it = sampleInfoAndIDs.begin(); it != sampleInfoAndIDs.end(); ++it) {
        const auto& sample   = *it;
        int         sampleID = (int)(it - sampleInfoAndIDs.begin());

        if (triangleSampler->isCornerIndex(sample.coord[0], sample.coord[1])) {
            UTriKey key =
                triangleSampler->getFeatureKeyFromSampleKey(UTriKey(sample.id[0], sample.id[1], sample.id[2]));
            int vtxID = std::max({key[0], key[1], key[2]});
            if (vtxID >= 0 && vtxID < (int)vertices.size()) {
                vertexID2SampleIDs.emplace(vtxID, sampleID);
                vertexID2SampleIDsLinear[vtxID] = sampleID;
            } else {
                fprintf(stderr, "vtxID %d out of bounds for arrays sized %d!\n", vtxID, (int)vertices.size());
            }

            // Vec3d p = vertices[vtxID];
            // Vec3d q = sample.pos;
            // maxDist = std::max(len(p - q), maxDist);
        }
    }
    // LGI << maxDist;

    if (subdivideTriangle <= 1) {
        for (auto it = vertexID2SampleIDs.begin(); it != vertexID2SampleIDs.end(); ++it) {
            PGO_ALOG(it->first == it->second);
        }
    }

    SPDLOG_LOGGER_INFO(Logging::lgr(), "Sampling done:");
    SPDLOG_LOGGER_INFO(Logging::lgr(), "  #samples (non-dup): {}", sampleIDQueryTable.size());
    SPDLOG_LOGGER_INFO(Logging::lgr(), "  #samples: {}", count.load());
    SPDLOG_LOGGER_INFO(Logging::lgr(), "  #vtx: {}", vertices.size());

    // if embedding weights and indices are given
    if (vertexEmbeddingIndices && vertexEmbeddingWeights) {
        const int embeddingArity =
            validateAndGetVertexEmbeddingArity(*vertexEmbeddingIndices, *vertexEmbeddingWeights,
                                               static_cast<int>(vertices.size()), "self-contact");

        // we need to first compute an interpolation matrix
        tbb::concurrent_vector<ES::TripletD> entries;

        // for (auto it = sampleIDQueryTable.begin(); it != sampleIDQueryTable.end(); ++it) {
        tbb::parallel_for(0, (int)sampleInfoAndIDs.size(), [&](int si) {
            int triID = sampleInfoAndIDs[si].triangleID;

            for (int vj = 0; vj < 3; vj++) {
                int vid = triangles[triID][vj];

                const int embedOffset = vid * embeddingArity;
                if (embedOffset + embeddingArity > (int)vertexEmbeddingIndices->size()) {
                    continue;
                }

                for (int j = 0; j < embeddingArity; j++) {
                    int    embedVid = (*vertexEmbeddingIndices)[embedOffset + j];
                    double embedW   = (*vertexEmbeddingWeights)[embedOffset + j];

                    double wfinal = sampleInfoAndIDs[si].w[vj] * embedW;
                    if (std::abs(wfinal) < 1e-16)
                        continue;

                    for (int dofi = 0; dofi < 3; dofi++) {
                        entries.emplace_back(si * 3 + dofi, embedVid * 3 + dofi, wfinal);
                    }
                }
            }
        });

        interpolationMatrix.resize(sampleInfoAndIDs.size() * 3, nDOFs);
        interpolationMatrix.setFromTriplets(entries.begin(), entries.end());

        tbb::parallel_for(0, (int)interpolationMatrix.rows(), [&](int rowi) {
            double wAll = 0;
            for (ES::SpMatD::InnerIterator it(interpolationMatrix, rowi); it; ++it) {
                wAll += it.value();
            }

            for (ES::SpMatD::InnerIterator it(interpolationMatrix, rowi); it; ++it) {
                it.valueRef() /= wAll;
            }
        });
    }
    // if there is not embedding given, but there is still more samples than just vertices
    else if (subdivideTriangle > 1) {
        tbb::concurrent_vector<ES::TripletD> entries;

        // for (auto it = sampleInfoAndID.begin(); it != sampleInfoAndID.end(); ++it) {
        tbb::parallel_for(0, (int)sampleInfoAndIDs.size(), [&](int si) {
            int triID = sampleInfoAndIDs[si].triangleID;

            Vec3i tri = triangles[triID];
            for (int j = 0; j < 3; j++) {
                entries.emplace_back(si * 3, tri[j] * 3, sampleInfoAndIDs[si].w[j]);
                entries.emplace_back(si * 3 + 1, tri[j] * 3 + 1, sampleInfoAndIDs[si].w[j]);
                entries.emplace_back(si * 3 + 2, tri[j] * 3 + 2, sampleInfoAndIDs[si].w[j]);
            }
        });

        interpolationMatrix.resize(sampleInfoAndIDs.size() * 3, nDOFs);
        interpolationMatrix.setFromTriplets(entries.begin(), entries.end());
    }
    else {
        std::vector<ES::TripletD> entries;
        entries.reserve(sampleInfoAndIDs.size() * 3);

        for (int si = 0; si < static_cast<int>(sampleInfoAndIDs.size()); ++si) {
            entries.emplace_back(si * 3, si * 3, 1.0);
            entries.emplace_back(si * 3 + 1, si * 3 + 1, 1.0);
            entries.emplace_back(si * 3 + 2, si * 3 + 2, 1.0);
        }

        interpolationMatrix.resize(sampleInfoAndIDs.size() * 3, nDOFs);
        interpolationMatrix.setFromTriplets(entries.begin(), entries.end());
    }

    n3 = (int)vertices.size() * 3;
    restP.setZero(n3);
    for (int i = 0; i < static_cast<int>(vertices.size()); i++)
        restP.segment<3>(i * 3) = vertices[i];

    curP0 = restP;
    curP1 = restP;
    lastP = restP;

    selfContactDetection = std::make_shared<TriangleMeshSelfContactDetection>(surfaceMeshRef);

    samplen3    = (int)sampleInfoAndIDs.size() * 3;
    sampleRestP = ES::VXd::Zero(samplen3);

    computeSamplePosition(restP, sampleRestP);
    sampleCurP0 = sampleRestP;
    sampleCurP1 = sampleRestP;
    sampleLastP = sampleRestP;

    SPDLOG_LOGGER_INFO(Logging::lgr(), "Computing vertex weights...");

    sampleWeights.assign(sampleIDQueryTable.size(), 0);
    for (size_t trii = 0; trii < triangleSamples.size(); trii++) {
        for (const auto& sample : triangleSamples[trii]) {
            auto it = sampleIDQueryTable.find(sample);
            PGO_ALOG(it != sampleIDQueryTable.end());

            sampleWeights[it->second] += sample.sampleWeight;
        }
    }

    if (!sampleWeights.empty()) {
        double maxW = *std::max_element(sampleWeights.begin(), sampleWeights.end());
        if (maxW > 1e-12) {
            for (auto& w : sampleWeights)
                w /= maxW;
        }
    }
}

void TriangleMeshSelfContactHandler::clearActivePairData() {
    contactedTrianglePairs.clear();
    contactedTriangleIDs.clear();
    activeClosestPoints.resize(0);
    activeNormals.resize(0);
    activeDistances.resize(0);
    activeBarycentricWeights.resize(0);
    contactEnergyObjIDs.clear();
}

void TriangleMeshSelfContactHandler::finalizeActivePairsFromSeeds(int maxSearchingNumTriangles) {
    if (rd->activeSeedRecords.empty()) {
        clearActivePairData();
        SPDLOG_LOGGER_INFO(Logging::lgr(), "# contacted point-triangle pairs (final): 0");
        return;
    }

    rd->contactedPointTrianglePairs.resize(rd->activeSeedRecords.size());
    rd->contactedPointTrianglePairsMask.assign(rd->activeSeedRecords.size(), 0);

    std::vector<ES::V3d> closestPoints(rd->activeSeedRecords.size(), ES::V3d::Zero());
    std::vector<ES::V3d> closestNormals(rd->activeSeedRecords.size(), ES::V3d::Zero());
    std::vector<ES::V3d> closestBarycentricWeights(rd->activeSeedRecords.size(), ES::V3d::Zero());
    std::vector<double>  closestDistances(rd->activeSeedRecords.size(), 0.0);

    std::atomic<int> counter(0);
    tbb::parallel_for(0, (int)rd->activeSeedRecords.size(), [&](int activeIdx) {
        const auto& seedRecord = rd->activeSeedRecords[activeIdx];
        const int   sampleIdx  = seedRecord.sampleID;
        const int   triIdx     = seedRecord.seedTriangleID;

        if (sampleIdx < 0 || triIdx < 0) {
            return;
        }

        auto& localSearchBuf    = rd->localSearchBuf.local();
        localSearchBuf.sampleID = sampleIdx;

        std::priority_queue<std::pair<double, int>>& Q            = localSearchBuf.Q;
        std::unordered_set<int>&                     searchFilter = localSearchBuf.searchFilter;

        while (!Q.empty()) {
            Q.pop();
        }
        searchFilter.clear();

        const ES::V3d p = sampleCurP0.segment<3>(sampleIdx * 3);

        auto evaluateTriangle = [&](int candidateTriID, TriMeshBVTree::ClosestTriangleQueryResult& bestRet,
                                    bool initializeBest) -> double {
            const ES::V3d va = curP0.segment<3>(triangles[candidateTriID][0] * 3);
            const ES::V3d vb = curP0.segment<3>(triangles[candidateTriID][1] * 3);
            const ES::V3d vc = curP0.segment<3>(triangles[candidateTriID][2] * 3);

            TriangleWithCollisionInfo triInfo(va, vb, vc);

            int     feature = -1;
            ES::V3d barycentricWeight;
            const double dist2 =
                triInfo.distanceToPoint2(p, &feature, barycentricWeight.data(), barycentricWeight.data() + 1,
                                         barycentricWeight.data() + 2);
            const ES::V3d closestPt = computeClosestPointFromBarycentric(va, vb, vc, barycentricWeight);

            if (initializeBest || dist2 < bestRet.dist2) {
                bestRet.closestPosition = closestPt;
                bestRet.dist2           = dist2;
                bestRet.feature         = feature;
                bestRet.triBaryWeight   = barycentricWeight;
                bestRet.triID           = candidateTriID;
            }

            return dist2;
        };

        TriMeshBVTree::ClosestTriangleQueryResult closestSite;
        closestSite.dist2 = std::numeric_limits<double>::infinity();
        const double seedDist2 = evaluateTriangle(triIdx, closestSite, true);

        Q.emplace(-seedDist2, triIdx);
        searchFilter.emplace(triIdx);

        while ((int)searchFilter.size() < maxSearchingNumTriangles && !Q.empty()) {
            const double negDistSq   = Q.top().first;
            const int    curTriangle = Q.top().second;
            Q.pop();

            if (-negDistSq > closestSite.dist2 * 5.0) {
                break;
            }

            const Vec3i neighbor = selfContactDetection->getTriMeshNeighbor().getTriangleNeighbors(curTriangle);
            for (int ni = 0; ni < 3; ni++) {
                const int neighborTriID = neighbor[ni];
                if (neighborTriID < 0) {
                    continue;
                }
                if (searchFilter.find(neighborTriID) != searchFilter.end()) {
                    continue;
                }

                const double candidateDist2 = evaluateTriangle(neighborTriID, closestSite, false);
                searchFilter.emplace(neighborTriID);
                Q.emplace(-candidateDist2, neighborTriID);
            }
        }

        if (closestSite.triID < 0) {
            return;
        }

        bool goodSample = true;
        for (const auto& triangleSample : triangleSamples[closestSite.triID]) {
            auto it = sampleIDQueryTable.find(triangleSample);
            PGO_ALOG(it != sampleIDQueryTable.end());
            if (sampleIdx == it->second) {
                goodSample = false;
                break;
            }
        }

        if (excludedTriangles.size()) {
            for (int tri : sampleTriangleIDs[sampleIdx]) {
                if (std::binary_search(excludedTriangles.begin(), excludedTriangles.end(), tri)) {
                    goodSample = false;
                    break;
                }
            }
            if (goodSample &&
                std::binary_search(excludedTriangles.begin(), excludedTriangles.end(), closestSite.triID)) {
                goodSample = false;
            }
        }

        if (!goodSample) {
            return;
        }

        rd->contactedPointTrianglePairs[activeIdx][0]  = sampleIdx;
        rd->contactedPointTrianglePairs[activeIdx][1]  = closestSite.triID;
        rd->contactedPointTrianglePairsMask[activeIdx] = 1;

        closestPoints[activeIdx] = closestSite.closestPosition;
        closestBarycentricWeights[activeIdx] = closestSite.triBaryWeight;
        closestDistances[activeIdx] = std::sqrt(std::max(closestSite.dist2, 0.0));

        ES::V3d normal = p - closestSite.closestPosition;
        const double normalNorm = normal.norm();
        if (normalNorm > 1e-12) {
            normal /= normalNorm;
        } else {
            const ES::V3d va = curP0.segment<3>(triangles[closestSite.triID][0] * 3);
            const ES::V3d vb = curP0.segment<3>(triangles[closestSite.triID][1] * 3);
            const ES::V3d vc = curP0.segment<3>(triangles[closestSite.triID][2] * 3);
            normal           = computeTriangleUnitNormal(va, vb, vc);
        }
        closestNormals[activeIdx] = normal;

        counter.fetch_add(1);
    });

    clearActivePairData();
    const int finalPairCount = counter.load();
    if (finalPairCount == 0) {
        SPDLOG_LOGGER_INFO(Logging::lgr(), "# contacted point-triangle pairs (final): 0");
        return;
    }

    contactedTrianglePairs.reserve(finalPairCount);
    contactedTriangleIDs.reserve(finalPairCount);
    activeClosestPoints.resize(finalPairCount * 3);
    activeNormals.resize(finalPairCount * 3);
    activeDistances.resize(finalPairCount);
    activeBarycentricWeights.resize(finalPairCount * 3);

    int outIdx = 0;
    for (size_t activeIdx = 0; activeIdx < rd->contactedPointTrianglePairs.size(); activeIdx++) {
        if (rd->contactedPointTrianglePairsMask[activeIdx] == 0) {
            continue;
        }

        const int triID = rd->contactedPointTrianglePairs[activeIdx][1];
        std::array<int, 4> sids{rd->contactedPointTrianglePairs[activeIdx][0],
                                vertexID2SampleIDsLinear[triangles[triID][0]],
                                vertexID2SampleIDsLinear[triangles[triID][1]],
                                vertexID2SampleIDsLinear[triangles[triID][2]]};
        contactedTrianglePairs.emplace_back(sids);
        contactedTriangleIDs.emplace_back(triID);
        activeClosestPoints.segment<3>(outIdx * 3) = closestPoints[activeIdx];
        activeNormals.segment<3>(outIdx * 3)       = closestNormals[activeIdx];
        activeDistances[outIdx]                    = closestDistances[activeIdx];
        activeBarycentricWeights.segment<3>(outIdx * 3) = closestBarycentricWeights[activeIdx];
        outIdx++;
    }

    SPDLOG_LOGGER_INFO(Logging::lgr(), "# contacted point-triangle pairs (final): {}", contactedTrianglePairs.size());
}

void TriangleMeshSelfContactHandler::handleContactDCD(double distThreshold, int maxSearchingNumTriangles) {
    if (collidingTrianglePairs.empty()) {
        clearActivePairData();
        return;
    }

    auto& sampleSeedRecords = rd->sampleSeedRecords;
    auto& sampleVisited     = rd->sampleVisited;

    if (sampleSeedRecords.size() != sampleInfoAndIDs.size()) {
        sampleSeedRecords.resize(sampleInfoAndIDs.size());
    }
    std::fill(sampleSeedRecords.begin(), sampleSeedRecords.end(),
              TriangleMeshSelfContactHandlerRuntimeData::SampleSeedRecord{});

    if (sampleVisited.size() != sampleInfoAndIDs.size()) {
        sampleVisited.assign(sampleInfoAndIDs.size(), 0);
    } else {
        std::fill(sampleVisited.begin(), sampleVisited.end(), 0);
    }

    std::vector<tbb::spin_mutex> sampleLocks(sampleInfoAndIDs.size());

    auto tryUpdateSeed = [&](int sampleID, int seedTriangleID, const ES::V3d& sampleP, const ES::V3d& va,
                             const ES::V3d& vb, const ES::V3d& vc, const ES::V3d& normal) {
        const double depth = (sampleP - va).dot(normal);
        if (depth >= distThreshold) {
            return;
        }

        const double dist2 = getSquaredDistanceToTriangle(sampleP, va, vb, vc);
        tbb::spin_mutex::scoped_lock lock(sampleLocks[sampleID]);
        if (sampleVisited[sampleID] == 0 || dist2 < sampleSeedRecords[sampleID].bestDist2) {
            sampleVisited[sampleID]                = 1;
            sampleSeedRecords[sampleID].seedTriangleID = seedTriangleID;
            sampleSeedRecords[sampleID].bestDist2      = dist2;
        }
    };

    tbb::parallel_for(0, (int)collidingTrianglePairs.size(), [&](int ci) {
        const int triA = collidingTrianglePairs[ci].triA;
        const int triB = collidingTrianglePairs[ci].triB;

        const ES::V3d vtxA[3] = {curP0.segment<3>(triangles[triA][0] * 3), curP0.segment<3>(triangles[triA][1] * 3),
                                 curP0.segment<3>(triangles[triA][2] * 3)};
        const ES::V3d vtxB[3] = {curP0.segment<3>(triangles[triB][0] * 3), curP0.segment<3>(triangles[triB][1] * 3),
                                 curP0.segment<3>(triangles[triB][2] * 3)};

        const ES::V3d nA = computeTriangleUnitNormal(vtxA[0], vtxA[1], vtxA[2]);
        const ES::V3d nB = computeTriangleUnitNormal(vtxB[0], vtxB[1], vtxB[2]);

        for (const auto& sinfo : triangleSamples[triB]) {
            auto it = sampleIDQueryTable.find(sinfo);
            PGO_ALOG(it != sampleIDQueryTable.end());
            const int sampleID   = it->second;
            const ES::V3d sampleP = sampleCurP0.segment<3>(sampleID * 3);
            tryUpdateSeed(sampleID, triA, sampleP, vtxA[0], vtxA[1], vtxA[2], nA);
        }

        for (const auto& sinfo : triangleSamples[triA]) {
            auto it = sampleIDQueryTable.find(sinfo);
            PGO_ALOG(it != sampleIDQueryTable.end());
            const int sampleID   = it->second;
            const ES::V3d sampleP = sampleCurP0.segment<3>(sampleID * 3);
            tryUpdateSeed(sampleID, triB, sampleP, vtxB[0], vtxB[1], vtxB[2], nB);
        }
    });

    rd->activeSeedRecords.clear();
    for (size_t sampleID = 0; sampleID < sampleVisited.size(); sampleID++) {
        if (sampleVisited[sampleID] == 0) {
            continue;
        }

        rd->activeSeedRecords.push_back(
            {static_cast<int>(sampleID), sampleSeedRecords[sampleID].seedTriangleID, sampleSeedRecords[sampleID].bestDist2});
    }

    finalizeActivePairsFromSeeds(maxSearchingNumTriangles);
}

std::vector<std::array<ES::V3d, 3>> TriangleMeshSelfContactHandler::getCollidingTriangles(const double* u) const {
    std::vector<std::array<ES::V3d, 3>> collidingTriangles;

    for (const auto& info : collidingTrianglePairs) {
        std::array<ES::V3d, 3> tri;
        Vec3i                  triAi = triangles[info.triA];
        Vec3i                  triBi = triangles[info.triB];

        for (int i = 0; i < 3; i++) {
            tri[i] = vertices[triAi[i]] + asVec3d(u + triAi[i] * 3);
        }
        collidingTriangles.emplace_back(tri);

        for (int i = 0; i < 3; i++) {
            tri[i] = vertices[triBi[i]] + asVec3d(u + triBi[i] * 3);
        }
        collidingTriangles.emplace_back(tri);
    }

    return collidingTriangles;
}

std::vector<ES::V3d> TriangleMeshSelfContactHandler::getSamplePoints(const double* u) const {
    ES::VXd P  = restP + Eigen::Map<const ES::VXd>(u, n3);
    ES::VXd SP = sampleRestP;
    computeSamplePosition(P, SP);

    std::vector<ES::V3d> pts;
    for (ES::IDX i = 0; i < SP.size() / 3; i++) {
        pts.push_back(asVec3d(SP.data() + i * 3));
    }

    return pts;
}

std::vector<ES::V3d> TriangleMeshSelfContactHandler::getCollidedSamplePoints(const double* u) const {
    ES::VXd P  = restP + Eigen::Map<const ES::VXd>(u, n3);
    ES::VXd SP = sampleRestP;
    computeSamplePosition(P, SP);

    std::vector<ES::V3d> pts;

    for (const auto& quad : contactedTrianglePairs) {
        for (int j = 0; j < 4; j++)
            pts.push_back(asVec3d(SP.data() + quad[j] * 3));
    }

    return pts;
}

std::vector<int> TriangleMeshSelfContactHandler::getCollidedSampleAffectedVertices() const {
    std::vector<int> vertexIndices;

    for (const auto& quad : contactedTrianglePairs) {
        for (int j = 0; j < 4; j++) {
            for (ES::SpMatD::InnerIterator it(interpolationMatrix, quad[j] * 3); it; ++it) {
                vertexIndices.push_back((int)it.col() / 3);
            }
        }
    }

    sortAndDeduplicateWithErase(vertexIndices);

    return vertexIndices;
}

double TriangleMeshSelfContactHandler::computeEmbeddedAlphaUpperBound(ES::ConstRefVecXd currentU,
                                                                      ES::ConstRefVecXd du, double alphaSafety,
                                                                      double dSafe) const {
    if (currentU.size() != totaln3 || du.size() != totaln3) {
        throw std::invalid_argument("computeEmbeddedAlphaUpperBound expects full embedded displacement vectors.");
    }

    ES::VXd sampleCurrentU = interpolationMatrix * currentU;
    ES::VXd sampleDu       = interpolationMatrix * du;

    return computeSampleAlphaUpperBound(sampleCurrentU, sampleDu, alphaSafety, dSafe);
}

double TriangleMeshSelfContactHandler::computeSampleAlphaUpperBound(ES::ConstRefVecXd sampleCurrentU,
                                                                    ES::ConstRefVecXd sampleDu, double alphaSafety,
                                                                    double dSafe) const {
    if (sampleCurrentU.size() != samplen3 || sampleDu.size() != samplen3) {
        throw std::invalid_argument("computeSampleAlphaUpperBound expects sample-sized displacement vectors.");
    }

    if (contactedTrianglePairs.empty()) {
        return 1.0;
    }

    double alphaUpper         = 1.0;
    bool   hasClosingDirection = false;

    for (int pi = 0; pi < static_cast<int>(contactedTrianglePairs.size()); ++pi) {
        const auto& pair = contactedTrianglePairs[pi];

        const ES::V3d n = activeNormals.segment<3>(pi * 3);
        const ES::V3d w = activeBarycentricWeights.segment<3>(pi * 3);

        const ES::V3d p  = sampleRestP.segment<3>(pair[0] * 3) + sampleCurrentU.segment<3>(pair[0] * 3);
        const ES::V3d dp = sampleDu.segment<3>(pair[0] * 3);

        const ES::V3d t0  = sampleRestP.segment<3>(pair[1] * 3) + sampleCurrentU.segment<3>(pair[1] * 3);
        const ES::V3d t1  = sampleRestP.segment<3>(pair[2] * 3) + sampleCurrentU.segment<3>(pair[2] * 3);
        const ES::V3d t2  = sampleRestP.segment<3>(pair[3] * 3) + sampleCurrentU.segment<3>(pair[3] * 3);
        const ES::V3d dt0 = sampleDu.segment<3>(pair[1] * 3);
        const ES::V3d dt1 = sampleDu.segment<3>(pair[2] * 3);
        const ES::V3d dt2 = sampleDu.segment<3>(pair[3] * 3);

        const ES::V3d closestPoint = t0 * w[0] + t1 * w[1] + t2 * w[2];
        const ES::V3d closestDelta = dt0 * w[0] + dt1 * w[1] + dt2 * w[2];

        const double d0   = (p - closestPoint).dot(n);
        const double dDot = (dp - closestDelta).dot(n);

        if (d0 <= dSafe) {
            if (dDot > 0.0) {
                continue;
            }

            return 0.0;
        }

        if (dDot < 0.0) {
            hasClosingDirection = true;
            alphaUpper          = std::min(alphaUpper, (d0 - dSafe) / (-dDot));
            if (alphaUpper <= 0.0) {
                return 0.0;
            }
        }
    }

    if (!hasClosingDirection) {
        return 1.0;
    }

    return std::clamp(alphaUpper * std::clamp(alphaSafety, 0.0, 1.0), 0.0, 1.0);
}

void TriangleMeshSelfContactHandler::execute(const std::vector<ES::V3d>& p0, const std::vector<ES::V3d>& p1) {
    for (int i = 0; i < static_cast<int>(p0.size()); i++)
        curP0.segment<3>(i * 3) = ES::V3d(p0[i][0], p0[i][1], p0[i][2]);

    for (int i = 0; i < static_cast<int>(p1.size()); i++)
        curP1.segment<3>(i * 3) = ES::V3d(p1[i][0], p1[i][1], p1[i][2]);

    computeSamplePosition(curP0, sampleCurP0);
    computeSamplePosition(curP1, sampleCurP1);

    executeCCD();
}

void TriangleMeshSelfContactHandler::execute(const double* u0, const double* u1) {
    curP0.noalias() = restP + Eigen::Map<const ES::VXd>(u0, n3);
    curP1.noalias() = restP + Eigen::Map<const ES::VXd>(u1, n3);

    computeSamplePosition(curP0, sampleCurP0);
    computeSamplePosition(curP1, sampleCurP1);

    executeCCD();
}

void TriangleMeshSelfContactHandler::execute(const std::vector<ES::V3d>& p0) {
    for (int i = 0; i < static_cast<int>(p0.size()); i++)
        curP0.segment<3>(i * 3) = ES::V3d(p0[i][0], p0[i][1], p0[i][2]);

    computeSamplePosition(curP0, sampleCurP0);

    executeDCD();
}

void TriangleMeshSelfContactHandler::execute(const double* u0) {
    curP0.noalias() = restP + Eigen::Map<const ES::VXd>(u0, n3);

    computeSamplePosition(curP0, sampleCurP0);

    executeDCD();
}

void TriangleMeshSelfContactHandler::execute(const std::vector<ES::V3d>& p0, double activationDistance) {
    for (int i = 0; i < static_cast<int>(p0.size()); i++) {
        curP0.segment<3>(i * 3) = ES::V3d(p0[i][0], p0[i][1], p0[i][2]);
    }

    computeSamplePosition(curP0, sampleCurP0);

    if (activationDistance <= 0.0) {
        executeDCD();
        if (collidingTrianglePairs.empty()) {
            clearActivePairData();
            return;
        }

        handleContactDCD(0.0, 100);
        return;
    }

    executeNearContactDCD(activationDistance, 100);
}

void TriangleMeshSelfContactHandler::execute(const double* u0, double activationDistance) {
    curP0.noalias() = restP + Eigen::Map<const ES::VXd>(u0, n3);

    computeSamplePosition(curP0, sampleCurP0);

    if (activationDistance <= 0.0) {
        executeDCD();
        if (collidingTrianglePairs.empty()) {
            clearActivePairData();
            return;
        }

        handleContactDCD(0.0, 100);
        return;
    }

    executeNearContactDCD(activationDistance, 100);
}

void TriangleMeshSelfContactHandler::executeNearContactDCD(double activationDistance, int maxSearchingNumTriangles) {
    const hclock::time_point t1 = hclock::now();

    auto& sampleSeedRecords = rd->sampleSeedRecords;
    auto& sampleVisited     = rd->sampleVisited;

    if (sampleSeedRecords.size() != sampleInfoAndIDs.size()) {
        sampleSeedRecords.resize(sampleInfoAndIDs.size());
    }
    std::fill(sampleSeedRecords.begin(), sampleSeedRecords.end(),
              TriangleMeshSelfContactHandlerRuntimeData::SampleSeedRecord{});

    if (sampleVisited.size() != sampleInfoAndIDs.size()) {
        sampleVisited.assign(sampleInfoAndIDs.size(), 0);
    } else {
        std::fill(sampleVisited.begin(), sampleVisited.end(), 0);
    }

    selfContactDetection->execute(curP0.data(), activationDistance);
    const auto& candidateTrianglePairs = selfContactDetection->getCandidateTrianglePairs();

    collidingTrianglePairs.clear();
    if (candidateTrianglePairs.empty()) {
        clearActivePairData();
        lastCDTime = dura(t1, hclock::now());
        SPDLOG_LOGGER_INFO(Logging::lgr(), "Self CD time: {}", lastCDTime);
        return;
    }

    std::vector<tbb::spin_mutex> sampleLocks(sampleInfoAndIDs.size());
    const double activationDistance2 = activationDistance * activationDistance;

    auto tryUpdateSeed = [&](int sampleID, int seedTriangleID, const ES::V3d& sampleP, const ES::V3d& va,
                             const ES::V3d& vb, const ES::V3d& vc) {
        const double dist2 = getSquaredDistanceToTriangle(sampleP, va, vb, vc);
        if (dist2 >= activationDistance2) {
            return;
        }

        tbb::spin_mutex::scoped_lock lock(sampleLocks[sampleID]);
        if (sampleVisited[sampleID] == 0 || dist2 < sampleSeedRecords[sampleID].bestDist2) {
            sampleVisited[sampleID]                = 1;
            sampleSeedRecords[sampleID].seedTriangleID = seedTriangleID;
            sampleSeedRecords[sampleID].bestDist2      = dist2;
        }
    };

    tbb::parallel_for(0, (int)candidateTrianglePairs.size(), [&](int pairIdx) {
        const int triA = candidateTrianglePairs[pairIdx].first;
        const int triB = candidateTrianglePairs[pairIdx].second;

        const ES::V3d vtxA[3] = {curP0.segment<3>(triangles[triA][0] * 3), curP0.segment<3>(triangles[triA][1] * 3),
                                 curP0.segment<3>(triangles[triA][2] * 3)};
        const ES::V3d vtxB[3] = {curP0.segment<3>(triangles[triB][0] * 3), curP0.segment<3>(triangles[triB][1] * 3),
                                 curP0.segment<3>(triangles[triB][2] * 3)};

        for (const auto& sinfo : triangleSamples[triB]) {
            auto it = sampleIDQueryTable.find(sinfo);
            PGO_ALOG(it != sampleIDQueryTable.end());
            const int sampleID   = it->second;
            const ES::V3d sampleP = sampleCurP0.segment<3>(sampleID * 3);
            tryUpdateSeed(sampleID, triA, sampleP, vtxA[0], vtxA[1], vtxA[2]);
        }

        for (const auto& sinfo : triangleSamples[triA]) {
            auto it = sampleIDQueryTable.find(sinfo);
            PGO_ALOG(it != sampleIDQueryTable.end());
            const int sampleID   = it->second;
            const ES::V3d sampleP = sampleCurP0.segment<3>(sampleID * 3);
            tryUpdateSeed(sampleID, triB, sampleP, vtxB[0], vtxB[1], vtxB[2]);
        }
    });

    rd->activeSeedRecords.clear();
    for (size_t sampleID = 0; sampleID < sampleVisited.size(); sampleID++) {
        if (sampleVisited[sampleID] == 0) {
            continue;
        }

        rd->activeSeedRecords.push_back(
            {static_cast<int>(sampleID), sampleSeedRecords[sampleID].seedTriangleID, sampleSeedRecords[sampleID].bestDist2});
    }

    finalizeActivePairsFromSeeds(maxSearchingNumTriangles);

    lastCDTime = dura(t1, hclock::now());
    SPDLOG_LOGGER_INFO(Logging::lgr(), "Self CD time: {}", lastCDTime);
}

void TriangleMeshSelfContactHandler::executeCCD() {
    hclock::time_point t1 = hclock::now();

    const std::vector<std::pair<int, int>>* potentialColliingTrianglePairsPtr;

    selfContactDetection->execute(curP0.data(), curP1.data());
    potentialColliingTrianglePairsPtr = &selfContactDetection->getPotentialCollidingTrianglePairs();

    for (auto it = rd->colliingTrianglePairsTLS.begin(); it != rd->colliingTrianglePairsTLS.end(); ++it)
        it->clear();

    // for (const auto &triPair : potentialColliingTrianglePairs) {
    tbb::parallel_for((size_t)0, potentialColliingTrianglePairsPtr->size(),
                      [&](size_t tritriID) {
                          const auto& triPair = potentialColliingTrianglePairsPtr->at(tritriID);

                          CCDKernel::TriangleCCD::CCDData ccdData;
                          bool ret = CCDKernel::TriangleCCD::CCDTest(static_cast<int>(vertices.size()), curP0.data(),
                                                                     curP1.data(), static_cast<int>(triangles.size()),
                                                                     triangles.data(), triPair.first, triPair.second,
                                                                     CCDKernel::TriangleCCD::CCDM_3RDPARTY, &ccdData);

                          if (ret) {
                              rd->colliingTrianglePairsTLS.local().push_back(ccdData);
                          }
                      },
                      tbb::static_partitioner());

    if (keepPrevious == 0) {
        collidingTrianglePairs.clear();
        for (auto it = rd->colliingTrianglePairsTLS.begin(); it != rd->colliingTrianglePairsTLS.end(); ++it)
            collidingTrianglePairs.insert(collidingTrianglePairs.end(), it->begin(), it->end());
    } else {
        std::set<std::pair<int, int>> visited;
        for (const auto& pr : collidingTrianglePairs) {
            std::pair<int, int> newTriPair{pr.triA, pr.triB};
            if (newTriPair.first > newTriPair.second) {
                std::swap(newTriPair.first, newTriPair.second);
            }

            visited.emplace(newTriPair);
        }

        for (auto it = rd->colliingTrianglePairsTLS.begin(); it != rd->colliingTrianglePairsTLS.end(); ++it) {
            for (const auto& pr : *it) {
                std::pair<int, int> newTriPair{pr.triA, pr.triB};
                if (newTriPair.first > newTriPair.second) {
                    std::swap(newTriPair.first, newTriPair.second);
                }

                if (visited.find(newTriPair) != visited.end()) {
                    continue;
                }

                visited.emplace(newTriPair);
                collidingTrianglePairs.emplace_back(pr);
            }
        }
    }

    hclock::time_point t2 = hclock::now();
    lastCDTime            = dura(t1, t2);
}

void TriangleMeshSelfContactHandler::executeDCD() {
    hclock::time_point t1 = hclock::now();

    const std::vector<std::pair<int, int>>* potentialColliingTrianglePairsPtr;

    selfContactDetection->execute(curP0.data(), nullptr);
    potentialColliingTrianglePairsPtr = &selfContactDetection->getPotentialCollidingTrianglePairs();

    for (auto it = rd->colliingTrianglePairsTLS.begin(); it != rd->colliingTrianglePairsTLS.end(); ++it)
        it->clear();

    tbb::parallel_for((size_t)0, potentialColliingTrianglePairsPtr->size(),
                      [&](size_t tritriID) {
                          const auto& triPair = potentialColliingTrianglePairsPtr->at(tritriID);
                          Vec3i       triA    = triangles[triPair.first];
                          Vec3i       triB    = triangles[triPair.second];

                          bool ret = intersectTriTri(curP0.data() + triA[0] * 3, curP0.data() + triA[1] * 3,
                                                     curP0.data() + triA[2] * 3, curP0.data() + triB[0] * 3,
                                                     curP0.data() + triB[1] * 3, curP0.data() + triB[2] * 3);

                          if (ret) {
                              CCDKernel::TriangleCCD::CCDData ccdData;
                              ccdData.triA    = triPair.first;
                              ccdData.triB    = triPair.second;
                              ccdData.ccdCase = CCDKernel::TriangleCCD::CCDC_COLLIDED;
                              ccdData.t       = 0;

                              rd->colliingTrianglePairsTLS.local().push_back(ccdData);
                          }
                      },
                      tbb::static_partitioner());

    collidingTrianglePairs.clear();
    for (auto it = rd->colliingTrianglePairsTLS.begin(); it != rd->colliingTrianglePairsTLS.end(); ++it)
        collidingTrianglePairs.insert(collidingTrianglePairs.end(), it->begin(), it->end());

    hclock::time_point t2 = hclock::now();

    lastCDTime = dura(t1, t2);

    SPDLOG_LOGGER_INFO(Logging::lgr(), "# colliding triangles: {}", getCollidingTrianglePair().size());
    SPDLOG_LOGGER_INFO(Logging::lgr(), "Self CD time: {}", lastCDTime);
}

void TriangleMeshSelfContactHandler::setExcludedVertices(const std::vector<int>& excludedVertices) {
    std::vector<int> vtx = excludedVertices;
    sortAndDeduplicate(vtx);

    excludedTriangles.clear();

    for (int tri = 0; tri < (int)triangles.size(); tri++) {
        bool excluded = true;
        for (int j = 0; j < 3; j++) {
            if (std::binary_search(vtx.begin(), vtx.end(), triangles[tri][j]) == false) {
                excluded = false;
                break;
            }
        }

        if (excluded)
            excludedTriangles.push_back(tri);
    }
    SPDLOG_LOGGER_INFO(Logging::lgr(), "# excluded triangles: {}", excludedTriangles.size());
}

void TriangleMeshSelfContactHandler::computeSamplePosition(const ES::VXd& P, ES::VXd& SP) const {
    tbb::parallel_for(
        0, (int)sampleInfoAndIDs.size(),
        [&](int si) {
            int tri = sampleInfoAndIDs[si].triangleID;

            ES::V3d vtx[3] = {P.segment<3>(triangles[tri][0] * 3), P.segment<3>(triangles[tri][1] * 3),
                              P.segment<3>(triangles[tri][2] * 3)};

            // compute sample position
            const auto& sinfo     = sampleInfoAndIDs[si];
            ES::V3d     sampleP   = vtx[0] * sinfo.w[0] + vtx[1] * sinfo.w[1] + vtx[2] * sinfo.w[2];
            SP.segment<3>(si * 3) = sampleP;
        },
        tbb::static_partitioner());
}

std::shared_ptr<PointTrianglePairCouplingEnergyWithCollision> TriangleMeshSelfContactHandler::buildContactEnergy(
    int checkingNeighboringContact, int changingTriangle) {
    contactEnergyObjIDs.assign(contactedTrianglePairs.size(), std::array<int, 4>{0, 0, 0, 0});

    contactEnergyObjectDOFOffsets[0] = 0;
    contactEnergyObjectDOFOffsets[1] = totaln3;

    return std::make_shared<PointTrianglePairCouplingEnergyWithCollision>(
        (int)contactedTrianglePairs.size(), 1, contactEnergyObjIDs.data(), contactedTrianglePairs.data(),
        contactedTriangleIDs.data(), contactEnergyObjectDOFOffsets.data(), &surfaceMeshRef,
        &selfContactDetection->getTriMeshNeighbor(), &vertexID2SampleIDsLinear, &interpolationMatrix, &sampleWeights);
}

std::shared_ptr<PointTrianglePairBarrierEnergy> TriangleMeshSelfContactHandler::buildBarrierEnergy(
    double dhat, double positiveZero) {
    contactEnergyObjIDs.assign(contactedTrianglePairs.size(), std::array<int, 4>{0, 0, 0, 0});

    contactEnergyObjectDOFOffsets[0] = 0;
    contactEnergyObjectDOFOffsets[1] = totaln3;

    return std::make_shared<PointTrianglePairBarrierEnergy>(
        static_cast<int>(contactedTrianglePairs.size()), 1, contactEnergyObjIDs.data(), contactedTrianglePairs.data(),
        contactEnergyObjectDOFOffsets.data(), &interpolationMatrix, &sampleWeights, activeNormals.data(),
        activeBarycentricWeights.data(), dhat, positiveZero);
}
