/*
author: Bohan Wang
copyright to USC
*/

#include "triangleMeshSelfContactHandler.h"
#include "triangleMeshSelfContactDetection.h"
#include "pointTrianglePairCouplingEnergyWithCollision.h"

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
#include <tuple>
#include <queue>
#include <unordered_set>
#include <chrono>

using namespace pgo;
using namespace pgo::Contact;
using namespace pgo::Mesh;
using namespace pgo::BasicAlgorithms;

namespace ES = pgo::EigenSupport;
using hclock = std::chrono::high_resolution_clock;

inline double dura(const hclock::time_point &t1, const hclock::time_point &t2)
{
  return std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1e6;
}

namespace pgo::Contact
{
class TriangleMeshSelfContactHandlerRuntimeData
{
public:
  using CCDD = CCDKernel::TriangleCCD::CCDData;
  tbb::enumerable_thread_specific<std::vector<CCDD, tbb::cache_aligned_allocator<CCDD>>> colliingTrianglePairsTLS;

  std::vector<std::tuple<int, double>> sampleTriDepthAll;
  std::vector<int> sampleVisited;

  std::vector<std::tuple<int, int, double>> sampleTriDepthActive;

  std::vector<std::array<int, 2>> contactedPointTrianglePairs;
  std::vector<int> contactedPointTrianglePairsMask;

  struct LocalSearchBuffer
  {
    int sampleID;
    std::priority_queue<std::pair<double, int>> Q;
    std::unordered_set<int> searchFilter;
  };

  tbb::enumerable_thread_specific<LocalSearchBuffer> localSearchBuf;
};
}  // namespace pgo::Contact

TriangleMeshSelfContactHandler::TriangleMeshSelfContactHandler(const std::vector<Vec3d> &V, const std::vector<Vec3i> &T, int nDOFs,
  int subdivideTriangle, const std::vector<int> *vertexEmbeddingIndices, const std::vector<double> *vertexEmbeddingWeights):
  vertices(V),
  triangles(T), surfaceMeshRef(vertices, triangles),
  totaln3(nDOFs)
{
  rd = std::make_shared<TriangleMeshSelfContactHandlerRuntimeData>();

  SPDLOG_LOGGER_INFO(Logging::lgr(), "Computing triangle area...");

  std::vector<double> triangleAreas(triangles.size(), 0);
  tbb::parallel_for(
    0, (int)triangles.size(), [&](int trii) {
      Vec3i tri = triangles[trii];
      Vec3d vtx[3] = {
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

  auto sortThree = [](const SampleID &idIn) -> SampleID {
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
    0, (int)triangles.size(), [&](int trii) {
      Vec3i tri = triangles[trii];
      Vec3d vtx[3] = {
        vertices[tri[0]],
        vertices[tri[1]],
        vertices[tri[2]],
      };

      triangleSampler->visitSample(vtx[0], vtx[1], vtx[2], triangleAreas[trii],
        [&](int i, int j, double area, const Vec3d &baryWeight, const Vec3d &pos) {
          SampleInfo sampleInfo;

          UTriKey key = triangleSampler->getSampleKey(tri, i, j);
          sampleInfo.id[0] = key[0];
          sampleInfo.id[1] = key[1];
          sampleInfo.id[2] = key[2];
          sampleInfo.idSorted = sortThree(sampleInfo.id);

          sampleInfo.pos = pos;
          sampleInfo.w = baryWeight;
          sampleInfo.sampleWeight = area;
          sampleInfo.coord = Vec2i(i, j);
          sampleInfo.triangleID = trii;
          triangleSamples[trii].push_back(sampleInfo);
        });
    },
    tbb::static_partitioner());

  // std::size_t count = 0;
  // sampleIDQueryTable.reserve(triangles.size());
  std::atomic<int> count(0);
  tbb::concurrent_unordered_map<SampleInfo, int, SampleInfoHash, SampleInfoEqual> sampleIDQueryTableCC;

  if (subdivideTriangle > 1) {
    tbb::parallel_for((int)0, (int)triangleSamples.size(), [&](int trii) {
      count.fetch_add((int)triangleSamples[trii].size());

      for (const auto &sample : triangleSamples[trii]) {
        sampleIDQueryTableCC.emplace(sample, 0);
      }

      if (trii % 1000 == 0)
        std::cout << trii << ' ' << std::flush;
    });
    std::cout << std::endl;

    int inc = 0;
    for (auto &pr : sampleIDQueryTableCC) {
      pr.second = inc++;
    }
  }
  // we want to keep the original vertex index
  else {
    for (auto it = triangleSamples.begin(); it != triangleSamples.end(); ++it) {
      count += (int)it->size();
      for (const auto &sample : *it) {
        // there must be two negative and 1 positive IDs
        int maxID = std::max({ sample.id[0], sample.id[1], sample.id[2] });
        PGO_ALOG((int64_t)sample.id[0] + (int64_t)sample.id[1] + (int64_t)sample.id[2] == (int64_t)maxID - 2);

        sampleIDQueryTableCC.emplace(sample, maxID / 2);
      }
    }
  }

  sampleIDQueryTable.reserve(sampleIDQueryTableCC.size());
  for (const auto &pr : sampleIDQueryTableCC) {
    sampleIDQueryTable.emplace(pr.first, pr.second);
  }

  SPDLOG_LOGGER_INFO(Logging::lgr(), "Sampling done.");

  std::vector<tbb::spin_mutex> sampleLocks(sampleIDQueryTable.size());

  sampleTriangleIDs.assign(sampleIDQueryTable.size(), std::vector<int>());
  tbb::parallel_for(0, (int)triangleSamples.size(), [&](int tri) {
    for (const auto &sinfo : triangleSamples[tri]) {
      int sid = sampleIDQueryTable[sinfo];
      sampleLocks[sid].lock();
      sampleTriangleIDs[sid].push_back(tri);
      sampleLocks[sid].unlock();
    }
  });

  // sampleInfoAndIDs.assign(sampleIDQueryTable.begin(), sampleIDQueryTable.end());
  sampleInfoAndIDs.resize(sampleIDQueryTable.size());
  for (const auto &pr : sampleIDQueryTable) {
    sampleInfoAndIDs[pr.second] = pr.first;
  }

  vertexID2SampleIDsLinear.resize(vertices.size());
  // double maxDist = 0;
  for (auto it = sampleInfoAndIDs.begin(); it != sampleInfoAndIDs.end(); ++it) {
    const auto &sample = *it;
    int sampleID = (int)(it - sampleInfoAndIDs.begin());

    if (triangleSampler->isCornerIndex(sample.coord[0], sample.coord[1])) {
      UTriKey key = triangleSampler->getFeatureKeyFromSampleKey(UTriKey(sample.id[0], sample.id[1], sample.id[2]));
      int vtxID = std::max({ key[0], key[1], key[2] });
      vertexID2SampleIDs.emplace(vtxID, sampleID);
      vertexID2SampleIDsLinear[vtxID] = sampleID;

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
    // we need to first compute an interpolation matrix
    tbb::concurrent_vector<ES::TripletD> entries;

    // for (auto it = sampleIDQueryTable.begin(); it != sampleIDQueryTable.end(); ++it) {
    tbb::parallel_for(0, (int)sampleInfoAndIDs.size(), [&](int si) {
      const Vec3d &p = sampleInfoAndIDs[si].pos;
      int triID = sampleInfoAndIDs[si].triangleID;

      for (int vj = 0; vj < 3; vj++) {
        int vid = triangles[triID][vj];

        for (int j = 0; j < 4; j++) {
          int tetVid = vertexEmbeddingIndices->at(vid * 4 + j);
          double tetw = vertexEmbeddingWeights->at(vid * 4 + j);

          double wfinal = sampleInfoAndIDs[si].w[vj] * tetw;
          if (std::abs(wfinal) < 1e-16)
            continue;

          for (int dofi = 0; dofi < 3; dofi++) {
            entries.emplace_back(si * 3 + dofi, tetVid * 3 + dofi, wfinal);
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

    interpolationMatrix.resize(sampleInfoAndIDs.size() * 3, vertices.size());
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

  samplen3 = (int)sampleInfoAndIDs.size() * 3;
  sampleRestP = ES::VXd::Zero(samplen3);

  computeSamplePosition(restP, sampleRestP);
  sampleCurP0 = sampleRestP;
  sampleCurP1 = sampleRestP;
  sampleLastP = sampleRestP;

  SPDLOG_LOGGER_INFO(Logging::lgr(), "Computing vertex weights...");

  sampleWeights.assign(sampleIDQueryTable.size(), 0);
  for (size_t trii = 0; trii < triangleSamples.size(); trii++) {
    for (const auto &sample : triangleSamples[trii]) {
      auto it = sampleIDQueryTable.find(sample);
      PGO_ALOG(it != sampleIDQueryTable.end());

      sampleWeights[it->second] += sample.sampleWeight;
    }
  }

  double maxW = *std::max_element(sampleWeights.begin(), sampleWeights.end());
  for (auto &w : sampleWeights)
    w /= maxW;
}

void TriangleMeshSelfContactHandler::handleContactDCD(double distThreshold, int maxSearchingNumTriangles)
{
  if (collidingTrianglePairs.size()) {
    auto &sampleTriDepth = rd->sampleTriDepthAll;
    auto &sampleVisited = rd->sampleVisited;

    if (sampleTriDepth.size() == 0) {
      sampleTriDepth.resize(sampleInfoAndIDs.size());
    }

    if (sampleVisited.size() == 0) {
      sampleVisited.assign(sampleInfoAndIDs.size(), 0);
    }

    memset(sampleVisited.data(), 0, sampleVisited.size() * sizeof(int));

    // for (size_t ci = 0; ci < collidingTrianglePairs.size(); ci++) {
    tbb::parallel_for(0, (int)collidingTrianglePairs.size(), [&](int ci) {
      int triA = collidingTrianglePairs[ci].triA;
      int triB = collidingTrianglePairs[ci].triB;

      ES::V3d vtxA[3] = {
        curP0.segment<3>(triangles[triA][0] * 3),
        curP0.segment<3>(triangles[triA][1] * 3),
        curP0.segment<3>(triangles[triA][2] * 3)
      };

      ES::V3d vtxB[3] = {
        curP0.segment<3>(triangles[triB][0] * 3),
        curP0.segment<3>(triangles[triB][1] * 3),
        curP0.segment<3>(triangles[triB][2] * 3)
      };

      ES::V3d nA = (vtxA[1] - vtxA[0]).cross(vtxA[2] - vtxA[0]);
      nA.normalize();

      ES::V3d nB = (vtxB[1] - vtxB[0]).cross(vtxB[2] - vtxB[0]);
      nB.normalize();

      // for each sample on the triangle B
      for (int si = 0; si < triangleSamples[triB].size(); si++) {
        const SampleInfo &sinfo = triangleSamples[triB][si];
        auto it = sampleIDQueryTable.find(sinfo);
        PGO_ALOG(it != sampleIDQueryTable.end());
        int sampleID = it->second;

        // compute sample position
        // ES::V3d sampleP = vtxB[0] * sinfo.w[0] + vtxB[1] * sinfo.w[1] + vtxB[2] * sinfo.w[2];
        ES::VXd sampleP = sampleCurP0.segment<3>(sampleID * 3);

        // if a point is in contact
        double depth = (sampleP - vtxA[0]).dot(nA);

        // if the point is under the surface or
        //    above surface but with dist < distThreshold
        if (depth < distThreshold) {
          // we find the distance between the point and the triangle
          Vec3d p(sampleP.data());
          Vec3d va(vtxA[0].data()), vb(vtxA[1].data()), vc(vtxA[2].data());

          double dist2 = getSquaredDistanceToTriangle(p, va, vb, vc);

          // if we encounter this point before
          if (sampleVisited[sampleID]) {
            // if the old one has bigger dist
            if (std::get<1>(sampleTriDepth[sampleID]) > dist2) {
              std::get<0>(sampleTriDepth[sampleID]) = triA;
              std::get<1>(sampleTriDepth[sampleID]) = dist2;
            }

            //// if the old one has bigger depth
            // if (itt->second.second > fabs(depth)) {
            //   itt->second.first = triA;
            //   itt->second.second = fabs(depth);
            // }
          }
          else {
            sampleTriDepth[sampleID] = std::make_tuple(triA, dist2);
            sampleVisited[sampleID] = 1;
          }
        }  // end if depth < 0
      }

      // for each sample on the triangle A
      for (int si = 0; si < triangleSamples[triA].size(); si++) {
        const SampleInfo &sinfo = triangleSamples[triA][si];
        auto it = sampleIDQueryTable.find(sinfo);
        PGO_ALOG(it != sampleIDQueryTable.end());
        int sampleID = it->second;

        // compute sample position
        // ES::V3d sampleP = vtxA[0] * sinfo.w[0] + vtxA[1] * sinfo.w[1] + vtxA[2] * sinfo.w[2];
        ES::V3d sampleP = sampleCurP0.segment<3>(sampleID * 3);

        // if a point is in contact
        double depth = (sampleP - vtxB[0]).dot(nB);
        if (depth < distThreshold) {
          // we find the distance between the point and the triangle
          Vec3d va(vtxB[0].data()), vb(vtxB[1].data()), vc(vtxB[2].data());
          Vec3d p(sampleP.data());

          double dist2 = getSquaredDistanceToTriangle(p, va, vb, vc);

          // if we encounter this point before
          if (sampleVisited[sampleID]) {
            // if the old one has bigger depth
            if (std::get<1>(sampleTriDepth[sampleID]) > dist2) {
              std::get<0>(sampleTriDepth[sampleID]) = triB;
              std::get<1>(sampleTriDepth[sampleID]) = dist2;
            }

            //// if the old one has bigger depth
            // if (itt->second.second > fabs(depth)) {
            //   itt->second.first = triB;
            //   itt->second.second = fabs(depth);
            // }
          }
          else {
            sampleTriDepth[sampleID] = std::make_tuple(triB, dist2);
            sampleVisited[sampleID] = 1;
          }
        }  // end if depth < 0
      }
    });    // end for

    // gather all samples
    rd->sampleTriDepthActive.clear();
    for (size_t i = 0; i < sampleVisited.size(); i++) {
      if (sampleVisited[i]) {
        rd->sampleTriDepthActive.emplace_back((int)i, std::get<0>(sampleTriDepth[i]), std::get<1>(sampleTriDepth[i]));
      }
    }

    rd->contactedPointTrianglePairs.resize(rd->sampleTriDepthActive.size());
    rd->contactedPointTrianglePairsMask.assign(rd->sampleTriDepthActive.size(), 0);

    std::atomic<int> counter(0);
    tbb::parallel_for(0, (int)rd->sampleTriDepthActive.size(), [&](int si) {
      // for (int si = 0; si < (int)rd->sampleTriDepthActive.size(); si++) {
      int sampleIdx = std::get<0>(rd->sampleTriDepthActive[si]);
      int triIdx = std::get<1>(rd->sampleTriDepthActive[si]);

      auto &localSearchBuf = rd->localSearchBuf.local();
      localSearchBuf.sampleID = sampleIdx;

      std::priority_queue<std::pair<double, int>> &Q = localSearchBuf.Q;
      std::unordered_set<int> &searchFilter = localSearchBuf.searchFilter;

      while (!Q.empty())
        Q.pop();

      searchFilter.clear();

      // get sample position and triangle positions
      Vec3d p(sampleCurP0.data() + sampleIdx * 3);
      Vec3d va(curP0.data() + triangles[triIdx][0] * 3), vb(curP0.data() + triangles[triIdx][1] * 3), vc(curP0.data() + triangles[triIdx][2] * 3);

      int feature;
      Vec3d w;
      Vec3d closestPt;

      // double dist2 = getSquaredDistanceToTriangle(p, va, vb, vc, feature, closestPt, w);
      TriangleWithCollisionInfo triInfo(va, vb, vc);

      // evaluate the distance
      double dist2 = triInfo.distanceToPoint2(p, &feature, w.data(), w.data() + 1, w.data() + 2);

      // calculate the min distance in the deformed shape
      TriMeshBVTree::ClosestTriangleQueryResult closestSite;
      closestSite.closestPosition = closestPt;
      closestSite.dist2 = dist2;
      closestSite.feature = feature;
      closestSite.triBaryWeight = w;
      closestSite.triID = triIdx;

      Q.emplace(-closestSite.dist2, closestSite.triID);
      searchFilter.emplace(closestSite.triID);

      while ((int)searchFilter.size() < maxSearchingNumTriangles && !Q.empty()) {
        double negDistSq = Q.top().first;
        int curTriangle = Q.top().second;

        Q.pop();

        // std::cout << "(" << -negDistSq << ',' << curTriangle << ") ";

        if (-negDistSq > closestSite.dist2 * 5.0)
          break;

        Vec3i neighbor = selfContactDetection->getTriMeshNeighbor().getTriangleNeighbors(curTriangle);
        for (int ni = 0; ni < 3; ni++) {
          int triId = neighbor[ni];

          // if no neighbors
          if (triId < 0)
            continue;

          // if this triangle is searched
          if (searchFilter.find(triId) != searchFilter.end())
            continue;

          // compute distance to the triangle
          va = asVec3d(curP0.data() + triangles[triId][0] * 3);
          vb = asVec3d(curP0.data() + triangles[triId][1] * 3);
          vc = asVec3d(curP0.data() + triangles[triId][2] * 3);

          TriangleWithCollisionInfo triInfo(va, vb, vc);

          // evaluate the distance
          double dist2 = triInfo.distanceToPoint2(p, &feature, w.data(), w.data() + 1, w.data() + 2);
          // dist2 = getSquaredDistanceToTriangle(p, va, vb, vc, feature, closestPt, w);

          if (dist2 < closestSite.dist2) {
            closestSite.closestPosition = closestPt;
            closestSite.dist2 = dist2;
            closestSite.feature = feature;
            closestSite.triBaryWeight = w;
            closestSite.triID = triId;
          }

          searchFilter.emplace(triId);
          Q.emplace(-dist2, triId);
        }
      }
      // std::cout << std::endl;

      // if the closest triangle is not itself
      Vec3i selectedTriangleVertices = surfaceMeshRef.tri(closestSite.triID);
      bool goodSample = true;

      for (size_t si = 0; si < triangleSamples[closestSite.triID].size(); si++) {
        auto it = sampleIDQueryTable.find(triangleSamples[closestSite.triID][si]);
        PGO_ALOG(it != sampleIDQueryTable.end());

        if (sampleIdx == it->second) {
          goodSample = false;
          break;
        }
      }

      for (int tri : sampleTriangleIDs[sampleIdx])
        if (excludedTriangles.size() && std::binary_search(excludedTriangles.begin(), excludedTriangles.end(), tri) == true)
          goodSample = false;

      if (excludedTriangles.size() && std::binary_search(excludedTriangles.begin(), excludedTriangles.end(), closestSite.triID) == true)
        goodSample = false;

      if (goodSample) {
        rd->contactedPointTrianglePairs[si][0] = sampleIdx;
        rd->contactedPointTrianglePairs[si][1] = closestSite.triID;
        rd->contactedPointTrianglePairsMask[si] = 1;

        counter.fetch_add(1);
      }
    });

    contactedTrianglePairs.clear();
    contactedTrianglePairs.reserve(counter);

    contactedTriangleIDs.clear();
    contactedTriangleIDs.reserve(counter);

    for (size_t si = 0; si < rd->contactedPointTrianglePairs.size(); si++) {
      if (rd->contactedPointTrianglePairsMask[si] == 0)
        continue;

      std::array<int, 4> sids{
        rd->contactedPointTrianglePairs[si][0],
        vertexID2SampleIDs[triangles[rd->contactedPointTrianglePairs[si][1]][0]],
        vertexID2SampleIDs[triangles[rd->contactedPointTrianglePairs[si][1]][1]],
        vertexID2SampleIDs[triangles[rd->contactedPointTrianglePairs[si][1]][2]]
      };
      contactedTrianglePairs.emplace_back(sids);
      contactedTriangleIDs.emplace_back(rd->contactedPointTrianglePairs[si][1]);
    }

    SPDLOG_LOGGER_INFO(Logging::lgr(), "# contacted point-triangle pairs (final): {}", contactedTrianglePairs.size());
  }
  else {
    contactedTrianglePairs.clear();
  }
}

std::vector<std::array<Vec3d, 3>> TriangleMeshSelfContactHandler::getCollidingTriangles(const double *u) const
{
  std::vector<std::array<Vec3d, 3>> collidingTriangles;

  for (const auto &info : collidingTrianglePairs) {
    std::array<Vec3d, 3> tri;
    Vec3i triAi = triangles[info.triA];
    Vec3i triBi = triangles[info.triB];

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

std::vector<Vec3d> TriangleMeshSelfContactHandler::getSamplePoints(const double *u) const
{
  ES::VXd P = restP + Eigen::Map<const ES::VXd>(u, n3);
  ES::VXd SP = sampleRestP;
  computeSamplePosition(P, SP);

  std::vector<Vec3d> pts;
  for (ES::IDX i = 0; i < SP.size() / 3; i++) {
    pts.push_back(asVec3d(SP.data() + i * 3));
  }

  return pts;
}

std::vector<Vec3d> TriangleMeshSelfContactHandler::getCollidedSamplePoints(const double *u) const
{
  ES::VXd P = restP + Eigen::Map<const ES::VXd>(u, n3);
  ES::VXd SP = sampleRestP;
  computeSamplePosition(P, SP);

  std::vector<Vec3d> pts;

  for (const auto &quad : contactedTrianglePairs) {
    for (int j = 0; j < 4; j++)
      pts.push_back(asVec3d(SP.data() + quad[j] * 3));
  }

  return pts;
}

std::vector<int> TriangleMeshSelfContactHandler::getCollidedSampleAffectedVertices() const
{
  std::vector<int> vertexIndices;

  for (const auto &quad : contactedTrianglePairs) {
    for (int j = 0; j < 4; j++) {
      for (ES::SpMatD::InnerIterator it(interpolationMatrix, quad[j] * 3); it; ++it) {
        vertexIndices.push_back((int)it.col() / 3);
      }
    }
  }

  sortAndDeduplicateWithErase(vertexIndices);

  return vertexIndices;
}

void TriangleMeshSelfContactHandler::execute(const std::vector<Vec3d> &p0, const std::vector<Vec3d> &p1)
{
  for (int i = 0; i < static_cast<int>(p0.size()); i++)
    curP0.segment<3>(i * 3) = ES::V3d(p0[i][0], p0[i][1], p0[i][2]);

  for (int i = 0; i < static_cast<int>(p1.size()); i++)
    curP1.segment<3>(i * 3) = ES::V3d(p1[i][0], p1[i][1], p1[i][2]);

  computeSamplePosition(curP0, sampleCurP0);
  computeSamplePosition(curP1, sampleCurP1);

  executeCCD();
}

void TriangleMeshSelfContactHandler::execute(const double *u0, const double *u1)
{
  curP0.noalias() = restP + Eigen::Map<const ES::VXd>(u0, n3);
  curP1.noalias() = restP + Eigen::Map<const ES::VXd>(u1, n3);

  computeSamplePosition(curP0, sampleCurP0);
  computeSamplePosition(curP1, sampleCurP1);

  executeCCD();
}

void TriangleMeshSelfContactHandler::execute(const std::vector<Vec3d> &p0)
{
  for (int i = 0; i < static_cast<int>(p0.size()); i++)
    curP0.segment<3>(i * 3) = ES::V3d(p0[i][0], p0[i][1], p0[i][2]);

  computeSamplePosition(curP0, sampleCurP0);

  executeDCD();
}

void TriangleMeshSelfContactHandler::execute(const double *u0)
{
  curP0.noalias() = restP + Eigen::Map<const ES::VXd>(u0, n3);

  computeSamplePosition(curP0, sampleCurP0);

  executeDCD();
}

void TriangleMeshSelfContactHandler::executeCCD()
{
  hclock::time_point t1 = hclock::now();

  const std::vector<std::pair<int, int>> *potentialColliingTrianglePairsPtr;

  selfContactDetection->execute(curP0.data(), curP1.data());
  potentialColliingTrianglePairsPtr = &selfContactDetection->getPotentialCollidingTrianglePairs();

  for (auto it = rd->colliingTrianglePairsTLS.begin(); it != rd->colliingTrianglePairsTLS.end(); ++it)
    it->clear();

  // for (const auto &triPair : potentialColliingTrianglePairs) {
  tbb::parallel_for((size_t)0, potentialColliingTrianglePairsPtr->size(), [&](size_t tritriID) {
    const auto &triPair = potentialColliingTrianglePairsPtr->at(tritriID);

    CCDKernel::TriangleCCD::CCDData ccdData;
    bool ret = CCDKernel::TriangleCCD::CCDTest(static_cast<int>(vertices.size()), curP0.data(), curP1.data(),
      static_cast<int>(triangles.size()), triangles.data(),
      triPair.first, triPair.second, CCDKernel::TriangleCCD::CCDM_3RDPARTY, &ccdData);

    if (ret) {
      rd->colliingTrianglePairsTLS.local().push_back(ccdData);
    }
  },
    tbb::static_partitioner());

  if (keepPrevious == 0) {
    collidingTrianglePairs.clear();
    for (auto it = rd->colliingTrianglePairsTLS.begin(); it != rd->colliingTrianglePairsTLS.end(); ++it)
      collidingTrianglePairs.insert(collidingTrianglePairs.end(), it->begin(), it->end());
  }
  else {
    std::set<std::pair<int, int>> visited;
    for (const auto &pr : collidingTrianglePairs) {
      std::pair<int, int> newTriPair{ pr.triA, pr.triB };
      if (newTriPair.first > newTriPair.second) {
        std::swap(newTriPair.first, newTriPair.second);
      }

      visited.emplace(newTriPair);
    }

    for (auto it = rd->colliingTrianglePairsTLS.begin(); it != rd->colliingTrianglePairsTLS.end(); ++it) {
      for (const auto &pr : *it) {
        std::pair<int, int> newTriPair{ pr.triA, pr.triB };
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
  lastCDTime = dura(t1, t2);
}

void TriangleMeshSelfContactHandler::executeDCD()
{
  hclock::time_point t1 = hclock::now();

  const std::vector<std::pair<int, int>> *potentialColliingTrianglePairsPtr;

  selfContactDetection->execute(curP0.data(), nullptr);
  potentialColliingTrianglePairsPtr = &selfContactDetection->getPotentialCollidingTrianglePairs();

  for (auto it = rd->colliingTrianglePairsTLS.begin(); it != rd->colliingTrianglePairsTLS.end(); ++it)
    it->clear();

  tbb::parallel_for((size_t)0, potentialColliingTrianglePairsPtr->size(), [&](size_t tritriID) {
    const auto &triPair = potentialColliingTrianglePairsPtr->at(tritriID);
    Vec3i triA = triangles[triPair.first];
    Vec3i triB = triangles[triPair.second];

    bool ret = intersectTriTri(curP0.data() + triA[0] * 3, curP0.data() + triA[1] * 3, curP0.data() + triA[2] * 3,
      curP0.data() + triB[0] * 3, curP0.data() + triB[1] * 3, curP0.data() + triB[2] * 3);

    if (ret) {
      CCDKernel::TriangleCCD::CCDData ccdData;
      ccdData.triA = triPair.first;
      ccdData.triB = triPair.second;
      ccdData.ccdCase = CCDKernel::TriangleCCD::CCDC_COLLIDED;
      ccdData.t = 0;

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

void TriangleMeshSelfContactHandler::setExcludedVertices(const std::vector<int> &excludedVertices)
{
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

void TriangleMeshSelfContactHandler::computeSamplePosition(const ES::VXd &P, ES::VXd &SP) const
{
  tbb::parallel_for(
    0, (int)sampleInfoAndIDs.size(), [&](int si) {
      int tri = sampleInfoAndIDs[si].triangleID;

      ES::V3d vtx[3] = {
        P.segment<3>(triangles[tri][0] * 3),
        P.segment<3>(triangles[tri][1] * 3),
        P.segment<3>(triangles[tri][2] * 3)
      };

      // compute sample position
      const auto &sinfo = sampleInfoAndIDs[si];
      ES::V3d sampleP = vtx[0] * sinfo.w[0] + vtx[1] * sinfo.w[1] + vtx[2] * sinfo.w[2];
      SP.segment<3>(si * 3) = sampleP;
    },
    tbb::static_partitioner());
}

std::shared_ptr<PointTrianglePairCouplingEnergyWithCollision> TriangleMeshSelfContactHandler::buildContactEnergy(int checkingNeighboringContact, int changingTriangle)
{
  contactEnergyObjIDs.assign(contactedTrianglePairs.size(), std::array<int, 4>{ 0, 0, 0, 0 });

  contactEnergyObjectDOFOffsets[0] = 0;
  contactEnergyObjectDOFOffsets[1] = totaln3;

  return std::make_shared<PointTrianglePairCouplingEnergyWithCollision>((int)contactedTrianglePairs.size(), 1,
    contactEnergyObjIDs.data(),
    contactedTrianglePairs.data(),
    contactedTriangleIDs.data(),
    contactEnergyObjectDOFOffsets.data(),
    &surfaceMeshRef,
    &selfContactDetection->getTriMeshNeighbor(),
    &vertexID2SampleIDsLinear,
    &interpolationMatrix,
    &sampleWeights);
}