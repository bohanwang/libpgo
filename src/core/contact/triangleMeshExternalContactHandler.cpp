#include "triangleMeshExternalContactHandler.h"
#include "pointPenetrationEnergy.h"

#include "pgoLogging.h"
#include "geometryQuery.h"
#include "triangleSampler.h"
#include "basicAlgorithms.h"

#include <tbb/parallel_for.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/spin_mutex.h>

#include <queue>
#include <unordered_map>
#include <unordered_set>

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

TriangleMeshExternalContactHandler::TriangleMeshExternalContactHandler(const std::vector<ES::V3d> &V, const std::vector<ES::V3i> &T, int nDOFs,
  const std::vector<TriMeshRef> &esurf, int subdivideTriangle, const std::vector<int> *vertexEmbeddingIndices, const std::vector<double> *vertexEmbeddingWeights):
  vertices(V),
  triangles(T), surfaceMeshRef(vertices, triangles),
  totaln3(nDOFs)
{
  SPDLOG_LOGGER_INFO(Logging::lgr(), "Computing triangle area...");

  std::vector<double> triangleAreas(triangles.size(), 0);
  tbb::parallel_for(
    0, (int)triangles.size(), [&](int trii) {
      ES::V3i tri = triangles[trii];
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
      // for (int trii = 0; trii < (int)triangles.size(); trii++) {
      ES::V3i tri = triangles[trii];
      ES::V3d vtx[3] = {
        vertices[tri[0]],
        vertices[tri[1]],
        vertices[tri[2]],
      };

      triangleSampler->visitSample(vtx[0], vtx[1], vtx[2], triangleAreas[trii],
        [&](int i, int j, double area, const ES::V3d &baryWeight, const ES::V3d &pos) {
          SampleInfo sampleInfo;

          UTriKey key = triangleSampler->getSampleKey(tri, i, j);
          sampleInfo.id[0] = key[0];
          sampleInfo.id[1] = key[1];
          sampleInfo.id[2] = key[2];
          sampleInfo.idSorted = sortThree(sampleInfo.id);

          sampleInfo.pos = pos;
          sampleInfo.w = baryWeight;
          sampleInfo.sampleWeight = area;
          sampleInfo.coord = ES::V2i(i, j);
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
    // for (auto it = triangleSamples.begin(); it != triangleSamples.end(); ++it) {
    //   count += it->size();
    //   for (const auto &sample : *it) {
    //     sampleIDQueryTable.emplace(sample, 0);
    //   }
    // }

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

      // ES::V3d p = vertices[vtxID];
      // ES::V3d q = sample.pos;
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
  SPDLOG_LOGGER_INFO(Logging::lgr(), "  #vtx : {}", vertices.size());

  // if embedding weights and indices are given
  if (vertexEmbeddingIndices && vertexEmbeddingWeights) {
    // we need to first compute an interpolation matrix
    tbb::concurrent_vector<ES::TripletD> entries;

    // for (auto it = sampleIDQueryTable.begin(); it != sampleIDQueryTable.end(); ++it) {
    tbb::parallel_for(0, (int)sampleInfoAndIDs.size(), [&](int si) {
      const ES::V3d &p = sampleInfoAndIDs[si].pos;
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

      ES::V3i tri = triangles[triID];
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
  restP = ES::VXd::Zero(n3);
  for (int i = 0; i < static_cast<int>(vertices.size()); i++)
    restP.segment<3>(i * 3) = ES::V3d(vertices[i][0], vertices[i][1], vertices[i][2]);

  curP = restP;
  lastP = restP;

  samplen3 = (int)sampleInfoAndIDs.size() * 3;
  sampleRestP = ES::VXd::Zero(samplen3);

  computeSamplePosition(restP, sampleRestP);
  sampleCurP = sampleRestP;
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

  surfaceMeshRuntime = surfaceMeshRef;
  surfaceMeshNormals.buildPseudoNormals(surfaceMeshRuntime);

  for (const auto &surf : esurf) {
    externalSurfaces.emplace_back(surf);

    externalSurfaceBVTrees.emplace_back();
    externalSurfaceBVTrees.back().buildByInertiaPartition(externalSurfaces.back());

    externalSurfaceNormals.emplace_back();
    externalSurfaceNormals.back().buildPseudoNormals(externalSurfaces.back());
  }
}

void TriangleMeshExternalContactHandler::updateExternalSurface(int idx, TriMeshRef esurf)
{
  externalSurfaces[idx] = esurf;
  externalSurfaceBVTrees[idx].clear();
  externalSurfaceBVTrees[idx].buildByInertiaPartition(externalSurfaces[idx]);

  externalSurfaceNormals[idx] = TriMeshPseudoNormal();
  externalSurfaceNormals[idx].buildPseudoNormals(externalSurfaces[idx]);
}

void TriangleMeshExternalContactHandler::execute(const std::vector<ES::V3d> &p0)
{
  for (int i = 0; i < static_cast<int>(p0.size()); i++) {
    curP.segment<3>(i * 3) = ES::V3d(p0[i][0], p0[i][1], p0[i][2]);
    surfaceMeshRuntime.pos(i) = p0[i];
  }

  computeSamplePosition(curP, sampleCurP);

  execute();
}

void TriangleMeshExternalContactHandler::execute(const double *u0)
{
  curP.noalias() = restP + Eigen::Map<const ES::VXd>(u0, n3);

  for (int i = 0; i < surfaceMeshRuntime.numVertices(); i++) {
    surfaceMeshRuntime.pos(i) = asVec3d(curP.data() + i * 3);
  }
  // BoundingBox bb(surfaceMeshRuntime.positions());

  computeSamplePosition(curP, sampleCurP);

  //Mesh::TriMeshGeo z;
  //for (int i = 0; i < (int)sampleCurP.size() / 3; i++) {
  //  z.addPos(sampleCurP.segment<3>(i * 3));
  //}
  //z.save("k.obj");

  execute();
}

void TriangleMeshExternalContactHandler::execute()
{
  hclock::time_point t1 = hclock::now();

  surfaceMeshNormals.updateVertexPositions(surfaceMeshRuntime);

  struct ContactInfo
  {
    int sId;
    int objId;
    ES::V3d closestPt;
    ES::V3d tgtNormal;
  };

  tbb::enumerable_thread_specific<std::vector<ContactInfo, tbb::cache_aligned_allocator<ContactInfo>>> contactInfoTLS;

  // ===========================================================
  // check bone contact
  for (auto it = contactInfoTLS.begin(); it != contactInfoTLS.end(); ++it) {
    it->clear();
  }

  // for each vertices
  tbb::parallel_for(0, (int)sampleInfoAndIDs.size(), [&](int si) {
    // for (int vi = 0; vi < (int)vertices.size(); vi++) {
    // for (int si = 0; si < (int)sampleInfoAndIDs.size(); si++) {
    ES::V3d srcPos = asVec3d(sampleCurP.data() + si * 3);
    int triID = sampleInfoAndIDs[si].triangleID;
    if (std::binary_search(excludedTriangles.begin(), excludedTriangles.end(), triID))
      return;

    ES::V3d baryW = sampleInfoAndIDs[si].w;
    ES::V3d srcNormal = surfaceMeshNormals.vtxNormal(triangles[triID][0]) * baryW[0] +
      surfaceMeshNormals.vtxNormal(triangles[triID][1]) * baryW[1] +
      surfaceMeshNormals.vtxNormal(triangles[triID][2]) * baryW[2];

    // find the furthest bone that it penetrates
    TriMeshBVTree::ClosestTriangleQueryResult maxRet;
    maxRet.dist2 = 0;
    int closestObjectID = -1;
    ES::V3d closestTgtNormal;

    thread_local std::vector<std::tuple<double, double, int>> nodeStack;

    for (int oi = 0; oi < (int)externalSurfaces.size(); oi++) {
      const TriMeshGeo &mesh = externalSurfaces[oi];
      const TriMeshBVTree &meshBVTree = externalSurfaceBVTrees[oi];
      const TriMeshPseudoNormal &meshNormals = externalSurfaceNormals[oi];

      nodeStack.clear();
      auto ret = meshBVTree.closestTriangleQuery(mesh, srcPos, nodeStack);
      ES::V3d tgtNormal = meshNormals.getPseudoNormal(mesh.triangles().data(), ret.triID, ret.feature);

      if (tgtNormal.dot(srcNormal) > 0)
        continue;

      ES::V3d tgtPos = ret.closestPosition;
      ES::V3d diff = srcPos - tgtPos;

      if (diff.dot(tgtNormal) > 0)
        continue;

      if (ret.dist2 > maxRet.dist2) {
        maxRet = ret;
        closestObjectID = oi;
        closestTgtNormal = tgtNormal;
      }
    }

    if (closestObjectID >= 0) {
      ContactInfo info;
      info.objId = closestObjectID;
      info.closestPt = maxRet.closestPosition;
      info.tgtNormal = closestTgtNormal;
      info.sId = si;

      contactInfoTLS.local().push_back(info);
    }
  });

  int count = 0;
  for (auto it = contactInfoTLS.begin(); it != contactInfoTLS.end(); ++it) {
    count += (int)it->size();
  }

  SPDLOG_LOGGER_INFO(Logging::lgr(), "# external contacts: {}", count);

  constraintCoeffs.resize(count);
  constraintNormals.resize(count * 3);
  constraintTargetPositions.resize(count * 3);

  barycentricIdx.clear();
  barycentricIdx.resize(count);

  barycentricWeights.clear();
  barycentricWeights.resize(count);

  contactedTriangles.clear();
  contactedTriangles.resize(count);

  contactedSamples.clear();
  contactedSamples.resize(count);

  count = 0;
  for (auto it = contactInfoTLS.begin(); it != contactInfoTLS.end(); ++it) {
    for (const auto &info : *it) {
      int vi = 0;
      for (ES::SpMatD::InnerIterator it(interpolationMatrix, info.sId * 3); it; ++it) {
        barycentricIdx[count].emplace_back((int)it.col() / 3);
        barycentricWeights[count].emplace_back(it.value());
        vi++;
      }

      constraintCoeffs[count] = sampleWeights[info.sId];
      constraintNormals.segment<3>(count * 3) = info.tgtNormal;
      constraintTargetPositions.segment<3>(count * 3) = info.closestPt;

      contactedTriangles[count] = sampleTriangleIDs[info.sId][0];
      contactedSamples[count] = info.sId;

      count++;
    }
  }

  hclock::time_point t2 = hclock::now();

  lastCDTime = dura(t1, t2);

  SPDLOG_LOGGER_INFO(Logging::lgr(), "CD time: {}", lastCDTime);
}

std::shared_ptr<PointPenetrationEnergy> TriangleMeshExternalContactHandler::buildContactEnergy(int checkingContactRuntime)
{
  if (checkingContactRuntime) {
    return std::make_shared<PointPenetrationEnergy>((int)barycentricIdx.size(), (int)interpolationMatrix.cols(),
      constraintCoeffs.data(), constraintTargetPositions.data(), constraintNormals.data(),
      barycentricIdx, barycentricWeights, 1.0,
      1, 1, &externalSurfaces, &externalSurfaceBVTrees, &externalSurfaceNormals);
  }
  else {
    return std::make_shared<PointPenetrationEnergy>((int)barycentricIdx.size(), (int)interpolationMatrix.cols(),
      constraintCoeffs.data(), constraintTargetPositions.data(), constraintNormals.data(),
      barycentricIdx, barycentricWeights, 1.0,
      1, 1, nullptr, nullptr, nullptr);
  }
}

std::vector<ES::V3d> TriangleMeshExternalContactHandler::getSamplePoints(const double *u) const
{
  ES::VXd P = restP + Eigen::Map<const ES::VXd>(u, n3);
  ES::VXd SP = sampleRestP;
  computeSamplePosition(P, SP);

  std::vector<ES::V3d> pts;
  for (ES::IDX i = 0; i < SP.size() / 3; i++) {
    pts.push_back(ES::V3d(SP.data() + i * 3));
  }

  return pts;
}

std::vector<ES::V3d> TriangleMeshExternalContactHandler::getCollidedSamplePoints(const double *u) const
{
  ES::VXd P = restP + Eigen::Map<const ES::VXd>(u, n3);
  ES::VXd SP = sampleRestP;
  computeSamplePosition(P, SP);

  std::vector<ES::V3d> pts;

  for (const auto &sid : contactedSamples) {
    pts.push_back(ES::V3d(SP.data() + sid * 3));
  }

  return pts;
}

std::vector<int> TriangleMeshExternalContactHandler::getCollidedSampleAffectedVertices() const
{
  std::vector<int> vertexIndices;

  for (const auto &sid : contactedSamples) {
    for (ES::SpMatD::InnerIterator it(interpolationMatrix, sid * 3); it; ++it) {
      vertexIndices.push_back((int)it.col() / 3);
    }
  }

  sortAndDeduplicateWithErase(vertexIndices);

  return vertexIndices;
}

void TriangleMeshExternalContactHandler::setExcludedVertices(const std::vector<int> &excludedVertices)
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

void TriangleMeshExternalContactHandler::computeSamplePosition(const ES::VXd &P, ES::VXd &SP) const
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