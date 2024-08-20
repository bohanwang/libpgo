#include "triMeshSampler.h"
#include "triangleSampler.h"
#include "geometryQuery.h"

#include <tbb/parallel_for.h>
#include <tbb/concurrent_unordered_map.h>

#include <cstring>
#include <cstdlib>
#include <iostream>

namespace pgo
{
namespace Mesh
{

namespace triMeshSamplerInternal
{
using SampleID = std::array<int, 3>;

struct SampleInfo
{
  SampleID id, idSorted;
  Vec3d pos;
  Vec3d w;
  Vec2i coord;

  double sampleWeight;
  int triangleID;
};

struct SampleInfoEqual
{
  bool operator()(const SampleInfo &s1, const SampleInfo &s2) const
  {
    return memcmp(s1.idSorted.data(), s2.idSorted.data(), sizeof(SampleID)) == 0;
  }
};

struct SampleInfoHash
{
  void hc(int v, std::size_t &seed) const
  {
    seed ^= std::hash<int>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  }

  std::size_t operator()(const SampleInfo &s) const
  {
    std::size_t seed = std::hash<int>()(s.idSorted[0]);
    hc(s.idSorted[1], seed);
    hc(s.idSorted[2], seed);

    return seed;
  }
};

struct ArrayCmp4Less
{
  bool operator()(const std::array<int, 4> &e1, const std::array<int, 4> &e2) const
  {
    int ret = memcmp(e1.data(), e2.data(), sizeof(int) * e1.size());
    return ret < 0;
  }
};

inline SampleID sortThree(const SampleID &idIn)
{
  SampleID id = idIn;

  if (id[0] > id[1])
    std::swap(id[0], id[1]);

  if (id[0] > id[2])
    std::swap(id[0], id[2]);

  if (id[1] > id[2])
    std::swap(id[1], id[2]);

  return id;
};
}  // namespace triMeshSamplerInternal

int sampleTriangle(const TriMeshGeo &mesh, int subdivideTriangle,
  std::vector<int> &sampleTriangleIDs, std::vector<double> &sampleWeights, std::vector<std::array<double, 3>> &sampleBarycentricWeights)
{
  using namespace triMeshSamplerInternal;

  int ntri = mesh.numTriangles();
  std::vector<double> triangleAreas(ntri, 0);
  tbb::parallel_for(
    0, ntri, [&](int trii) {
      Vec3d vtx[3] = {
        mesh.pos(trii, 0),
        mesh.pos(trii, 1),
        mesh.pos(trii, 2),
      };

      triangleAreas[trii] = getTriangleArea(vtx[0], vtx[1], vtx[2]);
    },
    tbb::static_partitioner());

  std::cout << "Sampling surface mesh..." << std::endl;

  TriangleSampler triangleSampler(std::max(subdivideTriangle, 1));
  std::vector<std::vector<SampleInfo>> triangleSamples(ntri, std::vector<SampleInfo>());

  // Sample surface
  tbb::parallel_for(
    0, ntri, [&](int trii) {
      triangleSamples[trii].reserve(100);

      Vec3d vtx[3] = {
        mesh.pos(trii, 0),
        mesh.pos(trii, 1),
        mesh.pos(trii, 2),
      };

      triangleSampler.visitSample(vtx[0], vtx[1], vtx[2], triangleAreas[trii],
        [&](int i, int j, double area, const Vec3d &baryWeight, const Vec3d &pos) {
          SampleInfo sampleInfo;

          UTriKey key = triangleSampler.getSampleKey(mesh.tri(trii), i, j);
          sampleInfo.id[0] = key[0];
          sampleInfo.id[1] = key[1];
          sampleInfo.id[2] = key[2];
          sampleInfo.idSorted = sortThree(sampleInfo.id);

          sampleInfo.pos = pos;
          sampleInfo.w = baryWeight;
          sampleInfo.sampleWeight = area;
          sampleInfo.coord = Vec2i(i, j);
          sampleInfo.triangleID = trii;
          triangleSamples[trii].emplace_back(sampleInfo);
        });
    },
    tbb::static_partitioner());

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
        if (!((int64_t)sample.id[0] + (int64_t)sample.id[1] + (int64_t)sample.id[2] == (int64_t)maxID - 2)) {
          std::cerr << "Impossible sample ID" << std::endl;
          return 1;
        }

        sampleIDQueryTableCC.emplace(sample, maxID / 2);
      }
    }
  }

  std::cout << "Sampling done." << std::endl;
  std::cout << "#vtx: " << mesh.numVertices() << "; #tri: " << mesh.numTriangles() << "; #samples: " << sampleIDQueryTableCC.size() << std::endl;

  if (subdivideTriangle <= 1) {
    // double maxDist = 0;
    for (auto it = sampleIDQueryTableCC.begin(); it != sampleIDQueryTableCC.end(); ++it) {
      const auto &sample = it->first;
      int sampleID = it->second;

      if (triangleSampler.isCornerIndex(sample.coord[0], sample.coord[1])) {
        UTriKey key = triangleSampler.getFeatureKeyFromSampleKey(UTriKey(sample.id[0], sample.id[1], sample.id[2]));
        int vtxID = std::max({ key[0], key[1], key[2] });
        if (vtxID != sampleID) {
          std::cerr << "Sample ID != vtx ID" << std::endl;
          return 1;
        }
      }
    }
  }

  sampleWeights.assign(sampleIDQueryTableCC.size(), 0.0);
  for (size_t trii = 0; trii < triangleSamples.size(); trii++) {
    for (const auto &sample : triangleSamples[trii]) {
      auto it = sampleIDQueryTableCC.find(sample);
      if (it == sampleIDQueryTableCC.end()) {
        std::cerr << "Impossible error happens." << std::endl;
        return 1;
      }

      sampleWeights[it->second] += sample.sampleWeight;
    }
  }

  sampleTriangleIDs.assign(sampleIDQueryTableCC.size(), 0);
  sampleBarycentricWeights.assign(sampleIDQueryTableCC.size(), std::array<double, 3>());
  for (const auto &pr : sampleIDQueryTableCC) {
    sampleTriangleIDs[pr.second] = pr.first.triangleID;

    sampleBarycentricWeights[pr.second][0] = pr.first.w[0];
    sampleBarycentricWeights[pr.second][1] = pr.first.w[1];
    sampleBarycentricWeights[pr.second][2] = pr.first.w[2];
  }

  return 0;
}
}  // namespace Mesh
}  // namespace pgo