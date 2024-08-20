/*
author: Bohan Wang
copyright to USC
*/

#pragma once

#include "CCDKernel.h"
#include "EigenSupport.h"
#include "triMeshGeo.h"
#include "triangleSampler.h"

#include <unordered_map>
#include <memory>

namespace pgo
{
namespace Contact
{
class TriangleMeshSelfContactDetection;
class PointTrianglePairCouplingEnergyWithCollision;
class TriangleMeshSelfContactHandlerRuntimeData;

class TriangleMeshSelfContactHandler
{
public:
  TriangleMeshSelfContactHandler(const std::vector<EigenSupport::V3d> &vertices, const std::vector<EigenSupport::V3i> &triangles, int nDOFs,
    int subdivideTriangle = 1, const std::vector<int> *vertexEmbeddingIndices = nullptr,
    const std::vector<double> *vertexEmbeddingWeights = nullptr);
  virtual ~TriangleMeshSelfContactHandler() {}

  // dcd
  virtual void execute(const std::vector<EigenSupport::V3d> &p0);
  virtual void execute(const double *u0);

  // ccd
  virtual void execute(const std::vector<EigenSupport::V3d> &p0, const std::vector<EigenSupport::V3d> &p1);
  virtual void execute(const double *u0, const double *u1);

  // dcd handling
  virtual void handleContactDCD(double distThreshold = 0.0, int maxSearchingNumTriangles = 100);
  std::shared_ptr<PointTrianglePairCouplingEnergyWithCollision> buildContactEnergy(int checkingNeighboringContact = 0, int changingTriangle = 0);

  void setExcludedVertices(const std::vector<int> &excludedVertices);
  void setExcludedTriangles(const std::vector<int> &excludedTriangles_) { excludedTriangles = excludedTriangles_; }
  void setKeepPrevious(int val) { keepPrevious = val; }

  // dump colliding triangle pair geometry
  std::vector<std::array<EigenSupport::V3d, 3>> getCollidingTriangles(const double *u) const;
  // dump colliding point-triangle pair sample ID
  const std::vector<std::array<int, 4>> &getCollidingVtxTriPairs() const { return contactedTrianglePairs; }
  // dump colliding triangle pair info
  const std::vector<CCDKernel::TriangleCCD::CCDData> &getCollidingTrianglePair() const { return collidingTrianglePairs; }
  // get embedding matrix
  const EigenSupport::SpMatD &getSampleEmbeddingMatrix() const { return interpolationMatrix; }
  const std::vector<double> &getSampleWeights() const { return sampleWeights; }

  std::vector<EigenSupport::V3d> getSamplePoints(const double *u) const;
  std::vector<EigenSupport::V3d> getCollidedSamplePoints(const double *u) const;
  std::vector<int> getCollidedSampleAffectedVertices() const;

  double getLastCDTime() const { return lastCDTime; }

protected:
  void executeCCD();
  void executeDCD();

  void computeSamplePosition(const EigenSupport::VXd &P, EigenSupport::VXd &SP) const;

protected:
  typedef std::array<int, 3> SampleID;

  struct SampleInfo
  {
    SampleID id, idSorted;
    EigenSupport::V3d pos;
    EigenSupport::V3d w;
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

protected:
  std::vector<EigenSupport::V3d> vertices;
  std::vector<EigenSupport::V3i> triangles;
  Mesh::TriMeshRef surfaceMeshRef;

  int n3;
  EigenSupport::VXd restP;
  EigenSupport::VXd curP0, curP1, lastP;

  int samplen3;
  EigenSupport::VXd sampleRestP;
  EigenSupport::VXd sampleCurP0, sampleCurP1, sampleLastP;

  int totaln3;

  EigenSupport::SpMatD interpolationMatrix;
  std::shared_ptr<Mesh::TriangleSampler> triangleSampler;
  std::vector<std::vector<SampleInfo>> triangleSamples;
  std::unordered_map<SampleInfo, int, SampleInfoHash, SampleInfoEqual> sampleIDQueryTable;
  std::unordered_map<int, int> vertexID2SampleIDs;
  std::vector<int> vertexID2SampleIDsLinear;
  std::vector<SampleInfo> sampleInfoAndIDs;
  std::vector<std::vector<int>> sampleTriangleIDs;
  std::vector<double> sampleWeights;
  std::vector<int> excludedTriangles;

  std::shared_ptr<TriangleMeshSelfContactHandlerRuntimeData> rd;

  // new CCD BVH using BFS and tbb multi-threading,
  // generally 2x fast, the more node it has, the faster it is
  std::shared_ptr<TriangleMeshSelfContactDetection> selfContactDetection;

  using CCDD = CCDKernel::TriangleCCD::CCDData;
  std::vector<CCDD> collidingTrianglePairs;
  std::vector<std::array<int, 4>> contactedTrianglePairs;
  std::vector<int> contactedTriangleIDs;
  std::vector<std::array<int, 4>> contactEnergyObjIDs;
  std::array<int, 2> contactEnergyObjectDOFOffsets;
  double lastCDTime = 0.0;
  int keepPrevious = 0;
};
}  // namespace Contact

}  // namespace pgo
