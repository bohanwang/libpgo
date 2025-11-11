#pragma once

#include "CCDKernel.h"
#include "EigenSupport.h"
#include "triMeshGeo.h"
#include "boundingVolumeTree.h"
#include "triMeshPseudoNormal.h"
#include "triangleSampler.h"

#include <unordered_map>
#include <memory>

namespace pgo
{
namespace Contact
{
class PointPenetrationEnergy;

class TriangleMeshExternalContactHandler
{
public:
  TriangleMeshExternalContactHandler(const std::vector<EigenSupport::V3d> &vertices, const std::vector<EigenSupport::V3i> &triangles, int nDOFs,
    const std::vector<Mesh::TriMeshRef> &externalSurfaces, int subdivideTriangle = 1,
    const std::vector<int> *vertexEmbeddingIndices = nullptr, const std::vector<double> *vertexEmbeddingWeights = nullptr);
  virtual ~TriangleMeshExternalContactHandler() {}

  void updateExternalSurface(int idx, Mesh::TriMeshRef esurf);

  virtual void execute(const std::vector<EigenSupport::V3d> &p0);
  virtual void execute(const double *u0);

  std::shared_ptr<PointPenetrationEnergy> buildContactEnergy(int checkingContactRuntime = 0);

  void setExcludedVertices(const std::vector<int> &excludedVertices);
  void setExcludedTriangles(const std::vector<int> &excludedTriangles_) { excludedTriangles = excludedTriangles_; }

  // get embedding matrix
  const EigenSupport::SpMatD &getSampleEmbeddingMatrix() const { return interpolationMatrix; }
  const std::vector<double> &getSampleWeights() const { return sampleWeights; }

  int getNumCollidingSamples() const { return (int)constraintCoeffs.size(); }

  std::vector<EigenSupport::V3d> getSamplePoints(const double *u) const;
  std::vector<EigenSupport::V3d> getCollidedSamplePoints(const double *u) const;
  std::vector<int> getCollidedSampleAffectedVertices() const;

  const EigenSupport::VXd &getLastCollisionPoints() const { return constraintTargetPositions; }
  const EigenSupport::VXd &getLastCollisionNormals() const { return constraintNormals; }

  double getLastCDTime() const { return lastCDTime; }

protected:
  void execute();

  void computeSamplePosition(const EigenSupport::VXd &P, EigenSupport::VXd &SP) const;

protected:
  typedef std::array<int, 3> SampleID;

  struct SampleInfo
  {
    SampleID id, idSorted;
    EigenSupport::V3d pos;
    EigenSupport::V3d w;
    EigenSupport::V2i coord;

    double sampleWeight;
    int triangleID;
  };

  struct SampleInfoEqual
  {
    bool operator()(const SampleInfo &s1, const SampleInfo &s2) const
    {
      return std::memcmp(s1.idSorted.data(), s2.idSorted.data(), sizeof(SampleID)) == 0;
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
      int ret = std::memcmp(e1.data(), e2.data(), sizeof(int) * e1.size());
      return ret < 0;
    }
  };

protected:
  std::vector<EigenSupport::V3d> vertices;
  std::vector<EigenSupport::V3i> triangles;
  Mesh::TriMeshRef surfaceMeshRef;
  Mesh::TriMeshGeo surfaceMeshRuntime;
  Mesh::TriMeshPseudoNormal surfaceMeshNormals;
  std::vector<int> excludedTriangles;

  std::vector<Mesh::TriMeshGeo> externalSurfaces;
  std::vector<Mesh::TriMeshBVTree> externalSurfaceBVTrees;
  std::vector<Mesh::TriMeshPseudoNormal> externalSurfaceNormals;

  int n3;
  EigenSupport::VXd restP;
  EigenSupport::VXd curP, lastP;

  int samplen3;
  EigenSupport::VXd sampleRestP;
  EigenSupport::VXd sampleCurP, sampleLastP;

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

  EigenSupport::VXd constraintCoeffs;
  EigenSupport::VXd constraintTargetPositions;
  EigenSupport::VXd constraintNormals;
  std::vector<int> contactedTriangles, contactedSamples;
  std::vector<std::vector<int>> barycentricIdx;
  std::vector<std::vector<double>> barycentricWeights;

  double lastCDTime = 0.0;
};
}  // namespace Contact

}  // namespace VegaFEM
