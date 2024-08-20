/*
author: Bohan Wang
copyright to USC
*/

#pragma once

#include "triMeshGeo.h"
#include "boundingVolumeTree.h"
#include "triMeshNeighbor.h"

#include <tbb/parallel_for.h>
#include <tbb/enumerable_thread_specific.h>

#include <vector>
#include <utility>

namespace pgo
{
namespace Contact
{
class TriangleMeshSelfContactDetection
{
public:
  TriangleMeshSelfContactDetection(Mesh::TriMeshRef mesh);

  void execute(const double *positions0, const double *positions1);
  const std::vector<std::pair<int, int>> &getPotentialCollidingTrianglePairs() const { return potentialCollidingTrianglePairs; }
  const Mesh::TriMeshNeighbor &getTriMeshNeighbor() const { return triMeshNeighbor; }

protected:
  static bool CCDIntersect(const Mesh::BoundingBox &AABB1Start, const Mesh::BoundingBox &AABB1End,
    const Mesh::BoundingBox &AABB2Start, const Mesh::BoundingBox &AABB2End);

  static bool DCDIntersect(const Mesh::BoundingBox &AABB1, const Mesh::BoundingBox &AABB2);

  Mesh::TriMeshRef triangleMeshRef;
  Mesh::TriMeshBVTree bvTree;
  Mesh::TriMeshNeighbor triMeshNeighbor;

  std::vector<Mesh::BoundingBox> lastAABBs, curAABBs;
  std::vector<int> triangles;
  std::vector<std::pair<int, int>> potentialCollidingTrianglePairs;

  struct BVTTNode
  {
    int nodeA, nodeB;
    BVTTNode(int na, int nb):
      nodeA(na), nodeB(nb) {}
  };

  using TrianglePairBuffer = std::vector<std::pair<int, int>, tbb::cache_aligned_allocator<std::pair<int, int>>>;
  using BVTTNodeBuffer = std::vector<BVTTNode, tbb::cache_aligned_allocator<BVTTNode>>;

  void addressSingleNode(const BVTTNode &node, int ccd);
  void assignTriangles(int nodeA, int nodeB, int ccd, TrianglePairBuffer &buf);
  void expandBothNodes(int nodeA, int nodeB, int ccd, BVTTNodeBuffer &buf);
  void expandFirstNodes(int nodeA, int nodeB, int ccd, BVTTNodeBuffer &buf);

  bool shareVertex(int triA, int triB);

  std::vector<BVTTNode> travseralFrontier0, travseralFrontier1;
  std::vector<BVTTNode> *lastFrontier, *currentFrontier;
  size_t parallel_threshold = 1000ull;

  tbb::enumerable_thread_specific<TrianglePairBuffer> trianglePairBufferTLS;
  tbb::enumerable_thread_specific<BVTTNodeBuffer> bvttNodeBufferTLS;
};

}  // namespace Contact
}  // namespace pgo
