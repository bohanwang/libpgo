/*
copyright to Bohan Wang
*/

#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace pgo
{
namespace Contact
{
namespace CIPC
{

struct PTPair
{
  int p, t0, t1, t2;  // vertex indices
  double weight;      // area(point) * area(triangle)
};

struct EEPair
{
  int ea0, ea1, eb0, eb1;  // vertex indices
  double weight;           // length(edgeA) * length(edgeB)
};

// Phase 2 external pairs: "dynamic side fields written first".
// Index-space convention: dynamic-side fields use dyn-surface-global indices;
// obstacle-side fields use obstacle-local indices. No runtime orientation tag.

struct ExternalPTPair
{                                     // dyn vertex x obs triangle
  int32_t            obstacleSlot;
  int                dynVertex;       // dyn-surface-global index
  std::array<int, 3> obsTri;          // obstacle-local indices
  double             weight;          // area(point) * area(triangle)
};

struct ExternalTPPair
{                                     // dyn triangle x obs vertex
  int32_t            obstacleSlot;
  std::array<int, 3> dynTri;          // dyn-surface-global indices
  int                obsVertex;       // obstacle-local index
  double             weight;
};

struct ExternalEEPair
{                                     // dyn edge x obs edge
  int32_t            obstacleSlot;
  std::array<int, 2> dynEdge;         // dyn-surface-global indices (row order in unique_edges)
  std::array<int, 2> obsEdge;         // obstacle-local indices (row order in ObstacleSurface::uniqueEdges)
  double             weight;          // length(edgeA) * length(edgeB)
};

struct ExternalPairSet
{
  std::vector<ExternalPTPair> ptPairs;
  std::vector<ExternalTPPair> tpPairs;
  std::vector<ExternalEEPair> eePairs;

  void clear()
  {
    ptPairs.clear();
    tpPairs.clear();
    eePairs.clear();
  }

  std::size_t size() const
  {
    return ptPairs.size() + tpPairs.size() + eePairs.size();
  }
};

struct SelfPairSet
{
  std::vector<PTPair> ptPairs;
  std::vector<EEPair> eePairs;

  void clear()
  {
    ptPairs.clear();
    eePairs.clear();
  }

  std::size_t size() const
  {
    return ptPairs.size() + eePairs.size();
  }
};

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
