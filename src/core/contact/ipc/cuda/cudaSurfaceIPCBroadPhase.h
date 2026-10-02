/*
copyright to Bohan Wang
*/

#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>

namespace pgo::Contact::CIPC
{

// Device-only topology view. The caller retains the arrays for this object's
// lifetime. Vertex positions use interleaved xyz, triangles/edges row-major.
struct CudaSurfaceIPCTopologyView
{
  const int *triangles;
  const int *edges;
  const std::uint8_t *deformable;
  int numVertices, numTriangles, numEdges;
};

struct CudaSurfaceIPCCandidate
{
  int kind;  // 0: point-triangle; 1: edge-edge, first < second
  int first, second;
};

// GPU Morton BVHs with double outward-rounded AABBs. Candidate generation
// excludes incident primitives and fixed-fixed pairs. For swept queries the
// bounds include x and x + dx. Padding applies to both primitive bounds.
// Optional per-coordinate error intervals enlarge both endpoints for CCD.
class CudaSurfaceIPCBroadPhase
{
public:
  explicit CudaSurfaceIPCBroadPhase(CudaSurfaceIPCTopologyView topology);
  ~CudaSurfaceIPCBroadPhase();
  void detect(const double *devicePositions, const double *deviceMotion, double padding,
    const double *deviceCoordinateErrors = nullptr);
  const CudaSurfaceIPCCandidate *getDeviceCandidates() const;
  std::size_t getNumCandidates() const;

private:
  class Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace pgo::Contact::CIPC
