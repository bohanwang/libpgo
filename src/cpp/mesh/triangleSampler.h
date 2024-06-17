#pragma once

#include "meshLinearAlgebra.h"
#include "triKey.h"

#include <cassert>
#include <cmath>

namespace pgo
{
namespace Mesh
{

//      /\         Left:
//     /__\        2D illustration of subdivision
//    /\  /\       2D triangle subdivision = 3:
//   /__\/__\      3D tets are the same
//  /\  /\  /\     //
// /__\/__\/__\    //
// In 2D, for subdivision k, the number of sub-triangles are: 1 + 3 + 5 + ... + (2k-1), where the first row has 1 sub-trianlge, second row 3, ...
//                                                           = k*k
//                           the number of sub-tri vtx are: 1 + 2 + 3 + ... + (k+1) = (k+2)(k+1)/2
// Each sub-tri vtx is surrounded by either 1, 3, or 6
// So assume the area of the triangle is 1, then each sub-tri has area 1/k^2,
// there will be 3 sub-tri vtx has 1 nbring sub-tri, 3(k-1) sub-tri vtx has 3 nbring sub-tri,
// and ((k+2)(k+1)/2 - 3k) sub-tri vtx has 6 nbring sub-tri
// If we sum up all area belonging to each sub-tri vtx:
// 3 * 1/(3kk) + 3(k-1) * 3/(3kk) + ((k+2)(k+1)/2 - 3k) * 6/(3kk) = 1, which is correct
class TriangleSampler
{
public:
  inline TriangleSampler(int numSubdivisionsOnEdge);

  // i, j: sample indices inside the triangle, i: [0, N], j: [0, N]
  // area: is the part of the triangle area belonging to this sample
  // baryWeight: barycentric weight of this sample w.r.t. the triangle vertices
  // pos: position of the sample
  using SamplePositionKernel = std::function<void(int i, int j, double area, const Vec3d &baryWeight, const Vec3d &pos)>;
  using SampleWeightKernel = std::function<void(int i, int j, double area, const Vec3d &baryWeight)>;

  // Kernel must be a SamplePositionKernel
  template<class Kernel>
  void visitSample(const Vec3d &triPosA, const Vec3d &triPosB, const Vec3d &triPosC, double triangleSurfaceArea, Kernel f) const;

  // Kernel must be a SampleIterationKernel
  template<class Kernel>
  void visitSample(double triangleSurfaceArea, Kernel f) const;

  inline int getNumEdgeSubdivisions() const { return N; }

  // get the sample index representing a triangle vertex (which is also a corner of the triangle area)
  const Vec2i &getCornerID(int i) const { return cornerIDs[i]; }

  // whether the sample index represents a triangle vertex
  inline bool isCornerIndex(int i, int j) const { return (i == 0 && j == 0) || (i == N && (j == 0 || j == N)); }

  // get an unique key to each sample (samples representing the same coordinate on the triangle mesh count as one)
  // values in sampleKey can be -1, which means the sample is on the edge of a triangle and it is interpolated by
  // less than three vertices
  inline UTriKey getSampleKey(const Vec3i &tri, int i, int j) const;

  // sampleKey is from getSampleKey()
  // values in sampleKey can be -1, in which case, the corresponding slot in vtxIndex is also invalidated.
  // e.g. N = 3, if sampleKey is (-1,-1, 43), then vtxIndex is (-1, -1, 10), barycentricWeight is (0.0, 0.0, 1.0)
  //             if sampleKey is (-1, 22, 41), then vtxIndex is (-1, 5, 10), barycentricWeight is (0.0, 2/3., 1/3.)
  inline void getBarycentricWeightFromKey(const UTriKey &sampleKey, Vec3i &vtxIndex, Vec3d &barycentricWeight) const;

  // get a feature key used to find the feature that this sample lies
  // a feature key can be:
  // (-1, -1, vtxID) representing a sample on a vtx, or
  // (-1, edgeVtxID0, edgeVtxID1) representing a sample on an edge, or
  // (triVtxID0, triVtxID1, triVtxID1) representing a sample inside a triangle
  inline UTriKey getFeatureKeyFromSampleKey(const UTriKey &sampleKey) const;

protected:
  int N = 1;
  double invN = 1.0;
  double edgeAreaWeight = 1.0;
  double cornerAreaWeight = (1.0 / 3.0);
  double interiorAreaWeight = 2;
  Vec2i cornerIDs[3];
};

inline TriangleSampler::TriangleSampler(int N):
  N(N), cornerIDs{ Vec2i(0, 0), Vec2i(N, 0), Vec2i(N, N) }
{
  invN = 1.0 / N;
  edgeAreaWeight = invN * invN;
  cornerAreaWeight = (1.0 / 3.0) * edgeAreaWeight;
  interiorAreaWeight = 2 * edgeAreaWeight;
}

template<class Kernel>
void TriangleSampler::visitSample(const Vec3d &triPosA, const Vec3d &triPosB, const Vec3d &triPosC, double triangleSurfaceArea, Kernel kernel) const
{
  static_assert(std::is_convertible<Kernel, SamplePositionKernel>::value, "TriangleSampler wrong kernel");
  const double cornerWeight = cornerAreaWeight * triangleSurfaceArea;
  const double edgeWeight = edgeAreaWeight * triangleSurfaceArea;
  const double interiorWeight = interiorAreaWeight * triangleSurfaceArea;

  const Vec3d start = triPosA;
  const Vec3d edge1 = triPosB - start;
  const Vec3d edge2 = triPosC - triPosB;

  const Vec3d vec1 = invN * edge1;
  const Vec3d vec2 = invN * edge2;

  Vec3d basei = start;

  for (int i = 0; i <= N; i++, basei += vec1) {
    Vec3d p = basei;                                                   // the sample pos, which is also the vtx of the sub-triangles
    for (int j = 0; j <= i; j++, p += vec2) {
      const Vec3d baryWeight(1 - i * invN, (i - j) * invN, j * invN);  // this wieght is also correct, no matter how the triangle deformes

      // compute the surface area associated with this sample point
      double area = 0.0;
      if (isCornerIndex(i, j)) {
        area = cornerWeight;
        //        numCornerCases++;
      }
      else if (j == 0 || i == N || i == j) {
        area = edgeWeight;
        //        numEdgeCases++;
      }
      else {
        area = interiorWeight;
        //        numInteriorCases++;
      }

      kernel(i, j, area, baryWeight, p);
    }  // end sample kernel
  }    // end i
}

template<class Kernel>
void TriangleSampler::visitSample(double triangleSurfaceArea, Kernel kernel) const
{
  static_assert(std::is_convertible<Kernel, SampleWeightKernel>::value, "TriangleSampler wrong kernel");
  const double cornerWeight = cornerAreaWeight * triangleSurfaceArea;
  const double edgeWeight = edgeAreaWeight * triangleSurfaceArea;
  const double interiorWeight = interiorAreaWeight * triangleSurfaceArea;

  for (int i = 0; i <= N; i++) {
    for (int j = 0; j <= i; j++) {
      const Vec3d baryWeight(1 - i * invN, (i - j) * invN, j * invN);  // this wieght is also correct, no matter how the triangle deformes

      // compute the surface area associated with this sample point
      double area = 0.0;
      if (isCornerIndex(i, j)) {
        area = cornerWeight;
      }
      else if (j == 0 || i == N || i == j) {
        area = edgeWeight;
      }
      else {
        area = interiorWeight;
      }

      kernel(i, j, area, baryWeight);
    }  // end sample kernel
  }    // end i
}

// the baryWeight is (1-i*invN, (i-j)*invN, j*invN)
// Remove the scale of invN, it is : (N-i, i-j, j), which are integers
// So we code triangle vtx ID with the bary weight by: vtxID * N + weight
inline UTriKey TriangleSampler::getSampleKey(const Vec3i &tri, int i, int j) const
{
  int a = (i < N ? (tri[0] * (N + 1) + (N - i)) : -1);
  int b = (i != j ? (tri[1] * (N + 1) + (i - j)) : -1);
  int c = (j > 0 ? (tri[2] * (N + 1) + j) : -1);
  return UTriKey(a, b, c);
}

inline void TriangleSampler::getBarycentricWeightFromKey(const UTriKey &sampleKey, Vec3i &vtxIndex, Vec3d &barycentricWeight) const
{
  Vec3d sampleWeight(0.0, 0.0, 0.0);
  Vec3i sampleVtxID(-1, -1, -1);
  for (int i = 0; i < 3; i++) {
    if (sampleKey[i] < 0)
      continue;
    sampleVtxID[i] = sampleKey[i] / (N + 1);
    sampleWeight[i] = invN * (sampleKey[i] % (N + 1));
  }

  assert(fabs(sampleWeight[0] + sampleWeight[1] + sampleWeight[2] - 1.0) < 1e-6);
  vtxIndex = sampleVtxID;
  barycentricWeight = sampleWeight;
}

inline UTriKey TriangleSampler::getFeatureKeyFromSampleKey(const UTriKey &sampleKey) const
{
  int fk[3];
  for (int i = 0; i < 3; i++) {
    if (sampleKey[i] < 0)
      fk[i] = -1;
    else
      fk[i] = sampleKey[i] / (N + 1);
  }
  return UTriKey(fk);
}

}  // namespace Mesh
}  // namespace pgo
