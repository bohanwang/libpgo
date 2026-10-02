/*
copyright to Bohan Wang
*/

#include "cudaEmbeddedSurfaceIPCPotentialEnergy.h"
#include "cudaIPCGeometry.h"
#include "cudaIPCContinuousCollision.h"
#include "cudaDeviceBuffer.h"
#include "cudaSurfaceIPCBroadPhase.h"
#include "ipc/topology/surfaceIPCTopology.h"

#include <libuipc/distance/distance_flagged.h>
#include <libuipc/distance/edge_edge_mollifier.h>
#include <libuipc/distance/ccd.h>
#include <libuipc/math/make_spd.h>

#include <cuda_runtime.h>
#include <math_constants.h>
#include <algorithm>
#include <array>
#include <climits>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace pgo::Contact::CIPC
{
namespace
{
using CudaUtilities::DeviceBuffer;
using CudaUtilities::checkCuda;
namespace gpuDistance = uipc::backend::cuda::distance;
using Vec3 = Eigen::Matrix<double, 3, 1>;
using Vec12 = Eigen::Matrix<double, 12, 1>;
using Mat12 = Eigen::Matrix<double, 12, 12>;

constexpr int blockSize = 128;
int blocks(int n)
{
  return n / blockSize + (n % blockSize != 0);
}

struct Pair
{
  int vertices[4];
  int kind;  // 0: point-triangle, 1: edge-edge
  double weight;
};

struct TopologyView
{
  const int *triangles;
  const int *edges;
  const uint8_t *deformable;
  const double *vertexArea;
  const double *triangleArea;
  const double *edgeLength;
  int numVertices, numTriangles, numEdges;
};

struct MapView
{
  const int *offsets;
  const int *columns;
  const double *values;
  int rows;
};

__device__ Vec3 position(const double *x, int v)
{
  return Eigen::Map<const Vec3>(x + 3 * v);
}

__global__ void mapKernel(MapView map, const double *x, const double *rest, double *out)
{
  const int row = blockIdx.x * blockDim.x + threadIdx.x;
  if (row >= map.rows)
    return;
  double value = rest ? rest[row] : 0.0;
  for (int k = map.offsets[row]; k < map.offsets[row + 1]; ++k)
    value += map.values[k] * x[map.columns[k]];
  out[row] = value;
}

// Bound every arithmetic operation used to map the chord and independently
// remap its endpoint. The absolute sums retain rest/displacement cancellation
// scale; an error estimate based only on the small world position is unsafe.
__global__ void mapCcdKernel(MapView map, const double *x, const double *dx,
  const double *rest, double *positions, double *motion, double *error)
{
  const int row = blockIdx.x * blockDim.x + threadIdx.x;
  if (row >= map.rows)
    return;
  double p = rest[row], v = 0.0;
  double positionSum = fabs(p), motionSum = 0.0, subtractionSum = 0.0;
  const int terms = map.offsets[row + 1] - map.offsets[row];
  for (int k = map.offsets[row]; k < map.offsets[row + 1]; ++k) {
    const double weight = map.values[k], value = x[map.columns[k]], step = dx[map.columns[k]];
    p += weight * value;
    v += weight * step;
    positionSum = __dadd_ru(positionSum, __dmul_ru(fabs(weight), fabs(value)));
    motionSum = __dadd_ru(motionSum, __dmul_ru(fabs(weight), fabs(step)));
    subtractionSum = __dadd_ru(subtractionSum,
      __dmul_ru(fabs(weight), __dadd_ru(fabs(value), fabs(step))));
  }
  // gamma_n = n*u/(1-n*u). Counting a multiply and add separately also
  // bounds fused evaluation. The final term covers represented endpoint
  // subtraction; the other three cover base, motion and endpoint remapping.
  const double unit = 0.5 * std::numeric_limits<double>::epsilon();
  const double nu = (2.0 * terms + 1.0) * unit;
  const double gamma = __ddiv_ru(nu, __dsub_rd(1.0, nu));
  const double sums = __dadd_ru(__dadd_ru(positionSum, motionSum), __dadd_ru(positionSum, motionSum));
  error[row] = terms ? __dadd_ru(__dmul_ru(gamma, sums), __dmul_ru(unit, subtractionSum)) : 0.0;
  positions[row] = p;
  motion[row] = v;
}

__global__ void transposeMapKernel(MapView map, const double *x, double *out)
{
  const int row = blockIdx.x * blockDim.x + threadIdx.x;
  if (row >= map.rows)
    return;
  for (int k = map.offsets[row]; k < map.offsets[row + 1]; ++k)
    atomicAdd(out + map.columns[k], map.values[k] * x[row]);
}

__device__ double squaredDistance(const Pair &pair, const double *x)
{
  const Vec3 a = position(x, pair.vertices[0]);
  const Vec3 b = position(x, pair.vertices[1]);
  const Vec3 c = position(x, pair.vertices[2]);
  const Vec3 d = position(x, pair.vertices[3]);
  double result;
  if (pair.kind == 0) {
    const auto flag = gpuDistance::point_triangle_distance_flag(a, b, c, d);
    gpuDistance::point_triangle_distance2(flag, a, b, c, d, result);
  }
  else {
    const auto flag = gpuDistance::edge_edge_distance_flag(a, b, c, d);
    gpuDistance::edge_edge_distance2(flag, a, b, c, d, result);
  }
  return result;
}

// The BVH has already applied topology and ownership filters. Compact only
// active barrier pairs; CCD retains every swept AABB candidate.
__global__ void filterPairsKernel(TopologyView topology, const CudaSurfaceIPCCandidate *candidates,
  int total, const double *x, bool swept, double dhat, Pair *pairs, int *count, int *invalid)
{
  const int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index >= total)
    return;
  const auto candidate = candidates[index];
  Pair pair;
  pair.kind = candidate.kind;
  if (pair.kind == 0) {
    pair.vertices[0] = candidate.first;
    for (int i = 0; i < 3; ++i)
      pair.vertices[i + 1] = topology.triangles[3 * candidate.second + i];
    pair.weight = topology.vertexArea[candidate.first] * topology.triangleArea[candidate.second];
  }
  else {
    for (int i = 0; i < 2; ++i) {
      pair.vertices[i] = topology.edges[2 * candidate.first + i];
      pair.vertices[i + 2] = topology.edges[2 * candidate.second + i];
    }
    pair.weight = topology.edgeLength[candidate.first] * topology.edgeLength[candidate.second];
  }
  if (!swept) {
    const double d2 = squaredDistance(pair, x);
    if (!isfinite(d2) || d2 <= 0.0) {
      atomicExch(invalid, 1);
      return;
    }
    if (d2 >= dhat * dhat)
      return;
  }
  const int slot = swept ? index : atomicAdd(count, 1);
  pairs[slot] = pair;
}

// These formulas intentionally match ipc/geometry/ipcBarrier.cpp, rather than
// libuipc's dt-scaled log-squared contact law.
__device__ void barrier(double s, double shat, double &value, double &first, double &second)
{
  const double r = s / shat;
  const double rm1 = r - 1.0;
  const double invr = 1.0 / r;
  const double lr = log(r);
  value = -rm1 * rm1 * lr;
  first = -(rm1 * (2.0 * lr + rm1 * invr)) / shat;
  second = -(2.0 * lr + 4.0 * rm1 * invr - rm1 * rm1 * invr * invr) / (shat * shat);
}

__global__ void mapDifferenceKernel(MapView map, const double *base, const double *trial,
  const double *positions, double *motion, double *trialPositions)
{
  const int row = blockIdx.x * blockDim.x + threadIdx.x;
  if (row >= map.rows)
    return;
  double difference = 0.0;
  for (int k = map.offsets[row]; k < map.offsets[row + 1]; ++k)
    difference += map.values[k] * (trial[map.columns[k]] - base[map.columns[k]]);
  motion[row] = difference;
  trialPositions[row] = positions[row] + difference;
}

__device__ double squaredNormChange(const Vec3 &base, const Vec3 &change)
{
  return (2.0 * base + change).dot(change);
}

// Squared primitive distances are polynomial quotients. Taking their changes
// before division preserves work from sub-ULP world-coordinate displacements.
__device__ double distanceChange(const Pair &pair, const double *base, const double *trial,
  const double *motion, double &s0, double &s1)
{
  Vec3 p[4], dp[4], q[4];
  bool stationary = true;
  for (int i = 0; i < 4; ++i) {
    p[i] = position(base, pair.vertices[i]);
    dp[i] = position(motion, pair.vertices[i]);
    q[i] = position(trial, pair.vertices[i]);
    stationary = stationary && dp[i].isZero(0.0);
  }
  if (stationary) {
    // Evaluate identical endpoints only once, including at feature boundaries.
    // Different compiler inlining/FMA choices must not manufacture motion.
    s0 = s1 = squaredDistance(pair, base);
    return 0.0;
  }
  const auto f0 = pair.kind == 0 ? gpuDistance::point_triangle_distance_flag(p[0], p[1], p[2], p[3]) :
                                   gpuDistance::edge_edge_distance_flag(p[0], p[1], p[2], p[3]);
  const auto f1 = pair.kind == 0 ? gpuDistance::point_triangle_distance_flag(q[0], q[1], q[2], q[3]) :
                                   gpuDistance::edge_edge_distance_flag(q[0], q[1], q[2], q[3]);
  if (f0 != f1) {
    s0 = squaredDistance(pair, base);
    s1 = squaredDistance(pair, trial);
    return s1 - s0;
  }
  uipc::Vector4i offsets;
  const int dimension = pair.kind == 0 ? gpuDistance::degenerate_point_triangle(f0, offsets) :
                                         gpuDistance::degenerate_edge_edge(f0, offsets);
  double change;
  if (dimension == 2) {
    const Vec3 w = p[offsets[0]] - p[offsets[1]], dw = dp[offsets[0]] - dp[offsets[1]];
    s0 = w.squaredNorm();
    change = squaredNormChange(w, dw);
  }
  else if (dimension == 3) {
    const Vec3 u = p[offsets[2]] - p[offsets[1]], du = dp[offsets[2]] - dp[offsets[1]];
    const Vec3 w = p[offsets[0]] - p[offsets[1]], dw = dp[offsets[0]] - dp[offsets[1]];
    const Vec3 n = u.cross(w), dn = du.cross(w) + u.cross(dw) + du.cross(dw);
    const double denominator = u.squaredNorm(), deltaDenominator = squaredNormChange(u, du);
    s0 = n.squaredNorm() / denominator;
    change = (squaredNormChange(n, dn) - s0 * deltaDenominator) / (denominator + deltaDenominator);
  }
  else {
    const int origin = pair.kind == 0 ? 1 : 0;
    const Vec3 u = p[pair.kind == 0 ? 2 : 1] - p[origin];
    const Vec3 du = dp[pair.kind == 0 ? 2 : 1] - dp[origin];
    const Vec3 v = p[3] - p[pair.kind == 0 ? 1 : 2];
    const Vec3 dv = dp[3] - dp[pair.kind == 0 ? 1 : 2];
    const Vec3 w = p[0] - p[pair.kind == 0 ? 1 : 2];
    const Vec3 dw = dp[0] - dp[pair.kind == 0 ? 1 : 2];
    const Vec3 n = u.cross(v), dn = du.cross(v) + u.cross(dv) + du.cross(dv);
    const double dot = w.dot(n), deltaDot = dw.dot(n) + w.dot(dn) + dw.dot(dn);
    const double denominator = n.squaredNorm(), deltaDenominator = squaredNormChange(n, dn);
    s0 = dot * dot / denominator;
    change = ((2.0 * dot + deltaDot) * deltaDot - s0 * deltaDenominator) / (denominator + deltaDenominator);
  }
  s1 = s0 + change;
  return change;
}

__device__ double barrierValue(double squaredDistance, double dhat2)
{
  if (squaredDistance >= dhat2)
    return 0.0;
  const double r = squaredDistance / dhat2;
  return -(r - 1.0) * (r - 1.0) * log(r);
}

__global__ void energyChangeKernel(const Pair *pairs, int count, const double *base, const double *trial, const double *motion,
  double dhat2, double kappa, double epsEE, double *change, int *invalid)
{
  const int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index >= count)
    return;
  const Pair pair = pairs[index];
  // Feature selection must be independent at the two endpoints. Swept AABB
  // candidates contain every pair active at either endpoint, including pairs
  // that enter or leave barrier support during a line search.
  double s0, s1;
  const double ds = distanceChange(pair, base, trial, motion, s0, s1);
  if (!(s0 > 0.0) || !(s1 > 0.0) || !isfinite(s0) || !isfinite(s1)) {
    atomicExch(invalid, 1);
    return;
  }
  if (s0 >= dhat2 && s1 >= dhat2)
    return;
  const double b0 = barrierValue(s0, dhat2), b1 = barrierValue(s1, dhat2);
  double delta = b1 - b0;
  if (s0 < dhat2 && s1 < dhat2) {
    // Exact algebra for -(r-1)^2 log(r), arranged to preserve small changes.
    // Work with the squared-distance increment before normalizing; subtracting
    // rounded normalized ratios can erase a representable distance change.
    const double t0 = s0 / dhat2 - 1.0;
    const double dr = ds / dhat2;
    const double t1 = t0 + dr;
    delta = -(2.0 * t0 + dr) * dr * log(s0 / dhat2) - t1 * t1 * log1p(ds / s0);
  }
  double m0 = 1.0, m1 = 1.0, dm = 0.0;
  if (pair.kind == 1 && epsEE > 0.0) {
    gpuDistance::edge_edge_mollifier(position(base, pair.vertices[0]), position(base, pair.vertices[1]),
      position(base, pair.vertices[2]), position(base, pair.vertices[3]), epsEE, m0);
    gpuDistance::edge_edge_mollifier(position(trial, pair.vertices[0]), position(trial, pair.vertices[1]),
      position(trial, pair.vertices[2]), position(trial, pair.vertices[3]), epsEE, m1);
    const Vec3 u = position(base, pair.vertices[1]) - position(base, pair.vertices[0]);
    const Vec3 v = position(base, pair.vertices[3]) - position(base, pair.vertices[2]);
    const Vec3 du = position(motion, pair.vertices[1]) - position(motion, pair.vertices[0]);
    const Vec3 dv = position(motion, pair.vertices[3]) - position(motion, pair.vertices[2]);
    const Vec3 n = u.cross(v), dn = du.cross(v) + u.cross(dv) + du.cross(dv);
    const double cross0 = n.squaredNorm(), deltaCross = squaredNormChange(n, dn);
    if (cross0 < epsEE && cross0 + deltaCross < epsEE)
      dm = (deltaCross / epsEE) * (2.0 - (2.0 * cross0 + deltaCross) / epsEE);
    else
      dm = m1 - m0;
  }
  const double weight = kappa * pair.weight;
  const double value = weight * (m0 * delta + dm * b1);
  const double scale = weight * (fabs(m0 * b0) + fabs(m1 * b1));
  if (!isfinite(value) || !isfinite(scale)) {
    atomicExch(invalid, 1);
    return;
  }
  atomicAdd(change, value);
  atomicAdd(change + 1, scale);
}

template<int Mode>
__global__ void evaluateKernel(const Pair *pairs, int count, const double *x,
  double dhat2, double kappa, double epsEE, bool projectPSD,
  double *energy, double *gradient, double *hessians)
{
  const int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index >= count)
    return;
  const Pair pair = pairs[index];
  const Vec3 a = position(x, pair.vertices[0]);
  const Vec3 b = position(x, pair.vertices[1]);
  const Vec3 c = position(x, pair.vertices[2]);
  const Vec3 d = position(x, pair.vertices[3]);
  double d2;
  Vec12 gd;
  Mat12 hd;
  if (pair.kind == 0) {
    const auto flag = gpuDistance::point_triangle_distance_flag(a, b, c, d);
    gpuDistance::point_triangle_distance2(flag, a, b, c, d, d2);
    if constexpr (Mode > 0)
      gpuDistance::point_triangle_distance2_gradient(flag, a, b, c, d, gd);
    if constexpr (Mode == 2)
      gpuDistance::point_triangle_distance2_hessian(flag, a, b, c, d, hd);
  }
  else {
    const auto flag = gpuDistance::edge_edge_distance_flag(a, b, c, d);
    if (flag.sum() == 4) {
      CudaIPCGeometry::interiorEdgeEdge(a, b, c, d, d2, Mode > 0 ? &gd : nullptr, Mode == 2 ? &hd : nullptr);
    }
    else {
      gpuDistance::edge_edge_distance2(flag, a, b, c, d, d2);
      if constexpr (Mode > 0)
        gpuDistance::edge_edge_distance2_gradient(flag, a, b, c, d, gd);
      if constexpr (Mode == 2)
        gpuDistance::edge_edge_distance2_hessian(flag, a, b, c, d, hd);
    }
  }
  double bv, bp, bpp;
  barrier(d2, dhat2, bv, bp, bpp);
  const double wk = kappa * pair.weight;
  double mollifier = 1.0;
  Vec12 gm = Vec12::Zero();
  Mat12 hm = Mat12::Zero();
  if (pair.kind == 1 && epsEE > 0.0) {
    gpuDistance::edge_edge_mollifier(a, b, c, d, epsEE, mollifier);
    if constexpr (Mode > 0)
      gpuDistance::edge_edge_mollifier_gradient(a, b, c, d, epsEE, gm);
    if constexpr (Mode == 2)
      gpuDistance::edge_edge_mollifier_hessian(a, b, c, d, epsEE, hm);
  }
  if constexpr (Mode == 0) {
    atomicAdd(energy, wk * mollifier * bv);
  }
  if constexpr (Mode == 1) {
    const Vec12 g = wk * (gm * bv + mollifier * bp * gd);
    for (int i = 0; i < 12; ++i)
      atomicAdd(gradient + 3 * pair.vertices[i / 3] + i % 3, g[i]);
  }
  if constexpr (Mode == 2) {
    const Vec12 gb = bp * gd;
    Mat12 h = wk * (bv * hm + gm * gb.transpose() + gb * gm.transpose() + mollifier * (bpp * gd * gd.transpose() + bp * hd));
    if (projectPSD)
      uipc::backend::cuda::make_spd<12>(h);
    Eigen::Map<Mat12>(hessians + 144 * index) = h;
  }
}

__global__ void hessianVectorKernel(const Pair *pairs, const double *hessians, int count,
  const double *v, double *out)
{
  const int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index >= count)
    return;
  const Pair pair = pairs[index];
  Vec12 localV;
  for (int i = 0; i < 12; ++i)
    localV[i] = v[3 * pair.vertices[i / 3] + i % 3];
  const Vec12 result = Eigen::Map<const Mat12>(hessians + 144 * index) * localV;
  for (int i = 0; i < 12; ++i)
    atomicAdd(out + 3 * pair.vertices[i / 3] + i % 3, result[i]);
}

__global__ void diagonalKernel(const Pair *pairs, const double *hessians, int count,
  MapView map, double *diagonal)
{
  const int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index >= count)
    return;
  const Pair pair = pairs[index];
  const Eigen::Map<const Mat12> h(hessians + 144 * index);
  // Compute diag(W^T H W) exactly, including cross terms when several surface
  // coordinates depend on the same simulation DOF.
  for (int i = 0; i < 12; ++i) {
    const int rowI = 3 * pair.vertices[i / 3] + i % 3;
    for (int ki = map.offsets[rowI]; ki < map.offsets[rowI + 1]; ++ki) {
      double sum = 0.0;
      for (int j = 0; j < 12; ++j) {
        const int rowJ = 3 * pair.vertices[j / 3] + j % 3;
        for (int kj = map.offsets[rowJ]; kj < map.offsets[rowJ + 1]; ++kj)
          if (map.columns[ki] == map.columns[kj])
            sum += map.values[kj] * h(i, j);
      }
      atomicAdd(diagonal + map.columns[ki], map.values[ki] * sum);
    }
  }
}

__device__ void atomicMinPositive(double *address, double value)
{
  auto *bits = reinterpret_cast<unsigned long long *>(address);
  unsigned long long old = *bits;
  while (value < __longlong_as_double(old)) {
    const auto expected = old;
    old = atomicCAS(bits, expected, __double_as_longlong(value));
    if (old == expected)
      break;
  }
}

__global__ void ccdKernel(const Pair *pairs, int count, const double *x, const double *dx,
  const double *errors, double slackness, double *alpha)
{
  const int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index >= count)
    return;
  const Pair pair = pairs[index];
  Vec3 points[4], motion[4], error[4];
  for (int i = 0; i < 4; ++i) {
    points[i] = position(x, pair.vertices[i]);
    motion[i] = position(dx, pair.vertices[i]);
    error[i] = position(errors, pair.vertices[i]);
  }
  const double toi = CudaIPCContinuousCollision::maximumStep(points, motion, error, pair.kind);
  if (toi < 1.0)
    atomicMinPositive(alpha, fmax(0.0, __dmul_rd(toi, slackness)));
}

template<typename T>
void upload(DeviceBuffer<T> &device, const std::vector<T> &host)
{
  device.upload(host.data(), host.size());
}
}  // namespace

class CudaEmbeddedSurfaceIPCPotentialEnergy::Impl
{
public:
  SurfaceIPCCore::Parameters parameters;
  int numDOFs = 0, numVertices = 0, numTriangles = 0, numEdges = 0;
  int activeCount = 0, hessianCount = 0;
  bool hessianReady = false;
  DeviceBuffer<int> triangles, edges, mapOffsets, mapColumns;
  DeviceBuffer<uint8_t> deformable;
  DeviceBuffer<double> vertexArea, triangleArea, edgeLength, rest, mapValues;
  DeviceBuffer<double> positions, trialPositions, displacements, ccdErrors, surfaceScratch, scalar, changeScalar, hessians;
  DeviceBuffer<int> count, invalid;
  DeviceBuffer<Pair> pairs, hessianPairs;
  std::unique_ptr<CudaSurfaceIPCBroadPhase> broadPhase;

  MapView map() const
  {
    return { mapOffsets.data(), mapColumns.data(), mapValues.data(), 3 * numVertices };
  }
  TopologyView topology() const
  {
    return { triangles.data(), edges.data(), deformable.data(), vertexArea.data(),
      triangleArea.data(), edgeLength.data(), numVertices, numTriangles, numEdges };
  }
  void updatePositions(const double *x)
  {
    mapKernel<<<blocks(3 * numVertices), blockSize>>>(map(), x, rest.data(), positions.data());
    checkCuda(cudaGetLastError(), "IPC surface mapping");
  }
  int findPairs(bool swept, bool rejectInvalid = true, const double *errors = nullptr)
  {
    try {
      broadPhase->detect(positions.data(), swept ? displacements.data() : nullptr, swept ? 0.0 : parameters.dhat, errors);
    }
    catch (const std::domain_error &) {
      if (rejectInvalid)
        throw;
      return -1;
    }
    const int candidates = static_cast<int>(broadPhase->getNumCandidates());
    pairs.resize(candidates);
    count.zero();
    invalid.zero();
    if (candidates) {
      filterPairsKernel<<<blocks(candidates), blockSize>>>(topology(), broadPhase->getDeviceCandidates(),
        candidates, positions.data(), swept, parameters.dhat, pairs.data(), count.data(), invalid.data());
      checkCuda(cudaGetLastError(), "IPC GPU narrow phase");
    }
    if (swept)
      return candidates;
    int required = 0, bad = 0;
    count.download(&required, 1);
    invalid.download(&bad, 1);
    if (bad) {
      if (rejectInvalid)
        throw std::runtime_error("CUDA IPC encountered a nonpositive or nonfinite primitive distance.");
      return -1;
    }
    return required;
  }
  void prepare(const double *x)
  {
    updatePositions(x);
    activeCount = findPairs(false);
  }
};

CudaEmbeddedSurfaceIPCPotentialEnergy::CudaEmbeddedSurfaceIPCPotentialEnergy(
  const EigenSupport::MXd &vertices, const EigenSupport::MXi &faces,
  const EigenSupport::SpMatD &mapping, const SurfaceIPCCore::Parameters &parameters):
  CudaEmbeddedSurfaceIPCPotentialEnergy(vertices, faces, mapping,
    std::vector<uint8_t>(static_cast<std::size_t>(vertices.rows()), uint8_t{ 1 }), parameters)
{
}

CudaEmbeddedSurfaceIPCPotentialEnergy::CudaEmbeddedSurfaceIPCPotentialEnergy(
  const EigenSupport::MXd &vertices, const EigenSupport::MXi &faces,
  const EigenSupport::SpMatD &mapping, const std::vector<uint8_t> &mask,
  const SurfaceIPCCore::Parameters &parameters):
  impl_(std::make_unique<Impl>())
{
  // Share CPU topology validation and quadrature measures at initialization.
  SurfaceIPCCore validator(parameters);
  validator.setMesh(vertices, faces, mask);
  if (parameters.frictionCoeff != 0.0)
    throw std::invalid_argument("CUDA IPC currently supports frictionCoeff=0 only.");
  if (mapping.rows() != 3 * vertices.rows() || mapping.cols() <= 0 ||
    mapping.rows() > INT_MAX || mapping.cols() > INT_MAX || mapping.nonZeros() > INT_MAX)
    throw std::invalid_argument("CUDA IPC displacement map has invalid dimensions or exceeds 32-bit indexing.");
  SurfaceIPCTopology topology;
  topology.setMesh(vertices, faces, mask);
  auto &data = *impl_;
  data.parameters = parameters;
  data.numDOFs = static_cast<int>(mapping.cols());
  data.numVertices = static_cast<int>(vertices.rows());
  data.numTriangles = static_cast<int>(topology.triangles.size());
  data.numEdges = static_cast<int>(topology.edges.size());
  std::vector<int> indices;
  for (const auto &triangle : topology.triangles)
    indices.insert(indices.end(), triangle.begin(), triangle.end());
  upload(data.triangles, indices);
  indices.clear();
  for (const auto &edge : topology.edges)
    indices.insert(indices.end(), edge.begin(), edge.end());
  upload(data.edges, indices);
  upload(data.deformable, topology.vertexIsDeformable);
  upload(data.vertexArea, topology.vertexArea);
  upload(data.triangleArea, topology.triArea);
  upload(data.edgeLength, topology.edgeLength);
  std::vector<double> rest(static_cast<std::size_t>(3) * vertices.rows());
  for (int vertex = 0; vertex < data.numVertices; ++vertex)
    for (int axis = 0; axis < 3; ++axis)
      rest[3 * vertex + axis] = vertices(vertex, axis);
  upload(data.rest, rest);
  Eigen::SparseMatrix<double, Eigen::RowMajor, int> rowMap = mapping;
  rowMap.makeCompressed();
  for (int row = 0; row < rowMap.rows(); ++row) {
    for (decltype(rowMap)::InnerIterator entry(rowMap, row); entry; ++entry) {
      if (!std::isfinite(entry.value()))
        throw std::invalid_argument("CUDA IPC displacement map contains a nonfinite value.");
      if (!mask[static_cast<std::size_t>(row / 3)] && entry.value() != 0.0)
        throw std::invalid_argument("CUDA IPC external vertices require zero displacement-map rows.");
    }
  }
  data.mapOffsets.upload(rowMap.outerIndexPtr(), rowMap.rows() + 1);
  data.mapColumns.upload(rowMap.innerIndexPtr(), rowMap.nonZeros());
  data.mapValues.upload(rowMap.valuePtr(), rowMap.nonZeros());
  data.positions.resize(rest.size());
  data.trialPositions.resize(rest.size());
  data.displacements.resize(rest.size());
  data.ccdErrors.resize(rest.size());
  data.surfaceScratch.resize(rest.size());
  data.scalar.resize(1);
  data.changeScalar.resize(2);
  data.count.resize(1);
  data.invalid.resize(1);
  data.broadPhase = std::make_unique<CudaSurfaceIPCBroadPhase>(CudaSurfaceIPCTopologyView{
    data.triangles.data(), data.edges.data(), data.deformable.data(), data.numVertices, data.numTriangles, data.numEdges });
}

CudaEmbeddedSurfaceIPCPotentialEnergy::~CudaEmbeddedSurfaceIPCPotentialEnergy() = default;

int CudaEmbeddedSurfaceIPCPotentialEnergy::getNumDOFs() const
{
  return impl_->numDOFs;
}
std::size_t CudaEmbeddedSurfaceIPCPotentialEnergy::getNumActivePairs() const
{
  return impl_->activeCount;
}

double CudaEmbeddedSurfaceIPCPotentialEnergy::func(const double *x)
{
  auto &data = *impl_;
  data.updatePositions(x);
  data.activeCount = data.findPairs(false, false);
  if (data.activeCount < 0) {
    data.activeCount = 0;
    return std::numeric_limits<double>::infinity();
  }
  data.scalar.zero();
  if (data.activeCount) {
    evaluateKernel<0><<<blocks(data.activeCount), blockSize>>>(data.pairs.data(), data.activeCount,
      data.positions.data(), data.parameters.dhat * data.parameters.dhat, data.parameters.kappa,
      data.parameters.eps_ee, false, data.scalar.data(), nullptr, nullptr);
    checkCuda(cudaGetLastError(), "IPC GPU energy");
  }
  double energy;
  data.scalar.download(&energy, 1);
  return energy;
}

NonlinearOptimization::CudaPotentialEnergy::EnergyChange
CudaEmbeddedSurfaceIPCPotentialEnergy::energyChange(const double *base, const double *trial)
{
  auto &data = *impl_;
  data.updatePositions(base);
  mapDifferenceKernel<<<blocks(3 * data.numVertices), blockSize>>>(data.map(), base, trial,
    data.positions.data(), data.displacements.data(), data.trialPositions.data());
  checkCuda(cudaGetLastError(), "IPC energy-change surface mapping");
  try {
    data.broadPhase->detect(data.positions.data(), data.displacements.data(), data.parameters.dhat);
  }
  catch (const std::domain_error &) {
    return { std::numeric_limits<double>::infinity(), 0.0 };
  }
  const int candidates = static_cast<int>(data.broadPhase->getNumCandidates());
  data.pairs.resize(candidates);
  data.invalid.zero();
  data.changeScalar.zero();
  if (candidates) {
    filterPairsKernel<<<blocks(candidates), blockSize>>>(data.topology(), data.broadPhase->getDeviceCandidates(),
      candidates, data.positions.data(), true, data.parameters.dhat, data.pairs.data(), data.count.data(), data.invalid.data());
    energyChangeKernel<<<blocks(candidates), blockSize>>>(data.pairs.data(), candidates, data.positions.data(), data.trialPositions.data(), data.displacements.data(),
      data.parameters.dhat * data.parameters.dhat, data.parameters.kappa, data.parameters.eps_ee, data.changeScalar.data(), data.invalid.data());
    checkCuda(cudaGetLastError(), "IPC local energy changes");
  }
  double result[2];
  int invalid = 0;
  data.changeScalar.download(result, 2);
  data.invalid.download(&invalid, 1);
  if (invalid)
    return { std::numeric_limits<double>::infinity(), 0.0 };
  return { result[0], result[1] };
}

void CudaEmbeddedSurfaceIPCPotentialEnergy::addGradient(const double *x, double *gradient)
{
  auto &data = *impl_;
  data.prepare(x);
  data.surfaceScratch.zero();
  if (data.activeCount) {
    evaluateKernel<1><<<blocks(data.activeCount), blockSize>>>(data.pairs.data(), data.activeCount,
      data.positions.data(), data.parameters.dhat * data.parameters.dhat, data.parameters.kappa,
      data.parameters.eps_ee, false, nullptr, data.surfaceScratch.data(), nullptr);
    checkCuda(cudaGetLastError(), "IPC GPU gradient");
  }
  transposeMapKernel<<<blocks(3 * data.numVertices), blockSize>>>(data.map(), data.surfaceScratch.data(), gradient);
  checkCuda(cudaGetLastError(), "IPC gradient transpose mapping");
}

void CudaEmbeddedSurfaceIPCPotentialEnergy::prepareHessian(const double *x)
{
  prepareHessianImpl(x, impl_->parameters.projectHessianToPSD);
}

void CudaEmbeddedSurfaceIPCPotentialEnergy::prepareExactHessian(const double *x)
{
  prepareHessianImpl(x, false);
}

void CudaEmbeddedSurfaceIPCPotentialEnergy::prepareHessianImpl(const double *x, bool projectPSD)
{
  auto &data = *impl_;
  data.hessianReady = false;
  data.prepare(x);
  data.hessianCount = data.activeCount;
  data.hessianPairs.resize(data.hessianCount);
  data.hessians.resize(static_cast<std::size_t>(144) * data.hessianCount);
  if (data.hessianCount) {
    checkCuda(cudaMemcpyAsync(data.hessianPairs.data(), data.pairs.data(),
                sizeof(Pair) * data.hessianCount, cudaMemcpyDeviceToDevice),
      "IPC cache Hessian pairs");
    evaluateKernel<2><<<blocks(data.hessianCount), blockSize>>>(data.hessianPairs.data(), data.hessianCount,
      data.positions.data(), data.parameters.dhat * data.parameters.dhat, data.parameters.kappa,
      data.parameters.eps_ee, projectPSD, nullptr, nullptr, data.hessians.data());
    checkCuda(cudaGetLastError(), "IPC GPU Hessian");
  }
  data.hessianReady = true;
}

void CudaEmbeddedSurfaceIPCPotentialEnergy::addHessianVector(const double *v, double *out) const
{
  auto &data = *impl_;
  if (!data.hessianReady)
    throw std::logic_error("CUDA IPC Hessian was not prepared.");
  if (!data.hessianCount)
    return;
  mapKernel<<<blocks(3 * data.numVertices), blockSize>>>(data.map(), v, nullptr, data.displacements.data());
  data.surfaceScratch.zero();
  hessianVectorKernel<<<blocks(data.hessianCount), blockSize>>>(data.hessianPairs.data(), data.hessians.data(),
    data.hessianCount, data.displacements.data(), data.surfaceScratch.data());
  transposeMapKernel<<<blocks(3 * data.numVertices), blockSize>>>(data.map(), data.surfaceScratch.data(), out);
  checkCuda(cudaGetLastError(), "IPC mapped Hessian-vector product");
}

void CudaEmbeddedSurfaceIPCPotentialEnergy::addDiagonal(double *diagonal) const
{
  const auto &data = *impl_;
  if (!data.hessianReady)
    throw std::logic_error("CUDA IPC Hessian was not prepared.");
  if (!data.hessianCount)
    return;
  diagonalKernel<<<blocks(data.hessianCount), blockSize>>>(data.hessianPairs.data(), data.hessians.data(),
    data.hessianCount, data.map(), diagonal);
  checkCuda(cudaGetLastError(), "IPC mapped Hessian diagonal");
}

double CudaEmbeddedSurfaceIPCPotentialEnergy::computeMaxStepSize(const double *x, const double *dx) const
{
  auto &data = *impl_;
  mapCcdKernel<<<blocks(3 * data.numVertices), blockSize>>>(data.map(), x, dx, data.rest.data(),
    data.positions.data(), data.displacements.data(), data.ccdErrors.data());
  checkCuda(cudaGetLastError(), "IPC bounded CCD surface mapping");
  const int candidates = data.findPairs(true, true, data.ccdErrors.data());
  double alpha = 1.0;
  data.scalar.upload(&alpha, 1);
  if (candidates) {
    ccdKernel<<<blocks(candidates), blockSize>>>(data.pairs.data(), candidates, data.positions.data(),
      data.displacements.data(), data.ccdErrors.data(), data.parameters.slackness, data.scalar.data());
    checkCuda(cudaGetLastError(), "IPC GPU continuous collision detection");
  }
  data.scalar.download(&alpha, 1);
  return alpha;
}
}  // namespace pgo::Contact::CIPC
