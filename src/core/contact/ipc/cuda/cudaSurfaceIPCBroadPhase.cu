/*
copyright to Bohan Wang
*/

#include "cudaSurfaceIPCBroadPhase.h"
#include "cudaDeviceBuffer.h"

#include <libuipc/collision_detection/linear_bvh_kernels.cuh>
#include <cub/device/device_radix_sort.cuh>
#include <cub/device/device_reduce.cuh>
#include <cuda/atomic>
#include <cuda_runtime.h>
#include <math_constants.h>

#include <algorithm>
#include <climits>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace pgo::Contact::CIPC
{
namespace
{
using CudaUtilities::DeviceBuffer;
using CudaUtilities::checkCuda;
using Node = uipc::backend::cuda::lbvh::Node;
constexpr int blockSize = 128;
int blocks(int n)
{
  return n / blockSize + (n % blockSize != 0);
}

struct Bounds
{
  double lower[3], upper[3];
};

struct MergeBounds
{
  __host__ __device__ Bounds operator()(const Bounds &a, const Bounds &b) const
  {
    Bounds result;
    for (int axis = 0; axis < 3; ++axis) {
      result.lower[axis] = fmin(a.lower[axis], b.lower[axis]);
      result.upper[axis] = fmax(a.upper[axis], b.upper[axis]);
    }
    return result;
  }
};

__device__ bool overlaps(const Bounds &a, const Bounds &b)
{
  for (int axis = 0; axis < 3; ++axis)
    if (a.lower[axis] > b.upper[axis] || b.lower[axis] > a.upper[axis])
      return false;
  return true;
}

template<int Vertices>
__global__ void primitiveBoundsKernel(const double *x, const double *dx, const int *indices,
  const std::uint8_t *deformable, double padding, const double *errors, Bounds *bounds, int *moving, int *invalid, int count)
{
  const int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index >= count)
    return;
  Bounds result;
  int anyDeformable = 0;
  bool finite = true;
  for (int axis = 0; axis < 3; ++axis) {
    double minimum = CUDART_INF, maximum = -CUDART_INF;
    for (int i = 0; i < Vertices; ++i) {
      const int vertex = Vertices == 1 ? index : indices[Vertices * index + i];
      const double p = x[3 * vertex + axis];
      const double qLower = dx ? __dadd_rd(p, dx[3 * vertex + axis]) : p;
      const double qUpper = dx ? __dadd_ru(p, dx[3 * vertex + axis]) : p;
      finite &= isfinite(p) && isfinite(qLower) && isfinite(qUpper);
      const double error = errors ? errors[3 * vertex + axis] : 0.0;
      finite &= isfinite(error) && error >= 0.0;
      minimum = fmin(minimum, __dsub_rd(fmin(p, qLower), error));
      maximum = fmax(maximum, __dadd_ru(fmax(p, qUpper), error));
      anyDeformable |= deformable[vertex] != 0;
    }
    result.lower[axis] = nextafter(minimum - padding, -CUDART_INF);
    result.upper[axis] = nextafter(maximum + padding, CUDART_INF);
    finite &= isfinite(result.lower[axis]) && isfinite(result.upper[axis]);
  }
  if (!finite) {
    atomicExch(invalid, 1);
    result = {};  // Keep tree construction well-defined before host rejection.
  }
  bounds[index] = result;
  moving[index] = anyDeformable;
}

__global__ void mortonKeysKernel(const Bounds *bounds, const Bounds *scene,
  std::uint64_t *keys, int count)
{
  const int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index >= count)
    return;
  double coordinate[3];
  for (int axis = 0; axis < 3; ++axis) {
    const double length = scene->upper[axis] - scene->lower[axis];
    const double center = 0.5 * bounds[index].lower[axis] + 0.5 * bounds[index].upper[axis];
    coordinate[axis] = length > 0 && isfinite(length) ? (center - scene->lower[axis]) / length : 0.5;
  }
  const auto morton = uipc::backend::cuda::lbvh::mortonCode(coordinate[0], coordinate[1], coordinate[2]);
  keys[index] = (static_cast<std::uint64_t>(morton) << 32) | static_cast<unsigned>(index);
}

__global__ void initializeNodesKernel(Node *nodes, const std::uint64_t *keys,
  const Bounds *primitiveBounds, const int *primitiveMoving, Bounds *bounds, int *moving, int count)
{
  const int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index >= 2 * count - 1)
    return;
  const int object = index >= count - 1 ? static_cast<int>(keys[index - count + 1] & 0xffffffffULL) : -1;
  nodes[index] = { -1, -1, -1, object };
  if (object >= 0) {
    bounds[index] = primitiveBounds[object];
    moving[index] = primitiveMoving[object];
  }
}

__global__ void mergeNodesKernel(const Node *nodes, Bounds *bounds, int *moving, int *arrivals, int count)
{
  const int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index >= count)
    return;
  int parent = nodes[index + count - 1].parent;
  while (parent >= 0) {
    // Release publishes the completed child bounds; acquire makes the first
    // child's bounds visible to the second arrival. No polling/spin loop.
    ::cuda::atomic_ref<int, ::cuda::thread_scope_device> flag(arrivals[parent]);
    if (flag.fetch_add(1, ::cuda::memory_order_acq_rel) == 0)
      return;
    const auto node = nodes[parent];
    bounds[parent] = MergeBounds{}(bounds[node.left], bounds[node.right]);
    moving[parent] = moving[node.left] | moving[node.right];
    parent = node.parent;
  }
}

__global__ void escapeNodesKernel(const Node *nodes, int *escape, int count)
{
  const int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index >= 2 * count - 1)
    return;
  int node = index, parent = nodes[node].parent;
  while (parent >= 0 && nodes[parent].right == node) {
    node = parent;
    parent = nodes[node].parent;
  }
  escape[index] = parent >= 0 ? nodes[parent].right : -1;
}

struct TreeView
{
  const Node *nodes;
  const Bounds *bounds;
  const int *moving;
  const int *escape;
  int count;
};

template<bool EdgeEdge>
__global__ void queryKernel(CudaSurfaceIPCTopologyView topology, TreeView tree,
  const Bounds *queries, const int *queryMoving, int queryCount,
  CudaSurfaceIPCCandidate *candidates, unsigned long long capacity, unsigned long long *count, const int *invalid)
{
  const int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index >= queryCount || tree.count == 0 || *invalid)
    return;
  const auto query = queries[index];
  const bool moving = queryMoving[index] != 0;
  int node = 0;
  while (node >= 0) {
    const auto current = tree.nodes[node];
    if ((!moving && !tree.moving[node]) || !overlaps(query, tree.bounds[node])) {
      node = tree.escape[node];
      continue;
    }
    if (current.object < 0) {
      node = current.left;
      continue;
    }
    const int other = current.object;
    bool adjacent = false;
    if constexpr (EdgeEdge) {
      adjacent = other <= index;
      for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
          adjacent |= topology.edges[2 * index + i] == topology.edges[2 * other + j];
    }
    else {
      for (int i = 0; i < 3; ++i)
        adjacent |= index == topology.triangles[3 * other + i];
    }
    if (!adjacent) {
      const auto slot = atomicAdd(count, 1ULL);
      if (slot < capacity)
        candidates[slot] = { EdgeEdge ? 1 : 0, index, other };
    }
    node = tree.escape[node];
  }
}

class Tree
{
public:
  DeviceBuffer<Bounds> primitiveBounds;
  DeviceBuffer<int> primitiveMoving;

  void build(int count)
  {
    count_ = count;
    if (count == 0)
      return;
    const int nodeCount = 2 * count - 1;
    nodes_.resize(nodeCount);
    bounds_.resize(nodeCount);
    moving_.resize(nodeCount);
    escape_.resize(nodeCount);
    keys_.resize(count);
    sortedKeys_.resize(count);
    scene_.resize(1);
    arrivals_.resize(count - 1);
    const double infinity = std::numeric_limits<double>::infinity();
    Bounds empty{ { infinity, infinity, infinity }, { -infinity, -infinity, -infinity } };
    std::size_t bytes = 0;
    checkCuda(cub::DeviceReduce::Reduce(nullptr, bytes, primitiveBounds.data(), scene_.data(), count, MergeBounds{}, empty), "IPC BVH reduction size");
    scratch_.resize(bytes);
    checkCuda(cub::DeviceReduce::Reduce(scratch_.data(), bytes, primitiveBounds.data(), scene_.data(), count, MergeBounds{}, empty), "IPC BVH scene bounds");
    mortonKeysKernel<<<blocks(count), blockSize>>>(primitiveBounds.data(), scene_.data(), keys_.data(), count);
    bytes = 0;
    checkCuda(cub::DeviceRadixSort::SortKeys(nullptr, bytes, keys_.data(), sortedKeys_.data(), count), "IPC BVH sort size");
    scratch_.resize(bytes);
    checkCuda(cub::DeviceRadixSort::SortKeys(scratch_.data(), bytes, keys_.data(), sortedKeys_.data(), count), "IPC BVH Morton sort");
    initializeNodesKernel<<<blocks(nodeCount), blockSize>>>(nodes_.data(), sortedKeys_.data(), primitiveBounds.data(),
      primitiveMoving.data(), bounds_.data(), moving_.data(), count);
    if (count > 1) {
      uipc::backend::cuda::lbvh::buildInternalNodes<<<blocks(count - 1), blockSize>>>(nodes_.data(), sortedKeys_.data(), count);
      arrivals_.zero();
      mergeNodesKernel<<<blocks(count), blockSize>>>(nodes_.data(), bounds_.data(), moving_.data(), arrivals_.data(), count);
    }
    escapeNodesKernel<<<blocks(nodeCount), blockSize>>>(nodes_.data(), escape_.data(), count);
    checkCuda(cudaGetLastError(), "IPC BVH construction");
  }

  TreeView view() const { return { nodes_.data(), bounds_.data(), moving_.data(), escape_.data(), count_ }; }

private:
  int count_ = 0;
  DeviceBuffer<Node> nodes_;
  DeviceBuffer<Bounds> bounds_, scene_;
  DeviceBuffer<int> moving_, escape_, arrivals_;
  DeviceBuffer<std::uint64_t> keys_, sortedKeys_;
  DeviceBuffer<char> scratch_;
};
}  // namespace

class CudaSurfaceIPCBroadPhase::Impl
{
public:
  CudaSurfaceIPCTopologyView topology;
  Tree triangles, edges;
  DeviceBuffer<Bounds> points;
  DeviceBuffer<int> pointMoving, invalid;
  DeviceBuffer<unsigned long long> count;
  DeviceBuffer<CudaSurfaceIPCCandidate> candidates;
  std::size_t candidateCount = 0;
};

CudaSurfaceIPCBroadPhase::CudaSurfaceIPCBroadPhase(CudaSurfaceIPCTopologyView topology):
  impl_(std::make_unique<Impl>())
{
  if (topology.numVertices < 0 || topology.numTriangles < 0 || topology.numEdges < 0 ||
    topology.numVertices > INT_MAX / 3 || topology.numTriangles > INT_MAX / 3 || topology.numEdges > INT_MAX / 2)
    throw std::length_error("CUDA IPC topology exceeds BVH node indexing.");
  auto &data = *impl_;
  data.topology = topology;
  data.points.resize(topology.numVertices);
  data.pointMoving.resize(topology.numVertices);
  data.triangles.primitiveBounds.resize(topology.numTriangles);
  data.triangles.primitiveMoving.resize(topology.numTriangles);
  data.edges.primitiveBounds.resize(topology.numEdges);
  data.edges.primitiveMoving.resize(topology.numEdges);
  data.invalid.resize(1);
  data.count.resize(1);
  data.candidates.resize(1024);
}

CudaSurfaceIPCBroadPhase::~CudaSurfaceIPCBroadPhase() = default;

void CudaSurfaceIPCBroadPhase::detect(const double *x, const double *dx, double padding, const double *errors)
{
  if (!(padding >= 0.0) || !std::isfinite(padding))
    throw std::invalid_argument("CUDA IPC BVH padding must be finite and nonnegative.");
  auto &data = *impl_;
  const auto topology = data.topology;
  data.candidateCount = 0;
  data.invalid.zero();
  if (topology.numVertices)
    primitiveBoundsKernel<1><<<blocks(topology.numVertices), blockSize>>>(x, dx, nullptr, topology.deformable, padding, errors,
      data.points.data(), data.pointMoving.data(), data.invalid.data(), topology.numVertices);
  if (topology.numTriangles)
    primitiveBoundsKernel<3><<<blocks(topology.numTriangles), blockSize>>>(x, dx, topology.triangles, topology.deformable, padding, errors,
      data.triangles.primitiveBounds.data(), data.triangles.primitiveMoving.data(), data.invalid.data(), topology.numTriangles);
  if (topology.numEdges)
    primitiveBoundsKernel<2><<<blocks(topology.numEdges), blockSize>>>(x, dx, topology.edges, topology.deformable, padding, errors,
      data.edges.primitiveBounds.data(), data.edges.primitiveMoving.data(), data.invalid.data(), topology.numEdges);
  checkCuda(cudaGetLastError(), "IPC primitive bounds");
  data.triangles.build(topology.numTriangles);
  data.edges.build(topology.numEdges);
  for (;;) {
    data.count.zero();
    if (topology.numVertices && topology.numTriangles)
      queryKernel<false><<<blocks(topology.numVertices), blockSize>>>(topology, data.triangles.view(), data.points.data(),
        data.pointMoving.data(), topology.numVertices, data.candidates.data(), data.candidates.size(), data.count.data(), data.invalid.data());
    if (topology.numEdges)
      queryKernel<true><<<blocks(topology.numEdges), blockSize>>>(topology, data.edges.view(), data.edges.primitiveBounds.data(),
        data.edges.primitiveMoving.data(), topology.numEdges, data.candidates.data(), data.candidates.size(), data.count.data(), data.invalid.data());
    checkCuda(cudaGetLastError(), "IPC BVH query");
    unsigned long long required = 0;
    int invalid = 0;
    data.count.download(&required, 1);
    data.invalid.download(&invalid, 1);
    if (invalid)
      throw std::domain_error("CUDA IPC BVH encountered a nonfinite position, motion, or bound.");
    if (required > INT_MAX)
      throw std::length_error("CUDA IPC candidates exceed supported pair indexing.");
    if (required <= data.candidates.size()) {
      data.candidateCount = static_cast<std::size_t>(required);
      return;
    }
    data.candidates.resize(std::min<std::size_t>(required + required / 2 + 1, INT_MAX));
  }
}

const CudaSurfaceIPCCandidate *CudaSurfaceIPCBroadPhase::getDeviceCandidates() const
{
  return impl_->candidates.data();
}

std::size_t CudaSurfaceIPCBroadPhase::getNumCandidates() const
{
  return impl_->candidateCount;
}
}  // namespace pgo::Contact::CIPC
