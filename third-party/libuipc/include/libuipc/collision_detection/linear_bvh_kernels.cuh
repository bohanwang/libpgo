// Selected and adapted from libuipc, Apache-2.0.
// Upstream 9c748a759bb9952a7aeedd59076a70591dc61d6c:
// src/backends/cuda/collision_detection/linear_bvh.h (node layout, MuGdxy)
// src/backends/cuda/collision_detection/details/linear_bvh.inl
// Morton expansion, radix-tree range/split, and internal-node construction.
// Raw pointers replace cuda_tool views; indices are signed; Morton coordinates
// are double and clamp to 1023. Keys include a unique primitive ID. Bounds,
// allocation, reduction, sorting and stackless traversal live in PGO.
#pragma once

#include <cuda_runtime.h>
#include <cstdint>

namespace uipc::backend::cuda::lbvh
{
struct Node
{
  int parent, left, right, object;
};

__device__ inline unsigned expandBits(unsigned value)
{
  value = (value * 0x00010001u) & 0xFF0000FFu;
  value = (value * 0x00000101u) & 0x0F00F00Fu;
  value = (value * 0x00000011u) & 0xC30C30C3u;
  value = (value * 0x00000005u) & 0x49249249u;
  return value;
}

__device__ inline unsigned mortonCode(double x, double y, double z)
{
  const unsigned xx = expandBits(static_cast<unsigned>(fmin(1023.0, fmax(0.0, x * 1024.0))));
  const unsigned yy = expandBits(static_cast<unsigned>(fmin(1023.0, fmax(0.0, y * 1024.0))));
  const unsigned zz = expandBits(static_cast<unsigned>(fmin(1023.0, fmax(0.0, z * 1024.0))));
  return xx * 4 + yy * 2 + zz;
}

__device__ inline int commonUpperBits(std::uint64_t a, std::uint64_t b)
{
  return __clzll(a ^ b);
}

__device__ inline int2 determineRange(const std::uint64_t *keys, int count, int index)
{
  if (index == 0)
    return make_int2(0, count - 1);
  const auto key = keys[index];
  const int leftDelta = commonUpperBits(key, keys[index - 1]);
  const int rightDelta = commonUpperBits(key, keys[index + 1]);
  const int direction = rightDelta > leftDelta ? 1 : -1;
  const int minimumDelta = min(leftDelta, rightDelta);
  // Wide search arithmetic avoids overflow near the supported index limit.
  long long maximum = 2;
  for (;;) {
    const long long other = index + direction * maximum;
    if (other < 0 || other >= count || commonUpperBits(key, keys[other]) <= minimumDelta)
      break;
    maximum *= 2;
  }
  long long length = 0;
  for (long long step = maximum / 2; step >= 1; step /= 2) {
    const long long other = index + direction * (length + step);
    if (other >= 0 && other < count && commonUpperBits(key, keys[other]) > minimumDelta)
      length += step;
  }
  const int other = static_cast<int>(index + direction * length);
  return make_int2(min(index, other), max(index, other));
}

__device__ inline int findSplit(const std::uint64_t *keys, int first, int last)
{
  const auto firstCode = keys[first];
  const auto lastCode = keys[last];
  if (firstCode == lastCode)
    return first + (last - first) / 2;
  const int delta = commonUpperBits(firstCode, lastCode);
  int split = first, stride = last - first;
  do {
    stride = (stride + 1) / 2;
    const int middle = split + stride;
    if (middle < last && commonUpperBits(firstCode, keys[middle]) > delta)
      split = middle;
  } while (stride > 1);
  return split;
}

__global__ void buildInternalNodes(Node *nodes, const std::uint64_t *keys, int count)
{
  const int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index >= count - 1)
    return;
  const int2 range = determineRange(keys, count, index);
  const int split = findSplit(keys, range.x, range.y);
  const int left = split + (range.x == split ? count - 1 : 0);
  const int right = split + 1 + (range.y == split + 1 ? count - 1 : 0);
  nodes[index].left = left;
  nodes[index].right = right;
  nodes[left].parent = index;
  nodes[right].parent = index;
}
}  // namespace uipc::backend::cuda::lbvh
