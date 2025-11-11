#pragma once

#include <vector>
#include <array>

namespace pgo
{
namespace AnimationIO
{
void dumpABC(const char *filename, const char *name,
  const std::vector<float> &positions,
  const std::vector<std::vector<float>> &displacements,
  const std::vector<std::vector<int>> &triangles,
  const std::vector<float> *uv = nullptr);
}
}  // namespace pgo