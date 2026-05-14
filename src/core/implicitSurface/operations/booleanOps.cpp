#include "operations/booleanOps.h"

#include <algorithm>
#include <stdexcept>

namespace pgo::ImplicitSurface {

void applyBoolean(const DenseGrid &a, const DenseGrid &b, BooleanOp op, DenseGrid &out)
{
  if (a.gridSpec() != b.gridSpec())
    throw std::runtime_error("applyBoolean: grid spec mismatch between inputs");

  if (a.gridSpec() != out.gridSpec())
    throw std::runtime_error("applyBoolean: out grid spec does not match inputs");

  if (&out == &a || &out == &b)
    throw std::runtime_error("applyBoolean: out must not alias a or b");

  const int total = a.size();

  switch (op) {
    case BooleanOp::Union:
      for (int i = 0; i < total; ++i)
        out[i] = std::min(a[i], b[i]);
      break;

    case BooleanOp::Intersection:
      for (int i = 0; i < total; ++i)
        out[i] = std::max(a[i], b[i]);
      break;

    case BooleanOp::Difference:
      for (int i = 0; i < total; ++i)
        out[i] = std::max(a[i], -b[i]);
      break;
  }
}

}  // namespace pgo::ImplicitSurface
