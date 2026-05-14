#pragma once

#include "field/gridSpec.h"
#include "EigenSupport.h"

#include <vector>

namespace pgo::ImplicitSurface {

class DenseGrid {
public:
  explicit DenseGrid(const GridSpec &spec);

  const GridSpec &gridSpec() const { return spec_; }
  [[nodiscard]] int resolution() const { return spec_.resolution; }
  [[nodiscard]] int size() const { return static_cast<int>(values_.size()); }

  double &operator[](int index) { return values_[index]; }
  const double &operator[](int index) const { return values_[index]; }

  double &at(int x, int y, int z);
  const double &at(int x, int y, int z) const;

  const double *data() const { return values_.data(); }
  double *data() { return values_.data(); }

  void fill(double value);
  void setZero();

private:
  GridSpec spec_;
  std::vector<double> values_;
};

}  // namespace pgo::ImplicitSurface
