/*************************************************************************
 *                                                                       *
 *                                                                       *
 * "basicAlgorithms" library , Copyright (C) 2018 USC                    *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Yijing Li, Jernej Barbic                                *
 * http://www.jernejbarbic.com/vega                                      *
 *                                                                       *
 * Research: Jernej Barbic, Hongyi Xu, Yijing Li,                        *
 *           Danyong Zhao, Bohan Wang,                                   *
 *           Fun Shing Sin, Daniel Schroeder,                            *
 *           Doug L. James, Jovan Popovic                                *
 *                                                                       *
 * Funding: National Science Foundation, Link Foundation,                *
 *          Singapore-MIT GAMBIT Game Lab,                               *
 *          Zumberge Research and Innovation Fund at USC,                *
 *          Sloan Foundation, Okawa Foundation,                          *
 *          USC Annenberg Foundation                                     *
 *                                                                       *
 * This library is free software; you can redistribute it and/or         *
 * modify it under the terms of the BSD-style license that is            *
 * included with this library in the file LICENSE.txt                    *
 *                                                                       *
 * This library is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the file     *
 * LICENSE.TXT for more details.                                         *
 *                                                                       *
 *************************************************************************/

#pragma once

#include <vector>
#include <array>

namespace pgo
{
namespace BasicAlgorithms
{
// a template class to refer to an external array of basic values
template<typename T>
class ArrayRef
{
public:
  typedef T value_type;

  ArrayRef() {}  // empty array
  ArrayRef(int numElements, const T *elements):
    n(numElements), v(elements) {}
  template<class Allocator>
  ArrayRef(const std::vector<T, Allocator> &vec):
    n(static_cast<int>(vec.size())), v(vec.data())
  {
  }
  template<std::size_t N>
  ArrayRef(const std::array<T, N> &arr):
    n(arr.size()), v(arr.data())
  {
  }

  int size() const { return n; }
  const T *data() const { return v; }
  const T *begin() const { return v; }
  const T *end() const { return v + n; }

  const T &operator[](int ID) const { return v[ID]; }

protected:
  int n = 0;
  const T *v = nullptr;
};

template<typename T>
ArrayRef<T> makeArrayRef(int numElements, const T *elements)
{
  return ArrayRef<T>(numElements, elements);
}

template<typename T, class Allocator>
ArrayRef<T> makeArrayRef(const std::vector<T, Allocator> &vec)
{
  return ArrayRef<T>(vec);
}

template<typename T, std::size_t N>
ArrayRef<T> makeArrayRef(const std::array<T, N> &arr)
{
  return ArrayRef<T>(arr);
}
}  // namespace BasicAlgorithms
}  // namespace pgo