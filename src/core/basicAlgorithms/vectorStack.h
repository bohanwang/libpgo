/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
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

namespace pgo
{
namespace BasicAlgorithms
{

// a stack/vector combination that holds multiple items in order and has an index
// that points to the current item
// useful for displaying history of an editing process
template<class T>
class VectorStack
{
public:
  VectorStack():
    i(-1) {}
  void clear()
  {
    vec.clear();
    i = -1;
  }

  void push(const T &newItem)
  {
    i = vec.size();
    vec.push_back(newItem);
  }
  void pop()
  {
    vec.pop_back();
    i = int(vec.size()) - 1;
  }

  void toBegin() { i = 0; }
  void toEnd() { i = int(vec.size()) - 1; }

  int size() const { return vec.size(); }
  bool empty() const { return vec.empty(); }

  const T &current() const
  {
    assert(i >= 0 && i < (int)vec.size());
    return vec[i];
  }
  T &current()
  {
    assert(i >= 0 && i < (int)vec.size());
    return vec[i];
  }
  int index() const { return this->i; }

  const T &top() const
  {
    assert(vec.size() > 0);
    return vec.back();
  }
  T &top()
  {
    assert(vec.size() > 0);
    return vec.back();
  }

  // ++index
  int &operator++()
  {
    if (i + 1 == (int)vec.size())
      return i;
    return ++i;
  }
  // --index
  int &operator--()
  {
    if (i <= 0)
      return i;
    return --i;
  }
  // index++
  int operator++(int)
  {
    if (i + 1 == (int)vec.size())
      return i;
    return i++;
  }
  // index--
  int operator--(int)
  {
    if (i <= 0)
      return i;
    return i--;
  }

protected:
  std::vector<T> vec;
  int i;
};

}  // namespace BasicAlgorithms
}  // namespace pgo
