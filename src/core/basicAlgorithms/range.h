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

#include <utility>

// used when you don't have a container but begin(), end() are needed for a range for loop

namespace pgo
{
namespace BasicAlgorithms
{

template<class InputIt>
class Range
{
public:
  template<class Container>
  Range(Container &container):
    a(container.begin()), b(container.end())
  {
  }
  Range(InputIt a, InputIt b):
    a(a), b(b) {}
  const InputIt &begin() const { return a; }
  const InputIt &end() const { return b; }

protected:
  InputIt a, b;
};

template<class InputIt>
Range<InputIt> makeRange(InputIt a, InputIt b)
{
  return Range<InputIt>(a, b);
}

template<class T>
Range<decltype(std::declval<T>().begin())> makeRange(T &container)
{
  return Range<decltype(std::declval<T>().begin())>(container);
}

}  // namespace BasicAlgorithms
}  // namespace pgo