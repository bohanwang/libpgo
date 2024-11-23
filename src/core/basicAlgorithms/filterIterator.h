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

#include <iterator>
#include <functional>

namespace pgo
{
namespace BasicAlgorithms
{

template<class BaseIterator>
struct FilterIterator
{
  typedef typename std::iterator_traits<BaseIterator>::value_type value_type;
  typedef std::function<bool(const value_type &)> Filter;

  FilterIterator() = default;
  FilterIterator(Filter filter, BaseIterator base, BaseIterator end = {}):
    iter_(base), end_(end), filter_(filter)
  {
    while (iter_ != end_ && !filter_(*iter_)) {
      iter_++;
    }
  }

  FilterIterator &operator++()
  {
    do {
      iter_++;
    } while (iter_ != end_ && !filter_(*iter_));
    return *this;
  }

  FilterIterator operator++(int)
  {
    FilterIterator copy = *this;
    ++(*this);
    return copy;
  }

  value_type &operator*() { return *iter_; }
  const value_type &operator*() const { return *iter_; }

  bool operator==(const FilterIterator &it2) const { return iter_ == it2.iter_; }
  bool operator!=(const FilterIterator &it2) const { return !(*this == it2); }

  bool operator==(const BaseIterator &it2) const { return iter_ == it2; }
  bool operator!=(const BaseIterator &it2) const { return !(*this == it2); }

  inline friend bool operator==(const BaseIterator &it, const FilterIterator &it2) { return it == it2.iter_; }
  inline friend bool operator!=(const BaseIterator &it, const FilterIterator &it2) { return it != it2.iter_; }

private:
  BaseIterator iter_, end_;
  Filter filter_;
};

template<class BaseIterator>
FilterIterator<BaseIterator> makeFilterIterator(typename FilterIterator<BaseIterator>::Filter f,
  BaseIterator base, BaseIterator end = {})
{
  return { f, base, end };
}

}  // namespace BasicAlgorithms
}  // namespace pgo
