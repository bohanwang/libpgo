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

#include "range.h"

#include <algorithm>
#include <functional>
#include <iterator>
#include <vector>
#include <cmath>
#include <map>

namespace pgo
{
namespace BasicAlgorithms
{

template<class Container>
int sizei(const Container &c)
{
  return static_cast<int>(c.size());
}

// return true if all elements in [first, last] == value
template<class InputIt, class T>
bool allOf(InputIt first, InputIt last, const T &value)
{
  auto l = [&](const typename std::iterator_traits<InputIt>::value_type &it) {
    return it == value;
  };
  return std::find_if_not(first, last, l) == last;
}

template<class Container, class Range>
void insertRange(Container &a, const Range &b)
{
  a.insert(b.begin(), b.end());
}

template<class Container, class Range>
void insertRangeToEnd(Container &a, const Range &b)
{
  a.insert(a.end(), b.begin(), b.end());
}

template<class Container, class Map>
void insertMapValues(Container &a, const Map &b)
{
  for (const auto &p : b)
    a.insert(p.second);
}

template<class Container, class Map>
void insertMapValuesToEnd(Container &a, const Map &b)
{
  for (const auto &p : b)
    a.insert(a.end(), p.second);
}

// remove duplicate elements in vec and sort them as a side effect
// elements of type T need a default constructor
template<class T>
void sortAndDeduplicate(std::vector<T> &vec)
{
  std::sort(vec.begin(), vec.end());
  auto newEnd = std::unique(vec.begin(), vec.end());
  vec.resize(std::distance(vec.begin(), newEnd));
}

template<class T>
void sortAndDeduplicate(std::vector<T> &vec, std::size_t startingIndex)
{
  std::sort(vec.begin() + startingIndex, vec.end());
  auto newEnd = std::unique(vec.begin() + startingIndex, vec.end());
  vec.resize(std::distance(vec.begin(), newEnd));
}

// remove duplicate elements in vec and sort them as a side effect
// use vector::erase to remove duplicate elements so that they don't need a default constructor
template<class T>
void sortAndDeduplicateWithErase(std::vector<T> &vec)
{
  std::sort(vec.begin(), vec.end());
  auto newEnd = std::unique(vec.begin(), vec.end());
  vec.erase(newEnd, vec.end());
}

// go through each index in range [first, last) and process the index by calling f(Index ID)
// different from std::for_each, which goes through an iterator range for(iter : range), and calls f(*iter)
template<class Index, class UnaryFunction>
void forEachIndex(Index first, Index last, UnaryFunction f)
{
  for (; first != last; ++first) {
    f(first);
  }
}

template<class Array>
auto convertArrayToFunction(const Array &arr) -> std::function<decltype(arr[0])(int)>
{
  return [&arr](int index) -> decltype(arr[index]) {
    return arr[index];
  };
}

// run binary search on a range to find whehter a val is in the range
template<class RamdomAccessRange, class T>
bool binarySearchFound(const RamdomAccessRange &sortedVector, const T &val)
{
  return std::binary_search(sortedVector.begin(), sortedVector.end(), val);
}

// if all values inside range can be found in the sorted vector by binary search
template<class RamdomAccessRange, class Range>
bool binarySearchFoundRange(const RamdomAccessRange &sortedVector, const Range &range)
{
  for (const auto &value : range)
    if (binarySearchFound(sortedVector, value) == false)
      return false;
  return true;
}

template<class ForwardIterator, class T>
ForwardIterator binarySearchFind(ForwardIterator first, ForwardIterator last, const T &val)
{
  first = std::lower_bound(first, last, val);
  return ((first != last && !(val < *first)) ? first : last);
}

template<class ForwardIterator, class T, class Compare>
ForwardIterator binarySearchFind(ForwardIterator first, ForwardIterator last, const T &val, Compare comp)
{
  first = std::lower_bound(first, last, val, comp);
  return ((first != last && !(comp(val, *first))) ? first : last);
}

// run binary search on a range to find the index of val in sortedVector
// return sortedVector.size() if not found
template<class RamdomAccessRange, class T>
std::size_t binarySearchGetIndex(const RamdomAccessRange &sortedVector, const T &val)
{
  auto it = std::lower_bound(sortedVector.begin(), sortedVector.end(), val);
  return ((it != sortedVector.end() && !(val < *it)) ? std::distance(sortedVector.begin(), it) : sortedVector.size());
}
template<class RamdomAccessRange, class T, class Compare>
std::size_t binarySearchGetIndex(const RamdomAccessRange &sortedVector, const T &val, Compare comp)
{
  auto it = std::lower_bound(sortedVector.begin(), sortedVector.end(), val, comp);
  return ((it != sortedVector.end() && !(comp(val, *it))) ? std::distance(sortedVector.begin(), it) : sortedVector.size());
}

// assume distance between keys can be measured using distance(const Key & a, const Key & b)
template<class Set, class T, class Distance>
typename Set::const_iterator findClosestElementInSet(const Set &set, const T &val, Distance dist)
{
  auto after = set.lower_bound(val);
  if (after == set.begin() || *after == val)
    return after;
  auto before = after;
  before--;
  if (dist(*after, val) < dist(val, *before))
    return after;
  return before;
}

template<class Map, class T, class Distance>
typename Map::const_iterator findClosestElementInMap(const Map &map, const T &val, Distance dist)
{
  auto after = map.lower_bound(val);
  if (after == map.begin() || after->first == val)
    return after;
  auto before = after;
  before--;
  if (dist(after->first, val) < dist(val, before->first))
    return after;
  return before;
}

template<class Set, class T>
typename Set::const_iterator findClosestRealInSet(const Set &set, const T &val)
{
  return findClosestElementInSet(set, val, [](const typename Set::key_type &va, const typename Set::key_type &vb) {
    auto dist = va - vb;
    return (dist < 0 ? -dist : dist);
  });
}

template<class Map, class T>
typename Map::const_iterator findClosestRealInMap(const Map &map, const T &val)
{
  return findClosestElementInMap(map, val, [](const typename Map::key_type &va, const typename Map::key_type &vb) {
    auto dist = va - vb;
    return (dist < 0 ? -dist : dist);
  });
}

// given a sorted vector of indices, find the indices that are not in the vector
template<class RamdomAccessRange, class Vector2>
void pushBackIntegersNotInSortedVector(int numElements, const RamdomAccessRange &vec, Vector2 &result)
{
  for (int i = 0; i < numElements; i++)
    if (binarySearchFound(vec, i) == false)
      result.push_back(i);
}

#if __cplusplus > 201402L
using std::clamp;
#else
template<class T>
T clamp(const T &a, const T &low, const T &high)
{
  if (a < low) {
    return low;
  }
  if (a > high) {
    return high;
  }
  return a;
}
#endif

// clamp a into [low, high]
template<class T>
void clampSelf(T &a, const T &low, const T &high)
{
  if (a < low) {
    a = low;
  }
  else if (a > high) {
    a = high;
  }
}

// whether a is in range of [first, last)
template<class T>
bool inRange(const T &a, const T &first, const T &last)
{
  return (first <= a && a < last);
}

// whether a is in closed range of [low, high]
template<class T>
bool inClosedRange(const T &a, const T &low, const T &high)
{
  return (low <= a && a <= high);
}

// return 1 if 0 < v, -1 if v < 0, 0 if v == 0
template<class T>
int signAsInt(const T &v)
{
  if (0 < v)
    return 1;
  else if (v < 0)
    return -1;
  else
    return 0;
}

// return 0.0 if v < 0.0
// useful when you compute a value which should be >= 0.0 but might gives negative values due to numerical errors
inline double sqrtSafe(double v)
{
  if (v < 0.0) {
    return 0.0;
  }
  return std::sqrt(v);
}

inline double acosSafe(double v)
{
  if (v < -1.0)
    return std::acos(-1.0);
  if (v > 1.0)
    return std::acos(1.0);
  return acos(v);
}

inline double asinSafe(double v)
{
  if (v < -1.0)
    return std::asin(-1.0);
  if (v > 1.0)
    return std::asin(1.0);
  return asin(v);
}

// for a range of elements [first, last) and a unary func on the element,
// find the element that gives max func(element)
template<class InputIt, class UnaryFunction>
InputIt maxFunctionValue(InputIt first, InputIt last, UnaryFunction func)
{
  if (first == last) {
    return last;
  }
  InputIt largest = first;
  auto value = func(*first);
  first++;
  for (; first != last; first++) {
    auto newValue = func(*first);
    if (value < newValue) {
      largest = first;
      value = newValue;
    }
  }
  return largest;
}

// scale data stored in inputRange to be in range of [low, high] and save it to outputRange
// inputRange and outputRange can be the same
template<class Range, typename T>
void normalizeRangeIntoBounds(const Range &inputRange, const T &low, const T &high, Range &outputRange)
{
  auto minMaxPair = std::minmax_element(inputRange.begin(), inputRange.end());
  T lowValue = *(minMaxPair.first);
  T span = *(minMaxPair.second) - *(minMaxPair.first);
  if (span != T(0)) {
    T scale = (high - low) / span;
    auto it = inputRange.begin();
    auto it2 = outputRange.begin();
    for (; it != inputRange.end(); ++it, ++it2) {
      *it2 = (scale * (*it - lowValue)) + low;
    }
  }
}

// compute median for [first, last)
// if first == last, return default value: value_type()
// computing median involves averaging two values, this is done here by dividing by value_type(2)
template<class RandomIt>
typename std::iterator_traits<RandomIt>::value_type median(RandomIt first, RandomIt last)
{
  typedef typename std::iterator_traits<RandomIt>::value_type value;
  if (first == last) {
    return value();
  }
  auto size = distance(first, last);
  auto nth = first + size / 2;
  nth_element(first, nth, last);
  if (size % 2 == 1) {
    return *nth;
  }

  // size % 2 == 0
  auto v0 = *nth;
  nth_element(first, nth - 1, last);
  return (v0 + *(nth - 1)) / value(2);
}

// remove elements from vector "inputVector" by given indices stored in "indices"
// input indices are sorted
template<class InputVector, class IndexRange>
void removeByIndices(InputVector &inputVector, const IndexRange &indices)
{
  int vecSize = static_cast<int>(inputVector.size());
  //  int indSize = distance(indices.begin(), indices.end());
  auto ID = indices.begin();
  int newEnd = 0;
  for (int i = 0; i < vecSize; i++) {
    if (ID != indices.end() && i == *ID) {
      ID++;
      continue;
    }
    if (i != newEnd) {
      inputVector[newEnd] = std::move(inputVector[i]);
    }
    newEnd++;
  }
  inputVector.resize(newEnd);
}

// The template function must be run on a sorted range from first to last.
// It functions like std::unique, which operates on a sorted range and moves duplicated elements to the back of the range,
// and finally returns the iterator pointing to the end of the deduplicated elements.
// However, reduceDuplicates not only moves duplicate elements to the back, but also calls a reduce operator on
// those duplicated elements.
// reduce is a binary function: void reduce(T & entryA, T & entryB), reduces data in entryA and entryB and store the result into entryA
// e.g. input buffer is [0, 1, 2.1, 2.2, 3, 4.1, 4.2, 4.3], and we only compare with the integer part of each number (e.g. 2.1 == 2.2)
// and reduce two elements by adding their fractional parts together (e.g. (2.1, 2.2) => 2.3)
// then the result is a buffer: [0, 1, 2.3, 3, 4.6, x, x, x], the function returns the iterator on the first x
// and reduce is called on (2.1,2.2) and (4.1,4.2), (4.3,4.3)
template<class ForwardIt, class BinaryReduce, class BinaryPredicate>
ForwardIt reduceDuplicates(ForwardIt first, ForwardIt last, BinaryReduce reduce, BinaryPredicate equal)
{
  if (first == last)
    return last;

  ForwardIt result = first;
  while (++first != last) {
    if (equal(*result, *first) == false)  // if *result and *first are not equal
    {
      ++result;                           // move result forward
      if (result != first)                // if result and first does not point to the same location, then
        *result = std::move(*first);      // move *first to *result
    }
    else                                  // *result and *first are equal
    {
      reduce(*result, *first);
    }
  }
  return ++result;
}

// the version of reduceDuplicates using std::equal_to as the binary equal compare functor
template<class ForwardIt, class BinaryReduce>
ForwardIt reduceDuplicates(ForwardIt first, ForwardIt last, BinaryReduce reduce)
{
  using T = decltype(*first);
  return reduceDuplicates(first, last, reduce, std::equal_to<const T &>());
}

inline void buildMapFromPairs(int numPairs, const std::pair<int, int> *pairs, std::map<int, std::vector<int>> &IDMap)
{
  for (int i = 0; i < numPairs; i++) {
    const auto &p = pairs[i];
    IDMap[p.first].push_back(p.second);
    IDMap[p.second].push_back(p.first);
  }
}

}  // namespace BasicAlgorithms
}  // namespace pgo
