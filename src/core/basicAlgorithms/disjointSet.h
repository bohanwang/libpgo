/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "basicAlgorithms" library , Copyright (C) 2007 CMU, 2009 MIT, 2018 USC*
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Yili Zhao, Jernej Barbic                                *
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
/*
  "Disjoint set", a union-find data structure
  that keeps track of a set of elements partitioned into a number of disjoint (non-overlapping) subsets.

  Operations:
  1) MakeSet: Make each set containing only a given element.
  2) Find: Determine which set a particular element is in. Also useful for determining if two elements are in the same set.
  3) Union: Combine or merge two sets into a single set.

  Heuristics such as path compression etc. have been implemented.
*/

#include <vector>
#include <map>

namespace pgo
{
namespace BasicAlgorithms
{

// implementation of a disjoint set on a fixed array of elements with continuous IDs
class DisjointSet
{
public:
  enum
  {
    NO_SUCCESSOR = -1
  };

  // constructor
  DisjointSet() {}
  // initialize num elements and each elemeent is in its own set
  // 'num' is the total number of elements
  DisjointSet(int num) { extendSizeTo(num); }

  void clear()
  {
    parent.clear();
    depth.clear();
  }

  int size() const { return static_cast<int>(parent.size()); }

  // extend the data structure to have size 'num'
  // previous set relationship is maintained
  void extendSizeTo(int num);

  // makes each element be its own set (already done in the constructor)
  void makeSet();

  // returns the representative of the set that x belongs to (path compression is implemented in this function as well)
  int findSet(int x) const;

  // merge two sets (the smaller one will be absorbed by the larger one)
  // x and y can be arbitrary elements, they need not be representatives
  void unionSet(int x, int y);
  template<class IntIterator>
  void unionRange(IntIterator itBegin, IntIterator itEnd);
  template<class IntRange>
  void unionRange(IntRange range)
  {
    unionRange(range.begin(), range.end());
  }

  // create a mapping from elementID to new continuous IDs for each set
  // returned NewIDMapping is: [0, size()) -> [0, #set)
  std::vector<int> createOldToNewIDMapping() const;

protected:
  mutable std::vector<int> parent;
  std::vector<int> depth;
};

// implementation of a disjoint set on a dynamic set of elements with possibly incontinuous IDs
class DisjointSetDynamic
{
public:
  DisjointSetDynamic() {}

  void clear()
  {
    dset.clear();
    dsetInputID.clear();
  }

  int size() const { return static_cast<int>(dsetInputID.size()); }
  int numElements() const { return static_cast<int>(dsetInputID.size()); }

  // return the setID of the set this element belongs to
  // if this element has not been added before, it will be added and form an individual set
  int findSet(int elementID);

  // union the sets the two input elements (x, y) belong to
  // if any one of the input element is not added before, it will be added and form an individual set
  void unionSet(int elementIDx, int elementIDy);

  // return all remaining sets with their elementIDs inside
  std::vector<std::vector<int>> getAllSets() const;

  // return a mapping: elementID -> new continuous set ID
  //                   [0, numElements()) -> [0, #set)
  void createElementToNewSetIDMapping(std::map<int, int> &element2setID) const;

protected:
  int getSetID(int elementID);

  DisjointSet dset;                // dsetInoutID -> setID
  std::map<int, int> dsetInputID;  // elementID -> dsetInoutID
};

template<class IntIterator>
void DisjointSet::unionRange(IntIterator itBegin, IntIterator itEnd)
{
  auto jt = itBegin;
  jt++;
  for (; jt != itEnd; jt++)
    unionSet(*itBegin, *jt);
}

}  // namespace BasicAlgorithms
}  // namespace pgo