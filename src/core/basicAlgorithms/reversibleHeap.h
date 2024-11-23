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

#include <map>
#include <set>
#include <cassert>
#include <cstddef>

namespace pgo
{
namespace BasicAlgorithms
{

// a heap that finds the Key with the smallest Value
// the input Key must be unique
// has the ability to find the Key inside the heap (therefore named Reversible)
template<class Value, class Key>
class ReversibleHeap
{
public:
  ReversibleHeap() {}
  virtual ~ReversibleHeap() {}

  typedef std::pair<Value, Key> ValueKeyPair;

  // add (value, key) into the heap, each key must be unique
  void push(const Value &value, const Key &key);

  void pop();

  // get the top (value, key) pair
  const ValueKeyPair &top() const { return *setQueue.begin(); }

  void erase(const Key &key);

  // add or update (value, key)
  void update(const Value &value, const Key &key);

  void clear();

  // get the value associated with key; return nullptr if not found
  const Value *find(const Key &key) const;

  size_t size() const { return setQueue.size(); }

protected:
  typedef std::set<ValueKeyPair> ValueKeySet;
  typedef typename ValueKeySet::iterator ValueKeySetIter;
  typedef std::map<Key, ValueKeySetIter> ReverseMap;
  typedef typename ReverseMap::iterator ReverseMapIter;

  ValueKeySet setQueue;
  ReverseMap reverseMap;
};

template<class Value, class Key>
void ReversibleHeap<Value, Key>::push(const Value &value, const Key &key)
{
  std::pair<ValueKeySetIter, bool> p = setQueue.insert(ValueKeyPair(value, key));
  assert(p.second == true);  // assert no same ValueKeyPair stored before
  ReverseMapIter it = reverseMap.find(key);
  assert(it == reverseMap.end());
  reverseMap.insert(std::pair<Key, ValueKeySetIter>(key, p.first));
}

template<class Value, class Key>
void ReversibleHeap<Value, Key>::update(const Value &value, const Key &key)
{
  ReverseMapIter it = reverseMap.find(key);
  if (it != reverseMap.end())  // this key has been added before
  {
    setQueue.erase(it->second);
    std::pair<ValueKeySetIter, bool> p = setQueue.insert(ValueKeyPair(value, key));
    assert(p.second == true);  // assert no same ValueKeyPair stored before
    it->second = p.first;      // update the ValueKeySetIter
  }
  else {
    std::pair<ValueKeySetIter, bool> p = setQueue.insert(ValueKeyPair(value, key));
    assert(p.second == true);  // assert no same ValueKeyPair stored before
    reverseMap.insert(std::pair<Key, ValueKeySetIter>(key, p.first));
  }
}

template<class Value, class Key>
void ReversibleHeap<Value, Key>::pop()
{
  ValueKeySetIter it = setQueue.begin();
  reverseMap.erase(it->second);
  setQueue.erase(it);
}

template<class Value, class Key>
void ReversibleHeap<Value, Key>::erase(const Key &key)
{
  ReverseMapIter it = reverseMap.find(key);
  if (it != reverseMap.end()) {
    setQueue.erase(it->second);
    reverseMap.erase(it);
  }
}

template<class Value, class Key>
void ReversibleHeap<Value, Key>::clear()
{
  reverseMap.clear();
  setQueue.clear();
}

template<class Value, class Key>
const Value *ReversibleHeap<Value, Key>::find(const Key &key) const
{
  typename std::map<Key, ValueKeySetIter>::const_iterator it = reverseMap.find(key);
  if (it != reverseMap.end())
    return &it->second->first;  // it->second is a ValueKeySetIter, it->second->first is the value pointed by the ValueKeySetIter
  return nullptr;
}

}  // namespace BasicAlgorithms
}  // namespace pgo