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

#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <unordered_set>

namespace pgo
{
namespace BasicAlgorithms
{

/////////////////////////////////////////////////////
//            Check element in set/map             //
/////////////////////////////////////////////////////

template<class T, class C, class A>
bool setFind(const std::set<T, C, A> &s, const T &x)
{
  return s.find(x) != s.end();
}
template<class T, class C, class A>
bool setNotFind(const std::set<T, C, A> &s, const T &x)
{
  return s.find(x) == s.end();
}

template<class T, class U, class C, class A>
bool mapFind(const std::map<T, U, C, A> &s, const T &x)
{
  return s.find(x) != s.end();
}
template<class T, class U, class C, class A>
bool mapNotFind(const std::map<T, U, C, A> &s, const T &x)
{
  return s.find(x) == s.end();
}

template<class T, class H, class P, class A>
bool usetFind(const std::unordered_set<T, H, P, A> &s, const T &x)
{
  return s.find(x) != s.end();
}

template<class T, class H, class P, class A>
bool usetNotFind(const std::unordered_set<T, H, P, A> &s, const T &x)
{
  return s.find(x) == s.end();
}

/////////////////////////////////////////////////////
//            insert Data in Container             //
/////////////////////////////////////////////////////

template<class T, class A, class Container>
void vectorInsertRangeBack(std::vector<T, A> &v, const Container &c)
{
  v.insert(v.end(), c.begin(), c.end());
}

template<class K, class V, class C, class A>
typename std::map<K, V, C, A>::iterator assignOrAddToMap(std::map<K, V, C, A> &m, const K &key, const V &value)
{
  auto it = m.lower_bound(key);
  if (it != m.end() && it->first == key)
    it->second += value;
  else
    it = m.emplace_hint(it, key, value);
  return it;
}

/////////////////////////////////////////////////////
//                Memory operations                //
/////////////////////////////////////////////////////

template<class T, class A>
inline void memset(std::vector<T, A> &v, unsigned char value = 0)
{
  memset(v.data(), value, sizeof(T) * v.size());
}

template<class T, class A>
inline void memcpy(std::vector<T, A> &v, const T *src)
{
  memcpy(v.data(), src, sizeof(T) * v.size());
}

// free any memory allocated in the vector
template<class T, class A>
inline void freeMemory(std::vector<T, A> &v)
{
  std::vector<T> tmp;
  v.swap(tmp);
}

template<class T, class A>
inline bool inVectorRange(int index, const std::vector<T, A> &vec)
{
  return index >= 0 && index < (int)vec.size();
}

/////////////////////////////////////////////////////
//               Mapping operations                //
/////////////////////////////////////////////////////

// create a reverse mapping from each value of vec to its index
// assuming the values in vec are unique
template<class T, class A>
std::map<int, int> getReverseMapping(const std::vector<T, A> &vec)
{
  std::map<int, int> ret;
  for (size_t i = 0; i < vec.size(); i++)
    ret[vec[i]] = i;
  return ret;
}

// return a vector containing unique elements
// inpID2unqID: optionally return the mapping: input ID -> output unique ID
// unqID2inpID: optionally return the mapping: output unique ID -> input ID
template<class T, class A>
std::vector<T, A> findUniqueElements(const std::vector<T, A> vec, std::vector<int> *inpID2unqID = nullptr, std::vector<int> *unqID2inpID = nullptr);

/////////////////////////////////////////////////////
//           Euclidean norm / distance             //
/////////////////////////////////////////////////////

// the squared Euclidean norm of the vector
template<class T>
T squaredEuclideanNorm(const std::vector<T> &v);
template<class T>
T squaredEuclideanNorm(size_t r, const T *v);

// the squared Euclidean distance between the two vectors of the same size
template<class T>
T squaredEuclideanDistance(const std::vector<T> &v1, const std::vector<T> &v2);
template<class T>
T squaredEuclideanDistance(size_t r, const T *v1, const T *v2);

// compute Euclidean norm / distance of double arrays
template<class T>
inline double EuclideanNorm(const std::vector<T> &v)
{
  return sqrt(squaredEuclideanNorm(v));
}

template<class T>
inline double EuclideanNorm(size_t r, const T *v)
{
  return sqrt(squaredEuclideanNorm(r, v));
}

template<class T>
inline double EuclideanDistance(const std::vector<T> &v1, const std::vector<T> &v2)
{
  return sqrt(squaredEuclideanDistance(v1, v2));
}

template<class T>
inline double EuclideanDistance(size_t r, const T *v1, const T *v2)
{
  return sqrt(squaredEuclideanDistance(r, v1, v2));
}

/////////////////////////////////////////////////////
//                  Set operation                  //
/////////////////////////////////////////////////////

// return whether two sets intersect
template<class T, class C, class A>
bool intersect(const std::set<T, C, A> &s1, const std::set<T, C, A> &s2);

// compute the intersection of two vectors and store the result in out (previous data in out is removed)
template<class T>
void intersect(const std::vector<T> &v1, const std::vector<T> &v2, std::vector<T> &out);

// return whether two sorted vectors intersect
template<class T>
bool intersectSorted(const std::vector<T> &v1, const std::vector<T> &v2);

template<class T, class C, class A>
int getIntersectionSize(const std::set<T, C, A> &s1, const std::set<T, C, A> &s2);

// return true if A is a subset of B
template<class T, class C, class A>
bool isSubset(const std::set<T, C, A> &setA, const std::set<T, C, A> &setB);

// compute s1 / s2
template<class T, class C, class A>
std::set<T, C, A> setMinus(const std::set<T, C, A> &s1, const std::set<T, C, A> &s2);

/////////////////////////////////////////////////////
//                 Print elements                  //
/////////////////////////////////////////////////////

// print (a->x, b->y, c->z, ...) to ostream
template<class T, class U>
void printMap(std::ostream &os, const std::map<T, U> &s);

/////////////////////////////////////////////////////
//         SAVE & LOAD IMPLEMENTATIONS             //
/////////////////////////////////////////////////////
// for use in saving/loading a class which includes vector as member var.
// Same save/load file format as used in vega mayaPlugin
template<class T>
bool saveToAscii(const std::vector<T> &v, std::ostream &out);

template<class T>
bool loadFromAscii(std::vector<T> &v, std::istream &in);

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

// Below are template implementation

template<class T, class C, class A>
bool intersect(const std::set<T, C, A> &s1, const std::set<T, C, A> &s2)
{
  if (s1.size() == 0 || s2.size() == 0)
    return false;
  const std::set<T> *ps1 = &s1, *ps2 = &s2;
  if (s1.size() > s2.size()) {
    ps1 = &s2;
    ps2 = &s1;
  }
  for (typename std::set<T, C, A>::const_iterator it = ps1->begin(); it != ps1->end(); it++)
    if (ps2->find(*it) != ps2->end())
      return true;
  return false;
}

template<class T>
void intersect(const std::vector<T> &v1, const std::vector<T> &v2, std::vector<T> &out)
{
  out.clear();
  if (v1.size() == 0 || v2.size() == 0)
    return;
  const std::vector<T> *pv1 = &v1, *pv2 = &v2;
  if (v1.size() > v2.size()) {
    pv1 = &v2;  // pv1 points to the smaller set
    pv2 = &v1;
  }
  std::set<T> s2;
  s2.insert(pv2->begin(), pv2->end());
  for (size_t i = 0; i < pv1->size(); i++) {
    const T &x = (*pv1)[i];
    if (s2.find(x) != s2.end())
      out.push_back(x);
  }
}

template<class T>
bool intersectSorted(const std::vector<T> &v1, const std::vector<T> &v2)
{
  if (v1.size() == 0 || v2.size() == 0)
    return false;
  const std::vector<T> *pv1 = &v1, *pv2 = &v2;
  if (v1.size() > v2.size()) {
    pv1 = &v2;  // pv1 points to the smaller set
    pv2 = &v1;
  }

  for (size_t i = 0; i < pv1->size(); i++)
    if (binary_search(pv2->begin(), pv2->end(), (*pv1)[i]))
      return true;

  return false;
}

template<class T, class C, class A>
int getIntersectionSize(const std::set<T, C, A> &s1, const std::set<T, C, A> &s2)
{
  if (s1.size() == 0 || s2.size() == 0)
    return 0;
  const std::set<T, C, A> *ps1 = &s1, *ps2 = &s2;
  if (s1.size() > s2.size())  // make ps1->size() <= ps2->size()
  {
    ps1 = &s2;
    ps2 = &s1;
  }
  int count = 0;
  for (typename std::set<T, C, A>::const_iterator it = ps1->begin(); it != ps1->end(); it++)
    if (ps2->find(*it) != ps2->end())
      count++;
  return count;
}

template<class T, class C, class A>
bool isSubset(const std::set<T, C, A> &setA, const std::set<T, C, A> &setB)
{
  if (setA.size() > setB.size())
    return false;
  for (typename std::set<T, C, A>::const_iterator it = setA.begin(); it != setA.end(); it++)
    if (setB.find(*it) == setB.end())
      return false;
  return true;
}

template<class T, class C, class A>
std::set<T, C, A> setMinus(const std::set<T, C, A> &s1, const std::set<T, C, A> &s2)
{
  std::set<T, C, A> ret;
  for (const auto &e : s1) {
    if (setNotFind(s2, e))
      ret.insert(e);
  }
  return ret;
}

template<class T>
T squaredEuclideanNorm(const std::vector<T> &v)
{
  T sum = 0;
  for (size_t i = 0; i < v.size(); i++)
    sum += v[i] * v[i];
  return sum;
}

template<class T>
T squaredEuclideanNorm(size_t r, const T *v)
{
  T sum = 0;
  for (size_t i = 0; i < r; i++)
    sum += v[i] * v[i];
  return sum;
}

template<class T>
T squaredEuclideanDistance(const std::vector<T> &v0, const std::vector<T> &v1)
{
  T sum = 0;
  for (size_t i = 0; i < v0.size(); i++) {
    T value = (v0[i] - v1[i]);
    sum += value * value;
  }
  return sum;
}

template<class T>
T squaredEuclideanDistance(size_t r, const T *v1, const T *v2)
{
  T sum = 0;
  for (size_t i = 0; i < r; i++) {
    T value = (v1[i] - v2[i]);
    sum += value * value;
  }
  return sum;
}

// helper template class to stream a range with begin(), end() implemented
template<class Range>
class RangeStreamer
{
public:
  RangeStreamer(const Range &input):
    range(input) {}

  friend std::ostream &operator<<(std::ostream &os, const RangeStreamer &r)
  {
    os << "{";
    bool first = true;
    for (const auto &v : r.range) {
      if (first)
        first = false;
      else
        os << ", ";
      os << v;
    }
    os << "}";
    return os;
  }

protected:
  const Range &range;
};

// use this function to create RangeStreamer so that you can use << to output a range (e.g. a container)
// usage:
// vector<int> a = {1,2,3};
// cout << streamRange(a) << endl; // print {1, 2, 3}
template<class Range>
RangeStreamer<Range> streamRange(const Range &input)
{
  return RangeStreamer<Range>(input);
}

template<class T, class U>
void printMap(std::ostream &os, const std::map<T, U> &s)
{
  os << "(";
  for (typename std::map<T, U>::const_iterator it = s.begin(); it != s.end(); it++) {
    if (it != s.begin())
      os << ", ";
    os << "<" << it->first << " -> " << it->second << ">";
  }
  os << ")";
}

template<class T, class A>
std::vector<T, A> findUniqueElements(const std::vector<T, A> vec, std::vector<int> *inpID2unqID, std::vector<int> *unqID2inpID)
{
  std::vector<T, A> ret;
  if (inpID2unqID)
    inpID2unqID->resize(vec.size());
  if (unqID2inpID)
    unqID2inpID->clear();

  if (vec.size() == 0)
    return ret;

  std::vector<int> buffer(vec.size());
  std::iota(buffer.begin(), buffer.end(), 0);

  auto cmp = [&](const int &p0, const int &p1) {
    return vec[p0] < vec[p1];
  };
  std::sort(buffer.begin(), buffer.end(), cmp);

  ret.push_back(vec[buffer[0]]);
  if (unqID2inpID)
    unqID2inpID->push_back(buffer[0]);
  if (inpID2unqID)
    inpID2unqID->at(buffer[0]) = 0;
  for (std::size_t i = 1; i < vec.size(); i++) {
    int vecID = buffer[i];
    if ((ret.back() != vec[vecID])) {
      ret.push_back(vec[vecID]);
      if (unqID2inpID)
        unqID2inpID->push_back(vecID);
    }
    if (inpID2unqID)
      inpID2unqID->at(vecID) = ret.size() - 1;
  }

  return ret;
}

}  // namespace BasicAlgorithms
}  // namespace pgo