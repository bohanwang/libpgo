/*
  This code is based on code from the Geometric Tools library,
  which is licensed under a boost license.
  Such usage is permitted by the boost license; for details,
  please see the boost license below.
*/

// Geometric Tools, LLC
// Copyright (c) 1998-2014
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt

/*************************************************************************
 *                                                                       *
 * We release our improvements to the wildMagic code under our standard  *
 * Vega FEM license, as follows:                                         *
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "improvements to the wildMagic library" , Copyright (C) 2018 USC      *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code author: Yijing Li                                                *
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

#include "edgeKey.h"
#include "meshLinearAlgebra.h"

#include <algorithm>
#include <ostream>

namespace pgo
{
namespace Mesh
{
// unoriented triangle key based on vtx indices
// v[0], v[1] and v[2] will be sorted to ensure v[0] <= v[1] <= v[2]
struct UTriKey
{
  inline UTriKey();  // creates an invalid key with v = {-1,-1,-1}
  inline UTriKey(int v0, int v1, int v2);
  inline UTriKey(const int v[3]):
    UTriKey(v[0], v[1], v[2]) {}
  inline UTriKey(const Vec3i &v):
    UTriKey(v[0], v[1], v[2]) {}

  inline bool operator<(const UTriKey &o) const { return v < o.v; }
  inline bool operator==(const UTriKey &o) const { return v == o.v; }
  inline bool operator!=(const UTriKey &o) const { return v != o.v; }
  const int &operator[](int index) const { return v[index]; }
  const int *indices() const { return &v[0]; }
  const int *data() const { return v.data(); }

  // return the unordered edge opposite to v[i], i:[0,3)
  inline UEdgeKey uEdgeKey(int i) const { return UEdgeKey(v[triEdgeIndex[i][0]], v[triEdgeIndex[i][1]]); }
  // given the global vtx index, return its first index in v [0,2]
  // else return -1
  inline int getInvertedIndex(const int globalVtxIndex) const { return pgo::getInvertedIndex(v, globalVtxIndex); }

  inline bool shareUEdge(const UTriKey &nbr) const;
  inline UEdgeKey getSharedUEdge(const UTriKey &nbr) const;  // return the first shared UEdge; return (-1,-1) if no shared UEdge

  inline bool hasUEdge(const UEdgeKey &edge) const;
  inline bool hasIndex(int vtxIndex) const { return v[0] == vtxIndex || v[1] == vtxIndex || v[2] == vtxIndex; }

  // return -1 if edge not found
  inline int getVertexOppositeEdge(const UEdgeKey &edge) const;

  // if all indices are >= 0
  inline bool allVerticesValid() const { return v[0] >= 0; }
  // if this triangle is valid: v[i] >=0 && v[i] != v[j]
  inline bool isValidTriangle() const;

  // index: { { 1, 2 }, { 0, 2 }, { 0, 1 } }
  static const int triEdgeIndex[3][2];

protected:
  Vec3i v;
};

inline std::ostream &operator<<(std::ostream &s, const UTriKey &v);

// oriented triangle key that remember its orientation
// v[0],v[1] and v[2] will be reordered to ensure v[0] <= v[1] and v[0] <= v[2] and (v[0], v[1], v[2]) remain the same orientation as the input
struct OTriKey
{
  inline OTriKey();  // creates an invalid key with v = {-1,-1,-1}
  inline OTriKey(int v0, int v1, int v2);
  inline OTriKey(const int v[3]):
    OTriKey(v[0], v[1], v[2]) {}
  inline OTriKey(const Vec3i &v):
    OTriKey(v[0], v[1], v[2]) {}

  inline bool operator<(const OTriKey &o) const { return v < o.v; }
  inline bool operator==(const OTriKey &o) const { return v == o.v; }
  inline bool operator!=(const OTriKey &o) const { return v != o.v; }
  const int &operator[](int index) const { return v[index]; }
  const int *indices() const { return &v[0]; }
  const int *data() const { return v.data(); }

  inline UTriKey uTriKey() const { return UTriKey(&v[0]); }

  // return the ordered/unordered edge opposite to v[i], i:[0,3)
  inline UEdgeKey uEdgeKey(int i) const { return UEdgeKey(v[triEdgeIndex[i][0]], v[triEdgeIndex[i][1]]); }
  inline OEdgeKey oEdgeKey(int i) const { return OEdgeKey(v[triEdgeIndex[i][0]], v[triEdgeIndex[i][1]]); }

  // given the global vtx index, return the index in v [0,2]; otherwise return -1
  inline int getInvertedIndex(const int globalVtxIndex) const { return pgo::getInvertedIndex(v, globalVtxIndex); }

  inline bool shareUEdge(const OTriKey &nbr) const;
  inline bool shareOEdge(const OTriKey &nbr) const;
  inline bool shareOppositeOEdge(const OTriKey &nbr) const;
  inline OEdgeKey getSharedOppositeOEdge(const OTriKey &nbr) const;  // return first owned OEdge whose reverse is owned by nbr; return (-1,-1) otherwise
  inline UEdgeKey getSharedUEdge(const OTriKey &nbr) const;          // return the first shared UEdge; return (-1,-1) if no shared UEdge

  inline int getInvertedOEdgeIndex(const OEdgeKey &edge) const;

  bool hasOEdge(const OEdgeKey &edge) const;
  bool hasUEdge(const UEdgeKey &edge) const;

  inline void reverse() { std::swap(v[1], v[2]); }                                // reverse the orientation
  inline OTriKey getReversedTriKey() const { return OTriKey(v[0], v[2], v[1]); }  // return reversed tri key

  // if all indices are >= 0
  inline bool allVerticesValid() const { return v[0] >= 0; }
  // if this triangle is valid: v[i] >=0 && v[i] != v[j]
  inline bool isValidTriangle() const;

  // permute v0, v1, v2 and store into r0, r1, r2 so that they share the same orientation but r0 = min(v0,v1,v2)
  static void permute(int v0, int v1, int v2, int &r0, int &r1, int &r2);

  // index: { { 1, 2 }, { 2, 0 }, { 0, 1 } }
  static const int triEdgeIndex[3][2];

protected:
  Vec3i v;
};

inline std::ostream &operator<<(std::ostream &s, const OTriKey &v);

///////////////////////////////////////////////////////////////////////////////
//                             IMPLEMENTATION                                //
///////////////////////////////////////////////////////////////////////////////

inline UTriKey::UTriKey()
{
  v[0] = v[1] = v[2] = -1;
}

inline UTriKey::UTriKey(int v0, int v1, int v2)
{
  v[0] = v0;
  v[1] = v1;
  v[2] = v2;
  std::sort(v.begin(), v.end());
}

// find whether two elemetns are shared among the two tuples: (v0, v1, v2) and (n0, n1, n2)
// we utilize the fact that vertices are pre-sorted in UTriKey
inline bool UTriKey::shareUEdge(const UTriKey &nbr) const
{
  for (int i = 0; i < 3; i++) {
    UEdgeKey key = uEdgeKey(i);
    for (int j = 0; j < 3; j++)
      if (key == nbr.uEdgeKey(j))
        return true;
  }
  return false;

  // const int * n = &nbr.v[0];

  // if (v[0] == n[0])
  // {
  //   if (v[1] == n[1] || v[1] == n[2] || v[2] == n[1] || v[2] == n[2])
  //     return true;
  //   return false;
  // }
  // else if (v[0] == n[1])
  // {
  //   // if enthering this branch (v[0] == n[1]), then n[0] <= v[1]. we can prove v[1] != n[0], otherwise v[0]==v[1]==n[0]==n[1], should enter the first branch already
  //   // same case for v[2] != n[0]
  //   if (v[1] == n[2] || v[2] == n[2])
  //     return true;
  //   return false;
  // }
  // // if v[0] == n[2], then n[1] <= v[1]
  // // if n[1] < v[1], these means (v1, v2) will never have one same element as (n0, n1), cannot share UEdge
  // // else, n[1] == v[1], then v0 == n1, should go to the second branch already
  // // so we don't consider the case where v0 == v2
  // if (v[1] == n[1] && v[2] == n[2])
  //   return true;
  // return false;
}

inline UEdgeKey UTriKey::getSharedUEdge(const UTriKey &nbr) const
{
  for (int i = 0; i < 3; i++) {
    UEdgeKey key = uEdgeKey(i);
    for (int j = 0; j < 3; j++)
      if (key == nbr.uEdgeKey(j))
        return key;
  }
  return UEdgeKey();  // return a default invalid UEdgeKey
}

inline bool UTriKey::isValidTriangle() const
{
  return v[0] >= 0 && v[0] != v[1] && v[1] != v[2];  // because v[0] <= v[1] <= v[2]
}

inline bool UTriKey::hasUEdge(const UEdgeKey &edge) const
{
  for (int i = 0; i < 3; i++)
    if (edge == uEdgeKey(i))
      return true;
  return false;
}

inline int UTriKey::getVertexOppositeEdge(const UEdgeKey &edge) const
{
  for (int i = 0; i < 3; i++)
    if (edge == uEdgeKey(i))
      return v[i];
  return -1;
}

////////////////////////////////////////////////////////////////////
//                         OTriKey
////////////////////////////////////////////////////////////////////

inline OTriKey::OTriKey()
{
  v[0] = v[1] = v[2] = -1;
}

inline OTriKey::OTriKey(int v0, int v1, int v2)
{
  permute(v0, v1, v2, v[0], v[1], v[2]);
}

// find whether two elemetns are shared among the two tuples: (v0, v1, v2) and (n0, n1, n2)
// we utilize the fact that v0 <= min(v1, v2) and n0 <= min(n1, n2)
inline bool OTriKey::shareUEdge(const OTriKey &nbr) const
{
  for (int i = 0; i < 3; i++) {
    UEdgeKey key = uEdgeKey(i);
    for (int j = 0; j < 3; j++)
      if (key == nbr.uEdgeKey(j))
        return true;
  }
  return false;

  // const int * n = &nbr.v[0];
  // if (v[0] == n[0])
  // {
  //   if (v[1] == n[1] || v[1] == n[2] || v[2] == n[1] || v[2] == n[2])
  //     return true;
  //   return false;
  // }
  // else if (v[0] == n[1])
  // {
  //   // if enthering this branch (v[0] == n[1]), then v[1] != n[0], otherwise v[0]==v[1]==n[0]==n[1], should enter the first branch already
  //   // same case for v[2] != n[0]
  //   if (v[1] == n[2] || v[2] == n[2])
  //     return true;
  //   return false;
  // }
  // else if (v[0] == n[2])
  // {
  //   // if enthering this branch (v[0] == n[2]), then v[1] != n[0], otherwise v[0]==v[1]==n[0]==n[1], should enter the first branch already
  //   // same case for v[2] != n[0]
  //   if (v[1] == n[1] || v[2] == n[1])
  //     return true;
  //   return false;
  // }

  // if (v[1] == n[1] && v[2] == n[2])
  //   return true;
  // return false;
}

// find whether two TriKeys share the same OEdge
// we utilize the fact that v0 <= min(v1, v2) and n0 <= min(n1, n2)
// only nine cases
// case 1, (v1, v2) == (n1, n2)
// case 2, (v1, v2) == (n2, n0)
// case 3, (v1, v2) == (n0, n1)
// case 4, (v2, v0) == (n1, n2)
// case 5, (v2, v0) == (n2, n0)
// case 6, (v2, v0) == (n0, n1)
// case 7, (v0, v1) == (n1, n2)
// case 8, (v0, v1) == (n2, n0)
// case 9, (v0, v1) == (n0, n1)
inline bool OTriKey::shareOEdge(const OTriKey &nbr) const
{
  // const int * n = nbr.v;
  // if (v[0] == n[0]) // case 5, 9
  // {
  //   if (v[1] == n[1] || v[1] == n[2])
  //     return true;
  // }
  // else if (v[0] == n[1]) // case 6, 7
  // {
  //   // if enthering this branch (v[0] == n[1]), then v[2] != n[0], otherwise v[0]==v[2]==n[0]==n[1], should enter the first branch already
  //   if (v[1] == n[2])
  //     return true;
  // }
  // else if (v[0] == n[2]) // case 4, 8
  // {
  //   // if enthering this branch (v[0] == n[2]), then v[1] != n[0], otherwise v[0]==v[1]==n[0]==n[1], should enter the first branch already
  //   if (v[2] == n[1])
  //     return true;
  // }
  // // check case 1, 2, 3
  // if (v[1] == n[1] && v[2] == n[2])
  //   return true;
  // if (v[1] == n[2] && v[2] == n[0])
  //   return true;
  // if (v[1] == n[0] && v[2] == n[1])
  //   return true;

  for (int i = 0; i < 3; i++) {
    OEdgeKey key = oEdgeKey(i);
    for (int j = 0; j < 3; j++)
      if (key == nbr.oEdgeKey(j))
        return true;
  }

  return false;
}

inline bool OTriKey::shareOppositeOEdge(const OTriKey &nbr) const
{
  for (int i = 0; i < 3; i++) {
    OEdgeKey key = oEdgeKey(i);
    key.reverse();
    for (int j = 0; j < 3; j++)
      if (key == nbr.oEdgeKey(j))
        return true;
  }
  return false;
}

inline OEdgeKey OTriKey::getSharedOppositeOEdge(const OTriKey &nbr) const
{
  for (int i = 0; i < 3; i++) {
    OEdgeKey key = oEdgeKey(i);
    key.reverse();
    for (int j = 0; j < 3; j++)
      if (key == nbr.oEdgeKey(j))
        return key.getReversedEdgeKey();
  }
  return OEdgeKey();
}

inline UEdgeKey OTriKey::getSharedUEdge(const OTriKey &nbr) const
{
  for (int i = 0; i < 3; i++) {
    UEdgeKey key = uEdgeKey(i);
    for (int j = 0; j < 3; j++)
      if (key == nbr.uEdgeKey(j))
        return key;
  }
  return UEdgeKey();  // return a default invalid UEdgeKey
}

inline bool OTriKey::hasOEdge(const OEdgeKey &edge) const
{
  for (int i = 0; i < 3; i++)
    if (edge == oEdgeKey(i))
      return true;
  return false;
}

inline bool OTriKey::hasUEdge(const UEdgeKey &edge) const
{
  for (int i = 0; i < 3; i++)
    if (edge == uEdgeKey(i))
      return true;
  return false;
}

inline int OTriKey::getInvertedOEdgeIndex(const OEdgeKey &edge) const
{
  for (int i = 0; i < 3; i++)
    if (edge == oEdgeKey(i))
      return i;
  return -1;
}

inline bool OTriKey::isValidTriangle() const
{
  return v[0] >= 0 && v[0] != v[1] && v[1] != v[2] && v[0] != v[2];  // because v[0] <= v[1], v[0] <= v[2]
}

inline std::ostream &operator<<(std::ostream &s, const UTriKey &v)
{
  return s << '(' << v[0] << ' ' << v[1] << ' ' << v[2] << ')';
}

inline std::ostream &operator<<(std::ostream &s, const OTriKey &v)
{
  return s << '(' << v[0] << ' ' << v[1] << ' ' << v[2] << ')';
}
}  // namespace Mesh
}  // namespace pgo

namespace std
{
template<>
struct hash<pgo::Mesh::UTriKey>
{
  size_t operator()(const pgo::Mesh::UTriKey &k) const
  {
    static_assert(sizeof(int) * 2 == sizeof(uint64_t), "uint64_t is not twice the same size as int");
    static_assert(sizeof(uint64_t) == sizeof(std::size_t), "uint64_t is not the same size as std::size_t");
    
    std::size_t v = (((std::size_t)k[0]) + ((std::size_t)(k[1]) << (sizeof(int) * 8)));
    pgo::EigenSupport::hashCombine(v, k[2]);
    return v;
  }
};

template<>
struct hash<pgo::Mesh::OTriKey>
{
  size_t operator()(const pgo::Mesh::OTriKey &k) const
  {
    static_assert(sizeof(int) * 2 == sizeof(uint64_t), "uint64_t is not twice the same size as int");
    static_assert(sizeof(uint64_t) == sizeof(std::size_t), "uint64_t is not the same size as std::size_t");
    
    std::size_t v = (((std::size_t)k[0]) + ((std::size_t)(k[1]) << (sizeof(int) * 8)));
    pgo::EigenSupport::hashCombine(v, k[2]);
    return v;
  }
};
}  // namespace std
