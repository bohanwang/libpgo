/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "mesh" library , Copyright (C) 2018 USC                               *
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

#include "meshLinearAlgebra.h"
#include "edgeKey.h"

#include <algorithm>
#include <ostream>
#include <cstring>

namespace pgo
{
namespace Mesh
{

// unoriented rectangle key based on vtx indices
// v[0], v[1], v[2] and v[3] will be sorted to ensure v[0] <= v[1] <= v[2] <= v[3]
struct URectKey
{
  inline URectKey(const int v[4]);
  inline URectKey(int v0, int v1, int v2, int v3);
  inline URectKey();  // creates an invalid key with v = {-1,-1,-1, -1}
  inline bool operator<(const URectKey &o) const { return v < o.v; }
  inline bool operator==(const URectKey &o) const { return v == o.v; }
  inline bool operator!=(const URectKey &o) const { return v != o.v; }
  const int &operator[](int index) const { return v[index]; }
  const int *indices() const { return &v[0]; }
  // return the unordered edge opposite to v[i], i: [0,4)
  inline UEdgeKey uEdgeKey(int i) const { return UEdgeKey(v[rectEdgeIndex[i][0]], v[rectEdgeIndex[i][1]]); }
  // given the global vtx index, return its first index in v [0,2]
  // else return -1
  inline int getInvertedIndex(const int globalVtxIndex) const { return pgo::getInvertedIndex(v, globalVtxIndex); }

  inline bool shareUEdge(const URectKey &nbr) const;

  static const int rectEdgeIndex[4][2];

protected:
  Vec4i v;
};

inline std::ostream &operator<<(std::ostream &s, const URectKey &v);

// oriented rectangle key that remember its orientation
// v[0],v[1] and v[2] will be reordered to ensure v[0] <= v[1] and v[0] <= v[2] and (v[0], v[1], v[2]) remain the same orientation as the input
struct ORectKey
{
  inline ORectKey(const int v[4]);
  inline ORectKey(int v0, int v1, int v2, int v3);
  inline ORectKey();  // creates an invalid key with v = {-1,-1,-1,-1}
  inline bool operator<(const ORectKey &o) const { return v < o.v; }
  inline bool operator==(const ORectKey &o) const { return v == o.v; }
  inline bool operator!=(const ORectKey &o) const { return v != o.v; }
  const int &operator[](int i) const { return v[i]; }
  const int *indices() const { return &v[0]; }

  inline URectKey uRectKey() const { return URectKey(&v[0]); }

  // return the ordered/unordered edge opposite to v[ind]
  inline UEdgeKey uEdgeKey(int i) const { return UEdgeKey(v[rectEdgeIndex[i][0]], v[rectEdgeIndex[i][1]]); }
  inline OEdgeKey oEdgeKey(int i) const { return OEdgeKey(v[rectEdgeIndex[i][0]], v[rectEdgeIndex[i][1]]); }

  // given the global vtx index, return the index in v [0,2]; otherwise return -1
  inline int getInvertedIndex(const int globalVtxIndex) const { return pgo::getInvertedIndex(v, globalVtxIndex); }

  inline bool shareUEdge(const ORectKey &nbr) const;
  inline bool shareOEdge(const ORectKey &nbr) const;
  inline bool shareOppositeOEdge(const ORectKey &nbr) const;
  inline OEdgeKey getSharedOppositeOEdge(const ORectKey &nbr) const;  // return first owned OEdge whose reverse is owned by nbr; return (-1,-1) otherwise
  inline UEdgeKey getSharedUEdge(const ORectKey &nbr) const;          // return the first shared UEdge; return (-1,-1) if no shared UEdge

  inline int getInvertedOEdgeIndex(const OEdgeKey &edge) const;

  bool hasOEdge(const OEdgeKey &edge) const;
  bool hasUEdge(const UEdgeKey &edge) const;

  inline void reverse() { std::swap(v[1], v[2]); }                                         // reverse the orientation
  inline ORectKey getReversedRectKey() const { return ORectKey(v[0], v[3], v[2], v[1]); }  // return reversed rectKey

  // rotate v0-v3 so that they share the same cyclic order but v0 = min(v0,v1,v2,v3)
  inline static void rotate(Vec4i &v);

  static const int rectEdgeIndex[4][2];

protected:
  Vec4i v;
};

inline std::ostream &operator<<(std::ostream &s, const ORectKey &v);

///////////////////////////////////////////////////////////////////////////////
//                             IMPLEMENTATION                                //
///////////////////////////////////////////////////////////////////////////////

inline URectKey::URectKey(int v0, int v1, int v2, int v3)
{
  v[0] = v0;
  v[1] = v1;
  v[2] = v2;
  v[3] = v3;
  std::sort(&v[0], &v[0] + 4);
}

inline URectKey::URectKey(const int vtx[4])
{
  memcpy(&v[0], vtx, sizeof(int) * 4);
  std::sort(&v[0], &v[0] + 4);
}

inline URectKey::URectKey()
{
  v[0] = v[1] = v[2] = v[3] = -1;
}

inline ORectKey::ORectKey(int v0, int v1, int v2, int v3)
{
  v[0] = v0;
  v[1] = v1;
  v[2] = v2;
  v[3] = v3;
  rotate(v);
}

inline ORectKey::ORectKey(const int vtx[4])
{
  memcpy(&v[0], vtx, sizeof(int) * 4);
  rotate(v);
}

inline ORectKey::ORectKey()
{
  v[0] = v[1] = v[2] = v[3] = -1;
}

inline void ORectKey::rotate(Vec4i &v)
{
  int *vmin = std::min_element(&v[0], &v[0] + 4);
  int minIdx = static_cast<int>(vmin - &v[0]);
  pgo::rotate(v, minIdx);
}

inline bool URectKey::shareUEdge(const URectKey &nbr) const
{
  for (int i = 0; i < 4; i++) {
    UEdgeKey key = uEdgeKey(i);
    for (int j = 0; j < 4; j++)
      if (key == nbr.uEdgeKey(j))
        return true;
  }
  return false;
}

inline bool ORectKey::shareUEdge(const ORectKey &nbr) const
{
  for (int i = 0; i < 4; i++) {
    UEdgeKey key = uEdgeKey(i);
    for (int j = 0; j < 4; j++)
      if (key == nbr.uEdgeKey(j))
        return true;
  }
  return false;
}

inline bool ORectKey::shareOEdge(const ORectKey &nbr) const
{
  for (int i = 0; i < 4; i++) {
    OEdgeKey key = oEdgeKey(i);
    for (int j = 0; j < 4; j++)
      if (key == nbr.oEdgeKey(j))
        return true;
  }
  return false;
}

inline bool ORectKey::shareOppositeOEdge(const ORectKey &nbr) const
{
  for (int i = 0; i < 4; i++) {
    OEdgeKey key = oEdgeKey(i);
    key.reverse();
    for (int j = 0; j < 4; j++)
      if (key == nbr.oEdgeKey(j))
        return true;
  }
  return false;
}

inline OEdgeKey ORectKey::getSharedOppositeOEdge(const ORectKey &nbr) const
{
  for (int i = 0; i < 4; i++) {
    OEdgeKey key = oEdgeKey(i);
    key.reverse();
    for (int j = 0; j < 4; j++)
      if (key == nbr.oEdgeKey(j))
        return key.getReversedEdgeKey();
  }
  return OEdgeKey();
}

inline UEdgeKey ORectKey::getSharedUEdge(const ORectKey &nbr) const
{
  for (int i = 0; i < 4; i++) {
    UEdgeKey key = uEdgeKey(i);
    for (int j = 0; j < 4; j++)
      if (key == nbr.uEdgeKey(j))
        return key;
  }
  return UEdgeKey();  // return a default invalid UEdgeKey
}

inline bool ORectKey::hasOEdge(const OEdgeKey &edge) const
{
  for (int i = 0; i < 4; i++)
    if (edge == oEdgeKey(i))
      return true;
  return false;
}

inline bool ORectKey::hasUEdge(const UEdgeKey &edge) const
{
  for (int i = 0; i < 4; i++)
    if (edge == uEdgeKey(i))
      return true;
  return false;
}

inline int ORectKey::getInvertedOEdgeIndex(const OEdgeKey &edge) const
{
  for (int i = 0; i < 4; i++)
    if (edge == oEdgeKey(i))
      return i;
  return -1;
}

inline std::ostream &operator<<(std::ostream &s, const URectKey &v)
{
  return s << '(' << v[0] << ' ' << v[1] << ' ' << v[2] << ' ' << v[3] << ')';
}

inline std::ostream &operator<<(std::ostream &s, const ORectKey &v)
{
  return s << '(' << v[0] << ' ' << v[1] << ' ' << v[2] << ' ' << v[3] << ')';
}

}  // namespace Mesh
}  // namespace pgo

namespace std
{
template<>
struct hash<pgo::Mesh::URectKey>
{
  size_t operator()(const pgo::Mesh::URectKey &k) const
  {
    static_assert(sizeof(int) * 2 == sizeof(uint64_t), "uint64_t is not twice the same size as int");
    static_assert(sizeof(uint64_t) == sizeof(std::size_t), "uint64_t is not the same size as std::size_t");

    std::size_t v = (((std::size_t)k[0]) + ((std::size_t)(k[1]) << (sizeof(int) * 8)));
    std::size_t v2 = (((std::size_t)k[2]) + ((std::size_t)(k[3]) << (sizeof(int) * 8)));
    pgo::EigenSupport::hashCombine(v, v2);

    return v;
  }
};

template<>
struct hash<pgo::Mesh::ORectKey>
{
  size_t operator()(const pgo::Mesh::ORectKey &k) const
  {
    static_assert(sizeof(int) * 2 == sizeof(uint64_t), "uint64_t is not twice the same size as int");
    static_assert(sizeof(uint64_t) == sizeof(std::size_t), "uint64_t is not the same size as std::size_t");
    
    std::size_t v = (((std::size_t)k[0]) + ((std::size_t)(k[1]) << (sizeof(int) * 8)));
    std::size_t v2 = (((std::size_t)k[2]) + ((std::size_t)(k[3]) << (sizeof(int) * 8)));
    pgo::EigenSupport::hashCombine(v, v2);

    return v;
  }
};
}  // namespace std
