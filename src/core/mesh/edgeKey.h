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

#include <ostream>
#include <functional>
#include <utility>
#include <cstdint>

namespace pgo
{
namespace Mesh
{
// unoriented edge
struct UEdgeKey
{
  inline UEdgeKey(int v0, int v1);
  inline UEdgeKey(const std::pair<int, int> &p):
    UEdgeKey(p.first, p.second) {}
  inline UEdgeKey();  // creates an invalid key with v = {-1,-1}
  const int &operator[](int index) const { return v[index]; }
  inline bool operator<(const UEdgeKey &) const;
  inline bool operator==(const UEdgeKey &o) const { return v[0] == o.v[0] && v[1] == o.v[1]; }
  inline bool operator!=(const UEdgeKey &o) const { return v[0] != o.v[0] || v[1] != o.v[1]; }

  inline bool shareIndex(const UEdgeKey &o) const { return v[0] == o.v[0] || v[0] == o.v[1] || v[1] == o.v[0] || v[1] == o.v[1]; }

protected:
  int v[2];
};

inline std::ostream &operator<<(std::ostream &s, const UEdgeKey &v);

// oriented edge
struct OEdgeKey
{
  inline OEdgeKey(int v0, int v1);
  inline OEdgeKey(const std::pair<int, int> &p):
    OEdgeKey(p.first, p.second) {}
  inline OEdgeKey();  // creates an invalid key with v = {-1,-1}
  const int &operator[](int index) const { return v[index]; }
  inline bool operator<(const OEdgeKey &) const;
  inline bool operator==(const OEdgeKey &o) const { return v[0] == o.v[0] && v[1] == o.v[1]; }
  inline bool operator!=(const OEdgeKey &o) const { return v[0] != o.v[0] || v[1] != o.v[1]; }
  inline void reverse() { std::swap(v[0], v[1]); }                             // reverse the orientation of this edge
  inline OEdgeKey getReversedEdgeKey() const { return OEdgeKey(v[1], v[0]); }  // return reversed edge key

protected:
  int v[2];
};

inline std::ostream &operator<<(std::ostream &s, const OEdgeKey &v);

///////////////////////////////////////////////////////////////////////////////
//                             IMPLEMENTATION                                //
///////////////////////////////////////////////////////////////////////////////

inline UEdgeKey::UEdgeKey(int v0, int v1)
{
  if (v0 < v1) {
    v[0] = v0;
    v[1] = v1;
  }
  else {
    v[0] = v1;
    v[1] = v0;
  }
}

inline UEdgeKey::UEdgeKey()
{
  v[0] = v[1] = -1;
}

inline bool UEdgeKey::operator<(const UEdgeKey &other) const
{
  if (v[0] < other.v[0])
    return true;
  if (v[0] > other.v[0])
    return false;
  return v[1] < other.v[1];
}

inline OEdgeKey::OEdgeKey(int v0, int v1)
{
  v[0] = v0;
  v[1] = v1;
}

inline OEdgeKey::OEdgeKey()
{
  v[0] = v[1] = -1;
}

inline bool OEdgeKey::operator<(const OEdgeKey &other) const
{
  if (v[0] < other.v[0])
    return true;
  if (v[0] > other.v[0])
    return false;
  return v[1] < other.v[1];
}

inline std::ostream &operator<<(std::ostream &s, const UEdgeKey &v)
{
  return s << '(' << v[0] << ' ' << v[1] << ')';
}

inline std::ostream &operator<<(std::ostream &s, const OEdgeKey &v)
{
  return s << '(' << v[0] << ' ' << v[1] << ')';
}

}  // namespace Mesh
}  // namespace pgo

namespace std
{
template<>
struct hash<pgo::Mesh::UEdgeKey>
{
  size_t operator()(const pgo::Mesh::UEdgeKey &k) const
  {
    static_assert(sizeof(int) * 2 == sizeof(uint64_t), "uint64_t is not twice the same size as int");
    uint64_t v = (((uint64_t)k[0]) + ((uint64_t)(k[1]) << (sizeof(int) * 8)));
    return std::hash<uint64_t>()(v);
  }
};

template<>
struct hash<pgo::Mesh::OEdgeKey>
{
  size_t operator()(const pgo::Mesh::OEdgeKey &k) const
  {
    static_assert(sizeof(int) * 2 == sizeof(uint64_t), "uint64_t is not twice the same size as int");
    uint64_t v = (((uint64_t)k[0]) + ((uint64_t)(k[1]) << (sizeof(int) * 8)));
    return std::hash<uint64_t>()(v);
  }
};
}  // namespace std