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

#include "triKey.h"

#include <algorithm>
#include <ostream>

namespace pgo
{
namespace Mesh
{

// unoriented tet key based on vtx indices
struct UTetKey
{
  inline UTetKey();  // creates an invalid key with v = {-1,-1,-1,-1}
  inline UTetKey(int v0, int v1, int v2, int v3);
  inline UTetKey(const int v[4]):
    UTetKey(v[0], v[1], v[2], v[3]) {}
  inline UTetKey(const Vec4i &v):
    UTetKey(v[0], v[1], v[2], v[3]) {}

  inline bool operator<(const UTetKey &o) const { return v < o.v; }
  inline bool operator==(const UTetKey &o) const { return v == o.v; }
  inline bool operator!=(const UTetKey &o) const { return v != o.v; }
  const int &operator[](int index) const { return v[index]; }
  const int *indices() const { return &v[0]; }
  UTriKey uFaceKey(int ind) const;

  inline bool shareUFace(const UTetKey &nbr) const;
  inline UTriKey getSharedUFace(const UTetKey &nbr) const;  // return the first shared UTriKey; return (-1,-1, -1) if no shared UTriKey

  // if all indices are >= 0
  inline bool allVerticesValid() const { return v[0] >= 0; }
  // if this tet is valid: v[i] >=0 && v[i] != v[j]
  inline bool isValidTet() const;
  inline bool hasIndex(int vtxIndex) const { return v[0] == vtxIndex || v[1] == vtxIndex || v[2] == vtxIndex || v[3] == vtxIndex; }

  // opposite face for each vtx in a tet. The faces are ordered so that its normals pointing outside the tet if tet has positive orientation
  // static const int tetFaceIndex[4][3];
protected:
  Vec4i v;
};

inline std::ostream &operator<<(std::ostream &s, const UTetKey &v);

// oriented tet key based on vtx indices
struct OTetKey
{
  inline OTetKey();  // creates an invalid key with v = {-1,-1,-1,-1}
  inline OTetKey(int v0, int v1, int v2, int v3);
  inline OTetKey(const int v[4]):
    OTetKey(v[0], v[1], v[2], v[3]) {}
  inline OTetKey(const Vec4i &v):
    OTetKey(v[0], v[1], v[2], v[3]) {}

  inline bool operator<(const OTetKey &o) const { return v < o.v; }
  inline bool operator==(const OTetKey &o) const { return v == o.v; }
  inline bool operator!=(const OTetKey &o) const { return v != o.v; }
  const int &operator[](int index) const { return v[index]; }
  const int *indices() const { return &v[0]; }

  // get the TriKey opposite to i-th vtx,  // 0 <= i < 4
  OTriKey oFaceKey(int i) const;
  UTriKey uFaceKey(int i) const;
  // get the i-th EdgeKey, // 0 <= i < 6
  //  OEdgeKey oEdgeKey(int i) const; // there's no meaning in producing an OEdgeKey, so it's commented out
  UEdgeKey uEdgeKey(int i) const;
  UTetKey uTetKey() const { return UTetKey(&v[0]); }
  // given the global vtx index, return its first index in v [0,3]
  // else return -1
  inline int getInvertedIndex(const int globalVtxIndex) const { return pgo::getInvertedIndex(v, globalVtxIndex); }

  inline int getInvertedEdgeIndex(const UEdgeKey &edge) const;
  inline int getInvertedTriIndex(const OTriKey &tri) const;
  inline int getInvertedTriIndex(const UTriKey &tri) const;

  inline bool hasIndex(int vtxIndex) const { return v[0] == vtxIndex || v[1] == vtxIndex || v[2] == vtxIndex || v[3] == vtxIndex; }
  // permute v0, v1, v2, v3 and store into r0, r1, r2, r3 so that they share the same orientation but r0 = min(v0,v1,v2,v3)
  static void permute(int v0, int v1, int v2, int v3, int &r0, int &r1, int &r2, int &r3);
  // opposite face for each vtx in a tet. The faces are ordered so that its normals pointing outside the tet if tet has positive orientation
  static const int tetFaceIndex[4][3];
  static const int tetEdgeIndex[6][2];

  Vec4i v;
};

inline std::ostream &operator<<(std::ostream &s, const OTetKey &v);

///////////////////////////////////////////////////////////////////////////////
//                             IMPLEMENTATION                                //
///////////////////////////////////////////////////////////////////////////////

inline UTetKey::UTetKey()
{
  v[0] = v[1] = v[2] = v[3] = -1;
}

inline UTetKey::UTetKey(int v0, int v1, int v2, int v3)
{
  v[0] = v0;
  v[1] = v1;
  v[2] = v2;
  v[3] = v3;
  std::sort(&v[0], &v[4]);
}

inline UTriKey UTetKey::uFaceKey(int ind) const
{
  return UTriKey(v[OTetKey::tetFaceIndex[ind][0]], v[OTetKey::tetFaceIndex[ind][1]], v[OTetKey::tetFaceIndex[ind][2]]);
}

inline bool UTetKey::shareUFace(const UTetKey &nbr) const
{
  for (int i = 0; i < 4; i++) {
    UTriKey key = uFaceKey(i);
    for (int j = 0; j < 4; j++)
      if (key == nbr.uFaceKey(j))
        return true;
  }
  return false;
}

inline UTriKey UTetKey::getSharedUFace(const UTetKey &nbr) const
{
  for (int i = 0; i < 4; i++) {
    UTriKey key = uFaceKey(i);
    for (int j = 0; j < 4; j++)
      if (key == nbr.uFaceKey(j))
        return key;
  }
  return UTriKey();  // return a default invalid UTriKey
}

inline bool UTetKey::isValidTet() const
{
  return v[0] >= 0 && v[0] != v[1] && v[1] != v[2] && v[2] != v[3];
}

inline OTetKey::OTetKey()
{
  v[0] = v[1] = v[2] = v[3] = -1;
}

inline OTetKey::OTetKey(int v0, int v1, int v2, int v3)
{
  permute(v0, v1, v2, v3, v[0], v[1], v[2], v[3]);
}

inline OTriKey OTetKey::oFaceKey(int ind) const
{
  return OTriKey(v[tetFaceIndex[ind][0]], v[tetFaceIndex[ind][1]], v[tetFaceIndex[ind][2]]);
}

inline UTriKey OTetKey::uFaceKey(int ind) const
{
  return UTriKey(v[tetFaceIndex[ind][0]], v[tetFaceIndex[ind][1]], v[tetFaceIndex[ind][2]]);
}

// inline OEdgeKey OTetKey::oEdgeKey(int i) const
//{
//   return OEdgeKey(v[tetEdgeIndex[i][0]], v[tetEdgeIndex[i][1]]);
// }

inline UEdgeKey OTetKey::uEdgeKey(int i) const
{
  return UEdgeKey(v[tetEdgeIndex[i][0]], v[tetEdgeIndex[i][1]]);
}

inline std::ostream &operator<<(std::ostream &s, const UTetKey &v)
{
  return s << '(' << v[0] << ' ' << v[1] << ' ' << v[2] << ' ' << v[3] << ')';
}

inline std::ostream &operator<<(std::ostream &s, const OTetKey &v)
{
  return s << '(' << v[0] << ' ' << v[1] << ' ' << v[2] << ' ' << v[3] << ')';
}

inline int OTetKey::getInvertedEdgeIndex(const UEdgeKey &edge) const
{
  for (int i = 0; i < 6; i++)
    if (edge == uEdgeKey(i))
      return i;
  return -1;
}

inline int OTetKey::getInvertedTriIndex(const OTriKey &tri) const
{
  for (int i = 0; i < 4; i++)
    if (oFaceKey(i) == tri)
      return i;
  return -1;
}

inline int OTetKey::getInvertedTriIndex(const UTriKey &tri) const
{
  for (int i = 0; i < 4; i++)
    if (uFaceKey(i) == tri)
      return i;
  return -1;
}

}  // namespace Mesh
}  // namespace pgo

namespace std
{
template<>
struct hash<pgo::Mesh::UTetKey>
{
  size_t operator()(const pgo::Mesh::UTetKey &k) const
  {
    static_assert(sizeof(int) * 2 == sizeof(size_t), "size_t is not twice the same size as int");
    size_t u = (((size_t)k[0]) + ((size_t)(k[1]) << (sizeof(int) * 8)));
    size_t v = (((size_t)k[2]) + ((size_t)(k[3]) << (sizeof(int) * 8)));
    pgo::EigenSupport::hashCombine(u, v);
    return u;
  }
};

template<>
struct hash<pgo::Mesh::OTetKey>
{
  size_t operator()(const pgo::Mesh::OTetKey &k) const
  {
    static_assert(sizeof(int) * 2 == sizeof(size_t), "size_t is not twice the same size as int");
    size_t u = (((size_t)k[0]) + ((size_t)(k[1]) << (sizeof(int) * 8)));
    size_t v = (((size_t)k[2]) + ((size_t)(k[3]) << (sizeof(int) * 8)));
    pgo::EigenSupport::hashCombine(u, v);
    return u;
  }
};
}  // namespace std
