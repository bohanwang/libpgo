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

#include "tetKey.h"

#include <map>
#include <set>
#include <unordered_map>

// maintain a tet mesh to enfore manifoldness on tet-face connection
// it can also query neighboring relationship on tets
// note that it only detects non-manifoldness created by one face shared by more than two tets,
// does not consider the case where two tets sharing no face but one vertex

// define tet(v0,v1,v2,v3) orientation:
//                                                  |(p1 - p0)^T|
// ((p1 - p0) x (p2 - p0)) dot (p3 - p0) > 0, <==>  |(p2 - p0)^T| > 0
//                                                  |(p3 - p0)^T|
// in struct UTetKey, we order v[4] to be v0 < v1 < v2 < v3,
// in struct OTetKey, we order v[4] to be v0 < v1 < (v2, v3) and the tet orientation is preserved
// we also preserve the input tet vertex order in class TetMeshManifold
// so that we can get correct tet face orientation from Tetrahedron::getOppositeFace(size_t)
//
// in a UTriKey, we order v[3] to be v0 < v1 < v2,
// no ordering is preserved in class TetMeshManifold::Triangle

namespace pgo
{
namespace Mesh
{

class TetMeshManifold
{
public:
  TetMeshManifold();
  virtual ~TetMeshManifold();

  void clear();

  class Triangle;
  class Tetrahedron;
  typedef std::unordered_map<UTetKey, Tetrahedron *> TetMap;
  typedef std::unordered_map<UTriKey, Triangle *> TriMap;
  typedef std::unordered_map<OTriKey, Tetrahedron *> SurfaceMap;

  // add a new tet (v0,v1,v2, v3) into the mesh, return the new built Tetrahedron
  // if the tet already exists, it just return the tet
  // return NULL if the new tet can cause one face shared by more than two tets (non-manifoldness)
  // an id can also be given to the input tet
  const Tetrahedron *add(int v0, int v1, int v2, int v3, int id = -1);
  const Tetrahedron *add(const int vtx[4], int id = -1);

  // remove a tet from the mesh, return false if it cannot be found
  bool remove(int v0, int v1, int v2, int v3);
  bool remove(const int vtx[4]);

  // get the connection structure for tets and triangles
  const TetMap &getTetMap() const { return tets; }
  const TriMap &getTriMap() const { return triangles; }
  // get the surface (boundary faces) of the mesh
  const SurfaceMap &getSurfaceMap() const { return surface; }
  void getSurface(std::set<OTriKey> &surface) const;

  class Triangle : public UTriKey
  {
  public:
    Triangle(const UTriKey &key);
    const Tetrahedron *getTet(std::size_t i) const { return tet[i]; }
    int getVtx(std::size_t i) const { return v[static_cast<int>(i)]; }

  protected:
    Tetrahedron *tet[2];
    friend class TetMeshManifold;
  };

  class Tetrahedron : public OTetKey
  {
  public:
    Tetrahedron(int v0, int v1, int v2, int v3);
    int getVtx(std::size_t i) const { return v[static_cast<int>(i)]; }
    // get the index of vtx in v[4], return -1 if not found
    int getVertexInvertedIndex(int vtx) const { return getInvertedIndex(vtx); }
    const Tetrahedron *getNeighbor(std::size_t i) const { return nbr[i]; }
    bool isNeighbor(const Tetrahedron *other) const;
    int getNeighborInvertedIndex(const Tetrahedron *other) const;  // get the index of other in this->nbr[4], return -1 if not found
    // get oriented opposite face for each vtx in the tet, vtx: 0-3
    void getOppositeFace(int vtx, int face[3]) const;
    OTriKey getOppositeOFace(int vtx) const;  // vtx: 0-3
    UTriKey getOppositeUFace(int vtx) const;  // vtx: 0-3
    int getOppositeVtx(const Triangle *face) const;
    int getID() const { return id; }

  protected:
    int id;
    Triangle *face[4];
    Tetrahedron *nbr[4];
    friend class TetMeshManifold;
  };

protected:
  TetMap tets;
  TriMap triangles;
  SurfaceMap surface;
};

}  // namespace Mesh
}  // namespace pgo
