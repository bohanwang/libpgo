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

#include "triKey.h"

#include <map>
#include <set>

namespace pgo
{
namespace Mesh
{
// Maintain a triangle mesh to enfore manifoldness on face-edge connections.
// It can also query neighboring relationship on faces.
// Note that it only detects non-manifoldness created by one edge shared by more than two faces,
// Does not consider the case where two faces share no edge but one vertex.

class TriMeshManifold
{
public:
  TriMeshManifold();
  virtual ~TriMeshManifold();

  // clear internal data
  void clear();

  class Triangle;
  class Edge;

  // add a new triangle key or (v0,v1,v2) into the mesh, return the new built Triangle
  // if the triangle already exists, it just return the triangle
  // return NULL if the new triangle can cause one edge shared by more than two faces (non-manifoldness)
  const Triangle *add(int v0, int v1, int v2);
  const Triangle *add(const OTriKey &key);

  // return a triangle, NULL if not found
  const Triangle *getTriangle(const OTriKey &key) const;

  // remove a triangle from the mesh, return false if it cannot be found
  bool remove(int v0, int v1, int v2);
  bool remove(const OTriKey &key);

  typedef std::map<OTriKey, Triangle *> TriMap;
  typedef std::map<UEdgeKey, Edge *> EdgeMap;

  // get the connection structure for triangles and edges
  const TriMap &getTriMap() const { return triangles; }
  const EdgeMap &getEdgeMap() const { return edges; }
  // get the boundary edges on the mesh
  const std::set<OEdgeKey> &getBoundary() const { return boundary; }

  // unoriented edge, storing its neighboring triangles
  class Edge : UEdgeKey
  {
  public:
    Edge(const UEdgeKey &key);
    const Triangle *getFace(std::size_t ind) const { return face[ind]; }

  protected:
    Triangle *face[2];
    friend class TriMeshManifold;
  };

  // oriented triangle, storing triangle indices, neighboring edges and neighboring triangles
  class Triangle : public OTriKey
  {
  public:
    Triangle(const OTriKey &key);
    const Triangle *getNeighbor(int ind) const { return nbr[ind]; }
    int getVtx(std::size_t i) const { return v[static_cast<int>(i)]; }

  protected:
    Edge *edge[3];
    Triangle *nbr[3];
    friend class TriMeshManifold;
  };

  typedef TriMap::iterator TriIter;
  typedef TriMap::const_iterator TriCIter;
  typedef EdgeMap::iterator EdgeIter;
  typedef EdgeMap::const_iterator EdgeCIter;

protected:
  TriMap triangles;
  EdgeMap edges;
  std::set<OEdgeKey> boundary;
};
}  // namespace Mesh
}  // namespace pgo
