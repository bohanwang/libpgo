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

#include "triMeshManifold.h"

#include "pgoLogging.h"

#include <algorithm>
#include <iostream>
#include <cstring>
#include <cassert>

namespace pgo
{
namespace Mesh
{

TriMeshManifold::TriMeshManifold()
{
}

TriMeshManifold::~TriMeshManifold()
{
  clear();
}

void TriMeshManifold::clear()
{
  for (TriIter it = triangles.begin(); it != triangles.end(); it++)
    delete it->second;

  for (EdgeIter it = edges.begin(); it != edges.end(); it++)
    delete it->second;

  triangles.clear();
  edges.clear();
  boundary.clear();
}

const TriMeshManifold::Triangle *TriMeshManifold::add(int v0, int v1, int v2)
{
  OTriKey trikey(v0, v1, v2);
  return add(trikey);
}

const TriMeshManifold::Triangle *TriMeshManifold::add(const OTriKey &trikey)
{
  //  cout << "adding OTriKey: " << trikey[0] << " " << trikey[1] << " " << trikey[2] << endl;
  TriIter it = triangles.find(trikey);
  // found this triangle in the manifold
  if (it != triangles.end()) {
    //    cout << "Already exists" << endl;
    return it->second;
  }

  Triangle *tri = new Triangle(trikey);
  UEdgeKey edgekeys[3];
  Edge *ownEdges[3];

  bool fail = false;
  // find all three edges in the manifold and test potential non-manifold adding
  for (int i = 0; i < 3; i++) {
    edgekeys[i] = tri->uEdgeKey(i);
    EdgeIter edgeit = edges.find(edgekeys[i]);
    if (edgeit == edges.end())  // this edge is new
    {
      ownEdges[i] = nullptr;
      PGO_ALOG(boundary.find(tri->oEdgeKey(i)) == boundary.end());  // PGO_ALOG this edge is not inside boundary
    }
    else {
      ownEdges[i] = edgeit->second;
      PGO_ALOG(ownEdges[i]->face[0] != nullptr);
      if (ownEdges[i]->face[1] != nullptr) {
        // this edge has two faces already. Mesh would become non-manifold if this triangle is added
        fail = true;
        break;
      }
      // test orientation of this edge is correct or not.
      if (boundary.find(tri->oEdgeKey(i).getReversedEdgeKey()) == boundary.end()) {
        fail = true;
        break;
      }
    }
  }
  if (fail) {
    delete tri;
    return nullptr;
  }

  // Now it's safe. Add this triangle to the manifold
  triangles[trikey] = tri;

  for (int i = 0; i < 3; i++) {
    Edge *edge = ownEdges[i];
    UEdgeKey &edgekey = edgekeys[i];
    OEdgeKey okey = tri->oEdgeKey(i);
    PGO_ALOG(edgekey == UEdgeKey(okey[0], okey[1]));
    if (edge == nullptr) {
      // this edge does not exist yet. We add a new Edge
      edge = new Edge(edgekey);
      PGO_ALOG(edge->v[0] <= edge->v[1]);
      edges[edgekey] = edge;
      edge->face[0] = tri;
      tri->edge[i] = edge;

      // modify boundary
      // PGO_ALOG(boundary.find(okey) == boundary.end()); // this PGO_ALOG has been checked above; so it's commented out
      //      cout << "boundary insertion: " << okey.v[0] << " " << okey.v[1] << endl;
      boundary.insert(okey);
    }
    else  // triangle has a neighbor on this edge
    {
      PGO_ALOG(boundary.find(okey) == boundary.end());
      edge->face[1] = tri;
      tri->edge[i] = edge;
      Triangle *nbr = edge->face[0];
      for (int j = 0; j < 3; j++)
        if (nbr->edge[j] == edge) {
          PGO_ALOG(nbr->nbr[j] == nullptr);
          nbr->nbr[j] = tri;
          tri->nbr[i] = nbr;
          break;
        }

      // update boundary
      okey.reverse();  // now okey is the neighbor's edge
                       //      cout << "removing boundary edge: " << okey.v[0] << " " << okey.v[1] << endl;
      PGO_ALOG(boundary.find(okey) != boundary.end());
      boundary.erase(okey);
    }
  }
  return tri;
}

const TriMeshManifold::Triangle *TriMeshManifold::getTriangle(const OTriKey &key) const
{
  TriCIter it = triangles.find(key);
  if (it != triangles.end())
    return it->second;
  return nullptr;
}

bool TriMeshManifold::remove(int v0, int v1, int v2)
{
  OTriKey trikey(v0, v1, v2);
  return remove(trikey);
}

bool TriMeshManifold::remove(const OTriKey &trikey)
{
  std::cout << "remove OTriKey: " << trikey[0] << " " << trikey[1] << " " << trikey[2] << std::endl;
  TriIter it = triangles.find(trikey);
  // can't find this tet in the manifold
  if (it == triangles.end())
    return false;

  Triangle *tri = it->second;
  for (int i = 0; i < 3; i++) {
    UEdgeKey edgekey = tri->uEdgeKey(i);
    EdgeIter edgeit = edges.find(edgekey);
    PGO_ALOG(edgeit != edges.end());
    Edge *edge = edgeit->second;
    PGO_ALOG(edge->face[0] == tri || edge->face[1] == tri);
    PGO_ALOG(!(edge->face[1] == tri && edge->face[0] == nullptr));
    PGO_ALOG(!(edge->face[0] == tri && edge->face[1] == tri));

    OEdgeKey okey = tri->oEdgeKey(i);
    Triangle *nbr = nullptr;
    if (edge->face[0] != tri)
      nbr = edge->face[0];
    else
      nbr = edge->face[1];

    if (nbr) {
      PGO_ALOG(boundary.find(OEdgeKey(okey[0], okey[1])) == boundary.end());
      PGO_ALOG(boundary.find(OEdgeKey(okey[1], okey[0])) == boundary.end());
      // we have a neighbor for this triangle
      // we'll reset the corresponding neighbor var
      PGO_ALOG(tri->nbr[i] == nbr);
      for (int j = 0; j < 3; j++)
        if (nbr->edge[j] == edge) {
          nbr->nbr[j] = nullptr;
          break;
        }
      // we remove this tri's reference in edge
      // and move neighbor to face[0] so that edge->face[0] is always not nullptr
      edge->face[0] = nbr;
      edge->face[1] = nullptr;

      // update boundary
      okey.reverse();  // now it's neighbor's edge
                       //      cout << "boundary insertion: " << okey.v[0] << " " << okey.v[1] << endl;
      boundary.insert(okey);
    }
    else {
      // only this tri shares this edge
      // we'll delete this edge
      edges.erase(edgeit);
      delete edge;
      // update boundary
      //      cout << "removing boundary edge: " << okey.v[0] << " " << okey.v[1] << endl;
      PGO_ALOG(boundary.find(okey) != boundary.end());
      boundary.erase(okey);
    }
  }
  triangles.erase(it);
  delete tri;
  return true;
}

TriMeshManifold::Edge::Edge(const UEdgeKey &key):
  UEdgeKey(key)
{
  memset(face, 0, sizeof(face));
}

TriMeshManifold::Triangle::Triangle(const OTriKey &key):
  OTriKey(key)
{
  memset(nbr, 0, sizeof(nbr));
  memset(edge, 0, sizeof(edge));
}

}  // namespace Mesh
}  // namespace pgo
