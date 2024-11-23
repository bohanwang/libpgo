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

#include "tetMeshManifold.h"

#include "pgoLogging.h"

#include <algorithm>
#include <iostream>
#include <cstring>

namespace pgo::Mesh
{

TetMeshManifold::TetMeshManifold()
{
}

TetMeshManifold::~TetMeshManifold()
{
  clear();
}

void TetMeshManifold::clear()
{
  for (auto &p : tets)
    delete p.second;
  for (auto &p : triangles)
    delete p.second;

  tets.clear();
  triangles.clear();
  surface.clear();
}

const TetMeshManifold::Tetrahedron *TetMeshManifold::add(const int vtx[4], int id)
{
  return add(vtx[0], vtx[1], vtx[2], vtx[3], id);
}

const TetMeshManifold::Tetrahedron *TetMeshManifold::add(int v0, int v1, int v2, int v3, int id)
{
  UTetKey tetkey(v0, v1, v2, v3);
  auto it = tets.find(tetkey);
  // found this tet in the manifold
  if (it != tets.end())
    return it->second;

  Tetrahedron *tet = new Tetrahedron(v0, v1, v2, v3);
  tet->id = id;

  UTriKey utrikeys[4];
  OTriKey otrikeys[4];
  Triangle *faces[4];
  bool fail = false;
  // find all four faces in the manifold and test potential non-manifold adding
  for (int i = 0; i < 4; i++) {
    utrikeys[i] = tet->uFaceKey(i);
    otrikeys[i] = tet->oFaceKey(i);
    auto triit = triangles.find(utrikeys[i]);
    if (triit == triangles.end())  // this triangle face is not inside triangles; it's new
    {
      faces[i] = nullptr;
      PGO_ALOG(surface.find(otrikeys[i]) == surface.end());  // it should also not inside surface
    }
    else {
      faces[i] = triit->second;
      PGO_ALOG(faces[i]->tet[0] != nullptr);
      if (faces[i]->tet[1] != nullptr) {
        // this triangle has two tets already. Mesh would become non-manifold if the tet is added
        fail = true;
        break;
      }
      if (surface.find(otrikeys[i].getReversedTriKey()) == surface.end())  // the orientation of this face is wrong
      {
        fail = true;
        break;
      }
    }
  }
  if (fail) {
    delete tet;
    return nullptr;
  }

  // Now it's safe. Add this tet to the manifold
  tets[tetkey] = tet;
  // OTriKey surfaceFaceToAdd[4];
  // int numSurfaceFaceToAdd = 0;

  for (int i = 0; i < 4; i++) {
    Triangle *tri = faces[i];
    UTriKey &facekey = utrikeys[i];
    OTriKey &okey = otrikeys[i];

    if (tri == nullptr) {
      // no this triangle exist. We add a new triangle
      tri = new Triangle(facekey);
      triangles[facekey] = tri;
      tri->tet[0] = tet;
      tet->face[i] = tri;
      // PGO_ALOG(surface.find(okey) == surface.end());
      surface.insert(std::pair<OTriKey, Tetrahedron *>(okey, tet));
      //      cout << "Begin to add surface triangle: " << endl;
      //
      // surfaceFaceToAdd[numSurfaceFaceToAdd++] = okey;
      // we don't add okey to surfaceManifold here now because this might cause the surfaceManifold to become non-manifold
      //      surfaceManifold.add(okey);
    }
    else {  // tet has a neighbor on this face
      PGO_ALOG(surface.find(okey) == surface.end());
      tri->tet[1] = tet;
      tet->face[i] = tri;
      Tetrahedron *nbr = tri->tet[0];
      for (int j = 0; j < 4; j++)
        if (nbr->face[j] == tri) {
          PGO_ALOG(nbr->nbr[j] == nullptr);
          nbr->nbr[j] = tet;
          tet->nbr[i] = nbr;
          break;
        }
      okey.reverse();  // now okey is the neighbor's face
      PGO_ALOG(surface.find(okey) != surface.end());
      surface.erase(okey);
      // bool success = surfaceManifold.remove(okey);
      // PGO_ALOG(success == true);
    }
  }

  // for(int i = 0; i < numSurfaceFaceToAdd; i++) {
  //   const TriMeshManifold::Triangle * ret = surfaceManifold.add(surfaceFaceToAdd[i]);
  //   PGO_ALOG(ret != nullptr);
  // }
  return tet;
}

bool TetMeshManifold::remove(const int vtx[4])
{
  return remove(vtx[0], vtx[1], vtx[2], vtx[3]);
}

bool TetMeshManifold::remove(int v0, int v1, int v2, int v3)
{
  UTetKey tetkey(v0, v1, v2, v3);
  auto it = tets.find(tetkey);
  // can't find this tet in the manifold
  if (it == tets.end())
    return false;

  // OTriKey surfaceFaceToAdd[4];
  // int numSurfaceFaceToAdd = 0;

  Tetrahedron *tet = it->second;
  for (int i = 0; i < 4; i++) {
    UTriKey trikey = tet->uFaceKey(i);
    auto triit = triangles.find(trikey);
    PGO_ALOG(triit != triangles.end());
    Triangle *tri = triit->second;
    PGO_ALOG(tri->tet[0] == tet || tri->tet[1] == tet);
    PGO_ALOG(!(tri->tet[1] == tet && tri->tet[0] == nullptr));
    PGO_ALOG(!(tri->tet[0] == tet && tri->tet[1] == tet));

    OTriKey okey = tet->oFaceKey(i);
    Tetrahedron *nbr = nullptr;
    if (tri->tet[0] != tet)
      nbr = tri->tet[0];
    else
      nbr = tri->tet[1];

    if (nbr) {
      // we have a neighbor for this tet
      // we'll reset the corresponding neighbor var
      PGO_ALOG(tet->nbr[i] == nbr);
      for (int j = 0; j < 4; j++)
        if (nbr->face[j] == tri) {
          nbr->nbr[j] = nullptr;
          break;
        }
      // we remove this tet's reference in tri
      // and move neighbor to tet[0] so that tri->tet[0] is always not nullptr
      tri->tet[0] = nbr;
      tri->tet[1] = nullptr;

      okey.reverse();  // now it's neighbor's face
      PGO_ALOG(surface.find(okey) == surface.end());
      surface.insert(std::pair<OTriKey, Tetrahedron *>(okey, tet));

      // surfaceFaceToAdd[numSurfaceFaceToAdd++] = okey;
      // we don't add okey to surfaceManifold here now because this might cause the surfaceManifold to become non-manifold
      //      surfaceManifold.add(okey);
    }
    else {
      // only this tet shares this face
      // we'll delete this face
      triangles.erase(triit);
      delete tri;
      PGO_ALOG(surface.find(okey) != surface.end());
      surface.erase(okey);
      // bool success = surfaceManifold.remove(okey);
      // PGO_ALOG(success == true);
    }
  }

  // for(int i = 0; i < numSurfaceFaceToAdd; i++) {
  //   const TriMeshManifold::Triangle * ret = surfaceManifold.add(surfaceFaceToAdd[i]);
  //   PGO_ALOG(ret != nullptr);
  // }

  tets.erase(it);
  delete tet;
  return true;
}

void TetMeshManifold::getSurface(std::set<OTriKey> &outputSurface) const
{
  outputSurface.clear();
  for (const auto &p : surface)
    outputSurface.insert(p.first);
}

TetMeshManifold::Tetrahedron::Tetrahedron(int v0, int v1, int v2, int v3):
  OTetKey(v0, v1, v2, v3), id(-1)
{
  //  v[0] = v0; v[1] = v1; v[2] = v2; v[3] = v3;
  memset(face, 0, sizeof(face));
  memset(nbr, 0, sizeof(nbr));
  // cout << "tetrahedron created: " << v[0] << " " << v[1] << " " << v[2] << " " << v[3] << endl;
}

int TetMeshManifold::Tetrahedron::getOppositeVtx(const Triangle *tri) const
{
  for (int i = 0; i < 4; i++)
    if (face[i] == tri)
      return v[i];
  return -1;
}

bool TetMeshManifold::Tetrahedron::isNeighbor(const TetMeshManifold::Tetrahedron *other) const
{
  for (int i = 0; i < 4; i++)
    if (nbr[i] == other)
      return true;
  return false;
}

int TetMeshManifold::Tetrahedron::getNeighborInvertedIndex(const TetMeshManifold::Tetrahedron *other) const
{
  for (int i = 0; i < 4; i++)
    if (nbr[i] == other)
      return i;
  return -1;
}

TetMeshManifold::Triangle::Triangle(const UTriKey &key):
  UTriKey(key)
{
  memset(tet, 0, sizeof(tet));
}

}  // namespace pgo::Mesh