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

#include "tetMeshGeo.h"
#include "tetKey.h"
#include "geometryQuery.h"

#include "basicAlgorithms.h"
#include "containerHelper.h"
#include "pgoLogging.h"

#include <cstring>
#include <cassert>
#include <iostream>
#include <string>

namespace pgo::Mesh
{

TetMeshRef::TetMeshRef(int numVertices, const double *vertices, int numTets, const int *t):
  numVertices_(numVertices), numTets_(numTets)
{
  positions_ = (const Vec3d *)vertices;
  tets_ = (const Vec4i *)t;
}

TetMeshRef::TetMeshRef(int numVertices, const Vec3d *vertices, int numTets, const int *t):
  numVertices_(numVertices), numTets_(numTets)
{
  positions_ = vertices;
  tets_ = (const Vec4i *)t;
}

TetMeshRef::TetMeshRef(int numVertices, const Vec3d *vertices, int numTets, const Vec4i *t):
  numVertices_(numVertices), numTets_(numTets)
{
  positions_ = vertices;
  tets_ = t;
}

TetMeshRef::TetMeshRef(const std::vector<Vec3d> &vertices, const std::vector<Vec4i> &tets):
  numVertices_(static_cast<int>(vertices.size())), numTets_(static_cast<int>(tets.size()))
{
  positions_ = vertices.data();
  tets_ = tets.data();
}

TetMeshRef::TetMeshRef(const std::vector<Vec3d> &vertices, const std::vector<int> &tets):
  numVertices_(static_cast<int>(vertices.size())), numTets_(static_cast<int>(tets.size()) / 4)
{
  positions_ = vertices.data();
  tets_ = (const Vec4i *)tets.data();
}

Vec3d TetMeshRef::computeTetCenter(int tetID) const
{
  return (pos(tetID, 0) + pos(tetID, 1) + pos(tetID, 2) + pos(tetID, 3)) / 4.0;
}

std::vector<BoundingBox> TetMeshRef::getTetBoundingBoxes() const
{
  std::vector<BoundingBox> bbs;
  bbs.reserve(numTets());
  for (int i = 0; i < numTets(); i++) {
    bbs.emplace_back(positions_, tets_[i]);
  }
  return bbs;
}

Tetrahedron TetMeshRef::getTetrahedron(int tetID) const
{
  return { pos(tetID, 0), pos(tetID, 1), pos(tetID, 2), pos(tetID, 3) };
}

// save to obj mesh
bool TetMeshRef::save(const std::string &filename) const
{
  // open file
  FILE *fout = fopen(filename.c_str(), "w");

  if (!fout) {
    std::cout << "Error: could not write to file " << filename << std::endl;
    return false;
  }

  fprintf(fout, "# Vega mesh file.\n");
  fprintf(fout, "# %d vertices, %d elements\n", numVertices(), numTets());
  fprintf(fout, "\n");

  // write vertices
  fprintf(fout, "*VERTICES\n");
  fprintf(fout, "%d 3 0 0\n", numVertices());

  for (int i = 0; i < numVertices(); i++) {
    const Vec3d &v = pos(i);
    fprintf(fout, "%d %.15G %.15G %.15G\n", i + 1, v[0], v[1], v[2]);
  }
  fprintf(fout, "\n");

  // write elements
  fprintf(fout, "*ELEMENTS\n");

  fprintf(fout, "TET\n");

  fprintf(fout, "%d %d 0\n", numTets(), 4);

  for (int el = 0; el < numTets(); el++) {
    fprintf(fout, "%d ", el + 1);
    for (int j = 0; j < 4; j++) {
      fprintf(fout, "%d", tet(el)[j] + 1);
      if (j != 3)
        fprintf(fout, " ");
    }
    fprintf(fout, "\n");
  }
  fprintf(fout, "\n");
  fclose(fout);

  std::cout << "Saved mesh (#v: " << numVertices() << ", #t: " << numTets() << ") to " << filename << "." << std::endl;
  return true;
}

TetMeshGeo::TetMeshGeo(int numVertices, const double *vertices, int numTets, const int *t):
  positions_(numVertices), tets_(numTets)
{
  memcpy(&positions_[0][0], vertices, sizeof(double) * numVertices * 3);
  memcpy(&tets_[0][0], t, sizeof(int) * numTets * 4);
}

TetMeshGeo::TetMeshGeo(std::vector<Vec3d> vertices, std::vector<Vec4i> tets):
  positions_(std::move(vertices)), tets_(std::move(tets))
{
}

TetMeshGeo::TetMeshGeo(int numVertices, const Vec3d *vertices, std::vector<Vec4i> tets):
  positions_(numVertices), tets_(std::move(tets))
{
  memcpy(&positions_[0][0], vertices, sizeof(double) * numVertices * 3);
}

/////////////////////////////////////////////////////////////
//                    TetNeighbor
/////////////////////////////////////////////////////////////

TetNeighbor::TetNeighbor(const std::vector<Vec4i> &tets):
  nTets(static_cast<int>(tets.size())),
  tetNbrs(tets.size(), Vec4i(-1))
{
  for (int tetID = 0; tetID < nTets; tetID++) {
    OTetKey tetkey(tets[tetID]);
    for (int j = 0; j < 4; j++) {
      int v0 = tets[tetID][j];
      OTriKey key = tetkey.oFaceKey(tetkey.getInvertedIndex(v0));
      if (key.isValidTriangle() == false)  // // error, TetMesh is not valid!
        throw 1;

      if (otriTet.find(key) != otriTet.end()) {
        std::cout << "Error, the tet neighbor is not manifold, tetID " << tetID << " tet " << tets[tetID] << " "
                  << "key " << key << " prevTetID: " << otriTet[key] << " prevTet " << tets[otriTet[key]] << std::endl;
        throw 1;  // this ordered face appears twice, not manifold!
      }
      otriTet.emplace(key, tetID);
    }
  }

  for (int tetID = 0; tetID < nTets; tetID++) {
    OTetKey tetkey(tets[tetID]);
    for (int j = 0; j < 4; j++) {
      if (tetNbrs[tetID][j] < 0) {
        int v0 = tets[tetID][j];
        OTriKey key = tetkey.oFaceKey(tetkey.getInvertedIndex(v0));
        auto iter = otriTet.find(key.getReversedTriKey());  // get the revserse key
        if (iter != otriTet.end()) {
          int tetID2 = iter->second;
          tetNbrs[tetID][j] = tetID2;
        }
      }
    }
  }
}

std::vector<std::pair<int, OTriKey>> TetNeighbor::findTetBoundaries(int numTets, const Vec4i *tets) const
{
  PGO_ALOG(tetNbrs.size() == (size_t)numTets);
  std::vector<std::pair<int, OTriKey>> ret;
  for (size_t tetID = 0; tetID < tetNbrs.size(); tetID++) {
    const auto &n = tetNbrs[tetID];
    const auto &t = tets[tetID];
    for (int i = 0; i < 4; i++) {
      if (n[i] < 0)  // find a bounary
      {
        ret.emplace_back((int)tetID, OTriKey(t[i], t[(i + 1) % 4], t[(i + 2) % 4]));
      }
    }
  }
  return ret;
}

NonManifoldTetNeighbor::NonManifoldTetNeighbor(const std::vector<Vec4i> &tets):
  nTets(static_cast<int>(tets.size()))
{
  for (int tetID = 0; tetID < nTets; tetID++) {
    OTetKey tetkey(tets[tetID]);
    for (int j = 0; j < 4; j++) {
      int v0 = tets[tetID][j];
      OTriKey key = tetkey.oFaceKey(tetkey.getInvertedIndex(v0));
      if (key.isValidTriangle() == false)  // // error, TetMesh is not valid!
        throw 1;
      otriTets[key].push_back(tetID);
    }
  }
}

std::vector<int> NonManifoldTetNeighbor::getTetNeighbors(const std::vector<Vec4i> &tets, int tetID) const
{
  PGO_ALOG(tets.size() == (size_t)nTets);
  std::vector<int> ret;
  OTetKey tetkey(tets[tetID]);
  for (int i = 0; i < 4; i++) {
    OTriKey nbrKey = tetkey.oFaceKey(i).getReversedTriKey();
    auto it = otriTets.find(nbrKey);
    if (it != otriTets.end()) {
      BasicAlgorithms::vectorInsertRangeBack(ret, it->second);
    }
  }
  BasicAlgorithms::sortAndDeduplicate(ret);
  return ret;
}

TetMeshGeo getSubTetMesh(const TetMeshRef &tetMesh, const std::vector<int> &subTetIDs,
  std::vector<int> *subVtxID2FullVtxID, std::map<int, int> *fullTetID2SubTetID)
{
  std::vector<int> usedVtxIDs;
  for (int tetID : subTetIDs) {
    for (int i = 0; i < 4; i++)
      usedVtxIDs.push_back(tetMesh.tet(tetID)[i]);
  }
  BasicAlgorithms::sortAndDeduplicate(usedVtxIDs);

  std::map<int, int> oldToNewVtxMap;
  for (size_t i = 0; i < usedVtxIDs.size(); i++) {
    oldToNewVtxMap[usedVtxIDs[i]] = static_cast<int>(i);
  }

  std::vector<Vec4i> newTets(subTetIDs.size());
  int newTetID = 0;
  for (int tetID : subTetIDs) {
    for (int i = 0; i < 4; i++) {
      int newID = oldToNewVtxMap[tetMesh.tet(tetID)[i]];
      newTets[newTetID][i] = newID;
    }
    newTetID++;
  }

  std::vector<Vec3d> newPos(usedVtxIDs.size());
  for (size_t i = 0; i < usedVtxIDs.size(); i++)
    newPos[i] = tetMesh.pos(usedVtxIDs[i]);
  if (subVtxID2FullVtxID) {
    *subVtxID2FullVtxID = std::move(usedVtxIDs);
  }
  if (fullTetID2SubTetID) {
    *fullTetID2SubTetID = std::move(oldToNewVtxMap);
  }

  return TetMeshGeo(std::move(newPos), std::move(newTets));
}

double TetMeshRef::computeTetDeterminant(int tetID) const
{
  return getTetDeterminant(pos(tetID, 0), pos(tetID, 1), pos(tetID, 2), pos(tetID, 3));
}

void TetMeshRef::computeTetBarycentricWeights(int tetID, const Vec3d &queryPos, double w[4]) const
{
  return getTetBarycentricWeights(queryPos, pos(tetID, 0), pos(tetID, 1), pos(tetID, 2), pos(tetID, 3), w);
}

double TetMeshRef::computeSquaredDistanceToTet(int tetID, const Vec3d &queryPos) const
{
  return getSquaredDistanceToTet(queryPos, pos(tetID, 0), pos(tetID, 1), pos(tetID, 2), pos(tetID, 3));
}

}  // namespace pgo::Mesh