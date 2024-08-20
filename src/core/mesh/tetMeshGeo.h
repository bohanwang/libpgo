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

#include "boundingBox.h"
#include "tetrahedron.h"
#include "triKey.h"

#include <string>
#include <set>
#include <map>
#include <cmath>

namespace pgo
{
namespace Mesh
{
// class to reference an external triangle mesh
class TetMeshRef
{
public:
  TetMeshRef() {}  // empty mesh
  TetMeshRef(int numVertices, const double *vertices, int numTets, const int *tets);
  TetMeshRef(int numVertices, const Vec3d *vertices, int numTets, const int *tets);
  TetMeshRef(int numVertices, const Vec3d *vertices, int numTets, const Vec4i *tets);
  TetMeshRef(const std::vector<Vec3d> &vertices, const std::vector<Vec4i> &tets);
  TetMeshRef(const std::vector<Vec3d> &vertices, const std::vector<int> &tets);

  int numVertices() const { return numVertices_; }
  int numTets() const { return numTets_; }

  const Vec3d &pos(int vtxID) const { return positions_[vtxID]; }
  const Vec4i &tet(int tetID) const { return tets_[tetID]; }

  int tetVtxID(int tetID, int i) const { return tets_[tetID][i]; }
  const Vec3d &pos(int tetID, int i) const { return positions_[tets_[tetID][i]]; }

  Vec3d computeTetCenter(int tetID) const;

  const Vec3d *positions() const { return positions_; }
  const Vec4i *tets() const { return tets_; }

  Tetrahedron getTetrahedron(int tetID) const;

  template<class TetIDContainer>
  void getVerticesInTets(const TetIDContainer &tetIDs, std::vector<int> &vertexIDs) const;

  // get bounding box for every tet
  std::vector<BoundingBox> getTetBoundingBoxes() const;

  double computeTetDeterminant(int tetID) const;

  void computeTetBarycentricWeights(int tetID, const Vec3d &queryPos, double weight[4]) const;

  double computeSquaredDistanceToTet(int tetID, const Vec3d &queryPos) const;
  double computeDistanceToTet(int tetID, const Vec3d &queryPos) const { return std::sqrt(computeSquaredDistanceToTet(tetID, queryPos)); }

  // save to obj mesh
  bool save(const std::string &filename) const;

protected:
  int numVertices_ = 0, numTets_ = 0;
  const Vec3d *positions_ = nullptr;
  const Vec4i *tets_ = nullptr;
};

class TetMeshGeo
{
public:
  TetMeshGeo() {}  // empty mesh
  TetMeshGeo(int numVertices, const double *vertices, int numTets, const int *tets);
  TetMeshGeo(std::vector<Vec3d> vertices, std::vector<Vec4i> tets);
  TetMeshGeo(int numVertices, const Vec3d *vertices, std::vector<Vec4i> tets);

  int numVertices() const { return static_cast<int>(positions_.size()); }
  int numTets() const { return static_cast<int>(tets_.size()); }

  const Vec3d &pos(int vtxID) const { return positions_[vtxID]; }
  Vec3d &pos(int vtxID) { return positions_[vtxID]; }
  const Vec4i &tet(int tetID) const { return tets_[tetID]; }
  Vec4i &tet(int tetID) { return tets_[tetID]; }

  int tetVtxID(int tetID, int i) const { return tets_[tetID][i]; }
  const Vec3d &pos(int tetID, int i) const { return positions_[tets_[tetID][i]]; }
  Vec3d &pos(int tetID, int i) { return positions_[tets_[tetID][i]]; }

  const std::vector<Vec3d> &positions() const { return positions_; }
  const std::vector<Vec4i> &tets() const { return tets_; }

  std::vector<Vec3d> &positions() { return positions_; }
  std::vector<Vec4i> &tets() { return tets_; }

  TetMeshRef ref() const { return { positions_, tets_ }; }
  //  implicit conversion
  operator TetMeshRef() const { return ref(); }

  // save to obj mesh
  bool save(const std::string &filename) const { return ref().save(filename); }

protected:
  std::vector<Vec3d> positions_;
  std::vector<Vec4i> tets_;
};

class TetNeighbor
{
public:
  TetNeighbor() {}
  TetNeighbor(const std::vector<Vec4i> &tets);
  int numTets() const { return nTets; }

  const Vec4i &getTetNeighbors(int tetID) const { return tetNbrs[tetID]; }

  std::vector<std::pair<int, OTriKey>> findTetBoundaries(int numTets, const Vec4i *tets) const;
  inline std::vector<std::pair<int, OTriKey>> findTetBoundaries(const std::vector<Vec4i> &tets) const { return findTetBoundaries((int)(tets.size()), tets.data()); }

protected:
  int nTets = 0;
  std::vector<Vec4i> tetNbrs;
  std::map<OTriKey, int> otriTet;  // <v0,v1,v2> -> tet with ordered face <v0, v1, v2>
};

class NonManifoldTetNeighbor
{
public:
  NonManifoldTetNeighbor() {}  // empty
  NonManifoldTetNeighbor(const std::vector<Vec4i> &tets);

  // return sorted buffer of all tets neighboring tetID
  std::vector<int> getTetNeighbors(const std::vector<Vec4i> &tets, int tetID) const;

protected:
  int nTets = 0;
  std::map<OTriKey, std::vector<int>> otriTets;
};

// return a sub mesh which are triangles from selected triangleIDs
// isolated vertices in the sub mesh ARE NOT removed
// call removeIsolatedVertices to remove them from the result TriMesh
template<class TetIDContainer>
TetMeshGeo getSubTetMeshWithSameVertices(const TetMeshRef &tetMesh, const TetIDContainer &tetIDs);

TetMeshGeo getSubTetMesh(const TetMeshRef &tetMesh, const std::vector<int> &subTetIDs,
  std::vector<int> *subVtxID2FullVtxID = nullptr,
  std::map<int, int> *fullTetID2SubTetID = nullptr);

//////////////////////////////////////////////////////////
//                  Implementation
//////////////////////////////////////////////////////////

template<class TetIDContainer>
void TetMeshRef::getVerticesInTets(const TetIDContainer &tetIDs, std::vector<int> &vertexIDs) const
{
  int numPrevIDs = vertexIDs.size();
  for (int tetID : tetIDs) {
    for (int i = 0; i < 4; i++) {
      vertexIDs.push_back(tets_[tetID][i]);
    }
  }
  // remove duplicate vertex IDs
  std::sort(vertexIDs.begin() + numPrevIDs, vertexIDs.end());
  auto newEnd = std::unique(vertexIDs.begin() + numPrevIDs, vertexIDs.end());
  vertexIDs.resize(std::distance(vertexIDs.begin(), newEnd));
}

template<class TetIDContainer>
TetMeshGeo getSubTetMeshWithSameVertices(const TetMeshRef &tetMesh, const TetIDContainer &tetIDs)
{
  std::vector<Vec4i> subTets;
  for (int tetID : tetIDs) {
    subTets.push_back(tetMesh.tet(tetID));
  }
  return TetMeshGeo(tetMesh.numVertices(), tetMesh.positions(), std::move(subTets));
}
}  // namespace Mesh
}  // namespace pgo
