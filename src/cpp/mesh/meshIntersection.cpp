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

#include "meshIntersection.h"
#include "boundingVolumeTree.h"
#include "predicates.h"

#include "basicAlgorithms.h"

#include <tbb/parallel_for.h>

std::vector<std::vector<int>> pgo::Mesh::computeTrianglesIntersectingEachTetExact(const TetMeshRef tetMesh, const TriMeshRef triMesh, const TriMeshBVTree &triMeshBVTree)
{
  std::vector<std::vector<int>> tetEmbedTri(tetMesh.numTets());

  tbb::parallel_for(0, tetMesh.numTets(), [&](int tetID) {
    std::array<Vec3d, 4> tet;
    for (int j = 0; j < 4; j++)
      tet[j] = tetMesh.pos(tetID, j);
    BoundingBox tetbb(tet);

    auto toBB = [&](const BoundingBox &bb) {
      return (tetbb.intersect(bb));
    };
    auto toTri = [&](int tri) {
      return intersectTriTet(triMesh.pos(tri, 0).data(), triMesh.pos(tri, 1).data(), triMesh.pos(tri, 2).data(),
        tet[0].data(), tet[1].data(), tet[2].data(), tet[3].data());
    };
    triMeshBVTree.rangeQuery(toBB, toTri, tetEmbedTri[tetID]);
    BasicAlgorithms::sortAndDeduplicate(tetEmbedTri[tetID]);
  });
  return tetEmbedTri;
}

std::vector<std::vector<int>> pgo::Mesh::computeTrianglesIntersectingEachTetExact(const TetMeshRef tetMesh, const TriMeshRef triMesh)
{
  TriMeshBVTree triMeshBVTree;
  triMeshBVTree.buildByInertiaPartition(triMesh);
  return computeTrianglesIntersectingEachTetExact(tetMesh, triMesh, triMeshBVTree);
}
