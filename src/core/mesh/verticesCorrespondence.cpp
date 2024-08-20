/*************************************************************************
 *                                                                       *
 *                                                                       *
 * "mesh" library , Copyright (C) 2018 USC                               *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Yijing Li, Jernej Barbic                                *
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

#include "verticesCorrespondence.h"
#include "boundingVolumeTree.h"
#include "valueIndex.h"

std::pair<int, double> pgo::Mesh::findVertexCorrespondences(BasicAlgorithms::ArrayRef<Vec3d> fullGroup, BasicAlgorithms::ArrayRef<Vec3d> subGroup, int *sub2fullMap)
{
  VertexBVTree fullTree;
  fullTree.buildByInertiaPartition(fullGroup);

  MaxValueIndex mvi;
  for (int subID = 0; subID < subGroup.size(); subID++) {
    Vec3d pos = subGroup[subID];
    int fullID = fullTree.getClosestVertex(fullGroup, pos);
    sub2fullMap[subID] = fullID;
    mvi.update((pos - fullGroup[fullID]).squaredNorm(), subID);
  }
  return std::make_pair(mvi.index, sqrt(mvi.value));
}
