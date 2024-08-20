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

#pragma once

#include <vector>

namespace pgo
{
namespace BasicAlgorithms
{

// low level BFS template routines for efficiency

// breadth-first search template
// IntVector: should have member function: .size() and .push_back()
// getNeighbor(int nodeID) -> range of nbring nodeIDs
// nodeVisited: a boolean array to check whether nodes have been visited
// foundNodeIDs: both as input (seed nodeIDs) and output (all the nodes visited by this BFS)
template<class GetNeighborRange, class BoolArray, class IntVector>
void doBFS(GetNeighborRange getNeighbor, BoolArray &nodeVisited, IntVector &foundNodeIDs)
{
  std::size_t candidateBegin = 0, candidateEnd = foundNodeIDs.size();
  while (candidateBegin != candidateEnd) {
    for (std::size_t i = candidateBegin; i < candidateEnd; i++) {
      int node = foundNodeIDs[i];
      for (int nbr : getNeighbor(node)) {
        if (nbr < 0 || nodeVisited[nbr])
          continue;
        nodeVisited[nbr] = true;
        foundNodeIDs.push_back(nbr);
      }
    }
    candidateBegin = candidateEnd;
    candidateEnd = foundNodeIDs.size();
  }
}

// another BFS template, except instead of using BoolArray to keep track of nodes visited, this one uses:
// ShouldSkipNode(int nodeID) -> whether to skip nodeID because it is visited or forbidden
// MarkNodeVisited(int nodeID) -> mark node as visited
template<class GetNeighborRange, class ShouldSkipNode, class MarkNodeVisited, class IntVector>
void doBFS(GetNeighborRange getNeighbor, ShouldSkipNode skipNode, MarkNodeVisited markNode, IntVector &foundNodeIDs)
{
  std::size_t candidateBegin = 0, candidateEnd = foundNodeIDs.size();
  while (candidateBegin != candidateEnd) {
    for (std::size_t i = candidateBegin; i < candidateEnd; i++) {
      int node = foundNodeIDs[i];
      for (int nbr : getNeighbor(node)) {
        if (nbr < 0 || skipNode(nbr))
          continue;
        markNode(nbr);
        foundNodeIDs.push_back(nbr);
      }
    }
    candidateBegin = candidateEnd;
    candidateEnd = foundNodeIDs.size();
  }
}

// breadth-first search template on directed tree
// since it's a directed tree, nodes will not be visited twice
template<class GetNeighborRange, class IntVector>
void doBFSOnDirectedTree(GetNeighborRange getNeighbor, IntVector &foundNodeIDs)
{
  std::size_t candidateBegin = 0, candidateEnd = foundNodeIDs.size();
  while (candidateBegin != candidateEnd) {
    for (std::size_t i = candidateBegin; i < candidateEnd; i++) {
      int node = foundNodeIDs[i];
      for (int nbr : getNeighbor(node))
        foundNodeIDs.push_back(nbr);
    }
    candidateBegin = candidateEnd;
    candidateEnd = foundNodeIDs.size();
  }
}

template<class GetNeighborRange, class BoolContainer>
std::vector<int> findNodeIDsFromSeedBFS(GetNeighborRange getNeighbor, BoolContainer &nodeVisited, int seed)
{
  if (nodeVisited[seed])
    return {};

  std::vector<int> foundNodeIDs = { seed };
  doBFS(getNeighbor, nodeVisited, foundNodeIDs);

  return foundNodeIDs;
}

}  // namespace Mesh
}  // namespace pgo
