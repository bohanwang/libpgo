/*************************************************************************
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

#include "boundingBox.h"

#include <vector>
#include <functional>
#include <tuple>
#include <cfloat>
#include <cassert>
#include <stack>

namespace pgo
{
namespace Mesh
{
class BoundingVolumeTreeBase
{
public:
  int getDepth() const { return depth; }  // get depth of the tree, which equals to the maximum Node::depth in all nodes
  int getNumNodes() const { return static_cast<int>(nodes.size()); }
  const BoundingBox &getRootBoundingBox() const
  {
    assert(nodes.size() > 0);
    return nodes[0].bb;
  }

  int getRootNodeID() const { return nodes.size() > 0 ? 0 : -1; }
  int getNumChildNodes(int nodeID) const { return static_cast<int>(nodes[nodeID].childrenIDs.size()); }
  int getChildNodeID(int nodeID, int j) const { return nodes[nodeID].childrenIDs[j]; }

  int getNumElements(int nodeID) const { return static_cast<int>(nodes[nodeID].indices.size()); }
  int getElementID(int nodeID, int j) const { return static_cast<int>(nodes[nodeID].indices[j]); }

  const BoundingBox &getNodeBoundingBox(int nodeID) const { return nodes[nodeID].bb; }
  //  const std::vector<int> & getChildNodeIDs(int nodeID) const { return nodes[nodeID].childrenIDs; }

  // return true if the queried object intersects bb
  using BBFilter = std::function<bool(const BoundingBox &bb)>;
  // return true if the queried object intersects the element
  using ElementFilter = std::function<bool(int eleID)>;
  // process each element, usually is defined to check if the queried object intersects the element
  // return false if aborting the entire hierarchical process
  using ElementProcess = std::function<void(int eleID)>;

  // Query the bounding volume tree:
  // Starting at the root, if the bounding box query toBB returns true,
  // recursively visit each children node, until reaching a leaf node,
  // where the element query toEle is called to each elementIDs stored in this leaf node.
  void rangeQuery(BBFilter toBB, ElementProcess toEle, std::stack<int> *externalStack = nullptr) const;

  // if element filter returns true, stop the entire hierarchical query process
  // return true only if element filter returns true
  bool rangeQueryWithEarlyStop(BBFilter toBB, ElementFilter toEle) const;

  // Query the bounding volume tree and return hit elementID list:
  // Starting at the root, if the bounding box query toBB returns true,
  // recursively visit each children node, until reaching a leaf node,
  // where if the element query toEle returns true on one elementID stored in this leaf node,
  // the elementID is then pushed into elementIDList.
  // Note: elementIDList is not sorted or deduplicated.
  inline void rangeQuery(BBFilter toBB, ElementFilter toEle, std::vector<int> &elementIDList) const;

  // return the near/far distance of the bounding box bb towards the queried object
  // return < DBL_MAX, --- > if this node should be skipped
  using DistanceToBB = std::function<std::pair<double, double>(const BoundingBox &bb)>;
  // return the distance of the element with elementID towards the queried object
  using DistanceToElement = std::function<double(int elementID, double minDistanceSoFar)>;

  // Nearest object query:
  // Starting at the root, call toBB to check the near/far distance of the bounding box towards the queried object.
  // If the near distance is farther than the closest distance found so far, then skip.
  // Otherwise, recursively visit each children node, until reaching a leaf node,
  // where if calling toElement to compute the distance between the element and the queried object gives
  // a distance smaller than the closest distance found so far. In this case, the closest distance will be updated.
  // distanceHigherBound: the value to initialze the "closest distance found so far" value in the algorithm.
  // Usually it is set to be DBL_MAX. But it can be other values to cull away elements that are too far away as a priori.
  // distanceLowerBound: when the value is found related to an element, then the search is done.
  // Usually it is set to be 0.0. But it can be other >0 values if you just need a close-enough element.
  // return the minimum distance found
  template<class DistanceToBB, class DistanceToElement>
  double nearestQuery(std::vector<std::tuple<double, double, int>> &nodeStack, DistanceToBB toBB, DistanceToElement toElement, double distanceHigherBound = DBL_MAX,
    double distanceLowerBound = 0.0) const;
  template<class DistanceToBB, class DistanceToElement>
  double nearestQuery(DistanceToBB toBB, DistanceToElement toElement, double distanceHigherBound = DBL_MAX,
    double distanceLowerBound = 0.0) const;

  using ElementPairProcess = std::function<void(int elementID, int otherElementID)>;

  void intersectionQuery(const BoundingVolumeTreeBase &otherTree, ElementPairProcess) const;

  void selfIntersectionQuery(ElementPairProcess) const;

protected:
  BoundingVolumeTreeBase() {}
  ~BoundingVolumeTreeBase() {}

  using Partition = std::function<void(const std::vector<int> &elements, const BoundingBox &bb,
    std::vector<std::vector<int>> &childIDs, std::vector<BoundingBox> &subBBs)>;

  // derived class calls this function to build the bounding volume tree
  void buildFromRoot(Partition divideNode, int maxDepth, int maxNumElementsPerNode);

  // whether a bounding box intersects an element
  using ElementBBIntersect = std::function<bool(int eleID, const BoundingBox &bb)>;
  // given element ID list, create a bounding box that covers them
  using BBFromElementList = std::function<BoundingBox(const std::vector<int> &elementList)>;

  // an implementation to divide the elements within the bounding box bb into eight sub groups, store them
  // in childIDs and their bounding boxes in subBBs
  // BBFromElementList: create a bounding box from elements
  // ElementBBIntersect: whether a bounding box intersect with an element
  // shrinkBoundingBoxes: whether to shrink the bounding box to tightly fit the elements after bounding box subdivision
  //                      Shrinking reduce the size of bounding boxes, improving performance.
  //                      But you may also want to not shrink when building an inexact octree to prevent elements lying on the bounding box faces.
  //                      In an inexact octree, elements lying on the bounding box faces may counted as not intersecting with the bounding box
  //                      Due to numerical issues.
  static void boundingBoxOctreePartition(const std::vector<int> &elements, const BoundingBox &bb,
    std::vector<std::vector<int>> &childIDs, std::vector<BoundingBox> &subBBs,
    const std::vector<BoundingBox> &elementBBs, BBFromElementList, ElementBBIntersect, bool shrinkBoundingBoxes = true);

  using ElementCenterOfMass = std::function<Vec3d(int eleID)>;
  using ElementMass = std::function<double(int eleID)>;
  using ElementInertiaTensor = std::function<Mat3d(int eleID, const Vec3d &referencePoint)>;
  // partition elements based on its inertia distribution
  static void elementsInertiaPartition(const std::vector<int> &elements, const BoundingBox &bb,
    int numPartitions, std::vector<std::vector<int>> &childIDs, ElementMass, ElementCenterOfMass, ElementInertiaTensor);
  // binary partition
  // if children size differ larger than childSizeBalanceThreshold, forcely rebalance them
  static void elementsInertiaBinaryPartition(const std::vector<int> &elements, const BoundingBox &bb,
    std::vector<std::vector<int>> &childIDs, int childSizeBalanceThreshold,
    ElementMass, ElementCenterOfMass, ElementInertiaTensor);

  // whether elementID is covered by nodeID. This functions calls itself recursively.
  bool findElementRecursive(int nodeID, int elementID) const;

  // get all elementIDs under the nodeID (including all its descendent nodes)
  // returned vector is sorted and deduplicated
  std::vector<int> getNodeCoveredElements(int nodeID) const;

  // elementIDs are not sorted or deduplicated
  void getCoveredElementsRecursive(int nodeID, std::vector<int> &elementIDs) const;

  struct Node
  {
    int depth;  // for the root node of the tree, depth = 0
    std::vector<int> indices;
    BoundingBox bb;
    std::vector<int> childrenIDs;

    Node(int depth, std::vector<int> indices):
      depth(depth), indices(std::move(indices)) {}
    bool isLeaf() const { return indices.size() > 0; }
    // get the next available child index in childrenIDs starting at nextChildIndex
    // return 8 if not available
    // return nextChildIndex if childrenIDs[nextChildIndex] is available
    // nextChildIndex: [0, 8]
    int getNextChild(int nextChildIndex) const { return (nextChildIndex < numChildren() ? nextChildIndex + 1 : nextChildIndex); }
    int numChildren() const { return static_cast<int>(childrenIDs.size()); }
  };

  int depth{ 0 };  // how deep in the hierarchy is this bounding volume tree

  std::vector<Node> nodes;
};

// =====================================================
// ============= BELOW ARE IMPLEMENTATION ==============
// =====================================================

inline void BoundingVolumeTreeBase::rangeQuery(BBFilter toBB, ElementFilter toEle, std::vector<int> &elementIDList) const
{
  auto elementProcess = [&](int eleID) {
    if (toEle(eleID))
      elementIDList.push_back(eleID);
  };
  rangeQuery(toBB, elementProcess);
}

}  // namespace Mesh
}  // namespace pgo
