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

#include "boundingVolumeTreeBase.h"

#include "triple.h"
#include "basicAlgorithms.h"
#include "containerHelper.h"
#include "pgoLogging.h"

#include <cassert>
#include <stack>

/////////////////////////////////////////////////////////////
//                   BoundingVolumeHierarchyBase
/////////////////////////////////////////////////////////////

namespace pgo::Mesh
{
void eigen_sym(const Mat3d &A, Vec3d &eigval, Vec3d eigvec[3])
{
  Eigen::SelfAdjointEigenSolver<Mat3d> eig(A, Eigen::ComputeEigenvectors);
  eigval = eig.eigenvalues();

  Mat3d V = eig.eigenvectors();
  eigvec[0] = V.col(0);
  eigvec[1] = V.col(1);
  eigvec[2] = V.col(2);
}

void BoundingVolumeTreeBase::rangeQuery(BBFilter toBB, ElementProcess toEle, std::stack<int> *externalStack) const
{
  PGO_ALOG(nodes.empty() == false);
  if (toBB(nodes[0].bb) == false)
    return;

  std::stack<int> *nodeStack;
  if (externalStack) {
    nodeStack = externalStack;
    while (!nodeStack->empty())
      nodeStack->pop();
  }
  else {
    nodeStack = new std::stack<int>;
  }

  nodeStack->push(0);
  while (!nodeStack->empty()) {
    int nodeID = nodeStack->top();
    nodeStack->pop();
    const auto &node = nodes[nodeID];

    const auto &indices = node.indices;
    if (indices.size() > 0) {  // leaf node
      for (int eleID : indices)
        toEle(eleID);
    }
    else {
      for (int childID : node.childrenIDs)
        if (childID >= 0 && toBB(nodes[childID].bb) == true)
          nodeStack->push(childID);
    }
  }

  if (externalStack == nullptr) {
    delete nodeStack;
  }
}

bool BoundingVolumeTreeBase::rangeQueryWithEarlyStop(BBFilter toBB, ElementFilter toEle) const
{
  PGO_ALOG(nodes.empty() == false);
  if (toBB(nodes[0].bb) == false)
    return false;
  std::stack<int> nodeStack;
  nodeStack.push(0);
  while (!nodeStack.empty()) {
    int nodeID = nodeStack.top();
    nodeStack.pop();
    const auto &node = nodes[nodeID];
    const auto &indices = node.indices;
    if (indices.size() > 0) {      // leaf node
      for (int eleID : indices)
        if (toEle(eleID) == true)  // stop
          return true;
    }
    else {
      for (int childID : node.childrenIDs)
        if (childID >= 0 && toBB(nodes[childID].bb) == true)
          nodeStack.push(childID);
    }
  }
  return false;
}

void BoundingVolumeTreeBase::intersectionQuery(const BoundingVolumeTreeBase &otherTree, ElementPairProcess elementsIntersect) const
{
  // we have two bounding volume hierarchy, A and B
  // we compare bounding boxes at the same level,
  //   or at immediately adjacent levels, with recursing first into the first node

  // define StackEntry: <// node from A, node from B, index of nodeA/B in its parent's childIDs array>
  // the third value in StackEntry seems ambiguous: nodeA/B?
  // in fact, for a particular node A/B pair, we can determine which node to recurse by using the following function:
  // chooseWhichNodeToRecurse -> return NODEA_CHOSEN or false
  const bool NODEA_CHOSEN = true;
  auto chooseWhichNodeToRecurse = [](const Node *nodeA, const Node *nodeB) {
    // we always try to recurse on B first
    // we recurse on A only if B cannot be recursed, due to B being a leaf or B too deep comparing to A
    return (nodeB->isLeaf() || ((nodeA->depth <= nodeB->depth) && nodeA->isLeaf() == false));
  };

  using StackEntry = triple<int, int, int>;
  std::stack<StackEntry> traversalStack;
  traversalStack.push(StackEntry(0, 0, 0));

  // int iter = 0;
  while (true) {
    StackEntry stackEntry = traversalStack.top();

    // cout << "====== " << iter++ << " ======" << endl;
    // printf("Top off stack: %d vs %d. Child index: %x . Stack size: %d\n", stackEntry.first, stackEntry.second, stackEntry.third, (int)traversalStack.size());

    int nodeAID = stackEntry.first;
    int nodeBID = stackEntry.second;
    PGO_ALOG(nodeAID >= 0 && nodeAID < BasicAlgorithms::sizei(nodes));
    PGO_ALOG(nodeBID >= 0 && nodeBID < BasicAlgorithms::sizei(otherTree.nodes));
    const Node *nodeA = &nodes[nodeAID];
    const Node *nodeB = &otherTree.nodes[nodeBID];

    // check intersection
    bool potentialCollision = nodeA->bb.intersect(nodeB->bb);

    // are the nodes on the deepest level
    // negative sign by convention indicates leaf
    bool deepestLevelA = nodeA->isLeaf();
    bool deepestLevelB = nodeB->isLeaf();
    bool bothLeaves = (deepestLevelA && deepestLevelB);

    // if both nodeA and nodeB are leaves and their bounding boxes intersect
    if (bothLeaves && potentialCollision)  // check intersection on elements
    {
      for (int eleA : nodeA->indices) {
        for (int eleB : nodeB->indices) {
          elementsIntersect(eleA, eleB);
        }
      }
    }

    // if there is no further intersection from this node pair
    // then we go to the stack and get the next node pair
    if ((potentialCollision == false) || bothLeaves) {
      int childIndex = 0;  // index of this stack entry's node in its parent's childIDs array

      while (true) {
        childIndex = traversalStack.top().third;
        traversalStack.pop();
        if (traversalStack.empty()) {
          return;  // end traversal, return the function intersectionQuery()
        }
        StackEntry parentEntry = traversalStack.top();
        nodeAID = parentEntry.first;   // now it stores the parent of the previous NodeAID
        nodeBID = parentEntry.second;  // now it stores the parent of the previous NodeBID
                                       //        cout << "  back to parent " << nodeAID << " " << nodeBID << endl;
        PGO_ALOG(nodeAID >= 0 && nodeAID < BasicAlgorithms::sizei(nodes));
        PGO_ALOG(nodeBID >= 0 && nodeBID < BasicAlgorithms::sizei(otherTree.nodes));
        nodeA = &nodes[nodeAID];                  // now it refers to parent nodeA
        nodeB = &otherTree.nodes[nodeBID];        // now it refers to parent nodeB

        childIndex += 1;                          // get next child
        if (chooseWhichNodeToRecurse(nodeA, nodeB) == NODEA_CHOSEN) {
          if (childIndex < nodeA->numChildren())  // we found a valid child node, so break
          {
            nodeAID = nodeA->childrenIDs[childIndex];
            break;
          }
        }
        else {
          if (childIndex < nodeB->numChildren()) {
            nodeBID = nodeB->childrenIDs[childIndex];
            break;
          }
        }
      }

      // now nodeAID and nodeBID stores the node IDs for the next stack entry
      //      cout << "found child " << nodeAID << " " << nodeBID << " " << endl;
      traversalStack.push(triple<int, int, int>(nodeAID, nodeBID, childIndex));
    }
    else  // we can recurse on this node pair
    {
      int childIndex = 0;

      // push the new node
      if (chooseWhichNodeToRecurse(nodeA, nodeB) == NODEA_CHOSEN) {
        nodeAID = nodeA->childrenIDs[0];
      }
      else {
        nodeBID = nodeB->childrenIDs[0];
      }
      //      cout << "at a new entry " << nodeAID << " " << nodeBID << endl;
      traversalStack.emplace(nodeAID, nodeBID, childIndex);
    }
  }  // end traversal loop
}

void BoundingVolumeTreeBase::selfIntersectionQuery(ElementPairProcess elementsIntersect) const
{
  // we traverse a bounding volume
  // we compare bounding boxes of two nodes at the same level,
  //   or at immediately adjacent levels, with recursing first into the first node

  // define StackEntry: <// node from A, node from B, index of nodeA/B in its parent's childIDs array>
  // we make sure NodeA <= NodeB to avoid redundant comparison
  // the third value in StackEntry seems ambiguous: nodeA/B?
  // in fact, for a particular node A/B pair, we can determine which node to recurse by using the following function:
  // chooseWhichNodeToRecurse -> return NODEA_CHOSEN or false
  auto chooseWhichNodeToRecurse = [](const Node *nodeA, const Node *nodeB) {
    // we always try to recurse on B first
    // we recurse on A only if B cannot be recursed, due to B being a leaf or B too deep comparing to A
    return (nodeB->isLeaf() || ((nodeA->depth <= nodeB->depth) && nodeA->isLeaf() == false));
  };

  using StackEntry = triple<int, int, int>;
  std::stack<StackEntry> traversalStack;
  traversalStack.push(StackEntry(0, 0, 0));

  while (true) {
    StackEntry stackEntry = traversalStack.top();

    int nodeAID = stackEntry.first;
    int nodeBID = stackEntry.second;

    PGO_ALOG(nodeAID >= 0 && nodeAID < (int)nodes.size());
    PGO_ALOG(nodeBID >= 0 && nodeBID < (int)nodes.size());
    const auto *nodeA = &nodes[nodeAID];
    const auto *nodeB = &nodes[nodeBID];

    //    cout << "now " << nodeAID << " " << nodeBID << " " << traversalStack.size() << endl;

    bool potentialCollision = true;

    if (nodeAID != nodeBID) {
      potentialCollision = nodeA->bb.intersect(nodeB->bb);
    }

    // are the nodes on deepest level
    // negative sign by convention indicates leaf
    bool deepestLevelA = nodeA->isLeaf();
    bool deepestLevelB = nodeB->isLeaf();

    bool bothLeaves = deepestLevelA && deepestLevelB;  // if both nodes are leaves

    // if both nodeA and nodeB are leaves and their bounding boxes intersect
    if (bothLeaves && potentialCollision)  // check itnersections between elements
    {
      // printf("Checking all triangle pairs...\n");
      //  check all element pairs; note: nodeAp and nodeBp could be the same

      if (nodeA == nodeB) {
        const auto &indices = nodeA->indices;
        for (size_t i = 0; i < indices.size(); i++)
          for (size_t j = i + 1; j < indices.size(); j++) {
            elementsIntersect(indices[i], indices[j]);
          }
      }
      else {
        for (int eleA : nodeA->indices)
          for (int eleB : nodeB->indices) {
            if (eleA == eleB)
              continue;  // same node, continue
            elementsIntersect(eleA, eleB);
          }
      }
    }  // end checking itnersections between elements

    // modify childIndex to get the next available child index
    // return true if we can find one
    auto goToValidChild = [&](int childIndex) {
      bool hasChild = true;
      if (nodeA == nodeB) {
        // for intra-node tests, we store (i,j) indices of current child pair in upper and lower parts of childIndex;
        // we only traverse children with i <= j; this avoids checking pairs of nodes twice
        short int *childIndexB = (short int *)&childIndex;
        short int *childIndexA = childIndexB + 1;

        // find an available childIndexA
        for (; *childIndexA < nodeA->numChildren() && nodeA->childrenIDs[*childIndexA] < 0; (*childIndexB) = ++(*childIndexA)) {
        }

        // find an available childIndexB
        for (; *childIndexB < nodeB->numChildren() && nodeB->childrenIDs[*childIndexB] < 0; (*childIndexB)++) {
        }

        if (*childIndexB >= nodeB->numChildren()) {
          // find next available childIndexA
          for ((*childIndexA)++; *childIndexA < nodeA->numChildren() && nodeA->childrenIDs[*childIndexA] < 0; (*childIndexA)++) {
          }
          (*childIndexB) = (*childIndexA);
          hasChild = (*childIndexA < nodeA->numChildren());
        }
        else {
          PGO_ALOG(nodeA->childrenIDs[*childIndexA] >= 0);
        }

        if (hasChild) {
          nodeAID = nodeA->childrenIDs[*childIndexA];
          nodeBID = nodeB->childrenIDs[*childIndexB];
          traversalStack.push(triple<int, int, int>(nodeAID, nodeBID, childIndex));
        }
      }
      else  // A and B are different nodes
      {
        bool recurseOnFirstNode = chooseWhichNodeToRecurse(nodeA, nodeB);
        const auto recurseNode = (recurseOnFirstNode ? nodeA : nodeB);

        for (; childIndex < recurseNode->numChildren() && recurseNode->childrenIDs[childIndex] < 0; childIndex++) {
        }

        hasChild = (childIndex < recurseNode->numChildren());
        if (hasChild)  // go to this node pair by pushing it into the stack
        {
          if (recurseOnFirstNode)
            nodeAID = nodeA->childrenIDs[childIndex];
          else
            nodeBID = nodeB->childrenIDs[childIndex];
          //          cout << "Find a new node pair on different nodes: " << nodeAID << " " << nodeBID << " " << childIndex << endl;
          traversalStack.push(triple<int, int, int>(nodeAID, nodeBID, childIndex));
        }
      }

      return hasChild;
    };

    // if there is no further intersection from this node pair
    // then we go to the stack and get the next node pair
    if (potentialCollision == false || bothLeaves) {
      //      cout << "pop stack and go to next node pair" << endl;
      int parentChildIndex = stackEntry.third;
      while (true) {
        int childIndex = parentChildIndex;
        traversalStack.pop();
        if (traversalStack.empty()) {
          return;  // end traversal, return function::selfIntersectionQuery
        }

        triple<int, int, int> parentEntry = traversalStack.top();
        parentChildIndex = parentEntry.third;
        nodeAID = parentEntry.first;
        nodeBID = parentEntry.second;
        nodeA = &nodes[nodeAID];
        nodeB = &nodes[nodeBID];

        // try to increment childIndex
        //        cout << "found a new stack top " << nodeAID << " " << nodeBID << " previous childIndex " << childIndex << endl;
        if (nodeA == nodeB) {
          short int *childIndexB = (short int *)&childIndex;
          (*childIndexB)++;
        }
        else
          childIndex++;

        if (goToValidChild(childIndex))  // if we successfully find the next children pair to recursive from this pair as in parentEntry, then break
          break;                         // otherwise, we have to pop out more
      }
    }
    else {
      // cout << "recurse into next node pair" << endl;
      // Here, nodeA/B cannot both be leaves
      // recursively go down one level
      int childIndex = 0;
      bool success = goToValidChild(childIndex);
      PGO_ALOG(success);
    }
  }
}

void BoundingVolumeTreeBase::buildFromRoot(Partition divideNode, int maxDepth, int maxNumElementsPerNode)
{
  depth = 0;

  std::stack<int> nodeStack;
  nodeStack.push(0);
  while (nodeStack.empty() == false) {
    int nodeID = nodeStack.top();
    nodeStack.pop();

    int nodeDepth = nodes[nodeID].depth;
    if (nodeDepth > depth)
      depth = nodeDepth;

    int numElements = (int)nodes[nodeID].indices.size();

    // if there are fewer elements than the threshold value, or if max depth has been reached,
    // this node becomes a leaf which stores the element indices
    // max depth checking is necessary, otherwise the box might get split forever
    if ((numElements <= maxNumElementsPerNode) || (nodeDepth >= maxDepth))
      continue;

    std::vector<int> indices = std::move(nodes[nodeID].indices);
    PGO_ALOG(nodes[nodeID].indices.empty());  // non-leaf nodes should not store vertex indices

    std::vector<std::vector<int>> childIDs;
    std::vector<BoundingBox> subBBs;

    // PGO_ALOG(nodes[nodeID].bb.verifyBox());
    divideNode(indices, nodes[nodeID].bb, childIDs, subBBs);

    for (int i = 0; i < BasicAlgorithms::sizei(childIDs); i++) {
      if (childIDs[i].size() == 0)
        continue;
      int childrenID = int(nodes.size());
      nodes[nodeID].childrenIDs.push_back(childrenID);
      nodes.emplace_back(nodeDepth + 1, std::move(childIDs[i]));
      auto &childNode = nodes.back();
      childNode.bb = subBBs[i];
      // PGO_ALOG(subBBs[i].verifyBox());
      nodeStack.push(childrenID);
    }
  }
}

void BoundingVolumeTreeBase::boundingBoxOctreePartition(const std::vector<int> &elementList, const BoundingBox &bb,
  std::vector<std::vector<int>> &childIDs, std::vector<BoundingBox> &subBBs, const std::vector<BoundingBox> &elementBBs,
  BBFromElementList createBBFromElements, ElementBBIntersect intersectElementBB, bool shrink)
{
  Vec3d center = bb.center();

  BoundingBox spatialChildrenBBs[8];
  // partition node.bb into 8 children spatially
  bb.createChildBoundingBoxes(spatialChildrenBBs);

  childIDs.assign(8, {});
  subBBs.resize(8);

  // for each child, generate a list of vertices that intersect the child bounding box
  for (int eleID : elementList) {
    const auto &triBB = elementBBs[eleID];

    unsigned char lowerChildrenBBCode[3], higherChildrenBBCode[3];  // stores [0, 1]
    // lowerChildrenBBCCode stores which side bmin of triBB is in the childrenBBs on all dimenstions
    for (int dim = 0; dim < 3; dim++) {
      lowerChildrenBBCode[dim] = (triBB.bmin()[dim] > center[dim]);
      higherChildrenBBCode[dim] = (triBB.bmax()[dim] > center[dim]);
    }

    for (unsigned char ii = lowerChildrenBBCode[0]; ii <= higherChildrenBBCode[0]; ii++)
      for (unsigned char jj = lowerChildrenBBCode[1]; jj <= higherChildrenBBCode[1]; jj++)
        for (unsigned char kk = lowerChildrenBBCode[2]; kk <= higherChildrenBBCode[2]; kk++) {
          unsigned char code = ii + (jj << 1) + (kk << 2);
          // cout << spatialChildrenBBs[code] << endl;
          if (intersectElementBB(eleID, spatialChildrenBBs[code]))
            childIDs[code].push_back(eleID);
          //          PGO_ALOG(spatialChildrenBBs[code].intersect(triBB));
        }
  }

  for (int i = 0; i < 8; i++) {
    if (childIDs[i].size() > 0) {
      BoundingBox bbFromElements = createBBFromElements(childIDs[i]);
      if (shrink) {
        // we can use spatialChildrenBBs[i] as childNode.bb
        // but it might be larger than what triangles occupy
        // we can shrink this bb by doing intersection with allTriVtxBB
        subBBs[i] = spatialChildrenBBs[i].getIntersection(bbFromElements);
      }
      else
        subBBs[i] = spatialChildrenBBs[i];
    }
  }
}

void BoundingVolumeTreeBase::elementsInertiaPartition(const std::vector<int> &elements, const BoundingBox & /*bb*/,
  int numPartitions, std::vector<std::vector<int>> &childIDs,
  ElementMass getEleMass, ElementCenterOfMass getEleCOM, ElementInertiaTensor getEleInertia)
{
  Vec3d centerOfMass = asVec3d(0.0);
  Vec3d unweightedCenterOfMass = asVec3d(0.0);
  double totalMass = 0.0;
  for (int eleID : elements) {
    double eleMass = getEleMass(eleID);
    Vec3d eleCOM = getEleCOM(eleID);
    totalMass += eleMass;
    centerOfMass += eleMass * eleCOM;
    unweightedCenterOfMass += eleCOM;
  }

  Mat3d IT = Mat3d::Zero();
  if (totalMass > 0.0) {
    centerOfMass /= totalMass;
    for (int eleID : elements)
      IT += getEleInertia(eleID, centerOfMass);
  }
  else {
    unweightedCenterOfMass /= static_cast<double>(elements.size());
    centerOfMass = unweightedCenterOfMass;
    IT.setIdentity();
  }

  // do eig on IT; eigenvector with smallest eigenvalue is the normal of the partitioning plane
  Vec3d eig_val;
  Vec3d eig_vec[3];

  eigen_sym(IT, eig_val, eig_vec);

  const int partitionDim = 0;
  Vec3d partitioningPlaneNormal = eig_vec[partitionDim];  // coefficients of partitioning plane (normal vector)

  // compute data span along the normal vector
  double minZ = DBL_MAX;
  double maxZ = -DBL_MAX;
  std::vector<double> eleZValues(elements.size());
  for (int i = 0; i < BasicAlgorithms::sizei(elements); i++) {
    int eleID = elements[i];
    // compute the projection of the element center along the partitioning plane normal
    double zValue = partitioningPlaneNormal.dot(getEleCOM(eleID) - centerOfMass);
    eleZValues[i] = zValue;
    if (zValue < minZ)
      minZ = zValue;
    if (zValue > maxZ)
      maxZ = zValue;
  }

  std::vector<double> D(numPartitions - 1);  // partitioning thresholds
  // we create (numPartitions-1) equally spaced partitioning planes, partitioning the data range along the normal vector
  for (int i = 0; i < numPartitions - 1; i++) {
    D[i] = (minZ + (maxZ - minZ) * (i + 1.0) / numPartitions);
    // printf("  D[%d] = %G\n", i, D[i]);
  }

  // do the actual partition
  childIDs.assign(numPartitions, {});
  for (int i = 0; i < BasicAlgorithms::sizei(elements); i++) {
    int eleID = elements[i];
    double zValue = eleZValues[i];
    bool inserted = false;
    for (int i = 0; i < numPartitions - 1; i++) {
      // printf("DTriCom = %G\n", DTriCom);
      if (zValue < D[i])  // D[i] is ascending sequence with i
      {
        childIDs[i].push_back(eleID);
        inserted = true;
        break;
      }
    }
    if (!inserted)
      childIDs[numPartitions - 1].push_back(eleID);
  }
}

void BoundingVolumeTreeBase::elementsInertiaBinaryPartition(const std::vector<int> &elements, const BoundingBox & /*bb*/,
  std::vector<std::vector<int>> &childIDs, int childSizeBalanceThreshold,
  ElementMass getEleMass, ElementCenterOfMass getEleCOM, ElementInertiaTensor getEleInertia)
{
  Vec3d centerOfMass = asVec3d(0.0);
  Vec3d unweightedCenterOfMass = asVec3d(0.0);
  double totalMass = 0.0;
  for (int eleID : elements) {
    double eleMass = getEleMass(eleID);
    Vec3d eleCOM = getEleCOM(eleID);
    totalMass += eleMass;
    centerOfMass += eleMass * eleCOM;
    unweightedCenterOfMass += eleCOM;
  }

  Mat3d IT;
  IT.setZero();
  if (totalMass > 0.0) {
    centerOfMass /= totalMass;
    for (int eleID : elements)
      IT += getEleInertia(eleID, centerOfMass);
  }
  else {
    unweightedCenterOfMass /= static_cast<double>(elements.size());
    centerOfMass = unweightedCenterOfMass;
    IT.setIdentity();
  }

  // do eig on IT; eigenvector with smallest eigenvalue is the normal of the partitioning plane
  Vec3d eig_val;
  Vec3d eig_vec[3];

  eigen_sym(IT, eig_val, eig_vec);

  const int partitionDim = 0;
  Vec3d partitioningPlaneNormal = eig_vec[partitionDim];  // coefficients of partitioning plane (normal vector)

  // compute data span along the normal vector
  double minZ = DBL_MAX;
  double maxZ = -DBL_MAX;
  std::vector<double> eleZValues(elements.size());
  for (int i = 0; i < BasicAlgorithms::sizei(elements); i++) {
    int eleID = elements[i];
    // compute the projection of the element center along the partitioning plane normal
    double zValue = partitioningPlaneNormal.dot(getEleCOM(eleID) - centerOfMass);
    eleZValues[i] = zValue;
    if (zValue < minZ)
      minZ = zValue;
    if (zValue > maxZ)
      maxZ = zValue;
  }

  // we create a partitioning plane, partitioning the data range along the normal vector
  double D = (minZ + (maxZ - minZ) / 2.0);

  // do the actual partition
  childIDs.assign(2, {});
  for (int i = 0; i < BasicAlgorithms::sizei(elements); i++) {
    int eleID = elements[i];
    double zValue = eleZValues[i];

    if (zValue < D) {
      childIDs[0].push_back(eleID);
    }
    else
      childIDs[1].push_back(eleID);
  }

  if (std::abs(BasicAlgorithms::sizei(childIDs[0]) - BasicAlgorithms::sizei(childIDs[1])) > childSizeBalanceThreshold ||
    BasicAlgorithms::sizei(childIDs[0]) == 0 || BasicAlgorithms::sizei(childIDs[1]) == 0) {
    std::vector<std::pair<double, int>> zValueSortBuffer(elements.size());
    for (int i = 0; i < BasicAlgorithms::sizei(elements); i++)
      zValueSortBuffer[i] = std::make_pair(eleZValues[i], elements[i]);
    sort(zValueSortBuffer.begin(), zValueSortBuffer.end());
    auto middle = zValueSortBuffer.begin() + BasicAlgorithms::sizei(elements) / 2;
    childIDs[0].clear();
    childIDs[1].clear();
    for (auto it = zValueSortBuffer.begin(); it != middle; it++)
      childIDs[0].push_back(it->second);
    for (auto it = middle; it != zValueSortBuffer.end(); it++)
      childIDs[1].push_back(it->second);
  }
}
bool BoundingVolumeTreeBase::findElementRecursive(int nodeID, int elementID) const
{
  if (nodes[nodeID].isLeaf()) {
    for (int eleID : nodes[nodeID].indices)
      if (eleID == elementID)
        return true;
  }
  else {
    for (int childID : nodes[nodeID].childrenIDs) {
      if (childID < 0)
        continue;
      if (findElementRecursive(childID, elementID))
        return true;
    }
  }
  return false;
}

std::vector<int> BoundingVolumeTreeBase::getNodeCoveredElements(int nodeID) const
{
  std::vector<int> ret;
  getCoveredElementsRecursive(nodeID, ret);
  BasicAlgorithms::sortAndDeduplicate(ret);
  return ret;
}

void BoundingVolumeTreeBase::getCoveredElementsRecursive(int nodeID, std::vector<int> &elementIDs) const
{
  if (nodes[nodeID].isLeaf()) {
    BasicAlgorithms::insertRangeToEnd(elementIDs, nodes[nodeID].indices);
    return;
  }

  for (int childID : nodes[nodeID].childrenIDs) {
    if (childID < 0)
      continue;
    getCoveredElementsRecursive(childID, elementIDs);
  }
}

}  // namespace pgo::Mesh
