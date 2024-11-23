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

#include "windingNumberTree.h"
#include "triMeshGeo.h"
#include "triMeshNeighbor.h"
#include "basicAlgorithms.h"
#include "containerHelper.h"
#include "pgoLogging.h"

#include <stack>
#include <cassert>
#include <algorithm>

namespace pgo
{
namespace Mesh
{

static std::vector<Vec3i> buildBoundary(const std::vector<Vec3i> &triangles)
{
  std::vector<Vec3i> boundaryFaces;
  auto exEdges = getExteriorEdges(triangles);
  if (exEdges.size() != 0) {
    int v = exEdges[0][0];
    for (size_t edgeID = 0; edgeID < exEdges.size(); edgeID++) {
      boundaryFaces.emplace_back(v, exEdges[edgeID][0], exEdges[edgeID][1]);
    }
  }
  return boundaryFaces;
}

WindingNumberTree::Node::Node(int d, const TriMeshRef mesh):
  depth(d),
  // We don't build a bounding box out of all triMesh's positions because some vertices might not be used by triangles
  triangles(mesh.exportTriangles()), bouFaces(buildBoundary(triangles)), bb(mesh.computeTriangleBoundingBox())
{
}

void WindingNumberTree::build(const TriMeshRef mesh, int maxDepth, int maxNumTrianglesPerNode)
{
  clear();

  numMeshVertices = mesh.numVertices();
  numMeshTriangles = mesh.numTriangles();
  PGO_ALOG(numMeshTriangles > 0);

  //  vector<int> rootIndices(triMesh.numTriangles());
  //  iota(rootIndices.begin(), rootIndices.end(), 0);

  nodes.emplace_back(0, mesh);

  std::stack<int> nodeStack;
  nodeStack.push(0);
  while (nodeStack.empty() == false) {
    int nodeID = nodeStack.top();
    nodeStack.pop();

    int nodeDepth = nodes[nodeID].depth;
    if (nodeDepth > depth)
      depth = nodeDepth;

    // if there are fewer triangles than the threshold value, or cap value is larger than #triangles,
    // or if max depth has been reached,
    // this node becomes a leaf
    size_t numNodeTri = nodes[nodeID].triangles.size();
    if (numNodeTri <= (size_t)maxNumTrianglesPerNode || (nodes[nodeID].bouFaces.size() >= numNodeTri) || (nodeDepth >= maxDepth)) {
      continue;
    }

    // compute the longest side (0,1,2) of the node mesh
    int dim = nodes[nodeID].bb.longestSide().first;

    std::vector<double> triCen(numNodeTri);
    TriMeshRef nodeMesh = nodes[nodeID].meshInNode(mesh);
    for (int i = 0; i < nodeMesh.numTriangles(); i++)
      triCen[i] = nodeMesh.computeTriangleCentroid(i)[dim];

    double middle = BasicAlgorithms::median(triCen.begin(), triCen.end());

    std::vector<Vec3i> subTriangles[2];
    for (size_t i = 0; i < numNodeTri; i++) {
      if (triCen[i] <= middle)  // left side
        subTriangles[0].push_back(nodes[nodeID].triangles[i]);
      else
        subTriangles[1].push_back(nodes[nodeID].triangles[i]);
    }

    if (subTriangles[0].size() == 0 || subTriangles[1].size() == 0)
      continue;

    for (int i = 0; i < 2; i++) {
      int childrenID = int(nodes.size());
      nodes[nodeID].childrenIDs[i] = childrenID;
      nodes.emplace_back(nodeDepth + 1, TriMeshRef(mesh.numVertices(), mesh.positions(), subTriangles[i]));
      nodeStack.push(childrenID);
    }
  }
}

void WindingNumberTree::clear()
{
  nodes.clear();
  depth = 0;
  numMeshVertices = numMeshTriangles = 0;
}

double WindingNumberTree::windingNumber(const TriMeshRef mesh, const Vec3d &p) const
{
  PGO_ALOG(mesh.numVertices() == numMeshVertices);
  PGO_ALOG(mesh.numTriangles() == numMeshTriangles);
  PGO_ALOG(nodes.size() > 0);

  std::stack<int> nodeStack;
  nodeStack.push(0);
  double sum = 0.0;
  while (!nodeStack.empty()) {
    int nodeID = nodeStack.top();
    nodeStack.pop();
    const auto &node = nodes[nodeID];

    if (node.bb.checkInside(p)) {
      if (node.childrenIDs[0] >= 0) {
        nodeStack.push(node.childrenIDs[0]);
        nodeStack.push(node.childrenIDs[1]);
      }
      else {
        sum += node.meshInNode(mesh).computeWindingNumber(p);
      }
    }
    else  // p is outside bounding box of the node
    {
      if (node.bouFaces.size() < node.triangles.size()) {
        sum += node.boundaryInNode(mesh).computeWindingNumber(p);
      }
      else {
        sum += node.meshInNode(mesh).computeWindingNumber(p);
      }
    }
  }

  return sum;
}

void WindingNumberTree::print() const
{
  std::stack<int> nodeStack;
  nodeStack.push(0);
  while (!nodeStack.empty()) {
    int nodeID = nodeStack.top();
    nodeStack.pop();
    const auto &node = nodes[nodeID];
    for (int i = 0; i < node.depth; i++)
      std::cout << "  ";
    std::cout << "nodeID: " << nodeID << " #tri " << node.triangles.size() << " #bou " << node.bouFaces.size() << " ch " << node.childrenIDs[0] << " " << node.childrenIDs[1] << std::endl;
    if (node.childrenIDs[0] >= 0) {
      nodeStack.push(node.childrenIDs[0]);
      nodeStack.push(node.childrenIDs[1]);
    }
  }
}

}  // namespace Mesh
}  // namespace pgo
