/*
author: Bohan Wang
copyright to MIT, USC
*/

#pragma once

#include "EigenDef.h"
#include "triMeshGeo.h"

namespace pgo
{
namespace CGALInterface
{
Mesh::TriMeshGeo isotropicRemeshing(const Mesh::TriMeshGeo &mesh, double targetEdgeLength, int numIter, double angleInDegree,
  std::vector<int> *constrainedVertices = nullptr, std::vector<std::pair<int, int>> *constrainedEdges = nullptr,
  std::vector<int> *faceSubset = nullptr);
Mesh::TriMeshGeo smoothMesh(const Mesh::TriMeshGeo &mesh, int numIter, double angleInDegree,
  std::vector<int> *constrainedVertices = nullptr, std::vector<std::pair<int, int>> *constrainedEdges = nullptr);
Mesh::TriMeshGeo smoothShape(const Mesh::TriMeshGeo &mesh, double time, int numIter);

Mesh::TriMeshGeo simplifyMesh(const Mesh::TriMeshGeo &meshIn, double edgeStoppingRatio);
Mesh::TriMeshGeo simplifyMeshGH(const Mesh::TriMeshGeo &meshIn, const std::string &method, double edgeStoppingRatio);

Mesh::TriMeshGeo refineMesh(const Mesh::TriMeshGeo &mesh, double density);
Mesh::TriMeshGeo refineSharpRegionOnMesh(const Mesh::TriMeshGeo &mesh, double density, double angle);
Mesh::TriMeshGeo subdivideMesh(const Mesh::TriMeshGeo &meshIn, int nIter, double smallSize = -1);

bool corefineAndComputeUnion(const Mesh::TriMeshGeo &mesh1, const Mesh::TriMeshGeo &mesh2, Mesh::TriMeshGeo &unionMesh);
void corefineOnly(Mesh::TriMeshGeo &mesh1, Mesh::TriMeshGeo &mesh2, bool noModify1 = false, bool noModify2 = false, double edgeLengthThreshold = 1e-5);

Mesh::TriMeshGeo triangulateHolePolyline(const std::vector<Vec3d> &polyline);
Mesh::TriMeshGeo triangulateRefineFairHole(const Mesh::TriMeshGeo &mesh);
Mesh::TriMeshGeo triangulate(const std::vector<Vec3d> &vtx, const std::vector<std::vector<int>> &faces);

Mesh::TriMeshGeo clipMesh(const Mesh::TriMeshGeo &meshIn, const Vec3d &planeN, const Vec3d &planeP);
Mesh::TriMeshGeo clipMesh(const Mesh::TriMeshGeo &meshIn, const Mesh::TriMeshGeo &volMesh);

void segmentMesh(const Mesh::TriMeshGeo &meshIn, int nClusters, std::vector<int> &classID);

void convexHullMesh(const std::vector<Mesh::TriMeshGeo> &meshes, Mesh::TriMeshGeo &meshOut);

bool isSelfIntersected(const Mesh::TriMeshGeo &meshIn);
bool isManifold(const Mesh::TriMeshGeo &meshIn);
void getLargestCC(const Mesh::TriMeshGeo &meshIn, Mesh::TriMeshGeo &meshOut);

}  // namespace CGALInterface
}  // namespace pgo