/*
author: Bohan Wang
copyright to MIT, USC
*/

#pragma once

#include "EigenDef.h"
#include "triMeshGeo.h"

namespace MedialAxisRepresentation
{
namespace CGALUtilities
{
TriMeshGeo isotropicRemeshing(const TriMeshGeo &mesh, double targetEdgeLength, int numIter, double angleInDegree,
  std::vector<int> *constrainedVertices = nullptr, std::vector<std::pair<int, int>> *constrainedEdges = nullptr,
  std::vector<int> *faceSubset = nullptr);
TriMeshGeo smoothMesh(const TriMeshGeo &mesh, int numIter, double angleInDegree,
  std::vector<int> *constrainedVertices = nullptr, std::vector<std::pair<int, int>> *constrainedEdges = nullptr);
TriMeshGeo smoothShape(const TriMeshGeo &mesh, double time, int numIter);

TriMeshGeo refineMesh(const TriMeshGeo &mesh, double density);
TriMeshGeo refineSharpRegionOnMesh(const TriMeshGeo &mesh, double density, double angle);

TriMeshGeo fixOrientitation(const TriMeshGeo &mesh);

bool corefineAndComputeUnion(const TriMeshGeo &mesh1, const TriMeshGeo &mesh2, TriMeshGeo &unionMesh);
void corefineOnly(TriMeshGeo &mesh1, TriMeshGeo &mesh2, bool noModify1 = false, bool noModify2 = false, double edgeLengthThreshold = 1e-5);

TriMeshGeo triangulateHolePolyline(const std::vector<Vec3d> &polyline);
TriMeshGeo triangulateRefineFairHole(const TriMeshGeo &mesh);
TriMeshGeo triangulate(const std::vector<Vec3d> &vtx, const std::vector<std::vector<int>> &faces);

bool segmentPlaneIntersection(const Vec3d &planePnt, const Vec3d &planeNorm, const Vec3d &segP1, const Vec3d &segP2, Vec3d &intersection);

TriMeshGeo clipMesh(const TriMeshGeo &meshIn, const Vec3d &planeN, const Vec3d &planeP);
TriMeshGeo clipMesh(const TriMeshGeo &meshIn, const TriMeshGeo &volMesh);

TriMeshGeo cleanMesh(const TriMeshGeo &meshIn);

TriMeshGeo subdivideMesh(const TriMeshGeo &meshIn, int nIter, double smallSize = -1);

void convexHullMesh(const std::vector<TriMeshGeo> &meshes, TriMeshGeo &meshOut);

TriMeshGeo cdt2D(const std::vector<Vec3d> &polyline);

bool isSelfIntersected(const TriMeshGeo &meshIn);
bool isManifold(const TriMeshGeo &meshIn);
void getLargestCC(const TriMeshGeo &meshIn, TriMeshGeo &meshOut);

}  // namespace CGALUtilities
}  // namespace MedialAxisRepresentation