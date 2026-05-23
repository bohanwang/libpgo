#include "surfaceQuality.h"

#include "arrayRef.h"
#include "boundingVolumeTree.h"
#include "triMeshGeo.h"
#include "triMeshNeighbor.h"

#if MESH_QUALITY_CHECK_HAS_CGAL
#include "cgalInterface.h"
#endif

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mesh_quality_check
{
namespace
{

using pgo::Vec3d;
using pgo::Vec3i;

const char *selfIntersectionBackendName(SelfIntersectionBackend backend)
{
  switch (backend) {
    case SelfIntersectionBackend::CgalBool:
      return "cgal-bool";
    case SelfIntersectionBackend::ExactCount:
      return "exact-count";
  }

  return "unknown";
}

bool canUseCgalSelfIntersectionBackend(const SurfaceQualityReport &report)
{
  if (report.isManifold.has_value() && *report.isManifold == false)
    return false;

  if (report.boundaryOrExteriorEdges.has_value() && *report.boundaryOrExteriorEdges == 0 &&
    report.isWindingConsistent.has_value() && *report.isWindingConsistent == false)
    return false;

  return true;
}

bool hasDegenerateGeometry(const pgo::Mesh::TriMeshGeo &mesh, const Vec3i &tri)
{
  if (pgo::Mesh::isTriangleInvalid(tri))
    return true;

  const Vec3d e0 = mesh.pos(tri[1]) - mesh.pos(tri[0]);
  const Vec3d e1 = mesh.pos(tri[2]) - mesh.pos(tri[0]);
  return e0.cross(e1).squaredNorm() <= std::numeric_limits<double>::epsilon();
}

int countInvalidTriangles(const pgo::Mesh::TriMeshGeo &mesh)
{
  int count = 0;
  for (const Vec3i &tri : mesh.triangles()) {
    if (hasDegenerateGeometry(mesh, tri))
      ++count;
  }
  return count;
}

EdgeLengthStats computeEdgeLengthStats(const pgo::Mesh::TriMeshGeo &mesh)
{
  std::set<std::pair<int, int>> uniqueEdges;
  for (const Vec3i &tri : mesh.triangles()) {
    for (int i = 0; i < 3; ++i) {
      const int a = tri[i];
      const int b = tri[(i + 1) % 3];
      uniqueEdges.emplace(std::min(a, b), std::max(a, b));
    }
  }

  EdgeLengthStats stats;
  if (uniqueEdges.empty())
    return stats;

  stats.min = std::numeric_limits<double>::max();
  double sum = 0.0;
  for (const auto &edge : uniqueEdges) {
    const double length = (mesh.pos(edge.first) - mesh.pos(edge.second)).norm();
    stats.min = std::min(stats.min, length);
    stats.max = std::max(stats.max, length);
    sum += length;
  }

  stats.mean = sum / (double)uniqueEdges.size();
  return stats;
}

void fillBoundingBox(const pgo::Mesh::TriMeshGeo &mesh, SurfaceQualityReport &report)
{
  if (mesh.numTriangles() == 0)
    return;

  const pgo::Mesh::BoundingBox bbox = mesh.ref().computeTriangleBoundingBox();
  for (int i = 0; i < 3; ++i) {
    report.bboxMin[i] = bbox.bmin()[i];
    report.bboxMax[i] = bbox.bmax()[i];
  }
}

double computeEnclosedVolume(const pgo::Mesh::TriMeshGeo &mesh)
{
  double v6 = 0.0;
  for (int t = 0; t < mesh.numTriangles(); ++t) {
    const Vec3d &a = mesh.pos(t, 0);
    const Vec3d &b = mesh.pos(t, 1);
    const Vec3d &c = mesh.pos(t, 2);
    v6 += a.dot(b.cross(c));
  }
  return std::abs(v6) / 6.0;
}

void appendFailureReasons(SurfaceQualityReport &report, const SurfaceQualityOptions &options)
{
  if (report.vertices == 0)
    report.errors.push_back("Input triangle mesh contains no vertices");
  if (report.triangles == 0)
    report.errors.push_back("Input triangle mesh contains no triangles");
  if (options.failOnInvalidTriangles && report.invalidTriangles > 0)
    report.errors.push_back("Input triangle mesh contains invalid or degenerate triangles");
  if (options.failOnInvalidTriangles == false && report.invalidTriangles > 0)
    report.warnings.push_back("Input triangle mesh contains invalid or degenerate triangles");
  if (options.expectedComponents.has_value() && report.componentsByEdge.has_value() &&
    *report.componentsByEdge != *options.expectedComponents) {
    report.errors.push_back(
      "Expected " + std::to_string(*options.expectedComponents) +
      " edge-connected components, found " + std::to_string(*report.componentsByEdge));
  }
  if (report.selfIntersections.has_value() && *report.selfIntersections > 0)
    report.errors.push_back("Input triangle mesh has self intersections");
  if (report.isManifold.has_value() && *report.isManifold == false)
    report.errors.push_back("Input triangle mesh is not manifold");
  if (report.boundaryOrExteriorEdges.has_value() && *report.boundaryOrExteriorEdges > 0)
    report.errors.push_back("Input triangle mesh is not closed");
  if (report.isWindingConsistent.has_value() && *report.isWindingConsistent == false)
    report.errors.push_back("Input triangle mesh has inconsistent winding");
}

template<class T>
nlohmann::json optionalToJson(const std::optional<T> &value)
{
  if (value.has_value())
    return *value;
  return nullptr;
}

}  // namespace

SurfaceQualityReport checkSurfaceMesh(const std::string &inputMesh, const SurfaceQualityOptions &options)
{
  SurfaceQualityReport report;
  report.input = inputMesh;

  pgo::Mesh::TriMeshGeo mesh;
  if (mesh.load(inputMesh) != true) {
    report.errors.push_back("Failed to load triangle mesh: " + inputMesh);
    appendFailureReasons(report, options);
    return report;
  }

  report.vertices = mesh.numVertices();
  report.triangles = mesh.numTriangles();

  if (report.triangles > 0)
    fillBoundingBox(mesh, report);

  report.invalidTriangles = countInvalidTriangles(mesh);
  report.edgeLength = computeEdgeLengthStats(mesh);

  const bool basicGeometryPassed =
    report.vertices > 0 &&
    report.triangles > 0 &&
    (!options.failOnInvalidTriangles || report.invalidTriangles == 0);
  report.basicGeometryStatus = basicGeometryPassed ? "passed" : "failed";

  const bool canAnalyzeMesh = report.vertices > 0 && report.triangles > 0;
  if (canAnalyzeMesh) {
    const auto triangleRef = pgo::BasicAlgorithms::makeArrayRef(mesh.triangles());

    if (options.checkTopology) {
      const pgo::Mesh::TriangleEdgeConnectivityStats topology = pgo::Mesh::computeTriangleEdgeConnectivityStats(triangleRef);
      report.componentsByEdge = topology.componentsByEdge;
      report.componentTriangleCountsByEdge = topology.componentTriangleCountsByEdge;
      report.isManifold = topology.isManifold;
      report.boundaryOrExteriorEdges = topology.boundaryOrNonManifoldEdges;
      const bool expectedComponentsPassed =
        options.expectedComponents.has_value() == false ||
        *report.componentsByEdge == *options.expectedComponents;
      report.topologyStatus = (*report.isManifold && *report.boundaryOrExteriorEdges == 0 && expectedComponentsPassed) ? "passed" : "failed";
    }

    if (options.checkWinding) {
      report.orientedBoundaryOrExteriorEdges = (int)pgo::Mesh::getExteriorEdges(triangleRef).size();
      report.isWindingConsistent = pgo::Mesh::areTrianglesManifold(triangleRef) && *report.orientedBoundaryOrExteriorEdges == 0;
      report.windingStatus = *report.isWindingConsistent ? "passed" : "failed";
    }

    if (options.checkSelfIntersection && report.triangles <= options.selfIntersectionTriangleLimit) {
      SelfIntersectionBackend backend = options.selfIntersectionBackend;
      if (backend == SelfIntersectionBackend::CgalBool && canUseCgalSelfIntersectionBackend(report) == false) {
        backend = SelfIntersectionBackend::ExactCount;
        report.warnings.push_back("CGAL self-intersection backend requires a polygon mesh; falling back to exact-count");
      }

      report.selfIntersectionBackend = selfIntersectionBackendName(backend);

      if (backend == SelfIntersectionBackend::CgalBool) {
#if MESH_QUALITY_CHECK_HAS_CGAL
        const bool hasSelfIntersection = pgo::CGALInterface::isSelfIntersected(mesh);
        report.selfIntersections = hasSelfIntersection ? 1 : 0;
        report.selfIntersectionsExact = false;
#else
        throw std::runtime_error("CGAL self-intersection backend is unavailable in this build");
#endif
      }
      else if (backend == SelfIntersectionBackend::ExactCount) {
        pgo::Mesh::TriMeshBVTree bvTree;
        bvTree.buildByInertiaPartition(mesh.ref());

        std::vector<std::pair<int, int>> selfIntersections;
        bvTree.selfIntersectionExact(mesh.ref(), selfIntersections);
        report.selfIntersections = (int)selfIntersections.size();
        report.selfIntersectionsExact = true;
      }
      else {
        throw std::runtime_error("Unknown self-intersection backend");
      }

      report.selfIntersectionStatus = (*report.selfIntersections == 0) ? "passed" : "failed";
    }
    else if (options.checkSelfIntersection) {
      report.selfIntersectionStatus = "skipped";
      report.errors.push_back("Input triangle mesh is too large for self-intersection checking in this tool run");
    }

    const bool closedManifold =
      report.isManifold.has_value() && *report.isManifold &&
      report.boundaryOrExteriorEdges.has_value() && *report.boundaryOrExteriorEdges == 0;
    const bool windingKnownBad =
      report.isWindingConsistent.has_value() && *report.isWindingConsistent == false;
    if (closedManifold && !windingKnownBad)
      report.enclosedVolume = computeEnclosedVolume(mesh);
  }

  report.passed =
    report.basicGeometryStatus == "passed" &&
    (!options.checkTopology || report.topologyStatus == "passed") &&
    (!options.checkWinding || report.windingStatus == "passed") &&
    (!options.checkSelfIntersection || report.selfIntersectionStatus == "passed");

  appendFailureReasons(report, options);
  return report;
}

nlohmann::json surfaceQualityReportToJson(const SurfaceQualityReport &report)
{
  return nlohmann::json{
    { "type", report.type },
    { "input", report.input },
    { "passed", report.passed },
    { "check_status", {
        { "basic_geometry", report.basicGeometryStatus },
        { "topology", report.topologyStatus },
        { "winding", report.windingStatus },
        { "self_intersection", report.selfIntersectionStatus },
      } },
    { "vertices", report.vertices },
    { "triangles", report.triangles },
    { "components_by_edge", optionalToJson(report.componentsByEdge) },
    { "component_triangle_counts_by_edge", optionalToJson(report.componentTriangleCountsByEdge) },
    { "boundary_or_exterior_edges", optionalToJson(report.boundaryOrExteriorEdges) },
    { "is_manifold", optionalToJson(report.isManifold) },
    { "is_winding_consistent", optionalToJson(report.isWindingConsistent) },
    { "oriented_boundary_or_exterior_edges", optionalToJson(report.orientedBoundaryOrExteriorEdges) },
    { "self_intersections", optionalToJson(report.selfIntersections) },
    { "self_intersections_exact", optionalToJson(report.selfIntersectionsExact) },
    { "self_intersection_backend", optionalToJson(report.selfIntersectionBackend) },
    { "enclosed_volume", optionalToJson(report.enclosedVolume) },
    { "invalid_triangles", report.invalidTriangles },
    { "edge_length", {
        { "min", report.edgeLength.min },
        { "max", report.edgeLength.max },
        { "mean", report.edgeLength.mean },
      } },
    { "bbox", {
        { "min", report.bboxMin },
        { "max", report.bboxMax },
      } },
    { "errors", report.errors },
    { "warnings", report.warnings },
  };
}

}  // namespace mesh_quality_check
