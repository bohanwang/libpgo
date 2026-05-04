#pragma once

#include <nlohmann/json_fwd.hpp>

#include <array>
#include <optional>
#include <string>
#include <vector>

namespace mesh_quality_check
{

struct EdgeLengthStats
{
  double min = 0.0;
  double max = 0.0;
  double mean = 0.0;
};

struct SurfaceQualityReport
{
  std::string type = "surface";
  std::string input;
  bool passed = false;
  std::string basicGeometryStatus = "failed";
  std::string topologyStatus = "skipped";
  std::string windingStatus = "skipped";
  std::string selfIntersectionStatus = "skipped";
  int vertices = 0;
  int triangles = 0;
  std::optional<int> componentsByEdge;
  std::optional<std::vector<int>> componentTriangleCountsByEdge;
  std::optional<int> boundaryOrExteriorEdges;
  std::optional<bool> isManifold;
  std::optional<bool> isWindingConsistent;
  std::optional<int> orientedBoundaryOrExteriorEdges;
  std::optional<int> selfIntersections;
  int invalidTriangles = 0;
  EdgeLengthStats edgeLength;
  std::array<double, 3> bboxMin{ { 0.0, 0.0, 0.0 } };
  std::array<double, 3> bboxMax{ { 0.0, 0.0, 0.0 } };
  std::vector<std::string> errors;
  std::vector<std::string> warnings;
};

struct SurfaceQualityOptions
{
  bool checkTopology = true;
  bool checkWinding = true;
  bool checkSelfIntersection = true;
  bool failOnInvalidTriangles = true;
  std::optional<int> expectedComponents;
  int selfIntersectionTriangleLimit = 200000;
};

SurfaceQualityReport checkSurfaceMesh(const std::string &inputMesh, const SurfaceQualityOptions &options = {});
nlohmann::json surfaceQualityReportToJson(const SurfaceQualityReport &report);

}  // namespace mesh_quality_check
