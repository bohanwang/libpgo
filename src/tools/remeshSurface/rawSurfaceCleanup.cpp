#include "arrayRef.h"
#include "triMeshGeo.h"
#include "triMeshNeighbor.h"

#include <argparse/argparse.hpp>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

using pgo::Vec3d;
using pgo::Vec3i;

struct MeshStats
{
  int vertices = 0;
  int triangles = 0;
  int invalidTriangles = 0;
  int components = 0;
  int boundaryOrNonmanifoldEdges = 0;
  bool isManifold = false;
};

struct CleanupReport
{
  std::string input;
  std::string output;
  bool dryRun = false;
  int expectedComponents = -1;
  double shortEdgeThreshold = 0.0;
  int maxPasses = 0;
  int maxCollapses = 0;

  MeshStats before;
  MeshStats after;
  bool topologyPreserved = false;
  bool cleanupComplete = false;

  int attemptedDeletions = 0;
  int acceptedDeletions = 0;
  int attemptedCollapses = 0;
  int acceptedCollapses = 0;
  int rejectedByTopology = 0;
  int rejectedByInvalidCount = 0;
};

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

std::vector<int> collectInvalidTriangleIDs(const pgo::Mesh::TriMeshGeo &mesh)
{
  std::vector<int> invalidIDs;
  for (int triID = 0; triID < mesh.numTriangles(); ++triID) {
    if (hasDegenerateGeometry(mesh, mesh.tri(triID)))
      invalidIDs.push_back(triID);
  }
  return invalidIDs;
}

MeshStats computeStats(const pgo::Mesh::TriMeshGeo &mesh)
{
  MeshStats stats;
  stats.vertices = mesh.numVertices();
  stats.triangles = mesh.numTriangles();
  stats.invalidTriangles = countInvalidTriangles(mesh);

  if (mesh.numTriangles() > 0) {
    const pgo::Mesh::TriangleEdgeConnectivityStats topology =
      pgo::Mesh::computeTriangleEdgeConnectivityStats(pgo::BasicAlgorithms::makeArrayRef(mesh.triangles()));
    stats.components = topology.componentsByEdge;
    stats.boundaryOrNonmanifoldEdges = topology.boundaryOrNonManifoldEdges;
    stats.isManifold = topology.isManifold;
  }

  return stats;
}

bool topologyGatePasses(const MeshStats &stats, int expectedComponents)
{
  if (stats.triangles <= 0)
    return false;
  if (expectedComponents >= 0 && stats.components != expectedComponents)
    return false;
  return stats.isManifold && stats.boundaryOrNonmanifoldEdges == 0;
}

pgo::Mesh::TriMeshGeo removeTrianglesAndCompact(const pgo::Mesh::TriMeshGeo &mesh, const std::vector<int> &triIDs)
{
  return pgo::Mesh::removeIsolatedVertices(pgo::Mesh::removeTriangles(mesh.ref(), triIDs).ref());
}

pgo::Mesh::TriMeshGeo collapseEdgesAndCompact(const pgo::Mesh::TriMeshGeo &mesh, const std::vector<std::pair<int, int>> &edges)
{
  std::vector<int> parent(mesh.numVertices());
  std::iota(parent.begin(), parent.end(), 0);

  auto findRoot = [&](int v) {
    int root = v;
    while (parent[root] != root)
      root = parent[root];
    while (parent[v] != v) {
      const int next = parent[v];
      parent[v] = root;
      v = next;
    }
    return root;
  };

  auto unite = [&](int a, int b) {
    int ra = findRoot(a);
    int rb = findRoot(b);
    if (ra == rb)
      return;
    if (ra > rb)
      std::swap(ra, rb);
    parent[rb] = ra;
  };

  for (const auto &[v0, v1] : edges)
    unite(v0, v1);

  std::vector<Vec3d> positions(mesh.numVertices(), Vec3d(0.0, 0.0, 0.0));
  std::vector<int> counts(mesh.numVertices(), 0);
  for (int v = 0; v < mesh.numVertices(); ++v) {
    const int root = findRoot(v);
    positions[root] += mesh.pos(v);
    counts[root] += 1;
  }
  for (int v = 0; v < mesh.numVertices(); ++v) {
    if (counts[v] > 0)
      positions[v] /= (double)counts[v];
  }

  std::vector<Vec3i> keptTriangles;
  keptTriangles.reserve(mesh.numTriangles());
  for (Vec3i tri : mesh.triangles()) {
    for (int i = 0; i < 3; ++i)
      tri[i] = findRoot(tri[i]);
    if (pgo::Mesh::isTriangleInvalid(tri) == false)
      keptTriangles.push_back(tri);
  }

  return pgo::Mesh::removeIsolatedVertices(pgo::Mesh::TriMeshGeo(std::move(positions), std::move(keptTriangles)).ref());
}

std::vector<std::pair<int, int>> collectShortEdgesOnInvalidTriangles(const pgo::Mesh::TriMeshGeo &mesh, double shortEdgeThreshold)
{
  std::set<std::pair<int, int>> edgeSet;
  const double threshold2 = shortEdgeThreshold * shortEdgeThreshold;

  for (const int triID : collectInvalidTriangleIDs(mesh)) {
    const Vec3i &tri = mesh.tri(triID);
    for (int i = 0; i < 3; ++i) {
      const int a = tri[i];
      const int b = tri[(i + 1) % 3];
      if (a < 0 || b < 0 || a == b)
        continue;
      const double length2 = (mesh.pos(a) - mesh.pos(b)).squaredNorm();
      if (length2 <= threshold2)
        edgeSet.emplace(std::min(a, b), std::max(a, b));
    }
  }

  return std::vector<std::pair<int, int>>(edgeSet.begin(), edgeSet.end());
}

bool tryAcceptCandidate(
  const pgo::Mesh::TriMeshGeo &candidate,
  pgo::Mesh::TriMeshGeo &mesh,
  int expectedComponents,
  CleanupReport &report)
{
  const MeshStats currentStats = computeStats(mesh);
  const MeshStats candidateStats = computeStats(candidate);

  if (candidateStats.invalidTriangles >= currentStats.invalidTriangles) {
    ++report.rejectedByInvalidCount;
    return false;
  }

  if (topologyGatePasses(candidateStats, expectedComponents) == false) {
    ++report.rejectedByTopology;
    return false;
  }

  mesh = candidate;
  return true;
}

void cleanupMesh(pgo::Mesh::TriMeshGeo &mesh, CleanupReport &report)
{
  for (int pass = 0; pass < report.maxPasses; ++pass) {
    bool changed = false;

    const std::vector<int> invalidIDs = collectInvalidTriangleIDs(mesh);
    if (invalidIDs.empty())
      break;

    report.attemptedDeletions += (int)invalidIDs.size();
    {
      pgo::Mesh::TriMeshGeo candidate = removeTrianglesAndCompact(mesh, invalidIDs);
      if (tryAcceptCandidate(candidate, mesh, report.expectedComponents, report)) {
        report.acceptedDeletions += (int)invalidIDs.size();
        changed = true;
      }
    }

    if (report.acceptedCollapses < report.maxCollapses) {
      const std::vector<std::pair<int, int>> shortEdges = collectShortEdgesOnInvalidTriangles(mesh, report.shortEdgeThreshold);
      const int collapseCount = std::min((int)shortEdges.size(), report.maxCollapses - report.acceptedCollapses);
      report.attemptedCollapses += collapseCount;

      if (collapseCount > 0) {
        std::vector<std::pair<int, int>> collapseEdges(shortEdges.begin(), shortEdges.begin() + collapseCount);
        pgo::Mesh::TriMeshGeo candidate = collapseEdgesAndCompact(mesh, collapseEdges);
        if (tryAcceptCandidate(candidate, mesh, report.expectedComponents, report)) {
          report.acceptedCollapses += collapseCount;
          changed = true;
        }
      }
    }

    if (changed == false)
      break;
  }
}

void writeJsonReport(const std::string &path, const nlohmann::json &report)
{
  const std::filesystem::path outputPath(path);
  if (outputPath.has_parent_path())
    std::filesystem::create_directories(outputPath.parent_path());

  std::ofstream out(outputPath);
  if (out.is_open() == false)
    throw std::runtime_error("Failed to open JSON report for writing: " + path);

  out << report.dump(2) << '\n';
}

nlohmann::json statsToJson(const MeshStats &stats)
{
  return {
    { "vertices", stats.vertices },
    { "triangles", stats.triangles },
    { "invalid_triangles", stats.invalidTriangles },
    { "components", stats.components },
    { "boundary_or_nonmanifold_edges", stats.boundaryOrNonmanifoldEdges },
    { "is_manifold", stats.isManifold },
  };
}

nlohmann::json reportToJson(const CleanupReport &report)
{
  return {
    { "type", "raw_surface_cleanup" },
    { "input", report.input },
    { "output", report.output },
    { "dry_run", report.dryRun },
    { "expected_components", report.expectedComponents >= 0 ? nlohmann::json(report.expectedComponents) : nlohmann::json(nullptr) },
    { "short_edge_threshold", report.shortEdgeThreshold },
    { "max_passes", report.maxPasses },
    { "max_collapses", report.maxCollapses },
    { "before", statsToJson(report.before) },
    { "after", statsToJson(report.after) },
    { "vertices_before", report.before.vertices },
    { "vertices_after", report.after.vertices },
    { "triangles_before", report.before.triangles },
    { "triangles_after", report.after.triangles },
    { "invalid_triangles_before", report.before.invalidTriangles },
    { "invalid_triangles_after", report.after.invalidTriangles },
    { "components_before", report.before.components },
    { "components_after", report.after.components },
    { "boundary_or_nonmanifold_edges_before", report.before.boundaryOrNonmanifoldEdges },
    { "boundary_or_nonmanifold_edges_after", report.after.boundaryOrNonmanifoldEdges },
    { "is_manifold_before", report.before.isManifold },
    { "is_manifold_after", report.after.isManifold },
    { "topology_preserved", report.topologyPreserved },
    { "cleanup_complete", report.cleanupComplete },
    { "attempted_deletions", report.attemptedDeletions },
    { "accepted_deletions", report.acceptedDeletions },
    { "attempted_collapses", report.attemptedCollapses },
    { "accepted_collapses", report.acceptedCollapses },
    { "rejected_by_topology", report.rejectedByTopology },
    { "rejected_by_invalid_count", report.rejectedByInvalidCount },
  };
}

}  // namespace

int main(int argc, char *argv[])
{
  argparse::ArgumentParser program("rawSurfaceCleanup");
  program.add_description("Conservative diagnostic cleanup for raw triangle surfaces");
  program.add_argument("--input")
    .help("Input surface mesh filename")
    .required()
    .metavar("PATH");
  program.add_argument("--output")
    .help("Output surface mesh filename")
    .required()
    .metavar("PATH");
  program.add_argument("--json")
    .help("Output JSON report filename")
    .required()
    .metavar("PATH");
  program.add_argument("--expected-components")
    .help("Expected edge-connected component count; omit to preserve the input count")
    .default_value(-1)
    .metavar("INT")
    .scan<'i', int>();
  program.add_argument("--short-edge-threshold")
    .help("Only existing edges at or below this length are eligible for collapse")
    .default_value(1e-5)
    .metavar("FLOAT")
    .scan<'g', double>();
  program.add_argument("--max-passes")
    .help("Maximum cleanup passes")
    .default_value(3)
    .metavar("INT")
    .scan<'i', int>();
  program.add_argument("--max-collapses")
    .help("Maximum accepted edge collapses")
    .default_value(10000)
    .metavar("INT")
    .scan<'i', int>();
  program.add_argument("--dry-run")
    .help("Write the report without writing the cleaned mesh")
    .default_value(false)
    .implicit_value(true);

  try {
    program.parse_args(argc, argv);
  }
  catch (const std::exception &err) {
    std::cerr << err.what() << std::endl;
    std::cerr << program;
    return 1;
  }

  try {
    CleanupReport report;
    report.input = program.get<std::string>("--input");
    report.output = program.get<std::string>("--output");
    const std::string jsonPath = program.get<std::string>("--json");
    report.expectedComponents = program.get<int>("--expected-components");
    report.shortEdgeThreshold = program.get<double>("--short-edge-threshold");
    report.maxPasses = program.get<int>("--max-passes");
    report.maxCollapses = program.get<int>("--max-collapses");
    report.dryRun = program.get<bool>("--dry-run");

    if (report.expectedComponents < -1)
      throw std::runtime_error("--expected-components must be non-negative");
    if (report.shortEdgeThreshold < 0.0)
      throw std::runtime_error("--short-edge-threshold must be non-negative");
    if (report.maxPasses < 0)
      throw std::runtime_error("--max-passes must be non-negative");
    if (report.maxCollapses < 0)
      throw std::runtime_error("--max-collapses must be non-negative");

    pgo::Mesh::TriMeshGeo mesh;
    if (mesh.load(report.input) == false)
      throw std::runtime_error("Failed to load triangle mesh: " + report.input);

    report.before = computeStats(mesh);
    if (report.expectedComponents < 0)
      report.expectedComponents = report.before.components;

    cleanupMesh(mesh, report);

    report.after = computeStats(mesh);
    report.topologyPreserved = topologyGatePasses(report.after, report.expectedComponents);
    report.cleanupComplete = report.after.invalidTriangles == 0;

    if (report.dryRun == false) {
      const std::filesystem::path outputPath(report.output);
      if (outputPath.has_parent_path())
        std::filesystem::create_directories(outputPath.parent_path());
      if (mesh.save(report.output) == false)
        throw std::runtime_error("Failed to write cleaned mesh: " + report.output);
    }

    writeJsonReport(jsonPath, reportToJson(report));
    return 0;
  }
  catch (const std::exception &err) {
    std::cerr << err.what() << std::endl;
    return 1;
  }
}
