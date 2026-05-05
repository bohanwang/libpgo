#include "boundingVolumeTree.h"
#include "cgalInterface.h"
#include "geogramInterface.h"
#include "pgoLogging.h"
#include "triMeshPseudoNormal.h"
#include "triMeshNeighbor.h"

#include <argparse/argparse.hpp>
#include <geogram/basic/process.h>

#include <cmath>

namespace
{

double triangleSignedVolume(const pgo::EigenSupport::V3d &a,
  const pgo::EigenSupport::V3d &b, const pgo::EigenSupport::V3d &c)
{
  return a.dot(b.cross(c)) / 6.0;
}

int orientNestedComponents(pgo::Mesh::TriMeshGeo &mesh)
{
  std::vector<int> componentTriangleCounts;
  const std::vector<int> componentIDs =
    pgo::Mesh::computeTriangleEdgeComponentIDs(mesh.ref().trianglesRef(), &componentTriangleCounts);
  if (componentTriangleCounts.size() <= 1)
    return 0;

  std::vector<double> componentVolumes(componentTriangleCounts.size(), 0.0);
  for (int triID = 0; triID < mesh.numTriangles(); ++triID) {
    const int componentID = componentIDs[triID];
    const pgo::Vec3i &tri = mesh.tri(triID);
    componentVolumes[componentID] += triangleSignedVolume(
      mesh.pos(tri[0]), mesh.pos(tri[1]), mesh.pos(tri[2]));
  }

  int outerComponent = 0;
  for (int componentID = 1; componentID < (int)componentVolumes.size(); ++componentID) {
    if (std::abs(componentVolumes[componentID]) > std::abs(componentVolumes[outerComponent]))
      outerComponent = componentID;
  }

  int flippedTriangles = 0;
  for (int triID = 0; triID < mesh.numTriangles(); ++triID) {
    const int componentID = componentIDs[triID];
    const double desiredSign = (componentID == outerComponent) ? 1.0 : -1.0;
    if (componentVolumes[componentID] * desiredSign < 0.0) {
      std::swap(mesh.tri(triID)[0], mesh.tri(triID)[1]);
      ++flippedTriangles;
    }
  }

  return flippedTriangles;
}

}  // namespace

int main(int argc, char *argv[])
{
  argparse::ArgumentParser program("Remesh surface");

  // git add subparser
  argparse::ArgumentParser cgal_smooth_cmd("cgal_smooth");
  cgal_smooth_cmd.add_description(
    "Use CGAL angle_and_area_smoothing() to remesh surface");
  cgal_smooth_cmd.add_argument("-i", "--input-mesh")
    .help("Input surface mesh filename")
    .required()
    .metavar("PATH");
  cgal_smooth_cmd.add_argument("-o", "--output-mesh")
    .help("Output surface mesh filename")
    .required()
    .metavar("PATH");
  cgal_smooth_cmd.add_argument("-s", "--sharp-edge-angle")
    .help("The angle threshold of the edge to be considered as a sharp edge")
    .default_value(180.0)
    .metavar("DEG")
    .scan<'g', double>();

  argparse::ArgumentParser cgal_iso_cmd("cgal_iso");
  cgal_iso_cmd.add_description(
    "Use CGAL isotropic_remeshing() to remesh surface");
  cgal_iso_cmd.add_argument("-i", "--input-mesh")
    .help("Input surface mesh filename")
    .required()
    .metavar("PATH");
  cgal_iso_cmd.add_argument("-o", "--output-mesh")
    .help("Output surface mesh filename")
    .required()
    .metavar("PATH");
  cgal_iso_cmd.add_argument("-l", "--edge-length")
    .help("The target edge length")
    .required()
    .metavar("FLOAT")
    .scan<'g', double>();
  cgal_iso_cmd.add_argument("-s", "--sharp-edge-angle")
    .help("The angle threshold of the edge to be considered as a sharp edge")
    .default_value(180.0)
    .metavar("DEG")
    .scan<'g', double>();
  cgal_iso_cmd.add_argument("--iterations")
    .help("Number of CGAL isotropic remeshing iterations")
    .default_value(10)
    .metavar("INT")
    .scan<'i', int>();

  argparse::ArgumentParser cgal_repair_self_intersections_cmd("cgal_repair_self_intersections");
  cgal_repair_self_intersections_cmd.add_description(
    "Use CGAL autorefine_and_remove_self_intersections() to repair local self-intersections");
  cgal_repair_self_intersections_cmd.add_argument("-i", "--input-mesh")
    .help("Input surface mesh filename")
    .required()
    .metavar("PATH");
  cgal_repair_self_intersections_cmd.add_argument("-o", "--output-mesh")
    .help("Output surface mesh filename")
    .required()
    .metavar("PATH");
  cgal_repair_self_intersections_cmd.add_argument("--method")
    .help("CGAL repair method: autorefine, autorefine-only, or remove")
    .default_value(std::string("autorefine"))
    .metavar("NAME");

  argparse::ArgumentParser cgal_simplify_cmd("cgal_simplify");
  cgal_simplify_cmd.add_description(
    "Use CGAL edge_collapse() to simplify surface");
  cgal_simplify_cmd.add_argument("-i", "--input-mesh")
    .help("Input surface mesh filename")
    .required()
    .metavar("PATH");
  cgal_simplify_cmd.add_argument("-o", "--output-mesh")
    .help("Output surface mesh filename")
    .required()
    .metavar("PATH");
  cgal_simplify_cmd.add_argument("-t", "--target-ratio")
    .help("The target ratio")
    .required()
    .metavar("FLOAT")
    .scan<'g', double>();

  argparse::ArgumentParser geogram_cmd("geogram");
  geogram_cmd.add_description("Use geogram remesher to remesh surface");
  geogram_cmd.add_argument("-i", "--input-mesh")
    .help("Input surface mesh filename")
    .required()
    .metavar("PATH");
  geogram_cmd.add_argument("-o", "--output-mesh")
    .help("Output surface mesh filename")
    .required()
    .metavar("PATH");
  geogram_cmd.add_argument("-t", "--target-num-vertices")
    .help("The target number vertices")
    .required()
    .metavar("FLOAT")
    .scan<'i', int>();
  geogram_cmd.add_argument("--threads")
    .help("Number of Geogram remeshing threads")
    .default_value(1)
    .metavar("INT")
    .scan<'i', int>();
  geogram_cmd.add_argument("--orient-nested-components")
    .help("Orient the largest closed component outward and all other components inward; useful for nested shell/cavity meshes before TetWild")
    .default_value(false)
    .implicit_value(true);

  program.add_subparser(cgal_smooth_cmd);
  program.add_subparser(cgal_iso_cmd);
  program.add_subparser(cgal_repair_self_intersections_cmd);
  program.add_subparser(cgal_simplify_cmd);
  program.add_subparser(geogram_cmd);

  try {
    program.parse_args(argc, argv);  // Example: ./main --color orange
  }
  catch (const std::exception &err) {
    std::cerr << err.what() << std::endl;
    std::cerr << program;
    return 1;
  }

  pgo::Logging::init();

  if (program.is_subcommand_used(cgal_smooth_cmd)) {
    pgo::Mesh::TriMeshGeo inputMesh;
    if (inputMesh.load(cgal_smooth_cmd.get<std::string>("--input-mesh")) != true)
      return 1;

    double angleThreshold = cgal_smooth_cmd.get<double>("--sharp-edge-angle");
    std::cout << "Sharp edge angle threshold: " << angleThreshold << std::endl;

    pgo::Mesh::TriMeshGeo meshOut =
      pgo::CGALInterface::smoothMesh(inputMesh, 10, angleThreshold);

    meshOut.save(cgal_smooth_cmd.get<std::string>("--output-mesh"));
  }
  else if (program.is_subcommand_used(cgal_iso_cmd)) {
    pgo::Mesh::TriMeshGeo inputMesh;
    if (inputMesh.load(cgal_iso_cmd.get<std::string>("--input-mesh")) != true)
      return 1;

    double avgLength = 0.0;
    for (int ti = 0; ti < inputMesh.numTriangles(); ti++) {
      for (int vi = 0; vi < 3; vi++) {
        pgo::EigenSupport::V3d diff =
          inputMesh.pos(ti, (vi + 1) % 3) - inputMesh.pos(ti, vi);
        avgLength += diff.norm();
      }
    }
    avgLength /= (inputMesh.numTriangles() * 3);
    std::cout << "Average length: " << avgLength << std::endl;

    double targetScale = cgal_iso_cmd.get<double>("--edge-length");
    double tgtLength = targetScale * avgLength;

    std::cout << "Target length: " << tgtLength << std::endl;

    double angleThreshold = cgal_iso_cmd.get<double>("--sharp-edge-angle");
    std::cout << "Sharp edge angle threshold: " << angleThreshold << std::endl;
    int iterations = cgal_iso_cmd.get<int>("--iterations");
    if (iterations < 0) {
      std::cerr << "--iterations must be non-negative" << std::endl;
      return 1;
    }
    std::cout << "Iterations: " << iterations << std::endl;

    pgo::Mesh::TriMeshGeo meshOut =
      pgo::CGALInterface::isotropicRemeshing(inputMesh, tgtLength, iterations, angleThreshold);

    meshOut.save(cgal_iso_cmd.get<std::string>("--output-mesh"));
  }
  else if (program.is_subcommand_used(cgal_repair_self_intersections_cmd)) {
    pgo::Mesh::TriMeshGeo inputMesh;
    if (inputMesh.load(cgal_repair_self_intersections_cmd.get<std::string>("--input-mesh")) != true)
      return 1;

    bool allFixed = false;
    const std::string method = cgal_repair_self_intersections_cmd.get<std::string>("--method");
    if (method != "autorefine" && method != "autorefine-only" && method != "remove") {
      std::cerr << "--method must be autorefine, autorefine-only, or remove" << std::endl;
      return 1;
    }
    std::cout << "Repair method: " << method << std::endl;

    pgo::Mesh::TriMeshGeo meshOut =
      pgo::CGALInterface::repairSelfIntersections(inputMesh, method, &allFixed);

    std::cout << "All self-intersections fixed: " << (allFixed ? "true" : "false") << std::endl;

    meshOut.save(cgal_repair_self_intersections_cmd.get<std::string>("--output-mesh"));
    if (allFixed == false)
      return 2;
  }
  else if (program.is_subcommand_used(cgal_simplify_cmd)) {
    pgo::Mesh::TriMeshGeo inputMesh;
    if (inputMesh.load(cgal_simplify_cmd.get<std::string>("--input-mesh")) != true)
      return 1;

    double tgtRatio = cgal_simplify_cmd.get<double>("--target-ratio");
    std::cout << "Target ratio: " << tgtRatio << std::endl;

    pgo::Mesh::TriMeshGeo meshOut =
      pgo::CGALInterface::simplifyMeshGH(inputMesh, "ptri", tgtRatio);

    meshOut.save(cgal_simplify_cmd.get<std::string>("--output-mesh"));
  }
  else if (program.is_subcommand_used(geogram_cmd)) {
    pgo::GeogramInterface::initGEO();

    int numThreads = geogram_cmd.get<int>("--threads");
    if (numThreads < 1) {
      std::cerr << "--threads must be positive" << std::endl;
      return 1;
    }
    GEO::Process::set_max_threads(numThreads);
    std::cout << "Geogram threads: " << numThreads << std::endl;

    int numInputTgtPts = geogram_cmd.get<int>("--target-num-vertices");
    int numTargetPoints = std::max(numInputTgtPts, 20);
    pgo::Mesh::TriMeshGeo outputMesh =
      pgo::GeogramInterface::remesh(geogram_cmd.get<std::string>("--input-mesh").c_str(), numTargetPoints, 0.8, -1.0);

    if constexpr (1) {
      pgo::Mesh::TriMeshBVTree outputMeshBVTree;
      outputMeshBVTree.buildByInertiaPartition(outputMesh);

      pgo::Mesh::TriMeshPseudoNormal outputMeshNormal;
      outputMeshNormal.buildPseudoNormals(outputMesh);

      pgo::Mesh::BoundingBox bb(outputMesh.positions());
      auto [idx, length] = bb.longestSide();

      pgo::EigenSupport::V3d pts[6] = {
        pgo::EigenSupport::V3d(length * 10, 0, 0) + bb.center(),
        pgo::EigenSupport::V3d(-length * 10, 0, 0) + bb.center(),

        pgo::EigenSupport::V3d(0, length * 10, 0) + bb.center(),
        pgo::EigenSupport::V3d(0, -length * 10, 0) + bb.center(),

        pgo::EigenSupport::V3d(0, 0, length * 10) + bb.center(),
        pgo::EigenSupport::V3d(0, 0, -length * 10) + bb.center(),
      };

      double maxRatio = 0;
      for (int viewi = 0; viewi < 6; viewi++) {
        int passThroughCount = 0;
        int correctNormalCount = 0;

        for (int vi = 0; vi < outputMesh.numVertices(); vi++) {
          pgo::EigenSupport::V3d rayDir = outputMesh.pos(vi) - pts[viewi];
          double dist = rayDir.norm();
          rayDir /= dist;

          pgo::EigenSupport::V3d segEnd = pts[viewi] + rayDir * (dist - 1e-6);
          pgo::EigenSupport::V3d segBeg = pts[viewi];

          if (outputMeshBVTree.hasLineSegmentIntersectionExact(outputMesh,
                segBeg, segEnd))
            continue;

          passThroughCount++;
          pgo::EigenSupport::V3d n = outputMeshNormal.vtxNormal(vi);
          if (n.dot(rayDir) < 0) {
            correctNormalCount++;
          }
        }

        std::cout << viewi << ":" << correctNormalCount << "/"
                  << passThroughCount << std::endl;
        maxRatio = std::max(maxRatio, correctNormalCount * 1.0 / passThroughCount);
      }

      if (maxRatio < 0.5) {
        std::cout << "Mesh is reversed." << std::endl;

        for (int ti = 0; ti < outputMesh.numTriangles(); ti++) {
          std::swap(outputMesh.tri(ti)[0], outputMesh.tri(ti)[1]);
        }
      }
    }

    if (geogram_cmd.get<bool>("--orient-nested-components")) {
      const int flippedTriangles = orientNestedComponents(outputMesh);
      std::cout << "Nested component orientation flipped triangles: " << flippedTriangles << std::endl;
    }

    pgo::Mesh::TriMeshGeo inputMesh;
    if (inputMesh.load(geogram_cmd.get<std::string>("--input-mesh")) != true)
      return 1;

    pgo::Mesh::TriMeshBVTree inputMeshBVTree;
    inputMeshBVTree.buildByInertiaPartition(inputMesh);
    for (int vi = 0; vi < outputMesh.numVertices(); vi++) {
      auto ret = inputMeshBVTree.closestTriangleQuery(inputMesh, outputMesh.pos(vi));
      outputMesh.pos(vi) = ret.closestPosition;
    }

    outputMesh.save(geogram_cmd.get<std::string>("--output-mesh"));
  }
  else {
    std::cerr << program;
    return 1;
  }

  return 0;
}
