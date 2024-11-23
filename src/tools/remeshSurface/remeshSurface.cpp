#include "boundingVolumeTree.h"
#include "cgalInterface.h"
#include "geogramInterface.h"
#include "pgoLogging.h"
#include "triMeshPseudoNormal.h"

#include <argparse/argparse.hpp>

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

  program.add_subparser(cgal_smooth_cmd);
  program.add_subparser(cgal_iso_cmd);
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

    pgo::Mesh::TriMeshGeo meshOut =
      pgo::CGALInterface::isotropicRemeshing(inputMesh, tgtLength, 10, angleThreshold);

    meshOut.save(cgal_iso_cmd.get<std::string>("--output-mesh"));
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