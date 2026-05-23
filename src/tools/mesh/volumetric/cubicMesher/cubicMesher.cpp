#include "cubicMesherIO.h"
#include "pgoLogging.h"
#include "triangleMeshVoxelizer.h"

#include <argparse/argparse.hpp>

#include <iostream>
#include <memory>
#include <string>

int main(int argc, char *argv[])
{
  argparse::ArgumentParser program("cubicMesher");
  program.add_description("Voxelize a closed triangle mesh into a cubic .veg mesh");
  program.add_argument("--input-mesh")
    .help("Input triangle mesh filename (.obj)")
    .required()
    .metavar("PATH");
  program.add_argument("--resolution")
    .help("Number of cubic elements along the shortest input AABB side")
    .required()
    .scan<'i', int>()
    .metavar("N");
  program.add_argument("--output-mesh")
    .help("Output cubic mesh filename (.veg)")
    .required()
    .metavar("PATH");
  program.add_argument("--output-surface")
    .help("Optional output surface mesh filename (.obj)")
    .metavar("PATH");
  program.add_argument("--E")
    .help("Young's modulus")
    .required()
    .scan<'g', double>();
  program.add_argument("--nu")
    .help("Poisson's ratio")
    .required()
    .scan<'g', double>();
  program.add_argument("--density")
    .help("Material density")
    .required()
    .scan<'g', double>();

  try {
    program.parse_args(argc, argv);
  }
  catch (const std::exception &err) {
    std::cerr << err.what() << std::endl;
    std::cerr << program;
    return 1;
  }

  pgo::Logging::init();

  try {
    cubic_mesher::TriangleMeshVoxelizerOptions options;
    options.inputMesh = program.get<std::string>("--input-mesh");
    options.resolution = program.get<int>("--resolution");
    options.E = program.get<double>("--E");
    options.nu = program.get<double>("--nu");
    options.density = program.get<double>("--density");

    std::unique_ptr<pgo::VolumetricMeshes::CubicMesh> cubicMesh = cubic_mesher::createTriangleMeshCubicMesh(options);
    cubic_mesher::saveCubicMesh(*cubicMesh, program.get<std::string>("--output-mesh"));

    if (program.is_used("--output-surface"))
      cubic_mesher::writeSurfaceMesh(*cubicMesh, program.get<std::string>("--output-surface"));

    std::cout << "Generated cubic mesh with " << cubicMesh->getNumVertices() << " vertices and "
              << cubicMesh->getNumElements() << " elements." << std::endl;
  }
  catch (const std::exception &err) {
    std::cerr << err.what() << std::endl;
    return 1;
  }

  return 0;
}
