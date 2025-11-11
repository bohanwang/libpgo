#include "tetgenInterface.h"
#include "pgoLogging.h"
#include "triMeshGeo.h"
#include "tetMesh.h"

#include <argparse/argparse.hpp>

int main(int argc, char *argv[])
{
  argparse::ArgumentParser program("Tetrahedralization surface");

  // git add subparser
  argparse::ArgumentParser tetgen_cmd("tetgen");
  tetgen_cmd.add_description(
    "Use tetgen to tetrahedralize surface");
  tetgen_cmd.add_argument("-i", "--input-mesh")
    .help("Input surface mesh filename")
    .required()
    .metavar("PATH");
  tetgen_cmd.add_argument("-o", "--output-mesh")
    .help("Output surface mesh filename")
    .required()
    .metavar("PATH");
  tetgen_cmd.add_argument("-c", "--command")
    .help("The angle threshold of the edge to be considered as a sharp edge")
    .required()
    .metavar("CMD");

  program.add_subparser(tetgen_cmd);

  try {
    program.parse_args(argc, argv);  // Example: ./main --color orange
  }
  catch (const std::exception &err) {
    std::cerr << err.what() << std::endl;
    std::cerr << program;
    return 1;
  }

  pgo::Logging::init();

  if (program.is_subcommand_used(tetgen_cmd)) {
    pgo::Mesh::TriMeshGeo inputMesh;
    if (inputMesh.load(tetgen_cmd.get<std::string>("--input-mesh")) != true)
      return 1;

    std::string tetgenCommand = tetgen_cmd.get<std::string>("--command");
    std::cout << "Tetgen command: " << tetgenCommand << std::endl;

    pgo::EigenSupport::MXd V;
    pgo::EigenSupport::MXi T;
    pgo::TetgenInterface::computeTetMesh(inputMesh, tetgenCommand, V, T);

    std::vector<pgo::Vec3d> positions;
    std::vector<pgo::Vec4i> tets;

    for (int i = 0; i < V.rows(); ++i) {
      positions.push_back(pgo::Vec3d(V(i, 0), V(i, 1), V(i, 2)));
    }

    for (int i = 0; i < T.rows(); ++i) {
      tets.push_back(pgo::Vec4i(T(i, 0), T(i, 1), T(i, 2), T(i, 3)));
    }

    pgo::VolumetricMeshes::TetMesh tetMesh(positions, tets);
    std::string outputMeshFilename = tetgen_cmd.get<std::string>("--output-mesh");
    if (tetMesh.save(outputMeshFilename.c_str()) != 0) {
      SPDLOG_LOGGER_ERROR(pgo::Logging::lgr(), "Failed to save tet mesh to file {}\n", outputMeshFilename);
      return 1;
    }
  }

  return 0;
}