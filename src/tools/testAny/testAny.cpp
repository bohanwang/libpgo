#include "pgo_c.h"
#include "pgoLogging.h"
#include <vector>
#include <fstream>

using namespace pgo;

int main(int argc, char const *argv[])
{
  pgo::Logging::init();
  std::string tetMeshFile = argv[1];
  std::string fixedVtxFile = argv[2];
  std::string saveFile = argv[3];
  // std::vector<int> fixedVertices;

  // std::ifstream ifs(fixedVtxFile);
  // int v;
  // while (ifs >> v)
  // {
  //   fixedVertices.push_back(v);
  // }

  // /////
  // pgoTetMeshStructHandle tetmesh = pgo_create_tetmesh_from_file(tetMeshFile.c_str());

  // int n = pgo_tetmesh_get_num_vertices(tetmesh);
  // int nele = pgo_tetmesh_get_num_tets(tetmesh);

  // ///// 1. quasi-static simulation on initShapeMeshFile
  // std::vector<double> xStaticEqRes(n * 3 + nele * 6);
  // pgo_create_quastic_static_sim(tetmesh, fixedVertices.data(), (int)fixedVertices.size(), xStaticEqRes.data(), nullptr, true);

  // pgoTetMeshStructHandle tetMeshStaticEqRes = pgo_tetmesh_update_vertices(tetmesh, xStaticEqRes.data());
  // pgo_save_tetmesh_to_file(tetMeshStaticEqRes, saveFile.c_str());
  double stepSize = 0.2;
  int verbose = 0;
  pgo_inverse_plasticity_opt(tetMeshFile.c_str(), fixedVtxFile.c_str(), saveFile.c_str(), stepSize, verbose);

  return 0;
}