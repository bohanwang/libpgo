#include "tetMesherBackend.h"

#include <floattetwild/FloatTetwild.h>
#include <floattetwild/Logger.hpp>
#include <floattetwild/MeshIO.hpp>

#include <geogram/basic/common.h>
#include <geogram/mesh/mesh.h>

#include <Eigen/Dense>

#include <array>
#include <cstdio>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

namespace tet_mesher
{
namespace
{

void initializeTetWild()
{
  static std::once_flag initFlag;
  std::call_once(initFlag, []() {
    GEO::initialize();
  });
}

void validateTetwildOptions(const TetwildOptions &options)
{
  if (options.hasLa) {
    if (options.la <= 0.0)
      throw std::runtime_error("tetwild.la must be positive");
  }
  else if (options.lr <= 0.0) {
    throw std::runtime_error("tetwild.lr must be positive");
  }

  if (options.epsr <= 0.0)
    throw std::runtime_error("tetwild.epsr must be positive");

  if (options.stopEnergy <= 0.0)
    throw std::runtime_error("tetwild.stop_energy must be positive");

  if (options.maxThreads < 0)
    throw std::runtime_error("tetwild.max_threads must be non-negative");
}

void cleanupTetwildScratchFiles(const std::string &scratchPrefix)
{
  const std::array<std::string, 5> suffixes = {
    "_pgo_tmp_tracked_surface.stl",
    "_pgo_tmp_opt.stl",
    "_pgo_tmp_cutting.stl",
    "_pgo_tmp_simplify.off",
    "_pgo_tmp_b_vs.xyz",
  };

  for (const std::string &suffix : suffixes)
    std::remove((scratchPrefix + suffix).c_str());
}

std::unique_ptr<pgo::VolumetricMeshes::TetMesh> convertTetwildOutput(const Eigen::MatrixXd &V, const Eigen::MatrixXi &T)
{
  if (V.cols() != 3)
    throw std::runtime_error("fTetWild returned a vertex matrix that does not have 3 columns");

  if (T.cols() != 4)
    throw std::runtime_error("fTetWild returned an element matrix that does not have 4 columns");

  std::vector<pgo::Vec3d> positions;
  positions.reserve(V.rows());
  for (int i = 0; i < V.rows(); ++i)
    positions.emplace_back(V(i, 0), V(i, 1), V(i, 2));

  std::vector<pgo::Vec4i> tets;
  tets.reserve(T.rows());
  for (int i = 0; i < T.rows(); ++i) {
    for (int j = 0; j < 4; ++j) {
      if (T(i, j) < 0 || T(i, j) >= V.rows())
        throw std::runtime_error("fTetWild returned a tetrahedron with an out-of-range vertex index");
    }

    tets.emplace_back(T(i, 0), T(i, 1), T(i, 2), T(i, 3));
  }

  auto tetMesh = std::make_unique<pgo::VolumetricMeshes::TetMesh>(positions, tets);
  tetMesh->orient();
  return tetMesh;
}

}  // namespace

std::unique_ptr<pgo::VolumetricMeshes::TetMesh> generateTetwildMesh(const TetwildOptions &options)
{
  validateTetwildOptions(options);
  initializeTetWild();

  floatTetWild::Logger::init(options.common.quiet == false, "");

  GEO::Mesh sfMesh;
  std::vector<floatTetWild::Vector3> inputVertices;
  std::vector<floatTetWild::Vector3i> inputFaces;
  std::vector<int> inputTags;
  if (floatTetWild::MeshIO::load_mesh(options.common.inputMesh, inputVertices, inputFaces, sfMesh, inputTags) == false)
    throw std::runtime_error("fTetWild failed to load input surface mesh " + options.common.inputMesh);

  if (inputVertices.empty() || inputFaces.empty())
    throw std::runtime_error("fTetWild input surface mesh is empty: " + options.common.inputMesh);

  floatTetWild::Parameters params;
  params.input_path = options.common.inputMesh;
  params.output_path = options.common.outputMesh + ".ftetwild";
  params.postfix = "pgo_tmp";
  params.is_quiet = options.common.quiet;
  params.log_level = options.common.quiet ? 6 : 3;
  params.eps_rel = options.epsr;
  params.stop_energy = options.stopEnergy;
  if (options.hasLa)
    params.ideal_edge_length_abs = options.la;
  else
    params.ideal_edge_length_rel = options.lr;

  if (options.maxThreads > 0)
    params.num_threads = static_cast<unsigned int>(options.maxThreads);

  Eigen::MatrixXd V;
  Eigen::MatrixXi T;
  const int ret = floatTetWild::tetrahedralization(sfMesh, params, V, T);
  cleanupTetwildScratchFiles(params.output_path);
  if (ret != 0)
    throw std::runtime_error("fTetWild tetrahedralization failed for " + options.common.inputMesh);

  return convertTetwildOutput(V, T);
}

}  // namespace tet_mesher
