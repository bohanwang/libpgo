#include "tetMesherBackend.h"

#include "tetgenInterface.h"
#include "triMeshGeo.h"

#include <iostream>
#include <stdexcept>
#include <vector>

namespace tet_mesher
{

std::unique_ptr<pgo::VolumetricMeshes::TetMesh> generateTetgenMesh(const TetgenOptions &options)
{
  pgo::Mesh::TriMeshGeo inputMesh;
  if (inputMesh.load(options.common.inputMesh) != true)
    throw std::runtime_error("Failed to load input surface mesh " + options.common.inputMesh);

  if (options.common.quiet == false)
    std::cout << "TetGen command: " << options.command << std::endl;

  pgo::EigenSupport::MXd V;
  pgo::EigenSupport::MXi T;
  pgo::TetgenInterface::computeTetMesh(inputMesh, options.command, V, T);

  std::vector<pgo::Vec3d> positions;
  positions.reserve(V.rows());
  for (int i = 0; i < V.rows(); ++i)
    positions.emplace_back(V(i, 0), V(i, 1), V(i, 2));

  std::vector<pgo::Vec4i> tets;
  tets.reserve(T.rows());
  for (int i = 0; i < T.rows(); ++i)
    tets.emplace_back(T(i, 0), T(i, 1), T(i, 2), T(i, 3));

  auto tetMesh = std::make_unique<pgo::VolumetricMeshes::TetMesh>(positions, tets);
  tetMesh->orient();
  return tetMesh;
}

}  // namespace tet_mesher
