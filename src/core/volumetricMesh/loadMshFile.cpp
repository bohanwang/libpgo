#include "tetMesh.h"

#include <gmsh.h>

#include <iostream>
#include <vector>

namespace pgo::VolumetricMeshes
{
TetMesh loadMshFile(const char *filename)
{
  // 1. Initialize the Gmsh API
  gmsh::initialize();

  // 2. Open your .msh file
  gmsh::open(filename);

  // 3. Retrieve mesh nodes
  std::vector<std::size_t> nodeTags;
  std::vector<double> nodeCoords, nodeParams;
  gmsh::model::mesh::getNodes(nodeTags, nodeCoords, nodeParams);

  std::cout << "Loaded " << nodeTags.size() << " nodes\n";

  std::vector<pgo::Vec3d> vertices(nodeTags.size());
  std::vector<pgo::Vec4i> elements;

  // Each node has (x, y, z)
  for (size_t i = 0; i < nodeTags.size(); ++i) {
    double x = nodeCoords[3 * i + 0];
    double y = nodeCoords[3 * i + 1];
    double z = nodeCoords[3 * i + 2];

    vertices[i] = pgo::Vec3d(x, y, z);
  }

  // 4. Retrieve element connectivity
  std::vector<std::size_t> elemTags, elemNodes;
  gmsh::model::mesh::getElementsByType(4, elemTags, elemNodes);

  std::cout << "Read " << elemTags.size() << " tetrahedra\n";
  for (size_t i = 0; i < elemTags.size(); ++i) {
    pgo::Vec4i tet;
    for (int j = 0; j < 4; ++j) {
      tet[j] = static_cast<int>(elemNodes[4 * i + j] - 1);  // convert to 0-based index
    }
    elements.push_back(tet);
  }

  // 5. Finalize API
  gmsh::finalize();

  return pgo::VolumetricMeshes::TetMesh(vertices, elements);
}
}  // namespace pgo::VolumetricMeshes
