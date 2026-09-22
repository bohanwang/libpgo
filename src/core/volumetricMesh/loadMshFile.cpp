#include "tetMesh.h"

#include <gmsh.h>

#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
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

  // Node tags are unique but need not be contiguous or ordered, so element
  // connectivity must be mapped through the tag -> array index table rather
  // than assumed to be tag - 1.
  std::unordered_map<std::size_t, int> tagToIndex;
  tagToIndex.reserve(nodeTags.size());

  // Each node has (x, y, z)
  for (size_t i = 0; i < nodeTags.size(); ++i) {
    double x = nodeCoords[3 * i + 0];
    double y = nodeCoords[3 * i + 1];
    double z = nodeCoords[3 * i + 2];

    vertices[i] = pgo::Vec3d(x, y, z);
    tagToIndex[nodeTags[i]] = static_cast<int>(i);
  }

  // 4. Retrieve element connectivity
  std::vector<std::size_t> elemTags, elemNodes;
  gmsh::model::mesh::getElementsByType(4, elemTags, elemNodes);

  std::cout << "Read " << elemTags.size() << " tetrahedra\n";
  for (size_t i = 0; i < elemTags.size(); ++i) {
    pgo::Vec4i tet;
    for (int j = 0; j < 4; ++j) {
      const std::size_t tag = elemNodes[4 * i + j];
      auto it = tagToIndex.find(tag);
      if (it == tagToIndex.end()) {
        gmsh::finalize();
        throw std::runtime_error("loadMshFile: tetrahedron " + std::to_string(elemTags[i]) + " references unknown node tag " + std::to_string(tag));
      }
      tet[j] = it->second;
    }
    elements.push_back(tet);
  }

  // 5. Finalize API
  gmsh::finalize();

  return pgo::VolumetricMeshes::TetMesh(vertices, elements);
}
}  // namespace pgo::VolumetricMeshes
