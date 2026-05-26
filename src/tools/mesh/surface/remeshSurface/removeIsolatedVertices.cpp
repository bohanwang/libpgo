#include "triMeshGeo.h"

int main(int argc, char **argv)
{
  if (argc < 3) {
    std::cout << "Usage: removeIsolatedVertices inputMesh outputMesh" << std::endl;
    return 0;
  }

  std::string inputMesh = argv[1];
  std::string outputMesh = argv[2];

  pgo::Mesh::TriMeshGeo mesh;
  if (mesh.load(inputMesh) == false) {
    std::cout << "Failed to load mesh from " << inputMesh << std::endl;
    return 1;
  }

  pgo::Mesh::TriMeshGeo newMesh = pgo::Mesh::removeIsolatedVertices(mesh);
  newMesh.save(outputMesh);

  return 0;
}