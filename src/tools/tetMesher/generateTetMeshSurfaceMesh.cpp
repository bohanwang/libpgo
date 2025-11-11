#include "tetMesh.h"
#include "generateSurfaceMesh.h"

#include <fstream>
#include <iomanip>

int main(int argc, char *argv[])
{
  pgo::VolumetricMeshes::TetMesh tetMesh(argv[1]);
  std::vector<pgo::EigenSupport::V3d> vertices;
  std::vector<std::vector<int>> faces;

  bool allEles = false;
  if (argc > 3 && argv[3][0] == 'A') {
    allEles = true;
  }

  pgo::VolumetricMeshes::GenerateSurfaceMesh::computeMesh(&tetMesh, vertices, faces, true, allEles);

  std::ofstream outFile(argv[2]);
  outFile << std::setprecision(15);
  for (size_t i = 0; i < vertices.size(); ++i) {
    outFile << "v " << vertices[i][0] << " " << vertices[i][1] << " " << vertices[i][2] << std::endl;
  }

  for (size_t i = 0; i < faces.size(); ++i) {
    outFile << "f ";
    for (size_t j = 0; j < faces[i].size(); ++j) {
      outFile << (faces[i][j] + 1) << " ";
    }
    outFile << std::endl;
  }

  outFile.close();

  return 0;
}