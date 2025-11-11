#include "loadMshFile.h"

int main(int argc, char **argv)
{
  if (argc != 3) {
    std::cerr << "Usage: mshFileToVegFile <input.msh> <output.veg>\n";
    return 1;
  }

  const char *inputMshFile = argv[1];
  const char *outputVegFile = argv[2];

  pgo::VolumetricMeshes::TetMesh tetMesh = pgo::VolumetricMeshes::loadMshFile(inputMshFile);

  if (tetMesh.save(outputVegFile) != 0) {
    std::cerr << "Failed to write VEG file: " << outputVegFile << "\n";
    return 1;
  }

  std::cout << "Successfully converted " << inputMshFile << " to " << outputVegFile << "\n";
  return 0;
}