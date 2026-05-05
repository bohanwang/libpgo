#pragma once

#include <string>

namespace pgo::Tools::FBMSUnionSurface
{

struct Options
{
  std::string fbmsPath;
  std::string spherePath;
  std::string outputSurfacePath;
  std::string debugFieldMode = "union";
  std::string extractionBackend = "marching-cubes";
  std::string isoOffsetMode = "zero";
  double fbmsThickness = 0.0;
  double sphereThickness = 0.0;
  double isoOffset = 0.0;
  double paddingRatio = 0.0;
  double vdbVoxelSize = 0.0;
  double vdbHalfWidth = 3.0;
  double vdbAdaptivity = 0.0;
  int resolution = 0;
  int vdbSmoothSteps = 0;
  int minComponentTriangles = 100;
  int keepLargestComponents = -1;
  bool enableTruncating = false;
  bool filterSmallComponents = false;
  bool projectFBMSBoundaryToSphere = false;
};

Options parseOptions(int argc, char *argv[]);
void validateOptions(const Options &options);

}  // namespace pgo::Tools::FBMSUnionSurface
