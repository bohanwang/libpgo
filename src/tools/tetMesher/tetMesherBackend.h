#pragma once

#include "tetMesh.h"

#include <memory>
#include <string>

namespace tet_mesher
{

struct CommonOptions
{
  std::string inputMesh;
  std::string outputMesh;
  std::string outputSurface;
  bool printStats = false;
  bool quiet = false;
};

struct TetgenOptions
{
  CommonOptions common;
  std::string command;
};

struct TetwildOptions
{
  CommonOptions common;
  double lr = 0.05;
  double la = 0.0;
  bool hasLa = false;
  double epsr = 0.001;
  double stopEnergy = 10.0;
  int maxThreads = 0;
};

std::unique_ptr<pgo::VolumetricMeshes::TetMesh> generateTetgenMesh(const TetgenOptions &options);
std::unique_ptr<pgo::VolumetricMeshes::TetMesh> generateTetwildMesh(const TetwildOptions &options);

void saveTetMeshOutputs(const pgo::VolumetricMeshes::TetMesh &tetMesh, const CommonOptions &options);

}  // namespace tet_mesher
