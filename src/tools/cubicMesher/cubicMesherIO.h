#pragma once

#include "cubicMesh.h"

#include <string>

namespace cubic_mesher
{

void saveCubicMesh(const pgo::VolumetricMeshes::CubicMesh &cubicMesh, const std::string &outputMesh);

void writeSurfaceMesh(const pgo::VolumetricMeshes::CubicMesh &cubicMesh, const std::string &outputSurface);

}  // namespace cubic_mesher
