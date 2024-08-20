/*
author: Bohan Wang
copyright to USC,MIT
*/

#pragma once

#include "EigenSupport.h"
#include "tetMeshGeo.h"

namespace pgo
{
namespace ES = EigenSupport;

namespace SolidDeformationModel
{
// F will be column-major
namespace TetMeshMatrix
{
void generateGradientMatrix(const Mesh::TetMeshRef &tetMesh, ES::SpMatD &G);
void generateElementGradientMatrix(const Mesh::TetMeshRef &tetMesh, int tetID, ES::M9x12d &G);

void generateBasicElementLaplacianMatrix(const Mesh::TetMeshRef &tetMesh, ES::SpMatD &L, int faceNeighbor, int scale);
void generateVertexLapacianMatrix(const Mesh::TetMeshRef &tetMesh, ES::SpMatD &L);

}  // namespace TetMeshMatrix

}  // namespace SolidDeformationModel
}  // namespace pgo
