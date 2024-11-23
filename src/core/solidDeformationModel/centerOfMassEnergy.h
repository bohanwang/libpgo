#pragma once

#include "tetMeshGeo.h"
#include "EigenDef.h"
#include <vector>

namespace pgo
{
namespace SolidDeformationModel
{
class CenterOfMassEnergy
{
public:
  CenterOfMassEnergy(const Mesh::TetMeshGeo &tetMesh, const EigenSupport::V3d &tgtCoM, const EigenSupport::M3d &projMat);

  void func_gradient(EigenSupport::ConstRefVecXd x, EigenSupport::RefVecXd grad, double &energy, EigenSupport::V3d &com);

  EigenSupport::MXi m_tet;
  EigenSupport::V3d m_tgtCoM;
  EigenSupport::M3d m_projMat;
};
}  // namespace SolidDeformationModel
}  // namespace pgo
