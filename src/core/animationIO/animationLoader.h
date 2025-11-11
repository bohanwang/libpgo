#pragma once

#include "EigenDef.h"
#include "triMeshGeo.h"
#include "tetMesh.h"
#include "barycentricCoordinates.h"

#include <string>
#include <vector>
#include <memory>

namespace pgo
{
namespace AnimationIO
{
struct AnimationSequence
{
  std::string name;
  std::string drivingMeshFilename;
  std::string displayMeshFilename;
  std::string sequenceName;
  std::string sequenceType;
  std::string scaleString;
  std::vector<int> sequenceRange;

  //
  std::shared_ptr<Mesh::TriMeshGeo> triMesh;
  std::shared_ptr<VolumetricMeshes::TetMesh> tetMesh;
  EigenSupport::MXd drivingDisplacements;

  std::shared_ptr<InterpolationCoordinates::BarycentricCoordinates> bary;
  EigenSupport::VXd displayMeshV;
  std::vector<std::vector<int>> displayMeshF;
  EigenSupport::MXd displayDisplacements;
};

class AnimationLoader
{
public:
  AnimationLoader() {}
  ~AnimationLoader() {}

  int load(const char *filename);
  int saveABC(const char *prefix);

private:
  std::vector<AnimationSequence> seqs;
};
}  // namespace AnimationIO
}  // namespace pgo