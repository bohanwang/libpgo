#pragma once

#include "EigenSupport.h"
#include "embeddedSurfaceIPCPotentialEnergy.h"
#include "ipc/assembly/ipcCollisionSurfaceAssembler.h"

#include <filesystem>
#include <vector>

namespace pgo
{
class ConfigFileJSON;

namespace RunIPCSim
{

struct IPCExternalObjectSpec
{
  std::filesystem::path path;
  double scale = 1.0;
  EigenSupport::V3d initialTranslation = EigenSupport::V3d::Zero();
  EigenSupport::V3d movement = EigenSupport::V3d::Zero();
};

std::vector<IPCExternalObjectSpec> parseIPCExternalObjectSpecs(
  const ConfigFileJSON &config);

std::vector<Contact::CIPC::IPCExternalSurface> loadIPCExternalSurfaces(
  const std::vector<IPCExternalObjectSpec> &specs);

// A kinematic external object, handled like the sampled runner's
// `external-objects`: frozen inside every frame and placed between frames.
// The geometry after scale and initial translation is kept as the base, so the
// placement for frame k is base + movement * k / (numSteps - 1) whether or not
// the run was restarted. The object never gains DOFs.
struct IPCKinematicExternalObject
{
  EigenSupport::MXd baseVertices;
  EigenSupport::V3d movement = EigenSupport::V3d::Zero();
  int firstSurfaceVertex = 0;  // offset inside the combined collision surface

  bool moves() const { return movement.cwiseAbs().maxCoeff() != 0.0; }
  EigenSupport::MXd positionsAt(double ratio) const;
};

struct IPCPreparedCollisionSurface
{
  Contact::CIPC::IPCCollisionSurfaceData surface;
  std::vector<IPCKinematicExternalObject> kinematicObjects;
};

IPCPreparedCollisionSurface prepareIPCCollisionSurface(
  const EigenSupport::MXd &deformableVertices,
  const EigenSupport::MXi &deformableTriangles,
  const EigenSupport::SpMatD &deformableDisplacementMap,
  const ConfigFileJSON &config);

// Places every moving object for `frame` of a run with `numSimSteps` frames,
// using the sampled runner's convention: the object reaches its full movement
// at the last frame, and a one-frame run keeps it at the base placement.
// Objects with zero movement are left untouched.
void placeIPCKinematicObjects(
  Contact::CIPC::EmbeddedSurfaceIPCPotentialEnergy &collisionHandler,
  const std::vector<IPCKinematicExternalObject> &objects,
  int frame, int numSimSteps);

}  // namespace RunIPCSim
}  // namespace pgo
