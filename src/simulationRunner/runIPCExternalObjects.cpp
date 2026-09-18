#include "runIPCExternalObjects.h"

#include "configFileJSON.h"
#include "geometryQuery.h"
#include "triMeshGeo.h"

#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

namespace pgo::RunIPCSim
{
namespace
{
namespace ES = pgo::EigenSupport;

[[noreturn]] void throwConfigError(const std::string &message)
{
  throw std::invalid_argument(message);
}

std::string externalObjectLabel(std::size_t objectIndex, const std::filesystem::path &path)
{
  return "ipc-external-objects[" + std::to_string(objectIndex) + "] (`" + path.string() + "`)";
}

ES::V3d parseOptionalFiniteVector(
  const nlohmann::json &object,
  const char *field,
  const ES::V3d &defaultValue,
  const std::string &label)
{
  if (!object.contains(field))
    return defaultValue;
  const auto &value = object[field];
  if (!value.is_array() || value.size() != 3)
    throwConfigError(label + " requires `" + field + "` with exactly three numeric values.");

  ES::V3d result;
  for (std::size_t component = 0; component < 3; ++component) {
    if (!value[component].is_number())
      throwConfigError(label + " requires numeric `" + field + "` components.");
    result[static_cast<Eigen::Index>(component)] = value[component].get<double>();
  }
  if (!result.allFinite())
    throwConfigError(label + " requires finite `" + field + "` components.");
  return result;
}

void validateExternalMesh(
  const pgo::Mesh::TriMeshGeo &mesh,
  std::size_t objectIndex,
  const std::filesystem::path &path)
{
  const std::string label = externalObjectLabel(objectIndex, path);
  if (mesh.numVertices() <= 0)
    throwConfigError(label + " must contain at least one vertex.");
  if (mesh.numTriangles() <= 0)
    throwConfigError(label + " must contain at least one triangle.");

  for (int vi = 0; vi < mesh.numVertices(); ++vi) {
    if (!mesh.pos(vi).allFinite())
      throwConfigError(label + " contains a non-finite vertex.");
  }
  for (int fi = 0; fi < mesh.numTriangles(); ++fi) {
    const auto &tri = mesh.tri(fi);
    for (int j = 0; j < 3; ++j) {
      if (tri[j] < 0 || tri[j] >= mesh.numVertices())
        throwConfigError(label + " contains an out-of-range triangle index.");
    }
    if (tri[0] == tri[1] || tri[1] == tri[2] || tri[2] == tri[0])
      throwConfigError(label + " contains a repeated-index triangle.");

    const ES::V3d e01 = mesh.pos(tri[1]) - mesh.pos(tri[0]);
    const ES::V3d e02 = mesh.pos(tri[2]) - mesh.pos(tri[0]);
    const double twiceArea = e01.cross(e02).norm();
    if (!e01.allFinite() || !e02.allFinite() || !std::isfinite(twiceArea))
      throwConfigError(label + " contains a triangle with a non-finite geometric measure.");
    if (twiceArea <= 0.0)
      throwConfigError(label + " contains a zero-area triangle.");
  }
}

}  // namespace

std::vector<IPCExternalObjectSpec> parseIPCExternalObjectSpecs(
  const pgo::ConfigFileJSON &config)
{
  if (config.exist("external-objects")) {
    throwConfigError(
      "IPC external objects use `ipc-external-objects`; `external-objects` is reserved for sampled contact.");
  }

  std::vector<IPCExternalObjectSpec> specs;
  if (!config.exist("ipc-external-objects"))
    return specs;

  const auto &objects = config.handle()["ipc-external-objects"];
  if (!objects.is_array())
    throwConfigError("`ipc-external-objects` must be an array.");

  specs.reserve(objects.size());
  for (std::size_t objectIndex = 0; objectIndex < objects.size(); ++objectIndex) {
    const auto &object = objects[objectIndex];
    if (!object.is_object())
      throwConfigError("ipc-external-objects[" + std::to_string(objectIndex) + "] must be an object.");
    if (!object.contains("filename") || !object["filename"].is_string()) {
      throwConfigError(
        "ipc-external-objects[" + std::to_string(objectIndex) + "] requires a string `filename`.");
    }

    IPCExternalObjectSpec spec;
    spec.path = std::filesystem::path(config.resolvePath(object["filename"].get<std::string>()));
    const std::string label = externalObjectLabel(objectIndex, spec.path);

    if (object.contains("scale")) {
      if (!object["scale"].is_number())
        throwConfigError(label + " requires a numeric `scale`.");
      spec.scale = object["scale"].get<double>();
    }
    if (!std::isfinite(spec.scale) || spec.scale <= 0.0)
      throwConfigError(label + " requires `scale` to be finite and strictly positive.");

    if (object.contains("init-disp"))
      throwConfigError(label + " uses `initial-translation`; per-object `init-disp` is not supported.");
    spec.initialTranslation = parseOptionalFiniteVector(
      object, "initial-translation", ES::V3d::Zero(), label);
    spec.movement = parseOptionalFiniteVector(object, "movement", ES::V3d::Zero(), label);
    if (spec.movement.cwiseAbs().maxCoeff() != 0.0)
      throwConfigError(label + " currently supports only `movement == [0, 0, 0]`.");

    specs.push_back(std::move(spec));
  }
  return specs;
}

std::vector<Contact::CIPC::IPCExternalSurface> loadIPCExternalSurfaces(
  const std::vector<IPCExternalObjectSpec> &specs)
{
  std::vector<Contact::CIPC::IPCExternalSurface> externalSurfaces;
  externalSurfaces.reserve(specs.size());
  for (std::size_t objectIndex = 0; objectIndex < specs.size(); ++objectIndex) {
    const IPCExternalObjectSpec &spec = specs[objectIndex];
    const std::string label = externalObjectLabel(objectIndex, spec.path);

    pgo::Mesh::TriMeshGeo mesh;
    if (!mesh.load(spec.path.string()))
      throwConfigError("Failed to load " + label + ".");
    validateExternalMesh(mesh, objectIndex, spec.path);
    mesh = pgo::Mesh::removeIsolatedVertices(mesh.ref());
    for (int vi = 0; vi < mesh.numVertices(); ++vi)
      mesh.pos(vi) = spec.scale * mesh.pos(vi) + spec.initialTranslation;
    validateExternalMesh(mesh, objectIndex, spec.path);

    Contact::CIPC::IPCExternalSurface surface;
    pgo::Mesh::triMeshGeoToMatrices(mesh, surface.vertices, surface.triangles);
    externalSurfaces.push_back(std::move(surface));
  }
  return externalSurfaces;
}

Contact::CIPC::IPCCollisionSurfaceData prepareIPCCollisionSurface(
  const ES::MXd &deformableVertices,
  const ES::MXi &deformableTriangles,
  const ES::SpMatD &deformableDisplacementMap,
  const pgo::ConfigFileJSON &config)
{
  const std::vector<IPCExternalObjectSpec> specs = parseIPCExternalObjectSpecs(config);
  const std::vector<Contact::CIPC::IPCExternalSurface> surfaces = loadIPCExternalSurfaces(specs);
  return Contact::CIPC::assembleIPCCollisionSurface(
    deformableVertices, deformableTriangles, deformableDisplacementMap, surfaces);
}

}  // namespace pgo::RunIPCSim
