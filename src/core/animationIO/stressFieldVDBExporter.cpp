#include "stressFieldVDBExporter.h"

#include "EigenSupport.h"
#include "initPredicates.h"
#include "pgoLogging.h"
#include "predicates.h"

#include <fmt/format.h>
#include <nlohmann/json.hpp>

#include <openvdb/io/File.h>
#include <openvdb/openvdb.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>

using namespace pgo;
using namespace pgo::AnimationIO;
namespace ES = EigenSupport;

namespace
{
ES::V3d tetVertexPosition(const VolumetricMeshes::TetMesh &mesh, const ES::VXd &u, int element, int corner)
{
  const int vi = mesh.getVertexIndex(element, corner);
  const ES::V3d &rest = mesh.getVertex(vi);
  return rest + u.segment<3>(vi * 3);
}

double averageRestEdgeLength(const VolumetricMeshes::TetMesh &mesh)
{
  static constexpr int edges[6][2] = { { 0, 1 }, { 0, 2 }, { 0, 3 }, { 1, 2 }, { 1, 3 }, { 2, 3 } };
  double sum = 0.0;
  std::int64_t count = 0;
  const int numElements = mesh.getNumElements();
  for (int el = 0; el < numElements; ++el) {
    for (const auto &e : edges) {
      const ES::V3d &p0 = mesh.getVertex(el, e[0]);
      const ES::V3d &p1 = mesh.getVertex(el, e[1]);
      sum += (p1 - p0).norm();
      ++count;
    }
  }
  return count > 0 ? sum / static_cast<double>(count) : 1.0;
}
}  // namespace

int StressFieldVDBExporter::loadTetMesh(const char *vegFilename)
{
  if (!vegFilename || !*vegFilename) {
    SPDLOG_LOGGER_ERROR(Logging::lgr(), "StressFieldVDBExporter: empty .veg path.");
    return 1;
  }

  if (!std::filesystem::exists(vegFilename)) {
    SPDLOG_LOGGER_ERROR(Logging::lgr(), "StressFieldVDBExporter: .veg file not found: {}", vegFilename);
    return 1;
  }

  tetMesh = std::make_shared<VolumetricMeshes::TetMesh>(vegFilename);
  SPDLOG_LOGGER_INFO(Logging::lgr(),
    "StressFieldVDBExporter: loaded {} ({} vertices, {} tets)",
    vegFilename, tetMesh->getNumVertices(), tetMesh->getNumElements());
  return 0;
}

int StressFieldVDBExporter::loadDeformationSequence(const char *folderPath, const char *pattern,
  int frameStart_, int frameEnd_)
{
  if (!tetMesh) {
    SPDLOG_LOGGER_ERROR(Logging::lgr(), "StressFieldVDBExporter: tet mesh must be loaded before deformation sequence.");
    return 1;
  }
  if (frameEnd_ <= frameStart_) {
    SPDLOG_LOGGER_ERROR(Logging::lgr(), "StressFieldVDBExporter: invalid frame range [{}, {}).", frameStart_, frameEnd_);
    return 1;
  }

  const int numFrames = frameEnd_ - frameStart_;
  const int n3 = tetMesh->getNumVertices() * 3;

  displacements.assign(numFrames, ES::VXd::Zero(n3));
  this->frameStart = frameStart_;

  for (int i = 0; i < numFrames; ++i) {
    const int frame = frameStart_ + i;
    const std::filesystem::path framePath =
      std::filesystem::path(folderPath) / fmt::format(fmt::runtime(pattern), frame);

    if (!std::filesystem::exists(framePath)) {
      SPDLOG_LOGGER_WARN(Logging::lgr(),
        "StressFieldVDBExporter: missing deformation frame {}: {}", frame, framePath.string());
      continue;
    }

    ES::MXd uMat;
    if (ES::readMatrix(framePath.string().c_str(), uMat) != 0) {
      SPDLOG_LOGGER_ERROR(Logging::lgr(),
        "StressFieldVDBExporter: failed to read {}", framePath.string());
      return 1;
    }

    if (static_cast<int>(uMat.rows()) != n3) {
      SPDLOG_LOGGER_ERROR(Logging::lgr(),
        "StressFieldVDBExporter: {} has {} rows, expected {}", framePath.string(), uMat.rows(), n3);
      return 1;
    }

    displacements[i] = uMat.col(0);
  }

  SPDLOG_LOGGER_INFO(Logging::lgr(),
    "StressFieldVDBExporter: loaded {} deformation frames from {}", numFrames, folderPath);
  return 0;
}

int StressFieldVDBExporter::loadVonMisesSequence(const char *folderPath, const char *pattern,
  int frameStart_, int frameEnd_)
{
  if (!tetMesh) {
    SPDLOG_LOGGER_ERROR(Logging::lgr(), "StressFieldVDBExporter: tet mesh must be loaded before stress sequence.");
    return 1;
  }
  if (frameEnd_ <= frameStart_) {
    SPDLOG_LOGGER_ERROR(Logging::lgr(), "StressFieldVDBExporter: invalid frame range [{}, {}).", frameStart_, frameEnd_);
    return 1;
  }

  const int numFrames = frameEnd_ - frameStart_;
  const int numElements = tetMesh->getNumElements();
  stresses.assign(numFrames, ES::VXd::Zero(numElements));

  for (int i = 0; i < numFrames; ++i) {
    const int frame = frameStart_ + i;
    const std::filesystem::path framePath =
      std::filesystem::path(folderPath) / fmt::format(fmt::runtime(pattern), frame);

    if (!std::filesystem::exists(framePath)) {
      SPDLOG_LOGGER_WARN(Logging::lgr(),
        "StressFieldVDBExporter: missing stress frame {}: {}", frame, framePath.string());
      continue;
    }

    std::ifstream in(framePath);
    if (!in.is_open()) {
      SPDLOG_LOGGER_ERROR(Logging::lgr(), "StressFieldVDBExporter: cannot open {}", framePath.string());
      return 1;
    }

    nlohmann::json doc;
    try {
      in >> doc;
    }
    catch (const std::exception &e) {
      SPDLOG_LOGGER_ERROR(Logging::lgr(),
        "StressFieldVDBExporter: failed to parse {}: {}", framePath.string(), e.what());
      return 1;
    }

    if (!doc.contains("values")) {
      SPDLOG_LOGGER_ERROR(Logging::lgr(),
        "StressFieldVDBExporter: {} has no \"values\" array", framePath.string());
      return 1;
    }

    const auto values = doc.at("values").get<std::vector<double>>();
    if (static_cast<int>(values.size()) != numElements) {
      SPDLOG_LOGGER_ERROR(Logging::lgr(),
        "StressFieldVDBExporter: {} has {} values, expected {}", framePath.string(), values.size(), numElements);
      return 1;
    }

    for (int e = 0; e < numElements; ++e)
      stresses[i][e] = values[e];
  }

  SPDLOG_LOGGER_INFO(Logging::lgr(),
    "StressFieldVDBExporter: loaded {} stress frames from {}", numFrames, folderPath);
  return 0;
}

int StressFieldVDBExporter::exportAnimationVDB(const char *outputDirectory, const char *prefix,
  double voxelSize) const
{
  if (!tetMesh) {
    SPDLOG_LOGGER_ERROR(Logging::lgr(), "StressFieldVDBExporter: no tet mesh.");
    return 1;
  }
  if (displacements.size() != stresses.size() || displacements.empty()) {
    SPDLOG_LOGGER_ERROR(Logging::lgr(),
      "StressFieldVDBExporter: displacement frames ({}) must match stress frames ({}) and be non-empty.",
      displacements.size(), stresses.size());
    return 1;
  }

  std::error_code ec;
  std::filesystem::create_directories(outputDirectory, ec);
  if (ec) {
    SPDLOG_LOGGER_ERROR(Logging::lgr(),
      "StressFieldVDBExporter: cannot create {}: {}", outputDirectory, ec.message());
    return 1;
  }

  const double resolvedVoxelSize = voxelSize > 0.0
    ? voxelSize
    : 0.5 * averageRestEdgeLength(*tetMesh);
  PGO_ALOG(resolvedVoxelSize > 0.0);

  Mesh::initPredicates();
  openvdb::initialize();

  const int numElements = tetMesh->getNumElements();
  const int numFrames = static_cast<int>(displacements.size());

  for (int fi = 0; fi < numFrames; ++fi) {
    const int frame = frameStart + fi;
    const ES::VXd &u = displacements[fi];
    const ES::VXd &stressPerElement = stresses[fi];

    auto grid = openvdb::FloatGrid::create(0.0f);
    grid->setTransform(openvdb::math::Transform::createLinearTransform(resolvedVoxelSize));
    grid->setName("von_mises");
    grid->setGridClass(openvdb::GRID_FOG_VOLUME);

    auto accessor = grid->getAccessor();

    for (int el = 0; el < numElements; ++el) {
      const ES::V3d a = tetVertexPosition(*tetMesh, u, el, 0);
      const ES::V3d b = tetVertexPosition(*tetMesh, u, el, 1);
      const ES::V3d c = tetVertexPosition(*tetMesh, u, el, 2);
      const ES::V3d d = tetVertexPosition(*tetMesh, u, el, 3);

      const ES::V3d bbMin = a.cwiseMin(b).cwiseMin(c).cwiseMin(d);
      const ES::V3d bbMax = a.cwiseMax(b).cwiseMax(c).cwiseMax(d);

      const openvdb::Vec3d ijkMinD = grid->transform().worldToIndex(openvdb::Vec3d(bbMin.x(), bbMin.y(), bbMin.z()));
      const openvdb::Vec3d ijkMaxD = grid->transform().worldToIndex(openvdb::Vec3d(bbMax.x(), bbMax.y(), bbMax.z()));

      const openvdb::Coord ijkMin(
        static_cast<int>(std::floor(ijkMinD.x())),
        static_cast<int>(std::floor(ijkMinD.y())),
        static_cast<int>(std::floor(ijkMinD.z())));
      const openvdb::Coord ijkMax(
        static_cast<int>(std::ceil(ijkMaxD.x())),
        static_cast<int>(std::ceil(ijkMaxD.y())),
        static_cast<int>(std::ceil(ijkMaxD.z())));

      const float stress = static_cast<float>(stressPerElement[el]);

      for (int k = ijkMin.z(); k <= ijkMax.z(); ++k) {
        for (int j = ijkMin.y(); j <= ijkMax.y(); ++j) {
          for (int i = ijkMin.x(); i <= ijkMax.x(); ++i) {
            const openvdb::Coord ijk(i, j, k);
            const openvdb::Vec3d worldP = grid->transform().indexToWorld(ijk);
            const double point[3] = { worldP.x(), worldP.y(), worldP.z() };
            if (!Mesh::pointInTet(point, a.data(), b.data(), c.data(), d.data()))
              continue;

            const float existing = accessor.getValue(ijk);
            // Voxels along shared faces may be hit by multiple tets; keep the
            // largest stress so face-shared voxels look continuous.
            if (stress > existing)
              accessor.setValue(ijk, stress);
          }
        }
      }
    }

    const std::filesystem::path outPath =
      std::filesystem::path(outputDirectory) / fmt::format("{}{:04d}.vdb", prefix, frame);

    openvdb::GridPtrVec grids;
    grids.push_back(grid);
    openvdb::io::File file(outPath.string());
    try {
      file.write(grids);
      file.close();
    }
    catch (const std::exception &e) {
      SPDLOG_LOGGER_ERROR(Logging::lgr(),
        "StressFieldVDBExporter: failed to write {}: {}", outPath.string(), e.what());
      return 1;
    }

    SPDLOG_LOGGER_INFO(Logging::lgr(),
      "StressFieldVDBExporter: wrote frame {} -> {}", frame, outPath.string());
  }

  return 0;
}
