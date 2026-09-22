#include "animationLoader.h"
#include "abcWriter.h"

#include "EigenSupport.h"
#include "configFileJSON.h"
#include "pgoLogging.h"
#include "generateSurfaceMesh.h"
#include "boundingVolumeTree.h"

#include <filesystem>
#include <iostream>

using namespace pgo;
using namespace pgo::AnimationIO;
namespace ES = EigenSupport;

int AnimationLoader::load(const char *filename)
{
  ConfigFileJSON jconfig;
  if (jconfig.open(filename) == false) {
    return 1;
  }

  int saveCacheGlobal = 1;
  if (jconfig.exist("save-cache"))
    saveCacheGlobal = jconfig.getInt("save-cache");

  auto jconfig_root = jconfig.handle();
  try {
    for (auto &jmesh : jconfig_root["meshes"]) {
      AnimationSequence aseq;
      aseq.name = jmesh.at("name").get<std::string>();
      aseq.drivingMeshFilename = jconfig.resolvePath(jmesh.at("driving-mesh").get<std::string>());
      aseq.displayMeshFilename = jconfig.resolvePath(jmesh.value("display-mesh", ""));
      aseq.sequenceName = jconfig.resolvePath(jmesh.at("sequence").get<std::string>());
      aseq.sequenceType = jmesh.at("sequence-type").get<std::string>();
      aseq.scaleString = jmesh.value("scale", "1,1,1");
      aseq.sequenceRange = jmesh.at("sequence-range").get<std::vector<int>>();
      seqs.push_back(aseq);
    }
  }
  catch (std::exception &e) {
    SPDLOG_LOGGER_ERROR(Logging::lgr(), "Failed to parse animation config file: {}", e.what());
    return 1;
  }

  for (auto &aseq : seqs) {
    double embeddedScale = 1, embeddingScale = 1, displacementScale = 1;
    if (aseq.scaleString.length()) {
      int ret = std::sscanf(aseq.scaleString.c_str(), "%lf,%lf,%lf", &embeddedScale, &embeddingScale, &displacementScale);
      PGO_ALOG(ret == 3);
    }

    if (aseq.drivingMeshFilename.find(".obj") != std::string::npos) {
      aseq.triMesh = std::make_shared<Mesh::TriMeshGeo>();
      if (aseq.triMesh->load(aseq.drivingMeshFilename) != true) {
        SPDLOG_LOGGER_ERROR(Logging::lgr(), "Failed to load mesh file: {}", aseq.drivingMeshFilename);
        return 1;
      }

      for (auto &p : aseq.triMesh->positions()) {
        p *= embeddingScale;
      }
    }
    else if (aseq.drivingMeshFilename.find(".veg") != std::string::npos) {
      aseq.tetMesh = std::make_shared<VolumetricMeshes::TetMesh>(aseq.drivingMeshFilename.c_str());

      for (int vi = 0; vi < aseq.tetMesh->getNumVertices(); vi++) {
        Vec3d p = aseq.tetMesh->getVertex(vi);
        p *= embeddingScale;
        aseq.tetMesh->setVertex(vi, p);
      }
    }

    if (aseq.displayMeshFilename.length()) {
      // The display (render) mesh is a separate OBJ embedded into the driving
      // mesh below; it is scaled by the embedded scale.
      Mesh::TriMeshGeo displayMesh;
      if (displayMesh.load(aseq.displayMeshFilename) != true) {
        SPDLOG_LOGGER_ERROR(Logging::lgr(), "Failed to load display mesh file: {}", aseq.displayMeshFilename);
        return 1;
      }
      if (displayMesh.numVertices() <= 0 || displayMesh.numTriangles() <= 0) {
        SPDLOG_LOGGER_ERROR(Logging::lgr(), "Display mesh {} has no geometry.", aseq.displayMeshFilename);
        return 1;
      }

      aseq.displayMeshV.resize(displayMesh.numVertices() * 3);
      for (int vi = 0; vi < displayMesh.numVertices(); vi++)
        aseq.displayMeshV.segment<3>(vi * 3) = displayMesh.pos(vi) * embeddedScale;

      aseq.displayMeshF.resize(displayMesh.numTriangles());
      for (int fi = 0; fi < displayMesh.numTriangles(); fi++) {
        const auto &tri = displayMesh.tri(fi);
        aseq.displayMeshF[fi] = { tri[0], tri[1], tri[2] };
      }
    }
    else {
      if (aseq.drivingMeshFilename.find(".obj") != std::string::npos) {
        ES::MXd V;
        ES::MXi F;
        Mesh::triMeshGeoToMatrices(*aseq.triMesh, V, F);

        aseq.displayMeshV.resize(V.rows() * 3);
        for (int vi = 0; vi < V.rows(); vi++) {
          aseq.displayMeshV.segment<3>(vi * 3) = V.row(vi);
        }

        aseq.displayMeshF.resize(F.rows());
        for (int i = 0; i < F.rows(); i++) {
          aseq.displayMeshF[i] = { F(i, 0), F(i, 1), F(i, 2) };
        }
      }
      else if (aseq.drivingMeshFilename.find(".veg") != std::string::npos) {
        std::vector<Vec3d> v;
        VolumetricMeshes::GenerateSurfaceMesh::computeMesh(aseq.tetMesh.get(), v, aseq.displayMeshF, false);
        aseq.displayMeshV.resize(v.size() * 3);
        for (size_t i = 0; i < v.size(); i++) {
          aseq.displayMeshV.segment<3>(i * 3) = v[i];
        }
      }
    }

    const int numDrivingVertices = aseq.triMesh ? aseq.triMesh->numVertices() : (aseq.tetMesh ? aseq.tetMesh->getNumVertices() : 0);
    if (numDrivingVertices <= 0) {
      SPDLOG_LOGGER_ERROR(Logging::lgr(), "Driving mesh {} must be an .obj or .veg file with vertices.", aseq.drivingMeshFilename);
      return 1;
    }
    const int n3 = numDrivingVertices * 3;

    if (aseq.sequenceRange.size() != 2 || aseq.sequenceRange[1] <= aseq.sequenceRange[0]) {
      SPDLOG_LOGGER_ERROR(Logging::lgr(), "sequence-range for {} must be [first, last) with last > first.", aseq.name);
      return 1;
    }
    const int numFrames = aseq.sequenceRange[1] - aseq.sequenceRange[0];

    // <name>-uAll.u next to the config caches the driving displacements
    // (rows: driving DOFs, columns: frames, before displacement scaling).
    // save-cache (default 1) enables both writing and reusing it; names
    // containing "temp" never use the cache.
    std::string uCacheFilename = jconfig.resolvePath(fmt::format("{}-uAll.u", aseq.name));
    size_t foundTemp = aseq.name.find("temp");
    int saveCache = saveCacheGlobal;
    if (saveCache) {
      if (foundTemp != std::string::npos)
        saveCache = 0;
    }

    if (saveCache && std::filesystem::exists(uCacheFilename)) {
      ES::MXd cached;
      if (ES::readMatrix(uCacheFilename.c_str(), cached) != 0) {
        SPDLOG_LOGGER_ERROR(Logging::lgr(), "Failed to read displacement cache {}; delete it to rebuild it from the sequence.", uCacheFilename);
        return 1;
      }
      if (cached.rows() != n3 || cached.cols() != numFrames) {
        SPDLOG_LOGGER_ERROR(Logging::lgr(),
          "Displacement cache {} is {}x{} but the sequence needs {}x{}; delete it to rebuild it from the sequence.",
          uCacheFilename, cached.rows(), cached.cols(), n3, numFrames);
        return 1;
      }
      SPDLOG_LOGGER_INFO(Logging::lgr(), "Using cached displacements from {}", uCacheFilename);
      aseq.drivingDisplacements = cached;
    }
    else {
      aseq.drivingDisplacements.resize(n3, numFrames);

      if (aseq.sequenceType == "objmesh") {
        if (!aseq.triMesh) {
          SPDLOG_LOGGER_ERROR(Logging::lgr(), "sequence-type objmesh requires an .obj driving mesh ({}).", aseq.name);
          return 1;
        }
        ES::VXd restPositions(n3);
        for (int vi = 0; vi < aseq.triMesh->numVertices(); vi++) {
          restPositions.segment<3>(vi * 3) = aseq.triMesh->pos(vi);
        }

        // Every frame in the range must exist: a skipped frame would shift all
        // later samples in the fixed-rate Alembic time sampling.
        for (int i = aseq.sequenceRange[0]; i < aseq.sequenceRange[1]; i++) {
          const std::string frameFilename = fmt::format(fmt::runtime(aseq.sequenceName), i);
          Mesh::TriMeshGeo frameMesh;
          if (frameMesh.load(frameFilename) != true) {
            SPDLOG_LOGGER_ERROR(Logging::lgr(), "Failed to load frame mesh: {}", frameFilename);
            return 1;
          }
          if (frameMesh.numVertices() != numDrivingVertices) {
            SPDLOG_LOGGER_ERROR(Logging::lgr(), "Frame mesh {} has {} vertices but the driving mesh has {}.",
              frameFilename, frameMesh.numVertices(), numDrivingVertices);
            return 1;
          }

          ES::VXd positions(n3);
          for (int vi = 0; vi < frameMesh.numVertices(); vi++) {
            positions.segment<3>(vi * 3) = frameMesh.pos(vi);
          }

          aseq.drivingDisplacements.col(i - aseq.sequenceRange[0]) = positions - restPositions;
          std::cout << i << ' ' << std::flush;
        }
        std::cout << std::endl;
      }
      else if (aseq.sequenceType == "u") {
        // One displacement matrix file per frame (driving DOFs x 1).
        for (int i = aseq.sequenceRange[0]; i < aseq.sequenceRange[1]; i++) {
          const std::string uFilename = fmt::format(fmt::runtime(aseq.sequenceName), i);
          ES::MXd u;
          if (ES::readMatrix(uFilename.c_str(), u) != 0) {
            SPDLOG_LOGGER_ERROR(Logging::lgr(), "Failed to load displacement file: {}", uFilename);
            return 1;
          }
          if (u.rows() != n3 || u.cols() < 1) {
            SPDLOG_LOGGER_ERROR(Logging::lgr(), "Displacement file {} is {}x{} but {} rows are required.", uFilename, u.rows(), u.cols(), n3);
            return 1;
          }
          aseq.drivingDisplacements.col(i - aseq.sequenceRange[0]) = u.col(0);
        }
      }
      else if (aseq.sequenceType == "uall") {
        // One matrix file holding every frame as a column (driving DOFs x frames).
        ES::MXd uAll;
        if (ES::readMatrix(aseq.sequenceName.c_str(), uAll) != 0) {
          SPDLOG_LOGGER_ERROR(Logging::lgr(), "Failed to load displacement file: {}", aseq.sequenceName);
          return 1;
        }
        if (uAll.rows() != n3 || uAll.cols() < aseq.sequenceRange[1]) {
          SPDLOG_LOGGER_ERROR(Logging::lgr(), "Displacement file {} is {}x{} but {} rows and at least {} columns are required.",
            aseq.sequenceName, uAll.rows(), uAll.cols(), n3, aseq.sequenceRange[1]);
          return 1;
        }
        aseq.drivingDisplacements = uAll.middleCols(aseq.sequenceRange[0], numFrames);
      }
      else {
        SPDLOG_LOGGER_ERROR(Logging::lgr(), "Unsupported sequence-type `{}` for {}; expected objmesh, u or uall.", aseq.sequenceType, aseq.name);
        return 1;
      }

      if (saveCache) {
        if (ES::writeMatrix(uCacheFilename.c_str(), aseq.drivingDisplacements) != 0)
          SPDLOG_LOGGER_WARN(Logging::lgr(), "Failed to write displacement cache {}", uCacheFilename);
      }
    }

    aseq.drivingDisplacements *= displacementScale;

    if (aseq.displayMeshFilename.length()) {
      // if the source mesh is a tet mesh,
      // we embed it into the tet mesh
      if (aseq.tetMesh) {
        SPDLOG_LOGGER_INFO(Logging::lgr(), "Rendering mesh embedded into the tet mesh");

        aseq.bary = std::make_shared<InterpolationCoordinates::BarycentricCoordinates>(
          (int)aseq.displayMeshV.size() / 3, aseq.displayMeshV.data(), aseq.tetMesh.get());

        aseq.displayDisplacements.resize(aseq.displayMeshV.rows(), aseq.drivingDisplacements.cols());
        for (int i = 0; i < (int)aseq.drivingDisplacements.cols(); i++) {
          aseq.bary->deform(aseq.drivingDisplacements.data() + i * aseq.drivingDisplacements.rows(),
            aseq.displayDisplacements.data() + i * aseq.displayDisplacements.rows());
        }
      }
      // else we embed it into the triangle surface mesh mesh
      else {
        SPDLOG_LOGGER_INFO(Logging::lgr(), "Rendering mesh embedded into the triangle mesh");
        PGO_ALOG(aseq.triMesh != nullptr);

        Mesh::TriMeshBVTree bvTree;
        bvTree.buildByInertiaPartition(*aseq.triMesh);

        std::vector<Vec3i> barycentricVertexIndices;
        std::vector<Vec3d> barycentricVertexWeights;

        double maxDist2 = 0;
        for (int vi = 0; vi < (int)aseq.displayMeshV.size() / 3; vi++) {
          ES::V3d p = aseq.displayMeshV.segment<3>(vi * 3);
          auto ret = bvTree.closestTriangleQuery(*aseq.triMesh, p);
          maxDist2 = std::max(ret.dist2, maxDist2);

          barycentricVertexIndices.emplace_back(aseq.triMesh->tri(ret.triID));
          barycentricVertexWeights.emplace_back(ret.triBaryWeight);
        }
        SPDLOG_LOGGER_INFO(Logging::lgr(), "Max dist: {}", std::sqrt(maxDist2));

        aseq.displayDisplacements.resize(aseq.displayMeshV.rows(), aseq.drivingDisplacements.cols());
        for (int i = 0; i < (int)aseq.drivingDisplacements.cols(); i++) {
          for (int vi = 0; vi < (int)aseq.displayMeshV.size() / 3; vi++) {
            const Vec3d &w = barycentricVertexWeights[vi];
            const Vec3i &idx = barycentricVertexIndices[vi];

            Vec3d ulocal[3] = {
              aseq.drivingDisplacements.col(i).segment<3>(idx[0] * 3),
              aseq.drivingDisplacements.col(i).segment<3>(idx[1] * 3),
              aseq.drivingDisplacements.col(i).segment<3>(idx[2] * 3),
            };

            Vec3d embededu = ulocal[0] * w[0] + ulocal[1] * w[1] + ulocal[2] * w[2];
            aseq.displayDisplacements.col(i).segment<3>(vi * 3) = embededu;
          }
        }
      }
    }
    else {
      aseq.displayDisplacements = aseq.drivingDisplacements;
    }
  }

  return 0;
}

int AnimationLoader::saveABC(const char *prefix)
{
  if (seqs.empty()) {
    SPDLOG_LOGGER_ERROR(Logging::lgr(), "No animation sequences loaded; nothing to write.");
    return 1;
  }

  for (const auto &aseq : seqs) {
    if (aseq.displayMeshV.size() == 0 || aseq.displayMeshF.empty() || aseq.displayDisplacements.cols() == 0) {
      SPDLOG_LOGGER_ERROR(Logging::lgr(), "Sequence {} has no geometry or no frames; refusing to write an empty archive.", aseq.name);
      return 1;
    }
    if (aseq.displayDisplacements.rows() != aseq.displayMeshV.size()) {
      SPDLOG_LOGGER_ERROR(Logging::lgr(), "Sequence {} displacements ({} rows) do not match its display mesh ({} DOFs).",
        aseq.name, aseq.displayDisplacements.rows(), aseq.displayMeshV.size());
      return 1;
    }

    std::vector<std::vector<float>> disp(aseq.displayDisplacements.cols(), std::vector<float>());
    std::vector<float> pos(aseq.displayMeshV.size());

    for (int vi = 0; vi < (int)aseq.displayMeshV.size() / 3; vi++) {
      pos[vi * 3] = float(aseq.displayMeshV[vi * 3]);
      pos[vi * 3 + 1] = float(aseq.displayMeshV[vi * 3 + 1]);
      pos[vi * 3 + 2] = float(aseq.displayMeshV[vi * 3 + 2]);
    }

    for (ES::IDX i = 0; i < aseq.displayDisplacements.cols(); i++) {
      disp[i].resize(aseq.displayDisplacements.rows());
      for (ES::IDX j = 0; j < aseq.displayDisplacements.rows(); j++) {
        disp[i][j] = float(aseq.displayDisplacements(j, i));
      }
    }

    dumpABC(fmt::format("{}/{}.abc", prefix, aseq.name).c_str(), aseq.name.c_str(), pos, disp, aseq.displayMeshF);
  }

  return 0;
}
