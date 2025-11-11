#include "animationLoader.h"
#include "abcWriter.h"

#include "configFileJSON.h"
#include "pgoLogging.h"
#include "generateSurfaceMesh.h"
#include "boundingVolumeTree.h"

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
      aseq.drivingMeshFilename = jmesh.at("driving-mesh").get<std::string>();
      aseq.displayMeshFilename = jmesh.value("display-mesh", "");
      aseq.sequenceName = jmesh.at("sequence").get<std::string>();
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
      // load obj mesh
      // scale by embedded scale
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

    std::string uCacheFilename = fmt::format("{}-uAll.u", aseq.name);
    size_t foundTemp = aseq.name.find("temp");
    int saveCache = saveCacheGlobal;
    if (saveCache) {
      if (foundTemp != std::string::npos)
        saveCache = 0;
    }

    if (saveCache == 0 || std::filesystem::exists(uCacheFilename) == false) {
      if (aseq.sequenceType == "objmesh") {
        ES::VXd restPositions(aseq.triMesh->numVertices() * 3);
        for (int vi = 0; vi < aseq.triMesh->numVertices(); vi++) {
          restPositions.segment<3>(vi * 3) = aseq.triMesh->pos(vi);
        }

        std::vector<ES::VXd> disp;
        for (int i = aseq.sequenceRange[0]; i < aseq.sequenceRange[1]; i++) {
          Mesh::TriMeshGeo frameMesh;
          if (frameMesh.load(fmt::format(fmt::runtime(aseq.sequenceName), i)) != true) {
            SPDLOG_LOGGER_WARN(Logging::lgr(), "Failed to load frame mesh: {}", fmt::format(fmt::runtime(aseq.sequenceName), i));
            continue;
          }

          ES::VXd positions(frameMesh.numVertices() * 3);
          for (int vi = 0; vi < frameMesh.numVertices(); vi++) {
            positions.segment<3>(vi * 3) = frameMesh.pos(vi);
          }

          disp.emplace_back(positions - restPositions);
          std::cout << i << ' ' << std::flush;
        }
        aseq.drivingDisplacements.resize(disp[0].size(), disp.size());
        for (size_t i = 0; i < disp.size(); i++) {
          aseq.drivingDisplacements.col(i) = disp[i];
        }
      }
      else if (aseq.sequenceType == "u") {
        // for (int i = m.sequenceRange[0]; i < m.sequenceRange[1]; i++) {
        //   std::string uFilename = fmt::format(m.sequenceName, i);

        //   std::vector<double> mat;
        //   int nr, nc;
        //   int ret = ReadMatrixFromDisk(uFilename.c_str(), nr, nc, mat);
        //   if (ret != 0)
        //     continue;

        //   ALOG(static_cast<int>(m.objMesh->getNumVertices()) * 3 == nr && nc == 1);

        //   m.displacements.emplace_back(std::move(mat));

        //   LG_ << i << ' ';
        // }

        // LG_ << '\n';
      }
      else if (aseq.sequenceType == "uall") {
        // std::string uFilename = m.sequenceName;
        // std::vector<double> mat;
        // int nr, nc;
        // int ret = ReadMatrixFromDisk(uFilename.c_str(), nr, nc, mat);
        // if (ret != 0)
        //   continue;

        // ALOG(static_cast<int>(m.objMesh->getNumVertices()) * 3 == nr);

        // for (int i = m.sequenceRange[0]; i < m.sequenceRange[1] && i < nc; i++) {
        //   std::vector<double> thisMat;
        //   thisMat.assign(mat.data() + i * nr, mat.data() + (i + 1) * nr);
        //   m.displacements.emplace_back(std::move(thisMat));
        //   LG_ << i << ' ';
        // }
        // LG_ << '\n';
      }

      // std::vector<double> uAll;
      // for (const auto &uframe : m.displacements) {
      //   uAll.insert(uAll.end(), uframe.begin(), uframe.end());
      // }

      // if (foundTemp == std::string::npos) {
      //   std::filesystem::path uAllPath = std::filesystem::path(outputDir) / std::filesystem::path(fmt::format("{}-uAll.u", m.name));
      //   WriteMatrixToDisk(uAllPath.string().c_str(), (int)m.objMesh->getNumVertices() * 3, (int)m.displacements.size(), uAll.data());
      // }
    }
    else {
      // std::vector<double> uAll;
      // int nr, nc;
      // int ret = ReadMatrixFromDisk(uAllPath.string().c_str(), nr, nc, uAll);
      // ALOG(ret == 0 && nr == static_cast<int>(m.objMesh->getNumVertices()) * 3);

      // int n3 = (int)m.objMesh->getNumVertices() * 3;
      // for (int i = 0; i < nc; i++) {
      //   m.displacements.emplace_back(uAll.data() + i * n3, uAll.data() + (i + 1) * n3);
      // }
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
  for (const auto &aseq : seqs) {
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