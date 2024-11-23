#include "triMeshGeo.h"
#include "pgoLogging.h"

#include <tiny_obj_loader.h>

#include <filesystem>

bool pgo::Mesh::TriMeshGeo::load(const std::string &filename)
{
  tinyobj::ObjReaderConfig reader_config;
  reader_config.mtl_search_path = std::filesystem::path(filename).parent_path().string();
  reader_config.triangulate = true;

  tinyobj::ObjReader reader;

  if (!reader.ParseFromFile(filename, reader_config)) {
    if (!reader.Error().empty()) {
      std::cerr << "TinyObjReader: " << reader.Error() << std::endl;
    }

    return false;
  }

  if (!reader.Warning().empty()) {
    std::cout << "TinyObjReader: " << reader.Warning() << std::endl;
  }

  static_assert(std::is_same<double, tinyobj::real_t>::value);

  auto &attrib = reader.GetAttrib();
  auto &shapes = reader.GetShapes();
  auto &materials = reader.GetMaterials();

  clear();

  for (int vi = 0; vi < (int)attrib.vertices.size() / 3; vi++) {
    positions_.emplace_back(asVec3d(attrib.vertices.data() + vi * 3));
  }

  // Loop over shapes
  for (size_t s = 0; s < shapes.size(); s++) {
    // Loop over faces(polygon)
    size_t index_offset = 0;
    for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
      size_t fv = size_t(shapes[s].mesh.num_face_vertices[f]);
      PGO_ALOG(fv == 3ull);

      Vec3i tri;
      // Loop over vertices in the face.
      for (size_t v = 0; v < fv; v++) {
        // access to vertex
        tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
        tri[v] = int(idx.vertex_index);
      }
      index_offset += fv;

      triangles_.emplace_back(tri);
    }
  }

  return true;
}