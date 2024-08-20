#include "geogramInterface.h"

#include "initPredicates.h"
#include "pgoLogging.h"

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_preprocessing.h>
#include <geogram/mesh/mesh_reorder.h>
#include <geogram/mesh/mesh_frame_field.h>
#include <geogram/mesh/mesh_tetrahedralize.h>
#include <geogram/mesh/mesh_decimate.h>
#include <geogram/mesh/mesh_remesh.h>
#include <geogram/basic/command_line_args.h>

#include <mutex>

namespace pgo
{
namespace GeogramInterface
{
std::once_flag geo_init;
}
}  // namespace pgo

void pgo::GeogramInterface::initGEO()
{
  std::call_once(geo_init, []() {
    GEO::initialize();
    GEO::CmdLine::import_arg_group("standard");
    GEO::CmdLine::import_arg_group("remesh");
    GEO::CmdLine::import_arg_group("algo");
  });
}

pgo::Mesh::TriMeshGeo pgo::GeogramInterface::remesh(const std::string &meshFilename, int targetNumPoints, double sizeFactor, double anisotropy)
{
  GEO::Mesh inputMesh;
  GEO::Mesh remeshedMesh;
  // GEO::mesh_load(meshFilename, inputMesh);

  if (!GEO::mesh_load(meshFilename, inputMesh)) {
    SPDLOG_LOGGER_ERROR(pgo::Logging::lgr(), "Cannot open file {}", meshFilename);
    return Mesh::TriMeshGeo();
  }

  if (sizeFactor > 0) {
    GEO::compute_sizing_field(inputMesh, sizeFactor);
  }

  if (anisotropy > 0) {
    GEO::compute_normals(inputMesh);
    GEO::simple_Laplacian_smooth(inputMesh, 3, true);  // true: smooth normals
    GEO::set_anisotropy(inputMesh, anisotropy * 0.02);
  }

  GEO::remesh_smooth(inputMesh, remeshedMesh, targetNumPoints);

  Mesh::TriMeshGeo meshOut;
  for (GEO::index_t vi = 0; vi < remeshedMesh.vertices.nb(); ++vi) {
    const double *p = remeshedMesh.vertices.point_ptr(vi);
    Vec3d pos(p);
    meshOut.addPos(pos);
  }

  for (GEO::index_t fi = 0; fi < remeshedMesh.facets.nb(); ++fi) {
    GEO::index_t nvtx = remeshedMesh.facets.nb_vertices(fi);
    PGO_ALOG(nvtx == 3);

    Vec3i tri(remeshedMesh.facets.vertex(fi, 0),
      remeshedMesh.facets.vertex(fi, 1),
      remeshedMesh.facets.vertex(fi, 2));

    meshOut.addTri(tri);
  }

  return meshOut;
}

pgo::Mesh::TriMeshGeo pgo::GeogramInterface::remesh(const Mesh::TriMeshGeo &inMesh, int targetNumPoints, double sizeFactor, double anisotropy)
{
  GEO::Mesh inputMesh;
  GEO::Mesh remeshedMesh;

  std::vector<double> positions(inMesh.numVertices() * 3);
  for (int vi = 0; vi < inMesh.numVertices(); vi++) {
    const auto &p = inMesh.pos(vi);
    positions[vi * 3] = p[0];
    positions[vi * 3 + 1] = p[1];
    positions[vi * 3 + 2] = p[2];
  }

  inputMesh.vertices.assign_points(positions.data(), (GEO::index_t)3, (GEO::index_t)inMesh.numVertices());

  GEO::vector<GEO::index_t> triangles(inMesh.numTriangles() * 3);
  for (int ti = 0; ti < inMesh.numTriangles(); ti++) {
    const auto &t = inMesh.tri(ti);
    triangles[ti * 3] = t[0];
    triangles[ti * 3 + 1] = t[1];
    triangles[ti * 3 + 2] = t[2];
  }

  inputMesh.facets.assign_triangle_mesh(triangles, false);

  if (sizeFactor > 0) {
    GEO::compute_sizing_field(inputMesh, sizeFactor);
  }

  if (anisotropy > 0) {
    GEO::compute_normals(inputMesh);
    GEO::simple_Laplacian_smooth(inputMesh, 3, true);  // true: smooth normals
    GEO::set_anisotropy(inputMesh, anisotropy * 0.02);
  }

  GEO::remesh_smooth(inputMesh, remeshedMesh, targetNumPoints);

  Mesh::TriMeshGeo meshOut;
  for (GEO::index_t vi = 0; vi < remeshedMesh.vertices.nb(); ++vi) {
    const double *p = remeshedMesh.vertices.point_ptr(vi);
    Vec3d pos(p);
    meshOut.addPos(pos);
  }

  for (GEO::index_t fi = 0; fi < remeshedMesh.facets.nb(); ++fi) {
    GEO::index_t nvtx = remeshedMesh.facets.nb_vertices(fi);
    PGO_ALOG(nvtx == 3);

    Vec3i tri(remeshedMesh.facets.vertex(fi, 0),
      remeshedMesh.facets.vertex(fi, 1),
      remeshedMesh.facets.vertex(fi, 2));

    meshOut.addTri(tri);
  }

  return meshOut;
}
