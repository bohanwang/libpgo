#include "abcWriter.h"

#include <Alembic/AbcGeom/All.h>
#include <Alembic/AbcCoreHDF5/All.h>
#include <Alembic/AbcCoreOgawa/All.h>

using namespace Alembic::AbcGeom;  // Contains Abc, AbcCoreAbstract

void pgo::AnimationIO::dumpABC(const char *filename, const char *name,
  const std::vector<float> &positions,
  const std::vector<std::vector<float>> &displacements,
  const std::vector<std::vector<int>> &triangles,
  const std::vector<float> *uv)
{
  OArchive archive(
    // The hard link to the implementation.
    Alembic::AbcCoreOgawa::WriteArchive(),
    // The file name.
    // Because we're an OArchive, this is creating (or clobbering)
    // the archive with this filename.
    filename);

  TimeSampling ts(1.0 / 24, 0.0);
  Alembic::Util::uint32_t tsidx = archive.addTimeSampling(ts);
  Alembic::AbcGeom::OXform xfobj(archive.getTop(), name, tsidx);

  std::vector<int> faceCount(triangles.size(), 3);

  std::vector<int> indices;
  for (const auto &tri : triangles) {
    for (int i = 0; i < (int)tri.size(); i++)
      indices.push_back(tri[i]);
  }

  char shapeName[512];
  sprintf(shapeName, "%sShape", name);
  OPolyMesh meshobj(xfobj, shapeName, tsidx);
  OPolyMeshSchema::Sample mesh_samp(
    V3fArraySample((const V3f *)positions.data(), positions.size() / 3),
    Int32ArraySample(indices.data(), indices.size()),
    Int32ArraySample(faceCount.data(), faceCount.size()));

  if (uv) {
    OV2fGeomParam::Sample uvsamp(V2fArraySample((const V2f *)uv->data(), uv->size() / 2), kFacevaryingScope);
    mesh_samp.setUVs(uvsamp);
  }

  std::vector<float> positionsCur = positions;

  XformSample xf_samp;
  for (int i = 0; i < (int)displacements.size(); i++) {
    xfobj.getSchema().set(xf_samp);

    for (int j = 0; j < (int)positionsCur.size(); j++) {
      positionsCur[j] = positions[j] + displacements[i][j];
    }

    mesh_samp.setPositions(P3fArraySample((const V3f *)positionsCur.data(), positionsCur.size() / 3));
    meshobj.getSchema().set(mesh_samp);
  }

  std::cout << "Writing: " << archive.getName() << std::endl;
}