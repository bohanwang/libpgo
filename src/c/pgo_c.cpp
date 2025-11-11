#include "pgo_c.h"

#include "basicIO.h"
#include "tetMeshGeo.h"
#include "triMeshGeo.h"
#include "generateTetMeshMatrix.h"
#include "tetMesh.h"
#include "pgoLogging.h"
#include "geometryQuery.h"
#include "boundingVolumeTree.h"
#include "initPredicates.h"
#include "EigenSupport.h"
#include "simulationMesh.h"
#include "deformationModelManager.h"
#include "tetMeshDeformationModel.h"
#include "basicIO.h"
#include "deformationModelAssembler.h"
#include "deformationModelEnergy.h"
#include "plasticModel.h"
#include "plasticModel3DDeformationGradient.h"
#include "multiVertexPullingSoftConstraints.h"
#include "implicitBackwardEulerTimeIntegrator.h"
#include "TRBDF2TimeIntegrator.h"
#include "generateMassMatrix.h"
#include "generateSurfaceMesh.h"
#include "barycentricCoordinates.h"
#include "triangleMeshExternalContactHandler.h"
#include "configFileJSON.h"
#include "pointPenetrationEnergy.h"
#include "triangleMeshSelfContactHandler.h"
#include "pointTrianglePairCouplingEnergyWithCollision.h"
#include "linearPotentialEnergy.h"
#include "NewtonRaphsonSolver.h"

#if defined(PGO_HAS_MKL)
#  include "smoothRSEnergy.h"
#endif

#if defined(PGO_HAS_CGAL)
#  include "cgalInterface.h"
#endif

#include <tbb/parallel_for.h>

#include <fmt/format.h>

#include <filesystem>

pgoTetMeshGeoStructHandle pgo_create_tetmeshgeo(int nv, double *vertices, int ntet, int *tets)
{
  pgo::Mesh::TetMeshGeo *tetmesh = new pgo::Mesh::TetMeshGeo(nv, vertices, ntet, tets);
  return reinterpret_cast<pgoTetMeshGeoStructHandle>(tetmesh);
}

pgoTetMeshGeoStructHandle pgo_create_tetmeshgeo_from_file(char *filename)
{
  try {
    pgo::VolumetricMeshes::TetMesh mesh(filename);
    pgo::Mesh::TetMeshGeo *tetmesh = new pgo::Mesh::TetMeshGeo;
    mesh.exportMeshGeometry(*tetmesh);

    std::cout << "#vtx:" << tetmesh->numVertices() << std::endl;
    std::cout << "#tets:" << tetmesh->numTets() << std::endl;

    return reinterpret_cast<pgoTetMeshGeoStructHandle>(tetmesh);
  }
  catch (...) {
    return nullptr;
  }
}

int pgo_tetmeshgeo_get_num_vertices(pgoTetMeshGeoStructHandle m)
{
  pgo::Mesh::TetMeshGeo *tetmesh = reinterpret_cast<pgo::Mesh::TetMeshGeo *>(m);

  return tetmesh->numVertices();
}

int pgo_tetmeshgeo_get_num_tets(pgoTetMeshGeoStructHandle m)
{
  pgo::Mesh::TetMeshGeo *tetmesh = reinterpret_cast<pgo::Mesh::TetMeshGeo *>(m);

  return tetmesh->numTets();
}

void pgo_tetmeshgeo_get_vertices(pgoTetMeshGeoStructHandle m, double *vertices)
{
  pgo::Mesh::TetMeshGeo *tetmesh = reinterpret_cast<pgo::Mesh::TetMeshGeo *>(m);

  for (int vi = 0; vi < tetmesh->numVertices(); vi++) {
    vertices[vi * 3] = tetmesh->pos(vi)[0];
    vertices[vi * 3 + 1] = tetmesh->pos(vi)[1];
    vertices[vi * 3 + 2] = tetmesh->pos(vi)[2];
  }
}

void pgo_tetmeshgeo_get_tets(pgoTetMeshGeoStructHandle m, int *tets)
{
  pgo::Mesh::TetMeshGeo *tetmesh = reinterpret_cast<pgo::Mesh::TetMeshGeo *>(m);

  for (int ti = 0; ti < tetmesh->numTets(); ti++) {
    tets[ti * 4] = tetmesh->tet(ti)[0];
    tets[ti * 4 + 1] = tetmesh->tet(ti)[1];
    tets[ti * 4 + 2] = tetmesh->tet(ti)[2];
    tets[ti * 4 + 3] = tetmesh->tet(ti)[3];
  }

  std::cout << tets[0] << ',' << tets[1] << ',' << tets[2] << ',' << tets[3] << std::endl;
}

void pgo_destroy_tetmeshgeo(pgoTetMeshGeoStructHandle m)
{
  pgo::Mesh::TetMeshGeo *tetmesh = reinterpret_cast<pgo::Mesh::TetMeshGeo *>(m);
  delete tetmesh;
}

pgoTetMeshStructHandle pgo_create_tetmesh(int nv, double *vertices, int ntet, int *tets, double E, double nu, double density)
{
  try {
    pgo::VolumetricMeshes::TetMesh *mesh = new pgo::VolumetricMeshes::TetMesh(nv, vertices, ntet, tets, E, nu, density);
    return reinterpret_cast<pgoTetMeshStructHandle>(mesh);
  }
  catch (...) {
    return nullptr;
  }
}

pgoTetMeshStructHandle pgo_create_tetmesh_from_file(const char *filename)
{
  try {
    pgo::VolumetricMeshes::TetMesh *mesh = new pgo::VolumetricMeshes::TetMesh(filename);
    return reinterpret_cast<pgoTetMeshStructHandle>(mesh);
  }
  catch (...) {
    return nullptr;
  }
}

void pgo_save_tetmesh_to_file(pgoTetMeshStructHandle tetMeshHandle, const char *filename)
{
  pgo::VolumetricMeshes::TetMesh *mesh = reinterpret_cast<pgo::VolumetricMeshes::TetMesh *>(tetMeshHandle);
  mesh->save(filename);
}

int pgo_tetmesh_get_num_vertices(pgoTetMeshStructHandle m)
{
  pgo::VolumetricMeshes::TetMesh *mesh = reinterpret_cast<pgo::VolumetricMeshes::TetMesh *>(m);
  pgo::Mesh::TetMeshGeo tetmeshGeo;
  mesh->exportMeshGeometry(tetmeshGeo);
  return tetmeshGeo.numVertices();
}

int pgo_tetmesh_get_num_tets(pgoTetMeshStructHandle m)
{
  pgo::VolumetricMeshes::TetMesh *mesh = reinterpret_cast<pgo::VolumetricMeshes::TetMesh *>(m);
  pgo::Mesh::TetMeshGeo tetmeshGeo;
  mesh->exportMeshGeometry(tetmeshGeo);

  return tetmeshGeo.numTets();
}

void pgo_tetmesh_get_vertices(pgoTetMeshStructHandle m, double *vertices)
{
  pgo::VolumetricMeshes::TetMesh *mesh = reinterpret_cast<pgo::VolumetricMeshes::TetMesh *>(m);
  pgo::Mesh::TetMeshGeo tetmeshGeo;
  mesh->exportMeshGeometry(tetmeshGeo);

  for (int vi = 0; vi < tetmeshGeo.numVertices(); vi++) {
    vertices[vi * 3] = tetmeshGeo.pos(vi)[0];
    vertices[vi * 3 + 1] = tetmeshGeo.pos(vi)[1];
    vertices[vi * 3 + 2] = tetmeshGeo.pos(vi)[2];
  }
}

pgoTetMeshStructHandle pgo_tetmesh_update_vertices(pgoTetMeshStructHandle m, double *vertices)
{
  pgo::VolumetricMeshes::TetMesh *tetMesh = reinterpret_cast<pgo::VolumetricMeshes::TetMesh *>(m);

  pgo::VolumetricMeshes::TetMesh *tetMeshNew = new pgo::VolumetricMeshes::TetMesh(*tetMesh);

  for (int vi = 0; vi < tetMeshNew->getNumVertices(); vi++) {
    pgo::Vec3d p(vertices + vi * 3);
    tetMeshNew->setVertex(vi, p);
  }
  return reinterpret_cast<pgoTetMeshStructHandle>(tetMeshNew);
}

void pgo_tetmesh_get_elements(pgoTetMeshStructHandle m, int *elements)
{
  pgo::VolumetricMeshes::TetMesh *mesh = reinterpret_cast<pgo::VolumetricMeshes::TetMesh *>(m);
  pgo::Mesh::TetMeshGeo tetmeshGeo;
  mesh->exportMeshGeometry(tetmeshGeo);

  for (int ti = 0; ti < tetmeshGeo.numTets(); ti++) {
    for (int i = 0; i < 4; i++) {
      elements[ti * 4 + i] = tetmeshGeo.tet(ti)[i];
    }
  }
}

pgoSmoothRSEnergyStructHandle pgo_create_smooth_rs_energy(pgoTetMeshGeoStructHandle m, double coeffR, double coeffS)
{
#if defined(PGO_HAS_MKL)
  pgo::Mesh::TetMeshGeo *tetmesh = reinterpret_cast<pgo::Mesh::TetMeshGeo *>(m);
  std::cout << tetmesh->numVertices() << '\n'
            << tetmesh->numTets() << std::endl;

  pgo::EigenSupport::SpMatD G;
  pgo::SolidDeformationModel::TetMeshMatrix::generateGradientMatrix(*tetmesh, G);

  std::cout << "G done." << std::endl;

  pgo::EigenSupport::SpMatD L;
  pgo::SolidDeformationModel::TetMeshMatrix::generateBasicElementLaplacianMatrix(*tetmesh, L, 0, 0);

  std::cout << "L done." << std::endl;

  double coeffs[2] = { coeffR, coeffS };
  pgo::PredefinedPotentialEnergies::SmoothRSEnergy *energy = new pgo::PredefinedPotentialEnergies::SmoothRSEnergy(*tetmesh, G, L, coeffs);

  std::cout << "energy done." << std::endl;

  return reinterpret_cast<pgoSmoothRSEnergyStructHandle>(energy);
#else
  std::cerr << "No MKL. Unsupported." << std::endl;
  return nullptr;
#endif
}

void pgo_destroy_smooth_rs_energy(pgoSmoothRSEnergyStructHandle energy)
{
#if defined(PGO_HAS_MKL)
  pgo::PredefinedPotentialEnergies::SmoothRSEnergy *eng = reinterpret_cast<pgo::PredefinedPotentialEnergies::SmoothRSEnergy *>(energy);
  delete eng;
#else
  std::cerr << "No MKL. Unsupported." << std::endl;
#endif
}

double pgo_smooth_rs_energy_func(pgoSmoothRSEnergyStructHandle energy, double *x)
{
#if defined(PGO_HAS_MKL)
  namespace ES = pgo::EigenSupport;
  pgo::PredefinedPotentialEnergies::SmoothRSEnergy *eng = reinterpret_cast<pgo::PredefinedPotentialEnergies::SmoothRSEnergy *>(energy);
  return eng->func(ES::Mp<const ES::VXd>(x, eng->getNumDOFs()));
#else
  std::cerr << "No MKL. Unsupported." << std::endl;
  return 0;
#endif
}

void pgo_smooth_rs_energy_grad(pgoSmoothRSEnergyStructHandle energy, double *x, double *grad)
{
#if defined(PGO_HAS_MKL)
  namespace ES = pgo::EigenSupport;
  pgo::PredefinedPotentialEnergies::SmoothRSEnergy *eng = reinterpret_cast<pgo::PredefinedPotentialEnergies::SmoothRSEnergy *>(energy);
  eng->gradient(ES::Mp<const ES::VXd>(x, eng->getNumDOFs()), ES::Mp<ES::VXd>(grad, eng->getNumDOFs()));
#else
  std::cerr << "No MKL. Unsupported." << std::endl;
#endif
}

int64_t pgo_smooth_rs_energy_hess_num_entries(pgoSmoothRSEnergyStructHandle energy)
{
#if defined(PGO_HAS_MKL)
  namespace ES = pgo::EigenSupport;
  pgo::PredefinedPotentialEnergies::SmoothRSEnergy *eng = reinterpret_cast<pgo::PredefinedPotentialEnergies::SmoothRSEnergy *>(energy);

  ES::SpMatD H;
  eng->createHessian(H);

  return (int64_t)H.nonZeros();
#else
  std::cerr << "No MKL. Unsupported." << std::endl;
  return 0;
#endif
}

void pgo_smooth_rs_energy_hess(pgoSmoothRSEnergyStructHandle energy, double *x, int *rows, int *cols, double *values)
{
#if defined(PGO_HAS_MKL)
  namespace ES = pgo::EigenSupport;
  pgo::PredefinedPotentialEnergies::SmoothRSEnergy *eng = reinterpret_cast<pgo::PredefinedPotentialEnergies::SmoothRSEnergy *>(energy);

  ES::SpMatD H;
  eng->createHessian(H);

  memset(H.valuePtr(), 0, H.nonZeros() * sizeof(double));
  eng->hessian(ES::Mp<const ES::VXd>(x, eng->getNumDOFs()), H);

  int64_t inc = 0;
  for (ES::IDX rowi = 0; rowi < H.rows(); rowi++) {
    for (ES::SpMatD::InnerIterator it(H, rowi); it; ++it) {
      rows[inc] = (int)it.row();
      cols[inc] = (int)it.col();
      values[inc] = (int)it.value();
      inc++;
    }
  }
#else
  std::cerr << "No MKL. Unsupported." << std::endl;
#endif
}

int64_t pgo_sparse_matrix_get_num_entries(pgoSparseMatrixStructHandle m)
{
  namespace ES = pgo::EigenSupport;
  ES::SpMatD *mat = reinterpret_cast<ES::SpMatD *>(m);

  return mat->nonZeros();
}

void pgo_sparse_matrix_get_row_indices(pgoSparseMatrixStructHandle m, int *rows)
{
  namespace ES = pgo::EigenSupport;

  ES::SpMatD *mat = reinterpret_cast<ES::SpMatD *>(m);

  int64_t inc = 0;
  for (ES::IDX rowi = 0; rowi < mat->rows(); rowi++) {
    for (ES::SpMatD::InnerIterator it(*mat, rowi); it; ++it) {
      rows[inc] = (int)it.row();
      inc++;
    }
  }
}

void pgo_sparse_matrix_get_col_indices(pgoSparseMatrixStructHandle m, int *cols)
{
  namespace ES = pgo::EigenSupport;
  ES::SpMatD *mat = reinterpret_cast<ES::SpMatD *>(m);

  int64_t inc = 0;
  for (ES::IDX rowi = 0; rowi < mat->rows(); rowi++) {
    for (ES::SpMatD::InnerIterator it(*mat, rowi); it; ++it) {
      cols[inc] = (int)it.col();
      inc++;
    }
  }
}

void pgo_sparse_matrix_get_values(pgoSparseMatrixStructHandle m, double *values)
{
  namespace ES = pgo::EigenSupport;
  ES::SpMatD *mat = reinterpret_cast<ES::SpMatD *>(m);

  int64_t inc = 0;
  for (ES::IDX rowi = 0; rowi < mat->rows(); rowi++) {
    for (ES::SpMatD::InnerIterator it(*mat, rowi); it; ++it) {
      values[inc] = it.value();
      inc++;
    }
  }
}

void pgo_destroy_sparse_matrix(pgoSparseMatrixStructHandle m)
{
  namespace ES = pgo::EigenSupport;
  ES::SpMatD *mat = reinterpret_cast<ES::SpMatD *>(m);

  delete mat;
}

pgoSparseMatrixStructHandle pgo_create_tet_laplacian_matrix(pgoTetMeshGeoStructHandle m, int faceNeighbor, int repeat, int scale)
{
  namespace ES = pgo::EigenSupport;
  pgo::Mesh::TetMeshGeo *tetmesh = reinterpret_cast<pgo::Mesh::TetMeshGeo *>(m);

  ES::SpMatD L;
  pgo::SolidDeformationModel::TetMeshMatrix::generateBasicElementLaplacianMatrix(*tetmesh, L, faceNeighbor, scale);

  if (repeat > 1) {
    ES::SpMatD Ln;
    ES::expandN(L, Ln, repeat);

    L = Ln;
  }

  ES::SpMatD *ret = new ES::SpMatD(L);
  return reinterpret_cast<pgoSparseMatrixStructHandle>(ret);
}

pgoSparseMatrixStructHandle pgo_create_tet_gradient_matrix(pgoTetMeshGeoStructHandle m)
{
  namespace ES = pgo::EigenSupport;
  pgo::Mesh::TetMeshGeo *tetmesh = reinterpret_cast<pgo::Mesh::TetMeshGeo *>(m);

  ES::SpMatD G;
  pgo::SolidDeformationModel::TetMeshMatrix::generateGradientMatrix(*tetmesh, G);

  return reinterpret_cast<pgoSparseMatrixStructHandle>(new ES::SpMatD(G));
}

void pgo_create_tet_gradient_per_element_matrix(pgoTetMeshGeoStructHandle m, double *smallMats)
{
  namespace ES = pgo::EigenSupport;
  pgo::Mesh::TetMeshGeo *tetmesh = reinterpret_cast<pgo::Mesh::TetMeshGeo *>(m);

  for (int i = 0; i < tetmesh->numTets(); i++) {
    ES::M9x12d G;
    pgo::SolidDeformationModel::TetMeshMatrix::generateElementGradientMatrix(*tetmesh, i, G);

    ES::Mp<ES::M9x12d>(smallMats + i * 9 * 12) = G;
  }
}

pgoSparseMatrixStructHandle pgo_create_tet_biharmonic_gradient_matrix(pgoTetMeshGeoStructHandle m, int faceNeighbor, int scale)
{
  namespace ES = pgo::EigenSupport;
  pgo::Mesh::TetMeshGeo *tetmesh = reinterpret_cast<pgo::Mesh::TetMeshGeo *>(m);

  ES::SpMatD L;
  pgo::SolidDeformationModel::TetMeshMatrix::generateBasicElementLaplacianMatrix(*tetmesh, L, faceNeighbor, 0);

  ES::SpMatD Ln;
  ES::expandN(L, Ln, 9);

  ES::SpMatD G;
  pgo::SolidDeformationModel::TetMeshMatrix::generateGradientMatrix(*tetmesh, G);

  ES::SpMatD LG;
  ES::mm(Ln, G, LG);

  ES::SpMatD *ret = new ES::SpMatD;
  ES::SpMatD &GTLTLG = *ret;
  if (scale) {
    std::vector<double> eleVols(tetmesh->numTets());
    for (int ei = 0; ei < tetmesh->numTets(); ei++) {
      eleVols[ei] = std::abs(pgo::Mesh::getTetDeterminant(tetmesh->pos(ei, 0), tetmesh->pos(ei, 1), tetmesh->pos(ei, 2), tetmesh->pos(ei, 3)));
    }

    double maxVol = *std::max_element(eleVols.begin(), eleVols.end());
    for (int ei = 0; ei < tetmesh->numTets(); ei++) {
      eleVols[ei] /= maxVol;
    }

    ES::SpMatD M;
    std::vector<ES::TripletD> entries;
    for (int i = 0; i < (int)tetmesh->numTets(); i++) {
      for (int j = 0; j < 9; j++) {
        entries.emplace_back(i * 9 + j, i * 9 + j, 1.0 / eleVols[i]);
      }
    }
    M.resize(Ln.rows(), Ln.cols());
    M.setFromTriplets(entries.begin(), entries.end());

    ES::SpMatD MLG;
    ES::mm(M, LG, MLG);

    ES::mm(LG, MLG, GTLTLG, 1);
  }
  else {
    ES::mm(LG, LG, GTLTLG, 1);
  }

  return reinterpret_cast<pgoSparseMatrixStructHandle>(ret);
}

double pgo_conjugate_mv(pgoSparseMatrixStructHandle m, double *v)
{
  namespace ES = pgo::EigenSupport;
  ES::SpMatD *M = reinterpret_cast<ES::SpMatD *>(m);
  ES::Mp<ES::VXd> vmap(v, M->rows());

  return vmap.dot(*M * vmap) * 0.5;
}

void pgo_sp_mv(pgoSparseMatrixStructHandle m, double *v, double *vout)
{
  namespace ES = pgo::EigenSupport;
  ES::SpMatD *M = reinterpret_cast<ES::SpMatD *>(m);
  ES::Mp<ES::VXd> vmap(v, M->rows());
  ES::Mp<ES::VXd> vmap_out(vout, M->rows());

  ES::mv(*M, vmap, vmap_out);
}

void pgo_init()
{
  pgo::Logging::init();
}

pgoTriMeshGeoStructHandle pgo_create_trimeshgeo(int nv, double *vertices, int ntri, int *tris)
{
  pgo::Mesh::TriMeshGeo *mesh = new pgo::Mesh::TriMeshGeo(nv, vertices, ntri, tris);
  return reinterpret_cast<pgoTriMeshGeoStructHandle>(mesh);
}

int pgo_trimeshgeo_get_num_vertices(pgoTriMeshGeoStructHandle trimesh)
{
  pgo::Mesh::TriMeshGeo *mesh = reinterpret_cast<pgo::Mesh::TriMeshGeo *>(trimesh);
  return mesh->numVertices();
}

int pgo_trimeshgeo_get_num_triangles(pgoTriMeshGeoStructHandle trimesh)
{
  pgo::Mesh::TriMeshGeo *mesh = reinterpret_cast<pgo::Mesh::TriMeshGeo *>(trimesh);
  return mesh->numTriangles();
}

void pgo_trimeshgeo_get_vertices(pgoTriMeshGeoStructHandle trimesh, double *vertices)
{
  pgo::Mesh::TriMeshGeo *mesh = reinterpret_cast<pgo::Mesh::TriMeshGeo *>(trimesh);
  memcpy(vertices, mesh->positions()[0].data(), sizeof(double) * mesh->positions().size() * 3);
}

void pgo_trimeshgeo_get_triangles(pgoTriMeshGeoStructHandle trimesh, int *tris)
{
  pgo::Mesh::TriMeshGeo *mesh = reinterpret_cast<pgo::Mesh::TriMeshGeo *>(trimesh);
  memcpy(tris, mesh->triangles()[0].data(), sizeof(int) * mesh->triangles().size() * 3);
}

void pgo_destroy_trimeshgeo(pgoTriMeshGeoStructHandle trimesh)
{
  pgo::Mesh::TriMeshGeo *mesh = reinterpret_cast<pgo::Mesh::TriMeshGeo *>(trimesh);
  delete mesh;
}

void pgo_mesh_segmentation(pgoTriMeshGeoStructHandle trimesh, int numSegs, int *classIDs)
{
#if defined(PGO_HAS_CGAL)
  pgo::Mesh::TriMeshGeo *mesh = reinterpret_cast<pgo::Mesh::TriMeshGeo *>(trimesh);
  std::vector<int> cids;
// pgo::CGALInterface::segmentMesh(*mesh, numSegs, cids);
// memcpy(classIDs, cids.data(), cids.size() * sizeof(int));
#else
  std::cout << "No CGAL. Unsupported." << std::endl;
#endif
}

pgoTriMeshGeoStructHandle pgo_mesh_isotropic_remeshing(pgoTriMeshGeoStructHandle trimesh, double targetEdgeLength, int nIter, double angleThreshold)
{
#if defined(PGO_HAS_CGAL)
  pgo::Mesh::TriMeshGeo *mesh = reinterpret_cast<pgo::Mesh::TriMeshGeo *>(trimesh);
  pgo::Mesh::TriMeshGeo meshOut = pgo::CGALInterface::isotropicRemeshing(*mesh, targetEdgeLength, nIter, angleThreshold);
  return reinterpret_cast<pgoTriMeshGeoStructHandle>(new pgo::Mesh::TriMeshGeo(meshOut));
#else
  std::cout << "No CGAL. Unsupported." << std::endl;
  return nullptr;
#endif
}

void pgo_trimesh_closest_distances(pgoTriMeshGeoStructHandle trimesh, int n, double *queryPos, double *queryDistance, int *queryTri)
{
  pgo::Mesh::initPredicates();

  pgo::Mesh::TriMeshGeo *mesh = reinterpret_cast<pgo::Mesh::TriMeshGeo *>(trimesh);
  pgo::Mesh::TriMeshBVTree bvTree;
  bvTree.buildByInertiaPartition(*mesh);

  tbb::parallel_for(0, n, [&](int i) {
    pgo::Vec3d pt(queryPos + i * 3);
    auto ret = bvTree.closestTriangleQuery(*mesh, pt);
    queryDistance[i] = ret.dist2;

    if (queryTri) {
      queryTri[i] = ret.triID;
    }
  });
}

void pgo_tetmesh_barycentric_weights(pgoTetMeshGeoStructHandle tetmesh, int n, double *queryPos, double *queryW, int *queryEle)
{
  pgo::Mesh::initPredicates();

  pgo::Mesh::TetMeshGeo *mesh = reinterpret_cast<pgo::Mesh::TetMeshGeo *>(tetmesh);
  pgo::Mesh::TetMeshBVTree bvTree;
  bvTree.buildByInertiaPartition(*mesh);

  tbb::parallel_for(0, n, [&](int i) {
    pgo::Vec3d pt(queryPos + i * 3);
    int ele = bvTree.getClosestTet(*mesh, pt);
    queryEle[i] = ele;
    pgo::Mesh::getTetBarycentricWeights(pt, mesh->pos(ele, 0), mesh->pos(ele, 1), mesh->pos(ele, 2), mesh->pos(ele, 3), queryW + i * 4);
  });
}

int pgo_run_sim_from_config(const char *configFileName)
{
  using namespace pgo;
  using namespace pgo::EigenSupport;
  namespace ES = EigenSupport;

  pgo::Mesh::initPredicates();

  ConfigFileJSON jconfig;
  if (jconfig.open(configFileName) != true) {
    return 0;
  }

  // tet mesh filename
  std::string tetMeshFilename = jconfig.getString("tet-mesh", 1);

  // surface mesh filename
  std::string surfaceMeshFilename = jconfig.getString("surface-mesh", 1);

  // external acceleration
  ES::V3d extAcc = ES::Mp<ES::V3d>(jconfig.getValue<std::array<double, 3>>("g", 1).data());

  // initial velocity
  ES::V3d initialVel = ES::Mp<ES::V3d>(jconfig.getValue<std::array<double, 3>>("init-vel", 1).data());

  ES::V3d initialDisp = ES::Mp<ES::V3d>(jconfig.getValue<std::array<double, 3>>("init-disp", 1).data());

  // scene scale
  double scale = jconfig.getDouble("scale", 1);

  // timestep
  double timestep = jconfig.getDouble("timestep", 1);

  // contact param
  double contactK = jconfig.getDouble("contact-stiffness", 1);
  int contactSamples = jconfig.getInt("contact-sample", 1);
  double fricCoeff = jconfig.getDouble("contact-friction-coeff", 1);
  double velEps = jconfig.getDouble("contact-vel-eps", 1);

  // solver param
  double solverEps = jconfig.getDouble("solver-eps", 1);
  int solverMaxIter = jconfig.getInt("solver-max-iter", 1);

  // damping
  std::array<double, 2> dampingParams = jconfig.getValue<std::array<double, 2>>("damping-params", 1);

  // material
  std::string material = jconfig.getString("elastic-material");
  pgo::SolidDeformationModel::DeformationModelElasticMaterial elasticMat;
  if (material == "stable-neo") {
    elasticMat = pgo::SolidDeformationModel::DeformationModelElasticMaterial::STABLE_NEO;
  }
  else if (material == "stvk-vol") {
    elasticMat = pgo::SolidDeformationModel::DeformationModelElasticMaterial::STVK_VOL;
  }
  else {
    SPDLOG_LOGGER_ERROR(Logging::lgr(), "Unsupported elastic material: {}", material);
    return 1;
  }

  int numSimSteps = jconfig.getInt("num-timestep", 1);
  int frameGap = jconfig.getInt("dump-interval", 1);

  // sim type
  std::string simType = jconfig.getString("sim-type");

  // output
  std::string outputFolder = jconfig.getString("output", 1);

  VolumetricMeshes::TetMesh tetMesh(tetMeshFilename.c_str());
  for (int vi = 0; vi < tetMesh.getNumVertices(); vi++) {
    tetMesh.setVertex(vi, tetMesh.getVertex(vi) * scale);
  }

  Mesh::TriMeshGeo surfaceMesh;
  if (surfaceMesh.load(surfaceMeshFilename) != true)
    return 1;

  for (int vi = 0; vi < surfaceMesh.numVertices(); vi++) {
    surfaceMesh.pos(vi) *= scale;
  }

  int surfn = surfaceMesh.numVertices();
  int surfn3 = surfn * 3;
  ES::VXd surfaceRestPositions(surfn3);
  for (int vi = 0; vi < surfn; vi++) {
    surfaceRestPositions.segment<3>(vi * 3) = surfaceMesh.pos(vi);
  }

  // initialize interpolation weights
  InterpolationCoordinates::BarycentricCoordinates bc(surfaceMesh.numVertices(), surfaceRestPositions.data(), &tetMesh);
  ES::SpMatD W = bc.generateInterpolationMatrix();

  // initialize fem
  std::shared_ptr<SolidDeformationModel::SimulationMesh> simMesh(SolidDeformationModel::loadTetMesh(&tetMesh));
  std::shared_ptr<SolidDeformationModel::DeformationModelManager> dmm = std::make_shared<SolidDeformationModel::DeformationModelManager>();

  dmm->setMesh(simMesh.get(), nullptr, nullptr);
  dmm->init(pgo::SolidDeformationModel::DeformationModelPlasticMaterial::VOLUMETRIC_DOF6, elasticMat, 1);

  std::vector<double> elementWeights(simMesh->getNumElements(), 1.0);
  std::shared_ptr<SolidDeformationModel::DeformationModelAssembler> assembler =
    std::make_shared<SolidDeformationModel::DeformationModelAssembler>(dmm, elementWeights.data());

  int n = simMesh->getNumVertices();
  int n3 = n * 3;
  int nele = simMesh->getNumElements();

  ES::VXd plasticity(nele * 6);
  ES::M3d I = ES::M3d::Identity();
  for (int ei = 0; ei < nele; ei++) {
    const SolidDeformationModel::PlasticModel3DDeformationGradient *pm =
      dynamic_cast<const SolidDeformationModel::PlasticModel3DDeformationGradient *>(dmm->getDeformationModel(ei)->getPlasticModel());
    if (!pm) {
      SPDLOG_LOGGER_ERROR(Logging::lgr(), "Plastic model is not of type PlasticModel3DDeformationGradient.");
      return 1;
    }
    pm->toParam(I.data(), plasticity.data() + ei * dmm->getNumPlasticParameters());
  }

  ES::VXd restPosition(n3);
  for (int vi = 0; vi < n; vi++) {
    double p[3];
    simMesh->getVertex(vi, p);
    restPosition.segment<3>(vi * 3) = ES::V3d(p[0], p[1], p[2]);
  }

  std::shared_ptr<SolidDeformationModel::DeformationModelEnergy> elasticEnergy =
    std::make_shared<SolidDeformationModel::DeformationModelEnergy>(assembler, &restPosition, 0);
  elasticEnergy->setPlasticParams(plasticity);

  ES::VXd zero(n3);
  zero.setZero();

  ES::SpMatD K;
  elasticEnergy->createHessian(K);
  elasticEnergy->hessian(zero, K);

  // attachments
  std::vector<std::shared_ptr<ConstraintPotentialEnergies::MultipleVertexPulling>> pullingEnergies;
  std::vector<ES::VXd> pullingTargets, pullingTargetRests;
  for (const auto &fv : jconfig.handle()["fixed-vertices"]) {
    std::string filename = fv["filename"].get<std::string>();
    std::array<double, 3> movement = fv["movement"].get<std::array<double, 3>>();
    double attachmentCoeff = fv["coeff"].get<double>();

    std::vector<int> fixedVertices;
    if (BasicIO::read1DText(filename.c_str(), std::back_inserter(fixedVertices)) != 0) {
      return 1;
    }
    std::sort(fixedVertices.begin(), fixedVertices.end());

    ES::VXd tgtVertexPositions(fixedVertices.size() * 3);
    ES::VXd tgtVertexRests(fixedVertices.size() * 3);
    for (int vi = 0; vi < (int)fixedVertices.size(); vi++) {
      tgtVertexPositions.segment<3>(vi * 3) = restPosition.segment<3>(fixedVertices[vi] * 3) + ES::Mp<ES::V3d>(movement.data());
      tgtVertexRests.segment<3>(vi * 3) = restPosition.segment<3>(fixedVertices[vi] * 3);
    }

    // initialize fixed constraints
    auto pullingEnergy = std::make_shared<ConstraintPotentialEnergies::MultipleVertexPulling>(K, restPosition.data(),
      (int)fixedVertices.size(), fixedVertices.data(), tgtVertexPositions.data(), nullptr, 1);
    pullingEnergy->setCoeff(attachmentCoeff);
    pullingEnergies.push_back(pullingEnergy);
    pullingTargets.push_back(tgtVertexPositions);
    pullingTargetRests.push_back(tgtVertexRests);
  }

  std::vector<std::string> kinematicObjectFilenames;
  std::vector<ES::V3d> kinematicObjectMovements;
  if (jconfig.exist("external-objects")) {
    auto jkinObjects = jconfig.handle()["external-objects"];
    for (const auto &jko : jkinObjects) {
      std::string koFilename = jko["filename"].get<std::string>();
      kinematicObjectFilenames.push_back(koFilename);
      kinematicObjectMovements.push_back(ES::Mp<ES::V3d>(jko["movement"].get<std::array<double, 3>>().data()));
    }
  }

  ES::SpMatD M;
  VolumetricMeshes::GenerateMassMatrix::computeMassMatrix(&tetMesh, M, true);

  // initialize gravity
  ES::VXd g(n3);
  for (int vi = 0; vi < n; vi++) {
    g.segment<3>(vi * 3) = extAcc;
  }

  ES::VXd fext(n3);
  ES::mv(M, g, fext);

  if (simType == "dynamic") {
    // initialize contact
    std::vector<Mesh::TriMeshGeo> kinematicObjects;
    for (const auto &filename : kinematicObjectFilenames) {
      kinematicObjects.emplace_back();
      if (kinematicObjects.back().load(filename) != true) {
        return 1;
      }

      for (int vi = 0; vi < kinematicObjects.back().numVertices(); vi++) {
        kinematicObjects.back().pos(vi) *= scale;
      }
    }

    std::vector<Mesh::TriMeshRef> kinematicObjectsRef;
    for (int i = 0; i < (int)kinematicObjects.size(); i++) {
      kinematicObjectsRef.emplace_back(kinematicObjects[i]);
    }

    std::shared_ptr<Contact::TriangleMeshExternalContactHandler> externalContactHandler;
    std::shared_ptr<Contact::TriangleMeshSelfContactHandler> selfCD;
    if (kinematicObjectsRef.size() && contactK > 0) {
      externalContactHandler = std::make_shared<Contact::TriangleMeshExternalContactHandler>(surfaceMesh.positions(), surfaceMesh.triangles(), n3,
        kinematicObjectsRef, contactSamples, &bc.getEmbeddingVertexIndices(), &bc.getEmbeddingWeights());
    }

    if (contactK > 0) {
      selfCD = std::make_shared<Contact::TriangleMeshSelfContactHandler>(
        surfaceMesh.positions(), surfaceMesh.triangles(), n3, contactSamples, &bc.getEmbeddingVertexIndices(), &bc.getEmbeddingWeights());
    }

    std::shared_ptr<Simulation::ImplicitBackwardEulerTimeIntegrator> intg =
      std::make_shared<Simulation::ImplicitBackwardEulerTimeIntegrator>(M, elasticEnergy,
        dampingParams[0], dampingParams[1], timestep, solverMaxIter, solverEps);

#if defined(PGO_HAS_KNITRO)
    intg->setSolverOption(Simulation::TimeIntegratorSolverOption::SO_KNITRO);
    intg->setSolverConfigFile("config.opt");
#endif

    for (auto pullingEnergy : pullingEnergies)
      intg->addImplicitForceModel(pullingEnergy, 0, 0);

    intg->setExternalForce(fext.data());

    ES::VXd curSurfacePos = surfaceRestPositions;
    ES::VXd x = restPosition, u(n3);
    ES::VXd uvel(n3), uacc(n3), usurf(surfn3);

    usurf.setZero();
    u.setZero();
    uvel.setZero();
    uacc.setZero();

    for (int i = 0; i < n; i++) {
      uvel.segment<3>(i * 3) = ES::V3d(initialVel[0], initialVel[1], initialVel[2]);
    }

    if (!std::filesystem::exists(outputFolder)) {
      std::filesystem::create_directories(outputFolder);
    }

    for (int framei = 0; framei < numSimSteps; framei++) {
      intg->clearGeneralImplicitForceModel();

      double ratio = (double)framei / (numSimSteps - 1);
      for (size_t pi = 0; pi < pullingEnergies.size(); pi++) {
        ES::VXd restTgt = pullingTargetRests[pi];
        ES::VXd curTgt = restTgt * (1 - ratio) + pullingTargets[pi] * ratio;
        pullingEnergies[pi]->setTargetPos(curTgt.data());

        std::cout << "Frame " << framei << ", attachment " << pi << " target: " << curTgt.transpose().head(3) << std::endl;
      }

      std::shared_ptr<Contact::PointPenetrationEnergy> extContactEnergy;
      Contact::PointPenetrationEnergyBuffer *extContactBuffer = nullptr;
      if (externalContactHandler) {
        externalContactHandler->execute(usurf.data());

        if (externalContactHandler->getNumCollidingSamples()) {
          extContactEnergy = externalContactHandler->buildContactEnergy();
          extContactBuffer = extContactEnergy->allocateBuffer();

          auto posFunc = [&restPosition](const EigenSupport::V3d &u, EigenSupport::V3d &p, int dofStart) {
            p = u + restPosition.segment<3>(dofStart);
          };

          auto lastPosFunc = [&restPosition, &u](const EigenSupport::V3d &x, EigenSupport::V3d &p, int dofStart) {
            p = u.segment<3>(dofStart) + restPosition.segment<3>(dofStart);
          };

          extContactEnergy->setComputePosFunction(posFunc);
          extContactEnergy->setBuffer(extContactBuffer);
          extContactEnergy->setCoeff(contactK);

          extContactEnergy->setFrictionCoeff(fricCoeff);
          extContactEnergy->setComputeLastPosFunction(lastPosFunc);
          extContactEnergy->setVelEps(velEps);
          extContactEnergy->setTimestep(timestep);

          intg->addGeneralImplicitForceModel(extContactEnergy, 0, 0);
        }
      }

      std::shared_ptr<Contact::PointTrianglePairCouplingEnergyWithCollision> selfContactEnergy;
      Contact::PointTrianglePairCouplingEnergyWithCollisionBuffer *selfContactEnergyBuf = nullptr;

      if (selfCD) {
        selfCD->execute(usurf.data());

        if (selfCD->getCollidingTrianglePair().size() > 0) {
          selfCD->handleContactDCD(0, 100);

          selfContactEnergy = selfCD->buildContactEnergy();
          selfContactEnergy->setToPosFunction([&restPosition](const ES::V3d &x, ES::V3d &p, int offset) {
            p = x + restPosition.segment<3>(offset);
          });

          selfContactEnergy->setToLastPosFunction([&restPosition, &u](const ES::V3d &x, ES::V3d &p, int offset) {
            p = restPosition.segment<3>(offset) + u.segment<3>(offset);
          });

          selfContactEnergyBuf = selfContactEnergy->allocateBuffer();
          selfContactEnergy->setBuffer(selfContactEnergyBuf);
          selfContactEnergy->setCoeff(contactK);
          selfContactEnergy->computeClosestPosition(u.data());

          selfContactEnergy->setFrictionCoeff(fricCoeff);
          selfContactEnergy->setTimestep(timestep);
          selfContactEnergy->setVelEps(velEps);

          intg->addGeneralImplicitForceModel(selfContactEnergy, 0, 0);
        }
      }

      intg->setqState(u, uvel, uacc);

      intg->doTimestep(1, 2, 1);

      intg->getq(u);
      intg->getqvel(uvel);
      intg->getqacc(uacc);

      if (extContactBuffer && extContactEnergy) {
        extContactEnergy->freeBuffer(extContactBuffer);
      }

      if (selfContactEnergyBuf && selfContactEnergy) {
        selfContactEnergy->freeBuffer(selfContactEnergyBuf);
      }

      ES::mv(W, u, usurf);

      if (framei % frameGap == 0) {
        ES::VXd psurf = surfaceRestPositions + usurf;

        Mesh::TriMeshGeo mesh = surfaceMesh;
        for (int vi = 0; vi < mesh.numVertices(); vi++) {
          mesh.pos(vi) = psurf.segment<3>(vi * 3) / scale;
        }
        mesh.save(fmt::format("{}/ret{:04d}.obj", outputFolder, framei / frameGap));
      }
    }
  }
  else if (simType == "static") {
    std::shared_ptr<PredefinedPotentialEnergies::LinearPotentialEnergy> externalForcesEnergy = std::make_shared<PredefinedPotentialEnergies::LinearPotentialEnergy>(fext);

    std::shared_ptr<NonlinearOptimization::PotentialEnergies> energyAll = std::make_shared<NonlinearOptimization::PotentialEnergies>(n3);
    energyAll->addPotentialEnergy(elasticEnergy);
    for (auto eng : pullingEnergies)
      energyAll->addPotentialEnergy(eng, 1.0);
    energyAll->addPotentialEnergy(externalForcesEnergy, -1.0);
    energyAll->init();

    NonlinearOptimization::NewtonRaphsonSolver::SolverParam solverParam;

    ES::VXd u(n3);
    u.setZero();

    energyAll->printEnergy(u);

    NonlinearOptimization::NewtonRaphsonSolver solver(u.data(), solverParam, energyAll, std::vector<int>(), nullptr);
    solver.solve(u.data(), solverMaxIter, solverEps, 2);

    ES::VXd x = restPosition + u;

    ES::VXd xsurf(surfn3);
    ES::mv(W, x, xsurf);

    Mesh::TriMeshGeo meshOut = surfaceMesh;
    for (int vi = 0; vi < meshOut.numVertices(); vi++) {
      meshOut.pos(vi) = xsurf.segment<3>(vi * 3);
    }
    meshOut.save(outputFolder);
  }

  return 0;
}