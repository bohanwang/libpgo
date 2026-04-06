#include "pgo_c.h"

#include "fileService.h"
#include "basicIO.h"
#include "runSimCore.h"
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
#include "configFileJSON.h"

#if defined(PGO_HAS_ANIMATIONIO)
#include "animationLoader.h"
#endif

#if defined(PGO_HAS_MKL)
#include "smoothRSEnergy.h"
#endif

#if defined(PGO_HAS_CGAL)
#include "cgalInterface.h"
#endif

#include <tbb/parallel_for.h>

#include <fmt/format.h>

pgoTetMeshGeoStructHandle pgo_create_tetmeshgeo(int nv, double* vertices, int ntet, int* tets) {
    pgo::Mesh::TetMeshGeo* tetmesh = new pgo::Mesh::TetMeshGeo(nv, vertices, ntet, tets);
    return reinterpret_cast<pgoTetMeshGeoStructHandle>(tetmesh);
}

pgoTetMeshGeoStructHandle pgo_create_tetmeshgeo_from_file(char* filename) {
    try {
        pgo::VolumetricMeshes::TetMesh mesh(filename);
        pgo::Mesh::TetMeshGeo*         tetmesh = new pgo::Mesh::TetMeshGeo;
        mesh.exportMeshGeometry(*tetmesh);

        std::cout << "#vtx:" << tetmesh->numVertices() << std::endl;
        std::cout << "#tets:" << tetmesh->numTets() << std::endl;

        return reinterpret_cast<pgoTetMeshGeoStructHandle>(tetmesh);
    } catch (...) {
        return nullptr;
    }
}

int pgo_tetmeshgeo_get_num_vertices(pgoTetMeshGeoStructHandle m) {
    pgo::Mesh::TetMeshGeo* tetmesh = reinterpret_cast<pgo::Mesh::TetMeshGeo*>(m);

    return tetmesh->numVertices();
}

int pgo_tetmeshgeo_get_num_tets(pgoTetMeshGeoStructHandle m) {
    pgo::Mesh::TetMeshGeo* tetmesh = reinterpret_cast<pgo::Mesh::TetMeshGeo*>(m);

    return tetmesh->numTets();
}

void pgo_tetmeshgeo_get_vertices(pgoTetMeshGeoStructHandle m, double* vertices) {
    pgo::Mesh::TetMeshGeo* tetmesh = reinterpret_cast<pgo::Mesh::TetMeshGeo*>(m);

    for (int vi = 0; vi < tetmesh->numVertices(); vi++) {
        vertices[vi * 3]     = tetmesh->pos(vi)[0];
        vertices[vi * 3 + 1] = tetmesh->pos(vi)[1];
        vertices[vi * 3 + 2] = tetmesh->pos(vi)[2];
    }
}

void pgo_tetmeshgeo_get_tets(pgoTetMeshGeoStructHandle m, int* tets) {
    pgo::Mesh::TetMeshGeo* tetmesh = reinterpret_cast<pgo::Mesh::TetMeshGeo*>(m);

    for (int ti = 0; ti < tetmesh->numTets(); ti++) {
        tets[ti * 4]     = tetmesh->tet(ti)[0];
        tets[ti * 4 + 1] = tetmesh->tet(ti)[1];
        tets[ti * 4 + 2] = tetmesh->tet(ti)[2];
        tets[ti * 4 + 3] = tetmesh->tet(ti)[3];
    }

    std::cout << tets[0] << ',' << tets[1] << ',' << tets[2] << ',' << tets[3] << std::endl;
}

void pgo_destroy_tetmeshgeo(pgoTetMeshGeoStructHandle m) {
    pgo::Mesh::TetMeshGeo* tetmesh = reinterpret_cast<pgo::Mesh::TetMeshGeo*>(m);
    delete tetmesh;
}

pgoTetMeshStructHandle pgo_create_tetmesh(int nv, double* vertices, int ntet, int* tets, double E, double nu,
                                          double density) {
    try {
        pgo::VolumetricMeshes::TetMesh* mesh =
            new pgo::VolumetricMeshes::TetMesh(nv, vertices, ntet, tets, E, nu, density);
        return reinterpret_cast<pgoTetMeshStructHandle>(mesh);
    } catch (...) {
        return nullptr;
    }
}

pgoTetMeshStructHandle pgo_create_tetmesh_from_file(const char* filename) {
    try {
        pgo::VolumetricMeshes::TetMesh* mesh = new pgo::VolumetricMeshes::TetMesh(filename);
        return reinterpret_cast<pgoTetMeshStructHandle>(mesh);
    } catch (...) {
        return nullptr;
    }
}

void pgo_save_tetmesh_to_file(pgoTetMeshStructHandle tetMeshHandle, const char* filename) {
    pgo::VolumetricMeshes::TetMesh* mesh = reinterpret_cast<pgo::VolumetricMeshes::TetMesh*>(tetMeshHandle);
    mesh->save(filename);
}

int pgo_tetmesh_get_num_vertices(pgoTetMeshStructHandle m) {
    pgo::VolumetricMeshes::TetMesh* mesh = reinterpret_cast<pgo::VolumetricMeshes::TetMesh*>(m);
    pgo::Mesh::TetMeshGeo           tetmeshGeo;
    mesh->exportMeshGeometry(tetmeshGeo);
    return tetmeshGeo.numVertices();
}

int pgo_tetmesh_get_num_tets(pgoTetMeshStructHandle m) {
    pgo::VolumetricMeshes::TetMesh* mesh = reinterpret_cast<pgo::VolumetricMeshes::TetMesh*>(m);
    pgo::Mesh::TetMeshGeo           tetmeshGeo;
    mesh->exportMeshGeometry(tetmeshGeo);

    return tetmeshGeo.numTets();
}

void pgo_tetmesh_get_vertices(pgoTetMeshStructHandle m, double* vertices) {
    pgo::VolumetricMeshes::TetMesh* mesh = reinterpret_cast<pgo::VolumetricMeshes::TetMesh*>(m);
    pgo::Mesh::TetMeshGeo           tetmeshGeo;
    mesh->exportMeshGeometry(tetmeshGeo);

    for (int vi = 0; vi < tetmeshGeo.numVertices(); vi++) {
        vertices[vi * 3]     = tetmeshGeo.pos(vi)[0];
        vertices[vi * 3 + 1] = tetmeshGeo.pos(vi)[1];
        vertices[vi * 3 + 2] = tetmeshGeo.pos(vi)[2];
    }
}

pgoTetMeshStructHandle pgo_tetmesh_update_vertices(pgoTetMeshStructHandle m, double* vertices) {
    pgo::VolumetricMeshes::TetMesh* tetMesh = reinterpret_cast<pgo::VolumetricMeshes::TetMesh*>(m);

    pgo::VolumetricMeshes::TetMesh* tetMeshNew = new pgo::VolumetricMeshes::TetMesh(*tetMesh);

    for (int vi = 0; vi < tetMeshNew->getNumVertices(); vi++) {
        pgo::Vec3d p(vertices + vi * 3);
        tetMeshNew->setVertex(vi, p);
    }
    return reinterpret_cast<pgoTetMeshStructHandle>(tetMeshNew);
}

void pgo_tetmesh_get_elements(pgoTetMeshStructHandle m, int* elements) {
    pgo::VolumetricMeshes::TetMesh* mesh = reinterpret_cast<pgo::VolumetricMeshes::TetMesh*>(m);
    pgo::Mesh::TetMeshGeo           tetmeshGeo;
    mesh->exportMeshGeometry(tetmeshGeo);

    for (int ti = 0; ti < tetmeshGeo.numTets(); ti++) {
        for (int i = 0; i < 4; i++) {
            elements[ti * 4 + i] = tetmeshGeo.tet(ti)[i];
        }
    }
}

pgoSmoothRSEnergyStructHandle pgo_create_smooth_rs_energy(pgoTetMeshGeoStructHandle m, double coeffR, double coeffS) {
#if defined(PGO_HAS_MKL)
    pgo::Mesh::TetMeshGeo* tetmesh = reinterpret_cast<pgo::Mesh::TetMeshGeo*>(m);
    std::cout << tetmesh->numVertices() << '\n' << tetmesh->numTets() << std::endl;

    pgo::EigenSupport::SpMatD G;
    pgo::SolidDeformationModel::TetMeshMatrix::generateGradientMatrix(*tetmesh, G);

    std::cout << "G done." << std::endl;

    pgo::EigenSupport::SpMatD L;
    pgo::SolidDeformationModel::TetMeshMatrix::generateBasicElementLaplacianMatrix(*tetmesh, L, 0, 0);

    std::cout << "L done." << std::endl;

    double                                            coeffs[2] = {coeffR, coeffS};
    pgo::PredefinedPotentialEnergies::SmoothRSEnergy* energy =
        new pgo::PredefinedPotentialEnergies::SmoothRSEnergy(*tetmesh, G, L, coeffs);

    std::cout << "energy done." << std::endl;

    return reinterpret_cast<pgoSmoothRSEnergyStructHandle>(energy);
#else
    std::cerr << "No MKL. Unsupported." << std::endl;
    return nullptr;
#endif
}

void pgo_destroy_smooth_rs_energy(pgoSmoothRSEnergyStructHandle energy) {
#if defined(PGO_HAS_MKL)
    pgo::PredefinedPotentialEnergies::SmoothRSEnergy* eng =
        reinterpret_cast<pgo::PredefinedPotentialEnergies::SmoothRSEnergy*>(energy);
    delete eng;
#else
    std::cerr << "No MKL. Unsupported." << std::endl;
#endif
}

double pgo_smooth_rs_energy_func(pgoSmoothRSEnergyStructHandle energy, double* x) {
#if defined(PGO_HAS_MKL)
    namespace ES = pgo::EigenSupport;
    pgo::PredefinedPotentialEnergies::SmoothRSEnergy* eng =
        reinterpret_cast<pgo::PredefinedPotentialEnergies::SmoothRSEnergy*>(energy);
    return eng->func(ES::Mp<const ES::VXd>(x, eng->getNumDOFs()));
#else
    std::cerr << "No MKL. Unsupported." << std::endl;
    return 0;
#endif
}

void pgo_smooth_rs_energy_grad(pgoSmoothRSEnergyStructHandle energy, double* x, double* grad) {
#if defined(PGO_HAS_MKL)
    namespace ES = pgo::EigenSupport;
    pgo::PredefinedPotentialEnergies::SmoothRSEnergy* eng =
        reinterpret_cast<pgo::PredefinedPotentialEnergies::SmoothRSEnergy*>(energy);
    eng->gradient(ES::Mp<const ES::VXd>(x, eng->getNumDOFs()), ES::Mp<ES::VXd>(grad, eng->getNumDOFs()));
#else
    std::cerr << "No MKL. Unsupported." << std::endl;
#endif
}

int64_t pgo_smooth_rs_energy_hess_num_entries(pgoSmoothRSEnergyStructHandle energy) {
#if defined(PGO_HAS_MKL)
    namespace ES = pgo::EigenSupport;
    pgo::PredefinedPotentialEnergies::SmoothRSEnergy* eng =
        reinterpret_cast<pgo::PredefinedPotentialEnergies::SmoothRSEnergy*>(energy);

    ES::SpMatD H;
    eng->createHessian(H);

    return (int64_t)H.nonZeros();
#else
    std::cerr << "No MKL. Unsupported." << std::endl;
    return 0;
#endif
}

void pgo_smooth_rs_energy_hess(pgoSmoothRSEnergyStructHandle energy, double* x, int* rows, int* cols, double* values) {
#if defined(PGO_HAS_MKL)
    namespace ES = pgo::EigenSupport;
    pgo::PredefinedPotentialEnergies::SmoothRSEnergy* eng =
        reinterpret_cast<pgo::PredefinedPotentialEnergies::SmoothRSEnergy*>(energy);

    ES::SpMatD H;
    eng->createHessian(H);

    memset(H.valuePtr(), 0, H.nonZeros() * sizeof(double));
    eng->hessian(ES::Mp<const ES::VXd>(x, eng->getNumDOFs()), H);

    int64_t inc = 0;
    for (ES::IDX rowi = 0; rowi < H.rows(); rowi++) {
        for (ES::SpMatD::InnerIterator it(H, rowi); it; ++it) {
            rows[inc]   = (int)it.row();
            cols[inc]   = (int)it.col();
            values[inc] = (int)it.value();
            inc++;
        }
    }
#else
    std::cerr << "No MKL. Unsupported." << std::endl;
#endif
}

int64_t pgo_sparse_matrix_get_num_entries(pgoSparseMatrixStructHandle m) {
    namespace ES    = pgo::EigenSupport;
    ES::SpMatD* mat = reinterpret_cast<ES::SpMatD*>(m);

    return mat->nonZeros();
}

void pgo_sparse_matrix_get_row_indices(pgoSparseMatrixStructHandle m, int* rows) {
    namespace ES = pgo::EigenSupport;

    ES::SpMatD* mat = reinterpret_cast<ES::SpMatD*>(m);

    int64_t inc = 0;
    for (ES::IDX rowi = 0; rowi < mat->rows(); rowi++) {
        for (ES::SpMatD::InnerIterator it(*mat, rowi); it; ++it) {
            rows[inc] = (int)it.row();
            inc++;
        }
    }
}

void pgo_sparse_matrix_get_col_indices(pgoSparseMatrixStructHandle m, int* cols) {
    namespace ES    = pgo::EigenSupport;
    ES::SpMatD* mat = reinterpret_cast<ES::SpMatD*>(m);

    int64_t inc = 0;
    for (ES::IDX rowi = 0; rowi < mat->rows(); rowi++) {
        for (ES::SpMatD::InnerIterator it(*mat, rowi); it; ++it) {
            cols[inc] = (int)it.col();
            inc++;
        }
    }
}

void pgo_sparse_matrix_get_values(pgoSparseMatrixStructHandle m, double* values) {
    namespace ES    = pgo::EigenSupport;
    ES::SpMatD* mat = reinterpret_cast<ES::SpMatD*>(m);

    int64_t inc = 0;
    for (ES::IDX rowi = 0; rowi < mat->rows(); rowi++) {
        for (ES::SpMatD::InnerIterator it(*mat, rowi); it; ++it) {
            values[inc] = it.value();
            inc++;
        }
    }
}

void pgo_destroy_sparse_matrix(pgoSparseMatrixStructHandle m) {
    namespace ES    = pgo::EigenSupport;
    ES::SpMatD* mat = reinterpret_cast<ES::SpMatD*>(m);

    delete mat;
}

pgoSparseMatrixStructHandle pgo_create_tet_laplacian_matrix(pgoTetMeshGeoStructHandle m, int faceNeighbor, int repeat,
                                                            int scale) {
    namespace ES                   = pgo::EigenSupport;
    pgo::Mesh::TetMeshGeo* tetmesh = reinterpret_cast<pgo::Mesh::TetMeshGeo*>(m);

    ES::SpMatD L;
    pgo::SolidDeformationModel::TetMeshMatrix::generateBasicElementLaplacianMatrix(*tetmesh, L, faceNeighbor, scale);

    if (repeat > 1) {
        ES::SpMatD Ln;
        ES::expandN(L, Ln, repeat);

        L = Ln;
    }

    ES::SpMatD* ret = new ES::SpMatD(L);
    return reinterpret_cast<pgoSparseMatrixStructHandle>(ret);
}

pgoSparseMatrixStructHandle pgo_create_tet_gradient_matrix(pgoTetMeshGeoStructHandle m) {
    namespace ES                   = pgo::EigenSupport;
    pgo::Mesh::TetMeshGeo* tetmesh = reinterpret_cast<pgo::Mesh::TetMeshGeo*>(m);

    ES::SpMatD G;
    pgo::SolidDeformationModel::TetMeshMatrix::generateGradientMatrix(*tetmesh, G);

    return reinterpret_cast<pgoSparseMatrixStructHandle>(new ES::SpMatD(G));
}

void pgo_create_tet_gradient_per_element_matrix(pgoTetMeshGeoStructHandle m, double* smallMats) {
    namespace ES                   = pgo::EigenSupport;
    pgo::Mesh::TetMeshGeo* tetmesh = reinterpret_cast<pgo::Mesh::TetMeshGeo*>(m);

    for (int i = 0; i < tetmesh->numTets(); i++) {
        ES::M9x12d G;
        pgo::SolidDeformationModel::TetMeshMatrix::generateElementGradientMatrix(*tetmesh, i, G);

        ES::Mp<ES::M9x12d>(smallMats + i * 9 * 12) = G;
    }
}

pgoSparseMatrixStructHandle pgo_create_tet_biharmonic_gradient_matrix(pgoTetMeshGeoStructHandle m, int faceNeighbor,
                                                                      int scale) {
    namespace ES                   = pgo::EigenSupport;
    pgo::Mesh::TetMeshGeo* tetmesh = reinterpret_cast<pgo::Mesh::TetMeshGeo*>(m);

    ES::SpMatD L;
    pgo::SolidDeformationModel::TetMeshMatrix::generateBasicElementLaplacianMatrix(*tetmesh, L, faceNeighbor, 0);

    ES::SpMatD Ln;
    ES::expandN(L, Ln, 9);

    ES::SpMatD G;
    pgo::SolidDeformationModel::TetMeshMatrix::generateGradientMatrix(*tetmesh, G);

    ES::SpMatD LG;
    ES::mm(Ln, G, LG);

    ES::SpMatD* ret    = new ES::SpMatD;
    ES::SpMatD& GTLTLG = *ret;
    if (scale) {
        std::vector<double> eleVols(tetmesh->numTets());
        for (int ei = 0; ei < tetmesh->numTets(); ei++) {
            eleVols[ei] = std::abs(pgo::Mesh::getTetDeterminant(tetmesh->pos(ei, 0), tetmesh->pos(ei, 1),
                                                                tetmesh->pos(ei, 2), tetmesh->pos(ei, 3)));
        }

        double maxVol = *std::max_element(eleVols.begin(), eleVols.end());
        for (int ei = 0; ei < tetmesh->numTets(); ei++) {
            eleVols[ei] /= maxVol;
        }

        ES::SpMatD                M;
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
    } else {
        ES::mm(LG, LG, GTLTLG, 1);
    }

    return reinterpret_cast<pgoSparseMatrixStructHandle>(ret);
}

double pgo_conjugate_mv(pgoSparseMatrixStructHandle m, double* v) {
    namespace ES      = pgo::EigenSupport;
    ES::SpMatD*     M = reinterpret_cast<ES::SpMatD*>(m);
    ES::Mp<ES::VXd> vmap(v, M->rows());

    return vmap.dot(*M * vmap) * 0.5;
}

void pgo_sp_mv(pgoSparseMatrixStructHandle m, double* v, double* vout) {
    namespace ES      = pgo::EigenSupport;
    ES::SpMatD*     M = reinterpret_cast<ES::SpMatD*>(m);
    ES::Mp<ES::VXd> vmap(v, M->rows());
    ES::Mp<ES::VXd> vmap_out(vout, M->rows());

    ES::mv(*M, vmap, vmap_out);
}

void pgo_init() {
    pgo::Logging::init();
}

pgoTriMeshGeoStructHandle pgo_create_trimeshgeo(int nv, double* vertices, int ntri, int* tris) {
    pgo::Mesh::TriMeshGeo* mesh = new pgo::Mesh::TriMeshGeo(nv, vertices, ntri, tris);
    return reinterpret_cast<pgoTriMeshGeoStructHandle>(mesh);
}

int pgo_trimeshgeo_get_num_vertices(pgoTriMeshGeoStructHandle trimesh) {
    pgo::Mesh::TriMeshGeo* mesh = reinterpret_cast<pgo::Mesh::TriMeshGeo*>(trimesh);
    return mesh->numVertices();
}

int pgo_trimeshgeo_get_num_triangles(pgoTriMeshGeoStructHandle trimesh) {
    pgo::Mesh::TriMeshGeo* mesh = reinterpret_cast<pgo::Mesh::TriMeshGeo*>(trimesh);
    return mesh->numTriangles();
}

void pgo_trimeshgeo_get_vertices(pgoTriMeshGeoStructHandle trimesh, double* vertices) {
    pgo::Mesh::TriMeshGeo* mesh = reinterpret_cast<pgo::Mesh::TriMeshGeo*>(trimesh);
    memcpy(vertices, mesh->positions()[0].data(), sizeof(double) * mesh->positions().size() * 3);
}

void pgo_trimeshgeo_get_triangles(pgoTriMeshGeoStructHandle trimesh, int* tris) {
    pgo::Mesh::TriMeshGeo* mesh = reinterpret_cast<pgo::Mesh::TriMeshGeo*>(trimesh);
    memcpy(tris, mesh->triangles()[0].data(), sizeof(int) * mesh->triangles().size() * 3);
}

void pgo_destroy_trimeshgeo(pgoTriMeshGeoStructHandle trimesh) {
    pgo::Mesh::TriMeshGeo* mesh = reinterpret_cast<pgo::Mesh::TriMeshGeo*>(trimesh);
    delete mesh;
}

void pgo_mesh_segmentation(pgoTriMeshGeoStructHandle trimesh, int numSegs, int* classIDs) {
#if defined(PGO_HAS_CGAL)
    pgo::Mesh::TriMeshGeo* mesh = reinterpret_cast<pgo::Mesh::TriMeshGeo*>(trimesh);
    std::vector<int>       cids;
// pgo::CGALInterface::segmentMesh(*mesh, numSegs, cids);
// memcpy(classIDs, cids.data(), cids.size() * sizeof(int));
#else
    std::cout << "No CGAL. Unsupported." << std::endl;
#endif
}

pgoTriMeshGeoStructHandle pgo_mesh_isotropic_remeshing(pgoTriMeshGeoStructHandle trimesh, double targetEdgeLength,
                                                       int nIter, double angleThreshold) {
#if defined(PGO_HAS_CGAL)
    pgo::Mesh::TriMeshGeo* mesh = reinterpret_cast<pgo::Mesh::TriMeshGeo*>(trimesh);
    pgo::Mesh::TriMeshGeo  meshOut =
        pgo::CGALInterface::isotropicRemeshing(*mesh, targetEdgeLength, nIter, angleThreshold);
    return reinterpret_cast<pgoTriMeshGeoStructHandle>(new pgo::Mesh::TriMeshGeo(meshOut));
#else
    std::cout << "No CGAL. Unsupported." << std::endl;
    return nullptr;
#endif
}

void pgo_trimesh_closest_distances(pgoTriMeshGeoStructHandle trimesh, int n, double* queryPos, double* queryDistance,
                                   int* queryTri) {
    pgo::Mesh::initPredicates();

    pgo::Mesh::TriMeshGeo*   mesh = reinterpret_cast<pgo::Mesh::TriMeshGeo*>(trimesh);
    pgo::Mesh::TriMeshBVTree bvTree;
    bvTree.buildByInertiaPartition(*mesh);

    tbb::parallel_for(0, n, [&](int i) {
        pgo::Vec3d pt(queryPos + i * 3);
        auto       ret   = bvTree.closestTriangleQuery(*mesh, pt);
        queryDistance[i] = ret.dist2;

        if (queryTri) {
            queryTri[i] = ret.triID;
        }
    });
}

void pgo_tetmesh_barycentric_weights(pgoTetMeshGeoStructHandle tetmesh, int n, double* queryPos, double* queryW,
                                     int* queryEle) {
    pgo::Mesh::initPredicates();

    pgo::Mesh::TetMeshGeo*   mesh = reinterpret_cast<pgo::Mesh::TetMeshGeo*>(tetmesh);
    pgo::Mesh::TetMeshBVTree bvTree;
    bvTree.buildByInertiaPartition(*mesh);

    tbb::parallel_for(0, n, [&](int i) {
        pgo::Vec3d pt(queryPos + i * 3);
        int        ele = bvTree.getClosestTet(*mesh, pt);
        queryEle[i]    = ele;
        pgo::Mesh::getTetBarycentricWeights(pt, mesh->pos(ele, 0), mesh->pos(ele, 1), mesh->pos(ele, 2),
                                            mesh->pos(ele, 3), queryW + i * 4);
    });
}

namespace {
int runSimFromConfigImpl(const char* configFileName, int deterministicOverride) {
    using namespace pgo;

    pgo::Logging::init();
    pgo::Mesh::initPredicates();

    const pgo::FileService::ConfigPathResolver configPathResolver(configFileName);

    ConfigFileJSON jconfig;
    if (jconfig.open(configPathResolver.filePath().c_str()) != true) {
        return 0;
    }

    try {
        pgo::api::RunSimConfig config = pgo::api::parseRunSimConfig(jconfig, configPathResolver.filePath());
        if (deterministicOverride >= 0) {
            config.runtime.deterministicMode = (deterministicOverride != 0);
        }
        return pgo::api::runSimFromConfig(config);
    } catch (const std::exception& err) {
        SPDLOG_LOGGER_ERROR(Logging::lgr(), "{}", err.what());
        return 1;
    }
}
}  // namespace

int pgo_run_sim_from_config(const char* configFileName) {
    return runSimFromConfigImpl(configFileName, -1);
}

int pgo_run_sim_from_config_with_deterministic(const char* configFileName, int deterministic) {
    return runSimFromConfigImpl(configFileName, deterministic ? 1 : 0);
}

int pgo_convert_animation_to_abc(const char* configFileName, const char* outputFolder) {
#if defined(PGO_HAS_ANIMATIONIO)
    pgo::Mesh::initPredicates();
    pgo::AnimationIO::AnimationLoader loader;
    if (loader.load(configFileName) != 0) {
        return 1;
    }

    return loader.saveABC(outputFolder);
#else
    (void)configFileName;
    (void)outputFolder;
    return 1;
#endif
}
