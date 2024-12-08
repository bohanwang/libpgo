#pragma once

#ifndef _LIBPGO_C_H_
#  define _LIBPGO_C_H_

#  include "pgo_c_def.h"

#  include <stdint.h>

LIBPGO_EXTERN_C_BEG

struct pgoTriMeshGeoStruct;
struct pgoTetMeshGeoStruct;
struct pgoTetMeshStruct;
struct pgoSmoothRSEnergyStruct;
struct pgoSparseMatrixStruct;

typedef struct pgoTriMeshGeoStruct *pgoTriMeshGeoStructHandle;
typedef struct pgoTetMeshGeoStruct *pgoTetMeshGeoStructHandle;
typedef struct pgoTetMeshStruct *pgoTetMeshStructHandle;
typedef struct pgoSmoothRSEnergyStruct *pgoSmoothRSEnergyStructHandle;
typedef struct pgoSparseMatrixStruct *pgoSparseMatrixStructHandle;

LIBPGO_C_EXPORT void pgo_init();

LIBPGO_C_EXPORT int pgo_load_1d_int_text(const char *filename, int **ptr, int sorted);

LIBPGO_C_EXPORT pgoTetMeshGeoStructHandle pgo_create_tetmeshgeo(int nv, double *vertices, int ntet, int *tets);
LIBPGO_C_EXPORT pgoTetMeshGeoStructHandle pgo_create_tetmeshgeo_from_file(char *filename);
LIBPGO_C_EXPORT int pgo_tetmeshgeo_get_num_vertices(pgoTetMeshGeoStructHandle tetmesh);
LIBPGO_C_EXPORT int pgo_tetmeshgeo_get_num_tets(pgoTetMeshGeoStructHandle tetmesh);
LIBPGO_C_EXPORT void pgo_tetmeshgeo_get_vertices(pgoTetMeshGeoStructHandle tetmesh, double *vertices);
LIBPGO_C_EXPORT void pgo_tetmeshgeo_get_tets(pgoTetMeshGeoStructHandle tetmesh, int *tets);
LIBPGO_C_EXPORT void pgo_destroy_tetmeshgeo(pgoTetMeshGeoStructHandle tetmesh);

LIBPGO_C_EXPORT pgoTetMeshStructHandle pgo_create_tetmesh_from_file(const char *filename);
LIBPGO_C_EXPORT pgoTetMeshStructHandle pgo_create_tetmesh(int nv, double *vertices, int ntet, int *tets, double E, double nu, double density);
LIBPGO_C_EXPORT void pgo_save_tetmesh_to_file(pgoTetMeshStructHandle tetMeshHandle, const char *filename);
LIBPGO_C_EXPORT int pgo_tetmesh_get_num_vertices(pgoTetMeshStructHandle m);
LIBPGO_C_EXPORT pgoTetMeshStructHandle pgo_tetmesh_update_vertices(pgoTetMeshStructHandle m, double *vertices);
LIBPGO_C_EXPORT int pgo_tetmesh_get_num_tets(pgoTetMeshStructHandle m);
LIBPGO_C_EXPORT void pgo_tetmesh_get_vertices(pgoTetMeshStructHandle m, double *vertices);
LIBPGO_C_EXPORT void pgo_tetmesh_get_elements(pgoTetMeshStructHandle m, int *elements);

LIBPGO_C_EXPORT int64_t pgo_sparse_matrix_get_num_entries(pgoSparseMatrixStructHandle m);
LIBPGO_C_EXPORT void pgo_sparse_matrix_get_row_indices(pgoSparseMatrixStructHandle m, int *rows);
LIBPGO_C_EXPORT void pgo_sparse_matrix_get_col_indices(pgoSparseMatrixStructHandle m, int *cols);
LIBPGO_C_EXPORT void pgo_sparse_matrix_get_values(pgoSparseMatrixStructHandle m, double *values);
LIBPGO_C_EXPORT void pgo_destroy_sparse_matrix(pgoSparseMatrixStructHandle m);

LIBPGO_C_EXPORT pgoSparseMatrixStructHandle pgo_create_tet_laplacian_matrix(pgoTetMeshGeoStructHandle tetmesh, int faceNeighbor, int repeat, int scale);
LIBPGO_C_EXPORT pgoSparseMatrixStructHandle pgo_create_tet_gradient_matrix(pgoTetMeshGeoStructHandle tetmesh);
LIBPGO_C_EXPORT void pgo_create_tet_gradient_per_element_matrix(pgoTetMeshGeoStructHandle tetmesh, double *smallMats);
LIBPGO_C_EXPORT pgoSparseMatrixStructHandle pgo_create_tet_biharmonic_gradient_matrix(pgoTetMeshGeoStructHandle tetmesh, int faceNeighbor, int scale);

LIBPGO_C_EXPORT pgoSmoothRSEnergyStructHandle pgo_create_smooth_rs_energy(pgoTetMeshGeoStructHandle tetmesh, double coeffR, double coeffS);
LIBPGO_C_EXPORT void pgo_destroy_smooth_rs_energy(pgoSmoothRSEnergyStructHandle energy);

LIBPGO_C_EXPORT double pgo_smooth_rs_energy_func(pgoSmoothRSEnergyStructHandle energy, double *x);
LIBPGO_C_EXPORT void pgo_smooth_rs_energy_grad(pgoSmoothRSEnergyStructHandle energy, double *x, double *grad);
LIBPGO_C_EXPORT int64_t pgo_smooth_rs_energy_hess_num_entries(pgoSmoothRSEnergyStructHandle energy);
LIBPGO_C_EXPORT void pgo_smooth_rs_energy_hess(pgoSmoothRSEnergyStructHandle energy, double *x, int *rows, int *cols, double *values);

LIBPGO_C_EXPORT double pgo_conjugate_mv(pgoSparseMatrixStructHandle m, double *v);
LIBPGO_C_EXPORT void pgo_sp_mv(pgoSparseMatrixStructHandle m, double *v, double *vout);

LIBPGO_C_EXPORT pgoTriMeshGeoStructHandle pgo_create_trimeshgeo(int nv, double *vertices, int ntri, int *tris);
LIBPGO_C_EXPORT int pgo_trimeshgeo_get_num_vertices(pgoTriMeshGeoStructHandle trimesh);
LIBPGO_C_EXPORT int pgo_trimeshgeo_get_num_triangles(pgoTriMeshGeoStructHandle trimesh);
LIBPGO_C_EXPORT void pgo_trimeshgeo_get_vertices(pgoTriMeshGeoStructHandle trimesh, double *vertices);
LIBPGO_C_EXPORT void pgo_trimeshgeo_get_triangles(pgoTriMeshGeoStructHandle trimesh, int *tris);
LIBPGO_C_EXPORT void pgo_destroy_trimeshgeo(pgoTriMeshGeoStructHandle trimesh);

LIBPGO_C_EXPORT void pgo_mesh_segmentation(pgoTriMeshGeoStructHandle trimesh, int numSegs, int *classIDs);
LIBPGO_C_EXPORT pgoTriMeshGeoStructHandle pgo_mesh_isotropic_remeshing(pgoTriMeshGeoStructHandle trimesh, double targetEdgeLength, int nIter, double angleThreshold);

LIBPGO_C_EXPORT void pgo_trimesh_closest_distances(pgoTriMeshGeoStructHandle trimesh, int n, double *queryPos, double *queryDistance, int *queryTri);
LIBPGO_C_EXPORT void pgo_tetmesh_barycentric_weights(pgoTetMeshGeoStructHandle tetmesh, int n, double *queryPos, double *queryW, int *queryEle);

LIBPGO_C_EXPORT int pgo_run_sim_from_config(const char *configFileName);

LIBPGO_EXTERN_C_END

#endif