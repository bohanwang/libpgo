import pypgo
import numpy as np

tetmesh = pypgo.create_tetmeshgeo_from_file("torus.veg")
v_load = pypgo.tetmeshgeo_get_vertices(tetmesh)
t_load = pypgo.tetmeshgeo_get_tets(tetmesh)
pypgo.destroy_tetmeshgeo(tetmesh)

v = np.reshape(v_load, (-1, 3))
f = np.reshape(t_load, (-1, 4))

all_v_param = []
all_surface_f = []
all_surface_vid = []
all_tet_faces = []

for sp_i in range(10):
    v_trans = v
    all_v_param.append(v_trans)
    all_tet_faces.append(f + sp_i * len(v))

all_v_param = np.concatenate(all_v_param, axis=0)
all_tet_faces = np.concatenate(all_tet_faces, axis=0)

v = all_v_param
f = all_tet_faces

v_flat = v.flatten().astype(np.float32)
f_flat = f.flatten().astype(np.int32)
tetmesh = pypgo.create_tetmeshgeo(v_flat, f_flat)

L = pypgo.create_element_laplacian_matrix(tetmesh, 0, 9, 0)

ridx = pypgo.sparse_matrix_get_row_indices(L)
cidx = pypgo.sparse_matrix_get_col_indices(L)
entries = pypgo.sparse_matrix_get_values(L)

print("L Info:")
print(pypgo.sparse_matrix_get_num_entries(L))
print(ridx.shape)
print(cidx.shape)
print(np.amax(entries))

pypgo.destroy_sparse_matrix(L)

GTLTLG = pypgo.create_tet_biharmonic_gradient_matrix(tetmesh, 1, 0)

ridx = pypgo.sparse_matrix_get_row_indices(GTLTLG)
cidx = pypgo.sparse_matrix_get_col_indices(GTLTLG)
entries = pypgo.sparse_matrix_get_values(GTLTLG)

print("GTLTLG Info:")
print(pypgo.sparse_matrix_get_num_entries(GTLTLG))
print(ridx.shape)
print(cidx.shape)
print(np.amax(entries))

pypgo.destroy_sparse_matrix(GTLTLG)

G_mats = pypgo.create_tet_gradient_per_element_matrix(tetmesh)

G0 = np.transpose(G_mats[0, :, :])
print(G0)

pypgo.destroy_tetmeshgeo(tetmesh)
