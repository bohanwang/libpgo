import os
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

for sp_i in range(20):
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
print(pypgo.sparse_matrix_get_num_entries(L))

ridx = pypgo.sparse_matrix_get_row_indices(L)
cidx = pypgo.sparse_matrix_get_col_indices(L)
entries = pypgo.sparse_matrix_get_values(L)

print(ridx.shape)
print(cidx.shape)
print(entries.shape)

G_mats = pypgo.create_tet_gradient_per_element_matrix(tetmesh)

G0 = np.transpose(G_mats[0, :, :])
print(G0)

print(pypgo.run_sim_from_config("box.json"))

'''
x = pypgo.tetmeshgeo_get_vertices(tetmesh)
u = np.loadtxt("a.u")
x1 = u + x

ele = pypgo.tetmeshgeo_get_tets(tetmesh)

for ei in range(pypgo.tetmeshgeo_get_num_tets(tetmesh)):
    xlocal = np.zeros((12))
    for j in range(4):
        vid = ele[ei * 4 + j]
        print(vid)
        xlocal[j * 3 : (j + 1) * 3] = x1[vid * 3 : (vid + 1) * 3]

    print(xlocal)
    F = np.matmul(G0, xlocal).reshape((3, 3))

    print(F)
    exit(1)

row_idx = pypgo.pog_sparse_matrix_get_row_indices(L)
col_idx = pypgo.pog_sparse_matrix_get_col_indices(L)
values = pypgo.pog_sparse_matrix_get_values(L)

print(row_idx)
print(col_idx)
print(values)
'''