import os
import platform
import time

if "Windows" in platform.platform():
    os.add_dll_directory("c:/code/libpog-install")

import pypog
import numpy as np
import torch
from tet_spheres import tet_spheres_ext
import scipy

device = "cuda:0"
torch.manual_seed(0)

tetmesh = pypog.create_tetmeshgeo_from_file("a.veg")
v_load = pypog.tetmeshgeo_get_vertices(tetmesh)
t_load = pypog.tetmeshgeo_get_tets(tetmesh)

v = np.reshape(v_load, (-1, 3))
f = np.reshape(t_load, (-1, 4))

tet_sp = tet_spheres_ext.TetSpheres("a.veg")
rand_x = tet_spheres_ext.random_x(tet_sp)
rand_x = rand_x.to(device=device)

print(rand_x.shape)

v_rest = torch.from_numpy(v).to(device=device)
gradH = torch.tensor(1.0, dtype=torch.float32).to(device=device)


def func(x):
    x_v = torch.from_numpy(x).view(-1, 3).to(device=device, dtype=torch.float32)
    ret = tet_spheres_ext.forward(x_v, tet_sp, 1.0, 0.0, 2)

    return ret.item()


def grad(x):
    x_v = torch.from_numpy(x).view(-1, 3).to(device=device, dtype=torch.float32)
    g = tet_spheres_ext.backward(gradH, x_v, tet_sp, 1.0, 0.0, 2)

    return g.cpu().numpy().flatten()


# x0 = rand_x.cpu().numpy().flatten()
# print(x0[0])

# g = tet_spheres_ext.backward(gradH, rand_x, tet_sp, 1.0, 0.0, 2)
# eps = 4e-5
# x1 = np.copy(x0)
# x1[0] += eps

# x2 = np.copy(x0)
# x2[0] -= eps
# v1 = func(x1)
# v2 = func(x2)

# print(f"{v1},{v2},{(v1 - v2) * 0.5 / eps}")
# print(g[0])


GTLTLG = pypog.create_tet_biharmonic_gradient_matrix(tetmesh, 1, 0)
row_idx = pypog.sparse_matrix_get_row_indices(GTLTLG)
col_idx = pypog.sparse_matrix_get_col_indices(GTLTLG)
values = pypog.sparse_matrix_get_values(GTLTLG)

n_ele = pypog.tetmeshgeo_get_num_tets(tetmesh)
n = pypog.tetmeshgeo_get_num_vertices(tetmesh)
GTLTLG_torch = (
    torch.sparse_coo_tensor(np.array([row_idx, col_idx]), values, (n * 3, n * 3))
    .coalesce()
    .to(device=device)
)

G = pypog.create_tet_gradient_matrix(tetmesh)
row_idx = pypog.sparse_matrix_get_row_indices(G)
col_idx = pypog.sparse_matrix_get_col_indices(G)
values = pypog.sparse_matrix_get_values(G)

G_torch = (
    torch.sparse_coo_tensor(np.array([row_idx, col_idx]), values, (n_ele * 9, n * 3))
    .coalesce()
    .to(device=device)
)

c1 = 1.0
c2 = 2.0

t1 = time.time()

for i in range(1000):
    f0 = tet_spheres_ext.forward(rand_x, tet_sp, c1, c2, 2)
    g0 = tet_spheres_ext.backward(gradH, rand_x, tet_sp, c1, c2, 2)

t2 = time.time()

print((t2 - t1) / 1000)

print("me:")
print(f0)
print(g0)

x = rand_x.flatten().requires_grad_()

t1 = time.time()

for i in range(1000):
    xtemp = torch.mv(GTLTLG_torch, x)
    f1 = torch.dot(xtemp, x) * 0.5 * c1

    Ftemp = torch.mv(G_torch, x)
    F_mat = Ftemp.reshape(-1, 3, 3)
    raw_J = torch.det(F_mat)
    J = torch.nn.functional.relu(-raw_J)
    f1 += (J * J).sum() * c2
    f1.backward()

    x.grad.data.zero_()
    print(x.grad)
    print(x.grad.data)

t2 = time.time()

print((t2 - t1) / 1000)

print("torch:")
print(x.grad)
print(f1)


exit(1)

print("=====================================")


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

tet_sp = tet_spheres_ext.TetSpheres(v_flat, f_flat)
rand_x = tet_spheres_ext.random_x(tet_sp)
rand_x = rand_x.to(device=device)

print(rand_x.shape)

# t1 = time.time()
# for i in range(1000):
#     tet_spheres_ext.forward(rand_x, tet_sp, 1.0, 1.0)

# t2 = time.time()

# print((t2 - t1) / 1000.0)

v_rest = torch.from_numpy(v).to(device=device)
tet_spheres_ext.forward(v_rest, tet_sp, 1.0, 1.0, 2)
