#include "pypgo.h"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <Eigen/Dense>

#include <iostream>
#include <chrono>
#include <fstream>

namespace py = pybind11;

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

void pypgo_init(py::module &m)
{
  py::class_<TetMeshGeo> TetMeshGeo_class(m, "TetMeshGeo");
  TetMeshGeo_class.def(py::init<>());

  py::class_<TriMeshGeo> TriMeshGeo_class(m, "TriMeshGeo");
  TriMeshGeo_class.def(py::init<>());

  py::class_<TetMesh> TetMesh_class(m, "TetMesh");
  TetMesh_class.def(py::init<>());

  py::class_<SmoothRSEnergy> SmoothRSEnergy_class(m, "SmoothRSEnergy");
  SmoothRSEnergy_class.def(py::init<>());

  py::class_<SparseMatrix> SparseMatrix_class(m, "SparseMatrix");
  SparseMatrix_class.def(py::init<>());

  using pyArrayFloat = py::array_t<float, py::array::c_style | py::array::forcecast>;
  using pyArrayInt = py::array_t<int, py::array::c_style | py::array::forcecast>;

  m.def(
    "create_tetmeshgeo", [](pyArrayFloat vertices, pyArrayInt elements) -> TetMeshGeo {
      py::buffer_info vtxInfo = vertices.request();
      py::buffer_info tetInfo = elements.request();

      if (vtxInfo.ndim != (py::ssize_t)1 || vtxInfo.format != py::format_descriptor<float>::format()) {
        std::cerr << "Wrong vertex type:" << vtxInfo.ndim << ',' << vtxInfo.format << std::endl;
        return TetMeshGeo();
      }

#if defined(_WIN32)
      if (tetInfo.ndim != (py::ssize_t)1) {
#else
      if (tetInfo.ndim != (py::ssize_t)1 || (tetInfo.format != py::format_descriptor<int>::format())) {
#endif
        std::cerr << "Wrong tet type:" << tetInfo.ndim << ',' << tetInfo.format << ',' << py::format_descriptor<int64_t>::format() << std::endl;
        return TetMeshGeo();
      }

      Eigen::VectorXd vertexPosDouble = Eigen::Map<const Eigen::VectorXf>((float *)vtxInfo.ptr, vtxInfo.size).cast<double>();
      Eigen::VectorXi tets = Eigen::Map<const Eigen::VectorXi>((int *)tetInfo.ptr, tetInfo.size);

      pgoTetMeshGeoStructHandle tetmesh = pgo_create_tetmeshgeo((int)vertexPosDouble.size() / 3, vertexPosDouble.data(), (int)tets.size() / 4, tets.data());
      return TetMeshGeo(tetmesh);
    });

  m.def(
    "create_tetmeshgeo_from_file", [](const std::string &filename) -> TetMeshGeo {
      pgoTetMeshGeoStructHandle tetmesh = pgo_create_tetmeshgeo_from_file(const_cast<char *>(filename.c_str()));
      return TetMeshGeo(tetmesh);
    });

  m.def("tetmeshgeo_get_num_vertices", [](TetMeshGeo m) -> int {
    return pgo_tetmeshgeo_get_num_vertices(m.handle);
  });

  m.def("tetmeshgeo_get_num_tets", [](TetMeshGeo m) -> int {
    return pgo_tetmeshgeo_get_num_tets(m.handle);
  });

  m.def("tetmeshgeo_get_vertices", [](TetMeshGeo m) -> py::array_t<float> {
    Eigen::VectorXd vtx(pgo_tetmeshgeo_get_num_vertices(m.handle) * 3);
    pgo_tetmeshgeo_get_vertices(m.handle, vtx.data());

    auto ret = py::array_t<float>(vtx.size());
    py::buffer_info binfo = ret.request();
    (Eigen::Map<Eigen::VectorXf>((float *)binfo.ptr, vtx.size())) = vtx.cast<float>();
    return ret;
  });

  m.def("tetmeshgeo_get_tets", [](TetMeshGeo m) -> py::array_t<int> {
    auto ret = py::array_t<int>(pgo_tetmeshgeo_get_num_tets(m.handle) * 4);
    py::buffer_info binfo = ret.request();

    pgo_tetmeshgeo_get_tets(m.handle, (int *)binfo.ptr);

    return ret;
  });

  m.def("destroy_tetmeshgeo", [](TetMeshGeo tetmesh) {
    pgo_destroy_tetmeshgeo(tetmesh.handle);
  });

  m.def("create_tetmesh", [](pyArrayFloat vertices, pyArrayInt elements, float E, float nu, float density) -> TetMesh {
    py::buffer_info vtxInfo = vertices.request();
    py::buffer_info tetInfo = elements.request();

    if (vtxInfo.ndim != (py::ssize_t)1 || vtxInfo.format != py::format_descriptor<float>::format()) {
      std::cerr << "Wrong vertex type:" << vtxInfo.ndim << ',' << vtxInfo.format << std::endl;
      return TetMesh();
    }

    if (tetInfo.ndim != (py::ssize_t)1 || (tetInfo.format != py::format_descriptor<int>::format())) {
      std::cerr << "Wrong tet type:" << tetInfo.ndim << ',' << tetInfo.format << ',' << py::format_descriptor<int64_t>::format() << std::endl;
      return TetMesh();
    }

    Eigen::VectorXd vertexPosDouble = Eigen::Map<const Eigen::VectorXf>((float *)vtxInfo.ptr, vtxInfo.size).cast<double>();
    Eigen::VectorXi tets = Eigen::Map<const Eigen::VectorXi>((int *)tetInfo.ptr, tetInfo.size);

    pgoTetMeshStructHandle tetmesh = pgo_create_tetmesh((int)vertexPosDouble.size() / 3, vertexPosDouble.data(), (int)tets.size() / 4, tets.data(), (double)E, (double)nu, (double)density);
    return TetMesh(tetmesh);
  });

  m.def("create_tetmesh_from_file", [](const std::string &tetmeshFilename) -> TetMesh {
    pgoTetMeshStructHandle tetmesh = pgo_create_tetmesh_from_file(const_cast<char *>(tetmeshFilename.c_str()));
    return TetMesh(tetmesh);
  });

  m.def("save_tetmesh_to_file", [](const TetMesh &tetmesh, const std::string &tetmeshFilename) {
    pgo_save_tetmesh_to_file(tetmesh.handle, const_cast<char *>(tetmeshFilename.c_str()));
  });

  m.def("get_tetmesh_vertex_positions", [](TetMesh &tetmesh) -> py::array_t<float> {
    int nVertices = pgo_tetmesh_get_num_vertices(tetmesh.handle);
    Eigen::VectorXd vertices(nVertices * 3);
    pgo_tetmesh_get_vertices(tetmesh.handle, vertices.data());
    auto ret = py::array_t<float>(nVertices * 3);
    py::buffer_info binfo = ret.request();
    (Eigen::Map<Eigen::VectorXf>((float *)binfo.ptr, nVertices * 3)) = vertices.cast<float>();
    ret.resize({ nVertices, 3 });
    return ret;
  });

  m.def("get_tetmesh_element_indices", [](TetMesh &tetmesh) -> py::array_t<int> {
    int nElements = pgo_tetmesh_get_num_tets(tetmesh.handle);
    Eigen::VectorXi elements(nElements * 4);
    pgo_tetmesh_get_elements(tetmesh.handle, elements.data());
    auto ret = py::array_t<int>(nElements * 4);
    py::buffer_info binfo = ret.request();
    (Eigen::Map<Eigen::VectorXi>((int *)binfo.ptr, nElements * 4)) = elements.cast<int>();
    ret.resize({ nElements, 4 });
    return ret;
  });

  m.def("update_tetmesh_vertices", [](TetMesh &tetmesh, pyArrayFloat vtxNew) -> TetMesh {
    py::buffer_info vtxInfo = vtxNew.request();
    if (vtxInfo.ndim != (py::ssize_t)2 || vtxInfo.shape[1] != 3 || vtxInfo.format != py::format_descriptor<float>::format()) {
      std::cerr << "Wrong vertex type:" << vtxInfo.ndim << ',' << vtxInfo.format << std::endl;
      return tetmesh;
    }

    vtxNew.resize({ vtxInfo.shape[0] * 3 });

    Eigen::VectorXd vtxNewDouble = Eigen::Map<const Eigen::VectorXf>((float *)vtxInfo.ptr, vtxInfo.shape[0] * 3).cast<double>();
    pgoTetMeshStructHandle tetmeshNewHandle = pgo_tetmesh_update_vertices(tetmesh.handle, vtxNewDouble.data());
    TetMesh tetmeshNew = TetMesh(tetmeshNewHandle);
    return tetmeshNew;
  });

  m.def(
    "create_element_laplacian_matrix", [](TetMeshGeo tetmesh, int face_neighbor, int n_repreat, int scale) -> SparseMatrix {
      pgoSparseMatrixStructHandle mat = pgo_create_tet_laplacian_matrix(tetmesh.handle, face_neighbor, n_repreat, scale);
      return SparseMatrix(mat);
    });

  m.def(
    "create_tet_biharmonic_gradient_matrix", [](TetMeshGeo tetmesh, int face_neighbor, int scale) -> SparseMatrix {
      pgoSparseMatrixStructHandle mat = pgo_create_tet_biharmonic_gradient_matrix(tetmesh.handle, face_neighbor, scale);
      return SparseMatrix(mat);
    });

  m.def(
    "conjugate_mv", [](SparseMatrix m, pyArrayFloat v) -> double {
      py::buffer_info binfo = v.request();
      Eigen::VectorXd vd = (Eigen::Map<Eigen::VectorXf>((float *)binfo.ptr, binfo.size)).cast<double>();

      double ret = pgo_conjugate_mv(m.handle, vd.data());
      return ret;
    });

  // m.def(
  //   "sp_mv", [](SparseMatrix m, pyArrayFloat v) -> py::array_t<float> {
  //     py::buffer_info binfo = v.request();
  //     Eigen::VectorXd vd = (Eigen::Map<Eigen::VectorXf>((float *)binfo.ptr, binfo.size)).cast<double>();

  //     auto ret = py::array_t<float>(binfo.size);
  //     py::buffer_info binfo_out = ret.request();

  //     Eigen::VectorXd voutd;
  //     pgo_sp_mv(m.handle, vd.data(), voutd.data());
  //     return ret;
  //   });

  m.def("sparse_matrix_get_num_entries", [](SparseMatrix m) -> int64_t {
    pgoSparseMatrixStructHandle mat = m.handle;
    return pgo_sparse_matrix_get_num_entries(mat);
  });

  m.def("sparse_matrix_get_row_indices", [](SparseMatrix m) -> py::array_t<int> {
    pgoSparseMatrixStructHandle mat = m.handle;
    int64_t nonZeros = pgo_sparse_matrix_get_num_entries(mat);
    auto ret = py::array_t<int>(nonZeros);
    py::buffer_info bInfo = ret.request();
    pgo_sparse_matrix_get_row_indices(mat, (int *)bInfo.ptr);
    return ret;
  });

  m.def("sparse_matrix_get_col_indices", [](SparseMatrix m) -> py::array_t<int> {
    pgoSparseMatrixStructHandle mat = m.handle;
    int64_t nonZeros = pgo_sparse_matrix_get_num_entries(mat);
    auto ret = py::array_t<int>(nonZeros);
    py::buffer_info bInfo = ret.request();
    pgo_sparse_matrix_get_col_indices(mat, (int *)bInfo.ptr);
    return ret;
  });

  m.def("sparse_matrix_get_values", [](SparseMatrix m) -> py::array_t<float> {
    pgoSparseMatrixStructHandle mat = m.handle;
    int64_t nonZeros = pgo_sparse_matrix_get_num_entries(mat);
    auto ret = py::array_t<float>(nonZeros);
    py::buffer_info bInfo = ret.request();
    Eigen::VectorXd values(nonZeros);
    pgo_sparse_matrix_get_values(mat, values.data());
    (Eigen::Map<Eigen::VectorXf>((float *)bInfo.ptr, nonZeros)) = values.cast<float>();
    return ret;
  });

  m.def("destroy_sparse_matrix", [](SparseMatrix m) {
    pgo_destroy_sparse_matrix(m.handle);
  });

  m.def("create_tet_gradient_matrix", [](TetMeshGeo m) -> SparseMatrix {
    return SparseMatrix(pgo_create_tet_gradient_matrix(m.handle));
  });

  m.def("create_tet_gradient_per_element_matrix", [](TetMeshGeo m) -> py::array_t<float> {
    int nTets = pgo_tetmeshgeo_get_num_tets(m.handle);
    Eigen::VectorXd v(nTets * 9 * 12);
    pgo_create_tet_gradient_per_element_matrix(m.handle, v.data());

    auto ret = py::array_t<float>({ nTets, 12, 9 });
    py::buffer_info bInfo = ret.request();
    (Eigen::Map<Eigen::VectorXf>((float *)bInfo.ptr, v.size())) = v.cast<float>();

    return ret;
  });

  m.def("create_smooth_rs_energy", [](TetMeshGeo tetmesh, double coeffR, double coeffS) -> SmoothRSEnergy {
    std::cout << "zz" << coeffR << ',' << coeffS << std::endl;
    pgoSmoothRSEnergyStructHandle h = pgo_create_smooth_rs_energy(tetmesh.handle, coeffR, coeffS);
    return SmoothRSEnergy(h);
  });

  m.def("destroy_smooth_rs_energy", [](SmoothRSEnergy energy) {
    pgo_destroy_smooth_rs_energy(energy.handle);
  });

  m.def("smooth_rs_energy_func", [](SmoothRSEnergy energy, py::array_t<float> x) -> double {
    py::buffer_info xInfo = x.request();

    if (xInfo.ndim != (py::ssize_t)1 || xInfo.format != py::format_descriptor<float>::format()) {
      std::cerr << "Wrong x type:" << xInfo.format << std::endl;
      return 1e20;
    }

    Eigen::VectorXd xDouble = Eigen::Map<const Eigen::VectorXf>((float *)xInfo.ptr, xInfo.size).cast<double>();
    return pgo_smooth_rs_energy_func(energy.handle, xDouble.data());
  });

  m.def("smooth_rs_energy_grad", [](SmoothRSEnergy energy, py::array_t<float> x) -> py::array_t<float> {
    py::buffer_info xInfo = x.request();
    auto grad = py::array_t<float>(xInfo.size);

    py::buffer_info gradInfo = grad.request();
    memset(gradInfo.ptr, 0, sizeof(float) * gradInfo.size);

    if (xInfo.ndim != (py::ssize_t)1 || xInfo.format != py::format_descriptor<float>::format()) {
      std::cerr << "Wrong x type:" << xInfo.format << std::endl;
      return grad;
    }

    Eigen::VectorXd xDouble = Eigen::Map<const Eigen::VectorXf>((float *)xInfo.ptr, xInfo.size).cast<double>();
    Eigen::VectorXd gradDouble(xInfo.size);
    gradDouble.setZero();
    pgo_smooth_rs_energy_grad(energy.handle, xDouble.data(), gradDouble.data());

    (Eigen::Map<Eigen::VectorXf>((float *)gradInfo.ptr, gradInfo.size)) = gradDouble.cast<float>();

    return grad;
  });

  m.def(
    "create_trimeshgeo", [](pyArrayFloat vertices, pyArrayInt triangles) -> TriMeshGeo {
      py::buffer_info vtxInfo = vertices.request();
      py::buffer_info triInfo = triangles.request();

      if (vtxInfo.ndim != (py::ssize_t)1 || vtxInfo.format != py::format_descriptor<float>::format()) {
        std::cerr << "Wrong vertex type:" << vtxInfo.ndim << ',' << vtxInfo.format << std::endl;
        return TriMeshGeo();
      }

#if defined(_WIN32)
      if (triInfo.ndim != (py::ssize_t)1) {
#else
        if (triInfo.ndim != (py::ssize_t)1 || (triInfo.format != py::format_descriptor<int>::format())) {
#endif
        std::cerr << "Wrong tri type:" << triInfo.ndim << ',' << triInfo.format << ',' << py::format_descriptor<int64_t>::format() << std::endl;
        return TriMeshGeo();
      }

      Eigen::VectorXd vertexPosDouble = Eigen::Map<const Eigen::VectorXf>((float *)vtxInfo.ptr, vtxInfo.size).cast<double>();
      Eigen::VectorXi tris = Eigen::Map<const Eigen::VectorXi>((int *)triInfo.ptr, triInfo.size);

      pgoTriMeshGeoStructHandle trimesh = pgo_create_trimeshgeo((int)vertexPosDouble.size() / 3, vertexPosDouble.data(), (int)tris.size() / 3, tris.data());
      return TriMeshGeo(trimesh);
    });

  m.def("trimeshgeo_get_num_vertices", [](TriMeshGeo m) -> int {
    return pgo_trimeshgeo_get_num_vertices(m.handle);
  });

  m.def("trimeshgeo_get_num_triangles", [](TriMeshGeo m) -> int {
    return pgo_trimeshgeo_get_num_triangles(m.handle);
  });

  m.def("trimeshgeo_get_vertices", [](TriMeshGeo m) -> py::array_t<float> {
    Eigen::VectorXd vtx(pgo_trimeshgeo_get_num_vertices(m.handle) * 3);
    pgo_trimeshgeo_get_vertices(m.handle, vtx.data());

    auto ret = py::array_t<float>(vtx.size());
    py::buffer_info binfo = ret.request();
    (Eigen::Map<Eigen::VectorXf>((float *)binfo.ptr, vtx.size())) = vtx.cast<float>();
    return ret;
  });

  m.def("trimeshgeo_get_triangles", [](TriMeshGeo m) -> py::array_t<int> {
    auto ret = py::array_t<int>(pgo_trimeshgeo_get_num_triangles(m.handle) * 3);
    py::buffer_info binfo = ret.request();

    pgo_trimeshgeo_get_triangles(m.handle, (int *)binfo.ptr);

    return ret;
  });

  m.def("destroy_trimeshgeo", [](TriMeshGeo m) {
    pgo_destroy_trimeshgeo(m.handle);
  });

  m.def("mesh_segmentation", [](TriMeshGeo m, int nSegs) -> py::array_t<int> {
    auto ret = py::array_t<int>(pgo_trimeshgeo_get_num_triangles(m.handle));
    py::buffer_info binfo = ret.request();

    pgo_mesh_segmentation(m.handle, nSegs, (int *)binfo.ptr);

    return ret;
  });

  m.def("mesh_isotropic_remeshing", [](TriMeshGeo m, double targetEdgeLength, int nIter, double angleThreshold) -> TriMeshGeo {
    return TriMeshGeo(pgo_mesh_isotropic_remeshing(m.handle, targetEdgeLength, nIter, angleThreshold));
  });

  m.def("trimesh_closest_distances", [](TriMeshGeo m, pyArrayFloat queryPt) -> py::array_t<float> {
    py::buffer_info quertPtInfo = queryPt.request();
    if (quertPtInfo.ndim != (py::ssize_t)1 || quertPtInfo.format != py::format_descriptor<float>::format()) {
      std::cerr << "Wrong vertex type:" << quertPtInfo.ndim << ',' << quertPtInfo.format << std::endl;
      return py::array_t<float>();
    }

    auto ret = py::array_t<float>(quertPtInfo.size / 3);

    Eigen::VectorXd queryPtDouble = Eigen::Map<const Eigen::VectorXf>((float *)quertPtInfo.ptr, quertPtInfo.size).cast<double>();
    Eigen::VectorXd queryDist(quertPtInfo.size / 3);
    pgo_trimesh_closest_distances(m.handle, (int)quertPtInfo.size / 3, queryPtDouble.data(), queryDist.data(), nullptr);

    py::buffer_info binfo = ret.request();
    (Eigen::Map<Eigen::VectorXf>((float *)binfo.ptr, binfo.size)) = queryDist.cast<float>();

    return ret;
  });

  m.def("tetmesh_barycentric_weights", [](TetMeshGeo m, pyArrayFloat queryPt) -> std::tuple<py::array_t<float>, py::array_t<int>> {
    py::buffer_info quertPtInfo = queryPt.request();
    if (quertPtInfo.ndim != (py::ssize_t)1 || quertPtInfo.format != py::format_descriptor<float>::format()) {
      std::cerr << "Wrong vertex type:" << quertPtInfo.ndim << ',' << quertPtInfo.format << std::endl;
      return std::tuple<py::array_t<float>, py::array_t<int>>();
    }

    auto ret_w = py::array_t<float>(quertPtInfo.size / 3 * 4);
    auto ret_i = py::array_t<int>(quertPtInfo.size / 3);

    py::buffer_info winfo = ret_w.request();
    py::buffer_info iinfo = ret_i.request();

    Eigen::VectorXd queryPtDouble = Eigen::Map<const Eigen::VectorXf>((float *)quertPtInfo.ptr, quertPtInfo.size).cast<double>();
    Eigen::VectorXd queryW(quertPtInfo.size / 3 * 4);

    pgo_tetmesh_barycentric_weights(m.handle, (int)quertPtInfo.size / 3, queryPtDouble.data(), queryW.data(), (int *)iinfo.ptr);

    (Eigen::Map<Eigen::VectorXf>((float *)winfo.ptr, winfo.size)) = queryW.cast<float>();

    return std::make_tuple(ret_w, ret_i);
  });

  m.def(
    "run_sim_from_config", [](const std::string &filename) -> int {
      return pgo_run_sim_from_config(filename.c_str());
    });

  m.def("debug", []() {
    std::cout << "Test!!!" << std::endl;
  });

#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}
