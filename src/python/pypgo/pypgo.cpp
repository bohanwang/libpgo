#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/eigen/dense.h>
#include <nanobind/eigen/sparse.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/trampoline.h>

#include "EigenSupport.h"
#include "generateTetMeshMatrix.h"
#include "geometryQuery.h"
#include "solveDiagnostics.h"
#include "solverResult.h"
#include "potentialEnergy.h"
#include "potentialEnergies.h"
#include "minimizeEnergy.h"
#include "NewtonSolver.h"
#include "pgoLogging.h"
#include "tetMesh.h"
#include "tetMeshGeo.h"

#include "ipc/core/surfaceIPCCore.h"
#include "embeddedSurfaceFloorPotentialEnergy.h"
#include "ipc/embeddedSurfaceIPCPotentialEnergy.h"

#if defined(PGO_HAS_ANIMATION_IO)
#  include "animationLoader.h"
#endif

namespace nb = nanobind;
using namespace nb::literals;

using pgo::NonlinearOptimization::PotentialEnergy;
using pgo::NonlinearOptimization::MaxStepResult;
using pgo::NonlinearOptimization::SolveStatus;
using pgo::NonlinearOptimization::SolverResult;
using pgo::NonlinearOptimization::PotentialEnergies;
using pgo::NonlinearOptimization::NewtonSolver;
namespace EnergyOptimizer = pgo::NonlinearOptimization::EnergyOptimizer;

using CRefXd = pgo::EigenSupport::ConstRefVecXd;
using RefXd  = pgo::EigenSupport::RefVecXd;
using SMat = pgo::EigenSupport::SpMatD;
using TetMeshGeoPtr = std::shared_ptr<pgo::Mesh::TetMeshGeo>;

struct SparseMatrix {
  SparseMatrix() = default;
  explicit SparseMatrix(pgo::EigenSupport::SpMatD m): mat(std::move(m)) {}
  pgo::EigenSupport::SpMatD mat;
};

using SparseMatrixPtr = std::shared_ptr<SparseMatrix>;

namespace {

void require_flat(const nb::ndarray<nb::ro, nb::c_contig, nb::device::cpu> &a, size_t multiple, const char *name) {
  if (a.ndim() != 1) {
    throw std::runtime_error(std::string(name) + " must be a flat 1D array");
  }
  if (a.size() % multiple != 0) {
    throw std::runtime_error(std::string(name) + " size is not divisible by " + std::to_string(multiple));
  }
}

std::vector<double> ndarray_to_double_vector(const nb::ndarray<nb::ro, nb::c_contig, nb::device::cpu> &a, const char *name) {
  require_flat(a, 1, name);
  std::vector<double> values(a.size());

  if (a.dtype() == nb::dtype<double>()) {
    const double *ptr = static_cast<const double *>(a.data());
    std::copy(ptr, ptr + a.size(), values.begin());
  }
  else if (a.dtype() == nb::dtype<float>()) {
    const float *ptr = static_cast<const float *>(a.data());
    std::transform(ptr, ptr + a.size(), values.begin(), [](float v) { return static_cast<double>(v); });
  }
  else {
    throw std::runtime_error(std::string(name) + " must have dtype float32 or float64");
  }

  return values;
}

std::vector<int> ndarray_to_int_vector(const nb::ndarray<nb::ro, nb::c_contig, nb::device::cpu> &a, const char *name) {
  require_flat(a, 1, name);
  std::vector<int> values(a.size());

  if (a.dtype() == nb::dtype<int>()) {
    const int *ptr = static_cast<const int *>(a.data());
    std::copy(ptr, ptr + a.size(), values.begin());
  }
  else if (a.dtype() == nb::dtype<int64_t>()) {
    const int64_t *ptr = static_cast<const int64_t *>(a.data());
    std::transform(ptr, ptr + a.size(), values.begin(), [](int64_t v) { return static_cast<int>(v); });
  }
  else {
    throw std::runtime_error(std::string(name) + " must have dtype int32 or int64");
  }

  return values;
}

Eigen::VectorXi sparse_rows(const pgo::EigenSupport::SpMatD &mat) {
  Eigen::VectorXi rows(static_cast<Eigen::Index>(mat.nonZeros()));
  Eigen::Index inc = 0;
  for (int outer = 0; outer < mat.outerSize(); ++outer) {
    for (pgo::EigenSupport::SpMatD::InnerIterator it(mat, outer); it; ++it) {
      rows[inc++] = static_cast<int>(it.row());
    }
  }
  return rows;
}

Eigen::VectorXi sparse_cols(const pgo::EigenSupport::SpMatD &mat) {
  Eigen::VectorXi cols(static_cast<Eigen::Index>(mat.nonZeros()));
  Eigen::Index inc = 0;
  for (int outer = 0; outer < mat.outerSize(); ++outer) {
    for (pgo::EigenSupport::SpMatD::InnerIterator it(mat, outer); it; ++it) {
      cols[inc++] = static_cast<int>(it.col());
    }
  }
  return cols;
}

Eigen::VectorXd sparse_values(const pgo::EigenSupport::SpMatD &mat) {
  Eigen::VectorXd values(static_cast<Eigen::Index>(mat.nonZeros()));
  Eigen::Index inc = 0;
  for (int outer = 0; outer < mat.outerSize(); ++outer) {
    for (pgo::EigenSupport::SpMatD::InnerIterator it(mat, outer); it; ++it) {
      values[inc++] = it.value();
    }
  }
  return values;
}

}  // namespace

// ── Trampoline for PotentialEnergy ──
struct PyPotentialEnergy : PotentialEnergy {
  NB_TRAMPOLINE(PotentialEnergy, 7);

  double func(CRefXd x) const override {
    NB_OVERRIDE_PURE(func, x);
  }
  void gradient(CRefXd x, RefXd grad) const override {
    NB_OVERRIDE_PURE(gradient, x, grad);
  }
  void hessian(CRefXd x, SMat &hess) const override {
    NB_OVERRIDE_PURE(hessian, x, hess);
  }
  void createHessian(SMat &hess) const override {
    NB_OVERRIDE_PURE(createHessian, hess);
  }
  void getDOFs(std::vector<int> &dofs) const override {
    NB_OVERRIDE_PURE(getDOFs, dofs);
  }
  int getNumDOFs() const override {
    NB_OVERRIDE_PURE(getNumDOFs);
  }
  MaxStepResult computeMaxStepLimit(CRefXd x, CRefXd dx) const override {
    NB_OVERRIDE(computeMaxStepLimit, x, dx);
  }
};

NB_MODULE(pypgo, m) {
  m.doc() = "libpgo nanobind bindings";

  nb::class_<pgo::Mesh::TetMeshGeo>(m, "TetMeshGeo");
  nb::class_<SparseMatrix>(m, "SparseMatrix");

  m.def("init", []() {
    pgo::Logging::init();
  });

  m.def("convert_animation_to_abc", [](const std::string &configFileName, const std::string &outputFolder) {
#if defined(PGO_HAS_ANIMATION_IO)
    pgo::AnimationIO::AnimationLoader loader;
    if (loader.load(configFileName.c_str()) != 0) {
      return 1;
    }
    return loader.saveABC(outputFolder.c_str());
#else
    (void)configFileName;
    (void)outputFolder;
    return 1;
#endif
  }, "config_file_name"_a, "output_folder"_a);

  m.def("create_tetmeshgeo",
    [](const nb::ndarray<nb::ro, nb::c_contig, nb::device::cpu> &vertices,
       const nb::ndarray<nb::ro, nb::c_contig, nb::device::cpu> &tets) {
      require_flat(vertices, 3, "vertices");
      require_flat(tets, 4, "tets");
      std::vector<double> vertexData = ndarray_to_double_vector(vertices, "vertices");
      std::vector<int> tetData = ndarray_to_int_vector(tets, "tets");
      return std::make_shared<pgo::Mesh::TetMeshGeo>(
        static_cast<int>(vertexData.size() / 3), vertexData.data(),
        static_cast<int>(tetData.size() / 4), tetData.data());
    },
    "vertices"_a, "tets"_a);

  m.def("create_tetmeshgeo_from_file", [](const std::string &filename) -> TetMeshGeoPtr {
    pgo::VolumetricMeshes::TetMesh mesh(filename.c_str());
    auto tetmesh = std::make_shared<pgo::Mesh::TetMeshGeo>();
    mesh.exportMeshGeometry(*tetmesh);
    return tetmesh;
  }, "filename"_a);

  m.def("tetmeshgeo_get_num_vertices", [](const TetMeshGeoPtr &tetmesh) {
    return tetmesh->numVertices();
  }, "tetmesh"_a);

  m.def("tetmeshgeo_get_num_tets", [](const TetMeshGeoPtr &tetmesh) {
    return tetmesh->numTets();
  }, "tetmesh"_a);

  m.def("tetmeshgeo_get_vertices", [](const TetMeshGeoPtr &tetmesh) {
    Eigen::VectorXd vertices(tetmesh->numVertices() * 3);
    for (int vi = 0; vi < tetmesh->numVertices(); ++vi) {
      vertices.segment<3>(vi * 3) = tetmesh->pos(vi);
    }
    return vertices;
  }, "tetmesh"_a);

  m.def("tetmeshgeo_get_tets", [](const TetMeshGeoPtr &tetmesh) {
    Eigen::VectorXi tets(tetmesh->numTets() * 4);
    for (int ti = 0; ti < tetmesh->numTets(); ++ti) {
      tets[ti * 4] = tetmesh->tet(ti)[0];
      tets[ti * 4 + 1] = tetmesh->tet(ti)[1];
      tets[ti * 4 + 2] = tetmesh->tet(ti)[2];
      tets[ti * 4 + 3] = tetmesh->tet(ti)[3];
    }
    return tets;
  }, "tetmesh"_a);

  m.def("destroy_tetmeshgeo", [](TetMeshGeoPtr &) {}, "tetmesh"_a);

  m.def("sparse_matrix_get_num_entries", [](const SparseMatrixPtr &mat) {
    return static_cast<int64_t>(mat->mat.nonZeros());
  }, "matrix"_a);

  m.def("sparse_matrix_get_row_indices", [](const SparseMatrixPtr &mat) {
    return sparse_rows(mat->mat);
  }, "matrix"_a);

  m.def("sparse_matrix_get_col_indices", [](const SparseMatrixPtr &mat) {
    return sparse_cols(mat->mat);
  }, "matrix"_a);

  m.def("sparse_matrix_get_values", [](const SparseMatrixPtr &mat) {
    return sparse_values(mat->mat);
  }, "matrix"_a);

  m.def("destroy_sparse_matrix", [](SparseMatrixPtr &) {}, "matrix"_a);

  m.def("create_tet_laplacian_matrix",
    [](const TetMeshGeoPtr &tetmesh, int faceNeighbor, int repeat, int scale) {
      namespace ES = pgo::EigenSupport;
      ES::SpMatD L;
      pgo::SolidDeformationModel::TetMeshMatrix::generateBasicElementLaplacianMatrix(*tetmesh, L, faceNeighbor, scale);
      if (repeat > 1) {
        ES::SpMatD Ln;
        ES::expandN(L, Ln, repeat);
        L = std::move(Ln);
      }
      return std::make_shared<SparseMatrix>(std::move(L));
    },
    "tetmesh"_a, "face_neighbor"_a, "repeat"_a, "scale"_a);

  m.def("create_element_laplacian_matrix",
    [](const TetMeshGeoPtr &tetmesh, int faceNeighbor, int repeat, int scale) {
      namespace ES = pgo::EigenSupport;
      ES::SpMatD L;
      pgo::SolidDeformationModel::TetMeshMatrix::generateBasicElementLaplacianMatrix(*tetmesh, L, faceNeighbor, scale);
      if (repeat > 1) {
        ES::SpMatD Ln;
        ES::expandN(L, Ln, repeat);
        L = std::move(Ln);
      }
      return std::make_shared<SparseMatrix>(std::move(L));
    },
    "tetmesh"_a, "face_neighbor"_a, "repeat"_a, "scale"_a);

  m.def("create_tet_gradient_matrix", [](const TetMeshGeoPtr &tetmesh) {
    namespace ES = pgo::EigenSupport;
    ES::SpMatD G;
    pgo::SolidDeformationModel::TetMeshMatrix::generateGradientMatrix(*tetmesh, G);
    return std::make_shared<SparseMatrix>(std::move(G));
  }, "tetmesh"_a);

  m.def("create_tet_gradient_per_element_matrix", [](const TetMeshGeoPtr &tetmesh) {
    namespace ES = pgo::EigenSupport;
    auto values = new std::vector<double>(static_cast<size_t>(tetmesh->numTets()) * 9 * 12);
    for (int i = 0; i < tetmesh->numTets(); ++i) {
      ES::M9x12d G;
      pgo::SolidDeformationModel::TetMeshMatrix::generateElementGradientMatrix(*tetmesh, i, G);
      ES::Mp<ES::M9x12d>(values->data() + static_cast<size_t>(i) * 9 * 12) = G;
    }
    nb::capsule owner(values, [](void *p) noexcept {
      delete static_cast<std::vector<double> *>(p);
    });
    return nb::ndarray<nb::numpy, double>(
      values->data(),
      { static_cast<size_t>(tetmesh->numTets()), static_cast<size_t>(12), static_cast<size_t>(9) },
      owner);
  }, "tetmesh"_a);

  m.def("create_tet_biharmonic_gradient_matrix",
    [](const TetMeshGeoPtr &tetmesh, int faceNeighbor, int scale) {
      namespace ES = pgo::EigenSupport;
      ES::SpMatD L;
      pgo::SolidDeformationModel::TetMeshMatrix::generateBasicElementLaplacianMatrix(*tetmesh, L, faceNeighbor, 0);

      ES::SpMatD Ln;
      ES::expandN(L, Ln, 9);

      ES::SpMatD G;
      pgo::SolidDeformationModel::TetMeshMatrix::generateGradientMatrix(*tetmesh, G);

      ES::SpMatD LG;
      ES::mm(Ln, G, LG);

      ES::SpMatD GTLTLG;
      if (scale) {
        std::vector<double> eleVols(tetmesh->numTets());
        for (int ei = 0; ei < tetmesh->numTets(); ++ei) {
          eleVols[ei] = std::abs(pgo::Mesh::getTetDeterminant(
            tetmesh->pos(ei, 0), tetmesh->pos(ei, 1), tetmesh->pos(ei, 2), tetmesh->pos(ei, 3)));
        }

        double maxVol = *std::max_element(eleVols.begin(), eleVols.end());
        for (int ei = 0; ei < tetmesh->numTets(); ++ei) {
          eleVols[ei] /= maxVol;
        }

        ES::SpMatD M;
        std::vector<ES::TripletD> entries;
        for (int i = 0; i < tetmesh->numTets(); ++i) {
          for (int j = 0; j < 9; ++j) {
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

      return std::make_shared<SparseMatrix>(std::move(GTLTLG));
    },
    "tetmesh"_a, "face_neighbor"_a, "scale"_a);

  // ── MaxStepResult ──
  nb::class_<MaxStepResult>(m, "MaxStepResult")
    .def(nb::init<>())
    .def_rw("alpha", &MaxStepResult::alpha)
    .def_rw("material_alpha", &MaxStepResult::materialAlpha)
    .def_rw("contact_alpha", &MaxStepResult::contactAlpha)
    .def_rw("material_clamped", &MaxStepResult::materialClamped)
    .def_rw("contact_clamped", &MaxStepResult::contactClamped)
    .def_static("unconstrained", &MaxStepResult::unconstrained);

  // ── SolveStatus ──
  nb::enum_<SolveStatus>(m, "SolveStatus")
    .value("Converged", SolveStatus::Converged)
    .value("MaxIterations", SolveStatus::MaxIterations)
    .value("LineSearchFailed", SolveStatus::LineSearchFailed)
    .value("StepTooSmall", SolveStatus::StepTooSmall)
    .value("NonFinite", SolveStatus::NonFinite);

  // ── SolverResult ──
  nb::class_<SolverResult>(m, "SolverResult")
    .def(nb::init<>())
    .def_ro("status", &SolverResult::status)
    .def_ro("iterations", &SolverResult::iterations);

  // ── NewtonSolver params ──
  nb::enum_<NewtonSolver::LineSearchMethod>(m, "LineSearchMethod")
    .value("GOLDEN", NewtonSolver::LSM_GOLDEN)
    .value("BRENTS", NewtonSolver::LSM_BRENTS)
    .value("BACKTRACK", NewtonSolver::LSM_BACKTRACK)
    .value("SIMPLE", NewtonSolver::LSM_SIMPLE);

  nb::class_<NewtonSolver::SolverParam>(m, "NewtonSolverParams")
    .def(nb::init<>())
    .def_rw("alpha", &NewtonSolver::SolverParam::alpha)
    .def_rw("line_search_method", &NewtonSolver::SolverParam::lsm)
    .def_rw("stop_after_increase", &NewtonSolver::SolverParam::stopAfterIncrease)
    .def_rw("add_damping", &NewtonSolver::SolverParam::addDamping);

  // ── PotentialEnergy (base + trampoline) ──
  nb::class_<PotentialEnergy, PyPotentialEnergy>(m, "PotentialEnergy")
    .def(nb::init<>())
    .def("func", &PotentialEnergy::func)
    .def("gradient", &PotentialEnergy::gradient)
    .def("hessian", &PotentialEnergy::hessian)
    .def("create_hessian", &PotentialEnergy::createHessian)
    .def("get_num_dofs", &PotentialEnergy::getNumDOFs)
    .def("compute_max_step_limit", &PotentialEnergy::computeMaxStepLimit);

  // ── PotentialEnergies ──
  nb::class_<PotentialEnergies, PotentialEnergy>(m, "PotentialEnergies")
    .def(nb::init<int>())
    .def("add", &PotentialEnergies::addPotentialEnergy, "energy"_a, "coeff"_a = 1.0)
    .def("init", &PotentialEnergies::init);

  // ── SolverType ──
  nb::enum_<pgo::NonlinearOptimization::EnergyOptimizer::SolverType>(m, "SolverType")
    .value("IPOPT", pgo::NonlinearOptimization::EnergyOptimizer::SolverType::ST_IPOPT)
    .value("KNITRO", pgo::NonlinearOptimization::EnergyOptimizer::SolverType::ST_KNITRO)
    .value("NEWTON", pgo::NonlinearOptimization::EnergyOptimizer::SolverType::ST_NEWTON);

  // ── minimize ──
  m.def("minimize", [](Eigen::Ref<Eigen::VectorXd> x,
                        std::shared_ptr<const PotentialEnergy> energy,
                        pgo::NonlinearOptimization::EnergyOptimizer::SolverType solver,
                        int maxIter, double eps, int verbose) {
    Eigen::VectorXd xlow = Eigen::VectorXd::Constant(x.size(), -1e100);
    Eigen::VectorXd xhi  = Eigen::VectorXd::Constant(x.size(),  1e100);
    return EnergyOptimizer::minimize(x, energy, xlow, xhi, solver, maxIter, eps, verbose);
  }, "x"_a, "energy"_a, "solver"_a = pgo::NonlinearOptimization::EnergyOptimizer::SolverType::ST_NEWTON,
     "max_iter"_a = 100, "eps"_a = 1e-6, "verbose"_a = 0);

  // ── NewtonSolver ──
  nb::class_<NewtonSolver>(m, "NewtonSolver")
    .def(nb::init<const double*, NewtonSolver::SolverParam,
                   std::shared_ptr<const PotentialEnergy>,
                   const std::vector<int>&, const double*>(),
         "x0"_a, "params"_a, "energy"_a, "fixed_dofs"_a, "fixed_values"_a = nb::none())
    .def("solve", [](NewtonSolver &self, Eigen::Ref<Eigen::VectorXd> x,
                      int numIter, double eps, int verbose) {
      return self.solve(x.data(), numIter, eps, verbose);
    }, "x"_a, "num_iter"_a = 100, "eps"_a = 1e-6, "verbose"_a = 0);

  // ── IPC parameters ──
  nb::class_<pgo::Contact::IPC::SurfaceIPCCore::Parameters>(m, "IPCCoreParams")
    .def(nb::init<>())
    .def_rw("dhat", &pgo::Contact::IPC::SurfaceIPCCore::Parameters::dhat)
    .def_rw("kappa", &pgo::Contact::IPC::SurfaceIPCCore::Parameters::kappa)
    .def_rw("eps_ee", &pgo::Contact::IPC::SurfaceIPCCore::Parameters::eps_ee)
    .def_rw("slackness", &pgo::Contact::IPC::SurfaceIPCCore::Parameters::slackness);

  nb::enum_<pgo::Contact::IPC::FloorAxis>(m, "FloorAxis")
    .value("INVALID", pgo::Contact::IPC::FloorAxis::INVALID)
    .value("X", pgo::Contact::IPC::FloorAxis::X)
    .value("Y", pgo::Contact::IPC::FloorAxis::Y)
    .value("Z", pgo::Contact::IPC::FloorAxis::Z);

  nb::enum_<pgo::Contact::IPC::FloorSide>(m, "FloorSide")
    .value("KEEP_ABOVE", pgo::Contact::IPC::FloorSide::KEEP_ABOVE)
    .value("KEEP_BELOW", pgo::Contact::IPC::FloorSide::KEEP_BELOW);

  nb::class_<pgo::Contact::IPC::FloorPenaltyParameters>(m, "FloorPenaltyParams")
    .def(nb::init<>())
    .def_rw("axis", &pgo::Contact::IPC::FloorPenaltyParameters::floorAxis)
    .def_rw("side", &pgo::Contact::IPC::FloorPenaltyParameters::floorSide)
    .def_rw("height", &pgo::Contact::IPC::FloorPenaltyParameters::floorHeight)
    .def_rw("kappa", &pgo::Contact::IPC::FloorPenaltyParameters::floorKappa);

  // ── Contact energies ──
  nb::class_<pgo::Contact::IPC::EmbeddedSurfaceIPCPotentialEnergy, PotentialEnergy>(m, "IPCEnergy")
    .def(nb::init<const Eigen::MatrixXd&, const Eigen::MatrixXi&,
                   const Eigen::SparseMatrix<double>&,
                   const pgo::Contact::IPC::SurfaceIPCCore::Parameters&>(),
         "surface_rest"_a, "surface_triangles"_a, "embedding"_a,
         "ipc_params"_a = pgo::Contact::IPC::SurfaceIPCCore::Parameters{})
    .def("mark_obstacle_static",
         &pgo::Contact::IPC::EmbeddedSurfaceIPCPotentialEnergy::markObstacleStatic);

  nb::class_<pgo::Contact::IPC::EmbeddedSurfaceFloorPotentialEnergy, PotentialEnergy>(m, "FloorEnergy")
    .def(nb::init<const Eigen::MatrixXd&, const Eigen::SparseMatrix<double>&,
                   const pgo::Contact::IPC::FloorPenaltyParameters&>(),
         "surface_rest"_a, "embedding"_a, "params"_a);
}
