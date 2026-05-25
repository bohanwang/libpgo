#include <nanobind/nanobind.h>
#include <nanobind/eigen/dense.h>
#include <nanobind/eigen/sparse.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/shared_ptr.h>

#include "solveDiagnostics.h"
#include "solverResult.h"
#include "potentialEnergy.h"
#include "potentialEnergies.h"
#include "minimizeEnergy.h"
#include "NewtonSolver.h"

#include "ipc/core/surfaceIPCCore.h"
#include "embeddedSurfaceFloorPotentialEnergy.h"
#include "ipc/embeddedSurfaceIPCPotentialEnergy.h"

namespace nb = nanobind;
using namespace nb::literals;

using pgo::NonlinearOptimization::PotentialEnergy;
using pgo::NonlinearOptimization::MaxStepResult;
using pgo::NonlinearOptimization::SolveStatus;
using pgo::NonlinearOptimization::SolverResult;
using pgo::NonlinearOptimization::PotentialEnergies;
using pgo::NonlinearOptimization::NewtonSolver;

using CRefXd = pgo::EigenSupport::ConstRefVecXd;
using RefXd  = pgo::EigenSupport::RefVecXd;
using SMat = pgo::EigenSupport::SpMatD;

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
    NB_OVERRIDE_PURE(getNumDOFs, );
  }
  MaxStepResult computeMaxStepLimit(CRefXd x, CRefXd dx) const override {
    NB_OVERRIDE(computeMaxStepLimit, x, dx);
  }
};

NB_MODULE(pypgo, m) {
  m.doc() = "libpgo nanobind bindings";

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
