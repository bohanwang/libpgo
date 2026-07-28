# libpgo 0.0.4 Implementation Plan

The short operational gate list and evidence record live in
[`RELEASE_0.0.4_CHECKLIST.md`](RELEASE_0.0.4_CHECKLIST.md). This document contains the detailed implementation contract.

| Field | Value |
| --- | --- |
| Status | Scope frozen; implementation not started |
| Target version | `0.0.4` |
| Target repository | `/Users/jinceyang/Desktop/codebase/merge/libpgo` |
| Target baseline branch | `upstream-temp-0.0.4` |
| Target baseline commit | `1670f6258d941d80cd9d718fc033fccc193d5e20` |

This is an internal release document. The absolute paths are intentional: they identify the three local code snapshots used for the 0.0.4 audit and prevent “current branch”, “temp”, and “main” from being confused during implementation.

## 1. Repository roles and immutable references

| Role | Absolute path | Branch / revision | How it may be used |
| --- | --- | --- | --- |
| 0.0.4 implementation target | `/Users/jinceyang/Desktop/codebase/merge/libpgo` | `upstream-temp-0.0.4` at `1670f6258d941d80cd9d718fc033fccc193d5e20` | All 0.0.4 implementation work lands here on `release/0.0.4`, created from this revision. |
| 0.0.4 integration target | `/Users/jinceyang/Desktop/codebase/merge/libpgo` | `upstream-0.0.4` | Merge completed release work here, then run every final gate against the exact integration commit before owner review. |
| Current experimental repository | `/Users/jinceyang/Desktop/codebase/libpgo` | `feat/generalized-neohookean` at `3bc08226c10e3a7dde862693904b49718a766dfa` | Reference implementation only. Selectively migrate build/CI/docs ideas and isolated fixes; do not merge the branch wholesale. |
| 0.0.3/main reference | `/Users/jinceyang/Desktop/codebase/main/libpgo` | detached worktree of `main` at `ed9675d56403b73a0cf32084e5da9b0011f88e2c` | Compatibility baseline and canonical example-asset source. Treat as read-only. |
| 0.0.4 docs target | `/Users/jinceyang/Desktop/codebase/merge/libpgo/docs` git submodule | `annajcy/libpgo-doc` branch `v0.0.4` at `ef0834c67a280d16836579942246129087c99e44` | Initialize/check out the canonical docs repository inside this release repository, preserve site infrastructure, and rewrite content to match 0.0.4. |
| Standalone noncanonical docs checkout | `/Users/jinceyang/Desktop/codebase/libpgo-doc` | different repository remote | Do not use as the 0.0.4 docs target unless its remote is explicitly changed and verified. |

Repository rules:

- [ ] Create the code release branch from exactly `/Users/jinceyang/Desktop/codebase/merge/libpgo@1670f6258d941d80cd9d718fc033fccc193d5e20`.
- [ ] Record the final code commit, docs commit, dependency revisions, toolchain versions, and artifact checksums in the release record.
- [ ] Merge completed work into `upstream-0.0.4` and bind all final test evidence to its exact integration commit.
- [ ] Stop for owner review before creating or merging an `upstream-0.0.4` to `main` pull request.
- [ ] Never merge `/Users/jinceyang/Desktop/codebase/libpgo` wholesale into `/Users/jinceyang/Desktop/codebase/merge/libpgo`.
- [ ] Never copy all CMake, Python package, benchmark, experiment, or asset changes from `/Users/jinceyang/Desktop/codebase/libpgo`.
- [ ] Use `/Users/jinceyang/Desktop/codebase/main/libpgo` as the 0.0.3 compatibility and example-asset baseline.
- [ ] Derive the documentation changelog from the code difference between the 0.0.3 compatibility baseline and final 0.0.4 code; the empty docs `main` branch is not a content reference.
- [ ] Preserve unrelated user changes if any worktree becomes dirty.

### 1.1 How an implementation agent must execute this document

This file is an implementation brief and checklist, not a shell script.

- [ ] Start commands from `/Users/jinceyang/Desktop/codebase/merge/libpgo`.
- [ ] Treat `/Users/jinceyang/Desktop/codebase/libpgo` and `/Users/jinceyang/Desktop/codebase/main/libpgo` as read-only references.
- [ ] Implement in the phase and commit order below; do not attempt all changes as one unreviewable patch.
- [ ] After each phase, update this checklist with the focused tests actually run and their results.
- [ ] Do not mark a server-only item complete without a Linux server result tied to an exact commit.
- [ ] Treat the Linux server matrix as an owner-run manual release gate, not a GitHub Actions workflow or required CI status check.
- [ ] Do not infer authorization to push, merge, tag, run final artifact-delivery workflows, or deploy docs merely because those actions appear in the final release checklist.
- [ ] Default stopping point without additional release authorization: the work is merged into `upstream-0.0.4`, all final gates pass on its exact integration commit, and the code, docs, artifacts, and Linux server evidence are ready for owner review.
- [ ] The release remains blocked until the owner supplies or confirms the Linux server matrix report.
- [ ] Final tag, artifact-delivery, and docs-deployment actions in Section 11 require a separate explicit go-ahead.

## 2. Frozen 0.0.4 release contract

### 2.1 Supported release features

- [ ] Retain the cubic volumetric mesh implementation.
- [ ] Retain the cubic FEM/deformation model.
- [ ] Retain the cubic mesher.
- [ ] Retain sampled-penalty contact.
- [ ] Retain the new IPC core, IPC energy, CCD maximum-step calculation, self-contact, and floor-contact path.
- [ ] Guarantee the public config-runner matrix in both dynamic and static modes:

  | Simulation domain | Contact model | Required |
  | --- | --- | --- |
  | tet | sampled | Yes |
  | tet | IPC | Yes |
  | cubic | sampled | Yes |
  | cubic | IPC | Yes |
  | shell | IPC | Yes |

- [ ] Make shell + IPC the only required shell case and migrate its configs/assets to the new separated example layout; sampled shell is not implemented or tested as part of the 0.0.4 release contract.
- [ ] Do not add IPC external-object contact if it is not already implemented. The current `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/tools/runSim/runIPCSimSetup.cpp` rejects `external-objects`; document IPC 0.0.4 as supporting self-contact plus the configured floor, not arbitrary external objects.
- [ ] Preserve source compatibility with the 0.0.3 APIs wherever the agreed fixes do not require a change.
- [ ] Do not preserve 0.0.3 example file paths. The 0.0.4 examples switch directly to the new config/asset-separated layout without compatibility JSON copies.

### 2.2 Public API contract

- [ ] Use one config-driven C/Python entry point for the five required combinations: tet/cubic × sampled/IPC plus shell + IPC.
- [ ] Replace the mixed CLI with two explicit executables:

  ```text
  runDynamicSim <config> [--resume]
  runStaticSim <config>
  ```

- [ ] `runDynamicSim` rejects static configs; `runStaticSim` rejects dynamic configs and does not accept `--resume`.
- [ ] Each executable supports all five required domain/contact cases; simulation mode is orthogonal to domain and contact selection.
- [ ] Both executables share the same config parser, domain classification, contact-model selection, and simulation implementation used by C and Python.
- [ ] Keep the C entry point:

  ```c
  int pgo_run_sim_from_config(const char *configFileName);
  ```

- [ ] Keep the Python entry point:

  ```python
  pypgo.run_sim_from_config(config_filename)
  ```

- [ ] Dispatch using config data:

  - `contact-model: "sampled"` selects the sampled runner.
  - `contact-model: "ipc"` selects the IPC runner.
  - A missing `contact-model` selects sampled contact for 0.0.3 compatibility.
  - Exactly one supported domain is selected: `tet-mesh`, `cubic-mesh`, or the existing shell schema.

- [ ] Keep the current 0.0.4 `pypgo` extension target.
- [ ] Do not migrate the `pypgo_core` target, abi3 packaging, nanobind conversion, or modular Python package architecture from `/Users/jinceyang/Desktop/codebase/libpgo`.
- [ ] Direct object-level Python wrappers for every new cubic/IPC class are not required for 0.0.4; the unified config runner is the required public surface.

### 2.3 Explicitly out of scope

- [ ] No new simulation algorithms beyond what exists in the 0.0.4 baseline.
- [ ] Wiring the existing IPC energy, derivatives, CCD/max-step, and solver components into `runStaticSim` is in scope; do not redesign IPC or introduce a different static-contact algorithm.
- [ ] No broad solver rewrite.
- [ ] No contact-system redesign.
- [ ] No mesh-system redesign.
- [ ] No package-architecture rewrite.
- [ ] No abi3 migration.
- [ ] No full CMake modularization copied from `/Users/jinceyang/Desktop/codebase/libpgo`.
- [ ] No generalized Neo-Hookean or other material work from the current experimental branch. If a minimal correctness dependency is discovered, stop and report it before expanding scope.
- [ ] No tricubic Hermite feature migration.
- [ ] No ARPACK restoration, option, dependency work, or CI job.
- [ ] No benchmark or research-experiment migration.
- [ ] No repository-wide formatting.
- [ ] No generated media expansion.
- [ ] No new permanent profiling framework.

### 2.4 Release parameters that do not block implementation

The implementation can proceed with the defaults below, but the owner must confirm the final values before the release-candidate gates:

- Wheel matrix: CPython 3.12 with `linux_x86_64`, `macosx` arm64, and `win_amd64` platform tags.
- Linux and Windows wheels use MKL via `PGO_USE_MKL=ON`, with `MKL_THREADING=tbb_thread`, `PGO_HAS_ORIG_PARDISO=OFF`, and runtime linkage verified to use TBB rather than Intel OpenMP.
- macOS wheels use `PGO_USE_MKL=OFF`; libpgo/SuiteSparse selects the system Accelerate framework with `BLA_VENDOR=Apple`.
- GMP, MPFR, Imath, fmt, TBB, and other non-system build/runtime dependencies come from the pinned Conda environment on all three platforms; do not install GMP/MPFR/Imath/fmt with `apt` or Homebrew.
- Wheels depend on the documented Conda runtime and do not vendor MKL, TBB, GMP, MPFR, Imath, fmt, or their transitive Conda runtime libraries.
- 0.0.4 is wheel-only: do not build or publish an sdist and do not upload to TestPyPI or PyPI.
- Linux server backend: MKL Pardiso via `PGO_USE_MKL=ON` and `PGO_HAS_ORIG_PARDISO=OFF`.
- Linux server duration/timeout: parameterized in the server driver and selected by the owner when the matrix is run.
- Linux server MKL thread count: explicitly supplied or accepted from a documented driver default; always recorded in the report.
- Cubic `.veg` policy: check in a deterministic Python generation script, not generated cubic `.veg` files. Generate meshes into a gitignored or explicitly supplied output directory for tests and examples.
- Delivery: the three platform CI workflows upload wheels and evidence as GitHub Actions artifacts; tagging, artifact retention/public linking, and docs deployment require separate explicit approval.

## 3. Release gates

A release candidate may be tagged only when every blocking gate is satisfied.

### Gate A — Correctness

- [ ] Native unit and integration tests pass on Linux, macOS, and Windows, excluding the explicitly server-only long matrix.
- [ ] Dynamic-Hessian tests pass.
- [ ] Fixed-topology solver symbolic-factorization lifecycle regression test passes.
- [ ] Tet and cubic maximum-step tests pass.
- [ ] IPC parameter, topology, state-vector, and spatial-hash validation tests pass.
- [ ] Cubic element rest-state validation passes for every cubic `.veg` generated by the release script during the release test workflow.
- [ ] Dynamic-output-directory and static-output-file safety tests prove that existing paths are rejected without deletion or modification.

### Gate B — API matrix

- [ ] Fast config/parser/API smoke passes in normal CI without running the long five-case dynamic simulation matrix.
- [ ] On the Linux server, C API passes tet + sampled smoke.
- [ ] On the Linux server, C API passes tet + IPC smoke.
- [ ] On the Linux server, C API passes cubic + sampled smoke.
- [ ] On the Linux server, C API passes cubic + IPC smoke.
- [ ] On the Linux server, C API passes shell + IPC smoke.
- [ ] On the Linux server, Python API passes the same five short API cases.
- [ ] On the Linux server, the representative longer tet/cubic sampled/IPC plus shell IPC matrix passes through `runDynamicSim`.
- [ ] The Linux server build demonstrably uses MKL Pardiso, not original/native Pardiso and not the Eigen fallback.
- [ ] Linux server and Linux/Windows wheel evidence demonstrate the MKL TBB threading layer; macOS evidence demonstrates no MKL linkage.
- [ ] Missing `contact-model` is tested and selects sampled contact.
- [ ] Invalid mesh/contact combinations fail with a useful error and a nonzero result.
- [ ] `runStaticSim` passes all five required short cases and rejects `--resume`.

### Gate C — Packaging

- [ ] Linux `linux_x86_64` wheel builds, is audited without repair/vendoring, and installs in a fresh documented Conda GMP/MPFR/Imath/fmt/MKL/TBB runtime environment.
- [ ] macOS arm64 wheel builds, is audited without vendoring Conda runtimes, and installs in a fresh documented Conda GMP/MPFR/Imath/fmt/TBB/non-MKL runtime environment.
- [ ] Windows `win_amd64` wheel builds, is audited without vendoring Conda runtimes, and installs in a fresh documented Conda GMP/MPFR/Imath/fmt/MKL/TBB runtime environment.
- [ ] Installed artifacts report `pypgo.__version__ == "0.0.4"`.
- [ ] `python -m pip check` passes.
- [ ] `python -m twine check` passes for all release artifacts.
- [ ] Wheels contain no test output, example output, cache, old wheel, build tree, or unapproved runtime library.
- [ ] Every platform workflow uploads exactly one wheel plus source-commit, dependency, linkage, test, checksum, and environment evidence as GitHub Actions artifacts.

### Gate D — Assets and examples

- [ ] Canonical assets inherited from `/Users/jinceyang/Desktop/codebase/main/libpgo/examples` and the retained shell inputs are migrated into one deduplicated `examples/assets/` tree.
- [ ] All configs are migrated into `examples/configs/`; no old-path compatibility configs are retained.
- [ ] Extra duplicate volume objects, generated media, and superseded scenes are removed without deleting either retained shell example tree.
- [ ] Generated cubic meshes pass deterministic topology, rest-state, scale, and smoke checks.
- [ ] New example configs resolve canonical `examples/assets/` and generated cubic paths from the config file’s directory after the generation step.
- [ ] Short smoke variants run in CI without writing into the source tree.

### Gate E — Documentation and release metadata

- [ ] Main repository release notes are complete.
- [ ] Detailed module changelogs are complete and build successfully from `annajcy/libpgo-doc` branch `v0.0.4`; deployment happens only after explicit authorization.
- [ ] The docs submodule is pinned to an exact docs commit.
- [ ] Known IPC limitations are documented.
- [ ] All version surfaces say `0.0.4`.
- [ ] The release tag points to the exact tested commit.

## 4. Implementation phase 0 — freeze, inventory, and branch setup

### P0.1 Create the release branches

- [ ] Create `release/0.0.4` in `/Users/jinceyang/Desktop/codebase/merge/libpgo`; if it already exists, verify its base commit before making changes.
- [ ] Add or initialize `/Users/jinceyang/Desktop/codebase/merge/libpgo/docs` as the `git@github.com:annajcy/libpgo-doc.git` git submodule; do not use a separate docs worktree.
- [ ] Inside that submodule, check out branch `v0.0.4` at the frozen starting commit `ef0834c67a280d16836579942246129087c99e44` before editing content.
- [ ] Record the SSH URL and `branch = v0.0.4` in the superproject `.gitmodules`.
- [ ] Prepare `upstream-0.0.4` as the integration branch; do not merge into it until the working branch is reviewable.
- [ ] Prepare release-tracking issue text if useful; create or modify a GitHub issue only with explicit authorization.

Acceptance:

- Code work starts from `1670f6258d941d80cd9d718fc033fccc193d5e20`.
- Docs work is performed inside `/Users/jinceyang/Desktop/codebase/merge/libpgo/docs`, on the canonical `annajcy/libpgo-doc` `v0.0.4` branch, even though its starting content currently matches the experimental docs commit.
- 0.0.3 (main) code checkout is used as the changelog baseline in docs;

### P0.2 Capture baselines before edits

- [ ] Run the existing native build and tests in `/Users/jinceyang/Desktop/codebase/merge/libpgo`.
- [ ] Run the existing `tests/pypgo` suite in `/Users/jinceyang/Desktop/codebase/merge/libpgo`.
- [ ] Record existing failures separately from regressions introduced during cleanup.
- [ ] Inventory all tracked files under `/Users/jinceyang/Desktop/codebase/merge/libpgo/dist`.
- [ ] Record current example asset sizes and element counts.
- [ ] Record the current example inputs and the cubic generator parameters needed to reproduce test meshes.

### P0.3 Resolve high-risk integration paths with bounded spikes

Before broad refactoring:

- [ ] Prove one tiny tet + IPC static solve can be assembled from the existing IPC energy/derivative/CCD/max-step and static Newton infrastructure, produces a finite result, and does not require a new contact algorithm.
- [ ] If that spike exposes a missing core capability, stop and record the exact gap before claiming the five-case static matrix or expanding the implementation.
- [ ] On each platform, build one minimal wheel from the clean checkout and prove a fresh Conda environment containing GMP/MPFR/Imath/fmt plus the platform runtime can load it without vendored libraries, build-prefix paths, or undeclared host-package-manager paths.
- [ ] Locally build the `v0.0.4` docs content after replacing the experimental `environment.yml`/`pypgo-ci` assumptions; prove the simplified workflow has a viable source checkout and Pages deployment route.
- [ ] Record spike commands/results in the release evidence so later implementation does not silently change the accepted contracts.

## 5. Implementation phase 1 — library correctness and cleanup

### P1.1 Delete the obsolete `CIPCPotentialEnergy`

Target files:

- `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/core/contact/CIPC.h`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/core/contact/CIPC.cpp`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/tests/src/core/contact/cipcPotentialEnergy_gtest.cpp`
- Corresponding entries in `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/core/contact/CMakeLists.txt`
- Corresponding entries in `/Users/jinceyang/Desktop/codebase/merge/libpgo/tests/src/core/contact/CMakeLists.txt`

Checklist:

- [ ] Remove the legacy wrapper completely.
- [ ] Keep `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/core/contact/embeddedSurfaceIPCPotentialEnergy.*`.
- [ ] Keep `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/core/contact/mappedSurfacePotentialEnergy.*`.
- [ ] Keep `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/core/contact/ipc/`.
- [ ] Move any still-useful behavior test to `SurfaceIPCCore` or `EmbeddedSurfaceIPCPotentialEnergy`.
- [ ] Confirm there are no remaining includes, forward declarations, test helpers, or CMake source entries referring to `CIPCPotentialEnergy`.

Acceptance:

```bash
rg "CIPCPotentialEnergy|CIPC\\.h" /Users/jinceyang/Desktop/codebase/merge/libpgo/src /Users/jinceyang/Desktop/codebase/merge/libpgo/tests
```

returns no obsolete-wrapper references.

### P1.2 Remove scoped profiling

Remove:

- `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/core/profiling/`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/tests/src/core/profiling/`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/core/contact/ipc/profiling/surfaceIPCProfiling.h`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/tests/src/core/contact/cipcProfiling_gtest.cpp`
- Profiling source/test CMake entries.
- `ScopedProfileSection` includes and local scope objects.
- Runner config fields and logging related only to profiling.

Checklist:

- [ ] Remove profiling calls from IPC, solver, energy aggregation, deformation, and runner code.
- [ ] Remove `profiling` config parsing from `runIPCSim`.
- [ ] Remove profile summary emission and profile reset/enable control flow.
- [ ] Preserve computation and error handling after mechanically removing scopes.
- [ ] Do not replace this with a different always-on profiler.

Acceptance:

```bash
rg "ScopedProfileSection|scopedProfile|surfaceIPCProfiling|profiling summary" /Users/jinceyang/Desktop/codebase/merge/libpgo/src /Users/jinceyang/Desktop/codebase/merge/libpgo/tests
```

returns no release-code profiling framework references.

### P1.3 Fix the Hessian interface contract

Primary files:

- `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/core/nonlinearOptimization/potentialEnergy.h`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/core/nonlinearOptimization/potentialEnergies.cpp`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/core/contact/mappedSurfacePotentialEnergy.cpp`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/core/contact/mappedSurfacePotentialEnergy.h`

Contract:

- Fixed topology:
  - `createHessian()` creates the reusable sparse pattern.
  - `hessian()` updates values in that pattern.
  - `hessianDirect()` may use the base implementation.
- Dynamic topology:
  - `createHessian()` returns an empty matrix with the correct dimensions.
  - `hessianDirect()` produces the complete current sparse matrix.
  - Normal dynamic-Hessian execution must not rely on an exception.
- Aggregate:
  - `PotentialEnergies::init()` calls `createHessian()` only for fixed-topology terms.
  - Fixed patterns are mapped into the aggregate template once.
  - Dynamic matrices are computed and scattered on each evaluation.

Checklist:

- [ ] Stop unconditionally invoking `createHessian()` before checking `isHessianTopologyFixed()`.
- [ ] Make the dynamic mapped-surface implementation return a correctly sized empty template rather than throw.
- [ ] Ensure `PotentialEnergies::hessianDirect()` includes fixed and dynamic contributions exactly once.
- [ ] Verify zero-energy-coefficient terms do not force unnecessary dynamic evaluation.
- [ ] Add a mixed fixed + dynamic aggregate regression test.
- [ ] Test dynamic contact-pair pattern changes between two states.
- [ ] Test matrix dimensions, symmetry, and finite values.

### P1.4 Fix Newton symbolic-factorization lifetime

Primary file:

- `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/core/nonlinearOptimization/NewtonSolver.cpp`

Checklist:

- [ ] Remove the current fixed-topology `solver.reset()` after each factorization/solve.
- [ ] For fixed topology, analyze once and reuse symbolic analysis across Newton iterations.
- [ ] For dynamic topology, re-run analysis when the sparse pattern may have changed.
- [ ] Continue numerical factorization on every iteration.
- [ ] Preserve the MKL Pardiso, original Pardiso, and Eigen fallback branches.
- [ ] Add a regression test that executes multiple Newton iterations for fixed topology.
- [ ] Add or retain a dynamic-topology test whose Hessian pattern changes.
- [ ] Verify numerical results are unchanged within tolerance.

Release-blocking acceptance:

- Fixed-topology solves no longer discard their symbolic solver immediately.
- Dynamic-topology solves never reuse an invalid symbolic pattern.

### P1.5 Simplify maximum-step behavior

Primary files:

- `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/core/nonlinearOptimization/potentialEnergy.h`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/core/nonlinearOptimization/potentialEnergies.cpp`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/core/solidDeformationModel/deformationModelEnergy.*`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/core/contact/embeddedSurfaceIPCPotentialEnergy.*`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/core/contact/ipc/core/surfaceIPCCore.*`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/core/simulation/`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/tools/runSim/`

Checklist:

- [ ] Give `PotentialEnergy::computeMaxStepSize()` a default implementation returning `1.0`.
- [ ] Remove boilerplate overrides that only return `1.0`.
- [ ] Keep meaningful overrides for deformation inversion prevention and IPC CCD.
- [ ] Aggregate by taking the minimum scalar feasible step.
- [ ] Do not introduce `FeasibleStepKind`.
- [ ] Remove material/contact classification reporting.
- [ ] Remove clamp counters, per-solve minimum reports, line-search classification callbacks, and concrete-type casts that exist only to print those diagnostics.
- [ ] Retain useful exceptional warnings only when they identify invalid numerical input, not ordinary step clamping.
- [ ] Test no-limit energy returns `1.0`.
- [ ] Test tet inversion constrains the step.
- [ ] Test cubic inversion constrains the step.
- [ ] Test IPC CCD constrains the step.
- [ ] Test a mixed energy returns the minimum of material and IPC limits.

### P1.6 Harden `CubicMeshDeformationModel`

Primary files:

- `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/core/solidDeformationModel/cubicMeshDeformationModel.h`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/core/solidDeformationModel/cubicMeshDeformationModel.cpp`

Checklist:

- [ ] Replace the raw internal pointer with `std::unique_ptr<CubicMeshDeformationModelInternal>`.
- [ ] Explicitly delete copy construction and copy assignment.
- [ ] Support move only if the base class and owning call sites make it safe; otherwise explicitly delete move too.
- [ ] Replace repeated `8`, `24`, and integration-point counts with named internal `constexpr` constants.
- [ ] At construction, validate every element and Gauss point.
- [ ] Require `det(Dm)` to be finite and greater than a scale-aware positive threshold.
- [ ] Do not apply `abs(det(Dm))` to hide inverted rest elements.
- [ ] Include element index and Gauss-point index in validation errors.
- [ ] Validate rest derivatives and inverse matrices are finite.
- [ ] Add tests for valid cube, inverted element, degenerate element, and non-finite input.
- [ ] Retain finite-difference energy/gradient/Hessian tests.

### P1.7 Merge tet/cubic `SimulationMesh` conversion

Primary files:

- `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/core/solidDeformationModel/simulationMesh.h`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/core/solidDeformationModel/simulationMesh.cpp`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/tools/runSim/runSimFEMSetup.cpp`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/tools/runSim/runIPCSimSetup.cpp`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/c/pgo_c.cpp`

Implementation:

```cpp
std::unique_ptr<SimulationMesh>
loadVolumeMesh(const VolumetricMeshes::VolumetricMesh &mesh);
```

Checklist:

- [ ] Implement the common conversion using the volumetric mesh element type, vertex count, element-vertex count, connectivity, and materials.
- [ ] Map tet to `SimulationMeshType::TET`.
- [ ] Map cubic to `SimulationMeshType::CUBIC`.
- [ ] Reject unsupported element types explicitly.
- [ ] Preserve `loadTetMesh()` and `loadCubicMesh()` as thin compatibility wrappers.
- [ ] Have legacy raw-pointer wrappers delegate to `loadVolumeMesh(...).release()` only at the compatibility boundary.
- [ ] Update new code to prefer `std::unique_ptr`.
- [ ] Verify material values and per-element material mapping are unchanged.
- [ ] Add equivalence tests for wrapper versus generic loader.

### P1.8 Put IPC validation in core

Primary files:

- `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/core/contact/ipc/core/surfaceIPCCore.cpp`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/core/contact/ipc/core/surfaceIPCCore.h`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/core/contact/ipc/topology/surfaceIPCTopology.cpp`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/core/contact/ipc/broadPhase/spatialHashGrid.cpp`

Required validation:

- [ ] `dhat` is finite and greater than zero.
- [ ] `kappa` is finite and greater than zero.
- [ ] `eps_ee` is finite and nonnegative.
- [ ] `slackness` is finite and satisfies `0 < slackness <= 1`.
- [ ] `V` has shape `N x 3`.
- [ ] `F` has shape `M x 3`.
- [ ] All vertex coordinates are finite.
- [ ] All triangle indices are in `[0, N)`.
- [ ] Each triangle uses three distinct vertex indices.
- [ ] Triangle areas are finite and above a scale-aware nondegeneracy threshold.
- [ ] Position and direction vectors have length `3N`.
- [ ] Position and direction vectors contain no NaN or infinity.
- [ ] The spatial-hash cell size is finite and greater than zero.
- [ ] State is invalidated after parameter or topology changes.
- [ ] Failures use `std::invalid_argument` for invalid public input and `std::logic_error` for missing required preparation state.

Tests:

- [ ] One negative test per rule above.
- [ ] Boundary tests for `eps_ee == 0` and `slackness == 1`.
- [ ] Rejected parameter updates do not partially mutate core state.
- [ ] Valid updates invalidate cached active pairs.

### P1.9 Retain optional components without expanding them

#### backward-cpp

- [ ] Keep `/Users/jinceyang/Desktop/codebase/merge/libpgo/CMakeModules/third-party/backward.cmake`.
- [ ] Replace `GIT_TAG master` with an exact tag or commit.
- [ ] Do not add new backward-cpp instrumentation.

#### CUDA utilities

- [ ] Keep `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/core/cudaUtilities`.
- [ ] Keep CUDA disabled in default release and wheel presets.
- [ ] Ensure a non-CUDA build never fetches or requires CUDA-only dependencies.
- [ ] If a CUDA runner/toolkit is available, run a non-blocking compile check; lack of a GPU runner does not block 0.0.4.

#### Native Pardiso

- [ ] Keep `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/core/eigenSupport/EigenOrigPardisoSupport.*`.
- [ ] Keep original Pardiso disabled unless explicitly configured.
- [ ] Do not bundle proprietary Pardiso libraries in wheels.
- [ ] Preserve its solver branch while fixing symbolic-factorization lifetime.

### P1.10 Small source and CMake cleanup

- [ ] Change `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/c/pgo_c.h` `pgo_create_tetmeshgeo_from_file(char *)` to `const char *`.
- [ ] Remove the corresponding `const_cast` from Python binding code.
- [ ] Add a small NumPy buffer-validation helper in `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/python/pypgo/` for dtype, ndim, shape, divisibility, and contiguity checks.
- [ ] Use the helper to replace duplicated ad hoc checks in `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/python/pypgo/pypgo.cpp`.
- [ ] Replace stderr-only validation failures with Python exceptions.
- [ ] Correct internal spelling `n_repreat` to `n_repeat`.
- [ ] Correct internal spelling `quertPtInfo` to `queryPtInfo`.
- [ ] Remove the duplicate `basicIO.h` include in `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/c/pgo_c.cpp`.
- [ ] Fix `pgo_smooth_rs_energy_hess()` so `values[inc]` receives the original `double` Hessian value rather than an `int` cast.
- [ ] Add a C API Hessian regression whose expected matrix includes finite noninteger entries.
- [ ] Correct the `pog_c` spelling in `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/c/CMakeLists.txt`.
- [ ] Remove trailing whitespace only in touched source and CMake files.
- [ ] Do not reformat generated IPC code or generated `.veg` files.

## 6. Implementation phase 2 — unified runner, C API, Python API, and config

### P2.1 Build one shared simulation-runner layer

#### Why the current layout must change

- `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/tools/runSim/runSim.cpp`,
  `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/tools/runSim/runIPCSim.cpp`, and
  `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/tools/runSim/runShellSim.cpp` are executables that each own CLI parsing, config access, setup, simulation, and output mutation.
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/c/pgo_c.cpp` contains a fourth, older tet + sampled implementation.
- The executable split is currently by contact/domain, while the 0.0.4 public split is by simulation mode. Adding two more large `main()` functions would multiply implementations instead of unifying them.
- A reusable library placed under `tools` would create an avoidable CMake ordering problem because `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/CMakeLists.txt` adds `c` before `tools`.

#### Target ownership

Create a non-installed application-layer library named `simulationRunner` outside the
CLI directory:

```text
src/simulationRunner/
  CMakeLists.txt
  simulationRunner.h/.cpp          # typed entry point and top-level pipeline
  simulationConfig.h/.cpp          # parse once, normalize, validate, classify
  sampledSimulationBackend.h/.cpp  # sampled setup and execution
  ipcSimulationBackend.h/.cpp      # IPC volume/shell setup and execution
  simulationOutput.h/.cpp          # new-run/resume output transaction

src/tools/runSim/
  CMakeLists.txt
  runDynamicSim.cpp                # argparse + result reporting only
  runStaticSim.cpp                 # argparse + result reporting only
```

The exact private filenames may be adjusted while implementing, but these ownership
boundaries are required:

- [ ] `src/tools/runSim` owns only CLI syntax, optional CLI-only logging setup, human-readable result reporting, and process exit codes.
- [ ] `simulationRunner` owns config loading, typed validation, route selection, setup, execution, output safety, and conversion of internal exceptions to a runner result.
- [ ] Sampled and IPC backend code owns solver/contact-specific setup and stepping, but does not parse CLI arguments, select another backend, clear output, or implement public API wrappers.
- [ ] `simulationOutput` implements the filesystem and resume contract in P2.5; backends receive an acquired output session and cannot bypass it with recursive deletion or direct overwrite.
- [ ] C and Python wrappers own ABI/language-boundary conversion only; they contain no mesh/contact/solver implementation.

#### One typed entry point and result

Use this internal contract:

```cpp
enum class SimulationMode { Dynamic, Static };
enum class RunDisposition { NewRun, Resume };

struct RunSimulationOptions
{
  std::optional<SimulationMode> expectedMode; // unset for C/Python
  RunDisposition disposition = RunDisposition::NewRun;
};

struct RunSimulationResult
{
  int code;
  std::string message;
};

RunSimulationResult runSimulationFromConfig(
  const std::filesystem::path &configPath,
  const RunSimulationOptions &options = {});
```

This is an internal C++ API, not a new installed public ABI. `RunSimulationResult`
centralizes an actionable failure message and integer status; callers must not need
to catch backend exceptions.

- [ ] `runDynamicSim <config>` passes `expectedMode = Dynamic` and `NewRun`.
- [ ] `runDynamicSim <config> --resume` passes `expectedMode = Dynamic` and `Resume`.
- [ ] `runStaticSim <config>` passes `expectedMode = Static` and `NewRun`; its argument parser rejects `--resume`.
- [ ] C and Python call the same function with no `expectedMode` and `NewRun`, so `sim-type` in the config selects static or dynamic while both language APIs remain new-run-only.
- [ ] Reject a CLI/config mode mismatch before creating or modifying output.
- [ ] Return `code == 0` only for a completed run; config, setup, solver, and output failures return nonzero with one actionable `message`.
- [ ] CLI callers print the result message and return its code; the C wrapper logs the message and returns its code; Python preserves its existing integer result contract through the C wrapper.

#### Parse once, then dispatch

The shared path must be visibly linear:

```text
CLI / C / Python
        |
        v
runSimulationFromConfig
  -> parse and validate one typed SimulationRequest
  -> validate expected mode and new-run/resume combination
  -> select backend from contact model
  -> prepare the selected domain without output mutation
  -> acquire the output transaction
  -> execute dynamic or static mode
  -> commit output and return RunSimulationResult
```

The typed request records at least the normalized config path, mode, domain
(`tet`, `cubic`, or `shell`), contact model (`sampled` or `ipc`), output path, and
validated common/backend settings. It is the only routing source after parsing.

- [ ] Parse the JSON once; do not reopen it independently in the dispatcher and backend.
- [ ] Classify mode, domain, and contact before expensive setup or output mutation.
- [ ] Reject zero/multiple domains, unsupported shell + sampled, unknown contact, and unsupported IPC external objects during typed validation.
- [ ] Select sampled versus IPC exactly once in the dispatcher.
- [ ] Within the selected backend, use the typed domain to choose volume or shell setup and the typed mode to choose dynamic or static execution.
- [ ] Keep the existing detailed FEM, sampled-contact, IPC, and shell setup logic where practical, but make it consume the typed request instead of reaching back into CLI state.
- [ ] Do not create ten independent functions for the five domain/contact cases in two modes; share preparation and mode execution wherever the underlying implementation is common.

#### Migrate, do not duplicate

- [ ] Extract the sampled flow from `runSim.cpp` into the sampled backend.
- [ ] Migrate the dynamic IPC flow from `runIPCSim.cpp` into the IPC backend; leave no second IPC implementation.
- [ ] Migrate shell IPC setup/execution from `runShellSim.cpp` into the same IPC backend; do not add sampled-shell execution.
- [ ] Add static IPC by composing the existing IPC energy/gradient/Hessian, CCD/maximum-step, constraints, and Newton solver infrastructure without the dynamic time integrator.
- [ ] Make static IPC support tet, cubic, and shell while preserving the 0.0.4 limitation on arbitrary `external-objects` and the supported self-contact/floor configuration.
- [ ] Add focused static IPC finite-output, convergence, no-NaN, maximum-step, and result-file tests before counting the five static matrix cells as supported.
- [ ] Delete the duplicated simulation body from `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/c/pgo_c.cpp`.
- [ ] Remove Release-build `std::cin.get()` pauses so every new entry point is noninteractive.
- [ ] Remove the `runSim`, `runIPCSim`, and `runShellSim` public targets after their implementations have been migrated; executable-name compatibility is not required for 0.0.4.

#### CMake target graph

- [ ] Add `src/simulationRunner` after all required core libraries and before `c`, `python`, and `tools` in `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/CMakeLists.txt`.
- [ ] Link `simulationRunner` only to the specific core targets it uses; do not add it to `PGO_GLOBAL_LIBRARY_TARGETS` and thereby link it into unrelated tools/tests.
- [ ] Link `pgo_c`, `pgo_c_static`, `runDynamicSim`, and `runStaticSim` explicitly to `simulationRunner`.
- [ ] Keep `argparse` linked only to the two executable targets.
- [ ] Verify the graph is one-way: core libraries → `simulationRunner` → C/Python/CLI; no runner → `pgo_c` edge and no dependency cycle.

#### Completion evidence

- [ ] A routing unit test covers every supported and rejected `(mode, domain, contact, disposition)` combination without running a long simulation.
- [ ] Mode mismatch, invalid config, and output-preflight tests prove failure occurs before output mutation.
- [ ] Focused integration tests prove the CLI, C, and Python surfaces reach the same sampled/IPC backend implementations.
- [ ] Source/CMake inspection confirms there is one config parser path, one dispatcher, one sampled backend, one IPC backend, and no simulation body in a public wrapper or `main()`.

This application-layer extraction is the only runner-level structural adjustment
planned for 0.0.4; it does not authorize a core solver, contact, mesh, or package
redesign.

### P2.2 C API alignment

Files:

- `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/c/pgo_c.h`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/c/pgo_c.cpp`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/c/pgo_c.symbol.txt`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/c/pgo_c.version`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/tests/pgo_c_gtest.cpp`

Checklist:

- [ ] Preserve the signature and symbol `pgo_run_sim_from_config`.
- [ ] Route it to the unified dispatcher.
- [ ] Return `0` on success and nonzero on config/setup/simulation failure.
- [ ] Treat a missing or unreadable config as failure; remove the current success return on `ConfigFileJSON::open()` failure from the C API and both CLI paths.
- [ ] Do not let C++ exceptions cross the C ABI boundary.
- [ ] Log or retain an actionable error message before returning failure.
- [ ] Verify the symbol remains exported on Linux, macOS, and Windows.
- [ ] Preserve the `LIBPGO_C_ABI_1.0` namespace because no binary-breaking C signature is introduced.
- [ ] Test all five required domain/contact combinations through the shared C library in Linux server Tier 1.
- [ ] Test repeated successful calls, a failure followed by success, and two different sequential configs in one process.

### P2.3 Python API alignment

Files:

- `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/python/pypgo/pypgo.cpp`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/python/pypgo/pypgo_main.cpp`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/tests/pypgo/`

Checklist:

- [ ] Preserve `pypgo.run_sim_from_config`.
- [ ] Continue routing through the C API or the exact same reusable C++ dispatcher.
- [ ] Preserve the existing integer result contract: `0` on success and nonzero on failure; do not silently convert runner failures into Python success.
- [ ] Keep the 0.0.4 Python binding new-run-only; do not infer resume from an existing output directory or the config file.
- [ ] Test path-like and string input as supported by the existing binding surface.
- [ ] Test all five required matrix cases through Python in Linux server Tier 1.
- [ ] Test invalid config behavior.
- [ ] Test missing config, existing output path, failure-then-success, and sequential multi-run behavior.
- [ ] Verify `pypgo.__version__`.
- [ ] Do not migrate `/Users/jinceyang/Desktop/codebase/libpgo/pypgo` or `/Users/jinceyang/Desktop/codebase/libpgo/src/python/pypgo` wholesale.

### P2.4 Config schema and compatibility

Checklist:

- [ ] Add `contact-model` with supported values `sampled` and `ipc`.
- [ ] For tet/cubic volume configs, missing `contact-model` means `sampled`; shell configs must select IPC explicitly.
- [ ] Reject unknown values.
- [ ] Require exactly one supported domain: `tet-mesh`, `cubic-mesh`, or a valid shell schema.
- [ ] Require `sim-type` to be exactly `dynamic` or `static`; validate it against the selected CLI mode.
- [ ] Detect a mismatched `.veg` element type and report the config key, expected type, actual type, and filename.
- [ ] Keep relative paths resolved relative to the config file.
- [ ] Require a finite positive timestep.
- [ ] Require finite scale and reject zero scale.
- [ ] Require `num-timestep >= 1`.
- [ ] Require `dump-interval >= 1`.
- [ ] Require finite positive solver tolerance and positive solver iteration limit.
- [ ] Validate gravity, initial velocity/displacement, damping, and movement vectors for expected length and finite values.
- [ ] Validate IPC-only parameters only when IPC is selected.
- [ ] Reject `external-objects` for IPC with the documented 0.0.4 limitation.
- [ ] Reject the legacy `restart-from-u` config key with an actionable message directing CLI users to `--resume`.
- [ ] Require shell configs to select `contact-model: "ipc"` explicitly; reject shell + sampled and missing `contact-model` for shell with an actionable error.
- [ ] Preserve the existing shell IPC setup and its `koiter-stvk`/`ipc-heuristic` contract.
- [ ] Keep IPC floor configuration explicit (`use-floor`, axis, height, kappa).
- [ ] Keep sampled contact parameters on the sampled path.
- [ ] Do not silently reinterpret IPC keys as sampled parameters or vice versa.

Required config parser tests:

- [ ] Explicit sampled.
- [ ] Explicit IPC.
- [ ] Unknown contact model.
- [ ] Volume config with missing contact model falls back to sampled.
- [ ] Shell config with missing contact model or explicit sampled contact is rejected with guidance to select IPC.
- [ ] Both tet and cubic keys present.
- [ ] Neither tet nor cubic key present with a valid shell config.
- [ ] Neither volume keys nor a valid shell config present.
- [ ] Tet key pointing to cubic mesh.
- [ ] Cubic key pointing to tet mesh.
- [ ] IPC with unsupported `external-objects`.
- [ ] IPC floor validation.

### P2.5 Dynamic-directory, static-file, and resume safety contract

The shared runner library, C API, Python API, `runDynamicSim`, and `runStaticSim` must never recursively delete a configured output path. Static and dynamic simulations have intentionally different output types.

Shared new-run checks:

- [ ] Resolve `output` relative to the config file and normalize it before any filesystem mutation.
- [ ] Require the output parent to exist and be a real directory.
- [ ] Require the final output path to be absent, including a symlink or broken symlink at that path.
- [ ] If the output path already exists as a file, directory, or symlink, fail with a nonzero result and an actionable message.
- [ ] Remove `clearOutputDirectory()` and every recursive-delete call from simulation execution.
- [ ] Create output only after config and setup validation have succeeded as far as practical.
- [ ] C and Python remain new-run-only because their frozen entry points do not carry CLI flags.

New dynamic run:

- [ ] `runDynamicSim` and a dynamic C/Python config create exactly one new final output directory and require directory creation to report success.
- [ ] Write frames, checkpoints, and a versioned run manifest only inside that directory.
- [ ] If execution fails after creating the new directory, leave it for diagnosis; never recursively clean it.

New static run:

- [ ] `runStaticSim` and a static C/Python config treat `output` as one new result-file path, preserving the existing static-config meaning.
- [ ] Write the result to a temporary sibling file, validate it, and atomically rename it to the requested absent final path.
- [ ] On failure, remove only the uniquely named temporary file created by this invocation; never delete or modify any pre-existing path.
- [ ] Reject an output path whose parent is absent, whose final component exists, or whose final component is a directory/symlink.

Dynamic CLI resume:

- [ ] Enter resume mode only when `runDynamicSim` receives `--resume`; `runStaticSim --resume` is a usage error.
- [ ] Require the configured output path to exist as a real directory, not a file or symlink.
- [ ] Require a versioned run manifest written by a new dynamic run that records the config hash, input identities/hashes, domain/contact mode, state dimensions, checkpoint schema, and last completed checkpoint/frame.
- [ ] Verify the requested config and required inputs match the run manifest before loading state.
- [ ] Locate the latest complete valid checkpoint and continue from the next frame.
- [ ] Write every checkpoint/frame to a temporary sibling, validate it, and atomically rename it before atomically advancing the manifest.
- [ ] Fail if no valid checkpoint exists, the manifest/config/input identity does not match, or the next output name would overwrite an existing file.
- [ ] Never remove, truncate, rename, or overwrite an existing frame/checkpoint during resume.
- [ ] Reject `restart-from-u` in config rather than maintaining two resume controls; report that `runDynamicSim --resume` is required.

Required tests:

- [ ] A new dynamic child directory under an existing parent succeeds.
- [ ] A new static result file under an existing parent succeeds and is atomically installed.
- [ ] Existing empty/nonempty directories, files, and symlinks fail and remain unchanged in both applicable modes.
- [ ] Filesystem root, home, repository root, and config directory fail because they already exist or are invalid result-file targets.
- [ ] Normalization such as `new/../existing` cannot bypass the existence check.
- [ ] A competing creation between validation and final creation/rename fails safely.
- [ ] Sampled, IPC, tet, cubic, shell, C, Python, and CLI surfaces return consistent failure results.
- [ ] `--resume` with valid sampled and IPC dynamic runs continues after the latest checkpoint without modifying earlier output.
- [ ] `--resume` coverage includes tet, cubic, and shell domains.
- [ ] `--resume` fails for a missing directory, static config, `runStaticSim`, non-libpgo directory, symlink, manifest/config mismatch, state-dimension mismatch, missing/corrupt checkpoint, and next-output collision.
- [ ] The C and Python APIs fail on an existing output path rather than implicitly resuming.

## 7. Implementation phase 3 — examples and mesh assets

### P3.1 Final config/asset-separated layout

Use `/Users/jinceyang/Desktop/codebase/main/libpgo/examples` as the canonical source for 0.0.3 volume inputs and the current 0.0.4 target as the source for retained shell/IPC inputs. Migrate the selected content into this single final layout:

```text
examples/
├── README.md
├── assets/
│   ├── common/
│   ├── media/
│   ├── volume/
│   │   ├── box/
│   │   ├── box-with-sphere/
│   │   ├── bunny/
│   │   ├── dragon/
│   │   └── dragon-dyn/
│   └── shell/
├── configs/
│   ├── volume/
│   │   ├── box/
│   │   ├── box-with-sphere/
│   │   ├── bunny/
│   │   ├── dragon/
│   │   └── dragon-dyn/
│   └── shell/
├── generated/
│   └── cubic/
└── scripts/
    └── generate_cubic_veg.py
```

Separation rules:

- [ ] Store each checked-in OBJ, tet `.veg`, fixed-vertex list, animation input, and other static input exactly once under `examples/assets/`.
- [ ] Put shared floor/external-object inputs such as `bottom.obj` under `examples/assets/common/`.
- [ ] Put canonical retained GIF/PNG media under `examples/assets/media/` without duplicating it by contact or mesh type.
- [ ] Put every simulation/animation JSON config under `examples/configs/`; JSON files are not assets.
- [ ] Put sampled and IPC volume configs beside each other under `examples/configs/volume/<scene>/`.
- [ ] Put dynamic/static shell IPC configs under `examples/configs/shell/`; both reference the one canonical `examples/assets/shell/` tree.
- [ ] Keep shell IPC functionality and content, but do not preserve its old directory paths or retain a sampled-shell config as a release example.
- [ ] Remove old duplicate trees such as `examples/cubic/`, `examples/ipc/`, and old scene directories after all selected content has been migrated and references have been updated.
- [ ] Do not retain compatibility JSON files, symlinks, forwarding files, or duplicate assets at the old paths.
- [ ] Run a content-hash duplicate audit over `examples/assets/` and justify or remove every duplicate.
- [ ] Update tests, README commands, docs, and CMake test definitions to use only the new paths.
- [ ] Add directory-level ignore rules for `examples/generated/` and keep generated cubic meshes out of git.

### P3.2 Generated cubic `.veg` layout

Add a checked-in Python generator, preferably:

```text
examples/scripts/generate_cubic_veg.py
```

The generator contract:

- [ ] Accept the cubic mesher executable/build directory, input surface OBJ, output directory, resolution/voxel parameters, and optional selected scenes as command-line arguments.
- [ ] Produce clearly named `<scene>-cubic.veg` files under an explicit output directory.
- [ ] Use `examples/generated/cubic/` as the documented local default and ignore that directory in git.
- [ ] Refuse to overwrite an existing output file or directory unless the caller explicitly chooses a different new path.
- [ ] Use canonical surface assets directly instead of copying OBJ files into a second tree.
- [ ] Print or write a machine-readable manifest containing input hashes, command arguments, mesher/code revision, element/vertex counts, bounds, and validation results.
- [ ] Keep existing tet `.veg` files.
- [ ] Do not commit generated cubic `.veg` files.
- [ ] Make cubic config filenames identify both mesh and contact choices and document that generation must run first.

Recommended dynamic config naming:

```text
<scene>-tet-sampled.json
<scene>-tet-ipc.json
<scene>-cubic-sampled.json
<scene>-cubic-ipc.json
shell-dynamic-ipc.json
shell-static-ipc.json
```

Static configs use the same domain/contact suffixes with an explicit `static` component, for example `<scene>-static-cubic-ipc.json`; dynamic configs include `dynamic` where needed to prevent ambiguity.

Path migration contract:

- [ ] Switch directly to the new paths in 0.0.4; old example paths are intentionally unsupported.
- [ ] Do not keep `examples/box/box.json`, `examples/shell/shell.json`, `examples/ipc/shell/shell-ipc.json`, or equivalent compatibility copies.
- [ ] Resolve checked-in assets relative to the config file, for example from `examples/configs/volume/box/` to `examples/assets/volume/box/`.
- [ ] Treat checked-in cubic JSON files as templates: after mesh generation, the CI/server driver rewrites generated-mesh and output paths to absolute paths and saves the exact materialized config under the report directory.
- [ ] Run and hash the materialized config rather than relying on an ambient generated-asset root or current working directory.
- [ ] Add a concise old-to-new path table to the example README and release notes; this is migration documentation, not backward compatibility.

### P3.3 Cubic generator validation

Useful implementation references, not wholesale migration sources:

- `/Users/jinceyang/Desktop/codebase/libpgo/examples/assets/generate_veg_assets.py`
- `/Users/jinceyang/Desktop/codebase/libpgo/examples/experiments/tricubic_hermite_fem/mesh/mesh_quality.py`
- `/Users/jinceyang/Desktop/codebase/libpgo/examples/experiments/tricubic_hermite_fem/mesh/generate_cubic_mesh.py`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/tools/cubicMesher/`

Checklist for meshes generated during bounded automated CI tests or the owner-run server validation:

- [ ] Generate cubic meshes from the canonical surface OBJ files into a fresh directory outside the source tree.
- [ ] Store the exact generation command, script and mesher revisions, resolution, occupancy rule, scale, input hashes, and output statistics.
- [ ] Ensure every element index is valid.
- [ ] Ensure every cubic element has eight distinct vertex indices.
- [ ] Ensure every Gauss-point `det(Dm)` is finite and above the model threshold.
- [ ] Ensure there are no isolated vertices or unintended disconnected components.
- [ ] Match the canonical source bounding box and scale.
- [ ] Run at least one short static or dynamic deformation smoke per regenerated cubic mesh.
- [ ] Verify no rest inversion, runtime inversion, NaN, or solver failure occurs.
- [ ] Verify a second clean run with the same inputs produces the same topology/content hash, or document and test the narrowest deterministic invariant the mesher can guarantee.
- [ ] Keep generated `.veg` files and reports as CI artifacts or in the owner’s manual server result bundle, not git files.

### P3.4 Example matrix acceptance

For each dynamic canonical scene selected for the public matrix:

- [ ] Tet + sampled config resolves and starts.
- [ ] Tet + IPC config resolves and starts.
- [ ] Generate the cubic `.veg` into a fresh test output directory.
- [ ] Cubic + sampled config resolves the generated mesh and starts.
- [ ] Cubic + IPC config resolves the generated mesh and starts.

For IPC scenes that previously used `external-objects`:

- [ ] Replace the external plane object with the supported IPC floor config rather than pretending external IPC is supported.
- [ ] Keep the sampled config’s original external-object behavior where appropriate.
- [ ] Explain this difference in the example README and release notes.

For retained shell cases:

- [ ] `examples/configs/shell/shell-dynamic-ipc.json` and `shell-static-ipc.json` resolve the canonical shell inputs and execute the retained shell IPC path.
- [ ] Dynamic shell IPC executes through `runDynamicSim`, C, and Python; static shell IPC executes through `runStaticSim`, C, and Python.
- [ ] No shell OBJ, fixed-vertex list, or animation input is duplicated between dynamic and static configs.

Short smoke config design:

- [ ] Use tiny generated meshes or override timestep count/output path for CI.
- [ ] Materialize every dynamic smoke config with a unique nonexistent output directory; an existing output path must make the run fail.
- [ ] Materialize every static smoke config with a unique nonexistent output file; an existing output path must make the run fail.
- [ ] Run enough steps to execute setup, energy/gradient/Hessian, maximum-step, factorization, and output.
- [ ] Do not run the full 2,000–4,000-frame showcase configs in CI.
- [ ] Hosted CI runs only a bounded subset; the Linux server driver runs the complete short API matrix and representative longer matrix.

## 8. Implementation phase 4 — build system, dependencies, and wheels

### P4.1 Dedicated pypgo configure presets

Reference implementation:

- `/Users/jinceyang/Desktop/codebase/libpgo/CMakePresets.json`
- `/Users/jinceyang/Desktop/codebase/libpgo/setup.py`

Target implementation:

- `/Users/jinceyang/Desktop/codebase/merge/libpgo/CMakePresets.json`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/setup.py`

Add configure presets:

```text
pypgo-wheel-no-mkl
pypgo-wheel-mkl
```

Checklist:

- [ ] Factor common Python options into a hidden preset.
- [ ] Enable `PGO_ENABLE_PYTHON`.
- [ ] Enable only the subprojects required to build `pypgo`.
- [ ] Keep the target name `pypgo`.
- [ ] Make MKL choice explicit in the two presets.
- [ ] Set `PGO_USE_MKL=ON`, `PGO_HAS_ORIG_PARDISO=OFF`, and dynamic MKL linkage in `pypgo-wheel-mkl`.
- [ ] Set `PGO_CHECK_CONDA=ON` in release CI presets and resolve MKL/TBB from the active Conda prefix rather than from an undeclared machine installation.
- [ ] Resolve GMP/GMPXX, MPFR, Imath, fmt, and TBB from the active `$CONDA_PREFIX` on every release platform. Do not resolve any of these non-system dependencies from `/opt/homebrew`, `/usr/local`, a system package manager, a developer cache, or another undeclared host prefix.
- [ ] Before compilation, record each of GMP/GMPXX, MPFR, Imath, fmt, and TBB's resolved include path, library path, and CMake package directory. Release configuration fails unless every non-system dependency path is under the active Conda prefix.
- [ ] Treat fmt as an explicit Conda build and runtime dependency unless it is changed to a fully static, non-exported dependency; record and audit the chosen policy in every wheel job.
- [ ] Release wheel configuration must not inherit a developer's `CMAKE_PREFIX_PATH`, `*_DIR`, compiler/linker flags, or cached dependency locations. CI supplies approved Conda paths and CMake options explicitly.
- [ ] Preserve and surface the existing `MKL_THREADING=tbb_thread` selection from `CMakeModules/third-party/mkl.cmake`; fail configuration if the selected MKL target does not use the TBB threading layer.
- [ ] Set `PGO_USE_MKL=OFF` in `pypgo-wheel-no-mkl`.
- [ ] Extend the Conda TBB discovery path to macOS (`$CONDA_PREFIX/lib/cmake/TBB`) or pass an explicit equivalent such as `TBB_DIR`; do not silently fetch a second oneTBB when the CI environment already provides TBB.
- [ ] Add an explicit macOS BLAS policy before SuiteSparse is configured: force `BLA_VENDOR=Apple`, use the LP64 interface, and make configuration fail if BLAS/LAPACK do not resolve to the system Accelerate framework.
- [ ] Replace the unconditional `set(BLA_VENDOR "" ... FORCE)` in `CMakeModules/third-party/suitesparse.cmake` with a platform-aware setting. On macOS it must preserve `Apple`; on Linux/Windows MKL builds it must preserve the approved MKL vendor/interface instead of resetting the caller's choice.
- [ ] Ensure SuiteSparse and downstream targets consume the same resolved `BLAS::BLAS` and `LAPACK::LAPACK` targets; do not allow a nested SuiteSparse configure to discover a different provider.
- [ ] Keep CUDA and native Pardiso off in wheel presets.
- [ ] Let `setup.py` read `PYPGO_CMAKE_PRESET`.
- [ ] Make CI set `PYPGO_CMAKE_PRESET` explicitly.
- [ ] Preserve a sensible local default compatible with the old platform behavior: no MKL on macOS, MKL on supported Linux/Windows environments.
- [ ] Configure into the setuptools-provided temporary build directory rather than a shared, stale source-tree cache.
- [ ] Continue supporting `CMAKE_ARGS`, parsed safely with `shlex.split`.
- [ ] Preserve macOS `ARCHFLAGS`.
- [ ] Preserve MSVC architecture/config handling.
- [ ] Preserve `CMAKE_BUILD_PARALLEL_LEVEL`.
- [ ] Do not depend on a developer’s existing `build/` cache.

### P4.2 Version and package metadata

Files requiring audit:

- `/Users/jinceyang/Desktop/codebase/merge/libpgo/setup.py`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/pyproject.toml`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/c/CMakeLists.txt`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/python/pypgo/CMakeLists.txt`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/python/pypgo/pypgo_main.cpp`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/README.md`

Checklist:

- [ ] Change Python distribution version from `0.0.3` to `0.0.4`.
- [ ] Change the C shared-library `VERSION` from `0.0.2` to `0.0.4`.
- [ ] Keep C `SOVERSION 0` unless an actual binary-incompatible change is introduced.
- [ ] Install both public headers, `pgo_c.h` and the included `pgo_c_def.h`.
- [ ] Generate and install `pgoConfigVersion.cmake` with compatibility appropriate for the retained `SOVERSION`.
- [ ] Ensure the compiled Python module receives `PYPGO_VERSION_INFO=0.0.4`.
- [ ] Ensure installed `pypgo.__version__` is exactly `0.0.4`.
- [ ] Ensure wheel filenames contain `0.0.4`.
- [ ] Add/verify package description, README content type, license file, project URL, and supported Python metadata.
- [ ] Audit and declare actual Python runtime dependencies; in particular, verify whether NumPy is required at import time or API-call time instead of leaving `install_requires=[]` without evidence.
- [ ] Define the supported runtime as a documented Conda/Miniforge environment followed by installation of the downloaded CI wheel: all platforms install GMP, MPFR, Imath, fmt, and TBB from `conda-forge`; Linux/Windows additionally install MKL-selected BLAS/LAPACK; on macOS libpgo links system Accelerate while the environment installs `_newaccelerate` BLAS/LAPACK packages so Python/NumPy use the same backend family.
- [ ] Remove release installation instructions that use `apt install libgmp-dev libmpfr-dev` or `brew install gmp mpfr imath`; provide one platform-adjusted Conda environment recipe instead.
- [ ] State explicitly that pip cannot install Conda packages, so these wheels are not claimed to work in an otherwise empty `venv`.
- [ ] Keep Conda-provided GMP, MPFR, Imath, fmt, MKL, and TBB external to the wheels, reject accidental vendoring during package audit, and document the exact Conda packages/channels required before installing the wheel.
- [ ] Define the release wheel guarantee explicitly: CPython 3.12 wheels tagged `linux_x86_64`, macOS arm64, and `win_amd64`.
- [ ] Keep `python_requires >= 3.9` only if source-build CI verifies the claimed range; otherwise narrow the metadata to the range actually tested.

### P4.3 Wheel-only source provenance

- [ ] Do not build, upload, or document an sdist for 0.0.4.
- [ ] Build each wheel directly from the exact checked-out commit after that platform's native and Python tests pass.
- [ ] Record `git rev-parse HEAD`, workflow run ID, runner image, CMake cache, dependency evidence, and wheel SHA-256 beside each artifact.
- [ ] Require a clean tracked worktree before packaging so uncommitted source changes cannot enter a wheel.

### P4.4 Pin third-party dependencies

Mandatory floating-reference fixes in `/Users/jinceyang/Desktop/codebase/merge/libpgo/CMakeModules/third-party/`:

- [ ] `backward.cmake`: replace `master`.
- [ ] `libigl.cmake`: replace `main`.
- [ ] `cuCollections.cmake`: replace `dev`.
Additional audit:

- [ ] Verify every other FetchContent URL is a fixed release URL or exact commit.
- [ ] Do not upgrade dependency versions without a release need.
- [ ] Record the final dependency table in docs.
- [ ] Audit dependency licenses and retain required notices for source and vendored wheel contents.
- [ ] Add archive hashes where practical.
- [ ] Ensure optional dependencies are fetched only when their feature is enabled.
- [ ] Pin the Miniforge installer version, GitHub Actions revisions, Python version, build tools, and critical Conda dependency versions sufficiently to avoid `latest` changing the release build unexpectedly.
- [ ] Use `conda-forge` with strict channel priority and remove/disable default channels in release CI.
- [ ] Retain `conda list --explicit`, `conda list`, channel configuration, and relevant CMake cache entries as release-build evidence.
- [ ] Explicitly pass every important CMake option in CI.

### P4.5 Remove tracked binary release artifacts

Current tracked artifacts are under:

- `/Users/jinceyang/Desktop/codebase/merge/libpgo/dist/macosx15_arm/`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/dist/ubuntu22.04/`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/dist/ubuntu24.04/`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/dist/win11/`

Checklist:

- [ ] Remove all tracked 0.0.2 and 0.0.3 wheel files.
- [ ] Ignore `dist/`.
- [ ] Ignore `wheelhouse/`.
- [ ] Ignore standard package build output.
- [ ] Replace the per-wheel `.gitignore` entries with directory rules.
- [ ] Update README installation instructions to download the appropriate wheel from the recorded GitHub Actions artifact and install it in the documented Conda environment.
- [ ] Keep release wheels and evidence in GitHub Actions artifacts, not in git; do not claim PyPI or sdist availability.

### P4.6 Platform CI

Reference workflows:

- `/Users/jinceyang/Desktop/codebase/libpgo/.github/workflows/linux-ci.yml`
- `/Users/jinceyang/Desktop/codebase/libpgo/.github/workflows/macos-ci.yml`
- `/Users/jinceyang/Desktop/codebase/libpgo/.github/workflows/windows-ci.yml`

Target workflows:

- `/Users/jinceyang/Desktop/codebase/merge/libpgo/.github/workflows/linux-ci.yml`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/.github/workflows/macos-ci.yml`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/.github/workflows/windows-ci.yml`

Rules:

- Selectively migrate workflow structure, dependency-audit logic, and clean-install checks.
- Keep dependency installation, target names, test paths, CMake options, and package layout aligned with `/Users/jinceyang/Desktop/codebase/merge/libpgo`.
- Do not copy experimental branch tests, generated research assets, `pypgo_core`, or modular package assumptions.
- Keep exactly these three independent workflow files. Each contains one primary `build-test-package` job for its platform; do not add an orchestrator workflow, reusable workflow, cross-workflow dependency, or separate hosted release-gate workflow.
- Give each workflow the same simple triggers (`push` for the release/main branches, `pull_request`, and `workflow_dispatch`), read-only contents permission, concurrency cancellation, a timeout, and one pinned runner.
- Use `conda-incubator/setup-miniconda@v3` with a pinned Miniforge release, `conda-forge`, strict channel priority, and no default channels on Linux, macOS, and Windows.
- Do not reuse the experimental repository's full `environment.yml`. Keep each platform's minimal Conda build/runtime package lists directly in its workflow so the complete CI contract is visible in one file.
- Keep hosted CI bounded: run unit tests, config validation, package checks, and at most a tiny representative simulation smoke. Do not run the long tet/cubic sampled/IPC plus shell IPC matrix in hosted CI.
- GitHub CI owns only build, automated tests, wheel packaging, dependency/tag audit, clean-install verification, and artifact upload.
- Hosted CI success does not replace the owner-run Linux server matrix gate.
- Do not add the Linux server matrix as a GitHub Actions job. The repository provides only the driver/config/report tooling needed for the owner to run it manually.

#### Three-workflow structure and Conda environment contract

Each of the three workflows is self-contained and follows the same step order:

```text
checkout
  -> set up compiler and pinned Miniforge
  -> create release build environment
  -> verify selected BLAS/MKL/TBB provider
  -> configure/build
  -> generate bounded test assets
  -> CTest and Python smoke
  -> build wheel directly from the tested clean checkout
  -> audit dependencies/tags without repairing or vendoring Conda runtimes
  -> create a separate clean Conda runtime environment
  -> install the wheel and run linkage/API smoke outside the checkout
  -> upload artifact and dependency evidence
```

The inline build package lists contain only 0.0.4 requirements such as Python 3.12, CMake, Ninja, GMP, MPFR, Imath, fmt, zlib, Eigen, Boost, TBB, NumPy, pip, setuptools, wheel/build, test tools, and the platform provider. GMP, MPFR, Imath, and fmt must come from Conda on Ubuntu, macOS, and Windows; release CI must not install them through `apt`, Homebrew, or a system-wide package manager.

- Linux/Windows add `mkl-devel` and MKL-selected `libblas`/`liblapack`.
- macOS adds `_newaccelerate` BLAS/LAPACK packages for Python/NumPy and no MKL package; the libpgo/SuiteSparse CMake build itself links the system Accelerate framework.

PyVista, Torch, research tests, notebook dependencies, and experimental package modules are not part of release CI.

Every platform job creates a second environment under the runner's temporary directory. It must not inherit the build environment:

- Linux/Windows runtime: Python 3.12, pip, NumPy, GMP, MPFR, Imath, fmt, TBB, MKL, and MKL-selected BLAS/LAPACK.
- macOS runtime: arm64 Python 3.12, pip, NumPy, GMP, MPFR, Imath, fmt, TBB, and `_newaccelerate` BLAS/LAPACK packages; no MKL. The installed libpgo extension must continue to resolve BLAS/LAPACK through the system Accelerate framework.

The workflow must set `MKL_THREADING_LAYER=TBB` globally for Linux/Windows build and runtime checks. `PGO_CHECK_CONDA=ON` is required on all platforms; Linux/Windows additionally use `PGO_USE_MKL=ON`, `PGO_MKL_LINK_DYNAMIC=ON`, and `PGO_HAS_ORIG_PARDISO=OFF`, while macOS uses `PGO_USE_MKL=OFF`, `BLA_VENDOR=Apple`, and the LP64 BLAS/LAPACK interface.

Before compiling, record and verify:

- active Python and Conda prefix;
- resolved Conda channels and package list;
- CMake/compiler versions;
- MKL/TBB/BLAS provider and versions;
- Linux `MKLROOT=$CONDA_PREFIX`;
- Windows `MKLROOT=%CONDA_PREFIX%\Library`;
- the resolved `MKL_DIR` and `TBB_DIR`/CMake target locations.

Port the structure of `/Users/jinceyang/Desktop/codebase/libpgo/tests/check_mkl_tbb_runtime.py`, but adapt it to import the 0.0.4 `pypgo` extension rather than the experimental `pypgo._core` layout. The runtime checker must execute a real MKL-backed operation before inspecting loaded libraries.

To keep the workflows independent, each workflow builds its wheel directly from its own exact, clean checked-out commit after native and Python tests pass. Each workflow uploads exactly one platform wheel plus its commit, checksum, dependency, linkage, test, and environment evidence. No workflow builds or uploads an sdist.

#### Linux

- [ ] In `linux-ci.yml`, create the minimal build environment from the inline Linux Conda package list.
- [ ] Verify Conda supplied GMP, MPFR, Imath, fmt, `mkl-devel`, `tbb-devel`, MKL BLAS/LAPACK, and `MKL_THREADING_LAYER=TBB` before configuration.
- [ ] Configure and build native Release with tests.
- [ ] Run CTest.
- [ ] Build the release wheel directly from the tested clean checkout using `pypgo-wheel-mkl`.
- [ ] Require the intentional `linux_x86_64` platform tag and document that the wheel is supported only in the specified Conda runtime, not as a generic manylinux wheel.
- [ ] Inspect required GLIBC/GLIBCXX symbol versions against the declared compatibility baseline.
- [ ] Run `auditwheel show` for dependency evidence only; do not run `auditwheel repair`. Inspect the extension extracted from the final wheel with `readelf -d`, `patchelf --print-rpath`, and `ldd`.
- [ ] Allow only system ABI libraries and approved external Conda runtime libraries: GMP/GMPXX, MPFR, Imath, fmt, TBB, and the approved MKL/BLAS/LAPACK runtime.
- [ ] Fail if wheel contents, ELF `NEEDED`, `RPATH`, or `RUNPATH` records a build Conda prefix, source/build-tree path, `/usr/local`, `/opt/homebrew`, or another undeclared host path. A required runtime rpath must be `$ORIGIN`-relative, never an absolute build-environment path.
- [ ] Inspect the installed extension and require `mkl_tbb_thread` plus TBB while rejecting `mkl_intel_thread`, `mkl_gnu_thread`, `libiomp5`, or another non-TBB threading runtime.
- [ ] Install in a fresh Conda GMP/MPFR/Imath/fmt/MKL/TBB runtime environment outside the checkout; it must not inherit build-environment library paths.
- [ ] Verify version and run package/API smoke.

#### macOS

- [ ] In `macos-ci.yml`, create the minimal arm64 build environment from the inline macOS Conda package list.
- [ ] Verify Conda supplied GMP, MPFR, Imath, fmt, TBB, and `_newaccelerate` BLAS/LAPACK for Python/NumPy, and verify that no MKL package or linkage is selected.
- [ ] Fix or explicitly configure macOS Conda TBB discovery so CMake resolves `$CONDA_PREFIX/lib/cmake/TBB`.
- [ ] Add the platform CMake policy that sets `BLA_VENDOR=Apple` and the LP64 interface before any SuiteSparse BLAS/LAPACK lookup.
- [ ] Remove the target repository's unconditional clearing of `BLA_VENDOR` and make SuiteSparse reuse the top-level `BLAS::BLAS` and `LAPACK::LAPACK` provider.
- [ ] Fail configure unless the resolved libpgo/SuiteSparse BLAS and LAPACK link interface contains the system Accelerate framework and excludes OpenBLAS, generic Conda BLAS libraries, and MKL.
- [ ] Configure and build native Release with tests.
- [ ] Run CTest.
- [ ] Assert the runner and built wheel architecture are `arm64`; select an arm64 runner or documented cross-build path rather than trusting the runner label.
- [ ] Set and verify an intentional `CMAKE_OSX_DEPLOYMENT_TARGET`; do not inherit a host-only tag such as `macosx_26_0`.
- [ ] Build the release wheel directly from the tested clean checkout using `pypgo-wheel-no-mkl`.
- [ ] Audit with `otool` and `delocate-listdeps`; do not run a repair step that vendors Conda runtimes.
- [ ] Run `otool -L`, `otool -l`, and `delocate-listdeps` on the extension extracted from the final wheel, not only on the build-tree extension.
- [ ] Allow only macOS system libraries/frameworks (including Accelerate), the standard C++ runtime, and `@rpath` names for approved Conda runtime libraries: GMP/GMPXX, MPFR, Imath, fmt, and TBB.
- [ ] Fail if the wheel contains an external Conda runtime or an install name/`LC_RPATH` containing `/opt/homebrew`, `/usr/local`, the build Conda prefix, a source/build-tree path, or another absolute non-system host path. A required runtime rpath must be loader-relative (for example `@loader_path/...`), never an absolute build-environment path.
- [ ] Inspect the linked libraries and final wheel tag with `otool`/`delocate-listdeps`.
- [ ] Assert that the extension has no MKL linkage.
- [ ] Install in a fresh arm64 Conda GMP/MPFR/Imath/fmt/TBB/non-MKL runtime environment outside the checkout. It must not receive Homebrew paths, manually copied libraries, or the build environment.
- [ ] Verify version and run package/API smoke.

#### Windows

- [ ] In `windows-ci.yml`, create the minimal build environment from the inline Windows Conda package list.
- [ ] Verify Conda supplied GMP, MPFR, Imath, fmt, `mkl-devel`, `tbb-devel`, MKL BLAS/LAPACK, and `MKL_THREADING_LAYER=TBB` before configuration.
- [ ] Verify CMake resolves MKL/TBB from `%CONDA_PREFIX%\Library`, not from a system installation.
- [ ] Configure and build native Release with tests under MSVC.
- [ ] Run CTest.
- [ ] Build the release wheel directly from the tested clean checkout using `pypgo-wheel-mkl`.
- [ ] Audit with `dumpbin` and `delvewheel show`; do not run a repair step that vendors Conda runtimes.
- [ ] Fail if the wheel contains Conda-provided GMP/GMPXX, MPFR, Imath, fmt, MKL, BLAS/LAPACK, or TBB DLLs.
- [ ] Fail if `dumpbin /DEPENDENTS`, `delvewheel show`, an embedded manifest, package metadata, or loader-search configuration refers to the build Conda prefix, source/build tree, a developer package-manager directory, or another undeclared absolute host path.
- [ ] Inspect dependencies with `dumpbin`/`delvewheel` and require the TBB MKL threading layer while rejecting the Intel OpenMP threading layer.
- [ ] Install in a fresh Conda GMP/MPFR/Imath/fmt/MKL/TBB runtime environment outside the checkout. It must not inherit build-environment `PATH` entries.
- [ ] Verify version and run package/API smoke.

#### Common package verification

Every wheel job must run the following inside its separate supported Conda runtime environment:

```bash
python -m pip install --no-deps <wheel>
python -c "import pypgo; assert pypgo.__version__ == '0.0.4'"
python -m pip check
pytest -q <installed-package smoke tests>
```

Also:

- [ ] The clean environment is created from a minimal runtime package list, never cloned from the build environment.
- [ ] On Linux/Windows, set `MKL_THREADING_LAYER=TBB`, execute an MKL-backed operation, require `mkl_tbb_thread` and TBB in loaded-library evidence, and reject `mkl_intel_thread`, `mkl_gnu_thread`, and `libiomp5`.
- [ ] On macOS, execute representative BLAS GEMM and LAPACK solve operations; verify the installed extension links system Accelerate, the clean Python environment loads the `_newaccelerate` shim/system Accelerate stack, and neither side loads MKL or OpenBLAS.
- [ ] Verify wheel contents do not include the external Conda GMP/MPFR/Imath/fmt/MKL/TBB runtime set.
- [ ] Run imports from a temporary directory outside `/Users/jinceyang/Desktop/codebase/merge/libpgo` so source-tree imports cannot mask a broken wheel.
- [ ] List wheel contents and reject forbidden paths.
- [ ] Upload exactly one final tested wheel per platform job.
- [ ] Retain source commit, wheel SHA-256, build logs, test results, and dependency/version/linkage reports in the same CI artifact bundle.
- [ ] Configure `upload-artifact` with an explicit repository-supported `retention-days` value and document that Actions artifacts expire; the release record must not imply permanent package hosting.
- [ ] Keep the expensive full matrix out of hosted CI and document the owner-run Linux server command in the release checklist.

## 9. Implementation phase 5 — tests

### P5.1 Native unit tests

- [ ] Cubic energy/gradient/Hessian finite differences.
- [ ] Cubic invalid rest-element validation.
- [ ] Cubic mesher deterministic topology/basic geometry.
- [ ] Cubic Python generator argument, no-overwrite, manifest, and repeatability tests.
- [ ] Tet material max-step.
- [ ] Cubic material max-step.
- [ ] IPC CCD max-step.
- [ ] IPC barrier energy/gradient/Hessian finite differences.
- [ ] IPC core validation matrix.
- [ ] Mixed fixed/dynamic Hessian assembly.
- [ ] Newton fixed-pattern reuse.
- [ ] Newton dynamic-pattern safety.
- [ ] Generic volume loader equivalence.
- [ ] Config dispatch and fallback.
- [ ] Missing/unreadable config returns failure through CLI, C, and Python.
- [ ] Dynamic-output-directory, static-output-file, no-overwrite, and no-recursive-delete contracts.
- [ ] `runDynamicSim --resume` validation, manifest matching, atomic checkpoint/frame installation, atomic manifest update, and no-overwrite behavior.
- [ ] `runStaticSim` atomic output-file behavior and rejection of `--resume`.
- [ ] C ABI runner matrix.
- [ ] C API Hessian export preserves finite fractional `double` values.
- [ ] Shell IPC regression tests through the shared dispatcher; shell + sampled is rejected and is not part of the release matrix.

### P5.2 Linux server API smoke matrix

Use small dedicated test configs and assets under `/Users/jinceyang/Desktop/codebase/merge/libpgo/tests/`, not large showcase runs.

These five short dynamic API cases run on the Linux server after implementation. Hosted CI only needs parser/unit coverage and a bounded representative simulation smoke.

| Case | C API | Python API | CLI | Expected |
| --- | --- | --- | --- | --- |
| tet + sampled | Required | Required | Required | Setup, one or more solves, finite output, return 0 |
| tet + IPC | Required | Required | Required | Dynamic Hessian + CCD path executed, return 0 |
| cubic + sampled | Required | Required | Required | Cubic FEM path executed, return 0 |
| cubic + IPC | Required | Required | Required | Cubic + dynamic IPC Hessian + CCD path executed, return 0 |
| shell + IPC | Required | Required | Required | Shell IPC path executes dynamic Hessian/CCD as applicable, return 0 |
| missing contact model | Required | Required | Required | Same behavior as tet/cubic sampled config |
| bad contact model | Required | Required | Required | Useful error, nonzero result |
| all five static domain/contact cases | Required | Required | `runStaticSim` | Shared sampled/IPC backend executes and writes one new result file atomically, return 0 |
| sampled `--resume` | Not supported | Not supported | `runDynamicSim` | Continues after latest valid checkpoint without modifying earlier output |
| IPC `--resume` | Not supported | Not supported | `runDynamicSim` | Continues after latest valid checkpoint without modifying earlier output |
| static `--resume` | Not supported | Not supported | Required failure | `runStaticSim` rejects the flag before output mutation |

### P5.3 Compatibility tests

- [ ] Selected main-era tet/sample scenarios are migrated to the new `examples/configs/` and `examples/assets/` paths and preserve their simulation semantics.
- [ ] Old example filesystem paths are absent and are not treated as compatibility entry points.
- [ ] Run one canonical 0.0.3 tet + sampled case in the read-only main checkout and the final 0.0.4 implementation, then compare selected displacement/energy/output observables within recorded tolerances.
- [ ] Existing public C symbols still export.
- [ ] Diff the Linux/macOS/Windows exported C symbol lists against the 0.0.3 reference and reject an unapproved removal or incompatible change.
- [ ] Existing Python smoke tests from `/Users/jinceyang/Desktop/codebase/main/libpgo` that apply to 0.0.4 still pass.
- [ ] `const char *` cleanup does not change binary calling convention.
- [ ] Document the intentional CLI migration from `runSim`/`runIPCSim`/`runShellSim` to `runDynamicSim` and `runStaticSim`; legacy executable-name compatibility is not required.
- [ ] Run `cmake --install` into a temporary prefix.
- [ ] Confirm the install contains both `pgo_c.h` and `pgo_c_def.h`.
- [ ] Compile and run a minimal external pure-C consumer against the installed headers and shared library.
- [ ] Verify `find_package(pgo 0.0.4 CONFIG REQUIRED)` and the installed CMake export work without paths from the source or build tree.

### P5.4 Sanitizer and numerical checks

- [ ] Run at least one Linux ASan/UBSan native test job if compatible with dependencies.
- [ ] Check every matrix smoke output for NaN and infinity.
- [ ] Check sparse Hessians for correct dimensions and finite entries.
- [ ] Verify no cubic rest determinant is zero or negative.
- [ ] Run a clean Release build after sanitizer/debug validation.

### P5.5 Owner-run Linux server validation with MKL Pardiso

This is the release-blocking home for expensive simulation testing. The owner runs and confirms it manually on the Linux server. It is intentionally outside hosted CI and is not represented as a GitHub required check.

#### Build contract

- [ ] Add a dedicated configure/build preset, preferably `linux-server-mkl-pardiso`.
- [ ] Configure `CMAKE_BUILD_TYPE=Release`.
- [ ] Configure `PGO_USE_MKL=ON`.
- [ ] Configure `PGO_HAS_ORIG_PARDISO=OFF`.
- [ ] Configure and record `MKL_THREADING=tbb_thread`.
- [ ] Enable the native tests and targets required by the runner, C API, and pypgo.
- [ ] Use the normal Newton solver path.
- [ ] Confirm the build defines `PGO_HAS_MKL`.
- [ ] Confirm `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/core/nonlinearOptimization/NewtonSolver.cpp` therefore selects `EigenMKLPardisoSupport`.
- [ ] Confirm `PGO_HAS_ORIG_PARDISO` is absent so the proprietary/original Pardiso branch cannot shadow MKL Pardiso.
- [ ] Capture `CMakeCache.txt`, configure output, compiler version, MKL version, TBB version, Python version, and `ldd` output for the runner and Python extension.
- [ ] Require runtime dependency evidence for `mkl_tbb_thread` and TBB, and reject the Intel OpenMP threading layer.
- [ ] Set and record a fixed MKL thread policy for the run, including `MKL_DYNAMIC` and `MKL_NUM_THREADS`; do not compare runs made with undocumented thread settings.
- [ ] Set `MKL_THREADING_LAYER=TBB` in the server test environment and record it in the report.

The name “MKL Pardiso” in this section means the `PGO_HAS_MKL && !PGO_HAS_ORIG_PARDISO` branch, not `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/core/eigenSupport/EigenOrigPardisoSupport.*`.

#### Reproducible manual server runner

- [ ] Add a release-test driver under `/Users/jinceyang/Desktop/codebase/merge/libpgo/tests/release/`.
- [ ] Keep the driver manually invocable; do not add a workflow that runs the server matrix.
- [ ] Make it accept the build directory, a new report/output root, timeout, generated-cubic asset directory, and selected cases as arguments.
- [ ] Run cases serially by default so one failure is attributable and output directories do not collide.
- [ ] Record the exact git commit and reject a dirty worktree unless explicitly overridden for diagnosis.
- [ ] Run the checked-in cubic generation script first and hash or record the exact configs, source assets, generator arguments, manifest, and generated meshes used.
- [ ] Write all simulation output outside `/Users/jinceyang/Desktop/codebase/merge/libpgo`.
- [ ] Materialize configs into the report directory with absolute generated-mesh/output paths and hash the exact materialized JSON used.
- [ ] Give every dynamic case a unique nonexistent output directory and every static case a unique nonexistent output file; verify that rerunning without a new path fails safely.
- [ ] Produce a machine-readable JSON or JUnit report plus a concise Markdown summary.
- [ ] Record command, return code, wall time, solver iterations where available, output file count, and failure reason per case.
- [ ] Support rerunning one failed cell without rerunning the entire matrix during diagnosis.

#### Tier 1 — short API matrix

- [ ] Run tet + sampled through C API.
- [ ] Run tet + IPC through C API.
- [ ] Run cubic + sampled through C API.
- [ ] Run cubic + IPC through C API.
- [ ] Run shell + IPC through C API.
- [ ] Repeat all five short cases through `pypgo.run_sim_from_config`.
- [ ] Repeat all five short cases through `runDynamicSim`.
- [ ] Run all five required static cases through `runStaticSim`, C, and Python.
- [ ] Keep these cases short enough to diagnose API/config integration, while still executing energy, gradient, Hessian, maximum-step, factorization, and output.

#### Tier 2 — representative numerical matrix

- [ ] Run one representative tet + sampled simulation at the agreed server duration.
- [ ] Run the comparable tet + IPC simulation.
- [ ] Run the comparable cubic + sampled simulation.
- [ ] Run the comparable cubic + IPC simulation.
- [ ] Run the comparable shell + IPC simulation.
- [ ] Use the same physical scene, scale, material, timestep policy, boundary conditions, and recorded cubic-generator parameters wherever the two contact models permit.
- [ ] Where sampled external-object contact and IPC floor contact differ by design, record that difference instead of claiming numerically identical contact geometry.
- [ ] Require zero process failures, zero uncaught exceptions, and no NaN or infinity.
- [ ] Require solver convergence under the configured iteration/tolerance policy.
- [ ] Require no invalid or inverted rest element.
- [ ] Check output frame count and basic displacement bounds.
- [ ] Treat timings as diagnostic information, not a 0.0.4 performance-regression promise.

#### Evidence and rerun rules

- [ ] Tie the server report to the exact release candidate commit.
- [ ] Store the preset name, CMake cache, dependency versions, environment variables, commands, logs, and report together.
- [ ] If code, configs, cubic generator code/parameters, generated-mesh hashes, solver logic, IPC logic, or dependency revisions change after the run, invalidate the affected results.
- [ ] After a correctness fix, rerun the affected cell and then rerun the complete five-case Tier 2 matrix before release.
- [ ] Do not mark Gate B complete based only on local or hosted-CI smoke tests.
- [ ] Attach or link the owner-confirmed server report in the final release record.

## 10. Implementation phase 6 — documentation and release notes

### P6.1 Docs repository branch and submodule

Canonical docs repository and branch:

- `https://github.com/annajcy/libpgo-doc.git`
- Branch `v0.0.4`
- Frozen starting commit `ef0834c67a280d16836579942246129087c99e44`
- Submodule checkout `/Users/jinceyang/Desktop/codebase/merge/libpgo/docs`

Checklist:

- [ ] Initialize/update the superproject's `docs` submodule, then check out the canonical `v0.0.4` branch from the frozen starting commit inside it; do not create a separate docs worktree or substitute a `release/0.0.4` docs branch.
- [ ] Preserve the docs-site design, theme/layout, navigation machinery, and reusable site assets, but simplify and repair the build/deployment workflow where the current experimental assumptions do not match 0.0.4.
- [ ] Treat the starting documentation content as experimental and nonauthoritative.
- [ ] Derive the 0.0.4 content and changelog from the code difference between `/Users/jinceyang/Desktop/codebase/main/libpgo@ed9675d56403b73a0cf32084e5da9b0011f88e2c` and the final 0.0.4 code.
- [ ] Do not use the empty docs `main` branch as a content or changelog reference.
- [ ] Delete or rewrite pages, sections, examples, API descriptions, and navigation entries that do not exist in the actual 0.0.4 code.
- [ ] Remove references to generalized Neo-Hookean, tricubic Hermite, experimental package architecture, ARPACK restoration, and other unreleased features.
- [ ] Replace the workflow's checkout of the libpgo default branch with the exact authorized `v0.0.4` tag/commit.
- [ ] Remove the workflow's dependency on experimental-only `libpgo-src/environment.yml` and the nonexistent `pypgo-ci` preset; use the actual 0.0.4 docs/source build path and only the dependencies needed to render the site.
- [ ] If GitHub requires the Pages workflow to exist on the default docs branch, keep only the minimal deployment workflow on docs `main` and make it check out the `v0.0.4` content branch explicitly; do not use empty `main` as a content source.
- [ ] Keep Pages permissions only in the docs deployment workflow; the three libpgo platform CI workflows remain read-only.
- [ ] Use a public HTTPS submodule URL in `/Users/jinceyang/Desktop/codebase/merge/libpgo/.gitmodules`.
- [ ] Optionally record `branch = v0.0.4`, but pin the code repository gitlink to the exact final docs commit.
- [ ] Commit docs content changes in the docs submodule first; after that exact docs commit exists, update and commit the superproject gitlink.
- [ ] Make CI checkout submodules recursively.
- [ ] Build the exact `v0.0.4` docs content in CI and deploy only after the local/CI link and snippet checks pass.
- [ ] Check internal links and code snippets.

### P6.2 Detailed changelog structure

Create detailed pages covering:

- [ ] 0.0.4 overview.
- [ ] Cubic mesh and cubic FEM.
- [ ] Cubic mesher and the Python `.veg` generation script; state that generated cubic `.veg` files are not stored in git.
- [ ] IPC core, barrier, self-contact, floor contact, and CCD.
- [ ] Hessian topology contract.
- [ ] Newton symbolic-factorization fix.
- [ ] Material and contact maximum-step behavior.
- [ ] Shared config dispatcher and the `runDynamicSim`/`runStaticSim` split.
- [ ] Dynamic-output-directory safety, static-output-file safety, and `runDynamicSim --resume` behavior.
- [ ] C API.
- [ ] Python API.
- [ ] Five-case config matrix: tet/cubic × sampled/IPC plus shell + IPC.
- [ ] Build presets.
- [ ] Wheel installation from GitHub Actions artifacts; explicitly state that 0.0.4 has no sdist and is not on PyPI.
- [ ] Linux/Windows MKL Pardiso wheel policy, TBB threading layer, macOS no-MKL policy, and Linux server validation environment/matrix summary.
- [ ] Dependency versions and optional features.
- [ ] Migration from 0.0.3.
- [ ] Old-to-new example config and asset path map, explicitly stating that old paths are not supported.
- [ ] Known limitations.

Known limitations must say:

- [ ] IPC external objects are not supported in 0.0.4 unless implementation and tests prove otherwise.
- [ ] IPC supports the implemented self-contact and explicit floor path.
- [ ] CUDA utilities are retained but not part of default wheels.
- [ ] Native Pardiso is retained but not bundled.
- [ ] Guaranteed wheel platforms/Python version are listed explicitly.
- [ ] Resume is CLI-only in 0.0.4; the frozen C and Python config-runner entry points remain new-run-only.

### P6.3 Main repository release notes

Files:

- `/Users/jinceyang/Desktop/codebase/merge/libpgo/CHANGELOG.md`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/README.md`

Checklist:

- [ ] Add a concise `0.0.4` section.
- [ ] Link to the detailed docs changelog.
- [ ] List cubic FEM/mesh/mesher.
- [ ] List IPC and CCD maximum-step.
- [ ] List the five-case unified config matrix, the volume missing-`contact-model` default, and shell's explicit IPC requirement.
- [ ] List C/Python unified entry.
- [ ] List the `runDynamicSim`/`runStaticSim` split, no-delete output policies, and dynamic-only `--resume` workflow.
- [ ] List correctness fixes.
- [ ] List packaging changes.
- [ ] List the config/asset-separated example layout and the intentional example-path migration.
- [ ] List known limitations.
- [ ] Remove instructions that point users to wheel files committed under `dist/`.

## 11. Implementation phase 7 — release-candidate, tag, and CI artifact flow

### RC1 — Local and PR readiness

- [ ] Worktree is clean except intentional release changes.
- [ ] Diff contains no benchmark/research/package-architecture migration.
- [ ] Diff contains no generated cache, output, old wheel, or local absolute path outside this internal plan.
- [ ] All release checklist items are linked to commits or CI results.
- [ ] Native and Python tests pass locally on at least one platform.
- [ ] Documentation build passes.

### RC2 — `upstream-0.0.4` integration candidate

- [ ] Merge the completed working branch into `upstream-0.0.4`.
- [ ] Record the exact integration merge commit and require a clean worktree.
- [ ] All platform native jobs pass.
- [ ] Hosted CI’s bounded unit/config/package smoke passes.
- [ ] The owner-confirmed Linux server Tier 1 API matrix passes.
- [ ] The owner-confirmed Linux server Tier 2 five-case dynamic numerical matrix passes with MKL Pardiso.
- [ ] The Linux server report identifies the exact `upstream-0.0.4` integration commit and clean build configuration.
- [ ] All three wheel-only workflows pass and upload one tested wheel each.
- [ ] Wheels pass validation in separate clean supported Conda runtime environments without vendored Conda libraries.
- [ ] Artifact contents are audited.
- [ ] Dependency versions and build tool versions are captured.
- [ ] Cubic generator manifest and generated-mesh validation report are attached.
- [ ] Prepare a concise owner handoff with the code diff, docs diff, gate evidence, artifacts, checksums, and known limitations.
- [ ] Stop for owner review. Do not create or merge the `upstream-0.0.4` to `main` pull request without explicit authorization.

### RC3 — Owner-authorized main integration and final artifact rerun

- [ ] Create or merge the `upstream-0.0.4` to `main` pull request only after owner approval.
- [ ] If the main merge changes the commit, bind the tag and final provenance to the exact resulting main commit and rerun the affected final checks.
- [ ] Run the three platform workflows on the exact final commit.
- [ ] Download each wheel/evidence bundle by workflow run and artifact ID.
- [ ] Verify every recorded SHA-256 and repeat clean-Conda version/import/API smoke tests against the downloaded wheel files.
- [ ] Record the exact final workflow URLs, run IDs, artifact IDs, source commit, and checksums in the release record.

### Final tag and artifact delivery

- [ ] Tag the exact tested commit as `v0.0.4`.
- [ ] Verify the tag includes the exact docs submodule commit.
- [ ] Ensure all three platform workflows have completed on the tagged commit and retained their wheel/evidence artifacts.
- [ ] Document how to download each platform wheel from its recorded GitHub Actions artifact and install the local file in the supported Conda environment.
- [ ] Do not build or publish an sdist and do not upload artifacts to TestPyPI or PyPI.
- [ ] Verify docs deployment.
- [ ] Announce the release only after the tag, all three workflow artifact bundles, checksums, installation verification, and docs are available.

### Delivered-artifact incident handling

If a serious defect is found after artifact delivery:

- [ ] Mark affected workflow artifacts and deployed docs with a concise warning in the release record.
- [ ] Disable or remove public links to known-bad artifacts where GitHub permits it while retaining their hashes and incident evidence.
- [ ] Fix forward in `0.0.5`; never silently substitute a different wheel for a recorded `0.0.4` artifact ID/checksum.

### Post-release

- [ ] Create a `0.0.4` maintenance tag/branch policy if patches are expected.
- [ ] Record any non-blocking follow-ups separately; do not silently fold new features into the release.
- [ ] Confirm the release record identifies one immutable-by-checksum wheel set and does not silently redirect to a later rebuild.
- [ ] Close the release checklist with final commit IDs and artifact checksums.

## 12. Recommended implementation/commit order

Use small, independently reviewable commits in this order:

1. `chore(release): establish 0.0.4 baseline and remove tracked wheels`
2. `fix(ipc): remove legacy CIPC wrapper and validate core inputs`
3. `fix(optimization): define dynamic Hessian contract`
4. `fix(solver): reuse fixed-topology symbolic factorization`
5. `refactor(energy): simplify feasible max-step interface`
6. `refactor(cubic): harden cubic deformation model ownership and rest validation`
7. `refactor(mesh): unify volumetric SimulationMesh conversion`
8. `chore(profiling): remove scoped profiling infrastructure`
9. `refactor(runner): add dynamic/static CLIs and unified domain/contact dispatcher`
10. `feat(api): align C and Python config-runner entry`
11. `refactor(examples): separate configs, deduplicate assets, and add cubic veg generator`
12. `test(release): add dynamic/static five-case domain-contact coverage`
13. `test(server): add reproducible Linux MKL Pardiso matrix driver`
14. `build(deps): pin release dependency revisions`
15. `chore(release): finalize version and wheel metadata`
16. `build(pypgo): add explicit MKL/no-MKL wheel presets`
17. `ci(release): build, audit, and isolate-test platform wheels`
18. `docs(release): rewrite v0.0.4 docs content and repair deployment`
19. `docs(release): finalize changelog and release notes`

After each correctness commit:

- [ ] Build affected native targets.
- [ ] Run focused tests.
- [ ] Run the full native suite before moving to packaging.
- [ ] Do not combine asset churn with solver or API correctness changes.

## 13. Final definition of done

0.0.4 is done only when:

- The code release commit is based on `/Users/jinceyang/Desktop/codebase/merge/libpgo`, not the experimental branch.
- Compatibility has been checked against `/Users/jinceyang/Desktop/codebase/main/libpgo`.
- Only explicitly selected reference work has been migrated from `/Users/jinceyang/Desktop/codebase/libpgo`.
- Cubic FEM, cubic mesher, IPC, sampled contact, and maximum-step behavior are correct.
- `runDynamicSim`, `runStaticSim`, C, and Python use one config parser and run tet/cubic × sampled/IPC plus shell + IPC through shared implementations; the old `runIPCSim` implementation has been migrated rather than retained in parallel.
- New dynamic runs create only new output directories; static runs atomically create only new result files; `runDynamicSim --resume` validates and continues an existing dynamic run without deleting or overwriting prior output.
- Existing shell assets and the shell IPC path are retained and regression-tested; shell + sampled is not part of the 0.0.4 contract.
- Example configs and static assets are separated under `examples/configs/` and `examples/assets/`, with no legacy-path compatibility copies or unintended duplicate assets.
- The obsolete CIPC wrapper and scoped profiling are gone.
- Dynamic Hessian and solver symbolic-analysis lifecycles are correct.
- A checked-in Python script reproducibly generates validated cubic `.veg` files into an untracked output directory; generated cubic `.veg` files are not stored in git.
- Optional backward/CUDA/native-Pardiso code remains, with release dependencies pinned and optional features off by default where required.
- Linux and Windows wheels use MKL Pardiso with the TBB threading layer; the macOS wheel does not use MKL and links libpgo/SuiteSparse to the system Accelerate framework; all external GMP/MPFR/Imath/fmt/TBB/MKL runtimes come from the documented Conda environment.
- Exactly one wheel per platform is built, audited, isolated-tested, checksummed, and uploaded with evidence as a GitHub Actions artifact; no sdist or PyPI release is produced.
- The exact `upstream-0.0.4` integration commit passes the owner-confirmed Linux server matrix using MKL Pardiso with TBB.
- The docs submodule is pinned to the exact final commit of `annajcy/libpgo-doc` branch `v0.0.4`.
- The exact workflow run IDs, artifact IDs, wheels, and checksums validated for the tagged commit are recorded as the 0.0.4 artifact set.
