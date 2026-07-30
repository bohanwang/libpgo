# libpgo 0.0.4 Implementation Plan

The short operational gate list and evidence record live in
[`RELEASE_0.0.4_CHECKLIST.md`](RELEASE_0.0.4_CHECKLIST.md). This document contains the detailed implementation contract.

| Field | Value |
| --- | --- |
| Status | Phases 1–4 implemented; Phase 5 validation is next; Phase 6 content is in progress |
| Target version | `0.0.4` |
| Target repository | `/Users/jinceyang/Desktop/codebase/merge/libpgo` |
| Target baseline branch | `upstream-temp-0.0.4` |
| Target baseline commit | `1670f6258d941d80cd9d718fc033fccc193d5e20` |

This is an internal release document. The absolute paths are intentional: they identify the three local code snapshots used for the 0.0.4 audit and prevent “current branch”, “temp”, and “main” from being confused during implementation.

## Current execution path

This section is the authoritative operational summary. Any later section
explicitly labelled historical is retained only to explain earlier decisions
and must not be used as implementation guidance.

| Phase | State | Current result / next gate |
| --- | --- | --- |
| Phase 1 — correctness | Complete | IPC, Hessian, Newton lifetime, maximum-step, cubic, and validation work is covered by the native suite. |
| Phase 2 — reusable runners | Complete | CLI implementations live in `simulationRunner`; C and Python reuse the same sampled/IPC paths without a backend hierarchy. |
| Phase 3 — examples/assets | Complete | Configs/assets are separated and lightweight tet/cubic assets are generated for bounded smoke tests. |
| Phase 4 — build/package | Implemented locally | No Conda; `uv` drives one locked environment; repaired standalone wheels are produced by three platform workflows. macOS local build, 154 native tests, build-tree Python smoke, wheel build, and docs build pass. |
| Phase 5 — validation | **Next** | Run hosted Linux/macOS/Windows `build-test` and `package` jobs, close compatibility/sanitizer gaps, then run the owner-controlled Linux MKL Pardiso matrix. |
| Phase 6 — docs/release notes | **In progress alongside Phase 5** | Handwritten 0.0.4 content, repository release notes, limitations, dependency table, migration guide, and repaired Pages workflow are drafted. Final CI/server evidence and the docs gitlink remain pending. |
| Phase 7 — release candidate | Requires explicit authorization | Integrate, rerun exact-commit workflows, record artifacts/checksums, then stop for owner review before tag/delivery. |

Current build and packaging rules:

- `pyproject.toml` is the single project-version source and the single pinned
  `uv` version source.
- `uv sync --locked` prepares `.venv`; it does not install pypgo.
- Developer builds use `uv run cmake --preset pypgo-wheel`, whose preset owns
  Ninja, `build/pypgo`, tests, and portable-wheel defaults.
- Developers import the extension directly from
  `build/pypgo/src/python/pypgo` through `PYTHONPATH`; ordinary C++ changes
  require only an incremental CMake build.
- Wheel packaging is a separate path: `uv build --wheel
  --no-build-isolation`, platform repair/audit, and installation into a fresh
  ordinary virtual environment.
- Linux and Windows default to MKL with its matching TBB runtime. macOS
  forcibly disables MKL and uses Accelerate.
- Linux/macOS acquire GMP and MPFR from the platform package manager for the
  build and bundle them during wheel repair. Windows uses the approved
  repository prebuilt GMP/MPFR files and bundles them with `delvewheel`.
- Each platform workflow has two jobs: `build-test` uploads native evidence;
  a fresh dependent `package` job downloads that evidence, builds and repairs
  the wheel, performs clean-install tests, finalizes provenance, and uploads
  exactly one wheel plus evidence.

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
- [ ] Derive the documentation release notes from the code difference between the 0.0.3 compatibility baseline and final 0.0.4 code; the empty docs `main` branch is not a content reference.
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
- [ ] Preserve the existing dynamic config-runner coverage:

  | Simulation domain | Contact model | Required |
  | --- | --- | --- |
  | tet | sampled | Yes |
  | tet | IPC | Yes |
  | cubic | sampled | Yes |
  | cubic | IPC | Yes |
  | shell | sampled | Yes, through `runShellSim` |
  | shell | IPC | Yes |

- [ ] Preserve both existing shell runners and migrate their configs/assets to the new separated example layout.
- [ ] Do not add IPC external-object contact if it is not already implemented. The current `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/tools/runSim/runIPCSimSetup.cpp` rejects `external-objects`; document IPC 0.0.4 as supporting self-contact plus the configured floor, not arbitrary external objects.
- [ ] Preserve source compatibility with the 0.0.3 APIs wherever the agreed fixes do not require a change.
- [ ] Do not preserve 0.0.3 example file paths. The 0.0.4 examples switch directly to the new config/asset-separated layout without compatibility JSON copies.

### 2.2 Public API contract

- [ ] Keep one config-driven C/Python entry point that dispatches sampled versus IPC.
- [ ] Preserve the established CLI split:

  ```text
  runSim <config>
  runIPCSim <config>
  runShellSim <config>
  ```

- [ ] Keep each executable focused on its existing sampled/IPC/domain responsibility.
- [ ] Extract the tool implementations into `simulationRunner` so C/Python reuse them without introducing a new backend hierarchy.
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
- [ ] Direct object-level Python wrappers for every new cubic/IPC class are not required for 0.0.4; the existing config runner is the required public surface.

### 2.3 Explicitly out of scope

- [ ] No new simulation algorithms beyond what exists in the 0.0.4 baseline.
- [x] Support static IPC through the existing `runIPCSim`/shared dispatcher
  without adding separate static/dynamic executables or a backend hierarchy.
- [ ] No new resume workflow or output transaction layer.
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

- Wheel matrix: CPython 3.12 with `manylinux_2_28_x86_64`,
  `macosx_26_0_arm64`, and `win_amd64` platform tags.
- Linux and Windows wheels use MKL via `PGO_USE_MKL=ON`, with `MKL_THREADING=tbb_thread`, `PGO_HAS_ORIG_PARDISO=OFF`, and runtime linkage verified to use TBB rather than Intel OpenMP.
- macOS wheels force `PGO_USE_MKL=OFF`; libpgo/SuiteSparse selects the system
  Accelerate framework.
- Ordinary source dependencies use pinned FetchContent archives. Imath is
  built statically. TBB follows the pinned MKL dependency on Linux/Windows and
  comes from Homebrew on macOS.
- Linux and macOS build against system-package GMP/MPFR; Windows uses the
  approved repository prebuilt files. Wheel repair bundles the required
  non-system runtime closure.
- Repaired wheels are standalone apart from declared Python requirements and
  platform system libraries; users do not need Conda or separately installed
  MKL, TBB, GMP, or MPFR.
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

### Gate B — API matrix

- [ ] Fast config/parser/API smoke passes in normal CI without running the long simulation matrix.
- [ ] On the Linux server, C API passes tet + sampled smoke.
- [ ] On the Linux server, C API passes tet + IPC smoke.
- [ ] On the Linux server, C API passes cubic + sampled smoke.
- [ ] On the Linux server, C API passes cubic + IPC smoke.
- [ ] On the Linux server, C API passes shell + IPC smoke.
- [ ] On the Linux server, Python API passes the same five short API cases.
- [ ] On the Linux server, representative sampled cases pass through `runSim`, IPC cases through `runIPCSim`, and sampled shell through `runShellSim`.
- [ ] The Linux server build demonstrably uses MKL Pardiso, not original/native Pardiso and not the Eigen fallback.
- [ ] Linux server and Linux/Windows wheel evidence demonstrate the MKL TBB threading layer; macOS evidence demonstrates no MKL linkage.
- [ ] Missing `contact-model` is tested and selects sampled contact.
- [ ] Invalid mesh/contact combinations fail with a useful error and a nonzero result.

### Gate C — Packaging

- [ ] Linux `manylinux_2_28_x86_64` wheel builds in the manylinux container,
  is repaired/audited with `auditwheel`, and installs in a fresh ordinary
  virtual environment with no inherited library path.
- [ ] macOS 26 arm64 raw wheel builds locally; final hosted CI must repair/audit it
  with `delocate` and install it in a fresh ordinary virtual environment.
- [ ] Windows `win_amd64` wheel builds, is repaired/audited with `delvewheel`,
  and installs in a fresh ordinary virtual environment.
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
- [ ] Detailed module release notes are complete and build successfully from `annajcy/libpgo-doc` branch `v0.0.4`; deployment happens only after explicit authorization.
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
- 0.0.3 (main) code checkout is used as the release-note baseline in docs;

### P0.2 Capture baselines before edits

- [ ] Run the existing native build and tests in `/Users/jinceyang/Desktop/codebase/merge/libpgo`.
- [ ] Run the existing `tests/pypgo` suite in `/Users/jinceyang/Desktop/codebase/merge/libpgo`.
- [ ] Record existing failures separately from regressions introduced during cleanup.
- [ ] Inventory all tracked files under `/Users/jinceyang/Desktop/codebase/merge/libpgo/dist`.
- [ ] Record current example asset sizes and element counts.
- [ ] Record the current example inputs and the cubic generator parameters needed to reproduce test meshes.

### P0.3 Resolve high-risk integration paths with bounded spikes

Before broad refactoring:

- [ ] Confirm a tiny tet + IPC dynamic smoke reaches the existing IPC energy/derivative/CCD/max-step path through the extracted runner.
- [ ] On each platform, build and repair one wheel from the clean checkout and
  prove a fresh ordinary virtual environment can load it without inherited
  build-prefix, source-tree, or host-package-manager search paths.
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

## 6. Implementation phase 2 — extract reusable simulation runners

Phase 2 deliberately preserves the existing executable split:

- `runSim` owns the sampled-penalty volume workflow.
- `runIPCSim` owns the IPC volume/shell workflow.
- `runShellSim` remains the sampled shell workflow.
- No `runDynamicSim` / `runStaticSim`, resume transaction, typed backend hierarchy, or solver architecture rewrite is part of 0.0.4.

The problem addressed here is implementation duplication between the tools and the
C/Python APIs, not the names or number of executables.

### P2.1 Move implementations into `simulationRunner`

- [x] Add `src/simulationRunner` after core libraries and before C/Python/tools.
- [x] Move the sampled and IPC simulation bodies plus their existing setup helpers into the module.
- [x] Keep `runSim.cpp` and `runIPCSim.cpp` as small CLI argument/result wrappers.
- [x] Preserve `runShellSim` and its existing sampled-shell behavior.
- [x] Avoid backend interfaces, prepared-request object graphs, and unified static/dynamic execution abstractions.
- [x] Preserve relative config-path resolution and CLI-only logging behavior.

The reusable entry points are intentionally small:

```cpp
int runSampledSimulationFromConfig(
  const std::filesystem::path &configPath,
  bool enableCliLog = false);

int runIPCSimulationFromConfig(
  const std::filesystem::path &configPath,
  bool enableCliLog = false);

int runSimulationFromConfig(const std::filesystem::path &configPath);
```

The third function is only the language-API dispatcher. Missing `contact-model`
keeps the 0.0.3 sampled default; `sampled` and `ipc` select the corresponding
existing runner, and unknown values fail.

### P2.2 C and Python reuse

- [x] Preserve the C symbol and signature `pgo_run_sim_from_config`.
- [x] Remove its duplicated tet + sampled simulation body.
- [x] Route it through `runSimulationFromConfig`.
- [x] Return nonzero for missing/unreadable configs, unknown contact models, setup failures, and simulation failures.
- [x] Prevent C++ exceptions from crossing the C ABI.
- [x] Keep Python routed through the C API and preserve its integer result contract.
- [x] Do not add a second Python-side dispatcher or simulation implementation.

### P2.3 Build and regression evidence

- [x] Link `simulationRunner` explicitly into `pgo_c`, `pgo_c_static`, `runSim`, and `runIPCSim`.
- [x] Keep argparse in executable targets rather than the reusable module.
- [x] Test missing config, sampled default routing, explicit IPC routing, and unknown-contact failure.
- [x] Run the existing sampled and IPC setup/integration suites after extraction.
- [x] Confirm the tool wrappers no longer contain solver/setup bodies.
- [x] Add focused tet/cubic/shell static IPC runner coverage using the existing
  Newton, dynamic-Hessian, material max-step, and CCD contracts.
- [x] Keep unrelated output-safety, resume, sampled-static contact redesign,
  and public CLI redesign work deferred beyond 0.0.4.

Acceptance: tools, C, and Python reach the same sampled/IPC implementations without
changing the established runner split or introducing a new application architecture.

## 7. Implementation phase 3 — examples and generated cubic assets

### P3.1 Separate configs from canonical assets

Final layout:

```text
examples/
├── README.md
├── assets/
│   ├── common/
│   ├── media/
│   ├── shell/
│   └── volume/{box,box-with-sphere,bunny,dragon}/
├── configs/
│   ├── shell/
│   └── volume/{box,box-with-sphere,bunny,dragon,dragon-dyn}/
├── generated/                 # gitignored
│   ├── cubic/
│   └── output/
└── scripts/generate_cubic_veg.py
```

`dragon-dyn` intentionally references the canonical `dragon` OBJ/tet VEG instead
of storing identical files in a second asset directory.

- [x] Move every retained static mesh, constraint input, and preview into `examples/assets/`.
- [x] Keep each retained asset once; distinguish tet/cubic constraint lists where their vertex indices differ.
- [x] Move all retained simulation and animation JSON files into `examples/configs/`.
- [x] Give simulation configs explicit mesh/contact names and an explicit `contact-model`.
- [x] Preserve sampled volume, IPC volume, sampled shell, and IPC shell examples.
- [x] Point generated output below `examples/generated/output/`.
- [x] Remove legacy scene, `examples/cubic`, and `examples/ipc` trees with no compatibility copies.
- [x] Remove superseded generated surfaces, checked-in cubic VEG files, and noncanonical cubic preview media.
- [x] Verify a SHA-256 audit reports no duplicate retained assets.
- [x] Update tests, CMake definitions, Python smoke paths, and README commands.
- [x] Ignore the complete `examples/generated/` tree.

### P3.2 Deterministic cubic generation

`examples/scripts/generate_cubic_veg.py`:

- [x] Accept either a `cubicMesher` executable or a build directory.
- [x] Support selected presets and a custom input surface/resolution/material invocation.
- [x] Keep separate full sampled and lightweight IPC presets for `box-with-sphere`.
- [x] Generate `<scene>-cubic.veg` into a caller-selected fresh directory.
- [x] Default to `examples/generated/cubic/`.
- [x] Refuse to overwrite an existing output directory.
- [x] Generate through a temporary sibling directory and publish only after all selected meshes validate.
- [x] Record input/output hashes, mesher path, source revision, arguments, counts, bounds, cell-size range, and validation results in `manifest.json`.
- [x] Validate finite vertices, cubic element arity/index ranges, and positive axis-aligned rest cells.
- [x] Reproduce the previous box cubic VEG byte-for-byte (matching SHA-256).
- [x] Keep generated cubic VEG and manifests out of git.

### P3.3 Test integration and documentation

- [x] Generate the cubic box fixture in the CMake build tree with `cubicMesher`.
- [x] Make every cubic-consuming native test depend on that generated fixture.
- [x] Materialize test-only configs in the build tree; tests do not write generated meshes into the source tree.
- [x] Keep tet and shell tests on canonical checked-in assets.
- [x] Provide a bounded all-config smoke driver backed by tetMesher/cubicMesher-generated lightweight tet and cubic assets; execute setup, assembly, solve, finite-output validation, and animation conversion for every retained config.
- [x] Add `examples/README.md` with generation, runner, and migration commands.
- [x] Add an old-to-new path table and state that old paths are intentionally unsupported.
- [x] Update the root README to the separated layout.

Acceptance: the source tree contains no generated cubic VEG, no unintended duplicate
example assets, and no operational references to removed paths; old paths appear only
in migration documentation and release checks.

## 8. Implementation phase 4 — build system, dependencies, and wheels

> **Current `add-no-conda` decision (2026-07-29):** release wheels are
> standalone pip-installable artifacts. Conda is not part of configure,
> packaging, CI, or runtime. The older Conda-based Phase 4 text retained below
> is historical design context only and must not be used as implementation
> guidance.

The active dependency contract is:

- CMake keeps the existing public `PGO_*` feature flags, except
  `PGO_CHECK_CONDA`, which is removed and rejected with a migration error.
- TBB is supplied by the platform packaging environment. Linux and Windows use
  the TBB runtime selected by Intel's pinned `mkl-devel` package; macOS uses
  Homebrew TBB because Intel does not publish a macOS arm64 package.
- Linux and Windows use pinned Intel PyPI MKL with the TBB threading layer.
  macOS uses the system Accelerate framework and does not link MKL.
- GMP and MPFR come from the manylinux system packages or Homebrew and are
  bundled by `auditwheel`/`delocate`. Windows keeps the repository's approved
  prebuilt GMP/MPFR DLLs and bundles them with `delvewheel`.
- Imath is built statically from its pinned upstream source archive. Other
  ordinary source dependencies continue to use pinned FetchContent inputs.
- CI repairs each raw wheel, audits its native dependency closure, installs it
  in a fresh ordinary `venv`, and runs import/API smoke without inherited
  library search paths.
- Release artifacts are produced only by the three platform workflows and are
  not tracked in git.

Acceptance: a user can install the repaired wheel plus its declared Python
requirements into a fresh CPython 3.12 virtual environment without Conda or
separately installed native runtime libraries.

### Active P4.1 — one clear CMake development path

- [x] Keep one public `pypgo-wheel` configure preset.
- [x] Make MKL default to ON and force it OFF on macOS.
- [x] Put the Ninja generator and `build/pypgo` binary directory in the
  preset; callers do not pass `-B` or `-G`.
- [x] Enable `BUILD_TESTING`, Python, subprojects, Alembic, and portable build
  behavior in the preset.
- [x] Defer gtest discovery to CTest with `PRE_TEST`, avoiding false
  post-link discovery timeouts during unrestricted parallel builds.
- [x] Remove obsolete no-MKL and `base_no_mkl` presets.
- [x] Preserve the existing public `PGO_*` feature flags other than the
  intentionally removed Conda check.
- [x] Keep `CMAKE_ARGS`, macOS `ARCHFLAGS`, and setuptools parallel-build
  compatibility. The packaging backend may override only the binary directory
  to use setuptools' isolated temporary tree.

Developer acceptance:

```bash
uv sync --locked
uv run cmake --preset pypgo-wheel
uv run cmake --build build/pypgo
uv run ctest --test-dir build/pypgo --output-on-failure
PYTHONPATH="$PWD/build/pypgo/src/python/pypgo" \
  uv run python -m pytest -q tests/pypgo/test_pgo_smoke.py
```

### Active P4.2 — uv environment and single version sources

- [x] Use `.python-version`, `pyproject.toml`, and one cross-platform
  `uv.lock`; no Conda environment is part of the supported path.
- [x] Pin the required `uv` version once in `pyproject.toml`; CI bootstraps
  that version through `scripts/uv_version.py`.
- [x] Keep the project version once in `pyproject.toml`; CMake and CI artifact
  naming read it rather than duplicating `0.0.4`.
- [x] Let CMake discover the active Python prefix and standard platform
  package locations without requiring `MKL_DIR`, `TBB_DIR`, `MKLROOT`, or
  `CMAKE_PREFIX_PATH` in normal commands.
- [x] Do not install pypgo during ordinary development. Import the extension
  from the build tree and rebuild incrementally after C++ edits.
- [x] Keep wheel construction explicit and separate through
  `uv build --wheel --no-build-isolation`.

### Active P4.3 — dependency acquisition and pinning

- [x] Remove `PGO_CHECK_CONDA` and all Conda discovery paths.
- [x] Pin FetchContent archives/commits and hashes where practical.
- [x] Build Imath statically from pinned source.
- [x] Use Intel PyPI `mkl-devel` and its matching TBB packages on
  Linux/Windows; do not independently build TBB there.
- [x] Use Homebrew TBB/GMP/MPFR on macOS and system GMP/MPFR packages in the
  manylinux container.
- [x] Use approved repository prebuilt GMP/MPFR files on Windows.
- [x] Retain third-party notices and dependency audit scripts.
- [ ] Resolve the remaining CGAL/SuiteSparse binary-distribution obligations
  and approve the final license/source-offer package before release.

### Active P4.4 — wheel repair, audit, and provenance

- [x] Keep 0.0.4 wheel-only; do not build an sdist or commit binary artifacts.
- [x] Require a clean checkout and exact commit before packaging.
- [x] Build raw wheels with `uv build`.
- [x] Configure Linux to repair with `auditwheel`, macOS with `delocate`, and
  Windows with `delvewheel`.
- [x] Configure CI to audit the repaired native dependency closure and reject
  source/build paths or undeclared external libraries.
- [x] Configure CI to install the repaired wheel into a new CPython 3.12
  virtual environment with no inherited native library search path.
- [x] Configure provenance to record CMake cache, dependency/toolchain
  evidence, native CTest,
  build-tree pytest, installed smoke, installed pytest, source commit, wheel
  metadata, and SHA-256.

### Active P4.5 — platform CI structure

- [x] Keep exactly three workflows: Linux, macOS, and Windows.
- [x] Give each workflow a `build-test` job and a fresh dependent `package`
  job.
- [x] Build/test directly from `build/pypgo` without installing pypgo.
- [x] Transfer native evidence between jobs through a short-lived artifact.
- [x] Build, repair, audit, clean-install-test, and upload exactly one wheel
  plus evidence in the package job.
- [x] Keep Actions pinned, permissions read-only, concurrency cancellation,
  explicit runner/container targets, and bounded timeouts.
- [x] Do not impose `CMAKE_BUILD_PARALLEL_LEVEL`, MKL runtime tuning
  environment variables, or redundant preset-selection variables.
- [ ] Run all six hosted jobs on the exact candidate commit and attach their
  run/artifact IDs to the release record.

### Active P4.6 — completed local evidence and remaining gate

Completed on macOS arm64:

- [x] `uv lock --check` and `uv sync --locked --check`.
- [x] Configure through only `cmake --preset pypgo-wheel`.
- [x] Unrestricted-parallel full build.
- [x] Native CTest: 154/154 passed.
- [x] Build-tree Python smoke: 2/2 passed.
- [x] Raw macOS 26 arm64 wheel build.
- [x] VitePress documentation build and sidebar validation.
- [x] Workflow YAML parsing and `git diff --check`.

Phase 4 exit condition: the implementation is complete and local macOS
evidence passes. Hosted Linux/macOS/Windows job results remain Phase 5/RC
evidence and are not implied by local success.

### Historical Phase 4 plan (superseded; non-operative)

#### Historical P4.1 Dedicated pypgo configure presets

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

- Scope rule: these presets are aliases for the existing `setup.py` build
  options. P4.1 does not replace dependency discovery, enforce a new BLAS
  provider, remove the existing Windows DLL handling, or otherwise change the
  native build architecture. Dependency and linkage auditing belongs in the
  platform CI phases.
- [x] Factor common Python options into a hidden preset.
- [x] Enable `PGO_ENABLE_PYTHON`.
- [x] Preserve the existing `PGO_BUILD_SUBPROJECTS=ON`,
  `PGO_ENABLE_ALEMBIC=ON`, and `PGO_CHECK_CONDA=ON` behavior.
- [x] Keep the target name `pypgo`.
- [x] Make MKL choice explicit in the two presets.
- [x] Set `PGO_USE_MKL=ON` in `pypgo-wheel-mkl`.
- [x] Set `PGO_USE_MKL=OFF` in `pypgo-wheel-no-mkl`.
- [x] Let `setup.py` read `PYPGO_CMAKE_PRESET`.
- [ ] Make CI set `PYPGO_CMAKE_PRESET` explicitly.
- [x] When `PYPGO_CMAKE_PRESET` is unset, preserve the original platform
  behavior and direct CMake invocation unchanged.
- [x] Continue configuring into the existing setuptools temporary build directory.
- [x] Continue supporting `CMAKE_ARGS` with its existing parsing behavior.
- [x] Preserve macOS `ARCHFLAGS`.
- [x] Preserve MSVC architecture/config handling.
- [x] Preserve `CMAKE_BUILD_PARALLEL_LEVEL`.
- [x] Preserve the existing MKLROOT setup and Windows GMP/MPFR DLL copy behavior.

#### Historical P4.2 Version and package metadata

Files requiring audit:

- `/Users/jinceyang/Desktop/codebase/merge/libpgo/setup.py`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/pyproject.toml`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/c/CMakeLists.txt`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/python/pypgo/CMakeLists.txt`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/src/python/pypgo/pypgo_main.cpp`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/README.md`

Checklist:

- [x] Change Python distribution version from `0.0.3` to `0.0.4`.
- [x] Change the C shared-library `VERSION` from `0.0.2` to `0.0.4`.
- [x] Keep C `SOVERSION 0` unless an actual binary-incompatible change is introduced.
- [x] Install both public headers, `pgo_c.h` and the included `pgo_c_def.h`.
- [x] Generate and install `pgoConfigVersion.cmake` with compatibility appropriate for the retained `SOVERSION`.
- [x] Ensure the compiled Python module receives `PYPGO_VERSION_INFO=0.0.4`.
- [x] Ensure installed `pypgo.__version__` is exactly `0.0.4`.
- [x] Ensure wheel filenames contain `0.0.4`.
- [x] Add/verify package description, README content type, license file, project URL, and supported Python metadata.
- [x] Audit and declare actual Python runtime dependencies; in particular, verify whether NumPy is required at import time or API-call time instead of leaving `install_requires=[]` without evidence.
- [x] Define the supported runtime as a documented Conda/Miniforge environment followed by installation of the downloaded CI wheel: all platforms install GMP, MPFR, Imath, fmt, and TBB from `conda-forge`; Linux/Windows additionally install MKL-selected BLAS/LAPACK; on macOS libpgo links system Accelerate while the environment installs `_newaccelerate` BLAS/LAPACK packages so Python/NumPy use the same backend family.
- [x] Remove release installation instructions that use `apt install libgmp-dev libmpfr-dev` or `brew install gmp mpfr imath`; provide one platform-adjusted Conda environment recipe instead.
- [x] State explicitly that pip cannot install Conda packages, so these wheels are not claimed to work in an otherwise empty `venv`.
- [x] Keep Conda-provided GMP, MPFR, Imath, fmt, MKL, and TBB external to the wheels, reject accidental vendoring during package audit, and document the exact Conda packages/channels required before installing the wheel.
- [x] Define the release wheel guarantee explicitly: CPython 3.12 wheels tagged `linux_x86_64`, macOS arm64, and `win_amd64`.
- [x] Keep `python_requires >= 3.9` only if source-build CI verifies the claimed range; otherwise narrow the metadata to the range actually tested.

#### Historical P4.3 Wheel-only source provenance

- Scope rule: `scripts/release_wheel_provenance.py` defines and enforces the
  wheel-only source/evidence contract without changing the existing setuptools
  build. P4.6 invokes it after each platform's tests and supplies the
  platform-specific dependency and test evidence.
- [x] Do not build, upload, or document an sdist for 0.0.4.
- [x] Build each wheel directly from the exact checked-out commit after that platform's native and Python tests pass.
- [x] Record `git rev-parse HEAD`, workflow run ID, runner image, CMake cache, dependency evidence, and wheel SHA-256 beside each artifact.
- [x] Require a clean tracked worktree before packaging so uncommitted source changes cannot enter a wheel.

#### Historical P4.4 Pin third-party dependencies

Mandatory floating-reference fixes in `/Users/jinceyang/Desktop/codebase/merge/libpgo/CMakeModules/third-party/`:

- [x] `backward.cmake`: replace `master`.
- [x] `libigl.cmake`: replace `main`.
- [x] `cuCollections.cmake`: replace `dev`.
Additional audit:

- [x] Verify every other FetchContent URL is a fixed release URL or exact commit.
- [x] Do not upgrade dependency versions without a release need.
- [x] Record the final dependency table in docs.
- [x] Audit dependency licenses and retain required notices for source and vendored wheel contents.
- [x] Remove the unnecessary `igl_copyleft::cgal` link from the default wheel path.
- [ ] Resolve the remaining CGAL/SuiteSparse binary-distribution obligations and approve the final wheel license/source-offer package before release.
- [x] Add archive hashes where practical.
- [x] Ensure optional dependencies are fetched only when their feature is enabled.
- [ ] Pin the Miniforge installer version, GitHub Actions revisions, Python version, build tools, and critical Conda dependency versions sufficiently to avoid `latest` changing the release build unexpectedly.
- [ ] Use `conda-forge` with strict channel priority and remove/disable default channels in release CI.
- [ ] Retain `conda list --explicit`, `conda list`, channel configuration, and relevant CMake cache entries as release-build evidence.
- [ ] Explicitly pass every important CMake option in CI.

#### Historical P4.5 Remove tracked binary release artifacts

Current tracked artifacts are under:

- `/Users/jinceyang/Desktop/codebase/merge/libpgo/dist/macosx15_arm/`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/dist/ubuntu22.04/`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/dist/ubuntu24.04/`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/dist/win11/`

Checklist:

- [x] Remove all tracked 0.0.2 and 0.0.3 wheel files.
- [x] Ignore `dist/`.
- [x] Ignore `wheelhouse/`.
- [x] Ignore standard package build output.
- [x] Replace the per-wheel `.gitignore` entries with directory rules.
- [x] Update README installation instructions to download the appropriate wheel from the recorded GitHub Actions artifact and install it in the documented Conda environment.
- [x] Keep release wheels and evidence in GitHub Actions artifacts, not in git; do not claim PyPI or sdist availability.

#### Historical P4.6 Platform CI

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

Phase 5 starts now. Execute it in this order:

1. Commit the Phase 4 implementation so provenance can enforce a clean,
   immutable source revision.
2. Run all three hosted workflows and require both `build-test` and `package`
   to pass on that exact revision.
3. Classify every P5.1 item as already covered by the 154-test suite or add the
   smallest missing regression; do not duplicate tests merely to rename them
   as release tests.
4. Complete the external C/CMake consumer and 0.0.3 compatibility checks.
5. Run the Linux sanitizer/numerical gate.
6. Prepare and execute the owner-controlled Linux MKL Pardiso Tier 1/Tier 2
   matrix only after the candidate source and generated assets are stable.

Local macOS results are a development baseline, not final cross-platform
evidence. Unchecked items below remain unchecked until their result is tied to
the exact candidate commit or explicitly audited as covered by such a result.

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
- [ ] C ABI runner matrix.
- [ ] C API Hessian export preserves finite fractional `double` values.
- [ ] Shell sampled and IPC regression tests keep using their established runners and canonical assets.

### P5.2 Linux server API smoke matrix

Use small dedicated test configs and assets under `/Users/jinceyang/Desktop/codebase/merge/libpgo/tests/`, not large showcase runs.

These five short dynamic API cases run on the Linux server after implementation. Hosted CI only needs parser/unit coverage and a bounded representative simulation smoke.

| Case | C API | Python API | CLI | Expected |
| --- | --- | --- | --- | --- |
| tet + sampled | Required | Required | `runSim` | Setup, one or more solves, finite output, return 0 |
| tet + IPC | Required | Required | `runIPCSim` | Dynamic Hessian + CCD path executed, return 0 |
| cubic + sampled | Required | Required | `runSim` | Cubic FEM path executed, return 0 |
| cubic + IPC | Required | Required | `runIPCSim` | Cubic + dynamic IPC Hessian + CCD path executed, return 0 |
| shell + IPC | Required | Required | `runIPCSim` | Shell IPC path executes dynamic Hessian/CCD as applicable, return 0 |
| shell + sampled | Not routed | Not routed | `runShellSim` | Existing sampled shell path remains functional |
| missing contact model | Required | Required | `runSim` | Same behavior as tet/cubic sampled config |
| bad contact model | Required | Required | N/A | Useful C/Python dispatcher error, nonzero result |

### P5.3 Compatibility tests

- [ ] Selected main-era tet/sample scenarios are migrated to the new `examples/configs/` and `examples/assets/` paths and preserve their simulation semantics.
- [ ] Old example filesystem paths are absent and are not treated as compatibility entry points.
- [ ] Run one canonical 0.0.3 tet + sampled case in the read-only main checkout and the final 0.0.4 implementation, then compare selected displacement/energy/output observables within recorded tolerances.
- [ ] Existing public C symbols still export.
- [ ] Diff the Linux/macOS/Windows exported C symbol lists against the 0.0.3 reference and reject an unapproved removal or incompatible change.
- [ ] Existing Python smoke tests from `/Users/jinceyang/Desktop/codebase/main/libpgo` that apply to 0.0.4 still pass.
- [ ] `const char *` cleanup does not change binary calling convention.
- [ ] Document that `runSim`, `runIPCSim`, and `runShellSim` remain the public CLI entry points.
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
- [ ] Give every case an isolated output path outside the source tree.
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
- [ ] Repeat sampled short cases through `runSim`, IPC short cases through `runIPCSim`, and sampled shell through `runShellSim`.
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

- `git@github.com:annajcy/libpgo-doc.git`
- Branch `v0.0.4`
- Frozen starting commit `ef0834c67a280d16836579942246129087c99e44`
- Submodule checkout `/Users/jinceyang/Desktop/codebase/merge/libpgo/docs`

Checklist:

- [x] Initialize/update the superproject's `docs` submodule, then check out the canonical `v0.0.4` branch from the frozen starting commit inside it; do not create a separate docs worktree or substitute a `release/0.0.4` docs branch.
- [x] Preserve the docs-site design, theme/layout, navigation machinery, and reusable site assets, but simplify and repair the build/deployment workflow where the current experimental assumptions do not match 0.0.4.
- [x] Treat the starting documentation content as experimental and nonauthoritative.
- [x] Derive the 0.0.4 content and release notes from the code difference between `/Users/jinceyang/Desktop/codebase/main/libpgo@ed9675d56403b73a0cf32084e5da9b0011f88e2c` and the final 0.0.4 code.
- [x] Do not use the empty docs `main` branch as a content or release-note reference.
- [x] Delete or rewrite pages, sections, examples, API descriptions, and navigation entries that do not exist in the actual 0.0.4 code.
- [x] Remove references to generalized Neo-Hookean, tricubic Hermite, experimental package architecture, ARPACK restoration, and other unreleased features.
- [x] Replace the workflow's checkout of the libpgo default branch with the exact authorized `v0.0.4` tag/commit.
- [x] Remove the workflow's dependency on experimental-only `libpgo-src/environment.yml` and the nonexistent `pypgo-ci` preset; use the actual 0.0.4 docs/source build path and only the dependencies needed to render the site.
- [ ] If GitHub requires the Pages workflow to exist on the default docs branch, keep only the minimal deployment workflow on docs `main` and make it check out the `v0.0.4` content branch explicitly; do not use empty `main` as a content source.
- [x] Keep Pages permissions only in the docs deployment workflow; the three libpgo platform CI workflows remain read-only.
- [x] Keep the SSH submodule URL
  `git@github.com:annajcy/libpgo-doc.git` in
  `/Users/jinceyang/Desktop/codebase/merge/libpgo/.gitmodules`.
- [ ] Optionally record `branch = v0.0.4`, but pin the code repository gitlink to the exact final docs commit.
- [ ] Commit docs content changes in the docs submodule first; after that exact docs commit exists, update and commit the superproject gitlink.
- [ ] Make CI checkout submodules recursively.
- [ ] Build the exact `v0.0.4` docs content in CI and deploy only after the local/CI link and snippet checks pass.
- [x] Check internal links and code snippets.

### P6.2 Detailed release-note structure

Create detailed pages covering:

- [x] 0.0.4 overview.
- [x] Cubic mesh and cubic FEM.
- [x] Cubic mesher and the Python `.veg` generation script; state that generated cubic `.veg` files are not stored in git.
- [x] IPC core, barrier, self-contact, floor contact, and CCD.
- [x] Hessian topology contract.
- [x] Newton symbolic-factorization fix.
- [x] Material and contact maximum-step behavior.
- [x] Extracted sampled/IPC runners and the preserved `runSim`/`runIPCSim`/`runShellSim` split.
- [x] C API.
- [x] Python API.
- [x] Config routing: tet/cubic × sampled/IPC plus the retained shell runners.
- [x] Build presets.
- [x] Wheel installation from GitHub Actions artifacts; explicitly state that 0.0.4 has no sdist and is not on PyPI.
- [x] Linux/Windows MKL Pardiso wheel policy, TBB threading layer, macOS no-MKL policy, and Linux server validation environment/matrix summary.
- [x] Dependency versions and optional features.
- [x] Migration from 0.0.3.
- [x] Old-to-new example config and asset path map, explicitly stating that old paths are not supported.
- [x] Known limitations.

Known limitations must say:

- [x] IPC external objects are not supported in 0.0.4 unless implementation and tests prove otherwise.
- [x] IPC supports the implemented self-contact and explicit floor path.
- [x] CUDA utilities are retained but not part of default wheels.
- [x] Native Pardiso is retained but not bundled.
- [x] Guaranteed wheel platforms/Python version are listed explicitly.

### P6.3 Main repository release notes

Files:

- `/Users/jinceyang/Desktop/codebase/merge/libpgo/release-notes.txt`
- `/Users/jinceyang/Desktop/codebase/merge/libpgo/README.md`

Checklist:

- [x] Add a concise `0.0.4` section.
- [x] Link to the detailed docs release notes.
- [x] List cubic FEM/mesh/mesher.
- [x] List IPC and CCD maximum-step.
- [x] List sampled/IPC config routing and the volume missing-`contact-model` default.
- [x] List the shared C/Python entry.
- [x] List the preserved `runSim`/`runIPCSim`/`runShellSim` split.
- [x] List correctness fixes.
- [x] List packaging changes.
- [x] List the config/asset-separated example layout and the intentional example-path migration.
- [x] List known limitations.
- [x] Remove instructions that point users to wheel files committed under `dist/`.

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
- [ ] The owner-confirmed Linux server Tier 2 numerical matrix passes with MKL Pardiso.
- [ ] The Linux server report identifies the exact `upstream-0.0.4` integration commit and clean build configuration.
- [ ] All three wheel-only workflows pass and upload one tested wheel each.
- [ ] Repaired wheels pass validation in separate fresh ordinary virtual
  environments without inherited build-time native library paths.
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
- [ ] Verify every recorded SHA-256 and repeat fresh-venv version/import/API
  smoke tests against the downloaded wheel files.
- [ ] Record the exact final workflow URLs, run IDs, artifact IDs, source commit, and checksums in the release record.

### Final tag and artifact delivery

- [ ] Tag the exact tested commit as `v0.0.4`.
- [ ] Verify the tag includes the exact docs submodule commit.
- [ ] Ensure all three platform workflows have completed on the tagged commit and retained their wheel/evidence artifacts.
- [ ] Document how to download each platform wheel from its recorded GitHub
  Actions artifact and install the local file in a CPython 3.12 virtual
  environment.
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
9. `refactor(runner): extract sampled and IPC implementations for API reuse`
10. `feat(api): align C and Python config-runner entry`
11. `refactor(examples): separate configs, deduplicate assets, and add cubic veg generator`
12. `test(release): cover sampled and IPC runner routing`
13. `test(server): add reproducible Linux MKL Pardiso matrix driver`
14. `build(deps): pin release dependency revisions`
15. `chore(release): finalize version and wheel metadata`
16. `build(pypgo): add unified uv and pypgo-wheel build path`
17. `ci(release): split native validation from repaired wheel packaging`
18. `docs(release): rewrite v0.0.4 docs content and repair deployment`
19. `docs(release): finalize release notes`

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
- `runSim` and `runIPCSim` are thin wrappers over reusable implementations; C and Python dispatch to those same implementations.
- Existing executable boundaries are preserved; static IPC reuses
  `runIPCSim` and the shared dispatcher without a new backend or resume architecture.
- Existing shell assets plus sampled and IPC shell paths are retained and regression-tested.
- Example configs and static assets are separated under `examples/configs/` and `examples/assets/`, with no legacy-path compatibility copies or unintended duplicate assets.
- The obsolete CIPC wrapper and scoped profiling are gone.
- Dynamic Hessian and solver symbolic-analysis lifecycles are correct.
- A checked-in Python script reproducibly generates validated cubic `.veg` files into an untracked output directory; generated cubic `.veg` files are not stored in git.
- Optional backward/CUDA/native-Pardiso code remains, with release dependencies pinned and optional features off by default where required.
- Linux and Windows wheels use MKL Pardiso with the matching TBB threading
  layer; the macOS wheel does not use MKL and links libpgo/SuiteSparse to the
  system Accelerate framework. Required non-system runtime libraries are
  bundled and audited by the platform wheel-repair tool.
- Exactly one wheel per platform is built, audited, isolated-tested, checksummed, and uploaded with evidence as a GitHub Actions artifact; no sdist or PyPI release is produced.
- The exact `upstream-0.0.4` integration commit passes the owner-confirmed Linux server matrix using MKL Pardiso with TBB.
- The docs submodule is pinned to the exact final commit of `annajcy/libpgo-doc` branch `v0.0.4`.
- The exact workflow run IDs, artifact IDs, wheels, and checksums validated for the tagged commit are recorded as the 0.0.4 artifact set.
