# libpgo 0.0.4 Release Checklist

This is the short operational checklist for the 0.0.4 release. Detailed implementation instructions live in
[`RELEASE_0.0.4_IMPLEMENTATION_PLAN.md`](RELEASE_0.0.4_IMPLEMENTATION_PLAN.md).

## 1. Frozen release identity

| Field | Value |
| --- | --- |
| Status | Scope frozen; implementation not started |
| Target version | `0.0.4` |
| Implementation repository | `/Users/jinceyang/Desktop/codebase/merge/libpgo` |
| Baseline branch | `upstream-temp-0.0.4` |
| Baseline commit | `1670f6258d941d80cd9d718fc033fccc193d5e20` |
| Working branch | `release/0.0.4`, created from the baseline commit |
| Integration branch | `upstream-0.0.4` |
| 0.0.3 compatibility reference | `/Users/jinceyang/Desktop/codebase/main/libpgo@ed9675d56403b73a0cf32084e5da9b0011f88e2c` |
| Experimental code reference | `/Users/jinceyang/Desktop/codebase/libpgo@3bc08226c10e3a7dde862693904b49718a766dfa` |
| Docs content/infrastructure base | `https://github.com/annajcy/libpgo-doc/tree/v0.0.4` at `ef0834c67a280d16836579942246129087c99e44` |

Release decisions:

- [ ] Do not merge the experimental code branch wholesale.
- [ ] Do not restore or enable ARPACK for 0.0.4.
- [ ] Retain the established sampled and IPC shell paths and migrate their assets/configs into the new separated layout.
- [ ] Do not preserve legacy example paths or compatibility config copies; 0.0.4 switches directly to the new `examples/assets/` and `examples/configs/` paths.
- [ ] Generate cubic `.veg` files with a checked-in Python script into a gitignored output directory; do not commit generated cubic `.veg` files.
- [ ] Linux `linux_x86_64` and Windows `win_amd64` wheels use MKL Pardiso with the TBB threading layer.
- [ ] macOS arm64 wheels do not use MKL.
- [ ] Release wheels are supported inside documented Conda/Miniforge environments: Conda supplies GMP, MPFR, Imath, fmt, TBB, and the other non-system runtime dependencies on all three platforms; Linux and Windows additionally obtain MKL-selected BLAS/LAPACK from `conda-forge`; on macOS, libpgo/SuiteSparse links the system Accelerate framework while Conda supplies `_newaccelerate`-selected BLAS/LAPACK for the Python environment.
- [ ] MKL and TBB are external Conda runtime dependencies and are not duplicated inside the wheels.
- [ ] Do not publish an sdist or upload this release to PyPI; each platform CI workflow uploads its tested wheel and evidence as GitHub Actions artifacts.
- [ ] Hosted CI consists of exactly three independent workflow files: `linux-ci.yml`, `macos-ci.yml`, and `windows-ci.yml`; there is no orchestration, reusable-workflow, or hosted Linux-server-matrix workflow.
- [ ] Preserve `runSim`, `runIPCSim`, and `runShellSim`; extract sampled/IPC implementations into `simulationRunner` for C/Python reuse.
- [ ] Work is merged into `upstream-0.0.4` and all final gates run against its exact merge commit.
- [ ] Opening or merging an `upstream-0.0.4` to `main` PR is an owner decision after review.

## 2. Scope gate

- [ ] Cubic volumetric mesh, cubic FEM, and cubic mesher are retained.
- [ ] Sampled contact is retained.
- [ ] IPC core, barrier energy, CCD maximum-step, self-contact, floor contact, and the existing shell IPC path are retained.
- [ ] The existing tools plus the C/Python dispatcher retain:
  - [ ] tet + sampled
  - [ ] tet + IPC
  - [ ] cubic + sampled
  - [ ] cubic + IPC
  - [ ] shell + IPC
  - [ ] shell + sampled through `runShellSim`
- [ ] Static IPC passes focused tet/cubic/shell runner tests through the
  existing `runIPCSim` and shared dispatcher paths.
- [ ] No separate static/dynamic runner architecture, resume workflow, or
  output transaction layer is included.
- [ ] Missing `contact-model` selects sampled contact for 0.0.3 compatibility.
- [ ] IPC external-object contact is not claimed as supported.
- [ ] No generalized Neo-Hookean, tricubic Hermite, benchmark, experiment, package-architecture, abi3, or ARPACK work is included.
- [ ] Nonessential cleanup is separated from release-blocking work and does not delay the release.

## 3. Correctness and API gate

- [ ] Dynamic-Hessian assembly tests pass.
- [ ] Fixed-topology symbolic-factorization reuse and dynamic-topology safety tests pass.
- [ ] Tet, cubic, and IPC maximum-step tests pass.
- [ ] Cubic rest-state and IPC public-input validation tests pass.
- [ ] `runSim` and `runIPCSim` are thin wrappers over the same implementations used by the C/Python dispatcher.
- [ ] CLI and C return `0` only on success and nonzero on every config/setup/simulation failure.
- [ ] No C++ exception crosses the C ABI boundary.
- [ ] Python failure behavior is documented and tested consistently.
- [ ] A missing or unreadable config is a failure, not success.
- [ ] The C API Hessian export preserves fractional `double` values.
- [ ] Repeated success, failure-then-success, and sequential multi-run calls work in one process.
- [ ] The 0.0.3 tet + sampled compatibility case matches recorded numerical observables within tolerance.
- [ ] Exported C symbols are compared against 0.0.3; no unapproved symbol is removed or incompatibly changed.

## 4. Runner and asset gate

- [ ] Tool wrappers contain CLI handling only; simulation setup/execution lives in `simulationRunner`.
- [ ] C and Python contain no duplicate simulation body.
- [ ] Missing `contact-model` selects sampled; explicit `ipc` reaches the IPC implementation; unknown values fail.
- [ ] Every checked-in static asset has one canonical copy under `examples/assets/`.
- [ ] A content-hash audit finds no unintended duplicate checked-in assets.
- [ ] Every checked-in example config lives under `examples/configs/` and references canonical shared assets.
- [ ] The shell sampled and IPC configs reference the one canonical `examples/assets/shell/` tree.
- [ ] The retained shell configs live under `examples/configs/shell/`.
- [ ] No compatibility copies remain at old paths such as `examples/box/box.json`, `examples/shell/shell.json`, or `examples/ipc/shell/shell-ipc.json`.
- [ ] Shell regression tests use the new separated paths.
- [ ] A checked-in Python script generates cubic `.veg` files deterministically from canonical surface assets.
- [ ] Generated cubic meshes go to a dedicated gitignored directory or an explicit caller-provided output directory.
- [ ] Hosted CI generates only the cubic meshes needed by bounded automated tests; the owner-run Linux server matrix generates its own recorded matrix inputs.
- [ ] The generation report records command, code revision, inputs, parameters, mesh statistics, and quality metrics.
- [ ] Cubic quality uses deterministic validity checks and recorded generation statistics; no generated cubic `.veg` is committed.

## 5. Platform and packaging gate

### Linux x86_64

- [ ] CI creates a pinned Miniforge build environment with Python 3.12, GMP, MPFR, Imath, fmt, `mkl-devel`, `tbb-devel`, and MKL-selected BLAS/LAPACK packages from `conda-forge`; do not install GMP/MPFR/Imath/fmt with `apt`.
- [ ] Native Release build and CTest pass.
- [ ] Wheel is built from the tested checkout with `PGO_CHECK_CONDA=ON`, `PGO_USE_MKL=ON`, `PGO_MKL_LINK_DYNAMIC=ON`, and `PGO_HAS_ORIG_PARDISO=OFF`.
- [ ] Configure evidence shows `MKL_THREADING=tbb_thread`.
- [ ] The installed wheel loads `mkl_tbb_thread` and TBB, and does not load the Intel OpenMP threading layer.
- [ ] The final wheel intentionally uses the `linux_x86_64` platform tag; `auditwheel show` and direct ELF inspection document external Conda libraries and reject vendored MKL/TBB or build-prefix RPATH/RUNPATH.
- [ ] Installation passes in a fresh Conda runtime environment containing the documented GMP/MPFR/Imath/fmt/MKL/TBB runtime packages.

### Windows amd64

- [ ] CI creates a pinned Miniforge build environment with Python 3.12, GMP, MPFR, Imath, fmt, `mkl-devel`, `tbb-devel`, and MKL-selected BLAS/LAPACK packages from `conda-forge`.
- [ ] Native Release build and CTest pass under MSVC.
- [ ] Wheel is built from the tested checkout with `PGO_CHECK_CONDA=ON`, `PGO_USE_MKL=ON`, `PGO_MKL_LINK_DYNAMIC=ON`, and `PGO_HAS_ORIG_PARDISO=OFF`.
- [ ] Configure evidence shows `MKL_THREADING=tbb_thread`.
- [ ] `dumpbin`/`delvewheel` evidence shows the TBB MKL threading layer and excludes the Intel OpenMP threading layer.
- [ ] Dependency inspection rejects vendored GMP/MPFR/Imath/fmt/MKL/TBB libraries and build-prefix or undeclared host paths; installation passes in a fresh Conda runtime environment containing the documented GMP/MPFR/Imath/fmt/MKL/TBB runtime packages.

### macOS arm64

- [ ] CI creates a pinned arm64 Miniforge build environment with Python 3.12, GMP, MPFR, Imath, fmt, `tbb-devel`, and `_newaccelerate` BLAS/LAPACK packages from `conda-forge`; do not install GMP/MPFR/Imath/fmt with Homebrew.
- [ ] Native Release build and CTest pass.
- [ ] The CMake dependency logic finds Conda TBB on macOS through `$CONDA_PREFIX/lib/cmake/TBB` or an explicit equivalent.
- [ ] Wheel is built from the tested checkout with `PGO_CHECK_CONDA=ON`, `PGO_USE_MKL=OFF`, `BLA_VENDOR=Apple`, and LP64 BLAS/LAPACK integers.
- [ ] CMake configure evidence shows that libpgo and SuiteSparse selected Accelerate rather than OpenBLAS, generic BLAS, or MKL.
- [ ] Linked-library inspection confirms system Accelerate linkage and no MKL, OpenBLAS, or non-system LAPACK dependency.
- [ ] `otool`/`delocate-listdeps` inspection rejects vendored external Conda runtimes, absolute build-prefix paths, `/opt/homebrew`, `/usr/local`, and other undeclared host install names/RPATHs; installation passes in a fresh Conda runtime environment containing the documented GMP/MPFR/Imath/fmt/TBB packages.

### Common artifact checks

- [ ] CI uses pinned Miniforge plus strict `conda-forge` channel priority on all three platforms; the resolved package list is retained as evidence.
- [ ] Every wheel job uses separate build and clean runtime environments so the build environment cannot mask missing runtime requirements.
- [ ] In the fresh supported Conda runtime environment, `python -m pip install --no-deps <wheel>` followed by import and API smoke tests succeeds.
- [ ] Linux/Windows runtime checks execute an MKL-backed operation, require `MKL_THREADING_LAYER=TBB` and `mkl_tbb_thread`, and reject `mkl_intel_thread`, `mkl_gnu_thread`, and Intel OpenMP.
- [ ] macOS checks exercise BLAS and LAPACK operations, confirm libpgo links system Accelerate, and confirm the clean Conda Python environment selects `_newaccelerate`.
- [ ] Wheel-content checks reject accidentally vendored MKL/TBB libraries and any unapproved runtime library.
- [ ] Installation documentation explains how to download the platform wheel from the GitHub Actions artifacts, create the required Conda environment, and install the local wheel path; it does not claim PyPI availability.
- [ ] Ubuntu/macOS installation docs remove the `apt install libgmp-dev libmpfr-dev` and `brew install gmp mpfr imath fmt` release paths; the documented release environment installs GMP, MPFR, Imath, and fmt through Conda on every platform.
- [ ] `pypgo.__version__ == "0.0.4"`.
- [ ] `python -m pip check` passes.
- [ ] `python -m twine check` passes.
- [ ] Wheel contents, linked libraries, platform tags, minimum OS/glibc target, and forbidden paths are audited.
- [ ] Each platform workflow uploads exactly one tested wheel plus source-commit, dependency, linkage, test, and SHA-256 evidence as GitHub Actions artifacts.
- [ ] `upload-artifact` uses an explicit repository-supported retention period, and the release record states that GitHub Actions artifacts expire rather than presenting them as permanent package hosting.

## 6. C/CMake installation gate

- [ ] `pgo_c.h` and `pgo_c_def.h` are both installed.
- [ ] `pgoConfig.cmake`, `pgoConfigVersion.cmake`, and exported targets are installed.
- [ ] `find_package(pgo 0.0.4 CONFIG REQUIRED)` works outside the source/build tree.
- [ ] A minimal external pure-C consumer compiles, links, and runs against the installed shared library.
- [ ] The installed package contains no source-tree or build-tree paths.
- [ ] C shared-library `VERSION` is `0.0.4`; `SOVERSION` remains `0` unless an ABI-breaking change is approved.

## 7. Owner-run Linux server release-candidate gate

This matrix is run and confirmed manually by the owner. It is not a GitHub Actions job or a required CI status check. The repository supplies the reproducible driver, configs, generation script, and report format.

- [ ] The exact `upstream-0.0.4` merge commit is recorded and the worktree is clean.
- [ ] The build uses MKL Pardiso, not original Pardiso or the Eigen fallback.
- [ ] The build and runtime use the TBB MKL threading layer.
- [ ] The short C API matrix passes the selected sampled/IPC cases.
- [ ] The short Python API matrix passes the same selected sampled/IPC cases.
- [ ] Sampled cases pass through `runSim`, IPC cases through `runIPCSim`, and sampled shell through `runShellSim`.
- [ ] The representative longer tet/cubic sampled/IPC plus shell IPC matrix passes.
- [ ] Reports include commit, configs, generated-mesh manifest, CMake cache, toolchain/dependency versions, environment, commands, logs, return codes, convergence, and numerical checks.
- [ ] Any affected result is rerun after code, config, generator, dependency, or test-input changes.

## 8. Documentation gate

- [ ] Use `annajcy/libpgo-doc` branch `v0.0.4`, based on `ef0834c67a280d16836579942246129087c99e44`.
- [ ] Perform docs work inside this repository's `/Users/jinceyang/Desktop/codebase/merge/libpgo/docs` git submodule; do not use a separate docs worktree.
- [ ] Preserve the documentation-site design and content tooling, but simplify and repair the deployment workflow as needed for 0.0.4.
- [ ] Delete or rewrite content that describes features not present in the actual 0.0.4 code.
- [ ] Derive the release notes from the code difference between `/Users/jinceyang/Desktop/codebase/main/libpgo@ed9675d56403b73a0cf32084e5da9b0011f88e2c` and the final 0.0.4 code; do not use the empty docs `main` branch as the content reference.
- [ ] Document the new non-backward-compatible example paths, public matrix, shell support, IPC limitations, cubic mesh generation workflow, MKL/TBB wheel policy, installation, and migration notes.
- [ ] Build docs and check internal links and code snippets.
- [ ] Pin the code repository submodule to the exact final docs commit and keep
  its SSH URL `git@github.com:annajcy/libpgo-doc.git`.

## 9. Integration and owner handoff

- [ ] Merge the completed working branch into `upstream-0.0.4`.
- [ ] Record the exact integration merge commit.
- [ ] Run every blocking gate against that commit.
- [ ] Produce one concise release record linking tests, server report, docs commit, artifacts, and checksums.
- [ ] Stop for owner review; do not create or merge the `upstream-0.0.4` to `main` PR without explicit authorization.
- [ ] If a later main merge changes the commit, bind the tag and final artifact evidence to the exact final commit and rerun affected final checks.

## 10. Tag, CI artifact delivery, and incident handling

Tagging, final workflow runs, artifact delivery, and docs deployment require separate explicit authorization.

- [ ] Tag the exact authorized and tested final commit as `v0.0.4`.
- [ ] Run all three platform workflows on the exact tagged commit and retain their wheel/evidence artifacts.
- [ ] Download each CI wheel by artifact ID, verify its recorded SHA-256 checksum, and install it in the documented clean Conda runtime environment.
- [ ] Verify the tag, workflow run links, artifact IDs/checksums, and docs deployment; do not publish an sdist or upload to PyPI.

If a serious defect is discovered after artifact delivery:

- [ ] Mark the affected workflow artifacts and docs with a clear warning in the release record.
- [ ] Disable or remove public links to known-bad artifacts where GitHub permits it, while retaining their hashes and incident record.
- [ ] Build a corrected `0.0.5` rather than silently replacing the recorded `0.0.4` artifact set.

## 11. Evidence record

| Evidence | Commit/artifact | Result/link | Owner confirmation |
| --- | --- | --- | --- |
| Linux native/CTest |  |  |  |
| macOS native/CTest |  |  |  |
| Windows native/CTest |  |  |  |
| Linux server Tier 1 |  |  |  |
| Linux server Tier 2 |  |  |  |
| Linux wheel |  |  |  |
| macOS wheel |  |  |  |
| Windows wheel |  |  |  |
| C ABI/install test |  |  |  |
| Cubic generator/validation report |  |  |  |
| Docs build/final commit |  |  |  |
