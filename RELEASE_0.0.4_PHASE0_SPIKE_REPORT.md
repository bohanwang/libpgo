# libpgo 0.0.4 Phase 0 Spike and Handover Report

Date: 2026-07-28

Repository: `/Users/jinceyang/Desktop/codebase/merge/libpgo`

Working branch: `release/0.0.4`

This report records the Phase 0 freeze, baseline, and bounded-spike work performed
before Phase 1. It is intended to be the starting context for a new implementation
session. The detailed target contract remains
[`RELEASE_0.0.4_IMPLEMENTATION_PLAN.md`](RELEASE_0.0.4_IMPLEMENTATION_PLAN.md);
the short release gates remain
[`RELEASE_0.0.4_CHECKLIST.md`](RELEASE_0.0.4_CHECKLIST.md).

## 1. Executive handover

| Area | Phase 0 result | Meaning for the next session |
| --- | --- | --- |
| Branch freeze | Complete locally | Implement on `release/0.0.4`; integrate into `upstream-0.0.4` only after review |
| Native baseline | Pass on macOS arm64 | 143/143 CTest tests pass in the current baseline build |
| C API baseline | Pass | Both current `PgoCTest` tests pass |
| Python source-build baseline | Pass | 2/2 `tests/pypgo` tests pass |
| Static IPC spike | Stopped at a precise core contract gap | Phase 1 must fix dynamic-Hessian initialization before claiming a static IPC solve |
| macOS wheel spike | Wheel builds and imports, but audit fails | The spike is evidence of dependency contamination, not a release candidate |
| Linux/Windows wheel spikes | Not run locally | They remain required in Phase 4 CI; do not infer success from the macOS result |
| Documentation spike | Pass locally | Doxygen, generated pdoc reference, and VitePress build complete successfully |
| Asset/package inventory | Captured below | Phase 3 can use this as the pre-migration baseline |

Phase 0 did not implement the Phase 1 library refactors. Its main outcome is a
known-good macOS native baseline plus explicit contracts for the static IPC,
packaging, and documentation work.

## 2. Frozen identities and branch roles

### Main repository

All of the following local branches currently point to the frozen code commit:

```text
1670f6258d941d80cd9d718fc033fccc193d5e20
Merge pull request #7 from annajcy/pre-merge/temp-0.0.4
```

| Branch | Role |
| --- | --- |
| `release/0.0.4` | Active implementation branch |
| `upstream-0.0.4` | Future reviewed integration branch |
| `upstream-temp-0.0.4` | Historical/source branch for the frozen commit; it is not needed for ongoing implementation |

The 0.0.3 compatibility reference remains:

```text
/Users/jinceyang/Desktop/codebase/main/libpgo
ed9675d56403b73a0cf32084e5da9b0011f88e2c
```

The experimental reference remains:

```text
/Users/jinceyang/Desktop/codebase/libpgo
3bc08226c10e3a7dde862693904b49718a766dfa
```

Do not merge the experimental checkout wholesale. Use it only as a source for
explicitly selected changes described by the implementation plan.

### Documentation repository

The main repository now contains `docs` as a submodule:

```ini
[submodule "docs"]
	path = docs
	url = git@github.com:annajcy/libpgo-doc.git
	branch = v0.0.4
```

The docs working branch and frozen starting commit are:

```text
branch: v0.0.4
commit: ef0834c67a280d16836579942246129087c99e44
```

The superproject currently stages the initial gitlink at that frozen commit.
The uncommitted docs edits described later are inside the submodule and therefore
are not yet represented by a newer superproject gitlink.

## 3. Local environments and baseline build

### Environment roles

| Conda environment | Role |
| --- | --- |
| `libpgo-dev` | Persistent development, native build, tests, pypgo build, and docs generation |
| `libpgo-wheel-runtime` | Separate wheel-install smoke environment; currently contains the Phase 0 test wheel |

The later release wheel verification must create fresh runtime environments
again. `libpgo-wheel-runtime` is evidence from this spike, not the final clean
release environment.

### Host and tools

```text
macOS 26.5.1 (25F80), arm64
Apple clang 21.0.0
CMake 4.3.1
Ninja 1.13.2
Python 3.12.13
```

### Native baseline configuration

Build directory:

```text
build/p0-baseline-conda-tbb
```

Important cache values:

```text
CMAKE_GENERATOR=Ninja
CMAKE_BUILD_TYPE=Release
CMAKE_PREFIX_PATH=/Users/jinceyang/miniconda3/envs/libpgo-dev
CMAKE_DISABLE_FIND_PACKAGE_OpenMP=TRUE
PGO_CHECK_CONDA=ON
PGO_ENABLE_FULL=ON
PGO_ENABLE_PYTHON=ON
PGO_USE_MKL=OFF
PGO_HAS_ORIG_PARDISO=OFF
PYTHON_EXECUTABLE=/Users/jinceyang/miniconda3/envs/libpgo-dev/bin/python
```

The build uses system Accelerate for BLAS/LAPACK. GMP, GMPXX, MPFR, Imath, and
TBB resolve from `libpgo-dev`, except that the current cache still resolves
`fmt_DIR` from `/opt/homebrew/lib/cmake/fmt`. This is acceptable only as a
recorded baseline defect; Phase 4 release configuration must reject it.

The `build/p0-baseline-conda-tbb` cache also contains the stale, unrecognized entry
`PGO_BUILD_P0_STATIC_IPC_SPIKE=ON` from the temporary investigation. There is no
corresponding source option or retained spike target. Phase 1 may continue in the
existing preset-owned `build/base_no_mkl` tree; do not use the Phase 0 evidence
tree as the Phase 1 build.

### Network setup used for FetchContent

SuiteSparse FetchContent initially appeared to hang because the direct GitHub
download was not progressing. The local proxy made the fetch complete:

```bash
export http_proxy=http://127.0.0.1:17891
export https_proxy=http://127.0.0.1:17891
export all_proxy=http://127.0.0.1:17891
```

These values are machine-local conveniences, not release configuration and not
values to commit.

### macOS Conda TBB discovery

`CMakeModules/third-party/tbb.cmake` was changed so macOS participates in the
existing Conda TBB lookup:

```cmake
if(CMAKE_SYSTEM_NAME STREQUAL "Linux" OR APPLE)
  set(CANDIDATE_PATH "$ENV{CONDA_PREFIX}/lib/cmake/TBB")
```

This prevents macOS from unnecessarily falling through to the oneTBB GitHub
FetchContent path when Conda TBB is present.

## 4. Baseline test evidence

The following results were rechecked on 2026-07-28.

### Native CTest

```bash
conda run --no-capture-output -n libpgo-dev \
  ctest --test-dir build/p0-baseline-conda-tbb \
  --output-on-failure -j4
```

Result:

```text
100% tests passed, 0 tests failed out of 143
```

The original run had 144 tests and one baseline-only failure:
`compile_flags_regression`. That test assumed Makefile-generator output and
required `flags.make`, which Ninja does not generate. Per owner direction, the
test registration and `tests/cmake/compile_flags_regression.cmake` were removed.
The resulting supported baseline is 143/143.

### C API

```bash
conda run --no-capture-output -n libpgo-dev \
  ctest --test-dir build/p0-baseline-conda-tbb \
  -R '^PgoCTest\.' --output-on-failure
```

Result:

```text
2/2 passed
PgoCTest.LoadsTetMeshFromExampleFile
PgoCTest.CreatesSmoothRSEnergyWhenSupported
```

This proves the existing C API smoke surface is usable. It does not satisfy the
larger Phase 2 ABI, error-boundary, repeated-call, and five-case matrix.

### Python API from the native build

```bash
PYTHONPATH=/Users/jinceyang/Desktop/codebase/merge/libpgo/build/p0-baseline-conda-tbb/src/python/pypgo \
  conda run --no-capture-output -n libpgo-dev \
  python -m pytest -q tests/pypgo
```

Result:

```text
2 passed
```

The built extension imports from:

```text
build/p0-baseline-conda-tbb/src/python/pypgo/pypgo.cpython-312-darwin.so
```

It exposes 45 public names. This proves the existing pypgo mesh/matrix smoke
surface is usable, but it does not prove the future runner matrix.

## 5. P0.3 spike A — static IPC integration

### Goal

Determine whether existing IPC energy, derivatives, CCD maximum-step, and static
Newton infrastructure can be connected without inventing another contact
algorithm.

### Result

The bounded spike did not reach a finite static solve. It stopped at a precise
core interface mismatch, as required by the Phase 0 stop rule:

1. `MappedSurfacePotentialEnergy::isHessianTopologyFixed()` returns false.
2. `PotentialEnergies::init()` nevertheless calls `createHessian()` on every
   child energy before checking whether its topology is fixed.
3. `MappedSurfacePotentialEnergy::createHessian()` intentionally throws and
   instructs callers to use `hessianDirect()`.

Therefore, placing mapped IPC contact into the aggregate used by static Newton
fails during aggregate initialization, before Newton can solve.

Relevant source locations:

```text
src/core/nonlinearOptimization/potentialEnergies.cpp:50
src/core/contact/mappedSurfacePotentialEnergy.cpp:90
src/core/contact/mappedSurfacePotentialEnergy.h:42
```

The existing aggregate already has a dynamic branch in
`PotentialEnergies::hessianDirect()` that evaluates a current local Hessian and
scatters it into the global matrix. The existing Newton solver also rebuilds and
reanalyzes its solver for a dynamic-topology energy. The immediate blocker is
the inconsistent initialization contract; Phase 1 must then lock all three
behaviors down with regression tests instead of relying on incidental behavior.

### Phase 1 implementation contract

P1.3 must:

- call `createHessian()` only for fixed-topology children during aggregate
  initialization;
- define a non-throwing, correctly sized empty-template behavior for a dynamic
  energy's `createHessian()`;
- preserve complete fixed-plus-dynamic assembly exactly once in
  `hessianDirect()`;
- add mixed fixed/dynamic and changing-contact-pattern tests.

P1.4 must:

- reuse symbolic analysis for fixed topology;
- retain reanalysis when the dynamic sparse pattern may change;
- remove the current fixed-topology `solver.reset()` after each solve.

After P1.3 and P1.4, rerun the tiny tet + IPC static spike and require a finite
result before expanding to the Phase 2 five-case static matrix.

### Scope conclusion

No new contact algorithm is indicated. The missing capability is the aggregate
dynamic-Hessian contract and its solver lifetime semantics. Temporary spike code
was removed; no exploratory executable or test target remains in the source
tree.

## 6. P0.3 spike B — Python wheel

### Build command

The macOS arm64 no-MKL wheel was built with the development environment:

```bash
mkdir -p build/p0-wheel-spike

http_proxy=http://127.0.0.1:17891 \
https_proxy=http://127.0.0.1:17891 \
all_proxy=http://127.0.0.1:17891 \
CMAKE_ARGS='-DCMAKE_DISABLE_FIND_PACKAGE_OpenMP=TRUE -DPGO_USE_MKL=OFF' \
CMAKE_BUILD_PARALLEL_LEVEL=4 \
conda run --no-capture-output -n libpgo-dev \
  python -m pip wheel --no-build-isolation --no-deps \
  --wheel-dir build/p0-wheel-spike .
```

Artifact:

```text
build/p0-wheel-spike/pypgo-0.0.3-cp312-cp312-macosx_26_0_arm64.whl
```

The wheel contains the extension plus metadata and no copied runtime dylibs.

### Runtime smoke

Installation into the separate `libpgo-wheel-runtime` environment with
`pip install --no-deps` succeeded. The installed extension imports, exposes 45
public names, and `pip check` reports no broken Python requirements.

This import success is not a clean dependency proof: absolute Homebrew install
names allow the host machine to satisfy dependencies that the runtime Conda
environment did not provide.

### Audit findings

The wheel is not a valid 0.0.4 release candidate:

1. Package metadata still says `0.0.3`; `setup.py` has
   `version="0.0.3"`.
2. The C shared library target still has `VERSION 0.0.2`.
3. `otool -L` records Homebrew dependencies:

   ```text
   /opt/homebrew/opt/fmt/lib/libfmt.12.dylib
   /opt/homebrew/opt/gmp/lib/libgmpxx.4.dylib
   /opt/homebrew/opt/mpfr/lib/libmpfr.6.dylib
   /opt/homebrew/opt/gmp/lib/libgmp.10.dylib
   ```

4. The extension also uses approved dependency families through `@rpath`:

   ```text
   @rpath/libtbbmalloc.2.dylib
   @rpath/libtbb.12.dylib
   @rpath/libImath.30.dylib
   ```

5. It correctly links the system Accelerate framework and has no MKL linkage.
6. `LC_RPATH` contains the absolute build-environment path:

   ```text
   /Users/jinceyang/miniconda3/envs/libpgo-dev/lib
   ```

7. The wheel build cache confirms mixed discovery: GMP/MPFR/Imath came from
   Conda in the native baseline, while `fmt_DIR` came from Homebrew; the wheel
   itself resolved fmt and GMP/MPFR to Homebrew install names.

### Contract carried into Phase 4

The implementation plan and checklist now explicitly require:

- GMP/GMPXX, MPFR, Imath, fmt, and TBB discovery from the active Conda prefix;
- failure on `/opt/homebrew`, `/usr/local`, source/build-tree paths, a build
  Conda prefix, or another undeclared host prefix;
- separate build and runtime Conda environments;
- no repair step that silently vendors the external Conda runtime;
- direct binary inspection on every platform;
- fresh runtime installation without inherited library-search paths.

Platform audit tools are:

| Platform | Required evidence |
| --- | --- |
| Linux | `auditwheel show`, `readelf -d`, `patchelf --print-rpath`, and `ldd` |
| macOS | `otool -L`, `otool -l`, and `delocate-listdeps` |
| Windows | `dumpbin /DEPENDENTS`, `delvewheel show`, manifest/metadata, and loader-search inspection |

Linux and Windows spikes were not executed on this macOS host. Their first
complete build/audit/install proof remains Phase 4 work.

## 7. P0.3 spike C — documentation

### Content reduction

The docs branch was reduced from the experimental architecture/parallelism and
large hand-maintained pypgo site to a small 0.0.4-oriented structure:

- landing page;
- source-build and wheel-build guides;
- Doxygen C++ reference;
- generated pdoc Python reference;
- minimal VitePress theme and navigation.

Obsolete parallel-runtime pages, duplicated pypgo prose pages, PDF export
scripts, Playwright/PDF dependencies, and the broken custom PDF icon were
removed.

### Python reference generation

A pybind11 extension has no normal Python source for pdoc to parse. The new
`docs/scripts/generate-pypgo-reference.py` imports the built extension,
enumerates its public classes and functions, emits a temporary pdoc-readable
Python module under `.generated/`, and then runs pdoc.

Discovery order is:

1. explicit `PYPGO_SOURCE_DIR`;
2. an already installed `pypgo`;
3. `../build/*/src/python/pypgo`.

The generator currently produces reference entries for 5 public classes and 40
public functions.

### Local proof

```bash
cd docs

LIBPGO_SOURCE_DIR=.. \
PYPGO_SOURCE_DIR=/Users/jinceyang/Desktop/codebase/merge/libpgo/build/p0-baseline-conda-tbb/src/python/pypgo \
PYPGO_CONDA_ENV=libpgo-dev \
pnpm build:full
```

Result on 2026-07-28:

```text
Doxygen: success
pdoc generation: success
VitePress 1.6.4: success
build complete
```

The C++ and Python API buttons use explicit static `index.html` links and open
in a new tab. This avoids VitePress SPA interception of generated reference
directories.

### Pages workflow

`docs/.github/workflows/pages.yml` is manual-only:

```yaml
on:
  workflow_dispatch:
```

It checks out:

```yaml
repository: annajcy/libpgo
ref: upstream-0.0.4
path: libpgo-src
```

That ref is deliberate: Pages should build from the reviewed integration branch,
not the temporary `merge/v0.0.4` name. Deployment still requires explicit owner
authorization by manually dispatching the workflow.

The workflow has been locally syntax-reviewed and its site build path has been
proven locally, but the GitHub-hosted job has not been dispatched.

## 8. Pre-edit package and example inventory

### Tracked `dist/` wheels

| Bytes | SHA-256 | Path |
| ---: | --- | --- |
| 1,249,566 | `4c1c9c8e0f11bb3982ee8d1bbba584f7bc5290e2ab5bfe01660d5477b318082b` | `dist/macosx15_arm/pypgo-0.0.2-cp311-cp311-macosx_15_0_arm64.whl` |
| 1,519,608 | `e8f60848c951e8ddcc6d1323f2aba9efe406dc72b9cca978260964e4bf483576` | `dist/macosx15_arm/pypgo-0.0.3-cp312-cp312-macosx_26_0_arm64.whl` |
| 5,810,742 | `6a7b25886cf52672901fadd188c5a548bfaf705159f4a1eaddfb99e495ce46c` | `dist/ubuntu22.04/pypgo-0.0.2-cp311-cp311-linux_x86_64.whl` |
| 5,792,761 | `2c867093e3cee87c6896100da70fc6c5f5b123bb4dc6cba1a59c59026954ee6b` | `dist/ubuntu24.04/pypgo-0.0.2-cp311-cp311-linux_x86_64.whl` |
| 6,495,098 | `75c133c5bc8d6c7606570e4eaf8f3bf4c131dd6b06aa0269bfbcda061e192f53` | `dist/ubuntu24.04/pypgo-0.0.3-cp312-cp312-linux_x86_64.whl` |
| 1,344,535 | `62ffc6a6960fdc07eae0dcc2c6dc391f70a872e54ef7aae3d4be9b250f93a85b` | `dist/win11/pypgo-0.0.2-cp311-cp311-win_amd64.whl` |
| 1,504,929 | `c3326b830350fecee4b00fd208bb44157318f1a6fc2cdfdc5456884732c16134` | `dist/win11/pypgo-0.0.3-cp312-cp312-win_amd64.whl` |

These are baseline inventory only. The 0.0.4 contract does not publish an sdist
or commit new release wheels into `dist/`; CI uploads tested wheels and evidence
as Actions artifacts.

### Tracked example summary

```text
120 files
65,348,035 bytes
39 JSON
32 OBJ
20 VEG
17 TXT
8 GIF
2 PNG
2 Markdown
```

### Volumetric mesh baseline

| Vertices | Elements | Bytes | Path |
| ---: | ---: | ---: | --- |
| 243 | 827 | 23,215 | `examples/box-hang/box.veg` |
| 243 | 827 | 23,215 | `examples/box-squash/box.veg` |
| 5,266 | 21,570 | 818,542 | `examples/box-with-sphere/box-with-sphere.veg` |
| 243 | 827 | 23,215 | `examples/box/box.veg` |
| 6,171 | 25,746 | 1,010,774 | `examples/bunny/bunny.veg` |
| 125 | 64 | 4,763 | `examples/cubic/box-hang/box.veg` |
| 125 | 64 | 4,763 | `examples/cubic/box-squash/box.veg` |
| 174 | 80 | 11,453 | `examples/cubic/box-with-sphere-xlite/box-with-sphere-xlite.veg` |
| 57,366 | 51,845 | 6,204,801 | `examples/cubic/box-with-sphere/box-with-sphere.veg` |
| 125 | 64 | 4,763 | `examples/cubic/box/box.veg` |
| 6,084 | 4,695 | 430,972 | `examples/cubic/bunny/bunny.veg` |
| 10,495 | 7,503 | 966,449 | `examples/cubic/dragon-dyn/dragon.veg` |
| 39,979 | 186,736 | 7,922,609 | `examples/dragon-dyn/dragon.veg` |
| 39,979 | 186,736 | 7,922,609 | `examples/dragon/dragon.veg` |
| 125 | 64 | 4,763 | `examples/ipc/cubic/box-hang/box.veg` |
| 125 | 64 | 4,763 | `examples/ipc/cubic/box-squash/box.veg` |
| 174 | 80 | 11,453 | `examples/ipc/cubic/box-with-sphere/box-with-sphere.veg` |
| 243 | 827 | 23,215 | `examples/ipc/tet/box-hang/box.veg` |
| 243 | 827 | 23,215 | `examples/ipc/tet/box-squash/box.veg` |
| 564 | 1,950 | 73,314 | `examples/torus.veg` |

The current tree contains many byte-identical copies, including box tet meshes,
small cubic box meshes, dragon meshes, shell meshes, source OBJ files, and fixed
vertex lists. This confirms the Phase 3 deduplication requirement. Preserve this
inventory until the canonical `examples/assets/` and `examples/configs/`
migration has its own hash/count comparison.

Historical commands and exact generator parameters for every existing cubic
`.veg` file are not encoded in the repository, so they cannot be reconstructed
reliably from Phase 0 evidence. Phase 3 must introduce a checked-in deterministic
generator and record its inputs, resolution, material values, output statistics,
and quality checks. Do not treat the current generated cubic files as having
complete provenance.

## 9. Current uncommitted worktree

At handover, no Phase 0 commit has been created.

Main repository status:

```text
A  .gitmodules
 M CMakeModules/third-party/tbb.cmake
Am docs
 M tests/CMakeLists.txt
 D tests/cmake/compile_flags_regression.cmake
?? RELEASE_0.0.4_CHECKLIST.md
?? RELEASE_0.0.4_IMPLEMENTATION_PLAN.md
?? RELEASE_0.0.4_PHASE0_SPIKE_REPORT.md
```

Meaning:

- `.gitmodules` and the initial docs gitlink are staged.
- TBB discovery and test deletion are unstaged.
- The docs submodule has its own unstaged content/workflow changes.
- The release plan, checklist, and this report are untracked.
- `build/` is ignored and contains local evidence only.

Do not use a destructive reset or checkout when starting Phase 1. Inspect and
preserve these changes. Decide the eventual commit split explicitly; a sensible
split is:

1. branch/submodule and Phase 0 evidence;
2. macOS Conda TBB plus baseline test cleanup;
3. docs submodule commit, followed by a superproject gitlink update.

No commit, push, merge, Pages deployment, release upload, or PR has been
performed.

## 10. Phase 1 start instructions

1. Read this report, the implementation plan, and the checklist before editing.
2. Keep implementation on `release/0.0.4`.
3. Preserve the current dirty worktree and use the existing
   `build/base_no_mkl` tree through the `base_no_mkl` configure preset and
   `base_no_mkl_release` build preset. Its current cache has
   `PGO_CHECK_CONDA=OFF`, so the first Phase 1 configure must explicitly set
   `PGO_CHECK_CONDA=ON`, the `libpgo-dev` prefix, and the macOS no-OpenMP
   baseline options.
4. Implement P1.1 through P1.10 in the plan's order unless a focused dependency
   requires a smaller reorder.
5. For P1.3/P1.4, add regression tests before declaring static IPC unblocked:
   mixed fixed/dynamic aggregation, changing dynamic sparsity, fixed symbolic
   reuse, and dynamic symbolic reanalysis.
6. Immediately rerun the tiny static tet + IPC proof after P1.3/P1.4. A finite
   result is the gate for the Phase 2 five-case static runner/API matrix.
7. Keep the macOS wheel spike classified as failed audit evidence. Packaging
   fixes belong to Phase 4 and must include Conda-only fmt/GMP/MPFR discovery.
8. Configure, build, and test with:

   ```bash
   conda run --no-capture-output -n libpgo-dev \
     cmake --preset base_no_mkl \
     -DPGO_CHECK_CONDA=ON \
     -DCMAKE_PREFIX_PATH=/Users/jinceyang/miniconda3/envs/libpgo-dev \
     -DCMAKE_DISABLE_FIND_PACKAGE_OpenMP=TRUE

   conda run --no-capture-output -n libpgo-dev \
     cmake --build --preset base_no_mkl_release

   conda run --no-capture-output -n libpgo-dev \
     ctest --test-dir build/base_no_mkl --output-on-failure -j4

   PYTHONPATH=/Users/jinceyang/Desktop/codebase/merge/libpgo/build/base_no_mkl/src/python/pypgo \
     conda run --no-capture-output -n libpgo-dev \
     python -m pytest -q tests/pypgo
   ```

9. Do not update `upstream-0.0.4`, commit the docs gitlink, dispatch Pages, or
   publish artifacts until the relevant work is reviewable and explicitly
   authorized.

## 11. Phase 0 decisions that should not be reopened accidentally

- The development environment is named `libpgo-dev`.
- Wheel runtime verification uses a separate fresh Conda environment.
- macOS obtains TBB from Conda and uses system Accelerate, not MKL.
- Linux and Windows use MKL with the TBB threading layer.
- GMP, MPFR, Imath, fmt, and TBB are explicit Conda dependencies on all three
  platforms.
- Release wheels must not be masked by Homebrew, apt, `/usr/local`, build-prefix
  RPATHs, or inherited loader paths.
- The Ninja-specific `compile_flags_regression` test was removed by owner
  direction.
- Docs use the SSH submodule URL and branch `v0.0.4`.
- Docs Pages checks out `annajcy/libpgo@upstream-0.0.4`.
- Pages deployment is manual and requires explicit authorization.
- Generated C++ and Python API references must continue to build successfully.
- The static IPC path reuses the existing IPC algorithm; Phase 1 repairs the
  dynamic-Hessian/solver contract rather than adding another contact method.
