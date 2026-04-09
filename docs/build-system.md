# Build System Documentation

This page describes how libpgo's build toolchain works end to end. It is the
canonical reference for how presets, feature flags, Conan bootstrap, generated
toolchain files, target organization, and Python packaging fit together.

For a Chinese, control-flow-oriented companion series that walks through the
current implementation file by file, see the
[Build System Walkthrough](walkthrough/build_system/index.md).

## Toolchain Layers

| Layer | Files / entry points | Responsibility |
|---|---|---|
| User entry points | `cmake --preset ...`, `uv sync`, `uv build` | Select the native C++ or Python-oriented workflow |
| Preset layer | `CMakePresets.json` | Defines build directories, `Ninja`, build type, and feature-oriented cache variables |
| Configure bootstrap | `cmake/BootstrapConan.cmake` | Normalizes features, exports local recipes, runs `conan install`, and loads the generated Conan toolchain |
| Dependency graph | `conanfile.py`, `conan/recipes/`, `conan/profiles/` | Defines package requirements, repository-local recipes, and host/build compiler profiles |
| Dependency exposure | `cmake/LibpgoConanDeps.cmake` | Resolves Conan-generated packages with `find_package()` and creates compatibility aliases for legacy target names |
| Project build graph | `CMakeLists.txt`, `src/core`, `src/api`, `src/tools` | Builds the core libraries, public APIs, and command-line tools |
| Python packaging | `pyproject.toml`, `scikit-build-core`, `uv` | Drives editable installs and wheel builds on top of the same CMake + Conan backend |

Two design choices shape the whole system:

- Native builds and Python builds share the same C++ dependency bootstrap.
- Feature flags are the main control surface; presets and Python packaging only provide different entry points into the same feature model.

## Configure-Time Control Flow

For a native build such as `cmake --preset core-release`, the configure path is:

1. `CMakePresets.json` selects `Ninja`, the build directory, `CMAKE_BUILD_TYPE`, and `CMAKE_TOOLCHAIN_FILE=cmake/BootstrapConan.cmake`.
2. `cmake/BootstrapConan.cmake` runs before the main project configuration and computes the effective feature set.
3. The bootstrap applies feature implications:
   - `PGO_PROFILE_FULL=ON` enables Python, animation I/O, geometry stack, libigl, geogram, and gmsh.
   - `PGO_FEATURE_PYTHON=ON`, `PGO_FEATURE_LIBIGL=ON`, or `PGO_FEATURE_GEOGRAM=ON` implies `PGO_FEATURE_GEOMETRY_STACK=ON`.
4. The bootstrap exports repository-local Conan recipes by running `conan/recipes/export_recipes.py`.
5. The bootstrap runs `conan install` against `conanfile.py`, passing the normalized feature set as Conan options such as `with_python`, `with_animation_io`, and `with_geometry_stack`.
6. Conan generates `conan_toolchain.cmake` and CMake package metadata under the preset build tree, inside `build/<preset>/.conan-build/<BuildType>/.../generators/`.
7. The root `CMakeLists.txt` loads the generated toolchain and then includes `cmake/LibpgoConanDeps.cmake` to resolve dependency targets.
8. The project build graph fans out into `src/core`, `src/api`, and `src/tools`.
9. `cmake --build --preset ...` then delegates actual compilation to `Ninja`; Conan is not normally rerun at build time.

The main consequence is that dependency resolution is a configure-time concern. When feature inputs and profiles stay unchanged, repeated builds reuse the same generated Conan metadata and the same CMake/Ninja build tree.

## Feature Model and Preset Mapping

libpgo exposes a feature-oriented CMake interface instead of separate hand-written build entry points. The primary options live in the root `CMakeLists.txt` and are mirrored in `conanfile.py`.

### Core feature relationships

- `PGO_PROFILE_FULL` is a convenience switch that expands to the full feature set.
- `PGO_FEATURE_PYTHON`, `PGO_FEATURE_LIBIGL`, and `PGO_FEATURE_GEOGRAM` all imply `PGO_FEATURE_GEOMETRY_STACK`.
- `PGO_FEATURE_MKL` and `PGO_FEATURE_ARPACK` stay independent from the full profile.

### Preset mapping

| Preset | Build type | Effective feature set |
|---|---|---|
| `core-release` / `core-debug` | Release / Debug | Core-only native build |
| `geometry-release` / `geometry-debug` | Release / Debug | Core + geometry stack |
| `animation-release` / `animation-debug` | Release / Debug | Core + animation I/O |
| `python-minimal-release` / `python-minimal-debug` | Release / Debug | Python bindings, geometry stack, tests off |
| `python-release` / `python-debug` | Release / Debug | Python bindings + animation I/O + geometry stack, tests off |
| `full-release` / `full-debug` | Release / Debug | Full profile: Python + animation I/O + geometry stack + libigl + geogram + gmsh |

The presets are thin wrappers. They mostly encode common feature combinations and consistent build directories so contributors do not need to remember long lists of `-D...=ON` flags.

## Conan Integration Details

Conan is the single C++ dependency manager for both the native and Python paths.

### What lives where

- `conanfile.py` declares the dependency graph and generates `CMakeToolchain` and `CMakeDeps` output.
- `conan/recipes/` contains repository-local recipes such as `tetgen`, `ccd-safe`, `ccd-exact`, `asa`, `geogram`, `gmsh`, `alembic`, and `suitesparse`.
- `conan/profiles/` contains repository-selected host/build profiles for Linux, macOS, and Windows CI.

### Why the bootstrap exports local recipes first

Some libpgo dependencies do not come directly from ConanCenter in the exact form this repository expects. The bootstrap script therefore exports the checked-in recipes into the local Conan cache before calling `conan install`. That keeps the whole dependency graph reproducible from the repository itself.

### When Conan reruns

`cmake/BootstrapConan.cmake` writes a feature signature file under the preset build tree, next to the Conan output directory. The signature records:

- build type
- C++ standard
- selected host/build profiles
- normalized feature options

If that signature matches the current configuration and the generated toolchain file already exists, configure skips `conan install`. If any of those inputs change, configure reruns Conan automatically.

### How generated packages become CMake targets

After Conan generates package config files, `cmake/LibpgoConanDeps.cmake` does two things:

1. Calls `find_package()` for the packages enabled by the effective feature set.
2. Creates compatibility aliases such as `TBB::tbb` <-> `onetbb::onetbb`, `boost::boost` <-> `Boost::boost`, and `Ceres::Ceres` <-> `Ceres::ceres`.

That alias layer lets the current source tree keep older target names while still consuming Conan-generated package metadata.

## Target Organization and Outputs

The top-level source graph is split into three CMake subtrees:

- `src/core` builds the internal static libraries such as `mesh`, `simulation`, `contact`, and `solidDeformationModel`.
- `src/api` builds `simulationApi`, the C API shared library `pgo_c`, and, when Python is enabled, the static helper target `pgo_c_static` plus the Python module.
- `src/tools` builds command-line executables such as `runSim`, `cubicMesher`, `tetMesher`, `remeshSurface`, and `animation`.

### Native preset output layout

For preset-based native builds, the root `CMakeLists.txt` sets `PGO_BINARY_ROOT` to the preset build directory:

- static libraries built through `add_libpgo_lib(...)` go to `build/<preset>/lib` on GNU/Clang and `build/<preset>/lib/<Config>` on MSVC
- executables built through `add_libpgo_tools(...)` go to `build/<preset>/bin` on GNU/Clang and `build/<preset>/bin/<Config>` on MSVC
- the shared C API library `pgo_c` also lands in the `bin` output directory for native builds

This layout makes the preset build directory self-contained: compiled libraries, runnable tools, and the Conan metadata all live under the same `build/<preset>/` root.

## Python Build Path

The Python path uses the same C++ backend, but enters it through packaging metadata instead of `CMakePresets.json`.

### Editable development build

`uv sync` follows this chain:

1. `uv` resolves the Python environment from `pyproject.toml`.
2. `scikit-build-core` becomes the build backend for the `pypgo` package.
3. `tool.scikit-build.cmake.args` passes `-DCMAKE_TOOLCHAIN_FILE=../../cmake/BootstrapConan.cmake` plus libpgo feature flags into CMake.
4. The same Conan bootstrap runs during configure, so the Python build reuses the same recipe export and dependency resolution logic as native builds.
5. CMake builds the `pypgo` extension module from `src/api/python/pypgo`.
6. `scikit-build-core` reuses the build tree under `build/scikit-build`, which keeps incremental `uv sync` runs fast.

The current `pyproject.toml` default Python path enables:

- `PGO_ENABLE_TESTS=OFF`
- `PGO_FEATURE_PYTHON=ON`
- `PGO_FEATURE_ANIMATION_IO=ON`
- `PGO_FEATURE_GEOMETRY_STACK=ON`

### CMake-only Python presets

If you want the Python bindings without going through Python packaging, use the native presets:

- `python-minimal-release` / `python-minimal-debug` for bindings without animation I/O
- `python-release` / `python-debug` for bindings with animation I/O and geometry stack

These presets still use the same bootstrap and dependency model; they simply stop at the native CMake layer instead of producing a wheel or editable install.

## CI Parity

The CI workflow in `.github/workflows/ci-cd.yml` intentionally mirrors the local build model instead of introducing a second build system.

- The C++ matrix runs `cmake --preset <preset>` and `cmake --build --preset <preset>` on Linux, macOS, and Windows.
- CI passes repository-selected Conan host/build profiles from `conan/profiles/` into the same configure entry point used locally.
- The Python job uses `uv sync`, with `CMAKE_ARGS` pointing back to the same repository Conan profiles.
- The release job uses `uv build`, so editable builds, tests, and distributable wheels all share the same backend assumptions.

This parity matters because it keeps local debugging and CI debugging aligned: if a preset or profile works locally, the same configuration shape is what CI executes.

## Practical Mental Model

Use these entry points depending on what you are trying to do:

- `cmake --preset ...` when you are working on native C++ code, command-line tools, or preset-based test runs
- `uv sync` when you are iterating on the Python bindings in an editable development environment
- `uv build` when you need distributable Python artifacts

If you remember only one thing, remember this: libpgo has one build backend and multiple front doors. Presets and Python packaging are different ways to drive the same CMake + Conan feature model.
