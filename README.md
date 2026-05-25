## libpgo: Library for Physically based Simulation (P), Geometric Shape Modeling (G), and Optimization (O)

The library is designed to primarily focus on physically based simulations, geometric shape modeling, and optimization.
The source code extends [VegaFEM](https://viterbi-web.usc.edu/~jbarbic/vega/) and is designed for academic research purposes.

---

## Prebuilt (Experimental)
The wheel package of following platform have been provided for ease of use. They are in `./dist` folder:
- Ubuntu 24.04: `ubuntu24.04/pypgo-0.0.3-cp312-cp312-linux_x86_64.whl`. Note that you still need to install `gmp` and `mpfr` as suggested in the prerequisites. You may `apt install` them if needed.
- Ubuntu 22.04: `ubuntu22.04/pypgo-0.0.2-cp311-cp311-linux_x86_64.whl`. Note that you still need to install `gmp` and `mpfr` as suggested in the prerequisites. You may `apt install` them if needed. This version depends on a lower version of the libc, so it should be more compatible.
- Windows: `win11/pypgo-0.0.3-cp312-cp312-win_amd64.whl`. The package is built under Windows 11, Visual Studio 2022. In theory, it supports other windows platforms.
- MacOS Arm: `pypgo-0.0.3-cp312-cp312-macosx_26_0_arm64.whl`. The package is built under Tahoe 26.0.1 on Apple M3.

Do `pip install ./dist/your-chosen.whl` to install the package. Note that the packages are experimental.

---

## Prerequisites

1. CMake >= **3.29**\
    We use several functionalities that are only supported by 3.29+. 
    > In most cases, both system's CMake and Conda Environment's CMake have a lower version of CMake unfortunately. In this sitation, please install a new CMake into your system. The latest CMake, either pre-built binaries or source files, can be obtained directly from the [official](https://cmake.org/download/) website. Once installed, hook `cmake` to the newly installed one, either by adding the `your-new-cmake/bin` to the front of the `PATH` or by replacing the existing `cmake` executable with the new one.

2. Compilers
    1. GCC **11, 12, 13** for Ubuntu\
        We use C++20, so only GCC 11, 12, and 13 are supported. You can get new gcc using `apt` or compile a new one from its source code.

    2. Apple Clang (We tested on 15.0.0, Mac OS 14.5)\
        Earlier versions might work if it supports C++20.

    3. Visual Studio 2022 (We tested on 17.9.5, Windows)\
        Earlier Visual Studio 2022 versions might work.

3. GMP, MPFR, and TBB for **Ubuntu** and **Mac OS**\
    This can be installed on Ubuntu by

    ```bash
    sudo apt install libgmp-dev libmpfr-dev libtbb-dev
    ```

    Or it can be installed on Mac OS by

    ```bash
    brew install gmp mpfr imath tbb
    ```

4. (Optional) Ninja\
    It can be installed by

    ```bash
    pip install ninja
    ```

    for better compilation performance

5. (Optional) numpy\
    This is used for running tests.

6. (Optional) Blender and ffmpeg\
    These are used by `scripts/render_abc_preview.py` when rendering Alembic
    `.abc` animations to GIF previews. The script finds Blender on `PATH`, via
    the `BLENDER` environment variable, or at
    `/Applications/Blender.app/Contents/MacOS/Blender` on macOS. It finds ffmpeg
    on `PATH` or via the `FFMPEG` environment variable.

## Compilation

### Dependency ownership

The build uses three dependency layers:

- FetchContent-managed C++ dependencies are always built from source by this repository: Eigen, fmt, spdlog, nlohmann_json, SuiteSparse, Ceres, Boost, CGAL, geogram, libigl, Alembic, and nanobind.
- External native SDKs and toolchain packages come from the system, Homebrew, apt, conda, or vendor installers: compilers, CMake, Ninja, GMP, MPFR, Imath, BLAS/LAPACK, TBB, OpenVDB, MKL, CUDA, Gmsh, Knitro, and Pardiso.
- Conda owns the Python API build environment: Python, pytest, numpy, setuptools/wheel, CMake/Ninja, and native runtime packages such as TBB, MKL, and Imath when needed.

For `pypgo`, use one conda environment for both Python packages and native build/runtime packages. This keeps Python, MKL, Imath, and runtime library lookup in the same prefix.

For Gmsh support on Linux, macOS, and Windows, install the conda-forge package into the active build environment before configuring with `PGO_ENABLE_GMSH=ON`:

```bash
mamba install -c conda-forge gmsh
# or
conda install -c conda-forge gmsh
```

### CMake Presets

Native C++ builds in this repository are preset-driven:

```bash
cmake --list-presets
cmake --preset <configure-preset>
cmake --build --preset <build-preset> [--target <target>...]
```

Configure presets define feature flags and build directories. Build presets map to configure presets. Build parallelism is controlled by `CMAKE_BUILD_PARALLEL_LEVEL` or command-line `--parallel`, not by shared presets.

#### Configure Presets

`base` is the default preset for CI and local development. It enables MKL, Alembic, Gmsh, TetWild, and OpenVDB. CUDA, Knitro, and Pardiso are enabled only by explicit presets or cache overrides.

| Configure preset | Binary directory | Purpose / key options |
| --- | --- | --- |
| `base` | `build/base` | Release baseline with MKL (`PGO_USE_MKL=ON`), Alembic, Gmsh, TetWild, and OpenVDB. |
| `base_debug` | `build/base_debug` | `base` in Debug mode. |
| `debug` | *(hidden fragment)* | Inheritance fragment that sets `CMAKE_BUILD_TYPE=Debug`. |
| `knitro` | *(hidden fragment)* | Inheritance fragment enabling Knitro (`PGO_OPT_USE_KNITRO=ON`). |
| `pardiso` | *(hidden fragment)* | Inheritance fragment enabling original Pardiso (`PGO_HAS_ORIG_PARDISO=ON`). |
| `cuda` | *(hidden fragment)* | Inheritance fragment enabling CUDA (`PGO_ENABLE_CUDA=ON`). |
| `base_knitro` | `build/base_knitro` | `base` + `knitro` (Linux). |
| `base_knitro_cuda` | `build/base_knitro_cuda` | `base` + `knitro` + `cuda` (Linux). |
| `all` | `build/all` | `base` + `knitro` + `pardiso` + `cuda` in Release mode (Linux). |
| `all_debug` | `build/all_debug` | `base` + `knitro` + `pardiso` + `cuda` in Debug mode (Linux). |
| `base_cuda` | `build/base_cuda` | `base` + `cuda` in Release mode (Linux/Windows). |
| `base_cuda_debug` | `build/base_cuda_debug` | `base` + `cuda` in Debug mode (Linux/Windows). |

The `debug`, `knitro`, `pardiso`, and `cuda` presets are hidden inheritance fragments: they are meant to be composed into other presets and are not selectable directly.

Machine-specific SDK paths belong in the untracked `CMakeUserPresets.json`, not in the shared presets. For example, local profiles can inherit `all` or `base_cuda` and provide `KNITRO_LIBRARY_HINT`, `PARDISO_LIBRARY_HINT`, or `cudss_DIR` there:

```json
{
    "version": 3,
    "configurePresets": [
        {
            "name": "local-all",
            "inherits": "all",
            "cacheVariables": {
                "KNITRO_LIBRARY_HINT": "/opt/artelys/knitro-15.0.1-Linux64",
                "PARDISO_LIBRARY_HINT": "/opt/panua-pardiso-20240229-linux"
            }
        },
        {
            "name": "local-base-cuda",
            "inherits": "base_cuda",
            "cacheVariables": {
                "cudss_DIR": "C:/Program Files/NVIDIA cuDSS/v0.7/lib/13/cmake/cudss"
            }
        }
    ],
    "buildPresets": [
        {
            "name": "local-all",
            "configurePreset": "local-all",
            "jobs": 32
        },
        {
            "name": "local-base-cuda",
            "configurePreset": "local-base-cuda",
            "jobs": 32
        }
    ]
}
```

```bash
cmake --preset local-all
cmake --build --preset local-all
```

The shared `all` preset remains path-free; use it directly when the required SDKs are already discoverable from the environment.

#### Build Presets

| Build preset | Configure preset | Typical use |
| --- | --- | --- |
| `base` | `base` | Release build with MKL and Alembic. |
| `base_debug` | `base_debug` | Debug build. |
| `base_knitro` | `base_knitro` | Release build with Knitro. |
| `base_knitro_cuda` | `base_knitro_cuda` | Release build with Knitro and CUDA. |
| `all` | `all` | Release build with all optional solvers/features (Linux). |
| `all_debug` | `all_debug` | Debug build with all optional solvers/features (Linux). |
| `base_cuda` | `base_cuda` | Release build with CUDA (Linux/Windows). |
| `base_cuda_debug` | `base_cuda_debug` | Debug build with CUDA (Linux/Windows). |

Platform auto-disable: on macOS, configuring with `base` (or any MKL/CUDA preset) emits a warning and forces `PGO_USE_MKL=OFF` and `PGO_ENABLE_CUDA=OFF`. This means the same `base` preset works across Linux, Windows, and macOS.

### Install pypgo into a conda environment

The Python package follows the CI shape: create one conda environment, install native/Python build requirements into it, then build the extension in place.

Linux with conda-provided MKL:

```bash
conda create -n libpgo -c conda-forge python=3.12 "cmake>=3.29" ninja mkl-devel tbb-devel numpy pytest setuptools wheel
conda activate libpgo

sudo apt-get install -y build-essential libblas-dev libgmp-dev libimath-dev liblapack-dev libmpfr-dev pkg-config zlib1g-dev

python setup.py build_ext --inplace
python -m pytest -q tests/pypgo
```

macOS without MKL:

```bash
conda create -n libpgo -c conda-forge python=3.12 "cmake>=3.29" ninja tbb-devel numpy pytest setuptools wheel
conda activate libpgo

brew install gmp mpfr imath

python setup.py build_ext --inplace
python -m pytest -q tests/pypgo
```

Windows with conda-provided MKL should run from an x64 MSVC developer shell:

```bash
conda create -n libpgo -c conda-forge python=3.12 "cmake>=3.29" ninja mkl-devel tbb-devel imath numpy pytest setuptools wheel
conda activate libpgo

python setup.py build_ext --inplace
python -m pytest -q tests/pypgo
```

For a minimal local install, after native prerequisites are available:

```bash
cd libpgo
conda activate libpgo
python setup.py build_ext --inplace
python -c "import pypgo; print(pypgo.__doc__)"
```

`pypgo` uses `PGO_PYTHON_USE_MKL=auto` by default. The accepted values are:

- `auto`: enable MKL only when `MKLROOT`, `CONDA_PREFIX`, or `CMAKE_PREFIX_PATH` points to an MKL-capable native prefix.
- `on`: require MKL and fail early if no MKL hint is available.
- `off`: always configure `pypgo` with `PGO_USE_MKL=OFF`.

The Python package build enables nanobind bindings directly and does not build the C API by default (`PGO_BUILD_C_API=OFF`). Alembic is enabled on Linux/macOS for `convert_animation_to_abc` and disabled on Windows; Gmsh, TetWild, and OpenVDB are off by default for `pypgo` and can be overridden through `CMAKE_ARGS`. When enabling Gmsh, install `gmsh` into the same conda environment. When enabling OpenVDB, install `openvdb`, `libboost-devel`, and `tbb-devel` into the same conda environment.

`setup.py` automatically adds the active conda prefix to CMake's package search path. Use `CMAKE_ARGS` only when you need extra local SDK paths or feature overrides.

The editable conda build uses a persistent CMake build directory under `build/pypgo-conda-base-<platform>-<python-tag>-<config>/`, so repeated builds reuse downloaded FetchContent dependencies. Override it with `PGO_PYTHON_BUILD_DIR` when you want a separate Python build tree.

After changing C++ binding or library sources, rebuild the extension
incrementally before running Python code:

```bash
python setup.py build_ext --inplace
python src/python/pypgo/pgo_test_01.py
```

Running Python uses the already built extension and does not rebuild it by itself. Use `build_ext --inplace` for the normal edit-build-run loop.

### Install pypgo with pip

```bash
cd libpgo
pip install .
```

To build a wheel package for the current platform and Python version:

```bash
cd libpgo
python setup.py bdist_wheel
```

The generated wheel is written to `dist/pypgo-*.whl` and can be installed with `pip install dist/pypgo-*.whl`.

If `ninja` has been installed, it will compile source files in parallel. If it is not installed,
set `CMAKE_BUILD_PARALLEL_LEVEL` to `n`, where `n` is the number of threads for compilation, to control the parallel compilation.

### Setup without Python

If you want to use the library with your C++ code or modify the source code, you may build it without python.

### Windows & Ubuntu

The default `base` preset enables MKL, Alembic, Gmsh, TetWild, and OpenVDB. Install [MKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html) and any enabled external packages first, then:

```bash
cd libpgo
cmake --preset base
cmake --build --preset base
```

To build without MKL, override the option: `cmake --preset base -DPGO_USE_MKL=OFF`.

> On Windows, a few extra steps are need before running the preset commands above. First, the library should be configured in "x64 Native Tools Command Prompt for VS 2022". In addition, before running the commands above, run `c:\Program Files (x86)\Intel\oneAPI\setvars.bat` to setup the environments for MKL, where `c:\Program Files (x86)\Intel\oneAPI` is the path to the oneAPI installation. Once setup, run above commands. OpenVDB is available when a compatible external package, such as conda-forge `openvdb`, is installed and `PGO_ENABLE_OPENVDB=ON` is set.

> On Ubuntu, a similar procedure is needed. Before configuring the library with presets, run `bash /opt/intel/oneapi/setvars.sh` to setup the MKL environments for the subsequent CMake configuration.

### Mac OS

macOS uses the same `base` preset. MKL and CUDA are not supported there, so configuring `base` emits a warning and forces `PGO_USE_MKL=OFF` and `PGO_ENABLE_CUDA=OFF` automatically.

Release build:

```bash
cd libpgo
cmake --preset base
cmake --build --preset base
```

Debug build:

```bash
cmake --preset base_debug
cmake --build --preset base_debug
```

The `base` preset keeps Alembic, Gmsh, TetWild, and OpenVDB enabled. For Gmsh, install the conda-forge package first, for example `mamba install -c conda-forge gmsh`, then configure with `-DPGO_ENABLE_GMSH=ON` when using custom presets or overrides. For OpenVDB, install the external packages first, for example `mamba install -c conda-forge openvdb libboost-devel tbb-devel`, then configure with `-DPGO_ENABLE_OPENVDB=ON`.

---

## Usage & Test

The primary runnable examples in this repository are now IPC examples driven by `runIPCSim` under `examples/ipc/`.

Build the IPC tools:

```bash
cmake --preset base
cmake --build --preset base --target runIPCSim convertAnimation
```

Run named IPC batches from the JSON config:

```bash
scripts/run_sim_batch.py --config examples/ipc/ipc_batch.json --job squash_regression --dry-run
scripts/run_sim_batch.py --config examples/ipc/ipc_batch.json --job squash_regression --skip-existing
scripts/run_sim_batch.py --config examples/ipc/ipc_batch.json --case cubic_box_with_sphere_lite --overwrite
scripts/run_sim_batch.py --config examples/ipc/ipc_batch.json --job all_ipc_abc
```

The generic batch runner reads [`examples/ipc/ipc_batch.json`](./examples/ipc/ipc_batch.json), runs the stages declared by each job, and defaults jobs without a `stages` field to `runIPCSim` with `--log` followed by `convertAnimation` with the matching per-case `anim.json`. Use [`examples/ipc/README.md`](./examples/ipc/README.md) for the full case list, job definitions, and output-overwrite policy.

The same runner also supports Alembic preview rendering when a case supplies
`render_config`. The render stage calls
[`scripts/render_abc_preview.py`](./scripts/render_abc_preview.py), which uses
Blender to render `.abc` frames and ffmpeg to encode a GIF:

```bash
scripts/render_abc_preview.py --config my_render_config.json --overwrite

scripts/run_sim_batch.py --config examples/ipc/ipc_batch.json --job sim --overwrite
```

Run representative IPC cases from the repo root:

```bash
build/base/bin/runIPCSim examples/ipc/shell/shell-hang/shell-ipc.json
build/base/bin/runIPCSim examples/ipc/shell/shell-drop/shell-ipc.json
build/base/bin/runIPCSim examples/ipc/tet/box-hang/box-ipc.json
build/base/bin/runIPCSim examples/ipc/cubic/box-with-sphere/box-ipc.json
```

Convert dumped frame sequences to Alembic:

```bash
build/base/bin/convertAnimation examples/ipc/shell/shell-hang/anim.json
build/base/bin/convertAnimation examples/ipc/shell/shell-drop/anim.json
build/base/bin/convertAnimation examples/ipc/tet/box-hang/anim.json
build/base/bin/convertAnimation examples/ipc/cubic/box-with-sphere/anim.json
```

For the full IPC case list and per-case notes, see [`examples/ipc/README.md`](./examples/ipc/README.md).

Legacy penalty-based volume contact is available through the same entrypoint:

```bash
build/base/bin/runIPCSim --legacy path/to/legacy-volume-config.json
```

`runIPCSim` accepts both `"sim-type": "dynamic"` and `"sim-type": "static"`. Static mode performs a one-shot Newton solve from the rest state, writes the same unified `states/deform0000.u` and `surface/ret0000.obj` layout as dynamic mode, and does not support `restart-from-u`. Static output is written only after Newton convergence; unconstrained gravity-only static drops, including legacy penalty-contact drops without attachments, are expected to fail instead of producing a partial state.

`--legacy` accepts the old volume JSON shape with either `tet-mesh` or `cubic-mesh` and uses the penalty contact model instead of IPC contact. Legacy static mode preserves the old volume static semantics: it solves elastic, attachment, and external-force energies without adding the legacy penalty contact energies. Shell legacy configs are no longer supported; use the IPC shell examples above for shell simulations.

Solver status is reported through the shared `SolverResult` / `SolveStatus` API used by `NewtonSolver`, `EnergyOptimizer`, and `TimeIntegratorSolver`. Static runs require `Converged` before writing output. Dynamic implicit Euler currently preserves the legacy timestep policy: `Converged`, `MaxIterations`, and `StepTooSmall` are accepted timestep statuses, while other statuses are failures. External solvers keep backend-specific return codes in `SolverResult::rawStatusCode`.

Optional Python API smoke test:

```bash
python src/python/pypgo/pgo_test_01.py
```

## Tools

### Cubic Mesher

`cubicMesher` converts a closed triangle surface mesh in `.obj` format into a cubic volumetric `.veg` mesh and can optionally export the extracted cubic surface as `.obj`.

Build the tool:

```bash
cmake --preset base
cmake --build --preset base --target cubicMesher
```

Basic usage:

```bash
build/base/bin/cubicMesher \
--input-mesh examples/ipc/cubic/box/box.obj \
--resolution 4 \
--output-mesh /tmp/libpgo-box.veg \
--output-surface /tmp/libpgo-box-surface.obj \
--E 10000000 \
--nu 0.45 \
--density 1000
```

Main arguments:

- `--input-mesh`: input closed triangle mesh in `.obj`
- `--resolution`: number of cubic cells along the shortest input AABB edge
- `--output-mesh`: output cubic `.veg`
- `--output-surface`: optional extracted surface `.obj`
- `--E`, `--nu`, `--density`: isotropic material parameters written into the output mesh

### Tet Mesher

`tetMesher` converts a closed triangle surface mesh into a tetrahedral `.veg` simulation mesh from a JSON job config. The JSON selects the backend, backend parameters, input/output paths, and optional generated boundary surface export. Paths inside the config are resolved relative to the config file.

Build `tetMesher` in the default no-MKL preset:

```bash
cmake --preset base
cmake --build --preset base --target tetMesher
```

Run a tet meshing job:

```bash
build/base/bin/tetMesher --config path/to/tetmesh.json
```

Basic TetGen config:

```json
{
  "version": 1,
  "backend": "tetgen",
  "input_mesh": "union_shell_remesh.obj",
  "output_mesh": "union_shell_tetgen.veg",
  "output_surface": "union_shell_tetgen_surface.obj",
  "print_stats": true,
  "tetgen": {
    "command": "pq1.414a0.01"
  }
}
```

The fTetWild backend is optional. Enable it when configuring, then build it in the preset build tree:

```bash
cmake --preset base -DPGO_TET_MESHER_USE_TET_WILD=ON
cmake --build --preset base --target tetMesher
```

Basic fTetWild config:

```json
{
  "version": 1,
  "backend": "tetwild",
  "input_mesh": "union_shell_remesh.obj",
  "output_mesh": "union_shell.veg",
  "output_surface": "union_shell_tet_surface.obj",
  "print_stats": true,
  "quiet": true,
  "tetwild": {
    "lr": 0.05,
    "epsr": 0.001,
    "stop_energy": 10,
    "max_threads": 8
  }
}
```

Config fields:

- `backend`: `tetgen` or `tetwild`
- `input_mesh`: input closed triangle mesh, typically `.obj`
- `output_mesh`: output tetrahedral `.veg`
- `output_surface`: optional generated tet boundary surface `.obj`
- `print_stats`, `quiet`: optional shared booleans
- `tetgen.command`: TetGen command string
- `tetwild.lr` / `tetwild.la`: relative or absolute fTetWild target edge length
- `tetwild.epsr`: fTetWild relative envelope tolerance
- `tetwild.stop_energy`, `tetwild.max_threads`: fTetWild optimization controls

## Third-party libraries

This library use the following third-party libraries:<br>
alembic, argparse, autodiff, boost, ceres, cgal, fmt, geogram, gmesh, json, knitro, libigl, mkl, spdlog, suitesparse, tbb, tinyobj-loader

---

## Licence

This library is developed using [VegaFEM](https://viterbi-web.usc.edu/~jbarbic/vega/) along with various third-party libraries, each governed by their respective licenses. Detailed copyright and license information is included within the majority of the source files.

In instances where specific licensing details are not provided within a source file, the copyright remains with the author. The licensing for those source files adhere to the principles of the pre-existing license framework. For instance, if a source file without licensing details incorporates components that fall under the GPL parts of CGAL, then that file will adhere to the GPL. All other source files default to the MIT License unless stated otherwise.

---

## TODO

- [x] Functional and compilable on three major platforms.
- [ ] Documentation
- [ ] More python interface
- [ ] Cleanup source code with non-MIT/non-FreeBSD licence.
- [ ] GUI
