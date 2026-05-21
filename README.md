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

3. GMP and MPFR for **Ubuntu** and **Mac OS**\
    This can be installed on Ubuntu by

    ```bash
    sudo apt install libgmp-dev libmpfr-dev
    ```

    Or it can be installed on Mac OS by

    ```bash
    brew install gmp mpfr imath
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

Going forward, it is assumed that all specified prerequisites are installed and that a Conda environment is used for python.

Install prerequisites:

```bash
conda install tbb tbb-devel mkl mkl-devel
conda install conda-forge::imath
```

### CMake Presets

Native C++ builds in this repository are preset-driven:

```bash
cmake --list-presets
cmake --preset <configure-preset>
cmake --build --preset <build-preset> [--target <target>...]
```

Configure presets define feature flags and build directories. Build presets define parallel build options and map to a configure preset.

#### Configure Presets

| Configure preset | Binary directory | Purpose / key options |
| --- | --- | --- |
| `base_no_mkl` | `build/base_no_mkl` | Release baseline without MKL. Full stack on: `PGO_ENABLE_FULL=ON`, Alembic/Gmsh/TetWild enabled. |
| `base` | `build/base` | Release baseline with MKL (`PGO_USE_MKL=ON`) and full stack enabled (non-macOS). |
| `base_win` | `build/base_win` | Windows-oriented release baseline: MKL on, Alembic/Gmsh/TetWild off. |
| `debug` | *(fragment preset)* | Inheritance fragment that sets `CMAKE_BUILD_TYPE=Debug`. |
| `knitro` | *(fragment preset)* | Inheritance fragment enabling Knitro (`PGO_OPT_USE_KNITRO=ON`) with `KNITRO_LIBRARY_HINT`. |
| `pardiso` | *(fragment preset)* | Inheritance fragment enabling original Pardiso (`PGO_HAS_ORIG_PARDISO=ON`) with `PARDISO_LIBRARY_HINT`. |
| `cuda` | *(fragment preset)* | Inheritance fragment enabling CUDA (`PGO_ENABLE_CUDA=ON`). |
| `base_knitro` | `build/base_knitro` | `base` + `knitro` (non-macOS). |
| `base_knitro_cuda` | `build/base_knitro_cuda` | `base` + `knitro` + `cuda` (non-macOS). |
| `all_debug` | `build/all_debug` | `base` + `knitro` + `pardiso` + `cuda` in Debug mode (non-macOS). |
| `all_release` | `build/all_release` | `base` + `knitro` + `pardiso` + `cuda` in Release mode (non-macOS). |
| `base_cuda_debug` | `build/base_cuda_debug` | `base` + `cuda` in Debug mode (non-macOS). |
| `base_no_mkl_debug` | `build/base_no_mkl_debug` | `base_no_mkl` in Debug mode. |
| `base_cuda_release` | `build/base_cuda_release` | `base` + `cuda` in Release mode (non-macOS). |
| `base_cuda_win` | `build/base_cuda_win` | `base_win` + `cuda`, with Windows `cudss_DIR` hint. |

#### Build Presets

| Build preset | Configure preset | Typical use |
| --- | --- | --- |
| `base` | `base` | Release build with MKL/full stack (non-macOS). |
| `all_debug` | `all_debug` | Debug build with all optional solvers/features (non-macOS). |
| `all_release` | `all_release` | Release build with all optional solvers/features (non-macOS). |
| `base_cuda_debug` | `base_cuda_debug` | Debug build with CUDA (non-macOS). |
| `base_no_mkl_debug` | `base_no_mkl_debug` | Debug build without MKL. |
| `base_no_mkl_release` | `base_no_mkl` | Release build without MKL. |
| `base_cuda_release` | `base_cuda_release` | Release build with CUDA (non-macOS). |

Some configure presets are composition-oriented and currently have no dedicated build preset (for example: `base_win`, `base_cuda_win`, `base_knitro`, `base_knitro_cuda`).

On macOS, use the no-MKL preset family only: `base_no_mkl`, `base_no_mkl_release`, and `base_no_mkl_debug`.

### Install libpgo

```bash
cd libpgo
pip install .
```

If `ninja` has been installed, it will compile source files in parallel. If it is not installed,
set `CMAKE_BUILD_PARALLEL_LEVEL` to `n`, where `n` is the number of threads for compilation, to control the parallel compilation.

### Setup without Python

If you want to use the library with your C++ code or modify the source code, you may build it without python.

### Windows & Ubuntu

To compile the lib with basic functionality (no MKL):

```bash
cd libpgo
cmake --preset base_no_mkl
cmake --build --preset base_no_mkl_release
```

To enable the MKL/full-feature stack, install [MKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html). Then,

```bash
cd libpgo
cmake --preset base
cmake --build --preset base
```

> On Windows, a few extra steps are need before running the preset commands above. First, the library should be configured in "x64 Native Tools Command Prompt for VS 2022". In addition, before running the commands above, run `c:\Program Files (x86)\Intel\oneAPI\setvars.bat` to setup the environments for MKL, where `c:\Program Files (x86)\Intel\oneAPI` is the path to the oneAPI installation. Once setup, run above commands.

> On Ubuntu, a similar procedure is needed. Before configuring the library with presets, run `bash /opt/intel/oneapi/setvars.sh` to setup the MKL environments for the subsequent CMake configuration.

### Mac OS

On macOS, MKL is not supported. Use only the no-MKL presets.

Release build:

```bash
cd libpgo
cmake --preset base_no_mkl
cmake --build --preset base_no_mkl_release
```

Debug build:

```bash
cmake --preset base_no_mkl_debug
cmake --build --preset base_no_mkl_debug
```

The `base_no_mkl` preset already keeps the full non-MKL feature stack enabled (including Alembic/Gmsh/TetWild). Alembic and Gmsh related features still depend on local third-party libraries (such as imath and gmsh).

---

## Usage & Test

The primary runnable examples in this repository are now IPC examples driven by `runIPCSim` under `examples/ipc/`.

Build the IPC tools:

```bash
cmake --preset base_no_mkl
cmake --build --preset base_no_mkl_release --target runIPCSim convertAnimation
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
build/base_no_mkl/bin/runIPCSim examples/ipc/shell/shell-hang/shell-ipc.json
build/base_no_mkl/bin/runIPCSim examples/ipc/shell/shell-drop/shell-ipc.json
build/base_no_mkl/bin/runIPCSim examples/ipc/tet/box-hang/box-ipc.json
build/base_no_mkl/bin/runIPCSim examples/ipc/cubic/box-with-sphere/box-ipc.json
```

Convert dumped frame sequences to Alembic:

```bash
build/base_no_mkl/bin/convertAnimation examples/ipc/shell/shell-hang/anim.json
build/base_no_mkl/bin/convertAnimation examples/ipc/shell/shell-drop/anim.json
build/base_no_mkl/bin/convertAnimation examples/ipc/tet/box-hang/anim.json
build/base_no_mkl/bin/convertAnimation examples/ipc/cubic/box-with-sphere/anim.json
```

For the full IPC case list and per-case notes, see [`examples/ipc/README.md`](./examples/ipc/README.md).

Legacy penalty-based volume contact is available through the same entrypoint:

```bash
build/base_no_mkl/bin/runIPCSim --legacy path/to/legacy-volume-config.json
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
cmake --preset base_no_mkl
cmake --build --preset base_no_mkl_release --target cubicMesher
```

Basic usage:

```bash
build/base_no_mkl/bin/cubicMesher \
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
cmake --preset base_no_mkl
cmake --build --preset base_no_mkl_release --target tetMesher
```

Run a tet meshing job:

```bash
build/base_no_mkl/bin/tetMesher --config path/to/tetmesh.json
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

The fTetWild backend is enabled by default in the main presets on macOS/Linux. Build it in the preset build tree:

```bash
cmake --preset base_no_mkl
cmake --build --preset base_no_mkl_release --target tetMesher
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
alembic, argparse, autodiff, boost, ceres, cgal, fmt, geogram, gmesh, json, knitro, libigl, mkl, pybind11, spdlog, suitesparse, tbb, tinyobj-loader

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
