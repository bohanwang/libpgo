## libpgo: Library for Physically based Simulation (P), Geometric Shape Modeling (G), and Optimization (O)

The library is designed to primarily focus on physically based simulations, geometric shape modeling, and optimization.
The source code extends [VegaFEM](https://viterbi-web.usc.edu/~jbarbic/vega/) and is designed for academic research purposes.

Release information:

- [Release notes](release-notes.txt)
- [Detailed 0.0.4 release notes](docs/release/0.0.4.md)

---

## Install a prebuilt wheel

Release 0.0.4 is distributed as standalone CPython 3.12 wheels through the
recorded [GitHub Actions artifacts](https://github.com/annajcy/libpgo/actions).
It is not published to PyPI, no source distribution is produced, and wheels
are not committed to this repository.

| Platform | Supported target | Artifact name |
| --- | --- | --- |
| Linux | `manylinux_2_28`, x86-64 | `pypgo-0.0.4-manylinux_2_28-x86_64` |
| macOS | macOS 26+, arm64 | `pypgo-0.0.4-macos-arm64` |
| Windows | Windows x86-64 | `pypgo-0.0.4-windows-x86_64` |

Download the artifact for the release commit. A prebuilt wheel does not
require Conda or a system installation of MKL, TBB, GMP, or MPFR.
Install [uv](https://docs.astral.sh/uv/getting-started/installation/) first;
`uv venv` downloads Python 3.12 when it is not already available.

Linux:

```bash
uv venv --python 3.12 .venv
uv pip install --python .venv/bin/python "numpy==2.0.2"
uv pip install --python .venv/bin/python --no-deps \
  /path/to/pypgo-0.0.4-cp312-cp312-manylinux_2_28_x86_64.whl
.venv/bin/python -c "import pypgo; print(pypgo)"
```

macOS 26 arm64:

```bash
uv venv --python 3.12 .venv
uv pip install --python .venv/bin/python "numpy==2.0.2"
uv pip install --python .venv/bin/python --no-deps \
  /path/to/pypgo-0.0.4-cp312-cp312-macosx_26_0_arm64.whl
.venv/bin/python -c "import pypgo; print(pypgo)"
```

Windows PowerShell:

```powershell
uv venv --python 3.12 .venv
uv pip install --python .venv\Scripts\python.exe "numpy==2.0.2"
uv pip install --python .venv\Scripts\python.exe --no-deps `
  C:\path\to\pypgo-0.0.4-cp312-cp312-win_amd64.whl
.\.venv\Scripts\python.exe -c "import pypgo; print(pypgo)"
```

Linux and Windows wheels bundle oneMKL, the matching oneTBB runtime, GMP, and
MPFR. The macOS wheel bundles oneTBB, GMP, and MPFR and links the system
Accelerate framework. See the
[prebuilt-wheel guide](docs/guide/build/build-from-wheel.md) for artifact
download and checksum commands.

## Build from source

Source builds require CMake 3.29 or newer, a C++20 compiler, Python 3.12, and
network access for pinned FetchContent archives. Unlike a repaired release
wheel, a source installation may depend on native libraries installed on the
build machine.

The repository tracks `.python-version` and a cross-platform `uv.lock`.
`uv sync --locked` installs the common build/test tools plus the current
platform's wheel-repair tool; Linux and Windows additionally receive the
matching Intel MKL/TBB packages through environment markers. It prepares the
environment but intentionally does not compile or install pypgo. The default
developer workflow builds the extension in `build/pypgo` and imports it
directly through `PYTHONPATH`, so an incremental C++ rebuild does not require
another `pip install`.

### Linux x86-64

Ubuntu or Debian:

```bash
sudo apt update
sudo apt install -y build-essential git libgmp-dev libmpfr-dev

uv sync --locked

uv run cmake --preset pypgo-wheel
uv run cmake --build build/pypgo
uv run ctest --test-dir build/pypgo --output-on-failure

export PYTHONPATH="$PWD/build/pypgo/src/python/pypgo"
export LD_LIBRARY_PATH="$PWD/.venv/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
uv run python -c "import pypgo; print(pypgo)"
uv run python -m pytest -q tests/pypgo/test_pgo_smoke.py
```

Use `gmp-devel` and `mpfr-devel` instead of `libgmp-dev` and `libmpfr-dev` on
Fedora/RHEL. Intel's `mkl-devel` package installs its matching `tbb-devel`
dependency into the same virtual environment. CMake discovers native
dependencies from the active Python environment and standard platform paths.

### macOS 26 arm64

```bash
xcode-select --install
brew install gmp mpfr tbb

uv sync --locked

uv run cmake --preset pypgo-wheel
uv run cmake --build build/pypgo
uv run ctest --test-dir build/pypgo --output-on-failure

export PYTHONPATH="$PWD/build/pypgo/src/python/pypgo"
uv run python -c "import pypgo; print(pypgo)"
uv run python -m pytest -q tests/pypgo/test_pgo_smoke.py
```

This source-built extension links the Homebrew TBB/GMP/MPFR libraries
directly. Keep those packages installed while using the environment.

### Windows x86-64

Install Python 3.12 and the Visual Studio 2022 **Desktop development with C++**
workload. Run the following from an x64 Native Tools PowerShell:

```powershell
uv sync --locked

uv run cmake --preset pypgo-wheel
uv run cmake --build build\pypgo
uv run ctest --test-dir build\pypgo --output-on-failure

$env:PYTHONPATH = "$PWD\build\pypgo\src\python\pypgo"
uv run python -c "import pypgo; print(pypgo)"
uv run python -m pytest -q tests\pypgo\test_pgo_smoke.py
```

Windows source builds use the approved GMP/MPFR files under `third-party`.
CMake stages the required GMP/MPFR, oneTBB, and oneMKL runtime DLLs beside
source-built executables and `pypgo.pyd`; no dependency-specific `PATH`
configuration is required.
After the first configure, ordinary C++ edits only require
`uv run cmake --build build/pypgo` (or `--target pypgo` when only the Python
module is needed). Detailed commands and the separate release-wheel workflow
are in the
[source-build guide](docs/guide/build/build-from-source.md).

CI separates these concerns into a source `build-test` job and a fresh
`package` job. Release automation runs
`scripts/release_wheel_provenance.py preflight` after downloading the native
test evidence. The command rejects a dirty checkout and records the exact
commit before `uv build --wheel --no-build-isolation` runs. After the repaired
wheel passes its installed smoke test and Python tests, the `record` command
requires exactly one wheel, rejects source archives, and records the wheel
checksum together with the CMake, dependency, test, and runner evidence.
Evidence and wheel output directories must be outside the source checkout.

## Usage & Test

We provide three Python scripts to exercise the installation. From a source
checkout, run them through the locked project environment created by
`uv sync --locked`.

1. `pgo_test_01.py`. It runs a few basic pgo APIs.

    ```bash
    cd examples
    uv run python ../src/python/pypgo/pgo_test_01.py
    ```

    The expected result will look like

    ```text
    Opening file torus.veg.
    #vtx:564
    #tets:1950
    164,134,506,563
    L Info:
    10067040
    (10067040,)
    (10067040,)
    125.0
    GTLTLG Info:
    503400
    (503400,)
    (503400,)
    9695578.0
    [[  6.958279    0.          0.        -17.495821    0.          0.
       13.10052     0.          0.         -2.5629783   0.          0.       ]
     [  0.          6.958279    0.          0.        -17.495821    0.
        0.         13.10052     0.          0.         -2.5629783   0.       ]
     [  0.          0.          6.958279    0.          0.        -17.495821
        0.          0.         13.10052     0.          0.         -2.5629783]
     [ -5.1109824   0.          0.         10.111505    0.          0.
        8.160282    0.          0.        -13.160804    0.          0.       ]
     [  0.         -5.1109824   0.          0.         10.111505    0.
        0.          8.160282    0.          0.        -13.160804    0.       ]
     [  0.          0.         -5.1109824   0.          0.         10.111505
        0.          0.          8.160282    0.          0.        -13.160804 ]
     [ 23.97409     0.          0.         -6.634346    0.          0.
       -1.4866991   0.          0.        -15.853046    0.          0.       ]
     [  0.         23.97409     0.          0.         -6.634346    0.
        0.         -1.4866991   0.          0.        -15.853046    0.       ]
     [  0.          0.         23.97409     0.          0.         -6.634346
        0.          0.         -1.4866991   0.          0.        -15.853046 ]]
    ```

2. `pgo_run_sim.py`. It reads input config file and run simulation. You can try `box`, `box-with-sphere`, `dragon`, and `dragon-dyn` to test different simulation results. Take the box example for illustration. You can run the box example using the following commands.
   
    ```bash
    uv run python src/python/pypgo/pgo_run_sim.py \
        examples/configs/volume/box/box-tet-sampled.json
    ```

    The expected result will look like the first image. The time integrator is hard-coded as implicit backward Euler (BE). You are free to change it to implicit Newmark (NW) or TR-BDF2 integrator (not support friction).
    <table style="width: 100%; table-layout: fixed; border-collapse: collapse;">
        <tr>
            <th style="width: 50%;text-align:center; border-top: 1px solid #ddd;">Box (NM)</th>
            <th style="width: 50%;text-align:center; border-top: 1px solid #ddd;">Box with Sphere (NM)</th>
        </tr>
        <tr>
            <td style="text-align: center; border-bottom: 1px solid #ddd;"><img src="./examples/assets/media/box.gif" alt="box"></td>
            <td style="text-align: center; border-bottom: 1px solid #ddd;"><img src="./examples/assets/media/box-with-sphere.gif" alt="box with sphere"></td>
        </tr>
        <tr>
            <th style="width: 50%;text-align:center;">Dragon (BE)</th>
            <th style="width: 50%;text-align:center;">Bunny (BE)</th>
        </tr>
        <tr>
            <td style="text-align: center; border-bottom: 1px solid #ddd;"><img src="./examples/assets/media/dragon-dynamic.gif" alt="dragon"></td>
            <td style="text-align: center; border-bottom: 1px solid #ddd;"><img src="./examples/assets/media/bunny.gif" alt="bunny"></td>
        </tr>
        <tr>
            <th style="width: 50%;text-align:center;">Rest Dragon</th>
            <th style="width: 50%;text-align:center;">Deformed Dragon</th>           
        </tr>
        <tr>
            <td style="text-align: center; border-bottom: 1px solid #ddd;"><img src="./examples/assets/media/dragon-rest.png" alt="dragon rest shape"></td>
            <td style="text-align: center; border-bottom: 1px solid #ddd;"><img src="./examples/assets/media/dragon-deformed.png" alt="dragon deformed shape"></td>
        </tr>
    </table>

3. `pgo_dump_abc.py`. It creates the abc file that can be used for blender/maya from config file `anim.json`. Essentially, it takes the simulation output `.obj` sequences and output a `.abc` file.

    ```bash
    uv run python src/python/pypgo/pgo_dump_abc.py \
        examples/configs/volume/box/box-tet-sampled-animation.json
    ```

    The `convertAnimation` tool provides the same conversion on the CLI:

    ```bash
        convertAnimation \
            examples/configs/volume/box/box-tet-sampled-animation.json
    ```

    If the optional second argument is omitted, the tool writes `.abc` files into the folder containing `anim.json`, and each output filename uses the mesh `name` field from the config.

## Tools

### Cubic Mesher

`cubicMesher` converts a closed triangle surface mesh in `.obj` format into a cubic volumetric `.veg` mesh and can optionally export the extracted cubic surface as `.obj`.

Build the tool:

```bash
    uv run cmake --preset pypgo-wheel
    uv run cmake --build build/pypgo --target cubicMesher
```

Generate the documented presets with the checked-in helper:

```bash
uv run python examples/scripts/generate_cubic_veg.py \
    --build-dir build/pypgo \
    --scene box
```

Main arguments:

- `--input-mesh`: input closed triangle mesh in `.obj`
- `--resolution`: number of cubic cells along the shortest input AABB edge
- `--output-mesh`: output cubic `.veg`
- `--output-surface`: optional extracted surface `.obj`
- `--E`, `--nu`, `--density`: isotropic material parameters written into the output mesh

Generated cubic meshes are written under the gitignored
`examples/generated/cubic/` directory. See
[`examples/README.md`](./examples/README.md) for presets, custom input options,
runner commands, and the old-to-new path table.

### Shell Simulation

Build the shell simulation CLI:

```bash
    uv run cmake --preset pypgo-wheel
    uv run cmake --build build/pypgo --target runShellSim
```

Run the bundled shell example:

```bash
    build/pypgo/bin/runShellSim \
        examples/configs/shell/shell-dynamic-sampled.json
```

To also write command-line output to a log next to the config, add `--log`:

```bash
    build/pypgo/bin/runShellSim \
        examples/configs/shell/shell-dynamic-sampled.json --log
```

---

## Build the C++ library without Python (optional)

After installing the platform dependencies from the source-build section,
configure a native tree instead of running `pip`.

Linux:

```bash
uv run cmake -S . -B build/native -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DPGO_ENABLE_FULL=ON
uv run cmake --build build/native
uv run ctest --test-dir build/native --output-on-failure
```

macOS 26 arm64:

```bash
uv run cmake -S . -B build/native -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DPGO_ENABLE_FULL=ON
uv run cmake --build build/native
uv run ctest --test-dir build/native --output-on-failure
```

Windows x64 Native Tools PowerShell:

```powershell
uv run cmake -S . -B build\native -G Ninja `
  -DCMAKE_BUILD_TYPE=Release `
  -DPGO_ENABLE_FULL=ON
uv run cmake --build build\native
uv run ctest --test-dir build\native --output-on-failure
```

---

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
