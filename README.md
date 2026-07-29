## libpgo: Library for Physically based Simulation (P), Geometric Shape Modeling (G), and Optimization (O)

The library is designed to primarily focus on physically based simulations, geometric shape modeling, and optimization.
The source code extends [VegaFEM](https://viterbi-web.usc.edu/~jbarbic/vega/) and is designed for academic research purposes.

---

## Python wheels

The 0.0.4 release supports CPython 3.12 on these targets:

- Linux x86-64 (`linux_x86_64`)
- macOS arm64
- Windows x86-64 (`win_amd64`)

Release wheels are downloaded from the recorded
[GitHub Actions artifacts](https://github.com/annajcy/libpgo/actions) and
installed from the local wheel file. They are designed for the documented
Conda/Miniforge runtime and are not generic wheels for an otherwise empty
virtual environment.

Version 0.0.4 is distributed as wheel artifacts only. It is not published on
PyPI, no source distribution is produced, and release wheels are not committed
to this repository.

Pip cannot install Conda packages. Before installing a wheel, create a runtime
environment containing NumPy, GMP, MPFR, Imath, fmt, and TBB:

```bash
conda create -n pypgo-004 -c conda-forge \
  python=3.12 "numpy>=1.26" gmp mpfr imath fmt tbb tbb-devel
conda activate pypgo-004
conda config --env --set channel_priority strict
```

On Linux and Windows, add the MKL runtime:

```bash
conda install -c conda-forge mkl mkl-devel \
  "libblas=*=*mkl" "liblapack=*=*mkl"
```

On macOS, libpgo links the system Accelerate framework. Select the matching
Accelerate family for NumPy:

```bash
conda install -c conda-forge \
  "libblas=*=*newaccelerate" "liblapack=*=*newaccelerate"
```

Then install the downloaded artifact:

```bash
python -m pip install /path/to/pypgo-0.0.4-<python>-<abi>-<platform>.whl
```

GMP, MPFR, Imath, fmt, MKL, and TBB remain external runtime libraries and are
not bundled into the wheel.

## Building from source

Source builds require CMake 3.28 or newer, a C++20 compiler, Python 3.12,
NumPy, and the native dependencies listed above. GCC 11–13, Apple Clang, and
Visual Studio 2022 are the currently supported compiler families.

Install CMake and Ninja in the active Conda environment, then build with the
existing setup entry point:

```bash
conda install -c conda-forge cmake ninja
python -m pip install .
```

The default build behavior remains platform-compatible with 0.0.3. Release
automation can explicitly select `pypgo-wheel-no-mkl` or `pypgo-wheel-mkl`
through `PYPGO_CMAKE_PRESET`. Set `CMAKE_BUILD_PARALLEL_LEVEL` to control
parallel compilation.

Release automation runs `scripts/release_wheel_provenance.py preflight` after
native and Python tests pass. The command rejects a dirty checkout and records
the exact commit before `python -m build --wheel --no-isolation` runs. Its
`record` command then requires exactly one wheel, rejects source archives, and
records the wheel checksum together with the CMake, dependency, test, and
runner evidence. Evidence and wheel output directories must be outside the
source checkout.

## Usage & Test

We provide three python scripts to test the installation.

1. `pgo_test_01.py`. It runs a few basic pgo APIs.

    ```bash
        cd examples
        python ../src/python/pypgo/pgo_test_01.py
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
        python src/python/pypgo/pgo_run_sim.py \
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
        python src/python/pypgo/pgo_dump_abc.py \
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
    cmake --preset base_no_mkl
    cmake --build build/base_no_mkl --target cubicMesher
```

Generate the documented presets with the checked-in helper:

```bash
    python3 examples/scripts/generate_cubic_veg.py \
        --build-dir build/base_no_mkl \
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
    cmake --preset base_no_mkl
    cmake --build build/base_no_mkl --target runShellSim
```

Run the bundled shell example:

```bash
    build/base_no_mkl/bin/runShellSim \
        examples/configs/shell/shell-dynamic-sampled.json
```

To also write command-line output to a log next to the config, add `--log`:

```bash
    build/base_no_mkl/bin/runShellSim \
        examples/configs/shell/shell-dynamic-sampled.json --log
```

---

## Setup without Python (Optional)

If you want to use the library with your C++ code or modify the source code, you may build it without python.

### Windows & Ubuntu

To compile the lib with a basic functionality,

```bash
    cd libpgo
    mkdir build
    cd build
    cmake ..
```

To enable all functionalities, Install [MKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html). Then,

```bash
    cd libpgo
    mkdir build
    cd build
    cmake .. -DPGO_USE_MKL=1 -DPGO_ENABLE_FULL=1
```

> On Windows, a few extra steps are need before running the above commands. First, the library should be configured in "x64 Native Tools Command Prompt for VS 2022". In addition, before running the commands above, run `c:\Program Files (x86)\Intel\oneAPI\setvars.bat` to setup the environments for MKL, where `c:\Program Files (x86)\Intel\oneAPI` is the path to the oneAPI installation. Once setup, run above commands.

> On Ubuntu, a similar procedure is needed. Before configuring the library, run `bash /opt/intel/oneapi/setvars.sh` to setup the MKL environments for the subsequent cmake configuration.

### Mac OS

To have a basic functionality, use CMake to compile it like on Windows & Ubuntu.

To enable all functionalities,

```bash
    cd libpgo
    mkdir build
    cd build
    cmake .. -DPGO_ENABLE_FULL=1 -DDPGO_ENABLE_ALEMBIC=1 -DPGO_ENABLE_GMSH=1
```
The last two flags work only if you have imath and gmesh libs.

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
