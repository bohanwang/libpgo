## libpgo: Library for Physically based Simulation (P), Geometric Shape Modeling (G), and Optimization (O)

The library is designed to primarily focus on physically based simulations, geometric shape modeling, and optimization.
The source code extends [VegaFEM](https://viterbi-web.usc.edu/~jbarbic/vega/) and is designed for academic research purposes.

---

## Prerequisites (Important!)

1. CMake >= **3.28**\
    We use several functionalities that are only supported by 3.28+. 
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
        brew install gmp mpfr
    ```

4. (Optional) Ninja\
    It can be installed by

    ```bash
        pip install ninja
    ```

    for better compilation performance

5. (Optional) numpy\
    This is used for running tests.

## Compilation

Going forward, it is assumed that all specified prerequisites are installed and that a Conda environment is used for python.

### Windows & Ubuntu

Install prerequisites:

```bash
    conda install tbb tbb-devel mkl mkl-devel
```

Install libpgo:

```bash
    cd libpgo
    pip install .
```

If `ninja` has been installed, it will compile source files in parallel. If it is not installed,
set `CMAKE_BUILD_PARALLEL_LEVEL` to `n`, where `n` is the number of threads for compilation, to control the parallel compilation.

### Mac OS

Install libpgo

```bash
    cd libpgo
    pip install .
```

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

2. `pgo_test_02.py`. It runs a quasi-static simulation for a dragon under attachments and gravity.

    ```bash
        cd examples
        python ../src/python/pypgo/pgo_test_02.py
    ```

    The expected result will look like the right image.

    | Rest Shape    | Deformed Shape |
    | ------------- | -------------- |
    | ![Rest!](/examples/dragon-img/rest.png) | ![Deformed!](/examples/dragon-img/deformed.png) |

3. `pgo_test_03.py`. It runs a dynamic simulation with contact handling.

    ```bash
        cd examples
        python ../src/python/pypgo/pgo_test_03.py
    ```

    In this scene, a box is dropped on the curved ground. The running speed of this demo is slow, so you may run it only if you are interested in the results. The result is a sequence of `.obj` files, each of which is the deformed shape of ecah timestep. These files are saved in `ret-box` folder. The result will look like:
    |   Box    |      |
    | -------| ---- |
    | ![box!](/examples/box-img/output.gif) | ![white!](/examples/box-img/white.png) |

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
    cmake .. -DPGO_ENABLE_FULL=1
```

---

## Third-party libraries

This library use the following third-party libraries:<br>
autodiff, boost, ceres, cgal, fmt, knitro, mkl, json, pybind11, spdlog, suitesparse, tbb, tinyobj-loader

---

## Licence

This library is developed using [VegaFEM](https://viterbi-web.usc.edu/~jbarbic/vega/) along with various third-party libraries, each governed by their respective licenses. Detailed copyright and license information is included within the majority of the source files.

In instances where specific licensing details are not provided within a source file, the copyright remains with the author. The licensing for those source files adhere to the principles of the pre-existing license framework. For instance, if a source file without licensing details incorporates components that fall under the GPL parts of CGAL, then that file will adhere to the GPL. All other source files default to the MIT License unless stated otherwise.

---

## TODO

- [x] Functional and compilable on three major platforms.
- [ ] Support cmake of lower version
- [ ] Documentation
- [ ] More python interface
- [ ] Cleanup source code with non-MIT/non-FreeBSD licence.
- [ ] GUI
