# Architecture Overview

libpgo is organized into three top-level source directories under `src/`.

```text
src/
├── core/          # Core library modules
│   ├── energy/    # Energy formulations (elasticity, contact, etc.)
│   ├── external/  # Vendored or adapted third-party code
│   ├── scene/     # Scene graph, mesh I/O, and simulation setup
│   ├── solve/     # Solvers and time integrators
│   └── utils/     # Shared utilities (math, containers, etc.)
├── api/           # Public API layer (Python bindings)
└── tools/         # Standalone command-line tools
    ├── animation/ # Animation export utilities
    ├── cubicMesher/   # Cubic/hexahedral mesh generator
    ├── remeshSurface/ # Surface remeshing tool
    ├── runSim/        # Simulation runner
    └── tetMesher/     # Tetrahedral mesh generator
```

## Build System

libpgo uses **CMake** with **Conan 2.x** for native dependency resolution and
**scikit-build-core** for Python packaging.

See the [Build System Documentation](../build-system.md) for the canonical
explanation of presets, feature normalization, Conan bootstrap, generated
toolchain files, and Python packaging flow.

## Feature Model

The build is feature-oriented. See the [C++ Library Guide](../guide/cpp-library.md) for the full option table.

Key features: Python bindings, animation I/O, geometry stack (Boost/CGAL/Ceres/NLopt/SuiteSparse), libigl, geogram, gmsh, MKL, ARPACK.
