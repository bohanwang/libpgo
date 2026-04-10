# C++ Library Guide

For the full build internals, see the [Build System Documentation](../build-system.md).

## Feature Options

| Option | Default | Description |
|---|---|---|
| `PGO_FEATURE_PYTHON` | `OFF` | Build Python bindings |
| `PGO_FEATURE_ANIMATION_IO` | `OFF` | Enable Alembic/Imath animation I/O |
| `PGO_FEATURE_GEOMETRY_STACK` | `OFF` | Enable Boost/CGAL/Ceres/NLopt/SuiteSparse |
| `PGO_FEATURE_LIBIGL` | `OFF` | Enable libigl integration |
| `PGO_FEATURE_GEOGRAM` | `OFF` | Enable geogram integration |
| `PGO_FEATURE_GMSH` | `OFF` | Enable gmsh integration |
| `PGO_FEATURE_MKL` | `OFF` | Enable MKL support |
| `PGO_FEATURE_ARPACK` | `OFF` | Enable ARPACK support |
| `PGO_PROFILE_FULL` | `OFF` | Convenience: enables Python + animation I/O + geometry stack + libigl + geogram + gmsh |

Behavior:

- `PGO_FEATURE_PYTHON=ON` implies `PGO_FEATURE_GEOMETRY_STACK=ON`
- `PGO_PROFILE_FULL=ON` expands to the full feature set

## Build Presets

### Core Build

```bash
cmake --preset core-release
cmake --build --preset core-release
ctest --preset core-release
```

### Geometry Stack Build

```bash
cmake --preset geometry-release
cmake --build --preset geometry-release
ctest --preset geometry-release
```

### Animation I/O (without Python)

```bash
cmake --preset animation-release
cmake --build --preset animation-release
```

### Full Build

```bash
cmake --preset full-release
cmake --build --preset full-release
ctest --preset full-release
```

## Conan Dependency Management

The build uses `cmake/BootstrapConan.cmake` which automatically:

1. Exports private Conan recipes from `conan/recipes/`
2. Runs `conan install` during the CMake configure step
3. Tracks a feature signature to re-run Conan only when options change

Conan options mirror the CMake feature model: `with_python`, `with_animation_io`, `with_geometry_stack`, etc.

See the [Build System Documentation](../build-system.md) for the configure-time
control flow, preset mapping, and generated Conan metadata layout.

## MKL Support

MKL is system-managed (not Conan-managed). Install via Intel oneAPI or Conda, then enable `PGO_FEATURE_MKL=ON`.

- **Linux**: `source /opt/intel/oneapi/setvars.sh` before configure
- **Windows**: run `setvars.bat` from the oneAPI install directory
