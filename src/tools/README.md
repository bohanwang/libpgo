# libpgo command-line tools

This directory contains the standalone C++ executables shipped by a full
libpgo source build. Installed Python users should normally prefer the `pypgo`
console commands for simulation and animation conversion:

```bash
pgo-run-sim <simulation-config.json>
pgo-run-cases <case-manifest.json> [CASE ...]
pgo-dump-abc <animation-config.json> <output-directory>
```

The C++ tools remain useful for native workflows, mesh preparation, and direct
access to individual backends.

## Build

Configure a full native build, then build every available tool:

```bash
uv run cmake --preset base
uv run cmake --build --preset base
```

The `pypgo-wheel` preset is intended for the portable Python wheel. It enables
Python bindings and selected optional dependencies, but leaves
`PGO_ENABLE_FULL` disabled, so it is not the canonical build for the complete
tool set documented here.

Build one tool by target name when only a single executable is needed:

```bash
uv run cmake --build build/base --target cubicMesher
```

Executables are written to `build/base/bin/` for the documented preset.

## Tool index

| Category | Tool | Purpose |
| --- | --- | --- |
| Simulation | `runSim` | Run sampled-contact tet or cubic simulations |
| Simulation | `runIPCSim` | Run IPC tet, cubic, or shell simulations |
| Animation | `convertAnimation` | Convert configured mesh sequences to Alembic |
| Cubic meshing | `cubicMesher` | Voxelize a closed surface into a cubic VEG mesh |
| Tet meshing | `tetMesher` | Tetrahedralize a closed surface through TetGen |
| Tet utilities | `generateTetMeshSurfaceMesh` | Extract a tet VEG surface as OBJ |
| Tet utilities | `mshFileToVegFile` | Convert a Gmsh MSH file to VEG |
| Surface processing | `remeshSurface` | Smooth, simplify, or remesh a surface |
| Surface processing | `mergeCloseVertices` | Merge nearby surface vertices and repair faces |
| Surface processing | `removeIsolatedVertices` | Remove vertices unused by any triangle |

`cubicMesherCore` is an internal static library used by `cubicMesher`, not a
command-line executable.

## Simulation

### `runSim`

Runs the sampled-contact simulation path. Tet and cubic volume configs are
supported. Static sampled configs run elasticity, gravity, and attachments,
but do not assemble sampled contact into the static Newton energy.

```bash
build/base/bin/runSim [--log] <simulation-config.json>
```

`--log` writes command output to a log file beside the config.

### `runIPCSim`

Runs the IPC simulation path for tet, cubic, and shell meshes in dynamic or
static mode.

```bash
build/base/bin/runIPCSim [--log] <simulation-config.json>
```

`--log` writes command output to a log file beside the config. New workflows
can use `pgo-run-sim` instead; it dispatches from the config's `contact-model`.

## Animation

### `convertAnimation`

Loads an animation JSON file and writes one Alembic file per configured mesh.
If the output directory is omitted, files are written beside the animation
config. The destination directory must already exist.

```bash
build/base/bin/convertAnimation <animation-config.json> [output-directory]
```

The installed Python equivalent is `pgo-dump-abc`.

## Cubic meshing

### `cubicMesher`

Voxelizes a closed triangle OBJ into an axis-aligned cubic VEG mesh. Material
parameters are embedded in the VEG file. The extracted cubic boundary surface
can be written at the same time.

```bash
build/base/bin/cubicMesher \
  --input-mesh input.obj \
  --resolution 20 \
  --output-mesh output.veg \
  --output-surface output-surface.obj \
  --E 100000 \
  --nu 0.45 \
  --density 1000
```

`--resolution` is the number of cubic cells along the shortest input bounding
box axis. `--output-surface` is optional.

The checked-in helper wraps this tool, validates both outputs, and records
hashes and mesh statistics in a manifest:

```bash
uv run python examples/scripts/generate_cubic_veg.py \
  --build-dir build/base
```

Its default output directory is `examples/assets/generated/cubic/`.

## Tet meshing and conversion

### `tetMesher`

Tetrahedralizes a closed triangle surface through TetGen and writes a tet VEG
mesh. The `--command` value is passed to TetGen as its switch string.

```bash
build/base/bin/tetMesher tetgen \
  --input-mesh input.obj \
  --output-mesh output.veg \
  --command pq1.2
```

### `generateTetMeshSurfaceMesh`

Extracts the boundary triangles of a tet VEG mesh as OBJ:

```bash
build/base/bin/generateTetMeshSurfaceMesh input.veg output.obj
```

An optional third argument beginning with `A` writes all element faces instead
of only the exterior boundary:

```bash
build/base/bin/generateTetMeshSurfaceMesh input.veg output.obj A
```

This legacy positional tool does not implement `--help`; invoke it only with
the arguments shown above.

### `mshFileToVegFile`

Converts a Gmsh MSH volume mesh to tet VEG:

```bash
build/base/bin/mshFileToVegFile input.msh output.veg
```

This target is available only when the build provides `Gmsh::Gmsh`.

## Surface processing

### `remeshSurface`

Provides four backends through subcommands:

```bash
# Isotropic CGAL remeshing. The edge-length value multiplies the input's
# average edge length.
build/base/bin/remeshSurface cgal_iso \
  --input-mesh input.obj --output-mesh output.obj \
  --edge-length 0.5 [--sharp-edge-angle 180]

# CGAL edge-collapse simplification.
build/base/bin/remeshSurface cgal_simplify \
  --input-mesh input.obj --output-mesh output.obj \
  --target-ratio 0.5

# CGAL angle-and-area smoothing.
build/base/bin/remeshSurface cgal_smooth \
  --input-mesh input.obj --output-mesh output.obj \
  [--sharp-edge-angle 180]

# Geogram remeshing followed by projection to the input surface.
build/base/bin/remeshSurface geogram \
  --input-mesh input.obj --output-mesh output.obj \
  --target-num-vertices 1000
```

The `remeshSurface` target is omitted when CGAL or Geogram support is missing.
Run `remeshSurface <subcommand> --help` for backend-specific options.

### `mergeCloseVertices`

Merges vertices within a distance tolerance, drops degenerate triangles, and
runs CGAL surface repair. If `epsilon` is omitted, the tool uses half of the
minimum input edge length.

```bash
build/base/bin/mergeCloseVertices input.obj output.obj [epsilon]
```

The input and output formats are those supported by the configured CGAL polygon
mesh I/O implementation.

### `removeIsolatedVertices`

Removes vertices that are not referenced by any triangle:

```bash
build/base/bin/removeIsolatedVertices input.obj output.obj
```

## Conditional availability

The exact target set depends on configured optional libraries:

- `convertAnimation` requires the `animationIO` target;
- `remeshSurface`, `mergeCloseVertices`, and `removeIsolatedVertices` require
  both CGAL and Geogram integration;
- `mshFileToVegFile` requires the imported `Gmsh::Gmsh` target.

Use the generated build target list when scripting across configurations:

```bash
uv run cmake --build build/base --target help
```
