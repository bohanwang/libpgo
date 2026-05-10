# FBMS Asset Pipeline

This directory contains the JSON-driven FBMS asset pipeline and packaged
simulation assets. The current production run is:

- `fbms_r512`: OpenVDB raw extraction at resolution 512, CGAL isotropic remesh,
  TetWild tetrahedralization, and post-TetWild repaired boundary surfaces when
  needed.

The generated assets are meant to provide one FBMS simulation asset and five
volume-matched TPMS baseline assets for each `g0_b*` case.

## Pipeline Overview

The pipeline starts from prepared shell meshes, generates raw union surfaces,
regularizes those surfaces, tetrahedralizes them, validates the TetWild
boundary, repairs failed boundary surfaces if needed, and finally packages only
the simulation-ready `.veg` and final `.veg.obj` files.

```text
Source FBMS OBJ files
        |
        v
Prepare assets
  - copy FBMS shells
  - generate/copy bounding spheres
  - generate/copy TPMS baseline shells
        |
        v
+-----------------------------+
| FBMS branch                 |
| fixed fbms thickness = 0.02 |
+-----------------------------+
        |
        v
Generate FBMS raw union surfaces
  generateFBMSUnionSurface openvdb
        |
        v
Raw quality gate
  - 3 components
  - invalid triangle check
  - enclosed_volume recorded
        |
        +-------------------------------+
        |                               |
        v                               v
Use enclosed_volume as TPMS budget      CGAL isotropic remesh
        |                               target edge ~= 0.02
        v                               |
+-----------------------------+         |
| TPMS baseline branch        |         |
| volume matched thickness    |         |
+-----------------------------+         |
        |                               |
        v                               |
Generate TPMS raw union surfaces        |
  generateFBMSUnionSurface openvdb      |
  --fbms-thickness = -volume            |
        |                               |
        v                               |
CGAL isotropic remesh                   |
  target edge ~= 0.02                   |
        |                               |
        +---------------+---------------+
                        |
                        v
Remesh quality check
  - full surface check
  - self-intersections are warnings
                        |
                        v
TetWild tetrahedralization
  tetMesher, lr = 0.008, epsr = 0.001
                        |
                        v
Final boundary quality check
  meshQualityCheck full
                        |
          +-------------+-------------+
          |                           |
          v                           v
   boundary passed             boundary failed
          |                           |
          v                           v
Use original .veg.obj          Repair boundary surface
                               repair_tpms_mesh.py
                               PyMeshFix per component
                                      |
                                      v
                               Repaired quality check
                                      |
                         +------------+------------+
                         |                         |
                         v                         v
                  repaired passed           repaired failed
                         |                         |
                         v                         v
        Use .veg.repaired.obj        Keep asset failed/not packageable
              as final surface
                         |
                         v
Dump 3 components and write volume info
                         |
                         v
Package simulation assets
  - .veg from TetWild
  - final .veg.obj from original or repaired surface
                         |
                         v
simulation_package/<case>/<asset>
```

For each `g0_b*`, the FBMS branch uses fixed shell thickness. The TPMS baseline
branch uses volume matching against that FBMS raw union volume, so every case
gets one FBMS asset plus five comparable TPMS baseline assets.

## Directory Layout

`fbms_r512/asset/`

Prepared source assets used by the pipeline. It contains:

- `asset/fbms/{case}.obj`
- `asset/fbms/{case}_bounding_sphere.obj`
- `asset/baseline/{shape}_fbms.obj`
- `asset/baseline/{shape}_fbms_bounding_sphere.obj`

`fbms_r512/raw_mesh_vdb_r512_t002/`

Raw FBMS union surfaces generated with fixed FBMS shell thickness `0.02`.

`fbms_r512/raw_mesh_vdb_r512_vmatch/`

Raw TPMS baseline union surfaces generated in volume-match mode. For each
`g0_b*`, the target volume is the corresponding FBMS raw union enclosed volume.

`fbms_r512/remesh_cgal_iso_r512_t002_e002/` and
`fbms_r512/remesh_cgal_iso_r512_vmatch_e002/`

CGAL isotropic remeshed surfaces. The physical target edge length is `0.02`;
the script converts that to the relative `remeshSurface cgal_iso --edge-length`
argument using `0.02 / raw_quality.edge_length.mean`.

`fbms_r512/tetmesh_tetwild_r512_t002_lr0008/` and
`fbms_r512/tetmesh_tetwild_r512_vmatch_lr0008/`

TetWild `.veg` outputs, exported boundary `.veg.obj`, quality JSONs, logs, and
optional repaired boundary surfaces.

`fbms_r512/simulation_package/`

The final packaged simulation assets. This is the directory to consume from
simulation code. Each case contains:

```text
simulation_package/
  g0_b5/
    fbms/
      g0_b5.veg
      g0_b5.veg.obj
    baseline/
      tpms_schwarz_p/
        tpms_schwarz_p.veg
        tpms_schwarz_p.veg.obj
      ...
  package_manifest.json
```

If a boundary was repaired, the package still exposes it as `{name}.veg.obj`;
the manifest records that its source was `{name}.veg.repaired.obj`.

## Pipeline Parameters

The r512 configuration lives in:

```bash
examples/fbms/fbms_r512/pipeline_r512.json
```

Current locked parameters:

- Raw extraction: `generateFBMSUnionSurface openvdb`
- Raw resolution: `512`
- FBMS thickness: `0.02`
- Sphere thickness: `0.03`
- Padding ratio: `0.08`
- Truncation, sphere-boundary projection, and small-component filtering enabled
- Remesh: `remeshSurface cgal_iso`
- Remesh physical target edge length: `0.02`
- Remesh sharp edge angle: `180`
- Remesh iterations: `10`
- TetWild: `lr=0.008`, `epsr=0.001`
- Material: density `1000`, Young's modulus `10000000`, Poisson ratio `0.45`

The configured TPMS baselines are:

- `tpms_schwarz_p`
- `tpms_schwarz_d`
- `tpms_gyroid`
- `tpms_iwp`
- `tpms_neovius`

## Tool Overview

The pipeline is a thin orchestration layer around several local tools.

`generate_bounding_sphere.py`

Creates a unit sphere OBJ at the origin (`--bounds-method unit`). Both FBMS
and TPMS use the same unit sphere so that raw union surfaces are comparable
across cases.

`generate_tpms_unit_ball.py`

Generates the five TPMS baseline shell meshes and matching unit spheres. With
`--flat-output`, it writes files directly as `{shape}_fbms.obj` and
`{shape}_fbms_bounding_sphere.obj`. The TPMS shells are sampled in `[-1,1]^3`.

`generateFBMSUnionSurface openvdb`

Builds a CSG union between a thickened shell mesh and a thickened sphere shell
using OpenVDB. Positive `--fbms-thickness` selects fixed shell thickness.
Negative `--fbms-thickness -V` triggers volume-budget mode: the tool binary-
searches for a shell thickness whose raw union volume is close to `V`.
Truncation clips the FBMS structure to the sphere interior. The raw output is
filtered to keep only the 3 largest edge-connected components
(`--filter-small-components --keep-largest-components 3`).

`meshQualityCheck surface`

Checks triangle surface quality and writes a JSON report. Raw checks are used
as hard gates for generated raw union meshes. Full checks are used for remesh
and final TetWild boundary surfaces. The JSON records component count,
invalid triangles, manifold/closed/winding status, self-intersections, edge
length statistics, and enclosed volume when meaningful.

`remeshSurface cgal_iso` / `remeshSurface geogram`

CGAL isotropic remeshing is the default. For noshell targets, if the CGAL
remesh has self-intersections the runner automatically falls back to geogram
at the same vertex count. The `--edge-length` argument to CGAL is relative to
the raw mesh mean edge length, so the runner computes a scale from the desired
physical target edge length.

`tetMesher`

Reads a `tetmesher.json` config, runs the TetWild backend, writes a `.veg`
volume mesh, and exports a `.veg.obj` boundary surface.

`volumetricMeshInfo`

Reads a `.veg` mesh and prints vertex count, element count, and total volume.
The runner saves this as `{name}.veg.info.txt`.

`repair_tpms_mesh.py`

Repairs a failed TetWild boundary OBJ using PyMeshFix per connected component.
It never overwrites the original `.veg` or `.veg.obj`; repaired surfaces are
written as `{name}.veg.repaired.obj`. The runner compares mtime of `veg.obj`
and `repaired.obj` to avoid reusing a stale repair from a previous run.

`dump_obj_components.py`

Splits the final boundary surface into connected components and writes a
`components.json` manifest plus per-component OBJ files. The pipeline requires
three final components.

`run_fbms_noshell_from_union.py`

Runs an independent no-shell pipeline from completed union outputs, extracting
the `union-minus-sphere` surface (FBMS structure with the sphere shell removed).

`package_fbms_sim_assets.py`

Copies the final simulation assets into `simulation_package`. Chooses the
original `.veg.obj` when it passed final quality, or `.veg.repaired.obj` when
repair was needed.

## How The Pipeline Works

### 1. Prepare Assets

The runner discovers all `g0_b*` or `g0b*` OBJ files under
`source_asset_directory`, copies them into `working_directory/asset/fbms`, and
generates a unit sphere for each:

```bash
python scripts/generate_bounding_sphere.py \
  --input examples/fbms/fbms_r512/asset/fbms/{case}.obj \
  --output examples/fbms/fbms_r512/asset/fbms/{case}_bounding_sphere.obj \
  --bounds-method unit \
  --method icosphere \
  --subdivisions 8
```

If TPMS baseline shell/sphere files are missing:

```bash
python scripts/generate_tpms_unit_ball.py \
  --cells 1.0 \
  --resolution 512 \
  --sphere-subdivisions 8 \
  --out examples/fbms/fbms_r512/asset/baseline \
  --flat-output
```

### 2. Generate FBMS Raw Union Surfaces

For each `g0_b*`, the runner creates a raw union surface using fixed thickness:

```bash
build/base_no_mkl/bin/generateFBMSUnionSurface openvdb \
  --fbms examples/fbms/fbms_r512/asset/fbms/{case}.obj \
  --sphere examples/fbms/fbms_r512/asset/fbms/{case}_bounding_sphere.obj \
  --fbms-thickness 0.02 \
  --sphere-thickness 0.03 \
  --resolution 512 \
  --padding-ratio 0.08 \
  --output-surface examples/fbms/fbms_r512/raw_mesh_vdb_r512_t002/fbms/{case}_union.obj \
  --enable-truncating \
  --project-fbms-boundary-to-sphere \
  --filter-small-components \
  --min-component-triangles 100 \
  --keep-largest-components 3
```

Raw quality is a hard gate. The raw quality JSON records `enclosed_volume`,
which becomes the volume budget for that case's TPMS baselines.

### 3. Generate Volume-Matched TPMS Baselines

For each `g0_b*` and TPMS shape, the runner generates a baseline raw union
surface in volume-budget mode (negative `--fbms-thickness`):

```bash
build/base_no_mkl/bin/generateFBMSUnionSurface openvdb \
  --fbms examples/fbms/fbms_r512/asset/baseline/{shape}_fbms.obj \
  --sphere examples/fbms/fbms_r512/asset/baseline/{shape}_fbms_bounding_sphere.obj \
  --fbms-thickness -{fbms_raw_enclosed_volume} \
  --sphere-thickness 0.03 \
  --resolution 512 \
  --padding-ratio 0.08 \
  --output-surface examples/fbms/fbms_r512/raw_mesh_vdb_r512_vmatch/baseline/{case}/{shape}_union.obj \
  --enable-truncating \
  --project-fbms-boundary-to-sphere \
  --filter-small-components \
  --min-component-triangles 100 \
  --keep-largest-components 3
```

### 4. Remesh Surfaces

Each raw union surface is remeshed with CGAL isotropic remeshing. The configured
physical target edge length is `0.02`, converted to a relative scale:

```text
edge_length_argument = 0.02 / raw_quality.edge_length.mean
```

For noshell targets, if the CGAL remesh has self-intersections the runner
automatically falls back to geogram at the same vertex count.

### 5. Generate Tet Meshes

The runner writes one `tetmesher.json` per asset and invokes
`build/base_no_mkl/bin/tetMesher --config {tetmesher.json}`. Each output
directory contains `.veg`, `.veg.obj`, `.veg.info.txt`, `.veg.quality.json`,
and `_tetwild.log`.

### 6. Validate And Repair Boundary Surfaces

Final boundary quality requires: 3 edge-connected components, no invalid
triangles, no self-intersections, closed/manifold topology, consistent winding.

If the original `.veg.obj` fails, the runner invokes `repair_tpms_mesh.py`.
Before running repair, it compares mtime of `veg.obj` and any existing
`veg.repaired.obj` — if the repaired file is older than the current `veg.obj`,
it is treated as stale and repair is re-run even without `--overwrite`.

The runner also captures volume mesh info and dumps final components.

### 7. No-Shell Pipeline (Separate Entry Point)

After the main pipeline completes, the no-shell pipeline generates surfaces
with the sphere shell removed (`union-minus-sphere` mode):

```bash
python scripts/run_fbms_noshell_from_union.py \
  --config examples/fbms/fbms_r512/pipeline_r512.json \
  --case g0_b15 \
  --include all \
  --overwrite
```

### 8. Package Simulation Assets

```bash
python scripts/package_fbms_sim_assets.py \
  --config examples/fbms/fbms_r512/pipeline_r512.json \
  --overwrite
```

The package script copies only `.veg` and the final `.veg.obj` (original or
repaired) for each passing asset and writes `package_manifest.json`.

## Common Commands

All scripts use `_resolve_python()` to prefer the project `.venv/bin/python`
over the launching interpreter, so child scripts always have their dependencies.

```bash
# Full pipeline
python scripts/run_fbms_all_asset_pipeline.py \
  --config examples/fbms/fbms_r512/pipeline_r512.json \
  --overwrite --update-readme

# Single case, from raw onward
python scripts/run_fbms_all_asset_pipeline.py \
  --config examples/fbms/fbms_r512/pipeline_r512.json \
  --case g0_b15 --from-stage raw --overwrite

# Validate + repair only
python scripts/run_fbms_all_asset_pipeline.py \
  --config examples/fbms/fbms_r512/pipeline_r512.json \
  --case g0_b15 --from-stage validate --overwrite

# Noshell
python scripts/run_fbms_noshell_from_union.py \
  --config examples/fbms/fbms_r512/pipeline_r512.json \
  --case g0_b15 --include all --overwrite

# Package
python scripts/package_fbms_sim_assets.py \
  --config examples/fbms/fbms_r512/pipeline_r512.json \
  --case g0_b10 --overwrite

# Refresh summary + results table without re-running
python scripts/run_fbms_all_asset_pipeline.py \
  --config examples/fbms/fbms_r512/pipeline_r512.json \
  --from-stage summary --update-readme
```

## Current r512 Results

The detailed result table is in `fbms_r512/README.md`.

- 18 final simulation assets are packaged across `g0_b5`, `g0_b10`, `g0_b15`.
- Each case has 1 FBMS asset + 5 volume-matched TPMS baselines.
- `g0_b5/tpms_schwarz_p` and `g0_b15/tpms_schwarz_p` use repaired final boundary surfaces.
- FBMS and TPMS baselines use the same unit sphere (center at origin, radius 1).
- Raw extraction keeps only the 3 largest components per asset.
- Noshell remesh falls back from CGAL to geogram when self-intersections are detected.

Use `simulation_package/package_manifest.json` for exact provenance of any packaged asset.
