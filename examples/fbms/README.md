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

Creates a bounding sphere OBJ for an input shell mesh. The FBMS pipeline uses it
when a `{case}_bounding_sphere.obj` is not already present. The current config
uses an AABB-based sphere and an icosphere mesh.

`generate_tpms_unit_ball.py`

Generates the five TPMS baseline shell meshes and their bounding spheres. With
`--flat-output`, it writes files directly as
`{shape}_fbms.obj` and `{shape}_fbms_bounding_sphere.obj`.

`generateFBMSUnionSurface openvdb`

Builds a surface union between a shell mesh and its bounding sphere shell using
OpenVDB. Positive `--fbms-thickness` means fixed shell thickness. Negative
`--fbms-thickness -V` means volume-budget mode: the tool searches for a shell
thickness whose raw union volume is close to `V`.

`meshQualityCheck surface`

Checks triangle surface quality and writes a JSON report. Raw checks are used
as hard gates for generated raw union meshes. Full checks are used for remesh
and final TetWild boundary surfaces. The JSON records component count,
invalid triangles, manifold/closed/winding status, self-intersections, edge
length statistics, and enclosed volume when meaningful.

`remeshSurface cgal_iso`

Runs CGAL isotropic surface remeshing. In this pipeline, it regularizes raw
OpenVDB triangle meshes before TetWild. The tool's `--edge-length` argument is
relative to the raw mesh mean edge length, so the runner computes a scale from
the desired physical target edge length.

`tetMesher`

Reads a `tetmesher.json` config, runs the TetWild backend, writes a `.veg`
volume mesh, and exports a `.veg.obj` boundary surface. The `.veg` is the
simulation mesh; the `.veg.obj` is the surface used for final quality checks
and packaging.

`volumetricMeshInfo`

Reads a `.veg` mesh and prints vertex count, element count, and total volume.
The runner saves this as `{name}.veg.info.txt` and uses it in the summary.

`repair_tpms_mesh.py`

Repairs a failed TetWild boundary OBJ after TetWild has completed. It uses
PyMeshFix per connected component, preserving components as much as possible.
It never overwrites the original `.veg` or `.veg.obj`; repaired surfaces are
written as `{name}.veg.repaired.obj`.

`dump_obj_components.py`

Splits the final boundary surface into connected components and writes a
`components.json` manifest plus per-component OBJ files. The pipeline requires
three final components.

`package_fbms_sim_assets.py`

Copies the final simulation assets into `simulation_package`. It chooses the
original `.veg.obj` when it passed final quality, or `.veg.repaired.obj` when
repair was needed, and records the source in `package_manifest.json`.

## How The Pipeline Works

### 1. Prepare Assets

The runner discovers all `g0_b*` or `g0b*` OBJ files under
`source_asset_directory`, copies them into `working_directory/asset/fbms`, and
ensures each has a bounding sphere.

If TPMS baseline shell/sphere files are missing, the runner can generate them
with:

- `scripts/generate_bounding_sphere.py`
- `scripts/generate_tpms_unit_ball.py --flat-output`

The commands generated by the runner look like:

```bash
uv run python scripts/generate_bounding_sphere.py \
  --input examples/fbms/fbms_r512/asset/fbms/{case}.obj \
  --output examples/fbms/fbms_r512/asset/fbms/{case}_bounding_sphere.obj \
  --bounds-method aabb \
  --method icosphere \
  --subdivisions 8 \
  --padding 1e-09
```

```bash
uv run python scripts/generate_tpms_unit_ball.py \
  --cells 1.0 \
  --resolution 512 \
  --sphere-subdivisions 8 \
  --out examples/fbms/fbms_r512/asset/baseline \
  --flat-output
```

For the current checked-in r512 run, the prepared assets already live under
`examples/fbms/fbms_r512/asset`.

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
  --min-component-triangles 100
```

Immediately after raw generation, the runner checks quality:

```bash
build/base_no_mkl/bin/meshQualityCheck surface \
  --check-level raw \
  --expected-components 3 \
  --invalid-triangles-policy warn \
  --input examples/fbms/fbms_r512/raw_mesh_vdb_r512_t002/fbms/{case}_union.obj \
  --json examples/fbms/fbms_r512/raw_mesh_vdb_r512_t002/fbms/quality/{case}_union.quality.json
```

Raw quality is a hard gate. The raw quality JSON also records
`enclosed_volume`, which becomes the volume budget for that case's TPMS
baselines.

### 3. Generate Volume-Matched TPMS Baselines

For each `g0_b*` and each TPMS shape, the runner generates a baseline raw union
surface with negative `--fbms-thickness`:

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
  --min-component-triangles 100
```

This tells `generateFBMSUnionSurface openvdb` to search for a shell thickness
that matches the target raw union volume.

### 4. Remesh Surfaces

Each raw union surface is remeshed with CGAL isotropic remeshing. The configured
physical target edge length is `0.02`, but `remeshSurface cgal_iso` takes a
relative edge-length scale, so the runner computes:

```text
edge_length_argument = 0.02 / raw_quality.edge_length.mean
```

The remesh command template is:

```bash
build/base_no_mkl/bin/remeshSurface cgal_iso \
  --input-mesh {raw_union_obj} \
  --output-mesh {remesh_obj} \
  --edge-length {edge_length_argument} \
  --sharp-edge-angle 180 \
  --iterations 10
```

The remesh quality command is:

```bash
build/base_no_mkl/bin/meshQualityCheck surface \
  --check-level full \
  --expected-components 3 \
  --invalid-triangles-policy warn \
  --input {remesh_obj} \
  --json {remesh_quality_json} \
  --self-intersection-backend exact-count \
  --self-intersection-triangle-limit 900000
```

Remesh quality is recorded, but self-intersections at this stage are treated as
warnings so that TetWild can still attempt recovery.

### 5. Generate Tet Meshes

The runner writes one `tetmesher.json` per asset and invokes:

```bash
build/base_no_mkl/bin/tetMesher --config {tet_output_dir}/tetmesher.json
```

The generated `tetmesher.json` has this shape:

```json
{
  "version": 1,
  "backend": "tetwild",
  "input_mesh": "../../../remesh_cgal_iso_r512_vmatch_e002/baseline/g0_b5/tpms_schwarz_p_remesh.obj",
  "output_mesh": "tpms_schwarz_p.veg",
  "output_surface": "tpms_schwarz_p.veg.obj",
  "print_stats": true,
  "tetwild": {
    "lr": 0.008,
    "epsr": 0.001
  },
  "material": {
    "density": 1000,
    "young_modulus": 10000000,
    "poisson_ratio": 0.45
  }
}
```

Each TetWild output directory contains:

- `{name}.veg`
- `{name}.veg.obj`
- `{name}.veg.info.txt`
- `{name}.veg.quality.json`
- `{name}_tetwild.log`
- `tetmesher.json`

The `.veg` file is the simulation volume mesh. The `.veg.obj` file is the
surface exported from the TetWild output.

### 6. Validate And Repair Boundary Surfaces

Final boundary quality is checked with `meshQualityCheck surface --check-level
full`. The final packaged surface must have:

- 3 edge-connected components
- no invalid triangles
- no self-intersections
- closed/manifold topology
- consistent winding

The final boundary check is:

```bash
build/base_no_mkl/bin/meshQualityCheck surface \
  --check-level full \
  --expected-components 3 \
  --invalid-triangles-policy warn \
  --input {tet_output_dir}/{name}.veg.obj \
  --json {tet_output_dir}/{name}.veg.quality.json \
  --self-intersection-backend exact-count \
  --self-intersection-triangle-limit 900000
```

If the original `.veg.obj` fails, the runner uses:

```bash
uv run python scripts/repair_tpms_mesh.py \
  {tet_output_dir}/{name}.veg.obj \
  --out {tet_output_dir}/{name}.veg.repaired.obj \
  --report {tet_output_dir}/{name}.veg.repair_report.json \
  --orient-positive
```

Then it checks the repaired surface with the same full quality gate:

```bash
build/base_no_mkl/bin/meshQualityCheck surface \
  --check-level full \
  --expected-components 3 \
  --invalid-triangles-policy warn \
  --input {tet_output_dir}/{name}.veg.repaired.obj \
  --json {tet_output_dir}/{name}.veg.repaired.quality.json \
  --self-intersection-backend exact-count \
  --self-intersection-triangle-limit 900000
```

This writes:

- `{name}.veg.repaired.obj`
- `{name}.veg.repaired.quality.json`
- `{name}.veg.repair_report.json`
- `{name}.veg.repaired.bad_edges.csv`

The original `.veg` and `.veg.obj` are never overwritten. If the repaired
surface passes final quality, it becomes the final surface for packaging. This
is why the r512 result table may show `pass(repaired)`.

The runner also captures volume mesh info and dumps final components:

```bash
build/base_no_mkl/bin/volumetricMeshInfo \
  {tet_output_dir}/{name}.veg \
  > {tet_output_dir}/{name}.veg.info.txt
```

```bash
uv run python scripts/dump_obj_components.py \
  --input {final_boundary_obj} \
  --output-dir {component_output_dir} \
  --overwrite
```

### 7. Package Simulation Assets

After validation, package the final simulation assets with:

```bash
uv run python scripts/package_fbms_sim_assets.py \
  --config examples/fbms/fbms_r512/pipeline_r512.json \
  --overwrite
```

The package script copies only the final simulation files:

- `.veg`: always from TetWild
- `.veg.obj`: original boundary surface if it passed, repaired boundary surface
  if repair was needed

It writes a `package_manifest.json` with source paths, package paths, repair
status, and repair volume drift.

## Common Commands

Dry-run the full pipeline:

```bash
uv run python scripts/run_fbms_all_asset_pipeline.py \
  --config examples/fbms/fbms_r512/pipeline_r512.json \
  --dry-run
```

Run or rerun the full pipeline:

```bash
uv run python scripts/run_fbms_all_asset_pipeline.py \
  --config examples/fbms/fbms_r512/pipeline_r512.json \
  --overwrite \
  --update-readme
```

Rerun only validation and README/summary refresh:

```bash
uv run python scripts/run_fbms_all_asset_pipeline.py \
  --config examples/fbms/fbms_r512/pipeline_r512.json \
  --from-stage validate \
  --overwrite \
  --update-readme
```

Rerun a single case and shape:

```bash
uv run python scripts/run_fbms_all_asset_pipeline.py \
  --config examples/fbms/fbms_r512/pipeline_r512.json \
  --case g0_b15 \
  --shape tpms_schwarz_p \
  --from-stage validate \
  --overwrite \
  --update-readme
```

Refresh summary only:

```bash
uv run python scripts/run_fbms_all_asset_pipeline.py \
  --config examples/fbms/fbms_r512/pipeline_r512.json \
  --from-stage summary \
  --update-readme
```

Package final simulation assets:

```bash
uv run python scripts/package_fbms_sim_assets.py \
  --config examples/fbms/fbms_r512/pipeline_r512.json \
  --overwrite
```

Package only one case:

```bash
uv run python scripts/package_fbms_sim_assets.py \
  --config examples/fbms/fbms_r512/pipeline_r512.json \
  --case g0_b10 \
  --overwrite
```

## Current r512 Results

The detailed result table is in:

```bash
examples/fbms/fbms_r512/README.md
```

Current summary:

- 18 final simulation assets are packaged.
- Each `g0_b*` has 1 FBMS asset and 5 volume-matched TPMS baseline assets.
- `g0_b5/tpms_schwarz_p` and `g0_b15/tpms_schwarz_p` use repaired final
  boundary surfaces.
- The packaged assets live in `examples/fbms/fbms_r512/simulation_package`.

Use `simulation_package/package_manifest.json` when you need exact provenance
for any packaged `.veg` or `.veg.obj`.
