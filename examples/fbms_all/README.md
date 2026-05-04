# FBMS All Asset Pipeline

This directory is a clean, reproducible FBMS simulation asset pipeline.

- Assets: `asset/{baseline,fbms}`
- Raw OpenVDB unions: `raw_mesh_vdb/{baseline,fbms,quality,logs}`
- Geogram remeshes: `remesh_geogram/{baseline,fbms,quality,logs}`
- TetWild simulation meshes: `tetmesh_tetwild/{baseline,fbms}`

The raw and remesh quality gates require `components=3`, manifold, and
watertight surfaces. `invalid_triangles` is treated as a warning for raw mesh
generation, because the current OpenVDB extraction can leave a small number of
near-degenerate triangles while still producing a valid three-component closed
surface. The final remesh and TetWild boundary checks use the full surface
checker and require zero self-intersections.

CGAL isotropic remesh and TetGen were tested in this directory first. That path
did not hold up: CGAL introduced one self-intersecting remesh case, and TetGen
either rejected some inputs or collapsed the three shell components into a
two-component boundary. The current production path is therefore:

```text
asset generation -> OpenVDB union -> Geogram remesh -> TetWild tetmesh
```

## Setup

```bash
uv sync
cmake --build build/base_no_mkl --target generateFBMSUnionSurface meshQualityCheck remeshSurface tetMesher volumetricMeshInfo --parallel
```

## Step 1: Assets

Generate the five baseline TPMS shells and unit bounding spheres:

```bash
uv run python scripts/generate_tpms_unit_ball.py --cells 1.0 --resolution 256 --sphere-subdivisions 4 --out examples/fbms_all/asset/baseline --flat-output
```

Copy the three FBMS source meshes:

```bash
cp examples/fbms_/asset/fbms/g0_b5.obj examples/fbms_all/asset/fbms/g0_b5.obj
cp examples/fbms_/asset/fbms/g0_b10.obj examples/fbms_all/asset/fbms/g0_b10.obj
cp examples/fbms_/asset/fbms/g0_b15.obj examples/fbms_all/asset/fbms/g0_b15.obj
```

Generate the three FBMS bounding spheres:

```bash
uv run python scripts/generate_bounding_sphere.py --input examples/fbms_all/asset/fbms/g0_b5.obj --output examples/fbms_all/asset/fbms/g0_b5_bounding_sphere.obj
uv run python scripts/generate_bounding_sphere.py --input examples/fbms_all/asset/fbms/g0_b10.obj --output examples/fbms_all/asset/fbms/g0_b10_bounding_sphere.obj
uv run python scripts/generate_bounding_sphere.py --input examples/fbms_all/asset/fbms/g0_b15.obj --output examples/fbms_all/asset/fbms/g0_b15_bounding_sphere.obj
```

## Step 2: Raw OpenVDB Meshes

All raw unions use `fbms-thickness=0.03`, `sphere-thickness=0.03`,
`resolution=256`, `padding-ratio=0.08`, boundary projection to the sphere, and
small component filtering.

```bash
build/base_no_mkl/bin/generateFBMSUnionSurface openvdb --fbms examples/fbms_all/asset/baseline/tpms_schwarz_p_fbms.obj --sphere examples/fbms_all/asset/baseline/tpms_schwarz_p_fbms_bounding_sphere.obj --fbms-thickness 0.03 --sphere-thickness 0.03 --resolution 256 --padding-ratio 0.08 --enable-truncating --project-fbms-boundary-to-sphere --filter-small-components --min-component-triangles 100 --output-surface examples/fbms_all/raw_mesh_vdb/baseline/tpms_schwarz_p_union.obj
build/base_no_mkl/bin/generateFBMSUnionSurface openvdb --fbms examples/fbms_all/asset/baseline/tpms_schwarz_d_fbms.obj --sphere examples/fbms_all/asset/baseline/tpms_schwarz_d_fbms_bounding_sphere.obj --fbms-thickness 0.03 --sphere-thickness 0.03 --resolution 256 --padding-ratio 0.08 --enable-truncating --project-fbms-boundary-to-sphere --filter-small-components --min-component-triangles 100 --output-surface examples/fbms_all/raw_mesh_vdb/baseline/tpms_schwarz_d_union.obj
build/base_no_mkl/bin/generateFBMSUnionSurface openvdb --fbms examples/fbms_all/asset/baseline/tpms_gyroid_fbms.obj --sphere examples/fbms_all/asset/baseline/tpms_gyroid_fbms_bounding_sphere.obj --fbms-thickness 0.03 --sphere-thickness 0.03 --resolution 256 --padding-ratio 0.08 --enable-truncating --project-fbms-boundary-to-sphere --filter-small-components --min-component-triangles 100 --output-surface examples/fbms_all/raw_mesh_vdb/baseline/tpms_gyroid_union.obj
build/base_no_mkl/bin/generateFBMSUnionSurface openvdb --fbms examples/fbms_all/asset/baseline/tpms_iwp_fbms.obj --sphere examples/fbms_all/asset/baseline/tpms_iwp_fbms_bounding_sphere.obj --fbms-thickness 0.03 --sphere-thickness 0.03 --resolution 256 --padding-ratio 0.08 --enable-truncating --project-fbms-boundary-to-sphere --filter-small-components --min-component-triangles 100 --output-surface examples/fbms_all/raw_mesh_vdb/baseline/tpms_iwp_union.obj
build/base_no_mkl/bin/generateFBMSUnionSurface openvdb --fbms examples/fbms_all/asset/baseline/tpms_neovius_fbms.obj --sphere examples/fbms_all/asset/baseline/tpms_neovius_fbms_bounding_sphere.obj --fbms-thickness 0.03 --sphere-thickness 0.03 --resolution 256 --padding-ratio 0.08 --enable-truncating --project-fbms-boundary-to-sphere --filter-small-components --min-component-triangles 100 --output-surface examples/fbms_all/raw_mesh_vdb/baseline/tpms_neovius_union.obj
build/base_no_mkl/bin/generateFBMSUnionSurface openvdb --fbms examples/fbms_all/asset/fbms/g0_b5.obj --sphere examples/fbms_all/asset/fbms/g0_b5_bounding_sphere.obj --fbms-thickness 0.03 --sphere-thickness 0.03 --resolution 256 --padding-ratio 0.08 --enable-truncating --project-fbms-boundary-to-sphere --filter-small-components --min-component-triangles 100 --output-surface examples/fbms_all/raw_mesh_vdb/fbms/g0_b5_union.obj
build/base_no_mkl/bin/generateFBMSUnionSurface openvdb --fbms examples/fbms_all/asset/fbms/g0_b10.obj --sphere examples/fbms_all/asset/fbms/g0_b10_bounding_sphere.obj --fbms-thickness 0.03 --sphere-thickness 0.03 --resolution 256 --padding-ratio 0.08 --enable-truncating --project-fbms-boundary-to-sphere --filter-small-components --min-component-triangles 100 --output-surface examples/fbms_all/raw_mesh_vdb/fbms/g0_b10_union.obj
build/base_no_mkl/bin/generateFBMSUnionSurface openvdb --fbms examples/fbms_all/asset/fbms/g0_b15.obj --sphere examples/fbms_all/asset/fbms/g0_b15_bounding_sphere.obj --fbms-thickness 0.03 --sphere-thickness 0.03 --resolution 256 --padding-ratio 0.08 --enable-truncating --project-fbms-boundary-to-sphere --filter-small-components --min-component-triangles 100 --output-surface examples/fbms_all/raw_mesh_vdb/fbms/g0_b15_union.obj
```

## Raw Quality

```bash
build/base_no_mkl/bin/meshQualityCheck surface --check-level raw --expected-components 3 --invalid-triangles-policy warn --input examples/fbms_all/raw_mesh_vdb/baseline/tpms_schwarz_p_union.obj --json examples/fbms_all/raw_mesh_vdb/quality/baseline/tpms_schwarz_p_union.quality.json
build/base_no_mkl/bin/meshQualityCheck surface --check-level raw --expected-components 3 --invalid-triangles-policy warn --input examples/fbms_all/raw_mesh_vdb/baseline/tpms_schwarz_d_union.obj --json examples/fbms_all/raw_mesh_vdb/quality/baseline/tpms_schwarz_d_union.quality.json
build/base_no_mkl/bin/meshQualityCheck surface --check-level raw --expected-components 3 --invalid-triangles-policy warn --input examples/fbms_all/raw_mesh_vdb/baseline/tpms_gyroid_union.obj --json examples/fbms_all/raw_mesh_vdb/quality/baseline/tpms_gyroid_union.quality.json
build/base_no_mkl/bin/meshQualityCheck surface --check-level raw --expected-components 3 --invalid-triangles-policy warn --input examples/fbms_all/raw_mesh_vdb/baseline/tpms_iwp_union.obj --json examples/fbms_all/raw_mesh_vdb/quality/baseline/tpms_iwp_union.quality.json
build/base_no_mkl/bin/meshQualityCheck surface --check-level raw --expected-components 3 --invalid-triangles-policy warn --input examples/fbms_all/raw_mesh_vdb/baseline/tpms_neovius_union.obj --json examples/fbms_all/raw_mesh_vdb/quality/baseline/tpms_neovius_union.quality.json
build/base_no_mkl/bin/meshQualityCheck surface --check-level raw --expected-components 3 --invalid-triangles-policy warn --input examples/fbms_all/raw_mesh_vdb/fbms/g0_b5_union.obj --json examples/fbms_all/raw_mesh_vdb/quality/fbms/g0_b5_union.quality.json
build/base_no_mkl/bin/meshQualityCheck surface --check-level raw --expected-components 3 --invalid-triangles-policy warn --input examples/fbms_all/raw_mesh_vdb/fbms/g0_b10_union.obj --json examples/fbms_all/raw_mesh_vdb/quality/fbms/g0_b10_union.quality.json
build/base_no_mkl/bin/meshQualityCheck surface --check-level raw --expected-components 3 --invalid-triangles-policy warn --input examples/fbms_all/raw_mesh_vdb/fbms/g0_b15_union.obj --json examples/fbms_all/raw_mesh_vdb/quality/fbms/g0_b15_union.quality.json
```

## Step 3: Geogram Remesh

The target vertex counts use one density rule for all cases:

```text
target-num-vertices = round_to_nearest_1000(563 * raw_surface_area)
```

Geogram is run with one thread for reproducibility and
`--orient-nested-components` so the outer sphere remains outward and the two
inner shell/cavity components are oriented inward before TetWild.

```bash
build/base_no_mkl/bin/remeshSurface geogram --threads 1 --orient-nested-components --input-mesh examples/fbms_all/raw_mesh_vdb/baseline/tpms_schwarz_p_union.obj --output-mesh examples/fbms_all/remesh_geogram/baseline/tpms_schwarz_p_remesh.obj --target-num-vertices 20000
build/base_no_mkl/bin/remeshSurface geogram --threads 1 --orient-nested-components --input-mesh examples/fbms_all/raw_mesh_vdb/baseline/tpms_schwarz_d_union.obj --output-mesh examples/fbms_all/remesh_geogram/baseline/tpms_schwarz_d_remesh.obj --target-num-vertices 22000
build/base_no_mkl/bin/remeshSurface geogram --threads 1 --orient-nested-components --input-mesh examples/fbms_all/raw_mesh_vdb/baseline/tpms_gyroid_union.obj --output-mesh examples/fbms_all/remesh_geogram/baseline/tpms_gyroid_remesh.obj --target-num-vertices 21000
build/base_no_mkl/bin/remeshSurface geogram --threads 1 --orient-nested-components --input-mesh examples/fbms_all/raw_mesh_vdb/baseline/tpms_iwp_union.obj --output-mesh examples/fbms_all/remesh_geogram/baseline/tpms_iwp_remesh.obj --target-num-vertices 23000
build/base_no_mkl/bin/remeshSurface geogram --threads 1 --orient-nested-components --input-mesh examples/fbms_all/raw_mesh_vdb/baseline/tpms_neovius_union.obj --output-mesh examples/fbms_all/remesh_geogram/baseline/tpms_neovius_remesh.obj --target-num-vertices 23000
build/base_no_mkl/bin/remeshSurface geogram --threads 1 --orient-nested-components --input-mesh examples/fbms_all/raw_mesh_vdb/fbms/g0_b5_union.obj --output-mesh examples/fbms_all/remesh_geogram/fbms/g0_b5_remesh.obj --target-num-vertices 22000
build/base_no_mkl/bin/remeshSurface geogram --threads 1 --orient-nested-components --input-mesh examples/fbms_all/raw_mesh_vdb/fbms/g0_b10_union.obj --output-mesh examples/fbms_all/remesh_geogram/fbms/g0_b10_remesh.obj --target-num-vertices 24000
build/base_no_mkl/bin/remeshSurface geogram --threads 1 --orient-nested-components --input-mesh examples/fbms_all/raw_mesh_vdb/fbms/g0_b15_union.obj --output-mesh examples/fbms_all/remesh_geogram/fbms/g0_b15_remesh.obj --target-num-vertices 24000
```

## Remesh Quality

```bash
build/base_no_mkl/bin/meshQualityCheck surface --check-level full --expected-components 3 --invalid-triangles-policy warn --input examples/fbms_all/remesh_geogram/baseline/tpms_schwarz_p_remesh.obj --json examples/fbms_all/remesh_geogram/quality/baseline/tpms_schwarz_p_remesh.quality.json
build/base_no_mkl/bin/meshQualityCheck surface --check-level full --expected-components 3 --invalid-triangles-policy warn --input examples/fbms_all/remesh_geogram/baseline/tpms_schwarz_d_remesh.obj --json examples/fbms_all/remesh_geogram/quality/baseline/tpms_schwarz_d_remesh.quality.json
build/base_no_mkl/bin/meshQualityCheck surface --check-level full --expected-components 3 --invalid-triangles-policy warn --input examples/fbms_all/remesh_geogram/baseline/tpms_gyroid_remesh.obj --json examples/fbms_all/remesh_geogram/quality/baseline/tpms_gyroid_remesh.quality.json
build/base_no_mkl/bin/meshQualityCheck surface --check-level full --expected-components 3 --invalid-triangles-policy warn --input examples/fbms_all/remesh_geogram/baseline/tpms_iwp_remesh.obj --json examples/fbms_all/remesh_geogram/quality/baseline/tpms_iwp_remesh.quality.json
build/base_no_mkl/bin/meshQualityCheck surface --check-level full --expected-components 3 --invalid-triangles-policy warn --input examples/fbms_all/remesh_geogram/baseline/tpms_neovius_remesh.obj --json examples/fbms_all/remesh_geogram/quality/baseline/tpms_neovius_remesh.quality.json
build/base_no_mkl/bin/meshQualityCheck surface --check-level full --expected-components 3 --invalid-triangles-policy warn --input examples/fbms_all/remesh_geogram/fbms/g0_b5_remesh.obj --json examples/fbms_all/remesh_geogram/quality/fbms/g0_b5_remesh.quality.json
build/base_no_mkl/bin/meshQualityCheck surface --check-level full --expected-components 3 --invalid-triangles-policy warn --input examples/fbms_all/remesh_geogram/fbms/g0_b10_remesh.obj --json examples/fbms_all/remesh_geogram/quality/fbms/g0_b10_remesh.quality.json
build/base_no_mkl/bin/meshQualityCheck surface --check-level full --expected-components 3 --invalid-triangles-policy warn --input examples/fbms_all/remesh_geogram/fbms/g0_b15_remesh.obj --json examples/fbms_all/remesh_geogram/quality/fbms/g0_b15_remesh.quality.json
```

## Step 4: TetWild

Each case has a local `tetmesher.json` using:

```json
{
  "backend": "tetwild",
  "material": {
    "density": 1000,
    "young_modulus": 10000000,
    "poisson_ratio": 0.45
  },
  "tetwild": {
    "lr": 0.013,
    "epsr": 0.001,
    "stop_energy": 10
  }
}
```

The `material` block is optional in `tetMesher`; if omitted, `tetMesher` writes
the Vega default material. These assets intentionally override only Young's
modulus to `10000000` while keeping density `1000` and Poisson ratio `0.45`.

```bash
build/base_no_mkl/bin/tetMesher --config examples/fbms_all/tetmesh_tetwild/baseline/tpms_schwarz_p/tetmesher.json
build/base_no_mkl/bin/tetMesher --config examples/fbms_all/tetmesh_tetwild/baseline/tpms_schwarz_d/tetmesher.json
build/base_no_mkl/bin/tetMesher --config examples/fbms_all/tetmesh_tetwild/baseline/tpms_gyroid/tetmesher.json
build/base_no_mkl/bin/tetMesher --config examples/fbms_all/tetmesh_tetwild/baseline/tpms_iwp/tetmesher.json
build/base_no_mkl/bin/tetMesher --config examples/fbms_all/tetmesh_tetwild/baseline/tpms_neovius/tetmesher.json
build/base_no_mkl/bin/tetMesher --config examples/fbms_all/tetmesh_tetwild/fbms/g0_b5/tetmesher.json
build/base_no_mkl/bin/tetMesher --config examples/fbms_all/tetmesh_tetwild/fbms/g0_b10/tetmesher.json
build/base_no_mkl/bin/tetMesher --config examples/fbms_all/tetmesh_tetwild/fbms/g0_b15/tetmesher.json
```

## Tet Boundary Quality

```bash
build/base_no_mkl/bin/meshQualityCheck surface --check-level full --expected-components 3 --invalid-triangles-policy warn --input examples/fbms_all/tetmesh_tetwild/baseline/tpms_schwarz_p/tpms_schwarz_p.veg.obj --json examples/fbms_all/tetmesh_tetwild/baseline/tpms_schwarz_p/tpms_schwarz_p.veg.quality.json
build/base_no_mkl/bin/meshQualityCheck surface --check-level full --expected-components 3 --invalid-triangles-policy warn --input examples/fbms_all/tetmesh_tetwild/baseline/tpms_schwarz_d/tpms_schwarz_d.veg.obj --json examples/fbms_all/tetmesh_tetwild/baseline/tpms_schwarz_d/tpms_schwarz_d.veg.quality.json
build/base_no_mkl/bin/meshQualityCheck surface --check-level full --expected-components 3 --invalid-triangles-policy warn --input examples/fbms_all/tetmesh_tetwild/baseline/tpms_gyroid/tpms_gyroid.veg.obj --json examples/fbms_all/tetmesh_tetwild/baseline/tpms_gyroid/tpms_gyroid.veg.quality.json
build/base_no_mkl/bin/meshQualityCheck surface --check-level full --expected-components 3 --invalid-triangles-policy warn --input examples/fbms_all/tetmesh_tetwild/baseline/tpms_iwp/tpms_iwp.veg.obj --json examples/fbms_all/tetmesh_tetwild/baseline/tpms_iwp/tpms_iwp.veg.quality.json
build/base_no_mkl/bin/meshQualityCheck surface --check-level full --expected-components 3 --invalid-triangles-policy warn --input examples/fbms_all/tetmesh_tetwild/baseline/tpms_neovius/tpms_neovius.veg.obj --json examples/fbms_all/tetmesh_tetwild/baseline/tpms_neovius/tpms_neovius.veg.quality.json
build/base_no_mkl/bin/meshQualityCheck surface --check-level full --expected-components 3 --invalid-triangles-policy warn --input examples/fbms_all/tetmesh_tetwild/fbms/g0_b5/g0_b5.veg.obj --json examples/fbms_all/tetmesh_tetwild/fbms/g0_b5/g0_b5.veg.quality.json
build/base_no_mkl/bin/meshQualityCheck surface --check-level full --expected-components 3 --invalid-triangles-policy warn --input examples/fbms_all/tetmesh_tetwild/fbms/g0_b10/g0_b10.veg.obj --json examples/fbms_all/tetmesh_tetwild/fbms/g0_b10/g0_b10.veg.quality.json
build/base_no_mkl/bin/meshQualityCheck surface --check-level full --expected-components 3 --invalid-triangles-policy warn --input examples/fbms_all/tetmesh_tetwild/fbms/g0_b15/g0_b15.veg.obj --json examples/fbms_all/tetmesh_tetwild/fbms/g0_b15/g0_b15.veg.quality.json
```

## Volumetric Mesh Info

```bash
build/base_no_mkl/bin/volumetricMeshInfo examples/fbms_all/tetmesh_tetwild/baseline/tpms_schwarz_p/tpms_schwarz_p.veg
build/base_no_mkl/bin/volumetricMeshInfo examples/fbms_all/tetmesh_tetwild/baseline/tpms_schwarz_d/tpms_schwarz_d.veg
build/base_no_mkl/bin/volumetricMeshInfo examples/fbms_all/tetmesh_tetwild/baseline/tpms_gyroid/tpms_gyroid.veg
build/base_no_mkl/bin/volumetricMeshInfo examples/fbms_all/tetmesh_tetwild/baseline/tpms_iwp/tpms_iwp.veg
build/base_no_mkl/bin/volumetricMeshInfo examples/fbms_all/tetmesh_tetwild/baseline/tpms_neovius/tpms_neovius.veg
build/base_no_mkl/bin/volumetricMeshInfo examples/fbms_all/tetmesh_tetwild/fbms/g0_b5/g0_b5.veg
build/base_no_mkl/bin/volumetricMeshInfo examples/fbms_all/tetmesh_tetwild/fbms/g0_b10/g0_b10.veg
build/base_no_mkl/bin/volumetricMeshInfo examples/fbms_all/tetmesh_tetwild/fbms/g0_b15/g0_b15.veg
```

## Tet Boundary Component Dumps

Dump each edge-connected component from a TetWild `.veg.obj` boundary. The
script only splits one OBJ at a time and writes component metadata; it does not
run quality checks.

```bash
uv run python scripts/dump_obj_components.py --input examples/fbms_all/tetmesh_tetwild/baseline/tpms_schwarz_p/tpms_schwarz_p.veg.obj --output-dir examples/fbms_all/tetmesh_tetwild_components/baseline/tpms_schwarz_p --overwrite
uv run python scripts/dump_obj_components.py --input examples/fbms_all/tetmesh_tetwild/baseline/tpms_schwarz_d/tpms_schwarz_d.veg.obj --output-dir examples/fbms_all/tetmesh_tetwild_components/baseline/tpms_schwarz_d --overwrite
uv run python scripts/dump_obj_components.py --input examples/fbms_all/tetmesh_tetwild/baseline/tpms_gyroid/tpms_gyroid.veg.obj --output-dir examples/fbms_all/tetmesh_tetwild_components/baseline/tpms_gyroid --overwrite
uv run python scripts/dump_obj_components.py --input examples/fbms_all/tetmesh_tetwild/baseline/tpms_iwp/tpms_iwp.veg.obj --output-dir examples/fbms_all/tetmesh_tetwild_components/baseline/tpms_iwp --overwrite
uv run python scripts/dump_obj_components.py --input examples/fbms_all/tetmesh_tetwild/baseline/tpms_neovius/tpms_neovius.veg.obj --output-dir examples/fbms_all/tetmesh_tetwild_components/baseline/tpms_neovius --overwrite
uv run python scripts/dump_obj_components.py --input examples/fbms_all/tetmesh_tetwild/fbms/g0_b5/g0_b5.veg.obj --output-dir examples/fbms_all/tetmesh_tetwild_components/fbms/g0_b5 --overwrite
uv run python scripts/dump_obj_components.py --input examples/fbms_all/tetmesh_tetwild/fbms/g0_b10/g0_b10.veg.obj --output-dir examples/fbms_all/tetmesh_tetwild_components/fbms/g0_b10 --overwrite
uv run python scripts/dump_obj_components.py --input examples/fbms_all/tetmesh_tetwild/fbms/g0_b15/g0_b15.veg.obj --output-dir examples/fbms_all/tetmesh_tetwild_components/fbms/g0_b15 --overwrite
```

The output goes to `tetmesh_tetwild_components/{baseline,fbms}/{case}`. Each
case has three OBJ files sorted by descending triangle count:

```text
{case}_component_0.obj
{case}_component_1.obj
{case}_component_2.obj
```

The script also writes `components.json` per case.

## Current Results

These results were regenerated with the commands above using
`ENU, 1000, 10000000, 0.45`.

| case | raw comps / invalid | remesh vertices / tris | remesh quality | tet vertices | tet elements | tet volume | tet boundary quality |
| --- | ---: | ---: | --- | ---: | ---: | ---: | --- |
| tpms_schwarz_p | 3 / 0 | 22235 / 44478 | pass | 22842 | 72010 | 0.567128 | pass |
| tpms_schwarz_d | 3 / 0 | 22829 / 45658 | pass | 23175 | 71298 | 0.592264 | pass |
| tpms_gyroid | 3 / 36 | 22079 / 44162 | pass | 22556 | 69371 | 0.587040 | pass |
| tpms_iwp | 3 / 0 | 23598 / 47212 | pass | 23811 | 73200 | 0.623365 | pass |
| tpms_neovius | 3 / 0 | 24706 / 49444 | pass | 25968 | 81405 | 0.637320 | pass |
| g0_b5 | 3 / 6 | 22709 / 45422 | pass | 22633 | 68996 | 0.621550 | pass |
| g0_b10 | 3 / 2 | 24994 / 50012 | pass | 24680 | 75470 | 0.679406 | pass |
| g0_b15 | 3 / 6 | 26040 / 52124 | pass | 25408 | 77774 | 0.710981 | pass |

The final TetWild boundary quality reports all have:

```text
components_by_edge=3
boundary_or_exterior_edges=0
is_manifold=true
invalid_triangles=0
self_intersections=0
passed=true
```
