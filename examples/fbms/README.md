# FBMS Asset Generation

This directory contains FBMS surface inputs and generated simulation-ready geometry assets.

The current b3/b8 shell assets use SDF thickening plus marching cubes at resolutions 64, 128, and 256, followed by CGAL isotropic remeshing.

## Batch Script

Use `shell_assets.json` to define raw SDF-union and remeshing jobs. The script command line is only for execution control: selecting jobs, dry-run mode, and overwrite/skip behavior.

Preview the lightweight thick-shell job without generating assets:

```bash
examples/fbms/generate_shell_assets.py \
  --job r64_thick \
  --dry-run
```

Generate the lightweight thick-shell job:

```bash
examples/fbms/generate_shell_assets.py \
  --job r64_thick \
  --skip-existing
```

Generate multiple jobs in one run:

```bash
examples/fbms/generate_shell_assets.py \
  --job r64_thick \
  --job r128_default \
  --skip-existing
```

Run every job from the JSON config explicitly:

```bash
examples/fbms/generate_shell_assets.py \
  --all-jobs \
  --skip-existing
```

Each job writes to its own folder under `examples/fbms/generated/<job>/<case>/`:

```text
generated/r64_thick/g0_b3/union_shell_raw.obj
generated/r64_thick/g0_b3/union_shell_remesh.obj
generated/r64_thick/g0_b3/stats.json
generated/r64_thick/g0_b8/union_shell_raw.obj
generated/r64_thick/g0_b8/union_shell_remesh.obj
generated/r64_thick/g0_b8/stats.json
```

Edit `shell_assets.json` to change `resolution`, `fbms_thickness`, `sphere_thickness`, `padding_ratio`, `edge_length`, `sharp_edge_angle`, or per-job `output_dir`. Pass `--overwrite` to replace existing outputs, or `--skip-existing` to leave them untouched.

The sections below document the existing flat assets and the equivalent manual commands. New experiment outputs should prefer the JSON job layout under `generated/`.

## Inputs

```text
g0_b3/g0_b3_fbms.obj
g0_b3/g0_b3_fbms_bounding_sphere.obj
g0_b8/g0_b8_fbms.obj
g0_b8/g0_b8_fbms_bounding_sphere.obj
```

## SDF Union Raw Surfaces

Generate the b3 raw shell:

```bash
build/base_no_mkl/bin/generateFBMSUnionSurface \
  --fbms examples/fbms/g0_b3/g0_b3_fbms.obj \
  --sphere examples/fbms/g0_b3/g0_b3_fbms_bounding_sphere.obj \
  --fbms-thickness 0.02 \
  --sphere-thickness 0.02 \
  --resolution 256 \
  --padding-ratio 0.08 \
  --output-surface examples/fbms/g0_b3/g0_b3_union_shell_raw.obj
```

Generate the b8 raw shell:

```bash
build/base_no_mkl/bin/generateFBMSUnionSurface \
  --fbms examples/fbms/g0_b8/g0_b8_fbms.obj \
  --sphere examples/fbms/g0_b8/g0_b8_fbms_bounding_sphere.obj \
  --fbms-thickness 0.02 \
  --sphere-thickness 0.02 \
  --resolution 256 \
  --padding-ratio 0.08 \
  --output-surface examples/fbms/g0_b8/g0_b8_union_shell_raw.obj
```

## Remeshing

Remesh the b3 raw shell:

```bash
build/base_no_mkl/bin/remeshSurface cgal_iso \
  --input-mesh examples/fbms/g0_b3/g0_b3_union_shell_raw.obj \
  --output-mesh examples/fbms/g0_b3/g0_b3_union_shell_remesh.obj \
  --edge-length 0.75 \
  --sharp-edge-angle 180
```

Remesh the b8 raw shell:

```bash
build/base_no_mkl/bin/remeshSurface cgal_iso \
  --input-mesh examples/fbms/g0_b8/g0_b8_union_shell_raw.obj \
  --output-mesh examples/fbms/g0_b8/g0_b8_union_shell_remesh.obj \
  --edge-length 0.75 \
  --sharp-edge-angle 180
```

`remeshSurface cgal_iso --edge-length` is a scale relative to the input average edge length. The value `0.75` slightly refines the marching-cubes mesh.

## Additional Resolutions

Resolution-specific assets use a filename suffix:

```text
g0_b*/g0_b*_union_shell_r128_raw.obj
g0_b*/g0_b*_union_shell_r128_remesh.obj
g0_b*/g0_b*_union_shell_r64_raw.obj
g0_b*/g0_b*_union_shell_r64_remesh.obj
```

Generate b3 at resolution 128:

```bash
build/base_no_mkl/bin/generateFBMSUnionSurface \
  --fbms examples/fbms/g0_b3/g0_b3_fbms.obj \
  --sphere examples/fbms/g0_b3/g0_b3_fbms_bounding_sphere.obj \
  --fbms-thickness 0.02 \
  --sphere-thickness 0.02 \
  --resolution 128 \
  --padding-ratio 0.08 \
  --output-surface examples/fbms/g0_b3/g0_b3_union_shell_r128_raw.obj

build/base_no_mkl/bin/remeshSurface cgal_iso \
  --input-mesh examples/fbms/g0_b3/g0_b3_union_shell_r128_raw.obj \
  --output-mesh examples/fbms/g0_b3/g0_b3_union_shell_r128_remesh.obj \
  --edge-length 0.75 \
  --sharp-edge-angle 180
```

Generate b8 at resolution 128:

```bash
build/base_no_mkl/bin/generateFBMSUnionSurface \
  --fbms examples/fbms/g0_b8/g0_b8_fbms.obj \
  --sphere examples/fbms/g0_b8/g0_b8_fbms_bounding_sphere.obj \
  --fbms-thickness 0.02 \
  --sphere-thickness 0.02 \
  --resolution 128 \
  --padding-ratio 0.08 \
  --output-surface examples/fbms/g0_b8/g0_b8_union_shell_r128_raw.obj

build/base_no_mkl/bin/remeshSurface cgal_iso \
  --input-mesh examples/fbms/g0_b8/g0_b8_union_shell_r128_raw.obj \
  --output-mesh examples/fbms/g0_b8/g0_b8_union_shell_r128_remesh.obj \
  --edge-length 0.75 \
  --sharp-edge-angle 180
```

Generate b3 at resolution 64:

```bash
build/base_no_mkl/bin/generateFBMSUnionSurface \
  --fbms examples/fbms/g0_b3/g0_b3_fbms.obj \
  --sphere examples/fbms/g0_b3/g0_b3_fbms_bounding_sphere.obj \
  --fbms-thickness 0.02 \
  --sphere-thickness 0.02 \
  --resolution 64 \
  --padding-ratio 0.08 \
  --output-surface examples/fbms/g0_b3/g0_b3_union_shell_r64_raw.obj

build/base_no_mkl/bin/remeshSurface cgal_iso \
  --input-mesh examples/fbms/g0_b3/g0_b3_union_shell_r64_raw.obj \
  --output-mesh examples/fbms/g0_b3/g0_b3_union_shell_r64_remesh.obj \
  --edge-length 0.75 \
  --sharp-edge-angle 180
```

Generate b8 at resolution 64:

```bash
build/base_no_mkl/bin/generateFBMSUnionSurface \
  --fbms examples/fbms/g0_b8/g0_b8_fbms.obj \
  --sphere examples/fbms/g0_b8/g0_b8_fbms_bounding_sphere.obj \
  --fbms-thickness 0.02 \
  --sphere-thickness 0.02 \
  --resolution 64 \
  --padding-ratio 0.08 \
  --output-surface examples/fbms/g0_b8/g0_b8_union_shell_r64_raw.obj

build/base_no_mkl/bin/remeshSurface cgal_iso \
  --input-mesh examples/fbms/g0_b8/g0_b8_union_shell_r64_raw.obj \
  --output-mesh examples/fbms/g0_b8/g0_b8_union_shell_r64_remesh.obj \
  --edge-length 0.75 \
  --sharp-edge-angle 180
```

## Lightweight Thickened Assets

For weaker machines, prefer resolution 64 with thicker shells. These assets use `fbms-thickness=0.06` and `sphere-thickness=0.06`, which gives the SDF about 1.4 grid cells of half-thickness at resolution 64 and is much less fragmented than the `0.02` shell at the same resolution.

Generate b3 lightweight thickened assets:

```bash
build/base_no_mkl/bin/generateFBMSUnionSurface \
  --fbms examples/fbms/g0_b3/g0_b3_fbms.obj \
  --sphere examples/fbms/g0_b3/g0_b3_fbms_bounding_sphere.obj \
  --fbms-thickness 0.06 \
  --sphere-thickness 0.06 \
  --resolution 64 \
  --padding-ratio 0.08 \
  --output-surface examples/fbms/g0_b3/g0_b3_union_shell_r64_t006_raw.obj

build/base_no_mkl/bin/remeshSurface cgal_iso \
  --input-mesh examples/fbms/g0_b3/g0_b3_union_shell_r64_t006_raw.obj \
  --output-mesh examples/fbms/g0_b3/g0_b3_union_shell_r64_t006_remesh.obj \
  --edge-length 0.75 \
  --sharp-edge-angle 180
```

Generate b8 lightweight thickened assets:

```bash
build/base_no_mkl/bin/generateFBMSUnionSurface \
  --fbms examples/fbms/g0_b8/g0_b8_fbms.obj \
  --sphere examples/fbms/g0_b8/g0_b8_fbms_bounding_sphere.obj \
  --fbms-thickness 0.06 \
  --sphere-thickness 0.06 \
  --resolution 64 \
  --padding-ratio 0.08 \
  --output-surface examples/fbms/g0_b8/g0_b8_union_shell_r64_t006_raw.obj

build/base_no_mkl/bin/remeshSurface cgal_iso \
  --input-mesh examples/fbms/g0_b8/g0_b8_union_shell_r64_t006_raw.obj \
  --output-mesh examples/fbms/g0_b8/g0_b8_union_shell_r64_t006_remesh.obj \
  --edge-length 0.75 \
  --sharp-edge-angle 180
```

## Generated Asset Stats

Generated on 2026-04-24 with `build/base_no_mkl`.

| Asset | Size | Vertices | Faces | BBox min | BBox max |
| --- | ---: | ---: | ---: | --- | --- |
| `g0_b3/g0_b3_union_shell_raw.obj` | 54 MB | 520426 | 1040848 | `(-1.013398, -1.013399, -1.013399)` | `(1.013399, 1.013399, 1.013399)` |
| `g0_b3/g0_b3_union_shell_remesh.obj` | 105 MB | 1011764 | 2023524 | `(-1.013398, -1.013399, -1.013399)` | `(1.013399, 1.013399, 1.013399)` |
| `g0_b8/g0_b8_union_shell_raw.obj` | 60 MB | 584338 | 1168692 | `(-1.012541, -1.012328, -1.012476)` | `(1.012452, 1.012666, 1.012517)` |
| `g0_b8/g0_b8_union_shell_remesh.obj` | 118 MB | 1125209 | 2250434 | `(-1.012541, -1.012328, -1.012476)` | `(1.012452, 1.012666, 1.012517)` |
| `g0_b3/g0_b3_union_shell_r128_raw.obj` | 13 MB | 128682 | 257432 | `(-1.004152, -1.004152, -1.004153)` | `(1.004153, 1.004153, 1.004153)` |
| `g0_b3/g0_b3_union_shell_r128_remesh.obj` | 25 MB | 246177 | 492422 | `(-1.004145, -1.004147, -1.004146)` | `(1.004148, 1.004146, 1.004146)` |
| `g0_b8/g0_b8_union_shell_r128_raw.obj` | 14 MB | 144464 | 289040 | `(-1.003309, -1.003095, -1.003243)` | `(1.003219, 1.003433, 1.003285)` |
| `g0_b8/g0_b8_union_shell_r128_remesh.obj` | 28 MB | 274536 | 549184 | `(-1.003303, -1.003091, -1.003241)` | `(1.003215, 1.003431, 1.003279)` |
| `g0_b3/g0_b3_union_shell_r64_raw.obj` | 2.2 MB | 22572 | 46940 | `(-1.013003, -1.013003, -1.013003)` | `(1.013004, 1.013004, 1.013003)` |
| `g0_b3/g0_b3_union_shell_r64_remesh.obj` | 3.2 MB | 32594 | 66984 | `(-1.013003, -1.013003, -1.013003)` | `(1.013004, 1.013004, 1.013003)` |
| `g0_b8/g0_b8_union_shell_r64_raw.obj` | 2.5 MB | 25512 | 53220 | `(-1.012146, -1.011932, -1.012081)` | `(1.012057, 1.012270, 1.012122)` |
| `g0_b8/g0_b8_union_shell_r64_remesh.obj` | 3.6 MB | 36890 | 75976 | `(-1.012146, -1.011932, -1.012081)` | `(1.012057, 1.012270, 1.012122)` |
| `g0_b3/g0_b3_union_shell_r64_t006_raw.obj` | 2.8 MB | 28866 | 57728 | `(-1.031804, -1.031805, -1.031805)` | `(1.031806, 1.031805, 1.031805)` |
| `g0_b3/g0_b3_union_shell_r64_t006_remesh.obj` | 5.4 MB | 55157 | 110310 | `(-1.031804, -1.031805, -1.031805)` | `(1.031806, 1.031805, 1.031805)` |
| `g0_b8/g0_b8_union_shell_r64_t006_raw.obj` | 3.1 MB | 32154 | 64324 | `(-1.030972, -1.030758, -1.030906)` | `(1.030882, 1.031096, 1.030948)` |
| `g0_b8/g0_b8_union_shell_r64_t006_remesh.obj` | 6.0 MB | 60989 | 121994 | `(-1.030972, -1.030758, -1.030906)` | `(1.030882, 1.031096, 1.030948)` |

The remeshed OBJ files are intended as TetWild/fTetWild inputs for later `.veg` generation.
