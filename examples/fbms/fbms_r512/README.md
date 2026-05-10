# FBMS All Asset Pipeline

Generates tet meshes for FBMS and volume-matched TPMS baselines, ready for simulation.

## Configuration

All parameters in `pipeline_r512.json`. Each stage can be re-run independently with `--from-stage` and `--overwrite`.

```json
"bounding_sphere_bounds_method": "unit"   // unit sphere at origin for both FBMS and TPMS
"keep_largest_components": 3              // keep only the 3 largest edge-connected components after raw extraction
"fbms_thickness": 0.02                    // fixed for FBMS; negative value triggers volume-matching search for TPMS
```

## Pipeline

```
prepare → raw → remesh → tet → validate (→ package)
                           ├── noshell: raw → remesh (CGAL → geogram fallback)
```

### prepare
- Generates a unit sphere per FBMS case (`generate_bounding_sphere.py --bounds-method unit`)
- Generates TPMS isosurfaces with matching unit spheres (`generate_tpms_unit_ball.py`)

### raw
- OpenVDB CSG intersection of thickened FBMS shell and sphere shell → union surface
- TPMS thickness searched to match FBMS raw volume
- Filters to keep only the 3 largest edge-connected components

### remesh
- CGAL isotropic remeshing
- **Noshell**: if CGAL remesh has self-intersections, automatically falls back to geogram at the same vertex count

### tet
- fTetWild volumetric meshing

### validate
- Checks boundary surface quality (topology, winding, self-intersections)
- Auto-repairs failing boundaries via `repair_tpms_mesh.py`
- Compares `veg.obj` and `repaired.obj` mtime to avoid reusing stale repair outputs

### noshell
- Runs from completed union pipeline outputs, extracts the `union-minus-sphere` surface (FBMS structure with sphere shell removed)
- Separate entry point: `run_fbms_noshell_from_union.py`

### package
- Copies passing `.veg` and surface `.obj` into `simulation_package/`
- Uses repaired surface when `boundary_source=repaired`

## Usage

```bash
# Full pipeline
python scripts/run_fbms_all_asset_pipeline.py --config examples/fbms/fbms_r512/pipeline_r512.json --overwrite

# Single case, raw through tet
python scripts/run_fbms_all_asset_pipeline.py --config ... --case g0_b15 --from-stage raw --overwrite

# Noshell
python scripts/run_fbms_noshell_from_union.py --config ... --case g0_b15 --include all --overwrite

# Package
python scripts/package_fbms_sim_assets.py --config ... --case g0_b15 --overwrite
```

Child scripts are invoked via `_resolve_python()`, which prefers `.venv/bin/python` over `sys.executable` to ensure dependencies are available.

## Notes

- Changing the sphere or component filter requires re-running from raw (`--from-stage raw --overwrite`)
- Changing remesh parameters requires re-running from remesh
- g0_b5 noshell FBMS and I-WP remesh targets trigger the geogram fallback (slower but self-intersection-free)

<!-- r512-pipeline-results:start -->
| asset | raw volume | selected thickness | remesh mean edge | remesh status | tet vertices | tet elements | tet volume | boundary | components |
| --- | ---: | ---: | ---: | --- | ---: | ---: | ---: | --- | ---: |
| g0_b5 | 0.517673306 | 0.0200000 | 0.018445 | pass | 74063 | 259219 | 0.519454 | pass | 3 |
| g0_b5 / tpms_schwarz_p | 0.517668314 | 0.0227347 | 0.018475 | warn: 32 self-intersections | 73296 | 263930 | 0.520424 | pass(repaired) | 3 |
| g0_b5 / tpms_schwarz_d | 0.517655600 | 0.0199523 | 0.018457 | pass | 73586 | 257948 | 0.519742 | pass | 3 |
| g0_b5 / tpms_gyroid | 0.517662403 | 0.0204869 | 0.018457 | pass | 74043 | 261344 | 0.519884 | pass | 3 |
| g0_b5 / tpms_iwp | 0.517585839 | 0.0174371 | 0.018453 | pass | 76131 | 264696 | 0.519992 | pass | 3 |
| g0_b5 / tpms_neovius | 0.517633508 | 0.0165015 | 0.018466 | pass | 78922 | 276095 | 0.520960 | pass | 3 |
| g0_b10 | 0.545841892 | 0.0200000 | 0.018450 | pass | 78949 | 274877 | 0.548030 | pass | 3 |
| g0_b10 / tpms_schwarz_p | 0.545741224 | 0.0273154 | 0.018489 | warn: 1 self-intersections | 73273 | 268371 | 0.548594 | pass | 3 |
| g0_b10 / tpms_schwarz_d | 0.545799922 | 0.0239619 | 0.018458 | pass | 75940 | 269557 | 0.547925 | pass | 3 |
| g0_b10 / tpms_gyroid | 0.545702813 | 0.0245937 | 0.018454 | pass | 75346 | 268515 | 0.547995 | pass | 3 |
| g0_b10 / tpms_iwp | 0.545746313 | 0.0209486 | 0.018461 | pass | 78079 | 273237 | 0.547952 | pass | 3 |
| g0_b10 / tpms_neovius | 0.545746739 | 0.0198186 | 0.018468 | pass | 78313 | 274253 | 0.548947 | pass | 3 |
| g0_b15 | 0.559187290 | 0.0200000 | 0.018455 | pass | 80694 | 280093 | 0.561605 | pass | 3 |
| g0_b15 / tpms_schwarz_p | 0.559044202 | 0.0294782 | 0.018454 | warn: 51 self-intersections | 76753 | 286997 | 0.561785 | pass(repaired) | 3 |
| g0_b15 / tpms_schwarz_d | 0.559122678 | 0.0258573 | 0.018454 | pass | 76508 | 273546 | 0.561211 | pass | 3 |
| g0_b15 / tpms_gyroid | 0.559169268 | 0.0265621 | 0.018464 | pass | 75644 | 272372 | 0.561386 | pass | 3 |
| g0_b15 / tpms_iwp | 0.559151992 | 0.0226253 | 0.018463 | pass | 78748 | 276354 | 0.561429 | pass | 3 |
| g0_b15 / tpms_neovius | 0.559102337 | 0.0213981 | 0.018472 | pass | 80018 | 281948 | 0.562312 | pass | 3 |
<!-- r512-pipeline-results:end -->
