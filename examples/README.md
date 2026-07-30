# libpgo examples

Example configuration and input data are intentionally separated:

- `assets/` contains each checked-in mesh, constraint file, and preview once.
- `configs/` contains sampled, IPC, and animation JSON files.
- `scripts/generate_cubic_veg.py` creates cubic VEG meshes.
- `scripts/generate_lite_tet_cubic_asset.py` creates dedicated lightweight
  tet/cubic smoke assets.
- `generated/` is ignored and contains local cubic meshes and simulation output.

Relative paths are resolved from the JSON file, so commands can be run from the
repository root.

## Generate cubic meshes

Build `cubicMesher`, then generate one or more presets:

```bash
cmake --build build/base --target cubicMesher
python3 examples/scripts/generate_cubic_veg.py \
  --build-dir build/base \
  --scene box \
  --scene bunny
```

Omit `--scene` to generate all five presets. `box-with-sphere` is the
full-resolution sampled showcase; `box-with-sphere-lite` preserves the small
resolution-5 IPC case. The script refuses to reuse an
existing output directory. It writes `<scene>-cubic.veg` and `manifest.json` to
`examples/generated/cubic/` by default. Use `--output-dir` for another fresh
directory.

For a custom surface, supply `--input-surface`, `--name`, `--resolution`, and
`--E`; `--nu` and `--density` are optional.

Run every retained simulation config for two timesteps, then convert two frames
through every retained animation config:

```bash
python3 examples/scripts/smoke_test_configs.py \
  --build-dir build/base
```

The smoke driver invokes both `tetMesher` and `cubicMesher` through
`generate_lite_tet_cubic_asset.py`. For each volume scene, that generator takes
the canonical surface bounds, writes an eight-vertex box surface, creates a
6-element tet mesh and a resolution-2 cubic mesh, and regenerates constraint
indices against those meshes. The materialized smoke configs use the lightweight
surface consistently for simulation and animation conversion. This retains each
config's domain, contact model, material, constraint, solver, and output route
without paying for showcase geometry. Simulation output must contain finite OBJ
vertices for the case to pass.

To generate only the smoke assets in a fresh directory:

```bash
python3 examples/scripts/generate_lite_tet_cubic_asset.py \
  --build-dir build/base \
  --output-dir /tmp/libpgo-smoke-assets
```

The standalone showcase generator retains each preset's full configured
resolution; it is not used by the bounded smoke matrix.

## Run examples

Sampled and IPC volume simulations remain separate entry points:

```bash
mkdir -p examples/generated/output

build/base/bin/runSim \
  examples/configs/volume/box/box-tet-sampled.json

build/base/bin/runIPCSim \
  examples/configs/volume/box/box-hang-tet-ipc.json

build/base/bin/runIPCSim \
  examples/configs/volume/box/box-hang-cubic-ipc.json
```

Cubic configs require the matching generated mesh first. The shell examples are:

```bash
build/base/bin/runShellSim \
  examples/configs/shell/shell-dynamic-sampled.json

build/base/bin/runIPCSim \
  examples/configs/shell/shell-dynamic-ipc.json
```

Simulation output is written below `examples/generated/output/`.

## Path migration

| Old path | New path |
| --- | --- |
| `examples/box/box.json` | `examples/configs/volume/box/box-tet-sampled.json` |
| `examples/cubic/box/box.json` | `examples/configs/volume/box/box-cubic-sampled.json` |
| `examples/ipc/tet/box-hang/box-ipc.json` | `examples/configs/volume/box/box-hang-tet-ipc.json` |
| `examples/ipc/cubic/box-hang/box-ipc.json` | `examples/configs/volume/box/box-hang-cubic-ipc.json` |
| `examples/ipc/cubic/box-with-sphere/box-ipc.json` | `examples/configs/volume/box-with-sphere/box-with-sphere-lite-cubic-ipc.json` |
| `examples/shell/shell.json` | `examples/configs/shell/shell-dynamic-sampled.json` |
| `examples/ipc/shell/shell-ipc.json` | `examples/configs/shell/shell-dynamic-ipc.json` |
| `examples/cubic/<scene>/<scene>.veg` | `examples/generated/cubic/<scene>-cubic.veg` |

There are no compatibility copies at the old paths.
