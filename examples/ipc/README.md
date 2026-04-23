# IPC Examples

This directory contains runnable examples for `runIPCSim`, along with the per-case animation configs used by `convertAnimation`.

The current examples are meant to cover the three IPC entry points that exist in this repo today:

1. shell self-contact IPC
2. tetrahedral volume unified IPC
3. cubic volume unified IPC
4. conversion of dumped frame sequences into Alembic caches

## Quick Start

Build the tools:

```bash
cmake --preset base_no_mkl_debug
cmake --build build/base_no_mkl_debug --target runIPCSim convertAnimation
```

Run the three shipped cases from the repo root:

```bash
build/base_no_mkl_debug/bin/runIPCSim examples/ipc/shell/shell-ipc.json
build/base_no_mkl_debug/bin/runIPCSim examples/ipc/tet/box-hang/box-ipc.json
build/base_no_mkl_debug/bin/runIPCSim examples/ipc/cubic/box-hang/box-ipc.json
```

Convert the dumped frames into Alembic:

```bash
build/base_no_mkl_debug/bin/convertAnimation examples/ipc/shell/anim.json
build/base_no_mkl_debug/bin/convertAnimation examples/ipc/tet/box-hang/anim.json
build/base_no_mkl_debug/bin/convertAnimation examples/ipc/cubic/box-hang/anim.json
```

The JSON configs use paths relative to the config file, so they can be launched from the repo root without first changing into the case directory.

## What Each Case Contains

The directory currently ships three runnable IPC inputs:

- `examples/ipc/shell/shell-ipc.json`
- `examples/ipc/tet/box-hang/box-ipc.json`
- `examples/ipc/cubic/box-hang/box-ipc.json`

The shell case lives directly under `examples/ipc/shell/`. The two volume cases currently live under:

- `examples/ipc/tet/box-hang/`
- `examples/ipc/cubic/box-hang/`

The two volume cases each contain:

- the display surface mesh: `box.obj`
- the volumetric mesh: `box.veg`
- the fixed-vertex list: `box-fixed.txt`
- the `runIPCSim` config: `box-ipc.json`
- the `convertAnimation` config: `anim.json`

The shell case contains:

- the display mesh: `shell.obj`
- the fixed-vertex list: `shell-fixed0.txt`
- the `runIPCSim` config: `shell-ipc.json`
- the `convertAnimation` config: `anim.json`
- previously generated output folders and a committed `shell-ipc.abc`

Current config convention:

- shell uses `surface-mesh` together with `elastic-material = koiter-stvk`
- tet uses `tet-mesh` together with `surface-mesh`
- cubic uses `cubic-mesh` together with `surface-mesh`
- tet and cubic provide explicit `ipc-dhat` and `ipc-kappa`
- shell currently uses `ipc-heuristic: true`
- in tet and cubic, `fixed-vertices` refers to volume simulation vertex indices, not surface vertex indices

## Case Guide

### `shell`

Existing shell self-contact IPC example. This is the legacy shell-only setup under `examples/ipc`, and it is the only case in this directory that already includes committed historical output folders and a generated `.abc` cache.

- files: `shell.obj`, `shell-fixed0.txt`, `shell-ipc.json`, `anim.json`
- material: `koiter-stvk`
- IPC params: shell heuristic via `ipc-heuristic: true`
- run: `build/base_no_mkl_debug/bin/runIPCSim examples/ipc/shell/shell-ipc.json`
- output: `examples/ipc/shell/ret-shell-ipc/`
- animation: `build/base_no_mkl_debug/bin/convertAnimation examples/ipc/shell/anim.json`
- Alembic: `examples/ipc/shell/shell-ipc.abc`

### `tet/box-hang`

Minimal tetrahedral unified IPC hanging-box case. This is the compact tet sanity check for the `runIPCSim` volume path: it pins the vertices listed in `box-fixed.txt`, uses `stable-neo`, and writes a display-surface OBJ sequence through the embedded IPC mapping.

- files: `box.obj`, `box.veg`, `box-fixed.txt`, `box-ipc.json`, `anim.json`
- material: `stable-neo`
- IPC params: explicit `ipc-dhat = 0.002`, `ipc-kappa = 3000.0`
- config note: `fixed-vertices` means tet simulation vertex indices
- run: `build/base_no_mkl_debug/bin/runIPCSim examples/ipc/tet/box-hang/box-ipc.json`
- output: `examples/ipc/tet/box-hang/ret-box-ipc/`
- animation: `build/base_no_mkl_debug/bin/convertAnimation examples/ipc/tet/box-hang/anim.json`
- Alembic: `examples/ipc/tet/box-hang/box-hang-ipc-tet.abc`

### `cubic/box-hang`

Minimal cubic unified IPC hanging-box case. This case mirrors the tet setup, but drives contact and output from a cubic volumetric simulation mesh.

- files: `box.obj`, `box.veg`, `box-fixed.txt`, `box-ipc.json`, `anim.json`
- material: `stable-neo`
- IPC params: explicit `ipc-dhat = 0.002`, `ipc-kappa = 3000.0`
- config note: `fixed-vertices` means cubic volume simulation vertex indices
- run: `build/base_no_mkl_debug/bin/runIPCSim examples/ipc/cubic/box-hang/box-ipc.json`
- output: `examples/ipc/cubic/box-hang/ret-box-ipc/`
- animation: `build/base_no_mkl_debug/bin/convertAnimation examples/ipc/cubic/box-hang/anim.json`
- Alembic: `examples/ipc/cubic/box-hang/box-hang-ipc-cubic.abc`

## Rebuilding Animation Outputs

After `runIPCSim` dumps the OBJ sequence, use the per-case animation config to convert it:

```bash
build/base_no_mkl_debug/bin/convertAnimation examples/ipc/shell/anim.json
build/base_no_mkl_debug/bin/convertAnimation examples/ipc/tet/box-hang/anim.json
build/base_no_mkl_debug/bin/convertAnimation examples/ipc/cubic/box-hang/anim.json
```

If the optional output path is omitted, `convertAnimation` writes the `.abc` file next to the animation config and uses the mesh `name` field from that config as the output filename stem.
