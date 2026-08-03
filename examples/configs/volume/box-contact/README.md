# Box contact comparison

This scene compresses the tet box against a static plane. It provides three
configurations with the same deformable mesh, top attachment target, and plane:

| Configuration | Contact model | Solve mode | Purpose |
| --- | --- | --- | --- |
| `box-contact-static-tet-ipc.json` | IPC | static | Solve the final constrained state directly |
| `box-contact-dynamic-tet-ipc.json` | IPC | dynamic | Slowly approach the same state and compare against Static IPC |
| `box-contact-dynamic-tet-sampled.json` | sampled penalty | dynamic | Compare penalty contact against Dynamic IPC |

The box starts above the plane at `y = 0.36`. Its four top corners are driven
downward by `0.1`. The two dynamic configurations use 1001 time steps and dump
every tenth step, producing 101 OBJ frames. The external plane does not move.

## Run the simulations

Run these commands from the repository root:

```bash
uv run pgo-run-sim \
  examples/configs/volume/box-contact/box-contact-static-tet-ipc.json

uv run pgo-run-sim \
  examples/configs/volume/box-contact/box-contact-dynamic-tet-ipc.json

uv run pgo-run-sim \
  examples/configs/volume/box-contact/box-contact-dynamic-tet-sampled.json
```

The OBJ results are written to:

| Configuration | OBJ output |
| --- | --- |
| Static IPC | `examples/generated/output/box-contact-static-tet-ipc/ret0000.obj` |
| Dynamic IPC | `examples/generated/output/box-contact-dynamic-tet-ipc/ret0000.obj` through `ret0100.obj` |
| Dynamic sampled | `examples/generated/output/box-contact-dynamic-tet-sampled/ret0000.obj` through `ret0100.obj` |

## Reference results

The reference run produced the following final surface results:

| Result | Static IPC | Dynamic IPC | Dynamic sampled |
| --- | ---: | ---: | ---: |
| Surface vertices | 194 | 194 | 194 |
| Minimum `y` | `0.36192673` | `0.36192658` | `0.35996874` |
| Relation to the `y = 0.36` plane | separated by `1.93e-3` | separated by `1.93e-3` | penetrated by `3.13e-5` |

Comparing corresponding vertices in the final Static IPC and Dynamic IPC
meshes gives a maximum position difference of `2.40e-4` and an RMS difference
of `4.23e-5`. The dynamic result is expected to be close rather than bitwise
identical because it retains small inertial and damping effects. The sampled
penalty result permits a small penetration, as expected for a penalty model.

## Export Alembic

Each simulation config has a sibling animation config. After generating the OBJ
sequences, export the deformable box and static plane to Alembic with:

```bash
mkdir -p \
  examples/generated/abc/box-contact/static-ipc \
  examples/generated/abc/box-contact/dynamic-ipc \
  examples/generated/abc/box-contact/dynamic-sampled

uv run pgo-dump-abc \
  examples/configs/volume/box-contact/box-contact-static-tet-ipc-animation.json \
  examples/generated/abc/box-contact/static-ipc

uv run pgo-dump-abc \
  examples/configs/volume/box-contact/box-contact-dynamic-tet-ipc-animation.json \
  examples/generated/abc/box-contact/dynamic-ipc

uv run pgo-dump-abc \
  examples/configs/volume/box-contact/box-contact-dynamic-tet-sampled-animation.json \
  examples/generated/abc/box-contact/dynamic-sampled
```

Each output directory contains one ABC for the deformable box and one for the
static plane. All files under `examples/generated/` are local generated
artifacts and are intentionally ignored by Git.
