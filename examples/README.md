# libpgo examples

This directory is the entry point for libpgo simulation examples. It summarizes
the supported mesh, contact, and solve-mode combinations and links to the
available scenes. Scene-specific setup, commands, outputs, and reference
results belong in the corresponding config directory.

## Capability coverage

| Mesh | Sampled dynamic | Sampled static | IPC dynamic | IPC static |
| --- | --- | --- | --- | --- |
| Tet | supported | solver only; contact ignored | supported | supported |
| Cubic | supported | solver only; contact ignored | supported | supported |
| Shell | unsupported | unsupported | supported | supported |

Static sampled volume configurations exercise the elastic, attachment, and
gravity solve, but sampled contact is not added to their Newton energy. The
runner prints an explicit warning. Sampled shell simulation is not supported.

## Example index

| Scene | Meshes | Contact models | Modes | Coverage |
| --- | --- | --- | --- | --- |
| [`box`](configs/volume/box/) | tet, cubic | sampled, IPC | dynamic | Small baseline volume simulation |
| [`box-contact`](configs/volume/box-contact/README.md) | tet, cubic | sampled, IPC | dynamic, static | External plane contact and Static IPC/Dynamic IPC comparison |
| [`box-with-sphere`](configs/volume/box-with-sphere/) | tet, cubic | sampled, IPC | dynamic | External-object contact |
| [`bunny`](configs/volume/bunny/) | tet, cubic | sampled, IPC | dynamic | Nonzero initial displacement |
| [`dragon-dyn`](configs/volume/dragon-dyn/) | tet, cubic | sampled, IPC | dynamic | Translated falling object |
| [`dragon`](configs/volume/dragon/README.md) | tet, cubic | sampled | static | Static volume solve without contact |
| [`shell_hang`](configs/shell/shell_hang/README.md) | shell | IPC | dynamic, static | Shell FEM through the IPC runner |

Every simulation JSON has a sibling `*-animation.json`. Cubic configurations
depend on meshes generated locally with
[`scripts/generate_cubic_veg.py`](scripts/generate_cubic_veg.py).

## Directory map

- [`assets/`](assets/) contains checked-in meshes and constraint files.
- [`configs/volume/`](configs/volume/) contains tet and cubic configurations.
- [`configs/shell/`](configs/shell/) contains shell configurations.
- [`scripts/`](scripts/) contains example asset-generation utilities.
- `assets/generated/cubic/` contains locally generated cubic VEG meshes and
  their extracted surface OBJ meshes.
- `generated/` contains local OBJ sequences and Alembic output.

Both generated directories are ignored by Git. Relative paths in a JSON file
are resolved from that JSON file, so commands can be launched from the
repository root.

## Common commands

The checked-in [`cases.json`](cases.json) manifest and
[`BATCH_RUNNER.md`](BATCH_RUNNER.md) provide one-command execution for all
examples or a selected set of named cases.

```bash
uv run pgo-run-cases examples/cases.json --list
uv run pgo-run-cases examples/cases.json --list-groups
uv run pgo-run-cases examples/cases.json box-contact-static-tet-ipc
uv run pgo-run-cases examples/cases.json --group box-contact
uv run pgo-run-cases examples/cases.json
```

Run any supported sampled or IPC simulation through the `pypgo` CLI:

```bash
uv run pgo-run-sim <simulation-config.json>
```

After the OBJ sequence has been generated, convert a sibling animation config
to Alembic. Create the destination directory first:

```bash
mkdir -p <abc-output-directory>
uv run pgo-dump-abc <animation-config.json> <abc-output-directory>
```

Follow the README inside a scene config directory when it provides additional
setup instructions or reference results.
