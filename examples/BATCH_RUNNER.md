# Example batch runner

[`cases.json`](cases.json) registers every checked-in simulation case and its
matching animation configuration. All paths are resolved relative to the batch
manifest, not the current working directory.

## Manifest format

```json
{
  "version": 1,
  "abc-output-root": "generated/abc",
  "groups": {
    "box": ["box-tet-ipc", "box-tet-sampled"]
  },
  "cases": [
    {
      "name": "box-tet-ipc",
      "simulation-config": "configs/volume/box/box-tet-ipc.json",
      "animation-config": "configs/volume/box/box-tet-ipc-animation.json"
    },
    {
      "name": "custom-output-example",
      "simulation-config": "configs/volume/box/box-tet-sampled.json",
      "animation-config": "configs/volume/box/box-tet-sampled-animation.json",
      "abc-output": "generated/abc/custom-directory"
    }
  ]
}
```

`version` and `cases` are required. Case names must be unique. `groups` maps a
group name to one or more registered case names. An optional per-case
`abc-output` overrides `abc-output-root/<case-name>`.

## Commands

List all registered case names:

```bash
uv run pgo-run-cases examples/cases.json --list
uv run pgo-run-cases examples/cases.json --list-groups
```

For every case, run its simulation and then export its animation:

```bash
uv run pgo-run-cases examples/cases.json
```

Run selected cases in the order given on the command line:

```bash
uv run pgo-run-cases examples/cases.json \
  box-contact-static-tet-ipc \
  box-contact-dynamic-tet-ipc
```

Run every case in a scene group:

```bash
uv run pgo-run-cases examples/cases.json --group box-contact
```

`--group` can be repeated, and groups can be combined with explicit case
names. Cases appearing more than once are run only once:

```bash
uv run pgo-run-cases examples/cases.json \
  shell-dynamic-ipc \
  --group box \
  --group bunny
```

Run or export only one stage:

```bash
uv run pgo-run-cases examples/cases.json --stage simulation
uv run pgo-run-cases examples/cases.json --stage animation
```

Use `--dry-run` to print the equivalent `pgo-run-sim` and `pgo-dump-abc`
commands. By default the batch stops at the first failed case. Add
`--keep-going` to run the remaining cases and return a failing status after the
batch summary.

The animation stage creates each configured ABC output directory. Simulation
output locations remain owned by the individual simulation JSON files.

The cubic cases require their generated VEG assets. Generate them before an
all-case run:

```bash
python3 examples/scripts/generate_cubic_veg.py --build-dir build/base
```
