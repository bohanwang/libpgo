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

## IPC friction

IPC supports regularized Coulomb friction for self-contact and stationary
external triangle meshes, on all three simulation mesh types. Set:

```json
"contact-friction-coeff": 0.3,
"ipc-friction-epsv": 0.001,
"ipc-friction-iterations": 1
```

The coefficient defaults to zero, preserving frictionless behavior in existing
configurations. `ipc-friction-epsv` is the positive tangential speed threshold
in simulation length units per second (default `0.001`). The displacement
threshold used by backward Euler is `timestep * ipc-friction-epsv`.
The shirt-drop example enables friction with the values above.

`ipc-friction-iterations` is a positive number of Newton solves per time step
(default `1`). Closest points, tangent planes, and barrier normal forces are
frozen for each solve. With one iteration they come from the previous time
step; newly formed contacts contribute friction starting in the following
step. Extra iterations refresh these quantities using the preceding solution,
while keeping the time-step reference displacement and inertia fixed.
Two to four iterations can improve the lagged approximation at extra cost.

Static friction is evaluated over one load increment, measured from
`init-disp`, using the configured positive `timestep` for regularization.
Use at least two iterations to include contacts first formed during the static
solve. This does not implement a history of quasi-static load increments.

The implementation follows the [IPC paper, section 5](https://ipc-sim.github.io/file/IPC-paper-350ppi.pdf),
[technical supplement, sections 8–9](https://ipc-sim.github.io/file/IPC-supplement-A-technical.pdf),
and [IPC Toolkit's smoothing functions](https://github.com/ipc-sim/ipc-toolkit/blob/main/src/ipc/friction/smooth_friction_mollifier.cpp).
It uses the same cubic smoothing law as this repository's sampled contact,
with normal force obtained from the weighted IPC barrier. Point–triangle and
edge–edge closest-feature cases are supported. The zero-slip Hessian uses its
analytic finite limit. Friction is mapped through the existing surface embedding;
external vertices remain fixed. CCD continues to limit collision-free steps.

Backward Euler evaluates its inertial/damping quadratic in the displacement
increment from the previous state. This is equivalent to the original energy
up to an additive constant, with the same forces and Hessian, but avoids
subtracting large nearly equal terms when resolving small frictional steps.

To run the full manifest with friction enabled and validate every state
and OBJ file, first build `runIPCSim`, `runSim`, and `cubicMesher`, then run:

```bash
python3 examples/scripts/generate_cubic_veg.py --build-dir build/<your-build>
python3 examples/scripts/test_contact_demos.py \
  --build-dir build/<your-build> \
  --output-dir examples/generated/friction-validation \
  --cores 64 --jobs 8 --friction-coeff 0.3 --friction-epsv 1e-4 \
  --static-sampled-ipc-companions
```

This Linux runner uses disjoint CPU affinities and writes isolated configs,
logs, per-case `result.json` and `newton-solves.csv` files, and a live `REPORT.md`
and `progress.json`. The coefficient applies to both IPC and sampled demos;
`--friction-epsv` sets their velocity regularization using the appropriate
configuration field. Newton's convergence tolerance is unchanged unless
`--solver-eps` is also supplied. `--max-steps N` explicitly requests a shorter
smoke run; the default runs every configured step. Use a fresh output directory
for each suite.

The solver reports its stop reason and recomputes the gradient at the returned
state. A solve has converged when the infinity norm of that gradient is below
`solver-eps`. `NewtonSolver::solve` returns `SOLVE_CONVERGED`,
`SOLVE_NOT_CONVERGED` (iteration limit, tiny step, or failed line search with a
finite state), or `SOLVE_NUMERICAL_FAILURE` (non-finite energy, gradient, or
step). The runners stop with a nonzero exit code on a numerical failure and log
a warning for each solve that stops short of the tolerance, so process success
still does not establish convergence of every step.
The two original static sampled dragon demos exercise no contact, as explained
in the capability table. `--static-sampled-ipc-companions` adds two separately
named IPC tests of those scenes using the dynamic dragon example's barrier
parameters, while retaining the original sampled runs.

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
| [`shirt-drop`](configs/shell/shirt_drop/) | shell | IPC | dynamic | Garment falling onto a ground plane |

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
