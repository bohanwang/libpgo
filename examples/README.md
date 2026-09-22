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

IPC supports regularized Coulomb friction for self-contact and kinematic
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
external vertices are fixed within a frame, so friction acts relative to the
object's current placement. CCD continues to limit collision-free steps.

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

## IPC barrier distance

Without `ipc-heuristic`, the barrier distance comes from exactly one of

```json
"ipc-dhat": 0.002
"ipc-dhat-relative": 0.002
```

`ipc-dhat` is an absolute distance in simulation units. `ipc-dhat-relative`
multiplies the diagonal of the deformable surface's bounding box (measured
after `scale`), which keeps dhat below the mesh's own feature size regardless
of the model's size. That matters: IPC only excludes primitives that share a
vertex from self-contact, so a dhat comparable to the local edge length turns
neighboring triangles into permanent barrier pairs. The bunny example, whose
bounding-box diagonal is about 0.25 with triangle heights down to 1.2 mm, had
6037 such pairs at frame 0 with `ipc-dhat = 0.002`; it now uses
`ipc-dhat-relative = 0.002` (dhat about 5e-4) and starts contact-free.
`ipc-kappa` is always given explicitly and is not rescaled with dhat.

## IPC external objects

`ipc-external-objects` lists kinematic triangle meshes that the simulated
surface collides with. Each entry takes `filename`, an optional positive
`scale`, an optional `initial-translation`, and an optional `movement`:

```json
"ipc-external-objects": [
  { "filename": "../../../assets/common/bottom.obj", "scale": 1.0,
    "initial-translation": [0, 0, 0], "movement": [0, 0.05, 0] }
]
```

The objects follow the same model as `runSim`'s `external-objects`. They have
no degrees of freedom and take no forces. Inside a frame they are frozen; the
runner places them between frames at `initial + movement * k / (N - 1)` for
frame `k` of an `N`-frame run, so the full movement is reached at the last
frame, a one-frame run keeps the initial placement, and a restarted run sees
the same placements. Static solves keep the initial placement. The barrier,
CCD, and friction all act between the simulated surface and the objects.
Because the object is static within a frame, friction resists motion relative
to its current placement and does not drag the simulated surface along.
Objects are never tested against themselves or each other, so an external
mesh may self-intersect or overlap another one.

Each placement must leave the simulated surface outside the objects. A jump
that lands an object inside the simulated surface is reported through the
start-state check below and stops the run; keep per-frame movement below the
current gap. A jump that passes completely through the surface cannot be
detected by a state check, since the object is not swept between frames.
Deformation of an object with fixed triangle topology is supported through
`EmbeddedSurfaceIPCPotentialEnergy::updateExternalSurface`; the config exposes
translation only.

## Prescribed vertex trajectories

A `fixed-vertices` entry can follow a recorded trajectory instead of a
constant `movement`:

```json
"fixed-vertices": [
  { "trajectory": "../../../assets/shell/towel_grasp/grasp-trajectory.json", "coeff": 1e4 }
]
```

An entry takes either `filename` plus `movement` or `trajectory`, never both;
`coeff` is the attachment stiffness in either case. The trajectory file is
JSON:

```json
{
  "format": "pgo-vertex-trajectory",
  "version": 1,
  "space": "position",
  "interpolation": "linear",
  "tracks": [
    { "name": "left-gripper", "vertices": [1221, 1222], "release-time": 11.4,
      "keys": [
        { "t": 0.00, "p": [[0.5678, 0.1394, 0.8329], [0.5647, 0.1548, 0.8308]] },
        { "t": 0.02, "p": [[0.5677, 0.1395, 0.8328], [0.5646, 0.1549, 0.8307]] } ] }
  ]
}
```

A track groups vertices that share keyframes, so a gripper holding two
vertices is one rigid clamp. Each key has the simulation time `t` in seconds
and one position per track vertex, in `vertices` order; the speed between keys
is therefore explicit and is logged per track. `space` is `position`
(coordinates in the mesh file's frame before `scale`, like the surface mesh)
or `displacement` (offsets from rest, also before `scale`). Only linear
interpolation exists. Times must increase strictly, a vertex may appear in one
track only, and vertex indices follow the fixed-vertex text files.

At frame k the target is the trajectory at t = k * `timestep`, the same rule
as the interpolated `movement`. Sampling holds the first key before it and the
last key after it. From `release-time` on, the track's vertices are free; the
attachment stays in the system with zero stiffness, so a restarted run makes
the same decisions. Static solves use the last key and ignore release times.
`runSim` rejects trajectory entries; they are an IPC runner feature.

## Initial surface state

`initial-surface-mesh` names an OBJ with the topology of `surface-mesh` whose
positions are the deformable surface at t = 0. The runner turns it into a
per-vertex initial displacement (after `scale`) and adds the uniform
`init-disp`. It is available for shell simulations; tet and cubic simulations
start from `init-disp` or a restart state. The
[`towel_grasp`](configs/shell/towel_grasp/README.md) example starts from a
folded towel this way.

## IPC start-state check

Before the first step, when restarting from a state file, and after every
placement of a moving external object, the IPC runner validates the collision
surface at the current displacement. It
rejects any point-triangle or edge-edge pair at zero distance, and any edge
that passes through a triangle it does not share a vertex with, using the exact
segment-triangle predicate. The second test is necessary: two triangles that
pierce each other have strictly positive point-triangle and edge-edge
distances, so the barrier energy is zero there and CCD then prevents the
surfaces from ever separating. States produced by a sampled-contact run, by an
overlapping external object, or by a self-intersecting input mesh are reported
with the offending edge and triangle instead of being simulated. Pairs whose
vertices are all static are ignored.

## Deterministic assembly

Parallel assembly is not bit-reproducible by default: element and contact-pair
contributions are accumulated with atomic adds and per-thread partial sums, and
the IPC broad phase lists pairs in the order the worker threads produced them,
so two identical runs differ in the last digits and the differences grow over a
simulation. Set

```json
"deterministic-assembly": true
```

to make the elastic assembly (tet, cubic and shell elements) and the IPC
contact assembly bit-reproducible across runs, independent of thread count and
scheduling. Contributions are still computed in parallel; they are written to
per-element or per-pair buffers and then summed per destination in a fixed
order, the IPC pairs are sorted canonically after every broad phase, and the
friction and Hessian terms follow that order. The cost is the extra buffers and
one gather pass per assembly. With this flag the runIPCSim tests compare
restarted and uninterrupted runs exactly, and two runs of the same
configuration write identical state files. The flag is accepted by both
runners; the sampled contact energies are not covered by it.

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
| [`towel_grasp`](configs/shell/towel_grasp/README.md) | shell | IPC | dynamic | Recorded two-arm grasp trajectory, folded initial state, release onto a table |

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

The converter fails with a nonzero exit code, instead of writing an empty or
mistimed archive, when a frame in `sequence-range` is missing, a frame's vertex
count differs from the driving mesh, the `display-mesh` cannot be loaded, or
the `sequence-type` is unknown. Supported types are `objmesh` (one OBJ per
frame), `u` (one displacement matrix file per frame, `.u` format) and `uall`
(a single matrix whose columns are the frames). With `save-cache` (default 1)
the driving displacements are stored next to the config as `<name>-uAll.u` and
reused on later conversions; delete that file after regenerating the frames,
because a cache whose size does not match the sequence is rejected.

Follow the README inside a scene config directory when it provides additional
setup instructions or reference results.
