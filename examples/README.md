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

## CUDA IPC

The CUDA backend supports tetrahedral and eight-node cubic `stable-neo` and `stvk-vol` solids,
`koiter-stvk` shells, regularized Coulomb friction, embedded collision surfaces, static or backward-Euler solves,
gravity, pulling trajectories, prescribed external objects, Rayleigh mass/stiffness damping, and
the existing state/restart files. CPU remains the default. Select CUDA with
`"simulation-backend": "cuda"` and optionally `"cuda-device": 0`.

The [CUDA parameter reference](CUDA_PARAMETERS.md) lists every exposed Newton,
PCG, IPC and material-validity setting, including defaults, units, allowed ranges,
C++-only controls and fixed implementation choices.

Build with an NVIDIA CUDA toolkit and run the box example:

```bash
cmake --preset base_cuda_release -DCMAKE_CUDA_ARCHITECTURES=native
cmake --build --preset base_cuda_release --target runIPCSim
build/base_cuda_release/bin/runIPCSim examples/configs/volume/box/box-tet-ipc-cuda.json --log
```

For a build without GPU visibility, set an explicit CUDA architecture supported
by the intended GPU (for example `90` for H200). CUDA is optional; builds with
`PGO_ENABLE_CUDA=OFF` have no dependency on the selected libuipc sources.

The CUDA box config uses the original box mesh, embedded surface, obstacle,
material and solver tolerance, and runs 750 timesteps through impact. It writes
states and meshes to `examples/generated/output/box-tet-ipc-cuda/`.

Material evaluation, projected Hessian products, contact geometry and assembly,
CCD, and Jacobi-PCG run on the device. Host code controls Newton/backtracking,
updates prescribed inputs, and handles checkpoints. The contact law retains
libpgo's normalized log barrier and geometric weights. Only the required
libuipc distance/CCD, LBVH topology, eigenprojection and PCG kernels are imported; see
[the import manifest](../third-party/libuipc/IMPORTS.md).

When material max-step protection is enabled, cubic feasibility defaults to
`cuda-cubic-validity: "gauss"`, matching the CPU's eight FEM quadrature points.
This preserves the original material and checks positive determinants at those
points throughout each displacement chord. It does not guarantee positivity
between the quadrature points: corners can invert while Newton converges and
the embedded contact surface remains collision-free. These whole-cell failures
are reported separately by the trajectory checker.
Set `cuda-cubic-validity: "whole-cell"` to additionally require positive Jacobians
throughout each trilinear cell. Directed interval bounds and bounded subdivision
certify this stronger condition; an inconclusive bound shortens or rejects the
step. This adds a feasibility restriction, without adding restoring forces to
the material. Initial states must pass the selected enabled feasibility checks
before Newton can report convergence. Each run writes `material-validity.json`
to record its criterion.
Restarts retain this metadata and reject changes to the criterion. A cubic
restart imported from an older run without metadata requires an explicit
`cuda-cubic-validity` setting; earlier states must be checked independently.

The backend uses GPU Morton LBVH queries with outward-rounded double AABBs,
checked candidate-buffer growth, and cached element/contact Hessians. Stronger
preconditioners remain follow-up work. Other materials and deterministic assembly
are rejected explicitly. Use the CPU backend for those options. Optional `cuda-pcg-max-iter` and
`cuda-pcg-tolerance` control the inner solver; `solver-eps` still measures the
Newton gradient infinity norm. The final log reports contact counts, convergence,
Newton/PCG iterations and elapsed time.
Line search evaluates component energy changes directly so small corrections
remain measurable after large rigid motion. After 20 Newton iterations, the
solver can try the unprojected tangent; PCG curvature and descent checks fall
back to the projected tangent when needed. The per-frame CSV records these
attempts, accepted directions, and fallbacks.
`cuda-exact-hessian-start-iteration` changes when these attempts begin (default
`20`; a negative value disables them). This changes the search direction,
while the original energy, gradient tolerance and collision checks still govern
acceptance.
`cuda-pcg-batch-size` optionally keeps PCG reduction scalars on the GPU and polls
after a bounded batch (`1`–`64`; default `0` uses host scalars). Each iteration
still checks curvature, finite values and residual tolerance. A terminal result
freezes the PCG state, and the final true residual is checked before Newton uses
the direction.
`cuda-stagnation-window` optionally detects a sequence of rounded, tiny updates
whose gradient norms stagnate after physical-chord clipping (`0` disables it).
The solver first tries bounded residual recovery; otherwise it reports a failed
solve so the dynamic runner can retry a smaller physical step. This early-failure
heuristic never declares convergence. The attempt CSV records its count and
full-precision step and residual ratios.
`cuda-stagnation-tolerance` sets the relative-change threshold for both gradient
norms in that optional detector (default `64 * epsilon(double)`). It must be
finite and strictly between zero and one. A larger value can detect plateaus
whose residual evaluations vary slightly despite tiny represented motion. This
setting affects early failure and recovery only; `solver-eps` remains the
convergence threshold. Record any explicit value with the solver strategy when
comparing trajectories.
If rounding prevents a near-stationary Newton direction from reducing the
residual, a bounded coordinate correction minimizes the linearized residual
using exact Hessian columns. It is accepted only within the existing energy
roundoff bound, with a measured residual reduction and certified motion. The
`residual_correction_*` CSV counters expose these attempts; the original Newton
tolerance still determines convergence.
Set `cuda-solver-verbose` to `2` to include per-iteration step and line-search
diagnostics when investigating convergence.

For timing a bounded invocation, `runIPCSim` accepts `--max-frames N` and
`--timing-csv PATH` on either backend. The frame limit preserves the configured
`num-timestep` and loading schedules, including after a checkpoint restart.
Timing rows separate compute from state/mesh exports and record convergence and
commit status. Compute includes runner work and rejected attempts; startup and
final invocation checks remain part of process wall time. The paired benchmark
helper is `scripts/benchmark_cuda_examples.py`.

GPU tests are labeled `cuda`; run them on a machine with a visible GPU:

```bash
cmake --build --preset base_cuda_release --target cudaNewtonSolver_gtest cudaDeformationModelEnergy_gtest cudaKoiterDeformationModelEnergy_gtest cudaEmbeddedSurfaceIPCPotentialEnergy_gtest cudaSurfaceIPCFrictionPotentialEnergy_gtest cudaIPCSim_gtest
ctest --test-dir build/base_cuda_release -L cuda --output-on-failure
```

### Validating complete CUDA examples

The validation criterion is that **every timestep reaches its configured Newton
gradient tolerance and the trajectory remains collision-free**. A successful
process exit, finite states, or a plausible animation is insufficient. The CUDA
runner writes full-precision per-timestep diagnostics to `newton-solves.csv` in
the output directory; this file covers the current invocation, including original
frame indices after a restart. Verbosity `2` also records individual Newton
iterations and line-search steps in the run log.

For dynamic simulation, the CUDA runner also checks motion from the physical
step's starting state to each Newton trial. If a dynamic interval cannot converge with certified motion,
it retries that interval with twice as many backward Euler substeps. The
configured timestep remains the output interval and maximum integration step;
`cuda-max-timestep-refinements` limits the number of halvings (default `8`, range
`0`–`12`). Every accepted physical substep must be certified. An exhausted retry
budget fails without committing that frame. The actual trajectory follows all
recorded substeps; a direct shortcut between coarse frame endpoints may cross
the mesh even when that physical path is clear. The frame CSV retains the
separate `macroframe_chord_certified` diagnostic. Do not interpolate such a
shortcut in a downstream consumer.
`newton-attempts.csv` records accepted and rejected solves, including each friction
lagging iteration and its frozen pair count; `newton-substeps.csv`
records accepted physical times and intermediate state files. Matching surfaces
are exported for every accepted step as `substeps/frame####-step####.obj`, in
addition to the configured coarse mesh dumps. The renderer selects these actual
states at uniform physical times without interpolating their geometry. The frame CSV's
gradient is the maximum over all accepted solves, including friction lagging
iterations, and iteration counts sum those solves. Energy and final
linear-solver/step diagnostics refer to the last solve.
`max_substep_newton_iterations` is the largest individual solve count. Every
lagging iteration keeps the physical step's reference, inertia, and damping
fixed; only the final solution advances time. Each solve retains the configured
tolerance and iteration cap.
The independent checker validates every recorded physical segment, endpoint,
and prescribed-object movement. It reports coarse shortcuts separately. A
missing or invalid row in a present physical manifest fails validation; legacy
outputs without a manifest use their macroframe chords as the trajectory.

Static equilibrium keeps the full configured load throughout Newton. Its
`static-newton-path.json` records the input and every accepted optimization
iterate for each friction lagging solve, with states in `newton-path/`.
Every actual rounded Newton transition is certified on the GPU and independently
checked by the CPU validator. Every terminal solve must meet the configured
tolerance; intermediate iterates are not physical timesteps or converged states.
This explicit path is authoritative for static validation, with the direct
initial-to-final chord reported separately as a diagnostic. The static
`newton-substeps.csv` row is terminal-state metadata. Missing or conflicting
entries in either present manifest fail validation. Static videos show the
final equilibrium and do not assign physical times to Newton iterates.

`scripts/validate_cuda_examples.py` prepares isolated copies of the original box,
box-with-sphere, and bunny configurations by default; `--volume-ipc` selects all
11 canonical volumetric IPC configurations, including cubic, dragon, static
contact and damped dynamic contact. `--shell-ipc` selects the four canonical
hanging-cloth, shirt and towel configurations. Preparation preserves their physics,
tolerances, timestep counts, and mesh dump intervals. Its `prepare`, `run`, and
`summarize` commands accept `--destination` and `--cases`; prepared configurations
contain an explicit `cuda-device` that can be adjusted for the available GPUs.
Explicit PCG accuracy/batching and exact-tangent start options are recorded as numerical
strategies in each report; the source material, contact law, Newton tolerance,
iteration cap and time schedule must still match.
Summarizing produces an HTML report, CSV tables, convergence plots, and links to
the raw logs and videos. It exits nonzero while any selected case fails or is
incomplete.
Frame iteration totals count accepted solves. `newton-attempt-work.csv` also
records completed rejected attempts from the verbose log, so their Newton and
PCG work is included when assessing runtime. Interrupted attempts remain in the
raw log and are excluded from these completed-attempt totals.

Build `validateIPCTrajectory` alongside `runIPCSim`, then independently check each
completed trajectory on the CPU:

```bash
build/base_cuda_release/bin/validateIPCTrajectory path/to/config.json --swept --report path/to/collision-validation.json
```

The checker examines every saved state for self/obstacle intersections and
nonpositive element Jacobians. Tetrahedra use determinant extrema. Cubic checks
follow the explicitly configured or recorded validity criterion: Gauss-point
validity matches the CPU model, while whole-cell validity additionally uses
outward-rounded interval bounds and tensor Bernstein subdivision over the
entire cell. Whole-cell diagnostics are retained when Gauss-point validity is
selected. Use `--cubic-validity gauss` or `--cubic-validity whole-cell` to select
the checking criterion explicitly. Historical runs without a configured or
recorded criterion retain whole-cell checking. With `--swept`, these
checks cover the actual piecewise linear trajectory along with independent CPU
CCD. A physical-substep manifest always requires these swept checks, including
when `--swept` is omitted.
In Gauss mode, the report separately counts unique elements with confirmed
negative Jacobian witnesses and lists their IDs and affected frames. Cells
whose whole-cell check is inconclusive are reported separately; a failed
positivity certificate alone does not establish an inversion. These counts
describe the recorded path, excluding diagnostic coarse-frame shortcuts.
Passing cubic determinant certificates report conservative lower bounds; failed
certificates are labeled as uncertified diagnostics. A segment
on the actual path that CCD cannot certify is reported as a validation failure, separately from a
detected intersection. `scripts/render_cuda_validation.py`, run with a Python
environment containing Blender's `bpy`, renders the saved OBJ frames and configured
obstacles to MP4; its overlays state simulation time and playback speed.
For a partial diagnostic, use the checker's `--first-frame` and `--last-frame`
options. These retain the original prescribed-motion timing; editing
`num-timestep` to inspect a prefix can change that timing.

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

## Shell material

Shell simulations use `elastic-material: koiter-stvk` with these keys:

| Key | Meaning | Default |
| --- | --- | --- |
| `shell-youngs-modulus` | membrane modulus E, Pa | 1e6 |
| `shell-poisson-ratio` | membrane Poisson ratio, in (-1, 0.5) | 0.4 |
| `shell-bending-youngs-modulus` | bending modulus E_b, Pa | E |
| `shell-bending-poisson-ratio` | bending Poisson ratio, in (-1, 0.5) | `shell-poisson-ratio` |
| `shell-thickness` | h, metres | 3e-3 |
| `shell-areal-density` | kg/m2, scales the lumped triangle mass matrix | 100 |
| `shell-zero-rest-curvature` | flat rest curvature instead of the loaded mesh's | false |

The model converts each modulus and Poisson ratio to 3D Lame constants and
evaluates a Koiter energy on twice the Green strain of the first fundamental
form, and on the change of the second fundamental form scaled by h^3 / 12.
The config values are therefore nominal: the effective in-plane Young's
modulus is 4 E / (1 - nu^2), the effective in-plane Poisson ratio is
nu / (1 - nu), and the bending rigidity of a cylindrical bend is
E_b (1 - nu_b) / ((1 + nu_b)(1 - 2 nu_b)) h^3 / 12. The
[`towel_grasp`](configs/shell/towel_grasp/README.md) README works an example
from measured cloth data.

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

## Linear solver

Both runners accept optional solver keys:

```json
"linear-solver": "default",
"pardiso-deterministic": true,
"pardiso-thread-limit": 64
```

`linear-solver` selects the direct solver behind every Newton solve:
`default` prefers Panua PARDISO when the build includes it, then Intel MKL
PARDISO, then Eigen's LDLT; `pardiso`, `mkl-pardiso` and `eigen` request one
explicitly and fail with a message when it is not part of the build. A build
with Panua PARDISO has no MKL PARDISO option, because both libraries export
the same `pardiso` entry point.

With Panua PARDISO, `pardiso-deterministic` (default true) applies the
user guide's recipe for an identical solution independent of the number of
processors (IPARM(34) = 1 with METIS reordering, the two-level parallel
factorization and solve, and sequential METIS reordering), and
`pardiso-thread-limit` (default 64) is the guide's IPARM(3), the upper limit of
processors for which the results coincide. Any `OMP_NUM_THREADS` up to that
limit then produces bit-identical solves; if the runtime would use more
threads than the limit, the limit is raised to the team size, since a lower
value corrupts the factorization. Together with `deterministic-assembly`,
runs with different thread counts write identical state files, which the
runIPCSim tests check.

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
