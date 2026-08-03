# IPC Static External Object MVP Implementation Plan

## 1. Goal

Add static kinematic `external-objects` to the IPC runner while preserving the
existing deformable self-contact path.

The MVP shall:

- merge deformable and external triangle meshes into one collision surface;
- keep the existing `PTPair` and `EEPair` representations, without adding a
  separate `TPPair` type;
- use a dynamic/static vertex mask to retain deformable-deformable and
  deformable-external pairs while rejecting external-external pairs;
- include retained pairs in the existing IPC energy, gradient, Hessian, and CCD
  maximum-step computations;
- keep external vertices out of the simulation unknowns by assigning them zero
  rows in the surface displacement map;
- accept only finite `external-objects[].movement == [0, 0, 0]`;
- support nonzero finite `init-disp`, including dynamic-restart fallback
  semantics; and
- migrate the bunny, dynamic dragon, and box-with-sphere IPC examples from the
  temporary analytic floor to their sampled-case external OBJ assets.

## 2. Non-goals

The MVP will not implement:

- moving kinematic objects;
- interpolation of external-object transforms over frames;
- CCD against nonzero external-object velocity;
- rotation, rigid transforms, or per-object scale;
- friction;
- solving external-object degrees of freedom;
- contact between two external objects; or
- a distinct external-contact energy or pair type.

`restart-from-u` is a dynamic-simulation feature in this MVP. Static solves use
`init-disp` only as their Newton initial guess and do not load a prior `.u`
state. A static config with `restart-from-u=true` is rejected with a clear
configuration error rather than silently ignoring the field.

The existing analytic floor remains supported. A config may contain a floor and
external objects at the same time, although the migrated showcase cases will
use only their external OBJ floor geometry.

## 3. Mathematical model and invariants

Let the deformable collision surface be `(Vd, Fd)` and its mapping from
simulation displacement `u` be `Wd`. Let all external meshes be concatenated as
`(Ve, Fe)`.

Build one collision surface:

```text
Vcollision = [ Vd ]
             [ Ve ]

Fcollision = [ Fd                  ]
             [ Fe + vertex offsets ]

Wcollision = [ Wd ]
             [  0 ]
```

The collision positions and trial displacements are:

```text
xcollision(u)  = Vcollision + Wcollision * u
dxcollision(du) = Wcollision * du
```

Consequently, external vertices remain fixed and the existing mapped energy
produces:

```text
gradient_u = Wcollision^T * gradient_collision
Hessian_u  = Wcollision^T * Hessian_collision * Wcollision
```

Required invariants:

1. The deformable vertices occupy the first `Vd.rows()` collision vertices.
2. Every face index is valid after concatenation.
3. Every triangle and edge belongs entirely to either deformable or external
   geometry; the builder never creates mixed primitives.
4. `Wcollision.rows() == 3 * Vcollision.rows()`.
5. `Wcollision.cols() == simulationDOFs`.
6. Rows corresponding to external vertices have no nonzero coefficients.
7. Output/render geometry continues to contain only the deformable surface.
8. Existing callers that do not provide a mask behave as if every vertex is
   dynamic.
9. The mask-aware embedded-energy API rejects a vertex marked static when any
   coefficient in the corresponding three rows of `Wcollision` is nonzero.

## 4. Pair policy

Do not add `TPPair`. A point-triangle pair is represented by the existing
`PTPair` regardless of which mesh owns the point or triangle.

Store a per-vertex dynamic mask in `SurfaceIPCTopology` and derive primitive
activity from it:

```text
triangleIsDynamic(fi) = any vertex of triangle fi is dynamic
edgeIsDynamic(ei)     = any endpoint of edge ei is dynamic
```

Apply the same policy in the barrier broad phase and CCD broad phase:

```text
accept PT(p, fi)  iff vertexIsDynamic[p] || triangleIsDynamic(fi)
accept EE(ei, ej) iff edgeIsDynamic(ei) || edgeIsDynamic(ej)
```

This gives the following behavior:

| Pair ownership | Result |
| --- | --- |
| deformable-deformable | retained as self-contact |
| deformable point-external triangle | retained as `PTPair` |
| external point-deformable triangle | retained as `PTPair` |
| deformable edge-external edge | retained as `EEPair` |
| external-external, same object | rejected |
| external-external, different objects | rejected |

Shared-vertex adjacency rejection and all existing distance/AABB tests remain
unchanged and run after the ownership filter.

## 5. Planned code changes

### 5.1 Topology and mask

Files:

- `src/core/contact/ipc/topology/surfaceIPCTopology.h`
- `src/core/contact/ipc/topology/surfaceIPCTopology.cpp`

Changes:

- Add a compact per-vertex dynamic mask.
- Preserve `setMesh(V, F)` and make it initialize an all-dynamic mask.
- Add `setMesh(V, F, vertexIsDynamic)`.
- Validate that the mask length equals `V.rows()`.
- Provide small query helpers for vertex, triangle, and edge activity.
- Rebuild all derived activity state whenever `setMesh` is called.

The default overload is required so all current self-contact tests and users
retain their behavior without source changes.

### 5.2 Barrier broad-phase filtering

File:

- `src/core/contact/ipc/broadPhase/surfaceIPCSelfBroadPhase.cpp`

Changes:

- Before evaluating a PT candidate, reject it when both the point and triangle
  are external.
- Before evaluating an EE candidate, reject it when both edges are external.
- Do not split candidate lists by contact direction and do not add new pair
  containers.

The spatial hash can continue to index the combined geometry. The first MVP
optimizes for a small, reviewable change; separate dynamic/external hash tables
are a later performance option only if profiling requires them.

### 5.3 CCD maximum-step filtering

File:

- `src/core/contact/ipc/core/surfaceIPCMaxStep.cpp`

Changes:

- Apply exactly the same PT and EE ownership filters as the barrier broad
  phase.
- Keep external trial displacement at zero through `Wcollision * du`.
- Ensure external-external candidates cannot reduce the Newton feasible step.

Using the same topology helpers in both paths prevents the active-pair and CCD
policies from diverging.

### 5.4 Core and embedded-energy API plumbing

Files:

- `src/core/contact/ipc/core/surfaceIPCCore.h`
- `src/core/contact/ipc/core/surfaceIPCCore.cpp`
- `src/core/contact/embeddedSurfaceIPCPotentialEnergy.h`
- `src/core/contact/embeddedSurfaceIPCPotentialEnergy.cpp`

Changes:

- Preserve the current all-dynamic `setMesh(V, F)` path.
- Add a mask-aware `SurfaceIPCCore::setMesh` overload.
- Add a mask-aware `EmbeddedSurfaceIPCPotentialEnergy` constructor overload.
- Pass the mask into `SurfaceIPCTopology`.
- In the mask-aware embedded-energy constructor, validate the cross-object
  invariant that every static vertex has three identically zero rows in the
  supplied displacement map. Reject an inconsistent mask/map pair before any
  energy or CCD evaluation. Dynamic vertices are allowed to have zero rows.
- Keep the barrier assembler unchanged.
- Add only the minimal read-only inspection needed by tests, if existing public
  state is insufficient to verify the combined surface.

`MappedSurfacePotentialEnergy` requires no algorithm change: it already maps
surface energy, gradient, Hessian, and maximum-step inputs through an arbitrary
sparse `W` with the correct formulas.

### 5.5 Collision-surface builder and config validation

Files:

- `src/simulationRunner/runIPCSimSetup.cpp`
- `src/simulationRunner/runIPCSimSetup.h` only if a reusable test-visible data
  structure is justified

Add an internal builder with responsibilities equivalent to:

```cpp
CollisionSurfaceData buildCollisionSurface(
  const MXd &deformableV,
  const MXi &deformableF,
  const SpMatD &deformableW,
  const ConfigFileJSON &config,
  double scale);
```

The returned data contains combined vertices, offset faces, sparse displacement
map, and dynamic mask.

Builder procedure:

1. Copy deformable `V/F` and mark its vertices dynamic.
2. Resolve each `external-objects[].filename` relative to the config file.
3. Parse `movement` as exactly three numeric values.
4. Reject non-finite values.
5. Reject any value other than zero, including a useful object index and path in
   the error.
6. Require the global simulation `scale` to be finite and strictly positive.
7. Load each OBJ through `TriMeshGeo`, remove vertices not referenced by any
   triangle, and apply the global simulation `scale`, matching the sampled
   runner.
8. Reject an external mesh that has no vertices or no triangles after cleanup,
   contains an invalid/repeated face index, contains a non-finite vertex, or
   contains a zero-area triangle after scaling. Include the object index and
   resolved path in the error.
9. Count total vertices/faces and allocate combined dense matrices once.
10. Append vertices and append faces with checked vertex offsets.
11. Construct `Wcollision` sparsely by copying only `Wd` triplets into its top
   rows; do not densify `W`.
12. Mark appended external vertices static.
13. Validate final dimensions, indices, and the static-mask/zero-row invariant
    before constructing the IPC energy.

Both shell and volume IPC setup should use the same builder. With no
`external-objects` field or with an empty array, it must reproduce the current
all-dynamic collision surface.

Both setup paths support finite positive non-unit `scale`. Geometry transforms
must be ordered as follows:

```text
Vd_scaled = scale * Vd_input
Ve_scaled = scale * Ve_input
Vcollision_scaled = [ Vd_scaled ]
                    [ Ve_scaled ]
xcollision(initialU) = Vcollision_scaled + Wcollision * initialU
```

In particular, scale the shell rest/simulation mesh about the input origin
before constructing its mass, elastic model, embedding, and collision surface;
then apply `init-disp` as a simulation-space translation. Do not multiply
`init-disp` by `scale`, and never apply it to external vertices. Remove the
current shell-only `scale ~= 1` assertion.

The `IpcSimulationContext` keeps its current deformable-only:

- `surfaceMesh`;
- `surfaceRestPositions`; and
- `surfaceFromSimulationDispMap`.

Those fields drive output generation. Only `collisionHandler` receives the
combined `Vcollision/Fcollision/Wcollision/mask`, so `retXXXX.obj` never gains
external geometry.

### 5.6 Nonzero `init-disp`

Files:

- `src/simulationRunner/runIPCSimSetup.cpp`
- `src/simulationRunner/ipcSimulationRunner.cpp`

Changes:

- Remove the zero-only `init-disp` rejection.
- Parse exactly three finite components.
- Expand the translation to every simulation vertex:

  ```text
  initialU.segment<3>(3 * vi) = initDisp
  ```

- Use `initialU` as the dynamic simulation's initial displacement.
- Use `initialU` as the static Newton initial guess.
- For a dynamic simulation with `restart-from-u=true`, scan restart candidates
  from newest to oldest. A valid dynamic state has exactly `simulationDOFs`
  rows, at least three columns, and finite displacement, velocity, and
  acceleration columns. The first valid candidate overrides `initialU`,
  configured initial velocity, and zero initial acceleration.
- If dynamic restart is requested but no state file exists, fall back to
  `initialU`, configured initial velocity, and zero initial acceleration.
- If a restart file exists but is malformed, has incompatible dimensions, or
  contains non-finite state, fail with a path-specific error instead of silently
  treating data corruption as a missing restart.
- Static solves never load `deform0000.u`; reject `restart-from-u=true` for
  `sim-type=static`.
- Never apply `init-disp` to external vertices; they are not simulation DOFs.

Initialization must occur before the first IPC energy/CCD evaluation so the
deformable-external relative placement is correct from the first solve.

### 5.7 Example config migration

Replace temporary floor fields in these six IPC simulation configs:

- `examples/configs/volume/bunny/bunny-tet-ipc.json`
- `examples/configs/volume/bunny/bunny-cubic-ipc.json`
- `examples/configs/volume/dragon-dyn/dragon-dynamic-tet-ipc.json`
- `examples/configs/volume/dragon-dyn/dragon-dynamic-cubic-ipc.json`
- `examples/configs/volume/box-with-sphere/box-with-sphere-tet-ipc.json`
- `examples/configs/volume/box-with-sphere/box-with-sphere-cubic-ipc.json`

Use the sampled-case values:

| Case | External asset | `init-disp` |
| --- | --- | --- |
| bunny tet/cubic | `bottom-large.obj` | `[0, -0.3, 0]` |
| dragon-dyn tet/cubic | `bottom.obj` | `[-1.2, 0.6, 0]` |
| box-with-sphere tet/cubic | `bottom-large.obj` | `[0, 0, 0]` |

For each config:

- remove `use-floor`, `floor-axis`, `floor-height`, and `floor-kappa`;
- add `external-objects` with zero movement; and
- restore the sampled case's `init-disp`.

The IPC animation configs reference deformable output sequences and require no
schema change. Validate their sequence paths after the simulation config edits.
The static dragon IPC case remains unchanged because its external-object list is
empty.

## 6. Test plan

### 6.1 Topology validation

File:

- `tests/src/core/contact/surfaceIPCTopology_gtest.cpp`

Add coverage for:

- omitted mask means all vertices dynamic;
- explicit mixed dynamic/static vertex mask is stored correctly;
- invalid mask length is rejected; and
- calling `setMesh` again rebuilds mask-derived state.

### 6.2 Broad-phase pair policy

File:

- `tests/src/core/contact/surfaceIPCSelfBroadPhase_gtest.cpp`

Use two close, disconnected triangles so ownership is unambiguous:

- first triangle dynamic and second external: retain PT pairs in both point
  ownership directions and retain cross-mesh EE pairs;
- both triangles external: produce no PT or EE pairs;
- both triangles dynamic: match the existing all-dynamic pair set.

Assertions should inspect pair vertex ownership rather than introducing a TP
label.

### 6.3 CCD ownership policy

File:

- `tests/src/core/contact/surfaceIPCMaxStep_gtest.cpp`

Add coverage for:

- a dynamic triangle moving through a static external triangle produces a
  feasible step strictly between zero and one;
- the same geometric motion marked entirely external does not reduce the step;
  and
- helper and `SurfaceIPCCore` maximum-step results still agree with a mask.

### 6.4 Mapped IPC derivatives with zero external rows

File:

- `tests/src/core/contact/embeddedSurfaceIPCPotentialEnergy_gtest.cpp`

Construct a small combined surface with:

- one dynamic triangle;
- one nearby external triangle;
- simulation DOFs only for the dynamic triangle; and
- zero `W` rows for the external triangle.

Verify:

- finite nonzero contact energy;
- simulation-space gradient against central finite differences;
- simulation-space Hessian against gradient differences, using the
  non-projected Hessian mode where required for consistency;
- Hessian dimensions equal the dynamic simulation DOF count; and
- external vertices remain unchanged when simulation displacement changes.

Also construct an inconsistent mask/map pair with a nonzero coefficient in a
static vertex row and verify that the mask-aware
`EmbeddedSurfaceIPCPotentialEnergy` constructor rejects it. Keep a companion
case showing that a dynamic vertex is allowed to have zero mapping rows.

### 6.5 Runner setup and config errors

File:

- `tests/src/tools/runIPCSim_gtest.cpp`

Add setup tests for shell and volume paths where practical:

- no external field retains current behavior;
- one static external OBJ loads successfully;
- multiple external OBJ meshes concatenate successfully;
- global `scale` is applied to deformable and external geometry consistently;
- shell and volume geometry are scaled about the input origin before
  `init-disp` is applied;
- movement `[0, 0, 0]` is accepted;
- nonzero movement is rejected with a useful message;
- NaN/infinite movement is rejected;
- missing or unreadable external OBJ is rejected;
- empty/face-free external OBJ is rejected;
- external OBJ isolated vertices are removed and cannot participate in CCD;
- repeated-index and zero-area external triangles are rejected;
- deformable output mapping dimensions remain unchanged; and
- combined collision mapping has the original simulation column count.

Add initialization tests for:

- nonzero dynamic `init-disp` reaches the first saved state/output;
- nonzero static `init-disp` is accepted as the initial guess;
- a valid dynamic restart state overrides `init-disp` and configured velocity;
- missing dynamic restart state falls back to `init-disp` and configured
  velocity;
- malformed, wrong-sized, and non-finite dynamic restart states fail with a
  useful error; and
- static `restart-from-u=true` is rejected.

### 6.6 Integration and pypgo smoke tests

Modify `tests/pypgo/test_pgo_smoke.py` so the materialized tet and cubic IPC
fixtures include a small static external OBJ, for example a two-triangle plane
written into `tmp_path` at the existing analytic floor height. Keep the current
`use-floor`/floor-energy fields in these smoke fixtures so analytic-floor
compatibility continues to be exercised; the focused native tests remain
responsible for proving external-only contact response.

Run those tet and cubic IPC fixtures through the existing shared runner paths.
Require:

- process/API return value is success;
- IPC energy, Hessian, and CCD paths execute;
- state matrices contain only finite values;
- output OBJ contains the deformable vertex/face counts, not combined collision
  counts; and
- at least one nonzero-`init-disp` fixture is exercised.

Additionally verify from the fixture/config and output that the external OBJ
was accepted and was not appended to `retXXXX.obj`. This prevents the existing
analytic floor from allowing the smoke matrix to pass without exercising the
external-object loading and combined-mapping path.

Then run the six migrated showcase configs only at deliberately shortened
timesteps using temporary materialized copies. Do not edit showcase timestep
counts for smoke testing.

## 7. Implementation order

1. Add topology mask storage, validation, and all-dynamic compatibility.
2. Apply mask filtering to barrier broad phase and add focused pair tests.
3. Apply the same filter to CCD and add maximum-step tests.
4. Plumb the mask through `SurfaceIPCCore` and
   `EmbeddedSurfaceIPCPotentialEnergy`.
5. Add the collision-surface builder and zero-movement validation.
6. Connect the builder to both shell and volume IPC setup without changing
   deformable output state.
7. Add zero-row mapping derivative tests and setup error tests.
8. Implement scale-before-translation initialization, nonzero `init-disp`, and
   dynamic-only restart precedence tests.
9. Migrate the six IPC configs and validate all referenced files.
10. Run focused tests, the complete IPC native suite, pypgo smoke tests, and
    shortened showcase runs.

Each step should keep the tree buildable and its focused tests passing before
moving to the next step.

## 8. Validation commands

Use the repository's configured build directory and exact target names. The
expected validation layers are:

1. Build the contact library and affected test targets.
2. Run topology, broad-phase, max-step, embedded-energy, and runner setup
   gtests.
3. Run all `RunIPC*` CTest entries.
4. Run the existing pypgo tet/cubic IPC smoke matrix.
5. Materialize shortened copies of the six showcase configs in a temporary
   directory and run them through pypgo.
6. Parse every edited JSON file and resolve every referenced asset.
7. Run `git diff --check`.

Long showcase simulations remain a manual validation tier and are not required
for the local MVP development loop.

## 9. Acceptance criteria

The MVP is complete only when all of the following hold:

- Existing all-dynamic IPC self-contact tests still pass unchanged.
- Combined collision topology produces deformable-deformable and
  deformable-external PT/EE pairs, with no external-external pairs.
- No `TPPair` type or external-specific barrier assembler is introduced.
- External vertices have zero displacement for every simulation state.
- The mask-aware embedded-energy API rejects static vertices with nonzero map
  rows.
- External contact contributes finite energy, simulation-space gradient, and
  simulation-space Hessian.
- CCD limits deformable motion against static external geometry.
- Nonzero or non-finite external movement is rejected before simulation.
- Shell and volume rest geometry are scaled about the input origin before
  `init-disp` is applied.
- Nonzero finite `init-disp` works with the documented dynamic-only restart
  precedence, and static restart requests are rejected.
- Output OBJ files contain only deformable geometry.
- The six target IPC configs use the same external assets and initial
  displacement as their sampled counterparts.
- Focused native tests, the IPC regression suite, and the pypgo smoke matrix
  pass.
- No unrelated user changes are overwritten.

## 10. Future extension point

Moving external objects should build on this representation rather than replace
it. A later phase can make the external base positions time-dependent and feed
external swept displacement into CCD. That work must define frame interpolation,
Newton/line-search time semantics, restart reconstruction, and tunneling
prevention before allowing nonzero `movement`.
