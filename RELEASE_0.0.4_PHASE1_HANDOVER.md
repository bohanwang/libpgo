# libpgo 0.0.4 Phase 1 Handover

## Snapshot

- Branch: `release/0.0.4`
- Handover commit: `727a037 clean up ipc core`
- Worktree status at handover: clean
- Phase 1 implementation is represented by the following commits:

  | Commit | Scope |
  | --- | --- |
  | `d9884ae` | P1.1: remove obsolete `CIPCPotentialEnergy` |
  | `ffb496a` | P1.2: remove scoped profiling |
  | `128de09` | P1.3/P1.4: Hessian contract and solver symbolic-analysis handling |
  | `f97421f` | P1.5: simplify max-step contract |
  | `727a037` | P1.8–P1.10, IPC/deformation FD coverage, IPC cleanup |

P1.6 and P1.7 were explicitly left unchanged by product decision; see the table below.

## Important decisions and contracts

These are deliberate choices made during Phase 1. A follow-up should preserve them unless the API contract is intentionally changed.

1. **Dynamic-topology `createHessian()` is invalid and throws.**  Dynamic IPC energy must not create a reusable symbolic Hessian. `PotentialEnergies::init()` skips Hessian creation for dynamic terms and the aggregate/direct path rebuilds as needed. This intentionally differs from the earlier proposal to return an empty matrix.
2. **IPC Hessian PSD projection is now optional.**  `SurfaceIPCCore::Parameters::projectHessianToPSD` defaults to `true`, preserving previous solver-facing behavior. Set it to `false` for the raw, energy-consistent barrier Hessian used by exact Hessian-vs-gradient FD tests.
3. **Mixed aggregate topology is dynamic.**  A `PotentialEnergies` aggregate with any dynamic child, including fixed deformation plus IPC, must be handled as dynamic for symbolic-factorization reuse.
4. **P1.6/P1.7 are intentionally not implementation work.**  Cubic hardening changes were reverted when it was decided that matching Tet behavior was not worthwhile; Tet/Cubic model unification was also declined.
5. **P1.9 is an audit, not a feature expansion.**  Existing optional MKL/original-Pardiso CMake branches were retained. No ARPACK, benchmark, or additional experimental path was introduced.

## P1.1–P1.10 status

| Item | Status | Result |
| --- | --- | --- |
| P1.1 Remove obsolete CIPC energy | Done | Deleted `CIPCPotentialEnergy`; IPC goes through Embedded/Map/core layers. |
| P1.2 Remove scoped profiling | Done | Removed the scoped profiling instrumentation. |
| P1.3 Fix Hessian interface contract | Done | Fixed and dynamic topology are distinguished; dynamic `createHessian()` remains an exception path by design. |
| P1.4 Fixed/dynamic symbolic analysis | Done | Fixed patterns reuse solver symbolic state; dynamic/mixed patterns trigger the correct reanalysis. Removed the incorrect fixed-path solver reset. |
| P1.5 Simplify max step | Done | Default energy max step is `1`; aggregate takes the minimum; deformation inversion and IPC CCD constraints remain. |
| P1.6 Cubic hardening | Intentionally not changed | Initial change was reverted. Existing cubic behavior/tests stay as before. |
| P1.7 Tet/Cubic merge | Intentionally not changed | No model hierarchy/API refactor. |
| P1.8 Put IPC validation in core | Done | Core owns mesh, parameter, state, and gradient-buffer validation. |
| P1.9 Optional-path retention audit | Done | Optional linear-solver branches remain; no unrelated optional feature added. |
| P1.10 Small cleanup | Done | Header cleanup plus focused tests; no broad refactor. |

## Main implementation details

### Hessian topology and Newton behavior

- `PotentialEnergies::init()` builds a symbolic Hessian only for fixed-topology children.
- A dynamic IPC adapter reports dynamic topology. Directly asking it to `createHessian()` is a contract violation and throws.
- The fixed topology Newton path no longer resets the solver after each solve, allowing symbolic factorization reuse.
- A fixed deformation energy combined with dynamic IPC is treated as dynamic at the aggregate level.

Relevant tests include dynamic-create failure, fixed+dynamic aggregation, zero-coefficient behavior, and Newton solver behavior in `tests/src/core/newtonSolver_gtest.cpp`.

### Max-step behavior

- `PotentialEnergy::computeMaxStepSize()` has the neutral default `1`.
- `PotentialEnergies` computes the minimum scalar step across enabled terms.
- The only material constraints kept are element inversion protection and IPC continuous-collision detection (CCD).

### IPC core validation and raw Hessian option

`SurfaceIPCCore` now validates:

- parameters: finite positive `dhat`/`kappa`, finite nonnegative `eps_ee`, finite `slackness` in `(0, 1]`;
- mesh: non-empty finite `N x 3` vertices, `M x 3` triangles, valid distinct indices;
- lifecycle: `setMesh()` before evaluation;
- state: correct `3N` size and finite surface positions/displacements;
- output: correctly sized gradient buffers.

The IPC adapter delegates mesh validation to the core. The mapped-surface layer additionally rejects non-finite rest vertices and non-finite embedding coefficients.

`projectHessianToPSD` is propagated from `SurfaceIPCCore::Parameters` into `SurfaceIPCBarrierAssembler::{computeHessian,computeAll}`. Keep it enabled for normal Newton/solver usage; disable it only when an exact raw Hessian is required for derivative verification.

### FD and transition coverage added

The commit `727a037` adds or expands the following tests.

- Tet deformation model: Neo-Hookean positional gradient/Hessian and plastic mixed derivatives; Hill material derivative coverage in `tests/src/core/solidDeformationModel/tetMeshDeformationModel_gtest.cpp`.
- Deformation assembler: single-tet and single-cubic directional energy-to-gradient and Hessian-vector-to-gradient FD checks.
- IPC geometry: all point-triangle dispatch classes and all edge-edge dispatch classes receive gradient/Hessian FD coverage.
- IPC barrier assembly: raw (unprojected) Hessian FD checks, `eps_ee > 0` mollifier coverage, PT face-to-edge/vertex and EE interior-to-endpoint mode-transition cases.
- IPC core/adapter: raw full Hessian-vs-gradient FD, including an identity and a coupled non-identity embedding map; default PSD behavior; `computeAll` consistency.
- Dynamic pair/pattern behavior: pair/pattern rebuilding around the cutoff, direct recomputation for changing states, and cutoff continuity checks.
- IPC max step: PT/EE CCD safe-step behavior in addition to existing core/adapter consistency coverage.

At an exact feature-mode switch, gradient FD is the required contract. Hessian FD is evaluated on each smooth side of the switch; the two one-sided Hessians are not required to match.

## Verification already performed

The following Phase-1-relevant targets were rebuilt in `build/base_no_mkl`, and the matching CTest group passed:

```bash
conda run --no-capture-output -n libpgo-dev cmake --build build/base_no_mkl --target \
  eigenSupport contact newtonSolver_gtest tetMeshDeformationModel_gtest \
  cubicMeshDeformationModel_gtest deformationModelAssembler_gtest \
  deformationModelEnergyMaxStep_gtest ipcGeometry_gtest surfaceIPCMaxStep_gtest \
  surfaceIPCBarrierAssembler_gtest surfaceIPCCore_gtest \
  embeddedSurfaceIPCPotentialEnergy_gtest embeddedSurfaceFloorPotentialEnergy_gtest \
  surfaceIPCTopology_gtest surfaceIPCSelfBroadPhase_gtest -- -j8

conda run --no-capture-output -n libpgo-dev ctest --test-dir build/base_no_mkl -R \
  '(NewtonSolverGTest|TetMeshDeformationModelGTest|CubicMeshDeformationModelGTest|DeformationModelAssemblerGTest|DeformationModelEnergyMaxStepGTest|IPCGeometryGTest|SurfaceIPCMaxStepGTest|SurfaceIPCBarrierAssemblerGTest|SurfaceIPCCoreGTest|EmbeddedSurfaceIPCPotentialEnergyGTest|EmbeddedSurfaceFloorPotentialEnergyGTest|surfaceIPCTopology_gtest|surfaceIPCSelfBroadPhase_gtest)' \
  --output-on-failure

git diff --check
```

Focused IPC CTest runs also passed: 27/27 for core/embedded IPC/floor validation and 39/39 for the wider IPC collection.

This is not a claim that the full repository, Python package, or C API test suite was rerun in this final pass. Run the full project matrix before a release or after toolchain/dependency changes.

## Files most likely to matter next

### IPC implementation

- `src/core/contact/ipc/core/surfaceIPCCore.h`
- `src/core/contact/ipc/core/surfaceIPCCore.cpp`
- `src/core/contact/ipc/core/surfaceIPCBarrierAssembler.h`
- `src/core/contact/ipc/core/surfaceIPCBarrierAssembler.cpp`
- `src/core/contact/embeddedSurfaceIPCPotentialEnergy.cpp`
- `src/core/contact/mappedSurfacePotentialEnergy.cpp`

### Tests

- `tests/src/core/contact/surfaceIPCCore_gtest.cpp`
- `tests/src/core/contact/surfaceIPCBarrierAssembler_gtest.cpp`
- `tests/src/core/contact/ipcGeometry_gtest.cpp`
- `tests/src/core/contact/surfaceIPCMaxStep_gtest.cpp`
- `tests/src/core/contact/embeddedSurfaceIPCPotentialEnergy_gtest.cpp`
- `tests/src/core/solidDeformationModel/tetMeshDeformationModel_gtest.cpp`
- `tests/src/core/solidDeformationModel/deformationModelAssembler_gtest.cpp`

## Recommended next-session entry point

1. Start from the clean `727a037` worktree and read this file plus `RELEASE_0.0.4_PHASE0_SPIKE_REPORT.md`.
2. If continuing IPC work, preserve the two Hessian modes: projected for robust solves, raw for mathematical FD checks.
3. Before changing dynamic topology behavior, run the Newton and IPC aggregate tests; the distinction between pure-fixed and mixed fixed+dynamic is intentional.
4. Do not restart P1.6/P1.7 without an explicit decision to broaden scope.
5. For release confidence, run the full configured test suite and any Python/C API checks used by the release pipeline.

The documentation submodule/worktree was not modified as part of Phase 1.
