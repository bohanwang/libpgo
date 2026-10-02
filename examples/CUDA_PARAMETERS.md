# CUDA Newton and IPC parameter reference

Source-audited on 2026-10-02. This describes the implemented libpgo CUDA route, not the upstream libuipc scene API. Defaults below are compiled defaults; values chosen in individual examples or benchmark configurations are not new defaults. All scalar numerical inputs described as finite exclude NaN and infinity.

The tables distinguish **model** inputs (geometry, forces, material/contact regularization), **discretization** inputs (time integration and friction lagging), **strategy** controls (how the same objective is solved), **validity** controls (which configurations are admissible), and **I/O** controls. Strategy changes can still select different admissible solutions or cause different adaptive retries; they are not promises of identical trajectories.

## JSON: selecting and configuring a run

Build with `PGO_ENABLE_CUDA=ON`; `CMAKE_CUDA_ARCHITECTURES` selects the target GPU architecture at build time. CUDA is optional and defaults off. The general C/Python simulation dispatcher reads `contact-model` (select `"ipc"`), whereas the `runIPCSim` executable already selects IPC. IPC defaults to the CPU backend unless `simulation-backend` is `"cuda"`. The legacy `runSim` executable directly invokes sampled CPU simulation; changing these keys does not turn that executable into the CUDA IPC entry point.

Sources: [backend dispatch](../src/simulationRunner/simulationRunner.cpp#L18), [IPC backend selection](../src/simulationRunner/ipcSimulationRunner.cpp#L74), [CUDA build option](../CMakeLists.txt#L20), [CUDA input validation](../src/simulationRunner/cudaIpcSimulationRunner.cpp#L172), [mesh input selection](../src/simulationRunner/runSimVolumeMeshIO.cpp#L32).

| Key | Default / allowed values | Effect |
| --- | --- | --- |
| `contact-model` | General C/Python dispatcher default `"sampled"`; use `"ipc"` | Route selection. `runIPCSim` already selects IPC; legacy `runSim` stays sampled CPU. CUDA is unavailable for sampled contact. |
| `simulation-backend` | `"cpu"`; `"cpu"` or `"cuda"` | Backend selection. |
| `cuda-device` | `0`; available logical CUDA device index | I/O/resource selection. Indices are after `CUDA_VISIBLE_DEVICES` filtering. |
| `sim-type` | **Required**; `"static"` or `"dynamic"` | Model/solve mode. Static solves the full-load equilibrium; dynamic uses backward Euler. |
| `elastic-material` | **Required** | Model. `"stable-neo"` or `"stvk-vol"` for volumes; `"koiter-stvk"` for shells. |
| `tet-mesh` | Required for tetrahedral volume; mutually exclusive with `cubic-mesh` | Model. Volume mesh path; four-node tets. |
| `cubic-mesh` | Required for cubic volume; mutually exclusive with `tet-mesh` | Model. Volume mesh path; eight-node trilinear hexes. |
| `surface-mesh` | **Required** | Model. Triangular collision surface, or the shell itself if neither volume-mesh key is present. Volume displacement embedding is computed from the meshes; there is no JSON embedding-matrix option. |
| `scale` | **Required**; finite `> 0` | Model. Multiplies deformable rest-mesh coordinates. Does not automatically scale every other dimensional parameter. |
| `init-disp` | **Required**; three finite numbers | Model. Uniform initial displacement, already in scaled simulation length units. |
| `init-vel` | **Required for dynamic**; three finite numbers | Model. Uniform initial velocity in length/time. Static uses zero. |
| `g` | **Required**; three finite numbers | Model. Acceleration; external force is `M*g`, including static runs. |
| `fixed-vertices` | **Required**; may be `[]` | Model. Soft quadratic attachment entries, detailed below; these do not remove degrees of freedom. |
| `ipc-external-objects` | Absent means none | Model. Prescribed triangular obstacles, detailed below; no rigid-body unknowns. |
| `output` | **Required**; path | I/O. Created if needed. Relative paths resolve against the configuration directory. |
| `restart-from-u` | `false`; boolean; dynamic only | I/O. Resume from the latest eligible `deform####.u` in `output`, including displacement, velocity and acceleration. |
| `dump-interval` | **Required for dynamic**, integer `>= 1`; static uses `1` | I/O. Macroframe OBJ cadence. Every committed macro state and every accepted physical-substep state/surface are still exported. |
| `loglevel` | `"info"`; `"trace"`, `"debug"`, `"info"`, `"warn"`, `"error"` | I/O. Logger level, separate from CUDA Newton verbosity. |
| `cuda-solver-verbose` | `0`; integer `0..3` | I/O. Increasing solver diagnostics; level 3 includes trial/component details and rejected-attempt state dumps. Does not change the objective or convergence criterion. |

Mesh material records supply volumetric Young's modulus, Poisson ratio and density. There are no CUDA JSON overrides for these. The CUDA FEM copies the first `SimulationMeshENuMaterial` per element and requires finite `E > 0`, `-1 < nu < 0.5`. Shell material parameters are separately exposed below. Units are user-selected but must be consistent: length, time, mass, modulus and density are not automatically converted to SI.

Sources: [initial displacement and scaling](../src/simulationRunner/runIPCSimSetup.cpp#L59), [volume material initialization](../src/simulationRunner/runIPCSimSetup.cpp#L521), [CUDA element material checks](../src/core/solidDeformationModel/cuda/cudaDeformationModelEnergy.cu#L728), [log-level parser](../src/simulationRunner/runSimCliLogging.cpp#L81), [restart](../src/simulationRunner/cudaIpcSimulationRunner.cpp#L331).

## JSON: Newton, PCG and adaptive stepping

Sources: [runner settings](../src/simulationRunner/cudaIpcSimulationRunner.cpp#L196), [CUDA strategy settings](../src/simulationRunner/cudaIpcSimulationRunner.cpp#L286), [compiled solver defaults](../src/core/nonlinearOptimization/cuda/cudaNewtonSolver.h#L25), [solver range validation](../src/core/nonlinearOptimization/cuda/cudaNewtonSolver.cu#L181).

| Key | Default / allowed values | Effect |
| --- | --- | --- |
| `solver-max-iter` | **Required**; integer `>= 1` | Strategy. Newton iteration cap for **each** solve, including each physical substep and each friction lagging iteration. |
| `solver-eps` | **Required**; finite `> 0` | Accuracy. Convergence requires the complete, unscaled objective gradient infinity norm to be **strictly less** than this value. Absolute, not relative, not a displacement or energy-change tolerance. |
| `cuda-pcg-max-iter` | `0`; integer `>= 0` | Strategy. `0` means `min(10000, 2*nDOF)` PCG iterations. Positive values set the cap per linear-system attempt, including regularized attempts. |
| `cuda-pcg-tolerance` | `1e-6`; finite `(0,1)` | Strategy. Relative **L2 linear residual** tolerance. Different quantity from `solver-eps`; formula below. |
| `cuda-pcg-batch-size` | `0`; integer `0..64` | Strategy/performance. `0` reads scalar checks on the host each iteration. Positive values perform the same checks on the GPU and read sticky completion/failure status once per batch. No convergence checks are skipped. |
| `cuda-exact-hessian-start-iteration` | `20`; integer | Strategy. Starting at this accepted Newton iteration, try the unprojected tangent first. Negative disables; `0` tries from the beginning. Invalid curvature, nonfinite results or non-descent cause projected fallback. Resets for each solve. |
| `cuda-stagnation-window` | `0`; integer `0..64` | Strategy. `0` disables early stagnation detection. Positive values bound a streak of input-chord-clipped, representationally tiny steps with insufficient measured residual progress. Relevant to dynamic solves, which enable the input-chord check. |
| `cuda-stagnation-tolerance` | `64*double_epsilon`, approximately `1.4210854715202004e-14`; finite `(0,1)` | Strategy. Relative-change threshold for both anchored gradient L2 and infinity norms in that window. May trigger residual recovery or a **not-converged** return so the runner retries; never relaxes `solver-eps`. Ineffective when the window is zero. |
| `timestep` | **Required for dynamic or friction-enabled static**; finite `> 0`; frictionless static internally uses `1` | Discretization. Macrostep duration `h`; also sets friction regularization displacement `ipc-friction-epsv*h`. Static friction uses it as a regularization timescale, not a simulated equilibrium duration. |
| `num-timestep` | **Required for dynamic**; integer `>= 0`; static uses `1` | Discretization. Number of macrosteps; also determines prescribed movement schedules. Zero is accepted but performs no dynamic steps. |
| `damping-params` | **Required for dynamic**; `[alpha,beta]`, both finite `>= 0`; static uses `[0,0]` | Model. Rayleigh mass and frozen-stiffness damping. `alpha` has inverse-time units; `beta` has time units. |
| `cuda-max-timestep-refinements` | `8`; integer `0..12` | Discretization/retry strategy. On a failed macrostep, restore its input and retry using `2^r` equal physical substeps, `r=0..limit`. Default allows 256 substeps; maximum allows 4096. `0` permits only the original step. Static does not subdivide. |

For a linear solve with right-hand side `-g`, PCG uses

```text
target_L2 = max(min(pcgAbsoluteTolerance, 0.1*solver-eps),
                cuda-pcg-tolerance * ||g||_2)
```

`pcgAbsoluteTolerance` is a C++-only field, default `1e-12`. The actual residual `H*dx + g` is recomputed before reporting linear convergence. An iteration-limited PCG direction can still be usable if finite and descending; it does not imply Newton convergence. The only successful Newton status is the strict final gradient test. There is no configurable “accept after cap,” relative Newton tolerance, or KKT convergence substitute. [PCG implementation](../src/core/nonlinearOptimization/cuda/cudaNewtonSolver.cu#L430), [Newton acceptance](../src/core/nonlinearOptimization/cuda/cudaNewtonSolver.cu#L552).

The dynamic quadratic Hessian is `M*(1/h_sub^2 + alpha/h_sub)`, plus attachment terms. Stiffness damping adds `beta/h_sub` times the elastic Hessian frozen at that physical substep's input. Elastic and contact energies are not multiplied by `h_sub^2`. Smaller timesteps therefore change the discretized objective and its raw-gradient scale. [Incremental objective](../src/simulationRunner/cudaIpcSimulationRunner.cpp#L485).

Dynamic acceptance requires convergence and a certificate for each actual input-to-output physical chord. Static acceptance records and certifies successive Newton-iterate chords; it does not require visibility from the initial configuration to every iterate. The coarse shortcut between macroframe endpoints is reported separately when recorded physical substeps exist. Failed refinements are rolled back completely. A converged intermediate solve in a subsequently rejected attempt is not a committed state. [Path controls](../src/simulationRunner/cudaIpcSimulationRunner.cpp#L286), [adaptive attempts](../src/simulationRunner/cudaIpcSimulationRunner.cpp#L452).

## JSON: contact and friction

Sources: [distance configuration](../src/simulationRunner/runIPCSimSetup.cpp#L133), [shell parameters](../src/simulationRunner/runIPCSimSetup.cpp#L160), [volume parameters](../src/simulationRunner/runIPCSimSetup.cpp#L185), [parameter validation](../src/core/contact/ipc/core/surfaceIPCCore.cpp#L109), [friction iterations](../src/simulationRunner/cudaIpcSimulationRunner.cpp#L184).

| Key | Default / allowed values | Effect |
| --- | --- | --- |
| `ipc-dhat` | **Required choice**, finite `> 0`, unless shell heuristic is enabled | Model/contact regularization. Absolute activation **distance**, in scaled simulation length units; not squared distance or collision thickness. Mutually exclusive with `ipc-dhat-relative`. |
| `ipc-dhat-relative` | Alternative to `ipc-dhat`; finite `> 0` | Model/contact regularization. Multiplied once by the scaled deformable rest surface's bounding-box diagonal, excluding external objects. |
| `ipc-kappa` | **Required**, finite `> 0`, unless shell heuristic is enabled | Model/contact regularization. Constant multiplier of this repository's weighted, normalized barrier; not automatically adapted during Newton. |
| `ipc-heuristic` | `false`; boolean; `true` supported only for shells | Model/contact setup. For shells, sets `dhat = 0.001*bbox_diagonal` and `kappa = shell-youngs-modulus*shell-thickness`. Explicit distance/stiffness keys are ignored in this branch. `true` on a volume is rejected. |
| `contact-friction-coeff` | `0`; finite `>= 0` | Model. Dimensionless Coulomb coefficient `mu`. Zero disables the separate friction energy. |
| `ipc-friction-epsv` | `1e-3`; finite `> 0` | Model/friction regularization. Tangential speed in simulation length/time. Smoothed slip threshold is `epsilon = epsv*h_sub`; validated even if friction is disabled. |
| `ipc-friction-iterations` | `1`; integer `>= 1` | Discretization/lagging. Fixed number of lagged-friction Newton solves per physical step or static solve. Validated even at `mu=0`, but then only one Newton solve is performed. There is no outer residual tolerance or automatic outer iteration cap. |

For squared primitive distance `s`, `r = s/dhat^2`, the active barrier is

```text
b(s) = -(r - 1)^2 * log(r),        0 < r < 1
E_pair = kappa * geometric_weight * mollifier * b(s)
```

It is inactive at `s >= dhat^2`. Nonpositive contact distances are outside the valid solve domain. Point–triangle weight is `vertexArea*triangleArea`, with `vertexArea` the sum of adjacent rest triangle areas divided by three. Edge–edge weight is the product of the two rest edge lengths. Consequently, the weights have different length powers; `kappa` is this implementation's quadrature multiplier, not a universal Young's-modulus quantity. Do not copy upstream libuipc stiffness/timestep/thickness settings assuming identical scaling. There is no extra `dhat^4` or timestep-squared multiplier here. [Barrier definition](../src/core/contact/ipc/geometry/ipcBarrier.cpp#L17), [GPU pair weights](../src/core/contact/ipc/cuda/cudaEmbeddedSurfaceIPCPotentialEnergy.cu#L150), [rest measures](../src/core/contact/ipc/topology/surfaceIPCTopology.cpp#L96).

Friction freezes closest-point weights, tangent projector and normal force at the lagged state, while retaining the **original physical-step reference** across outer iterations. Only the last outer solve advances velocity/time. Prescribed obstacle vertices use their current placement in both reference and lagged states, matching the CPU convention; this is not an independently integrated rigid friction model. [Friction interface](../src/core/contact/ipc/cuda/cudaSurfaceIPCFrictionPotentialEnergy.h#L15), [runner lifecycle](../src/simulationRunner/cudaIpcSimulationRunner.cpp#L505).

## JSON: material validity and shells

Sources: [material guard switch](../src/simulationRunner/runIPCSimSetup.cpp#L337), [cubic selection and restart metadata](../src/simulationRunner/cudaIpcSimulationRunner.cpp#L172), [CUDA material contract](../src/core/solidDeformationModel/cuda/cudaDeformationModelEnergy.h#L12), [shell material parsing](../src/simulationRunner/runIPCSimSetup.cpp#L101).

| Key | Default / allowed values | Effect |
| --- | --- | --- |
| `enable-material-max-step` | `true`; boolean | Validity. Enables volumetric determinant step bounds, including initial zero-direction feasibility. `false` disables both Gauss/tet and additional whole-cell **motion** guards; it does not change the energy. Whole-cell mode still requires certified rest cells at construction. No area max-step guard is implemented by CUDA Koiter, so this flag does not create one for shells. |
| `cuda-cubic-validity` | `"gauss"`; `"gauss"` or `"whole-cell"`; explicit key requires `cubic-mesh` | Validity. Gauss mode follows the CPU's eight quadrature locations. Whole-cell additionally certifies the entire trilinear cell along the chord, failing closed on uncertainty. Neither mode changes material energy, quadrature, derivatives or `solver-eps`. |
| `shell-thickness` | `0.003`; finite `> 0` | Model. Koiter membrane/bending thickness. **Not** IPC collision thickness or shell mass density. |
| `shell-areal-density` | `100`; finite `> 0` | Model. Shell mass per unit rest area; independent of thickness. |
| `shell-youngs-modulus` | `1e6`; finite `> 0` | Model. Membrane Young's modulus. |
| `shell-bending-youngs-modulus` | Membrane Young's modulus; finite `> 0` | Model. Bending modulus. |
| `shell-poisson-ratio` | `0.4`; finite `(-1,0.5)` | Model. Membrane Poisson ratio. |
| `shell-bending-poisson-ratio` | Membrane Poisson ratio; finite `(-1,0.5)` | Model. Bending Poisson ratio. |
| `shell-zero-rest-curvature` | `false`; boolean | Model. Use zero rest curvature instead of the rest mesh's curvature. |
| `initial-surface-mesh` | Absent; shell only | Model/initial condition. Same vertex count and triangle topology as the shell rest mesh; initializes its deformed shape. Scaled initial-minus-rest positions are added to `init-disp`. Rejected for volume meshes. |

The fixed enabled volumetric floor is **`det(F) > 1e-8`**, not a configurable JSON value. For a tet, `F` is constant inside the element. For a cubic cell, the eight Gauss positions are all combinations of `0.5 +/- 0.5/sqrt(3)` in the unit parameter cube. Gauss positivity does **not** certify corners or the interior between samples: a converged, surface-collision-free configuration may have whole-cell inversions. Independent reports retain those inversions as separate diagnostics in Gauss mode. Whole-cell mode imposes a stronger feasible domain without adding a restoring force, so it can prevent reaching a CPU Gauss-valid stationary point. [Gauss locations](../src/core/solidDeformationModel/cuda/cudaDeformationModelEnergy.cu#L705), [determinant constants](../src/core/solidDeformationModel/cuda/cudaDeformationModelEnergy.cu#L27).

`material-validity.json` records the selected criterion, guard flag and floor. Restart rejects a change to recorded validity settings. A cubic restart without legacy validity metadata requires explicit `cuda-cubic-validity`; it does not retroactively validate earlier states. Shell-specific material keys are consumed only for shells. [Restart checks](../src/simulationRunner/cudaIpcSimulationRunner.cpp#L246), [shell initial state](../src/simulationRunner/runIPCSimSetup.cpp#L452).

## JSON: prescribed obstacles and attachments

`ipc-external-objects` is an array of these objects. All vector components must be finite. Meshes must be valid triangular surfaces. Sources: [object parser](../src/simulationRunner/runIPCExternalObjects.cpp#L87), [defaults and movement convention](../src/simulationRunner/runIPCExternalObjects.h#L17).

| Object field | Default / range | Meaning |
| --- | --- | --- |
| `filename` | Required path | Obstacle mesh; relative to simulation configuration directory. |
| `scale` | `1`; finite `> 0` | Separate obstacle scale; global deformable `scale` is not applied. |
| `initial-translation` | `[0,0,0]`; finite vector | Added after obstacle scaling. |
| `movement` | `[0,0,0]`; finite vector | Total prescribed translation in simulation units. At macroframe `k` of `N>1`, placement is base plus `movement*k/(N-1)`. A one-frame run stays at base. Static ignores movement. |

Obstacles are fixed within a macroframe, including its physical refinements and Newton solves, and moved between macroframes. The independent trajectory checker checks this separate prescribed-motion segment. There is no JSON obstacle rotation, dynamic rigid mass or contact group mask.

Each `fixed-vertices` entry requires finite nonnegative `coeff` and exactly one of the following forms. `coeff` scales a soft positional penalty; zero disables its force. Source: [attachment parser](../src/simulationRunner/runIPCSimSetup.cpp#L234).

| Form | Fields | Meaning |
| --- | --- | --- |
| Vertex-list movement | `filename`, `movement: [x,y,z]`, `coeff` | File of simulation vertex indices; target is scaled rest position plus the supplied movement. The movement itself is not multiplied by global scale. Dynamic interpolates the target over the configured frame schedule; static uses the full target. Supply valid indices and finite movement: this legacy branch does not explicitly validate them as thoroughly as trajectory tracks. |
| Keyed trajectory | `trajectory`, `coeff` | Path to the trajectory JSON below. Cannot also contain `filename` or `movement`. Dynamic samples the prescribed physical-step target time; static samples the last key time and keeps the attachment active. |

Trajectory file fields ([parser](../src/simulationRunner/runSimVertexTrajectory.cpp#L34), [file-level contract](../src/simulationRunner/runSimVertexTrajectory.cpp#L150), [target sampling](../src/simulationRunner/cudaIpcSimulationRunner.cpp#L463)):

| Field | Default / range | Meaning |
| --- | --- | --- |
| `format`, `version` | Required `"pgo-vertex-trajectory"`, `1` | File format. |
| `space` | `"position"`; `"position"` or `"displacement"` | Key coordinates before global scale; displacement keys are offsets from rest. |
| `interpolation` | `"linear"`; only supported value | Piecewise linear; clamp to first/last key outside the key interval. |
| `tracks` | Required nonempty array | Vertices cannot repeat across tracks within the same file. |
| `tracks[].name` | Generated `tracks[i]` | Diagnostic label. |
| `tracks[].vertices` | Required nonempty list of unique, in-range, nonnegative integers | Simulation vertex indices. |
| `tracks[].keys` | Required nonempty array | Every entry has finite time `t`, strictly increasing, and `p`, one finite XYZ vector per listed vertex. Times need not start at zero. |
| `tracks[].release-time` | Absent means never; finite if present | Dynamic attachment active only while sampled time is less than release time. Static does not release. |

## C++ controls not all exposed in JSON

### `CudaNewtonSolver::SolverParam`

The type is named `SolverParam`, not `Options`. All 14 fields are listed here, including the JSON mappings. Set these before constructing the solver. Sources: [declaration/defaults](../src/core/nonlinearOptimization/cuda/cudaNewtonSolver.h#L25), [validation](../src/core/nonlinearOptimization/cuda/cudaNewtonSolver.cu#L181).

| C++ field | Default | Range and JSON mapping / effect |
| --- | --- | --- |
| `maxPCGIterations` | `0` | Integer `>=0`; `cuda-pcg-max-iter`. |
| `pcgReductionBatchSize` | `0` | Integer `0..64`; `cuda-pcg-batch-size`. |
| `stagnationWindow` | `0` | Integer `0..64`; `cuda-stagnation-window`. |
| `stagnationRelativeTolerance` | `64*double_epsilon` | Finite `(0,1)`; `cuda-stagnation-tolerance`. |
| `pcgRelativeTolerance` | `1e-6` | Finite `(0,1)`; `cuda-pcg-tolerance`. |
| `pcgAbsoluteTolerance` | `1e-12` | **C++ only**, finite `>=0`; effective tolerance is capped by `0.1*NewtonTolerance` as above. |
| `maxLineSearchIterations` | `60` | **C++ only**, integer `>=1`; maximum backtracking trials per direction. |
| `armijoCoeff` | `1e-4` | **C++ only**, finite `(0,1)`; sufficient-decrease coefficient. |
| `maxRegularizationAttempts` | `12` | **C++ only**, integer `>=0`; shifted retries after the unshifted projected attempt. Zero disables shifted retries. |
| `initialRegularization` | `1e-8` | **C++ only**, finite `>0`; initial diagonal shift is this times `max(1e-12,maxAbs(diag(H)))`; subsequent shifts multiply by ten. |
| `exactHessianStartIteration` | `20` | Integer; `cuda-exact-hessian-start-iteration`. Negative disables. |
| `maxResidualCorrectionCandidates` | `32` | **C++ only**, integer `>=0`; bounded coordinate residual-recovery candidates, limited by DOF count. Zero disables this recovery. Energy, feasibility and strict full-residual safeguards remain active. |
| `checkInputChord` | `false` | **C++ only** boolean. Runner sets true for dynamic, false for static. Certifies each trial's whole chord from the solve input. |
| `checkTrialChord` | `false` | **C++ only** boolean. Runner sets true for static, false for dynamic. Explicitly certifies current-to-represented-trial endpoints after rounding. |

`solve(hostX, maxIterations, gradientTolerance, verbose=0)` supplies the Newton cap/tolerance/verbosity separately. The C++ cap allows zero (JSON requires at least one). The state must be finite and feasible even if its initial gradient already meets tolerance. `setAcceptedIterateObserver(...)` is an output callback for input iteration zero and all accepted iterates, including recovery; it is not a solver tolerance. Status values are `0` converged, `1` not converged, `2` numerical failure; failure retains the last finite accepted state. [Public API](../src/core/nonlinearOptimization/cuda/cudaNewtonSolver.h#L101).

### `SurfaceIPCCore::Parameters` used by CUDA contact

Sources: [parameter declaration](../src/core/contact/ipc/core/surfaceIPCCore.h#L39), [validation](../src/core/contact/ipc/core/surfaceIPCCore.cpp#L109), [GPU constructors](../src/core/contact/ipc/cuda/cudaEmbeddedSurfaceIPCPotentialEnergy.h#L22).

| C++ field | C++ default | Range and runner behavior |
| --- | --- | --- |
| `dhat` | `0.1` | Finite `>0`. JSON setup overrides it using required distance configuration or shell heuristic. |
| `kappa` | `0.1` | Finite `>0`. JSON setup overrides it using required stiffness or shell heuristic. |
| `eps_ee` | `0` | Finite `>=0`, **C++ only**. Edge–edge mollifier threshold for squared cross-product norm, with length-to-the-fourth units. Zero disables mollification. Passed directly, not converted from a relative rest-angle coefficient. JSON setup fixes zero. Changing this changes the contact law. |
| `slackness` | `1` | Finite `(0,1]`, **C++ only**. Multiplies partial CCD bounds to shorten a collision-limited Newton step. A full certificate stays `1`. JSON setup fixes one. |
| `projectHessianToPSD` | `true` | Boolean, **C++ only**. Project local barrier Hessians in normal `prepareHessian`; `prepareExactHessian` explicitly bypasses this projection. Changes the search tangent, not barrier energy/gradient. |
| `frictionCoeff` | `0` | Finite `>=0`; maps to JSON `contact-friction-coeff`. The barrier-only CUDA class rejects nonzero friction; the runner passes zero to it and creates a separate `CudaSurfaceIPCFrictionPotentialEnergy` with the original coefficient. |
| `frictionEpsV` | `1e-3` | Finite `>0`; maps to `ipc-friction-epsv`. |

The C++ constructor additionally accepts surface rest vertices, triangles, sparse displacement map `W`, and optionally a deformable-vertex mask. These define the geometry/embedding rather than tuning constants. Prescribed vertices have no simulation DOFs; fixed–fixed and incident primitive pairs are excluded. Self contact is included. Contact gradients and Hessian products are pulled back through `W`. The friction class takes the same geometry/parameters, and its `updateFriction(reference, lagged, timestep)` requires the caller to preserve the physical reference across lag iterations.

### C++ material and quadratic energy inputs

`CudaDeformationModelEnergy(mesh, enableMaterialMaxStep=true, CubicValidity::Gauss)` defaults to Stable Neo-Hookean. The overload accepts `DeformationModelElasticMaterial::STABLE_NEO` or `STVK_VOL`. `CubicValidity::WholeCell` selects the stronger guard. Mesh/material data are copied at construction; modifying the original mesh afterward is not a live parameter update. [Constructor](../src/core/solidDeformationModel/cuda/cudaDeformationModelEnergy.cu#L813).

The per-element `SimulationMeshENuMaterial` has `E`, `nu`, and a StVK volume parameter named `J` (`setE`, `setNu`, `setJ`). Its bare defaults are `E=6000`, `nu=0.4`, `J=10000`; loaded volume meshes replace `E` and `nu`. The standard loader leaves `J=10000`; there is no JSON compression-ratio setting. For StVK-vol, this parameter has **two roles**: it directly supplies the coefficient of `0.5*J*(det(F)-1)^2`, and, when positive, sets the separate compression response coefficient `J*E/(1-2*nu)`. Nonpositive `J` disables only the latter response; a negative value still supplies a negative quadratic volume coefficient. CUDA checks finiteness, not positivity, of this parameter and its derived coefficients. Thus it should not be interpreted solely as a dimensionless multiplier. The material parameter `J` is distinct from the deformation determinant `det(F)`. Stable Neo-Hookean ignores it. [Material data](../src/core/solidDeformationModel/simulationMesh.h#L34), [volume energy](../src/core/solidDeformationModel/cuda/cudaDeformationModelEnergy.cu#L171), [CUDA conversion](../src/core/solidDeformationModel/cuda/cudaDeformationModelEnergy.cu#L770).

For StVK-vol, the normal Hessian path projects the invariant-material tangent and then adds the unprojected volume tangent, matching the CPU's combined-material convention. The resulting tangent is not guaranteed positive semidefinite; PCG curvature checks and regularization still matter. The exact path keeps both tangents unprojected. [Tangent assembly](../src/core/solidDeformationModel/cuda/cudaDeformationModelEnergy.cu#L438).

`CudaKoiterDeformationModelEnergy(mesh, parameters, zeroRestCurvature=false)` takes five finite coefficients per patch: membrane `E,nu`, bending `E,nu`, thickness. Runner-level physical ranges are in the shell table. `CudaQuadraticEnergy` accepts its sparse matrix, reference, linear term and constant through C++; `CudaFrozenHessianEnergy::setReference` takes a reference and coefficient. The runner derives these from inertia, damping and attachments; there are no separate JSON Hessian, force-offset or energy-offset knobs. The volume mass matrix is consistent and the shell mass matrix is lumped; neither is selected by a CUDA JSON mass-lumping switch.

## Fixed implementation choices: not tunable parameters

These values explain behavior near difficult cases. Changing them requires implementation work; adding similarly named JSON fields has no effect.

| Area | Current fixed behavior / limits | Source |
| --- | --- | --- |
| Linear algebra | Double precision; Hessian-vector operators with Jacobi preconditioning; no selectable direct solver, multigrid/MAS, mixed precision or preconditioner switch. | [PCG](../src/core/nonlinearOptimization/cuda/cudaNewtonSolver.cu#L430) |
| Jacobi diagonal | Floor `max(1e-30, diagonalScale*1e-12)`; diagonal scale at least `1e-12`. | [preconditioning](../src/core/nonlinearOptimization/cuda/cudaNewtonSolver.cu#L443) |
| Backtracking | Step halving by `0.5`; component-wise stable energy differences. Roundoff allowance `16*double_epsilon*roundoffScale`; roundoff-limited acceptance requires gradient infinity norm below tolerance or reduced to at most `0.9` of the old value. Recovery additionally requires strict full L2 residual improvement. | [line search](../src/core/nonlinearOptimization/cuda/cudaNewtonSolver.cu#L674) |
| Stagnation step scale | Represented step/position scale at most `16*double_epsilon`; user threshold controls residual changes only. | [stagnation](../src/core/nonlinearOptimization/cuda/cudaNewtonSolver.cu#L785) |
| Broadphase | GPU Morton LBVH, outward-rounded AABBs, automatic checked candidate-buffer growth. No JSON BVH/candidate-capacity/contact-count truncation setting. | [broadphase](../src/core/contact/ipc/cuda/cudaSurfaceIPCBroadPhase.cu) |
| Contact CCD | Zero-thickness surfaces; original-chord conservative advancement, at most 1000 iterations, then interval fallback with depth 32/work 4096. Uncertainty returns a certified prefix, never a full pass. | [CCD](../src/core/contact/ipc/cuda/cudaIPCContinuousCollision.h#L193) |
| Collision-limited prefix | Retains a conservative 10% initial-gap reserve; advancement safety factor `0.9`. Full clear chords still return `1`. This is separate from C++ `slackness`. | [prefix reserve](../src/core/contact/ipc/cuda/cudaIPCContinuousCollision.h#L235) |
| CCD arithmetic uncertainty | Directed bounds include mapping/pose roundoff; no configurable absolute CCD epsilon, contact thickness or tolerance bypass. | [CCD evaluation](../src/core/contact/ipc/cuda/cudaIPCContinuousCollision.h) |
| Volumetric Gauss/tet bound | Determinant floor `1e-8`; boundary tolerance `1e-12`; polynomial root bisection up to 64 iterations; partial-step safety `0.99`. | [material step bound](../src/core/solidDeformationModel/cuda/cudaDeformationModelEnergy.cu#L512) |
| Whole-cell bound | Directed interval proof with depth 96/work 4096; 16 prefix-bisection iterations and `0.99` safety. Unresolved boxes fail closed. | [hex proof](../src/core/solidDeformationModel/cuda/cudaHexJacobian.h#L185), [prefix search](../src/core/solidDeformationModel/cuda/cudaDeformationModelEnergy.cu#L645) |
| Cubic integration | Eight Gauss samples; per-sample material tangent projection. No quadrature-order option. | [quadrature](../src/core/solidDeformationModel/cuda/cudaDeformationModelEnergy.cu#L705) |
| Material protection for shells | CUDA Koiter `computeMaxStepSize` returns `1`; no runtime triangle-area guard. Independent trajectory validation checks shell geometry separately. | [Koiter interface](../src/core/solidDeformationModel/cuda/cudaKoiterDeformationModelEnergy.h) |
| CUDA launch/storage choices | Kernel block sizes, allocation growth, BVH sorting, reduction ordering and local eigensolver implementation are internal; no JSON launch-block or assembly-order controls. | [CUDA modules](../src/core/contact/ipc/cuda/) |

## Rejected settings, ineffective settings and aliases

There are no alternate spellings/aliases for the CUDA solver keys above. `ipc-dhat-relative` is an alternative distance representation, not an additional distance multiplier. JSON `cuda-pcg-tolerance` maps to C++ `pcgRelativeTolerance`; it is not the Newton tolerance.

Explicitly rejected in this route:

- `deterministic-assembly` set to `true` (default false; GPU assembly is not deterministic).
- Any `linear-solver` field, even one selecting a CPU solver that is otherwise available.
- Legacy `external-objects`; use `ipc-external-objects`.
- Any `use-floor`, `floor-axis`, `floor-height` or `floor-kappa` field; use an external triangle mesh.
- Per-obstacle `init-disp`; use `initial-translation`.
- Volumetric `ipc-heuristic: true`, volumetric `initial-surface-mesh`, `cuda-cubic-validity` without a cubic mesh, or static `restart-from-u: true`.
- Unsupported material names and sampled contact with `simulation-backend: "cuda"`.

The parser does **not** implement a catch-all unknown-key rejection. CPU-only `pardiso-deterministic` and `pardiso-thread-limit` are not consumed by the CUDA route. Invented keys such as `ipc-eps-ee`, `ipc-slackness`, `cuda-armijo`, `cuda-pcg-absolute-tolerance`, a JSON preconditioner selector or a determinant-floor override do not configure their similarly named C++ quantities. An unrecognized key can therefore be silently ineffective. Parameters for inactive branches can also be ignored: static ignores dynamic damping/initial velocity/count settings; volume ignores shell material keys; shell heuristic ignores explicit `ipc-dhat`, `ipc-dhat-relative` and `ipc-kappa`.

Sources: [CUDA rejections](../src/simulationRunner/cudaIpcSimulationRunner.cpp#L172), [removed floor options](../src/simulationRunner/runIPCSimSetup.cpp#L342), [external-object checks](../src/simulationRunner/runIPCExternalObjects.cpp#L87), [CPU-only solver parser](../src/simulationRunner/runSimCliLogging.cpp#L103), [JSON access semantics](../src/core/configFileJSON/configFileJSON.h#L68).

## Invocation and measurement controls

The `runIPCSim` executable also accepts these flags; they are **not JSON keys or solver strategy parameters**. Sources: [CLI](../src/tools/runSim/runIPCSim.cpp#L13), [C++ invocation options](../src/simulationRunner/simulationRunner.h#L7), [timing schema](../src/simulationRunner/runSimTiming.h#L12).

| CLI flag | Default / range | Effect |
| --- | --- | --- |
| `--log` | Off; flag | Save CLI output to `output/runIPCSim.log` through the runner's log redirect. |
| `--max-frames N` | `-1` for all remaining frames; otherwise integer `>0` | Limit this invocation, including a checkpoint restart, without changing configured `num-timestep`, timestep or prescribed movement/force schedules. This is the bounded-checkpoint benchmark control. |
| `--timing-csv PATH` | Empty, disabled | Write per-frame compute, explicit output and total wall times, plus terminal/accepted/all-solve status and solve count. Compute includes runner/logging and rejected attempts; it is not a GPU-kernel-only timer. Process startup/setup are outside the frame timer. |

The matching C++ type is `IPCSimulationRunOptions`, with `maxFrames=-1` and empty `timingOutput`, passed to `runIPCSimulationFromConfig`. These controls apply to both CPU and CUDA IPC. A bounded successful invocation is not evidence that all originally configured frames have run. For an isolated step, seed a separate output directory with the exact preceding displacement/velocity/acceleration checkpoint, retain the original configuration schedule, enable restart, and use `--max-frames 1`; do not shorten `num-timestep` to one because that changes scheduled inputs.

Validation/checker CLI options, video-renderer camera settings, CUDA toolkit environment variables and CPU thread counts are separate tools/runtime controls. They do not change the CUDA Newton/contact parameter set documented here; changing a checker criterion does not change the simulated material or retroactively certify a trajectory.
