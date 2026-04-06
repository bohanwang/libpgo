# Pulled Cubic Box Self IPC

This is the formal self-contact `ipc-barrier` smoke example.

Visual result:

![IPC Pull Example](ipc-pull.gif)

![IPC vs Penalty Comparison](ipc_penalty_compare.gif)

- The volumetric mesh and collision surface stay consistent.
- The surface mesh is the local boundary mesh in [`cubic-box.obj`](cubic-box.obj).
- The cubic mesh is a softened copy of the same cube geometry in [`soft-cubic-box.veg`](soft-cubic-box.veg).
- Two pulling constraints keep the bottom face fixed and drive the top face downward, so self near-contact appears naturally from the simulated deformation instead of a hand-edited folded rest surface.
- The current tuning prioritizes a steadier visual result over maximal compression: the top target now stops slightly above the bottom face, the pull stiffness is lower, and the motion is spread over more solver steps before each dumped frame.

Directory contents:

- [`pulled-cubic-box-self-ipc.json`](pulled-cubic-box-self-ipc.json): simulation config
- [`soft-cubic-box.veg`](soft-cubic-box.veg): softened volumetric mesh used by the smoke
- [`cubic-box.obj`](cubic-box.obj): matching boundary mesh for simulation and playback
- [`bottom-face.txt`](bottom-face.txt) and [`top-face.txt`](top-face.txt): constrained vertex sets
- [`anim.json`](anim.json): viewer sequence config for the dumped OBJ frames
- [`ipc-pull.gif`](ipc-pull.gif): IPC-only result preview
- [`ipc_penalty_compare.gif`](ipc_penalty_compare.gif): IPC vs penalty comparison preview

What this example is for:

- Regression coverage for the self `ipc-barrier` runtime path.
- A checked-in scene that reaches self near-contact without the obvious late-frame over-compression of the earlier debug setup.
- A longer default run that takes 40 solver steps and still dumps 20 OBJ frames for playback through [`anim.json`](anim.json).

What the current IPC update added around this example:

- External and self `ipc-barrier` runtime branches on the dynamic path
- Self near-contact active-set construction
- Self barrier energy
- External / self / merged feasible-alpha runtime filtering
- API smoke coverage for external-only, self-only, and merged ext+self paths

What it is not:

- Not a polished physics showcase.
- Not a full primitive self-IPC benchmark or paper-equivalent implementation.
- Not a frictional IPC benchmark.

Current modeling scope:

- The self path is still sample-based (point-triangle samples), not full primitive PT/EE constraints.
- Runtime feasibility uses a linearized feasible-alpha upper bound callback.
- This scene should be read as an integration/validation example for the current repo-aligned IPC path.

The intentionally softened material only exists to make the smoke scene short and repeatable while keeping the surface and volume meshes aligned.

Runtime note:

- The checked-in example runs for 41 timesteps with `dump-interval = 2`, so it still dumps frames `ret0001.obj` through `ret0020.obj`.
- It uses slightly larger `external-ipc-dhat` / `self-ipc-dhat` and `external-ipc-kappa` / `self-ipc-kappa`, a softer top pull, and a higher solver iteration cap than the earlier smoke version so the barrier activates earlier while the compression arrives more gradually.
- The API smoke test overrides this back to a short run so CI does not pay the full example cost.

How to run:

```bash
cd libpgo
cmake --preset core-debug
cmake --build --preset core-debug --target runSim
./build/core-debug/bin/runSim examples/pulled-cubic-box-self-ipc/pulled-cubic-box-self-ipc.json
```

How to run the focused IPC regression tests:

```bash
cd libpgo/build/core-debug/bin
./core_scene_self_contact_handler_test
./api_runSim_config_parse_test --gtest_filter='RunSimConfigParseTest.RunSimFromConfigCubicDynamicIpcNearContactActivatesExternalBarrier:RunSimConfigParseTest.RunSimFromConfigCubicDynamicIpcSelfBarrierSmokeTest:RunSimConfigParseTest.RunSimFromConfigCubicDynamicIpcMergedBarrierSmokeTest:RunSimConfigParseTest.RunSimFromConfigCubicDynamicIpcMergedDeterministicSmokeTest'
```
