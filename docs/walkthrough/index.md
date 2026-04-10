# Walkthrough

This section contains code-oriented walkthroughs for important runtime flows.

- [Build System Walkthrough](build_system/index.md) - a Chinese, contributor-oriented series on how Conan, CMake presets, scikit-build-core, root CMake logic, local recipes, and CI/CD fit together in the current repo.
- [Simulation Tick: Pre-Tick Setup](simulation-tick-pre-tick.md) - volumetric mesh loading, simulation mesh construction, surface embedding, interpolation matrix assembly, and deformation model setup before stepping the solver.
- [Cubic Mesh Walkthrough](cubic_mesh/index.md) - a repo-aligned series on how cubic/hexahedral meshes move from asset generation into `SimulationMesh`, element models, global energy assembly, and shared runtime examples.
- [IPC](ipc/index.md) - a repo-aligned series on how the IPC collision detection and energy evaluation is structured and implemented in the current repo.
