# libpgo 0.0.4 Phase 2.0 Static-contact Spike Report

## Result

Phase 2.0 establishes two implementation facts before the unified-runner
extraction begins.

1. Static IPC is viable with the current core contracts. A two-tet,
   self-contact test combines deformation, a downward load, dynamic IPC
   Hessian assembly, the static Newton solver, and CCD maximum-step handling.
   It reaches a finite lower-energy result while preserving fixed DOFs.
2. The existing sampled-contact energy has a fixed active set. A contact
   energy built before penetration contains no constraints and remains zero
   after the surface subsequently penetrates an external mesh. Re-running the
   handler and rebuilding the energy detects the penetration and yields a
   positive energy.

## Evidence

- `NewtonSolverGTest.StaticIpcSelfContactUsesCcdMaxStepAndProducesFiniteResult`
  in `tests/src/core/newtonSolver_gtest.cpp` uses two disconnected tet
  surfaces. The lower tet and the upper tet base are fixed; a load drives the
  upper apex toward the lower tet. The test verifies dynamic aggregate
  topology, a sub-unit CCD step, finite energy/gradient, energy decrease,
  fixed-DOF preservation, and calls to the IPC maximum-step path during
  Newton line search.
- `SampledStaticContactGTest.ContactEnergyBuiltBeforePenetrationMissesLaterContact`
  in `tests/src/core/contact/sampledStaticContact_gtest.cpp` first constructs
  a sampled external-contact energy with zero active samples, then moves the
  surface into a box. The old energy remains zero; a rebuilt energy is
  positive.

## Consequence for Phase 2

The IPC static backend can keep one persistent
`EmbeddedSurfaceIPCPotentialEnergy`: it prepares contact pairs from the
current state and reports a CCD maximum step.

The sampled static backend cannot use one contact energy for an entire Newton
solve. Its implementation must either:

- rebuild sampled contact energies after each accepted static step, while
  destroying the old aggregate and buffers before mutating the handlers; or
- receive an explicit product decision that narrows static sampled-contact
  support.

The spike does not establish protection against a Newton step tunneling
through an initially separated sampled-contact obstacle. That requires a
separate product/algorithm decision; current sampled contact has no IPC-like
CCD maximum-step contract. Static sampled friction is also not defined by the
current components and should remain unsupported unless a concrete model is
approved.

## Scope deliberately not started

This spike does not add `simulationRunner`, change the public C/Python APIs,
or alter existing runner output behavior. Those changes belong to the next
Phase 2 milestone after the typed config and output-safety contracts are in
place.
