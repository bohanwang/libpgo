# libpgo 0.0.4 Remaining Test Plan

This document tracks only the validation work that remains before the 0.0.4
release. Implementation history and superseded build/package designs are not
release guidance and are intentionally omitted.

## Verified baseline

The current baseline is commit
`ff654876b997489cc4a38f90b5a7bb6c5b780371` on branch `add-no-conda`.

The build-tree tests, repaired-wheel audit, installed-wheel smoke tests, and
installed `tests/pypgo/test_pgo_smoke.py` suite pass in all three hosted
workflows:

- [Linux CI](https://github.com/annajcy/libpgo/actions/runs/30682955487)
- [macOS CI](https://github.com/annajcy/libpgo/actions/runs/30682955507)
- [Windows CI](https://github.com/annajcy/libpgo/actions/runs/30682955497)

These runs establish the packaging baseline. They do not by themselves close
the compatibility, sanitizer, or owner-run Linux server gates below.

## Test policy

- Bind every result to an exact commit and a clean worktree.
- Reuse existing tests when they already cover a requirement; do not duplicate
  a test only to give it a release-specific name.
- Use lightweight assets and write generated meshes and simulation output
  outside the source tree.
- A simulation smoke passes only when it returns success, produces finite
  state/output, and exercises the intended mesh, contact, and solve path.
- Keep long numerical runs off hosted CI. Run them manually on the Linux server.
- After a relevant code, config, asset-generator, solver, or dependency change,
  invalidate and rerun the affected evidence.

## P5.1 Audit native regression coverage

First map each requirement below to an existing test and CI result. Add only the
smallest missing regression.

- [ ] Cubic energy, gradient, and Hessian finite differences.
- [ ] Cubic invalid rest-element validation.
- [ ] Cubic mesher deterministic topology and basic geometry.
- [ ] Lightweight cubic generator arguments, no-overwrite behavior, manifest,
      and repeatability.
- [ ] Tet and cubic material maximum-step behavior.
- [ ] IPC CCD maximum-step behavior.
- [ ] IPC barrier energy, gradient, and Hessian finite differences.
- [ ] IPC parameter, topology, state-vector, and spatial-hash validation.
- [ ] Mixed fixed/dynamic Hessian assembly.
- [ ] Newton fixed-pattern reuse and dynamic-pattern safety.
- [ ] Generic tet/cubic volume-loader equivalence.
- [ ] Config dispatch, missing `contact-model` fallback, and invalid routing.
- [ ] Missing or unreadable config failure through CLI, C, and Python.
- [ ] C API Hessian export preserves finite fractional `double` values.
- [ ] Sampled and IPC shell regressions retain their established runners.

Record the source test and command beside each checked item or in the final test
report.

## P5.2 Complete the lightweight API and runner matrix

Use dedicated fixtures under `tests/`. Each supported cell must be exercised
through the shared simulation runner from both the C and Python APIs. The CLI
column verifies that the public executable still routes to the same path.

| Case | C API | Python API | CLI | Required evidence |
| --- | --- | --- | --- | --- |
| tet + sampled + dynamic | [ ] | [ ] | `runSim` [ ] | finite output, return 0 |
| tet + sampled + static | [ ] | [ ] | `runSim` [ ] | finite output, return 0 |
| tet + IPC + dynamic | [ ] | [ ] | `runIPCSim` [ ] | IPC Hessian and CCD path |
| tet + IPC + static | [ ] | [ ] | `runIPCSim` [ ] | finite static solve |
| cubic + sampled + dynamic | [ ] | [ ] | `runSim` [ ] | cubic FEM path |
| cubic + sampled + static | [ ] | [ ] | `runSim` [ ] | warning that sampled contact energy is not included |
| cubic + IPC + dynamic | [ ] | [ ] | `runIPCSim` [ ] | cubic, IPC Hessian, and CCD path |
| cubic + IPC + static | [ ] | [ ] | `runIPCSim` [ ] | finite static solve |
| shell + IPC | [ ] | [ ] | `runIPCSim` [ ] | shell IPC path |
| shell + sampled | N/A | N/A | `runShellSim` [ ] | existing sampled shell path |
| missing `contact-model` | [ ] | [ ] | `runSim` [ ] | sampled fallback |
| invalid routing | [ ] | [ ] | N/A | useful error and nonzero result |

For sampled + static, the test must verify the explicit warning; it must not
claim that sampled contact contributes to the Newton energy.

## P5.3 Compatibility and installed-consumer checks

- [ ] Compare one canonical 0.0.3 tet + sampled case with 0.0.4 using recorded
      displacement, energy, or output tolerances.
- [ ] Compare exported C symbols on Linux, macOS, and Windows with the 0.0.3
      reference and explain every intentional difference.
- [ ] Run applicable 0.0.3 Python smoke tests against 0.0.4.
- [ ] Verify the `const char *` cleanup did not change the C ABI.
- [ ] Install into a temporary prefix and verify `pgo_c.h`, `pgo_c_def.h`, and
      the CMake package export are present.
- [ ] Build and run a minimal external pure-C consumer using only the installed
      prefix.
- [ ] Verify `find_package(pgo 0.0.4 CONFIG REQUIRED)` without source-tree or
      build-tree paths.
- [ ] Import and test every wheel from outside the checkout in a clean runtime
      environment.

## P5.4 Sanitizer and numerical checks

- [ ] Run the native suite in a clean Linux ASan/UBSan build where supported by
      the dependencies.
- [ ] Reject NaN or infinity in every matrix smoke output.
- [ ] Check sparse Hessian dimensions and finite entries.
- [ ] Verify every generated cubic rest determinant is finite and positive.
- [ ] Finish with a clean Release build and rerun the affected smoke tests.

## P5.5 Owner-run Linux MKL Pardiso validation

This manual gate is intentionally outside GitHub Actions.

### Environment evidence

- [ ] Build Release with `PGO_USE_MKL=ON`, `PGO_HAS_ORIG_PARDISO=OFF`, and the
      TBB threading layer.
- [ ] Confirm the Newton solver selects `EigenMKLPardisoSupport`.
- [ ] Confirm loaded-library evidence contains `mkl_tbb_thread` and TBB and
      excludes the Intel OpenMP threading layer.
- [ ] Record the commit, CMake cache, compiler, MKL/TBB/Python versions,
      environment variables, thread policy, and `ldd` output.

### Tier 1: short integration matrix

- [ ] Run the supported P5.2 tet, cubic, and shell cases through C, Python, and
      their public CLIs.
- [ ] Record command, return code, wall time, output count, and failure reason
      for every cell.

### Tier 2: representative numerical matrix

- [ ] Run representative tet + sampled, tet + IPC, cubic + sampled,
      cubic + IPC, and shell + IPC simulations.
- [ ] Require successful processes, finite output, solver convergence, valid
      rest elements, expected frame counts, and bounded displacement.
- [ ] Treat timing as diagnostic information, not a performance promise.

### Report

- [ ] Store the exact materialized configs, generated-asset manifest and hashes,
      commands, logs, environment evidence, and results together.
- [ ] Produce a machine-readable report plus a concise Markdown summary.
- [ ] Tie the report to the exact release-candidate commit.
- [ ] After a correctness change, rerun the affected cell and the complete Tier
      2 matrix.

## Completion criteria

Phase 5 is complete when:

- every item above is checked with evidence or explicitly marked not applicable
  with a reason;
- the three hosted workflows pass on the final candidate commit;
- the external C/CMake consumer passes;
- sanitizer and clean Release results pass; and
- the owner confirms the Linux MKL Pardiso report for that same commit.

Completing this plan does not authorize merging, tagging, publishing artifacts,
or deploying documentation.
