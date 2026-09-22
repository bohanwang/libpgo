# Towel grasp (ABC-130k episode replay)

A two-arm robot lifts a folded towel by four grasped vertices, carries it, and
releases it onto a table. The towel is a shell simulated with IPC contact
against a static table plane and itself; the grippers are prescribed vertex
trajectories recorded from the source episode. There is no robot geometry in
the source package, so the arms are not collision objects here either.

## Source

`examples/scripts/convert_abc130k_towel_grasp.py` converts the public
`abc130k-uipc-organized-v1` package (episode
`episode_0104361e-4e0e-4254-83f2-80c59a7038cd`, task
`fold_and_stack_the_towels`; robot replayed in MuJoCo, cloth in libuipc) into
`examples/assets/shell/towel_grasp/`. See `SOURCE.md` there for checksums.
The source marks its material as uncalibrated and its runs as not physics
qualified; the per-frame cloth states are not part of the public package, so
nothing here is a ground-truth comparison.

| Quantity | Value |
| --- | --- |
| Towel grid | 49 x 25 vertices, 2304 triangles, 0.85 x 0.38 m rest |
| Mass, thickness | 65 g (areal density 0.2012 kg/m2), 1.6 mm |
| Table | plane at z = 0.75; towel rests 4 mm above it |
| Initial state | `towel-initial.obj`, folded hairpin at the first recorded frame |
| Grasp | two vertices per hand on the far edge: 1221, 1222 (left), 1202, 1203 (right) |
| Trajectory | 826 keys at 50 Hz, peak 0.9 m/s, lift to z about 1.49 |
| Release | right hand at 11.32 s, left at 11.40 s; free fall to the table afterwards |
| Frame | MuJoCo world, z up, metres |

## Config notes

- `initial-surface-mesh` supplies the folded start state as a per-vertex
  displacement from the planar rest mesh.
- The `fixed-vertices` entry uses a `trajectory` file. Frame k is the state at
  t = k * 0.002 s; 8251 steps cover the 16.5 s timeline and `dump-interval`
  10 writes 826 OBJ frames aligned with the recorded 50 Hz frames.
- Material: the source package lists stretch and shear moduli of 4000 Pa for
  its visual reference, 1e5 Pa for its rejected candidate, and a bending
  stiffness of 1.42e-6 N m, and marks all of them as uncalibrated. Measured
  cotton fabrics of similar weight give a different picture. Wang, O'Brien and
  Ramamoorthi (SIGGRAPH 2011) list a 60 percent cotton knit at 0.187 kg/m2,
  close to this towel's 0.201 kg/m2, whose ARCSim material file has bending
  stiffness coefficients of 7e-6 to 39e-6 N m and stretching coefficients of
  tens to a few hundred N/m in the linear regime; other cloth-simulation
  measurements put cotton stretch stiffness at about 1000 to 2000 N/m and
  Kawabata bending rigidity of light woven cotton at 0.05 to 0.2 gf cm2/cm,
  which is 5e-6 to 2e-5 N m. Measured pseudo Poisson ratios of woven cotton
  are 0.2 to 0.5. This config targets a membrane stiffness E h of 1000 N/m,
  a bending rigidity of 1e-5 N m and an in-plane Poisson ratio of 0.4 (the
  source's bending value is about seven times softer than that range).
- Mapping to the config keys: libpgo's Koiter St. Venant-Kirchhoff model
  converts each modulus to 3D Lame constants with `shell-poisson-ratio` and
  applies them to twice the Green strain, so the config values are not the
  physical ones. With config values E, nu and E_b the effective in-plane
  Young's modulus is 4 E / (1 - nu^2), the effective in-plane Poisson ratio is
  nu / (1 - nu), and the bending rigidity for a cylindrical bend is
  E_b (1 - nu) / ((1 + nu)(1 - 2 nu)) h^3 / 12. Inverting these for the
  targets above with h = 1.6 mm gives `shell-poisson-ratio` 0.29,
  `shell-youngs-modulus` 1.4e5 and `shell-bending-youngs-modulus` 2.3e4.
  To reproduce the source's bending stiffness instead, use 3200 for the
  bending modulus.
- Contact: `ipc-heuristic` picks dhat as 1e-3 of the bounding-box diagonal
  (about 0.9 mm, below the 4 mm table clearance and the 22 mm layer gap) and
  kappa from the membrane stiffness. Friction 0.65 applies to the table and
  to self-contact, as in the source.
- `coeff` 1e4 N/m holds each grasped vertex to within about 20 um under the
  towel's weight.

## Run

```bash
build/all_release/bin/runIPCSim examples/configs/shell/towel_grasp/towel-grasp-dynamic-ipc.json
mkdir -p examples/generated/output/towel-grasp-dynamic-ipc-abc
uv run pgo-dump-abc examples/configs/shell/towel_grasp/towel-grasp-dynamic-ipc-animation.json \
  examples/generated/output/towel-grasp-dynamic-ipc-abc
```
