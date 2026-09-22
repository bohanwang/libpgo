# Towel grasp assets

Converted from the `abc130k-uipc-organized-v1` package, episode `episode_0104361e-4e0e-4254-83f2-80c59a7038cd`, task `fold_and_stack_the_towels`, run `baseline_fine_legacy`, by `examples/scripts/convert_abc130k_towel_grasp.py`.

| Output | Input | Input sha256 |
| --- | --- | --- |
| `towel-rest.obj` | `towel_rest.obj` | `607921bfae82cdbf2b2bd90a58cc164a8d02b08f2f216e695f11b60eae8cfada` |
| `towel-initial.obj` | `towel_initial.obj` | `b4d319ee3440a551ee87311559b88ce203397d716647f9497fba28cffdfeb352` |
| `grasp-trajectory.json` | `baseline_fine_legacy_robot_grasp_50hz.npz` | `ab0deb510b685f747007d6a32d947fdf3d5c586b371130a6e4e92e545e523d6e` |
| `table-plane.obj` | manifest `table_z_m = 0.75` | |

Trajectory times are recorded time minus 10.5 s; 826 keys per track at 50 Hz.
- `left-gripper`: vertices [1221, 1222], release at 11.4 s
- `right-gripper`: vertices [1202, 1203], release at 11.32 s
