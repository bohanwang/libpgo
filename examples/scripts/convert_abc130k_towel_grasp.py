#!/usr/bin/env python3
"""Convert the ABC-130k towel-grasp package into libpgo example assets.

Reads the public data directory (manifest.json, grasp_config.json, the 50 Hz
robot/grasp NPZ and the two towel OBJ files) and writes:

  towel-rest.obj         planar rest metric (copied)
  towel-initial.obj      estimated state at the first recorded frame (copied)
  table-plane.obj        z-up table plane at the manifest's table height
  grasp-trajectory.json  pgo-vertex-trajectory v1: one track per hand, the
                         recorded pin targets keyed by simulation time
                         (recorded time minus the first recorded time), and
                         the release time of each hand from the active mask
  SOURCE.md              provenance and checksums of the inputs

Requires numpy. Example:

  python3 examples/scripts/convert_abc130k_towel_grasp.py \
      --data-dir examples/abc130k-uipc-3d-trajectory-downloads-20260921/data \
      --out-dir examples/assets/shell/towel_grasp
"""
import argparse
import hashlib
import json
import os
import shutil

import numpy as np


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def dump_trajectory(trajectory):
    """JSON text with one keyframe per line, so the file stays readable."""
    head = {k: v for k, v in trajectory.items() if k != "tracks"}
    lines = [json.dumps(head, indent=1)[:-2].rstrip() + ","]
    lines.append(' "tracks": [')
    for ti, track in enumerate(trajectory["tracks"]):
        lines.append("  {")
        lines.append(f'   "name": {json.dumps(track["name"])},')
        lines.append(f'   "vertices": {json.dumps(track["vertices"])},')
        if "release-time" in track:
            lines.append(f'   "release-time": {json.dumps(track["release-time"])},')
        lines.append('   "keys": [')
        for ki, key in enumerate(track["keys"]):
            comma = "," if ki + 1 < len(track["keys"]) else ""
            lines.append(f'    {{"t": {key["t"]}, "p": {json.dumps(key["p"])}}}{comma}')
        lines.append("   ]")
        lines.append("  }" + ("," if ti + 1 < len(trajectory["tracks"]) else ""))
    lines.append(" ]")
    lines.append("}")
    return "\n".join(lines) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--run", default="baseline_fine_legacy",
                    help="manifest run whose robot/grasp NPZ is converted (the two runs share identical grasp data)")
    ap.add_argument("--table-size", type=float, default=4.0, help="side length of the square table plane, metres")
    ap.add_argument("--decimals", type=int, default=6, help="rounding of trajectory coordinates (6 = micrometres)")
    args = ap.parse_args()

    data = args.data_dir
    manifest = json.load(open(os.path.join(data, "manifest.json")))
    grasp = json.load(open(os.path.join(data, "grasp_config.json")))
    run = manifest["runs"][args.run]
    npz_name = os.path.basename(run["compact_trajectory"])
    npz_path = os.path.join(data, npz_name)
    z = np.load(npz_path)

    timeline = z["timeline_s"]
    targets = z["grasp_target_position_m"]  # (frames, 4, 3)
    active = z["active_grasp_mask"]  # (frames, 4)
    pins = [int(v) for v in z["pin_indices"]]
    names = grasp["pin_order"]
    assert pins == grasp["pin_indices_zero_based"], "pin order in NPZ and grasp_config.json differ"
    assert len(names) == 4 and targets.shape[1] == 4

    t0 = float(timeline[0])
    times = timeline - t0

    # group pins by hand: names are left_0, left_1, right_0, right_1
    hands = {}
    for j, name in enumerate(names):
        hands.setdefault(name.split("_")[0], []).append(j)

    tracks = []
    for hand, cols in hands.items():
        keys = []
        for k in range(len(times)):
            keys.append({
                "t": round(float(times[k]), 6),
                "p": [[round(float(c), args.decimals) for c in targets[k, j]] for j in cols],
            })
        track = {
            "name": f"{hand}-gripper",
            "vertices": [pins[j] for j in cols],
            "keys": keys,
        }
        inactive = np.where(~active[:, cols].all(axis=1))[0]
        if len(inactive):
            track["release-time"] = round(float(times[inactive[0]]), 6)
        # keys must come after release-time in the JSON for readability
        track = {k: track[k] for k in ("name", "vertices", "release-time", "keys") if k in track}
        tracks.append(track)

    trajectory = {
        "format": "pgo-vertex-trajectory",
        "version": 1,
        "space": "position",
        "interpolation": "linear",
        "units": "metre",
        "source": {
            "package": manifest["schema_version"],
            "episode": manifest["source_episode"],
            "run": args.run,
            "file": npz_name,
            "time-origin-s": t0,
            "record-hz": run["record_hz"],
        },
        "tracks": tracks,
    }

    os.makedirs(args.out_dir, exist_ok=True)
    with open(os.path.join(args.out_dir, "grasp-trajectory.json"), "w") as f:
        f.write(dump_trajectory(trajectory))

    shutil.copyfile(os.path.join(data, "towel_rest.obj"), os.path.join(args.out_dir, "towel-rest.obj"))
    shutil.copyfile(os.path.join(data, "towel_initial.obj"), os.path.join(args.out_dir, "towel-initial.obj"))

    # z-up table plane centred under the towel's rest footprint
    table_z = float(manifest["mesh"]["table_z_m"])
    cx, cy = 0.7, 0.0  # under the towel's rest footprint (x in [0.3, 1.15], y in [-0.19, 0.19])
    h = 0.5 * args.table_size
    with open(os.path.join(args.out_dir, "table-plane.obj"), "w") as f:
        f.write(f"# table plane, z-up, z = {table_z} m (manifest table_z_m)\n")
        for x, y in ((cx - h, cy - h), (cx + h, cy - h), (cx + h, cy + h), (cx - h, cy + h)):
            f.write(f"v {x:.6f} {y:.6f} {table_z:.6f}\n")
        f.write("f 1 2 3\nf 1 3 4\n")

    with open(os.path.join(args.out_dir, "SOURCE.md"), "w") as f:
        f.write("# Towel grasp assets\n\n")
        f.write(f"Converted from the `{manifest['schema_version']}` package, episode `{manifest['source_episode']}`, ")
        f.write(f"task `{manifest['source_task']}`, run `{args.run}`, by `examples/scripts/convert_abc130k_towel_grasp.py`.\n\n")
        f.write("| Output | Input | Input sha256 |\n| --- | --- | --- |\n")
        for out, src in (("towel-rest.obj", "towel_rest.obj"), ("towel-initial.obj", "towel_initial.obj"),
                         ("grasp-trajectory.json", npz_name)):
            f.write(f"| `{out}` | `{src}` | `{sha256(os.path.join(data, src))}` |\n")
        f.write(f"| `table-plane.obj` | manifest `table_z_m = {table_z}` | |\n\n")
        f.write(f"Trajectory times are recorded time minus {t0} s; {len(times)} keys per track at {run['record_hz']} Hz.\n")
        for tr in tracks:
            f.write(f"- `{tr['name']}`: vertices {tr['vertices']}, release at {tr.get('release-time', 'never')} s\n")

    print(f"wrote {len(tracks)} tracks x {len(times)} keys to {args.out_dir}")


if __name__ == "__main__":
    main()
