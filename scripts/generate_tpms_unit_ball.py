#!/usr/bin/env python3
"""Generate classical TPMS shells in [-1, 1]^3 for use with generateFBMSUnionSurface.

For each TPMS (Schwarz P, Schwarz D, Gyroid, I-WP, Neovius) we sample
f_tpms(x, y, z) on a uniform grid over [-1, 1]^3 and march out the open
isosurface f_tpms == 0 -- the "shell". The shell extends to the cube
boundary; clipping it down to a ball is the job of generateFBMSUnionSurface
with --enable-truncating, which CSG-intersects the thickened shell against
the outer surface of a thickened bounding sphere.

For convenience we also drop a matching unit-sphere mesh into each folder,
with a header that loadSphereParameters() in generateFBMSUnionSurface.cpp
can parse (# center = ..., # radius = ...).

By default, each TPMS lands in its own subdirectory under examples/fbms with
files named per the existing fbms convention:
    <prefix>_fbms.obj
    <prefix>_fbms_bounding_sphere.obj

Pass --flat-output to write all files directly into the output folder.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np
import trimesh
from skimage import measure


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = REPO_ROOT / "examples" / "fbms"


@dataclass(frozen=True)
class TPMS:
    name: str
    folder: str
    field: Callable[[np.ndarray, np.ndarray, np.ndarray, float], np.ndarray]


# w = omega is the angular frequency. Period along each axis is 2*pi/w; for
# `cells` periods across the [-1, 1] domain (length 2), use w = pi * cells.
def schwarz_p(x, y, z, w):
    return np.cos(w * x) + np.cos(w * y) + np.cos(w * z)


def schwarz_d(x, y, z, w):
    return (
        np.cos(w * x) * np.cos(w * y) * np.cos(w * z)
        - np.sin(w * x) * np.sin(w * y) * np.sin(w * z)
    )


def gyroid(x, y, z, w):
    return (
        np.sin(w * x) * np.cos(w * y)
        + np.sin(w * y) * np.cos(w * z)
        + np.sin(w * z) * np.cos(w * x)
    )


def iwp(x, y, z, w):
    return 2.0 * (
        np.cos(w * x) * np.cos(w * y)
        + np.cos(w * y) * np.cos(w * z)
        + np.cos(w * z) * np.cos(w * x)
    ) - (np.cos(2.0 * w * x) + np.cos(2.0 * w * y) + np.cos(2.0 * w * z))


def neovius(x, y, z, w):
    return 3.0 * (np.cos(w * x) + np.cos(w * y) + np.cos(w * z)) + 4.0 * (
        np.cos(w * x) * np.cos(w * y) * np.cos(w * z)
    )


TPMS_LIST: list[TPMS] = [
    TPMS("Schwarz P", "tpms_schwarz_p", schwarz_p),
    TPMS("Schwarz D", "tpms_schwarz_d", schwarz_d),
    TPMS("Gyroid", "tpms_gyroid", gyroid),
    TPMS("I-WP", "tpms_iwp", iwp),
    TPMS("Neovius", "tpms_neovius", neovius),
]


def build_tpms_field(
    tpms_field: Callable[[np.ndarray, np.ndarray, np.ndarray, float], np.ndarray],
    resolution: int,
    cells: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (field, spacing, origin) on a uniform grid over [-1, 1]^3."""
    axis = np.linspace(-1.0, 1.0, resolution, dtype=np.float64)
    spacing = np.full(3, axis[1] - axis[0])
    origin = np.array([axis[0], axis[0], axis[0]])

    # ij indexing keeps axes in (x, y, z) order to match marching_cubes below.
    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")

    omega = math.pi * cells
    return tpms_field(x, y, z, omega), spacing, origin


def marching_cubes_to_mesh(
    field: np.ndarray, spacing: np.ndarray, origin: np.ndarray
) -> trimesh.Trimesh:
    verts, faces, _, _ = measure.marching_cubes(
        field,
        level=0.0,
        spacing=tuple(spacing),
        allow_degenerate=False,
    )
    verts = verts + origin  # shift from grid coords to world coords
    return trimesh.Trimesh(vertices=verts, faces=faces, process=True)


def write_obj_with_header(mesh: trimesh.Trimesh, path: Path, header_lines: list[str]) -> None:
    obj_text = trimesh.exchange.obj.export_obj(mesh, include_normals=False, include_color=False)
    with path.open("w", encoding="utf-8") as fp:
        for line in header_lines:
            fp.write(f"# {line}\n")
        fp.write(obj_text)


def write_unit_sphere_obj(path: Path, subdivisions: int) -> tuple[int, int]:
    """Write a unit sphere mesh with parseable center/radius header."""
    sphere = trimesh.creation.icosphere(subdivisions=subdivisions, radius=1.0)
    vertices = np.asarray(sphere.vertices, dtype=np.float64)
    faces = np.asarray(sphere.faces, dtype=np.int64)

    with path.open("w", encoding="utf-8") as fp:
        fp.write("# unit sphere for TPMS+sphere CSG via generateFBMSUnionSurface\n")
        fp.write("# center = (0, 0, 0)\n")
        fp.write("# radius = 1\n")
        fp.write(f"# sphere_resolution = icosphere subdivisions {subdivisions}\n")
        fp.write(f"# sphere_vertices = {len(vertices)}\n")
        fp.write(f"# sphere_faces = {len(faces)}\n")
        fp.write("o unit_sphere\n")
        for x, y, z in vertices:
            fp.write(f"v {x:.17g} {y:.17g} {z:.17g}\n")
        for i, j, k in faces:
            fp.write(f"f {i + 1} {j + 1} {k + 1}\n")
    return len(vertices), len(faces)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cells", type=float, default=1.0,
                        help="Number of TPMS unit cells across the [-1, 1] domain (default: 1).")
    parser.add_argument("--resolution", type=int, default=256,
                        help="Grid samples per axis for marching cubes (default: 256).")
    parser.add_argument("--sphere-subdivisions", type=int, default=4,
                        help="Icosphere subdivisions for the bounding-sphere mesh (default: 4).")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT,
                        help="Parent output folder (default: examples/fbms).")
    parser.add_argument("--flat-output", action="store_true",
                        help="Write all TPMS assets directly under --out instead of one subfolder per TPMS.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.resolution < 8:
        raise SystemExit("--resolution must be at least 8")
    if args.cells <= 0:
        raise SystemExit("--cells must be positive")

    args.out = args.out.resolve()
    args.out.mkdir(parents=True, exist_ok=True)

    print(f"[tpms] resolution={args.resolution} cells={args.cells} out={args.out}")
    for tpms in TPMS_LIST:
        sub_dir = args.out if args.flat_output else args.out / tpms.folder
        sub_dir.mkdir(parents=True, exist_ok=True)

        field, spacing, origin = build_tpms_field(tpms.field, args.resolution, args.cells)
        f_min, f_max = float(field.min()), float(field.max())
        if not (f_min <= 0.0 <= f_max):
            print(f"[tpms] {tpms.name}: zero level set absent (range {f_min:g}..{f_max:g}); skipping.")
            continue

        mesh = marching_cubes_to_mesh(field, spacing, origin)
        shell_path = sub_dir / f"{tpms.folder}_fbms.obj"
        write_obj_with_header(
            mesh,
            shell_path,
            [
                f"TPMS = {tpms.name} open shell (f == 0) in [-1, 1]^3",
                f"cells = {args.cells} (omega = pi*{args.cells})",
                f"resolution = {args.resolution}",
                f"vertices = {len(mesh.vertices)}",
                f"faces = {len(mesh.faces)}",
                "intended use: feed --fbms into generateFBMSUnionSurface;",
                "  pair with this folder's *_fbms_bounding_sphere.obj as --sphere",
                "  and pass --enable-truncating to clip to the unit ball.",
            ],
        )

        sphere_path = sub_dir / f"{tpms.folder}_fbms_bounding_sphere.obj"
        sphere_v, sphere_f = write_unit_sphere_obj(sphere_path, args.sphere_subdivisions)

        print(
            f"[tpms] {tpms.name:>10s} -> {shell_path.relative_to(REPO_ROOT)} "
            f"(V={len(mesh.vertices)}, F={len(mesh.faces)}); "
            f"sphere V={sphere_v}, F={sphere_f}"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
