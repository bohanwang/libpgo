#!/usr/bin/env python3
"""Generate genuinely lightweight tet/cubic assets for example smoke tests."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import shutil
import subprocess
import tempfile
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT_DIR = REPO_ROOT / "examples" / "generated" / "smoke-assets"
SCENE_SURFACES = {
    "box": REPO_ROOT / "examples/assets/volume/box/box.obj",
    "box-with-sphere": (
        REPO_ROOT / "examples/assets/volume/box-with-sphere/box-with-sphere.obj"
    ),
    "bunny": REPO_ROOT / "examples/assets/volume/bunny/bunny.obj",
    "dragon": REPO_ROOT / "examples/assets/volume/dragon/dragon.obj",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def find_binary(explicit: Path | None, build_dir: Path | None, name: str) -> Path:
    if explicit is not None:
        candidates = [explicit]
    else:
        assert build_dir is not None
        suffixes = (f"{name}.exe", name) if os.name == "nt" else (name,)
        candidates = [
            build_dir / directory / suffix
            for directory in ("bin", "src/tools", f"src/tools/{name}", "")
            for suffix in suffixes
        ]
    for candidate in candidates:
        resolved = candidate.expanduser().resolve()
        if resolved.is_file():
            return resolved
    raise FileNotFoundError(
        f"Could not find {name}; pass --{name.replace('Mesher', '-mesher')} PATH "
        "or --build-dir DIR."
    )


def obj_bounds(path: Path) -> tuple[list[float], list[float]]:
    vertices = []
    for line in path.read_text().splitlines():
        if line.startswith("v "):
            fields = line.split()
            vertices.append(tuple(float(value) for value in fields[1:4]))
    if not vertices:
        raise ValueError(f"{path} contains no OBJ vertices")
    bounds_min = [min(point[axis] for point in vertices) for axis in range(3)]
    bounds_max = [max(point[axis] for point in vertices) for axis in range(3)]
    if any(
        not math.isfinite(value)
        for value in bounds_min + bounds_max
    ) or any(bounds_min[axis] >= bounds_max[axis] for axis in range(3)):
        raise ValueError(f"{path} has invalid bounds")
    return bounds_min, bounds_max


def write_box_surface(
    path: Path, bounds_min: list[float], bounds_max: list[float]
) -> None:
    x0, y0, z0 = bounds_min
    x1, y1, z1 = bounds_max
    vertices = (
        (x0, y0, z0),
        (x1, y0, z0),
        (x1, y1, z0),
        (x0, y1, z0),
        (x0, y0, z1),
        (x1, y0, z1),
        (x1, y1, z1),
        (x0, y1, z1),
    )
    faces = (
        (1, 3, 2), (1, 4, 3),
        (5, 6, 7), (5, 7, 8),
        (1, 2, 6), (1, 6, 5),
        (4, 8, 7), (4, 7, 3),
        (1, 5, 8), (1, 8, 4),
        (2, 3, 7), (2, 7, 6),
    )
    lines = ["# Eight-vertex bounding-box surface for libpgo smoke tests."]
    lines.extend(f"v {x:.17g} {y:.17g} {z:.17g}" for x, y, z in vertices)
    lines.extend(f"f {a} {b} {c}" for a, b, c in faces)
    path.write_text("\n".join(lines) + "\n")


def parse_veg(path: Path, element_type: str, arity: int) -> dict:
    lines = [line.strip() for line in path.read_text().splitlines()]
    try:
        vertices_marker = lines.index("*VERTICES")
        elements_marker = lines.index("*ELEMENTS")
    except ValueError as error:
        raise ValueError(f"{path} is missing a required VEG section") from error

    vertex_count = int(lines[vertices_marker + 1].split()[0])
    vertices = []
    for line in lines[vertices_marker + 2 : vertices_marker + 2 + vertex_count]:
        fields = line.split()
        if len(fields) != 4:
            raise ValueError(f"Malformed VEG vertex row: {line}")
        point = tuple(float(value) for value in fields[1:])
        if not all(math.isfinite(value) for value in point):
            raise ValueError(f"{path} contains a non-finite vertex")
        vertices.append(point)

    if lines[elements_marker + 1] != element_type:
        raise ValueError(f"{path} does not contain {element_type} elements")
    element_count = int(lines[elements_marker + 2].split()[0])
    for line in lines[elements_marker + 3 : elements_marker + 3 + element_count]:
        fields = line.split()
        if len(fields) != arity + 1:
            raise ValueError(f"Malformed VEG element row: {line}")
        indices = [int(value) - 1 for value in fields[1:]]
        if len(set(indices)) != arity or min(indices) < 0 or max(indices) >= vertex_count:
            raise ValueError(f"{path} contains invalid element indices")

    if vertex_count <= 0 or element_count <= 0:
        raise ValueError(f"{path} must contain vertices and elements")
    return {
        "vertices": vertices,
        "vertex_count": vertex_count,
        "element_count": element_count,
        "element_arity": arity,
        "bounds": {
            "min": [min(point[axis] for point in vertices) for axis in range(3)],
            "max": [max(point[axis] for point in vertices) for axis in range(3)],
        },
    }


def boundary_indices(vertices: list[tuple[float, float, float]], axis: int, high: bool):
    values = [point[axis] for point in vertices]
    target = max(values) if high else min(values)
    span = max(values) - min(values)
    tolerance = max(1.0, span) * 1.0e-9
    return [
        index for index, point in enumerate(vertices)
        if abs(point[axis] - target) <= tolerance
    ]


def write_indices(path: Path, indices: list[int]) -> None:
    if not indices:
        raise ValueError(f"Cannot write empty constraint set: {path}")
    path.write_text("".join(f"{index}\n" for index in indices))


def revision() -> str:
    result = subprocess.run(
        ["git", "-C", str(REPO_ROOT), "rev-parse", "HEAD"],
        check=False,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip() if result.returncode == 0 else "unknown"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", type=Path)
    parser.add_argument("--tet-mesher", type=Path)
    parser.add_argument("--cubic-mesher", type=Path)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"New output directory (default: {DEFAULT_OUTPUT_DIR})",
    )
    parser.add_argument(
        "--cubic-resolution",
        type=int,
        default=2,
        help="Cubic cells along the shortest bounding-box side (default: 2)",
    )
    parser.add_argument(
        "--tetgen-command",
        default="pq1.2Y",
        help="TetGen command passed through tetMesher (default: pq1.2Y)",
    )
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    try:
        if args.build_dir is None and (
            args.tet_mesher is None or args.cubic_mesher is None
        ):
            raise ValueError(
                "Pass --build-dir, or pass both --tet-mesher and --cubic-mesher"
            )
        if args.cubic_resolution <= 0:
            raise ValueError("--cubic-resolution must be positive")
        build_dir = args.build_dir.expanduser().resolve() if args.build_dir else None
        tet_mesher = find_binary(args.tet_mesher, build_dir, "tetMesher")
        cubic_mesher = find_binary(args.cubic_mesher, build_dir, "cubicMesher")
        output_dir = args.output_dir.expanduser().resolve()
        if output_dir.exists():
            raise FileExistsError(
                f"Refusing to overwrite existing output directory: {output_dir}"
            )
        for surface in SCENE_SURFACES.values():
            if not surface.is_file():
                raise FileNotFoundError(f"Missing canonical surface: {surface}")

        output_dir.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(
            prefix=f".{output_dir.name}-", dir=output_dir.parent
        ) as temporary:
            staging = Path(temporary) / output_dir.name
            staging.mkdir()
            records = []
            for scene, canonical_surface in SCENE_SURFACES.items():
                bounds_min, bounds_max = obj_bounds(canonical_surface)
                surface = staging / f"{scene}-surface.obj"
                tet = staging / f"{scene}-tet.veg"
                cubic = staging / f"{scene}-cubic.veg"
                write_box_surface(surface, bounds_min, bounds_max)

                tet_command = [
                    str(tet_mesher), "tetgen",
                    "--input-mesh", str(surface),
                    "--output-mesh", str(tet),
                    "--command", args.tetgen_command,
                ]
                cubic_command = [
                    str(cubic_mesher),
                    "--input-mesh", str(surface),
                    "--resolution", str(args.cubic_resolution),
                    "--output-mesh", str(cubic),
                    "--E", "1000000",
                    "--nu", "0.45",
                    "--density", "1000",
                ]
                print(" ".join(tet_command), flush=True)
                subprocess.run(tet_command, check=True)
                print(" ".join(cubic_command), flush=True)
                subprocess.run(cubic_command, check=True)

                meshes = {}
                for domain, mesh, element_type, arity in (
                    ("tet", tet, "TET", 4),
                    ("cubic", cubic, "CUBIC", 8),
                ):
                    metadata = parse_veg(mesh, element_type, arity)
                    vertices = metadata.pop("vertices")
                    constraint_names = {"hang": [0]}
                    if scene == "box":
                        constraint_names.update({
                            "zmin": boundary_indices(vertices, 2, False),
                            "zmax": boundary_indices(vertices, 2, True),
                        })
                    constraints = {}
                    for kind, indices in constraint_names.items():
                        filename = f"{scene}-{domain}-{kind}.txt"
                        write_indices(staging / filename, indices)
                        constraints[kind] = {
                            "file": filename,
                            "count": len(indices),
                            "sha256": sha256(staging / filename),
                        }
                    metadata["constraints"] = constraints
                    metadata["file"] = mesh.name
                    metadata["sha256"] = sha256(mesh)
                    meshes[domain] = metadata

                records.append({
                    "scene": scene,
                    "canonical_surface": str(canonical_surface),
                    "canonical_surface_sha256": sha256(canonical_surface),
                    "smoke_surface": surface.name,
                    "smoke_surface_sha256": sha256(surface),
                    "canonical_bounds": {"min": bounds_min, "max": bounds_max},
                    "meshes": meshes,
                })

            manifest = {
                "schema_version": 1,
                "purpose": "lightweight generated assets for all-config smoke tests",
                "generator": str(Path(__file__).resolve()),
                "source_revision": revision(),
                "tet_mesher": str(tet_mesher),
                "cubic_mesher": str(cubic_mesher),
                "parameters": {
                    "surface": "eight-vertex canonical AABB",
                    "tetgen_command": args.tetgen_command,
                    "cubic_resolution": args.cubic_resolution,
                },
                "assets": records,
            }
            (staging / "manifest.json").write_text(
                json.dumps(manifest, indent=2) + "\n"
            )
            shutil.move(str(staging), str(output_dir))

        print(f"Wrote {len(records)} tet/cubic smoke asset pairs to {output_dir}")
        return 0
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        parser.exit(1, f"error: {error}\n")


if __name__ == "__main__":
    raise SystemExit(main())
