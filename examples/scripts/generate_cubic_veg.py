#!/usr/bin/env python3
"""Generate validated cubic VEG meshes from the canonical example surfaces."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT_DIR = REPO_ROOT / "examples" / "generated" / "cubic"


@dataclass(frozen=True)
class Scene:
    surface: Path
    resolution: int
    youngs_modulus: float
    poisson_ratio: float = 0.45
    density: float = 1000.0


SCENES = {
    "box": Scene(REPO_ROOT / "examples/assets/volume/box/box.obj", 4, 1.0e7),
    "box-with-sphere": Scene(
        REPO_ROOT / "examples/assets/volume/box-with-sphere/box-with-sphere.obj",
        50,
        1.0e6,
    ),
    "box-with-sphere-lite": Scene(
        REPO_ROOT / "examples/assets/volume/box-with-sphere/box-with-sphere.obj",
        5,
        1.0e6,
    ),
    "bunny": Scene(REPO_ROOT / "examples/assets/volume/bunny/bunny.obj", 20, 1.0e5),
    "dragon": Scene(REPO_ROOT / "examples/assets/volume/dragon/dragon.obj", 20, 1.0e6),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def find_mesher(executable: Path | None, build_dir: Path | None) -> Path:
    if executable is not None:
        candidates = [executable]
    else:
        assert build_dir is not None
        names = ("cubicMesher.exe", "cubicMesher") if os.name == "nt" else ("cubicMesher",)
        candidates = [
            build_dir / directory / name
            for directory in ("bin", "src/tools", "src/tools/cubicMesher", "")
            for name in names
        ]

    for candidate in candidates:
        resolved = candidate.expanduser().resolve()
        if resolved.is_file():
            return resolved
    raise FileNotFoundError(
        "Could not find cubicMesher; pass --mesher PATH or --build-dir DIR."
    )


def parse_cubic_veg(path: Path) -> dict:
    lines = [line.strip() for line in path.read_text().splitlines()]
    try:
        vertices_marker = lines.index("*VERTICES")
        elements_marker = lines.index("*ELEMENTS")
    except ValueError as error:
        raise ValueError(f"{path} is missing a required VEG section") from error

    vertex_count = int(lines[vertices_marker + 1].split()[0])
    if vertex_count <= 0:
        raise ValueError("VEG must contain at least one vertex")
    vertices = []
    for line in lines[vertices_marker + 2 : vertices_marker + 2 + vertex_count]:
        fields = line.split()
        if len(fields) != 4:
            raise ValueError(f"Malformed VEG vertex row: {line}")
        point = tuple(float(value) for value in fields[1:])
        if not all(math.isfinite(value) for value in point):
            raise ValueError("VEG contains a non-finite vertex")
        vertices.append(point)

    if lines[elements_marker + 1] != "CUBIC":
        raise ValueError(f"{path} does not contain CUBIC elements")
    element_count = int(lines[elements_marker + 2].split()[0])
    if element_count <= 0:
        raise ValueError("VEG must contain at least one cubic element")
    elements = []
    for line in lines[elements_marker + 3 : elements_marker + 3 + element_count]:
        fields = line.split()
        if len(fields) != 9:
            raise ValueError(f"Malformed VEG cubic element row: {line}")
        indices = tuple(int(value) - 1 for value in fields[1:])
        if len(set(indices)) != 8 or min(indices) < 0 or max(indices) >= vertex_count:
            raise ValueError("VEG contains an invalid cubic element index set")
        elements.append(indices)

    tolerance = 1.0e-9
    min_cell_size = math.inf
    max_cell_size = 0.0
    for indices in elements:
        points = [vertices[index] for index in indices]
        coordinates = [
            sorted({round(point[axis], 12) for point in points}) for axis in range(3)
        ]
        if any(len(values) != 2 for values in coordinates):
            raise ValueError("A cubic element is not an axis-aligned eight-corner cell")
        sizes = [values[1] - values[0] for values in coordinates]
        if min(sizes) <= 0.0 or max(sizes) - min(sizes) > tolerance * max(1.0, max(sizes)):
            raise ValueError("A cubic element has an invalid rest-state cell shape")
        min_cell_size = min(min_cell_size, *sizes)
        max_cell_size = max(max_cell_size, *sizes)

    bounds = {
        "min": [min(point[axis] for point in vertices) for axis in range(3)],
        "max": [max(point[axis] for point in vertices) for axis in range(3)],
    }
    return {
        "vertex_count": vertex_count,
        "element_count": element_count,
        "element_arity": 8,
        "bounds": bounds,
        "cell_size_range": [min_cell_size, max_cell_size],
        "validation": {
            "finite_vertices": True,
            "valid_indices": True,
            "axis_aligned_positive_rest_cells": True,
        },
    }


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
    mesher = parser.add_mutually_exclusive_group(required=True)
    mesher.add_argument("--mesher", type=Path, help="Path to the cubicMesher executable")
    mesher.add_argument("--build-dir", type=Path, help="Build tree containing cubicMesher")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"New output directory (default: {DEFAULT_OUTPUT_DIR})",
    )
    parser.add_argument(
        "--scene",
        action="append",
        choices=sorted(SCENES),
        help="Preset scene to generate; repeat to select multiple (default: all)",
    )
    parser.add_argument("--input-surface", type=Path, help="Generate one custom surface")
    parser.add_argument("--name", help="Output name for --input-surface")
    parser.add_argument("--resolution", type=int, help="Resolution for --input-surface")
    parser.add_argument("--E", type=float, help="Young's modulus for --input-surface")
    parser.add_argument("--nu", type=float, default=0.45)
    parser.add_argument("--density", type=float, default=1000.0)
    return parser


def selected_scenes(args: argparse.Namespace) -> dict[str, Scene]:
    if args.input_surface is not None:
        if args.scene:
            raise ValueError("--input-surface cannot be combined with --scene")
        if not args.name or args.resolution is None or args.E is None:
            raise ValueError("--input-surface requires --name, --resolution, and --E")
        if Path(args.name).name != args.name or args.name in {".", ".."}:
            raise ValueError("--name must be a single filename-safe component")
        return {
            args.name: Scene(
                args.input_surface.expanduser().resolve(),
                args.resolution,
                args.E,
                args.nu,
                args.density,
            )
        }
    if args.name or args.resolution is not None or args.E is not None:
        raise ValueError("--name, --resolution, and --E require --input-surface")
    names = args.scene or list(SCENES)
    return {name: SCENES[name] for name in names}


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    try:
        mesher = find_mesher(args.mesher, args.build_dir)
        scenes = selected_scenes(args)
        output_dir = args.output_dir.expanduser().resolve()
        if output_dir.exists():
            raise FileExistsError(
                f"Refusing to overwrite existing output directory: {output_dir}"
            )
        for name, scene in scenes.items():
            if not scene.surface.is_file():
                raise FileNotFoundError(f"Missing input surface for {name}: {scene.surface}")
            if scene.resolution <= 0:
                raise ValueError(f"Resolution for {name} must be positive")

        output_dir.parent.mkdir(parents=True, exist_ok=True)
        records = []
        with tempfile.TemporaryDirectory(
            prefix=f".{output_dir.name}-", dir=output_dir.parent
        ) as temporary:
            temporary_dir = Path(temporary)
            for name, scene in scenes.items():
                output = temporary_dir / f"{name}-cubic.veg"
                command = [
                    str(mesher),
                    "--input-mesh",
                    str(scene.surface),
                    "--resolution",
                    str(scene.resolution),
                    "--output-mesh",
                    str(output),
                    "--E",
                    str(scene.youngs_modulus),
                    "--nu",
                    str(scene.poisson_ratio),
                    "--density",
                    str(scene.density),
                ]
                print(" ".join(command), flush=True)
                subprocess.run(command, check=True)
                record = {
                    "scene": name,
                    "input_surface": str(scene.surface),
                    "input_sha256": sha256(scene.surface),
                    "output": f"{name}-cubic.veg",
                    "output_sha256": sha256(output),
                    "arguments": {
                        "input_mesh": str(scene.surface),
                        "resolution": scene.resolution,
                        "output_mesh": str(output_dir / f"{name}-cubic.veg"),
                        "E": scene.youngs_modulus,
                        "nu": scene.poisson_ratio,
                        "density": scene.density,
                    },
                    **parse_cubic_veg(output),
                }
                records.append(record)

            manifest = {
                "schema_version": 1,
                "generator": str(Path(__file__).resolve()),
                "mesher": str(mesher),
                "source_revision": revision(),
                "assets": records,
            }
            (temporary_dir / "manifest.json").write_text(
                json.dumps(manifest, indent=2) + "\n"
            )
            shutil.move(str(temporary_dir), str(output_dir))
        print(f"Wrote {len(records)} mesh(es) and manifest to {output_dir}")
        return 0
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        parser.exit(1, f"error: {error}\n")


if __name__ == "__main__":
    raise SystemExit(main())
