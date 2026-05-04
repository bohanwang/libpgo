#!/usr/bin/env python3
"""Dump edge-connected triangle components from OBJ triangle meshes.

This script only splits meshes and writes component metadata. It does not run
mesh quality checks.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class ComponentRecord:
    source: str
    output: str
    component: int
    vertices: int
    triangles: int


class UnionFind:
    def __init__(self, size: int) -> None:
        self.parent = list(range(size))
        self.size = [1] * size

    def find(self, value: int) -> int:
        while self.parent[value] != value:
            self.parent[value] = self.parent[self.parent[value]]
            value = self.parent[value]
        return value

    def union(self, a: int, b: int) -> None:
        root_a = self.find(a)
        root_b = self.find(b)
        if root_a == root_b:
            return
        if self.size[root_a] < self.size[root_b]:
            root_a, root_b = root_b, root_a
        self.parent[root_b] = root_a
        self.size[root_a] += self.size[root_b]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Dump edge-connected triangle components from OBJ meshes without running quality checks."
    )
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Input OBJ path.",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Directory where component OBJ files and components.json are written.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing component OBJ and JSON files.",
    )
    return parser.parse_args()


def component_prefix(input_path: Path) -> str:
    name = input_path.name
    if name.endswith(".veg.obj"):
        return name[: -len(".veg.obj")]
    if name.endswith(".obj"):
        return name[: -len(".obj")]
    return input_path.stem


def load_obj_triangles(path: Path) -> tuple[list[str], list[tuple[int, int, int]]]:
    vertices: list[str] = []
    triangles: list[tuple[int, int, int]] = []
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("v "):
                vertices.append(line.rstrip())
            elif line.startswith("f "):
                ids: list[int] = []
                for token in line.split()[1:]:
                    vertex_token = token.split("/")[0]
                    vertex_id = int(vertex_token)
                    if vertex_id < 0:
                        vertex_id = len(vertices) + 1 + vertex_id
                    ids.append(vertex_id - 1)
                if len(ids) != 3:
                    raise ValueError(
                        f"{path}:{line_number}: expected triangular face, got {line.rstrip()}"
                    )
                for vertex_id in ids:
                    if vertex_id < 0 or vertex_id >= len(vertices):
                        raise ValueError(
                            f"{path}:{line_number}: face references invalid vertex {vertex_id + 1}"
                        )
                triangles.append((ids[0], ids[1], ids[2]))
    if not vertices:
        raise ValueError(f"{path}: no vertices found")
    if not triangles:
        raise ValueError(f"{path}: no triangles found")
    return vertices, triangles


def edge_connected_components(triangles: list[tuple[int, int, int]]) -> list[list[int]]:
    components = UnionFind(len(triangles))
    edge_owner: dict[tuple[int, int], int] = {}
    for tri_id, tri in enumerate(triangles):
        for edge in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])):
            a, b = sorted(edge)
            previous = edge_owner.get((a, b))
            if previous is None:
                edge_owner[(a, b)] = tri_id
            else:
                components.union(previous, tri_id)

    grouped: dict[int, list[int]] = defaultdict(list)
    for tri_id in range(len(triangles)):
        grouped[components.find(tri_id)].append(tri_id)
    ret = list(grouped.values())
    ret.sort(key=lambda component: (-len(component), component[0]))
    return ret


def write_component_obj(
    path: Path,
    source: Path,
    name: str,
    vertices: list[str],
    triangles: list[tuple[int, int, int]],
    triangle_ids: list[int],
    overwrite: bool,
) -> tuple[int, int]:
    if path.exists() and not overwrite:
        raise FileExistsError(f"{path} exists; pass --overwrite to replace it")

    used_vertices: list[int] = []
    seen: set[int] = set()
    for tri_id in triangle_ids:
        for vertex_id in triangles[tri_id]:
            if vertex_id not in seen:
                seen.add(vertex_id)
                used_vertices.append(vertex_id)
    remap = {old: new + 1 for new, old in enumerate(used_vertices)}

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        handle.write(f"# component dump from {source}\n")
        handle.write(f"o {name}\n")
        for vertex_id in used_vertices:
            handle.write(vertices[vertex_id] + "\n")
        for tri_id in triangle_ids:
            tri = triangles[tri_id]
            handle.write(f"f {remap[tri[0]]} {remap[tri[1]]} {remap[tri[2]]}\n")

    return len(used_vertices), len(triangle_ids)


def dump_components_for_file(
    input_path: Path,
    output_dir: Path,
    overwrite: bool,
) -> list[ComponentRecord]:
    if not input_path.is_file():
        raise FileNotFoundError(input_path)
    vertices, triangles = load_obj_triangles(input_path)
    components = edge_connected_components(triangles)
    prefix = component_prefix(input_path)

    records: list[ComponentRecord] = []
    for component_id, triangle_ids in enumerate(components):
        output_path = output_dir / f"{prefix}_component_{component_id}.obj"
        vertex_count, triangle_count = write_component_obj(
            output_path,
            input_path,
            f"{prefix}_component_{component_id}",
            vertices,
            triangles,
            triangle_ids,
            overwrite,
        )
        records.append(
            ComponentRecord(
                source=str(input_path),
                output=str(output_path),
                component=component_id,
                vertices=vertex_count,
                triangles=triangle_count,
            )
        )

    case_summary_path = output_dir / "components.json"
    if case_summary_path.exists() and not overwrite:
        raise FileExistsError(f"{case_summary_path} exists; pass --overwrite to replace it")
    case_summary_path.parent.mkdir(parents=True, exist_ok=True)
    case_summary_path.write_text(
        json.dumps([record.__dict__ for record in records], indent=2) + "\n"
    )
    return records


def main() -> None:
    args = parse_args()
    records = dump_components_for_file(
        input_path=args.input.resolve(),
        output_dir=args.output_dir,
        overwrite=args.overwrite,
    )

    print(f"wrote {len(records)} components under {args.output_dir}")
    for record in records:
        print(
            f"{record.output}: component={record.component} "
            f"vertices={record.vertices} triangles={record.triangles}"
        )


if __name__ == "__main__":
    main()
