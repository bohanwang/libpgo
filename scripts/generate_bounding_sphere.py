#!/usr/bin/env python3
"""Generate an OBJ bounding sphere for an OBJ mesh."""

from __future__ import annotations

import argparse
import math
import random
from itertools import combinations
from pathlib import Path


Vec3 = tuple[float, float, float]
Face = tuple[int, int, int]
Sphere = tuple[Vec3, float]


def add(a: Vec3, b: Vec3) -> Vec3:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def sub(a: Vec3, b: Vec3) -> Vec3:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def scale(v: Vec3, s: float) -> Vec3:
    return (v[0] * s, v[1] * s, v[2] * s)


def dot(a: Vec3, b: Vec3) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross(a: Vec3, b: Vec3) -> Vec3:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def squared_length(v: Vec3) -> float:
    return dot(v, v)


def distance(a: Vec3, b: Vec3) -> float:
    return math.sqrt(squared_length(sub(a, b)))


def normalize(v: Vec3) -> Vec3:
    length = math.sqrt(squared_length(v))
    if length == 0.0:
        raise ValueError("cannot normalize zero-length vector")
    return (v[0] / length, v[1] / length, v[2] / length)


def determinant_3x3(rows: tuple[Vec3, Vec3, Vec3]) -> float:
    a, b, c = rows
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def solve_3x3(rows: tuple[Vec3, Vec3, Vec3], rhs: Vec3) -> Vec3 | None:
    det = determinant_3x3(rows)
    if abs(det) < 1e-15:
        return None

    col0 = (
        (rhs[0], rows[0][1], rows[0][2]),
        (rhs[1], rows[1][1], rows[1][2]),
        (rhs[2], rows[2][1], rows[2][2]),
    )
    col1 = (
        (rows[0][0], rhs[0], rows[0][2]),
        (rows[1][0], rhs[1], rows[1][2]),
        (rows[2][0], rhs[2], rows[2][2]),
    )
    col2 = (
        (rows[0][0], rows[0][1], rhs[0]),
        (rows[1][0], rows[1][1], rhs[1]),
        (rows[2][0], rows[2][1], rhs[2]),
    )
    return (
        determinant_3x3(col0) / det,
        determinant_3x3(col1) / det,
        determinant_3x3(col2) / det,
    )


def read_obj_vertices(path: Path) -> list[Vec3]:
    vertices: list[Vec3] = []
    with path.open("r", encoding="utf-8") as fin:
        for line_number, line in enumerate(fin, start=1):
            fields = line.strip().split()
            if not fields or fields[0] != "v":
                continue
            if len(fields) < 4:
                raise ValueError(f"{path}:{line_number}: vertex line has fewer than 3 coordinates")
            try:
                vertices.append((float(fields[1]), float(fields[2]), float(fields[3])))
            except ValueError as exc:
                raise ValueError(f"{path}:{line_number}: invalid vertex coordinates") from exc

    if not vertices:
        raise ValueError(f"{path}: no OBJ vertices found")
    return vertices


def compute_aabb_center_bounding_sphere(vertices: list[Vec3], padding: float) -> tuple[Vec3, float]:
    mins = [min(v[axis] for v in vertices) for axis in range(3)]
    maxs = [max(v[axis] for v in vertices) for axis in range(3)]
    center = tuple((mins[axis] + maxs[axis]) * 0.5 for axis in range(3))

    radius = 0.0
    for vertex in vertices:
        dx = vertex[0] - center[0]
        dy = vertex[1] - center[1]
        dz = vertex[2] - center[2]
        radius = max(radius, math.sqrt(dx * dx + dy * dy + dz * dz))

    return center, radius * (1.0 + padding)


def sphere_contains_point(sphere: Sphere, point: Vec3, tolerance: float = 1e-12) -> bool:
    center, radius = sphere
    return distance(center, point) <= radius + tolerance * max(1.0, radius)


def sphere_from_2_points(a: Vec3, b: Vec3) -> Sphere:
    center = scale(add(a, b), 0.5)
    return center, distance(a, b) * 0.5


def sphere_from_3_points(a: Vec3, b: Vec3, c: Vec3) -> Sphere | None:
    ab = sub(b, a)
    ac = sub(c, a)
    normal = cross(ab, ac)
    denom = 2.0 * squared_length(normal)
    if denom < 1e-30:
        return None

    ab2 = squared_length(ab)
    ac2 = squared_length(ac)
    offset = scale(
        add(scale(cross(normal, ab), ac2), scale(cross(ac, normal), ab2)),
        1.0 / denom,
    )
    center = add(a, offset)
    return center, distance(center, a)


def sphere_from_4_points(a: Vec3, b: Vec3, c: Vec3, d: Vec3) -> Sphere | None:
    rows = (
        scale(sub(b, a), 2.0),
        scale(sub(c, a), 2.0),
        scale(sub(d, a), 2.0),
    )
    rhs = (
        squared_length(b) - squared_length(a),
        squared_length(c) - squared_length(a),
        squared_length(d) - squared_length(a),
    )
    center = solve_3x3(rows, rhs)
    if center is None:
        return None
    return center, distance(center, a)


def sphere_from_boundary(boundary: list[Vec3]) -> Sphere:
    if not boundary:
        return (0.0, 0.0, 0.0), 0.0
    if len(boundary) == 1:
        return boundary[0], 0.0

    candidates: list[Sphere] = []
    for pair in combinations(boundary, 2):
        sphere = sphere_from_2_points(pair[0], pair[1])
        if all(sphere_contains_point(sphere, point) for point in boundary):
            candidates.append(sphere)

    if len(boundary) >= 3:
        for triple in combinations(boundary, 3):
            sphere = sphere_from_3_points(triple[0], triple[1], triple[2])
            if sphere is not None and all(
                sphere_contains_point(sphere, point) for point in boundary
            ):
                candidates.append(sphere)

    if len(boundary) == 4:
        sphere = sphere_from_4_points(boundary[0], boundary[1], boundary[2], boundary[3])
        if sphere is not None and all(sphere_contains_point(sphere, point) for point in boundary):
            candidates.append(sphere)

    if not candidates:
        raise ValueError("could not construct minimal sphere from boundary points")
    return min(candidates, key=lambda sphere: sphere[1])


def compute_minimal_bounding_sphere(vertices: list[Vec3], padding: float) -> tuple[Vec3, float]:
    points = list(vertices)
    random.Random(0).shuffle(points)
    sphere: Sphere = ((0.0, 0.0, 0.0), -1.0)

    for i, point in enumerate(points):
        if sphere[1] >= 0.0 and sphere_contains_point(sphere, point):
            continue
        sphere = (point, 0.0)
        for j in range(i):
            point_j = points[j]
            if sphere_contains_point(sphere, point_j):
                continue
            sphere = sphere_from_boundary([point, point_j])
            for k in range(j):
                point_k = points[k]
                if sphere_contains_point(sphere, point_k):
                    continue
                sphere = sphere_from_boundary([point, point_j, point_k])
                for l in range(k):
                    point_l = points[l]
                    if sphere_contains_point(sphere, point_l):
                        continue
                    sphere = sphere_from_boundary([point, point_j, point_k, point_l])

    center, radius = sphere
    return center, radius * (1.0 + padding)


def create_uv_sphere_mesh(center: Vec3, radius: float, lat_segments: int, lon_segments: int) -> tuple[list[Vec3], list[Face]]:
    if lat_segments < 2:
        raise ValueError("--lat must be at least 2")
    if lon_segments < 3:
        raise ValueError("--lon must be at least 3")

    cx, cy, cz = center
    vertices: list[Vec3] = [(cx, cy, cz + radius)]

    for lat in range(1, lat_segments):
        theta = math.pi * lat / lat_segments
        z = cz + radius * math.cos(theta)
        ring_radius = radius * math.sin(theta)
        for lon in range(lon_segments):
            phi = 2.0 * math.pi * lon / lon_segments
            x = cx + ring_radius * math.cos(phi)
            y = cy + ring_radius * math.sin(phi)
            vertices.append((x, y, z))

    vertices.append((cx, cy, cz - radius))

    bottom_vertex = len(vertices)
    faces: list[Face] = []

    first_ring = 2
    for lon in range(lon_segments):
        faces.append((1, first_ring + lon, first_ring + (lon + 1) % lon_segments))

    for lat in range(lat_segments - 2):
        upper_start = 2 + lat * lon_segments
        lower_start = upper_start + lon_segments
        for lon in range(lon_segments):
            upper0 = upper_start + lon
            upper1 = upper_start + (lon + 1) % lon_segments
            lower0 = lower_start + lon
            lower1 = lower_start + (lon + 1) % lon_segments
            faces.append((upper0, lower0, lower1))
            faces.append((upper0, lower1, upper1))

    last_ring = 2 + (lat_segments - 2) * lon_segments
    for lon in range(lon_segments):
        faces.append((last_ring + lon, bottom_vertex, last_ring + (lon + 1) % lon_segments))

    return vertices, faces


def create_icosphere_mesh(center: Vec3, radius: float, subdivisions: int) -> tuple[list[Vec3], list[Face]]:
    if subdivisions < 0:
        raise ValueError("--subdivisions must be non-negative")

    phi = (1.0 + math.sqrt(5.0)) * 0.5
    vertices: list[Vec3] = [
        normalize((-1.0, phi, 0.0)),
        normalize((1.0, phi, 0.0)),
        normalize((-1.0, -phi, 0.0)),
        normalize((1.0, -phi, 0.0)),
        normalize((0.0, -1.0, phi)),
        normalize((0.0, 1.0, phi)),
        normalize((0.0, -1.0, -phi)),
        normalize((0.0, 1.0, -phi)),
        normalize((phi, 0.0, -1.0)),
        normalize((phi, 0.0, 1.0)),
        normalize((-phi, 0.0, -1.0)),
        normalize((-phi, 0.0, 1.0)),
    ]
    faces: list[Face] = [
        (0, 11, 5),
        (0, 5, 1),
        (0, 1, 7),
        (0, 7, 10),
        (0, 10, 11),
        (1, 5, 9),
        (5, 11, 4),
        (11, 10, 2),
        (10, 7, 6),
        (7, 1, 8),
        (3, 9, 4),
        (3, 4, 2),
        (3, 2, 6),
        (3, 6, 8),
        (3, 8, 9),
        (4, 9, 5),
        (2, 4, 11),
        (6, 2, 10),
        (8, 6, 7),
        (9, 8, 1),
    ]

    for _ in range(subdivisions):
        midpoint_cache: dict[tuple[int, int], int] = {}

        def midpoint_index(i: int, j: int) -> int:
            key = (min(i, j), max(i, j))
            if key in midpoint_cache:
                return midpoint_cache[key]
            midpoint = normalize(scale(add(vertices[i], vertices[j]), 0.5))
            vertices.append(midpoint)
            midpoint_cache[key] = len(vertices) - 1
            return midpoint_cache[key]

        refined_faces: list[Face] = []
        for i, j, k in faces:
            a = midpoint_index(i, j)
            b = midpoint_index(j, k)
            c = midpoint_index(k, i)
            refined_faces.extend([
                (i, a, c),
                (j, b, a),
                (k, c, b),
                (a, b, c),
            ])
        faces = refined_faces

    sphere_vertices = [add(center, scale(v, radius)) for v in vertices]
    sphere_faces = [(i + 1, j + 1, k + 1) for i, j, k in faces]
    return sphere_vertices, sphere_faces


def write_sphere_obj(
    path: Path,
    source_path: Path,
    source_vertex_count: int,
    center: Vec3,
    radius: float,
    bounds_method: str,
    sphere_method: str,
    lat_segments: int | None,
    lon_segments: int | None,
    subdivisions: int | None,
    padding: float,
    vertices: list[Vec3],
    faces: list[Face],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8") as fout:
        fout.write(f"# Bounding sphere generated for {source_path.name}\n")
        if bounds_method == "aabb":
            fout.write("# bounds_method = axis-aligned bounding-box center plus maximal vertex distance\n")
        else:
            fout.write("# bounds_method = minimal enclosing sphere\n")
        fout.write(f"# source_vertex_count = {source_vertex_count}\n")
        fout.write(f"# center = ({center[0]:.17g}, {center[1]:.17g}, {center[2]:.17g})\n")
        fout.write(f"# radius = {radius:.17g}\n")
        fout.write(f"# relative_padding = {padding:.17g}\n")
        fout.write(f"# sphere_parameterization = {sphere_method}\n")
        if sphere_method == "uv":
            fout.write(f"# sphere_resolution = latitude {lat_segments}, longitude {lon_segments}\n")
        else:
            fout.write(f"# sphere_resolution = icosphere subdivisions {subdivisions}\n")
        fout.write(f"# sphere_vertices = {len(vertices)}\n")
        fout.write(f"# sphere_faces = {len(faces)}\n")
        object_name = "minimal_bounding_sphere" if bounds_method == "minimal" else "bounding_sphere"
        fout.write(f"o {object_name}\n")

        for x, y, z in vertices:
            fout.write(f"v {x:.17g} {y:.17g} {z:.17g}\n")
        for i, j, k in faces:
            fout.write(f"f {i} {j} {k}\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate an OBJ sphere enclosing the vertices of an OBJ mesh.",
    )
    parser.add_argument("--input", required=True, type=Path, help="Input OBJ mesh path")
    parser.add_argument("--output", required=True, type=Path, help="Output bounding sphere OBJ path")
    parser.add_argument(
        "--bounds-method",
        choices=("aabb", "minimal"),
        default="aabb",
        help="Bounding sphere computation. Default: aabb, matching the original g0_b8 asset.",
    )
    parser.add_argument(
        "--method",
        choices=("icosphere", "uv"),
        default="icosphere",
        help="Sphere triangulation. Default: icosphere for more uniform triangles.",
    )
    parser.add_argument("--subdivisions", default=5, type=int, help="Icosphere subdivision count, default: 5")
    parser.add_argument("--lat", default=64, type=int, help="Latitude segments, default: 64")
    parser.add_argument("--lon", default=128, type=int, help="Longitude segments, default: 128")
    parser.add_argument(
        "--padding",
        default=1e-9,
        type=float,
        help="Relative radius padding. Default: 1e-9, matching the existing g0_b8 sphere.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    try:
        if args.padding < 0.0:
            raise ValueError("--padding must be non-negative")
        source_vertices = read_obj_vertices(args.input)
        if args.bounds_method == "minimal":
            center, radius = compute_minimal_bounding_sphere(source_vertices, args.padding)
        else:
            center, radius = compute_aabb_center_bounding_sphere(source_vertices, args.padding)
        if args.method == "uv":
            sphere_vertices, sphere_faces = create_uv_sphere_mesh(center, radius, args.lat, args.lon)
        else:
            sphere_vertices, sphere_faces = create_icosphere_mesh(center, radius, args.subdivisions)
        write_sphere_obj(
            args.output,
            args.input,
            len(source_vertices),
            center,
            radius,
            args.bounds_method,
            args.method,
            args.lat if args.method == "uv" else None,
            args.lon if args.method == "uv" else None,
            args.subdivisions if args.method == "icosphere" else None,
            args.padding,
            sphere_vertices,
            sphere_faces,
        )
    except OSError as exc:
        print(f"error: {exc}")
        return 1
    except ValueError as exc:
        print(f"error: {exc}")
        return 1

    print(
        f"Saved bounding sphere to {args.output} "
        f"(bounds_method={args.bounds_method}, "
        f"center=({center[0]:.17g}, {center[1]:.17g}, {center[2]:.17g}), "
        f"radius={radius:.17g}, vertices={len(sphere_vertices)}, faces={len(sphere_faces)})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
