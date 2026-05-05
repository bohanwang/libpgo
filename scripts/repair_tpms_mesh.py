#!/usr/bin/env python3
"""
Repair a problematic OBJ triangle mesh, preserving connected components.

Designed for TPMS-like meshes where running MeshFix on the whole file may drop
smaller components. This script:
  1. loads and triangulates OBJ faces,
  2. reports boundary / non-manifold / winding issues,
  3. splits the mesh into connected components,
  4. runs PyMeshFix on each component separately,
  5. merges the repaired components,
  6. writes repaired OBJ and JSON report.

Install dependency:
  pip install pymeshfix numpy

Example:
  python repair_tpms_mesh.py tpms_schwarz_p.veg.obj --out tpms_schwarz_p.fixed.obj
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pymeshfix

Edge = Tuple[int, int]


def load_obj(path: str) -> Tuple[np.ndarray, np.ndarray]:
    vertices: List[List[float]] = []
    faces: List[List[int]] = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("v "):
                p = line.strip().split()
                if len(p) >= 4:
                    vertices.append([float(p[1]), float(p[2]), float(p[3])])
            elif line.startswith("f "):
                idx: List[int] = []
                for token in line.strip().split()[1:]:
                    s = token.split("/")[0]
                    if not s:
                        continue
                    i = int(s)
                    idx.append(i - 1 if i > 0 else len(vertices) + i)
                if len(idx) == 3:
                    faces.append(idx)
                elif len(idx) > 3:
                    for k in range(1, len(idx) - 1):
                        faces.append([idx[0], idx[k], idx[k + 1]])
    if len(vertices) == 0 or len(faces) == 0:
        raise ValueError("OBJ has no vertices or no faces")
    return np.asarray(vertices, dtype=np.float64), np.asarray(faces, dtype=np.int64)


def write_obj(path: str, vertices: np.ndarray, faces: np.ndarray, comment: str = "") -> None:
    with open(path, "w", encoding="utf-8") as f:
        if comment:
            for line in comment.splitlines():
                f.write(f"# {line}\n")
        for x, y, z in vertices:
            f.write(f"v {x:.17g} {y:.17g} {z:.17g}\n")
        for a, b, c in faces:
            f.write(f"f {int(a) + 1} {int(b) + 1} {int(c) + 1}\n")


def edge_key(a: int, b: int) -> Edge:
    return (a, b) if a < b else (b, a)


def face_edges(face: Iterable[int]) -> Iterable[Tuple[int, int]]:
    a, b, c = map(int, face)
    yield a, b
    yield b, c
    yield c, a


def build_edge_faces(faces: np.ndarray) -> Dict[Edge, List[int]]:
    edge_faces: Dict[Edge, List[int]] = defaultdict(list)
    for fi, face in enumerate(faces):
        for a, b in face_edges(face):
            edge_faces[edge_key(a, b)].append(fi)
    return edge_faces


def remove_degenerate_duplicate_faces(vertices: np.ndarray, faces: np.ndarray) -> Tuple[np.ndarray, dict]:
    keep = []
    seen = set()
    removed_degenerate = 0
    removed_duplicate = 0
    for i, face in enumerate(faces):
        a, b, c = map(int, face)
        if a == b or b == c or c == a:
            removed_degenerate += 1
            continue
        area2 = np.linalg.norm(np.cross(vertices[b] - vertices[a], vertices[c] - vertices[a]))
        if area2 <= 0.0:
            removed_degenerate += 1
            continue
        key = tuple(sorted((a, b, c)))
        if key in seen:
            removed_duplicate += 1
            continue
        seen.add(key)
        keep.append(i)
    return faces[keep].copy(), {
        "removed_degenerate_faces": removed_degenerate,
        "removed_duplicate_faces": removed_duplicate,
    }


def connected_face_components(faces: np.ndarray) -> List[List[int]]:
    edge_faces = build_edge_faces(faces)
    adj: List[List[int]] = [[] for _ in range(len(faces))]
    for fs in edge_faces.values():
        if len(fs) >= 2:
            for i in range(len(fs)):
                for j in range(i + 1, len(fs)):
                    adj[fs[i]].append(fs[j])
                    adj[fs[j]].append(fs[i])
    seen = np.zeros(len(faces), dtype=bool)
    comps: List[List[int]] = []
    for i in range(len(faces)):
        if seen[i]:
            continue
        stack = [i]
        seen[i] = True
        comp: List[int] = []
        while stack:
            cur = stack.pop()
            comp.append(cur)
            for nb in adj[cur]:
                if not seen[nb]:
                    seen[nb] = True
                    stack.append(nb)
        comps.append(comp)
    return comps


def compact_component(vertices: np.ndarray, faces: np.ndarray, face_indices: List[int]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    sub_faces = faces[face_indices]
    used_vertices = np.unique(sub_faces)
    remap = {int(old): i for i, old in enumerate(used_vertices)}
    new_faces = np.vectorize(lambda x: remap[int(x)])(sub_faces).astype(np.int64)
    return vertices[used_vertices], new_faces, used_vertices


def winding_conflict_edges(faces: np.ndarray) -> int:
    edge_dirs: Dict[Edge, List[int]] = defaultdict(list)
    for face in faces:
        for a, b in face_edges(face):
            e = edge_key(a, b)
            edge_dirs[e].append(+1 if (a, b) == e else -1)
    return sum(1 for dirs in edge_dirs.values() if len(dirs) == 2 and dirs[0] == dirs[1])


def nonmanifold_vertex_count(faces: np.ndarray) -> int:
    edge_faces = build_edge_faces(faces)
    vertex_faces: Dict[int, List[int]] = defaultdict(list)
    for fi, face in enumerate(faces):
        for v in map(int, face):
            vertex_faces[v].append(fi)
    bad = 0
    for v, inc in vertex_faces.items():
        if len(inc) <= 1:
            continue
        inc_set = set(inc)
        adj: Dict[int, List[int]] = defaultdict(list)
        for fi in inc:
            for a, b in face_edges(faces[fi]):
                if a != v and b != v:
                    continue
                fs = edge_faces[edge_key(a, b)]
                if len(fs) == 2:
                    other = fs[0] if fs[1] == fi else fs[1] if fs[0] == fi else None
                    if other is not None and other in inc_set:
                        adj[fi].append(other)
        unseen = set(inc)
        fan_count = 0
        while unseen:
            fan_count += 1
            stack = [unseen.pop()]
            while stack:
                cur = stack.pop()
                for nb in adj.get(cur, []):
                    if nb in unseen:
                        unseen.remove(nb)
                        stack.append(nb)
        if fan_count > 1:
            bad += 1
    return bad


def signed_volume(vertices: np.ndarray, faces: np.ndarray, comp: List[int] | None = None) -> float:
    tri = faces if comp is None else faces[comp]
    return float(np.einsum("ij,ij->i", vertices[tri[:, 0]], np.cross(vertices[tri[:, 1]], vertices[tri[:, 2]])).sum() / 6.0)


def stats(vertices: np.ndarray, faces: np.ndarray) -> dict:
    edge_faces = build_edge_faces(faces)
    hist = defaultdict(int)
    for fs in edge_faces.values():
        hist[len(fs)] += 1
    comps = connected_face_components(faces)
    return {
        "vertices": int(len(vertices)),
        "faces": int(len(faces)),
        "edges": int(len(edge_faces)),
        "boundary_edges": int(sum(1 for fs in edge_faces.values() if len(fs) == 1)),
        "nonmanifold_edges": int(sum(1 for fs in edge_faces.values() if len(fs) > 2)),
        "winding_conflict_edges": int(winding_conflict_edges(faces)),
        "edge_face_count_histogram": {str(k): int(v) for k, v in sorted(hist.items())},
        "connected_components": int(len(comps)),
        "component_face_counts": [int(len(c)) for c in sorted(comps, key=len, reverse=True)],
        "component_signed_volumes": [signed_volume(vertices, faces, c) for c in comps],
    }


def write_bad_edges_csv(path: str, vertices: np.ndarray, faces: np.ndarray) -> None:
    edge_faces = build_edge_faces(faces)
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["edge_v0", "edge_v1", "face_count", "v0_x", "v0_y", "v0_z", "v1_x", "v1_y", "v1_z", "faces"])
        for (a, b), fs in sorted(edge_faces.items(), key=lambda item: (-len(item[1]), item[0])):
            if len(fs) != 2:
                writer.writerow([a, b, len(fs), *vertices[a], *vertices[b], " ".join(map(str, fs))])


def repair_by_component(vertices: np.ndarray, faces: np.ndarray, remove_smallest_components: bool = False) -> Tuple[np.ndarray, np.ndarray, list]:
    comps = connected_face_components(faces)
    all_vertices = []
    all_faces = []
    offset = 0
    component_reports = []
    for ci, comp in enumerate(comps):
        cv, cf, _ = compact_component(vertices, faces, comp)
        before = stats(cv, cf)
        meshfix = pymeshfix.MeshFix(cv, cf)
        meshfix.repair(joincomp=False, remove_smallest_components=remove_smallest_components)
        rv = np.asarray(meshfix.points, dtype=np.float64)
        rf = np.asarray(meshfix.faces, dtype=np.int64)
        after = stats(rv, rf)
        all_vertices.append(rv)
        all_faces.append(rf + offset)
        offset += len(rv)
        component_reports.append({"component": ci, "before": before, "after": after})
    if all_vertices:
        return np.vstack(all_vertices), np.vstack(all_faces), component_reports
    return vertices, faces, component_reports


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="input OBJ file")
    parser.add_argument("--out", default=None, help="output OBJ file")
    parser.add_argument("--report", default=None, help="output JSON report")
    parser.add_argument("--remove-smallest-per-component", action="store_true", help="let MeshFix remove small islands inside each component")
    args = parser.parse_args()

    in_path = Path(args.input)
    out_path = Path(args.out) if args.out else in_path.with_name(in_path.stem + ".fixed.obj")
    report_path = Path(args.report) if args.report else out_path.with_suffix(".repair_report.json")
    bad_csv_path = out_path.with_suffix(".bad_edges.csv")

    vertices, faces = load_obj(str(in_path))
    report = {"input": str(in_path), "output": str(out_path)}
    report["before"] = stats(vertices, faces)

    faces, cleanup = remove_degenerate_duplicate_faces(vertices, faces)
    report["cleanup"] = cleanup
    report["after_cleanup"] = stats(vertices, faces)

    vertices, faces, component_reports = repair_by_component(vertices, faces, args.remove_smallest_per_component)
    report["component_repairs"] = component_reports
    report["after_meshfix_by_component"] = stats(vertices, faces)

    report["final"] = stats(vertices, faces)

    write_obj(str(out_path), vertices, faces, comment="repaired by repair_tpms_mesh.py using PyMeshFix per connected component")
    write_bad_edges_csv(str(bad_csv_path), vertices, faces)
    with open(report_path, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)

    print(json.dumps({
        "output": str(out_path),
        "report": str(report_path),
        "bad_edges_csv": str(bad_csv_path),
        "before": report["before"],
        "final": report["final"],
    }, indent=2))


if __name__ == "__main__":
    main()
