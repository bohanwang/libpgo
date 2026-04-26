#!/usr/bin/env python3
"""Export FBMS tet stress frames to ParaView-readable VTU/PVD files."""

from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import struct
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from xml.sax.saxutils import escape


@dataclass(frozen=True)
class ExportConfig:
    veg: Path
    states: Path
    stress: Path
    output: Path
    frame_start: int
    frame_end: int


def resolve_config_path(config_dir: Path, value: Any, key: str) -> Path:
    if not isinstance(value, str) or not value:
        raise ValueError(f"config field '{key}' must be a non-empty string")
    path = Path(value)
    if not path.is_absolute():
        path = config_dir / path
    return Path(os.path.normpath(str(path)))


def load_config(config_path: Path | str) -> ExportConfig:
    config_path = Path(config_path).expanduser().absolute()
    with config_path.open("r", encoding="utf-8") as f:
        data = json.load(f)

    if not isinstance(data, dict):
        raise ValueError("export config must be a JSON object")

    required_paths = ("veg", "states", "stress", "output")
    missing = [key for key in required_paths if key not in data]
    if missing:
        raise ValueError(f"missing required config fields: {', '.join(missing)}")

    if "frame-end" not in data:
        raise ValueError("missing required config field: frame-end")
    frame_start = int(data.get("frame-start", 0))
    frame_end = int(data["frame-end"])
    if frame_start < 0 or frame_end <= frame_start:
        raise ValueError("frame range must satisfy 0 <= frame-start < frame-end")

    config_dir = config_path.parent
    return ExportConfig(
        veg=resolve_config_path(config_dir, data["veg"], "veg"),
        states=resolve_config_path(config_dir, data["states"], "states"),
        stress=resolve_config_path(config_dir, data["stress"], "stress"),
        output=resolve_config_path(config_dir, data["output"], "output"),
        frame_start=frame_start,
        frame_end=frame_end,
    )


def _next_content_line(lines: list[str], index: int) -> tuple[str, int]:
    while index < len(lines):
        line = lines[index].strip()
        index += 1
        if line and not line.startswith("#"):
            return line, index
    raise ValueError("unexpected end of .veg file")


def parse_veg(path: Path | str) -> tuple[list[tuple[float, float, float]], list[tuple[int, int, int, int]]]:
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    index = 0

    while index < len(lines) and lines[index].strip().upper() != "*VERTICES":
        index += 1
    if index >= len(lines):
        raise ValueError(f"{path}: missing *VERTICES section")
    index += 1

    header, index = _next_content_line(lines, index)
    header_parts = header.split()
    if len(header_parts) < 2:
        raise ValueError(f"{path}: invalid vertex header")
    vertex_count = int(header_parts[0])

    vertices: list[tuple[float, float, float]] = []
    for _ in range(vertex_count):
        line, index = _next_content_line(lines, index)
        parts = line.split()
        if len(parts) < 4:
            raise ValueError(f"{path}: invalid vertex row: {line}")
        vertices.append((float(parts[1]), float(parts[2]), float(parts[3])))

    while index < len(lines) and lines[index].strip().upper() != "*ELEMENTS":
        index += 1
    if index >= len(lines):
        raise ValueError(f"{path}: missing *ELEMENTS section")
    index += 1

    element_type, index = _next_content_line(lines, index)
    if element_type.upper() != "TET":
        raise ValueError(f"{path}: only TET .veg elements are supported")

    header, index = _next_content_line(lines, index)
    header_parts = header.split()
    if len(header_parts) < 2:
        raise ValueError(f"{path}: invalid element header")
    tet_count = int(header_parts[0])

    tets: list[tuple[int, int, int, int]] = []
    for _ in range(tet_count):
        line, index = _next_content_line(lines, index)
        parts = line.split()
        if len(parts) < 5:
            raise ValueError(f"{path}: invalid tet row: {line}")
        tet = tuple(int(value) - 1 for value in parts[1:5])
        if any(vertex_index < 0 or vertex_index >= vertex_count for vertex_index in tet):
            raise ValueError(f"{path}: tet references a vertex outside 1..{vertex_count}")
        tets.append(tet)

    return vertices, tets


def read_state_displacement(path: Path | str, vertex_count: int) -> list[tuple[float, float, float]]:
    path = Path(path)
    data = path.read_bytes()
    if len(data) < 12:
        raise ValueError(f"{path}: state file is too small")

    row_count, column_count, entry_size = struct.unpack_from("<iii", data, 0)
    expected_rows = vertex_count * 3
    if row_count != expected_rows:
        raise ValueError(f"{path}: expected {expected_rows} rows, found {row_count}")
    if column_count < 1:
        raise ValueError(f"{path}: expected at least one matrix column")
    if entry_size != 8:
        raise ValueError(f"{path}: only Float64 state matrices are supported")

    value_count = row_count * column_count
    expected_bytes = 12 + value_count * entry_size
    if len(data) != expected_bytes:
        raise ValueError(f"{path}: expected {expected_bytes} bytes, found {len(data)}")

    flat_displacement = struct.unpack_from(f"<{row_count}d", data, 12)
    return [
        (
            flat_displacement[vertex_index * 3],
            flat_displacement[vertex_index * 3 + 1],
            flat_displacement[vertex_index * 3 + 2],
        )
        for vertex_index in range(vertex_count)
    ]


def read_stress(path: Path | str, tet_count: int) -> tuple[float | int, list[float]]:
    path = Path(path)
    with path.open("r", encoding="utf-8") as f:
        data = json.load(f)

    if not isinstance(data, dict):
        raise ValueError(f"{path}: stress file must be a JSON object")
    if data.get("location") != "tet_element":
        raise ValueError(f"{path}: expected stress location 'tet_element'")
    values = data.get("values")
    if not isinstance(values, list):
        raise ValueError(f"{path}: stress field 'values' must be a list")
    if len(values) != tet_count:
        raise ValueError(f"{path}: expected {tet_count} stress values, found {len(values)}")

    timestep = data.get("time", data.get("frame", 0))
    return timestep, [float(value) for value in values]


def format_numbers(values: list[float] | tuple[float, ...]) -> str:
    return " ".join(f"{value:.17g}" for value in values)


def clamped_percentile(values: list[float], percentile: float) -> list[float]:
    if not values:
        return []
    sorted_values = sorted(values)
    index = int(percentile * (len(sorted_values) - 1))
    threshold = sorted_values[index]
    return [min(value, threshold) for value in values]


def stress_cell_data_arrays(stress_values: list[float]) -> list[tuple[str, list[float]]]:
    log_floor = 1.0e-30
    return [
        ("von_mises", stress_values),
        ("von_mises_log10", [math.log10(max(value, log_floor)) for value in stress_values]),
        ("von_mises_clamped_99", clamped_percentile(stress_values, 0.99)),
    ]


def cell_data_array_xml(name: str, values: list[float]) -> str:
    return f"""        <DataArray type="Float64" Name="{escape(name)}" format="ascii">
          {format_numbers(tuple(values))}
        </DataArray>"""


def write_vtu(
    path: Path,
    points: list[tuple[float, float, float]],
    tets: list[tuple[int, int, int, int]],
    stress_values: list[float],
) -> None:
    point_text = "\n          ".join(format_numbers(point) for point in points)
    connectivity = " ".join(str(index) for tet in tets for index in tet)
    offsets = " ".join(str((index + 1) * 4) for index in range(len(tets)))
    types = " ".join("10" for _ in tets)
    cell_data_text = "\n".join(cell_data_array_xml(name, values) for name, values in stress_cell_data_arrays(stress_values))

    path.write_text(
        f"""<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">
  <UnstructuredGrid>
    <Piece NumberOfPoints="{len(points)}" NumberOfCells="{len(tets)}">
      <Points>
        <DataArray type="Float64" NumberOfComponents="3" format="ascii">
          {point_text}
        </DataArray>
      </Points>
      <Cells>
        <DataArray type="Int32" Name="connectivity" format="ascii">
          {connectivity}
        </DataArray>
        <DataArray type="Int32" Name="offsets" format="ascii">
          {offsets}
        </DataArray>
        <DataArray type="UInt8" Name="types" format="ascii">
          {types}
        </DataArray>
      </Cells>
      <CellData Scalars="von_mises">
{cell_data_text}
      </CellData>
    </Piece>
  </UnstructuredGrid>
</VTKFile>
""",
        encoding="utf-8",
    )


def write_pvd(path: Path, entries: list[tuple[float | int, Path]]) -> None:
    dataset_lines = []
    for timestep, frame_path in entries:
        dataset_lines.append(
            f'    <DataSet timestep="{escape(str(timestep))}" group="" part="0" file="{escape(frame_path.name)}"/>'
        )
    path.write_text(
        f"""<?xml version="1.0"?>
<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">
  <Collection>
{chr(10).join(dataset_lines)}
  </Collection>
</VTKFile>
""",
        encoding="utf-8",
    )


def prepare_output_dir(output: Path, overwrite: bool) -> None:
    if output.exists():
        if any(output.iterdir()):
            if not overwrite:
                raise FileExistsError(f"{output} already contains files; pass --overwrite to replace it")
            shutil.rmtree(output)
        elif not output.is_dir():
            raise FileExistsError(f"{output} exists and is not a directory")
    output.mkdir(parents=True, exist_ok=True)


def export(config: ExportConfig, overwrite: bool = False) -> list[Path]:
    vertices, tets = parse_veg(config.veg)
    prepare_output_dir(config.output, overwrite)

    written_frames: list[Path] = []
    pvd_entries: list[tuple[float | int, Path]] = []
    for frame in range(config.frame_start, config.frame_end):
        state_path = config.states / f"deform{frame:04d}.u"
        stress_path = config.stress / f"von_mises{frame:04d}.json"
        displacements = read_state_displacement(state_path, len(vertices))
        timestep, stress_values = read_stress(stress_path, len(tets))
        points = [
            (
                vertex[0] + displacement[0],
                vertex[1] + displacement[1],
                vertex[2] + displacement[2],
            )
            for vertex, displacement in zip(vertices, displacements)
        ]
        frame_path = config.output / f"frame{frame:04d}.vtu"
        write_vtu(frame_path, points, tets, stress_values)
        written_frames.append(frame_path)
        pvd_entries.append((timestep, frame_path))

    pvd_path = config.output / "series.pvd"
    write_pvd(pvd_path, pvd_entries)
    return written_frames + [pvd_path]


def export_from_config(config_path: Path | str, overwrite: bool = False) -> list[Path]:
    return export(load_config(config_path), overwrite=overwrite)


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path, help="JSON export config path")
    parser.add_argument("--overwrite", action="store_true", help="replace existing output files")
    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    try:
        written = export_from_config(args.config, overwrite=args.overwrite)
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    print(f"wrote {len(written)} files to {written[-1].parent}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
