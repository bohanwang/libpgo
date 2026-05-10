#!/usr/bin/env python3
"""Generate no-shell FBMS assets from a completed union pipeline run."""

from __future__ import annotations

import argparse
import json
import re
import shlex
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import run_fbms_all_asset_pipeline as base  # noqa: E402


DEFAULT_CONFIG = base.REPO_ROOT / "examples" / "fbms" / "fbms_r512" / "pipeline_r512.json"
CommandSpec = base.CommandSpec
PipelineConfig = base.PipelineConfig
SourceCase = base.SourceCase
load_config = base.load_config
read_json = base.read_json
quality_passed = base.quality_passed


@dataclass(frozen=True)
class NoShellTarget:
    kind: str
    case: str
    name: str
    source_obj: Path
    sphere_obj: Path
    thickness: float
    union_raw_log: Path
    union_raw_quality: Path
    raw_obj: Path
    raw_quality: Path
    raw_log: Path
    raw_quality_log: Path
    remesh_obj: Path
    remesh_quality: Path
    remesh_log: Path
    remesh_quality_log: Path
    expected_components: int | None


def raw_resolution_tag(config: PipelineConfig) -> str:
    return base.raw_resolution_tag(config)


def fixed_thickness_tag(config: PipelineConfig) -> str:
    return base.fixed_thickness_tag(config)


def remesh_edge_tag(config: PipelineConfig) -> str:
    return base.remesh_edge_tag(config)


def tool_path(config: PipelineConfig, tool_name: str) -> Path:
    return base.tool_path(config, tool_name)


def target_key(target: NoShellTarget) -> str:
    return f"{target.kind}/{target.case}/{target.name}"


def progress_line(event: str, **fields: Any) -> str:
    parts = [f"[fbms-noshell] {event}"]
    for key, value in fields.items():
        if value is not None:
            parts.append(f"{key}={value}")
    return " ".join(parts)


def log_progress(event: str, **fields: Any) -> None:
    print(progress_line(event, **fields), flush=True)


def elapsed_text(start: float) -> str:
    return f"{time.monotonic() - start:.2f}s"


def print_command(argv: list[str]) -> None:
    print(f"+ {shlex.join(argv)}", flush=True)


def run_command(command: CommandSpec, dry_run: bool) -> int:
    start = time.monotonic()
    log_progress(
        "command",
        label=command.label,
        status="start",
        log=command.log_path,
        dry_run=str(dry_run).lower(),
    )
    print_command(command.argv)
    if dry_run:
        log_progress("command", label=command.label, status="dry-run", elapsed=elapsed_text(start))
        return 0

    if command.log_path is not None:
        command.log_path.parent.mkdir(parents=True, exist_ok=True)
    proc = subprocess.run(
        command.argv,
        cwd=REPO_ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    if command.log_path is not None:
        command.log_path.write_text(f"$ {shlex.join(command.argv)}\n{proc.stdout}", encoding="utf-8")
    status = "done" if proc.returncode == 0 else "failed"
    log_progress(
        "command",
        label=command.label,
        status=status,
        exit=proc.returncode,
        elapsed=elapsed_text(start),
        log=command.log_path,
    )
    if proc.returncode != 0:
        if command.allow_failure:
            log_progress("warning", label=command.label, reason="allowed-nonzero-exit")
        else:
            raise subprocess.CalledProcessError(proc.returncode, command.argv, output=proc.stdout)
    return proc.returncode


def parse_selected_thickness(log_text: str) -> float:
    match = re.search(r"\[volume-search\]\s+selected thickness\s*=\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)", log_text)
    if match:
        return float(match.group(1))

    tokens = shlex.split(log_text.splitlines()[0].lstrip("$ ")) if log_text.splitlines() else []
    for i, token in enumerate(tokens[:-1]):
        if token == "--fbms-thickness":
            thickness = float(tokens[i + 1])
            if thickness <= 0.0:
                raise ValueError("fixed-thickness fallback found a non-positive --fbms-thickness")
            return thickness

    raise ValueError("could not find selected thickness in union generate log")


def read_selected_thickness(log_path: Path) -> float:
    return parse_selected_thickness(log_path.read_text(encoding="utf-8"))


def noshell_fbms_target(config: PipelineConfig, case_name: str) -> NoShellTarget:
    union = base.fbms_target(config, case_name)
    resolution = raw_resolution_tag(config)
    thickness_tag = fixed_thickness_tag(config)
    edge = remesh_edge_tag(config)
    raw_dir = config.working_directory / f"raw_mesh_vdb_{resolution}_{thickness_tag}_noshell" / "fbms"
    remesh_dir = config.working_directory / f"remesh_cgal_iso_{resolution}_{thickness_tag}_noshell_{edge}" / "fbms"
    output_name = f"{case_name}_union_minus_sphere"
    return NoShellTarget(
        kind="fbms",
        case=case_name,
        name=case_name,
        source_obj=union.fbms_obj,
        sphere_obj=union.sphere_obj,
        thickness=float(config.parameters["raw"]["fbms_thickness"]),
        union_raw_log=union.raw_log,
        union_raw_quality=union.raw_quality,
        raw_obj=raw_dir / f"{output_name}.obj",
        raw_quality=raw_dir / "quality" / f"{output_name}.quality.json",
        raw_log=raw_dir / "logs" / f"{output_name}_generate.log",
        raw_quality_log=raw_dir / "logs" / f"{output_name}_raw_quality.log",
        remesh_obj=remesh_dir / f"{output_name}_remesh.obj",
        remesh_quality=remesh_dir / "quality" / f"{output_name}_remesh.quality.json",
        remesh_log=remesh_dir / "logs" / f"{output_name}_remesh.log",
        remesh_quality_log=remesh_dir / "logs" / f"{output_name}_remesh_quality.log",
        expected_components=1,
    )


def noshell_baseline_target(config: PipelineConfig, case_name: str, shape: str) -> NoShellTarget:
    union = base.baseline_target(config, case_name, shape, target_volume=1.0)
    resolution = raw_resolution_tag(config)
    edge = remesh_edge_tag(config)
    raw_dir = config.working_directory / f"raw_mesh_vdb_{resolution}_vmatch_noshell" / "baseline" / case_name
    remesh_dir = config.working_directory / f"remesh_cgal_iso_{resolution}_vmatch_noshell_{edge}" / "baseline" / case_name
    output_name = f"{shape}_union_minus_sphere"
    return NoShellTarget(
        kind="baseline",
        case=case_name,
        name=shape,
        source_obj=union.fbms_obj,
        sphere_obj=union.sphere_obj,
        thickness=read_selected_thickness(union.raw_log),
        union_raw_log=union.raw_log,
        union_raw_quality=union.raw_quality,
        raw_obj=raw_dir / f"{output_name}.obj",
        raw_quality=raw_dir / "quality" / f"{output_name}.quality.json",
        raw_log=raw_dir / "logs" / f"{output_name}_generate.log",
        raw_quality_log=raw_dir / "logs" / f"{output_name}_raw_quality.log",
        remesh_obj=remesh_dir / f"{output_name}_remesh.obj",
        remesh_quality=remesh_dir / "quality" / f"{output_name}_remesh.quality.json",
        remesh_log=remesh_dir / "logs" / f"{output_name}_remesh.log",
        remesh_quality_log=remesh_dir / "logs" / f"{output_name}_remesh_quality.log",
        expected_components=None,
    )


def raw_command(config: PipelineConfig, target: NoShellTarget) -> CommandSpec:
    raw = config.parameters["raw"]
    argv = [
        str(tool_path(config, "generate")),
        "openvdb",
        "--fbms",
        str(target.source_obj),
        "--sphere",
        str(target.sphere_obj),
        "--fbms-thickness",
        f"{target.thickness:.17g}",
        "--sphere-thickness",
        str(raw["sphere_thickness"]),
        "--resolution",
        str(raw["resolution"]),
        "--padding-ratio",
        str(raw["padding_ratio"]),
        "--output-surface",
        str(target.raw_obj),
        "--surface-mode",
        "union-minus-sphere",
    ]
    if raw.get("enable_truncating", True):
        argv.append("--enable-truncating")
    if raw.get("project_fbms_boundary_to_sphere", True):
        argv.append("--project-fbms-boundary-to-sphere")
    if raw.get("filter_small_components", True):
        argv.extend(["--filter-small-components", "--min-component-triangles", str(raw["min_component_triangles"])])
    if raw.get("keep_largest_components", -1) > 0:
        argv.extend(["--keep-largest-components", str(raw["keep_largest_components"])])
    return CommandSpec(f"raw-noshell:{target.kind}:{target.case}:{target.name}", argv, target.raw_log)


def quality_command(
    config: PipelineConfig,
    label: str,
    input_mesh: Path,
    output_json: Path,
    log_path: Path,
    check_level: str,
    invalid_policy: str,
    allow_failure: bool,
    expected_components: int | None,
    self_intersection_backend: str | None = None,
    self_intersection_triangle_limit: int | None = None,
) -> CommandSpec:
    argv = [
        str(tool_path(config, "quality")),
        "surface",
        "--check-level",
        check_level,
        "--invalid-triangles-policy",
        invalid_policy,
        "--input",
        str(input_mesh),
        "--json",
        str(output_json),
    ]
    if expected_components is not None:
        argv.extend(["--expected-components", str(expected_components)])
    if self_intersection_backend is not None:
        argv.extend(["--self-intersection-backend", self_intersection_backend])
    if self_intersection_triangle_limit is not None:
        argv.extend(["--self-intersection-triangle-limit", str(self_intersection_triangle_limit)])
    return CommandSpec(label, argv, log_path, allow_failure=allow_failure)


def _remesh_has_self_intersections(quality_path: Path) -> bool:
    if not quality_path.exists():
        return False
    q = read_json(quality_path)
    return q.get("self_intersections", 0) > 0


def remesh_geogram_command(config: PipelineConfig, target: NoShellTarget, vertex_count: int) -> CommandSpec:
    return CommandSpec(
        f"remesh-geogram-noshell:{target.kind}:{target.case}:{target.name}",
        [
            str(tool_path(config, "remesh")),
            "geogram",
            "--input-mesh",
            str(target.raw_obj),
            "--output-mesh",
            str(target.remesh_obj),
            "--target-num-vertices",
            str(vertex_count),
        ],
        target.remesh_log,
    )


def compute_remesh_edge_length_arg(raw_quality: dict[str, Any], target_edge_length: float) -> float:
    return base.compute_remesh_edge_length_arg(raw_quality, target_edge_length)


def remesh_edge_length_text(config: PipelineConfig, target: NoShellTarget, dry_run: bool) -> str:
    if target.raw_quality.exists():
        raw_quality = read_json(target.raw_quality)
        return f"{compute_remesh_edge_length_arg(raw_quality, float(config.parameters['remesh']['target_edge_length'])):.17g}"
    if dry_run:
        return f"computed:{config.parameters['remesh']['target_edge_length']}/noshell_raw_mean_edge"
    raise FileNotFoundError(target.raw_quality)


def remesh_command(config: PipelineConfig, target: NoShellTarget, dry_run: bool) -> CommandSpec:
    remesh = config.parameters["remesh"]
    return CommandSpec(
        f"remesh-noshell:{target.kind}:{target.case}:{target.name}",
        [
            str(tool_path(config, "remesh")),
            "cgal_iso",
            "--input-mesh",
            str(target.raw_obj),
            "--output-mesh",
            str(target.remesh_obj),
            "--edge-length",
            remesh_edge_length_text(config, target, dry_run=dry_run),
            "--sharp-edge-angle",
            str(remesh["sharp_edge_angle"]),
            "--iterations",
            str(remesh["iterations"]),
        ],
        target.remesh_log,
    )


def ensure_union_dependency(target: NoShellTarget, dry_run: bool) -> None:
    if dry_run:
        return
    for path in (target.source_obj, target.sphere_obj, target.union_raw_log, target.union_raw_quality):
        if not path.exists():
            raise FileNotFoundError(path)
    if not quality_passed(target.union_raw_quality):
        raise RuntimeError(f"union raw quality did not pass: {target.union_raw_quality}")


def run_noshell_target(config: PipelineConfig, target: NoShellTarget, dry_run: bool, overwrite: bool) -> None:
    ensure_union_dependency(target, dry_run=dry_run)
    raw = config.parameters["raw"]
    quality = config.parameters["quality"]

    if overwrite or not quality_passed(target.raw_quality):
        target.raw_obj.parent.mkdir(parents=True, exist_ok=True) if not dry_run else None
        run_command(raw_command(config, target), dry_run=dry_run)
        run_command(
            quality_command(
                config,
                f"raw-quality-noshell:{target.kind}:{target.case}:{target.name}",
                target.raw_obj,
                target.raw_quality,
                target.raw_quality_log,
                "raw",
                str(quality["raw_invalid_triangles_policy"]),
                allow_failure=False,
                expected_components=target.expected_components,
            ),
            dry_run=dry_run,
        )
    if not dry_run and not quality_passed(target.raw_quality):
        raise RuntimeError(f"no-shell raw quality failed for {target_key(target)}")

    if overwrite or not target.remesh_obj.exists() or not target.remesh_quality.exists():
        target.remesh_obj.parent.mkdir(parents=True, exist_ok=True) if not dry_run else None
        run_command(remesh_command(config, target, dry_run=dry_run), dry_run=dry_run)
        run_command(
            quality_command(
                config,
                f"remesh-quality-noshell:{target.kind}:{target.case}:{target.name}",
                target.remesh_obj,
                target.remesh_quality,
                target.remesh_quality_log,
                "full",
                str(quality["remesh_invalid_triangles_policy"]),
                allow_failure=True,
                expected_components=target.expected_components,
                self_intersection_backend=str(quality["remesh_self_intersection_backend"]),
                self_intersection_triangle_limit=int(quality["remesh_self_intersection_triangle_limit"]),
            ),
            dry_run=dry_run,
        )
        if not dry_run and _remesh_has_self_intersections(target.remesh_quality):
            cgal_quality = read_json(target.remesh_quality)
            vertex_count = int(cgal_quality["vertices"])
            log_progress("remesh", target=target_key(target), status="fallback-geogram",
                         reason="cgal-self-intersections", target_vertices=vertex_count)
            run_command(remesh_geogram_command(config, target, vertex_count), dry_run=dry_run)
            run_command(
                quality_command(
                    config,
                    f"remesh-quality-noshell:{target.kind}:{target.case}:{target.name}",
                    target.remesh_obj,
                    target.remesh_quality,
                    target.remesh_quality_log,
                    "full",
                    str(quality["remesh_invalid_triangles_policy"]),
                    allow_failure=True,
                    expected_components=target.expected_components,
                    self_intersection_backend=str(quality["remesh_self_intersection_backend"]),
                    self_intersection_triangle_limit=int(quality["remesh_self_intersection_triangle_limit"]),
                ),
                dry_run=dry_run,
            )
    if not dry_run and not target.remesh_quality.exists():
        raise RuntimeError(f"no-shell remesh quality JSON was not produced for {target_key(target)}")


def filter_cases(cases: tuple[SourceCase, ...], selected: list[str]) -> tuple[SourceCase, ...]:
    return base.filter_cases(cases, selected)


def filter_shapes(config: PipelineConfig, selected: list[str]) -> tuple[str, ...]:
    return base.filter_shapes(config, selected)


def requested_targets(config: PipelineConfig, cases: tuple[SourceCase, ...], shapes: tuple[str, ...], include: str) -> list[NoShellTarget]:
    targets: list[NoShellTarget] = []
    for case in cases:
        if include in ("fbms", "all"):
            targets.append(noshell_fbms_target(config, case.name))
        if include in ("baseline", "all"):
            for shape in shapes:
                targets.append(noshell_baseline_target(config, case.name, shape))
    return targets


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate no-shell FBMS assets from completed union outputs.")
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG, help="Pipeline JSON config.")
    parser.add_argument("--case", action="append", default=[], help="Restrict to one FBMS case. Can be repeated.")
    parser.add_argument("--shape", action="append", default=[], help="Restrict to one TPMS baseline shape. Can be repeated.")
    parser.add_argument("--include", choices=("fbms", "baseline", "all"), default="all", help="Which no-shell targets to generate.")
    parser.add_argument("--dry-run", action="store_true", help="Print planned commands without writing assets.")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing no-shell raw/remesh outputs.")
    return parser.parse_args()


def run_pipeline(args: argparse.Namespace) -> list[NoShellTarget]:
    start = time.monotonic()
    config = load_config(args.config)
    cases = filter_cases(base.discover_fbms_cases(config.source_asset_directory), args.case)
    shapes = filter_shapes(config, args.shape)
    targets = requested_targets(config, cases, shapes, args.include)
    log_progress(
        "pipeline",
        status="start",
        config=config.path,
        working_directory=config.working_directory,
        include=args.include,
        targets=len(targets),
        overwrite=str(args.overwrite).lower(),
        dry_run=str(args.dry_run).lower(),
    )
    for target in targets:
        log_progress("target", target=target_key(target), status="start")
        run_noshell_target(config, target, dry_run=args.dry_run, overwrite=args.overwrite)
        log_progress("target", target=target_key(target), status="done")
    log_progress("pipeline", status="done", elapsed=elapsed_text(start))
    return targets


def main() -> int:
    args = parse_args()
    try:
        run_pipeline(args)
        return 0
    except (
        FileNotFoundError,
        ValueError,
        RuntimeError,
        json.JSONDecodeError,
        subprocess.CalledProcessError,
    ) as exc:
        print(f"error: {exc}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
