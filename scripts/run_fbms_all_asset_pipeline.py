#!/usr/bin/env python3
"""Run the JSON-driven FBMS-all r512 asset pipeline.

The pipeline starts from a source asset directory containing FBMS shell OBJ
files named g0_b*.obj or g0b*.obj. It prepares a working asset folder, generates
FBMS assets at fixed thickness 0.02, then generates volume-matched TPMS
baselines using each FBMS raw union volume as the budget.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = REPO_ROOT / "examples" / "fbms_all" / "pipeline_r512.json"
STAGES = ("prepare", "raw", "remesh", "tet", "validate", "summary")
README_RESULTS_START = "<!-- r512-pipeline-results:start -->"
README_RESULTS_END = "<!-- r512-pipeline-results:end -->"

TPMS_SHAPES = (
    "tpms_schwarz_p",
    "tpms_schwarz_d",
    "tpms_gyroid",
    "tpms_iwp",
    "tpms_neovius",
)

DEFAULT_PARAMETERS: dict[str, Any] = {
    "asset": {
        "tpms_cells": 1.0,
        "tpms_resolution": 256,
        "tpms_sphere_subdivisions": 4,
        "bounding_sphere_bounds_method": "aabb",
        "bounding_sphere_method": "icosphere",
        "bounding_sphere_subdivisions": 5,
        "bounding_sphere_padding": 1e-9,
    },
    "raw": {
        "resolution": 512,
        "fbms_thickness": 0.02,
        "sphere_thickness": 0.03,
        "padding_ratio": 0.08,
        "enable_truncating": True,
        "project_fbms_boundary_to_sphere": True,
        "filter_small_components": True,
        "min_component_triangles": 100,
    },
    "remesh": {
        "target_edge_length": 0.02,
        "sharp_edge_angle": 180,
        "iterations": 10,
    },
    "tetwild": {
        "lr": 0.008,
        "epsr": 0.001,
    },
    "material": {
        "density": 1000,
        "young_modulus": 10000000,
        "poisson_ratio": 0.45,
    },
    "quality": {
        "expected_components": 3,
        "raw_invalid_triangles_policy": "warn",
        "remesh_invalid_triangles_policy": "warn",
        "remesh_self_intersection_backend": "exact-count",
        "remesh_self_intersection_triangle_limit": 900000,
        "final_invalid_triangles_policy": "fail",
        "final_self_intersection_backend": "exact-count",
        "final_self_intersection_triangle_limit": 900000,
    },
}


@dataclass(frozen=True)
class CommandSpec:
    label: str
    argv: list[str]
    log_path: Path | None = None
    allow_failure: bool = False


@dataclass(frozen=True)
class PrepareAction:
    label: str
    argv: list[str]
    output: Path


@dataclass(frozen=True)
class SourceCase:
    name: str
    source_obj: Path


@dataclass(frozen=True)
class PipelineConfig:
    path: Path
    source_asset_directory: Path
    working_directory: Path
    build_dir: Path
    parameters: dict[str, Any]
    tpms_shapes: tuple[str, ...]


@dataclass(frozen=True)
class PipelineTarget:
    kind: str
    case: str
    name: str
    fbms_obj: Path
    sphere_obj: Path
    raw_thickness: float
    raw_obj: Path
    raw_quality: Path
    raw_log: Path
    raw_quality_log: Path
    remesh_obj: Path
    remesh_quality: Path
    remesh_log: Path
    remesh_quality_log: Path
    tet_dir: Path
    component_dir: Path

    @property
    def tet_config(self) -> Path:
        return self.tet_dir / "tetmesher.json"

    @property
    def veg(self) -> Path:
        return self.tet_dir / f"{self.name}.veg"

    @property
    def veg_obj(self) -> Path:
        return self.tet_dir / f"{self.name}.veg.obj"

    @property
    def veg_quality(self) -> Path:
        return self.tet_dir / f"{self.name}.veg.quality.json"

    @property
    def repaired_veg_obj(self) -> Path:
        return self.tet_dir / f"{self.name}.veg.repaired.obj"

    @property
    def repaired_veg_quality(self) -> Path:
        return self.tet_dir / f"{self.name}.veg.repaired.quality.json"

    @property
    def repair_report(self) -> Path:
        return self.tet_dir / f"{self.name}.veg.repair_report.json"

    @property
    def veg_info(self) -> Path:
        return self.tet_dir / f"{self.name}.veg.info.txt"

    @property
    def tet_log(self) -> Path:
        return self.tet_dir / f"{self.name}_tetwild.log"

    @property
    def boundary_quality_log(self) -> Path:
        return self.tet_dir / f"{self.name}_boundary_quality.log"

    @property
    def boundary_repair_log(self) -> Path:
        return self.tet_dir / f"{self.name}_boundary_repair.log"

    @property
    def boundary_repaired_quality_log(self) -> Path:
        return self.tet_dir / f"{self.name}_boundary_repaired_quality.log"

    @property
    def components_log(self) -> Path:
        return self.tet_dir / f"{self.name}_components.log"


def require_mapping(value: Any, context: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ValueError(f"{context} must be an object")
    return value


def require_string(value: Any, context: str) -> str:
    if not isinstance(value, str) or not value:
        raise ValueError(f"{context} must be a non-empty string")
    return value


def repo_path(path_text: str | Path) -> Path:
    path = Path(path_text)
    return path if path.is_absolute() else REPO_ROOT / path


def deep_merge(defaults: dict[str, Any], overrides: dict[str, Any]) -> dict[str, Any]:
    merged: dict[str, Any] = {}
    for key, value in defaults.items():
        if isinstance(value, dict):
            merged[key] = deep_merge(value, {})
        else:
            merged[key] = value
    for key, value in overrides.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = deep_merge(merged[key], value)
        else:
            merged[key] = value
    return merged


def load_config(config_path: Path) -> PipelineConfig:
    config_path = repo_path(config_path)
    with config_path.open("r", encoding="utf-8") as fin:
        config = require_mapping(json.load(fin), "config")

    shapes_value = config.get("tpms_shapes", list(TPMS_SHAPES))
    if not isinstance(shapes_value, list) or not shapes_value:
        raise ValueError("tpms_shapes must be a non-empty array")
    shapes = tuple(require_string(value, "tpms_shapes[]") for value in shapes_value)

    return PipelineConfig(
        path=config_path,
        source_asset_directory=repo_path(
            require_string(config.get("source_asset_directory"), "source_asset_directory")
        ),
        working_directory=repo_path(require_string(config.get("working_directory"), "working_directory")),
        build_dir=repo_path(require_string(config.get("build_dir", "build/base_no_mkl"), "build_dir")),
        parameters=deep_merge(DEFAULT_PARAMETERS, require_mapping(config.get("parameters", {}), "parameters")),
        tpms_shapes=shapes,
    )


def is_fbms_case_obj(path: Path) -> bool:
    if path.suffix.lower() != ".obj":
        return False
    if path.name.endswith("_bounding_sphere.obj"):
        return False
    return re.fullmatch(r"g0_?b.+", path.stem) is not None


def case_sort_key(case: SourceCase) -> tuple[int, str]:
    match = re.search(r"b(\d+)", case.name)
    return (int(match.group(1)) if match else 10**9, case.name)


def discover_fbms_cases(source_asset_directory: Path) -> tuple[SourceCase, ...]:
    if not source_asset_directory.is_dir():
        raise FileNotFoundError(source_asset_directory)

    by_name: dict[str, Path] = {}
    for path in sorted(source_asset_directory.rglob("*.obj")):
        if not is_fbms_case_obj(path):
            continue
        existing = by_name.get(path.stem)
        if existing is not None and existing.resolve() != path.resolve():
            raise ValueError(f"duplicate FBMS case {path.stem}: {existing} and {path}")
        by_name[path.stem] = path

    if not by_name:
        raise ValueError(f"no g0_b*/g0b* OBJ files found under {source_asset_directory}")
    return tuple(sorted((SourceCase(name, path) for name, path in by_name.items()), key=case_sort_key))


def find_source_asset(source_asset_directory: Path, filename: str) -> Path | None:
    matches = sorted(source_asset_directory.rglob(filename))
    return matches[0] if matches else None


def tool_path(config: PipelineConfig, tool_name: str) -> Path:
    if tool_name == "generate":
        return config.build_dir / "bin" / "generateFBMSUnionSurface"
    if tool_name == "quality":
        return config.build_dir / "bin" / "meshQualityCheck"
    if tool_name == "remesh":
        return config.build_dir / "bin" / "remeshSurface"
    if tool_name == "tetmesher":
        return config.build_dir / "bin" / "tetMesher"
    if tool_name == "volume_info":
        return config.build_dir / "bin" / "volumetricMeshInfo"
    raise ValueError(f"unknown tool name: {tool_name}")


def ensure_tool(path: Path, dry_run: bool) -> None:
    if not dry_run and not path.exists():
        raise FileNotFoundError(path)


def progress_line(event: str, **fields: Any) -> str:
    parts = [f"[fbms-all] {event}"]
    for key, value in fields.items():
        if value is not None:
            parts.append(f"{key}={value}")
    return " ".join(parts)


def log_progress(event: str, **fields: Any) -> None:
    print(progress_line(event, **fields), flush=True)


def target_key(target: PipelineTarget) -> str:
    return f"{target.kind}/{target.case}/{target.name}"


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
        command.log_path.write_text(
            f"$ {shlex.join(command.argv)}\n{proc.stdout}",
            encoding="utf-8",
        )
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


def run_command_capture(argv: list[str], output_path: Path, dry_run: bool) -> None:
    start = time.monotonic()
    label = f"capture:{output_path.name}"
    log_progress("command", label=label, status="start", output=output_path, dry_run=str(dry_run).lower())
    print_command(argv + [">", str(output_path)])
    if dry_run:
        log_progress("command", label=label, status="dry-run", elapsed=elapsed_text(start))
        return
    output_path.parent.mkdir(parents=True, exist_ok=True)
    proc = subprocess.run(
        argv,
        cwd=REPO_ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    output_path.write_text(proc.stdout, encoding="utf-8")
    log_progress(
        "command",
        label=label,
        status="done" if proc.returncode == 0 else "failed",
        exit=proc.returncode,
        elapsed=elapsed_text(start),
        output=output_path,
    )
    if proc.returncode != 0:
        raise subprocess.CalledProcessError(proc.returncode, argv, output=proc.stdout)


def copy_file(src: Path, dst: Path, dry_run: bool, overwrite: bool) -> PrepareAction | None:
    action = PrepareAction("copy", ["copy", str(src), str(dst)], dst)
    if dst.exists() and not overwrite:
        log_progress("prepare", action="copy", status="skip", output=dst, reason="exists")
        return None
    if src.exists() and dst.exists() and src.resolve() == dst.resolve():
        log_progress("prepare", action="copy", status="skip", output=dst, reason="same-file")
        return None
    log_progress("prepare", action="copy", status="start", source=src, output=dst)
    print(f"+ copy {src} {dst}", flush=True)
    if not dry_run:
        dst.parent.mkdir(parents=True, exist_ok=True)
        if overwrite and dst.exists():
            dst.unlink()
        shutil.copy2(src, dst)
    log_progress("prepare", action="copy", status="dry-run" if dry_run else "done", output=dst)
    return action


def expected_tpms_paths(asset_baseline_dir: Path, shape: str) -> tuple[Path, Path]:
    return (
        asset_baseline_dir / f"{shape}_fbms.obj",
        asset_baseline_dir / f"{shape}_fbms_bounding_sphere.obj",
    )


def prepare_assets(
    config: PipelineConfig,
    cases: tuple[SourceCase, ...],
    dry_run: bool,
    overwrite: bool,
) -> list[PrepareAction | CommandSpec]:
    start = time.monotonic()
    log_progress(
        "stage",
        stage="prepare",
        status="start",
        cases=len(cases),
        shapes=len(config.tpms_shapes),
        overwrite=str(overwrite).lower(),
    )
    params = config.parameters["asset"]
    actions: list[PrepareAction | CommandSpec] = []
    planned_outputs: set[Path] = set()
    asset_fbms_dir = config.working_directory / "asset" / "fbms"
    asset_baseline_dir = config.working_directory / "asset" / "baseline"
    if not dry_run:
        asset_fbms_dir.mkdir(parents=True, exist_ok=True)
        asset_baseline_dir.mkdir(parents=True, exist_ok=True)

    for case in cases:
        dst_obj = asset_fbms_dir / f"{case.name}.obj"
        copy_action = copy_file(case.source_obj, dst_obj, dry_run=dry_run, overwrite=overwrite)
        if copy_action is not None:
            actions.append(copy_action)
            planned_outputs.add(copy_action.output)

        dst_sphere = asset_fbms_dir / f"{case.name}_bounding_sphere.obj"
        source_sphere = case.source_obj.with_name(f"{case.name}_bounding_sphere.obj")
        if source_sphere.exists():
            sphere_action = copy_file(source_sphere, dst_sphere, dry_run=dry_run, overwrite=overwrite)
            if sphere_action is not None:
                actions.append(sphere_action)
                planned_outputs.add(sphere_action.output)
        elif overwrite or not dst_sphere.exists():
            cmd = CommandSpec(
                label=f"prepare-sphere:{case.name}",
                argv=[
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "generate_bounding_sphere.py"),
                    "--input",
                    str(dst_obj),
                    "--output",
                    str(dst_sphere),
                    "--bounds-method",
                    str(params["bounding_sphere_bounds_method"]),
                    "--method",
                    str(params["bounding_sphere_method"]),
                    "--subdivisions",
                    str(params["bounding_sphere_subdivisions"]),
                    "--padding",
                    str(params["bounding_sphere_padding"]),
                ],
            )
            actions.append(cmd)
            run_command(cmd, dry_run=dry_run)

    missing_tpms = []
    def prepared_or_existing(path: Path) -> bool:
        return path.exists() or path in planned_outputs

    for shape in config.tpms_shapes:
        shell, sphere = expected_tpms_paths(asset_baseline_dir, shape)
        source_shell = find_source_asset(config.source_asset_directory, shell.name)
        source_sphere = find_source_asset(config.source_asset_directory, sphere.name)
        if source_shell is not None:
            copy_action = copy_file(source_shell, shell, dry_run=dry_run, overwrite=overwrite)
            if copy_action is not None:
                actions.append(copy_action)
                planned_outputs.add(copy_action.output)
        if source_sphere is not None:
            copy_action = copy_file(source_sphere, sphere, dry_run=dry_run, overwrite=overwrite)
            if copy_action is not None:
                actions.append(copy_action)
                planned_outputs.add(copy_action.output)
        if not prepared_or_existing(shell) or not prepared_or_existing(sphere):
            missing_tpms.append(shape)

    if missing_tpms:
        log_progress("prepare", action="generate-tpms", status="needed", shapes=",".join(missing_tpms))
        cmd = CommandSpec(
            label="prepare-tpms",
            argv=[
                sys.executable,
                str(REPO_ROOT / "scripts" / "generate_tpms_unit_ball.py"),
                "--cells",
                str(params["tpms_cells"]),
                "--resolution",
                str(params["tpms_resolution"]),
                "--sphere-subdivisions",
                str(params["tpms_sphere_subdivisions"]),
                "--out",
                str(asset_baseline_dir),
                "--flat-output",
            ],
        )
        actions.append(cmd)
        run_command(cmd, dry_run=dry_run)
    else:
        log_progress("prepare", action="generate-tpms", status="skip", reason="baseline-assets-present")

    log_progress("stage", stage="prepare", status="done", actions=len(actions), elapsed=elapsed_text(start))
    return actions


def stage_enabled(from_stage: str, stage: str) -> bool:
    return STAGES.index(stage) >= STAGES.index(from_stage)


def asset_fbms_obj(config: PipelineConfig, case_name: str) -> Path:
    return config.working_directory / "asset" / "fbms" / f"{case_name}.obj"


def asset_fbms_sphere(config: PipelineConfig, case_name: str) -> Path:
    return config.working_directory / "asset" / "fbms" / f"{case_name}_bounding_sphere.obj"


def asset_tpms_obj(config: PipelineConfig, shape: str) -> Path:
    return config.working_directory / "asset" / "baseline" / f"{shape}_fbms.obj"


def asset_tpms_sphere(config: PipelineConfig, shape: str) -> Path:
    return config.working_directory / "asset" / "baseline" / f"{shape}_fbms_bounding_sphere.obj"


def raw_resolution_tag(config: PipelineConfig) -> str:
    return f"r{int(config.parameters['raw']['resolution'])}"


def fixed_thickness_tag(config: PipelineConfig) -> str:
    return f"t{int(round(abs(float(config.parameters['raw']['fbms_thickness'])) * 100)):03d}"


def remesh_edge_tag(config: PipelineConfig) -> str:
    return f"e{int(round(float(config.parameters['remesh']['target_edge_length']) * 100)):03d}"


def tetwild_lr_tag(config: PipelineConfig) -> str:
    return f"lr{int(round(float(config.parameters['tetwild']['lr']) * 1000)):04d}"


def summary_path(config: PipelineConfig) -> Path:
    return config.working_directory / f"pipeline_{raw_resolution_tag(config)}_summary.json"


def fbms_target(config: PipelineConfig, case_name: str) -> PipelineTarget:
    resolution = raw_resolution_tag(config)
    thickness = fixed_thickness_tag(config)
    edge = remesh_edge_tag(config)
    lr = tetwild_lr_tag(config)
    raw_dir = config.working_directory / f"raw_mesh_vdb_{resolution}_{thickness}" / "fbms"
    remesh_dir = config.working_directory / f"remesh_cgal_iso_{resolution}_{thickness}_{edge}" / "fbms"
    tet_dir = config.working_directory / f"tetmesh_tetwild_{resolution}_{thickness}_{lr}" / "fbms" / case_name
    component_dir = config.working_directory / f"tetmesh_tetwild_{resolution}_{thickness}_{lr}_components" / "fbms" / case_name
    raw = config.parameters["raw"]
    return PipelineTarget(
        kind="fbms",
        case=case_name,
        name=case_name,
        fbms_obj=asset_fbms_obj(config, case_name),
        sphere_obj=asset_fbms_sphere(config, case_name),
        raw_thickness=float(raw["fbms_thickness"]),
        raw_obj=raw_dir / f"{case_name}_union.obj",
        raw_quality=raw_dir / "quality" / f"{case_name}_union.quality.json",
        raw_log=raw_dir / "logs" / f"{case_name}_generate.log",
        raw_quality_log=raw_dir / "logs" / f"{case_name}_raw_quality.log",
        remesh_obj=remesh_dir / f"{case_name}_remesh.obj",
        remesh_quality=remesh_dir / "quality" / f"{case_name}_remesh.quality.json",
        remesh_log=remesh_dir / "logs" / f"{case_name}_remesh.log",
        remesh_quality_log=remesh_dir / "logs" / f"{case_name}_remesh_quality.log",
        tet_dir=tet_dir,
        component_dir=component_dir,
    )


def baseline_target(config: PipelineConfig, case_name: str, shape: str, target_volume: float) -> PipelineTarget:
    resolution = raw_resolution_tag(config)
    edge = remesh_edge_tag(config)
    lr = tetwild_lr_tag(config)
    raw_dir = config.working_directory / f"raw_mesh_vdb_{resolution}_vmatch" / "baseline" / case_name
    remesh_dir = config.working_directory / f"remesh_cgal_iso_{resolution}_vmatch_{edge}" / "baseline" / case_name
    tet_dir = (
        config.working_directory
        / f"tetmesh_tetwild_{resolution}_vmatch_{lr}"
        / "baseline"
        / case_name
        / shape
    )
    component_dir = (
        config.working_directory
        / f"tetmesh_tetwild_{resolution}_vmatch_{lr}_components"
        / "baseline"
        / case_name
        / shape
    )
    return PipelineTarget(
        kind="baseline",
        case=case_name,
        name=shape,
        fbms_obj=asset_tpms_obj(config, shape),
        sphere_obj=asset_tpms_sphere(config, shape),
        raw_thickness=-float(target_volume),
        raw_obj=raw_dir / f"{shape}_union.obj",
        raw_quality=raw_dir / "quality" / f"{shape}_union.quality.json",
        raw_log=raw_dir / "logs" / f"{shape}_generate.log",
        raw_quality_log=raw_dir / "logs" / f"{shape}_raw_quality.log",
        remesh_obj=remesh_dir / f"{shape}_remesh.obj",
        remesh_quality=remesh_dir / "quality" / f"{shape}_remesh.quality.json",
        remesh_log=remesh_dir / "logs" / f"{shape}_remesh.log",
        remesh_quality_log=remesh_dir / "logs" / f"{shape}_remesh_quality.log",
        tet_dir=tet_dir,
        component_dir=component_dir,
    )


def read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as fin:
        return require_mapping(json.load(fin), str(path))


def quality_passed(path: Path) -> bool:
    if not path.exists():
        return False
    try:
        return bool(read_json(path).get("passed"))
    except (OSError, json.JSONDecodeError, ValueError):
        return False


def boundary_quality_path(target: PipelineTarget) -> Path:
    if quality_passed(target.veg_quality):
        return target.veg_quality
    if quality_passed(target.repaired_veg_quality):
        return target.repaired_veg_quality
    if target.repaired_veg_quality.exists():
        return target.repaired_veg_quality
    return target.veg_quality


def boundary_surface_path(target: PipelineTarget) -> Path:
    if quality_passed(target.repaired_veg_quality):
        return target.repaired_veg_obj
    if target.repaired_veg_obj.exists() and target.repaired_veg_quality.exists():
        return target.repaired_veg_obj
    return target.veg_obj


def boundary_source(target: PipelineTarget) -> str:
    if quality_passed(target.veg_quality):
        return "original"
    if quality_passed(target.repaired_veg_quality):
        return "repaired"
    if target.repaired_veg_quality.exists():
        return "repaired"
    if target.veg_quality.exists():
        return "original"
    return "missing"


def final_gate_passed(target: PipelineTarget) -> bool:
    if not quality_passed(boundary_quality_path(target)):
        return False
    components = target.component_dir / "components.json"
    if not components.exists():
        return False
    try:
        return len(json.loads(components.read_text(encoding="utf-8"))) == 3
    except (OSError, json.JSONDecodeError):
        return False


def raw_command(config: PipelineConfig, target: PipelineTarget) -> CommandSpec:
    raw = config.parameters["raw"]
    argv = [
        str(tool_path(config, "generate")),
        "openvdb",
        "--fbms",
        str(target.fbms_obj),
        "--sphere",
        str(target.sphere_obj),
        "--fbms-thickness",
        f"{target.raw_thickness:.17g}",
        "--sphere-thickness",
        str(raw["sphere_thickness"]),
        "--resolution",
        str(raw["resolution"]),
        "--padding-ratio",
        str(raw["padding_ratio"]),
        "--output-surface",
        str(target.raw_obj),
    ]
    if raw.get("enable_truncating", True):
        argv.append("--enable-truncating")
    if raw.get("project_fbms_boundary_to_sphere", True):
        argv.append("--project-fbms-boundary-to-sphere")
    if raw.get("filter_small_components", True):
        argv.extend(["--filter-small-components", "--min-component-triangles", str(raw["min_component_triangles"])])
    return CommandSpec(f"raw:{target.case}:{target.name}", argv, target.raw_log)


def quality_command(
    config: PipelineConfig,
    label: str,
    input_mesh: Path,
    output_json: Path,
    log_path: Path,
    check_level: str,
    invalid_policy: str,
    allow_failure: bool,
    self_intersection_backend: str | None = None,
    self_intersection_triangle_limit: int | None = None,
) -> CommandSpec:
    quality = config.parameters["quality"]
    argv = [
        str(tool_path(config, "quality")),
        "surface",
        "--check-level",
        check_level,
        "--expected-components",
        str(quality["expected_components"]),
        "--invalid-triangles-policy",
        invalid_policy,
        "--input",
        str(input_mesh),
        "--json",
        str(output_json),
    ]
    if self_intersection_backend is not None:
        argv.extend(["--self-intersection-backend", self_intersection_backend])
    if self_intersection_triangle_limit is not None:
        argv.extend(["--self-intersection-triangle-limit", str(self_intersection_triangle_limit)])
    return CommandSpec(label, argv, log_path, allow_failure=allow_failure)


def repair_command(target: PipelineTarget) -> CommandSpec:
    return CommandSpec(
        f"boundary-repair:{target.case}:{target.name}",
        [
            sys.executable,
            str(REPO_ROOT / "scripts" / "repair_tpms_mesh.py"),
            str(target.veg_obj),
            "--out",
            str(target.repaired_veg_obj),
            "--report",
            str(target.repair_report),
            "--orient-positive",
        ],
        target.boundary_repair_log,
        allow_failure=True,
    )


def compute_remesh_edge_length_arg(raw_quality: dict[str, Any], target_edge_length: float) -> float:
    edge_length = require_mapping(raw_quality.get("edge_length"), "raw_quality.edge_length")
    mean = edge_length.get("mean")
    if not isinstance(mean, (int, float)) or mean <= 0:
        raise ValueError("raw quality edge_length.mean must be positive")
    return float(target_edge_length) / float(mean)


def remesh_edge_length_text(config: PipelineConfig, target: PipelineTarget, dry_run: bool) -> str:
    remesh = config.parameters["remesh"]
    if target.raw_quality.exists():
        raw_quality = read_json(target.raw_quality)
        return f"{compute_remesh_edge_length_arg(raw_quality, float(remesh['target_edge_length'])):.17g}"
    if dry_run:
        return f"computed:{remesh['target_edge_length']}/raw_mean_edge"
    raise FileNotFoundError(target.raw_quality)


def remesh_command(config: PipelineConfig, target: PipelineTarget, dry_run: bool) -> CommandSpec:
    remesh = config.parameters["remesh"]
    return CommandSpec(
        f"remesh:{target.case}:{target.name}",
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


def relative_path_from(path: Path, base: Path) -> str:
    return os.path.relpath(path, start=base)


def write_tetmesher_config(config: PipelineConfig, target: PipelineTarget, dry_run: bool) -> dict[str, Any]:
    data = {
        "version": 1,
        "backend": "tetwild",
        "input_mesh": relative_path_from(target.remesh_obj, target.tet_dir),
        "output_mesh": f"{target.name}.veg",
        "output_surface": f"{target.name}.veg.obj",
        "print_stats": True,
        "tetwild": dict(config.parameters["tetwild"]),
        "material": dict(config.parameters["material"]),
    }
    log_progress("write", target=target_key(target), kind="tet-config", status="start", path=target.tet_config)
    print(f"+ write {target.tet_config}", flush=True)
    if not dry_run:
        target.tet_config.parent.mkdir(parents=True, exist_ok=True)
        target.tet_config.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    log_progress(
        "write",
        target=target_key(target),
        kind="tet-config",
        status="dry-run" if dry_run else "done",
        path=target.tet_config,
    )
    return data


def tet_command(config: PipelineConfig, target: PipelineTarget) -> CommandSpec:
    return CommandSpec(
        f"tet:{target.case}:{target.name}",
        [str(tool_path(config, "tetmesher")), "--config", str(target.tet_config)],
        target.tet_log,
    )


def run_raw_stage(config: PipelineConfig, target: PipelineTarget, dry_run: bool, overwrite: bool) -> None:
    ensure_tool(tool_path(config, "generate"), dry_run)
    ensure_tool(tool_path(config, "quality"), dry_run)
    start = time.monotonic()
    if overwrite or not quality_passed(target.raw_quality):
        log_progress("stage", stage="raw", target=target_key(target), status="start")
        target.raw_obj.parent.mkdir(parents=True, exist_ok=True) if not dry_run else None
        run_command(raw_command(config, target), dry_run=dry_run)
        quality = config.parameters["quality"]
        run_command(
            quality_command(
                config,
                f"raw-quality:{target.case}:{target.name}",
                target.raw_obj,
                target.raw_quality,
                target.raw_quality_log,
                "raw",
                str(quality["raw_invalid_triangles_policy"]),
                allow_failure=False,
            ),
            dry_run=dry_run,
        )
        log_progress("stage", stage="raw", target=target_key(target), status="done", elapsed=elapsed_text(start))
    else:
        log_progress("stage", stage="raw", target=target_key(target), status="skip", reason="raw-quality-passed")
    if not dry_run and not quality_passed(target.raw_quality):
        raise RuntimeError(f"raw quality failed for {target.kind}/{target.case}/{target.name}")


def run_remesh_stage(config: PipelineConfig, target: PipelineTarget, dry_run: bool, overwrite: bool) -> None:
    ensure_tool(tool_path(config, "remesh"), dry_run)
    ensure_tool(tool_path(config, "quality"), dry_run)
    start = time.monotonic()
    if overwrite or not target.remesh_obj.exists() or not target.remesh_quality.exists():
        log_progress("stage", stage="remesh", target=target_key(target), status="start")
        target.remesh_obj.parent.mkdir(parents=True, exist_ok=True) if not dry_run else None
        run_command(remesh_command(config, target, dry_run=dry_run), dry_run=dry_run)
        quality = config.parameters["quality"]
        run_command(
            quality_command(
                config,
                f"remesh-quality:{target.case}:{target.name}",
                target.remesh_obj,
                target.remesh_quality,
                target.remesh_quality_log,
                "full",
                str(quality["remesh_invalid_triangles_policy"]),
                allow_failure=True,
                self_intersection_backend=str(quality["remesh_self_intersection_backend"]),
                self_intersection_triangle_limit=int(quality["remesh_self_intersection_triangle_limit"]),
            ),
            dry_run=dry_run,
        )
        remesh_status = "unknown"
        if target.remesh_quality.exists():
            remesh_status = "pass" if quality_passed(target.remesh_quality) else "warning"
        log_progress(
            "stage",
            stage="remesh",
            target=target_key(target),
            status="done",
            quality=remesh_status,
            elapsed=elapsed_text(start),
        )
    else:
        log_progress("stage", stage="remesh", target=target_key(target), status="skip", reason="outputs-present")
    if not dry_run and not target.remesh_quality.exists():
        raise RuntimeError(f"remesh quality JSON was not produced for {target.kind}/{target.case}/{target.name}")


def run_tet_stage(config: PipelineConfig, target: PipelineTarget, dry_run: bool, overwrite: bool) -> None:
    ensure_tool(tool_path(config, "tetmesher"), dry_run)
    start = time.monotonic()
    if overwrite or not target.veg.exists() or not target.veg_obj.exists():
        log_progress("stage", stage="tet", target=target_key(target), status="start")
        write_tetmesher_config(config, target, dry_run=dry_run)
        run_command(tet_command(config, target), dry_run=dry_run)
        log_progress("stage", stage="tet", target=target_key(target), status="done", elapsed=elapsed_text(start))
    else:
        log_progress("stage", stage="tet", target=target_key(target), status="skip", reason="outputs-present")


def run_validate_stage(config: PipelineConfig, target: PipelineTarget, dry_run: bool, overwrite: bool) -> None:
    ensure_tool(tool_path(config, "quality"), dry_run)
    ensure_tool(tool_path(config, "volume_info"), dry_run)
    start = time.monotonic()
    log_progress("stage", stage="validate", target=target_key(target), status="start")
    quality = config.parameters["quality"]
    if overwrite or not quality_passed(target.veg_quality):
        run_command(
            quality_command(
                config,
                f"boundary-quality:{target.case}:{target.name}",
                target.veg_obj,
                target.veg_quality,
                target.boundary_quality_log,
                "full",
                str(quality["final_invalid_triangles_policy"]),
                allow_failure=True,
                self_intersection_backend=str(quality["final_self_intersection_backend"]),
                self_intersection_triangle_limit=int(quality["final_self_intersection_triangle_limit"]),
            ),
            dry_run=dry_run,
        )
    else:
        log_progress("validate", target=target_key(target), check="boundary-quality", status="skip")
    if not dry_run and not quality_passed(target.veg_quality):
        if overwrite or not target.repaired_veg_obj.exists() or not target.repaired_veg_quality.exists():
            run_command(repair_command(target), dry_run=dry_run)
            if target.repaired_veg_obj.exists():
                run_command(
                    quality_command(
                        config,
                        f"boundary-repaired-quality:{target.case}:{target.name}",
                        target.repaired_veg_obj,
                        target.repaired_veg_quality,
                        target.boundary_repaired_quality_log,
                        "full",
                        str(quality["final_invalid_triangles_policy"]),
                        allow_failure=True,
                        self_intersection_backend=str(quality["final_self_intersection_backend"]),
                        self_intersection_triangle_limit=int(quality["final_self_intersection_triangle_limit"]),
                    ),
                    dry_run=dry_run,
                )
        else:
            log_progress("validate", target=target_key(target), check="boundary-repair", status="skip")
    elif dry_run:
        log_progress("validate", target=target_key(target), check="boundary-repair", status="dry-run-if-needed")
    if overwrite or not target.veg_info.exists():
        run_command_capture([str(tool_path(config, "volume_info")), str(target.veg)], target.veg_info, dry_run=dry_run)
    else:
        log_progress("validate", target=target_key(target), check="volume-info", status="skip")
    component_input = boundary_surface_path(target)
    if overwrite or component_input == target.repaired_veg_obj or not (target.component_dir / "components.json").exists():
        run_command(
            CommandSpec(
                f"components:{target.case}:{target.name}",
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "dump_obj_components.py"),
                    "--input",
                    str(component_input),
                    "--output-dir",
                    str(target.component_dir),
                    "--overwrite",
                ],
                target.components_log,
            ),
            dry_run=dry_run,
        )
    else:
        log_progress("validate", target=target_key(target), check="components", status="skip")
    boundary_status = "unknown"
    final_quality = boundary_quality_path(target)
    if final_quality.exists():
        boundary_status = "pass" if quality_passed(final_quality) else "warning"
    if not dry_run and not final_gate_passed(target):
        log_progress(
            "warning",
            label=f"boundary-gate:{target.case}:{target.name}",
            reason="final-gate-not-passed",
        )
    log_progress(
        "stage",
        stage="validate",
        target=target_key(target),
        status="done",
        quality=boundary_status,
        elapsed=elapsed_text(start),
    )


def run_target(
    config: PipelineConfig,
    target: PipelineTarget,
    from_stage: str,
    dry_run: bool,
    overwrite: bool,
) -> None:
    if not overwrite and from_stage == "prepare" and final_gate_passed(target):
        log_progress("target", target=target_key(target), status="skip", reason="final-boundary-passed")
        return
    if stage_enabled(from_stage, "raw"):
        run_raw_stage(config, target, dry_run=dry_run, overwrite=overwrite)
    if stage_enabled(from_stage, "remesh"):
        run_remesh_stage(config, target, dry_run=dry_run, overwrite=overwrite)
    if stage_enabled(from_stage, "tet"):
        run_tet_stage(config, target, dry_run=dry_run, overwrite=overwrite)
    if stage_enabled(from_stage, "validate"):
        run_validate_stage(config, target, dry_run=dry_run, overwrite=overwrite)


SELECTED_THICKNESS_RE = re.compile(r"selected thickness = ([0-9eE+\-.]+)")
INFO_INT_RE = re.compile(r"^#(vtx|elements):\s+(\d+)", re.MULTILINE)
INFO_VOLUME_RE = re.compile(r"^total volume:\s+([0-9eE+\-.]+)", re.MULTILINE)


def parse_selected_thickness(target: PipelineTarget) -> float | None:
    if target.kind == "fbms":
        return target.raw_thickness
    if not target.raw_log.exists():
        return None
    match = SELECTED_THICKNESS_RE.search(target.raw_log.read_text(encoding="utf-8", errors="replace"))
    return float(match.group(1)) if match else None


def parse_volume_info(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"vertices": None, "elements": None, "volume": None}
    text = path.read_text(encoding="utf-8", errors="replace")
    ints = {match.group(1): int(match.group(2)) for match in INFO_INT_RE.finditer(text)}
    volume_match = INFO_VOLUME_RE.search(text)
    return {
        "vertices": ints.get("vtx"),
        "elements": ints.get("elements"),
        "volume": float(volume_match.group(1)) if volume_match else None,
    }


def component_count(path: Path) -> int | None:
    if not path.exists():
        return None
    try:
        return len(json.loads(path.read_text(encoding="utf-8")))
    except (OSError, json.JSONDecodeError):
        return None


def repair_report_abs_volume(report: dict[str, Any], key: str) -> float | None:
    section = report.get(key)
    if not isinstance(section, dict):
        return None
    volumes = section.get("component_signed_volumes")
    if not isinstance(volumes, list):
        return None
    total = 0.0
    for value in volumes:
        if not isinstance(value, (int, float)):
            return None
        total += abs(float(value))
    return total


def repair_volume_drift(
    original_boundary_quality: dict[str, Any] | None,
    repaired_boundary_quality: dict[str, Any] | None,
    repair_report_path: Path,
) -> float | None:
    original_volume = original_boundary_quality.get("enclosed_volume") if original_boundary_quality else None
    repaired_volume = repaired_boundary_quality.get("enclosed_volume") if repaired_boundary_quality else None
    if not (
        isinstance(original_volume, (int, float))
        and isinstance(repaired_volume, (int, float))
        and abs(float(original_volume)) > 0
    ):
        if not repair_report_path.exists():
            return None
        try:
            report = read_json(repair_report_path)
        except (OSError, json.JSONDecodeError, ValueError):
            return None
        original_volume = repair_report_abs_volume(report, "before")
        repaired_volume = repair_report_abs_volume(report, "final")
    if (
        isinstance(original_volume, (int, float))
        and isinstance(repaired_volume, (int, float))
        and abs(float(original_volume)) > 0
    ):
        return abs(float(repaired_volume) - float(original_volume)) / abs(float(original_volume))
    return None


def warning_summary(quality: dict[str, Any] | None) -> str:
    if quality is None:
        return "missing"
    self_intersections = quality.get("self_intersections")
    if quality.get("passed"):
        return "pass"
    if isinstance(self_intersections, int) and self_intersections > 0:
        return f"warn: {self_intersections} self-intersections"
    return "warn"


def collect_target_summary(target: PipelineTarget) -> dict[str, Any]:
    raw_quality = read_json(target.raw_quality) if target.raw_quality.exists() else None
    remesh_quality = read_json(target.remesh_quality) if target.remesh_quality.exists() else None
    original_boundary_quality = read_json(target.veg_quality) if target.veg_quality.exists() else None
    repaired_boundary_quality = read_json(target.repaired_veg_quality) if target.repaired_veg_quality.exists() else None
    boundary_quality = read_json(boundary_quality_path(target)) if boundary_quality_path(target).exists() else None
    info = parse_volume_info(target.veg_info)
    components = component_count(target.component_dir / "components.json")
    volume_drift = repair_volume_drift(original_boundary_quality, repaired_boundary_quality, target.repair_report)
    return {
        "kind": target.kind,
        "case": target.case,
        "name": target.name,
        "selected_thickness": parse_selected_thickness(target),
        "raw_volume": raw_quality.get("enclosed_volume") if raw_quality else None,
        "raw_passed": raw_quality.get("passed") if raw_quality else None,
        "raw_invalid_triangles": raw_quality.get("invalid_triangles") if raw_quality else None,
        "remesh_mean_edge": (
            remesh_quality.get("edge_length", {}).get("mean") if remesh_quality else None
        ),
        "remesh_passed": remesh_quality.get("passed") if remesh_quality else None,
        "remesh_warning": warning_summary(remesh_quality),
        "tet_vertices": info["vertices"],
        "tet_elements": info["elements"],
        "tet_volume": info["volume"],
        "boundary_passed": boundary_quality.get("passed") if boundary_quality else None,
        "boundary_source": boundary_source(target),
        "boundary_repair_attempted": target.repair_report.exists() or target.repaired_veg_quality.exists(),
        "boundary_repair_passed": (
            repaired_boundary_quality.get("passed") if repaired_boundary_quality else None
        ),
        "boundary_repair_volume_drift": volume_drift,
        "component_count": components,
    }


def collect_summary(config: PipelineConfig, cases: tuple[SourceCase, ...]) -> dict[str, Any]:
    results: list[dict[str, Any]] = []
    for case in cases:
        fbms = fbms_target(config, case.name)
        results.append(collect_target_summary(fbms))
        budget_quality = read_json(fbms.raw_quality) if fbms.raw_quality.exists() else {}
        budget = budget_quality.get("enclosed_volume")
        if isinstance(budget, (int, float)):
            for shape in config.tpms_shapes:
                results.append(collect_target_summary(baseline_target(config, case.name, shape, float(budget))))

    return {
        "generated_at_unix": int(time.time()),
        "working_directory": str(config.working_directory),
        "results": results,
    }


def format_float(value: Any, digits: int = 6) -> str:
    if isinstance(value, (int, float)):
        return f"{float(value):.{digits}f}"
    return "-"


def markdown_results_table(summary: dict[str, Any]) -> str:
    lines = [
        "| asset | raw volume | selected thickness | remesh mean edge | remesh status | tet vertices | tet elements | tet volume | boundary | components |",
        "| --- | ---: | ---: | ---: | --- | ---: | ---: | ---: | --- | ---: |",
    ]
    for result in summary.get("results", []):
        asset = result["case"] if result["kind"] == "fbms" else f"{result['case']} / {result['name']}"
        if result.get("boundary_passed") and result.get("boundary_source") == "repaired":
            boundary = "pass(repaired)"
        else:
            boundary = "pass" if result.get("boundary_passed") else "fail/missing"
        lines.append(
            "| "
            + " | ".join(
                [
                    asset,
                    format_float(result.get("raw_volume"), 9),
                    format_float(result.get("selected_thickness"), 7),
                    format_float(result.get("remesh_mean_edge"), 6),
                    str(result.get("remesh_warning") or "-"),
                    str(result.get("tet_vertices") or "-"),
                    str(result.get("tet_elements") or "-"),
                    format_float(result.get("tet_volume"), 6),
                    boundary,
                    str(result.get("component_count") or "-"),
                ]
            )
            + " |"
        )
    return "\n".join(lines)


def replace_readme_marker_block(text: str, replacement: str) -> str:
    if README_RESULTS_START not in text or README_RESULTS_END not in text:
        block = f"\n\n{README_RESULTS_START}\n{replacement}\n{README_RESULTS_END}\n"
        return text.rstrip() + block
    start = text.index(README_RESULTS_START) + len(README_RESULTS_START)
    end = text.index(README_RESULTS_END)
    return text[:start] + "\n" + replacement + "\n" + text[end:]


def update_readme(config: PipelineConfig, summary: dict[str, Any], dry_run: bool) -> None:
    readme = config.working_directory / "README.md"
    replacement = markdown_results_table(summary)
    text = readme.read_text(encoding="utf-8") if readme.exists() else "# FBMS All Asset Pipeline\n"
    updated = replace_readme_marker_block(text, replacement)
    log_progress("write", kind="readme", status="start", path=readme)
    print(f"+ update {readme}", flush=True)
    if not dry_run:
        readme.write_text(updated, encoding="utf-8")
    log_progress("write", kind="readme", status="dry-run" if dry_run else "done", path=readme)


def write_summary(config: PipelineConfig, summary: dict[str, Any], dry_run: bool) -> Path:
    path = summary_path(config)
    log_progress("write", kind="summary", status="start", path=path)
    print(f"+ write {path}", flush=True)
    if not dry_run:
        path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    log_progress("write", kind="summary", status="dry-run" if dry_run else "done", path=path)
    return path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the FBMS-all r512 asset pipeline from JSON.")
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG, help="Pipeline JSON config.")
    parser.add_argument("--case", action="append", default=[], help="Restrict to one FBMS case. Can be repeated.")
    parser.add_argument("--shape", action="append", default=[], help="Restrict to one TPMS shape. Can be repeated.")
    parser.add_argument("--from-stage", choices=STAGES, default="prepare", help="Start from this stage.")
    parser.add_argument("--dry-run", action="store_true", help="Print planned commands without writing assets.")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing stage outputs.")
    parser.add_argument("--update-readme", action="store_true", help="Refresh the README results marker block.")
    return parser.parse_args()


def filter_cases(cases: tuple[SourceCase, ...], selected: list[str]) -> tuple[SourceCase, ...]:
    if not selected:
        return cases
    requested = set(selected)
    filtered = tuple(case for case in cases if case.name in requested)
    missing = requested.difference(case.name for case in filtered)
    if missing:
        raise ValueError(f"unknown cases: {', '.join(sorted(missing))}")
    return filtered


def filter_shapes(config: PipelineConfig, selected: list[str]) -> tuple[str, ...]:
    if not selected:
        return config.tpms_shapes
    requested = set(selected)
    missing = requested.difference(config.tpms_shapes)
    if missing:
        raise ValueError(f"unknown TPMS shapes: {', '.join(sorted(missing))}")
    return tuple(shape for shape in config.tpms_shapes if shape in requested)


def run_pipeline(args: argparse.Namespace) -> dict[str, Any]:
    pipeline_start = time.monotonic()
    config = load_config(args.config)
    cases = filter_cases(discover_fbms_cases(config.source_asset_directory), args.case)
    shapes = filter_shapes(config, args.shape)
    config = PipelineConfig(
        path=config.path,
        source_asset_directory=config.source_asset_directory,
        working_directory=config.working_directory,
        build_dir=config.build_dir,
        parameters=config.parameters,
        tpms_shapes=shapes,
    )
    log_progress(
        "pipeline",
        status="start",
        config=config.path,
        working_directory=config.working_directory,
        cases=",".join(case.name for case in cases),
        shapes=",".join(shapes),
        from_stage=args.from_stage,
        overwrite=str(args.overwrite).lower(),
        dry_run=str(args.dry_run).lower(),
    )

    if stage_enabled(args.from_stage, "prepare"):
        prepare_assets(config, cases, dry_run=args.dry_run, overwrite=args.overwrite)

    for case in cases:
        log_progress("target", target=f"fbms/{case.name}/{case.name}", status="start")
        print(f"== FBMS {case.name} ==", flush=True)
        fbms = fbms_target(config, case.name)
        if args.from_stage != "summary":
            run_target(config, fbms, args.from_stage, dry_run=args.dry_run, overwrite=args.overwrite)

        if fbms.raw_quality.exists():
            budget_quality = read_json(fbms.raw_quality)
            budget = budget_quality.get("enclosed_volume")
            if not isinstance(budget, (int, float)) or budget <= 0:
                raise ValueError(f"{fbms.raw_quality} must contain positive enclosed_volume")
        elif args.dry_run:
            budget = 1.0
        else:
            budget_quality = read_json(fbms.raw_quality)
            budget = budget_quality.get("enclosed_volume")
            if not isinstance(budget, (int, float)) or budget <= 0:
                raise ValueError(f"{fbms.raw_quality} must contain positive enclosed_volume")

        for shape in shapes:
            log_progress("target", target=f"baseline/{case.name}/{shape}", status="start")
            print(f"== Baseline {case.name} / {shape} ==", flush=True)
            target = baseline_target(config, case.name, shape, float(budget))
            if args.from_stage != "summary":
                run_target(config, target, args.from_stage, dry_run=args.dry_run, overwrite=args.overwrite)

    summary = collect_summary(config, cases)
    if stage_enabled(args.from_stage, "summary"):
        write_summary(config, summary, dry_run=args.dry_run)
    if args.update_readme:
        update_readme(config, summary, dry_run=args.dry_run)
    log_progress("pipeline", status="done", elapsed=elapsed_text(pipeline_start))
    return summary


def main() -> int:
    args = parse_args()
    try:
        run_pipeline(args)
        return 0
    except (
        FileNotFoundError,
        FileExistsError,
        ValueError,
        RuntimeError,
        json.JSONDecodeError,
        subprocess.CalledProcessError,
    ) as exc:
        print(f"error: {exc}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
