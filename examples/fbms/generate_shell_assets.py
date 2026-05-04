#!/usr/bin/env python3
"""Generate FBMS SDF-union raw surfaces and remeshed shell assets from JSON jobs."""

from __future__ import annotations

import argparse
import json
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[2]
FBMS_ROOT = REPO_ROOT / "examples" / "fbms"
DEFAULT_CONFIG = FBMS_ROOT / "shell_assets.json"


@dataclass(frozen=True)
class AssetCase:
    name: str
    directory: Path
    prefix: str

    @property
    def fbms_obj(self) -> Path:
        return self.directory / f"{self.prefix}_fbms.obj"

    @property
    def sphere_obj(self) -> Path:
        return self.directory / f"{self.prefix}_fbms_bounding_sphere.obj"


@dataclass(frozen=True)
class RawParams:
    resolution: int
    fbms_thickness: float
    sphere_thickness: float
    padding_ratio: float
    enable_truncating: bool


@dataclass(frozen=True)
class RemeshParams:
    edge_length: float
    sharp_edge_angle: float


@dataclass(frozen=True)
class AssetJob:
    name: str
    cases: tuple[str, ...]
    raw: RawParams
    remesh: RemeshParams
    output_dir: str
    raw_filename: str
    remesh_filename: str
    stats_filename: str


CASES = {
    "g0_b3": AssetCase("g0_b3", FBMS_ROOT / "g0_b3", "g0_b3"),
    "g0_b8": AssetCase("g0_b8", FBMS_ROOT / "g0_b8", "g0_b8"),
}


def compact_float(value: float) -> str:
    return f"{value:.12g}"


def require_mapping(value: Any, context: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ValueError(f"{context} must be an object")
    return value


def require_string(value: Any, context: str) -> str:
    if not isinstance(value, str) or not value:
        raise ValueError(f"{context} must be a non-empty string")
    return value


def require_number(value: Any, context: str) -> float:
    if not isinstance(value, (int, float)):
        raise ValueError(f"{context} must be a number")
    return float(value)


def require_int(value: Any, context: str) -> int:
    if not isinstance(value, int):
        raise ValueError(f"{context} must be an integer")
    return value


def require_bool(value: Any, context: str) -> bool:
    if not isinstance(value, bool):
        raise ValueError(f"{context} must be a boolean")
    return value


def repo_path(path_text: str) -> Path:
    path = Path(path_text)
    return path if path.is_absolute() else REPO_ROOT / path


def fbms_path(path_text: str) -> Path:
    path = Path(path_text)
    return path if path.is_absolute() else FBMS_ROOT / path


def merge_mapping(defaults: dict[str, Any], overrides: dict[str, Any]) -> dict[str, Any]:
    merged = dict(defaults)
    merged.update(overrides)
    return merged


def parse_cases(value: Any, context: str) -> tuple[str, ...]:
    if value == "all":
        return tuple(CASES)
    if not isinstance(value, list) or not value:
        raise ValueError(f"{context} must be a non-empty list or \"all\"")
    names = tuple(require_string(case, f"{context}[]") for case in value)
    unknown = [name for name in names if name not in CASES]
    if unknown:
        raise ValueError(f"{context} contains unknown cases: {', '.join(unknown)}")
    return names


def parse_job(job_config: dict[str, Any], defaults: dict[str, Any]) -> AssetJob:
    name = require_string(job_config.get("name"), "job.name")
    raw_defaults = require_mapping(defaults.get("raw", {}), "defaults.raw")
    remesh_defaults = require_mapping(defaults.get("remesh", {}), "defaults.remesh")
    raw_config = merge_mapping(raw_defaults, require_mapping(job_config.get("raw"), f"job {name}.raw"))
    remesh_config = merge_mapping(
        remesh_defaults,
        require_mapping(job_config.get("remesh", {}), f"job {name}.remesh"),
    )

    resolution = require_int(raw_config.get("resolution"), f"job {name}.raw.resolution")
    fbms_thickness = require_number(raw_config.get("fbms_thickness"), f"job {name}.raw.fbms_thickness")
    sphere_thickness = require_number(
        raw_config.get("sphere_thickness", fbms_thickness),
        f"job {name}.raw.sphere_thickness",
    )
    padding_ratio = require_number(raw_config.get("padding_ratio"), f"job {name}.raw.padding_ratio")
    enable_truncating = require_bool(
        raw_config.get("enable_truncating", False),
        f"job {name}.raw.enable_truncating",
    )
    edge_length = require_number(remesh_config.get("edge_length"), f"job {name}.remesh.edge_length")
    sharp_edge_angle = require_number(
        remesh_config.get("sharp_edge_angle"),
        f"job {name}.remesh.sharp_edge_angle",
    )

    if resolution < 2:
        raise ValueError(f"job {name}.raw.resolution must be at least 2")
    if fbms_thickness == 0.0:
        raise ValueError(f"job {name}.raw.fbms_thickness must be non-zero")
    if sphere_thickness <= 0.0:
        raise ValueError(f"job {name}.raw.sphere_thickness must be positive")
    if padding_ratio < 0.0:
        raise ValueError(f"job {name}.raw.padding_ratio must be non-negative")
    if edge_length <= 0.0:
        raise ValueError(f"job {name}.remesh.edge_length must be positive")

    return AssetJob(
        name=name,
        cases=parse_cases(job_config.get("cases", defaults.get("cases", "all")), f"job {name}.cases"),
        raw=RawParams(
            resolution=resolution,
            fbms_thickness=fbms_thickness,
            sphere_thickness=sphere_thickness,
            padding_ratio=padding_ratio,
            enable_truncating=enable_truncating,
        ),
        remesh=RemeshParams(edge_length=edge_length, sharp_edge_angle=sharp_edge_angle),
        output_dir=require_string(
            job_config.get("output_dir", defaults.get("output_dir", "generated/{job}/{case}")),
            f"job {name}.output_dir",
        ),
        raw_filename=require_string(
            job_config.get("raw_filename", defaults.get("raw_filename", "union_shell_raw.obj")),
            f"job {name}.raw_filename",
        ),
        remesh_filename=require_string(
            job_config.get("remesh_filename", defaults.get("remesh_filename", "union_shell_remesh.obj")),
            f"job {name}.remesh_filename",
        ),
        stats_filename=require_string(
            job_config.get("stats_filename", defaults.get("stats_filename", "stats.json")),
            f"job {name}.stats_filename",
        ),
    )


def load_config(config_path: Path) -> tuple[Path, dict[str, AssetJob]]:
    with config_path.open("r", encoding="utf-8") as fin:
        config = require_mapping(json.load(fin), "config")

    defaults = require_mapping(config.get("defaults", {}), "defaults")
    build_dir = repo_path(require_string(config.get("build_dir", "build/base_no_mkl"), "build_dir"))
    jobs_config = config.get("jobs")
    if not isinstance(jobs_config, list) or not jobs_config:
        raise ValueError("jobs must be a non-empty array")

    jobs: dict[str, AssetJob] = {}
    for index, job_value in enumerate(jobs_config):
        job = parse_job(require_mapping(job_value, f"jobs[{index}]"), defaults)
        if job.name in jobs:
            raise ValueError(f"duplicate job name: {job.name}")
        jobs[job.name] = job

    return build_dir, jobs


def render_output_dir(job: AssetJob, case: AssetCase) -> Path:
    text = job.output_dir.format(job=job.name, case=case.name)
    return fbms_path(text)


def output_paths(job: AssetJob, case: AssetCase) -> tuple[Path, Path, Path]:
    directory = render_output_dir(job, case)
    return directory / job.raw_filename, directory / job.remesh_filename, directory / job.stats_filename


def run_command(command: list[str], dry_run: bool) -> None:
    printable = " ".join(command)
    print(f"+ {printable}")
    if dry_run:
        return
    subprocess.run(command, check=True, cwd=REPO_ROOT)


def count_obj(path: Path) -> tuple[int, int]:
    vertices = 0
    faces = 0
    with path.open("r", encoding="utf-8", errors="replace") as fin:
        for line in fin:
            if line.startswith("v "):
                vertices += 1
            elif line.startswith("f "):
                faces += 1
    return vertices, faces


def check_inputs(case: AssetCase) -> None:
    for path in (case.fbms_obj, case.sphere_obj):
        if not path.exists():
            raise FileNotFoundError(path)


def check_output_policy(paths: tuple[Path, ...], overwrite: bool, skip_existing: bool) -> bool:
    existing = [path for path in paths if path.exists()]
    if not existing:
        return False
    if overwrite:
        return False
    if skip_existing:
        print("Skipping existing outputs:")
        for path in existing:
            print(f"  {path.relative_to(REPO_ROOT)}")
        return True
    joined = ", ".join(str(path.relative_to(REPO_ROOT)) for path in existing)
    raise FileExistsError(f"output already exists: {joined}; pass --overwrite or --skip-existing")


def write_stats(job: AssetJob, case: AssetCase, raw_path: Path, remesh_path: Path, stats_path: Path) -> None:
    assets = {}
    for label, path in (("raw", raw_path), ("remesh", remesh_path)):
        if path.exists():
            vertices, faces = count_obj(path)
            assets[label] = {
                "path": str(path.relative_to(REPO_ROOT)),
                "vertices": vertices,
                "faces": faces,
            }
            print(f"{job.name} {case.name} {label}: {path.relative_to(REPO_ROOT)} vertices={vertices} faces={faces}")

    stats = {
        "job": job.name,
        "case": case.name,
        "raw": {
            "resolution": job.raw.resolution,
            "fbms_thickness": job.raw.fbms_thickness,
            "sphere_thickness": job.raw.sphere_thickness,
            "padding_ratio": job.raw.padding_ratio,
            "enable_truncating": job.raw.enable_truncating,
        },
        "remesh": {
            "edge_length": job.remesh.edge_length,
            "sharp_edge_angle": job.remesh.sharp_edge_angle,
        },
        "assets": assets,
    }
    with stats_path.open("w", encoding="utf-8") as fout:
        json.dump(stats, fout, indent=2)
        fout.write("\n")


def generate_asset(args: argparse.Namespace, build_dir: Path, job: AssetJob, case: AssetCase) -> None:
    check_inputs(case)

    raw_path, remesh_path, stats_path = output_paths(job, case)
    policy_paths = (remesh_path,) if args.remesh_only else (raw_path,) if args.raw_only else (raw_path, remesh_path)
    if check_output_policy(policy_paths, args.overwrite, args.skip_existing):
        return
    if args.remesh_only and not args.dry_run and not raw_path.exists():
        raise FileNotFoundError(raw_path)

    if not args.dry_run:
        raw_path.parent.mkdir(parents=True, exist_ok=True)

    generate_bin = build_dir / "bin" / "generateFBMSUnionSurface"
    remesh_bin = build_dir / "bin" / "remeshSurface"
    if not args.dry_run:
        tools = (remesh_bin,) if args.remesh_only else (generate_bin,) if args.raw_only else (generate_bin, remesh_bin)
        for tool in tools:
            if not tool.exists():
                raise FileNotFoundError(tool)

    if not args.remesh_only:
        command = [
            str(generate_bin),
            "--fbms",
            str(case.fbms_obj),
            "--sphere",
            str(case.sphere_obj),
            "--fbms-thickness",
            compact_float(job.raw.fbms_thickness),
            "--sphere-thickness",
            compact_float(job.raw.sphere_thickness),
            "--resolution",
            str(job.raw.resolution),
            "--padding-ratio",
            compact_float(job.raw.padding_ratio),
            "--output-surface",
            str(raw_path),
        ]
        if job.raw.enable_truncating:
            command.append("--enable-truncating")
        run_command(command, args.dry_run)

    if not args.raw_only:
        run_command(
            [
                str(remesh_bin),
                "cgal_iso",
                "--input-mesh",
                str(raw_path),
                "--output-mesh",
                str(remesh_path),
                "--edge-length",
                compact_float(job.remesh.edge_length),
                "--sharp-edge-angle",
                compact_float(job.remesh.sharp_edge_angle),
            ],
            args.dry_run,
        )

    if not args.dry_run:
        write_stats(job, case, raw_path, remesh_path, stats_path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate FBMS raw shell OBJ files and remeshed OBJ files from JSON jobs."
    )
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG, help="JSON job config file.")
    parser.add_argument("--job", action="append", default=[], help="Job name to run. Can be repeated.")
    parser.add_argument("--all-jobs", action="store_true", help="Run every job from the JSON config.")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without running them.")
    parser.add_argument("--overwrite", action="store_true", help="Allow existing outputs to be overwritten.")
    parser.add_argument("--skip-existing", action="store_true", help="Skip jobs whose output files already exist.")
    parser.add_argument("--raw-only", action="store_true", help="Generate raw OBJ only.")
    parser.add_argument("--remesh-only", action="store_true", help="Run remesh only, assuming raw OBJ exists.")
    args = parser.parse_args()

    if args.overwrite and args.skip_existing:
        parser.error("--overwrite and --skip-existing are mutually exclusive")
    if args.raw_only and args.remesh_only:
        parser.error("--raw-only and --remesh-only are mutually exclusive")
    if args.all_jobs and args.job:
        parser.error("--all-jobs cannot be combined with --job")
    if not args.all_jobs and not args.job:
        parser.error("pass at least one --job, or pass --all-jobs explicitly")

    args.config = args.config if args.config.is_absolute() else REPO_ROOT / args.config
    return args


def main() -> int:
    args = parse_args()
    try:
        build_dir, jobs = load_config(args.config)
        selected_jobs = list(jobs) if args.all_jobs else args.job
        unknown_jobs = [name for name in selected_jobs if name not in jobs]
        if unknown_jobs:
            raise ValueError(f"unknown jobs: {', '.join(unknown_jobs)}")

        for job_name in selected_jobs:
            job = jobs[job_name]
            for case_name in job.cases:
                case = CASES[case_name]
                print(
                    f"== {job.name}/{case.name}: resolution={job.raw.resolution}, "
                    f"fbms_thickness={job.raw.fbms_thickness}, sphere_thickness={job.raw.sphere_thickness}, "
                    f"enable_truncating={job.raw.enable_truncating}, edge_length={job.remesh.edge_length} =="
                )
                generate_asset(args, build_dir, job, case)
        return 0
    except (FileNotFoundError, FileExistsError, ValueError, json.JSONDecodeError, subprocess.CalledProcessError) as exc:
        print(f"error: {exc}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
