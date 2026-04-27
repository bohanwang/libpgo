#!/usr/bin/env python3
"""Run runIPCSim and animation-conversion batches from a JSON config."""

from __future__ import annotations

import argparse
import json
import shlex
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = REPO_ROOT / "examples" / "ipc" / "ipc_batch.json"
DEFAULT_STAGES = ("sim", "abc")
STAGE_ORDER = ("sim", "abc", "vtu")


@dataclass(frozen=True)
class CaseConfig:
    name: str
    sim_config: Path
    anim_config: Path
    vtu_config: Path | None
    log: bool


@dataclass(frozen=True)
class BatchJob:
    name: str
    cases: tuple[str, ...]
    stages: tuple[str, ...]


@dataclass(frozen=True)
class CommandSpec:
    label: str
    argv: list[str]


@dataclass(frozen=True)
class CaseResult:
    name: str
    status: str
    elapsed_seconds: float


def require_mapping(value: Any, context: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ValueError(f"{context} must be an object")
    return value


def require_string(value: Any, context: str) -> str:
    if not isinstance(value, str) or not value:
        raise ValueError(f"{context} must be a non-empty string")
    return value


def require_bool(value: Any, context: str) -> bool:
    if not isinstance(value, bool):
        raise ValueError(f"{context} must be a boolean")
    return value


def repo_path(path_text: str) -> Path:
    path = Path(path_text)
    return path if path.is_absolute() else REPO_ROOT / path


def resolve_anim_path(anim_text: str, sim_config: Path) -> Path:
    path = Path(anim_text)
    if path.is_absolute():
        return path

    case_relative = sim_config.parent / path
    if case_relative.exists():
        return case_relative
    return repo_path(anim_text)


def resolve_case_relative_path(path_text: str, sim_config: Path) -> Path:
    path = Path(path_text)
    if path.is_absolute():
        return path

    case_relative = sim_config.parent / path
    if case_relative.exists():
        return case_relative
    return repo_path(path_text)


def parse_stages(value: Any, context: str) -> tuple[str, ...]:
    if value is None:
        return DEFAULT_STAGES
    if not isinstance(value, list) or not value:
        raise ValueError(f"{context} must be a non-empty list")

    requested = tuple(require_string(stage, f"{context}[]") for stage in value)
    unknown = [stage for stage in requested if stage not in STAGE_ORDER]
    if unknown:
        raise ValueError(f"{context} contains unknown stages: {', '.join(unknown)}")
    duplicates = sorted({stage for stage in requested if requested.count(stage) > 1})
    if duplicates:
        raise ValueError(f"{context} contains duplicate stages: {', '.join(duplicates)}")
    return tuple(stage for stage in STAGE_ORDER if stage in requested)


def parse_cases(value: Any, known_cases: dict[str, CaseConfig], context: str) -> tuple[str, ...]:
    if value == "all":
        return tuple(known_cases)
    if not isinstance(value, list) or not value:
        raise ValueError(f"{context} must be a non-empty list or \"all\"")

    case_names = tuple(require_string(case, f"{context}[]") for case in value)
    unknown = [case for case in case_names if case not in known_cases]
    if unknown:
        raise ValueError(f"{context} contains unknown cases: {', '.join(unknown)}")
    return case_names


def parse_case(name: str, case_config: dict[str, Any], defaults: dict[str, Any]) -> CaseConfig:
    sim_config = repo_path(require_string(case_config.get("sim_config"), f"case {name}.sim_config"))
    anim_config = resolve_anim_path(
        require_string(case_config.get("anim_config", defaults.get("anim_config", "anim.json")), f"case {name}.anim_config"),
        sim_config,
    )
    vtu_config_value = case_config.get("vtu_config", defaults.get("vtu_config"))
    vtu_config = None
    if vtu_config_value is not None:
        vtu_config = resolve_case_relative_path(require_string(vtu_config_value, f"case {name}.vtu_config"), sim_config)
    log = require_bool(case_config.get("log", defaults.get("log", True)), f"case {name}.log")
    return CaseConfig(name=name, sim_config=sim_config, anim_config=anim_config, vtu_config=vtu_config, log=log)


def load_config(config_path: Path) -> tuple[Path, dict[str, CaseConfig], dict[str, BatchJob]]:
    with config_path.open("r", encoding="utf-8") as fin:
        config = require_mapping(json.load(fin), "config")

    defaults = require_mapping(config.get("defaults", {}), "defaults")
    build_dir = repo_path(require_string(config.get("build_dir", "build/base_no_mkl"), "build_dir"))

    cases_config = require_mapping(config.get("cases"), "cases")
    cases: dict[str, CaseConfig] = {}
    for name, case_value in cases_config.items():
        case_name = require_string(name, "case name")
        cases[case_name] = parse_case(case_name, require_mapping(case_value, f"case {case_name}"), defaults)

    jobs_config = config.get("jobs")
    if not isinstance(jobs_config, list) or not jobs_config:
        raise ValueError("jobs must be a non-empty array")

    jobs: dict[str, BatchJob] = {}
    for index, job_value in enumerate(jobs_config):
        job_config = require_mapping(job_value, f"jobs[{index}]")
        name = require_string(job_config.get("name"), f"jobs[{index}].name")
        if name in jobs:
            raise ValueError(f"duplicate job name: {name}")
        jobs[name] = BatchJob(
            name=name,
            cases=parse_cases(job_config.get("cases", "all"), cases, f"job {name}.cases"),
            stages=parse_stages(job_config.get("stages"), f"job {name}.stages"),
        )

    return build_dir, cases, jobs


def build_commands(build_dir: Path, case: CaseConfig, stages: tuple[str, ...], overwrite: bool) -> list[CommandSpec]:
    commands: list[CommandSpec] = []
    if "sim" in stages:
        argv = [str(build_dir / "bin" / "runIPCSim"), str(case.sim_config)]
        if case.log:
            argv.append("--log")
        commands.append(CommandSpec("sim", argv))
    if "abc" in stages:
        commands.append(CommandSpec("abc", [str(build_dir / "bin" / "convertAnimation"), str(case.anim_config)]))
    if "vtu" in stages:
        if case.vtu_config is None:
            raise ValueError(f"case {case.name} needs vtu_config for vtu stage")
        argv = [str(REPO_ROOT / "scripts" / "export_fbms_stress_vtu.py"), "--config", str(case.vtu_config)]
        if overwrite:
            argv.append("--overwrite")
        commands.append(CommandSpec("vtu", argv))
    return commands


def sim_output_dir(case: CaseConfig) -> Path:
    with case.sim_config.open("r", encoding="utf-8") as fin:
        config = require_mapping(json.load(fin), str(case.sim_config))
    output = Path(require_string(config.get("output"), f"{case.sim_config}.output"))
    return output if output.is_absolute() else case.sim_config.parent / output


def check_case_inputs(case: CaseConfig, stages: tuple[str, ...]) -> None:
    if "sim" in stages and not case.sim_config.exists():
        raise FileNotFoundError(case.sim_config)
    if "abc" in stages and not case.anim_config.exists():
        raise FileNotFoundError(case.anim_config)
    if "vtu" in stages:
        if case.vtu_config is None:
            raise ValueError(f"case {case.name} needs vtu_config for vtu stage")
        if not case.vtu_config.exists():
            raise FileNotFoundError(case.vtu_config)


def check_tools(commands: list[CommandSpec]) -> None:
    for command in commands:
        tool = Path(command.argv[0])
        if not tool.exists():
            raise FileNotFoundError(tool)


def check_output_policy(case: CaseConfig, overwrite: bool, skip_existing: bool, stages: tuple[str, ...]) -> bool:
    if "sim" not in stages:
        return False

    output_dir = sim_output_dir(case)
    if not output_dir.exists():
        return False
    if overwrite:
        return False
    if skip_existing:
        print(f"Skipping {case.name}: output exists at {output_dir.relative_to(REPO_ROOT)}")
        return True
    raise FileExistsError(
        f"output already exists: {output_dir.relative_to(REPO_ROOT)}; pass --overwrite or --skip-existing"
    )


def run_command(command: CommandSpec, dry_run: bool) -> None:
    print(f"+ {shlex.join(command.argv)}")
    if dry_run:
        return
    subprocess.run(command.argv, check=True, cwd=REPO_ROOT)


def run_case_stages(args: argparse.Namespace, build_dir: Path, case: CaseConfig, stages: tuple[str, ...]) -> CaseResult:
    start = time.monotonic()

    check_case_inputs(case, stages)
    if not args.dry_run and check_output_policy(case, args.overwrite, args.skip_existing, stages):
        return CaseResult(case.name, "skipped", time.monotonic() - start)

    commands = build_commands(build_dir, case, stages, overwrite=args.overwrite)
    if not args.dry_run:
        check_tools(commands)

    for command in commands:
        run_command(command, args.dry_run)
    return CaseResult(case.name, "ok", time.monotonic() - start)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run runIPCSim and convertAnimation batches from JSON jobs.")
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG, help="JSON batch config file.")
    parser.add_argument("--job", action="append", default=[], help="Job name to run. Can be repeated.")
    parser.add_argument("--all-jobs", action="store_true", help="Run every job from the JSON config.")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without running them.")
    parser.add_argument("--overwrite", action="store_true", help="Allow runIPCSim to replace existing output folders.")
    parser.add_argument("--skip-existing", action="store_true", help="Skip cases whose runIPCSim output folder exists.")
    args = parser.parse_args()

    if args.overwrite and args.skip_existing:
        parser.error("--overwrite and --skip-existing are mutually exclusive")
    if args.all_jobs and args.job:
        parser.error("--all-jobs cannot be combined with --job")
    if not args.all_jobs and not args.job:
        parser.error("pass at least one --job, or pass --all-jobs explicitly")

    args.config = args.config if args.config.is_absolute() else REPO_ROOT / args.config
    return args


def print_summary(results: list[CaseResult]) -> None:
    print("== Summary ==")
    for result in results:
        print(f"{result.status:7} {result.name:24} {result.elapsed_seconds:.2f}s")


def main() -> int:
    args = parse_args()
    try:
        build_dir, cases, jobs = load_config(args.config)
        selected_jobs = list(jobs) if args.all_jobs else args.job
        unknown_jobs = [name for name in selected_jobs if name not in jobs]
        if unknown_jobs:
            raise ValueError(f"unknown jobs: {', '.join(unknown_jobs)}")

        results: list[CaseResult] = []
        for job_name in selected_jobs:
            job = jobs[job_name]
            print(f"== job {job.name} ==")
            for case_name in job.cases:
                case = cases[case_name]
                print(f"-- case {case.name} --")
                args.stages = job.stages
                results.append(run_case_stages(args, build_dir, case, job.stages))

        print_summary(results)
        return 0
    except (FileNotFoundError, FileExistsError, ValueError, json.JSONDecodeError, subprocess.CalledProcessError) as exc:
        print(f"error: {exc}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
