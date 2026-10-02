#!/usr/bin/env python3
"""Prepare, run and summarize the original full-length CUDA IPC examples.

No physical parameter, solver tolerance, frame count or dump interval is changed.
The runner owns output/newton-solves.csv; this tool reads that file and writes
separate summary.json and newton-iterations.csv artifacts in the case directory.
"""

from __future__ import annotations

import argparse
import array
import csv
import hashlib
import html
import json
import math
import os
from pathlib import Path
import re
import struct
import subprocess
import sys
import time
from collections import Counter


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DESTINATION = ROOT / "examples/generated/cuda-validation-20261002"
VOLUME_IPC_CASES = {
    "box-tet": ("box/box-tet-ipc.json", 0),
    "box-cubic": ("box/box-cubic-ipc.json", 1),
    "box-with-sphere-tet": ("box-with-sphere/box-with-sphere-tet-ipc.json", 1),
    "box-with-sphere-cubic": ("box-with-sphere/box-with-sphere-cubic-ipc.json", 0),
    "bunny-tet": ("bunny/bunny-tet-ipc.json", 0),
    "bunny-cubic": ("bunny/bunny-cubic-ipc.json", 1),
    "dragon-dynamic-tet": ("dragon-dyn/dragon-dynamic-tet-ipc.json", 1),
    "dragon-dynamic-cubic": ("dragon-dyn/dragon-dynamic-cubic-ipc.json", 0),
    "box-contact-static-tet": ("box-contact/box-contact-static-tet-ipc.json", 0),
    "box-contact-static-cubic": ("box-contact/box-contact-static-cubic-ipc.json", 1),
    "box-contact-dynamic-tet": ("box-contact/box-contact-dynamic-tet-ipc.json", 0),
}
SHELL_IPC_CASES = {
    "shell-hang-static": ("shell_hang/shell-static-ipc.json", 0),
    "shell-hang-dynamic": ("shell_hang/shell-dynamic-ipc.json", 0),
    "shirt-drop": ("shirt_drop/shirt-drop-dynamic-ipc.json", 1),
    "towel-grasp": ("towel_grasp/towel-grasp-dynamic-ipc.json", 1),
}
CASE_ALIASES = {"box": "box-tet", "box-with-sphere": "box-with-sphere-tet",
                "bunny": "bunny-tet", "dragon": "dragon-dynamic-tet"}
CASES = {**VOLUME_IPC_CASES, **SHELL_IPC_CASES, **{alias: VOLUME_IPC_CASES[name] for alias, name in CASE_ALIASES.items()}}
PRIMARY_CASES = ("box", "box-with-sphere", "bunny")
PATH_KEYS = {"tet-mesh", "surface-mesh", "initial-surface-mesh", "cubic-mesh", "filename", "trajectory"}
CUDA_STRATEGY_KEYS = {"cuda-pcg-batch-size", "cuda-pcg-tolerance", "cuda-exact-hessian-start-iteration", "cuda-stagnation-window", "cuda-stagnation-tolerance"}
NONPHYSICAL_KEYS = {"output", "simulation-backend", "cuda-device", "cuda-solver-verbose", "cuda-cubic-validity"} | CUDA_STRATEGY_KEYS
INT_FIELDS = {"frame", "newton_iterations", "pcg_iterations", "active_pairs", "status", "regularizations", "roundoff_steps"}
OPTIONAL_INT_FIELDS = {"exact_hessian_attempts", "exact_hessian_directions", "exact_hessian_fallbacks",
                       "substeps", "refinement_retries", "max_substep_newton_iterations",
                       "rejected_newton_solves", "input_chord_rejections", "macroframe_chord_certified", "trial_chord_rejections",
                       "residual_correction_attempts", "residual_correction_candidates", "residual_correction_steps", "friction_solves", "stagnation_triggers"}
FLOAT_FIELDS = {"gradient_inf", "tolerance", "energy", "linear_residual", "last_feasible_alpha", "last_accepted_alpha"}
NUMBER = r"[-+]?(?:\d*\.?\d+(?:[eE][-+]?\d+)?|inf|nan)"
ITERATION = re.compile(rf"CUDA Newton iter=(\d+); E=({NUMBER}); \|\|g\|\|_inf=({NUMBER})")
FRAME = re.compile(r"CUDA frame (\d+):")
SOLVE = re.compile(r"CUDA solve frame (\d+)(?:, refinement=(\d+), substep=(\d+)/(\d+))?(?:, friction-iteration=(\d+))?")
SOLVE_WORK = re.compile(rf"CUDA Newton status=(-?\d+); stop=([^;]+); \|\|g\|\|_inf=({NUMBER}); Newton iterations=(\d+); PCG iterations=(\d+)")
COMPLETE = re.compile(rf"CUDA IPC complete: frames=(\d+), peak-contacts=(\d+), Newton=(\d+), PCG=(\d+), nonconverged=(\d+), elapsed=({NUMBER}) s")
UNCOMMITTED_FAILURE = re.compile(r"CUDA frame (\d+) could not converge with certified physical motion after (\d+) refinements; no state was committed")


def digest(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def absolute_paths(value, parent: Path):
    if isinstance(value, list):
        return [absolute_paths(item, parent) for item in value]
    if isinstance(value, dict):
        return {
            key: str((parent / item).resolve()) if key in PATH_KEYS and isinstance(item, str)
            else absolute_paths(item, parent)
            for key, item in value.items()
        }
    return value


def write_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def original_config(name: str):
    source = ROOT / "examples/configs" / ("shell" if name in SHELL_IPC_CASES else "volume") / CASES[name][0]
    return source, absolute_paths(json.loads(source.read_text()), source.parent)


def unchanged_physics(name: str, config: dict) -> bool:
    _, original = original_config(name)
    return {k: v for k, v in config.items() if k not in NONPHYSICAL_KEYS} == {
        k: v for k, v in original.items() if k not in NONPHYSICAL_KEYS
    }


def expected_frames(config: dict) -> int:
    return 1 if config.get("sim-type") == "static" else config["num-timestep"]


def prepare(destination: Path, selected: list[str]) -> None:
    manifest_path = destination / "manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {"schema_version": 1, "cases": {}}
    for name in selected:
        source, config = original_config(name)
        directory = destination / name
        config.update({"output": str(directory / "output"), "simulation-backend": "cuda",
                       "cuda-device": CASES[name][1], "cuda-solver-verbose": 2})
        if "cubic-mesh" in config:
            config["cuda-cubic-validity"] = "gauss"
        assert unchanged_physics(name, config)
        filename = directory / "config.json"
        if filename.exists() and json.loads(filename.read_text()) != config:
            raise ValueError(f"Refusing to overwrite a different existing config: {filename}")
        if not filename.exists():
            write_json(filename, config)
        manifest["cases"][name] = {
            "source": str(source), "source_sha256": digest(source), "config": str(filename),
            "config_sha256": digest(filename), "output": config["output"], "device": config["cuda-device"],
            "num_timestep": expected_frames(config), "timestep": config.get("timestep"),
            "solver_eps": config["solver-eps"], "planned": True,
            "physics_unchanged": True,
            "cubic_validity": config.get("cuda-cubic-validity"),
            "not_run_reason": None,
        }
    write_json(manifest_path, manifest)


def run_case(destination: Path, name: str, binary: Path, threads: int, cpu_set: str | None) -> int:
    directory = destination / name
    config_path = directory / "config.json"
    config = json.loads(config_path.read_text())
    if not unchanged_physics(name, config):
        raise ValueError(f"Physical parameters differ from the source example: {name}")
    if (Path(config["output"]) / "newton-solves.csv").exists() or (directory / "run.log").exists():
        raise ValueError(f"Refusing to overwrite an existing run: {directory}")
    environment = os.environ.copy()
    for key in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
        environment[key] = str(threads)
    command = [str(binary), str(config_path)]
    if cpu_set:
        # Affinity also constrains oneTBB, whose current runner has no thread
        # limit option and does not consume a TBB_NUM_THREADS environment key.
        command = ["taskset", "--cpu-list", cpu_set, *command]
    record = {"command": command, "binary_sha256": digest(binary), "config_sha256": digest(config_path),
              "source_config_sha256": digest(original_config(name)[0]), "thread_environment": {
                  key: environment[key] for key in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS")},
              "cpu_affinity": cpu_set, "started_unix": time.time()}
    write_json(directory / "run-record.json", record)
    started = time.monotonic()
    with (directory / "run.log").open("w") as log:
        result = subprocess.run(command, cwd=ROOT, env=environment, stdout=log, stderr=subprocess.STDOUT, check=False)
    record.update({"exit_code": result.returncode, "finished_unix": time.time(), "wall_seconds": time.monotonic() - started})
    write_json(directory / "run-record.json", record)
    return result.returncode


def read_solver_rows(path: Path) -> tuple[list[dict], int]:
    if not path.exists():
        return [], 0
    # Snapshot complete lines: the simulator may currently be appending a row.
    snapshot = path.read_text()
    lines = snapshot.splitlines(keepends=True)
    incomplete = int(bool(lines) and not lines[-1].endswith("\n"))
    if incomplete:
        lines.pop()
    rows = list(csv.DictReader(lines))
    for row in rows:
        for key in INT_FIELDS:
            row[key] = int(row[key])
        for key in OPTIONAL_INT_FIELDS:
            if key in row:
                row[key] = int(row[key])
        for key in FLOAT_FIELDS:
            row[key] = float(row[key])
    return rows, incomplete


def iteration_rows(log: str) -> list[dict]:
    output = []
    pending = []
    context = {"frame": None, "refinement": None, "substep": None, "substeps": None, "friction_iteration": None}
    for line in log.splitlines():
        if match := SOLVE.search(line):
            context = {"frame": int(match[1]), "refinement": int(match[2]) if match[2] else None,
                       "substep": int(match[3]) - 1 if match[3] else None,
                       "substeps": int(match[4]) if match[4] else None,
                       "friction_iteration": int(match[5]) if match[5] else 0}
        if match := ITERATION.search(line):
            pending.append({**context, "iteration": int(match[1]), "energy_logged": match[2],
                            "gradient_inf_logged": match[3]})
        if match := FRAME.search(line):
            for row in pending:
                if row["frame"] is None:
                    row["frame"] = int(match[1])
            output.extend(pending)
            pending = []
            context = {"frame": None, "refinement": None, "substep": None, "substeps": None, "friction_iteration": None}
    # Preserve diagnostics from an interrupted final solve without assigning
    # them to an unconfirmed completed frame.
    output.extend(pending)
    return output


def inspect_states(output: Path, expected: int) -> dict:
    states = sorted(output.glob("deform*.u"))
    indices = [int(path.stem.removeprefix("deform")) for path in states]
    rows = None
    failures = []
    for path in states:
        with path.open("rb") as stream:
            header = stream.read(12)
            if len(header) != 12:
                failures.append(f"Invalid header: {path.name}")
                continue
            count, columns, scalar_bytes = struct.unpack("=iii", header)
            data = array.array("d")
            try:
                data.frombytes(stream.read())
            except ValueError:
                failures.append(f"Invalid double data: {path.name}")
                continue
        if count <= 0 or count % 3 or columns != 3 or scalar_bytes != 8 or len(data) != count * columns:
            failures.append(f"Invalid dimensions: {path.name}")
        elif (rows is not None and rows != count) or not all(map(math.isfinite, data)):
            failures.append(f"Nonfinite or inconsistent state: {path.name}")
        rows = count
    return {"state_files": len(states), "expected_state_files": expected, "dofs": rows,
            "state_frames_complete": indices == list(range(expected)), "all_saved_states_finite": not failures,
            "state_errors": failures}


def solver_work_rows(log: str) -> list[dict]:
    rows = []
    context = {"frame": None, "refinement": None, "substep": None, "substeps": None, "friction_iteration": None}
    for line in log.splitlines():
        if match := SOLVE.search(line):
            context = {"frame": int(match[1]), "refinement": int(match[2]) if match[2] else None,
                       "substep": int(match[3]) - 1 if match[3] else None,
                       "substeps": int(match[4]) if match[4] else None,
                       "friction_iteration": int(match[5]) if match[5] else 0}
        if match := SOLVE_WORK.search(line):
            rows.append({**context, "status": int(match[1]), "stop": match[2],
                         "gradient_inf_logged": match[3], "newton_iterations": int(match[4]),
                         "pcg_iterations": int(match[5])})
    return rows


def summarize(destination: Path, name: str, check_states: bool) -> dict:
    directory = destination / name
    config = json.loads((directory / "config.json").read_text())
    output = Path(config["output"])
    log_candidates = (output / "runIPCSim.log", directory / "run.log", directory / "console.log")
    logfile = next((path for path in log_candidates if path.exists()), log_candidates[0])
    log = logfile.read_text(errors="replace") if logfile.exists() else ""
    rows, incomplete_rows = read_solver_rows(output / "newton-solves.csv")
    completed = COMPLETE.findall(log)
    uncommitted_failures = UNCOMMITTED_FAILURE.findall(log)
    expected = expected_frames(config)
    failures = []
    for row in rows:
        finite = all(math.isfinite(row[key]) for key in FLOAT_FIELDS)
        if not finite or row["status"] != 0 or not row["gradient_inf"] < row["tolerance"] or row["tolerance"] != config["solver-eps"]:
            failures.append({key: row[key] if not isinstance(row[key], float) or math.isfinite(row[key]) else str(row[key])
                             for key in ("frame", "status", "stop", "gradient_inf", "tolerance")})
    frame_sequence = [row["frame"] for row in rows]
    worst = max(rows, key=lambda row: row["gradient_inf"] / row["tolerance"], default=None)
    result = {
        "case": name, "source": str(original_config(name)[0]), "config": str(directory / "config.json"),
        "simulation_type": config.get("sim-type"),
        "physics_unchanged": unchanged_physics(name, config), "expected_frames": expected,
        "uncommitted_failed_frame": int(uncommitted_failures[-1][0]) if uncommitted_failures else None,
        "uncommitted_failed_refinements": int(uncommitted_failures[-1][1]) if uncommitted_failures else None,
        "cuda_solver_strategy": {key: config[key] for key in sorted(CUDA_STRATEGY_KEYS) if key in config},
        "configured_cubic_validity": config.get("cuda-cubic-validity"),
        "frames_recorded": len(rows), "frame_sequence_complete": frame_sequence == list(range(expected)),
        "incomplete_csv_rows_ignored": incomplete_rows,
        "completion_marker_present": bool(completed), "solver_tolerance": config["solver-eps"],
        "all_frames_converged": bool(rows) and not failures and len(rows) == expected,
        "failed_frames": failures, "stop_counts": dict(Counter(row["stop"] for row in rows)),
        "failed_stop_counts": dict(Counter(row["stop"] for row in failures)),
        "first_failed_frame": failures[0]["frame"] if failures else None,
        "worst_gradient_frame": worst["frame"] if worst else None,
        "worst_gradient_inf": worst["gradient_inf"] if worst and math.isfinite(worst["gradient_inf"]) else None,
        "worst_gradient_tolerance_ratio": worst["gradient_inf"] / worst["tolerance"]
            if worst and math.isfinite(worst["gradient_inf"]) and worst["tolerance"] > 0 else None,
        "peak_active_pairs": max((row["active_pairs"] for row in rows), default=0),
        "contact_expected": name not in ("shell-hang-static", "shell-hang-dynamic"),
        "first_contact_frame": next((row["frame"] for row in rows if row["active_pairs"] > 0), None),
        "contact_frames": sum(row["active_pairs"] > 0 for row in rows),
        "total_newton_iterations": sum(row["newton_iterations"] for row in rows),
        "total_pcg_iterations": sum(row["pcg_iterations"] for row in rows),
        "frame_iteration_totals_scope": "accepted physical solves, or accepted static equilibrium solves; excludes rejected attempts",
        "max_newton_iterations": max((row["newton_iterations"] for row in rows), default=0),
        "max_newton_iterations_per_physical_solve": max((row.get("max_substep_newton_iterations", row["newton_iterations"]) for row in rows), default=0),
        "max_pcg_iterations": max((row["pcg_iterations"] for row in rows), default=0),
        "total_regularizations": sum(row["regularizations"] for row in rows),
        "total_roundoff_steps": sum(row["roundoff_steps"] for row in rows),
        "total_exact_hessian_attempts": sum(row.get("exact_hessian_attempts", 0) for row in rows),
        "total_exact_hessian_directions": sum(row.get("exact_hessian_directions", 0) for row in rows),
        "total_exact_hessian_fallbacks": sum(row.get("exact_hessian_fallbacks", 0) for row in rows),
        "accepted_physical_substeps": sum(row.get("substeps", 1) for row in rows) if config.get("sim-type") != "static" else 0,
        "static_terminal_records": len(rows) if config.get("sim-type") == "static" else 0,
        "refined_frames": sum(row.get("substeps", 1) > 1 for row in rows),
        "max_substeps_per_frame": max((row.get("substeps", 1) for row in rows), default=0),
        "rejected_newton_solves": sum(row.get("rejected_newton_solves", 0) for row in rows),
        "input_chord_rejections": sum(row.get("input_chord_rejections", 0) for row in rows),
        "trial_chord_rejections": sum(row.get("trial_chord_rejections", 0) for row in rows),
        "accepted_stagnation_triggers": sum(row.get("stagnation_triggers", 0) for row in rows),
        "uncertified_coarse_chords": sum(row.get("macroframe_chord_certified", 1) == 0 for row in rows),
        "friction_solves": sum(row.get("friction_solves", 0) for row in rows),
        "residual_correction_attempts": sum(row.get("residual_correction_attempts", 0) for row in rows),
        "residual_correction_candidates": sum(row.get("residual_correction_candidates", 0) for row in rows),
        "residual_correction_steps": sum(row.get("residual_correction_steps", 0) for row in rows),
        "runner_csv": str(output / "newton-solves.csv"), "log": str(logfile),
    }
    if completed:
        frames, contacts, newton, pcg, unconverged, elapsed = completed[-1]
        result["runner_completion"] = {"frames": int(frames), "peak_contacts": int(contacts), "newton": int(newton),
                                       "pcg": int(pcg), "nonconverged": int(unconverged), "elapsed_seconds": float(elapsed)}
    record = directory / "run-record.json"
    if record.exists():
        result["run_record"] = json.loads(record.read_text())
    execution_path = directory / "execution.json"
    result["execution_record"] = str(execution_path) if execution_path.exists() else str(record) if record.exists() else None
    status_path = directory / "run-status.json"
    if status_path.exists():
        result["run_status"] = json.loads(status_path.read_text())
    result["exit_code"] = result.get("run_status", result.get("run_record", {})).get("exit_code")
    result["completed"] = result["frame_sequence_complete"] and bool(completed) and result["exit_code"] == 0
    collision_path = directory / "collision-validation.json"
    collision = json.loads(collision_path.read_text()) if collision_path.exists() else None
    collision_execution_path = directory / "collision-execution.json"
    collision_execution = json.loads(collision_execution_path.read_text()) if collision_execution_path.exists() else {}
    result["collision_execution_record"] = str(collision_execution_path) if collision_execution else None
    result["collision_checks_provisional"] = bool(collision_execution.get("checker_results_provisional"))
    result["collision_report"] = str(collision_path) if collision is not None else None
    result["minimum_tet_determinant"] = collision.get("minimum_tet_determinant") if collision is not None else None
    result["minimum_element_determinant"] = collision.get("minimum_element_determinant",
        collision.get("minimum_tet_determinant")) if collision is not None else None
    result["minimum_triangle_area_squared_ratio"] = collision.get("minimum_triangle_area_squared_ratio") if collision is not None else None
    result["element_validity_criterion"] = collision.get("element_validity_criterion",
        "whole-cell" if collision.get("element_type") == "cubic" else
        "positive-triangle-area" if collision.get("element_type") == "shell" else "positive-tet-determinant") if collision is not None else None
    whole_cell = collision.get("whole_cell_diagnostics", {}) if collision else {}
    strict_volume_pass = bool(collision and collision.get("passed") and collision.get("element_type") in ("tet", "tetrahedral", "cubic") and not whole_cell)
    result["confirmed_inverted_element_count"] = whole_cell.get("confirmed_inverted_element_count", 0 if strict_volume_pass else None)
    result["uncertified_without_negative_witness_element_count"] = whole_cell.get("uncertified_without_negative_witness_element_count", 0 if strict_volume_pass else None)
    result["collision_validation_passed"] = bool(collision is not None and collision.get("passed") is True
        and collision.get("frames_checked") == collision.get("frames_expected") == expected
        and collision.get("frames_failed") == 0
        and collision.get("swept_checks_performed", collision.get("swept_checks_requested")) is True
        and collision.get("first_frame", 0) == 0 and collision.get("full_frame_range_requested", True) is True
        and not result["collision_checks_provisional"])
    if collision is not None:
        result["collision_summary"] = {key: collision.get(key) for key in (
            "passed", "frames_expected", "frames_checked", "frames_failed", "self_edge_triangle_intersections",
            "external_edge_triangle_intersections", "minimum_tet_determinant", "minimum_element_determinant",
            "element_type", "element_validity_criterion", "whole_cell_diagnostics", "swept_checks_requested", "minimum_is_certified_lower_bound", "elapsed_seconds",
            "minimum_scope", "minimum_triangle_area_squared_ratio", "physical_substeps", "static_newton_path", "trajectory_path", "swept_checks_performed", "coarse_chord_diagnostics",
            "uncertified_sweep_segments", "sweep_midpoint_self_edge_triangle_intersections",
            "sweep_midpoint_external_edge_triangle_intersections")}
        result["collision_summary"]["minimum_element_determinant"] = result["minimum_element_determinant"]
    static_path = output / "static-newton-path.json"
    result["static_newton_path_manifest"] = str(static_path) if static_path.exists() else None
    iterations = iteration_rows(log)
    work = solver_work_rows(log)
    work_path = directory / "newton-attempt-work.csv"
    result["solver_work_csv"] = str(work_path) if work else None
    result["completed_attempts_logged"] = len(work)
    result["total_completed_attempt_newton_iterations"] = sum(row["newton_iterations"] for row in work) if work else None
    result["total_completed_attempt_pcg_iterations"] = sum(row["pcg_iterations"] for row in work) if work else None
    result["attempt_work_scope"] = "all completed Newton attempts in the verbose log, including rejected attempts; excludes any interrupted in-progress solve"
    if work:
        with work_path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=tuple(work[0]))
            writer.writeheader()
            writer.writerows(work)
    attempts_path = output / "newton-attempts.csv"
    result["attempts_csv"] = str(attempts_path) if attempts_path.exists() else None
    substeps_path = output / "newton-substeps.csv"
    result["substeps_csv"] = str(substeps_path) if substeps_path.exists() else None
    if attempts_path.exists():
        attempts = list(csv.DictReader(attempts_path.open()))
        accepted = {(int(row["frame"]), int(row["refinement"]), int(row["substep"]), int(row.get("friction_iteration", 0))): row["accepted"] for row in attempts}
        for row in iterations:
            row["accepted"] = accepted.get((row["frame"], row["refinement"], row["substep"], row["friction_iteration"]))
    with (directory / "newton-iterations.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=("frame", "refinement", "substep", "substeps", "friction_iteration", "accepted", "iteration", "energy_logged", "gradient_inf_logged"))
        writer.writeheader()
        writer.writerows(iterations)
    result["iteration_rows"] = len(iterations)
    result["iteration_precision_note"] = "Iteration values come from verbose text; terminal frame values use the full-precision runner CSV."
    if check_states and result["completed"]:
        # Immutable completed trajectories can be checked once. A changed file
        # count, size, or modification time invalidates the cached full scan.
        states = sorted(output.glob("deform*.u"))
        signature = [[path.name, path.stat().st_size, path.stat().st_mtime_ns] for path in states]
        cache_path = directory / "state-validation.json"
        cache = json.loads(cache_path.read_text()) if cache_path.exists() else {}
        if cache.get("signature") != signature or cache.get("expected_frames") != expected:
            cache = {"signature": signature, "expected_frames": expected, "result": inspect_states(output, expected)}
            write_json(cache_path, cache)
        result.update(cache["result"])
        result["state_check"] = "completed trajectory checked"
    else:
        result["state_files"] = sum(1 for _ in output.glob("deform*.u"))
        result["state_check"] = "pending completion" if not result["completed"] else "skipped"
    if rows:
        result["convergence_plot"] = str(make_plot(directory, name, rows))
    # Prefix diagnostics have separate, explicitly labeled report links.
    # They must not stand in for the complete case's simulation video.
    videos = [path for path in sorted(directory.glob("*.mp4")) + sorted(output.glob("*.mp4"))
              if path.name != "partial-prefix.mp4"]
    result["video"] = str(videos[0]) if videos else None
    preview = directory / "preview.png"
    result["preview"] = str(preview) if preview.exists() else None
    result["validation_passed"] = bool(result["physics_unchanged"] and result["frame_sequence_complete"] and
        result["completion_marker_present"] and result["all_frames_converged"] and
        (not result["contact_expected"] or result["peak_active_pairs"] > 0) and
        result.get("state_frames_complete", False) and result.get("all_saved_states_finite", False) and
        result["exit_code"] == 0 and result["collision_validation_passed"])
    if result["validation_passed"]:
        result["assessment"] = "passed"
    elif result["exit_code"] not in (None, 0):
        result["assessment"] = "stopped" if result["exit_code"] < 0 else "process failed"
    elif failures or (collision is not None and collision.get("passed") is False):
        result["assessment"] = "failed" if result["completed"] else "running; failures observed"
    else:
        result["assessment"] = "checks pending" if result["completed"] else "running" if rows else "not run"
    write_json(directory / "summary.json", result)
    return result


def make_plot(directory: Path, name: str, rows: list[dict]) -> Path:
    os.environ.setdefault("MPLCONFIGDIR", str(directory.parent / ".matplotlib"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    frames = [row["frame"] for row in rows]
    ratio = [max(row["gradient_inf"] / row["tolerance"], 1e-16) for row in rows]
    fig, axes = plt.subplots(3, 1, sharex=True, figsize=(11, 7), constrained_layout=True)
    axes[0].semilogy(frames, ratio, linewidth=0.85, color="#166a9a")
    axes[0].axhline(1, color="#ba3333", linestyle="--", linewidth=1, label="Required: gradient / tolerance < 1")
    failed = [row for row in rows if row["status"] != 0 or row["gradient_inf"] >= row["tolerance"]]
    axes[0].scatter([row["frame"] for row in failed],
                    [max(row["gradient_inf"] / row["tolerance"], 1e-16) for row in failed],
                    s=7, color="#ba3333", label="Nonconverged frames", zorder=3)
    axes[0].set_ylabel("Gradient / tolerance")
    axes[0].legend(loc="best", fontsize=8)
    axes[0].set_title(f"{name}: {len(rows):,} recorded frames; {len(failed):,} nonconverged")
    axes[1].plot(frames, [row["newton_iterations"] for row in rows], color="#256f55", linewidth=0.85)
    axes[1].set_ylabel("Newton iterations\n(accepted solves)")
    axes[2].plot(frames, [row["pcg_iterations"] for row in rows], color="#7d55a0", linewidth=0.85)
    axes[2].set_ylabel("PCG iterations\n(accepted solves)")
    axes[2].set_xlabel("Simulation frame (zero based)")
    for axis in axes:
        axis.grid(alpha=0.18)
    filename = directory / "newton-convergence.png"
    fig.savefig(filename, dpi=150)
    plt.close(fig)
    return filename


def write_report(destination: Path, results: list[dict]) -> None:
    write_json(destination / "summary.json", {"schema_version": 1, "cases": results})
    fields = ("case", "assessment", "expected_frames", "frames_recorded", "failed_frame_count",
              "uncommitted_failed_frame", "uncommitted_failed_refinements",
              "solver_tolerance", "worst_gradient_inf", "worst_gradient_tolerance_ratio", "first_failed_frame",
              "worst_gradient_frame", "peak_active_pairs", "contact_frames", "total_newton_iterations",
              "total_pcg_iterations", "total_roundoff_steps", "total_exact_hessian_attempts",
              "total_completed_attempt_newton_iterations", "total_completed_attempt_pcg_iterations",
              "total_exact_hessian_directions", "total_exact_hessian_fallbacks", "accepted_physical_substeps",
              "max_newton_iterations_per_physical_solve",
              "refined_frames", "max_substeps_per_frame", "rejected_newton_solves", "input_chord_rejections", "trial_chord_rejections", "uncertified_coarse_chords",
              "residual_correction_attempts", "residual_correction_candidates", "residual_correction_steps", "friction_solves",
              "accepted_stagnation_triggers", "minimum_tet_determinant", "minimum_element_determinant", "minimum_triangle_area_squared_ratio",
              "configured_cubic_validity", "element_validity_criterion", "confirmed_inverted_element_count",
              "uncertified_without_negative_witness_element_count", "collision_validation_passed", "validation_passed")
    with (destination / "summary-all-cases.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for result in results:
            writer.writerow({key: len(result["failed_frames"]) if key == "failed_frame_count" else result.get(key) for key in fields})

    def relative(path):
        return html.escape(os.path.relpath(path, destination), quote=True)

    def link(path, label):
        return f'<a href="{relative(path)}">{html.escape(label)}</a>' if path and Path(path).exists() else html.escape(label + " pending")

    table_rows, cards = [], []
    for result in results:
        name = html.escape(result["case"])
        failed = len(result["failed_frames"])
        cls = "pass" if result["validation_passed"] else "fail" if failed or result.get("collision_summary", {}).get("passed") is False or result["exit_code"] not in (None, 0) else "pending"
        collisions = "Passed" if result["collision_validation_passed"] else "Failed" if result.get("collision_summary", {}).get("passed") is False else "Provisional" if result["collision_checks_provisional"] else "Pending"
        ratio = result["worst_gradient_tolerance_ratio"]
        ratio_text = f"{ratio:.6g}" if ratio is not None else "—"
        table_rows.append(f'<tr><td><a href="#{name}">{name}</a></td><td class="{cls}">{html.escape(result["assessment"])}</td>'
                          f'<td>{result["frames_recorded"]:,} / {result["expected_frames"]:,}</td><td>{failed:,}</td>'
                          f'<td>{result["solver_tolerance"]:g}</td><td>{ratio_text}</td><td>{result["peak_active_pairs"]:,}</td><td>{collisions}</td></tr>')
        directory = destination / result["case"]
        links = " · ".join((link(directory / "config.json", "Config"), link(directory / "summary.json", "JSON summary"),
                            link(result["execution_record"], "Execution provenance"),
                            link(result["runner_csv"], "Newton solves CSV"), link(directory / "newton-iterations.csv", "Newton iterations CSV"),
                            link(result["attempts_csv"], "Accepted and rejected attempts"), link(result["substeps_csv"], "Physical substeps / static terminal state"),
                            link(result["log"], "Solver log"), link(directory / "console.log", "Console log"),
                            link(result["collision_report"], "CPU trajectory report"), link(directory / "collision-frames.csv", "Trajectory checks CSV"),
                            link(result["collision_execution_record"], "Checker provenance"),
                            link(result["video"], "Video")))
        if result.get("static_newton_path_manifest"):
            links += " · " + link(result["static_newton_path_manifest"], "Static Newton path")
        if result.get("solver_work_csv"):
            links += " · " + link(result["solver_work_csv"], "Newton/PCG work including retries")
        for filename, label in (("render-record.json", "Render provenance"),
                                ("partial-prefix.mp4", "Validated prefix video (incomplete run)"),
                                ("partial-prefix-render-record.json", "Prefix video provenance"),
                                ("independent-final-audit.json", "Independent final artifact audit"),
                                ("final-run-artifact-audit.json", "Final run and manifest audit"),
                                ("frozen-build-provenance.json", "Frozen build provenance"),
                                ("frozen-build-record.json", "Build record"),
                                ("frozen-source.tar.gz", "Source snapshot"),
                                ("previous-attempt.html", "Prior run and diagnosis"),
                                ("current-run-diagnosis.html", "Current cubic investigation"),
                                ("criterion-decision.json", "Selected validity criterion"),
                                ("verification.json", "Build verification"),
                                ("whole-cell-diagnostics.json", "Separate whole-cell diagnostics"),
                                ("inversion-timeline/inversion-counts.csv", "Inverted-element counts by frame"),
                                ("inversion-timeline/inversion-timeline.png", "Inversion timeline"),
                                ("inversion-timeline/inversion-audit.json", "Inversion count audit"),
                                ("frame413-plateau-audit.json", "Contact plateau diagnostic"),
                                ("interruption-record.json", "Interrupted run record"),
                                ("material-feasibility-audit.json", "Final GPU material guard audit"),
                                ("material-objective-comparison.json", "CPU/GPU material comparison"),
                                ("element-validity-history.png", "Element Jacobian history"),
                                ("element-validity-history.csv", "Element Jacobian samples"),
                                ("render-metadata-correction.json", "Static metadata correction"),
                                ("numerical-strategy-selection.json", "Solver strategy evidence"),
                                ("resource-scheduling.json", "GPU scheduling")):
            if (directory / filename).exists():
                links += " · " + link(directory / filename, label)
        plot = f'<img class="plot" src="{relative(result["convergence_plot"])}" alt="{name} Newton convergence and iteration counts">' if result.get("convergence_plot") else ""
        video = f'<video controls preload="metadata" src="{relative(result["video"])}"></video>' if result["video"] else '<p class="muted">Video pending.</p>'
        failure = (f'First nonconverged frame: {result["first_failed_frame"]}; worst gradient: {result["worst_gradient_inf"]:.9g} '
                   f'at frame {result["worst_gradient_frame"]}.' if failed and result["worst_gradient_inf"] is not None else "No nonconverged frames in the recorded CSV.")
        if result.get("uncommitted_failed_frame") is not None:
            failure += (f' Frame {result["uncommitted_failed_frame"]} exhausted its refinement budget '
                        f'(maximum level {result["uncommitted_failed_refinements"]}) without a converged, certified result; '
                        'no state was committed for that frame.')
        minimum_determinant = result["minimum_element_determinant"]
        determinant_label = ("Certified lower bound on det(F)" if
            result.get("collision_summary", {}).get("minimum_is_certified_lower_bound") else
            "Uncertified Jacobian diagnostic" if result.get("collision_summary", {}).get("element_type") == "cubic" else
            "Minimum det(F)")
        if result.get("collision_summary", {}).get("element_type") == "shell":
            minimum_determinant = result["minimum_triangle_area_squared_ratio"]
            determinant_label = ("Certified lower bound on squared area / squared rest area" if
                result.get("collision_summary", {}).get("minimum_is_certified_lower_bound") else "Uncertified shell-area diagnostic")
        determinant_note = (f' {determinant_label}: {minimum_determinant:.9g}.'
                            if minimum_determinant is not None else "")
        collision_summary = result.get("collision_summary", {})
        if collision_summary.get("element_validity_criterion"):
            determinant_note += f' Element validity criterion: {html.escape(collision_summary["element_validity_criterion"])}.'
        whole_cell = collision_summary.get("whole_cell_diagnostics")
        if whole_cell and whole_cell.get("required_for_validation") is False:
            determinant_note += (' Whole-cell positivity is a separate diagnostic under the CPU-compatible criterion; '
                                 'see the independent check report for inverted or uncertified cells.')
        if result.get("confirmed_inverted_element_count") is not None:
            determinant_note += f' Confirmed inverted elements: {result["confirmed_inverted_element_count"]:,} unique elements.'
        if result.get("uncertified_without_negative_witness_element_count"):
            determinant_note += (f' Elements with inconclusive whole-cell checks and no negative witness: '
                                 f'{result["uncertified_without_negative_witness_element_count"]:,}; these are not counted as confirmed inversions.')
        sweep_note = (" Piecewise linear sweeps were checked." if
            collision_summary.get("swept_checks_performed", collision_summary.get("swept_checks_requested")) else "")
        if collision_summary.get("trajectory_path") == "recorded_physical_substeps":
            if result.get("simulation_type") == "static":
                sweep_note += (" The static initial-to-equilibrium straight path was checked. "
                               "This geometric feasibility path has no physical timing.")
            else:
                sweep_note += (f' The recorded motion follows {result["accepted_physical_substeps"]:,} physical substeps. '
                               f'{result["uncertified_coarse_chords"]:,} direct coarse-frame shortcuts were not certified; '
                               'use the recorded substep path for interpolation.')
        if collision_summary.get("trajectory_path") == "recorded_static_newton_iterates":
            sweep_note += (" Static geometry follows the recorded accepted Newton path. "
                           "Intermediate knots are optimization iterates; every terminal equilibrium solve "
                           "must meet the original tolerance. This path has no physical timing.")
        cpu_replay = directory / "repro-frame1106-cpu/output/runIPCSim.log"
        comparison = ""
        if result["case"] == "box" and cpu_replay.exists():
            replay_text = cpu_replay.read_text(errors="replace")
            cpu_result = re.findall(r"Newton result: stop=([^,]+), iteration=(\d+), gradient_inf=([^,]+), tolerance=([^,]+), converged=(\w+)", replay_text)
            gpu_result = next((row for row in result["failed_frames"] if row["frame"] == 1106), None)
            if cpu_result and gpu_result:
                stop, iteration, gradient, tolerance, converged = cpu_result[-1]
                comparison = (f'<p>Checkpoint replay diagnostic: restarting the CPU solver from the saved GPU frame 1105 state '
                              f'and solving frame 1106 at the original tolerance {html.escape(tolerance)} reached '
                              f'gradient {html.escape(gradient)} in {iteration} iterations ({html.escape(stop)}). '
                              f'The GPU trajectory recorded gradient {gpu_result["gradient_inf"]:.12g} for that solve. '
                              f'This is a single checkpoint replay, not a complete CPU trajectory. '
                              f'{link(cpu_replay, "CPU replay log")}.</p>')
        strategy_note = (f'<p>Explicit CUDA solver strategy: {html.escape(json.dumps(result["cuda_solver_strategy"], sort_keys=True))}.</p>'
                         if result.get("cuda_solver_strategy") else "")
        cards.append(f'<section id="{name}"><h2>{name}</h2><p>{links}</p><p>{failure} '
                     f'Stops: {html.escape(json.dumps(result["stop_counts"]))}. CPU trajectory checks: {collisions.lower()}. '
                     f'Saved states: {html.escape(result["state_check"])}.{determinant_note}{sweep_note}</p>{strategy_note}{comparison}{plot}{video}</section>')
    document = f'''<!doctype html>
<html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>CUDA IPC example validation</title><style>
body{{font:15px/1.5 system-ui,sans-serif;color:#243447;background:#f4f6f8;max-width:1200px;margin:auto;padding:32px}}
h1,h2{{line-height:1.2;color:#12253a}}a{{color:#145f95}}table{{border-collapse:collapse;width:100%;font-variant-numeric:tabular-nums}}
th,td{{padding:10px 12px;text-align:left;border-bottom:1px solid #dae1e8}}th{{background:#eaf0f5}}section{{background:white;border:1px solid #dae1e8;border-radius:10px;padding:24px;margin:24px 0}}
.pass{{color:#136f43;font-weight:600}}.fail{{color:#a52b32;font-weight:600}}.pending,.muted{{color:#657285}}.plot,video{{display:block;width:100%;height:auto;margin-top:20px}}.scroll{{overflow:auto}}
</style><h1>CUDA IPC example validation</h1>
<p>The original physics, stopping tolerances, and configured frame counts are retained. A pass requires every accepted physical solve to converge below its original tolerance, complete finite saved states, successful process completion, and passing independent CPU checks for collisions and the recorded element validity criterion. For cubic runs explicitly using the CPU-compatible Gauss criterion, whole-cell corner inversions are reported separately and are not contact failures. Historical strict whole-cell checks remain labeled as such. Adaptive runs follow every recorded physical substep; a direct shortcut between coarse frame endpoints can intersect even when that actual path is clear. Those shortcuts are reported separately. Static runs with an explicit Newton path validate each accepted optimization transition and require convergence of every terminal equilibrium solve; these iterates have no physical timing. Process exit 0 alone does not establish convergence.</p>
<p>{link(destination / 'summary-all-cases.csv', 'All cases CSV')} · {link(destination / 'summary.json', 'All cases JSON')} · {link(destination / 'manifest.json', 'Input manifest')} · {link(destination / 'validator-controls.json', 'Validator positive and negative controls')}</p>
{f'<p>{link(destination / "hardware.json", "CUDA validation hardware")}</p>' if (destination / 'hardware.json').exists() else ''}
<div class="scroll"><table><thead><tr><th>Case</th><th>Assessment</th><th>Frames</th><th>Nonconverged</th><th>Tolerance</th><th>Worst gradient / tolerance</th><th>Peak pairs</th><th>CPU trajectory checks</th></tr></thead><tbody>{''.join(table_rows)}</tbody></table></div>
<p class="muted">Partial runs are snapshots. Contact counts include self contact and contact with obstacles. Frame indices are zero based. Frame plots count accepted solves; the Newton/PCG work CSV includes completed rejected attempts as well. Per-case execution provenance records binary and config hashes and run conditions; historical source_sha256 fields refer to the original example config. Frozen code snapshots are linked separately. Overlapping jobs mean elapsed times are not isolated performance measurements.</p>
{''.join(cards)}</html>'''
    (destination / "index.html").write_text(document)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "run", "summarize"), nargs="?")
    parser.add_argument("--summarize", action="store_true", help="Alias for the summarize action.")
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument("--cases", nargs="+", choices=tuple(CASES), help="Selected canonical cases or legacy aliases; default: box, box-with-sphere, bunny.")
    selection.add_argument("--volume-ipc", action="store_true", help="Select all eleven canonical tetrahedral and cubic volume IPC examples.")
    selection.add_argument("--shell-ipc", action="store_true", help="Select the four canonical Koiter shell IPC examples.")
    parser.add_argument("--destination", "--output-dir", type=Path, default=DEFAULT_DESTINATION)
    parser.add_argument("--binary", type=Path, default=ROOT / "build/all_release/bin/runIPCSim")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--cpu-set", help="Optional Linux taskset CPU list; also constrains oneTBB CPU preparation.")
    parser.add_argument("--skip-state-check", action="store_true", help="Allow a quick incomplete summary; it cannot mark validation passed.")
    args = parser.parse_args()
    args.cases = list(VOLUME_IPC_CASES) if args.volume_ipc else list(SHELL_IPC_CASES) if args.shell_ipc else args.cases or list(PRIMARY_CASES)
    if args.summarize:
        if args.action not in (None, "summarize"):
            parser.error("--summarize conflicts with this action")
        args.action = "summarize"
    if args.action is None:
        parser.error("Specify prepare, run, or summarize / --summarize")
    args.destination = args.destination.resolve()
    if args.threads < 1:
        parser.error("--threads must be positive")
    if args.action == "prepare":
        prepare(args.destination, args.cases)
        print(args.destination / "manifest.json")
        return 0
    if args.action == "run":
        for name in args.cases:
            code = run_case(args.destination, name, args.binary.resolve(), args.threads, args.cpu_set)
            print(f"{name}: process exit {code}", flush=True)
            if code:
                return code
        return 0
    results = [summarize(args.destination, name, not args.skip_state_check) for name in args.cases]
    write_report(args.destination, results)
    for result in results:
        print(f"{result['case']}: frames={result['frames_recorded']}/{result['expected_frames']}, "
              f"failed_frames={len(result['failed_frames'])}, peak_pairs={result['peak_active_pairs']}, "
              f"passed={result['validation_passed']}")
    return 0 if all(result["validation_passed"] for result in results) else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, KeyError) as error:
        print(f"Validation error: {error}", file=sys.stderr)
        sys.exit(2)
