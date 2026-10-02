#!/usr/bin/env python3
"""Measure matched CPU/CUDA IPC checkpoint solves without changing their schedule.

Full-trajectory correctness evidence is linked separately. These measurements
are isolated checkpoint-step costs, not extrapolated full-simulation runtimes.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import math
import os
from pathlib import Path
import re
import shutil
import statistics
import subprocess
import time

import validate_cuda_examples as validation


ROOT = validation.ROOT
DEFAULT = ROOT / "examples/generated/cuda-performance-20261002"
REFERENCE = ROOT / "examples/generated/cuda-ipc-validation-20261002"
VARIANTS = ROOT / "examples/generated/cuda-material-variants-20261002"
NATIVE_RESULT = re.compile(
    r"Newton result: stop=([^,]+), iteration=(\d+), gradient_inf=([^,]+), "
    r"tolerance=([^,]+), converged=(true|false), status=([^\s]+)"
)


def digest(path: Path) -> str:
    return validation.digest(path)


def save(path: Path, data) -> None:
    validation.write_json(path, data)


def physical_config(config: dict) -> dict:
    excluded = validation.NONPHYSICAL_KEYS | {"restart-from-u"}
    return {key: value for key, value in config.items() if key not in excluded}


def object_digest(value) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def input_assets(value) -> set[Path]:
    result = set()
    if isinstance(value, list):
        for item in value:
            result.update(input_assets(item))
    elif isinstance(value, dict):
        for key, item in value.items():
            if key in validation.PATH_KEYS and isinstance(item, str):
                result.add(Path(item))
            else:
                result.update(input_assets(item))
    return result


def prepare_trial(directory: Path, config: dict, checkpoint: Path | None) -> Path:
    output = Path(config["output"])
    allowed = {"config.json", "output"}
    if directory.exists() and {p.name for p in directory.iterdir()} - allowed:
        raise ValueError(f"Trial directory contains existing run products: {directory}")
    config_path = directory / "config.json"
    if config_path.exists() and json.loads(config_path.read_text()) != config:
        raise ValueError(f"Trial configuration differs: {config_path}")
    output.mkdir(parents=True, exist_ok=True)
    expected = {checkpoint.name} if checkpoint else set()
    if {p.name for p in output.iterdir()} - expected:
        raise ValueError(f"Trial output contains more than its input checkpoint: {output}")
    if checkpoint:
        target = output / checkpoint.name
        if target.exists():
            if digest(target) != digest(checkpoint):
                raise ValueError(f"Trial input checkpoint differs: {target}")
        else:
            shutil.copyfile(checkpoint, target)
    if not config_path.exists():
        save(config_path, config)
    return config_path


def prepare(destination: Path, repeats: int) -> None:
    plan_path = destination / "benchmark-plan.json"
    if plan_path.exists():
        raise FileExistsError(f"Benchmark plan already exists: {plan_path}")
    cases = []
    for name in (*validation.VOLUME_IPC_CASES, *validation.SHELL_IPC_CASES):
        source, original = validation.original_config(name)
        cases.append((name, source, original, (REFERENCE / name).resolve(), "canonical IPC"))
    for name in ("box-tet", "box-cubic"):
        reference = VARIANTS / name
        source, original = validation.original_config(name)
        original["elastic-material"] = "stvk-vol"
        cases.append((name + "-stvk-vol", source, original, reference, "explicit StVK material variant"))

    trials = []
    entries = []
    pinned_inputs = {}
    for name, source, original, reference, kind in cases:
        summary_path = reference / "summary.json"
        if not summary_path.exists():
            summary_path = reference / "material-variant-summary.json"
        summary = json.loads(summary_path.read_text())
        if summary.get("validation_passed", summary.get("material_variant_passed")) is not True:
            raise ValueError(f"Incomplete reference trajectory: {reference}")
        previous_config = json.loads((reference / "config.json").read_text())
        for path in input_assets(original) | {source}:
            pinned_inputs[str(path)] = digest(path)
        with (reference / "output/newton-solves.csv").open() as stream:
            rows = list(csv.DictReader(stream))
        static = original["sim-type"] == "static"
        if static:
            phases = {0: ["static equilibrium"]}
        else:
            peak = max(rows, key=lambda row: int(row["active_pairs"]))
            phases = {}
            for phase, frame in (("initial", 0), ("peak contacts", int(peak["frame"])), ("late", int(rows[-1]["frame"]))):
                phases.setdefault(frame, []).append(phase)
        strategy = {key: previous_config[key] for key in validation.CUDA_STRATEGY_KEYS if key in previous_config}
        entry = dict(case=name, kind=kind, source=str(source), source_sha256=digest(source),
                     reference_directory=str(reference), reference_report_sha256=digest(reference / "collision-validation.json"),
                     original_frames=validation.expected_frames(original), original_tolerance=original["solver-eps"],
                     original_iteration_cap=original["solver-max-iter"], material=original["elastic-material"],
                     model="tet" if "tet-mesh" in original else "cubic" if "cubic-mesh" in original else "koiter",
                     numerical_strategy=strategy, full_validation_passed=True, phases=phases,
                     physical_config=physical_config(original), physical_config_sha256=object_digest(physical_config(original)))
        entries.append(entry)
        for frame, labels in phases.items():
            checkpoint = reference / "output" / f"deform{frame - 1:04d}.u" if frame > 0 else None
            if checkpoint and not checkpoint.is_file():
                raise FileNotFoundError(checkpoint)
            for repeat in range(1 if static else repeats):
                # Alternate order to avoid always measuring one backend first.
                for backend in (("cpu", "cuda") if repeat % 2 == 0 else ("cuda", "cpu")):
                    trial_id = f"{name}/frame{frame:04d}/{backend}-r{repeat + 1}"
                    directory = destination / "solver" / trial_id
                    output = directory / "output"
                    config = dict(original, output=str(output), **{"simulation-backend": backend})
                    if frame > 0:
                        config["restart-from-u"] = True
                    if backend == "cuda":
                        config.update(strategy)
                        config.update({"cuda-device": 0, "cuda-solver-verbose": 2})
                        if "cubic-mesh" in config:
                            config["cuda-cubic-validity"] = "gauss"
                    config_path = prepare_trial(directory, config, checkpoint)
                    assert physical_config(config) == physical_config(original)
                    trials.append(dict(id=trial_id, case=name, frame=frame, phases=labels, backend=backend,
                                       repeat=repeat + 1, config=str(config_path), config_sha256=digest(config_path),
                                       checkpoint=str(checkpoint) if checkpoint else None,
                                       checkpoint_sha256=digest(checkpoint) if checkpoint else None,
                                       expected_frames=1, max_frames=1, tolerance=original["solver-eps"],
                                       iteration_cap=original["solver-max-iter"], static=static,
                                       physical_config_sha256=object_digest(physical_config(original)),
                                       expected_cpu_solves=original.get("ipc-friction-iterations", 1) if original.get("contact-friction-coeff", 0) > 0 else 1))

    # This checked-in CUDA demonstration is a separate, shorter full run.
    source = ROOT / "examples/configs/volume/box/box-tet-ipc-cuda.json"
    original = validation.absolute_paths(json.loads(source.read_text()), source.parent)
    for path in input_assets(original) | {source}:
        pinned_inputs[str(path)] = digest(path)
    for backend in ("cpu", "cuda"):
        directory = destination / "solver" / "box-tet-cuda-750" / backend
        output = directory / "output"
        config = {key: value for key, value in original.items() if not key.startswith("cuda-")}
        config.update({"simulation-backend": backend, "output": str(output)})
        if backend == "cuda":
            config.update({"cuda-device": 0, "cuda-solver-verbose": 2})
        prepare_trial(directory, config, None)
        assert physical_config(config) == physical_config(original)
        trials.append(dict(id=f"box-tet-cuda-750/{backend}", case="box-tet-cuda-750", frame=0,
                           phases=["full 750-frame CUDA demonstration"], backend=backend, repeat=1,
                           config=str(directory / "config.json"), config_sha256=digest(directory / "config.json"),
                           checkpoint=None, checkpoint_sha256=None, expected_frames=750, max_frames=-1,
                           tolerance=config["solver-eps"], iteration_cap=config["solver-max-iter"], static=False,
                           physical_config_sha256=object_digest(physical_config(original)), expected_cpu_solves=750))
    plan = dict(schema_version=1, cases=entries, trials=trials, dynamic_repeats=repeats, static_repeats=1,
                pinned_inputs_sha256=pinned_inputs,
                methodology={
                    "comparison": "Original inputs and identical u/v/a checkpoint for each paired solve; original num-timestep retained.",
                    "sample_scope": "Initial, maximum-reference-contact and final macroframe, not a full runtime extrapolation.",
                    "warmup_scope": "Fresh-process checkpoint latency; no in-process warmup. Separate elastic operator timings use warmup and repeated resident-vector batches.",
                    "static_scope": "One complete original equilibrium solve per backend; no repeat-based variability estimate.",
                    "cpu_resources": "8 physical CPU cores; one process at a time, including all CPU and CUDA trials.",
                    "gpu_resources": "One H200 NVL; no simultaneous simulation, render, or equivalence workload.",
                    "timings": "Frame compute and export timers plus process wall time; process includes setup and final checks.",
                    "validation": "Original Newton tolerance and cap; independent CPU checks run after timing, including recorded physical motion.",
                    "baseline_limit": "Native CPU may accept nonconvergence and lacks CUDA adaptive physical-chord enforcement; such trials are explicitly marked, not counted as valid speedups.",
                    "stvk_scope": "Two additional explicit material variants, separate from the unchanged canonical suite.",
                })
    save(plan_path, plan)
    print(f"Prepared {len(trials)} isolated trials in {plan_path}")


def native_results(log: str) -> list[dict]:
    return [dict(stop=m[1], iterations=int(m[2]), gradient_inf=float(m[3]), tolerance=float(m[4]),
                 converged=m[5] == "true", status=m[6]) for m in NATIVE_RESULT.finditer(log)]


def finite_states(directory: Path, first: int, count: int) -> bool:
    import array
    import struct
    for frame in range(first, first + count):
        path = directory / f"deform{frame:04d}.u"
        if not path.is_file():
            return False
        data = path.read_bytes()
        if len(data) < 12:
            return False
        rows, columns, width = struct.unpack("=iii", data[:12])
        if rows <= 0 or columns != 3 or width != 8 or len(data) != 12 + rows * columns * width:
            return False
        values = array.array("d")
        values.frombytes(data[12:])
        if not all(map(math.isfinite, values)):
            return False
    return True


def valid_frame_timings(rows: list[dict], trial: dict) -> bool:
    try:
        if [int(row["frame"]) for row in rows] != list(range(trial["frame"], trial["frame"] + trial["expected_frames"])):
            return False
        for row in rows:
            compute, output, total = (float(row[key]) for key in ("compute_seconds", "output_seconds", "total_seconds"))
            if not all(math.isfinite(value) and value >= 0 for value in (compute, output, total)):
                return False
            if not math.isclose(compute + output, total, rel_tol=1e-12, abs_tol=1e-12):
                return False
            if row["committed"] != "1" or row["accepted_solvers_converged"] != "1" or int(row["terminal_status"]) != 0 or int(row["newton_solves"]) < 1:
                return False
            if trial["backend"] == "cpu" and row["all_solver_status_converged"] != "1":
                return False
        if trial["backend"] == "cpu" and sum(int(row["newton_solves"]) for row in rows) != trial["expected_cpu_solves"]:
            return False
        return True
    except (KeyError, TypeError, ValueError):
        return False


def run(destination: Path, binary: Path, checker: Path, cpu_set: str, timeout: float,
        selected: set[str] | None, limit: int | None) -> None:
    plan = json.loads((destination / "benchmark-plan.json").read_text())
    for filename, expected in plan["pinned_inputs_sha256"].items():
        if digest(Path(filename)) != expected:
            raise ValueError(f"Pinned input changed: {filename}")
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS="8", MKL_NUM_THREADS="8", OPENBLAS_NUM_THREADS="8", CUDA_VISIBLE_DEVICES="0")
    hashes = {"runner_sha256": digest(binary), "checker_sha256": digest(checker), "harness_sha256": digest(Path(__file__))}
    completed = 0
    for trial in plan["trials"]:
        if selected and trial["case"] not in selected:
            continue
        directory = Path(trial["config"]).parent
        record_path = directory / "run-record.json"
        if record_path.exists():
            previous = json.loads(record_path.read_text())
            if previous.get("complete") is True and all(previous.get(k) == v for k, v in hashes.items()):
                continue
            raise RuntimeError(f"Existing incomplete or different-binary trial: {directory}")
        if limit is not None and completed >= limit:
            break
        assert digest(Path(trial["config"])) == trial["config_sha256"]
        config = json.loads(Path(trial["config"]).read_text())
        assert object_digest(physical_config(config)) == trial["physical_config_sha256"]
        output = Path(config["output"])
        expected_seed = {Path(trial["checkpoint"]).name} if trial["checkpoint"] else set()
        if {p.name for p in output.iterdir()} != expected_seed:
            raise ValueError(f"Unexpected stale output before trial: {output}")
        if {p.name for p in directory.iterdir()} - {"output", "config.json"}:
            raise ValueError(f"Unexpected stale artifacts before trial: {directory}")
        if trial["checkpoint"]:
            seed = directory / "output" / Path(trial["checkpoint"]).name
            assert digest(seed) == digest(Path(trial["checkpoint"])) == trial["checkpoint_sha256"]
        activity = subprocess.run(["nvidia-smi", "--query-compute-apps=pid,process_name", "--format=csv,noheader"],
                                  capture_output=True, text=True, check=True).stdout.strip()
        if activity:
            raise RuntimeError("Isolated measurement requires idle GPUs; observed: " + activity)
        timing_path = directory / "frame-timing.csv"
        command = ["taskset", "-c", cpu_set, str(binary), trial["config"], "--timing-csv", str(timing_path)]
        if trial["max_frames"] > 0:
            command += ["--max-frames", str(trial["max_frames"])]
        record = dict(trial, **hashes, command=command, cpu_affinity=cpu_set,
                      thread_environment={key: env[key] for key in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "CUDA_VISIBLE_DEVICES")},
                      started_unix=time.time(), complete=False, timed_out=False)
        save(record_path, record)
        print("START", trial["id"], flush=True)
        start = time.monotonic()
        with (directory / "run.log").open("w") as log:
            try:
                record["exit_code"] = subprocess.run(command, cwd=ROOT, env=env, stdout=log,
                                                      stderr=subprocess.STDOUT, timeout=timeout).returncode
            except subprocess.TimeoutExpired:
                record.update(exit_code=None, timed_out=True)
        record.update(wall_seconds=time.monotonic() - start, finished_unix=time.time())
        if timing_path.exists():
            with timing_path.open() as stream:
                record["frame_timings"] = list(csv.DictReader(stream))
        else:
            record["frame_timings"] = []
        record["valid_frame_timings"] = valid_frame_timings(record["frame_timings"], trial)
        record["frame_timings_sha256"] = digest(timing_path) if timing_path.exists() else None
        if trial["backend"] == "cuda":
            rows, incomplete = validation.read_solver_rows(output / "newton-solves.csv")
            record["newton_results"] = rows
            requested_frames = list(range(trial["frame"], trial["frame"] + trial["expected_frames"]))
            strict = [row["frame"] for row in rows] == requested_frames and not incomplete and all(
                row["status"] == 0 and math.isfinite(row["gradient_inf"]) and row["gradient_inf"] < trial["tolerance"]
                and row["tolerance"] == trial["tolerance"]
                and row.get("max_substep_newton_iterations", row["newton_iterations"]) <= trial["iteration_cap"] for row in rows)
        else:
            rows = native_results((directory / "run.log").read_text(errors="replace"))
            record["newton_results"] = rows
            strict = len(rows) == trial["expected_cpu_solves"] and all(
                row["converged"] and row["status"] == "converged" and math.isfinite(row["gradient_inf"])
                and row["tolerance"] == trial["tolerance"]
                and row["gradient_inf"] < trial["tolerance"] and row["iterations"] <= trial["iteration_cap"] for row in rows)
        record["finite_complete_states"] = finite_states(output, trial["frame"], trial["expected_frames"])
        record["strict_convergence"] = bool(strict and record["exit_code"] == 0 and record["finite_complete_states"])
        save(record_path, record)
        if record["finite_complete_states"]:
            report = directory / "collision-validation.json"
            command = ["taskset", "-c", cpu_set, str(checker), trial["config"], "--first-frame", str(trial["frame"]),
                       "--last-frame", str(trial["frame"] + trial["expected_frames"] - 1), "--swept",
                       "--report", str(report), "--csv", str(directory / "collision-frames.csv")]
            if "cubic-mesh" in config:
                command += ["--cubic-validity", "gauss"]
            checked = dict(command=command, started_unix=time.time(), excluded_from_timing=True)
            with (directory / "checker.log").open("w") as log:
                try:
                    checked["exit_code"] = subprocess.run(command, cwd=ROOT, env=env, stdout=log,
                                                          stderr=subprocess.STDOUT, timeout=timeout).returncode
                except subprocess.TimeoutExpired:
                    checked.update(exit_code=None, timed_out=True)
            checked["finished_unix"] = time.time()
            if report.exists():
                proof = json.loads(report.read_text())
                checks = {
                    "checker_exit_zero": checked["exit_code"] == 0,
                    "report_passed": proof.get("passed") is True,
                    "frame_count": proof.get("frames_checked") == proof.get("frames_expected") == trial["expected_frames"],
                    "no_failed_frames": proof.get("frames_failed") == 0,
                    "first_frame": proof.get("first_frame") == trial["frame"],
                    "last_frame": proof.get("last_frame") == trial["frame"] + trial["expected_frames"] - 1,
                    "swept": proof.get("swept_checks_performed") is True,
                    "certified": proof.get("uncertified_sweep_segments") == 0,
                    "config": proof.get("config") == trial["config"],
                    "states": proof.get("state_directory") == str(output),
                    "criterion": "cubic-mesh" not in config or proof.get("element_validity_criterion") == "gauss",
                }
                checked.update(report_sha256=digest(report), passed=all(checks.values()), checks=checks)
            else:
                checked["passed"] = False
            save(directory / "checker-record.json", checked)
            record["geometry_passed"] = checked["passed"]
        else:
            record["geometry_passed"] = False
        record.update(complete=True, valid_comparison_result=record["strict_convergence"] and record["geometry_passed"] and record["valid_frame_timings"],
                      run_log_sha256=digest(directory / "run.log"))
        save(record_path, record)
        print("DONE", trial["id"], "wall", round(record["wall_seconds"], 6), "strict", record["strict_convergence"],
              "geometry", record["geometry_passed"], flush=True)
        completed += 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=["prepare", "run"])
    parser.add_argument("--destination", type=Path, default=DEFAULT)
    parser.add_argument("--binary", type=Path, default=ROOT / "build/all_release/bin/runIPCSim")
    parser.add_argument("--checker", type=Path, default=ROOT / "build/all_release/bin/validateIPCTrajectory")
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--cpu-set", default="8-15")
    parser.add_argument("--timeout", type=float, default=1800)
    parser.add_argument("--cases", nargs="+")
    parser.add_argument("--limit", type=int)
    args = parser.parse_args()
    if args.repeats < 1 or args.timeout <= 0:
        parser.error("repeats and timeout must be positive")
    if args.action == "prepare":
        prepare(args.destination.resolve(), args.repeats)
    else:
        run(args.destination.resolve(), args.binary.resolve(), args.checker.resolve(), args.cpu_set, args.timeout,
            set(args.cases) if args.cases else None, args.limit)


if __name__ == "__main__":
    main()
