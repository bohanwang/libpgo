from __future__ import annotations

import argparse
import csv
import json
import re
import statistics
from pathlib import Path


ACTIVE_SAMPLES_RE = re.compile(r"# external active samples: (\d+)")
SOLVER_TIME_RE = re.compile(r"Newton solve time: ([0-9.eE+-]+)")
STEP_RE = re.compile(r"T(\d+): \|\|g\|\|=([0-9.eE+-]+); Solver Ret: (-?\d+)")
ENERGY_MAIN_RE = re.compile(r"\s*main:\s*([0-9.eE+-]+)")
ENERGY_SUB_RE = re.compile(r"\s*sub (\d+):\s*([0-9.eE+-]+)")


def parse_log(log_path: Path, timestep: float) -> list[dict[str, float | int | None]]:
    rows: list[dict[str, float | int | None]] = []
    pending_active_samples = 0
    pending_solver_time = 0.0
    current: dict[str, float | int | None] | None = None

    for line in log_path.read_text().splitlines():
        match = ACTIVE_SAMPLES_RE.search(line)
        if match:
            pending_active_samples = int(match.group(1))
            continue

        match = SOLVER_TIME_RE.search(line)
        if match:
            pending_solver_time = float(match.group(1))
            continue

        match = STEP_RE.search(line)
        if match:
            if current is not None:
                rows.append(current)

            step = int(match.group(1))
            current = {
                "step": step,
                "time_sec": step * timestep,
                "active_samples": pending_active_samples,
                "solver_time_sec": pending_solver_time,
                "grad_norm": float(match.group(2)),
                "solver_ret": int(match.group(3)),
                "energy_main": None,
            }

            pending_active_samples = 0
            pending_solver_time = 0.0
            continue

        if current is None:
            continue

        match = ENERGY_MAIN_RE.search(line)
        if match:
            current["energy_main"] = float(match.group(1))
            continue

        match = ENERGY_SUB_RE.search(line)
        if match:
            current[f"energy_sub{int(match.group(1))}"] = float(match.group(2))

    if current is not None:
        rows.append(current)

    return rows


def compute_contact_intervals(rows: list[dict[str, float | int | None]]) -> list[dict[str, int]]:
    intervals: list[dict[str, int]] = []
    start_idx: int | None = None
    peak_active = -1
    peak_step = -1

    for idx, row in enumerate(rows):
        active = int(row["active_samples"])
        step = int(row["step"])

        if active > 0:
            if start_idx is None:
                start_idx = idx
                peak_active = active
                peak_step = step
            elif active > peak_active:
                peak_active = active
                peak_step = step
        elif start_idx is not None:
            start_step = int(rows[start_idx]["step"])
            end_step = int(rows[idx - 1]["step"])
            intervals.append(
                {
                    "start_step": start_step,
                    "end_step": end_step,
                    "duration_steps": end_step - start_step + 1,
                    "peak_active": peak_active,
                    "peak_step": peak_step,
                }
            )
            start_idx = None

    if start_idx is not None:
        start_step = int(rows[start_idx]["step"])
        end_step = int(rows[-1]["step"])
        intervals.append(
            {
                "start_step": start_step,
                "end_step": end_step,
                "duration_steps": end_step - start_step + 1,
                "peak_active": peak_active,
                "peak_step": peak_step,
            }
        )

    return intervals


def write_curve_csv(rows: list[dict[str, float | int | None]], curve_csv: Path) -> None:
    sub_indices = sorted(
        {
            int(key[len("energy_sub"):])
            for row in rows
            for key in row.keys()
            if key.startswith("energy_sub")
        }
    )
    fieldnames = [
        "step",
        "time_sec",
        "active_samples",
        "solver_time_sec",
        "grad_norm",
        "solver_ret",
        "energy_main",
    ] + [f"energy_sub{sub_idx}" for sub_idx in sub_indices]

    with curve_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            out_row = {key: row.get(key) for key in fieldnames}
            writer.writerow(
                {
                    key: ("" if value is None else value)
                    for key, value in out_row.items()
                }
            )


def build_summary(log_path: Path, curve_csv: Path, rows: list[dict[str, float | int | None]]) -> dict:
    active_rows = [row for row in rows if int(row["active_samples"]) > 0]
    solver_times = [float(row["solver_time_sec"]) for row in rows]
    first_nonzero = active_rows[0] if active_rows else None
    max_active_row = max(rows, key=lambda row: int(row["active_samples"]), default=None)

    return {
        "log_path": str(log_path),
        "curve_csv": str(curve_csv),
        "num_steps_recorded": len(rows),
        "first_step": int(rows[0]["step"]) if rows else None,
        "last_step": int(rows[-1]["step"]) if rows else None,
        "first_nonzero_step": int(first_nonzero["step"]) if first_nonzero else None,
        "first_nonzero_time_sec": float(first_nonzero["time_sec"]) if first_nonzero else None,
        "first_nonzero_active_samples": int(first_nonzero["active_samples"]) if first_nonzero else None,
        "last_nonzero_step": int(active_rows[-1]["step"]) if active_rows else None,
        "max_active_samples": int(max_active_row["active_samples"]) if max_active_row else None,
        "max_active_step": int(max_active_row["step"]) if max_active_row else None,
        "mean_active_nonzero": statistics.mean(int(row["active_samples"]) for row in active_rows) if active_rows else 0.0,
        "median_active_nonzero": statistics.median(int(row["active_samples"]) for row in active_rows) if active_rows else 0.0,
        "num_nonzero_steps": len(active_rows),
        "contact_intervals": compute_contact_intervals(rows),
        "num_solver_failures": sum(int(row["solver_ret"]) != 0 for row in rows),
        "max_solver_time_sec": max(solver_times) if solver_times else 0.0,
        "mean_solver_time_sec": statistics.mean(solver_times) if solver_times else 0.0,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract contact curves and summary statistics from a runSim log.")
    parser.add_argument("--log", type=Path, required=True)
    parser.add_argument("--timestep", type=float, required=True)
    parser.add_argument("--curve-csv", type=Path, required=True)
    parser.add_argument("--summary-json", type=Path, required=True)
    args = parser.parse_args()

    rows = parse_log(args.log, args.timestep)
    write_curve_csv(rows, args.curve_csv)
    summary = build_summary(args.log, args.curve_csv, rows)
    args.summary_json.write_text(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
