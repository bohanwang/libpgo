from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt


def load_curve(csv_path: Path) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    with csv_path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(
                {
                    "step": float(row["step"]),
                    "time_sec": float(row["time_sec"]),
                    "active_samples": float(row["active_samples"]),
                    "solver_time_sec": float(row["solver_time_sec"]),
                    "grad_norm": float(row["grad_norm"]),
                    "solver_ret": float(row["solver_ret"]),
                    "energy_main": float(row["energy_main"]) if row["energy_main"] else 0.0,
                    "energy_sub0": float(row["energy_sub0"]) if row["energy_sub0"] else 0.0,
                }
            )
    return rows


def load_summary(json_path: Path) -> dict:
    return json.loads(json_path.read_text())


def plot_series(ax, rows: list[dict[str, float]], key: str, title: str, color: str, ylabel: str) -> None:
    ax.plot([r["step"] for r in rows], [r[key] for r in rows], color=color, linewidth=1.6)
    ax.set_title(title)
    ax.set_xlabel("Step")
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.25)


def annotate_contact(ax, summary: dict, color: str) -> None:
    first_nonzero = summary.get("first_nonzero_step")
    max_active_step = summary.get("max_active_step")
    max_active = summary.get("max_active_samples")

    if first_nonzero is not None:
        ax.axvline(first_nonzero, color=color, linestyle="--", alpha=0.55, linewidth=1.0)
        ax.text(first_nonzero, ax.get_ylim()[1] * 0.96, f"first: {first_nonzero}", color=color, rotation=90,
                va="top", ha="right", fontsize=9)

    if max_active_step is not None and max_active is not None:
        ax.scatter([max_active_step], [max_active], color=color, s=22, zorder=3)
        ax.text(max_active_step, max_active, f" peak {int(max_active)} @ {max_active_step}", color=color,
                va="bottom", ha="left", fontsize=9)


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare IPC barrier and penalty contact curves.")
    parser.add_argument("--barrier-csv", type=Path, required=True)
    parser.add_argument("--barrier-summary", type=Path, required=True)
    parser.add_argument("--penalty-csv", type=Path, required=True)
    parser.add_argument("--penalty-summary", type=Path, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    args = parser.parse_args()

    barrier_rows = load_curve(args.barrier_csv)
    barrier_summary = load_summary(args.barrier_summary)
    penalty_rows = load_curve(args.penalty_csv)
    penalty_summary = load_summary(args.penalty_summary)

    fig, axes = plt.subplots(2, 2, figsize=(14, 8), dpi=180, constrained_layout=True)
    plot_series(axes[0, 0], barrier_rows, "active_samples", "IPC Barrier: Active Samples", "#1f77b4", "samples")
    plot_series(axes[0, 1], penalty_rows, "active_samples", "Penalty: Active Samples", "#ff7f0e", "samples")
    plot_series(axes[1, 0], barrier_rows, "solver_time_sec", "IPC Barrier: Newton Solve Time", "#1f77b4", "seconds")
    plot_series(axes[1, 1], penalty_rows, "solver_time_sec", "Penalty: Newton Solve Time", "#ff7f0e", "seconds")

    annotate_contact(axes[0, 0], barrier_summary, "#1f77b4")
    annotate_contact(axes[0, 1], penalty_summary, "#ff7f0e")

    fig.suptitle("cubic-box-ipc: IPC Barrier vs Penalty (core-release)", fontsize=16, fontweight="bold")

    png_path = args.output_prefix.with_suffix(".png")
    svg_path = args.output_prefix.with_suffix(".svg")
    fig.savefig(png_path)
    fig.savefig(svg_path)
    plt.close(fig)


if __name__ == "__main__":
    main()
