from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt


def slugify(label: str) -> str:
    slug = re.sub(r"[^a-zA-Z0-9]+", "_", label.strip().lower()).strip("_")
    return slug or "series"


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
                    "energy_main": float(row["energy_main"]) if row.get("energy_main") else 0.0,
                    "energy_sub0": float(row["energy_sub0"]) if row.get("energy_sub0") else 0.0,
                }
            )
    return rows


def load_summary(json_path: Path) -> dict:
    return json.loads(json_path.read_text())


def annotate_contact(ax, summary: dict, color: str) -> None:
    first_nonzero = summary.get("first_nonzero_step")
    max_active_step = summary.get("max_active_step")
    max_active = summary.get("max_active_samples")

    if first_nonzero is not None:
        ax.axvline(first_nonzero, color=color, linestyle="--", alpha=0.3, linewidth=1.0)

    if max_active_step is not None and max_active is not None:
        ax.scatter([max_active_step], [max_active], color=color, s=18, zorder=3)


def build_combined_rows(series: list[dict]) -> list[dict[str, float]]:
    step_maps: list[tuple[str, dict[float, dict[str, float]]]] = []
    all_steps: set[float] = set()
    for item in series:
        slug = item["slug"]
        row_map = {row["step"]: row for row in item["rows"]}
        step_maps.append((slug, row_map))
        all_steps.update(row_map.keys())

    combined_rows: list[dict[str, float]] = []
    for step in sorted(all_steps):
        row: dict[str, float] = {"step": step, "time_sec": step_maps[0][1][step]["time_sec"]}
        for slug, row_map in step_maps:
            source = row_map.get(step)
            row[f"{slug}_active_samples"] = source["active_samples"] if source else 0.0
            row[f"{slug}_solver_time_sec"] = source["solver_time_sec"] if source else 0.0
            row[f"{slug}_energy_main"] = source["energy_main"] if source else 0.0
            row[f"{slug}_energy_sub0"] = source["energy_sub0"] if source else 0.0
        combined_rows.append(row)

    return combined_rows


def write_combined_csv(output_path: Path, series: list[dict], rows: list[dict[str, float]]) -> None:
    fieldnames = ["step", "time_sec"]
    for item in series:
        slug = item["slug"]
        fieldnames.extend(
            [
                f"{slug}_active_samples",
                f"{slug}_solver_time_sec",
                f"{slug}_energy_main",
                f"{slug}_energy_sub0",
            ]
        )

    with output_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def build_comparison_summary(series: list[dict], combined_csv: Path) -> dict:
    models = {}
    for item in series:
        models[item["label"]] = item["summary"]

    ratios = {}
    for lhs in series:
        for rhs in series:
            if lhs["label"] == rhs["label"]:
                continue

            lhs_mean = lhs["summary"].get("mean_solver_time_sec", 0.0) or 0.0
            rhs_mean = rhs["summary"].get("mean_solver_time_sec", 0.0) or 0.0
            if rhs_mean > 0.0:
                ratios[f"{lhs['slug']}_vs_{rhs['slug']}_mean_solver_time_ratio"] = lhs_mean / rhs_mean

    return {
        "models": models,
        "ratios": ratios,
        "comparison_csv": str(combined_csv),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare multiple contact-model curves on one figure.")
    parser.add_argument(
        "--series",
        action="append",
        nargs=4,
        metavar=("LABEL", "COLOR", "CSV", "SUMMARY"),
        required=True,
        help="Add one series as: LABEL COLOR CURVE_CSV SUMMARY_JSON",
    )
    parser.add_argument("--output-prefix", type=Path, required=True)
    args = parser.parse_args()

    series = []
    for label, color, csv_path, summary_path in args.series:
        series.append(
            {
                "label": label,
                "slug": slugify(label),
                "color": color,
                "rows": load_curve(Path(csv_path)),
                "summary": load_summary(Path(summary_path)),
            }
        )

    fig, axes = plt.subplots(2, 2, figsize=(14, 8), dpi=180, constrained_layout=True)

    for item in series:
        steps = [row["step"] for row in item["rows"]]
        active = [row["active_samples"] for row in item["rows"]]
        solver_time = [row["solver_time_sec"] for row in item["rows"]]
        cumulative_solver_time = []
        cumulative = 0.0
        for value in solver_time:
            cumulative += value
            cumulative_solver_time.append(cumulative)

        axes[0, 0].plot(steps, active, label=item["label"], color=item["color"], linewidth=1.6)
        axes[0, 1].plot(steps, solver_time, label=item["label"], color=item["color"], linewidth=1.3)
        axes[1, 0].plot(steps, cumulative_solver_time, label=item["label"], color=item["color"], linewidth=1.6)
        annotate_contact(axes[0, 0], item["summary"], item["color"])

    axes[0, 0].set_title("Active Samples")
    axes[0, 0].set_xlabel("Step")
    axes[0, 0].set_ylabel("samples")
    axes[0, 0].grid(True, alpha=0.25)
    axes[0, 0].legend()

    axes[0, 1].set_title("Newton Solve Time per Step")
    axes[0, 1].set_xlabel("Step")
    axes[0, 1].set_ylabel("seconds")
    axes[0, 1].grid(True, alpha=0.25)
    axes[0, 1].legend()

    axes[1, 0].set_title("Cumulative Newton Solve Time")
    axes[1, 0].set_xlabel("Step")
    axes[1, 0].set_ylabel("seconds")
    axes[1, 0].grid(True, alpha=0.25)
    axes[1, 0].legend()

    axes[1, 1].axis("off")
    summary_lines = []
    for item in series:
        summary = item["summary"]
        summary_lines.extend(
            [
                f"{item['label']}:",
                f"  first contact: {summary.get('first_nonzero_step')}",
                f"  peak active: {summary.get('max_active_samples')} @ {summary.get('max_active_step')}",
                f"  mean solve: {summary.get('mean_solver_time_sec', 0.0):.4f}s",
                f"  max solve: {summary.get('max_solver_time_sec', 0.0):.4f}s",
                "",
            ]
        )
    axes[1, 1].text(0.0, 1.0, "\n".join(summary_lines).rstrip(), va="top", ha="left", family="monospace")

    fig.suptitle("cubic-box-ipc: IPC vs IPC+LS vs Penalty (core-release)", fontsize=16, fontweight="bold")

    png_path = args.output_prefix.with_suffix(".png")
    svg_path = args.output_prefix.with_suffix(".svg")
    fig.savefig(png_path)
    fig.savefig(svg_path)
    plt.close(fig)

    combined_rows = build_combined_rows(series)
    combined_csv = args.output_prefix.with_suffix(".csv")
    write_combined_csv(combined_csv, series, combined_rows)

    comparison_summary = build_comparison_summary(series, combined_csv)
    args.output_prefix.with_name(args.output_prefix.name + "-summary").with_suffix(".json").write_text(
        json.dumps(comparison_summary, indent=2)
    )


if __name__ == "__main__":
    main()
