#!/usr/bin/env python3
"""Package final FBMS simulation assets from a generated pipeline."""

from __future__ import annotations

import argparse
import json
import shutil
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import run_fbms_all_asset_pipeline as pipeline


DEFAULT_CONFIG = pipeline.DEFAULT_CONFIG


@dataclass(frozen=True)
class PackageItem:
    kind: str
    case: str
    name: str
    source_veg: Path
    source_surface: Path
    package_veg: Path
    package_surface: Path
    summary: dict[str, Any]


def read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as fin:
        data = json.load(fin)
    if not isinstance(data, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return data


def summary_path_for_config(config: pipeline.PipelineConfig) -> Path:
    return pipeline.summary_path(config)


def package_output_dir(config: pipeline.PipelineConfig) -> Path:
    return config.working_directory / "simulation_package"


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Package final FBMS simulation assets.")
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG, help="Pipeline JSON config.")
    parser.add_argument("--summary", type=Path, default=None, help="Pipeline summary JSON.")
    parser.add_argument("--out", type=Path, default=None, help="Output package directory.")
    parser.add_argument("--case", action="append", default=[], help="Restrict to one FBMS case. Can be repeated.")
    parser.add_argument("--shape", action="append", default=[], help="Restrict to one TPMS shape. Can be repeated.")
    parser.add_argument("--dry-run", action="store_true", help="Plan package without copying files.")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing package files.")
    return parser.parse_args(argv)


def result_key(result: dict[str, Any]) -> tuple[str, str, str]:
    return (str(result.get("kind")), str(result.get("case")), str(result.get("name")))


def summary_lookup(summary: dict[str, Any]) -> dict[tuple[str, str, str], dict[str, Any]]:
    results = summary.get("results")
    if not isinstance(results, list):
        raise ValueError("summary.results must be an array")
    lookup: dict[tuple[str, str, str], dict[str, Any]] = {}
    for result in results:
        if not isinstance(result, dict):
            raise ValueError("summary.results[] must be objects")
        key = result_key(result)
        if key in lookup:
            raise ValueError(f"duplicate summary result: {key}")
        lookup[key] = result
    return lookup


def selected_cases(config: pipeline.PipelineConfig, requested: list[str]) -> tuple[pipeline.SourceCase, ...]:
    cases = pipeline.discover_fbms_cases(config.source_asset_directory)
    return pipeline.filter_cases(cases, requested)


def selected_shapes(config: pipeline.PipelineConfig, requested: list[str]) -> tuple[str, ...]:
    return pipeline.filter_shapes(config, requested)


def target_for_result(
    config: pipeline.PipelineConfig,
    kind: str,
    case_name: str,
    name: str,
    lookup: dict[tuple[str, str, str], dict[str, Any]],
) -> pipeline.PipelineTarget:
    if kind == "fbms":
        return pipeline.fbms_target(config, case_name)
    fbms_result = lookup.get(("fbms", case_name, case_name), {})
    budget = fbms_result.get("raw_volume")
    return pipeline.baseline_target(config, case_name, name, float(budget) if isinstance(budget, (int, float)) else 0.0)


def package_subdir(out: Path, kind: str, case_name: str, name: str) -> Path:
    if kind == "fbms":
        return out / case_name / "fbms"
    return out / case_name / "baseline" / name


def validate_result(result: dict[str, Any]) -> None:
    if result.get("boundary_passed") is not True or result.get("component_count") != 3:
        raise RuntimeError(
            f"asset is not packageable: {result.get('kind')}/{result.get('case')}/{result.get('name')}"
        )


def source_surface_for_result(target: pipeline.PipelineTarget, result: dict[str, Any]) -> Path:
    if result.get("boundary_source") == "repaired":
        return target.repaired_veg_obj
    return target.veg_obj


def build_package_items(
    config: pipeline.PipelineConfig,
    summary: dict[str, Any],
    out: Path,
    cases: tuple[pipeline.SourceCase, ...],
    shapes: tuple[str, ...],
) -> list[PackageItem]:
    lookup = summary_lookup(summary)
    items: list[PackageItem] = []
    for case in cases:
        result = lookup.get(("fbms", case.name, case.name))
        if result is not None:
            items.append(build_item(config, out, "fbms", case.name, case.name, result, lookup))
        for shape in shapes:
            result = lookup.get(("baseline", case.name, shape))
            if result is not None:
                items.append(build_item(config, out, "baseline", case.name, shape, result, lookup))
    return items


def build_item(
    config: pipeline.PipelineConfig,
    out: Path,
    kind: str,
    case_name: str,
    name: str,
    result: dict[str, Any],
    lookup: dict[tuple[str, str, str], dict[str, Any]],
) -> PackageItem:
    validate_result(result)
    target = target_for_result(config, kind, case_name, name, lookup)
    source_surface = source_surface_for_result(target, result)
    dest_dir = package_subdir(out, kind, case_name, name)
    return PackageItem(
        kind=kind,
        case=case_name,
        name=name,
        source_veg=target.veg,
        source_surface=source_surface,
        package_veg=dest_dir / f"{name}.veg",
        package_surface=dest_dir / f"{name}.veg.obj",
        summary=result,
    )


def ensure_sources_exist(items: list[PackageItem]) -> None:
    for item in items:
        for path in (item.source_veg, item.source_surface):
            if not path.exists():
                raise FileNotFoundError(path)


def copy_one(src: Path, dst: Path, overwrite: bool, dry_run: bool) -> None:
    if dst.exists() and not overwrite:
        raise FileExistsError(f"{dst} exists; pass --overwrite to replace it")
    print(f"+ copy {src} {dst}", flush=True)
    if dry_run:
        return
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        dst.unlink()
    shutil.copy2(src, dst)


def write_manifest(out: Path, manifest: dict[str, Any], overwrite: bool, dry_run: bool) -> None:
    manifest_path = out / "package_manifest.json"
    if manifest_path.exists() and not overwrite:
        raise FileExistsError(f"{manifest_path} exists; pass --overwrite to replace it")
    print(f"+ write {manifest_path}", flush=True)
    if dry_run:
        return
    out.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def manifest_for_items(
    config_path: Path,
    summary_path: Path,
    out: Path,
    items: list[PackageItem],
) -> dict[str, Any]:
    return {
        "generated_at_unix": int(time.time()),
        "config": str(config_path),
        "summary": str(summary_path),
        "output_directory": str(out),
        "assets": [
            {
                "kind": item.kind,
                "case": item.case,
                "name": item.name,
                "boundary_source": item.summary.get("boundary_source"),
                "boundary_repair_attempted": item.summary.get("boundary_repair_attempted"),
                "boundary_repair_volume_drift": item.summary.get("boundary_repair_volume_drift"),
                "source_veg": str(item.source_veg),
                "source_surface": str(item.source_surface),
                "package_veg": str(item.package_veg),
                "package_surface": str(item.package_surface),
            }
            for item in items
        ],
    }


def package_assets_from_args(argv: list[str] | None = None) -> dict[str, Any]:
    args = parse_args(argv)
    config = pipeline.load_config(args.config)
    summary_path = args.summary if args.summary is not None else summary_path_for_config(config)
    out = pipeline.repo_path(str(args.out)) if args.out is not None else package_output_dir(config)
    summary = read_json(summary_path)
    cases = selected_cases(config, args.case)
    shapes = selected_shapes(config, args.shape)
    items = build_package_items(config, summary, out, cases, shapes)
    ensure_sources_exist(items)
    manifest = manifest_for_items(args.config, summary_path, out, items)

    for item in items:
        copy_one(item.source_veg, item.package_veg, overwrite=args.overwrite, dry_run=args.dry_run)
        copy_one(item.source_surface, item.package_surface, overwrite=args.overwrite, dry_run=args.dry_run)
    write_manifest(out, manifest, overwrite=args.overwrite, dry_run=args.dry_run)
    print(f"packaged {len(items)} assets into {out}", flush=True)
    return manifest


def main() -> int:
    try:
        package_assets_from_args()
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
