#!/usr/bin/env python3
"""Short-run every retained simulation and animation example config."""

from __future__ import annotations

import argparse
import json
import math
import shutil
import subprocess
import tempfile
import time
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
CONFIG_ROOT = REPO_ROOT / "examples" / "configs"
GENERATOR = Path(__file__).with_name("generate_lite_tet_cubic_asset.py")


def find_binary(build_dir: Path, name: str) -> Path:
    suffixes = (f"{name}.exe", name)
    candidates = [
        build_dir / directory / suffix
        for directory in ("bin", "src/tools", "src/tools/runSim", "")
        for suffix in suffixes
    ]
    for candidate in candidates:
        resolved = candidate.expanduser().resolve()
        if resolved.is_file():
            return resolved
    raise FileNotFoundError(f"Could not find {name} below {build_dir}")


def resolve_input(config_path: Path, value: str) -> str:
    path = Path(value)
    if path.is_absolute():
        return str(path)
    return str((config_path.parent / path).resolve())


def scene_name(source: Path) -> str | None:
    parts = source.relative_to(CONFIG_ROOT).parts
    if parts[0] == "shell":
        return None
    directory = parts[1]
    if directory == "dragon-dyn":
        return "dragon"
    return "box-with-sphere" if directory == "box-with-sphere" else directory


def lite_constraint(
    generated_assets: Path, scene: str, domain: str, original: str
) -> str:
    name = Path(original).name
    kind = "zmin" if "zmin" in name else "zmax" if "zmax" in name else "hang"
    return str(generated_assets / f"{scene}-{domain}-{kind}.txt")


def bound_nonzero_movement(item: dict, maximum: float = 1.0e-4) -> None:
    movement = item.get("movement")
    if not movement:
        return
    magnitude = max(abs(float(value)) for value in movement)
    if magnitude > maximum:
        item["movement"] = [float(value) * maximum / magnitude for value in movement]


def materialize_simulation(
    source: Path, destination: Path, generated_assets: Path, output_root: Path
) -> tuple[dict, Path]:
    config = json.loads(source.read_text())

    scene = scene_name(source)
    domain = "cubic" if "cubic-mesh" in config else "tet" if "tet-mesh" in config else None
    if scene is not None:
        assert domain is not None
        config[f"{domain}-mesh"] = str(generated_assets / f"{scene}-{domain}.veg")
        config["surface-mesh"] = str(generated_assets / f"{scene}-surface.obj")
        for item in config.get("fixed-vertices", []):
            item["filename"] = lite_constraint(
                generated_assets, scene, domain, item["filename"]
            )
            bound_nonzero_movement(item)
    else:
        config["surface-mesh"] = resolve_input(source, config["surface-mesh"])
        for item in config.get("fixed-vertices", []):
            item["filename"] = resolve_input(source, item["filename"])
            bound_nonzero_movement(item)
    for item in config.get("external-objects", []):
        item["filename"] = resolve_input(source, item["filename"])
        bound_nonzero_movement(item)

    stem = source.stem
    config["solver-max-iter"] = min(config["solver-max-iter"], 5)
    if config["sim-type"] == "dynamic":
        config["num-timestep"] = 2
        config["dump-interval"] = 1
        output = output_root / stem
    else:
        output = output_root / f"{stem}.obj"
    config["output"] = str(output)

    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(config, indent=2) + "\n")
    return config, output


def materialize_animation(
    source: Path, destination: Path, generated_assets: Path, output_root: Path
) -> Path:
    config = json.loads(source.read_text())
    simulation_stem = source.stem.removesuffix("-animation")
    scene = scene_name(source)
    for mesh in config["meshes"]:
        mesh["driving-mesh"] = (
            str(generated_assets / f"{scene}-surface.obj")
            if scene is not None
            else resolve_input(source, mesh["driving-mesh"])
        )
        mesh["sequence"] = str(output_root / simulation_stem / "ret{:04d}.obj")
        mesh["sequence-range"] = [0, 1]

    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(config, indent=2) + "\n")
    return destination


def run_case(
    name: str,
    command: list[str],
    log_path: Path,
    timeout: int,
) -> dict:
    print(f"[smoke] {name}", flush=True)
    start = time.monotonic()
    try:
        with log_path.open("w") as log:
            result = subprocess.run(
                command,
                cwd=REPO_ROOT,
                stdout=log,
                stderr=subprocess.STDOUT,
                text=True,
                timeout=timeout,
                check=False,
            )
        return_code = result.returncode
    except subprocess.TimeoutExpired:
        return_code = -1
        with log_path.open("a") as log:
            log.write(f"\nTimed out after {timeout} seconds.\n")
    elapsed = time.monotonic() - start
    return {
        "name": name,
        "command": command,
        "return_code": return_code,
        "elapsed_seconds": elapsed,
        "log": str(log_path),
    }


def obj_is_finite(path: Path) -> bool:
    vertex_count = 0
    try:
        for line in path.read_text().splitlines():
            if not line.startswith("v "):
                continue
            fields = line.split()
            if len(fields) < 4 or not all(
                math.isfinite(float(value)) for value in fields[1:4]
            ):
                return False
            vertex_count += 1
    except (OSError, ValueError):
        return False
    return vertex_count > 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument(
        "--work-dir",
        type=Path,
        help="Fresh directory for generated meshes, configs, outputs, and logs",
    )
    parser.add_argument("--timeout", type=int, default=600, help="Seconds per case")
    parser.add_argument(
        "--cubic-resolution",
        type=int,
        default=2,
        help="Cubic resolution for generated smoke assets (default: 2)",
    )
    parser.add_argument(
        "--keep",
        action="store_true",
        help="Keep an automatically created work directory",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    build_dir = args.build_dir.expanduser().resolve()
    binaries = {
        name: find_binary(build_dir, name)
        for name in (
            "tetMesher",
            "cubicMesher",
            "runSim",
            "runIPCSim",
            "runShellSim",
            "convertAnimation",
        )
    }

    temporary = None
    if args.work_dir:
        work_dir = args.work_dir.expanduser().resolve()
        if work_dir.exists():
            raise FileExistsError(f"Refusing existing work directory: {work_dir}")
        work_dir.mkdir(parents=True)
    else:
        temporary = tempfile.TemporaryDirectory(prefix="libpgo-example-smoke-")
        work_dir = Path(temporary.name)

    generated_assets = work_dir / "generated" / "smoke-assets"
    output_root = work_dir / "output"
    config_root = work_dir / "configs"
    log_root = work_dir / "logs"
    log_root.mkdir(parents=True)
    output_root.mkdir(parents=True)

    generator_command = [
        str(GENERATOR),
        "--tet-mesher",
        str(binaries["tetMesher"]),
        "--cubic-mesher",
        str(binaries["cubicMesher"]),
        "--output-dir",
        str(generated_assets),
        "--cubic-resolution",
        str(args.cubic_resolution),
    ]
    subprocess.run(generator_command, cwd=REPO_ROOT, check=True)

    simulation_sources = []
    animation_sources = []
    for source in sorted(CONFIG_ROOT.rglob("*.json")):
        config = json.loads(source.read_text())
        (animation_sources if "meshes" in config else simulation_sources).append(
            source
        )

    results = []
    failed = False
    for source in simulation_sources:
        relative = source.relative_to(CONFIG_ROOT)
        destination = config_root / relative
        config, output = materialize_simulation(
            source, destination, generated_assets, output_root
        )
        contact = config["contact-model"]
        is_shell = "tet-mesh" not in config and "cubic-mesh" not in config
        runner = (
            "runIPCSim"
            if contact == "ipc"
            else "runShellSim"
            if is_shell
            else "runSim"
        )
        result = run_case(
            str(relative),
            [str(binaries[runner]), str(destination)],
            log_root / f"{source.stem}.log",
            args.timeout,
        )

        output_files = (
            [output]
            if config["sim-type"] == "static" and output.is_file()
            else sorted(output.glob("ret*.obj"))
            if config["sim-type"] == "dynamic"
            else []
        )
        output_ok = bool(output_files) and all(obj_is_finite(path) for path in output_files)
        result["output_files"] = [str(path) for path in output_files]
        result["finite_output"] = output_ok
        result["output_ok"] = output_ok
        result["passed"] = result["return_code"] == 0 and output_ok
        failed |= not result["passed"]
        results.append(result)

    for source in animation_sources:
        relative = source.relative_to(CONFIG_ROOT)
        destination = materialize_animation(
            source, config_root / relative, generated_assets, output_root
        )
        outputs_before = set(destination.parent.glob("*.abc"))
        result = run_case(
            str(relative),
            [str(binaries["convertAnimation"]), str(destination)],
            log_root / f"{source.stem}.log",
            args.timeout,
        )
        output_files = sorted(set(destination.parent.glob("*.abc")) - outputs_before)
        result["output_files"] = [str(path) for path in output_files]
        result["output_ok"] = bool(output_files)
        result["passed"] = result["return_code"] == 0 and result["output_ok"]
        failed |= not result["passed"]
        results.append(result)

    asset_manifest = json.loads((generated_assets / "manifest.json").read_text())
    report = {
        "work_dir": str(work_dir),
        "asset_parameters": asset_manifest["parameters"],
        "asset_sizes": {
            asset["scene"]: {
                domain: {
                    "vertices": asset["meshes"][domain]["vertex_count"],
                    "elements": asset["meshes"][domain]["element_count"],
                }
                for domain in ("tet", "cubic")
            }
            for asset in asset_manifest["assets"]
        },
        "simulation_count": len(simulation_sources),
        "animation_count": len(animation_sources),
        "passed": sum(result["passed"] for result in results),
        "failed": sum(not result["passed"] for result in results),
        "results": results,
    }
    report_path = work_dir / "report.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({key: report[key] for key in report if key != "results"}, indent=2))

    if temporary is not None and (args.keep or failed):
        retained = Path(tempfile.mkdtemp(prefix="libpgo-example-smoke-results-"))
        shutil.copytree(work_dir, retained, dirs_exist_ok=True)
        print(f"Retained smoke results at {retained}")
    if temporary is not None:
        temporary.cleanup()
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
