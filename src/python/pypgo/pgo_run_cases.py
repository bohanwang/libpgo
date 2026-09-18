import argparse
import json
import re
import shlex
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pypgo


_CASE_NAME_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
_MANIFEST_KEYS = {"version", "abc-output-root", "cases", "groups"}
_CASE_KEYS = {
    "name",
    "simulation-config",
    "animation-config",
    "abc-output",
}


class ManifestError(ValueError):
    pass


@dataclass(frozen=True)
class Case:
    name: str
    simulation_config: Path
    animation_config: Path
    abc_output: Path


@dataclass(frozen=True)
class Manifest:
    cases: tuple[Case, ...]
    groups: dict[str, tuple[str, ...]]


def _expect_object(value: Any, context: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ManifestError(f"{context} must be a JSON object")
    return value


def _expect_string(value: Any, context: str) -> str:
    if not isinstance(value, str) or not value:
        raise ManifestError(f"{context} must be a non-empty string")
    return value


def _reject_unknown_keys(
    value: dict[str, Any], allowed: set[str], context: str
) -> None:
    unknown = sorted(set(value) - allowed)
    if unknown:
        raise ManifestError(f"{context} has unknown field(s): {', '.join(unknown)}")


def _resolve_path(manifest_directory: Path, value: Any, context: str) -> Path:
    configured_path = Path(_expect_string(value, context)).expanduser()
    if not configured_path.is_absolute():
        configured_path = manifest_directory / configured_path
    return configured_path.resolve()


def load_manifest(filename: str | Path) -> Manifest:
    manifest_path = Path(filename).expanduser().resolve()
    try:
        manifest_value = json.loads(manifest_path.read_text(encoding="utf-8"))
    except OSError as error:
        raise ManifestError(f"cannot read manifest {manifest_path}: {error}") from error
    except json.JSONDecodeError as error:
        raise ManifestError(f"invalid JSON in {manifest_path}: {error}") from error

    manifest = _expect_object(manifest_value, "manifest")
    _reject_unknown_keys(manifest, _MANIFEST_KEYS, "manifest")
    if manifest.get("version") != 1:
        raise ManifestError("manifest.version must be 1")

    raw_cases = manifest.get("cases")
    if not isinstance(raw_cases, list) or not raw_cases:
        raise ManifestError("manifest.cases must be a non-empty array")

    manifest_directory = manifest_path.parent
    abc_output_root = _resolve_path(
        manifest_directory,
        manifest.get("abc-output-root", "generated/abc"),
        "manifest.abc-output-root",
    )
    cases = []
    names = set()
    for index, raw_case in enumerate(raw_cases):
        context = f"manifest.cases[{index}]"
        case_value = _expect_object(raw_case, context)
        _reject_unknown_keys(case_value, _CASE_KEYS, context)

        name = _expect_string(case_value.get("name"), f"{context}.name")
        if not _CASE_NAME_PATTERN.fullmatch(name):
            raise ManifestError(
                f"{context}.name must contain only letters, digits, '.', '_', or '-'"
            )
        if name in names:
            raise ManifestError(f"duplicate case name: {name}")
        names.add(name)

        abc_output_value = case_value.get("abc-output")
        abc_output = (
            _resolve_path(
                manifest_directory,
                abc_output_value,
                f"{context}.abc-output",
            )
            if abc_output_value is not None
            else abc_output_root / name
        )
        cases.append(
            Case(
                name=name,
                simulation_config=_resolve_path(
                    manifest_directory,
                    case_value.get("simulation-config"),
                    f"{context}.simulation-config",
                ),
                animation_config=_resolve_path(
                    manifest_directory,
                    case_value.get("animation-config"),
                    f"{context}.animation-config",
                ),
                abc_output=abc_output,
            )
        )
    raw_groups = manifest.get("groups", {})
    if not isinstance(raw_groups, dict):
        raise ManifestError("manifest.groups must be a JSON object")

    groups = {}
    for name, raw_members in raw_groups.items():
        context = f"manifest.groups.{name}"
        if not isinstance(name, str) or not _CASE_NAME_PATTERN.fullmatch(name):
            raise ManifestError(
                "group names must contain only letters, digits, '.', '_', or '-'"
            )
        if not isinstance(raw_members, list) or not raw_members:
            raise ManifestError(f"{context} must be a non-empty array")
        if not all(isinstance(member, str) for member in raw_members):
            raise ManifestError(f"{context} must contain only case names")
        if len(raw_members) != len(set(raw_members)):
            raise ManifestError(f"{context} contains duplicate case names")
        unknown_members = sorted(set(raw_members) - names)
        if unknown_members:
            raise ManifestError(
                f"{context} references unknown case(s): {', '.join(unknown_members)}"
            )
        groups[name] = tuple(raw_members)

    return Manifest(cases=tuple(cases), groups=groups)


def _select_cases(
    manifest: Manifest,
    requested_names: Sequence[str],
    requested_groups: Sequence[str],
) -> list[Case]:
    cases = manifest.cases
    if not requested_names and not requested_groups:
        return list(cases)

    cases_by_name = {case.name: case for case in cases}
    unknown = sorted(set(requested_names) - set(cases_by_name))
    if unknown:
        raise ManifestError(f"unknown case(s): {', '.join(unknown)}")

    unknown_groups = sorted(set(requested_groups) - set(manifest.groups))
    if unknown_groups:
        raise ManifestError(f"unknown group(s): {', '.join(unknown_groups)}")

    expanded_names = list(requested_names)
    for group in requested_groups:
        expanded_names.extend(manifest.groups[group])

    selected = []
    seen = set()
    for name in expanded_names:
        if name not in seen:
            selected.append(cases_by_name[name])
            seen.add(name)
    return selected


def _validate_inputs(cases: Sequence[Case], stage: str) -> None:
    missing = []
    for case in cases:
        if stage in {"all", "simulation"} and not case.simulation_config.is_file():
            missing.append(f"{case.name}: {case.simulation_config}")
        if stage in {"all", "animation"} and not case.animation_config.is_file():
            missing.append(f"{case.name}: {case.animation_config}")
    if missing:
        raise ManifestError("missing config file(s):\n  " + "\n  ".join(missing))


def _display_commands(case: Case, stage: str) -> None:
    if stage in {"all", "simulation"}:
        print(shlex.join(["pgo-run-sim", str(case.simulation_config)]))
    if stage in {"all", "animation"}:
        print(f"mkdir -p {shlex.quote(str(case.abc_output))}")
        print(
            shlex.join(
                [
                    "pgo-dump-abc",
                    str(case.animation_config),
                    str(case.abc_output),
                ]
            )
        )


def _run_case(case: Case, stage: str) -> tuple[bool, str | None]:
    try:
        if stage in {"all", "simulation"}:
            print(
                f"[{case.name}] running simulation: {case.simulation_config}",
                flush=True,
            )
            status = pypgo.run_sim_from_config(str(case.simulation_config))
            if status != 0:
                return False, f"simulation returned {status}"

        if stage in {"all", "animation"}:
            case.abc_output.mkdir(parents=True, exist_ok=True)
            print(
                f"[{case.name}] exporting animation: {case.abc_output}",
                flush=True,
            )
            status = pypgo.convert_animation_to_abc(
                str(case.animation_config), str(case.abc_output)
            )
            if status != 0:
                return False, f"animation export returned {status}"
    except Exception as error:  # Native bindings surface failures as exceptions.
        return False, str(error)
    return True, None


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="pgo-run-cases",
        description="Run simulation and animation cases from a batch manifest.",
    )
    parser.add_argument("manifest", help="batch manifest JSON file")
    parser.add_argument(
        "cases",
        nargs="*",
        metavar="CASE",
        help="case names to run; omit to run every case",
    )
    parser.add_argument(
        "--stage",
        choices=("all", "simulation", "animation"),
        default="all",
        help="run simulations, export animations, or do both (default: all)",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="list selected cases without running them",
    )
    parser.add_argument(
        "--list-groups",
        action="store_true",
        help="list groups and their member cases without running them",
    )
    parser.add_argument(
        "--group",
        action="append",
        default=[],
        metavar="GROUP",
        help="run a named group; may be provided more than once",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="print the equivalent pypgo CLI commands without running them",
    )
    parser.add_argument(
        "--keep-going",
        action="store_true",
        help="continue with later cases after a failure",
    )
    arguments = parser.parse_args(argv)

    try:
        manifest = load_manifest(arguments.manifest)
        cases = _select_cases(manifest, arguments.cases, arguments.group)
        if not arguments.list and not arguments.list_groups:
            _validate_inputs(cases, arguments.stage)
    except ManifestError as error:
        parser.error(str(error))

    if arguments.list_groups:
        for name, members in manifest.groups.items():
            print(f"{name}: {', '.join(members)}")
        return 0

    if arguments.list:
        for case in cases:
            print(case.name)
        return 0

    if arguments.dry_run:
        for case in cases:
            _display_commands(case, arguments.stage)
        return 0

    failures = []
    succeeded = 0
    for case in cases:
        success, error = _run_case(case, arguments.stage)
        if success:
            succeeded += 1
            print(f"[{case.name}] completed")
            continue

        failures.append((case.name, error or "unknown failure"))
        print(f"[{case.name}] failed: {failures[-1][1]}")
        if not arguments.keep_going:
            break

    print(f"Batch summary: {succeeded} succeeded, {len(failures)} failed")
    if failures:
        for name, error in failures:
            print(f"  {name}: {error}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
