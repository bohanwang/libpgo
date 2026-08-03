#!/usr/bin/env python3
"""Guard and record provenance for wheel-only libpgo releases."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import zipfile
from datetime import datetime, timezone
from email.parser import BytesParser
from pathlib import Path

from project_version import read_project_version
from pypgo_wheel.contracts import SUPPORTED_PLATFORMS
from pypgo_wheel.contracts.common import RELEASE_DISTRIBUTION

SOURCE_RECORD = "source-provenance.json"
FINAL_RECORD = "wheel-provenance.json"
SDIST_SUFFIXES = (".tar.gz", ".tar.bz2", ".tar.xz", ".tgz", ".zip")


class ProvenanceError(RuntimeError):
    pass


def run_git(source_dir: Path, *args: str) -> str:
    result = subprocess.run(
        ["git", "-c", f"safe.directory={source_dir}", *args],
        cwd=source_dir,
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def ensure_outside_source(path: Path, source_dir: Path, label: str) -> None:
    try:
        path.relative_to(source_dir)
    except ValueError:
        return
    raise ProvenanceError(f"{label} must be outside the source checkout: {path}")


def require_clean_checkout(source_dir: Path) -> dict[str, str]:
    source_dir = source_dir.resolve()
    git_root = Path(run_git(source_dir, "rev-parse", "--show-toplevel")).resolve()
    if not source_dir.samefile(git_root):
        raise ProvenanceError(
            f"source directory is not the Git root: {source_dir} "
            f"(Git reported {git_root})"
        )

    status = run_git(
        source_dir,
        "status",
        "--porcelain=v1",
        "--untracked-files=all",
    )
    if status:
        raise ProvenanceError(
            "release packaging requires a completely clean checkout; "
            f"Git reported:\n{status}"
        )

    commit = run_git(source_dir, "rev-parse", "HEAD")
    github_sha = os.environ.get("GITHUB_SHA")
    if github_sha and github_sha != commit:
        raise ProvenanceError(
            f"GITHUB_SHA ({github_sha}) does not match checked-out HEAD ({commit})"
        )

    return {
        "commit": commit,
        "tree": run_git(source_dir, "rev-parse", "HEAD^{tree}"),
    }


def write_json(path: Path, payload: dict[str, object]) -> None:
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def command_preflight(args: argparse.Namespace) -> None:
    source_dir = args.source_dir.resolve()
    release_version = read_project_version(source_dir)
    evidence_dir = args.evidence_dir.resolve()
    ensure_outside_source(evidence_dir, source_dir, "evidence directory")
    checkout = require_clean_checkout(source_dir)

    evidence_dir.mkdir(parents=True, exist_ok=True)
    record = {
        "schema_version": 1,
        "release_version": release_version,
        "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
        "source": checkout,
    }
    output = evidence_dir / SOURCE_RECORD
    write_json(output, record)
    print(output)


def wheel_metadata(
    wheel: Path,
    release_version: str,
) -> dict[str, object]:
    with zipfile.ZipFile(wheel) as archive:
        names = archive.namelist()
        metadata_names = [
            name for name in names if name.endswith(".dist-info/METADATA")
        ]
        wheel_names = [name for name in names if name.endswith(".dist-info/WHEEL")]
        if len(metadata_names) != 1 or len(wheel_names) != 1:
            raise ProvenanceError(
                f"{wheel.name} must contain exactly one METADATA and one WHEEL file"
            )

        metadata = BytesParser().parsebytes(archive.read(metadata_names[0]))
        distribution = metadata.get("Name")
        if distribution != RELEASE_DISTRIBUTION:
            raise ProvenanceError(
                f"{wheel.name} contains distribution {distribution!r}, "
                f"expected {RELEASE_DISTRIBUTION!r}"
            )
        version = metadata.get("Version")
        if version != release_version:
            raise ProvenanceError(
                f"{wheel.name} contains version {version!r}, expected {release_version}"
            )

        wheel_text = archive.read(wheel_names[0]).decode("utf-8")
        tags = [
            line.removeprefix("Tag: ").strip()
            for line in wheel_text.splitlines()
            if line.startswith("Tag: ")
        ]
        return {
            "distribution": distribution,
            "version": version,
            "tags": tags,
            "members": len(names),
        }


def evidence_entry(path: Path) -> dict[str, object]:
    return {
        "name": path.name,
        "sha256": sha256(path),
        "size": path.stat().st_size,
    }


def command_record(args: argparse.Namespace) -> None:
    source_dir = args.source_dir.resolve()
    release_version = read_project_version(source_dir)
    evidence_dir = args.evidence_dir.resolve()
    wheel_dir = args.wheel_dir.resolve()
    ensure_outside_source(evidence_dir, source_dir, "evidence directory")
    ensure_outside_source(wheel_dir, source_dir, "wheel directory")
    source_record_path = evidence_dir / SOURCE_RECORD
    if not source_record_path.is_file():
        raise ProvenanceError(
            f"missing {SOURCE_RECORD}; run the preflight command before packaging"
        )

    source_record = json.loads(source_record_path.read_text(encoding="utf-8"))
    checkout = require_clean_checkout(source_dir)
    if source_record.get("source") != checkout:
        raise ProvenanceError("checkout changed after the release preflight")

    wheels = sorted(wheel_dir.glob("*.whl"))
    if len(wheels) != 1:
        raise ProvenanceError(
            f"expected exactly one wheel in {wheel_dir}, found {len(wheels)}"
        )

    sdists = [
        path
        for path in wheel_dir.iterdir()
        if path.is_file() and path.name.endswith(SDIST_SUFFIXES)
    ]
    if sdists:
        raise ProvenanceError(
            "wheel-only release directory contains source archives: "
            + ", ".join(path.name for path in sdists)
        )

    cmake_cache = args.cmake_cache.resolve()
    dependency_evidence = [path.resolve() for path in args.dependency_evidence]
    audit_evidence = [path.resolve() for path in args.audit_evidence]
    test_evidence = [path.resolve() for path in args.test_evidence]
    required_files = [
        cmake_cache,
        *dependency_evidence,
        *audit_evidence,
        *test_evidence,
    ]
    missing = [str(path) for path in required_files if not path.is_file()]
    if missing:
        raise ProvenanceError("missing evidence files: " + ", ".join(missing))

    wheel = wheels[0]
    expected_prefix = f"{RELEASE_DISTRIBUTION}-{release_version}-"
    if not wheel.name.startswith(expected_prefix):
        raise ProvenanceError(
            f"wheel filename must start with {expected_prefix!r}: {wheel.name}"
        )
    payload = {
        "schema_version": 1,
        "release_version": release_version,
        "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
        "platform": args.platform,
        "source": checkout,
        "wheel": {
            **evidence_entry(wheel),
            **wheel_metadata(wheel, release_version),
        },
        "cmake_cache": evidence_entry(cmake_cache),
        "dependency_evidence": [
            evidence_entry(path)
            for path in sorted(dependency_evidence, key=lambda path: path.name)
        ],
        "audit_evidence": [
            evidence_entry(path)
            for path in sorted(audit_evidence, key=lambda path: path.name)
        ],
        "test_evidence": [
            evidence_entry(path)
            for path in sorted(test_evidence, key=lambda path: path.name)
        ],
        "ci": {
            "github_run_id": os.environ.get("GITHUB_RUN_ID"),
            "github_run_attempt": os.environ.get("GITHUB_RUN_ATTEMPT"),
            "runner_os": os.environ.get("RUNNER_OS"),
            "runner_arch": os.environ.get("RUNNER_ARCH"),
            "runner_image_os": os.environ.get("ImageOS"),
            "runner_image_version": os.environ.get("ImageVersion"),
        },
    }
    output = evidence_dir / FINAL_RECORD
    write_json(output, payload)
    print(output)


def parser() -> argparse.ArgumentParser:
    source_default = Path(__file__).resolve().parents[1]
    root = argparse.ArgumentParser(description=__doc__)
    subparsers = root.add_subparsers(dest="command", required=True)

    preflight = subparsers.add_parser(
        "preflight",
        help="require a clean checkout and record the exact source commit",
    )
    preflight.add_argument("--source-dir", type=Path, default=source_default)
    preflight.add_argument("--evidence-dir", type=Path, required=True)
    preflight.set_defaults(handler=command_preflight)

    record = subparsers.add_parser(
        "record",
        help="validate the wheel-only output and write its evidence manifest",
    )
    record.add_argument("--source-dir", type=Path, default=source_default)
    record.add_argument("--evidence-dir", type=Path, required=True)
    record.add_argument("--wheel-dir", type=Path, required=True)
    record.add_argument(
        "--platform",
        choices=SUPPORTED_PLATFORMS,
        required=True,
    )
    record.add_argument("--cmake-cache", type=Path, required=True)
    record.add_argument(
        "--dependency-evidence",
        type=Path,
        action="append",
        required=True,
    )
    record.add_argument(
        "--audit-evidence",
        type=Path,
        action="append",
        required=True,
    )
    record.add_argument(
        "--test-evidence",
        type=Path,
        action="append",
        required=True,
    )
    record.set_defaults(handler=command_record)
    return root


def main() -> int:
    args = parser().parse_args()
    try:
        args.handler(args)
    except (OSError, subprocess.CalledProcessError, ProvenanceError) as error:
        print(f"release provenance error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
