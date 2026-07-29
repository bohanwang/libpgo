#!/usr/bin/env python3
"""Audit a final pypgo wheel without repairing or vendoring dependencies."""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path


class AuditError(RuntimeError):
    pass


def run(*command: str) -> str:
    result = subprocess.run(
        command,
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AuditError(message)


def dependency_report(native: Path, target: str) -> str:
    if target == "macos":
        return run("otool", "-L", str(native)) + run("otool", "-l", str(native))
    if target == "linux":
        return run("readelf", "-d", str(native)) + run("ldd", str(native))
    return run("dumpbin", "/DEPENDENTS", str(native))


def audit(args: argparse.Namespace) -> None:
    wheel = args.wheel.resolve()
    require(wheel.is_file(), f"wheel does not exist: {wheel}")
    require(wheel.name.startswith("pypgo-0.0.4-"), f"unexpected wheel: {wheel.name}")

    with zipfile.ZipFile(wheel) as archive:
        members = archive.namelist()
        forbidden_members = [
            name
            for name in members
            if "/tests/" in f"/{name.lower()}/"
            or "__pycache__" in name.lower()
            or name.lower().endswith((".pyc", ".whl", ".tar.gz"))
        ]
        require(not forbidden_members, f"forbidden wheel members: {forbidden_members}")

        native_members = [
            name
            for name in members
            if name.lower().endswith((".so", ".pyd", ".dylib", ".dll"))
        ]
        require(
            len(native_members) == 1,
            f"wheel must contain only the pypgo extension, found: {native_members}",
        )
        require(
            Path(native_members[0]).name.startswith("pypgo."),
            f"unexpected native wheel member: {native_members[0]}",
        )

        with tempfile.TemporaryDirectory(prefix="pypgo-wheel-audit-") as directory:
            native = Path(directory) / Path(native_members[0]).name
            native.write_bytes(archive.read(native_members[0]))
            report = dependency_report(native, args.platform)

            forbidden_paths = {
                str(args.source_dir.resolve()),
                str(args.conda_prefix.resolve()),
                "/opt/homebrew",
                "/usr/local",
            }
            leaked = sorted(path for path in forbidden_paths if path and path in report)
            require(not leaked, f"forbidden dependency or RPATH entries: {leaked}")

            lowered = report.lower()
            if args.platform == "macos":
                architectures = run("lipo", "-archs", str(native)).split()
                require(architectures == ["arm64"], f"unexpected architectures: {architectures}")
                require("accelerate.framework" in lowered, "system Accelerate is not linked")
                require("mkl" not in lowered and "openblas" not in lowered, "unexpected macOS BLAS provider")
                build_info = run("otool", "-l", str(native))
                min_versions = re.findall(r"\\bminos ([0-9.]+)", build_info)
                require(min_versions, "missing LC_BUILD_VERSION minimum OS")
                require(
                    all(tuple(map(int, version.split("."))) <= (14, 0) for version in min_versions),
                    f"wheel requires newer than macOS 14.0: {min_versions}",
                )
            elif args.platform == "linux":
                require("x86-64" in run("file", str(native)).lower(), "extension is not x86-64")
                require("mkl" in lowered, "MKL is not present in Linux linkage evidence")
                require("tbb" in lowered, "TBB is not present in Linux linkage evidence")
                require("libiomp" not in lowered, "Intel OpenMP is linked")
            else:
                headers = run("dumpbin", "/HEADERS", str(native)).lower()
                require("machine (x64)" in headers or "8664 machine" in headers, "extension is not x64")
                require("mkl" in lowered, "MKL is not present in Windows linkage evidence")
                require("tbb" in lowered, "TBB is not present in Windows linkage evidence")
                require("libiomp" not in lowered, "Intel OpenMP is linked")

            print(report)

    print(f"audited {wheel.name} for {args.platform}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--platform", choices=("linux", "macos", "windows"), required=True)
    parser.add_argument("--wheel", type=Path, required=True)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument(
        "--conda-prefix",
        type=Path,
        default=Path(os.environ.get("CONDA_PREFIX", sys.prefix)),
    )
    args = parser.parse_args()
    try:
        audit(args)
    except (AuditError, OSError, subprocess.CalledProcessError, zipfile.BadZipFile) as error:
        print(f"wheel audit error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
