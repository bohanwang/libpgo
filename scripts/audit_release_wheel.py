#!/usr/bin/env python3
"""Audit a repaired, standalone pypgo wheel."""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

from project_version import read_project_version


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
    project_version = read_project_version(args.source_dir.resolve())
    require(wheel.is_file(), f"wheel does not exist: {wheel}")
    require(
        wheel.name.startswith(f"pypgo-{project_version}-"),
        f"unexpected wheel: {wheel.name}",
    )

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
            or ".so." in Path(name).name.lower()
        ]
        extension_members = [
            name
            for name in native_members
            if Path(name).name.startswith("pypgo.")
        ]
        require(
            len(extension_members) == 1,
            f"wheel must contain exactly one pypgo extension, found: {extension_members}",
        )

        lowered_members = "\n".join(native_members).lower()
        require("libiomp" not in lowered_members, "Intel OpenMP must not be bundled")
        require("mkl_intel_thread" not in lowered_members, "MKL Intel threading runtime must not be bundled")
        require("mkl_gnu_thread" not in lowered_members, "MKL GNU threading runtime must not be bundled")

        if args.platform in {"linux", "windows"}:
            require("mkl" in lowered_members, "repaired wheel does not bundle MKL")
            require("tbb" in lowered_members, "repaired wheel does not bundle TBB")
            require(
                "mkl_tbb_thread" in lowered_members,
                "repaired wheel does not bundle the MKL TBB threading layer",
            )
            dispatch_components = (
                "mkl_def",
                "mkl_mc3",
                "mkl_avx2",
                "mkl_avx512",
                "mkl_vml_def",
                "mkl_vml_cmpt",
                "mkl_vml_mc3",
                "mkl_vml_avx2",
                "mkl_vml_avx512",
            )
            missing_dispatch = [
                component
                for component in dispatch_components
                if component not in lowered_members
            ]
            require(
                not missing_dispatch,
                f"repaired wheel is missing MKL dispatch libraries: {missing_dispatch}",
            )
            if args.platform == "windows":
                native_filenames = {Path(member).name.lower() for member in native_members}
                required_original_mkl_names = {
                    "mkl_core.2.dll",
                    "mkl_tbb_thread.2.dll",
                }
                required_original_mkl_names.update(
                    f"{component}.2.dll" for component in dispatch_components
                )
                missing_original_mkl_names = sorted(
                    required_original_mkl_names - native_filenames
                )
                require(
                    not missing_original_mkl_names,
                    "Windows wheel is missing original oneMKL runtime names: "
                    f"{missing_original_mkl_names}",
                )
        else:
            require("tbb" in lowered_members, "repaired macOS wheel does not bundle TBB")

        with tempfile.TemporaryDirectory(prefix="pypgo-wheel-audit-") as directory:
            native = Path(directory) / Path(extension_members[0]).name
            native.write_bytes(archive.read(extension_members[0]))
            report = dependency_report(native, args.platform)

            forbidden_paths = {
                str(args.source_dir.resolve()),
                str(sys.prefix),
                "/opt/homebrew",
                "/usr/local",
            }
            leaked = sorted(path for path in forbidden_paths if path and path in report)
            require(not leaked, f"forbidden dependency or RPATH entries: {leaked}")

            lowered = report.lower()
            if args.platform == "macos":
                macos_natives = []
                for member in native_members:
                    extracted = Path(directory) / Path(member).name
                    extracted.write_bytes(archive.read(member))
                    macos_natives.append(extracted)

                for macos_native in macos_natives:
                    architectures = run("lipo", "-archs", str(macos_native)).split()
                    require(
                        architectures == ["arm64"],
                        f"unexpected architectures in {macos_native.name}: {architectures}",
                    )
                    build_info = run("otool", "-l", str(macos_native))
                    min_versions = re.findall(r"\bminos ([0-9.]+)", build_info)
                    require(
                        min_versions,
                        f"missing LC_BUILD_VERSION minimum OS in {macos_native.name}",
                    )
                    require(
                        all(
                            tuple(map(int, version.split("."))) <= (26, 0)
                            for version in min_versions
                        ),
                        f"{macos_native.name} requires newer than macOS 26.0: "
                        f"{min_versions}",
                    )

                require("accelerate.framework" in lowered, "system Accelerate is not linked")
                require("mkl" not in lowered and "openblas" not in lowered, "unexpected macOS BLAS provider")
            elif args.platform == "linux":
                require("x86-64" in run("file", str(native)).lower(), "extension is not x86-64")
                require("mkl" in lowered, "MKL is not present in Linux linkage evidence")
                require("tbb" in lowered, "TBB is not present in Linux linkage evidence")
                require("mkl_tbb_thread" in lowered, "MKL TBB threading layer is not linked")
                require("libiomp" not in lowered, "Intel OpenMP is linked")
            else:
                headers = run("dumpbin", "/HEADERS", str(native)).lower()
                require("machine (x64)" in headers or "8664 machine" in headers, "extension is not x64")
                require("mkl" in lowered, "MKL is not present in Windows linkage evidence")
                require("tbb" in lowered, "TBB is not present in Windows linkage evidence")
                require("mkl_tbb_thread" in lowered, "MKL TBB threading layer is not linked")
                require("libiomp" not in lowered, "Intel OpenMP is linked")

            print(report)

    print(f"audited {wheel.name} for {args.platform}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--platform", choices=("linux", "macos", "windows"), required=True)
    parser.add_argument("--wheel", type=Path, required=True)
    parser.add_argument("--source-dir", type=Path, required=True)
    args = parser.parse_args()
    try:
        audit(args)
    except (AuditError, OSError, subprocess.CalledProcessError, zipfile.BadZipFile) as error:
        print(f"wheel audit error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
