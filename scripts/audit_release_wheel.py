#!/usr/bin/env python3
"""Audit a repaired, standalone pypgo wheel."""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
import tempfile
import zipfile
from dataclasses import dataclass
from email.parser import BytesParser
from pathlib import Path, PurePosixPath

from packaging.tags import Tag, parse_tag
from packaging.utils import (
    InvalidWheelFilename,
    canonicalize_name,
    parse_wheel_filename,
)
from packaging.version import InvalidVersion, Version

from project_version import read_project_version
from release_wheel_contract import (
    RELEASE_DISTRIBUTION,
    SUPPORTED_PLATFORMS,
    get_platform_contract,
)


class AuditError(RuntimeError):
    pass


@dataclass(frozen=True)
class WheelArchiveInspection:
    members: tuple[str, ...]
    native_members: tuple[str, ...]
    extension_member: str
    filename_tags: frozenset[Tag]
    metadata_tags: frozenset[Tag]


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


def format_tags(tags: frozenset[Tag]) -> str:
    return ", ".join(sorted(str(tag) for tag in tags)) or "<none>"


def native_component_pattern(platform: str, component: str) -> re.Pattern[str]:
    escaped = re.escape(component)
    if platform == "macos":
        return re.compile(rf"^lib{escaped}(?:\.[0-9]+)*\.dylib$")
    if platform == "linux":
        return re.compile(rf"^lib{escaped}(?:-[0-9a-f]+)?\.so(?:\.[0-9]+)*$")

    windows_basenames = {
        "tbb": "tbb12",
        "gmp": "gmp-10",
        "gmpxx": "gmpxx-4",
        "mpfr": "mpfr-6",
    }
    basename = windows_basenames.get(component, f"{component}.2")
    return re.compile(rf"^{re.escape(basename)}(?:-[0-9a-f]+)?\.dll$")


def missing_native_components(
    platform: str,
    native_members: tuple[str, ...],
) -> list[str]:
    filenames = tuple(
        PurePosixPath(member).name.lower() for member in native_members
    )
    contract = get_platform_contract(platform)
    return [
        component
        for component in contract.required_native_components
        if not any(
            native_component_pattern(platform, component).fullmatch(filename)
            for filename in filenames
        )
    ]


def audit_wheel_archive(
    wheel: Path,
    archive: zipfile.ZipFile,
    platform: str,
    project_version: str,
) -> WheelArchiveInspection:
    contract = get_platform_contract(platform)
    try:
        distribution, version, _, filename_tags = parse_wheel_filename(wheel.name)
        expected_version = Version(project_version)
    except (InvalidVersion, InvalidWheelFilename) as error:
        raise AuditError(f"invalid wheel filename: {wheel.name}: {error}") from error

    require(
        distribution == canonicalize_name(RELEASE_DISTRIBUTION),
        f"unexpected wheel distribution: {distribution}",
    )
    require(
        version == expected_version,
        f"unexpected wheel version: {version}; expected {expected_version}",
    )

    members = tuple(archive.namelist())
    wheel_metadata_members = [
        member for member in members if member.endswith(".dist-info/WHEEL")
    ]
    require(
        len(wheel_metadata_members) == 1,
        "wheel must contain exactly one .dist-info/WHEEL file, "
        f"found: {wheel_metadata_members}",
    )
    wheel_metadata = BytesParser().parsebytes(
        archive.read(wheel_metadata_members[0])
    )
    root_is_purelib = wheel_metadata.get("Root-Is-Purelib")
    require(
        root_is_purelib is not None and root_is_purelib.strip().lower() == "false",
        "repaired wheel must declare Root-Is-Purelib: false",
    )

    metadata_tags: set[Tag] = set()
    try:
        for value in wheel_metadata.get_all("Tag", []):
            metadata_tags.update(parse_tag(value.strip()))
        expected_tags = parse_tag(contract.expected_tag)
    except ValueError as error:
        raise AuditError(f"invalid wheel compatibility tag: {error}") from error
    frozen_metadata_tags = frozenset(metadata_tags)
    require(
        filename_tags == frozen_metadata_tags,
        "wheel filename and WHEEL metadata tags differ: "
        f"filename=[{format_tags(filename_tags)}], "
        f"metadata=[{format_tags(frozen_metadata_tags)}]",
    )
    require(
        filename_tags == expected_tags,
        f"unexpected compatibility tags for {platform}: "
        f"expected=[{format_tags(expected_tags)}], "
        f"actual=[{format_tags(filename_tags)}]",
    )

    forbidden_members = [
        name
        for name in members
        if "/tests/" in f"/{name.lower()}/"
        or "__pycache__" in name.lower()
        or name.lower().endswith((".pyc", ".whl", ".tar.gz"))
    ]
    require(not forbidden_members, f"forbidden wheel members: {forbidden_members}")

    native_members = tuple(
        name
        for name in members
        if name.lower().endswith((".so", ".pyd", ".dylib", ".dll"))
        or ".so." in PurePosixPath(name).name.lower()
    )
    extension_members = [
        name
        for name in native_members
        if PurePosixPath(name).parent.as_posix() == "pypgo"
        and PurePosixPath(name).name.startswith("_pypgo.")
    ]
    require(
        len(extension_members) == 1,
        "wheel must contain exactly one pypgo/_pypgo extension, "
        f"found: {extension_members}",
    )
    require(
        "pypgo/__init__.py" in members,
        "wheel does not contain the pypgo package initializer",
    )
    if platform == "windows":
        initializer = archive.read("pypgo/__init__.py").decode("utf-8")
        require(
            "delvewheel" in initializer and "add_dll_directory" in initializer,
            "Windows wheel initializer does not register its bundled DLL directory",
        )

    unexpected_roots = {
        PurePosixPath(name).parts[0]
        for name in members
        if PurePosixPath(name).parts
        and PurePosixPath(name).parts[0]
        in {"c", "core", "python", "simulationRunner", "tests", "tools"}
    }
    require(
        not unexpected_roots,
        f"wheel contains unintended namespace package roots: {sorted(unexpected_roots)}",
    )

    lowered_members = "\n".join(native_members).lower()
    require("libiomp" not in lowered_members, "Intel OpenMP must not be bundled")
    require(
        "mkl_intel_thread" not in lowered_members,
        "MKL Intel threading runtime must not be bundled",
    )
    require(
        "mkl_gnu_thread" not in lowered_members,
        "MKL GNU threading runtime must not be bundled",
    )

    missing_components = missing_native_components(platform, native_members)
    require(
        not missing_components,
        f"repaired wheel is missing bundled native components: {missing_components}",
    )
    if platform == "windows":
        native_filenames = {
            PurePosixPath(member).name.lower() for member in native_members
        }
        missing_original_runtime_names = sorted(
            set(contract.required_original_runtime_names) - native_filenames
        )
        require(
            not missing_original_runtime_names,
            "Windows wheel is missing original oneMKL/oneTBB runtime names: "
            f"{missing_original_runtime_names}",
        )

    return WheelArchiveInspection(
        members=members,
        native_members=native_members,
        extension_member=extension_members[0],
        filename_tags=filename_tags,
        metadata_tags=frozen_metadata_tags,
    )


def audit(args: argparse.Namespace) -> None:
    wheel = args.wheel.resolve()
    project_version = read_project_version(args.source_dir.resolve())
    require(wheel.is_file(), f"wheel does not exist: {wheel}")

    with zipfile.ZipFile(wheel) as archive:
        inspection = audit_wheel_archive(
            wheel,
            archive,
            args.platform,
            project_version,
        )
        native_members = inspection.native_members

        with tempfile.TemporaryDirectory(prefix="pypgo-wheel-audit-") as directory:
            native = Path(directory) / PurePosixPath(inspection.extension_member).name
            native.write_bytes(archive.read(inspection.extension_member))
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

                require(
                    "accelerate.framework" in lowered,
                    "system Accelerate is not linked",
                )
                require(
                    "mkl" not in lowered and "openblas" not in lowered,
                    "unexpected macOS BLAS provider",
                )
            elif args.platform == "linux":
                require(
                    "x86-64" in run("file", str(native)).lower(),
                    "extension is not x86-64",
                )
                require("mkl" in lowered, "MKL is not present in Linux linkage evidence")
                require("tbb" in lowered, "TBB is not present in Linux linkage evidence")
                require(
                    "mkl_tbb_thread" in lowered,
                    "MKL TBB threading layer is not linked",
                )
                require("libiomp" not in lowered, "Intel OpenMP is linked")
            else:
                headers = run("dumpbin", "/HEADERS", str(native)).lower()
                require(
                    "machine (x64)" in headers or "8664 machine" in headers,
                    "extension is not x64",
                )
                require("mkl" in lowered, "MKL is not present in Windows linkage evidence")
                require("tbb" in lowered, "TBB is not present in Windows linkage evidence")
                require(
                    "tbb12.dll" in lowered and "tbb12-" not in lowered,
                    "Windows extension does not link the original tbb12.dll name",
                )
                require(
                    "mkl_tbb_thread" in lowered,
                    "MKL TBB threading layer is not linked",
                )
                require("libiomp" not in lowered, "Intel OpenMP is linked")

            print(report)

    print(f"audited {wheel.name} for {args.platform}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--platform", choices=SUPPORTED_PLATFORMS, required=True)
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
