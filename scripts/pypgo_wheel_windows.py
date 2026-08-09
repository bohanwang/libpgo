#!/usr/bin/env python3
"""Package or audit self-contained Windows pypgo wheels."""

from __future__ import annotations

import os
import re
import zipfile
from collections.abc import Mapping
from pathlib import Path, PurePosixPath
from typing import Final

from pypgo_wheel.audit import WheelArchiveInspection, require
from pypgo_wheel.common import (
    COMMON_NATIVE_COMPONENTS,
    MKL_DISPATCH_COMPONENTS,
    MKL_NATIVE_COMPONENTS,
    require_tool,
    run_capture,
)

from pypgo_wheel.platform import PypgoWheelPlatform


WINDOWS_UNMANGLED_RUNTIME_DLLS: Final = (
    "tbb12.dll",
    "mkl_core.2.dll",
    "mkl_tbb_thread.2.dll",
    *(f"mkl_{component}.2.dll" for component in MKL_DISPATCH_COMPONENTS),
)

WINDOWS_FORCED_INCLUDE_DLLS: Final = (
    "gmpxx-4.dll",
    "mpfr-6.dll",
    *WINDOWS_UNMANGLED_RUNTIME_DLLS,
)


class WindowsPypgoWheelPlatform(PypgoWheelPlatform):
    platform = "windows"
    platform_tag = "win_amd64"
    repair_report = "delvewheel.txt"
    required_native_components = (*COMMON_NATIVE_COMPONENTS, *MKL_NATIVE_COMPONENTS)
    required_original_runtime_names = WINDOWS_UNMANGLED_RUNTIME_DLLS
    supported_system = "win32"
    supported_machines = frozenset({"amd64", "x86_64"})

    def __init__(self) -> None:
        self.native_directories: tuple[Path, ...] = ()

    def validate_package_environment(self, source_dir: Path) -> None:
        require_tool("cl")
        self.native_directories = (
            Path(os.sys.prefix).resolve() / "Library" / "bin",
            source_dir / "third-party" / "gmp-msvc" / "release",
            source_dir / "third-party" / "mpfr-msvc" / "release",
        )
        missing = [
            directory for directory in self.native_directories if not directory.is_dir()
        ]
        if missing:
            raise FileNotFoundError(
                f"missing Windows native runtime directories: {missing}"
            )

    def repair_wheel(
        self,
        raw_wheel: Path,
        repair_dir: Path,
        source_dir: Path,
        environment: Mapping[str, str],
    ) -> None:
        native_paths = os.pathsep.join(
            str(directory) for directory in self.native_directories
        )
        included = os.pathsep.join(WINDOWS_FORCED_INCLUDE_DLLS)
        unmangled = os.pathsep.join(WINDOWS_UNMANGLED_RUNTIME_DLLS)
        run_capture(
            (
                os.sys.executable,
                "-m",
                "delvewheel",
                "repair",
                "--add-path",
                native_paths,
                "--include",
                included,
                "--no-mangle",
                unmangled,
                "--wheel-dir",
                repair_dir,
                raw_wheel,
            ),
            cwd=source_dir,
            env=environment,
        )

    def collect_repair_report(self, repaired_wheel: Path, source_dir: Path) -> str:
        return run_capture(
            (os.sys.executable, "-m", "delvewheel", "show", repaired_wheel),
            cwd=source_dir,
        )

    def native_component_pattern(self, component: str) -> re.Pattern[str]:
        basenames = {
            "tbb": "tbb12",
            "gmp": "gmp-10",
            "gmpxx": "gmpxx-4",
            "mpfr": "mpfr-6",
        }
        basename = basenames.get(component, f"{component}.2")
        return re.compile(rf"^{re.escape(basename)}(?:-[0-9a-f]+)?\.dll$")

    def audit_platform_archive(
        self,
        archive: zipfile.ZipFile,
        inspection: WheelArchiveInspection,
    ) -> None:
        initializer = archive.read("pypgo/__init__.py").decode("utf-8")
        require(
            "delvewheel" in initializer and "add_dll_directory" in initializer,
            "Windows wheel initializer does not register its bundled DLL directory",
        )
        native_filenames = {
            PurePosixPath(member).name.lower()
            for member in inspection.native_members
        }
        missing_original_names = sorted(
            set(self.required_original_runtime_names)
            - native_filenames
        )
        require(
            not missing_original_names,
            "Windows wheel is missing original oneMKL/oneTBB runtime names: "
            f"{missing_original_names}",
        )

    def audit_native_binaries(
        self,
        extension: Path,
        natives: tuple[Path, ...],
    ) -> str:
        del natives
        dumpbin = require_tool("dumpbin")
        report = run_capture((dumpbin, "/DEPENDENTS", extension))
        headers = run_capture((dumpbin, "/HEADERS", extension))
        lowered = report.lower()
        lowered_headers = headers.lower()
        require(
            "machine (x64)" in lowered_headers or "8664 machine" in lowered_headers,
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
        return headers + report


if __name__ == "__main__":
    raise SystemExit(WindowsPypgoWheelPlatform().main())
