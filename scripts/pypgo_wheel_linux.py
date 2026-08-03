#!/usr/bin/env python3
"""Package or audit self-contained Linux pypgo wheels."""

from __future__ import annotations

import os
import re
from collections.abc import Mapping
from pathlib import Path

from pypgo_wheel.audit import require
from pypgo_wheel.common import PackagingError, require_tool, run_capture
from pypgo_wheel.contracts.linux import LINUX_CONTRACT, MANYLINUX_TAG
from pypgo_wheel.linux_mkl import restore
from pypgo_wheel.platform import PypgoWheelPlatform


class LinuxPypgoWheelPlatform(PypgoWheelPlatform):
    contract = LINUX_CONTRACT
    supported_system = "linux"
    supported_machines = frozenset({"x86_64", "amd64"})

    def validate_package_environment(self, source_dir: Path) -> None:
        del source_dir
        actual_tag = os.environ.get("AUDITWHEEL_PLAT")
        if actual_tag != MANYLINUX_TAG:
            raise PackagingError(
                f"Linux packaging requires the {MANYLINUX_TAG} container; "
                f"AUDITWHEEL_PLAT={actual_tag!r}"
            )
        require_tool("auditwheel")

    def configure_build_environment(self, source_dir: Path) -> dict[str, str]:
        del source_dir
        environment = os.environ.copy()
        python_library_dir = str(Path(os.sys.prefix).resolve() / "lib")
        current = environment.get("LD_LIBRARY_PATH")
        environment["LD_LIBRARY_PATH"] = (
            f"{python_library_dir}{os.pathsep}{current}"
            if current
            else python_library_dir
        )
        return environment

    def repair_wheel(
        self,
        raw_wheel: Path,
        repair_dir: Path,
        source_dir: Path,
        environment: Mapping[str, str],
    ) -> None:
        auditwheel = require_tool("auditwheel")
        run_capture(
            (
                auditwheel,
                "repair",
                "--plat",
                MANYLINUX_TAG,
                "--only-plat",
                "--wheel-dir",
                repair_dir,
                raw_wheel,
            ),
            cwd=source_dir,
            env=environment,
        )

    def post_repair(self, repaired_wheel: Path, source_dir: Path) -> None:
        del source_dir
        restore(repaired_wheel)

    def collect_repair_report(self, repaired_wheel: Path, source_dir: Path) -> str:
        return run_capture(
            (require_tool("auditwheel"), "show", repaired_wheel),
            cwd=source_dir,
        )

    def native_component_pattern(self, component: str) -> re.Pattern[str]:
        return re.compile(
            rf"^lib{re.escape(component)}(?:-[0-9a-f]+)?\.so(?:\.[0-9]+)*$"
        )

    def audit_native_binaries(
        self,
        extension: Path,
        natives: tuple[Path, ...],
    ) -> str:
        del natives
        report = run_capture((require_tool("readelf"), "-d", extension))
        report += run_capture((require_tool("ldd"), extension))
        file_report = run_capture((require_tool("file"), extension))
        lowered = report.lower()
        require("x86-64" in file_report.lower(), "extension is not x86-64")
        require("mkl" in lowered, "MKL is not present in Linux linkage evidence")
        require("tbb" in lowered, "TBB is not present in Linux linkage evidence")
        require(
            "mkl_tbb_thread" in lowered,
            "MKL TBB threading layer is not linked",
        )
        require("libiomp" not in lowered, "Intel OpenMP is linked")
        return file_report + report


if __name__ == "__main__":
    raise SystemExit(LinuxPypgoWheelPlatform().main())
