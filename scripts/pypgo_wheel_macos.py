#!/usr/bin/env python3
"""Package or audit self-contained macOS pypgo wheels."""

from __future__ import annotations

import os
import re
from collections.abc import Mapping
from pathlib import Path

from packaging.version import Version

from pypgo_wheel.audit import require
from pypgo_wheel.common import require_tool, run_capture
from pypgo_wheel.contracts.macos import (
    MACOS_ARCHITECTURE,
    MACOS_CONTRACT,
    MACOS_DEPLOYMENT_TARGET,
)
from pypgo_wheel.platform import PypgoWheelPlatform


class MacOSPypgoWheelPlatform(PypgoWheelPlatform):
    contract = MACOS_CONTRACT
    supported_system = "darwin"
    supported_machines = frozenset({"arm64", "aarch64"})

    def validate_package_environment(self, source_dir: Path) -> None:
        del source_dir
        require_tool("delocate-wheel")
        require_tool("delocate-listdeps")

    def configure_build_environment(self, source_dir: Path) -> dict[str, str]:
        del source_dir
        environment = os.environ.copy()
        environment.update(
            {
                "ARCHFLAGS": f"-arch {MACOS_ARCHITECTURE}",
                "MACOSX_DEPLOYMENT_TARGET": MACOS_DEPLOYMENT_TARGET,
                "_PYTHON_HOST_PLATFORM": (
                    f"macosx-{MACOS_DEPLOYMENT_TARGET}-{MACOS_ARCHITECTURE}"
                ),
            }
        )
        return environment

    def repair_wheel(
        self,
        raw_wheel: Path,
        repair_dir: Path,
        source_dir: Path,
        environment: Mapping[str, str],
    ) -> None:
        run_capture(
            (
                require_tool("delocate-wheel"),
                "--require-archs",
                MACOS_ARCHITECTURE,
                "--wheel-dir",
                repair_dir,
                raw_wheel,
            ),
            cwd=source_dir,
            env=environment,
        )

    def collect_repair_report(self, repaired_wheel: Path, source_dir: Path) -> str:
        return run_capture(
            (require_tool("delocate-listdeps"), "--all", repaired_wheel),
            cwd=source_dir,
        )

    def native_component_pattern(self, component: str) -> re.Pattern[str]:
        return re.compile(rf"^lib{re.escape(component)}(?:\.[0-9]+)*\.dylib$")

    def audit_native_binaries(
        self,
        extension: Path,
        natives: tuple[Path, ...],
    ) -> str:
        otool = require_tool("otool")
        lipo = require_tool("lipo")
        report = run_capture((otool, "-L", extension))
        report += run_capture((otool, "-l", extension))
        for native in natives:
            architectures = run_capture((lipo, "-archs", native)).split()
            require(
                architectures == [MACOS_ARCHITECTURE],
                f"unexpected architectures in {native.name}: {architectures}",
            )
            build_info = run_capture((otool, "-l", native))
            report += build_info
            min_versions = re.findall(r"\bminos ([0-9.]+)", build_info)
            require(
                min_versions,
                f"missing LC_BUILD_VERSION minimum OS in {native.name}",
            )
            require(
                all(
                    Version(version) <= Version(MACOS_DEPLOYMENT_TARGET)
                    for version in min_versions
                ),
                f"{native.name} requires newer than macOS "
                f"{MACOS_DEPLOYMENT_TARGET}: {min_versions}",
            )

        lowered = report.lower()
        require("accelerate.framework" in lowered, "system Accelerate is not linked")
        require(
            "mkl" not in lowered and "openblas" not in lowered,
            "unexpected macOS BLAS provider",
        )
        return report


if __name__ == "__main__":
    raise SystemExit(MacOSPypgoWheelPlatform().main())
