"""Template-method base class for pypgo wheel packaging and auditing."""

from __future__ import annotations

import argparse
import os
import platform as host_platform
import re
import subprocess
import sys
import tempfile
import zipfile
from abc import ABC, abstractmethod
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import ClassVar

from .audit import AuditError, WheelArchiveInspection, audit_repaired_wheel
from .common import (
    ABI_TAG,
    PYTHON_TAG,
    PackagingError,
    build_raw_wheel,
    exactly_one_wheel,
    prepare_wheel_dir,
    publish_wheel,
    write_report,
)


SOURCE_ROOT = Path(__file__).resolve().parents[2]


class PypgoWheelPlatform(ABC):
    """Define platform hooks around one shared package/audit workflow."""

    platform: ClassVar[str]
    platform_tag: ClassVar[str]
    repair_report: ClassVar[str]
    required_native_components: ClassVar[tuple[str, ...]]
    required_original_runtime_names: ClassVar[tuple[str, ...]] = ()
    supported_system: ClassVar[str]
    supported_machines: ClassVar[frozenset[str]]

    @property
    def expected_tag(self) -> str:
        return f"{PYTHON_TAG}-{ABI_TAG}-{self.platform_tag}"

    def validate_host(self) -> None:
        if sys.version_info[:2] != (3, 12):
            raise PackagingError(
                "pypgo wheels require CPython 3.12; "
                f"running {sys.version.split()[0]}"
            )
        machine = host_platform.machine().lower()
        if sys.platform != self.supported_system or machine not in self.supported_machines:
            raise PackagingError(
                f"{self.platform} wheel tooling requires {self.supported_system} "
                f"{sorted(self.supported_machines)}; running {sys.platform} {machine}"
            )

    def validate_package_environment(self, source_dir: Path) -> None:
        """Perform packaging-only environment checks after host validation."""

    def configure_build_environment(self, source_dir: Path) -> dict[str, str]:
        return os.environ.copy()

    @abstractmethod
    def repair_wheel(
        self,
        raw_wheel: Path,
        repair_dir: Path,
        source_dir: Path,
        environment: Mapping[str, str],
    ) -> None:
        """Repair one raw wheel into repair_dir."""

    def post_repair(self, repaired_wheel: Path, source_dir: Path) -> None:
        """Apply platform-specific changes after the repair tool finishes."""

    @abstractmethod
    def collect_repair_report(self, repaired_wheel: Path, source_dir: Path) -> str:
        """Return the platform repair tool's report."""

    @abstractmethod
    def native_component_pattern(self, component: str) -> re.Pattern[str]:
        """Return the accepted archive filename pattern for a native component."""

    def audit_platform_archive(
        self,
        archive: zipfile.ZipFile,
        inspection: WheelArchiveInspection,
    ) -> None:
        """Validate additional platform-specific archive properties."""

    @abstractmethod
    def audit_native_binaries(
        self,
        extension: Path,
        natives: tuple[Path, ...],
    ) -> str:
        """Audit extracted native binaries and return linkage evidence."""

    def audit(self, wheel: Path, source_dir: Path = SOURCE_ROOT) -> str:
        self.validate_host()
        return audit_repaired_wheel(wheel, source_dir, self)

    def package(
        self,
        source_dir: Path,
        wheel_dir: Path,
        report_dir: Path | None = None,
    ) -> Path:
        self.validate_host()
        source_dir = source_dir.resolve()
        self.validate_package_environment(source_dir)
        wheel_dir = prepare_wheel_dir(wheel_dir)
        environment = self.configure_build_environment(source_dir)

        with tempfile.TemporaryDirectory(prefix="pypgo-wheel-package-") as directory:
            temporary = Path(directory)
            raw_wheel = build_raw_wheel(source_dir, temporary / "raw", environment)
            repair_dir = temporary / "repaired"
            repair_dir.mkdir()
            self.repair_wheel(
                raw_wheel,
                repair_dir,
                source_dir,
                environment,
            )
            repaired_wheel = exactly_one_wheel(repair_dir)
            self.post_repair(repaired_wheel, source_dir)
            repair_report = self.collect_repair_report(repaired_wheel, source_dir)
            write_report(
                report_dir,
                self.repair_report,
                repair_report,
            )
            audit_report = self.audit(repaired_wheel, source_dir)
            write_report(report_dir, "wheel-linkage.txt", audit_report)
            return publish_wheel(repaired_wheel, wheel_dir)

    def main(self, argv: Sequence[str] | None = None) -> int:
        parser = argparse.ArgumentParser(
            description=f"Package or audit self-contained pypgo {self.platform} wheels."
        )
        subparsers = parser.add_subparsers(dest="command", required=True)

        package_parser = subparsers.add_parser(
            "package",
            help="build, repair, audit, and publish one wheel",
        )
        package_parser.add_argument("--wheel-dir", type=Path, required=True)
        package_parser.add_argument("--report-dir", type=Path)
        package_parser.add_argument("--source-dir", type=Path, default=SOURCE_ROOT)

        audit_parser = subparsers.add_parser(
            "audit",
            help="audit an existing repaired wheel",
        )
        audit_parser.add_argument("--wheel", type=Path, required=True)
        audit_parser.add_argument("--source-dir", type=Path, default=SOURCE_ROOT)

        arguments = parser.parse_args(argv)
        try:
            if arguments.command == "package":
                wheel = self.package(
                    arguments.source_dir,
                    arguments.wheel_dir,
                    arguments.report_dir,
                )
                print(f"self-contained wheel: {wheel}")
            else:
                self.audit(arguments.wheel, arguments.source_dir)
                print(f"audited {arguments.wheel.name} for {self.platform}")
        except (
            AuditError,
            PackagingError,
            OSError,
            subprocess.CalledProcessError,
            zipfile.BadZipFile,
        ) as error:
            print(f"wheel {arguments.command} error: {error}", file=sys.stderr)
            return 1
        return 0
