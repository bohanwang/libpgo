"""Platform-independent pypgo wheel archive and native linkage auditing."""

from __future__ import annotations

import sys
import tempfile
import zipfile
from dataclasses import dataclass
from email.parser import BytesParser
from pathlib import Path, PurePosixPath
from typing import Protocol

from packaging.tags import Tag, parse_tag
from packaging.utils import (
    InvalidWheelFilename,
    canonicalize_name,
    parse_wheel_filename,
)
from packaging.version import InvalidVersion, Version

from project_version import read_project_version
from .common import RELEASE_DISTRIBUTION


class AuditError(RuntimeError):
    """Raised when a repaired wheel violates the release contract."""


@dataclass(frozen=True)
class WheelArchiveInspection:
    members: tuple[str, ...]
    native_members: tuple[str, ...]
    extension_member: str
    filename_tags: frozenset[Tag]
    metadata_tags: frozenset[Tag]


class WheelAuditPlatform(Protocol):
    platform: str
    platform_tag: str
    required_native_components: tuple[str, ...]
    required_original_runtime_names: tuple[str, ...]

    @property
    def expected_tag(self) -> str: ...

    def native_component_pattern(self, component: str): ...

    def audit_platform_archive(
        self,
        archive: zipfile.ZipFile,
        inspection: WheelArchiveInspection,
    ) -> None: ...

    def audit_native_binaries(
        self,
        extension: Path,
        natives: tuple[Path, ...],
    ) -> str: ...


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AuditError(message)


def format_tags(tags: frozenset[Tag]) -> str:
    return ", ".join(sorted(str(tag) for tag in tags)) or "<none>"


def missing_native_components(
    platform: WheelAuditPlatform,
    native_members: tuple[str, ...],
) -> list[str]:
    filenames = tuple(PurePosixPath(member).name.lower() for member in native_members)
    return [
        component
        for component in platform.required_native_components
        if not any(
            platform.native_component_pattern(component).fullmatch(filename)
            for filename in filenames
        )
    ]


def audit_wheel_archive(
    wheel: Path,
    archive: zipfile.ZipFile,
    platform: WheelAuditPlatform,
    project_version: str,
) -> WheelArchiveInspection:
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
    wheel_metadata = BytesParser().parsebytes(archive.read(wheel_metadata_members[0]))
    root_is_purelib = wheel_metadata.get("Root-Is-Purelib")
    require(
        root_is_purelib is not None and root_is_purelib.strip().lower() == "false",
        "repaired wheel must declare Root-Is-Purelib: false",
    )

    metadata_tags: set[Tag] = set()
    try:
        for value in wheel_metadata.get_all("Tag", []):
            metadata_tags.update(parse_tag(value.strip()))
        expected_tags = parse_tag(platform.expected_tag)
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
        f"unexpected compatibility tags for {platform.platform}: "
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

    inspection = WheelArchiveInspection(
        members=members,
        native_members=native_members,
        extension_member=extension_members[0],
        filename_tags=filename_tags,
        metadata_tags=frozen_metadata_tags,
    )
    platform.audit_platform_archive(archive, inspection)
    return inspection


def audit_repaired_wheel(
    wheel: Path,
    source_dir: Path,
    platform: WheelAuditPlatform,
) -> str:
    wheel = wheel.resolve()
    source_dir = source_dir.resolve()
    require(wheel.is_file(), f"wheel does not exist: {wheel}")
    project_version = read_project_version(source_dir)

    with zipfile.ZipFile(wheel) as archive:
        inspection = audit_wheel_archive(
            wheel,
            archive,
            platform,
            project_version,
        )
        with tempfile.TemporaryDirectory(prefix="pypgo-wheel-audit-") as directory:
            temporary = Path(directory)
            extracted_natives = []
            for index, member in enumerate(inspection.native_members):
                extracted = temporary / f"{index}-{PurePosixPath(member).name}"
                extracted.write_bytes(archive.read(member))
                extracted_natives.append(extracted)
            extension_index = inspection.native_members.index(
                inspection.extension_member
            )
            extension = extracted_natives[extension_index]
            report = platform.audit_native_binaries(
                extension,
                tuple(extracted_natives),
            )

    forbidden_paths = {
        str(source_dir),
        str(Path(sys.prefix).resolve()),
        "/opt/homebrew",
        "/usr/local",
    }
    leaked = sorted(path for path in forbidden_paths if path and path in report)
    require(not leaked, f"forbidden dependency or RPATH entries: {leaked}")
    return report
