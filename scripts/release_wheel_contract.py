#!/usr/bin/env python3
"""Shared release-wheel policy for packaging, auditing, and provenance."""

from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Final, Mapping


RELEASE_DISTRIBUTION: Final = "pypgo"
PYTHON_TAG: Final = "cp312"
ABI_TAG: Final = "cp312"

SUPPORTED_PLATFORMS: Final = ("linux", "macos", "windows")
EVIDENCE_CATEGORIES: Final = ("environment", "audit", "test")

MKL_DISPATCH_COMPONENTS: Final = (
    "def",
    "mc3",
    "avx2",
    "avx512",
    "vml_def",
    "vml_cmpt",
    "vml_mc3",
    "vml_avx2",
    "vml_avx512",
)
MKL_DISPATCH_LIBRARY_COMPONENTS: Final = tuple(
    f"mkl_{component}" for component in MKL_DISPATCH_COMPONENTS
)

WINDOWS_UNMANGLED_RUNTIME_DLLS: Final = (
    "tbb12.dll",
    "mkl_core.2.dll",
    "mkl_tbb_thread.2.dll",
    *(f"mkl_{component}.2.dll" for component in MKL_DISPATCH_COMPONENTS),
)


@dataclass(frozen=True)
class WheelPlatformContract:
    """Immutable release requirements that differ by wheel platform."""

    platform: str
    platform_tag: str
    repair_report: str
    required_native_components: tuple[str, ...]
    required_original_runtime_names: tuple[str, ...] = ()

    @property
    def expected_tag(self) -> str:
        return f"{PYTHON_TAG}-{ABI_TAG}-{self.platform_tag}"

    @property
    def required_audit_evidence(self) -> tuple[str, ...]:
        return (self.repair_report, "wheel-linkage.txt")


_COMMON_NATIVE_COMPONENTS: Final = ("tbb", "gmp", "gmpxx", "mpfr")
_MKL_NATIVE_COMPONENTS: Final = (
    "mkl_core",
    "mkl_tbb_thread",
    *MKL_DISPATCH_LIBRARY_COMPONENTS,
)

PLATFORM_CONTRACTS: Final[Mapping[str, WheelPlatformContract]] = MappingProxyType(
    {
        "macos": WheelPlatformContract(
            platform="macos",
            platform_tag="macosx_26_0_arm64",
            repair_report="delocate-listdeps.txt",
            required_native_components=_COMMON_NATIVE_COMPONENTS,
        ),
        "linux": WheelPlatformContract(
            platform="linux",
            platform_tag="manylinux_2_28_x86_64",
            repair_report="auditwheel.txt",
            required_native_components=(
                *_COMMON_NATIVE_COMPONENTS,
                *_MKL_NATIVE_COMPONENTS,
            ),
        ),
        "windows": WheelPlatformContract(
            platform="windows",
            platform_tag="win_amd64",
            repair_report="delvewheel.txt",
            required_native_components=(
                *_COMMON_NATIVE_COMPONENTS,
                *_MKL_NATIVE_COMPONENTS,
            ),
            required_original_runtime_names=WINDOWS_UNMANGLED_RUNTIME_DLLS,
        ),
    }
)


def get_platform_contract(platform: str) -> WheelPlatformContract:
    """Return the contract for a CLI platform name."""

    try:
        return PLATFORM_CONTRACTS[platform]
    except KeyError as error:
        supported = ", ".join(SUPPORTED_PLATFORMS)
        raise ValueError(
            f"unsupported release-wheel platform {platform!r}; expected one of: "
            f"{supported}"
        ) from error
