"""macOS release-wheel contract."""

from __future__ import annotations

from typing import Final

from .common import COMMON_NATIVE_COMPONENTS, WheelPlatformContract


MACOS_DEPLOYMENT_TARGET: Final = "26.0"
MACOS_ARCHITECTURE: Final = "arm64"

MACOS_CONTRACT: Final = WheelPlatformContract(
    platform="macos",
    platform_tag=(
        f"macosx_{MACOS_DEPLOYMENT_TARGET.replace('.', '_')}_{MACOS_ARCHITECTURE}"
    ),
    repair_report="delocate-listdeps.txt",
    required_native_components=COMMON_NATIVE_COMPONENTS,
)
