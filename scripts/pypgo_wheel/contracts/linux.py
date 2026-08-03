"""Linux release-wheel contract."""

from __future__ import annotations

from typing import Final

from .common import (
    COMMON_NATIVE_COMPONENTS,
    WheelPlatformContract,
)
from .mkl import MKL_NATIVE_COMPONENTS


MANYLINUX_TAG: Final = "manylinux_2_28_x86_64"

LINUX_CONTRACT: Final = WheelPlatformContract(
    platform="linux",
    platform_tag=MANYLINUX_TAG,
    repair_report="auditwheel.txt",
    required_native_components=(
        *COMMON_NATIVE_COMPONENTS,
        *MKL_NATIVE_COMPONENTS,
    ),
)
