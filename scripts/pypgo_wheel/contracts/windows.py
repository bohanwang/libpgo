"""Windows release-wheel contract and DLL filename policy."""

from __future__ import annotations

from typing import Final

from .common import (
    COMMON_NATIVE_COMPONENTS,
    WheelPlatformContract,
)
from .mkl import MKL_DISPATCH_COMPONENTS, MKL_NATIVE_COMPONENTS


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

WINDOWS_CONTRACT: Final = WheelPlatformContract(
    platform="windows",
    platform_tag="win_amd64",
    repair_report="delvewheel.txt",
    required_native_components=(
        *COMMON_NATIVE_COMPONENTS,
        *MKL_NATIVE_COMPONENTS,
    ),
    required_original_runtime_names=WINDOWS_UNMANGLED_RUNTIME_DLLS,
)
