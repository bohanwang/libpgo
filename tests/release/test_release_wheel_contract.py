import sys
from pathlib import Path

import pytest


SCRIPTS_DIR = Path(__file__).resolve().parents[2] / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

from pypgo_wheel.common import (  # noqa: E402
    ABI_TAG,
    COMMON_NATIVE_COMPONENTS,
    MKL_DISPATCH_COMPONENTS,
    MKL_DISPATCH_LIBRARY_COMPONENTS,
    PYTHON_TAG,
    RELEASE_DISTRIBUTION,
    SUPPORTED_PLATFORMS,
)
from pypgo_wheel_linux import LinuxPypgoWheelPlatform  # noqa: E402
from pypgo_wheel_macos import MacOSPypgoWheelPlatform  # noqa: E402
from pypgo_wheel_windows import (  # noqa: E402
    WINDOWS_FORCED_INCLUDE_DLLS,
    WINDOWS_UNMANGLED_RUNTIME_DLLS,
    WindowsPypgoWheelPlatform,
)


PLATFORM_CLASSES = {
    "linux": LinuxPypgoWheelPlatform,
    "macos": MacOSPypgoWheelPlatform,
    "windows": WindowsPypgoWheelPlatform,
}


def test_shared_release_identity() -> None:
    assert RELEASE_DISTRIBUTION == "pypgo"
    assert PYTHON_TAG == "cp312"
    assert ABI_TAG == "cp312"


@pytest.mark.parametrize(
    ("platform", "expected_tag", "repair_report"),
    [
        ("macos", "cp312-cp312-macosx_26_0_arm64", "delocate-listdeps.txt"),
        ("linux", "cp312-cp312-manylinux_2_28_x86_64", "auditwheel.txt"),
        ("windows", "cp312-cp312-win_amd64", "delvewheel.txt"),
    ],
)
def test_platform_release_facts(
    platform: str,
    expected_tag: str,
    repair_report: str,
) -> None:
    instance = PLATFORM_CLASSES[platform]()

    assert instance.platform == platform
    assert instance.expected_tag == expected_tag
    assert instance.repair_report == repair_report
    assert {"tbb", "gmp", "gmpxx", "mpfr"}.issubset(
        instance.required_native_components
    )


def test_supported_platforms_match_platform_classes() -> None:
    assert set(SUPPORTED_PLATFORMS) == set(PLATFORM_CLASSES)


def test_every_platform_bundles_common_components() -> None:
    assert all(
        tuple(PLATFORM_CLASSES[name]().required_native_components[:4])
        == COMMON_NATIVE_COMPONENTS
        for name in SUPPORTED_PLATFORMS
    )


def test_mkl_dispatch_names_are_derived_once() -> None:
    assert MKL_DISPATCH_LIBRARY_COMPONENTS == tuple(
        f"mkl_{component}" for component in MKL_DISPATCH_COMPONENTS
    )
    assert WINDOWS_UNMANGLED_RUNTIME_DLLS == (
        "tbb12.dll",
        "mkl_core.2.dll",
        "mkl_tbb_thread.2.dll",
        *(f"mkl_{component}.2.dll" for component in MKL_DISPATCH_COMPONENTS),
    )


def test_windows_original_runtime_names_are_declared_on_the_platform() -> None:
    windows = PLATFORM_CLASSES["windows"]()

    assert windows.required_original_runtime_names == WINDOWS_UNMANGLED_RUNTIME_DLLS
    assert not PLATFORM_CLASSES["linux"]().required_original_runtime_names
    assert not PLATFORM_CLASSES["macos"]().required_original_runtime_names


def test_windows_forced_includes_cover_indirect_and_unmangled_runtimes() -> None:
    assert WINDOWS_FORCED_INCLUDE_DLLS == (
        "gmpxx-4.dll",
        "mpfr-6.dll",
        *WINDOWS_UNMANGLED_RUNTIME_DLLS,
    )


def test_expected_tag_is_derived_from_shared_and_platform_tags() -> None:
    linux = PLATFORM_CLASSES["linux"]()

    assert linux.expected_tag == f"{PYTHON_TAG}-{ABI_TAG}-{linux.platform_tag}"
