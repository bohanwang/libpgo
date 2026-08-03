import sys
from dataclasses import FrozenInstanceError
from pathlib import Path

import pytest


SCRIPTS_DIR = Path(__file__).resolve().parents[2] / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

from release_wheel_contract import (  # noqa: E402
    ABI_TAG,
    EVIDENCE_CATEGORIES,
    MKL_DISPATCH_COMPONENTS,
    MKL_DISPATCH_LIBRARY_COMPONENTS,
    PLATFORM_CONTRACTS,
    PYTHON_TAG,
    RELEASE_DISTRIBUTION,
    SUPPORTED_PLATFORMS,
    WINDOWS_UNMANGLED_RUNTIME_DLLS,
    get_platform_contract,
)


def test_shared_release_identity() -> None:
    assert RELEASE_DISTRIBUTION == "pypgo"
    assert PYTHON_TAG == "cp312"
    assert ABI_TAG == "cp312"
    assert EVIDENCE_CATEGORIES == ("environment", "audit", "test")


@pytest.mark.parametrize(
    ("platform", "expected_tag", "repair_report"),
    [
        ("macos", "cp312-cp312-macosx_26_0_arm64", "delocate-listdeps.txt"),
        ("linux", "cp312-cp312-manylinux_2_28_x86_64", "auditwheel.txt"),
        ("windows", "cp312-cp312-win_amd64", "delvewheel.txt"),
    ],
)
def test_platform_contracts(
    platform: str,
    expected_tag: str,
    repair_report: str,
) -> None:
    contract = get_platform_contract(platform)

    assert contract.platform == platform
    assert contract.expected_tag == expected_tag
    assert contract.repair_report == repair_report
    assert contract.required_audit_evidence == (
        repair_report,
        "wheel-linkage.txt",
    )
    assert {"tbb", "gmp", "gmpxx", "mpfr"}.issubset(
        contract.required_native_components
    )


def test_supported_platforms_match_contract_keys() -> None:
    assert set(SUPPORTED_PLATFORMS) == set(PLATFORM_CONTRACTS)


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


def test_windows_original_runtime_names_are_part_of_its_contract() -> None:
    windows = get_platform_contract("windows")

    assert windows.required_original_runtime_names == WINDOWS_UNMANGLED_RUNTIME_DLLS
    assert not get_platform_contract("linux").required_original_runtime_names
    assert not get_platform_contract("macos").required_original_runtime_names


def test_platform_contract_is_immutable() -> None:
    with pytest.raises(FrozenInstanceError):
        get_platform_contract("linux").platform_tag = "linux_x86_64"  # type: ignore[misc]


def test_unknown_platform_is_rejected() -> None:
    with pytest.raises(ValueError, match="unsupported release-wheel platform"):
        get_platform_contract("freebsd")
