import sys
import zipfile
from pathlib import Path, PurePosixPath

import pytest


SCRIPTS_DIR = Path(__file__).resolve().parents[2] / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

from audit_release_wheel import AuditError, audit_wheel_archive  # noqa: E402
from release_wheel_contract import get_platform_contract  # noqa: E402


PROJECT_VERSION = "0.0.4"
REQUIRED_COMPONENT_CASES = tuple(
    (platform, component)
    for platform in ("macos", "linux", "windows")
    for component in get_platform_contract(platform).required_native_components
)


def native_members(platform: str) -> dict[str, str]:
    if platform == "macos":
        return {
            "extension": "pypgo/_pypgo.cpython-312-darwin.so",
            "tbb": "pypgo/.dylibs/libtbb.12.dylib",
            "gmp": "pypgo/.dylibs/libgmp.10.dylib",
            "gmpxx": "pypgo/.dylibs/libgmpxx.4.dylib",
            "mpfr": "pypgo/.dylibs/libmpfr.6.dylib",
        }
    if platform == "linux":
        members = {
            "extension": "pypgo/_pypgo.cpython-312-x86_64-linux-gnu.so",
            "tbb": "pypgo.libs/libtbb-deadbeef.so.12",
            "gmp": "pypgo.libs/libgmp-a1b2c3d4.so.10",
            "gmpxx": "pypgo.libs/libgmpxx-b2c3d4e5.so.4",
            "mpfr": "pypgo.libs/libmpfr-c3d4e5f6.so.6",
            "mkl_core": "pypgo.libs/libmkl_core-deadbeef.so.2",
            "mkl_tbb_thread": "pypgo.libs/libmkl_tbb_thread-deadbeef.so.2",
        }
        for component in get_platform_contract(platform).required_native_components:
            if component.startswith("mkl_") and component not in members:
                members[component] = f"pypgo.libs/lib{component}.so.2"
        return members

    members = {
        "extension": "pypgo/_pypgo.cp312-win_amd64.pyd",
        "tbb": "pypgo.libs/tbb12.dll",
        "gmp": "pypgo.libs/gmp-10-a1b2c3d4.dll",
        "gmpxx": "pypgo.libs/gmpxx-4-b2c3d4e5.dll",
        "mpfr": "pypgo.libs/mpfr-6-c3d4e5f6.dll",
    }
    for component in get_platform_contract(platform).required_native_components:
        if component.startswith("mkl_"):
            members[component] = f"pypgo.libs/{component}.2.dll"
    return members


def write_synthetic_wheel(
    tmp_path: Path,
    platform: str,
    *,
    filename_tag: str | None = None,
    metadata_tags: tuple[str, ...] | None = None,
    root_is_purelib: bool = False,
    remove_component: str | None = None,
    mangle_original_component: str | None = None,
    extra_members: tuple[str, ...] = (),
) -> Path:
    contract = get_platform_contract(platform)
    filename_tag = filename_tag or contract.expected_tag
    metadata_tags = metadata_tags or (contract.expected_tag,)
    wheel = tmp_path / f"pypgo-{PROJECT_VERSION}-{filename_tag}.whl"
    members = native_members(platform)
    if remove_component is not None:
        members.pop(remove_component)
    if mangle_original_component is not None:
        original = PurePosixPath(members[mangle_original_component])
        members[mangle_original_component] = original.with_name(
            f"{original.stem}-deadbeef{original.suffix}"
        ).as_posix()

    initializer = (
        "# patched by delvewheel\n"
        "import os\n"
        "os.add_dll_directory('pypgo.libs')\n"
        if platform == "windows"
        else "from . import _pypgo\n"
    )
    wheel_metadata = (
        "Wheel-Version: 1.0\n"
        "Generator: release-wheel-test\n"
        f"Root-Is-Purelib: {'true' if root_is_purelib else 'false'}\n"
        + "".join(f"Tag: {tag}\n" for tag in metadata_tags)
    )
    with zipfile.ZipFile(wheel, "w") as archive:
        archive.writestr("pypgo/__init__.py", initializer)
        for member in members.values():
            archive.writestr(member, b"native-placeholder")
        for member in extra_members:
            archive.writestr(member, b"placeholder")
        archive.writestr(
            f"pypgo-{PROJECT_VERSION}.dist-info/WHEEL",
            wheel_metadata,
        )
    return wheel


def inspect(wheel: Path, platform: str):
    with zipfile.ZipFile(wheel) as archive:
        return audit_wheel_archive(
            wheel,
            archive,
            platform,
            PROJECT_VERSION,
        )


@pytest.mark.parametrize("platform", ["macos", "linux", "windows"])
def test_valid_repaired_wheel_satisfies_archive_contract(
    tmp_path: Path,
    platform: str,
) -> None:
    wheel = write_synthetic_wheel(tmp_path, platform)

    result = inspect(wheel, platform)

    assert result.filename_tags == result.metadata_tags
    assert result.extension_member.startswith("pypgo/_pypgo.")


def test_filename_and_metadata_tags_must_match(tmp_path: Path) -> None:
    wheel = write_synthetic_wheel(
        tmp_path,
        "linux",
        filename_tag="cp312-cp312-linux_x86_64",
    )

    with pytest.raises(AuditError, match="filename and WHEEL metadata tags differ"):
        inspect(wheel, "linux")


def test_linux_tag_must_include_contract_tag(tmp_path: Path) -> None:
    wrong_tag = "cp312-cp312-manylinux_2_27_x86_64"
    wheel = write_synthetic_wheel(
        tmp_path,
        "linux",
        filename_tag=wrong_tag,
        metadata_tags=(wrong_tag,),
    )

    with pytest.raises(AuditError, match="unexpected compatibility tags for linux"):
        inspect(wheel, "linux")


def test_linux_accepts_more_compatible_auditwheel_tag(tmp_path: Path) -> None:
    more_compatible_tag = "cp312-cp312-manylinux_2_27_x86_64"
    contract_tag = get_platform_contract("linux").expected_tag
    wheel = write_synthetic_wheel(
        tmp_path,
        "linux",
        filename_tag=(
            "cp312-cp312-"
            "manylinux_2_27_x86_64.manylinux_2_28_x86_64"
        ),
        metadata_tags=(more_compatible_tag, contract_tag),
    )

    result = inspect(wheel, "linux")

    assert result.filename_tags == result.metadata_tags


def test_wheel_must_not_be_purelib(tmp_path: Path) -> None:
    wheel = write_synthetic_wheel(
        tmp_path,
        "macos",
        root_is_purelib=True,
    )

    with pytest.raises(AuditError, match="Root-Is-Purelib: false"):
        inspect(wheel, "macos")


@pytest.mark.parametrize(("platform", "component"), REQUIRED_COMPONENT_CASES)
def test_required_native_component_must_be_bundled(
    tmp_path: Path,
    platform: str,
    component: str,
) -> None:
    wheel = write_synthetic_wheel(
        tmp_path,
        platform,
        remove_component=component,
    )

    with pytest.raises(
        AuditError,
        match=rf"missing bundled native components: .*{component}",
    ):
        inspect(wheel, platform)


@pytest.mark.parametrize("component", ["tbb", "mkl_core", "mkl_avx2"])
def test_windows_required_original_runtime_name_cannot_be_mangled(
    tmp_path: Path,
    component: str,
) -> None:
    wheel = write_synthetic_wheel(
        tmp_path,
        "windows",
        mangle_original_component=component,
    )

    with pytest.raises(AuditError, match="missing original oneMKL/oneTBB"):
        inspect(wheel, "windows")


def test_multiple_extensions_are_rejected(tmp_path: Path) -> None:
    wheel = write_synthetic_wheel(
        tmp_path,
        "linux",
        extra_members=("pypgo/_pypgo.second.so",),
    )

    with pytest.raises(AuditError, match="exactly one pypgo/_pypgo extension"):
        inspect(wheel, "linux")


def test_unintended_source_root_is_rejected(tmp_path: Path) -> None:
    wheel = write_synthetic_wheel(
        tmp_path,
        "macos",
        extra_members=("core/accidental.py",),
    )

    with pytest.raises(AuditError, match="unintended namespace package roots"):
        inspect(wheel, "macos")
