import json
import os
import re
import subprocess
import sys
from collections.abc import Mapping
from pathlib import Path

import pytest


SOURCE_DIR = Path(__file__).resolve().parents[2]
SCRIPTS_DIR = SOURCE_DIR / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

import pypgo_wheel.platform as platform_module  # noqa: E402
import pypgo_wheel_windows as windows_module  # noqa: E402
from pypgo_wheel.audit import AuditError  # noqa: E402
from pypgo_wheel.common import (  # noqa: E402
    PackagingError,
    clean_setuptools_wheel_staging,
    exactly_one_wheel,
    prepare_wheel_dir,
)
from pypgo_wheel.contracts.common import WheelPlatformContract  # noqa: E402
from pypgo_wheel.contracts.linux import LINUX_CONTRACT  # noqa: E402
from pypgo_wheel.contracts.macos import MACOS_CONTRACT  # noqa: E402
from pypgo_wheel.contracts.windows import (  # noqa: E402
    WINDOWS_CONTRACT,
    WINDOWS_FORCED_INCLUDE_DLLS,
    WINDOWS_UNMANGLED_RUNTIME_DLLS,
)
from pypgo_wheel_linux import LinuxPypgoWheelPlatform  # noqa: E402
from pypgo_wheel_macos import MacOSPypgoWheelPlatform  # noqa: E402
from pypgo_wheel.platform import PypgoWheelPlatform  # noqa: E402
from pypgo_wheel_windows import WindowsPypgoWheelPlatform  # noqa: E402


ENTRY_POINTS = (
    "pypgo_wheel_linux.py",
    "pypgo_wheel_macos.py",
    "pypgo_wheel_windows.py",
)

PLATFORM_IMPLEMENTATIONS = (
    (LinuxPypgoWheelPlatform, LINUX_CONTRACT),
    (MacOSPypgoWheelPlatform, MACOS_CONTRACT),
    (WindowsPypgoWheelPlatform, WINDOWS_CONTRACT),
)


def test_cmake_presets_separate_portability_and_runtime_layout() -> None:
    presets = json.loads((SOURCE_DIR / "CMakePresets.json").read_text(encoding="utf-8"))
    configure_presets = {
        preset["name"]: preset for preset in presets["configurePresets"]
    }

    base = configure_presets["base"]["cacheVariables"]
    wheel = configure_presets["pypgo-wheel"]["cacheVariables"]

    assert base["PGO_PORTABLE_BUILD"] == "OFF"
    assert base["PGO_ENABLE_RELEASE_DEBUG_INFO"] == "ON"
    assert base["PGO_RUNTIME_LAYOUT"] == "SOURCE"
    assert wheel["PGO_PORTABLE_BUILD"] == "ON"
    assert wheel["PGO_ENABLE_RELEASE_DEBUG_INFO"] == "OFF"
    assert wheel["PGO_RUNTIME_LAYOUT"] == "WHEEL"
    assert wheel["PGO_BUILD_TESTS"] == "OFF"
    assert wheel["PGO_BUILD_SUBPROJECTS"] == "OFF"


def test_wheel_runtime_layout_configures_without_source_staging(
    tmp_path: Path,
) -> None:
    result = subprocess.run(
        [
            "cmake",
            "-S",
            str(SOURCE_DIR / "tests" / "cmake" / "runtimeDependencies"),
            "-B",
            str(tmp_path / "build"),
            f"-DLIBPGO_SOURCE_DIR={SOURCE_DIR}",
            "-DPGO_RUNTIME_LAYOUT=WHEEL",
            "-DPGO_ENABLE_PYTHON=ON",
        ],
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    assert "Runtime dependency layout: WHEEL" in result.stdout


def test_active_environment_replaces_stale_cmake_python_cache(
    tmp_path: Path,
) -> None:
    result = subprocess.run(
        [
            "cmake",
            "-S",
            str(SOURCE_DIR / "tests" / "cmake" / "pythonDependencies"),
            "-B",
            str(tmp_path / "build"),
            f"-DLIBPGO_SOURCE_DIR={SOURCE_DIR}",
            f"-DEXPECTED_PYTHON_PREFIX={Path(sys.prefix).resolve()}",
        ],
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    assert "Using Python executable:" in result.stdout
    assert "Using Python dependency prefix:" in result.stdout


@pytest.mark.parametrize(
    ("release_debug_info", "expected_strip"),
    (("ON", "OFF"), ("OFF", "ON")),
)
def test_release_debug_info_controls_strip_policy(
    tmp_path: Path,
    release_debug_info: str,
    expected_strip: str,
) -> None:
    result = subprocess.run(
        [
            "cmake",
            "-S",
            str(SOURCE_DIR / "tests" / "cmake" / "releaseBuildPolicy"),
            "-B",
            str(tmp_path / "build"),
            "-DCMAKE_BUILD_TYPE=Release",
            f"-DLIBPGO_SOURCE_DIR={SOURCE_DIR}",
            f"-DPGO_ENABLE_RELEASE_DEBUG_INFO={release_debug_info}",
            f"-DEXPECTED_STRIP={expected_strip}",
        ],
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stdout + result.stderr


class FakeWheelPlatform(PypgoWheelPlatform):
    contract = LINUX_CONTRACT
    supported_system = sys.platform
    supported_machines = frozenset({"test"})

    def __init__(self, events: list[str], *, fail_audit: bool = False) -> None:
        self.events = events
        self.fail_audit = fail_audit

    def validate_host(self) -> None:
        self.events.append("validate")

    def validate_package_environment(self, source_dir: Path) -> None:
        del source_dir
        self.events.append("validate-package")

    def repair_wheel(
        self,
        raw_wheel: Path,
        repair_dir: Path,
        source_dir: Path,
        environment: Mapping[str, str],
    ) -> None:
        del raw_wheel, source_dir, environment
        self.events.append("repair")
        (repair_dir / "pypgo-repaired.whl").touch()

    def post_repair(self, repaired_wheel: Path, source_dir: Path) -> None:
        del repaired_wheel, source_dir
        self.events.append("post-repair")

    def collect_repair_report(self, repaired_wheel: Path, source_dir: Path) -> str:
        del repaired_wheel, source_dir
        self.events.append("repair-report")
        return "repair report\n"

    def native_component_pattern(self, component: str) -> re.Pattern[str]:
        return re.compile(re.escape(component))

    def audit_native_binaries(
        self,
        extension: Path,
        natives: tuple[Path, ...],
    ) -> str:
        del extension, natives
        return "unused"

    def audit(self, wheel: Path, source_dir: Path = SOURCE_DIR) -> str:
        del wheel, source_dir
        self.events.append("audit")
        if self.fail_audit:
            raise AuditError("synthetic audit failure")
        return "linkage report\n"


@pytest.mark.parametrize(("implementation", "contract"), PLATFORM_IMPLEMENTATIONS)
def test_platform_implementation_binds_its_contract(
    implementation: type[PypgoWheelPlatform],
    contract: WheelPlatformContract,
) -> None:
    platform = implementation()

    assert platform.contract is contract
    assert platform.platform == contract.platform


def test_prepare_wheel_dir_rejects_stale_wheels(tmp_path: Path) -> None:
    wheel_dir = tmp_path / "wheelhouse"
    wheel_dir.mkdir()
    (wheel_dir / "stale.whl").touch()

    with pytest.raises(PackagingError, match="already contains wheel"):
        prepare_wheel_dir(wheel_dir)


def test_exactly_one_wheel_requires_one_output(tmp_path: Path) -> None:
    with pytest.raises(PackagingError, match="expected exactly one wheel"):
        exactly_one_wheel(tmp_path)

    expected = tmp_path / "pypgo.whl"
    expected.touch()
    assert exactly_one_wheel(tmp_path) == expected.resolve()

    (tmp_path / "another.whl").touch()
    with pytest.raises(PackagingError, match="expected exactly one wheel"):
        exactly_one_wheel(tmp_path)


def test_clean_setuptools_wheel_staging_preserves_build_caches(
    tmp_path: Path,
) -> None:
    build_dir = tmp_path / "build"
    stale_lib = build_dir / "lib.macosx-26.0-arm64-cpython-312"
    stale_bdist = build_dir / "bdist.macosx-26.0-arm64"
    cmake_cache = build_dir / "base"
    extension_cache = build_dir / "temp.macosx-26.0-arm64-cpython-312"
    for directory in (stale_lib, stale_bdist, cmake_cache, extension_cache):
        directory.mkdir(parents=True)

    clean_setuptools_wheel_staging(tmp_path)

    assert not stale_lib.exists()
    assert not stale_bdist.exists()
    assert cmake_cache.is_dir()
    assert extension_cache.is_dir()


@pytest.mark.parametrize("entry_point", ENTRY_POINTS)
@pytest.mark.parametrize("arguments", (("--help",), ("package", "--help"), ("audit", "--help")))
def test_platform_entry_point_help_is_portable(
    entry_point: str,
    arguments: tuple[str, ...],
) -> None:
    result = subprocess.run(
        [sys.executable, str(SCRIPTS_DIR / entry_point), *arguments],
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    assert "usage:" in result.stdout


def test_base_template_repairs_audits_then_publishes(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    events: list[str] = []
    implementation = FakeWheelPlatform(events)

    def fake_build(_source: Path, raw_dir: Path, _env: Mapping[str, str]) -> Path:
        events.append("build")
        raw_dir.mkdir(parents=True)
        wheel = raw_dir / "pypgo-raw.whl"
        wheel.touch()
        return wheel

    monkeypatch.setattr(platform_module, "build_raw_wheel", fake_build)
    wheel_dir = tmp_path / "wheelhouse"
    report_dir = tmp_path / "reports"

    wheel = implementation.package(tmp_path, wheel_dir, report_dir)

    assert events == [
        "validate",
        "validate-package",
        "build",
        "repair",
        "post-repair",
        "repair-report",
        "audit",
    ]
    assert wheel == (wheel_dir / "pypgo-repaired.whl").resolve()
    assert wheel.is_file()
    assert (report_dir / "auditwheel.txt").read_text() == "repair report\n"
    assert (report_dir / "wheel-linkage.txt").read_text() == "linkage report\n"


def test_base_template_does_not_publish_failed_audit(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    implementation = FakeWheelPlatform([], fail_audit=True)

    def fake_build(_source: Path, raw_dir: Path, _env: Mapping[str, str]) -> Path:
        raw_dir.mkdir(parents=True)
        wheel = raw_dir / "pypgo-raw.whl"
        wheel.touch()
        return wheel

    monkeypatch.setattr(platform_module, "build_raw_wheel", fake_build)
    wheel_dir = tmp_path / "wheelhouse"

    with pytest.raises(AuditError, match="synthetic audit failure"):
        implementation.package(tmp_path, wheel_dir)

    assert not list(wheel_dir.glob("*.whl"))


def test_windows_repair_uses_contract_dll_lists(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    implementation = WindowsPypgoWheelPlatform()
    implementation.native_directories = (tmp_path / "native",)
    commands: list[tuple[str, ...]] = []

    def fake_run(command, **_kwargs) -> str:
        normalized = tuple(str(item) for item in command)
        commands.append(normalized)
        return "repair output\n"

    monkeypatch.setattr(windows_module, "run_capture", fake_run)
    raw_wheel = tmp_path / "raw.whl"
    raw_wheel.touch()
    repair_dir = tmp_path / "repaired"
    repair_dir.mkdir()

    implementation.repair_wheel(
        raw_wheel,
        repair_dir,
        tmp_path,
        os.environ,
    )
    repair_command = commands[0]

    assert repair_command[repair_command.index("--include") + 1] == os.pathsep.join(
        WINDOWS_FORCED_INCLUDE_DLLS
    )
    assert repair_command[
        repair_command.index("--no-mangle") + 1
    ] == os.pathsep.join(WINDOWS_UNMANGLED_RUNTIME_DLLS)
