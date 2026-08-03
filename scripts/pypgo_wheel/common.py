"""Shared subprocess and filesystem helpers for pypgo wheel packaging."""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path
from typing import Mapping, Sequence


class PackagingError(RuntimeError):
    """Raised when a wheel packaging step cannot satisfy its contract."""


def run(
    command: Sequence[str | os.PathLike[str]],
    *,
    cwd: Path | None = None,
    env: Mapping[str, str] | None = None,
) -> None:
    printable = tuple(str(item) for item in command)
    print("+", " ".join(printable), flush=True)
    subprocess.run(printable, cwd=cwd, env=env, check=True)


def run_capture(
    command: Sequence[str | os.PathLike[str]],
    *,
    cwd: Path | None = None,
    env: Mapping[str, str] | None = None,
) -> str:
    printable = tuple(str(item) for item in command)
    print("+", " ".join(printable), flush=True)
    result = subprocess.run(
        printable,
        cwd=cwd,
        env=env,
        check=True,
        capture_output=True,
        text=True,
    )
    output = result.stdout + result.stderr
    print(output, end="" if output.endswith("\n") else "\n")
    return output


def require_tool(name: str) -> str:
    executable = shutil.which(name)
    if executable is None:
        raise PackagingError(f"required executable is unavailable: {name}")
    return executable


def prepare_wheel_dir(wheel_dir: Path) -> Path:
    wheel_dir = wheel_dir.resolve()
    wheel_dir.mkdir(parents=True, exist_ok=True)
    existing = sorted(wheel_dir.glob("*.whl"))
    if existing:
        raise PackagingError(
            f"wheel output directory already contains wheel(s): {existing}"
        )
    return wheel_dir


def exactly_one_wheel(directory: Path) -> Path:
    wheels = sorted(directory.glob("*.whl"))
    if len(wheels) != 1:
        raise PackagingError(f"expected exactly one wheel in {directory}, found {wheels}")
    return wheels[0].resolve()


def publish_wheel(repaired_wheel: Path, wheel_dir: Path) -> Path:
    destination = wheel_dir / repaired_wheel.name
    if destination.exists():
        raise PackagingError(f"refusing to overwrite wheel: {destination}")
    shutil.copy2(repaired_wheel, destination)
    return destination.resolve()


def clean_setuptools_wheel_staging(source_dir: Path) -> None:
    """Remove stale setuptools install trees without touching CMake caches."""

    build_dir = (source_dir / "build").resolve()
    if not build_dir.is_dir():
        return
    for pattern in ("lib.*", "bdist.*"):
        for staging_dir in build_dir.glob(pattern):
            if staging_dir.is_dir():
                shutil.rmtree(staging_dir)


def build_raw_wheel(source_dir: Path, raw_dir: Path, env: Mapping[str, str]) -> Path:
    require_tool("uv")
    clean_setuptools_wheel_staging(source_dir)
    raw_dir.mkdir(parents=True, exist_ok=True)
    run(
        (
            "uv",
            "build",
            "--wheel",
            "--no-build-isolation",
            "--out-dir",
            str(raw_dir),
        ),
        cwd=source_dir,
        env=env,
    )
    return exactly_one_wheel(raw_dir)


def write_report(report_dir: Path | None, filename: str, contents: str) -> None:
    if report_dir is None:
        return
    report_dir = report_dir.resolve()
    report_dir.mkdir(parents=True, exist_ok=True)
    (report_dir / filename).write_text(contents, encoding="utf-8")
