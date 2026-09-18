#!/usr/bin/env python3
"""Read the package version from pyproject.toml."""

from __future__ import annotations

import tomllib
from pathlib import Path


def read_project_version(source_dir: Path) -> str:
    with (source_dir / "pyproject.toml").open("rb") as stream:
        return tomllib.load(stream)["project"]["version"]


if __name__ == "__main__":
    print(read_project_version(Path(__file__).resolve().parents[1]))
