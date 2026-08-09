#!/usr/bin/env python3
"""Read the required uv version from pyproject.toml."""

from __future__ import annotations

import re
import tomllib
from pathlib import Path


def read_uv_version(source_dir: Path) -> str:
    with (source_dir / "pyproject.toml").open("rb") as stream:
        requirement = tomllib.load(stream)["tool"]["uv"]["required-version"]

    match = re.fullmatch(r"==([0-9]+(?:\.[0-9]+){2})", requirement)
    if not match:
        raise ValueError(
            "tool.uv.required-version must be an exact version such as ==0.11.6"
        )
    return match.group(1)


if __name__ == "__main__":
    print(read_uv_version(Path(__file__).resolve().parents[1]))
