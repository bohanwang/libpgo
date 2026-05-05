#!/usr/bin/env python3
"""Unit tests for the TPMS boundary repair helper."""

from __future__ import annotations

import subprocess
import sys
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = REPO_ROOT / "scripts" / "repair_tpms_mesh.py"


class RepairTpmsMeshTest(unittest.TestCase):
    def test_cli_help_does_not_expose_orient_positive(self) -> None:
        result = subprocess.run(
            [sys.executable, str(SCRIPT_PATH), "--help"],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

        self.assertNotIn("--orient-positive", result.stdout)


if __name__ == "__main__":
    unittest.main()
