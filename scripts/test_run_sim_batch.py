#!/usr/bin/env python3
"""Unit tests for the generic runIPCSim batch runner."""

from __future__ import annotations

import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
RUNNER_PATH = REPO_ROOT / "scripts" / "run_sim_batch.py"


def load_runner_module():
    spec = importlib.util.spec_from_file_location("run_sim_batch", RUNNER_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {RUNNER_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class RunSimBatchRunnerTest(unittest.TestCase):
    def test_load_ipc_config_applies_defaults_and_builds_commands(self) -> None:
        runner = load_runner_module()
        build_dir, cases, jobs = runner.load_config(REPO_ROOT / "examples" / "ipc" / "ipc_batch.json")

        case = cases["cubic_box_squash"]

        self.assertEqual(build_dir, REPO_ROOT / "build" / "base_no_mkl")
        self.assertEqual(
            case.sim_config,
            REPO_ROOT / "examples" / "ipc" / "cubic" / "box-squash" / "box-ipc.json",
        )
        self.assertEqual(
            case.anim_config,
            REPO_ROOT / "examples" / "ipc" / "cubic" / "box-squash" / "anim.json",
        )
        self.assertTrue(case.log)
        self.assertEqual(jobs["squash_regression"].cases, ("tet_box_squash", "cubic_box_squash"))
        self.assertEqual(jobs["squash_regression"].stages, ("sim", "abc"))
        self.assertEqual(
            cases["shell_hang"].sim_config,
            REPO_ROOT / "examples" / "ipc" / "shell" / "shell-hang" / "shell-ipc.json",
        )
        self.assertEqual(
            cases["shell_drop"].sim_config,
            REPO_ROOT / "examples" / "ipc" / "shell" / "shell-drop" / "shell-ipc.json",
        )
        self.assertEqual(jobs["shell"].cases, ("shell_hang", "shell_drop"))

        commands = runner.build_commands(build_dir, case, jobs["squash_regression"].stages, overwrite=False)
        self.assertEqual(commands[0].label, "sim")
        self.assertEqual(
            commands[0].argv,
            [
                str(REPO_ROOT / "build" / "base_no_mkl" / "bin" / "runIPCSim"),
                str(case.sim_config),
                "--log",
            ],
        )
        self.assertEqual(commands[1].label, "abc")
        self.assertEqual(
            commands[1].argv,
            [
                str(REPO_ROOT / "build" / "base_no_mkl" / "bin" / "convertAnimation"),
                str(case.anim_config),
            ],
        )

    def test_sim_output_dir_comes_from_sim_config(self) -> None:
        runner = load_runner_module()
        _, cases, _ = runner.load_config(REPO_ROOT / "examples" / "ipc" / "ipc_batch.json")

        self.assertEqual(
            runner.sim_output_dir(cases["tet_box_squash"]),
            REPO_ROOT / "examples" / "ipc" / "tet" / "box-squash" / "ret-box-squash-ipc",
        )

    def test_load_config_rejects_unknown_job_case(self) -> None:
        runner = load_runner_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            config = Path(tmpdir) / "bad.json"
            config.write_text(
                """
{
  "cases": {
    "known": {
      "sim_config": "examples/ipc/shell/shell-hang/shell-ipc.json"
    }
  },
  "jobs": [
    {
      "name": "broken",
      "cases": ["missing"]
    }
  ]
}
""",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "unknown cases"):
                runner.load_config(config)

    def test_load_config_rejects_unknown_job_stage(self) -> None:
        runner = load_runner_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            config = Path(tmpdir) / "bad-stage.json"
            config.write_text(
                """
{
  "cases": {
    "known": {
      "sim_config": "examples/ipc/shell/shell-hang/shell-ipc.json"
    }
  },
  "jobs": [
    {
      "name": "broken",
      "stages": ["sim", "movie"],
      "cases": ["known"]
    }
  ]
}
""",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "unknown stages"):
                runner.load_config(config)

    def test_render_stage_requires_render_config(self) -> None:
        runner = load_runner_module()
        build_dir, cases, _ = runner.load_config(REPO_ROOT / "examples" / "ipc" / "ipc_batch.json")

        with self.assertRaisesRegex(ValueError, "render_config"):
            runner.build_commands(build_dir, cases["tet_box_squash"], ("render",), overwrite=False)

    def test_cli_runs_single_case_by_name(self) -> None:
        result = subprocess.run(
            [
                sys.executable,
                str(RUNNER_PATH),
                "--config",
                str(REPO_ROOT / "examples" / "ipc" / "ipc_batch.json"),
                "--case",
                "cubic_box_with_sphere_lite",
                "--dry-run",
            ],
            cwd=REPO_ROOT,
            check=False,
            capture_output=True,
            text=True,
        )

        self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
        self.assertIn("== case cubic_box_with_sphere_lite ==", result.stdout)
        self.assertIn("box-with-sphere-lite/box-ipc.json --log", result.stdout)
        self.assertIn("box-with-sphere-lite/anim.json", result.stdout)

    def test_cli_runs_shell_drop_case_by_name(self) -> None:
        result = subprocess.run(
            [
                sys.executable,
                str(RUNNER_PATH),
                "--config",
                str(REPO_ROOT / "examples" / "ipc" / "ipc_batch.json"),
                "--case",
                "shell_drop",
                "--dry-run",
            ],
            cwd=REPO_ROOT,
            check=False,
            capture_output=True,
            text=True,
        )

        self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
        self.assertIn("== case shell_drop ==", result.stdout)
        self.assertIn("shell/shell-drop/shell-ipc.json --log", result.stdout)
        self.assertIn("shell/shell-drop/anim.json", result.stdout)


if __name__ == "__main__":
    unittest.main()
