#!/usr/bin/env python3
"""Unit tests for the generic runIPCSim batch runner."""

from __future__ import annotations

import importlib.util
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
        self.assertIsNone(case.vtu_config)
        self.assertEqual(jobs["squash_regression"].cases, ("tet_box_squash", "cubic_box_squash"))
        self.assertEqual(jobs["squash_regression"].stages, ("sim", "abc"))

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

    def test_load_fbms_config_uses_generated_case_assets(self) -> None:
        runner = load_runner_module()
        build_dir, cases, jobs = runner.load_config(REPO_ROOT / "examples" / "fbms" / "fbms_batch.json")

        case = cases["g0_b8_case1_pressure"]

        self.assertEqual(build_dir, REPO_ROOT / "build" / "base_no_mkl")
        self.assertEqual(
            case.sim_config,
            REPO_ROOT / "examples" / "fbms" / "generated" / "r128_default" / "g0_b8" / "g0_b8_case1_pressure-ipc.json",
        )
        self.assertEqual(
            case.anim_config,
            REPO_ROOT
            / "examples"
            / "fbms"
            / "generated"
            / "r128_default"
            / "g0_b8"
            / "g0_b8_case1_pressure-anim.json",
        )
        self.assertTrue(case.log)
        self.assertEqual(
            case.vtu_config,
            REPO_ROOT
            / "examples"
            / "fbms"
            / "generated"
            / "r128_default"
            / "g0_b8"
            / "g0_b8_case1_pressure-stress-vtu.json",
        )
        self.assertEqual(jobs["case1_pressure"].stages, ("sim", "abc", "vtu"))
        self.assertEqual(jobs["case1_pressure"].cases, ("g0_b8_case1_pressure", "g0_b3_case1_pressure"))
        self.assertEqual(
            jobs["g0_b8"].cases,
            (
                "g0_b8_case1_pressure",
                "g0_b8_case2_squash_floor_prototype",
                "g0_b8_case3_wall_impact_floor_prototype",
            ),
        )
        self.assertEqual(jobs["g0_b8"].stages, ("sim", "abc", "vtu"))
        self.assertEqual(jobs["g0_b8_post"].stages, ("abc", "vtu"))
        self.assertEqual(jobs["g0_b8_vtu"].stages, ("vtu",))
        self.assertEqual(
            jobs["g0_b3"].cases,
            (
                "g0_b3_case1_pressure",
                "g0_b3_case2_squash_floor_prototype",
                "g0_b3_case3_wall_impact_floor_prototype",
            ),
        )
        self.assertEqual(
            cases["g0_b3_case2_squash_floor_prototype"].anim_config,
            REPO_ROOT
            / "examples"
            / "fbms"
            / "generated"
            / "r128_default"
            / "g0_b3"
            / "g0_b3_case2_squash_floor-prototype-anim.json",
        )
        commands = runner.build_commands(build_dir, case, jobs["g0_b8"].stages, overwrite=True)
        self.assertEqual([command.label for command in commands], ["sim", "abc", "vtu"])
        self.assertEqual(
            commands[2].argv,
            [
                str(REPO_ROOT / "scripts" / "export_fbms_stress_vtu.py"),
                "--config",
                str(case.vtu_config),
                "--overwrite",
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
      "sim_config": "examples/ipc/shell/shell-ipc.json"
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
      "sim_config": "examples/ipc/shell/shell-ipc.json"
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

    def test_vtu_stage_requires_vtu_config(self) -> None:
        runner = load_runner_module()
        build_dir, cases, _ = runner.load_config(REPO_ROOT / "examples" / "ipc" / "ipc_batch.json")

        with self.assertRaisesRegex(ValueError, "vtu_config"):
            runner.build_commands(build_dir, cases["tet_box_squash"], ("vtu",), overwrite=False)


if __name__ == "__main__":
    unittest.main()
