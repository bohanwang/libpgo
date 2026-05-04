import argparse
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = REPO_ROOT / "examples" / "fbms" / "generate_shell_assets.py"


def load_module():
    spec = importlib.util.spec_from_file_location("generate_shell_assets", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class GenerateShellAssetsTest(unittest.TestCase):
    def setUp(self):
        self.module = load_module()

    def test_parse_job_accepts_volume_budget_and_truncating(self):
        job = self.module.parse_job(
            {
                "name": "volume_budget",
                "raw": {
                    "resolution": 64,
                    "fbms_thickness": -0.2,
                    "sphere_thickness": 0.02,
                    "enable_truncating": True,
                },
            },
            {
                "cases": ["g0_b3"],
                "raw": {"padding_ratio": 0.08},
                "remesh": {"edge_length": 1.0, "sharp_edge_angle": 180.0},
            },
        )

        self.assertEqual(job.raw.fbms_thickness, -0.2)
        self.assertTrue(job.raw.enable_truncating)

    def test_parse_job_rejects_zero_fbms_thickness(self):
        with self.assertRaisesRegex(ValueError, "fbms_thickness"):
            self.module.parse_job(
                {
                    "name": "zero",
                    "raw": {
                        "resolution": 64,
                        "fbms_thickness": 0.0,
                        "sphere_thickness": 0.02,
                    },
                },
                {
                    "cases": ["g0_b3"],
                    "raw": {"padding_ratio": 0.08},
                    "remesh": {"edge_length": 1.0, "sharp_edge_angle": 180.0},
                },
            )

    def test_generate_asset_passes_truncating_flag(self):
        case = self.module.AssetCase("unit", Path("case-dir"), "case")
        job = self.module.AssetJob(
            name="volume_budget",
            cases=("unit",),
            raw=self.module.RawParams(
                resolution=64,
                fbms_thickness=-0.2,
                sphere_thickness=0.02,
                padding_ratio=0.08,
                enable_truncating=True,
            ),
            remesh=self.module.RemeshParams(edge_length=1.0, sharp_edge_angle=180.0),
            output_dir=str(Path(tempfile.mkdtemp()) / "{job}" / "{case}"),
            raw_filename="raw.obj",
            remesh_filename="remesh.obj",
            stats_filename="stats.json",
        )
        args = argparse.Namespace(
            dry_run=True,
            raw_only=True,
            remesh_only=False,
            overwrite=True,
            skip_existing=False,
        )

        commands = []
        with mock.patch.object(self.module, "check_inputs"), mock.patch.object(
            self.module, "run_command", side_effect=lambda command, dry_run: commands.append(command)
        ):
            self.module.generate_asset(args, Path("build/base_no_mkl"), job, case)

        self.assertEqual(len(commands), 1)
        self.assertIn("--enable-truncating", commands[0])
        self.assertIn("-0.2", commands[0])


if __name__ == "__main__":
    unittest.main()
