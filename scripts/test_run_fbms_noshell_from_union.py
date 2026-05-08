#!/usr/bin/env python3
"""Unit tests for the FBMS no-shell asset generator."""

from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = REPO_ROOT / "scripts" / "run_fbms_noshell_from_union.py"


def load_module():
    spec = importlib.util.spec_from_file_location("run_fbms_noshell_from_union", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class RunFbmsNoShellFromUnionTest(unittest.TestCase):
    def setUp(self) -> None:
        self.module = load_module()
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)

    def tearDown(self) -> None:
        self.tmp.cleanup()

    def write_config(self, source: Path, working: Path) -> Path:
        config = self.root / "pipeline.json"
        config.write_text(
            json.dumps(
                {
                    "source_asset_directory": str(source),
                    "working_directory": str(working),
                    "build_dir": "build/base_no_mkl",
                    "tpms_shapes": ["tpms_gyroid"],
                    "parameters": {
                        "raw": {"resolution": 512, "fbms_thickness": 0.02, "sphere_thickness": 0.03},
                        "remesh": {"target_edge_length": 0.02},
                    },
                }
            ),
            encoding="utf-8",
        )
        return config

    def test_parse_selected_thickness_prefers_volume_search_summary(self) -> None:
        thickness = self.module.parse_selected_thickness(
            "$ generateFBMSUnionSurface --fbms-thickness -0.6\n"
            "[volume-search] selected thickness = 0.0333906 volume = 0.606279 iterations = 11\n"
        )

        self.assertAlmostEqual(thickness, 0.0333906)

    def test_parse_selected_thickness_falls_back_to_command_line_for_fixed_fbms(self) -> None:
        thickness = self.module.parse_selected_thickness(
            "$ generateFBMSUnionSurface openvdb --fbms-thickness 0.02 --sphere-thickness 0.03\n"
        )

        self.assertEqual(thickness, 0.02)

    def test_fbms_target_paths_use_noshell_directories_and_config_thickness(self) -> None:
        source = self.root / "source"
        working = self.root / "work"
        source.mkdir()
        config = self.module.load_config(self.write_config(source, working))

        target = self.module.noshell_fbms_target(config, "g0_b15")

        self.assertEqual(target.thickness, 0.02)
        self.assertEqual(target.expected_components, 1)
        self.assertEqual(
            target.raw_obj,
            working / "raw_mesh_vdb_r512_t002_noshell" / "fbms" / "g0_b15_union_minus_sphere.obj",
        )
        self.assertEqual(
            target.remesh_obj,
            working
            / "remesh_cgal_iso_r512_t002_noshell_e002"
            / "fbms"
            / "g0_b15_union_minus_sphere_remesh.obj",
        )
        self.assertEqual(target.union_raw_quality.name, "g0_b15_union.quality.json")

    def test_baseline_target_parses_selected_thickness_from_union_log(self) -> None:
        source = self.root / "source"
        working = self.root / "work"
        source.mkdir()
        config = self.module.load_config(self.write_config(source, working))
        log = (
            working
            / "raw_mesh_vdb_r512_vmatch"
            / "baseline"
            / "g0_b15"
            / "logs"
            / "tpms_gyroid_generate.log"
        )
        log.parent.mkdir(parents=True)
        log.write_text("[volume-search] selected thickness = 0.0333906 volume = 0.606279\n", encoding="utf-8")

        target = self.module.noshell_baseline_target(config, "g0_b15", "tpms_gyroid")

        self.assertAlmostEqual(target.thickness, 0.0333906)
        self.assertIsNone(target.expected_components)
        self.assertEqual(
            target.raw_obj,
            working
            / "raw_mesh_vdb_r512_vmatch_noshell"
            / "baseline"
            / "g0_b15"
            / "tpms_gyroid_union_minus_sphere.obj",
        )

    def test_raw_command_uses_union_minus_sphere_surface_mode(self) -> None:
        source = self.root / "source"
        working = self.root / "work"
        source.mkdir()
        config = self.module.load_config(self.write_config(source, working))
        target = self.module.noshell_fbms_target(config, "g0_b15")

        command = self.module.raw_command(config, target)

        self.assertEqual(command.label, "raw-noshell:fbms:g0_b15:g0_b15")
        self.assertIn("--surface-mode", command.argv)
        self.assertIn("union-minus-sphere", command.argv)
        self.assertIn("--fbms-thickness", command.argv)
        self.assertIn("0.02", command.argv)
        self.assertEqual(command.argv[command.argv.index("--output-surface") + 1], str(target.raw_obj))

    def test_run_noshell_target_writes_raw_quality_before_remesh_scale(self) -> None:
        source = self.root / "source"
        working = self.root / "work"
        source.mkdir()
        config = self.module.load_config(self.write_config(source, working))
        target = self.module.noshell_fbms_target(config, "g0_b15")
        target.source_obj.parent.mkdir(parents=True)
        target.source_obj.write_text("v 0 0 0\n", encoding="utf-8")
        target.sphere_obj.write_text("v 0 0 0\n", encoding="utf-8")
        target.union_raw_log.parent.mkdir(parents=True)
        target.union_raw_log.write_text("$ generate --fbms-thickness 0.02\n", encoding="utf-8")
        target.union_raw_quality.parent.mkdir(parents=True)
        target.union_raw_quality.write_text(json.dumps({"passed": True}), encoding="utf-8")

        commands = []

        def fake_run_command(command, dry_run):
            commands.append(command)
            if command.label.startswith("raw-noshell:"):
                target.raw_obj.parent.mkdir(parents=True, exist_ok=True)
                target.raw_obj.write_text("v 0 0 0\n", encoding="utf-8")
            elif command.label.startswith("raw-quality-noshell:"):
                target.raw_quality.parent.mkdir(parents=True, exist_ok=True)
                target.raw_quality.write_text(
                    json.dumps({"passed": True, "edge_length": {"mean": 0.005}}),
                    encoding="utf-8",
                )
            elif command.label.startswith("remesh-noshell:"):
                target.remesh_obj.parent.mkdir(parents=True, exist_ok=True)
                target.remesh_obj.write_text("v 0 0 0\n", encoding="utf-8")
            elif command.label.startswith("remesh-quality-noshell:"):
                target.remesh_quality.parent.mkdir(parents=True, exist_ok=True)
                target.remesh_quality.write_text(json.dumps({"passed": True}), encoding="utf-8")
            return 0

        self.module.run_command = fake_run_command

        self.module.run_noshell_target(config, target, dry_run=False, overwrite=True)

        labels = [command.label for command in commands]
        self.assertEqual(
            labels,
            [
                "raw-noshell:fbms:g0_b15:g0_b15",
                "raw-quality-noshell:fbms:g0_b15:g0_b15",
                "remesh-noshell:fbms:g0_b15:g0_b15",
                "remesh-quality-noshell:fbms:g0_b15:g0_b15",
            ],
        )
        raw_quality = commands[1]
        self.assertEqual(raw_quality.argv[raw_quality.argv.index("--expected-components") + 1], "1")
        remesh = commands[2]
        self.assertEqual(remesh.argv[remesh.argv.index("--edge-length") + 1], "4")

    def test_baseline_quality_command_does_not_force_component_count(self) -> None:
        source = self.root / "source"
        working = self.root / "work"
        source.mkdir()
        config = self.module.load_config(self.write_config(source, working))

        command = self.module.quality_command(
            config,
            "raw-quality-noshell:baseline:g0_b15:tpms_gyroid",
            working / "in.obj",
            working / "quality.json",
            working / "quality.log",
            "raw",
            "warn",
            allow_failure=False,
            expected_components=None,
        )

        self.assertNotIn("--expected-components", command.argv)


if __name__ == "__main__":
    unittest.main()
