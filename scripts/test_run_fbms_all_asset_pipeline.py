#!/usr/bin/env python3
"""Unit tests for the JSON-driven FBMS-all asset pipeline runner."""

from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = REPO_ROOT / "scripts" / "run_fbms_all_asset_pipeline.py"


def load_module():
    spec = importlib.util.spec_from_file_location("run_fbms_all_asset_pipeline", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class RunFbmsAllAssetPipelineTest(unittest.TestCase):
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
                    "tpms_shapes": ["tpms_schwarz_p"],
                    "parameters": {
                        "remesh": {"target_edge_length": 0.02},
                    },
                }
            ),
            encoding="utf-8",
        )
        return config

    def write_r384_config(self, source: Path, working: Path) -> Path:
        config = self.root / "pipeline_r384.json"
        config.write_text(
            json.dumps(
                {
                    "source_asset_directory": str(source),
                    "working_directory": str(working),
                    "build_dir": "build/base_no_mkl",
                    "tpms_shapes": ["tpms_schwarz_p"],
                    "parameters": {
                        "raw": {"resolution": 384},
                        "tetwild": {"lr": 0.008},
                    },
                }
            ),
            encoding="utf-8",
        )
        return config

    def test_load_config_applies_defaults(self) -> None:
        source = self.root / "source"
        working = self.root / "work"
        source.mkdir()
        config = self.module.load_config(self.write_config(source, working))

        self.assertEqual(config.source_asset_directory, source)
        self.assertEqual(config.working_directory, working)
        self.assertEqual(config.parameters["raw"]["resolution"], 512)
        self.assertEqual(config.parameters["remesh"]["target_edge_length"], 0.02)
        self.assertEqual(config.parameters["tetwild"]["lr"], 0.008)
        self.assertEqual(config.tpms_shapes, ("tpms_schwarz_p",))

    def test_discover_fbms_cases_supports_both_name_forms_and_excludes_spheres(self) -> None:
        source = self.root / "source"
        nested = source / "fbms"
        nested.mkdir(parents=True)
        (nested / "g0_b5.obj").write_text("v 0 0 0\n", encoding="utf-8")
        (nested / "g0b10.obj").write_text("v 0 0 0\n", encoding="utf-8")
        (nested / "g0_b5_bounding_sphere.obj").write_text("v 0 0 0\n", encoding="utf-8")
        (nested / "tpms_schwarz_p_fbms.obj").write_text("v 0 0 0\n", encoding="utf-8")

        cases = self.module.discover_fbms_cases(source)

        self.assertEqual([case.name for case in cases], ["g0_b5", "g0b10"])

    def test_prepare_assets_dry_run_builds_generation_commands_without_outputs(self) -> None:
        source = self.root / "source"
        source.mkdir()
        source_obj = source / "g0_b5.obj"
        source_obj.write_text("v 0 0 0\n", encoding="utf-8")
        working = self.root / "work"
        config = self.module.load_config(self.write_config(source, working))
        cases = self.module.discover_fbms_cases(source)

        actions = self.module.prepare_assets(config, cases, dry_run=True, overwrite=False)

        labels = [action.label for action in actions]
        self.assertIn("copy", labels)
        self.assertIn("prepare-sphere:g0_b5", labels)
        self.assertIn("prepare-tpms", labels)
        self.assertFalse((working / "asset").exists())

    def test_compute_remesh_edge_length_arg_uses_raw_mean_edge(self) -> None:
        scale = self.module.compute_remesh_edge_length_arg(
            {"edge_length": {"mean": 0.005}},
            0.02,
        )

        self.assertEqual(scale, 4.0)

    def test_write_tetmesher_config_matches_tetwild_shape(self) -> None:
        source = self.root / "source"
        source.mkdir()
        working = self.root / "work"
        config = self.module.load_config(self.write_config(source, working))
        target = self.module.fbms_target(config, "g0_b5")

        data = self.module.write_tetmesher_config(config, target, dry_run=False)

        self.assertEqual(data["backend"], "tetwild")
        self.assertEqual(data["tetwild"], {"lr": 0.008, "epsr": 0.001})
        self.assertEqual(data["material"]["young_modulus"], 10000000)
        self.assertEqual(data["output_mesh"], "g0_b5.veg")
        self.assertTrue(target.tet_config.exists())

    def test_target_directories_follow_config_resolution(self) -> None:
        source = self.root / "source"
        source.mkdir()
        working = self.root / "work"
        config = self.module.load_config(self.write_r384_config(source, working))
        fbms = self.module.fbms_target(config, "g0_b15")
        baseline = self.module.baseline_target(config, "g0_b15", "tpms_schwarz_p", 0.5)

        self.assertIn("raw_mesh_vdb_r384_t002", str(fbms.raw_obj))
        self.assertIn("remesh_cgal_iso_r384_t002_e002", str(fbms.remesh_obj))
        self.assertIn("tetmesh_tetwild_r384_t002_lr0008", str(fbms.tet_dir))
        self.assertIn("raw_mesh_vdb_r384_vmatch", str(baseline.raw_obj))
        self.assertIn("tetmesh_tetwild_r384_vmatch_lr0008", str(baseline.tet_dir))
        self.assertEqual(self.module.summary_path(config).name, "pipeline_r384_summary.json")

    def test_replace_readme_marker_block_only_updates_marker_content(self) -> None:
        original = "\n".join(
            [
                "# Title",
                "keep before",
                self.module.README_RESULTS_START,
                "old table",
                self.module.README_RESULTS_END,
                "keep after",
                "",
            ]
        )

        updated = self.module.replace_readme_marker_block(original, "new table")

        self.assertIn("keep before", updated)
        self.assertIn("new table", updated)
        self.assertNotIn("old table", updated)
        self.assertIn("keep after", updated)

    def test_progress_line_includes_stage_target_and_status(self) -> None:
        line = self.module.progress_line(
            "stage",
            stage="tet",
            target="baseline/g0_b5/tpms_schwarz_p",
            status="start",
        )

        self.assertEqual(
            line,
            "[fbms-all] stage stage=tet target=baseline/g0_b5/tpms_schwarz_p status=start",
        )

    def test_repair_paths_are_derived_from_veg_obj_without_overwriting_original(self) -> None:
        source = self.root / "source"
        source.mkdir()
        working = self.root / "work"
        config = self.module.load_config(self.write_config(source, working))
        target = self.module.baseline_target(config, "g0_b5", "tpms_schwarz_p", 0.5)

        self.assertEqual(target.repaired_veg_obj.name, "tpms_schwarz_p.veg.repaired.obj")
        self.assertEqual(target.repaired_veg_quality.name, "tpms_schwarz_p.veg.repaired.quality.json")
        self.assertEqual(target.repair_report.name, "tpms_schwarz_p.veg.repair_report.json")
        self.assertNotEqual(target.repaired_veg_obj, target.veg_obj)

    def test_validation_repairs_failed_boundary_and_dumps_repaired_components(self) -> None:
        source = self.root / "source"
        source.mkdir()
        working = self.root / "work"
        config = self.module.load_config(self.write_config(source, working))
        target = self.module.baseline_target(config, "g0_b5", "tpms_schwarz_p", 0.5)
        target.tet_dir.mkdir(parents=True)
        target.veg.write_text("veg\n", encoding="utf-8")
        target.veg_obj.write_text("v 0 0 0\nf 1 1 1\n", encoding="utf-8")

        commands = []

        def fake_run_command(command, dry_run):
            commands.append(command)
            if command.label.startswith("boundary-quality:"):
                target.veg_quality.parent.mkdir(parents=True, exist_ok=True)
                target.veg_quality.write_text(json.dumps({"passed": False, "enclosed_volume": 1.0}), encoding="utf-8")
            elif command.label.startswith("boundary-repair:"):
                target.repaired_veg_obj.write_text("v 0 0 0\nf 1 1 1\n", encoding="utf-8")
                target.repair_report.write_text(json.dumps({"final": {"connected_components": 3}}), encoding="utf-8")
            elif command.label.startswith("boundary-repaired-quality:"):
                target.repaired_veg_quality.write_text(
                    json.dumps({"passed": True, "enclosed_volume": 1.01}),
                    encoding="utf-8",
                )
            elif command.label.startswith("components:"):
                target.component_dir.mkdir(parents=True, exist_ok=True)
                (target.component_dir / "components.json").write_text(json.dumps([{}, {}, {}]), encoding="utf-8")
            return 0

        captures = []

        def fake_capture(argv, output_path, dry_run):
            captures.append((argv, output_path))
            output_path.parent.mkdir(parents=True, exist_ok=True)
            output_path.write_text("#vtx: 1\n#elements: 1\ntotal volume: 1\n", encoding="utf-8")

        self.module.run_command = fake_run_command
        self.module.run_command_capture = fake_capture

        self.module.run_validate_stage(config, target, dry_run=False, overwrite=True)

        labels = [command.label for command in commands]
        self.assertIn("boundary-quality:g0_b5:tpms_schwarz_p", labels)
        self.assertIn("boundary-repair:g0_b5:tpms_schwarz_p", labels)
        self.assertIn("boundary-repaired-quality:g0_b5:tpms_schwarz_p", labels)
        component_command = next(command for command in commands if command.label.startswith("components:"))
        self.assertIn(str(target.repaired_veg_obj), component_command.argv)

    def test_validation_rewrites_existing_components_after_repair_passes(self) -> None:
        source = self.root / "source"
        source.mkdir()
        working = self.root / "work"
        config = self.module.load_config(self.write_config(source, working))
        target = self.module.baseline_target(config, "g0_b5", "tpms_schwarz_p", 0.5)
        target.tet_dir.mkdir(parents=True)
        target.veg.write_text("veg\n", encoding="utf-8")
        target.veg_obj.write_text("v 0 0 0\nf 1 1 1\n", encoding="utf-8")
        target.component_dir.mkdir(parents=True, exist_ok=True)
        (target.component_dir / "components.json").write_text(json.dumps([{}, {}, {}]), encoding="utf-8")

        commands = []

        def fake_run_command(command, dry_run):
            commands.append(command)
            if command.label.startswith("boundary-quality:"):
                target.veg_quality.parent.mkdir(parents=True, exist_ok=True)
                target.veg_quality.write_text(json.dumps({"passed": False}), encoding="utf-8")
            elif command.label.startswith("boundary-repair:"):
                target.repaired_veg_obj.write_text("v 0 0 0\nf 1 1 1\n", encoding="utf-8")
                target.repair_report.write_text(json.dumps({}), encoding="utf-8")
            elif command.label.startswith("boundary-repaired-quality:"):
                target.repaired_veg_quality.write_text(json.dumps({"passed": True}), encoding="utf-8")
            elif command.label.startswith("components:"):
                (target.component_dir / "components.json").write_text(json.dumps([{"repaired": True}, {}, {}]), encoding="utf-8")
            return 0

        def fake_capture(argv, output_path, dry_run):
            output_path.parent.mkdir(parents=True, exist_ok=True)
            output_path.write_text("#vtx: 1\n#elements: 1\ntotal volume: 1\n", encoding="utf-8")

        self.module.run_command = fake_run_command
        self.module.run_command_capture = fake_capture

        self.module.run_validate_stage(config, target, dry_run=False, overwrite=False)

        component_command = next(command for command in commands if command.label.startswith("components:"))
        self.assertIn(str(target.repaired_veg_obj), component_command.argv)

    def test_summary_uses_repaired_quality_when_original_boundary_failed(self) -> None:
        source = self.root / "source"
        source.mkdir()
        working = self.root / "work"
        config = self.module.load_config(self.write_config(source, working))
        target = self.module.baseline_target(config, "g0_b5", "tpms_schwarz_p", 0.5)
        for path, data in [
            (target.raw_quality, {"passed": True, "enclosed_volume": 1.0, "invalid_triangles": 0}),
            (target.remesh_quality, {"passed": False, "self_intersections": 3, "edge_length": {"mean": 0.02}}),
            (target.veg_quality, {"passed": False, "enclosed_volume": 1.0}),
            (target.repaired_veg_quality, {"passed": True, "enclosed_volume": 1.02}),
        ]:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(data), encoding="utf-8")
        target.veg_info.parent.mkdir(parents=True, exist_ok=True)
        target.veg_info.write_text("#vtx: 10\n#elements: 20\ntotal volume: 1.5\n", encoding="utf-8")
        target.component_dir.mkdir(parents=True, exist_ok=True)
        (target.component_dir / "components.json").write_text(json.dumps([{}, {}, {}]), encoding="utf-8")
        target.repair_report.write_text(json.dumps({"input": str(target.veg_obj)}), encoding="utf-8")

        summary = self.module.collect_target_summary(target)

        self.assertTrue(summary["boundary_passed"])
        self.assertEqual(summary["boundary_source"], "repaired")
        self.assertTrue(summary["boundary_repair_attempted"])
        self.assertTrue(summary["boundary_repair_passed"])
        self.assertAlmostEqual(summary["boundary_repair_volume_drift"], 0.02)
        self.assertIn("pass(repaired)", self.module.markdown_results_table({"results": [summary]}))

    def test_summary_uses_repair_report_volume_drift_when_quality_volume_is_missing(self) -> None:
        source = self.root / "source"
        source.mkdir()
        working = self.root / "work"
        config = self.module.load_config(self.write_config(source, working))
        target = self.module.baseline_target(config, "g0_b5", "tpms_schwarz_p", 0.5)
        for path, data in [
            (target.raw_quality, {"passed": True, "enclosed_volume": 1.0, "invalid_triangles": 0}),
            (target.remesh_quality, {"passed": True, "edge_length": {"mean": 0.02}}),
            (target.veg_quality, {"passed": False, "enclosed_volume": None}),
            (target.repaired_veg_quality, {"passed": True, "enclosed_volume": None}),
        ]:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(data), encoding="utf-8")
        target.veg_info.parent.mkdir(parents=True, exist_ok=True)
        target.veg_info.write_text("#vtx: 10\n#elements: 20\ntotal volume: 1.5\n", encoding="utf-8")
        target.component_dir.mkdir(parents=True, exist_ok=True)
        (target.component_dir / "components.json").write_text(json.dumps([{}, {}, {}]), encoding="utf-8")
        target.repair_report.write_text(
            json.dumps(
                {
                    "before": {"component_signed_volumes": [2.0, -3.0]},
                    "final": {"component_signed_volumes": [2.0, 3.1]},
                }
            ),
            encoding="utf-8",
        )

        summary = self.module.collect_target_summary(target)

        self.assertAlmostEqual(summary["boundary_repair_volume_drift"], 0.02)


if __name__ == "__main__":
    unittest.main()
