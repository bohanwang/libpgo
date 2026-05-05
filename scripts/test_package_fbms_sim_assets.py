#!/usr/bin/env python3
"""Unit tests for the FBMS simulation asset packaging script."""

from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = REPO_ROOT / "scripts" / "package_fbms_sim_assets.py"


def load_module():
    spec = importlib.util.spec_from_file_location("package_fbms_sim_assets", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class PackageFbmsSimAssetsTest(unittest.TestCase):
    def setUp(self) -> None:
        self.module = load_module()
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.source = self.root / "source"
        self.working = self.root / "work"
        self.source.mkdir()
        self.working.mkdir()
        (self.source / "g0_b5.obj").write_text("v 0 0 0\n", encoding="utf-8")
        self.config = self.root / "pipeline_r512.json"
        self.config.write_text(
            json.dumps(
                {
                    "source_asset_directory": str(self.source),
                    "working_directory": str(self.working),
                    "build_dir": "build/base_no_mkl",
                    "tpms_shapes": [
                        "tpms_schwarz_p",
                        "tpms_schwarz_d",
                        "tpms_gyroid",
                        "tpms_iwp",
                        "tpms_neovius",
                    ],
                    "parameters": {"raw": {"resolution": 512}, "tetwild": {"lr": 0.008}},
                }
            ),
            encoding="utf-8",
        )

    def tearDown(self) -> None:
        self.tmp.cleanup()

    def asset_dir(self, kind: str, case: str, name: str) -> Path:
        if kind == "fbms":
            return self.working / "tetmesh_tetwild_r512_t002_lr0008" / "fbms" / case
        return self.working / "tetmesh_tetwild_r512_vmatch_lr0008" / "baseline" / case / name

    def write_asset_files(self, kind: str, case: str, name: str, repaired: bool = False) -> None:
        asset_dir = self.asset_dir(kind, case, name)
        asset_dir.mkdir(parents=True, exist_ok=True)
        (asset_dir / f"{name}.veg").write_text(f"{kind}:{case}:{name}:veg\n", encoding="utf-8")
        (asset_dir / f"{name}.veg.obj").write_text(f"{kind}:{case}:{name}:original\n", encoding="utf-8")
        if repaired:
            (asset_dir / f"{name}.veg.repaired.obj").write_text(
                f"{kind}:{case}:{name}:repaired\n",
                encoding="utf-8",
            )

    def write_summary(self, results: list[dict]) -> Path:
        summary = self.working / "pipeline_r512_summary.json"
        summary.write_text(
            json.dumps(
                {
                    "generated_at_unix": 1,
                    "working_directory": str(self.working),
                    "results": results,
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        return summary

    def result(self, kind: str, case: str, name: str, repaired: bool = False) -> dict:
        return {
            "kind": kind,
            "case": case,
            "name": name,
            "boundary_passed": True,
            "boundary_source": "repaired" if repaired else "original",
            "boundary_repair_attempted": repaired,
            "boundary_repair_volume_drift": 1e-6 if repaired else None,
            "component_count": 3,
        }

    def test_packages_original_and_repaired_surfaces_with_uniform_names_and_manifest(self) -> None:
        shapes = ["tpms_schwarz_p", "tpms_schwarz_d", "tpms_gyroid", "tpms_iwp", "tpms_neovius"]
        results = [self.result("fbms", "g0_b5", "g0_b5")]
        self.write_asset_files("fbms", "g0_b5", "g0_b5")
        for shape in shapes:
            repaired = shape == "tpms_schwarz_p"
            results.append(self.result("baseline", "g0_b5", shape, repaired=repaired))
            self.write_asset_files("baseline", "g0_b5", shape, repaired=repaired)
        self.write_summary(results)

        manifest = self.module.package_assets_from_args(
            [
                "--config",
                str(self.config),
                "--overwrite",
            ]
        )

        out = self.working / "simulation_package"
        self.assertEqual((out / "g0_b5" / "fbms" / "g0_b5.veg").read_text(encoding="utf-8"), "fbms:g0_b5:g0_b5:veg\n")
        self.assertEqual(
            (out / "g0_b5" / "baseline" / "tpms_schwarz_p" / "tpms_schwarz_p.veg.obj").read_text(
                encoding="utf-8"
            ),
            "baseline:g0_b5:tpms_schwarz_p:repaired\n",
        )
        self.assertEqual(
            (out / "g0_b5" / "baseline" / "tpms_schwarz_d" / "tpms_schwarz_d.veg.obj").read_text(
                encoding="utf-8"
            ),
            "baseline:g0_b5:tpms_schwarz_d:original\n",
        )
        self.assertEqual(len(manifest["assets"]), 6)
        manifest_path = out / "package_manifest.json"
        self.assertTrue(manifest_path.exists())
        written_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        repaired_entries = [asset for asset in written_manifest["assets"] if asset["boundary_source"] == "repaired"]
        self.assertEqual(len(repaired_entries), 1)
        self.assertTrue(repaired_entries[0]["source_surface"].endswith(".veg.repaired.obj"))
        self.assertTrue(repaired_entries[0]["package_surface"].endswith("tpms_schwarz_p.veg.obj"))

    def test_rejects_failed_boundary_or_wrong_component_count(self) -> None:
        self.write_asset_files("fbms", "g0_b5", "g0_b5")
        self.write_summary(
            [
                {
                    **self.result("fbms", "g0_b5", "g0_b5"),
                    "boundary_passed": False,
                }
            ]
        )

        with self.assertRaisesRegex(RuntimeError, "not packageable"):
            self.module.package_assets_from_args(["--config", str(self.config), "--case", "g0_b5"])

    def test_case_and_shape_filters_package_subset(self) -> None:
        self.write_asset_files("fbms", "g0_b5", "g0_b5")
        self.write_asset_files("baseline", "g0_b5", "tpms_schwarz_p")
        self.write_asset_files("baseline", "g0_b5", "tpms_schwarz_d")
        self.write_summary(
            [
                self.result("fbms", "g0_b5", "g0_b5"),
                self.result("baseline", "g0_b5", "tpms_schwarz_p"),
                self.result("baseline", "g0_b5", "tpms_schwarz_d"),
            ]
        )

        manifest = self.module.package_assets_from_args(
            ["--config", str(self.config), "--shape", "tpms_schwarz_p", "--overwrite"]
        )

        self.assertEqual([asset["name"] for asset in manifest["assets"]], ["g0_b5", "tpms_schwarz_p"])
        self.assertFalse(
            (self.working / "simulation_package" / "g0_b5" / "baseline" / "tpms_schwarz_d").exists()
        )

    def test_dry_run_does_not_create_package(self) -> None:
        self.write_asset_files("fbms", "g0_b5", "g0_b5")
        self.write_summary([self.result("fbms", "g0_b5", "g0_b5")])

        manifest = self.module.package_assets_from_args(["--config", str(self.config), "--dry-run"])

        self.assertEqual(len(manifest["assets"]), 1)
        self.assertFalse((self.working / "simulation_package").exists())

    def test_refuses_to_overwrite_existing_files_without_flag(self) -> None:
        self.write_asset_files("fbms", "g0_b5", "g0_b5")
        self.write_summary([self.result("fbms", "g0_b5", "g0_b5")])
        out_file = self.working / "simulation_package" / "g0_b5" / "fbms" / "g0_b5.veg"
        out_file.parent.mkdir(parents=True)
        out_file.write_text("existing\n", encoding="utf-8")

        with self.assertRaisesRegex(FileExistsError, "exists"):
            self.module.package_assets_from_args(["--config", str(self.config)])


if __name__ == "__main__":
    unittest.main()
