#!/usr/bin/env python3
"""Unit tests for the FBMS stress VTU exporter."""

from __future__ import annotations

import importlib.util
import json
import struct
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
EXPORTER_PATH = REPO_ROOT / "scripts" / "export_fbms_stress_vtu.py"


def load_exporter_module():
    spec = importlib.util.spec_from_file_location("export_fbms_stress_vtu", EXPORTER_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {EXPORTER_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def write_state(path: Path, displacement: list[float]) -> None:
    row_count = len(displacement)
    velocity = [0.0] * row_count
    acceleration = [0.0] * row_count
    payload = displacement + velocity + acceleration
    path.write_bytes(struct.pack("<iii", row_count, 3, 8) + struct.pack(f"<{len(payload)}d", *payload))


class ExportFbmsStressVtuTest(unittest.TestCase):
    def test_load_config_resolves_paths_relative_to_config(self) -> None:
        exporter = load_exporter_module()

        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            output_dir = root / "case_output"
            output_dir.mkdir()
            config_path = output_dir / "stress_vtu.json"
            config_path.write_text(
                json.dumps(
                    {
                        "veg": "../mesh.veg",
                        "states": "states",
                        "stress": "stress",
                        "output": "vtu",
                        "frame-end": 3,
                    }
                ),
                encoding="utf-8",
            )

            config = exporter.load_config(config_path)

            self.assertEqual(config.veg, root / "mesh.veg")
            self.assertEqual(config.states, output_dir / "states")
            self.assertEqual(config.stress, output_dir / "stress")
            self.assertEqual(config.output, output_dir / "vtu")
            self.assertEqual(config.frame_start, 0)
            self.assertEqual(config.frame_end, 3)

    def test_exports_tiny_vtu_series_from_json_config(self) -> None:
        exporter = load_exporter_module()

        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            case_dir = root / "case_output"
            states_dir = case_dir / "states"
            stress_dir = case_dir / "stress"
            states_dir.mkdir(parents=True)
            stress_dir.mkdir()

            (root / "mesh.veg").write_text(
                """# Vega mesh file.
# 4 vertices, 1 elements

*VERTICES
4 3 0 0
1 0 0 0
2 1 0 0
3 0 1 0
4 0 0 1

*ELEMENTS
TET
1 4 0
1 1 2 3 4
""",
                encoding="utf-8",
            )
            write_state(
                states_dir / "deform0010.u",
                [
                    0.0,
                    0.0,
                    0.0,
                    0.25,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    -0.5,
                ],
            )
            (stress_dir / "von_mises0010.json").write_text(
                json.dumps(
                    {
                        "frame": 10,
                        "time": 0.005,
                        "stress_type": "von_mises",
                        "location": "tet_element",
                        "values": [123.5],
                    }
                ),
                encoding="utf-8",
            )
            config_path = case_dir / "stress_vtu.json"
            config_path.write_text(
                json.dumps(
                    {
                        "veg": "../mesh.veg",
                        "states": "states",
                        "stress": "stress",
                        "output": "vtu",
                        "frame-start": 10,
                        "frame-end": 11,
                    }
                ),
                encoding="utf-8",
            )

            written = exporter.export_from_config(config_path, overwrite=False)

            frame_path = case_dir / "vtu" / "frame0010.vtu"
            pvd_path = case_dir / "vtu" / "series.pvd"
            self.assertEqual(written, [frame_path, pvd_path])
            self.assertTrue(frame_path.exists())
            self.assertTrue(pvd_path.exists())

            vtu_root = ET.parse(frame_path).getroot()
            self.assertEqual(vtu_root.tag, "VTKFile")
            piece = vtu_root.find("./UnstructuredGrid/Piece")
            self.assertIsNotNone(piece)
            self.assertEqual(piece.attrib["NumberOfPoints"], "4")
            self.assertEqual(piece.attrib["NumberOfCells"], "1")
            points = piece.find("./Points/DataArray")
            self.assertIsNotNone(points)
            self.assertIn("1.25 0 0", points.text)
            self.assertIn("0 0 0.5", points.text)
            stress = piece.find("./CellData/DataArray[@Name='von_mises']")
            self.assertIsNotNone(stress)
            self.assertEqual(stress.text.strip(), "123.5")
            log_stress = piece.find("./CellData/DataArray[@Name='von_mises_log10']")
            self.assertIsNotNone(log_stress)
            self.assertEqual(log_stress.text.strip(), "2.0916669575956846")
            clamped_stress = piece.find("./CellData/DataArray[@Name='von_mises_clamped_99']")
            self.assertIsNotNone(clamped_stress)
            self.assertEqual(clamped_stress.text.strip(), "123.5")

            pvd_root = ET.parse(pvd_path).getroot()
            dataset = pvd_root.find("./Collection/DataSet")
            self.assertIsNotNone(dataset)
            self.assertEqual(dataset.attrib["timestep"], "0.005")
            self.assertEqual(dataset.attrib["file"], "frame0010.vtu")

    def test_rejects_stress_value_count_mismatch(self) -> None:
        exporter = load_exporter_module()

        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            stress_path = root / "von_mises0000.json"
            stress_path.write_text(
                json.dumps(
                    {
                        "location": "tet_element",
                        "stress_type": "von_mises",
                        "values": [],
                    }
                ),
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "1 stress values"):
                exporter.read_stress(stress_path, tet_count=1)


if __name__ == "__main__":
    unittest.main()
