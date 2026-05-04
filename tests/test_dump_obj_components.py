import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = REPO_ROOT / "scripts" / "dump_obj_components.py"


def load_module():
    spec = importlib.util.spec_from_file_location("dump_obj_components", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class DumpObjComponentsTest(unittest.TestCase):
    def setUp(self):
        self.module = load_module()
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def write_two_component_mesh(self) -> Path:
        mesh = self.root / "mesh.veg.obj"
        mesh.write_text(
            "\n".join(
                [
                    "v 0 0 0",
                    "v 1 0 0",
                    "v 0 1 0",
                    "v 3 0 0",
                    "v 4 0 0",
                    "v 3 1 0",
                    "f 1 2 3",
                    "f 4 5 6",
                    "",
                ]
            )
        )
        return mesh

    def test_dump_components_sorts_by_size_and_writes_local_indices(self):
        mesh = self.write_two_component_mesh()
        output_root = self.root / "components"

        records = self.module.dump_components_for_file(
            input_path=mesh,
            output_dir=output_root,
            overwrite=False,
        )

        self.assertEqual(len(records), 2)
        self.assertEqual([record.triangles for record in records], [1, 1])
        first_obj = output_root / "mesh_component_0.obj"
        self.assertTrue(first_obj.exists())
        text = first_obj.read_text()
        self.assertIn("v 0 0 0", text)
        self.assertIn("f 1 2 3", text)
        components_json = output_root / "components.json"
        self.assertTrue(components_json.exists())


if __name__ == "__main__":
    unittest.main()
