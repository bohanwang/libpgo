#!/usr/bin/env python3
"""Unit tests for the Alembic preview renderer wrapper."""

from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
RENDERER_PATH = REPO_ROOT / "scripts" / "render_abc_preview.py"


def load_renderer_module():
    spec = importlib.util.spec_from_file_location("render_abc_preview", RENDERER_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {RENDERER_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class RenderAbcPreviewTest(unittest.TestCase):
    def test_load_config_resolves_paths_relative_to_config(self) -> None:
        renderer = load_renderer_module()

        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            case_dir = root / "case"
            case_dir.mkdir()
            config_path = case_dir / "preview-render.json"
            config_path.write_text(
                json.dumps(
                    {
                        "abc": "output/abc/case.abc",
                        "frames_dir": "output/render_frames",
                        "output_gif": "../../previews/case.gif",
                        "frame_start": 5,
                        "frame_end": 25,
                        "frame_step": 2,
                        "fps": 30,
                        "gif_fps": 10,
                        "width": 640,
                        "height": 360,
                        "camera": {
                            "mode": "auto",
                            "view": [0.0, -1.0, 0.35],
                            "ortho_scale_multiplier": 2.1,
                        },
                    }
                ),
                encoding="utf-8",
            )

            config = renderer.load_config(config_path)

            self.assertEqual(config.abc, (case_dir / "output" / "abc" / "case.abc").resolve())
            self.assertEqual(config.frames_dir, (case_dir / "output" / "render_frames").resolve())
            self.assertEqual(config.output_gif, (root.parent / "previews" / "case.gif").resolve())
            self.assertEqual(config.frame_start, 5)
            self.assertEqual(config.frame_end, 25)
            self.assertEqual(config.frame_step, 2)
            self.assertEqual(config.fps, 30)
            self.assertEqual(config.gif_fps, 10)
            self.assertEqual(config.width, 640)
            self.assertEqual(config.height, 360)
            self.assertEqual(config.camera.view, (0.0, -1.0, 0.35))
            self.assertEqual(config.camera.ortho_scale_multiplier, 2.1)

    def test_build_commands_use_blender_script_and_ffmpeg_palette(self) -> None:
        renderer = load_renderer_module()

        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            config = renderer.RenderConfig(
                config_path=root / "case-render.json",
                abc=root / "case.abc",
                frames_dir=root / "frames",
                output_gif=root / "case.gif",
                frame_start=1,
                frame_end=9,
                frame_step=3,
                fps=30,
                gif_fps=10,
                width=960,
                height=540,
                samples=32,
                background=(1.0, 1.0, 1.0),
                material_color=(0.72, 0.78, 0.86, 1.0),
                camera=renderer.CameraConfig(mode="auto", view=(0.0, -1.0, 0.4), ortho_scale_multiplier=2.4),
            )

            commands = renderer.build_commands(config, Path("/Applications/Blender.app/Contents/MacOS/Blender"), Path("/usr/bin/ffmpeg"))

            self.assertEqual(len(commands), 3)
            self.assertEqual(commands[0].label, "blender")
            self.assertEqual(commands[0].argv[:3], ["/Applications/Blender.app/Contents/MacOS/Blender", "-b", "--python"])
            self.assertIn(str(REPO_ROOT / "scripts" / "blender_render_abc.py"), commands[0].argv)
            self.assertIn("--frame-step", commands[0].argv)
            self.assertIn("3", commands[0].argv)
            self.assertEqual(commands[1].label, "palette")
            self.assertEqual(commands[1].argv[0], "/usr/bin/ffmpeg")
            self.assertIn("-pattern_type", commands[1].argv)
            self.assertIn("glob", commands[1].argv)
            self.assertIn(str(root / "frames" / "frame_*.png"), commands[1].argv)
            self.assertIn("palettegen=stats_mode=diff", commands[1].argv)
            self.assertEqual(commands[2].label, "gif")
            self.assertIn("-pattern_type", commands[2].argv)
            self.assertIn("glob", commands[2].argv)
            self.assertIn("paletteuse=dither=bayer:bayer_scale=5:diff_mode=rectangle", commands[2].argv)
            self.assertEqual(commands[2].argv[-1], str(root / "case.gif"))


if __name__ == "__main__":
    unittest.main()
