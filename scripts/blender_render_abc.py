#!/usr/bin/env python3
"""Blender-side Alembic frame renderer used by render_abc_preview.py."""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import bpy
from mathutils import Vector


def parse_vec(text: str, length: int) -> tuple[float, ...]:
    values = tuple(float(part) for part in text.split(","))
    if len(values) != length:
        raise ValueError(f"expected {length} comma-separated values, got {text}")
    return values


def parse_args() -> argparse.Namespace:
    argv = sys.argv
    argv = argv[argv.index("--") + 1 :] if "--" in argv else []
    parser = argparse.ArgumentParser(description="Render an Alembic animation to PNG frames.")
    parser.add_argument("--abc", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--frame-start", type=int, required=True)
    parser.add_argument("--frame-end", type=int, required=True)
    parser.add_argument("--frame-step", type=int, default=1)
    parser.add_argument("--fps", type=int, default=30)
    parser.add_argument("--width", type=int, default=960)
    parser.add_argument("--height", type=int, default=540)
    parser.add_argument("--samples", type=int, default=64)
    parser.add_argument("--background", default="1,1,1")
    parser.add_argument("--material-color", default="0.72,0.78,0.86,1")
    parser.add_argument("--camera-view", default="0,-1,0.35")
    parser.add_argument("--ortho-scale-multiplier", type=float, default=2.4)
    args = parser.parse_args(argv)
    args.background = parse_vec(args.background, 3)
    args.material_color = parse_vec(args.material_color, 4)
    args.camera_view = parse_vec(args.camera_view, 3)
    return args


def clear_scene() -> None:
    bpy.ops.object.select_all(action="SELECT")
    bpy.ops.object.delete()


def import_abc(path: Path) -> None:
    bpy.ops.wm.alembic_import(filepath=str(path), set_frame_range=False)


def mesh_objects() -> list[bpy.types.Object]:
    return [obj for obj in bpy.context.scene.objects if obj.type == "MESH"]


def scene_bounds() -> tuple[Vector, float]:
    mins = Vector((math.inf, math.inf, math.inf))
    maxs = Vector((-math.inf, -math.inf, -math.inf))

    found = False
    for obj in mesh_objects():
        found = True
        for corner in obj.bound_box:
            point = obj.matrix_world @ Vector(corner)
            mins.x = min(mins.x, point.x)
            mins.y = min(mins.y, point.y)
            mins.z = min(mins.z, point.z)
            maxs.x = max(maxs.x, point.x)
            maxs.y = max(maxs.y, point.y)
            maxs.z = max(maxs.z, point.z)

    if not found:
        raise RuntimeError("Alembic import did not create any mesh objects")

    center = (mins + maxs) * 0.5
    radius = max((maxs - mins).length * 0.5, 1e-6)
    return center, radius


def assign_material(color: tuple[float, float, float, float]) -> None:
    material = bpy.data.materials.new("PreviewMaterial")
    material.diffuse_color = color
    for obj in mesh_objects():
        obj.data.materials.clear()
        obj.data.materials.append(material)


def setup_camera(center: Vector, radius: float, view: tuple[float, float, float], scale_multiplier: float) -> None:
    view_vector = Vector(view)
    if view_vector.length == 0:
        raise ValueError("camera view vector must be non-zero")
    view_vector.normalize()

    camera_data = bpy.data.cameras.new("Camera")
    camera = bpy.data.objects.new("Camera", camera_data)
    bpy.context.collection.objects.link(camera)
    camera.location = center + view_vector * (3.0 * radius)
    direction = center - camera.location
    camera.rotation_euler = direction.to_track_quat("-Z", "Y").to_euler()
    camera_data.type = "ORTHO"
    camera_data.ortho_scale = radius * scale_multiplier
    bpy.context.scene.camera = camera


def setup_light(center: Vector, radius: float) -> None:
    light_data = bpy.data.lights.new("KeyArea", "AREA")
    light_data.energy = 500.0
    light_data.size = radius * 3.0
    light = bpy.data.objects.new("KeyArea", light_data)
    bpy.context.collection.objects.link(light)
    light.location = center + Vector((0.0, -2.0 * radius, 3.0 * radius))


def setup_render(args: argparse.Namespace) -> None:
    scene = bpy.context.scene
    scene.frame_start = args.frame_start
    scene.frame_end = args.frame_end
    scene.frame_step = args.frame_step
    scene.render.fps = args.fps
    scene.render.resolution_x = args.width
    scene.render.resolution_y = args.height
    scene.render.film_transparent = False
    scene.render.image_settings.file_format = "PNG"
    scene.render.filepath = str(args.out / "frame_")

    world = scene.world or bpy.data.worlds.new("World")
    scene.world = world
    world.color = args.background

    scene.render.engine = "BLENDER_EEVEE_NEXT" if "BLENDER_EEVEE_NEXT" in bpy.types.RenderSettings.bl_rna.properties["engine"].enum_items else "BLENDER_EEVEE"
    if hasattr(scene, "eevee"):
        scene.eevee.taa_render_samples = args.samples


def main() -> None:
    args = parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    clear_scene()
    import_abc(args.abc)
    assign_material(args.material_color)
    center, radius = scene_bounds()
    setup_camera(center, radius, args.camera_view, args.ortho_scale_multiplier)
    setup_light(center, radius)
    setup_render(args)
    bpy.ops.render.render(animation=True)


if __name__ == "__main__":
    main()
