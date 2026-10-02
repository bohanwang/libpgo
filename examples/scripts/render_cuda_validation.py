#!/usr/bin/env python3
"""Render saved simulation OBJ frames and their configured IPC obstacles.

Run with a Python environment containing Blender's bpy, NumPy and Pillow:

  python render_cuda_validation.py case.json --device 1 --duration 8

The meshes are the actual simulator output. This script does not simulate,
interpolate geometry, or modify the input files. It fits a fixed or tracking
camera to the saved trajectory, then encodes images with ffmpeg/libx264.
Obstacle motion follows the simulation's frame/(num-timestep-1) convention.
"""

import argparse
import csv
import json
import math
import re
import shutil
import subprocess
from pathlib import Path

import bpy
import numpy as np
from mathutils import Vector
from PIL import Image, ImageDraw, ImageFont


coordinate_basis = np.asarray([[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]])


def configure_up_direction(gravity):
    global coordinate_basis
    gravity = np.asarray(gravity, dtype=float)
    if gravity.shape != (3,) or not np.isfinite(gravity).all():
        raise ValueError("Gravity must be a finite three-vector")
    norm = np.linalg.norm(gravity)
    up = -gravity / norm if norm > 0 else np.asarray([0.0, 1.0, 0.0])
    right = np.asarray([1.0, 0.0, 0.0])
    right -= up * np.dot(right, up)
    if np.linalg.norm(right) < 1e-8:
        right = np.asarray([0.0, 1.0, 0.0])
        right -= up * np.dot(right, up)
    right /= np.linalg.norm(right)
    coordinate_basis = np.asarray([right, np.cross(up, right), up])
    return up.tolist()


def read_obj(path):
    vertices, faces = [], []
    with Path(path).open() as source:
        for line in source:
            fields = line.split()
            if not fields:
                continue
            if fields[0] == "v":
                vertices.append([float(x) for x in fields[1:4]])
            elif fields[0] == "f":
                indices = [int(x.split("/")[0]) for x in fields[1:]]
                indices = [x - 1 if x > 0 else len(vertices) + x for x in indices]
                for j in range(1, len(indices) - 1):
                    faces.append([indices[0], indices[j], indices[j + 1]])
    if not vertices or not faces:
        raise ValueError(f"No triangle mesh in {path}")
    return np.asarray(vertices, dtype=float), np.asarray(faces, dtype=np.int32)


def blender_coordinates(points):
    """Rotate the gravity-defined simulation up to Blender +Z, preserving orientation."""
    return points @ coordinate_basis.T


def make_material(name, color, roughness):
    material = bpy.data.materials.new(name)
    material.use_nodes = True
    bsdf = material.node_tree.nodes["Principled BSDF"]
    bsdf.inputs["Base Color"].default_value = (*color, 1.0)
    bsdf.inputs["Roughness"].default_value = roughness
    bsdf.inputs["Metallic"].default_value = 0.0
    bsdf.inputs["Specular IOR Level"].default_value = 0.35
    return material


def make_mesh(name, vertices, faces, material, smooth=False):
    mesh = bpy.data.meshes.new(name)
    mesh.from_pydata(blender_coordinates(vertices).tolist(), [], faces.tolist())
    mesh.update()
    for face in mesh.polygons:
        face.use_smooth = smooth
    if smooth:
        mesh.set_sharp_from_angle(angle=math.radians(35))
    obj = bpy.data.objects.new(name, mesh)
    bpy.context.collection.objects.link(obj)
    mesh.materials.append(material)
    return obj


def set_vertices(obj, vertices):
    if len(vertices) != len(obj.data.vertices):
        raise ValueError("OBJ trajectory changes the number of vertices")
    obj.data.vertices.foreach_set("co", blender_coordinates(vertices).astype(np.float32).ravel())
    obj.data.update()


def set_render_device(scene, device_index, backend):
    if backend == "CPU":
        scene.cycles.device = "CPU"
        return "CPU"
    preferences = bpy.context.preferences.addons["cycles"].preferences
    preferences.compute_device_type = backend
    preferences.get_devices()
    candidates = [device for device in preferences.devices if device.type == backend]
    if not candidates or device_index < 0 or device_index >= len(candidates):
        raise RuntimeError(f"Requested {backend} render device {device_index}; available devices: {[(d.name, d.type) for d in preferences.devices]}")
    selected = candidates[device_index]
    for device in preferences.devices:
        device.use = device.id == selected.id
    scene.cycles.device = "GPU"
    return f"{backend}: {selected.name} ({device_index})"


def area_light(name, location, target, energy, size):
    light = bpy.data.lights.new(name, "AREA")
    light.energy, light.shape, light.size = energy, "DISK", size
    obj = bpy.data.objects.new(name, light)
    bpy.context.collection.objects.link(obj)
    obj.location = location
    obj.rotation_euler = (Vector(target) - Vector(location)).to_track_quat("-Z", "Y").to_euler()


def setup_scene(args, lower, upper):
    scene = bpy.context.scene
    scene.render.engine = "CYCLES"
    scene.cycles.samples = args.samples
    scene.cycles.use_denoising = True
    scene.cycles.max_bounces = 5
    scene.cycles.diffuse_bounces = 3
    scene.cycles.glossy_bounces = 2
    scene.render.use_persistent_data = True
    scene.render.resolution_x, scene.render.resolution_y = args.width, args.height
    scene.render.resolution_percentage = 100
    scene.render.image_settings.file_format = "PNG"
    scene.render.image_settings.color_mode = "RGB"
    scene.render.film_transparent = False
    scene.view_settings.view_transform = "AgX"
    device = set_render_device(scene, args.device, args.backend)
    world = bpy.data.worlds.new("studio background")
    world.use_nodes = True
    world.node_tree.nodes["Background"].inputs["Color"].default_value = (0.83, 0.88, 0.94, 1)
    world.node_tree.nodes["Background"].inputs["Strength"].default_value = 0.65
    scene.world = world

    corners = blender_coordinates(np.asarray([(x, y, z) for x in (lower[0], upper[0]) for y in (lower[1], upper[1]) for z in (lower[2], upper[2])]))
    center = Vector(corners.mean(axis=0))
    extent = max(float(np.ptp(corners, axis=0).max()), 1e-6)
    azimuth, elevation = math.radians(args.azimuth), math.radians(args.elevation)
    view = Vector((math.sin(azimuth) * math.cos(elevation), -math.cos(azimuth) * math.cos(elevation), math.sin(elevation)))
    camera_data = bpy.data.cameras.new("trajectory camera")
    camera_data.type = "ORTHO"
    camera = bpy.data.objects.new("trajectory camera", camera_data)
    bpy.context.collection.objects.link(camera)
    camera.location = center + view * (5 * extent)
    camera.rotation_euler = (-view).to_track_quat("-Z", "Y").to_euler()
    rotation = camera.rotation_euler.to_matrix()
    right, up = rotation.col[0], rotation.col[1]
    half_width = max(abs((Vector(corner) - center).dot(right)) for corner in corners)
    half_height = max(abs((Vector(corner) - center).dot(up)) for corner in corners)
    # Reserve space for the header and the time strip so neither covers meshes.
    camera_data.ortho_scale = 2 * max(half_width / 0.90, half_height * args.width / args.height / 0.72)
    camera.location += up * (camera_data.ortho_scale * 0.025)
    camera_data.clip_start = extent * 0.001
    camera_data.clip_end = extent * 30
    scene.camera = camera
    area_light("soft key", center + Vector((2, -3, 4)) * extent, center, 700 * extent ** 2, 3.0 * extent)
    area_light("soft fill", center + Vector((-3, -1, 2)) * extent, center, 400 * extent ** 2, 3.0 * extent)
    area_light("rim", center + Vector((1, 3, 3)) * extent, center, 600 * extent ** 2, 2.0 * extent)
    return scene, device


def font(size, bold=False):
    name = "DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf"
    path = Path("/usr/share/fonts/truetype/dejavu") / name
    return ImageFont.truetype(str(path), size) if path.exists() else ImageFont.load_default()


def annotate(path, title, subtitle, frame, simulation_time, contacts, speed, static=False):
    with Image.open(path) as source:
        image = source.convert("RGBA")
    width, height = image.size
    scale = width / 960
    canvas = Image.new("RGBA", image.size)
    draw = ImageDraw.Draw(canvas)
    pad = int(22 * scale)
    draw.rounded_rectangle((pad, pad, width - pad, pad + int(80 * scale)), radius=int(10 * scale), fill=(249, 251, 253, 231))
    draw.text((pad + int(18 * scale), pad + int(10 * scale)), title, font=font(int(25 * scale), True), fill=(24, 41, 54))
    draw.text((pad + int(18 * scale), pad + int(47 * scale)), subtitle, font=font(int(16 * scale)), fill=(77, 94, 106))
    draw.rounded_rectangle((pad, height - pad - int(45 * scale), width - pad, height - pad), radius=int(8 * scale), fill=(249, 251, 253, 233))
    contact_label = f"  |  contact pairs {contacts}" if contacts is not None else ""
    label = f"Static equilibrium{contact_label}" if static else f"t = {simulation_time:.3f} s  |  frame {frame}{contact_label}  |  playback {speed:.2f}x"
    draw.text((pad + int(16 * scale), height - pad - int(32 * scale)), label, font=font(int(15 * scale)), fill=(35, 53, 67))
    Image.alpha_composite(image, canvas).convert("RGB").save(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path)
    parser.add_argument("--output-dir", type=Path, help="Override directory containing ret*.obj")
    parser.add_argument("--video", type=Path, help="Default: output-dir/cuda_validation.mp4")
    parser.add_argument("--title", help="Caption, defaults to config filename")
    parser.add_argument("--subtitle", help="Status or context shown below the title")
    parser.add_argument("--preview", type=Path, help="Override the preview PNG path")
    parser.add_argument("--width", type=int, default=960)
    parser.add_argument("--height", type=int, default=720)
    parser.add_argument("--samples", type=int, default=32)
    parser.add_argument("--backend", choices=["CUDA", "OPTIX", "CPU"], default="CUDA")
    parser.add_argument("--device", type=int, default=0)
    parser.add_argument("--duration", type=float, default=8.0, help="Video playback duration in seconds")
    parser.add_argument("--max-frames", type=int, default=200, help="Evenly select at most this many saved OBJ frames")
    parser.add_argument("--track-body", action="store_true", help="Follow the body center with a fixed scale fitted to all saved body extents")
    parser.add_argument("--azimuth", type=float, default=32.0)
    parser.add_argument("--elevation", type=float, default=18.0)
    parser.add_argument("--preview-only", action="store_true", help="Render only the selected representative frame")
    parser.add_argument("--preview-frame", type=int, help="Simulation frame closest to this index, otherwise peak contact")
    parser.add_argument("--skip-existing", action="store_true", help="Resume rendering a previously fixed input trajectory")
    args = parser.parse_args()
    if args.width < 64 or args.height < 64 or args.width % 2 or args.height % 2 or args.samples < 1 or args.duration <= 0 or args.max_frames < 1:
        parser.error("Resolution must be even and >=64; samples, duration and max-frames must be positive")
    config_path = args.config.resolve()
    config = json.loads(config_path.read_text())
    if args.subtitle is None:
        body_type = "Koiter shell" if config.get("elastic-material") == "koiter-stvk" else "deformable solid"
        contact_type = "IPC contact + friction" if config.get("contact-friction-coeff", 0) > 0 else "IPC contact"
        args.subtitle = f"CUDA  /  {body_type}  /  {contact_type}"
    simulation_up = configure_up_direction(config.get("g", [0, -1, 0]))
    static = config.get("sim-type") == "static"
    output_dir = args.output_dir.resolve() if args.output_dir else (config_path.parent / config["output"]).resolve()
    scale = float(config.get("scale", 1.0))
    timestep = float(config.get("timestep", 1.0))
    total_steps = int(config.get("num-timestep", 1))
    manifest = output_dir / "newton-substeps.csv"
    records = []
    if manifest.exists():
        with manifest.open() as source:
            for row in csv.DictReader(source):
                frame, step = int(row["frame"]), int(row["substep"])
                path = output_dir / "substeps" / f"frame{frame:04d}-step{step:04d}.obj"
                time = float(row["physical_time"])
                if not path.is_file() or not math.isfinite(time) or (records and time <= records[-1]["time"]):
                    raise ValueError(f"Incomplete or unordered physical surface trajectory: {path}")
                records.append(dict(path=path, frame=frame, substep=step, time=time, contacts=int(row["active_pairs"])))
        trajectory_path = "recorded_physical_substeps"
    else:
        for path in sorted(output_dir.glob("ret[0-9]*.obj")):
            frame = int(re.search(r"(\d+)\.obj$", path.name).group(1)) * config.get("dump-interval", 1)
            records.append(dict(path=path, frame=frame, substep=None, time=0.0 if static else (frame + 1) * timestep, contacts=None))
        trajectory_path = "legacy_macroframe_snapshots"
    if static:
        trajectory_path = "static_equilibrium"
    if not records:
        raise FileNotFoundError(f"No recorded surface trajectory in {output_dir}")
    all_files = [row["path"] for row in records]
    # Select held snapshots at uniform physical times. Uniform row indices
    # would slow down portions with adaptive substeps and distort playback.
    times = np.asarray([row["time"] for row in records])
    sample_times = np.linspace(times[0], times[-1], min(len(records), args.max_frames))
    selected_indices = np.maximum(0, np.searchsorted(times, sample_times, side="right") - 1)
    samples = [records[index] for index in selected_indices]
    files = [row["path"] for row in samples]
    frame_indices = [row["frame"] for row in samples]
    lower, upper = np.full(3, np.inf), np.full(3, -np.inf)
    maximum_body_extent = np.zeros(3)
    body_centers = {}
    for path in all_files:
        vertices, _ = read_obj(path)
        vertices *= scale
        frame_lower, frame_upper = vertices.min(axis=0), vertices.max(axis=0)
        lower, upper = np.minimum(lower, frame_lower), np.maximum(upper, frame_upper)
        maximum_body_extent = np.maximum(maximum_body_extent, frame_upper - frame_lower)
        body_centers[path] = (frame_lower + frame_upper) / 2
    body_lower, body_upper = lower.copy(), upper.copy()
    padding = max(float(np.max(upper - lower)), 1e-6) * 0.15
    obstacles = []
    for index, spec in enumerate(config.get("ipc-external-objects", [])):
        vertices, faces = read_obj(config_path.parent / spec["filename"])
        vertices = vertices * float(spec.get("scale", 1.0)) + np.asarray(spec.get("initial-translation", [0, 0, 0]))
        movement = np.asarray(spec.get("movement", [0, 0, 0]), dtype=float)
        obstacles.append((vertices, faces, movement))
        # Large environments remain present, but only their local neighborhood
        # participates in camera fitting, so a huge floor does not hide the body.
        obstacle_lower = np.minimum(vertices.min(axis=0), (vertices + movement).min(axis=0))
        obstacle_upper = np.maximum(vertices.max(axis=0), (vertices + movement).max(axis=0))
        clipped_lower = np.maximum(obstacle_lower, body_lower - padding)
        clipped_upper = np.minimum(obstacle_upper, body_upper + padding)
        if np.all(clipped_lower <= clipped_upper):
            lower, upper = np.minimum(lower, clipped_lower), np.maximum(upper, clipped_upper)

    contact_counts = {}
    solver_csv = output_dir / "newton-solves.csv"
    if solver_csv.exists():
        with solver_csv.open() as source:
            contact_counts = {int(row["frame"]): int(row["active_pairs"]) for row in csv.DictReader(source)}
    else:
        log_path = next((path for path in (output_dir / "runIPCSim.log", config_path.parent / "run.log") if path.exists()), None)
        if log_path:
            contact_counts = {int(frame): int(count) for frame, count in re.findall(r"CUDA frame (\d+):.*?contacts=(\d+)", log_path.read_text())}
    preview_frame = args.preview_frame
    if preview_frame is None:
        preview_frame = max(frame_indices, key=lambda frame: contact_counts.get(frame, 0)) if contact_counts and max(contact_counts.values()) > 0 else frame_indices[len(frame_indices) // 2]
    preview_index = min(range(len(files)), key=lambda i: abs(frame_indices[i] - preview_frame))
    video_path = args.video.resolve() if args.video else output_dir / "cuda_validation.mp4"
    video_path.parent.mkdir(parents=True, exist_ok=True)
    preview_path = args.preview.resolve() if args.preview else video_path.with_name(video_path.stem + "_preview.png")
    preview_path.parent.mkdir(parents=True, exist_ok=True)
    image_dir = output_dir / ".cuda_validation_render"
    image_dir.mkdir(parents=True, exist_ok=True)
    title = args.title or config_path.stem.replace("-", " ")
    simulated_duration = times[-1] - times[0]
    playback_speed = simulated_duration / args.duration

    bpy.ops.wm.read_factory_settings(use_empty=True)
    first_vertices, first_faces = read_obj(files[0])
    body_material = make_material("deformable solid", (0.035, 0.27, 0.48), 0.32)
    body = make_mesh("simulated mesh", first_vertices * scale, first_faces, body_material, smooth=len(first_faces) > 1500)
    obstacle_material = make_material("IPC environment", (0.43, 0.27, 0.12), 0.68)
    obstacle_objects = [make_mesh(f"IPC obstacle {i}", vertices, faces, obstacle_material, smooth=len(faces) > 800) for i, (vertices, faces, _) in enumerate(obstacles)]
    framing_lower, framing_upper = lower.copy(), upper.copy()
    if args.track_body:
        framing_lower = body_centers[files[0]] - maximum_body_extent * 0.70
        framing_upper = body_centers[files[0]] + maximum_body_extent * 0.70
    scene, device = setup_scene(args, framing_lower, framing_upper)
    camera_rig = [(obj, obj.location.copy()) for obj in scene.objects if obj.type in {"CAMERA", "LIGHT"}]
    print(f"Render device: {device}; saved frames: {len(all_files)}; selected: {len(files)}; camera bounds: {framing_lower.tolist()} to {framing_upper.tolist()}", flush=True)
    selected = [preview_index] if args.preview_only else range(len(files))
    for index in selected:
        simulation_frame = frame_indices[index]
        image_path = image_dir / f"frame{index:05d}.png"
        if not (args.skip_existing and image_path.exists()):
            vertices, faces = read_obj(files[index])
            if not np.array_equal(faces, first_faces):
                raise ValueError(f"OBJ topology changes in {files[index]}")
            set_vertices(body, vertices * scale)
            if args.track_body:
                offset = Vector(blender_coordinates(np.asarray([body_centers[files[index]] - body_centers[files[0]]]))[0])
                for obj, initial_location in camera_rig:
                    obj.location = initial_location + offset
            ratio = simulation_frame / (total_steps - 1) if total_steps > 1 else 0.0
            for obj, (base, _, movement) in zip(obstacle_objects, obstacles):
                set_vertices(obj, base + ratio * movement)
            scene.render.filepath = str(image_path)
            bpy.ops.render.render(write_still=True)
            snapshot = samples[index]
            frame_label = simulation_frame if snapshot["substep"] is None else f"{simulation_frame}:{snapshot['substep']}"
            contacts = snapshot["contacts"] if snapshot["contacts"] is not None else contact_counts.get(simulation_frame)
            annotate(image_path, title, args.subtitle, frame_label, snapshot["time"], contacts, playback_speed, static=static)
        if index == preview_index:
            shutil.copyfile(image_path, preview_path)
        print(f"Rendered {index + 1}/{len(files)}: simulation frame {simulation_frame}", flush=True)
    if not args.preview_only:
        ffmpeg = shutil.which("ffmpeg")
        if not ffmpeg:
            raise RuntimeError("ffmpeg is required to encode the video")
        video_path.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run([ffmpeg, "-y", "-hide_banner", "-loglevel", "warning", "-framerate", f"{len(files) / args.duration:.8f}",
                        "-i", str(image_dir / "frame%05d.png"), "-frames:v", str(math.ceil(args.duration * 30)), "-r", "30",
                        "-c:v", "libx264", "-preset", "medium", "-crf", "19", "-pix_fmt", "yuv420p", "-movflags", "+faststart", str(video_path)], check=True)
    metadata = dict(simulation_type="static" if static else "dynamic", config=str(config_path), output_directory=str(output_dir), title=title, subtitle=args.subtitle, source_frames=len(all_files), rendered_frames=len(selected),
                    selected_simulation_frames=frame_indices, selected_physical_times=[] if static else [row["time"] for row in samples],
                    selected_physical_substeps=[] if static else [row["substep"] for row in samples], trajectory_path=trajectory_path,
                    physical_time=not static,
                    geometry_sampling="final static equilibrium; no physical timing" if static else "held recorded snapshots at uniform physical times; no geometry interpolation",
                    simulation_up_direction=simulation_up,
                    render_device=device, resolution=[args.width, args.height], samples=args.samples,
                    playback_duration_seconds=args.duration, playback_speed=playback_speed, preview_simulation_frame=frame_indices[preview_index],
                    camera_mode="track-body" if args.track_body else "fixed", camera_bounds=[framing_lower.tolist(), framing_upper.tolist()],
                    trajectory_bounds=[lower.tolist(), upper.tolist()], preview=str(preview_path), video=None if args.preview_only else str(video_path))
    video_path.with_name(video_path.stem + "_render.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(metadata, indent=2), flush=True)


if __name__ == "__main__":
    main()
