#!/usr/bin/env python3
"""Exercise an installed pypgo wheel and verify the loaded BLAS runtime."""

from __future__ import annotations

import argparse
import ctypes
import importlib
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

from project_version import read_project_version


def windows_modules() -> str:
    from ctypes import wintypes

    kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
    psapi = ctypes.WinDLL("psapi", use_last_error=True)

    get_current_process = kernel32.GetCurrentProcess
    get_current_process.argtypes = ()
    get_current_process.restype = wintypes.HANDLE

    enum_process_modules = psapi.EnumProcessModules
    enum_process_modules.argtypes = (
        wintypes.HANDLE,
        ctypes.POINTER(wintypes.HMODULE),
        wintypes.DWORD,
        ctypes.POINTER(wintypes.DWORD),
    )
    enum_process_modules.restype = wintypes.BOOL

    get_module_file_name = psapi.GetModuleFileNameExW
    get_module_file_name.argtypes = (
        wintypes.HANDLE,
        wintypes.HMODULE,
        wintypes.LPWSTR,
        wintypes.DWORD,
    )
    get_module_file_name.restype = wintypes.DWORD

    process = get_current_process()
    modules = (wintypes.HMODULE * 4096)()
    needed = wintypes.DWORD()
    if not enum_process_modules(
        process, modules, ctypes.sizeof(modules), ctypes.byref(needed)
    ):
        raise ctypes.WinError(ctypes.get_last_error())
    paths = []
    module_count = min(needed.value // ctypes.sizeof(wintypes.HMODULE), len(modules))
    for module in modules[:module_count]:
        buffer = ctypes.create_unicode_buffer(32768)
        if get_module_file_name(process, module, buffer, len(buffer)):
            paths.append(buffer.value)
    return "\n".join(paths)


def loaded_modules(target: str) -> str:
    if target == "linux":
        return Path("/proc/self/maps").read_text(encoding="utf-8")
    if target == "macos":
        return subprocess.run(
            ["vmmap", str(os.getpid())],
            check=True,
            capture_output=True,
            text=True,
        ).stdout
    return windows_modules()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--platform", choices=("linux", "macos", "windows"), required=True)
    parser.add_argument("--source-dir", type=Path, required=True)
    args = parser.parse_args()

    import pypgo

    package_path = Path(pypgo.__file__).resolve()
    native_module = importlib.import_module("pypgo._pypgo")
    module_path = Path(native_module.__file__).resolve()
    source_dir = args.source_dir.resolve()
    for loaded_path in (package_path, module_path):
        if source_dir == loaded_path or source_dir in loaded_path.parents:
            raise RuntimeError(f"pypgo was imported from the source checkout: {loaded_path}")
    expected_version = read_project_version(source_dir)
    if pypgo.__version__ != expected_version:
        raise RuntimeError(f"unexpected pypgo version: {pypgo.__version__}")

    matrix = np.arange(256 * 256, dtype=np.float64).reshape(256, 256)
    product = matrix @ matrix.T
    if not np.isfinite(product).all():
        raise RuntimeError("BLAS GEMM produced non-finite values")
    solution = np.linalg.solve(
        np.array([[4.0, 1.0], [2.0, 3.0]]),
        np.array([1.0, 1.0]),
    )
    if not np.allclose(solution, [0.2, 0.2]):
        raise RuntimeError(f"LAPACK solve produced {solution}")

    vertices = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        dtype=np.float32,
    )
    tets = np.array([[0, 1, 2, 3]], dtype=np.int32)
    mesh = pypgo.create_tetmeshgeo(vertices.ravel(), tets.ravel())
    try:
        matrix_handle = pypgo.create_element_laplacian_matrix(mesh, 0, 9, 0)
        try:
            values = pypgo.sparse_matrix_get_values(matrix_handle)
            if values.size == 0 or not np.isfinite(values).all():
                raise RuntimeError("pypgo sparse operation failed")
        finally:
            pypgo.destroy_sparse_matrix(matrix_handle)
    finally:
        pypgo.destroy_tetmeshgeo(mesh)

    modules = loaded_modules(args.platform).lower()
    if args.platform == "macos":
        if "accelerate.framework" not in modules:
            raise RuntimeError("system Accelerate is not loaded")
        if "mkl" in modules or "openblas" in modules:
            raise RuntimeError("unexpected macOS BLAS runtime")
    else:
        if "mkl_tbb_thread" not in modules or "tbb" not in modules:
            raise RuntimeError("MKL TBB threading runtime is not loaded")
        if any(name in modules for name in ("mkl_intel_thread", "mkl_gnu_thread", "libiomp5")):
            raise RuntimeError("an unapproved MKL threading runtime is loaded")

    print(f"pypgo {pypgo.__version__}: {package_path} ({module_path.name})")
    print(f"{args.platform} BLAS/LAPACK runtime smoke passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
