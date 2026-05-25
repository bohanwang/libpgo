"""CMake-backed setuptools entry point for pypgo."""

import os
import re
import shlex
import subprocess
import sys
import sysconfig
import platform
from pathlib import Path
import shutil
from typing import List, Optional, Union

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext

install_requires = ["numpy"]

# Convert distutils Windows platform specifiers to CMake -A arguments
PLAT_TO_CMAKE = {
    "win32": "Win32",
    "win-amd64": "x64",
    "win-arm32": "ARM",
    "win-arm64": "ARM64",
}


def cmake_bool(value: bool) -> str:
    return "ON" if value else "OFF"


def conda_native_prefix(is_windows: bool) -> Optional[Path]:
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if not conda_prefix:
        return None
    prefix = Path(conda_prefix)
    return prefix / "Library" if is_windows else prefix


def cmake_prefix_path_entries(raw: Optional[str]) -> List[str]:
    if not raw:
        return []
    return [entry for entry in raw.split(";") if entry]


def add_unique_path(paths: List[str], path: Optional[Union[Path, str]]) -> None:
    if not path:
        return
    path_str = os.fspath(path)
    if path_str and path_str not in paths:
        paths.append(path_str)


def path_has_mkl(prefix: Path) -> bool:
    include_dir = prefix / "include"
    return (include_dir / "mkl.h").exists() or (include_dir / "mkl_version.h").exists()


def has_mkl_hint(is_windows: bool, cmake_prefix_paths: List[str]) -> bool:
    if os.environ.get("MKLROOT") and path_has_mkl(Path(os.environ["MKLROOT"])):
        return True

    candidates = [Path(path) for path in cmake_prefix_paths]
    conda_prefix = conda_native_prefix(is_windows)
    if conda_prefix is not None:
        candidates.append(conda_prefix)

    for candidate in candidates:
        if path_has_mkl(candidate):
            return True

    return False


def python_mkl_enabled(is_macos: bool, is_windows: bool, cmake_prefix_paths: List[str]) -> bool:
    mode = os.environ.get("PGO_PYTHON_USE_MKL", "auto").strip().lower()
    if mode not in {"auto", "on", "off"}:
        raise RuntimeError("PGO_PYTHON_USE_MKL must be one of: auto, on, off")

    if is_macos:
        if mode == "on":
            raise RuntimeError("PGO_PYTHON_USE_MKL=on is not supported on macOS; use PGO_PYTHON_USE_MKL=off or auto.")
        return False

    if mode == "off":
        return False

    found_hint = has_mkl_hint(is_windows, cmake_prefix_paths)
    if mode == "on" and not found_hint:
        raise RuntimeError(
            "PGO_PYTHON_USE_MKL=on requires an MKL hint. Install mkl-devel and set "
            "MKLROOT, CONDA_PREFIX, or CMAKE_PREFIX_PATH to the native MKL prefix."
        )
    return found_hint


# A CMakeExtension needs a sourcedir instead of a file list.
# The name must be the _single_ output extension from the CMake build.
# If you need multiple extensions, see scikit-build.
class CMakeExtension(Extension):
    def __init__(self, name: str, sourcedir: str = "") -> None:
        super().__init__(name, sources=[])
        self.sourcedir = os.fspath(Path(sourcedir).resolve())


class CMakeBuild(build_ext):
    def build_extension(self, ext: CMakeExtension) -> None:
        # Must be in this form due to bug in .resolve() only fixed in Python 3.10+
        ext_fullpath = Path.cwd() / self.get_ext_fullpath(ext.name)
        extdir = ext_fullpath.parent.resolve()

        # Using this requires trailing slash for auto-detection & inclusion of
        # auxiliary "native" libs

        debug = int(os.environ.get("DEBUG", 0)) if self.debug is None else self.debug
        cfg = "Debug" if debug else "Release"
        is_windows = sys.platform.startswith("win")
        is_macos = sys.platform == "darwin"
        preset_like = "base"

        # CMake lets you override the generator - we need to check this.
        # Can be set with Conda-Build, for example.
        cmake_generator = os.environ.get("CMAKE_GENERATOR", "")

        cmake_prefix_paths = cmake_prefix_path_entries(os.environ.get("CMAKE_PREFIX_PATH"))
        conda_prefix = conda_native_prefix(is_windows)
        if conda_prefix is not None:
            add_unique_path(cmake_prefix_paths, conda_prefix)
            if "MKLROOT" not in os.environ and path_has_mkl(conda_prefix):
                os.environ["MKLROOT"] = os.fspath(conda_prefix)

        use_mkl = python_mkl_enabled(is_macos, is_windows, cmake_prefix_paths)

        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}{os.sep}",
            f"-DPython_EXECUTABLE={sys.executable}",
            f"-DCMAKE_BUILD_TYPE={cfg}",  # not used on MSVC, but no harm
            f"-DPGO_USE_MKL={cmake_bool(use_mkl)}",
            "-DPGO_CHECK_CONDA=ON",
            "-DPGO_ENABLE_FULL=ON",
            "-DPGO_BUILD_SUBPROJECTS=ON",
            "-DPGO_ENABLE_PYTHON=ON",
            "-DPGO_BUILD_C_API=OFF",
            "-DPGO_ENABLE_ALEMBIC=ON",
            "-DPGO_ENABLE_GMSH=OFF",
            "-DPGO_TET_MESHER_USE_TET_WILD=OFF",
            "-DPGO_ENABLE_OPENVDB=OFF",
        ]

        if cmake_prefix_paths:
            cmake_args.append(f"-DCMAKE_PREFIX_PATH={';'.join(cmake_prefix_paths)}")

        build_args = []
        # Adding CMake arguments set as environment variable
        # (needed e.g. to build for ARM OSx on conda-forge)
        if "CMAKE_ARGS" in os.environ:
            cmake_args += [item for item in shlex.split(os.environ["CMAKE_ARGS"]) if item]

        # In this example, we pass in the version to C++. You might not need to.
        cmake_args += [f"-DPYPGO_VERSION_INFO={self.distribution.get_version()}"]

        if self.compiler.compiler_type != "msvc":
            # Using Ninja-build since it a) is available as a wheel and b)
            # multithreads automatically. MSVC would require all variables be
            # exported for Ninja to pick it up, which is a little tricky to do.
            # Users can override the generator with CMAKE_GENERATOR in CMake
            # 3.15+.
            if not cmake_generator or cmake_generator == "Ninja":
                try:
                    import ninja

                    ninja_executable_path = Path(ninja.BIN_DIR) / "ninja"
                    cmake_args += [
                        "-GNinja",
                        f"-DCMAKE_MAKE_PROGRAM:FILEPATH={ninja_executable_path}",
                    ]
                except ImportError:
                    pass

        else:
            # Single config generators are handled "normally"
            single_config = any(x in cmake_generator for x in {"NMake", "Ninja"})

            # CMake allows an arch-in-generator style for backward compatibility
            contains_arch = any(x in cmake_generator for x in {"ARM", "Win64"})

            # Specify the arch if using MSVC generator, but only if it doesn't
            # contain a backward-compatibility arch spec already in the
            # generator name.
            if not single_config and not contains_arch:
                cmake_args += ["-A", PLAT_TO_CMAKE[self.plat_name]]

            # Multi-config generators have a different way to specify configs
            if not single_config:
                cmake_args += [f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{cfg.upper()}={extdir}"]
                build_args += ["--config", cfg]

        if sys.platform.startswith("darwin"):
            # Cross-compile support for macOS - respect ARCHFLAGS if set
            archs = re.findall(r"-arch (\S+)", os.environ.get("ARCHFLAGS", ""))
            if archs:
                cmake_args += ["-DCMAKE_OSX_ARCHITECTURES={}".format(";".join(archs))]

        # Set CMAKE_BUILD_PARALLEL_LEVEL to control the parallel build level
        # across all generators.
        if "CMAKE_BUILD_PARALLEL_LEVEL" not in os.environ:
            # self.parallel is a Python 3 only way to set parallel jobs by hand
            # using -j in the build_ext call, not supported by pip or PyPA-build.
            if hasattr(self, "parallel") and self.parallel:
                # CMake 3.12+ only.
                build_args += [f"-j{self.parallel}"]

        default_build_dir = (
            Path.cwd()
            / "build"
            / f"pypgo-conda-{preset_like}-{sysconfig.get_platform()}-{sys.implementation.cache_tag}-{cfg.lower()}"
        )
        build_temp = Path(os.environ.get("PGO_PYTHON_BUILD_DIR", default_build_dir))
        if not build_temp.exists():
            build_temp.mkdir(parents=True)

        # e.g., copy a known DLL into the same folder as the built .pyd/.so
        # for the extension named mypackage._example
        if "Windows" in platform.platform():
            ext_build_path = self.get_ext_fullpath(ext.name)
            ext_dir = os.path.dirname(os.path.abspath(ext_build_path))
            if not os.path.exists(ext_dir):
                os.makedirs(ext_dir, exist_ok=True)

            third_party_folder = (Path.cwd() / "third-party").resolve()

            for folder in [f"{third_party_folder}/gmp-msvc/release", f"{third_party_folder}/mpfr-msvc/release"]:
                for file in os.listdir(folder):
                    full_filename = os.path.join(folder, file)
                    if os.path.isfile(full_filename) and file.lower().endswith(".dll"):
                        dest_dll = os.path.join(ext_dir, file)
                        print(f"copying {full_filename} to {dest_dll}")
                        shutil.copyfile(full_filename, dest_dll)

        subprocess.run(["cmake", ext.sourcedir, *cmake_args], cwd=build_temp, check=True)
        subprocess.run(["cmake", "--build", ".", "--target", "pypgo", *build_args], cwd=build_temp, check=True)

# The information here can also be placed in setup.cfg - better separation of
# logic and declaration, and simpler if you include description/version in a file.
setup(
    name="pypgo",
    version="0.0.3",
    author="Bohan Wang",
    author_email="wangbh11@gmail.com",
    description="build pypgo",
    long_description="",
    ext_modules=[CMakeExtension("pypgo")],
    cmdclass={"build_ext": CMakeBuild},
    zip_safe=False,
    install_requires=install_requires,
    extras_require={"test": ["pytest>=6.0"]},
    python_requires=">=3.9",
    # Tell setuptools to include extra non-Python files in the wheel
    # include_package_data=include_package_data,  # needs a MANIFEST.in or package_data below
    # package_dir=package_dir,
    # One way: use package_data
    # package_data=package_data,
)
