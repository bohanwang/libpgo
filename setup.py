"""CMake-backed setuptools entry point for pypgo."""

import os
import re
import shlex
import subprocess
import sys
import sysconfig
from pathlib import Path
import shutil
from typing import List, Optional, Union

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext

install_requires = ["numpy"]


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
        preset_name = os.environ.get("PGO_PYTHON_CMAKE_PRESET", "python-build")

        cmake_prefix_paths = cmake_prefix_path_entries(os.environ.get("CMAKE_PREFIX_PATH"))
        conda_prefix = conda_native_prefix(is_windows)
        if conda_prefix is not None:
            add_unique_path(cmake_prefix_paths, conda_prefix)

        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}{os.sep}",
            f"-DPython_EXECUTABLE={sys.executable}",
            f"-DCMAKE_BUILD_TYPE={cfg}",  # not used on MSVC, but no harm
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

        try:
            import ninja

            cmake_args.append(f"-DCMAKE_MAKE_PROGRAM:FILEPATH={Path(ninja.BIN_DIR) / 'ninja'}")
        except ImportError:
            pass

        if sys.platform.startswith("darwin"):
            # Cross-compile support for macOS - respect ARCHFLAGS if set
            archs = re.findall(r"-arch (\S+)", os.environ.get("ARCHFLAGS", ""))
            if archs:
                cmake_args += ["-DCMAKE_OSX_ARCHITECTURES={}".format(";".join(archs))]

        # Set CMAKE_BUILD_PARALLEL_LEVEL to control the parallel build level
        # across all generators.
        parallel = getattr(self, "parallel", None)
        if "CMAKE_BUILD_PARALLEL_LEVEL" not in os.environ and parallel:
            build_args += [f"-j{parallel}"]

        default_build_dir = (
            Path.cwd()
            / "build"
            / f"pypgo-conda-{preset_name}-{sysconfig.get_platform()}-{sys.implementation.cache_tag}-{cfg.lower()}"
        )
        build_temp = Path(os.environ.get("PGO_PYTHON_BUILD_DIR", default_build_dir))
        build_temp = build_temp.resolve()
        if not build_temp.exists():
            build_temp.mkdir(parents=True)

        if is_windows:
            ext_dir = (Path.cwd() / self.get_ext_fullpath(ext.name)).resolve().parent
            ext_dir.mkdir(parents=True, exist_ok=True)

            third_party_folder = Path.cwd() / "third-party"
            for folder in ["gmp-msvc/release", "mpfr-msvc/release"]:
                for dll in (third_party_folder / folder).glob("*.dll"):
                    dest_dll = ext_dir / dll.name
                    print(f"copying {dll} to {dest_dll}")
                    shutil.copyfile(dll, dest_dll)

        subprocess.run(
            ["cmake", "--preset", preset_name, "-B", os.fspath(build_temp), *cmake_args],
            cwd=ext.sourcedir,
            check=True,
        )
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
