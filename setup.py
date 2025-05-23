"""
modified from pybind11 example
"""

import os
import re
import subprocess
import sys
import platform
from pathlib import Path
import shutil

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext

if "macOS" not in platform.platform():
    install_requires=["tbb", "tbb-devel", "mkl", "mkl-devel", "mkl-include"]
else:
    install_requires=[]


# Convert distutils Windows platform specifiers to CMake -A arguments
PLAT_TO_CMAKE = {
    "win32": "Win32",
    "win-amd64": "x64",
    "win-arm32": "ARM",
    "win-arm64": "ARM64",
}


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

        # CMake lets you override the generator - we need to check this.
        # Can be set with Conda-Build, for example.
        cmake_generator = os.environ.get("CMAKE_GENERATOR", "")

        # Set Python_EXECUTABLE instead if you use PYBIND11_FINDPYTHON
        # EXAMPLE_VERSION_INFO shows you how to pass a value into the C++ code
        # from Python.
        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}{os.sep}",
            f"-DPYTHON_EXECUTABLE={sys.executable}",
            f"-DCMAKE_BUILD_TYPE={cfg}",  # not used on MSVC, but no harm
        ]

        cmake_args += [f"-DPGO_ENABLE_PYTHON=1", f"-DPGO_BUILD_SUBPROJECTS=1"]

        if "macOS" in platform.platform():
            cmake_args += [
                f"-DPGO_USE_MKL=0",
            ]
        else:
            cmake_args += [
                f"-DPGO_USE_MKL=1",
            ]

        # enable mkl
        if "CONDA_PREFIX" in os.environ:
            if "Windows" in platform.platform():
                os.environ["MKLROOT"] = os.path.join(os.environ["CONDA_PREFIX"], "Library")
            else:
                os.environ["MKLROOT"] = os.environ["CONDA_PREFIX"]

        build_args = []
        # Adding CMake arguments set as environment variable
        # (needed e.g. to build for ARM OSx on conda-forge)
        if "CMAKE_ARGS" in os.environ:
            cmake_args += [item for item in os.environ["CMAKE_ARGS"].split(" ") if item]

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

        build_temp = Path(self.build_temp) / ext.name
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
        # subprocess.run(["cmake", "--build", ".", "--target", "pgo_c", *build_args], cwd=build_temp, check=True)

        # if "CONDA_PREFIX" in os.environ:
        #     if "Windows" in platform.platform():
        #         install_prefix = os.path.join(os.environ["CONDA_PREFIX"], "Library")
        #     else:
        #         install_prefix = os.environ["CONDA_PREFIX"]

        #     subprocess.run(["cmake", "--install", ".", "--prefix", install_prefix, *build_args], cwd=build_temp, check=True)


# The information here can also be placed in setup.cfg - better separation of
# logic and declaration, and simpler if you include description/version in a file.
setup(
    name="pypgo",
    version="0.0.2",
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
