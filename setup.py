import os
import re
import shlex
import subprocess
import sys
from pathlib import Path

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext


class CMakeExtension(Extension):
    def __init__(self, name: str, sourcedir: str = "") -> None:
        super().__init__(name, sources=[])
        self.sourcedir = os.fspath(Path(sourcedir).resolve())


class CMakeBuild(build_ext):
    def build_extension(self, ext: CMakeExtension) -> None:
        ext_fullpath = Path.cwd() / self.get_ext_fullpath(ext.name)
        extdir = ext_fullpath.parent.resolve()

        debug = int(os.environ.get("DEBUG", 0)) if self.debug is None else self.debug
        cfg = "Debug" if debug else "Release"
        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}{os.sep}",
            f"-DPYTHON_EXECUTABLE={sys.executable}",
            f"-DCMAKE_BUILD_TYPE={cfg}",
        ]
        if sys.platform == "win32":
            # The release workflow lets delvewheel vendor the runtime DLLs.
            # Direct CMake builds keep the default source-build staging enabled.
            cmake_args += ["-DPGO_STAGE_WINDOWS_RUNTIME=OFF"]

        build_args = []
        if "CMAKE_ARGS" in os.environ:
            cmake_args += shlex.split(os.environ["CMAKE_ARGS"], posix=sys.platform != "win32")

        if sys.platform.startswith("darwin"):
            # Cross-compile support for macOS - respect ARCHFLAGS if set
            archs = re.findall(r"-arch (\S+)", os.environ.get("ARCHFLAGS", ""))
            if archs:
                cmake_args += ["-DCMAKE_OSX_ARCHITECTURES={}".format(";".join(archs))]

        if "CMAKE_BUILD_PARALLEL_LEVEL" not in os.environ:
            if hasattr(self, "parallel") and self.parallel:
                build_args += [f"-j{self.parallel}"]

        build_temp = Path(self.build_temp) / f"{ext.name}-ninja"
        build_temp.mkdir(parents=True, exist_ok=True)

        subprocess.run(
            [
                "cmake",
                "--preset",
                "pypgo-wheel",
                "-S",
                ext.sourcedir,
                "-B",
                str(build_temp),
                *cmake_args,
            ],
            cwd=ext.sourcedir,
            check=True,
        )
        subprocess.run(["cmake", "--build", ".", "--target", "pypgo", *build_args], cwd=build_temp, check=True)

setup(
    ext_modules=[CMakeExtension("pypgo")],
    cmdclass={"build_ext": CMakeBuild},
    zip_safe=False,
)
