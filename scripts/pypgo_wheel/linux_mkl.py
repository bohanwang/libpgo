"""Restore oneMKL dispatch filenames after auditwheel repair.

oneMKL opens its CPU-specific kernels by their original filenames at runtime.
auditwheel cannot rewrite those embedded dlopen() strings, so its normal
hash-mangling makes an otherwise self-contained wheel fail at the first
dispatched MKL operation.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from .common import MKL_DISPATCH_COMPONENTS


def run(*command: str) -> None:
    subprocess.run(command, check=True)


def restore(wheel: Path) -> None:
    wheel = wheel.resolve()
    if not wheel.is_file():
        raise FileNotFoundError(wheel)

    with tempfile.TemporaryDirectory(prefix="pypgo-mkl-dispatch-") as directory:
        temporary = Path(directory)
        unpack_dir = temporary / "unpacked"
        output_dir = temporary / "repacked"
        run(
            sys.executable,
            "-m",
            "wheel",
            "unpack",
            "-d",
            str(unpack_dir),
            str(wheel),
        )

        roots = [path for path in unpack_dir.iterdir() if path.is_dir()]
        if len(roots) != 1:
            raise RuntimeError(f"expected one unpacked wheel root, found {roots}")
        root = roots[0]

        extensions = list(root.glob("pypgo/_pypgo*.so"))
        if len(extensions) != 1:
            raise RuntimeError(f"expected one pypgo extension, found {extensions}")
        extension = extensions[0]

        library_dirs = list(root.glob("pypgo.libs"))
        if len(library_dirs) != 1:
            raise RuntimeError(
                f"expected one pypgo.libs directory, found {library_dirs}"
            )
        library_dir = library_dirs[0]

        for component in MKL_DISPATCH_COMPONENTS:
            pattern = re.compile(
                rf"^libmkl_{re.escape(component)}-[0-9a-f]+\.so\.2$"
            )
            candidates = [
                path for path in library_dir.iterdir() if pattern.match(path.name)
            ]
            if len(candidates) != 1:
                raise RuntimeError(
                    f"expected one repaired libmkl_{component} library, "
                    f"found {candidates}"
                )

            repaired = candidates[0]
            original_name = f"libmkl_{component}.so.2"
            restored = library_dir / original_name
            repaired_name = repaired.name
            repaired.rename(restored)
            run("patchelf", "--set-soname", original_name, str(restored))
            run(
                "patchelf",
                "--replace-needed",
                repaired_name,
                original_name,
                str(extension),
            )

        output_dir.mkdir()
        run(
            sys.executable,
            "-m",
            "wheel",
            "pack",
            "-d",
            str(output_dir),
            str(root),
        )
        packed = list(output_dir.glob("*.whl"))
        if len(packed) != 1:
            raise RuntimeError(f"expected one repacked wheel, found {packed}")
        staged = wheel.with_suffix(f"{wheel.suffix}.tmp")
        shutil.copyfile(packed[0], staged)
        os.replace(staged, wheel)

    print(f"restored oneMKL dispatch filenames in {wheel.name}")
