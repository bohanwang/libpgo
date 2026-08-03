"""oneMKL runtime components shared by Linux and Windows wheel contracts."""

from __future__ import annotations

from typing import Final


MKL_DISPATCH_COMPONENTS: Final = (
    "def",
    "mc3",
    "avx2",
    "avx512",
    "vml_def",
    "vml_cmpt",
    "vml_mc3",
    "vml_avx2",
    "vml_avx512",
)
MKL_DISPATCH_LIBRARY_COMPONENTS: Final = tuple(
    f"mkl_{component}" for component in MKL_DISPATCH_COMPONENTS
)
MKL_NATIVE_COMPONENTS: Final = (
    "mkl_core",
    "mkl_tbb_thread",
    *MKL_DISPATCH_LIBRARY_COMPONENTS,
)
