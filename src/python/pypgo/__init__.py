"""Python bindings for libpgo."""

import os as _os

if _os.name == "nt":
    # Keep the handle alive so source-build runtime DLLs staged beside the
    # extension remain available to oneMKL/oneTBB runtime loading.
    _dll_directory_handle = _os.add_dll_directory(_os.path.dirname(__file__))

from . import _pypgo

__all__ = [name for name in dir(_pypgo) if not name.startswith("_")]

from ._pypgo import *

__version__ = _pypgo.__version__
__all__.append("__version__")
