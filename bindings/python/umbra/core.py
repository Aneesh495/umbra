"""
umbra: Python Foreign Function Interface (FFI) bindings
Loads the compiled libumbra shared library and binds C ABI entrypoints.
"""

import os
import sys
import ctypes
import numpy as np
from typing import Optional

_LIB_UMBRA: Optional[ctypes.CDLL] = None

def find_library() -> Optional[str]:
    """Search common build and installation locations for libumbra."""
    candidates = [
        os.path.join(os.path.dirname(__file__), "..", "..", "..", "build", "libumbra.dylib"),
        os.path.join(os.path.dirname(__file__), "..", "..", "..", "build", "libumbra.so"),
        os.path.join(os.path.dirname(__file__), "..", "..", "..", "lib", "libumbra.dylib"),
        os.path.join(os.path.dirname(__file__), "..", "..", "..", "lib", "libumbra.so"),
        "libumbra.dylib",
        "libumbra.so",
    ]
    for p in candidates:
        abs_p = os.path.abspath(p)
        if os.path.isfile(abs_p):
            return abs_p
    return None

def get_lib() -> Optional[ctypes.CDLL]:
    """Returns the loaded ctypes CDLL handle or None if not yet compiled."""
    global _LIB_UMBRA
    if _LIB_UMBRA is not None:
        return _LIB_UMBRA
    path = find_library()
    if path:
        try:
            _LIB_UMBRA = ctypes.CDLL(path)
            _bind_signatures(_LIB_UMBRA)
        except OSError:
            _LIB_UMBRA = None
    return _LIB_UMBRA

def _bind_signatures(lib: ctypes.CDLL) -> None:
    """Configures C types and signatures for all exported ABI symbols."""
    # Version
    lib.umbra_version_major.argtypes = []
    lib.umbra_version_major.restype = ctypes.c_int32

    lib.umbra_version_minor.argtypes = []
    lib.umbra_version_minor.restype = ctypes.c_int32

    lib.umbra_version_patch.argtypes = []
    lib.umbra_version_patch.restype = ctypes.c_int32

    # Schwarzschild shadow radius
    lib.umbra_kerr_shadow_schwarzschild_radius.argtypes = []
    lib.umbra_kerr_shadow_schwarzschild_radius.restype = ctypes.c_double

    # Kerr shadow analytic boundary
    lib.umbra_kerr_shadow_analytic.argtypes = [
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_int32,
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ]
    lib.umbra_kerr_shadow_analytic.restype = None

    # con2prim 1d
    lib.umbra_c_con2prim_1d.argtypes = [
        ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
        ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
        ctypes.c_double,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
    ]
    lib.umbra_c_con2prim_1d.restype = ctypes.c_int32

    # con2prim 2d
    lib.umbra_c_con2prim_2d.argtypes = [
        ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
        ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
        ctypes.c_double,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
    ]
    lib.umbra_c_con2prim_2d.restype = ctypes.c_int32
