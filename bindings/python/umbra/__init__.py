"""
umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
Python scientific interface for relativistic magnetohydrodynamics and black hole raytracing.
"""

from .core import get_lib, find_library
from .spacetime import KerrSpacetime
from .raytracer import BlackHoleRaytracer
from .simulation import ConservativeToPrimitive, GRMHDSolver
from .analysis import compute_error_norms, compute_convergence_rate

__version__ = "1.0.0"
__author__ = "umbra development team"

__all__ = [
    "KerrSpacetime",
    "BlackHoleRaytracer",
    "ConservativeToPrimitive",
    "GRMHDSolver",
    "compute_error_norms",
    "compute_convergence_rate",
    "get_lib",
    "find_library",
    "__version__",
]
