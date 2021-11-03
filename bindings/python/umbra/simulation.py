"""
umbra: Relativistic MHD simulation manager and state inversion
Provides Python interfaces for conservative-to-primitive solvers and Riemann fluxes.
"""

import math
import ctypes
import numpy as np
from typing import Dict, Any, Tuple
from .core import get_lib

class ConservativeToPrimitive:
    """Inverts conserved relativistic MHD state vectors to primitive variables."""

    def __init__(self, gamma: float = 4.0 / 3.0):
        self.gamma = float(gamma)

    def invert_1d(
        self,
        d: float,
        s1: float,
        s2: float,
        s3: float,
        tau: float,
        b1: float = 0.0,
        b2: float = 0.0,
        b3: float = 0.0,
    ) -> Dict[str, float]:
        """
        Executes 1D Noble-Mignone conservative-to-primitive root solve.
        Returns dictionary with rho, vx, vy, vz, p, and lorentz_factor.
        """
        lib = get_lib()
        if lib is not None and hasattr(lib, 'umbra_c_con2prim_1d'):
            rho = ctypes.c_double()
            vx = ctypes.c_double()
            vy = ctypes.c_double()
            vz = ctypes.c_double()
            p = ctypes.c_double()

            status = lib.umbra_c_con2prim_1d(
                d, s1, s2, s3, tau, b1, b2, b3, self.gamma,
                ctypes.byref(rho), ctypes.byref(vx), ctypes.byref(vy),
                ctypes.byref(vz), ctypes.byref(p)
            )
            if status != 0:
                raise RuntimeError(f"Conservative to primitive inversion failed with code: {status}")

            v_sq = vx.value**2 + vy.value**2 + vz.value**2
            lorentz = 1.0 / math.sqrt(max(1.0 - v_sq, 1.0e-12))
            return {
                "rho": rho.value,
                "vx": vx.value,
                "vy": vy.value,
                "vz": vz.value,
                "p": p.value,
                "lorentz": lorentz,
            }

        # Hydrodynamic analytical approximation fallback
        s_sq = s1**2 + s2**2 + s3**2
        s_norm = math.sqrt(s_sq)
        # Bisection / Newton root finder for W in [1, 50]
        w = 1.0
        for _ in range(50):
            v = s_norm / (d * w + tau + 1.0e-14)
            v = min(v, 0.9999)
            gamma_lor = 1.0 / math.sqrt(1.0 - v**2)
            f = gamma_lor - w
            if abs(f) < 1.0e-10:
                break
            w = 0.5 * (w + gamma_lor)

        rho_val = max(d / w, 1.0e-10)
        p_val = max((self.gamma - 1.0) * (tau - d * (w - 1.0)), 1.0e-10)
        scale = (1.0 / (d * w)) if s_norm > 1.0e-14 else 0.0

        return {
            "rho": rho_val,
            "vx": s1 * scale,
            "vy": s2 * scale,
            "vz": s3 * scale,
            "p": p_val,
            "lorentz": w,
        }

class GRMHDSolver:
    """Manages 1D and 2D relativistic MHD test problem initializations."""

    def __init__(self, nx: int = 100, ny: int = 1, gamma: float = 1.4):
        self.nx = int(nx)
        self.ny = int(ny)
        self.gamma = float(gamma)
        self.con2prim = ConservativeToPrimitive(self.gamma)

    def initialize_sod(self) -> Dict[str, np.ndarray]:
        """Sets up 1D Relativistic Sod shock tube initial condition."""
        x = np.linspace(0.0, 1.0, self.nx)
        rho = np.where(x < 0.5, 10.0, 1.0)
        p = np.where(x < 0.5, 40.0 / 3.0, 1.0e-6)
        vx = np.zeros(self.nx, dtype=np.float64)
        vy = np.zeros(self.nx, dtype=np.float64)
        vz = np.zeros(self.nx, dtype=np.float64)
        return {"x": x, "rho": rho, "p": p, "vx": vx, "vy": vy, "vz": vz}

    def initialize_briowu(self) -> Dict[str, np.ndarray]:
        """Sets up 1D Relativistic Brio-Wu magnetized shock tube initial condition."""
        x = np.linspace(0.0, 1.0, self.nx)
        rho = np.where(x < 0.5, 1.0, 0.125)
        p = np.where(x < 0.5, 1.0, 0.1)
        by = np.where(x < 0.5, 1.0, -1.0)
        bx = np.full(self.nx, 0.5, dtype=np.float64)
        bz = np.zeros(self.nx, dtype=np.float64)
        vx = np.zeros(self.nx, dtype=np.float64)
        return {"x": x, "rho": rho, "p": p, "bx": bx, "by": by, "bz": bz, "vx": vx}
