"""
umbra: Kerr spacetime geometry and black hole horizons
Provides analytical horizon radii, ergospheres, and Bardeen shadow contours.
"""

import math
import numpy as np
from typing import Tuple, Dict, Any
from .core import get_lib

class KerrSpacetime:
    """Represents a rotating Kerr black hole spacetime with mass M=1."""

    def __init__(self, spin: float = 0.0):
        if abs(spin) > 1.0:
            raise ValueError(f"Extreme Kerr limit exceeded: |a| = {abs(spin)} > 1.0")
        self.spin = float(spin)

    @property
    def horizon_outer(self) -> float:
        """Event horizon radius r+ = 1 + sqrt(1 - a^2)."""
        return 1.0 + math.sqrt(max(1.0 - self.spin**2, 0.0))

    @property
    def horizon_inner(self) -> float:
        """Cauchy horizon radius r- = 1 - sqrt(1 - a^2)."""
        return 1.0 - math.sqrt(max(1.0 - self.spin**2, 0.0))

    def ergosphere_radius(self, theta_rad: float) -> float:
        """Ergosphere boundary radius r_E(theta) = 1 + sqrt(1 - a^2 cos^2 theta)."""
        term = max(1.0 - (self.spin * math.cos(theta_rad))**2, 0.0)
        return 1.0 + math.sqrt(term)

    def photon_orbit_radii(self) -> Tuple[float, float]:
        """Prograde and retrograde circular equatorial photon orbit radii (Bardeen 1972)."""
        a_abs = abs(self.spin)
        r_pro = 2.0 * (1.0 + math.cos(2.0 / 3.0 * math.acos(-a_abs)))
        r_ret = 2.0 * (1.0 + math.cos(2.0 / 3.0 * math.acos(a_abs)))
        return (r_pro, r_ret)

    def shadow_boundary(self, theta_obs_deg: float = 90.0, n_points: int = 200) -> Dict[str, np.ndarray]:
        """
        Computes celestial coordinates (alpha, beta) of the apparent black hole shadow boundary.
        Uses compiled Fortran kernel when available with exact NumPy fallback.
        """
        theta_rad = math.radians(theta_obs_deg)
        lib = get_lib()

        if lib is not None and hasattr(lib, 'umbra_kerr_shadow_analytic'):
            alpha = np.empty(n_points, dtype=np.float64)
            beta_plus = np.empty(n_points, dtype=np.float64)
            beta_minus = np.empty(n_points, dtype=np.float64)
            lib.umbra_kerr_shadow_analytic(self.spin, theta_rad, n_points, alpha, beta_plus, beta_minus)
            return {"alpha": alpha, "beta_plus": beta_plus, "beta_minus": beta_minus}

        # Exact analytical fallback
        a = self.spin
        if abs(a) < 1.0e-7:
            # Schwarzschild limit: circular profile of radius 3*sqrt(3)
            r_shad = 3.0 * math.sqrt(3.0)
            alpha = np.linspace(-r_shad, r_shad, n_points)
            beta_sq = np.maximum(r_shad**2 - alpha**2, 0.0)
            beta_plus = np.sqrt(beta_sq)
            beta_minus = -beta_plus
            return {"alpha": alpha, "beta_plus": beta_plus, "beta_minus": beta_minus}

        sin_th = math.sin(theta_rad)
        cos_th = math.cos(theta_rad)
        cot_th = cos_th / max(abs(sin_th), 1.0e-12)

        r_pro, r_ret = self.photon_orbit_radii()
        r_p = np.linspace(r_pro, r_ret, n_points)
        delta = r_p**2 - 2.0 * r_p + a**2

        xi = ((r_p**2 - a**2) - r_p * delta) / (a * (r_p - 1.0))
        eta = (r_p**3 * (4.0 * delta - r_p * (r_p - 1.0)**2)) / (a**2 * (r_p - 1.0)**2)

        alpha = -xi / sin_th
        term_beta = eta + a**2 * cos_th**2 - xi**2 * cot_th**2
        beta_plus = np.sqrt(np.maximum(term_beta, 0.0))
        beta_minus = -beta_plus

        return {"alpha": alpha, "beta_plus": beta_plus, "beta_minus": beta_minus}
