"""
umbra: Post-processing and scientific error analysis
Computes L1, L2, L_inf error norms and convergence rates across spatial grid refinements.
"""

import numpy as np
from typing import Dict, Any, List

def compute_error_norms(numerical: np.ndarray, exact: np.ndarray) -> Dict[str, float]:
    """
    Computes discrete L1, L2, and L_infinity relative error norms between numerical and exact fields.
    """
    diff = numerical - exact
    norm_exact_l1 = float(np.sum(np.abs(exact))) + 1.0e-14
    norm_exact_l2 = float(np.sqrt(np.sum(exact**2))) + 1.0e-14
    norm_exact_linf = float(np.max(np.abs(exact))) + 1.0e-14

    l1 = float(np.sum(np.abs(diff))) / norm_exact_l1
    l2 = float(np.sqrt(np.sum(diff**2))) / norm_exact_l2
    linf = float(np.max(np.abs(diff))) / norm_exact_linf

    return {
        "l1": l1,
        "l2": l2,
        "linf": linf,
        "max_abs": float(np.max(np.abs(diff))),
    }

def compute_convergence_rate(resolutions: List[int], errors: List[float]) -> float:
    """
    Calculates numerical order of convergence:
    p = -log(E2 / E1) / log(N2 / N1) via linear regression in log-log space.
    """
    if len(resolutions) < 2 or len(errors) < 2:
        return 0.0

    log_n = np.log(resolutions)
    log_e = np.log(errors)

    # Slope of log(error) vs log(N)
    poly = np.polyfit(log_n, log_e, 1)
    rate = -poly[0] # Order of accuracy is negative of slope
    return float(rate)
