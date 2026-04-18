"""
First Chern number of the Hopf bundle.

c₁ = (1/2π) ∫_{S²} F = 1

The topological origin of the charge quantum: you cannot have half a unit
of charge because a bundle with c₁ = ½ does not exist on S³.
"""

import numpy as np
from scipy.integrate import trapezoid


def compute_c1(N_chi: int = 32_000) -> dict:
    """Numerically integrate the Hopf curvature over S².

    The φ dimension is integrated analytically (F is φ-independent → ×2π).
    The χ dimension is integrated numerically with the trapezoidal rule.

    Returns
    -------
    dict with keys:
        c1_chiphi   : integrated value in dχ∧dφ orientation
        c1_abs      : |c₁|  (should be 1.0)
        err_abs     : |  |c₁| − 1  |
        analytic    : exact analytic result (−1)
        method      : integration method description
    """
    chi = np.linspace(0.0, np.pi, N_chi)
    F_chi = -0.5 * np.sin(chi)  # F integrated over φ gives F_chi × 2π
    c1_chiphi = float(trapezoid(F_chi, chi))

    return dict(
        c1_chiphi=c1_chiphi,
        c1_phichi=-c1_chiphi,
        c1_abs=abs(c1_chiphi),
        err_abs=abs(abs(c1_chiphi) - 1.0),
        analytic=-1.0,
        method=f"trapezoid_1D_N{N_chi}",
    )
