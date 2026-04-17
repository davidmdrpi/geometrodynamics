"""
Hopf connection on S³ — the electromagnetic potential from pure geometry.

The Hopf bundle  S¹ → S³ → S²  carries a canonical connection 1-form

    A = (1/2) cos(χ) dφ

whose curvature is the electromagnetic field strength 2-form

    F = dA = -(1/2) sin(χ) dχ ∧ dφ

and whose holonomy around a full fibre at hyper-latitude χ is

    ∮ A = π cos(χ)

giving e^{iπ} = −1 at the poles (spin-½ sign flip) and trivial phase
at the equator χ = π/2 (stable orbit, zero self-energy).
"""

import numpy as np


def hopf_connection(chi: float | np.ndarray) -> float | np.ndarray:
    """φ-component of the Hopf connection at hyper-latitude χ.

    A_φ(χ) = ½ cos(χ).  Vanishes at χ = π/2 (equatorial orbit = zero
    electromagnetic self-energy).
    """
    return 0.5 * np.cos(chi)


def hopf_curvature(chi: float | np.ndarray) -> float | np.ndarray:
    """Magnitude of the field-strength 2-form |F|(χ) = ½ sin(χ).

    Maximum at χ = π/2: the equatorial orbit sits at the peak of the
    electromagnetic field emerging from pure S³ geometry.
    """
    return 0.5 * np.sin(chi)


def hopf_holonomy(chi: float | np.ndarray) -> float | np.ndarray:
    """Phase accumulated by parallel transport around a full Hopf fibre.

    ∮ A = (½ cos χ) · 2π = π cos(χ).

    At χ = 0: holonomy = π  → e^{iπ} = −1 (spinor sign flip).
    At χ = π/2: holonomy = 0 (trivial, stable orbit).
    """
    return np.pi * np.cos(chi)


def hopf_circle(
    base_theta: float,
    base_phi: float,
    N: int = 120,
    twist_offset: float = 0.0,
    scale: float = 1.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Stereographically projected Hopf circle in R³.

    Parameters
    ----------
    base_theta, base_phi : float
        Polar coordinates of the base point on S².
    N : int
        Number of sample points around the fibre.
    twist_offset : float
        Phase offset along the fibre.
    scale : float
        Overall scaling of the stereographic projection.

    Returns
    -------
    x, y, z : ndarray
        Cartesian coordinates of the projected circle.
    """
    psi = np.linspace(0.0, 2.0 * np.pi, N) + twist_offset
    chi2 = base_theta / 2.0
    x1 = np.cos(chi2) * np.cos((psi + base_phi) / 2.0)
    x2 = np.cos(chi2) * np.sin((psi + base_phi) / 2.0)
    x3 = np.sin(chi2) * np.cos((psi - base_phi) / 2.0)
    x4 = np.sin(chi2) * np.sin((psi - base_phi) / 2.0)
    d = np.maximum(1.0 - x4, 1e-8)
    return scale * x1 / d, scale * x2 / d, scale * x3 / d
