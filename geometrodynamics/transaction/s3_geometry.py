"""
S³ geometry utilities.

Geodesic distance, antipodal map, Green function for the Laplacian
on the 3-sphere, and the Hopf-sphere parameterisation.
"""

import numpy as np

from geometrodynamics.constants import S3_RADIUS, S3_GREEN_EPS


def nrm4(v) -> np.ndarray:
    """Normalise a 4-vector to unit length on S³."""
    v = np.asarray(v, dtype=float)
    n = np.linalg.norm(v)
    if n < 1e-15:
        raise ValueError("zero vector cannot be normalised")
    return v / n


def antipode4(p4) -> np.ndarray:
    """Antipodal point on S³: −p."""
    return -np.asarray(p4, dtype=float)


def geo4(a, b) -> float:
    """Geodesic distance (angle) between two points on S³."""
    a = nrm4(a)
    b = nrm4(b)
    return float(np.arccos(np.clip(np.dot(a, b), -1.0, 1.0)))


def hsp(chi: float, th: float, ph: float) -> np.ndarray:
    """Hopf-sphere parameterisation: (χ, θ, φ) → unit 4-vector."""
    return np.array(
        [
            np.sin(chi) * np.sin(th) * np.cos(ph),
            np.sin(chi) * np.sin(th) * np.sin(ph),
            np.sin(chi) * np.cos(th),
            np.cos(chi),
        ],
        dtype=float,
    )


def s3_tangent_direction(src, dst) -> np.ndarray:
    """Unit tangent vector at *src* pointing toward *dst* on S³."""
    src = nrm4(src)
    dst = nrm4(dst)
    tang = dst - np.dot(dst, src) * src
    n = np.linalg.norm(tang)
    if n < 1e-12:
        e = np.array([1.0, 0.0, 0.0, 0.0])
        tang = e - np.dot(e, src) * src
        n = np.linalg.norm(tang)
    return tang / max(n, 1e-12)


def s3_green_potential(
    psi: float,
    radius: float = S3_RADIUS,
    eps: float = S3_GREEN_EPS,
) -> float:
    """Zero-mean Green function for the Laplacian on S³.

    G(ψ) = ((π − ψ) cot ψ − ½) / (4π² R)

    reproducing the Euclidean 1/(4πs) singularity near the source while
    remaining globally well-defined on compact S³.
    """
    psi_eff = float(np.clip(psi, eps, np.pi - eps))
    return float(
        (((np.pi - psi_eff) / np.tan(psi_eff)) - 0.5)
        / (4.0 * np.pi ** 2 * radius)
    )


def s3_green_field_kernel(
    psi: float,
    radius: float = S3_RADIUS,
    eps: float = S3_GREEN_EPS,
) -> float:
    """Magnitude of the geodesic electric field |dG/ds| on S³.

    For potential G(ψ) = ((π−ψ)cotψ − ½) / (4π²R), the derivative is

        dG/dψ = −[cotψ + (π−ψ)csc²ψ] / (4π²R)

    and the field magnitude is |E| = |dG/dψ| / R.
    """
    psi_eff = float(np.clip(psi, eps, np.pi - eps))
    sin_p = np.sin(psi_eff)
    cos_p = np.cos(psi_eff)
    cot_p = cos_p / max(sin_p, 1e-12)
    csc2_p = 1.0 / max(sin_p ** 2, 1e-12)
    dGdpsi = -(cot_p + (np.pi - psi_eff) * csc2_p) / (
        4.0 * np.pi ** 2 * radius
    )
    return float(abs(dGdpsi) / radius)
