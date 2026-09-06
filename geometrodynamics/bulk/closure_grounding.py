"""Conditional classical parents and prior controls, not a BAM derivation."""

import math

import numpy as np
from scipy.integrate import quad
from scipy.linalg import expm, solve_banded

from geometrodynamics.waves.finite_throat import dtn_matrix


def tube_frame_energy(rotation, length=1.3, area=.7, n_interior=31):
    """Minimize nine scalar gradient energies with endpoint matrices I,R.

    The interior matrix is UNCONSTRAINED, not restricted to SO(3). Vector
    channels and endpoint/holonomy matching are additions to the scalar
    monopole implementation. The finite-k DtN call only checks its static
    limit; static elimination is not a dynamical or causal construction.
    """
    R = np.asarray(rotation, dtype=float)
    if (R.shape != (3, 3) or not np.allclose(R.T @ R, np.eye(3), atol=1e-10)
            or not np.isclose(np.linalg.det(R), 1., atol=1e-10)):
        raise ValueError("endpoint rotation must be in SO(3)")
    if length <= 0 or area <= 0 or n_interior < 1:
        raise ValueError("positive tube dimensions and interior size required")
    ab = np.zeros((3, n_interior))
    ab[0, 1:], ab[1, :], ab[2, :-1] = -1., 2., -1.
    rhs = np.zeros((n_interior, 9))
    rhs[0] += np.eye(3).ravel()
    rhs[-1] += R.ravel()
    interior = solve_banded((1, 1), ab, rhs)
    field = np.vstack([np.eye(3).ravel(), interior, R.ravel()])
    dx = length / (n_interior + 1)
    fd = area / (2 * dx) * np.sum(np.diff(field, axis=0) ** 2)
    ends = np.vstack([np.eye(3).ravel(), R.ravel()])
    # Existing dtn_matrix has a removable 0/0 at k=0. Approach, do not replace it.
    dtn = dtn_matrix(1e-6 / length, length, area).real
    dtn_energy = .5 * np.sum(ends * (dtn @ ends))
    predicted = area / (2 * length) * np.linalg.norm(R - np.eye(3)) ** 2
    return {"finite_difference_energy": float(fd), "dtn_energy": float(dtn_energy),
            "predicted_energy": float(predicted), "K": 8 * area / length,
            "fd_error": float(abs(fd - predicted)),
            "dtn_error": float(abs(dtn_energy - predicted))}


def momentum_integrals():
    """Distinct radial momentum laws giving the same uniform position law."""
    gaussian = 2 * math.pi * quad(lambda r: r * math.exp(-r * r / 2), 0, np.inf)[0]
    quartic = 2 * math.pi * quad(lambda r: r * math.exp(-r ** 4 / 4), 0, np.inf)[0]
    return {"gaussian": gaussian, "quartic": quartic}


def inertia_controls(x, epsilon=.5):
    """Smooth metrics on S2: variable volume vs anisotropy with unit volume.

    Embedded normal eigenvalue is one, so the 3D determinant equals the
    intrinsic tangent determinant. The anisotropic field is smooth at poles.
    """
    x = np.asarray(x, dtype=float)
    x = x / np.linalg.norm(x)
    if abs(epsilon) >= 1:
        raise ValueError("abs(epsilon) must be below one")
    normal = np.outer(x, x)
    tangent = np.eye(3) - normal
    inertia = 1 + epsilon * x[2]
    variable = inertia * tangent + normal
    v = tangent @ np.array([0., 0., 1.])
    traceless = np.outer(v, v) - .5 * (v @ v) * tangent
    anisotropic = expm(epsilon * traceless)
    return {"variable_volume": float(np.sqrt(np.linalg.det(variable))),
            "expected_variable_volume": float(inertia),
            "anisotropic_volume": float(np.sqrt(np.linalg.det(anisotropic))),
            "anisotropy": float(np.linalg.norm(anisotropic - np.eye(3)))}
