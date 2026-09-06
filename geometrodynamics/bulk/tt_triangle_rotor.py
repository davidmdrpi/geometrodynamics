"""Test the triangle director inside the repository's existing ESU TT mode.

Public freeze 0eb684b. The ambient field is a five-component symmetric
trace-free tensor, not a rigid rotor. Restricting its action to a uniaxial
cone does not in general restrict its solutions to that cone.
"""

from dataclasses import dataclass
from functools import lru_cache
import math

import numpy as np


PUBLIC_PREREG = "0eb684b86a57be9cfbd2d41d89377cc11f6b76cf"
SEED = 20260907
IDENTITY = np.eye(3)


@dataclass(frozen=True)
class TensorModel:
    radius: float = 1.
    kappa: float = 1.

    def __post_init__(self):
        if not all(np.isfinite(v) and v > 0 for v in (self.radius, self.kappa)):
            raise ValueError("radius and gravitational coupling must be positive")

    @property
    def omega2(self):
        return 8. / self.radius ** 2

    @property
    def volume(self):
        return 2 * math.pi ** 2 * self.radius ** 3

    @property
    def normalization(self):
        """C = Vol(S3)/kappa in L = C/2 (tr Bdot^2 - omega^2 tr B^2)."""
        return self.volume / self.kappa


def _director(n):
    n = np.asarray(n, dtype=float)
    if n.shape != (3,) or not np.all(np.isfinite(n)) or abs(n @ n - 1) > 1e-10:
        raise ValueError("director must be a finite unit three-vector")
    return n


def stf(matrix):
    matrix = np.asarray(matrix, dtype=float)
    symmetric = (matrix + np.swapaxes(matrix, -1, -2)) / 2
    return symmetric - np.trace(symmetric, axis1=-2, axis2=-1)[..., None, None] * IDENTITY / 3


def axis_tensor(n):
    n = _director(n)
    return np.outer(n, n) - IDENTITY / 3


def embedding(amplitude, n):
    return float(amplitude) * axis_tensor(n)


def field_velocity(amplitude, amplitude_dot, n, v):
    n, v = _director(n), np.asarray(v, dtype=float)
    if v.shape != (3,) or not np.all(np.isfinite(v)) or abs(n @ v) > 1e-10:
        raise ValueError("director velocity must be finite and tangent")
    return amplitude_dot * axis_tensor(n) + amplitude * (np.outer(v, n) + np.outer(n, v))


def field_acceleration(amplitude, amplitude_dot, amplitude_ddot, n, v, n_ddot):
    return (amplitude_ddot * axis_tensor(n)
            + 2 * amplitude_dot * (np.outer(v, n) + np.outer(n, v))
            + amplitude * (np.outer(n_ddot, n) + np.outer(n, n_ddot) + 2 * np.outer(v, v)))


def normal_projection(matrix, n):
    """The two omitted biaxial components, orthogonal to the cone tangent."""
    n = _director(n)
    P = IDENTITY - np.outer(n, n)
    transverse = P @ np.asarray(matrix) @ P
    return transverse - np.trace(transverse) * P / 2


def tensor_components(matrix, n):
    """Orthogonal radial, director and biaxial projections of an STF tensor."""
    n = _director(n)
    matrix = np.asarray(matrix)
    radial = 1.5 * (n @ matrix @ n) * axis_tensor(n)
    v = (IDENTITY - np.outer(n, n)) @ matrix @ n
    angular = np.outer(v, n) + np.outer(n, v)
    return {"radial": radial, "angular": angular, "normal": normal_projection(matrix, n)}


def restricted_acceleration(amplitude, amplitude_dot, n, v, source=None,
                            model=TensorModel()):
    """Euler-Lagrange acceleration of the substituted action, NOT full flow."""
    if amplitude == 0:
        raise ValueError("amplitude-zero apex has no director chart")
    n, v = _director(n), np.asarray(v, dtype=float)
    source = np.zeros((3, 3)) if source is None else np.asarray(source)
    P = IDENTITY - np.outer(n, n)
    addot = 3 * amplitude * (v @ v) - model.omega2 * amplitude + 1.5 * (n @ source @ n)
    nddot = -(v @ v) * n - 2 * amplitude_dot * v / amplitude + P @ source @ n / amplitude
    return float(addot), nddot


def omitted_equation(amplitude, n, v, source=None):
    """Exact normal residual predicted by the embedding's second derivative."""
    n, v = _director(n), np.asarray(v)
    P = IDENTITY - np.outer(n, n)
    residual = 2 * amplitude * (np.outer(v, v) - (v @ v) * P / 2)
    if source is not None:
        residual = residual - normal_projection(source, n)
    return residual


def pullback_metric(amplitude, polar_angle, model=TensorModel()):
    """Coordinates (A, theta, phi) on the RP2 chart; singular at A=0."""
    C = model.normalization
    return np.diag([2 * C / 3, 2 * C * amplitude ** 2,
                    2 * C * amplitude ** 2 * math.sin(polar_angle) ** 2])


def normalized_source(integrated_stress, model=TensorModel()):
    """S=kappa/Vol * integral T_TT dVol; backreaction.shear_sources integrates."""
    return np.asarray(integrated_stress) / model.normalization


def lagrangian(beta, beta_dot, source=None, model=TensorModel()):
    source = np.zeros((3, 3)) if source is None else np.asarray(source)
    return float(model.normalization * (
        .5 * np.sum(beta_dot ** 2) - .5 * model.omega2 * np.sum(beta ** 2)
        + np.sum(beta * source)))


def restricted_lagrangian(amplitude, amplitude_dot, n, v, source=None,
                         model=TensorModel()):
    source = np.zeros((3, 3)) if source is None else np.asarray(source)
    return float(model.normalization * (amplitude_dot ** 2 / 3
        + amplitude ** 2 * (v @ v) - model.omega2 * amplitude ** 2 / 3
        + amplitude * (n @ source @ n)))


def exact_free_flow(beta0, velocity0, times, model=TensorModel()):
    times = np.asarray(times)
    omega = math.sqrt(model.omega2)
    cosine = np.cos(omega * times)[..., None, None]
    sine = np.sin(omega * times)[..., None, None]
    return (cosine * beta0 + sine * velocity0 / omega,
            -omega * sine * beta0 + cosine * velocity0)


def free_energy(beta, beta_dot, model=TensorModel()):
    return .5 * model.normalization * (np.sum(beta_dot ** 2, axis=(-2, -1))
                                      + model.omega2 * np.sum(beta ** 2, axis=(-2, -1)))


def rotational_charge(beta, beta_dot, model=TensorModel()):
    """Conserved antisymmetric matrix C [beta, beta_dot]."""
    return model.normalization * (beta @ beta_dot - beta_dot @ beta)


def nearest_uniaxial(beta):
    """Global Frobenius minimum over all signed amplitudes and directors.

    For a fixed n, A*=3 n.B.n/2. The best n is an eigenvector with largest
    absolute eigenvalue. The apex/ties have no unique director; no physical
    axis is assigned at the zero field.
    """
    beta = np.asarray(beta, dtype=float)
    if beta.shape != (3, 3) or not np.all(np.isfinite(beta)):
        raise ValueError("expected a finite three by three tensor")
    if np.linalg.norm(beta - stf(beta)) > 1e-10:
        raise ValueError("tensor must be symmetric and trace-free")
    values, vectors = np.linalg.eigh(beta)
    index = int(np.argmax(np.abs(values)))
    amplitude = 1.5 * values[index]
    closest = embedding(amplitude, vectors[:, index])
    # Independent eigenvalue-distance identity avoids subtracting large squares.
    pair = np.delete(values, index)
    return {"distance": float(np.linalg.norm(beta - closest)),
            "eigenvalue_distance": float(abs(pair[0] - pair[1]) / math.sqrt(2)),
            "amplitude": float(amplitude), "closest": closest,
            "axis_defined": bool(np.linalg.norm(beta) > 0
                                 and np.sum(np.isclose(abs(values), abs(values[index]),
                                                      rtol=1e-12, atol=0)) == 1)}


def uniform_rotation(time, amplitude=.01, speed=.4, model=TensorModel()):
    """A manufactured rotating director and the source it REQUIRES."""
    n = np.array([math.sin(speed * time), 0., math.cos(speed * time)])
    v = speed * np.array([math.cos(speed * time), 0., -math.sin(speed * time)])
    acceleration = -speed ** 2 * n
    beta = embedding(amplitude, n)
    bdot = field_velocity(amplitude, 0., n, v)
    bddot = field_acceleration(amplitude, 0., 0., n, v, acceleration)
    return {"n": n, "v": v, "beta": beta, "beta_dot": bdot,
            "beta_ddot": bddot, "source": bddot + model.omega2 * beta}


def closed_form_free_spectrum(times, amplitude=.01, speed=.4, model=TensorModel()):
    """Post-review exact spectrum/distance for the frozen A_dot0=0 family.

    Coordinates are (n0, vhat, n0 cross vhat). This is the entire free
    trajectory, including periodic returns, not a small-time approximation.
    """
    times = np.asarray(times, dtype=float)
    omega = math.sqrt(model.omega2)
    c, s = np.cos(omega * times), speed * np.sin(omega * times) / omega
    radius = np.hypot(c / 2, s)
    values = amplitude * np.stack((c / 6 - radius, c / 6 + radius, -c / 3), axis=-1)
    denominator = radius + abs(c) / 2
    # Rationalize the small gap to avoid cancellation near a return.
    gap = np.divide(s * s, denominator, out=np.zeros_like(s), where=denominator > 0)
    return {"eigenvalues": np.sort(values, axis=-1),
            "distance": abs(amplitude) * gap / math.sqrt(2), "c": c, "s": s}


def continuous_eigenline(times, speed=.4, model=TensorModel()):
    """A continuously advancing eigenline of the biaxial field, not a rotor.

    n0=e_z, vhat=e_x, A0>0. At repeated eigenvalues this continuation is
    a choice; the field alone does not distinguish a line in that eigenspace.
    The nearest-uniaxial axis instead switches branch at cos(omega t)=0.
    """
    if not np.isfinite(speed) or speed <= 0:
        raise ValueError("positive director speed required for this branch")
    phase = math.sqrt(model.omega2) * np.asarray(times, dtype=float)
    winding = np.floor((phase + math.pi) / (2 * math.pi))
    reduced = phase - 2 * math.pi * winding
    c = np.cos(reduced)
    s = speed * np.sin(reduced) / math.sqrt(model.omega2)
    angle = .5 * np.arctan2(2 * s, c) + math.pi * winding
    axes = np.stack((np.sin(angle), np.zeros_like(angle), np.cos(angle)), axis=-1)
    return {"angle": angle, "axis": axes, "angular_speed": speed / (c * c + 4 * s * s)}


@lru_cache(maxsize=1)
def adm_quadratic_derivation():
    """Independent spatial Koszul curvature and ADM kinetic/action expansion.

    Uses [e_i,e_j]=2 epsilon_ijk a_k/(a_i a_j) e_k on SU(2), not the
    repository's four-dimensional Cartan/Einstein routine or its frequency.
    Isotropy extends the diagonal quadratic coefficient to the STF space.
    """
    import sympy as sp
    lengths = sp.symbols("l1 l2 l3", positive=True)
    c = sp.MutableDenseNDimArray.zeros(3, 3, 3)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                c[i, j, k] = 2 * sp.LeviCivita(i, j, k) * lengths[k] / (lengths[i] * lengths[j])
    connection = sp.MutableDenseNDimArray.zeros(3, 3, 3)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                connection[i, j, k] = (c[i, j, k] - c[j, k, i] + c[k, i, j]) / 2
    scalar = sp.simplify(sum(
        connection[j, j, m] * connection[i, m, i]
        - connection[i, j, m] * connection[j, m, i]
        - c[i, j, m] * connection[m, j, i]
        for i in range(3) for j in range(3) for m in range(3)))
    a, eps = sp.symbols("a epsilon", positive=True)
    b1, b2 = sp.symbols("b1 b2", real=True)
    b = (b1, b2, -b1 - b2)
    anisotropic = scalar.subs({lengths[i]: a * sp.exp(eps * b[i]) for i in range(3)})
    round_value = sp.simplify(anisotropic.subs(eps, 0))
    linear = sp.simplify(sp.diff(anisotropic, eps).subs(eps, 0))
    quadratic = sp.simplify(sp.diff(anisotropic, eps, 2).subs(eps, 0) / 2)
    norm2 = sum(v ** 2 for v in b)
    coefficient = sp.simplify(quadratic / norm2)
    # N=1, fixed a, sum b_i=0: K^i_j=epsilon bdot_i delta^i_j.
    d1, d2 = sp.symbols("d1 d2", real=True)
    rates = (d1, d2, -d1 - d2)
    kinetic = sp.expand(sum(d ** 2 for d in rates) - sum(rates) ** 2)
    return {"spatial_scalar": str(scalar), "round_scalar": str(round_value),
            "linear_scalar": str(linear), "curvature_quadratic": str(quadratic),
            "curvature_coefficient_per_tr_beta2": str(coefficient),
            "adm_kinetic_quadratic": str(kinetic),
            "omega2_radius2": float(sp.simplify(-coefficient * a ** 2)),
            "round_anchor": bool(sp.simplify(round_value - 6 / a ** 2) == 0),
            "linear_vanishes": bool(linear == 0),
            "kinetic_is_frobenius": bool(sp.simplify(kinetic - sum(d ** 2 for d in rates)) == 0)}


def physical_verdict(checks, analytic_obstruction):
    if not checks or not all(checks.values()) or not analytic_obstruction:
        return "UNRESOLVED"
    return "FREE_ROTATING_UNIAXIAL_TT_FAMILY_NOT_INVARIANT"
