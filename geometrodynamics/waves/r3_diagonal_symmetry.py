"""Diagonal half-clock and spatial actions, additive to frozen R3 maps.

Coordinates are (A,p_A,x1,p1,x2,p2), beta=x1 E0+x2 E1.
The eight-state diagonal reduction is checked against the full matrix system.
"""
from math import lcm
import numpy as np
from scipy.integrate import solve_ivp
from . import r3_family as rf
from . import nonlinear_supported_tt as full

BASIS = np.diagonal(rf.E[:2], axis1=1, axis2=2)
CYCLE = np.eye(3)[[1, 2, 0]]
REFLECTION = np.eye(3)[[1, 0, 2]]


def point(z):
    z = np.asarray(z, float)
    if z.shape != (6,) or not np.isfinite(z).all():
        raise ValueError('finite six-dimensional diagonal point required')
    return z


def spatial(z, matrix=CYCLE):
    z = point(z)
    representation = np.einsum('aij,bji->ab', rf.E[:2],
                              matrix @ rf.E[:2] @ matrix.T)
    out = z.copy()
    out[2::2] = representation @ z[2::2]
    out[3::2] = representation @ z[3::2]
    return out


def reduced_initial(z):
    z = point(z)
    A, pA = z[:2]
    x, v = z[2::2], z[3::2]/A**2
    m = np.exp(2*x @ BASIS)
    r = 2*(2*np.sum(1/m)-np.sum(m*m))
    qp2 = pA*pA/6-A*A*(v@v)+A*A*r-3*A**4
    if A <= 0 or qp2 <= 0:
        raise ValueError('nondegenerate positive chart and clock required')
    return np.r_[A, -pA/6, 0., -np.sqrt(qp2), x, v]


def reduced_rhs(t, y):
    A, Ap, q, qp = y[:4]
    x, v = y[4:6], y[6:8]
    h = A*A-q*q/6
    if A <= 0 or h <= 0:
        raise ValueError('left positive chart')
    hp = 2*A*Ap-q*qp/3
    m = np.exp(2*x @ BASIS)
    inv = 1/m
    r, ell = 2*(2*inv.sum()-(m*m).sum()), v@v
    force = BASIS @ ((-4+q*q/h)*inv-4*m*m)
    return np.r_[Ap, -A*(r+ell)/6+A**3, qp,
                 -(inv.sum()+(r+ell)/6)*q, v, -hp/h*v+force]


def clock_map(z, half=True, method='DOP853', matrix=False, tol=(1e-13, 1e-15)):
    """First noninitial crossing, with sign identification for a half return."""
    z = point(z)
    initial = rf.to_state(np.r_[z, np.zeros(6)]) if matrix else reduced_initial(z)
    clock_index = 3 if matrix else 2
    def event(t, y):
        return y[clock_index]
    event.direction = 1 if half else -1
    # Full return starts on a downward zero: ignore the initial event explicitly.
    event.terminal = half
    rhs = (lambda t, y: full.conformal_rhs(y)) if matrix else reduced_rhs
    sol = solve_ivp(rhs, (0., 2.2 if half else 3.95), initial,
                    events=event, method=method, rtol=tol[0], atol=tol[1], max_step=.05)
    hits = np.flatnonzero(sol.t_events[0] > 1e-8)
    if not sol.success or not len(hits):
        raise ArithmeticError('clock crossing failed: '+sol.message)
    k = hits[0]
    y = sol.y_events[0][k].copy()
    if half:
        if matrix:
            y[3:11] *= -1  # all four q and qprime components
        else:
            y[2:4] *= -1
    if matrix:
        out = rf.to_section(y)[:6]
        constraint = abs(full.constraints(y)['residual'][0])
    else:
        A, Ap, q, qp = y[:4]
        h = A*A-q*q/6
        out = np.r_[A, -6*Ap, y[4], h*y[6], y[5], h*y[7]]
        m = np.exp(2*y[4:6] @ BASIS)
        r = 2*(2*np.sum(1/m)-np.sum(m*m))
        constraint = abs(-3*Ap**2+qp**2/2+h*(y[6:8]@y[6:8])/2-h*r/2+q*q*np.sum(1/m)/2+1.5*A**4)
    return dict(z=out, time=float(sol.t_events[0][k]), constraint=float(constraint))


def angle(z):
    z = np.asarray(z)
    return np.arctan2(z[..., 4], z[..., 2])


def chirality(z):
    z = np.asarray(z)
    return z[..., 2]*z[..., 5]-z[..., 4]*z[..., 3]


def fourier_fit(nodes, degree=20):
    """Unsymmetrized fit in geometric shear azimuth, not a normal-form angle."""
    nodes = np.asarray(nodes, float)
    k = np.arange(-degree, degree+1)
    design = np.exp(1j*np.outer(angle(nodes), k))
    coefficients = np.linalg.lstsq(design, nodes, rcond=None)[0]
    def curve(theta):
        return (np.exp(1j*np.asarray(theta)[..., None]*k) @ coefficients).real
    return curve, float(np.linalg.cond(design))


def first_harmonic(half_order, spatial_order):
    """Conditional on faithful commuting circle rotations in a common angle."""
    if half_order < 1 or spatial_order < 1:
        raise ValueError('positive orders required')
    return lcm(half_order, spatial_order)
