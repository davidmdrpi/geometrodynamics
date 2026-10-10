"""Clock-sign quotient of the LRS section map; see r3_half_clock_prereg.md.

The exact involution S reverses (q,q') and leaves the geometry unchanged.
H = S o flow_to_next_upward_zero acts on the downward clock section. P=H^2.
This module uses the scalar reduced equations; the existing esu_map uses
the full homogeneous matrix equations and is the independent check.
"""
from fractions import Fraction

import numpy as np
from scipy.integrate import solve_ivp

from . import r3_return_map as rm


def clock_sign(state):
    out = np.asarray(state, float).copy()
    out[2:4] *= -1
    return out


def lifted_rotation(p, q):
    """Branch established on the existing LRS circles: rho_H=(rho_P+1)/2."""
    return (Fraction(p, q) + 1) / 2


def half_map(z, method='DOP853', tol=(1e-13, 1e-15), full=False):
    """Next upward q=0 crossing, identified by S with the downward section.

    Section coordinates do not contain q or q', so applying S is implicit
    in the returned coordinates. It is explicit in the returned full state.
    """
    z = np.asarray(z, float)
    if z.shape != (4,) or not np.isfinite(z).all():
        raise ValueError('finite four-dimensional section point required')
    initial = np.asarray(rm.section_to_state(z), float)
    if not np.isfinite(initial).all() or initial[3] >= 0:
        raise ValueError('nondegenerate downward clock required')

    def crossing(t, y):
        return y[2]
    crossing.direction = 1
    crossing.terminal = True
    sol = solve_ivp(lambda t, y: rm.rhs(y), (0., 2.2), initial,
                    method=method, rtol=tol[0], atol=tol[1], events=crossing)
    if not sol.success or len(sol.t_events[0]) != 1:
        raise ArithmeticError('half-clock crossing unresolved: '+sol.message)
    y = clock_sign(sol.y_events[0][0])
    if y[3] >= 0 or not np.isfinite(y).all():
        raise ArithmeticError('invalid half-clock endpoint')
    result = np.asarray(rm.state_to_section(y), float)
    if full:
        return dict(z=result, state=y, time=float(sol.t_events[0][0]),
                    constraint=float(abs(rm.constraint(y))))
    return result


def squared_map(z):
    return half_map(half_map(z))
