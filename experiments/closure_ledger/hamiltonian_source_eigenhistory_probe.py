"""
The Hamiltonian source: the imposed phase law replaced by a minimal
conservative source model, and the eigenhistory solved as
U(X) X = X together with total-energy closure (PR #219).

> Framing: QFT on the *fixed classical* throat geometry (geometry -> fields),
> not quantum gravity.

WHAT #218 IMPOSED
-----------------
#218 established the eigenhistory existence theorem with a source whose
phase law phi_s(I) = beta I/(1 + I/I_sat) was IMPOSED - a reasonable
caricature of internal-state pulling, but not derived.  This probe
replaces it with the minimal conservative source that has an explicit
Hamiltonian:

    H = p^2/2 + w0^2 q^2/2 + (mu/4) q^4 + g q u(0)

- a side-coupled Duffing oscillator (the quartic is the minimal
anharmonicity; the cubic is absent by q -> -q symmetry), linearly
coupled to the loop field at the crossing.  Harmonic balance at
frequency w (q = Re[a e^{-iwt}], u(0) = Re[U0 e^{-iwt}], gauge U0 real):

    D(a) a = -g U0 ,     D(a) = w0^2 - w^2 + (3/4) mu a^2 ,

and the field jump [u'] = g a makes the source a point scatterer of
DERIVED strength kappa_eff = g a / U0 (-> -g^2/D linearly):

    t_s = 2iw/(2iw - kappa) ,   r_s = kappa/(2iw - kappa) ,

UNITARY BY CONSTRUCTION (|t|^2 + |r|^2 = 1 for real kappa: the
Hamiltonian is conservative, so the scattering must be - nothing
imposed), REACTIVE (zero net power: a and U0 in phase on the real
branch), with the amplitude-dependent phase now DERIVED - and reducing
to the #218 form in the weak-coupling limit (phi -> -kappa/2w =
g^2/(2 w D(I)), the amplitude dependence carried by D).

THE RING EIGENPROBLEM WITH THE SOURCE INSIDE.  Because the Hamiltonian
source also REFLECTS, the loop is a genuine two-direction ring:
[source] -pi- [mouth A] -tau- [mouth B] -pi- [source], with the state
X = (right-mover, left-mover; source amplitude a).  The homogeneous
condition U(X) X = X is the ring monodromy eigenproblem with the
source's S-matrix depending on the state itself.  Its solutions form a
BRANCH w*(A) (the nonlinear ring modes at the gap edges - and the fine
scan resolves the tr > 2 gap segment that #217's coarse ring scan
missed), so the amplitude is NOT fixed by the homogeneous condition
alone: TOTAL-ENERGY CLOSURE fixes it -

    tr T_ring(w, A) = 2   AND   E_field + E_source + <g q u(0)> = E_0

solved jointly by 2D Newton - the ledger CORRECTED to include the
time-averaged interaction energy <g q u(0)> = g a U0/2 (negative:
a and u(0) antiphased).  Reported, as requested: the fixed-point
residual (|U(X)X - X| ~ 1e-14) and the FULL HAMILTONIAN STABILITY
SPECTRUM - the Duffing (q, p) evolved as INDEPENDENT variables in the
4x4 variational monodromy of the reduced 2-dof Hamiltonian about its
shooting-refined periodic orbit, not the algebraically slaved
harmonic-balance spectrum: the Floquet-trivial pair at 1 plus the
SOURCE pair rotating at its dressed frequency (3.210 vs bare 3.2),
all on the unit circle, symplectic to machine precision.

Tests:
  T1. Goal.
  T2. The Hamiltonian source, derived (unitary by construction;
      reactive; the cubic branch; the #218 law as the weak-coupling
      limit).
  T3. The ring with the source (real trace; the resolved gap; the
      nonlinear mode branch w*(A)).
  T4. The joint solve U(X)X = X + total-energy closure (2D Newton;
      the fixed-point residuals; the state including the source; the
      energy partition).
  T5. Energy closure element-wise and dynamically (flux constant
      around the ring; zero source power; 1e4-pass persistence).
  T6. The FULL Hamiltonian stability spectrum ((q, p) independent in
      the variational monodromy; the source pair at its dressed
      frequency; symplectic; the slaved HB spectrum retained for
      comparison).
  T7. Honest scope.
  T8. Assessment.

Verdict:
  THE_SOURCE_PHASE_LAW_IS_DERIVED_FROM_A_HAMILTONIAN_THE_EIGENHISTORY_
  SOLVES_UXX_EQUALS_X_WITH_TOTAL_ENERGY_CLOSURE_AND_ITS_STABILITY_
  SPECTRUM_IS_MARGINAL
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

_CACHE: dict = {}
_RH = 1.0
_W0 = 3.2                # source natural frequency
_MU = 0.5                # Duffing quartic (minimal anharmonicity)
_G = 0.6                 # linear coupling
_TAU = 0.8               # throat interior transit
_L1 = math.pi            # source -> mouth A
_L2 = math.pi            # mouth B -> source

# ========================================================================
# SECTION A - the two-sided greybody and the unitarized port tables
# ========================================================================


def _v_of_r(r):
    f = 1.0 - (_RH / r) ** 2
    return f * (0.75 / r ** 2 + 2.25 * _RH ** 2 / r ** 4)


def _x_of_r(r):
    return r + (_RH / 2) * np.log((r - _RH) / (r + _RH))


def _rhs(w):
    def rhs(x, y):
        pr, pi_, qr, qi, r = y
        vv = _v_of_r(r)
        f = 1 - (_RH / r) ** 2
        return [qr, qi, (vv - w ** 2) * pr, (vv - w ** 2) * pi_, f]
    return rhs


def _solve_and_match(w, x0, x_far, y0, span, basis_sign, rtol=1e-9):
    sol = solve_ivp(_rhs(w), span, y0, rtol=rtol, atol=1e-12,
                    dense_output=True, method="DOP853")
    if basis_sign > 0:
        xs = x_far - np.linspace(0.0, 3.7, 8)
        basis = lambda xm: [np.exp(-1j * w * xm), np.exp(1j * w * xm)]
    else:
        xs = x0 + np.linspace(0.0, 2.0, 8)
        basis = lambda xm: [np.exp(1j * w * xm), np.exp(-1j * w * xm)]
    rows, vals = [], []
    for xm in xs:
        pr, pi_, *_ = sol.sol(xm)
        rows.append(basis(xm))
        vals.append(pr + 1j * pi_)
    (alpha, beta), *_ = np.linalg.lstsq(np.array(rows), np.array(vals),
                                        rcond=None)
    return complex(1.0 / alpha), complex(beta / alpha)


def _greybody(w):
    key = ('g', w)
    if key in _CACHE:
        return _CACHE[key]
    r0 = _RH * (1 + 1e-7)
    x0 = float(_x_of_r(r0))
    r_far = max(60.0 / w, 50.0 * w, 30.0)
    x_far = float(_x_of_r(r_far))
    y0 = [math.cos(-w * x0), math.sin(-w * x0),
          w * math.sin(-w * x0), -w * math.cos(-w * x0), r0]
    t, r_out = _solve_and_match(w, x0, x_far, y0, (x0, x_far), +1)
    p0 = np.exp(1j * w * x_far)
    q0 = 1j * w * p0
    y1 = [p0.real, p0.imag, q0.real, q0.imag, r_far]
    _, r_in = _solve_and_match(w, x0, x_far, y1, (x_far, x0), -1)
    out = (t, r_out, r_in)
    _CACHE[key] = out
    return out


def _exact_port_raw(w):
    t, r_out, r_in = _greybody(w)
    n = math.sqrt(abs(t) ** 2 + abs(r_out) ** 2)
    t_u = t / n
    r_u = math.sqrt(max(1 - abs(t_u) ** 2, 0.0)) * np.exp(
        1j * np.angle(r_in))
    return complex(t_u), complex(r_u)


_PGRID = np.linspace(2.4, 3.8, 57)


def _splines():
    if 'SPL' not in _CACHE:
        vals = [_exact_port_raw(float(w)) for w in _PGRID]
        tv = np.array([v[0] for v in vals])
        rv = np.array([v[1] for v in vals])
        _CACHE['SPL'] = (CubicSpline(_PGRID, tv.real),
                         CubicSpline(_PGRID, tv.imag),
                         CubicSpline(_PGRID, rv.real),
                         CubicSpline(_PGRID, rv.imag))
    return _CACHE['SPL']


def exact_port(w):
    s = _splines()
    t_u = complex(s[0](w), s[1](w))
    r_u = complex(s[2](w), s[3](w))
    n = math.sqrt(abs(t_u) ** 2 + abs(r_u) ** 2)
    return t_u / n, r_u / n          # re-unitarized interpolant


def barrier_M(w, rev=False):
    t_u, r_u = exact_port(w)
    r_out_u = -np.conj(r_u) * t_u / np.conj(t_u)
    ro, ri = (r_u, r_out_u) if rev else (r_out_u, r_u)
    return np.array([[(t_u * t_u - ro * ri) / t_u, ri / t_u],
                     [-ro / t_u, 1 / t_u]])


# ── the Hamiltonian source ──────────────────────────────────────────────

def source_a(w, U0):
    """The continuous real Duffing branch of D(a) a = -g U0."""
    if abs(U0) < 1e-30:
        return 0.0
    roots = np.roots([0.75 * _MU, 0.0, _W0 ** 2 - w ** 2, _G * U0])
    real = [r.real for r in roots if abs(r.imag) < 1e-9]
    lin = -_G * U0 / (_W0 ** 2 - w ** 2)
    return float(min(real, key=lambda r: abs(r - lin)))


def kappa_eff(w, U0):
    if abs(U0) < 1e-30:
        return -_G ** 2 / (_W0 ** 2 - w ** 2)
    return _G * source_a(w, U0) / U0


def source_scattering(w, U0):
    k = kappa_eff(w, U0)
    t_s = 2j * w / (2j * w - k)
    r_s = k / (2j * w - k)
    return complex(t_s), complex(r_s)


def source_M(w, U0):
    t_s, r_s = source_scattering(w, U0)
    return np.array([[(t_s * t_s - r_s * r_s) / t_s, r_s / t_s],
                     [-r_s / t_s, 1 / t_s]])


def P(w, L):
    return np.diag([np.exp(1j * w * L), np.exp(-1j * w * L)])


def ring_T(w, U0):
    return (P(w, _L2) @ barrier_M(w, rev=True) @ P(w, _TAU)
            @ barrier_M(w) @ P(w, _L1) @ source_M(w, U0))


def tr_ring(w, U0):
    return float(complex(np.trace(ring_T(w, U0))).real)


def eigvec(w, U0):
    T = ring_T(w, U0)
    _, s, vh = np.linalg.svd(T - np.eye(2))
    return vh[-1].conj(), float(s[-1])


def energies(w, U0):
    """Total energy of the U0-normalized eigenhistory: field on the
    three segments + the time-averaged Duffing energy + the
    time-averaged INTERACTION energy <g q u(0)> = g a U0 / 2 (both
    real in the gauge; negative here since a and U0 are antiphased) -
    the full Hamiltonian ledger."""
    v, resid = eigvec(w, U0)
    segs = [( _L1, v)]
    x = barrier_M(w) @ (P(w, _L1) @ v)
    segs.append((_TAU, x))
    x = barrier_M(w, rev=True) @ (P(w, _TAU) @ x)
    segs.append((_L2, x))
    scale = U0 / abs(v[0] + v[1])
    E_f = sum(L * w ** 2 * (abs(s0[0]) ** 2 + abs(s0[1]) ** 2) / 2
              for L, s0 in segs) * scale ** 2
    a = source_a(w, U0)
    E_s = (w ** 2 * a ** 2 / 4 + _W0 ** 2 * a ** 2 / 4
           + (3 / 32) * _MU * a ** 4)
    E_int = 0.5 * _G * a * U0
    return float(E_f), float(E_s), float(E_int), resid


def joint_solve(E0, w_init, U_init):
    """2D Newton: tr T_ring = 2 (homogeneous condition) together with
    the CORRECTED total-energy closure
    E_field + E_source + <g q u(0)> = E0."""
    x = np.array([w_init, U_init])
    F = np.array([1.0, 1.0])
    for _ in range(60):
        w, U0 = x
        Ef, Es, Ei, _ = energies(w, U0)
        F = np.array([tr_ring(w, U0) - 2, Ef + Es + Ei - E0])
        if np.abs(F).max() < 1e-12:
            break
        J = np.zeros((2, 2))
        h = 1e-6
        for j, dx in enumerate(([h, 0], [0, h])):
            wp, Up = x + np.array(dx)
            Efp, Esp, Eip, _ = energies(wp, Up)
            J[0, j] = (tr_ring(wp, Up) - 2 - F[0]) / h
            J[1, j] = (Efp + Esp + Eip - E0 - F[1]) / h
        x = x - np.linalg.solve(J, F)
    return x, F


def scaled_state(w, U0):
    """The physical eigenstate: (v scaled and gauged, a)."""
    v, _ = eigvec(w, U0)
    v = v * (U0 / abs(v[0] + v[1]))
    v = v * np.exp(-1j * np.angle(v[0] + v[1]))
    return v, source_a(w, abs(v[0] + v[1]))


# ========================================================================
# T1. Goal
# ========================================================================


def test_T1_goal() -> dict:
    return {
        'name': 'T1_goal',
        'description': (
            'Replace #218\'s imposed source phase law with a minimal '
            'conservative source derived from an explicit Hamiltonian '
            '(side-coupled Duffing oscillator, harmonic balance); '
            'include the source state and energy in the loop; solve '
            'the homogeneous condition U(X)X = X together with '
            'total-energy closure; report the fixed-point residual '
            'and the stability eigenvalues.'
        ),
        'pass': True,
    }


# ========================================================================
# T2. The Hamiltonian source, derived
# ========================================================================


def test_T2_source_derived() -> dict:
    # (i) unitarity is a CONSEQUENCE (conservative H), not an input
    unit_dev, power_dev = 0.0, 0.0
    for w, U0 in ((2.8, 0.5), (3.0, 2.0), (3.4, 1.0), (2.73, 0.9)):
        t_s, r_s = source_scattering(w, U0)
        unit_dev = max(unit_dev, abs(abs(t_s) ** 2 + abs(r_s) ** 2 - 1))
        a = source_a(w, U0)
        power_dev = max(power_dev,
                        abs(0.5 * _G * w * np.imag(a * np.conj(U0))))

    # (ii) the cubic branch is continuous from the linear response
    U_grid = np.linspace(1e-4, 2.0, 100)
    a_grid = [source_a(2.9, U) for U in U_grid]
    jumps = max(abs(a_grid[i + 1] - a_grid[i]) for i in range(99))

    # (iii) the #218 imposed law is the weak-coupling limit:
    # phi_s -> -kappa/(2w) = g^2/(2 w D(I)), amplitude dependence in D
    w = 2.8
    lim_ratio = []
    for gg in (0.05, 0.02):
        k = -gg ** 2 / (_W0 ** 2 - w ** 2)
        t_s = 2j * w / (2j * w - k)
        lim_ratio.append(float(np.angle(t_s) / (-k / (2 * w))))

    ok = (unit_dev < 1e-12 and power_dev < 1e-14 and jumps < 0.05
          and all(abs(r - 1) < 1e-4 for r in lim_ratio))
    return {
        'name': 'T2_source_derived',
        'description': (
            'H = p^2/2 + w0^2 q^2/2 + mu q^4/4 + g q u(0): the '
            'scattering is unitary BECAUSE the Hamiltonian is '
            'conservative, reactive (zero net power) on the real '
            'branch, and reduces to the #218 law at weak coupling'
        ),
        'unitarity_deviation': float(unit_dev),
        'net_power_deviation': float(power_dev),
        'branch_max_jump': float(jumps),
        'weak_coupling_limit_ratios': lim_ratio,
        'pass': bool(ok),
    }


# ========================================================================
# T3. The ring with the Hamiltonian source
# ========================================================================


def test_T3_ring_modes() -> dict:
    if 'T3' in _CACHE:
        return _CACHE['T3']
    # (i) time-reversal: the ring trace is real
    tr_imag = max(abs(complex(np.trace(ring_T(w, 0.8))).imag)
                  for w in (2.7, 3.0, 3.5))

    # (ii) the resolved gap: tr > 2 on a finite segment (the split
    # mode pair #217's coarse scan reported as a tangency)
    ws = np.linspace(2.70, 2.78, 801)
    trs = np.array([tr_ring(float(w), 0.8) for w in ws])
    above = trs > 2
    gap_found = bool(above.any())
    gap_lo = float(ws[above][0]) if gap_found else None
    gap_hi = float(ws[above][-1]) if gap_found else None
    gap_max = float(trs.max())

    # (iii) the nonlinear mode branch w*(A): the homogeneous condition
    # alone defines a CURVE, not a point
    branch = []
    for U0 in (0.2, 0.5, 0.8, 1.2, 1.6):
        wlo, whi = 2.70, 2.78
        wg = np.linspace(wlo, whi, 200)
        tg = [tr_ring(float(w), U0) for w in wg]
        wstar = None
        for i in range(len(wg) - 1):
            if (tg[i] - 2) * (tg[i + 1] - 2) < 0:
                wstar = brentq(lambda w: tr_ring(w, U0) - 2,
                               wg[i], wg[i + 1], xtol=1e-11)
                break
        _, _, _, resid = energies(wstar, U0)
        branch.append({'U0': U0, 'w_star': float(wstar),
                       'eig_residual': resid})
    w_span = branch[-1]['w_star'] - branch[0]['w_star']

    ok = (tr_imag < 1e-10 and gap_found and gap_max > 2
          and all(b['eig_residual'] < 1e-12 for b in branch)
          and abs(w_span) > 1e-6)
    out = {
        'name': 'T3_ring_modes',
        'description': (
            'the two-direction ring with the source inside: real '
            'trace (time reversal); the tr > 2 gap segment resolved '
            '(the #217 tangency was a coarse-grid artifact); the '
            'homogeneous condition defines a nonlinear mode BRANCH '
            'w*(A) - the amplitude is not fixed by homogeneity alone'
        ),
        'trace_imag_max': float(tr_imag),
        'gap_segment': [gap_lo, gap_hi],
        'gap_max_trace': gap_max,
        'mode_branch': branch,
        'branch_frequency_span': float(w_span),
        'pass': bool(ok),
    }
    _CACHE['T3'] = out
    return out


# ========================================================================
# T4. The joint solve: U(X)X = X + total-energy closure
# ========================================================================


def test_T4_joint_solve() -> dict:
    if 'T4' in _CACHE:
        return _CACHE['T4']
    E0 = 11.177
    x, F = joint_solve(E0, 2.7324, 0.9)
    w_s, U_s = float(x[0]), float(x[1])
    Ef, Es, Ei, eig_resid = energies(w_s, U_s)
    v, a_s = scaled_state(w_s, U_s)
    T = ring_T(w_s, U_s)
    ux_resid = float(np.linalg.norm(T @ v - v) / np.linalg.norm(v))
    # the source's slaving equation residual
    D = _W0 ** 2 - w_s ** 2 + 0.75 * _MU * a_s ** 2
    slave_resid = abs(D * a_s + _G * abs(v[0] + v[1]))

    # the interpolated-port systematic: re-evaluate the homogeneous
    # condition with RAW (non-interpolated) ports at the fixed point
    def raw_tr(w, U0):
        t_u, r_u = _exact_port_raw(w)
        r_o = -np.conj(r_u) * t_u / np.conj(t_u)
        MA = np.array([[(t_u * t_u - r_o * r_u) / t_u, r_u / t_u],
                       [-r_o / t_u, 1 / t_u]])
        MB = np.array([[(t_u * t_u - r_u * r_o) / t_u, r_o / t_u],
                       [-r_u / t_u, 1 / t_u]])
        T2 = (P(w, _L2) @ MB @ P(w, _TAU) @ MA @ P(w, _L1)
              @ source_M(w, U0))
        return float(complex(np.trace(T2)).real)
    interp_syst = abs(raw_tr(w_s, U_s) - 2.0)

    ok = (np.abs(F).max() < 1e-11 and ux_resid < 1e-12
          and eig_resid < 1e-12 and slave_resid < 1e-10
          and abs(Ef + Es + Ei - E0) < 1e-9
          and Ei < 0
          and interp_syst < 1e-3)
    out = {
        'name': 'T4_joint_solve',
        'description': (
            'the homogeneous condition solved TOGETHER with the '
            'CORRECTED total-energy closure (E_field + E_source + '
            '<g q u(0)> = E0): 2D Newton; the state X* includes the '
            'source amplitude; residuals reported as requested'
        ),
        'E0': E0,
        'w_star': w_s,
        'U0_star': U_s,
        'a_star': float(a_s),
        'newton_residuals': [float(F[0]), float(F[1])],
        'UX_equals_X_residual': ux_resid,
        'eigen_residual': eig_resid,
        'source_slaving_residual': float(slave_resid),
        'E_field': Ef,
        'E_source': Es,
        'E_interaction': Ei,
        'energy_closure_residual': float(abs(Ef + Es + Ei - E0)),
        'raw_port_systematic': float(interp_syst),
        'pass': bool(ok),
    }
    _CACHE['T4'] = out
    return out


# ========================================================================
# T5. Energy closure, element-wise and dynamically
# ========================================================================


def test_T5_energy_closure() -> dict:
    t4 = test_T4_joint_solve()
    w_s, U_s = t4['w_star'], t4['U0_star']
    v, a_s = scaled_state(w_s, U_s)

    # (i) the flux form |R|^2 - |L|^2 is constant around the ring
    def flux(s):
        return abs(s[0]) ** 2 - abs(s[1]) ** 2
    f1 = flux(v)
    x = barrier_M(w_s) @ (P(w_s, _L1) @ v)
    f2 = flux(x)
    x = barrier_M(w_s, rev=True) @ (P(w_s, _TAU) @ x)
    f3 = flux(x)
    flux_dev = max(abs(f2 - f1), abs(f3 - f1))

    # (ii) zero net power into the source at the fixed point
    p_src = abs(0.5 * _G * w_s * np.imag(a_s * np.conj(v[0] + v[1])))

    # (iii) 1e4-pass persistence of the full nonlinear state
    z = v.copy()
    for _ in range(10000):
        U0 = abs(z[0] + z[1])
        z = ring_T(w_s, U0) @ z
    norm_drift = abs(np.linalg.norm(z) / np.linalg.norm(v) - 1)
    U_drift = abs(abs(z[0] + z[1]) - U_s) / U_s

    ok = (flux_dev < 1e-12 and p_src < 1e-12
          and norm_drift < 1e-8 and U_drift < 1e-8)
    return {
        'name': 'T5_energy_closure',
        'description': (
            'the conserved flux form is constant around the ring; the '
            'source absorbs zero net power; the full nonlinear state '
            'persists 1e4 passes'
        ),
        'flux_constancy': float(flux_dev),
        'source_net_power': float(p_src),
        'norm_drift_1e4': float(norm_drift),
        'U0_drift_1e4': float(U_drift),
        'pass': bool(ok),
    }


# ── the reduced 2-dof Hamiltonian for the full Floquet analysis ─────────
#
# The bare ring (source removed) has a mode at w_r with normalized mode
# function psi(x); the reduced Hamiltonian keeps that one field mode
# with the Duffing source as INDEPENDENT (q, p):
#
#   H_red = (P^2 + w_r^2 Q^2)/2 + (p^2 + w0^2 q^2)/2 + mu q^4/4
#           + g_eff q Q ,        g_eff = g psi(0) ,
#
# whose periodic orbit (nonlinear normal mode) is the eigenhistory in
# reduced form, and whose 4x4 variational monodromy over one period is
# the FULL Hamiltonian stability spectrum - (q, p) evolved, not slaved.


def bare_ring_mode():
    if 'BARE' in _CACHE:
        return _CACHE['BARE']

    def tr_bare(w):
        T = (P(w, _L2) @ barrier_M(w, rev=True) @ P(w, _TAU)
             @ barrier_M(w) @ P(w, _L1))
        return float(complex(np.trace(T)).real)

    ws = np.linspace(2.70, 2.78, 400)
    tb = [tr_bare(w) for w in ws]
    w_r = None
    for i in range(len(ws) - 1):
        if (tb[i] - 2) * (tb[i + 1] - 2) < 0:
            w_r = brentq(lambda w: tr_bare(w) - 2, ws[i], ws[i + 1],
                         xtol=1e-12)
            break
    T = (P(w_r, _L2) @ barrier_M(w_r, rev=True) @ P(w_r, _TAU)
         @ barrier_M(w_r) @ P(w_r, _L1))
    _, sv, vh = np.linalg.svd(T - np.eye(2))
    v = vh[-1].conj()
    segs = [(_L1, v)]
    x = barrier_M(w_r) @ (P(w_r, _L1) @ v)
    segs.append((_TAU, x))
    x = barrier_M(w_r, rev=True) @ (P(w_r, _TAU) @ x)
    segs.append((_L2, x))
    norm2 = 0.0
    for Lg, s0 in segs:
        xs = np.linspace(0, Lg, 2001)
        prof = (s0[0] * np.exp(1j * w_r * xs)
                + s0[1] * np.exp(-1j * w_r * xs))
        norm2 += np.trapezoid(np.abs(prof) ** 2, xs)
    psi0 = abs(v[0] + v[1]) / math.sqrt(norm2)
    out = {'w_r': float(w_r), 'psi0': float(psi0),
           'g_eff': float(_G * psi0), 'bare_resid': float(sv[-1])}
    _CACHE['BARE'] = out
    return out


def _h_red(y, w_r, g_eff):
    Q, Pq, q, pp = y
    return ((Pq ** 2 + w_r ** 2 * Q ** 2) / 2
            + (pp ** 2 + _W0 ** 2 * q ** 2) / 2
            + _MU * q ** 4 / 4 + g_eff * q * Q)


def _rhs_red(w_r, g_eff):
    def rhs(t, y):
        Q, Pq, q, pp = y
        return [Pq, -w_r ** 2 * Q - g_eff * q,
                pp, -_W0 ** 2 * q - _MU * q ** 3 - g_eff * Q]
    return rhs


def _rhs_var(w_r, g_eff):
    base = _rhs_red(w_r, g_eff)

    def rhs(t, y):
        q = y[2]
        A = np.array([[0, 1, 0, 0],
                      [-w_r ** 2, 0, -g_eff, 0],
                      [0, 0, 0, 1],
                      [-g_eff, 0, -(_W0 ** 2 + 3 * _MU * q ** 2), 0]])
        M = y[4:].reshape(4, 4)
        return np.concatenate([base(t, y[:4]), (A @ M).ravel()])
    return rhs


# ========================================================================
# T6. The full Hamiltonian stability spectrum
# ========================================================================


def test_T6_stability() -> dict:
    t4 = test_T4_joint_solve()
    w_s, U_s = t4['w_star'], t4['U0_star']
    bare = bare_ring_mode()
    w_r, psi0, g_eff = bare['w_r'], bare['psi0'], bare['g_eff']

    # the reduced-model initial guess from harmonic balance
    Q_amp = U_s / psi0
    a_hb = source_a(w_s, U_s)
    E_red = _h_red([Q_amp, 0, a_hb, 0], w_r, g_eff)
    rhs = _rhs_red(w_r, g_eff)

    # (i) the periodic orbit (nonlinear normal mode) by shooting:
    # turning-point start (P = p = 0); unknowns (Q0, q0, T); equations
    # P(T) = 0, p(T) = 0, H_red = E_red
    def shot(Q0, q0, Tp):
        sol = solve_ivp(rhs, (0, Tp), [Q0, 0, q0, 0], rtol=1e-12,
                        atol=1e-14, method="DOP853")
        y = sol.y[:, -1]
        return np.array([y[1], y[3]])

    xk = np.array([Q_amp, a_hb, 2 * math.pi / w_s])
    F3 = np.array([1.0, 1.0, 1.0])
    for _ in range(40):
        Q0, q0, Tp = xk
        F3 = np.concatenate([shot(Q0, q0, Tp),
                             [_h_red([Q0, 0, q0, 0], w_r, g_eff)
                              - E_red]])
        if np.abs(F3).max() < 1e-11:
            break
        J = np.zeros((3, 3))
        h = 1e-7
        for j in range(3):
            d = np.zeros(3)
            d[j] = h
            Q1, q1, T1 = xk + d
            F1 = np.concatenate([shot(Q1, q1, T1),
                                 [_h_red([Q1, 0, q1, 0], w_r, g_eff)
                                  - E_red]])
            J[:, j] = (F1 - F3) / h
        xk = xk - np.linalg.solve(J, F3)
    orbit_resid = float(np.abs(F3).max())
    w_nnm = 2 * math.pi / xk[2]
    freq_consistency = abs(w_nnm - w_s) / w_s

    # (ii) the 4x4 variational monodromy over one period: the FULL
    # Hamiltonian spectrum with (q, p) evolved as independent variables
    y0 = np.concatenate([[xk[0], 0, xk[1], 0], np.eye(4).ravel()])
    sol = solve_ivp(_rhs_var(w_r, g_eff), (0, xk[2]), y0, rtol=1e-12,
                    atol=1e-14, method="DOP853")
    Mf = sol.y[4:, -1].reshape(4, 4)
    ev = np.linalg.eigvals(Mf)
    mods = np.abs(ev)
    det_M = float(np.linalg.det(Mf).real)
    Jsym = np.array([[0, 1, 0, 0], [-1, 0, 0, 0],
                     [0, 0, 0, 1], [0, 0, -1, 0]], dtype=float)
    sympl_resid = float(np.max(np.abs(Mf.T @ Jsym @ Mf - Jsym)))

    # identify the pairs: the Floquet-trivial pair at 1 (along-flow +
    # energy) and the source pair rotating at its dressed frequency
    idx = np.argsort(np.abs(ev - 1))
    trivial = ev[idx[:2]]
    nontriv = ev[idx[2:]]
    theta = float(abs(np.angle(nontriv[0])))
    w_pair = (theta + 2 * math.pi) / xk[2]
    dressed_dev = abs(w_pair - _W0) / _W0

    # (iii) energy conservation along the orbit
    solE = solve_ivp(rhs, (0, xk[2]), [xk[0], 0, xk[1], 0], rtol=1e-12,
                     atol=1e-14, dense_output=True, method="DOP853")
    Hs = [_h_red(solE.sol(t), w_r, g_eff)
          for t in np.linspace(0, xk[2], 50)]
    e_drift = float(max(Hs) - min(Hs))

    # (iv) the slaved harmonic-balance winding-map spectrum, retained
    # for comparison (this is what v1 reported; it cannot see the
    # source pair)
    v, _ = scaled_state(w_s, U_s)

    def themap(vr):
        vv = vr[:2] + 1j * vr[2:]
        return np.concatenate([
            (ring_T(w_s, abs(vv[0] + vv[1])) @ vv).real,
            (ring_T(w_s, abs(vv[0] + vv[1])) @ vv).imag])

    X0 = np.concatenate([v.real, v.imag])
    Jm = np.zeros((4, 4))
    h = 1e-7
    for j in range(4):
        dp = X0.copy()
        dp[j] += h
        dm = X0.copy()
        dm[j] -= h
        Jm[:, j] = (themap(dp) - themap(dm)) / (2 * h)
    ev_slaved = np.linalg.eigvals(Jm)

    ok = (orbit_resid < 1e-10
          and freq_consistency < 1e-3
          and np.all(np.abs(mods - 1.0) < 1e-9)
          and abs(det_M - 1.0) < 1e-10
          and sympl_resid < 1e-10
          and np.max(np.abs(trivial - 1.0)) < 1e-6
          and dressed_dev < 0.05
          and e_drift < 1e-9)
    return {
        'name': 'T6_stability',
        'description': (
            'the FULL Hamiltonian stability spectrum: the Duffing '
            '(q, p) evolved as independent variables in the 4x4 '
            'variational monodromy of the reduced 2-dof Hamiltonian '
            'about its shooting-refined periodic orbit - the '
            'Floquet-trivial pair at 1 plus the SOURCE pair rotating '
            'at its dressed frequency (invisible to the slaved '
            'harmonic-balance map); all on the unit circle, '
            'symplectic to machine precision'
        ),
        'bare_ring_mode': bare,
        'orbit': {'Q0': float(xk[0]), 'q0': float(xk[1]),
                  'period': float(xk[2]),
                  'shooting_residual': orbit_resid},
        'nnm_frequency': float(w_nnm),
        'ring_frequency': float(w_s),
        'reduction_consistency': float(freq_consistency),
        'floquet_eigenvalues': [[float(e.real), float(e.imag)]
                                for e in ev],
        'eigenvalue_moduli': [float(m) for m in mods],
        'max_modulus_deviation': float(np.max(np.abs(mods - 1))),
        'trivial_pair_deviation': float(np.max(np.abs(trivial - 1.0))),
        'source_pair_frequency': float(w_pair),
        'source_bare_frequency': _W0,
        'dressed_frequency_deviation': float(dressed_dev),
        'det_monodromy': det_M,
        'symplectic_residual': sympl_resid,
        'energy_drift_on_orbit': e_drift,
        'slaved_hb_spectrum': [[float(e.real), float(e.imag)]
                               for e in ev_slaved],
        'pass': bool(ok),
    }


# ========================================================================
# T7. Honest scope
# ========================================================================


def test_T7_honest_scope() -> dict:
    t4 = test_T4_joint_solve()
    a = abs(t4['a_star'])
    third_harm = 0.75 * _MU * a ** 2 / (8 * t4['w_star'] ** 2)
    scope = [
        'Harmonic balance keeps the fundamental only; the neglected '
        f'third-harmonic weight is O(mu a^2/8w^2) ~ {third_harm:.1e} '
        'at the fixed point - small, and computable if needed.',
        'The Duffing quartic is the MINIMAL conservative '
        'anharmonicity (the cubic is absent by q -> -q symmetry); '
        'any conservative source with a nonlinear frequency pull '
        'plays the same role.',
        'E_0 is a specified energy budget; its physical value (e.g. '
        'the #58 nucleation quantum) is program-level input, not '
        'derived here.  hbar is still not derived: the scales are '
        '(w0, mu, g, E_0).',
        'The corrected ledger includes the time-averaged '
        'interaction energy <g q u(0)> (negative, ~0.5% of E_0 '
        'here); with it the working point shifts (U0*: 0.9122 -> '
        '0.9144) and closure holds to 1e-13.',
        'The full stability spectrum is computed on the reduced '
        '2-dof Hamiltonian (the bare ring mode + the Duffing source, '
        'g_eff = g psi(0) derived from the bare eigenvector); higher '
        'ring modes are dropped - the reduction consistency metric '
        'is the NNM-vs-ring frequency match (8.8e-5).  The slaved '
        'harmonic-balance map iterates WINDINGS; the Floquet '
        'monodromy evolves TIME over one carrier period - the two '
        'spectra answer different questions, and only the latter '
        'sees the source pair.',
        'The mouth ports are cubic-spline interpolants of the '
        'unitarized Tangherlini greybody (re-unitarized pointwise); '
        'the raw-port systematic at the fixed point is '
        f"{t4['raw_port_systematic']:.1e} in tr - the working point "
        'is defined in the interpolated model, anchored to the '
        'physical solver at that tolerance.',
        'The ring is the monochromatic (carrier-closed) skeleton; '
        'packet eigenhistories need the group condition (named '
        'successor since #218).',
        'Classical, zonal scalar, frozen geometry; MTY network '
        'history posited.',
    ]
    return {
        'name': 'T7_honest_scope',
        'description': 'what this PR does and does not establish',
        'scope': scope,
        'pass': True,
    }


# ========================================================================
# T8. Assessment
# ========================================================================


def test_T8_assessment() -> dict:
    t2 = test_T2_source_derived()
    t3 = test_T3_ring_modes()
    t4 = test_T4_joint_solve()
    t5 = test_T5_energy_closure()
    t6 = test_T6_stability()
    core = all(t['pass'] for t in (t2, t3, t4, t5, t6))
    assessment = (
        'Nothing about the source is imposed anymore, and nothing '
        'about its dynamics is frozen: the scattering is unitary '
        'because the Hamiltonian is conservative, #218\'s law is its '
        'weak-coupling shadow, the ledger now carries the interaction '
        'energy <g q u(0)> explicitly, and the amplitude is fixed by '
        'the CORRECTED total-energy closure on the nonlinear mode '
        'branch.  The stability question is answered at the full '
        'Hamiltonian level: with the Duffing (q, p) evolved as '
        'independent variables in the variational monodromy about the '
        'shooting-refined periodic orbit, the spectrum is the '
        'Floquet-trivial pair at 1 plus the source pair rotating at '
        'its dressed frequency - the degree of freedom the slaved '
        'harmonic-balance map could not see - all on the unit circle, '
        'symplectic to machine precision: the conservative '
        'eigenhistory is marginal in the full Hamiltonian sense, '
        'Novikov-passive, with no runaway direction anywhere in its '
        'phase space.'
        if core else
        'INCONCLUSIVE - a core check failed; do not quote.'
    )
    return {
        'name': 'T8_assessment',
        'description': 'the standing of the Hamiltonian-source solve',
        'core_green': bool(core),
        'assessment': assessment,
        'pass': bool(core),
    }


# ========================================================================
# Probe runner
# ========================================================================


def run_probe() -> dict:
    tests = [
        test_T1_goal(),
        test_T2_source_derived(),
        test_T3_ring_modes(),
        test_T4_joint_solve(),
        test_T5_energy_closure(),
        test_T6_stability(),
        test_T7_honest_scope(),
        test_T8_assessment(),
    ]
    t2, t3, t4, t5, t6 = tests[1], tests[2], tests[3], tests[4], tests[5]
    all_ok = all(t['pass'] for t in tests)
    if all_ok:
        verdict_class = (
            "THE_SOURCE_PHASE_LAW_IS_DERIVED_FROM_A_HAMILTONIAN_THE_"
            "EIGENHISTORY_SOLVES_UXX_EQUALS_X_WITH_TOTAL_ENERGY_"
            "CLOSURE_AND_ITS_STABILITY_SPECTRUM_IS_MARGINAL"
        )
        verdict = (
            "ESTABLISHED (the argument is in "
            "docs/hamiltonian_source_eigenhistory.md).\n\n"
            "THE SOURCE, DERIVED. H = p^2/2 + w0^2 q^2/2 + mu q^4/4 "
            "+ g q u(0): unitary scattering as a CONSEQUENCE "
            f"({t2['unitarity_deviation']:.0e}), reactive "
            f"({t2['net_power_deviation']:.0e}), and #218's imposed "
            "law recovered as the weak-coupling limit (ratios "
            f"{t2['weak_coupling_limit_ratios']}).\n\n"
            "THE RING. The reflecting source makes the loop a true "
            "two-direction ring: real trace, the tr > 2 gap segment "
            f"resolved ({t3['gap_segment'][0]:.4f}-"
            f"{t3['gap_segment'][1]:.4f}, max {t3['gap_max_trace']:.4f} "
            "- cleaning up #217's coarse-scan tangency), and the "
            "homogeneous condition defining a nonlinear mode BRANCH "
            f"w*(A) (span {t3['branch_frequency_span']:.1e}): "
            "homogeneity alone cannot fix the amplitude.\n\n"
            "THE JOINT SOLVE, CORRECTED LEDGER. U(X)X = X TOGETHER "
            "WITH E_field + E_source + <g q u(0)> = E0: (w*, U0*) = "
            f"({t4['w_star']:.6f}, {t4['U0_star']:.6f}), a* = "
            f"{t4['a_star']:.4f}; FIXED-POINT RESIDUALS: Newton "
            f"({t4['newton_residuals'][0]:.0e}, "
            f"{t4['newton_residuals'][1]:.0e}), U(X)X = X "
            f"{t4['UX_equals_X_residual']:.0e}, source slaving "
            f"{t4['source_slaving_residual']:.0e}; energy partition "
            f"E_field = {t4['E_field']:.4f} + E_source = "
            f"{t4['E_source']:.4f} + E_int = {t4['E_interaction']:.5f} "
            f"= E0 to {t4['energy_closure_residual']:.0e}. Flux "
            f"constant around the ring ({t5['flux_constancy']:.0e}); "
            f"1e4-pass persistence (drift {t5['norm_drift_1e4']:.0e})."
            "\n\n"
            "THE FULL HAMILTONIAN STABILITY SPECTRUM. The Duffing "
            "(q, p) evolved as INDEPENDENT variables in the 4x4 "
            "variational monodromy about the shooting-refined "
            f"periodic orbit (residual {t6['orbit']['shooting_residual']:.0e}; "
            "NNM frequency consistent with the full ring to "
            f"{t6['reduction_consistency']:.1e}): eigenvalues "
            f"{t6['floquet_eigenvalues'][2]} (double - the Floquet-"
            f"trivial pair, to {t6['trivial_pair_deviation']:.0e}) and "
            f"{t6['floquet_eigenvalues'][0]} - THE SOURCE PAIR, "
            f"rotating at its dressed frequency "
            f"{t6['source_pair_frequency']:.4f} (bare w0 = "
            f"{t6['source_bare_frequency']}), invisible to the slaved "
            "harmonic-balance map; ALL moduli within "
            f"{t6['max_modulus_deviation']:.0e} of the unit circle; "
            f"det M = 1 to {abs(t6['det_monodromy']-1):.0e}, "
            f"symplectic to {t6['symplectic_residual']:.0e}, energy "
            f"drift {t6['energy_drift_on_orbit']:.0e}: the "
            "conservative eigenhistory is marginal in the FULL "
            "Hamiltonian sense - Novikov-passive, no runaway."
        )
    else:
        verdict_class = "HAMILTONIAN_SOURCE_INCONCLUSIVE"
        verdict = (
            "INCONCLUSIVE. A core check failed; re-examine before "
            "quoting the solve."
        )

    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "identification": (
            "The Hamiltonian source: #218's imposed phase law replaced "
            "by a minimal conservative source (side-coupled Duffing, "
            "explicit Hamiltonian, unitary-by-construction harmonic-"
            "balance scattering); the eigenhistory solved as U(X)X = X "
            "together with total-energy closure on the two-direction "
            "ring; fixed-point residuals at machine precision and a "
            "marginal (all-unit-circle) stability spectrum reported"
        ),
        "executes": (
            "the requested model replacement and joint solve"
        ),
        "tests": tests,
        "n_passed": sum(1 for t in tests if t["pass"]),
        "n_total": len(tests),
        "verdict_class": verdict_class,
        "verdict": verdict,
    }


def _json_default(o):
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.ndarray):
        return "<array>"
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, complex):
        return [o.real, o.imag]
    return str(o)


def render_markdown(s: dict) -> str:
    out = []
    out.append("# The Hamiltonian source eigenhistory (PR #219)")
    out.append("")
    out.append(f"**Run:** {s['timestamp_utc']}")
    out.append("")
    out.append(
        "The deliverable is `docs/hamiltonian_source_eigenhistory.md` "
        "- the imposed source law replaced by an explicit Hamiltonian, "
        "and U(X)X = X solved with total-energy closure. *(QFT on the "
        "fixed classical throat geometry, not quantum gravity.)*"
    )
    out.append("")
    out.append("## Test summary")
    out.append("")
    out.append("| # | Test | Key finding | PASS? |")
    out.append("|---|---|---|---|")
    labels = {
        "T1": "the goal: derive the source, solve jointly",
        "T2": "the source from H: unitary/reactive by construction",
        "T3": "the ring: the gap resolved; a mode BRANCH w*(A)",
        "T4": "U(X)X = X + energy closure: residuals at 1e-14",
        "T5": "flux constant; zero source power; 1e4-pass persistence",
        "T6": "the FULL Hamiltonian Floquet spectrum ((q,p) evolved)",
        "T7": "honest scope",
        "T8": "assessment",
    }
    for t in s["tests"]:
        p = "**PASS**" if t["pass"] else "**FAIL**"
        pre = t["name"][:2]
        out.append(f"| {pre} | `{t['name']}` | {labels.get(pre,'-')} | {p} |")
    out.append("")
    out.append("## Verdict")
    out.append("")
    out.append(f"**Class:** `{s['verdict_class']}`")
    out.append("")
    out.append(s["verdict"])
    out.append("")
    return "\n".join(out)


def main(argv: Optional[list] = None) -> int:
    summary = run_probe()
    md = render_markdown(summary)
    print(md)
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = here / "runs" / f"{ts}_hamiltonian_source_eigenhistory_probe"
    out.mkdir(parents=True, exist_ok=True)
    (out / "probe.json").write_text(
        json.dumps(summary, indent=2, default=_json_default))
    (out / "probe.md").write_text(md)
    print(f"\nWrote: {out / 'probe.json'}")
    print(f"Wrote: {out / 'probe.md'}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
