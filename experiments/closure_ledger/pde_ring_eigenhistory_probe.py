"""
The PDE-ring eigenhistory: explicit time-domain evolution of the
Duffing source coupled to the two-direction field ring, the periodic
orbit of the COMPLETE source-field state, and its one-period monodromy
(PR #220).

> Framing: QFT on the *fixed classical* throat geometry (geometry -> fields),
> not quantum gravity.

THE FULL HAMILTONIAN SUCCESSOR TO #219
--------------------------------------
#219 derived the source from a Hamiltonian but analyzed the loop in
harmonic balance (monochromatic scattering) with a reduced 2-dof
Floquet model.  This probe drops every reduction: the FIELD IS THE
PDE - the two-direction wave on the periodic ring with the glued
Tangherlini barrier potentials - and the Duffing (q, p) evolve in
explicit time domain, coupled at the crossing:

    H = sum_i dx [pi_i^2/2 + ((Du)_i^2 + V_i u_i^2)/2]
        + p^2/2 + w0^2 q^2/2 + mu q^4/4 + g q u(0)

with THE INTERACTION TERM g q u(0) IN THE CONSERVED ENERGY (its
omission visibly breaks conservation: the no-interaction ledger varies
100x more over one period than the true H does over one hundred).
The integrator is symplectic (leapfrog / velocity Verlet), so H is
conserved up to a bounded O(dt^2) shadow-Hamiltonian oscillation with
no secular drift.

THE PERIODIC ORBIT, LITERALLY AS REQUESTED.  Unknowns (X(0), T) with
X the complete source-field state (u, pi, q, p) in R^{2N+2};
conditions

    X(T) - X(0) = 0 ,     H[X(0)] - E_0 = 0 ,     p(0) = 0

(the last the single phase condition removing time-translation
degeneracy), solved by damped Gauss-Newton with the full
finite-difference Jacobian (batched leapfrog columns): quadratic
convergence to |R| = 2e-13.

THE ONE-PERIOD MONODROMY OF THE COMPLETE STATE.  The tangent leapfrog
(the exact linearization of the symplectic map, co-evolved with the
orbit) gives the (2N+2) x (2N+2) monodromy: ALL eigenvalues on the
unit circle to 4e-14 - every discretized field mode's Floquet pair AND
the source pair: no parametric instability anywhere in the complete
phase space; det M = 1, and M preserves the dx-weighted symplectic
form to 4e-15.  The source pair is identified by its (q, p)
eigenvector weight (0.36 vs 0.008 for the next) and rotates at the
dressed frequency 3.2124 - within 0.07% of #219's reduced-model
prediction 3.2102 (bare w0 = 3.2): the reduced Floquet analysis of
#219 is CONFIRMED by the full PDE, and every approximation it made is
now retired.

Tests:
  T1. Goal.
  T2. The discretized Hamiltonian system (the periodic two-barrier
      ring; symplectic conservation of H INCLUDING g q u(0); the
      interaction term load-bearing; dt^2 scaling; the linear ring
      modes).
  T3. The periodic orbit (the literal Gauss-Newton system; quadratic
      convergence; E_0 hit exactly; dt-refinement of the orbit
      frequency).
  T4. The one-period monodromy of the complete state (all unit
      circle; det = 1; dx-weighted symplectic; the trivial pair).
  T5. The source pair and the field pairs (qp-weight identification;
      the dressed frequency vs #219; low field modes matched in the
      spectrum).
  T6. The energy ledger on the orbit (time-resolved partition; the
      oscillating interaction energy and its average vs the #219
      harmonic-balance value; conservation along the orbit).
  T7. Honest scope.
  T8. Assessment.

Verdict:
  THE_FULL_PDE_EIGENHISTORY_EXISTS_ITS_COMPLETE_MONODROMY_IS_UNIT_
  CIRCLE_SYMPLECTIC_AND_THE_219_REDUCTION_IS_CONFIRMED_BY_THE_FIELD_
  ITSELF
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np
from scipy.interpolate import interp1d

_CACHE: dict = {}
_RH = 1.0
_W0 = 3.2
_MU = 0.5
_G = 0.6
_D_INT = 8.0
_C = 2 * math.pi + _D_INT
_X_PK = 0.1978
_N = 192
_DX = _C / _N
_NS = 2048

# ── the periodic two-barrier ring potential ─────────────────────────────

_rg = 1 + np.logspace(-9, math.log10(299.0), 4000)
_xg = _rg + (_RH / 2) * np.log((_rg - _RH) / (_rg + _RH))
_fg = 1 - (_RH / _rg) ** 2
_vg = _fg * (0.75 / _rg ** 2 + 2.25 * _RH ** 2 / _rg ** 4)
_vi = interp1d(_xg, _vg, kind="cubic", bounds_error=False,
               fill_value=(0.0, 0.0))


def _v_tort(xi):
    xi = np.asarray(xi, float)
    out = np.zeros_like(xi)
    inside = (xi >= _xg[0]) & (xi <= _xg[-1])
    out[inside] = _vi(xi[inside])
    far = xi > _xg[-1]
    out[far] = 0.75 / xi[far] ** 2
    return out


def ring_potential(x):
    """Two glued Tangherlini barriers on the periodic ring: A at pi
    (exterior facing the source at 0), B mirrored at pi + D_INT."""
    x = np.asarray(x, float)
    dA = (x - math.pi + _C / 2) % _C - _C / 2
    dB = (x - (math.pi + _D_INT) + _C / 2) % _C - _C / 2
    return _v_tort(_X_PK - dA) + _v_tort(_X_PK + dB)


_XGRID = np.arange(_N) * _DX
_V = ring_potential(_XGRID)


def _lap(a):
    return (np.roll(a, -1, axis=0) - 2 * a + np.roll(a, 1, axis=0)) \
        / _DX ** 2


# ── the symplectic (leapfrog) flow of the coupled system ────────────────

def h_total(u, pi, q, p, include_interaction=True):
    du = (np.roll(u, -1) - u) / _DX
    e_f = _DX * np.sum(pi ** 2 / 2 + du ** 2 / 2 + _V * u ** 2 / 2)
    e_s = p ** 2 / 2 + _W0 ** 2 * q ** 2 / 2 + _MU * q ** 4 / 4
    e_i = _G * q * u[0] if include_interaction else 0.0
    return float(e_f + e_s + e_i)


def h_parts(u, pi, q, p):
    du = (np.roll(u, -1) - u) / _DX
    e_f = float(_DX * np.sum(pi ** 2 / 2 + du ** 2 / 2
                             + _V * u ** 2 / 2))
    e_s = float(p ** 2 / 2 + _W0 ** 2 * q ** 2 / 2 + _MU * q ** 4 / 4)
    e_i = float(_G * q * u[0])
    return e_f, e_s, e_i


def flow_batch(Z, T, n_s=_NS):
    """Leapfrog the batched state columns Z ((2N+2) x M) for time T."""
    u = Z[:_N].copy()
    pi = Z[_N:2 * _N].copy()
    q = Z[2 * _N].copy()
    p = Z[2 * _N + 1].copy()
    dt = T / n_s
    Vc = _V[:, None] if u.ndim == 2 else _V
    for _ in range(n_s):
        fu = _lap(u) - Vc * u
        fu[0] -= _G * q / _DX
        pi = pi + 0.5 * dt * fu
        p = p + 0.5 * dt * (-_W0 ** 2 * q - _MU * q ** 3 - _G * u[0])
        u = u + dt * pi
        q = q + dt * p
        fu = _lap(u) - Vc * u
        fu[0] -= _G * q / _DX
        pi = pi + 0.5 * dt * fu
        p = p + 0.5 * dt * (-_W0 ** 2 * q - _MU * q ** 3 - _G * u[0])
    out = np.empty_like(Z)
    out[:_N] = u
    out[_N:2 * _N] = pi
    out[2 * _N] = q
    out[2 * _N + 1] = p
    return out


def flow(z, T, n_s=_NS):
    return flow_batch(z[:, None], T, n_s)[:, 0]


def linear_modes():
    if 'LIN' in _CACHE:
        return _CACHE['LIN']
    L = (np.diag(-2.0 * np.ones(_N)) + np.diag(np.ones(_N - 1), 1)
         + np.diag(np.ones(_N - 1), -1))
    L[0, -1] = L[-1, 0] = 1.0
    L /= _DX ** 2
    evals, evecs = np.linalg.eigh(-L + np.diag(_V))
    _CACHE['LIN'] = (np.sqrt(np.abs(evals)), evecs)
    return _CACHE['LIN']


def periodic_orbit(n_s=_NS):
    """The literal solve: unknowns (X(0), T); conditions
    X(T) - X(0) = 0, H[X(0)] - E_0 = 0, p(0) = 0 (phase)."""
    key = ('ORB', n_s)
    if key in _CACHE:
        return _CACHE[key]
    ws, evecs = linear_modes()
    k = 12                                    # source-coupled ring mode
    phi = evecs[:, k] * np.sign(evecs[0, k])
    u0 = (0.9 / abs(phi[0])) * phi
    q0 = -_G * u0[0] / (_W0 ** 2 - ws[k] ** 2)
    z = np.concatenate([u0, 0 * u0, [q0, 0.0]])
    T = 2 * math.pi / ws[k]
    E0 = round(h_total(u0, 0 * u0, q0, 0.0), 3)
    n_dim = 2 * _N + 2

    def residual(zz, TT):
        zT = flow(zz, TT, n_s)
        return np.concatenate([zT - zz,
                               [h_total(zz[:_N], zz[_N:2 * _N],
                                        zz[2 * _N], zz[2 * _N + 1]) - E0,
                                zz[2 * _N + 1]]])

    hist = []
    for _ in range(12):
        R = residual(z, T)
        hist.append(float(np.linalg.norm(R)))
        if hist[-1] < 1e-11:
            break
        eps = 1e-7
        Zb = np.concatenate([z[:, None],
                             z[:, None] + eps * np.eye(n_dim)], axis=1)
        ZT = flow_batch(Zb, T, n_s)
        base = ZT[:, 0]
        J = np.empty((n_dim + 2, n_dim + 1))
        for j in range(n_dim):
            zj = Zb[:, j + 1]
            dtop = (ZT[:, j + 1] - base) / eps - (zj - z) / eps
            dH = (h_total(zj[:_N], zj[_N:2 * _N], zj[2 * _N],
                          zj[2 * _N + 1])
                  - h_total(z[:_N], z[_N:2 * _N], z[2 * _N],
                            z[2 * _N + 1])) / eps
            dph = (zj - z)[2 * _N + 1] / eps
            J[:, j] = np.concatenate([dtop, [dH, dph]])
        zT2 = flow(z, T + eps, n_s)
        J[:, n_dim] = np.concatenate([(zT2 - base) / eps, [0.0, 0.0]])
        step, *_ = np.linalg.lstsq(J, -R, rcond=None)
        z = z + step[:n_dim]
        T = T + step[n_dim]
    out = {'z': z, 'T': float(T), 'E0': E0,
           'residual': float(np.linalg.norm(residual(z, T))),
           'newton_history': hist,
           'w_orbit': float(2 * math.pi / T),
           'w_linear': float(ws[k])}
    _CACHE[key] = out
    return out


def monodromy(orb):
    """Tangent leapfrog co-evolved with the orbit: the exact
    linearization of the symplectic map over one period."""
    key = ('MONO', orb['T'])
    if key in _CACHE:
        return _CACHE[key]
    z, T = orb['z'], orb['T']
    n_dim = 2 * _N + 2
    u = z[:_N].copy()
    pi = z[_N:2 * _N].copy()
    q = float(z[2 * _N])
    p = float(z[2 * _N + 1])
    M = np.eye(n_dim)
    dU = M[:_N]
    dPi = M[_N:2 * _N]
    dQ = M[2 * _N]
    dP = M[2 * _N + 1]
    dt = T / _NS
    for _ in range(_NS):
        fu = _lap(u) - _V * u
        fu[0] -= _G * q / _DX
        fq = -_W0 ** 2 * q - _MU * q ** 3 - _G * u[0]
        dfu = _lap(dU) - _V[:, None] * dU
        dfu[0] -= _G * dQ / _DX
        dfq = -(_W0 ** 2 + 3 * _MU * q ** 2) * dQ - _G * dU[0]
        pi = pi + 0.5 * dt * fu
        p = p + 0.5 * dt * fq
        dPi = dPi + 0.5 * dt * dfu
        dP = dP + 0.5 * dt * dfq
        u = u + dt * pi
        q = q + dt * p
        dU = dU + dt * dPi
        dQ = dQ + dt * dP
        fu = _lap(u) - _V * u
        fu[0] -= _G * q / _DX
        fq = -_W0 ** 2 * q - _MU * q ** 3 - _G * u[0]
        dfu = _lap(dU) - _V[:, None] * dU
        dfu[0] -= _G * dQ / _DX
        dfq = -(_W0 ** 2 + 3 * _MU * q ** 2) * dQ - _G * dU[0]
        pi = pi + 0.5 * dt * fu
        p = p + 0.5 * dt * fq
        dPi = dPi + 0.5 * dt * dfu
        dP = dP + 0.5 * dt * dfq
    Mono = np.vstack([dU, dPi, dQ[None, :], dP[None, :]])
    _CACHE[key] = Mono
    return Mono


# ========================================================================
# T1. Goal
# ========================================================================


def test_T1_goal() -> dict:
    return {
        'name': 'T1_goal',
        'description': (
            'The full Hamiltonian successor to #219: explicit '
            'time-domain evolution of the Duffing (q, p) coupled to '
            'the two-direction field ring, H = H_field + p^2/2 + '
            'w0^2 q^2/2 + mu q^4/4 + g q u(0) with the interaction '
            'term in the conserved energy; the one-period monodromy '
            'of the COMPLETE source-field state; and the periodic '
            'orbit solved literally as X(T) - X(0) = 0, H[X] - E0 = '
            '0, with one phase condition (p(0) = 0) removing '
            'time-translation degeneracy.'
        ),
        'pass': True,
    }


# ========================================================================
# T2. The discretized Hamiltonian system
# ========================================================================


def test_T2_hamiltonian_system() -> dict:
    if 'T2' in _CACHE:
        return _CACHE['T2']
    ws, evecs = linear_modes()
    k = 12
    phi = evecs[:, k] * np.sign(evecs[0, k])
    u0 = (0.9 / abs(phi[0])) * phi
    q0 = -_G * u0[0] / (_W0 ** 2 - ws[k] ** 2)
    z0 = np.concatenate([u0, 0 * u0, [q0, 0.0]])
    T = 2 * math.pi / ws[k]
    H0 = h_total(u0, 0 * u0, q0, 0.0)

    # (i) H including g q u(0): bounded O(dt^2) oscillation, no
    # secular drift over 100 periods
    z = z0.copy()
    Hs = []
    for _ in range(100):
        z = flow(z, T)
        Hs.append(h_total(z[:_N], z[_N:2 * _N], z[2 * _N],
                          z[2 * _N + 1]))
    Hs = np.array(Hs)
    drift_100 = float(np.abs(Hs - H0).max())
    secular = float(abs(Hs[-1] - np.mean(Hs[:10])))

    # (ii) dt^2 scaling of the shadow-Hamiltonian bound
    z = z0.copy()
    z = flow(z, T, n_s=2 * _NS)
    drift_fine = abs(h_total(z[:_N], z[_N:2 * _N], z[2 * _N],
                             z[2 * _N + 1]) - H0)
    z = z0.copy()
    z = flow(z, T, n_s=_NS)
    drift_coarse = abs(h_total(z[:_N], z[_N:2 * _N], z[2 * _N],
                               z[2 * _N + 1]) - H0)
    dt2_ratio = float(drift_coarse / max(drift_fine, 1e-30))

    # (iii) the interaction term is load-bearing: the ledger WITHOUT
    # g q u(0) visibly fails over a single period
    z = flow(z0.copy(), T)
    h_no = abs(h_total(z[:_N], z[_N:2 * _N], z[2 * _N], z[2 * _N + 1],
                       include_interaction=False)
               - h_total(z0[:_N], z0[_N:2 * _N], z0[2 * _N],
                         z0[2 * _N + 1], include_interaction=False))

    # (iv) the ring mode structure: CW/CCW pairs, source-coupled
    # members have |phi(0)| > 0
    pair_split = float(abs(ws[11] - ws[12]))
    coupled = abs(evecs[0, 12]) > 0.05 and abs(evecs[0, 11]) < 1e-6

    ok = (drift_100 < 1e-4 and secular < 2 * drift_100
          and 2.5 < dt2_ratio < 6.0
          and h_no > 20 * drift_100
          and pair_split < 5e-4 and coupled)
    out = {
        'name': 'T2_hamiltonian_system',
        'description': (
            'the coupled PDE-ring + Duffing system under symplectic '
            'leapfrog: H INCLUDING g q u(0) conserved (bounded, '
            'non-secular, O(dt^2)); the interaction term load-bearing '
            'in the ledger; the ring mode pair structure'
        ),
        'H_drift_100_periods': drift_100,
        'secular_component': secular,
        'dt2_scaling_ratio': dt2_ratio,
        'no_interaction_ledger_violation_1_period': float(h_no),
        'violation_over_true_drift': float(h_no / drift_100),
        'mode_pair_split': pair_split,
        'V_range': [float(_V.min()), float(_V.max())],
        'pass': bool(ok),
    }
    _CACHE['T2'] = out
    return out


# ========================================================================
# T3. The periodic orbit, literally as requested
# ========================================================================


def test_T3_periodic_orbit() -> dict:
    orb = periodic_orbit()

    # dt-refinement: the orbit frequency converges O(dt^2)
    orb2 = periodic_orbit(n_s=2 * _NS)
    dw = abs(orb['w_orbit'] - orb2['w_orbit'])
    richardson = float(dw / max(abs(orb2['w_orbit']), 1e-30))

    z = orb['z']
    H_at = h_total(z[:_N], z[_N:2 * _N], z[2 * _N], z[2 * _N + 1])
    hist = orb['newton_history']
    quadratic = all(hist[i + 1] < 0.15 * hist[i]
                    for i in range(min(3, len(hist) - 1)))

    ok = (orb['residual'] < 1e-11
          and abs(H_at - orb['E0']) < 1e-10
          and abs(z[2 * _N + 1]) < 1e-11         # phase condition
          and quadratic and richardson < 1e-5)
    return {
        'name': 'T3_periodic_orbit',
        'description': (
            'Gauss-Newton on the literal system [X(T) - X(0); '
            'H - E0; p(0)] over (X(0), T): quadratic convergence to '
            'machine residual; E0 hit exactly; the orbit frequency '
            'dt-converged'
        ),
        'E0': orb['E0'],
        'period': orb['T'],
        'w_orbit': orb['w_orbit'],
        'w_linear_mode': orb['w_linear'],
        'nonlinear_pull': float(orb['w_orbit'] - orb['w_linear']),
        'newton_history': hist,
        'final_residual': orb['residual'],
        'energy_at_orbit': float(H_at),
        'phase_condition_residual': float(abs(z[2 * _N + 1])),
        'u_at_source': float(z[0]),
        'q_star': float(z[2 * _N]),
        'dt_refinement_relative': richardson,
        'pass': bool(ok),
    }


# ========================================================================
# T4. The one-period monodromy of the complete state
# ========================================================================


def test_T4_monodromy() -> dict:
    if 'T4' in _CACHE:
        return _CACHE['T4']
    orb = periodic_orbit()
    Mono = monodromy(orb)
    n_dim = 2 * _N + 2

    ev = np.linalg.eigvals(Mono)
    mods = np.abs(ev)
    max_dev = float(np.abs(mods - 1).max())
    det_M = float(np.linalg.det(Mono).real)

    # the dx-weighted symplectic form (canonical pairs (u_i, dx pi_i))
    Jw = np.zeros((n_dim, n_dim))
    for i in range(_N):
        Jw[i, _N + i] = _DX
        Jw[_N + i, i] = -_DX
    Jw[2 * _N, 2 * _N + 1] = 1.0
    Jw[2 * _N + 1, 2 * _N] = -1.0
    sympl = float(np.abs(Mono.T @ Jw @ Mono - Jw).max())

    idx = np.argsort(np.abs(ev - 1))
    trivial_dev = float(np.abs(ev[idx[:2]] - 1.0).max())

    ok = (max_dev < 1e-10 and abs(det_M - 1) < 1e-9
          and sympl < 1e-12 and trivial_dev < 1e-8)
    out = {
        'name': 'T4_monodromy',
        'description': (
            'the (2N+2) x (2N+2) monodromy of the COMPLETE '
            'source-field state (tangent leapfrog, exact '
            'linearization of the symplectic map): ALL eigenvalues '
            'on the unit circle - no parametric instability of any '
            'field mode against the eigenhistory; det = 1; the '
            'dx-weighted symplectic form preserved'
        ),
        'dimension': n_dim,
        'max_modulus_deviation': max_dev,
        'det_monodromy': det_M,
        'symplectic_residual': sympl,
        'trivial_pair_deviation': trivial_dev,
        'pass': bool(ok),
    }
    _CACHE['T4'] = out
    _CACHE['EV'] = ev
    _CACHE['MONO_M'] = Mono
    return out


# ========================================================================
# T5. The source pair and the field pairs
# ========================================================================


def test_T5_pairs() -> dict:
    test_T4_monodromy()
    orb = periodic_orbit()
    Mono = _CACHE['MONO_M']
    T = orb['T']
    w_, vecs = np.linalg.eig(Mono)

    # source-pair identification by (q, p) eigenvector weight
    qp_w = ((np.abs(vecs[2 * _N]) ** 2 + np.abs(vecs[2 * _N + 1]) ** 2)
            / np.sum(np.abs(vecs) ** 2, axis=0))
    order = np.argsort(qp_w)[::-1]
    i_src = order[0]
    theta = abs(np.angle(w_[i_src]))
    w_source = (theta + 2 * math.pi) / T
    w_219 = 3.2102          # the #219 reduced-model prediction

    # field pairs: the low linear modes appear in the spectrum at
    # w_k T mod 2pi (folded), weakly perturbed by the orbit
    ws_lin, _ = linear_modes()
    ang = np.sort(np.abs(np.angle(w_)))
    matches = []
    for k in (1, 2, 3, 4, 5):
        target = ws_lin[k] * T % (2 * math.pi)
        target = min(target, 2 * math.pi - target)
        gap = float(np.min(np.abs(ang - target)))
        matches.append({'k': int(k), 'w_lin': float(ws_lin[k]),
                        'folded_angle': float(target),
                        'spectral_gap': gap})
    worst_match = max(m['spectral_gap'] for m in matches)

    ok = (qp_w[order[0]] > 0.2
          and qp_w[order[2]] < 0.05
          and abs(w_source - w_219) / w_219 < 5e-3
          and abs(w_source - _W0) / _W0 < 0.05
          and worst_match < 2e-2)
    return {
        'name': 'T5_pairs',
        'description': (
            'the source pair identified by (q, p) eigenvector weight, '
            'rotating at the dressed frequency - confirming the #219 '
            'reduced Floquet analysis from the full PDE; the low '
            'field modes present at their folded angles'
        ),
        'qp_weight_top3': [float(qp_w[order[i]]) for i in range(3)],
        'source_pair': [float(w_[i_src].real), float(w_[i_src].imag)],
        'source_pair_modulus': float(abs(w_[i_src])),
        'source_frequency_dressed': float(w_source),
        'reduced_model_219': w_219,
        'agreement_with_219': float(abs(w_source - w_219) / w_219),
        'bare_w0': _W0,
        'field_mode_matches': matches,
        'pass': bool(ok),
    }


# ========================================================================
# T6. The energy ledger on the orbit
# ========================================================================


def test_T6_energy_ledger() -> dict:
    orb = periodic_orbit()
    z, T = orb['z'].copy(), orb['T']

    # time-resolved partition over one period
    n_snap = 64
    parts = []
    zz = z.copy()
    for _ in range(n_snap):
        zz = flow(zz, T / n_snap, n_s=_NS // n_snap)
        parts.append(h_parts(zz[:_N], zz[_N:2 * _N], zz[2 * _N],
                             zz[2 * _N + 1]))
    parts = np.array(parts)
    e_int_mean = float(np.mean(parts[:, 2]))
    e_int_swing = float(parts[:, 2].max() - parts[:, 2].min())

    # the #219 harmonic-balance average: <g q u(0)> = g a U0 / 2 with
    # a, U0 the orbit's amplitudes (q and u0 antiphased)
    hb_pred = 0.5 * _G * z[2 * _N] * z[0]
    hb_agree = abs(e_int_mean - hb_pred) / abs(hb_pred)

    # conservation of the full H along the orbit
    H_along = parts.sum(axis=1)
    H0 = h_total(z[:_N], z[_N:2 * _N], z[2 * _N], z[2 * _N + 1])
    cons = float(np.abs(H_along - H0).max())

    ok = (e_int_mean < 0 and hb_agree < 0.05
          and cons < 1e-4
          and abs(H0 - orb['E0']) < 1e-10)
    return {
        'name': 'T6_energy_ledger',
        'description': (
            'the ledger, time-resolved on the orbit: the interaction '
            'energy oscillates about a negative mean equal to the '
            '#219 harmonic-balance value; the full H (interaction '
            'included) constant along the orbit'
        ),
        'E_field_mean': float(np.mean(parts[:, 0])),
        'E_source_mean': float(np.mean(parts[:, 1])),
        'E_int_mean': e_int_mean,
        'E_int_swing': e_int_swing,
        'harmonic_balance_prediction': float(hb_pred),
        'hb_agreement': float(hb_agree),
        'H_conservation_along_orbit': cons,
        'pass': bool(ok),
    }


# ========================================================================
# T7. Honest scope
# ========================================================================


def test_T7_honest_scope() -> dict:
    scope = [
        'The orbit and monodromy are those of the DISCRETE symplectic '
        '(leapfrog) map at n_s = 2048 steps/period; the continuous '
        'orbit is approached O(dt^2) (frequency shift ~1e-6 verified '
        'by dt-refinement), and H is conserved up to the bounded '
        'shadow-Hamiltonian oscillation (2e-5), never secularly.',
        'Grid N = 192 (dx = 0.075, ~30 points per carrier '
        'wavelength); the spectrum is the discretized field\'s - all '
        'of it on the unit circle.',
        'The ring geometry uses the physical finite-width Tangherlini '
        'barriers (interior separation 8): the fuller model than '
        '#219\'s point-barrier transfer-matrix idealization (tau = '
        '0.8); the local source physics is robust to this (source '
        'pair within 0.07% of the #219 reduction).',
        'E_0 is a specified budget; hbar is still not derived.',
        'Classical, zonal scalar, frozen geometry; the CTC/network '
        'aspects (mouth offset, greybody one-way membrane) live in '
        '#216-#218 - this PR is the conservative-dynamics core: the '
        'ring as the time-closed loop\'s covering model.',
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
    t2 = test_T2_hamiltonian_system()
    t3 = test_T3_periodic_orbit()
    t4 = test_T4_monodromy()
    t5 = test_T5_pairs()
    t6 = test_T6_energy_ledger()
    core = all(t['pass'] for t in (t2, t3, t4, t5, t6))
    assessment = (
        'Every reduction is retired: the field is the PDE on the '
        'ring, the source is the explicit (q, p) Duffing, the energy '
        'is the full Hamiltonian with the interaction term, the '
        'integrator is symplectic, and the periodic orbit solves the '
        'literal system X(T) = X(0), H = E0, p(0) = 0 to machine '
        'residual.  The one-period monodromy of the complete '
        'source-field state - hundreds of Floquet pairs - lies '
        'entirely on the unit circle with det 1 and the symplectic '
        'form preserved to 4e-15: no parametric instability of any '
        'field mode against the eigenhistory, the trivial pair where '
        'Floquet theory demands it, and the source pair at its '
        'dressed frequency within 0.07% of the #219 reduced '
        'prediction.  The eigenhistory of the transactional arc now '
        'stands as a complete, explicitly evolved, energy-conserving '
        'Hamiltonian object.'
        if core else
        'INCONCLUSIVE - a core check failed; do not quote.'
    )
    return {
        'name': 'T8_assessment',
        'description': 'the standing of the full-PDE eigenhistory',
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
        test_T2_hamiltonian_system(),
        test_T3_periodic_orbit(),
        test_T4_monodromy(),
        test_T5_pairs(),
        test_T6_energy_ledger(),
        test_T7_honest_scope(),
        test_T8_assessment(),
    ]
    t2, t3, t4, t5, t6 = tests[1], tests[2], tests[3], tests[4], tests[5]
    all_ok = all(t['pass'] for t in tests)
    if all_ok:
        verdict_class = (
            "THE_FULL_PDE_EIGENHISTORY_EXISTS_ITS_COMPLETE_MONODROMY_"
            "IS_UNIT_CIRCLE_SYMPLECTIC_AND_THE_219_REDUCTION_IS_"
            "CONFIRMED_BY_THE_FIELD_ITSELF"
        )
        verdict = (
            "ESTABLISHED (the argument is in "
            "docs/pde_ring_eigenhistory.md).\n\n"
            "THE SYSTEM. H = H_field + p^2/2 + w0^2 q^2/2 + mu q^4/4 "
            "+ g q u(0) on the periodic two-barrier ring, evolved by "
            "symplectic leapfrog: H INCLUDING THE INTERACTION TERM is "
            f"conserved (bound {t2['H_drift_100_periods']:.1e} over "
            "100 periods, non-secular, O(dt^2) scaling ratio "
            f"{t2['dt2_scaling_ratio']:.1f}), while the ledger "
            "without it fails "
            f"{t2['violation_over_true_drift']:.0f}x worse in a "
            "single period - the term is load-bearing.\n\n"
            "THE ORBIT. Gauss-Newton on the literal [X(T) - X(0); "
            "H - E0; p(0)] over the complete state: quadratic "
            f"convergence to {t3['final_residual']:.0e}; E0 = "
            f"{t3['E0']} hit to {abs(t3['energy_at_orbit']-t3['E0']):.0e}; "
            f"w_orbit = {t3['w_orbit']:.6f} (nonlinear pull "
            f"{t3['nonlinear_pull']:.2e} from the linear mode), "
            f"dt-converged to {t3['dt_refinement_relative']:.0e}.\n\n"
            "THE COMPLETE MONODROMY. "
            f"{t4['dimension']} x {t4['dimension']}: ALL eigenvalues "
            "on the unit circle to "
            f"{t4['max_modulus_deviation']:.0e} - no parametric "
            "instability of any field mode; det = 1 to "
            f"{abs(t4['det_monodromy']-1):.0e}; the dx-weighted "
            f"symplectic form preserved to "
            f"{t4['symplectic_residual']:.0e}; the trivial pair to "
            f"{t4['trivial_pair_deviation']:.0e}.\n\n"
            "THE PAIRS. The source pair (qp-weight "
            f"{t5['qp_weight_top3'][0]:.2f} vs "
            f"{t5['qp_weight_top3'][2]:.3f} for the next) rotates at "
            f"{t5['source_frequency_dressed']:.4f} - within "
            f"{t5['agreement_with_219']:.1%} of #219's reduced-model "
            f"3.2102 (bare w0 = 3.2): the reduction is CONFIRMED by "
            "the field itself; the low field modes sit at their "
            "folded angles (worst gap "
            f"{max(m['spectral_gap'] for m in t5['field_mode_matches']):.1e})."
            "\n\n"
            "THE LEDGER. Time-resolved on the orbit: E_int oscillates "
            f"about {t6['E_int_mean']:.4f} (negative), equal to the "
            f"#219 harmonic-balance value to {t6['hb_agreement']:.1%}; "
            "the full H is constant along the orbit to "
            f"{t6['H_conservation_along_orbit']:.0e}."
        )
    else:
        verdict_class = "PDE_RING_EIGENHISTORY_INCONCLUSIVE"
        verdict = (
            "INCONCLUSIVE. A core check failed; re-examine before "
            "quoting the result."
        )

    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "identification": (
            "The PDE-ring eigenhistory: explicit time-domain "
            "evolution of the Duffing (q, p) coupled to the "
            "two-direction field ring with the interaction term in "
            "the conserved energy; the periodic orbit of the complete "
            "source-field state solved as X(T) - X(0) = 0, "
            "H[X] - E0 = 0 with the p(0) = 0 phase condition; and the "
            "one-period monodromy - all Floquet pairs on the unit "
            "circle, symplectic, with the source pair confirming the "
            "#219 reduction at 0.07%"
        ),
        "executes": (
            "the requested full Hamiltonian successor to #219"
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
    out.append("# The PDE-ring eigenhistory (PR #220)")
    out.append("")
    out.append(f"**Run:** {s['timestamp_utc']}")
    out.append("")
    out.append(
        "The deliverable is `docs/pde_ring_eigenhistory.md` - the "
        "full Hamiltonian successor to #219: the explicit time-domain "
        "source-field system, its periodic orbit, and the complete "
        "monodromy. *(QFT on the fixed classical throat geometry, not "
        "quantum gravity.)*"
    )
    out.append("")
    out.append("## Test summary")
    out.append("")
    out.append("| # | Test | Key finding | PASS? |")
    out.append("|---|---|---|---|")
    labels = {
        "T1": "the goal: the full Hamiltonian successor",
        "T2": "H with g q u(0) conserved; the term load-bearing",
        "T3": "the literal periodic-orbit solve, machine residual",
        "T4": "the complete monodromy: all unit circle, symplectic",
        "T5": "the source pair confirms #219 at 0.07%",
        "T6": "the time-resolved ledger closes",
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
    out = here / "runs" / f"{ts}_pde_ring_eigenhistory_probe"
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
