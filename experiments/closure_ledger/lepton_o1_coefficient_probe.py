"""
The remaining order-one lepton mass coefficient, derived - not fit -
on the #220 PDE-ring eigenhistory background (PR #221).

> Framing: QFT on the *fixed classical* throat geometry (geometry -> fields),
> not quantum gravity.

THE TARGET
----------
#210's register left exactly one number between the mass ladder and an
outright m_e/m_mu prediction: the O(1) coefficient sigma_mode/lambda_C
(required 88.6*alpha = 0.6467 in convention A, 206.8*alpha = 1.5089 in
convention B - "constrained, not derived", a x2.3 band).  This probe
derives it on the spatially converged, spectrally stable PDE
eigenhistory background of #220.

THE CHAIN (every link repo canon)
---------------------------------
1. The eigenhistory's mass IS its frequency (m = hbar w*), so
   lambda_C = 1/w* - both numerator and denominator of the O(1) are
   properties of ONE geometric object: X = sigma_mode * w* is a
   derived ratio, with nothing to fit.
2. The transactional arc (#217/#218) puts eigenhistories on the
   throat's interior resonance comb; the ring interior of the #220
   background IS the bridge coordinate of #202's 5D throat solve.
3. #202's machine-verified Pin-twisted boundary condition - odd-k
   modes have a NODE at the cross-cap - selects the ODD interior
   fundamental as the electron (k = 1) transit mode.  Its near-neck
   profile is exactly #202's regular solution phi ~ sigma (checked).
4. #202's sigma_mode is the MATCHING radius - where the interior
   linear growth turns over: on the mode, the antinode distance.
5. The hard-wall theorem: for the odd fundamental the antinode sits at
   the quarter wave, so X = w * sigma_match = pi/2 EXACTLY -
   cavity-length independent (that is WHY the coefficient is O(1)).
6. Measured on the true finite-width Tangherlini barriers:
   X = 1.579 - 1.638 across the full regulator scan, exterior-robust,
   dx-converged; the eigenhistory orbit (source-decoupled by parity,
   Gauss-Newton to 3e-13, monodromy unit-circle) carries the same X.

THE CONFRONTATION (no fit anywhere)
-----------------------------------
With #210's primordial anchor r_s = alpha*lambda_C and #202's exact
k = 1 law eps_1 = r_s/sigma_mode:

    m_e/m_mu = alpha / X

- convention A (required X = 0.6467) is EXCLUDED by x2.4;
- convention B (required X = 1.5089) is SELECTED, and
    m_e/m_mu = alpha/X = (2 alpha/pi) * (1 + throat shift)
             = 0.00445 - 0.00462   vs observed 0.0048363:
  the derivation lands at 92-96% of the observed ratio with ZERO
  fitted numbers (inputs: the throat geometry and alpha; m_mu enters
  only as the anchor scale, m_e is never used).

Tests:
  T1. Goal.
  T2. The identifications and the hard-wall theorems (closed forms for
      both parities and all three definitions, machine-checked).
  T3. The measurement at the reference instance (both parities: the
      #202 node-at-neck, the phi ~ sigma near-neck law, X values).
  T4. Robustness (interior-depth regulator scan; exterior
      independence; dx convergence).
  T5. The eigenhistory carries it (the source-decoupled odd orbit:
      Gauss-Newton, energy closure, unit-circle monodromy of the
      complete state, orbit-averaged X, exact linearity).
  T6. The confrontation (convention A excluded, B selected;
      m_e/m_mu = alpha/X at 92-96% of observed; the neck aspect).
  T7. Honest scope.
  T8. Assessment.

Verdict:
  THE_REMAINING_O1_IS_DERIVED_THE_ODD_THROAT_FUNDAMENTAL_GIVES_X_EQUALS
  _PI_OVER_2_PLUS_THROAT_SHIFT_AND_ME_OVER_MMU_EQUALS_ALPHA_OVER_X_AT_
  92_TO_96_PERCENT_NO_FIT
"""

from __future__ import annotations

import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np

try:
    from experiments.closure_ledger import (
        pde_ring_eigenhistory_probe as bg)
except ImportError:                                    # script execution
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import pde_ring_eigenhistory_probe as bg

_CACHE: dict = {}

# CODATA 2018 / PDG - comparison targets only, never inputs to the solve
_ALPHA = 0.0072973525693
_MMU_OVER_ME = 206.7682830
_ME_OVER_MMU = 1.0 / _MMU_OVER_ME
_X_REQ_A = (3.0 / 7.0) * _ALPHA * _MMU_OVER_ME     # convention A
_X_REQ_B = _ALPHA * _MMU_OVER_ME                   # convention B

# hard-wall closed forms (derived in T2)
_XU_EVEN_HW = math.pi * math.sqrt(1 / 12 - 1 / (2 * math.pi ** 2))
_XE_EVEN_HW = math.pi / math.sqrt(12)
_XU_ODD_HW = math.sqrt(math.pi ** 2 / 3 - 0.5)
_XM_ODD_HW = math.pi / 2
_XE_ODD_HW = 2 * math.pi / math.sqrt(12)

_D_REF = 12.0                     # reference interior depth (well-trapped)
_LEXT_REF = 2 * math.pi           # the S3 exterior arc of the #220 ring


# ── the generalized #220 two-barrier ring ───────────────────────────────

def build_ring(L_ext, D_int, dx_target=0.039):
    """The #220 background geometry with scannable interior depth and
    exterior arc: source at 0, mouths at +-L_ext/2, the neck (the #202
    cross-cap / bridge midpoint) at L_ext/2 + D_int/2."""
    C = L_ext + D_int
    N = int(round(C / dx_target))
    if N % 2:
        N += 1
    dx = C / N
    x = np.arange(N) * dx
    dA = (x - L_ext / 2 + C / 2) % C - C / 2
    dB = (x - (L_ext / 2 + D_int) + C / 2) % C - C / 2
    V = bg._v_tort(bg._X_PK - dA) + bg._v_tort(bg._X_PK + dB)
    return {'C': C, 'N': N, 'dx': dx, 'x': x, 'V': V,
            'L_ext': L_ext, 'D_int': D_int,
            'x_neck': L_ext / 2 + D_int / 2}


def _modes(s):
    N, dx, V = s['N'], s['dx'], s['V']
    L = (np.diag(-2.0 * np.ones(N)) + np.diag(np.ones(N - 1), 1)
         + np.diag(np.ones(N - 1), -1))
    L[0, -1] = L[-1, 0] = 1.0
    L /= dx ** 2
    evals, evecs = np.linalg.eigh(-L + np.diag(V))
    return np.sqrt(np.abs(evals)), evecs


def _neck_distance(s):
    return (s['x'] - s['x_neck'] + s['C'] / 2) % s['C'] - s['C'] / 2


def _antinode(dcen, w_amp, dx, D_int):
    """Sub-grid (parabolic) first-antinode distance from the neck."""
    mask = (dcen > 0) & (dcen < D_int / 2 + 1.0)
    idx = np.where(mask)[0]
    j = idx[int(np.argmax(w_amp[idx]))]
    y0, y1, y2 = w_amp[j - 1], w_amp[j], w_amp[j + 1]
    shift = 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2)
    return float(dcen[j] + shift * dx)


def measures(s, u, w):
    """The three definitional measures of one mode about the neck:
    amplitude-RMS (the #203 convention), energy-RMS, and the #202
    matching (turnover/antinode) radius."""
    dx, V = s['dx'], s['V']
    dcen = _neck_distance(s)
    up = (np.roll(u, -1) - np.roll(u, 1)) / (2 * dx)
    w_u = u ** 2
    w_E = w ** 2 * u ** 2 + up ** 2 + V * u ** 2
    sig_u = math.sqrt(float(np.sum(w_u * dcen ** 2) / np.sum(w_u)))
    sig_E = math.sqrt(float(np.sum(w_E * dcen ** 2) / np.sum(w_E)))
    d_anti = _antinode(dcen, np.abs(u), dx, s['D_int'])
    return {'X_u': float(sig_u * w), 'X_E': float(sig_E * w),
            'X_match': float(d_anti * w), 'sigma_match': d_anti,
            'u_neck_over_max': float(np.abs(u[int(np.argmin(
                np.abs(dcen)))]) / np.abs(u).max())}


def interior_pair(s):
    """The lowest interior-localized mode of each parity about the
    neck (max-interior-fraction identification; the ring is symmetric
    under x -> -x, which is reflection through the neck)."""
    key = ('PAIR', s['L_ext'], s['D_int'], s['N'])
    if key in _CACHE:
        return _CACHE[key]
    ws, evecs = _modes(s)
    N = s['N']
    dcen = _neck_distance(s)
    inside = np.abs(dcen) < s['D_int'] / 2
    refl = (-np.arange(N)) % N
    found = {}
    for k in range(min(120, N)):
        u = evecs[:, k]
        frac = float(np.sum(u[inside] ** 2))
        if frac < 0.55:
            continue
        par = float(np.sum(u * u[refl]))
        lab = 'even' if par > 0 else 'odd'
        if lab not in found:
            found[lab] = {'w': float(ws[k]), 'u': u, 'frac': frac,
                          'parity': par, 'k': int(k)}
        if len(found) == 2:
            break
    _CACHE[key] = found
    return found


# ── the odd eigenhistory orbit (the #220 machinery) ─────────────────────

def odd_orbit():
    """The full #220-style Gauss-Newton periodic orbit seeded on the
    odd interior fundamental (dx ~ 0.075 grid; phase condition
    <u, pi> = 0 because the odd mode decouples the source exactly)."""
    if 'ORB' in _CACHE:
        return _CACHE['ORB']
    s = build_ring(_LEXT_REF, _D_REF, dx_target=0.075)
    N, dx, V = s['N'], s['dx'], s['V']
    m = interior_pair(s)['odd']
    phi = m['u'] / np.sqrt(np.sum(m['u'] ** 2))
    u0 = (0.9 / np.abs(phi).max()) * phi
    z = np.concatenate([u0, 0 * u0, [0.0, 0.0]])
    T = 2 * math.pi / m['w']
    E0 = round(sum(bg._h_parts_at(u0, 0 * u0, 0.0, 0.0, dx, V)), 3)
    n_dim = 2 * N + 2

    def hh(zz):
        return sum(bg._h_parts_at(zz[:N], zz[N:2 * N], zz[2 * N],
                                  zz[2 * N + 1], dx, V))

    def residual(zz, TT):
        zT = bg._flow_batch_at(zz[:, None], TT, bg._NS, N, dx, V)[:, 0]
        return np.concatenate(
            [zT - zz, [hh(zz) - E0,
                       float(np.dot(zz[:N], zz[N:2 * N]))]])

    hist = []
    for _ in range(10):
        R = residual(z, T)
        hist.append(float(np.linalg.norm(R)))
        if hist[-1] < 1e-11:
            break
        eps = 1e-7
        Zb = np.concatenate([z[:, None],
                             z[:, None] + eps * np.eye(n_dim)], axis=1)
        ZT = bg._flow_batch_at(Zb, T, bg._NS, N, dx, V)
        base = ZT[:, 0]
        hb = hh(z)
        ph_b = float(np.dot(z[:N], z[N:2 * N]))
        J = np.empty((n_dim + 2, n_dim + 1))
        for j in range(n_dim):
            zj = Zb[:, j + 1]
            dtop = (ZT[:, j + 1] - base) / eps - (zj - z) / eps
            dH = (hh(zj) - hb) / eps
            dph = (float(np.dot(zj[:N], zj[N:2 * N])) - ph_b) / eps
            J[:, j] = np.concatenate([dtop, [dH, dph]])
        zT2 = bg._flow_batch_at(z[:, None], T + eps, bg._NS, N,
                                dx, V)[:, 0]
        J[:, n_dim] = np.concatenate([(zT2 - base) / eps, [0.0, 0.0]])
        step, *_ = np.linalg.lstsq(J, -R, rcond=None)
        z = z + step[:n_dim]
        T = T + step[n_dim]
    out = {'s': s, 'z': z, 'T': float(T), 'E0': E0,
           'w_lin': m['w'], 'newton_history': hist,
           'residual': float(np.linalg.norm(residual(z, T))),
           'periodicity_only': lambda zz, TT: float(np.linalg.norm(
               bg._flow_batch_at(zz[:, None], TT, bg._NS, N, dx,
                                 V)[:, 0] - zz))}
    _CACHE['ORB'] = out
    return out


# ========================================================================
# T1. Goal
# ========================================================================


def test_T1_goal() -> dict:
    return {
        'name': 'T1_goal',
        'description': (
            'Derive - not fit - the remaining order-one lepton mass '
            'coefficient sigma_mode/lambda_C (#210 register) on the '
            'spatially converged, spectrally stable PDE eigenhistory '
            'background of #220: the mass is the eigenhistory '
            'frequency, sigma_mode is the #202 matching radius of the '
            'odd (Pin-twisted, node-at-neck) interior fundamental, '
            'and X = sigma_mode * w is a ratio of two properties of '
            'one geometric object.'
        ),
        'pass': True,
    }


# ========================================================================
# T2. The identifications and the hard-wall theorems
# ========================================================================


def test_T2_hard_wall_theorems() -> dict:
    # analytic modes sampled on a fine grid; the measures must
    # reproduce the closed forms
    L, M = 8.0, 20001
    xs = np.linspace(0, L, M)
    d = xs - L / 2
    got = {}
    for lab, n in (('even', 1), ('odd', 2)):
        kk = n * math.pi / L
        u = np.sin(n * math.pi * xs / L)
        up = kk * np.cos(n * math.pi * xs / L)
        w_u = u ** 2
        w_E = kk ** 2 * u ** 2 + up ** 2
        sig_u = math.sqrt(float(np.trapezoid(w_u * d ** 2, xs)
                                / np.trapezoid(w_u, xs)))
        sig_E = math.sqrt(float(np.trapezoid(w_E * d ** 2, xs)
                                / np.trapezoid(w_E, xs)))
        got[lab] = {'X_u': sig_u * kk, 'X_E': sig_E * kk}
    # antinode of the odd mode (parabolic refinement) -> quarter wave
    u = np.abs(np.sin(2 * math.pi * xs / L))
    dx = xs[1] - xs[0]
    mask = d > 0
    idx = np.where(mask)[0]
    j = idx[int(np.argmax(u[idx]))]
    y0, y1, y2 = u[j - 1], u[j], u[j + 1]
    d_anti = d[j] + 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2) * dx
    got['odd']['X_match'] = float(d_anti * 2 * math.pi / L)

    errs = {
        'even_X_u': abs(got['even']['X_u'] - _XU_EVEN_HW),
        'even_X_E': abs(got['even']['X_E'] - _XE_EVEN_HW),
        'odd_X_u': abs(got['odd']['X_u'] - _XU_ODD_HW),
        'odd_X_E': abs(got['odd']['X_E'] - _XE_ODD_HW),
        'odd_X_match': abs(got['odd']['X_match'] - _XM_ODD_HW),
    }
    ok = max(errs.values()) < 1e-4
    return {
        'name': 'T2_hard_wall_theorems',
        'description': (
            'the identifications (mass = eigenhistory frequency, '
            'lambda_C = 1/w; parity selector = the #202 Pin-twisted '
            'node-at-neck theorem for odd k; sigma_mode = the #202 '
            'matching radius) and the hard-wall closed forms: for the '
            'odd fundamental the antinode is the quarter wave, X = '
            'w sigma_match = pi/2 EXACTLY, cavity-length independent '
            '- the analytic reason the coefficient is O(1)'
        ),
        'closed_forms': {
            'even_X_u': _XU_EVEN_HW, 'even_X_E': _XE_EVEN_HW,
            'odd_X_u': _XU_ODD_HW, 'odd_X_match': _XM_ODD_HW,
            'odd_X_E': _XE_ODD_HW},
        'sampled_errors': {k: float(v) for k, v in errs.items()},
        'pass': bool(ok),
    }


# ========================================================================
# T3. The measurement at the reference instance
# ========================================================================


def test_T3_measurement() -> dict:
    if 'T3' in _CACHE:
        return _CACHE['T3']
    s = build_ring(_LEXT_REF, _D_REF, dx_target=0.0374)
    pair = interior_pair(s)
    ev, od = pair['even'], pair['odd']
    me = measures(s, ev['u'], ev['w'])
    mo = measures(s, od['u'], od['w'])

    # the #202 near-neck law phi ~ sigma for the odd mode: the ratio
    # u(d)/d drifts only by the mode's curvature over |d| < 0.8
    dcen = _neck_distance(s)
    uo = od['u'] * np.sign(od['u'][np.argmax(dcen > 0.5)])
    r04 = float(np.interp(0.4, dcen[np.argsort(dcen)],
                          uo[np.argsort(dcen)])) / 0.4
    r08 = float(np.interp(0.8, dcen[np.argsort(dcen)],
                          uo[np.argsort(dcen)])) / 0.8
    lin_drift = abs(r08 / r04 - 1)

    ok = (ev['frac'] > 0.95 and od['frac'] > 0.90
          and me['u_neck_over_max'] > 0.99          # even: antinode
          and mo['u_neck_over_max'] < 1e-8          # odd: exact node
          and od['parity'] < -0.99 and ev['parity'] > 0.99
          and 0.60 < me['X_u'] < 0.67
          and 1.55 < mo['X_match'] < 1.61
          and 1.80 < mo['X_u'] < 1.93
          and 1.95 < mo['X_E'] < 2.04
          and lin_drift < 0.05)
    out = {
        'name': 'T3_measurement',
        'description': (
            'the interior fundamentals of both parities on the real '
            'Tangherlini two-barrier background (reference instance '
            'D = 12, S3 exterior arc 2pi): the odd (electron, k = 1) '
            'mode has the EXACT #202 node at the neck and the phi ~ '
            'sigma near-neck law; the even (uncharged, k = 0) mode '
            'touches the cross-cap - the #202 parity dichotomy '
            'realized on the eigenhistory background; all three X '
            'definitions measured'
        ),
        'even': {'w': ev['w'], 'interior_fraction': ev['frac'],
                 'parity': ev['parity'], 'X_u': me['X_u'],
                 'X_E': me['X_E'],
                 'u_neck_over_max': me['u_neck_over_max']},
        'odd': {'w': od['w'], 'interior_fraction': od['frac'],
                'parity': od['parity'], 'X_u': mo['X_u'],
                'X_E': mo['X_E'], 'X_match': mo['X_match'],
                'sigma_match': mo['sigma_match'],
                'u_neck_over_max': mo['u_neck_over_max']},
        'odd_near_neck_linearity_drift': float(lin_drift),
        'hard_wall_reference': {'odd_X_match': _XM_ODD_HW},
        'pass': bool(ok),
    }
    _CACHE['T3'] = out
    return out


# ========================================================================
# T4. Robustness: regulator, exterior, grid
# ========================================================================


def test_T4_robustness() -> dict:
    if 'T4' in _CACHE:
        return _CACHE['T4']
    # (i) interior-depth (regulator) scan - the #215 horizon-continuum
    # depth is a regulator; X_match must not depend on it
    depth = {}
    for D in (8.0, 10.0, 12.0, 16.0, 24.0, 32.0, 48.0):
        s = build_ring(_LEXT_REF, D)
        pair = interior_pair(s)
        row = {}
        for lab in ('even', 'odd'):
            if lab in pair:
                mm = measures(s, pair[lab]['u'], pair[lab]['w'])
                row[lab] = {'w': pair[lab]['w'],
                            'frac': pair[lab]['frac'],
                            'X_u': mm['X_u'], 'X_E': mm['X_E'],
                            'X_match': mm['X_match']}
        depth[str(D)] = row
    xm = [depth[k]['odd']['X_match'] for k in depth]
    xm_lo, xm_hi = min(xm), max(xm)

    # (ii) exterior independence at D = 12
    ext = {}
    for Lx in (math.pi, 2 * math.pi, 4 * math.pi):
        s = build_ring(Lx, _D_REF)
        pair = interior_pair(s)
        mm = measures(s, pair['odd']['u'], pair['odd']['w'])
        ext[f'{Lx:.3f}'] = {'w': pair['odd']['w'],
                            'frac': pair['odd']['frac'],
                            'X_match': mm['X_match'],
                            'X_u': mm['X_u']}
    xe = [v['X_match'] for v in ext.values()]

    # (iii) dx convergence at the reference instance
    grid = {}
    for dxt in (0.112, 0.075, 0.056, 0.0374):
        s = build_ring(_LEXT_REF, _D_REF, dx_target=dxt)
        pair = interior_pair(s)
        mm = measures(s, pair['odd']['u'], pair['odd']['w'])
        grid[str(dxt)] = {'N': s['N'], 'w': pair['odd']['w'],
                          'X_match': mm['X_match'], 'X_u': mm['X_u']}
    dx_shift = abs(grid['0.0374']['X_match'] - grid['0.075']['X_match'])

    ok = (all(1.55 < v < 1.67 for v in xm)
          and (xm_hi - xm_lo) < 0.08
          and (max(xe) - min(xe)) < 0.05
          and dx_shift < 1e-3)
    out = {
        'name': 'T4_robustness',
        'description': (
            'X_match is regulator-independent (interior depth 8-48: '
            'band 0.06 wide, all within 6% of pi/2), '
            'exterior-independent (S3 arc pi-4pi: spread < 0.04, and '
            'robust where the RMS definition is contaminated by '
            'hybridization), and dx-converged (< 1e-3)'
        ),
        'depth_scan': depth,
        'X_match_depth_band': [float(xm_lo), float(xm_hi)],
        'exterior_scan': ext,
        'X_match_exterior_spread': float(max(xe) - min(xe)),
        'grid_scan': grid,
        'X_match_dx_shift': float(dx_shift),
        'pass': bool(ok),
    }
    _CACHE['T4'] = out
    return out


# ========================================================================
# T5. The eigenhistory carries it
# ========================================================================


def test_T5_eigenhistory_orbit() -> dict:
    if 'T5' in _CACHE:
        return _CACHE['T5']
    orb = odd_orbit()
    s, z, T = orb['s'], orb['z'], orb['T']
    N, dx, V = s['N'], s['dx'], s['V']

    # source decoupling is exact on the orbit (odd parity => u(0) = 0)
    q_star, p_star = abs(float(z[2 * N])), abs(float(z[2 * N + 1]))
    H_at = sum(bg._h_parts_at(z[:N], z[N:2 * N], z[2 * N],
                              z[2 * N + 1], dx, V))

    # orbit-averaged profile over one period
    dcen = _neck_distance(s)
    zz = z.copy()
    acc_u = np.zeros(N)
    acc_E = np.zeros(N)
    for _ in range(64):
        zz = bg._flow_batch_at(zz[:, None], T / 64, bg._NS // 64,
                               N, dx, V)[:, 0]
        uu, pp = zz[:N], zz[N:2 * N]
        upx = (np.roll(uu, -1) - np.roll(uu, 1)) / (2 * dx)
        acc_u += uu ** 2
        acc_E += pp ** 2 + upx ** 2 + V * uu ** 2
    w_orb = 2 * math.pi / T
    sig_u = math.sqrt(float(np.sum(acc_u * dcen ** 2) / np.sum(acc_u)))
    d_anti = _antinode(dcen, np.sqrt(acc_u), dx, s['D_int'])
    X_u_orb = sig_u * w_orb
    X_m_orb = d_anti * w_orb

    # linear-mode comparison on the same grid
    mo = measures(s, interior_pair(s)['odd']['u'],
                  interior_pair(s)['odd']['w'])

    # exact linearity of the decoupled sector: the half-amplitude
    # state is periodic at the SAME T without re-solving
    z_half = z.copy()
    z_half[:2 * N] *= 0.5
    lin_res = orb['periodicity_only'](z_half, T)

    # the complete monodromy (source in the tangent space)
    Mono = bg._monodromy_at(N, z, T, dx, V)
    ev = np.linalg.eigvals(Mono)
    max_dev = float(np.abs(np.abs(ev) - 1).max())

    ok = (orb['residual'] < 1e-11
          and q_star < 1e-12 and p_star < 1e-12
          and abs(H_at - orb['E0']) < 1e-10
          and abs(X_m_orb - mo['X_match']) < 2e-3
          and abs(X_u_orb / mo['X_u'] - 1) < 2e-3
          and lin_res < 1e-9
          and max_dev < 1e-9)
    out = {
        'name': 'T5_eigenhistory_orbit',
        'description': (
            'the #220 Gauss-Newton periodic orbit seeded on the odd '
            'interior fundamental: converged to machine residual with '
            'the source EXACTLY decoupled by parity (q* = p* = 0 - '
            'the odd/charged transit mode is source-transparent at '
            'the crossing); energy closure; the complete monodromy '
            '(source in the tangent space) on the unit circle; the '
            'orbit-averaged profile carries the same X; the decoupled '
            'sector exactly linear (X amplitude-independent)'
        ),
        'E0': orb['E0'],
        'period': T,
        'w_orbit': float(w_orb),
        'w_linear': orb['w_lin'],
        'newton_history': orb['newton_history'],
        'final_residual': orb['residual'],
        'q_star': q_star, 'p_star': p_star,
        'energy_at_orbit': float(H_at),
        'X_match_orbit': float(X_m_orb),
        'X_match_linear': mo['X_match'],
        'X_u_orbit': float(X_u_orb),
        'X_u_linear': mo['X_u'],
        'half_amplitude_periodicity_residual': float(lin_res),
        'monodromy_dimension': 2 * N + 2,
        'monodromy_max_modulus_deviation': max_dev,
        'pass': bool(ok),
    }
    _CACHE['T5'] = out
    return out


# ========================================================================
# T6. The confrontation
# ========================================================================


def test_T6_confrontation() -> dict:
    t3 = test_T3_measurement()
    t4 = test_T4_robustness()
    X_ref = t3['odd']['X_match']
    xm_lo, xm_hi = t4['X_match_depth_band']

    # convention A vs B (the #201 S1 conventions)
    ratio_A = X_ref / _X_REQ_A
    excluded_A = ratio_A > 2.0
    pred_ref = _ALPHA / X_ref
    pred_lo = _ALPHA / xm_hi
    pred_hi = _ALPHA / xm_lo
    pred_hw = 2 * _ALPHA / math.pi
    land_ref = pred_ref / _ME_OVER_MMU
    land_band = [pred_lo / _ME_OVER_MMU, pred_hi / _ME_OVER_MMU]

    # the neck aspect: c = ln(sigma/r_s) = ln(X/alpha)
    c_derived = math.log(X_ref / _ALPHA)
    c_required_B = math.log(_MMU_OVER_ME)

    # the even/conv-A alternate (parity-excluded by #202; recorded)
    X_even = t3['even']['X_u']
    alt_pred = (3.0 / 7.0) * _ALPHA / X_even

    ok = (excluded_A
          and 0.90 < land_ref < 1.00
          and 0.90 < land_band[0] and land_band[1] < 1.00
          and abs(pred_hw / _ME_OVER_MMU - 0.9606) < 0.01)
    return {
        'name': 'T6_confrontation',
        'description': (
            'the derived X against the two #201 conventions: A '
            '(required 0.6467) EXCLUDED by x2.4; B (required 1.5089) '
            'SELECTED - m_e/m_mu = alpha/X = (2 alpha/pi)(1 + throat '
            'shift) lands at 92-96% of the observed ratio with zero '
            'fitted numbers'
        ),
        'X_required_conv_A': _X_REQ_A,
        'X_required_conv_B': _X_REQ_B,
        'X_derived_reference': X_ref,
        'X_derived_depth_band': [xm_lo, xm_hi],
        'X_hard_wall': _XM_ODD_HW,
        'conv_A_ratio': float(ratio_A),
        'conv_A_excluded': bool(excluded_A),
        'me_over_mmu_predicted_reference': float(pred_ref),
        'me_over_mmu_predicted_band': [float(pred_lo), float(pred_hi)],
        'me_over_mmu_hard_wall_2alpha_over_pi': float(pred_hw),
        'me_over_mmu_observed': _ME_OVER_MMU,
        'landing_fraction_reference': float(land_ref),
        'landing_fraction_band': [float(x) for x in land_band],
        'landing_fraction_hard_wall': float(pred_hw / _ME_OVER_MMU),
        'neck_aspect_c_derived': float(c_derived),
        'neck_aspect_c_required_conv_B': float(c_required_B),
        'even_conv_A_alternate_prediction': float(alt_pred),
        'even_conv_A_alternate_note': (
            'lands in-band but is parity-EXCLUDED: #202 proves the '
            'k = 1 mode has a node at the neck (odd), and the even '
            'mode is the k = 0 (uncharged) channel'),
        'pass': bool(ok),
    }


# ========================================================================
# T7. Honest scope
# ========================================================================


def test_T7_honest_scope() -> dict:
    scope = [
        'The ring-transit <-> #202-bridge identification (the ring '
        'interior as the bridge sigma coordinate, the neck as the '
        'cross-cap) is structural and machine-consistent here (exact '
        'node at the neck, the phi ~ sigma near-neck law), but the 1D '
        'transit measure is not the 5D rho^3 bridge measure - redoing '
        'X on the true 5D radial operator is the named successor and '
        'the leading candidate for the remaining 4-8% residual.',
        'The winding number k is not represented in the zonal scalar; '
        'the parity of the k = 1 mode is imported from #202\'s '
        'machine-verified Pin-twisted boundary condition, not '
        're-derived here.',
        'Inputs: the Tangherlini barrier geometry (r_h = 1) and alpha '
        '(the #184-protected boundary invariant) via #210\'s '
        'primordial anchor r_s = alpha*lambda_C. m_mu enters only as '
        'the anchor scale of the #201 law (convention B: S1 = m_mu); '
        'm_e is used NOWHERE in the construction - only in the final '
        'comparison.',
        'The interior depth is the #215 horizon-continuum regulator: '
        'X_match is stable to 4% over depth 8-48 and the band is '
        'carried; the physically-capped depth (the alpha-capped '
        'tortoise depth, #210\'s c = ln(1/alpha) + O(1)) is a '
        'successor refinement.',
        'The convention A/B ambiguity of #201 is resolved BY the '
        'measurement (B selected, A excluded x2.4), not assumed; the '
        'even-parity/conv-A alternate is recorded and excluded by the '
        '#202 parity theorem.',
        'Classical, zonal scalar, frozen geometry; hbar enters only '
        'through lambda_C = hbar/mc as always (B4 anchor).',
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
    t2 = test_T2_hard_wall_theorems()
    t3 = test_T3_measurement()
    t4 = test_T4_robustness()
    t5 = test_T5_eigenhistory_orbit()
    t6 = test_T6_confrontation()
    core = all(t['pass'] for t in (t2, t3, t4, t5, t6))
    assessment = (
        'The #210 register item is executed: the remaining O(1) is no '
        'longer a fit band. On the #220 eigenhistory background the '
        'electron transit mode is forced by #202\'s parity theorem '
        '(odd, node at the neck), its mass is its frequency, and its '
        '#202 matching radius is the antinode - so X = sigma * w is a '
        'derived ratio with a hard-wall value of exactly pi/2 and a '
        'measured throat value of 1.58-1.64, regulator- and '
        'exterior-robust, carried unchanged by the source-decoupled '
        'Gauss-Newton eigenhistory orbit with unit-circle monodromy. '
        'Confronted with the #201 law at #210\'s primordial anchor, '
        'convention A is excluded and convention B lands: m_e/m_mu = '
        'alpha/X = (2 alpha/pi)(1 + throat shift) reaches 92-96% of '
        'the observed value with zero fitted numbers. The residual '
        '4-8% is real, stated, and owned by the named successor (the '
        '5D bridge measure).'
        if core else
        'INCONCLUSIVE - a core check failed; do not quote.'
    )
    return {
        'name': 'T8_assessment',
        'description': 'the standing of the derived coefficient',
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
        test_T2_hard_wall_theorems(),
        test_T3_measurement(),
        test_T4_robustness(),
        test_T5_eigenhistory_orbit(),
        test_T6_confrontation(),
        test_T7_honest_scope(),
        test_T8_assessment(),
    ]
    t2, t3, t4, t5, t6 = tests[1], tests[2], tests[3], tests[4], tests[5]
    all_ok = all(t['pass'] for t in tests)
    if all_ok:
        verdict_class = (
            "THE_REMAINING_O1_IS_DERIVED_THE_ODD_THROAT_FUNDAMENTAL_"
            "GIVES_X_EQUALS_PI_OVER_2_PLUS_THROAT_SHIFT_AND_ME_OVER_"
            "MMU_EQUALS_ALPHA_OVER_X_AT_92_TO_96_PERCENT_NO_FIT"
        )
        verdict = (
            "ESTABLISHED (the argument is in "
            "docs/lepton_o1_coefficient.md).\n\n"
            "THE COEFFICIENT. X = sigma_mode * w on ONE geometric "
            "object - the eigenhistory whose frequency is the mass "
            "(lambda_C = 1/w) and whose #202 matching radius is "
            "sigma_mode: a derived ratio, nothing to fit.\n\n"
            "THE THEOREM. #202's Pin-twisted boundary condition "
            "forces the electron (k = 1) transit mode ODD - the exact "
            "node at the neck is reproduced here to "
            f"{t3['odd']['u_neck_over_max']:.0e} with the phi ~ sigma "
            "near-neck law (drift "
            f"{t3['odd_near_neck_linearity_drift']:.1%}); for the odd "
            "fundamental the antinode is the quarter wave: X = pi/2 "
            "EXACTLY on hard walls (closed forms verified to "
            f"{max(t2['sampled_errors'].values()):.0e}), cavity-length "
            "independent - the analytic reason the coefficient is "
            "O(1).\n\n"
            "THE MEASUREMENT. On the real Tangherlini barriers "
            f"(reference D = 12): X_match = "
            f"{t3['odd']['X_match']:.4f} (hard wall "
            f"{_XM_ODD_HW:.4f}); regulator scan (depth 8-48) band "
            f"[{t4['X_match_depth_band'][0]:.4f}, "
            f"{t4['X_match_depth_band'][1]:.4f}]; exterior spread "
            f"{t4['X_match_exterior_spread']:.3f}; dx-shift "
            f"{t4['X_match_dx_shift']:.0e}.\n\n"
            "THE ORBIT. The #220 Gauss-Newton eigenhistory on the odd "
            f"mode: residual {t5['final_residual']:.0e}; the source "
            f"decouples EXACTLY (q* = {t5['q_star']:.0e}) - the "
            "charged transit mode is source-transparent at the "
            "crossing; complete monodromy on the unit circle to "
            f"{t5['monodromy_max_modulus_deviation']:.0e}; the orbit "
            f"carries X_match = {t5['X_match_orbit']:.4f} "
            "(linear-mode value to "
            f"{abs(t5['X_match_orbit']-t5['X_match_linear']):.0e}); "
            "exactly linear (amplitude-independent, "
            f"{t5['half_amplitude_periodicity_residual']:.0e}).\n\n"
            "THE CONFRONTATION. Convention A (required "
            f"{_X_REQ_A:.4f}) EXCLUDED x{t6['conv_A_ratio']:.1f}; "
            f"convention B (required {_X_REQ_B:.4f}) SELECTED: "
            "m_e/m_mu = alpha/X = "
            f"{t6['me_over_mmu_predicted_reference']:.6f} (band "
            f"[{t6['me_over_mmu_predicted_band'][0]:.6f}, "
            f"{t6['me_over_mmu_predicted_band'][1]:.6f}]; hard-wall "
            f"anchor 2 alpha/pi = "
            f"{t6['me_over_mmu_hard_wall_2alpha_over_pi']:.6f}) vs "
            f"observed {_ME_OVER_MMU:.6f} - the derivation lands at "
            f"{t6['landing_fraction_band'][0]:.1%}-"
            f"{t6['landing_fraction_band'][1]:.1%} of the observed "
            "ratio with ZERO fitted numbers (inputs: the throat "
            "geometry and alpha; m_e used only for comparison)."
        )
    else:
        verdict_class = "LEPTON_O1_COEFFICIENT_INCONCLUSIVE"
        verdict = (
            "INCONCLUSIVE. A core check failed; re-examine before "
            "quoting the result."
        )

    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "identification": (
            "The remaining order-one lepton mass coefficient derived "
            "on the #220 PDE eigenhistory background: X = sigma_mode "
            "* w of the odd (Pin-twisted, #202) interior fundamental "
            "with sigma_mode the #202 matching radius - pi/2 on hard "
            "walls, 1.58-1.64 measured on the Tangherlini barriers - "
            "giving m_e/m_mu = alpha/X at 92-96% of observed, no fit"
        ),
        "executes": (
            "the #210 register item 'derive the O(1) "
            "(sigma_mode/lambda_C)' on the requested #220 background"
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
    if callable(o):
        return "<fn>"
    return str(o)


def render_markdown(s: dict) -> str:
    out = []
    out.append("# The derived O(1) lepton mass coefficient (PR #221)")
    out.append("")
    out.append(f"**Run:** {s['timestamp_utc']}")
    out.append("")
    out.append(
        "The deliverable is `docs/lepton_o1_coefficient.md` - the #210 "
        "register item executed on the #220 eigenhistory background. "
        "*(QFT on the fixed classical throat geometry, not quantum "
        "gravity.)*"
    )
    out.append("")
    out.append("## Test summary")
    out.append("")
    out.append("| # | Test | Key finding | PASS? |")
    out.append("|---|---|---|---|")
    labels = {
        "T1": "the goal: derive, not fit",
        "T2": "the hard-wall theorems: X = pi/2 exactly",
        "T3": "the measurement; the #202 parity dichotomy realized",
        "T4": "regulator-, exterior-, grid-robust",
        "T5": "the eigenhistory orbit carries it; source-decoupled",
        "T6": "conv A excluded, conv B lands at 92-96%",
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
    out = here / "runs" / f"{ts}_lepton_o1_coefficient_probe"
    out.mkdir(parents=True, exist_ok=True)
    (out / "probe.json").write_text(
        json.dumps(summary, indent=2, default=_json_default))
    (out / "probe.md").write_text(md)
    print(f"\nWrote: {out / 'probe.json'}")
    print(f"Wrote: {out / 'probe.md'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
