"""
The dynamical absorber: S = S_field + S_absorber + S_coupling, and the
epsilon derived as a damping rate (PR #214).

> Framing: QFT on the *fixed classical* throat geometry (geometry -> fields),
> not quantum gravity.

THE DECISIVE SUCCESSOR TO #213
------------------------------
#213 derived the Feynman propagator from the frozen bulk using complete
absorption as a BOUNDARY CONDITION: no free remnants on complete
histories.  The decisive step is to make the absorber an actual degree
of freedom with its own action,

    S = S_field + S_absorber + S_coupling ,

    S_field    = sum_n int (udot_n^2 - w_n^2 u_n^2)/2 dt     (the tower
                 w_n = n/R, the conformal S^3 modes)
    S_absorber = sum_j int (qdot_j^2 - nu_j^2 q_j^2)/2 dt    (a bank of
                 oscillators - the throat's internal continuum; the
                 engine's MouthState.modes are exactly this d.o.f.,
                 which until now had no action)
    S_coupling = - sum_j g_j int q_j Phi(psi_a, t) dt ,
                 Phi(psi_a) = sum_n Y_n(psi_a) u_n ,
                 Y_n(psi) = sin(n psi)/(sin(psi) sqrt(2 pi^2)) ,

and show that ABSORPTION IS AN OUTCOME OF THE DYNAMICS, with the
i*epsilon now a COMPUTED function of the coupling:

1. EXACT ELIMINATION.  The action is Gaussian, so integrating out the
   absorber is exact: the field's effective resolvent is
   G_eff(W)^-1 = w0^2 - W^2 - Sigma(W), Sigma(W) = sum_j g_j^2/(nu_j^2 - W^2)
   - verified against direct matrix inversion of the coupled system.

2. THE DERIVED DAMPING.  The resolvent pole moves BELOW the real axis
   by a computed amount: Om* = w_tilde - i*gamma/2 with
   gamma = pi g^2 rho / (2 w0^2) (flat bank).  Three independent
   measurements agree: the complex pole, the live time-domain energy
   decay of the coupled system, and the golden-rule formula; and
   gamma scales exactly as g^2 (log-log slope 2).  #213's imposed
   epsilon is now eps_derived = w_tilde * gamma - the absorber's
   response rate, with g -> 0 recovering #213's kernel exactly.

3. THE GEOMETRY ADJUDICATES PLACEMENT (antipodal vs distributed).
   Coupling of a point absorber at psi_a to tower mode n is
   Y_n(psi_a); the per-mode damping rate is
   gamma_n = pi (g Y_n)^2 rho / (2 w_n^2):
     - ANTIPODE: Y_n(pi) = n(-1)^{n+1}/sqrt(2 pi^2) - the coupling
       GROWS as n and never vanishes, exactly cancelling the 1/w_n^2
       kinematic suppression: gamma_n = g^2/(4 pi), THE SAME RATE FOR
       EVERY MODE.  The antipode is the impedance-matched absorber -
       the #166 focusing caustic as absorber matching.
     - GENERIC POINT: Y_n oscillates with near-zeros (|sin 22| = 0.009
       at psi_a = 1) and gamma_n falls as 1/n^2 - long-lived remnants.
     - UNIFORM DISTRIBUTED: the exact selection rule
       int sin(n psi) sin(psi) dpsi = (pi/2) delta_{n1} - a uniformly
       distributed absorber couples ONLY to the tower ground mode and
       absorbs nothing else.
   Live adjudication (kicked multi-mode field, exact normal-mode
   evolution): residual field energy 0.014 (antipode, uniform across
   modes as the flat rate predicts) vs 0.95 (generic) vs 1.00
   (uniform).

4. RECURRENCE AND THE ORDER OF LIMITS.  A finite bank revives near
   T_rec = 2 pi / d(nu) (linear in N): on the closed bulk "epsilon
   finite" means Poincare recurrence; eps -> 0 is N -> infinity FIRST
   - the continuum limit of the throat's internal spectrum.

Tests:
  T1. Goal.
  T2. The action and the exact elimination (resolvent identity vs
      matrix inversion; absorber stability).
  T3. The derived damping (pole / live decay / golden rule agree;
      exact g^2 scaling; continuum Sigma verified).
  T4. The geometry adjudicates placement (flat antipodal rate law;
      selection rules; the live three-way adjudication).
  T5. Recurrence and the order of limits.
  T6. The #213 closure: eps_derived = w_tilde*gamma fills in the
      orderings' i*epsilon; the Delta^2 pole form agrees to O(gamma^2);
      g -> 0 recovers #213.
  T7. Honest scope.
  T8. Assessment.

Verdict:
  EPSILON_IS_THE_ABSORBER_RESPONSE_RATE_AND_THE_ANTIPODE_IS_THE_
  IMPEDANCE_MATCHED_ABSORBER_OF_THE_CLOSED_BULK
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np

_CACHE: dict = {}
_SQ = math.sqrt(2 * math.pi ** 2)

# ========================================================================
# SECTION A - the coupled Gaussian system and its exact elimination
# ========================================================================

_W0 = 3.0                 # reference tower mode for the single-mode model
_LO, _HI = 0.5, 6.0       # bank frequency span (single-mode model)


def flat_bank(n_bath: int, lo: float = _LO, hi: float = _HI):
    nu = np.linspace(lo, hi, n_bath)
    dnu = nu[1] - nu[0]
    return nu, dnu


def stiffness_matrix(w0: float, nu: np.ndarray, gj: np.ndarray) -> np.ndarray:
    n = len(nu)
    K = np.zeros((n + 1, n + 1))
    K[0, 0] = w0 ** 2
    K[np.arange(1, n + 1), np.arange(1, n + 1)] = nu ** 2
    K[0, 1:] = gj
    K[1:, 0] = gj
    return K


def sigma_discrete(z: complex, nu: np.ndarray, gj: np.ndarray) -> complex:
    return complex(np.sum(gj ** 2 / (nu ** 2 - z ** 2)))


def sigma_continuum(z: complex, kappa: float,
                    lo: float = _LO, hi: float = _HI) -> complex:
    """N -> infinity limit of the flat bank, g_j^2 = kappa * dnu:
    Sigma(z) = kappa * int_lo^hi dnu/(nu^2 - z^2)."""
    return (kappa / (2 * z)) * (np.log((hi - z) / (hi + z))
                                - np.log((lo - z) / (lo + z)))


def resolvent_pole(kappa: float, w0: float = _W0) -> complex:
    """Newton iteration for the complex pole of the continuum resolvent."""
    d = lambda z: w0 ** 2 - z ** 2 - sigma_continuum(z, kappa)
    om = complex(w0, -1e-3)
    for _ in range(100):
        h = 1e-8
        dd = (d(om + h) - d(om - h)) / (2 * h)
        om = om - d(om) / dd
    return om


def normal_mode_state(K: np.ndarray, x0: np.ndarray, v0: np.ndarray,
                      t: float):
    """Exact evolution of xddot = -K x (no time-stepping error)."""
    lam, V = np.linalg.eigh(K)
    rt = np.sqrt(lam)
    cx, cv = V.T @ x0, V.T @ v0
    x = V @ (np.cos(rt * t) * cx + np.sin(rt * t) / rt * cv)
    xd = V @ (-rt * np.sin(rt * t) * cx + np.cos(rt * t) * cv)
    return x, xd, lam


# ========================================================================
# T1. Goal
# ========================================================================


def test_T1_goal() -> dict:
    return {
        'name': 'T1_goal',
        'description': (
            'Promote the absorber from a boundary condition (#213) to '
            'a degree of freedom: S = S_field + S_absorber + '
            'S_coupling, all Gaussian, all classical.  Show absorption '
            'is an outcome of the dynamics, epsilon is the computed '
            'damping rate of the absorber response, and the geometry '
            'itself adjudicates where the absorber must live - the '
            'antipodal point is impedance-matched to every tower mode '
            'while a uniformly distributed absorber couples only to '
            'the ground mode.'
        ),
        'pass': True,
    }


# ========================================================================
# T2. The action and the exact elimination
# ========================================================================


def test_T2_exact_elimination() -> dict:
    w0 = _W0
    nu, dnu = flat_bank(400)
    kappa = 0.28
    gj = np.sqrt(kappa * dnu) * np.ones(len(nu))
    K = stiffness_matrix(w0, nu, gj)

    # (i) resolvent identity: [(K - W^2 I)^-1]_{00} = 1/(w0^2 - W^2 - Sigma)
    ident_err = 0.0
    for z in (2.0 + 0.3j, 3.5 - 0.7j, 0.9 + 0.1j, 5.1 - 0.2j):
        A = K - z ** 2 * np.eye(len(nu) + 1)
        direct = np.linalg.inv(A.astype(complex))[0, 0]
        formula = 1.0 / (w0 ** 2 - z ** 2 - sigma_discrete(z, nu, gj))
        ident_err = max(ident_err, abs(direct - formula) / abs(formula))

    # (ii) absorber stability: positive-definite stiffness (no runaway,
    # the coupled action is bounded below)
    min_eig = float(np.linalg.eigvalsh(K).min())

    # (iii) continuum Sigma is the N -> infinity limit of the bank
    z = 3.0 - 0.1j
    nu_l, dnu_l = flat_bank(200000)
    gj_l = np.sqrt(kappa * dnu_l) * np.ones(len(nu_l))
    cont_err = abs(sigma_discrete(z, nu_l, gj_l)
                   - sigma_continuum(z, kappa))

    ok = ident_err < 1e-10 and min_eig > 0 and cont_err < 1e-5
    return {
        'name': 'T2_exact_elimination',
        'description': (
            'the Gaussian absorber integrates out exactly: '
            'G_eff^-1 = w0^2 - W^2 - Sigma(W), verified against '
            'direct matrix inversion; the coupled action is stable'
        ),
        'resolvent_identity_error': float(ident_err),
        'min_stiffness_eigenvalue': min_eig,
        'continuum_sigma_error': float(cont_err),
        'pass': bool(ok),
    }


# ========================================================================
# T3. The derived damping
# ========================================================================


def test_T3_derived_damping() -> dict:
    if 'T3' in _CACHE:
        return _CACHE['T3']
    w0 = _W0
    gamma_set = 0.05                      # design point
    kappa = 2 * gamma_set * w0 ** 2 / math.pi

    # (a) the finite-bank complex pole (analytic continuation of the
    # discrete sum; bank dense enough that dnu << gamma)
    nu, dnu = flat_bank(1200)
    gj = np.sqrt(kappa * dnu) * np.ones(len(nu))
    d = lambda z: w0 ** 2 - z ** 2 - sigma_discrete(z, nu, gj)
    om = complex(w0, -0.02)
    for _ in range(80):
        h = 1e-7
        om = om - d(om) / ((d(om + h) - d(om - h)) / (2 * h))
    gamma_pole = -2 * om.imag

    # (b) live time-domain energy decay of the coupled system
    K = stiffness_matrix(w0, nu, gj)
    lam, V = np.linalg.eigh(K)
    rt = np.sqrt(lam)
    v0 = np.zeros(len(nu) + 1)
    v0[0] = 1.0
    cv = V.T @ v0
    ts = np.linspace(10, 80, 200)
    e_field = []
    for t in ts:
        x = V @ (np.sin(rt * t) / rt * cv)
        xd = V @ (np.cos(rt * t) * cv)
        e_field.append(0.5 * (xd[0] ** 2 + w0 ** 2 * x[0] ** 2))
    gamma_live = -float(np.polyfit(ts, np.log(e_field), 1)[0])

    # (c) the golden-rule formula gamma = pi*kappa/(2 w0^2) (the design)
    gamma_rule = math.pi * kappa / (2 * w0 ** 2)

    # (d) exact g^2 scaling on the continuum pole
    gs = np.array([0.05, 0.1, 0.2, 0.4])
    gams = np.array([-2 * resolvent_pole(g ** 2).imag for g in gs])
    slope = float(np.polyfit(np.log(gs), np.log(gams), 1)[0])
    rule_ratio = gams / (math.pi * gs ** 2 / (2 * w0 ** 2))

    ok = (abs(gamma_pole - gamma_set) / gamma_set < 0.01
          and abs(gamma_live - gamma_set) / gamma_set < 0.01
          and abs(gamma_rule - gamma_set) / gamma_set < 1e-12
          and abs(slope - 2.0) < 0.01
          and np.all(np.abs(rule_ratio - 1) < 0.01))
    out = {
        'name': 'T3_derived_damping',
        'description': (
            'the pole moves below the axis by a computed amount: '
            'complex pole, live decay, and golden rule agree; '
            'gamma scales exactly as g^2'
        ),
        'gamma_design': gamma_set,
        'gamma_from_pole': float(gamma_pole),
        'gamma_from_live_decay': gamma_live,
        'gamma_from_golden_rule': gamma_rule,
        'g2_scaling_slope': slope,
        'golden_rule_ratio_over_g_grid': [float(r) for r in rule_ratio],
        'pass': bool(ok),
    }
    _CACHE['T3'] = out
    return out


# ========================================================================
# T4. The geometry adjudicates placement
# ========================================================================


def _placement_couplings(placement: str, wn: np.ndarray) -> np.ndarray:
    if placement == 'antipode':
        # Y_n(pi) = lim sin(n psi)/sin(psi) / sqrt(2 pi^2)
        return wn * (-1.0) ** (wn + 1) / _SQ
    if placement == 'generic':
        psi0 = 1.0
        return np.sin(wn * psi0) / (math.sin(psi0) * _SQ)
    if placement == 'uniform':
        # exact selection rule: int sin(n psi) sin psi dpsi = (pi/2) d_{n1}
        c = np.zeros(len(wn))
        c[0] = 1.0 / _SQ
        return c
    raise ValueError(placement)


def test_T4_placement_adjudication() -> dict:
    if 'T4' in _CACHE:
        return _CACHE['T4']
    M, N, T, g = 24, 3000, 600.0, 0.3
    wn = np.arange(1, M + 1).astype(float)
    nu = np.linspace(0.5, 26.0, N)
    dnu = nu[1] - nu[0]
    gj = g * np.sqrt(dnu) * np.ones(N)

    # (i) the antipodal flat-rate law: gamma_n = pi (g Y_n)^2 rho/(2 w_n^2)
    # with Y_n = n/sqrt(2 pi^2) gives gamma = g^2/(4 pi) for EVERY n
    gam_flat = g ** 2 / (4 * math.pi)
    resid_pred = math.exp(-gam_flat * T)

    # (ii) the antipodal coupling is the psi -> pi limit of Y_n, exact
    psi = math.pi - 1e-6
    lim_err = float(np.max(np.abs(
        np.sin(wn * psi) / (math.sin(psi) * _SQ)
        - _placement_couplings('antipode', wn))))

    # (iii) the uniform selection rule by quadrature
    pgrid = np.linspace(0, math.pi, 20001)
    sel = [np.trapezoid(np.sin(n * pgrid) * np.sin(pgrid), pgrid)
           for n in (1, 2, 3, 7)]
    sel_ok = (abs(sel[0] - math.pi / 2) < 1e-8
              and all(abs(s) < 1e-8 for s in sel[1:]))

    # (iv) the generic-point near-zero
    min_sin = float(np.min(np.abs(np.sin(wn * 1.0))))

    # (v) the live three-way adjudication: kicked multi-mode field,
    # exact normal-mode evolution to t = T
    sig = 0.1
    v_field = wn * np.exp(-0.5 * (wn * sig) ** 2)
    v_field /= np.sqrt(np.sum(v_field ** 2))
    residuals, per_mode = {}, {}
    for placement in ('antipode', 'generic', 'uniform'):
        cn = _placement_couplings(placement, wn)
        K = np.zeros((M + N, M + N))
        K[np.arange(M), np.arange(M)] = wn ** 2
        K[np.arange(M, M + N), np.arange(M, M + N)] = nu ** 2
        K[:M, M:] = np.outer(cn, gj)
        K[M:, :M] = np.outer(gj, cn)
        x0 = np.zeros(M + N)
        v0 = np.zeros(M + N)
        v0[:M] = v_field
        x, xd, lam = normal_mode_state(K, x0, v0, T)
        if lam.min() <= 0:
            residuals[placement] = float('nan')
            continue
        em = 0.5 * (xd[:M] ** 2 + wn ** 2 * x[:M] ** 2)
        residuals[placement] = float(np.sum(em) / 0.5)
        per_mode[placement] = em / (0.5 * v_field ** 2)

    anti_range = (float(per_mode['antipode'].min()),
                  float(per_mode['antipode'].max()))
    ok = (abs(residuals['antipode'] - resid_pred) < 0.01
          and anti_range[1] < 2.0 * resid_pred
          and anti_range[0] > 0.3 * resid_pred
          and residuals['generic'] > 0.85
          and residuals['uniform'] > 0.98
          and lim_err < 1e-4 and sel_ok and min_sin < 0.01)
    out = {
        'name': 'T4_placement_adjudication',
        'description': (
            'antipodal point absorber: flat rate g^2/(4 pi) for every '
            'mode (impedance matched); generic point: near-zeros, '
            '1/n^2 rates; uniform distributed: couples only to n = 1'
        ),
        'flat_rate_g2_over_4pi': gam_flat,
        'predicted_antipodal_residual': resid_pred,
        'live_residuals': residuals,
        'antipodal_per_mode_residual_range': anti_range,
        'antipodal_limit_error': lim_err,
        'uniform_selection_integrals': [float(s) for s in sel],
        'generic_min_coupling': min_sin,
        'pass': bool(ok),
    }
    _CACHE['T4'] = out
    return out


# ========================================================================
# T5. Recurrence and the order of limits
# ========================================================================


def test_T5_recurrence() -> dict:
    if 'T5' in _CACHE:
        return _CACHE['T5']
    w0 = _W0
    gamma_set = 0.05
    kappa = 2 * gamma_set * w0 ** 2 / math.pi
    revivals = {}
    for n_bath in (300, 600):
        nu, dnu = flat_bank(n_bath)
        gj = np.sqrt(kappa * dnu) * np.ones(n_bath)
        K = stiffness_matrix(w0, nu, gj)
        lam, V = np.linalg.eigh(K)
        rt = np.sqrt(lam)
        v0 = np.zeros(n_bath + 1)
        v0[0] = 1.0
        cv = V.T @ v0
        t_pred = 2 * math.pi / dnu
        tt = np.linspace(0.5 * t_pred, 1.5 * t_pred, 600)
        e = []
        for t in tt:
            x = V @ (np.sin(rt * t) / rt * cv)
            xd = V @ (np.cos(rt * t) * cv)
            e.append(0.5 * (xd[0] ** 2 + w0 ** 2 * x[0] ** 2))
        e = np.array(e)
        revivals[n_bath] = {
            't_rec_predicted': float(t_pred),
            't_revival_measured': float(tt[np.argmax(e)]),
            'revival_energy': float(np.max(e)),
        }
    r3, r6 = revivals[300], revivals[600]
    ratio = r6['t_revival_measured'] / r3['t_revival_measured']
    ok = (abs(r3['t_revival_measured'] / r3['t_rec_predicted'] - 1) < 0.15
          and abs(r6['t_revival_measured'] / r6['t_rec_predicted'] - 1) < 0.15
          and abs(ratio - 2.0) < 0.15
          and r3['revival_energy'] > 0.1)
    out = {
        'name': 'T5_recurrence',
        'description': (
            'finite bank revives near T_rec = 2 pi/dnu, linear in N: '
            'epsilon finite = Poincare recurrence on the closed bulk; '
            'eps -> 0 is the continuum limit N -> infinity FIRST'
        ),
        'revivals': revivals,
        'revival_time_scaling_ratio': float(ratio),
        'pass': bool(ok),
    }
    _CACHE['T5'] = out
    return out


# ========================================================================
# T6. The #213 closure: epsilon filled in by the absorber
# ========================================================================


def test_T6_213_closure() -> dict:
    w0 = _W0

    # (i) eps_derived = w_tilde * gamma: the Delta^2 pole form
    # 1/(Delta^2 - w_t^2 + i eps_d) has its pole at sqrt(w_t^2 - i eps_d),
    # agreeing with the true pole w_t - i gamma/2 to O(gamma^2)
    pole_form_gap, gamma_sq = [], []
    for g in (0.2, 0.4):
        om = resolvent_pole(g ** 2)
        wt, gm = om.real, -2 * om.imag
        eps_d = wt * gm
        pole2 = np.sqrt(wt ** 2 - 1j * eps_d)
        pole_form_gap.append(abs(pole2 - (wt - 1j * gm / 2)))
        gamma_sq.append(gm ** 2)
    o2_ok = all(gap < 0.1 * g2 for gap, g2 in zip(pole_form_gap, gamma_sq))

    # (ii) the ordering integral with the derived epsilon: the damped
    # offer segment (i/2 w_t) e^{-i w_t t - gamma t/2} half-line
    # transforms to -(1/2 w_t)/(Delta - w_t + i gamma/2) - the i*epsilon
    # of #213's orderings, now a physical rate (numerical quadrature)
    om = resolvent_pole(0.2 ** 2)
    wt, gm = om.real, -2 * om.imag
    delta = 2.0
    t = np.linspace(0, 40 / gm, 4000001)
    f = (1j / (2 * wt)) * np.exp(1j * delta * t - 1j * wt * t - gm * t / 2)
    i_quad = np.trapezoid(f, t)
    i_closed = -(1 / (2 * wt)) / (delta - wt + 1j * gm / 2)
    ordering_err = abs(i_quad - i_closed)

    # (iii) g -> 0 recovers #213: the damped kernel converges to
    # (i/2 w0) e^{-i w0 |t|} pointwise as gamma -> 0
    tg = np.linspace(-20, 20, 2001)
    recover = []
    for g in (0.2, 0.1, 0.05):
        omg = resolvent_pole(g ** 2)
        wtg, gmg = omg.real, -2 * omg.imag
        kern = (1j / (2 * wtg)) * np.exp(-1j * wtg * np.abs(tg)
                                         - gmg * np.abs(tg) / 2)
        kern0 = (1j / (2 * w0)) * np.exp(-1j * w0 * np.abs(tg))
        recover.append(float(np.max(np.abs(kern - kern0))))
    monotone = recover[0] > recover[1] > recover[2]

    ok = o2_ok and ordering_err < 1e-4 and monotone and recover[-1] < 0.02
    return {
        'name': 'T6_213_closure',
        'description': (
            'eps_derived = w_tilde*gamma fills in the orderings\' '
            'i*epsilon as a physical rate; the pole form agrees to '
            'O(gamma^2); g -> 0 recovers the #213 kernel'
        ),
        'pole_form_gap_vs_gamma_sq': [
            [float(a), float(b)] for a, b in zip(pole_form_gap, gamma_sq)],
        'ordering_integral_error': float(ordering_err),
        'kernel_recovery_sequence': recover,
        'pass': bool(ok),
    }


# ========================================================================
# T7. Honest scope
# ========================================================================


def test_T7_honest_scope() -> dict:
    scope = [
        'The absorber is Gaussian/bilinear BY DESIGN - that is what '
        'makes the elimination exact; anharmonic absorber dynamics '
        '(real registration/irreversibility) is the standing #209 '
        'open, not claimed here.',
        'The flat bank density is a modeling choice; the physical '
        'absorber is the throat\'s internal continuum - the engine\'s '
        'MouthState.modes are exactly this degree of freedom (this PR '
        'gives it an action; deriving its spectral density from the '
        'throat geometry is the successor).',
        'The source point psi = 0 shares the antipode\'s matching '
        'property (Y_n(0) = n/sqrt(2 pi^2)): source and antipode are '
        'the two distinguished points; the antipode is the unique '
        'OTHER one - recorded, not hidden.',
        'Classical throughout; per-mode field; frozen geometry, no '
        'backreaction.',
        'The adjudication is at fixed absorber budget (one point vs '
        'one distribution); a dense random cloud of point absorbers '
        'also works generically - the sharp statements are the flat '
        'antipodal rate law and the exact uniform selection rule.',
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
    t2 = test_T2_exact_elimination()
    t3 = test_T3_derived_damping()
    t4 = test_T4_placement_adjudication()
    t5 = test_T5_recurrence()
    t6 = test_T6_213_closure()
    core = all(t['pass'] for t in (t2, t3, t4, t5, t6))
    assessment = (
        'The absorber is now a degree of freedom, and everything #213 '
        'imposed is an outcome: absorption happens dynamically (live '
        'energy decay matching the derived rate), the i*epsilon is the '
        'computed absorber response rate (pole, live decay, and golden '
        'rule agreeing, scaling exactly as g^2, vanishing with the '
        'coupling to recover #213), the finite-bank recurrences locate '
        'the physical meaning of epsilon on a closed bulk, and the '
        'geometry itself picks the absorber\'s address: the antipodal '
        'point is impedance-matched to every tower mode - the #166 '
        'caustic in its absorber role - while a uniformly distributed '
        'absorber is forbidden by an exact selection rule from '
        'absorbing anything above the ground mode.'
        if core else
        'INCONCLUSIVE - a core check failed; do not quote.'
    )
    return {
        'name': 'T8_assessment',
        'description': 'the standing of the dynamical-absorber derivation',
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
        test_T2_exact_elimination(),
        test_T3_derived_damping(),
        test_T4_placement_adjudication(),
        test_T5_recurrence(),
        test_T6_213_closure(),
        test_T7_honest_scope(),
        test_T8_assessment(),
    ]
    t3, t4, t5, t6 = tests[2], tests[3], tests[4], tests[5]
    all_ok = all(t['pass'] for t in tests)
    if all_ok:
        verdict_class = (
            "EPSILON_IS_THE_ABSORBER_RESPONSE_RATE_AND_THE_ANTIPODE_"
            "IS_THE_IMPEDANCE_MATCHED_ABSORBER_OF_THE_CLOSED_BULK"
        )
        r = t4['live_residuals']
        verdict = (
            "DERIVED (the argument is in "
            "docs/dynamical_absorber_propagator.md).\n\n"
            "THE PROMOTION. S = S_field + S_absorber + S_coupling, all "
            "Gaussian: the elimination is exact (resolvent identity to "
            "machine precision), and #213's imposed epsilon is now the "
            "computed damping rate of the absorber response - the "
            "complex pole, the live time-domain decay, and the golden "
            f"rule agree ({t3['gamma_from_pole']:.4f} / "
            f"{t3['gamma_from_live_decay']:.4f} / "
            f"{t3['gamma_from_golden_rule']:.4f}), scaling exactly as "
            f"g^2 (slope {t3['g2_scaling_slope']:.4f}), with g -> 0 "
            "recovering the #213 kernel.\n\n"
            "THE ADDRESS. The geometry adjudicates placement: the "
            "antipodal point absorber couples as Y_n ~ n - exactly "
            "cancelling the kinematic suppression - so EVERY tower "
            f"mode damps at the same rate g^2/(4 pi); live residuals "
            f"{r['antipode']:.3f} (antipode, uniform across modes) vs "
            f"{r['generic']:.2f} (generic point, near-zero couplings) "
            f"vs {r['uniform']:.2f} (uniform distribution - an exact "
            "selection rule confines it to the ground mode). The #166 "
            "focusing caustic is the impedance matching.\n\n"
            "THE MEANING OF EPSILON. A finite bank revives at "
            "T_rec = 2 pi/dnu (linear in N, measured ratio "
            f"{t5['revival_time_scaling_ratio']:.2f}): on the closed "
            "bulk, finite epsilon is Poincare recurrence, and the "
            "eps -> 0 limit is the continuum limit of the throat's "
            "internal spectrum taken first. The orderings' i*epsilon "
            "is filled in as a physical rate (quadrature error "
            f"{t6['ordering_integral_error']:.0e}), and the covariant "
            "pole form agrees to O(gamma^2)."
        )
    else:
        verdict_class = "DYNAMICAL_ABSORBER_INCONCLUSIVE"
        verdict = (
            "INCONCLUSIVE. A core check failed; re-examine before "
            "quoting the derivation."
        )

    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "identification": (
            "The dynamical absorber: S = S_field + S_absorber + "
            "S_coupling on the frozen bulk - absorption as an outcome, "
            "epsilon as the computed absorber response rate, and the "
            "antipodal point as the geometry's own impedance-matched "
            "absorber (vs generic-point remnants and the exact "
            "uniform-distribution selection rule)"
        ),
        "executes": (
            "the decisive #213 successor: the absorber boundary "
            "condition promoted to a degree of freedom"
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
    return str(o)


def render_markdown(s: dict) -> str:
    out = []
    out.append("# The dynamical absorber (PR #214)")
    out.append("")
    out.append(f"**Run:** {s['timestamp_utc']}")
    out.append("")
    out.append(
        "The deliverable is `docs/dynamical_absorber_propagator.md` - "
        "the absorber promoted from a boundary condition to a degree "
        "of freedom. *(QFT on the fixed classical throat geometry, "
        "not quantum gravity.)*"
    )
    out.append("")
    out.append("## Test summary")
    out.append("")
    out.append("| # | Test | Key finding | PASS? |")
    out.append("|---|---|---|---|")
    labels = {
        "T1": "the goal: absorption as an outcome",
        "T2": "exact elimination of the Gaussian absorber",
        "T3": "the damping derived: pole = live decay = rule",
        "T4": "the antipode is the impedance-matched absorber",
        "T5": "finite bank = recurrence; the order of limits",
        "T6": "the #213 epsilon filled in as a physical rate",
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
    out = here / "runs" / f"{ts}_dynamical_absorber_propagator_probe"
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
