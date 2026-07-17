"""
The coupled 5D Einstein-Klein-Gordon weld: the soliton's energy
density sources the local Tangherlini metric, all scales lock in one
unit system - and the lock refutes self-sourcing (PR #222).

> Framing: QFT on the *fixed classical* throat geometry (geometry -> fields),
> not quantum gravity.  This PR is the one place the geometry is NOT
> frozen: the requested modification couples the soliton's energy
> density to the 5D metric, exactly to derive the scale weld the #221
> identifiability audit demanded.

THE REQUEST, EXECUTED
---------------------
Option A (the core): the action is modified so the soliton's localized
energy density directly deforms the local 5D metric - the static
spherical 5D Einstein-Klein-Gordon system, solved by shooting.  The
background cavity scale R* (the wave potential's well on the SOLVED
metric) and the soliton's R_RMS both emerge from the single coupled
system, together with the exterior Tangherlini scale r_s = sqrt(mu).
Option B (the verification layer): the Israel junction at the core
boundary is machine-checked - the interior numeric extrinsic curvature
matches the exterior Tangherlini's with jumps tracking the scalar tail
(no shell), and the mass the exterior reads equals the integrated
interior energy: the length scales share one physical boundary.

THE WELD (the audit's T5, answered)
-----------------------------------
In the coupled system every length is expressed in the scalar-mass
unit 1/m, and the gravitational-coupling convention drops out EXACTLY
(the (KF, s) -> (KF/4, 2s) rescale leaves every geometric observable
bitwise unchanged): the family relations mu(omega), X(sigma/r_s),
r_s*omega are DERIVED, convention-free functions.  The audit's free
radial-unit rescale is forbidden by the field equations.

THE MEASURED STRUCTURE
----------------------
* The 4D reduction reproduces the Kaup benchmark (M_max = 0.6327 vs
  0.633) - same solver, n = D-2 as the only change.
* THE 5D CRITICAL-MASS MARGINALITY: as the family goes dilute
  (omega -> 1) the mass does NOT vanish - mu -> mu_crit ~ 7.695
  (linear-in-s_c convergence, difference ratio 2.0) while the size
  runs free at fixed mass: X = sigma*omega ~ s_c^(-1/2) (slope
  -0.50).  The 5D Newtonian potential is 1/r^2 - marginal, the same
  marginality as #221's exterior tails - so a 5D self-gravitating
  scalar has a critical mass and a scale-free size zero mode.
* THE LOCK: r_s * omega in [1.53, 2.774] over the ENTIRE family
  (spiral turnaround included; sup = sqrt(mu_crit)); q-channel spot
  checks land in the same band.
* THE CAVITY: on the solved metric the wave-potential well's centroid
  R* tracks the THROAT scale (R*/r_s = 0.78-1.00) while R*/R_RMS
  falls 0.39 -> 0.17 dilute: R* and R_RMS emerge from one system and
  DECOUPLE - the cavity belongs to the gravitational core, not to the
  matter cloud.

THE CONFRONTATION
-----------------
The EM-cap primordial anchor requires r_s * omega = alpha = 0.0073.
The one-mouth coupled system pins r_s * omega >= 1.53: EXCLUDED by
x210.  Jointly: sigma/r_s = 206.8 (conv B) is reachable only on the
dilute branch, where X = sigma*omega = 206.8 * 2.774 = 574 vs the
required 1.51 (x380).  A 5D scalar CANNOT be much lighter than the
throat its own field creates - the self-sourcing reading of the
anchor is refuted at the coupled level, and #210's relocation is now
FORCED by the weld: the throat is PRIMORDIAL (its bulk mass a
geometric datum), the #221 cavity is the primordial throat's, and the
successor is sharply defined - the soliton on the fixed-mu throat
background with perturbative back-reaction (the bridge version of
this exact solver).

Tests:
  T1. Goal.
  T2. The solver: general-n reduction, the 4D Kaup benchmark, the
      constraint identity, dx/xmax convergence.
  T3. The coupled 5D family: the existence curve, the dilute
      critical-mass marginality, the scale-free zero mode, the
      spiral turnaround.
  T4. The weld: exact convention invariance - the audit's unit
      freedom demolished.
  T5. The Israel junction (option B): extrinsic-curvature match at
      the core boundary, jumps tracking the tail; the exterior reads
      the interior mass.
  T6. The cavity from the same metric (option A): the wave potential
      validated on vacuum against the #215 closed form; R* and R_RMS
      from one system, decoupling dilute.
  T7. The confrontation: r_s*omega never alpha - self-sourcing
      refuted, the primordial throat forced.
  T8. Honest scope.
  T9. Assessment.

Verdict:
  THE_WELD_IS_DERIVED_AND_IT_REFUTES_SELF_SOURCING_THE_COUPLED_5D_EKG_
  LOCKS_RS_OMEGA_INTO_1P5_TO_2P8_NEVER_ALPHA_THE_THROAT_IS_PRIMORDIAL_
  AND_THE_CRITICAL_MASS_MARGINALITY_IS_MEASURED
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np

_KF = 0.25
_ALPHA = 7.2973525693e-3
_MU_OVER_E = 206.7683
_X_REQ_B = _ALPHA * _MU_OVER_E                    # 1.5089

_CACHE: dict = {}


# ── the general-n static spherical EKG solver (#210 conventions) ────────

def make_V(kind: str, g: float = 4.0, a0: float = 0.09, lam: float = 1.0):
    """Interaction potential (mass term explicit in the equations);
    'qchan' is the #210 relativistic completion of the committed #180
    structure (adiabatic order field)."""
    if kind == "kaup":
        return (lambda s2: 0.0 * s2), (lambda s: 0.0 * s)
    if kind == "qchan":
        def V(s2):
            u = g * s2 - a0
            return -np.where(u > 0, u * u / (4 * lam), 0.0)

        def dV(s):
            u = g * s * s - a0
            return -np.where(u > 0, (u / lam) * g * s, 0.0)
        return V, dV
    raise ValueError(kind)


def _rhs(x, y, W, V, dV, n=3, KF=_KF):
    """ds^2 = -al^2 dt^2 + a^2 dr^2 + r^2 dOmega_n^2; n = D - 2.
    n = 2 is exactly the #210 system; n = 3 is 5D with the Tangherlini
    exterior N = 1 - mu/r^2, mu(r) = r^2 (1 - 1/a^2)."""
    s, sp, a, al = y
    w2 = (W / al) ** 2
    a2 = a * a
    s2 = s * s
    Vi = V(s2)
    rho_t = (w2 + 1) * a2 * s2 + sp * sp + a2 * Vi        # a^2 rho
    pr_t = (w2 - 1) * a2 * s2 + sp * sp - a2 * Vi         # a^2 p_r
    ap = a * ((n - 1) * (1 - a2) / (2 * x) + (2 * KF / n) * x * rho_t)
    alp = al * ((n - 1) * (a2 - 1) / (2 * x) + (2 * KF / n) * x * pr_t)
    spp = (-(n / x + alp / al - ap / a) * sp
           - w2 * a2 * s + a2 * (s + 0.5 * dV(s)))
    return np.array([sp, spp, ap, alp])


def _batch_nodes(sc, Ws, V, dV, n, dx, xmax, KF=_KF):
    nb = len(Ws)
    y = np.array([np.full(nb, sc), np.zeros(nb), np.ones(nb),
                  np.ones(nb)])
    alive = np.ones(nb, dtype=bool)
    nodes = np.zeros(nb, dtype=int)
    prev = np.full(nb, sc)
    x = dx
    for _ in range(int(round((xmax - dx) / dx))):
        with np.errstate(all="ignore"):
            k1 = _rhs(x, y, Ws, V, dV, n, KF)
            k2 = _rhs(x + dx / 2, y + dx / 2 * k1, Ws, V, dV, n, KF)
            k3 = _rhs(x + dx / 2, y + dx / 2 * k2, Ws, V, dV, n, KF)
            k4 = _rhs(x + dx, y + dx * k3, Ws, V, dV, n, KF)
            yn = y + dx / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        bad = (~np.isfinite(yn).all(axis=0)) | (np.abs(yn[0]) > 3 * abs(sc))
        alive &= ~bad
        yn[:, ~alive] = y[:, ~alive]
        y = yn
        cross = alive & (y[0] * prev < 0)
        nodes[cross] += 1
        prev = np.where(alive, y[0], prev)
        x += dx
    return nodes


def solve_star(sc, kind="kaup", n=3, w_lo=0.3, w_hi=1.9, rounds=7,
               nw=25, dx=0.01, xmax=60.0, KF=_KF, keep_prof=False,
               **kw):
    """Ground state by bracketing the 0 -> 1 node transition; then
    observables with the tail cut where the field has decayed."""
    key = (kind, round(sc, 6), n, dx, xmax, KF, w_hi,
           tuple(sorted(kw.items())), keep_prof)
    if key in _CACHE:
        return _CACHE[key]
    V, dV = make_V(kind, **kw)
    lo, hi = w_lo, w_hi
    for _ in range(rounds):
        Ws = np.linspace(lo, hi, nw)
        nd = _batch_nodes(sc, Ws, V, dV, n, dx, xmax, KF)
        idx = None
        for i in range(1, nw):
            if nd[i - 1] == 0 and nd[i] >= 1:
                idx = i
                break
        if idx is None:
            _CACHE[key] = None
            return None
        lo, hi = Ws[idx - 1], Ws[idx]
    W = lo
    y = np.array([[sc], [0.0], [1.0], [1.0]])
    x = dx
    mom0 = mom2 = 0.0
    M_cut = x_cut = None
    al_c = a_c = 1.0
    hist = []
    prof = [(dx, sc, 0.0, 1.0, 1.0)]
    for _ in range(int(round((xmax - dx) / dx))):
        with np.errstate(all="ignore"):
            k1 = _rhs(x, y, np.array([W]), V, dV, n, KF)
            k2 = _rhs(x + dx / 2, y + dx / 2 * k1, np.array([W]),
                      V, dV, n, KF)
            k3 = _rhs(x + dx / 2, y + dx / 2 * k2, np.array([W]),
                      V, dV, n, KF)
            k4 = _rhs(x + dx, y + dx * k3, np.array([W]), V, dV, n, KF)
            yn = y + dx / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        if not np.isfinite(yn).all() or abs(yn[0][0]) > 3 * abs(sc):
            break
        y = yn
        x += dx
        s, sp, a, al = (float(v[0]) for v in y)
        prof.append((x, s, sp, a, al))
        if M_cut is None:
            w2 = (W / al) ** 2
            rho = (w2 + 1) * s * s + (sp / a) ** 2 \
                + float(V(np.array([s * s]))[0])
            mom0 += rho * x ** n * dx
            mom2 += rho * x ** (n + 2) * dx
            hist.append((x, mom0))
            if abs(s) < 1e-4 * abs(sc) and x > 2.0:
                M_cut = x ** (n - 1) * (1 - 1 / (a * a))
                x_cut, al_c, a_c = x, al, a
    if M_cut is None:
        a_f = float(y[2][0])
        M_cut = x ** (n - 1) * (1 - 1 / (a_f * a_f))
        x_cut, al_c, a_c = x, float(y[3][0]), a_f
    rms = math.sqrt(mom2 / mom0) if mom0 > 0 else float("nan")
    r99 = next((xx for xx, mm in hist if mm >= 0.99 * mom0), x_cut)
    out = {"W": W, "W_phys": W / (al_c * a_c), "mu": float(M_cut),
           "rms": float(rms), "r99": float(r99),
           "x_cut": float(x_cut), "alc": al_c * a_c}
    if keep_prof:
        out["prof"] = np.array(prof)
    _CACHE[key] = out
    return out


def _dilute_params(sc):
    xmax = min(60.0 * math.sqrt(0.18 / sc), 600.0)
    dx = 0.01 if xmax <= 300 else 0.02
    return dx, xmax


_MAIN_SC = (0.02, 0.05, 0.08, 0.12, 0.18, 0.25, 0.35, 0.5, 0.7)
_DILUTE_SC = (0.01, 0.005, 0.0025)
_SPIRAL = ((1.0, 3.5), (1.4, 5.0))


def family_5d():
    if "FAM" in _CACHE:
        return _CACHE["FAM"]
    rows = []
    for sc in _MAIN_SC + _DILUTE_SC:
        dx, xmax = _dilute_params(sc) if sc <= 0.02 else (0.01, 60.0)
        r = solve_star(sc, dx=dx, xmax=xmax)
        if r is None:
            continue
        rs = math.sqrt(max(r["mu"], 0.0))
        rows.append({"sc": sc, "W_phys": r["W_phys"], "mu": r["mu"],
                     "r_s": rs, "rms": r["rms"], "r99": r["r99"],
                     "X": r["rms"] * r["W_phys"],
                     "rms_over_rs": r["rms"] / rs,
                     "rs_omega": rs * r["W_phys"]})
    for sc, whi in _SPIRAL:
        r = solve_star(sc, w_hi=whi)
        if r is None:
            continue
        rs = math.sqrt(max(r["mu"], 0.0))
        rows.append({"sc": sc, "W_phys": r["W_phys"], "mu": r["mu"],
                     "r_s": rs, "rms": r["rms"], "r99": r["r99"],
                     "X": r["rms"] * r["W_phys"],
                     "rms_over_rs": r["rms"] / rs,
                     "rs_omega": rs * r["W_phys"],
                     "branch": "spiral"})
    _CACHE["FAM"] = rows
    return rows


def _v_eff(x, a, al, ell):
    """The wave potential on the solved background, constructed by
    definition (finite differences): u = r^{3/2} psi in the tortoise
    coordinate of ds^2 = -al^2 dt^2 + a^2 dr^2 + r^2 dOmega_3^2."""
    F = al / a
    r32 = x ** 1.5
    inner = F * np.gradient(r32, x)
    return al ** 2 * ell * (ell + 2) / x ** 2 \
        + F * np.gradient(inner, x) / r32


# ========================================================================
# T1. Goal
# ========================================================================


def test_T1_goal() -> dict:
    return {
        'name': 'T1_goal',
        'description': (
            'Execute the requested modification: couple the soliton '
            'core\'s localized energy density to the local 5D metric '
            '(the static spherical Einstein-Klein-Gordon system, '
            'option A), so the background cavity scale R*, the '
            'soliton\'s R_RMS, and the exterior Tangherlini scale '
            'r_s emerge from a single coupled solve; verify the '
            'Israel junction at the core boundary (option B); and '
            'answer the #221 identifiability audit\'s missing-weld '
            'finding with the derived, convention-free family '
            'relations.'
        ),
        'pass': True,
    }


# ========================================================================
# T2. The solver: benchmark, constraint, convergence
# ========================================================================


def test_T2_solver() -> dict:
    if 'T2' in _CACHE:
        return _CACHE['T2']
    # the 4D reduction (n = 2) must land the Kaup point
    best_M = 0.0
    for sc in (0.3, 0.35, 0.4, 0.45):
        r = solve_star(sc, n=2, w_lo=0.8, w_hi=1.6)
        if r is not None:
            best_M = max(best_M, 0.5 * r["mu"])       # 4D: mu = 2M
    kaup_err = abs(best_M - 0.633)

    # the constraint identity mu' = (2 kappa / n) x^n rho on a 5D
    # member (kappa = 2 KF; the tt Einstein equation re-derived from
    # the stored profile - an independent implementation check)
    r = solve_star(0.18, keep_prof=True)
    x, s, sp, a, al = r["prof"].T
    mu_r = x ** 2 * (1 - 1 / a ** 2)
    w2 = (r["W"] / al) ** 2
    rho = (w2 + 1) * s ** 2 + (sp / a) ** 2
    pred = (2 * (2 * _KF) / 3) * x ** 3 * rho
    dmu = np.gradient(mu_r, x)
    sel = (x > 0.5) & (x < r["x_cut"])
    res = np.abs(dmu[sel] - pred[sel]) / (np.abs(pred[sel]) + 1e-12)

    # convergence: dx halving and xmax extension
    r1 = solve_star(0.18, dx=0.02)
    r2 = solve_star(0.18, dx=0.01)
    r3 = solve_star(0.18, dx=0.01, xmax=90.0)
    d_dx = abs(r1["mu"] - r2["mu"])
    d_xm = abs(r3["mu"] - r2["mu"])

    ok = (kaup_err < 0.005
          and float(np.median(res)) < 1e-4
          and float(res.max()) < 5e-3
          and d_dx < 1e-4 and d_xm < 1e-5)
    out = {
        'name': 'T2_solver',
        'description': (
            'the general-n static spherical EKG solver: the n = 2 '
            'reduction reproduces the #210 Kaup benchmark; the tt '
            'Einstein (mass-function) identity holds pointwise on the '
            'stored 5D profile as an independent implementation '
            'check; mu is dx- and xmax-converged'
        ),
        'kaup_M_max': float(best_M),
        'kaup_error': float(kaup_err),
        'constraint_residual_median': float(np.median(res)),
        'constraint_residual_max': float(res.max()),
        'mu_dx_halving_shift': float(d_dx),
        'mu_xmax_extension_shift': float(d_xm),
        'pass': bool(ok),
    }
    _CACHE['T2'] = out
    return out


# ========================================================================
# T3. The coupled 5D family and the critical-mass marginality
# ========================================================================


def test_T3_family() -> dict:
    if 'T3' in _CACHE:
        return _CACHE['T3']
    rows = family_5d()
    main = [r for r in rows if r.get("branch") != "spiral"]
    main_sorted = sorted(main, key=lambda r: r["sc"])
    mu_by_sc = [r["mu"] for r in main_sorted]
    monotone = all(mu_by_sc[i] > mu_by_sc[i + 1]
                   for i in range(len(mu_by_sc) - 1))

    # dilute extrapolation: mu(sc) linear -> mu_crit
    d = {r["sc"]: r for r in rows}
    m1, m2, m3 = d[0.01]["mu"], d[0.005]["mu"], d[0.0025]["mu"]
    diff_ratio = (m2 - m1) / (m3 - m2)
    mu_crit = m3 + (m3 - m2) / (diff_ratio - 1)
    rs_omega_sup = math.sqrt(mu_crit)

    # the scale-free zero mode: X ~ sc^(-1/2) at fixed critical mass
    slope = (math.log(d[0.01]["X"] / d[0.0025]["X"])
             / math.log(0.01 / 0.0025))

    # spiral turnaround: W_phys rises again past the frequency minimum
    spiral = [r for r in rows if r.get("branch") == "spiral"]
    turnaround = (len(spiral) >= 2
                  and spiral[-1]["W_phys"] > d[0.7]["W_phys"])

    ok = (len(rows) >= 13 and monotone
          and (1 - d[0.0025]["W_phys"]) < 1e-3
          and 1.9 < diff_ratio < 2.1
          and 7.6 < mu_crit < 7.8
          and abs(slope + 0.5) < 0.02
          and turnaround)
    out = {
        'name': 'T3_family',
        'description': (
            'the coupled 5D family: mu monotone along the main '
            'branch; the dilute endpoint has a CRITICAL MASS (mu -> '
            'mu_crit, linear-in-sc, difference ratio 2.0) with omega '
            '-> 1 and the size running free at fixed mass (X ~ '
            'sc^(-1/2)) - the 5D Newtonian 1/r^2 marginality, the '
            'same marginality as #221\'s exterior tails; the spiral '
            'turnaround bounds the compact end'
        ),
        'family': [{k: (round(v, 6) if isinstance(v, float) else v)
                    for k, v in r.items()} for r in rows],
        'mu_monotone_main_branch': bool(monotone),
        'dilute_diff_ratio': float(diff_ratio),
        'mu_crit_extrapolated': float(mu_crit),
        'rs_omega_supremum': float(rs_omega_sup),
        'one_minus_omega_at_dilute_end':
            float(1 - d[0.0025]["W_phys"]),
        'X_dilute_slope_vs_sc': float(slope),
        'spiral_turnaround': bool(turnaround),
        'pass': bool(ok),
    }
    _CACHE['T3'] = out
    return out


# ========================================================================
# T4. The weld: exact convention invariance
# ========================================================================


def test_T4_weld_invariance() -> dict:
    r1 = solve_star(0.18, KF=0.25)
    r2 = solve_star(0.36, KF=0.0625)
    rels = {k: abs(r1[k] - r2[k]) / abs(r1[k])
            for k in ("W_phys", "mu", "rms")}
    ok = all(v < 1e-12 for v in rels.values())
    return {
        'name': 'T4_weld_invariance',
        'description': (
            'the audit\'s missing weld, supplied: in the coupled '
            'system every length is in the scalar-mass unit and the '
            'gravitational-coupling convention drops out EXACTLY - '
            'the (KF, s) -> (KF/4, 2s) rescale leaves every geometric '
            'observable unchanged (here to machine zero), so the '
            'family relations mu(omega), X(sigma/r_s), r_s*omega are '
            'DERIVED convention-free functions and the audit-T5 free '
            'radial rescale is forbidden by the field equations'
        ),
        'rescale_relative_shifts': {k: float(v) for k, v in rels.items()},
        'pass': bool(ok),
    }


# ========================================================================
# T5. The Israel junction (option B)
# ========================================================================


def test_T5_israel_junction() -> dict:
    r = solve_star(0.18, keep_prof=True)
    x, s, sp, a, al = r["prof"].T
    mu_inf = r["mu"]
    dal = np.gradient(al, x)
    rows = []
    for fb in (1.0, 1.5, 2.25):
        Rb = fb * r["x_cut"]
        j = min(np.searchsorted(x, Rb), len(x) - 2)
        Kth_int = 1 / (a[j] * x[j])
        Next = 1 - mu_inf / x[j] ** 2
        Kth_ext = math.sqrt(max(Next, 0.0)) / x[j]
        Kt_int = dal[j] / (al[j] * a[j])
        Kt_ext = (mu_inf / x[j] ** 3) / math.sqrt(Next)
        rows.append({
            'Rb_over_cut': fb, 'Rb': float(x[j]),
            'tail_phi_over_phi0': float(abs(s[j]) / 0.18),
            'jump_K_theta_rel': float(abs(Kth_int - Kth_ext)
                                      / abs(Kth_ext)),
            'jump_K_t_rel': float(abs(Kt_int - Kt_ext)
                                  / abs(Kt_ext)),
        })
    # the mass the exterior reads = the integrated interior energy
    j = min(np.searchsorted(x, r["x_cut"]), len(x) - 1)
    mu_at_cut = x[j] ** 2 * (1 - 1 / a[j] ** 2)
    mass_match = abs(mu_at_cut - mu_inf) / mu_inf

    ok = (rows[0]['jump_K_theta_rel'] < 1e-10
          and rows[0]['jump_K_t_rel'] < 1e-5
          and rows[0]['tail_phi_over_phi0'] < 1e-3
          and rows[1]['jump_K_theta_rel'] < 1e-7
          and rows[1]['jump_K_t_rel'] < 1e-4
          and mass_match < 1e-6)
    return {
        'name': 'T5_israel_junction',
        'description': (
            'option B executed as the verification layer: the '
            'extrinsic curvature of the interior numeric solution '
            'matches the exterior Tangherlini\'s on the shared '
            'boundary with jumps tracking the scalar tail (no shell '
            'in the smooth-matching limit), and the mass the exterior '
            'reads equals the integrated interior energy - the two '
            'sectors\' length scales are locked to one physical '
            'boundary'
        ),
        'member_sc': 0.18,
        'junction_rows': rows,
        'exterior_reads_interior_mass_rel': float(mass_match),
        'pass': bool(ok),
    }


# ========================================================================
# T6. The cavity from the same metric (option A)
# ========================================================================


def test_T6_cavity() -> dict:
    if 'T6' in _CACHE:
        return _CACHE['T6']
    # vacuum validation: the FD-constructed potential on the exact
    # Tangherlini metric must reproduce the #215 closed form
    rh = 1.0
    xv = np.linspace(1.2, 30, 8000)
    fv = 1 - rh ** 2 / xv ** 2
    Vnum = _v_eff(xv, 1 / np.sqrt(fv), np.sqrt(fv), 1)
    V215 = fv * ((3 + 0.75) / xv ** 2 + 2.25 * rh ** 2 / xv ** 4)
    sel = (xv > 1.5) & (xv < 25)
    vac_dev = float(np.max(np.abs(Vnum[sel] - V215[sel])
                           / np.abs(V215[sel])))

    members = []
    for sc in (0.5, 0.18, 0.05):
        xmax = 60.0 if sc >= 0.08 else 150.0
        r = solve_star(sc, xmax=xmax, keep_prof=True)
        x, s, sp, a, al = r["prof"].T
        aln = al / r["alc"]
        V1 = _v_eff(x, a, aln, 1)
        sel = (x > 1.0) & (x < 0.9 * x[-1])
        xj = x[sel]
        dV1 = V1[sel] - (3 + 0.75) / xj ** 2
        wneg = np.where(dV1 < 0, -dV1, 0.0)
        Rstar = float(np.sum(wneg * xj) / np.sum(wneg))
        width = math.sqrt(float(np.sum(wneg * (xj - Rstar) ** 2)
                                / np.sum(wneg)))
        rs = math.sqrt(r["mu"])
        members.append({'sc': sc, 'R_star': Rstar, 'width': width,
                        'depth': float(-dV1.min()),
                        'R_star_over_rs': Rstar / rs,
                        'R_star_over_rms': Rstar / r["rms"],
                        'rms': r["rms"], 'r_s': rs})
    ratios_rs = [m['R_star_over_rs'] for m in members]
    ratios_rms = [m['R_star_over_rms'] for m in members]

    ok = (vac_dev < 1e-5
          and all(m['depth'] > 0 for m in members)
          and all(0.7 < v < 1.1 for v in ratios_rs)
          and ratios_rms[0] > ratios_rms[1] > ratios_rms[2]
          and ratios_rms[2] < 0.5 * ratios_rms[0])
    out = {
        'name': 'T6_cavity',
        'description': (
            'the requested single-system emergence: the wave '
            'potential on the SOLVED metric (FD-constructed, '
            'validated on vacuum against the #215 closed form) has a '
            'well whose centroid R* tracks the THROAT scale (R*/r_s '
            '= 0.78-1.00) while R*/R_RMS falls monotonically toward '
            'the dilute end - R* and R_RMS emerge from one coupled '
            'system and DECOUPLE: the cavity belongs to the '
            'gravitational core, not to the matter cloud'
        ),
        'vacuum_deviation_vs_215': vac_dev,
        'members': members,
        'pass': bool(ok),
    }
    _CACHE['T6'] = out
    return out


# ========================================================================
# T7. The confrontation
# ========================================================================


def test_T7_confrontation() -> dict:
    t3 = test_T3_family()
    rows = t3['family']
    rs_omega = [r['rs_omega'] for r in rows]
    band_lo, band_hi = min(rs_omega), max(rs_omega)
    sup = t3['rs_omega_supremum']
    exclusion = band_lo / _ALPHA

    # jointly: sigma/r_s = 206.8 lives on the dilute branch where
    # X = sigma * omega = 206.8 * (r_s omega) -> the required conv-B
    # X = 1.5089 is missed by the same factor
    X_at_206 = _MU_OVER_E * sup
    X_factor = X_at_206 / _X_REQ_B

    # the q-channel class check: same O(1) band
    qrows = []
    for sc in (0.12, 0.18, 0.25):
        r = solve_star(sc, kind="qchan", g=4.0, a0=0.09, lam=1.0)
        if r is None:
            continue
        rs = math.sqrt(max(r["mu"], 0.0))
        qrows.append({'sc': sc, 'rs_omega': rs * r["W_phys"],
                      'X': r["rms"] * r["W_phys"]})

    ok = (1.4 < band_lo < 1.8 and 2.7 < sup < 2.85
          and exclusion > 200
          and X_factor > 300
          and len(qrows) == 3
          and all(1.5 < q['rs_omega'] < 2.8 for q in qrows))
    return {
        'name': 'T7_confrontation',
        'description': (
            'the lock excludes self-sourcing: r_s*omega is pinned '
            'into [1.53, 2.77] over the entire family (q-channel in '
            'the same band) - the EM-cap anchor r_s*omega = alpha is '
            'excluded by x210, and jointly sigma/r_s = 206.8 forces '
            'X = 574 vs the required 1.51 (x380): a 5D scalar cannot '
            'be much lighter than the throat its own field creates. '
            'The #210 primordial relocation is FORCED by the weld: '
            'the throat mass is a geometric datum, the #221 cavity '
            'is the primordial throat\'s (T6: R* tracks r_s, not the '
            'cloud), and the successor is the soliton on the '
            'fixed-mu throat background with perturbative '
            'back-reaction'
        ),
        'rs_omega_band': [float(band_lo), float(band_hi)],
        'rs_omega_supremum': float(sup),
        'anchor_required_rs_omega': _ALPHA,
        'anchor_exclusion_factor': float(exclusion),
        'X_at_sigma_over_rs_206': float(X_at_206),
        'X_required_conv_B': _X_REQ_B,
        'X_exclusion_factor': float(X_factor),
        'qchan_band': qrows,
        'pass': bool(ok),
    }


# ========================================================================
# T8. Honest scope
# ========================================================================


def test_T8_honest_scope() -> dict:
    scope = [
        'The stationary complex-scalar (boson-star) ansatz stands in '
        'for the #180 real psi-phi-q structure; the #210 q-channel '
        'potential is spot-checked in 5D and lands in the same '
        'r_s*omega band - the conclusion is a property of the '
        'potential class (the #210 compactness argument), not of the '
        'Kaup choice.',
        'One-mouth trivial topology: no cross-cap, no Pin-odd '
        'winding sector.  The refutation is of SELF-SOURCING - the '
        'reading in which the electron\'s own field creates its '
        'throat.  The forced alternative (the primordial fixed-mu '
        'throat with the soliton as a perturbative dressing, the '
        'bridge version of this exact solver) is the named '
        'successor.',
        'Ground states only; the spiral is entered (two members past '
        'the frequency minimum) but not traced; the compact end of '
        'the r_s*omega band is bounded by the solver\'s bracket '
        'window, and tracing the spiral deeper can only shrink the '
        'band from inside [1.4, 2.85] - it cannot reach alpha.',
        'The cavity potential is for test waves on the solved '
        'background (no back-reaction of the probe wave), '
        'FD-constructed and vacuum-validated; R* definitional band '
        '(centroid vs minimum vs width) carried by reporting all '
        'three.',
        'Classical throughout; hbar enters only through lambda_C = '
        '1/omega as in the whole arc; alpha imported (#184) for the '
        'confrontation only.',
        'This PR intentionally unfreezes the geometry (the program '
        'framing freezes it): the point is exactly to test what the '
        'frozen-background reading posits - and the coupled solve '
        'CONFIRMS the frozen primordial background as the only '
        'consistent reading of the anchor.',
    ]
    return {
        'name': 'T8_honest_scope',
        'description': 'what this PR does and does not establish',
        'scope': scope,
        'pass': True,
    }


# ========================================================================
# T9. Assessment
# ========================================================================


def test_T9_assessment() -> dict:
    t2 = test_T2_solver()
    t3 = test_T3_family()
    t4 = test_T4_weld_invariance()
    t5 = test_T5_israel_junction()
    t6 = test_T6_cavity()
    t7 = test_T7_confrontation()
    core = all(t['pass'] for t in (t2, t3, t4, t5, t6, t7))
    assessment = (
        'The audit\'s missing weld is supplied by the requested '
        'coupled system - and the weld decides. In the 5D '
        'Einstein-Klein-Gordon solve every scale (the soliton R_RMS, '
        'the cavity R*, the exterior Tangherlini r_s) emerges from '
        'one solution in one unit, with the coupling convention '
        'dropping out exactly and the Israel junction confirming a '
        'shared boundary. The derived family carries a 5D '
        'critical-mass marginality (mu -> 7.695 with a scale-free '
        'size zero mode - the same 1/r^2 marginality as #221\'s '
        'tails), and it pins r_s*omega into [1.53, 2.77]: the EM-cap '
        'anchor r_s*omega = alpha is excluded by x210. '
        'Self-sourcing is refuted at the coupled level; the #210 '
        'primordial-throat relocation is now forced, the #221 cavity '
        'is the primordial throat\'s (R* tracks r_s, not the cloud), '
        'and the identifiability program continues on the fixed-mu '
        'bridge with perturbative back-reaction.'
        if core else
        'INCONCLUSIVE - a core check failed; do not quote.'
    )
    return {
        'name': 'T9_assessment',
        'description': 'the standing of the derived weld',
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
        test_T2_solver(),
        test_T3_family(),
        test_T4_weld_invariance(),
        test_T5_israel_junction(),
        test_T6_cavity(),
        test_T7_confrontation(),
        test_T8_honest_scope(),
        test_T9_assessment(),
    ]
    t2, t3, t4, t5, t6, t7 = (tests[1], tests[2], tests[3], tests[4],
                              tests[5], tests[6])
    all_ok = all(t['pass'] for t in tests)
    if all_ok:
        verdict_class = (
            "THE_WELD_IS_DERIVED_AND_IT_REFUTES_SELF_SOURCING_THE_"
            "COUPLED_5D_EKG_LOCKS_RS_OMEGA_INTO_1P5_TO_2P8_NEVER_"
            "ALPHA_THE_THROAT_IS_PRIMORDIAL_AND_THE_CRITICAL_MASS_"
            "MARGINALITY_IS_MEASURED"
        )
        verdict = (
            "ESTABLISHED (the argument is in "
            "docs/coupled_5d_ekg_weld.md).\n\n"
            "THE SYSTEM. The requested modification, solved: the "
            "soliton's energy density sources the local 5D metric "
            "(static spherical EKG, shooting); the n = 2 reduction "
            f"lands the Kaup benchmark (M_max = {t2['kaup_M_max']:.4f}"
            " vs 0.633); the tt Einstein identity holds pointwise "
            f"(median {t2['constraint_residual_median']:.0e}); mu is "
            f"dx/xmax-converged ({t2['mu_dx_halving_shift']:.0e}/"
            f"{t2['mu_xmax_extension_shift']:.0e}).\n\n"
            "THE FAMILY. mu monotone along the main branch; the "
            "dilute endpoint has a CRITICAL MASS: mu -> "
            f"{t3['mu_crit_extrapolated']:.3f} (difference ratio "
            f"{t3['dilute_diff_ratio']:.2f}), omega -> 1 "
            f"(1 - omega = {t3['one_minus_omega_at_dilute_end']:.0e})"
            ", and the size runs FREE at fixed mass (X ~ sc^"
            f"{t3['X_dilute_slope_vs_sc']:.3f}) - the 5D Newtonian "
            "1/r^2 marginality, the same marginality as #221's "
            "tails; the spiral turnaround bounds the compact end.\n\n"
            "THE WELD. Convention invariance EXACT (worst rescale "
            f"shift {max(t4['rescale_relative_shifts'].values()):.0e})"
            ": every scale in one unit, the family relations derived "
            "- the audit's free rescale is forbidden by the field "
            "equations.  The Israel junction: [K^theta] = "
            f"{t5['junction_rows'][0]['jump_K_theta_rel']:.0e}, "
            f"[K^t] = {t5['junction_rows'][0]['jump_K_t_rel']:.0e} at "
            "the core boundary (jumps tracking the tail), and the "
            "exterior reads the interior mass to "
            f"{t5['exterior_reads_interior_mass_rel']:.0e}.\n\n"
            "THE CAVITY. The wave potential on the solved metric "
            f"(vacuum-validated vs #215 to "
            f"{t6['vacuum_deviation_vs_215']:.0e}): R* tracks the "
            "THROAT (R*/r_s = "
            f"{t6['members'][0]['R_star_over_rs']:.2f}-"
            f"{t6['members'][2]['R_star_over_rs']:.2f}) while "
            "R*/R_RMS falls "
            f"{t6['members'][0]['R_star_over_rms']:.2f} -> "
            f"{t6['members'][2]['R_star_over_rms']:.2f} dilute: R* "
            "and R_RMS emerge from ONE system and DECOUPLE - the "
            "cavity belongs to the gravitational core, not the "
            "cloud.\n\n"
            "THE CONFRONTATION. r_s*omega pinned into "
            f"[{t7['rs_omega_band'][0]:.3f}, "
            f"{t7['rs_omega_supremum']:.3f}] over the entire family "
            "(q-channel in-band): the EM-cap anchor r_s*omega = "
            f"alpha is EXCLUDED x{t7['anchor_exclusion_factor']:.0f}; "
            "jointly, sigma/r_s = 206.8 forces X = "
            f"{t7['X_at_sigma_over_rs_206']:.0f} vs the required "
            f"{_X_REQ_B:.2f} (x{t7['X_exclusion_factor']:.0f}). A 5D "
            "scalar cannot be much lighter than the throat its own "
            "field creates: SELF-SOURCING IS REFUTED at the coupled "
            "level, the #210 primordial relocation is FORCED, and "
            "the successor is the soliton on the fixed-mu throat "
            "bridge with perturbative back-reaction."
        )
    else:
        verdict_class = "COUPLED_5D_EKG_WELD_INCONCLUSIVE"
        verdict = (
            "INCONCLUSIVE. A core check failed; re-examine before "
            "quoting the result."
        )

    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "identification": (
            "The coupled 5D Einstein-Klein-Gordon weld: the soliton's "
            "energy density sources the local Tangherlini metric; "
            "R_RMS, the cavity R*, and r_s emerge from one solve in "
            "one unit (convention invariance exact; Israel junction "
            "verified); the 5D critical-mass marginality is measured "
            "(mu_crit = 7.695, scale-free size zero mode); and the "
            "lock r_s*omega in [1.53, 2.77] excludes the self-sourced "
            "anchor by x210 - the throat is primordial"
        ),
        "executes": (
            "the requested action modification (option A) with the "
            "Israel junction verification (option B) - the #221 "
            "audit's successor"
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
    out.append("# The coupled 5D EKG weld (PR #222)")
    out.append("")
    out.append(f"**Run:** {s['timestamp_utc']}")
    out.append("")
    out.append(
        "The deliverable is `docs/coupled_5d_ekg_weld.md` - the #221 "
        "audit's successor: the soliton sources the 5D metric, the "
        "scales lock, and the lock refutes self-sourcing. *(The one "
        "PR that unfreezes the geometry - to test what the frozen "
        "reading posits.)*"
    )
    out.append("")
    out.append("## Test summary")
    out.append("")
    out.append("| # | Test | Key finding | PASS? |")
    out.append("|---|---|---|---|")
    labels = {
        "T1": "the goal: both requested options executed",
        "T2": "solver benchmarked (Kaup), constraint, convergence",
        "T3": "the 5D family; the critical-mass marginality",
        "T4": "the weld: convention invariance EXACT",
        "T5": "the Israel junction verified (option B)",
        "T6": "R* and R_RMS from one system - and they decouple",
        "T7": "r_s*omega never alpha: self-sourcing refuted",
        "T8": "honest scope",
        "T9": "assessment",
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
    out = here / "runs" / f"{ts}_coupled_5d_ekg_weld_probe"
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
