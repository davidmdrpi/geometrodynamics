"""
The throat's spectral density: the flat bank replaced by the geometry's
own transmission spectrum (PR #215).

> Framing: QFT on the *fixed classical* throat geometry (geometry -> fields),
> not quantum gravity.

THE NAMED SUCCESSOR OF #214
---------------------------
#214 promoted the absorber to a degree of freedom but modelled its
spectrum as a FLAT oscillator bank - a shape chosen by hand.  The
physical absorber is the throat itself: a 5D Tangherlini mouth
(f(r) = 1 - r_h^2/r^2, the #168 identification) whose horizon is a
one-way membrane with a COMPUTABLE frequency response.  This probe
derives the absorber's spectral density from the throat geometry and
retires the bank:

1. THE SCATTERING PROBLEM.  A massless scalar mode phi(r) Y_l(S^3)
   e^{-iwt} on the Tangherlini exterior obeys, with psi = r^{3/2} phi
   and the tortoise coordinate dx = dr/f,
       psi_xx + [w^2 - V] psi = 0 ,
       V(r) = f * [ (l(l+2) + 3/4)/r^2 + (9/4) r_h^2/r^4 ] ,
       x(r) = r + (r_h/2) ln((r - r_h)/(r + r_h))  (analytic).
   Ingoing at the horizon, Hankel-matched at large r; the greybody
   transmission T_l(w) is flux-conserving to 1e-5 (1e-3 at the
   highest frequencies).

2. THE IR IS FIXED BY THE AREA THEOREM.  The universal low-frequency
   absorption theorem (any dimension): sigma_s(w -> 0) = A_horizon
   = 2 pi^2 r_h^3.  Measured: sigma/A_h = 1.012 at w r_h = 0.04,
   monotonically -> 1, and T_0 grows exactly as w^3 (slope 3.02).
   The spectral density's IR wing is A_h w^3/(4 pi) - GEOMETRY, not
   a model choice.

3. THE UV IS FIXED BY THE PHOTON SPHERE.  Null geodesics: b(r)^2 =
   r^2/f maximizes at r_c = sqrt(2) r_h with b_c = 2 r_h.  The
   T = 1/2 crossings sit at the eikonal w = (l+1)/b_c (1-4%), and
   above the barrier T -> 1: the throat is UV-BLACK - a perfect
   absorber for every mode that resolves it.

4. THE HORIZON IS THE CONTINUUM.  The tortoise coordinate diverges
   logarithmically at the horizon: a regulated interior of depth
   delta has length L ~ -(r_h/2) ln(delta) -> infinity, level spacing
   pi/L -> 0, recurrence time 2L -> infinity.  #214's order of limits
   (bank size N -> infinity FIRST) is not a limit to be taken here -
   the geometry has already taken it.  Classically the throat never
   revives.

5. THE WELD (parameter-free).  A cavity (the universe stand-in)
   terminated by the throat: the complex quasimode frequencies
   w_q - i gamma/2 (ingoing at the horizon, node at the wall) obey
   the transit law
       gamma = T(w_q) / (2 L_cav) ,   L_cav = x_wall - x_barrier ,
   measured at 0.6% / 1.3% / 8% over the first three modes of the
   ladder, with the ladder spacing pi/L_cav (4%).  The mode's epsilon
   is fully geometric: eps_n = w_n T(w_n)/(2 L_cav).  Applied to the
   S^3 tower via the #213 refocusing transit (one antipodal delivery
   per half-period pi R): gamma_n = T(w_n r_h)/(pi R) - IR-transparent
   as n^3, UV-black above n ~ b_c... the bank of #214 is retired with
   zero shape freedom left.

Tests:
  T1. Goal.
  T2. The scattering solver (flux conservation; regulator
      independence).
  T3. The IR anchor: the area theorem + the w^3 law.
  T4. The UV anchor: the photon sphere, the eikonal crossings,
      T -> 1 above the barrier.
  T5. The horizon is the continuum (tortoise log law; the #214
      order of limits taken by the geometry).
  T6. The weld: the quasimode ladder obeys the parameter-free
      transit law; the tower spectral density assembled.
  T7. Honest scope.
  T8. Assessment.

Verdict:
  THE_BANK_RETIRED_THE_SPECTRAL_DENSITY_IS_THE_THROATS_GREYBODY_
  IR_TRANSPARENT_BY_THE_AREA_THEOREM_UV_BLACK_ABOVE_THE_PHOTON_
  SPHERE_AND_THE_HORIZON_IS_THE_CONTINUUM
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar
from scipy.special import hankel1, hankel2

_CACHE: dict = {}
_RH = 1.0                       # throat radius (units)

# ========================================================================
# SECTION A - the Tangherlini scattering problem
# ========================================================================


def v_of_r(r, l: int):
    """psi_xx + [w^2 - V] psi = 0 potential, psi = r^{3/2} phi."""
    f = 1.0 - (_RH / r) ** 2
    return f * ((l * (l + 2) + 0.75) / r ** 2 + 2.25 * _RH ** 2 / r ** 4)


def x_of_r(r):
    """Tortoise coordinate, analytic: diverges -log at the horizon."""
    return r + (_RH / 2) * np.log((r - _RH) / (r + _RH))


def transmission(w: float, l: int, delta: float = 1e-7,
                 rtol: float = 1e-11) -> tuple:
    """Greybody transmission T_l(w): ingoing at the horizon, Hankel
    least-squares matched at large r (the far radius scales with w to
    keep the metric's long-range phase correction below the matching
    tolerance).  Returns (T, R)."""
    key = ('T', w, l, delta)
    if key in _CACHE:
        return _CACHE[key]
    r0 = _RH * (1 + delta)
    x0 = float(x_of_r(r0))
    r_far = max(60.0 / w, 50.0 * w, 30.0)
    x_far = float(x_of_r(r_far))

    def rhs(x, y):
        pr, pi_, qr, qi, r = y
        vv = v_of_r(r, l)
        f = 1 - (_RH / r) ** 2
        return [qr, qi, (vv - w ** 2) * pr, (vv - w ** 2) * pi_, f]

    y0 = [math.cos(-w * x0), math.sin(-w * x0),
          w * math.sin(-w * x0), -w * math.cos(-w * x0), r0]
    sol = solve_ivp(rhs, (x0, x_far), y0, rtol=rtol, atol=1e-13,
                    dense_output=True, method="DOP853")
    # least-squares match over spread points (a 2-point solve is
    # phase-degenerate when w*dx ~ n*pi)
    xs = x_far - np.linspace(0.0, 3.7, 8)
    rows, vals = [], []
    for xm in xs:
        pr, pi_, _, _, r = sol.sol(xm)
        z = w * r
        rows.append([math.sqrt(r) * hankel2(l + 1, z),
                     math.sqrt(r) * hankel1(l + 1, z)])
        vals.append(pr + 1j * pi_)
    (a, b), *_ = np.linalg.lstsq(np.array(rows), np.array(vals),
                                 rcond=None)
    T = float(np.pi * w / (2 * abs(a) ** 2))
    R = float(abs(b / a) ** 2)
    _CACHE[key] = (T, R)
    return T, R


def quasimode_shoot(w: complex, l: int = 0,
                    r_wall: float = 30.0) -> complex:
    """psi(x_wall) for ingoing-at-horizon data; zeros = quasimodes."""
    w = complex(w)
    r0 = _RH * (1 + 1e-7)
    x0 = float(x_of_r(r0))
    x_wall = float(x_of_r(r_wall))

    def rhs(x, y):
        psi = y[0] + 1j * y[1]
        q = y[2] + 1j * y[3]
        r = y[4]
        dq = (v_of_r(r, l) - w ** 2) * psi
        f = 1 - (_RH / r) ** 2
        return [q.real, q.imag, dq.real, dq.imag, f]

    p0 = np.exp(-1j * w * x0)
    q0 = -1j * w * p0
    sol = solve_ivp(rhs, (x0, x_wall),
                    [p0.real, p0.imag, q0.real, q0.imag, r0],
                    rtol=1e-11, atol=1e-13, method="DOP853")
    return complex(sol.y[0][-1], sol.y[1][-1])


def find_quasimode(w_guess: float, l: int = 0) -> complex:
    """Secant iteration in complex w for a quasimode near w_guess."""
    w0 = complex(w_guess)
    w1 = complex(w_guess) - 2e-4j
    f0, f1 = quasimode_shoot(w0, l), quasimode_shoot(w1, l)
    for _ in range(60):
        if abs(f1 - f0) < 1e-14:
            break
        w2 = w1 - f1 * (w1 - w0) / (f1 - f0)
        w0, f0, w1, f1 = w1, f1, w2, quasimode_shoot(w2, l)
        if abs(w1 - w0) < 1e-12:
            break
    return w1


# ========================================================================
# T1. Goal
# ========================================================================


def test_T1_goal() -> dict:
    return {
        'name': 'T1_goal',
        'description': (
            'Retire #214\'s flat bank: derive the absorber\'s spectral '
            'density from the throat geometry itself - the greybody '
            'transmission T_l(w) of the 5D Tangherlini mouth, with the '
            'IR fixed by the universal area theorem, the UV fixed by '
            'the photon sphere, the continuum limit taken by the '
            'horizon\'s infinite tortoise depth, and the whole welded '
            'to live cavity quasimodes by a parameter-free transit '
            'law.'
        ),
        'pass': True,
    }


# ========================================================================
# T2. The scattering solver
# ========================================================================


def test_T2_solver() -> dict:
    # (i) flux conservation R + T = 1 across the grid
    worst_low, worst_high = 0.0, 0.0
    for l in (0, 2):
        for w in (0.05, 0.2, 1.0):
            T, R = transmission(w, l)
            worst_low = max(worst_low, abs(R + T - 1))
        for w in (2.0, 6.0):
            T, R = transmission(w, l)
            worst_high = max(worst_high, abs(R + T - 1))

    # (ii) horizon-regulator independence
    t6, _ = transmission(0.5, 0, delta=1e-6)
    t8, _ = transmission(0.5, 0, delta=1e-8)
    reg_dep = abs(t6 - t8) / t8

    ok = worst_low < 5e-5 and worst_high < 5e-3 and reg_dep < 1e-4
    return {
        'name': 'T2_solver',
        'description': (
            'greybody solver: ingoing horizon data, Hankel-matched '
            'asymptotics; flux conservation and regulator independence'
        ),
        'flux_conservation_low_freq': float(worst_low),
        'flux_conservation_high_freq': float(worst_high),
        'horizon_regulator_dependence': float(reg_dep),
        'pass': bool(ok),
    }


# ========================================================================
# T3. The IR anchor: the area theorem
# ========================================================================


def test_T3_area_theorem() -> dict:
    if 'T3' in _CACHE:
        return _CACHE['T3']
    A_h = 2 * math.pi ** 2 * _RH ** 3

    # sigma_s(w) = 4 pi T_0 / w^3  ->  A_h  as  w -> 0
    ws = np.array([0.04, 0.06, 0.08])
    ratios = []
    for w in ws:
        T, _ = transmission(float(w), 0)
        ratios.append(4 * math.pi * T / w ** 3 / A_h)
    intercept = float(np.polyfit(ws, ratios, 1)[1])

    # the w^3 law
    ws2 = np.array([0.02, 0.03, 0.045, 0.06])
    Ts = np.array([transmission(float(w), 0)[0] for w in ws2])
    slope = float(np.polyfit(np.log(ws2), np.log(Ts), 1)[0])

    ok = (abs(ratios[0] - 1) < 0.02
          and ratios[0] < ratios[1] < ratios[2]
          and abs(intercept - 1) < 0.03
          and abs(slope - 3.0) < 0.05)
    out = {
        'name': 'T3_area_theorem',
        'description': (
            'the universal low-frequency theorem: sigma_s(w->0) = '
            'A_horizon = 2 pi^2 r_h^3, with T_0 ~ w^3 - the spectral '
            'density\'s IR wing is pure geometry'
        ),
        'sigma_over_area': {f"{w:.2f}": float(r)
                            for w, r in zip(ws, ratios)},
        'extrapolated_intercept': intercept,
        'omega_cubed_slope': slope,
        'pass': bool(ok),
    }
    _CACHE['T3'] = out
    return out


# ========================================================================
# T4. The UV anchor: the photon sphere
# ========================================================================


def test_T4_photon_sphere() -> dict:
    if 'T4' in _CACHE:
        return _CACHE['T4']
    # (i) null geodesics: b(r)^2 = r^2/f maximizes the potential at
    # r_c = sqrt(2) r_h with critical impact parameter b_c = 2 r_h
    res = minimize_scalar(lambda r: -(1 - (_RH / r) ** 2) / r ** 2,
                          bounds=(1.001, 10.0), method="bounded")
    r_c = float(res.x)
    b_c = r_c / math.sqrt(1 - (_RH / r_c) ** 2)

    # (ii) eikonal crossings: T = 1/2 at w ~ (l+1)/b_c; WKB ratio -> 1
    eik, wkb = {}, {}
    for l in (1, 2, 3):
        vres = minimize_scalar(lambda r: -v_of_r(r, l),
                               bounds=(1.001, 10.0), method="bounded")
        w_wkb = math.sqrt(-vres.fun)
        w_half = brentq(lambda w: transmission(w, l)[0] - 0.5,
                        0.4, 2.5 * w_wkb, xtol=1e-3)
        eik[l] = w_half / ((l + 1) / b_c)
        wkb[l] = w_half / w_wkb

    # (iii) UV black: T -> 1 above the barrier
    t3, _ = transmission(3.0, 0)
    t6, _ = transmission(6.0, 0)

    ok = (abs(r_c - math.sqrt(2)) < 1e-6 and abs(b_c - 2.0) < 1e-6
          and all(abs(e - 1) < 0.05 for e in eik.values())
          and wkb[1] < wkb[2] < wkb[3] < 1.0
          and abs(t3 - 1) < 2e-3 and abs(t6 - 1) < 5e-3)
    out = {
        'name': 'T4_photon_sphere',
        'description': (
            'r_c = sqrt(2) r_h, b_c = 2 r_h; T = 1/2 at the eikonal '
            '(l+1)/b_c; T -> 1 above the barrier: the throat is '
            'UV-black for every mode that resolves it'
        ),
        'photon_sphere_radius': r_c,
        'critical_impact_parameter': b_c,
        'eikonal_ratio': {str(k): float(v) for k, v in eik.items()},
        'wkb_ratio': {str(k): float(v) for k, v in wkb.items()},
        'T_above_barrier': {'3.0': t3, '6.0': t6},
        'pass': bool(ok),
    }
    _CACHE['T4'] = out
    return out


# ========================================================================
# T5. The horizon is the continuum
# ========================================================================


def test_T5_horizon_continuum() -> dict:
    # regulated interior depth: L(delta) = x(r_wall) - x(r_h(1+delta))
    r_wall = 30.0
    deltas = np.array([1e-2, 1e-4, 1e-6, 1e-8])
    Ls = np.array([x_of_r(r_wall) - x_of_r(_RH * (1 + d))
                   for d in deltas])
    # analytic: dL/d(ln delta) = -r_h/2 exactly as delta -> 0
    coef = float(np.polyfit(np.log(deltas), Ls, 1)[0])
    spacing = math.pi / Ls[-1]
    t_rec = 2 * Ls[-1]

    ok = (abs(coef + _RH / 2) < 1e-3
          and np.all(np.diff(Ls) > 0)
          and spacing < math.pi / Ls[0])
    return {
        'name': 'T5_horizon_continuum',
        'description': (
            'the tortoise depth diverges as -(r_h/2) ln(delta): level '
            'spacing pi/L -> 0, recurrence time 2L -> infinity - the '
            '#214 order of limits (N -> infinity first) is taken by '
            'the geometry, and classically the throat never revives'
        ),
        'log_coefficient': coef,
        'log_coefficient_analytic': -_RH / 2,
        'depths': {f"{d:.0e}": float(L) for d, L in zip(deltas, Ls)},
        'spacing_at_finest_regulator': float(spacing),
        'recurrence_time_at_finest_regulator': float(t_rec),
        'pass': bool(ok),
    }


# ========================================================================
# T6. The weld: quasimodes obey the parameter-free transit law
# ========================================================================


def test_T6_quasimode_weld() -> dict:
    if 'T6' in _CACHE:
        return _CACHE['T6']
    r_wall = 30.0
    # cavity length: wall to the barrier peak (the sub-barrier region
    # is the drain, not the cavity)
    vres = minimize_scalar(lambda r: -v_of_r(r, 0),
                           bounds=(1.001, 5.0), method="bounded")
    L_cav = float(x_of_r(r_wall) - x_of_r(vres.x))

    # locate the ladder by scanning |shoot| on the real axis
    ws = np.linspace(0.2, 0.46, 105)
    vals = np.abs([quasimode_shoot(w) for w in ws])
    guesses = [float(ws[i]) for i in range(1, len(ws) - 1)
               if vals[i] < vals[i - 1] and vals[i] < vals[i + 1]][:3]

    modes = []
    for g in guesses:
        wq = find_quasimode(g)
        gam = -2 * wq.imag
        T, _ = transmission(wq.real, 0)
        transit = T / (2 * L_cav)
        modes.append({
            'w': float(wq.real),
            'gamma': float(gam),
            'T_at_w': float(T),
            'transit_prediction': float(transit),
            'ratio': float(gam / transit),
            'eps_geometric': float(wq.real * gam),
        })
    spacings = np.diff([m['w'] for m in modes])
    spacing_ratio = [float(s / (math.pi / L_cav)) for s in spacings]

    # the tower spectral density (r_h/R = 0.1; per-pass absorption
    # T(w_n r_h) delivered once per half-period by the #213 refocus)
    rh_over_R = 0.1
    tower = {}
    for n in (1, 2, 4, 8, 12, 16, 24):
        T, _ = transmission(n * rh_over_R, 0)
        tower[n] = float(min(T, 1.0))
    tvals = [tower[n] for n in (1, 2, 4, 8, 12, 16, 24)]
    # IR: n^3 growth between n = 1 and 2
    ir_ratio = tower[2] / tower[1]

    ok = (len(modes) == 3
          and abs(modes[0]['ratio'] - 1) < 0.02
          and all(abs(m['ratio'] - 1) < 0.12 for m in modes)
          and all(abs(s - 1) < 0.06 for s in spacing_ratio)
          and all(tvals[i] < tvals[i + 1] + 1e-3
                  for i in range(len(tvals) - 1))
          and 7.5 < ir_ratio < 9.5
          and tower[1] < 0.01 and tower[24] > 0.99)
    out = {
        'name': 'T6_quasimode_weld',
        'description': (
            'the cavity quasimode widths obey gamma = T(w)/(2 L_cav) '
            'with NO free parameter; the ladder spacing is pi/L_cav; '
            'the tower density is IR-transparent (n^3) and UV-black'
        ),
        'L_cav': L_cav,
        'modes': modes,
        'ladder_spacing_over_pi_L': spacing_ratio,
        'tower_per_pass_absorption': {str(k): v
                                      for k, v in tower.items()},
        'tower_ir_cubed_ratio': float(ir_ratio),
        'pass': bool(ok),
    }
    _CACHE['T6'] = out
    return out


# ========================================================================
# T7. Honest scope
# ========================================================================


def test_T7_honest_scope() -> dict:
    scope = [
        'Classical scalar greybody: no Hawking flux, no spin/tensor '
        'channels (the photon\'s Lorentz structure rides on #46; its '
        'greybody differs at O(1) and is a follow-up).',
        'The tower rate uses the #213 refocusing transit (one '
        'antipodal delivery per half-period, s-channel): the zonal 1D '
        'reduction the #214 cavity weld validates exactly; the full '
        'S^3-with-throat matched asymptotics (3D flux spreading, '
        'l-mixing at the caustic) is the named successor.',
        'r_h/R is a parameter here (0.1 in the table); the #210 '
        'anchor r_s ~ alpha lambda_C supplies the physical value - '
        'laboratory modes sit FAR into the IR-transparent wing '
        '(consistency: everyday fields propagate freely; the throat '
        'blackens only at its own Compton-edge scale).',
        'Absorption is classical one-way membrane physics; what the '
        'horizon does with the energy (the interior, re-emission, '
        'the antipodal identification of the mouths) is #168/#200 '
        'territory, untouched.',
        'Matching accuracy degrades to ~1e-3 at the highest '
        'frequencies scanned (the long-range metric phase); all '
        'quoted UV numbers carry that error bar.',
        'Frozen geometry, no backreaction, throughout.',
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
    t2 = test_T2_solver()
    t3 = test_T3_area_theorem()
    t4 = test_T4_photon_sphere()
    t5 = test_T5_horizon_continuum()
    t6 = test_T6_quasimode_weld()
    core = all(t['pass'] for t in (t2, t3, t4, t5, t6))
    assessment = (
        'The flat bank is retired: the absorber\'s spectral density '
        'is now the throat\'s own greybody transmission, with both '
        'wings pinned by theorems of the geometry - the IR by the '
        'universal area theorem (sigma_s -> A_h, T ~ w^3), the UV by '
        'the photon sphere (T -> 1 above (l+1)/b_c, b_c = 2 r_h) - '
        'the continuum limit #214 had to take by hand supplied by the '
        'horizon\'s infinite tortoise depth, and the whole welded to '
        'live dynamics by a parameter-free transit law that the '
        'quasimode ladder obeys at the percent level. Nothing about '
        'the absorber is chosen anymore: its spectrum, its address '
        '(#214), and its epsilon are all read off the frozen bulk.'
        if core else
        'INCONCLUSIVE - a core check failed; do not quote.'
    )
    return {
        'name': 'T8_assessment',
        'description': 'the standing of the geometric spectral density',
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
        test_T3_area_theorem(),
        test_T4_photon_sphere(),
        test_T5_horizon_continuum(),
        test_T6_quasimode_weld(),
        test_T7_honest_scope(),
        test_T8_assessment(),
    ]
    t3, t4, t5, t6 = tests[2], tests[3], tests[4], tests[5]
    all_ok = all(t['pass'] for t in tests)
    if all_ok:
        verdict_class = (
            "THE_BANK_RETIRED_THE_SPECTRAL_DENSITY_IS_THE_THROATS_"
            "GREYBODY_IR_TRANSPARENT_BY_THE_AREA_THEOREM_UV_BLACK_"
            "ABOVE_THE_PHOTON_SPHERE_AND_THE_HORIZON_IS_THE_CONTINUUM"
        )
        m = t6['modes']
        verdict = (
            "DERIVED (the argument is in "
            "docs/throat_spectral_density.md).\n\n"
            "THE SPECTRUM FROM GEOMETRY. The 5D Tangherlini greybody "
            "T_l(w), flux-conserving across the band: the IR "
            "wing is the universal area theorem - sigma_s/A_h = "
            f"{t3['sigma_over_area']['0.04']:.3f} at w r_h = 0.04, "
            f"extrapolating to {t3['extrapolated_intercept']:.3f}, "
            f"with T ~ w^3 (slope {t3['omega_cubed_slope']:.2f}) - "
            "and the UV wing is the photon sphere: r_c = sqrt(2) r_h, "
            "b_c = 2 r_h, T = 1/2 at the eikonal (l+1)/b_c (ratios "
            f"{t4['eikonal_ratio']['1']:.2f}/"
            f"{t4['eikonal_ratio']['2']:.2f}/"
            f"{t4['eikonal_ratio']['3']:.2f}), T -> 1 above.\n\n"
            "THE CONTINUUM. The tortoise depth grows as -(r_h/2) ln "
            f"delta (coefficient {t5['log_coefficient']:.4f}): level "
            "spacing -> 0, recurrence -> infinity - the #214 order of "
            "limits is a property of the horizon, not a choice.\n\n"
            "THE WELD. Cavity quasimodes obey the parameter-free "
            "transit law gamma = T(w)/(2 L_cav): ratios "
            f"{m[0]['ratio']:.3f}/{m[1]['ratio']:.3f}/"
            f"{m[2]['ratio']:.3f} down the ladder, spacing pi/L_cav "
            f"to {max(abs(s - 1) for s in t6['ladder_spacing_over_pi_L']):.0%}. "
            "The tower density (r_h/R = 0.1): per-pass absorption "
            f"{t6['tower_per_pass_absorption']['1']:.1e} at n = 1 "
            "rising as n^3 (ratio "
            f"{t6['tower_ir_cubed_ratio']:.1f}/8 expected) to "
            f"{t6['tower_per_pass_absorption']['24']:.3f} at n = 24: "
            "IR-transparent, UV-black. Nothing about the absorber is "
            "chosen anymore - spectrum, address, and epsilon are all "
            "read off the frozen bulk."
        )
    else:
        verdict_class = "THROAT_SPECTRAL_DENSITY_INCONCLUSIVE"
        verdict = (
            "INCONCLUSIVE. A core check failed; re-examine before "
            "quoting the derivation."
        )

    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "identification": (
            "The throat's spectral density: #214's flat bank replaced "
            "by the 5D Tangherlini greybody transmission - IR pinned "
            "by the universal area theorem, UV by the photon sphere, "
            "the continuum by the horizon's infinite tortoise depth, "
            "and the quasimode ladder obeying the parameter-free "
            "transit law gamma = T(w)/(2 L_cav)"
        ),
        "executes": (
            "the named #214 successor: the absorber spectrum derived "
            "from the throat geometry itself"
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
    out.append("# The throat's spectral density (PR #215)")
    out.append("")
    out.append(f"**Run:** {s['timestamp_utc']}")
    out.append("")
    out.append(
        "The deliverable is `docs/throat_spectral_density.md` - the "
        "flat bank replaced by the geometry's own transmission "
        "spectrum. *(QFT on the fixed classical throat geometry, not "
        "quantum gravity.)*"
    )
    out.append("")
    out.append("## Test summary")
    out.append("")
    out.append("| # | Test | Key finding | PASS? |")
    out.append("|---|---|---|---|")
    labels = {
        "T1": "the goal: the bank retired",
        "T2": "the greybody solver, flux-conserving",
        "T3": "IR pinned by the universal area theorem",
        "T4": "UV pinned by the photon sphere",
        "T5": "the horizon is the continuum",
        "T6": "quasimodes obey the parameter-free transit law",
        "T7": "honest scope",
        "T8": "assessment",
    }
    for t in s["tests"]:
        p = "**PASS**" if t["pass"] else "**FAIL**"
        pre = t["name"][:2]
        out.append(f"| {pre} | `{t['name']}` | {labels.get(pre,'-')} | {p} |")
    out.append("")
    t6 = s["tests"][5]
    if 'modes' in t6:
        out.append("## The quasimode ladder vs the transit law")
        out.append("")
        out.append("| w_q | gamma (measured) | T(w)/(2L) (predicted) | ratio |")
        out.append("|---:|---:|---:|---:|")
        for m in t6['modes']:
            out.append(f"| {m['w']:.5f} | {m['gamma']:.3e} | "
                       f"{m['transit_prediction']:.3e} | {m['ratio']:.3f} |")
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
    out = here / "runs" / f"{ts}_throat_spectral_density_probe"
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
