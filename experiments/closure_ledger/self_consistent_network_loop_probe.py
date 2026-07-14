"""
The self-consistent network loop: the full two-mouth transfer system
solved, and the effective Green function DERIVED before comparison
with I+- (PR #217).

> Framing: QFT on the *fixed classical* throat geometry (geometry -> fields),
> not quantum gravity.

WHAT #216 LEFT OPEN
-------------------
#216 built the two-port throat and showed the network's projected
kernel is advanced - but treated the loop at first order: one
traversal, closure conditions read off the engine's bookkeeping, the
I+- assembly postulated with the network weight.  The decisive step is
to construct the FULL two-mouth transfer system and solve the loop
SELF-CONSISTENTLY: the returned confirmation re-enters the network
(the closed universe refocuses it onto mouth A again), so the field at
the crossing obeys the CTC constraint

    F = F_0 + Lambda F      =>     G_eff(w) = g / (1 - Lambda(w)) ,

with everything in Lambda derived.  Doing this carefully forces three
corrections/completions of #216:

1. CLOCK-RATE-CORRECT TRAVERSAL.  The throat is static in its own
   proper time, so ports and interior loops run at w_tau = w/rate_A,
   and the emergent global frequency is w * rate_B / rate_A - a
   slow-clock (deep-well) exit mouth REDSHIFTS the wave.  (This
   corrects an inverted rate ratio shipped in #216, exercised only at
   equal rates.)  Elastic confirmation - the return matching the
   source mode - then REQUIRES rate_A = rate_B at traversal epoch:
   state closure selects rate-matched networks, only the aged-in
   offset Delta_BA surviving.

2. THE VALUE-TRANSPORT LOOP EIGENVALUE.  The self-consistency of the
   field is governed by the returned VALUE per emitted value.  Summing
   the echo train with its k-dependent arrival displacements,

       Lambda(w) = t_net(w_tau) * deco * e^{+i w D_loop} ,
       D_loop = d_A + d_B + tau_glob + Delta_BA ,

   validated from first principles against a RING spectrum (transfer-
   matrix ring modes at 2.7429/3.6000 vs the Lambda = 1 prediction
   2.7390/3.6071 - 0.2%, the residual the O(|r|) backscatter
   splitting; the opposite phase convention misses by 0.16 and is
   excluded).  This CORRECTS the #216 closure bookkeeping: at time
   closure the carrier closes on the throat's own scattering phase,
   arg t_net = 0 (mod 2 pi) - the engine-convention comb of #216
   tracked retro_phase_match's time-phase amplitude, which coincides
   with value transport on the tower's exterior legs (e^{-2 pi i n}=1)
   but not for the transit phase.  The deform knob (phi = w * delta)
   is unchanged.

3. WIGNER-DELAY-CORRECT CLOSURE and THE TRANSACTION POINTS.  Group
   closure demands D_loop = -tau_W (the composite Wigner delay), and
   carrier closure arg Lambda = 0 (mod 2 pi).  Solving BOTH:
   above the barrier the throat phase is nearly flat, and the mouths'
   geometric transfer phase theta supplies the carrier tuning (the
   packet then lands ON the crossing: peak 0.006 vs 0.175
   uncorrected); below the barrier the equations are solved by the
   throat alone - and the solution sits ON the interior Fabry-Perot
   resonance (|t_net|^2 = 0.969): COMPLETED TRANSACTIONS LIVE ON THE
   THROAT'S RESONANCE COMB, derived rather than observed.

THE EFFECTIVE GREEN FUNCTION, THEN I+-.  G_eff = g/(1 - Lambda) is
derived as the resolvent of the explicit transfer system (machine
identity) and only THEN compared with the transactional assembly:
its first-winding truncation is exactly #216's K1 = I+ + Lambda I-;
the resummation renormalizes the confirmation weight to
Lambda/(1 - Lambda); the line shape at a completed transaction is
QUARTIC (group closure makes the loop phase stationary - an
anomalously flat resonance), Lorentzian when only the carrier closes,
with widths predicted by (1 - |Lambda|); and passivity |Lambda| <= 1
(two-port unitarity) makes every self-consistent solution stable or
marginal - the Novikov fixed point cannot run away, and the completed
transaction (Lambda -> 1) is the #213 epsilon -> 0 coherent limit,
with epsilon_eff = confirmation deficit per loop time = the #214
absorber damping in its transactional form.

Tests:
  T1. Goal.
  T2. Clock-rate-correct traversal (the corrected redshift; the exact
      clock composition of the exit time; the elastic selection rule).
  T3. The full transfer system (matrix resolvent = closed form; the
      nested interior x winding double resummation; the ring-spectrum
      validation of the value-transport eigenvalue; the #216
      convention correction quantified).
  T4. Wigner-correct closure and the transaction points (above-barrier
      point via the mouth phase; below-barrier point ON the interior
      resonance; packets land on the crossing, phase-aligned).
  T5. Source and mouth state evolution (input-output/coupled-mode:
      kappa_tot vs the t_net linewidth and the resonant Wigner delay;
      the source fixed-point iteration converging at rate |Lambda|;
      the cavity build-up ledger).
  T6. The effective Green function vs I+- (quartic and Lorentzian
      line shapes with predicted widths; passivity sweep; the O(Lambda)
      truncation = the #216 assembly; the resummed weight; the
      epsilon unification).
  T7. Honest scope.
  T8. Assessment.

Verdict:
  THE_LOOP_SOLVED_SELF_CONSISTENTLY_G_EFF_DERIVED_COMPLETED_
  TRANSACTIONS_SIT_ON_THE_THROAT_RESONANCES_AND_THE_NOVIKOV_FIXED_
  POINT_IS_PASSIVE
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from geometrodynamics.transaction import (
    MouthPort,
    NetworkMouth,
    NetworkThroat,
    closure_offset,
    effective_green,
    emergent_frequency,
    loop_eigenvalue,
    network_confirmation,
    traverse_throat,
)

_CACHE: dict = {}
_RH = 1.0
_D_A = math.pi
_D_B = math.pi

# ========================================================================
# SECTION A - the two-sided greybody (PR #215/#216 solver)
# ========================================================================


def _v_of_r(r, l: int = 0):
    f = 1.0 - (_RH / r) ** 2
    return f * ((l * (l + 2) + 0.75) / r ** 2 + 2.25 * _RH ** 2 / r ** 4)


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


def greybody_exterior(w: float) -> tuple:
    key = ('ge', w)
    if key in _CACHE:
        return _CACHE[key]
    r0 = _RH * (1 + 1e-7)
    x0 = float(_x_of_r(r0))
    r_far = max(60.0 / w, 50.0 * w, 30.0)
    x_far = float(_x_of_r(r_far))
    y0 = [math.cos(-w * x0), math.sin(-w * x0),
          w * math.sin(-w * x0), -w * math.cos(-w * x0), r0]
    out = _solve_and_match(w, x0, x_far, y0, (x0, x_far), +1)
    _CACHE[key] = out
    return out


def greybody_interior(w: float) -> tuple:
    key = ('gi', w)
    if key in _CACHE:
        return _CACHE[key]
    r0 = _RH * (1 + 1e-7)
    x0 = float(_x_of_r(r0))
    r_far = max(60.0 / w, 50.0 * w, 30.0)
    x_far = float(_x_of_r(r_far))
    p0 = np.exp(1j * w * x_far)
    q0 = 1j * w * p0
    y0 = [p0.real, p0.imag, q0.real, q0.imag, r_far]
    out = _solve_and_match(w, x0, x_far, y0, (x_far, x0), -1)
    _CACHE[key] = out
    return out


def tangherlini_port() -> MouthPort:
    return MouthPort(
        t=lambda w: greybody_exterior(w)[0],
        r_out=lambda w: greybody_exterior(w)[1],
        r_in=lambda w: greybody_interior(w)[1],
    )


def make_throat(tau_th: float, delta: float,
                theta: float = 0.0) -> NetworkThroat:
    A = NetworkMouth("A", psi=math.pi, link_id="L1", clock_offset=0.0)
    B = NetworkMouth("B", psi=math.pi, link_id="L1", clock_offset=delta,
                     transfer_phase=theta)
    return NetworkThroat(A, B, tau_th=tau_th,
                         port_A=tangherlini_port(),
                         port_B=tangherlini_port())


def tau_wigner(w: float, tau_th: float, h: float = 1e-3) -> float:
    """Composite Wigner delay d arg t_net / dw."""
    th = make_throat(tau_th, 0.0)
    return float((np.angle(th.t_AB(w + h))
                  - np.angle(th.t_AB(w - h))) / (2 * h))


# ========================================================================
# T1. Goal
# ========================================================================


def test_T1_goal() -> dict:
    return {
        'name': 'T1_goal',
        'description': (
            'Construct the full two-mouth transfer system and solve '
            'the loop self-consistently: both mouth barriers, '
            'clock-rate-correct traversal, Wigner-delay-correct packet '
            'closure, repeated winding resummation, source and mouth '
            'state evolution - and the resulting effective Green '
            'function DERIVED (as the resolvent of the transfer '
            'system, validated against a first-principles ring '
            'spectrum) before any comparison with I+-.'
        ),
        'pass': True,
    }


# ========================================================================
# T2. Clock-rate-correct traversal
# ========================================================================


def test_T2_clock_rates() -> dict:
    # (i) the corrected redshift: a slow-clock exit mouth redshifts
    A = NetworkMouth("A", psi=math.pi, link_id="L")
    B_slow = NetworkMouth("B", psi=math.pi, link_id="L", clock_rate=0.5)
    th_slow = NetworkThroat(A, B_slow, tau_th=0.8,
                            port_A=tangherlini_port(),
                            port_B=tangherlini_port())
    w_out = emergent_frequency(th_slow, 2.4)
    redshift_ok = abs(w_out - 1.2) < 1e-12          # NOT 4.8 (the #216 bug)

    # (ii) the exact clock composition of the exit time (equal rates
    # rho: global transit tau_th / rho)
    A2 = NetworkMouth("A", psi=math.pi, link_id="L", clock_rate=0.5)
    B2 = NetworkMouth("B", psi=math.pi, link_id="L", clock_rate=0.5,
                      clock_offset=-9.0)
    th2 = NetworkThroat(A2, B2, tau_th=0.8,
                        port_A=tangherlini_port(),
                        port_B=tangherlini_port())
    leg = traverse_throat(th2, 1.3, t_entry=2.0)
    exit_ok = (abs(leg.t_end - (2.0 + 0.8 / 0.5 - 9.0)) < 1e-12
               and abs(leg.local_duration - 0.8) < 1e-12)

    # (iii) the elastic selection rule: unequal rates break state
    # closure - the return frequency misses the source mode
    mismatch = abs(emergent_frequency(th_slow, 2.4) - 2.4)
    elastic = abs(emergent_frequency(th2, 2.4) - 2.4)   # equal rates

    ok = redshift_ok and exit_ok and mismatch > 1.0 and elastic < 1e-12
    return {
        'name': 'T2_clock_rates',
        'description': (
            'the throat runs at w_tau = w/rate_A; the emergent '
            'frequency is w rate_B/rate_A (redshift out of the '
            'slow-clock well - correcting the inverted #216 ratio); '
            'elastic confirmation forces rate-matched mouths'
        ),
        'emergent_frequency_slow_exit': float(w_out),
        'corrected_from_216_value': 4.8,
        'exit_time_composition_ok': bool(exit_ok),
        'state_closure_mismatch_unequal_rates': float(mismatch),
        'state_closure_mismatch_equal_rates': float(elastic),
        'pass': bool(ok),
    }


# ========================================================================
# T3. The full transfer system and the loop eigenvalue
# ========================================================================


def test_T3_transfer_system() -> dict:
    if 'T3' in _CACHE:
        return _CACHE['T3']
    tau = 0.8
    delta_off = closure_offset(_D_A, _D_B, tau) + 0.13
    th = make_throat(tau, delta_off)

    # (i) the signal-flow matrix (both barriers explicit) resolvent
    # equals the closed form 1/(1 - Lambda)
    w = 2.7
    t, _ = greybody_exterior(w)
    _, r_in = greybody_interior(w)
    M = np.zeros((4, 4), dtype=complex)
    M[1, 0] = t
    M[2, 1] = r_in ** 2 * np.exp(2j * w * tau)
    M[1, 2] = 1.0
    M[3, 1] = t
    M[0, 3] = np.exp(1j * w * (_D_A + _D_B + tau + delta_off))
    x = np.linalg.solve(np.eye(4) - M, np.array([1, 0, 0, 0], complex))
    lam = loop_eigenvalue(th, w, _D_A, _D_B)
    matrix_err = abs(x[0] - 1.0 / (1.0 - lam))

    # (ii) the nested double resummation: interior loops x windings
    w2 = 0.5
    th_lo = make_throat(tau, delta_off)
    lam_lo = loop_eigenvalue(th_lo, w2, _D_A, _D_B)
    t2, _ = greybody_exterior(w2)
    _, ri2 = greybody_interior(w2)
    inner = sum((ri2 ** 2 * np.exp(2j * w2 * tau)) ** k
                for k in range(80))
    lam_trunc = (t2 * t2 * inner
                 * np.exp(1j * w2 * (_D_A + _D_B + tau + delta_off)))
    double = sum(lam_trunc ** j for j in range(80))
    nested_err = abs(double - 1.0 / (1.0 - lam_lo))

    # (iii) the ring-spectrum validation of the value-transport
    # eigenvalue: transfer-matrix ring modes vs the Lambda = 1 comb
    def barrier_M(wq, rev=False):
        tq, ro = greybody_exterior(wq)
        _, riq = greybody_interior(wq)
        if rev:
            ro, riq = riq, ro
        return np.array([[(tq * tq - ro * riq) / tq, riq / tq],
                         [-ro / tq, 1 / tq]])

    L_ext = _D_A + _D_B

    def ring_tr(wq):
        P_ext = np.diag([np.exp(1j * wq * L_ext),
                         np.exp(-1j * wq * L_ext)])
        P_int = np.diag([np.exp(1j * wq * tau), np.exp(-1j * wq * tau)])
        return complex(np.trace(
            P_ext @ barrier_M(wq) @ P_int @ barrier_M(wq, True))).real

    ws = np.linspace(2.5, 4.2, 69)
    trs = np.array([ring_tr(wq) for wq in ws])
    modes = []
    for i in range(1, len(ws) - 1):
        if trs[i] > trs[i - 1] and trs[i] > trs[i + 1]:
            # parabolic sub-grid refinement of the trace maximum
            d = (trs[i - 1] - trs[i + 1]) / (
                2 * (trs[i - 1] - 2 * trs[i] + trs[i + 1]))
            modes.append(float(ws[i] + d * (ws[1] - ws[0])))

    def pred_phase(wq, sign):
        th0 = make_throat(tau, 0.0)
        return float(np.angle(th0.t_AB(wq)
                              * np.exp(sign * 1j * wq * (L_ext + tau))))

    def comb_roots(sign):
        vals = [pred_phase(wq, sign) for wq in ws]
        roots = []
        for i in range(len(ws) - 1):
            if (vals[i] * vals[i + 1] < 0
                    and abs(vals[i + 1] - vals[i]) < math.pi):
                roots.append(brentq(lambda q: pred_phase(q, sign),
                                    ws[i], ws[i + 1], xtol=1e-8))
        return roots

    good = comb_roots(+1)
    bad = comb_roots(-1)
    ring_err = max(min(abs(m - r) for r in good) for m in modes)
    wrong_conv = min(min(abs(m - r) for r in bad) for m in modes)

    # (iv) the #216 convention correction, quantified: the engine's
    # loop amplitude differs from the value-transport eigenvalue by
    # exactly e^{-i w (d_A + d_B + tau_glob)}
    tr = network_confirmation(th, w, 0.0, _D_A, _D_B)
    d_loop = tr.t_return - tr.t_emit
    engine_lam = tr.amp * np.exp(1j * w * d_loop)
    conv_err = abs(engine_lam
                   - lam * np.exp(-1j * w * (_D_A + _D_B + tau)))

    ok = (matrix_err < 1e-12 and nested_err < 1e-10
          and len(modes) >= 2 and ring_err < 0.01 and wrong_conv > 0.1
          and conv_err < 1e-12)
    out = {
        'name': 'T3_transfer_system',
        'description': (
            'the explicit transfer system resolves to G = 1/(1-Lambda) '
            'exactly; interior x winding loops resum nested; the ring '
            'spectrum validates Lambda = t_net e^{+i w D_loop} from '
            'first principles and excludes the opposite convention'
        ),
        'matrix_resolvent_error': float(matrix_err),
        'nested_resummation_error': float(nested_err),
        'ring_modes': modes,
        'ring_prediction_error': float(ring_err),
        'wrong_convention_miss': float(wrong_conv),
        'engine_convention_offset_error': float(conv_err),
        'pass': bool(ok),
    }
    _CACHE['T3'] = out
    return out


# ========================================================================
# T4. Wigner-correct closure and the transaction points
# ========================================================================


def test_T4_transaction_points() -> dict:
    if 'T4' in _CACHE:
        return _CACHE['T4']
    # (a) ABOVE THE BARRIER (w0 = 3, a tower mode): group closure
    # fixes Delta* = -(dA+dB+tau) - tau_W; the throat phase is nearly
    # flat there, so the carrier is closed by the mouths' geometric
    # transfer phase theta*
    w0, tau = 3.0, 0.8
    tw = tau_wigner(w0, tau)
    th0 = make_throat(tau, 0.0)
    theta_star = float(-np.angle(th0.t_AB(w0) * np.exp(-1j * w0 * tw)))
    delta_star = closure_offset(_D_A, _D_B, tau) - tw
    th_star = make_throat(tau, delta_star, theta_star)
    lam_star = loop_eigenvalue(th_star, w0, _D_A, _D_B)

    def packet_peak(th, w_c, sig=0.3, npts=41):
        wsg = np.linspace(w_c - 4 * sig, w_c + 4 * sig, npts)
        c = np.exp(-0.5 * ((wsg - w_c) / sig) ** 2)
        lams = np.array([loop_eigenvalue(th, float(wq), _D_A, _D_B)
                         for wq in wsg])
        tg = np.linspace(-4, 4, 8001)
        g = np.sum(c[None, :] * lams[None, :]
                   * np.exp(-1j * np.outer(tg, wsg)), axis=1)
        return (float(tg[np.argmax(np.abs(g))]),
                float(np.angle(g[np.argmin(np.abs(tg))])))

    peak_corr, phase_corr = packet_peak(th_star, w0)
    th_unc = make_throat(tau, closure_offset(_D_A, _D_B, tau), theta_star)
    peak_unc, _ = packet_peak(th_unc, w0)

    # (b) BELOW THE BARRIER (w0 = 0.5): both closures solved by the
    # throat alone - the root sits ON the interior resonance
    w1 = 0.5

    def carrier_residual(tau_q):
        twq = tau_wigner(w1, tau_q)
        thq = make_throat(tau_q, 0.0)
        return float(np.angle(thq.t_AB(w1) * np.exp(-1j * w1 * twq)))

    taus = np.linspace(4.9, 5.4, 26)
    hs = [carrier_residual(tq) for tq in taus]
    tau_star = None
    for i in range(len(taus) - 1):
        if hs[i] * hs[i + 1] < 0 and abs(hs[i + 1] - hs[i]) < math.pi:
            tau_star = brentq(carrier_residual, taus[i], taus[i + 1],
                              xtol=1e-9)
            break
    th_res = make_throat(tau_star, 0.0)
    T_at_point = float(abs(th_res.t_AB(w1)) ** 2)
    resid = abs(carrier_residual(tau_star))
    tw_res = tau_wigner(w1, tau_star)

    ok = (abs(lam_star - 1) < 1e-4
          and abs(np.angle(lam_star)) < 1e-6
          and abs(peak_corr) < 0.02
          and abs(peak_unc - tw) < 0.02
          and abs(phase_corr) < 0.02
          and tau_star is not None and resid < 1e-8
          and T_at_point > 0.9)
    out = {
        'name': 'T4_transaction_points',
        'description': (
            'Wigner-correct closure lands the packet ON the crossing, '
            'phase-aligned; the above-barrier point needs the mouth '
            'phase; the below-barrier point is solved by the throat '
            'alone and sits ON the interior Fabry-Perot resonance'
        ),
        'above_barrier': {
            'w0': w0, 'tau': tau, 'theta_star': theta_star,
            'delta_star': float(delta_star),
            'Lambda': [float(lam_star.real), float(lam_star.imag)],
            'packet_peak_corrected': peak_corr,
            'packet_peak_uncorrected': peak_unc,
            'wigner': float(tw),
            'phase_at_crossing': phase_corr,
        },
        'below_barrier': {
            'w0': w1, 'tau_star': float(tau_star),
            'carrier_residual': float(resid),
            'T_net_at_point': T_at_point,
            'wigner_at_point': float(tw_res),
        },
        'pass': bool(ok),
    }
    _CACHE['T4'] = out
    return out


# ========================================================================
# T5. Source and mouth state evolution
# ========================================================================


def test_T5_state_evolution() -> dict:
    if 'T5' in _CACHE:
        return _CACHE['T5']
    t4 = test_T4_transaction_points()
    w1 = 0.5
    tau_star = t4['below_barrier']['tau_star']

    # (i) input-output / coupled-mode rates: kappa_tot = (T_A+T_B)/
    # (2 tau) vs the measured linewidth of |t_net|^2 and the resonant
    # Wigner delay 2/kappa_tot
    T1 = float(abs(greybody_exterior(w1)[0]) ** 2)
    kappa_tot = 2 * T1 / (2 * tau_star)
    th_res = make_throat(tau_star, 0.0)
    ws = np.linspace(w1 - 0.08, w1 + 0.08, 41)
    prof = np.array([abs(th_res.t_AB(float(wq))) ** 2 for wq in ws])
    i0 = int(np.argmax(prof))
    above = prof > prof[i0] / 2
    fwhm = float(ws[above][-1] - ws[above][0])
    tw_res = t4['below_barrier']['wigner_at_point']
    cmt_ok = (abs(fwhm - kappa_tot) < 0.25 * kappa_tot
              and abs(tw_res - 2 / kappa_tot) < 0.25 * (2 / kappa_tot))

    # (ii) the source fixed-point iteration x <- 1 + Lambda x converges
    # to G_eff at rate |Lambda| (the winding/generation picture of the
    # self-consistent constraint)
    tau, dlt = 0.8, closure_offset(_D_A, _D_B, 0.8) + 0.13
    th = make_throat(tau, dlt)
    lam = loop_eigenvalue(th, w1, _D_A, _D_B)
    g_exact = effective_green(th, w1, _D_A, _D_B)
    x, errs = 0.0 + 0.0j, []
    for _ in range(40):
        x = 1.0 + lam * x
        errs.append(abs(x - g_exact))
    rate = float(np.mean([errs[i + 1] / errs[i] for i in range(8, 16)]))
    iter_ok = (errs[-1] < 1e-12
               and abs(rate - abs(lam)) < 0.02 * abs(lam))

    # (iii) the mouth-state ledger: cavity amplitude build-up through
    # the transfer matrix iteration saturates at the steady state
    t2, _ = greybody_exterior(w1)
    _, ri2 = greybody_interior(w1)
    M = np.zeros((4, 4), dtype=complex)
    M[1, 0] = t2
    M[2, 1] = ri2 ** 2 * np.exp(2j * w1 * tau)
    M[1, 2] = 1.0
    M[3, 1] = t2
    M[0, 3] = np.exp(1j * w1 * (_D_A + _D_B + tau + dlt))
    v = np.array([1, 0, 0, 0], complex)
    acc = v.copy()
    cavity = []
    for _ in range(30000):
        v = M @ v
        acc += v
        cavity.append(abs(acc[1]) ** 2)
        if len(cavity) > 10 and abs(v).max() < 1e-14:
            break
    x_ss = np.linalg.solve(np.eye(4) - M, np.array([1, 0, 0, 0], complex))
    build_ok = (abs(cavity[-1] - abs(x_ss[1]) ** 2)
                < 1e-3 * abs(x_ss[1]) ** 2
                and cavity[5] < cavity[-1])

    ok = cmt_ok and iter_ok and build_ok
    out = {
        'name': 'T5_state_evolution',
        'description': (
            'input-output rates match the derived linewidth and the '
            'resonant storage time; the source fixed-point iteration '
            'converges at rate |Lambda| (marginal at the completed '
            'transaction); the cavity state builds to steady state'
        ),
        'kappa_tot': float(kappa_tot),
        'linewidth_fwhm': fwhm,
        'wigner_at_resonance': float(tw_res),
        'storage_time_2_over_kappa': float(2 / kappa_tot),
        'iteration_rate': rate,
        'Lambda_magnitude': float(abs(lam)),
        'cavity_steady_state_reached': bool(build_ok),
        'pass': bool(ok),
    }
    _CACHE['T5'] = out
    return out


# ========================================================================
# T6. The effective Green function, then I+-
# ========================================================================


def _fwhm_interp(ws, G2):
    """Half-max width with linear interpolation of the crossings."""
    i0 = int(np.argmax(G2))
    half = G2[i0] / 2
    lo = hi = None
    for i in range(i0, 0, -1):
        if G2[i - 1] < half <= G2[i]:
            lo = ws[i - 1] + (half - G2[i - 1]) / (G2[i] - G2[i - 1]) \
                * (ws[i] - ws[i - 1])
            break
    for i in range(i0, len(ws) - 1):
        if G2[i + 1] < half <= G2[i]:
            hi = ws[i] + (G2[i] - half) / (G2[i] - G2[i + 1]) \
                * (ws[i + 1] - ws[i])
            break
    if lo is None or hi is None:
        return float(ws[-1] - ws[0])
    return float(hi - lo)


def test_T6_geff_vs_orderings() -> dict:
    t4 = test_T4_transaction_points()
    ab = t4['above_barrier']
    w0, tau = ab['w0'], ab['tau']
    th_star = make_throat(tau, ab['delta_star'], ab['theta_star'])

    # (i) the QUARTIC line at the completed transaction (group closure
    # makes the loop phase stationary - an anomalously flat resonance)
    ws = np.linspace(w0 - 0.02, w0 + 0.02, 41)
    lams = np.array([loop_eigenvalue(th_star, float(wq), _D_A, _D_B)
                     for wq in ws])
    G2 = 1.0 / np.abs(1 - lams) ** 2
    i0 = int(np.argmax(G2))
    fwhm_q = _fwhm_interp(ws, G2)
    h = 1e-3
    ph = lambda wq: np.angle(loop_eigenvalue(th_star, wq, _D_A, _D_B))
    thpp = (ph(w0 + h) - 2 * ph(w0) + ph(w0 - h)) / h ** 2
    one_m = 1 - abs(lams[i0])
    fwhm_q_pred = float(2 * math.sqrt(2 * one_m / abs(thpp)))

    # (ii) the LORENTZIAN line when only the carrier closes
    w2 = 0.9
    th0 = make_throat(tau, 0.0)
    theta2 = float(-np.angle(th0.t_AB(w2) * np.exp(1j * w2 * 2.0)))
    th_l = make_throat(tau, closure_offset(_D_A, _D_B, tau) + 2.0, theta2)
    ws2 = np.linspace(w2 - 0.12, w2 + 0.12, 41)
    lams2 = np.array([loop_eigenvalue(th_l, float(wq), _D_A, _D_B)
                      for wq in ws2])
    G22 = 1.0 / np.abs(1 - lams2) ** 2
    i2 = int(np.argmax(G22))
    fwhm_l = _fwhm_interp(ws2, G22)
    ph2 = lambda wq: np.angle(loop_eigenvalue(th_l, wq, _D_A, _D_B))
    thp = (ph2(w2 + h) - ph2(w2 - h)) / (2 * h)
    fwhm_l_pred = float(2 * (1 - abs(lams2[i2])) / abs(thp))

    # (iii) passivity sweep: |Lambda| <= 1 everywhere (Novikov)
    th_p = make_throat(tau, closure_offset(_D_A, _D_B, tau))
    lmax = max(abs(loop_eigenvalue(th_p, float(wq), _D_A, _D_B))
               for wq in np.linspace(0.3, 4.0, 25))

    # (iv) the I+- comparison, only now: at a PARTIAL-confirmation
    # point the winding series resums to the renormalized weight
    # Lambda/(1 - Lambda) exactly (its O(Lambda) truncation is #216's
    # K1 = I+ + Lambda I-); at the COMPLETED point the weight
    # diverges - that divergence IS the transaction pole
    w_p = 0.5
    th_p2 = make_throat(tau, closure_offset(_D_A, _D_B, tau) + 0.13)
    lam_p = loop_eigenvalue(th_p2, w_p, _D_A, _D_B)
    eps = 1e-3
    dgrid = np.linspace(-5, 5, 401)
    Ip = -(1 / (2 * w_p)) / (dgrid - w_p + 1j * eps)
    Im_ = +(1 / (2 * w_p)) / (dgrid + w_p - 1j * eps)
    weight_full = lam_p / (1 - lam_p)
    K_series = Ip + sum(lam_p ** j for j in range(1, 200)) * Im_
    K_closed = Ip + weight_full * Im_
    resum_err = float(np.max(np.abs(K_series - K_closed)))
    lam_c = loop_eigenvalue(th_star, w0, _D_A, _D_B)
    weight_completed = float(abs(lam_c / (1 - lam_c)))

    # (v) the epsilon unification: the confirmation deficit
    # 1 - |Lambda| -> 0 at completion (the #213 coherent limit; the
    # #214 damping in transactional form), and the transaction pole
    # dominates: peak G^2 at completion >> partial by the deficit
    # ratio squared.  The line WIDTH at completion is anomalously
    # broad for its deficit (quartic ~ sqrt(deficit)): group closure
    # flattens the resonance - completed transactions are robust to
    # detuning.
    deficit_completed = float(1 - abs(lam_c))
    deficit_partial = float(1 - abs(lams2[i2]))
    peak_ratio = float((1.0 / deficit_completed ** 2)
                       / np.max(G22))

    ok = (abs(fwhm_q - fwhm_q_pred) < 0.15 * fwhm_q_pred
          and abs(fwhm_l - fwhm_l_pred) < 0.15 * fwhm_l_pred
          and lmax <= 1.0 + 1e-9
          and resum_err < 1e-10
          and weight_completed > 1e4
          and deficit_completed < 1e-4
          and deficit_completed < 1e-3 * deficit_partial
          and peak_ratio > 1e6)
    return {
        'name': 'T6_geff_vs_orderings',
        'description': (
            'G_eff derived first: quartic line at the completed '
            'transaction, Lorentzian at carrier-only closure, widths '
            'predicted by 1-|Lambda|; passivity (Novikov); THEN the '
            'I+- comparison: O(Lambda) = the #216 assembly, resummed '
            'weight Lambda/(1-Lambda), epsilon -> 0 at completion'
        ),
        'quartic_fwhm': fwhm_q,
        'quartic_fwhm_predicted': fwhm_q_pred,
        'lorentzian_fwhm': fwhm_l,
        'lorentzian_fwhm_predicted': fwhm_l_pred,
        'passivity_max_Lambda': float(lmax),
        'resummation_identity_error': resum_err,
        'confirmation_weight_at_completion': weight_completed,
        'deficit_completed': deficit_completed,
        'deficit_partial': deficit_partial,
        'transaction_pole_peak_ratio': peak_ratio,
        'pass': bool(ok),
    }


# ========================================================================
# T7. Honest scope
# ========================================================================


def test_T7_honest_scope() -> dict:
    scope = [
        'The convention correction to #216: the closure comb there '
        'tracked the engine\'s time-phase loop amplitude '
        '(retro_phase_match\'s bookkeeping), which coincides with '
        'value transport on the tower\'s exterior legs but not for '
        'the transit phase; the self-consistent field is governed by '
        'Lambda = t_net e^{+i w D_loop} (ring-validated).  The #216 '
        'deform knob (phi = w delta) and all magnitude/advanced-'
        'projection results are unchanged.',
        'Continuous-w analysis interpolates between the physical '
        'tower modes; the above-barrier transaction point uses the '
        'mouths\' geometric transfer phase as the carrier tuning (a '
        'free parameter of the frozen network), while the '
        'below-barrier point is solved by the throat alone.',
        'The coupled-mode (input-output) rates are valid at moderate '
        'finesse (25% tolerances here); they are a consistency check, '
        'not the derivation.',
        'The elastic selection rule fixes rate_A = rate_B at '
        'traversal epoch; the differential aging that builds Delta_BA '
        'must have happened earlier - the network history is posited '
        '(MTY), not solved.',
        'The winding resummation assumes an unchanged network per '
        'winding (frozen geometry, linear field); back-reaction of '
        'the stored cavity energy on the throat is out of scope.',
        'Classical, zonal scalar, l = 0 greybody throughout; the '
        'engine migration remains staged.',
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
    t2 = test_T2_clock_rates()
    t3 = test_T3_transfer_system()
    t4 = test_T4_transaction_points()
    t5 = test_T5_state_evolution()
    t6 = test_T6_geff_vs_orderings()
    core = all(t['pass'] for t in (t2, t3, t4, t5, t6))
    assessment = (
        'The loop is solved, not sampled: the full two-mouth transfer '
        'system resolves to G_eff = g/(1 - Lambda) with the '
        'value-transport eigenvalue validated against a first-'
        'principles ring spectrum (and the #216 bookkeeping corrected '
        'in the process); clock rates enter exactly (the redshift '
        'inversion fixed, elastic confirmation forcing rate-matched '
        'mouths); Wigner-correct closure lands the packet on the '
        'crossing phase-aligned, and the two closure equations solved '
        'together put the completed transaction ON the throat\'s '
        'interior resonance below the barrier; the source and cavity '
        'states evolve to the steady state at rate |Lambda| with '
        'input-output rates matching the derived linewidth; and the '
        'effective Green function - derived before any I+- comparison '
        '- reproduces #216\'s assembly at first winding order, '
        'renormalizes the confirmation weight to Lambda/(1-Lambda), '
        'shows an anomalously flat (quartic) resonance exactly at '
        'completed transactions, and is passive everywhere: the '
        'Novikov fixed point cannot run away, and the completed '
        'transaction is the #213 coherent epsilon -> 0 limit with the '
        'width given by the confirmation deficit - #214\'s absorber '
        'damping in its transactional form.'
        if core else
        'INCONCLUSIVE - a core check failed; do not quote.'
    )
    return {
        'name': 'T8_assessment',
        'description': 'the standing of the self-consistent solve',
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
        test_T2_clock_rates(),
        test_T3_transfer_system(),
        test_T4_transaction_points(),
        test_T5_state_evolution(),
        test_T6_geff_vs_orderings(),
        test_T7_honest_scope(),
        test_T8_assessment(),
    ]
    t3, t4, t5, t6 = tests[2], tests[3], tests[4], tests[5]
    all_ok = all(t['pass'] for t in tests)
    if all_ok:
        verdict_class = (
            "THE_LOOP_SOLVED_SELF_CONSISTENTLY_G_EFF_DERIVED_"
            "COMPLETED_TRANSACTIONS_SIT_ON_THE_THROAT_RESONANCES_AND_"
            "THE_NOVIKOV_FIXED_POINT_IS_PASSIVE"
        )
        ab, bb = t4['above_barrier'], t4['below_barrier']
        verdict = (
            "DERIVED (the argument is in "
            "docs/self_consistent_network_loop.md).\n\n"
            "THE SYSTEM. Both barriers in one signal-flow system: the "
            "resolvent equals g/(1 - Lambda) exactly "
            f"({t3['matrix_resolvent_error']:.0e}); the interior x "
            "winding double series resums nested "
            f"({t3['nested_resummation_error']:.0e}); and the "
            "value-transport eigenvalue Lambda = t_net e^{i w D_loop} "
            "is validated against a first-principles ring spectrum "
            f"(modes to {t3['ring_prediction_error']:.3f}; the "
            f"opposite convention misses by "
            f"{t3['wrong_convention_miss']:.2f}) - correcting the "
            "#216 closure bookkeeping.\n\n"
            "THE CLOSURES. Clock rates exact (redshift fixed; elastic "
            "confirmation forces rate-matched mouths). Wigner-correct "
            f"closure lands the packet at {ab['packet_peak_corrected']:.3f} "
            f"(uncorrected: {ab['packet_peak_uncorrected']:.3f} = "
            "tau_W), phase-aligned "
            f"({ab['phase_at_crossing']:.3f}). Above the barrier the "
            "mouth phase closes the carrier (Lambda = "
            f"{ab['Lambda'][0]:.6f}); below it BOTH closures are "
            "solved by the throat alone and the completed transaction "
            f"sits ON the interior resonance (|t_net|^2 = "
            f"{bb['T_net_at_point']:.3f} at tau* = "
            f"{bb['tau_star']:.4f}, residual "
            f"{bb['carrier_residual']:.0e}).\n\n"
            "THE STATES. Input-output rates match the derived "
            f"linewidth (kappa = {t5['kappa_tot']:.4f} vs FWHM "
            f"{t5['linewidth_fwhm']:.4f}) and the resonant storage "
            f"time ({t5['wigner_at_resonance']:.1f} vs "
            f"{t5['storage_time_2_over_kappa']:.1f}); the source "
            f"iteration converges at rate |Lambda| "
            f"({t5['iteration_rate']:.4f} vs "
            f"{t5['Lambda_magnitude']:.4f}) - marginal at completion; "
            "the cavity builds to steady state.\n\n"
            "G_EFF, THEN I+-. Quartic line at the completed "
            f"transaction (FWHM {t6['quartic_fwhm']:.4f} vs "
            f"{t6['quartic_fwhm_predicted']:.4f} - group closure "
            "makes the resonance anomalously flat), Lorentzian at "
            f"carrier-only closure ({t6['lorentzian_fwhm']:.4f} vs "
            f"{t6['lorentzian_fwhm_predicted']:.4f}); passivity "
            f"max|Lambda| = {t6['passivity_max_Lambda']:.6f} <= 1 "
            "(Novikov, no runaway); the O(Lambda) truncation IS the "
            "#216 assembly, the resummation renormalizes the "
            "confirmation weight to Lambda/(1-Lambda) "
            f"({t6['resummation_identity_error']:.0e}; at completion "
            f"the weight reaches {t6['confirmation_weight_at_completion']:.1e} "
            "- the divergence IS the transaction pole), and the "
            "confirmation deficit 1-|Lambda| -> 0 at completion "
            f"({t6['deficit_completed']:.1e} vs partial "
            f"{t6['deficit_partial']:.2f}): the completed transaction "
            "is the #213 coherent limit, its deficit the #214 "
            "absorber damping in transactional form."
        )
    else:
        verdict_class = "SELF_CONSISTENT_LOOP_INCONCLUSIVE"
        verdict = (
            "INCONCLUSIVE. A core check failed; re-examine before "
            "quoting the derivation."
        )

    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "identification": (
            "The self-consistent network loop: the full two-mouth "
            "transfer system solved - clock-rate-correct traversal, "
            "Wigner-delay-correct closure, nested winding resummation, "
            "source and mouth state evolution - and the effective "
            "Green function G_eff = g/(1 - Lambda) derived (ring-"
            "validated value-transport eigenvalue) before comparison "
            "with I+-: completed transactions sit on the throat's "
            "interior resonances and the Novikov fixed point is "
            "passive"
        ),
        "executes": (
            "the decisive #216 successor: the loop solved "
            "self-consistently"
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
    out.append("# The self-consistent network loop (PR #217)")
    out.append("")
    out.append(f"**Run:** {s['timestamp_utc']}")
    out.append("")
    out.append(
        "The deliverable is `docs/self_consistent_network_loop.md` - "
        "the full two-mouth transfer system solved and the effective "
        "Green function derived before comparison with I+-. *(QFT on "
        "the fixed classical throat geometry, not quantum gravity.)*"
    )
    out.append("")
    out.append("## Test summary")
    out.append("")
    out.append("| # | Test | Key finding | PASS? |")
    out.append("|---|---|---|---|")
    labels = {
        "T1": "the goal: solve the loop, then compare",
        "T2": "clock-rate-correct traversal (redshift fixed)",
        "T3": "the transfer system; ring-validated eigenvalue",
        "T4": "Wigner-correct closure; the transaction points",
        "T5": "source and mouth state evolution",
        "T6": "G_eff first, then I+-: lines, passivity, epsilon",
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
    out = here / "runs" / f"{ts}_self_consistent_network_loop_probe"
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
