"""
The wormhole-network confirmation: the advanced half of the Wheeler
transaction as a globally causal, everywhere-future-directed network
traversal - through a TWO-PORT throat with repeated-loop physics
(PR #216).

> Framing: QFT on the *fixed classical* throat geometry (geometry -> fields),
> not quantum gravity.

THE PROBLEM WITH "ADVANCED"
---------------------------
Wheeler-Feynman transactions are retarded + advanced, and #213 showed
the frozen bulk forces the combination.  But waves propagate forward in
time.  The engine's ``advanced_confirm_amplitude`` assigns the
confirmation's conjugated phase BY FIAT.  This probe supplies the
mechanism: a wormhole NETWORK can carry a retarded wave into the global
past with every segment locally future-directed -

    emit one retarded C-wave                        (t = 0, psi = 0)
    -> propagate forward to the future antipode     (t = pi)
    -> transmit through mouth A's interface         (#215 greybody)
    -> traverse the throat forward; LOOP k times    ((2k+1) tau_th > 0
       between the two interior barrier faces        local each)
    -> exit through mouth B's OWN interface         (second barrier)
    -> emerge through the clock-offset mouth        (t < 0: differential
                                                     aging, Morris-
                                                     Thorne-Yurtsever)
    -> propagate forward again on S^3               (back to psi = 0)
    -> intersect the original particle crossing     (t_return = t_emit)

THE TWO-PORT THROAT.  Each mouth carries its own barrier: entry
transmission t_A, interior faces r_inA/r_inB, exit transmission t_B.
The composite transmission is DERIVED by summing the interior loops,

    t_net(w) = t_A t_B / (1 - r_inA r_inB e^{2 i w tau_th}) ,

validated against a DIRECT solve of the glued two-barrier potential
(3e-4 in |T|^2 across a sweep through resonances).  New physics:
RESONANT CONFIRMATION - identical mouths transmit perfectly on the
interior resonance comb even deep below the barrier (T_net = 1 at
single-port T = 0.33), puncturing the #215 IR transparency; off
resonance soft modes confirm at T^2/(1+R)^2.

THE CENTRAL IDENTITY.  Clock continuity forces the throat's
global-frame transfer factor to be t_net(w) e^{i w Delta_BA}: the
projected kernel between absorption and return is

    K = Lambda(w) * e^{-i w D'} ,   D' = tau_th + Delta_BA + d_B < 0 ,
    Lambda(w) = t_net(w) e^{i w Delta_BA} * decorations ,

i.e. the retarded phase rule analytically continued to a NEGATIVE
exterior interval - THE ADVANCED KERNEL - with greybody-composite
weight |Lambda| = |t_net|.  #213's coherence requirement (relative
phase phi = 0) becomes the network's PHASE-CLOSURE condition, satisfied
on a discrete comb of frequencies: the transaction selects its modes.
The two Compton completion families are two globally causal path
classes coherent exactly at closure.

Tests:
  T1. Goal.
  T2. The traversal ledger + the echo train: local future-
      directedness of every leg INCLUDING the repeated loops, global
      reach into the past, time closure at the crossing, the derived
      transfer factor.
  T3. The two-port mouth S-matrix: exterior AND interior greybody
      (match to #215, flux closure both sides, the free-interface
      null test, reciprocity t_in = t, the unitarity relation
      r_in = -conj(r) t / conj(t), the Wigner delay).
  T4. The composite throat: transparent reduction, loop-series
      convergence, the DIRECT two-barrier validation, resonant
      tunneling below the barrier, the Airy average identity.
  T5. The advanced projection: K = Lambda e^{-i w D'} with D' < 0 and
      |Lambda| = |t_net|; the phase-closure comb; the engine weld; the
      live packet (emerges in the past, intersects the crossing within
      the composite Wigner delay; below-barrier packets confirm
      strongly ON the interior resonance and weakly off - resonant
      confirmation).
  T6. The pole structure: I+ + Lambda I- reconstructs the covariant
      pole at closure; network detuning reproduces the #213 deform
      diagnostic; partial confirmation splits by |t_net|.
  T7. Honest scope.
  T8. Assessment.

Verdict:
  THE_ADVANCED_HALF_IS_A_FUTURE_DIRECTED_NETWORK_TRAVERSAL_THE_CLOCK_
  OFFSET_MOUTH_IS_THE_ANALYTIC_CONTINUATION_AND_PHASE_CLOSURE_IS_THE_
  COHERENCE_CONDITION_NOW_THROUGH_A_TWO_PORT_THROAT
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

from geometrodynamics.transaction import (
    MouthPort,
    NetworkMouth,
    NetworkThroat,
    closure_offset,
    emergence_train,
    network_confirmation,
    projected_kernel,
    retro_phase_match,
    transparent_port,
    traverse_throat,
)

_CACHE: dict = {}
_RH = 1.0

# ========================================================================
# SECTION A - the two-sided greybody of a Tangherlini mouth (PR #215
# solver, plane-wave referenced so a switched-off barrier gives t = 1)
# ========================================================================


def _v_of_r(r, l: int = 0, strength: float = 1.0):
    f = 1.0 - (_RH / r) ** 2
    return strength * f * ((l * (l + 2) + 0.75) / r ** 2
                           + 2.25 * _RH ** 2 / r ** 4)


def _x_of_r(r):
    return r + (_RH / 2) * np.log((r - _RH) / (r + _RH))


def _rhs_factory(w: float, strength: float = 1.0):
    def rhs(x, y):
        pr, pi_, qr, qi, r = y
        vv = _v_of_r(r, 0, strength)
        f = 1 - (_RH / r) ** 2
        return [qr, qi, (vv - w ** 2) * pr, (vv - w ** 2) * pi_, f]
    return rhs


def greybody_exterior(w: float, strength: float = 1.0,
                      rtol: float = 3e-10) -> tuple:
    """Exterior incidence: (t, r_out)."""
    key = ('ge', w, strength)
    if key in _CACHE:
        return _CACHE[key]
    r0 = _RH * (1 + 1e-7)
    x0 = float(_x_of_r(r0))
    r_far = max(60.0 / w, 50.0 * w, 30.0)
    x_far = float(_x_of_r(r_far))
    y0 = [math.cos(-w * x0), math.sin(-w * x0),
          w * math.sin(-w * x0), -w * math.cos(-w * x0), r0]
    sol = solve_ivp(_rhs_factory(w, strength), (x0, x_far), y0,
                    rtol=rtol, atol=1e-13, dense_output=True,
                    method="DOP853")
    xs = x_far - np.linspace(0.0, 3.7, 8)
    rows, vals = [], []
    for xm in xs:
        pr, pi_, *_ = sol.sol(xm)
        rows.append([np.exp(-1j * w * xm), np.exp(1j * w * xm)])
        vals.append(pr + 1j * pi_)
    (alpha, beta), *_ = np.linalg.lstsq(np.array(rows), np.array(vals),
                                        rcond=None)
    out = (complex(1.0 / alpha), complex(beta / alpha))
    _CACHE[key] = out
    return out


def greybody_interior(w: float, rtol: float = 3e-10) -> tuple:
    """Interior (horizon-side) incidence: (t_in, r_in) - the face the
    repeated loops bounce off.  Outgoing-only at large x, integrated
    inward, decomposed near the horizon where V is exponentially
    small."""
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
    sol = solve_ivp(_rhs_factory(w), (x_far, x0), y0, rtol=rtol,
                    atol=1e-13, dense_output=True, method="DOP853")
    xs = x0 + np.linspace(0.0, 2.0, 8)
    rows, vals = [], []
    for xm in xs:
        pr, pi_, *_ = sol.sol(xm)
        rows.append([np.exp(1j * w * xm), np.exp(-1j * w * xm)])
        vals.append(pr + 1j * pi_)
    (alpha, beta), *_ = np.linalg.lstsq(np.array(rows), np.array(vals),
                                        rcond=None)
    out = (complex(1.0 / alpha), complex(beta / alpha))
    _CACHE[key] = out
    return out


def tangherlini_port() -> MouthPort:
    """The full one-sided S-matrix of a Tangherlini mouth."""
    return MouthPort(
        t=lambda w: greybody_exterior(w)[0],
        r_out=lambda w: greybody_exterior(w)[1],
        r_in=lambda w: greybody_interior(w)[1],
    )


# V as a function of tortoise x, for the glued two-barrier direct solve
_rg = 1 + np.logspace(-9, math.log10(299.0), 4000)
_xg = _x_of_r(_rg)
_vg = _v_of_r(_rg)
_v_interp = interp1d(_xg, _vg, kind="cubic", bounds_error=False,
                     fill_value=(0.0, 0.0))


def _v_of_x(x):
    x = float(x)
    if x < _xg[0]:
        return 0.0
    if x > _xg[-1]:
        return 0.75 / x ** 2
    return float(_v_interp(x))


_X_PK = 0.1978          # barrier-peak tortoise position (#215)


def two_barrier_direct(w: float, c: float = 4.0,
                       rtol: float = 3e-10) -> complex:
    """DIRECT transmission through the glued two-mouth potential
    V(y) = V_x(x_pk - (y + c)) + V_x(x_pk + (y - c)) - the ground
    truth the Fabry-Perot composition must reproduce."""
    y_far = c + max(60.0 / w, 50.0 * w, 30.0)

    def rhs(y, s):
        pr, pi_, qr, qi = s
        vv = _v_of_x(_X_PK - (y + c)) + _v_of_x(_X_PK + (y - c))
        return [qr, qi, (vv - w ** 2) * pr, (vv - w ** 2) * pi_]

    y0 = -y_far
    p0 = np.exp(-1j * w * y0)
    q0 = -1j * w * p0
    sol = solve_ivp(rhs, (y0, y_far), [p0.real, p0.imag, q0.real, q0.imag],
                    rtol=rtol, atol=1e-13, dense_output=True,
                    method="DOP853")
    ys = y_far - np.linspace(0.0, 3.7, 8)
    rows, vals = [], []
    for ym in ys:
        pr, pi_, *_ = sol.sol(ym)
        rows.append([np.exp(-1j * w * ym), np.exp(1j * w * ym)])
        vals.append(pr + 1j * pi_)
    (alpha, beta), *_ = np.linalg.lstsq(np.array(rows), np.array(vals),
                                        rcond=None)
    return complex(1.0 / alpha)


# ── the standard network of this probe ──────────────────────────────────

_D_A = math.pi          # source -> future antipodal mouth A
_D_B = math.pi          # past mouth B -> source
_TAU = 0.8              # interior transit between the port planes


def standard_throat(delta: Optional[float] = None,
                    tau_th: float = _TAU) -> NetworkThroat:
    if delta is None:
        delta = closure_offset(_D_A, _D_B, tau_th)
    A = NetworkMouth("A", psi=math.pi, link_id="L1", clock_offset=0.0)
    B = NetworkMouth("B", psi=math.pi, link_id="L1", clock_offset=delta)
    return NetworkThroat(A, B, tau_th=tau_th,
                         port_A=tangherlini_port(),
                         port_B=tangherlini_port())


def resonant_tau(w: float, k: int = 2) -> float:
    """Interior transit for which the throat's Fabry-Perot resonance
    sits exactly at frequency w: arg(r_in^2) + 2 w tau = 2 pi k."""
    _, r_in = greybody_interior(w)
    return (2 * math.pi * k - 2 * np.angle(r_in)) / (2 * w)


# ========================================================================
# T1. Goal
# ========================================================================


def test_T1_goal() -> dict:
    return {
        'name': 'T1_goal',
        'description': (
            'Replace advanced_confirm_amplitude with an explicit '
            'network traversal through a TWO-PORT throat: greybody '
            'entry at mouth A, repeated interior loops between the two '
            'barrier faces, exit through mouth B\'s own interface, '
            'emergence through the clock-offset mouth in the global '
            'past, and return to the original crossing - showing the '
            'projected response is advanced, every segment (including '
            'every loop) is locally future-directed, the closures '
            'hold, and the effective kernel has the #213 pole '
            'structure. The advanced half of a Wheeler transaction '
            'gets a classical geometric mechanism with the throat\'s '
            'full interior dynamics.'
        ),
        'pass': True,
    }


# ========================================================================
# T2. The traversal ledger + the echo train
# ========================================================================


def test_T2_traversal_ledger() -> dict:
    th = standard_throat()
    w = 3.0
    tr = network_confirmation(th, w, t_emit=0.0, d_A=_D_A, d_B=_D_B)

    ledger = [{
        'leg': leg.name,
        't_start_global': float(leg.t_start),
        't_end_global': float(leg.t_end),
        'local_duration': float(leg.local_duration),
    } for leg in tr.legs]

    # the echo train at a below-barrier frequency (visible loops):
    # emergences at t_absorb + (2k+1) tau + Delta, geometrically damped,
    # every echo future-directed in the throat clock
    w_lo = 0.5
    train = emergence_train(th, w_lo, t_entry=math.pi, kmax=3)
    _, r_in = greybody_interior(w_lo)
    damp_pred = abs(r_in) ** 2
    train_ok = True
    echo_ledger = []
    for k, leg in enumerate(train):
        echo_ledger.append({
            'k': k,
            'local_duration': float(leg.local_duration),
            't_emerge_global': float(leg.t_end),
            'amplitude': float(abs(leg.factor)),
        })
        train_ok = train_ok and leg.local_duration > 0
        train_ok = train_ok and abs(
            leg.local_duration - (2 * k + 1) * _TAU) < 1e-12
        if k:
            ratio = abs(leg.factor) / abs(train[k - 1].factor)
            train_ok = train_ok and abs(ratio - damp_pred) < 1e-9

    # the derived transfer factor at every frequency
    transfer_err = 0.0
    for wt in (0.7, 1.9, 3.3):
        leg = traverse_throat(th, wt, t_entry=math.pi)
        gframe = leg.factor * np.exp(1j * wt * (leg.t_end - leg.t_start))
        expected = th.t_AB(wt) * np.exp(1j * wt * th.delta_BA)
        transfer_err = max(transfer_err, abs(gframe - expected))

    ok = (tr.locally_forward
          and tr.globally_advanced
          and tr.t_emerge < tr.t_emit < tr.t_absorb
          and abs(tr.t_return - tr.t_emit) < 1e-12
          and transfer_err < 1e-12
          and train_ok)
    return {
        'name': 'T2_traversal_ledger',
        'description': (
            'every leg future-directed in its own clock - including '
            'every interior loop of the echo train; the global past '
            'reached only through the frozen mouth offset; time '
            'closure at the crossing; U_BA = e^{i w Delta_BA} derived'
        ),
        'ledger': ledger,
        'echo_train': echo_ledger,
        'echo_damping_predicted': float(damp_pred),
        't_emit': float(tr.t_emit),
        't_absorb': float(tr.t_absorb),
        't_emerge': float(tr.t_emerge),
        't_return': float(tr.t_return),
        'delta_BA': float(th.delta_BA),
        'locally_forward': bool(tr.locally_forward),
        'globally_advanced': bool(tr.globally_advanced),
        'transfer_factor_identity_error': float(transfer_err),
        'pass': bool(ok),
    }


# ========================================================================
# T3. The two-port mouth S-matrix
# ========================================================================


def test_T3_mouth_smatrix() -> dict:
    if 'T3' in _CACHE:
        return _CACHE['T3']
    # (i) exterior side against the #215 committed values
    ref = {0.1: 1.663731e-03, 0.5: 3.331217e-01,
           1.0: 9.629316e-01, 2.0: 0.9998469}
    rel_err, flux_out = 0.0, 0.0
    for w, Tref in ref.items():
        t, r = greybody_exterior(w)
        rel_err = max(rel_err, abs(abs(t) ** 2 - Tref) / Tref)
        flux_out = max(flux_out, abs(abs(t) ** 2 + abs(r) ** 2 - 1))

    # (ii) interior side: flux, reciprocity, the unitarity relation
    recip_err, flux_in, unit_rel = 0.0, 0.0, 0.0
    for w in (0.5, 1.0, 2.0):
        t, r = greybody_exterior(w)
        t_in, r_in = greybody_interior(w)
        recip_err = max(recip_err, abs(t_in - t) / abs(t))
        flux_in = max(flux_in, abs(abs(t_in) ** 2 + abs(r_in) ** 2 - 1))
        unit_rel = max(unit_rel,
                       abs(r_in - (-np.conj(r) * t / np.conj(t))))

    # (iii) the free-interface null test
    t_null, _ = greybody_exterior(1.3, strength=0.0)
    null_err = abs(t_null - 1.0)

    # (iv) the Wigner delay (exterior amplitude)
    def wig(w, h=1e-3):
        return float((np.angle(greybody_exterior(w + h)[0])
                      - np.angle(greybody_exterior(w - h)[0])) / (2 * h))
    tw = {w: wig(w) for w in (0.5, 1.0, 3.0)}

    ok = (rel_err < 1e-3 and flux_out < 3e-4 and null_err < 1e-8
          and recip_err < 1e-3 and flux_in < 1e-4 and unit_rel < 1e-3
          and tw[0.5] > tw[1.0] > tw[3.0] > 0)
    out = {
        'name': 'T3_mouth_smatrix',
        'description': (
            'the full one-sided S-matrix of a mouth: exterior t/r_out '
            '(#215 match), interior r_in (the loop face) with '
            'reciprocity t_in = t and the unitarity relation '
            'r_in = -conj(r) t/conj(t); V = 0 gives t = 1 exactly'
        ),
        'match_to_215_relative': float(rel_err),
        'flux_closure_exterior': float(flux_out),
        'flux_closure_interior': float(flux_in),
        'reciprocity_error': float(recip_err),
        'unitarity_relation_error': float(unit_rel),
        'free_interface_null_error': float(null_err),
        'wigner_delay': {str(k): float(v) for k, v in tw.items()},
        'pass': bool(ok),
    }
    _CACHE['T3'] = out
    return out


# ========================================================================
# T4. The composite throat: repeated loops, validated directly
# ========================================================================


def test_T4_composite_throat() -> dict:
    if 'T4' in _CACHE:
        return _CACHE['T4']
    th = standard_throat()

    # (i) transparent-port reduction (the v1 single-interface limit)
    A = NetworkMouth("A", psi=math.pi, link_id="L")
    B = NetworkMouth("B", psi=math.pi, link_id="L", clock_offset=-7.0)
    th_open = NetworkThroat(A, B, tau_th=_TAU,
                            port_A=transparent_port(),
                            port_B=transparent_port())
    open_err = max(abs(th_open.t_AB(w) - 1.0) for w in (0.5, 2.0))

    # (ii) loop-series convergence to the closed form (real ports)
    conv_err = 0.0
    for w in (0.5, 1.0):
        partial = sum(th.loop_expansion(w, 80))
        conv_err = max(conv_err, abs(partial - th.t_AB(w)))

    # (iii) THE DIRECT VALIDATION: the glued two-barrier potential vs
    # the Fabry-Perot composition, |T|^2 across a sweep through
    # resonances (interior separation D = 2(c - x_pk))
    c = 4.0
    D = 2 * (c - _X_PK)
    th_D = standard_throat(tau_th=D)
    direct_errs = []
    for w in np.linspace(0.42, 0.68, 10):
        comp = th_D.t_AB(float(w))
        direct = two_barrier_direct(float(w), c=c)
        direct_errs.append(abs(abs(direct) ** 2 - abs(comp) ** 2))
    direct_err = float(max(direct_errs))

    # (iv) resonant tunneling: identical ports transmit PERFECTLY on
    # the interior resonance even deep below the barrier
    w_r = 0.5
    T_single = abs(greybody_exterior(w_r)[0]) ** 2
    th_res = standard_throat(tau_th=resonant_tau(w_r))
    T_on = abs(th_res.t_AB(w_r)) ** 2
    R1 = 1 - T_single
    T_off_pred = T_single ** 2 / (1 + R1) ** 2
    # off resonance: shift tau by a quarter FP period
    th_off = standard_throat(tau_th=resonant_tau(w_r)
                             + math.pi / (2 * w_r))
    T_off = abs(th_off.t_AB(w_r)) ** 2

    # (v) the Airy average identity: mean T over one FP period of tau
    # equals the incoherent echo sum T^2/(1 - R^2)
    taus = resonant_tau(w_r) + np.linspace(0, math.pi / w_r, 2001)
    _, r_in = greybody_interior(w_r)
    t1 = greybody_exterior(w_r)[0]
    Tnet = np.abs(t1 * t1 / (1 - r_in ** 2
                             * np.exp(2j * w_r * taus))) ** 2
    airy = float(np.mean(Tnet[:-1]))
    airy_pred = T_single ** 2 / (1 - R1 ** 2)

    ok = (open_err < 1e-14 and conv_err < 1e-10
          and direct_err < 1e-3
          and abs(T_on - 1.0) < 2e-3
          and abs(T_off - T_off_pred) < 2e-3
          and T_on / T_off > 20
          and abs(airy - airy_pred) < 2e-3)
    out = {
        'name': 'T4_composite_throat',
        'description': (
            'the two-port composition validated against a direct '
            'glued-barrier solve; resonant tunneling: perfect '
            'transmission on the interior comb below the barrier; '
            'the Airy average = the incoherent echo sum'
        ),
        'transparent_reduction_error': float(open_err),
        'loop_series_convergence': float(conv_err),
        'direct_two_barrier_validation': direct_err,
        'single_port_T': float(T_single),
        'resonant_T_net': float(T_on),
        'off_resonance_T_net': float(T_off),
        'off_resonance_predicted': float(T_off_pred),
        'resonant_enhancement': float(T_on / T_off),
        'airy_average': airy,
        'airy_predicted': float(airy_pred),
        'pass': bool(ok),
    }
    _CACHE['T4'] = out
    return out


# ========================================================================
# T5. The advanced projection, phase closure, and the live packet
# ========================================================================


def test_T5_advanced_projection() -> dict:
    if 'T5' in _CACHE:
        return _CACHE['T5']
    th = standard_throat()

    # (i) the projected kernel: pastward interval, composite weight
    interval_ok, weight_err = True, 0.0
    for w in (0.8, 2.0, 3.2):
        pk = projected_kernel(th, w, _D_B)
        interval_ok = interval_ok and pk['interval'] < 0
        weight_err = max(weight_err,
                         abs(abs(pk['confirmation_weight'])
                             - pk['greybody_magnitude']))

    # (ii) phase closure on a discrete comb: arg Lambda(w) = 0 (mod 2pi)
    ws = np.linspace(2.2, 4.2, 81)
    argl = np.array([np.angle(th.t_AB(float(w))
                              * np.exp(1j * w * th.delta_BA))
                     for w in ws])
    roots = []
    for i in range(len(ws) - 1):
        d = argl[i + 1] - argl[i]
        if argl[i] * argl[i + 1] < 0 and abs(d) < math.pi:
            roots.append(float(ws[i] - argl[i]
                               * (ws[i + 1] - ws[i]) / d))

    def wig_net(w, h=1e-3):
        return float((np.angle(th.t_AB(w + h))
                      - np.angle(th.t_AB(w - h))) / (2 * h))
    spacing_errs = []
    for i, s in enumerate(np.diff(roots)):
        w_mid = 0.5 * (roots[i] + roots[i + 1])
        pred = 2 * math.pi / abs(th.delta_BA + wig_net(w_mid))
        spacing_errs.append(abs(s / pred - 1))
    spacing_err = float(max(spacing_errs))

    # (iii) the engine weld
    w_close = roots[0]
    tr = network_confirmation(th, w_close, 0.0, _D_A, _D_B)
    _, weight_close = retro_phase_match(1.0 + 0.0j, tr.amp)
    w_off = w_close + 0.25 * (roots[1] - roots[0])
    tr_off = network_confirmation(th, w_off, 0.0, _D_A, _D_B)
    _, weight_off = retro_phase_match(1.0 + 0.0j, tr_off.amp)

    # (iv) the live packet, above the barrier (loops negligible there)
    hi = _packet_run(th, 3.0, 0.3, npts=41)
    # (v) resonant confirmation below the barrier: a narrow packet
    # centered ON an interior resonance confirms strongly - stored in
    # the cavity and released over the storage time (the composite
    # Wigner delay) - while off resonance it barely confirms.  The
    # k = 1 resonance keeps the line broad; the time grid must span
    # the ringdown.
    w_r = 0.5
    th_res = standard_throat(tau_th=resonant_tau(w_r, k=1))
    on = _packet_run(th_res, w_r, 0.008, npts=33, t_span=250.0)
    th_off = standard_throat(tau_th=resonant_tau(w_r, k=1)
                             + math.pi / (2 * w_r))
    off = _packet_run(th_off, w_r, 0.008, npts=33, t_span=250.0)

    ok = (interval_ok and weight_err < 1e-12
          and len(roots) >= 2 and spacing_err < 0.1
          and weight_close > 0.99 and weight_off < 0.5
          and hi['emergence_peak'] < 0.0
          and abs(hi['return_peak'] - hi['wigner_net']) < 0.05
          and abs(hi['energy_fraction'] - hi['mean_T']) < 1e-3
          and hi['energy_fraction'] > 0.99
          and on['energy_fraction'] > 0.7
          and abs(on['energy_fraction'] - on['mean_T']) < 0.05
          and off['energy_fraction'] < 0.08
          and on['energy_fraction'] / off['energy_fraction'] > 10
          and on['return_peak'] > 3.0)      # stored, delayed release
    out = {
        'name': 'T5_advanced_projection',
        'description': (
            'K = Lambda e^{-iwD\'} with D\' < 0, |Lambda| = |t_net|: '
            'the advanced kernel; the closure comb; the engine weld; '
            'live packets - and RESONANT CONFIRMATION: below-barrier '
            'packets confirm on the interior comb, delayed by the '
            'cavity storage time'
        ),
        'pastward_interval': bool(interval_ok),
        'weight_equals_composite_error': float(weight_err),
        'closure_comb_roots': [float(r) for r in roots],
        'comb_spacing_max_error': spacing_err,
        'engine_weight_at_closure': float(weight_close),
        'engine_weight_at_quarter_detune': float(weight_off),
        'packet_above_barrier': hi,
        'packet_on_resonance': on,
        'packet_off_resonance': off,
        'resonant_confirmation_contrast': float(
            on['energy_fraction'] / off['energy_fraction']),
        'pass': bool(ok),
    }
    _CACHE['T5'] = out
    return out


def _packet_run(th: NetworkThroat, w0: float, sig: float,
                npts: int = 41, t_span: float = 40.0) -> dict:
    D = _D_A + th.tau_th + th.delta_BA + _D_B     # = 0 at time closure
    t_e = _D_A + th.tau_th + th.delta_BA
    ws = np.linspace(w0 - 4 * sig, w0 + 4 * sig, npts)
    c = np.exp(-0.5 * ((ws - w0) / sig) ** 2)
    t_amp = np.array([th.t_AB(float(w)) for w in ws])

    tgrid = np.linspace(-t_span, t_span, int(200 * t_span) + 1)
    f0 = np.sum(c[None, :] * np.exp(-1j * np.outer(tgrid, ws)), axis=1)
    g = np.sum(c[None, :] * t_amp[None, :]
               * np.exp(-1j * np.outer(tgrid - D, ws)), axis=1)
    sgrid = t_e + tgrid
    h = np.sum(c[None, :] * t_amp[None, :]
               * np.exp(-1j * np.outer(sgrid - t_e, ws)), axis=1)

    def wig_net(w, hh=1e-3):
        return float((np.angle(th.t_AB(w + hh))
                      - np.angle(th.t_AB(w - hh))) / (2 * hh))

    e_frac = float(np.trapezoid(np.abs(g) ** 2, tgrid)
                   / np.trapezoid(np.abs(f0) ** 2, tgrid))
    mean_T = float(np.sum(c ** 2 * np.abs(t_amp) ** 2) / np.sum(c ** 2))
    return {
        'return_peak': float(tgrid[np.argmax(np.abs(g))]),
        'emergence_peak': float(sgrid[np.argmax(np.abs(h))]),
        't_emergence_geometric': float(t_e),
        'wigner_net': wig_net(w0),
        'energy_fraction': e_frac,
        'mean_T': mean_T,
    }


# ========================================================================
# T6. The pole structure
# ========================================================================


def test_T6_pole_structure() -> dict:
    th = standard_throat()
    t5 = test_T5_advanced_projection()
    w_close = t5['closure_comb_roots'][1]     # a closing mode, T ~ 1
    eps = 1e-3

    def ordering_plus(dl, w):
        return -(1 / (2 * w)) / (dl - w + 1j * eps)

    def ordering_minus(dl, w):
        return +(1 / (2 * w)) / (dl + w - 1j * eps)

    dgrid = np.linspace(-5, 5, 801)
    mask = np.abs(np.abs(dgrid) - w_close) > 0.3

    lam = th.t_AB(w_close) * np.exp(1j * w_close * th.delta_BA)
    K = ordering_plus(dgrid, w_close) + lam * ordering_minus(dgrid, w_close)
    pole_exact = -(w_close - 1j * eps) / (
        w_close * (dgrid ** 2 - (w_close - 1j * eps) ** 2))
    dev_closed = float(np.max(np.abs((K - pole_exact) / pole_exact)[mask]))
    lam_phase = float(abs(np.angle(lam)))
    lam_mag = float(abs(lam))

    deform = {}
    for delta in (0.1, 0.5):
        th_d = standard_throat(delta=th.delta_BA + delta)
        lam_d = th_d.t_AB(w_close) * np.exp(1j * w_close * th_d.delta_BA)
        phi_geo = float(np.angle(lam_d / lam))
        phi_pred = ((w_close * delta + math.pi) % (2 * math.pi)) - math.pi
        K_d = (ordering_plus(dgrid, w_close)
               + lam_d * ordering_minus(dgrid, w_close))
        dev_d = float(np.max(np.abs((K_d - pole_exact) / pole_exact)[mask]))
        K_213 = (ordering_plus(dgrid, w_close)
                 + np.exp(1j * phi_geo)
                 * ordering_minus(dgrid, w_close))
        dev_213 = float(np.max(np.abs((K_213 - pole_exact)
                                      / pole_exact)[mask]))
        deform[delta] = {'phi_geometric': phi_geo,
                         'phi_predicted_w_delta': phi_pred,
                         'pole_deviation': dev_d,
                         'pole_deviation_213_diag': dev_213}

    # partial confirmation off the interior comb: |t_net| < 1 splits
    # into coherent pole share and ordering-asymmetric deficit
    w_low = 0.5
    s = abs(standard_throat().t_AB(w_low))
    Kp = ordering_plus(dgrid, w_low) + s * ordering_minus(dgrid, w_low)
    share = (1 + s) / 2 * (ordering_plus(dgrid, w_low)
                           + ordering_minus(dgrid, w_low))
    deficit = (1 - s) / 2 * (ordering_plus(dgrid, w_low)
                             - ordering_minus(dgrid, w_low))
    split_err = float(np.max(np.abs(Kp - share - deficit)))

    ok = (dev_closed < 1e-4
          and lam_phase < 0.03 and lam_mag > 0.995
          and all(abs(d['phi_geometric'] - d['phi_predicted_w_delta'])
                  < 1e-9 for d in deform.values())
          and all(abs(d['pole_deviation'] - d['pole_deviation_213_diag'])
                  < 0.05 * d['pole_deviation_213_diag'] + 1e-6
                  for d in deform.values())
          and split_err < 1e-12)
    return {
        'name': 'T6_pole_structure',
        'description': (
            'the two completion families as path classes: at network '
            'closure I+ + Lambda I- is the covariant pole; network '
            'detuning IS the #213 deform test (phi = w delta); partial '
            'confirmation splits coherent share/deficit by |t_net|'
        ),
        'closing_mode': float(w_close),
        'Lambda_magnitude': lam_mag,
        'Lambda_phase': lam_phase,
        'pole_deviation_at_closure': dev_closed,
        'deform_knob': {str(k): v for k, v in deform.items()},
        'partial_confirmation_split_error': split_err,
        'pass': bool(ok),
    }


# ========================================================================
# T7. Honest scope
# ========================================================================


def test_T7_honest_scope() -> dict:
    scope = [
        'The network breaks global hyperbolicity by construction (a '
        'time-machine geometry); the closure condition is exactly the '
        'Novikov self-consistency fixed point - transactions are its '
        'solutions.  No chronology-protection dynamics is addressed: '
        'the network is frozen background, consistent with the '
        'program framing.',
        'The network geometry is posited, not solved: the mouth pair '
        'from the #58/#200 nucleation channel, the offset Delta_BA '
        'from differential aging (Morris-Thorne-Yurtsever mechanism) '
        '- cited as mechanism, not derived from the 5D field '
        'equations.',
        'The interior transit tau_th and the mapping of the two '
        'mouths\' near-horizon regions onto a single 1D channel are '
        'model-level (the glued-potential direct solve validates the '
        'composition on exactly that model); the true 5D interior '
        'is #168/#200 territory.',
        'Single-mode scalar, zonal (1D transit) reduction throughout; '
        'the greybody is the l = 0 channel.',
        'The elastic case (equal mouth clock rates) is analyzed; '
        'unequal rates red/blueshift the traversal and are only '
        'unit-tested.',
        'Time closure is tuned to the PRIMARY (k = 0) emergence; the '
        'k-th echo returns 2k tau_th after the crossing (a network '
        'can instead close on any single echo).  The composite t_net '
        'is the coherent monochromatic sum; packet-level echo '
        'resolution would need bandwidth > 1/tau_th.',
        'The engine\'s advanced_confirm_amplitude is retained; '
        'network_confirmation is the mechanism-level replacement and '
        'the T5 weld shows the engine\'s phase closure accepts it - '
        'full engine migration is staged, not performed.',
        'Which tower modes land on the closure comb (and on the '
        'interior resonance comb) depends on network parameters - '
        'selectivity is demonstrated, not matched to a spectrum.',
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
    t2 = test_T2_traversal_ledger()
    t3 = test_T3_mouth_smatrix()
    t4 = test_T4_composite_throat()
    t5 = test_T5_advanced_projection()
    t6 = test_T6_pole_structure()
    core = all(t['pass'] for t in (t2, t3, t4, t5, t6))
    assessment = (
        'The advanced half of the Wheeler transaction now has a '
        'classical geometric mechanism with the throat\'s full '
        'interior dynamics: a retarded wave, absorbed at the future '
        'antipodal mouth, transmitted through TWO interfaces with '
        'repeated interior loops between them (every loop locally '
        'future-directed), and re-emitted through the clock-offset '
        'mouth in the global past, projects onto the exterior as '
        'EXACTLY the advanced kernel with the composite two-port '
        'weight |t_net| - validated against a direct glued-barrier '
        'solve.  The repeated-loop physics adds resonant '
        'confirmation: soft modes, IR-transparent at a single '
        'interface, confirm perfectly on the throat\'s interior '
        'resonance comb.  #213\'s coherent relative phase is '
        're-derived as the network\'s phase-closure condition, its '
        'deform test becomes a geometric detuning knob, and the two '
        'Compton completion families become two globally causal '
        'wormhole-network path classes.'
        if core else
        'INCONCLUSIVE - a core check failed; do not quote.'
    )
    return {
        'name': 'T8_assessment',
        'description': 'the standing of the network-confirmation model',
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
        test_T2_traversal_ledger(),
        test_T3_mouth_smatrix(),
        test_T4_composite_throat(),
        test_T5_advanced_projection(),
        test_T6_pole_structure(),
        test_T7_honest_scope(),
        test_T8_assessment(),
    ]
    t2, t4, t5, t6 = tests[1], tests[3], tests[4], tests[5]
    all_ok = all(t['pass'] for t in tests)
    if all_ok:
        verdict_class = (
            "THE_ADVANCED_HALF_IS_A_FUTURE_DIRECTED_NETWORK_TRAVERSAL"
            "_THE_CLOCK_OFFSET_MOUTH_IS_THE_ANALYTIC_CONTINUATION_AND"
            "_PHASE_CLOSURE_IS_THE_COHERENCE_CONDITION_NOW_THROUGH_A_"
            "TWO_PORT_THROAT"
        )
        hi = t5['packet_above_barrier']
        verdict = (
            "DERIVED (the argument is in "
            "docs/wormhole_network_confirmation.md).\n\n"
            "THE MECHANISM. One retarded C-wave: source -> future "
            f"antipode (t = {t2['t_absorb']:.3f}) -> mouth A's "
            "interface -> interior loops (echo train damped by "
            f"|r_in|^2 = {t2['echo_damping_predicted']:.3f}, every "
            "loop future-directed) -> mouth B's OWN interface -> "
            f"emergence at t = {t2['t_emerge']:.3f} (the global PAST, "
            f"via Delta_BA = {t2['delta_BA']:.3f}) -> forward return "
            f"-> the crossing (t_return = {t2['t_return']:.0e}). "
            "U_BA = e^{i w Delta_BA} derived (identity "
            f"{t2['transfer_factor_identity_error']:.0e}).\n\n"
            "THE TWO-PORT THROAT. t_net = t_A t_B/(1 - r_inA r_inB "
            "e^{2 i w tau}) - validated against the DIRECT glued-"
            f"barrier solve to {t4['direct_two_barrier_validation']:.0e} "
            "in |T|^2 across resonances. RESONANT CONFIRMATION: at "
            f"single-port T = {t4['single_port_T']:.2f}, the interior "
            f"comb transmits {t4['resonant_T_net']:.4f} (perfect) vs "
            f"{t4['off_resonance_T_net']:.4f} off resonance "
            f"({t4['resonant_enhancement']:.0f}x); the Airy average "
            "equals the incoherent echo sum to "
            f"{abs(t4['airy_average'] - t4['airy_predicted']):.0e}.\n\n"
            "THE PROJECTION IS ADVANCED. K = Lambda e^{-i w D'} with "
            "D' < 0, |Lambda| = |t_net|; the closure comb (spacing to "
            f"{t5['comb_spacing_max_error']:.0%}); the engine accepts "
            f"the loop at the comb ({t5['engine_weight_at_closure']:.3f}) "
            f"and rejects it detuned "
            f"({t5['engine_weight_at_quarter_detune']:.2f}). Packets: "
            f"above-barrier emerges at {hi['emergence_peak']:.2f} < 0 "
            f"and lands on the crossing ({hi['return_peak']:.3f} vs "
            f"Wigner {hi['wigner_net']:.3f}); below-barrier confirms "
            f"{t5['packet_on_resonance']['energy_fraction']:.2f} ON "
            "the interior resonance vs "
            f"{t5['packet_off_resonance']['energy_fraction']:.3f} off "
            f"({t5['resonant_confirmation_contrast']:.0f}x contrast), "
            "delayed by the cavity storage time.\n\n"
            "THE POLE. At closure I+ + Lambda I- is the covariant "
            f"pole (|Lambda| = {t6['Lambda_magnitude']:.4f}, deviation "
            f"{t6['pole_deviation_at_closure']:.0e}); detuning the "
            "network IS the #213 deform test (phi = w delta exactly). "
            "The two Compton completion families are two globally "
            "causal wormhole-network path classes."
        )
    else:
        verdict_class = "NETWORK_CONFIRMATION_INCONCLUSIVE"
        verdict = (
            "INCONCLUSIVE. A core check failed; re-examine before "
            "quoting the mechanism."
        )

    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "identification": (
            "The wormhole-network confirmation through a two-port "
            "throat: the advanced half of the Wheeler transaction as "
            "an everywhere-future-directed network traversal - "
            "greybody entry, repeated interior loops, exit through the "
            "far mouth's own interface, emergence in the global past, "
            "return to the crossing; the projected kernel is the "
            "advanced kernel with the composite weight |t_net|, phase "
            "closure is #213's coherence condition, and the interior "
            "comb adds resonant confirmation of soft modes"
        ),
        "executes": (
            "the requested model plus the review items: the second "
            "interface and the repeated-loop physics"
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
    out.append("# The wormhole-network confirmation (PR #216)")
    out.append("")
    out.append(f"**Run:** {s['timestamp_utc']}")
    out.append("")
    out.append(
        "The deliverable is `docs/wormhole_network_confirmation.md` "
        "plus `geometrodynamics/transaction/network.py` - the advanced "
        "confirmation as an explicit, globally causal network "
        "traversal through a two-port throat. *(QFT on the fixed "
        "classical throat geometry, not quantum gravity.)*"
    )
    out.append("")
    out.append("## Test summary")
    out.append("")
    out.append("| # | Test | Key finding | PASS? |")
    out.append("|---|---|---|---|")
    labels = {
        "T1": "the goal: a mechanism for the advanced half",
        "T2": "the ledger + echo train: locally forward, globally past",
        "T3": "the two-port mouth S-matrix (both faces)",
        "T4": "the composite throat, validated directly",
        "T5": "the projection is advanced; resonant confirmation",
        "T6": "the pole at closure; detuning = the #213 deform",
        "T7": "honest scope",
        "T8": "assessment",
    }
    for t in s["tests"]:
        p = "**PASS**" if t["pass"] else "**FAIL**"
        pre = t["name"][:2]
        out.append(f"| {pre} | `{t['name']}` | {labels.get(pre,'-')} | {p} |")
    out.append("")
    t2 = s["tests"][1]
    if 'ledger' in t2:
        out.append("## The traversal ledger")
        out.append("")
        out.append("| leg | global start | global end | local duration |")
        out.append("|---|---:|---:|---:|")
        for leg in t2['ledger']:
            out.append(f"| {leg['leg']} | {leg['t_start_global']:.3f} | "
                       f"{leg['t_end_global']:.3f} | "
                       f"{leg['local_duration']:.3f} |")
        out.append("")
        out.append("## The echo train (below-barrier, w = 0.5)")
        out.append("")
        out.append("| k | local duration | global emergence | amplitude |")
        out.append("|---:|---:|---:|---:|")
        for e in t2['echo_train']:
            out.append(f"| {e['k']} | {e['local_duration']:.2f} | "
                       f"{e['t_emerge_global']:.3f} | "
                       f"{e['amplitude']:.4f} |")
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
    out = here / "runs" / f"{ts}_wormhole_network_confirm_probe"
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
