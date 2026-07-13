"""
The wormhole-network confirmation: the advanced half of the Wheeler
transaction as a globally causal, everywhere-future-directed network
traversal (PR #216).

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
    -> transmit into the throat                     (#215 greybody)
    -> traverse the throat forward                  (tau_th > 0 local)
    -> emerge through a clock-offset mouth          (t = pi + tau_th
                                                     + Delta_BA < 0:
                                                     differential aging,
                                                     Morris-Thorne-
                                                     Yurtsever)
    -> propagate forward again on S^3               (back to psi = 0)
    -> intersect the original particle crossing     (t_return = t_emit)

THE CENTRAL IDENTITY.  Clock continuity forces the throat's
global-frame transfer factor to be t_AB(w) e^{i w Delta_BA}: the
projected kernel between absorption and return is

    K = Lambda(w) * e^{-i w D'} ,   D' = tau_th + Delta_BA + d_B < 0 ,
    Lambda(w) = t_AB(w) e^{i w Delta_BA} * decorations ,

i.e. the retarded phase rule analytically continued to a NEGATIVE
exterior interval - THE ADVANCED KERNEL - with greybody weight
|Lambda| = sqrt(T).  #213's coherence requirement (relative phase
phi = 0) becomes the network's PHASE-CLOSURE condition
arg t_AB(w) + w Delta_BA = 0 (mod 2 pi), satisfied on a discrete comb
of frequencies: the transaction selects its modes.  The two Compton
completion families are two globally causal path classes: the direct
exterior path (the offer ordering I+) and the network path (the
confirmation ordering I-), coherent exactly at closure.

Tests:
  T1. Goal.
  T2. The ten-step traversal ledger: local future-directedness
      everywhere, global reach into the past, time closure at the
      crossing, the derived transfer factor U_BA = e^{i w Delta_BA}.
  T3. The greybody amplitude t_AB(w): |t|^2 matches #215, flux
      closure, the free-throat null test (V = 0 => t = 1), the Wigner
      delay.
  T4. The advanced projection: K = Lambda e^{-i w D'} with D' < 0 and
      |Lambda| = sqrt(T); phase closure on a discrete comb (spacing
      2 pi/|Delta + tau_W|); the engine weld (retro_phase_match
      accepts the network loop at closure).
  T5. The live packet: emerges in the global past, intersects the
      crossing within the Wigner delay, returns the greybody energy
      fraction; the below-barrier packet confirms weakly (partial
      confirmation = the #215 IR transparency).
  T6. The pole structure: I+ + Lambda I- reconstructs the covariant
      pole exactly at closure; detuning the network reproduces the
      #213 deform diagnostic (the refutation edge now has a geometric
      knob); partial confirmation splits into coherent-pole share
      (1+sqrt(T))/2 and deficit (1-sqrt(T))/2.
  T7. Honest scope.
  T8. Assessment.

Verdict:
  THE_ADVANCED_HALF_IS_A_FUTURE_DIRECTED_NETWORK_TRAVERSAL_THE_CLOCK_
  OFFSET_MOUTH_IS_THE_ANALYTIC_CONTINUATION_AND_PHASE_CLOSURE_IS_THE_
  COHERENCE_CONDITION
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np
from scipy.integrate import solve_ivp

from geometrodynamics.transaction import (
    NetworkMouth,
    NetworkThroat,
    closure_offset,
    network_confirmation,
    projected_kernel,
    retro_phase_match,
    traverse_throat,
)

_CACHE: dict = {}
_RH = 1.0

# ========================================================================
# SECTION A - the complex greybody amplitude (PR #215 solver, plane-wave
# referenced in the tortoise coordinate so that V = 0 gives t = 1)
# ========================================================================


def _v_of_r(r, l: int, strength: float = 1.0):
    f = 1.0 - (_RH / r) ** 2
    return strength * f * ((l * (l + 2) + 0.75) / r ** 2
                           + 2.25 * _RH ** 2 / r ** 4)


def _x_of_r(r):
    return r + (_RH / 2) * np.log((r - _RH) / (r + _RH))


def greybody_amplitude(w: float, l: int = 0, strength: float = 1.0,
                       rtol: float = 3e-10) -> tuple:
    """Complex transmission amplitude t(w): psi -> t e^{-iwx} toward the
    horizon, psi -> e^{-iwx} + r e^{+iwx} outside.  Referenced to plane
    waves in x, so a switched-off potential (strength = 0) gives t = 1
    exactly."""
    key = ('g', w, l, strength)
    if key in _CACHE:
        return _CACHE[key]
    r0 = _RH * (1 + 1e-7)
    x0 = float(_x_of_r(r0))
    r_far = max(60.0 / w, 50.0 * w, 30.0)
    x_far = float(_x_of_r(r_far))

    def rhs(x, y):
        pr, pi_, qr, qi, r = y
        vv = _v_of_r(r, l, strength)
        f = 1 - (_RH / r) ** 2
        return [qr, qi, (vv - w ** 2) * pr, (vv - w ** 2) * pi_, f]

    y0 = [math.cos(-w * x0), math.sin(-w * x0),
          w * math.sin(-w * x0), -w * math.cos(-w * x0), r0]
    sol = solve_ivp(rhs, (x0, x_far), y0, rtol=rtol, atol=1e-13,
                    dense_output=True, method="DOP853")
    xs = x_far - np.linspace(0.0, 3.7, 8)
    rows, vals = [], []
    for xm in xs:
        pr, pi_, _, _, _ = sol.sol(xm)
        rows.append([np.exp(-1j * w * xm), np.exp(1j * w * xm)])
        vals.append(pr + 1j * pi_)
    (alpha, beta), *_ = np.linalg.lstsq(np.array(rows), np.array(vals),
                                        rcond=None)
    out = (complex(1.0 / alpha), complex(beta / alpha))
    _CACHE[key] = out
    return out


def wigner_delay(w: float, l: int = 0, h: float = 1e-3) -> float:
    tp, _ = greybody_amplitude(w + h, l)
    tm, _ = greybody_amplitude(w - h, l)
    return float((np.angle(tp) - np.angle(tm)) / (2 * h))


# ── the standard network of this probe ──────────────────────────────────

_D_A = math.pi          # source -> future antipodal mouth A
_D_B = math.pi          # past mouth B -> source
_TAU = 0.8              # throat traversal proper time


def standard_throat(delta: Optional[float] = None) -> NetworkThroat:
    if delta is None:
        delta = closure_offset(_D_A, _D_B, _TAU)
    A = NetworkMouth("A", psi=math.pi, link_id="L1", clock_offset=0.0)
    B = NetworkMouth("B", psi=math.pi, link_id="L1", clock_offset=delta)
    return NetworkThroat(A, B, tau_th=_TAU,
                         t_AB=lambda w: greybody_amplitude(w)[0])


# ========================================================================
# T1. Goal
# ========================================================================


def test_T1_goal() -> dict:
    return {
        'name': 'T1_goal',
        'description': (
            'Replace advanced_confirm_amplitude with an explicit '
            'network traversal: one retarded C-wave to the future '
            'antipode, greybody transmission into the throat, forward '
            'traversal, emergence through a clock-offset mouth deep '
            'in the global past, forward propagation back to the '
            'original crossing - showing the projected response is '
            'advanced, every segment is locally future-directed, the '
            'closures hold, and the effective kernel has the #213 '
            'pole structure. The advanced half of a Wheeler '
            'transaction gets a classical geometric mechanism.'
        ),
        'pass': True,
    }


# ========================================================================
# T2. The ten-step traversal ledger
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

    # the derived transfer factor: global-frame factor of the throat
    # leg equals t_AB(w) * e^{i w Delta_BA} at every frequency
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
          and transfer_err < 1e-12)
    return {
        'name': 'T2_traversal_ledger',
        'description': (
            'every leg future-directed in its own clock; the global '
            'past reached only through the frozen mouth offset; time '
            'closure at the crossing; U_BA = e^{i w Delta_BA} derived'
        ),
        'ledger': ledger,
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
# T3. The greybody amplitude
# ========================================================================


def test_T3_greybody_amplitude() -> dict:
    if 'T3' in _CACHE:
        return _CACHE['T3']
    # (i) |t|^2 against the #215 committed values (Hankel-matched)
    ref = {0.1: 1.663731e-03, 0.5: 3.331217e-01,
           1.0: 9.629316e-01, 2.0: 0.9998469}
    rel_err, flux_err = 0.0, 0.0
    for w, Tref in ref.items():
        t, r = greybody_amplitude(w)
        rel_err = max(rel_err, abs(abs(t) ** 2 - Tref) / Tref)
        flux_err = max(flux_err, abs(abs(t) ** 2 + abs(r) ** 2 - 1))

    # (ii) the free-throat null test: strength = 0 => t = 1 exactly
    t_null, _ = greybody_amplitude(1.3, strength=0.0)
    null_err = abs(t_null - 1.0)

    # (iii) the Wigner delay: large below the barrier, small above
    tw = {w: wigner_delay(w) for w in (0.5, 1.0, 3.0)}

    ok = (rel_err < 1e-3 and flux_err < 3e-4 and null_err < 1e-8
          and tw[0.5] > tw[1.0] > tw[3.0] > 0)
    out = {
        'name': 'T3_greybody_amplitude',
        'description': (
            'the complex t_AB(w) from the #215 scattering problem, '
            'plane-wave referenced (V = 0 gives t = 1 exactly); '
            'flux-closed; the Wigner delay of the throat'
        ),
        'match_to_215_relative': float(rel_err),
        'flux_closure': float(flux_err),
        'free_throat_null_error': float(null_err),
        'wigner_delay': {str(k): float(v) for k, v in tw.items()},
        'pass': bool(ok),
    }
    _CACHE['T3'] = out
    return out


# ========================================================================
# T4. The advanced projection and phase closure
# ========================================================================


def test_T4_advanced_projection() -> dict:
    if 'T4' in _CACHE:
        return _CACHE['T4']
    th = standard_throat()

    # (i) the projected kernel: pastward interval, greybody weight
    interval_ok, weight_err = True, 0.0
    for w in (0.8, 2.0, 3.2):
        pk = projected_kernel(th, w, _D_B)
        interval_ok = interval_ok and pk['interval'] < 0
        weight_err = max(weight_err,
                         abs(abs(pk['confirmation_weight'])
                             - pk['greybody_magnitude']))

    # (ii) phase closure on a discrete comb: arg Lambda(w) = arg t(w)
    # + w Delta_BA = 0 (mod 2 pi); roots and their spacing
    ws = np.linspace(2.2, 4.2, 81)
    argl = np.array([np.angle(greybody_amplitude(float(w))[0]
                              * np.exp(1j * w * th.delta_BA))
                     for w in ws])
    roots = []
    for i in range(len(ws) - 1):
        d = argl[i + 1] - argl[i]
        if argl[i] * argl[i + 1] < 0 and abs(d) < math.pi:  # true crossing
            roots.append(float(ws[i] - argl[i]
                               * (ws[i + 1] - ws[i]) / d))
    spacings = np.diff(roots)
    tau_w_mid = wigner_delay(3.0)
    spacing_pred = 2 * math.pi / abs(th.delta_BA + tau_w_mid)
    spacing_err = float(np.max(np.abs(spacings / spacing_pred - 1)))

    # (iii) the engine weld: at time closure the full loop amplitude
    # equals Lambda(w); retro_phase_match accepts it at a closing root
    # and rejects it at quarter detuning
    w_close = roots[0]
    tr = network_confirmation(th, w_close, 0.0, _D_A, _D_B)
    _, weight_close = retro_phase_match(1.0 + 0.0j, tr.amp)
    w_off = w_close + spacing_pred / 4        # arg Lambda ~ pi/2
    tr_off = network_confirmation(th, w_off, 0.0, _D_A, _D_B)
    _, weight_off = retro_phase_match(1.0 + 0.0j, tr_off.amp)

    ok = (interval_ok and weight_err < 1e-12
          and len(roots) >= 2 and spacing_err < 0.1
          and weight_close > 0.99 and weight_off < 0.5)
    out = {
        'name': 'T4_advanced_projection',
        'description': (
            'K = Lambda e^{-iwD\'} with D\' < 0 and |Lambda| = '
            'sqrt(T): the advanced kernel with greybody weight; phase '
            'closure selects a discrete frequency comb; the engine\'s '
            'phase closure accepts the network loop at the comb'
        ),
        'pastward_interval': bool(interval_ok),
        'weight_equals_greybody_error': float(weight_err),
        'closure_comb_roots': [float(r) for r in roots],
        'comb_spacing_predicted': float(spacing_pred),
        'comb_spacing_max_error': spacing_err,
        'engine_weight_at_closure': float(weight_close),
        'engine_weight_at_quarter_detune': float(weight_off),
        'pass': bool(ok),
    }
    _CACHE['T4'] = out
    return out


# ========================================================================
# T5. The live packet
# ========================================================================


def _packet_run(w0: float, sig: float, npts: int = 61):
    th = standard_throat()
    D = _D_A + _TAU + th.delta_BA + _D_B          # = 0 at time closure
    t_e = _D_A + _TAU + th.delta_BA               # emergence time
    ws = np.linspace(w0 - 4 * sig, w0 + 4 * sig, npts)
    c = np.exp(-0.5 * ((ws - w0) / sig) ** 2)
    t_amp = np.array([greybody_amplitude(float(w))[0] for w in ws])

    tgrid = np.linspace(-4, 4, 4001)
    # incident packet at the source (peak t = 0)
    f0 = np.sum(c[None, :] * np.exp(-1j * np.outer(tgrid, ws)), axis=1)
    # returned packet: filtered and shifted by the closed network
    g = np.sum(c[None, :] * t_amp[None, :]
               * np.exp(-1j * np.outer(tgrid - D, ws)), axis=1)
    # emergence envelope at mouth B (peak expected at t_e + tau_W)
    sgrid = np.linspace(t_e - 4, t_e + 4, 4001)
    h = np.sum(c[None, :] * t_amp[None, :]
               * np.exp(-1j * np.outer(sgrid - t_e, ws)), axis=1)

    e_frac = float(np.trapezoid(np.abs(g) ** 2, tgrid)
                   / np.trapezoid(np.abs(f0) ** 2, tgrid))
    mean_T = float(np.sum(c ** 2 * np.abs(t_amp) ** 2) / np.sum(c ** 2))
    return {
        'return_peak': float(tgrid[np.argmax(np.abs(g))]),
        'emergence_peak': float(sgrid[np.argmax(np.abs(h))]),
        't_emergence_geometric': float(t_e),
        'wigner': wigner_delay(w0),
        'energy_fraction': e_frac,
        'mean_T': mean_T,
    }


def test_T5_live_packet() -> dict:
    if 'T5' in _CACHE:
        return _CACHE['T5']
    hi = _packet_run(3.0, 0.3)          # above the barrier
    lo = _packet_run(0.5, 0.1, npts=41)  # below the barrier

    ok = (hi['emergence_peak'] < 0.0                       # in the past
          and abs(hi['emergence_peak']
                  - (hi['t_emergence_geometric'] + hi['wigner'])) < 0.05
          and abs(hi['return_peak'] - hi['wigner']) < 0.05  # the crossing
          and abs(hi['energy_fraction'] - hi['mean_T']) < 1e-3
          and hi['energy_fraction'] > 0.99
          and 0.2 < lo['energy_fraction'] < 0.5)            # weak confirm
    out = {
        'name': 'T5_live_packet',
        'description': (
            'a Gaussian packet run through the closed network: emerges '
            'in the global past, intersects the original crossing '
            'within the Wigner delay, returns the greybody energy '
            'fraction; a below-barrier packet confirms weakly (the '
            '#215 IR transparency as partial confirmation)'
        ),
        'above_barrier': hi,
        'below_barrier': lo,
        'pass': bool(ok),
    }
    _CACHE['T5'] = out
    return out


# ========================================================================
# T6. The pole structure
# ========================================================================


def test_T6_pole_structure() -> dict:
    th = standard_throat()
    t4 = test_T4_advanced_projection()
    w_close = t4['closure_comb_roots'][1]     # a closing mode, T ~ 1
    eps = 1e-3

    def ordering_plus(dl, w):
        return -(1 / (2 * w)) / (dl - w + 1j * eps)

    def ordering_minus(dl, w):
        return +(1 / (2 * w)) / (dl + w - 1j * eps)

    dgrid = np.linspace(-5, 5, 801)
    mask = np.abs(np.abs(dgrid) - w_close) > 0.3

    # (i) at closure: Lambda real positive, |Lambda| = sqrt(T) ~ 1;
    # K = I+ + Lambda I- is the covariant pole to the Lambda deficit
    lam = (greybody_amplitude(w_close)[0]
           * np.exp(1j * w_close * th.delta_BA))
    K = ordering_plus(dgrid, w_close) + lam * ordering_minus(dgrid, w_close)
    pole_exact = -(w_close - 1j * eps) / (
        w_close * (dgrid ** 2 - (w_close - 1j * eps) ** 2))
    dev_closed = float(np.max(np.abs((K - pole_exact) / pole_exact)[mask]))
    lam_phase = float(abs(np.angle(lam)))
    lam_mag = float(abs(lam))

    # (ii) the geometric deform knob: detune Delta_BA by delta ->
    # relative phase phi = w*delta; the pole-form deviation matches the
    # #213 T5 diagnostic at the same phi
    deform = {}
    for delta in (0.1, 0.5):
        th_d = standard_throat(delta=th.delta_BA + delta)
        lam_d = (greybody_amplitude(w_close)[0]
                 * np.exp(1j * w_close * th_d.delta_BA))
        phi_geo = float(np.angle(lam_d / lam))
        phi_pred = ((w_close * delta + math.pi) % (2 * math.pi)) - math.pi
        K_d = (ordering_plus(dgrid, w_close)
               + lam_d * ordering_minus(dgrid, w_close))
        dev_d = float(np.max(np.abs((K_d - pole_exact) / pole_exact)[mask]))
        # the #213 diagnostic at the same phi (unit-magnitude deform)
        K_213 = (ordering_plus(dgrid, w_close)
                 + np.exp(1j * phi_geo)
                 * ordering_minus(dgrid, w_close))
        dev_213 = float(np.max(np.abs((K_213 - pole_exact)
                                      / pole_exact)[mask]))
        deform[delta] = {'phi_geometric': phi_geo,
                         'phi_predicted_w_delta': phi_pred,
                         'pole_deviation': dev_d,
                         'pole_deviation_213_diag': dev_213}

    # (iii) partial confirmation below the barrier: K = I+ + sqrt(T) I-
    # splits into coherent-pole share (1+sqrt T)/2 + deficit (1-sqrt T)/2
    w_low = 0.5
    s = math.sqrt(abs(greybody_amplitude(w_low)[0]) ** 2)
    Kp = (ordering_plus(dgrid, w_low)
          + s * ordering_minus(dgrid, w_low))
    share = (1 + s) / 2 * (ordering_plus(dgrid, w_low)
                           + ordering_minus(dgrid, w_low))
    deficit = (1 - s) / 2 * (ordering_plus(dgrid, w_low)
                             - ordering_minus(dgrid, w_low))
    split_err = float(np.max(np.abs(Kp - share - deficit)))

    ok = (dev_closed < 1e-5
          and lam_phase < 0.02 and lam_mag > 0.999
          and all(abs(d['phi_geometric'] - d['phi_predicted_w_delta'])
                  < 1e-9 for d in deform.values())
          and all(abs(d['pole_deviation'] - d['pole_deviation_213_diag'])
                  < 0.02 * d['pole_deviation_213_diag'] + 1e-6
                  for d in deform.values())
          and split_err < 1e-12)
    return {
        'name': 'T6_pole_structure',
        'description': (
            'the two completion families as path classes: at network '
            'closure I+ + Lambda I- is the covariant pole; network '
            'detuning IS the #213 deform test (phi = w delta); partial '
            'confirmation splits coherent share/deficit by sqrt(T)'
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
        'Single-mode scalar, zonal (1D transit) reduction throughout; '
        'the greybody is the l = 0 channel.',
        'The elastic case (equal mouth clock rates) is analyzed; '
        'unequal rates red/blueshift the traversal '
        '(emergent_frequency) and are only unit-tested.',
        'The engine\'s advanced_confirm_amplitude is retained; '
        'network_confirmation is the mechanism-level replacement and '
        'the T4 weld shows the engine\'s phase closure accepts it - '
        'full engine migration is staged, not performed.',
        'Phase closure selects a frequency comb; which tower modes '
        'land on the comb depends on the network parameters - mode '
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
    t3 = test_T3_greybody_amplitude()
    t4 = test_T4_advanced_projection()
    t5 = test_T5_live_packet()
    t6 = test_T6_pole_structure()
    core = all(t['pass'] for t in (t2, t3, t4, t5, t6))
    assessment = (
        'The advanced half of the Wheeler transaction now has a '
        'classical geometric mechanism: a retarded wave, absorbed at '
        'the future antipodal mouth, transmitted with the throat\'s '
        'own greybody amplitude, traversed forward, and re-emitted '
        'through a clock-offset mouth in the global past, projects '
        'onto the exterior as EXACTLY the advanced kernel with weight '
        'sqrt(T) - every segment locally future-directed, the past '
        'reached only through the frozen differential aging of the '
        'mouths.  #213\'s coherent relative phase is re-derived as the '
        'network\'s phase-closure condition, its deform test becomes a '
        'geometric detuning knob, and the two Compton completion '
        'families become two globally causal path classes through the '
        'wormhole network.'
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
        test_T3_greybody_amplitude(),
        test_T4_advanced_projection(),
        test_T5_live_packet(),
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
            "_PHASE_CLOSURE_IS_THE_COHERENCE_CONDITION"
        )
        hi = t5['above_barrier']
        verdict = (
            "DERIVED (the argument is in "
            "docs/wormhole_network_confirmation.md).\n\n"
            "THE MECHANISM. One retarded C-wave: source -> future "
            f"antipode (t = {t2['t_absorb']:.3f}) -> greybody "
            "transmission -> forward throat traversal -> emergence at "
            f"t = {t2['t_emerge']:.3f} (the global PAST, via the "
            f"frozen mouth offset Delta_BA = {t2['delta_BA']:.3f}) -> "
            "forward return -> the original crossing (t_return = "
            f"{t2['t_return']:.0e}). Every leg locally future-"
            "directed; the derived transfer factor U_BA = "
            "e^{i w Delta_BA} (identity to "
            f"{t2['transfer_factor_identity_error']:.0e}).\n\n"
            "THE PROJECTION IS ADVANCED. The absorption -> return "
            "kernel is Lambda(w) e^{-i w D'} with D' < 0 and |Lambda| "
            "= sqrt(T): the retarded rule continued to a negative "
            "exterior interval - the advanced kernel with the #215 "
            "greybody weight. Phase closure arg t + w Delta = 0 (mod "
            f"2 pi) selects a comb (spacing 2 pi/|Delta + tau_W| to "
            f"{t4['comb_spacing_max_error']:.0%}); the engine's "
            "retro_phase_match accepts the loop at the comb "
            f"({t4['engine_weight_at_closure']:.3f}) and rejects it "
            f"detuned ({t4['engine_weight_at_quarter_detune']:.2f}). "
            "Live packet: emerges at "
            f"{hi['emergence_peak']:.2f} (< 0), intersects the "
            f"crossing at {hi['return_peak']:.3f} (the Wigner delay "
            f"{hi['wigner']:.3f}), energy fraction "
            f"{hi['energy_fraction']:.4f} = mean T.\n\n"
            "THE POLE. At closure I+ + Lambda I- is the covariant "
            f"pole (|Lambda| = {t6['Lambda_magnitude']:.4f}, phase "
            f"{t6['Lambda_phase']:.3f}); detuning the network IS the "
            "#213 deform test (phi = w delta exactly; deviations match "
            "the #213 diagnostic to 2%); below the barrier the "
            "confirmation is partial with the coherent share "
            "(1+sqrt T)/2. The two Compton completion families are "
            "two globally causal wormhole-network path classes."
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
            "The wormhole-network confirmation: the advanced half of "
            "the Wheeler transaction as an everywhere-future-directed "
            "network traversal - greybody transmission at the future "
            "antipodal mouth, forward traversal, emergence through a "
            "clock-offset mouth in the global past, and return to the "
            "original crossing; the projected kernel is the advanced "
            "kernel with weight sqrt(T), and phase closure is #213's "
            "coherence condition"
        ),
        "executes": (
            "the requested model: advanced_confirm_amplitude replaced "
            "by an explicit network traversal"
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
        "traversal. *(QFT on the fixed classical throat geometry, not "
        "quantum gravity.)*"
    )
    out.append("")
    out.append("## Test summary")
    out.append("")
    out.append("| # | Test | Key finding | PASS? |")
    out.append("|---|---|---|---|")
    labels = {
        "T1": "the goal: a mechanism for the advanced half",
        "T2": "the ten-step ledger: locally forward, globally past",
        "T3": "the complex greybody amplitude + Wigner delay",
        "T4": "the projection is advanced; phase closure = comb",
        "T5": "the live packet intersects the crossing",
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
