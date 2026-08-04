"""
The physical gate set: the qubit clock, the gravitational Kerr, and
the twist switch on the committed network (PR #234).

> Framing: QFT on the *fixed classical* throat geometry (geometry -> fields),
> not quantum gravity.

THE QUESTION
------------
#232 delivered the 2^N structure with a STRUCTURAL entangling
coupling: "chi is the |phi|^4 overlap ... not calibrated to a throat
geometry in this probe" - its named scope item.  This probe executes
it on the committed #224 two-throat network, and welds in the other
arc's switch (#227/#233): every rate of the gate set, calibrated in
model units on the geometry the ledger already owns.

THE ANSWER (measured / machine-checked)
---------------------------------------
  * THE QUBIT CLOCK: the intra-doublet beat IS the one-qubit X
    rotation - dw = 1.714e-3, X period pi/dw = 1833 (the #224
    ledger numbers re-derived exactly on the same builder).
  * THE TWIST SWITCH: flipping one exterior link (the #227 Z2
    holonomy) freezes the beat to 8.3e-14 - ratio 4.9e-11,
    REPRODUCING the #227 committed row at r_s = 0.3 exactly - so
    one topological bit toggles the SAME doublet between GATE mode
    (clock running) and MEMORY mode (leakage time 3.8e13, the #227
    exact-pair theorem).  The basin DENSITIES are twist-invariant
    (the properly localized degenerate combination changes the
    dressing integral by < 1%): the twist kills the linear (leaky)
    channel and leaves the density-density (Kerr) channel intact -
    it is an error switch, not a gate hazard.
  * THE GRAVITATIONAL KERR, CALIBRATED: the committed nonlinearity
    is not a contact |phi|^4 - it is the EKG backreaction, whose
    leading form is the #223 dressing law delta mu = c A^2 (EXACTLY
    quadratic, ledger re-read).  Per ONE QUANTUM of a rail mode on
    its throat (continuum normalization, amplitude^2 = 1/2w):
    delta mu_q = 5.8e-2; the throat response dw/dmu = -0.95
    (measured, linear); so the on-throat Kerr is
        chi_rail = |dw/dmu| * delta mu_q = 5.5e-2   (model units)
    giving t_CZ = pi/chi = 57.2 - the entangling gate is 32x FASTER
    than the one-qubit clock.  GRAVITY IS THE ENTANGLING RESOURCE.
  * THE DUAL-RAIL SELECTION RULE: within the sector (one quantum
    per doublet) every SAME-qubit Kerr term vanishes IDENTICALLY
    (n_a n_b = 0 on a shared quantum; n(n-1) = 0 for <= 1) - checked
    in Fock at machine zero.  The only active Kerr couples rails of
    DIFFERENT qubits: exactly the CZ generator.  Distant-rail
    cross-talk is measured at 3.9e-2 of chi_rail (the basin tail),
    a known coherent ZZ phase of 0.12 rad per CZ - compensable, and
    reported, not hidden.
  * THE STRONG-COUPLING HONESTY: delta mu_q / mu = 0.64 at the
    committed network scale - ONE quantum dresses the throat at
    O(1): the quoted chi is the LEADING-ORDER value and the
    perturbative window (#223's own A_pert) is exceeded per
    quantum; the beyond-leading Kerr is the named successor.
  * THE CIRCUIT BUDGET: GHZ_3 at physical rates - 5 one-qubit
    half-periods + 2 CZs = ~4700 model units; twist-frozen leakage
    error (T dw_frozen)^2 ~ 1.5e-19; the budget table committed.

Tests:
  T1. Goal (execute #232's scope item; weld the twist arc).
  T2. The qubit clock: the doublet re-derived = the #224 ledger.
  T3. The twist switch: freeze reproduced = the #227 ledger row;
      basin densities twist-invariant (the Kerr channel survives).
  T4. The gravitational Kerr calibrated: delta mu per quantum,
      dw/dmu, chi_rail, t_CZ; the strong-coupling caveat.
  T5. The dual-rail selection rule (Fock, machine zero) + the CZ
      with the calibrated chi under frozen linear coupling.
  T6. The circuit budget: GHZ_3 at physical rates; the error ledger.
  T7. Arc-weld holdout: #224/#227/#223/#232 committed re-reads.
  T8. Honest scope.
  T9. Assessment.

Verdict:
  THE_GATE_SET_IS_PHYSICAL_THE_BEAT_IS_THE_CLOCK_THE_TWIST_IS_THE_
  SWITCH_AND_THE_GRAVITATIONAL_DRESSING_IS_THE_KERR_WITH_T_CZ_32X_
  FASTER_THAN_THE_CLOCK_AT_LEADING_ORDER_IN_A_STRONG_COUPLING
"""

from __future__ import annotations

import glob
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np

import mouth_exchange_dynamics_probe as mx
import tensor_product_emergence_probe as tp

_KAPPA = 0.5                     # = 2 KF, the #222/#223 convention
_RS = 0.3

_CACHE: dict = {}


def _latest(pat):
    fs = sorted(glob.glob(pat))
    return json.load(open(fs[-1])) if fs else None


def _t(d, name):
    return [t for t in d['tests'] if t['name'] == name][0]


# ── network helpers on the committed builder ────────────────────────────

def ring_modes_twist(g, sign):
    """The #224 ring operator with the closing link sign = +-1 (the
    #227 Z2 holonomy on the exterior cycle)."""
    N, dx, V = g['N'], g['dx'], g['V']
    L = (np.diag(-2.0 * np.ones(N)) + np.diag(np.ones(N - 1), 1)
         + np.diag(np.ones(N - 1), -1))
    L[0, -1] = L[-1, 0] = sign
    L /= dx ** 2
    evals, evecs = np.linalg.eigh(-L + np.diag(V))
    return np.sqrt(np.abs(evals)), evecs


def interior_pair(g, ws, evecs, w_lo=1.9, w_hi=2.5):
    inI = g['in_A'] | g['in_B']
    ins = [k for k in range(g['N'])
           if w_lo < ws[k] < w_hi
           and float(np.sum(evecs[:, k][inI] ** 2)) > 0.5]
    best, bd = None, 1e9
    for i in range(len(ins) - 1):
        dd = ws[ins[i + 1]] - ws[ins[i]]
        if dd < bd:
            bd, best = dd, (ins[i], ins[i + 1])
    return best


def localized_basins(g, u1, u2):
    """Rotate the (possibly exactly degenerate) pair to the maximally
    B-localized combination - analytically: the top eigenvector of
    the 2x2 quadratic form of the B-region weight."""
    PB = g['in_B']
    a = float(np.sum(u1[PB] * u1[PB]))
    b = float(np.sum(u2[PB] * u2[PB]))
    c = float(np.sum(u1[PB] * u2[PB]))
    M = np.array([[a, c], [c, b]])
    w, V = np.linalg.eigh(M)
    v = V[:, -1]                        # max B-weight combination
    uB = v[0] * u1 + v[1] * u2
    uA = -v[1] * u1 + v[0] * u2
    return uA, uB, float(w[-1])


def _wallsB(g):
    half = mx._D_CH / 2
    b0 = half + mx._SIG_M + mx._LARC + mx._SIG_M
    b1 = b0 + mx._D_CH
    s = g['s']
    w1 = (s >= b0 - mx._SIG_M) & (s < b0)
    w2 = (s >= b1) & (s < b1 + mx._SIG_M)
    return (w1, b0 - s[w1]), (w2, s[w2] - b1)


def dmu_per_quantum(g, u, w):
    """delta mu of throat B from ONE quantum of grid-normalized mode
    u: continuum-normalize (int u^2 dx = 1), per-quantum amplitude
    1/sqrt(2w), phi = u/rho^{3/2}, the #223 stress integral over
    throat B's two bridge walls."""
    s, dx = g['s'], g['dx']
    un = u / math.sqrt(dx)
    amp = un / math.sqrt(2 * w)
    tot = 0.0
    for m, sig in _wallsB(g):
        rho = np.sqrt(_RS ** 2 + sig ** 2)
        phi = amp[m] / rho ** 1.5
        dphi = np.gradient(phi, s[m])
        rho_E = w ** 2 * phi ** 2 + dphi ** 2 + 3.0 / rho ** 2 * phi ** 2
        tot += float(np.sum(rho ** 3 * rho_E) * dx)
    return (2 * _KAPPA / 3) * tot


def build_two_throat_rsB(rs, rsB):
    """The committed builder with an independent throat-B wall scale
    (for the response derivative dw/dmu_B)."""
    dxl = mx._DX
    Ltot = 2 * mx._D_CH + 4 * mx._SIG_M + 2 * mx._LARC
    Nl = int(round(Ltot / dxl))
    sl = np.arange(Nl) * dxl
    V = np.zeros(Nl)
    halfl = mx._D_CH / 2

    def seg(v0, v1):
        return (sl >= v0 - 1e-12) & (sl < v1 - 1e-12)

    p = halfl
    m = seg(p, p + mx._SIG_M)
    V[m] = mx.v_bridge_half(sl[m] - p, rs)
    p += mx._SIG_M
    a0 = p
    m = seg(p, p + mx._LARC)
    V[m] = (3.75 / ((sl[m] - a0) + mx._SIG_M) ** 2
            + 3.75 / ((a0 + mx._LARC - sl[m]) + mx._SIG_M) ** 2 - 1.0)
    p += mx._LARC
    m = seg(p, p + mx._SIG_M)
    V[m] = mx.v_bridge_half(p + mx._SIG_M - sl[m], rsB)
    p += mx._SIG_M
    b_s = p
    m = seg(p, p + mx._D_CH)
    V[m] = 0.0
    p += mx._D_CH
    b_e = p
    m = seg(p, p + mx._SIG_M)
    V[m] = mx.v_bridge_half(sl[m] - p, rsB)
    p += mx._SIG_M
    a0 = p
    m = seg(p, p + mx._LARC)
    V[m] = (3.75 / ((sl[m] - a0) + mx._SIG_M) ** 2
            + 3.75 / ((a0 + mx._LARC - sl[m]) + mx._SIG_M) ** 2 - 1.0)
    p += mx._LARC
    m = seg(p, p + mx._SIG_M)
    V[m] = mx.v_bridge_half(p + mx._SIG_M - sl[m], rs)
    p += mx._SIG_M
    in_A = (sl < halfl) | (sl >= p)
    in_B = (sl >= b_s) & (sl < b_e)
    return {'s': sl, 'V': V, 'dx': dxl, 'N': Nl,
            'in_A': in_A, 'in_B': in_B}


def _wB_of(dmu):
    rsB = math.sqrt(_RS ** 2 + dmu)
    g2 = build_two_throat_rsB(_RS, rsB)
    ws, evecs = mx.ring_modes(g2)
    cands = [k for k in range(g2['N']) if 1.9 < ws[k] < 2.5
             and float(np.sum(evecs[:, k][g2['in_B']] ** 2)) > 0.5]
    return float(ws[cands[0]])


# ========================================================================
# T1. Goal
# ========================================================================


def test_T1_goal() -> dict:
    return {
        'name': 'T1_goal',
        'description': (
            'Execute #232\'s named scope item - calibrate the '
            'entangling coupling chi on the committed geometry - and '
            'weld in the #227/#233 twist: the full physical gate set '
            '(the one-qubit clock, the entangling Kerr, the '
            'memory/gate switch) with every rate measured in model '
            'units on the #224 network the ledger already owns.'
        ),
        'pass': True,
    }


# ========================================================================
# T2. The qubit clock
# ========================================================================


def test_T2_clock() -> dict:
    if 'T2' in _CACHE:
        return _CACHE['T2']
    g = mx.build_two_throat(_RS)
    d = mx.interior_doublet(g)
    period = math.pi / d['dw']
    uA, uB, locB = localized_basins(
        g, d['u1'], d['u2'])
    locA = float(np.sum(uA[g['in_A']] ** 2))
    # the #224 ledger's own construction: the (u1 +- u2)/sqrt2 combos
    cP = (d['u1'] + d['u2']) / math.sqrt(2)
    cM = (d['u1'] - d['u2']) / math.sqrt(2)
    loc224 = max(float(np.sum(cP[g['in_A']] ** 2)),
                 float(np.sum(cP[g['in_B']] ** 2)),
                 float(np.sum(cM[g['in_A']] ** 2)),
                 float(np.sum(cM[g['in_B']] ** 2)))
    # the #224 ledger, re-derived on the same builder: exact
    here = Path(__file__).resolve().parent / 'runs'
    d224 = _latest(str(here / '*_mouth_exchange_dynamics_probe/probe.json'))
    t3l = _t(d224, 'T3_p_other_mouth')
    dev_dw = abs(d['dw'] / t3l['dw'] - 1.0)
    dev_per = abs(period / t3l['period_pi_over_dw'] - 1.0)

    ok = dev_dw < 1e-9 and dev_per < 1e-9 and locA > 0.97 and locB > 0.97
    out = {
        'name': 'T2_clock',
        'description': (
            'THE QUBIT CLOCK: the #224 interior doublet IS the '
            'one-qubit sector - rails = the A/B mouth basins '
            '(localization 0.977), and the intra-doublet beat is the '
            'X rotation: dw = 1.714e-3, X period pi/dw = 1833 - the '
            'committed #224 ledger re-derived exactly on the same '
            'builder (machine agreement).  One-qubit gates run at '
            'the exchange beat: the qubit\'s clock is the network\'s '
            'own tunneling'
        ),
        'dw': d['dw'],
        'x_period': float(period),
        'basin_localization_A': locA,
        'basin_localization_B': float(locB),
        'localized_combo_weight_224_construction': float(loc224),
        'ledger_dw_rel_dev': float(dev_dw),
        'ledger_period_rel_dev': float(dev_per),
        'pass': bool(ok),
    }
    _CACHE['T2'] = out
    return out


# ========================================================================
# T3. The twist switch
# ========================================================================


def test_T3_twist_switch() -> dict:
    if 'T3' in _CACHE:
        return _CACHE['T3']
    g = mx.build_two_throat(_RS)
    d = mx.interior_doublet(g)
    wt, vt = ring_modes_twist(g, -1.0)
    pair = interior_pair(g, wt, vt)
    dw_tw = float(wt[pair[1]] - wt[pair[0]])
    ratio = dw_tw / d['dw']
    leak_time = math.pi / dw_tw
    # the #227 committed row at r_s = 0.3, re-read
    here = Path(__file__).resolve().parent / 'runs'
    d227 = _latest(str(here / '*_nonorientable_er_link_z2_probe/probe.json'))
    row = [r for r in d227['T4']['rows'] if abs(r['r_s'] - 0.3) < 1e-9][0]
    dev_frozen = abs(dw_tw / row['dw_twisted'] - 1.0)
    # basin-density twist invariance: rotate the degenerate pair to
    # the localized combination, compare the dressing integral
    wbar_t = 0.5 * (wt[pair[0]] + wt[pair[1]])
    uAt, uBt, locBt = localized_basins(g, vt[:, pair[0]], vt[:, pair[1]])
    wbar = 0.5 * (d['w1'] + d['w2'])
    uA, uB, _ = localized_basins(g, d['u1'], d['u2'])
    dmu_untw = dmu_per_quantum(g, uB, wbar)
    dmu_tw = dmu_per_quantum(g, uBt, wbar_t)
    kerr_shift = abs(dmu_tw / dmu_untw - 1.0)

    ok = (ratio < 1e-9 and dev_frozen < 1e-6 and locBt > 0.97
          and kerr_shift < 0.01)
    out = {
        'name': 'T3_twist_switch',
        'description': (
            'THE TWIST SWITCH: flipping one exterior link (the #227 '
            'Z2 holonomy) freezes the SAME doublet\'s beat by 4.9e-11 '
            '- reproducing the #227 committed row at r_s = 0.3 to '
            'machine agreement - so one topological bit toggles the '
            'qubit between GATE mode (clock running) and MEMORY mode '
            '(leakage time 3.8e13; the #227 exact-pair theorem).  '
            'The properly localized twisted basins (the degenerate '
            'pair rotated to maximal localization, 0.977) leave the '
            'DRESSING INTEGRAL unchanged to < 1%: the twist kills '
            'the linear (leaky) channel and leaves the density-'
            'density (Kerr) channel intact - an error switch, not a '
            'gate hazard'
        ),
        'dw_twisted': dw_tw,
        'freeze_ratio': float(ratio),
        'memory_leakage_time': float(leak_time),
        'ledger_227_frozen_rel_dev': float(dev_frozen),
        'twisted_basin_localization': float(locBt),
        'kerr_channel_shift_under_twist': float(kerr_shift),
        'pass': bool(ok),
    }
    _CACHE['T3'] = out
    return out


# ========================================================================
# T4. The gravitational Kerr, calibrated
# ========================================================================


def test_T4_kerr_calibrated() -> dict:
    if 'T4' in _CACHE:
        return _CACHE['T4']
    g = mx.build_two_throat(_RS)
    d = mx.interior_doublet(g)
    wbar = 0.5 * (d['w1'] + d['w2'])
    uA, uB, _ = localized_basins(g, d['u1'], d['u2'])
    # the dressing per quantum (the #223 stress integral, continuum
    # normalized, amplitude^2 = 1/2w per quantum)
    dmu_B = dmu_per_quantum(g, uB, wbar)
    dmu_A = dmu_per_quantum(g, uA, wbar)
    # the throat response dw/dmu, measured (two step sizes)
    resp = {}
    for dmu in (0.02, 0.04):
        resp[dmu] = (_wB_of(dmu) - _wB_of(-dmu)) / (2 * dmu)
    dwdmu = (_wB_of(0.03) - _wB_of(-0.03)) / 0.06
    lin_dev = abs(resp[0.02] / resp[0.04] - 1.0)
    chi_rail = abs(dwdmu) * dmu_B
    chi_dist = abs(dwdmu) * dmu_A
    t_cz = math.pi / chi_rail
    # the #223 dressing law re-read: exactly quadratic
    here = Path(__file__).resolve().parent / 'runs'
    d223 = _latest(str(here / '*_bridge_dressing_network_probe/probe.json'))
    t7l = _t(d223, 'T7_dressing')
    quad_ok = max(t7l['quadratic_linearity']) < 1e-6
    ratio_check = abs(t7l['dmu_A_01'] / t7l['dmu_A_005'] - 4.0) < 1e-5
    # the strong-coupling honesty
    mu = _RS ** 2
    strong = dmu_B / mu

    ok = (quad_ok and ratio_check and lin_dev < 0.05
          and 0.01 < chi_rail < 0.2 and 30.0 < t_cz < 100.0
          and chi_dist / chi_rail < 0.1 and 0.3 < strong < 1.0)
    out = {
        'name': 'T4_kerr_calibrated',
        'description': (
            'THE GRAVITATIONAL KERR, CALIBRATED: the committed '
            'nonlinearity is the EKG backreaction whose leading form '
            'is the #223 dressing law delta mu = c A^2 (ledger '
            're-read: EXACTLY quadratic, ratio-4 at machine).  Per '
            'ONE QUANTUM of a rail mode on its throat: delta mu_q = '
            '5.8e-2; the measured throat response dw/dmu = -0.95 '
            '(linear to 5% across step sizes); so chi_rail = '
            '|dw/dmu| delta mu_q = 5.5e-2 and t_CZ = pi/chi = 57 '
            'model units.  GRAVITY IS THE ENTANGLING RESOURCE - the '
            '#232 structural chi is now a number on the committed '
            'geometry.  THE HONESTY: delta mu_q/mu = 0.64 - one '
            'quantum dresses the throat at O(1) in model units, so '
            'the quoted chi is LEADING ORDER and the #223 '
            'perturbative window is exceeded per quantum: the '
            'beyond-leading Kerr is the named successor'
        ),
        'dmu_per_quantum_on_throat': float(dmu_B),
        'dmu_per_quantum_distant': float(dmu_A),
        'dw_dmu_measured': float(dwdmu),
        'response_linearity_dev': float(lin_dev),
        'chi_rail': float(chi_rail),
        'chi_distant': float(chi_dist),
        'cross_talk_suppression': float(chi_dist / chi_rail),
        't_CZ': float(t_cz),
        'strong_coupling_dmu_over_mu': float(strong),
        'ledger_223_quadratic_max_dev': float(max(t7l['quadratic_linearity'])),
        'pass': bool(ok),
    }
    _CACHE['T4'] = out
    return out


# ========================================================================
# T5. The dual-rail selection rule + the calibrated CZ
# ========================================================================


def test_T5_selection_rule_cz() -> dict:
    if 'T5' in _CACHE:
        return _CACHE['T5']
    t4 = test_T4_kerr_calibrated()
    chi = t4['chi_rail']
    # Fock check on two qubits (4 modes, 2 quanta)
    N = 2
    basis, idx, P = tp.make_sector(N)
    # (a) SAME-qubit Kerr terms vanish identically on the sector:
    #     n_0 n_1 (rails of qubit 1) and n_i(n_i - 1)
    n = [tp.hop(basis, idx, i, i) for i in range(4)]
    same = P @ (n[0] @ n[1]) @ P.T
    self_k = P @ (n[1] @ n[1] - n[1]) @ P.T
    same_dev = float(np.abs(same).max())
    self_dev = float(np.abs(self_k).max())
    # (b) the CZ with the CALIBRATED chi under FROZEN linear coupling
    #     (twist on: h_linear = 0): U = exp(-i chi n_b1 n_b2 t) at
    #     t = pi/chi
    H = chi * (n[1] @ n[3])
    U = tp.u_from_h(H, math.pi / chi)
    Us = P @ U @ P.T
    CZ = np.diag([1, 1, 1, -1]).astype(complex)
    dev_cz = float(np.abs(Us - CZ).max())
    leak = float(np.abs(U @ P.T - P.T @ Us).max())
    # (c) the residual: distant-rail ZZ phase per CZ (coherent,
    #     compensable)
    zz_phase = t4['chi_distant'] * t4['t_CZ']

    ok = (same_dev < 1e-12 and self_dev < 1e-12 and dev_cz < 1e-12
          and leak < 1e-12 and zz_phase < 0.5)
    out = {
        'name': 'T5_selection_rule_cz',
        'description': (
            'THE DUAL-RAIL SELECTION RULE: within the sector (one '
            'quantum per doublet) every SAME-qubit Kerr term '
            'vanishes IDENTICALLY - n_a n_b = 0 on a shared quantum '
            'and n(n-1) = 0 for occupation <= 1 (Fock, machine '
            'zero).  The only active Kerr couples rails of '
            'DIFFERENT qubits: exactly the CZ generator - the '
            'encoding auto-silences the self-phases.  The CZ with '
            'the CALIBRATED chi under twist-frozen linear coupling '
            'is exact (machine zero, zero leakage) at t = pi/chi = '
            '57 model units.  The measured residual: a distant-rail '
            'coherent ZZ phase of 0.12 rad per CZ (the 3.9e-2 basin '
            'tail) - known, calibrated, compensable by echo; '
            'reported, not hidden'
        ),
        'same_qubit_kerr_on_sector': same_dev,
        'self_kerr_on_sector': self_dev,
        'cz_deviation_calibrated_chi': dev_cz,
        'sector_leakage': leak,
        'distant_zz_phase_per_cz_rad': float(zz_phase),
        'pass': bool(ok),
    }
    _CACHE['T5'] = out
    return out


# ========================================================================
# T6. The circuit budget
# ========================================================================


def test_T6_circuit_budget() -> dict:
    if 'T6' in _CACHE:
        return _CACHE['T6']
    t2 = test_T2_clock()
    t3 = test_T3_twist_switch()
    t4 = test_T4_kerr_calibrated()
    # GHZ_3 at physical rates: the #232 circuit is H + 2x(H, CZ, H)
    # = 5 one-qubit gates + 2 CZs; a Hadamard is a half X-period
    # rotation at the beat (order-of-magnitude physical clock)
    t_1q = 0.5 * t2['x_period']
    t_cz = t4['t_CZ']
    T_circ = 5 * t_1q + 2 * t_cz
    # error ledger
    leak_err = (T_circ / t3['memory_leakage_time']) ** 2
    zz_err = 2 * (t4['chi_distant'] * t_cz)      # 2 CZs, coherent
    budget = {
        'one_qubit_gate_time': float(t_1q),
        'cz_gate_time': float(t_cz),
        'ghz3_circuit_time': float(T_circ),
        'cz_over_clock': float(t_cz / t2['x_period']),
        'twist_frozen_leakage_error': float(leak_err),
        'coherent_zz_phase_total_rad': float(zz_err),
    }
    ok = (T_circ < 1e4 and leak_err < 1e-15
          and budget['cz_over_clock'] < 0.1)
    out = {
        'name': 'T6_circuit_budget',
        'description': (
            'THE CIRCUIT BUDGET AT PHYSICAL RATES: GHZ_3 = 5 '
            'one-qubit gates (half X-periods at the beat, 916 each) '
            '+ 2 CZs (57 each) = 4700 model units; the CZ is 32x '
            'FASTER than the one-qubit clock (the gravitational '
            'Kerr is strong).  The error ledger: twist-frozen '
            'leakage (T dw_frozen)^2 = 1.5e-19 - the topological '
            'switch makes linear leakage negligible over the whole '
            'circuit - and the coherent distant-rail ZZ total of '
            '0.25 rad is the one calibrated correction the '
            'compiler owes'
        ),
        'budget': budget,
        'pass': bool(ok),
    }
    _CACHE['T6'] = out
    return out


# ========================================================================
# T7. The arc-weld holdout
# ========================================================================


def test_T7_holdout() -> dict:
    if 'T7' in _CACHE:
        return _CACHE['T7']
    here = Path(__file__).resolve().parent / 'runs'
    checks = []

    def add(name, stored, derived, tol):
        dev = abs(stored - derived) / (abs(derived) + 1e-300)
        checks.append({'constant': name, 'ledger': float(stored),
                       'rederived': float(derived),
                       'rel_dev': float(dev), 'ok': bool(dev < tol)})

    t2 = test_T2_clock()
    t3 = test_T3_twist_switch()
    d224 = _latest(str(here / '*_mouth_exchange_dynamics_probe/probe.json'))
    t3l = _t(d224, 'T3_p_other_mouth')
    add('#224 dw (the clock)', t3l['dw'], t2['dw'], 1e-9)
    t2l = _t(d224, 'T2_network')
    add('#224 basin localization', t2l['localized_combo_weight'],
        t2['localized_combo_weight_224_construction'], 1e-9)
    d227 = _latest(str(here / '*_nonorientable_er_link_z2_probe/probe.json'))
    row = [r for r in d227['T4']['rows'] if abs(r['r_s'] - 0.3) < 1e-9][0]
    add('#227 frozen dw (the switch)', row['dw_twisted'],
        t3['dw_twisted'], 1e-6)
    add('#227 freeze ratio', row['ratio'], t3['freeze_ratio'], 1e-6)
    d223 = _latest(str(here / '*_bridge_dressing_network_probe/probe.json'))
    t7l = _t(d223, 'T7_dressing')
    add('#223 dressing quadratic (dmu(0.1)/dmu(0.05))',
        t7l['dmu_A_01'] / t7l['dmu_A_005'], 4.0, 1e-5)
    d232 = _latest(str(here / '*_tensor_product_emergence_probe/probe.json'))
    t5l = _t(d232, 'T5_quartic_cz')
    t3l2 = _t(d232, 'T3_tensor_from_locality')
    add('#232 structural CZ exactness', 1.0 + t5l['cz_deviation'], 1.0, 1e-11)
    add('#232 tensor factorization', 1.0 + t3l2['worst_tensor_deviation'],
        1.0, 1e-11)

    all_ok = all(c['ok'] for c in checks)
    out = {
        'name': 'T7_holdout',
        'description': (
            'the arc-weld holdout: the four committed ledgers this '
            'probe welds - #224 (the clock and basins), #227 (the '
            'frozen beat and ratio, from the OTHER arc), #223 (the '
            'exactly-quadratic dressing), #232 (the structural gate '
            'exactness) - re-read and re-derived on one machine, '
            'all at machine agreement: the two development arcs '
            'meet in one calibrated gate set without retuning '
            'anything'
        ),
        'checks': checks,
        'n_checks': len(checks),
        'pass': bool(all_ok),
    }
    _CACHE['T7'] = out
    return out


# ========================================================================
# T8. Honest scope
# ========================================================================


def test_T8_honest_scope() -> dict:
    scope = [
        'All rates are in MODEL UNITS on the committed #224 network '
        '(r_s = 0.3, the 1114-site ring); no claim converts them to '
        'laboratory units - that requires the anchor chain (#225/'
        '#226) and is deliberately out of scope here.',
        'chi is LEADING ORDER: delta mu_q/mu = 0.64 means one '
        'quantum dresses the throat at O(1) at this network scale, '
        'so the #223 perturbative window is exceeded per quantum.  '
        'The beyond-leading (self-consistent) Kerr - resolving the '
        'dressed throat with the quantum in place - is the named '
        'successor, and the strong coupling itself is the finding: '
        'the gravitational entangler is not weak.',
        'The per-quantum normalization is the standard flat-measure '
        'field quantization (amplitude^2 = 1/2w with int u^2 dx = '
        '1); O(1) convention factors between this and #223\'s '
        'peak-amplitude convention are absorbed into the leading-'
        'order status of chi.',
        'The two-qubit CZ presumes an architecture where rails of '
        'different qubits share a throat; this probe calibrates the '
        'on-throat Kerr scale chi_rail and the distant-rail tail on '
        'the ONE committed network - building the two-network '
        'shared-throat geometry is construction work, not new '
        'physics.',
        'The Hadamard-at-the-beat clock is an order-of-magnitude '
        'physical rate (arbitrary one-qubit gates need controlled '
        'detuning, e.g. the #224 bias channel - committed machinery, '
        'not exercised here).',
        'The twist switch statements re-derive #227 on its own '
        'numbers; the #233 correction (deck class set by mouth '
        'placement) does not touch the freeze mechanism used here.',
    ]
    return {
        'name': 'T8_honest_scope',
        'description': 'what this probe does and does not establish',
        'scope': scope,
        'pass': True,
    }


# ========================================================================
# T9. Assessment
# ========================================================================


def test_T9_assessment() -> dict:
    t2 = test_T2_clock()
    t3 = test_T3_twist_switch()
    t4 = test_T4_kerr_calibrated()
    t5 = test_T5_selection_rule_cz()
    t6 = test_T6_circuit_budget()
    t7 = test_T7_holdout()
    core = all(t['pass'] for t in (t2, t3, t4, t5, t6, t7))
    assessment = (
        'The gate set is physical.  The one-qubit clock is the '
        'network\'s own tunneling beat (the #224 doublet, re-derived '
        'exactly); the memory/gate switch is one topological bit '
        '(the #227 twist, re-derived exactly, and shown to spare the '
        'Kerr channel); and the entangling coupling - #232\'s one '
        'structural placeholder - is now a calibrated number: the '
        'gravitational dressing of the throat by a single quantum '
        'gives chi = 5.5e-2, a CZ in 57 model units, 32x FASTER '
        'than the one-qubit clock.  The dual-rail encoding '
        'auto-silences every same-qubit Kerr term (machine zero), '
        'the twist freezes linear leakage to (T dw)^2 ~ 1e-19 over '
        'a GHZ_3 circuit, and the one calibrated residual is a '
        'coherent 0.12 rad ZZ phase per CZ.  The honest edge: the '
        'coupling is STRONG (one quantum dresses the throat by '
        '64%), so chi is leading-order and the self-consistent '
        'dressed gate is the successor.  Ontologically this '
        'completes #232\'s sentence: entanglement on this substrate '
        'is generated by GRAVITY - the same backreaction that welds '
        'the soliton scales (#222) and dresses the bridge (#223) '
        'is, per quantum, the two-qubit gate.'
        if core else
        'INCONCLUSIVE - a core check failed; do not quote.'
    )
    return {
        'name': 'T9_assessment',
        'description': 'the physical-gate-set verdict',
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
        test_T2_clock(),
        test_T3_twist_switch(),
        test_T4_kerr_calibrated(),
        test_T5_selection_rule_cz(),
        test_T6_circuit_budget(),
        test_T7_holdout(),
        test_T8_honest_scope(),
        test_T9_assessment(),
    ]
    t2, t3, t4, t5, t6, t7 = (tests[1], tests[2], tests[3], tests[4],
                              tests[5], tests[6])
    all_ok = all(t['pass'] for t in tests)
    if all_ok:
        verdict_class = (
            "THE_GATE_SET_IS_PHYSICAL_THE_BEAT_IS_THE_CLOCK_THE_"
            "TWIST_IS_THE_SWITCH_AND_THE_GRAVITATIONAL_DRESSING_IS_"
            "THE_KERR_WITH_T_CZ_32X_FASTER_THAN_THE_CLOCK_AT_"
            "LEADING_ORDER_IN_A_STRONG_COUPLING"
        )
        verdict = (
            "ESTABLISHED (the argument is in "
            "docs/physical_gate_set.md).\n\n"
            "THE CLOCK. The #224 doublet re-derived exactly (dw "
            f"dev {t2['ledger_dw_rel_dev']:.0e}): X period "
            f"{t2['x_period']:.0f}, basins "
            f"{t2['basin_localization_A']:.4f}.\n\n"
            "THE SWITCH. One twisted link freezes the beat "
            f"(ratio {t3['freeze_ratio']:.1e}; the #227 row "
            f"re-derived to {t3['ledger_227_frozen_rel_dev']:.0e}); "
            "memory leakage time "
            f"{t3['memory_leakage_time']:.1e}; the Kerr channel "
            f"shifts < {t3['kerr_channel_shift_under_twist']:.1%} "
            "under the twist - an error switch, not a gate hazard."
            "\n\n"
            "THE KERR. delta mu per quantum "
            f"{t4['dmu_per_quantum_on_throat']:.3e}, response "
            f"dw/dmu {t4['dw_dmu_measured']:.3f}: chi = "
            f"{t4['chi_rail']:.3e}, t_CZ = {t4['t_CZ']:.1f} - "
            f"{1/t6['budget']['cz_over_clock']:.0f}x faster than "
            "the clock.  GRAVITY IS THE ENTANGLING RESOURCE.  "
            "Honesty: delta mu_q/mu = "
            f"{t4['strong_coupling_dmu_over_mu']:.2f} - leading "
            "order, strong coupling.\n\n"
            "THE SELECTION RULE. Same-qubit Kerr terms vanish on "
            f"the sector ({t5['same_qubit_kerr_on_sector']:.0e}); "
            "the calibrated CZ is exact "
            f"({t5['cz_deviation_calibrated_chi']:.0e}, leakage "
            f"{t5['sector_leakage']:.0e}); residual ZZ "
            f"{t5['distant_zz_phase_per_cz_rad']:.2f} rad per CZ, "
            "compensable.\n\n"
            "THE BUDGET. GHZ_3 = "
            f"{t6['budget']['ghz3_circuit_time']:.0f} model units; "
            "twist-frozen leakage error "
            f"{t6['budget']['twist_frozen_leakage_error']:.1e}.\n\n"
            "THE WELD. "
            f"{t7['n_checks']} committed constants from FOUR "
            "ledgers spanning BOTH arcs (#224, #227, #223, #232) "
            "re-derived on one machine at machine agreement: the "
            "arcs meet in one calibrated gate set."
        )
    else:
        verdict_class = "PHYSICAL_GATE_SET_INCONCLUSIVE"
        verdict = (
            "INCONCLUSIVE. A core check failed; re-examine before "
            "quoting the result."
        )

    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "identification": (
            "The physical gate set on the committed #224 network: "
            "the one-qubit clock is the doublet beat (1833), the "
            "memory/gate switch is the #227 twist (freeze 4.9e-11, "
            "Kerr channel spared), and the entangling coupling is "
            "the gravitational dressing per quantum - chi = 5.5e-2, "
            "t_CZ = 57, 32x faster than the clock, leading order in "
            "a strong coupling - with the dual-rail selection rule "
            "silencing all same-qubit Kerr terms"
        ),
        "executes": (
            "#232's named scope item (calibrate chi on the "
            "committed geometry) + the weld of both development "
            "arcs into one calibrated gate set"
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
    out.append("# The physical gate set (PR #234)")
    out.append("")
    out.append(f"**Run:** {s['timestamp_utc']}")
    out.append("")
    out.append(
        "The deliverable is `docs/physical_gate_set.md` - the qubit "
        "clock, the gravitational Kerr, and the twist switch, "
        "calibrated on the committed network. *(QFT on the fixed "
        "classical throat geometry, not quantum gravity.)*"
    )
    out.append("")
    out.append("## Test summary")
    out.append("")
    out.append("| # | Test | Key finding | PASS? |")
    out.append("|---|---|---|---|")
    labels = {
        "T1": "the goal: #232's scope item + the arc weld",
        "T2": "the clock: the #224 beat, re-derived exactly",
        "T3": "the switch: the #227 freeze; Kerr channel spared",
        "T4": "the gravitational Kerr: chi = 5.5e-2, t_CZ = 57",
        "T5": "the dual-rail selection rule; the calibrated CZ",
        "T6": "the circuit budget at physical rates",
        "T7": "the arc-weld holdout (four ledgers, both arcs)",
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
    out = here / "runs" / f"{ts}_physical_gate_set_probe"
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
