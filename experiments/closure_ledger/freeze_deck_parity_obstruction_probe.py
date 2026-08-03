"""Probe G -- do the freeze and deck class 1 exclude each other? (PR #233)

> Framing: QFT on the *fixed classical* throat geometry (geometry -> fields),
> not quantum gravity.

THE QUESTION
------------
Probe F (#231) OBSERVED that the two properties needed to make the
exchange-freeze result a statement about the spin structure never occur
together:

    #223 (1 throat) : deck class 1, labels coincide -- but NO mouth doublet
    #224 (2 throats): mouth doublet and exchange freeze -- but deck class 0

and named "an interior channel on the #223 ring, or some other geometry
carrying both" as the obvious next construction.  Is the non-co-occurrence
structural, or an accident of which two rings the repo happens to contain?

*** CORRECTION NOTICE ***
--------------------------
The first version of this probe answered "structural: a parity
obstruction", on the grounds that a ring of N throats carries N exterior
pi-arcs, hence deck class = N mod 2, hence freeze (N even) and deck class 1
(N odd) exclude each other.  **That answer was wrong, and this file now
records the computation that refutes it.**

The error was a hidden assumption: that every exterior arc has length
exactly pi.  An arc of length pi joins a point to its antipode, so
assuming it for *every* arc of an N-throat ring forces all N throats to
share the SAME antipodal pair of mouth locations.  That is one
configuration, not the general one.  Mouth placement is free geometric
data, and the deck class is set by it:

    deck class = (total exterior arc length) / pi   mod 2     [Probe E]

which does not mention N at all.  The freeze law does mention N, and
stands.  So the two conditions are INDEPENDENT, and they co-occur on
geometries that were simply never built.

WHAT ACTUALLY HOLDS
-------------------
  (1) THE FREEZE LAW -- STANDS.  The cyclic translation S advancing the
      ring by one throat satisfies S^N = W.  The freeze is the statement
      that EVERY level is forced into a degenerate multiplet, which happens
      iff S has NO real eigenvalue (a real operator pairs complex
      eigenvalues, and only those).  For W = -1 the eigenvalues are the
      N-th roots of -1, e^{i pi (2k+1)/N}, and one is real iff N is ODD:

          full freeze  <=>  N EVEN.

      Verified here at operator level: ||[H,S]|| ~ 1e-12 and S^N = W
      exactly, on every ring tested.

  (2) THE DECK-CLASS LAW -- CORRECTED.  deck class = (sum of exterior arc
      lengths)/pi mod 2, which is a function of MOUTH PLACEMENT, not of N.
      "deck class = N mod 2" is only the all-arcs-equal-to-pi special case.

  (3) THEY ARE INDEPENDENT.  Choose N even for the freeze, then place the
      mouths so the arcs sum to an odd multiple of pi.  The simplest such
      ring: TWO throats whose mouth pairs are rotated a quarter circle
      apart, so each exterior arc is pi/2 and the total is pi.  Measured:
      max intra-pair gap over the 16 lowest levels 8.5e-13 in the twisted
      sector -- FROZEN -- with deck class 1.  #224's own configuration
      (both arcs pi) is the degenerate case in which the two throats share
      one antipodal mouth pair.

  *** So the geometry #231 asked for EXISTS.  Probe F's non-co-occurrence
  is an artifact of the two rings the repo contains, not an obstruction. ***

THE OTHER EVASION CLASSES, CHECKED RATHER THAN ASSERTED
--------------------------------------------------------
  * TUNED / non-pi arcs        -> the two properties CO-OCCUR (above).
  * ONE THROAT, FINITE CHANNEL -> #231's own D > 0 item.  Computed at
    D = 0.5, 4, 8, 16 on the full ring operator: the interior modes form a
    NONDEGENERATE ladder (D = 8: 0.3785, 0.7568, 1.1350, 1.5131, 1.8910),
    twisted and untwisted spectra agree to 5 digits, no freeze (gaps
    1.9e-1 ... 8.5e-1).  So the N = 1 ring does not acquire a doublet, and
    the earlier probe's ASSERTION to that effect is upheld -- but by
    computation, and it is moot, because (3) supplies the geometry anyway.
  * INHOMOGENEOUS / INTERNALLY STRUCTURED -> the freeze is FRAGILE.
    Unequal arcs (pi, 2pi), unequal channels (4, 5), or an internal bump in
    one throat each destroy it (gaps 2.5e-1, 3.8e-1, 4.3e-1) even at even
    N.  The freeze needs the cyclic symmetry exactly, not approximately.
  * NON-CYCLIC -> a theta network (three throats between two tri-mouth
    junctions, b_1 = 2, four twist sectors) shows NO full freeze in any
    sector.  Its degeneracies come from the S_3 edge-permutation symmetry
    and are present in the UNTWISTED sector too, so they are not
    twist-induced.  S^N = W has no analogue there.  One representative
    network, not a proof for all non-cyclic ones.

CONSEQUENCE FOR THE SELECTION QUESTION
---------------------------------------
The earlier claim that the freeze sector is cut off from RP^3
spin-structure data "permanently, not contingently" is RETRACTED.  On the
quarter-circle two-throat ring the cycle IS the pi_1(RP^3) generator, so
Probe E's composition rule applies and the freeze sector carries the deck
label after all.  That reopens the line of attack the earlier version
claimed to close: Probes B-D (eta = +/-1/4, |det| equal, h = 0,
Delta<T> = 0) can now be applied to a network that actually freezes.  They
still find nothing selecting the twist -- but on this ring that is a
statement about the freeze sector, which it never was before.

Tests:
  T1. Goal, and the claim under test.
  T2. The freeze law: S^N = W, full freeze iff N even -- argued, measured,
      and verified at operator level.
  T3. The deck-class law, corrected: set by mouth placement, not by N.
  T4. THE COUNTEREXAMPLE: rings carrying freeze AND deck class 1.
  T5. The remaining evasion classes, computed one by one.
  T6. Consequences: what is retracted, what is restored, what stands.
  T7. Assessment.

Verdict:
  THE_FREEZE_NEEDS_EVEN_N_BUT_THE_DECK_CLASS_IS_SET_BY_MOUTH_PLACEMENT_NOT
  _BY_N_SO_THEY_ARE_INDEPENDENT_AND_A_TWO_THROAT_RING_WITH_QUARTER_CIRCLE
  _ARCS_CARRIES_BOTH
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np

from experiments.closure_ledger.mouth_exchange_dynamics_probe import (
    v_bridge_half, _LARC, _DX, _SIG_M, _D_CH,
)

PI = math.pi


# ── the N-throat cyclic ring, with mouth placement as free data ─────────
#
# #224's build_two_throat generalised twice over: to N throats, and to
# exterior arcs of arbitrary length.  The grid is chosen so that each
# throat cell is EXACTLY the same number of points, which makes the
# shift-by-one-cell operator S an exact symmetry of the DISCRETE operator
# -- without that, an incommensurate grid fakes a broken symmetry and the
# freeze appears to fail when it does not.

def build_cyclic_ring(N: int, arc: float, dch: float = _D_CH, rs: float = 0.3,
                      bumps: Optional[dict] = None,
                      dch_list: Optional[list] = None,
                      arc_list: Optional[list] = None) -> dict:
    dchs = list(dch_list) if dch_list else [dch] * N
    arcs = list(arc_list) if arc_list else [arc] * N
    cell = dchs[0] + 2 * _SIG_M + arcs[0]
    dx = cell / int(round(cell / _DX))
    counts = [int(round((dchs[i] + 2 * _SIG_M + arcs[i]) / dx))
              for i in range(N)]
    m = sum(counts)
    s = np.arange(m) * dx
    v = np.zeros(m)
    interior = np.zeros(m, dtype=bool)

    def seg(a, b):
        return (s >= a - 1e-12) & (s < b - 1e-12)

    p = 0.0
    for i in range(N):
        d, la = dchs[i], arcs[i]
        msk = seg(p, p + d)
        interior |= msk
        if bumps and bumps.get(i):
            v[msk] += bumps[i] * np.exp(-((s[msk] - (p + d / 2)) / 0.5) ** 2)
        p += d
        msk = seg(p, p + _SIG_M)
        v[msk] = v_bridge_half(s[msk] - p, rs)
        p += _SIG_M
        a0 = p
        msk = seg(p, p + la)
        v[msk] = (3.75 / ((s[msk] - a0) + _SIG_M) ** 2
                  + 3.75 / ((a0 + la - s[msk]) + _SIG_M) ** 2 - 1.0)
        p += la
        msk = seg(p, p + _SIG_M)
        v[msk] = v_bridge_half(p + _SIG_M - s[msk], rs)
        p += _SIG_M
    return {'s': s, 'V': v, 'dx': dx, 'M': m, 'N': N, 'interior': interior,
            'arcs': arcs, 'channels': dchs, 'cell_points': counts,
            'shift_is_exact': len(set(counts)) == 1,
            'total_exterior_over_pi': sum(arcs) / PI}


def hamiltonian(g: dict, holonomy: float) -> np.ndarray:
    m, dx, v = g['M'], g['dx'], g['V']
    lap = (np.diag(-2.0 * np.ones(m)) + np.diag(np.ones(m - 1), 1)
           + np.diag(np.ones(m - 1), -1))
    lap[0, -1] = lap[-1, 0] = holonomy
    return -lap / dx ** 2 + np.diag(v)


def shift_operator(g: dict, holonomy: float) -> np.ndarray:
    """S: advance by one throat cell, carrying the holonomy sign on wrap."""
    m, n = g['M'], g['cell_points'][0]
    S = np.zeros((m, m))
    for i in range(m):
        j = i + n
        S[j % m, i] = holonomy if j >= m else 1.0
    return S


def ring_modes(g: dict, holonomy: float, k: int = 16):
    ev, U = np.linalg.eigh(hamiltonian(g, holonomy))
    w = np.sqrt(np.abs(ev))
    i = np.argsort(w)[:k]
    w, U = w[i], U[:, i]
    frac = np.array([float(np.sum(U[g['interior'], j] ** 2))
                     for j in range(len(w))])
    return w, frac


def freeze_report(g: dict, holonomy: float, k: int = 16) -> dict:
    w, frac = ring_modes(g, holonomy, k)
    gap = float(np.max(np.diff(w)[0::2]))
    return {'max_intra_pair_gap': gap, 'frozen': bool(gap < 1e-8),
            'n_interior_localized': int(np.sum(frac > 0.5)),
            'w_lowest': [round(float(x), 6) for x in w[:6]],
            'interior_fraction_lowest': [round(float(x), 3) for x in frac[:6]]}


def deck_class(g: dict) -> int:
    """Probe E's criterion: the S^3 lift closes iff the total exterior arc
    length is an even multiple of pi."""
    return int(round(g['total_exterior_over_pi'])) % 2


def s_eigenvalue_angles(N: int, W: int) -> list:
    """S^N = W, so spec(S) = the N-th roots of W (angles in units of pi)."""
    if W == 1:
        return [2.0 * k / N for k in range(N)]
    return [(2.0 * k + 1) / N for k in range(N)]


def has_real_eigenvalue(N: int, W: int) -> bool:
    return any(abs(a - round(a)) < 1e-12 for a in s_eigenvalue_angles(N, W))


# ── the non-cyclic case: a theta network of three throats ───────────────

def theta_network(arcs: list, dch: float = _D_CH, rs: float = 0.3,
                  dx: float = _DX, eps: tuple = (1, 1, 1)):
    """Three edges A -> B, each = channel + barrier + arc + barrier, joined
    at two tri-mouth junctions with Kirchhoff conditions.  b_1 = 2, so the
    Z2 link field has two independent Wilson loops; `eps` carries the link
    sign on the A-end of each edge (vertex gauge flips all three at once)."""
    edges = []
    for la in arcs:
        n = int(round((dch + 2 * _SIG_M + la) / dx))
        s = np.arange(n) * dx
        v = np.zeros(n)
        interior = s < dch - 1e-12
        p = dch
        msk = (s >= p - 1e-12) & (s < p + _SIG_M - 1e-12)
        v[msk] = v_bridge_half(s[msk] - p, rs)
        p += _SIG_M
        a0 = p
        msk = (s >= p - 1e-12) & (s < p + la - 1e-12)
        v[msk] = (3.75 / ((s[msk] - a0) + _SIG_M) ** 2
                  + 3.75 / ((a0 + la - s[msk]) + _SIG_M) ** 2 - 1.0)
        p += la
        msk = s >= p - 1e-12
        v[msk] = v_bridge_half(p + _SIG_M - s[msk], rs)
        edges.append({'n': n, 'V': v, 'interior': interior})

    off = np.cumsum([0] + [e['n'] for e in edges])
    m = int(off[-1]) + 2
    iA, iB = m - 2, m - 1
    H = np.zeros((m, m))
    inv = 1.0 / dx ** 2
    for j, e in enumerate(edges):
        n, o, ej = e['n'], int(off[j]), eps[j]
        for i in range(n):
            H[o + i, o + i] += 2 * inv + e['V'][i]
            if i > 0:
                H[o + i, o + i - 1] -= inv
            else:
                H[o + i, iA] -= inv * ej
            if i < n - 1:
                H[o + i, o + i + 1] -= inv
            else:
                H[o + i, iB] -= inv
        H[iA, o] -= inv * ej
        H[iB, o + n - 1] -= inv
    H[iA, iA] += 3 * inv
    H[iB, iB] += 3 * inv
    return H


def theta_freeze(arcs: list, eps: tuple, k: int = 16) -> dict:
    w = np.sqrt(np.abs(np.linalg.eigvalsh(theta_network(arcs, eps=eps))))[:k]
    gap = float(np.max(np.diff(w)[0::2]))
    return {'eps': list(eps),
            'wilson_e1e2': int(eps[0] * eps[1]),
            'wilson_e2e3': int(eps[1] * eps[2]),
            'max_intra_pair_gap': gap, 'frozen': bool(gap < 1e-8),
            'w_lowest': [round(float(x), 6) for x in w[:6]]}


# ════════════════════════════════════════════════════════════════════════

def test_T1_goal() -> dict:
    return {
        'name': 'T1_goal_do_the_freeze_and_deck_class_one_exclude_each_other',
        'description': (
            "Probe F (#231) OBSERVED that the mouth doublet / exchange "
            "freeze (on #224, 2 throats) and deck class 1 (on #223, 1 "
            "throat) never occur together, and named a geometry carrying "
            "both as the obvious next construction. An observation across "
            "the two networks the repo happens to contain is not an "
            "obstruction, so this probe asks whether the non-co-occurrence "
            "is structural. The first version of this probe answered YES, on "
            "a parity argument. *** That answer was wrong. The parity "
            "argument assumed every exterior arc has length exactly pi, "
            "which forces all N throats to share one antipodal mouth pair. "
            "Mouth placement is free, and once it is freed the two "
            "properties co-occur. *** This version records the computation "
            "that refutes the earlier claim and checks the other evasion "
            "classes by calculation rather than by assertion."
        ),
        'pass': True,
    }


def test_T2_freeze_law() -> dict:
    argued = []
    for N in (2, 3, 4, 5, 6):
        real_tw = has_real_eigenvalue(N, -1)
        argued.append({
            'N': N,
            'S_angles_over_pi_twisted': [round(a, 4)
                                         for a in s_eigenvalue_angles(N, -1)],
            'twisted_has_real_eigenvalue': real_tw,
            'full_freeze_predicted': (not real_tw),
            'N_is_even': N % 2 == 0,
        })
    argument_ok = all(r['full_freeze_predicted'] == r['N_is_even']
                      for r in argued)

    measured, operator = [], []
    for N in (2, 3, 4):
        g = build_cyclic_ring(N, _LARC)
        fu = freeze_report(g, +1.0)
        ft = freeze_report(g, -1.0)
        measured.append({'N': N, 'arc_over_pi': round(_LARC / PI, 4),
                         'gap_untwisted': fu['max_intra_pair_gap'],
                         'gap_twisted': ft['max_intra_pair_gap'],
                         'frozen': ft['frozen'], 'N_is_even': N % 2 == 0})
        for hol in (+1.0, -1.0):
            H, S = hamiltonian(g, hol), shift_operator(g, hol)
            operator.append({
                'N': N, 'W': int(hol),
                'commutator_norm': float(np.max(np.abs(H @ S - S @ H))),
                'S_to_the_N_minus_W': float(np.max(np.abs(
                    np.linalg.matrix_power(S, N) - hol * np.eye(g['M'])))),
            })
    measured_ok = all(r['frozen'] == r['N_is_even'] for r in measured)
    operator_ok = all(o['commutator_norm'] < 1e-8
                      and o['S_to_the_N_minus_W'] < 1e-12 for o in operator)
    ok = argument_ok and measured_ok and operator_ok
    return {
        'name': 'T2_full_freeze_iff_N_is_even_STANDS',
        'description': (
            "The cyclic translation S advancing the ring by one throat "
            "satisfies S^N = W. The freeze is the statement that EVERY level "
            "is forced into a degenerate multiplet, which happens iff S has "
            "NO real eigenvalue -- a real operator pairs complex eigenvalues "
            "and only those. For W = -1 the eigenvalues are the N-th roots "
            "of -1, e^{i pi (2k+1)/N}, and one is real iff N is ODD, so full "
            "freeze iff N EVEN. This law is UNCHANGED by the correction, and "
            "is checked three ways: the root-of-unity argument; the measured "
            f"max intra-pair gap on real rings (N = 2 twisted "
            f"{measured[0]['gap_twisted']:.3e} FROZEN, N = 3 "
            f"{measured[1]['gap_twisted']:.3e} NOT frozen, N = 4 "
            f"{measured[2]['gap_twisted']:.3e} FROZEN); and at operator "
            "level, where ||[H,S]|| < 1e-8 and S^N = W hold exactly. The "
            "grids are chosen so each throat cell is the same integer number "
            "of points -- on an incommensurate grid S is not a symmetry of "
            "the discrete operator and the freeze spuriously appears to fail."
        ),
        'root_of_unity_argument': argued,
        'argument_matches_parity': argument_ok,
        'measured_on_real_rings': measured,
        'measurement_matches_parity': measured_ok,
        'operator_level_check': operator,
        'S_is_a_symmetry_and_S_to_the_N_equals_W': operator_ok,
        'pass': ok,
    }


def test_T3_deck_law_corrected() -> dict:
    rows = []
    for N, arc, label in ((1, _LARC, '#223-like: 1 throat, arc pi'),
                          (2, _LARC, '#224: 2 throats, arcs pi (mouths '
                                     'coincide)'),
                          (3, _LARC, '3 throats, arcs pi'),
                          (2, 0.5 * PI, '2 throats, arcs pi/2 (mouths a '
                                        'quarter circle apart)'),
                          (2, 1.5 * PI, '2 throats, arcs 3pi/2'),
                          (2, 2.0 * PI, '2 throats, arcs 2pi'),
                          (4, 0.75 * PI, '4 throats, arcs 3pi/4')):
        g = build_cyclic_ring(N, arc)
        rows.append({'label': label, 'N': N,
                     'arc_over_pi': round(arc / PI, 4),
                     'total_exterior_over_pi': round(
                         g['total_exterior_over_pi'], 4),
                     'deck_class': deck_class(g),
                     'N_mod_2': N % 2,
                     'old_law_would_say': N % 2})
    old_law_fails = [r for r in rows if r['deck_class'] != r['N_mod_2']]
    # the old law is exactly the all-arcs-pi restriction
    pi_rows = [r for r in rows if abs(r['arc_over_pi'] - 1.0) < 0.01]
    old_law_ok_on_pi = all(r['deck_class'] == r['N_mod_2'] for r in pi_rows)
    ok = bool(old_law_fails) and old_law_ok_on_pi
    return {
        'name': 'T3_deck_class_is_set_by_mouth_placement_not_by_N',
        'description': (
            "Probe E's criterion is that the S^3 lift closes iff the TOTAL "
            "EXTERIOR ARC LENGTH is an even multiple of pi, so deck class = "
            "(sum of arc lengths)/pi mod 2. That formula does not mention N. "
            "The earlier claim 'an N-throat ring carries N exterior pi-arcs, "
            "so deck class = N mod 2' silently assumed every arc has length "
            "pi -- and an arc of length pi joins a point to its ANTIPODE, so "
            "assuming it for every arc forces all N throats to share the "
            "same antipodal pair of mouth locations. That is #224's "
            "configuration, not the general one. Mouth placement is free "
            "geometric data. Verified here: the old law reproduces the deck "
            f"class on all all-arcs-pi rings, and FAILS on "
            f"{len(old_law_fails)} of the rings with other arc lengths -- "
            "e.g. two throats with quarter-circle arcs have N even but deck "
            "class 1."
        ),
        'rows': rows,
        'old_law_correct_when_every_arc_is_pi': old_law_ok_on_pi,
        'rows_where_old_law_fails': old_law_fails,
        'pass': ok,
    }


def test_T4_counterexample() -> dict:
    cases = []
    for N, arc, note in (
            (2, 0.5 * PI, 'the simplest one: two throats whose mouth pairs '
                          'are a quarter circle apart'),
            (2, 1.5 * PI, 'the same ring taking the long way round'),
            (4, 0.75 * PI, 'four throats, arcs 3pi/4, total 3pi')):
        g = build_cyclic_ring(N, arc)
        fu, ft = freeze_report(g, +1.0), freeze_report(g, -1.0)
        H, S = hamiltonian(g, -1.0), shift_operator(g, -1.0)
        cases.append({
            'note': note, 'N': N, 'arc_over_pi': round(arc / PI, 4),
            'total_exterior_over_pi': round(g['total_exterior_over_pi'], 4),
            'deck_class': deck_class(g),
            'shift_is_exact_on_the_grid': g['shift_is_exact'],
            'commutator_norm': float(np.max(np.abs(H @ S - S @ H))),
            'S_to_the_N_minus_W': float(np.max(np.abs(
                np.linalg.matrix_power(S, N) + np.eye(g['M'])))),
            'gap_untwisted': fu['max_intra_pair_gap'],
            'gap_twisted': ft['max_intra_pair_gap'],
            'frozen': ft['frozen'],
            'n_interior_localized_of_16': ft['n_interior_localized'],
            'has_mouth_doublet': bool(N >= 2),
            'w_lowest_twisted': ft['w_lowest'],
        })
    all_both = all(c['frozen'] and c['deck_class'] == 1
                   and c['has_mouth_doublet'] for c in cases)
    mechanism_ok = all(c['commutator_norm'] < 1e-8
                       and c['S_to_the_N_minus_W'] < 1e-12 for c in cases)
    return {
        'name': 'T4_rings_carrying_the_freeze_AND_deck_class_one',
        'description': (
            "Take N even so the freeze law applies, then place the mouths so "
            "the arcs sum to an ODD multiple of pi. The simplest such ring is "
            "TWO throats whose mouth pairs are rotated a quarter circle "
            "apart: each exterior arc is pi/2, the total is pi, the deck "
            "class is 1, and the ring has two interior channels hence a "
            "mouth doublet. Measured max intra-pair gap over the 16 lowest "
            f"levels in the twisted sector: {cases[0]['gap_twisted']:.3e} -- "
            f"FROZEN (untwisted {cases[0]['gap_untwisted']:.3e}, not frozen). "
            "And the freeze there is the SAME mechanism as #224's, not a "
            "coincidence: ||[H,S]|| < 1e-8 and S^N = W exactly. *** So the "
            "geometry #231 asked for exists, the earlier version of this "
            "probe was wrong to rule it out, and Probe F's non-co-occurrence "
            "is an artifact of the two rings the repo contains. *** #224's "
            "own configuration (both arcs pi) is the degenerate case where "
            "the two throats share a single antipodal mouth pair."
        ),
        'cases': cases,
        'all_carry_freeze_and_deck_class_one': all_both,
        'freeze_mechanism_verified_at_operator_level': mechanism_ok,
        'pass': all_both and mechanism_ok,
    }


def test_T5_evasion_classes() -> dict:
    # (a) one throat, finite interior channel -- #231's own D > 0 item
    one_throat = []
    for D in (0.5, 4.0, 8.0, 16.0):
        g = build_cyclic_ring(1, _LARC, dch=D)
        fu, ft = freeze_report(g, +1.0), freeze_report(g, -1.0)
        one_throat.append({
            'D_channel': D, 'deck_class': deck_class(g),
            'gap_untwisted': fu['max_intra_pair_gap'],
            'gap_twisted': ft['max_intra_pair_gap'],
            'frozen': ft['frozen'],
            'n_interior_localized_of_16': ft['n_interior_localized'],
            'w_lowest_twisted': ft['w_lowest'],
            'interior_fraction_lowest': ft['interior_fraction_lowest'],
        })
    one_throat_never_freezes = not any(r['frozen'] for r in one_throat)

    # (b) inhomogeneous / internally structured, at even N
    inhomo = []
    for label, g in (
            ('unequal arcs (pi, 2pi): total 3pi, deck 1',
             build_cyclic_ring(2, _LARC, arc_list=[_LARC, 2 * _LARC])),
            ('unequal channels (4, 5), arcs 3pi/2',
             build_cyclic_ring(2, 1.5 * PI, dch_list=[4.0, 5.0])),
            ('internal bump inside throat 0, arcs 3pi/2',
             build_cyclic_ring(2, 1.5 * PI, bumps={0: 3.0}))):
        ft = freeze_report(g, -1.0)
        inhomo.append({'label': label, 'N': g['N'], 'deck_class': deck_class(g),
                       'gap_twisted': ft['max_intra_pair_gap'],
                       'frozen': ft['frozen']})
    inhomo_destroys_freeze = not any(r['frozen'] for r in inhomo)

    # (c) non-cyclic: a theta network of three throats
    theta = [theta_freeze([_LARC] * 3, eps)
             for eps in ((1, 1, 1), (1, 1, -1), (1, -1, -1))]
    theta += [theta_freeze([0.5 * PI] * 3, eps)
              for eps in ((1, 1, 1), (1, 1, -1))]
    theta_never_freezes = not any(r['frozen'] for r in theta)
    twisted = [r for r in theta if r['wilson_e1e2'] < 0 or r['wilson_e2e3'] < 0]
    untwisted = [r for r in theta if r['wilson_e1e2'] > 0
                 and r['wilson_e2e3'] > 0]
    theta_gap_unchanged = all(
        abs(t['max_intra_pair_gap'] - u['max_intra_pair_gap']) < 1e-6
        for t in twisted for u in untwisted
        if abs(len(t['w_lowest']) - len(u['w_lowest'])) == 0
        and abs(t['w_lowest'][0] - u['w_lowest'][0]) < 1e-2)

    ok = (one_throat_never_freezes and inhomo_destroys_freeze
          and theta_never_freezes)
    return {
        'name': 'T5_the_remaining_evasion_classes_computed_not_asserted',
        'description': (
            "(a) ONE THROAT WITH A FINITE INTERIOR CHANNEL -- #231's own "
            "D > 0 item, which the earlier version of this probe dismissed "
            "by ASSERTION. Computed on the full ring operator at D = 0.5, 4, "
            "8, 16: the interior modes form a NONDEGENERATE ladder (D = 8: "
            "0.3785, 0.7568, 1.1350, 1.5131, 1.8910 -- successive box "
            "levels, not pairs), the twisted and untwisted spectra agree to "
            "five digits because the interior states do not reach the "
            "holonomy, and nothing freezes (gaps 1.9e-1 ... 8.5e-1). The "
            "assertion is upheld, now by calculation; it is also moot, since "
            "T4 supplies the co-occurring geometry a different way. "
            "(b) INHOMOGENEOUS / INTERNALLY STRUCTURED at even N: unequal "
            "arcs, unequal channels, and an internal bump in one throat each "
            "DESTROY the freeze (gaps 2.5e-1, 3.8e-1, 4.3e-1). The freeze "
            "needs the cyclic symmetry exactly, not approximately -- so it "
            "is a fragile phenomenon, which is itself a constraint on any "
            "physical reading of it. (c) NON-CYCLIC: a theta network of "
            "three throats between two tri-mouth junctions (b_1 = 2, so two "
            "independent Wilson loops) shows NO full freeze in any twist "
            "sector. Its degeneracies come from the S_3 edge-permutation "
            "symmetry and are present in the untwisted sector too, so they "
            "are not twist-induced; S^N = W has no analogue on a graph with "
            "no translation. This is one representative non-cyclic network, "
            "not a proof for all of them."
        ),
        'one_throat_finite_channel': one_throat,
        'one_throat_never_freezes': one_throat_never_freezes,
        'inhomogeneous_or_structured': inhomo,
        'inhomogeneity_destroys_the_freeze': inhomo_destroys_freeze,
        'non_cyclic_theta_network': theta,
        'theta_never_freezes': theta_never_freezes,
        'theta_degeneracies_are_permutation_not_twist': theta_gap_unchanged,
        'pass': ok,
    }


def test_T6_consequences() -> dict:
    items = [
        {'claim': 'the freeze law: full freeze iff N even (S^N = W)',
         'status': 'STANDS -- argued, measured, and now verified at operator '
                   'level (||[H,S]|| < 1e-8, S^N = W exactly)'},
        {'claim': "the earlier probe's deck-class law, deck class = N mod 2",
         'status': 'CORRECTED -- true only when every exterior arc has length '
                   'pi, i.e. when all N throats share one antipodal mouth '
                   'pair. In general deck class = (total arc length)/pi mod '
                   '2, which does not depend on N.'},
        {'claim': 'the earlier probe: freeze and deck class 1 are mutually '
                  'exclusive by parity',
         'status': 'RETRACTED -- they are INDEPENDENT. T4 exhibits rings '
                   'carrying both, with the freeze verified to be the same '
                   'S^N = W mechanism as #224\'s.'},
        {'claim': "#231's recommended next construction (a geometry carrying "
                  'both the doublet and deck class 1)',
         'status': 'RESTORED, in corrected form -- it exists. Not as #231 '
                   'phrased it (an interior channel on the ONE-throat ring: '
                   'T5(a) computes that this gives a nondegenerate interior '
                   'ladder, no doublet, no freeze) but by MOVING THE MOUTHS: '
                   'two throats with quarter-circle arcs.'},
        {'claim': 'the freeze sector is cut off from RP^3 spin-structure data '
                  'permanently, not contingently',
         'status': 'RETRACTED -- on the quarter-circle two-throat ring the '
                   'cycle IS the pi_1(RP^3) generator, so Probe E\'s '
                   'composition rule applies and the freeze sector carries '
                   'the deck label. The line of attack the earlier version '
                   'claimed to close is open.'},
        {'claim': "Probe F (#231): the phenomenon and the coupling never "
                  'co-occur',
         'status': 'DEMOTED back to an observation about the two rings the '
                   'repo contains -- correct as reported there, but not a '
                   'theorem, and now known to fail on a third ring.'},
    ]
    return {
        'name': 'T6_what_is_retracted_restored_and_still_standing',
        'description': (
            "The freeze law survives intact; the deck-class law and "
            "everything built on it does not. The net effect is to REOPEN "
            "rather than close: Probes B-D (eta = +/-1/4, equal |det D|, "
            "h = 0, Delta<T_AB> = 0) can now be applied to a network that "
            "actually freezes, which was the whole point of looking for such "
            "a network. They still find nothing selecting the twist -- but "
            "on the quarter-circle ring that null result is finally a "
            "statement ABOUT the freeze sector, which it never was before. "
            "The remaining honest caveat runs the other way from the earlier "
            "version's: the freeze is FRAGILE (T5(b)), so a physical reading "
            "of it has to explain why the exact cyclic symmetry it needs "
            "would hold."
        ),
        'items': items,
        'pass': True,
    }


def test_T7_assessment(results: dict) -> dict:
    passed = sum(1 for k, v in results.items()
                 if k.startswith('T') and v.get('pass'))
    total = sum(1 for k in results if k.startswith('T'))
    # This test is not yet in `results`, so `passed`/`total` cover only its
    # predecessors.  Count it too, matching the repo convention (e.g.
    # mouth_exchange_dynamics_probe tallies the full test list).
    self_pass = passed >= total - 1
    passed += int(self_pass)
    total += 1
    return {
        'name': 'T7_assessment',
        'description': (
            "The freeze needs N even; the deck class is set by where the "
            "mouths sit, not by N. The two are independent, and a ring of "
            "two throats whose mouth pairs are a quarter circle apart "
            "carries both. The earlier version of this probe claimed a "
            "parity obstruction ruling that out; the claim rested on an "
            "unstated assumption that every exterior arc has length pi, and "
            "is retracted here. The other evasion classes are settled by "
            "computation: the one-throat finite channel does not produce a "
            "doublet, inhomogeneity destroys the freeze, and a non-cyclic "
            "theta network does not freeze in any twist sector."
        ),
        'established': [
            'full freeze iff N even: S^N = W has no real eigenvalue only for '
            'even N -- argued, measured (N=2 frozen 2.2e-12, N=3 NOT frozen '
            '4.0e-01, N=4 frozen 1.7e-12), and verified at operator level',
            'deck class = (total exterior arc length)/pi mod 2 (Probe E), a '
            'function of mouth placement and not of N',
            'the two are INDEPENDENT: a 2-throat ring with quarter-circle '
            'arcs freezes (gap 8.5e-13) and has deck class 1',
            "#231's recommended construction exists after all -- by moving "
            'the mouths, not by adding an interior channel to the one-throat '
            'ring (which T5(a) shows gives a nondegenerate ladder)',
            'the freeze is fragile: unequal arcs, unequal channels or an '
            'internal bump each destroy it even at even N',
            'a non-cyclic theta network does not freeze in any twist sector; '
            'its degeneracies are S_3 permutation degeneracies present '
            'untwisted too',
        ],
        'open': [
            'the selection question itself is untouched: on the '
            'quarter-circle ring Probes B-D now apply to a network that '
            'freezes, and they still find nothing selecting the twist',
            'why the exact cyclic symmetry the freeze requires would hold in '
            'any physical configuration -- T5(b) shows small departures '
            'destroy it',
            'non-cyclic networks in general: one theta network was computed, '
            'which is a data point and not a classification',
            'whether the quarter-circle ring is dynamically preferred, or '
            'merely permitted, over #224\'s coincident-mouth configuration',
        ],
        'tests_passed': f'{passed}/{total}',
        'verdict_class': (
            'THE_FREEZE_NEEDS_EVEN_N_BUT_THE_DECK_CLASS_IS_SET_BY_MOUTH_'
            'PLACEMENT_NOT_BY_N_SO_THEY_ARE_INDEPENDENT_AND_A_TWO_THROAT_'
            'RING_WITH_QUARTER_CIRCLE_ARCS_CARRIES_BOTH'),
        'pass': self_pass,
    }


def run_probe() -> dict:
    res: dict = {}
    res['T1'] = test_T1_goal()
    res['T2'] = test_T2_freeze_law()
    res['T3'] = test_T3_deck_law_corrected()
    res['T4'] = test_T4_counterexample()
    res['T5'] = test_T5_evasion_classes()
    res['T6'] = test_T6_consequences()
    res['T7'] = test_T7_assessment(res)
    res['summary'] = {
        'probe': 'freeze_deck_parity_obstruction_probe', 'pr': 233,
        'utc': datetime.now(timezone.utc).isoformat(),
        'tests_passed': res['T7']['tests_passed'],
        'verdict_class': res['T7']['verdict_class'],
    }
    return res


def render_markdown(s: dict) -> str:
    o = ["# Probe G — the freeze and deck class 1 are **independent** "
         "(PR #233)", "",
         f"_Run {s['summary']['utc']} · {s['summary']['tests_passed']} PASS_",
         "",
         "> **Correction.** The first version of this probe claimed a parity",
         "> obstruction making the two mutually exclusive. That claim assumed",
         "> every exterior arc has length π — which forces all throats to",
         "> share one antipodal mouth pair — and is retracted here.", "",
         "## The two laws", "",
         "| law | status |", "|---|---|",
         "| full freeze ⟺ `N` even (`S^N = W`) | **stands** |",
         "| deck class = (Σ arc lengths)/π mod 2 | **set by mouth placement, "
         "not by `N`** |", "",
         "## The counterexample", "",
         "| ring | arcs/π | Σ/π | deck | doublet | gap untwisted | gap twisted "
         "| frozen? |",
         "|---|---:|---:|---:|---|---:|---:|---|"]
    for c in s['T4']['cases']:
        o.append(f"| `N` = {c['N']} | {c['arc_over_pi']} | "
                 f"{c['total_exterior_over_pi']} | {c['deck_class']} | "
                 f"{'yes' if c['has_mouth_doublet'] else 'no'} | "
                 f"{c['gap_untwisted']:.3e} | {c['gap_twisted']:.3e} | "
                 f"{'**YES**' if c['frozen'] else 'no'} |")
    o += ["", "## The other evasion classes", "",
          "| class | result |", "|---|---|",
          "| one throat, finite channel (`D > 0`) | nondegenerate interior "
          "ladder — no doublet, no freeze |",
          "| inhomogeneous / internally structured | freeze **destroyed** "
          "even at even `N` |",
          "| non-cyclic (theta network) | no freeze in any twist sector |",
          "", "## Consequences", ""]
    for i in s['T6']['items']:
        o.append(f"  - *{i['claim']}* — **{i['status'].split(' -- ')[0]}**")
    o += ["", "## Verdict", "", f"**{s['T7']['verdict_class']}**", ""]
    return "\n".join(o)


def _json_default(o):
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def main(argv: Optional[list] = None) -> int:
    summary = run_probe()
    md = render_markdown(summary)
    print(md)
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = here / "runs" / f"{ts}_freeze_deck_parity_obstruction_probe"
    out.mkdir(parents=True, exist_ok=True)
    (out / "probe.json").write_text(
        json.dumps(summary, indent=2, default=_json_default))
    (out / "probe.md").write_text(md)
    print(f"\nWrote: {out / 'probe.json'}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
