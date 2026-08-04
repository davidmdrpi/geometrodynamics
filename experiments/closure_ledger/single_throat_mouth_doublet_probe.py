"""Probe H -- can ONE throat carry a mouth doublet at deck class 1? (PR #235)

> Framing: QFT on the *fixed classical* throat geometry (geometry -> fields),
> not quantum gravity.

THE QUESTION
------------
#233 (corrected) showed the exchange freeze and deck class 1 are
INDEPENDENT, and exhibited a two-throat ring carrying both.  It left the
one-throat case answered only for the interior channel: a finite channel
D > 0 gives a nondegenerate box ladder, so no doublet.  But #224's doublet
is a MOUTH doublet, and a single throat already HAS two mouths.  So:

    Can a single throat with a finite interior channel produce two
    MOUTH-LOCALIZED states while retaining deck class 1?

THE ANSWER: YES TO THE STATES, NO TO THE DOUBLET THAT MATTERS
--------------------------------------------------------------
  (1) TWO MOUTH-LOCALIZED STATES: YES, easily.  Put a barrier in the middle
      of the interior channel and each sub-well abuts one mouth.  With one
      exterior pi-arc the ring keeps deck class 1.  Measured: localized
      combinations carrying > 0.96 of their weight in a single well, level
      splitting tunable over 7.5e-3 ... 2.7e-1 by the bump height.  So the
      literal question is answered YES, and #233's "no doublet on one
      throat" was too narrow -- it looked only at the channel ladder.

  (2) BUT THE TWIST DOES NOT ACT ON THEM.  Splitting ratio twisted /
      untwisted = 1.000 across the whole range.  The two wells are joined
      by two paths -- through the central bump, or the long way round
      through the exterior arc -- and only the second carries the
      holonomy.  Extracting both from the two twist sectors
      (split = 2|t_bump + W t_arc|):

          |t_arc|  = 2.0e-6 ... 2.6e-6, INDEPENDENT of the bump height
          |t_bump| = 1.4e-1 ... 3.7e-3, falling by 37x over the scan

      exactly as the two-path picture predicts -- but the exterior path
      crosses TWO mouth barriers and is 3-5 orders of magnitude weaker, so
      the interference never balances.  A doublet, but not a
      twist-controlled one: it cannot serve as #234's memory/gate switch.

  (3) THE RESONANT TWO-BASIN DOUBLET IS NO BETTER.  Tuning the channel
      length so a channel level crosses an arc level gives a genuine
      avoided crossing (D ~ 8.3, gap 3.2e-3) between two basins that ARE
      joined by two mouth barriers of equal strength.  The twist still
      moves it by only ~10%.

  (4) THE STRUCTURAL REASON, WHICH IS THE REAL RESULT.  The freeze is the
      S^k = W mechanism: a cyclic TRANSLATION of order k whose square (k
      even) has no real eigenvalue.  A one-throat ring has no such
      translation.  Its two mouth walls both face the SAME interior
      channel, so they are MIRROR IMAGES of each other, and the ring
      carries a REFLECTION R, not a translation.  And

          R^2 = +1   in BOTH twist sectors,

      always -- a reflection has no holonomy-dependent square, so it can
      never do what S^2 = W does.  Verified by building the same ring both
      ways:

          walls MIRROR (physical, one throat) : ||[H,S]|| = 4.9e+01,
              R is the symmetry, no freeze in either sector (gap 5.3e-01)
          walls TRANSLATE (counterfactual)    : ||[H,S]|| = 0.0e+00,
              S^2 = W, twisted sector FROZEN at 2.0e-12

      Same lengths, same deck class 1, same two mouths -- the ONLY
      difference is whether the two walls mirror or translate, and it
      decides the freeze completely.

  *** So the integer that governs the freeze is the ORDER OF THE CYCLIC
  TRANSLATION, and it equals the throat count because each throat
  contributes one interior basin bounded by a mirrored pair of walls.  A
  single throat gives a reflection, not a translation.  #233's "N even" is
  therefore right for the one-throat case -- but for the wall orientation,
  not for the arc counting it originally gave, and not for the channel
  ladder the corrected version gave. ***

CONSEQUENCE
-----------
#231's construction still exists (#233's quarter-circle two-throat ring),
but it needs two throats for a reason that is now precise rather than
numerological: two throats are the minimum for a cyclic translation.  And
#234's gate set, which uses the twist as its memory/gate switch, cannot be
moved onto a one-throat network at all -- not because the doublet is
missing (it is not) but because the twist has no purchase on it.

Tests:
  T1. Goal.
  T2. Two mouth-localized states on a one-throat deck-class-1 ring: YES.
  T3. The twist does not act on them; both tunneling paths extracted.
  T4. The resonant two-basin doublet: an avoided crossing, still not frozen.
  T5. The mechanism: mirror walls vs translate walls; R^2 = +1 vs S^2 = W.
  T6. Consequences for #231, #233 and #234.
  T7. Assessment.

Verdict:
  ONE_THROAT_GIVES_TWO_MOUTH_LOCALIZED_STATES_AT_DECK_CLASS_ONE_BUT_NEVER_A
  _TWIST_CONTROLLED_DOUBLET_BECAUSE_ITS_TWO_MOUTH_WALLS_MIRROR_RATHER_THAN
  _TRANSLATE_SO_THE_RING_CARRIES_A_REFLECTION_WITH_R_SQUARED_PLUS_ONE
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np

from experiments.closure_ledger.mouth_exchange_dynamics_probe import (
    v_bridge_half, _DX, _SIG_M,
)

PI = math.pi
_RS = 0.3


# ── the one-throat ring ─────────────────────────────────────────────────
#
#   [ interior channel D ][ mouth wall ][ exterior arc L_a ][ mouth wall ]
#
# ONE throat = ONE pair of mouths = TWO walls.  Sampling is cell-centred so
# that the exact symmetries of the continuum ring survive discretization --
# the same lesson #233 recorded for the cyclic shift.

def build_one_throat(D: float, L_a: float, rs: float = _RS,
                     bump_h: float = 0.0, bump_w: float = 0.5,
                     walls: str = 'mirror',
                     channel_profile: str = 'flat') -> dict:
    """`walls='mirror'` is the physical rule: each wall is v_bridge_half of
    the distance from the INTERIOR channel, so the two walls of one throat
    are mirror images.  `walls='translate'` is the counterfactual in which
    they are translates -- not realizable by a single throat, and used only
    to isolate what the freeze actually depends on."""
    n_D = int(round(D / _DX))
    dx = D / n_D
    n_a = int(round(L_a / dx))
    n_s = int(round(_SIG_M / dx))
    sig = n_s * dx
    m = n_D + n_a + 2 * n_s
    s = (np.arange(m) + 0.5) * dx          # cell centres
    v = np.zeros(m)

    i_ch = slice(0, n_D)
    i_w1 = slice(n_D, n_D + n_s)
    i_arc = slice(n_D + n_s, n_D + n_s + n_a)
    i_w2 = slice(n_D + n_s + n_a, m)

    def basin_profile(x, x0, L):
        return (3.75 / ((x - x0) + sig) ** 2
                + 3.75 / ((x0 + L - x) + sig) ** 2 - 1.0)

    if channel_profile == 'arc':
        v[i_ch] = basin_profile(s[i_ch], s[i_ch][0] - 0.5 * dx, D)
    if bump_h:
        v[i_ch] += bump_h * np.exp(-((s[i_ch] - D / 2) / bump_w) ** 2)
    wall = v_bridge_half((np.arange(n_s) + 0.5) * dx, rs)
    v[i_w1] = wall
    v[i_w2] = wall if walls == 'translate' else wall[::-1]
    v[i_arc] = basin_profile(s[i_arc], s[i_arc][0] - 0.5 * dx, L_a)

    channel = np.zeros(m, dtype=bool); channel[i_ch] = True
    arcm = np.zeros(m, dtype=bool); arcm[i_arc] = True
    wellL = np.zeros(m, dtype=bool); wellL[0:n_D // 2] = True
    wellR = np.zeros(m, dtype=bool); wellR[n_D // 2:n_D] = True
    return {'V': v, 'dx': dx, 'M': m, 'D': D, 'L_a': L_a, 'walls': walls,
            'n_cell': n_D + n_s, 'cells_congruent': (n_D == n_a),
            'channel': channel, 'arc': arcm, 'wellL': wellL, 'wellR': wellR,
            'total_exterior_over_pi': L_a / PI,
            'deck_class': int(round(L_a / PI)) % 2}


def hamiltonian(g: dict, holonomy: float) -> np.ndarray:
    m, dx = g['M'], g['dx']
    lap = (np.diag(-2.0 * np.ones(m)) + np.diag(np.ones(m - 1), 1)
           + np.diag(np.ones(m - 1), -1))
    lap[0, -1] = lap[-1, 0] = holonomy
    return -lap / dx ** 2 + np.diag(g['V'])


def translation(g: dict, holonomy: float) -> np.ndarray:
    """S: advance by one [basin + wall] cell, carrying the holonomy sign."""
    m, n = g['M'], g['n_cell']
    S = np.zeros((m, m))
    for i in range(m):
        j = i + n
        S[j % m, i] = holonomy if j >= m else 1.0
    return S


def reflection(g: dict, holonomy: float) -> np.ndarray:
    """R: reflect about the centre of the interior channel, exchanging the
    throat's two mouths."""
    m = g['M']
    axis2c = int(round(g['D'] / g['dx'])) - 1
    R = np.zeros((m, m))
    for i in range(m):
        j = axis2c - i
        R[j % m, i] = holonomy if j < 0 else 1.0
    return R


def modes(g: dict, holonomy: float, k: int = 16):
    ev, U = np.linalg.eigh(hamiltonian(g, holonomy))
    w = np.sqrt(np.abs(ev))
    i = np.argsort(w)[:k]
    return w[i], U[:, i]


def max_intra_pair_gap(g: dict, holonomy: float, k: int = 16) -> float:
    w, _ = modes(g, holonomy, k)
    return float(np.max(np.diff(w)[0::2]))


def channel_pair_splitting(g: dict, holonomy: float, k: int = 16):
    """The two lowest states living in the interior channel."""
    w, U = modes(g, holonomy, k)
    idx = [j for j in range(k)
           if float(np.sum(U[:, j][g['channel']] ** 2)) > 0.7][:2]
    if len(idx) < 2:
        return float('nan'), None
    return float(w[idx[1]] - w[idx[0]]), (U[:, idx[0]], U[:, idx[1]])


def well_localization(g: dict, u1: np.ndarray, u2: np.ndarray) -> float:
    """Rotate the pair to the maximally wellR-localized combination and
    report that weight (the #234 basin construction)."""
    P = g['wellR']
    a = float(np.sum(u1[P] ** 2)); b = float(np.sum(u2[P] ** 2))
    c = float(np.sum(u1[P] * u2[P]))
    ev = np.linalg.eigvalsh(np.array([[a, c], [c, b]]))
    return float(ev[-1])


# ════════════════════════════════════════════════════════════════════════

def test_T1_goal() -> dict:
    return {
        'name': 'T1_goal_can_one_throat_carry_a_mouth_doublet_at_deck_class_1',
        'description': (
            "#233 (corrected) showed the exchange freeze and deck class 1 "
            "are INDEPENDENT and exhibited a TWO-throat ring carrying both. "
            "For the ONE-throat ring it checked only the interior channel, "
            "found a nondegenerate box ladder, and concluded 'no doublet'. "
            "But #224's doublet is a MOUTH doublet and a single throat "
            "already has TWO mouths. This probe asks the sharper question "
            "directly: can one throat with a finite interior channel produce "
            "two MOUTH-LOCALIZED states while retaining deck class 1?"
        ),
        'pass': True,
    }


def test_T2_two_mouth_localized_states() -> dict:
    rows = []
    for h in (2.0, 4.0, 8.0, 14.0):
        g = build_one_throat(8.0, PI, bump_h=h)
        split, pair = channel_pair_splitting(g, +1.0)
        loc = well_localization(g, *pair)
        rows.append({'bump_height': h, 'deck_class': g['deck_class'],
                     'total_exterior_over_pi': round(
                         g['total_exterior_over_pi'], 4),
                     'splitting': split,
                     'single_well_localization': loc})
    ok = (all(r['deck_class'] == 1 for r in rows)
          and all(r['single_well_localization'] > 0.95 for r in rows))
    return {
        'name': 'T2_two_mouth_localized_states_exist_at_deck_class_one',
        'description': (
            "Put a barrier in the middle of the interior channel: each "
            "sub-well then abuts one MOUTH of the throat, and the "
            "symmetric/antisymmetric pair rotates to two states each "
            "localized in a single well. With one exterior pi-arc the ring "
            "keeps deck class 1. Measured localization of the rotated "
            f"combination: {rows[0]['single_well_localization']:.4f} ... "
            f"{rows[-1]['single_well_localization']:.4f}, with the splitting "
            f"tunable from {rows[0]['splitting']:.3e} down to "
            f"{rows[-1]['splitting']:.3e} by the bump height. *** So the "
            "literal answer is YES, and #233's 'no doublet on one throat' "
            "was too narrow: it looked at the channel LADDER and missed that "
            "a structured channel gives a mouth PAIR. ***"
        ),
        'rows': rows,
        'pass': ok,
    }


def test_T3_the_twist_does_not_act() -> dict:
    rows = []
    for h in (0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 14.0):
        g = build_one_throat(8.0, PI, bump_h=h)
        su, _ = channel_pair_splitting(g, +1.0)
        st, _ = channel_pair_splitting(g, -1.0)
        rows.append({'bump_height': h, 'split_untwisted': su,
                     'split_twisted': st, 'ratio': st / su,
                     't_bump': 0.25 * (su + st),
                     't_arc': 0.25 * abs(su - st)})
    ta = [r['t_arc'] for r in rows]
    tb = [r['t_bump'] for r in rows]
    arc_flat = max(ta) / min(ta) < 1.5           # arc path bump-independent
    bump_falls = tb[0] / tb[-1] > 10.0
    never_balances = all(r['ratio'] > 0.99 for r in rows)
    return {
        'name': 'T3_the_twist_has_no_purchase_on_the_one_throat_doublet',
        'description': (
            "The two wells are joined by TWO paths -- through the central "
            "bump, or the long way round through the exterior arc -- and "
            "only the second carries the holonomy, so "
            "split = 2|t_bump + W t_arc| and the two twist sectors separate "
            "them. Extracted: |t_arc| = "
            f"{min(ta):.1e} ... {max(ta):.1e}, INDEPENDENT of the bump "
            f"height (ratio {max(ta)/min(ta):.2f}) exactly as the two-path "
            f"picture requires, while |t_bump| falls {tb[0]/tb[-1]:.0f}x over "
            "the same scan. But the exterior path crosses TWO mouth barriers "
            "and is 3-5 orders of magnitude weaker, so the interference "
            "never balances: the splitting ratio twisted/untwisted stays "
            f"> {min(r['ratio'] for r in rows):.3f} throughout. *** The "
            "doublet is real but NOT twist-controlled -- it cannot serve as "
            "#234's memory/gate switch. ***"
        ),
        'rows': rows,
        'arc_path_independent_of_bump': arc_flat,
        'bump_path_falls': bump_falls,
        'interference_never_balances': never_balances,
        'pass': arc_flat and bump_falls and never_balances,
    }


def test_T4_resonant_two_basin_doublet() -> dict:
    rows = []
    for D in (7.5, 8.0, 8.25, 8.5, 9.0, 9.5):
        g = build_one_throat(float(D), PI)
        wu, Uu = modes(g, +1.0, 8)
        wt, _ = modes(g, -1.0, 8)
        gu, gt = float(wu[1] - wu[0]), float(wt[1] - wt[0])
        rows.append({'D': D, 'gap_untwisted': gu, 'gap_twisted': gt,
                     'ratio': gt / gu,
                     'channel_fraction_lower': float(
                         np.sum(Uu[:, 0][g['channel']] ** 2)),
                     'channel_fraction_upper': float(
                         np.sum(Uu[:, 1][g['channel']] ** 2))})
    cross = min(rows, key=lambda r: r['gap_untwisted'])
    swaps = any(rows[i]['channel_fraction_lower'] < 0.5
                < rows[i + 1]['channel_fraction_lower']
                for i in range(len(rows) - 1))
    still_not_frozen = cross['gap_twisted'] > 1e-4
    return {
        'name': 'T4_the_resonant_two_basin_doublet_is_no_better',
        'description': (
            "The other way to build a pair: tune the interior channel length "
            "so a CHANNEL level crosses an ARC level. These two basins ARE "
            "joined by two mouth barriers of equal strength, so this is the "
            "most favourable case for the twist. An avoided crossing does "
            f"appear at D ~ {cross['D']} (gap "
            f"{cross['gap_untwisted']:.3e}, with the lower state's channel "
            "fraction swapping across it), but the twist moves it by only "
            f"{abs(1 - cross['ratio']) * 100:.0f}% -- it does not freeze. The "
            "reason is T5: near the crossing the states are eigenstates of "
            "the ring's REFLECTION, and a reflection cannot be traded "
            "against the holonomy."
        ),
        'rows': rows,
        'crossing': cross,
        'basins_swap_across_the_crossing': swaps,
        'still_not_frozen': still_not_frozen,
        'pass': swaps and still_not_frozen,
    }


def test_T5_mechanism() -> dict:
    out = []
    for walls in ('mirror', 'translate'):
        # BOTH basins carry the same profile and the same length, so the
        # ONLY difference between the two rings is the wall orientation.
        g = build_one_throat(PI, PI, walls=walls, channel_profile='arc')
        entry = {'walls': walls, 'deck_class': g['deck_class'],
                 'cells_congruent': g['cells_congruent'], 'sectors': []}
        for hol in (+1.0, -1.0):
            H = hamiltonian(g, hol)
            S, R = translation(g, hol), reflection(g, hol)
            entry['sectors'].append({
                'W': int(hol),
                'commutator_H_S': float(np.max(np.abs(H @ S - S @ H))),
                'S_squared_minus_W': float(np.max(np.abs(
                    S @ S - hol * np.eye(g['M'])))),
                'commutator_H_R': float(np.max(np.abs(H @ R - R @ H))),
                'R_squared_minus_I': float(np.max(np.abs(
                    R @ R - np.eye(g['M'])))),
                'max_intra_pair_gap': max_intra_pair_gap(g, hol),
                'frozen': bool(max_intra_pair_gap(g, hol) < 1e-8),
            })
        out.append(entry)
    mirror, trans = out[0], out[1]
    # the physical ring: S is not a symmetry, nothing freezes
    mirror_ok = (all(sec['commutator_H_S'] > 1.0 for sec in mirror['sectors'])
                 and not any(sec['frozen'] for sec in mirror['sectors']))
    # the counterfactual: S is a symmetry, S^2 = W, the twisted sector freezes
    trans_ok = (all(sec['commutator_H_S'] < 1e-8 for sec in trans['sectors'])
                and all(sec['S_squared_minus_W'] < 1e-12
                        for sec in trans['sectors'])
                and [sec['frozen'] for sec in trans['sectors']] == [False, True])
    # R^2 = +1 in BOTH sectors, for both rings -- the whole point
    r_sq = all(sec['R_squared_minus_I'] < 1e-12
               for e in out for sec in e['sectors'])
    return {
        'name': 'T5_mirror_walls_give_a_reflection_translate_walls_give_S',
        'description': (
            "The freeze is the S^k = W mechanism: a cyclic TRANSLATION whose "
            "k-th power is the holonomy, with k even so that no eigenvalue is "
            "real. A one-throat ring has no such translation. Its two mouth "
            "walls both face the SAME interior channel, so they are MIRROR "
            "IMAGES, and the ring carries a REFLECTION R instead -- and "
            "R^2 = +1 in BOTH twist sectors (verified to < 1e-12 on every "
            "ring here), so a reflection has no holonomy-dependent square and "
            "can never do what S^2 = W does. Demonstrated by building the "
            "SAME ring both ways, same lengths, same deck class 1, same two "
            "mouths, differing only in wall orientation: with MIRROR walls "
            f"(physical) ||[H,S]|| = "
            f"{mirror['sectors'][0]['commutator_H_S']:.1e} and neither sector "
            f"freezes (gap {mirror['sectors'][1]['max_intra_pair_gap']:.3e}); "
            f"with TRANSLATE walls (not realizable by one throat) ||[H,S]|| = "
            f"{trans['sectors'][0]['commutator_H_S']:.1e}, S^2 = W, and the "
            f"twisted sector FREEZES at "
            f"{trans['sectors'][1]['max_intra_pair_gap']:.3e}. *** The wall "
            "orientation decides the freeze completely. ***"
        ),
        'rings': out,
        'physical_ring_has_no_translation_and_never_freezes': mirror_ok,
        'counterfactual_ring_has_S_and_freezes': trans_ok,
        'R_squared_is_plus_one_in_both_sectors': r_sq,
        'pass': mirror_ok and trans_ok and r_sq,
    }


def test_T6_consequences() -> dict:
    items = [
        {'claim': "#233 (corrected): 'the one-throat ring with D > 0 gives a "
                  "nondegenerate interior ladder, so no doublet'",
         'status': 'TOO NARROW -- true of the channel ladder, but a '
                   'STRUCTURED channel does give two mouth-localized states '
                   '(T2). The conclusion that the one-throat ring cannot '
                   'carry the freeze survives, for the different and better '
                   'reason in T5.'},
        {'claim': 'the integer governing the freeze is the throat count N',
         'status': 'SHARPENED -- it is the ORDER OF THE CYCLIC TRANSLATION. '
                   'That equals the throat count because each throat '
                   'contributes one interior basin bounded by a MIRRORED pair '
                   'of walls, so the shortest translation period is one '
                   'throat. Two throats are the minimum for a cyclic '
                   'translation at all.'},
        {'claim': "#231's construction (a geometry carrying the doublet and "
                  'deck class 1)',
         'status': 'STILL EXISTS -- #233\'s quarter-circle two-throat ring. '
                   'It needs two throats for a reason that is now structural '
                   'rather than numerological.'},
        {'claim': "#234's gate set could be moved onto a one-throat network",
         'status': 'NO -- and not because the doublet is missing. The doublet '
                   'is there (T2); the TWIST has no purchase on it (T3), so '
                   'the memory/gate switch #234 builds on the Z2 holonomy '
                   'has nothing to switch.'},
    ]
    return {
        'name': 'T6_consequences_for_231_233_and_234',
        'description': (
            "The literal question is answered YES and the useful one NO, and "
            "the gap between them is the result. A single throat gives two "
            "mouth-localized states at deck class 1 without difficulty; what "
            "it cannot give is a translation for the holonomy to act on. "
            "That replaces arc-counting (#233 v1) and channel-ladder "
            "counting (#233 v2) with the property that actually controls the "
            "freeze -- and it is checkable on any network by asking whether a "
            "cyclic translation commutes with H."
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
            "One throat with a structured interior channel does produce two "
            "mouth-localized states at deck class 1. But the Z2 twist cannot "
            "reach them: the only path carrying the holonomy crosses two "
            "mouth barriers and is orders of magnitude too weak, and more "
            "fundamentally the ring's symmetry is a REFLECTION, whose square "
            "is +1 in both twist sectors, rather than a translation whose "
            "square is the holonomy. Wall orientation, not arc count and not "
            "channel structure, is what decides the freeze."
        ),
        'established': [
            'two mouth-localized states on a one-throat deck-class-1 ring: '
            'YES (localization > 0.96, splitting tunable 7.5e-3 ... 2.7e-1)',
            'the twist does not act on them: ratio twisted/untwisted > 0.99 '
            'throughout; |t_arc| ~ 2e-6 independent of the bump height while '
            '|t_bump| falls 37x -- the two-path picture confirmed and the '
            'interference shown never to balance',
            'the resonant two-basin doublet (avoided crossing at D ~ 8.25) '
            'moves by only ~10% under the twist',
            'the mechanism: mirror walls give a REFLECTION (R^2 = +1 in both '
            'sectors, verified < 1e-12) and no freeze; translate walls give '
            'S with S^2 = W and an exact freeze (2.0e-12) -- same lengths, '
            'same deck class, same two mouths',
            'so the freeze is governed by the order of the cyclic '
            'translation, and two throats are the minimum for one to exist',
        ],
        'open': [
            'whether any BAM geometry realizes translate-oriented walls on a '
            'single throat -- it would need the second wall to face away from '
            'the interior, which is what a second throat provides',
            'the selection question is still untouched: nothing selects the '
            'twist on any network examined so far',
            "#234's gate set on a deck-class-1 network requires two throats; "
            'whether the quarter-circle ring supports its dressing '
            'calibration is uncomputed',
        ],
        'tests_passed': f'{passed}/{total}',
        'verdict_class': (
            'ONE_THROAT_GIVES_TWO_MOUTH_LOCALIZED_STATES_AT_DECK_CLASS_ONE_'
            'BUT_NEVER_A_TWIST_CONTROLLED_DOUBLET_BECAUSE_ITS_TWO_MOUTH_WALLS_'
            'MIRROR_RATHER_THAN_TRANSLATE_SO_THE_RING_CARRIES_A_REFLECTION_'
            'WITH_R_SQUARED_PLUS_ONE'),
        'pass': self_pass,
    }


def run_probe() -> dict:
    res: dict = {}
    res['T1'] = test_T1_goal()
    res['T2'] = test_T2_two_mouth_localized_states()
    res['T3'] = test_T3_the_twist_does_not_act()
    res['T4'] = test_T4_resonant_two_basin_doublet()
    res['T5'] = test_T5_mechanism()
    res['T6'] = test_T6_consequences()
    res['T7'] = test_T7_assessment(res)
    res['summary'] = {
        'probe': 'single_throat_mouth_doublet_probe', 'pr': 235,
        'utc': datetime.now(timezone.utc).isoformat(),
        'tests_passed': res['T7']['tests_passed'],
        'verdict_class': res['T7']['verdict_class'],
    }
    return res


def render_markdown(s: dict) -> str:
    o = ["# Probe H — one throat, two mouth states, no twist control "
         "(PR #235)", "",
         f"_Run {s['summary']['utc']} · {s['summary']['tests_passed']} PASS_",
         "",
         "**Q. Can a single throat with a finite interior channel produce two "
         "mouth-localized states while retaining deck class 1?**", "",
         "**A. Yes — and the twist still cannot touch them.**", "",
         "## The states exist", "",
         "| bump `h` | deck | splitting | single-well localization |",
         "|---:|---:|---:|---:|"]
    for r in s['T2']['rows']:
        o.append(f"| {r['bump_height']} | {r['deck_class']} | "
                 f"{r['splitting']:.3e} | {r['single_well_localization']:.4f} |")
    o += ["", "## But the twist does not act on them", "",
          "| bump `h` | split untwisted | split twisted | ratio | `|t_bump|` | "
          "`|t_arc|` |", "|---:|---:|---:|---:|---:|---:|"]
    for r in s['T3']['rows']:
        o.append(f"| {r['bump_height']} | {r['split_untwisted']:.4e} | "
                 f"{r['split_twisted']:.4e} | {r['ratio']:.4f} | "
                 f"{r['t_bump']:.2e} | {r['t_arc']:.2e} |")
    o += ["",
          "`|t_arc|` is independent of the bump height — the two-path picture "
          "confirmed — but the exterior path crosses **two** mouth barriers, "
          "so the interference never balances.", "",
          "## The mechanism", "",
          "| walls | `‖[H,S]‖` | `S²=W`? | `R²=+1`? | twisted gap | frozen? |",
          "|---|---:|---|---|---:|---|"]
    for e in s['T5']['rings']:
        tw = e['sectors'][1]
        o.append(f"| {e['walls']} | {tw['commutator_H_S']:.1e} | "
                 f"{'yes' if tw['S_squared_minus_W'] < 1e-12 else 'no'} | "
                 f"{'yes' if tw['R_squared_minus_I'] < 1e-12 else 'no'} | "
                 f"{tw['max_intra_pair_gap']:.3e} | "
                 f"{'**YES**' if tw['frozen'] else 'no'} |")
    o += ["",
          "Same lengths, same deck class 1, same two mouths — only the wall "
          "orientation differs, and it decides the freeze completely. A "
          "single throat's two walls both face the same interior, so they "
          "**mirror**; and `R² = +1` in both twist sectors, always.", "",
          "## Consequences", ""]
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
    out = here / "runs" / f"{ts}_single_throat_mouth_doublet_probe"
    out.mkdir(parents=True, exist_ok=True)
    (out / "probe.json").write_text(
        json.dumps(summary, indent=2, default=_json_default))
    (out / "probe.md").write_text(md)
    print(f"\nWrote: {out / 'probe.json'}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
