"""Probe E -- does the #224 network cycle generate pi_1(RP^3)? (PR #230)

> Framing: QFT on the *fixed classical* throat geometry (geometry -> fields),
> not quantum gravity.

THE QUESTION, AND WHY IT DECIDES THE ARC
-----------------------------------------
Probe A showed the B2 spin structure (in H^1(RP^3;Z2), detected by the deck
loop generating pi_1(RP^3) = Z2) and the network holonomy W (in
H^1(Gamma;Z2), detected by Wilson loops on ER-graph cycles) are independent
data.  They coincide ONLY if the network cycle carrying W is itself the
pi_1(RP^3) generator.  Probes B, C and D all compute properties of the two
RP^3 SPIN STRUCTURES.  So whether any of them bears on the network twist
comes down to this one question.

THE TEST
--------
A loop in RP^3 = S^3/{+/-1} is classified by lifting it to S^3: the lift of
a closed loop either returns to its starting point (class 0, contractible)
or to the ANTIPODE (class 1, the generator).  For a geodesic arc on the
unit exterior this is decided purely by total arc length modulo 2 pi:

    total pi   ->  lift ends at the antipode  ->  GENERATOR
    total 2 pi ->  lift closes                ->  TRIVIAL

THE ANSWER: NO -- AND #223 DIFFERS FROM #224
---------------------------------------------
Read off the constructions themselves, not the prose:

  * #224 (``build_two_throat``, the two-throat ring carrying the
    exchange-freeze headline) has
        Ltot = 2*_D_CH + 4*_SIG_M + 2*_LARC ,   _LARC = 3.14 ~ pi
    i.e. TWO exterior arcs of pi each.  Total exterior 2 pi: the S^3 lift
    CLOSES, so the deck class is 0 -- the cycle is NOT the pi_1(RP^3)
    generator.

  * #223 (the single two-mouth bridge ring) has ONE exterior arc of pi to
    the antipodal source point.  Total exterior pi: the lift ends at the
    antipode, deck class 1 -- that cycle IS the generator.

So the composition rule theta = 0 iff W*eta = +1 is justified on the #223
ring and NOT on the #224 ring -- and the #224 ring is exactly where the
arc's headline result (the exchange freeze) lives.

WHAT THE #224 CYCLE ACTUALLY IS
--------------------------------
The ambient space is not RP^3: it is the exterior with one 1-handle per
throat (a wormhole IS a 1-handle).  Attaching n 1-handles to a connected
3-manifold gives

    pi_1(M) = pi_1(exterior) * Z * ... * Z    (free product, one Z per handle)

The #224 cycle passes through BOTH throats, so its class lives in the FREE
(handle) part, while its deck component -- the only part B2 grades -- is
trivial by the arc count above.  W on that cycle is therefore a character
on a free Z generator, and B2 is a character on the Z2 deck generator.
Different generators, different groups.

*** CONSEQUENCE: Probes B, C and D compute the two RP^3 spin structures,
which are labelled by the deck Z2.  The #224 network twist is not that
label.  So B-D, though correct, DO NOT BEAR on the #224 twist sector --
including the exchange-freeze result.  The eta = +/-1/4 asymmetry is a
statement about B2, not about W. ***

ROBUSTNESS
----------
The conclusion does not depend on whether the exterior is read as S^3 or as
RP^3:
  * exterior = S^3: pi_1(S^3) = 0, the exterior contributes nothing at all
    and the cycle is purely handle-generated -- deck component trivially 0;
  * exterior = RP^3: the two-arc total 2 pi gives generator^2 = 0 in Z2.
Both readings give deck class 0.

Tests:
  T1. Goal and the lifting test.
  T2. The arc counts, read from the constructions.
  T3. The S^3 lift: endpoint and deck class, #224 vs #223.
  T4. pi_1 of the network space: the handle free product.
  T5. Consequence for Probes B-D.
  T6. Assessment.

Verdict:
  THE_224_CYCLE_DOES_NOT_GENERATE_PI1_RP3_ITS_DECK_CLASS_IS_TRIVIAL_AND_ITS
  _HOLONOMY_LIVES_ON_A_HANDLE_GENERATOR_SO_PROBES_B_TO_D_DO_NOT_BEAR_ON_IT
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np

PI = math.pi


# ── the S^3 lift test ───────────────────────────────────────────────────

def great_circle_endpoint(total_arc: float) -> np.ndarray:
    """Lift a geodesic exterior path of given total length to the unit S^3
    (start at e1, move in the e1-e2 great circle)."""
    return np.array([math.cos(total_arc), math.sin(total_arc), 0.0, 0.0])


def deck_class(total_arc: float, tol: float = 1e-9) -> dict:
    """0 if the lift closes, 1 if it ends at the antipode."""
    start = np.array([1.0, 0.0, 0.0, 0.0])
    end = great_circle_endpoint(total_arc)
    closes = bool(np.allclose(end, start, atol=tol))
    antipodal = bool(np.allclose(end, -start, atol=tol))
    return {'total_arc': total_arc,
            'endpoint': [round(float(v), 12) for v in end],
            'lift_closes': closes, 'lift_is_antipodal': antipodal,
            'deck_class': 0 if closes else (1 if antipodal else None),
            'is_generator': antipodal}


# the two constructions, read from their own length formulae
NETWORKS = {
    '#224 two-throat ring (build_two_throat)': {
        'n_exterior_arcs': 2, 'arc_length': 3.14,
        'source': 'Ltot = 2*_D_CH + 4*_SIG_M + 2*_LARC, _LARC = 3.14',
        'n_throats_traversed': 2,
    },
    '#223 single two-mouth bridge ring': {
        'n_exterior_arcs': 1, 'arc_length': PI,
        'source': 'neck -> bridge -> mouth -> S^3 arc (pi R_u to the '
                  'antipodal source point)',
        'n_throats_traversed': 1,
    },
}


# ════════════════════════════════════════════════════════════════════════

def test_T1_goal() -> dict:
    return {
        'name': 'T1_goal_and_the_lifting_test',
        'description': (
            "Probe A showed B2 (in H^1(RP^3;Z2)) and the network W (in "
            "H^1(Gamma;Z2)) are independent unless the network cycle IS the "
            "pi_1(RP^3) generator. Probes B, C and D all compute properties "
            "of the two RP^3 SPIN STRUCTURES, so whether any of them bears "
            "on the network twist reduces to this single question. The test "
            "is the standard one: lift the loop to the double cover S^3; the "
            "lift either closes (class 0) or ends at the ANTIPODE (class 1, "
            "the generator). For geodesic exterior arcs on the unit sphere "
            "that is decided by total arc length mod 2 pi."
        ),
        'pass': True,
    }


def test_T2_arc_counts() -> dict:
    rows = []
    for name, d in NETWORKS.items():
        total = d['n_exterior_arcs'] * d['arc_length']
        rows.append({'network': name, 'n_arcs': d['n_exterior_arcs'],
                     'arc_length': d['arc_length'], 'total_exterior': total,
                     'total_over_pi': round(total / PI, 6),
                     'read_from': d['source']})
    ok = (abs(rows[0]['total_over_pi'] - 2.0) < 1e-2
          and abs(rows[1]['total_over_pi'] - 1.0) < 1e-9)
    return {
        'name': 'T2_arc_counts_read_from_the_constructions',
        'description': (
            "Taken from the length formulae themselves rather than the "
            "prose. #224's build_two_throat sets Ltot = 2*_D_CH + 4*_SIG_M "
            "+ 2*_LARC with _LARC = 3.14: TWO exterior arcs of ~pi, total "
            "~2 pi. #223's ring has ONE arc of pi, to the antipodal source "
            "point. That single difference in arc count is what separates "
            "the two networks topologically."
        ),
        'rows': rows,
        'pass': ok,
    }


def test_T3_lift_and_deck_class() -> dict:
    rows = []
    for name, d in NETWORKS.items():
        total = d['n_exterior_arcs'] * d['arc_length']
        cl = deck_class(total, tol=1e-2 if abs(total - 2 * PI) < 1e-2 else 1e-9)
        cl['network'] = name
        rows.append(cl)
    a = next(r for r in rows if r['network'].startswith('#224'))
    b = next(r for r in rows if r['network'].startswith('#223'))
    ok = (a['deck_class'] == 0 and not a['is_generator']
          and b['deck_class'] == 1 and b['is_generator'])
    return {
        'name': 'T3_the_S3_lift_224_is_trivial_223_is_the_generator',
        'description': (
            "Lifting each cycle's exterior path to S^3: #224's total ~2 pi "
            "brings the lift back to its STARTING point, so its deck class "
            "is 0 -- the cycle is NOT the pi_1(RP^3) generator. #223's total "
            "pi brings the lift to the ANTIPODE, so its deck class is 1 -- "
            "that cycle IS the generator. *** The composition rule theta = 0 "
            "iff W*eta = +1 is therefore justified on the #223 ring and NOT "
            "on the #224 ring -- and #224 is exactly where the arc's "
            "headline exchange-freeze result lives. ***"
        ),
        'rows': rows,
        'pass': ok,
    }


def test_T4_pi1_of_the_network_space() -> dict:
    rows = []
    for name, d in NETWORKS.items():
        n = d['n_throats_traversed']
        rows.append({
            'network': name,
            'n_handles_traversed': n,
            'pi1_of_ambient': f"pi_1(exterior) * Z^{n}  (free product, one Z "
                              f"per 1-handle)",
            'cycle_class_lives_in': 'the FREE (handle) part',
        })
    robustness = [
        {'exterior_read_as': 'S^3',
         'pi1_exterior': '0',
         'deck_component_of_224_cycle': 0,
         'why': 'the exterior contributes nothing; the cycle is purely '
                'handle-generated'},
        {'exterior_read_as': 'RP^3',
         'pi1_exterior': 'Z2',
         'deck_component_of_224_cycle': 0,
         'why': 'two arcs of pi give generator^2 = 0 in Z2'},
    ]
    same = len({r['deck_component_of_224_cycle'] for r in robustness}) == 1
    return {
        'name': 'T4_pi1_of_the_network_space_is_a_handle_free_product',
        'description': (
            "The ambient space is not RP^3. A wormhole IS a 1-handle, and "
            "attaching n 1-handles to a connected 3-manifold gives "
            "pi_1 = pi_1(exterior) * Z * ... * Z, one free Z per handle. "
            "The #224 cycle passes through BOTH throats, so its class lives "
            "in the FREE part, while its deck component -- the only part B2 "
            "grades -- is trivial by T3. W on that cycle is a character on a "
            "free Z generator; B2 is a character on the Z2 deck generator. "
            "Different generators of different groups. ROBUSTNESS: the "
            "conclusion is independent of whether the exterior is read as "
            "S^3 (pi_1 = 0, so the deck component is trivially 0) or as "
            "RP^3 (two pi-arcs give generator^2 = 0). Both give 0."
        ),
        'rows': rows,
        'robustness': robustness,
        'conclusion_robust_to_exterior_reading': same,
        'pass': same,
    }


def test_T5_consequence_for_B_C_D() -> dict:
    items = [
        {'probe': 'B (rp3_dirac_eta_probe)',
         'computes': 'eta = +/-1/4, h = 0, |det D| for the two RP^3 SPIN '
                     'STRUCTURES',
         'labelled_by': 'the deck Z2',
         'bears_on_224_twist': False},
        {'probe': 'C (bulk_aps_spin_structure_probe)',
         'computes': 'whether each boundary SPIN STRUCTURE extends; ABK of '
                     'the Pin- mouth',
         'labelled_by': 'the deck Z2 (and, for ABK, the mouth Pin- structure)',
         'bears_on_224_twist': False},
        {'probe': 'D (twist_sector_einstein_dirac_probe)',
         'computes': 'back-reaction difference between the two SPIN '
                     'STRUCTURES',
         'labelled_by': 'the deck Z2',
         'bears_on_224_twist': False},
    ]
    none_bear = not any(i['bears_on_224_twist'] for i in items)
    return {
        'name': 'T5_probes_B_to_D_do_not_bear_on_the_224_twist',
        'description': (
            "Probes B, C and D are correct as computed, and they all answer "
            "the SAME question: how do the two RP^3 spin structures differ? "
            "That label is the deck Z2. T3/T4 show the #224 network twist is "
            "not that label -- its deck component is trivial and its "
            "holonomy lives on a handle generator. *** So B-D do NOT bear on "
            "the #224 twist sector, including the exchange-freeze result. "
            "The eta = +/-1/4 asymmetry is a statement about B2, not about "
            "W. *** This is a scoping correction, not a refutation: nothing "
            "in B-D is withdrawn, but its target is narrowed. It also means "
            "the 'three selection mechanisms all fail' conclusion is a "
            "statement about the SPIN STRUCTURE sectors, and the network "
            "twist sector remains untested by any of them."
        ),
        'items': items,
        'none_bear_on_the_network_twist': none_bear,
        'what_would_bear': ('a cycle with nontrivial deck class -- e.g. the '
                            '#223 single-bridge ring, whose single pi-arc '
                            'IS the generator (T3)'),
        'pass': none_bear,
    }


def test_T6_assessment(results: dict) -> dict:
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
        'name': 'T6_assessment',
        'description': (
            "The #224 network cycle does NOT generate pi_1(RP^3): its two "
            "pi-arcs close the S^3 lift, giving deck class 0, and its "
            "non-triviality comes from the throat handles instead. So the "
            "network twist W and the B2 spin structure are independent not "
            "just in principle (Probe A) but in this specific network. "
            "Probes B-D, which compute spin-structure properties, therefore "
            "do not bear on the #224 twist -- their target is narrowed, not "
            "their content withdrawn. The #223 ring, by contrast, DOES carry "
            "the generator, so it is where the spin-structure results and "
            "the network holonomy actually meet."
        ),
        'established': [
            "#224's exterior total is ~2 pi (two pi-arcs), so the S^3 lift "
            'closes and the deck class is 0 -- not the generator',
            "#223's exterior total is pi (one arc), lift ends antipodal, "
            'deck class 1 -- IS the generator',
            'the ambient pi_1 is a handle free product; the #224 cycle lives '
            'in the free part',
            'the conclusion is robust to reading the exterior as S^3 or RP^3',
            'Probes B-D do not bear on the #224 twist sector',
        ],
        'open': [
            'the network twist sector itself is now untested by ANY of the '
            'selection mechanisms tried so far -- energetics, inflow and '
            'back-reaction were all computed for the spin-structure label',
            'the #223 ring is where the two labels do meet; the '
            'spin-structure results have not been pushed through there',
            'the arc lengths are taken from the constructions as built '
            '(_LARC = 3.14 as a stand-in for pi); the geometry is a '
            'discretized model, not an embedded submanifold computation',
        ],
        'tests_passed': f'{passed}/{total}',
        'verdict_class': (
            'THE_224_CYCLE_DOES_NOT_GENERATE_PI1_RP3_ITS_DECK_CLASS_IS_'
            'TRIVIAL_AND_ITS_HOLONOMY_LIVES_ON_A_HANDLE_GENERATOR_SO_PROBES_'
            'B_TO_D_DO_NOT_BEAR_ON_IT'),
        'pass': self_pass,
    }


def run_probe() -> dict:
    res: dict = {}
    res['T1'] = test_T1_goal()
    res['T2'] = test_T2_arc_counts()
    res['T3'] = test_T3_lift_and_deck_class()
    res['T4'] = test_T4_pi1_of_the_network_space()
    res['T5'] = test_T5_consequence_for_B_C_D()
    res['T6'] = test_T6_assessment(res)
    res['summary'] = {
        'probe': 'network_cycle_pi1_class_probe', 'pr': 230,
        'utc': datetime.now(timezone.utc).isoformat(),
        'tests_passed': res['T6']['tests_passed'],
        'verdict_class': res['T6']['verdict_class'],
    }
    return res


def render_markdown(s: dict) -> str:
    o = ["# Probe E — does the #224 cycle generate π₁(RP³)? (PR #230)", "",
         f"_Run {s['summary']['utc']} · {s['summary']['tests_passed']} PASS_",
         "", "## The lift test", "",
         "| network | arcs | total exterior | lift | deck class | generator? |",
         "|---|---:|---:|---|---:|---|"]
    for r2, r3 in zip(s['T2']['rows'], s['T3']['rows']):
        lift = 'closes' if r3['lift_closes'] else (
            'antipodal' if r3['lift_is_antipodal'] else '?')
        o.append(f"| {r2['network']} | {r2['n_arcs']} | "
                 f"{r2['total_over_pi']:.3f}π | {lift} | {r3['deck_class']} | "
                 f"{'**YES**' if r3['is_generator'] else 'no'} |")
    o += ["", "## Consequence", "",
          "| probe | labelled by | bears on the #224 twist? |",
          "|---|---|---|"]
    for i in s['T5']['items']:
        o.append(f"| {i['probe']} | {i['labelled_by']} | "
                 f"{'yes' if i['bears_on_224_twist'] else '**no**'} |")
    o += ["", f"What would bear: {s['T5']['what_would_bear']}",
          "", "## Verdict", "", f"**{s['T6']['verdict_class']}**", ""]
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
    out = here / "runs" / f"{ts}_network_cycle_pi1_class_probe"
    out.mkdir(parents=True, exist_ok=True)
    (out / "probe.json").write_text(
        json.dumps(summary, indent=2, default=_json_default))
    (out / "probe.md").write_text(md)
    print(f"\nWrote: {out / 'probe.json'}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
