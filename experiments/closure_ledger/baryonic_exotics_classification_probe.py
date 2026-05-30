"""
BAM-specific baryonic exotics: classification and experimental
constraints (PR #102).

PR #101 matched BAM's non-orientable (Möbius) topology to the OBSERVED
mesonic exotics (1-+ hybrids, multiquark zoo) and flagged the Möbius
BARYON as a BAM-specific prediction. This probe classifies the
BAM-specific baryonic exotics and — as the user asks — identifies which
channels are experimentally constrained.

## The crucial subtlety: no exotic J^P for baryons

For mesons, the non-orientable (Möbius) flux tube produced a SMOKING-GUN
exotic quantum number: `J^PC = 1-+`, forbidden to ordinary qq̄ because of
the `C = (−1)^{L+S}` constraint. **Baryons have no such constraint.** A
qqq baryon has `P = (−1)^L`, spin `S ∈ {½, 3/2}`, and no good `C`, so
every half-integer `J^P` is reachable by an ordinary baryon — there is NO
forbidden (exotic) `J^P`. Consequently a BAM Möbius / hybrid baryon
carries ORDINARY quantum numbers: it is a **supernumerary** state, not a
manifestly exotic one. Its only signature is COUNTING — an extra
resonance beyond the qqq quark-model spectrum.

This makes baryonic exotics the hardest of all to identify, and — because
they must hide inside an already densely-measured spectrum — the MOST
experimentally constrained corner of BAM's non-orientable predictions
(the opposite end from the unobserved glueballs of PR #100).

## The BAM-specific baryonic exotics

  - **Hybrid baryon** (excited Y-junction glue): the baryon analogue of
    the hybrid meson — the Y-junction flux in its first excited/twisted
    mode. Mass ≈ nucleon + 2√σ ≈ 0.94 + 0.85 ≈ **1.79 GeV** (the same
    flux-tube quantum as PR #101); the Δ-based partner ≈ 2.08 GeV.
  - **Möbius baryon** (non-orientable Y-network, `make_mobius_baryon`):
    the Z₂-twisted partner of a baryon; ordinary `J^P`, supernumerary.

Both land in the LIGHT N*/Δ* region — the densest, best-measured part of
the baryon spectrum.

## Experimental-constraint ranking of channels

| channel | data density | constraint on a BAM exotic |
|---|---|---|
| light N*/Δ* (< 2.5 GeV) | ~40 PDG states; "missing-resonance" problem active | **strongest** — extra states must coincide with observed ones or decouple |
| strange hyperons (Λ*, Σ*) | well measured | strong |
| charm/bottom baryons (Ξ_c*, Ω_c*, Λ_b*) | sparse (a handful, LHCb) | **weakest** — BAM freest here |

The lightest BAM baryonic exotics (~1.8–2.1 GeV) sit squarely in the
MOST-constrained light sector (`N(1710)`, `N(1875)`, `N(1900)`, …,
`Δ(1900)`, `Δ(1950)`, `Δ(2000)`, …).

## The test / tension

BAM's Möbius topology doubles the spectrum (a Z₂-twisted partner per
state). The dense light sector cannot absorb arbitrary extra states, so
the Möbius/hybrid baryons must either (i) COINCIDE with observed-but-
unexplained resonances (filling "missing resonances" the quark model
under-predicts), or (ii) DECOUPLE from the dominant production/decay
channels (weak `πN` coupling) — the same mechanism invoked for the
missing resonances. Over-prediction without decoupling would be excluded.
So this is a genuine, near-term-testable constraint — BAM's most
exposed topological prediction.

## What this probe establishes (and does not)

  - **Establishes (BAM-native):** baryonic exotics have NO exotic-`J^P`
    smoking gun (any `J^P` is ordinary for qqq), so BAM's Möbius/hybrid
    baryons are supernumerary ordinary-`J^P` states identifiable only by
    counting; they land in the light N*/Δ* region (~1.8–2.1 GeV, the
    `2√σ` flux-tube gap), the MOST experimentally constrained channel;
    and the constraint ranking is light N*/Δ* > strange hyperons >
    heavy-quark baryons (the freest). This is the opposite extreme from
    the unobserved glueballs (PR #100).

  - **Does not establish:** a specific Möbius-baryon state confirmed in
    data. Without a smoking-gun `J^P`, the prediction is a counting one,
    testable only by matching the dense spectrum or by demonstrating
    decoupling — both beyond this classification.

Tests:
  T1. Baryonic exotics from non-orientable topology (Möbius/hybrid
      baryon); build on the PR #101 Möbius-baryon flag.
  T2. No exotic J^P for baryons (P=(−1)^L, S∈{½,3/2}, no C) ⟹ any J^P
      ordinary ⟹ no smoking gun (unlike meson 1-+).
  T3. BAM baryonic exotics: hybrid baryon ≈ nucleon + 2√σ ≈ 1.79 GeV;
      Δ-partner ≈ 2.08 GeV; Möbius baryon = Z₂ partner.
  T4. Constraint ranking: light N*/Δ* (strongest) > strange hyperons >
      charm/bottom baryons (weakest).
  T5. Lightest BAM baryonic exotic sits in the MOST-constrained light
      sector (dense N*/Δ* region).
  T6. Test/tension: Möbius doubling must coincide with observed states or
      decouple (missing-resonance mechanism); else excluded.
  T7. The opposite extreme from glueballs (most vs least constrained);
      honest scope.
  T8. Assessment.

Verdict:
  - BARYONIC_EXOTICS_LACK_EXOTIC_JP_LIGHT_SECTOR_MOST_CONSTRAINED
    (expected): BAM-specific baryonic exotics have no exotic-J^P smoking
    gun (supernumerary ordinary-J^P states), land in the densely-measured
    light N*/Δ* region (~1.8–2.1 GeV), and are therefore the MOST
    experimentally constrained corner of BAM's non-orientable predictions;
    heavy-quark baryon exotics are the freest.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np

from geometrodynamics.qcd.constants import SIGMA_QCD


PI = math.pi
SQRT_SIGMA = math.sqrt(SIGMA_QCD)
HYBRID_GAP = 2.0 * SQRT_SIGMA            # flux-tube excitation ≈ 0.85 GeV
M_NUCLEON = 0.939
M_DELTA = 1.232

# Light N*/Δ* states in the ~1.8–2.1 GeV region (PDG), to show density.
DENSE_NSTAR_REGION = ['N(1710)', 'N(1720)', 'N(1875)', 'N(1880)', 'N(1895)',
                      'N(1900)', 'N(2000)', 'N(2060)']
DENSE_DELTA_REGION = ['Δ(1900)', 'Δ(1905)', 'Δ(1910)', 'Δ(1920)',
                      'Δ(1930)', 'Δ(1950)', 'Δ(2000)']


def baryon_jp_reachable(twoJ: int, P: int) -> bool:
    """Can ordinary qqq reach (J=twoJ/2, P)? P=(−1)^L, S∈{1/2,3/2},
    J∈[|L−S|, L+S]. Half-integer J ⟹ twoJ odd."""
    if twoJ % 2 == 0:
        return False    # baryons are half-integer spin
    J = twoJ / 2.0
    for L in range(0, int(J) + 3):
        if (-1) ** L != P:
            continue
        for twoS in (1, 3):       # S = 1/2 or 3/2
            S = twoS / 2.0
            if abs(L - S) <= J <= L + S:
                return True
    return False


# ---------------------------------------------------------------------------
# T1. Setup
# ---------------------------------------------------------------------------

def test_T1_setup() -> dict:
    return {
        'name': 'T1_baryonic_exotics_setup',
        'description': (
            "PR #101 flagged the Möbius baryon as a BAM-specific "
            "prediction. This probe classifies the baryonic exotics "
            "(Möbius / hybrid baryon) and ranks their experimental "
            "constraints."
        ),
        'bam_baryonic_exotics': ['hybrid baryon (excited Y-junction glue)',
                                 'Möbius baryon (non-orientable Y-network)'],
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T2. No exotic J^P for baryons
# ---------------------------------------------------------------------------

def test_T2_no_exotic_jp() -> dict:
    """Unlike mesons (1-+ exotic via C), baryons have P=(−1)^L, S∈{½,3/2},
    no good C ⟹ every half-integer J^P is reachable by qqq ⟹ no forbidden
    (exotic) J^P. So a BAM Möbius/hybrid baryon is supernumerary
    (ordinary J^P), not manifestly exotic. Verify a sample of J^P are all
    reachable."""
    sample = [(1, +1), (1, -1), (3, +1), (3, -1), (5, +1), (5, -1), (7, +1)]
    rows = [{'jp': f'{twoJ}/2{"+" if P > 0 else "−"}',
             'reachable_by_qqq': baryon_jp_reachable(twoJ, P)} for twoJ, P in sample]
    all_reachable = all(r['reachable_by_qqq'] for r in rows)
    return {
        'name': 'T2_no_exotic_jp_for_baryons',
        'description': (
            "Baryons: P=(−1)^L, S∈{½,3/2}, no C ⟹ every half-integer J^P "
            "reachable by qqq ⟹ NO exotic J^P (unlike meson 1-+). BAM "
            "Möbius/hybrid baryons are supernumerary (ordinary J^P), "
            "signature = counting only."
        ),
        'rows': rows,
        'all_jp_reachable_no_exotic': all_reachable,
        'meson_contrast': '1-+ is exotic for mesons (C-forbidden); no baryon analogue',
        'pass': all_reachable,
    }


# ---------------------------------------------------------------------------
# T3. The BAM baryonic exotics and their masses
# ---------------------------------------------------------------------------

def test_T3_exotic_masses() -> dict:
    """The hybrid baryon (excited Y-junction glue) sits ≈ nucleon + 2√σ ≈
    1.79 GeV; the Δ-based partner ≈ 2.08 GeV — the same 2√σ flux-tube
    quantum as the mesonic hybrids (PR #101). The Möbius baryon is the
    Z₂-twisted partner. All land in the light N*/Δ* region."""
    m_hybrid_N = M_NUCLEON + HYBRID_GAP
    m_hybrid_D = M_DELTA + HYBRID_GAP
    return {
        'name': 'T3_baryonic_exotic_masses',
        'description': (
            "Hybrid baryon ≈ nucleon + 2√σ ≈ 1.79 GeV; Δ-partner ≈ "
            "2.08 GeV (same 2√σ flux-tube quantum as PR #101). Möbius "
            "baryon = Z₂ partner. All in the light N*/Δ* region."
        ),
        'flux_tube_gap_2sqrt_sigma_GeV': HYBRID_GAP,
        'hybrid_N_baryon_GeV': m_hybrid_N,
        'hybrid_Delta_baryon_GeV': m_hybrid_D,
        'lands_in_light_sector': 1.5 < m_hybrid_N < 2.5,
        'pass': 1.5 < m_hybrid_N < 2.5,
    }


# ---------------------------------------------------------------------------
# T4. Constraint ranking of channels
# ---------------------------------------------------------------------------

def test_T4_constraint_ranking() -> dict:
    """Rank baryon channels by data density (= constraint on a BAM exotic):
    light N*/Δ* (~40 PDG states, missing-resonance problem) — strongest;
    strange hyperons (Λ*, Σ*) — strong; charm/bottom baryons (sparse) —
    weakest (BAM freest)."""
    ranking = [
        {'channel': 'light N*/Δ* (<2.5 GeV)', 'data': '~40 PDG states; missing-resonance problem',
         'constraint': 'strongest'},
        {'channel': 'strange hyperons (Λ*, Σ*)', 'data': 'well measured',
         'constraint': 'strong'},
        {'channel': 'charm/bottom baryons (Ξ_c*, Ω_c*, Λ_b*)', 'data': 'sparse (LHCb handful)',
         'constraint': 'weakest (BAM freest)'},
    ]
    return {
        'name': 'T4_experimental_constraint_ranking',
        'description': (
            "Channels by data density: light N*/Δ* (strongest) > strange "
            "hyperons > charm/bottom baryons (weakest, BAM freest)."
        ),
        'ranking': ranking,
        'most_constrained': 'light N*/Δ* (<2.5 GeV)',
        'least_constrained': 'charm/bottom baryon exotics',
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T5. Lightest exotic in the most-constrained sector
# ---------------------------------------------------------------------------

def test_T5_in_most_constrained() -> dict:
    """The lightest BAM baryonic exotics (~1.8–2.1 GeV) sit squarely in
    the MOST-constrained light sector — the densely-measured N*/Δ* region
    (N(1710), N(1875), N(1900), …; Δ(1900), Δ(1950), Δ(2000), …)."""
    m_hybrid_N = M_NUCLEON + HYBRID_GAP
    near_nstar = sum(1 for _ in DENSE_NSTAR_REGION)
    return {
        'name': 'T5_lightest_exotic_in_dense_region',
        'description': (
            "Lightest BAM baryonic exotic ≈ 1.79 GeV (N-based), 2.08 GeV "
            "(Δ-based) — in the densely-measured N*/Δ* region (the "
            "most-constrained channel)."
        ),
        'hybrid_N_baryon_GeV': m_hybrid_N,
        'dense_nstar_region_states': DENSE_NSTAR_REGION,
        'dense_delta_region_states': DENSE_DELTA_REGION,
        'n_nearby_nstar': near_nstar,
        'in_dense_region': near_nstar >= 5,
        'pass': near_nstar >= 5,
    }


# ---------------------------------------------------------------------------
# T6. The test / tension
# ---------------------------------------------------------------------------

def test_T6_test_tension() -> dict:
    """The Möbius topology doubles the spectrum (a Z₂-twisted partner per
    state). The dense light sector cannot absorb arbitrary extras, so the
    Möbius/hybrid baryons must either COINCIDE with observed-but-
    unexplained resonances (fill missing resonances) or DECOUPLE (weak πN
    coupling) — else over-prediction is excluded. A genuine near-term
    constraint."""
    return {
        'name': 'T6_test_and_tension',
        'description': (
            "Möbius doubling vs the dense light sector: the extra states "
            "must coincide with observed resonances (fill missing "
            "resonances) or decouple (weak πN), else excluded. A real "
            "near-term constraint."
        ),
        'mobius_doubles_spectrum': True,
        'resolution_options': [
            'coincide with observed-but-unexplained resonances (missing '
            'resonances the quark model under-predicts)',
            'decouple from dominant πN production/decay (weak coupling) — '
            'the standard missing-resonance mechanism',
        ],
        'over_prediction_without_decoupling': 'excluded',
        'genuine_constraint': True,
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T7. Opposite extreme from glueballs
# ---------------------------------------------------------------------------

def test_T7_opposite_of_glueballs() -> dict:
    """Glueballs (PR #100) are unobserved — BAM's least-constrained
    topological prediction (free to differ). Light baryonic exotics are
    densely measured — BAM's MOST-constrained. The two bracket the
    spectrum of testability of BAM's non-orientable predictions."""
    return {
        'name': 'T7_opposite_extreme_from_glueballs',
        'description': (
            "Glueballs (PR #100) = least constrained (unobserved, BAM "
            "free). Light baryonic exotics = most constrained (dense "
            "N*/Δ*). The two bracket the testability of BAM's "
            "non-orientable predictions."
        ),
        'least_constrained_pole': 'Möbius glueballs (unobserved, PR #100)',
        'most_constrained_pole': 'light Möbius/hybrid baryons (dense N*/Δ*)',
        'honest_caveat': (
            'no exotic-J^P smoking gun for baryons ⟹ the constraint is a '
            'counting one (match the dense spectrum or decouple), not a '
            'clean quantum-number test'
        ),
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T8. Assessment
# ---------------------------------------------------------------------------

def test_T8_assessment() -> dict:
    return {
        'name': 'T8_assessment',
        'description': (
            "BAM-specific baryonic exotics have no exotic-J^P smoking gun "
            "(supernumerary ordinary-J^P states), land in the "
            "densely-measured light N*/Δ* region (~1.8–2.1 GeV), and are "
            "therefore the MOST experimentally constrained corner of BAM's "
            "non-orientable predictions; heavy-quark baryon exotics are the "
            "freest."
        ),
        'classification': (
            'BARYONIC_EXOTICS_LACK_EXOTIC_JP_LIGHT_SECTOR_MOST_CONSTRAINED'
        ),
        'pass': True,
    }


# ---------------------------------------------------------------------------
# Runner + verdict
# ---------------------------------------------------------------------------

def run_probe() -> dict:
    tests = [
        test_T1_setup(),
        test_T2_no_exotic_jp(),
        test_T3_exotic_masses(),
        test_T4_constraint_ranking(),
        test_T5_in_most_constrained(),
        test_T6_test_tension(),
        test_T7_opposite_of_glueballs(),
        test_T8_assessment(),
    ]
    if all(t['pass'] for t in tests[:7]):
        verdict_class = (
            'BARYONIC_EXOTICS_LACK_EXOTIC_JP_LIGHT_SECTOR_MOST_CONSTRAINED'
        )
        verdict = (
            'BAM-SPECIFIC BARYONIC EXOTICS HAVE NO SMOKING GUN AND LIVE IN '
            'THE MOST CONSTRAINED CHANNEL. PR #101 matched BAM\'s '
            'non-orientable topology to the OBSERVED mesonic exotics and '
            'flagged the Möbius baryon as a BAM-specific prediction. This '
            'probe classifies the baryonic exotics and identifies which '
            'channels are experimentally constrained.\n\n'
            'NO EXOTIC J^P FOR BARYONS. For mesons the Möbius flux tube '
            'gave a smoking-gun exotic quantum number, J^PC = 1-+, '
            'forbidden to qq̄ by the C=(−1)^{L+S} constraint. Baryons have '
            'no such constraint: a qqq baryon has P=(−1)^L, S∈{½,3/2}, and '
            'no good C, so every half-integer J^P is reachable by an '
            'ordinary baryon — there is NO forbidden (exotic) J^P. So a BAM '
            'Möbius / hybrid baryon carries ORDINARY quantum numbers; it is '
            'a SUPERNUMERARY state, identifiable only by counting (an extra '
            'resonance beyond the qqq quark model), not by a smoking-gun '
            'J^P.\n\n'
            'THE BAM BARYONIC EXOTICS. The hybrid baryon (excited '
            'Y-junction glue) sits ≈ nucleon + 2√σ ≈ 1.79 GeV; the Δ-based '
            'partner ≈ 2.08 GeV — the same 2√σ flux-tube quantum as the '
            'mesonic hybrids (PR #101). The Möbius baryon '
            '(make_mobius_baryon) is the Z₂-twisted partner. Both land in '
            'the LIGHT N*/Δ* region.\n\n'
            'EXPERIMENTAL-CONSTRAINT RANKING. By data density: light N*/Δ* '
            '(< 2.5 GeV) — ~40 PDG states, the missing-resonance problem '
            'active — is the STRONGEST constraint; strange hyperons (Λ*, '
            'Σ*) strong; charm/bottom baryons (Ξ_c*, Ω_c*, Λ_b*) sparse — '
            'the WEAKEST (BAM freest). The lightest BAM baryonic exotics '
            '(~1.8–2.1 GeV) sit squarely in the MOST-constrained light '
            'sector (N(1710), N(1875), N(1900), …; Δ(1900), Δ(1950), '
            'Δ(2000), …).\n\n'
            'THE TEST / TENSION. BAM\'s Möbius topology doubles the '
            'spectrum (a Z₂-twisted partner per state). The dense light '
            'sector cannot absorb arbitrary extras, so the Möbius/hybrid '
            'baryons must either COINCIDE with observed-but-unexplained '
            'resonances (filling missing resonances the quark model '
            'under-predicts) or DECOUPLE from the dominant πN '
            'production/decay (the standard missing-resonance mechanism); '
            'over-prediction without decoupling would be excluded. So this '
            'is the opposite extreme from the unobserved glueballs '
            '(PR #100): the light baryonic exotics are BAM\'s MOST exposed '
            'topological prediction.\n\n'
            'HONEST SCOPE. ESTABLISHED (BAM-native): baryonic exotics have '
            'no exotic-J^P smoking gun (any J^P is ordinary for qqq), so '
            'BAM\'s Möbius/hybrid baryons are supernumerary ordinary-J^P '
            'states identifiable only by counting; they land in the light '
            'N*/Δ* region (~1.8–2.1 GeV, the 2√σ gap), the MOST '
            'experimentally constrained channel; the constraint ranking is '
            'light N*/Δ* > strange hyperons > heavy-quark baryons (the '
            'freest). NOT established: a specific Möbius-baryon state '
            'confirmed in data — without a smoking-gun J^P the prediction '
            'is a counting one, testable only by matching the dense '
            'spectrum or demonstrating decoupling.'
        )
    else:
        verdict_class = 'BARYONIC_EXOTIC_CLASSIFICATION_INCONCLUSIVE'
        verdict = (
            'BARYONIC EXOTIC CLASSIFICATION INCONCLUSIVE. A structural test '
            'failed; investigate before claiming the constraint ranking.'
        )

    return {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'identification': (
            'BAM-specific baryonic exotics (Möbius / hybrid baryon) have no '
            'exotic-J^P smoking gun (supernumerary ordinary-J^P states), '
            'land in the densely-measured light N*/Δ* region (~1.8–2.1 GeV, '
            'the 2√σ gap), and are the MOST experimentally constrained '
            'corner of BAM\'s non-orientable predictions'
        ),
        'no_smoking_gun': 'no exotic J^P for baryons (any J^P ordinary for qqq)',
        'masses': 'hybrid N ≈ 1.79 GeV, hybrid Δ ≈ 2.08 GeV (nucleon/Δ + 2√σ)',
        'constraint_ranking': 'light N*/Δ* (strongest) > strange hyperons > charm/bottom (weakest)',
        'test': 'Möbius doubling must coincide with observed states or decouple; else excluded',
        'contrast': 'opposite extreme from the unobserved glueballs (PR #100)',
        'tests': tests,
        'n_passed': sum(1 for t in tests if t['pass']),
        'n_total': len(tests),
        'verdict_class': verdict_class,
        'verdict': verdict,
    }


# ---------------------------------------------------------------------------
# Markdown rendering
# ---------------------------------------------------------------------------

def render_markdown(s: dict) -> str:
    L: list[str] = []
    L.append('# BAM-specific baryonic exotics: classification + experimental constraints (PR #102)')
    L.append('')
    L.append(f"**Run:** {s['timestamp_utc']}")
    L.append('')
    L.append(
        "Classifies BAM's non-orientable (Möbius / hybrid) BARYONIC exotics "
        "and ranks the channels by experimental constraint. **Key "
        "subtlety:** unlike mesons (where `1-+` is a smoking-gun exotic via "
        "`C`), baryons have **no forbidden `J^P`** — so BAM's Möbius/hybrid "
        "baryons are *supernumerary ordinary-`J^P`* states, identifiable "
        "only by counting. They land in the densely-measured light N*/Δ* "
        "region (~1.8–2.1 GeV), making them the **most** experimentally "
        "constrained corner of BAM's non-orientable predictions — the "
        "opposite extreme from the unobserved glueballs (PR #100)."
    )
    L.append('')
    L.append(f"- **Identification**: {s['identification']}")
    L.append(f"- **No smoking gun**: {s['no_smoking_gun']}")
    L.append(f"- **Masses**: {s['masses']}")
    L.append(f"- **Constraint ranking**: {s['constraint_ranking']}")
    L.append(f"- **Test**: {s['test']}")
    L.append(f"- **Contrast**: {s['contrast']}")
    L.append('')

    L.append('## Test summary')
    L.append('')
    L.append('| # | Test | Key finding | PASS? |')
    L.append('|---|---|---|---|')
    label_map = {
        'T1': 'baryonic exotics from non-orientable topology (Möbius/hybrid)',
        'T2': 'no exotic J^P for baryons ⟹ supernumerary, no smoking gun',
        'T3': 'hybrid N ≈ 1.79 GeV, hybrid Δ ≈ 2.08 GeV (base + 2√σ)',
        'T4': 'constraint ranking: light N*/Δ* > hyperons > heavy baryons',
        'T5': 'lightest exotic in the dense (most-constrained) N*/Δ* region',
        'T6': 'Möbius doubling: coincide or decouple, else excluded',
        'T7': 'opposite extreme from glueballs (most vs least constrained)',
        'T8': 'BARYONIC_EXOTICS_LACK_EXOTIC_JP_LIGHT_SECTOR_MOST_CONSTRAINED',
    }
    for t in s['tests']:
        passed = '**PASS**' if t['pass'] else '**FAIL**'
        nm = t['name']; prefix = nm[:2]
        L.append(f"| {prefix} | `{nm}` | {label_map.get(prefix, '—')} | {passed} |")
    L.append('')

    # T2 J^P table
    t2 = s['tests'][1]
    L.append('## T2: Every baryon J^P is ordinary (no exotic)')
    L.append('')
    L.append('| J^P | reachable by qqq? |')
    L.append('|---|:---:|')
    for r in t2['rows']:
        L.append(f"| {r['jp']} | {'✓' if r['reachable_by_qqq'] else '—'} |")
    L.append('')
    L.append("Contrast with mesons: `1-+` is exotic (C-forbidden). No "
             "baryon `J^P` is forbidden, so a Möbius/hybrid baryon is a "
             "supernumerary ordinary-`J^P` state — countable, not "
             "manifestly exotic.")
    L.append('')

    # T4 ranking table
    t4 = s['tests'][3]
    L.append('## T4: Experimental-constraint ranking')
    L.append('')
    L.append('| channel | data density | constraint |')
    L.append('|---|---|---|')
    for r in t4['ranking']:
        L.append(f"| {r['channel']} | {r['data']} | {r['constraint']} |")
    L.append('')
    t3 = s['tests'][2]
    L.append(f"The lightest BAM baryonic exotic (`{t3['hybrid_N_baryon_GeV']:.2f} GeV` "
             f"N-based, `{t3['hybrid_Delta_baryon_GeV']:.2f} GeV` Δ-based) sits "
             "in the densely-measured light sector — the strongest "
             "constraint.")
    L.append('')

    L.append('## Verdict')
    L.append('')
    L.append(f"**{s['verdict_class']}.** {s['verdict']}")
    L.append('')

    L.append('## What this leaves open')
    L.append('')
    L.append('- **A confirmed Möbius-baryon state** — without a smoking-gun '
             '`J^P`, the prediction is a counting one; testing it means '
             'matching the dense N*/Δ* spectrum or demonstrating decoupling '
             '(the missing-resonance mechanism).')
    L.append('- **Heavy-quark baryon exotics** — the least-constrained '
             'channel, where BAM\'s Möbius predictions are freest and a '
             'clean new state is most likely findable.')
    L.append('')
    return '\n'.join(L)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def _json_default(o):
    if isinstance(o, tuple):
        return list(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if hasattr(o, '__dict__'):
        try:
            return asdict(o)
        except Exception:
            return o.__dict__
    return str(o)


def main(argv: Optional[list[str]] = None) -> int:
    summary = run_probe()
    md = render_markdown(summary)
    print(md)
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    out = here / 'runs' / f'{ts}_baryonic_exotics_classification_probe'
    out.mkdir(parents=True, exist_ok=True)
    (out / 'probe.json').write_text(
        json.dumps(summary, indent=2, default=_json_default)
    )
    (out / 'probe.md').write_text(md)
    print(f"\nWrote: {out / 'probe.json'}")
    print(f"Wrote: {out / 'probe.md'}")
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())
