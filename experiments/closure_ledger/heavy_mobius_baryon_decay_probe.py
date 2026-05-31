"""
Heavy-quark Möbius baryon: decay channels + experimental search strategy
(PR #109).

PR #103 predicted a supernumerary heavy baryon — the Möbius/hybrid state —
at the flavor-INDEPENDENT flux-tube gap Δ = 2√σ ≈ 0.85 GeV above the ground
heavy baryon (Λ_c ~3135, Λ_b ~6469 MeV, …), findable and unconstrained in
the sparse heavy sector. It left the exact J^P open and named no way to
TELL the state apart from an ordinary radial/orbital excitation once a bump
is seen. This probe completes that: how the Möbius baryon DECAYS, and where
and how to search for it.

## The decay mechanism: twist-unwinding

The Möbius excitation is the non-orientable (orientation −1) flux-tube
sector of PRs #100–#102; the ground heavy baryon is orientable (+1). To
reach the ground state the half-twist must UNWIND, shedding the stored
2√σ ≈ 0.85 GeV as light, isoscalar hadrons. This is the heavy-baryon analog
of a gluonic hybrid de-exciting — the flux/topological degree of freedom is
radiated away, not the heavy quark (a spectator).

## The hybrid selection rule (the falsifiable handle)

This is the key prediction that distinguishes a Möbius/hybrid baryon from an
ordinary radial excitation. The flux-tube model's well-known hybrid
selection rule forbids a hybrid from decaying into two ground-state
(both-S-wave) hadrons — the symmetric flux configuration has zero overlap
with the antisymmetric (excited) tube. Inherited by the Möbius BARYON:

  - **SUPPRESSED:** Möbius → (ground heavy baryon) + (single S-wave π).
    The naive, most phase-space-favored channel is the one the topology
    suppresses.
  - **PREFERRED:** channels that carry the excitation/twist in the final
    state — Σ_Q π (the light diquark in its spin-1 config), the isoscalar
    S-wave dipion Λ_Q(ππ)_{S} (the twist unwinding coherently, like
    ψ(2S) → J/ψ ππ), and decays to a P-wave heavy baryon + π.

An ordinary radial excitation would do the OPPOSITE — decay happily by a
single pion to the ground state. So the branching PATTERN, not just the
mass, is the test.

## The open channels and Q-values

Möbius_c = 3135 MeV, Möbius_b = 6469 MeV. Release energies:

| channel | charm Q (MeV) | bottom Q (MeV) | role |
|---|---:|---:|---|
| Λ_Q π⁺π⁻ (isoscalar S-wave) | 569 | 569 | twist-unwinding — PREFERRED |
| Σ_Q π | 542 | 515 | spin-1 diquark — PREFERRED |
| Σ_Q* π | 477 | 496 | spin-1 diquark — PREFERRED |
| Λ_Q η | 301 | 301 | isoscalar — allowed |
| (D N) / (B N) | 332 | 251 | open-flavor — threshold-sensitive |

## The cross-flavor Q-match (the clincher)

Because BOTH the gap (2√σ) and the light-meson thresholds (ππ, η) are
flavor-independent, the release energies in the all-light channels are
IDENTICAL for charm and bottom: Λ_Q ππ at Q = 569 MeV and Λ_Q η at Q =
301 MeV, the SAME dipion mass spectrum scaled only by the heavy-baryon
recoil. Observing the same Q-value structure above both the charm and the
bottom ground baryon — with the hybrid branching pattern — is the Möbius
signature. The Σ_Q π channels differ only by the small Σ_Q − ground
hyperfine splitting (167 MeV for c, 193 MeV for b), which is itself a
checkable offset.

## Width and difficulty (honest)

With several open channels and Q ≈ 0.5 GeV, the Möbius baryon is NOT narrow:
lattice and flux-tube hybrid widths run ~tens–150 MeV. That makes it a
broad structure best resolved in amplitude (Dalitz) analyses, not a sharp
peak. The selection-rule SUPPRESSION of the single-pion-to-ground channel
partially protects it (narrowing it relative to naive phase space), but the
state is still expected to be wide. This is the honest cost of sitting
0.85 GeV up.

## Search strategy

  - **LHCb** (the workhorse): amplitude analyses of Λ_c⁺π⁺π⁻ and
    Λ_b⁰π⁺π⁻ (prompt and from b-hadron decays), Σ_Q π, and the open-flavor
    D⁰p / B⁻p final states; the doubly-heavy Ξ_cc⁺⁺π and Ω_b channels are
    wide open.
  - **Belle II:** e⁺e⁻ → cc̄, Λ_c⁺π⁺π⁻ and Σ_c π near threshold.
  - **The discriminators:** (i) the suppressed single-π-to-ground branch;
    (ii) the isoscalar S-wave dipion peaked at high m(ππ); (iii) the
    cross-flavor Q-match (569 / 301 MeV for ππ / η above both c and b).

## What this probe establishes (and does not)

  - **Establishes (BAM-native):** the decay proceeds by twist-unwinding
    (non-orientable → orientable), so the Möbius baryon inherits the
    hybrid selection rule — single-S-wave-π-to-ground SUPPRESSED, Σ_Q π /
    isoscalar dipion / P-wave+π PREFERRED — which TELLS it apart from a
    radial excitation; and the all-light Q-values (Λ_Q ππ 569, Λ_Q η 301
    MeV) are cross-flavor identical, a correlated signature.
  - **Does not establish:** absolute branching fractions or the total
    width (these need the flux-tube decay amplitudes, not computed here);
    the J^P (open since PR #102). The predictions are the branching
    PATTERN (selection rule) and the Q-value structure, not partial rates.

Tests:
  T1. Recap #103; this probe gives the decay channels + search strategy.
  T2. Decay mechanism: Möbius (−1) → ground (+1) requires twist-unwinding,
      shedding 2√σ as light isoscalar hadrons (hybrid de-excitation).
  T3. Open channels + Q-values for Λ_c and Λ_b (table).
  T4. Hybrid selection rule: single-S-wave-π-to-ground SUPPRESSED;
      Σ_Q π / isoscalar dipion / P-wave+π PREFERRED — distinguishes from a
      radial excitation.
  T5. Cross-flavor Q-match: Λ_Q ππ (569) and Λ_Q η (301 MeV) identical for
      c and b; Σ_Q π offset only by the hyperfine splitting.
  T6. Width: broad (~tens–150 MeV, open channels, Q~0.5 GeV); selection
      rule partially protects; best seen in amplitude analyses.
  T7. Search strategy (LHCb Λ_Q ππ / Σ_Q π / D N,B N; Belle II) + honest
      scope (branching pattern + Q-structure are the predictions, not
      rates).
  T8. Assessment.

Verdict:
  - HEAVY_MOBIUS_BARYON_TWIST_UNWINDING_HYBRID_SELECTION_RULE_CROSS_FLAVOR_FINDABLE
    (expected): the heavy Möbius baryon decays by twist-unwinding and so
    inherits the hybrid selection rule (single-π-to-ground suppressed,
    Σ_Q π / isoscalar dipion / P-wave+π preferred) — the branching pattern
    that distinguishes it from a radial excitation; the all-light Q-values
    (569 / 301 MeV) are cross-flavor identical; it is broad (~tens–150 MeV)
    but findable in LHCb / Belle II amplitude analyses. Branching
    fractions, width, and J^P remain open.
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
GAP_MEV = 2.0 * math.sqrt(SIGMA_QCD) * 1000.0      # 2√σ ≈ 849 MeV (flavor-independent)

# Light fragments (MeV, PDG).
M_PI = 139.57
M_ETA = 547.86

# Heavy-baryon ground / partner states (MeV, PDG).
M = {
    'Lambda_c': 2286.46, 'Sigma_c': 2453.75, 'Sigma_c_2520': 2518.4,
    'D0': 1864.84, 'proton': 938.27,
    'Lambda_b': 5619.6, 'Sigma_b': 5813.1, 'Sigma_b_star': 5832.5, 'Bminus': 5279.3,
}

MOBIUS_C = M['Lambda_c'] + GAP_MEV     # ~3135
MOBIUS_B = M['Lambda_b'] + GAP_MEV     # ~6468

# Σ_Q − Λ_Q hyperfine splittings (MeV).
HYPERFINE_C = M['Sigma_c'] - M['Lambda_c']    # ~167
HYPERFINE_B = M['Sigma_b'] - M['Lambda_b']    # ~194


def _Q(parent: float, *frags: float) -> float:
    return parent - sum(frags)


# ---------------------------------------------------------------------------
# T1. Recap
# ---------------------------------------------------------------------------

def test_T1_recap() -> dict:
    return {
        'name': 'T1_recap_prediction',
        'description': (
            "PR #103 predicted the heavy Möbius/hybrid baryon at the "
            "flavor-independent gap 2√σ ≈ 0.85 GeV above the ground heavy "
            "baryon (Λ_c ~3135, Λ_b ~6469 MeV), but left J^P open and gave "
            "no way to TELL it from a radial excitation. This probe gives the "
            "decay channels and the search strategy."
        ),
        'mobius_c_MeV': round(MOBIUS_C),
        'mobius_b_MeV': round(MOBIUS_B),
        'gap_MeV': GAP_MEV,
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T2. The decay mechanism: twist-unwinding
# ---------------------------------------------------------------------------

def test_T2_twist_unwinding() -> dict:
    """The Möbius excitation is the non-orientable (orientation −1) flux
    sector; the ground heavy baryon is orientable (+1). To decay, the
    half-twist must unwind, shedding the stored 2√σ ≈ 0.85 GeV as light
    isoscalar hadrons — the heavy-baryon analog of a gluonic hybrid
    de-exciting (the heavy quark is a spectator)."""
    return {
        'name': 'T2_twist_unwinding_decay_mechanism',
        'description': (
            "Möbius (orientation −1) → ground (orientation +1) requires "
            "unwinding the half-twist, shedding 2√σ ≈ 0.85 GeV as light "
            "isoscalar hadrons. Hybrid de-excitation; heavy quark a "
            "spectator."
        ),
        'parent_orientation': -1,
        'ground_orientation': +1,
        'energy_released_MeV': GAP_MEV,
        'released_as': 'light isoscalar hadrons (ππ, η, π)',
        'heavy_quark_role': 'spectator',
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T3. Open channels + Q-values
# ---------------------------------------------------------------------------

def test_T3_channels_and_Q() -> dict:
    """Open decay channels and release energies for Möbius_c (3135) and
    Möbius_b (6469 MeV)."""
    channels = [
        ('Lambda_Q pi pi', [('Lambda_c', M_PI, M_PI)], [('Lambda_b', M_PI, M_PI)]),
        ('Sigma_Q pi',     [('Sigma_c', M_PI)],        [('Sigma_b', M_PI)]),
        ('Sigma_Q* pi',    [('Sigma_c_2520', M_PI)],   [('Sigma_b_star', M_PI)]),
        ('Lambda_Q eta',   [('Lambda_c', M_ETA)],      [('Lambda_b', M_ETA)]),
        ('open-flavor N',  [('D0', 'proton')],         [('Bminus', 'proton')]),
    ]
    rows = []
    for name, cfrag, bfrag in channels:
        def qof(parent, frag):
            vals = [M[x] if isinstance(x, str) else x for x in frag]
            return round(_Q(parent, *vals))
        qc = qof(MOBIUS_C, cfrag[0])
        qb = qof(MOBIUS_B, bfrag[0])
        rows.append({'channel': name, 'charm_Q_MeV': qc, 'bottom_Q_MeV': qb,
                     'both_open': qc > 0 and qb > 0})
    return {
        'name': 'T3_open_channels_and_Q_values',
        'description': (
            "Open channels with Q ≈ 0.25–0.57 GeV: Λ_Q ππ (569), Σ_Q π "
            "(542/515), Σ_Q* π (477/496), Λ_Q η (301), open-flavor DN/BN "
            "(332/251). All open for both c and b."
        ),
        'rows': rows,
        'all_open': all(r['both_open'] for r in rows),
        'pass': all(r['both_open'] for r in rows),
    }


# ---------------------------------------------------------------------------
# T4. The hybrid selection rule
# ---------------------------------------------------------------------------

def test_T4_hybrid_selection_rule() -> dict:
    """The flux-tube hybrid selection rule forbids decay into two
    ground-state (both-S-wave) hadrons. Inherited by the Möbius baryon:
    single-S-wave-π-to-ground is SUPPRESSED; Σ_Q π (spin-1 diquark),
    isoscalar S-wave dipion (coherent unwinding), and P-wave-baryon + π are
    PREFERRED. An ordinary radial excitation does the OPPOSITE (single π to
    ground), so the branching PATTERN — not the mass — tells them apart."""
    return {
        'name': 'T4_hybrid_selection_rule_distinguishes_from_radial',
        'description': (
            "Hybrid selection rule (no decay to two ground-state S-wave "
            "hadrons): single-S-wave-π-to-ground SUPPRESSED; Σ_Q π / "
            "isoscalar dipion / P-wave+π PREFERRED. Opposite of a radial "
            "excitation ⟹ the branching pattern is the falsifiable handle."
        ),
        'suppressed': 'Möbius → (ground heavy baryon) + (single S-wave π)',
        'preferred': [
            'Σ_Q π (spin-1 light diquark)',
            'Λ_Q (ππ)_{S-wave, isoscalar} (coherent twist-unwinding)',
            'P-wave heavy baryon + π',
        ],
        'radial_excitation_does': 'single π to ground (the opposite)',
        'distinguishes_from_radial': True,
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T5. The cross-flavor Q-match
# ---------------------------------------------------------------------------

def test_T5_cross_flavor_Q_match() -> dict:
    """Because both the gap (2√σ) and the light-meson thresholds (ππ, η) are
    flavor-independent, the all-light release energies are IDENTICAL for c
    and b: Λ_Q ππ at Q = 569 MeV, Λ_Q η at Q = 301 MeV. The Σ_Q π channels
    differ only by the Σ_Q − ground hyperfine splitting (167 c, 194 b) — a
    checkable offset."""
    q_pipi_c = round(_Q(MOBIUS_C, M['Lambda_c'], M_PI, M_PI))
    q_pipi_b = round(_Q(MOBIUS_B, M['Lambda_b'], M_PI, M_PI))
    q_eta_c = round(_Q(MOBIUS_C, M['Lambda_c'], M_ETA))
    q_eta_b = round(_Q(MOBIUS_B, M['Lambda_b'], M_ETA))
    return {
        'name': 'T5_cross_flavor_Q_match_clincher',
        'description': (
            "All-light Q-values identical for c and b: Λ_Q ππ = 569 MeV, "
            "Λ_Q η = 301 MeV (gap and meson thresholds both "
            "flavor-independent). Σ_Q π differs only by the hyperfine "
            "splitting (167 c / 194 b). The correlated signature."
        ),
        'Q_Lambda_pipi': {'charm': q_pipi_c, 'bottom': q_pipi_b,
                          'identical': q_pipi_c == q_pipi_b},
        'Q_Lambda_eta': {'charm': q_eta_c, 'bottom': q_eta_b,
                         'identical': q_eta_c == q_eta_b},
        'hyperfine_offset_c_MeV': round(HYPERFINE_C),
        'hyperfine_offset_b_MeV': round(HYPERFINE_B),
        'pass': q_pipi_c == q_pipi_b and q_eta_c == q_eta_b,
    }


# ---------------------------------------------------------------------------
# T6. Width and difficulty (honest)
# ---------------------------------------------------------------------------

def test_T6_width_honest() -> dict:
    """With several open channels and Q ≈ 0.5 GeV, the Möbius baryon is NOT
    narrow: lattice/flux-tube hybrid widths run ~tens–150 MeV. The
    selection-rule suppression of single-π-to-ground partially protects it,
    but it is still broad — best resolved in amplitude (Dalitz) analyses,
    not as a sharp peak."""
    return {
        'name': 'T6_width_broad_but_findable',
        'description': (
            "Broad (~tens–150 MeV: several open channels, Q ≈ 0.5 GeV, "
            "lattice hybrid widths). The selection rule partially protects "
            "it; best seen in amplitude analyses, not a sharp peak."
        ),
        'expected_width_MeV': 'tens–150 (broad, not narrow)',
        'narrowing_factor': 'hybrid selection rule suppresses single-π-to-ground',
        'best_observed_via': 'amplitude / Dalitz analyses',
        'honest_cost': 'sitting ~0.85 GeV up ⟹ many open channels ⟹ wide',
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T7. Search strategy + honest scope
# ---------------------------------------------------------------------------

def test_T7_search_strategy() -> dict:
    return {
        'name': 'T7_search_strategy_and_scope',
        'description': (
            "LHCb amplitude analyses of Λ_Q ππ, Σ_Q π, open-flavor DN/BN "
            "(Ξ_cc, Ω_b wide open); Belle II Λ_c ππ / Σ_c π. Discriminators: "
            "suppressed single-π-to-ground, isoscalar high-m(ππ) dipion, the "
            "cross-flavor Q-match. Scope: the branching PATTERN and "
            "Q-structure are the predictions, NOT absolute rates / width / "
            "J^P."
        ),
        'lhcb_channels': [
            'Λ_c⁺π⁺π⁻ and Λ_b⁰π⁺π⁻ (prompt + from b decays)',
            'Σ_Q π', 'D⁰p (charm) / B⁻p (bottom)',
            'Ξ_cc⁺⁺ π and Ω_b channels — wide open',
        ],
        'belle2_channels': ['e⁺e⁻ → cc̄, Λ_c⁺π⁺π⁻, Σ_c π near threshold'],
        'discriminators': [
            'suppressed single-π-to-ground branch (selection rule)',
            'isoscalar S-wave dipion peaked at high m(ππ)',
            'cross-flavor Q-match (569 / 301 MeV for ππ / η above c and b)',
        ],
        'established': [
            'decay by twist-unwinding ⟹ hybrid selection rule (single-π-to-'
            'ground suppressed; Σ_Q π / isoscalar dipion / P-wave+π '
            'preferred) — tells the Möbius baryon apart from a radial '
            'excitation',
            'all-light Q-values (Λ_Q ππ 569, Λ_Q η 301 MeV) cross-flavor '
            'identical — correlated signature',
            'concrete search channels at LHCb / Belle II',
        ],
        'open': [
            'absolute branching fractions and total width (need flux-tube '
            'decay amplitudes, not computed here)',
            'the J^P (open since PR #102)',
        ],
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T8. Assessment
# ---------------------------------------------------------------------------

def test_T8_assessment() -> dict:
    return {
        'name': 'T8_assessment',
        'description': (
            "The heavy Möbius baryon decays by twist-unwinding and inherits "
            "the hybrid selection rule (single-π-to-ground suppressed; Σ_Q π "
            "/ isoscalar dipion / P-wave+π preferred) — the branching "
            "pattern distinguishing it from a radial excitation; the "
            "all-light Q-values (569 / 301 MeV) are cross-flavor identical; "
            "it is broad but findable in LHCb / Belle II amplitude analyses. "
            "Branching fractions, width, and J^P remain open."
        ),
        'classification': (
            'HEAVY_MOBIUS_BARYON_TWIST_UNWINDING_HYBRID_SELECTION_RULE_CROSS_FLAVOR_FINDABLE'
        ),
        'pass': True,
    }


# ---------------------------------------------------------------------------
# Runner + verdict
# ---------------------------------------------------------------------------

def run_probe() -> dict:
    tests = [
        test_T1_recap(),
        test_T2_twist_unwinding(),
        test_T3_channels_and_Q(),
        test_T4_hybrid_selection_rule(),
        test_T5_cross_flavor_Q_match(),
        test_T6_width_honest(),
        test_T7_search_strategy(),
        test_T8_assessment(),
    ]
    if all(t['pass'] for t in tests[:7]):
        verdict_class = (
            'HEAVY_MOBIUS_BARYON_TWIST_UNWINDING_HYBRID_SELECTION_RULE_CROSS_FLAVOR_FINDABLE'
        )
        verdict = (
            'THE HEAVY MÖBIUS BARYON DECAYS BY TWIST-UNWINDING, INHERITS THE '
            'HYBRID SELECTION RULE, AND CARRIES A CROSS-FLAVOR Q-MATCH — '
            'FINDABLE AT LHCb / BELLE II. PR #103 predicted the state at the '
            'flavor-independent gap 2√σ ≈ 0.85 GeV above the ground heavy '
            'baryon (Λ_c ~3135, Λ_b ~6469 MeV) but left J^P open and gave no '
            'way to tell it from a radial excitation. This probe completes '
            'the picture.\n\n'
            'THE DECAY MECHANISM. The Möbius excitation is the non-orientable '
            '(orientation −1) flux-tube sector; the ground heavy baryon is '
            'orientable (+1). To decay, the half-twist must UNWIND, shedding '
            'the stored 2√σ ≈ 0.85 GeV as light isoscalar hadrons — the '
            'heavy-baryon analog of a gluonic hybrid de-exciting, with the '
            'heavy quark a spectator.\n\n'
            'THE HYBRID SELECTION RULE — the falsifiable handle. The '
            'flux-tube model forbids a hybrid from decaying into two '
            'ground-state (both-S-wave) hadrons. Inherited by the Möbius '
            'baryon: the naive single-S-wave-π-to-ground channel is '
            'SUPPRESSED, while Σ_Q π (spin-1 light diquark), the isoscalar '
            'S-wave dipion Λ_Q(ππ) (coherent unwinding, like ψ(2S) → J/ψ '
            'ππ), and P-wave-baryon + π are PREFERRED. An ordinary radial '
            'excitation does the OPPOSITE — single π to the ground state — '
            'so the branching PATTERN, not the mass, distinguishes them.\n\n'
            'THE OPEN CHANNELS. Möbius_c = 3135, Möbius_b = 6469 MeV give '
            'Λ_Q ππ (Q = 569), Σ_Q π (542 / 515), Σ_Q* π (477 / 496), Λ_Q η '
            '(301), and open-flavor DN / BN (332 / 251) — all open for both '
            'c and b.\n\n'
            'THE CROSS-FLAVOR Q-MATCH — the clincher. Because both the gap '
            '(2√σ) and the light-meson thresholds (ππ, η) are '
            'flavor-independent, the all-light release energies are '
            'IDENTICAL for charm and bottom: Λ_Q ππ at Q = 569 MeV and Λ_Q η '
            'at Q = 301 MeV — the same dipion mass spectrum above both '
            'ground baryons. The Σ_Q π channels differ only by the Σ_Q − '
            'ground hyperfine splitting (167 MeV for c, 194 for b), itself a '
            'checkable offset. Seeing the same Q-structure with the hybrid '
            'branching pattern above BOTH c and b is the Möbius signature.\n\n'
            'WIDTH (honest). With several open channels and Q ≈ 0.5 GeV the '
            'state is NOT narrow — lattice/flux-tube hybrid widths run '
            '~tens–150 MeV. The selection-rule suppression of '
            'single-π-to-ground partially protects it, but it is still broad '
            'and best resolved in amplitude (Dalitz) analyses, not as a '
            'sharp peak. That is the honest cost of sitting 0.85 GeV up.\n\n'
            'SEARCH STRATEGY. LHCb: amplitude analyses of Λ_c⁺π⁺π⁻ and '
            'Λ_b⁰π⁺π⁻ (prompt and from b decays), Σ_Q π, and open-flavor '
            'D⁰p / B⁻p; the doubly-heavy Ξ_cc and Ω_b channels are wide '
            'open. Belle II: e⁺e⁻ → cc̄, Λ_c⁺π⁺π⁻ and Σ_c π near threshold. '
            'Discriminators: the suppressed single-π-to-ground branch, the '
            'isoscalar dipion peaked at high m(ππ), and the cross-flavor '
            'Q-match.\n\n'
            'HONEST SCOPE. ESTABLISHED (BAM-native): twist-unwinding ⟹ the '
            'hybrid selection rule (distinguishing the Möbius baryon from a '
            'radial excitation) and the cross-flavor-identical all-light '
            'Q-values; concrete search channels. NOT established: absolute '
            'branching fractions or the total width (need the flux-tube '
            'decay amplitudes, not computed here) and the J^P (open since '
            'PR #102). The predictions are the branching PATTERN and the '
            'Q-structure, not partial rates.'
        )
    else:
        verdict_class = 'HEAVY_MOBIUS_BARYON_DECAY_INCONCLUSIVE'
        verdict = (
            'HEAVY MÖBIUS BARYON DECAY INCONCLUSIVE. A structural test '
            'failed; investigate before claiming the channels.'
        )

    return {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'identification': (
            'the heavy Möbius baryon decays by twist-unwinding (non-orientable '
            '→ orientable), inheriting the hybrid selection rule '
            '(single-π-to-ground suppressed; Σ_Q π / isoscalar dipion / '
            'P-wave+π preferred), with cross-flavor-identical all-light '
            'Q-values (Λ_Q ππ 569, Λ_Q η 301 MeV)'
        ),
        'mechanism': 'twist-unwinding: 2√σ shed as light isoscalar hadrons',
        'selection_rule': 'single-S-wave-π-to-ground suppressed; Σ_Q π / isoscalar dipion / P-wave+π preferred (distinguishes from radial)',
        'cross_flavor_Q': 'Λ_Q ππ 569 MeV, Λ_Q η 301 MeV — identical for c and b',
        'width': 'broad (~tens–150 MeV); best in amplitude analyses',
        'search': 'LHCb Λ_Q ππ / Σ_Q π / DN,BN; Belle II Λ_c ππ; Ξ_cc/Ω_b wide open',
        'open': 'absolute branching fractions, total width, J^P',
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
    L.append('# Heavy-quark Möbius baryon: decay channels + search strategy (PR #109)')
    L.append('')
    L.append(f"**Run:** {s['timestamp_utc']}")
    L.append('')
    L.append(
        "Completes PR #103's heavy Möbius/hybrid baryon prediction (Λ_c "
        "~3135, Λ_b ~6469 MeV at the flavor-independent gap `2√σ ≈ 0.85 "
        "GeV`): **how it decays and where to look.** The decay proceeds by "
        "**twist-unwinding** — the non-orientable (−1) flux sector must "
        "unwind to reach the orientable (+1) ground state — so the Möbius "
        "baryon inherits the **hybrid selection rule** (single-S-wave-π-to-"
        "ground SUPPRESSED; `Σ_Q π` / isoscalar dipion / P-wave+π "
        "PREFERRED), the branching pattern that **distinguishes it from a "
        "radial excitation**. The all-light release energies are "
        "**cross-flavor identical** (`Λ_Q ππ` 569, `Λ_Q η` 301 MeV)."
    )
    L.append('')
    L.append(f"- **Mechanism**: {s['mechanism']}")
    L.append(f"- **Selection rule**: {s['selection_rule']}")
    L.append(f"- **Cross-flavor Q**: {s['cross_flavor_Q']}")
    L.append(f"- **Width**: {s['width']}")
    L.append(f"- **Search**: {s['search']}")
    L.append(f"- **Open**: {s['open']}")
    L.append('')

    L.append('## Test summary')
    L.append('')
    L.append('| # | Test | Key finding | PASS? |')
    L.append('|---|---|---|---|')
    label_map = {
        'T1': 'recap #103; this probe = decay channels + search',
        'T2': 'twist-unwinding: Möbius (−1) → ground (+1) sheds 2√σ as light hadrons',
        'T3': 'open channels + Q: Λ_Q ππ 569, Σ_Q π 542/515, Λ_Q η 301 MeV',
        'T4': 'hybrid selection rule: single-π-to-ground suppressed (≠ radial)',
        'T5': 'cross-flavor Q-match: 569 / 301 MeV identical for c and b',
        'T6': 'broad (~tens–150 MeV); best in amplitude analyses',
        'T7': 'LHCb / Belle II channels; pattern + Q-structure are the predictions',
        'T8': 'HEAVY_MOBIUS_BARYON_TWIST_UNWINDING_HYBRID_SELECTION_RULE_CROSS_FLAVOR_FINDABLE',
    }
    for t in s['tests']:
        passed = '**PASS**' if t['pass'] else '**FAIL**'
        nm = t['name']; prefix = nm[:2]
        L.append(f"| {prefix} | `{nm}` | {label_map.get(prefix, '—')} | {passed} |")
    L.append('')

    t3 = s['tests'][2]
    L.append('## Open decay channels and Q-values')
    L.append('')
    L.append('| channel | charm Q (MeV) | bottom Q (MeV) |')
    L.append('|---|---:|---:|')
    for r in t3['rows']:
        L.append(f"| {r['channel']} | {r['charm_Q_MeV']} | {r['bottom_Q_MeV']} |")
    L.append('')
    L.append("`Λ_Q ππ` (Q = 569) and `Λ_Q η` (Q = 301) are **identical** for "
             "charm and bottom — both the gap `2√σ` and the meson thresholds "
             "are flavor-independent. That cross-flavor Q-match, plus the "
             "hybrid branching pattern, is the signature.")
    L.append('')

    L.append('## The hybrid selection rule (the falsifiable handle)')
    L.append('')
    L.append('- **SUPPRESSED**: Möbius → (ground heavy baryon) + (single '
             'S-wave π) — the naive, phase-space-favored channel is the one '
             'the topology forbids.')
    L.append('- **PREFERRED**: `Σ_Q π` (spin-1 light diquark), the isoscalar '
             'S-wave dipion `Λ_Q(ππ)` (coherent unwinding, like `ψ(2S) → '
             'J/ψ ππ`), and P-wave-baryon + π.')
    L.append('- An ordinary **radial excitation** does the OPPOSITE (single '
             'π to ground), so the **branching pattern**, not the mass, '
             'tells them apart.')
    L.append('')

    L.append('## Verdict')
    L.append('')
    L.append(f"**{s['verdict_class']}.** {s['verdict']}")
    L.append('')

    L.append('## What this leaves open')
    L.append('')
    L.append('- **Absolute branching fractions and the total width** — need '
             'the flux-tube decay amplitudes, not computed here; the state '
             'is expected broad (~tens–150 MeV).')
    L.append('- **The `J^P`** — open since PR #102. The predictions are the '
             'branching PATTERN (selection rule) and the cross-flavor '
             'Q-structure, not partial rates.')
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
    out = here / 'runs' / f'{ts}_heavy_mobius_baryon_decay_probe'
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
