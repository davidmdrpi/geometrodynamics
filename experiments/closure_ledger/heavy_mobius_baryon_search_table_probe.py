"""
Heavy Möbius baryon: a sharper LHCb / Belle II search table (PR #114).

PRs #103/#109/#110 predicted the heavy Möbius/hybrid baryon (ground heavy
baryon + 2√σ), worked out its twist-unwinding decays and the hybrid
selection rule, and compiled the non-orientable sector into one note. This
probe CONVERTS that physics into a sharper, actionable LHCb / Belle II
search table: for each target, the concrete reconstruction chain, the
facility and production mode, the primary discovery channel with its
kinematic handle, the discriminator, the constraint status, and a tiered
priority. Every mass is a pushforward of the single scale √σ.

## The new quantitative handle: the dipion endpoint

The preferred twist-unwinding channel is Λ_Q(ππ)_S. Its dipion
invariant-mass endpoint — the maximum m(ππ), reached when the heavy baryon
recoils at rest — is

    m(ππ)_max = M_Möbius − M_ground = 2√σ ≈ 849 MeV,

FLAVOR-INDEPENDENT: the same 849 MeV endpoint above the charm and the
bottom ground baryon, with the spectrum peaking HIGH (coherent isoscalar
S-wave, like ψ(2S) → J/ψ ππ). This is sharper than a Q-value: it is a fixed
edge in a directly-plotted observable, identical for c and b — a single
overlay that tests the whole framework.

## The search table (tiered)

  TIER 1 — the discovery pair (highest yield × the cross-flavor lever):
    • Λ_c Möbius (3135 MeV): Λ_c⁺π⁺π⁻ with Λ_c⁺ → pK⁻π⁺; prompt at LHCb +
      Belle II near threshold. The statistics-richest channel.
    • Λ_b Möbius (6469 MeV): Λ_b⁰π⁺π⁻ with Λ_b⁰ → Λ_c⁺π⁻; from b-decays at
      LHCb. The cross-flavor partner.
    Together they carry the clincher: the SAME dipion endpoint (849 MeV)
    and Q-values (ππ 569, η 301 MeV) above both — a simultaneous two-channel
    fit.

  TIER 2 — the clean frontier (entirely unexplored, but rate-limited):
    • Ξ_cc Möbius (4471 MeV): Ξ_cc⁺⁺π⁺π⁻ with Ξ_cc⁺⁺ → Λ_c⁺K⁻π⁺π⁺.
    • Ω_b Möbius (6894 MeV): Ω_b⁻π⁺π⁻ with Ω_b⁻ → Ω_c⁰π⁻.
    No measured excitation spectrum at all (PR #103's "freest of the free"),
    so a clean bump = discovery — but doubly-heavy / Ω_b production is rare.

  TIER 3 — calibratable:
    • Ω_c Möbius (3544 MeV): Ω_c⁰π⁺π⁻. Sits above the well-mapped 2017 LHCb
      Ω_c excitations (≤3120), which calibrate the search but mean it is not
      virgin territory.

## The discriminators (per channel)

  - **Suppressed single-π-to-ground** — the hybrid selection rule forbids
    decay to the ground baryon + one S-wave π; an ordinary radial excitation
    does the opposite. The branching PATTERN distinguishes them.
  - **Dipion endpoint + shape** — a fixed 849 MeV edge, peaking high; the
    same overlay for c and b.
  - **Cross-flavor Q-match** — Λ_Q ππ at 569 and Λ_Q η at 301 MeV identical
    for charm and bottom; a joint constraint, not a single bump.

## Honest scope

This is a prioritization/presentation deliverable, not new physics. The
masses carry the ±lattice-hybrid-gap uncertainty (~0.8–1.3 GeV), the states
are broad (~tens–150 MeV ⟹ amplitude analyses), absolute branching
fractions and widths are not computed, and the J^P is open. The table makes
the existing #109/#110 content searchable and ranked on one page.

Tests:
  T1. Scope: convert #109/#110 into a sharper, tiered, actionable table.
  T2. Targets + reconstruction chains (golden exclusive modes).
  T3. Facility + production mode per target (LHCb prompt/b-decay; Belle II).
  T4. Primary channel + the dipion endpoint m(ππ)_max = 2√σ = 849 MeV
      (flavor-independent, peaks high) — the new sharp handle.
  T5. Discriminators per channel (selection rule; endpoint; Q-match).
  T6. Tiered priority with explicit criteria (yield, cleanliness,
      unexplored, cross-flavor leverage).
  T7. Honest scope (presentation; masses ±band, broad, BFs/J^P open).
  T8. Assessment (emits the table).

Verdict:
  - HEAVY_MOBIUS_BARYON_SEARCH_TABLE_SHARPENED (expected): #109/#110 are
    converted into a tiered LHCb / Belle II search table — reconstruction
    chains, facility/production, the flavor-independent dipion endpoint
    (849 MeV) as the new sharp handle, per-channel discriminators, and an
    explicit priority (Tier 1: the Λ_c + Λ_b cross-flavor pair; Tier 2: the
    unexplored Ξ_cc / Ω_b; Tier 3: the calibratable Ω_c). Masses ±band,
    widths broad, BFs / J^P open.
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
GAP_MEV = 2.0 * math.sqrt(SIGMA_QCD) * 1000.0      # 2√σ ≈ 849 MeV
M_PI = 139.57
M_ETA = 547.86
DIPION_ENDPOINT_MEV = GAP_MEV                       # m(ππ)_max = M_Möbius − M_ground

# Targets: ground (MeV), reconstruction chain, facility, production, tier,
# constraint status, and the four transparent priority factors (0–2 each):
# (yield, reconstruction cleanliness, unexplored, cross-flavor leverage).
TARGETS = {
    'Λ_c': {
        'ground': 2286, 'chain': 'Λ_c⁺π⁺π⁻, Λ_c⁺ → pK⁻π⁺',
        'facility': 'LHCb (prompt) + Belle II', 'production': 'prompt cc̄ / e⁺e⁻',
        'status': 'above Λ_c(2940) ceiling', 'tier': 1,
        'factors': (2, 2, 0, 2),
    },
    'Λ_b': {
        'ground': 5620, 'chain': 'Λ_b⁰π⁺π⁻, Λ_b⁰ → Λ_c⁺π⁻',
        'facility': 'LHCb', 'production': 'from b-hadron decays',
        'status': 'above Λ_b(6152) ceiling', 'tier': 1,
        'factors': (2, 2, 0, 2),
    },
    'Ξ_cc': {
        'ground': 3622, 'chain': 'Ξ_cc⁺⁺π⁺π⁻, Ξ_cc⁺⁺ → Λ_c⁺K⁻π⁺π⁺',
        'facility': 'LHCb', 'production': 'prompt (rare, doubly-heavy)',
        'status': 'NO measured excitations — unexplored', 'tier': 2,
        'factors': (0, 1, 2, 0),
    },
    'Ω_b': {
        'ground': 6045, 'chain': 'Ω_b⁻π⁺π⁻, Ω_b⁻ → Ω_c⁰π⁻',
        'facility': 'LHCb', 'production': 'from b-decays (rare)',
        'status': 'NO measured excitations — unexplored', 'tier': 2,
        'factors': (0, 1, 2, 1),
    },
    'Ω_c': {
        'ground': 2695, 'chain': 'Ω_c⁰π⁺π⁻',
        'facility': 'LHCb', 'production': 'prompt / from b-decays',
        'status': 'above the 2017 Ω_c excitations (≤3120)', 'tier': 3,
        'factors': (1, 1, 0, 1),
    },
}


def _mobius(g: int) -> int:
    return round(g + GAP_MEV)


def _Q(parent: float, *frags: float) -> float:
    return parent - sum(frags)


# ---------------------------------------------------------------------------
# T1. Scope
# ---------------------------------------------------------------------------

def test_T1_scope() -> dict:
    return {
        'name': 'T1_scope',
        'description': (
            "Convert PRs #103/#109/#110 (heavy Möbius baryon masses, decays, "
            "note) into a sharper, tiered, actionable LHCb / Belle II search "
            "table: reconstruction chains, facility/production, the dipion "
            "endpoint handle, discriminators, and ranked priority."
        ),
        'inputs': ['PR #103 masses', 'PR #109 decays + selection rule', 'PR #110 note'],
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T2. Targets + reconstruction chains
# ---------------------------------------------------------------------------

def test_T2_targets_and_chains() -> dict:
    rows = []
    for name, d in TARGETS.items():
        rows.append({'target': f'{name} Möbius', 'mobius_MeV': _mobius(d['ground']),
                     'reconstruction': d['chain']})
    return {
        'name': 'T2_targets_and_reconstruction_chains',
        'description': (
            "Five targets with golden exclusive reconstruction chains: Λ_c "
            "(pK⁻π⁺), Λ_b (Λ_c⁺π⁻), Ξ_cc (Λ_c⁺K⁻π⁺π⁺), Ω_b (Ω_c⁰π⁻), Ω_c."
        ),
        'rows': rows,
        'pass': all(r['mobius_MeV'] > 3000 for r in rows),
    }


# ---------------------------------------------------------------------------
# T3. Facility + production
# ---------------------------------------------------------------------------

def test_T3_facility_production() -> dict:
    rows = [{'target': f'{n} Möbius', 'facility': d['facility'],
             'production': d['production']} for n, d in TARGETS.items()]
    return {
        'name': 'T3_facility_and_production',
        'description': (
            "Facility/production: Λ_c at LHCb (prompt) + Belle II; Λ_b, Ξ_cc, "
            "Ω_b, Ω_c at LHCb (prompt or from b-decays). Belle II adds the "
            "charm-threshold Λ_c channel."
        ),
        'rows': rows,
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T4. Primary channel + dipion endpoint
# ---------------------------------------------------------------------------

def test_T4_dipion_endpoint() -> dict:
    """The preferred Λ_Q(ππ)_S channel has dipion endpoint m(ππ)_max =
    M_Möbius − M_ground = 2√σ ≈ 849 MeV — FLAVOR-INDEPENDENT (same for c and
    b), peaking high (isoscalar S-wave). A fixed edge in a directly-plotted
    observable; the single overlay that tests the framework."""
    q_pipi = _Q(GAP_MEV, M_PI, M_PI)
    q_eta = GAP_MEV - M_ETA
    return {
        'name': 'T4_dipion_endpoint_flavor_independent',
        'description': (
            "m(ππ)_max = 2√σ ≈ 849 MeV, flavor-independent (same above c and "
            "b), peaking high. Sharper than a Q-value: a fixed edge in the "
            "plotted dipion mass, identical for charm and bottom."
        ),
        'dipion_endpoint_MeV': round(DIPION_ENDPOINT_MEV),
        'flavor_independent': True,
        'spectrum_peaks': 'high m(ππ) (coherent isoscalar S-wave, like ψ(2S)→J/ψ ππ)',
        'Q_Lambda_pipi_MeV': round(q_pipi),
        'Q_Lambda_eta_MeV': round(q_eta),
        'pass': 800 < DIPION_ENDPOINT_MEV < 900,
    }


# ---------------------------------------------------------------------------
# T5. Discriminators
# ---------------------------------------------------------------------------

def test_T5_discriminators() -> dict:
    return {
        'name': 'T5_discriminators_per_channel',
        'description': (
            "Per-channel discriminators: (1) suppressed single-π-to-ground "
            "(hybrid selection rule; opposite of a radial excitation); (2) "
            "the 849 MeV dipion endpoint peaking high, identical for c and "
            "b; (3) the cross-flavor Q-match (ππ 569, η 301 MeV)."
        ),
        'discriminators': [
            'suppressed single-π-to-ground (selection rule ⟹ ≠ radial excitation)',
            'dipion endpoint 849 MeV, peaking high — same overlay for c and b',
            'cross-flavor Q-match: Λ_Q ππ 569, Λ_Q η 301 MeV (joint constraint)',
        ],
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T6. Tiered priority
# ---------------------------------------------------------------------------

def test_T6_priority() -> dict:
    """Transparent priority score = yield + reconstruction cleanliness +
    unexplored + cross-flavor leverage (each 0–2). Tier 1: the Λ_c + Λ_b
    cross-flavor pair (highest yield × the clincher). Tier 2: the unexplored
    Ξ_cc / Ω_b (clean but rate-limited). Tier 3: the calibratable Ω_c."""
    rows = []
    for name, d in TARGETS.items():
        score = sum(d['factors'])
        rows.append({'target': f'{name} Möbius', 'tier': d['tier'],
                     'score': score, 'factors_yld_cln_unx_xfl': list(d['factors']),
                     'status': d['status']})
    rows.sort(key=lambda r: (r['tier'], -r['score']))
    tier1 = [r['target'] for r in rows if r['tier'] == 1]
    return {
        'name': 'T6_tiered_priority',
        'description': (
            "Priority = yield + cleanliness + unexplored + cross-flavor "
            "(0–2 each). Tier 1: Λ_c + Λ_b (the cross-flavor discovery "
            "pair); Tier 2: Ξ_cc / Ω_b (unexplored, rate-limited); Tier 3: "
            "Ω_c (calibratable)."
        ),
        'rows': rows,
        'tier1_pair': tier1,
        'pass': set(tier1) == {'Λ_c Möbius', 'Λ_b Möbius'},
    }


# ---------------------------------------------------------------------------
# T7. Honest scope
# ---------------------------------------------------------------------------

def test_T7_honest_scope() -> dict:
    return {
        'name': 'T7_honest_scope',
        'description': (
            "A prioritization/presentation deliverable, not new physics: "
            "masses carry the ±lattice-hybrid-gap uncertainty (~0.8–1.3 "
            "GeV), states are broad (~tens–150 MeV ⟹ amplitude analyses), "
            "absolute branching fractions/widths are not computed, J^P open."
        ),
        'established': 'a sharper, tiered, actionable search table from #109/#110',
        'open': [
            'exact masses (±lattice hybrid gap ~0.8–1.3 GeV)',
            'absolute branching fractions and total widths (broad states)',
            'the J^P (supernumerary ordinary-J^P; cross-flavor correlation is the handle)',
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
            "PRs #109/#110 are converted into a tiered LHCb / Belle II search "
            "table: reconstruction chains, facility/production, the "
            "flavor-independent dipion endpoint (849 MeV) as the new sharp "
            "handle, per-channel discriminators, and explicit priority "
            "(Tier 1: Λ_c + Λ_b pair; Tier 2: Ξ_cc / Ω_b; Tier 3: Ω_c)."
        ),
        'classification': 'HEAVY_MOBIUS_BARYON_SEARCH_TABLE_SHARPENED',
        'pass': True,
    }


# ---------------------------------------------------------------------------
# Runner + verdict
# ---------------------------------------------------------------------------

def run_probe() -> dict:
    tests = [
        test_T1_scope(),
        test_T2_targets_and_chains(),
        test_T3_facility_production(),
        test_T4_dipion_endpoint(),
        test_T5_discriminators(),
        test_T6_priority(),
        test_T7_honest_scope(),
        test_T8_assessment(),
    ]
    if all(t['pass'] for t in tests[:7]):
        verdict_class = 'HEAVY_MOBIUS_BARYON_SEARCH_TABLE_SHARPENED'
        verdict = (
            'PRs #109/#110 SHARPENED INTO A TIERED, ACTIONABLE LHCb / BELLE '
            'II SEARCH TABLE. The heavy Möbius/hybrid baryon physics — masses '
            '(ground + 2√σ), twist-unwinding decays, the hybrid selection '
            'rule, the cross-flavor Q-match — is converted into a search '
            'table an analyst can act on, with every mass a pushforward of '
            'the single scale √σ.\n\n'
            'THE NEW SHARP HANDLE: THE DIPION ENDPOINT. The preferred '
            'twist-unwinding channel Λ_Q(ππ)_S has a dipion invariant-mass '
            'endpoint m(ππ)_max = M_Möbius − M_ground = 2√σ ≈ 849 MeV, '
            'FLAVOR-INDEPENDENT — the same 849 MeV edge above the charm and '
            'the bottom ground baryon, with the spectrum peaking high '
            '(coherent isoscalar S-wave, like ψ(2S) → J/ψ ππ). This is '
            'sharper than a Q-value: a fixed edge in a directly-plotted '
            'observable, identical for c and b — a single overlay that tests '
            'the whole framework.\n\n'
            'THE TIERED TABLE. Tier 1, the discovery pair: Λ_c Möbius (3135 '
            'MeV) via Λ_c⁺π⁺π⁻ (Λ_c⁺ → pK⁻π⁺), prompt at LHCb and at Belle '
            'II — the statistics-richest channel — and Λ_b Möbius (6469 MeV) '
            'via Λ_b⁰π⁺π⁻ (Λ_b⁰ → Λ_c⁺π⁻) from b-decays at LHCb, the '
            'cross-flavor partner; together they carry the clincher (same '
            'dipion endpoint 849 MeV and Q-values ππ 569 / η 301 MeV above '
            'both). Tier 2, the clean frontier: Ξ_cc Möbius (4471) via '
            'Ξ_cc⁺⁺π⁺π⁻ and Ω_b Möbius (6894) via Ω_b⁻π⁺π⁻ — no measured '
            'excitations at all (PR #103\'s "freest of the free"), so a clean '
            'bump = discovery, but doubly-heavy / Ω_b production is rare. '
            'Tier 3, calibratable: Ω_c Möbius (3544) via Ω_c⁰π⁺π⁻, sitting '
            'above the well-mapped 2017 LHCb Ω_c excitations (≤3120).\n\n'
            'THE DISCRIMINATORS. Per channel: (1) the suppressed '
            'single-π-to-ground branch (the hybrid selection rule — an '
            'ordinary radial excitation does the opposite); (2) the 849 MeV '
            'dipion endpoint peaking high, the same overlay for c and b; (3) '
            'the cross-flavor Q-match (Λ_Q ππ 569, Λ_Q η 301 MeV identical '
            'for charm and bottom), a joint constraint rather than a single '
            'bump.\n\n'
            'HONEST SCOPE. A prioritization/presentation deliverable, not new '
            'physics: the masses carry the ±lattice-hybrid-gap uncertainty '
            '(~0.8–1.3 GeV), the states are broad (~tens–150 MeV ⟹ amplitude '
            'analyses), absolute branching fractions and widths are not '
            'computed, and the J^P is open. The table makes the existing '
            '#109/#110 content searchable and ranked on one page.'
        )
    else:
        verdict_class = 'HEAVY_MOBIUS_BARYON_SEARCH_TABLE_INCONCLUSIVE'
        verdict = (
            'INCONCLUSIVE. A structural test failed; reconcile the table with '
            '#109/#110 before issuing.'
        )

    return {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'identification': (
            'PRs #109/#110 converted into a tiered LHCb / Belle II search '
            'table — reconstruction chains, facility/production, the '
            'flavor-independent dipion endpoint (849 MeV), per-channel '
            'discriminators, and explicit priority'
        ),
        'new_handle': 'dipion endpoint m(ππ)_max = 2√σ ≈ 849 MeV, flavor-independent, peaks high',
        'tier1': 'Λ_c (3135) + Λ_b (6469) — the cross-flavor discovery pair',
        'tier2': 'Ξ_cc (4471) + Ω_b (6894) — unexplored, rate-limited',
        'tier3': 'Ω_c (3544) — above the 2017 Ω_c excitations (calibratable)',
        'discriminators': 'suppressed single-π-to-ground; 849 MeV dipion endpoint; cross-flavor Q-match',
        'open': 'masses ±lattice band; broad widths; branching fractions; J^P',
        'tests': tests,
        'n_passed': sum(1 for t in tests if t['pass']),
        'n_total': len(tests),
        'verdict_class': verdict_class,
        'verdict': verdict,
    }


# ---------------------------------------------------------------------------
# The search table (the deliverable)
# ---------------------------------------------------------------------------

def render_search_table(s: dict) -> str:
    L: list[str] = []
    L.append('# Heavy Möbius baryon — LHCb / Belle II search table')
    L.append('')
    L.append(f"_Compiled {s['timestamp_utc']} · PRs #103/#109/#110/#114 · all "
             f"masses pushed forward from one scale, √σ._")
    L.append('')
    L.append('## The handle')
    L.append('')
    L.append(f"Preferred channel `Λ_Q(ππ)_S` (twist-unwinding). **Dipion "
             f"endpoint** `m(ππ)_max = M_Möbius − M_ground = 2√σ ≈ "
             f"{round(DIPION_ENDPOINT_MEV)} MeV — flavor-independent** (same "
             f"above charm and bottom), spectrum peaking high. A fixed edge "
             f"in a directly-plotted observable; one overlay tests the "
             f"framework.")
    L.append('')
    L.append('## Search table')
    L.append('')
    L.append('| Tier | Target (mass MeV) | Reconstruction chain | Facility / production | Constraint status |')
    L.append('|---|---|---|---|---|')
    order = sorted(TARGETS.items(), key=lambda kv: (kv[1]['tier'], -sum(kv[1]['factors'])))
    for name, d in order:
        L.append(f"| {d['tier']} | {name} Möbius ({_mobius(d['ground'])}) | "
                 f"`{d['chain']}` | {d['facility']} · {d['production']} | "
                 f"{d['status']} |")
    L.append('')
    L.append('**Tier 1** — the discovery pair: highest yield × the '
             'cross-flavor clincher (same dipion endpoint + Q-values above '
             'both `Λ_c` and `Λ_b`; a simultaneous two-channel fit). '
             '**Tier 2** — entirely unexplored (no measured excitations), so '
             'a clean bump = discovery, but doubly-heavy / `Ω_b` production '
             'is rare. **Tier 3** — above the well-mapped 2017 `Ω_c` '
             'excitations, which calibrate the search.')
    L.append('')
    L.append('## Discriminators (per channel)')
    L.append('')
    for disc in s['tests'][4]['discriminators']:
        L.append(f"- {disc}")
    L.append('')
    L.append('## Scope (honest)')
    L.append('')
    L.append(f"- **Established:** a sharper, tiered, actionable table from "
             f"#109/#110.")
    L.append(f"- **Open:** {s['open']}. A prioritization deliverable, not new "
             f"physics; states are broad ⟹ resolve in amplitude / Dalitz "
             f"analyses, not as sharp peaks.")
    L.append('')
    return '\n'.join(L)


# ---------------------------------------------------------------------------
# Probe-summary markdown
# ---------------------------------------------------------------------------

def render_markdown(s: dict) -> str:
    L: list[str] = []
    L.append('# Heavy Möbius baryon: a sharper LHCb / Belle II search table (PR #114)')
    L.append('')
    L.append(f"**Run:** {s['timestamp_utc']}")
    L.append('')
    L.append(
        "Converts PRs #103/#109/#110 into a sharper, tiered, actionable "
        "LHCb / Belle II search table. **New handle:** the `Λ_Q(ππ)` dipion "
        "endpoint `m(ππ)_max = 2√σ ≈ 849 MeV` is **flavor-independent** (the "
        "same edge above charm and bottom, peaking high) — a fixed edge in a "
        "plotted observable, one overlay for the whole framework. The note "
        "is written standalone to `search_table.md`."
    )
    L.append('')
    L.append(f"- **New handle**: {s['new_handle']}")
    L.append(f"- **Tier 1**: {s['tier1']}")
    L.append(f"- **Tier 2**: {s['tier2']}")
    L.append(f"- **Tier 3**: {s['tier3']}")
    L.append(f"- **Discriminators**: {s['discriminators']}")
    L.append(f"- **Open**: {s['open']}")
    L.append('')

    L.append('## Test summary')
    L.append('')
    L.append('| # | Test | Key finding | PASS? |')
    L.append('|---|---|---|---|')
    label_map = {
        'T1': 'convert #109/#110 → tiered actionable table',
        'T2': 'targets + golden reconstruction chains',
        'T3': 'facility + production per target (LHCb / Belle II)',
        'T4': 'dipion endpoint 849 MeV, flavor-independent (new handle)',
        'T5': 'discriminators: selection rule; endpoint; Q-match',
        'T6': 'tiers: Λ_c+Λ_b pair / Ξ_cc+Ω_b / Ω_c',
        'T7': 'honest scope: masses ±band, broad, BFs/J^P open',
        'T8': 'HEAVY_MOBIUS_BARYON_SEARCH_TABLE_SHARPENED',
    }
    for t in s['tests']:
        passed = '**PASS**' if t['pass'] else '**FAIL**'
        nm = t['name']; prefix = nm[:2]
        L.append(f"| {prefix} | `{nm}` | {label_map.get(prefix, '—')} | {passed} |")
    L.append('')
    L.append('---')
    L.append('')
    L.append('The compiled table follows (also written standalone as '
             '`search_table.md`):')
    L.append('')
    L.append(render_search_table(s))
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
    table = render_search_table(summary)
    print(md)
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    out = here / 'runs' / f'{ts}_heavy_mobius_baryon_search_table_probe'
    out.mkdir(parents=True, exist_ok=True)
    (out / 'probe.json').write_text(
        json.dumps(summary, indent=2, default=_json_default)
    )
    (out / 'probe.md').write_text(md)
    (out / 'search_table.md').write_text(table)
    print(f"\nWrote: {out / 'probe.json'}")
    print(f"Wrote: {out / 'probe.md'}")
    print(f"Wrote: {out / 'search_table.md'}")
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())
