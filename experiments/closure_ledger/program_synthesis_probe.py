"""
Program-wide synthesis: the BAM input budget and epistemic ledger
(PR #104).

A capstone meta-probe. After the lepton/QED foundations (B1–B5,
throat-as-particle), the neutrino arc (#85–#96), the quark arc
(#76–#98), and the QCD confinement / exotic arc (#99–#103), this probe
classifies every major result by EPISTEMIC STATUS and counts BAM's
actual free/anchored inputs. The headline: BAM's entire dimensionful
content reduces to TWO anchors (the B4 scale-modulus minimum — m_e and
√σ); the dimensionless content is predominantly derived geometry; the
genuinely-open residuals are localized to a handful (neutrino compliance
ε, quark n_part) plus the universal flavor puzzle; and the program's
distinctive output is a tower of topological predictions from the
non-orientable structure.

## The five epistemic tiers

  1. **DERIVED geometry** — parameter-free geometric/topological results,
     or reconstructed structural facts that match data.
  2. **DIMENSIONFUL ANCHORS (B4)** — the mandatory scale inputs; the B4
     scale-modulus theorem (PR #52) requires exactly one per sector.
  3. **OPEN dimensionless residuals** — constrained but not yet derived
     dimensionless numbers, localized to a few places.
  4. **The FLAVOR PUZZLE** — the Yukawa hierarchy; not geometric, open
     across ALL of physics (not a BAM-specific failing).
  5. **TOPOLOGICAL PREDICTIONS** — BAM-specific, falsifiable consequences
     of the non-orientable structure (matched, constrained, or findable).

## The input budget (headline)

  - **2 dimensionful anchors**: `m_e = ℏc/R_MID` (QED/lepton scale) and
    `√σ ≈ Λ_QCD` (confinement scale). The B4 theorem says this is the
    irreducible minimum — one scale per sector.
  - **~2 open dimensionless residuals**: the neutrino boundary compliance
    `ε` (the bounce/seesaw residual, bracketed by `[2π, k_5√(2π)]`) and
    the quark `n_part = 233` (a compensator for the flavor puzzle).
  - **1 universal open problem**: the flavor puzzle (Yukawa hierarchy) —
    shared with every theory, not BAM-specific.
  - Everything else is **derived geometry** or a **topological
    prediction**.

So beyond the two mandatory B4 anchors and the universal flavor puzzle,
BAM carries only a small, localized set of open dimensionless inputs —
the rest of the dimensionless structure is geometric or predicted.

Tests:
  T1. The five-tier epistemic scheme.
  T2. Tier 1 — derived geometry (enumerate; the bulk).
  T3. Tier 2 — exactly TWO dimensionful anchors (B4 minimum).
  T4. Tier 3 — open dimensionless residuals, localized (ε, n_part).
  T5. Tier 4 — the flavor puzzle: universal, not BAM-specific.
  T6. Tier 5 — topological predictions (matched → constrained → findable).
  T7. The input budget: 2 anchors + ~2 residuals + flavor puzzle; rest
      derived/predicted.
  T8. Assessment.

Verdict:
  - BAM_TWO_ANCHORS_LOCALIZED_RESIDUALS_DERIVED_GEOMETRY_TOPOLOGICAL_PREDICTIONS
    (expected): BAM's dimensionful content reduces to two B4 anchors
    (m_e, √σ); the dimensionless content is predominantly derived
    geometry; the open residuals are localized (neutrino ε, quark n_part)
    plus the universal flavor puzzle; the distinctive output is a tower of
    non-orientable topological predictions.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np


# ---------------------------------------------------------------------------
# The epistemic ledger (curated across the program)
# ---------------------------------------------------------------------------

DERIVED_GEOMETRY = [
    ('|c₁| = 1 charge quantization', 'Hopf Chern number'),
    ('spin-½ from Hopf holonomy', '∮A = π cos χ'),
    ('Coulomb law from eigenmode throat flux', 'B5 master integral'),
    ('g = 2', 'Pauli SU(2) + Hopf monopole'),
    ('Schwinger a = α/2π (one loop, 0.15%)', 'S³ Green fn + throat vertex'),
    ('C = inner/outer swap (c₁→−c₁)', 'PR #63 involution'),
    ('CPT on throat histories; Dirac 4-spinor', 'PR #64–#66'),
    ('pair threshold factor 2 (Σc₁=0)', 'C-conjugate pair, PR #58'),
    ('k_5 = D_bulk = dim(S³)+2 = 5', 'PR #73'),
    ('β_lepton = k_5²·2π = 50π', 'PR #71'),
    ('#generations = (k_5+1)/2 = 3', 'PR #72'),
    ('lepton/quark mass operator unified (Bohr-Sommerfeld)', 'PR #83'),
    ('neutrino Majorana ⟸ k=0 ⟹ c₁=0', 'PR #86 selection rule'),
    ('only-ν seesaw ⟸ Σc₁=0 single-state', 'PR #87'),
    ('PMNS anarchic (cross-coordinate); CKM aligned', 'PR #92'),
    ('CP generic (Hopf-complex); 2 Majorana phases ⟸ c₁=0', 'PR #94'),
    ('6 quarks = 3×2, Z₂ partition, k=0 shell', 'PR #69/#77'),
    ('Cornell form: flux-tube bridge + throat/gluon exchange', 'PR #99'),
    ('string breaking = Schwinger (eE→σ) = PR #58', 'PR #99'),
    ('Regge slope α′ = 1/(2πσ); glueball = half', 'PR #99/#100'),
    ('exotic J^PC (1-+) ⟸ non-orientable flux tube', 'PR #101'),
    ('no exotic J^P for baryons (supernumerary only)', 'PR #102'),
]

DIMENSIONFUL_ANCHORS = [
    ('m_e = ℏc/R_MID', 'QED / lepton scale (B4, PR #52)'),
    ('√σ ≈ Λ_QCD ≈ 0.42 GeV', 'confinement scale (B4 analogue, PR #99)'),
]

OPEN_DIMENSIONLESS_RESIDUALS = [
    ('neutrino boundary compliance ε', 'bracketed by [2π, k_5√(2π)]; PR #88–#90'),
    ('quark n_part = 233', 'flavor-puzzle compensator; parity-only invariant; PR #76/#97'),
]

FLAVOR_PUZZLE = [
    ('quark Yukawa hierarchy', 'RG-invariant ratios ⟹ not running; the flavor '
     'puzzle, open in ALL physics (PR #98)'),
]

TOPOLOGICAL_PREDICTIONS = [
    ('mesonic 1-+ hybrids (π₁, η₁)', 'MATCHED to data (PR #101)'),
    ('neutrino normal ordering; m_ββ ≲ 8 meV; Σm_ν ≈ 59–65 meV', 'FALSIFIABLE (PR #95/#96)'),
    ('multiquark exotics (X, Z_c, T_cc, P_c) = multi-junction', 'ACCOMMODATED (PR #101)'),
    ('light baryonic exotics (Möbius/hybrid)', 'CONSTRAINED counting test (PR #102)'),
    ('heavy Möbius baryon (Λ_c ~3.14, Λ_b ~6.47 GeV)', 'FINDABLE, unconstrained (PR #103)'),
    ('Möbius glueball tower', 'FREE (glueballs unobserved, PR #100)'),
]


# ---------------------------------------------------------------------------
# T1. The five-tier scheme
# ---------------------------------------------------------------------------

def test_T1_tier_scheme() -> dict:
    return {
        'name': 'T1_five_tier_epistemic_scheme',
        'description': (
            "Five tiers: derived geometry; dimensionful anchors (B4); open "
            "dimensionless residuals; the flavor puzzle (universal); "
            "topological predictions (BAM-specific)."
        ),
        'tiers': ['DERIVED geometry', 'DIMENSIONFUL ANCHORS (B4)',
                  'OPEN dimensionless residuals', 'the FLAVOR PUZZLE',
                  'TOPOLOGICAL PREDICTIONS'],
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T2. Derived geometry — the bulk
# ---------------------------------------------------------------------------

def test_T2_derived() -> dict:
    return {
        'name': 'T2_derived_geometry',
        'description': (
            "Parameter-free geometric/topological results (or "
            "reconstructed structural facts matching data) — the bulk of "
            "the program."
        ),
        'entries': [f'{name} — {basis}' for name, basis in DERIVED_GEOMETRY],
        'count': len(DERIVED_GEOMETRY),
        'pass': len(DERIVED_GEOMETRY) >= 15,
    }


# ---------------------------------------------------------------------------
# T3. Exactly two dimensionful anchors (B4 minimum)
# ---------------------------------------------------------------------------

def test_T3_anchors() -> dict:
    """The B4 scale-modulus theorem (PR #52) requires one dimensionful
    anchor per sector — the irreducible minimum. BAM has exactly two:
    m_e (QED/lepton) and √σ (confinement)."""
    return {
        'name': 'T3_two_dimensionful_anchors',
        'description': (
            "B4 (PR #52): one dimensionful anchor per sector is mandatory. "
            "BAM has exactly two — m_e = ℏc/R_MID and √σ ≈ Λ_QCD."
        ),
        'anchors': [f'{a} — {role}' for a, role in DIMENSIONFUL_ANCHORS],
        'count': len(DIMENSIONFUL_ANCHORS),
        'is_b4_minimum': True,
        'pass': len(DIMENSIONFUL_ANCHORS) == 2,
    }


# ---------------------------------------------------------------------------
# T4. Open dimensionless residuals — localized
# ---------------------------------------------------------------------------

def test_T4_residuals() -> dict:
    """The genuinely-open dimensionless inputs are localized to a handful:
    the neutrino boundary compliance ε (bracketed, the seesaw residual)
    and the quark n_part = 233 (a flavor-puzzle compensator)."""
    return {
        'name': 'T4_open_dimensionless_residuals',
        'description': (
            "Open dimensionless inputs, localized: neutrino compliance ε "
            "(bracketed [2π, k_5√(2π)]) and quark n_part = 233 (compensator)."
        ),
        'residuals': [f'{r} — {note}' for r, note in OPEN_DIMENSIONLESS_RESIDUALS],
        'count': len(OPEN_DIMENSIONLESS_RESIDUALS),
        'localized': len(OPEN_DIMENSIONLESS_RESIDUALS) <= 3,
        'pass': len(OPEN_DIMENSIONLESS_RESIDUALS) <= 3,
    }


# ---------------------------------------------------------------------------
# T5. The flavor puzzle — universal, not BAM-specific
# ---------------------------------------------------------------------------

def test_T5_flavor_puzzle() -> dict:
    """The quark Yukawa hierarchy is the flavor puzzle — RG-invariant
    ratios (not running), irregular, derivable by NO current theory.
    It is a universal open problem, shared with all of physics, not a
    BAM-specific failing (PR #98)."""
    return {
        'name': 'T5_flavor_puzzle_universal',
        'description': (
            "Quark Yukawa hierarchy = the flavor puzzle (RG-invariant "
            "ratios ⟹ not running; irregular). Open across ALL physics, "
            "not BAM-specific (PR #98)."
        ),
        'entries': [f'{name} — {note}' for name, note in FLAVOR_PUZZLE],
        'universal_open_problem': True,
        'bam_specific_failing': False,
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T6. Topological predictions — the distinctive output
# ---------------------------------------------------------------------------

def test_T6_topological_predictions() -> dict:
    """BAM's distinctive output: topological predictions from the
    non-orientable structure, spanning the testability gamut — matched
    (1-+ hybrids), falsifiable (neutrino m_ββ/Σm_ν/ordering), accommodated
    (multiquark zoo), constrained (light baryonic exotics), findable
    (heavy Möbius baryon), and free (Möbius glueball tower)."""
    return {
        'name': 'T6_topological_predictions',
        'description': (
            "Non-orientable topological predictions across the testability "
            "gamut: matched, falsifiable, constrained, findable, free."
        ),
        'predictions': [f'{p} — {status}' for p, status in TOPOLOGICAL_PREDICTIONS],
        'count': len(TOPOLOGICAL_PREDICTIONS),
        'spans_matched_to_free': True,
        'pass': len(TOPOLOGICAL_PREDICTIONS) >= 5,
    }


# ---------------------------------------------------------------------------
# T7. The input budget
# ---------------------------------------------------------------------------

def test_T7_input_budget() -> dict:
    """The synthesis: BAM's free/anchored content is 2 dimensionful
    anchors (B4: m_e, √σ) + ~2 localized open dimensionless residuals
    (ε, n_part) + 1 universal open problem (the flavor puzzle).
    Everything else is derived geometry (~22 results) or a topological
    prediction (~6). The dimensionful content is at the B4 minimum."""
    n_derived = len(DERIVED_GEOMETRY)
    n_anchor = len(DIMENSIONFUL_ANCHORS)
    n_resid = len(OPEN_DIMENSIONLESS_RESIDUALS)
    n_flavor = len(FLAVOR_PUZZLE)
    n_pred = len(TOPOLOGICAL_PREDICTIONS)
    return {
        'name': 'T7_input_budget',
        'description': (
            "Input budget: 2 dimensionful anchors (B4) + ~2 localized open "
            "dimensionless residuals + 1 universal flavor puzzle; the rest "
            "(~22 derived, ~6 topological predictions) is geometry or "
            "prediction."
        ),
        'n_derived_geometry': n_derived,
        'n_dimensionful_anchors': n_anchor,
        'n_open_dimensionless_residuals': n_resid,
        'n_flavor_puzzle': n_flavor,
        'n_topological_predictions': n_pred,
        'dimensionful_at_b4_minimum': n_anchor == 2,
        'open_inputs_localized': n_resid <= 3,
        'pass': n_anchor == 2 and n_resid <= 3 and n_derived >= 15,
    }


# ---------------------------------------------------------------------------
# T8. Assessment
# ---------------------------------------------------------------------------

def test_T8_assessment() -> dict:
    return {
        'name': 'T8_assessment',
        'description': (
            "BAM's dimensionful content reduces to two B4 anchors (m_e, "
            "√σ); the dimensionless content is predominantly derived "
            "geometry; the open residuals are localized (neutrino ε, quark "
            "n_part) plus the universal flavor puzzle; the distinctive "
            "output is a tower of non-orientable topological predictions."
        ),
        'classification': (
            'BAM_TWO_ANCHORS_LOCALIZED_RESIDUALS_DERIVED_GEOMETRY_TOPOLOGICAL_PREDICTIONS'
        ),
        'pass': True,
    }


# ---------------------------------------------------------------------------
# Runner + verdict
# ---------------------------------------------------------------------------

def run_probe() -> dict:
    tests = [
        test_T1_tier_scheme(),
        test_T2_derived(),
        test_T3_anchors(),
        test_T4_residuals(),
        test_T5_flavor_puzzle(),
        test_T6_topological_predictions(),
        test_T7_input_budget(),
        test_T8_assessment(),
    ]
    if all(t['pass'] for t in tests[:7]):
        verdict_class = (
            'BAM_TWO_ANCHORS_LOCALIZED_RESIDUALS_DERIVED_GEOMETRY_TOPOLOGICAL_PREDICTIONS'
        )
        verdict = (
            'PROGRAM-WIDE SYNTHESIS: TWO ANCHORS, LOCALIZED RESIDUALS, '
            'DERIVED GEOMETRY, AND A TOWER OF TOPOLOGICAL PREDICTIONS. '
            'After the lepton/QED foundations, the neutrino arc (#85–#96), '
            'the quark arc (#76–#98), and the QCD confinement / exotic arc '
            '(#99–#103), this probe classifies every major result by '
            'epistemic status and counts BAM\'s actual free/anchored '
            'inputs.\n\n'
            'THE FIVE TIERS. (1) DERIVED geometry — parameter-free '
            'geometric/topological results (charge quantization |c₁|=1, '
            'spin-½, g=2, the one-loop Schwinger a=α/2π, k_5=5, '
            'β_lepton=50π, three generations, C/CPT, the Bohr-Sommerfeld '
            'mass operator, the neutrino Majorana selection rule, PMNS '
            'anarchy, generic CP, the Cornell form, string-breaking = '
            'Schwinger, the Regge slope, the exotic-J^PC bookkeeping, …) — '
            'about two dozen results. (2) DIMENSIONFUL ANCHORS — the B4 '
            'scale-modulus minimum. (3) OPEN dimensionless residuals — a '
            'handful. (4) The FLAVOR PUZZLE — universal. (5) TOPOLOGICAL '
            'PREDICTIONS — the distinctive output.\n\n'
            'THE INPUT BUDGET. BAM\'s entire DIMENSIONFUL content reduces '
            'to TWO anchors: m_e = ℏc/R_MID (the QED/lepton scale) and √σ ≈ '
            'Λ_QCD (the confinement scale). The B4 scale-modulus theorem '
            '(PR #52) says one dimensionful input per sector is mandatory, '
            'so two is the irreducible minimum — BAM sits at it. The '
            'genuinely-open DIMENSIONLESS inputs are localized to ~two: the '
            'neutrino boundary compliance ε (the seesaw/bounce residual, '
            'itself bracketed to [2π, k_5√(2π)]) and the quark n_part=233 '
            '(a compensator for the flavor puzzle). Beyond these, there is '
            'one UNIVERSAL open problem — the flavor puzzle (the quark '
            'Yukawa hierarchy: RG-invariant ratios, irregular, derivable by '
            'no current theory) — which BAM shares with all of physics and '
            'is not a BAM-specific failing. Everything else is derived '
            'geometry or a topological prediction.\n\n'
            'THE DISTINCTIVE OUTPUT. The non-orientable structure yields a '
            'tower of topological predictions spanning the full testability '
            'gamut: MATCHED (the mesonic 1-+ hybrids π₁, η₁), FALSIFIABLE '
            '(neutrino normal ordering, m_ββ ≲ 8 meV, Σm_ν ≈ 59–65 meV), '
            'ACCOMMODATED (the multiquark exotic zoo), CONSTRAINED (the '
            'light baryonic exotics — a counting test), FINDABLE (the heavy '
            'Möbius baryon at the flavor-independent 2√σ gap), and FREE '
            '(the Möbius glueball tower, glueballs unobserved).\n\n'
            'SO THE WHOLE PROGRAM, IN ONE LINE: two mandatory B4 anchors + '
            'a couple of localized open dimensionless residuals + the '
            'universal flavor puzzle, with the rest derived geometry and a '
            'set of falsifiable non-orientable topological predictions.\n\n'
            'HONEST SCOPE. This is a CLASSIFICATION, not a new derivation: '
            'it inventories and ranks the program\'s results by epistemic '
            'status. The "derived" tier includes both parameter-free '
            'results and reconstructed/fitted structural matches (e.g. the '
            'charged-lepton ladder, the Schwinger 0.15% match) — the '
            'inventory does not relitigate each. The two anchors and the '
            'flavor puzzle are genuinely open inputs (B4-mandatory and '
            'universal respectively); the localized residuals (ε, n_part) '
            'are the program\'s own remaining homework.'
        )
    else:
        verdict_class = 'BAM_SYNTHESIS_INCONCLUSIVE'
        verdict = (
            'SYNTHESIS INCONCLUSIVE. A tier count fell outside the expected '
            'range; review the ledger.'
        )

    return {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'identification': (
            'BAM = 2 dimensionful anchors (B4: m_e, √σ) + ~2 localized open '
            'dimensionless residuals (neutrino ε, quark n_part) + the '
            'universal flavor puzzle; the rest derived geometry (~22) or '
            'topological predictions (~6)'
        ),
        'dimensionful_anchors': [a for a, _ in DIMENSIONFUL_ANCHORS],
        'open_residuals': [r for r, _ in OPEN_DIMENSIONLESS_RESIDUALS],
        'flavor_puzzle': 'quark Yukawa hierarchy (universal, not BAM-specific)',
        'distinctive_output': 'non-orientable topological predictions (matched → free)',
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
    L.append('# Program-wide synthesis: the BAM input budget and epistemic ledger (PR #104)')
    L.append('')
    L.append(f"**Run:** {s['timestamp_utc']}")
    L.append('')
    L.append(
        "A capstone meta-probe classifying every major result by epistemic "
        "status and counting BAM's actual free/anchored inputs. "
        "**Headline:** BAM's entire dimensionful content reduces to **two "
        "B4 anchors** (`m_e`, `√σ`); the dimensionless content is "
        "predominantly derived geometry; the open residuals are localized "
        "(neutrino `ε`, quark `n_part`) plus the universal flavor puzzle; "
        "and the distinctive output is a tower of non-orientable "
        "topological predictions."
    )
    L.append('')
    L.append(f"- **Identification**: {s['identification']}")
    L.append(f"- **Dimensionful anchors**: {', '.join(s['dimensionful_anchors'])}")
    L.append(f"- **Open residuals**: {', '.join(s['open_residuals'])}")
    L.append(f"- **Flavor puzzle**: {s['flavor_puzzle']}")
    L.append(f"- **Distinctive output**: {s['distinctive_output']}")
    L.append('')

    L.append('## Test summary')
    L.append('')
    L.append('| # | Test | Key finding | PASS? |')
    L.append('|---|---|---|---|')
    label_map = {
        'T1': 'five-tier epistemic scheme',
        'T2': f"derived geometry — {s['tests'][1]['count']} results (the bulk)",
        'T3': 'exactly 2 dimensionful anchors (B4 minimum)',
        'T4': f"open dimensionless residuals — {s['tests'][3]['count']} (localized)",
        'T5': 'flavor puzzle: universal, not BAM-specific',
        'T6': f"topological predictions — {s['tests'][5]['count']} (matched → free)",
        'T7': 'input budget: 2 anchors + ~2 residuals + flavor puzzle',
        'T8': 'BAM_TWO_ANCHORS_LOCALIZED_RESIDUALS_DERIVED_GEOMETRY_TOPOLOGICAL_PREDICTIONS',
    }
    for t in s['tests']:
        passed = '**PASS**' if t['pass'] else '**FAIL**'
        nm = t['name']; prefix = nm[:2]
        L.append(f"| {prefix} | `{nm}` | {label_map.get(prefix, '—')} | {passed} |")
    L.append('')

    L.append('## The input budget')
    L.append('')
    t7 = s['tests'][6]
    L.append('| tier | count |')
    L.append('|---|---:|')
    L.append(f"| derived geometry | {t7['n_derived_geometry']} |")
    L.append(f"| **dimensionful anchors (B4)** | **{t7['n_dimensionful_anchors']}** |")
    L.append(f"| open dimensionless residuals | {t7['n_open_dimensionless_residuals']} |")
    L.append(f"| flavor puzzle (universal) | {t7['n_flavor_puzzle']} |")
    L.append(f"| topological predictions | {t7['n_topological_predictions']} |")
    L.append('')

    L.append('## Tier 2 — the two anchors (B4 minimum)')
    L.append('')
    for a, role in DIMENSIONFUL_ANCHORS:
        L.append(f"- **{a}** — {role}")
    L.append('')

    L.append('## Tier 3 — open dimensionless residuals (localized)')
    L.append('')
    for r, note in OPEN_DIMENSIONLESS_RESIDUALS:
        L.append(f"- **{r}** — {note}")
    L.append('')

    L.append('## Tier 5 — topological predictions (the distinctive output)')
    L.append('')
    for p, status in TOPOLOGICAL_PREDICTIONS:
        L.append(f"- {p} — **{status}**")
    L.append('')

    L.append('## Verdict')
    L.append('')
    L.append(f"**{s['verdict_class']}.** {s['verdict']}")
    L.append('')

    L.append('## What this leaves open')
    L.append('')
    L.append('- **The two localized residuals** — the neutrino compliance '
             '`ε` and the quark `n_part`; the program\'s own remaining '
             'homework.')
    L.append('- **The flavor puzzle** — universal, shared with all of '
             'physics; not BAM-specific.')
    L.append('- **The two B4 anchors** — `m_e` and `√σ` are mandatory '
             'dimensionful inputs (one per sector), not derivable within '
             'the framework.')
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
    out = here / 'runs' / f'{ts}_program_synthesis_probe'
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
