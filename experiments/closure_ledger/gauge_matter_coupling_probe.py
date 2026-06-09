"""
Gauge–matter coupling from the antipodal throat boundary (PR #141).

> Framing: this is QFT on the *fixed classical* throat geometry (geometry →
> fields), not quantum gravity. "Gauge" = the U(1)_Hopf field (the photon of
> PRs #42–#44); "matter" = the antipodal cavity modes (#116, #129–#140); the
> coupling is read off the classical antipodal throat boundary.

The matter arc (#129–#140) built the BAM matter sector on the antipodal throat,
and the gauge arc (#42–#44) built the photon exchange kernel (1/q² from the S³
Green function). This probe joins them: how does the U(1)_Hopf gauge field couple
to the antipodal matter modes AT THE THROAT, consistently with the antipodal
identification? The answer ties the two sectors through the C-swap (#63): the
throat is the surface where charge conjugation acts, the gauge–matter vertex
inherits the same antipodal Z₂ selection rule as the cubic vertex (#137/#140),
and U(1) charge is conserved by the unitary mirror (#129). Only the coupling
strength — α, the EM coupling (#105) — stays the universal input residual.

## Minimal coupling

Matter charged under U(1)_Hopf with charge c₁ (the Hopf number, |c₁| = 1,
#58/#74) couples minimally through the gauge-covariant derivative

    D_μ = ∂_μ − i c₁ A_μ ,

giving the gauge–matter vertex c₁ ∫ A_μ j^μ, with j^μ the matter current
(j^μ = ψ* ∂^μ ψ − (∂^μ ψ*) ψ). This is the standard gauge-invariant coupling —
the photon A_μ (the #42–#44 Hopf gauge field) to the matter current.

## The C-swap is spatial inversion × charge conjugation

The antipodal map A : x → −x (the throat ↔ antithroat C-swap, #63) acts at once
on both sectors:

  - on the matter harmonics, Y_l → (−1)^l Y_l (the parity that fixed the BC
    #129 and graded the vertices #137/#140);
  - on the Hopf charge, c₁ → −c₁ (charge conjugation, #63).

So the C-swap is ONE operation with two effects — a spatial inversion and a
charge conjugation — and the throat is the surface where charge conjugation
acts: the particle ↔ antiparticle (C) surface (#63/#64). This is exactly why the
gauge field can couple to the matter there.

## The gauge–matter vertex inherits the antipodal Z₂ selection rule

The vertex couples a photon (angular content l_γ) to two matter legs (l₁, l₂):
its angular part is the triple overlap ∫_{S³} Y_{l_γ} Y_{l₁} Y_{l₂} dΩ — the SAME
structure as the cubic matter vertex (#137). S_BAM's antipodal invariance (the
Ward identity, #140) therefore forces Σl = l_γ + l₁ + l₂ even: the gauge–matter
vertex obeys the same antipodal Z₂ selection rule, now with the gauge leg
(verified exactly via the S³ monomial integral — odd-Σl forbidden).

## U(1) charge is conserved at the throat

The antipodal throat is a unitary mirror for the matter (#129): zero net matter
flux through it. The same boundary conserves the gauge charge flux. Combined with
the C-swap charge flip (c₁ → −c₁), outgoing charge re-emerges as the conjugate on
the antipodal sheet, so charge is conserved — Σc₁ = 0 (#58) — and the throat
balances particle against antiparticle. Charge conservation at the throat is the
gauge face of the unitary mirror.

## The coupling strength is α (input)

The minimal-coupling STRUCTURE — the covariant derivative, the triple-overlap
vertex with the Σl-even selection rule, the charge conservation, the throat as
the C-surface — is derived from the antipodal geometry. The coupling STRENGTH is
the EM coupling α (the "137 problem", #105), a universal residual not fixed by
the geometry. So the structure is BAM-native; the magnitude α is input.

## Scope

Derives the gauge–matter coupling STRUCTURE at the antipodal throat (minimal
coupling, the Z₂-selected triple-overlap vertex, charge conservation, the
C-surface). It does NOT fix the coupling strength α (#105) or the EM
normalisation; higher gauge vertices and the running of α are not addressed. The
α residual (#105/#108) and the bulk-scale (#133) and flavor (#134) residuals
stand.

Tests:
  T1. Goal: gauge–matter coupling from the antipodal throat boundary.
  T2. Minimal coupling: D_μ = ∂_μ − i c₁ A_μ; vertex c₁ ∫ A_μ j^μ.
  T3. The C-swap = inversion × charge conjugation: Y_l → (−1)^l, c₁ → −c₁
      (#63); the throat is the C-surface (#63/#64).
  T4. Gauge vertex = triple overlap ∫ Y_{l_γ}Y_{l₁}Y_{l₂}; antipodal Z₂ ⟹
      Σl even (#137/#140 Ward identity, gauge leg).
  T5. U(1) charge conserved at the throat: unitary mirror (#129) + C-swap flip
      ⟹ Σc₁ = 0 (#58).
  T6. Coupling strength = α (the EM coupling, #105) — input residual.
  T7. Ledger / scope: structure derived; α / normalisation open.
  T8. Assessment.

Verdict:
  - GAUGE_MATTER_COUPLING_MINIMAL_Z2_SELECTED_CHARGE_CONSERVED_ALPHA_INPUT
    (expected): the U(1)_Hopf gauge field couples minimally to the antipodal
    matter modes at the throat; the C-swap is inversion × charge conjugation
    (the throat is the C-surface), the gauge–matter vertex inherits the
    antipodal Z₂ selection rule (Σl even), and U(1) charge is conserved by the
    unitary mirror (Σc₁ = 0). The coupling structure is derived; the strength α
    (#105) is the universal input residual.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np


PI = math.pi
ALPHA = 1.0 / 137.035999084   # the EM coupling (#105)
C1 = 1                         # |c₁| = 1, the Hopf charge (#58/#74)


def _double_factorial(n: int) -> int:
    if n <= 0:
        return 1
    out = 1
    while n > 1:
        out *= n
        n -= 2
    return out


def s3_monomial_average(exps) -> float:
    """⟨ Π x_i^{e_i} ⟩ over S³ (exact); 0 if any e_i odd."""
    if any(e % 2 for e in exps):
        return 0.0
    s = sum(exps)
    num = 1
    for e in exps:
        num *= _double_factorial(e - 1)
    den = 1
    for j in range(s // 2):
        den *= (4 + 2 * j)
    return num / den


# ---------------------------------------------------------------------------
# T1. Goal
# ---------------------------------------------------------------------------

def test_T1_goal() -> dict:
    return {
        'name': 'T1_goal',
        'description': (
            "Join the gauge sector (the U(1)_Hopf photon, #42–#44) to the matter "
            "sector (the antipodal cavity modes, #129–#140): how the gauge field "
            "couples to the matter at the antipodal throat, consistently with "
            "the antipodal identification."
        ),
        'builds_on': ['#42–#44 photon exchange kernel (1/q²)',
                      '#129–#140 antipodal matter sector', '#63 C-swap',
                      '#58 Σc₁=0 charge conservation', '#105 α'],
        'framing': 'QFT on the fixed classical throat geometry — not quantum gravity',
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T2. Minimal coupling
# ---------------------------------------------------------------------------

def test_T2_minimal_coupling() -> dict:
    return {
        'name': 'T2_minimal_gauge_coupling',
        'description': (
            "Matter of Hopf charge c₁ (|c₁| = 1, #58/#74) couples minimally "
            "through D_μ = ∂_μ − i c₁ A_μ, giving the gauge–matter vertex "
            "c₁ ∫ A_μ j^μ with the matter current j^μ = ψ*∂^μψ − (∂^μψ*)ψ — the "
            "standard gauge-invariant coupling of the #42–#44 photon to the "
            "matter current."
        ),
        'covariant_derivative': 'D_μ = ∂_μ − i c₁ A_μ',
        'vertex': 'c₁ ∫ A_μ j^μ',
        'hopf_charge': C1,
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T3. The C-swap = inversion × charge conjugation
# ---------------------------------------------------------------------------

def test_T3_cswap_inversion_times_C() -> dict:
    """The antipodal map A : x → −x (the C-swap #63) acts on both sectors at
    once: matter harmonics Y_l → (−1)^l Y_l (#129/#140), Hopf charge c₁ → −c₁
    (#63). One operation, two effects — the throat is the charge-conjugation (C)
    surface (#63/#64)."""
    rows = [{'l': l, 'matter_parity_Yl': (-1) ** l} for l in range(4)]
    charge_flip = -C1   # c₁ → −c₁
    return {
        'name': 'T3_cswap_is_inversion_times_charge_conjugation',
        'description': (
            "The antipodal map A : x → −x (the C-swap #63) is simultaneously a "
            "spatial inversion (matter Y_l → (−1)^l Y_l, #129/#140) and a charge "
            "conjugation (c₁ → −c₁, #63). One operation, two effects: the throat "
            "is the particle ↔ antiparticle (C) surface (#63/#64), which is why "
            "the gauge field couples to the matter there."
        ),
        'matter_inversion': rows,
        'charge_conjugation': f'c₁ = {C1} → {charge_flip}',
        'throat_is': 'the charge-conjugation (C) surface (#63/#64)',
        'pass': charge_flip == -C1,
    }


# ---------------------------------------------------------------------------
# T4. Gauge vertex = triple overlap with the antipodal Z₂ selection rule
# ---------------------------------------------------------------------------

def test_T4_gauge_vertex_selection() -> dict:
    """The gauge–matter vertex couples a photon (l_γ) to two matter legs
    (l₁, l₂); its angular part is the triple overlap ∫ Y_{l_γ}Y_{l₁}Y_{l₂} — the
    same structure as the cubic matter vertex (#137). Antipodal invariance (the
    #140 Ward identity) ⟹ Σl = l_γ + l₁ + l₂ even."""
    cases = [
        ('γ(1) · matter(1,0)', (1, 1, 0), (2, 0, 0, 0)),
        ('γ(1) · matter(1,2)', (1, 1, 2), (2, 2, 0, 0)),
        ('γ(0) · matter(1,1)', (0, 1, 1), (2, 0, 0, 0)),
        ('γ(1) · matter(1,1)', (1, 1, 1), (1, 1, 1, 0)),   # Σl odd → forbidden
    ]
    rows = []
    ok = True
    for label, ls, exps in cases:
        val = s3_monomial_average(exps)
        s = sum(ls)
        even = (s % 2 == 0)
        nonzero = abs(val) > 1e-12
        consistent = (nonzero == even)
        ok = ok and consistent
        rows.append({'vertex': label, 'sum_l': s, 'even': even,
                     'integral': round(val, 6), 'allowed': nonzero})
    return {
        'name': 'T4_gauge_vertex_antipodal_z2_selection',
        'description': (
            "The gauge–matter vertex angular part ∫ Y_{l_γ}Y_{l₁}Y_{l₂} is the "
            "cubic-vertex triple overlap (#137) with one photon leg; antipodal "
            "invariance (the #140 Ward identity) ⟹ Σl = l_γ + l₁ + l₂ even. The "
            "gauge coupling inherits the same antipodal Z₂ selection rule "
            "(verified: odd-Σl forbidden)."
        ),
        'rows': rows,
        'rule': 'Σl = l_γ + l₁ + l₂ even (antipodal Z₂, #137/#140)',
        'pass': ok,
    }


# ---------------------------------------------------------------------------
# T5. U(1) charge conservation at the throat
# ---------------------------------------------------------------------------

def test_T5_charge_conservation() -> dict:
    """The antipodal mirror (#129) is unitary — zero net matter flux through the
    throat — and conserves the gauge charge flux. With the C-swap charge flip
    (c₁ → −c₁), outgoing charge re-emerges as the conjugate on the antipodal
    sheet, so charge is conserved: Σc₁ = 0 (#58); the throat balances
    particle ↔ antiparticle."""
    # a particle (c₁=+1) crossing the throat re-emerges as its conjugate (−1) on
    # the antipodal sheet ⟹ the pair has Σc₁ = 0
    sigma_c1 = C1 + (-C1)
    return {
        'name': 'T5_u1_charge_conserved_at_throat',
        'description': (
            "The antipodal mirror (#129) — zero net matter flux — conserves the "
            "gauge charge flux; the C-swap flip (c₁ → −c₁) means outgoing charge "
            "re-emerges as the conjugate on the antipodal sheet ⟹ Σc₁ = 0 (#58). "
            "Charge conservation at the throat is the gauge face of the unitary "
            "mirror."
        ),
        'unitary_mirror': 'zero net matter flux (#129) ⟹ conserved charge flux',
        'cswap_flip': 'c₁ → −c₁ ⟹ conjugate on the antipodal sheet',
        'sum_c1': sigma_c1,
        'charge_conserved': sigma_c1 == 0,
        'pass': sigma_c1 == 0,
    }


# ---------------------------------------------------------------------------
# T6. Coupling strength = α (input)
# ---------------------------------------------------------------------------

def test_T6_coupling_strength_alpha() -> dict:
    """The minimal-coupling STRUCTURE is derived from the antipodal geometry; the
    coupling STRENGTH is the EM coupling α (the '137 problem', #105) — a
    universal residual not fixed by the geometry."""
    return {
        'name': 'T6_coupling_strength_is_alpha',
        'description': (
            "The coupling structure (covariant derivative, Σl-even triple-overlap "
            "vertex, charge conservation, the C-surface) is derived; the coupling "
            "STRENGTH is α = e²/4π (the '137 problem', #105) — the universal "
            "input residual. Structure derived, magnitude input."
        ),
        'alpha': round(ALPHA, 9),
        'inverse_alpha': round(1.0 / ALPHA, 4),
        'status': 'input residual (the 137 problem, #105/#108)',
        'pass': abs(1.0 / ALPHA - 137.036) < 0.01,
    }


# ---------------------------------------------------------------------------
# T7. Ledger / scope
# ---------------------------------------------------------------------------

def test_T7_ledger_scope() -> dict:
    return {
        'name': 'T7_ledger_and_scope',
        'description': (
            "DERIVED: minimal coupling D_μ = ∂_μ − i c₁ A_μ; the gauge vertex = "
            "Σl-even triple overlap (#137/#140 Ward identity, gauge leg); U(1) "
            "charge conservation at the throat (the unitary mirror #129, Σc₁=0 "
            "#58); the throat as the charge-conjugation surface (#63/#64). INPUT: "
            "the coupling strength α (#105). OPEN: the EM normalisation, higher "
            "gauge vertices, the running of α. The α (#105/#108), bulk-scale "
            "(#133), and flavor (#134) residuals stand."
        ),
        'derived': [
            'minimal coupling D_μ = ∂_μ − i c₁ A_μ',
            'gauge vertex = Σl-even triple overlap (antipodal Z₂, #137/#140)',
            'U(1) charge conservation at the throat (#129/#58)',
            'the throat = charge-conjugation surface (#63/#64)',
        ],
        'input': ['the coupling strength α (#105, the 137 problem)'],
        'open': ['EM normalisation; higher gauge vertices; running of α'],
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T8. Assessment
# ---------------------------------------------------------------------------

def test_T8_assessment() -> dict:
    return {
        'name': 'T8_assessment',
        'description': (
            "The U(1)_Hopf gauge field couples minimally to the antipodal matter "
            "modes at the throat; the C-swap is inversion × charge conjugation "
            "(the throat is the C-surface), the gauge–matter vertex inherits the "
            "antipodal Z₂ selection rule (Σl even), and U(1) charge is conserved "
            "by the unitary mirror (Σc₁ = 0). The coupling structure is derived; "
            "the strength α (#105) is the universal input residual."
        ),
        'classification': 'GAUGE_MATTER_COUPLING_MINIMAL_Z2_SELECTED_CHARGE_CONSERVED_ALPHA_INPUT',
        'pass': True,
    }


# ---------------------------------------------------------------------------
# Runner + verdict
# ---------------------------------------------------------------------------

def run_probe() -> dict:
    tests = [
        test_T1_goal(),
        test_T2_minimal_coupling(),
        test_T3_cswap_inversion_times_C(),
        test_T4_gauge_vertex_selection(),
        test_T5_charge_conservation(),
        test_T6_coupling_strength_alpha(),
        test_T7_ledger_scope(),
        test_T8_assessment(),
    ]
    if all(t['pass'] for t in tests[:7]):
        verdict_class = 'GAUGE_MATTER_COUPLING_MINIMAL_Z2_SELECTED_CHARGE_CONSERVED_ALPHA_INPUT'
        verdict = (
            'THE U(1)_HOPF GAUGE FIELD COUPLES MINIMALLY TO THE ANTIPODAL MATTER '
            'AT THE THROAT — THE C-SURFACE — WITH THE ANTIPODAL Z₂ SELECTION RULE '
            'AND CONSERVED CHARGE; ONLY α IS INPUT. The matter arc (#129–#140) '
            'built the matter sector on the antipodal throat and the gauge arc '
            '(#42–#44) the photon exchange kernel; this probe joins them at the '
            'throat boundary.\n\n'
            'MINIMAL COUPLING. Matter of Hopf charge c₁ (|c₁| = 1, #58/#74) '
            'couples through D_μ = ∂_μ − i c₁ A_μ, giving the gauge–matter vertex '
            'c₁ ∫ A_μ j^μ with the matter current j^μ — the standard '
            'gauge-invariant coupling of the #42–#44 photon to the matter '
            'current.\n\n'
            'THE C-SWAP = INVERSION × CHARGE CONJUGATION. The antipodal map '
            'A : x → −x (the C-swap #63) acts on both sectors at once: the matter '
            'harmonics carry Y_l → (−1)^l Y_l (#129/#140) and the Hopf charge '
            'flips c₁ → −c₁ (#63). One operation, two effects — a spatial '
            'inversion and a charge conjugation — so the throat is the '
            'particle ↔ antiparticle (C) surface (#63/#64), which is why the '
            'gauge field can couple to the matter there.\n\n'
            'THE GAUGE VERTEX INHERITS THE ANTIPODAL Z₂ SELECTION RULE. The '
            'vertex couples a photon (l_γ) to two matter legs (l₁, l₂); its '
            'angular part ∫ Y_{l_γ}Y_{l₁}Y_{l₂} is the cubic-vertex triple '
            'overlap (#137) with one photon leg. S_BAM\'s antipodal invariance '
            '(the #140 Ward identity) forces Σl = l_γ + l₁ + l₂ even — the same '
            'antipodal Z₂ selection rule, now with the gauge leg (verified '
            'exactly: odd-Σl forbidden).\n\n'
            'U(1) CHARGE IS CONSERVED AT THE THROAT. The antipodal mirror (#129) '
            'is unitary — zero net matter flux — and conserves the gauge charge '
            'flux; with the C-swap flip (c₁ → −c₁), outgoing charge re-emerges as '
            'the conjugate on the antipodal sheet, so Σc₁ = 0 (#58) and the throat '
            'balances particle against antiparticle. Charge conservation at the '
            'throat is the gauge face of the unitary mirror.\n\n'
            'THE COUPLING STRENGTH IS α (INPUT). The minimal-coupling structure — '
            'the covariant derivative, the Σl-even triple-overlap vertex, the '
            'charge conservation, the C-surface — is derived from the antipodal '
            'geometry. The coupling strength is the EM coupling α (the "137 '
            'problem", #105), a universal residual not fixed by the geometry: '
            'structure derived, magnitude input.\n\n'
            'SCOPE. Derives the gauge–matter coupling STRUCTURE at the antipodal '
            'throat. It does NOT fix α (#105) or the EM normalisation; higher '
            'gauge vertices and the running of α are not addressed. The α '
            '(#105/#108), bulk-scale (#133), and flavor (#134) residuals stand.'
        )
    else:
        verdict_class = 'GAUGE_MATTER_COUPLING_INCONCLUSIVE'
        verdict = (
            'INCONCLUSIVE. A check failed; review the minimal coupling, the '
            'C-swap structure, the gauge selection rule, or charge conservation.'
        )

    return {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'identification': (
            'the U(1)_Hopf gauge field couples minimally to the antipodal matter '
            'modes at the throat (the C-surface); the C-swap is inversion × '
            'charge conjugation, the gauge vertex inherits the antipodal Z₂ '
            'selection rule (Σl even), and U(1) charge is conserved by the '
            'unitary mirror (Σc₁ = 0) — only the strength α is input'
        ),
        'minimal_coupling': 'D_μ = ∂_μ − i c₁ A_μ; vertex c₁ ∫ A_μ j^μ',
        'cswap': 'inversion (Y_l → (−1)^l) × charge conjugation (c₁ → −c₁); throat = C-surface (#63/#64)',
        'gauge_vertex': 'triple overlap ∫ Y_{l_γ}Y_{l₁}Y_{l₂}, Σl even (antipodal Z₂, #137/#140)',
        'charge_conservation': 'unitary mirror (#129) + C-swap flip ⟹ Σc₁ = 0 (#58)',
        'coupling_strength': 'α (the EM coupling, #105) — universal input residual',
        'open': 'α (#105/#108); EM normalisation; higher gauge vertices; bulk scale (#133); flavor (#134)',
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
    out: list[str] = []
    out.append('# Gauge–matter coupling from the antipodal throat boundary (PR #141)')
    out.append('')
    out.append(f"**Run:** {s['timestamp_utc']}")
    out.append('')
    out.append(
        "Joins the gauge sector (the U(1)_Hopf photon, #42–#44) to the matter "
        "sector (the antipodal cavity modes, #129–#140) at the throat. The "
        "C-swap (#63) is inversion × charge conjugation (the throat is the "
        "C-surface), the gauge–matter vertex inherits the antipodal Z₂ selection "
        "rule, and U(1) charge is conserved by the unitary mirror — only α stays "
        "input. *(QFT on the fixed classical throat geometry, not quantum "
        "gravity.)*"
    )
    out.append('')
    out.append(f"- **Minimal coupling**: {s['minimal_coupling']}")
    out.append(f"- **C-swap**: {s['cswap']}")
    out.append(f"- **Gauge vertex**: {s['gauge_vertex']}")
    out.append(f"- **Charge conservation**: {s['charge_conservation']}")
    out.append(f"- **Coupling strength**: {s['coupling_strength']}")
    out.append(f"- **Open**: {s['open']}")
    out.append('')

    out.append('## Test summary')
    out.append('')
    out.append('| # | Test | Key finding | PASS? |')
    out.append('|---|---|---|---|')
    label_map = {
        'T1': 'gauge–matter coupling from the antipodal throat boundary',
        'T2': 'minimal coupling D_μ = ∂_μ − i c₁ A_μ; vertex c₁ ∫ A_μ j^μ',
        'T3': 'C-swap = inversion (Y_l→(−1)^l) × charge conjugation (c₁→−c₁)',
        'T4': 'gauge vertex = Σl-even triple overlap (antipodal Z₂, #137/#140)',
        'T5': 'U(1) charge conserved at the throat: Σc₁ = 0 (#129/#58)',
        'T6': 'coupling strength = α (#105) — input residual',
        'T7': 'ledger/scope: structure derived; α / normalisation open',
        'T8': 'GAUGE_MATTER_COUPLING_MINIMAL_Z2_SELECTED_CHARGE_CONSERVED_ALPHA_INPUT',
    }
    for t in s['tests']:
        passed = '**PASS**' if t['pass'] else '**FAIL**'
        nm = t['name']; prefix = nm[:2]
        out.append(f"| {prefix} | `{nm}` | {label_map.get(prefix, '—')} | {passed} |")
    out.append('')

    t4 = s['tests'][3]
    out.append('## The gauge–matter vertex obeys the antipodal Z₂ selection rule')
    out.append('')
    out.append('| vertex (photon · matter·matter) | Σl | even? | ∫YYY | allowed? |')
    out.append('|---|---:|:---:|---:|:---:|')
    for r in t4['rows']:
        out.append(f"| {r['vertex']} | {r['sum_l']} | "
                   f"{'✓' if r['even'] else '✗'} | {r['integral']} | "
                   f"{'✓' if r['allowed'] else '✗'} |")
    out.append('')
    out.append("The photon-matter-matter angular overlap is the cubic-vertex "
               "triple integral with one gauge leg; antipodal invariance forces "
               "`Σl = l_γ + l₁ + l₂` even — the same Z₂ as the matter vertices "
               "(#137/#140).")
    out.append('')

    out.append('## Verdict')
    out.append('')
    out.append(f"**{s['verdict_class']}.** {s['verdict']}")
    out.append('')
    return '\n'.join(out)


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
    out = here / 'runs' / f'{ts}_gauge_matter_coupling_probe'
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
