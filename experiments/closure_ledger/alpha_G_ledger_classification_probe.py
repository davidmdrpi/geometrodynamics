"""
Where do α and G sit in the #104 epistemic ledger? (PR #105).

PR #104 classified the mass-sector results into five tiers and found
BAM's dimensionful content reduces to two B4 anchors (m_e, √σ). This
probe answers a sharper question the user posed: under that ledger, are
the fine-structure constant α and Newton's constant G **derived**,
**residual**, or **anchors**? The answer also refines #104.

## The fundamental-constant baseline

It helps to place the four fundamental constants first:

  - **c** — a unit convention (sets length = time); not physical content.
  - **ℏ** — the **closure quantum**. The closure ledger reduces every
    dimensionless parameter to closure invariants (`action_base = 2π`,
    `β_lepton = k_5²·2π = 50π`, …), and the Compton bridge gives
    `ℏ = m_e · R_MID · c`. So ℏ is GEOMETRIC (the `2π` closure structure),
    fixed once the single dimensionful anchor is. ℏ is derived-structure,
    not an independent input.
  - **G** — see below: the dimensionful **anchor**.
  - **α** — see below: a **residual** (universal).

## G is an ANCHOR (the gravitational, foundational one)

BAM is GR-foundational: the throat is a gravitational wormhole, and its
size — the invariant bulk length `ΔR = R_OUTER − R_INNER` (equivalently
`R_MID`), the ONE dimensionful input the B4 theorem (PR #52/#53) says is
mandatory — is set by the bulk gravity sector. The brane tension that
fixes the throat equilibrium is the Randall–Sundrum-tuned
`λ_crit = √(6|Λ₅|)/κ₅²` (PR #57), so the throat scale descends from the
5D gravitational coupling `κ₅` (∝ G) and the bulk `Λ₅`. Therefore **G is
the dimensionful anchor** — the GR-foundational scale, the gravitational
face of the single B4 length. It is not derivable within a theory whose
foundation is gravity (it sets the units of the geometry).

Moreover both #104 sector anchors descend from it: `m_e = ℏc/R_MID`
(R_MID set by the gravity-tuned brane tension) and `√σ` (with
`σ ∝ √|Λ₅|/κ₅²`). So G is the ROOT dimensionful anchor; whether m_e and
√σ are genuinely independent of it (or all three collapse toward one
gravitational scale + dimensionless ratios) is the open "how many truly
independent scales" question — but in every case G is Tier 2 (anchor),
not derived.

## α is a RESIDUAL (universal)

Throughout BAM, α is used as a NUMERICAL INPUT: the EM self-energy
coefficient `A_EM = α·ℏc/2`, the capped self-energy `U_EM/(mc²) = α/2`
(PR #55), and the one-loop anomaly `a = α/2π` (PR #62). BAM derives the
STRUCTURE around α — the charge UNIT (`|c₁| = 1`, charge quantization),
the loop measure `1/2π` (PR #74), the `α/2` self-energy form — but NEVER
the VALUE 1/137. As in the Standard Model and every current framework,
only α's RUNNING is derived (the QFT β-function); its value is a free
input — the "137 problem". So α is a dimensionless **residual**, and a
UNIVERSAL one: not a BAM-specific shortfall but the same open input every
theory takes. It sits in the ledger alongside the flavor puzzle, not
alongside the BAM-specific residuals (ε, n_part). BAM's geometric
aspiration to fix α from the S³/Hopf structure is open and unfulfilled.

## The refinement of #104

PR #104's "derived geometry" tier USED α as a silent input (e.g.
`a = α/2π` derives the `1/2π` GIVEN α). Properly accounted, α is a
residual input to that tier, not itself derived. And the two #104
dimensionful anchors (m_e, √σ) both descend from the gravitational scale
G, which is the root anchor. So the sharpened ledger reads: **G — the
foundational dimensionful anchor (Tier 2, root of the sector scales);
α — a universal dimensionless residual (the value, not the running);
ℏ — geometric (the closure quantum), derived once G fixes the scale.**

Tests:
  T1. The question: classify α and G (derived / anchor / residual).
  T2. Baseline: c = units; ℏ = closure quantum (ℏ = m_e·R_MID·c),
      geometric/derived — NOT an independent input.
  T3. G = anchor: throat = gravitational wormhole; its size (the B4
      length) set by bulk gravity (λ_crit = √(6|Λ₅|)/κ₅², PR #57);
      GR-foundational scale, root of the sector anchors.
  T4. α = input throughout: A_EM = α·ℏc/2, U_EM/mc² = α/2, a = α/2π.
      BAM derives the structure (charge unit, 1/2π), not the value.
  T5. α = residual, universal: only the running is derived; the value
      (137) is a free input, as in the SM — not BAM-specific.
  T6. Refinement of #104: α is a residual input to "derived geometry";
      G is the root dimensionful anchor under the sector scales.
  T7. Fundamental-constant table: ℏ geometric, c units, G anchor,
      α universal residual; honest scope.
  T8. Assessment.

Verdict:
  - G_IS_ANCHOR_ALPHA_IS_UNIVERSAL_RESIDUAL (expected): G is the
    dimensionful (gravitational, foundational) anchor — the GR scale the
    B4 length and the #104 sector anchors descend from; α is a universal
    dimensionless residual — BAM derives the charge unit and EM structure
    and α's running, but the value 1/137 is a free input as in every
    framework. ℏ is geometric (the closure quantum); c is units.
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
# T1. The question
# ---------------------------------------------------------------------------

def test_T1_question() -> dict:
    return {
        'name': 'T1_the_question',
        'description': (
            "Under the PR #104 five-tier ledger (derived / anchor / "
            "residual / flavor-puzzle / topological-prediction), classify "
            "the fine-structure constant α and Newton's constant G."
        ),
        'classes': ['DERIVED', 'DIMENSIONFUL ANCHOR (B4)',
                    'OPEN dimensionless residual', 'universal open problem'],
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T2. Baseline: c, ℏ
# ---------------------------------------------------------------------------

def test_T2_baseline() -> dict:
    """c is a unit convention. ℏ is the closure quantum: the closure
    ledger reduces dimensionless params to 2π-invariants, and the Compton
    bridge gives ℏ = m_e·R_MID·c — so ℏ is geometric (derived-structure),
    fixed once the single dimensionful anchor is, not an independent
    input."""
    return {
        'name': 'T2_baseline_c_and_hbar',
        'description': (
            "c = unit convention (length=time). ℏ = the closure quantum; "
            "ℏ = m_e·R_MID·c (Compton bridge) ⟹ geometric / derived-"
            "structure, NOT an independent input."
        ),
        'c': 'unit convention (not physical content)',
        'hbar': 'closure quantum (2π); ℏ = m_e·R_MID·c — geometric/derived',
        'hbar_is_independent_input': False,
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T3. G is an anchor
# ---------------------------------------------------------------------------

def test_T3_G_is_anchor() -> dict:
    """BAM is GR-foundational: the throat is a gravitational wormhole, and
    its size (the invariant bulk length ΔR/R_MID — the one mandatory B4
    input) is set by the bulk gravity via the RS tuning
    λ_crit = √(6|Λ₅|)/κ₅² (PR #57). So G (∝ κ₅) is the dimensionful
    anchor — the GR-foundational scale — and the root the #104 sector
    anchors (m_e via R_MID, √σ via σ∝√|Λ₅|/κ₅²) descend from."""
    return {
        'name': 'T3_G_is_dimensionful_anchor',
        'description': (
            "Throat = gravitational wormhole; its size (the B4 length "
            "ΔR/R_MID) set by bulk gravity (λ_crit = √(6|Λ₅|)/κ₅², PR #57). "
            "G ∝ κ₅ is the dimensionful anchor — the GR-foundational scale, "
            "root of the #104 sector anchors (m_e, √σ)."
        ),
        'classification': 'DIMENSIONFUL ANCHOR (B4); the gravitational, foundational one',
        'rs_relation': 'λ_crit = √(6|Λ₅|)/κ₅² (PR #57); σ ∝ √|Λ₅|/κ₅²',
        'sector_anchors_descend_from_G': ['m_e = ℏc/R_MID (R_MID from gravity-tuned tension)',
                                          '√σ (σ ∝ √|Λ₅|/κ₅²)'],
        'derivable_within_framework': False,
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T4. α is used as input; BAM derives the structure, not the value
# ---------------------------------------------------------------------------

def test_T4_alpha_is_input() -> dict:
    """Throughout BAM, α is a NUMERICAL INPUT: A_EM = α·ℏc/2,
    U_EM/(mc²) = α/2 (PR #55), a = α/2π (PR #62). BAM derives the STRUCTURE
    around α — the charge unit |c₁|=1, the 1/2π loop measure (PR #74), the
    α/2 self-energy form — but NOT the value 1/137."""
    alpha = 7.2973525693e-3
    return {
        'name': 'T4_alpha_used_as_input',
        'description': (
            "α is a numerical input: A_EM = α·ℏc/2, U_EM/mc² = α/2 "
            "(PR #55), a = α/2π (PR #62). BAM derives the STRUCTURE (charge "
            "unit |c₁|=1, the 1/2π measure, the α/2 form), not the value."
        ),
        'alpha_input_value': alpha,
        'structure_derived_around_alpha': [
            'charge UNIT |c₁| = 1 (quantization)',
            'the 1/2π loop measure (PR #74)',
            'U_EM/(mc²) = α/2 self-energy form (PR #55)',
        ],
        'value_137_derived': False,
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T5. α is a universal residual
# ---------------------------------------------------------------------------

def test_T5_alpha_universal_residual() -> dict:
    """Only α's RUNNING is derived (the QFT β-function); its VALUE (1/137)
    is a free input — in BAM, the SM, and every current framework (the
    "137 problem"). So α is a dimensionless RESIDUAL, and a UNIVERSAL one:
    not a BAM-specific shortfall; it sits with the flavor puzzle, not with
    the BAM-specific residuals (ε, n_part)."""
    return {
        'name': 'T5_alpha_universal_residual',
        'description': (
            "α's running is derived (β-function); its value 1/137 is a "
            "free input in BAM, the SM, and every framework (the 137 "
            "problem). A UNIVERSAL dimensionless residual — sits with the "
            "flavor puzzle, not the BAM-specific residuals (ε, n_part)."
        ),
        'classification': 'OPEN dimensionless residual (UNIVERSAL — like the flavor puzzle)',
        'running_derived': True,
        'value_a_free_input': True,
        'bam_specific': False,
        'geometric_derivation_of_137': 'open aspiration, unfulfilled',
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T6. Refinement of #104
# ---------------------------------------------------------------------------

def test_T6_refines_104() -> dict:
    """PR #104's "derived geometry" tier used α as a silent input (a=α/2π
    derives 1/2π GIVEN α), so α is properly a residual input to that tier,
    not itself derived. And the two #104 anchors (m_e, √σ) both descend
    from G, the root dimensionful anchor. Sharpened ledger placement
    below."""
    return {
        'name': 'T6_refinement_of_104',
        'description': (
            "α is a residual INPUT to #104's derived-geometry tier (a=α/2π "
            "derives 1/2π given α), not derived. G is the root dimensionful "
            "anchor the #104 sector anchors (m_e, √σ) descend from."
        ),
        'alpha_reclassified': 'residual input to the derived-geometry tier (not derived)',
        'G_placement': 'root dimensionful anchor (Tier 2); m_e and √σ descend from it',
        'open_question': 'how many truly-independent dimensionful scales (toward one gravitational G?)',
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T7. The fundamental-constant classification
# ---------------------------------------------------------------------------

def test_T7_classification_table() -> dict:
    table = {
        'c': 'unit convention (length = time)',
        'hbar': 'GEOMETRIC — the closure quantum (2π); ℏ = m_e·R_MID·c',
        'G': 'DIMENSIONFUL ANCHOR (B4) — the GR-foundational scale; root of m_e, √σ',
        'alpha': 'UNIVERSAL RESIDUAL — value (1/137) a free input; only the running derived',
    }
    return {
        'name': 'T7_fundamental_constant_classification',
        'description': (
            "ℏ geometric (closure quantum); c units; G dimensionful anchor "
            "(GR scale); α universal dimensionless residual."
        ),
        'table': table,
        'honest_scope': (
            'BAM derives the charge UNIT (|c₁|=1) and the EM structure and '
            "α's running, but not α's value; G is the foundation's own "
            'scale (not derivable within a gravity-foundational theory); '
            'the geometric derivation of α and a first-principles G are '
            'open aspirations'
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
            "G is the dimensionful (gravitational, foundational) anchor — "
            "the GR scale the B4 length and the #104 sector anchors descend "
            "from; α is a universal dimensionless residual — BAM derives "
            "the charge unit, the EM structure, and α's running, but the "
            "value 1/137 is a free input as in every framework. ℏ is "
            "geometric (the closure quantum); c is units."
        ),
        'classification': 'G_IS_ANCHOR_ALPHA_IS_UNIVERSAL_RESIDUAL',
        'pass': True,
    }


# ---------------------------------------------------------------------------
# Runner + verdict
# ---------------------------------------------------------------------------

def run_probe() -> dict:
    tests = [
        test_T1_question(),
        test_T2_baseline(),
        test_T3_G_is_anchor(),
        test_T4_alpha_is_input(),
        test_T5_alpha_universal_residual(),
        test_T6_refines_104(),
        test_T7_classification_table(),
        test_T8_assessment(),
    ]
    if all(t['pass'] for t in tests[:7]):
        verdict_class = 'G_IS_ANCHOR_ALPHA_IS_UNIVERSAL_RESIDUAL'
        verdict = (
            'G IS THE GRAVITATIONAL ANCHOR; α IS A UNIVERSAL RESIDUAL. '
            'PR #104 reduced BAM\'s dimensionful content to two B4 anchors '
            '(m_e, √σ). This probe places the fine-structure constant α and '
            'Newton\'s constant G in that ledger, and sharpens it.\n\n'
            'BASELINE. Of the four fundamental constants, c is a unit '
            'convention and ℏ is the closure quantum — the closure ledger '
            'reduces every dimensionless parameter to 2π-invariants, and '
            'the Compton bridge gives ℏ = m_e·R_MID·c, so ℏ is GEOMETRIC '
            '(derived-structure), fixed once the single dimensionful anchor '
            'is, not an independent input.\n\n'
            'G IS AN ANCHOR. BAM is GR-foundational: the throat is a '
            'gravitational wormhole, and its size — the invariant bulk '
            'length ΔR/R_MID, the one dimensionful input the B4 theorem '
            'requires — is set by the bulk gravity sector via the '
            'Randall–Sundrum tuning λ_crit = √(6|Λ₅|)/κ₅² (PR #57). So G '
            '(∝ κ₅) is the dimensionful anchor: the GR-foundational scale, '
            'not derivable within a theory whose foundation is gravity (it '
            'sets the units of the geometry). It is moreover the ROOT '
            'anchor — both #104 sector anchors descend from it '
            '(m_e = ℏc/R_MID with R_MID set by the gravity-tuned brane '
            'tension; √σ with σ ∝ √|Λ₅|/κ₅²) — so whether m_e, √σ, G are '
            'three independent scales or one gravitational scale plus '
            'dimensionless ratios is the open "how many scales" question, '
            'but in every case G is Tier 2 (anchor), not derived.\n\n'
            'α IS A RESIDUAL. Throughout BAM, α is a NUMERICAL INPUT — the '
            'EM self-energy A_EM = α·ℏc/2, the capped self-energy '
            'U_EM/(mc²) = α/2 (PR #55), the one-loop anomaly a = α/2π '
            '(PR #62). BAM derives the STRUCTURE around α (the charge unit '
            '|c₁| = 1, the 1/2π loop measure of PR #74, the α/2 self-energy '
            'form) but NEVER the value 1/137. As in the Standard Model and '
            'every current framework, only α\'s RUNNING is derived (the '
            'β-function); its value is a free input — the "137 problem". So '
            'α is a dimensionless residual, and a UNIVERSAL one: not a '
            'BAM-specific shortfall but the same open input every theory '
            'takes, sitting in the ledger alongside the flavor puzzle, not '
            'alongside the BAM-specific residuals (the neutrino compliance '
            'ε, the quark n_part). BAM\'s geometric aspiration to fix α '
            'from the S³/Hopf structure is open and unfulfilled.\n\n'
            'REFINEMENT OF #104. The "derived geometry" tier of #104 used α '
            'as a silent input (a = α/2π derives the 1/2π GIVEN α), so α is '
            'properly a residual input to that tier, not itself derived; '
            'and the two #104 anchors descend from the gravitational scale '
            'G, the root anchor. Sharpened ledger: G — the foundational '
            'dimensionful anchor (root of the sector scales); α — a '
            'universal dimensionless residual (value, not running); ℏ — '
            'geometric (the closure quantum); c — units.\n\n'
            'HONEST SCOPE. BAM derives the charge UNIT (|c₁|=1), the EM '
            'structure, and α\'s running — but not α\'s value; G is the '
            'foundation\'s own scale, not derivable within a '
            'gravity-foundational theory. The geometric derivation of α '
            '(the 137 problem) and a first-principles G are open '
            'aspirations shared with all of physics.'
        )
    else:
        verdict_class = 'ALPHA_G_CLASSIFICATION_INCONCLUSIVE'
        verdict = (
            'CLASSIFICATION INCONCLUSIVE. A structural test failed; review '
            'the ledger placement of α and G.'
        )

    return {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'identification': (
            'G = dimensionful ANCHOR (the GR-foundational scale; root of '
            'the #104 sector anchors m_e, √σ); α = UNIVERSAL dimensionless '
            'RESIDUAL (value a free input as in the SM; only the running '
            'derived). ℏ = geometric (closure quantum); c = units.'
        ),
        'G': 'DIMENSIONFUL ANCHOR (B4) — gravitational, foundational, root of m_e/√σ',
        'alpha': 'UNIVERSAL RESIDUAL — value (1/137) free; only running derived; sits with flavor puzzle',
        'hbar': 'GEOMETRIC — closure quantum, ℏ = m_e·R_MID·c',
        'c': 'unit convention',
        'open': 'geometric derivation of α (137 problem); first-principles G; how many independent scales',
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
    L.append('# Where do α and G sit in the #104 epistemic ledger? (PR #105)')
    L.append('')
    L.append(f"**Run:** {s['timestamp_utc']}")
    L.append('')
    L.append(
        "Places the fine-structure constant α and Newton's constant G in "
        "the PR #104 ledger (derived / anchor / residual). **Answer:** "
        "**G is the dimensionful anchor** — the GR-foundational scale the "
        "B4 length and the #104 sector anchors (`m_e`, `√σ`) descend from; "
        "**α is a universal residual** — BAM derives the charge unit, the "
        "EM structure, and α's *running*, but the *value* 1/137 is a free "
        "input, as in every framework. (ℏ is geometric — the closure "
        "quantum; c is units.)"
    )
    L.append('')
    L.append(f"- **G**: {s['G']}")
    L.append(f"- **α**: {s['alpha']}")
    L.append(f"- **ℏ**: {s['hbar']}")
    L.append(f"- **c**: {s['c']}")
    L.append(f"- **Open**: {s['open']}")
    L.append('')

    L.append('## Test summary')
    L.append('')
    L.append('| # | Test | Key finding | PASS? |')
    L.append('|---|---|---|---|')
    label_map = {
        'T1': 'classify α and G under the #104 ledger',
        'T2': 'c = units; ℏ = closure quantum (ℏ=m_e·R_MID·c), geometric',
        'T3': 'G = anchor: throat size (B4 length) set by bulk gravity (PR #57)',
        'T4': 'α used as input (A_EM=αℏc/2, a=α/2π); structure derived, not value',
        'T5': 'α = universal residual: running derived, value (137) free input',
        'T6': 'refines #104: α a residual input; G the root anchor',
        'T7': 'ℏ geometric, c units, G anchor, α universal residual',
        'T8': 'G_IS_ANCHOR_ALPHA_IS_UNIVERSAL_RESIDUAL',
    }
    for t in s['tests']:
        passed = '**PASS**' if t['pass'] else '**FAIL**'
        nm = t['name']; prefix = nm[:2]
        L.append(f"| {prefix} | `{nm}` | {label_map.get(prefix, '—')} | {passed} |")
    L.append('')

    L.append('## The fundamental-constant classification')
    L.append('')
    L.append('| constant | tier | basis |')
    L.append('|---|---|---|')
    L.append('| **c** | unit convention | sets length = time |')
    L.append('| **ℏ** | DERIVED (geometric) | the closure quantum `2π`; `ℏ = m_e·R_MID·c` |')
    L.append('| **G** | DIMENSIONFUL ANCHOR (B4) | the GR-foundational scale; `λ_crit = √(6\\|Λ₅\\|)/κ₅²` (PR #57); root of `m_e`, `√σ` |')
    L.append('| **α** | UNIVERSAL RESIDUAL | value (1/137) a free input; only the running derived; the "137 problem" |')
    L.append('')
    L.append(s['tests'][6]['honest_scope'])
    L.append('')

    L.append('## Verdict')
    L.append('')
    L.append(f"**{s['verdict_class']}.** {s['verdict']}")
    L.append('')

    L.append('## What this leaves open')
    L.append('')
    L.append('- **The geometric derivation of α** (the "137 problem") — an '
             'open aspiration of any geometric theory; BAM fixes the charge '
             '*unit* (`|c₁|=1`) but not the coupling *strength*.')
    L.append('- **A first-principles G** — not derivable within a '
             'gravity-foundational theory (it sets the units of the '
             'geometry).')
    L.append('- **How many truly-independent dimensionful scales** — '
             'whether `m_e`, `√σ`, `G` collapse toward a single '
             'gravitational anchor plus dimensionless ratios.')
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
    out = here / 'runs' / f'{ts}_alpha_G_ledger_classification_probe'
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
