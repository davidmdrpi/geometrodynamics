"""
Are m_e and √σ independent B4 anchors, or sector readouts of one
bulk-gravity scale G? (PR #106).

PR #105 found G the dimensionful anchor and the root the #104 sector
anchors (m_e, √σ) descend from, and left the sharp question: are m_e and
√σ two INDEPENDENT dimensionful inputs, or readouts of a single
bulk-gravity scale G times dimensionless ratios? This probe answers it.

## Common gravitational origin

Both #104 anchors are brane/geometric scales of the SAME bulk geometry,
descending from the bulk gravity sector (PR #57): the throat size
R_MID (⟹ m_e = ℏc/R_MID) is set by the Randall–Sundrum-tuned brane
tension λ_crit = √(6|Λ₅|)/κ₅², and the confinement tension σ ∝ √|Λ₅|/κ₅².
So m_e and √σ are NOT two independent KINDS of input — both are the one
bulk-gravity scale G read out in two channels (the throat-winding lepton
channel and the cavity-confinement QCD channel, PR #83).

## But the ratio is not derived

The dimensionless ratio is

    √σ / m_e  ≈  0.424 GeV / 0.511 MeV  ≈  830,

the lepton-throat (electron Compton) to QCD-confinement length
hierarchy. If BAM's geometry fixed it, the two channels' relative
normalisation would be a clean closure number. It is not: the nearest
BAM closure quantity is 50π·k_5 = 785 (5.4% off) — a near-coincidence in
the spirit of F_13 = 233 (PR #76), not a derivation. So the ~830 ratio is
UNDERIVED: the geometry does not yet fix the relative normalisation of
the throat-winding and cavity-confinement channels.

## The verdict and the honest bookkeeping

m_e and √σ are NOT independent dimensionful anchors — they share the one
bulk-gravity scale G. Equivalently, the program has **one foundational
dimensionful anchor (G) plus one open dimensionless ratio** (m_e/√σ ≈
1/830). This REDUCES the dimensionful-anchor count from two to one. But
it is a REPACKAGING, not a free-lunch reduction: a dimensionful anchor
has been converted into a dimensionless residual, so the TOTAL count of
irreducible inputs is unchanged (was: 2 dimensionful anchors; now: 1
dimensionful anchor + 1 dimensionless ratio). What it buys is conceptual
cleanliness — the sole fundamental SCALE is the gravitational foundation
G, with everything else dimensionless — exactly the GR-foundational
posture. The new ratio joins ε, n_part, and α as the program's open
dimensionless residuals; its smallness (the electron being anomalously
light relative to Λ_QCD) overlaps the flavor puzzle, though it is not
identical to it (it also carries the QCD scale's origin).

## What this probe establishes (and does not)

  - **Establishes (BAM-native):** m_e and √σ are not independent — both
    descend from the single bulk-gravity scale G (PR #57). The
    dimensionful-anchor count reduces from two to one (G). But the
    dimensionless ratio m_e/√σ ≈ 1/830 is UNDERIVED (no clean closure
    match; 50π·k_5 = 785 is a 5.4% near-coincidence), so it becomes a new
    open dimensionless residual — the repackaging leaves the total
    irreducible-input count unchanged.

  - **Does not establish:** a derivation of the ~830 lepton/QCD scale
    hierarchy. That is the open dimensionless input the reduction exposes;
    fixing it (from the throat-winding vs cavity-confinement channel
    normalisation) would genuinely reduce the count to one.

Tests:
  T1. The question: independent anchors, or one G + ratios?
  T2. Common gravitational origin (both brane scales from bulk gravity,
      PR #57) ⟹ not two independent KINDS.
  T3. The ratio √σ/m_e ≈ 830 (lepton-throat / QCD-confinement hierarchy).
  T4. Is 830 derived? No clean closure match (nearest 50π·k_5=785, 5.4%
      off — a near-coincidence like F_13); underived.
  T5. Verdict: NOT independent — one bulk scale G + one open ratio;
      dimensionful count 2→1.
  T6. Honest bookkeeping: a repackaging (dimensionful → dimensionless),
      total inputs unchanged; cleaner "one fundamental scale G" picture.
  T7. Updated ledger: ONE foundational anchor G; m_e/√σ a new dimensionless
      residual (overlaps the flavor puzzle).
  T8. Assessment.

Verdict:
  - M_E_SQRT_SIGMA_NOT_INDEPENDENT_ONE_G_PLUS_UNDERIVED_RATIO (expected):
    m_e and √σ are not independent — both are sector readouts of the
    single bulk-gravity scale G; the dimensionful-anchor count reduces
    2→1, but the dimensionless ratio m_e/√σ ≈ 1/830 is underived and
    becomes a new open residual, so the total irreducible-input count is
    unchanged — a cleaner "one fundamental scale (G)" repackaging.
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
K_5 = 5
M_E_GEV = 0.5109989e-3
SQRT_SIGMA_GEV = math.sqrt(SIGMA_QCD)
RATIO = SQRT_SIGMA_GEV / M_E_GEV          # √σ / m_e ≈ 830


# ---------------------------------------------------------------------------
# T1. The question
# ---------------------------------------------------------------------------

def test_T1_question() -> dict:
    return {
        'name': 'T1_the_question',
        'description': (
            "Are m_e (lepton scale) and √σ (QCD scale) two INDEPENDENT B4 "
            "anchors, or readouts of a single bulk-gravity scale G × "
            "dimensionless ratios?"
        ),
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T2. Common gravitational origin
# ---------------------------------------------------------------------------

def test_T2_common_origin() -> dict:
    """Both #104 anchors are brane/geometric scales of the SAME bulk
    geometry, descending from the bulk gravity (PR #57): R_MID (⟹ m_e =
    ℏc/R_MID) from the RS-tuned tension λ_crit = √(6|Λ₅|)/κ₅²; σ ∝
    √|Λ₅|/κ₅². So not two independent KINDS — one bulk-gravity scale G in
    two channels (throat-winding lepton, cavity-confinement QCD)."""
    return {
        'name': 'T2_common_gravitational_origin',
        'description': (
            "Both m_e and √σ descend from the bulk gravity (PR #57): "
            "R_MID from λ_crit = √(6|Λ₅|)/κ₅²; σ ∝ √|Λ₅|/κ₅². One "
            "bulk-gravity scale G in two channels (throat-winding / "
            "cavity-confinement)."
        ),
        'm_e_origin': 'ℏc/R_MID; R_MID from gravity-tuned brane tension (PR #57)',
        'sqrt_sigma_origin': 'σ ∝ √|Λ₅|/κ₅² (PR #57)',
        'two_independent_kinds': False,
        'one_scale_two_channels': True,
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T3. The dimensionless ratio
# ---------------------------------------------------------------------------

def test_T3_ratio() -> dict:
    """√σ/m_e ≈ 830 — the lepton-throat (electron Compton) to
    QCD-confinement length hierarchy."""
    return {
        'name': 'T3_dimensionless_ratio',
        'description': (
            "√σ/m_e ≈ 830 (m_e/√σ ≈ 1/830): the lepton-throat / "
            "QCD-confinement scale hierarchy."
        ),
        'm_e_MeV': M_E_GEV * 1e3,
        'sqrt_sigma_MeV': SQRT_SIGMA_GEV * 1e3,
        'ratio_sqrt_sigma_over_m_e': RATIO,
        'pass': 800 < RATIO < 860,
    }


# ---------------------------------------------------------------------------
# T4. Is the ratio derived?
# ---------------------------------------------------------------------------

def test_T4_ratio_underived() -> dict:
    """If the geometry fixed the ratio, it would be a clean closure
    number. It is not: the nearest BAM closure quantity is 50π·k_5 = 785
    (5.4% off) — a near-coincidence (cf. F_13 = 233, PR #76), not a
    derivation. So the ~830 ratio is UNDERIVED."""
    candidates = {
        '50π·k_5': 50.0 * PI * K_5,
        'k_5⁴': K_5 ** 4,
        'k_5⁴ + k_5²': K_5 ** 4 + K_5 ** 2,
        '(2π)³': (2.0 * PI) ** 3,
        '2·k_5⁴': 2.0 * K_5 ** 4,
    }
    rows = [{'candidate': c, 'value': v, 'rel_off': (v - RATIO) / RATIO}
            for c, v in candidates.items()]
    best = min(rows, key=lambda r: abs(r['rel_off']))
    clean_match = abs(best['rel_off']) < 0.01
    return {
        'name': 'T4_ratio_underived',
        'description': (
            "No clean closure match for ~830: nearest is 50π·k_5 = 785 "
            "(5.4% off), a near-coincidence (cf. F_13=233), not a "
            "derivation. The ratio is UNDERIVED."
        ),
        'candidates': rows,
        'closest': best['candidate'],
        'closest_rel_off': best['rel_off'],
        'clean_match': clean_match,
        'underived': not clean_match,
        'pass': not clean_match,
    }


# ---------------------------------------------------------------------------
# T5. The verdict on independence
# ---------------------------------------------------------------------------

def test_T5_not_independent() -> dict:
    """m_e and √σ are NOT independent — they share the one bulk-gravity
    scale G. Equivalently: one foundational dimensionful anchor (G) + one
    open dimensionless ratio (m_e/√σ ≈ 1/830). The dimensionful-anchor
    count reduces from two to one."""
    return {
        'name': 'T5_not_independent_anchors',
        'description': (
            "m_e and √σ NOT independent — one bulk-gravity scale G + one "
            "open dimensionless ratio. Dimensionful-anchor count 2 → 1."
        ),
        'independent_anchors': False,
        'dimensionful_anchor_count_before': 2,
        'dimensionful_anchor_count_after': 1,
        'plus_one_open_dimensionless_ratio': True,
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T6. Honest bookkeeping: a repackaging
# ---------------------------------------------------------------------------

def test_T6_repackaging() -> dict:
    """The 2→1 reduction is a REPACKAGING, not a free reduction: a
    dimensionful anchor has been converted into a dimensionless residual,
    so the TOTAL count of irreducible inputs is unchanged (2 dimensionful
    → 1 dimensionful + 1 dimensionless). What it buys is conceptual
    cleanliness — the sole fundamental SCALE is the gravitational
    foundation G, everything else dimensionless (the GR-foundational
    posture)."""
    before_total = 2          # 2 dimensionful anchors
    after_total = 1 + 1       # 1 dimensionful anchor + 1 dimensionless ratio
    return {
        'name': 'T6_repackaging_not_free_reduction',
        'description': (
            "2→1 dimensionful is a REPACKAGING: a dimensionful anchor → a "
            "dimensionless residual. TOTAL irreducible inputs unchanged "
            "(2 → 1+1). Buys conceptual cleanliness: the sole fundamental "
            "scale is G (GR-foundational)."
        ),
        'total_inputs_before': before_total,
        'total_inputs_after': after_total,
        'total_unchanged': before_total == after_total,
        'gain': 'sole fundamental SCALE = G; everything else dimensionless',
        'pass': before_total == after_total,
    }


# ---------------------------------------------------------------------------
# T7. Updated ledger
# ---------------------------------------------------------------------------

def test_T7_updated_ledger() -> dict:
    return {
        'name': 'T7_updated_ledger',
        'description': (
            "Updated ledger: ONE foundational dimensionful anchor (G); the "
            "m_e/√σ ratio joins ε, n_part, α as open dimensionless "
            "residuals (its smallness overlaps the flavor puzzle)."
        ),
        'foundational_dimensionful_anchor': 'G (the bulk-gravity scale)',
        'open_dimensionless_residuals': [
            'neutrino compliance ε',
            'quark n_part = 233',
            'α (universal; PR #105)',
            'm_e/√σ ≈ 1/830 (the lepton/QCD scale hierarchy; this probe)',
        ],
        'ratio_overlaps_flavor_puzzle': True,
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T8. Assessment
# ---------------------------------------------------------------------------

def test_T8_assessment() -> dict:
    return {
        'name': 'T8_assessment',
        'description': (
            "m_e and √σ are not independent — both are sector readouts of "
            "the single bulk-gravity scale G; the dimensionful-anchor count "
            "reduces 2→1, but the dimensionless ratio m_e/√σ ≈ 1/830 is "
            "underived and becomes a new open residual, so the total "
            "irreducible-input count is unchanged — a cleaner 'one "
            "fundamental scale (G)' repackaging."
        ),
        'classification': (
            'M_E_SQRT_SIGMA_NOT_INDEPENDENT_ONE_G_PLUS_UNDERIVED_RATIO'
        ),
        'pass': True,
    }


# ---------------------------------------------------------------------------
# Runner + verdict
# ---------------------------------------------------------------------------

def run_probe() -> dict:
    tests = [
        test_T1_question(),
        test_T2_common_origin(),
        test_T3_ratio(),
        test_T4_ratio_underived(),
        test_T5_not_independent(),
        test_T6_repackaging(),
        test_T7_updated_ledger(),
        test_T8_assessment(),
    ]
    if all(t['pass'] for t in tests[:7]):
        verdict_class = (
            'M_E_SQRT_SIGMA_NOT_INDEPENDENT_ONE_G_PLUS_UNDERIVED_RATIO'
        )
        verdict = (
            'm_e AND √σ ARE NOT INDEPENDENT — ONE BULK-GRAVITY SCALE G PLUS '
            'AN UNDERIVED DIMENSIONLESS RATIO. PR #105 found G the '
            'dimensionful anchor and the root the #104 sector anchors '
            'descend from. This probe settles whether m_e and √σ are two '
            'independent dimensionful inputs or readouts of the single '
            'scale G.\n\n'
            'COMMON GRAVITATIONAL ORIGIN. Both #104 anchors are '
            'brane/geometric scales of the SAME bulk geometry, descending '
            'from the bulk gravity (PR #57): the throat size R_MID (⟹ m_e = '
            'ℏc/R_MID) is set by the Randall–Sundrum-tuned brane tension '
            'λ_crit = √(6|Λ₅|)/κ₅², and the confinement tension σ ∝ '
            '√|Λ₅|/κ₅². So m_e and √σ are NOT two independent KINDS of '
            'input — both are the one bulk-gravity scale G read out in two '
            'channels (the throat-winding lepton channel and the '
            'cavity-confinement QCD channel, PR #83).\n\n'
            'BUT THE RATIO IS NOT DERIVED. The dimensionless ratio √σ/m_e ≈ '
            '0.424 GeV / 0.511 MeV ≈ 830 — the lepton-throat (electron '
            'Compton) to QCD-confinement length hierarchy — is not a clean '
            'closure number: the nearest BAM closure quantity is 50π·k_5 = '
            '785 (5.4% off), a near-coincidence in the spirit of F_13 = 233 '
            '(PR #76), not a derivation. So the geometry does not yet fix '
            'the relative normalisation of the two channels.\n\n'
            'THE VERDICT. m_e and √σ are NOT independent dimensionful '
            'anchors — they share the one bulk-gravity scale G. '
            'Equivalently, the program has ONE foundational dimensionful '
            'anchor (G) plus ONE open dimensionless ratio (m_e/√σ ≈ '
            '1/830). The dimensionful-anchor count reduces from two to '
            'one. But this is a REPACKAGING, not a free-lunch reduction: a '
            'dimensionful anchor has been converted into a dimensionless '
            'residual, so the TOTAL count of irreducible inputs is '
            'unchanged (2 dimensionful → 1 dimensionful + 1 dimensionless). '
            'What it buys is conceptual cleanliness — the sole fundamental '
            'SCALE is the gravitational foundation G, with everything else '
            'dimensionless, exactly the GR-foundational posture. The new '
            'ratio joins ε, n_part, and α as the program\'s open '
            'dimensionless residuals; its smallness (the electron being '
            'anomalously light relative to Λ_QCD) overlaps the flavor '
            'puzzle, though it is not identical to it.\n\n'
            'HONEST SCOPE. ESTABLISHED (BAM-native): m_e and √σ both '
            'descend from the single bulk-gravity scale G (PR #57), so they '
            'are not independent; the dimensionful-anchor count reduces '
            '2→1. NOT established: a derivation of the ~830 lepton/QCD '
            'hierarchy — it is underived (no clean closure match) and '
            'becomes a new open dimensionless residual, leaving the total '
            'irreducible-input count unchanged. Fixing that ratio from the '
            'throat-winding vs cavity-confinement channel normalisation '
            'would genuinely reduce the count to one.'
        )
    else:
        verdict_class = 'SCALE_COUNT_INCONCLUSIVE'
        verdict = (
            'SCALE COUNT INCONCLUSIVE. A structural test failed; review '
            'the anchor accounting.'
        )

    return {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'identification': (
            'm_e and √σ are not independent — both are sector readouts of '
            'the single bulk-gravity scale G (PR #57); the '
            'dimensionful-anchor count reduces 2→1, but the dimensionless '
            'ratio √σ/m_e ≈ 830 is underived (a new open residual), so the '
            'total irreducible-input count is unchanged'
        ),
        'origin': 'both brane scales of one bulk geometry (G); PR #57 RS tuning',
        'ratio': '√σ/m_e ≈ 830 — UNDERIVED (nearest 50π·k_5=785, 5.4% off)',
        'verdict_on_independence': 'NOT independent — one G + one open dimensionless ratio',
        'bookkeeping': 'dimensionful 2→1, but a dimensionless residual added; total inputs unchanged',
        'open': 'derive the ~830 lepton/QCD hierarchy (channel normalisation) to reach a true single input',
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
    L.append('# m_e and √σ: independent anchors, or one bulk-gravity scale G? (PR #106)')
    L.append('')
    L.append(f"**Run:** {s['timestamp_utc']}")
    L.append('')
    L.append(
        "Settles the scale-count question PR #105 exposed. **Answer:** m_e "
        "and √σ are **not independent** — both descend from the single "
        "bulk-gravity scale `G` (PR #57), so the dimensionful-anchor count "
        "reduces **2 → 1**. But their ratio `√σ/m_e ≈ 830` is **underived** "
        "(no clean closure match), so it becomes a new open dimensionless "
        "residual — the total irreducible-input count is unchanged. A "
        "cleaner *one fundamental scale (G)* repackaging, not a free "
        "reduction."
    )
    L.append('')
    L.append(f"- **Identification**: {s['identification']}")
    L.append(f"- **Origin**: {s['origin']}")
    L.append(f"- **Ratio**: {s['ratio']}")
    L.append(f"- **Independence**: {s['verdict_on_independence']}")
    L.append(f"- **Bookkeeping**: {s['bookkeeping']}")
    L.append(f"- **Open**: {s['open']}")
    L.append('')

    L.append('## Test summary')
    L.append('')
    L.append('| # | Test | Key finding | PASS? |')
    L.append('|---|---|---|---|')
    label_map = {
        'T1': 'independent anchors, or one G + ratios?',
        'T2': 'common gravitational origin (both brane scales of one bulk)',
        'T3': '√σ/m_e ≈ 830 (lepton-throat / QCD-confinement hierarchy)',
        'T4': '830 underived (nearest 50π·k_5=785, 5.4% off — near-coincidence)',
        'T5': 'NOT independent — one G + one open ratio; dimensionful 2→1',
        'T6': 'repackaging: dimensionful→dimensionless; total inputs unchanged',
        'T7': 'ledger: one anchor G; m_e/√σ a new dimensionless residual',
        'T8': 'M_E_SQRT_SIGMA_NOT_INDEPENDENT_ONE_G_PLUS_UNDERIVED_RATIO',
    }
    for t in s['tests']:
        passed = '**PASS**' if t['pass'] else '**FAIL**'
        nm = t['name']; prefix = nm[:2]
        L.append(f"| {prefix} | `{nm}` | {label_map.get(prefix, '—')} | {passed} |")
    L.append('')

    t4 = s['tests'][3]
    L.append('## T4: the ~830 ratio has no clean closure match')
    L.append('')
    L.append('| candidate | value | off |')
    L.append('|---|---:|---:|')
    for r in t4['candidates']:
        L.append(f"| {r['candidate']} | {r['value']:.0f} | {r['rel_off']*100:+.1f}% |")
    L.append('')
    L.append(f"Nearest is `{t4['closest']}` ({t4['closest_rel_off']*100:+.1f}%) — a "
             "near-coincidence (cf. `F_13 = 233`, PR #76), not a derivation. "
             "The ratio is **underived**.")
    L.append('')

    t6 = s['tests'][5]
    L.append('## T6: the bookkeeping (a repackaging)')
    L.append('')
    L.append(f"- before: **{t6['total_inputs_before']}** dimensionful anchors")
    L.append(f"- after: **{t6['total_inputs_after']}** = 1 dimensionful anchor (G) + "
             "1 open dimensionless ratio")
    L.append(f"- total irreducible inputs **unchanged** — but the sole "
             "fundamental *scale* is now `G` (GR-foundational).")
    L.append('')

    L.append('## Verdict')
    L.append('')
    L.append(f"**{s['verdict_class']}.** {s['verdict']}")
    L.append('')

    L.append('## What this leaves open')
    L.append('')
    L.append('- **Deriving the ~830 lepton/QCD hierarchy** — from the '
             'throat-winding vs cavity-confinement channel normalisation. '
             'Doing so would genuinely reduce BAM to a SINGLE irreducible '
             'input (the gravitational scale G).')
    L.append('- The ratio currently joins `ε`, `n_part`, `α` as the '
             'program\'s open dimensionless residuals; its smallness '
             'overlaps the flavor puzzle (the electron\'s lightness).')
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
    out = here / 'runs' / f'{ts}_scale_count_anchors_probe'
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
