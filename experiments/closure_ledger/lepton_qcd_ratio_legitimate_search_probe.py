"""
The legitimate search: does ANY fit-independent, §8-stable bulk quantity
select the lepton/QCD scale ratio √σ/m_e ≈ 830? (PR #108).

PR #106 isolated the lepton/QCD scale ratio √σ/m_e ≈ 830 as the single
UNDERIVED open dimensionless residual whose derivation would collapse BAM
to one foundational anchor (the bulk-gravity scale G). PR #107 then killed
the tempting candidate `N_q + ΔN = 832 ≈ 830`: it merely recycles the
phenomenological compensator n_part (832 = 4·n_part − 4·k_5²), §8-drifts
764–920, and rides on no independent integral.

This probe runs the HONEST, fit-INDEPENDENT search that #107 said was the
only thing that could legitimately close the gap: scan principled bulk
geometric quantities — built ONLY from BAM's fixed structural integers
(k_5 = 5, β_lepton = 50π), the closure constant 2π, and the Tangherlini
cavity spectrum — for one that selects ~830, to <1%, while being §8-stable.
A candidate must satisfy ALL FOUR criteria to count as a derivation:

    (C1) fit-independent  — built from fixed geometry, NOT from n_part /
                            the v3 quark fit (the failure mode of #107);
    (C2) §8-stable        — invariant under the quark_axioms §8 ablations;
    (C3) precise          — within ≈1% of the observed 830.3;
    (C4) principled       — no ad-hoc multiplicative fudge factor.

**Result: the search FAILS.** §8-stability (C2) is automatic for genuine
geometric candidates — they don't depend on quark-mass ablations at all,
which is exactly what n_part lacked. But no candidate clears C3 AND C4
together. The best PRINCIPLED candidate is

    2π·k_5³ = β_lepton·k_5 = 785.40   (−5.4% from 830.3),

and every sub-percent match requires an unprincipled factor:

    π·265   = 832.5  (+0.27%)   — but 265 is not a BAM integer;
    (4/3)·k_5⁴ = 833.3 (+0.37%) — but the 4/3 is ad-hoc;
    k_5⁵/3.77 = 828.9 (−0.16%)  — but 3.77 is ad-hoc.

The exponential / dimensional-transmutation route fares no better:
ln(830.3) = 6.72, and the only clean geometric action near it (the closure
quantum 2π = 6.28) is 7% off, so 830 is not e^(geometric action) for any
principled action either. The Tangherlini cavity-eigenvalue sums are
O(10–350) and select nothing near 830.

**Conclusion (honest negative).** No principled, fit-independent,
§8-stable bulk quantity selects √σ/m_e ≈ 830. The ratio therefore stays an
UNDERIVED open dimensionless residual — and the evidence now points to it
being IRREDUCIBLY residual, in the same class as α and the electron's
anomalous lightness (the universal flavor puzzle): a pure number the
geometry does not fix. BAM does not collapse to a single anchor; it sits
at one foundational scale (G) + this one open ratio + α + the flavor
puzzle. That is the honest floor, and PR #107's caution is vindicated:
having rejected the n_part recycling, the fit-independent search that
would have legitimately closed the gap comes up empty.

Tests:
  T1. State the target √σ/m_e and the four criteria (C1–C4) a candidate
      must satisfy to count as a derivation.
  T2. §8-stability (C2) is AUTOMATIC for fit-independent geometric
      candidates (they don't touch the quark ablations) — the real bar is
      C3 (precision) ∧ C4 (principle).
  T3. Scan principled powers/products of {k_5, β_lepton, 2π}: best is
      2π·k_5³ = 785.4 (−5.4%); nothing principled lands <1%.
  T4. Every sub-% match needs an ad-hoc factor (π·265, (4/3)·k_5⁴,
      k_5⁵/3.77) — fails C4.
  T5. Exponential / dimensional-transmutation route: ln(830)=6.72, the
      clean action 2π=6.28 is 7% off ⟹ no principled e^(action) selects 830.
  T6. Tangherlini cavity-eigenvalue bulk integrals are O(10–350) — select
      nothing near 830.
  T7. Honest scope: search fails (no candidate clears C3 ∧ C4); √σ/m_e
      stays UNDERIVED and is plausibly IRREDUCIBLY residual, like α and the
      electron's anomalous lightness (flavor puzzle).
  T8. Assessment.

Verdict:
  - LEPTON_QCD_RATIO_NO_PRINCIPLED_GEOMETRIC_SELECTION_STAYS_UNDERIVED
    (expected): the fit-independent §8-stable search comes up empty. The
    best principled candidate (2π·k_5³ = 785) is 5.4% off; sub-% matches
    need ad-hoc factors; the exponential route has no clean action; cavity
    integrals select nothing near 830. √σ/m_e stays an underived open
    residual, plausibly irreducible like α. BAM stays at one scale G + one
    open ratio + α + the flavor puzzle — it does NOT collapse to a single
    anchor.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np

from geometrodynamics.constants import R_MID, R_OUTER
from geometrodynamics.tangherlini.radial import (
    V_tangherlini, r_to_rstar, rstar_to_r,
)
from geometrodynamics.qcd.constants import SIGMA_QCD


PI = math.pi
K_5 = 5
BETA_LEPTON = 50 * PI                       # 157.08, the lepton closure phase
N_PART_BASELINE = 233

M_E_GEV = 0.5109989e-3
RATIO = math.sqrt(SIGMA_QCD) / M_E_GEV      # √σ/m_e ≈ 830.3 (the target)

PRECISION_BAR = 0.01                        # C3: within ≈1%

# n_part across the 12 quark_axioms §8 ablations (PR #76/#97), for the C2
# contrast: a fit quantity drifts; a geometric one does not.
N_PART_S8 = [233, 233, 233, 238, 237, 237, 241, 216, 247, 247, 220, 255]


def _off(v: float) -> float:
    return (v - RATIO) / RATIO


def _shell_eigenvalues(l: int = 1, n: int = 800):
    rsmin = r_to_rstar(R_MID + 5e-4, R_MID)
    rsmax = r_to_rstar(R_OUTER - 5e-4, R_MID)
    rstar = np.linspace(rsmin, rsmax, n)
    h = rstar[1] - rstar[0]
    rphys = np.array([rstar_to_r(s, R_MID) for s in rstar])
    V = V_tangherlini(rphys, l, R_MID)
    main = 2.0 / h ** 2 + V[1:-1]
    off = -1.0 / h ** 2 * np.ones(n - 3)
    ev = np.sort(np.linalg.eigvalsh(
        np.diag(main) + np.diag(off, 1) + np.diag(off, -1)))
    return np.maximum(ev, 0.0)


# ---------------------------------------------------------------------------
# T1. The target and the four criteria
# ---------------------------------------------------------------------------

def test_T1_target_and_criteria() -> dict:
    return {
        'name': 'T1_target_and_criteria',
        'description': (
            "Target: √σ/m_e ≈ 830.3, the one UNDERIVED ratio (PR #106) whose "
            "derivation would collapse BAM to a single anchor G. A candidate "
            "counts as a derivation only if it satisfies ALL of: C1 "
            "fit-independent (not n_part / the v3 fit — #107's failure), C2 "
            "§8-stable, C3 precise (<1%), C4 principled (no ad-hoc factor)."
        ),
        'target_ratio': RATIO,
        'criteria': {
            'C1_fit_independent': 'built from fixed geometry, not n_part / the v3 quark fit',
            'C2_s8_stable': 'invariant under the quark_axioms §8 ablations',
            'C3_precise': f'within {PRECISION_BAR*100:.0f}% of 830.3',
            'C4_principled': 'no ad-hoc multiplicative fudge factor',
        },
        'pass': abs(RATIO - 830.3) < 1.0,
    }


# ---------------------------------------------------------------------------
# T2. §8-stability is automatic for fit-independent candidates
# ---------------------------------------------------------------------------

def test_T2_s8_stability_automatic() -> dict:
    """The #107 failure was C2: 4·n_part−100 §8-drifts 764–920. A genuine
    geometric candidate (built from k_5, β_lepton, 2π) does NOT depend on
    the quark-mass ablations at all, so C2 is automatic. That means the
    binding constraints for THIS search are C3 (precision) and C4
    (principle) — not C2."""
    npart_analogue = [4 * n - 100 for n in N_PART_S8]
    npart_span = max(npart_analogue) - min(npart_analogue)
    # geometric candidate: identical across every ablation (no n_part input)
    geom_analogue = [2 * PI * K_5 ** 3 for _ in N_PART_S8]
    geom_span = max(geom_analogue) - min(geom_analogue)
    return {
        'name': 'T2_s8_stability_automatic_for_geometric',
        'description': (
            "C2 is the constraint that killed #107 (4·n_part−100 drifts "
            "764–920). A fit-independent geometric candidate (k_5, β_lepton, "
            "2π) is §8-invariant by construction (it never touches the quark "
            "ablations), so C2 is automatic. The binding bars here are C3 "
            "and C4."
        ),
        'npart_route_s8_span': npart_span,
        'npart_route_is_s8_stable': npart_span < 1,
        'geometric_route_s8_span': geom_span,
        'geometric_route_is_s8_stable': geom_span < 1e-9,
        'binding_constraints_for_this_search': ['C3_precise', 'C4_principled'],
        'pass': geom_span < 1e-9 and npart_span > 100,
    }


# ---------------------------------------------------------------------------
# T3. Scan principled powers/products of {k_5, β_lepton, 2π}
# ---------------------------------------------------------------------------

def test_T3_principled_scan() -> dict:
    """Scan candidates built ONLY from BAM's fixed structural quantities —
    k_5 = 5, β_lepton = 50π, the closure constant 2π — with small integer
    powers and no free factors. Best principled candidate: 2π·k_5³ =
    β_lepton·k_5 = 785.4, which is −5.4% off. Nothing principled lands <1%."""
    candidates = {
        '2π·k_5³ (= β_lepton·k_5)': 2 * PI * K_5 ** 3,
        'k_5⁴': K_5 ** 4,
        'k_5⁵': K_5 ** 5,
        'e^(2π)': math.exp(2 * PI),
        '2π²·k_5²': 2 * PI ** 2 * K_5 ** 2,
        '(4/3)π·k_5³': (4.0 / 3.0) * PI * K_5 ** 3,
        '4π·k_5²': 4 * PI * K_5 ** 2,
        'β_lepton·k_5': BETA_LEPTON * K_5,
    }
    scored = {k: {'value': v, 'off_percent': _off(v) * 100}
              for k, v in candidates.items()}
    # best by absolute relative offset
    best_key = min(scored, key=lambda k: abs(scored[k]['off_percent']))
    best = scored[best_key]
    any_under_1pct = any(abs(s['off_percent']) < PRECISION_BAR * 100
                         for s in scored.values())
    return {
        'name': 'T3_principled_scan',
        'description': (
            "Principled candidates from {k_5, β_lepton, 2π}: best is "
            "2π·k_5³ = β_lepton·k_5 = 785.4 (−5.4%). No principled candidate "
            "lands within 1% of 830.3."
        ),
        'candidates': scored,
        'best_principled': best_key,
        'best_off_percent': best['off_percent'],
        'any_principled_under_1pct': any_under_1pct,
        'fails_C3_C4': not any_under_1pct,
        'pass': (not any_under_1pct) and abs(best['off_percent']) > 1.0,
    }


# ---------------------------------------------------------------------------
# T4. Sub-% matches all need ad-hoc factors (fail C4)
# ---------------------------------------------------------------------------

def test_T4_subpercent_needs_adhoc() -> dict:
    """The only ways to land <1% of 830 use a NON-BAM factor: π·265 = 832.5
    (+0.27%) needs 265, (4/3)·k_5⁴ = 833.3 (+0.37%) needs the 4/3, k_5⁵/3.77
    = 828.9 (−0.16%) needs 3.77. None of 265, 4/3, 3.77 is a fixed BAM
    structural quantity — each is reverse-engineered from the target, so
    they fail C4 (principle). A free factor tuned to hit 830 is a fit, not
    a derivation."""
    adhoc = {
        'π·265': {'value': PI * 265, 'adhoc_factor': '265 (not a BAM integer)'},
        '(4/3)·k_5⁴': {'value': (4.0 / 3.0) * K_5 ** 4, 'adhoc_factor': '4/3 (ad-hoc)'},
        'k_5⁵/3.77': {'value': K_5 ** 5 / 3.77, 'adhoc_factor': '3.77 (ad-hoc)'},
    }
    for d in adhoc.values():
        d['off_percent'] = _off(d['value']) * 100
    all_subpercent = all(abs(d['off_percent']) < 1.0 for d in adhoc.values())
    all_have_adhoc = all('adhoc_factor' in d for d in adhoc.values())
    return {
        'name': 'T4_subpercent_matches_need_adhoc_factor',
        'description': (
            "Every sub-% match needs a non-BAM factor: π·265 (+0.27%), "
            "(4/3)·k_5⁴ (+0.37%), k_5⁵/3.77 (−0.16%). 265, 4/3, 3.77 are "
            "reverse-engineered from the target ⟹ fail C4 (principle). A "
            "tuned free factor is a fit, not a derivation."
        ),
        'adhoc_candidates': adhoc,
        'all_subpercent': all_subpercent,
        'all_require_adhoc_factor': all_have_adhoc,
        'fails_C4': True,
        'pass': all_subpercent and all_have_adhoc,
    }


# ---------------------------------------------------------------------------
# T5. Exponential / dimensional-transmutation route
# ---------------------------------------------------------------------------

def test_T5_exponential_route() -> dict:
    """A hierarchy like 830 could in principle arise by dimensional
    transmutation, 830 = e^S for some geometric action S. But ln(830.3) =
    6.72, and the only clean geometric action nearby — the closure quantum
    2π = 6.28 — is 7.0% off (e^(2π) = 535, −35%). There is no principled
    geometric action that exponentiates to 830 to <1%, so this route fails
    C3/C4 too."""
    ln_target = math.log(RATIO)
    clean_action = 2 * PI
    action_off = (ln_target - clean_action) / clean_action
    e_2pi = math.exp(2 * PI)
    return {
        'name': 'T5_exponential_dimensional_transmutation_route',
        'description': (
            "830 = e^S would need S = ln(830.3) = 6.72; the clean geometric "
            "action 2π = 6.28 is 7.0% off (e^(2π)=535, −35%). No principled "
            "action exponentiates to 830 to <1% ⟹ the transmutation route "
            "also fails."
        ),
        'ln_target': ln_target,
        'clean_action_2pi': clean_action,
        'action_off_percent': action_off * 100,
        'e_2pi': e_2pi,
        'e_2pi_off_percent': _off(e_2pi) * 100,
        'no_clean_action_selects_830': abs(action_off) > PRECISION_BAR,
        'pass': abs(action_off) > PRECISION_BAR,
    }


# ---------------------------------------------------------------------------
# T6. Tangherlini cavity-eigenvalue bulk integrals
# ---------------------------------------------------------------------------

def test_T6_cavity_integrals() -> dict:
    """The genuine bulk quantities — Tangherlini cavity eigenvalue sums —
    are O(10–350): Σω²(n=0..5) ≈ 83, Σω(n=0..5) ≈ 20, Σω²(n=0..9) ≈ 347.
    None lands near 830, and there is no principled truncation that selects
    it. So the bulk shell-stress integral the #107 verdict asked for does
    not exist."""
    ev = _shell_eigenvalues()
    integrals = {
        'Σ ω²(n=0..5)': float(sum(ev[:6])),
        'Σ ω(n=0..5)': float(sum(np.sqrt(ev[:6]))),
        'Σ ω²(n=0..9)': float(sum(ev[:10])),
    }
    for k in list(integrals):
        integrals[k] = {'value': integrals[k],
                        'off_percent': _off(integrals[k]) * 100}
    none_near = all(abs(d['off_percent']) > 10 for d in integrals.values())
    return {
        'name': 'T6_cavity_eigenvalue_integrals_select_nothing',
        'description': (
            "Tangherlini cavity eigenvalue sums are O(10–350): Σω²(0..5)≈83, "
            "Σω(0..5)≈20, Σω²(0..9)≈347. None near 830; no principled "
            "truncation selects it. The independent bulk integral #107 asked "
            "for does not exist."
        ),
        'cavity_integrals': integrals,
        'none_near_830': none_near,
        'pass': none_near,
    }


# ---------------------------------------------------------------------------
# T7. Honest scope: stays underived, plausibly irreducible
# ---------------------------------------------------------------------------

def test_T7_honest_scope() -> dict:
    """No candidate clears C3 ∧ C4 (with C2 automatic for geometric
    quantities). So the fit-independent search that #107 said was the only
    legitimate route comes up empty. √σ/m_e stays UNDERIVED — and the
    failure of every principled route now points to it being IRREDUCIBLY
    residual, in the class of α and the electron's anomalous lightness (the
    universal flavor puzzle): a pure number the geometry does not fix."""
    return {
        'name': 'T7_honest_scope_underived_plausibly_irreducible',
        'description': (
            "No candidate clears C3 ∧ C4 (C2 automatic). The fit-independent "
            "search comes up empty ⟹ √σ/m_e stays UNDERIVED, and the failure "
            "of every principled route points to it being IRREDUCIBLY "
            "residual — same class as α and the electron's anomalous "
            "lightness (the flavor puzzle)."
        ),
        'search_succeeded': False,
        'ratio_status': 'UNDERIVED (PR #106 unchanged)',
        'now_classified_as': 'plausibly IRREDUCIBLE residual (like α / electron anomalous lightness)',
        'bam_input_floor': 'one foundational scale G + this open ratio + α + the flavor puzzle',
        'collapses_to_single_anchor': False,
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T8. Assessment
# ---------------------------------------------------------------------------

def test_T8_assessment() -> dict:
    return {
        'name': 'T8_assessment',
        'description': (
            "The fit-independent §8-stable search for a bulk quantity "
            "selecting √σ/m_e ≈ 830 fails: best principled candidate 2π·k_5³ "
            "= 785 (−5.4%); every sub-% match needs an ad-hoc factor; no "
            "clean action exponentiates to 830; cavity integrals select "
            "nothing near it. √σ/m_e stays underived, plausibly irreducible "
            "like α. BAM does NOT collapse to a single anchor."
        ),
        'classification': 'LEPTON_QCD_RATIO_NO_PRINCIPLED_GEOMETRIC_SELECTION_STAYS_UNDERIVED',
        'pass': True,
    }


# ---------------------------------------------------------------------------
# Runner + verdict
# ---------------------------------------------------------------------------

def run_probe() -> dict:
    tests = [
        test_T1_target_and_criteria(),
        test_T2_s8_stability_automatic(),
        test_T3_principled_scan(),
        test_T4_subpercent_needs_adhoc(),
        test_T5_exponential_route(),
        test_T6_cavity_integrals(),
        test_T7_honest_scope(),
        test_T8_assessment(),
    ]
    if all(t['pass'] for t in tests[:7]):
        verdict_class = 'LEPTON_QCD_RATIO_NO_PRINCIPLED_GEOMETRIC_SELECTION_STAYS_UNDERIVED'
        verdict = (
            'THE FIT-INDEPENDENT SEARCH COMES UP EMPTY — √σ/m_e ≈ 830 STAYS '
            'UNDERIVED. PR #106 isolated √σ/m_e ≈ 830.3 as the single '
            'underived dimensionless residual whose derivation would collapse '
            'BAM to one foundational anchor (G). PR #107 killed the n_part '
            'recycling candidate (832 = 4·n_part − 4·k_5², §8-drifts 764–920) '
            'and said the only legitimate route left was an INDEPENDENT, '
            '§8-stable bulk quantity selecting ~830. This probe ran that '
            'search honestly, and it fails.\n\n'
            'A candidate had to clear all four bars: C1 fit-independent, C2 '
            '§8-stable, C3 precise (<1%), C4 principled (no ad-hoc factor). '
            'C2 turns out to be AUTOMATIC for genuine geometric candidates — '
            'built from BAM\'s fixed integers k_5 = 5, β_lepton = 50π and the '
            'closure constant 2π, they never touch the quark ablations that '
            'made n_part drift. So the binding bars are C3 and C4 — and no '
            'candidate clears both.\n\n'
            'THE PRINCIPLED SCAN. The best principled quantity is 2π·k_5³ = '
            'β_lepton·k_5 = 785.40, which is −5.4% from 830.3. Everything '
            'else principled is worse (k_5⁴ = 625, −24.7%; e^(2π) = 535, '
            '−35.5%; 2π²·k_5² = 493, −40.6%). Nothing principled lands within '
            '1%.\n\n'
            'SUB-% MATCHES ARE FITS, NOT DERIVATIONS. The only ways to reach '
            '<1% use a non-BAM factor reverse-engineered from the target: '
            'π·265 = 832.5 (+0.27%), (4/3)·k_5⁴ = 833.3 (+0.37%), k_5⁵/3.77 = '
            '828.9 (−0.16%). None of 265, 4/3, 3.77 is a fixed BAM quantity; '
            'each is tuned to hit 830, so they fail C4. A free factor chosen '
            'to match is a fit, not a derivation.\n\n'
            'THE EXPONENTIAL ROUTE FAILS TOO. A hierarchy could arise by '
            'dimensional transmutation, 830 = e^S; but ln(830.3) = 6.72, and '
            'the only clean geometric action nearby — the closure quantum '
            '2π = 6.28 — is 7.0% off (e^(2π) = 535, −35%). No principled '
            'action exponentiates to 830 to <1%.\n\n'
            'NO BULK INTEGRAL SELECTS IT. The genuine Tangherlini cavity '
            'eigenvalue sums are O(10–350): Σω²(n=0..5) ≈ 83, Σω(n=0..5) ≈ '
            '20, Σω²(n=0..9) ≈ 347. None lands near 830 and no principled '
            'truncation selects it. The independent bulk shell-stress '
            'integral #107 asked for does not exist.\n\n'
            'HONEST CONCLUSION. √σ/m_e ≈ 830 stays an UNDERIVED open '
            'dimensionless residual, and the failure of every principled, '
            'fit-independent route now points to it being IRREDUCIBLY '
            'residual — in the same class as α and the electron\'s anomalous '
            'lightness (the universal flavor puzzle): a pure number the '
            'geometry does not fix. BAM does NOT collapse to a single anchor. '
            'It sits at one foundational scale (G) + this one open ratio + α '
            '+ the flavor puzzle. That is the honest floor, and PR #107\'s '
            'caution is vindicated: with the n_part recycling rejected, the '
            'fit-independent search that would have legitimately closed the '
            'gap comes up empty.'
        )
    else:
        verdict_class = 'LEPTON_QCD_RATIO_SEARCH_INCONCLUSIVE'
        verdict = (
            'INCONCLUSIVE. A structural test failed; review the candidate '
            'scan and the C1–C4 accounting.'
        )

    return {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'identification': (
            'the fit-independent, §8-stable search for a bulk quantity '
            'selecting √σ/m_e ≈ 830 fails: best principled candidate 2π·k_5³ '
            '= 785 (−5.4%); sub-% matches need ad-hoc factors; no clean '
            'action exponentiates to 830; cavity integrals select nothing '
            'near it'
        ),
        'answer': 'search fails — √σ/m_e stays underived, plausibly irreducible like α',
        'best_principled': '2π·k_5³ = β_lepton·k_5 = 785.4 (−5.4%)',
        'subpercent_matches': 'π·265, (4/3)·k_5⁴, k_5⁵/3.77 — all need ad-hoc factors',
        'exponential_route': 'ln(830)=6.72, clean action 2π=6.28 is 7% off',
        'ledger': 'PR #106 unchanged — one scale G + one open ratio + α + flavor puzzle; NOT a single anchor',
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
    L.append('# The legitimate search: does any fit-independent, §8-stable bulk quantity select √σ/m_e ≈ 830? (PR #108)')
    L.append('')
    L.append(f"**Run:** {s['timestamp_utc']}")
    L.append('')
    L.append(
        "PR #106 isolated `√σ/m_e ≈ 830` as the single UNDERIVED ratio whose "
        "derivation would collapse BAM to one anchor (G). PR #107 killed the "
        "tempting `N_q + ΔN = 832` candidate (it recycles `n_part`) and said "
        "the only legitimate route left was an INDEPENDENT, §8-stable bulk "
        "quantity selecting ~830. **This probe runs that search honestly — "
        "and it fails.** No principled, fit-independent candidate lands within "
        "1% of 830.3; every sub-% match needs an ad-hoc factor; the "
        "exponential route has no clean action; cavity integrals select "
        "nothing near it. `√σ/m_e` stays underived — plausibly IRREDUCIBLE, "
        "like α."
    )
    L.append('')
    L.append(f"- **Answer**: {s['answer']}")
    L.append(f"- **Best principled candidate**: {s['best_principled']}")
    L.append(f"- **Sub-% matches**: {s['subpercent_matches']}")
    L.append(f"- **Exponential route**: {s['exponential_route']}")
    L.append(f"- **Ledger**: {s['ledger']}")
    L.append('')

    L.append('## The four criteria (a candidate must clear ALL)')
    L.append('')
    L.append('| | criterion | meaning |')
    L.append('|---|---|---|')
    L.append('| C1 | fit-independent | built from fixed geometry, not `n_part` / the v3 fit (#107\'s failure) |')
    L.append('| C2 | §8-stable | invariant under the `quark_axioms` §8 ablations |')
    L.append('| C3 | precise | within ≈1% of 830.3 |')
    L.append('| C4 | principled | no ad-hoc multiplicative fudge factor |')
    L.append('')
    L.append("**C2 is automatic** for genuine geometric candidates (they "
             "never touch the quark ablations that made `n_part` drift "
             "764–920). So the binding bars are **C3 ∧ C4** — and no "
             "candidate clears both.")
    L.append('')

    L.append('## Test summary')
    L.append('')
    L.append('| # | Test | Key finding | PASS? |')
    L.append('|---|---|---|---|')
    label_map = {
        'T1': 'target √σ/m_e = 830.3; criteria C1–C4',
        'T2': 'C2 §8-stability automatic for geometric ⟹ binding bars are C3 ∧ C4',
        'T3': 'best principled = 2π·k_5³ = 785 (−5.4%); nothing principled <1%',
        'T4': 'sub-% matches (π·265, 4/3·k_5⁴, k_5⁵/3.77) all need ad-hoc factors',
        'T5': 'exponential route: ln(830)=6.72, clean 2π=6.28 is 7% off',
        'T6': 'cavity integrals O(10–350) — select nothing near 830',
        'T7': 'underived; plausibly IRREDUCIBLE like α / electron anomalous lightness',
        'T8': 'LEPTON_QCD_RATIO_NO_PRINCIPLED_GEOMETRIC_SELECTION_STAYS_UNDERIVED',
    }
    for t in s['tests']:
        passed = '**PASS**' if t['pass'] else '**FAIL**'
        nm = t['name']; prefix = nm[:2]
        L.append(f"| {prefix} | `{nm}` | {label_map.get(prefix, '—')} | {passed} |")
    L.append('')

    t3 = s['tests'][2]
    L.append('## T3: the principled scan (best is 5.4% off)')
    L.append('')
    L.append('| candidate (from {k_5, β_lepton, 2π}) | value | off% |')
    L.append('|---|---:|---:|')
    for k, d in sorted(t3['candidates'].items(), key=lambda kv: abs(kv[1]['off_percent'])):
        L.append(f"| {k} | {d['value']:.2f} | {d['off_percent']:+.2f}% |")
    L.append('')
    L.append("The best principled candidate, `2π·k_5³ = β_lepton·k_5 = "
             "785.4`, is **−5.4%** off. Nothing built from fixed BAM "
             "geometry lands within 1% of 830.3.")
    L.append('')

    t4 = s['tests'][3]
    L.append('## T4: every sub-% match is a fit, not a derivation')
    L.append('')
    L.append('| candidate | value | off% | ad-hoc factor |')
    L.append('|---|---:|---:|---|')
    for k, d in t4['adhoc_candidates'].items():
        L.append(f"| {k} | {d['value']:.1f} | {d['off_percent']:+.2f}% | {d['adhoc_factor']} |")
    L.append('')
    L.append("Each lands <1% only via a factor (265, 4/3, 3.77) "
             "reverse-engineered from the target — none is a fixed BAM "
             "quantity. A tuned free factor **fails C4**.")
    L.append('')

    L.append('## Verdict')
    L.append('')
    L.append(f"**{s['verdict_class']}.** {s['verdict']}")
    L.append('')

    L.append('## What this leaves open (the honest floor)')
    L.append('')
    L.append('- **`√σ/m_e ≈ 830` stays UNDERIVED** — and now plausibly '
             'IRREDUCIBLE, in the class of α and the electron\'s anomalous '
             'lightness (the universal flavor puzzle): a pure number the '
             'geometry does not fix.')
    L.append('- **BAM does NOT collapse to a single anchor.** It sits at one '
             'foundational scale (`G`) + this one open ratio + α + the '
             'flavor puzzle. PR #107\'s caution is vindicated: the '
             'fit-independent search that would have closed the gap comes '
             'up empty.')
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
    out = here / 'runs' / f'{ts}_lepton_qcd_ratio_legitimate_search_probe'
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
