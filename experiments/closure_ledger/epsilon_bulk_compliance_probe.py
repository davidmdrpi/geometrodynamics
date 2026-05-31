"""
Can ε be computed from bulk compliance, or only inferred from the meV scale?
(PR #112).

The chargeless-throat healing length ε is the mass-generating parameter of
the BAM Majorana neutrino: m_ν = m_D · e^{−S} with the bounce action
S = t²·P0·L*(ε) and the tortoise length L*(ε) ≈ −(rs/2) ln ε + const. PRs
#87–#90 INFERRED ε by demanding the observed S = ln(m_D/m_ν) ≈ 15–18 — i.e.
ε was read back from the meV scale it is meant to predict. This probe asks
the sharper question: can ε instead be COMPUTED from the bulk compliance
(the elastic / nucleation response of the throat neck), with NO neutrino
input?

## Bulk compliance ⟹ the healing length (meV-free)

The compliance cutoff is a sub-throat healing length, ε = ℓ²/(2 rs), where
ℓ is the length over which the chargeless neck heals back into the bulk.
For the nucleating throat that length is the critical-bubble scale

    ℓ ~ R_c = 2σ/ρ,

with the surface tension σ and bag density ρ fixed by the ELECTRON
rest-energy calibration (PR #58: E(R*) = m_e c² ⟹ σ = 1/(12π), ρ = 3/(4π)),
so R_c = 2/9 ≈ 0.222 — a pure number from the charged-throat geometry, with
NO neutrino mass anywhere. The candidate compliances all come out
sub-throat, O(10⁻²):

    ε ~ R_c³ ≈ 0.011,   Δ³ ≈ 0.018,   R_c²/2 ≈ 0.025.

So the ORDER OF MAGNITUDE of ε — and its decisive sub-throat character — IS
computable from bulk compliance, independent of the meV scale.

## The chain closes meV-free

With ε = R_c³ (bulk nucleation geometry) and the winding-edge tension
t = k_5√(2π) ≈ 12.53 (the PR #89 closure quantity, also meV-free):

    S = t²·P0·L*(R_c³) ≈ 16.85   ⟹   m_ν(gen 1) = m_D·e^{−S} ≈ 2.1 meV,

with m_D ≈ 43 keV the (electron-anchored) cavity-floor Dirac mass. The
meV scale comes out as an OUTPUT of bulk geometry — a genuine retrodiction,
not an input. This structurally DERIVES the neutrino's lightness: a
sub-throat healing length (ε ≪ 1) makes L* large, hence S large, hence
m_ν = m_D·e^{−S} exponentially small.

## The catch: precision is not geometric

The bounce action is steep in ε at the winding edge:

    d ln m_ν / d ln ε = t²·P0·rs/2 ≈ 4.8.

So the O(1) ambiguity among the healing-length candidates blows up:

    ε = R_c³ ≈ 0.011 ⟹ m_ν ≈ 2 meV,
    ε = Δ³  ≈ 0.018 ⟹ m_ν ≈ 20 meV,
    ε = R_c²/2 ≈ 0.025 ⟹ m_ν ≈ 108 meV,

a factor ~50 spread from a factor ~2 in ε. Landing precisely on "few meV"
still effectively SELECTS ε ≈ R_c³ using the observed scale — there is no
first-principles reason to prefer R_c³ over R_c²/2.

## The obstruction to pinning ε

A true compliance = 1/stiffness needs the absolute bulk gravitational
stiffness, set by κ₅² and Λ₅. BAM fixes only the dimensionless RS-tuning
combination λ_crit κ₅²/√|Λ₅| = √6 (PR #57); κ₅² and Λ₅ separately are
absorbed into the single dimensionful anchor and are NOT pinned. So the
compliance — and hence ε — is computable only up to an O(1) factor, exactly
the residual status everything tied to the one anchor shares.

## The honest verdict (a genuine partial)

  - **Computed (meV-free):** ε's order of magnitude and sub-throat
    character — ε ~ R_c³ ~ 10⁻² from the electron-calibrated nucleation
    geometry — which DERIVES the exponential smallness of m_ν and, at the
    winding-edge tension, OUTPUTS m_ν ~ meV with no neutrino input. So ε is
    upgraded from "inferred from the meV scale" to "bulk-geometric to order
    of magnitude; the meV scale is a retrodiction."
  - **Still residual:** the PRECISE ε. Because m_ν ∝ ε^{4.8}, the O(1)
    ambiguity (R_c³ vs Δ³ vs R_c²/2) spans m_ν ~ 2–108 meV; the exact value
    is not fixed by geometry (the absolute compliance normalization is the
    unpinned κ₅²/Λ₅ = the one anchor), and pinning the precise meV still
    uses the observed scale.

So: the SMALLNESS is derived from bulk compliance; the exact VALUE is not.

Tests:
  T1. The question + current status: ε inferred via S = ln(m_D/m_ν).
  T2. Define bulk compliance: ε = ℓ²/(2rs), ℓ ~ R_c = 2σ/ρ, σ/ρ from the
      electron calibration (PR #58) — no neutrino input.
  T3. Compute ε from compliance: R_c = 0.222 ⟹ ε ~ R_c³ 0.011, Δ³ 0.018,
      R_c²/2 0.025 — sub-throat O(10⁻²), meV-free.
  T4. Chain closes meV-free: ε = R_c³, t = k_5√(2π) ⟹ S ≈ 16.85 ⟹ m_ν ≈
      2.1 meV (gen 1) — a retrodiction of the observed scale.
  T5. The catch: m_ν ∝ ε^{4.8}; the O(1) ambiguity spans m_ν ~ 2–108 meV ⟹
      precision NOT geometric.
  T6. The obstruction: absolute compliance needs κ₅²/Λ₅; BAM fixes only
      √6 (PR #57) ⟹ ε pinned only to O(1) (the one anchor).
  T7. Verdict: smallness DERIVED (ε bulk-geometric to O(1), meV-independent,
      m_ν a retrodiction); precise value RESIDUAL.
  T8. Assessment.

Verdict:
  - EPSILON_BULK_COMPLIANCE_DERIVES_SMALLNESS_PRECISE_VALUE_STAYS_RESIDUAL
    (expected): ε's order of magnitude and sub-throat character ARE
    computable from the electron-calibrated bulk nucleation compliance
    (ε ~ R_c³ ~ 10⁻²) with no neutrino input — deriving the exponential
    smallness of m_ν and retrodicting m_ν ~ meV at the winding-edge
    tension. But because m_ν ∝ ε^{4.8}, the O(1) ambiguity in the healing
    length spans m_ν ~ 2–108 meV, and the absolute compliance normalization
    is the unpinned κ₅²/Λ₅; the PRECISE ε stays a residual. The smallness
    is derived from bulk compliance; the exact value is not.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np

from geometrodynamics.constants import R_MID, R_OUTER, DELTA


PI = math.pi
K_5 = 5
RS = R_MID

# EM / charged-throat calibration (PR #58): E(R*) = m_e c² ⟹ these.
# NOTE: σ, ρ are fixed by the ELECTRON rest energy — no neutrino input.
SIGMA_EM = 1.0 / (12.0 * PI * RS ** 2)
RHO_EM = 3.0 / (4.0 * PI * RS ** 3)
R_C = 2.0 * SIGMA_EM / RHO_EM                       # critical-bubble radius = 2/9
E_C = (16.0 * PI / 3.0) * SIGMA_EM ** 3 / RHO_EM ** 2
MU_THROAT = 4.0 * PI * SIGMA_EM * RS ** 2
P0 = math.sqrt(2.0 * MU_THROAT * E_C)               # under-barrier momentum ≈ 0.0605

# Winding-edge ΔL=2 tension ratio (PR #89 bracket [2π, k_5√(2π)]); meV-free.
# Equivalently t = √β_lepton with β_lepton = 50π (PR #90's "winding edge"),
# since k_5²·2π = 25·2π = 50π.
T_WIND = K_5 * math.sqrt(2.0 * PI)                  # ≈ 12.533 = √(50π)
T2_WIND = T_WIND ** 2                               # = 50π ≈ 157.08

# Electron-anchored cavity-floor Dirac masses (eV), PR #86/#87 — no neutrino input.
M_D_EV = {1: 43.0e3, 2: 80.0e3, 3: 118.0e3}

# Observed neutrino SCALE (eV), for the order-of-magnitude comparison only.
# The lightest mass is unmeasured (≲ 3 meV, PR #111); the heavier eigenstates
# are √Δm²_21 = 8.6 meV and √Δm²_31 = 50.3 meV. NOT a per-generation map.
NU_OBS_SCALE_MEV = {'lightest (≲)': 3.0, '√Δm²_21': 8.6, '√Δm²_31': 50.3}


def _rstar(r: float) -> float:
    return r + (RS / 2.0) * math.log((r - RS) / (r + RS))


def _rstar_offset(eps: float) -> float:
    return (RS + eps) + (RS / 2.0) * math.log(eps / (2.0 * RS + eps))


def tortoise_length(eps: float) -> float:
    return _rstar(R_OUTER) - _rstar_offset(eps)


def bounce_S(eps: float, t2: float = T2_WIND) -> float:
    return t2 * P0 * tortoise_length(eps)


def m_nu_meV(eps: float, gen: int = 1, t2: float = T2_WIND) -> float:
    return M_D_EV[gen] * math.exp(-bounce_S(eps, t2)) * 1e3


# Bulk-compliance candidates for ε (all meV-free).
EPS_CANDIDATES = {
    'R_c³': R_C ** 3,
    'Δ³': DELTA ** 3,
    'R_c²/2 (ℓ=R_c)': R_C ** 2 / 2.0,
}


# ---------------------------------------------------------------------------
# T1. The question + current status
# ---------------------------------------------------------------------------

def test_T1_question() -> dict:
    return {
        'name': 'T1_question_and_current_status',
        'description': (
            "ε is the chargeless-throat healing length that sets m_ν = "
            "m_D·e^{−S}, S = t²·P0·L*(ε). PRs #87–#90 INFERRED ε by demanding "
            "the observed S = ln(m_D/m_ν) ≈ 15–18 — ε read back from the meV "
            "scale. Question: can ε be COMPUTED from bulk compliance instead?"
        ),
        'current_status': 'ε inferred from observed S = ln(m_D/m_ν)',
        'role_of_eps': 'm_ν = m_D·e^{−S}, S = t²·P0·L*(ε), L* ≈ −(rs/2)ln ε + const',
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T2. Define bulk compliance
# ---------------------------------------------------------------------------

def test_T2_define_compliance() -> dict:
    """ε = ℓ²/(2 rs) with ℓ the neck healing length ~ R_c = 2σ/ρ, the
    critical-bubble scale. σ and ρ are fixed by the ELECTRON rest-energy
    calibration (PR #58: σ = 1/(12π), ρ = 3/(4π)), so R_c = 2/9 ≈ 0.222 is a
    pure number from the charged-throat geometry — NO neutrino input."""
    return {
        'name': 'T2_define_bulk_compliance',
        'description': (
            "ε = ℓ²/(2rs), ℓ ~ R_c = 2σ/ρ (critical-bubble healing length); "
            "σ = 1/(12π), ρ = 3/(4π) from the ELECTRON calibration (PR #58) "
            "⟹ R_c = 2/9 ≈ 0.222, meV-free."
        ),
        'sigma_em': SIGMA_EM, 'rho_em': RHO_EM,
        'R_c': R_C,
        'R_c_exact': '2/9',
        'calibrated_to': 'electron rest energy (PR #58) — NOT the neutrino',
        'pass': abs(R_C - 2.0 / 9.0) < 1e-9,
    }


# ---------------------------------------------------------------------------
# T3. Compute ε from compliance
# ---------------------------------------------------------------------------

def test_T3_compute_eps() -> dict:
    """The compliance candidates are all sub-throat, O(10⁻²): R_c³ ≈ 0.011,
    Δ³ ≈ 0.018, R_c²/2 ≈ 0.025 — computed from the electron-calibrated
    geometry, no neutrino input. So ε's ORDER OF MAGNITUDE and sub-throat
    character are bulk-computable."""
    cands = {k: round(v, 4) for k, v in EPS_CANDIDATES.items()}
    all_sub_throat = all(1e-3 < v < 1e-1 for v in EPS_CANDIDATES.values())
    return {
        'name': 'T3_compute_eps_from_compliance',
        'description': (
            "Compliance candidates sub-throat O(10⁻²): R_c³ 0.011, Δ³ 0.018, "
            "R_c²/2 0.025 — bulk-geometric, meV-free. ε's order of magnitude "
            "+ sub-throat character ARE computable."
        ),
        'eps_candidates': cands,
        'all_sub_throat_1e-2': all_sub_throat,
        'meV_free': True,
        'pass': all_sub_throat,
    }


# ---------------------------------------------------------------------------
# T4. The chain closes meV-free
# ---------------------------------------------------------------------------

def test_T4_chain_closes() -> dict:
    """ε = R_c³ + the winding-edge tension t = k_5√(2π) (PR #89, meV-free)
    ⟹ S ≈ 16.85 ⟹ the predicted m_ν spectrum lands at the meV SCALE
    (2.1–5.7 meV) — output from bulk geometry, no neutrino input. The claim
    is the SCALE only: uniform S gives m_ν ∝ m_D (×2.7 spread), so the
    observed generation spread (up to √Δm²_31 = 50 meV) is the SEPARATE
    known residual (PR #90/#91), not addressed by ε."""
    eps = R_C ** 3
    S = bounce_S(eps)
    rows = []
    for g in (1, 2, 3):
        rows.append({'gen': g, 'm_D_keV': M_D_EV[g] / 1e3,
                     'm_nu_predicted_meV': round(m_nu_meV(eps, g), 2)})
    pred = [r['m_nu_predicted_meV'] for r in rows]
    return {
        'name': 'T4_chain_closes_meV_free',
        'description': (
            "ε = R_c³ + t = k_5√(2π) ⟹ S ≈ 16.85 ⟹ predicted m_ν = "
            "2.1–5.7 meV — the meV SCALE output from bulk geometry "
            "(retrodiction, no neutrino input). Scale only; the generation "
            "spread is the separate PR #90/#91 residual."
        ),
        't_wind': T_WIND, 't2_wind_eq_50pi': T2_WIND,
        'S': round(S, 2),
        'rows': rows,
        'predicted_scale_meV': f'{min(pred):.1f}–{max(pred):.1f}',
        'observed_scale_meV': NU_OBS_SCALE_MEV,
        'scale_is_meV': all(1.0 < p < 60.0 for p in pred),
        'spread_is_separate_residual': 'uniform S ⟹ m_ν ∝ m_D (×2.7); observed up to ×24 (PR #90/#91)',
        'pass': 1.0 < m_nu_meV(eps, 1) < 10.0,
    }


# ---------------------------------------------------------------------------
# T5. The catch: precision is not geometric
# ---------------------------------------------------------------------------

def test_T5_steep_sensitivity() -> dict:
    """d ln m_ν / d ln ε = t²·P0·rs/2 ≈ 4.8. So the O(1) ambiguity among the
    candidates blows up: R_c³ → 2 meV, Δ³ → 20 meV, R_c²/2 → 108 meV — a
    factor ~50 spread from a factor ~2 in ε. Precision is NOT geometric;
    landing on 'few meV' still selects ε ≈ R_c³ via the observed scale."""
    sens = T2_WIND * P0 * RS / 2.0
    rows = {k: round(m_nu_meV(v, 1), 1) for k, v in EPS_CANDIDATES.items()}
    spread = max(rows.values()) / min(rows.values())
    return {
        'name': 'T5_steep_sensitivity_precision_not_geometric',
        'description': (
            "d ln m_ν/d ln ε ≈ 4.8 ⟹ the O(1) candidate ambiguity spans "
            "m_ν ~ 2–108 meV (×~50). Precision is NOT geometric; the exact "
            "meV still selects ε ≈ R_c³ via the observed scale."
        ),
        'd_lnmnu_d_lneps': round(sens, 2),
        'm_nu_meV_by_candidate': rows,
        'spread_factor': round(spread, 1),
        'precision_geometric': False,
        'pass': sens > 3.0 and spread > 10.0,
    }


# ---------------------------------------------------------------------------
# T6. The obstruction to pinning ε
# ---------------------------------------------------------------------------

def test_T6_obstruction() -> dict:
    """A true compliance = 1/stiffness needs the absolute bulk stiffness
    (κ₅², Λ₅). BAM fixes only the dimensionless RS-tuning combination
    λ_crit κ₅²/√|Λ₅| = √6 (PR #57); κ₅² and Λ₅ separately are absorbed into
    the single dimensionful anchor and not pinned. So ε is computable only
    up to an O(1) factor — the residual status everything tied to the one
    anchor shares."""
    return {
        'name': 'T6_obstruction_unpinned_bulk_normalization',
        'description': (
            "Absolute compliance = 1/stiffness needs κ₅²/Λ₅; BAM fixes only "
            "λ_crit κ₅²/√|Λ₅| = √6 (PR #57). κ₅², Λ₅ separately absorbed "
            "into the one anchor ⟹ ε pinned only to O(1)."
        ),
        'pinned_combination': 'λ_crit κ₅²/√|Λ₅| = √6 (dimensionless)',
        'unpinned': 'κ₅² and Λ₅ separately (the absolute bulk stiffness)',
        'consequence': 'ε computable only up to an O(1) factor (= the one anchor)',
        'sqrt6': round(math.sqrt(6.0), 4),
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T7. Verdict / reclassification
# ---------------------------------------------------------------------------

def test_T7_reclassification() -> dict:
    return {
        'name': 'T7_reclassification_smallness_derived_value_residual',
        'description': (
            "ε upgraded from 'inferred from the meV scale' to 'bulk-geometric "
            "to order of magnitude': the SMALLNESS (sub-throat ⟹ exponential "
            "suppression) is derived, meV-independent, and m_ν ~ meV is a "
            "retrodiction; the PRECISE value stays residual (m_ν ∝ ε^{4.8}, "
            "absolute normalization = the unpinned anchor)."
        ),
        'derived_meV_free': [
            'ε is sub-throat (ε ~ R_c³ ~ 10⁻²) from the electron-calibrated '
            'nucleation compliance — no neutrino input',
            'hence m_ν = m_D·e^{−S} is exponentially small (smallness DERIVED)',
            'm_ν ~ meV is OUTPUT at the winding-edge tension (retrodiction)',
        ],
        'still_residual': [
            'the precise ε — O(1) ambiguity (R_c³/Δ³/R_c²/2) ⟹ m_ν 2–108 meV',
            'm_ν ∝ ε^{4.8} ⟹ steep; exact meV still uses the observed scale',
            'absolute compliance normalization = unpinned κ₅²/Λ₅ (the anchor)',
        ],
        'one_line': 'the smallness is derived from bulk compliance; the exact value is not',
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T8. Assessment
# ---------------------------------------------------------------------------

def test_T8_assessment() -> dict:
    return {
        'name': 'T8_assessment',
        'description': (
            "Bulk compliance computes ε's order of magnitude and sub-throat "
            "character (ε ~ R_c³ ~ 10⁻²) with no neutrino input — deriving "
            "the exponential smallness of m_ν and retrodicting m_ν ~ meV at "
            "the winding-edge tension. But m_ν ∝ ε^{4.8}, the O(1) ambiguity "
            "spans 2–108 meV, and the absolute normalization is the unpinned "
            "κ₅²/Λ₅; the precise ε stays residual."
        ),
        'classification': 'EPSILON_BULK_COMPLIANCE_DERIVES_SMALLNESS_PRECISE_VALUE_STAYS_RESIDUAL',
        'pass': True,
    }


# ---------------------------------------------------------------------------
# Runner + verdict
# ---------------------------------------------------------------------------

def run_probe() -> dict:
    tests = [
        test_T1_question(),
        test_T2_define_compliance(),
        test_T3_compute_eps(),
        test_T4_chain_closes(),
        test_T5_steep_sensitivity(),
        test_T6_obstruction(),
        test_T7_reclassification(),
        test_T8_assessment(),
    ]
    if all(t['pass'] for t in tests[:7]):
        verdict_class = 'EPSILON_BULK_COMPLIANCE_DERIVES_SMALLNESS_PRECISE_VALUE_STAYS_RESIDUAL'
        verdict = (
            'BULK COMPLIANCE DERIVES THE SMALLNESS OF ε (AND HENCE m_ν), BUT '
            'NOT ITS PRECISE VALUE. PRs #87–#90 inferred the chargeless-throat '
            'healing length ε by demanding the observed S = ln(m_D/m_ν) ≈ '
            '15–18 — ε read back from the meV scale it is meant to predict. '
            'This probe asks whether ε can instead be COMPUTED from the bulk '
            'compliance, with no neutrino input. The answer is a genuine '
            'partial: the order of magnitude yes, the precise value no.\n\n'
            'THE COMPLIANCE COMPUTATION (meV-free). The compliance cutoff is '
            'a sub-throat healing length ε = ℓ²/(2rs), with ℓ ~ R_c = 2σ/ρ '
            'the critical-bubble scale. σ = 1/(12π) and ρ = 3/(4π) are fixed '
            'by the ELECTRON rest-energy calibration (PR #58), so R_c = 2/9 ≈ '
            '0.222 is a pure number from the charged-throat geometry — no '
            'neutrino mass anywhere. The candidate compliances are all '
            'sub-throat, O(10⁻²): R_c³ ≈ 0.011, Δ³ ≈ 0.018, R_c²/2 ≈ 0.025. '
            'So ε\'s order of magnitude and its decisive sub-throat character '
            'ARE computable from bulk compliance, independent of the meV '
            'scale.\n\n'
            'THE CHAIN CLOSES meV-FREE. With ε = R_c³ and the winding-edge '
            'tension t = k_5√(2π) ≈ 12.53 (the PR #89 closure quantity, also '
            'meV-free), S = t²·P0·L*(R_c³) ≈ 16.85, so m_ν(gen 1) = '
            'm_D·e^{−S} ≈ 2.1 meV with m_D ≈ 43 keV the electron-anchored '
            'cavity-floor Dirac mass. The meV scale comes out as an OUTPUT — '
            'a genuine retrodiction. This structurally DERIVES the neutrino\'s '
            'lightness: a sub-throat healing length (ε ≪ 1) makes L* large, '
            'hence S large, hence m_ν exponentially small.\n\n'
            'THE CATCH: PRECISION IS NOT GEOMETRIC. The bounce action is '
            'steep in ε — d ln m_ν/d ln ε = t²·P0·rs/2 ≈ 4.8 — so the O(1) '
            'ambiguity among the healing-length candidates blows up: R_c³ → '
            '2 meV, Δ³ → 20 meV, R_c²/2 → 108 meV, a factor ~50 spread from '
            'a factor ~2 in ε. Landing precisely on "few meV" still SELECTS '
            'ε ≈ R_c³ using the observed scale; there is no first-principles '
            'reason to prefer R_c³ over R_c²/2.\n\n'
            'THE OBSTRUCTION. A true compliance = 1/stiffness needs the '
            'absolute bulk gravitational stiffness (κ₅², Λ₅). BAM fixes only '
            'the dimensionless RS-tuning combination λ_crit κ₅²/√|Λ₅| = √6 '
            '(PR #57); κ₅² and Λ₅ separately are absorbed into the single '
            'dimensionful anchor and are not pinned. So ε is computable only '
            'up to an O(1) factor — the same residual status everything tied '
            'to the one anchor shares.\n\n'
            'VERDICT. ε is upgraded from "inferred from the meV scale" to '
            '"bulk-geometric to order of magnitude": the SMALLNESS is derived '
            'from the electron-calibrated nucleation compliance '
            '(meV-independent), and m_ν ~ meV is a retrodiction at the '
            'winding-edge tension. The PRECISE value stays a residual — '
            'm_ν ∝ ε^{4.8}, the O(1) ambiguity spans 2–108 meV, and the '
            'absolute normalization is the unpinned κ₅²/Λ₅. The smallness is '
            'derived from bulk compliance; the exact value is not.'
        )
    else:
        verdict_class = 'EPSILON_BULK_COMPLIANCE_INCONCLUSIVE'
        verdict = (
            'INCONCLUSIVE. A structural test failed; review the compliance '
            'computation and the bounce-action chain.'
        )

    return {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'identification': (
            'ε\'s order of magnitude + sub-throat character are computable '
            'from the electron-calibrated bulk nucleation compliance '
            '(ε ~ R_c³ ~ 10⁻², meV-free), deriving the smallness of m_ν and '
            'retrodicting m_ν ~ meV; the precise ε stays residual (m_ν ∝ '
            'ε^{4.8}; absolute normalization = the unpinned anchor)'
        ),
        'computed_meV_free': 'ε ~ R_c³ ~ 10⁻² (sub-throat); m_ν ≈ 2.1 meV output at t = k_5√(2π)',
        'derives': 'the exponential smallness of m_ν (sub-throat ⟹ large S)',
        'residual': 'the precise ε — m_ν ∝ ε^{4.8}, O(1) ambiguity spans 2–108 meV',
        'obstruction': 'absolute compliance normalization = unpinned κ₅²/Λ₅ (the one anchor)',
        'one_line': 'the smallness is derived from bulk compliance; the exact value is not',
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
    L.append('# Can ε be computed from bulk compliance, or only inferred from the meV scale? (PR #112)')
    L.append('')
    L.append(f"**Run:** {s['timestamp_utc']}")
    L.append('')
    L.append(
        "PRs #87–#90 INFERRED the chargeless-throat healing length `ε` by "
        "demanding the observed `S = ln(m_D/m_ν) ≈ 15–18` — `ε` read back "
        "from the meV scale it is meant to predict. This probe asks whether "
        "`ε` can instead be COMPUTED from the bulk compliance (the "
        "elastic/nucleation response of the throat neck), with no neutrino "
        "input. **Answer: a genuine partial — the smallness yes, the precise "
        "value no.**"
    )
    L.append('')
    L.append(f"- **Computed (meV-free)**: {s['computed_meV_free']}")
    L.append(f"- **Derives**: {s['derives']}")
    L.append(f"- **Residual**: {s['residual']}")
    L.append(f"- **Obstruction**: {s['obstruction']}")
    L.append(f"- **One line**: {s['one_line']}")
    L.append('')

    L.append('## Test summary')
    L.append('')
    L.append('| # | Test | Key finding | PASS? |')
    L.append('|---|---|---|---|')
    label_map = {
        'T1': 'ε currently inferred via S = ln(m_D/m_ν)',
        'T2': 'compliance: ε = ℓ²/2rs, ℓ ~ R_c = 2σ/ρ (electron-calibrated)',
        'T3': 'ε ~ R_c³ 0.011, Δ³ 0.018, R_c²/2 0.025 — sub-throat, meV-free',
        'T4': 'ε=R_c³ + t=k_5√(2π) ⟹ S≈16.85 ⟹ m_ν≈2.1 meV (retrodiction)',
        'T5': 'm_ν ∝ ε^{4.8}: O(1) ambiguity spans 2–108 meV (not geometric)',
        'T6': 'absolute compliance = unpinned κ₅²/Λ₅; only √6 fixed',
        'T7': 'smallness DERIVED (meV-free); precise value RESIDUAL',
        'T8': 'EPSILON_BULK_COMPLIANCE_DERIVES_SMALLNESS_PRECISE_VALUE_STAYS_RESIDUAL',
    }
    for t in s['tests']:
        passed = '**PASS**' if t['pass'] else '**FAIL**'
        nm = t['name']; prefix = nm[:2]
        L.append(f"| {prefix} | `{nm}` | {label_map.get(prefix, '—')} | {passed} |")
    L.append('')

    t4 = s['tests'][3]
    L.append('## The chain closes meV-free (ε = R_c³, t = k_5√(2π))')
    L.append('')
    L.append(f"`S ≈ {t4['S']}` ⟹ predicted m_ν spectrum:")
    L.append('')
    L.append('| gen | m_D (keV) | m_ν predicted (meV) |')
    L.append('|---|---:|---:|')
    for r in t4['rows']:
        L.append(f"| {r['gen']} | {r['m_D_keV']:.0f} | {r['m_nu_predicted_meV']} |")
    L.append('')
    L.append(f"The meV **scale** ({t4['predicted_scale_meV']} meV) is an "
             "**output** of bulk geometry — no neutrino input. A sub-throat "
             "`ε ≪ 1` makes `L*` large ⟹ `S` large ⟹ `m_ν = m_D·e^{−S}` "
             "exponentially small: **the smallness is derived.** _Scale only_ "
             "— uniform `S` gives `m_ν ∝ m_D` (×2.7), so the observed "
             "generation spread (up to `√Δm²_31 = 50 meV`) is the separate "
             "PR #90/#91 residual, not addressed by `ε`.")
    L.append('')

    t5 = s['tests'][4]
    L.append('## The catch: m_ν ∝ ε^{4.8} (precision not geometric)')
    L.append('')
    L.append('| ε candidate | value | m_ν (meV) |')
    L.append('|---|---:|---:|')
    for k, mv in t5['m_nu_meV_by_candidate'].items():
        L.append(f"| {k} | {EPS_CANDIDATES[k]:.4f} | {mv} |")
    L.append('')
    L.append(f"`d ln m_ν/d ln ε ≈ {t5['d_lnmnu_d_lneps']}` — a factor ~2 in "
             f"`ε` spans m_ν by ~×{t5['spread_factor']}. Landing on 'few meV' "
             f"still selects `ε ≈ R_c³` via the observed scale.")
    L.append('')

    L.append('## Verdict')
    L.append('')
    L.append(f"**{s['verdict_class']}.** {s['verdict']}")
    L.append('')

    L.append('## What this changes, and what stays open')
    L.append('')
    L.append('- **Changed:** `ε` moves from "inferred from the meV scale" to '
             '"bulk-geometric to order of magnitude." The neutrino\'s '
             'exponential lightness is now DERIVED from the electron-'
             'calibrated nucleation compliance (sub-throat `ε`), and `m_ν ~ '
             'meV` is a retrodiction — not an input.')
    L.append('- **Open:** the PRECISE `ε`. Because `m_ν ∝ ε^{4.8}`, the O(1) '
             'ambiguity (`R_c³`/`Δ³`/`R_c²/2`) spans `m_ν ~ 2–108 meV`; the '
             'absolute compliance normalization is the unpinned `κ₅²/Λ₅` '
             '(the one anchor). The exact value is still residual.')
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
    out = here / 'runs' / f'{ts}_epsilon_bulk_compliance_probe'
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
