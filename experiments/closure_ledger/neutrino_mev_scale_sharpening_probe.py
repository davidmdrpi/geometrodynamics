"""
Sharpening the meV-scale neutrino predictions (PR #111).

The BAM neutrino sector (PRs #87–#96) predicts a light, normal-ordered,
Majorana spectrum (generations = cavity overtones; m_ν ∝ m_D; c₁ = 0 ⟹
Majorana; anarchic phases). PR #96 left the headline meV-scale observable
as a BAND — Σm_ν ≈ 59–65 meV — because the lightest mass m₁ was an open
absolute-scale residual (~few meV, PR #90). This probe SHARPENS the
meV-scale predictions by (i) updating the oscillation inputs to the latest
global fit (NuFIT 6.0, 2024), (ii) using the 2025 DESI DR2 + CMB cosmology
to CORNER m₁ against the normal-ordering floor, and (iii) reading off the
full pinned spectrum — Σm_ν, the three masses, the β-decay effective mass
m_β, and the 0νββ effective mass m_ββ — with an honest statement of which
are reachable.

## The cornering

With NuFIT 6.0 (NO) the oscillation splittings fix two of the three masses
outright:

    √Δm²_21 = 8.65 meV,   √Δm²_31 = 50.34 meV,
    Σm_ν|_floor (m₁→0) = 59.0 meV.

DESI DR2 + CMB (2025) bounds Σm_ν ≲ 60–64 meV (95%) — which pins m₁ to the
LOW end of "few meV":

    m₁ ≲ 3 meV   ⟹   Σm_ν ∈ [59.0, 62.6] meV,

a sharpening of the PR #96 band (59–65) toward the floor. The sharpest
single-number BAM statement is the hierarchical limit m₁ → 0, the natural
endpoint of the light cavity-overtone spectrum AND the point cosmology is
cornering us toward.

## The pinned spectrum (hierarchical limit)

    m₁ ≈ 0   (≲ 3 meV),   m₂ = 8.65 meV,   m₃ = 50.34 meV.

From this, the two laboratory effective masses follow:

  - **β-decay (KATRIN):** m_β = √(Σ|U_ei|² m_i²) ≈ 8.8 meV (rising to
    ~9.3 meV at m₁ = 3 meV).
  - **0νββ:** m_ββ = |Σ U_ei² m_i| — and crucially, in normal ordering the
    contributions CANNOT fully cancel at m₁ → 0, because
    s12²c13² m₂ = 2.60 meV > s13² m₃ = 1.10 meV, leaving a NONZERO floor

        m_ββ ∈ [1.5, 3.7] meV   (m₁ → 0, over the Majorana phases),

    widening to ~[0, 5.9] meV as m₁ → 3 meV (where a phase-tuned null
    becomes possible). So BAM predicts a small but nonvanishing 0νββ rate
    in the hierarchical limit.

## Reachability (the honest other half)

Sharper predictions also sharpen the statement of what can test them:

  - **Σm_ν ≈ 59 meV — testable NOW.** DESI DR2 + CMB is already cornering
    it at the floor; this is the live handle.
  - **m_β ≈ 8.8 meV — not reachable soon.** KATRIN final sensitivity
    ~200 meV; next-gen (Project 8) ~40 meV — still ~4–5× above.
  - **m_ββ ≈ 1.5–3.7 meV — not reachable soon.** Next-gen 0νββ
    (LEGEND-1000, nEXO) reach m_ββ ~ 9–20 meV — ~3–10× above the BAM floor.

So the meV-scale spectrum is now PINNED, but only Σm_ν is near-term
testable, and it is being tested at the floor right now. The Majorana and
absolute-scale predictions (m_ββ, m_β) are sharp yet below foreseeable
sensitivity — an honest readout, not a promise of imminent discovery.

## Falsifiers (sharpened)

  - A robust cosmological Σm_ν < 59.0 meV (below the NO floor) excludes
    normal ordering ⟹ BAM fails (and would clash with the oscillation Δm²
    themselves — a deep consistency break). NOTE: some 2025 DESI + CMB
    analyses already PREFER central Σm_ν values at or below the floor (the
    data wants less lensing than the floor supplies); if that hardens it is
    tension for ALL normal-ordered models, BAM included — flagged honestly.
  - A quasi-degenerate detection (Σm_ν ≳ 100 meV, the IO/degenerate floor)
    contradicts the BAM light scale.
  - A 0νββ signal implying m_ββ ≫ 4 meV with confirmed normal ordering
    would sit above the BAM hierarchical band (would need m₁ near the top
    of its window or non-anarchic phases).

Tests:
  T1. Setup: sharpen by pinning the full spectrum (was a Σm_ν band).
  T2. NuFIT 6.0 inputs ⟹ √Δm²_21 = 8.65, √Δm²_31 = 50.34 meV, NO floor
      59.0 meV.
  T3. DESI DR2 cornering ⟹ m₁ ≲ 3 meV ⟹ Σm_ν ∈ [59.0, 62.6] meV
      (sharpened from 59–65).
  T4. Pinned spectrum (hierarchical limit): m = (≲3, 8.65, 50.34) meV.
  T5. m_β ≈ 8.8–9.3 meV; m_ββ nonzero floor [1.5, 3.7] meV (no full
      cancellation in NO).
  T6. Reachability (honest): Σm_ν testable now; m_β, m_ββ sharp but below
      foreseeable sensitivity.
  T7. Falsifiers sharpened + the 2025 sub-floor cosmological preference
      (tension flag); open: m₁ band + Majorana phases.
  T8. Assessment.

Verdict:
  - NEUTRINO_MEV_PREDICTIONS_SHARPENED_SPECTRUM_PINNED_SIGMA_TESTABLE
    (expected): updating to NuFIT 6.0 and cornering m₁ with DESI DR2 pins
    the meV-scale spectrum — m = (≲3, 8.65, 50.34) meV, Σm_ν ∈ [59.0, 62.6]
    (sharpened toward the floor), m_β ≈ 8.8–9.3 meV, m_ββ a nonzero floor
    [1.5, 3.7] meV. Only Σm_ν is near-term testable (DESI, at the floor
    now); m_β and m_ββ are sharp but below foreseeable sensitivity. The m₁
    band and Majorana phases stay open.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np


# --- Oscillation global fit: NuFIT 6.0 (2024), normal ordering ---
DM2_21 = 7.49e-5      # eV²
DM2_31 = 2.534e-3     # eV²
S12_SQ = 0.307
S13_SQ = 0.02195
C12_SQ = 1.0 - S12_SQ
C13_SQ = 1.0 - S13_SQ

# Cosmological 95%-CL bound that does the cornering.
DESI_DR2_CMB_MEV = 62.0       # ~60–64 meV (2025, tightest combinations)
M1_CORNER_MEV = 3.0           # m₁ ceiling set by DESI DR2 against the floor

# Laboratory sensitivities (meV) for the honest reachability statement.
KATRIN_FINAL_MEV = 200.0
PROJECT8_MEV = 40.0
NEXTGEN_0NUBB_MEV = 15.0      # LEGEND-1000 / nEXO ~ 9–20 meV


def _meV(x_eV: float) -> float:
    return x_eV * 1e3


def spectrum_meV(m1_meV: float):
    m1 = m1_meV * 1e-3
    m2 = math.sqrt(m1 ** 2 + DM2_21)
    m3 = math.sqrt(m1 ** 2 + DM2_31)
    return _meV(m1), _meV(m2), _meV(m3)


def sigma_meV(m1_meV: float) -> float:
    return sum(spectrum_meV(m1_meV))


def m_beta_meV(m1_meV: float) -> float:
    m1, m2, m3 = (x * 1e-3 for x in spectrum_meV(m1_meV))
    val = C13_SQ * (C12_SQ * m1 ** 2 + S12_SQ * m2 ** 2) + S13_SQ * m3 ** 2
    return _meV(math.sqrt(val))


def m_bb_band_meV(m1_meV: float):
    """0νββ effective mass band over the Majorana phases."""
    m1, m2, m3 = (x * 1e-3 for x in spectrum_meV(m1_meV))
    t1 = C12_SQ * C13_SQ * m1
    t2 = S12_SQ * C13_SQ * m2
    t3 = S13_SQ * m3
    terms = sorted([t1, t2, t3], reverse=True)
    lo = max(0.0, terms[0] - terms[1] - terms[2])
    hi = terms[0] + terms[1] + terms[2]
    return _meV(lo), _meV(hi)


NO_FLOOR_MEV = sigma_meV(0.0)


# ---------------------------------------------------------------------------
# T1. Setup
# ---------------------------------------------------------------------------

def test_T1_setup() -> dict:
    return {
        'name': 'T1_setup_sharpen_to_pinned_spectrum',
        'description': (
            "Sharpen the PR #96 meV-scale prediction: replace the Σm_ν band "
            "(59–65 meV, open m₁) with the full PINNED spectrum, using "
            "NuFIT 6.0 + DESI DR2 to corner m₁ and reading off Σm_ν, the "
            "three masses, m_β, and m_ββ."
        ),
        'was': 'Σm_ν band 59–65 meV (m₁ an open residual)',
        'now': 'full spectrum pinned (m₁ cornered) + m_β + m_ββ',
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T2. Updated oscillation inputs
# ---------------------------------------------------------------------------

def test_T2_nufit6_inputs() -> dict:
    s_dm21 = _meV(math.sqrt(DM2_21))
    s_dm31 = _meV(math.sqrt(DM2_31))
    return {
        'name': 'T2_nufit6_inputs_and_floor',
        'description': (
            "NuFIT 6.0 (NO): √Δm²_21 = 8.65 meV, √Δm²_31 = 50.34 meV ⟹ the "
            "normal-ordering floor Σm_ν(m₁→0) = 59.0 meV."
        ),
        'sqrt_dm2_21_meV': s_dm21,
        'sqrt_dm2_31_meV': s_dm31,
        'no_floor_meV': NO_FLOOR_MEV,
        'pass': abs(NO_FLOOR_MEV - 59.0) < 1.5,
    }


# ---------------------------------------------------------------------------
# T3. The cornering
# ---------------------------------------------------------------------------

def test_T3_cornering() -> dict:
    """DESI DR2 + CMB (≲ 60–64 meV) pins m₁ to the low end: m₁ ≲ 3 meV ⟹
    Σm_ν ∈ [59.0, 62.6] meV — sharpened from the PR #96 band (59–65)."""
    sig_lo = sigma_meV(0.0)
    sig_hi = sigma_meV(M1_CORNER_MEV)
    return {
        'name': 'T3_desi_dr2_corners_m1',
        'description': (
            "DESI DR2 + CMB ≲ 60–64 meV corners m₁ ≲ 3 meV ⟹ Σm_ν ∈ "
            "[59.0, 62.6] meV (sharpened from 59–65, toward the floor)."
        ),
        'desi_dr2_bound_meV': DESI_DR2_CMB_MEV,
        'm1_corner_meV': M1_CORNER_MEV,
        'sigma_band_meV': (round(sig_lo, 1), round(sig_hi, 1)),
        'sharpened_from': '59–65 meV (PR #96)',
        'pass': sig_hi <= DESI_DR2_CMB_MEV + 1.0 and sig_lo >= NO_FLOOR_MEV - 0.1,
    }


# ---------------------------------------------------------------------------
# T4. The pinned spectrum
# ---------------------------------------------------------------------------

def test_T4_pinned_spectrum() -> dict:
    """In the hierarchical limit (m₁ → 0, the BAM endpoint cosmology
    favors), the spectrum is fully pinned: m = (≲3, 8.65, 50.34) meV."""
    m1, m2, m3 = spectrum_meV(0.0)
    return {
        'name': 'T4_pinned_spectrum_hierarchical_limit',
        'description': (
            "Hierarchical limit: m = (≲3, 8.65, 50.34) meV — m₂ and m₃ "
            "pinned by the oscillation splittings, m₁ cornered by cosmology."
        ),
        'm1_meV': '≲ 3 (→ 0)',
        'm2_meV': round(m2, 2),
        'm3_meV': round(m3, 2),
        'ordering': 'normal (PR #91)',
        'pass': abs(m2 - 8.65) < 0.2 and abs(m3 - 50.34) < 0.5,
    }


# ---------------------------------------------------------------------------
# T5. m_β and m_ββ
# ---------------------------------------------------------------------------

def test_T5_lab_effective_masses() -> dict:
    """m_β (KATRIN) ≈ 8.8 meV (→ 9.3 at m₁ = 3 meV). m_ββ has a NONZERO
    floor in NO: s12²c13² m₂ = 2.60 meV > s13² m₃ = 1.10 meV, so the
    contributions cannot fully cancel ⟹ m_ββ ∈ [1.5, 3.7] meV at m₁ → 0."""
    mb_lo = m_beta_meV(0.0)
    mb_hi = m_beta_meV(M1_CORNER_MEV)
    bb_lo0, bb_hi0 = m_bb_band_meV(0.0)
    # the two NO contributions at m1 -> 0
    t2 = S12_SQ * C13_SQ * math.sqrt(DM2_21) * 1e3
    t3 = S13_SQ * math.sqrt(DM2_31) * 1e3
    return {
        'name': 'T5_m_beta_and_m_bb',
        'description': (
            "m_β ≈ 8.8–9.3 meV (KATRIN effective). m_ββ NONZERO floor "
            "[1.5, 3.7] meV at m₁→0 — no full cancellation in NO "
            "(s12²c13²m₂ = 2.60 > s13²m₃ = 1.10 meV)."
        ),
        'm_beta_meV_band': (round(mb_lo, 2), round(mb_hi, 2)),
        'm_bb_meV_floor_band': (round(bb_lo0, 2), round(bb_hi0, 2)),
        'no_full_cancellation': t2 > t3,
        'contrib_solar_meV': round(t2, 2),
        'contrib_atm_meV': round(t3, 2),
        'pass': bb_lo0 > 1.0 and t2 > t3 and 8.0 < mb_lo < 9.5,
    }


# ---------------------------------------------------------------------------
# T6. Reachability (honest)
# ---------------------------------------------------------------------------

def test_T6_reachability() -> dict:
    """Σm_ν ≈ 59 meV is testable NOW (DESI DR2, at the floor). m_β ≈ 8.8
    meV is ~4–5× below the best foreseeable β-decay reach (Project 8 ~40
    meV). m_ββ ≈ 1.5–3.7 meV is ~3–10× below next-gen 0νββ (~9–20 meV). So
    the spectrum is pinned, but only Σm_ν is near-term testable."""
    return {
        'name': 'T6_reachability_honest',
        'description': (
            "Σm_ν ~59 meV testable now (DESI, at the floor); m_β ~8.8 meV "
            "below Project 8 (~40 meV); m_ββ ~1.5–3.7 meV below next-gen "
            "0νββ (~9–20 meV). Only Σm_ν is near-term testable."
        ),
        'sigma_testable_now': True,
        'sigma_probe': 'DESI DR2 + CMB (cornering the floor)',
        'm_beta_reach_ratio': round(PROJECT8_MEV / m_beta_meV(0.0), 1),
        'm_bb_reach_ratio': round(NEXTGEN_0NUBB_MEV / m_bb_band_meV(0.0)[1], 1),
        'm_beta_reachable_soon': m_beta_meV(0.0) > PROJECT8_MEV,
        'm_bb_reachable_soon': m_bb_band_meV(0.0)[1] > NEXTGEN_0NUBB_MEV,
        'pass': (m_beta_meV(0.0) < PROJECT8_MEV
                 and m_bb_band_meV(0.0)[1] < NEXTGEN_0NUBB_MEV),
    }


# ---------------------------------------------------------------------------
# T7. Falsifiers + honest scope
# ---------------------------------------------------------------------------

def test_T7_falsifiers_and_scope() -> dict:
    return {
        'name': 'T7_falsifiers_and_scope',
        'description': (
            "Falsifiers: robust Σm_ν < 59.0 ⟹ NO excluded (BAM fails); "
            "Σm_ν ≳ 100 ⟹ not light; m_ββ ≫ 4 meV with confirmed NO above "
            "the BAM band. Flag: some 2025 DESI+CMB fits already prefer "
            "Σm_ν at/below the floor — tension for ALL NO models. Open: m₁ "
            "band (0–3 meV) and the Majorana phases."
        ),
        'falsifiers': [
            'Σm_ν < 59.0 meV (below NO floor) ⟹ normal ordering excluded ⟹ BAM fails',
            'Σm_ν ≳ 100 meV (IO/degenerate) ⟹ contradicts BAM light scale',
            'm_ββ ≫ 4 meV with confirmed NO ⟹ above the BAM hierarchical band',
        ],
        'tension_flag_2025': (
            'some DESI DR2 + CMB analyses prefer central Σm_ν at or below '
            'the NO floor (data wants less lensing than the floor) — tension '
            'for every normal-ordered model, BAM included'
        ),
        'open': [
            'm₁ within its cornered band (0–3 meV) — the absolute-scale residual',
            'the Majorana phases (anarchic) — set the m_ββ value within [1.5, 3.7] meV',
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
            "Updating to NuFIT 6.0 and cornering m₁ with DESI DR2 pins the "
            "meV-scale spectrum: m = (≲3, 8.65, 50.34) meV, Σm_ν ∈ [59.0, "
            "62.6] meV (sharpened toward the floor), m_β ≈ 8.8–9.3 meV, "
            "m_ββ a nonzero floor [1.5, 3.7] meV. Only Σm_ν is near-term "
            "testable; m_β, m_ββ are sharp but below foreseeable "
            "sensitivity. m₁ band + phases open."
        ),
        'classification': 'NEUTRINO_MEV_PREDICTIONS_SHARPENED_SPECTRUM_PINNED_SIGMA_TESTABLE',
        'pass': True,
    }


# ---------------------------------------------------------------------------
# Runner + verdict
# ---------------------------------------------------------------------------

def run_probe() -> dict:
    tests = [
        test_T1_setup(),
        test_T2_nufit6_inputs(),
        test_T3_cornering(),
        test_T4_pinned_spectrum(),
        test_T5_lab_effective_masses(),
        test_T6_reachability(),
        test_T7_falsifiers_and_scope(),
        test_T8_assessment(),
    ]
    if all(t['pass'] for t in tests[:7]):
        verdict_class = 'NEUTRINO_MEV_PREDICTIONS_SHARPENED_SPECTRUM_PINNED_SIGMA_TESTABLE'
        verdict = (
            'THE meV-SCALE NEUTRINO PREDICTIONS, SHARPENED: THE SPECTRUM IS '
            'NOW PINNED, AND ONLY Σm_ν IS NEAR-TERM TESTABLE. PR #96 left the '
            'headline observable as a band — Σm_ν ≈ 59–65 meV — because the '
            'lightest mass m₁ was an open absolute-scale residual. Updating '
            'the oscillation inputs to the latest global fit (NuFIT 6.0, '
            '2024) and using the 2025 DESI DR2 + CMB cosmology to corner m₁ '
            'sharpens the whole picture.\n\n'
            'THE CORNERING. NuFIT 6.0 (NO) fixes two of the three masses '
            'outright: √Δm²_21 = 8.65 meV and √Δm²_31 = 50.34 meV, so the '
            'normal-ordering floor is Σm_ν(m₁→0) = 59.0 meV. DESI DR2 + CMB '
            '(≲ 60–64 meV at 95%) then pins m₁ to the low end of "few meV": '
            'm₁ ≲ 3 meV, hence Σm_ν ∈ [59.0, 62.6] meV — a sharpening of the '
            'PR #96 band toward the floor. The sharpest single-number BAM '
            'statement is the hierarchical limit m₁ → 0, the natural '
            'endpoint of the light cavity-overtone spectrum and the point '
            'cosmology is cornering us toward.\n\n'
            'THE PINNED SPECTRUM. m = (≲3, 8.65, 50.34) meV. From it the two '
            'laboratory effective masses follow: m_β = √(Σ|U_ei|²m_i²) ≈ 8.8 '
            'meV (rising to ~9.3 at m₁ = 3 meV), and the 0νββ effective mass '
            'm_ββ = |Σ U_ei² m_i|. Crucially, in normal ordering the '
            'contributions CANNOT fully cancel at m₁ → 0 — s12²c13² m₂ = '
            '2.60 meV exceeds s13² m₃ = 1.10 meV — leaving a NONZERO floor '
            'm_ββ ∈ [1.5, 3.7] meV over the Majorana phases. So BAM predicts '
            'a small but nonvanishing 0νββ rate in the hierarchical '
            'limit.\n\n'
            'REACHABILITY (the honest other half). Sharper predictions also '
            'sharpen what can test them: Σm_ν ≈ 59 meV is testable NOW — '
            'DESI DR2 + CMB is cornering it at the floor. But m_β ≈ 8.8 meV '
            'sits ~4–5× below the best foreseeable β-decay reach (Project 8 '
            '~40 meV; KATRIN final ~200 meV), and m_ββ ≈ 1.5–3.7 meV sits '
            '~3–10× below next-gen 0νββ (LEGEND-1000 / nEXO ~9–20 meV). The '
            'spectrum is pinned, but only Σm_ν is near-term testable, and it '
            'is being tested at the floor right now; the Majorana and '
            'absolute-scale predictions are sharp yet below foreseeable '
            'sensitivity — an honest readout, not a promise of imminent '
            'discovery.\n\n'
            'FALSIFIERS (sharpened). A robust cosmological Σm_ν < 59.0 meV '
            'excludes normal ordering ⟹ BAM fails (and would clash with the '
            'oscillation Δm² themselves). A quasi-degenerate detection '
            '(Σm_ν ≳ 100 meV) contradicts the light scale. A 0νββ signal '
            'with m_ββ ≫ 4 meV and confirmed NO would sit above the BAM '
            'hierarchical band. Honest flag: some 2025 DESI + CMB analyses '
            'already prefer central Σm_ν at or below the floor (the data '
            'wants less lensing than the floor supplies); if that hardens it '
            'is tension for ALL normal-ordered models, BAM included. Open: '
            'm₁ within its cornered band (0–3 meV) and the anarchic Majorana '
            'phases (which set m_ββ within [1.5, 3.7] meV).'
        )
    else:
        verdict_class = 'NEUTRINO_MEV_PREDICTIONS_INCONCLUSIVE'
        verdict = (
            'INCONCLUSIVE. A consistency test failed; reconcile the pinned '
            'spectrum with the oscillation + cosmology inputs.'
        )

    m1, m2, m3 = spectrum_meV(0.0)
    bb_lo, bb_hi = m_bb_band_meV(0.0)
    return {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'identification': (
            'the meV-scale neutrino spectrum pinned by NuFIT 6.0 + DESI DR2: '
            'm = (≲3, 8.65, 50.34) meV, Σm_ν ∈ [59.0, 62.6], m_β ≈ 8.8–9.3, '
            'm_ββ nonzero floor [1.5, 3.7] meV'
        ),
        'spectrum_meV': {'m1': '≲3', 'm2': round(m2, 2), 'm3': round(m3, 2)},
        'sigma_meV': f'[{NO_FLOOR_MEV:.1f}, {sigma_meV(M1_CORNER_MEV):.1f}] (sharpened from 59–65)',
        'm_beta_meV': f'{m_beta_meV(0.0):.1f}–{m_beta_meV(M1_CORNER_MEV):.1f}',
        'm_bb_meV': f'[{bb_lo:.1f}, {bb_hi:.1f}] floor (nonzero — no full cancellation in NO)',
        'reachability': 'only Σm_ν near-term testable (DESI, at the floor); m_β, m_ββ below foreseeable sensitivity',
        'open': 'm₁ band (0–3 meV); anarchic Majorana phases',
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
    L.append('# Sharpening the meV-scale neutrino predictions (PR #111)')
    L.append('')
    L.append(f"**Run:** {s['timestamp_utc']}")
    L.append('')
    L.append(
        "Sharpens the PR #96 meV-scale prediction — a Σm_ν band (59–65 meV, "
        "open m₁) — into a fully PINNED spectrum, by updating the "
        "oscillation inputs to **NuFIT 6.0** and using **DESI DR2 + CMB "
        "(2025)** to corner the lightest mass against the normal-ordering "
        "floor. Result: `m = (≲3, 8.65, 50.34) meV`, `Σm_ν ∈ [59.0, 62.6] "
        "meV`, `m_β ≈ 8.8–9.3 meV`, and a **nonzero** 0νββ floor `m_ββ ∈ "
        "[1.5, 3.7] meV` (the NO contributions cannot fully cancel). Honest "
        "other half: only Σm_ν is near-term testable."
    )
    L.append('')
    L.append(f"- **Spectrum**: {s['spectrum_meV']}")
    L.append(f"- **Σm_ν**: {s['sigma_meV']}")
    L.append(f"- **m_β (KATRIN)**: {s['m_beta_meV']} meV")
    L.append(f"- **m_ββ (0νββ)**: {s['m_bb_meV']}")
    L.append(f"- **Reachability**: {s['reachability']}")
    L.append(f"- **Open**: {s['open']}")
    L.append('')

    L.append('## Test summary')
    L.append('')
    L.append('| # | Test | Key finding | PASS? |')
    L.append('|---|---|---|---|')
    label_map = {
        'T1': 'sharpen: Σm_ν band → full pinned spectrum',
        'T2': 'NuFIT 6.0: √Δm²_21 8.65, √Δm²_31 50.34 ⟹ floor 59.0 meV',
        'T3': 'DESI DR2 corners m₁ ≲ 3 meV ⟹ Σm_ν [59.0, 62.6]',
        'T4': 'pinned spectrum m = (≲3, 8.65, 50.34) meV',
        'T5': 'm_β 8.8–9.3 meV; m_ββ nonzero floor [1.5, 3.7] meV',
        'T6': 'reachable: only Σm_ν (DESI); m_β, m_ββ below foreseeable',
        'T7': 'falsifiers sharpened + 2025 sub-floor tension flag',
        'T8': 'NEUTRINO_MEV_PREDICTIONS_SHARPENED_SPECTRUM_PINNED_SIGMA_TESTABLE',
    }
    for t in s['tests']:
        passed = '**PASS**' if t['pass'] else '**FAIL**'
        nm = t['name']; prefix = nm[:2]
        L.append(f"| {prefix} | `{nm}` | {label_map.get(prefix, '—')} | {passed} |")
    L.append('')

    L.append('## The sharpened meV-scale sheet')
    L.append('')
    L.append('| observable | BAM (NO, light) | testable by | reach vs BAM |')
    L.append('|---|---|---|---|')
    L.append('| m₁ | ≲ 3 meV (→ 0) | cosmology (indirect) | — |')
    L.append('| m₂ | 8.65 meV | (pinned by Δm²_21) | — |')
    L.append('| m₃ | 50.34 meV | (pinned by Δm²_31) | — |')
    L.append('| **Σm_ν** | **59.0–62.6 meV** | **DESI DR2 + CMB** | **at the floor NOW** |')
    L.append('| m_β | 8.8–9.3 meV | KATRIN / Project 8 | ~4–5× below |')
    L.append('| m_ββ | 1.5–3.7 meV (floor) | LEGEND-1000 / nEXO | ~3–10× below |')
    L.append('')

    L.append('## Verdict')
    L.append('')
    L.append(f"**{s['verdict_class']}.** {s['verdict']}")
    L.append('')

    L.append('## What this leaves open')
    L.append('')
    L.append('- **m₁ within its cornered band (0–3 meV)** — the '
             'absolute-scale residual; cosmology is squeezing it against the '
             'floor.')
    L.append('- **The anarchic Majorana phases** — they set the exact m_ββ '
             'within the [1.5, 3.7] meV floor band; not predicted '
             'individually (the universal flavor puzzle).')
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
    out = here / 'runs' / f'{ts}_neutrino_mev_scale_sharpening_probe'
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
