"""
Mixing/anarchy origin of the ε_n profile (PR #151).

> Framing: QFT on the *fixed classical* throat geometry (geometry → fields), not
> quantum gravity. The seesaw matrix lives in the overtone basis of the
> classical cavity; the audit asks whether mixing reshapes its spectrum.

PR #149 bracketed the ε_n overshoot and concluded the gentle required profile
"plausibly belongs to the mixing/anarchy sector (#92)". This probe TESTS that
hypothesis. The model: the Majorana seesaw matrix in the overtone basis, with
the DERIVED steep χ-driven compliances (ε_n ∝ 1/χ_n, #91/#113) setting the
per-channel bounce suppressions, the DERIVED cavity floors (#91) setting the
Dirac growth, anarchic O(1) cross-channel overlaps (the #91 large-PMNS
structure), and one structural choice β interpolating the pair-tunneling rule
from FACTORIZED (β = 0: each element carries the geometric mean of its two
channel suppressions) to CHANNEL-DOMINANT (β = 1: each element tunnels
through the most compliant neck available to the pair — the least-action
saddle).

## The result, in three parts

1. CHANNEL DOMINANCE RESOLVES THE MEASURED OVERSHOOT. The only measured mass
   ratio is m₃/m₂ ≈ 5.7 (both splittings; m₁ is unmeasured). The factorized
   rule (β = 0) reproduces the #113/#149 overshoot in matrix form: ensemble
   median m₃/m₂ ≈ 113, the observed value at the 0.1th percentile — EXCLUDED.
   Channel dominance (β = 1) compresses it to median ≈ 3.0 with the observed
   value at the ~75th percentile — NATURAL. The steep χ-driven hierarchy
   collapses out of the heavy pair because every pair element shares the
   widest neck.

2. THE SAME β PRODUCES LARGE MIXING. The mixing indicator (the
   second-largest component of the heaviest eigenvector) grows from 0.085
   (β = 0, aligned) to 0.43 (β = 1, large mixing ~ sin θ ≳ 0.4) — one
   mechanism moves both observables toward the data, consistent with the #91
   cross-channel large-PMNS identification.

3. THE UNMEASURED RATIO BECOMES A PREDICTION. m₂/m₁ stays large (median
   ≈ 200) at every β: the lightest (most throat-coupled) channel keeps its
   suppression. Since m₁ is NOT measured, this is not a misfit — it is a
   falsifiable PREDICTION: m₁ ≈ 0.04 meV (vs the #112 uniform-anchor
   2.08 meV), hence Σm_ν ≈ 58.8 meV (vs 61.1) — distinguishable by ~1 meV
   cosmology precision, with m_ββ shifted accordingly. Both scenarios keep
   normal ordering.

## The residual relocation

The #149 residual — a three-number fine-tuned ε_n profile — relocates to:
the derived χ-driven compliances + the derived cavity floors + ONE structural
saddle rule (channel dominance) + an anarchic O(1) draw. The price is that
the precise ratios become STATISTICAL (anarchy percentiles), not
deterministic — the flavor puzzle's BAM face, now localized. No new
continuous knob: β is a discrete structural choice between two saddle rules,
selected by the measured ratio.

## Scope

The pair-tunneling rule G_ij is MODELLED (motivated as the least-action
saddle through the widest available neck — the #136 posture); the
charged-lepton side of the PMNS matrix is not computed (the mixing indicator
is ν-side only); CP/Majorana phases open. The m₁ prediction inherits the
ensemble spread.

Tests:
  T1. Goal: test #149's mixing/anarchy hypothesis with a derived-input model.
  T2. The model: derived inputs verified (cavity floors = the #149 m_D
      endpoint to <1%; χ-driven suppressions); the β-interpolated rule.
  T3. The β scan on the measured ratio: β = 0 excluded (0.1th pct — the
      overshoot re-derived in matrix form); β = 1 natural (~75th pct).
  T4. Large mixing from the same β: indicator 0.085 → 0.43.
  T5. The m₁ prediction: r₂₁ ≈ 200 ⟹ m₁ ≈ 0.04 meV, Σm_ν ≈ 58.8 meV vs the
      uniform-anchor 61.1 — a falsifiable discriminator.
  T6. Residual relocation: three-number profile → one saddle rule + anarchic
      draw; naturalness percentiles; input budget unchanged.
  T7. Scope.
  T8. Assessment.

Verdict:
  - EPS_N_FROM_CHANNEL_DOMINANT_ANARCHY_R32_NATURAL_M1_PREDICTED_LIGHT
    (expected): channel-dominant anarchic mixing in the overtone basis
    resolves the measured overshoot (m₃/m₂ natural at the ~75th percentile),
    produces large mixing from the same mechanism, and converts the
    unmeasured m₂/m₁ into a falsifiable prediction (m₁ ≈ 0.04 meV,
    Σm_ν ≈ 58.8 meV) — relocating the #149 residual from a fine-tuned
    profile to one structural saddle rule plus an anarchic draw.
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
P_BOUNCE = 4.8                                  # #112
EPS_RATIO = np.array([1.0, 3.134, 7.787])       # χ-driven ε_n ∝ 1/χ_n (#113)
MD_FLOORS = np.array([1.055, 1.974, 2.894])     # cavity floors (#91)
MD_149_ENDPOINT = (1.879, 1.478)                # the #149 m_D attribution
LNF = P_BOUNCE * np.log(EPS_RATIO)              # ln(suppression release) per channel

# Oscillation data (#149 conventions)
DM21_SQ = 7.42e-5      # eV²
DM31_SQ = 2.514e-3     # eV²
M1_UNIFORM_ANCHOR = 2.08   # meV (#112)

N_ENSEMBLE = 2000
SEED = 42


def observed_r32() -> float:
    """The measured heavy-pair ratio (m₁-insensitive for small m₁)."""
    m2 = math.sqrt(M1_UNIFORM_ANCHOR**2 + DM21_SQ * 1e6)
    m3 = math.sqrt(M1_UNIFORM_ANCHOR**2 + DM31_SQ * 1e6)
    return m3 / m2


def seesaw_ensemble(beta: float, n: int = N_ENSEMBLE, seed: int = SEED):
    """Ensemble of overtone-basis seesaw matrices
    M_ij = m_D,i m_D,j · c_ij · G_ij(β), with c_ij symmetric anarchic O(1)
    (log-uniform in [1/3, 3], random sign) and
    G_ij = exp[(1−β)·(lnF_i + lnF_j)/2 + β·max(lnF_i, lnF_j)] — factorized
    (β = 0) to channel-dominant (β = 1) pair tunneling.
    Returns (r21, r32, mixing-indicator) arrays."""
    rng = np.random.default_rng(seed)
    G = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            G[i, j] = math.exp((1 - beta) * (LNF[i] + LNF[j]) / 2
                               + beta * max(LNF[i], LNF[j]))
    r21s, r32s, mixes = [], [], []
    for _ in range(n):
        c = np.exp(rng.uniform(-math.log(3), math.log(3), (3, 3)))
        c = (c + c.T) / 2
        sgn = rng.choice([-1.0, 1.0], (3, 3))
        sgn = np.triu(sgn) + np.triu(sgn, 1).T
        M = np.outer(MD_FLOORS, MD_FLOORS) * (c * sgn) * G
        w, V = np.linalg.eigh(M)
        m = np.sort(np.abs(w))
        if m[0] <= 0:
            continue
        r21s.append(m[1] / m[0])
        r32s.append(m[2] / m[1])
        v_heavy = np.abs(V[:, int(np.argmax(np.abs(w)))])
        mixes.append(float(np.sort(v_heavy)[1]))
    return np.array(r21s), np.array(r32s), np.array(mixes)


_SCAN = {}
for _beta in (0.0, 0.5, 1.0):
    _SCAN[_beta] = seesaw_ensemble(_beta)


# ---------------------------------------------------------------------------
# T1. Goal
# ---------------------------------------------------------------------------

def test_T1_goal() -> dict:
    return {
        'name': 'T1_goal',
        'description': (
            "Test #149's hypothesis that the gentle ε_n profile belongs to "
            "the mixing/anarchy sector: model the seesaw in the overtone "
            "basis with DERIVED inputs (χ-driven compliances, cavity-floor "
            "Dirac growth, anarchic cross-channel overlaps) and audit the "
            "pair-tunneling saddle rule from factorized to channel-dominant."
        ),
        'builds_on': ['#149 ε_n bracket (the hypothesis)', '#113 overshoot',
                      '#91 cavity floors + cross-channel PMNS', '#79 χ_n',
                      '#112 p = 4.8 bounce', '#92 mixing/anarchy sector'],
        'framing': 'QFT on the fixed classical throat geometry — not quantum gravity',
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T2. The model and its derived inputs
# ---------------------------------------------------------------------------

def test_T2_model_inputs() -> dict:
    """The Dirac growth is the DERIVED cavity floors (#91) — which reproduce
    the #149 m_D-attribution endpoint to <1%; the channel suppressions are
    the DERIVED χ-driven compliances through the p = 4.8 bounce."""
    md21 = MD_FLOORS[1] / MD_FLOORS[0]
    md32 = MD_FLOORS[2] / MD_FLOORS[1]
    d21 = abs(md21 / MD_149_ENDPOINT[0] - 1.0)
    d32 = abs(md32 / MD_149_ENDPOINT[1] - 1.0)
    F = np.exp(LNF)
    ok = d21 < 0.01 and d32 < 0.01
    return {
        'name': 'T2_model_and_derived_inputs',
        'description': (
            "M_ij = m_D,i m_D,j · c_ij · G_ij(β). Every hierarchical input "
            "is DERIVED: the Dirac growth is the #91 cavity floors — which "
            "equal the #149 m_D-attribution endpoint to <1% (an identity "
            "between two independently-introduced numbers, verified) — and "
            "the channel suppressions F_n = (ε_n/ε₁)^p are the χ-driven "
            "compliances through the #112 bounce. The anarchic c_ij are the "
            "#91 cross-channel O(1) overlaps; β interpolates the pair "
            "tunneling from factorized to channel-dominant (least-action "
            "saddle through the widest neck — modelled, the #136 posture)."
        ),
        'cavity_floor_ratios': [round(md21, 4), round(md32, 4)],
        'pr149_md_endpoint': list(MD_149_ENDPOINT),
        'agreement': [float(f'{d21:.2e}'), float(f'{d32:.2e}')],
        'channel_release_factors_F': [round(float(f), 1) for f in F],
        'pass': ok,
    }


# ---------------------------------------------------------------------------
# T3. The β scan on the measured ratio
# ---------------------------------------------------------------------------

def test_T3_beta_scan_measured_ratio() -> dict:
    """β = 0 reproduces the #113/#149 overshoot in matrix form (observed r₃₂
    at the ~0.1th percentile — excluded); β = 1 compresses the spread
    (observed at the ~75th percentile — natural)."""
    obs = observed_r32()
    rows, percentiles = [], {}
    for beta, (r21s, r32s, mixes) in _SCAN.items():
        pct = float(np.mean(r32s < obs) * 100)
        percentiles[beta] = pct
        rows.append({'beta': beta,
                     'median_r32': round(float(np.median(r32s)), 2),
                     'observed_r32': round(obs, 2),
                     'observed_percentile': round(pct, 1)})
    excluded_0 = percentiles[0.0] < 1.0
    natural_1 = 20.0 < percentiles[1.0] < 95.0
    return {
        'name': 'T3_beta_scan_measured_ratio',
        'description': (
            "The only measured mass ratio is the heavy pair m₃/m₂ ≈ 5.7 "
            "(both Δm²; m₁ unmeasured). FACTORIZED (β = 0): ensemble median "
            "r₃₂ ≈ 113 — the #113/#149 overshoot re-derived in matrix form; "
            "the observed value sits at the 0.1th percentile, EXCLUDED. "
            "CHANNEL-DOMINANT (β = 1): every pair element tunnels through "
            "the widest neck, the steep hierarchy collapses out of the "
            "heavy pair, median r₃₂ ≈ 3.0, observed at the ~75th percentile "
            "— NATURAL. Mixing resolves the measured overshoot."
        ),
        'rows': rows,
        'pass': excluded_0 and natural_1,
    }


# ---------------------------------------------------------------------------
# T4. Large mixing from the same mechanism
# ---------------------------------------------------------------------------

def test_T4_large_mixing_emerges() -> dict:
    """The mixing indicator grows 0.085 → 0.43 as β: 0 → 1 — the same saddle
    rule that fixes the spread produces the large mixing (#91 cross-channel
    consistency)."""
    rows = []
    for beta, (r21s, r32s, mixes) in _SCAN.items():
        rows.append({'beta': beta,
                     'median_mixing_indicator': round(float(np.median(mixes)), 3)})
    m0 = rows[0]['median_mixing_indicator']
    m1 = rows[-1]['median_mixing_indicator']
    ok = m0 < 0.15 and m1 > 0.35
    return {
        'name': 'T4_large_mixing_from_same_beta',
        'description': (
            "The mixing indicator (second-largest component of the heaviest "
            "mass eigenvector) grows from 0.085 (β = 0: the factorized "
            "hierarchy aligns the eigenvectors with the overtones — SMALL "
            "mixing, contradicting #91) to 0.43 (β = 1: large mixing, "
            "sin θ ≳ 0.4 — the anarchy/large-PMNS regime). ONE structural "
            "rule moves BOTH observables toward the data: the spread "
            "compresses and the mixing grows together — consistent with the "
            "#91 cross-channel large-PMNS identification."
        ),
        'rows': rows,
        'pass': ok,
    }


# ---------------------------------------------------------------------------
# T5. The m₁ prediction
# ---------------------------------------------------------------------------

def test_T5_m1_prediction() -> dict:
    """r₂₁ stays ≈ 200 at β = 1 — and m₁ is unmeasured, so this is a
    PREDICTION: m₁ ≈ 0.04 meV, Σm_ν ≈ 58.8 meV (vs the #112 uniform-anchor
    2.08 meV, Σ = 61.1) — a ~2 meV discriminator for cosmology."""
    r21s, r32s, _ = _SCAN[1.0]
    r21_med = float(np.median(r21s))
    m2 = math.sqrt(DM21_SQ) * 1e3       # m₁ ≈ 0 limit, meV
    m3 = math.sqrt(DM31_SQ) * 1e3
    m1_pred = m2 / r21_med
    sigma_pred = m1_pred + m2 + m3
    m2_u = math.sqrt(M1_UNIFORM_ANCHOR**2 + DM21_SQ * 1e6)
    m3_u = math.sqrt(M1_UNIFORM_ANCHOR**2 + DM31_SQ * 1e6)
    sigma_uniform = M1_UNIFORM_ANCHOR + m2_u + m3_u
    ok = r21_med > 50 and m1_pred < 0.2 and abs(sigma_pred - 58.8) < 1.0
    return {
        'name': 'T5_m1_prediction',
        'description': (
            "The lightest channel (most throat-coupled, least compliant) "
            "keeps its suppression under every β: median r₂₁ ≈ 200 at "
            "β = 1. m₁ is NOT measured (only the two Δm²), so this is not a "
            "misfit — it is a falsifiable PREDICTION distinguishing the "
            "mixing solution from the #112 uniform anchor: m₁ ≈ 0.04 meV "
            "(vs 2.08), Σm_ν ≈ 58.8 meV (vs 61.1) — a ~2 meV discriminator "
            "at next-generation cosmology precision, with m_ββ shifted "
            "accordingly. Both scenarios keep normal ordering."
        ),
        'median_r21_beta1': round(r21_med, 1),
        'm1_predicted_meV': round(m1_pred, 3),
        'sigma_mnu_mixing_meV': round(sigma_pred, 1),
        'm1_uniform_anchor_meV': M1_UNIFORM_ANCHOR,
        'sigma_mnu_uniform_meV': round(sigma_uniform, 1),
        'discriminator': 'Σm_ν at ~1–2 meV precision; m_ββ',
        'pass': ok,
    }


# ---------------------------------------------------------------------------
# T6. The residual relocation
# ---------------------------------------------------------------------------

def test_T6_residual_relocation() -> dict:
    """The #149 three-number profile residual relocates to one structural
    saddle rule + an anarchic O(1) draw; the ratios become statistical
    (percentile-natural), not deterministic; no new continuous knob."""
    obs = observed_r32()
    _, r32s, _ = _SCAN[1.0]
    pct = float(np.mean(r32s < obs) * 100)
    iqr = [round(float(np.percentile(r32s, 25)), 2),
           round(float(np.percentile(r32s, 75)), 2)]
    return {
        'name': 'T6_residual_relocation',
        'description': (
            "Before (#149): a three-number fine-tuned ε_n profile, bracketed "
            "but underived. After: the DERIVED χ-driven compliances + the "
            "DERIVED cavity floors + ONE discrete structural choice (the "
            "channel-dominant saddle rule, selected by the measured ratio) "
            "+ an anarchic O(1) draw. The price: the precise ratios are "
            "STATISTICAL — the observed r₃₂ sits at the ~75th percentile "
            "(natural, IQR shown) rather than being predicted exactly — the "
            "flavor puzzle's BAM face, now localized to the anarchic draw. "
            "No new continuous knob enters: β is a discrete saddle-rule "
            "choice, not a tuned parameter; the input budget (#150) is "
            "unchanged."
        ),
        'observed_r32': round(obs, 2),
        'beta1_r32_iqr': iqr,
        'observed_percentile': round(pct, 1),
        'relocated_to': 'channel-dominant saddle rule + anarchic O(1) draw',
        'budget_unchanged': True,
        'pass': 20.0 < pct < 95.0,
    }


# ---------------------------------------------------------------------------
# T7. Scope
# ---------------------------------------------------------------------------

def test_T7_scope() -> dict:
    return {
        'name': 'T7_scope',
        'description': (
            "The pair-tunneling rule G_ij is MODELLED (motivated as the "
            "least-action saddle through the widest neck available to the "
            "pair; deriving it from the bounce path integral is open — the "
            "#136 posture). The mixing indicator is ν-side only (the "
            "charged-lepton rotation and explicit PMNS angles are not "
            "computed); CP/Majorana phases open. The m₁ prediction carries "
            "the ensemble spread (median quoted). The hypothesis test is "
            "POSITIVE on the measured ratio and converts the unmeasured one "
            "into a discriminator — it does not yet derive the anarchic "
            "draw."
        ),
        'open': [
            'derive G_ij from the bounce path integral (saddle rule)',
            'charged-lepton rotation / explicit PMNS angles; CP phases',
            'the anarchic draw (the flavor puzzle, localized)',
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
            "Channel-dominant anarchic mixing in the overtone basis resolves "
            "the measured overshoot (r₃₂ natural at the ~75th percentile; "
            "factorized excluded at 0.1%), produces large mixing from the "
            "same rule, and converts the unmeasured m₂/m₁ into a falsifiable "
            "prediction (m₁ ≈ 0.04 meV, Σm_ν ≈ 58.8 meV vs the uniform "
            "anchor's 61.1). The #149 residual relocates from a fine-tuned "
            "profile to one saddle rule plus an anarchic draw."
        ),
        'classification': 'EPS_N_FROM_CHANNEL_DOMINANT_ANARCHY_R32_NATURAL_M1_PREDICTED_LIGHT',
        'pass': True,
    }


# ---------------------------------------------------------------------------
# Runner + verdict
# ---------------------------------------------------------------------------

def run_probe() -> dict:
    tests = [
        test_T1_goal(),
        test_T2_model_inputs(),
        test_T3_beta_scan_measured_ratio(),
        test_T4_large_mixing_emerges(),
        test_T5_m1_prediction(),
        test_T6_residual_relocation(),
        test_T7_scope(),
        test_T8_assessment(),
    ]
    t3, t5 = tests[2], tests[4]
    if all(t['pass'] for t in tests[:7]):
        verdict_class = 'EPS_N_FROM_CHANNEL_DOMINANT_ANARCHY_R32_NATURAL_M1_PREDICTED_LIGHT'
        verdict = (
            'THE #149 HYPOTHESIS TESTS POSITIVE: CHANNEL-DOMINANT ANARCHIC '
            'MIXING IN THE OVERTONE BASIS RESOLVES THE MEASURED OVERSHOOT, '
            'PRODUCES LARGE MIXING FROM THE SAME RULE, AND CONVERTS THE '
            'UNMEASURED m₂/m₁ INTO A FALSIFIABLE PREDICTION — THE ε_n '
            'PROFILE RESIDUAL RELOCATES TO ONE SADDLE RULE PLUS AN ANARCHIC '
            'DRAW. #149 bracketed the overshoot and pointed at the '
            'mixing/anarchy sector; this probe builds the test.\n\n'
            'THE MODEL, DERIVED INPUTS. M_ij = m_D,i m_D,j·c_ij·G_ij(β): '
            'the Dirac growth is the #91 cavity floors — which equal the '
            '#149 m_D-attribution endpoint to <1% (verified identity) — the '
            'channel suppressions are the χ-driven compliances through the '
            '#112 bounce, the c_ij are the #91 anarchic cross-channel '
            'overlaps, and β interpolates the pair-tunneling saddle from '
            'factorized to channel-dominant (the widest neck available to '
            'the pair).\n\n'
            'THE MEASURED RATIO SELECTS CHANNEL DOMINANCE. The only '
            'measured ratio is the heavy pair r₃₂ ≈ 5.7. Factorized '
            '(β = 0): ensemble median ≈ 113 — the #113/#149 overshoot '
            're-derived in matrix form; observed at the 0.1th percentile, '
            'excluded. Channel-dominant (β = 1): the steep hierarchy '
            'collapses out of the heavy pair (every element shares the '
            'widest neck); median ≈ 3.0; observed at the ~75th percentile — '
            'natural.\n\n'
            'ONE RULE, TWO OBSERVABLES. The same β that compresses the '
            'spread grows the mixing indicator 0.085 → 0.43 (small-mixing '
            'aligned → large-mixing anarchic) — consistent with the #91 '
            'cross-channel large-PMNS identification. The spread and the '
            'mixing are two faces of one saddle rule.\n\n'
            'THE m₁ PREDICTION. The lightest channel keeps its suppression: '
            f'median r₂₁ ≈ {t5["median_r21_beta1"]} at β = 1. m₁ is '
            'unmeasured, so this is a PREDICTION, not a misfit: '
            f'm₁ ≈ {t5["m1_predicted_meV"]} meV and '
            f'Σm_ν ≈ {t5["sigma_mnu_mixing_meV"]} meV, against the #112 '
            f'uniform anchor (m₁ = 2.08, Σ = {t5["sigma_mnu_uniform_meV"]}) '
            '— a ~2 meV discriminator for next-generation cosmology, with '
            'm_ββ shifted accordingly; both keep normal ordering.\n\n'
            'THE RESIDUAL RELOCATES. Before: a three-number fine-tuned '
            'profile (bracketed, #149). After: derived compliances + '
            'derived floors + one discrete saddle rule + an anarchic O(1) '
            'draw — the ratios become percentile-natural statistics rather '
            'than deterministic outputs (the flavor puzzle\'s BAM face, '
            'localized). No new continuous knob; the #150 budget is '
            'unchanged.\n\n'
            'SCOPE. The saddle rule is modelled, not derived from the '
            'bounce path integral; the charged-lepton rotation, explicit '
            'PMNS angles, and CP phases are open; the anarchic draw itself '
            'is the localized residual.'
        )
    else:
        verdict_class = 'EPS_N_MIXING_ANARCHY_INCONCLUSIVE'
        verdict = (
            'INCONCLUSIVE. A hypothesis-test check failed; review the '
            'ensemble scan, the derived-input identities, or the prediction '
            'extraction.'
        )

    return {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'identification': (
            'channel-dominant anarchic mixing resolves the measured '
            'overshoot (r₃₂ natural at ~75th percentile; factorized '
            'excluded), grows large mixing from the same saddle rule, and '
            'predicts m₁ ≈ 0.04 meV / Σm_ν ≈ 58.8 meV — the #149 residual '
            'relocates to one saddle rule plus an anarchic draw'
        ),
        'model': 'M_ij = m_D,i m_D,j·c_ij·G_ij(β); floors #91, χ-compliances #113, anarchic c #91',
        'selection': 'β = 0 excluded (0.1th pct); β = 1 natural (~75th pct) on the measured r₃₂',
        'mixing': 'indicator 0.085 → 0.43 with the same β (cross-channel #91)',
        'prediction': 'm₁ ≈ 0.04 meV; Σm_ν ≈ 58.8 meV (vs uniform-anchor 61.1) — falsifiable',
        'relocation': 'three-number profile → one saddle rule + anarchic O(1) draw; budget unchanged',
        'open': 'derive the saddle rule; PMNS angles/CP; the anarchic draw (flavor puzzle)',
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
    out.append('# Mixing/anarchy origin of the ε_n profile (PR #151)')
    out.append('')
    out.append(f"**Run:** {s['timestamp_utc']}")
    out.append('')
    out.append(
        "Tests #149's hypothesis that the gentle ε_n profile belongs to the "
        "mixing/anarchy sector. The seesaw matrix in the overtone basis — "
        "derived χ-driven compliances, derived cavity-floor Dirac growth, "
        "anarchic cross-channel overlaps — with the pair-tunneling saddle "
        "rule audited from factorized to channel-dominant. Result: channel "
        "dominance resolves the measured overshoot (natural at the ~75th "
        "percentile), grows large mixing from the same rule, and converts "
        "the unmeasured m₂/m₁ into a falsifiable m₁ prediction. *(QFT on "
        "the classical throat, not quantum gravity.)*"
    )
    out.append('')
    out.append(f"- **Model**: {s['model']}")
    out.append(f"- **Selection**: {s['selection']}")
    out.append(f"- **Mixing**: {s['mixing']}")
    out.append(f"- **Prediction**: {s['prediction']}")
    out.append(f"- **Relocation**: {s['relocation']}")
    out.append(f"- **Open**: {s['open']}")
    out.append('')

    out.append('## Test summary')
    out.append('')
    out.append('| # | Test | Key finding | PASS? |')
    out.append('|---|---|---|---|')
    label_map = {
        'T1': 'test the #149 mixing/anarchy hypothesis with derived inputs',
        'T2': 'cavity floors = the #149 m_D endpoint (<1%, verified identity)',
        'T3': 'measured r₃₂: β = 0 excluded (0.1th pct); β = 1 natural (~75th)',
        'T4': 'large mixing from the same β: indicator 0.085 → 0.43',
        'T5': 'prediction: m₁ ≈ 0.04 meV, Σm_ν ≈ 58.8 meV (vs 61.1 uniform)',
        'T6': 'residual relocates: profile → saddle rule + anarchic draw',
        'T7': 'scope: saddle rule modelled; PMNS angles/CP open',
        'T8': 'EPS_N_FROM_CHANNEL_DOMINANT_ANARCHY_R32_NATURAL_M1_PREDICTED_LIGHT',
    }
    for t in s['tests']:
        passed = '**PASS**' if t['pass'] else '**FAIL**'
        nm = t['name']; prefix = nm[:2]
        out.append(f"| {prefix} | `{nm}` | {label_map.get(prefix, '—')} | {passed} |")
    out.append('')

    t3 = s['tests'][2]
    out.append('## The β scan on the measured ratio (m₃/m₂)')
    out.append('')
    out.append('| β (saddle rule) | median r₃₂ | observed r₃₂ | observed percentile |')
    out.append('|---:|---:|---:|---:|')
    for r in t3['rows']:
        out.append(f"| {r['beta']} | {r['median_r32']} | {r['observed_r32']} | {r['observed_percentile']}% |")
    out.append('')
    out.append("β = 0 (factorized) re-derives the #113/#149 overshoot in "
               "matrix form — excluded; β = 1 (channel-dominant: every pair "
               "element tunnels through the widest neck) makes the observed "
               "ratio anarchy-typical.")
    out.append('')

    t4 = s['tests'][3]
    out.append('## Large mixing emerges from the same rule')
    out.append('')
    out.append('| β | median mixing indicator |')
    out.append('|---:|---:|')
    for r in t4['rows']:
        out.append(f"| {r['beta']} | {r['median_mixing_indicator']} |")
    out.append('')

    t5 = s['tests'][4]
    out.append('## The m₁ prediction (the unmeasured ratio becomes falsifiable)')
    out.append('')
    out.append('| quantity | mixing solution (β = 1) | uniform anchor (#112) |')
    out.append('|---|---:|---:|')
    out.append(f"| m₁ (meV) | {t5['m1_predicted_meV']} | {t5['m1_uniform_anchor_meV']} |")
    out.append(f"| Σm_ν (meV) | {t5['sigma_mnu_mixing_meV']} | {t5['sigma_mnu_uniform_meV']} |")
    out.append('')
    out.append(f"Median r₂₁ = {t5['median_r21_beta1']} at β = 1; the "
               "discriminator is Σm_ν at ~1–2 meV cosmology precision (and "
               "m_ββ); both scenarios keep normal ordering.")
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
    out = here / 'runs' / f'{ts}_eps_n_mixing_anarchy_origin_probe'
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
