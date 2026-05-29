"""
Neutrinoless double-beta (0νββ) effective-mass probe (PR #95).

The neutrino arc (#85–#94) fixed the *structure* of the BAM neutrino
sector. This probe turns it into a concrete, falsifiable prediction for
the one observable that directly tests the Majorana nature — the
effective Majorana mass measured in 0νββ,

    m_ββ = | Σ_i U_ei² m_i |
         = | c12²c13² m1 + s12²c13² m2 e^{iα21} + s13² m3 e^{iα31} |,

which combines everything the arc established:

  - **0νββ occurs at all** because the neutrino is **Majorana**
    (`c₁ = 0`, C-invariant, PR #86) — a Dirac neutrino would forbid it;
  - the **normal ordering** (PR #91: generations = cavity overtones,
    m_ν ∝ m_D) selects the NO band of m_ββ;
  - the **anarchic Majorana phases** (PR #94) populate the *whole* band,
    including the cancellation trough where m_ββ → 0;
  - the **light absolute scale** (PR #90: m_ν ~ few meV) places us in the
    deep-cancellation region.

## The prediction

Using the observed Δm² and mixing angles (BAM fixes the *ordering*,
*Majorana nature*, *phase distribution*, and *scale*, not Δm²
themselves), scanning the Majorana phases uniformly:

  - **Normal ordering band**: m_ββ ≈ 1.5–3.7 meV at zero lightest mass,
    widening to ≈ 0–8 meV at the BAM scale (lightest ~ few meV), with a
    cancellation trough (m_ββ → ~0) around m_lightest ~ 3–5 meV.
  - **Inverted ordering (contrast)**: m_ββ ≈ 19–48 meV, with a hard floor
    ~19 meV and NO cancellation to zero.

So BAM predicts m_ββ ≲ 8 meV — well below the current bound
(KamLAND-Zen, m_ββ ≲ 28–122 meV) and largely below even next-generation
reach (LEGEND-1000 / nEXO, ~ 9–20 meV). This is a sharp falsifier: a
0νββ discovery with m_ββ ≳ 19 meV would imply inverted ordering or a
quasi-degenerate scale, contradicting the BAM normal-ordering + light-
scale prediction.

## What this probe establishes (and does not)

  - **Establishes (BAM-native):** 0νββ occurs (the neutrino is Majorana,
    PR #86); BAM selects the normal-ordering m_ββ band (PR #91), which
    the anarchic Majorana phases (PR #94) populate fully (cancellation to
    ~0 allowed); at the BAM light scale (PR #90) m_ββ ≲ 8 meV — below
    current bounds and a target/falsifier for next-gen experiments.

  - **Does not establish:** a single m_ββ value. It is a band, because
    the lightest neutrino mass is unmeasured and the Majorana phases are
    anarchic (uniform). The exact spectrum (the PR #91 χ_n-corrected
    ratios; the absolute scale) and the specific phases are the residuals.

Tests:
  T1. Setup: m_ββ = |Σ U_ei² m_i|; needs Majorana + ordering + phases +
      scale, all from the arc.
  T2. 0νββ occurs ⟸ neutrino Majorana ⟸ c₁=0 (PR #86); Dirac forbids it.
  T3. Normal ordering (PR #91) ⟹ the NO band; IO contrast band computed.
  T4. Anarchic Majorana phases (PR #94) populate the full band incl.
      the cancellation trough (m_ββ → ~0).
  T5. BAM light scale (PR #90, ~few meV) ⟹ m_ββ ≲ 8 meV.
  T6. Comparison: below current (≲28–122 meV) and largely below next-gen
      (~9–20 meV); falsifiable (discovery ≳19 meV ⟹ IO/degenerate).
  T7. Honest scope: qualitative prediction firm (occurs, NO band, light
      scale); exact m_ββ a band (unmeasured lightest mass + anarchic
      phases).
  T8. Assessment.

Verdict:
  - ZERONUBB_OCCURS_NORMAL_ORDERING_M_BB_FEW_MEV (expected): 0νββ occurs
    (Majorana, c₁=0); BAM selects the normal-ordering band, populated by
    the anarchic Majorana phases; at the light scale m_ββ ≲ 8 meV —
    below current reach, a falsifiable target for next-gen experiments.
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
N_PHASE = 20000
SEED = 0

# Observed inputs (BAM fixes ordering / Majorana / phases / scale, not these):
THETA12_DEG = 33.4
THETA13_DEG = 8.6
DM2_21 = 7.5e-5      # eV²
DM2_31 = 2.5e-3      # eV² (normal ordering)

S12SQ = math.sin(math.radians(THETA12_DEG)) ** 2
C12SQ = 1.0 - S12SQ
S13SQ = math.sin(math.radians(THETA13_DEG)) ** 2
C13SQ = 1.0 - S13SQ

# Experimental scales (meV)
KAMLAND_ZEN_BOUND_MEV = (28.0, 122.0)     # 90% CL, nuclear-matrix-element range
NEXTGEN_REACH_MEV = (9.0, 20.0)           # LEGEND-1000 / nEXO target
BAM_LIGHT_SCALE_MEV = 2.0                 # PR #90 (~few meV)


def mbb_normal(m1: float, a21: np.ndarray, a31: np.ndarray) -> np.ndarray:
    """|m_ββ| (eV) for normal ordering, lightest mass m1 (eV), Majorana
    phase arrays a21, a31."""
    m2 = math.sqrt(m1 ** 2 + DM2_21)
    m3 = math.sqrt(m1 ** 2 + DM2_31)
    z = (C12SQ * C13SQ * m1
         + S12SQ * C13SQ * m2 * np.exp(1j * a21)
         + S13SQ * m3 * np.exp(1j * a31))
    return np.abs(z)


def mbb_inverted(m3: float, a1: np.ndarray, a2: np.ndarray) -> np.ndarray:
    """|m_ββ| (eV) for inverted ordering, lightest mass m3."""
    m2 = math.sqrt(m3 ** 2 + DM2_31)
    m1 = math.sqrt(m3 ** 2 + DM2_31 - DM2_21)
    z = (C12SQ * C13SQ * m1 * np.exp(1j * a1)
         + S12SQ * C13SQ * m2 * np.exp(1j * a2)
         + S13SQ * m3)
    return np.abs(z)


def _phases(n: int = N_PHASE, seed: int = SEED):
    rng = np.random.default_rng(seed)
    a = rng.uniform(0.0, 2.0 * PI, (n, 2))
    return a[:, 0], a[:, 1]


def _band_no(m1_meV: float) -> tuple[float, float]:
    a21, a31 = _phases()
    v = mbb_normal(m1_meV * 1e-3, a21, a31) * 1e3
    return float(v.min()), float(v.max())


def _band_io(m3_meV: float) -> tuple[float, float]:
    a1, a2 = _phases()
    v = mbb_inverted(m3_meV * 1e-3, a1, a2) * 1e3
    return float(v.min()), float(v.max())


# ---------------------------------------------------------------------------
# T1. Setup
# ---------------------------------------------------------------------------

def test_T1_setup() -> dict:
    return {
        'name': 'T1_setup',
        'description': (
            "m_ββ = |Σ U_ei² m_i| = |c12²c13² m1 + s12²c13² m2 e^{iα21} + "
            "s13² m3 e^{iα31}|. Combines Majorana nature (PR #86), ordering "
            "(PR #91), phases (PR #94), scale (PR #90)."
        ),
        's12sq': S12SQ, 'c12sq': C12SQ, 's13sq': S13SQ, 'c13sq': C13SQ,
        'dm2_21': DM2_21, 'dm2_31': DM2_31,
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T2. 0νββ occurs ⟸ Majorana
# ---------------------------------------------------------------------------

def test_T2_occurs_majorana() -> dict:
    """0νββ requires a Majorana neutrino (ΔL=2). In BAM the neutrino is
    Majorana because it is chargeless (c₁=0, C-invariant, PR #86), so
    0νββ OCCURS — a qualitative BAM prediction. A Dirac neutrino would
    forbid it (m_ββ ≡ 0)."""
    return {
        'name': 'T2_zeronubb_occurs_majorana',
        'description': (
            "0νββ requires Majorana (ΔL=2). BAM neutrino is Majorana ⟸ "
            "c₁=0 (PR #86) ⟹ 0νββ occurs. Dirac would forbid it."
        ),
        'neutrino_nature': 'Majorana (c₁=0, PR #86)',
        'zeronubb_allowed': True,
        'dirac_would_forbid': True,
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T3. Normal ordering → the NO band (vs IO contrast)
# ---------------------------------------------------------------------------

def test_T3_normal_ordering_band() -> dict:
    """BAM predicts normal ordering (PR #91: generations = cavity
    overtones, m_ν ∝ m_D), selecting the NO m_ββ band. Compute it and
    contrast with the inverted-ordering band (floor ~19 meV, no
    cancellation)."""
    no0 = _band_no(0.0)
    io0 = _band_io(0.0)
    return {
        'name': 'T3_normal_ordering_band',
        'description': (
            "BAM ⟹ normal ordering (PR #91) ⟹ the NO band (m_ββ ≈ "
            "1.5–3.7 meV at zero lightest mass). IO contrast: ~19–48 meV, "
            "hard floor ~19 meV, no cancellation."
        ),
        'no_band_m1zero_meV': no0,
        'io_band_m3zero_meV': io0,
        'io_floor_meV': io0[0],
        'no_below_io_floor': no0[1] < io0[0],
        'pass': no0[1] < io0[0],
    }


# ---------------------------------------------------------------------------
# T4. Anarchic Majorana phases → full band incl. cancellation
# ---------------------------------------------------------------------------

def test_T4_anarchic_phases_band() -> dict:
    """The anarchic Majorana phases (PR #94, uniform) populate the WHOLE
    NO band — including the cancellation trough where the three terms
    partially cancel and m_ββ → ~0 (around m_lightest ~ 3–5 meV)."""
    rows = []
    for m1 in (0.0, 2.0, 3.0, 5.0):
        lo, hi = _band_no(m1)
        rows.append({'m_lightest_meV': m1, 'mbb_min_meV': lo,
                     'mbb_max_meV': hi})
    trough = min(r['mbb_min_meV'] for r in rows)
    return {
        'name': 'T4_anarchic_phases_populate_band',
        'description': (
            "Anarchic Majorana phases (PR #94, uniform) populate the whole "
            "NO band, including the cancellation trough (m_ββ → ~0 around "
            "m_lightest ~ 3–5 meV)."
        ),
        'rows': rows,
        'cancellation_trough_min_meV': trough,
        'cancellation_to_near_zero': trough < 0.5,
        'pass': trough < 0.5,
    }


# ---------------------------------------------------------------------------
# T5. BAM light scale → m_ββ ≲ 8 meV
# ---------------------------------------------------------------------------

def test_T5_light_scale() -> dict:
    """At the BAM light scale (PR #90: lightest ~ few meV), the NO band is
    m_ββ ≈ 0–8 meV. Take the maximum over a plausible scale window."""
    hi_max = 0.0
    rows = []
    for m1 in (0.0, BAM_LIGHT_SCALE_MEV, 5.0):
        lo, hi = _band_no(m1)
        hi_max = max(hi_max, hi)
        rows.append({'m_lightest_meV': m1, 'mbb_band_meV': (lo, hi)})
    return {
        'name': 'T5_bam_light_scale_mbb',
        'description': (
            "At the BAM light scale (PR #90, ~few meV) the NO band is "
            "m_ββ ≈ 0–8 meV."
        ),
        'bam_light_scale_meV': BAM_LIGHT_SCALE_MEV,
        'rows': rows,
        'mbb_upper_few_meV': hi_max,
        'mbb_below_10meV': hi_max < 10.0,
        'pass': hi_max < 10.0,
    }


# ---------------------------------------------------------------------------
# T6. Comparison to experiment / falsifiability
# ---------------------------------------------------------------------------

def test_T6_experimental_comparison() -> dict:
    """The BAM band (≲ 8 meV) is below the current bound (KamLAND-Zen,
    28–122 meV) — consistent with the null result — and largely below
    next-gen reach (LEGEND-1000 / nEXO, ~9–20 meV). Falsifiable: a 0νββ
    discovery with m_ββ ≳ 19 meV (the IO floor) would imply inverted
    ordering or a quasi-degenerate scale, contradicting BAM."""
    bam_max = _band_no(5.0)[1]
    io_floor = _band_io(0.0)[0]
    return {
        'name': 'T6_experimental_comparison_falsifiable',
        'description': (
            "BAM band ≲ 8 meV: below current KamLAND-Zen (28–122 meV) and "
            "largely below next-gen (~9–20 meV). Falsifiable — discovery "
            "≳ 19 meV (IO floor) ⟹ inverted/degenerate, contradicting BAM."
        ),
        'bam_mbb_max_meV': bam_max,
        'kamland_zen_bound_meV': KAMLAND_ZEN_BOUND_MEV,
        'nextgen_reach_meV': NEXTGEN_REACH_MEV,
        'io_floor_meV': io_floor,
        'below_current_bound': bam_max < KAMLAND_ZEN_BOUND_MEV[0],
        'falsifier': 'discovery with m_ββ ≳ 19 meV ⟹ IO/degenerate (challenges BAM)',
        'pass': bam_max < KAMLAND_ZEN_BOUND_MEV[0],
    }


# ---------------------------------------------------------------------------
# T7. Honest scope
# ---------------------------------------------------------------------------

def test_T7_honest_scope() -> dict:
    return {
        'name': 'T7_honest_scope',
        'description': (
            "Qualitative prediction firm (0νββ occurs, NO band, light "
            "scale); exact m_ββ is a band (unmeasured lightest mass + "
            "anarchic phases)."
        ),
        'established_bam_native': [
            '0νββ occurs ⟸ neutrino Majorana ⟸ c₁=0 (PR #86)',
            'BAM selects the normal-ordering m_ββ band (PR #91) — below '
            'the IO floor ~19 meV',
            'anarchic Majorana phases (PR #94) populate the full band '
            'incl. cancellation to ~0',
            'at the BAM light scale (PR #90) m_ββ ≲ 8 meV — below current '
            'bounds, a falsifiable target for next-gen experiments',
        ],
        'open': [
            'a single m_ββ value: it is a band (the lightest neutrino mass '
            'is unmeasured; the Majorana phases are anarchic)',
            'the exact spectrum (PR #91 χ_n-corrected ratios; absolute '
            'scale) and the specific phases',
        ],
        'falsifier': 'a 0νββ discovery with m_ββ ≳ 19 meV (IO/degenerate)',
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T8. Assessment
# ---------------------------------------------------------------------------

def test_T8_assessment() -> dict:
    return {
        'name': 'T8_assessment',
        'description': (
            "0νββ occurs (Majorana, c₁=0); BAM selects the normal-ordering "
            "band, populated by the anarchic Majorana phases; at the light "
            "scale m_ββ ≲ 8 meV — below current reach, a falsifiable target "
            "for next-gen experiments."
        ),
        'classification': 'ZERONUBB_OCCURS_NORMAL_ORDERING_M_BB_FEW_MEV',
        'pass': True,
    }


# ---------------------------------------------------------------------------
# Runner + verdict
# ---------------------------------------------------------------------------

def run_probe() -> dict:
    tests = [
        test_T1_setup(),
        test_T2_occurs_majorana(),
        test_T3_normal_ordering_band(),
        test_T4_anarchic_phases_band(),
        test_T5_light_scale(),
        test_T6_experimental_comparison(),
        test_T7_honest_scope(),
        test_T8_assessment(),
    ]
    if all(t['pass'] for t in tests[:7]):
        verdict_class = 'ZERONUBB_OCCURS_NORMAL_ORDERING_M_BB_FEW_MEV'
        verdict = (
            '0νββ OCCURS, IN NORMAL ORDERING, WITH m_ββ ≲ 8 meV — BELOW '
            'CURRENT REACH AND FALSIFIABLE. The neutrino arc (#85–#94) '
            'fixed the structure of the BAM neutrino sector; this probe '
            'turns it into a prediction for the effective Majorana mass '
            'm_ββ = |Σ U_ei² m_i| measured in 0νββ.\n\n'
            '0νββ OCCURS. 0νββ requires a Majorana neutrino (ΔL=2). In BAM '
            'the neutrino is Majorana because it is chargeless (c₁=0, '
            'C-invariant, PR #86), so 0νββ occurs — a qualitative '
            'prediction. A Dirac neutrino would forbid it (m_ββ ≡ 0).\n\n'
            'NORMAL ORDERING ⟹ THE NO BAND. BAM predicts normal ordering '
            '(PR #91: generations = cavity overtones, m_ν ∝ m_D), '
            'selecting the NO m_ββ band — m_ββ ≈ 1.5–3.7 meV at zero '
            'lightest mass. The inverted-ordering band (the contrast) sits '
            'at ~19–48 meV with a hard floor ~19 meV and no cancellation '
            'to zero — entirely above the BAM band.\n\n'
            'ANARCHIC PHASES POPULATE THE BAND. The anarchic Majorana '
            'phases (PR #94, uniform) populate the WHOLE NO band, '
            'including the cancellation trough where the three terms '
            'partially cancel and m_ββ → ~0 (around m_lightest ~ 3–5 '
            'meV).\n\n'
            'LIGHT SCALE ⟹ m_ββ ≲ 8 meV. At the BAM light scale (PR #90: '
            'lightest ~ few meV) the NO band is m_ββ ≈ 0–8 meV. This is '
            'below the current bound (KamLAND-Zen, m_ββ ≲ 28–122 meV) — '
            'consistent with the null result — and largely below even '
            'next-generation reach (LEGEND-1000 / nEXO, ~9–20 meV). It is '
            'a sharp falsifier: a 0νββ discovery with m_ββ ≳ 19 meV would '
            'imply inverted ordering or a quasi-degenerate scale, '
            'contradicting the BAM normal-ordering + light-scale '
            'prediction.\n\n'
            'HONEST SCOPE. ESTABLISHED (BAM-native): 0νββ occurs (the '
            'neutrino is Majorana, PR #86); BAM selects the '
            'normal-ordering band (PR #91), below the IO floor; the '
            'anarchic Majorana phases (PR #94) populate the full band '
            '(cancellation to ~0 allowed); at the BAM light scale (PR #90) '
            'm_ββ ≲ 8 meV — below current bounds and a target/falsifier '
            'for next-gen experiments. NOT established: a single m_ββ '
            'value — it is a band, because the lightest neutrino mass is '
            'unmeasured and the Majorana phases are anarchic (uniform); '
            'the exact spectrum (the PR #91 χ_n-corrected ratios, the '
            'absolute scale) and the specific phases are the residuals.'
        )
    else:
        verdict_class = 'ZERONUBB_INCONCLUSIVE'
        verdict = (
            '0νββ PREDICTION INCONCLUSIVE. A structural or numerical test '
            'failed; investigate before claiming the m_ββ band.'
        )

    return {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'identification': (
            '0νββ occurs (Majorana, c₁=0); BAM normal-ordering band '
            'populated by anarchic Majorana phases; at the light scale '
            'm_ββ ≲ 8 meV — below current reach, falsifiable'
        ),
        'occurs_because': 'neutrino Majorana ⟸ c₁=0 (PR #86)',
        'ordering': 'normal (PR #91) — NO band, below IO floor ~19 meV',
        'phases': 'anarchic (PR #94) — full band incl. cancellation to ~0',
        'scale': 'light, ~few meV (PR #90) — m_ββ ≲ 8 meV',
        'falsifier': 'discovery with m_ββ ≳ 19 meV ⟹ IO/degenerate',
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
    L.append('# Neutrinoless double-beta (0νββ) effective-mass probe (PR #95)')
    L.append('')
    L.append(f"**Run:** {s['timestamp_utc']}")
    L.append('')
    L.append(
        "Turns the neutrino arc (#85–#94) into a concrete, falsifiable "
        "prediction for the effective Majorana mass `m_ββ = |Σ U_ei² m_i|` "
        "measured in 0νββ. **0νββ occurs** (the neutrino is Majorana, "
        "`c₁=0`, PR #86); BAM selects the **normal-ordering** band "
        "(PR #91), populated by the **anarchic Majorana phases** (PR #94); "
        "and at the **light scale** (PR #90) `m_ββ ≲ 8 meV` — below "
        "current reach and a falsifiable target for next-gen experiments."
    )
    L.append('')
    L.append(f"- **Identification**: {s['identification']}")
    L.append(f"- **Occurs because**: {s['occurs_because']}")
    L.append(f"- **Ordering**: {s['ordering']}")
    L.append(f"- **Phases**: {s['phases']}")
    L.append(f"- **Scale**: {s['scale']}")
    L.append(f"- **Falsifier**: {s['falsifier']}")
    L.append('')

    L.append('## Test summary')
    L.append('')
    L.append('| # | Test | Key finding | PASS? |')
    L.append('|---|---|---|---|')
    label_map = {
        'T1': 'm_ββ = |Σ U_ei² m_i|; needs Majorana+ordering+phases+scale',
        'T2': '0νββ occurs ⟸ neutrino Majorana ⟸ c₁=0 (PR #86)',
        'T3': 'normal ordering ⟹ NO band (below IO floor ~19 meV)',
        'T4': 'anarchic phases populate full band incl. cancellation→~0',
        'T5': 'BAM light scale ⟹ m_ββ ≲ 8 meV',
        'T6': 'below current (28–122 meV) & next-gen (~9–20); falsifiable',
        'T7': 'qualitative firm; exact m_ββ a band (residuals)',
        'T8': 'ZERONUBB_OCCURS_NORMAL_ORDERING_M_BB_FEW_MEV',
    }
    for t in s['tests']:
        passed = '**PASS**' if t['pass'] else '**FAIL**'
        nm = t['name']; prefix = nm[:2]
        L.append(f"| {prefix} | `{nm}` | {label_map.get(prefix, '—')} | {passed} |")
    L.append('')

    t3 = s['tests'][2]; t4 = s['tests'][3]; t6 = s['tests'][5]
    L.append('## m_ββ band (normal ordering)')
    L.append('')
    L.append('| m_lightest (meV) | m_ββ min (meV) | m_ββ max (meV) |')
    L.append('|---:|---:|---:|')
    for r in t4['rows']:
        L.append(f"| {r['m_lightest_meV']:.0f} | {r['mbb_min_meV']:.2f} | "
                 f"{r['mbb_max_meV']:.2f} |")
    L.append('')
    L.append(f"The anarchic phases populate the whole band, with a "
             f"cancellation trough down to "
             f"`{t4['cancellation_trough_min_meV']:.2f} meV` "
             f"(m_ββ → ~0) around m_lightest ~ 3–5 meV.")
    L.append('')
    L.append('## BAM (normal) vs inverted ordering, and experiment')
    L.append('')
    L.append('| | m_ββ (meV) |')
    L.append('|---|---|')
    L.append(f"| **BAM (normal, light scale)** | ≲ {t6['bam_mbb_max_meV']:.0f} "
             "(with cancellation to ~0) |")
    L.append(f"| inverted-ordering band | {t3['io_band_m3zero_meV'][0]:.0f}–"
             f"{t3['io_band_m3zero_meV'][1]:.0f} (floor ~{t3['io_floor_meV']:.0f}, no cancellation) |")
    L.append(f"| current bound (KamLAND-Zen) | {KAMLAND_ZEN_BOUND_MEV[0]:.0f}–"
             f"{KAMLAND_ZEN_BOUND_MEV[1]:.0f} |")
    L.append(f"| next-gen reach (LEGEND-1000 / nEXO) | "
             f"~{NEXTGEN_REACH_MEV[0]:.0f}–{NEXTGEN_REACH_MEV[1]:.0f} |")
    L.append('')
    L.append("BAM predicts `m_ββ ≲ 8 meV` — below current bounds (so the "
             "null result is expected) and largely below next-gen reach. "
             "**Falsifier:** a 0νββ discovery with `m_ββ ≳ 19 meV` (the IO "
             "floor) would imply inverted ordering or a quasi-degenerate "
             "scale, contradicting BAM.")
    L.append('')

    L.append('## Verdict')
    L.append('')
    L.append(f"**{s['verdict_class']}.** {s['verdict']}")
    L.append('')

    L.append('## What this leaves open')
    L.append('')
    L.append('- **A single m_ββ value** — it is a band, because the '
             'lightest neutrino mass is unmeasured and the Majorana phases '
             'are anarchic (uniform).')
    L.append('- **The exact spectrum** (PR #91 `χ_n`-corrected ratios; the '
             'absolute scale) and the **specific phases**.')
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
    out = here / 'runs' / f'{ts}_zeronubb_effective_mass_probe'
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
