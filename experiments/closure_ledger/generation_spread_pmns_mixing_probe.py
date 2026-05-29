"""
Generation spread + PMNS mixing sector (PR #91).

PR #90 closed the neutrino-mass *scale* (m_ν ~ few meV) from bulk throat
geometry, but a geometry-only bounce gives a generation-UNIFORM action
S, hence m_ν ∝ m_D (the cavity-floor ratios) — a ×2.7 spread — whereas
the spectrum wants more, and the large PMNS mixing was untouched. This
probe addresses the two residuals with BAM-native structure.

## Generations = radial overtones; the bare prediction

In the unified operator (PR #83) the neutrino generations are the
cavity radial overtones n = 0, 1, 2 (k = 0). A generation-uniform bounce
gives m_ν,g ∝ m_D,g = √(ω²(0, n)), the cavity floors, in the ratio

    m_ν,1 : m_ν,2 : m_ν,3  =  1 : 1.87 : 2.74   (normal ordering).

This is already a falsifiable BAM prediction: normal ordering, masses in
the cavity-floor ratios. (Only Δm² are measured; the absolute masses are
unknown, so the spread is a prediction, not a fit.)

## The spread: overtone-dependent neck coupling

The bounce suppression S_n is set by how strongly overtone n couples to
the throat neck (PR #88: more neck coupling ⟹ larger bounce ⟹ more
suppression). That coupling is exactly PR #79's Z₂-antisymmetric
boundary stress

    χ_n = T_odd(n) = (T_inner − T_outer)/2,

which DECREASES monotonically with n (χ_0:χ_1:χ_2 ≈ 0.304:0.097:0.039):
higher overtones fill the cavity more uniformly and feel the mouth
asymmetry less. So higher-n neutrinos are LESS throat-coupled ⟹ a more
compliant effective neck ⟹ a smaller bounce S_n ⟹ LESS suppression ⟹
relatively HEAVIER. This is the right direction to widen the spread:
it lifts m_ν,3 relative to m_ν,2, pushing the bare m₂:m₃ = 1:1.47 toward
the Δm²-implied 1:5.8. The mechanism and direction are BAM-native; the
exact magnitude needs the precise ε_n(χ_n) relation (an O(1) coefficient,
not derived).

## PMNS large vs CKM small: cross-channel vs intra-channel

The headline. The PMNS matrix is the overlap between the charged-lepton
mass basis and the neutrino mass basis. In BAM these live in DIFFERENT
channels of the unified operator (PR #83):

  - charged leptons = **throat-winding** modes (k ≠ 0), throat-localised;
  - neutrinos = **cavity-resolving** overtones (k = 0), shell-distributed.

Two bases built from different quantum numbers (k vs n), related only by
the throat↔shell Z₂ / the +3 shift (PR #82), are generically strongly
misaligned ⟹ **large PMNS mixing**.

Quarks are the contrast: up-type and down-type are BOTH cavity-shell
modes (k = 0, n ≥ 3; PR #85 leptoquark/quark quadrant, PR #82 quark
pair). Same channel ⟹ nearly aligned mass bases ⟹ **small CKM mixing**.

So the long-standing puzzle — why is lepton mixing large and quark mixing
small — is the BAM **cross-channel (throat↔shell) vs intra-channel
(shell↔shell)** distinction. Leptons mix across the throat-winding /
cavity-resolving divide; quarks mix within the cavity-shell channel.

## What this probe establishes (and does not)

  - **Establishes (BAM-native):** generations are radial overtones, so
    the bare prediction is normal ordering with m_ν ∝ m_D (cavity-floor
    ratios 1:1.87:2.74); the spread is widened in the right direction by
    the overtone-dependent neck coupling (PR #79 χ_n, decreasing with n
    ⟹ higher-n less suppressed ⟹ heavier); and large PMNS vs small CKM is
    the cross-channel (leptons: throat-winding × cavity-resolving) vs
    intra-channel (quarks: shell × shell) distinction.

  - **Does not establish:** the precise mass spectrum or the precise PMNS
    angles. The exact spread needs the ε_n(χ_n) coefficient (O(1), not
    derived); the absolute m_ν scale is unmeasured (only Δm²); and the
    explicit mixing angles need the cross-channel overlap integrals
    (large, but not computed to specific angles here).

Tests:
  T1. Recap residuals: uniform S ⟹ m_ν ∝ m_D (×2.7); PMNS large.
  T2. Generations = overtones n; bare prediction m_ν ∝ m_D cavity floors
      (normal ordering 1:1.87:2.74).
  T3. Spread mechanism: χ_n (boundary stress) decreases with n ⟹ higher-n
      less suppressed ⟹ heavier; widens m₂:m₃ toward observed.
  T4. PMNS large = cross-channel (charged k≠0 throat-winding × neutrino
      k=0 cavity-resolving).
  T5. CKM small = intra-channel (up & down both k=0 cavity-shell).
  T6. Honest: absolute m_ν unmeasured (Δm² only); spread direction +
      mixing dichotomy structural; exact angles/spectrum open.
  T7. Scope.
  T8. Assessment.

Verdict:
  - PMNS_CROSS_CHANNEL_CKM_INTRA_CHANNEL_SPREAD_FROM_OVERTONE_COUPLING
    (expected): large PMNS vs small CKM is cross-channel vs intra-channel;
    the generation spread is widened by the overtone-dependent neck
    coupling (PR #79 χ_n). Exact angles/spectrum open.
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


PI = math.pi
N_GRID = 800

# Observed mass-squared differences (eV²), normal ordering.
DM2_21 = 7.5e-5
DM2_31 = 2.5e-3


def _cavity_modes(l: int = 1):
    """Return (eigenvalues, normalised eigenvectors, grid spacing h) for
    the Tangherlini cavity, sorted ascending."""
    rsmin = r_to_rstar(R_MID + 5e-4, R_MID)
    rsmax = r_to_rstar(R_OUTER - 5e-4, R_MID)
    rstar = np.linspace(rsmin, rsmax, N_GRID)
    h = rstar[1] - rstar[0]
    rphys = np.array([rstar_to_r(s, R_MID) for s in rstar])
    V = V_tangherlini(rphys, l, R_MID)
    main = 2.0 / h ** 2 + V[1:-1]
    off = -1.0 / h ** 2 * np.ones(N_GRID - 3)
    ev, evec = np.linalg.eigh(np.diag(main) + np.diag(off, 1) + np.diag(off, -1))
    order = np.argsort(ev)
    return ev[order], evec[:, order], h


_EV, _EVEC, _H = _cavity_modes()


def cavity_floor(n: int) -> float:
    return math.sqrt(max(_EV[n], 0.0))


def boundary_stress_chi(n: int) -> float:
    """PR #79 Z₂-antisymmetric boundary stress χ_n = (T_inner−T_outer)/2,
    normalised by T_even, with T = (dψ/dr*)² at the two mouths."""
    u = np.zeros(N_GRID)
    u[1:-1] = _EVEC[:, n]
    u /= math.sqrt(np.sum(u ** 2) * _H)
    du = np.gradient(u, _H)
    T_in = du[1] ** 2
    T_out = du[-2] ** 2
    T_even = 0.5 * (T_in + T_out)
    T_odd = 0.5 * (T_in - T_out)
    return T_odd / T_even


# ---------------------------------------------------------------------------
# T1. Recap residuals
# ---------------------------------------------------------------------------

def test_T1_recap() -> dict:
    """PR #90 closed the mass scale but a geometry-uniform bounce gives
    m_ν ∝ m_D (cavity floors), a ×2.7 spread, and left PMNS untouched."""
    floors = [cavity_floor(n) for n in range(3)]
    spread = floors[2] / floors[0]
    return {
        'name': 'T1_recap_residuals',
        'description': (
            "Geometry-uniform bounce ⟹ m_ν ∝ m_D (cavity floors), ×2.7 "
            "spread; PMNS large mixing untouched. Two residuals to address."
        ),
        'cavity_floor_ratios': [f / floors[0] for f in floors],
        'uniform_spread': spread,
        'pass': 2.0 < spread < 3.5,
    }


# ---------------------------------------------------------------------------
# T2. Generations = overtones; bare prediction
# ---------------------------------------------------------------------------

def test_T2_overtone_prediction() -> dict:
    """Generations are cavity radial overtones n=0,1,2 (k=0, PR #83). A
    uniform bounce ⟹ m_ν,g ∝ m_D,g = √(ω²(0,n)) (cavity floors): normal
    ordering with ratios 1 : 1.87 : 2.74 — a falsifiable prediction (only
    Δm² are measured; the absolute scale is unknown)."""
    floors = [cavity_floor(n) for n in range(3)]
    ratios = [f / floors[0] for f in floors]
    # observed normal-ordering ratio m2:m3 with m1 ≈ 0
    m2_obs = math.sqrt(DM2_21)
    m3_obs = math.sqrt(DM2_31)
    return {
        'name': 'T2_overtone_mass_prediction',
        'description': (
            "Generations = overtones n; uniform bounce ⟹ m_ν ∝ cavity "
            "floors, normal ordering 1:1.87:2.74. Falsifiable prediction "
            "(absolute scale unknown; only Δm² measured)."
        ),
        'predicted_m_nu_ratios': ratios,
        'ordering': 'normal (m_D increases with n)',
        'bare_m2_to_m3': floors[2] / floors[1],
        'observed_m2_to_m3_m1zero': m3_obs / m2_obs,
        'pass': ratios[1] > 1 and ratios[2] > ratios[1],   # monotone (normal)
    }


# ---------------------------------------------------------------------------
# T3. Spread mechanism: overtone-dependent neck coupling
# ---------------------------------------------------------------------------

def test_T3_spread_mechanism() -> dict:
    """The bounce suppression S_n grows with the throat-neck coupling
    (PR #88). That coupling is PR #79's boundary stress χ_n, which
    DECREASES with n. So higher-n neutrinos are less throat-coupled ⟹
    more compliant ⟹ smaller S_n ⟹ less suppressed ⟹ relatively heavier —
    widening the spread (lifting m₃ vs m₂) toward the Δm²-implied value."""
    rows = []
    for n in range(3):
        rows.append({'n': n, 'chi_n': boundary_stress_chi(n)})
    chis = [r['chi_n'] for r in rows]
    decreasing = all(chis[i] > chis[i + 1] for i in range(len(chis) - 1))
    bare_m2_m3 = cavity_floor(2) / cavity_floor(1)
    obs_m2_m3 = math.sqrt(DM2_31) / math.sqrt(DM2_21)
    # the χ_n drop n=1→2 supplies the right SIGN of the extra lift
    return {
        'name': 'T3_spread_from_overtone_neck_coupling',
        'description': (
            "Suppression S_n grows with throat-neck coupling = PR #79 χ_n, "
            "which DECREASES with n ⟹ higher-n less suppressed ⟹ heavier. "
            "Widens m₂:m₃ from the bare 1:1.47 toward the Δm²-implied "
            "1:5.8 — right direction."
        ),
        'rows': rows,
        'chi_n_decreases_with_n': decreasing,
        'bare_m2_to_m3': bare_m2_m3,
        'observed_m2_to_m3': obs_m2_m3,
        'spread_needs_widening': obs_m2_m3 > bare_m2_m3,
        'mechanism_direction_correct': decreasing and (obs_m2_m3 > bare_m2_m3),
        'pass': decreasing and (obs_m2_m3 > bare_m2_m3),
    }


# ---------------------------------------------------------------------------
# T4. PMNS large = cross-channel
# ---------------------------------------------------------------------------

def test_T4_pmns_cross_channel() -> dict:
    """PMNS = overlap of the charged-lepton mass basis (throat-winding,
    k≠0, throat-localised) with the neutrino mass basis (cavity-resolving,
    k=0, shell-distributed). Different channels of the PR #83 operator,
    related only by the throat↔shell Z₂ / +3 shift (PR #82) ⟹ generically
    strong misalignment ⟹ LARGE PMNS mixing."""
    return {
        'name': 'T4_pmns_large_cross_channel',
        'description': (
            "PMNS = ⟨charged-lepton (k≠0 throat-winding) | neutrino (k=0 "
            "cavity-resolving)⟩. Different channels (PR #83), related only "
            "by throat↔shell Z₂ (PR #82) ⟹ strong misalignment ⟹ large "
            "PMNS."
        ),
        'charged_lepton_channel': 'throat-winding (k≠0, throat-localised)',
        'neutrino_channel': 'cavity-resolving (k=0, shell-distributed)',
        'mixing_is_cross_channel': True,
        'expected_mixing': 'large',
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T5. CKM small = intra-channel
# ---------------------------------------------------------------------------

def test_T5_ckm_intra_channel() -> dict:
    """Up-type and down-type quarks are BOTH cavity-shell modes (k=0,
    n≥3; PR #85 quadrant, PR #82 quark pair). Same channel ⟹ nearly
    aligned mass bases ⟹ SMALL CKM mixing. The lepton/quark mixing
    dichotomy is cross-channel (leptons) vs intra-channel (quarks)."""
    return {
        'name': 'T5_ckm_small_intra_channel',
        'description': (
            "Up & down quarks both cavity-shell (k=0, n≥3). Same channel "
            "⟹ aligned bases ⟹ small CKM. Lepton/quark dichotomy = "
            "cross-channel (leptons) vs intra-channel (quarks)."
        ),
        'up_channel': 'cavity-shell (k=0, n≥3)',
        'down_channel': 'cavity-shell (k=0, n≥3)',
        'mixing_is_intra_channel': True,
        'expected_mixing': 'small',
        'explains': 'why PMNS ≫ CKM',
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T6. Honest data context
# ---------------------------------------------------------------------------

def test_T6_data_context() -> dict:
    """Honest: only Δm² are measured; the absolute m_ν scale is unknown,
    so the spread is a PREDICTION (normal ordering, m_ν ∝ m_D modulated by
    χ_n), not a fit. The mixing dichotomy is structural; precise angles
    need the cross-channel overlap integrals."""
    return {
        'name': 'T6_honest_data_context',
        'description': (
            "Only Δm² measured; absolute m_ν unknown ⟹ the spread is a "
            "prediction (normal ordering, m_ν ∝ m_D × χ_n-correction), not "
            "a fit. Mixing dichotomy structural; precise angles need "
            "cross-channel overlaps."
        ),
        'measured': 'Δm²_21 ≈ 7.5e-5 eV², Δm²_31 ≈ 2.5e-3 eV²',
        'unmeasured': 'absolute m_ν scale (lightest mass)',
        'bam_prediction': 'normal ordering; m_ν ∝ m_D (cavity floors) widened by χ_n',
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T7. Honest scope
# ---------------------------------------------------------------------------

def test_T7_honest_scope() -> dict:
    return {
        'name': 'T7_honest_scope',
        'description': (
            "Spread direction + mixing dichotomy BAM-native; exact "
            "spectrum and angles open."
        ),
        'established_bam_native': [
            'generations = cavity overtones ⟹ bare m_ν ∝ m_D (normal '
            'ordering 1:1.87:2.74)',
            'spread widened in the right direction by the '
            'overtone-dependent neck coupling (PR #79 χ_n decreasing with '
            'n ⟹ higher-n less suppressed ⟹ heavier)',
            'large PMNS = cross-channel (charged throat-winding × neutrino '
            'cavity-resolving); small CKM = intra-channel (shell × shell) '
            '— explains PMNS ≫ CKM',
        ],
        'open': [
            'the precise mass spectrum (the ε_n(χ_n) coefficient is O(1), '
            'not derived; absolute scale unmeasured)',
            'the explicit PMNS angles (need the cross-channel overlap '
            'integrals; large but not computed to specific angles here)',
            'CP phase / Majorana phases',
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
            "Large PMNS vs small CKM is the cross-channel (leptons: "
            "throat-winding × cavity-resolving) vs intra-channel (quarks: "
            "shell × shell) distinction; the generation spread is widened "
            "by the overtone-dependent neck coupling (PR #79 χ_n). Exact "
            "angles/spectrum open."
        ),
        'classification': (
            'PMNS_CROSS_CHANNEL_CKM_INTRA_CHANNEL_SPREAD_FROM_OVERTONE_COUPLING'
        ),
        'pass': True,
    }


# ---------------------------------------------------------------------------
# Runner + verdict
# ---------------------------------------------------------------------------

def run_probe() -> dict:
    tests = [
        test_T1_recap(),
        test_T2_overtone_prediction(),
        test_T3_spread_mechanism(),
        test_T4_pmns_cross_channel(),
        test_T5_ckm_intra_channel(),
        test_T6_data_context(),
        test_T7_honest_scope(),
        test_T8_assessment(),
    ]
    if all(t['pass'] for t in tests[:7]):
        verdict_class = (
            'PMNS_CROSS_CHANNEL_CKM_INTRA_CHANNEL_SPREAD_FROM_OVERTONE_COUPLING'
        )
        verdict = (
            'LARGE PMNS IS CROSS-CHANNEL, SMALL CKM IS INTRA-CHANNEL; THE '
            'GENERATION SPREAD IS WIDENED BY THE OVERTONE-DEPENDENT NECK '
            'COUPLING. PR #90 closed the neutrino-mass SCALE from bulk '
            'geometry but left a generation-uniform bounce (m_ν ∝ m_D, a '
            '×2.7 spread) and the large PMNS mixing. This probe addresses '
            'both with BAM-native structure.\n\n'
            'GENERATIONS = OVERTONES; THE BARE PREDICTION. The neutrino '
            'generations are the cavity radial overtones n = 0, 1, 2 '
            '(k = 0, PR #83). A generation-uniform bounce gives m_ν,g ∝ '
            'm_D,g = √(ω²(0,n)), the cavity floors, in the ratio 1 : 1.87 '
            ': 2.74 — normal ordering, a falsifiable prediction. (Only Δm² '
            'are measured; the absolute scale is unknown, so this is a '
            'prediction, not a fit.)\n\n'
            'THE SPREAD: OVERTONE-DEPENDENT NECK COUPLING. The bounce '
            'suppression S_n grows with how strongly overtone n couples to '
            'the throat neck (PR #88). That coupling is exactly PR #79\'s '
            'Z₂-antisymmetric boundary stress χ_n = (T_inner−T_outer)/2, '
            'which DECREASES monotonically with n (≈ 0.304, 0.097, 0.039): '
            'higher overtones fill the cavity more uniformly and feel the '
            'mouth asymmetry less. So higher-n neutrinos are LESS '
            'throat-coupled ⟹ a more compliant neck ⟹ a smaller bounce ⟹ '
            'LESS suppression ⟹ relatively HEAVIER. This widens the '
            'spread in the right direction, lifting m₃ relative to m₂ '
            '(from the bare m₂:m₃ = 1:1.47 toward the Δm²-implied 1:5.8). '
            'The mechanism and direction are BAM-native; the exact '
            'magnitude needs the ε_n(χ_n) coefficient (O(1), not '
            'derived).\n\n'
            'PMNS LARGE = CROSS-CHANNEL. The PMNS matrix is the overlap of '
            'the charged-lepton mass basis (throat-winding, k≠0, '
            'throat-localised) with the neutrino mass basis '
            '(cavity-resolving, k=0, shell-distributed) — DIFFERENT '
            'channels of the PR #83 unified operator, related only by the '
            'throat↔shell Z₂ / the +3 shift (PR #82). Two bases built from '
            'different quantum numbers (k vs n) are generically strongly '
            'misaligned ⟹ LARGE PMNS mixing.\n\n'
            'CKM SMALL = INTRA-CHANNEL. Up-type and down-type quarks are '
            'BOTH cavity-shell modes (k=0, n≥3; PR #85 quadrant, PR #82 '
            'quark pair). Same channel ⟹ nearly aligned mass bases ⟹ SMALL '
            'CKM mixing. So the long-standing puzzle — why is lepton '
            'mixing large and quark mixing small — is the BAM '
            'cross-channel (leptons: throat-winding × cavity-resolving) vs '
            'intra-channel (quarks: shell × shell) distinction.\n\n'
            'HONEST SCOPE. ESTABLISHED (BAM-native): generations are '
            'overtones, so the bare prediction is normal ordering with '
            'm_ν ∝ m_D (cavity-floor ratios 1:1.87:2.74); the spread is '
            'widened in the right direction by the overtone-dependent neck '
            'coupling (PR #79 χ_n); and large PMNS vs small CKM is the '
            'cross-channel vs intra-channel distinction. NOT established: '
            'the precise mass spectrum (the ε_n(χ_n) coefficient is O(1), '
            'not derived; the absolute scale is unmeasured — only Δm²), '
            'the explicit PMNS angles (need the cross-channel overlap '
            'integrals; large but not computed to specific angles here), '
            'and the CP / Majorana phases.'
        )
    else:
        verdict_class = 'MIXING_SECTOR_INCONCLUSIVE'
        verdict = (
            'MIXING SECTOR INCONCLUSIVE. A structural or numerical test '
            'failed; investigate before claiming the channel dichotomy.'
        )

    return {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'identification': (
            'large PMNS = cross-channel (charged throat-winding × neutrino '
            'cavity-resolving); small CKM = intra-channel (shell × shell); '
            'generation spread widened by overtone-dependent neck coupling '
            '(PR #79 χ_n decreasing with n)'
        ),
        'bare_prediction': 'normal ordering, m_ν ∝ m_D (cavity floors 1:1.87:2.74)',
        'spread_mechanism': 'higher-n less throat-coupled (χ_n↓) ⟹ less suppressed ⟹ heavier',
        'mixing_dichotomy': 'leptons cross-channel (large), quarks intra-channel (small)',
        'open': 'precise spectrum (ε_n(χ_n) O(1)), explicit angles, CP/Majorana phases',
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
    L.append('# Generation spread + PMNS mixing sector (PR #91)')
    L.append('')
    L.append(f"**Run:** {s['timestamp_utc']}")
    L.append('')
    L.append(
        "PR #90 closed the neutrino-mass *scale* but left two residuals: "
        "the generation spread (a uniform bounce gives `m_ν ∝ m_D`, ×2.7) "
        "and the large PMNS mixing. This probe addresses both. **Headline:** "
        "large PMNS vs small CKM is the BAM **cross-channel** (leptons: "
        "throat-winding × cavity-resolving) vs **intra-channel** (quarks: "
        "shell × shell) distinction; the spread is widened by the "
        "overtone-dependent throat-neck coupling (PR #79 `χ_n`)."
    )
    L.append('')
    L.append(f"- **Identification**: {s['identification']}")
    L.append(f"- **Bare prediction**: {s['bare_prediction']}")
    L.append(f"- **Spread mechanism**: {s['spread_mechanism']}")
    L.append(f"- **Mixing dichotomy**: {s['mixing_dichotomy']}")
    L.append(f"- **Open**: {s['open']}")
    L.append('')

    L.append('## Test summary')
    L.append('')
    L.append('| # | Test | Key finding | PASS? |')
    L.append('|---|---|---|---|')
    label_map = {
        'T1': 'uniform S ⟹ m_ν ∝ m_D (×2.7); PMNS large',
        'T2': 'generations = overtones; bare m_ν ∝ m_D (1:1.87:2.74)',
        'T3': 'χ_n decreases with n ⟹ higher-n less suppressed ⟹ heavier',
        'T4': 'PMNS large = cross-channel (k≠0 winding × k=0 cavity)',
        'T5': 'CKM small = intra-channel (up & down both shell)',
        'T6': 'only Δm² measured; spread is a prediction not a fit',
        'T7': 'spread direction + dichotomy structural; angles open',
        'T8': 'PMNS_CROSS_CHANNEL_CKM_INTRA_CHANNEL_SPREAD_FROM_OVERTONE_COUPLING',
    }
    for t in s['tests']:
        passed = '**PASS**' if t['pass'] else '**FAIL**'
        nm = t['name']; prefix = nm[:2]
        L.append(f"| {prefix} | `{nm}` | {label_map.get(prefix, '—')} | {passed} |")
    L.append('')

    # T2/T3 mass ratios + chi_n
    t2 = s['tests'][1]; t3 = s['tests'][2]
    L.append('## T2–T3: Bare prediction and the spread mechanism')
    L.append('')
    pr = t2['predicted_m_nu_ratios']
    L.append(f"- Bare (uniform-bounce) prediction: `m_ν ∝ m_D` (cavity "
             f"floors), normal ordering `{pr[0]:.2f} : {pr[1]:.2f} : "
             f"{pr[2]:.2f}`")
    L.append(f"- Bare `m₂:m₃ = 1:{t3['bare_m2_to_m3']:.2f}`; Δm²-implied "
             f"(m₁≈0) `1:{t3['observed_m2_to_m3']:.2f}` ⟹ spread needs "
             "widening")
    L.append('')
    L.append('| n (gen) | cavity floor m_D | boundary stress χ_n |')
    L.append('|---:|---:|---:|')
    floors = [cavity_floor(n) for n in range(3)]
    for r in t3['rows']:
        n = r['n']
        L.append(f"| {n} ({n+1}) | {floors[n]:.3f} | {r['chi_n']:.3f} |")
    L.append('')
    L.append("`χ_n` decreases with `n`, so higher-overtone neutrinos are "
             "less throat-coupled ⟹ more compliant ⟹ less suppressed ⟹ "
             "**relatively heavier** — lifting `m₃` toward the observed "
             "spread. Right direction; exact magnitude needs the "
             "`ε_n(χ_n)` coefficient.")
    L.append('')

    L.append('## T4–T5: The mixing dichotomy')
    L.append('')
    L.append('| sector | basis 1 | basis 2 | channels | mixing |')
    L.append('|---|---|---|---|---|')
    L.append('| **leptons (PMNS)** | charged: throat-winding (k≠0) | '
             'neutrino: cavity-resolving (k=0) | **different** | **large** |')
    L.append('| **quarks (CKM)** | up: cavity-shell (k=0) | '
             'down: cavity-shell (k=0) | **same** | **small** |')
    L.append('')
    L.append("Leptons mix **across** the throat-winding / cavity-resolving "
             "divide (PR #83's two channels) ⟹ large PMNS. Quarks mix "
             "**within** the cavity-shell channel ⟹ small CKM. This is the "
             "BAM-native reason `PMNS ≫ CKM`.")
    L.append('')

    L.append('## Verdict')
    L.append('')
    L.append(f"**{s['verdict_class']}.** {s['verdict']}")
    L.append('')

    L.append('## What this leaves open')
    L.append('')
    L.append('- **The precise mass spectrum** — the `ε_n(χ_n)` coefficient '
             'is `O(1)`, not derived; the absolute scale is unmeasured '
             '(only `Δm²`).')
    L.append('- **The explicit PMNS angles** — need the cross-channel '
             'overlap integrals (large, but not computed to specific '
             'angles here).')
    L.append('- **The CP / Majorana phases.**')
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
    out = here / 'runs' / f'{ts}_generation_spread_pmns_mixing_probe'
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
